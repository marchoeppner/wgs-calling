#!/usr/bin/env nextflow

inputFile = file(params.samples)

params.outdir = "alignments"

OUTDIR = file(params.outdir)
params.report = true

// Specifies the underlying genome assembly
params.assembly = "hg38"

// *****************************************
// Assembly-specific variables and resources
// *****************************************

if (params.genomes.containsKey(params.assembly) == false) {
   exit 1, "Specified unknown genome assembly, please consult the documentation for valid assemblies."
}

REF = file(params.genomes[ params.assembly ].fasta)
DBSNP = file(params.genomes[ params.assembly ].dbsnp )
G1K = file(params.genomes[ params.assembly ].g1k )
MILLS = file(params.genomes[ params.assembly ].mills )
INTERVALS = file(params.genomes[ params.assembly ].intervals )

// ******************
// Misc
// ******************

params.format = "cram"

params.email = false

// Trimming parameters
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Define regular variables so that they can be overwritten
params.saveTrimmed = false

clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2

// Whether to use a local scratch disc
use_scratch = params.scratch

// Collect validated intervals for calling
// drastically increases parallelism
regions = Channel.fromPath(INTERVALS).splitText( by: 1)


// Make sure the Nextflow version is current enough
try {
    if( ! nextflow.version.matches(">= $params.nextflow_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

logParams(params, "nextflow_parameters-gatk4_alignment.txt")

VERSION = "0.2"

// Header log info
log.info "========================================="
log.info "GATK Best Practice for Genome-Seq Preprocessing v${VERSION}"
log.info "Section:             		Read alignment"
log.info "Commit hash:			$workflow.commitId"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version:		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "========================================="

Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set {  inputFastp }

process runFastp {

  tag "${indivID}|${sampleID}|${libraryID}"

  scratch true

  input:
  set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, run_date, fastqR1, fastqR2 from inputFastp

  output:
  set val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(center), val(run_date),file("*_R1*.trimmed.fastq.gz"),file("*_R2*.trimmed.fastq.gz") into outputTrimAndSplit
  set indivID, sampleID, libraryID, file(json),file(html) into outputReportTrimming

  script:
  left = file(fastqR1).getBaseName() + ".trimmed.fastq.gz"
  right = file(fastqR2).getBaseName() + ".trimmed.fastq.gz"
  json = indivID + "_" + sampleID + "_" + libraryID + ".fastp.json"
  html = indivID + "_" + sampleID + "_" + libraryID + ".fastp.html"

  """
	fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right -w ${task.cpus} -s ${task.cpus*3} -j $json -h $html
  """

}

inputBwa = outputTrimAndSplit.transpose( by: [9,10] )

// Run BWA on each trimmed chunk
process runBwa {

    tag "${indivID}|${sampleID}|${libraryID}|${rgID}|${this_chunk}|${params.assembly}"
    // publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/BWA/"

    scratch true

    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center,file(fastqR1),file(fastqR2) from inputBwa

    output:

    set indivID, sampleID, file(outfile) into runBWAOutput

    script:
    this_chunk = fastqR1.getName().split("-")[0].substring(0,4)
    outfile = sampleID + "_" + libraryID + "_" + rgID + "_" + this_chunk + ".aligned.bam"
    outfile_index = outfile + ".bai"
    dict_file = REF.getBaseName() + ".dict"

    """
	bwa mem -H $dict_file -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${REF}\\tCN:${center}" -t ${task.cpus} ${REF} $fastqR1 $fastqR2 | samtools sort -m 4G -@ 4 -o $outfile - 
	samtools index $outfile
    """
}

inputFixTags = runBWAOutput.groupTuple(by: [0,1])

process runFixTags {

	tag "${indivID}|${sampleID}|${params.assembly}"
	
	scratch true

	input:
    	set indivID, sampleID, file(aligned_bam_list) from inputFixTags

	output:
	set indivID,sampleID,file(bam_fixed),file(bam_fixed_bai) into inputMarkDuplicates

	script:
	bam_fixed = indivID + "_" + sampleID + ".fixed_tags.bam"
	bam_fixed_bai = bam_fixed + ".bai"

	"""
		gatk MergeSamFiles \
                    -I ${aligned_bam_list.join(' -I ')} \
                    -O merged.bam \
		    --USE_THREADING true \
                    --SORT_ORDER coordinate

		gatk SetNmMdAndUqTags \
		-I merged.bam \
		-O $bam_fixed \
		-R $REF \
		--IS_BISULFITE_SEQUENCE false

		samtools index $bam_fixed
		rm -f merged.bam

	"""
}

// Mark duplicate reads. This uses a discontinuted implementation of MD to fully leverage CRAM format
process runMarkDuplicates {

    tag "${indivID}|${sampleID}|${params.assembly}"
    // publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/Processing/MarkDuplicates", mode: 'copy'

    scratch true

    input:
    set indivID, sampleID, file(bam), file(bai) from inputMarkDuplicates
    
    output:
    set indivID, sampleID, file(outfile_bam),file(outfile_bai) into runMarkDuplicatesOutput
    
    file(outfile_metrics) into runMarkDuplicatesOutput_QC
    file(outfile_md5) into MarkDuplicatesMD5
    
    script:
    outfile_bam = indivID + "." + sampleID + ".dedup.bam"
    outfile_bai = indivID + "." + sampleID + ".dedup.bai"
    outfile_md5 = indivID + "." + sampleID + ".dedup.bam.md5"

    outfile_metrics = sampleID + "_duplicate_metrics.txt"	

    """
        gatk --java-options "-Xms4G -Xmx${task.memory.toGiga()-3}G" MarkDuplicates \
                -I ${bam} \
                -O ${outfile_bam} \
                -R ${REF} \
                -M ${outfile_metrics} \
                --CREATE_INDEX true \
                --TMP_DIR \$TMPDIR \
                --MAX_RECORDS_IN_RAM 50000 \
		--ASSUME_SORT_ORDER coordinate \
                --CREATE_MD5_FILE true
        """
}

// Generate a model for base recalibration within target intervals
process runBaseRecalibrator {

	tag "${indivID}|${sampleID}|${params.assembly}|batch: ${region_tag}"
    	// publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/Processing/BaseRecalibrator/batches", mode: 'copy'

    	scratch use_scratch
	    
    	input:
    	set indivID, sampleID, dedup_bam, dedup_bai from runMarkDuplicatesOutput
	each region from regions
    
    	output:
    	set indivID, sampleID, file(recal_table) into outputBaseRecalibrator
	set indivID, sampleID, dedup_bam into BamForBQSR
    
    	script:
        region_tag = region.getName().split("-")[0]
    	recal_table = indivID + "." + sampleID + "." + region_tag + ".recal_table.txt" 
       
    	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" BaseRecalibrator \
		--reference ${REF} \
		--input ${dedup_bam} \
		--known-sites ${MILLS} \
		--known-sites ${DBSNP} \
		--use-original-qualities \
		-L $region \
		-ip 150 \
		--output ${recal_table}
	"""
}

ReportsBySample = outputBaseRecalibrator.groupTuple(by: [0,1])

process runGatherBQSRReports {

	tag "${indivID}|${sampleID}|${params.assembly}|ALL"
        // publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/Processing/BaseRecalibrator/"

	input:
	set indivID,sampleID,file(reports) from ReportsBySample

	output:
	set indivID,sampleID,file(merged_report) into MergedReport

	script:
	sorted_reports = []
	regions.each { region ->
                region_tag = region.getName().split("-")[0]
                this_report = sampleID + "." + region_tag + ".recal_table.txt" 
                sorted_reports << reports.find { it =~ this_report }
        }	
	merged_report = indivID + "_" + sampleID + ".recal_table.txt"

	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" GatherBQSRReports \
		-I ${sorted_reports.join(' -I ')} \
		-O $merged_report 
	"""
}

inputForApplyBQSR = MergedReport.join(BamForBQSR, by: [0,1])

process runApplyBQSR {

	tag "${indivID}|${sampleID}|${params.assembly}"
	publishDir "${OUTDIR}/${params.assembly}/cram/", mode: 'copy'

	scratch use_scratch
	    
	input:
	set indivID, sampleID, file(recal_table), file(realign_bam) from inputForApplyBQSR

	output:
	set indivID, sampleID, file(outfile_bam), file(outfile_bai) into BamForWGSStats
	            
    	script:
    	outfile_bam = indivID + "." + sampleID + ".clean.cram"
    	outfile_bai = indivID + "." + sampleID + ".clean.cram.bai"
	outfile_md5 = indivID + "." + sampleID + ".clean.cram.md5"
          
    	"""
          gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyBQSR \
             --reference ${REF} \
             --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
             --use-original-qualities \
             --input ${realign_bam} \
	     -L $INTERVALS \
             -bqsr ${recal_table} \
             --output ${outfile_bam} \
             -OBM true \
	     -OVM true
    	"""	
} 


// ------------------------------------------------------------------------------------------------------------
//
// Perform a several tasks to assess QC:
// 1) Depth of coverage over targets
// 2) Generate alignment stats, insert size stats, quality score distribution
// 3) Generate hybrid capture stats
// 4) Run FASTQC to assess read quality
//
// ------------------------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------------------------
//
// Plot results with multiqc
//
// ------------------------------------------------------------------------------------------------------------

process runWgsCoverage {

	tag "${indivID}|${sampleID}|${params.assembly}"
        publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/picard_stats", mode: 'copy'

	input:
	set val(indivID),val(sampleID),val(bam),val(bai) from BamForWGSStats

	output:
	file(wgs_stats) into CoverageStats

	script:
	wgs_stats = indivID + "_" + sampleID + "_wgs_coverage.txt"

	"""
		picard CollectAlignmentSummaryMetrics \
		R=$REF \
		I=$bam \
		O=$wgs_stats
	"""
	
}

process runMultiQCLibrary {

    tag "Generating library level summary and QC plots"
	publishDir "${OUTDIR}/Summary/Library", mode: 'copy'
	    
    input:
    file('*') from runMarkDuplicatesOutput_QC.flatten().toList()
    file('*') from outputReportTrimming.flatten().toList()
    file('*') from CoverageStats.flatten().toList()

    output:
    file("library_multiqc*") into runMultiQCLibraryOutput

    script:

    """
    cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n library_multiqc *
    """
}


workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}

if (params.email) {
	workflow.onComplete {
	    def subject = 'WGS alignment finished.'
	    def recipient = params.email

	    ['mail', '-s', subject, recipient].execute() << """

	    Pipeline execution summary
	    ---------------------------
	    Completed at: ${workflow.complete}
	    Duration    : ${workflow.duration}
	    Success     : ${workflow.success}
	    workDir     : ${workflow.workDir}
	    exit status : ${workflow.exitStatus}
	    Error report: ${workflow.errorReport ?: '-'}
	    """
	}
}



//#############################################################################################################
//#############################################################################################################
//
// FUNCTIONS
//
//#############################################################################################################
//#############################################################################################################


// ------------------------------------------------------------------------------------------------------------
//
// Read input file and save it into list of lists
//
// ------------------------------------------------------------------------------------------------------------
def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}

