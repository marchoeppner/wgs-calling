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
INTERVALS = file(params.genomes[params.assembly ].intervals )
INTERVAL_CHUNKS = file(params.genomes[params.assembly ].interval_chunks )

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
regions = []
file(INTERVAL_CHUNKS).eachFile() { file ->
         regions << file
}
regions.sort()

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

process runTrimAndSplit {

   tag "${indivID}|${sampleID}|${libraryID}"

  input:
  set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, run_date, fastqR1, fastqR2 from inputFastp

  output:
  set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, run_date,file("*_R1*.trimmed.fastq.gz"),file("*_R2*.trimmed.fastq.gz") into outputTrimAndSplit
  set indivID, sampleID, libraryID, file(json),file(html) into outputReportTrimming

  script:
  left = file(fastqR1).getBaseName() + ".trimmed.fastq.gz"
  right = file(fastqR2).getBaseName() + ".trimmed.fastq.gz"
  json = indivID + "_" + sampleID + "_" + libraryID + ".fastp.json"
  html = indivID + "_" + sampleID + "_" + libraryID + ".fastp.html"

  """
	fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right -s ${task.cpus*3} -w ${task.cpus} -j $json -h $html
  """

}

inputBwa = outputTrimAndSplit.transpose( by: [0,1,2,3,4,5,6,7,8] )

// Run BWA on each trimmed chunk
process runBwa {

    tag "${indivID}|${sampleID}|${libraryID}|${rgID}|batch: ${chunk}|${params.assembly}"
    // publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/BWA/"

    //scratch use_scratch

    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center,file(left),file(right) from inputBwa

    output:

    set indivID, sampleID, file(outfile) into runBWAOutput

    script:
    outfile = sampleID + "_" + libraryID + "_" + rgID + "_" + chunk + ".aligned.cram"

    """
	bwa mem -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${REF}\\tCN:${center}" -t ${task.cpus} ${REF} $left $right | samtools view -h -f 0x2 - | samtools sort -O cram -m 7G --reference $REF - > $outfile
    """
}

runBWAOutput_grouped_by_sample = runBWAOutput.groupTuple(by: [0,1])

// Merge chunked alignments into single file per sample/individual
process runMergeCram {
    tag "${indivID}|${sampleID}|${params.assembly}"
    publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/Processing/MergeAlignments"

   //scratch use_scratch

   input:
   set indivID, sampleID, aligned_bam_list from runBWAOutput_grouped_by_sample

   output: 
   set indivID, sampleID, file(outfile_bam),file(outfile_bai) into inputMarkDuplicates 
   
   script:
   outfile_bam = sampleID + ".merged.cram"
   outfile_bai = sampleID + ".merged.cram.crai"

   """
	samtools merge -p -c -@ ${task.cpus} $outfile_bam ${aligned_bam_list.join(' ') }
	samtools index $outfile_bam
   """
   
}

// Mark duplicate reads. This uses a discontinuted implementation of MD to fully leverage CRAM format
process runMarkDuplicates {

    tag "${indivID}|${sampleID}|${params.assembly}"
    // publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/Processing/MarkDuplicates", mode: 'copy'

    //scratch use_scratch

    input:
    set indivID, sampleID, file(bam), file(bai) from inputMarkDuplicates
    
    output:
    set indivID, sampleID, file(outfile_bam),file(outfile_bai) into runMarkDuplicatesOutput
    
    file(outfile_metrics) into runMarkDuplicatesOutput_QC
    file(outfile_md5) into MarkDuplicatesMD5
    
    script:
    outfile_bam = sampleID + ".dedup.cram"
    outfile_bai = sampleID + ".dedup.cram.bai"
    outfile_md5 = sampleID + ".dedup.cram.md5"

    outfile_metrics = sampleID + "_duplicate_metrics.txt"	

    """
        gatk --java-options "-Xms32G -Xmx72G" MarkDuplicatesGATK \
                -I ${bam} \
                -O ${outfile_bam} \
                -R ${REF} \
                -M ${outfile_metrics} \
                --CREATE_INDEX true \
                --TMP_DIR \$TMPDIR \
                --MAX_RECORDS_IN_RAM 1000000 \
                --ASSUME_SORTED true \
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
	each region from  regions
    
    	output:
    	set indivID, sampleID, file(recal_table) into outputBaseRecalibrator
	set indivID, sampleID, dedup_bam into BamForBQSR
    
    	script:
        region_tag = region.getName().split("-")[0]
    	recal_table = sampleID + "." + region_tag + ".recal_table.txt" 
       
    	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" BaseRecalibrator \
		--reference ${REF} \
		--input ${dedup_bam} \
		--known-sites ${MILLS} \
		--known-sites ${DBSNP} \
                --known-sites ${G1K} \
		-L $region \
		-ip 500 \
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
	publishDir "${OUTDIR}/${params.assembly}/CRAM_FINAL", mode: 'copy'

	scratch use_scratch
	    
	input:
	set indivID, sampleID, file(recal_table), file(realign_bam) from inputForApplyBQSR

	output:
	set indivID, sampleID, file(outfile_bam), file(outfile_bai) into BamForWGSStats
	            
    	script:
    	outfile_bam = sampleID + ".clean.cram"
    	outfile_bai = sampleID + ".clean.cram.bai"
	outfile_md5 = sampleID + ".clean.cram.md5"
          
    	"""
          gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyBQSR \
             --reference ${REF} \
             --input ${realign_bam} \
	     -L $INTERVALS \
	     -ip 500 \
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

process runMultiQCLibrary {

    tag "Generating library level summary and QC plots"
	publishDir "${OUTDIR}/Summary/Library", mode: 'copy'
	    
    input:
    file('*') from runMarkDuplicatesOutput_QC.flatten().toList()
    file('*') from outputReportTrimming.flatten().toList()

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

