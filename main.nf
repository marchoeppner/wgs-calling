#!/usr/bin/env nextflow

inputFile = file(params.samples)

params.outdir = "results"

OUTDIR = file(params.outdir)

// Help message
def helpMessage() {
	log.info """
	===============================================================================
	IKMB WGS pipeline | version ${workflow.manifest.version}
	===============================================================================

	Usage: nextflow run marchoeppner/wgs-calling --samples Samples.csv --assembly hg38

	Required parameters:
	--samples	A CSV formatted sample sheet (see github documentation for details)
	--assembly 	A version of the human reference genome (see github documentation for details)

	Optional parameters:
	--email		An Email adress to which reports are sent

	Output:
	--outdir	Local directory to which all output is written (default: results)
	""".stripIndent()
}

// Show help message
if (params.help){
	helpMessage()
	exit 0
}
// *****************************************
// Assembly-specific variables and resources
// *****************************************

if (params.genomes.containsKey(params.assembly) == false) {
   exit 1, "Specified unknown genome assembly, please consult the documentation for valid assemblies."
}

REF = file(params.genomes[ params.assembly ].fasta)
DICT = REF.getBaseName() + ".dict"
DBSNP = file(params.genomes[ params.assembly ].dbsnp )
G1K = file(params.genomes[ params.assembly ].g1k )
MILLS = file(params.genomes[ params.assembly ].mills )
OMNI = file(params.genomes[ params.assembly ].omni )
HAPMAP = file(params.genomes[ params.assembly ].hapmap )
AXIOM = file(params.genomes[ params.assembly ].axiom )

INTERVALS = file(params.genomes[ params.assembly ].intervals )

// This is a bit simplistic and uses each interval , instead of pooling smaller intervals into one job. 
regions = []

INTERVALS.eachLine { str ->
        if(! str.startsWith("@") ) {
                regions << str.trim()
        }
}

summary = [:]

summary['Assembly'] = params.assembly
summary['SampleList']=  params.samples
summary['Intervals'] = INTERVALS

// ******************
// Misc
// ******************

run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

// Rules for hard filtering
SNP_RULES = params.snp_filter_rules
INDEL_RULES = params.indel_filter_rules

// Annotations to use for variant recalibration
snp_recalibration_values = params.snp_recalibration_values 
indel_recalbration_values = params.indel_recalbration_values


// Header log info
log.info "========================================="
log.info "GATK Best Practice for WGS variant calling: v${workflow.manifest.version}"
log.info "Section:             		5-dollar-genome (all-in-one)"
log.info "Commit hash:			$workflow.commitId"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version:		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "========================================="

Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set {  inputFastp }

process runFastp {

  scratch params.scratch

  input:
  set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, run_date, fastqR1, fastqR2 from inputFastp

  output:
  set val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(center), val(run_date),file("*${left}"),file("*${right}") into outputTrimAndSplit
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

    scratch params.scratch

    input:
    set val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(center), val(run_date),file(fastqR1),file(fastqR2) from inputBwa

    output:

    set indivID, sampleID, file(outfile) into runBWAOutput

    script:
    //this_chunk = fastqR1.getName()
    this_chunk = fastqR1.getName().substring(0,4)
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

	scratch params.scratch

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
	                -O /dev/stdout \
			--USE_THREADING true \
	                --SORT_ORDER coordinate | gatk SetNmMdAndUqTags \
			-I /dev/stdin \
			-O $bam_fixed \
			-R $REF \
			--IS_BISULFITE_SEQUENCE false

		samtools index $bam_fixed

	"""
}

// Mark duplicate reads. This uses a discontinuted implementation of MD to fully leverage CRAM format
process runMarkDuplicates {

    scratch params.scratch

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
                --MAX_RECORDS_IN_RAM 1000000 \
		--ASSUME_SORT_ORDER coordinate \
                --CREATE_MD5_FILE true \
		--TMP_DIR tmpdir 

		rm -Rf tmpdir
        """
}

// Generate a model for base recalibration within target intervals
process runBaseRecalibrator {

    	scratch params.scratch
	    
    	input:
    	set indivID, sampleID, dedup_bam, dedup_bai from runMarkDuplicatesOutput
	each region from regions
    
    	output:
    	set indivID, sampleID, file(recal_table) into outputBaseRecalibrator
	set indivID, sampleID, dedup_bam into BamForBQSR
    
    	script:
        region_tag = region.trim().replace(/:/, '_')
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

	input:
	set indivID,sampleID,file(reports) from ReportsBySample

	output:
	set indivID,sampleID,file(merged_report) into MergedReport

	script:
	sorted_reports = []
	regions.each { region ->
                region_tag = region.trim().replace(/:/, '_')
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

	publishDir "${OUTDIR}/${params.assembly}/cram/", mode: 'copy'

	scratch params.scratch
	    
	input:
	set indivID, sampleID, file(recal_table), file(realign_bam) from inputForApplyBQSR

	output:
	set indivID, sampleID, file(outfile_bam), file(outfile_bai) into BamForWGSStats,inputHCSample
	            
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
             -bqsr ${recal_table} \
             --output ${outfile_bam} \
             -OBM true \
	     -OVM true
    	"""	
} 

process runHCSample {

	scratch params.scratch

	input:
	set indivID,sampleID,bam,bai from inputHCSample
	each region from regions

	output:
	set val(region),file(vcf),file(vcf_index) into inputCombineVariants

	script:
	region_tag = region.trim().replace(/:/, '_')
	vcf = indivID + "_" + sampleID + "." + region_tag + ".raw_variants.g.vcf.gz"
	vcf_index = vcf + ".tbi"

	""" 
	gatk --java-options "-Xms16G -Xmx${task.memory.toGiga()}G" HaplotypeCaller \
		-R $REF \
		-I $bam \
		--intervals ${region.trim()} \
		-O $vcf \
		-OVI true \
		-ERC GVCF
	"""  

}

// gather all gvcfs for a given calling interval
VariantsPerRegion = inputCombineVariants.groupTuple()

process runGenomicsDBImport  {

        publishDir "${OUTDIR}/${params.assembly}/Variants/GenomicsDB"
	
	scratch params.scratch

	input:
	set val(region),file(vcf_list),file(indices) from VariantsPerRegion

	output:
        set region,file(genodb) into inputJoinedGenotyping

	script:
 	region_tag = region.trim().replace(/:/, '_')
	genodb = "genodb_${region_tag}"

	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" GenomicsDBImport  \
		--variant ${vcf_list.join(" --variant ")} \
		--reference $REF \
		-L ${region.trim()} \
		--genomicsdb-workspace-path $genodb \
	"""

}

// Perform genotyping on a per interval basis

process runGenotypeGVCFs {
  
	publishDir "${OUTDIR}/${params.assembly}/Variants/JointGenotypes/PerRegion"

        scratch params.scratch
  
	input:
	set region,file(genodb) from inputJoinedGenotyping
  
	output:
	file(gvcf) into inputCombineVariantsFromGenotyping
  
	script:
        region_tag = region.trim().replace(/:/, '_')
	gvcf = "genotypes.${region_tag}.g.vcf.gz"
  
	"""
 	gatk --java-options "-Xmx${task.memory.toGiga()}G" GenotypeGVCFs \
		--reference $REF \
		--only-output-calls-starting-in-intervals \
		--use-new-qual-calculator \
		--dbsnp $DBSNP \
		-V gendb://${genodb} \
               	--output $gvcf \
                -G StandardAnnotation \
		-L ${region.trim()} \
		-OVI true
  	"""
}

// Merging the scatter-gather VCF files into one file

process combineVariantsFromGenotyping {

	publishDir "${OUTDIR}/${params.assembly}/Variants/JointGenotypes", mode: 'copy'

	input:
	file(vcf_files) from inputCombineVariantsFromGenotyping.collect()

	output:
	set file(gvcf),file(gvcf_index) into (inputRecalSNP , inputRecalIndel, inputHardFilterSnp, inputHardFilterIndel )

	script:
	gvcf = "genotypes.merged.vcf.gz"
	gvcf_index = gvcf + ".tbi"

        def sorted_vcf = [ ]
	regions.each { region -> 
		region_tag = region.trim().replace(/:/, '_')
		this_vcf = "genotypes.${region_tag}.g.vcf.gz"
		sorted_vcf << vcf_files.find { it =~ this_vcf }
	}

	"""
		gatk GatherVcfsCloud \
			-I ${sorted_vcf.join(" -I ")} \
			--output $gvcf 
		gatk IndexFeatureFile -F $gvcf
	"""
}

process runHardFilterSNP {

	publishDir "${OUTDIR}/${params.assembly}/Variants/HardFilter", mode: 'copy'

	input:
	set file(vcf),file(index) from inputHardFilterSnp
	
	output:
	set file(vcf_filtered),file(vcf_filtered_index) into outputHardFilterSnp

	script:
	vcf_filtered = "genotypes.hard_filter.snp.vcf.gz"
	vcf_filtered_index = vcf_filtered + ".tbi"

	"""
		gatk SelectVariants \
                        --select-type SNP \
                        -V $vcf \
                        -O snps.vcf.gz \
                        -OVI
		 gatk VariantFiltration \
                      -R $REF \
                      -V snps.vcf.gz \
                      -O $vcf_filtered \
                      --filter-expression "${SNP_RULES}" \
                      --filter-name "hard_snp_filter" \
                      -OVI true
	"""
	
}

process runHardFilterIndel {

        publishDir "${OUTDIR}/${params.assembly}/Variants/HardFilter", mode: 'copy'

        input:
        set file(vcf),file(index) from inputHardFilterIndel

        output:
        set file(vcf_filtered),file(vcf_filtered_index) into outputHardFilterIndel

        script:
        vcf_filtered = "genotypes.hard_filter.indel.vcf.gz"
        vcf_filtered_index = vcf_filtered + ".tbi"

        """
		gatk SelectVariants \
			--select-type INDEL \
			-V $vcf \
			-O indels.vcf.gz \
			-OVI 
                gatk VariantFiltration \
                      -R $REF \
                      -V indels.vcf.gz \
                      -O $vcf_filtered \
                      --filter-expression "${INDEL_RULES}" \
                      --filter-name "hard_indel_filter" \
                       -OVI true
        """

}

process runMergeHardFilterVcf {

        publishDir "${OUTDIR}/${params.assembly}/Variants/HardFilter", mode: 'copy'

        input:
        set file(indels),file(indels_index) from outputHardFilterIndel
	set file(snps),file(snps_index) from outputHardFilterSnp


        output:
        set file(vcf_merged),file(vcf_merged_index) into outputHardFilter

	script:
	vcf_merged = "genotypes.hard_filter.merged.vcf.gz"
	vcf_merged_index = vcf_merged + ".tbi"


	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G"  MergeVcfs \
		-I $indels \
		-I $snps \
		-O $vcf_merged
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

    publishDir "${OUTDIR}/Summary/Library", mode: 'copy'
	    
    input:
    file('*') from runMarkDuplicatesOutput_QC.flatten().toList()
    file('*') from outputReportTrimming.flatten().toList()
    file('*') from CoverageStats.flatten().toList()

    output:
    file("library_multiqc*") into multiqc_report

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

	def email_fields = [:]
	email_fields['version'] = workflow.manifest.version
	email_fields['session'] = workflow.sessionId
	email_fields['runName'] = run_name
	email_fields['Samples'] = params.samples
	email_fields['success'] = workflow.success
	email_fields['dateStarted'] = workflow.start
	email_fields['dateComplete'] = workflow.complete
	email_fields['duration'] = workflow.duration
	email_fields['exitStatus'] = workflow.exitStatus
	email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
	email_fields['errorReport'] = (workflow.errorReport ?: 'None')
	email_fields['commandLine'] = workflow.commandLine
	email_fields['projectDir'] = workflow.projectDir
	email_fields['script_file'] = workflow.scriptFile
	email_fields['launchDir'] = workflow.launchDir
	email_fields['user'] = workflow.userName
	email_fields['Pipeline script hash ID'] = workflow.scriptId
	email_fields['kit'] = TARGETS
	email_fields['assembly'] = REF
	email_fields['manifest'] = workflow.manifest
	email_fields['summary'] = summary

	email_info = ""
  for (s in email_fields) {
	email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }	

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()
  
  def subject = "WGS analysis finished ($run_name)."

  if (params.email) {

  	def mqc_report = null
  	try {
        	if (workflow.success && !params.skip_multiqc) {
            		mqc_report = multiqc_report.getVal()
            		if (mqc_report.getClass() == ArrayList){
                		log.warn "[IKMB GenomeSeq] Found multiple reports from process 'multiqc', will use only one"
                		mqc_report = mqc_report[0]
                	}
        	}
    	} catch (all) {
        	log.warn "[IKMB GenomeSeq] Could not attach MultiQC report to summary email"
  	}

	def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
	def sf = new File("$baseDir/assets/sendmail_template.txt")	
    	def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    	def sendmail_html = sendmail_template.toString()

	try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}
