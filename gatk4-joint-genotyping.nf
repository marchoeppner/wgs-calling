#!/usr/bin/env nextflow

FOLDER = file(params.folder)

params.outdir = "genotypes"

OUTDIR = file(params.outdir)

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
OMNI = file(params.genomes[ params.assembly ].omni )
HAPMAP = file(params.genomes[ params.assembly ].hapmap )
AXIOM = file(params.genomes[ params.assembly ].axiom )
INTERVALS = file(params.genomes[params.assembly ].intervals )
INTERVAL_CHUNKS = file(params.genomes[params.assembly ].interval_chunks )

// Annotations to use for variant recalibration
snp_recalibration_values = params.snp_recalibration_values 
indel_recalbration_values = params.indel_recalbration_values

// ******************
// Misc
// ******************

params.email = false

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

logParams(params, "nextflow_parameters-gatk4_joint_genotyping.txt")

VERSION = "0.1"

// Header log info
log.info "========================================="
log.info "GATK Best Practice for Genome-Seq calling v${VERSION}"
log.info "Section: 			Joint Variant Calling"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version:		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "========================================="

GVCF = Channel.fromPath(FOLDER + "/*.g.vcf.gz")
GVCF_INDICES = Channel.fromPath(FOLDER + "/*.g.vcf.gz.tbi")

// ------------------------------------------------------------------------------------------------------------
// Haplotype Caller for raw genomic variants
// ------------------------------------------------------------------------------------------------------------


// Import individual vcf files into a GenomicsDB database on a per chromosome basis
// From here on all samples are in the same file
process runGenomicsDBImport  {

	tag "ALL|${params.assembly}|batch: ${region_tag}"
        publishDir "${OUTDIR}/${params.assembly}/Variants/GenomicsDB"
	
	//scratch use_scratch 

	input:
	file(vcf_list) from GVCF.collect()
	file(indices) from GVCF_INDICES.collect()

	each region from regions
	
	output:
        set region,file(genodb) into inputJoinedGenotyping

	script:
 	region_tag = region.getName().split("-")[0]	
	genodb = "genodb_${region_tag}"

	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" GenomicsDBImport  \
		--variant ${vcf_list.join(" --variant ")} \
		--reference $REF \
		-L $region \
		--genomicsdb-workspace-path $genodb
	"""

}

// Perform genotyping on a per interval basis

process runGenotypeGVCFs {
  
	tag "ALL|${params.assembly}|batch: ${region_tag}"
	publishDir "${OUTDIR}/${params.assembly}/Variants/JointGenotypes/PerRegion"
  
	input:
	set region,file(genodb) from inputJoinedGenotyping
  
	output:
	file(gvcf) into inputCombineVariantsFromGenotyping
  
	script:
        region_tag = region.getName().split("-")[0]  
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
		-L $region \
		-OVI true
  	"""
}

// Merging the scatter-gather VCF files into one file

process combineVariantsFromGenotyping {
	tag "ALL"
	publishDir "${OUTDIR}/${params.assembly}/Variants/JoinedGenotypes", mode: 'copy'

	input:
	file(vcf_files) from inputCombineVariantsFromGenotyping.collect()

	output:
	set file(gvcf),file(gvcf_index) into (inputRecalSNP , inputRecalIndel )

	script:
	gvcf = "genotypes.merged.vcf.gz"
	gvcf_index = gvcf + ".tbi"

        def sorted_vcf = [ ]
	regions.each { region -> 
		region_tag = region.getName().split("-")[0]
		this_vcf = "genotypes.${region_tag}.g.vcf.gz"
		sorted_vcf << vcf_files.find { it =~ this_vcf }
	}

	"""
		gatk GatherVcfsCloud \
			-I ${sorted_vcf.join(" -I ")} \
			--output $gvcf \

		gatk IndexFeatureFile -F $gvcf
	"""
}

process runRecalibrationModeSNP {

	tag "ALL"
	// publishDir "${OUTDIR}/${params.assembly}/Variants/Recal"

	input:
	set file(vcf),file(index) from inputRecalSNP

	output:
  	set file(recal_file),file(tranches) into inputRecalSNPApply

	script:
	recal_file = "genotypes.recal_SNP.recal"
  	tranches = "genotypes.recal_SNP.tranches"

  	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
		-R $REF \
		-V $vcf \
               	-O $recal_file \
                --tranches-file $tranches \
		-an ${snp_recalibration_values.join(' -an ')} \
                -mode SNP \
		--trust-all-polymorphic \
		-L $INTERVALS \
		--resource hapmap,known=false,training=true,truth=true,prior=15:$HAPMAP \
		--resource omni,known=false,training=true,truth=true,prior=12:$OMNI \
		--resource 1000G,known=false,training=true,truth=false,prior=10:$G1K \
		--resource dbsnp,known=true,training=false,truth=false,prior=7:$DBSNP \
                -tranche ${params.snp_recalibration_tranche_values.join(' -tranche ')} \
  	"""
}

process runRecalibrationModeIndel {

	tag "ALL"
	// publishDir "${OUTDIR}/${params.assembly}/Variants/Recal"

	input:
	set file(vcf),file(index) from inputRecalIndel

	output:
	set file(recal_file),file(tranches),file(vcf),file(index) into inputRecalIndelApply

	script:

	recal_file = "genotypes.recal_Indel.recal"
	tranches = "genotypes.recal_Indel.tranches"

	"""
        gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
       	        -R $REF \
                -V $vcf \
               	-O $recal_file \
       	        --tranches-file $tranches \
		-an ${indel_recalbration_values.join(' -an ')} \
       	        -mode INDEL \
		-L $INTERVALS \
                --resource mills,known=false,training=true,truth=true,prior=12.0:$MILLS \
		--resource axiomPoly,known=false,training=true,truth=false,prior=10:$AXIOM \
               	--resource dbsnp,known=true,training=false,truth=false,prior=2:$DBSNP \
		-tranche ${params.indel_recalibration_tranche_values.join(' -tranche ')} \
		-max-gaussians 3 \
		--trust-all-polymorphic
	"""
}

process runRecalIndelApply {
        tag "ALL"
        // publishDir "${OUTDIR}/${params.assembly}/Variants/Recal"

        input:
        set file(recal_file),file(tranches),file(gvcf),file(gvcf_index) from inputRecalIndelApply

        output:
        set file(vcf_indel),file(vcf_indel_index) into VcfRecalSNPApply

        script:

        vcf_indel = "genotypes.recal_Indel.vcf.gz"
        vcf_indel_index = vcf_indel + ".tbi"

        """
                gatk IndexFeatureFile -F $recal_file
                gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
                        -R $REF \
                        -V $gvcf \
                        --recal-file $recal_file \
                        --tranches-file $tranches \
                        -mode INDEL \
                        --ts-filter-level ${params.indel_filter_level} \
                        -O $vcf_indel \
                        -OVI true
          """
}

process runRecalSNPApply {

	tag "ALL"
	publishDir "${OUTDIR}/${params.assembly}/Variants/Filtered"

	input:
	set file(recal_file),file(tranches) from inputRecalSNPApply
        set file(vcf_indel),file(vcf_indel_index) from VcfRecalSNPApply

	output:
	set file(vcf_recalibrated),file(vcf_recalibrated_index) into inputFilterIndel

	script:
 
	vcf_recalibrated = "genotypes.recal_Indel.recal_SNP.vcf.gz"
	vcf_recalibrated_index = vcf_recalibrated + ".tbi"

	"""
		gatk IndexFeatureFile -F $recal_file
	 	gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
			-R $REF \
			-V $vcf_indel \
		        --recal-file $recal_file \
                	--tranches-file $tranches \
			-mode SNP \
			--ts-filter-level ${params.snp_filter_level} \
			-O $vcf_recalibrated \
			-OVI true
  	"""
}

process runVariantFiltrationIndel {

	tag "ALL"
	publishDir "${OUTDIR}/${params.assembly}/Variants/Final"

  	input:
	set file(gvcf),file(gvcf_index) from inputFilterIndel

  	output:
  	set file(filtered_gvcf),file(filtered_gvcf_index) into inputCollectMetrics

  	script:

  	filtered_gvcf = "genotypes.recal_Indel.recal_SNP.filtered.vcf.gz"
	filtered_gvcf_index = filtered_gvcf + ".tbi"

	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantFiltration \
                -R $REF \
                -V $gvcf \
		--filter-expression "QD < 2.0" \
		--filter-name "QDFilter" \
                -O $filtered_gvcf \
		-OVI true
  	"""
}

process runCollectVariantCallingMetrics {

	tag "ALL"
        publishDir "${OUTDIR}/${params.assembly}/Variants/QC"

	input:
        set file(filtered_gvcf),file(filtered_gvcf_index) from inputCollectMetrics

	output:
	set file(metrics_details),file(metrics_summary) into outputCallingMetrics

	script:
	metrics_details = "genotypes.variant_calling_detail_metrics"
	metrics_summary = "genotypes.variant_calling_summary_metrics"

	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" CollectVariantCallingMetrics \
		-I ${filtered_gvcf} \
		-O genotypes \
		-TI $INTERVALS \
		--DBSNP $DBSNP
	"""
	
}


workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}

if (params.email) {
	workflow.onComplete {
	    def subject = 'WGS joint calling finished.'
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

