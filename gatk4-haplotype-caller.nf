#!/usr/bin/env nextflow

inputFile = file(params.samples)

params.outdir = "gvcf"

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
INTERVALS = file(params.genomes[params.assemly].intervals)

regions = []

INTERVALS.eachLine { str ->
        if(! str.startsWith("@") ) {
                println str.trim()
                regions << str.trim()
        }
}

// ******************
// Misc
// ******************

params.email = false

// Whether to use a local scratch disc
use_scratch = params.scratch

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

logParams(params, "nextflow_parameters_gatk4-haplotypecaller.txt")

VERSION = "0.1"

// Header log info
log.info "========================================="
log.info "GATK Best Practice for Genome-Seq calling v${VERSION}"
log.info "Section:             		HaplotypeCaller"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version:		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "========================================="

Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set {  inputHCSample }

// ------------------------------------------------------------------------------------------------------------
// Haplotype Caller for raw genomic variants
// ------------------------------------------------------------------------------------------------------------

regions = Channel.fromPath(INTERVALS).splitText(by: 1)

process runHCSample {

  tag "${indivID}|${sampleID}|${params.assembly}|batch: ${region_tag}"
  publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/HaplotypeCaller", mode: 'copy'

  // scratch use_scratch

  input:
  set indivID,sampleID,bam,bai from inputHCSample
  each region from regions

  output:
  set indivID,sampleID,file(vcf),file(vcf_index) into inputCombineVariants

  script:
  region_tag = region.trim().replace(/:/, '_')
  vcf = indivID + "_" + sampleID + "." + region_tag + ".raw_variants.g.vcf.gz"
  vcf_index = vcf + ".tbi"

  """ 
	gatk --java-options "-Xms16G -Xmx${task.memory.toGiga()}G" HaplotypeCaller \
		-R $REF \
		-I $bam \
		--intervals $region \
		-O $vcf \
		-OVI true \
		-ERC GVCF
  """  

}

VariantsPerSample = inputCombineVariants.groupTuple(by: [0,1])

process combineVariants {
	tag "${indivID}|${params.assembly}"
	publishDir "${OUTDIR}/${params.assembly}/HC", mode: 'copy'

	input:
	set indivID,sampleID,file(vcf_files),file(indices) from VariantsPerSample

	output:
	set file(gvcf),file(gvcf_index) into outputCombineVariants

	script:
	gvcf = indivID + "_" + sampleID + ".g.vcf.gz"
	gvcf_index = gvcf + ".tbi"

	// The interval files are numbered based on their correct genomic order. We can use this
	// to sort the partial gvcfs into the correct order for merging
        def sorted_vcf = [ ]
	regions.each { region -> 
		region_tag = region.trim().replace(/:/, '_')
		this_vcf = indivID + "_" + sampleID + "." + region_tag + ".raw_variants.g.vcf.gz"
		sorted_vcf << vcf_files.find { it =~ this_vcf }
	}

	"""
		gatk GatherVcfsCloud \
			-I ${sorted_vcf.join(" -I ")} \
			--output $gvcf \

		gatk IndexFeatureFile -F $gvcf
	"""
}
workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}

if (params.email) {
	workflow.onComplete {
	    def subject = 'WGS HC step finished.'
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

