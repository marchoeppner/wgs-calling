#!/usr/bin/env nextflow

FOLDER = file(params.samples)

FORMATTER = "ruby $baseDir/bin/format_intervals.rb"
params.outdir = "output"

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
INTERVALS = file(params.genomes[params.assembly ].intervals )
INTERVAL_CHUNKS = file(params.genomes[params.assembly ].interval_chunks )

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

logParams(params, "nextflow_parameters_gatk4-haplotypecaller.txt")

VERSION = "0.1"

// Header log info
log.info "========================================="
log.info "GATK Best Practice for Genome-Seq calling v${VERSION}"
log.info "Section:             		Freebayes"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version:		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "========================================="

Channel.fromPath(FOLDER)
       .set {  inputFreebayes }

// ------------------------------------------------------------------------------------------------------------
// Haplotype Caller for raw genomic variants
// ------------------------------------------------------------------------------------------------------------

process runFreebayes {
  tag "ALL|batch: ${region_tag}"
  publishDir "${OUTDIR}/Freebayes/byRegion", mode: 'copy'

  // scratch use_scratch

  input:
  file(samples) from inputFreebayes
  each region from regions

  output:
  file(vcf) into outputFreebayes

  script:
  region_tag = region.getName().split("-")[0]
  vcf = region_tag + ".raw_variants.vcf.gz"
  vcf_index = vcf + ".tbi"

  """ 
	freebayes-parallel <( $FORMATTER  $region ) ${task.cpus} \
	-f $REF \
	--genotype-qualities \
	-L $samples > $vcf

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

