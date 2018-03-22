#!/usr/bin/env nextflow

inputFile = file(params.samples)

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
GOLD1 = file(params.genomes[ params.assembly ].gold )
OMNI = file(params.genomes[ params.assembly ].omni )
HAPMAP = file(params.genomes[ params.assembly ].hapmap )
INTERVALS = file(params.genomes[params.assembly ].intervals )
DICT = file(params.genomes[params.assembly ].dict )

// ******************
// Misc
// ******************

params.email = false

// Whether to use a local scratch disc
use_scratch = params.scratch

// Collect validated intervals for calling
// drastically increases parallelism
regions =  []
file(INTERVALS).eachLine { line ->
	regions << line.trim()
}

logParams(params, "nextflow_parameters.txt")

VERSION = "0.1"

// Header log info
log.info "========================================="
log.info "GATK Best Practice for Genome-Seq calling v${VERSION}"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version:		${params.assembly}"
log.info "Adapter sequence used:	${adapters}"
log.info "Command Line:			$workflow.commandLine"
log.info "========================================="

Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set {  inputSplitter }

// ------------------------------------------------------------------------------------------------------------
// Haplotype Caller for raw genomic variants
// ------------------------------------------------------------------------------------------------------------

process runHCSample {

  tag "${indivID}|${params.assembly}|${region}"
  publishDir "${OUTDIR}/${params.assembly}/${indivID}/${sampleID}/HaplotypeCaller/ByIntervalRegion", mode: 'copy'

  //scratch use_scratch

  input:
  set indivID,sampleID,file(bam),file(bai) from inputHCSample
  each region from regions

  output:
  set region,file(vcf) into outputHCSample
  set region,file(vcf_index) into outputHCSampleIndex
  set indivID,sampleID,file(vcf) into inputMergeSamples
  set indivID,sampleID,file(vcf_index) into inputMergeSamplesIndex

  script:

  vcf = indivID + "_" + sampleID + ".raw_variants.${region.replace(/:/, "_")}.g.vcf.gz"
  vcf_index = vcf + ".tbi"

  """ 
	gatk --java-options "-Xms32G -Xmx64G" HaplotypeCaller \
		-R $REF \
		-I $bam \
		-L $region \
		-O $vcf \
		-OVI true \
		-ERC GVCF
  """  

}

vcfByRegion = outputHCSample.groupTuple(by: 0)
vcfindexByRegion = outputHCSampleIndex.groupTuple(by: 0)

// Import individual vcf files into a GenomicsDB database on a per chromosome basis
// From here on all samples are in the same file
process runGenomicsDBImport  {

	tag "ALL|${params.assembly}|${region}"
        publishDir "${OUTDIR}/${params.assembly}/Variants/JoinedGenotypes/PerRegion"
	
	scratch use_scratch 

	input:
        set region,file(vcf_list) from vcfByRegion
	set region_index,file(index_list) from vcfindexByRegion

	output:
        set region,file(genodb) into inputJoinedGenotyping

	script:
	genodb = "genodb_${region.replace(/:/, "_")}"

	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" GenomicsDBImport  \
		--variant ${vcf_list.join(" --variant ")} \
		--reference $REF \
		--intervals $region \
		--genomicsdb-workspace-path $genodb
	"""

}

// Perform genotyping on a per chromosome basis

process runJoinedGenotyping {
  
	tag "${region}"
	// publishDir "${OUTDIR}/Variants/JoinedGenotypes/PerRegion"
  
	input:
	set region,file(genodb) from inputJoinedGenotyping
  
	output:
	file(gvcf) into inputCombineVariantsFromGenotyping
  
	script:
  
	gvcf = "genotypes.${region.replace(/:/, "_")}.g.vcf.gz"
  
	"""
 	gatk --java-options "-Xmx${task.memory.toGiga()}G" GenotypeGVCFs \
		--reference $REF \
		--dbsnp $DBSNP \
		-V gendb://${genodb} \
               	--output $gvcf \
                -G StandardAnnotation \
		-OVI true
  	"""
}

// Merging the scatter-gather VCF files into one file

process combineVariantsFromGenotyping {
	tag "ALL"
	publishDir "${OUTDIR}/Variants/JoinedGenotypes", mode: 'copy'

	input:
	file(vcf_files) from inputCombineVariantsFromGenotyping.collect()

	output:
	file(gvcf) into (inputRecalSNP , inputRecalIndel)
	file(gvcf_index) into (inputRecalSNPIndex, inputRecalIndelIndex)
	file(gvcf_index) into combinedVariantsIndex

	script:

	gvcf = "genotypes.merged.vcf.gz"
	gvcf_index = gvcf + ".tbi"

	"""
       		vcf-concat ${vcf_files.join(" ")} | vcf-sort | bgzip > $gvcf
		tabix $gvcf
	"""
}

process runRecalibrationModeSNP {

	tag "ALL"
	// publishDir "${OUTDIR}/Variants/Recal"

	input:
	file(vcf) from inputRecalSNP
	file(index) from inputRecalSNPIndex

	output:
  	set file(recal_file),file(tranches),file(rscript),file(snp_file),file(snp_index) into inputRecalSNPApply

	script:
	snp_file = "genotypes.merged.snps.vcf.gz"
	snp_index = snp_file + ".tbi"
	recal_file = "genotypes.recal_SNP.recal"
  	tranches = "genotypes.recal_SNP.tranches"
  	rscript = "genowtypes.recal_SNP.R"

  	"""

	gatk --java-options "-Xmx${task.memory.toGiga()}G" SelectVariants \
		-R $REF \
		-V $vcf \
		-select-type SNP \
		-O $snp_file

	gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
		-R $REF \
		-V $snp_file \
               	-O $recal_file \
                --tranches-file $tranches \
	        --rscript-file $rscript \
		-an MQ -an MQRankSum -an ReadPosRankSum -an FS -an DP -an ReadPosRankSum -an InbreedingCoeff -an QD \
                -mode SNP \
		--resource hapmap,known=false,training=true,truth=true,prior=15.0:$HAPMAP \
		--resource omni,known=false,training=true,truth=true,prior=12.0:$OMNI \
		--resource 1000G,known=false,training=true,truth=false,prior=10.0:$G1K \
		--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
                -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  	"""
}

process runRecalibrationModeIndel {

	tag "ALL"
	// publishDir "${OUTDIR}/Variants/Recal"

	input:
	file(vcf) from inputRecalIndel
	file(index) from inputRecalIndelIndex

	output:
	set file(recal_file),file(tranches),file(rscript),file(indel_file),file(indel_index) into inputRecalIndelApply

	script:
	indel_file = "genotypes.merged.indel.vcf.gz"
	indel_index = indel_file + ".tbi"
	recal_file = "genotypes.recal_Indel.recal"
	tranches = "genotypes.recal_Indel.tranches"
	rscript = "genotypes.recal_Indel.R"

	"""
		
	gatk --java-options "-Xmx${task.memory.toGiga()}G" SelectVariants \
		-R $REF \
		-V $vcf \
		-select-type INDEL \
		-select-type MIXED \
		-select-type MNP \
		-O $indel_file \
		-OVI true

        gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
       	        -R $REF \
                -V $indel_file \
               	-O $recal_file \
       	        --tranches-file $tranches \
                --rscript-file $rscript \
               	-an MQRankSum -an SOR -an ReadPosRankSum -an FS -an DP -an InbreedingCoeff -an MQ -an QD \
       	        -mode INDEL \
                --resource mills,known=false,training=true,truth=true,prior=15.0:$GOLD1 \
               	--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
                -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	"""
}

process runRecalSNPApply {

	tag "ALL"
	// publishDir "${OUTDIR}/Variants/Filtered"

	input:
	set file(recal_file),file(tranches),file(rscript),file(gvcf),file(gvcf_index) from inputRecalSNPApply

	output:
	file vcf_snp   into outputRecalSNPApply

	script:
 
	vcf_snp = "genotypes.recal_SNP.vcf.gz"

	"""
		gatk IndexFeatureFile -F $recal_file
	 	gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
			-R $REF \
			-V $gvcf \
		        --recal-file $recal_file \
                	--tranches-file $tranches \
			-mode SNP \
			--ts-filter-level 99.0 \
			-O $vcf_snp	
  	"""
}

process runRecalIndelApply {

	tag "ALL"
  	// publishDir "${OUTDIR}/Variants/Recal"

	input:
  	set file(recal_file),file(tranches),file(rscript),file(gvcf),file(gvcf_index) from inputRecalIndelApply

  	output:
  	set file(vcf_indel),file(gvcf_index) into outputRecalIndelApply

  	script:

	vcf_indel = "genotypes.recal_Indel.vcf.gz"
	gvcf_index = vcf_indel + ".tbi"

	"""
        	gatk IndexFeatureFile -F $recal_file
        	gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
                	-R $REF \
	                -V $gvcf \
        	        --recal-file $recal_file \
                	--tranches-file $tranches \
	                -mode INDEL \
        	        --ts-filter-level 99.0 \
                	-O $vcf_indel \
			-OVI true
	  """
}

process runVariantFiltrationIndel {

	tag "ALL"
	// publishDir "${OUTDIR}//Variants/Filtered"

  	input:
	set file(gvcf),file(gvcf_index) from outputRecalIndelApply

  	output:
  	file(filtered_gvcf) into outputVariantFiltrationIndel

  	script:

  	filtered_gvcf = "genotypes.recal_Indel.filtered.vcf.gz"

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

inputCombineVariants = outputVariantFiltrationIndel.mix(outputRecalSNPApply)

process runCombineVariants {

	tag "ALL"
  	publishDir "${OUTDIR}/Variants/Final", mode: 'copy'

	input: 
    	set file(indel),file(snp) from inputCombineVariants.collect()

	output:
	set file(merged_file),file(merged_file_index) into inputLeftNormalize

	script:
	merged_file = "merged_callset.vcf.gz"
	merged_file_index = merged_file + ".tbi"

	"""
		gatk IndexFeatureFile -F $indel
		gatk IndexFeatureFile -F $snp

		gatk SortVcf -I $indel -O indels.sorted.vcf.gz	
		gatk IndexFeatureFile -F indels.sorted.vcf.gz

		gatk SortVcf -I $snp -O snps.sorted.vcf.gz
		gatk IndexFeatureFile -F snps.sorted.vcf.gz

		gatk MergeVcfs \
		-I=indels.sorted.vcf.gz \
		-I=snps.sorted.vcf.gz \
		-O=merged.vcf.gz \
		-R=$REF \
	
		gatk IndexFeatureFile -F merged.vcf.gz

		gatk SelectVariants \
		-R $REF \
		-V merged.vcf.gz \
		-OVI true \
		-O $merged_file \
		--remove-unused-alternates true \
		--exclude-non-variants true
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

