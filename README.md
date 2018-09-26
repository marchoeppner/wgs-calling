![](images/ikmb_bfx_logo.png)

# IKMB WGS variant calling  pipeline 

This pipeline generates variant calls from short read files. 

Tools and versions:

GATK 4.0.9.0
Samtools 1.5

## Running the pipeline

First you must clone the repository to a location visible from your cluster:

`git clone http://git.ikmb.uni-kiel.de/bfx-core/NF-wgs-calling`

The pipeline consists, for practical reasons, of three independent workflows: read alignment (produces a WGS cram), haplotype calling (produces a GVCF) and finally joint genotyping (produces a final VCF). 

You can now run the pipeline like so:

`nextflow -c path_to_repo/nextflow.config run path_to_repo/gatk-alignment.nf --samples Samples.csv`

where `Samples.csv` contains information on the input data (see section "Sample sheet").

## Valid assemblies

This pipeline supports two reference assemblies:

`GRCh37` - the gold-standard reference originally used by the 1000 Genomes consortium. This assembly includes numerous decoy sequences and some non-reference haplotypes. 

`hg38` - the current version of the human genome assembly, including all ALT contigs and HLA regions. 

Variants for both assemblies are only called in pre-defined intervals, limited to the canonical chromosomes and skipping regions such as centromeres etc. 

Obviously, you should choose the same assembly for all steps of the pipeline. 

## Sample sheet

Each pipeline requires a specific samplesheet. Scripts are included under bin/ that can produce compliant inputs from a folder that contains relevant data.

### gatk4-alignment.nf
bin/samplesheet_from_folder.rb -f /path/to/fastq

### gatk4-haplotype-caller.nf
bin/samplesheet_from_cram.rb -f /path/to/cram

### gatk4-joint-genotyping.nf 
bin/samplesheet_from_gvcf.rb -f path/to/gvcf

## Outputs

Each pipeline create an output directory: alignment, gvcf and genotypes

## Reporting

By default, the pipeline writes relevant parameters used during pipeline execution into a log file `nextflow_parameters.txt`. 

If you wish to be alerted of a completed pipeline run, you can use the optional Email flag:

`--email 'your.name@somewhere.org'`

