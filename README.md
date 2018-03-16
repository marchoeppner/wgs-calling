![](images/ikmb_bfx_logo.png)

# IKMB WGS variant calling  pipeline 

This pipeline generates variant calls from genomic CRAM files. Please see the companion pipeline for generating the necessary read alignments here: http://git.ikmb.uni-kiel.de/bfx-core/NF-wgs-alignment

Tools and versions:

GATK 4.0.1.2
Samtools 1.5

## Running the pipeline

First you must clone the repository to a location visible from your cluster:

`git clone http://git.ikmb.uni-kiel.de/bfx-core/NF-wgs-calling`

You can now run the pipeline like so:

`nextflow -c path_to_repo/nextflow.config run path_to_repo/main.nf --samples samples.csv`

where `samples.csv` contains information on the input data (see section "Sample sheet").

## Valid assemblies

This pipeline supports two reference assemblies:

`GRCh37` - the gold-standard reference originally used by the 1000 Genomes consortium. This assembly includes numerous decoy sequences and some non-reference haplotypes. 

`hg38` - the current version of the human genome assembly, including all ALT contigs and HLA regions. 

Variants for both assemblies are only called in pre-defined intervals, limited to the canonical chromosomes and skipping regions such as centromeres etc. 

You MUST choose the same assembly as the one that was used for generating the alignments, obviously. 

## Sample sheet

In order to keep track of all relevant information (individual id, sample id, library id etc), the pipeline requires a CSV-formatted sample sheet as input:

  * INDIVIDUAl_ID - The ID of the individual from which the sample was derived.
  * SAMPLE_ID - The ID of the sample. Note that more than one sample can come from the same individual (e.g. tumor/normal pair)
  * BAM - Full path to read alignment
  * BAI - Full path to read alignment index

## Outputs

Unless specified otherwise (using the flag `--outdir SOME_LOCATION`), all data will be written to the folder "output". All relevant data inside will be sorted by Individual ID, sample ID and library ID (in that order). 

## Reporting

By default, the pipeline writes relevant parameters used during pipeline execution into a log file `nextflow_parameters.txt`. 

If you wish to be alerted of a completed pipeline run, you can use the optional Email flag:

`--email 'your.name@somewhere.org'`

