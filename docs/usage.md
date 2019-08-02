![](../images/ikmb_bfx_logo.png)

# Usage

## Sample sheet

This pipeline requires a sample sheet in CSV format providing information on which read files belng to which sample. 

The general format looks as follows:

```bash

IndivID;SampleID;libraryID;rgID;rgPU;platform;platform_model;Center;Date;R1;R2

```

To automatically generate such a sample sheet from a folder of PE reads, use the included script `bin/samplesheet_from_folder.rb`

```bash
ruby samplesheet_from_folder.rb --folder /path/to/reads > Samples.csv
```

Limitation: This assuems that each sample is represented by exactly one set of PE files. This will not be true if one sample was sequenced across multiple
lanes, as is the case for smaller sequencing instruments or BGI sequencers. Here, you will need to modify the first to columns so that all files belonging
to the same sample actually have the same sample name attached to them. 

## Execution

If your cluster profile is already part of this code base, you can simply run this pipeline like so:

`
nextflow run marchoeppner/wgs-calling --samples Samples.csv --assembly GRCh38 
`

The total run time depends on the size of your cluster and the speed of the underlying storage system/networking. Please be advised that, by design,
nextflows keeps all of the intermediate files in the work/ subdirectory. Running a just a dozen or so WGS data sets will therefore require a temporary
storage space of several TBs. Most of that will be removed when you delete work/ after the pipeline has finished. Still, each 30X data set will require around
300GB after all is said and done - and that is with the highly compressed CRAM format for the read alignments, instead of BAM. 


