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


