#!/bin/bash

file=$1
base=$(basename $1 .fastq.gz)_chunk
echo $base

zcat $file | split -d --lines 40000000 --filter='pigz -p8 > $FILE.fastq.gz' - $base

