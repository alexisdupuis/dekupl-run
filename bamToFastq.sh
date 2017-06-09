#!/bin/sh

# Converts BAM files into fastq files.

listBams=`ls chrY_bam/*.bam`

for i in $listBams

do

  name=`basename $i '.sorted.dup.recal.bam'`
  samtools bam2fq $i | gzip > chrY_compressed-fastq/$name.fastq.gz

done
