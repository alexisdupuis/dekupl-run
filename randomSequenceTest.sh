#!/bin/sh

listFastas=`ls test_samples/*.fasta`

for i in $listFastas

do

    name=`basename $i '.fasta'`
    ./wgsim/wgsim -1 100 -2 100 -e 0.1 -r 0 -d 250 -N 600 $i test_samples/$name-1.fastq test_samples/$name-2.fastq
    gzip test_samples/$name-1.fastq test_samples/$name-2.fastq
    mv test_samples/$name-1.fastq.gz data/
    mv test_samples/$name-2.fastq.gz data/

done
