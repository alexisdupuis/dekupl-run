#!/bin/sh

# Convert bam files into fastq files and compress them

if [ $# -ne 1 ]
then
    echo "This script needs one argument which is the directory where the bam files are."
    exit 1
fi

ls $1 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "The directory given in argument doesn't exist."
    exit 1
fi

# Checks if the directory chr21_compressed_fastq/ already exists, otherwise creates it.

ls chr21_compressed_fastq > /dev/null 2>&1
if [ $? -ne 0 ]
then
    mkdir chr21_compressed_fastq
fi

listBams=`ls $1/*.bam`

for i in $listBams
do
    name=`basename $i '.bam'`
    fastqFile=chr21_compressed_fastq/$name.fastq.gz
    ls $fastqFile > /dev/null 2>&1
    # Skip the creation of the file if it already exists.
    if [ $? -eq 0 ]
    then
        continue
    fi
    echo "Converting $name."
    samtools bam2fq $i | gzip > $fastqFile
done

exit 0