#!/bin/sh

# Convert bam files into fastq files and compress them.

if [ $# -ne 2 ]
then
    echo "This script needs two arguments :"
    echo "The first one must be the directory where the bam files are."
    echo "The second one must be the path of the directory for the fastq files."
    exit 1
fi

ls $1 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "The directory given in argument doesn't exist."
    exit 1
fi

# Checks if the second argument correspond to a directory which already exists, otherwise creates it.

ls $2 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    mkdir $2
fi

listBams=`ls $1/*.bam`

for bam in $listBams
do
    name=`basename $bam '.bam'`
    fastqFile=$2/$name.fastq.gz
    ls $fastqFile > /dev/null 2>&1
    # Skip the creation of the file if it already exists.
    if [ $? -eq 0 ]
    then
        continue
    fi
    echo "Converting $name."
    samtools bam2fq $bam | gzip > $fastqFile
done

exit 0