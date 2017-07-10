#!/bin/sh

# Select the reads from a bam file aligning in a specific region and copy them into a new bam file.

if [ $# -ne 2 ]
then
    echo "This script needs two arguments :"
    echo "First argument must be the directory where the bam files are."
    echo "Second argument must be the directory name for the new bam files."
    echo "Last argument is of type chrName:[intervalStart-intervalEnd]"
    exit 1
fi

ls $1 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "The directory given in argument doesn't exist."
    exit 1
fi

# Checks if the directory chr21_bam/ already exists, otherwise creates it.

ls $2 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    mkdir $2
fi

listBams=`ls $1/*.bam`

for bam in $listBams
do
    name=`basename $bam '.sorted.dup.recal.bam'`
    bamFile=./chr21_bam/$name-chr21.bam
    ls $bamFile > /dev/null 2>&1
    # Skip the creation of the file if it already exists.
    if [ $? -eq 0 ]
    then
        continue
    fi
    echo "Selecting aligned reads from $name."
    samtools view -b $bam $2 | samtools sort > $bamFile
done

exit 0