#!/bin/sh

# Select the reads from a bam file aligning in a specific region and copy them into a new bam file,
# then create an index for this bam file.

if [ $# -ne 2 ]
then
    echo "This script needs two arguments :"
    echo "First argument must be a directory name."
    echo "Second argument is of type chrName:[intervalStart-intervalEnd]"
    exit 1
fi

ls $1 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "The directory given in argument doesn't exist."
    exit 1
fi

# Checks if the directory chr21_bam/ already exists, otherwise creates it.

ls chr21_bam > /dev/null 2>&1
if [ $? -ne 0 ]
then
    mkdir chr21_bam
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
    # echo "Indexing $bamFile."
    # samtools index -b $bamFile
done

exit 0

# alignement sur gène TTTY14 GRCh38
# samtools view -b $1 Y:18872501-19077416 | samtools sort > $bamFile
# samtools index -b $bamFile

# alignement sur gène TTTY14 GRCh37
# samtools view -b $1 Y:21034387-21239302 | samtools sort > $bamFile
# samtools index -b $bamFile

# samtools view -b $1 Y:18500000-19500000 | samtools sort > $bamFileBis
# samtools index -b $bamFileBis

# alignement région contig T et A