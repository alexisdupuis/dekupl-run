#!/bin/sh

# Create a txt file containing the coverage for each base from a bam file, then
# select a specific region, and finally merge the coverage per base in intervals.

if [ $# -ne 2 ]
then
    if [ $# -ne 4 ]
    then
        echo "This script needs at least two arguments :"
        echo "First argument must be a directory name."
        echo "Second argument is the length of the interval for the merging."
        echo "Third and last arguments (optional) are the bounds of an interval."
        exit 1
    fi
fi

ls $1 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "The directory given in argument doesn't exist."
    exit 1
fi

[ $2 -eq 0 ] 2>/dev/null
test=$?
if [ $test -gt 1 ]
then
    echo "The second argument must be an integer."
    exit 1
fi

if [ $# -ne 2 ]
then
    [ $2 -eq 0 ] 2>/dev/null
    test=$?
    if [ $test -ne 0 ]
    then
        echo "The third argument must be an integer."
        exit 1
    fi
    [ $4 -gt $3 ] 2>/dev/null
    test=$?
    if [ $test -gt 1 ]
    then
        echo "The last argument must be an integer greater than the third argument."
        exit 1
    fi
fi

# Checks if the directory chrY-genomecov/ already exists, otherwise creates it.

ls chrY-genomecov/ > /dev/null 2>&1
if [ $? -ne 0 ]
then
    mkdir chrY-genomecov/
fi

mergingInterval=$2

if [ $# -eq 2 ]
then
    bInterval=$3
    eInterval=$4
fi

listBams=`ls $1/*.bam`

for bam in $listBams
do
    name=`basename $bam '.sorted.dup.recal.bam'`
    rawHistogramFile=./chrY-genomecov/$name-raw-genomecov.txt
    histogramFile=./chrY-genomecov/$name-genomecov.txt
    mergedFile=./chrY-genomecov/$name-merged-genomecov.txt
    ls $mergedFile > /dev/null 2>&1
    # Skip the creation of the file if it already exists.
    if [ $? -eq 0 ]
    then
        continue
    fi
    echo "Start of bedtools genomecov for $name."

    bedtools genomecov -d -ibam $bam > $rawHistogramFile

    echo "End of bedtools genomecov for $name."

    if [ $# -eq 2 ]
    then
        echo "Filtering with awk commands."
        awk '$2 > $bInterval' $rawHistogramFile | awk '$2 < $eInterval' > $histogramFile
        echo "End of filtering."
    fi

    python3 genomecov_merge.py $histogramFile $mergedFile $mergingInterval

    # rm $rawHistogramFile
    # rm $histogramFile
done

exit 0