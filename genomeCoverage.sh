#!/bin/sh

# Create a txt file containing the coverage for each base from a bam file, then
# select a specific region, and finally merge the coverage per base in intervals.

if [ $# -ne 4 ]
then
    echo "This script needs four arguments :"
    echo "First argument must be a directory name."
    echo "Second argument is the length of the interval for the merging."
    echo "Third and last arguments are the bounds of an interval."
    exit 1
fi

ls $1 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "The directory given in argument doesn't exist."
    exit 1
fi

[ $2 -eq 0 ] 2>/dev/null
if [ $? -gt 1 ]
then
    echo "The second argument must be an integer."
    exit 1
fi

[ $2 -eq 0 ] 2>/dev/null
if [ $? -gt 1 ]
then
    echo "The third argument must be an integer."
    exit 1
fi
[ $4 -gt $3 ] 2>/dev/null
if [ $? -ne 0 ]
then
    echo "The last argument must be an integer greater than the third argument."
    exit 1
fi

# Checks if the directory chrY-genomecov/ already exists, otherwise creates it.

ls ./chrY-genomecov/ > /dev/null 2>&1
if [ $? -ne 0 ]
then
    mkdir ./chrY-genomecov/
fi

mergingInterval=$2

listBams=`ls $1/*.bam`

for bam in $listBams
do
    name=`basename $bam '.sorted.dup.recal.bam'`
    rawHistogramFile=/nas/chrY-genomecov/$name-raw-genomecov.txt
    ls $rawHistogramFile > /dev/null 2>&1
    # Skip the creation of the raw genomecov file if it already exists.
    if [ $? -eq 0 ]
    then
        continue
    else
        echo "Start of bedtools genomecov for $name."
        bedtools genomecov -d -ibam $bam > $rawHistogramFile
        echo "End of bedtools genomecov for $name."
    fi
done

listRaws=`ls /nas/chrY-genomecov/*-raw-*`

for raw in $listRaws
do
    name=`basename $raw '-raw-genomecov.txt'`
    rawHistogramFile=/nas/chrY-genomecov/$name-raw-genomecov.txt
    histogramFile=/nas/chrY-genomecov/$name-genomecov.txt
    ls $histogramFile > /dev/null 2>&1
    # Skip the creation of the cut genomecov file if it already exists.
    if [ $? -eq 0 ]
    then
        continue
    else
        echo "Filtering with awk commands for $name."
        awk -v bInterval=$3 '$2 > bInterval' $rawHistogramFile | awk -v eInterval=$4 '$2 < eInterval' > $histogramFile
        echo "End of filtering for $name."
    fi
done

listHists=`ls /nas/chrY-genomecov/*Y-genomecov.txt`

for hist in $listHists
do
    name=`basename $hist '-genomecov.txt'`
    histogramFile=/nas/chrY-genomecov/$name-genomecov.txt
    mergedFile=./chrY-genomecov/$name-merged-genomecov.txt
    ls $mergedFile > /dev/null 2>&1
    # Skip the creation of the merged file if it already exists.
    if [ $? -eq 0 ]
    then
        continue
    else
        echo "Merging positions into intervals for $name."
        python3 genomecov_merge.py $histogramFile $mergedFile $mergingInterval
        echo "Merging done for $name."
    fi
done

exit 0