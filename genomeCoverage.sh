#!/bin/sh

name=`basename $1 '.sorted.dup.recal.bam'`
rawHistogramFile=./chrY-cTetA/$name-raw-genomecov.txt
histogramFile=./chrY-cTetA/$name-genomecov.txt

mergeFile=./chrY-cTetA/$name-cTetA.txt

echo "Start of bedtools genomecov."

bedtools genomecov -d -ibam $1 > $rawHistogramFile

echo "End of bedtools genomecov."
echo "Start of filtering with awk commands."

awk '$2 > 3897400' $rawHistogramFile | awk '$2 < 3907200' > $histogramFile

echo "End of filtering."

python3 genomecov_merge.py per-base $histogramFile $mergeFile $2

rm $rawHistogramFile
rm $histogramFile

# bedtools genomecov -bg -ibam $1 | awk '$2 > 18500000' | awk '$3 < 19500000' > $histogramFile

# ttty14
# awk '$2 > 21044260' $rawHistogramFile | awk '$2 < 21046260' > $histogramFile


exit