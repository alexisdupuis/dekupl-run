#!/bin/sh


if [ $# -ne 1 ]
then
    echo "This script needs one argument, which is 'condition' if you want to order the samples by their condition,"
    echo "or 'contig' if you want to order them by their contig expression."
    exit 1
fi

# Creates a file for each condition, then copies the first two columns of one of
# the files since those columns are common to all the files.

touch c1_merged_genomecov.txt
cat ./chrY-genomecov/002E_chrY-merged-genomecov.txt | awk 'BEGIN {OFS="\t"} ; {print $1, $2}' > c1_merged_genomecov.txt

touch c2_merged_genomecov.txt
cat ./chrY-genomecov/002E_chrY-merged-genomecov.txt | awk 'BEGIN {OFS="\t"} ; {print $1, $2}' > c2_merged_genomecov.txt

if [ $1 = "condition" ]
then
    listFiles=`ls ./chrY-genomecov/*-merged-*`
    echo "Copying coverage values from each samples."
    for file in $listFiles
    do
        name=`basename $file '-merged-genomecov.txt'`
        awk '{print $3}' $file > tmp1.txt
        condition=$(python3 total_kmers.py get-condition ../pheno/ADES_FR_20161209.pheno.WGS_males.txt $name)
        if [ $condition -eq 1 ]
        then
            paste c1_merged_genomecov.txt tmp1.txt > tmp2.txt
            cat tmp2.txt > c1_merged_genomecov.txt
        elif [ $condition -eq 2 ]
        then
            paste c2_merged_genomecov.txt tmp1.txt > tmp2.txt
            cat tmp2.txt > c2_merged_genomecov.txt
        else
            echo "$name was not found in the phenotype file."
            exit 1
        fi
    done
elif [ $1 = "contig" ]
then
    listNames=$(python3 total_kmers.py contig-or-not 0)
    echo "Copying coverage values from each samples."
    for name in $listNames
    do
        filename=`echo $name | tr -d "," | tr -d "'" | tr -d "[" | tr -d "]"`
        file=./chrY-genomecov/$filename-merged-genomecov.txt
        awk '{print $3}' $file > tmp1.txt
        paste c1_merged_genomecov.txt tmp1.txt > tmp2.txt
        cat tmp2.txt > c1_merged_genomecov.txt
    done

    listNames=$(python3 total_kmers.py contig-or-not 1)
    for name in $listNames
    do
        filename=`echo $name | tr -d "," | tr -d "'" | tr -d "[" | tr -d "]"`
        file=./chrY-genomecov/$filename-merged-genomecov.txt
        awk '{print $3}' $file > tmp1.txt
        paste c2_merged_genomecov.txt tmp1.txt > tmp2.txt
        cat tmp2.txt > c2_merged_genomecov.txt
    done
fi

rm tmp1.txt
rm tmp2.txt

# Deletes the ten last lines of the two files because they miss informations for some samples,
# then computes the mean for each lines in those two files.

head -n -10 c1_merged_genomecov.txt > temp.txt
mv temp.txt c1_merged_genomecov.txt

head -n -10 c2_merged_genomecov.txt > temp.txt
mv temp.txt c2_merged_genomecov.txt

echo "Computing the mean for each position and each condition."

if [ $1 = "condition" ]
then
    python3 total_kmers.py genomecov-mean c1_merged_genomecov.txt c1_final_genomecov.txt -healthy
    python3 total_kmers.py genomecov-mean c2_merged_genomecov.txt c2_final_genomecov.txt -sick
else
    python3 total_kmers.py genomecov-mean c1_merged_genomecov.txt c1_final_genomecov.txt -no-contig
    python3 total_kmers.py genomecov-mean c2_merged_genomecov.txt c2_final_genomecov.txt -contig
fi

cat c1_final_genomecov.txt > file_for_R.txt
cat c2_final_genomecov.txt >> file_for_R.txt

echo "Drawing the graph of coverage per position."

if [ $1 = "condition" ]
then
    Rscript two_curves.R file_for_R.txt conditions_genomecov.jpg
else
    Rscript two_curves.R file_for_R.txt contig_or_not_contig.jpg
fi

rm c1*
rm c2*

exit 0