#!/bin/sh

# Cut the file containing the contigs into multiple files containing thirty contigs at most.

if [ $# -ne 2 ]
then
    echo "This script needs two arguments :"
    echo "The first one must be the path of the file containing all the contigs"
    echo "The second one must be the path of the directory for all the files obtained from the first one."
    exit 1
fi

ls $1 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "The file given in argument doesn't exist."
    exit 1
fi

# Checks if the second argument correspond to a directory which already exists, otherwise creates it.

ls $2 > /dev/null 2>&1
if [ $? -ne 0 ]
then
    mkdir $2
fi

# Makes a copy of the file since it will be cut multiple times.

cp $1 copy-$1
copy=copy-$1
cpt=1

while [ true ]
do
    newFile=spread_contigs/thirty-contigs-$cpt.txt
    head -60 $copy > $newFile
    sed '1,60d' $copy > tmp.txt
    rm $copy
    mv tmp.txt $copy
    cpt=$(($cpt + 1))
    nbLines=`less $newFile | wc -l`
    if [ $nbLines -ne 60 ]
    then
        break
    fi
done

rm $copy

exit 0