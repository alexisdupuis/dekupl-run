#!/bin/sh

# Select the reads from a BAM file aligning in a specific region and copy them into a new BAM file,
# then create an index for this BAM file.

name=`basename $1 '.sorted.dup.recal.bam'`
bamFile=./chrY-cTetA/$name-cTetA.bam
# bamFileBis=./chrY_ttty14-alignment/$name-ttty14-37alignment.bam

# alignement sur gène TTTY14 GRCh38
# samtools view -b $1 Y:18872501-19077416 | samtools sort > $bamFile
# samtools index -b $bamFile

# alignement sur gène TTTY14 GRCh37
# samtools view -b $1 Y:21034387-21239302 | samtools sort > $bamFile
# samtools index -b $bamFile

#samtools view -b $1 Y:18500000-19500000 | samtools sort > $bamFileBis
# samtools index -b $bamFileBis

# alignement région contig T et A
samtools view -b $1 Y:3901400-3903200 | samtools sort > $bamFile
samtools index -b $bamFile

exit