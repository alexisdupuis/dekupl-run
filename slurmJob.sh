#!/bin/sh

#SBATCH --job-name=selectReads
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=2048M

./selectReads.sh /storage/ipl/collaboration/WGS 21