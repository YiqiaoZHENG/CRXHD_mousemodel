#!/bin/bash

#SBATCH --job-name=yq_homer
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=00:05:00
#SBATCH --output=/scratch/sclab/021420_ChIP/slurm_outputs/slurm_homer/homer_%A_%a.out

#set up the temporary directory
TMPDIR="/tmp/.yiqiao.zheng/${SLURM_ARRAY_JOB_ID}"
#TMPDIR="/tmp/.yiqiao.zheng/${SLURM_JOB_ID}"
mkdir -p $TMPDIR
echo "temporary working directory ${TMPDIR}"

# clear all module loaded in case of dependency conflicts
module purge
#load new modules required, using default
module load homer/4.8

#### Usage of this script ####
# This script takes in a list of regulatory region coordinates and performs genomic region annotation.
# Input:
#   Bed file containing coordinates of genomic regions. Each region should be assigned a unique name.
# Output:
#   Bed file with additional annotation fields.
#

# move into ChIPseq analysis directory
basedir="/scratch/sclab/021420_ChIP"
cd $basedir
echo "here i am !! $PWD"

#directory variables
fastadir="${basedir}DiffBind/difffasta/"
beddir="${basedir}DiffBind/diffbed/"
matrixdir="${basedir}deeptools/matrix/"
homerdir="${basedir}HOMER/"

#### process file names ####
query_bed="${beddir}hdmuts_chip_all_regions_forhome.tsv"
output_file="${beddir}hdmuts_chip_all_regions_annotation.tsv"

# annotate peaks
annotatePeaks.pl ${query_bed} mm10 > ${output_file}

echo "ha! tihs is the end of the script!"

# clear all module loaded
module purge

