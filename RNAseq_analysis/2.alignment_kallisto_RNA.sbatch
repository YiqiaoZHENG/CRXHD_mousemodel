#!/bin/bash

#SBATCH --job-name=yq_kallisto
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH  --array=14
#SBATCH --output=/scratch/sclab/CRX_RNAseq/slurm_outputs/slurm_array_kallisto/kallisto_%A_%a.out

# clear all module loaded in case of dependency conflicts
module purge
#load new modules required
module load mpich/3.3-python-2.7.15
module load kallisto/0.46.2-python-2.7.15

#### Usage of this script ####
# This script takes in a trimmed fastq.gz file, performs alignment and quantification of RNA transcripts.
# Input:
#   Trimmed fastq.gz file
# Output:
#   Transcript abundance matrix
#

# move into RNAseq analysis directory
basedir="/scratch/sclab/CRX_RNAseq"
cd $basedir
echo "here i am !! $PWD"
mkdir -p abundances

#directory variables
basedir="/scratch/sclab/CRX_RNAseq/"
indexdir="/scratch/sclab/CRX_RNAseq/annotation/Ensembl/Mus_musculus.GRCm38.cdna.all.idx"
trimdir="/scratch/sclab/CRX_RNAseq/trim/"
outputdir="/scratch/sclab/CRX_RNAseq/abundances/"

#### Processing file names ####
#retrieve all intput fastq file names from lookup file, one sample per line
#all but the last field contains the input fastq file names to be quantified
#the last field contains sample name
IFS=$'\n'
sample_list=($(<${basedir}kallisto_lookup.csv))
#sparse each line into array
IFS=',' read -a array <<< "${sample_list[${SLURM_ARRAY_TASK_ID}]}"
echo "total number of input fq.gz files $((${#array[@]}-1))"
echo "all input fq.fz files ${array[@]::$((${#array[@]}-1))}"
inputfiles=("${array[@]::$((${#array[@]}-1))}")
inputlist=("${inputfiles[@]/#/${trimdir}${array[${#array[@]}-1]}$'/'}")
echo "kallisto quant outputs will be written to abundances/${array[${#array[@]}-1]}"
abundancedir=${outputdir}${array[${#array[@]}-1]}.kallisto/

#### Main ####
#mm10 cdna fasta downloaded from Ensembl: https://nov2020.archive.ensembl.org/Mus_musculus/Info/Index
# ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/
# build a kallisto index of the transcriptome, Ensembl format
# this command only needs to be run once
#kallisto index -i /scratch/sclab/CRX_RNAseq/annotation/Ensembl/Mus_musculus.GRCm38.cdna.all.idx /scratch/sclab/CRX_RNAseq/annotation/Ensembl/Mus_musculus.GRCm38.cdna.all.fa.gz

#map and quantify reads against transcriptome index
kallisto quant -i $indexdir -o $abundancedir -t ${SLURM_CPUS_PER_TASK} ${inputlist[@]}

echo "ha! tihs is the end of the script!"

# clear all module loaded
module purge