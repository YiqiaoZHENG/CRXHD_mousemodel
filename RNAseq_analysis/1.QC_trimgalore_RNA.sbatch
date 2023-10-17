#!/bin/bash

#SBATCH --job-name=yq_trimRNA
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --output=/scratch/sclab/CRX_RNAseq/slurm_outputs/slurm_array_trimming/trimRNA_%A_%a.out
#SBATCH --array=

# clear all module loaded in case of dependency conflicts
module purge
#load new modules required
module load trimgalore/0.6.1-python-3.6.5-java-11

#### Usage of this script ####
# This script takes in a raw fastq.gz file, performs adapter trimming, quality checking before and after trimming.
# A trimmed fastq.gz file and two fastqc folders will be generated at the end of the script.
# Input:
#   Demultiplexed fastq.gz file
# Output:
#   Trimmed fastq.gz file
#   Fastqc folders with QC matrix before and after trimming
#

# move into RNAseq analysis directory
basedir="/scratch/sclab/CRX_RNAseq"
cd $basedir
echo "here i am !! $PWD"
mkdir -p trim/fastqc_files

#### Processing file names ####
IFS=$'\n'
sample_list=($(<${basedir}trimming_lookup.csv))
#sparse each line into array
IFS=',' read -a array <<< "${sample_list[${SLURM_ARRAY_TASK_ID}]}"
echo "all input fq.fz files ${array[@]::$((${#array[@]}-1))}"
inputfiles=("${array[@]::$((${#array[@]}-1))}")
inputlist=("${inputfiles[@]/#/$basedir}")
echo "trimgalore outputs will be written to trim/${array[${#array[@]}-1]}"
outputdir=${basedir}/trim/${array[${#array[@]}-1]}/

# output dir for fastqc, make sure to mkdir before running fastqc
qc1="${basedir}/fq_qc_files/fastqc_files/${sample}_fastqc"
mkdir -p $qc1
qc2="${basedir}/fq_qc_files/fastqc_files/${sample}_fastqctrim"
mkdir -p $qc2

#### Main ####
# move into the folder where i keep all my concatenated raw fastq.gz files
cd "${basedir}/raw_fq"
#fastqc before trim
fastqc -f fastq -o $qc1 ${inputlist[@]}
#trim and fastqc the output files
trim_galore --fastqc_args "-f fastq -o $qc2" --gzip -o $outputdir --cores ${SLURM_CPUS_PER_TASK} --paired ${inputlist[@]}
#qc after trim, does not need to call if fastqc already specified in trim_galore
#fastqc -f fastq -o $qc2 $r1_trim $r2_trim

echo "ha! tihs is the end of the script!"

# clear all module loaded
module purge