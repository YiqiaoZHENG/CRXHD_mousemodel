#!/bin/bash

#SBATCH --job-name=yq_trimChIP
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --output=/scratch/sclab/021420_ChIP/slurm_outputs/slurm_array_trimming/trimChIP_%A_%a.out
#SBATCH --array=0-7

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

# move into ChIPseq analysis directory
basedir="/scratch/sclab/021420_ChIP"
cd $basedir
echo "here i am !! $PWD"
mkdir -p trim/fastqc_files

#### process file names ####
# retrieve sample name from the {SLURM array ID}th line of metadata file
IFS=$'\n'
sample_list=($(<$BASEDIR/CRXHD_ChIP_lookup.txt))
sample=${sample_list[${SLURM_ARRAY_TASK_ID}]}
echo $sample

# input raw fastq.gz files 
r1="${basedir}/concat_fq/${sample}_r001.fastq.gz"
r2="${basedir}/concat_fq/${sample}_r002.fastq.gz"

# output fastq.gz files after adapter trimmer and filtering
r1_trim="${basedir}/trim/${sample}_r001_val_1.fq.gz"
r2_trim="${basedir}/trim/${sample}_r002_val_2.fq.gz"

# output dir for fastqc, make sure to mkdir before running fastqc
qc1="${basedir}/trim/fastqc_files/${sample}_fastqc"
mkdir -p $qc1
qc2="${basedir}/trim/fastqc_files/${sample}_fastqctrim"
mkdir -p $qc2

for elem in $r1 $r2 #$r1_trim $r2_trim $qc1 $qc2
do
    echo $elem
done

#### Main ####
# move into the folder where i keep all my concatenated raw fastq.gz files
cd "${basedir}/concat_fq"

# fastqc before trim
fastqc -f fastq -o $qc1 $r1 $r2
# trim and fastqc the output files
trim_galore --fastqc_args "-f fastq -o $qc2" --gzip -o "${basedir}/trim" --cores ${SLURM_CPUS_PER_TASK} --paired $r1 $r2
# qc after trim, does not need to call if fastqc already specified in trim_galore
#fastqc -f fastq -o $qc2 $r1_trim $r2_trim

echo "ha! tihs is the end of the script!"

# clear all module loaded
module purge