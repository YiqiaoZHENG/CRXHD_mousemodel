#!/bin/bash

#SBATCH --job-name=yq_macs2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/sclab/021420_ChIP/slurm_outputs/slurm_array_peakcall/macs2_%A_%a.out
#SBATCH --array=0-7 #because i am using sample pair, only need an array half total number of samples

# clear all module loaded in case of dependency conflicts
module purge
#load new modules required, using default
module load macs2/2.1.1.20160309

#### Usage of this script ####
# This script takes in bam files from a chip ip-input pair and perform peak calling.
# A new folder containing all peak calling outputs will be generated at the end of the script.
# Input:
#   bam file from ip sample (experimental)
#   bam file from input sample (control/background)
# Output:
#   Folder containing peak calling outputs from macs2
#

# move into ChIPseq analysis directory
basedir="/scratch/sclab/021420_ChIP"
cd $basedir
echo "here i am !! $PWD"
mkdir -p macs2peak

#### process file names ####
IFS=$'\n'
sample_list=($(<${basedir}/CRXHD_ChIP_lookup.txt))

input_sample=${sample_list[${SLURM_ARRAY_TASK_ID}*2]}
ip_sample=${sample_list[${SLURM_ARRAY_TASK_ID}*2+1]}

# input bam files
inputbam="${basedir}/alignedbam/${input_sample}.sorted.blk.bam"
ipbam="${basedir}/alignedbam/${ip_sample}.sorted.blk.bam"
# keep only protein name and replicate number as macs2 peak file prefix
sample_name=${ip_sample/_ip/}

for elem in $sample_name $ipbam $inputbam
do
    echo $elem
done

#### Main ####
# move into the directory where all the peak files will be stored
cd macs2peak
# call peak and write output to ${basedir}/macs2peak/${sample_name}_macs2/
macs2 callpeak -t $ipbam -c $inputbam -f BAMPE -g mm --keep-dup all \
               -n ${sample_name}_macs2 \
               --tempdir ${TMPDIR} \
               --outdir "${basedir}/macs2peak/${sample_name}_macs2"

echo "ha! tihs is the end of the script!"

# clear all module loaded
module purge