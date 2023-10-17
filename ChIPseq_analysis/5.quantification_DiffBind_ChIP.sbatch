#!/bin/bash

#SBATCH --job-name=yq_DiffBind
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/sclab/021420_ChIP/slurm_outputs/slurm_ChIPQC/DiffBind_%j.out

# clear all module loaded in case of dependency conflicts
module purge
#load new modules required
module load r/4.0.4-python-2.7.15-java-11

#### Usage of this script ####
# This is a wrapper script to run the R script 
# 1.analyze_diffpeak.R to perform differential binding analysis
# Input:
#   Metadata spreadsheet
# Output:
#   RData file containing the differential binding analysis object
#   rds file containing the consensus peakset object
# 2. extract_fasta.R to retrieve the fasta sequences based on genomic coordinates
# Input:
#   GRanges object or bed file
# Output:
#   Fasta file with each genomic region uniquely named
# 3. GREAT_analysis.R to run GREAT analysis on a set of genomic regions
# Input:
#   GRanges object or bed file
# Output:
#   GRanges object containing the GREAT analysis result
#   Spreadsheet containing the GREAT analysis result
#

# move into ChIPseq analysis directory
basedir="/scratch/sclab/021420_ChIP"
cd $basedir
echo "here i am !! $PWD"

#run differential peak analysis, the trailing parameters specifies parallel core numbers to run the analysis
#para: --mem-per-cpu=4G --cpus-per-task=16 --exclusive --time=01:00:00
Rscript analyze_diffpeak.R ${SLURM_CPUS_PER_TASK}

#retrieve FASTA sequences based on GRanges objects
#para: --mem=8G --nodes=1 --time=01:00:00
#Rscript extract_fasta.R

#run GREAT analysis
#para: --mem=8G --nodes=1 --time=01:00:00
#Rscript GREAT_analysis.R

echo "ha! tihs is the end of the script!"

# clear all module loaded
module purge