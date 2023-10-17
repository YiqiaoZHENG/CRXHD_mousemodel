#!/bin/bash

#SBATCH --job-name=yq_DEanalysis
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/sclab/CRX_RNAseq/slurm_outputs/slurm_DEanalysis/DEseq2_%j.out


# clear all module loaded in case of dependency conflicts
module purge
#load new modules required, using default
module load r/4.0.4-python-2.7.15-java-11

#### Usage of this script ####
# This is a wrapper script to run the R script
# 1. annotate_kallisto.R: to convert kallisto output transcript level abundance to gene level abundance.
#                       And compile all samples in this study into a single dataframe for downstream analysis.
# Input:
#   Kallisto output abundances h5 and tsv files
# Output:
#   Gene abundance matrix
# 2. DEseq2_RNA.R: to perform differential expression analysis.
# Input:
#   Summary count matrices for all samples
# Output:
#   rds files containing P10 or P21 DESeq2 analysis results
#   Summary tables containing DE analysis for all mutant genotypes in P10 or P21
#

# move into RNAseq analysis directory
basedir="/scratch/sclab/CRX_RNAseq"
cd $basedir
echo "here i am !! $PWD"

#convert transcript to genes for all kallisto outputs, by default extract count, can specify tpm
#Note: raw count matrix should be used as input to DESeq2 analysis
#para: --nodes=1 --mem=8G --time=00:20:00
Rscript annotate_kallisto.R count

#DE analysis with DEseq2
#para: --mem-per-cpu=4G --cpus-per-task=4 --time=00:15:00
Rscript DEseq2_RNA.R ${SLURM_CPUS_PER_TASK} kallisto_estcounts_gene.ensembl.txt


echo "ha! tihs is the end of the script!"

module purge