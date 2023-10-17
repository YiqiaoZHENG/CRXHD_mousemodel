#!/bin/bash

#SBATCH --job-name=yq_picard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --output=/scratch/sclab/021420_ChIP/slurm_outputs/slurm_array_samtobam/traimalign_%A_%a.out
#SBATCH --array=0-7

# clear all module loaded in case of dependency conflicts
module purge

#load new modules required
module load bowtie2/2.3.5-python-3.6.5
module load samtools/1.9
module load bedtools/2.27.1
module load picard/2.21.4-java-11

#### Usage of this script ####
# This script takes in a trimmed fastq.gz file and aligns it to the mouse genome (mm10), generating an initial sam file.
# After sam to bam file convertion, read pairs that are not uniquely mapped are then filtered out.
# Reads mapped to mitochondria genome and overlapped with ENCODE documented blacklist regions are then filtered out.
# Duplicated reads are then marked and removed.
# The final "clean" bam file is sorted and indexed.
# Input:
#   Trimmed fastq.gz file
# Output:
#   Filtered and indexed bam file
#


# move into ChIPseq analysis directory
basedir="/scratch/sclab/021420_ChIP"
cd $basedir
echo "here i am !! $PWD"
mkdir -p alignedbam/filtering_stats

#### process file names ####
# retrieve sample name from the {SLURM array ID}th line of metadata file
IFS=$'\n'
sample_list=($(<${basedir}/CRXHD_ChIP_lookup.txt))
sample=${sample_list[${SLURM_ARRAY_TASK_ID}]}
echo $sample

# input fastq.gz files after adapter trimmer and filtering
r1="${basedir}/trim/${sample}_r001_val_1.fq.gz"
r2="${basedir}/trim/${sample}_r002_val_2.fq.gz"

# intermediate files
outsam="${sample}.aligned.sam"
pairmappedbam="${sample}.pairmapped.bam"
pairmappedindex="${sample}.pairmapped.bam.bai"
nomito="${sample}.nomito.bam"
outstat="filtering_stats/${sample}.picard.metric.txt"
outbam="${sample}.picard.filtered.bam"

# final output bam file and flagstat file
outbamblk="${sample}.sorted.blk.bam"
flagstat="filtering_stats/${sample}.sorted.blk.flagstat.txt"

for elem in $r1 $r2 #$outsam $outstat $outbamblk
do
    echo $elem
done

#### Main ####
# move into the folder where all mapped files will be stored
cd ${basedir}/alignedbam

# first map trimmed and filtered raw sequences to a reference genome (using mm10)
# UCSC mm10 bowtie2 index location: /scratch/sclab/genomes/mm10/Sequence/Bowtie2Index
bowtie2 -t -p $SLURM_CPUS_PER_TASK -X 2000 --very-sensitive \
       -x /scratch/sclab/genomes/mm10/Sequence/Bowtie2Index/genome \
        -1 $r1 -2 $r2 -S $outsam

# processing mapped sam file using samtools
echo "filtering and keeping only uniquely mapped read pairs"
# convert sam to bam format, keep only uniquely mapped and properly paired reads, sort by leftmost coordinates and index
samtools view -b $outsam | samtools view -f 0x2 -q 30 -bu | samtools sort -o $pairmappedbam && samtools index $pairmappedbam
echo "removing chrM reads"
# then remove mitochondria reads and sort
samtools idxstats $pairmappedbam | cut -f 1 | egrep -v chrM | xargs samtools view -bu $pairmappedbam | samtools sort -o $nomito

rm $pairmappedbam
rm $pairmappedindex

echo "marking duplicated reads with picard"
# use picard to filter duplicate reads, very memory intensive, good to limit the amount of RAM allocated with -Xmx value
# read more http://broadinstitute.github.io/picard/faq.html
java -Xmx2g -jar $PICARD MarkDuplicates I=$nomito O=$outbam M=$outstat REMOVE_DUPLICATES=true TMP_DIR=${TMPDIR}

rm $nomito

echo "removing blacklist regions"
# remove blacklist regions of the genome, blacklist bed download from ENCODE, and index the output file for deeptools usage later
bedtools intersect -v -a $outbam -b /scratch/sclab/genomes/blacklist_files/mm10.blacklist.bed.gz > $outbamblk && samtools index $outbamblk
samtools flagstat $outbamblk > $flagstat

rm $outbam

echo "ha! tihs is the end of the script!"

# clear all module loaded
module purge