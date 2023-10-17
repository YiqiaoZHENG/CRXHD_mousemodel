#!/bin/bash

#SBATCH --job-name=yq_meme
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/sclab/021420_ChIP/slurm_outputs/slurm_meme/memechip_%A_%a.out
#SBATCH --array=0-2

# clear all module loaded in case of dependency conflicts
module purge
#load new modules required
module load meme/5.0.4
module load ghostscript/9.19

#### Usage of this script ####
# This script takes in a list of fasta sequences of ChIPseq peaks and runs de novo motif searching on these sequences.
# For the best performance, I have centered all peaks on their summit and takes the +/-200 bp region surrounding the summit.
# So all the input fasta sequences have a uniform length of 401bp.
# Input:
#   Fasta file containing the ChIPseq peak sequences.
# Output:
#   Folder of all MEME-ChIP output.
#

# move into ChIPseq analysis directory
basedir="/scratch/sclab/021420_ChIP/"
cd $basedir
echo "here i am !! $PWD"
mkdir -p MEME

#directory variables
fastadir="${basedir}DiffBind/difffasta/"
memedir="${basedir}MEME/"
refdb="/scratch/ref/meme/motif_databases/"
db1=${refdb}JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme
db2=${refdb}MOUSE/uniprobe_mouse.meme
db3=${refdb}EUKARYOTE/jolma2013.meme
db4=${refdb}EUKARYOTE/homeodomain.meme
#db5=${refdb}MOUSE/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme

# retrieve the names of all fasta files to be submitted to meme
IFS=$'\n'
sample_list=($(<${basedir}meme_meta.txt))
sample=${sample_list[${SLURM_ARRAY_TASK_ID}]}
query_fa="${fastadir}${sample}_centered.fa"
new_dir="${memedir}olap/${sample}"

echo "submitting sample ${sample} to MEME-CHIP"
echo "query ${query_fa}"
echo "meme outputs will be written to directory ${new_dir}" 
# run de novo motif searching with meme-chip
# only search motifs of length 6 to 14bp, this encompasses monomeric and dimeric homeodomain motifs
meme-chip -oc $new_dir -db $db1 -db $db2 -db $db3 -db $db4 -dna -order 1 -ccut 0 \
-meme-mod anr -meme-nmotifs 20 -meme-minw 6 -meme-maxw 14 \
-meme-p ${SLURM_NTASKS_PER_NODE} \
-spamo-skip -fimo-skip  \
$query_fa

echo "ha! tihs is the end of the script!"

# clear all module loaded
module purge