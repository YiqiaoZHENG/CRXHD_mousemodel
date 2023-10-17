#!/bin/bash

# Usage: Given FASTQ files, parse raw fastq.gz files and then count motifs in each file
if [ "$#" -ne 5 ]
then
    echo "Usage: fastqCountMotifs.sh specseq_meta.txt libraries consensus path/to/this/script/ ID"
    echo "\t- specseq_meta.txt is formatted as [sample name],[raw fq directory],[output directory prefix],[library path]."
    echo "\t- libraries is a comman deliminated list of specseq libraries in the samples, and used for pattern searching."
    echo "\t- consensus is a comman deliminated list of consensus sequences for each specseq library, separated by comma, and used in calculating relative binding energy."
    exit 1
fi

# parse bash command line arguments
fastqToSample=$1
IFS=',' read -a libraries <<< $2
IFS=',' read -a consensus <<< $3
scriptPath=$4

SLURM_ARRAY_TASK_ID=$5

# If SLURM_ARRAY_TASK_ID is not set, that means either (1) this script is being run locally or (2) there is no job array. In this case, it is assumed that there is only one sample to process
if [ -z "${SLURM_ARRAY_TASK_ID}" ]
then
    SLURM_ARRAY_TASK_ID=0
fi

# parse arguments in specseq_meta.txt
IFS=$'\n'
sample_list=($(<${fastqToSample}))
IFS=',' read -a array <<< "${sample_list[SLURM_ARRAY_TASK_ID]}"
sample_name=${array[0]}
rawfq_dir=${array[1]}
prefix=${array[2]}
libpath=${array[3]}
#short_name=${array[4]} # for samples where the names in fq.gz are abbreviation

# initialize a new specseq analysis directory if not already existed
mkdir -p "${prefix}"
mkdir -p "${prefix}/raw_fq"
mkdir -p "${prefix}/counts"
mkdir -p "${prefix}/RBE"
mkdir -p "${prefix}/RBE/avgRBE"
mkdir -p "${prefix}/plots"

# find all R1 read fastq.gz files for ${sample_name} in ${rawfq_dir} and parse
# then extract and pre-count testing sequences based on constant flanking sequences
fq_to_parse=($(find "${rawfq_dir}" -maxdepth  1 -type f  -name "*${sample_name}*_R1_*"))
for fq in "${fq_to_parse[@]}"; do
    # retrieve band name and preparse fq.gz files
    band_name=$(echo "${fq}" | sed -n 's/.*Chen_\(.*\)_SIC_.*/\1/p')
    echo "- preparsing R1 fastq.gz file for ${sample_name} ${band_name} band"
    seqtk seq -a ${fq} > "${prefix}"/raw_fq/"${sample_name}${band_name: -1}".fastq && echo "- successfully written to ${prefix}/raw_fq/${band_name}.fastq"
    # updatae the band name and point to the preparsed fq.gz file
    band_name="${sample_name}${band_name: -1}"
    # retrieve the total number of reads in the raw fastq
    raw_reads="$(grep ">" "${prefix}/raw_fq/${band_name}.fastq" | wc -l)"
    echo "total raw reads for ${band_name}: ${raw_reads}"

    # extract motifs bewteen constant flanking sequences and count occurence
    echo "- extracting and counting motifs"
    # Pull lines from the FASTQ that contain at least 50 consecutive letters in the set [ACGTN] ...
    grep -e "[ACGTN]\{50,\}" "${prefix}"/raw_fq/"${band_name}".fastq |
    # ...then extract the letters between constant flanking sequences, i.e. testing motifs ...
    sed -n 's/.*TACCGAG\(.*\)AGATCAG.*/\1/p' |
    # ...sort and count the motifs...
    sort | uniq -c |
    # ...reformat to get rid of leading spacer and change the delimiter to a tab
    sed -e 's/^ *//;s/ /\t/' > "${prefix}"/counts/"${band_name}"_counts.tmp && echo "- successfully written to ${prefix}/counts/${band_name}_counts.tmp"

done

# process all pre-count files and calculate relative binidng energy
python3 "${scriptPath}"/calculate_ratios.py --sample "${sample_name}" --lib "${libraries[@]}" --ref "${consensus[@]}" --libpath "${libpath}" --prefix "${prefix}"

# Cleanups
rm "${prefix}"/counts/"${sample_name}"*.tmp

echo "ha! this is end of the script!!"