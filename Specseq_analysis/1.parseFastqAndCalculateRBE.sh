#!bin/bash
"""
Author: Yiqiao Zheng
Email: yiqiao.zheng@wustl.edu
"""
# Usage: Wrapper script for running the spec-seq count processing pipeline
if [ "$#" -ne 2 ]
then
    echo "Usage: parseFastqAndCalculateRBE.sh path/to/current/script path/to/specseq/analysis/directory"
    echo "\t- should supply a full path as {basedir} to store all processed files"
    echo "\t- example uages: bash parseFastqAndCalculateRBE.sh /mnt/v/yqzheng/qiaoer/VSCode_yiqiao/SPEC-SEQ/scripts /mnt/v/yqzheng/qiaoer/VSCode_yiqiao/SPEC-SEQ/091021_HD_mono_pool9"
    exit 1
fi

scriptpath=$1
basedir=$2

mkdir -p "${basedir}"
mkdir -p "${basedir}/logs"

# process all samples in the metadata file
sample_list=($(<"${basedir}/specseq_meta.txt"))
for (( ID=0; ID<("${#sample_list[@]}"); ID++)); do
    IFS=',' read -a array <<< "${sample_list[ID]}"
    sample_name=${array[0]}
    # retrieve protein name only
    protein_name=$(echo "${sample_name}" | sed 's/_.*//')
    # use different consensus based on mutation
    case ${protein_name} in
    k88q|k88n) 
        echo "Using Q50 consensus: ${sample_name}"
        bash "${scriptpath}"/fastqCountMotifs.sh "${basedir}"/specseq_meta.txt M,Mrev TAATTA TAATTA "${scriptpath}" "${ID}" #> ${basedir}/logs/${sample_name}_log.txt

    ;;
    *)
        echo "Using K50 consensus: ${sample_name}"
        bash "${scriptpath}"/fastqCountMotifs.sh "${basedir}"/specseq_meta.txt M,Mrev TAATCC,GGATTA "${scriptpath}" "${ID}" #> ${basedir}/logs/${sample_name}_log.txt

    ;;
    esac
done

echo "ha! this is the end of the script!"

