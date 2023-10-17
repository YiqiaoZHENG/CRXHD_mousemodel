#!/bin/bash

#SBATCH --job-name=yq_deeptools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/sclab/021420_ChIP/slurm_outputs/slurm_array_deeptools/plotMx_%j.out

# clear all module loaded in case of dependency conflicts
module purge
#load new modules required, using default
module load deeptools/3.0.0


#### Usage of this script ####
# There are three deeptools functions contained in this script.
# These three functions should be used sequentially/individually.
# 1.bamCoverage: calculate coverage at 10kb bin windows from input bam file and output to bigWig format.
#               The output bigwig files are used to generate average coverage files.
# 2.bigwigCompare: calculate average coverage from two bigWig files and output to an avg.bigWig file
#               The output avg.bigWig files are used to generate intensity heatmaps.
# 3.computeMatrix and plotHeatmap: generate binding intensity matrix from bigwig and bed files if supplied.
#               Coverage heatmaps are then generated using the intensity matrix.
#


# move into ChIPseq analysis directory
basedir="/scratch/sclab/021420_ChIP"
cd $basedir
echo "here i am !! $PWD"

#directory variables
bamdir="${basedir}/alignedbam/"
bigwigdir="${basedir}/bigwigs/"
beddir="${basedir}/ChIPQC/diffbed/"
matrixdir="${basedir}/deeptools/matrix/"
plotdir="${basedir}/deeptools/plots/"

#### bamCoverage to calculate coverage from bam file and output to bigwig file ####
# retrieve sample name from the {SLURM array ID}th line of metadata file
IFS=$'\n'
sample_list=($(<${basedir}/CRXHD_ChIP_lookup.txt))
sample=${sample_list[${SLURM_ARRAY_TASK_ID}]}
echo $sample

echo ${SLURM_ARRAY_TASK_ID} 
echo $sample

# input bam file name for bamCoverage
# output bigwig file name
outbamblk=${bamdir}${SAMPLE}.sorted.blk.bam
coverage=${bigwigdir}${SAMPLE}.sorted.blk.bigWig

# convert filtered, sorted, indexed bam file to bigwig
bamCoverage -b $outbamblk -o $coverage --binSize 10 -e --normalizeUsing CPM -p ${SLURM_NTASKS_PER_NODE}

#### end of bamCoverage ####


#### bigwigCompare to calcualte the average coverage from replicates ####
#retrieve sample name from the metadata file
IFS=$'\n'
sample_list=($(<$basedir/021420_ChIP_lookup.txt))
rep1name="${sample_list[4*${SLURM_ARRAY_TASK_ID}]}"
rep2name="${sample_list[4*${SLURM_ARRAY_TASK_ID}+2]}"
outputname="${bigwigdir}$(echo ${rep1name} | sed 's/1_input//')input_avg.bigWig"

# input bigwig files
rep1coverage="${bigwigdir}${rep1name}.sorted.blk.bigWig"
rep2coverage="${bigwigdir}${rep2name}.sorted.blk.bigWig"

echo "calculating average coverage for ${rep1coverage} and ${rep2coverage}"
echo "will write output to ${outputname}"

# calculate mean coverage from replicate bigwig files
bigwigCompare -b1 ${rep1coverage} -b2 ${rep2coverage} --operation mean --binSize 10 -p ${SLURM_NTASKS_PER_NODE} -o ${outputname} -of bigwig

#### end of bigwigCompare ####


#### computeMatrix and plotHeatmap ####
#standardize output names
name_suffix=hdmuts_intensity_e80a_k88n
output_gz="${matrixdir}${name_suffix}.gz"
output_svg="${plotdir}${name_suffix}.svg"
ref_point=summit

## process and select input bigwig files
# Method1. reads automatically, might not be in order, so need to use --smartLabels if not specific requirement for the format
# the variables in the array will be a full path string /path/filename.suffix, but deeptools will remove the path
#score_list=($(find $bigwigdir -maxdepth 1 -type f  -name "*ip*" | sort ))
#score_files=("${score_list[@]##$bigwigdir}")
#score_names=("${score_files[@]%.sorted.blk.bigWig}")

# Method2. force sample order by storing the names in an array
# crx chip samples with replicates
#score_names=(wt1 wt2 e1 e2 k1 k2 r1 r2)
# force same order as lables array, prepanding path and appending suffix
#score_files=("${score_names[@]/%/_ip.sorted.blk.bigWig}")
#score_list=("${score_files[@]/#//scratch/sclab/021420_ChIP/bigwigs/}")
#score_names=(wt.1 wt.2 e80a.1 e80a.2 k88n.1 k88n.2 r90w.1 r90w.2)

# crx chip sample average coverages (as in final manuscript)
score_names=(wt e k r)
score_files=("${score_names[@]/%/_avg.bigWig}")
score_list=("${score_files[@]/#//scratch/sclab/021420_ChIP/bigwigs/}")
score_names=(wt.chip e80a.chip k88n.chip r90w.chip)

echo "all input bigwig files in the current run"
IFS=$'\n'; echo "${score_list[*]}" 

## process and select input region files
# two clusters of peaks based on intensity clustering
# directly retreive bed file names from a matedata text file
region_files=($(sed -n '1,2p' /scratch/sclab/021420_ChIP/deeptools_bed_meta.txt))
region_list=("${region_files[@]/#/${beddir}}")
region_names=("${region_files[@]%_regions.bed}")
region_names=("${region_names[@]##contrast_peakset/}")

echo "all regions in the current run"
IFS=$'\n'; echo "${region_files[*]}"

echo "running computeMatrix, output will be written to ${output_gz}"
# generate the intermediate matrix of scores per genome region for plotHeatmap
computeMatrix reference-point -S ${score_list[@]} \
-R ${region_list[@]} \
-o $output_gz \
--samplesLabel ${score_names[@]} --referencePoint center -b 2000 -a 2000 -bs 10 \
--sortRegions descend --sortUsing mean --sortUsingSamples 1 3 \
--missingDataAsZero --skipZeros \
-p ${SLURM_NTASKS_PER_NODE}

intensity_hm_colors=('#FFFFFF, #08306B' '#FFFFFF, #08306B' '#FFFFFF, #08306B' '#FFFFFF, #08306B')

echo "making deeptools heatmap, will be written to ${output_svg}"
# make the heatmap from matrix (svg format for publication)
plotHeatmap -m $output_gz -o ${plotdir}${name_suffix}.svg --dpi 600 \
--refPointLabel $ref_point \
--xAxisLabel 'dist. to summit (bp)' --labelRotation 0 \
--regionsLabel ${region_names[@]} --samplesLabel ${score_names[@]} \
--sortRegions descend --sortUsing mean --sortUsingSamples 1 3 \
--colorList ${intensity_hm_colors[@]} --colorNumber 444 \
--whatToShow 'heatmap and colorbar' \
--boxAroundHeatmaps yes --legendLocation best \
--zMin 0 --zMax 1.8 \
--plotFileFormat svg

# make the heatmap from matrix (png format for quick check on server)
plotHeatmap -m $output_gz -o ${plotdir}${name_suffix}.png --dpi 600 \
--refPointLabel $ref_point \
--xAxisLabel 'dist. to summit (bp)' --labelRotation 0 \
--regionsLabel ${region_names[@]} --samplesLabel ${score_names[@]} \
--sortRegions descend --sortUsing mean --sortUsingSamples 1 3 \
--colorList ${intensity_hm_colors[@]} --colorNumber 444 \
--whatToShow 'heatmap and colorbar' \
--boxAroundHeatmaps yes --legendLocation best \
--zMin 0 --zMax 1.8 \
--plotFileFormat png

echo "ha! tihs is the end of the script!"

# clear all module loaded
module purge