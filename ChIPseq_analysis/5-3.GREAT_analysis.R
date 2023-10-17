#!/bin/Rscript

print(paste("here i am!!",getwd()))

#load libraries
suppressPackageStartupMessages(library(DiffBind))
suppressPackageStartupMessages(library(rGREAT))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape))


#directory variables
basedir="/scratch/sclab/021420_ChIP"
greatdir="/scratch/sclab/021420_ChIP/DiffBind/great/"
beddir="/scratch/sclab/021420_ChIP/DiffBind/diffbed/"

#### run GREAT analysis for all genomic regions in the master consensus peakset ####
cons.peaks <- read.csv(paste0(beddir, "hdmuts_chip_all_regions.tsv"), sep="\t", header = T) %>% mutate(strand="*")
cons.peaks <- makeGRangesFromDataFrame(cons.peaks, keep.extra.columns=FALSE) 
# submit for great analysis, smaller the request_interval (in s) faster the processing speed
print(paste("submitting all peaks in CRX consensus peakset for great analysis"))
job <- submitGreatJob(gr=cons.peaks, species="mm10", adv_upstream=5, adv_downstream=1, version="4", request_interval=30)
# retrieve GO enrichment table
#print("retrieveing GO term enrichment tables")
#gotb <- getEnrichmentTables(job, category = "GO")
# retrieve GRange objects containing the gene-region associations without making plots
print("retrieving gene-region associations")
res <- plotRegionGeneAssociationGraphs(job, request_interval=30, plot=FALSE)
write.table(res, file=paste(beddir,"great_subset/","crx_chip_all_regions_great_gene.txt",sep=""), sep="\t", quote=F, row.names=F, col.names=T)
# save the job object in rds format
print(job)
saveRDS(job, file=.path(greatdir,"crx_hdmuts_consensus_peaks.great.rds"))

cat("\n")
sessionInfo()file
print(paste("any warnings",warnings(),sep=":"))