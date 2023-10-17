#!/bin/Rscript

print("here i am!!")
print(getwd())

args <- commandArgs(trailingOnly = TRUE)
cores_to_use <- as.numeric(args)
# tell if running the program in serial or parallel
if (length(args)==0) {
  print("running in serial, if not desired, cancel job and supply a core number")
  BiocParallel::register(BiocParallel::SerialParam()) #default
  run_para<-FALSE
} else {
  print(paste("running in parallel with ",cores_to_use," cores",sep=""))
  BiocParallel::register(BiocParallel::MulticoreParam(cores_to_use))
  run_para<-TRUE
}
# to check bpparam: 
print(BiocParallel::bpparam())
print(paste("cores detected by parallel",parallel::detectCores(logical = FALSE),sep=":")) #this is what the default use of cores by dba.blacklist()

#load DiffBind libraries
suppressPackageStartupMessages(library(DiffBind))
suppressPackageStartupMessages(library(tidyverse))

#directory variables
basedir="/scratch/sclab/021420_ChIP"
outputdir="/scratch/sclab/021420_ChIP/DiffBind/outputs"
beddir="/scratch/sclab/021420_ChIP/DiffBind/diffbed/"

# check metadata spreadsheet
data_sheet<-read.csv(file=file.path(basedir,"CRXHD_DiffBind_meta"))
#print(data_sheet)

#read in a set of peaksets and associated metadata
print("here\'s the metadata")
dbObj <- dba(sampleSheet=data_sheet, dir=basedir)
dbObj

# specify cpu/thread number
dbObj$config$cores<-cores_to_use #this is what i asked in the sbatch script
dbObj$config$RunParallel<-run_para

#### Make consensus peaksets from each set of replicates, then derive master consensus set ####
#print("generating consensus peaksets by factor")
consObj<-dba.peakset(dbObj, consensus=DBA_FACTOR, minOverlap=2) #overlap replicates
consObj
#print("retrieve consensus peakset by factor as GRanges object")
consObj <-  dba(consObj, mask=consObj$masks$Consensus, minOverlap=1) #need to specify minOverlap to include all peaks, default is 2
consObj
factor.consensus <- dba.peakset(consObj, bRetrieve=TRUE)

#### Run differential binding analysis ####
print("counting with consensus peakset by factor")
# if summits=TRUE, summits will be calculated and peaksets unaffected
# if summits>0 then all consensus peaks will be re-centered around a consensus summit with width 2*summit+1, default:200
dbObj<-dba.count(dbObj,peaks=factor.consensus,bUseSummarizeOverlaps=TRUE) 
dbObj
print("normalizing")
dbObj<-dba.normalize(dbObj,method=DBA_ALL_METHODS,normalize=DBA_NORM_LIB,library=DBA_LIBSIZE_FULL) #Use the full library size (total number of reads in BAM/SAM/BED file)
print("contrasting")
dbObj<-dba.contrast(dbObj,categories=DBA_FACTOR,minMembers=2)
print(paste("blacklisting cores = ",cores_to_use,sep=""))
dbObj<-dba.blacklist(dbObj)
print("analyzing")
dbObj<-dba.analyze(dbObj,method=DBA_ALL_METHODS,bRetrieveAnalysis=FALSE)
#note:dba.analyze() retrieves a DBA object with results of analysis added to DBA$contrasts with bRetrieveAnalysis=FALSE
print("analysis done, printing results")
contrast_res<-dba.show(dbObj, bContrasts=T)
print(contrast_res)
print("saving post analysis object")
dba.save(dbObj, file='crx_hdmuts_postanalyze_full_factor_obj', dir=outputdir, pre='dba_', ext='RData', 
         bRemoveAnalysis=FALSE, bRemoveBackground=FALSE,
         bCompress=FALSE)


# retreive peakset post analysis (now the peakset is already centered on summit)
print("retrieving consesus peakset post analysis")
consensus.peaks <- dba.peakset(dbObj, bRetrieve=TRUE)
# make sure to use UCSC conventions i.e. chr1 instead of 1
seqlevelsStyle(consensus.peaks) <- "UCSC"
# save the consensus peakset object
saveRDS(consensus.peaks , file=file.path(outputdir,"hdmuts_chip_peaks.consensus.centered.rds"))

print("the end")

cat("\n")
sessionInfo()
print(paste("any warnings",warnings(),sep=":"))