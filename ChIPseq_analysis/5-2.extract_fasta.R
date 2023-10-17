#!/bin/Rscript

print(paste("here i am!!",getwd()))

#load libraries
suppressPackageStartupMessages(library(DiffBind))
suppressPackageStartupMessages(library(tidyverse))
library("BSgenome.Mmusculus.UCSC.mm10")

#directory variables
basedir="/scratch/sclab/021420_ChIP/DiffBind/"
beddir="/scratch/sclab/021420_ChIP/DiffBind/diffbed/"
fastadir="/scratch/sclab/021420_ChIP/DiffBind/difffasta/"

#### I. high confidence peaks for each protein (replicate overlap) ####
## read metadata spreadsheet
data_sheet<-read.csv("CRXHD_DiffBind_meta")
dbObj <- dba(sampleSheet=data_sheet)

# so, let's specify cpu/thread number
dbObj$config$cores<-cores_to_use #this is what i asked in the sbatch script
dbObj$config$RunParallel<-run_para

# subset dba object by each protein and count reads from bam files
print("retrieving DBA objects by protein and count reads")
# specify minOverlap to include only peaks that are presented in both replicates, default 2
# and run dba.count to recenter summits +/-200bp around summit
wt.dba <- dba(dbObj, mask=dbObj$masks$WT) %>%  dba.count(minOverlap=2, bUseSummarizeOverlaps=TRUE) 
e.dba <- dba(dbObj, mask=dbObj$masks$E80A) %>%  dba.count(minOverlap=2, bUseSummarizeOverlaps=TRUE)
k.dba <- dba(dbObj, mask=dbObj$masks$K88N) %>%  dba.count(minOverlap=2, bUseSummarizeOverlaps=TRUE)

# retrieve the corresponding peakset as GRanges object and append to a GRangesList object
print("retrieving peaksets from each protein sample")
wt.olap <- dba.peakset(wt.dba,bRetrieve=TRUE)
e.olap <- dba.peakset(e.dba,bRetrieve=TRUE)
k.olap <- dba.peakset(k.dba,bRetrieve=TRUE)
grl.olap<-GRangesList("wt_olap"=wt.olap, "e80a_olap"=e.olap, "k88n_olap"=k.olap)

# save the peakset files in bed format, will be used in deeptools
for (name in names(grl.olap)){
  out<-as.data.frame(granges(grl.olap[[name]]))
  write.table(out, file=paste(beddir,name,"_peaks",".bed",sep=""), sep="\t", quote=F, row.names=F)
}

# give each region a number, formated as peakset_#, will be stored as the name for each fasta record
for (name in names(grl.olap)){
  names(grl.olap[[name]]) <- paste(name,1:length(grl.olap[[name]]),sep="_")
}

# check the names and number of elements in each peakset
elementNROWS(grl.olap)

# retrieve the DNA sequences based on regions in the GRangesList
seq.olap = getSeq(Mmusculus, grl.olap)
seq.olap

# output the DNA sequences to fasta files, named as protein_olap.fa
for (name in names(seq.olap)){
  file_name <- paste(name,"_200_centered.fa",sep="")
  print(paste("writing",file_name,sep=" "))
  writeXStringSet(seq.olap[[name]], filepath=paste(fastadir,file_name,sep=""), 
                    append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
}

#### II. consensus peakset for all possible CRX binding sites ####
## using bed file format as input
# retreive genomic coordinates dataframe
df <-  read.csv(paste0(beddir, "hdmuts_chip_all_regions.tsv"), sep="\t", header = T) %>% column_to_rownames(var="peak.id")
# convert dataframe to GRanges object
gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE) %>% sort()
# retreive fasta records
seq <- getSeq(genomedb, gr)
# write fasta records to output file
writeXStringSet(seq , filepath=paste0(fastadir,"hdmuts_chip_all_regions.fa"), 
                    append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
    
cat("\n")
sessionInfo()
print(paste("any warnings",warnings(),sep=":"))