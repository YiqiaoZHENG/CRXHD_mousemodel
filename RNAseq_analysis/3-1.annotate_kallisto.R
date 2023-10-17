#!/bin/Rscript

print("here i am!!")
print(getwd())

args <- commandArgs(trailingOnly = TRUE)
# tell if extract tpm or non-length normalized counts
if (args=="tpm") {
  print("will extract tpm values")
  count_col <- "tpm"
  col.suffix <- "tpm"
  name.suffix <- "tpm"
} else {
  print("will extract non-length normalized counts values")
  count_col <- "tpm"
  col.suffix <- "count"
  name.suffix <- "estcounts"
}

# Load libraries
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(tidyverse))

#directory variables
basedir="/scratch/sclab/CRX_RNAseq/"
abundancedir="/scratch/sclab/CRX_RNAseq/abundances"
annotationdir="/scratch/sclab/CRX_RNAseq/mapped_tpm"

# helper function to make the most updated TxDb object from Ensemble
# note this database will not have genename
#EnsDb <- makeTxDbFromEnsembl(organism="Mus musculus", release=102, server="ensembldb.ensembl.org")
# save and load
#saveDb(EnsDb, "/scratch/sclab/CRX_RNAseq/annotation/Ensembl/Ens.Mmus.GRCm38.102.db")
#EnsDb <- loadDb("/scratch/sclab/CRX_RNAseq/annotation/Ensembl/Ens.Mmus.GRCm38.102.db", packageName=NA)

#load EnsDb object made from bioconductor AnnotationHub
library(EnsDb.Mmusculus.v102)
Ensdb <- EnsDb.Mmusculus.v102

#### I. Prepare the input abundance files ####
# List all directories containing kallisto output abundance matrix 
samples <- list.files(path = "/scratch/sclab/CRX_RNAseq/abundances/", full.names = F)
# Obtain a vector of all filenames including the path
files <- file.path(abundancedir,samples, "abundance.tsv")
h5files <- file.path(abundancedir,samples, "abundance.h5")
# Since all abundance files have the same name it is useful to have names for each element
names(files) <- samples
names(h5files) <- samples


#### II. Retrieve the abundance from each sample and compile into a single summary dataframe #### 
# helper function to retrieve transcript ids and tpm
select_data <- function(sample){
    #print(paste("reading ",sample,sep=""))
    # read in the abundance.tsv file as data frame
    tx.counts <- read.delim(file.path(abundancedir,sample,"abundance.tsv"), sep="\t", 
                            header=T, stringsAsFactors = FALSE) %>% 
                            dplyr::select("target_id",count_col)
    # retrieve the transcript names as a list of characters
    txids <- as.character(tx.counts[,1])
    require(stringr)
    # remove version number from transcript IDs
    join.txids <- tibble("target_id"=txids,"tx_id"=str_replace(txids, "([.][0-9])", ""))
    tx.counts <- join.txids %>% left_join(tx.counts, by=c("target_id"="target_id")) %>% 
                                dplyr::select("tx_id",count_col)
    # rename tpm column to sample name
    names(tx.counts)[names(tx.counts) == count_col] <- paste(sample,col.suffix,sep=".")
    return(tx.counts)
}

# make a summary transcript x sample count table for all the samples in the abundances folder
if (args=="tpm") {
  # retrieve lengthed normalized tpm values
  all.tpm <- lapply(samples, select_data) %>% `names<-`(samples)
  # merge all count tables in the list by tx_id into a summary dataframe
  all.tx.matrix <- Reduce(function(x, y) merge(x, y, by = "tx_id", all = TRUE), all.tpm)
  print ("take a look at the dimension of transcript tpm matrix")
  dim(all.tx.matrix)
} else {
  # retrieve non-length normalized counts using tximport packages
  txi <- tximport(h5files, type="kallisto", txOut = TRUE, ignoreTxVersion=TRUE, countsFromAbundance="lengthScaledTPM")
  tx.counts <- txi$counts %>% data.frame() %>% `names<-`(samples) %>% rownames_to_column(var = "target_id")
  # remove the version suffix before mapping to gene
  txids <- as.character(tx.counts[,1])
  require(stringr)
  join.txids <- tibble("target_id"=txids,"tx_id"=str_replace(txids, "([.][0-9])", "")) 
  all.tx.matrix <- join.txids %>% left_join(tx.counts, by=c("target_id"="target_id")) %>% dplyr::select(-"target_id")
  print ("take a look at the dimension of transcript count matrix")
  dim(all.tx.matrix)
  # write tximport output count table to file
  #write.table(all.tx.matrix,file.path(annotationdir,paste0("kallisto_countsbytx.ensembl.txt")), sep="\t", quote=F, col.names=T, row.names=F)
}


#### III. Convert transcript level abundance to gene level abundance ####
# retrieve the gene ids and gene names for the set of transcript ids
txtogene <- AnnotationDbi::select(Ensdb,
                                  keys = all.tx.matrix$tx_id, 
                                  columns=c("GENEID","GENENAME"), 
                                  keytype="TXID")
#all_txs <- transcripts(EnsDb, columns=c("tx_id", "tx_name"), filter=NULL, use.names=FALSE)
# combine gene names with abundance data
gene_table <- txtogene %>% left_join(all.tx.matrix, by=c("TXID"="tx_id"))
ordered_table <- gene_table[order(gene_table$GENENAME),]
# handling na entries, either mask to 0 or remove (probably not a good idea for summary matrix)
ordered_table[is.na(ordered_table)] <- 0
#small_table <- gene_table %>% drop_na()

# write to file
write.table(ordered_table,file.path(annotationdir,paste0("kallisto_",name.suffix,"_all.ensembl.txt")), sep="\t", quote=F, col.names=T, row.names=F)
print(paste("results written to mapped_tpm/","kallisto_",name.suffix,"_all.ensembl.txt",sep=""))

# collapse transcript to gene level matrix, use gene id do not use genename
# convert numbers from character to numeric class
ordered_table[4:length(ordered_table)] <- lapply(ordered_table[4:length(ordered_table)], function(x) as.numeric(x))

tpm.names <- names(ordered_table)[4:length(ordered_table)]

sumbygene <- function(sample){
    res <- tapply(ordered_table[,sample], ordered_table[,"GENEID"], sum) %>% data.frame() %>% `names<-`(sample) 
    rownames_to_column(res, var = "GENEID")
}

all.sum <- lapply(tpm.names, sumbygene) %>% `names<-`(tpm.names)
length(all.sum)
names(all.sum)

# merge all sum count by gene tables in the list by GENEID into a summary dataframe
all.sum.matrix <- Reduce(function(x, y) merge(x, y, by = "GENEID", all = TRUE), all.sum)

write.table(all.sum.matrix,file.path(annotationdir,paste0("kallisto_",name.suffix,"_gene.ensembl.txt")), sep="\t", quote=F, col.names=T, row.names=F)
print(paste("results written to mapped_tpm/","kallisto_",name.suffix,"_gene.ensembl.txt",sep=""))

cat("\n")
sessionInfo()
print(paste("any warnings",warnings(),sep=":"))