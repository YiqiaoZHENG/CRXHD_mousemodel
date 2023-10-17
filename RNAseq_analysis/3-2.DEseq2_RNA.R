#!/bin/Rscript

print("here i am!!")
print(getwd())

args <- commandArgs(trailingOnly = TRUE)
count_file <- as.character(args[2])
cores_to_use <- as.numeric(args[1])

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
#print(paste("cores detected by parallel",parallel::detectCores(logical = FALSE),sep=":")) #this is what the default use of cores by dba.blacklist()
#print(parallel::detectCores(logical = FALSE))

# Load libraries for DEseq2 analysis
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(DEGreport))
suppressPackageStartupMessages(library(pheatmap))
#load EnsDb object made from bioconductor AnnotationHub
library(EnsDb.Mmusculus.v102)
Ensdb <- EnsDb.Mmusculus.v102

#directory variables
basedir="/scratch/sclab/CRX_RNAseq/deseq2/"
annotationdir="/scratch/sclab/CRX_RNAseq/mapped_tpm"
count_file <- "kallisto_estcounts_gene.ensembl.txt"
outputdir="/scratch/sclab/CRX_RNAseq/deseq2/output/hdmuts/"

print(paste0("at the end of task, workspace image will be saved at ",outputdir,"kallisto_hdmuts_DEseq2.RData"))

print(paste0("input count file: ",count_file))
print(paste0("type of counts: ",str_split(count_file,"_",simplify=TRUE)[,2]))

#### I.retrieve count matrix as dataframe ####
all.matrix <- read.table(file.path(annotationdir,count_file), 
                        sep="\t", header=T, row.names="GENEID")
countData <- all.matrix %>% dplyr::mutate_if(is.numeric, round)
rm(all.matrix)
rm(all.genes)
# check dimemsion of count matrix
print("take a quick look at the dimension of our count matrix")
dim(countData)

countData <- countData %>% rownames_to_column("GENEID")
# find genenames for all geneids and write to file
txtogene <- AnnotationDbi::select(Ensdb,
                                  keys = countData$GENEID, 
                                  columns="GENENAME", 
                                  keytype="GENEID")
write.table(txtogene,file.path(annotationdir,paste0("kallisto_hdmuts.allgenes.ensembl.txt")), sep="\t", quote=F, col.names=T, row.names=F)

#### II.prepare p10 and p21 input count matrix separately ####
p10_countData <- countData %>% select(contains("p10")) %>% select(-contains(c("crxko_hom","d_hom")))
p21_countData <- countData %>% select(contains("p21")) %>% select(-contains(c("crxko_hom","d_hom")))

# retain only genes that have at least 1 count in all samples
p10_countData <- p10_countData %>% filter(rowSums(p10_countData>1)==ncol(p10_countData)) #[1] 15772    18
p21_countData <- p21_countData %>% filter(rowSums(p21_countData>1)==ncol(p21_countData)) #[1] 16646    18

# generate metadata table
make_metaTable <- function(count_table) {
    # automatically process experimental design for colData based on the column names
    samples <- str_split(colnames(count_table),"\\.",simplify=TRUE)[,1]
    splitnames <- str_split(samples,"_",simplify=TRUE)
    mut.all <- c("e_het","e_hom", "k_het","k_hom","r_hom")
    genotype <- factor(paste(splitnames[,1],splitnames[,2],sep="_"),levels=c("wt_hom",mut.all))
    age <- factor(splitnames[,3])
    replicate <- factor(splitnames[,4])
    geno_age <- factor(paste(genotype, age,sep="."))
    metaTable <- data.frame(samples,genotype,age,replicate,geno_age,row.names="samples")
    #take a look at the metaTable
    print(metaTable)
    return(metaTable)
}
p10_metaTable <- make_metaTable(p10_countData)
p21_metaTable <- make_metaTable(p21_countData)

# initialize DEseq2 objects
p10.deObj <- DESeqDataSetFromMatrix(p10_countData, colData=p10_metaTable, design= ~ 1)
p21.deObj <- DESeqDataSetFromMatrix(p21_countData, colData=p21_metaTable, design= ~ 1)

# filter and keep only genes with at least 1 fpm count in at least three samples (which mostly likely come from the same genotype)
filtered_p10_countData <- p10_countData %>% filter(rowSums(fpm(p10.deObj, robust=TRUE)>1)>3)
filtered_p21_countData <- p21_countData %>% filter(rowSums(fpm(p21.deObj, robust=TRUE)>1)>3)

rm(p10.deObj)
rm(p21.deObj)

#### III. Run differential expression analysis ####
# update the DEseq2 objects with fpm filtered count matrix
p10.deObj <- DESeqDataSetFromMatrix(filtered_p10_countData, colData=p10_metaTable, design= ~ genotype)
p21.deObj <- DESeqDataSetFromMatrix(filtered_p21_countData, colData=p21_metaTable, design= ~ genotype)
# and run DE analysis by genotype
p10.deObj <- DESeq(p10.deObj)
p21.deObj <- DESeq(p21.deObj)

# retrieve DE analysis results by age and genotype contrast
shrink_by_age <- function(sub.deObj){
    # retrieve results table for given contrast
    # results returned are filtered by adjusted p value below a given FDR cutoff, alpha, default 0.1
    log.res <- sapply(mut.all, function(x)
                               results(sub.deObj,contrast=c("genotype",x,"wt_hom"),
                               pAdjustMethod="BH", alpha=0.1,
                               format='DataFrame', parallel=run_para))
    # perform Wald testing on the shrunken log2 foldchanges between specific conditions
    print("performing pair wise Wald testing")
    all.res <- mapply(function(x,y)
                               lfcShrink(sub.deObj, coef=x, res=y, type="apeglm", parallel=run_para),
                               paste("genotype",mut.all,"vs","wt_hom",sep="_"), log.res, SIMPLIFY = FALSE)
    return(all.res)
}
run_para=FALSE
p10.DEanalysis <- shrink_by_age(p10.deObj)
p21.DEanalysis <- shrink_by_age(p21.deObj)

# convert list of tables to tibbles for easy access
tbs_to_tibble <- function(tables){
  all.tibble <- lapply(tables, function(x) x %>% data.frame() %>% rownames_to_column(var = "GENEID") %>% as_tibble())
  return(all.tibble)
}
p10.DEanalysis <- tbs_to_tibble(p10.DEanalysis)
p21.DEanalysis <- tbs_to_tibble(p21.DEanalysis)

# save the differential analysis objects in rds format
saveRDS(p10.DEanalysis, paste0(outputdir,"kallisto_hdmuts_p10deAnalysis.rds"))
saveRDS(p21.DEanalysis, paste0(outputdir,"kallisto_hdmuts_p21deAnalysis.rds"))

# save individual DE analysis results in table format and write to a .tsv file
save_raw_res <- function(res,prefix){
  saving <- sapply(names(res), function (x)
                               write.table(file=paste0(outputdir,"contrasts/",prefix,".",str_replace(x, paste0(str_split(x, "_",simplify=TRUE)[,1],"_"),""),".tsv"),
                               res[[x]], sep="\t", quote=F, row.names=F))
}
save_raw_res(p10.DEanalysis,"p10")
save_raw_res(p21.DEanalysis,"p21")

#### IV. Compile DE analysis for each age into a summary table ####
#for each age, gather DE analysis results for all genotypes into a summary table 
DE_matrix_age <- function(age, count.matrix, sig.DE.genes){
    lfc.list <- lapply(names(sig.DE.genes), function(x) sig.DE.genes[[x]] %>% 
                        dplyr::select(GENEID,log2FoldChange,padj) %>% 
                        as_tibble() %>% 
                        `names<-`(c("GENEID",paste0(age,".",str_split(x, pattern="_", simplify=TRUE)[2],str_split(x, pattern="_", simplify=TRUE)[3],".lfc"),paste0(age,".",str_split(x, pattern="_", simplify=TRUE)[2],str_split(x, pattern="_", simplify=TRUE)[3],".padj"))))
    all.lfc <- Reduce(function(x, y) merge(x, y, by = "GENEID", all = TRUE), lfc.list)
    DE.matrix <- count.matrix %>% data.frame() %>% rownames_to_column(var="GENEID") %>%
                 dplyr::filter(GENEID %in% all.lfc$GENEID)

    # attach genenames using geneid as mapping
    all.lfc <- all.lfc %>% left_join(txtogene, by=c("GENEID"="GENEID")) 
    DE.matrix <- DE.matrix %>% left_join(txtogene, by=c("GENEID"="GENEID"))

    return(list("log2FoldChange"=all.lfc ,"normalized.counts"=DE.matrix))
}
p10.compiled.lfc <- DE_matrix_age(filtered_p10_countData,p10.DEanalysis)
p21.compiled.lfc <- DE_matrix_age(filtered_p21_countData,p21.DEanalysis)

# write compiled tables to .tsv files
write.table(p10.compiled.lfc[["log2FoldChange"]], file=paste0(outputdir,"compiled_matrix/", "p10_log2FoldChange.tsv"), sep="\t", quote=F, row.names=F)
write.table(p10.compiled.lfc[["normalized.counts"]], file=paste0(outputdir,"compiled_matrix/", "p10_normalizedCounts.tsv"), sep="\t", quote=F, row.names=F)
write.table(p21.compiled.lfc[["log2FoldChange"]], file=paste0(outputdir,"compiled_matrix/", "p21_log2FoldChange.tsv"), sep="\t", quote=F, row.names=F)
write.table(p21.compiled.lfc[["normalized.counts"]], file=paste0(outputdir,"compiled_matrix/", "p21_normalizedCounts.tsv"), sep="\t", quote=F, row.names=F)

# Save the entire workspace and individual DESeq2 object
save.image(file = paste0(outputdir,"kallisto_hdmuts_DEseq2.RData"))
saveRDS(p10.deObj, paste0(outputdir,"kallisto_hdmuts_p10deObj.rds"))
saveRDS(p21.deObj, paste0(outputdir,"kallisto_hdmuts_p21deObj.rds"))

cat("\n")
sessionInfo()
print(paste("any warnings",warnings(),sep=":"))