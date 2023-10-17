#!/bin/Rscript

#remotes::install_github("zeropin/TFCookbook")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(TFCookbook))

# check number and dtype of arguments passed to script
cmd_args <- commandArgs(TRUE)
if (length(cmd_args) < 5) stop("Not enough arguments. Please supply 3 arguments.")

inputFile <-  as.character(cmd_args[1])
library <-  as.character(cmd_args[2])
band <- as.character(cmd_args[3])
MMth <- as.numeric(cmd_args[4])
write_tofile_bool <- as.logical(cmd_args[5])

# infer output dir from input file
outputdir <- file.path(dirname(inputFile), "ewm_mlr")
ewm_id <- str_split(basename(inputFile), pattern = "_", simplify = TRUE)[1]

# create new outputdir if not already exist
dir.create(outputdir, showWarnings = FALSE)

# load RBE data
energy_matrix <- read.csv(file = inputFile, sep = "\t", header = TRUE)

name_to_search <- paste0("avg.", library, ".", band, "ddG")

# check if the specified lib x band exits in the given dataframe
if (!(all_of(name_to_search) %in% colnames(energy_matrix))) {
    stop(paste(name_to_search, "not found, Please check."))
}else{
    # bind the energy model
    energy_matrix %>%
        dplyr::filter(lib == library) %>%
        dplyr::filter(MMcount <= MMth) %>%
        dplyr::select(sequence,name_to_search) %>%
        dplyr::rename(Energy = `name_to_search`, Sequence = `sequence`) %>%
        TFCookbook::buildEnergyModel() -> core_model

    # check mlr model
    summary(core_model)

    # convert to matrix formatted as nucleotide x position
    core_model %>% 
        TFCookbook::getEnergyMatrix() %>%
        as.data.frame() -> core_ewm

    # convert colume names to position
    colnames(core_ewm) <- substring(colnames(core_ewm), 2)

    # format identifier
    ewm_id <- paste(ewm_id,library,band,sep=".")

    # remove columns where all entries are NAs and transpose the matrix to match the ewm reading function
    core_ewm <- core_ewm %>% dplyr::select_if(~all(!is.na(.))) %>% t() %>% as.data.frame() %>% rownames_to_column(var="pos")

    if (!write_tofile_bool) {
        # return identifier
        cat(ewm_id)
        cat("\n")
        # return matrix to stdout as input back to python script
        print(core_ewm, row.names=FALSE)
        
    }else {
        # make position index as column and write to text file
        # first write header with primiary motif identifier start with ">"
        outputFile <- file.path(outputdir,paste0(ewm_id, "_ewm_mlr.txt"))
        # write identifier as header followed by the full ewm
        writeLines(paste0(">", ewm_id), file(outputFile))
        write.table(core_ewm, file=outputFile, append = TRUE, quote = FALSE, sep = " ", row.names=FALSE, col.names=FALSE)
        close(file(outputFile))
        cat(outputFile)
        cat("\n")
    }
}

cat("ha! this is the end of the specseq_mlr.R script!")

# print warnings is any
if (!is.null(warnings())) {message(paste("any warnings",warnings(),sep = ":"))}
