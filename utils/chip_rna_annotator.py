"""
Author: Yiqiao Zheng
Email: yiqiao.zheng@wustl.edu
"""

import os
import warnings
import re
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels
import statsmodels.api as sm


#### I. Functions for data filtering ####

def chip_fpm_filter(chip_count_df, genotypes=["WT", "E80A", "K88N"], fpm_th=2):
    """
    Given a chip intensity/count matrix, filter regions above fpm cutoff in the specified samples (columns).

    Parameters
    ----------
    chip_count_df: pandas DataFrame
        ChIP binding intensity matrix.

    genotypes: list
        List of samples/columns to be based on for filtering.
    
    fpm_th: float
        CPM cutoff

    Returns
    -------
    chip_count_df: pandas DataFrame
        Filtered ChIP binding intensity matrix.
    """
    chip_count_df = chip_count_df.copy()

    # select by genotypes of interests
    genotypes = "|".join(genotypes)
    check_cols = [x for x in chip_count_df.columns if re.search(genotypes, x)]

    # filter genes based on age and fpm cutoff
    mask = chip_count_df.loc[:, check_cols] >= fpm_th
    chip_count_df =chip_count_df[mask.sum(axis=1) == len(mask.columns)].reset_index(drop=True) # keep rows where all samples pass fpm threshold

    return chip_count_df

def rna_fpm_filter(rna_count_df, rna_lfc_df, genotypes=["wt", "e_het", "e_hom", "k_het", "k_hom", "r_hom"], age="all", fpm_th=2):
    """
    Given a expression count matrix, filter genes above fpm cutoff in the specified samples (columns).
    Genotypes and ages are both required to identify samples.

    Parameters
    ----------
    rna_count_df: pandas DataFrame
        RNA expression count matrix.

    genotypes: list
        List of genotypes to be based on for filtering.

    ages: str
        Ages to be based on for filtering. Default: "all" to filter based on both p10 and p21 RNA data
    
    fpm_th: float
        CPM cutoff

    Returns
    -------
    rna_count_df: pandas DataFrame
        Filtered RNA expression count matrix.
    """
    rna_count_df = rna_count_df.copy()
    rna_lfc_df = rna_lfc_df.copy()

    # select by genotypes of interests
    genotypes = "|".join(genotypes)
    check_cols = [x for x in rna_count_df.columns if re.search(genotypes, x)]

    if age == "p10" or age == "p21":
        check_cols = [x for x in check_cols if re.search(age, x)]

    # filter genes based on age and fpm cutoff
    mask = rna_count_df.loc[:, check_cols] >= fpm_th
    rna_lfc_df = rna_lfc_df[mask.sum(axis=1) == len(mask.columns)].reset_index(drop=True) # keep rows where all samples pass fpm threshold

    return rna_lfc_df


# II. Functions for peakset compiling and gene association

def merge_peakset(chip_peakset, chip_lfc_matrix, fpm_th=2):
    """
    Given a list of genomic regions, retrieve the differential binding matrix.

    Parameters
    ----------
    chip_peakset: one or list of pandas DataFrame
        One or multiple dataframe of ChIP binding intensity matrix with each genomic region identified uniquely by "peak.id"
    
    chip_lfc_matrix: pandas DataFrame
        ChIP differential binding analysis results table
    
    fpm_th: float
        ChIPseq CPM cutoff

    Returns
    -------
    chip_peakset: pandas DataFrame
        Dataframe containing ChIP binding intensity and differential binding statistics
    """
    if type(chip_peakset) == list and len(chip_peakset) > 1: # find union peakset if multiple samples specified
        chip_peakset = pd.concat(chip_peakset, ignore_index=True).copy().drop_duplicates(subset=["peak.id"], keep="first", inplace=False)

    # retrieve chip lfc data
    chip_peakset = pd.merge(chip_peakset, chip_lfc_matrix, how="inner")
    # fpm filtering
    chip_peakset = chip_fpm_filter(chip_peakset, fpm_th=fpm_th)

    return chip_peakset


def is_nearby_peak(chip_data, DE_cutoff=100):
    """Filter all genomic regions within DE_cutoff kb to a TSS."""
    chip_data = chip_data.copy()
    
    # apply distance cutoff
    mask = abs(chip_data.loc[:, "distTSS"]) <= DE_cutoff*1e+3

    # call a gene associated peak within distance cutoff
    nearby_peak = chip_data[mask].reset_index(drop=True)
    
    return nearby_peak


def is_cloest_peak(chip_data, DE_cutoff=100):
    """Identify genomic regions within DE_cutoff kb and is the cloest to a TSS. Only unique genes with a cloest peak are retained."""
    chip_data = chip_data.copy()

    # apply distance cutoff
    nearby_peak = is_nearby_peak(chip_data, DE_cutoff=DE_cutoff)
    
    # order genes by distance to TSS
    nearby_peak = nearby_peak.sort_values(by='distTSS', ascending=True, key=abs, inplace=False)

    # call unique genes with a cloest peak
    cloest_peak = nearby_peak.drop_duplicates('gene', keep='first', inplace=False).reset_index(drop=True)

    return cloest_peak


def is_cloest_gene(chip_data, DE_cutoff=100):
    """Identify genomic regions within DE_cutoff kb and is the cloest to a TSS. Only unique peaks cloest to a gene are retained."""
    chip_data = chip_data.copy()

    # apply distance cutoff
    nearby_peak = is_nearby_peak(chip_data, DE_cutoff=DE_cutoff)
    
    # order genes by distance to TSS
    nearby_peak = nearby_peak.sort_values(by='distTSS', ascending=True, key=abs, inplace=False)
    
    # call unique peaks with a cloest genes
    cloest_gene = nearby_peak.drop_duplicates('peak.id', keep='first', inplace=False).reset_index(drop=True)

    return cloest_gene



def is_depedent_gene(chip_data, rna_data):
    """
    Split a list of genes into two based on whether a ChIP peak is associated not or.

    Parameters
    ----------
    chip_data: pandas DataFrame
        Peak to gene association matrix. Identifier of genes should be under column name "gene".

    rna_data: pandas DataFrame
        Gene name identified matrix. Identifier of genes should be under column name "gene".

    Returns
    -------
    genes_near_peak: pandas DataFrame
        Dataframe of all genes with an associated peak.

    genes_no_peak: pandas DataFrame
        Dataframe of all genes without an associated peak.
    """
    chip_data = chip_data.copy()
    rna_data = rna_data.copy()

    #  now split the rna data into genes with or without an associated peak
    genes_near_peak = pd.merge(rna_data, chip_data, how="right").reset_index(drop=True)
    mask = ~rna_data["gene"].isin(genes_near_peak["gene"]) # mask for genes with no associated peaks
    genes_no_peak = rna_data[mask].reset_index(drop=True) 

    return genes_near_peak, genes_no_peak


# III. functions for chip and rna data annotation

def annot_chip_lfc(data, chip_lfc, chip_fdr, chip_lfc_th=1.0, chip_fdr_th=0.1):
    """
    Given a differential binding analysis results matrix, annotate the state of Gained/Lost/No change of a genomic region. 
    A region with NaN log2FC value will be annotated as "No Peak".

    Parameters
    ----------
    data: pandas DataFrame
        Differential binding analysis results matrix. Need to contain both log2FC and FDR columns.

    chip_lfc: str
        Name of column containing log2FC numbers.

    chip_fdr: str
        Name of column containing fdr numbers.

    chip_lfc_th: float
        Cut-off log2FC value of calling differential binding. Default: 1.0 (2 fold).

    chip_fdr_th: float
        Cut-off FDR value of calling differential binding. Default: 0.1.
    
    Returns
    -------
    data: pandas DataFrame
        The original input dataframe with an added field of "chip_group" for annotation.
    """
    data = data.copy()
    data["chip_group"] = np.nan

    # annotate by chip signal gained/lost
    data[["chip_group"]] = data[["chip_group"]].mask((data[chip_lfc] > chip_lfc_th) & (data[chip_fdr] < chip_fdr_th), 'Gained', inplace=False)
    data[["chip_group"]] = data[["chip_group"]].mask((data[chip_lfc] < -chip_lfc_th) & (data[chip_fdr] < chip_fdr_th), 'Lost', inplace=False)
    data[["chip_group"]] = data[["chip_group"]].mask((abs(data[chip_lfc]) <= chip_lfc_th) | (data[chip_fdr] >= chip_fdr_th) | (data[chip_fdr].isna()), 'No Change', inplace=False)
    data[["chip_group"]] = data[["chip_group"]].mask(data[chip_lfc].isna(), 'No Peak', inplace=False) 

    return data


def annot_rna_lfc(data, rna_lfc, rna_padj, rna_lfc_th=1.5, rna_padj_th=1e-5, simplify=False):
    """
    Given a differential expression analysis results matrix, annotate the state of Gained/Lost/No change of a gene. 
    A gene with NaN log2FC value will be annotated as "NaN"

    Parameters
    ----------
    data: pandas DataFrame
        Differential expression analysis results matrix. Need to contain both log2FC and padj columns.
    
    rna_lfc: str
        Name of column containing log2FC numbers.

    rna_padj: str
        Nmae of column containing padj (adjusted p-value) numbers.

    rna_lfc_th: float
        Cut-off log2FC value of calling differential expression. Default: 1.5. (~3 fold).

    rna_padj_th: float
        Cut-off padj value of calling differential expression. Default: 1e-5.

    simplify: boolean
        If true, only DE or not DE will be annotated. If false, Gained/Lost/No change will be annotated. Default: False.

    Returns
    -------
    data: pandas DataFrame
        The original input dataframe with an added field of "rna_group" for annotation.
    """
    data = data.copy()
    data["rna_group"] = np.nan

    # annotate by gene expression gained/lost
    if simplify: # do not specify gained or lost
        data[["rna_group"]] = data[["rna_group"]].mask((abs(data[rna_lfc]) > rna_lfc_th) & (data[rna_padj] < rna_padj_th), 'Diff. Exp.', inplace=False) 
    else: # specify gained and lost, default
        data[["rna_group"]] = data[["rna_group"]].mask((data[rna_lfc] > rna_lfc_th) & (data[rna_padj] < rna_padj_th), 'Gained', inplace=False)
        data[["rna_group"]] = data[["rna_group"]].mask((data[rna_lfc] < -rna_lfc_th) & (data[rna_padj] < rna_padj_th), 'Lost', inplace=False)
    data[["rna_group"]] = data[["rna_group"]].mask((abs(data[rna_lfc]) <= rna_lfc_th) | (data[rna_padj] >= rna_padj_th), 'No Change', inplace=False)
    data[["rna_group"]] = data[["rna_group"]].mask(data[rna_lfc].isna(), 'NaN', inplace=False)

    return data


def annot_pro_enh(data, chip_lfc, chip_fdr, chip_lfc_th=1.0, chip_fdr_th=0.1, Pr_cutoff=.5, DE_cutoff=100):
    """
    For a list of peak-gene pairs, annotate whether a peak is located at proximal, distal, or not at the regulatory region of a gene based on distance cutoffs.
    Peaks with NaN in the log2FC field or outside of the distal regulatory region distance cutoff will both be annotated as "No CRE".

    Parameters
    ----------
    data: pandas DataFrame
        Dataframe containing peak to gene association and differential binding analysis results.

    chip_lfc: str
        Name of column containing log2FC numbers.

    chip_fdr: str
        Name of column containing fdr numbers.

    chip_lfc_th: float
        Cut-off log2FC value of calling differential binding. Default: 1.0 (2 fold).

    chip_fdr_th: float
        Cut-off FDR value of calling differential binding. Default: 0.1.

    Pr_cutoff: float
        Cut-off distance in kb for calling a peak to be located at proximal regulatory region of a gene. Default: 0.5.

    DE_cutoff: float
        Cut-off distance in kb for calling a peak to be located at distal regulatory region of a gene. Default: 100.

    Returns
    -------
    data: pandas DataFrame
        The original input dataframe with an added field of "CRE_group" of regulatory region annotation.
    """
    data = data.copy()
    data["CRE_group"] = np.nan

    if "chip_group" not in data.columns:
        data = annot_chip_lfc(data, chip_lfc, chip_fdr, chip_lfc_th, chip_fdr_th)
    
    # annotate as Promoter vs Enhancer
    data[["CRE_group"]] = data[["CRE_group"]].mask((abs(data["distTSS"]) <= Pr_cutoff*1e+3) & (data["chip_group"]!="No Peak") , 'Promoter', inplace=False)
    data[["CRE_group"]] = data[["CRE_group"]].mask((abs(data["distTSS"]) > Pr_cutoff*1e+3) & (abs(data["distTSS"]) <= DE_cutoff*1e+3) & (data["chip_group"]!="No Peak"), 'Distal Enhancer', inplace=False)
    data[["CRE_group"]] = data[["CRE_group"]].mask((abs(data["distTSS"]) > DE_cutoff*1e+3) | (data["distTSS"].isna()) | (data["chip_group"]=="No Peak"), 'No CRE', inplace=False)

    return data


def annot_distgroup(data, dist_cutoff=[0, .5, 5, 50, 100]):
    """
    For a list of peak-gene association pairs, assign a peak to different distance interval based on the cutoff values specified.

    Parameters
    ----------
    data: pandas DataFrame
        Dataframe containing peak to gene association and the distance to the TSS stored in "disTSS" column.
    
    dist_cutoff: list
        List of distance cut-off values in kb. Default: [0, .5, 5, 50, 100].

    Returns
    -------
    data: pandas DataFrame
        The original input dataframe with an added field of "dist_group" of distance interval group assigned.    
    """
    data = data.copy()
    data["dist_group"] = np.nan

    # annotate by distance to cloest TSS
    for i in range(len(dist_cutoff)-1):
        cutoff1 = dist_cutoff[i]
        cutoff2 = dist_cutoff[i+1]
        data[["dist_group"]] = data[["dist_group"]].mask((abs(data["distTSS"]) >= cutoff1*1e+3) &
                                                            (abs(data["distTSS"]) < cutoff2*1e+3),
                                                            ' '.join([str(cutoff1), "to", str(cutoff2)+"kb"]), inplace=False)
    data[["dist_group"]] = data[["dist_group"]].mask(abs(data["distTSS"]) >= dist_cutoff[-1]*1e+3, ">"+str(dist_cutoff[-1])+"kb", inplace=False)

    return data


def annotate_rnalfc_by_chip(chip_peakset, chip_lfc_df, peak_togene_df, rna_count_df, rna_lfc_df, chip_geno, rna_geno, age,
                                Pr_cutoff=0.5, DE_cutoff=100, map_method="cloest", dist_cutoff=[0, 0.5, 5, 50, 100], chip_fpm_th=2, rna_fpm_th=5, 
                                chip_lfc_th = 1.0, chip_fdr_th = 0.1, rna_lfc_th = 1.5, rna_padj_th = 1e-5,
                                annot_list=["ProEnh", "rnalfc", "chiplfc"]):
    """
    A wrapper function to compile and annotate a pair of differential binding and differential expression data. Example: CRX ChIP-seq and RNA-seq in CRX HD mutant retinas compared to WT retinas. By default, peaks will be assigned to genes if within 100kb to TSS.

    Parameters
    ----------
    chip_peakset: one or a list of pandas DataFrame
        One or multiple dataframe of ChIP binding intensity matrix with each genomic region identified uniquely by "peak.id"
    
    chip_lfc_df: pandas DataFrame
        ChIP differential binding analysis results matrix. Each peak is uniquely identified by a "peak.id".

    peak_togene_df: pandas DataFrame
        Peak to gene association with distance to TSS information under "distTSS" column, such as that obtained by GREAT analysis. Each peak is uniquely identified by a "peak.id".

    rna_count_df: pandas DataFrame
        Normalized expression matrix with both "GENEID" and "gene" name identifiers.

    rna_lfc_df: pandas DataFrame
        Differential expression alysis results matrix with both "GENEID" and "gene" identifiers.

    chip_geno: str
        Genotype identifier of ChIP data columns.
    
    rna_geno: str
        Genotype identifier of RNA data columns.

    age: str
        Age identifier of RNA data columns.
        
    Pr_cutoff: float
        Cut-off distance in kb for calling a peak to be located at proximal regulatory region of a gene. Default: 0.5.

    DE_cutoff: float
        Cut-off distance in kb for calling a peak to be located at distal regulatory region of a gene. Default: 100.
    
    map_method: str
        Rule to assign a peak to gene. Default is to assign all peaks within a distance cut off to a gene. If "cloest_peak" is specified, each gene is assigned to the cloest peak. If "cloest_gene" is specified, each peak is assigned to the cloest gene.
    
    dist_cutoff: list
        List of distance cut-off values in kb. Default: [0, .5, 5, 50, 100].
    
    chip_fpm_th: float
        ChIPseq CPM cutoff. Default: 2.
    
    rna_fpm_th: float 
        RNAseq CPM cutoff. Default: 5.

    chip_lfc_th: float
        Cut-off log2FC value of calling differential binding. Default: 1.0 (2 fold).

    chip_fdr_th: float
        Cut-off FDR value of calling differential binding. Default: 0.1.
    
    rna_lfc_th: float
        Cut-off log2FC value of calling differential expression. Default: 1.5. (~3 fold).

    rna_padj_th: float
        Cut-off padj value of calling differential expression. Default: 1e-5.

    annot_list: list
        Specify the type of additional annotations to run. If "rna_lfc" is specified, genes will be annotated as Gained/Lost/No Change contained in "rna_group" column. If "chip_lfc" is specified, peaks will be annotated as Gained/Lost/No Change contained in "chip_group" column. If "ProEnh" is specified, peaks will be annotaed as proximal or distal regulatory element contained in "CRE_group" column. If "dist" is specified, peaks will be assigned to different distance interval contained in "dist_group" column.

    Returns
    -------
    data_to_plot: pandas DataFrame
        Dataframe containing the compield ChIP and RNA analysis results with additional fields of annotation.
    """
    ## data preprocessing
    # associate peaks to nearby genes by mapping method specified
    peakset_df = merge_peakset(chip_peakset, chip_lfc_df, fpm_th=chip_fpm_th)
    # map peak to genename
    peakset_df = pd.merge(peak_togene_df, peakset_df, how="inner")
    if map_method == "cloest_peak": # keep unique genes with cloest peak assigned
        peakset_df = is_cloest_peak(peakset_df, DE_cutoff)
    elif map_method == "cloest_gene": # keep unique peaks with cloest gene assigned
        peakset_df = is_cloest_gene(peakset_df, DE_cutoff)
    else: # map_method == "dist" or map_method == "nearby": # keep all gene and peaks within distance cutoff, default
        peakset_df = is_nearby_peak(peakset_df, DE_cutoff)
    # compile chip and rna data
    rna_df = rna_fpm_filter(rna_count_df, rna_lfc_df, age=age, fpm_th=rna_fpm_th)
    genes_near_peak, genes_no_peak = is_depedent_gene(peakset_df, rna_df)

    # concatenate two gene lists and associated data (rbind)
    compiled_df = pd.concat([genes_near_peak, genes_no_peak], axis=0, ignore_index=True)
    
    ## slice out only data specificed by genotype and age
    chip_lfc_col = ".".join(["chip",chip_geno,"lfc"])
    chip_fdr_col = ".".join(["chip",chip_geno,"fdr"])
    rna_lfc_col = ".".join([age,rna_geno,"lfc"])
    rna_padj_col = ".".join([age,rna_geno,"padj"])
    columns_to_select = [chip_lfc_col, chip_fdr_col, rna_lfc_col, rna_padj_col]
    #print("selecting columns " + ", ".join(columns_to_select))
    data_to_plot = compiled_df.loc[:,["gene", "GENEID", "distTSS", "peak.id"] + columns_to_select].copy().reset_index(drop=True)

    ## apply annotations
    if "rnalfc" in annot_list:
        data_to_plot = annot_rna_lfc(data=data_to_plot, rna_lfc=rna_lfc_col, rna_padj=rna_padj_col,  rna_lfc_th=rna_lfc_th, rna_padj_th=rna_padj_th, simplify=False)
        data_to_plot = data_to_plot[data_to_plot.rna_group != "NaN"] # remove any row where there's no rna reads
    
    if "chiplfc" in annot_list:
        data_to_plot = annot_chip_lfc(data=data_to_plot, chip_lfc=chip_lfc_col, chip_fdr=chip_fdr_col, chip_lfc_th=chip_lfc_th, chip_fdr_th=chip_fdr_th)
    
    if "ProEnh" in annot_list:
        data_to_plot = annot_pro_enh(data=data_to_plot, chip_lfc=chip_lfc_col, chip_fdr=chip_fdr_col, chip_lfc_th=chip_lfc_th, chip_fdr_th=chip_fdr_th, Pr_cutoff=Pr_cutoff, DE_cutoff=DE_cutoff)

    if "dist" in annot_list:
        data_to_plot = annot_distgroup(data=data_to_plot, DE_cutoff=dist_cutoff)

    return data_to_plot

# annotate pair of rna data
def annot_paired_rna(rna_df1, rna_df2, diff_th=1.0, simplify=False):
    """
    For a pair of differential expression analysis results tables, identify genes that show concordant changes or sample-specific changes. The input data need to contain the "rna_group" field indicating the whether a gene is called differentially expressed in the sample. Otherwise, run the input dataframe through annot_rna_lfc function fisrt.

    Parameters
    ----------
    rna_df1: pandas DataFrame
        Differential expression analysis results from the first sample.

    rna_df2: pandas DataFrame
        Differential expression analysis results from the second sample.

    diff_th: float
        The difference in log2FC value to call sample-specific changes. Note the two sample should be compared to the same control samples, e.g. WT.

    simplify: boolean
        If True, only sample-specific ("Diff.") or not ("No Diff.") will be annotated. If False, the sample specific group will be further separated into four depending of the direction of changes in sample 1 vs 2. "Gp1": uniquely gained in sample2. "Gp2": uniquely gained in sample1. "Gp3": uniquely lost in sample2. "Gp4": unqiuely lost in sample1.

    Returns
    -------
    data: pandas DataFrame
        Merged differential expression matrix of genes that show differential expression in at least of one the two samples.

    notDE_data: pandas DataFrame
        Merged differential expression matrix of genes that do not show differential expression in both samples.
    """
    # input dataframes should be already be annotated by rna lfc
    rna_df1 = rna_df1.copy().rename(columns=dict(rna_group = "rna_group1"))
    rna_lfc_col1 = [x for x in rna_df1.columns if re.search("lfc", x) and not re.search("chip", x)][0]
    rna_df2 = rna_df2.copy().rename(columns=dict(rna_group = "rna_group2"))
    rna_lfc_col2 = [x for x in rna_df2.columns if re.search("lfc", x) and not re.search("chip", x)][0]

    data = pd.merge(rna_df1, rna_df2, how="inner").reset_index(drop=True)
    # take out genes that are not called D.E. in both samples
    noDE_data = data.loc[(data["rna_group1"] == "No Change") & (data["rna_group2"] == "No Change"), :]
    noDE_data = noDE_data.assign(rna_group="Ctrl")

    # keeps genes that are at least differentially expressed in one sample
    data = data.loc[(data["rna_group1"] != "No Change") | (data["rna_group2"] != "No Change"), :]
    # calculate lfc ratio
    data["lfc.diff"] = data[rna_lfc_col2]-data[rna_lfc_col1]
    data["rna_group"] = np.nan

    # annotate by concordance of change
    if simplify: # do not specify sample
        data[["rna_group"]] = data[["rna_group"]].mask(abs(data["lfc.diff"]) > diff_th , 'Diff.', inplace=False)
    else: # specify gained and lost, default
        data[["rna_group"]] = data[["rna_group"]].mask((data["lfc.diff"] > diff_th) & (data["rna_group2"] == "Gained"), 'Gp1', inplace=False)
        data[["rna_group"]] = data[["rna_group"]].mask((data["lfc.diff"] > diff_th) & (data["rna_group2"] == "No Change"), 'Gp2', inplace=False)
        data[["rna_group"]] = data[["rna_group"]].mask((data["lfc.diff"] < -diff_th) & (data["rna_group2"] == "Lost"), 'Gp3', inplace=False)
        data[["rna_group"]] = data[["rna_group"]].mask((data["lfc.diff"] < -diff_th) & (data["rna_group2"] == "No Change"), 'Gp4', inplace=False)
    
    data[["rna_group"]] = data[["rna_group"]].mask((data["rna_group"].isna()), 'No Diff.', inplace=False)


    keep_cols = [col for col in data.columns if col not in ["rna_group1", "rna_group2", "lfc.diff"]]
    data = data[keep_cols].sort_values(by="rna_group", ascending=False).drop_duplicates('gene', keep='first', inplace=False).reset_index(drop=True) # reorder columns and rows
    noDE_data = noDE_data[keep_cols].sort_values(by="rna_group", ascending=False).drop_duplicates('gene', keep='first', inplace=False).reset_index(drop=True) # reorder columns and rows
    
    return data, noDE_data


### IV. functions for pair-wise rna comparison

def build_contingency_table(data):
    """Helper function to build a 2x2 contigency table comparing peaks annotated as a proximal or distal reulatory element"""
    # prepare the data
    count_tb = data.groupby("CRE_group").rna_group.value_counts().unstack().fillna(0)
    count_tb.loc["with CRE",:] = count_tb.loc["Distal Enhancer",:] + count_tb.loc["Promoter",:]
    count_tb = count_tb.loc[["with CRE", "No CRE"]].T

    return count_tb


def rnaGp_exact_tests(data):
    """
    Wrapper function to run exact tests comparing peaks annotated as a proximal or distal reulatory element.

    Returns
    -------
    stats_summary: pandas DataFrame
        Results summary from two-sided Fisher's, Barnard's, Boschloo's exact tests with multiple correction by "BH" method.
    
    """
    # prepare the data
    summary_tb = build_contingency_table(data)

    v_idx = [idx for idx in summary_tb.index]
    stats_summary = np.empty((len(v_idx),3),dtype=object)

    for i, name in enumerate(v_idx):
        test_table = pd.DataFrame(summary_tb.loc[name, :]).T
        mask = [row for row in summary_tb.index if row != name]
        test_table.loc["None "+name] = summary_tb.loc[mask,:].sum(axis=0)

        # fisher's
        #print("fisher's two sided")
        res = stats.fisher_exact(test_table, alternative='two-sided')
        stats_summary[i,0] = res[1]

        # barnard's
        #print("barnard's two sided")
        res = stats.barnard_exact(test_table, alternative='two-sided')
        stats_summary[i,1] = res.pvalue
        
        # boschloo's
        #print("boschloo's two sided")
        res = stats.boschloo_exact(test_table, alternative='two-sided')
        stats_summary[i,2] = res.pvalue

    stats_summary = pd.DataFrame(data=stats_summary, index=v_idx, columns=["Fisher's p", "Barnard's p", "Boschloo's p"])

    # multiple comparison correction
    stats_summary["Fisher's padj"] = statsmodels.stats.multitest.multipletests(stats_summary["Fisher's p"], alpha=0.05, method='fdr_bh')[1]
    stats_summary["Barnard's padj"] = statsmodels.stats.multitest.multipletests(stats_summary["Barnard's p"], alpha=0.05, method='fdr_bh')[1]
    stats_summary["Boschloo's padj"] = statsmodels.stats.multitest.multipletests(stats_summary["Boschloo's p"], alpha=0.05, method='fdr_bh')[1]
    
    return stats_summary


def odds_ratio_by_Gp(data, uselog=True):
    """Calculate the odds ratio comparing peaks annotated as a proximal or distal reulatory element. If uselog is Ture, return log odds ratio.
    """
    # prepare the data
    summary_tb = build_contingency_table(data)

    v_idx = [idx for idx in summary_tb.index]
    odds_ratio_ser={}
    for i, name in enumerate(v_idx):
        test_table = pd.DataFrame(summary_tb.loc[name, :]).T
        mask = [row for row in summary_tb.index if row != name]
        test_table.loc["None "+name] = summary_tb.loc[mask,:].sum(axis=0)

        # create a statsmodels Table object
        or_table = sm.stats.Table(test_table)
        # retrieve odds ratio
        if uselog:
            odds_ratio = or_table.local_log_oddsratios.iloc[0,0]
        else:
            odds_ratio = or_table.local_oddsratios.iloc[0,0]
        # fill the odds ratio table
        odds_ratio_ser[name] = odds_ratio
    
    odds_ratio_ser=pd.Series(odds_ratio_ser)

    return odds_ratio_ser


def pval_to_star(pval, hide_ns=True):
    """Helper function to convert pvalues from number to asterisk format
    """
    pvalue_thresholds = pd.Series({1e-4: "****", 1e-3: "***", 1e-2: "**", 0.05: "*", 1: "ns"})
    if hide_ns:
         pvalue_thresholds[1] = ""
    for th in pvalue_thresholds.index:
        if pval <= th:
            pval_text = pvalue_thresholds[th]
            break

    return pval_text
