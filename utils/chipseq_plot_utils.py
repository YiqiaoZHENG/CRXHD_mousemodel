"""
Author: Yiqiao Zheng
Email: yiqiao.zheng@wustl.edu
"""

import os
import warnings

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import mannwhitneyu, normaltest
import fastcluster

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import matplotlib.font_manager
import seaborn as sns
from statannotations.Annotator import Annotator

# I. Miscellaneous functions for data processing and plot handling

# heatmap colors
div_heat_colors = mpl.colors.LinearSegmentedColormap.from_list(
                        "yq_divergent", [(0, "#CF8B03"), (0.5, "#FFFFFF"), (1, "#08306B")])
single_heat_colors = mpl.colors.LinearSegmentedColormap.from_list(
                        "yq_single", [(0, "#D5DCE6"), (1, "#08306B")])

def palette2hex(mpl_pal):
    """Given a matplotlib palette name, convert it to a list of hex color codes."""
    # check if requested palette is valid
    if mpl_pal in cm._cmap_registry.keys():
        # from a matplotlib palette, retireve color in hex
        cmap = cm.get_cmap(mpl_pal)
        cmap_hex = []
        for i in range(cmap.N):
            rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
            cmap_hex.append(matplotlib.colors.rgb2hex(rgb))
        return cmap_hex
    else:
        warnings.warn(mpl_pal + " is not a matplotlib palette")


def chip_intensity_heatmap(data,  hm_title=None, hm_xlabel=None, cb_title=None, cmap=div_heat_colors, paramdict=None):
    """
    A wrapper function to plot a matrix dataset as a hierarchically-clustered heatmap. By default, row z-score will be calculated and used for clustering.

    Parameters
    ----------
    data: pandas DataFrame
        Binding intensity data for clustering. Cannot contain NAs.

    hm_title: str or None
        If specified, will be set to title name of heatmap
    
    hm_xlabel: str or None
        If specified, will be set to the x-axis label name.

    cb_title: str or None
        If specified, will be set to the color bar name.

    cmap: Matplotlib colormap
        Color map used to make the heatmap.

    paramdict: dict or None
        All other keyword arguments are passed to seaborn.heatmap().

    Returns
    -------
    cg: Seaborn ClusterGrid instance.
    """
    # default seaborn clustermap parameters
    default_params = {
    # data, normalization, clustering
    'z_score': 0, # z normalized by row
    'metric': "euclidean", # plot with euclidean distance
    'method': "complete", # linkage method to use for calculating clusters: Farthest Point Algorithm/Voor Hees Algorithm
    'row_cluster': True,
    'col_cluster': False,

    # dendrogram, colorbar, colormap
    'cbar_pos': (1, .3, .03, .4), #tuple of (left, bottom, width, height),
    'robust': True,
    'center': 0.0,
    'cmap': cmap,
    'cbar_kws': { # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html
            'orientation': "vertical",
            'ticklocation': "right",
            'pad': 0.08 # padding between the color bar and new image axes
            },

    # figsize, axes ticks, axes lables
    'yticklabels': False,
    'figsize': mpl.rcParams["figure.figsize"]
    }

    # add additional parameters if specified
    if paramdict:
        plot_params = default_params | paramdict # for py3.9+
        # plot_params = {**default_params, **paramdict} # py3.5+
    else:
        plot_params = default_params

    # cluster data and make heatmap
    cg = sns.clustermap(data=data.copy(), **plot_params)
    # omit the dendrogram
    cg.ax_row_dendrogram.set_visible(False)

    # get the heatmap axes
    ax = cg.ax_heatmap
    # heatmap title
    if hm_title:
        ax.set_title(hm_title)
    # heatmap axes labels
    if hm_xlabel:
        ax.set_xlabel(hm_xlabel, fontsize=mpl.rcParams["axes.titlesize"])
    ax.get_yaxis().set_visible(False) # ax.set_ylabel("")
    # heatmap axis tick labels
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=mpl.rcParams["axes.titlesize"])

    # name the color bar
    if cb_title:
        cg.ax_cbar.set_title(cb_title, fontsize=mpl.rcParams["axes.labelsize"], pad=8.0)

    return cg


def nclust_heatmap(data, clust_col, nclust=1, cmap=div_heat_colors, paramdict=None):
    """
    A helper function to make a set of heatmaps for each cluster based on hierarchical clustering of the full matrix. By default, row z-score will be calculated and used for clustering.

    Parameters
    ----------
    data: pandas DataFrame
        Binding intensity data for clustering. Cannot contain NAs.

    clust_col: list
        List of columns to be used for hierarchical clustering.

    nclus: int
        Number of clusters to cut.

    cmap: Matplotlib colormap
        Color map used to make the heatmap.
    
    paramdict: dict or None
        All other keyword arguments are passed to seaborn.heatmap().

    Returns
    -------
    data: pandas DataFrame
        Original input data matrix with an added field of cluster number.
    cg_list: List of seaborn ClusterGrid instances for all clusters.
    """
    # first check if dataframe has all clust_col
    for col in clust_col:
        if col not in data.columns:
            warnings.warn("Column name " + col + " not found in dataframe")

    # default seaborn clustermap parameters
    default_params = {
    # data, normalization, clustering
    'z_score': 0, # z normalized by row
    'metric': "euclidean", # plot with euclidean distance
    'method': "complete", # linkage method to use for calculating clusters: Farthest Point Algorithm/Voor Hees Algorithm
    'row_cluster': True,
    'col_cluster': False,

    # dendrogram, colorbar, colormap
    'cbar_pos': (1, .3, .03, .4), #tuple of (left, bottom, width, height),
    'robust': True,
    'center': 0.0,
    'cmap': cmap,
    'cbar_kws': { # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html
            'orientation': "vertical",
            'ticklocation': "right",
            'pad': 0.08 # padding between the color bar and new image axes
            },

    # figsize, axes ticks, axes lables
    'yticklabels': False,
    'figsize': mpl.rcParams["figure.figsize"]
    }

    # add additional parameters if specified
    if paramdict:
        plot_params = default_params | paramdict # for py3.9+
    else:
        plot_params = default_params
    
    data = data.copy()
    # calculate linkage matrix
    row_linkage = fastcluster.linkage(data.loc[:,clust_col], method='complete', metric='euclidean', preserve_input='True')

    # retrieve flat clusters
    row_cluster = scipy.cluster.hierarchy.fcluster(row_linkage, nclust, criterion='maxclust')
    # assign cluster number to peak
    data["row_cluster"] = row_cluster

    cg_list = []

    # make a small heatmap for each cluster defined
    for i in range(1,nclust+1,1):
        data_to_plot = data.loc[data.row_cluster == i,clust_col].copy()
        print("number of peaks in cluster " + str(i) + " " + str(len(data_to_plot)))
        # make heatmap
        cg = sns.clustermap(data = data_to_plot, **plot_params)
        # omit the dendrogram
        cg.ax_row_dendrogram.set_visible(False)
        # name heatmap with cluster number
        cg.ax_heatmap.set_title("cluster " + str(i))
        # remove heatmap yaxis labels
        cg.ax_heatmap.get_yaxis().set_visible(False)

        cg_list.append(cg)
    
    return data, cg_list


def parse_clustered_peakset(df, cluster_col, prefix):
    """
    Helper function to parse the clustered chipseq matrix by cluster and save in .bed format.
    
    Parameters
    ----------
    df: pandas DataFrame
        Binding intensity data with cluster definition.
    
    cluster_col: str
        Column name of cluster definition.

    prefix: str
        Directory where the bed files will be saved to.

    Returns
    -------
    None
    """
    df = df.copy()
    for name in df[cluster_col].unique():
        small_df = df.loc[df.row_cluster == name, ["seqnames", "start", "end"]]
        small_df.to_csv(os.path.join(prefix, name+"_regions.bed"), sep="\t", header=False, index=False)


