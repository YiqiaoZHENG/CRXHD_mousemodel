"""
Author: Yiqiao Zheng
Email: yiqiao.zheng@wustl.edu
"""

import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager
import seaborn as sns
from statannotations.Annotator import Annotator

import chipseq_plot_utils

from utils import specseq_plot_utils

# color palettes for developmental boxplot
select_ages = ["e14.5","e17.5","p0","p3","p7","p10","p14","p21"]
age_pal = pd.Series({age:color for age,color in zip(select_ages, chipseq_plot_utils.palette2hex("Set2"))})


def annotate_test_stats(data, x_column, y_column, morder, figax, test_pairs, test='Mann-Whitney', test_format="star"):
    """
    A wrapper to annotate a categorized box/strip/violin plot. Default it to run Mann-Whitney U test and annotate in asterisks.

    Parameters
    ----------
    data: pandas DataFrame
        Data in the plot to be annotated.
    
    x_column: str
        Column name of x-axis data.

    y_column: str
        Column name of y-axis data.

    morder: str
        Order of  x-axis labels.
    
    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.

    test_pairs: list-like
        List of sample pairs to compare.

    test: str
        Name of statistical test to run. Availabel tests: Mann-Whitney, t-test (independent and paired), Welch's t-test, Levene test, Wilcoxon test, Kruskal-Wallis test. See statannotations package documentation for more information.
    
    test_format: "star" or "simple_text"
        How to format the star notation on the plot.

    Returns
    -------
    fig : Figure handle
    ax : Axes handle
    """
    # retrieve ax to add annotation
    fig, ax = figax

    valid_name = [x.get_text() for x in ax.get_xticklabels()]

    pairs = []
    # check if all pairs are valid
    for pair in test_pairs:
        if pair[0] in valid_name and pair[1] in valid_name:
            pairs.append(pair)
    
    if len(pairs) >= 1:
        annotator = Annotator(ax, pairs, data=data, x=x_column, y=y_column, order=morder)
        # adjust for multiple tests
        # ref: https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
        # test value should be a StatTest instance or one of the following strings: t-test_ind, t-test_welch, t-test_paired, Mann-Whitney, Mann-Whitney-gt, Mann-Whitney-ls, Levene, Wilcoxon, Kruskal.
        annotator.configure(test=test, text_format=test_format, comparisons_correction="fdr_bh", verbose=2, line_width=mpl.rcParams["lines.linewidth"], line_height=0.0)
        annotator.apply_and_annotate()
    
    return fig, ax


def box_by_category(data, x_column, y_column, morder, mpal, annot_pairs, annot_bool=True, xlabel=None, ylabel=None, patch_artist=False, figax=None):
    """
    A wrappter to make boxplot with specific order, box aesthetics, and optional statistical test annotations. 

    Parameters
    ----------
    data: pandas DataFrame
        Long format dataframe to make the boxplot.

    x_column: str
        Column name of x-axis data.

    y_column: str
        Column name of y-axis data.

    morder: list
        Order of x-axis data.
    
    mpal: dictionary
        Dictionary specifying the color of boxes.
    
    annot_pairs: list-like
        List of sample pairs to compare.
    
    annot_bool: boolean
        If True, results from Mann-Whitney U tests with FDR-BH correction will be annotated as asterisks on the plot.

    patch-artist: boolean
        If True, the box will be filled. If False, the box will be open.

    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.
    
    Returns
    -------
    fig : Figure handle
    ax : Axes handle
    """
    if figax:
        fig, ax = figax
    else:
        fig, ax = plt.subplots()

    data = data.copy()

    # check if all keys exist
    v_order = [od for od in morder if od in data[x_column].unique()]
    v_color = pd.Series({k:mpal[k] for k in v_order})

    data = data[data[x_column].isin(v_order)]

    
    if patch_artist:
        # fill box with color
        sns.boxplot(x=x_column, y=y_column, data=data, showfliers = False, ax=ax, order=v_order, palette=v_color)
    else:
        # remove fill color
        boxplot_kwargs = {'boxprops' : {'edgecolor': 'k', 'facecolor': 'none'}}
        # make the plot
        sns.boxplot(x=x_column, y=y_column, data=data, showfliers = False, ax=ax, order=v_order, **boxplot_kwargs)
        
        # color the lines of each box
        for i, artist in enumerate(ax.artists):
            col = v_color[i]
            # This sets the color for the main box
            artist.set_edgecolor(col)
            # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
            # if outlier is removed only 5 Line2D objects
            # automatically get the object number
            n_line2d_objs =  int(len(ax.lines)/len(v_order))
            # Loop over them here, and use the same colour as above
            for j in range(i*n_line2d_objs,i*n_line2d_objs+n_line2d_objs):
                line = ax.lines[j]
                line.set_color(col)
                line.set_mfc(col)
                line.set_mec(col)

    # set axis labels
    if xlabel:
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel("")
    if ylabel:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel("")
    
    # add annotation
    if annot_bool and len(annot_pairs)>0 and len(v_order)>1:
        fig, ax = annotate_test_stats(data, x_column, y_column, v_order, figax=(fig, ax), test_pairs=annot_pairs)

    return fig, ax


def strip_by_category(data, x_column, y_column, morder,  mpal, add_mean=False, showmean_th=15, mean_wd=.25, markersize=8, xlabel=None, ylabel=None, ylimits=None, annot_pairs=None, annot_bool=True, figax=None):
    """
    A wrappter to make stripplot with specific order, marker aesthetics, and optional statistical test annotations. 

    Parameters
    ----------
    data: pandas DataFrame
        Long format dataframe to make the boxplot.

    x_column: str
        Column name of x-axis data.

    y_column: str
        Column name of y-axis data.

    morder: list
        Order of x-axis data.
    
    mpal: dictionary
        Dictionary specifying the color of boxes.

    add_mean: boolean
        If True, a short bar showing the mean of data will be shown.
    
    showmean_th: int
        Only show mean value bar for categories with datapoint over this cut-off. 
    
    mean_wd: float
        Length of the short bar showing the mean.
    
    markersize: float
        Strip plot marker size. 

    ylimits: list
        Top and bottom limits of the y-aixs.
    
    annot_pairs: list-like
        List of sample pairs to compare.
    
    annot_bool: boolean
        If True, results from Mann-Whitney U tests with FDR-BH correction will be annotated as asterisks on the plot.

    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.
    
    Returns
    -------
    fig : Figure handle
    ax : Axes handle
    """
    if figax:
        fig, ax = figax
    else:
        fig, ax = plt.subplots()

    data = data.copy()

    idx_counts = data[x_column].value_counts()
    # check if all keys exist and keep only groups with at least five datapoints
    v_order = [od for od in morder if od in data[x_column].unique() and idx_counts[od] >=5]
    v_color = pd.Series({k:mpal[k] for k in v_order})
    data = data[data[x_column].isin(v_order)]

    # check if all keys exist
    #v_order = [od for od in morder if od in data[x_column].unique()]
    #v_color = pd.Series({k:mpal[k] for k in v_order})
    
    sns.stripplot(x=x_column, y=y_column, data=data, ax=ax, order=v_order, palette=v_color, size=markersize, zorder=1)
    
    # add mean line
    if add_mean:
        if not mean_wd:
            mean_wd=.25
        hline_paras = dict(linestyle='-', color = 'black', alpha=0.6)
        means = data.groupby(x_column, sort=False)[y_column].mean()
        _ = [ax.hlines(means[od], i-mean_wd, i+mean_wd, **hline_paras, zorder=2) for i,od in enumerate(v_order) if idx_counts[od] >= showmean_th]

    # set axis labels
    if xlabel:
        ax.set(xlabel=xlabel)
    if ylabel:
        ax.set(ylabel=ylabel)

    if ylimits:
        ax.set_ylim(ylimits)

    # add annotation if all conditions met
    if annot_bool and len(annot_pairs)>0 and len(v_order)>1:
        idx_counts = data[x_column].value_counts()
        # only annotate pairs with at least xx datapoint available
        pairs = [pair for pair in annot_pairs if pair[0] in idx_counts.index and pair[1] in idx_counts.index and idx_counts[pair[0]]>=showmean_th and idx_counts[pair[1]]>=showmean_th]
        fig, ax = annotate_test_stats(data, x_column, y_column, v_order, figax=(fig, ax), test_pairs=pairs)

    return fig, ax


def violin_by_category(data, x_column, y_column, morder, mpal, xlabel=None, ylabel=None, figax=None):
    """
    A wrappter to make violin plot with specific order.

    Parameters
    ----------
    data: pandas DataFrame
        Long format dataframe to make the boxplot.

    x_column: str
        Column name of x-axis data.

    y_column: str
        Column name of y-axis data.

    morder: list
        Order of x-axis data.
    
    mpal: dictionary
        Dictionary specifying the color of boxes.

    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.
    
    Returns
    -------
    fig : Figure handle
    ax : Axes handle
    """    

    if figax:
        fig, ax = figax
    else:
        fig, ax = plt.subplots()

    data = data.copy()

    # check if all keys exist
    v_order = [od for od in morder if od in data[x_column].unique()]
    v_color = pd.Series({k:mpal[k] for k in v_order})

    data = data[data[x_column].isin(v_order)]

    sns.violinplot(x=x_column, y=y_column, data=data, cut=0, ax=ax, order= v_order, palette=v_color)

    # set axis labels
    if xlabel:
        ax.set(xlabel=xlabel)
    if ylabel:
        ax.set(ylabel=ylabel)

    return fig, ax


# developmental line plot
def single_line_by_category(data, x_column, y_column, xlabel, ylabel, color, figax=None):
    """
    Draw a single line showing the relationship between x and y.

    Parameters
    ----------
    data: pandas DataFrame
        Long format dataframe.

    x_column: str
        Column name of identifier variables to plot on the x-axis.

    y_column: str
        Column name of measured variables to plot on the y-axis.

    xlabel: str
        X-axis label name.

    ylabel: str
        Y-axis label name.

    color: str
        Color of line to be drawn.

    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.
    
    Returns
    -------
    fig : Figure handle
    ax : Axes handle
    """
    if figax:
        fig, ax = figax
    else:
        fig, ax = plt.subplots()
    
    ax =  sns.lineplot(data=data, x=x_column, y=y_column, ax=ax, color=color)
    ax.set(xlabel=xlabel, ylabel=ylabel)

    return fig, ax


def dev_line_by_rnaGp(data, morder, mpal, ages=["p3","p7","p10","p14","p21"], xlabel=None, ylabel=None, uniform_y=False, figax=None):
    """
    Given a categoried gene name identified RNA expression data in wide format, draw a set of line plots for each category specified in the "rna_group" 
    
    Parameters
    ----------
    data: pandas DataFrame
        RNA expression matrix. Row z-score normalized is prefered.

    morder: list
        Order of rna_group category to draw the line plot.
    
    mpal: dictionary
        Dictionary specifying the color of lines.
    
    ages: list
        Column name list of measured varivables to plot.

    xlabel: str
        X-axis label name.

    ylabel: str
        Y-axis label name.

    uniform_y: boolean
        If True, all line plots will use the same y-axis limits.

    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.
    
    Returns
    -------
    fig : Figure handle
    ax_list: List of axes handles
    """

    data = data.copy()
    
    # retrieve unique groups to plot
    if figax:
        fig, ax_list = figax
    else:
        fig, ax_list = specseq_plot_utils.setup_multiplot(len(morder), len(morder), sharex=False, sharey=False) # note, do not sharey, it will use the smallest y limits
    ax_list = ax_list.flatten()
    
    # format axes labels:
    if xlabel:
        xlabel = xlabel
    else:
        xlabel = ""
    if ylabel:
        ylabel = ylabel
    else:
        ylabel = ""

    plot_parameters = dict(
        xlabel = xlabel,
        ylabel = ylabel,
        x_column="age",
        y_column="norm.counts"
    )

    for gp, ax in zip(morder, ax_list):
        # reshape data
        data_to_plot = data[data.rna_group == gp].drop(columns="rna_group").reset_index(drop=True)
        # convert wide to long dataframe
        data_to_plot = pd.melt(data_to_plot, id_vars=["gene"], value_vars=[age for age in ages if age in data.columns], var_name="age", value_name="norm.counts")
        fig, ax =  single_line_by_category(data=data_to_plot, color=mpal[gp], **plot_parameters, figax=(fig,ax))
        ax.set(title=gp)

    # adjust for y axes, make them the same
    yaxis_limits = list(zip(*[ax.get_ylim() for ax in ax_list.flat]))
    bottom, top = (min(yaxis_limits[0]), max(yaxis_limits[1]))

    for ax in fig.get_axes():
        ss = ax.get_subplotspec()
        if uniform_y:
            ax.set_ylim(bottom, top)
            if not ss.is_first_col():
                ax.yaxis.set_visible(False) # hide the entire axis, labels, ticks, title
        else:
            if not ss.is_first_col():
                ax.set_ylabel("") # high only the title

    return fig, ax_list