"""
Author: Yiqiao Zheng
Email: yiqiao.zheng@wustl.edu
"""

import os
from datetime import datetime
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import logomaker

from utils import specseq_ewm_utils

# Define the two colormaps: classic and colorblind safe
dna_classic = {
    # as in most of the pwm figures
    "A": "#003200",
    "C": "#000064",
    "G": "#644100",
    "T": "#640000",
    "filler": "#979797"
}
dna_safe = {
    # from the Stormo lab
    "A": "#0E927B",
    "C": "#59A9D8",
    "G": "#DC9514",
    "T": "#1A1A1A",
    "filler": "#979797"
}
# Define library design (update as new libraries added)
lib_designs = {
    # for automatically adding filter sequences
    "M":"TAANNN",
    "Mrev":"NNNTTA",
}


def set_manuscript_params():
    """
    Set the matplotlib rcParams to values for manuscript-size figures.

    """
    mpl.rcParams["figure.figsize"] = (2, 2)
    mpl.rcParams["figure.titlesize"] = 7
    mpl.rcParams["axes.titlesize"] = 7
    mpl.rcParams["axes.labelsize"] = 7
    mpl.rcParams["axes.titlepad"] = 4
    mpl.rcParams["xtick.labelsize"] = 6
    mpl.rcParams["ytick.labelsize"] = 6
    mpl.rcParams["axes.labelpad"] = 2
    mpl.rcParams["legend.fontsize"] = 6
    mpl.rcParams["image.cmap"] = "viridis"
    mpl.rcParams["lines.markersize"] = 1.0
    mpl.rcParams["lines.linewidth"] = 1.0
    mpl.rcParams["font.size"] = 5
    mpl.rcParams["savefig.dpi"] = 300
    mpl.rcParams['font.family'] = "sans-serif"
    mpl.rcParams['font.sans-serif'] = "Arial"
    mpl.rcParams['mathtext.default'] = "regular"



def set_color(values):
    """
    A wrapper for converting numbers into colors. Given a number between 0 and 1, convert it to the corresponding color in the color scheme.
    
    """
    my_cmap = mpl.cm.get_cmap()
    return my_cmap(values)


def add_letter(ax, x, y, letter):
    """
    Add a letter to label an axes as a panel of a larger figure.

    Parameters
    ----------
    ax : Axes object
        The panel to add the letter to.

    x : int
        x coordinate of the right side of the letter, in ax.transAxes coordinates

    y : int
        y coordinate of the top side of the letter, in ax.transAxes coordinates

    letter : str
        The letter to add

    Returns
    -------
    Text
        The created Text instance
    """
    return ax.text(x, y, letter, fontsize=mpl.rcParams["axes.labelsize"], fontweight="bold", 
                    ha="right", va="top", transform=ax.transAxes)


def save_fig(fig, prefix, tight_layout=True, timestamp=True, tight_pad=1.08):
    """
    Save a figure as a PNG and an SVG.
    
    """
    if tight_layout:
        fig.tight_layout(pad=tight_pad)
    if timestamp:
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        fig.text(0, 0, now, transform=fig.transFigure)

    # create directory if not already exists
    savedir = os.path.split(prefix)[0]
    if not os.path.exists(savedir):
        os.mkdir(savedir)

    fig.savefig(f"{prefix}.svg", bbox_inches="tight")
    fig.savefig(f"{prefix}.png", bbox_inches="tight")
    

def setup_multiplot(n_plots, n_cols=2, sharex=True, sharey=True, big_dimensions=True):
    """
    Setup a multiplot and hide any superfluous axes that may result.

    Parameters
    ----------
    n_plots : int
        Number of subplots to make
    n_cols : int
        Number of columns in the multiplot. Number of rows is inferred.
    sharex : bool
        Indicate if the x-axis should be shared.
    sharey : bool
        Indicate if the y-axis should be shared.
    big_dimensions : bool
        If True, then the size of the multiplot is the default figure size multiplied by the number of rows/columns.
        If False, then the entire figure is the default figure size.

    Returns
    -------
    fig : figure handle
    ax_list : list-like
        The list returned by plt.subplots(), but any superfluous axes are removed and replaced by None
    """
    n_rows = int(np.ceil(n_plots / n_cols))
    row_size, col_size = mpl.rcParams["figure.figsize"]

    if big_dimensions:
        row_size *= n_rows
        col_size *= n_cols

    fig, ax_list = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=(col_size, row_size), sharex=sharex, sharey=sharey)

    # The index corresponding to n_plots is the first subplot to be hidden
    for i in range(ax_list.size):
        coords = np.unravel_index(i, ax_list.shape)
        ax = ax_list[coords]
        if i >= n_plots:
            ax.remove()
            ax_list[coords] = None

    return fig, ax_list


def make_correlation_scatter(sample_tuple, sample_labels, lr=True, colors="black", marker='o', mrsize=8, 
                                xticks=None, yticks=None, stepsize=1, figax=None):
    """
    Wrapper function to make a scatter plot of relative binding energy and fit data with linear regression.


    Parameters
    ----------
    sample_tuple: tuple
        Pair of RBE data
    
    sample_labels: list
        Labels to be used as axis labels

    lr: boolean
        Boolean indicating whether draw the linear fit line. Default: True
    
    colors: str
        Marker color
    
    mrsize: float
        Marker size

    xticks : list-like
        If specified, set the x axis ticks to these values.
        
    yticks: list-like
        If specified, set the y axis ticks to these values.

    stepsize: float
        Stepsize of the major axis ticks. Default: 1
        Overrided when xticks and yticks are specified.

    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.

    Returns
    -------
    fig : Figure handle
    ax : Axes handle
    text: tuple
        Tuple of PCC, SCC, slope, intercept, mean_diffRBE values
    """

    # prepare data to plot
    x=sample_tuple[sample_labels[0]]
    y=sample_tuple[sample_labels[1]]
    
    # calculate mean dddG
    sample_tuple["diffRBE"] = abs(sample_tuple[sample_labels[0]]-sample_tuple[sample_labels[1]])
    mean_diffRBE=np.mean(sample_tuple.diffRBE)
  
    # regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    x_lin_reg = np.linspace(min(x), max(x), 100)
    y_lin_reg = slope*x_lin_reg+intercept

    # calculate Pearson correlation coefficient
    pearson_corr,_ = stats.pearsonr(x, y)
    # calculate Spearman rank correlation coefficient
    spearman_corr,_ = stats.spearmanr(x, y)

    # format disply text
    text = (f'r = {pearson_corr:.3f}', f'\u03C1 = {spearman_corr:.3f}', slope, intercept, mean_diffRBE)

    if figax:
        fig, ax = figax
    else:
        fig, ax = plt.subplots()

    ax.scatter(x, y, marker = marker, s=mrsize, facecolors='none', edgecolors = colors)
    ax.set_xlabel(str.replace(sample_labels[0].upper(),"_"," rep"))
    ax.set_ylabel(str.replace(sample_labels[1].upper(),"_"," rep"))

    # define tick positions, make x and y axes the same
    axis_lower_limit = round(min(min(x), min(y))-stepsize/2)
    axis_upper_limit = round(max(max(x), max(y))+stepsize/2)

    if xticks is not None:
        ax.set_xticks(xticks)
    else:
        ax.set_xticks(np.arange(axis_lower_limit, axis_upper_limit, step=stepsize))

    if yticks is not None:
        ax.set_yticks(yticks)
    else: 
        ax.set_yticks(np.arange(axis_lower_limit, axis_upper_limit, step=stepsize))

    # draw regression line
    if lr:
        ax.plot(x_lin_reg, y_lin_reg, '--', color = 'orange', alpha=0.8)
    # draw x=y line
    ax.plot(x_lin_reg, x_lin_reg, '--', color = 'grey', alpha=0.8)

    # draw axis lines
    ax.hlines(y=0, xmin=axis_lower_limit, xmax=axis_upper_limit, linestyle='--', color = 'grey', alpha=0.4)
    ax.vlines(x=0, ymin=axis_lower_limit, ymax=axis_upper_limit, linestyle='--', color = 'grey', alpha=0.4)

    # add plot title
    ax.set_title("Relative Binding Energy [kT]")

    # clean up the plot
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
                
    # force equal aspect
    ax.set_aspect("equal")

    return fig, ax, text


def create_filler_ewm(filler_seq, positions, height):
    """
    Helper function to generate a dummy ewm for filler sequences.

    Parameters
    ----------
    filler_seq: str
        DNA sequence to be filled in the full length motif

    positions: str
        Positions of the filler sequence in the full length motif

    height: float
        Dummy energy value of all bases in the filler sequences

    Returns
    -------
    filler_ewm: dataFrame
        Dummy ewm of filler sequences
    """
    filler_ewm = {}
    ewm_map={"A": [1,0,0,0],
             "C": [0,1,0,0],
             "G": [0,0,1,0],
             "T": [0,0,0,1]}
    for seq,pos in zip(filler_seq, positions):
        filler_ewm[pos] = [x*height for x in ewm_map[seq]]
    filler_ewm = pd.DataFrame(filler_ewm).T
    filler_ewm.columns = ['A', 'C', 'G', 'T']
    filler_ewm.index = filler_ewm.index.astype('int64')
    filler_ewm.index.name = "pos"

    return filler_ewm


def make_ewm_logo(energyMatrix, sequence_design, add_filler=True, colors=dna_safe, x_label=None, title=None, figax=None):
    """
    Make EWM logo based on energy model and sequence design.

    Parameters
    ----------
    energyMatrix: dataFrame
        EWM dataFrame, positions need to be 1 offset.
    
    sequence_design: str
        Full length library design. Variable regions are masked by N.
    
    add_filler: boolean
        Boolean indicating whether to plot the constant regions of the Spec-seq library

    colors: dictionary
        Color scheme used to make the EWM logo. Default: dna_safe

    x_label: str
        If specified, x-axis label will be set.
    
    title: str
        If specified, figure title will be set.

    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.
    
    Returns
    -------
    fig : Figure handle
    ax : Axes handle
    ewm_logo: Logomaker object
    """
    # normalize the energyMatrix
    energyMatrix = specseq_ewm_utils.normalize_ewm(energyMatrix.copy())
    # plotting -E for visualization purpose
    energyMatrix = -energyMatrix

    style_filler = False
    # check parameters before adding filler sequences
    if add_filler:
        # first retrieve filler sequence and positions according to seqeunce design
        N_pos = np.array([pos+1 for pos,nuc in enumerate(sequence_design) if nuc == "N"])
        nonN_pos = [pos+1 for pos,nuc in enumerate(sequence_design) if nuc != "N"]
        ewm_pos = np.array(energyMatrix.index)
        # check if given ewm matches that of the sequence design
        if bool(N_pos.all() == ewm_pos.all()):
            filler_seq = "".join(sequence_design[i-1] for i in nonN_pos)
            filler_pos = nonN_pos
            # add filler (constant) sequences for display purposes
            half_max_height = 0.4*max(list(energyMatrix[energyMatrix>0].sum(1)))
            filler_ewm = create_filler_ewm(filler_seq, filler_pos, half_max_height)
            energyMatrix = pd.concat([filler_ewm, energyMatrix]).sort_index()
            style_filler = True
    
    if figax:
        fig, ax = figax
    else:
        fig, ax = plt.subplots()

    # create Logo object
    ewm_logo = logomaker.Logo(energyMatrix,
                                ax=ax,
                                #shade_below=.5,
                                #fade_below=.5,
                                font_name="Arial",
                                color_scheme=colors,
                                flip_below = False)

    if style_filler:
        # styling the filler letters
        for letter,pos in zip(filler_seq, filler_pos):
            ewm_logo.style_single_glyph(pos, letter, ax=ax, color=dna_safe["filler"], ceiling=half_max_height)

    # style using Logo methods
    ewm_logo.style_spines(visible=False)
    ewm_logo.style_spines(spines=['left', 'bottom'], visible=True)
    ewm_logo.style_xticks(rotation=0, fmt='%d', anchor=0)

    # style using Axes methods
    ewm_logo.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    #ewm_logo.ax.set_ylabel("$-\Delta \Delta G$\n(kcal/mol)", fontsize=mpl.rcParams["ytick.labelsize"])
    ewm_logo.ax.set_ylabel("-Energy\n(kcal/mol)", fontsize=mpl.rcParams["ytick.labelsize"])
    ewm_logo.ax.tick_params(axis='x')
    ewm_logo.ax.tick_params(axis='y')
    ewm_logo.ax.xaxis.set_ticks_position('none')
    ewm_logo.ax.xaxis.set_tick_params(pad=-1)

    # optional styling
    if x_label:
        ewm_logo.ax.set_xlabel("position", fontsize=mpl.rcParams["xtick.labelsize"])
    if title:
        ewm_logo.ax.set_title(title)

    return fig, ax, ewm_logo


def make_pwm_logo(probabilityMatrix, sequence_design, add_filler=True, colors=dna_safe, x_label=None, title=None, figax=None):
    """
    Make PWM logo based on energy model and sequence design.

    Parameters
    ----------
    probabilityMatrix: dataFrame
        PWM dataFrame, positions need to be 1 offset.
    
    sequence_design: str
        Full length library design. Variable regions are masked by N.
    
    add_filler: boolean
        Boolean indicating whether to plot the constant regions of the Spec-seq library

    colors: dictionary
        Color scheme used to make the PWM logo. Default: dna_safe

    x_label: str
        If specified, x-axis label will be set.
    
    title: str
        If specified, figure title will be set.

    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.
    
    Returns
    -------
    fig : Figure handle
    ax : Axes handle
    pwm_logo: Logomaker object
    """
    probabilityMatrix = probabilityMatrix.copy()
    # convert from probability to information matrix
    pwm_information = logomaker.transform_matrix(probabilityMatrix, from_type="probability", to_type="information")

    style_filler = False
    # check parameters before adding filler sequences
    if add_filler:
        # first retrieve filler sequence and positions according to seqeunce design
        N_pos = np.array([pos+1 for pos,nuc in enumerate(sequence_design) if nuc == "N"])
        nonN_pos = [pos+1 for pos,nuc in enumerate(sequence_design) if nuc != "N"]
        ewm_pos = np.array(probabilityMatrix.index)
        # check if given ewm matches that of the sequence design
        if bool(N_pos.all() == ewm_pos.all()):
            filler_seq = "".join(sequence_design[i-1] for i in nonN_pos)
            filler_pos = nonN_pos
            # add filler (constant) sequences for display purposes
            half_max_height = 0.4*max(list(pwm_information[pwm_information>0].sum(1)))
            filler_ewm = create_filler_ewm(filler_seq, filler_pos, half_max_height)
            pwm_information = pd.concat([filler_ewm, pwm_information]).sort_index()
            style_filler = True

    if figax:
        fig, ax = figax
    else:
        fig, ax = plt.subplots()

    # create Logo object
    pwm_logo = logomaker.Logo(pwm_information,
                                ax=ax,
                                #shade_below=.5,
                                #fade_below=.5,
                                font_name="Arial",
                                color_scheme=colors)
    
    if style_filler:
        # styling the filler letters
        for letter,pos in zip(filler_seq, filler_pos):
            pwm_logo.style_single_glyph(pos, letter, ax=ax, color=dna_safe["filler"], ceiling=half_max_height)

    # style using Logo methods
    pwm_logo.style_spines(visible=False)
    pwm_logo.style_spines(spines=['left', 'bottom'], visible=True)
    pwm_logo.style_xticks(rotation=0, fmt='%d', anchor=0)

    # style using Axes methods
    pwm_logo.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    pwm_logo.ax.set_ylabel("Information\nContent (bits)", labelpad=3)
    pwm_logo.ax.tick_params(axis='x')
    pwm_logo.ax.tick_params(axis='y')
    pwm_logo.ax.xaxis.set_ticks_position('none')
    pwm_logo.ax.xaxis.set_tick_params(pad=-1)

    # optional styling
    if x_label:
        pwm_logo.ax.set_xlabel("position", labelpad=3)
    if title:
        pwm_logo.ax.set_title(title)

    return fig, ax, pwm_logo


def probBound_RBE_lines(sequence_design, ewm, consensus, mu_values, figax=None):
    """
    Helper function that given a ewm and sequence design, plot predicted occupancy v.s. several different mu values
    
    Parameters
    ----------
    sequence_design: str
        Full length library design. Variable regions are masked by N.

    ewm: dataFrame
        Energy weight matrix
    
    consensus: str
        Consensus DNA sequence
    
    mu_values: list-like
        List of chemical energies to make prediction on

    figax : (figure, axes) or None
        If specified, make the plot in the provided axes. Otherwise, generate a new axes.
    
    Returns
    -------
    fig : Figure handle
    ax : Axes handle
    """
    # first retireve predicted binding energy matrix
    predicted_energy_df = specseq_ewm_utils.predict_bindingEnergy(sequence_design=sequence_design, ewm=ewm)
    # offset by consensus binding energy
    predicted_energy_df[["pred.ddG"]] = predicted_energy_df[["pred.ddG"]] - predicted_energy_df.loc[consensus,"pred.ddG"]
    predicted_energy_df = predicted_energy_df.sort_values(by=["pred.ddG"])
    relative_energy = np.array(predicted_energy_df["pred.ddG"])
    mu_colors = set_color(np.arange(len(mu_values)) / len(mu_values))

    if figax:
        fig, ax = figax
    else:
        fig, ax = plt.subplots(figsize=(4, 4))
        
    for mu, color in zip(mu_values, mu_colors):
        occupancy_probability = 1.0 / (1.0 + np.exp(relative_energy - mu))
        #ax.semilogx(relative_energy, occupancy_probability, color=color, label=rf"$\mu = {mu}$") # plot energy in log scale
        ax.plot(relative_energy, occupancy_probability, color=color, label=rf"$\mu = {mu}$")
            
    ax.set_xlabel("Relative Binding Energy [kT]")
    ax.set_ylabel("Pr(bound)")
    ax.legend(loc=(1.02, 0), fontsize=11)

    # Show where the 50% cutoff is for mu of 0.1
    ax.axhline(0.5, color="k", linestyle="--")
    ax.axvline(0.1, color="k", linestyle="--")

    return fig, ax

def annotate_heatmap(ax, df, thresh, adjust_lower_triangle=False):
    """
    Display numbers on top of a heatmap to make it easier to view for a reader. If adjust_lower_triangle is True,
    then the lower triangle of the heatmap will display values in parentheses. This should only happen if the heatmap
    is symmetric. Assumes that low values are displayed as a light color and high values are a dark color.

    Parameters
    ----------
    ax : Axes object
        The plot containing the heatmap on which annotations should be made

    df : pd.DataFrame
        The data underlying the heatmap.

    thresh : float
        Cutoff for switching from dark to light colors. Values above the threshold will be displayed as white text,
        those below as black text.

    adjust_lower_triangle : bool
        If True, the lower triangle values will be shown in parentheses.

    Returns
    -------
    None
    """
    for row in range(df.shape[0]):
        for col in range(df.shape[1]):
            value = df.iloc[row, col]
            if value > thresh:
                color = "white"
            else:
                color = "black"

            # Format the value as text
            value = f"{value:.2f}"
            # Add parentheses if desired and in the lower triangle and the heatmap is square
            if adjust_lower_triangle and row < col and df.shape[0] == df.shape[1]:
                value = "(" + value + ")"

            ax.text(row, col, value, ha="center", va="center", color=color)


def motif_rank_order_dots(energy_df, annotate=None, y_column="pred.ddG", ticklabel_size=3, markersize=1, title=None, colors=["#1A1A1A", "#DC9514"], figax=None):
    """
    Helper function that given a list of sequence and their RBE, make a ranked order dot plot.

    Parameters
    ----------
    energy_df: dataFrame
        Dataframe storing sequences and their RBE
    
    annotate: list-like
        List of sequence variants to be highlighted

    y_column: str
        Values in this column will be plotted

    ticklabel_size: float
        Fontsize of sequence labels
    
    markersize: float
        Dot size

    colors: list-like
        List of two colors used for all and highlighted sequences

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

    energy_df = energy_df.copy()
    # initilize an array for y axis
    y_pos = np.arange(len(energy_df.index))

    # add selected sequences on ticks
    ax.set_yticks(y_pos)
    labels = [item.get_text() for item in ax.get_yticklabels()]
    for i,seq in enumerate(energy_df.sequence):
        if annotate: #highlight those interesting ones
            if seq in annotate:
                ax.scatter(y=i,x=energy_df.at[i,y_column], color=colors[1], marker="o", s=markersize)
            else:
                ax.scatter(y=i,x=energy_df.at[i,y_column], color=colors[0], marker="o", s=markersize) 
        else:
            ax.scatter(y=i,x=energy_df.at[i,y_column], color=colors[0], marker="o", s=markersize)
        labels[i] = seq

    ax.yaxis.set_ticks_position("right")
    ax.set_yticklabels(labels, fontsize=ticklabel_size, va="center", ha="left")
    ax.tick_params(axis='y', which='major', pad=2)

    # add color to the labels
    for i,seq in enumerate(energy_df.sequence):
        if annotate: #highlight those interesting ones
            if seq in annotate:
                ax.get_yticklabels()[i].set_color(colors[1])

    # add axis decorations
    ax.set_xlabel("$-\Delta \Delta G$ (kcal/mol)")

    if title:
        ax.set_title(title)

    return fig, ax
