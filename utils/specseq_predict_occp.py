
import os
import re
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn import linear_model

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.font_manager


### I. Scan a  single sequence and predict occupancy and relative affinity given a PWM or EWM

def ewm_to_dict(ewm):
    """Convert a DataFrame representation of an EWM to a dictionary for faster indexing.

    Parameters
    ----------
    ewm : pd.DataFrame

    Returns
    -------
    ewm_dict : {int: {str: float}}
        Dictionary of dictionaries, where the outer keys are positions, the inner keys are letters, and the values
        are values of the matrix
    """
    ewm_dict = ewm.to_dict(orient="index")
    return ewm_dict

# find reverse complement
def rev_comp(seq):
    """Take the reverse compliment of a sequence

    Parameters
    ----------
    seq : str
        The original sequence.

    Returns
    -------
    new_seq : str
        The reverse compliment.
    """
    compliment = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    rev_seq = seq[::-1]
    rev_comp = "".join([compliment[i] for i in rev_seq])
    return rev_comp

def calculate_energy_landscape(f_sequence, ewm, sequence_design=None, consensus=None):
    if type(ewm) is dict:
        ewm_pos = np.array(ewm.keys())
    else:
        ewm_pos = np.array(ewm.index)

    if sequence_design:
        # first check if the seqeunce design and ewm match in both length and positions
        motif_N_pos = np.array([pos+1 for pos,nuc in enumerate(sequence_design) if nuc == "N"])
        motif_nonN_pos = np.array([pos+1 for pos,nuc in enumerate(sequence_design) if nuc != "N"])
        motif_len = len(sequence_design)
    else:
        # if design not specified, infer from ewm, assuming all positions are given
        motif_N_pos = ewm_pos
        motif_nonN_pos = None
        motif_len = len(ewm_pos)

    if np.array_equal(motif_N_pos, ewm_pos):
        # total number of positions to scan per strand
        n_pos = len(f_sequence) - motif_len + 1

        # initialize arrays to store energy landscape
        f_energy = np.zeros(n_pos)
        r_energy = np.zeros(n_pos)

        # retrieve the reverse complement sequences
        r_sequence = rev_comp(f_sequence)

        # scan through DNA, calculate f. and r. energy
        for pos in range(n_pos):
            # get kmer starting at position pos
            f_motif = f_sequence[pos : pos+motif_len]
            r_motif = r_sequence[pos : pos+motif_len]

            """
            # about three times as slow
            f_score = sum([ewm.loc[i+1,nuc] for i,nuc in enumerate(f_motif) if i+1 in motif_N_pos])
            r_score = sum([ewm.loc[i+1,nuc] for i,nuc in enumerate(r_motif) if i+1 in motif_N_pos])
            
            """
            """
            # Initialize energy score
            f_score = 0
            r_score = 0
            
            # energy score f. and r. kmer, accounting of the position offset (+1) and partial ewm (as from specseq)
            for i in motif_N_pos: # for if ewm is dictionary
                fscore += ewm[i][f_kmer[i]]
                rscore += ewm[i][r_kmer[i]]
            """

            # skip sequences containing ambiguous nucleotide N
            if re.search("N", f_sequence):
                f_score=np.nan
                r_score=np.nan
            
            elif motif_nonN_pos: # for partial ewm, need perfect match of constant region
                # first check constant base positions, calcualte energy only if constant regions matches
                if sum(f_motif[i-1]!=sequence_design[i-1] for i in motif_nonN_pos) != 0:
                    f_score = np.nan
                else:
                    f_score = sum([ewm.at[i,f_motif[i-1]] for i in motif_N_pos])
                # do the same thing of motif on the opposite strand
                if sum(r_motif[i-1]!=sequence_design[i-1] for i in motif_nonN_pos) != 0:
                    r_score = np.nan
                else:
                    r_score = sum([ewm.at[i,r_motif[i-1]] for i in motif_N_pos])
            
            else: # for ewm where all positions are specified
                f_score = sum([ewm.at[i,f_motif[i-1]] for i in motif_N_pos])
                r_score = sum([ewm.at[i,r_motif[i-1]] for i in motif_N_pos])
            
            # update energy andscape array
            f_energy[pos] = f_score
            r_energy[pos] = r_score

        # reorder the reverse strand score to match forward strand position
        r_energy = r_energy[::-1]

        # offset by the energy of consensus sequence if given
        if consensus:
            if len(consensus) == len(sequence_design):
                consensus_score = sum([ewm.at[i,consensus[i-1]] for i in motif_N_pos])
                f_energy = f_energy - consensus_score
                r_energy = r_energy - consensus_score
        
    else:
        warnings.warn("ewm does not match sequence design, please check!")
    
    return f_energy, r_energy       

def calcualte_occupancy_landscape(seq, ewm, sequence_design, consensus, mu):
    # essentially converting relative energy landscape to occupancy (probability of bound)
    f_energy, r_energy = calculate_energy_landscape(seq, ewm, sequence_design, consensus)
    # Convert energy scores to occupancies
    f_occupancy = 1 / (1 + np.exp(f_energy - mu))
    r_occupancy = 1 / (1 + np.exp(r_energy - mu))
    return f_occupancy, r_occupancy

def calculate_relaffinity_landscape(seq, ewm, sequence_design, consensus, temp=0):
    # essentially convertting relative energy landscape to relative affinity
    f_energy, r_energy = calculate_energy_landscape(seq, ewm, sequence_design, consensus)
    # calcualte RT, in the unit of kJ/mol
    rt = (temp + 273.15)*8.3145/1000
    # convert energy scores to relative affinity
    f_relaffinity = np.exp(-f_energy/rt)
    r_relaffinity = np.exp(-r_energy/rt)

    return f_relaffinity, r_relaffinity

### II. compiled energy landscape of a list of sequences stored in dataframe, identifier stored as index

def total_landscape(seq, ewms, designs, refs, mu):
    """Compute the occupancy landscape for each TF and join it all together into a DataFrame. Pad the ends of the
    positional information so every TF occupancy landscape is the same length.

    Parameters
    ----------
    seq : str
        The DNA sequence.
    ewms : pd.Series or dict {str: {int: {str: float}}}
        Keys/index are TF names and values are dictionary representations of the EWMs.
    designs: pd.Series, dtype=pd.DataFrame
        Each value of the series is a sequence design string.
    refs: pd.Series, dypte=pd.DataFrame
        Each value of the series is a reference sequence according to ewm.
    mu : int or float
        TF chemical potential.

    Returns
    -------
    landscape : pd.DataFrame, dtype=float
        The occupancy of each TF at each position in each orientation. Rows are positions, columns are TFs and
        orientations, values indicate the predicted occupancy starting at the position.
    """
    # initialize list of None values if sequence design or consensus not specified
    if len(designs) == 0:
        designs = [None] * len(ewms)
    if len(refs) == 0:
        refs = [None] * len(ewms)
    
    landscape = {}
    seq_len = len(seq)
    if type(ewms) is dict:
        keys = ewms.keys()
    else:
        keys = ewms.index
    # For each TF
    for i,name in enumerate(keys):
        # Get the predicted occupancy and add it to the list
        fscores, rscores = calcualte_occupancy_landscape(seq, ewm=ewms[name], sequence_design=designs[i], consensus=refs[i], mu=mu)
        landscape[f"{name}_F"] = fscores
        landscape[f"{name}_R"] = rscores

    # Pad the ends of the lists to the length of the sequence
    for key, val in landscape.items():
        amount_to_add = seq_len - len(val)
        landscape[key] = np.pad(val, (0, amount_to_add), mode="constant", constant_values=0)

    landscape = pd.DataFrame(landscape)

    return landscape


def total_occupancy(seq, ewms, designs, refs, mu):
    """For each TF, calculate its predicted occupancy over the sequence given the energy matrix and chemical
    potential. Then, summarize the information as the total occupancy of each TF over the entire sequence.

    Parameters
    ----------
    seq : str
        The DNA sequence.
    ewms : pd.Series or dict {str: {int: {str: float}}}
        Keys/index are TF names and values are dictionary representations of the EWMs.
    designs: pd.Series, dtype=pd.DataFrame
        Each value of the series is a sequence design string.
    refs: pd.Series, dypte=pd.DataFrame
        Each value of the series is a reference sequence according to ewm.
    mu : int or float
        TF chemical potential.

    Returns
    -------
    occ_profile : pd.Series, dtype=float
        The total occupancy profile of each TF on the sequence.
    """
    occ_landscape = total_landscape(seq, ewms, designs, refs, mu)
    occ_profile = {}
    # Add together F and R strand
    if type(ewms) is dict:
        keys = ewms.keys()
    else:
        keys = ewms.index
    for tf in keys:
        occ_profile[tf] = occ_landscape[[f"{tf}_F", f"{tf}_R"]].sum().sum()

    occ_profile = pd.Series(occ_profile)
    
    return occ_profile


def all_seq_total_occupancy(seq_ser, ewm_ser, design_ser, cons_ser, mu, convert_ewm=False):
    """Calculate the total predicted occupancy of each TF over each sequence.

    Parameters
    ----------
    seq_ser : pd.Series, dtype=str
        Representation of FASTA file, where each value is a different sequence. Index is the FASTA header.
    ewm_ser : pd.Series, dtype=pd.DataFrame
        Each value of the series is an energy matrix for a different TF.
    design_ser: pd.Series, dtype=pd.DataFrame
        Each value of the series is a sequence design string.
    cons_ser: pd.Series, dypte=pd.DataFrame
        Each value of the series is a reference sequence according to ewm.
    mu : int or float
        TF chemical potential.

    Returns
    -------
    occ_df : pd.DataFrame, shape=[n_seq, n_tf]
        Total predicted occupancy of each TF over each sequence. Rows are sequences with same index as seq_ser,
        columns represent different TFs.
    """
    # Convert pd.Series to dictionary representations for speedups
    if convert_ewm:
        ewm_ser = {name: ewm_to_dict(ewm) for name, ewm in ewm_ser.iteritems()}
    #if len(design_ser)>1 and type(design_ser) != dict:
    #    design_ser = {name: seq for name, seq in zip(ewm_ser.keys(), design_ser.values)}
    #if len(cons_ser)>1 and type(cons_ser) != dict:
    #    cons_ser = {name: seq for name, seq in zip(ewm_ser.keys(), cons_ser.values)}

    seq_ser = seq_ser.str.upper()
    print("using mu equals " + str(mu) + " for calculation")
    occ_df = seq_ser.apply(lambda x: total_occupancy(x, ewm_ser, design_ser, cons_ser, mu))

    return occ_df

def save_df(df, outfile):
    """Save a DataFrame to file."""
    # check if directory already exists, create new if not
    if not os.path.exists(os.path.split(outfile)[0]):
            os.mkdir(os.path.split(outfile)[0])
    df.to_csv(outfile, sep="\t", na_rep="NaN")

