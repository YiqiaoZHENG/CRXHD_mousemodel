
import os
import re
import itertools
import math
import numpy as np
import pandas as pd

# I. I/O of fasta files
def fasta_iter(fin, sep=""):
    """A generator function to parse through one entry in a FASTA or FASTA-like file.

    Parameters
    ----------
    fin : file input stream
        Handle to the file to parse
    sep : str
        Delimiter for adjacent bases in the file

    Yields
    -------
    header : str
        Name of the sequence
    sequence : str
        The sequence
    """
    # Generator yields True if on a header line
    generator = itertools.groupby(fin, lambda x: len(x) > 0 and x[0] == ">")
    for _, header in generator:
        # Syntactic sugar to get the header string
        header = list(header)[0].strip()[1:]
        # Get all the lines for this sequence and concatenate together
        sequence = sep.join(i.strip() for i in generator.__next__()[1])
        yield header, sequence


def read_fasta(filename):
    """Parse through a FASTA file and store the sequences as a Series.

    Parameters
    ----------
    filename : str
        Name of the file.

    Returns
    -------
    seq_series : pd.Series, dtype=str
        Index is the FASTA header, values are the sequence strings.
    """
    seq_series = {}
    with open(filename) as fin:
        # Add each sequence to the series
        for header, sequence in fasta_iter(fin):
            sequence = sequence.upper()
            seq_series[header] = sequence
    
    seq_series = pd.Series(seq_series)
    seq_series.index.name = "label"
    return seq_series


def write_fasta(fasta_ser, filename):
    """Write the given series to a file in FASTA format.

    Parameters
    ----------
    fasta_ser : pd.Series
        Index is the FASAT header, values are the sequence strings.
    filename : str
        Name of the file to write to.

    Returns
    -------
    None
    """
    if not os.path.exists(os.path.split(filename)[0]):
        os.mkdir(os.path.split(filename)[0])

    with open(filename, "w") as fout:
        for header, seq in fasta_ser.iteritems():
            fout.write(f">{header}\n{seq}\n")

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
    compliment = {"A": "T", "C": "G", "G": "C", "T": "A"}
    new_seq = seq[::-1]
    new_seq = "".join([compliment[i] for i in new_seq])
    return new_seq





# II. Read and parse raw FIMO files
def peek(fin):
    """ Peek at the next line in a file.

    Parameters
    ----------
    fin : file input stream

    Returns
    -------
    line : str
    """
    pos = fin.tell()
    line = fin.readline()
    fin.seek(pos)
    return line


def gobble(fin):
    """Gobble up lines in the file until we have reached the start of a motif or EOF.

    Parameters
    ----------
    fin : file input stream

    Returns
    -------
    lines : str
        The lines that got gobbled, including newline characters.
    """
    lines = ""

    while True:
        line = peek(fin)
        if len(line) == 0:
            break
        else:
            lines += fin.readline()

    return lines

def build_fimo_table(filename):
    # initialize an empty table with all unique sequence identifiers as index and motif names as column names
    fimo_df = pd.read_csv(filename, sep="\t", header=0)

    # drop motif_alt_id column
    fimo_df = fimo_df.drop(columns='motif_alt_id', inplace=False)

    # drop the last three fimo stats rows
    fimo_df = fimo_df.dropna(axis='index', how='any', inplace=False)
    fimo_df = fimo_df[fimo_df.motif_id != np.NaN]

    # initialize a dataframe to store transformed fimo results
    # formatted as seq x motif, each cell will store the rest fields of the fimo result table if exists, any other cells will contain NaN or 0 (?)
    index_name = pd.Series(fimo_df.sequence_name.unique()).sort_values(ascending=True)
    column_name =  pd.Series(fimo_df.motif_id.unique())

    transformed_fimo = pd.DataFrame(np.zeros((len(index_name),len(column_name))), index=index_name, columns=column_name, dtype=object)

    transformed_fimo.index.name = "label"

    return transformed_fimo

def parse_raw_fimo(filename, sep="\t"):
    # generate empty dataframe
    fimo_table = build_fimo_table(filename)

    with open(filename) as fin:
        # skip header
        fin.readline()

        # Add each sequence to the series
        while True:
            newline = fin.readline()
            # skip annotation and empty lines
            if len(newline) > 0 and newline[0] != "#" and newline != "\n":
                splited_line = newline.strip().split(sep)
                identifer = splited_line[2]
                motif = splited_line[0]
                fimo_output = ",".join(splited_line[3:])

                if fimo_table.at[identifer, motif] == 0.0:
                    fimo_table.at[identifer, motif] = fimo_output
                else:
                    fimo_table.at[identifer, motif] = fimo_table.at[identifer, motif] + ";" + fimo_output
            
            # Check if EOF
            if len(peek(fin)) == 0:
                break
    
    return fimo_table

def split_string_toSeries(string):
    # helper function to split string into Series if the string is not zero
    idx_name = ["start", "end", "strand", "score", "p-value", "q-value", "matched_sequence"]
    score_ser = pd.Series(string).str.split(",", expand=True).set_axis(idx_name, axis='columns').squeeze('index')
    
    return score_ser

def read_fimo_file(filename):
    fimo_df = pd.read_csv(filename, sep="\t").set_index("label")
    # convert all '0.0' to float
    fimo_df = fimo_df.applymap(lambda x: float(x) if x == '0.0' else x)
    
    return fimo_df

# III. Mask matched motif to strings of N
def mask_motif_by_fimo(fasta_ser, fimo_df, motif_name, coremotif_dict=None):
    masked_fasta = {}
    match_count_ser = {}

    for identifier, seq in fasta_ser.iteritems():
        match_count = 0
        if identifier in fimo_df.index:
            # find the fimo result entry if exists
            fimo_matches = fimo_df.at[identifier, motif_name]
            # only extract information from sequence with motif match
            if fimo_matches != 0.0:
                # single match
                if fimo_matches.count(";") == 0:
                    fimo_matches = [fimo_matches]
                # multiple matches
                else:
                    fimo_matches = fimo_matches.split(";")
                # iterate through all matches
                for match in fimo_matches:
                    apply_mask = False
                    # parse the fimo results
                    fimo_ser = split_string_toSeries(match)
                    
                    match_seq = fimo_ser.matched_sequence
                    # if a core motif dictionary is specified, mask only sequences with perfect match of the core motif
                    # otherwise, mask the entire fimo found matched sequence
                    if coremotif_dict:
                        matched_core = "".join([match_seq[pos] for pos in coremotif_dict.keys()])
                        ref_core = "".join([nuc for nuc in coremotif_dict.values()])
                        if matched_core == ref_core:
                            apply_mask = True
                    else:
                        apply_mask = True
                    
                    start = math.trunc(float(fimo_ser.start))
                    end = math.trunc(float(fimo_ser.end))
                    #strand = str(fimo_ser.strand)
                    #score = float(fimo_ser.score)

                    # mask the matched substring into strings of N in the fasta record
                    if apply_mask:
                        masked_seq = seq.replace(seq[start-1:end], "N"*len(match_seq))
                        # update the fasta sequence
                        seq = masked_seq
                        match_count += 1

        # update the fasta and match count dictionary
        masked_fasta[identifier] = seq
        match_count_ser[identifier] = match_count

    # convert dictionary to Series
    masked_fasta = pd.Series(masked_fasta)
    masked_fasta.index.name = "label"

    match_count_ser = pd.Series(match_count_ser)
    match_count_ser.index.name = "label"

    return match_count_ser, masked_fasta



#### something sad...  i tried generator function to parse fimo table, but it is slow as hell hhhh
## %timeit sequence_annotator.parse_raw_fimo(f) 386 ms ± 11.7 ms per loop (mean ± std. dev. of 7 runs, 1 loop each
## %timeit read_raw_fimo(f) 2.69 s ± 22.3 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
# now use a generator function to iterate through the entire fimo result table
def read_raw_fimo(filename, sep="\t"):
    # initialize an empty table with all unique sequence identifiers as index and motif names as column names
    fimo_df = pd.read_csv(filename, sep="\t", header=0)
    # drop motif_alt_id column
    fimo_df = fimo_df.drop(columns='motif_alt_id', inplace=False)

    # initialize an empty multi index dataframe,
    mesh = np.array(np.meshgrid(np.array(fimo_df.sequence_name.unique()), np.array(fimo_df.motif_id.unique())))
    combinations = mesh.T.reshape(-1, 2)
    arrays = list(zip(*combinations))
    tuples = list(zip(*arrays))

    index = pd.MultiIndex.from_tuples(tuples, names=["label", "motif_id"])
    fimo_df = pd.DataFrame(np.zeros(len(tuples)), index=index, columns=["fimo_score"])

    with open(filename) as fin:
        # skip header
        fin.readline()

        # generator for the rest of the file
        fimo_gen = (row for row in fin)
        # get next line
        newline = next(fimo_gen)

        while True:
            # retrieve data
            splited_line = newline.strip().split(sep)
            index = (splited_line[2],  splited_line[0])
            score = ",".join(splited_line[3:])

            # update the fimo score dataframe
            #for index, score in fimo_iter(fin, sep="\t"):
            if fimo_df.loc[index,"fimo_score"] == 0.0:
                fimo_df.loc[index,"fimo_score"] = score
            else:
                fimo_df.loc[index,"fimo_score"] += ";" + score

            newline = next(fimo_gen)
            if len(newline) == 0 or newline[0] == "#" or newline == "\n":
                break
                #newline = next(fimo_gen) # skip unmeaningful lines

    return fimo_df


# IV. Functions for motif associated analysis

# fimo

def get_sum_fimo_score(fimo_scores):
    #sum = np.array([float(sequence_annotator.split_string_toSeries(s)["score"]) for s in fimo_scores]).sum()
    sum = np.array([float(s.split(",")[3]) for s in fimo_scores]).sum()
    return sum


def get_max_fimo_score(fimo_scores):
    #max = np.array([float(sequence_annotator.split_string_toSeries(s)["score"]) for s in fimo_scores]).max()
    max = np.array([float(s.split(",")[3]) for s in fimo_scores]).max()
    return max


def get_mean_fimo_score(fimo_scores):
    #avg = np.array([float(sequence_annotator.split_string_toSeries(s)["score"]) for s in fimo_scores]).mean()
    avg = np.array([float(s.split(",")[3]) for s in fimo_scores]).mean()
    return avg


def extract_fimo_score(fimo_scores, get_sum=True, get_max=False, get_mean=False):
     # multiple matches
    if fimo_scores.count(";") > 0:
        fimo_scores = fimo_scores.split(";")
        if not get_max and not get_mean: # default
            score = get_sum_fimo_score(fimo_scores)
        elif get_max:
            score = get_max_fimo_score(fimo_scores)
        else:
            score = get_mean_fimo_score(fimo_scores)
    # single match
    else:
        score = float(fimo_scores.split(",")[3])

    return score


def count_fimo_score(fimo_scores):
    count = fimo_scores.count(";") + 1
    return count


def annot_motif_count(data, fimo_table):
    fimo_table = fimo_table.copy()
    count_df = fimo_table.applymap(lambda x: x if x == 0.0 or x[0:4] == "peak" else count_fimo_score(fimo_scores=x))
    data = pd.merge(data, count_df, how="left", on="peak.id").reset_index(drop=True)
    
    return data


def annot_motif_score(data, fimo_table, **kwargs):
    fimo_table = fimo_table.copy()
    score_df = fimo_table.applymap(lambda x: x if x == 0.0 or x[0:4] == "peak" else extract_fimo_score(fimo_scores=x, **kwargs))
    data = pd.merge(data, score_df, how="left", on="peak.id").reset_index(drop=True)

    return data


def annot_pred_occp(data, occp_table, **kwargs):
    return data


def annotate_rnalfc_by_motif(annot_df, scoring_df, annot_list=["motif_score"], score_type=dict(get_mean=True)):
    # the input dataframe should at least contain peak.id column  
    data_to_plot = annot_df.copy()

    if "motif_count" in annot_list:
        data_to_plot = annot_motif_count(data=data_to_plot, fimo_table=scoring_df)
    
    if "motif_score" in annot_list:
        data_to_plot = annot_motif_score(data=data_to_plot, fimo_table=scoring_df, **score_type)
    
    if "occupancy" in annot_list:
        data_to_plot = annot_pred_occp(data=data_to_plot)

    # sort by peak.id
    data_to_plot = data_to_plot.sort_values(by="peak.id").reset_index(drop=True)
    
    return data_to_plot

