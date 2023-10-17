"""
Author: Yiqiao Zheng
Email: yiqiao.zheng@wustl.edu
"""

import os
import re
import itertools
import subprocess

import numpy as np
import pandas as pd
from scipy import stats
from sklearn import linear_model


# Define library design (update as new libraries added)
lib_designs = {
    # for automatically adding filter sequences
    "M":"TAANNN",
    "Mrev":"NNNTTA",
}


### I. Build energy model from Specseq RBE data ####
# essentially a python sklearn implementation of Zuo's TFCookbook R package main functions: https://github.com/zeropin/TFCookbook

def parse_by_library(seq_to_energy_df, lib="M", band ="m", update_bool=False):
    """
    Parse full RBE dataframe by library, find and use matched energy column name as identifier
    
    Parameters
    ----------
    seq_to_energy_df : dataframe
        Dataframe containing average RBE data. At least four columns are required: sequence, lib, MMcount, avg.{lib}.{band}.ddG
    
    lib: str
        Identifier of specseq library

    band: str
        Identifier of EMSA band
    
    update_bool: bool
        Boolean indicating whether the requested lib x band RBE data existed. Default: False

    Returns
    -------
    parsed_energy_matrix: dataframe
        A parsed dataframe specified by library and band 

    update_bool: bool
        Boolean indicating whether the requested lib x band RBE data existed
    """
    # check if requested band data exists, skip if not, format: 'avg.M.mddG'
    name_to_search = "avg." + lib + "." + band +"ddG"
    if name_to_search in seq_to_energy_df.columns:
        # only extract columns with matched lib x band energy data
        parsed_energy_matrix= seq_to_energy_df.loc[seq_to_energy_df["lib"]==lib, ["sequence","lib","MMcount",name_to_search]]
        # update the ewm matrix when at least one library is found and properly parsed
        update_bool = True
    else:
        print("\t- "+name_to_search + " not found in the given dataframe, skipping!")
        
    return parsed_energy_matrix, update_bool


def seq_to_byte(sequence):
    """
    Convert sequence character to bytes for regression
    
    Parameters
    ----------
    sequence: str
        A string of DNA characters [ACGTN-]
    
    Returns
    -------
    full_seq_inbyte: list
        A list of byte encoded DNA characters

    """
    # split sequence into list of individual characters
    nucleotides = list(sequence)
    # define nucleotide to byte mapping
    base_to_byte={  "A": [0,1,0],
                    "C": [0,0,0],
                    "G": [1,0,0],
                    "T": [0,0,1],
                    "N": [0.25,0.25,0.25],
                    "-": [0,0,0]}
    # format and generate full length byte map for the given sequence, length should be len(nucleotide)*3
    all_nucleotides_inbyte=[base_to_byte[k] for k in nucleotides]
    all_nucleotides_inbyte=[str(n)for m in all_nucleotides_inbyte for n in m]
    
    full_seq_inbyte = ""
    full_seq_inbyte = full_seq_inbyte.join(all_nucleotides_inbyte)

    return full_seq_inbyte


def buildEnergyModel(seq_to_energy_df, lib="M", band ="m", MM_th=2):
    """
    Build energy model using multiple linear regression
    
    Parameters
    ----------
    seq_to_energy_df : dataframe
        Dataframe containing average RBE data. At least four columns are required: sequence, lib, MMcount, avg.{lib}.{band}.ddG
    
    lib: str
        Identifier of specseq library

    band: str
        Identifier of EMSA band
    
    MM_th: int
        Mismatch cutoff of sequence variants to be included for generating EWM. Default: 2

    Returns
    -------
    scheme: string
        Energy model scheme

    regr: object
        Fitted estimator of the Ordinary least squares Linear Regression
        Will return 0 if the specified RBE data can not be found
    """
    num = len(seq_to_energy_df.loc[seq_to_energy_df.lib==lib,'sequence'].reset_index(drop=True)[0])

    # generate energy model scheme, len(sequence) * 3
    bases = ["CG", "CA", "CT"]*num
    positions = [str(k) for k in np.repeat(list(range(1,num+1,1)),3)]
    scheme = np.char.add(positions,bases)

    # parse the dataframe library and band
    parsed_df, update_bool = parse_by_library(seq_to_energy_df, lib=lib, band=band)

    # only proceed if the specifed lib x band energy data is found
    if update_bool:
        # select only sequences within 2 mismatches to consensus
        parsed_df = parsed_df.loc[parsed_df["MMcount"]<=MM_th] # drop all the rows not matching condition
        parsed_df.reset_index(drop=True, inplace=True)

        # retrieve the energy column name, should be the last column
        energy_colname = parsed_df.columns[-1]

        print("  - fitting energy data " + energy_colname + " with linear model")

        # convert sequence characters to byte encoding
        byte_df = parsed_df.sequence.apply(lambda s: pd.Series({'ByteCode':seq_to_byte(s)}))
        # concatenate with the sequence-energy dataframe
        byte_to_energy_df = pd.concat([parsed_df,byte_df], axis=1)
        # select only the energy and byte code columns
        byte_to_energy_df = byte_to_energy_df[[energy_colname ,"ByteCode"]]

        # divide the byte code column by specified scheme, each element in the scheme maps to a single code
        splited_byte = byte_to_energy_df.ByteCode.apply(lambda x: pd.Series(list(x)))
        # names the columns as specified in scheme
        splited_byte.columns = scheme
        # concatenate the two dataframes
        byte_to_energy_df = pd.concat([byte_to_energy_df.loc[:, byte_to_energy_df.columns != 'ByteCode'],splited_byte], axis=1)

        # select only necessary columns (should be written outside of this function)
        byte_to_energy_df = byte_to_energy_df.rename(columns={energy_colname : "energy"})
        byte_to_energy_df = byte_to_energy_df[["energy"]+list(scheme)]

        # Ordinary least squares Linear Regression
        regr = linear_model.LinearRegression()
        # fit linear model
        regr.fit(byte_to_energy_df[list(scheme)], byte_to_energy_df[["energy"]])

        return scheme,regr

    else:
        return scheme, 0
        
        
def retreieveEnergyModel(scheme, linear_model):
    """
    Retrieve coefficients from linear regression model and format it as a standardized energyMatrix
    
    Parameters
    ----------
    scheme: string
        Energy model scheme
    
    regr: object
        Fitted estimator of the Ordinary least squares Linear Regression

    Returns
    -------
    energyMatrix: dataframe
        Dataframe containing the EWM generated by the linear model
        Position numbers in the EWM is offset by 1 to match the exact position
    """
    # generate nucleotide to coefficients map
    coeff_map = pd.DataFrame(list(zip(scheme, linear_model.coef_[0])), columns=['pos', 'coeffs']).set_index("pos")

    # initialize an empty matrix of 4*len(sequence) to store the energy matrix
    matrix_length = int(len(coeff_map)/3)+1
    energyMatrix = pd.DataFrame(index=["A", "C", "G", "T"], columns=range(1,matrix_length,1))
    # fill the energy matrix with coefficients from the linear model
    for i in range(1,matrix_length,1):
        energyMatrix.loc["C", i] = 0
        energyMatrix.loc["G", i] = coeff_map.loc[str(i)+"CG","coeffs"]
        energyMatrix.loc["A", i] = coeff_map.loc[str(i)+"CA","coeffs"]
        energyMatrix.loc["T", i] = coeff_map.loc[str(i)+"CT","coeffs"]
    
    # drop columns with all zeros, since it is never absolute 0, use a very small number for comparison
    energyMatrix = energyMatrix.loc[:, (abs(energyMatrix)>=0.0001).any()].loc[:, (abs(energyMatrix)<=10).any()].T
    # name the index column
    energyMatrix.index.name = "pos"

    # convert all numbers to float dtype
    energyMatrix = energyMatrix.astype(float)

    return energyMatrix


def ewm_from_RBE(RBE_matrix, lib_list=["M"], find_band="m", MMth=2, normalize=False):
    """
    Wrapper function to calculate ewm from raw RBE data
    
    Parameters
    ----------
    RBE_matrix : list
        List of dataframe containing average RBE data. 
        At least four columns are required in the dataframe: sequence, lib, MMcount, avg.{lib}.{band}.ddG
    
    lib: str
        Identifier of specseq library

    band: str
        Identifier of EMSA band
    
    MM_th: int
        Mismatch cutoff of sequence variants to be included for generating EWM. Default: 2

    normalize: bool
        Boolean to indicate whether to normalize the EWM. Default: False

    Returns
    -------
    energy_models: list
        List of EWMs

    """
    energy_models = {}
    for i in RBE_matrix.index:
        for find_lib in lib_list:
            id = ".".join([find_lib,find_band])
            print(".".join([i,id]))
            if any(id in string for string in RBE_matrix[i].columns.tolist()):
                scheme, model = buildEnergyModel(RBE_matrix[i], find_lib, find_band, MMth)
                # update the ewm dictionary and only if a linear model has been properly generated
                if model != 0:
                    ewm = retreieveEnergyModel(scheme, model)
                    ewm = normalize_ewm(ewm)
                    if not normalize:
                        ewm = denormalize_ewm(ewm)    
                    energy_models[".".join([i,id])] = ewm

    # convert library to series if not empty
    if energy_models:
        energy_models = pd.Series(energy_models, dtype=object)

    return energy_models

def denormalize_ewm(energyMatrix):
    """Reverse normalize EWM so that the favourite base at each position has energy 0"""
    # check if row sum is not 0
    if energyMatrix.sum(axis=1).sum() <= 0.0001:
        # normalize to the lowest energy/highest affinity base
        energyMatrix = energyMatrix.apply(lambda x: x-min(x), axis=1)
    return energyMatrix

def normalize_ewm(energyMatrix):
    """Normalize EWM to make energy sum at each position 0"""
    # check if row sum is not 0
    if not energyMatrix.sum(axis=1).sum() <= 0.0001:
        # normalization to make energy sum 0 
        energyMatrix = energyMatrix.apply(lambda x: x-0.25*sum(x), axis=1)
    return energyMatrix


### II. Build energy model with specseq_mlr.R script, using lm() function in R ####

def get_motif_revcomp(pwm):
    """Get reverse complement of a sequence"""
    for_pwm = pwm.copy()
    for_pwm = for_pwm[::-1].reset_index(drop=True)
    pwm_rc = for_pwm.copy()
    pwm_rc["A"] = for_pwm["T"]
    pwm_rc["C"] = for_pwm["G"]
    pwm_rc["G"] = for_pwm["C"]
    pwm_rc["T"] = for_pwm["A"]
    pwm = pwm_rc

    # offset the index by 1 to match base position
    pwm.index = range(1, len(pwm.index)+1,1)
    pwm.index.name = "pos"

    return pwm


def read_pwm_motif(path_to_file, header=False, get_revcomp=False):
    """
    Read in one MEME pwm from file"""
    # infer identifier name from file name
    pwm_name = os.path.basename(path_to_file).strip(".txt")
    # if header exists, udpate the motif identifier
    if header:
        with open(path_to_file, "r") as f:  
            pwm_name = f.readline()[1:].strip()
            f.close()
        print("retrievig pwm for: "+pwm_name)
        pwm = pd.read_csv(path_to_file, sep='\s+', names=["A","C","G","T"], skiprows=[0]) # skip the header line
    else:
        print("retrievig pwm for: "+pwm_name)
        pwm = pd.read_csv(path_to_file, sep='\s+', names=["A","C","G","T"])
    
    # get the reverse complement if specified
    if get_revcomp:
        pwm = get_motif_revcomp(pwm)

    # offset the index by 1 to match base position
    pwm.index = range(1, len(pwm.index)+1,1)

    return pwm_name,pwm


def run_specseq_mlr(scriptdir, file, lib, band, MMth=2, write_file_bool=False):
    """
    Wrapper function to call the specseq_mlr.R script and get ewm back from stdout

    Parameters
    ----------
    scriptdir: string
        Path to the R script
    
    file: string
        Path to output file

    lib: str
        Identifier of specseq libary to find

    band: str
        Identified of EMSA band to find
    
    MM_th: int
        Mismatch cutoff of sequence variants to be included in generating EWM. Default: 2

    file_bool: boolean
        Boolean indicating if write output to file or return to stdout. Default: False

    Returns
    -------
    ewm_name: str
        Identifier of EWM. Return 0 if an error occurs.
    
    ewm: dataframe
        Dataframe containing of calcualted EWM. Return 0 if an error occurs.
    """
    # run multiple linear regression to fit energy data
    mlr_output = read_specseq_mlr_stdout(scriptdir, file, lib, band, MMth, write_file_bool)
    # retreive the stdout only if run sucessfully
    if not mlr_output.returncode:
        if write_file_bool:
            # ewm written to file, output file directory returned in CompletedProcess.stdout
            path_to_mlr_ewm = mlr_output.stdout.rstrip()
            ewm_name, ewm = read_pwm_motif(path_to_mlr_ewm, header=True)
        else:
            # ewm returned in CompletedProcess.stdout
             ewm_name, ewm = get_ewm_from_compProcess(mlr_output)

        return ewm_name, ewm
    else:
        print(" - specseq_mlr.R run failed, please check")
        return 0,0

 
def read_specseq_mlr_stdout(scriptdir, file, lib, band, MMth=2, write_file_bool=False):
    """
    Call specseq_mlr.R and capture the output from stdout

    Parameters
    ----------
    scriptdir: string
        Path to the R script
    
    file: string
        Path to output file

    lib: str
        Identifier of specseq libary to find

    band: str
        Identified of EMSA band to find
    
    MM_th: int
        Mismatch cutoff of sequence variants to be included in generating EWM. Default: 2

    file_bool: boolean
        Boolean indicating if write output to file or return to stdout. Default: False

    Returns
    -------
    mlr_output: object
        CompletedProcess instance from specseq_mlr.R
    """
    # inputfile, library to fit, band to fit, mismatch threshold, all in str dypte
    args = [file, lib, band, str(MMth), str(write_file_bool)]

    # Build subprocess command
    cmd = ['Rscript', os.path.join(scriptdir, "specseq_mlr.R")] + args

    # run the command and capture the output
    mlr_output = subprocess.run(cmd, universal_newlines=True, capture_output=True)

    return mlr_output


def get_ewm_from_compProcess(compPro_object):
    """
    Retrieve ewm identifier and ewm from CompletedProcess output when not write to file
    
    Parameters
    ----------
    compPro_object: object
        CompletedProcess instance from specseq_mlr.R
    
    Returns
    -------
    identifier: str
        Identifier/name of ewm

    ewm: dataframe
        Spec-seq ewm

    """
    df = compPro_object.stdout
    df = df.split("\n")
    identifier = df[0]
    col_names = list(re.split(' +', df[1].lstrip()))
    ewm = {}
    for entry in df[2:-1]:
        splited_entry = list(re.split(' +', entry)[1:])
        ewm[splited_entry[0]] = splited_entry[1:]
    # convert dictionary to dataframe and format
    ewm = pd.DataFrame.from_dict(ewm, orient='index')
    # add column names
    ewm.columns = col_names[1:]
    ewm.index.name = col_names[0]
    
    return identifier, ewm


def mlr_ewm_from_RBE(scriptdir, RBE_filename, lib_list=["M"], find_band="m"):
    """
    Wrapper function to generate ewm models from one or more RBE data files

    Parameters
    ----------
    scriptdir: str
        File path to specseq_mlr.R script

    RBE_filenames: list
        List of file paths to average RBE data

    lib_list: list
        List of Spec-seq libraries to generate ewm

    find_band: str
        Identifier of EMSA band to find
    
    Returns
    -------
    energy_models: pandas.Series
        1d array storing all calculated ewms.
        Indexes are automatically generated and formatted as {filename prefix}.{lib}.{band}
    """
    energy_models = {}
    for f in RBE_filename:
        for find_lib in lib_list:
            i = os.path.split(f)[1].split("_avgRBE.txt")[0]
            id = ".".join([find_lib,find_band])
            print(".".join([i,id]))
            fline=open(f).readline().rstrip().split()
            if any(id in string for string in fline):
                ewm_name, ewm = run_specseq_mlr(scriptdir, f, find_lib, find_band, MMth=2, write_file_bool=False)
                # store the identifier and ewm
                if(ewm_name):
                    print(" - mlr ewm successfully extracted for " + ewm_name)
                    energy_models[ewm_name] = ewm
        
    # convert library to series if not empty
    if energy_models:
        energy_models = pd.Series(energy_models, dtype=object)

    return energy_models

### III. Read EWMs from text file

def peek(fin):
    """
    Peek at the next line in a file.

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


def gobble(fin, wm_type="ewm"):
    """
    Gobble up lines in the file until we have reached the start of a motif or EOF. Note EWM and MEME PWM files are formatted differently and thus are handled separately.

    Parameters
    ----------
    fin : file input stream

    Returns
    -------
    lines : str
        The lines that got gobbled, including newline characters.
    """
    lines = ""
    if wm_type == "ewm":
        while True:
            line = peek(fin)
            if len(line) == 0 or line[0] == ">":
                break
            else:
                lines += fin.readline()
    elif wm_type == "pwm":
        while True:
            line = peek(fin)
            if len(line) == 0 or line[:5] == "MOTIF":
                break
            else:
                lines += fin.readline()

    return lines


def read_ewm_files(filename):
    """
    Given a summary energy model file, read in all EWMs. 
 
    Parameters
    ----------
    filename: str
        Path to EWM file

    Returns
    -------
    ewm_ser: pandas Series
        1d array storing EWMs as DataFrames and indexess are motif identifiers
    
    """
    ewm_ser = {}
    with open(filename) as fin:
        # Lines before the first motif is encountered
        gobble(fin, "ewm")

        # Do-while like behavior to read in the data
        # Do <read in motif> while not EOF
        while True:
            # ><motif name>
            motif_id = fin.readline().strip()[1:]
            print('reading ewm for ' + motif_id)

            # Every line that starts with a integer is a new position in the EWM, if the first character is not a integer
            # it is not part of the EWM.
            ewm = []
            while len(peek(fin)) != 0 and peek(fin)[0] != ">":
                ewm.append(fin.readline().split())

            # Make a DataFrame and add to the list
            ewm = pd.DataFrame(ewm, dtype=float, columns=["pos", "A", "C", "G", "T"])
            # make the position column as index
            ewm = ewm.astype({'pos': 'int'}).set_index('pos')
            ewm_ser[motif_id] = ewm
            # Read up any extra info such as the URL
            gobble(fin, "ewm")

            # Check if EOF
            if len(peek(fin)) == 0:
                break

    ewm_ser = pd.Series(ewm_ser)
    return ewm_ser

def read_meme_files(filename):
    """
    Given a summary meme motif file, read in all PWMs.

    Parameters
    ----------
    filename: str
        Path to PWM file

    Returns
    -------
    pwm_ser: pandas Series
        1d array storing PWMs as DataFrames and indexes are motif identifiers
    """
    pwm_ser = {}
    with open(filename) as fin:
        # Lines before the first motif is encountered
        gobble(fin, "pwm")

        # Do-while like behavior to read in the data
        # Do <read in motif> while not EOF
        while True:
            # MOTIF <motif name>
            motif_id = fin.readline().split()[1]
            print('reading pwm for ' + motif_id)
            
            # Empty line
            fin.readline()
            # "letter-probability matrix: [other info]"
            fin.readline()

            # Every line that starts with a space is a new position in the PWM, if the first character is not a space
            # it is not part of the PWM.
            pwm = []
            while peek(fin)[0] == " ":
                pwm.append(fin.readline().split())
            
            # Make a DataFrame and add to the list
            pwm = pd.DataFrame(pwm, dtype=float, columns=["A", "C", "G", "T"])

            # offset the index by 1 to match base position
            pwm.index = range(1, len(pwm.index)+1,1)

            pwm_ser[motif_id] = pwm
            # Read up any extra info
            gobble(fin, "pwm")

            # Check if EOF
            if len(peek(fin)) == 0:
                break

    pwm_ser = pd.Series(pwm_ser)
    return pwm_ser


def save_ewm(energyMatrices, output_file):
    """
    Write an array of EWMs to text file
    
    Parameters
    ----------
    energyMatrices: pandas Series
        1d array storing EWMs as DataFrames and indexes are motif identifiers
    
    output_file: str
        Path to output text file

    Returns
    -------
    None
    """
    fout = open(output_file,'w')
    for id,mtx in energyMatrices.items():
        # write the ewm identifier
        fout.write(">"+id+"\n")
        # write the entire ewm matrix with ewm positions using pandas
        for line in mtx.index:
            fout.write(str(line) + "\t" + "\t".join([str(i) for i in mtx.loc[line,:]])+"\n")
    fout.close()

### III. Compare predicted with observed relative binding energy

def EnergyModel_to_Coeffs(energyMatrix):
    """
    Convert energy matrix to energy model coefficients

    Parameters
    ----------
    energyMatrix: dataFrame
        Dataframe storing EWMs in ACGT format

    Returns
    -------
    coeff_map: dataFrame
        Dataframe storing the energy model as coefficients in 3L+1 encoding
    """
    # initialize an empty matrix of length 3*len(sequence in ewm) to store the coefficients
    coeff_map = pd.DataFrame([str(a)+b for a,b in zip(np.repeat(list(energyMatrix.index),3),["CA","CG","CT"]*3)], columns=["pos"])
    coeff_map["coeffs"] = 0
    coeff_map = coeff_map.set_index("pos")

    # now fill the coefficients column based on energy matrix
    for i in energyMatrix.index:
      coeff_map.loc[str(i) + "CG", "coeffs"] = energyMatrix.loc[i, "G"] - energyMatrix.loc[i, "C"]
      coeff_map.loc[str(i) + "CA", "coeffs"] = energyMatrix.loc[i, "A"] - energyMatrix.loc[i, "C"]
      coeff_map.loc[str(i) + "CT", "coeffs"] = energyMatrix.loc[i, "T"] - energyMatrix.loc[i, "C"]

    return coeff_map


def design_library(sequence_design):
    """
    Helper function to generate all k-mer sequences based on given seqeuence design

    Parameters
    sequence_design: str
        Spec-seq DNA library design. Only ACGTN are allowed.

    Returns
    -------
    libary: dataFrame
        Single column dataframe containing all k-mers based on library design
    """
    # first check if all the based in the sequence design is in [ACGTN]
    if bool(re.match("^[ACGTN]+$", sequence_design)):
        # count the number of variable nucleotides and generate kmer
        k = sequence_design.count("N")
        kmer_df = pd.DataFrame([list(p) for p in itertools.product(['A','T','G','C'], repeat=k)])
        
        # note the positions of variabe nucleotides
        N_pos = [pos for pos,nuc in enumerate(sequence_design) if nuc == "N"]
        # attach the N positions as column names to the kmer dataframe
        kmer_df.columns = N_pos

        # initialize an empty dataframe to store all sequence 
        library = pd.DataFrame(index=range(4**k),columns=range(len(sequence_design)))
        for pos,nuc in enumerate(sequence_design):
            # first fill columns corresponding to constant bases in the design
            if nuc != "N":
                library[pos] = nuc
            # then fill columns corresponding to variable based
            else:
                library[[pos]] = kmer_df[[pos]]
        
        # now concatenate all based into a single string
        library['sequence'] = library.apply(''.join, axis=1)
        # discard all columns except full length sequence
        library = library[["sequence"]]

        return library
    
    else:
        print("None ACGTN characters found. No library generated. Please check!")

        return 0


def old_predict_bindingEnergy(sequence_design, coeffs):
    """
    Predict relative binding energy for given sequnce design and ewm (deprecated, using 3L+1 encoding, a little bit complicated)
    
    Parameters
    ----------
    sequence_design: str
        Spec-seq DNA library design. Only ACGTN are allowed.

    coeffs: dataFrame
        Dataframe storing the energy model as coefficients in 3L+1 encoding

    Returns
    -------
    predicted_energy_df: dataFrame
        Two column dataframe storing all k-mers specified by library design each with a predicted RBE
    
    """
    # first check if the seqeunce design and coeffs match in both length and positions
    N_pos = np.array([pos+1 for pos,nuc in enumerate(sequence_design) if nuc == "N"])
    coeff_pos = np.unique([int(p[:-2]) for p in coeffs.index])

    if bool(N_pos.all() == coeff_pos.all()):
        # generate all oligos based on sequence design
        sequence_df = design_library(sequence_design)
        # convert sequence to byte encodeing
        sequence_df["ByteCode"] = sequence_df.sequence.apply(lambda s: seq_to_byte(s))
        sequence_df = sequence_df.set_index("sequence")

        # generate energy model scheme, len(sequence) * 3
        bases = ["CG", "CA", "CT"]*len(sequence_design)
        positions = [str(k) for k in np.repeat(list(range(1,len(sequence_design)+1,1)),3)]
        scheme = np.char.add(positions,bases)

        # divide the byte code column by specified scheme, each element in the scheme maps to a single code
        splited_byte = sequence_df.ByteCode.apply(lambda x: pd.Series(list(x)))
        # names the columns as specified in scheme
        splited_byte.columns = scheme
        # keep columns that are specified in the coefficient dataframe
        splited_byte = splited_byte[coeffs.index]

        # calcualte the binding energy by taking matrix product of the byte and coefficients matrix
        sequence_df["pred.ddG"] = 0.0001
        for seq in splited_byte.index:
            sequence_df.at[seq,"pred.ddG"] = sum(int(a)*b for a,b in zip(splited_byte.loc[seq,:], coeffs.coeffs))

        predicted_energy_df = sequence_df[["pred.ddG"]]

        return predicted_energy_df
    else:
        print("ewm does not match sequence design, please check!")

        return 0

# predict relative binding energy for given sequnce design and ewm
def predict_bindingEnergy(sequence_design, ewm):
    """
    Predict relative binding energy for given sequnce design and ewm.

    sequence_design: str
        Spec-seq DNA library design. Only ACGTN are allowed.

    ewm: dataFrame
        EWM model

    Returns
    -------
    sequence_df: dataFrame
        Two column dataframe storing all k-mers specified by library design each with a predicted RBE.
        The base position is offset by 1 to match the exact position.
    """
    # first check if the seqeunce design and ewm match in both length and positions
    N_pos = np.array([pos+1 for pos,nuc in enumerate(sequence_design) if nuc == "N"])
    ewm_pos = np.array(ewm.index)

    if bool(N_pos.all() == ewm_pos.all()):
        # generate all oligos based on sequence design
        sequence_df = design_library(sequence_design)
        # initialize pred.ddG column
        sequence_df["pred.ddG"] = 0.0001
        # set sequence as index for easy access
        sequence_df = sequence_df.set_index("sequence")
        # iterate through all sequences and calculate and update energy
        for seq in sequence_df.index:
            sequence_df.at[seq,"pred.ddG"] = sum([ewm.at[i,seq[i-1]] for i in ewm.index]) # account for the position offset
            
        return sequence_df
    else:
        print("ewm does not match sequence design, please check!")

        return 0


### IV. Convesion bewteen frequency and energy based pwm

def pwm_to_ewm(probabilityMatrix, pseudocount = 0.0001, temp=0, normalize=False):
    """
    Given a position frequency matrix matrix, convert it to energy weight matrix.

    Parameters
    ----------
    probabilityMatrix: dataFrame
        Position frequency matrix matrix
    
    pseudocount: float
        Pseudocount to be added to each number in the matrix
    
    temp: float
        Temperature in Celsius
    
    normalize: boolean
        Boolean to indicate whether to normalize the converted EWM. Default: False

    Returns
    -------
    energyMatrix: dataFrame
        EWM dataframe
    """
    probabilityMatrix = probabilityMatrix.copy()
    probabilityMatrix += pseudocount
    # Normalize each position by the most frequent letter to get relative Kd
    probabilityMatrix = probabilityMatrix.apply(lambda x: x / x.max(), axis=1)
    # calculate RT (in unit of kJ/mol)
    rt = (temp + 273.15)*8.3145/1000
    # Convert to EWM
    energyMatrix = -rt * np.log(probabilityMatrix)
    if normalize:
        # Normalization to make energy sum 0
        energyMatrix = normalize_ewm(energyMatrix)
    return energyMatrix

# by definition, probability of bound ~ occupancy, essentially, both are in relative terms, relative to the consensus
def ewm_to_pwm(energyMatrix, mu):
    """
    Given a energy weight matrix, convert it to position frequency matrix.

    Parameters
    ----------
    energyMatrix: dataFrame
        Energy weight matrix
    
    mu: float
        Chemical potential

    Returns
    -------
    probabilityMatrix: dataFrame
        Position frequency matrix
    
    """
    # Normalize each position by the least energy base
    energyMatrix = energyMatrix.apply(lambda x: x - x.min(), axis=1)
    # Convert EWM scores to occupancies
    probabilityMatrix = 1/(energyMatrix.sub(mu).applymap(np.exp)+1)
    # Normalize to probability sum 1
    probabilityMatrix = probabilityMatrix.apply(lambda x: x/x.sum(), axis=1)
    
    return probabilityMatrix
