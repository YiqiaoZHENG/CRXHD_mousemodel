#!/usr/bin/env python3

"""
read in raw Miseq R1 read files of monomeric and unbound fractions of the same lane 
calcualte the ratio of B/Ub for each band and relative binding energy (-ln(B/Ub)) in the unit of kT

Usage: python3 calculate_ratios.py --sample [sample_name] --lib ["M",...] --ref ["TAATCC",...] --prefix [output parent directory]

Args: 
        --sample = path to fasta file of dimer bound fraction
        --lib = path to fasta file of monomer bound fraction
        --ref = path to fasta file of unboudn fraction
        --prefix = 

Outputs:
        counts/{sample_name}_counts.txt  = summary text file of counts compiled from all bands
        RBE/{sample_name}_RBE.txt = summary text file of ratios and relative binding energy compiled from all bands

"""

#read in raw Miseq R1 read file and output the sequence reads into a new txt file with the name formatted as sampleID_index.txt
#each line correspond to one sequence read in the raw file

import os
import re
import argparse
import numpy as np
import pandas as pd
import glob2

def generate_motif_map(lib, libpath):
    """
    For generate a compiled list of all motifs specified in lib.

    Parameters
    ----------
    lib: list
        List of Spec-seq libraries in the current run to be counted
        Default would be to search for all libraries: [M, Mrev]
    
    Returns
    ----------
    motif_lib: dictionary {str: [library, int]}
        A compiled dictionary with all sequence variants as keys, assigned library and counts as values. 
    
    """
    motif_lib = {}

    for element in lib:
        lib_map = f"{libpath}/{element}_library.txt"
       
        # Read in BC counts
        lib_map = pd.read_csv(lib_map, header=None, names=["sequence"],squeeze=True)
        # initialize a library with sequence variants as keys, assigned library and count=0 as values
        small_lib = {key: element for key in lib_map}
        #print(small_lib)
        # update the motif lib
        motif_lib = {**motif_lib, **small_lib}
    # make sure TAATTA - M in case Mrev overrides M
    motif_lib["TAATTA"] = "M"
    # convert motif library to dataframe
    motif_lib = pd.DataFrame.from_dict(motif_lib, orient='index', columns=['lib'])

    return motif_lib


def count_motifs(sample, input, lib, libpath):
    """
    For each input fasta fil, count all motifs in specified libraries and calculate the b/ub ratios.

    Parameters
    ----------
    input : list 
        List of fasta files containing raw sequencing reads for one sample/lane
    lib : list
        List of Spec-seq libraries in the current run to be counted
        Default would be to search for all libraries: [M, Mrev]

    Returns
    -------
    motif_count: pd dataframe
        A compiled dataframe containing (1)sequences (2)assigned library (3)d.count (4)m.count (5)f.count 
    """

    print("\t- Searching for matches in library\t" +' '.join(lib))
    
    # Initialize dictionary of all sequence variants specified in lib
    motiflib = generate_motif_map(lib, libpath)
    # Initialize dataframe to store count dataframes of a single sample/lane
    motif_counts = motiflib
    
    for file in input:
        # Retrieve the band name
        band_name = re.split(sample,file)[1][0]
        print("\t- processing file "+file)
        # Read in the tentative motif counts
        putative_motif_count_df = pd.read_csv(file, sep="\t", header=None, names=["count", "sequence"],
                                       index_col=1, squeeze=True)
                                       
        # Join together the two DataFrames to determine motif counts and representation
        motif_to_lib_df = pd.merge(motiflib, putative_motif_count_df, left_index=True, right_index=True)
        # Remove reads not in the library mapping
        motif_to_lib_df = motif_to_lib_df[motif_to_lib_df["lib"].notnull()]
        motif_to_lib_df = motif_to_lib_df.rename (columns = {"count" : band_name+".count"})

        # Compile counts from all bands into one dataframe
        motif_counts = pd.merge(motif_counts, motif_to_lib_df[[band_name+".count"]], left_index=True, right_index=True)
        # give the axis a name
        motif_counts = motif_counts.rename_axis("sequence")
    
    # print the number of reads with perfect match
    print("summary of reads with perfect match")
    for col in motif_counts.columns:
        if re.search("count",col):
            print(motif_counts[[col]].sum())
        
    return motif_counts

def cpm_motifs(count_df):
    """Helper function to perform cpm normalization of a read count matrix."""
    print("converting raw counts to cpm")
    for col in count_df.columns:
        if re.search("count",col):
            count_df[col] = count_df[col]/float(count_df[[col]].sum()) * 1000000
    return count_df

def find_rev_complement(query):
    """Return the reversed complement of query sequence."""
    complement={'A':'T','C':'G','G':'C','T':'A'}
    revcomp_query = ''.join(complement.get(base, base) for base in reversed(query))
    return revcomp_query


def check_consensus(lib, consensus):
    """Checking given consensus for length and nuleotide match to ACGT."""
    print("\t- Checking given consensus for length and nuleotide match to ACGT")
    reference_list = []
    default_consensus={"M":"TAATCC", "Mrev":"GGATTA"}
    # length of consensus should match that of lib
    expected_length={"M":6, "Mrev":6}
    for i in range(len(lib)):       
        # check given motif length
        if len(consensus[i]) != expected_length[lib[i]]:
            reference_list.append(default_consensus[lib[i]])
            print("\t- Length of library "+ lib[i] + " consensus does not match expected. Using default (CRX WT).")
        # check any nucleotide not matched ACGT
        elif bool(re.compile(r'[^ACGT]').search(consensus[i])) == True:
            reference_list.append(default_consensus[lib[i]])
            print("\t- Library "+ lib[i] + " consensus contains non-ACGT characters. Using default (CRX WT).")
        else:
            reference_list.append(consensus[i])
    print("\t- Finalized consensus " +' '.join(reference_list))
    return reference_list


def count_to_ratios(count_df, lib, consensus):
    """Check if the length of given consensus matches expected, if not, use default and print warning message."""
    reference_list=check_consensus(lib, consensus)

    # common fields for monomer and free bands
    # caculate ratios
    count_df['m/u'] = (count_df['m.count']/count_df['f.count'])
    # adjust for TAATTA (palindrome by itself)
    count_df.loc['TAATTA','m/u']=float(count_df.loc['TAATTA','m/u'])/2
    
    # create a new entry for TAATTA Mrev if both M and Mrev are in lib
    if "M" in lib and "Mrev" in lib:
        tmp_df =count_df.loc[['TAATTA']]
        tmp_df['lib']='Mrev'
        count_df=pd.concat([count_df,tmp_df]).sort_values(by='lib')
        del tmp_df

    # temporarily make the sequence index a column in order to apply calculation row wise
    count_df = count_df.reset_index()
    
    # calcualte normalized ratios for each consensus in reference_list
    for i in range(len(lib)):
        ref_RBE = count_df.loc[(count_df["sequence"]==reference_list[i]) & (count_df["lib"]==lib[i]), 'm/u']
        count_df[lib[i] + 'refm/u'] =  count_df['m/u']/float(ref_RBE)

    # fields only when dimer band also presnet
    if "d.count" in count_df.columns:
        # calculate ratios
        count_df['d/u'] = (count_df['d.count']/count_df['f.count'])
        # calcualte normalized ratios for each consensus in reference_list
        for i in range(len(lib)):
            ref_RBE = count_df.loc[(count_df["sequence"]==reference_list[i]) & (count_df["lib"]==lib[i]), 'd/u']
            count_df[lib[i] + 'refd/u'] =  count_df['d/u']/float(ref_RBE)
        # calculate du/mm
        count_df['du/mm']=count_df['d.count']*count_df['f.count']/(count_df['m.count']*count_df['m.count'])

    # rebuild the sequence index
    count_df = count_df.set_index("sequence")
    
    print("\t- All fields in motif ratio dataframe\n\t- " + '  '.join(count_df.columns))

    return reference_list, count_df;


def count_mismatches(query, lib, consensus):
    """Count and record mismatches of a list of query variants to the consensus."""
    # query should include at least two columns: sequence, lib
    # build a reference library
    if type(lib)==list and len(lib) > 1:
        reference_lib={k:v for k,v in zip(lib,consensus)}
    else:
        reference_lib={lib:consensus}
    print("using references:\n"+("\n\t-").join([f"lib {k}: {v}" for k,v in reference_lib.items()]))
    # initialize a new column
    query["MMcount"] = 0
    # now update MMcount column with mismatches to corresponding consensus
    query = query.reset_index(drop=False)
    for i in query.index:
        mismatch_num=sum(a!=b for a,b in zip(query.loc[i,"sequence"],reference_lib[query.loc[i,"lib"]]))
        query.at[i,"MMcount"]=mismatch_num
    # reset sequence as index column
    query = query.set_index("sequence")

    return query


def calculate_RBE(ratio_df, lib, consensus):
    """Calculate relative binding energy from count matrix."""
    # count mismatches to each reference sequence
    ratio_df = count_mismatches(ratio_df, lib, consensus)
    # convert normalized ratios to relative binding energy
    for k in ratio_df.columns:
        if re.search("ref", k):
            library = re.split("ref",k)[0]
            band = re.split("ref",k)[1][0]
            ratio_df[library+"."+band+"ddG"] = -np.log(ratio_df[k])

    print("\t- All fields in motif RBE dataframe\n\t- " + '  '.join(ratio_df.columns))

    return ratio_df


####################
##  main program  ##
####################

def main(sample, input, lib, ref, cpm, libpath, prefix):
    # set directory parameters
    os.chdir(prefix)
    if len(input) > 1:
        # count motif occurrence for all inputs
        motif_counts = count_motifs(sample, input, lib, libpath)
        # convert raw counts to cpm if specified
        if cpm:
            motif_counts = cpm_motifs(motif_counts)
        motif_counts.reset_index().rename(columns={"index" : "sequence"}).to_csv(prefix + "/counts/" + sample + "_count.txt", sep='\t', index=False)
        print("\t- Counts written to file: " + prefix + "/RBE/" + sample + "_count.txt")

        #print(motif_counts.head(4))
        
        # compile counts from all bands and calculate ratios
        updated_ref, motif_ratios = count_to_ratios(motif_counts, lib, ref)
        # convert normalized ratios to relative binding energy
        motif_RBE = calculate_RBE(motif_ratios, lib, updated_ref)
        # write to file
        motif_RBE.to_csv(prefix + "/RBE/" + sample + "_RBE.txt", sep='\t')
        print("\t- Relative binding energies written to file: " + prefix + "/RBE/" + sample + "_RBE.txt")
    else:
        # terminate the program if only one fasta file is provided
        raise ValueError("\t- Less than two input fasta files are provided. Counting terminated.")

if __name__ == "__main__":
    print("This is the start of calculate_ratios.py!!")
    # Setup argument parsing
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--sample", '-s', type=str, required=True,
        help="Names of sample/lane"
    )
    parser.add_argument(
        "--lib", "-l", nargs='*', choices=['M','Mrev'],
        default=['M','Mrev'],
        help="A list of Spec-seq libraries in the current run to be counted"
             "Default would be to search for all libraries: [M, Mrev]"
    )
    parser.add_argument(
        "--ref", "-r", nargs="*", default=["TAATCC", "GGATTA"],
        help="A list of consensus sequences to be used for calcualting relative binding energy"
             "Only count ratios will be calculated if reference are not specified"
    )
    parser.add_argument(
        "--cpm", action='store_true',
        help="Specify whether to convert raw counts to cpm values, defualt is to use raw counts"
    )
    parser.add_argument(
        "--libpath", type=str, default=glob2.glob("*/calculate_ratios.py".replace("/calculate_ratios.py","")),
        help="Path to directory storing all library maps."
             "Default would be to search in the same directory as this script."
    )
    parser.add_argument(
        "--prefix", "-p", type=str, default=os.getcwd(),
        help="The global prefix to give to all output text files"
    )
    args = parser.parse_args()

    # automatically search for all bands matched with sample/lane
    #print("Searching for input count files at " + args.prefix + "/counts/" + args.sample + "*.tmp")
    input = glob2.glob(os.path.join(args.prefix, "counts", args.sample + "*.tmp"))
    print("\t- All input count files :\n\t- " + '\n\t- '.join(input))

    # Pass arguments to main
    main(args.sample, input, args.lib,  args.ref, args.cpm, libpath=args.libpath, prefix=args.prefix)