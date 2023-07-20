import os
import gzip
import shutil
from collections import defaultdict
import numpy as np
import pandas as pd

# Define the codon table
codontable = defaultdict(lambda: "-")
codontable.update({
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'})


def open_fasta(path, filter_species=None):
    if filter_species is None:
        filter_species = []
    outfile = {}
    ali_file = gzip.open(path, 'rt') if path.endswith(".gz") else open(path, 'r')
    for seq_id in ali_file:
        sp_id = seq_id.replace('>', '').strip()
        seq = ali_file.readline().strip()
        if (len(filter_species) > 0) and (sp_id not in filter_species):
            continue
        outfile[sp_id] = seq
    print(f"Read {len(outfile)} sequences from {path}")
    return outfile


def write_fasta(dico_fasta, output):
    outfile = gzip.open(output, 'wt') if output.endswith(".gz") else open(output, 'w')
    outfile.write("\n".join([f">{seq_id}\n{seq}" for seq_id, seq in dico_fasta.items()]))
    outfile.close()
    print(f"Written {len(dico_fasta)} sequences to {output}")


def filter_fasta(dico_fasta1, dico_fasta2):
    n = len(list(dico_fasta1.values())[0])

    bool_filter = []

    for nt in range(0, n, 3):
        triplets1 = [seq[nt:nt + 3] for seq in dico_fasta1.values()]
        triplets2 = [seq[nt:nt + 3] for seq in dico_fasta2.values()]

        if (triplets1.count("---") / len(triplets1)) > 0.5:
            bool_filter.extend([False] * 3)

        elif (triplets2.count("---") / len(triplets2)) > 0.5:
            bool_filter.extend([False] * 3)

        else:
            bool_filter.extend([True] * 3)

    # Filter the sequences
    for seq_id in dico_fasta1.keys():
        dico_fasta1[seq_id] = "".join([dico_fasta1[seq_id][i] for i in range(len(bool_filter)) if bool_filter[i]])
    for seq_id in dico_fasta2.keys():
        dico_fasta2[seq_id] = "".join([dico_fasta2[seq_id][i] for i in range(len(bool_filter)) if bool_filter[i]])

    return dico_fasta1, dico_fasta2


def write_ali(dico_fasta, output):
    outfile = gzip.open(output, 'wt') if output.endswith(".gz") else open(output, 'w')
    # Assert that all sequences have the same length
    set_seq_len = set([len(seq) for seq in dico_fasta.values()])
    assert len(set_seq_len) == 1, "All sequences must have the same length"
    seq_len = set_seq_len.pop()

    # Write the header (the number of sequences and the sequence length)
    outfile.write(f"{len(dico_fasta)} {seq_len}\n")
    # Write the names and sequences
    outfile.write("\n".join([f"{name} {seq}" for name, seq in dico_fasta.items()]))
    outfile.close()
    print(f"Written {len(dico_fasta)} sequences to {output}")


def zip_file(input_path, output_path=None):
    if not os.path.isfile(input_path):
        return
    if output_path is None:
        output_path = input_path + ".gz"
    with open(input_path, 'rb') as f_in:
        with gzip.open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(input_path)


def translate_cds(seq):
    return "".join([codontable[seq[i * 3:i * 3 + 3]] for i in range(len(seq) // 3)])


def read_profiles(path):
    return pd.read_csv(path, sep="\t", skiprows=1, header=None, names="site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y".split(","))


def js(dataframe1: pd.DataFrame, dataframe2: pd.DataFrame) -> pd.DataFrame:

    # Assert that both dataframes contain 20 columns
    assert dataframe1.shape[1] == 20 and dataframe2.shape[1] == 20, "Both dataframes must contain 20 columns"

    # Assert that both dataframes contain the same number of rows
    assert dataframe1.shape[0] == dataframe2.shape[0], "Both dataframes must contain the same number of rows"

    # Assert that each row sum up to 1
    assert np.allclose(dataframe1.sum(axis=1), 1), "Each row of dataframe 1 must sum up to 1"
    assert np.allclose(dataframe2.sum(axis=1), 1), "Each row of dataframe 2 must sum up to 1"

    # Average of the two profiles
    m = (dataframe1 + dataframe2) / 2

    # Compute the Kullback-Leibler divergence
    kl1 = np.sum(dataframe1 * np.log(dataframe1 / m), axis=1)
    kl2 = np.sum(dataframe2 * np.log(dataframe2 / m), axis=1)

    return (kl1 + kl2) / 2
