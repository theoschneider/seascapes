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


def kl(path1: str, path2: str):
    # Compute Kullback-Leibler divergence between two profiles
    df1 = pd.read_csv(path1, sep="\t", header=0)
    df2 = pd.read_csv(path2, sep="\t", header=0)

    kl1 = np.sum(df1 * np.log(df1 / df2), axis=1)
    kl2 = np.sum(df2 * np.log(df2 / df1), axis=1)

    return (kl1 + kl2) / 2
