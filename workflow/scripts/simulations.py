import numpy as np
import pandas as pd
from collections import defaultdict


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


def get_R(filepath):
    R = pd.DataFrame(np.zeros((4, 4)), index=["A", "C", "G", "T"], columns=["A", "C", "G", "T"])

    with open(filepath) as f:
        R_file = f.readlines()

    for line in R_file[1:]:
        line = line.strip().split("\t")
        R.loc[line[0].split("_")[1], line[0].split("_")[2]] = float(line[1])

    # Diagonal: 1 - sum of other elements
    for i in range(4):
        R.iloc[i, i] = 0 - R.iloc[i, :].sum()

    return R


# Determine sigma (vector) by this relation: sigma x R = 0 (null vector)
def get_steady_state(R_matrix):
    dimension = R_matrix.shape[0]
    M = np.vstack((R_matrix.transpose()[:-1], np.ones(dimension)))
    b = np.vstack((np.zeros((dimension - 1, 1)), [1]))
    return np.linalg.solve(M, b).transpose()[0]


def get_fitness(folderpath):
    df1 = pd.read_csv(f"{folderpath}/Euarchontoglires.Exclude.siteprofiles", sep="\t", header=0, index_col=None)
    df2 = pd.read_csv(f"{folderpath}/Euarchontoglires.Include.siteprofiles", sep="\t", header=0, index_col=None)

    # concatenate both files (keep column names)
    df = pd.concat([df1, df2], axis=0, ignore_index=True)
    df = df.iloc[:, 1:]

    # Take the col means, stored in a dictionary where key is the column name
    fitness = {aa: df[aa].mean() for aa in df.columns}

    return fitness


R = get_R("/Users/theo/THÉO/seascapes/data/Experiments/ENSG00000000003_TSPAN6_NT/sitemutsel_1.run.nucmatrix.tsv")
sigma = get_steady_state(R)
fitness = get_fitness("/Users/theo/THÉO/seascapes/processed/ENSG00000006715_VPS41")

print(R)
print(sigma)
print(fitness)




def get_Q(R, f):
    # All possible codons except stop codons
    codons = [c for c in codontable.keys() if c not in {"TAA", "TAG", "TGA"}]

    Q = pd.DataFrame(np.zeros((len(codons), len(codons))), index=codons, columns=codons)

    for i in codons:
        for j in codons:

            diff = [a+b for a, b in zip(i, j) if a != b]

            if len(diff) == 0 or len(diff) > 1:
                Q.loc[i, j] = 0

            elif codontable[i] == codontable[j]:
                Q.loc[i, j] = R.loc[diff[0][0], diff[0][1]]

            else:
                Q.loc[i, j] = R.loc[diff[0][0], diff[0][1]] * ((f[codontable[j]] - f[codontable[i]]) / (1 - np.exp(f[codontable[i]] - f[codontable[j]])))

    return Q


Q = get_Q(R, fitness)


def get_pi(sigma, f):
    return None


print(Q)

