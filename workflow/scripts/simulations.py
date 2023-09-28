import numpy as np
import pandas as pd
from collections import defaultdict, Counter


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

# All possible codons except stop codons
codons = [c for c in codontable.keys() if c not in {"TAA", "TAG", "TGA"}]


def get_R(filepath):
    R = pd.DataFrame(np.zeros((4, 4)), index=["A", "C", "G", "T"], columns=["A", "C", "G", "T"])

    with open(filepath) as f:
        R_file = f.readlines()

    for line in R_file[1:]:
        line = line.strip().split("\t")
        R.loc[line[0].split("_")[1], line[0].split("_")[2]] = float(line[1])

    # Diagonal: 1 - sum of other elements
    for i in range(4):
        R.iloc[i, i] = -R.iloc[i, :].sum()

    return R


# Determine sigma (vector) by this relation: sigma x R = 0 (null vector), using linear algebra
def get_steady_state(R_matrix):
    dimension = R_matrix.shape[0]
    M = np.vstack((R_matrix.transpose()[:-1], np.ones(dimension)))
    b = np.vstack((np.zeros((dimension - 1, 1)), [1]))
    solution = np.linalg.solve(M, b).transpose()[0]
    # Save solution in a dictionary where keys are the matrix column names
    sigma_dict = {R_matrix.columns[i]: solution[i] for i in range(dimension)}
    return sigma_dict


def get_fitness(filepath):
    df = pd.read_csv(filepath, sep="\t", header=0, index_col=None)

    # Take the first row of the df, save it in a dictionary where key is col name
    fitness = {aa: np.log(df[aa][0]) for aa in df.columns[1:]}

    return fitness


def get_Q(R, f):

    Q = pd.DataFrame(np.zeros((len(codons), len(codons))), index=codons, columns=codons)

    for i in codons:
        for j in codons:

            diff = [a+b for a, b in zip(i, j) if a != b]

            if len(diff) == 0 or len(diff) > 1:
                continue

            Ai = codontable[i]
            Aj = codontable[j]

            if Ai == Aj:
                Q.loc[i, j] = R.loc[diff[0][0], diff[0][1]]

            else:
                s = f[Aj] - f[Ai]

                if abs(s) < 1e-6:
                    pfix = 1 + s / 2
                else:
                    pfix = s / (1 - np.exp(-s))

                Q.loc[i, j] = R.loc[diff[0][0], diff[0][1]] * pfix

    # Diagonal: 1 - sum of other elements
    for i in range(len(codons)):
        Q.iloc[i, i] = -Q.iloc[i, :].sum()

    return Q


def get_pi(sigma, f):
    pi = {}

    for codon in codons:
        pi[codon] = (sigma[codon[0]] * sigma[codon[1]] * sigma[codon[2]] * np.exp(f[codontable[codon]]))

    denom = sum(pi.values())
    pi = {k: v / denom for k, v in pi.items()}

    return pi


R = get_R("/Users/theo/THÉO/seascapes/data/Experiments/ENSG00000000003_TSPAN6_NT/sitemutsel_1.run.nucmatrix.tsv")
sigma = get_steady_state(R)
fitness = get_fitness("/Users/theo/THÉO/seascapes/processed/ENSG00000006715_VPS41/Euarchontoglires.Exclude.siteprofiles")
Q = get_Q(R, fitness)
pi = get_pi(sigma, fitness)


print(R)
print(sigma)

print(fitness)

print(Q)
print(pi)


def run_simulation(Q, max_time=10000, seed=None):

    if seed is not None:
        np.random.seed(seed)

    dimension = Q.shape[0]
    state, clock = 0, 0
    history = defaultdict(lambda: 0)

    while clock < max_time:

        row = Q.iloc[state, :]
        potential_states = np.where(row > 0)
        potential_states = potential_states[0]
        rates = row[potential_states]
        samples = np.random.exponential(1 / rates)
        time = np.min(samples)
        clock += time

        history[state] += time

        state = potential_states[np.argmin(samples)]

    return Counter(history)


simulation = run_simulation(Q, seed=42)


print(simulation)
