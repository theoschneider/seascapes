import numpy as np
import pandas as pd
from collections import defaultdict
from ete3 import Tree


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


def get_R(filepath: str) -> pd.DataFrame:
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
def get_steady_state(R_matrix: pd.DataFrame) -> dict:
    dimension = R_matrix.shape[0]
    M = np.vstack((R_matrix.transpose()[:-1], np.ones(dimension)))
    b = np.vstack((np.zeros((dimension - 1, 1)), [1]))
    solution = np.linalg.solve(M, b).transpose()[0]
    # Save solution in a dictionary where keys are the matrix column names
    sigma_dict = {R_matrix.columns[i]: solution[i] for i in range(dimension)}
    return sigma_dict


def get_Q(R: pd.DataFrame, f: dict) -> pd.DataFrame:

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


def get_pi(sigma: dict, f: dict) -> list:
    pi = []

    for codon in codons:
        pi.append(sigma[codon[0]] * sigma[codon[1]] * sigma[codon[2]] * np.exp(f[codontable[codon]]))

    pi = pi / np.sum(pi)

    return pi


def run_simulation(Q: pd.DataFrame, initial_state: int, max_time: int, seed=None):

    if seed is not None:
        np.random.seed(seed)

    dimension = Q.shape[0]
    clock = 0
    history = defaultdict(lambda: 0)
    state = initial_state

    codon_list = [Q.columns[state]]
    time_list = []

    while True:

        row = Q.iloc[state, :]
        potential_states = np.where(row > 0)[0]
        rates = row[potential_states]
        samples = np.random.exponential(1 / rates)
        time = np.min(samples)
        time_list.append(clock)
        clock += time

        if clock > max_time:
            break

        history[state] += time

        state = potential_states[np.argmin(samples)]
        codon_list.append(Q.columns[state])

        # assert difference is only 1 mutation
        assert len([a+b for a, b in zip(codon_list[-1], codon_list[-2]) if a != b]) == 1, "More than 1 mutation"

    if len(time_list) > 1:
        assert time_list[-1] < max_time, "Clock exceeded max_time"

    return codon_list, time_list, state


folder = "/Users/theo/THÉO/seascapes"
gene = "ENSG00000006715_VPS41"

tree = Tree(f"{folder}/data/omm_RooTree.v10b_116/{gene}_NT.rootree", format=1)
tree_depth = sum([node.dist for node in tree.traverse()])


def get_seq(tree: Tree, fitness_path: str, R: str, filename: str, seed=None) -> None:

    # calculate and print tree total depth
    tree_depth = sum([node.dist for node in tree.traverse()])
    print(f"Tree depth: {tree_depth}")

    if seed is not None:
        np.random.seed(seed)

    R = get_R(R)
    sigma = get_steady_state(R)
    fitness_df = pd.read_csv(fitness_path, sep="\t", header=0, index_col=None)
    seq_length = fitness_df.shape[0]

    for node in tree.traverse("preorder"):
        node.add_feature("states", [None] * seq_length)

    for i in range(seq_length):

        fit_pos = {aa: np.log(fitness_df[aa][i]) for aa in fitness_df.columns[1:]}
        Q = get_Q(R, fit_pos)
        pi = get_pi(sigma, fit_pos)

        for node in tree.traverse("preorder"):
            if node.is_root():
                state = np.random.choice(range(len(codons)), p=pi)

            else:
                initial_state = node.up.states[i]
                assert initial_state is not None, "initial state is None"
                codon_list, _, state = run_simulation(Q, initial_state, max_time=node.dist)

            node.states[i] = state

    with open(filename, "w") as f:
        for node in tree.traverse("preorder"):
            if node.is_leaf():
                f.write(">" + node.name + "\n")
                f.write("".join([codons[j] for j in node.states]) + "\n")

    return None


get_seq(tree=tree,
        fitness_path=f"{folder}/processed/{gene}/Euarchontoglires.Exclude.siteprofiles",
        R=f"{folder}/data/Experiments/{gene}_NT/sitemutsel_1.run.nucmatrix.tsv",
        filename=f"{folder}/results/sim_seqs.fasta",
        seed=42)

