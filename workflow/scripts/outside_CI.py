import numpy as np
import pandas as pd
from libraries import getconfint, open_fasta

folder = "/Users/theo/THÉO/seascapes/"

# Read the data
fasta = open_fasta(folder + "data/omm_NT_fasta.v10c_116/ENSG00000150401_DCUN1D2_NT.fasta")
distance = pd.read_csv(folder + "processed/ENSG00000150401_DCUN1D2/distance.tsv", sep="\t", header=0, index_col=None)
omega = pd.read_csv(folder + "processed/ENSG00000150401_DCUN1D2/Euarchontoglires.Whole.omega.ci0.025.tsv", sep="\t", header=0, index_col=0)

# Run the linear regression
lower_ci, upper_ci = getconfint(x=omega.iloc[1:, 1].tolist(), y=distance.iloc[:, 1].tolist(), alpha=0.05)

for position in range(len(distance.iloc[:, 1])):
    print(distance.iloc[position, 1], omega.iloc[position + 1, 1], lower_ci[position], upper_ci[position])

