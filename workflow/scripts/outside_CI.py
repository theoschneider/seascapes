import numpy as np
import pandas as pd
from libraries import getconfint, open_ali, write_fasta

folder = "/Users/theo/THÉO/seascapes/"

# Read the data
sequences = open_ali(folder + "processed/ENSG00000006715_VPS41/Euarchontoglires.Whole.phy")
distance = pd.read_csv(folder + "processed/ENSG00000006715_VPS41/distance.tsv", sep="\t", header=0, index_col=None)
omega = pd.read_csv(folder + "processed/ENSG00000006715_VPS41/Euarchontoglires.Whole.omega.ci0.025.tsv", sep="\t", header=0, index_col=0)

# Run the linear regression
lower_ci, upper_ci = getconfint(x=omega.iloc[1:, 1].tolist(), y=distance.iloc[:, 1].tolist(), alpha=0.05)

sequences["conf_filter"] = ""

for position in range(len(distance.iloc[:, 1])):

    if (distance.iloc[position, 1] < lower_ci[position]) or (distance.iloc[position, 1] > upper_ci[position]):
        sequences["conf_filter"] += "AAA"
    else:
        sequences["conf_filter"] += "TTT"

write_fasta(sequences, folder + "processed/ENSG00000006715_VPS41/Euarchontoglires.Whole.conf_filter.fasta")
