import os
import pandas as pd
import matplotlib.pyplot as plt
from libraries import linearfit
import numpy as np

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Define the range (FROM and TO), the flag to process all folders (ALL) and the names of the genes (if needed)
CLADE = "Euarchontoglires"

# Get the list of all genes, remove DS_Store if present
all_genes = os.listdir(os.path.join(SOURCE_DIR, "processed"))
if ".DS_Store" in all_genes:
    all_genes.remove(".DS_Store")

# Initialise an empty df with distance, omega0, omegaA
df = pd.DataFrame(columns=["distance", "omega0", "omegaA"])

# Initialize colors
fill_reg = [0, 0, 0, 0.1]
lin_reg = [0, 0, 0, 1]

# Initialize the number of bins
n_bins = 8

# Iterate over the contents of SOURCE_DIR and process folders
for folder_name in all_genes:

    folder_path = os.path.join(SOURCE_DIR, "processed", folder_name)

    # Check if the folder contains distance.tsv, omega, omega0
    if not (os.path.exists(os.path.join(folder_path, "distance.tsv"))
    and os.path.exists(os.path.join(folder_path, CLADE + ".Whole.omega.ci0.025.tsv"))
    and os.path.exists(os.path.join(folder_path, CLADE + ".Whole.omega_0.ci0.025.tsv"))):
        continue

    # For every folder
    print(f"Processing folder: {folder_path}")

    name = folder_name.split("_")[1]

    # Open distance.tsv
    distance = pd.read_csv(os.path.join(folder_path, "distance.tsv"), sep="\t", header=0, index_col=None)

    # Open omega
    omega = pd.read_csv(os.path.join(folder_path, CLADE + ".Whole.omega.ci0.025.tsv"), sep="\t", header=0, index_col=0)

    # Open omega0
    omega0 = pd.read_csv(os.path.join(folder_path, CLADE + ".Whole.omega_0.ci0.025.tsv"), sep="\t", header=0, index_col=0)

    # Keep only the col of interest and convert to lists
    distance = distance.iloc[:, 1].tolist()
    omega = omega.iloc[1:, 1].tolist()
    omega0 = omega0.iloc[1:, 1].tolist()

    # Assert that the length of omega and omega0 are the same
    assert len(omega) == len(omega0) == len(distance), "Length of omega, omega0 and distance are not the same"

    for i in range(len(omega)):
        # Add all values to the df
        df = df._append(
            {"distance": distance[i], "omega0": omega0[i], "omegaA": omega[i] - omega0[i]},
            ignore_index=True)


# Plot distance vs omega and save it (for every bin)
for bin in range(n_bins):

    lower = min(df["omega0"] + bin / n_bins * (max(df["omega0"]) - min(df["omega0"])))
    upper = min(df["omega0"] + (bin + 1) / n_bins * (max(df["omega0"]) - min(df["omega0"])))

    # Filter the df
    df_bin = df[(df["omega0"] > lower) & (df["omega0"] <= upper)]

    # Plot distance vs omega and save it
    yfit, lower_ci, upper_ci = linearfit(df_bin["omegaA"], df_bin["distance"], np.linspace(min(df_bin["omegaA"]), max(df_bin["omegaA"]), 100))
    plt.figure(figsize=(7, 7))
    plt.scatter(df_bin["omegaA"], df_bin["distance"], s=2)
    plt.fill_between(x=np.linspace(min(df_bin["omegaA"]), max(df_bin["omegaA"]), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
    plt.plot(np.linspace(min(df_bin["omegaA"]), max(df_bin["omegaA"]), 100), yfit, linewidth=2, color=lin_reg)
    plt.xlabel("Omega")
    plt.ylabel("Distance")
    plt.title("Distance as a function of omega for bin " + str(bin))
    plt.xscale("log")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(SOURCE_DIR + "/results/distance-vs-omega_bin-" + str(bin+1) + ".pdf")

