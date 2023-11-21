import os
import pandas as pd
import matplotlib.pyplot as plt
from libraries import linearfit
import numpy as np

plt.rcParams['text.usetex'] = True

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
    and os.path.exists(os.path.join(folder_path, CLADE + ".Whole.omega_0.ci0.025.tsv"))
    and os.path.exists(os.path.join(folder_path, CLADE + ".Whole.mask.tsv"))):
        print(f"Skipping folder: {folder_path}")
        continue

    # For every folder
    print(f"Processing folder: {folder_path}")

    name = folder_name.split("_")[1]

    # Open the mask
    mask = pd.read_csv(os.path.join(folder_path, CLADE + ".Whole.mask.tsv"), sep="\t", header=0, index_col=None)

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

    # Filter the distance and omega lists
    distance = [distance[i] for i in range(len(distance)) if mask.iloc[i, 1] < 0.95]
    omega = [omega[i] for i in range(len(omega)) if mask.iloc[i, 1] < 0.95]
    omega0 = [omega0[i] for i in range(len(omega0)) if mask.iloc[i, 1] < 0.95]

    # Assert that the length of omega and omega0 are the same
    assert len(omega) == len(omega0) == len(distance), "Length of omega, omega0 and distance are not the same"

    for i in range(len(omega)):
        # Add all values to the df
        df = df._append(
            {"distance": distance[i], "omega0": omega0[i], "omegaA": omega[i] - omega0[i]},
            ignore_index=True)

print("Done with concatenating all the dataframes")

# Sort the df and calculate the step
df.sort_values(by=["omega0"], inplace=True)
step = len(df) // n_bins

# Iterate over the number of bins
for bin in range(n_bins):

    lower = df.iloc[bin*step, 1]

    if bin == n_bins - 1:
        # If it is the last bin, take the last value of omega0 as upper bound
        upper = max(df["omega0"])
    else:
        # Take the last value of the bin
        upper = df.iloc[(bin+1)*step, 1]

    df_bin = df[(df["omega0"] > lower) & (df["omega0"] <= upper)]

    # Plot distance vs omega and save it
    yfit, lower_ci, upper_ci = linearfit(df_bin["omegaA"], df_bin["distance"], np.linspace(min(df_bin["omegaA"]), max(df_bin["omegaA"]), 100))
    plt.figure(figsize=(7, 7))
    plt.scatter(df_bin["omegaA"], df_bin["distance"], s=2)
    plt.fill_between(x=np.linspace(min(df_bin["omegaA"]), max(df_bin["omegaA"]), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
    plt.plot(np.linspace(min(df_bin["omegaA"]), max(df_bin["omegaA"]), 100), yfit, linewidth=2, color=lin_reg)
    plt.xlabel("$\omega_A$")
    plt.ylabel("Distance")
    plt.title(f"Distance as a function of $\omega_A$ for bin {bin+1}: {lower:.2f} $< \omega_0 \le$ {upper:.2f}")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(f"{SOURCE_DIR}/results/distance-vs-omegaA_bin-{bin+1}.pdf")

