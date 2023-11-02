import os
import pandas as pd
import matplotlib.pyplot as plt
from libraries import linearfit
import numpy as np


# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Define the gene
CDS = "ENSG00000040199_PHLPP2"

# Get folder with processed data
folder_path = os.path.join(SOURCE_DIR, "processed", CDS)

print(f"Processing folder: {folder_path}")
name = CDS.split("_")[1]

# Open distance.tsv
distance = pd.read_csv(os.path.join(folder_path, "distance.tsv"), sep="\t", header=0, index_col=None)

# Open omega
omega = pd.read_csv(os.path.join(folder_path, "Euarchontoglires.Whole.omega.ci0.025.tsv"), sep="\t", header=0, index_col=0)

# Keep only the col of interest and convert to lists
distance = distance.iloc[:, 1].tolist()
omega = omega.iloc[1:, 1].tolist()

# Open the mask
mask = pd.read_csv(os.path.join(folder_path, "Euarchontoglires.Whole.mask.tsv"), sep="\t", header=0, index_col=None)

# Filter the distance and omega lists
distance = [distance[i] for i in range(len(distance)) if mask.iloc[i, 1] < 0.95]
omega = [omega[i] for i in range(len(omega)) if mask.iloc[i, 1] < 0.95]

# Initialize colors
fill_reg = [0, 0, 0, 0.1]
lin_reg = [0, 0, 0, 1]

# Plot distance vs omega and save it
yfit, lower_ci, upper_ci = linearfit(omega, distance, np.linspace(min(omega), max(omega), 100))
plt.figure(figsize=(7, 7))
plt.scatter(omega, distance, s=2)
plt.fill_between(x=np.linspace(min(omega), max(omega), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
plt.plot(np.linspace(min(omega), max(omega), 100), yfit, linewidth=2, color=lin_reg)
plt.xlabel("Omega")
plt.ylabel("Distance")
plt.title("Distance as a function of omega for gene " + name + " (only variable sites)")
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/" + CDS + "/filtered_distance-vs-omega.pdf")
