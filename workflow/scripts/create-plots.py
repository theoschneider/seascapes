import os
import pandas as pd
import matplotlib.pyplot as plt
from libraries import linearfit
import numpy as np

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Define the range (FROM and TO), the flag to process all folders (ALL) and the names of the genes (if needed)
FROM = 1
TO = 5
ALL = False
GENES = ["ENSG00000204120_GIGYF2",
         "ENSG00000150401_DCUN1D2",
         "ENSG00000151364_KCTD14",
         "ENSG00000110583_NAA40",
         "ENSG00000204536_CCHCR1",
         "ENSG00000198604_BAZ1A",
         "ENSG00000154146_NRGN",
         "ENSG00000108582_CPD",
         "ENSG00000163349_HIPK1",
         "ENSG00000111540_RAB5B",
         "ENSG00000147432_CHRNB3",
         "ENSG00000011198_ABHD5"]
CLADE = "Euarchontoglires"

# Check if ALL is true, then adjust FROM and TO to process all folders
if ALL:
    FROM = 1
    TO = len(os.listdir(SOURCE_DIR))

# Check if some genes have been specified
if GENES:
    genes_list = GENES
else:
    genes_list = os.listdir(os.path.join(SOURCE_DIR, "processed"))[FROM:TO + 1]

# Initialize empty dataframes for all genes
all_omega = []
all_omega0 = []
all_distance = []

# Iterate over the contents of SOURCE_DIR and process folders
for folder_name in genes_list:
    # For every folder
    folder_path = os.path.join(SOURCE_DIR, "processed", folder_name)
    print(f"Processing folder: {folder_path}")

    # Open distance.tsv
    distance = pd.read_csv(os.path.join(folder_path, "distance.tsv"), sep="\t", header=0, index_col=None)

    # Append distance to the all_distance dataframe
    all_distance.extend(distance.iloc[:, 1].tolist())

    # Open estimated omega and omega0, compute omegaA
    omega = pd.read_csv(SOURCE_DIR + "/processed/" + folder_name + "/" + CLADE + ".Whole.omega.ci0.025.tsv", sep="\t", header=0, index_col=0)
    omega0 = pd.read_csv(SOURCE_DIR + "/processed/" + folder_name + "/" + CLADE + ".Whole.omega_0.ci0.025.tsv", sep="\t", header=0, index_col=0)

    omega = omega.iloc[1:, 1].tolist()
    omega0 = omega0.iloc[1:, 1].tolist()
    omegaA = [omega[i] - omega0[i] for i in range(len(omega))]

    # Plot the distance and omega on the same plot
    fig, ax1 = plt.subplots()
    col1 = "#2674D9"
    col2 = "#D92674"
    width = len(distance) / 20
    fig.set_size_inches(width, 7)
    ax1.set_xlabel("Position")
    ax1.set_ylabel("Distance", color=col1)
    ax1.bar(x=distance.iloc[:, 0]-0.25, height=distance.iloc[:, 1], width=0.5, color=col1, alpha=1)
    ax1.tick_params(axis='y', labelcolor=col1)
    ax2 = ax1.twinx()  # instantiate a second axis that shares the same x-axis
    ax2.set_ylabel("omega", color=col2)
    ax2.bar(x=distance.iloc[:, 0]+0.25, height=omega, width=0.5, color=col2, alpha=1)
    ax2.tick_params(axis='y', labelcolor=col2)
    plt.title("Distance and omega per site")
    fig.tight_layout()
    os.makedirs(os.path.join(SOURCE_DIR, "results", folder_name), exist_ok=True)
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/distance.pdf")

    # Append omega and omega0 to the all_omega and all_omega0 dataframes
    all_omega.extend(omega)
    all_omega0.extend(omega0)

    # Plot distance vs omega and save it
    yfit, lower_ci, upper_ci = linearfit(omega, distance.iloc[:, 1].tolist(), np.linspace(min(omega), max(omega), 100))
    plt.figure(figsize=(7, 7))
    plt.scatter(omega, distance.iloc[:, 1], s=2)
    plt.fill_between(x=np.linspace(min(omega), max(omega), 100), y1=lower_ci, y2=upper_ci, color=[1, 0, 0, 0.15], edgecolor=None)
    plt.plot(np.linspace(min(omega), max(omega), 100), yfit, linewidth=1, color=[1, 0, 0, .8])
    plt.xlabel("Omega")
    plt.ylabel("Distance")
    plt.title("Distance vs omega")
    plt.tight_layout()
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/distance-vs-omega.pdf")

    # Plot distance vs omega0 and save it
    yfit, lower_ci, upper_ci = linearfit(omega0, distance.iloc[:, 1].tolist(), np.linspace(min(omega0), max(omega0), 100))
    plt.figure(figsize=(7, 7))
    plt.scatter(omega0, distance.iloc[:, 1], s=2)
    plt.fill_between(x=np.linspace(min(omega0), max(omega0), 100), y1=lower_ci, y2=upper_ci, color=[1, 0, 0, 0.15], edgecolor=None)
    plt.plot(np.linspace(min(omega0), max(omega0), 100), yfit, linewidth=1, color=[1, 0, 0, .8])
    plt.xlabel("Omega 0")
    plt.ylabel("Distance")
    plt.title("Distance vs omega 0")
    plt.tight_layout()
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/distance-vs-omega0.pdf")

    # Plot distance vs omegaA and save it
    yfit, lower_ci, upper_ci = linearfit(omegaA, distance.iloc[:, 1].tolist(), np.linspace(min(omegaA), max(omegaA), 100))
    plt.figure(figsize=(7, 7))
    plt.scatter(omegaA, distance.iloc[:, 1], s=2)
    plt.fill_between(x=np.linspace(min(omegaA), max(omegaA), 100), y1=lower_ci, y2=upper_ci, color=[1, 0, 0, 0.15], edgecolor=None)
    plt.plot(np.linspace(min(omegaA), max(omegaA), 100), yfit, linewidth=1, color=[1, 0, 0, .8])
    plt.xlabel("OmegaA")
    plt.ylabel("Distance")
    plt.title("Distance vs omegaA (omega - omega0)")
    plt.tight_layout()
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/distance-vs-omegaA.pdf")


all_omegaA = [all_omega[i] - all_omega0[i] for i in range(len(all_omega))]

# Plot distance vs total omega and save it
yfit, lower_ci, upper_ci = linearfit(all_omega, all_distance, np.linspace(min(all_omega), max(all_omega), 100))
plt.figure(figsize=(7, 7))
plt.scatter(all_omega, all_distance, s=2)
plt.fill_between(x=np.linspace(min(all_omega), max(all_omega), 100), y1=lower_ci, y2=upper_ci, color=[1, 0, 0, 0.15], edgecolor=None)
plt.plot(np.linspace(min(all_omega), max(all_omega), 100), yfit, linewidth=1, color=[1, 0, 0, .8])
plt.xlabel("Omega")
plt.ylabel("Distance")
plt.title("Distance vs omega")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/allgenes-distance-vs-omega.pdf")

# Plot distance vs total omega0 and save it
yfit, lower_ci, upper_ci = linearfit(all_omega0, all_distance, np.linspace(min(all_omega0), max(all_omega0), 100))
plt.figure(figsize=(7, 7))
plt.scatter(all_omega0, all_distance, s=2)
plt.fill_between(x=np.linspace(min(all_omega0), max(all_omega0), 100), y1=lower_ci, y2=upper_ci, color=[1, 0, 0, 0.15], edgecolor=None)
plt.plot(np.linspace(min(all_omega0), max(all_omega0), 100), yfit, linewidth=1, color=[1, 0, 0, .8])
plt.xlabel("Omega 0")
plt.ylabel("Distance")
plt.title("Distance vs omega 0")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/allgenes-distance-vs-omega0.pdf")

# Plot distance vs total omegaA and save it
yfit, lower_ci, upper_ci = linearfit(all_omegaA, all_distance, np.linspace(min(all_omegaA), max(all_omegaA), 100))
plt.figure(figsize=(7, 7))
plt.scatter(all_omegaA, all_distance, s=2)
plt.fill_between(x=np.linspace(min(all_omegaA), max(all_omegaA), 100), y1=lower_ci, y2=upper_ci, color=[1, 0, 0, 0.15], edgecolor=None)
plt.plot(np.linspace(min(all_omegaA), max(all_omegaA), 100), yfit, linewidth=1, color=[1, 0, 0, .8])
plt.xlabel("OmegaA")
plt.ylabel("Distance")
plt.title("Distance vs omegaA (omega - omega0)")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/allgenes-distance-vs-omegaA.pdf")
