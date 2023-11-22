import os
import pandas as pd
import matplotlib.pyplot as plt
from libraries import linearfit, set_size
import numpy as np

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Define the range (FROM and TO), the flag to process all folders (ALL) and the names of the genes (if needed)
FROM = 0
TO = 5
ALL = True
GENES = []
CLADE = "Euarchontoglires"

# Get the list of all genes, remove DS_Store if present
all_genes = os.listdir(os.path.join(SOURCE_DIR, "processed"))
if ".DS_Store" in all_genes:
    all_genes.remove(".DS_Store")

# Check if ALL is true, then adjust FROM and TO to process all folders
if ALL:
    FROM = 1
    TO = len(all_genes)

# Check if some genes have been specified
if GENES:
    genes_list = GENES
else:
    genes_list = all_genes[FROM:TO]


# Initialize empty dataframes for all genes
all_omega = []
all_omega0 = []
all_distance = []

# Initialize empty lists for 1 point per gene
omegaA_per_gene = []
omega_per_gene = []
distance_per_gene = []
all_names = []

# Initialize colors
fill_reg = [0, 0, 0, 0.1]
lin_reg = [0, 0, 0, 1]

# Iterate over the contents of SOURCE_DIR and process folders
for folder_name in genes_list:
    # For every folder
    folder_path = os.path.join(SOURCE_DIR, "processed", folder_name)

    # Check if the folder contains distance.tsv, omega, omega0
    if not (os.path.exists(os.path.join(folder_path, "distance.tsv"))
    and os.path.exists(os.path.join(folder_path, CLADE + ".Whole.omega.ci0.025.tsv"))
    and os.path.exists(os.path.join(folder_path, CLADE + ".Whole.omega_0.ci0.025.tsv"))
    and os.path.exists(os.path.join(folder_path, CLADE + ".Whole.mask.tsv"))):
        print(f"Skipping folder: {folder_path}")
        continue

    print(f"Processing folder: {folder_path}")
    name = folder_name.split("_")[1]
    all_names.append(name)

    # Check if the result folder exists, if not create it
    if not os.path.exists(os.path.join(SOURCE_DIR, "results", folder_name)):
        os.makedirs(os.path.join(SOURCE_DIR, "results", folder_name))

    # Open the mask
    mask = pd.read_csv(os.path.join(folder_path, CLADE + ".Whole.mask.tsv"), sep="\t", header=0, index_col=None)
    mask = mask.iloc[:, 1]

    # Open distance.tsv
    distance_df = pd.read_csv(os.path.join(folder_path, "distance.tsv"), sep="\t", header=0, index_col=None)

    # Filter distance but keep both columns
    filt_distance_df = distance_df[mask < 0.95]

    # Keep only the col of interest and convert to list
    distance = distance_df.iloc[:, 1].tolist()
    positions = distance_df.iloc[:, 0].tolist()
    filt_distance = filt_distance_df.iloc[:, 1].tolist()
    filt_positions = filt_distance_df.iloc[:, 0].tolist()

    # Open estimated omega and omega0, compute omegaA
    omega = pd.read_csv(SOURCE_DIR + "/processed/" + folder_name + "/" + CLADE + ".Whole.omega.ci0.025.tsv", sep="\t", header=0, index_col=0)
    omega0 = pd.read_csv(SOURCE_DIR + "/processed/" + folder_name + "/" + CLADE + ".Whole.omega_0.ci0.025.tsv", sep="\t", header=0, index_col=0)

    # Append distance to the all_distance dataframe
    all_distance.extend(distance)

    # Append omega and omega0 to the all_omega and all_omega0 dataframes
    all_omega.extend(omega)
    all_omega0.extend(omega0)

    # Append distance to the distance_per_gene list
    distance_per_gene.append(np.mean(distance))

    # Append omegaA to the omegaA_per_gene list
    omegaA_per_gene.append(omega.iloc[0, 1] - omega0.iloc[0, 1])
    omega_per_gene.append(omega.iloc[0, 1])

    # Keep only the col of interest and convert to lists
    omega = omega.iloc[1:, 1].tolist()
    omega0 = omega0.iloc[1:, 1].tolist()
    omegaA = [omega[i] - omega0[i] for i in range(len(omega))]
    filt_omega = [omega[i] for i in range(len(omega)) if mask[i] < 0.95]
    filt_omega0 = [omega0[i] for i in range(len(omega0)) if mask[i] < 0.95]
    filt_omegaA = [omegaA[i] for i in range(len(omegaA)) if mask[i] < 0.95]

    # Check that filt_omega are not empty
    if not filt_omega or not filt_omega0 or not filt_omegaA or not filt_distance:
        print(f"Skipping folder: {folder_path}")
        continue

    # Plot the distance and omega on the same plot
    fig, ax1 = plt.subplots()
    col1 = "#1AFF1A"
    col2 = "#4B0092"
    width = len(distance) / 20
    # fig.set_size_inches(width, 7)
    ax1.set_xlabel("Position")
    ax1.set_ylabel("Distance", color=col1)
    ax1.bar(x=[a-0.25 for a in positions], height=distance, width=0.5, color=col1, alpha=1)
    ax1.tick_params(axis='y', labelcolor=col1)
    ax2 = ax1.twinx() # instantiate a second axis that shares the same x-axis
    ax2.set_ylabel("omega", color=col2)
    ax2.bar(x=[a+0.25 for a in positions], height=omega, width=0.5, color=col2, alpha=1)
    ax2.tick_params(axis='y', labelcolor=col2)
    plt.title("Distance and omega per site for gene " + name)
    plt.xlim([0, len(distance)])
    fig.tight_layout()
    os.makedirs(os.path.join(SOURCE_DIR, "results", folder_name), exist_ok=True)
    set_size(width, 7)
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/distance.pdf", bbox_inches='tight')

    # Plot the filtered distance and omega on the same plot
    fig, ax1 = plt.subplots()
    col1 = "#1AFF1A"
    col2 = "#4B0092"
    width = len(distance) / 20
    # fig.set_size_inches(width, 7)
    ax1.set_xlabel("Position")
    ax1.set_ylabel("Distance", color=col1)
    ax1.bar(x=[a-0.25 for a in filt_positions], height=filt_distance, width=0.5, color=col1, alpha=1)
    ax1.tick_params(axis='y', labelcolor=col1)
    ax2 = ax1.twinx() # instantiate a second axis that shares the same x-axis
    ax2.set_ylabel("omega", color=col2)
    ax2.bar(x=[a+0.25 for a in filt_positions], height=filt_omega, width=0.5, color=col2, alpha=1)
    ax2.tick_params(axis='y', labelcolor=col2)
    plt.title("Distance and omega per site, for variable sites, for gene " + name)
    plt.xlim([0, len(distance)])
    fig.tight_layout()
    os.makedirs(os.path.join(SOURCE_DIR, "results", folder_name), exist_ok=True)
    set_size(width, 7)
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/filtered-distance.pdf", bbox_inches='tight')

    # Plot filtered distance vs omega and save it
    yfit, lower_ci, upper_ci = linearfit(filt_omega, filt_distance, np.linspace(min(filt_omega), max(filt_omega), 100))
    plt.figure(figsize=(7, 7))
    plt.scatter(filt_omega, filt_distance, s=2)
    plt.fill_between(x=np.linspace(min(filt_omega), max(filt_omega), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
    plt.plot(np.linspace(min(filt_omega), max(filt_omega), 100), yfit, linewidth=2, color=lin_reg)
    plt.xlabel("Omega")
    plt.ylabel("Distance")
    plt.title("Distance as a function of omega, for variable positions, for gene " + name)
    plt.xscale("log")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/filtered-distance-vs-omega.pdf")

    # Plot filtered distance vs omega0 and save it
    yfit, lower_ci, upper_ci = linearfit(filt_omega0, filt_distance, np.linspace(min(filt_omega0), max(filt_omega0), 100))
    plt.figure(figsize=(7, 7))
    plt.scatter(filt_omega0, filt_distance, s=2)
    plt.fill_between(x=np.linspace(min(filt_omega0), max(filt_omega0), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
    plt.plot(np.linspace(min(filt_omega0), max(filt_omega0), 100), yfit, linewidth=2, color=lin_reg)
    plt.xlabel("Omega0")
    plt.ylabel("Distance")
    plt.title("Distance as a function of omega0, for variable positions, for gene " + name)
    plt.xscale("log")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/filtered-distance-vs-omega0.pdf")

    # Plot distance vs omegaA and save it
    yfit, lower_ci, upper_ci = linearfit(filt_omegaA, filt_distance, np.linspace(min(filt_omegaA), max(filt_omegaA), 100))
    plt.figure(figsize=(7, 7))
    plt.scatter(filt_omegaA, filt_distance, s=2)
    plt.fill_between(x=np.linspace(min(filt_omegaA), max(filt_omegaA), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
    plt.plot(np.linspace(min(filt_omegaA), max(filt_omegaA), 100), yfit, linewidth=2, color=lin_reg)
    plt.xlabel("OmegaA")
    plt.ylabel("Distance")
    plt.title("Distance as a function of omegaA (omega - omega0), for filtered positions, for gene " + name)
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/filtered-distance-vs-omegaA.pdf")


all_omegaA = [all_omega[i] - all_omega0[i] for i in range(len(all_omega))]

# Plot distance vs total omega and save it
yfit, lower_ci, upper_ci = linearfit(all_omega, all_distance, np.linspace(min(all_omega), max(all_omega), 100))
plt.figure(figsize=(7, 7))
plt.scatter(all_omega, all_distance, s=2)
plt.fill_between(x=np.linspace(min(all_omega), max(all_omega), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
plt.plot(np.linspace(min(all_omega), max(all_omega), 100), yfit, linewidth=2, color=lin_reg)
plt.xlabel("Omega")
plt.ylabel("Distance")
plt.title("Distance as a function of omega")
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/allgenes-distance-vs-omega.pdf")

# Plot distance vs total omega0 and save it
yfit, lower_ci, upper_ci = linearfit(all_omega0, all_distance, np.linspace(min(all_omega0), max(all_omega0), 100))
plt.figure(figsize=(7, 7))
plt.scatter(all_omega0, all_distance, s=2)
plt.fill_between(x=np.linspace(min(all_omega0), max(all_omega0), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
plt.plot(np.linspace(min(all_omega0), max(all_omega0), 100), yfit, linewidth=2, color=lin_reg)
plt.xlabel("Omega0")
plt.ylabel("Distance")
plt.title("Distance as a function of omega 0")
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/allgenes-distance-vs-omega0.pdf")

# Plot distance vs total omegaA and save it
yfit, lower_ci, upper_ci = linearfit(all_omegaA, all_distance, np.linspace(min(all_omegaA), max(all_omegaA), 100))
plt.figure(figsize=(7, 7))
plt.scatter(all_omegaA, all_distance, s=2)
plt.fill_between(x=np.linspace(min(all_omegaA), max(all_omegaA), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
plt.plot(np.linspace(min(all_omegaA), max(all_omegaA), 100), yfit, linewidth=2, color=lin_reg)
plt.xlabel("OmegaA")
plt.ylabel("Distance")
plt.title("Distance as a function of omegaA (omega - omega0)")
plt.yscale("log")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/allgenes-distance-vs-omegaA.pdf")

# Plot omegaA vs distance per gene and save it
yfit, lower_ci, upper_ci = linearfit(omega_per_gene, distance_per_gene, np.linspace(min(omega_per_gene), max(omega_per_gene), 100))
plt.figure(figsize=(7, 7))
for i, txt in enumerate(all_names):
    plt.annotate(txt, (omega_per_gene[i] + 0.0005, distance_per_gene[i] + 0.0005), fontsize=6)
plt.scatter(omega_per_gene, distance_per_gene, s=5)
plt.fill_between(x=np.linspace(min(omega_per_gene), max(omega_per_gene), 100), y1=lower_ci, y2=upper_ci, color=fill_reg, edgecolor=None)
plt.plot(np.linspace(min(omega_per_gene), max(omega_per_gene), 100), yfit, linewidth=2, color=lin_reg)
plt.xlabel("Average omega")
plt.ylabel("Average distance")
plt.rcParams["axes.titlesize"] = 10
plt.title("Average distance as a function of average omega, per gene")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/pergene-distance-vs-omega.pdf")
