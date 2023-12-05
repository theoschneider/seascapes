import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Set the clade
CLADE = "Euarchontoglires"

# Get a list of folders that are both in processed and processed_sim
processed_folders = os.listdir(os.path.join(SOURCE_DIR, "processed"))
if ".DS_Store" in processed_folders:
    processed_folders.remove(".DS_Store")
processed_sim_folders = os.listdir(os.path.join(SOURCE_DIR, "processed_sim"))
if ".DS_Store" in processed_sim_folders:
    processed_sim_folders.remove(".DS_Store")
folders = [folder for folder in processed_folders if folder in processed_sim_folders]

print("Processing " + str(len(folders)) + " folders")

# Iterate for every folder
for folder_name in folders:

    # Check if both distance files exist and also the mask
    if not (os.path.exists(os.path.join(SOURCE_DIR, "processed", folder_name, "distance.tsv"))
    and os.path.exists(os.path.join(SOURCE_DIR, "processed_sim", folder_name, "distance.tsv"))
    and os.path.exists(os.path.join(SOURCE_DIR, "processed", folder_name, CLADE + ".Whole.mask.tsv"))):
        print(f"Skipping folder: {folder_name}")
        continue

    print(f"Processing folder: {folder_name}")

    empirical_path = os.path.join(SOURCE_DIR, "processed", folder_name)
    simulated_path = os.path.join(SOURCE_DIR, "processed_sim", folder_name)

    # Check if the result folder exists, if not create it
    if not os.path.exists(os.path.join(SOURCE_DIR, "results", folder_name)):
        os.makedirs(os.path.join(SOURCE_DIR, "results", folder_name))

    # Open the mask
    mask = pd.read_csv(os.path.join(empirical_path, CLADE + ".Whole.mask.tsv"), sep="\t", header=0, index_col=None)
    mask = mask.iloc[:, 1]

    # Open empirical distance
    distance_df = pd.read_csv(os.path.join(empirical_path, "distance.tsv"), sep="\t", header=0, index_col=None)

    # Filter distance, keep 1 col
    filt_dist = distance_df.iloc[:, 1][mask < 0.95]

    # Open simulated distance
    sim_distance_df = pd.read_csv(os.path.join(simulated_path, "distance.tsv"), sep="\t", header=0, index_col=None)

    # Filter distance, keep 1 col
    sim_filt_dist = sim_distance_df.iloc[:, 1][mask < 0.95]

    # Plot the 2 distances in a boxplot
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.boxplot([filt_dist, sim_filt_dist], labels=["Empirical", "Simulated"])
    ax.set_ylabel("Distance")
    ax.set_title(folder_name)
    plt.savefig(os.path.join(SOURCE_DIR, "results", folder_name, "emp-vs-sim.pdf"))
    plt.close()