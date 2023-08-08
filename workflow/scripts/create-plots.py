import os
import pandas as pd
import matplotlib.pyplot as plt

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Define the range (FROM and TO) and the flag to process all folders (ALL)
FROM = 1
TO = 5
ALL = False

# Check if ALL is true, then adjust FROM and TO to process all folders
if ALL:
    FROM = 1
    TO = len(os.listdir(SOURCE_DIR))

# Iterate over the contents of SOURCE_DIR and process folders
for folder_name in os.listdir(os.path.join(SOURCE_DIR, "processed"))[FROM:TO + 1]:
    folder_path = os.path.join(SOURCE_DIR, "processed", folder_name)
    print(f"Processing folder: {folder_path}")
    # In every folder, open distance.tsv
    df = pd.read_csv(os.path.join(folder_path, "distance.tsv"), sep="\t", header=0, index_col=None)
    # Plot the distance
    plt.figure(figsize=(12, 7))
    plt.plot(df.iloc[:, 0], df.iloc[:, 1])
    plt.xlabel("Position")
    plt.ylabel("Distance")
    plt.title("Jensen-Shannon divergence between two subsets, compared to base distance for each subset")
    plt.tight_layout()
    # Create directory and save the plot
    os.makedirs(os.path.join(SOURCE_DIR, "results", folder_name), exist_ok=True)
    plt.savefig(SOURCE_DIR + "/results/" + folder_name + "/distance.pdf")

