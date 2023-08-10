import os
import pandas as pd
import matplotlib.pyplot as plt

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

# Check if ALL is true, then adjust FROM and TO to process all folders
if ALL:
    FROM = 1
    TO = len(os.listdir(SOURCE_DIR))

# Check if some genes have been specified
if GENES:
    genes_list = GENES
else:
    genes_list = os.listdir(os.path.join(SOURCE_DIR, "processed"))[FROM:TO + 1]

# Iterate over the contents of SOURCE_DIR and process folders
for folder_name in genes_list:
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

