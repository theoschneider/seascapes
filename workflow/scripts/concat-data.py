import os
import pandas as pd


# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# List folders and sim-folders
emp_folders = set(os.listdir(os.path.join(SOURCE_DIR, "processed")))
emp_folders.discard(".DS_Store")

sim_folders = set(os.listdir(os.path.join(SOURCE_DIR, "processed_sim")))
sim_folders.discard(".DS_Store")

# Take the intersection of the two lists
folders = emp_folders & sim_folders

# Initialize count
count = 0

# Initialize dataframe
df = pd.DataFrame(columns=["gene", "emp_distance", "emp_omega", "emp_omega0", "sim_distance", "sim_omega", "sim_omega0"])

# Iterate over the contents of folders and process them
for folder in folders:
    # Get folder with processed data
    emp_path = os.path.join(SOURCE_DIR, "processed", folder)
    sim_path = os.path.join(SOURCE_DIR, "processed_sim", folder)

    # Check if the necessary files exist in the folder
    if not (os.path.exists(os.path.join(emp_path, "distance.tsv"))
    and os.path.exists(os.path.join(emp_path, "Euarchontoglires.Whole.omega.ci0.025.tsv"))
    and os.path.exists(os.path.join(emp_path, "Euarchontoglires.Whole.omega_0.ci0.025.tsv"))
    and os.path.exists(os.path.join(sim_path, "distance.tsv"))
    and os.path.exists(os.path.join(sim_path, "Euarchontoglires.Whole.omega.ci0.025.tsv"))
    and os.path.exists(os.path.join(sim_path, "Euarchontoglires.Whole.omega_0.ci0.025.tsv"))):
        continue

    # Increment count
    count += 1

    # Get gene name
    gene = folder.split("_")[1]

    # Open distance.tsv
    emp_distance = pd.read_csv(os.path.join(emp_path, "distance.tsv"), sep="\t", header=0, index_col=None)
    emp_distance = emp_distance.iloc[:, 1].tolist()
    sim_distance = pd.read_csv(os.path.join(sim_path, "distance.tsv"), sep="\t", header=0, index_col=None)
    sim_distance = sim_distance.iloc[:, 1].tolist()

    # Open omega
    emp_omega = pd.read_csv(os.path.join(emp_path, "Euarchontoglires.Whole.omega.ci0.025.tsv"), sep="\t", header=0, index_col=0)
    emp_omega = emp_omega.iloc[1:, 1].tolist()
    sim_omega = pd.read_csv(os.path.join(sim_path, "Euarchontoglires.Whole.omega.ci0.025.tsv"), sep="\t", header=0, index_col=0)
    sim_omega = sim_omega.iloc[1:, 1].tolist()

    # Open omega0
    emp_omega0 = pd.read_csv(os.path.join(emp_path, "Euarchontoglires.Whole.omega_0.ci0.025.tsv"), sep="\t", header=0, index_col=0)
    emp_omega0 = emp_omega0.iloc[1:, 1].tolist()
    sim_omega0 = pd.read_csv(os.path.join(sim_path, "Euarchontoglires.Whole.omega_0.ci0.025.tsv"), sep="\t", header=0, index_col=0)
    sim_omega0 = sim_omega0.iloc[1:, 1].tolist()

    # Check that all lists are the same length
    assert len(emp_distance) == len(emp_omega) == len(emp_omega0) == len(sim_distance) == len(sim_omega) == len(sim_omega0), "Files do not have the same length"
    n_row = len(emp_distance)

    # Create a list of dictionaries, each representing a row
    rows = [{"gene": gene, "emp_distance": emp_distance[i], "emp_omega": emp_omega[i],
             "emp_omega0": emp_omega0[i], "sim_distance": sim_distance[i], "sim_omega": sim_omega[i],
             "sim_omega0": sim_omega0[i]} for i in range(n_row)]

    # Create a new DataFrame for these rows
    new_rows = pd.DataFrame(rows)

    # Use pd.concat to append the new rows to the existing DataFrame
    df = pd.concat([df, new_rows], ignore_index=True)


# Save the file
df.to_csv(os.path.join(SOURCE_DIR, "results", "concat_data.tsv"), sep="\t", header=True, index=False)

# Print the number of genes processed
print(f"Processed {count} genes, saved in {os.path.join(SOURCE_DIR, 'results', 'concat_data.tsv')} ({df.shape[0]} rows)")
