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
df = pd.DataFrame(columns=["gene", "emp_distance", "emp_omega", "emp_omega0", "emp_omegaA", "sim_distance", "sim_omega", "sim_omega0", "sim_omegaA"])

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

    # Calculate omegaA
    emp_omegaA = [emp_omega[i] - emp_omega0[i] for i in range(len(emp_omega))]
    sim_omegaA = [sim_omega[i] - sim_omega0[i] for i in range(len(sim_omega))]

    # Check that all lists are the same length
    assert len(emp_distance) == len(emp_omega) == len(emp_omega0) == len(sim_distance) == len(sim_omega) == len(sim_omega0) == len(emp_omegaA) == len(sim_omegaA), "Files do not have the same length"
    n_row = len(emp_distance)

    # Create a list of dictionaries, each representing a row
    rows = [{"gene": gene, "emp_distance": emp_distance[i], "emp_omega": emp_omega[i], "emp_omega0": emp_omega0[i],
             "emp_omegaA": emp_omegaA[i], "sim_distance": sim_distance[i], "sim_omega": sim_omega[i],
             "sim_omega0": sim_omega0[i], "sim_omegaA": sim_omegaA[i]} for i in range(n_row)]

    # Create a new DataFrame for these rows
    new_rows = pd.DataFrame(rows)

    # Use pd.concat to append the new rows to the existing DataFrame
    df = pd.concat([df, new_rows], ignore_index=True)


# Same for empirical data only
# Initialize count
emp_count = 0

# Initialize dataframe
emp_df = pd.DataFrame(columns=["gene", "distance", "omega", "omega0", "omegaA", "mask"])

# Iterate over folders
for folder in emp_folders:
    # Get folder with processed data
    emp_path = os.path.join(SOURCE_DIR, "processed", folder)

    # Check if the necessary files exist in the folder
    if not (os.path.exists(os.path.join(emp_path, "distance.tsv"))
    and os.path.exists(os.path.join(emp_path, "Euarchontoglires.Whole.omega.ci0.025.tsv"))
    and os.path.exists(os.path.join(emp_path, "Euarchontoglires.Whole.omega_0.ci0.025.tsv"))
    and os.path.exists(os.path.join(emp_path, "Euarchontoglires.Whole.mask.tsv"))):
        continue

    # Increment count
    emp_count += 1

    # Get gene name
    gene = folder.split("_")[1]

    # Open distance.tsv
    distance = pd.read_csv(os.path.join(emp_path, "distance.tsv"), sep="\t", header=0, index_col=None)
    distance = distance.iloc[:, 1].tolist()

    # Open omega
    omega = pd.read_csv(os.path.join(emp_path, "Euarchontoglires.Whole.omega.ci0.025.tsv"), sep="\t", header=0, index_col=0)
    omega = omega.iloc[1:, 1].tolist()

    # Open omega0
    omega0 = pd.read_csv(os.path.join(emp_path, "Euarchontoglires.Whole.omega_0.ci0.025.tsv"), sep="\t", header=0, index_col=0)
    omega0 = omega0.iloc[1:, 1].tolist()

    # Calculate omegaA
    omegaA = [omega[i] - omega0[i] for i in range(len(omega))]

    # Open mask
    mask = pd.read_csv(os.path.join(emp_path, "Euarchontoglires.Whole.mask.tsv"), sep="\t", header=0, index_col=None)
    mask = mask.iloc[:, 1].tolist()

    # Check that all lists are the same length
    assert len(distance) == len(omega) == len(omega0) == len(omegaA) == len(mask), "Files do not have the same length"
    n_row = len(distance)

    # Create a list of dictionaries, each representing a row
    rows = [{"gene": gene, "distance": distance[i], "omega": omega[i], "omega0": omega0[i],
             "omegaA": omegaA[i], "mask": mask[i]} for i in range(n_row)]

    # Create a new DataFrame for these rows
    new_rows = pd.DataFrame(rows)

    # Use pd.concat to append the new rows to the existing DataFrame
    emp_df = pd.concat([emp_df, new_rows], ignore_index=True)


# Save the files
df.to_csv(os.path.join(SOURCE_DIR, "processed", "concat_data.tsv"), sep="\t", header=True, index=False)
emp_df.to_csv(os.path.join(SOURCE_DIR, "processed", "emp_data.tsv"), sep="\t", header=True, index=False)

# Print the number of genes processed
print(f"Emp+Sim: Processed {count} genes, saved in {os.path.join(SOURCE_DIR, 'processed', 'concat_data.tsv')} ({df.shape[0]} rows)")
print(f"Emp only: Processed {emp_count} genes, saved in {os.path.join(SOURCE_DIR, 'processed', 'emp_data.tsv')} ({emp_df.shape[0]} rows)")
