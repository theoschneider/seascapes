import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Read from concat_data.tsv
df = pd.read_csv(os.path.join(SOURCE_DIR, "processed", "emp_data.tsv"), sep="\t", header=0, index_col=None)

# Create empty df for plot
rows = [pd.DataFrame(columns=["gene", "mean_dist", "mean_omega", "mean_omega0", "mean_omegaA"])]
groupby = df.groupby("gene")
for gene, group in groupby:
    # Calculate the mean distance
    mean_dist = np.mean(group["distance"])
    # Calculate the mean omega
    mean_omega = np.mean(group["omega"])
    # Calculate the mean omega0
    mean_omega0 = np.mean(group["omega0"])
    # Calculate the mean omegaA
    mean_omegaA = np.mean(group["omegaA"])
    # Create a df with 1 row
    row = pd.DataFrame([[gene, mean_dist, mean_omega, mean_omega0, mean_omegaA]], columns=["gene", "mean_dist", "mean_omega", "mean_omega0", "mean_omegaA"])
    # Append to plot_df
    rows.append(row)
plot_df = pd.concat(rows, ignore_index=True)

# Plot the distances vs omega per gene in a scatterplot
fig, ax = plt.subplots(dpi=100, figsize=(10, 10))
# Gene names
# for i, row in plot_df.iterrows():
#     if row["mean_omega"] > 0.6 or row["mean_dist"] > 0.1:
#         txt = row["gene"]
#         ax.annotate(txt, (row["mean_omega"] + 0.0005, row["mean_dist"] + 0.0005), fontsize=12)
# Scatter
ax.scatter(plot_df["mean_omega"], plot_df["mean_dist"], alpha=0.2)
# Labels
plt.xlabel("Average $\omega$ (log-scaled)")
plt.ylabel("Average distance (log-scaled)")
plt.xscale("log")
plt.yscale("log")
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
ax.tick_params(axis='both', which='major', labelsize=16)
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/pergene-distance-vs-omega.pdf")
plt.show()
plt.close()

# Plot the distances vs omegaA per gene in a scatterplot
fig, ax = plt.subplots(dpi=100, figsize=(10, 10))
# Gene names
# for i, row in plot_df.iterrows():
#     if row["mean_omegaA"] > 0.1 or row["mean_dist"] > 0.1:
#         txt = row["gene"]
#         ax.annotate(txt, (row["mean_omegaA"] + 0.0005, row["mean_dist"] + 0.0005), fontsize=12)
# Scatter
ax.scatter(plot_df["mean_omegaA"], plot_df["mean_dist"], alpha=0.2)
# Labels
plt.xlabel("Average $\omega_A$")
plt.ylabel("Average distance (log-scaled)")
plt.yscale("log")
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
ax.tick_params(axis='both', which='major', labelsize=16)
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/pergene-distance-vs-omegaA.pdf")
plt.show()
plt.close()

