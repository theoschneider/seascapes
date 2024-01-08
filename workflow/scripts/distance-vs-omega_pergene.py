import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Read from concat_data.tsv
df = pd.read_csv(os.path.join(SOURCE_DIR, "processed", "emp_data.tsv"), sep="\t", header=0, index_col=None)

# Create empty df for plot
plot_df = pd.DataFrame(columns=["gene", "mean_dist", "mean_omega", "mean_omega0", "mean_omegaA"])

for gene in df["gene"].unique():
    # Calculate the mean distance
    mean_dist = np.mean(df[df["gene"] == gene]["distance"])
    # Calculate the mean omega
    mean_omega = np.mean(df[df["gene"] == gene]["omega"])
    # Calculate the mean omega0
    mean_omega0 = np.mean(df[df["gene"] == gene]["omega0"])
    # Calculate the mean omegaA
    mean_omegaA = np.mean(df[df["gene"] == gene]["omegaA"])
    # Create a df with 1 row
    row = pd.DataFrame([[gene, mean_dist, mean_omega, mean_omega0, mean_omegaA]], columns=["gene", "mean_dist", "mean_omega", "mean_omega0", "mean_omegaA"])
    # Append to plot_df
    plot_df = pd.concat([plot_df, row], ignore_index=True)

# Plot the distances vs omega per gene in a scatterplot
fig, ax = plt.subplots(figsize=(10, 10))
# Gene names
for i, row in plot_df.iterrows():
    if row["mean_omega"] > 0.6 or row["mean_dist"] > 0.1:
        txt = row["gene"]
        ax.annotate(txt, (row["mean_omega"] + 0.0005, row["mean_dist"] + 0.0005), fontsize=6)
# Scatter
ax.scatter(plot_df["mean_omega"], plot_df["mean_dist"])
# Labels
plt.xlabel("Average omega")
plt.ylabel("Average distance")
plt.rcParams["axes.titlesize"] = 10
plt.title("Average distance as a function of average omega, per gene")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/pergene-distance-vs-omega.pdf")
plt.close()

# Plot the distances vs omega0 per gene in a scatterplot
fig, ax = plt.subplots(figsize=(10, 10))
# Gene names
for i, row in plot_df.iterrows():
    if row["mean_omega0"] > 0.6 or row["mean_dist"] > 0.1:
        txt = row["gene"]
        ax.annotate(txt, (row["mean_omega0"] + 0.0005, row["mean_dist"] + 0.0005), fontsize=6)
# Scatter
ax.scatter(plot_df["mean_omega0"], plot_df["mean_dist"])
# Labels
plt.xlabel("Average omega0")
plt.ylabel("Average distance")
plt.rcParams["axes.titlesize"] = 10
plt.title("Average distance as a function of average omega0, per gene")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/pergene-distance-vs-omega0.pdf")
plt.close()

# Plot the distances vs omegaA per gene in a scatterplot
fig, ax = plt.subplots(figsize=(10, 10))
# Gene names
for i, row in plot_df.iterrows():
    if row["mean_omegaA"] > 0.6 or row["mean_dist"] > 0.1:
        txt = row["gene"]
        ax.annotate(txt, (row["mean_omegaA"] + 0.0005, row["mean_dist"] + 0.0005), fontsize=6)
# Scatter
ax.scatter(plot_df["mean_omegaA"], plot_df["mean_dist"])
# Labels
plt.xlabel("Average omegaA")
plt.ylabel("Average distance")
plt.rcParams["axes.titlesize"] = 10
plt.title("Average distance as a function of average omegaA, per gene")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/pergene-distance-vs-omegaA.pdf")
plt.close()


