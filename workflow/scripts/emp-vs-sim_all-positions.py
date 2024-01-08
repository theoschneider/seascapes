import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Set the clade
CLADE = "Euarchontoglires"

# Read from concat_data.tsv
df = pd.read_csv(os.path.join(SOURCE_DIR, "processed", "concat_data.tsv"), sep="\t", header=0, index_col=None)

# Plot the distances in a scatterplot with omega as color
fig, ax = plt.subplots(figsize=(12, 10))
# x=y line
x = np.linspace(min(df["emp_distance"]), max(df["emp_distance"]), 100)
ax.plot(x, x, color="red")
# scatter
scatter = ax.scatter(df["emp_distance"], df["sim_distance"], c=df["emp_omega"], cmap="viridis", zorder=10)
# colorbar, with label omega
fig.colorbar(scatter, ax=ax).ax.set_title("omega")
# Labels
plt.xlabel("Empirical distance")
plt.ylabel("Simulated distance")
plt.rcParams["axes.titlesize"] = 10
plt.title("Distance between two subsets, for all positions, for all genes (n = 256)")
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/allpositions_emp-vs-sim_scatter.pdf")
plt.close()

