import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Set the clade
CLADE = "Euarchontoglires"

# Read from concat_data.tsv
df = pd.read_csv(os.path.join(SOURCE_DIR, "processed", "concat_data.tsv"), sep="\t", header=0, index_col=None)

# Plot the distances in a scatterplot with omega as color
fig, ax = plt.subplots(dpi=130, figsize=(10, 10))
# x=y line
x = np.linspace(min(df["emp_distance"]), max(df["emp_distance"]), 100)
ax.plot(x, x, color="red")
# scatter
scatter = ax.scatter(df["emp_distance"], df["sim_distance"], zorder=10, alpha=0.05)
# Labels
plt.xlabel("Empirical distance", fontsize=16)
plt.ylabel("Simulated distance", fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=16)
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/allpositions_emp-vs-sim_scatter.pdf")
plt.show()
plt.close()

# Calculate r2 and slope
x = df["emp_distance"]
y = df["sim_distance"]
r2 = np.corrcoef(x, y)[0, 1]**2
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
print("r2:", r2)
print("slope:", slope)
print("p_value:", p_value)
