import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import ttest_rel, ks_2samp
from libraries import cohen_d

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Read from concat_data.tsv
df = pd.read_csv(os.path.join(SOURCE_DIR, "processed", "concat_data.tsv"), sep="\t", header=0, index_col=None)

# Number of obs
n_obs = len(df)

# Perform t-test
t, p = ttest_rel(df["emp_distance"], df["sim_distance"])
print("t:", t)
print("p:", p)
# Perform KS test
D, p = ks_2samp(df["emp_distance"], df["sim_distance"])
print("D:", D)
print("p:", p)
# Get a r2
x = df["emp_distance"]
y = df["sim_distance"]
# Calculate Cohen's D
d = cohen_d(x, y)
print("d:", d)

# Create 2 lists
distance_col = df["emp_distance"].tolist() + df["sim_distance"].tolist()
type_col = ["Empirical"]*n_obs + ["Simulated"]*n_obs

# Create plot df
plot_df = pd.DataFrame({"distance": distance_col, "type": type_col})

# Create violin plot
fig, ax = plt.subplots(dpi=130, figsize=(8, 8))
sns.set_style("whitegrid")
p = sns.violinplot(x="type", y="distance", data=plot_df, inner="quart", log_scale=True, ax=ax, color="#517db3")
p.set_xlabel("Type of data", fontsize=16)
p.set_ylabel("Distance (log-scaled)", fontsize=16)
p.tick_params(axis='both', which='major', labelsize=16)
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/violin-plot.pdf")
plt.show()
plt.close()
