import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Read from emp_data.tsv
df = pd.read_csv(os.path.join(SOURCE_DIR, "processed", "emp_data.tsv"), sep="\t", header=0, index_col=None)

# Keep filtered positions
df = df[df["mask"] < 0.95]

# Add 1 to omegaA col to avoid negative values
df["omegaA"] += 1

# Create a 2D density plot with log-transformed data
p = sns.jointplot(x='omegaA', y='distance', data=df, kind='kde', cmap='Blues', fill=True, log_scale=True, height=9)
p.fig.tight_layout()
p.fig.subplots_adjust(top=0.95)  # Reduce plot to make room
p.fig.suptitle("Distance as a function of omegaA", fontsize=16)

# Add labels and title
plt.xlabel('OmegaA (log-scaled with a constant shift)')
plt.ylabel('Distance (log-scaled)')

plt.savefig(SOURCE_DIR + "/results/density-plot-omegaA.pdf")
plt.show()
plt.close()

