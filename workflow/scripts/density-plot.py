import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Read from emp_data.tsv
df = pd.read_csv(os.path.join(SOURCE_DIR, "processed", "emp_data.tsv"), sep="\t", header=0, index_col=None)

# Keep filtered positions
# df = df[df["mask"] < 0.95]

# Clip omega and distance
df["omega"] = df["omega"].clip(3e-3, 2)
df["distance"] = df["distance"].clip(3e-3, 1)

# Create a 2D density plot with log-transformed distance
p = sns.jointplot(x='omega', y='distance', data=df, kind='kde', cmap='Blues', fill=True, height=6, log_scale=(True, True))

# Add labels and title
plt.xlabel('$\omega$ (clipped and log-scaled)', fontsize=14)
plt.ylabel('Distance (log-scaled)', fontsize=14)

# Increase labels and ticks size
plt.tick_params(axis='both', which='major', labelsize=14)
# plt.rcParams["axes.titlesize"] = 14
# plt.rcParams["axes.labelsize"] = 14
# plt.rcParams["xtick.labelsize"] = 14
# plt.rcParams["ytick.labelsize"] = 14

p.fig.tight_layout()
plt.savefig(SOURCE_DIR + "/results/density-plot-omega.pdf")
plt.show()
plt.close()


# Clip omega and distance
df["omegaA"] = df["omegaA"].clip(-0.5, 0.5)

# Create a 2D density plot with log-transformed distance
p = sns.jointplot(x='omegaA', y='distance', data=df, kind='kde', cmap='Blues', fill=True, height=6, log_scale=(False, True))

# Add labels and title
plt.xlabel('$\omega_A$ (clipped and log-scaled)', fontsize=14)
plt.ylabel('Distance (log-scaled)', fontsize=14)

# Increase labels and ticks size
plt.tick_params(axis='both', which='major', labelsize=14)
# plt.rcParams["axes.titlesize"] = 14
# plt.rcParams["axes.labelsize"] = 14
# plt.rcParams["xtick.labelsize"] = 14
# plt.rcParams["ytick.labelsize"] = 14

p.fig.tight_layout()
plt.savefig(SOURCE_DIR + "/results/density-plot-omegaA.pdf")
plt.show()
plt.close()
