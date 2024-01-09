import os
import pandas as pd
import matplotlib.pyplot as plt
from libraries import set_size
import numpy as np


# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"
df = pd.read_csv(os.path.join(SOURCE_DIR, "processed", "emp_data.tsv"), sep="\t", header=0, index_col=None)
gene = "NUDT18"

# Keep gene of interest
df = df[df["gene"] == gene]
# Where mask is < 0.95, set distance and omega to 0
df.loc[df["mask"] < 0.95, ["distance", "omega"]] = 0
positions = [i for i in range(1, len(df)+1)]

# Plot the distance and omega on the same plot
fig, ax1 = plt.subplots(dpi=130, figsize=(16, 7))
col1 = "#003F70"
col2 = "#FF8E00"
width = len(df) / 20
ax1.set_xlabel("Position")
ax1.set_ylabel("Distance", color=col1)
ax1.bar(x=[a - 0.25 for a in positions], height=df["distance"], width=0.5, color=col1, alpha=1)
ax1.tick_params(axis='y', labelcolor=col1)
ax2 = ax1.twinx()  # instantiate a second axis that shares the same x-axis
ax2.set_ylabel("$\omega$", color=col2)
ax2.bar(x=[a + 0.25 for a in positions], height=df["omega"], width=0.5, color=col2, alpha=1)
ax2.tick_params(axis='y', labelcolor=col2)
plt.xlim([0, len(df)])
# Increase ticks and size
ax1.tick_params(axis='both', which='major', labelsize=14)
ax2.tick_params(axis='both', which='major', labelsize=14)
# Increase labels size
ax1.xaxis.label.set_size(14)
ax1.yaxis.label.set_size(14)
ax2.yaxis.label.set_size(14)
# set_size(width, 7)
fig.tight_layout()
plt.savefig(SOURCE_DIR + "/results/" + gene + "_distance-omega.pdf")
plt.show()

