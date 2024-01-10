import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

# Set the source directory
SOURCE_DIR = "/Users/theo/THÉO/seascapes"

# Read from concat_data.tsv
df = pd.read_csv(os.path.join(SOURCE_DIR, "processed", "emp_data.tsv"), sep="\t", header=0, index_col=None)

# Create empty df for plot
rows = []
groupby = df.groupby("gene")
for gene, group in groupby:
    # Bootstrap distance, omega and omegaA
    dist_res = scipy.stats.bootstrap((group["distance"].tolist(),), statistic=np.mean, confidence_level=0.95, n_resamples=1000, axis=0)
    omega_res = scipy.stats.bootstrap((group["omega"].tolist(),), statistic=np.mean, confidence_level=0.95, n_resamples=1000, axis=0)
    omegaA_res = scipy.stats.bootstrap((group["omegaA"].tolist(),), statistic=np.mean, confidence_level=0.95, n_resamples=1000, axis=0)
    # Save higher and lower
    dist_high = dist_res.confidence_interval[1]
    dist_low = dist_res.confidence_interval[0]
    omega_high = omega_res.confidence_interval[1]
    omega_low = omega_res.confidence_interval[0]
    omegaA_high = omegaA_res.confidence_interval[1]
    omegaA_low = omegaA_res.confidence_interval[0]
    # Calculate the mean distance
    mean_dist = np.mean(group["distance"])
    # Calculate the mean omega
    mean_omega = np.mean(group["omega"])
    # Calculate the mean omega0
    mean_omega0 = np.mean(group["omega0"])
    # Calculate the mean omegaA
    mean_omegaA = np.mean(group["omegaA"])
    # If omegaA bar crosses 0
    if omegaA_low < 0 and omegaA_high > 0:
        c = "#FF8E00"
    else:
        c = "C0"
    # Create a df with 1 row
    row = pd.DataFrame(data=[[gene, mean_dist, dist_low, dist_high, mean_omega, omega_low, omega_high, mean_omegaA, omegaA_low, omegaA_high, c]],
                       columns=["gene", "mean_dist", "low_dist", "high_dist", "mean_omega", "low_omega", "high_omega", "mean_omegaA", "low_omegaA", "high_omegaA", "color"])
    # Append to plot_df
    rows.append(row)
plot_df = pd.concat(rows, ignore_index=True)

# Plot the distances vs omega per gene in a scatterplot
fig, ax = plt.subplots(dpi=100, figsize=(10, 10))
plt.xscale("log")
plt.yscale("log")
# Scatter
ax.scatter(plot_df["mean_omega"], plot_df["mean_dist"], alpha=0.2, zorder=10)
# Error bars
ax.errorbar(plot_df["mean_omega"], plot_df["mean_dist"],
            xerr=[plot_df["mean_omega"] - plot_df["low_omega"], plot_df["high_omega"] - plot_df["mean_omega"]],
            yerr=[plot_df["mean_dist"] - plot_df["low_dist"], plot_df["high_dist"] - plot_df["mean_dist"]],
            fmt="none", ecolor="C0", elinewidth=1, capsize=0, alpha=0.2)
# Labels
plt.xlabel("Average $\omega$ (log-scaled)")
plt.ylabel("Average distance (log-scaled)")
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
ax.tick_params(axis='both', which='major', labelsize=16)
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/pergene-distance-vs-omega.pdf")
plt.show()
plt.close()

# Plot the distances vs omegaA per gene in a scatterplot
fig, ax = plt.subplots(dpi=100, figsize=(10, 10))
plt.yscale("log")
# Scatter
ax.scatter(plot_df["mean_omegaA"], plot_df["mean_dist"], alpha=0.2, zorder=10, color=plot_df["color"])
# Error bars
ax.errorbar(plot_df["mean_omegaA"], plot_df["mean_dist"],
            xerr=[plot_df["mean_omegaA"] - plot_df["low_omegaA"], plot_df["high_omegaA"] - plot_df["mean_omegaA"]],
            yerr=[plot_df["mean_dist"] - plot_df["low_dist"], plot_df["high_dist"] - plot_df["mean_dist"]],
            fmt="none", ecolor=plot_df["color"], elinewidth=1, capsize=0, alpha=0.2)
# Labels
plt.xlabel("Average $\omega_A$")
plt.ylabel("Average distance (log-scaled)")
ax.xaxis.label.set_size(16)
ax.yaxis.label.set_size(16)
ax.tick_params(axis='both', which='major', labelsize=16)
plt.tight_layout()
plt.savefig(SOURCE_DIR + "/results/pergene-distance-vs-omegaA.pdf")
plt.show()
plt.close()

