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

# Calculate the error margin of distance per gene using bootstrap
rows = [pd.DataFrame(columns=["gene", "emp_dist", "emp_dist_low", "emp_dist_high", "sim_dist", "sim_dist_low", "sim_dist_high", "emp_omega"])]
groupby = df.groupby("gene")
for gene, group in groupby:
    # Get the distance values for this gene
    emp_gene_dist = group["emp_distance"].tolist()
    # Calculate the error margin
    emp_res = scipy.stats.bootstrap((emp_gene_dist,), statistic=np.mean, confidence_level=0.95,
                                n_resamples=1000, axis=0)
    # Same for simulated
    sim_gene_dist = group["sim_distance"].tolist()
    sim_res = scipy.stats.bootstrap((sim_gene_dist,), statistic=np.mean, confidence_level=0.95,
                                n_resamples=1000, axis=0)
    # Get omega
    emp_omega = np.mean(df[df["gene"] == gene]["emp_omega"])
    # Create a df with 1 row
    row = pd.DataFrame([[gene, np.mean(emp_gene_dist), emp_res.confidence_interval[0], emp_res.confidence_interval[1], np.mean(sim_gene_dist), sim_res.confidence_interval[0], sim_res.confidence_interval[1], emp_omega]],
                       columns=["gene", "emp_dist", "emp_dist_low", "emp_dist_high",
                                "sim_dist", "sim_dist_low", "sim_dist_high", "emp_omega"])
    rows.append(row)

# Create plot_df
plot_df = pd.concat(rows, ignore_index=True)

# Plot the distances per gene in a scatterplot, with confidence interval, with omega as color
fig, ax = plt.subplots(dpi=130, figsize=(12, 10))
# x=y line
x = np.linspace(min(plot_df["emp_dist"]), max(plot_df["emp_dist"]), 100)
ax.plot(x, x, color="red")
# error bars, colored for omega
ax.errorbar(plot_df["emp_dist"], plot_df["sim_dist"],
            xerr=[plot_df["emp_dist"] - plot_df["emp_dist_low"], plot_df["emp_dist_high"] - plot_df["emp_dist"]],
            yerr=[plot_df["sim_dist"] - plot_df["sim_dist_low"], plot_df["sim_dist_high"] - plot_df["sim_dist"]],
            fmt="none", ecolor="grey", elinewidth=1, capsize=2)
# scatter
scatter = ax.scatter(plot_df["emp_dist"], plot_df["sim_dist"], c=plot_df["emp_omega"], cmap="viridis", zorder=10)
# colorbar, with label omega
cbar = fig.colorbar(scatter, ax=ax)
cbar.ax.set_title("$\omega$", fontsize=16)
cbar.ax.tick_params(labelsize=16)
# add gene names if one of the error bars does not overlap with x=y line and count
count = 0
for i, row in plot_df.iterrows():
    if row["emp_dist_low"] > row["sim_dist"] or row["sim_dist_high"] < row["emp_dist"]:
        count += 1
        # lower right corner
        ax.annotate(row["gene"],  # this is the text
                    (row["emp_dist"], row["sim_dist"]),  # this is the point to label
                    textcoords="offset points",
                    xytext=(10, -40),  # distance from text to points (x,y)
                    ha='left',  # text horizontal alignment
                    arrowprops=dict(arrowstyle="-"),
                    fontsize=12)
    elif row["emp_dist_high"] < row["sim_dist"] or row["sim_dist_low"] > row["emp_dist"]:
        count += 1
        # upper left corner
        ax.annotate(row["gene"],  # this is the text
                    (row["emp_dist"], row["sim_dist"]),  # this is the point to label
                    textcoords="offset points",
                    xytext=(-100, 10),  # distance from text to points (x,y)
                    ha='left',  # text horizontal alignment
                    arrowprops=dict(arrowstyle="-"),
                    fontsize=12)
ax.set_xlabel("Empirical distance", fontsize=16)
ax.set_ylabel("Simulated distance", fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=16)
plt.tight_layout()
plt.savefig(os.path.join(SOURCE_DIR, "results", "pergene_emp-vs-sim_ci.pdf"), bbox_inches='tight')
plt.show()
plt.close()
print("Total number of genes:", len(plot_df))
print("Number of genes with non-overlapping error bars:", count)

# Calculate r2 and slope
r2 = np.corrcoef(plot_df["emp_dist"], plot_df["sim_dist"])[0, 1] ** 2
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(plot_df["emp_dist"], plot_df["sim_dist"])
print("r2:", r2)
print("slope:", slope)
print("p_value:", p_value)
