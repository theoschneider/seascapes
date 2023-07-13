import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

path = "/Users/theo/THÉO/seascapes/"

# Load data

distance = pd.read_csv(path + "results/rodents_subset/distance.tsv", sep="\t", header=None, index_col=0)
distance = distance.iloc[:, 0].tolist()

omega0 = pd.read_csv(path + "processed/mutsel_TSPAN6_1.ci0.025.tsv", sep="\t", header=0, index_col=0)
omega0 = omega0.iloc[1:, 1].tolist()

omega = pd.read_csv(path + "processed/classical_mutsel_TSPAN6.ci0.025.tsv", sep="\t", header=0, index_col=0)
omega = omega.iloc[1:, 1].tolist()

omegaA = [omega[i] - omega0[i] for i in range(len(omega))]


# Plots

plt.figure(figsize=(7, 7))
plt.scatter(omega0, distance, s=2)
plt.xlabel("Omega 0")
plt.ylabel("Distance")
plt.title("Distance vs omega 0")
plt.tight_layout()
plt.savefig(path + "results/rodents_subset/distance-vs-omega0.pdf")

plt.figure(figsize=(7, 7))
plt.scatter(omega, distance, s=2)
plt.xlabel("Omega")
plt.ylabel("Distance")
plt.title("Distance vs omega")
plt.tight_layout()
plt.savefig(path + "results/rodents_subset/distance-vs-omega.pdf")

plt.figure(figsize=(7, 7))
plt.scatter(omegaA, distance, s=2)
plt.xlabel("OmegaA")
plt.ylabel("Distance")
plt.title("Distance vs omegaA (omega - omega0)")
plt.tight_layout()
plt.savefig(path + "results/rodents_subset/distance-vs-omegaA.pdf")