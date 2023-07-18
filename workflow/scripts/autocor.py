import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from libraries import js


def main(path: str, outdir: str):

    filenames = os.listdir(path)
    filenames.sort(key=lambda x: int(x.split("_")[1]))
    print(filenames)
    dataframes1 = [pd.read_csv(path+filename, sep="\t", header=0, index_col=0) for filename in filenames
                   if filename.split("_")[0] == "chain1"]
    dataframes2 = [pd.read_csv(path+filename, sep="\t", header=0, index_col=0) for filename in filenames
                   if filename.split("_")[0] == "chain2"]

    assert len(dataframes1) == len(dataframes2), "Number of windows in chain1 and chain2 are not equal"

    matrix = np.zeros((len(dataframes1), len(dataframes2)))

    for i in range(len(dataframes1)):
        df1 = dataframes1[i]

        for j in range(len(dataframes2)):
            df2 = dataframes2[j]

            distance = js(df1, df2)

            matrix[i, j] = np.mean(distance)

    # Save the matrix
    pd.DataFrame(matrix).to_csv(outdir + "all_distances.tsv", sep="\t", header=False, index=False)

    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, cmap='hot', interpolation='nearest') # origin='lower'
    plt.colorbar()
    plt.xticks(np.arange(0, len(dataframes1)), np.arange(1, len(dataframes1)+1))
    plt.yticks(np.arange(0, len(dataframes2)), np.arange(1, len(dataframes2)+1))
    plt.xlabel("Window from chain 2 (of size delta t)")
    plt.ylabel("Window from chain 1 (of size delta t)")
    plt.title("Autocorrelation plot, showing Jensen–Shannon divergence between all pairs of windows")
    plt.savefig(f"{outdir}/heatmap_bothchains_{int(2000 / len(dataframes1))}.pdf")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--path', required=True, type=str, dest="path",
                        help="The path to the directory containing all the .siteprofiles files")
    parser.add_argument('--outdir', required=True, type=str, dest="outdir",
                        help="Path to the output directory")

    args = parser.parse_args()
    main(args.path, args.outdir)
