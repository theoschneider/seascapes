import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import os


def main(path: str, outdir: str):

    filenames = os.listdir(path)
    filenames.sort(key=lambda x: int(x.split("_")[0]))
    print(filenames)
    dataframes = [pd.read_csv(path+filename, sep="\t", header=0) for filename in filenames]
    matrix = np.zeros((len(filenames), len(filenames)))

    for i in range(len(filenames)):
        df1 = dataframes[i]

        for j in range(i+1, len(filenames)):
            df2 = dataframes[j]

            kl1 = np.sum(df1 * np.log(df1 / df2), axis=1)
            kl2 = np.sum(df2 * np.log(df2 / df1), axis=1)

            distance = (kl1 + kl2)/2

            matrix[i, j] = np.mean(distance)

    # Save the matrix
    pd.DataFrame(matrix).to_csv(outdir + "all_distances.tsv", sep="\t", header=False, index=False)

    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, cmap='hot', interpolation='nearest') # origin='lower'
    plt.colorbar()
    plt.xticks(np.arange(0, len(filenames)), np.arange(1, len(filenames)+1))
    plt.yticks(np.arange(0, len(filenames)), np.arange(1, len(filenames)+1))
    plt.xlabel("Second window (of size delta t)")
    plt.ylabel("First window (of size delta t)")
    plt.title("Autocorrelation plot, showing Jensen–Shannon divergence between all pairs of windows")
    plt.savefig(outdir + "heatmap.pdf")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--path', required=True, type=str, dest="path",
                        help="The path to the directory containing all the .siteprofiles files")
    parser.add_argument('--outdir', required=True, type=str, dest="outdir",
                        help="Path to the output directory")

    args = parser.parse_args()
    main(args.path, args.outdir)
