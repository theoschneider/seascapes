import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def main(path: str, outdir: str):

    filenames = os.listdir(path)
    matrix = np.zeros((len(filenames), len(filenames)))

    for i in range(len(filenames)):
        for j in range(i+1, len(filenames)):

            df1 = pd.read_csv(path+filenames[i], sep="\t", header=0)
            df2 = pd.read_csv(path+filenames[j], sep="\t", header=0)

            kl1 = np.sum(df1 * np.log(df1 / df2), axis=1)
            kl2 = np.sum(df2 * np.log(df2 / df1), axis=1)

            distance = (kl1 + kl2)/2

            matrix[i, j] = np.mean(distance)

    # Save the matrix
    pd.DataFrame(matrix).to_csv(outdir + "all_distances.tsv", sep="\t", header=False, index=False)

    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.xlabel("Second window (delta t)")
    plt.ylabel("First window (delta t)")
    plt.title("Autocorrelation between all pairs of site profiles")
    plt.savefig(outdir + "heatmap.pdf")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--path', required=True, type=str, dest="path",
                        help="The path to the directory containing all the .siteprofiles files")
    parser.add_argument('--outdir', required=True, type=str, dest="outdir",
                        help="Path to the output directory")

    args = parser.parse_args()
    main(args.path, args.outdir)
