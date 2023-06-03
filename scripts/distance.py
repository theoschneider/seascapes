import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main(file1: str, file2: str, output_dir: str):
    df1 = pd.read_csv(file1, sep="\t", header=0)
    df2 = pd.read_csv(file2, sep="\t", header=0)

    kl1 = np.sum(df1 * np.log(df1 / df2), axis=1)
    kl2 = np.sum(df2 * np.log(df2 / df1), axis=1)

    distance = (kl1 + kl2)/2

    distance.to_csv(output_dir+"/distance.tsv", sep="\t", header=False)

    plt.figure(figsize=(12, 7))
    plt.plot(distance)
    plt.xlabel("Position")
    plt.ylabel("Distance")
    plt.title("Distance between two site profiles")
    plt.tight_layout()
    plt.savefig(output_dir+"/distance.pdf")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--file1', required=True, type=str, dest="file1", help="First .siteprofiles file path")
    parser.add_argument('--file2', required=True, type=str, dest="file2", help="Second .siteprofiles file path")
    parser.add_argument('--output', required=True, type=str, dest="output_dir", help="Output file path")

    args = parser.parse_args()
    main(args.file1, args.file2, args.output_dir)
