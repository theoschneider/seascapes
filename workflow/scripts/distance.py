import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from libraries import js


def main(file1: str, file2: str, file3: str, file4: str, output_dir: str):

    df1 = pd.read_csv(file1, sep="\t", header=0, index_col=0)
    df2 = pd.read_csv(file2, sep="\t", header=0, index_col=0)
    df3 = pd.read_csv(file3, sep="\t", header=0, index_col=0)
    df4 = pd.read_csv(file4, sep="\t", header=0, index_col=0)

    # Calculate KL for base subset1; base subset2; subset1 vs subset2
    distance1 = js(df1, df2)
    distance2 = js(df3, df4)

    avg_sub_1 = (df1 + df2) / 2
    avg_sub_2 = (df3 + df4) / 2

    distance3 = js(avg_sub_1, avg_sub_2)

    distance1.to_csv(output_dir+"/base_sub1.tsv", sep="\t", header=False)
    distance2.to_csv(output_dir+"/base_sub2.tsv", sep="\t", header=False)
    distance3.to_csv(output_dir+"/distance.tsv", sep="\t", header=False)

    plt.figure(figsize=(12, 7))
    plt.plot(distance1)
    plt.plot(distance2)
    plt.plot(distance3)
    plt.xlabel("Position")
    plt.ylabel("Distance")
    plt.legend(["Base distance for subset 1", "Base distance for subset 2", "Distance between subsets"])
    plt.title("Jensen-Shannon divergence between two subsets, compared to base distance for each subset")
    plt.tight_layout()
    plt.savefig(output_dir+"/distance.pdf")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--file1', required=True, type=str, dest="file1",
                        help="First run of subset 1, path to .siteprofiles file")
    parser.add_argument('--file2', required=True, type=str, dest="file2",
                        help="Second run of subset 1, path to .siteprofiles file")
    parser.add_argument('--file3', required=True, type=str, dest="file3",
                        help="First run of subset 2, path to .siteprofiles file")
    parser.add_argument('--file4', required=True, type=str, dest="file4",
                        help="Second run of subset 2, path to .siteprofiles file")
    parser.add_argument('--outdir', required=True, type=str, dest="output_dir",
                        help="Output file path")

    args = parser.parse_args()
    main(args.file1, args.file2, args.file3, args.file4, args.output_dir)
