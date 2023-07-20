import argparse
import pandas as pd
from libraries import js


def main(path1, path2, outpath):
    df1 = pd.read_csv(path1, sep="\t", header=0, index_col=0)
    df2 = pd.read_csv(path2, sep="\t", header=0, index_col=0)
    df_distance = js(df1, df2)
    df_distance.to_csv(outpath, sep="\t")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--path1', required=True, type=str, dest="path1",
                        help="First profile path")
    parser.add_argument('--path2', required=True, type=str, dest="path2",
                        help="Second profile path")
    parser.add_argument('--outpath', required=True, type=str, dest="outpath",
                        help="Output file path")

    args = parser.parse_args()
    main(args.path1, args.path2, args.outpath)
