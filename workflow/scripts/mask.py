from libraries import open_ali
import pandas as pd
import argparse
from collections import defaultdict


def main(ali_path: str, out_path: str) -> None:
    ali_dic = open_ali(ali_path)

    # write in a file, for every position, the frequency of the minor amino-acid
    with open(out_path, "w") as f:
        f.write("site\tminor_aa_freq\n")
        l = len(list(ali_dic.keys())[0])
        for pos in range(l/3):
            freq = defaultdict(0)

            for sp in ali_dic:
                freq[ali_dic[sp][pos:pos+3]] += 1

            minor_aa = min(freq, key=freq.get)
            minor_aa_freq = freq[minor_aa] / l

            f.write(f"{pos}\t{minor_aa_freq}\n")

    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--alipath', required=True, type=str, dest="alipath",
                        help="Path to the input alignment file")
    parser.add_argument('--outpath', required=True, type=str, dest="outpath",
                        help="Path to the output text file")

    args = parser.parse_args()

    main(args.alipath, args.outpath)
