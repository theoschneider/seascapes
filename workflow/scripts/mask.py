from libraries import open_ali
import argparse
from collections import defaultdict


# Define the codon table
codontable = defaultdict(lambda: "-")
codontable.update({
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'})


def main(ali_path: str, out_path: str) -> None:
    ali_dic = open_ali(ali_path)

    # write in a file, for every position, the frequency of the minor amino-acid
    with open(out_path, "w") as f:
        f.write("site\tmajor_aa_freq\n")
        l = len(list(ali_dic.values())[0].rstrip()) / 3

        assert l.is_integer(), "Alignment length is not a multiple of 3"

        l = int(l)
        for pos in range(l):
            n_species = 0
            freq = defaultdict(lambda: 0)

            for sp in ali_dic.keys():
                aa = codontable[ali_dic[sp][pos*3:pos*3+3]]

                if aa != "-":
                    freq[aa] += 1
                    n_species += 1

            major_freq = max(freq.values()) / n_species

            f.write(f"{pos+1}\t{major_freq}\n")

    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--alipath', required=True, type=str, dest="alipath",
                        help="Path to the input alignment file")
    parser.add_argument('--outpath', required=True, type=str, dest="outpath",
                        help="Path to the output text file")

    args = parser.parse_args()

    main(args.alipath, args.outpath)
