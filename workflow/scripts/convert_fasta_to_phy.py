import argparse
from libraries import open_fasta
from libraries import write_ali


def main(input_path: str, output_path: str):
    fasta_seq = open_fasta(input_path)
    write_ali(fasta_seq, output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', required=True, type=str, dest="input_path", help="Input fasta file path")
    parser.add_argument('--output', required=True, type=str, dest="output_path", help="Output phy file path")
    args = parser.parse_args()
    main(args.input_path, args.output_path)
