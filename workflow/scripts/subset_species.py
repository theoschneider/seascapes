import argparse
from ete3 import Tree
import pandas as pd
from libraries import open_fasta, write_fasta, filter_fasta


def main(tree_path: str, fasta_path: str, subset_path: str, outdir: str):

    # Read the subset file
    subset = pd.read_csv(subset_path, header=None, index_col=None)
    subset = subset[0].tolist()

    # Read the tree
    tree1 = Tree(tree_path, format=1)
    tree2 = tree1.copy()

    # Split the tree in 2 subtrees (according to subset argument)
    tree1.prune(subset, preserve_branch_length=True)
    tree2.prune(list(set(tree2.get_leaf_names()) - set(subset)), preserve_branch_length=True)

    # Open the 2 fasta and split them
    fasta1 = open_fasta(fasta_path, filter_species=subset)
    fasta2 = open_fasta(fasta_path, filter_species=list(tree2.get_leaf_names()))

    fasta1, fasta2 = filter_fasta(fasta1, fasta2)

    # Check that there are no empty sequences
    for seq_id in fasta1.keys():
        assert fasta1[seq_id] != len(fasta1[seq_id]) * "-", f"The sequence {seq_id} is empty!"

    for seq_id in fasta2.keys():
        assert fasta2[seq_id] != len(fasta2[seq_id]) * "-", f"The sequence {seq_id} is empty!"

    # Check that there are at least 20 species in each subset
    assert len(tree1.get_leaf_names()) >= 20 and len(tree2.get_leaf_names()) >= 20, \
        "Both subsets must contain at least 20 species!"

    # Write the 2 subtrees and 2 fasta files
    with open(outdir + "subset1.rootree", 'w') as f:
        f.write(tree1.write(format=1))

    with open(outdir + "subset2.rootree", 'w') as f:
        f.write(tree2.write(format=1))

    write_fasta(fasta1, outdir + "subset1.fasta")
    write_fasta(fasta2, outdir + "subset2.fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--tree', required=True, type=str, dest="tree_path",
                        help="Path to the tree file")
    parser.add_argument('--fasta', required=True, type=str, dest="fasta_path",
                        help="Path to the fasta file")
    parser.add_argument('--outdir', required=True, type=str, dest="outdir",
                        help="Path to the output directory")
    parser.add_argument('--subset', required=True, nargs='*', dest="subset_path",
                        help="csv file containing the species to keep")

    args = parser.parse_args()
    main(args.tree_path, args.fasta_path, args.subset, args.outdir)
