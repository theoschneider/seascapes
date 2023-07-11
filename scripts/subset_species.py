import argparse
from ete3 import Tree
import pandas as pd


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

    # Split the fasta files in 2 (according to subset argument)
    with open(fasta_path, 'r') as f:
        fasta = f.readlines()

    fasta1 = ""
    fasta2 = ""


    for i in range(0, len(fasta), 2):

        sp = fasta[i][1:].strip()
        seq = fasta[i+1].strip()

        # Check that there are no empty sequences
        assert seq != len(seq) * "-", "One or more sequences are empty!"

        # Look for species name in the tree file and append the corresponding fasta sequence to the right file
        if sp in tree1.get_leaf_names():
            fasta1 += ">" + sp + "\n" + seq + "\n"
        else:
            fasta2 += ">" + sp + "\n" + seq + "\n"

    # Check that there are at least 20 species in each subset
    assert len(tree1.get_leaf_names()) >= 20 and len(tree2.get_leaf_names()) >= 20, \
        "Both subsets must contain at least 20 species!"

    # Write the 2 subtrees and 2 fasta files
    with open(outdir + "subset1.rootree", 'w') as f:
        f.write(tree1.write(format=1))

    with open(outdir + "subset2.rootree", 'w') as f:
        f.write(tree2.write(format=1))

    with open(outdir + "subset1.fasta", 'w') as f:
        f.write(fasta1)

    with open(outdir + "subset2.fasta", 'w') as f:
        f.write(fasta2)


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
