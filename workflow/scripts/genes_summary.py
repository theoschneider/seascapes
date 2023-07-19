import os
import pandas as pd
from ete3 import Tree

folder = os.path.expanduser("~/THÉO/seascapes/")


def find_tree_length(input_tree):
    out_tree_length = 0
    for node in input_tree.traverse():
        out_tree_length += node.dist
    return out_tree_length


# Read the species list
df = pd.read_csv(folder + "data/species.tsv", sep="\t", header=0, index_col=None)

# Convert the dataframe to a dictionary
sp_dic = {df.iloc[i, 0]: df.iloc[i, 1] for i in range(len(df))}
interest = {species for species in sp_dic.keys() if sp_dic[species] in {"Glires", "Scandentia", "Dermoptera", "Primates"}}

# Create empty out dictionary
out = {"gene": [],
       "tree_length": [],
       "tree_species": [],
       "s1_length": [],
       "s2_length": [],
       "s1_species": [],
       "s2_species": []}

# List all tree files
tree_files = os.listdir(folder + "data/omm_RooTree.v10b_116/")


for file in tree_files:
    tree = Tree(folder + "data/omm_RooTree.v10b_116/" + file, format=1)

    for species in tree.get_leaf_names():
        assert species in sp_dic.keys(), f"{species} not in the species list!"

    toKeep = set(tree.get_leaf_names()) & interest
    notKeep = set(tree.get_leaf_names()) - toKeep

    if len(toKeep) == 0 or len(notKeep) == 0:
        continue

    tree_length = find_tree_length(tree)

    subset1 = tree.copy()
    subset1.prune(toKeep, preserve_branch_length=True)

    subset1_length = find_tree_length(subset1)

    subset2 = tree.copy()
    subset2.prune(notKeep, preserve_branch_length=True)

    subset2_length = find_tree_length(subset2)

    out["gene"].append(file.split("_NT")[0])
    out["tree_length"].append(tree_length)
    out["tree_species"].append(len(tree.get_leaf_names()))
    out["s1_length"].append(subset1_length)
    out["s2_length"].append(subset2_length)
    out["s1_species"].append(len(subset1.get_leaf_names()))
    out["s2_species"].append(len(subset2.get_leaf_names()))

out_df = pd.DataFrame.from_dict(out)
out_df.to_csv(folder + "data/genes_summary.csv", header=True, index=False)
