#!/bin/bash
# Set directory
folder=~/THÉO/seascapes

# Read the 2 mutsel previously ran
${folder}/utils/bayescode/bin/readmutselomega --every 1 --until 2000 --burnin 1000 --ss ${folder}/processed/mutsel_TSPAN6_1
${folder}/utils/bayescode/bin/readmutselomega --every 1 --until 2000 --burnin 1000 --ss ${folder}/processed/mutsel_TSPAN6_2

# moyenne = premiere ligne commentées
# toute les valeurs en csv en dessous