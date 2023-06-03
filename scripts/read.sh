#!/bin/bash
# Set directory
folder=~/THÉO/seascapes

# Read the 2 mutsel previously ran
${folder}/utils/bayescode/bin/readmutselomega --every 1 --until 2000 --burnin 1000 --ss ${folder}/processed/mutsel_TSPAN6_1
${folder}/utils/bayescode/bin/readmutselomega --every 1 --until 2000 --burnin 1000 --ss ${folder}/processed/mutsel_TSPAN6_2

# Compute the distance between the 2 runs
python3 ${folder}/scripts/distance.py --file1 ${folder}/processed/mutsel_TSPAN6_1.siteprofiles \
                                      --file2 ${folder}/processed/mutsel_TSPAN6_2.siteprofiles \
                                      --output ${folder}/results/
