#!/bin/bash
# Set directory
folder=~/THÉO/seascapes

# Read MutSel for all 4 files
for i in 1_1 1_2 2_1 2_2
do

  ${folder}/utils/bayescode/bin/readmutselomega --every 1 \
                                                --burnin 1000 \
                                                --until 2000 \
                                                --ss "${folder}/processed/rodents_subset/mutsel_subset${i}"

done


# Compute the distance between all the .siteprofiles in the folder
python3 "${folder}/scripts/distance.py" --file1 "${folder}/processed/rodents_subset/mutsel_subset1_1.siteprofiles" \
                                        --file2 "${folder}/processed/rodents_subset/mutsel_subset1_2.siteprofiles" \
                                        --file3 "${folder}/processed/rodents_subset/mutsel_subset2_1.siteprofiles" \
                                        --file4 "${folder}/processed/rodents_subset/mutsel_subset2_2.siteprofiles" \
                                        --outdir "${folder}/results/rodents_subset/"
