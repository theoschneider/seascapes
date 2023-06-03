#!/bin/bash
# Set directory
folder=~/THÉO/seascapes

# Convert fasta to phy
python3 ${folder}/scripts/convert_fasta_to_phy.py \
        --input ${folder}/data/omm_NT_fasta.v10c_116/ENSG00000000003_TSPAN6_NT.fasta \
        --output ${folder}/processed/TSPAN6_NT.phy

for i in 1 2
do
  # Run mutsel
  ${folder}/utils/bayescode/bin/mutselomega --ncat 30 \
  -a ${folder}/processed/TSPAN6_NT.phy \
  -t ${folder}/data/omm_RooTree.v10b_116/ENSG00000000003_TSPAN6_NT.rootree \
  --until 2000 ${folder}/processed/mutsel_TSPAN6_${i}
done
