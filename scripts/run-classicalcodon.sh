#!/bin/bash

# Set directory
folder=~/THÉO/seascapes

${folder}/utils/bayescode/bin/mutselomega --freeomega --omegancat 30 --flatfitness \
                                          -a ${folder}/processed/TSPAN6_NT.phy \
                                          -t ${folder}/data/omm_RooTree.v10b_116/ENSG00000000003_TSPAN6_NT.rootree \
                                          -u 2000 ${folder}/classical_mutsel_TSPAN6
