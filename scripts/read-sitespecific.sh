#!/bin/bash

# Set directory
folder=~/THÉO/seascapes


# Compute site omega 0
${folder}/utils/bayescode/bin/readmutselomega --every 1 --until 2000 --burnin 1000 --confidence_interval 0.025 \
                                              --omega_0 ${folder}/processed/mutsel_TSPAN6_1

# Compute site omega 1

