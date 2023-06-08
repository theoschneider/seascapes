#!/bin/bash
# Set directory
folder=~/THÉO/seascapes

# Variables that can be changed
loop_until=2000 # This is the --until parameter in readmutselomega
delta=200

python3 "${folder}/scripts/autocor.py" --path "${folder}/processed/TSPAN6_${delta}/" \
                                       --outdir "${folder}/results/"
