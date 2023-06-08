#!/bin/bash
# Set directory
folder=~/THÉO/seascapes

# Variables that can be changed
loop_until=2000 # This is the --until parameter in readmutselomega
delta=200

# From 0 to "loop_until" with a step of "delta"
for ((burnin=0; burnin<=loop_until-delta; burnin+=delta))
do
  mkdir -p "${folder}/processed/TSPAN6_${delta}"

  ${folder}/utils/bayescode/bin/readmutselomega --every 1 \
                                                --burnin $burnin \
                                                --until $((burnin+delta)) \
                                                --ss "${folder}/processed/mutsel_TSPAN6_1"

  mv "${folder}/processed/mutsel_TSPAN6_1.siteprofiles" "${folder}/processed/TSPAN6_${delta}/${burnin}_$((burnin+delta)).siteprofiles"
done

# Compute the distance between all the .siteprofiles in the folder
python3 "${folder}/scripts/autocor.py" --path "${folder}/processed/TSPAN6_${delta}/" \
                                       --outdir "${folder}/results/"
