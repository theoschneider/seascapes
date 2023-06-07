#!/bin/bash
# Set directory
folder=~/THÉO/seascapes

until=2000
delta_t=100

for burnin in seq 0 ${until} ${delta_t}
do
  ${folder}/utils/bayescode/bin/readmutselomega --every 1 \
                                                --burnin ${burnin} \
                                                --until ${burnin}+${delta_t} \
                                                --ss ${folder}/processed/mutsel_TSPAN6_1

  ${folder}/utils/bayescode/bin/readmutselomega --every 1 \
                                                --burnin ${burnin} \
                                                --until ${burnin}+${delta_t} \
                                                --ss ${folder}/processed/mutsel_TSPAN6_2

  # Compute the distance between the 2 runs
  python3 ${folder}/scripts/distance.py --file1 ${folder}/processed/mutsel_TSPAN6_1.siteprofiles \
                                        --file2 ${folder}/processed/mutsel_TSPAN6_2.siteprofiles \
                                        --output ${folder}/results/
done
