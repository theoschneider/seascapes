#!/bin/bash

# Set directory
folder=~/THÉO/seascapes
outdir=${folder}/processed/rodents_subset/

mkdir -p ${outdir}

# Create subsets of species
python3 ${folder}/scripts/subset_species.py --tree ${folder}/data/omm_RooTree.v10b_116/ENSG00000000003_TSPAN6_NT.rootree \
                                            --fasta ${folder}/data/omm_NT_fasta.v10c_116/ENSG00000000003_TSPAN6_NT.fasta \
                                            --outdir ${outdir} \
                                            --subset "Dipodomys_ordii" "Cricetulus_griseus" "Cavia_porcellus" "Rattus_norvegicus" "Mus_musculus" "Mus_pahari" "Mus_spretus" "Mus_caroli" "Mesocricetus_auratus" "Peromyscus_maniculatus" "Meriones_unguiculatus" "Heterocephalus_glaber" "Octodon_degus" "Jaculus_jaculus" "Castor_canadensis" "Nannospalax_galili" "Marmota_marmota_marmota" "Ictidomys_tridecemlineatus" "Microtus_ochrogaster" "Chinchilla_lanigera" "Fukomys_damarensis"

# Convert fasta to phy and run mutsel 2 times for each subset
for subset in 1 2
do
  python3 ${folder}/scripts/convert_fasta_to_phy.py \
        --input ${outdir}/subset${subset}.fasta \
        --output ${outdir}/subset${subset}.phy
  for chain in 1 2
  do
    # Run mutsel
    ${folder}/utils/bayescode/bin/mutselomega --ncat 30 \
                                              -a ${outdir}/subset${subset}.phy \
                                              -t ${outdir}/subset${subset}.rootree \
                                              --until 2000 ${outdir}/mutsel_subset${subset}_${chain}
  done
done
