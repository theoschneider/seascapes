# Set directory
folder=~/Users/theo/desktop/uni/M/seascapes/

# Convert fasta to phy
python3 ${folder}/scripts/convert_fasta_to_phy.py --input ${folder}/data/omm_NT_fasta.v10c_116/ENSG00000000003_TSPAN6_NT.fasta --output ${folder}processed/TSPAN6_NT.phy

# Run mutsel
mutselomega --ncat 30 -a ${folder}/processed/TSPAN6_NT.phy -t ${folder}/data/omm_RooTree.v10b_116/ENSG00000000003_TSPAN6_NT.rootree --until 2000 run_mutsel_TSPAN6
