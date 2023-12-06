#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition cpu
#SBATCH --mem 6G
#SBATCH --time 2-23:59:59


snakemake -k -s /scratch/tschnei6/seascapes/workflow/SnakeSimu -j 10000 --cluster "sbatch -J {params.name} -N 1 -o ./slurm/%x.%j.out -e ./slurm/%x.%j.err --cpus-per-task={params.threads} --mem={params.mem} --time={params.time}"
