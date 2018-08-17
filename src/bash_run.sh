#!/usr/bin/env bash

#SBATCH -J disordered_protein
#SBATCH -p DPB
#SBATCH -c 24
#SBATCH --mem=10000

module load Python/3.6.0
srun python disorderly "$@"
