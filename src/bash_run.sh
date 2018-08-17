#!/usr/bin/env bash

#SBATCH -J disordered_protein
#SBATCH -p DPB
#SBATCH -c 4
#SBATCH --mem=6000

module load Python/3.6.0
srun python disorderly.py "$@"
