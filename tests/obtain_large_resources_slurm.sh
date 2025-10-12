#!/bin/bash
#SBATCH --partition batch
#SBATCH --time 48:00:00
#SBATCH -c 10
#SBATCH --mem 128G

srun bash obtain_large_resources.sh