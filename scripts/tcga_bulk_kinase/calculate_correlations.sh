#!/bin/bash
#SBATCH --job-name=tcga_erpos_esr1_correlations
#SBATCH --partition=allnodes
#SBATCH --time=UNLIMITED
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10g
#SBATCH --output=/path/to/TCGA_data/correlations_out.txt
#SBATCH --error=/path/to/TCGA_data/correlations_err.txt
#SBATCH --chdir=/path/to/TCGA_data/

set -euo pipefail

BASE_DIR="/path/to/TCGA_data"
python "${BASE_DIR}/compute_esr1_kinase_correlations.py" --base-dir "${BASE_DIR}"
