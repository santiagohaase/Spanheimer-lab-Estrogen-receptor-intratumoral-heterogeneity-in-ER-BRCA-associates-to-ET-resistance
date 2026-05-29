#!/bin/bash
#SBATCH --job-name=tcga_erpos_ssgsea
#SBATCH --partition=allnodes
#SBATCH --time=UNLIMITED
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=10g
#SBATCH --output=/path/to/TCGA_data/ssgsea_out.txt
#SBATCH --error=/path/to/TCGA_data/ssgsea_err.txt
#SBATCH --chdir=/path/to/TCGA_data/

set -euo pipefail

BASE_DIR="/path/to/TCGA_data"

python "${BASE_DIR}/filter_tcga_er_positive.py" --base-dir "${BASE_DIR}"
python "${BASE_DIR}/run_ssgsea_tcga_erpos.py" --base-dir "${BASE_DIR}"
