#!/bin/bash
#
#SBATCH --job-name=first-vl-test-pams
#SBATCH --output=../output/output-first-vl-test-pams.txt
#
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --mem=128G

set -euo pipefail

cd $VL_THESIS_DIR
source ~/work/miniforge3/etc/profile.d/conda.sh
conda activate r-paper2-me

time Rscript stats_tests_first_vl.R --group pams
