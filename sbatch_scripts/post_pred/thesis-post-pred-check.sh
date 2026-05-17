#!/bin/bash
#
#SBATCH --job-name=thesis-post-pred-check
#SBATCH --output=../output/thesis-post-pred-check.txt
#
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem=200G

set -euo pipefail

cd $VL_THESIS_DIR
source ~/work/miniforge3/etc/profile.d/conda.sh
conda activate r-paper2-me

time Rscript post_pred_check_all_models.R --sampling_iterations 2000

