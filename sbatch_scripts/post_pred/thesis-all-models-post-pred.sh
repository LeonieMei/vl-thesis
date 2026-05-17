#!/bin/bash
#
#SBATCH --job-name=thesis-all-models-post-pred
#SBATCH --output=../output/thesis-all-models-post-pred.txt
#
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem=250G

set -euo pipefail

cd $VL_THESIS_DIR
source ~/work/miniforge3/etc/profile.d/conda.sh
conda activate r-paper2-me

time Rscript post_pred_all_models.R --sampling_iterations 2000

