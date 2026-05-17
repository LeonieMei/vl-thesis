#!/bin/bash
#
#SBATCH --job-name=vl-tc-model1-mintests6-2000-samples-chains4-multi-threading-filter-by-T2-trim-neg-pcrs-exclude-lrt-samples-min-max-load5-post-pred
#SBATCH --output=../output/output-model1-mintests6-2000-samples-chains4-multi_threading-filter-by-T2-trim-neg-pcrs-exclude-lrt-samples-min-max-load5-post-pred.txt
#
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --mem=128G

set -euo pipefail

cd $VL_THESIS_DIR 
source ~/work/miniforge3/etc/profile.d/conda.sh
conda activate r-paper2-me

time Rscript post_pred.R --model_stan 1 --threading --min_tests 6 --min_pos_tests 2 --n_iter_sampling 2000 --n_chains 4 --filter_by_testname "T2" --trim_neg_pcrs --exclude_lrt_samples --min_max_load 5 --color_by "is.difficult"