#!/bin/bash
#
#SBATCH --job-name=vl-tc-model1-mintests6-min-pos-tests2-chains4-multi-threading-2000-samples-filter-by-T2-trim-neg-pcrs-exclude-lrt-samples-min-max-load5
#SBATCH --output=../output/output-model1-mintests6-min-pos-tests2-chains4-multi_threading-2000-samples-filter-by-T2-trim-neg-pcrs-exclude-lrt-samples-min-max-load5.txt
#
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --time=2-00
#SBATCH --mem=256G

set -euo pipefail

cd $VL_THESIS_DIR 
source ~/work/miniforge3/etc/profile.d/conda.sh
conda activate r-paper2-me

time Rscript bin/run_stan_model.R --model_stan 1 --threading --max_tree_depth 12 --min_tests 6 --min_pos_tests 2 --n_iter_warmup 1000 --n_iter_sampling 2000 --n_chains 4 --filter_by_testname "T2" --trim_neg_pcrs --exclude_lrt_samples --min_max_load 5
