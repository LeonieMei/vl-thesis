#!/bin/bash
#
#SBATCH --job-name=vl-tc-model1-chains4-multi-threading-sim-data-all-params-zero-filter-by-T2-trim-neg-pcrs-exclude-lrt-samples-min-max-load5
#SBATCH --output=../output/output-model1-chains4-multi_threading-sim-data-all-params-zero-filter-by-T2-trim-neg-pcrs-exclude-lrt-samples-min-max-load5.txt
#
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --time=4-00
#SBATCH --mem=128G

set -euo pipefail

cd $VL_THESIS_DIR 
source ~/work/miniforge3/etc/profile.d/conda.sh
conda activate r-paper2-me


for min_tests in {3..7}; do
    for min_pos_tests in {2..7}; do
      if [ "$min_pos_tests" -gt "$min_tests" ] ; then
        continue
      fi
      time Rscript bin/run_stan_model.R --model_stan 1 --threading --max_tree_depth 12 --min_tests ${min_tests} --min_pos_tests ${min_pos_tests} --n_iter_warmup 1000 --n_iter_sampling 1000 --n_chains 4 --simulate_data --simulate_param_sizes_all_zero --filter_by_testname "T2" --trim_neg_pcrs --exclude_lrt_samples --min_max_load 5
  done
done
