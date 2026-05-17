#!/bin/bash
#
#SBATCH --job-name=vl-tc-model3-chains4-multi-threading-sim-data-all-params-zero-filter-by-T2-trim-neg-pcrs-7-days-to-negative-exclude-lrt-samples-min-max-load5-post-pred
#SBATCH --output=../output/output-model3-chains4-multi_threading-sim-data-all-params-zero-filter-by-T2-trim-neg-pcrs-7-days-to-negative-exclude-lrt-samples-min-max-load5-post-pred.txt
#
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --time=1-00
#SBATCH --mem=128G

set -euo pipefail

cd $VL_THESIS_DIR 
source ~/work/miniforge3/etc/profile.d/conda.sh
conda activate r-paper2-me

for min_tests in {5..7}; do
  for min_pos_tests in {2..7}; do
    if [ "$min_pos_tests" -gt "$min_tests" ] ; then
      continue
    fi
    time Rscript post_pred.R --model_stan 3 --threading --min_tests ${min_tests} --min_pos_tests ${min_pos_tests} --n_iter_sampling 1000 --n_chains 4 --simulate_data --simulate_param_sizes_all_zero --filter_by_testname "T2" --trim_neg_pcrs --days_to_negative 7 --single_leading_trailing_neg_pcr --exclude_lrt_samples --min_max_load 5 --make_only_mean_figure --color_by "is.difficult"
  done
done
