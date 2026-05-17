# Code for viral load analysis as part of my PhD thesis

This repository contains code for analyses of log<sub>10</sub> viral load data from Labor Berlin
used in my PhD thesis. The data are not included.

# Data

Viral load and associated data (date of PCR testing, PCR target, test center, individual's gender and date of birth) are provided in two files: `viral-load-with-negatives-first-positive-paper2.tsv` contains data from negative PCR tests and the first positive PCR test per infection. `min-3-timeseries-paper2-lines.json` contains viral load and associated data from consecutive PCR testing, including only infections with at least 3 PCR tests performed on different days.

# Code

The directory `sbatch_scripts` contains bash scripts that handle resource allocation and execution of MCMC sampling (directory `fit_models`) and posterior analysis (directory `post_pred`). Figures, tables, MCMC results are stored in the `output` directory. The `stan` directory contains the model specifications.

The directory `markdown` contains R markdown files for

* creating tables describing the study population and viral loads (`make_tables.Rmd`)
* evaluating the results from model fitting with real data (`real_data_evaluation.Rmd`)
* analyzing Rhat and ESS values from MCMC runs with simulated and real data (`rhat_evaluation.Rmd`)
* evaluating the results from model fitting with simulated data (`simulation_evaluation.Rmd`)
* analyzing estimated differences in log<sub>10</sub> viral load from first positive PCR tests between different groups of the study population (`stats_tests_first_vl.Rmd`)

The remaining files contain functions for 

* data processing (`data_utils.R`)
* figure generation (`plot_utils.R`)
* processing of MCMC results (`post_analysis.R`)
* figure generation of MCMC results (`post_analysis_plot.R`)


, scripts for 

* making swarmplots of log<sub>10</sub> viral loads from first positive PCR tests (`make_swarmplots_first_vl.R`)
* posterior analysis (`post_pred.R`)
* making plots of posterior results (`post_pred_all_models.R`)
* making posterior predictive check plots (`post_pred_check_all_models.R`)
* making prior predictive check plots (`prior_pred_plots.R`)
* model fitting (`run_stan_model.R`)
* estimating differences in log<sub>10</sub> viral load from first positive PCR tests between different groups of the study population (`stats_tests_first_vl.R`)

and specifications of

* data-related parameters (`data_params.R`)
* plotting-related parameters (`plot_params.R`)


# License

This project is licensed under the terms of the MIT license.




