library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(optparse)
library(here)

TOP <- here()
source(file.path(TOP, "data_params.R"))
source(file.path(TOP, "post_analysis_plot.R"))

MAIN_PARAMS_MEAN_STATS <- file.path(DIR_PDATA, "params_stats.tsv")
MAIN_PARAMS_MEAN_STATS_POST_PRED_GROUPS2 <- file.path(DIR_PDATA, "params_stats_post_pred_groups2.tsv")
MAIN_PARAMS_MEAN_DIFF_STATS <- file.path(DIR_PDATA, "params_diff_stats.tsv")
MAIN_PARAMS_RHAT_STATS <- file.path(DIR_PDATA, "params_rhat_stats.tsv")

option_list <- list(
  make_option(c("-s", "--model_stan"),
    type = "integer", default = 1,
    help = "stan model number (see see .stan files in CP directory)"
  ),
  make_option(c("--min_tests", type = "integer"),
    default = 3,
    help = "Minimum number of tests per infection."
  ),
  make_option(c("--max_tests", type = "integer"),
    default = 20,
    help = "Maximum number of tests per infection."
  ),
  make_option(c("--min_pos_tests", type = "integer"),
    default = 2,
    help = "Minimum number of positive tests per infection."
  ),
  make_option(c("--n_iter_sampling", type = "integer"),
    default = 750,
    help = "Number of iterations in the MCMC sampling phase."
  ),
  make_option(c("--n_chains", type = "integer"),
    default = 8,
    help = "Number of parallel MCMC chains."
  ),
  make_option(c("--include_warmup", type = "logical"),
    action = "store_true", default = FALSE,
    help = "Include warmup samples when loading posterior samples."
  ),
  make_option(c("--testing", type = "logical"),
    action = "store_true", default = FALSE,
    help = "Testing (with a subset of samples)."
  ),
  make_option(c("--simulate_data"),
    default = FALSE, action = "store_true",
    help = "Simulated viral load data."
  ),
  make_option(c("--simulate_data_from_prior"),
    default = FALSE, action = "store_true",
    help = "Simulated viral load data from model's prior predictions was used."
  ),
  make_option(c("--simulate_param_sizes_all_zero"),
    default = FALSE, action = "store_true",
    help = "Simulate viral load data with all param sizes (except for main model params and infection random effects params) set to 0."
  ),
  make_option(c("--n_samples", type = "integer"),
    default = 0,
    help = "The number of infections viral load data was simulated for if option --simulate_data is set or the number of infections randomly drawn from the actual data."
  ),
  make_option(c("--simulate_mixed_courses"),
    default = FALSE, action = "store_true",
    help = "Simulate infection courses from one of two patterns: 1) a sharp peak and a relatively slow clearance or 2) a smooth transition at the peak and a relatively fast clearance."
  ),
  make_option(c("--remove_problematic_infections"),
    default = FALSE, action = "store_true",
    help = "Removed infections whose Rhat values in the time_shift parameter (b_shift) were above 1.01 when running model17_17 in May 2024."
  ),
  make_option(c("--filter_by_testname"),
    default = NULL,
    help = "Only use PCRs with the given testname (LC480, T2, T3, T6, T9, Ta)."
  ),
  make_option(c("--suffix", type = "character"),
    default = "",
    help = "Suffix added to the end of filenames that are created from the wanted run."
  ),
  make_option(c("--make_only_mean_figure"),
    default = FALSE, action = "store_true",
    help = "Make only a figure showing predicted mean viral load courses."
  ),
  make_option(c("--make_only_group_figures"),
    default = FALSE, action = "store_true",
    help = "Make only the figures showing predicted mean viral load courses by subgroup."
  ),
  make_option(c("--collect_only_rhat_stats"),
    default = FALSE, action = "store_true",
    help = "Collect only Rhat statistics as opposed to making figures."
  ),
  make_option(c("--threading"),
    default = FALSE, action = "store_true",
    help = "For model fitting multiple threads per chain were run."
  ),
  make_option(c("--use_true_vls"),
    default = FALSE, action = "store_true",
    help = "When data was simulated, the true viral loads (including negative viral loads) where used for fitting the model."
  ),
  make_option(c("--trim_neg_pcrs"),
    default = FALSE, action = "store_true",
    help = "Trim all negative PCRs that happened more than days_to_negative days before/after a positive PCR from the same infection."
  ),
  make_option(c("--days_to_negative", type = "integer"),
    default = 14,
    help = "Number of days between a negative and a first positive or a negative and a last positive PCR for the negative PCR to still be included in the infection. For now, this is only effective when simulating data."
  ),
  make_option(c("--single_leading_trailing_neg_pcr"),
    default = FALSE, action = "store_true",
    help = "Of all negative PCRs in an infection, keep only the last leading negative PCR (before the first positive PCR) and the first trailing negative PCR (after the last positive PCR)."
  ),
  make_option(c("--no_vl_error"),
    default = FALSE, action = "store_true",
    help = "No error was added to the observed viral loads (i.e. no jittering) when simulating data."
  ),
  make_option(c("--error_dist"),
    default = "normal", default = "normal",
    help = "When simulating data, add error from the specified distribution to the observed viral loads."
  ),
  make_option(c("--add_neg_pcrs"),
    default = FALSE, action = "store_true",
    help = "When simulating data, (as of now) 7 negative PCRs were added to the beginning and/or end of an infection, if there is already a negative PCR at the beginning/end of that infection."
  ),
  make_option(c("--no_shifts"),
    default = FALSE, action = "store_true",
    help = "When simulating data, the PCR time series was not shifted in time before fitting."
  ),
  make_option(c("--only_pos_pcrs"),
    default = FALSE, action = "store_true",
    help = "When simulating data, only positive PCR tests were used (then, model with truncated normal likelihood was fitted)."
  ),
  make_option(c("--color_by"),
    default = "imputed", type = "character",
    help = "Color shifted data points according to whether Rhat of the shift parameter for the corresponding infection is above 1.05."
  ),
  make_option(c("--adapt_delta", type = "real"),
    default = .8,
    help = "Set the adapt_delta (acceptance rate) for MCMC sampling."
  ),
  make_option(c("--use_paper1_data"),
    default = FALSE, action = "store_true",
    help = "Use only the data from the first viral load paper (before April 3rd 2021)."
  ),
  make_option(c("--exclude_lrt_samples"),
    default = FALSE, action = "store_true",
    help = "Exclude PCRs from lower respiratory tract (LRT) samples."
  ),
  make_option(c("--exclude_too_steep_courses"),
    default = FALSE, action = "store_true",
    help = "Infections with a slope above 5 (or below -5) between two PCRs were excluded."
  ),
  make_option(c("--max_diff_load_perday", type = "real"),
    default = 5,
    help = "The maximum increase in viral load between the first positive and its preceeding negative PCR."
  ),
  make_option(c("--min_duration", type = "integer"),
    default = 0,
    help = "The minimum number of days between the earliest and latest positive PCR in an infection."
  ),
  make_option(c("--min_max_load", type = "real"),
    default = 4,
    help = "The maximum log10 viral load of an infection needs to be at least min_max_load."
  ),
  make_option(c("--plot_lines"),
    default = FALSE, action = "store_true",
    help = "Plot lines instead of points for shifted infection data."
  ),
  make_option(c("--plot_n_infections", type = "integer"),
    default = 0,
    help = "The number of infections to sample for plotting."
  ),
  make_option(c("--summary_stat"),
    default = "median",
    help = "The summary statistic (median or mean) to use for summarizing the main model parameters (intercept, up-slope, down-slope)."
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

output_file <- MAIN_PARAMS_MEAN_STATS

assert_that(opt$model_stan %in% c(1, 2, 3, 4, 5), msg = "You can choose from stan models 1, 2, 3, 4 and 5.")

if (opt$color_by == "is.difficult") {
  assertthat::assert_that(opt$no_shifts == FALSE)
}

if (opt$n_samples == 0) {
  opt$n_samples <- NULL
}

if (opt$plot_n_infections == 0) {
  opt$plot_n_infections <- NULL
}

if (opt$use_true_vls == TRUE) {
  assertthat::assert_that(opt$simulate_data == TRUE,
    msg = "You can only set option --use_true_vls together with option --simulate_data"
  )
}

assert_that(opt$summary_stat %in% c("median", "mean"), msg = "You can only use median or mean as a summary statistic.")

model_suffix_stan <- paste0(opt$model_stan, ifelse(opt$threading == TRUE, "_threading", ""))
opt_suffix_samples <- ifelse(opt$simulate_data == TRUE, paste0("_simvl", opt$n_samples),
  ifelse(!is.null(opt$n_samples), paste0("_samples", opt$n_samples), "")
)
opt_suffix_simulate_param_sizes_all_zero <- ifelse(opt$simulate_param_sizes_all_zero == TRUE, "_param_sizes_all_zero", "")
opt_suffix_simulate_mixed_courses <- ifelse(opt$simulate_mixed_courses, "_mixed_courses", "")
opt_suffix <- ifelse(opt$suffix == "", "", paste0("_", opt$suffix))
opt_suffix_remove_problematic_infections <- ifelse(opt$remove_problematic_infections == TRUE, "_wo_probl_inf", "")
opt_testname_suffix <- ifelse(!is.null(opt$filter_by_testname), paste0("_", opt$filter_by_testname), "")
opt_use_true_vls_suffix <- ifelse(opt$use_true_vls == TRUE, "_true_vls", "")
opt_trim_neg_pcrs <- ifelse(opt$trim_neg_pcrs == TRUE, "_trim_neg_pcrs", "")
opt_days_to_negative <- ifelse((opt$trim_neg_pcrs == TRUE) | (opt$add_neg_pcrs == TRUE) | (opt$days_to_negative != 14), paste0("_", opt$days_to_negative, "_days_to_negative"), "")
opt_single_leading_trailing_neg_pcr <- ifelse(opt$single_leading_trailing_neg_pcr, "_single_surrounding_neg_pcrs", "")
opt_no_vl_error <- ifelse(opt$no_vl_error == TRUE, "_no_vl_error", "")
opt_error_dist <- ifelse(opt$error_dist == "normal" | opt_no_vl_error == TRUE, "", paste0("_", gsub("-", "_", opt$error_dist)))
opt_add_neg_pcrs <- ifelse(opt$add_neg_pcrs == TRUE, "_add_neg_pcrs", "")
opt_add_no_shifts <- ifelse(opt$no_shifts == TRUE, "_no_shifts", "")
opt_only_pos_pcrs <- ifelse(opt$only_pos_pcrs == TRUE, "_only_pos_pcrs", "")
opt_max_tests <- ifelse(opt$max_tests != 20, paste0("_maxtests", opt$max_tests), "")
opt_min_pos_tests <- ifelse(opt$min_pos_tests != 2, paste0("_min_pos_tests", opt$min_pos_tests), "")
opt_simulate_data_from_prior <- ifelse(opt$simulate_data_from_prior == TRUE, "_from_prior", "")
opt_adapt_delta <- ifelse(opt$adapt_delta == 0.8, "",
  paste0("_adapt_delta_", sub("\\.", "_", opt$adapt_delta))
)
opt_use_paper1_data <- ifelse(opt$use_paper1_data, "_paper1_data", "")
opt_exclude_lrt_samples <- ifelse(opt$exclude_lrt_samples, "_exclude_lrt_samples", "")
opt_exclude_too_steep_courses <- ifelse(opt$exclude_too_steep_courses, "_exclude_too_steep_courses", "")
opt_max_diff_load_perday <- ifelse(opt$exclude_too_steep_courses, opt$max_diff_load_perday, "")
opt_min_duration <- ifelse(opt$min_duration, paste0("_min_", opt$min_duration, "_days_duration"), "")
opt_min_max_load <- paste0("_min_max_load_", opt$min_max_load)
opt_summary_stat <- ifelse(opt$summary_stat == "median", "_median", "")


if (opt$model_stan %in% MODELS_DAY2) {
  day_var <- "day2"
} else {
  day_var <- "day"
}


model_dir <- paste0("model", opt$model_stan)
dir_figures <- file.path(DIR_FIGURES, "vl_trajectories", "real", model_dir)
dir_pdata <- file.path(DIR_PDATA, model_dir)

if (!dir.exists(dir_figures)) {
  dir.create(dir_figures)
}
if (!dir.exists(dir_pdata)) {
  dir.create(dir_pdata)
}

model <- paste0(
  "model",
  model_suffix_stan,
  "_sel",
  opt$min_tests,
  opt_max_tests,
  opt_min_pos_tests,
  "_chains",
  opt$n_chains,
  "_n_iter",
  opt$n_iter_sampling,
  opt_suffix_samples,
  opt_simulate_data_from_prior,
  opt_suffix_simulate_param_sizes_all_zero,
  opt_suffix_simulate_mixed_courses,
  opt_suffix_remove_problematic_infections,
  opt_testname_suffix,
  opt_use_true_vls_suffix,
  opt_trim_neg_pcrs,
  opt_days_to_negative,
  opt_single_leading_trailing_neg_pcr,
  opt_no_vl_error,
  opt_error_dist,
  opt_add_neg_pcrs,
  opt_add_no_shifts,
  opt_only_pos_pcrs,
  opt_adapt_delta,
  opt_use_paper1_data,
  opt_exclude_lrt_samples,
  opt_exclude_too_steep_courses,
  opt_max_diff_load_perday,
  opt_min_duration,
  opt_min_max_load,
  opt_suffix
)

run_data <- load_draws(
  model = model,
  model_dir = model_dir,
  model_no_stan = opt$model_stan,
  inc_warmup = opt$include_warmup,
  testing = opt$testing,
  test_thin = 1
)

pp_results <- make_post_pred_all(
  draws = run_data$draws,
  day_data = run_data$day_data,
  model = model,
  dir_pdata = dir_pdata,
  day_var = day_var,
  model_no_stan = opt$model_stan,
  ss = run_data$ss,
  stat = opt$summary_stat,
  testing = opt$testing
)


# Make plot showing predicted mean viral load courses.
if (opt$simulate_data) {
  sim_data <- run_data$sim_data
  output_file2 <- MAIN_PARAMS_MEAN_DIFF_STATS
} else {
  sim_data <- NULL
  output_file2 <- NULL
}


collect_rhat_stats(
  model = model,
  model_no_stan = opt$model_stan,
  ss = run_data$ss,
  output_file = MAIN_PARAMS_RHAT_STATS,
  min_tests = opt$min_tests,
  simulated = opt$simulate_data
)


if (!opt$make_only_group_figures && !opt$collect_only_rhat_stats) {
  make_post_pred_summary_stat(
    VLCP_by_day_draws = pp_results$VLCP_by_day_draws,
    VLCP_by_day_ID = pp_results$VLCP_by_day_ID,
    ss = run_data$ss,
    shift_draws = pp_results$shift_draws,
    shifted_data_by_draw_day_ID = pp_results$shifted_data_by_draw_day_ID,
    draws = run_data$draws,
    day_data = run_data$day_data,
    dir_pdata = dir_pdata,
    dir_figures = dir_figures,
    model = model,
    n_min_tests = opt$min_tests,
    n_samples = opt$n_samples,
    model_no_stan = opt$model_stan,
    day_var = day_var,
    sim_data = sim_data,
    testing = opt$testing,
    annotate_plot = FALSE,
    output_file = output_file,
    output_file2 = output_file2,
    color_by = opt$color_by,
    plot_lines = opt$plot_lines,
    plot_n_infections = opt$plot_n_infections,
    stat = opt$summary_stat
  )
}

if (!opt$make_only_mean_figure && !opt$collect_only_rhat_stats) {
  # Make predictions based on grouping variables.
  make_post_pred_by_group(
    VLCP_by_day_draws = pp_results$VLCP_by_day_draws,
    draws = run_data$draws,
    day_data = run_data$day_data,
    dir_figures = dir_figures,
    model_no_stan = opt$model_stan,
    model = model,
    output_file = MAIN_PARAMS_MEAN_STATS_POST_PRED_GROUPS2
  )

  ggsave(pp_results$ppc_timecourse,
    file = file.path(dir_figures, paste0("ppc_timecourse_", model, ".png")),
    width = COL_WIDTH * 2, height = COL_WIDTH * 1.5, units = "in", dpi = 300,
    type = "cairo"
  )
  ggsave(pp_results$ppc_timecourse,
    file = file.path(dir_figures, paste0("ppc_timecourse_", model, ".pdf")),
    width = COL_WIDTH * 2, height = COL_WIDTH * 1.5, units = "in"
  )
}
