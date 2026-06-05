library(data.table)
library(here)
library(optparse)
library(cmdstanr)
library(posterior)
setDTthreads(1)


TOP <- here()
source(file.path(TOP, "sars2vl", "data_params.R"))
source(file.path(TOP, "sars2vl", "data_funcs.R"))
CMD_STAN_DIR <- Sys.getenv("CMD_STAN_DIR")
set_cmdstan_path(CMD_STAN_DIR)


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
  make_option(c("--n_iter_warmup", type = "integer"),
    default = 1000,
    help = "Number of iterations in the MCMC warm-up phase."
  ),
  make_option(c("--n_iter_sampling", type = "integer"),
    default = 1000,
    help = "Number of iterations in the MCMC sampling phase."
  ),
  make_option(c("--n_chains", type = "integer"),
    default = 4,
    help = "Number of parallel MCMC chains."
  ),
  make_option(c("--max_tree_depth", type = "integer"),
    default = 10,
    help = "Maximum tree depth during MCMC sampling."
  ),
  make_option(c("--suffix", type = "character"),
    default = "",
    help = "Suffix to add to the end of filenames that are created from this run."
  ),
  make_option(c("--tc_file", type = "character"),
    default = TC_FILE_JSON,
    help = paste("File containing time course data. Note: only the filename needs to be provided, assuming that the file sits in", DIR_DATA, ".")
  ),
  make_option(c("--vl_first_pos_file", type = "character"),
    default = VL_FIRST_POS_FILE_TSV,
    help = paste("File containing first positive viral load data. Note: only the filename needs to be provided, assuming that the file sits in", DIR_DATA, ".")
  ),
  make_option(c("--simulate_data"),
    default = FALSE, action = "store_true",
    help = "Simulate viral load data."
  ),
  make_option(c("--simulate_data_from_prior"),
    default = FALSE, action = "store_true",
    help = "Simulate viral load data from model's prior predictions."
  ),
  make_option(c("--simulate_param_sizes_all_zero"),
    default = FALSE, action = "store_true",
    help = "Simulate viral load data with all param sizes (except for main model params and infection random effects params) set to 0."
  ),
  make_option(c("--n_samples", type = "integer"),
    default = 0, # Note: using NA or NULL here results in some weird behaviour
    help = "The number of infections to simulate viral load data for if option --simulate_data is set or the number of infections to randomly draw from the actual data."
  ),
  make_option(c("--simulate_mixed_courses"),
    default = FALSE, action = "store_true",
    help = "Simulate infection courses from one of two patterns: 1) a sharp peak and a relatively slow clearance or 2) a smooth transition at the peak and a relatively fast clearance."
  ),
  make_option(c("--remove_problematic_infections"),
    default = FALSE, action = "store_true",
    help = "Remove infections whose Rhat values in the time_shift parameter (b_shift) were above 1.01 when running model17_17 in May 2024."
  ),
  make_option(c("--filter_by_testname"),
    default = NULL,
    help = "Only use PCRs with the given testname LC480, T2, T3, T6, T9, Ta for the E-gene and T1, T4, T5, T9, Tb for the N-gene."
  ),
  make_option(c("--threading"),
    default = FALSE, action = "store_true",
    help = "For model fitting use multiple threads per chain."
  ),
  make_option(c("--use_true_vls"),
    default = FALSE, action = "store_true",
    help = "When simulating data, use the true viral loads (i.e. including negative log10 viral load values) when fitting the model."
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
    help = "When simulating data, don't add any error to the observed viral loads (i.e. no jittering)."
  ),
  make_option(c("--error_dist"),
    default = "normal",
    help = "When simulating data, add error from the specified distribution to the observed viral loads."
  ),
  make_option(c("--compile_only"),
              default = FALSE, action = "store_true",
              help = "Only compile stan models (don't run)."
  ),
  make_option(c("--add_neg_pcrs"),
    default = FALSE, action = "store_true",
    help = "When simulating data, add (as of now) 7 negative PCRs at the beginning and/or end of an infection, if there is already a negative PCR at the beginning/end of that infection."
  ),
  make_option(c("--no_shifts"),
    default = FALSE, action = "store_true",
    help = "When simulating data, don't shift the PCR time series before fitting."
  ),
  make_option(c("--only_pos_pcrs"),
    default = FALSE, action = "store_true",
    help = "When simulating data, use only positive PCR tests (then run model with truncated normal likelihood)."
  ),
  make_option(c("--adapt_delta", type = "real"),
    default = .8,
    help = "Set the adapt_delta (acceptance rate) for MCMC sampling."
  ),
  make_option(c("--overall_vl_error"),
    default = FALSE, action = "store_true",
    help = "Simulate data using the same viral load error (sigma) for every data point instead of infection-wise errors."
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
    help = "Exclude infections with a slope above max_diff_load_perday between two PCRs."
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
  make_option(c("--process_rds_file"),
    default = FALSE, action = "store_true",
    help = "Instead of fitting, load and process the existing .RDS file with fitting results."
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

assert_that(opt$model_stan %in% c(1, 2, 3, 4, 5), msg = "You can choose from stan models 1, 2, 3, 4 and 5.")

options(mc.cores = 4)

if (opt$n_samples == 0) {
  opt$n_samples <- NULL
}

if (opt$use_true_vls == TRUE) {
  assert_that(opt$simulate_data == TRUE,
    msg = "You can only set option --use_true_vls together with option --simulate_data"
  )
}


if (opt$simulate_data_from_prior == TRUE) {
  assert_that(opt$simulate_data == TRUE)
}

assert_that(opt$error_dist %in% c("normal", "skew-normal", "student-t"))


if (opt$tc_file == TC_FILE_JSON) {
  tc_file <- opt$tc_file
} else {
  tc_file <- file.path(DIR_DATA, opt$tc_file)
}

if (opt$vl_first_pos_file == VL_FIRST_POS_FILE_TSV) {
  vl_first_pos_file <- opt$vl_first_pos_file
} else {
  vl_first_pos_file <- file.path(DIR_DATA, opt$vl_first_pos_file)
}

if (!is.null(opt$filter_by_testname)) {
  if (opt$filter_by_testname == "T1") {
    tc_file <- TC_FILE_JSON_N_GENE
    vl_first_pos_file <- VL_FIRST_POS_FILE_TSV_N_GENE
  }
}

# Set the number of threads per chain.
if (opt$threading == TRUE) {
  threads_per_chain <- 4
} else {
  threads_per_chain <- NULL
}


if (opt$exclude_too_steep_courses) {
  max_diff_load_perday <- opt$max_diff_load_perday
} else {
  max_diff_load_perday <- NULL
}

print(paste("Simulate viral load data:", opt$simulate_data))
print(paste("Remove problematic infections:", opt$remove_problematic_infections))

model_suffix_stan <- paste0(opt$model_stan, ifelse(opt$threading == TRUE, "_threading", ""))
opt_suffix_samples <- ifelse(opt$simulate_data == TRUE, paste0("_simvl", opt$n_samples), ifelse(!is.null(opt$n_samples), paste0("_samples", opt$n_samples), ""))
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
opt_adapt_delta <- ifelse(opt$adapt_delta == 0.8, "", paste0("_adapt_delta_", sub("\\.", "_", opt$adapt_delta)))
opt_use_paper1_data <- ifelse(opt$use_paper1_data, "_paper1_data", "")
opt_exclude_lrt_samples <- ifelse(opt$exclude_lrt_samples, "_exclude_lrt_samples", "")
opt_exclude_too_steep_courses <- ifelse(opt$exclude_too_steep_courses, "_exclude_too_steep_courses", "")
opt_max_diff_load_perday <- ifelse(opt$exclude_too_steep_courses, opt$max_diff_load_perday, "")
opt_min_duration <- ifelse(opt$min_duration, paste0("_min_", opt$min_duration, "_days_duration"), "")
opt_min_max_load <- paste0("_min_max_load_", opt$min_max_load)


# minimum number of data points per subjects
my_selections <- opt$min_tests:opt$min_tests

if (opt$no_vl_error) {
  opt$error_dist <- NULL
}

# run analyses for subjects with different
# minimum number of data points per subjects
for (selection in my_selections) {
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

  basic_dirname <- file.path(DIR_FIT, paste0("model", opt$model_stan))
  fn <- file.path(basic_dirname, paste0(model, ".Rdata"))

  if (!file.exists(fn)) {
    basic_filename <- file.path(basic_dirname, model)
    if (!dir.exists(basic_dirname)) {
      dir.create(basic_dirname)
    }
    print(paste0("File to be created: ", basic_filename, ".RDS"))

    # Note: we are setting force_recompile to true because we use a different 
    # cmdstan path when running the model via sbatch than via an OOD session.
    sm <- cmdstan_model(file.path(DIR_STAN, paste0("model", model_suffix_stan, ".stan")),
      cpp_options = list(stan_threads = opt$threading),
      force_recompile = TRUE
    )

    if (!opt$compile_only){
      make_time_course_standata(
        selection = selection, # minimum number of data points per infection
        max_n_tests = opt$max_tests, # maximum number of data points per infection
        min_n_pos_tests = opt$min_pos_tests, # Minimum number of positive PCRs per infection.
        imputation_limit = 3, # maximum value of imputed viral load for negative tests
        model_no_stan = opt$model_stan,
        samples = opt$n_samples,
        filter_by_testname = opt$filter_by_testname,
        tc_file = tc_file,
        vl_first_pos_file = vl_first_pos_file,
        simulate_data = opt$simulate_data,
        simulate_param_sizes_all_zero = opt$simulate_param_sizes_all_zero,
        simulate_mixed_courses = opt$simulate_mixed_courses,
        use_true_vls = opt$use_true_vls,
        trim_neg_pcrs = opt$trim_neg_pcrs,
        days_to_negative = opt$days_to_negative,
        single_leading_trailing_neg_pcr = opt$single_leading_trailing_neg_pcr,
        error_dist_sim_data = opt$error_dist,
        add_neg_pcrs = opt$add_neg_pcrs,
        no_shifts = opt$no_shifts,
        only_pos_pcrs_sim_data = opt$only_pos_pcrs,
        remove_problematic_infections = opt$remove_problematic_infections,
        threading = opt$threading,
        infection_wise_vl_error = !opt$overall_vl_error,
        use_paper1_data = opt$use_paper1_data,
        exclude_lrt_samples = opt$exclude_lrt_samples,
        min_duration = opt$min_duration,
        min_max_load = opt$min_max_load,
        max_diff_load_perday = max_diff_load_perday
      ) %>%
        list2env(.GlobalEnv)
  
      if (opt$no_shifts == TRUE) {
        # When using the actual sampling days, we can assume that there are
        # sampling days before the day of peak viral load, i.e. days < 0
        assertthat::assert_that(min(datalist$X_DAY) < 0)
      }
      if (opt$only_pos_pcrs == TRUE) {
        assertthat::assert_that(min(datalist$Y_DAY) > 2)
      }
  
      assertthat::assert_that(length(sim_data$intercept) == uniqueN(sim_data$dt$ID))
      assertthat::assert_that(length(sim_data$upSlope) == uniqueN(sim_data$dt$ID))
      assertthat::assert_that(length(sim_data$downSlope) == uniqueN(sim_data$dt$ID))
  
      assert_that(length(unique(names(datalist))) == length(names(datalist)))
  
      if (opt$process_rds_file) {
        csf <- readRDS(file = paste0(basic_filename, ".RDS"))
      } else {
        csf <- sm$sample(
          data = datalist,
          iter_warmup = opt$n_iter_warmup,
          iter_sampling = opt$n_iter_sampling,
          chains = opt$n_chains,
          parallel_chains = opt$n_chains,
          threads_per_chain = threads_per_chain,
          save_warmup = FALSE,
          adapt_delta = opt$adapt_delta,
          max_treedepth = opt$max_tree_depth,
          init = lapply(1:opt$n_chains, function(x) make_TC_inits(datalist)),
          seed = SEED
        )
  
        csf$save_object(file = paste0(basic_filename, ".RDS"))
      }
  
      draws <- csf$draws()
  
      ss <-
        draws %>%
        summarise_draws() %>%
        data.table()
      ss_mcse <-
        draws %>%
        summarise_draws(default_mcse_measures()) %>%
        data.table()
      sampler_diags <-
        csf$sampler_diagnostics() %>%
        as_draws_df() %>%
        data.table()
      inits <- csf$init()
      if (opt$simulate_data == TRUE) {
        save(csf, ss, ss_mcse, day_data, datalist, sampler_diags, inits, sim_data,
          file = paste0(basic_filename, ".Rdata")
        )
      } else {
        save(csf, ss, ss_mcse, day_data, datalist, sampler_diags, inits,
          file = paste0(basic_filename, ".Rdata")
        )
      }
  
      csf$save_output_files(dir = basic_dirname, basename = "output")
    }
  }
}

warnings()
