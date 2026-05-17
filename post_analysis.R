library(rlang)
library(data.table)
library(posterior)
library(assertthat)
library(here)


TOP <- here()
source(file.path(TOP, "data_params.R"))
source(file.path(TOP, "plot_params.R"))

MODELS_DAY2 <- c(2, 5)

DINA4_HEIGHT <- 8.35

VL_PLOT_WIDTH <- DINA4_WIDTH / 1.5 # for thesis: divide by 1 instead of 1.5
VL_PLOT_HEIGHT <- DINA4_HEIGHT / 3

VL_PLOT_WIDTH_THESIS <- DINA4_WIDTH / 1
VL_PLOT_HEIGHT_THESIS <- DINA4_HEIGHT / 3

VL_PLOT_WIDTH_PANEL_3 <- VL_PLOT_WIDTH
VL_PLOT_HEIGHT_PANEL_3 <- VL_PLOT_HEIGHT / 3


#' Load input data and posterior draws from an MCMC run
#' #'
#' @param model A string specifying the model name (includes numeric model
#' identifiers, data filtering criteria, and run settings).
#' @param model_dir Path to the directory containing model fitting results.
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param inc_warmup Samples from the warm-up phase are included in the results.
#' @param testing A logical value. If `TRUE`, loads only a subset of posterior
#' draws (for testing/debugging).
#' @param test_thin Thinning interval for posterior draws (e.g., `thin = 10`
#' selects every 10th draw).
#'
#' @return A named list with
#'  -`draws`: A `draws` object from the `posterior` package (containing values 
#'            sampled in each draw, for each model parameter).
#'  -`day_data`: A data.table with PCR test and metadata used for fitting.
#'  -`sim_data`: A named list, containing a data.table with simulated data used 
#'               for model fitting, and true (simulated) parameter values.
#'  -`ss`: A data.table containing summary statistics (including Rhat and ESS 
#'         values) from a specific run.
#'  -`datalist`: List with data for the stan model.
load_draws <- function(model,
                       model_dir,
                       model_no_stan,
                       inc_warmup,
                       testing = FALSE,
                       test_thin = 10) {
  sim_data <- NULL
  load(file.path(
    DIR_FIT, model_dir,
    paste0(model, ".Rdata")
  ))
  draws <- csf$draws(inc_warmup = inc_warmup)
  if (testing == TRUE) draws <- draws %>% thin_draws(test_thin)


  day_data <- label_problematic_infections(ss, day_data)
  return(list(
    draws = draws,
    day_data = day_data,
    sim_data = sim_data,
    ss = ss,
    datalist = datalist
  ))
}


##### Compute summary statistics #####

#' Compute or validate parameter differences.
#'
#' @param param_true Numeric vector. True parameter values (optional if `diff`
#' is provided).
#' @param param_est Numeric vector. Estimated parameter values (optional if
#' `diff` is provided).
#' @param diff Numeric vector. Precomputed differences (optional).
#' @return Numeric vector of differences (`param_est - param_true` if not
#' provided).
handle_diff <- function(param_true,
                        param_est,
                        diff) {
  # If diff is not provided, compute it
  if (is.null(diff)) {
    if (is.null(param_true) || is.null(param_est)) {
      stop("Either 'diff' or both 'param_true' and 'param_est' must be 
           provided.")
    }
    diff <- param_est - param_true
  }
  return(diff)
}


#' Compute the 94% HDI (lower and upper quartiles) on a numeric vector.
#'
#' @param x A numeric vector.
calc_94_hdi <- function(x) {
  return(unname(quantile(x, probs = c(0.03, 0.97))))
}


#' Compute the 50% HDI (lower and upper quartiles) on a numeric vector.
#'
#' @param x A numeric vector.
calc_50_hdi <- function(x) {
  return(unname(quantile(x, probs = c(0.25, 0.75))))
}


#' Compute a summary statistic on a numeric vector.
#'
#' @param summary_stat The summary statistic that will be computed.
#' @param x A numeric vector.
calc_summary_stat <- function(summary_stat, x) {
  if (is.null(summary_stat)) {
    return(x)
  } else {
    assert_that(summary_stat %in% c("mean", "median", "hdi", "hdi50"))
    if (summary_stat == "mean") {
      return(mean(x))
    } else if (summary_stat == "median") {
      return(median(x))
    } else if (summary_stat == "hdi") {
      return(calc_94_hdi(x))
    } else if (summary_stat == "hdi50") {
      return(calc_50_hdi(x))
    }
  }
}


#' Compute a specified summary statistic for a given model parameter,
#' aggregated across infections and stratified by posterior draw.
#'
#' @param dt A `data.table` containing estimated parameter values.
#' @param param_name The name of the model parameter (column name in `dt`)
#' for which the summary statistic will be computed.
#' @param summary_stat The summary statistic to compute. Supported options:
#' "mean", "median", "hdi", "hdi50".
#'
#' @return A vector where each element is the summary statistic for the
#' parameter, computed across all infections for a specific posterior draw.
calc_stat_by_draw <- function(dt,
                              param_name,
                              summary_stat = "median") {
  summary_by_draw <- dt[, .(summary = calc_summary_stat(summary_stat, 
                                                        get(param_name))), 
                        by = .draw]
  return(summary_by_draw$summary)
}


#' Compute the differences between estimated and true values, return a summary
#' statistic (mean, median or HPDI) across posterior draws for each difference.
#'
#' @param param_true A vector of true parameter values (known because data for
#' fitting were simulated), or a single value.
#' @param param_est A vector of estimated parameter values, or a single value.
#' @param summary_stat The summary statistic (mean, median, 94% HDI, 50% HDI)
#' that will be computed across posterior draws.
#'
#' @return Summary statistics of the difference between estimated and true
#' value. Alternatively, just a single such value.
calc_difference <- function(param_true,
                            param_est,
                            summary_stat = "median") {
  diff <- param_est - param_true
  return(calc_summary_stat(summary_stat, diff))
}


#' Compute the multiplicative deviations between estimated and true parameter
#' values. A summary statistic (e.g., mean, median) across posterior draws will
#' be computed and returned for each deviation.
#'
#' @param param_true A vector of true parameter values (known because data for
#' fitting were simulated), or a single value, or NULL.
#' @param param_est A vector of estimated parameter values, or a single value,
#' or NULL.
#' @param diff The difference between estimated and true values, or a single
#' difference value, or NULL.
#' @param exponentiate A logical value indicating whether the multiplicative
#' deviation of the exponentiated parameter values should be computed (i.e.
#' they are provided in log-space).
#' @param summary_stat The summary statistic (mean, median, 94% HDI, 50% HDI)
#' that will be computed across posterior draws.
#'
#' @return Summary statistics of the multiplicative deviations between
#' estimated and true parameter values. Alternatively, just a single such value.
calc_multiplicative_deviation <- function(param_true = NULL,
                                          param_est = NULL,
                                          diff = NULL,
                                          exponentiate = FALSE,
                                          summary_stat = "median") {
  if (exponentiate) {
    diff <- handle_diff(param_true, param_est, diff)
    x <- exp(diff) - 1
  } else {
    assert_that(!is.null(param_true) & !is.null(param_est))
    assert_that(param_true != 0)
    x <- (param_est - param_true) / abs(param_true)
  }
  return(calc_summary_stat(summary_stat, x))
}


#' Compute the proportional relative differences between estimated and true
#' parameter values, defined as (estimated - true) / true.
#' A summary statistic (e.g., mean, median) across posterior draws will be
#' computed and returned for each difference.
#' Multiply the results by 100 to convert to percentage deviations.

#' @param param_true A vector of true parameter values (known because data for
#' fitting were simulated), or a single value, or NULL.
#' @param param_est A vector of estimated parameter values, or a single value,
#' or NULL.
#' @param diff The difference between estimated and true values, or a single
#' difference value, or NULL.
#' @param exponentiate A logical value indicating whether the relative
#' difference of the exponentiated parameter values should be computed (i.e.
#' they are provided in log-space).
#' @param summary_stat The summary statistic (mean, median, 94% HDI, 50% HDI)
#' that will be computed across posterior draws.
#'
#' @return Summary statistics (e.g., means, medians) of the proportional
#' relative differences between the estimated and true parameter values.
#' Alternatively, just a single such value.
calc_relative_difference <- function(param_true = NULL,
                                     param_est = NULL,
                                     diff = NULL,
                                     exponentiate = FALSE,
                                     summary_stat = "median") {
  return(calc_multiplicative_deviation(
    param_true = param_true,
    param_est = param_est,
    diff = diff,
    exponentiate = exponentiate,
    summary_stat = summary_stat
  ))
}


#' Compute the difference between estimated and true values for posterior 
#' predictions of infection-group-specific quantities. A summary statistic 
#' (e.g., mean, median) across posterior draws will be computed and returned 
#' for each difference.
#' @param dt A data.table containing estimated and true values for each 
#' infection ID and draw.
#' @param param_name_true The name of the column containing the true parameter 
#' values (known because data for fitting were simulated), or NULL.
#' @param param_name_est The name of the column containing the estimated 
#' parameter values, or NULL.
#' @param diff The differences between estimated and true values, or NULL.
#' @param exponentiate A logical value indicating whether the relative 
#' difference of the exponentiated parameter values should be computed (i.e. 
#' they are provided in log-space).
#' @param summary_stat The summary statistic (mean, median, 94% HDI, 50% HDI) 
#' that will be computed across posterior draws.
#'
#' @return Summary statistics of the differences between estimated and true 
#' quantities for infection subgroups.
calc_difference_pop_avg <- function(dt,
                                    param_name_true,
                                    param_name_est,
                                    summary_stat = "median",
                                    summary_stat2 = "median") {
  # Compute the difference for each row using get()
  dt[, diff := get(param_name_est) - get(param_name_true)]

  # Apply the summary_stat2 by .draw
  result <- calc_stat_by_draw(dt, "diff", summary_stat = summary_stat2)

  # Apply the summary_stat to the results (e.g., to get HPDI of the medians)
  result_summary <- calc_summary_stat(summary_stat, result)

  return(result_summary)
}

#' Compute the multiplicative deviations between estimated and true values for 
#' posterior predictions of infection-group-specific quantities. A summary 
#' statistic (e.g., mean, median) across posterior draws will be computed and 
#' returned for each deviation.
#'
#' @param dt A data.table containing estimated and true values.
#' @param param_name_true The name of the column containing the true parameter 
#' values (known because data for fitting were simulated), or NULL.
#' @param param_name_est The name of the column containing the estimated 
#' parameter values, or NULL.
#' @param diff The differences between estimated and true values, or NULL.
#' @param exponentiate A logical value indicating whether the relative 
#' difference of the exponentiated parameter values should be computed (i.e. 
#' they are provided in log-space).
#' @param summary_stat The summary statistic (mean, median, 94% HDI, 50% HDI) 
#' that will be computed across posterior draws.
#' @param summary_stat2 The summary statistic (mean or median) that will be used 
#' for averaging over infections within subgroups first, to then summarize over 
#' draws.
#'
#' @return Summary statistics of the multiplicative deviations between 
#' estimated and true quantities for infection subgroups.
calc_multiplicative_deviation_pop_avg <- function(dt,
                                                  param_name_true = NULL,
                                                  param_name_est = NULL,
                                                  diff = NULL,
                                                  exponentiate = FALSE,
                                                  summary_stat = "median",
                                                  summary_stat2 = "median") {
  if (exponentiate) {
    diff <- handle_diff(dt[[param_name_true]], dt[[param_name_est]], diff)
    dt[, x := exp(diff)] # Store x as a new column in dt
  } else {
    assert_that(!is.null(param_name_true) & !is.null(param_name_est))
    # Store x as a new column in dt
    dt[, x := get(param_name_est) / get(param_name_true)]
  }
  # Apply summary_stat2 by .draw
  result <- calc_stat_by_draw(dt, "x", summary_stat = summary_stat2)

  # Apply summary_stat to the results (e.g., to get HPDI of the medians)
  result_summary <- calc_summary_stat(summary_stat, result)


  return(result_summary)
}


#' Compute the proportional relative differences between the estimated and true 
#' values, defined as (estimated - true) / true, for posterior predictions of 
#' infection-group-specific quantities.
#' Summary statistics (e.g., mean, median) across posterior draws of these 
#' proportional deviations will be computed and returned. Multiply the results 
#' by 100 to convert to percentage deviations.
#'
#' @param dt A data.table containing estimated and true values.
#' @param param_name_true The name of the column containing the true parameter 
#' values (known because data for fitting were simulated), or NULL.
#' @param param_name_est The name of the column containing the estimated parameter 
#' values, or NULL.
#' @param diff The differences between estimated and true values, or NULL.
#' @param exponentiate A logical value indicating whether the relative 
#' difference of the exponentiated parameter values should be computed (i.e. 
#' they are provided in log-space).
#' @param summary_stat The summary statistic (mean, median, 94% HDI, 50% HDI) 
#' that will be computed across posterior draws.
#' @param summary_stat2 The summary statistic (mean or median) that will be used 
#' for averaging over infections within subgroups first, to then summarize over 
#' draws.
#'
#' @return Summary statistics (e.g., mean, median) of the proportional relative 
#' difference between the estimated and true quantities for infection subgroups.
calc_relative_difference_pop_avg <- function(dt,
                                             param_name_true = NULL,
                                             param_name_est = NULL,
                                             diff = NULL,
                                             exponentiate = FALSE,
                                             summary_stat = "median",
                                             summary_stat2 = "median") {
  # Call calc_multiplicative_deviation with the same arguments
  result <- calc_multiplicative_deviation_pop_avg(
    dt,
    param_name_true,
    param_name_est,
    diff,
    exponentiate,
    summary_stat,
    summary_stat2
  )

  # Subtract 1 to get the relative multiplicative deviation
  return(result - 1)
}


##### Manipulate draws and calculate statistics ####


" Generate a draws data.table.
#"
#' @param draws A draws object from the `posterior` package.
#' @param thin An optional parameter to select subset of posterior draws.
#' @return A data.table.
as_draws_dt <- function(draws, thin = 1) {
  dt <-
    as_draws_df(draws) %>%
    thin_draws(thin) %>%
    data.table() %>%
    .[, c(".chain", ".iteration") := NULL]
  return(dt)
}


#' Generate list with basic stats (mean 3% and 97% quantiles). 
#' Mainly for use in data.tables.
#'
#' @param x A vector with real numbers.
#' @return A list.
post_stats_list <- function(x) {
  m <- collapse::fmean(x)
  qs <- quantile(x, c(.03, .97), names = FALSE)
  return(list(m = m, q3 = qs[1], q97 = qs[2]))
}


#' Generate list with mean, median and a larger number of quantiles. 
#' Mainly for use in data.tables.
#'
#' @param x A vector with real numbers.
#' @param quantiles A vector with quantiles to be calculated 
#' (default : seq(.5,.95, by = .05)).
#' @return A list.
my_stats_list_long <- function(x, quantiles = seq(.5, .90, by = .05)) {
  m <- collapse::fmean(x)
  med <- collapse::fmedian(x)
  qs <- quantile(x, c((1 - quantiles) / 2, 1 - (1 - quantiles) / 2), 
                 names = FALSE)
  dstats <- as.list(c(m = m, med = med, qs))
  names(dstats) <- c("mean", "median", c(paste0("lower", quantiles * 100), 
                                         paste0("upper", quantiles * 100)))
  return(dstats)
}


#' Calculate posterior statistics for a number of variables.
#'
#' @param dt A long-form data.table with posterior samples and auxiliary 
#' variables.
#' @param var A variable in `dt` for which statistics are calculated.
#' @param by A grouping variable in `dt`.
#' @param quantiles A vector with quantiles to be calculated 
#' (default : seq(.5,.95, by = .05)).
#' @return A data.table.
get_stats <- function(dt, 
                      var = "value", 
                      by = NULL, 
                      quantiles = seq(.5, .95, by = .05)) {
  my_stats <-
    dt[, as.list(my_stats_list_long(get(var), quantiles = quantiles)), by = by]
  return(my_stats)
}


#' Compute log(1 + exp(x)) numerically stably.
#'
#' Avoids overflow/underflow by using `log1p` and conditional logic.
#' For x > 0, computes x + log1p(exp(-x)); otherwise, computes log1p(exp(x)).
#'
#' @param x Numeric vector. Input values.
#' @return Numeric vector. log(1 + exp(x)) computed stably for all x.
log1p_exp <- function(x) {
  ifelse(x > 0,
    x + log1p(exp(-x)), # Stable for large positive x
    log1p(exp(x))
  ) # Stable for other values of x
}


#' Determine the alpha variable name in a Stan/MCMC draws object.
#'
#' Checks for common alpha parameter naming conventions 
#' ("alpha", "alpha_mu", "alpha[1]") and returns the first match found. 
#' Returns NULL if no alpha parameter is present.
#'
#' @param draws Object. Stan/MCMC draws or similar object with `variables()`
#'  method.
#' @return Character or NULL. Name of the alpha parameter if found, otherwise 
#' NULL.
determine_alpha_var <- function(draws) {
  if ("alpha" %in% variables(draws)) {
    alpha_var <- "alpha"
  } else if ("alpha_mu" %in% variables(draws)) {
    alpha_var <- "alpha_mu"
  } else if ("alpha[1]" %in% variables(draws)) {
    alpha_var <- "alpha[]"
  } else {
    alpha_var <- NULL
  }
  return(alpha_var)
}


#' Extract posterior draws and group by infection identifiers (id).
#'
#' @param draws A draws object from the `posterior` package.
#' @param day_data A data.table with PCR test data.
#' @param params A parameter to be extracted from the draws object.
#' @param thin An optional parameter to select subset of posterior draws.
#' @return A data.table.
draws_by_id <- function(draws, day_data, params, thin = 1) {
  for (p in params) {
    tmp <-
      subset_draws(draws, p) %>%
      as_draws_dt(thin = thin) %>%
      melt(id.vars = ".draw", variable.name = "ID") %>%
      setnames("value", p)

    if (length(params) == 1 | p == params[1]) {
      draws.dt <- copy(tmp)
    } else {
      draws.dt <- cbind(draws.dt, tmp[, c(p), with = FALSE])
    }
  }
  draws.dt %>%
    .[, ID_stan := as.numeric(gsub("[^0-9]", "", ID, perl = TRUE))] %>%
    .[, ID := factor(ID_stan, labels = day_data[day == 0, ID])] %>%
    setkeyv(c(".draw", "ID"))
  return(draws.dt)
}


#' Calculate statics over posterior draws that were first aggregated.
#'
#' @param dt A data.table with an indicator variable `.draw`.
#' @param by A grouping variable.
#' @param target.var A parameter for which statistics are calculated.
#' @return A data.table.
summarise_draws_dt_by <- function(dt, by, target.var = "value", varname = NULL,
                                  stat = "median") {
  agg_fun <- match.fun(if (stat == "median") collapse::fmedian else collapse::fmean)
  summarised_draws <-
    dt %>%
    .[, list(m = agg_fun(get(target.var))), by = c(".draw", by)] %>%
    .[, as.list(post_stats_list(m)), by = by]
  if (!is.null(varname)) {
    setnames(summarised_draws, "m", varname)
  }
  return(summarised_draws)
}


#' Compute infection-group-specific posterior estimates of the main model 
#' parameters (intercept, up-slope etc.).
#'
#' @param draws A `draws` object from the `posterior` package (containing values 
#' sampled in each draw, for each model parameter).
#' @param grp_var Character. The variable to group by.
#' @param thin Thinning interval for posterior draws 
#' (e.g., `thin = 10` selects every 10th draw).
#' @param calc_times A logical value indicating whether to compute estimated 
#' proliferation and clearance times.
#' @param calc_bc A logical value indicating whether to compute smoothing 
#' parameters (`b` and `c`).
#' @param reference_variant Character. Specifies the reference SARS-CoV-2 
#' variant.
#' @param reference_hosp_status Character. Specifies the reference disease 
#' severity status.
#' @return A data.table with infection-group-specific estimates of the main 
#' model parameters per draw.
draws_by_grp <- function(draws,
                         grp_var,
                         thin = 1,
                         calc_times = FALSE,
                         calc_bc = TRUE,
                         reference_variant = "omicron",
                         reference_hosp_status = "non-PAMS, non-hosp.") {
  assert_that(reference_variant %in% c("wildtype", "omicron"))
  assert_that(reference_hosp_status %in% c(
    "hospitalized",
    "non-PAMS, non-hosp.",
    "PAMS"
  ))
  # The following has only been implemented for a 2-slopes model using the 
  # alpha var. Make sure we are not dealing with a three-slope model
  assert_that(!"beta2[1]" %in% variables(draws))
  # Make sure we are not dealing with a two-slope model using beta_sweight_mu
  assert_that(!"beta_sweight_mu" %in% variables(draws))

  mapping <- SUBGRP_VARS_MAPPING
  idx_mapping <- SUBGRPS
  label_mapping <- SUBGRP_CATS_LABEL_MAPPING2

  assert_that(grp_var %in% names(mapping))

  grp_cats <- mapping[[grp_var]]
  grp_idcs <- which(idx_mapping %in% grp_cats)
  # We are using Omicron infections as the reference category
  omi_variant_idx <- which(idx_mapping == "omicron")
  hosp_idx <- which(idx_mapping == "hospitalized")
  pams_idx <- which(idx_mapping == "PAMS")

  # If hospitalized is used as the reference category in PAMS3
  intercept_hosp <- paste0("betaPGH_intercept[", hosp_idx, "]")
  slope_up_hosp <- paste0("betaPGH_slope_up[", hosp_idx, "]")
  slope_down_hosp <- paste0("betaPGH_slope_down[", hosp_idx, "]")

  # If PAMS is used as the reference category in PAMS3
  intercept_pams <- paste0("betaPGH_intercept[", pams_idx, "]")
  slope_up_pams <- paste0("betaPGH_slope_up[", pams_idx, "]")
  slope_down_pams <- paste0("betaPGH_slope_down[", pams_idx, "]")

  intercept_omi <- paste0("betaPGH_intercept[", omi_variant_idx, "]")
  slope_up_omi <- paste0("betaPGH_slope_up[", omi_variant_idx, "]")
  slope_down_omi <- paste0("betaPGH_slope_down[", omi_variant_idx, "]")


  return_draws_dt <- function(draws, grp_idcs, thin) {
    grp_params_intercept <- paste0("betaPGH_intercept[", grp_idcs, "]")
    grp_params_slope_up <- paste0("betaPGH_slope_up[", grp_idcs, "]")
    grp_params_slope_down <- paste0("betaPGH_slope_down[", grp_idcs, "]")

    params <- c(
      "log_slope_up_mu", "log_slope_down_mu", "log_intercept_mu",
      grp_params_intercept, grp_params_slope_up,
      grp_params_slope_down
    )

    # Reference category draws
    draws_dt <-
      subset_draws(draws, params) %>%
      as_draws_dt(thin = thin)
    return(draws_dt)
  }

  alpha_hosp_idx <- which(idx_mapping == "hospitalized:alpha")
  delta_hosp_idx <- which(idx_mapping == "hospitalized:delta")
  omicron_hosp_idx <- which(idx_mapping == "hospitalized:omicron")
  grp_idcs <- c(
    grp_idcs, omi_variant_idx, hosp_idx, pams_idx,
    alpha_hosp_idx, delta_hosp_idx, omicron_hosp_idx
  )

  draws_dt <- return_draws_dt(draws = draws, grp_idcs = grp_idcs, thin = thin)

  intercept_alpha_hosp <- paste0("betaPGH_intercept[", alpha_hosp_idx, "]")
  slope_up_alpha_hosp <- paste0("betaPGH_slope_up[", alpha_hosp_idx, "]")
  slope_down_alpha_hosp <- paste0("betaPGH_slope_down[", alpha_hosp_idx, "]")

  intercept_delta_hosp <- paste0("betaPGH_intercept[", delta_hosp_idx, "]")
  slope_up_delta_hosp <- paste0("betaPGH_slope_up[", delta_hosp_idx, "]")
  slope_down_delta_hosp <- paste0("betaPGH_slope_down[", delta_hosp_idx, "]")

  intercept_omi_hosp <- paste0("betaPGH_intercept[", omicron_hosp_idx, "]")
  slope_up_omi_hosp <- paste0("betaPGH_slope_up[", omicron_hosp_idx, "]")
  slope_down_omi_hosp <- paste0("betaPGH_slope_down[", omicron_hosp_idx, "]")

  if ((grp_var == "PAMS3") ||
    (grp_var == "variant") ||
    (reference_hosp_status %in% c("PAMS", "non-PAMS, non-hosp.")) ||
    (reference_variant == "wildtype")) {
    intercept_add_omi_hosp <- 0
    slope_up_add_omi_hosp <- 0
    slope_down_add_omi_hosp <- 0
  } else {
    assert_that(reference_hosp_status == "hospitalized")
    assert_that(reference_variant == "omicron")

    intercept_add_omi_hosp <- draws_dt[[intercept_omi_hosp]]
    slope_up_add_omi_hosp <- draws_dt[[slope_up_omi_hosp]]
    slope_down_add_omi_hosp <- draws_dt[[slope_down_omi_hosp]]
  }

  if ((grp_var == "PAMS3") || 
      (reference_hosp_status == "non-PAMS, non-hosp.")) {
    intercept_add_pams <- 0
    slope_up_add_pams <- 0
    slope_down_add_pams <- 0

    intercept_add_hosp <- 0
    slope_up_add_hosp <- 0
    slope_down_add_hosp <- 0
  } else if (reference_hosp_status == "hospitalized") {
    intercept_add_pams <- 0
    slope_up_add_pams <- 0
    slope_down_add_pams <- 0

    intercept_add_hosp <- draws_dt[[intercept_hosp]]
    slope_up_add_hosp <- draws_dt[[slope_up_hosp]]
    slope_down_add_hosp <- draws_dt[[slope_down_hosp]]
  } else if (reference_hosp_status == "PAMS") {
    intercept_add_pams <- draws_dt[[intercept_pams]]
    slope_up_add_pams <- draws_dt[[slope_up_pams]]
    slope_down_add_pams <- draws_dt[[slope_down_pams]]

    intercept_add_hosp <- 0
    slope_up_add_hosp <- 0
    slope_down_add_hosp <- 0
  }

  # We are using Omicron infections as the reference category for all
  # variables other than variant
  if ((grp_var == "variant") || (reference_variant == "wildtype")) {
    intercept_add_omi <- 0
    slope_up_add_omi <- 0
    slope_down_add_omi <- 0
  } else {
    assert_that(reference_variant == "omicron")
    intercept_add_omi <- draws_dt[[intercept_omi]]
    slope_up_add_omi <- draws_dt[[slope_up_omi]]
    slope_down_add_omi <- draws_dt[[slope_down_omi]]
  }

  draws_dt_new <- data.table(
    intercept = exp(draws_dt[["log_intercept_mu"]] + intercept_add_omi + intercept_add_hosp + intercept_add_pams + intercept_add_omi_hosp),
    slope_up = exp(draws_dt[["log_slope_up_mu"]] + slope_up_add_omi + slope_up_add_hosp + slope_up_add_pams + slope_up_add_omi_hosp),
    slope_down = -exp(draws_dt[["log_slope_down_mu"]] + slope_down_add_omi + slope_down_add_hosp + slope_down_add_pams + slope_down_add_omi_hosp),
    .draw = draws_dt[[".draw"]]
  )

  draws_dt_new[, (grp_var) := "reference"]

  iterate_grp_cats <- function(grp_cats, draws_dt, draws_dt_new, idx_mapping) {
    for (grp_cat in grp_cats) {
      # Note: the following scenario can only happen if grp_var is equal to
      # "variant" or to "PAMS3", so add_omi_hosp variables will be set to 0
      # Add interaction effect between hospitalization and variant
      if (((grp_cat == "hospitalized") && (reference_variant == "omicron")) ||
        ((grp_cat == "omicron") && (reference_hosp_status == "hospitalized"))) {
        add_to_intercept <- draws_dt[[intercept_omi_hosp]]
        add_to_slope_up <- draws_dt[[slope_up_omi_hosp]]
        add_to_slope_down <- draws_dt[[slope_down_omi_hosp]]
        # We only allow wildtype and omicron as reference categories for variant
      } else if ((grp_cat == "delta") && (reference_hosp_status == "hospitalized")) {
        add_to_intercept <- draws_dt[[intercept_delta_hosp]]
        add_to_slope_up <- draws_dt[[slope_up_delta_hosp]]
        add_to_slope_down <- draws_dt[[slope_down_delta_hosp]]
      } else if ((grp_cat == "alpha") && (reference_hosp_status == "hospitalized")) {
        add_to_intercept <- draws_dt[[intercept_alpha_hosp]]
        add_to_slope_up <- draws_dt[[slope_up_alpha_hosp]]
        add_to_slope_down <- draws_dt[[slope_down_alpha_hosp]]
      } else {
        add_to_intercept <- 0
        add_to_slope_up <- 0
        add_to_slope_down <- 0
      }

      dt_tmp <- data.table()
      cat_idx <- which(idx_mapping %in% grp_cat)
      # intercept_param = paste0("intercept_", grp_cat)
      # slope_up_param = paste0("slope_up_", grp_cat)
      # slope_down_param = paste0("slope_down_", grp_cat)

      intercept_grp <- paste0("betaPGH_intercept[", cat_idx, "]")
      slope_up_grp <- paste0("betaPGH_slope_up[", cat_idx, "]")
      slope_down_grp <- paste0("betaPGH_slope_down[", cat_idx, "]")

      dt_tmp[, intercept := exp(draws_dt[["log_intercept_mu"]] +
        draws_dt[[intercept_grp]] +
        add_to_intercept +
        intercept_add_omi +
        intercept_add_pams +
        intercept_add_hosp +
        intercept_add_omi_hosp)]
      dt_tmp[, slope_up := exp(draws_dt[["log_slope_up_mu"]] +
        draws_dt[[slope_up_grp]] +
        add_to_slope_up +
        slope_up_add_omi +
        slope_up_add_pams +
        slope_up_add_hosp +
        slope_up_add_omi_hosp)]
      dt_tmp[, slope_down := -exp(draws_dt[["log_slope_down_mu"]] +
        draws_dt[[slope_down_grp]] +
        add_to_slope_down +
        slope_down_add_omi +
        slope_down_add_pams +
        slope_down_add_hosp +
        slope_down_add_omi_hosp)]
      dt_tmp[, (grp_var) := grp_cat]
      dt_tmp[, .draw := draws_dt[[".draw"]]]
      draws_dt_new <- rbind(draws_dt_new, dt_tmp)
    }

    return(draws_dt_new)
  }

  draws_dt <- iterate_grp_cats(
    grp_cats = grp_cats, draws_dt = draws_dt,
    draws_dt_new = draws_dt_new,
    idx_mapping = idx_mapping
  )
  if (calc_bc) {
    draws_dt[, b := (slope_up + slope_down) / 2]
    draws_dt[, c := (slope_down - slope_up) / 2]
  }
  if (calc_times) {
    draws_dt[, time2peak := intercept / slope_up]
    draws_dt[, time_from_peak := -intercept / slope_down]
    draws_dt[, time_from_peak2 := -(intercept - 3) / slope_down]
  }
  # Rename columns
  labels <- unlist(unname(label_mapping[grp_var]))
  val_name_mapping <- setNames(as.list(labels[2:length(labels)]), grp_cats)
  val_name_mapping["reference"] <- labels[1]
  draws_dt[, (grp_var) := factor(val_name_mapping[get(grp_var)],
    levels = labels,
    labels = labels
  )]

  return(draws_dt)
}


#' Generate predicted time courses of viral load and culture positivity for
#' specific infection subgroups (predictions are computed from population-level
#' parameter estimates).
#'
#' @param grp_var Character. Column name to group by.
#' @param draws A draws object from the time course model generated with the 
#' posterior package.
#' @param days Vector of days (can be fractions of days) for which to calculate 
#' predictions.
#' @param thin Thinning interval for posterior draws 
#' (e.g., `thin = 10` selects every 10th draw).
#' @return A data.table with time courses of viral load and culture positivity.
make_VLCP_by_grp <- function(grp_var,
                             draws,
                             days = seq(-10, 40, by = 1),
                             thin = 1,
                             reference_variant = "omicron",
                             reference_hosp_status = "non-PAMS, non-hosp.") {
  assert_that(reference_variant %in% c("wildtype", "omicron"))
  assert_that(reference_hosp_status %in% c(
    "PAMS", "hospitalized",
    "non-PAMS, non-hosp."
  ))

  draws_dt <- draws_by_grp(
    draws = draws, grp_var = grp_var, thin = thin,
    reference_variant = reference_variant,
    reference_hosp_status = reference_hosp_status
  )

  alpha_var <- determine_alpha_var(draws)
  alpha_draw <-
    subset_draws(draws, c(alpha_var)) %>%
    as_draws_dt(thin = thin)

  draws_dt <- merge(draws_dt,
    alpha_draw,
    by = c(".draw"),
    allow.cartesian = TRUE
  )

  # 2-slope model
  VLCP_by_draw_grp <-
    draws_dt[,
      as.list(yhat_smooth_2_slopes(days, intercept, b, c, get(alpha_var))),
      by = c(grp_var, ".draw")
    ] %>%
    melt(id.var = c(grp_var, ".draw"), value.name = "log10_load") %>%
    .[, day_shifted := days[variable]] %>%
    .[, variable := NULL]

  setkeyv(VLCP_by_draw_grp, ".draw")

  if ("alpha_CP" %in% variables(draws)) {
    CP_params <- subset_draws(draws, c("alpha_CP", "beta_CP")) %>%
      as_draws_dt(thin = thin) %>%
      setnames("beta_CP[1]", "beta_CP")
    setkeyv(CP_params, ".draw")
    VLCP_by_draw_grp <-
      VLCP_by_draw_grp[CP_params, 
                       CP := inv.logit(alpha_CP + log10_load * beta_CP)]
  }
  return(VLCP_by_draw_grp)
}


#' Generate predicted time courses of viral load and culture positivity for 
#' each individual and posterior sample.
#'
#' @param draws A draws object from the time course model generated with the 
#' `posterior` package.
#' @param days A vector of days (can be fractions of days) for which to 
#' calculate predictions.
#' @param thin An optional parameter to select subset of posterior draws.
#' @return aA data.table with time courses of viral load and culture positivity.
make_VLCP_by_draw_ID <- function(draws, 
                                 day_data, 
                                 days = seq(-10, 40, by = 1), 
                                 thin = 1) {
  alpha_var <- determine_alpha_var(draws)

  # If the sigmoid function is used for smoothing
  if (is.null(alpha_var)) {
    by_draw_ID <-
      draws_by_id(draws,
        day_data,
        c("slope_up", "slope_down", "intercept"),
        thin = thin
      )

    beta_sweight_draw <-
      subset_draws(draws, c("beta_sweight_mu")) %>%
      as_draws_dt(thin = thin)

    by_draw_ID <- merge(by_draw_ID,
      beta_sweight_draw,
      by = c(".draw"),
      allow.cartesian = TRUE
    )

    VLCP_by_draw_ID <-
      by_draw_ID[,
        as.list(yhat_2_slopes(days, 
                              intercept, 
                              slope_up, 
                              slope_down, 
                              beta_sweight_mu)),
        by = c("ID", ".draw")
      ] %>%
      melt(id.var = c("ID", ".draw"), value.name = "log10_load") %>%
      .[, variable := as.numeric(gsub("[^0-9]", "", variable))] %>%
      .[, day_shifted := days[variable]] %>%
      .[, variable := NULL]
  } else {
    # If instead the transition at the peak was smoothed according to 
    # https://springerplus.springeropen.com/articles/10.1186/s40064-016-3278-y#Equ1
    # 3-slope model
    three_slopes <- "beta2[1]" %in% variables(draws)
    if (alpha_var == "alpha[]") {
      if (three_slopes) {
        params <- c("a", "b", "c1", "c2", "beta2", "alpha")
      } else {
        params <- c("intercept", "b", "c", "alpha")
      }
      by_draw_ID <- draws_by_id(draws, day_data, params, thin = thin)
      alpha_var <- "alpha"
    } else {
      if (three_slopes) {
        params <- c("a", "b", "c1", "c2", "beta2")
      } else {
        params <- c("intercept", "b", "c")
      }
      by_draw_ID <- draws_by_id(draws, day_data, params, thin = thin)
      alpha_draw <- subset_draws(draws, c(alpha_var)) %>% as_draws_dt(thin = thin)
      by_draw_ID <- merge(by_draw_ID,
        alpha_draw,
        by = c(".draw"),
        allow.cartesian = TRUE
      )
    }
    if ("beta2[1]" %in% variables(draws)) {
      VLCP_by_draw_ID <-
        by_draw_ID[, as.list(yhat_smooth_3_slopes(days, 
                                                  a, 
                                                  b, 
                                                  c1, 
                                                  c2, 
                                                  beta2, 
                                                  get(alpha_var))),
          by = c("ID", ".draw")
        ] %>%
        melt(id.var = c("ID", ".draw"), value.name = "log10_load") %>%
        .[, variable := as.numeric(gsub("[^0-9]", "", variable))] %>%
        .[, day_shifted := days[variable]] %>%
        .[, variable := NULL]
    } else {
      # 2-slope model
      VLCP_by_draw_ID <-
        by_draw_ID[,
          as.list(yhat_smooth_2_slopes(days, intercept, b, c, get(alpha_var))),
          by = c("ID", ".draw")
        ] %>%
        melt(id.var = c("ID", ".draw"), value.name = "log10_load") %>%
        .[, variable := as.numeric(gsub("[^0-9]", "", variable))] %>%
        .[, day_shifted := days[variable]] %>%
        .[, variable := NULL]
    }
  }
  setkeyv(VLCP_by_draw_ID, ".draw")

  if ("alpha_CP" %in% variables(draws)) {
    CP_params <- subset_draws(draws, c("alpha_CP", "beta_CP")) %>%
      as_draws_dt(thin = thin) %>%
      setnames("beta_CP[1]", "beta_CP")
    setkeyv(CP_params, ".draw")
    VLCP_by_draw_ID <-
      VLCP_by_draw_ID[CP_params, CP := inv.logit(alpha_CP + log10_load * beta_CP)]
  }

  return(VLCP_by_draw_ID)
}


#' Compute parameters for a smooth two-slope model.
#'
#' Transforms intercept and two slopes (`slope_up`, `slope_down`) into parameters
#' for a smooth transition model: `a` (intercept), `b` (average slope), and 
#' `c` (slope difference/2).
#'
#' @param intercept Numeric. Intercept term.
#' @param slope_up Numeric. Slope for the "up" phase.
#' @param slope_down Numeric. Slope for the "down" phase.
#' @return Named list. Parameters `a`, `b`, and `c` for the smooth transition 
#' model.
compute_smooth_2_slope_params <- function(intercept, slope_up, slope_down) {
  a <- intercept
  b <- (slope_up + slope_down) / 2
  c <- (slope_down - slope_up) / 2

  return(list(a = a, b = b, c = c))
}


#' Compute parameters for a smooth three-slope model.
#'
#' Transforms intercept, two slopes (`slope_up`, `slope_down`), and a lower
#'  bound `L` into parameters for a smooth three-phase transition model:
#' `a` (adjusted intercept), 
#' `b` (half of `slope_up`), 
#' `c1` (half slope difference),
#' `c2` (half of `slope_down`), and 
#' `beta2` (second breakpoint).
#'
#' @param intercept Numeric. Intercept term.
#' @param slope_up Numeric. Slope for the initial "up" phase.
#' @param slope_down Numeric. Slope for the final "down" phase.
#' @param L Numeric. Lower bound for the infection range.
#' @return Named list. Parameters `a`, `b`, `c1`, `c2`, and `beta2` for the 
#' three-slope model.
compute_smooth_3_slope_params <- function(intercept, slope_up, slope_down, L) {
  range_infection <- intercept - L

  b <- slope_up / 2 # Final slope is 0
  c1 <- (slope_down - slope_up) / 2
  c2 <- -slope_down / 2
  beta2 <- -(range_infection / slope_down) # Second breakpoint
  a <- intercept - c2 * beta2

  return(list(a = a, b = b, c1 = c1, c2 = c2, beta2 = beta2))
}


#' Compute smooth two-slope prediction.
#'
#' Models a response as a smooth combination of two linear slopes, 
#' with transition smoothness controlled by `alpha`. Uses `log1p_exp` for 
#' numerical stability.
#'
#' @param days Numeric vector. Days relative to the day of peak viral load.
#' @param a Numeric. Intercept term.
#' @param b Numeric. Average slope.
#' @param c Numeric. Half the difference between slopes.
#' @param alpha Numeric. Transition rate for logistic weighting.
#' @return Numeric vector. Smooth prediction for each input `days`.
yhat_smooth_2_slopes <- function(days, a, b, c, alpha) {
  return((a + b * days + c * days + 2 / alpha * c * log1p_exp(-alpha * days)))
}


#' Compute smooth three-slope prediction.
#'
#' Models a response as a smooth combination of three linear phases, with 
#' transitions controlled by `alpha` and a breakpoint at `beta2`. 
#' Uses `log1p_exp` for numerical stability.
#'
#' @param days Numeric vector. Days relative to the day of peak viral load.
#' @param a Numeric. Adjusted intercept.
#' @param b Numeric. Half of the initial slope.
#' @param c1 Numeric. Half the difference between initial and middle slopes.
#' @param c2 Numeric. Half of the final slope.
#' @param beta2 Numeric. Second breakpoint for the transition.
#' @param alpha Numeric. Transition rate for logistic weighting.
#' @return Numeric vector. Smooth prediction for each input `days`.
yhat_smooth_3_slopes <- function(days, a, b, c1, c2, beta2, alpha) {
  return(a + b * days + c1 * days + 2 / alpha * c1 * log1p_exp(-alpha * days)
    + c2 * (days - beta2) + 2 / alpha * c2 * log1p_exp(-alpha * (days - beta2)))
}


#' Calculate viral load trajectory using a two- or three-slope model.
#'
#' Implements the smooth transition model described in
#' <https://springerplus.springeropen.com/articles/10.1186/s40064-016-3278-y#Equ1>.
#' For `n_slopes=2`, uses `yhat_smooth_2_slopes`; for `n_slopes=3`, 
#' uses `yhat_smooth_3_slopes`.
#'
#' @param days Numeric vector. Days relative to the day of peak viral load.
#' @param intercept Numeric. Intercept term.
#' @param slope_up Numeric. Initial slope.
#' @param slope_down Numeric. Final slope.
#' @param alpha Numeric. Transition rate.
#' @param n_slopes Integer. Number of slopes (2 or 3). Defaults to 3.
#' @return Numeric vector. Predicted viral load trajectory for each input 
#' `days`.
smooth_function <- function(days, 
                            intercept, 
                            slope_up, 
                            slope_down, 
                            alpha, 
                            n_slopes = 3) {
  assertthat::assert_that(n_slopes %in% c(2, 3))
  if (n_slopes == 2) {
    b <- (slope_up + slope_down) / 2
    c <- (slope_down - slope_up) / 2
    return(yhat_smooth_2_slopes(days, 
                                a = intercept, 
                                b = b, 
                                c = c, 
                                alpha = alpha))
  } else {
    range_infection <- intercept - 3
    b <- slope_up / 2
    c1 <- (slope_down - slope_up) / 2
    c2 <- -slope_down / 2
    beta2 <- -(range_infection / slope_down) # Second breakpoint
    a <- intercept - c2 * beta2
    return(yhat_smooth_3_slopes(days, a, b, c1, c2, beta2, alpha))
  }
}


#' Assign Numeric IDs to unique infection IDs
#'
#' Adds a new column `ID_numeric` to `day_data`, assigning a unique numeric 
#' group identifier to each unique value in the `ID` column. Useful for 
#' downstream analysis or plotting.
#'
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @return A data.table containing input data with an added column 
#' `ID_numeric` (integer).
map_IDs <- function(day_data) {
  day_data <- day_data %>%
    .[, ID_numeric := .GRP, by = .(ID)]
}


#' Label infections with poor convergence (R-hat ≥ 1.05) of the temporal shift 
#' parameter.
#'
#' @param ss A data.table. Posterior summary statistics (e.g., from `rstan` or 
#' `brms`), including `variable` and `rhat`.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @return data.table. Input `day_data` with an added logical column 
#' `is.difficult` (`TRUE` for infections with `rhat ≥ 1.05` for temporal shift 
#' parameter).
label_problematic_infections <- function(ss, day_data) {
  is.difficult.infections <-
    ss %>%
    .[grepl("^shift\\[", variable) & rhat >= 1.05, ] %>%
    .[, tmp_id := as.numeric(gsub("[^0-9.]", "", variable))] %>%
    .[, ID := day_data[day == 0]$ID[tmp_id]]

  day_data <- day_data %>%
    .[, is.difficult := FALSE] %>%
    .[ID %in% unique(is.difficult.infections$ID), is.difficult := TRUE]

  return(day_data)
}


#' Adjust Skew-Normal posterior draws for population-level correction
#'
#' @param draws A `draws` object from the `posterior` package 
#' (containing values sampled in each draw, for each model parameter).
#' @return A `draws_array` object with adjusted `log_intercept_mu` to account 
#' for skew-normal mean correction.
process_skew_normal_draws <- function(draws) {
  draws_dt <- as_draws_dt(draws)
  draws_dt[, intercept := intercept_mean]
  # Compute the "correction" term that is added to go from location to expected value
  # (mean) of the skew-normal distribution
  # This is without random effects (omega_mu as opposed to omega) as we want to
  # add it to log_intercept_mu
  draws_dt[, correction_pop_level := exp(omega_mu) * delta * sqrt(2 / pi)]
  # We are adding the correction on the log scale so we don't have to take care of that later on
  draws_dt[, log_intercept_mu := log_intercept_mu + log1p(correction_pop_level * exp(-log_intercept_mu))]
  draws <- as_draws_array(draws_dt)
  return(draws)
}


#' Compute summary statistics for temporal shifts.
#'
#' @param shifted_data_by_draw_day_ID A data.table with posterior draws of 
#' shifted days (i.e. estimated days of sampling relative to the day of peak 
#' viral load), and imputed log10 viral loads grouped by infection.
#' @param day_var A string specifying the day variable used in the model 
#' (e.g. "day" with day 0 being the first day with a log10 viral load > 3; 
#' or "day2" with day 0 being the day of the highest viral load measured in the 
#' infection).
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed on shift parameters.
#' @param color_by A column name to color data points by (e.g. imputed viral 
#' loads).
#' @return List with:
#'   - `tmp_shifted_data`: Aggregated data.table.
#'   - `alpha`: Numeric transparency value for plotting (value depends on the 
#'   number of data points).
compute_shifted_data <- function(shifted_data_by_draw_day_ID,
                                 day_var,
                                 stat,
                                 color_by) {
  compute_summary_stat <- TRUE
  agg_fun <- match.fun(if (stat == "median") median else mean)

  if (compute_summary_stat == TRUE) {
    if (day_var == "day") {
      tmp_shifted_data <-
        shifted_data_by_draw_day_ID %>%
        .[, .(
          day_shifted = agg_fun(day_shifted),
          log10_load = agg_fun(log10_load)
        ), by = .(ID, ID_stan, day, color_column = get(color_by))]
      alpha <- 0.1 # ifelse(uniqueN(day_data$infection_hash) < 4000, 0.2, 0.05)
    } else {
      tmp_shifted_data <-
        shifted_data_by_draw_day_ID %>%
        .[, .(
          day_shifted = agg_fun(day_shifted),
          log10_load = agg_fun(log10_load)
        ), by = .(ID, ID_stan, value = get(day_var), day, 
                  color_column = get(color_by))] %>%
        rename(!!day_var := value)
      alpha <- 0.1
    }
  } else {
    alpha <- .01
    tmp_shifted_data <- shifted_data_by_draw_day_ID
  }

  tmp <- gc(verbose = FALSE)

  return(list(tmp_shifted_data = tmp_shifted_data, alpha = alpha))
}

#' Compute the summary statistic of the estimated (or true) values of a specific
#'  parameter for an infection subgroup.
#'
#' @param dt A data.table of posterior draws, grouped by infection identifiers 
#' (ID).
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param param The parameter of interest (e.g. intercept).
#' @param category_ids The infection identifiers that belong to a specific 
#' infection subgroup.
#' @param param_suffix A string specifying the suffix of the parameter of 
#' interest (e.g. "_true" if we are summarizing over true parameter values).
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#'
#' @return A vector with summarized estimated (or true) values for a specified 
#' infection subgroup. If sim_data is NULL, estimates for each posterior draw 
#' are returned.
compute_avg_groups <- function(dt,
                               stat2,
                               param,
                               category_ids,
                               param_suffix = "",
                               sim_data = FALSE) {
  assert_that(stat2 %in% c("mean", "median"))
  stat_func <- ifelse(stat2 == "mean", mean, median)
  param <- paste0(param, param_suffix)

  # Loop over each category to compare with the reference
  if (sim_data) {
    # Compute medians for the current category for simulated data therefore not 
    # by draws (i.e. we don't have draws because these are the true values, 
    # not the estimated ones)
    summary <- dt[, .(summary_cat = stat_func(eval(as.name(param))[ID %in% category_ids]))]
  } else {
    # Compute medians for the current category, per draw
    summary <- dt[, .(summary_cat = stat_func(eval(as.name(param))[ID %in% category_ids])),
      by = .draw
    ]
  }
  estimated <- summary$summary_cat

  return(estimated)
}


#' Compute the summary statistic of the difference between the reference 
#' infection subgroup and another specified subgroup.
#'
#' @param dt A data.table of posterior draws, grouped by infection identifiers 
#' (ID).
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param param The parameter of interest (e.g. intercept).
#' @param ref_category_ids The infection identifiers that belong to the 
#' reference infection subgroup.
#' @param other_category_ids The infection identifiers that belong to the 
#' other infection subgroup.
#' @param param_suffix A string specifying the suffix of the parameter of 
#' interest (e.g. "_true" if we are summarizing over true
#' parameter values).
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#'
#' @return A vector with summarized estimated (or true) values for a specified 
#' infection subgroup. If sim_data is NULL, estimates for each posterior draw 
#' are returned.
compute_avg_group_differences <- function(dt,
                                          stat2,
                                          param,
                                          ref_category_ids,
                                          other_category_ids,
                                          param_suffix = "",
                                          sim_data = FALSE) {
  assert_that(stat2 %in% c("mean", "median"))
  stat_func <- ifelse(stat2 == "mean", mean, median)
  param <- paste0(param, param_suffix)

  # Loop over each category to compare with the reference
  if (sim_data) {
    # Compute medians for the reference and current category for simulated data
    # therefore not by draws
    summary <- dt[, .(
      summary_ref = stat_func(eval(as.name(param))[ID %in% ref_category_ids]),
      summary_cat = stat_func(eval(as.name(param))[ID %in% other_category_ids])
    )]
  } else {
    # Compute medians for the reference and current category, per draw
    summary <- dt[, .(
      summary_ref = stat_func(eval(as.name(param))[ID %in% ref_category_ids]),
      summary_cat = stat_func(eval(as.name(param))[ID %in% other_category_ids])
    ), by = .draw]
  }
  # Compute the difference for each draw
  diff_estimated <- summary$summary_cat - summary$summary_ref

  return(diff_estimated)
}


#' Create and return a list of estimated (or true) quantities (intercept, 
#' up-slope etc.) for infection subgroups.
#'
#' @param stat_group_params A nested named list containing subgroup-specific 
#' results.
#' @param PAMS3 A logical value indicating whether predictions were computed 
#' for different disease severity statuses according to variable "PAMS3" as 
#' opposed to "PAMS1". "PAMS1" only differentiates between non-PAMS and PAMS, 
#' while "PAMS3" has 3 categories (unknown status, hospitalized, PAMS).
#'
#'
#' @return A list of estimated (or true) quantities (intercept, up-slope etc.) 
#' for infection subgroups.
process_avg_group_params_list <- function(stat_group_params, PAMS3 = FALSE) {
  if (PAMS3) {
    group_var_mapping <- list(
      PAMS3 = c(3, 2, 1), gender = c(0, 1),
      prior_infection = c(0, 1),
      age_category3_code = c(1, 2, 3),
      variant = c(
        "Wildtype", "Alpha", "Delta", "Omicron",
        "Unknown"
      )
    )
  } else {
    group_var_mapping <- list(
      PAMS1 = c(0, 1), gender = c(0, 1),
      hospitalized = c(0, 1), prior_infection = c(0, 1),
      age_category3_code = c(1, 2, 3),
      variant = c(
        "Wildtype", "Alpha", "Delta", "Omicron",
        "Unknown"
      )
    )
  }
  param_list <- list()
  for (group_var in names(group_var_mapping)) {
    cats <- group_var_mapping[[group_var]]
    for (cat in cats) {
      cat_str <- paste0(cat, "_str")
      for (param in c(
        "intercept", "slope_up", "slope_down", "time2peak",
        "time_from_peak"
      )) {
        param_list <- append(param_list, 
                             list(stat_group_params[[group_var]][[cat_str]][[param]]))
      }
    }
  }
  return(param_list)
}


#' Create and return a list of estimated (or true) predicted differences 
#' (in intercept, up-slope etc.) between reference infection subgroups and 
#' remaining infection subgroups.
#'
#' @param stat_group_params A nested named list containing subgroup-specific 
#' results.
#' @param PAMS3 A logical value indicating whether predictions were computed 
#' for different disease severity statuses according to variable "PAMS3"
#' as opposed to "PAMS1". "PAMS1" only differentiates between non-PAMS and 
#' PAMS, while "PAMS3" has 3 categories (unknown status, hospitalized, PAMS).
#'
#' @return A list of estimated (or true) predicted differences between reference
#'  infection subgroups and remaining subgroups.
process_avg_group_params_diff_list <- function(stat_group_params,
                                               PAMS3 = FALSE) {
  if (PAMS3) {
    group_var_mapping <- list(
      PAMS3 = c(3, 2, 1), gender = c(0, 1),
      prior_infection = c(0, 1),
      age_category3_code = c(1, 2, 3),
      variant = c(
        "Wildtype", "Alpha", "Delta", "Omicron",
        "Unknown"
      )
    )
  } else {
    group_var_mapping <- list(
      PAMS1 = c(FALSE, TRUE), gender = c(0, 1),
      hospitalized = c(0, 1), prior_infection = c(0, 1),
      age_category3_code = c(1, 2, 3),
      variant = c(
        "Wildtype", "Alpha", "Delta", "Omicron",
        "Unknown"
      )
    )
  }
  param_list <- list()
  for (group_var in names(group_var_mapping)) {
    cats <- group_var_mapping[[group_var]]
    for (cat in cats[2:length(cats)]) {
      for (param in c(
        "intercept", "slope_up", "slope_down", "time2peak",
        "time_from_peak"
      )) {
        param_list <- append(param_list, 
                             list(stat_group_params[[group_var]][[cat]][[param]]))
      }
    }
  }
  return(param_list)
}

#' Compute summary statistics of differences in posterior estimates between the 
#' reference and the remaining infection subgroups (if sim_data is NULL). 
#' If sim_data is not NULL, compute these values using the known (true) 
#' parameter values. If diff is TRUE, compute the differences between estimated 
#' and true values.
#'
#' @param by_draw_ID A data.table of posterior draws, grouped by infection 
#' identifiers (ID).
#' @param stat_func The function to compute the summary statistic.
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#'
#' @return A nested named list containing subgroup-specific results.
return_main_params_diff <- function(dt,
                                    stat_func,
                                    sim_data) {
  # The following parameters are on the original scale in the model.
  slope_up_stat <- stat_func(dt, "slope_up_true", "slope_up")
  slope_down_stat <- stat_func(dt, "slope_down_true", "slope_down")
  intercept_stat <- stat_func(dt, "intercept_true", "intercept")

  dt[, time2peak_true := intercept_true / slope_up_true]
  dt[, time_from_peak_true := intercept_true / -slope_down_true]


  time2peak_stat <- stat_func(dt, "time2peak_true", "time2peak")
  time_from_peak_stat <- stat_func(dt, "time_from_peak_true", "time_from_peak")

  return(list(
    slope_up = slope_up_stat,
    slope_down = slope_down_stat,
    intercept = intercept_stat,
    time2peak = time2peak_stat,
    time_from_peak = time_from_peak_stat
  ))
}


#' Compute summary statistics of predicted differences between the reference 
#' and the remaining infection subgroups (if sim_data is NULL; predictions are 
#' computed from infection-level parameter estimates). Differences in the 
#' intercept, up-slope, down-slope and proliferation and clearance times 
#' between infection subgroups are computed. If sim_data is not NULL, compute 
#' these values using the known (true) parameter values. If diff is TRUE, 
#' compute the differences between estimated and true values.
#'
#' @param dt A data.table of posterior draws, grouped by infection identifiers 
#' (ID).
#' @param stat_func The function to compute the summary statistic.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param diff A logical value indicating whether to compute the differences 
#' between estimated and true values.
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#' @param PAMS3 A logical value indicating whether to compute predictions for 
#' different disease severity statuses according to variable "PAMS3"
#' as opposed to "PAMS1". "PAMS1" only differentiates between non-PAMS and PAMS, 
#' while "PAMS3" has 3 categories (unknown status, hospitalized, PAMS).
#'
#' @return A nested named list containing subgroup-specific results.
return_group_params_avg_diff <- function(dt,
                                         stat_func,
                                         stat2,
                                         day_data,
                                         diff = FALSE,
                                         sim_data = NULL,
                                         PAMS3 = FALSE) {
  if (diff) {
    assert_that(!is.null(sim_data))
    compute_est <- TRUE
    compute_sim <- TRUE
    param_suffix <- "_true"
  } else if (!is.null(sim_data)) {
    compute_est <- FALSE
    compute_sim <- TRUE
    # This is only to get the right parameter names (see bottom)
    param_suffix <- "" 
  } else {
    compute_est <- TRUE
    compute_sim <- FALSE
    param_suffix <- ""
  }
  # Variable names and categories
  # Define the reference categories first
  if (PAMS3) {
    group_var_mapping <- list(
      PAMS3 = c(3, 2, 1), gender = c(0, 1),
      prior_infection = c(0, 1),
      age_category3_code = c(1, 2, 3),
      variant = c(
        "Wildtype", "Alpha", "Delta", "Omicron",
        "Unknown"
      )
    )
  } else {
    group_var_mapping <- list(
      PAMS1 = c(FALSE, TRUE), gender = c(0, 1),
      hospitalized = c(0, 1), prior_infection = c(0, 1),
      age_category3_code = c(1, 2, 3),
      variant = c(
        "Wildtype", "Alpha", "Delta", "Omicron",
        "Unknown"
      )
    )
  }

  results_est <- list()
  results_sim <- list()
  results_diff <- list()
  for (group_var in names(group_var_mapping)) {
    results_est[[group_var]] <- list()
    results_sim[[group_var]] <- list()
    results_diff[[group_var]] <- list()

    cats <- group_var_mapping[[group_var]]
    ref_cat <- cats[1]
    ref_category_ids <- unique(day_data[eval(as.name(group_var)) == ref_cat]$ID)
    other_cats <- cats[2:length(cats)]
    for (other_cat in other_cats) {
      results_est[[group_var]][[other_cat]] <- list()
      results_sim[[group_var]][[other_cat]] <- list()
      results_diff[[group_var]][[other_cat]] <- list()

      other_category_ids <- unique(day_data[eval(as.name(group_var)) == other_cat]$ID)
      for (param in c("intercept", 
                      "slope_up", 
                      "slope_down", 
                      "time2peak", 
                      "time_from_peak")) {
        if (compute_est) {
          est_diffs <- compute_avg_group_differences(
            dt = dt,
            stat2 = stat2,
            param = param,
            ref_category_ids = ref_category_ids,
            other_category_ids = other_category_ids,
            param_suffix = "",
            sim_data = FALSE
          )
        }
        if (compute_sim) {
          sim_diffs <- compute_avg_group_differences(
            dt = dt,
            stat2 = stat2,
            param = param,
            ref_category_ids = ref_category_ids,
            other_category_ids = other_category_ids,
            param_suffix = param_suffix,
            sim_data = TRUE
          )
        }
        # Compute difference of differences
        # (i.e. estimated differences - true differences)
        if (diff) {
          results_diff[[group_var]][[other_cat]][[param]] <- stat_func(
            sim_diffs,
            est_diffs
          )
        } else if (compute_sim) {
          results_sim[[group_var]][[other_cat]][[param]] <- stat_func(sim_diffs)
        } else {
          results_est[[group_var]][[other_cat]][[param]] <- stat_func(est_diffs)
        }
      }
    }
  }
  if (diff) {
    results <- results_diff
  } else if (compute_sim) {
    results <- results_sim
  } else {
    results <- results_est
  }
  return(results)
}


#' Compute summary statistics of conditional posterior predictions (if sim_data 
#' is NULL; computed from infection-level parameter estimates) of intercept, 
#' up-slope, down-slope, and proliferation and clearance times for different 
#' infection subgroups. If sim_data is not NULL, compute these values using the 
#' known (true) parameter values. If diff is TRUE, compute the differences 
#' between estimated and true values.
#'
#' @param dt A data.table of posterior draws, grouped by infection identifiers 
#' (ID).
#' @param stat_func The function to compute the summary statistic.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param diff A logical value indicating whether to compute the differences 
#' between estimated and true values.
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#' @param PAMS3 A logical value indicating whether to compute predictions for 
#' different disease severity statuses according to variable PAMS3 as opposed 
#' to PAMS1. PAMS1 only differentiates between non-PAMS and PAMS, while PAMS3 
#' has 3 categories (unknown status, hospitalized, PAMS).
#'
#' @return A nested named list with values for the different groups and 
#' parameters (up-slope, intercept etc.).
return_group_params_avg <- function(dt,
                                    stat_func,
                                    stat2,
                                    day_data,
                                    diff = FALSE,
                                    sim_data = NULL,
                                    PAMS3 = FALSE) {
  if (diff) {
    assert_that(!is.null(sim_data))
    compute_est <- TRUE
    compute_sim <- TRUE
    param_suffix <- "_true"
  } else if (!is.null(sim_data)) {
    compute_est <- FALSE
    compute_sim <- TRUE
    # This is only to get the right parameter names (see bottom)
    param_suffix <- ""
  } else {
    compute_est <- TRUE
    compute_sim <- FALSE
    param_suffix <- ""
  }
  # Variable names and categories
  # Define the reference categories first
  if (PAMS3) {
    group_var_mapping <- list(
      PAMS3 = c(3, 2, 1), gender = c(0, 1),
      prior_infection = c(0, 1),
      age_category3_code = c(1, 2, 3),
      variant = c(
        "Wildtype", "Alpha", "Delta", "Omicron",
        "Unknown"
      )
    )
  } else {
    group_var_mapping <- list(
      PAMS1 = c(0, 1), gender = c(0, 1),
      hospitalized = c(0, 1), prior_infection = c(0, 1),
      age_category3_code = c(1, 2, 3),
      variant = c(
        "Wildtype", "Alpha", "Delta", "Omicron",
        "Unknown"
      )
    )
  }

  results_est <- list()
  results_sim <- list()
  results_diff <- list()
  for (group_var in names(group_var_mapping)) {
    results_est[[group_var]] <- list()
    results_sim[[group_var]] <- list()
    results_diff[[group_var]] <- list()

    cats <- group_var_mapping[[group_var]]
    for (cat in cats) {
      cat_str <- paste0(cat, "_str")
      results_est[[group_var]][[cat_str]] <- list()
      results_sim[[group_var]][[cat_str]] <- list()
      results_diff[[group_var]][[cat_str]] <- list()

      category_ids <- unique(day_data[eval(as.name(group_var)) == cat]$ID)
      for (param in c("intercept", 
                      "slope_up", 
                      "slope_down", 
                      "time2peak", 
                      "time_from_peak")) {
        if (compute_est) {
          est <- compute_avg_groups(
            dt = dt,
            stat2 = stat2,
            param = param,
            category_ids = category_ids,
            param_suffix = "",
            sim_data = FALSE
          )
        }
        if (compute_sim) {
          sim <- compute_avg_groups(
            dt = dt,
            stat2 = stat2,
            param = param,
            category_ids = category_ids,
            param_suffix = param_suffix,
            sim_data = TRUE
          )
        }
        # Compute difference of differences
        # (i.e. estimated differences - true differences)
        if (diff) {
          results_diff[[group_var]][[cat_str]][[param]] <- stat_func(
            sim,
            est
          )
        } else if (compute_sim) {
          results_sim[[group_var]][[cat_str]][[param]] <- stat_func(sim)
        } else {
          results_est[[group_var]][[cat_str]][[param]] <- stat_func(est)
        }
      }
    }
  }
  if (diff) {
    results <- results_diff
  } else if (compute_sim) {
    results <- results_sim
  } else {
    results <- results_est
  }
  return(results)
}


#' Return parameter values (estimates, true values, differences between 
#' estimated and true values) to be written to a common .tsv file.
#'
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param n_min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
#' @param n_samples A numeric (or NULL) specifying the number of infections 
#' that were simulated (if sim_data not NULL, else NULL).
#' @param sim A logical value indicating whether parameter values are true 
#' values (when fitting simulated data).
#' @param stat_params A named list containing values of the main model 
#' parameters.
#' @param intercept_dict A named list containing values of 
#' infection-group-specific intercept parameters.
#' @param up_dict A named list containing values of infection-group-specific 
#' up-slope parameters.
#' @param down_dict A named list containing values of infection-group-specific 
#' down-slope parameters.
#' @param earliest_sim_day A numeric specifying the earliest day (relative to 
#' day of peak viral load) that a viral load value was simulated for.
#' @param latest_sim_day A numeric specifying the latest day (relative to day 
#' of peak viral load) that a viral load value was simulated for.
#'
#' @return A vector of parameter values.
return_params_vec <- function(model,
                              model_no_stan,
                              n_min_tests,
                              n_samples,
                              sim,
                              stat_params,
                              intercept_dict,
                              up_dict,
                              down_dict,
                              earliest_sim_day,
                              latest_sim_day) {
  if (sim) {
    assert_that(!is.na(earliest_sim_day))
    assert_that(!is.na(latest_sim_day))
  }
  assert_that(!is.null(earliest_sim_day))
  assert_that(!is.null(latest_sim_day))
  sim_int <- int(sim)

  group_params <- process_avg_group_params_list(stat_params$avg_group_params)
  group_params_diff <- process_avg_group_params_diff_list(stat_params$avg_group_params_diff)

  list_top <- list(
    model,
    model_no_stan,
    n_min_tests,
    ifelse(is.null(n_samples), NA, n_samples),
    sim_int,
    stat_params$alpha_mu,
    stat_params$intercept_mu,
    stat_params$slope_up_mu,
    stat_params$slope_down_mu,
    stat_params$sigma_mu,
    stat_params$sigma_sigma,
    stat_params$sigma,
    stat_params$intercept,
    stat_params$slope_up,
    stat_params$slope_down,
    stat_params$intercept_sigma,
    stat_params$slope_up_sigma,
    stat_params$slope_down_sigma,
    stat_params$time2peak,
    stat_params$time_from_peak,
    earliest_sim_day,
    latest_sim_day,
    intercept_dict[["PAMS"]],
    up_dict[["PAMS"]],
    down_dict[["PAMS"]],
    intercept_dict[["gender"]],
    up_dict[["gender"]],
    down_dict[["gender"]],
    intercept_dict[["hospitalized"]],
    up_dict[["hospitalized"]],
    down_dict[["hospitalized"]],
    intercept_dict[["prior_infection"]],
    up_dict[["prior_infection"]],
    down_dict[["prior_infection"]],
    intercept_dict[["age_cat2"]],
    up_dict[["age_cat2"]],
    down_dict[["age_cat2"]],
    intercept_dict[["age_cat3"]],
    up_dict[["age_cat3"]],
    down_dict[["age_cat3"]],
    intercept_dict[["alpha"]],
    up_dict[["alpha"]],
    down_dict[["alpha"]],
    intercept_dict[["delta"]],
    up_dict[["delta"]],
    down_dict[["delta"]],
    intercept_dict[["omicron"]],
    up_dict[["omicron"]],
    down_dict[["omicron"]],
    intercept_dict[["unknown"]],
    up_dict[["unknown"]],
    down_dict[["unknown"]]
  )
  list_top2 <- list(
    intercept_dict[["hospitalized:alpha"]],
    up_dict[["hospitalized:alpha"]],
    down_dict[["hospitalized:alpha"]],
    intercept_dict[["hospitalized:delta"]],
    up_dict[["hospitalized:delta"]],
    down_dict[["hospitalized:delta"]],
    intercept_dict[["hospitalized:omicron"]],
    up_dict[["hospitalized:omicron"]],
    down_dict[["hospitalized:omicron"]],
    intercept_dict[["hospitalized:unknown"]],
    up_dict[["hospitalized:unknown"]],
    down_dict[["hospitalized:unknown"]]
  )

  final_list <- c(list_top, list_top2, group_params, group_params_diff)


  return(final_list)
}


#' Return the true parameter values used for simulating data (which were in 
#' turn used for fitting).
#'
#' @param model A string specifying the suffix of the model (includes numeric 
#' model identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param n_min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
#' @param n_samples A numeric (or NULL) specifying the number of infections 
#' that were simulated (if sim_data not NULL, else NULL).
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#' @param stat_params_sim A named list of the true parameter values used for 
#' simulating data.
#' @param ncols A numeric indicating the number of columns intended for the 
#' output file.
#' @param exponentiate A logical value specifying whether parameter estimates 
#' were exponentiated where it was appropriate (e.g. fixed-effects parameters).
#' @return A string with the true parameter values used for simulating data.
return_sim_params_string <- function(model,
                                     model_no_stan,
                                     n_min_tests,
                                     n_samples,
                                     sim_data,
                                     stat_params_sim,
                                     ncols,
                                     exponentiate) {
  sim_params_groups <- calc_sim_group_stats(sim_data$params, 
                                            exponentiate = exponentiate)
  intercept_dict <- sim_params_groups[["intercept"]]
  slope_up_dict <- sim_params_groups[["slope_up"]]
  slope_down_dict <- sim_params_groups[["slope_down"]]
  params_sim_vec <- return_params_vec(
    model = model,
    model_no_stan = model_no_stan,
    n_min_tests = n_min_tests,
    n_samples = n_samples,
    sim = 1,
    stat_params = stat_params_sim,
    intercept_dict = intercept_dict,
    up_dict = slope_up_dict,
    down_dict = slope_down_dict,
    earliest_sim_day = min(sim_data$dt$day_sampled),
    latest_sim_day = max(sim_data$dt$day_sampled)
  )

  assert_that(ncols == length(params_sim_vec))
  params_sim <- paste(params_sim_vec, collapse = "\t")
  return(params_sim)
}


#' Return the header of the tsv file where posterior estimates and predictions 
#' will be stored. This file contains either
#' A) posterior estimates and conditional predictions (computed from 
#' infection-level parameter estimates), if post_pred_groups_inf_only is set 
#' `FALSE`,
#' B) only posterior predictions (computed from population-level parameter 
#' estimates), if post_pred_groups_inf_only is TRUE
#'
#' @param exponentiate A logical value specifying whether parameter estimates 
#' were exponentiated where it is appropriate (e.g. fixed-effects parameters).
#' @param post_pred_groups_inf_only A logical value indicating whether to save 
#' only conditional posterior predictions (computed from infection-level 
#' parameter estimates), as opposed to saving both posterior estimates and 
#' posterior predictions.
return_params_header <- function(exponentiate = FALSE,
                                 post_pred_groups_inf_only = FALSE) {
  if (post_pred_groups_inf_only) {
    subgroup2_string_vec <- c(
      "non_pams_non_hosp_avg", "hospitalized_avg",
      "pams_avg", "female_avg", "male_avg",
      "no_prior_infection_avg", "prior_infection_avg",
      "age_cat1_avg", "age_cat2_avg", "age_cat3_avg",
      "wildtype_avg", "alpha_avg", "delta_avg",
      "omicron_avg", "unknown_avg"
    )
    subgroup3_string_vec <- c(
      "hospitalized_avg_diff",
      "pams_avg_diff", "gender_avg_diff",
      "prior_infection_avg_diff",
      "age_cat2_avg_diff", "age_cat3_avg_diff",
      "alpha_avg_diff", "delta_avg_diff",
      "omicron_avg_diff", "unknown_avg_diff"
    )

    header_subgroups2 <- paste(
      c(sapply(
        subgroup2_string_vec,
        function(subgroup_string) {
          c(
            paste0("intercept_", subgroup_string),
            paste0("slope_up_", subgroup_string),
            paste0("slope_down_", subgroup_string),
            paste0("time2peak_", subgroup_string),
            paste0("time_from_peak_", subgroup_string)
          )
        }
      )),
      collapse = "\t"
    )
    header_subgroups3 <- paste(
      c(sapply(
        subgroup3_string_vec,
        function(subgroup_string) {
          c(
            paste0("intercept_", subgroup_string),
            paste0("slope_up_", subgroup_string),
            paste0("slope_down_", subgroup_string),
            paste0("time2peak_", subgroup_string),
            paste0("time_from_peak_", subgroup_string)
          )
        }
      )),
      collapse = "\t"
    )
    header_base <- paste(c("model", "model_no_stan", "n_min_tests"), 
                         collapse = "\t")
    header <- paste(header_base, header_subgroups2, header_subgroups3, 
                    sep = "\t")
    ncols <- 4 + 5 * length(subgroup2_string_vec) + 5 * length(subgroup3_string_vec)
  } else {
    subgroup_string_vec <- c(
      "pams", "gender", "hospitalized", "prior_infection",
      "age_cat2", "age_cat3", "alpha", "delta", "omicron",
      "unknown", "hospitalized_alpha", "hospitalized_delta",
      "hospitalized_omicron", "hospitalized_unknown"
    )

    subgroup2_string_vec <- c(
      "non_pams_avg", "pams_avg", "female_avg", "male_avg",
      "non_hospitalized_avg", "hospitalized_avg",
      "no_prior_infection_avg", "prior_infection_avg",
      "age_cat1_avg", "age_cat2_avg", "age_cat3_avg",
      "wildtype_avg", "alpha_avg", "delta_avg", "omicron_avg",
      "unknown_avg"
    )
    subgroup3_string_vec <- c(
      "pams_avg_diff", "gender_avg_diff",
      "hospitalized_avg_diff", "prior_infection_avg_diff",
      "age_cat2_avg_diff", "age_cat3_avg_diff",
      "alpha_avg_diff", "delta_avg_diff",
      "omicron_avg_diff", "unknown_avg_diff"
    )

    exp_prefix <- ifelse(exponentiate, "", "log_")
    header_vec <- c(
      "model", "model_no_stan", "n_min_tests",
      "n_samples", "params_true", "alpha_mu",
      paste0(exp_prefix, "intercept_mu"),
      paste0(exp_prefix, "slope_up_mu"),
      paste0(exp_prefix, "slope_down_mu"),
      paste0(exp_prefix, "sigma_mu"),
      paste0(exp_prefix, "sigma_sigma"),
      "sigma", "intercept", "slope_up",
      "slope_down", paste0(exp_prefix, "intercept_sigma"),
      paste0(exp_prefix, "slope_up_sigma"),
      paste0(exp_prefix, "slope_down_sigma"),
      "time2peak", "time_from_peak",
      "earliest_day_sim", "latest_day_sim"
    )
    header <- paste(header_vec, collapse = "\t")
    intercept_prefix <- paste0(exp_prefix, "intercept_")
    up_prefix <- paste0(exp_prefix, "slope_up_")
    down_prefix <- paste0(exp_prefix, "slope_down_")
    header_subgroups <- paste(c(sapply(
      subgroup_string_vec,
      function(subgroup_string) {
        c(
          paste0(intercept_prefix, subgroup_string),
          paste0(up_prefix, subgroup_string),
          paste0(down_prefix, subgroup_string)
        )
      }
    )), collapse = "\t")


    header_subgroups2 <- paste(
      c(sapply(
        subgroup2_string_vec,
        function(subgroup_string) {
          c(
            paste0("intercept_", subgroup_string),
            paste0("slope_up_", subgroup_string),
            paste0("slope_down_", subgroup_string),
            paste0("time2peak_", subgroup_string),
            paste0("time_from_peak_", subgroup_string)
          )
        }
      )),
      collapse = "\t"
    )
    header_subgroups3 <- paste(
      c(sapply(
        subgroup3_string_vec,
        function(subgroup_string) {
          c(
            paste0("intercept_", subgroup_string),
            paste0("slope_up_", subgroup_string),
            paste0("slope_down_", subgroup_string),
            paste0("time2peak_", subgroup_string),
            paste0("time_from_peak_", subgroup_string)
          )
        }
      )),
      collapse = "\t"
    )

    header <- paste(header, 
                    header_subgroups, 
                    header_subgroups2, 
                    header_subgroups3, sep = "\t")
    ncols <- length(header_vec) + 3 * length(subgroup_string_vec) + 5 * length(subgroup2_string_vec) + 5 * length(subgroup3_string_vec)
  }
  return(list("header" = header, "ncols" = ncols))
}


#' Compute and return differences between estimated and true 
#' infection-subgroup-specific parameter/fixed-effects values (mean, median or 
#' HPDIs).
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param sim_params_groups A named list of the true infection-group-specific 
#' parameter values used for simulating data.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed.
#' @param exponentiate A logical value specifying whether parameter estimates 
#' should be exponentiated where appropriate (e.g. fixed-effects parameters).
#' @return A named list of differences between estimated and true main 
#' parameter values (mean, median or HPDIs).
calc_group_params_diff <- function(draws,
                                   day_data,
                                   sim_params_groups,
                                   stat = "median",
                                   exponentiate = FALSE) {
  assert_that(stat %in% c("median", "hdi", "hdi50", "mean"))
  if (exponentiate) {
    stat_func <- function(param_true, param_est) {
      return(calc_relative_difference(
        param_true = param_true,
        param_est = param_est,
        summary_stat = stat,
        exponentiate = TRUE
      ))
    }
  } else {
    stat_func <- function(param_true, param_est) {
      return(calc_difference(
        param_true = param_true,
        param_est = param_est,
        summary_stat = stat
      ))
    }
  }

  mapping <- SUBGRP_MAPPING
  # True parameter values
  intercept_dict_sim <- sim_params_groups[["intercept"]]
  slope_up_dict_sim <- sim_params_groups[["slope_up"]]
  slope_down_dict_sim <- sim_params_groups[["slope_down"]]
  # Estimated parameters
  # Convert the draws object to a data.table
  draws_dt <- as_draws_dt(draws, thin = 1)
  intercept_dict <- list()
  slope_up_dict <- list()
  slope_down_dict <- list()
  for (group_number in names(mapping)) {
    group_name <- mapping[[group_number]]
    intercept_param_current <- paste0("betaPGH_intercept[", group_number, "]")
    slope_up_param_current <- paste0("betaPGH_slope_up[", group_number, "]")
    slope_down_param_current <- paste0("betaPGH_slope_down[", group_number, "]")

    intercept_dict[[group_name]] <- stat_func(intercept_dict_sim[[group_name]], 
                                              draws_dt[[intercept_param_current]])
    slope_up_dict[[group_name]] <- stat_func(slope_up_dict_sim[[group_name]], 
                                             draws_dt[[slope_up_param_current]])
    slope_down_dict[[group_name]] <- stat_func(slope_down_dict_sim[[group_name]], 
                                               draws_dt[[slope_down_param_current]])
  }
  return(list("intercept" = intercept_dict, 
              "slope_up" = slope_up_dict, 
              "slope_down" = slope_down_dict))
}


#' Compute and return differences between estimated and true main parameter 
#' values (mean, median or HPDIs), and differences between 
#' infection-group-specific posterior predictions (computed from 
#' infection-level parameter estimates) and true values. This function is meant 
#' for runs with simulated data.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#' @param sim_data_main_params A named list of the true main parameter values 
#' used for simulating data.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param exponentiate A logical value specifying whether parameter estimates 
#' should be exponentiated where appropriate (e.g. fixed-effects
#' parameters).
#' @return A named list of differences between estimated and true main 
#' parameter values (mean, median or HPDIs) and differences between
#' infection-group-specific posterior predictions (computed from 
#' infection-level parameter estimates) and true values. This function is meant
#' for runs with simulated data.
calc_main_params_diff <- function(draws,
                                  day_data,
                                  model_no_stan,
                                  sim_data,
                                  sim_data_main_params,
                                  stat = "median",
                                  stat2 = "median",
                                  exponentiate = FALSE) {
  assert_that(stat %in% c("median", "hdi", "hdi50", "mean"))
  # stat2 is used for averaging over infections first (for calculating 
  # infection group averages), then computing stat over draws

  if (exponentiate) {
    stat_func <- function(param_true, param_est) {
      return(calc_relative_difference(
        param_true = param_true,
        param_est = param_est,
        summary_stat = stat,
        exponentiate = TRUE
      ))
    }
    # These are main model parameters like intercept, time2peak etc., 
    # therefore no need to exponentiate
    stat_func2 <- function(dt, param_name_true, param_name_est) {
      return(calc_relative_difference_pop_avg(dt,
        param_name_true = param_name_true,
        param_name_est = param_name_est,
        summary_stat = stat,
        summary_stat2 = stat2,
        exponentiate = FALSE
      ))
    }
    stat_func3 <- function(param_true, param_est) {
      return(calc_relative_difference(
        param_true = param_true,
        param_est = param_est,
        summary_stat = stat,
        exponentiate = FALSE
      ))
    }
  } else {
    stat_func <- function(param_true, param_est) {
      return(calc_difference(param_true, param_est, summary_stat = stat))
    }
    stat_func2 <- function(dt, param_name_true, param_name_est) {
      return(calc_difference_pop_avg(dt,
        param_name_true,
        param_name_est,
        summary_stat = stat,
        summary_stat2 = stat2
      ))
    }
    stat_func3 <- stat_func
  }

  by_draw_ID <-
    draws_by_id(draws,
      day_data,
      c("slope_up", 
        "slope_down", 
        "intercept", 
        "sigma", 
        "time2peak", 
        "time_from_peak"),
      thin = 1
    )
  # Convert the draws object to a data.table
  draws_dt <- as_draws_dt(draws, thin = 1)

  day_data_ID <- day_data[, .(ID, up_slope, down_slope, intercept, sigma)]
  # Keep only one row per ID (the first occurrence)
  day_data_ID <- unique(day_data_ID, by = "ID")
  # Add time2peak and time_from_peak columns
  L <- sim_data$dt$L

  day_data_ID[, time2peak_true := intercept / up_slope]
  day_data_ID[, time_from_peak_true := intercept / -down_slope]

  # Rename a single column
  setnames(day_data_ID, old = "intercept", new = "intercept_true")
  setnames(day_data_ID, old = "up_slope", new = "slope_up_true")
  setnames(day_data_ID, old = "down_slope", new = "slope_down_true")
  setnames(day_data_ID, old = "sigma", new = "sigma_true")


  stopifnot(all(day_data_ID$ID %in% by_draw_ID$ID)) # every ID in dt1 is in dt2
  stopifnot(all(by_draw_ID$ID %in% day_data_ID$ID)) # every ID in dt2 is in dt1

  by_draw_ID <- merge(by_draw_ID, day_data_ID, by = "ID")

  # The following parameters are on the original scale in the model.
  main_params_list <- return_main_params_diff(
    dt = by_draw_ID,
    stat_func = stat_func2,
    sim_data = sim_data
  )
  avg_group_params_list <- return_group_params_avg(
    dt = by_draw_ID,
    stat_func = stat_func3,
    stat2 = stat2,
    day_data = day_data,
    sim_data = sim_data,
    diff = TRUE
  )
  avg_group_params_diff_list <- return_group_params_avg_diff(
    dt = by_draw_ID,
    stat_func = stat_func3,
    stat2 = stat2,
    day_data = day_data,
    sim_data = sim_data,
    diff = TRUE
  )

  sigma_stat <- stat_func2(by_draw_ID, "sigma_true", "sigma")
  alpha_mu_stat <- stat_func3(sim_data$params$alpha_mu[1], draws_dt$alpha_mu)

  # Parameters on log-scale.
  intercept_mu_stat <- stat_func(sim_data_main_params$intercept_mu, 
                                 draws_dt$log_intercept_mu)
  slope_up_mu_stat <- stat_func(sim_data_main_params$slope_up_mu, 
                                draws_dt$log_slope_up_mu)
  slope_down_mu_stat <- stat_func(sim_data_main_params$slope_down_mu, 
                                  draws_dt$log_slope_down_mu)
  sigma_mu_stat <- stat_func(0, draws_dt$sigma_mu)
  sigma_sigma_stat <- stat_func(0.3, draws_dt$sigma_sigma)

  # Extract values for sigma params
  slope_up_sigma_stat <- stat_func(sim_data$params$sigma_up_slope, 
                                   draws_dt$slope_up_sigma)
  slope_down_sigma_stat <- stat_func(sim_data$params$sigma_down_slope, 
                                     draws_dt$slope_down_sigma)
  intercept_sigma_stat <- stat_func(sim_data$params$sigma_intercept, 
                                    draws_dt$intercept_sigma)


  return(list(
    slope_up = main_params_list$slope_up,
    slope_down = main_params_list$slope_down,
    intercept = main_params_list$intercept,
    sigma = sigma_stat,
    time2peak = main_params_list$time2peak,
    time_from_peak = main_params_list$time_from_peak,
    slope_up_sigma = slope_up_sigma_stat,
    slope_down_sigma = slope_down_sigma_stat,
    intercept_sigma = intercept_sigma_stat,
    intercept_mu = intercept_mu_stat,
    slope_up_mu = slope_up_mu_stat,
    slope_down_mu = slope_down_mu_stat,
    alpha_mu = alpha_mu_stat,
    sigma_mu = sigma_mu_stat,
    sigma_sigma = sigma_sigma_stat,
    avg_group_params = avg_group_params_list,
    avg_group_params_diff = avg_group_params_diff_list
  ))
}


#' Return a named list of infection-group-specific/fixed-effects parameter 
#' estimates (mean, median or HPDIs).
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed.
#' @param exponentiate A logical value specifying whether parameter estimates 
#' should be exponentiated where appropriate (e.g. fixed-effects parameters).
#'
#' @return A named list.
calc_group_params_summary_stat <- function(draws,
                                           stat = "median",
                                           exponentiate = FALSE) {
  assert_that(stat %in% c("median", "hdi", "hdi50", "mean"))
  if (stat == "median") {
    stat_func <- median
  } else if (stat == "mean") {
    stat_func <- mean
  } else if (stat == "hdi50") {
    stat_func <- calc_50_hdi
  } else {
    stat_func <- calc_94_hdi
  }
  if (exponentiate) {
    exp_func <- exp
  } else {
    # Do nothing.
    exp_func <- function(x) {
      return(x)
    }
  }

  # Convert the draws object to a data.table
  draws_dt <- as_draws_dt(draws, thin = 1)

  mapping <- c(
    "PAMS", "gender", "hospitalized", "prior_infection", "age_cat2", "age_cat3",
    "alpha", "delta", "omicron", "unknown", "hospitalized:alpha",
    "hospitalized:delta", "hospitalized:omicron", "hospitalized:unknown"
  )
  names(mapping) <- 1:14

  intercept_dict <- list()
  up_dict <- list()
  down_dict <- list()

  log_intercept_mu <- draws_dt[["log_intercept_mu"]]
  log_slope_up_mu <- draws_dt[["log_slope_up_mu"]]
  log_slope_down_mu <- draws_dt[["log_slope_down_mu"]]


  for (group_number in names(mapping)) {
    group_name <- mapping[[group_number]]
    intercept_param_current <- paste0("betaPGH_intercept[", group_number, "]")
    slope_up_param_current <- paste0("betaPGH_slope_up[", group_number, "]")
    slope_down_param_current <- paste0("betaPGH_slope_down[", group_number, "]")

    intercept_dict[[group_name]] <- exp_func(stat_func(draws_dt[[intercept_param_current]]))
    up_dict[[group_name]] <- exp_func(stat_func(draws_dt[[slope_up_param_current]]))
    down_dict[[group_name]] <- exp_func(stat_func(draws_dt[[slope_down_param_current]]))
  }

  return(list("intercept" = intercept_dict, "up" = up_dict, "down" = down_dict))
}


#' Return a named list of main parameter estimates (mean, median or HPDIs) and 
#' infection-group-specific posterior predictions (computed from infection-level 
#' parameter estimates).
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed for posterior estimates and predictions.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param exponentiate A logical value specifying whether parameter estimates 
#' should be exponentiated where appropriate (e.g. fixed-effects parameters).
#'
#' @return A named list of main parameter estimates (mean, median or HPDIs) and 
#' infection-group-specific posterior predictions (computed from infection-level 
#' parameter estimates).
calc_main_params_summary_stat <- function(draws,
                                          day_data,
                                          stat = "median",
                                          stat2 = "median",
                                          exponentiate = FALSE) {
  assert_that(stat %in% c("median", "hdi", "hdi50", "mean"))
  assert_that(stat2 %in% c("median", "mean"))
  if (stat == "median") {
    stat_func <- function(dt, param_name) {
      median(calc_stat_by_draw(dt, param_name, summary_stat = "median"))
    }
    stat_func2 <- median
  } else if (stat == "mean") {
    stat_func <- function(dt, param_name) {
      mean(calc_stat_by_draw(dt, param_name, summary_stat = "mean"))
    }
    stat_func2 <- mean
  } else if (stat == "hdi50") {
    stat_func <- function(dt, param_name) {
      calc_50_hdi(calc_stat_by_draw(dt, param_name, summary_stat = stat2))
    }
    stat_func2 <- calc_50_hdi
  } else {
    stat_func <- function(dt, param_name) {
      calc_94_hdi(calc_stat_by_draw(dt, param_name, summary_stat = stat2))
    }
    stat_func2 <- calc_94_hdi
  }
  if (exponentiate) {
    exp_func <- exp
  } else {
    # Do nothing.
    exp_func <- function(x) {
      return(x)
    }
  }

  by_draw_ID <-
    draws_by_id(draws,
      day_data,
      c(
        "slope_up", "slope_down", "intercept", "sigma", "time2peak",
        "time_from_peak"
      ),
      thin = 1
    )
  # Convert the draws object to a data.table
  draws_dt <- as_draws_dt(draws, thin = 1)

  slope_up_stat <- stat_func(by_draw_ID, "slope_up")
  slope_down_stat <- stat_func(by_draw_ID, "slope_down")
  intercept_stat <- stat_func(by_draw_ID, "intercept")
  sigma_stat <- stat_func(by_draw_ID, "sigma")
  time2peak_stat <- stat_func(by_draw_ID, "time2peak")
  time_from_peak_stat <- stat_func(by_draw_ID, "time_from_peak")
  alpha_mu_stat <- stat_func2(draws_dt$alpha_mu)

  avg_group_params_list <- return_group_params_avg(
    dt = by_draw_ID,
    stat_func = stat_func2,
    stat2 = stat2,
    day_data = day_data,
    sim_data = NULL,
    diff = FALSE
  )

  avg_group_params_diff_list <- return_group_params_avg_diff(
    dt = by_draw_ID,
    stat_func = stat_func2,
    stat2 = stat2,
    day_data = day_data,
    sim_data = NULL,
    diff = FALSE
  )

  log_intercept_mu_stat <- stat_func2(draws_dt$log_intercept_mu)
  log_slope_up_mu_stat <- stat_func2(draws_dt$log_slope_up_mu)
  log_slope_down_mu_stat <- stat_func2(draws_dt$log_slope_down_mu)
  sigma_mu_stat <- stat_func2(draws_dt$sigma_mu)
  sigma_sigma_stat <- stat_func2(draws_dt$sigma_sigma)

  # Extract values for sigma params, exponentiate, and calculate the stat_func
  slope_up_sigma_stat <- stat_func2(draws_dt$slope_up_sigma)
  slope_down_sigma_stat <- stat_func2(draws_dt$slope_down_sigma)
  intercept_sigma_stat <- stat_func2(draws_dt$intercept_sigma)


  return(list(
    slope_up = slope_up_stat,
    slope_down = slope_down_stat,
    intercept = intercept_stat,
    sigma = sigma_stat,
    time2peak = time2peak_stat,
    time_from_peak = time_from_peak_stat,
    alpha_mu = alpha_mu_stat,
    slope_up_sigma = exp_func(slope_up_sigma_stat),
    slope_down_sigma = exp_func(slope_down_sigma_stat),
    intercept_sigma = exp_func(intercept_sigma_stat),
    intercept_mu = exp_func(log_intercept_mu_stat),
    slope_up_mu = exp_func(log_slope_up_mu_stat),
    slope_down_mu = exp_func(log_slope_down_mu_stat),
    sigma_mu = exp_func(sigma_mu_stat),
    sigma_sigma = exp_func(sigma_sigma_stat),
    avg_group_params = avg_group_params_list,
    avg_group_params_diff = avg_group_params_diff_list
  ))
}


#' Return a named list of main parameter values used for simulating viral load 
#' data, as well as infection-group-specific quantities (up-slope, intercept, 
#' down-slope), computed by averaging over infection-level parameters of 
#' infections belonging to the respective group.
#'
#'
#' @param stat_params_sim A named list of the true parameter values used for 
#' simulating data.
#' @param exponentiate A logical value specifying whether parameter estimates 
#' should be exponentiated where appropriate (e.g. fixed-effects parameters).
#'
#' @return A named list of main parameter values and infection-group-specific 
#' quantities (up-slope, intercept, down-slope ) used for simulating data.
calc_main_params_summary_stat_sim_data <- function(sim_data,
                                                   day_data,
                                                   model_no_stan,
                                                   stat = "median",
                                                   stat2 = "median",
                                                   exponentiate = FALSE) {
  assert_that(stat %in% c("median", "hdi", "hdi50", "mean"))
  if (stat == "median") {
    stat_func <- median
  } else if (stat == "mean") {
    stat_func <- mean
  } else if (stat == "hdi50") {
    stat_func <- calc_50_hdi
  } else {
    stat_func <- calc_94_hdi
  }
  if (exponentiate) {
    exp_func <- exp
  } else {
    # Do nothing.
    exp_func <- function(x) {
      return(x)
    }
  }

  sim_data_dt <- copy(sim_data$dt)
  # Need to rename for computation in return_group_params_avg_diff
  sim_data_dt[, slope_up := up_slope]
  sim_data_dt[, slope_down := down_slope]
  L <- sim_data$dt$L

  sim_data_dt[, time2peak := intercept / up_slope]
  sim_data_dt[, time_from_peak := intercept / -down_slope]

  time2peak_stat <- stat_func(sim_data_dt$time2peak)
  time_from_peak_stat <- stat_func(sim_data_dt$time_from_peak)

  stat_up_slope <- stat_func(sim_data$up_slope)
  stat_down_slope <- stat_func(sim_data$down_slope)
  stat_intercept <- stat_func(sim_data$intercept)
  stat_sigma <- stat_func(sim_data$sigma)

  avg_group_params_list <- return_group_params_avg(
    dt = sim_data_dt,
    stat_func = stat_func,
    stat2 = stat2,
    day_data = day_data,
    sim_data = sim_data,
    diff = FALSE
  )
  avg_group_params_diff_list <- return_group_params_avg_diff(
    dt = sim_data_dt,
    stat_func = stat_func,
    stat2 = stat2,
    day_data = day_data,
    sim_data = sim_data,
    diff = FALSE
  )

  sim_params <- sim_data$params
  param_sizes_all_zero <- sim_params$param_sizes_all_zero

  # random effects sigmas are on the log-scale
  stat_up_slope_sigma <- sim_params$sigma_up_slope
  stat_down_slope_sigma <- sim_params$sigma_down_slope
  stat_intercept_sigma <- sim_params$sigma_intercept

  # We are using different reference categories in our stan model, so we
  # need to adjust accordingly (add these params, the rest that are not added,
  # are 0).
  add_to_intercept <- (ifelse(param_sizes_all_zero, 
                            0, sim_params$log_intercept_pams3[3]) +
    ifelse(param_sizes_all_zero, 0, sim_params$log_intercept_age[1]) +
    ifelse(param_sizes_all_zero, 0, sim_params$log_intercept_variant$Wildtype)
  )
  add_to_slope_up <- (ifelse(param_sizes_all_zero, 
                          0, sim_params$log_up_slope_pams3[3]) +
    ifelse(param_sizes_all_zero, 0, sim_params$log_up_slope_age[1]) +
    ifelse(param_sizes_all_zero, 0, sim_params$log_up_slope_variant$Wildtype)
  )
  add_to_slope_down <- (ifelse(param_sizes_all_zero, 
                            0, sim_params$log_down_slope_pams3[3]) +
    ifelse(param_sizes_all_zero, 0, sim_params$log_down_slope_age[1]) +
    ifelse(param_sizes_all_zero, 0, sim_params$log_down_slope_variant$Wildtype)
  )
  if (!param_sizes_all_zero) {
    # We are only using T2 for now (and I'm pretty sure this won't change as
    # there is really no time or will for an additional analysis), so we just
    # use the first (and only) position.
    # There is also a contribution from the test_centre but it's so small that
    # it probably doesn't matter.
    assert_that(length(sim_params$intercept_pcr[1]) == 1)
    add_to_intercept <- add_to_intercept + sim_params$intercept_pcr[1]
  }
  stat_log_intercept_mu <- sim_params$log_intercept_mu[1] + add_to_intercept
  stat_log_slope_up_mu <- sim_params$log_up_slope_mu[1] + add_to_slope_up
  stat_log_slope_down_mu <- sim_params$log_down_slope_mu[1] + add_to_slope_down
  stat_alpha_mu <- sim_params$alpha_mu[1]
  # Note: this in on the logarithmic scale
  stat_log_sigma_mu <- ifelse(is.null(sim_data$sigma_mu), 0, sim_data$sigma_mu)

  return(list(
    slope_up = stat_up_slope,
    slope_down = stat_down_slope,
    intercept = stat_intercept,
    sigma = stat_sigma,
    time2peak = time2peak_stat,
    time_from_peak = time_from_peak_stat,
    alpha_mu = stat_alpha_mu,
    slope_up_sigma = exp_func(stat_up_slope_sigma),
    slope_down_sigma = exp_func(stat_down_slope_sigma),
    intercept_sigma = exp_func(stat_intercept_sigma),
    intercept_mu = exp_func(stat_log_intercept_mu),
    slope_up_mu = exp_func(stat_log_slope_up_mu),
    slope_down_mu = exp_func(stat_log_slope_down_mu),
    sigma_mu = exp_func(stat_log_sigma_mu),
    sigma_sigma = exp_func(0.3),
    avg_group_params = avg_group_params_list,
    avg_group_params_diff = avg_group_params_diff_list
  ))
}


#' Return a named list of infection-group-specific/fixed-effects parameter 
#' values used for simulating viral load data.
#'
#' @param stat_params_sim A named list of the true parameter values used for 
#' simulating data.
#' @param exponentiate A logical value specifying whether parameter estimates 
#' should be exponentiated where appropriate (e.g. fixed-effects parameters).
#'
#' @return A named list of infection-group-specific/fixed-effects parameter 
#' values used for simulating data.
calc_sim_group_stats <- function(stat_params_sim, exponentiate = FALSE) {
  if (exponentiate) {
    exp_func <- exp
  } else {
    # Do nothing.
    exp_func <- function(x) {
      return(x)
    }
  }

  intercept_dict <- list()
  slope_up_dict <- list()
  slope_down_dict <- list()

  param_sizes_all_zero <- stat_params_sim$param_sizes_all_zero

  # in our model, PAMS1 is relative to the last category (i.e. not hospitalized, not PAMS)
  intercept_dict[["PAMS"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_intercept_pams3[1] - stat_params_sim$log_intercept_pams3[3]))
  slope_up_dict[["PAMS"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_up_slope_pams3[1] - stat_params_sim$log_up_slope_pams3[3]))
  slope_down_dict[["PAMS"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_down_slope_pams3[1] - stat_params_sim$log_down_slope_pams3[3]))
  # gender (no effect in simulated data)
  intercept_dict[["gender"]] <- exp_func(0)
  slope_up_dict[["gender"]] <- exp_func(0)
  slope_down_dict[["gender"]] <- exp_func(0)
  # in our model, hospitalized is relative to the last category (i.e. not hospitalized, not PAMS)
  intercept_dict[["hospitalized"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_intercept_pams3[2] - stat_params_sim$log_intercept_pams3[3]))
  slope_up_dict[["hospitalized"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_up_slope_pams3[2] - stat_params_sim$log_up_slope_pams3[3]))
  slope_down_dict[["hospitalized"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_down_slope_pams3[2] - stat_params_sim$log_down_slope_pams3[3]))
  intercept_dict[["prior_infection"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_intercept_prior_infection[2]))
  slope_up_dict[["prior_infection"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_up_slope_prior_infection[2]))
  slope_down_dict[["prior_infection"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_down_slope_prior_infection[2]))
  # in our model, age categories 2 and 3 are relative to age category 1, this is not the
  # case in the simulated data, adjust accordingly
  intercept_dict[["age_cat2"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_intercept_age[2] - stat_params_sim$log_intercept_age[1]))
  slope_up_dict[["age_cat2"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_up_slope_age[2] - stat_params_sim$log_up_slope_age[1]))
  slope_down_dict[["age_cat2"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_down_slope_age[2] - stat_params_sim$log_down_slope_age[1]))
  # Age category 3
  intercept_dict[["age_cat3"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_intercept_age[3] - stat_params_sim$log_intercept_age[1]))
  slope_up_dict[["age_cat3"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_up_slope_age[3] - stat_params_sim$log_up_slope_age[1]))
  slope_down_dict[["age_cat3"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_down_slope_age[3] - stat_params_sim$log_down_slope_age[1]))
  # For the SARS-CoV-2 variant, Wildtype is the reference category
  # Alpha
  intercept_dict[["alpha"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_intercept_variant$Alpha - stat_params_sim$log_intercept_variant$Wildtype))
  slope_up_dict[["alpha"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_up_slope_variant$Alpha - stat_params_sim$log_up_slope_variant$Wildtype))
  slope_down_dict[["alpha"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_down_slope_variant$Alpha - stat_params_sim$log_down_slope_variant$Wildtype))
  # Delta
  intercept_dict[["delta"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_intercept_variant$Delta - stat_params_sim$log_intercept_variant$Wildtype))
  slope_up_dict[["delta"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_up_slope_variant$Delta - stat_params_sim$log_up_slope_variant$Wildtype))
  slope_down_dict[["delta"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_down_slope_variant$Delta - stat_params_sim$log_down_slope_variant$Wildtype))
  # Omicron
  intercept_dict[["omicron"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_intercept_variant$Omicron - stat_params_sim$log_intercept_variant$Wildtype))
  slope_up_dict[["omicron"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_up_slope_variant$Omicron - stat_params_sim$log_up_slope_variant$Wildtype))
  slope_down_dict[["omicron"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_down_slope_variant$Omicron - stat_params_sim$log_down_slope_variant$Wildtype))
  # Unknown
  intercept_dict[["unknown"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_intercept_variant$Unknown - stat_params_sim$log_intercept_variant$Wildtype))
  slope_up_dict[["unknown"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_up_slope_variant$Unknown - stat_params_sim$log_up_slope_variant$Wildtype))
  slope_down_dict[["unknown"]] <- exp_func(ifelse(param_sizes_all_zero, 0, stat_params_sim$log_down_slope_variant$Unknown - stat_params_sim$log_down_slope_variant$Wildtype))
  # hospitalized:alpha interaction (no effect in simulated data)
  intercept_dict[["hospitalized:alpha"]] <- exp_func(0)
  slope_up_dict[["hospitalized:alpha"]] <- exp_func(0)
  slope_down_dict[["hospitalized:alpha"]] <- exp_func(0)
  # hospitalized:delta
  intercept_dict[["hospitalized:delta"]] <- exp_func(0)
  slope_up_dict[["hospitalized:delta"]] <- exp_func(0)
  slope_down_dict[["hospitalized:delta"]] <- exp_func(0)
  # hospitalized:omicron
  intercept_dict[["hospitalized:omicron"]] <- exp_func(0)
  slope_up_dict[["hospitalized:omicron"]] <- exp_func(0)
  slope_down_dict[["hospitalized:omicron"]] <- exp_func(0)
  # hospitalized:unknown
  intercept_dict[["hospitalized:unknown"]] <- exp_func(0)
  slope_up_dict[["hospitalized:unknown"]] <- exp_func(0)
  slope_down_dict[["hospitalized:unknown"]] <- exp_func(0)


  return(list(
    "intercept" = intercept_dict,
    "slope_up" = slope_up_dict,
    "slope_down" = slope_down_dict
  ))
}


#### Save (processed) posterior estimates ####

#' Compute and save differences (mean, median or HPDIs) between parameter 
#' estimates and true parameter values, and differences between posterior 
#' predictions (computed from infection-level parameter estimates) and true 
#' values to a common file. This function is only appropriate when simulated 
#' data were used for model fitting.
#'
#' @param output_file A path to the file where posterior estimates and 
#' predictions will be stored.
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param ss A data.table containing summary statistics (including Rhat and ESS 
#' values) from a specific run.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param n_min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
#' @param n_samples A numeric (or NULL) specifying the number of infections 
#' that were simulated (if sim_data not NULL, else NULL).
#' @param ncols A numeric indicating the number of columns intended for the 
#' output file.
#' @param exponentiate A logical value specifying whether to exponentiate 
#' parameter estimates where it is appropriate (e.g. fixed-effects parameters).
write_diff_param_stats <- function(output_file,
                                   draws,
                                   day_data,
                                   sim_data,
                                   model,
                                   model_no_stan,
                                   ss,
                                   stat,
                                   stat2,
                                   n_min_tests,
                                   n_samples,
                                   ncols,
                                   exponentiate) {
  assert_that(!is.null(sim_data))
  con <- file(output_file, open = "a")
  # Note: we are not exponentiating yet, this is done in calc_main_params_diff
  # and calc_group_params_diff
  stat_params_sim <- calc_main_params_summary_stat_sim_data(sim_data,
    day_data = day_data,
    model_no_stan = model_no_stan,
    stat = stat,
    stat2 = stat2,
    exponentiate = FALSE
  )
  sim_params_groups <- calc_sim_group_stats(sim_data$params,
    exponentiate = FALSE
  )

  main_params_diff <- calc_main_params_diff(
    draws = draws,
    day_data = day_data,
    model_no_stan = model_no_stan,
    sim_data = sim_data,
    sim_data_main_params = stat_params_sim,
    stat = stat,
    stat2 = stat2,
    exponentiate = exponentiate
  )
  group_params_diff <- calc_group_params_diff(
    draws = draws,
    day_data = day_data,
    sim_params_groups = sim_params_groups,
    stat = stat,
    exponentiate = exponentiate
  )

  params_vec <- return_params_vec(
    model = model,
    model_no_stan = model_no_stan,
    n_min_tests = n_min_tests,
    n_samples = n_samples, sim = 1,
    stat_params = main_params_diff,
    intercept_dict = group_params_diff$intercept,
    up_dict = group_params_diff$slope_up,
    down_dict = group_params_diff$slope_down,
    earliest_sim_day = min(sim_data$dt$day_sampled),
    latest_sim_day = max(sim_data$dt$day_sampled)
  )

  params_estimated <- paste(params_vec, collapse = "\t")
  assert_that(ncols == length(params_vec))
  writeLines(params_estimated, con)
  close(con)
}


#' Compute and save parameter estimates (mean, median or HPDIs) and posterior 
#' predictions (computed from infection-level parameter estimates) to a common 
#' file. This function saves only estimates, not differences between estimated 
#' and true values (in the case of simulated data).
#'
#' @param output_file A path to the file where posterior estimates and 
#' predictions will be stored.
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param sim_data A named list (or NULL), containing a data.table with 
#' simulated data used for model fitting, and true (simulated) parameter values.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param ss A data.table containing summary statistics (including Rhat and ESS 
#' values) from a specific run.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed for posterior estimates and predictions.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param n_min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
#' @param n_samples A numeric (or NULL) specifying the number of infections 
#' that were simulated (if sim_data not NULL, else NULL).
#' @param ncols A numeric indicating the number of columns intended for the 
#' output file.
#' @param exponentiate A logical value specifying whether to exponentiate 
#' parameter estimates where it is appropriate (e.g. fixed-effects parameters).
write_orig_param_stats <- function(output_file,
                                   draws,
                                   day_data,
                                   sim_data,
                                   model,
                                   model_no_stan,
                                   ss,
                                   stat,
                                   stat2,
                                   n_min_tests,
                                   n_samples,
                                   ncols,
                                   exponentiate = FALSE) {
  con <- file(output_file, open = "a")
  stat_params_estimated <- calc_main_params_summary_stat(draws,
    day_data = day_data,
    stat = stat,
    stat2 = stat2,
    exponentiate = exponentiate
  )

  subgroup_stat_params_estimated <- calc_group_params_summary_stat(
    draws = draws,
    stat = stat,
    exponentiate = exponentiate
  )


  params_estimated_vec <- return_params_vec(
    model = model,
    model_no_stan = model_no_stan,
    n_min_tests = n_min_tests,
    n_samples = n_samples,
    sim = 0,
    stat_params = stat_params_estimated,
    intercept_dict = subgroup_stat_params_estimated$intercept,
    up_dict = subgroup_stat_params_estimated$up,
    down_dict = subgroup_stat_params_estimated$down,
    earliest_sim_day = NA,
    latest_sim_day = NA
  )

  assert_that(ncols == length(params_estimated_vec))
  params_estimated <- paste(params_estimated_vec, collapse = "\t")
  writeLines(params_estimated, con)
  if (!is.null(sim_data)) {
    stat_params_sim <- calc_main_params_summary_stat_sim_data(sim_data,
      day_data = day_data,
      model_no_stan = model_no_stan,
      stat = stat,
      stat2 = stat2,
      exponentiate = exponentiate
    )
    params_sim <- return_sim_params_string(model,
      model_no_stan,
      n_min_tests,
      n_samples = n_samples,
      sim_data = sim_data,
      stat_params_sim = stat_params_sim,
      ncols = ncols, exponentiate = exponentiate
    )
    writeLines(params_sim, con)
  }
  close(con)
}


#' Wrapper to compute and save posterior parameter estimates (mean, median or 
#' HPDIs) and predictions (computed from infection-level parameter estimates) 
#' to a common file.
#'
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param ss A data.table containing summary statistics (including Rhat and 
#' ESS values) from a specific run.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed for posterior estimates and predictions.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#' @param output_file A path to the file where posterior estimates and 
#' predictions will be stored.
#' @param output_file2 A path to the file where posterior estimates of 
#' differences (between true and estimated values in the case of simulated data)
#' and predicted differences will be stored.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param n_min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
#' @param n_samples A numeric (or NULL) specifying the number of infections 
#' that were simulated (if sim_data not NULL, else NULL).
write_param_stats <- function(output_file,
                              stat,
                              stat2,
                              draws,
                              day_data,
                              sim_data,
                              ss,
                              model,
                              model_no_stan,
                              n_min_tests,
                              n_samples,
                              calc_diff = FALSE,
                              exponentiate = FALSE) {
  model_no_suffix <- ifelse(model_no_stan == 173, "_model24_173", "")
  exp_suffix <- ifelse(exponentiate, "_exp", "")
  dir_path <- dirname(output_file)
  filename_base <- tools::file_path_sans_ext(basename(output_file))
  output_file <- file.path(dir_path, paste0(
    filename_base, exp_suffix,
    model_no_suffix, ".tsv"
  ))
  output_file_hdi <- file.path(dir_path, paste0(
    filename_base, exp_suffix,
    model_no_suffix, "_hdi.tsv"
  ))
  output_file_hdi50 <- file.path(dir_path, paste0(
    filename_base, exp_suffix,
    model_no_suffix,
    "_hdi50.tsv"
  ))

  for (output_file_current in c(output_file, 
                                output_file_hdi, 
                                output_file_hdi50)) {
    if (output_file_current == output_file) {
      stat_current <- stat
    } else if (output_file_current == output_file_hdi) {
      stat_current <- "hdi"
    } else {
      stat_current <- "hdi50"
    }
    skip <- FALSE
    # Store the result in a separate file for later evaluation.
    if (!file.exists(output_file_current)) {
      header_params <- return_params_header(exponentiate = exponentiate)
      header <- header_params$header
      ncols <- header_params$ncols
      cat(paste0(header, "\n"), file = output_file_current)
    } else {
      dt_params <- fread(output_file_current, sep = "\t")
      ncols <- ncol(dt_params)
      skip <- model %in% dt_params$model
      remove(dt_params)
    }
    if (skip == FALSE) {
      if (calc_diff) {
        write_diff_param_stats(
          output_file = output_file_current,
          draws = draws,
          day_data = day_data,
          sim_data = sim_data,
          model = model,
          model_no_stan = model_no_stan,
          ss = ss,
          stat = stat_current,
          stat2 = stat2,
          n_min_tests = n_min_tests,
          n_samples = n_samples,
          ncols = ncols,
          exponentiate = exponentiate
        )
      } else {
        write_orig_param_stats(
          output_file = output_file_current,
          draws = draws,
          day_data = day_data,
          sim_data = sim_data,
          model = model,
          model_no_stan = model_no_stan,
          ss = ss,
          stat = stat_current,
          stat2 = stat2,
          n_min_tests = n_min_tests,
          n_samples = n_samples, ncols = ncols,
          exponentiate = exponentiate
        )
      }
    }
  }
}


#' Compute and save conditional posterior predictions computed from 
#' population-level parameter estimates (medians and HPDIs) of main model 
#' parameters for different infection subgroups.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param output_file A string specifying the path to the file where posterior 
#' predictions will be stored.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param reference_variant A string indicating the variant used for prediction.
#' @param reference_hosp_status A string indicating the disease severity 
#' status (PAMS or hospitalized) used for prediction.
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param thin Thinning interval for posterior draws (e.g., `thin = 10` selects 
#' every 10th draw).
write_post_pred_groups <- function(draws,
                                   output_file,
                                   model,
                                   reference_variant,
                                   reference_hosp_status,
                                   model_no_stan,
                                   thin = 1) {
  assert_that(reference_variant %in% c("wildtype", "omicron"))
  assert_that(reference_hosp_status %in% c("hospitalized", 
                                           "PAMS", 
                                           "non-PAMS, non-hosp."))
  dir_path <- dirname(output_file)
  filename_base <- tools::file_path_sans_ext(basename(output_file))
  hosp_suffix <- ifelse(reference_hosp_status == "hospitalized", "_hosp",
    ifelse(reference_hosp_status == "PAMS", "_pams",
      "_non_hosp_non_pams"
    )
  )
  output_file <- file.path(dir_path, paste0(
    filename_base, "_ref_variant_",
    reference_variant,
    hosp_suffix, ".tsv"
  ))

  params <- c(
    "slope_up", "intercept", "slope_down", "time2peak",
    "time_from_peak", "time_from_peak2"
  )
  if (!file.exists(output_file)) {
    header_base <- params
    hdi_header <- rbind(paste0(header_base, "_hdi_lower"), 
                        paste0(header_base, "_hdi_upper"))
    hdi_50_header <- rbind(paste0(header_base, "_hdi50_lower"), 
                           paste0(header_base, "_hdi50_upper"))
    hdi_header_diff <- rbind(paste0(header_base, "_hdi_lower_diff"), 
                             paste0(header_base, "_hdi_upper_diff"))
    hdi_50_header_diff <- rbind(paste0(header_base, "_hdi50_lower_diff"), 
                                paste0(header_base, "_hdi50_upper_diff"))
    header <- c(
      "model", "model_no_stan", "variable", "group",
      rbind(header_base, paste0(header_base, "_diff")), # Like zipping in python
      rbind(hdi_header, hdi_header_diff),
      rbind(hdi_50_header, hdi_50_header_diff)
    )
    ncols <- length(header)
    header <- paste(header, collapse = "\t")
    cat(paste0(header, "\n"), file = output_file)
    skip <- FALSE
  } else {
    dt_out <- fread(output_file, sep = "\t")
    ncols <- ncol(dt_out)
    skip <- model %in% dt_out$model
  }
  mapping <- SUBGRP_CATS_LABEL_MAPPING2
  grp_vars <- SUBGRPS_MAIN

  if (skip == FALSE) {
    con <- file(output_file, open = "a")
    for (group_var in grp_vars) {
      draws_dt <- draws_by_grp(
        draws = draws,
        grp_var = group_var,
        thin = thin,
        calc_times = TRUE,
        calc_bc = FALSE,
        reference_variant = reference_variant,
        reference_hosp_status = reference_hosp_status
      )
      reference_label <- mapping[[group_var]][1]
      counter <- 1
      for (group in mapping[[group_var]]) {
        line <- c(model, model_no_stan, group_var, group)
        for (summary_stat in c("median", "hdi", "hdi50")) {
          for (param in params) {
            reference_vals <- draws_dt[get(group_var) == eval(reference_label)][[param]]
            group_vals <- draws_dt[get(group_var) == eval(group)][[param]]

            grp_ss <- calc_summary_stat(summary_stat = summary_stat, group_vals)
            grp_ss_diff <- calc_difference(
              param_true = reference_vals,
              param_est = group_vals,
              summary_stat = summary_stat
            )
            line <- c(line, grp_ss, grp_ss_diff)
          }
        }
        assert_that(length(line) == ncols)
        counter <- counter + 1
        writeLines(paste(line, collapse = "\t"), con = con)
      }
    }
    close(con)
  }
}


#' Compute and save posterior predictions computed from infection-level 
#' parameter estimates for different infection subgroups.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed for posterior predictions.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param output_file A path to the file where posterior predictions will be 
#' stored.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param n_min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
write_param_stats_post_pred_inf <- function(output_file,
                                            stat,
                                            stat2,
                                            draws,
                                            day_data,
                                            model,
                                            model_no_stan,
                                            n_min_tests) {
  suffix <- "_post_pred_inf"
  dir_path <- dirname(output_file)
  filename_base <- tools::file_path_sans_ext(basename(output_file))
  output_file <- file.path(dir_path, 
                           paste0(filename_base, suffix, ".tsv"))
  output_file_hdi <- file.path(dir_path, 
                               paste0(filename_base, suffix, "_hdi.tsv"))
  output_file_hdi50 <- file.path(dir_path, 
                                 paste0(filename_base, suffix, "_hdi50.tsv"))

  for (output_file_current in c(output_file, 
                                output_file_hdi, 
                                output_file_hdi50)) {
    if (output_file_current == output_file) {
      stat_current <- stat
    } else if (output_file_current == output_file_hdi) {
      stat_current <- "hdi"
    } else {
      stat_current <- "hdi50"
    }
    skip <- FALSE
    # Store the result in a separate file for later evaluation.
    if (!file.exists(output_file_current)) {
      header_params <- return_params_header(
        exponentiate = FALSE,
        post_pred_groups_inf_only = TRUE
      )
      header <- header_params$header
      ncols <- header_params$ncols
      cat(paste0(header, "\n"), file = output_file_current)
    } else {
      dt_params <- fread(output_file_current, sep = "\t")
      ncols <- ncol(dt_params)
      skip <- model %in% dt_params$model
      remove(dt_params)
    }
    if (skip == FALSE) {
      con <- file(output_file_current, open = "a")
      assert_that(stat_current %in% c("median", "hdi", "hdi50", "mean"))
      if (stat_current == "median") {
        stat_func <- median
      } else if (stat_current == "mean") {
        stat_func <- mean
      } else if (stat_current == "hdi50") {
        stat_func <- calc_50_hdi
      } else {
        stat_func <- calc_94_hdi
      }

      by_draw_ID <-
        draws_by_id(draws,
          day_data,
          c("slope_up", 
            "slope_down", 
            "intercept", 
            "sigma", 
            "time2peak", 
            "time_from_peak"),
          thin = 1
        )

      avg_group_params_list <- return_group_params_avg(
        dt = by_draw_ID,
        stat_func = stat_func,
        stat2 = stat2,
        day_data = day_data,
        sim_data = NULL,
        diff = FALSE,
        PAMS3 = TRUE
      )
      avg_group_params_diff_list <- return_group_params_avg_diff(
        dt = by_draw_ID,
        stat_func = stat_func,
        stat2 = stat2,
        day_data = day_data,
        sim_data = NULL,
        diff = FALSE,
        PAMS3 = TRUE
      )

      stat_params <- list(
        avg_group_params = avg_group_params_list,
        avg_group_params_diff = avg_group_params_diff_list
      )

      group_params <- process_avg_group_params_list(stat_params$avg_group_params, 
                                                    PAMS3 = TRUE)
      group_params_diff <- process_avg_group_params_diff_list(stat_params$avg_group_params_diff,
        PAMS3 = TRUE
      )
      params_vec <- c(
        model,
        model_no_stan,
        n_min_tests,
        group_params,
        group_params_diff
      )
      params_estimated <- paste(params_vec, collapse = "\t")
      assert_that(ncols == length(params_vec))
      writeLines(params_estimated, con)
      close(con)
    }
  }
}


#' Second wrapper to compute and save posterior parameter estimates (mean, 
#' median or HPDIs) and predictions (computed from infection-level parameter 
#' estimates) to a common file.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param ss A data.table containing summary statistics (including Rhat and 
#' ESS values) from a specific run.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed for posterior estimates and predictions.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#' @param output_file A path to the file where posterior estimates and 
#' predictions will be stored.
#' @param output_file2 A path to the file where posterior estimates of 
#' differences (between true and estimated values in the case of simulated data)
#' and predicted differences will be stored.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param n_min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
#' @param n_samples A numeric (or NULL) specifying the number of infections 
#' that were simulated (if sim_data not NULL, else NULL).
#' @param post_pred_groups A logical value indicating whether to save 
#' conditional posterior predictions computed from population-level
#' parameter estimates instead of posterior estimates. I think this option 
#' should be removed eventually as posterior predictions are already computed
#' and saved in make_post_pred_by_group (to file 
#' "params_stats_post_pred_groups.tsv"). We are only using the values in
#' "params_stats_post_pred_groups.tsv" anyway.
save_param_stats <- function(draws,
                             day_data,
                             ss,
                             stat,
                             stat2,
                             sim_data,
                             output_file,
                             output_file2,
                             model,
                             model_no_stan,
                             n_min_tests,
                             n_samples) {
  for (exponentiate in c(FALSE, TRUE)) {
    if (!is.null(output_file)) {
      write_param_stats(
        output_file = output_file,
        stat = stat,
        stat2 = stat2,
        draws = draws,
        day_data = day_data,
        sim_data = sim_data,
        model = model,
        model_no_stan = model_no_stan,
        ss = ss,
        n_min_tests = n_min_tests,
        n_samples = n_samples,
        exponentiate = exponentiate
      )
    }
    if (!is.null(output_file2)) {
      # Write files specifying differences between true and estimated parameter
      # values.
      assert_that(!is.null(sim_data))
      write_param_stats(
        output_file = output_file2,
        stat = stat,
        stat2 = stat2,
        draws = draws,
        day_data = day_data,
        sim_data = sim_data,
        ss = ss,
        model_no_stan = model_no_stan,
        model = model,
        n_min_tests = n_min_tests,
        n_samples = n_samples,
        calc_diff = TRUE,
        exponentiate = exponentiate
      )
    }
  }
}


#' Wrapper function to compute and save posterior predictions computed from 
#' infection-level parameter estimates for different infection subgroups.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed for posterior predictions.
#' @param stat2 The summary statistic (mean or median) that will be used for 
#' averaging over infections within subgroups first, to then summarize over 
#' draws.
#' @param output_file A path to the file where posterior predictions will be 
#' stored.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param n_min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
save_param_stats_post_pred_inf <- function(draws,
                                           day_data,
                                           stat,
                                           stat2,
                                           output_file,
                                           model,
                                           model_no_stan,
                                           n_min_tests) {
  if (!is.null(output_file)) {
    write_param_stats_post_pred_inf(
      output_file = output_file,
      stat = stat,
      stat2 = stat2,
      draws = draws,
      day_data = day_data,
      model = model,
      model_no_stan = model_no_stan,
      n_min_tests = n_min_tests
    )
  }
}


#' Save Rhat values to a common file.
#'
#' @param model A string specifying the suffix of the model (includes numeric 
#' model identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param ss A data.table containing summary statistics (including Rhat and 
#' ESS values) from a specific run.
#' @param output_file A string specifying the path to the file where Rhat 
#' values will be stored.
#' @param min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
#' @param simulated A logical value indicating whether the run was performed on 
#' simulated data.
collect_rhat_stats <- function(model,
                               model_no_stan,
                               ss,
                               output_file,
                               min_tests,
                               simulated) {
  mapping <- SUBGRP_MAPPING

  params <- c(
    "alpha_mu", "log_intercept_mu", "log_slope_up_mu",
    "log_slope_down_mu", "sigma_mu", "sigma_sigma", "intercept_sigma",
    "slope_up_sigma", "slope_down_sigma"
  )
  params_colnames <- c(
    "alpha_mu", "log_intercept_mu", "log_slope_up_mu",
    "log_slope_down_mu", "log_sigma_mu", "log_sigma_sigma",
    "log_intercept_sigma", "log_slope_up_sigma",
    "log_slope_down_sigma"
  )
  assert_that(length(params) == length(params_colnames))

  intercept_group_param <- paste0("betaPGH_intercept[", names(mapping), "]")
  slope_up_group_param <- paste0("betaPGH_slope_up[", names(mapping), "]")
  slope_down_group_param <- paste0("betaPGH_slope_down[", names(mapping), "]")

  intercept_group_param_col <- paste0("intercept_", mapping)
  slope_up_group_param_col <- paste0("slope_up_", mapping)
  slope_down_group_param_col <- paste0("slope_down_", mapping)

  params <- c(
    params, intercept_group_param, slope_up_group_param,
    slope_down_group_param
  )
  params_colnames <- c(
    params_colnames, intercept_group_param_col,
    slope_up_group_param_col, slope_down_group_param_col
  )

  if (!file.exists(output_file)) {
    header <- paste(c(c(
      "model", "model_no_stan", "min_tests",
      "simulated"
    ), params_colnames), collapse = "\t")
    cat(paste0(header, "\n"), file = output_file)
    skip <- FALSE
  } else {
    dt_params <- fread(output_file, sep = "\t")
    skip <- model %in% dt_params$model
    remove(dt_params)
  }
  if (skip == FALSE) {
    con <- file(output_file, open = "a")
    params_estimated_1 <- c(model, model_no_stan, min_tests, simulated)
    params_estimated_2 <- c()
    for (p in params) {
      params_estimated_2 <- c(params_estimated_2, ss[variable == p]$rhat)
    }
    params_estimated <- paste(c(params_estimated_1, params_estimated_2),
      collapse = "\t"
    )
    writeLines(params_estimated, con)
    close(con)
  }
}
 