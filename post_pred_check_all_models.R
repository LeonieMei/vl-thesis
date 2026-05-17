library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(optparse)
library(here)

TOP <- here()
source(file.path(TOP, "data_params.R"))
source(file.path(TOP, "post_analysis_plot.R"))

DIR_FIGURES_POST_PRED_CHECK <- file.path(DIR_FIGURES_VL, "post_pred_checks")

option_list <- list(
  make_option(c("--summary_stat"),
    default = "median",
    help = "The summary statistic (median or mean) to use for summarizing the main model parameters (intercept, up-slope, down-slope)."
  ),
  make_option(c("--sampling_iterations"),
    default = 1000,
    help = "The number of MCMC sampling iterations run."
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

assert_that(opt$summary_stat %in% c("median", "mean"), msg = "You can only use median or mean as a summary statistic.")
assert_that(opt$sampling_iterations %in% c(1000, 2000))


#' Make plots of posterior estimates of log10 viral load trajectories and 
#' temporal placement of data points.
#'
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param VLCP_by_day_ID A data.table with infection-wise posterior predictions 
#' of median viral load trajectories.
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param shift_draws A data.table with posterior draws of temporal shifts, 
#' grouped by infection.
#' @param shifted_data_by_draw_day_ID A data.table with posterior draws of 
#' shifted days (i.e. estimated days of sampling relative to the day of peak 
#' viral load), and imputed log10 viral loads grouped by infection.
#' @param day_data A data.table with PCR test data.
#' @param dir_pdata The name of the directory where posterior predictions are 
#' stored.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param day_var A string specifying the day variable used in the model 
#' (e.g. "day" with day 0 being the first day with a log10 viral
#' @param testing A logical value. If `TRUE`, loads only a subset of posterior 
#' draws (for testing/debugging).
#' @param color_by A column name to color data points by (e.g. imputed viral 
#' loads).
#' @param plot_lines A logical value indicating whether to plot shifted data 
#' points or estimated infection-wise viral load trajectories.
#' @param plot_n_infections Randomly select n infections to make posterior 
#' prediction plots for.
#' @param stat The summary statistic (mean or median) applied to posterior 
#' estimates and predictions.
#' @param hide_x_label_p1 A logical value indicating whether to hide the x-axis 
#' label for the upper plot.
#'
#' @return A named list of ggplot objects, the first one showing estimated 
#' log10 viral load trajectories, the second one showing shifted data points 
#' (and optionally trajectories).
make_plots <- function(VLCP_by_day_draws,
                       VLCP_by_day_ID,
                       draws,
                       shift_draws,
                       shifted_data_by_draw_day_ID,
                       day_data,
                       dir_pdata,
                       model,
                       model_no_stan,
                       day_var,
                       testing = FALSE,
                       color_by = "imputed",
                       plot_lines = FALSE,
                       plot_n_infections = NULL,
                       stat = "median",
                       hide_x_label_p1 = FALSE) {
  assert_that(stat %in% c("mean", "median"))

  if (testing == TRUE) {
    tc_days <- c(-7:4, seq(5, 25, 2.5), 40)
    thin <- 10
  } else {
    tc_days <- c(
      seq(-11, -7, by = 1), seq(-6.5, -3, by = .5),
      seq(-2.75, 0, by = .25), seq(0.5, 5, by = .5),
      seq(6, 25, 1), 27, 30, 35, 40
    )
    thin <- 4
  }
  vl_by_day_filepath <- file.path(
    dir_pdata,
    paste0(
      "VL_by_day_", model,
      ifelse(stat == "median", "_median", ""),
      ".Rdata"
    )
  )
  if (file.exists(vl_by_day_filepath)) {
    load(vl_by_day_filepath)
  } else {
    VL_by_key_days_draws <-
      VLCP_by_day_draws %>%
      .[day_shifted %in% c(0, 5, 8)] %>%
      .[, list(mean = mean(log10_load), sd = sd(log10_load), 
               median = median(log10_load)),
        by = c(".draw", "day_shifted")
      ] %>%
      melt(id.vars = c(".draw", "day_shifted"), variable.name = "statistic")
    if (stat == "median") {
      VL_by_day <-
        VLCP_by_day_draws %>%
        .[, value := collapse::fmedian(log10_load), 
          by = .(day_shifted, .draw)] %>%
        get_stats(by = "day_shifted") %>%
        setnames("median", "log10_load") %>%
        .[, ID := 0]
      save(VL_by_day, VL_by_key_days_draws, file = vl_by_day_filepath)
    } else {
      VL_by_day <-
        VLCP_by_day_draws %>%
        .[, value := collapse::fmean(log10_load), 
          by = .(day_shifted, .draw)] %>%
        get_stats(by = "day_shifted") %>%
        setnames("mean", "log10_load") %>%
        .[, ID := 0]
      save(VL_by_day, VL_by_key_days_draws, file = vl_by_day_filepath)
    }
  }

  xlim <- c(-10, 40)
  ylim <- c(-2.5, 14)


  if (is.null(plot_n_infections)) {
    random_IDs <- unique(day_data$ID)
  } else {
    random_IDs <- sample(unique(day_data$ID),
      size = plot_n_infections
    )
  }

  p1 <- modified_plot(
    plot_by_day(VL_by_day,
      VLCP_by_day_ID %>% .[ID %in% random_IDs] %>% as.data.frame(),
      y.var = "log10_load",
      ylim = ylim, xlim = xlim,
      alpha = 0.05
    ) +
      geom_hline(yintercept = 0, lty = 3),
    hide_x_label = hide_x_label_p1
  )

  if (plot_lines) {
    p2 <- plot_shifted_data_lines(
      day_data = day_data,
      draws = draws,
      VL_by_day = VL_by_day,
      shift_draws = shift_draws,
      shifted_data_by_draw_day_ID = shifted_data_by_draw_day_ID,
      day_var = day_var,
      random_IDs = random_IDs,
      color_by = color_by,
      xlim = xlim,
      ylim = ylim,
      stat = stat,
      thin = 1
    )
  } else {
    p2 <- plot_shifted_data_points(
      day_data = day_data,
      draws = draws,
      VL_by_day = VL_by_day,
      shift_draws = shift_draws,
      shifted_data_by_draw_day_ID = shifted_data_by_draw_day_ID,
      day_var = day_var,
      random_IDs = random_IDs,
      color_by = color_by,
      xlim = xlim,
      ylim = ylim,
      stat = stat,
      thin = 1,
      alpha = 0.05
    )
    p2 <- modified_plot(p2, hide_legend = TRUE)
  }
  return(list(p1 = p1, p2 = p2))
}


#' Generate plots of posterior predictive checks.
#'
#' @param opt A named list of command-line arguments.
#' @param models A vector of strings specifying the model names (includes 
#' numeric model identifiers, data filtering criteria, and run settings).
#' @param model_dirs Paths to the directories containing model fitting results.
#' @param model_stan_nos The numbers/identifiers of the stan files used for 
#' fitting.
#' @param day_vars A vector of strings specifying the day variables used in the 
#' models (e.g. "day" with day 0 being the first day with a log10 viral
#'  load > 3; or "day2" with day 0 being the day of the highest viral load 
#'  measured in the infection).
#' @param dir_figures The directory to save figures to.
make_post_pred_check_plots <- function(opt,
                                       models,
                                       model_dirs,
                                       model_stan_nos,
                                       day_vars,
                                       dir_figures) {
  width <- VL_PLOT_WIDTH_THESIS
  height <- VL_PLOT_HEIGHT_THESIS * 2
  sampling_iteration_suffix <- ifelse(grep("iter2000", models[1]),
    "_iter2000", ""
  )

  counter_model <- 1
  for (model_list in models) {
    plot_list <- list()
    counter_tests <- 1
    for (model in model_list) {
      run_data <- load_draws(
        model = model,
        model_dir = model_dirs[counter_model],
        model_no_stan = model_stan_nos[counter_model],
        inc_warmup = FALSE, testing = FALSE,
        test_thin = 1
      )
      dir_pdata <- file.path(DIR_PDATA, 
                             paste0("model", model_stan_nos[counter_model]))
      pp_results <- make_post_pred_all(
        draws = run_data$draws,
        day_data = run_data$day_data,
        model = model,
        dir_pdata = dir_pdata,
        day_var = day_vars[counter_model],
        model_no_stan = model_stan_nos[counter_model],
        ss = run_data$ss,
        stat = opt$summary_stat,
        testing = FALSE
      )
      ps <- make_plots(
        VLCP_by_day_draws = pp_results$VLCP_by_day_draws,
        VLCP_by_day_ID = pp_results$VLCP_by_day_ID,
        shift_draws = pp_results$shift_draws,
        shifted_data_by_draw_day_ID = pp_results$shifted_data_by_draw_day_ID,
        draws = run_data$draws,
        day_data = run_data$day_data,
        dir_pdata = dir_pdata,
        model = model,
        model_no_stan = model_stan_nos[counter_model],
        day_var = day_vars[counter_model],
        testing = FALSE,
        color_by = "is.difficult",
        plot_lines = FALSE,
        plot_n_infections = NULL,
        stat = "median"
      )


      if (counter_tests == 1) {
        p1 <- ps$p1 + theme(axis.title.x = element_blank())
        p2 <- ps$p2 + theme(axis.title.x = element_blank(), 
                            axis.title.y = element_blank())
      } else {
        p1 <- ps$p1 + labs(x = DAYS_LABEL)
        p2 <- ps$p2 + labs(x = DAYS_LABEL) + theme(axis.title.y = element_blank())
      }
      plot_list <- c(plot_list, p1, p2)
      counter_tests <- counter_tests + 1

      rm(pp_results)
      rm(run_data)
      rm(ps)
    }
    p_all <- arrangeGrob(plot_annotation(plot_list[[1]], letter = "A"),
      plot_annotation(plot_list[[2]], letter = "B"),
      plot_annotation(plot_list[[3]], letter = "C"),
      plot_annotation(plot_list[[4]], letter = "D"),
      ncol = 2
    )
    filename_suffix <- paste0("_model", model_stan_nos[counter_model])
    ggsave(p_all,
      width = width, height = height,
      units = "in",
      filename = file.path(dir_figures, paste0(
        "TC",
        filename_suffix,
        "_post_pred_check",
        sampling_iteration_suffix,
        ".png"
      )),
      dpi = 300, bg = "white"
    )
    rm(p_all)
    counter_model <- counter_model + 1
  }
}


main <- function(opt) {
  models <- list(
    list(
      paste0("model1_threading_sel3_chains4_n_iter", opt$sampling_iterations, 
             "_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5"),
      paste0("model1_threading_sel6_chains4_n_iter", opt$sampling_iterations, 
             "_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5")
    ),
    list(
      paste0("model2_threading_sel3_chains4_n_iter", opt$sampling_iterations, 
             "_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5"),
      paste0("model2_threading_sel6_chains4_n_iter", opt$sampling_iterations, 
             "_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5")
    )
  )
  day_vars <- c("day", "day2")
  model_stan_nos <- c(1, 2)
  model_dirs <- c("model1", "model2")

  make_post_pred_check_plots(opt,
    models = models,
    model_dirs = model_dirs,
    model_stan_nos = model_stan_nos,
    day_vars = day_vars,
    dir_figures = DIR_FIGURES_POST_PRED_CHECK
  )
}

main(opt)
