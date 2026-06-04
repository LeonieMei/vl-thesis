library(ggplot2)
library(optparse)
library(here)

TOP <- here()
source(file.path(TOP, "sars2vl", "data_params.R"))
source(file.path(TOP, "sars2vl", "plot_params.R"))
source(file.path(TOP, "sars2vl", "post_analysis_plot.R"))
DIR_FIGURES_TRAJECTORIES <- file.path(DIR_FIGURES_VL, "all_models")

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


#' Create and save conditional posterior prediction plots (for disease 
#' severity, age, variant and prior infection status). Make two plots for each
#' of these variables, one with predictions computed from infection-level 
#' parameter estimates, one computed from population-level parameter estimates.
#'
#' @param opt A named list of command-line arguments.
#' @param models A vector of strings specifying the model names (includes 
#' numeric model identifiers, data filtering criteria, and run settings).
#' @param day_vars A vector of strings specifying the day variables used in the 
#' models (e.g. "day" with day 0 being the first day with a log10 viral
#'  load > 3; or "day2" with day 0 being the day of the highest viral load 
#'  measured in the infection).
#' @param model_nos_stan The numbers/identifiers of the stan files used for 
#' fitting.
#' @param model_dirs Paths to the directories containing model fitting results.
#' @param dir_figures A string specifying the directory path where figures will 
#' be saved.
#' @param width The figure width (in inches).
#' @param height The figure height (in inches).
create_multi_model_plots <- function(opt,
                                     models,
                                     day_vars,
                                     model_nos_stan,
                                     model_dirs,
                                     dir_figures,
                                     width = VL_PLOT_WIDTH_THESIS,
                                     height = VL_PLOT_HEIGHT_THESIS) {
  height_factor <- 1.333333

  post_pred_pop_funcs <- rev(c(
    make_pams_plot_post_pred_pop,
    make_age_plot_post_pred_pop,
    make_variant_plot_post_pred_pop,
    make_prior_infection_plot_post_pred_pop
  ))
  post_pred_inf_funcs <- rev(c(
    make_pams_plot,
    make_age_plot,
    make_variant_plot,
    make_prior_infection_plot
  ))

  annotations <- c("A", "B", "C", "D")
  legend_titles <- rev(list("", "Age (years)", "", "Prior infection"))
  legend_positions_inf <- rev(c(list(
    list(c(.45, .95), c(.45, .95), NULL, NULL),
    list(c(.75, .8), c(.75, .8), NULL, NULL),
    list(c(.7, .85), c(.7, .85), NULL, NULL),
    list(c(.7, .8), c(.7, .8), NULL, NULL)
  )))
  legend_positions_pop <- rev(c(list(
    list(c(.6, .95), c(.6, .95), NULL, NULL),
    list(c(.75, .8), c(.75, .8), NULL, NULL),
    list(c(.7, .85), c(.7, .85), NULL, NULL),
    list(c(.7, .8), c(.7, .8), NULL, NULL)
  )))

  hide_x_labels <- c(TRUE, TRUE, FALSE, FALSE)
  hide_y_labels <- c(FALSE, TRUE, FALSE, TRUE)
  hide_legends_inf <- c(FALSE, FALSE, TRUE, TRUE)
  hide_legends_pop <- c(FALSE, TRUE, TRUE, TRUE)
  filename_suffixes <- rev(c("PAMS3", "age", "variant", "prior_infection"))
  filename_suffixes2 <- c("_hosp", "_pams")
  filename_suffixes3 <- c("_omicron_pop", "_inf")
  sampling_iteration_suffix <- ifelse(grep("iter2000", models[1]),
    "_iter2000", ""
  )


  # Posterior predictions for subgroups, NOT based on observed data
  counter_hosp_status <- 1
  for (hosp_status_current in c("hospitalized", "PAMS")) {
    filename_suffix2 <- filename_suffixes2[counter_hosp_status]
    counter_pred_type <- 1
    for (plot_func_list in c(list(post_pred_pop_funcs), 
                             list(post_pred_inf_funcs))) {
      filename_suffix3 <- filename_suffixes3[counter_pred_type]
      if (counter_pred_type == 2) {
        legend_positions <- legend_positions_inf
        hide_legends <- hide_legends_inf
        if (counter_hosp_status == 1) {
          filename_suffix2 <- ""
        } else {
          # We don't need to run prediction twice when using infection-level 
          # parameters.
          break
        }
      } else {
        legend_positions <- legend_positions_pop
        hide_legends <- hide_legends_pop
      }

      counter_grp <- 1
      for (plot_func in plot_func_list) {
        legend_title <- legend_titles[counter_grp]
        filename_suffix <- filename_suffixes[counter_grp]
        counter_model <- 1
        plot_list <- list()
        for (model in models) {
          run_data <- load_draws(
            model = model,
            model_dir = model_dirs[counter_model],
            model_no_stan = model_nos_stan[counter_model],
            inc_warmup = FALSE, testing = FALSE,
            test_thin = 1
          )
          dir_pdata <- file.path(DIR_PDATA, paste0("model", model_nos_stan[counter_model]))

          # Predictions based on infection-level parameter estimates.
          if (counter_pred_type == 2) {
            pp_results <- make_post_pred_all(
              draws = run_data$draws,
              day_data = run_data$day_data,
              model = model,
              dir_pdata = dir_pdata,
              day_var = day_vars[counter_model],
              model_no_stan = model_nos_stan[counter_model],
              ss = run_data$ss,
              stat = opt$summary_stat,
              testing = FALSE
            )

            p <- plot_func(
              VLCP_by_day_draws = pp_results$VLCP_by_day_draws,
              day_data = run_data$day_data,
              legend_title = legend_title,
              legend_position = legend_positions[[counter_grp]][[counter_model]],
              hide_x_label = hide_x_labels[counter_model],
              hide_y_label = hide_y_labels[counter_model],
              hide_legend = hide_legends[counter_model],
              add_to_font_size = -2
            )
          } else {
            # Predictions based on population-level parameter estimates.
            p <- plot_func(
              draws = run_data$draws,
              day_data = run_data$day_data,
              legend_title = legend_title,
              legend_position = legend_positions[[counter_grp]][[counter_model]],
              hide_x_label = hide_x_labels[counter_model],
              hide_y_label = hide_y_labels[counter_model],
              hide_legend = hide_legends[counter_model],
              add_to_font_size = -2,
              reference_hosp_status = hosp_status_current,
              reference_variant = "omicron"
            )
          }
          if (!hide_legends[counter_model]) {
            p <- p + shrink_legend()
          }
          plot_list <- c(plot_list, p)

          ggsave(p,
            width = width * 2 / 3, height = height * height_factor * 2 / 3,
            units = "in",
            filename = file.path(dir_figures, paste0(
              "TC_by_",
              filename_suffix,
              filename_suffix2,
              filename_suffix3,
              "_",
              model,
              sampling_iteration_suffix,
              ".png"
            )),
            dpi = 300, bg = "white"
          )
          counter_model <- counter_model + 1
        }
        p_all <- arrangeGrob(plot_annotation(plot_list[[1]], letter = "A"),
          plot_annotation(plot_list[[2]], letter = "B"),
          plot_annotation(plot_list[[3]], letter = "C"),
          plot_annotation(plot_list[[4]], letter = "D"),
          ncol = 2
        )


        ggsave(p_all,
          width = width, height = height * height_factor,
          units = "in",
          filename = file.path(dir_figures, paste0(
            "TC_by_",
            filename_suffix,
            filename_suffix2,
            filename_suffix3,
            sampling_iteration_suffix,
            ".png"
          )),
          dpi = 300, bg = "white"
        )
        rm(p_all)
        counter_grp <- counter_grp + 1
      }
      counter_pred_type <- counter_pred_type + 1
    }
    counter_hosp_status <- counter_hosp_status + 1
  }
}


#' Create and save conditional posterior prediction plots for different 
#' SARS-CoV-2 variants. Make two plots with predictions, one computed with 
#' infection-level parameter estimates, one computed with population-level 
#' parameter estimates.
#'
#' @param opt A named list of command-line arguments.
#' @param models A vector of strings specifying the model names (includes 
#' numeric model identifiers, data filtering criteria, and run settings).
#' @param day_vars A vector of strings specifying the day variables used in the 
#' models (e.g. "day" with day 0 being the first day with a log10 viral 
#' load > 3; or "day2" with day 0 being the day of the highest viral load 
#' measured in the infection).
#' @param model_nos_stan The numbers/identifiers of the stan files used for 
#' fitting.
#' @param model_dirs Paths to the directories containing model fitting results.
#' @param dir_figures A string specifying the directory path where figures will 
#' be saved.
#' @param width The figure width (in inches).
#' @param height The figure height (in inches).
create_variant_plots <- function(opt,
                                 models,
                                 day_vars,
                                 model_nos_stan,
                                 model_dirs,
                                 dir_figures,
                                 width = VL_PLOT_WIDTH_THESIS,
                                 height = VL_PLOT_HEIGHT_THESIS) {
  height_factor <- 1.333333

  plot_funcs <- c(make_variant_plot_post_pred_pop, make_variant_plot)

  annotations <- c("A", "B", "C", "D")
  legend_title <- ""
  legend_position <- c(list(c(.7, .85), c(.7, .85), NULL, NULL))

  hide_x_labels <- c(TRUE, TRUE, FALSE, FALSE)
  hide_y_labels <- c(FALSE, TRUE, FALSE, TRUE)
  hide_legends <- c(FALSE, FALSE, TRUE, TRUE)
  filename_suffix <- "variant"
  filename_suffixes2 <- c("_pop", "_inf")
  filename_suffixes3 <- c("_hosp", "_pams")
  sampling_iteration_suffix <- ifelse(grep("iter2000", models[1]),
    "_iter2000", ""
  )

  # Posterior predictions for subgroups, NOT based on observed data
  counter_hosp_status <- 1
  for (hosp_status_current in c("hospitalized", "PAMS")) {
    filename_suffix3 <- filename_suffixes3[counter_hosp_status]
    counter_pred_type <- 1
    filename_suffix2 <- filename_suffixes2[counter_pred_type]
    for (plot_func in plot_funcs) {
      if (counter_pred_type == 2) {
        if (counter_hosp_status == 1) {
          filename_suffix3 <- ""
        } else {
          # We don't need to run prediction twice when using infection-level 
          # parameters.
          break
        }
      }
      counter_model <- 1
      plot_list <- list()
      for (model in models) {
        run_data <- load_draws(
          model = model,
          model_dir = model_dirs[counter_model],
          model_no_stan = model_nos_stan[counter_model],
          inc_warmup = FALSE,
          testing = FALSE,
          test_thin = 1
        )
        dir_pdata <- file.path(DIR_PDATA, paste0("model", 
                                                 model_nos_stan[counter_model]))

        # Make plots based on grouping variables.
        # Predictions computed with infection-level parameter estimates.
        if (counter_pred_type == 2) {
          pp_results <- make_post_pred_all(
            draws = run_data$draws,
            day_data = run_data$day_data,
            model = model,
            dir_pdata = dir_pdata,
            day_var = day_vars[counter_model],
            model_no_stan = model_nos_stan[counter_model],
            ss = run_data$ss,
            stat = opt$summary_stat,
            testing = FALSE
          )

          p <- plot_func(
            VLCP_by_day_draws = pp_results$VLCP_by_day_draws,
            day_data = run_data$day_data,
            legend_title = legend_title,
            legend_position = legend_position[[counter_model]],
            hide_x_label = hide_x_labels[counter_model],
            hide_y_label = hide_y_labels[counter_model],
            hide_legend = hide_legends[counter_model],
            add_to_font_size = -2
          )
        } else {
          # Predictions computed with population-level parameter estimates.
          p <- plot_func(
            draws = run_data$draws,
            day_data = run_data$day_data,
            legend_title = legend_title,
            legend_position = legend_position[[counter_model]],
            hide_x_label = hide_x_labels[counter_model],
            hide_y_label = hide_y_labels[counter_model],
            hide_legend = hide_legends[counter_model],
            add_to_font_size = -2,
            reference_hosp_status = hosp_status_current,
            reference_variant = "omicron"
          )
        }
        if (!hide_legends[counter_model]) {
          p <- p + shrink_legend()
        }
        plot_list <- c(plot_list, p)
        counter_model <- counter_model + 1
      }

      p_all <- arrangeGrob(plot_annotation(plot_list[[1]], letter = "A"),
        plot_annotation(plot_list[[2]], letter = "B"),
        plot_annotation(plot_list[[3]], letter = "C"),
        plot_annotation(plot_list[[4]], letter = "D"),
        ncol = 2
      )


      ggsave(p_all,
        width = width, height = height * height_factor,
        units = "in",
        filename = file.path(dir_figures, paste0(
          "TC_by_",
          filename_suffix,
          filename_suffix2,
          filename_suffix3,
          sampling_iteration_suffix,
          "_omicron.png"
        )),
        dpi = 300, bg = "white"
      )
      rm(p_all)
      counter_pred_type <- counter_pred_type + 1
    }
    counter_hosp_status <- counter_hosp_status + 1
  }
}


#' Create and save histograms to visualize the distribution of estimated 
#' sampling days.
#'
#' @param opt A named list of command-line arguments.
#' @param models A vector of strings specifying the model names (includes 
#' numeric model identifiers, data filtering criteria, and run settings).
#' @param day_vars A vector of strings specifying the day variables used in the 
#' models (e.g. "day" with day 0 being the first day with a log10 viral
#' load > 3; or "day2" with day 0 being the day of the highest viral load 
#' measured in the infection).
#' @param model_nos_stan The numbers/identifiers of the stan files used for 
#' fitting.
#' @param model_dirs Paths to the directories containing model fitting results.
#' @param dir_figures A string specifying the directory path where figures will 
#' be saved.
#' @param width The figure width (in inches).
#' @param height The figure height (in inches).
create_day_shifted_plots <- function(opt,
                                     models,
                                     day_vars,
                                     model_nos_stan,
                                     model_dirs,
                                     dir_figures,
                                     width = VL_PLOT_WIDTH_THESIS,
                                     height = VL_PLOT_HEIGHT_THESIS) {
  height_factor <- 1.333333

  plot_list <- list()
  annotations <- c("A", "B", "C", "D")
  hide_x_labels <- c(TRUE, TRUE, FALSE, FALSE)

  sampling_iteration_suffix <- ifelse(grep("iter2000", models[1]),
    "_iter2000", ""
  )

  # Plot distributions of estimated sampling days for all models
  counter <- 1
  for (model in models) {
    run_data <- load_draws(
      model = model,
      model_dir = model_dirs[counter],
      model_no_stan = model_nos_stan[counter],
      inc_warmup = FALSE,
      testing = FALSE,
      test_thin = 1
    )
    dir_pdata <- file.path(DIR_PDATA, paste0("model", model_nos_stan[counter]))
    pp_results <- make_post_pred_all(
      draws = run_data$draws,
      day_data = run_data$day_data,
      model = model,
      dir_pdata = dir_pdata,
      day_var = day_vars[counter],
      model_no_stan = model_nos_stan[counter],
      ss = run_data$ss,
      stat = opt$summary_stat,
      testing = FALSE
    )

    tmp_shifted_data <-
      pp_results$shifted_data_by_draw_day_ID %>%
      .[, .(
        day_shifted = median(day_shifted),
        log10_load = median(log10_load)
      ), by = .(ID, day)]

    # Create the ggplot object
    p <- ggplot(tmp_shifted_data, aes(x = day_shifted)) +
      geom_histogram(
        bins = 100, # Equivalent to 'breaks=100' in base R
        fill = "gray",
        color = "gray",
        alpha = 0.7
      ) +
      scale_x_continuous(limits = c(-30, 60)) +
      labs(
        x = "Estimated day of sampling",
        # ggplot2 automatically adds a y-axis label for histograms
        y = "Frequency"
      )

    p <- modified_plot(p, hide_x_label = hide_x_labels[counter])

    plot_list <- c(plot_list, p)
    counter <- counter + 1
    rm(tmp_shifted_data)
    rm(pp_results)
    rm(run_data)
  }

  p_all <- arrangeGrob(plot_annotation(plot_list[[1]], letter = "A"),
    plot_annotation(plot_list[[2]], letter = "B"),
    plot_annotation(plot_list[[3]], letter = "C"),
    plot_annotation(plot_list[[4]], letter = "D"),
    ncol = 2
  )


  ggsave(p_all,
    width = width, height = height * height_factor,
    units = "in",
    filename = file.path(dir_figures, paste0(
      "hist_days_shifted",
      sampling_iteration_suffix,
      ".png"
    )),
    dpi = 300, bg = "white"
  )
  rm(p_all)
}


#' Save summary statistics (ss) from a STAN model to a TSV file.
#'
#' @param model String specifying the model name (includes numeric model 
#' identifier, data filtering criteria, and run settings).
#' @param model_dir Path to the directory containing model fitting results.
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#
#' @return Save a TSV file named "{model}_ss_file.tsv" to the DIR_SS_DATA 
#' directory.
save_ss_file <- function(model,
                         model_dir,
                         model_no_stan) {
  run_data <- load_draws(
    model = model,
    model_dir = model_dir,
    model_no_stan = model_no_stan,
    inc_warmup = FALSE,
    testing = FALSE,
    test_thin = 1
  )

  # Save to TSV
  filepath <- file.path(DIR_SS_DATA, paste0(model, "_ss_file.tsv"))
  write.table(run_data$ss,
    file = filepath, sep = "\t", row.names = FALSE,
    quote = FALSE
  )
}


#' Save summary statistics (ss) from multiple STAN models to a TSV file.
#'
#' @param models A vector of strings specifying the model names (includes 
#' numeric model identifier, data filtering criteria, and run settings).
#' @param model_nos_stan The numbers/identifiers of the stan files used for 
#' fitting.
#' @param model_dirs Paths to the directories containing model fitting results.
#
#' @return Save TSV files named "{model}_ss_file.tsv" to the DIR_SS_DATA 
#' directory.
save_ss_files <- function(models,
                          model_nos_stan,
                          model_dirs) {
  counter <- 1
  for (model in models) {
    save_ss_file(
      model = model,
      model_dir = model_dirs[counter],
      model_no_stan = model_nos_stan[counter]
    )
    counter <- counter + 1
  }
}


main <- function(opt) {
  models <- c(
    paste0("model1_threading_sel3_chains4_n_iter", opt$sampling_iterations, 
           "_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5"),
    paste0("model1_threading_sel6_chains4_n_iter", opt$sampling_iterations, 
           "_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5"),
    paste0("model2_threading_sel3_chains4_n_iter", opt$sampling_iterations, 
           "_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5"),
    paste0("model2_threading_sel6_chains4_n_iter", opt$sampling_iterations, 
           "_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5")
  )
  day_vars <- c("day", "day", "day2", "day2")
  model_nos_stan <- c(1, 1, 2, 2)
  model_dirs <- c("model1", "model1", "model2", "model2")


  if (!dir.exists(DIR_FIGURES_TRAJECTORIES)) {
    dir.create(DIR_FIGURES_TRAJECTORIES)
  }

  save_ss_files(
    models = models,
    model_nos_stan = model_nos_stan,
    model_dirs = model_dirs
  )

  create_multi_model_plots(
    opt = opt,
    models = models,
    day_vars = day_vars,
    model_nos_stan = model_nos_stan,
    model_dirs = model_dirs,
    dir_figures = DIR_FIGURES_TRAJECTORIES
  )
}

main(opt)
