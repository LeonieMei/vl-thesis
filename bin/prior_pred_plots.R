library(ggplot2)
library(here)

TOP <- here()
source(file.path(TOP, "sars2vl", "plot_params.R"))
source(file.path(TOP, "sars2vl", "data_funcs.R"))
source(file.path(TOP, "sars2vl", "post_analysis.R"))


DIR_FIGURES_PRIOR_PRED <- file.path(DIR_FIGURES_VL, "prior_pred")


density_curve_theme <- function() {
  return(
    theme_thesis() +
      theme(
        panel.grid.major.y = element_blank(), # Remove major y-axis gridlines
        panel.grid.minor.y = element_blank(), # Remove minor y-axis gridlines
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  )
}


make_prior_preds <- function(datalist, model_no_stan) {
  model_prefix <- "model"
  model_stan <- paste0(model_prefix, model_no_stan, ".stan")
  sm <- cmdstan_model(file.path(DIR_STAN, model_stan), force_recompile = TRUE)

  datalist$condition_on_data <- 0
  my_inits <- lapply(1:4, function(x) make_TC_inits(datalist))
  priorpred <-
    sm$sample(
      data = datalist,
      iter_warmup = 250,
      iter_sampling = 500,
      chains = 4,
      init = my_inits,
      seed = 123
    )
  draws <- priorpred$draws()
  # save(draws, datalist, day_data, file = pp_file)
  datalist$condition_on_data <- 1
  return(draws)
}


make_prior_pred_plot <- function(draws, tcmp) {
  slope_up_param <- "slope_up_mu"
  slope_down_param <- "slope_down_mu"
  intercept_param <- "intercept_mu"
  alpha_param <- "alpha_mu"

  par(mar = c(2, 1, 0, 1), mgp = c(1, .25, 0), tck = -.01)
  layout(t(matrix(c(1:4, rep(5, 8)), nrow = 2)))

  ## prior distributions for group-level means
  # Example: Convert your density plots to ggplot2
  pd_up <- ggplot(data.frame(x = seq(0, 3, length.out = 100)), aes(x)) +
    stat_function(fun = dlnorm, args = list(meanlog = tcmp$log_slope_up_mu, 
                                            sdlog = tcmp$log_slope_up_sigma)) +
    labs(x = "Proliferation rate", y = "") +
    density_curve_theme()

  pd_int <- ggplot(data.frame(x = seq(4, 12, length.out = 100)), aes(x)) +
    stat_function(fun = dlnorm, args = list(meanlog = tcmp$log_intercept_mu, 
                                            sdlog = tcmp$log_intercept_sigma)) +
    labs(x = "Peak", y = "") +
    density_curve_theme()


  pd_down <- ggplot(data.frame(x = seq(0, 0.9, length.out = 100)), aes(x)) +
    stat_function(fun = dlnorm, args = list(meanlog = tcmp$log_slope_down_mu, 
                                            sdlog = tcmp$log_slope_down_sigma)) +
    labs(x = "-Clearance rate", y = "") +
    density_curve_theme()


  # For pd_intup, generate the density data and plot
  set.seed(123)
  time2peak_data <- abs(rlnorm(1000000, 
                               tcmp$log_intercept_mu, 
                               tcmp$log_intercept_sigma)) /
    rlnorm(1000000, tcmp$log_slope_up_mu, tcmp$log_slope_up_sigma)
  pd_time2peak <- ggplot(data.frame(x = time2peak_data), aes(x)) +
    geom_density() +
    labs(x = "Proliferation time", y = "") +
    xlim(1, 10) +
    density_curve_theme()


  ## prior predictions from group-level means
  plot_viral_trajectories_gg <- function(draws,
                                         intercept_param,
                                         slope_up_param,
                                         slope_down_param,
                                         alpha_param) {
    days <- seq(-10, 40, by = .1)
    model_pars <-
      draws %>%
      subset_draws(slope_up_param) %>%
      as_draws_dt() %>%
      merge(draws %>% subset_draws(slope_down_param) %>% as_draws_dt(), by = ".draw") %>%
      merge(draws %>% subset_draws(intercept_param) %>% as_draws_dt(), by = ".draw") %>%
      merge(draws %>% subset_draws("alpha_mu") %>% as_draws_dt(), by = ".draw")

    # Create a data frame to store trajectories
    trajectories <- tibble(draw = 1:nrow(model_pars)) %>%
      mutate(
        intercept = model_pars[, get(intercept_param)],
        slope_up = model_pars[, get(slope_up_param)],
        slope_down = -model_pars[, get(slope_down_param)],
        alpha_mu = model_pars[, get(alpha_param)]
      ) %>%
      # Generate predictions for each draw and day
      tidyr::expand_grid(days) %>%
      mutate(
        log_viral_load = smooth_function(days, 
                                         intercept, 
                                         slope_up, 
                                         slope_down, 
                                         alpha_mu, 
                                         n_slopes = 2)
      )

    # Create the plot
    p <- ggplot() +
      # Add trajectories
      geom_line(
        data = trajectories, aes(x = days, y = log_viral_load, group = draw),
        color = "black", alpha = 0.05
      ) +
      # Add vertical line at day 0
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
      # Labels and theme
      labs(x = "", y = expression(log[10] ~ viral ~ load)) +
      xlim(-10, 40) +
      ylim(0, 15) +
      ylab(VL_LABEL) +
      xlab(DAYS_LABEL) +
      theme_thesis()


    return(p)
  }

  p <- plot_viral_trajectories_gg(
    draws = draws,
    intercept_param = "intercept_mu",
    slope_up_param = "slope_up_mu",
    slope_down_param = "slope_down_mu",
    alpha_param = "alpha_mu"
  )

  # Arrange the plots in a grid:
  # - Top row: p1, p2, p3, p4 (equally sized)
  # - Bottom row: p5 spanning the width of all 4 top plots
  combined_plot <- (
    (pd_up + pd_int + pd_down + pd_time2peak) /
      p
  ) + plot_layout(heights = c(1.5, 2)) # Adjust heights to make p taller

  return(combined_plot)
}


make_prior_pred_infection_plot <- function(draws, day_data, tcmp) {
  x_limits <- list(
    "slope_up" = c(0, 6), "intercept" = c(3, 15), "slope_down" = c(0, 1),
    "time2peak" = c(0, 15)
  )
  prior_sample_densities <- vector(mode = "list", length = 4)
  names(prior_sample_densities) <-
    c("slope_down", "slope_up", "intercept", "time2peak")
  for (p in names(prior_sample_densities)) {
    prior_sample_densities[[p]] <- draws %>%
      draws_by_id(params = p, day_data = day_data) %>%
      setnames(p, "par") %>%
      .[par > min(0, quantile(par, .005)) & par < max(0, quantile(par, .995))] %>%
      .[, par] %>%
      density()
  }

  # Create ggplot objects for each plot
  x_labels <- c("Proliferation rate", "Peak", "Clearance rate", "Proliferation time")
  p_list <- list()
  counter <- 1
  for (param in c("slope_up", "intercept", "slope_down", "time2peak")) {
    x_factor <- ifelse(param == "slope_down", -1, 1)
    mu_param <- paste0("log_", param, "_mu")
    sigma_param <- paste0("log_", param, "_sigma")
    if (param == "time2peak") {
      time2peak_data <- abs(rlnorm(1000000, 
                                   tcmp$log_intercept_mu, 
                                   tcmp$log_intercept_sigma)) /
        rlnorm(1000000, tcmp$log_slope_up_mu, tcmp$log_slope_up_sigma)
      p <- ggplot(data.frame(x = time2peak_data), aes(x)) +
        geom_density(linetype = 2)
    } else {
      p <- ggplot() +
        stat_function(
          fun = dlnorm, args = list(
            meanlog = tcmp[[mu_param]],
            sdlog = tcmp[[sigma_param]]
          ),
          color = "black", linetype = 2
        )
    }
    p <- p + geom_line(
      data = data.frame(
        x = prior_sample_densities[[param]]$x * x_factor,
        y = prior_sample_densities[[param]]$y
      ),
      aes(x = x, y = y), color = "black"
    ) +
      labs(x = x_labels[counter], y = "") +
      xlim(x_limits[[counter]][1], x_limits[[counter]][2]) +
      density_curve_theme()
    # theme(axis.text.y = element_blank())
    p_list <- c(p_list, p)
    counter <- counter + 1
  }
  # Combine the plots in a 2x2 grid
  p <- (p_list[[1]] + p_list[[2]]) / (p_list[[3]] + p_list[[4]])
  return(p)
}


make_plots <- function(model, model_no_stan, model_dir,
                       tcmp_list) {
  run_data <- load_draws(
    model = model, model_dir = model_dir,
    model_no_stan = model_no_stan,
    inc_warmup = FALSE, testing = FALSE,
    test_thin = 1
  )
  draws <- make_prior_preds(
    datalist = run_data$datalist,
    model_no_stan = model_no_stan
  )
  p <- make_prior_pred_plot(
    draws = draws,
    tcmp = tcmp_list
  )
  p2 <- make_prior_pred_infection_plot(
    draws = draws,
    day_data = run_data$day_data,
    tcmp = tcmp_list
  )
  ggsave(p,
    width = COL_WIDTH * 2, height = HEIGHT_PRIOR_PRED_PLOT * 1.4,
    units = "in",
    filename = file.path(
      DIR_FIGURES_PRIOR_PRED,
      paste0(
        "prior_pred_pop_mean_model_",
        model_no_stan,
        ".png"
      )
    ),
    dpi = 300, bg = "white"
  )

  ggsave(p2,
    width = COL_WIDTH * 2, height = HEIGHT_PRIOR_PRED_PLOT * 0.75,
    units = "in",
    filename = file.path(
      DIR_FIGURES_PRIOR_PRED,
      paste0(
        "prior_pred_infection_model_",
        model_no_stan,
        ".png"
      )
    ),
    dpi = 300, bg = "white"
  )
}


main <- function() {
  tcmp1 <- list(
    # mean, sd
    log_intercept_mu = 2.2, log_intercept_sigma = .15,
    # mean, sd
    log_slope_up_mu = .5, log_slope_up_sigma = .3,
    # mean and sd of the log down-slope
    log_slope_down_mu = -0.8, log_slope_down_sigma = .5
  )

  tcmp_lists <- list(tcmp1)

  # Note: Models 2 and 1 have the same specifications (except for shift 
  # parameter) so we need to only present prior predictions from one of them
  model_stan_nos <- c(1)
  model_dirs <- c("model1")

  models <- c("model1_threading_sel6_chains4_n_iter1000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5")

  counter <- 1
  for (model in models) {
    make_plots(
      model = model,
      model_no_stan = model_stan_nos[counter],
      model_dir = model_dirs[counter],
      tcmp_list = tcmp_lists[[counter]]
    )
    counter <- counter + 1
  }
}

main()
