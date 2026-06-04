library(patchwork)
library(ggplot2)
library(bayesplot)
library(here)
library(grid)
library(gridExtra)


TOP <- here()
source(file.path(TOP, "plot_params.R"))
source(file.path(TOP, "data_params.R"))
source(file.path(TOP, "post_analysis.R"))
source(file.path(TOP, "plot_funcs.R"))


# Functions for making posterior prediction plots

#### Utility functions ####


#' Add progressively shaded credible interval ribbons to ggplot.
#'
#' @param data is a `data.table` or `data.frame` in wide format with mean 
#' and/or median and quantiles lower50, lower55, lower95 ,... upper50, upper55, 
#' upper95 in columns.
#' @param fill is a character indicating the fill colour or a variable in the 
#' data.
#' @param alpha is the transparency level for each ribbon.
#' @param conf_levels The quantiles to add as shaded ribbons.
#'
#' @return A ggplot layer with a progressively shaded confidence ribbon.
conf_ribbon <- function(data,
                        fill = "red",
                        alpha = .05,
                        conf_levels = seq(50, 90, 5)) {
  if (exists(fill, data)) {
    lapply(
      conf_levels,
      function(x) {
        geom_ribbon(
          alpha = alpha, color = NA,
          aes_string(
            ymin = paste0("lower", x),
            ymax = paste0("upper", x), fill = fill
          )
        )
      }
    )
  } else {
    lapply(
      conf_levels,
      function(x) {
        geom_ribbon(
          alpha = alpha, fill = fill, color = NA,
          aes_string(
            ymin = paste0("lower", x),
            ymax = paste0("upper", x)
          )
        )
      }
    )
  }
}


#' Add progressively shaded credible interval ribbons to ggplot.
#'
#' @param data is a `data.table` or `data.frame` in wide format with mean 
#' and/or median and quantiles lower50, lower55, lower95 ,... upper50, upper55, 
#' upper95 in columns.
#' @param fill is a character indicating the fill color or a variable in the 
#' data.
#' @param alpha is the transparency level for each ribbon.
#' @param conf_levels The quantiles to add as shaded ribbons.
#' @return A ggplot layer with a progressively shaded confidence ribbon.
conf_ribbon.hdi <- function(data, fill = "red",
                            alpha = .05,
                            conf_levels = seq(50, 90, 5)) {
  if (exists(fill, data)) {
    lapply(
      conf_levels,
      function(x) {
        geom_ribbon(
          data = data, alpha = alpha, color = NA,
          aes_string(
            ymin = paste0("lower", x[[1]][1]),
            ymax = paste0("upper", x[[1]][1]),
            fill = fill,
            group = ifelse(exists("hdi.group", data), "hdi.group", "NULL")
          )
        )
      }
    )
  } else {
    lapply(
      conf_levels,
      function(x) {
        geom_ribbon(
          data = data, alpha = alpha, fill = fill, color = NA,
          aes_string(
            ymin = paste0("lower", x),
            ymax = paste0("upper", x),
            group = ifelse(exists("hdi.group", data), "hdi.group", "NULL")
          )
        )
      }
    )
  }
}


#' Aggregate draws by day and group.
#'
#' Aggregates daily draws for each category of a given variable, optionally 
#' including CP values.
#'
#' @param by_day_draws A `data.table` with time courses for each posterior draw.
#' @param grp.dt A `data.table` with a column "ID" and a column with a grouping 
#' variable.
#' @param grp.var Character. A grouping variable.
#' @param summary_stat The summary statistic that will be computed for each 
#' category of the grouping variable.
#' @param include_CP A logical value indicating wheter to compute the culture 
#' probability trajectory in addition to the viral load trajectory.
#'
#' @return A data.table with aggregated statistics by day and group.
return_by_day_grp_VL <- function(by_day_draws,
                                 grp.dt,
                                 grp.var,
                                 summary_stat = "median",
                                 include_CP = FALSE) {
  assert_that(summary_stat %in% c("mean", "median"))
  agg_fun <- match.fun(if (summary_stat == "median") collapse::fmedian else collapse::fmean)

  if (include_CP) {
    by_day_grp <-
      do.call(
        rbind,
        lapply(
          unique(grp.dt[[grp.var]]),
          function(x) {
            idx <- grp.dt[get(grp.var) == x, ID]
            by_day_draws[ID %in% idx,
              list(
                log10_load = agg_fun(log10_load),
                CP = agg_fun(CP)
              ),
              by = .(.draw, day_shifted)
            ] %>%
              melt(id.vars = c(".draw", "day_shifted")) %>%
              get_stats(by = c("variable", "day_shifted")) %>%
              .[, (grp.var) := x]
          }
        )
      )
  } else {
    by_day_grp <-
      do.call(
        rbind,
        lapply(
          unique(grp.dt[[grp.var]]),
          function(x) {
            idx <- grp.dt[get(grp.var) == x, ID]
            by_day_draws[ID %in% idx,
              list(log10_load = agg_fun(log10_load)),
              by = .(.draw, day_shifted)
            ] %>%
              melt(id.vars = c(".draw", "day_shifted")) %>%
              get_stats(by = c("variable", "day_shifted")) %>%
              .[, (grp.var) := x]
          }
        )
      )
  }
  return(by_day_grp)
}


#' Plot time course of viral load or culture positivity.
#' Each individual is plotted as a blue line and group level mean and/or median
#' and credible intervals are plotted in red.
#'
#' @param by_day is a data.table with group level statistics.
#' @param by_day_id is a data.table or data.frame with individual level 
#' means/medians.
#' @param y.var is the name of the variable on the y-axis.
#' @param xlim Numeric vector of length 2, specifying the lower and upper 
#' x-axis limits.
#' @param ylim Numeric vector of length 2, specifying the lower and upper 
#' y-axis limits.
#' @param clr is the color for group level results.
#' @param alpha The transparency level for each infection trajectory if 
#' `by_day_id` is not NULL.
#'
#' @return A ggplot layer with a progressively shaded confidence ribbon.
plot_by_day <- function(by_day,
                        by_day_id = NULL,
                        y.var,
                        xlim = c(-11, 37.5),
                        ylim = NULL,
                        clr = "red",
                        alpha = .075) {
  if (exists(clr, by_day)) {
    p <- ggplot(by_day, aes_string(x = "day_shifted", y = y.var, color = clr))
    gl <- geom_line(size = 1.3, alpha = 0.8)
  } else {
    p <- ggplot(by_day, aes_string(x = "day_shifted", y = y.var))
    gl <- geom_line(color = "red", size = 0.75)
  }

  if (!is.null(by_day_id)) {
    p <-
      p + geom_line(data = by_day_id, aes_string(group = "ID"), alpha = alpha, 
                    color = "blue", size = .15)
  }

  p <-
    p +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    gl +
    conf_ribbon(by_day, fill = clr, alpha = 0.08) +
    xlab(DAYS_LABEL) +
    ylab(VL_LABEL) +
    gg_expand()

  return(p)
}


#### Posterior prediction plots ####


#' Create two plots showing simulated viral load trajectories and data points.
#'
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#' @param days Numeric vector. Days relative to the day of peak viral load.
#' @param model_no_stan The number/identifier of the stan file.
#' @param xlim Numeric vector of length 2, specifying the lower and upper 
#' x-axis limits.
#' @param ylim Numeric vector of length 2, specifying the lower and upper 
#' y-axis limits.
#' @param alpha Numeric. The transparency level for lines and data points.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed to plot the population viral load trajectory.

#'
#' @return A named list with two ggplot objects
#'  -`p1`: Simulated viral load trajectories.
#'  -`p2`: Simulated viral load data points.
plot_by_day_sim_data <- function(sim_data,
                                 days,
                                 model_no_stan,
                                 xlim = c(-10, 40),
                                 ylim = c(0, 14),
                                 alpha = .075,
                                 stat = "median") {
  agg_fun <- match.fun(if (stat == "median") median else mean)

  y_list <- list()
  days_list <- list()
  ID_list <- list()

  y_stat <- smooth_function(
    days = days,
    intercept = agg_fun(sim_data$intercept),
    slope_up = agg_fun(sim_data$up_slope),
    slope_down = agg_fun(sim_data$down_slope),
    alpha = 8,
    n_slopes = 2
  )

  y_list <- append(y_list, y_stat)
  days_list <- append(days_list, days)
  ID_list <- append(ID_list, list(rep("avg", length(y_stat))))

  for (d in 1:length(sim_data$intercept)) {
    y <- smooth_function(
      days = days,
      intercept = sim_data$intercept[d],
      slope_up = sim_data$up_slope[d],
      slope_down = sim_data$down_slope[d],
      alpha = 8,
      n_slopes = 2
    )
    y_list <- append(y_list, y)
    days_list <- append(days_list, days)
    ID_list <- append(ID_list, list(rep(d, length(y))))
  }

  # Plot shifted data lines
  VL_by_day_sim_data <- data.table(
    y = unlist(y_list),
    days = unlist(days_list),
    ID = unlist(ID_list)
  )
  VL_by_day_sim_data <- VL_by_day_sim_data %>%
    .[, is_avg := ifelse(ID == "avg", TRUE, FALSE)]

  p <- VL_by_day_sim_data[is_avg == FALSE] %>%
    ggplot(aes(x = days, y = y, group = ID)) +
    geom_line(alpha = alpha, color = "blue", size = .15) +
    labs(x = DAYS_LABEL, y = VL_LABEL) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    gg_expand()

  p <- p +
    geom_line(
      data = VL_by_day_sim_data[is_avg == TRUE],
      color = "red", linewidth = 0.75
    )

  p <- modified_plot(
    p + theme(
      axis.title.y = VL_LABEL,
      axis.title.x = DAYS_LABEL
    ) +
      geom_hline(yintercept = 0, lty = 3),
    hide_x_label = FALSE
  )

  # Plot shifted data points
  vl_var <- "log10_load_true"

  p2 <- modified_plot(
    sim_data$dt %>%
      ggplot(aes(x = day_sampled, y = !!sym(vl_var))) +
      geom_point(alpha = 0.1, shape = 21, stroke = 0, size = 2, fill = "blue") +
      coord_cartesian(xlim = xlim, ylim = ylim) +
      gg_expand() +
      geom_hline(yintercept = 0, lty = 3),
    hide_y_label = TRUE,
    hide_x_label = FALSE
  ) +
    xlab(DAYS_LABEL)
  p2 <- p2 +
    geom_line(
      data = VL_by_day_sim_data[is_avg == TRUE, ],
      aes(x = days, y = y),
      color = "red",
      linewidth = 0.75
    )
  return(list(p1 = p, p2 = p2))
}


#' Plot estimated viral load trajectories and viral load data points.
#'
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param VL_by_day A data.table containing data for the population "average" 
#' viral load trajectory (mean or median).
#' @param shift_draws A data.table with posterior draws of temporal shifts, 
#' grouped by infection.
#' @param shifted_data_by_draw_day_ID A data.table with posterior draws of 
#' shifted days (i.e. estimated days of sampling relative to the day of peak 
#' viral load), and imputed log10 viral loads grouped by infection.
#' @param day_var A string specifying the day variable used in the model (e.g. 
#' "day" with day 0 being the first day with a log10 viral load > 3; 
#' or "day2" with day 0 being the day of the highest viral load measured in the 
#' infection).
#' @param random_IDs A vector of infection identifiers whose viral loads will 
#' be plotted.
#' @param color_by A column name to color data points by (e.g. imputed viral 
#' loads).
#' @param xlim Numeric vector of length 2, specifying the lower and upper 
#' x-axis limits.
#' @param ylim Numeric vector of length 2, specifying the lower and upper 
#' y-axis limits.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed on shift parameters.
#' @param thin Thinning interval for posterior draws (e.g., `thin = 10` selects 
#' every 10th draw).
#'
#' @return A ggplot object.
plot_shifted_data_lines <- function(day_data,
                                    draws,
                                    VL_by_day,
                                    shift_draws,
                                    shifted_data_by_draw_day_ID,
                                    day_var,
                                    random_IDs,
                                    color_by = "imputed",
                                    xlim,
                                    ylim,
                                    stat = stat,
                                    thin = 1) {
  tmp_data <- compute_shifted_data(
    shifted_data_by_draw_day_ID, day_var, stat,
    color_by
  )


  p_shifted_data <-
    tmp_data$tmp_shifted_data %>%
    .[ID %in% random_IDs] %>%
    ggplot(aes(x = day_shifted, y = log10_load, group = ID)) +
    geom_hline(yintercept = 0, lty = 3) +
    geom_line(alpha = 0.1, aes(color = color_column)) +
    geom_point(alpha = tmp_data$alpha, shape = 21, stroke = 0, size = 1, 
               aes(fill = color_column)) +
    red_blue(c(3, 2)) +
    geom_point(data = tmp_data$tmp_shifted_data[ID %in% random_IDs][get(day_var) == 0], 
               shape = 16, size = 0.05) +
    geom_line(data = VL_by_day, size = 0.75, color = "red") +
    theme(legend.position = "none") +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    labs(x = DAYS_LABEL, y = VL_LABEL) +
    gg_expand()

  return(p_shifted_data)
}


#' Plot viral load data points with estimated temporal shifts.
#'
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param VL_by_day A data.table containing data for the population "average" 
#' viral load trajectory (mean or median).
#' @param shift_draws A data.table with posterior draws of temporal shifts, 
#' grouped by infection.
#' @param shifted_data_by_draw_day_ID A data.table with posterior draws of 
#' shifted days (i.e. estimated days of sampling relative to the day of peak 
#' viral load), and imputed log10 viral loads grouped by infection.
#' @param day_var A string specifying the day variable used in the model 
#' (e.g. "day" with day 0 being the first day with a log10 viral load > 3; 
#' or "day2" with day 0 being the day of the highest viral load measured in 
#' the infection).
#' @param random_IDs A vector of infection identifiers whose viral loads will 
#' be plotted.
#' @param color_by A column name to color data points by (e.g. imputed viral 
#' loads).
#' @param xlim Numeric vector of length 2, specifying the lower and upper 
#' x-axis limits.
#' @param ylim Numeric vector of length 2, specifying the lower and upper 
#' y-axis limits.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed on shift parameters.
#' @param thin Thinning interval for posterior draws (e.g., `thin = 10` selects 
#' every 10th draw).
#' @param alpha Transparency of plotted data points.
#'
#' @return A ggplot object.
plot_shifted_data_points <- function(day_data,
                                     draws,
                                     VL_by_day,
                                     shift_draws,
                                     shifted_data_by_draw_day_ID,
                                     day_var,
                                     random_IDs,
                                     color_by = "imputed",
                                     xlim,
                                     ylim,
                                     stat,
                                     thin = 1,
                                     alpha = 0.1) {
  tmp_data <- compute_shifted_data(
    shifted_data_by_draw_day_ID,
    day_var,
    stat,
    color_by
  )

  p_shifted_data <-
    tmp_data$tmp_shifted_data %>%
    .[ID %in% random_IDs] %>%
    ggplot(aes(x = day_shifted, y = log10_load, group = ID)) +
    geom_hline(yintercept = 0, lty = 3) +
    geom_point(
      alpha = alpha, shape = 21, stroke = 0, size = 2,
      aes(fill = color_column)
    ) +
    red_blue(c(3, 2)) +
    geom_line(data = VL_by_day, size = 0.75, color = "red") +
    theme(legend.position = "none") +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    labs(x = DAYS_LABEL, y = VL_LABEL) +
    gg_expand()

  return(p_shifted_data)
}


#' Plot individual viral load trajectories of infections for which shifting 
#' parameter had bad convergence.
#'
#' @param ss A data.table containing summary statistics (including Rhat and 
#' ESS values) from a specific run.
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param shifted_data_by_draw_day_ID A data.table with posterior draws of 
#' shifted days (i.e. estimated days of sampling relative
#' to the day of peak viral load), and imputed log10 viral loads grouped by 
#' infection.
#' @param day_var A string specifying the day variable used in the model 
#' (e.g. "day" with day 0 being the first day with a log10 viral load > 3; 
#' or "day2" with day 0 being the day of the highest viral load measured in the 
#' infection).
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed to plot viral load trajectories.
#'
#' @return A ggplot object.
plot_is_difficult <- function(ss,
                              VLCP_by_day_draws,
                              shifted_data_by_draw_day_ID,
                              day_var,
                              stat = "median") {
  draw_samples <- sample(max(shifted_data_by_draw_day_ID$.draw), 250)

  is.difficult.infections <-
    ss %>%
    .[grepl("^shift\\[", variable) & rhat >= 1.1, ] %>%
    .[, tmp_id := as.numeric(gsub("[^0-9.]", "", variable))] %>%
    .[, ID := day_data[day == 0]$ID[tmp_id]]

  is.difficult.infections <-
    is.difficult.infections %>%
    .[, parameter := tstrsplit(variable, "\\[")[[1]]] %>%
    .[, parameter := paste0(parameter, " (", round(rhat, 1), ")")] %>%
    .[, mean_rhat := mean(rhat), by = .(ID)] %>%
    .[, parameters := paste(parameter, collapse = ", "), by = .(ID)] %>%
    .[, parameter_lb := paste(parameter, collapse = "\n"), by = .(ID)] %>%
    .[, n_pars := .N, by = .(ID)] %>%
    .[, .(ID, tmp_id, parameters, n_pars, parameter_lb, mean_rhat)] %>%
    unique() %>%
    setkeyv("ID") %>%
    .[order(-mean_rhat)]

  difficult.ID <- is.difficult.infections$ID[1:20]

  VL_by_day_ID_difficult <-
    VLCP_by_day_draws[ID %in% difficult.ID] %>%
    summarise_draws_dt_by(
      by = c("day_shifted", "ID"),
      target.var = "log10_load",
      varname = "log10_load",
      stat = stat
    )

  ppc_timecourse <-
    ggplot(VL_by_day_ID_difficult, aes(x = day_shifted, y = log10_load)) +
    geom_ribbon(aes(ymin = q3, ymax = q97), fill = "blue", alpha = .3) +
    geom_line(col = "blue") +
    coord_cartesian(xlim = c(-15, 50), ylim = c(0, 12)) +
    facet_wrap(~ID, ncol = 5) +
    geom_jitter(
      data = shifted_data_by_draw_day_ID[ID %in% difficult.ID & .draw %in% draw_samples],
      height = .05, alpha = .01, color = "red",
      aes(x = day_shifted, y = log10_load)
    ) +
    geom_point(
      data = day_data[ID %in% difficult.ID], aes(x = !!sym(day_var)),
      pch = "x", size = 4
    ) +
    theme_thesis2() +
    xlab(DAYS_LABEL) +
    ylab(VL_LABEL) +
    scale_alpha(range = c(.01, .5))

  return(ppc_timecourse)
}


#' Plot individual viral load trajectories of unusual infections.
#'
#' @param ss A data.table containing summary statistics (including Rhat and ESS 
#' values) from a specific run.
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param shifted_data_by_draw_day_ID A data.table with posterior draws of 
#' shifted days (i.e. estimated days of sampling relative to the day of peak 
#' viral load), and imputed log10 viral loads grouped by infection.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param day_var A string specifying the day variable used in the model 
#' (e.g. "day" with day 0 being the first day with a log10 viral load > 3; 
#' or "day2" with day 0 being the day of the highest viral load measured in the 
#' infection).
#' @param late_positive_infections A logical value indicating whether to only 
#' plot infections with high viral loads late in infection.
#' @param high_positive_infections A logical value indicating whether to only 
#' plot infections with high peak viral load.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed to plot viral load trajectories.
#'
#' @return A ggplot object.
plot_is_weirdly_placed <- function(ss,
                                   VLCP_by_day_draws,
                                   shifted_data_by_draw_day_ID,
                                   day_data,
                                   day_var,
                                   late_positive_infections = FALSE,
                                   high_positive_infections = FALSE,
                                   stat = "median") {
  draw_samples <- sample(max(shifted_data_by_draw_day_ID$.draw), 350)

  weirdly.placed.infections <-
    ss %>%
    .[grepl("^shift\\[", variable)] %>%
    .[, tmp_id := as.numeric(gsub("[^0-9.]", "", variable))] %>%
    .[, ID := day_data[day == 0]$ID[tmp_id]] %>%
    rname(median_shift = median) %>%
    rename(mean_shift = mean) %>%
    merge(day_data[, .(ID, log10_load, value = get(day_var))], by = "ID") %>%
    rename(!!day_var := value) %>%
    .[, mean_day_shifted := get(day_var) + mean_shift] %>%
    .[, median_day_shifted := get(day_var) + median_shift]

  if (late_positive_infections) {
    weirdly_placed_infections_ids <- unique(weirdly.placed.infections[(log10_load > 9) &
      (median_day_shifted > 20)]$ID)
    print(weirdly_placed_infections_ids)
  } else if (high_positive_infections) {
    weirdly_placed_infections_ids <- unique(weirdly.placed.infections[(log10_load > 9.3) &
      (median_day_shifted < 5) &
      (median_day_shifted) > -3]$ID)
  } else {
    weirdly_placed_infections_ids <- unique(weirdly.placed.infections[(log10_load <= 7) &
      (median_day_shifted > 1) &
      (median_day_shifted < 7)]$ID)
  }

  weirdly.placed.infections <-
    weirdly.placed.infections[, .SD[1], by = ID] %>%
    .[ID %in% weirdly_placed_infections_ids] %>%
    .[, parameter := tstrsplit(variable, "\\[")[[1]]] %>%
    .[, parameter := paste0(parameter, " (", round(rhat, 1), ")")] %>%
    .[, mean_rhat := mean(rhat), by = .(ID)] %>%
    .[, parameters := paste(parameter, collapse = ", "), by = .(ID)] %>%
    .[, parameter_lb := paste(parameter, collapse = "\n"), by = .(ID)] %>%
    .[, n_pars := .N, by = .(ID)] %>%
    .[, .(ID, tmp_id, parameters, n_pars, parameter_lb, mean_rhat)] %>%
    unique() %>%
    setkeyv("ID") %>%
    .[order(-mean_rhat)]

  weirdly.placed.infections.ID <- weirdly.placed.infections$ID[1:20]

  VL_by_day_ID_weirdly_placed <-
    VLCP_by_day_draws[ID %in% weirdly.placed.infections.ID] %>%
    summarise_draws_dt_by(
      by = c("day_shifted", "ID"),
      target.var = "log10_load",
      varname = "log10_load",
      stat = stat
    )

  ppc_timecourse <-
    ggplot(VL_by_day_ID_weirdly_placed, aes(x = day_shifted, y = log10_load)) +
    geom_ribbon(aes(ymin = q3, ymax = q97), fill = "blue", alpha = .3) +
    geom_line(col = "blue") +
    coord_cartesian(xlim = c(-15, 50), ylim = c(0, 12)) +
    facet_wrap(~ID, ncol = 5) +
    geom_jitter(
      data = shifted_data_by_draw_day_ID[ID %in% weirdly.placed.infections.ID & .draw %in% draw_samples],
      height = .05, alpha = .01, color = "red",
      aes(x = day_shifted, y = log10_load)
    ) +
    geom_point(
      data = day_data[ID %in% weirdly.placed.infections.ID], 
      aes(x = !!sym(day_var)),
      pch = "x", size = 4
    ) +
    theme_thesis2() +
    xlab(DAYS_LABEL) +
    ylab(VL_LABEL) +
    scale_alpha(range = c(.01, .5))
  return(ppc_timecourse)
}


#### Conditional posterior prediction plots ####


#' Plot time courses of viral load or culture positivity by infection group.
#'
#' @param by_day_grp A data.table with summary statistics (e.g. median) for 
#' posterior predictions of log10 viral load by group and day (post peak load).
#' @param grp.var A string specifying the grouping variable.
#' @param grp.dt A data.table containing the data that were used for fitting.
#' @param y.var A string specifying the variable to plot on the y-axis (will be 
#' "median" if median of predicted log10 viral load was computed).
#' @param show_legend_counts A logical value specifying whether to show 
#' infection counts in the legend.
#'
#' @return A ggplot with time courses and credible intervals.
plot_by_day_grp <- function(by_day_grp,
                            grp.var,
                            grp.dt,
                            y.var = "median",
                            show_legend_counts = TRUE) {
  p <- plot_by_day(by_day_grp,
    clr = grp.var,
    y.var = y.var,
    alpha = 0.075
  )
  if (length(unique(grp.dt[[grp.var]])) == 2) {
    p <- p + red_blue()
  }

  N_table <-
    table(grp.dt[, get(grp.var)]) %>%
    data.table() %>%
    .[N > 0]
  if (show_legend_counts) {
    N_table[, label := paste0(V1, paste0(" (", N, ")"))]
  } else {
    N_table[, label := V1]
  }

  p <-
    p +
    theme(legend.position = c(.8, .8)) +
    gg_legend_size(1)

  legend_title <- paste(gsub("_", " ", grp.var), " (N)")
  if ((class(grp.dt[, get(grp.var)])[1] == "ordered") & (!grepl("age_category", 
                                                                grp.var))) {
    p <- p +
      scale_color_ordinal(name = legend_title, labels = N_table$label) +
      scale_fill_ordinal(name = legend_title, labels = N_table$label)
  } else {
    if (grp.var == "prior_infection") {
      color_values <- c("red", "blue")
    } else if (grp.var == "variant") {
      color_values <- VARIANT_COLORS
      N_table %>%
        .[, V1 := VARIANT_LEGEND_LABELS[V1]]
    } else if (grp.var == "variant2") {
      color_values <- VARIANT2_COLORS
      N_table %>%
        .[, V1 := VARIANT2_LEGEND_LABELS[V1]]
    } else if (grepl("age_category", grp.var)) {
      color_values <- AGE_CODE_COLORS
      N_table %>%
        .[, V1 := AGE_LEGEND_LABELS[V1]]
    } else if (nrow(N_table) == 2) {
      color_values <- c("red", "blue")
    } else if (nrow(N_table) == 3) {
      color_values <- c("gray", "red", "blue") # c("black","red","blue")
    } else if (nrow(N_table) == 4) {
      color_values <- c("#7E6148FF", "#00A087FF", "#FF3CCC", "#0033CC")
    } else {
      color_values <- colorspace::qualitative_hcl(nrow(N_table), 
                                                  palette = "Dark 3")
    }

    if (show_legend_counts) {
      N_table[, label := paste0(V1, paste0(" (", N, ")"))]
    } else {
      N_table[, label := V1]
    }

    p <- p +
      scale_color_manual(name = legend_title, labels = N_table$label, 
                         values = color_values) +
      scale_fill_manual(name = legend_title, labels = N_table$label, 
                        values = color_values)
  }
  return(p)
}


#' Plot time courses of viral load by infection group (predictions are computed 
#' from population-level parameter estimates).
#'
#' @param by_day_draws is a `data.table` with time courses for each infection 
#' group and posterior draw.
#' @param day_data A data.table with PCR test data.
#' @param grp.var The variable specifying the infection subgroups for which 
#' time courses should be plotted.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position of 
#' the legend in the plot (e.g., "top", "right", "bottom", "left", or a numeric 
#' vector like c(0.8, 0.2)).
#' @param xlim Numeric vector of length 2, specifying the lower and upper 
#' x-axis limits.
#' @param ylim Numeric vector of length 2, specifying the lower and upper 
#' y-axis limits.
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#'
#' @return A ggplot with time courses and credible intervals for different 
#' infection subgroups.
plot_by_day_grp_pop <- function(by_day_draws,
                                day_data,
                                grp.var,
                                legend_title = NULL,
                                legend_position = NULL,
                                xlim = c(-10, 40),
                                ylim = c(0, 11),
                                hide_x_label = FALSE,
                                hide_y_label = FALSE,
                                hide_legend = FALSE,
                                add_to_font_size = 0) {
  cats <- unique(by_day_draws[[grp.var]])

  by_day_grp <-
    do.call(
      rbind,
      lapply(
        cats,
        function(x) {
          by_day_draws[get(grp.var) == x, .(log10_load, .draw, day_shifted)] %>%
            melt(id.vars = c(".draw", "day_shifted")) %>%
            get_stats(by = c("variable", "day_shifted")) %>%
            .[, (grp.var) := x]
        }
      )
    )

  VL <- plot_by_day_grp(
    by_day_grp,
    grp.var = grp.var,
    grp.dt = day_data[day == 0],
    show_legend_counts = FALSE
  ) +
    ylab(expression(Log[10] ~ viral ~ load)) +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    labs(x = DAYS_LABEL, y = VL_LABEL) +
    theme(legend.position = c(.75, .825))
  tmp <- gc()

  VL <- modified_plot(VL,
    legend_title = legend_title,
    legend_position = legend_position,
    hide_x_label = hide_x_label,
    hide_y_label = hide_y_label,
    hide_legend = hide_legend,
    add_size = add_to_font_size
  )
  return(VL)
}


#' Combined plot of time courses of viral load and culture positivity by group 
#' and group differences.
#'
#' @param by_day_draws is a `data.table` with time courses for each infection 
#' and posterior draw.
#' @param grp.dt A `data.table` with a column "ID" and a column with a grouping 
#' variable.
#' @param stat is "RD" and "RR" for risk difference and risk ratio, 
#' respectively.
#' @param comparisons An optional matrix with pairwise comparisons.
#' @return A ggplot with three panels: viral load time course, culture 
#' positivity time course and group differences.
plot_by_day_grp_delta <- function(by_day_draws,
                                  grp.dt,
                                  stat = "RD",
                                  summary_stat = "median",
                                  comparisons = NULL,
                                  legend_title = NULL,
                                  include_CP = FALSE,
                                  legend_position = "right",
                                  xlim = c(-10, 40),
                                  ylim = c(0, 11),
                                  hide_x_label = FALSE,
                                  hide_y_label = FALSE,
                                  hide_legend = FALSE,
                                  add_to_font_size = 0) {
  assert_that(summary_stat %in% c("mean", "median"))
  agg_fun <- match.fun(if (summary_stat == "median") collapse::fmedian else collapse::fmean)

  grp.var <- grep("ID", names(grp.dt), value = TRUE, invert = TRUE)
  if (include_CP) {
    by_day_grp <-
      do.call(
        rbind,
        lapply(
          unique(grp.dt[[grp.var]]),
          function(x) {
            idx <- grp.dt[get(grp.var) == x, ID]
            by_day_draws[ID %in% idx,
              list(
                log10_load = agg_fun(log10_load),
                CP = agg_fun(CP)
              ),
              by = .(.draw, day_shifted)
            ] %>%
              melt(id.vars = c(".draw", "day_shifted")) %>%
              get_stats(by = c("variable", "day_shifted")) %>%
              .[, (grp.var) := x]
          }
        )
      )
  } else {
    by_day_grp <-
      do.call(
        rbind,
        lapply(unique(grp.dt[[grp.var]]), function(x) {
          idx <- grp.dt[get(grp.var) == x, ID]
          by_day_draws[ID %in% idx,
            list(log10_load = agg_fun(log10_load)),
            by = .(.draw, day_shifted)
          ] %>%
            get_stats(by = "day_shifted", var = "log10_load") %>%
            .[, (grp.var) := x]
        })
      )
    VL <- plot_by_day_grp(
      by_day_grp,
      grp.var = grp.var, grp.dt = grp.dt
    ) +
      ylab(expression(Log[10] ~ viral ~ load)) +
      coord_cartesian(ylim = ylim, xlim = xlim) +
      labs(x = DAYS_LABEL, y = VL_LABEL) +
      theme(legend.position = c(.75, .825))
    tmp <- gc()
  }

  if (include_CP) {
    VL <- plot_by_day_grp(
      by_day_grp[variable == "log10_load"],
      grp.var = grp.var,
      grp.dt = grp.dt
    ) +
      ylab(expression(Log[10] ~ viral ~ load)) +
      coord_cartesian(ylim = ylim, xlim = xlim) +
      theme(legend.position = c(.75, .825))
    tmp <- gc()
    CP <- plot_by_day_grp(
      by_day_grp[variable == "CP"],
      grp.var = grp.var,
      grp.dt = grp.dt
    ) +
      ylab("Probability of positive culture") +
      coord_cartesian(xlim = xlim, ylim = c(-.005, 1)) +
      theme(legend.position = c(.75, .825))
    tmp <- gc()
    Peak_delta <-
      plot_delta_grps(by_day_draws[day_shifted == 0],
        grp.dt = grp.dt,
        y.var = "CP",
        stat = stat,
        comparisons = comparisons
      ) +
      ylab("Difference of peak culture probability")
    tmp <- gc()
    CP <- modified_plot(CP,
      legend_title = legend_title,
      legend_position = legend_position,
      hide_x_label = hide_x_label,
      hide_y_label = hide_y_label,
      hide_legend = hide_legend,
      add_size = add_to_font_size
    )
    Peak_delta <- modified_plot(Peak_delta,
      legend_title = legend_title,
      legend_position = legend_position,
      hide_x_label = hide_x_label,
      hide_y_label = hide_y_label,
      hide_legend = hide_legend,
      add_size = add_to_font_size
    )
  }

  VL <- modified_plot(VL,
    legend_title = legend_title,
    legend_position = legend_position,
    hide_x_label = hide_x_label,
    hide_y_label = hide_y_label,
    hide_legend = hide_legend,
    add_size = add_to_font_size
  )

  if (include_CP) {
    final_plot <- plot_grid(VL, CP, Peak_delta, hjust = 1, vjust = 1, ncol = 1)
    # create common x labels
    x.grob <- textGrob("Days from peak viral load",
      gp = gpar(fontsize = TEXT_SIZE)
    )

    # add to plot
    final_plot <- grid.arrange(arrangeGrob(final_plot, bottom = x.grob))
    return(final_plot)
  } else {
    return(VL)
  }
}


#' Plot of time courses of viral load by group and group differences.
#' @param by_day_draws is a `data.table` with time courses for each infection 
#' and posterior draw.
#' @param grp.dt A `data.table` with a column "ID" and a column with a grouping 
#' variable.
#' @param stat is "RD" and "RR" for risk difference and risk ration, 
#' respectively.
#' @param comparisons An optional matrix with pairwise comparisons.
#' @return A ggplot with one panel: viral load time course.
plot_by_day_grp_VL <- function(by_day_draws,
                               grp.dt,
                               stat = "RD",
                               comparisons = NULL,
                               xlim = c(-11, 37.5),
                               ylim = c(0, 11)) {
  grp.var <- grep("ID", names(grp.dt), value = TRUE, invert = TRUE)
  by_day_grp <- return_by_day_grp_VL(by_day_draws, grp.dt, grp.var)

  VL <- plot_by_day_grp(
    by_day_grp[variable == "log10_load"],
    grp.var = grp.var,
    grp.dt = grp.dt
  ) +
    ylab(expression(Log[10] ~ viral ~ load)) +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    theme(legend.position = c(.75, .825))
  tmp <- gc()
  return(VL)
}


#' Plot group differences in viral load or culture positivity of time course 
#' analysis.
#'
#' @param by_day_draws is a `data.table` with time courses for each infection 
#' and posterior draw.
#' @param y.var The variable for which the difference is calculated.
#' @param grp.dt A `data.table` with a column "ID" and a column with a grouping 
#' variable.
#' @param stat is "RD" and "RR" for risk difference and risk ratio, 
#' respectively.
#' @return A ggplot with group differences and credible intervals.
plot_delta_grps <- function(by_day_draws,
                            y.var = "log10_load",
                            grp.dt,
                            stat = "RD",
                            comparisons = NULL,
                            plot = TRUE,
                            summary_stat = "median") {
  agg_fun <- match.fun(if (summary_stat == "median") median else mean)

  grp.var <- setdiff(names(grp.dt), "ID")
  cast_formula <- as.formula(paste0(".draw ~", grp.var))

  setkeyv(by_day_draws, "ID")
  setkeyv(grp.dt, "ID")
  tmp <-
    by_day_draws[day_shifted == 0] %>%
    .[grp.dt, c(grp.var) := get(grp.var)] %>%
    .[, .(value = agg_fun(get(y.var))), by = c(".draw", grp.var)] %>%
    .[!is.na(get(grp.var))] %>%
    dcast(cast_formula)
  if (is.null(comparisons)) {
    comparisons <- combn(names(tmp)[-1], 2)
  }

  comp_stats <- c()
  for (k in 1:ncol(comparisons)) {
    if (stat == "RR") { ## risk ratio
      delta <- tmp[, get(comparisons[1, k])] / tmp[, get(comparisons[2, k])]
    } else if (stat == "RD") { ## risk difference
      delta <- tmp[, get(comparisons[1, k])] - tmp[, get(comparisons[2, k])]
    }
    comp_stats <- rbind(
      comp_stats,
      get_stats(data.table(delta = delta), var = "delta")
    )
  }
  rm(tmp)
  tmp <- gc()
  comp_stats <-
    data.table(comp_stats) %>%
    .[, comparison := paste0(
      comparisons[1, ],
      ifelse(stat == "RD", "–", "/"),
      ifelse(ncol(comparisons) > 22, "\n", ""), comparisons[2, ]
    )] %>%
    .[, x := ncol(comparisons):1]

  if (stat == "RR") {
    ylim <- c(0.5, 1.1)
  } else if (stat == "RD") {
    offset <- (max(comp_stats[, upper90]) - min(comp_stats[, lower90])) * .5
    ylim <- c(min(comp_stats[, lower90]) - offset, max(comp_stats[, upper90]))
  }

  if (plot == TRUE) {
    return(
      ggplot(comp_stats, aes(x = x, y = .data[[summary_stat]])) +
        geom_hline(
          yintercept = ifelse(stat == "RR", 1, 0),
          col = "red", lty = 3, size = .5
        ) +
        geom_point(color = "black") +
        conf_linerange() +
        theme(
          axis.line.y = element_blank(),
          axis.text.y = element_blank()
        ) +
        geom_text(
          aes(
            y = min(comp_stats[, lower90]) - offset * .75,
            label = comparison
          ),
          size = 6,
          hjust = 0
        ) +
        xlab("Comparison") +
        ylab(stat) +
        coord_cartesian(ylim = ylim) +
        coord_flip()
    )
  } else {
    return(
      comp_stats
    )
  }
}


#' Create a plot showing conditional posterior predictions of viral load 
#' trajectories for infections with different viral variants.
#' Predictions are computed from infection-level parameter estimates.
#'
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position 
#' of the legend in the plot (e.g., "top", "right", "bottom", "left", or a 
#' numeric vector like c(0.8, 0.2)).
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#'
#' @return A ggplot object.
make_variant_plot <- function(VLCP_by_day_draws,
                              day_data,
                              legend_title = NULL,
                              legend_position = NULL,
                              hide_x_label = FALSE,
                              hide_y_label = FALSE,
                              hide_legend = FALSE,
                              add_to_font_size = 0) {
  # Group by variant.
  TC_by_variant <-
    plot_by_day_grp_delta(VLCP_by_day_draws,
      grp.dt = unique(day_data[, .(ID, variant)]),
      legend_title = legend_title,
      legend_position = legend_position,
      hide_x_label = hide_x_label,
      hide_y_label = hide_y_label,
      hide_legend = hide_legend,
      add_to_font_size = add_to_font_size
    )

  return(TC_by_variant)
}


#' Create a plot showing conditional posterior predictions computed from 
#' population-level parameter estimates of viral load trajectories for 
#' infections with different viral variants.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting, 
#' with results stored in draws.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position 
#' of the legend in the plot (e.g., "top", "right",
#' "bottom", "left", or a numeric vector like c(0.8, 0.2)).
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#' @param reference_variant A string indicating the variant used for prediction.
#' @param reference_hosp_status A string indicating the disease severity 
#' status (PAMS or hospitalized) used for prediction.
#'
#' @return A ggplot object.
make_variant_plot_post_pred_pop <- function(draws,
                                            day_data,
                                            legend_title = NULL,
                                            legend_position = NULL,
                                            hide_x_label = FALSE,
                                            hide_y_label = FALSE,
                                            hide_legend = FALSE,
                                            add_to_font_size = 0,
                                            reference_variant = "omicron",
                                            reference_hosp_status = "non-PAMS, non-hosp.") {
  # Population-level posterior prediction
  VLCP_by_day_draws_variant <- make_VLCP_by_grp("variant",
    draws = draws,
    days = seq(-10, 40, by = 1),
    thin = 1,
    reference_variant = reference_variant,
    reference_hosp_status = reference_hosp_status
  )
  TC_by_variant <- plot_by_day_grp_pop(VLCP_by_day_draws_variant,
    day_data = day_data,
    grp.var = "variant",
    legend_title = legend_title,
    legend_position = legend_position,
    hide_x_label = hide_x_label,
    hide_y_label = hide_y_label,
    hide_legend = hide_legend,
    add_to_font_size = add_to_font_size
  )
  return(TC_by_variant)
}


#' Create a plot showing posterior predictions of viral load trajectories for 
#' infections with different viral variants (including Omicron sublineages). 
#' Posterior predictions are computed from infection-level parameter estimates.
#'
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position 
#' of the legend in the plot (e.g., "top", "right", "bottom", "left", or a 
#' numeric vector like c(0.8, 0.2)).
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#'
#' @return A ggplot object.
make_variant2_plot <- function(VLCP_by_day_draws,
                               day_data,
                               legend_title = NULL,
                               legend_position = NULL,
                               hide_x_label = FALSE,
                               hide_y_label = FALSE,
                               hide_legend = FALSE,
                               add_to_font_size = 0) {
  TC_by_variant2 <-
    plot_by_day_grp_delta(VLCP_by_day_draws,
      grp.dt = unique(day_data[, .(ID, variant2)]),
      legend_title = legend_title,
      legend_position = legend_position,
      hide_x_label = hide_x_label,
      hide_y_label = hide_y_label,
      hide_legend = hide_legend,
      add_to_font_size = add_to_font_size
    )

  return(TC_by_variant2)
}


#' Create a plot showing posterior predictions of viral load trajectories for 
#' infections of different disease severity.
#' Predictions are computed from infection-level parameter estimates.
#'
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position of 
#' the legend in the plot (e.g., "top", "right", "bottom", "left", or a numeric 
#' vector like c(0.8, 0.2)).
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#'
#' @return A ggplot object.
make_pams_plot <- function(VLCP_by_day_draws,
                           day_data,
                           hospitalized_var = "hospitalized",
                           legend_title = NULL,
                           legend_position = NULL,
                           hide_x_label = FALSE,
                           hide_y_label = FALSE,
                           hide_legend = FALSE,
                           add_to_font_size = 0) {
  # Group by PAMS status.
  grp.dt.Group <-
    unique(day_data[, .(ID, PAMS1, get(hospitalized_var))])
  grp.dt.Group <- set_names(grp.dt.Group, c("ID", "PAMS1", "hospitalized"))
  grp.dt.Group %>%
    .[, Group := ifelse(PAMS1 == TRUE,
      "PAMS",
      ifelse(hospitalized == TRUE,
        "Hospitalized", "Non-PAMS, non-hosp."
      )
    )] %>%
    .[, Group := factor(Group, levels = c("Non-PAMS, non-hosp.", 
                                          "Hospitalized", 
                                          "PAMS"))] %>%
    .[, PAMS1 := NULL] %>%
    .[, hospitalized := NULL]

  TC_by_PAMSHosp <-
    plot_by_day_grp_delta(VLCP_by_day_draws,
      grp.dt = grp.dt.Group,
      legend_title = legend_title,
      legend_position = legend_position,
      hide_x_label = hide_x_label,
      hide_y_label = hide_y_label,
      hide_legend = hide_legend,
      add_to_font_size = add_to_font_size
    )

  return(TC_by_PAMSHosp)
}


#' Create a plot showing conditional posterior predictions computed from 
#' population-level parameter estimates of viral load trajectories for 
#' infections of different disease severity.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting, 
#' with results stored in draws.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position of 
#' the legend in the plot (e.g., "top", "right", "bottom", "left", or a numeric 
#' vector like c(0.8, 0.2)).
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#' @param reference_variant A string indicating the variant used for prediction.
#' @param reference_hosp_status A string indicating the disease severity status 
#' (PAMS or hospitalized) used for prediction.
#'
#' @return A ggplot object.
make_pams_plot_post_pred_pop <- function(draws,
                                         day_data,
                                         pams_var = "PAMS3",
                                         legend_title = NULL,
                                         legend_position = NULL,
                                         hide_x_label = FALSE,
                                         hide_y_label = FALSE,
                                         hide_legend = FALSE,
                                         add_to_font_size = 0,
                                         reference_variant = "omicron",
                                         reference_hosp_status = "non-PAMS, non-hosp.") {
  assert_that(pams_var %in% c("PAMS1", "PAMS3", "hospitalized"))

  # Make population-level posterior predictions (computed from population-level 
  # parameter estimates, not on infection-level parameter estimates).
  VLCP_by_day_draws_PAMS <- make_VLCP_by_grp(pams_var,
    draws = draws,
    days = seq(-10, 40, by = 1),
    thin = 1,
    reference_variant = reference_variant,
    reference_hosp_status = reference_hosp_status
  )
  if (pams_var == "PAMS3") {
    day_data <- copy(day_data)
    day_data[, PAMS3 := factor(PAMS3,
      levels = c(3, 2, 1),
      labels = c("Non-PAMS, non-hosp.", "Hospitalized", "PAMS")
    )]
  }
  TC_by_PAMS <- plot_by_day_grp_pop(
    by_day_draws = VLCP_by_day_draws_PAMS,
    day_data = day_data,
    grp.var = pams_var,
    legend_title = legend_title,
    legend_position = legend_position,
    hide_x_label = hide_x_label,
    hide_y_label = hide_y_label,
    hide_legend = hide_legend,
    add_to_font_size = add_to_font_size
  )

  return(TC_by_PAMS)
}


#' Create a plot showing conditional posterior predictions of viral load 
#' trajectories for infections of individuals in different age categories.
#' Predictions are computed from infection-level parameter estimates.
#'
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position of 
#' the legend in the plot (e.g., "top", "right", "bottom", "left", or a numeric 
#' vector like c(0.8, 0.2)).
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#'
#' @return A ggplot object.
make_age_plot <- function(VLCP_by_day_draws,
                          day_data,
                          legend_title = NULL,
                          legend_position = NULL,
                          hide_x_label = FALSE,
                          hide_y_label = FALSE,
                          hide_legend = FALSE,
                          add_to_font_size = 0) {
  # Group by age group.
  day_data %>%
    .[, age_category3 := factor(
      cut(age,
        breaks = AGE_BREAKS3,
        ordered_result = TRUE, right = FALSE
      ),
      levels = AGE_LEVELS3,
      labels = AGE_LABELS3,
      ordered = TRUE
    )]

  TC_by_age <-
    plot_by_day_grp_delta(VLCP_by_day_draws,
      grp.dt = unique(day_data[, .(ID, age_category3)]),
      legend_title = legend_title,
      legend_position = legend_position,
      hide_x_label = hide_x_label,
      hide_y_label = hide_y_label,
      hide_legend = hide_legend,
      add_to_font_size = add_to_font_size
    )

  return(TC_by_age)
}


#' Create a plot showing conditional posterior predictions computed from 
#' population-level parameter estimates of viral load trajectories for 
#' infections of individuals in different age categories.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting, 
#' with results stored in draws.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position 
#' of the legend in the plot (e.g., "top", "right",
#' "bottom", "left", or a numeric vector like c(0.8, 0.2)).
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#' @param reference_variant A string indicating the variant used for prediction.
#' @param reference_hosp_status A string indicating the disease severity status 
#' (PAMS or hospitalized) used for prediction.
#'
#' @return A ggplot object.
make_age_plot_post_pred_pop <- function(draws, day_data,
                                        legend_title = NULL,
                                        legend_position = "right",
                                        hide_x_label = FALSE,
                                        hide_y_label = FALSE,
                                        hide_legend = FALSE,
                                        add_to_font_size = 0,
                                        reference_variant = "omicron",
                                        reference_hosp_status = "non-PAMS, non-hosp.") {
  # Population-level posterior prediction
  VLCP_by_day_draws_age <- make_VLCP_by_grp("age_category3",
    draws = draws,
    days = seq(-10, 40, by = 1),
    thin = 1,
    reference_variant = reference_variant,
    reference_hosp_status = reference_hosp_status
  )
  TC_by_age <- plot_by_day_grp_pop(VLCP_by_day_draws_age,
    day_data = day_data,
    grp.var = "age_category3",
    legend_title = legend_title,
    legend_position = legend_position,
    hide_x_label = hide_x_label,
    hide_y_label = hide_y_label,
    hide_legend = hide_legend,
    add_to_font_size = add_to_font_size
  )
  return(TC_by_age)
}


#' Create a plot showing posterior predictions of viral load trajectories for 
#' infections of individuals that were immune naive (to our knowledge) or had a 
#' prior infection. Posterior predictions are computed from infection-level 
#' parameter estimates.
#'
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position of 
#' the legend in the plot (e.g., "top", "right", "bottom", "left", or a numeric 
#' vector like c(0.8, 0.2)).
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#'
#' @return A ggplot object.
make_prior_infection_plot <- function(VLCP_by_day_draws,
                                      day_data = day_data,
                                      legend_title = NULL,
                                      legend_position = "right",
                                      hide_x_label = FALSE,
                                      hide_y_label = FALSE,
                                      hide_legend = FALSE,
                                      add_to_font_size = 0) {
  TC_by_prior_infection <-
    plot_by_day_grp_delta(VLCP_by_day_draws,
      grp.dt = unique(day_data[, .(ID, prior_infection)]),
      legend_title = legend_title,
      legend_position = legend_position,
      hide_x_label = hide_x_label,
      hide_y_label = hide_y_label,
      hide_legend = hide_legend,
      add_to_font_size = add_to_font_size
    )
  return(TC_by_prior_infection)
}


#' Create a plot showing conditional posterior predictions computed from 
#' population-level parameter estimates of viral load trajectories for 
#' infections of individuals that were immune naive (to our knowledge) or had a 
#' prior infection.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting, 
#' with results stored in draws.
#' @param legend_title The title of the figure legend. Specify NULL for a 
#' legend with no title.
#' @param legend_position Character or numeric vector. Controls the position of 
#' the legend in the plot (e.g., "top", "right", "bottom", "left", or a numeric 
#' vector like c(0.8, 0.2)).
#' @param hide_x_label A logical value indicating whether to hide the x-axis 
#' label.
#' @param hide_y_label A logical value indicating whether to hide the y-axis 
#' label.
#' @param hide_legend A logical value indicating whether to hide the figure 
#' legend.
#' @param add_to_font_size Numeric. Adjusts the base font size for all text 
#' elements in the plot.
#' @param reference_variant A string indicating the variant used for prediction.
#' @param reference_hosp_status A string indicating the disease severity status 
#' (PAMS or hospitalized) used for prediction.
#'
#' @return A ggplot object.
make_prior_infection_plot_post_pred_pop <- function(draws,
                                                    day_data,
                                                    legend_title = NULL,
                                                    legend_position = "right",
                                                    hide_x_label = FALSE,
                                                    hide_y_label = FALSE,
                                                    hide_legend = FALSE,
                                                    add_to_font_size = 0,
                                                    reference_variant = "omicron",
                                                    reference_hosp_status = "non-PAMS, non-hosp.") {
  # Population-level posterior prediction
  VLCP_by_day_draws_prior_infection <- make_VLCP_by_grp("prior_infection",
    draws = draws,
    days = seq(-10, 40, by = 1),
    thin = 1,
    reference_variant = reference_variant,
    reference_hosp_status = reference_hosp_status
  )
  TC_by_prior_infection <- plot_by_day_grp_pop(VLCP_by_day_draws_prior_infection,
    day_data = day_data,
    grp.var = "prior_infection",
    legend_title = legend_title,
    legend_position = legend_position,
    hide_x_label = hide_x_label,
    hide_y_label = hide_y_label,
    hide_legend = hide_legend,
    add_to_font_size = add_to_font_size
  )
  return(TC_by_prior_infection)
}


#' Create and save plots showing predicted viral load trajectories for 
#' different infection subgroups. Predictions computed from infection-level 
#' parameter estimates, and those computed from population-level estimates are 
#' plotted.
#'
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting, 
#' with results stored in draws.
#' @param dir_figures A string specifying the directory path where figures will 
#' be saved.
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param model A string specifying the name of the model (includes numeric 
#' model identifiers, data filtering criteria, and run settings).
#' @param output_file A string specifying the path to the file where posterior 
#' predictions for different infection subgroups will be stored.
#' @param thesis A logical value specifying whether figures are intended for 
#' use in PhD thesis.
make_post_pred_by_group <- function(VLCP_by_day_draws,
                                    draws,
                                    day_data,
                                    dir_figures,
                                    model_no_stan,
                                    model,
                                    output_file,
                                    thesis = TRUE) {
  # Plotting parameters
  include_CP <- FALSE
  width <- ifelse(thesis, VL_PLOT_WIDTH_THESIS, VL_PLOT_WIDTH) # 30, 10
  height <- ifelse(thesis, VL_PLOT_HEIGHT_THESIS, VL_PLOT_HEIGHT)
  height3 <- height * 3

  dir_pdata <- file.path(DIR_PDATA, paste0("model", model_no_stan))

  # The following is to compute posterior predictions of main model parameters 
  # for different subgroups, computed from population-level parameter estimates.
  write_post_pred_groups(
    draws = draws,
    output_file = output_file,
    model = model,
    reference_variant = "wildtype",
    reference_hosp_status = "hospitalized",
    model_no_stan = model_no_stan
  )
  write_post_pred_groups(
    draws = draws,
    output_file = output_file,
    model = model,
    reference_variant = "omicron",
    reference_hosp_status = "hospitalized",
    model_no_stan = model_no_stan
  )
  write_post_pred_groups(
    draws = draws,
    output_file = output_file,
    model = model,
    reference_variant = "wildtype",
    reference_hosp_status = "PAMS",
    model_no_stan = model_no_stan
  )
  write_post_pred_groups(
    draws = draws,
    output_file = output_file,
    model = model,
    reference_variant = "omicron",
    reference_hosp_status = "PAMS",
    model_no_stan = model_no_stan
  )
}


#### Main functions ####


#' Create and save plots showing posterior predictions for viral load 
#' trajectories, computed from infection-level parameter estimates.
#' Save posterior estimates (mean, median, or HPDIs) or alternatively 
#' conditional predictions; save posterior predictions computed from 
#' infection-level parameter estimates.
#'
#' @param VLCP_by_day_draws A data.table with posterior predictions of viral 
#' load trajectories, computed from infection-level parameter estimates.
#' @param VLCP_by_day_ID A data.table with infection-wise posterior predictions 
#' of median viral load trajectories.
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param ss A data.table containing summary statistics (including Rhat and
#'  ESS values) from a specific run.
#' @param shift_draws A data.table with posterior draws of temporal shifts, 
#' grouped by infection.
#' @param shifted_data_by_draw_day_ID A data.table with posterior draws of 
#' shifted days (i.e. estimated days of sampling relative to the day of peak 
#' viral load), and imputed log10 viral loads grouped by infection.
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param dir_pdata A string specifying the directory path where posterior 
#' predictions are stored.
#' @param dir_figures A string specifying the directory path where figures 
#' will be saved.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param n_min_tests The minimum number of PCR days per infection that were 
#' required for inclusion in the run.
#' @param n_samples A numeric (or NULL) specifying the number of infections 
#' that were simulated (if sim_data not NULL, else NULL).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param day_var A string specifying the day variable used in the model 
#' (e.g. "day" with day 0 being the first day with a log10 viral load > 3; 
#' or "day2" with day 0 being the day of the highest viral load measured in the 
#' infection).
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#' @param testing A logical value. If `TRUE`, loads only a subset of posterior 
#' draws (for testing/debugging).
#' @param annotate_plot A logical value indicating whether to annotate the plot 
#' (with letters A, B, C, D).
#' @param output_file A path to the file where posterior estimates (medians and 
#' HPDIs), and posterior predictions will be stored. The file name will be 
#' extended in the case of posterior predictions, so that these are stored in a 
#' separate file.
#' @param output_file2 A path to the file where differences between estimated 
#' and true parameter values (medians and HPDIs), and differences between 
#' posterior predictions and true observations will be stored (if simulated 
#' data were used for model fitting). The file name will be extended in the 
#' case of posterior predictions, so that these are stored in a separate file. 
#' Can be NULL if real (not simulated) data were used.
#' @param color_by A column name to color data points by (e.g. imputed viral 
#' loads).
#' @param plot_lines A logical value indicating whether to plot shifted data 
#' points or estimated infection-wise viral load trajectories.
#' @param plot_n_infections Randomly select n infections to make posterior 
#' prediction plots for.
#' @param stat The summary statistic (mean or median) that will be computed 
#' for posterior estimates and predictions.
#' @param thesis A logical value specifying whether figures are intended for 
#' use in PhD thesis.
make_post_pred_summary_stat <- function(VLCP_by_day_draws,
                                        VLCP_by_day_ID,
                                        draws,
                                        ss,
                                        shift_draws,
                                        shifted_data_by_draw_day_ID,
                                        day_data,
                                        dir_pdata,
                                        dir_figures,
                                        model,
                                        n_min_tests,
                                        n_samples,
                                        model_no_stan,
                                        day_var,
                                        sim_data = NULL,
                                        testing = FALSE,
                                        annotate_plot = FALSE,
                                        output_file = NULL,
                                        output_file2 = NULL,
                                        color_by = "imputed",
                                        plot_lines = FALSE,
                                        plot_n_infections = NULL,
                                        stat = "median",
                                        thesis = TRUE) {
  width <- ifelse(thesis, VL_PLOT_WIDTH_THESIS, VL_PLOT_WIDTH) # 30, 10
  height <- ifelse(thesis, VL_PLOT_HEIGHT_THESIS, VL_PLOT_HEIGHT)
  height3 <- height * 3

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
      ylim = ylim, xlim = xlim
    ) +
      geom_hline(yintercept = 0, lty = 3),
    hide_x_label = TRUE
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
      color_by = color_by, xlim = xlim, ylim = ylim,
      stat = stat,
      thin = 1
    )
  } else {
    p2 <- plot_shifted_data_points(
      day_data = day_data, draws = draws,
      VL_by_day = VL_by_day,
      shift_draws = shift_draws,
      shifted_data_by_draw_day_ID = shifted_data_by_draw_day_ID,
      day_var = day_var,
      random_IDs = random_IDs,
      color_by = color_by, xlim = xlim, ylim = ylim,
      stat = stat, thin = 1
    )
  }
  # Plot underlying simulated viral load courses.
  if (!is.null(sim_data)) {
    if (testing == TRUE) {
      days <- c(-7:4, seq(5, 25, 2.5), 40)
    } else {
      days <- c(
        seq(-11, -7, by = 1), seq(-6.5, -3, by = .5),
        seq(-2.75, 0, by = .25), seq(0.5, 5, by = .5),
        seq(6, 25, 1), 27, 30, 35, 40
      )
    }
    p2 <- modified_plot(p2,
      hide_x_label = TRUE,
      hide_y_label = TRUE,
      hide_legend = TRUE
    )

    p_sim <- plot_by_day_sim_data(sim_data,
      days = tc_days,
      xlim = xlim,
      ylim = ylim,
      model_no_stan = model_no_stan
    )

    p3 <- p_sim$p1
    p4 <- p_sim$p2
  } else {
    p2 <- modified_plot(p2, hide_legend = TRUE)
  }
  if (annotate_plot == TRUE) {
    stat_params_estimated <- calc_main_params_summary_stat(draws,
      day_data = day_data,
      stat = stat,
      stat2 = stat
    )
    xpos_label <- ifelse(is.null(sim_data), 23, 17.5)
    ypos_label <- ifelse(is.null(sim_data), 11, 11.75)
    p1 <- annotate_summary_stat_plot(stat_params_estimated, p1,
      xpos_label = xpos_label,
      ypos_label = ypos_label,
      sim_data = !is.null(sim_data)
    )
    if (!is.null(sim_data)) {
      stat_params_sim <- calc_main_params_summary_stat_sim_data(sim_data,
        day_data = day_data,
        model_no_stan = model_no_stan,
        stat = stat, stat2 = stat
      )
      p3 <- annotate_summary_stat_plot(stat_params_sim, p3,
        xpos_label = xpos_label,
        ypos_label = ypos_label,
        sim_data = !is.null(sim_data)
      )
    }
  }
  save_param_stats(
    draws = draws, day_data = day_data, ss = ss, stat = stat,
    stat2 = stat, sim_data = sim_data, output_file = output_file,
    output_file2 = output_file2, model = model,
    model_no_stan = model_no_stan, n_min_tests = n_min_tests,
    n_samples = n_samples
  )
  # The following is to compute posterior predictions of subgroups, 
  # computed from infection-level parameter estimates.
  save_param_stats_post_pred_inf(
    draws = draws, day_data = day_data,
    stat = stat, stat2 = stat,
    output_file = output_file,
    model = model,
    model_no_stan = model_no_stan,
    n_min_tests = n_min_tests
  )

  if (!is.null(sim_data)) {
    ### combine all plots
    final_plot <- arrangeGrob(plot_annotation(p1, letter = "A"),
      plot_annotation(p2, letter = "B"),
      plot_annotation(p3, letter = "C"),
      plot_annotation(p4, letter = "D"),
      ncol = 2
    )
    plot_width <- width # 30
    plot_height <- height * 2 # 65
  } else {
    final_plot <- arrangeGrob(plot_annotation(p1, letter = "A"),
      plot_annotation(p2, letter = "B"),
      ncol = 1
    )
    plot_width <- width
    plot_height <- height * 2
  }

  lines_suffix <- ifelse(plot_lines, "_plot_lines", "")
  infection_samples_suffix <- ifelse(plot_n_infections, paste0(
    "_",
    plot_n_infections,
    "infections"
  ),
  ""
  )
  output_figure_filepath <- file.path(dir_figures, paste0(
    "ppc_timecourse_",
    stat, "_",
    model,
    lines_suffix,
    infection_samples_suffix,
    "_color_by_", color_by,
    ".png"
  ))
  ggsave(
    filename = output_figure_filepath, plot = final_plot,
    width = plot_width, height = plot_height, dpi = 300, bg = "white"
  )
}


#' Process posterior draws and make predictions.
#'
#' @param draws A `draws` object from the `posterior` package (containing 
#' values sampled in each draw, for each model parameter).
#' @param day_data A data.table with PCR test and metadata used for fitting.
#' @param model A string specifying the model name (includes numeric model 
#' identifiers, data filtering criteria, and run settings).
#' @param model_no_stan The number/identifier of the stan file used for fitting.
#' @param dir_pdata A string specifying the directory path where posterior 
#' predictions are stored.
#' @param day_var A string specifying the day variable used in the model 
#' (e.g. "day" with day 0 being the first day with a log10 viral load > 3; 
#' or "day2" with day 0 being the day of the highest viral load measured in the 
#' infection).
#' @param ss A data.table containing summary statistics (including Rhat and ESS 
#' values) from a specific run.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that will 
#' be computed on shift parameters.
#' @param testing A logical value. If `TRUE`, loads only a subset of posterior 
#' draws (for testing/debugging).
#'
#' @return A named list with:
#'  -`VLCP_by_day_draws`: A data.table with conditional posterior predictions of 
#'                        viral load trajectories
#'                       ,computed from infection-level parameter estimates.
#'  -`ppc_timecourse`: A ggplot object (with posterior predictions of viral 
#'                     load trajectories).
#'  -`VLCP_by_day_ID`: A data.table with infection-wise posterior predictions 
#'                     of median viral load trajectories.
#'  -`shift_draws`: A data.table with posterior draws of temporal shifts, 
#'                  grouped by infection.
#'  -`shifted_data_by_draw_day_ID`: A data.table with posterior draws of shifted 
#'                                  days (i.e. estimated days of sampling 
#'                                  relative to the day of peak viral load),
#'                                  and imputed log10 viral loads grouped by 
#'                                  infection.
make_post_pred_all <- function(draws,
                               day_data,
                               model,
                               model_no_stan,
                               dir_pdata,
                               day_var,
                               ss,
                               stat = "median",
                               testing = FALSE) {
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

  vlcp_filepath <- file.path(dir_pdata, paste0(
    "VLCP_by_day_draws_min",
    "_",
    ifelse(stat == "median", "median_", ""),
    model,
    ".Rdata"
  ))
  if (file.exists(vlcp_filepath)) {
    load(vlcp_filepath)
  } else {
    VLCP_by_day_draws <-
      make_VLCP_by_draw_ID(draws,
        day_data,
        days = tc_days,
        thin = thin
      )
    setkeyv(VLCP_by_day_draws, c("ID", ".draw"))
    if (stat == "median") {
      # calculate posterior expectation of viral load
      # by ID and day
      if ("alpha_CP" %in% variables(draws)) {
        VLCP_by_day_ID <-
          VLCP_by_day_draws %>%
          .[, list(
            CP = collapse::fmedian(CP),
            log10_load = collapse::fmedian(log10_load)
          ),
          by = c("ID", "day_shifted")
          ]
      } else {
        VLCP_by_day_ID <-
          VLCP_by_day_draws %>%
          .[, list(log10_load = collapse::fmedian(log10_load)),
            by = c("ID", "day_shifted")
          ]
      }
    } else {
      # calculate posterior expectation of viral load
      # by ID and day
      if ("alpha_CP" %in% variables(draws)) {
        VLCP_by_day_ID <-
          VLCP_by_day_draws %>%
          .[, list(
            CP = collapse::fmean(CP),
            log10_load = collapse::fmean(log10_load)
          ),
          by = c("ID", "day_shifted")
          ]
      } else {
        VLCP_by_day_ID <-
          VLCP_by_day_draws %>%
          .[, list(log10_load = collapse::fmean(log10_load)),
            by = c("ID", "day_shifted")
          ]
      }
    }

    # select sub-set of cases for plotting
    # Only sampling 20 infections here for smaller plot.
    ids <- sample(unique(VLCP_by_day_ID$ID), 
                  min(20, uniqueN(VLCP_by_day_ID$ID)))

    VL_by_dayID <-
      VLCP_by_day_draws[ID %in% ids] %>%
      summarise_draws_dt_by(
        by = c("day_shifted", "ID"),
        target.var = "log10_load",
        varname = "log10_load",
        stat = stat
      )
    save(VLCP_by_day_draws,
      VLCP_by_day_ID,
      VL_by_dayID,
      ids,
      file = vlcp_filepath
    )
  }

  # setkeyv(day_data,"ID") --> this leads to reordering which then messes with
  # later code that relies on day_data rows being in the original order
  tmp <- gc(verbose = FALSE)
  shift_var <- "shift"

  shift_draws <-
    draws_by_id(draws, day_data, c(shift_var), thin = thin)
  setkeyv(shift_draws, c("ID", ".draw"))

  day_data[log10_load <= VL_LOWER_LIMIT, 
           idx := 1:sum(day_data$log10_load <= VL_LOWER_LIMIT)]
  imputed_loads <-
    subset_draws(draws, "imp_neg") %>%
    thin_draws(thin = thin) %>%
    as_draws_dt() %>%
    melt(
      id.var = ".draw",
      variable.name = "par",
      value.name = "log10_load"
    ) %>%
    .[, idx := as.numeric(gsub("[^0-9]", "", par))] %>%
    .[, par := NULL]

  if (day_var == "day") {
    imputed_loads <- imputed_loads %>%
      merge(day_data[!is.na(idx), .(ID, day, idx)],
        by = "idx", all.x = TRUE, all.y = FALSE
      ) %>%
      .[, idx := NULL] %>%
      .[, imputed := "Yes"] %>%
      setkeyv("ID")
    day_data_to_merge <- day_data[log10_load > VL_LOWER_LIMIT, 
                                  .(ID, day, log10_load, material)]
  } else {
    imputed_loads <- imputed_loads %>%
      merge(day_data[!is.na(idx), .(ID, day, value = get(day_var), idx)],
        by = "idx", all.x = TRUE, all.y = FALSE
      ) %>%
      rename(!!day_var := value) %>%
      .[, idx := NULL] %>%
      .[, imputed := "Yes"] %>%
      setkeyv("ID")
    day_data_to_merge <- day_data[log10_load > VL_LOWER_LIMIT, 
                                  .(ID, 
                                    value = get(day_var), 
                                    day, 
                                    log10_load, 
                                    material)] %>% rename(!!day_var := value)
  }

  # assign imputed loads to correct ID and day
  # combine imputed and observed loads
  imputed_day_data_by_draw_ID <-
    shift_draws %>%
    .[, .(ID, .draw)] %>%
    unique() %>%
    merge(
      day_data_to_merge %>%
        .[, imputed := "No"],
      by = "ID",
      allow.cartesian = TRUE,
      all = TRUE
    ) %>%
    rbind(imputed_loads, fill = TRUE) %>%
    setkeyv(c("ID", ".draw", day_var))

  # 3 is the imputed maximum viral load (changed to 6 for now).
  assert_that(all(imputed_day_data_by_draw_ID[imputed == "Yes"]$log10_load <= 6))

  if (day_var == "day") {
    shifted_data_by_draw_day_ID <-
      shift_draws %>%
      merge(day_data[, .(day, ID, is.difficult, material)],
        by = "ID",
        allow.cartesian = TRUE
      ) %>%
      .[, day_shifted := day + get(shift_var)] %>%
      .[, c(shift_var) := NULL] %>%
      setkeyv(c("ID", ".draw", "day")) %>%
      .[imputed_day_data_by_draw_ID, log10_load := log10_load] %>%
      .[imputed_day_data_by_draw_ID, imputed := imputed]
  } else {
    shifted_data_by_draw_day_ID <-
      shift_draws %>%
      merge(day_data[, .(value = get(day_var), 
                         day, 
                         ID, 
                         is.difficult, 
                         material)],
        by = "ID",
        allow.cartesian = TRUE
      ) %>%
      rename(!!day_var := value) %>%
      .[, day_shifted := get(day_var) + get(shift_var)] %>%
      .[, c(shift_var) := NULL] %>%
      setkeyv(c("ID", ".draw", day_var)) %>%
      .[imputed_day_data_by_draw_ID, log10_load := log10_load] %>%
      .[imputed_day_data_by_draw_ID, imputed := imputed]
  }
  rm(imputed_loads)
  tmp <- gc(verbose = FALSE)

  draw_samples <- sample(max(shifted_data_by_draw_day_ID$.draw), 250)
  ppc_timecourse <-
    ggplot(VL_by_dayID[ID %in% ids], aes(x = day_shifted, y = log10_load)) +
    geom_ribbon(aes(ymin = q3, ymax = q97), fill = "blue", alpha = .3) +
    geom_line(col = "blue") +
    coord_cartesian(xlim = c(-10, 40), ylim = c(0, 10)) +
    facet_wrap(~ID, ncol = 5) +
    # geom_text(aes(x = 27, y = 8, label = ID)) +
    geom_jitter(
      data = shifted_data_by_draw_day_ID[ID %in% ids & .draw %in% draw_samples],
      height = .05, alpha = .01, color = "red",
      aes(x = day_shifted, y = log10_load)
    ) +
    geom_point(
      data = day_data[ID %in% ids], aes(x = get(day_var)),
      col = "black", pch = "x", size = 4
    ) +
    theme_thesis2() +
    xlab(DAYS_LABEL) +
    ylab(VL_LABEL)

  return(list(
    VLCP_by_day_draws = VLCP_by_day_draws,
    ppc_timecourse = ppc_timecourse,
    VLCP_by_day_ID = VLCP_by_day_ID,
    shift_draws = shift_draws,
    shifted_data_by_draw_day_ID = shifted_data_by_draw_day_ID
  ))
}
