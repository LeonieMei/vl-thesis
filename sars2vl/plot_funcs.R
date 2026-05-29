library(scales)
library(ggh4x)
library(ggstats)
library(here)

TOP <- here()
source(file.path(TOP, "sars2vl", "plot_params.R"))
source(file.path(TOP, "sars2vl", "data_funcs.R"))



#' Add manual scales with two colours red and blue, for fill and colour to 
#' ggplot plot.
#'
#' @return A `list` with manual scales for fill and color
red_blue <- function(o = 2:3) {
  clrs <- c("black", "red", "blue")
  list(
    scale_colour_manual(values = clrs[o]),
    scale_fill_manual(values = clrs[o])
  )
}


#' Adjust legend size in ggplots.
#'
#' @return A `list` with ggplot-theme for modified legend text and icons.
gg_legend_size <- function(ncol = 3) {
  return(
    list(
      theme(
        legend.text = element_text(size = 8), # 5),
        legend.title = element_text(size = 9), # 5),
        legend.key.size = unit(0.5, "cm")
      ),
      guides(
        fill = guide_legend(
          ncol = ncol,
          keywidth = .75,
          keyheight = .75
        ),
        color = guide_legend(
          ncol = ncol,
          keywidth = .75,
          keyheight = .75
        )
      )
    )
  )
}


#' Expand x and y axes in ggplot.
#'
#' @return A `list` with instructions to expand axes.
gg_expand <- function(x1 = 0, x2 = 0, y1 = .01, y2 = 0) {
  list(
    scale_x_continuous(expand = expansion(x1, x2)),
    scale_y_continuous(expand = expansion(y1, y2))
  )
}


#' Set the size of x- and y-axis tick labels in a ggplot object to 14.
#'
#' @param p A ggplot object.
#' @return The ggplot object with updated tick label sizes.
gg_set_tick_label_size <- function(p) {
  return(p + theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ))
}


#' Custom theme for PhD thesis plots.
#'
#' @param add Numeric. Adjusts the base font size for all text elements in the 
#' plot.
#' @return A ggplot theme object.
theme_thesis <- function(add = 0) {
  return(theme_minimal() + theme(
    axis.title = element_text(size = TEXT_SIZE + add),
    axis.text = element_text(size = TEXT_SIZE + add), # Size of axis tick labels
    legend.title = element_text(size = LEGEND_TITLE_SIZE + add), 
    legend.text = element_text(size = LEGEND_TEXT_SIZE + add),
    strip.text = element_text(size = TEXT_SIZE - 3 + add),
    plot.title = element_text(size = PLOT_TITLE_SIZE + add),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ))
}


#' Custom theme for PhD thesis plots (with facet strip backgrounds, text, and 
#' legends removed).
#'
#' @param add Numeric. Adjusts the base font size for all text elements in the 
#' plot.
#' @return A ggplot theme object.
theme_thesis2 <- function(add = 0) {
  return(theme_minimal() + theme(
    axis.title = element_text(size = TEXT_SIZE + add), 
    # Size of axis tick labels
    axis.text = element_text(size = TEXT_SIZE - 4 + add), 
    legend.title = element_text(size = LEGEND_TITLE_SIZE + add), 
    legend.text = element_text(size = LEGEND_TEXT_SIZE + add),
    plot.title = element_text(size = PLOT_TITLE_SIZE + add),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_blank(),
    legend.position = "none",
    strip.text = element_blank(), # Remove subplot titles
    strip.text.x = element_blank(),
    strip.text.x.top = element_blank(),
  ))
}


#' Customize a ggplot object with adjusted labels, legend, and axis visibility.
#'
#' Apply a theme for doctoral thesis document, adjust axis expansion, and 
#' control visibility of axis labels and legend. Optionally set a legend title 
#' and position.
#'
#' @param ggplot_object ggplot object. Input plot to modify.
#' @param legend_title Character. Title for the legend (NULL for no title).
#' @param hide_x_label Logical. If TRUE, hide the x-axis label.
#' @param hide_y_label Logical. If TRUE, hide the y-axis label.
#' @param add_size Numeric. Additional font size adjustment.
#' @param legend_position Character. Position of the legend ("right", "left", 
#' etc.).
#' @param hide_legend Logical. If TRUE, hide the legend entirely.
#'
#' @return ggplot object. Modified plot with updated theme and layout.
modified_plot <- function(ggplot_object,
                          legend_title = NULL,
                          hide_x_label = FALSE,
                          hide_y_label = FALSE,
                          add_size = 0,
                          legend_position = "right",
                          hide_legend = FALSE) {
  ggplot_object_new <- ggplot_object +
    # coord_cartesian(clip = "off", ylim = ylim, xlim = xlim) +
    theme_thesis(add_size) +
    # Expand limits to avoid clipping
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.02)),
      position = "bottom"
    ) + # move x-axis to top
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.02)),
      position = "left"
    ) # move y-axis to right

  if (hide_x_label) {
    ggplot_object_new <- ggplot_object_new + theme(axis.title.x = element_blank())
  }
  if (hide_y_label) {
    ggplot_object_new <- ggplot_object_new + theme(axis.title.y = element_blank())
  }
  if (hide_legend) {
    ggplot_object_new <- ggplot_object_new + theme(legend.position = "none")
  } else {
    ggplot_object_new <- ggplot_object_new + theme(legend.position = legend_position)
  }

  if (is.null(legend_title)) {
    return(ggplot_object_new + theme(legend.title = element_blank()))
  } else {
    return(ggplot_object_new + guides(
      fill = guide_legend(title = legend_title),
      color = guide_legend(title = legend_title)
    ))
  }
}


#' Theme for shrinking the legend in a ggplot object.
#' @return A ggplot theme object.
shrink_legend <- function() {
  return(
    theme(
      legend.key.size = unit(0.3, "cm"), # change legend key size
      legend.key.height = unit(0.3, "cm"), # change legend key height
      legend.key.width = unit(0.3, "cm"), # change legend key width
      # change legend title font size
      legend.title = element_text(size = LEGEND_TITLE_SIZE - 2),
      legend.text = element_text(size = LEGEND_TEXT_SIZE - 2)
    ) # change legend text font size
  )
}


#' Add annotation tag to a ggplot object.
#'
#' @param plot A ggplot2 plot object.
#' @param letter The annotation tag (e.g., "A", "B").
#' @param size Font size (numeric) of the annotation tag.
#' @param hjust Horizontal adjustment (numeric) for the tag position.
#' @param left_margin Left margin (numeric) of the plot.
#'
#' @return A ggplot2 plot with an annotation tag added.
plot_annotation <- function(plot, letter = "A", size = 14, hjust = 1.7, 
                           left_margin = 17) {
  return(plot + labs(tag = letter) +
    theme(
      plot.tag = element_text(
        size = size,
        hjust = hjust,
        vjust = 0
      ),
      plot.tag.position = c(0, 0.98),
      plot.margin = margin(l = left_margin, t = 10, r = left_margin)
    ))
}


#' Plot grouped viral load data as beeswarm.
#'
#' Creates a faceted beeswarm plot of log10 viral load by group.
#'
#' @param bdata A data.table
#' @param group_var Character. The name of the grouping variable.
#' @param x_var Character. The name of the x-axis variable.
#' @param xlabel Character. The label for the x-axis.
#'
#' @return A ggplot object with a faceted beeswarm plot.
bees_group_func <- function(bdata, group_var, x_var, xlabel) {
  bdata_var <- copy(bdata)
  bdata_var %>%
    .[, group := eval(parse(text = group_var))]

  bees_group_var <-
    bdata_var %>%
    ggplot(aes(x = get(x_var), y = log10_load, color = get(group_var))) +
    geom_quasirandom(method = "tukeyDense", size = .1) +
    ylab(expression(Log[10] ~ viral ~ load)) +
    facet_grid(group ~ .) +
    gg_add_grid("y") +
    theme(strip.text.x = element_blank(), legend.position = "none") +
    theme(axis.text.x = element_text(size = 10)) +
    xlab(xlabel)

  return(bees_group_var)
}


#' Annotate a plot with estimated (or true) viral load trajectories, 
#' indicating estimated (or true) peak viral load, proliferation time, and 
#' clearance time.
#'
#' @param stat_params A named list with estimated or true (in the case of 
#' simulated data) parameter values.
#' @param p A ggplot object.
#' @param stat The summary statistic (mean, median, 94% HDI, 50% HDI) that was 
#' computed for posterior estimates and predictions.
#' @param xpos_label The x-coordinate for the legend's position.
#' @param ypos_label The x-coordinate for the legend's position.
#' @param sim_data A named list, containing a data.table with simulated data 
#' used for model fitting, and true (simulated) parameter values.
#'
#' @return The annotated ggplot object.
annotate_summary_stat_plot <- function(stat_params,
                                       p,
                                       stat = "median",
                                       xpos_label = 23,
                                       ypos_label = 11,
                                       sim_data = FALSE) {
  stat_str <- ifelse(stat == "median", "Median ", "Mean ")
  label_text <- paste(paste0(stat_str, " peak viral load: ", 
                             format(round(stat_params$intercept, 1), 
                                    nsmall = 1)),
    paste0(stat_str, " proliferation time: ", 
           format(round(stat_params$time2peak, 1), nsmall = 1), " days"),
    paste0(stat_str, " clearance time: ", 
           format(round(stat_params$time_from_peak, 1), nsmall = 1), " days"),
    sep = "\n"
  )

  if (sim_data) {
    # Decrease label size
    return(p + geom_text(aes(x = xpos_label, y = ypos_label), 
                         label = label_text, size = 3))
  } else {
    return(p + geom_text(aes(x = xpos_label, y = ypos_label), 
                         label = label_text))
  }
}


#' Create and save plots showing differences between estimated and true values 
#' for runs with simulated data.
#'
#' @param param_sizes_all_zero A logical value indicating whether to make plots 
#' for simulation runs where all fixed effects (infection-group-specific 
#' parameters) were set to 0 or only those where these values were different 
#' from 0.
generate_sim_data_group_plots <- function(param_sizes_all_zero) {
  # Current implementation only for exponentiated values
  log_prefix <- ""
  nrow_legend <- ifelse(param_sizes_all_zero, 1, 2)

  if (param_sizes_all_zero) {
    param <- "hdi_contains_true_value"
    column_suffix <- "_hdi_contains_true_value"
    models <- c(1, 2, 3)
    plot_dir <- file.path(DIR_FIGURES_VL_SIM, "param_sizes_all_zero")
    shapes <- c(15, 16) # shapes in order
    shape_label <- "94% HPDI includes true value"
    filename_suffix <- ""
    contains_true_value_cols <- TRUE
    combination_cols <- FALSE
    shape_label_mapping <- c(
      "TRUE" = "True",
      "FALSE" = "False"
    )
  } else {
    param <- "hdi_combination"
    column_suffix <- "_hdi_combination"
    models <- c(1, 2)
    plot_dir <- file.path(DIR_FIGURES_VL_SIM, "param_sizes_not_zero")
    shapes <- c(15, 16, 17, 18)
    shape_label <- "94% HPDI"
    filename_suffix <- "_combination"
    contains_true_value_cols <- FALSE
    combination_cols <- TRUE
    shape_label_mapping <- c(
      "1" = "includes true value,\nexcludes 1",
      "2" = "includes true value,\nincludes 1",
      "3" = "excludes true value,\nincludes 1",
      "4" = "excludes true value,\nexcludes 1"
    )
  }

  group_cols <- return_group_cols(log_prefix)
  group_cols_hdi <- return_group_cols(log_prefix, hdi = TRUE)
  group_cols_hdi50 <- return_group_cols(log_prefix, hdi50 = TRUE)
  group_cols_mapping2 <- return_group_cols_hdi_mapping(log_prefix, 
                                                       column_suffix = column_suffix)


  name_mapping <- list(
    alpha = "Alpha",
    delta = "Delta",
    omicron = "Omicron",
    gender = "Male",
    hospitalized = "Hospitalized",
    pams = "PAMS",
    age_cat2 = "30-60 years old",
    age_cat3 = ">60 years old",
    prior_infection = "Prior infection",
    hospitalized_alpha = "Hospitalized,\nAlpha",
    hospitalized_delta = "Hospitalized,\nDelta",
    hospitalized_omicron = "Hospitalized,\nOmicron"
  )
  true_vals_intercept <- list(
    alpha = exp(0),
    delta = exp(0.1),
    omicron = exp(0.15),
    gender = exp(0),
    hospitalized = exp(0.05),
    pams = exp(-0.05),
    age_cat2 = exp(0.05),
    age_cat3 = exp(0.1),
    prior_infection = exp(-0.05),
    hospitalized_alpha = exp(0),
    hospitalized_delta = exp(0),
    hospitalized_omicron = exp(0)
  )
  true_vals_slope_up <- list(
    alpha = exp(0.2),
    delta = exp(0.3),
    omicron = exp(0.4),
    gender = exp(0),
    hospitalized = exp(0),
    pams = exp(0),
    age_cat2 = exp(0),
    age_cat3 = exp(0),
    prior_infection = exp(-0.05),
    hospitalized_alpha = exp(0),
    hospitalized_delta = exp(0),
    hospitalized_omicron = exp(0)
  )
  true_vals_slope_down <- list(
    alpha = exp(0.1),
    delta = exp(0.3),
    omicron = exp(0.4),
    hospitalized = exp(-0.1),
    gender = exp(0),
    pams = exp(0.1),
    age_cat2 = exp(0.025),
    age_cat3 = exp(-0.15),
    prior_infection = exp(0.2),
    hospitalized_alpha = exp(0),
    hospitalized_delta = exp(0),
    hospitalized_omicron = exp(0)
  )

  for (model_no in models) {
    dts_model <- get_sim_dts(model_no = model_no, 
                             param_sizes_all_zero = param_sizes_all_zero)
    dt_current <- dts_model[["dte"]]

    if ((model_no == 2) && (param_sizes_all_zero == TRUE)) {
      id_vars <- c("model_no_stan", "min_n_tests", "iteration")
    } else {
      id_vars <- c("model_no_stan", "min_n_tests")
    }
    id_vars_new <- c(id_vars, "min_n_pos_tests")
    xlim_mapping <- list(
      "intercept" = c(0.75, 1.28), "slope_up" = c(0.35, 2.15),
      "slope_down" = c(0.75, 1.9)
    )
    by_mapping <- list("intercept" = 0.2, "slope_up" = 0.4, "slope_down" = 0.2)
    for (statistic_current in c("intercept", "slope_up", "slope_down")) {
      x_lims <- xlim_mapping[[statistic_current]]
      by <- by_mapping[[statistic_current]]

      group_cols_current <- group_cols$group_cols_mapping[[statistic_current]]
      group_cols_hdi_current <- group_cols_hdi$group_cols_mapping[[statistic_current]]
      group_cols_hdi50_current <- group_cols_hdi50$group_cols_mapping[[statistic_current]]
      group_cols_current2 <- group_cols_mapping2[[statistic_current]]
      cols_reshaped <- c(id_vars_new, group_cols_hdi_current, 
                         group_cols_hdi50_current, group_cols_current2)

      if (param_sizes_all_zero) {
        true_data <- data.table(
          group = unname(unlist(name_mapping)),
          true_value = rep(exp(0), length(name_mapping))
        )
      } else {
        if (statistic_current == "intercept") {
          true_data <- data.table(
            statistic_label = factor(unname(unlist(name_mapping)),
              levels = unname(unlist(name_mapping))
            ),
            true_value = as.numeric(true_vals_intercept[names(name_mapping)])
          )
        } else if (statistic_current == "slope_up") {
          true_data <- data.table(
            statistic_label = factor(unname(unlist(name_mapping)), 
                                     levels = unname(unlist(name_mapping))),
            true_value = as.numeric(true_vals_slope_up[names(name_mapping)])
          )
        } else {
          true_data <- data.table(
            statistic_label = factor(unname(unlist(name_mapping)), 
                                     levels = unname(unlist(name_mapping))),
            true_value = as.numeric(true_vals_slope_down[names(name_mapping)])
          )
        }
      }

      name_mapping_new <- setNames(
        unname(unlist(name_mapping[group_cols$group_cols_base])),
        group_cols_current
      )
      dtr_reshaped <- reshape_dt(dt_current[, ..cols_reshaped],
        id_vars = id_vars_new,
        contains_true_value_cols = contains_true_value_cols,
        combination_cols = combination_cols
      )
      dtr_reshaped[, n_tests_combination := paste0(min_n_tests, 
                                                   "_", 
                                                   min_n_pos_tests)]
      dtr_reshaped[, statistic_label := factor(name_mapping_new[statistic],
        levels = unname(unlist(name_mapping))
      )]

      dodge_height <- 0.5
      p <- ggplot(dtr_reshaped, aes(
        y = min_n_tests, x = value, color = min_n_pos_tests,
        shape = !!sym(param), group = min_n_pos_tests
      )) +
        geom_errorbar(aes(xmin = lower, xmax = upper, group = min_n_pos_tests),
          position = ggstance::position_dodgev(height = dodge_height),
          width = 0,
          linewidth = 0.2
        ) +
        geom_errorbar(aes(xmin = lower50, xmax = upper50, 
                          group = min_n_pos_tests),
          position = ggstance::position_dodgev(height = dodge_height),
          width = 0,
          linewidth = 0.45
        ) +
        geom_point(
          position = ggstance::position_dodgev(height = dodge_height), 
          alpha = ALPHA_SIM_PLOTS,
          size = DATA_POINT_SIZE_SIM_PLOTS - 2
        ) +
        scale_color_viridis_d(guide = guide_legend(reverse = TRUE)) +
        geom_stripped_rows(color = NA) +
        geom_vline(xintercept = 1, color = "darkgray", linetype = "solid", 
                   linewidth = 0.4) +
        geom_vline(data = true_data, aes(xintercept = true_value), 
                   color = "darkred", linewidth = 0.2) +
        scale_x_continuous(
          limits = x_lims,
          breaks = seq(x_lims[1] + 0.05, x_lims[2], by = by)
        ) +
        labs(
          title = "",
          x = paste0("Parameter estimate"),
          shape = shape_label,
          y = MINIMUM_PCRS_LABEL,
          color = MINIMUM_POS_PCRS_LABEL
        ) +
        scale_shape_manual(
          values = shapes,
          labels = shape_label_mapping
        ) +
        # ylim(ylim_plot[1], ylim_plot[2]) +
        theme_thesis() +
        theme(
          # Remove major horizontal gridlines
          panel.grid.major.y = element_blank(),
          # Remove minor horizontal gridlines
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          legend.title = element_text(hjust = 0.5), # Center the title
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.box.just = "top",
          # Remove top margin inside legend
          legend.margin = margin(t = 0, unit = "pt"), 
          # Top, right, bottom, left margins
          plot.margin = margin(1, 4, 1, 1, "pt")
        ) +
        guides(
          color = guide_legend(nrow = 1, byrow = FALSE),
          shape = guide_legend(nrow = nrow_legend, byrow = FALSE)
        ) +
        facet_wrap(. ~ statistic_label, ncol = 4)

      ggsave(
        filename = file.path(plot_dir, paste0(
          "sim", str_to_sentence(statistic_current),
          "GroupDiffSummaryModel", filename_suffix,
          model_no, ".png"
        )),
        plot = p,
        width = COL_WIDTH * 2,
        height = HEIGHT_SIM_PLOT * 2.7,
        units = "in", # Units: inches ("in"), centimeters ("cm"), etc.
        dpi = DPI,
        limitsize = FALSE
      )
    }
  }
}


#' Create and save plots of posterior estimates (medians and 94% HPDIs) of 
#' infection-group-specific parameters/fixed effects (separately for 
#' proliferation time, peak, clearance rate). Estimates from Models 1 and 2 
#' (with thresholds of 3 and 5 PCR days) are plotted.
#'
#' @param dt A data.table containing posterior estimates of fixed effects.
#' @param plot_dir A string specifying the directory path where figures will be 
#' saved.
#' @param height_plot A numeric specifying the plot height (in inches).
#' @param sampling_iterations A numeric value specifying the number of sampling 
#' iterations run during model fitting. Show results only for the runs with the 
#' corresponding number of sampling iterations.
generate_group_plots <- function(dt,
                                 plot_dir,
                                 height_plot,
                                 sampling_iterations = 1000,
                                 model_nos = c(1, 2)) {
  assert_that(identical(model_nos, c(1, 2)) | identical(model_nos, c(4, 5)))
  dt_models <- dt[model_no_stan %in% model_nos]
  group_cols <- return_group_cols(log_prefix = "")$group_cols_mapping
  group_cols_labels <- return_group_cols(log_prefix = "")$group_cols_mapping_labels

  # Define patterns and corresponding groups
  patterns <- list(
    "SARS-CoV-2 variant" = "(?<!hospitalized_)wildtype|alpha|delta|omicron",
    "PAMS" = "pams",
    "Hospitalized" = "hospitalized$",
    "Prior infection" = "prior_infection",
    "Gender" = "female|male",
    "Age" = "age_cat",
    "Hospitalized, SARS-CoV-2 variant" = "(?<=hospitalized_)wildtype|alpha|delta|omicron"
  )

  # Function to find the first matching pattern
  find_group <- function(value) {
    for (group in names(patterns)) {
      if (grepl(patterns[[group]], value, perl = TRUE)) {
        return(group)
      }
    }
    return("unknown")
  }
  labels <- rep("Parameter estimate", 3)
  counter <- 1
  for (main_param in c("slope_up", "intercept", "slope_down")) {
    cols_current <- unlist(unname(group_cols[[main_param]]))
    cols_hdi_current <- paste0(cols_current, "_hdi")
    cols_hdi50_current <- paste0(cols_current, "_hdi50")
    cols_reshape <- c(
      "model_no_stan", "min_n_tests", cols_current,
      cols_hdi_current, cols_hdi50_current
    )
    name_mapping_current <- group_cols_labels[[main_param]]

    model_levels <- model_nos
    model_labels <- c("Model 1", "Model 2")
    model_n_levels <- rev(c("Model 1 3", "Model 1 6", "Model 2 3", "Model 2 6"))
    model_n_labels <- rev(c(
      "Model 1, 3 PCR days", "Model 1, 6 PCR days", "Model 2, 3 PCR days",
      "Model 2, 6 PCR days"
    ))

    dt_reshaped <- reshape_dt(dt_models[, ..cols_reshape], hdi_cols = TRUE)
    dt_reshaped[, model := factor(model_no_stan, levels = model_levels, 
                                  labels = model_labels)]
    dt_reshaped[, model_n := factor(paste(model, min_n_tests, sep = " "),
      levels = model_n_levels,
      labels = model_n_labels
    )]
    # Create the new column
    dt_reshaped[, group := factor(sapply(statistic, find_group),
      levels = c(
        "SARS-CoV-2 variant", "PAMS",
        "Hospitalized", "Gender", "Age",
        "Prior infection",
        "Hospitalized, SARS-CoV-2 variant"
      )
    )]
    dt_reshaped[, statistic := factor(statistic,
      levels = rev(cols_current),
      labels = rev(unname(name_mapping_current))
    )]
    dt_reshaped[, model_n_combination := factor(paste0(model_no_stan, 
                                                       min_n_tests))]

    dodge_height <- 0.7
    p <- ggplot(dt_reshaped, aes(
      x = value, y = statistic,
      color = model_n,
      group = model_n
    )) +
      geom_vline(xintercept = 1, color = "darkgray", linetype = "solid", 
                 linewidth = 0.4) +
      scale_color_manual(
        values = c(
          "Model 1, 3 PCR days" = "#FF6B6B",
          "Model 1, 6 PCR days" = "#C91818",
          "Model 2, 3 PCR days" = "#4ECDC4",
          "Model 2, 6 PCR days" = "#1E9AA8"
        ),
      ) +
      geom_stripped_rows(color = NA) +
      geom_errorbar(aes(xmin = lower, xmax = upper, group = model_n),
        position = ggstance::position_dodgev(height = dodge_height),
        width = 0,
        linewidth = 0.15
      ) +
      geom_errorbar(aes(xmin = lower50, xmax = upper50, group = model_n),
        position = ggstance::position_dodgev(height = dodge_height),
        width = 0,
        linewidth = 0.45
      ) +
      geom_point(
        position = ggstance::position_dodgev(height = dodge_height),
        alpha = ALPHA_SIM_PLOTS,
        size = DATA_POINT_SIZE_SIM_PLOTS - 1.7
      ) +
      labs(x = labels[counter], y = "") +
      theme_thesis() +
      theme(
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),# Remove major horizontal gridlines
        panel.grid.minor.y = element_blank(),# Remove minor horizontal gridlines
        legend.position = "top",
        plot.margin = margin(1, 4, 1, 1, "pt")# Top, right, bottom, left margins
      ) +
      guides(color = guide_legend(nrow = 2, byrow = FALSE, reverse = TRUE))

    ggsave(
      filename = file.path(plot_dir, paste0(
        "group_", main_param, "_iter",
        sampling_iterations, ".png"
      )),
      plot = p,
      width = COL_WIDTH * 2,
      height = height_plot * 2.3,
      units = "in", # Units: inches ("in"), centimeters ("cm"), etc.
      dpi = DPI
    )
    counter <- counter + 1
  }
}

#' Create and save plot of posterior predictions (medians and 94% HPDIs) for 
#' population-wide proliferation time, peak, and clearance rates. Predictions 
#' from Models 1 and 2 (with thresholds of 3 and 5 PCR days) are plotted.
#'
#' @param plot_dir A string specifying the directory path where figures will be 
#' saved.
#' @param sampling_iterations A numeric value specifying the number of sampling 
#' iterations run during model fitting. Show results only for the runs with the 
#' corresponding number of sampling iterations.
generate_main_params_plots <- function(plot_dir, sampling_iterations = 1000) {
  # Current implementation only for exponentiated values
  log_prefix <- ""
  models <- c(1, 2)
  model_levels <- c(1, 2)
  model_labels <- c("Model 1", "Model 2")
  model_nos <- c(1, 2)

  dts <- get_dts(models, sampling_iterations = sampling_iterations)
  # Use exponentiated as opposed to log values
  dt_current <- dts[["dte"]]
  id_vars <- c("model_no_stan", "min_n_tests")

  cols <- c("intercept", "time2peak", "slope_down")
  hdi_cols <- paste0(cols, "_hdi")
  hdi50_cols <- paste0(cols, "_hdi50")
  cols_reshaped <- c(id_vars, cols, hdi_cols, hdi50_cols)

  name_mapping <- c(
    time2peak = PROLIFERATION_LABEL2,
    intercept = PEAK_LABEL2,
    slope_down = DOWN_SLOPE_LABEL2_2
  )
  main_params <- c("time2peak", "intercept", "slope_down")
  labels <- c(PROLIFERATION_LABEL2, PEAK_LABEL2, DOWN_SLOPE_LABEL2_2)

  dt_reshaped <- reshape_dt(dt_current[, ..cols_reshaped],
    id_vars = id_vars,
    hdi_cols = TRUE
  )
  dt_reshaped[, model := factor(model_no_stan, levels = rev(model_levels), 
                                labels = rev(model_labels))]
  dt_reshaped[, statistic := factor(statistic, levels = names(name_mapping))]
  dt_reshaped[, statistic_label := factor(statistic,
    levels = main_params,
    labels = labels
  )]

  dodge_height <- 0.5
  p <- ggplot(dt_reshaped, aes(
    y = min_n_tests, x = value, color = model,
    group = model
  )) +
    geom_errorbar(aes(xmin = lower, xmax = upper, group = model),
      position = ggstance::position_dodgev(height = dodge_height),
      width = 0,
      linewidth = 0.15
    ) +
    geom_errorbar(aes(xmin = lower50, xmax = upper50, group = model),
      position = ggstance::position_dodgev(height = dodge_height),
      width = 0,
      linewidth = 0.45
    ) +
    geom_point(
      position = ggstance::position_dodgev(height = dodge_height),
      alpha = ALPHA_SIM_PLOTS,
      size = DATA_POINT_SIZE_SIM_PLOTS - 1.7
    ) +
    geom_stripped_rows(color = NA) +
    scale_color_manual(
      values = c(
        "Model 1" = "#FF6B6B",
        "Model 2" = "#4ECDC4"
      ),
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(
      title = "",
      y = MINIMUM_PCRS_LABEL,
      x = "",
      color = "Model"
    ) +
    theme_thesis() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_blank(), # Remove major horizontal gridlines
      panel.grid.minor.y = element_blank(), # Remove minor horizontal gridlines
      legend.title = element_blank(),
      strip.text = element_text(size = 10)
    ) +
    facet_wrap(. ~ statistic_label,
      label = label_parsed,
      ncol = 3, scales = "free_x"
    )

  # Specify different breaks for each facet
  p <- p +
    facetted_pos_scales(
      x = list(
        # Breaks for facet A
        scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)), 
        # Breaks for facet B
        scale_x_continuous(breaks = c(8.8, 8.9, 9.0)), 
        # Breaks for facet C
        scale_x_continuous(breaks = c(-0.25, -0.24, -0.23, -0.22)) 
      )
    )

  ggsave(
    filename = file.path(plot_dir, paste0(
      "paramsSummaryModelsIter",
      sampling_iterations, ".png"
    )),
    plot = p,
    width = COL_WIDTH * 2,
    height = HEIGHT_SIM_PLOT,
    units = "in", # Units: inches ("in"), centimeters ("cm"), etc.
    dpi = DPI
  )
}


#' Create and save an overview plot of R-hat values for selected models and 
#' parameters.
#'
#' @param dtrhat_reshaped A data.table containing Rhat values.
#' @param model_nos A vector of model numbers to include (refers to stan file 
#' numbers).
#' @param scales_plot Scaling for facet axes (e.g., "free", "fixed").
#' @param labeller_model Labeller function. Custom labels for model numbers.
#' @param filename_suffix Suffix for the output filename.
#' @param filename_suffix2 Additional suffix for the output filename.
#' @param height_plot Height of the plot in inches.
#' @param color_by_min_tests A logical value indicating whether to color 
#' histogram by minimum number of PCR days specified for a run.
#' @param sampling_iterations A numeric value specifying the number of sampling 
#' iterations run during model fitting. Show results only for the runs with the 
#' corresponding number of sampling iterations.
#'
#' @return Save a histogram plot of R-hat values, with facets for each model 
#' and main parameter (proliferation time, peak, clearance rate).
make_rhat_plot <- function(dtrhat_reshaped, model_nos, scales_plot,
                           labeller_model, filename_suffix, filename_suffix2,
                           height_plot, color_by_min_tests = TRUE,
                           sampling_iterations = 1000) {
  assert_that(sampling_iterations %in% c(1000, 2000))
  dtrhat_plot <- dtrhat_reshaped[iter_sampling == eval(sampling_iterations)]
  if (color_by_min_tests) {
    p_rhat_hist <- ggplot(
      data = dtrhat_plot[(model_no_stan %in% model_nos) &
        (main_param %in% c(
          "Peak",
          "Proliferation rate",
          "Clearance rate"
        ))],
      aes(x = value, fill = min_tests_factor)
    ) +
      scale_fill_viridis_d()
  } else {
    p_rhat_hist <- ggplot(
      data = dtrhat_plot[(model_no_stan %in% model_nos) &
        (main_param %in% c(
          "Peak",
          "Proliferation rate",
          "Clearance rate"
        ))],
      aes(x = value)
    )
  }

  p_rhat_hist <- p_rhat_hist +
    geom_histogram(position = "stack", bins = 30) +
    # coord_cartesian(xlim = c(min_val, max_val)) +
    geom_segment(
      aes(
        x = 1.01,
        xend = 1.01,
        y = 0,
        yend = Inf
      ),
      color = "black",
      linewidth = 0.2
    ) +
    geom_segment(
      aes(
        x = 1.1,
        xend = 1.1,
        y = 0,
        yend = Inf
      ),
      color = "red",
      linewidth = 0.2
    ) +
    facet_wrap(~ model_no_plot + main_param,
      ncol = 3,
      scales = scales_plot, labeller = labeller_model
    ) +
    labs(
      fill = MINIMUM_PCRS_LABEL2,
      x = expression(hat(R)), # hat on top of R
      y = "Count"
    ) +
    theme_thesis() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = TEXT_SIZE - 2),
      legend.position = "top"
    )

  ggsave(
    filename = file.path(plot_dir, paste0(
      "rhatsHistMainParamModel",
      filename_suffix, filename_suffix2,
      ".png"
    )),
    plot = p_rhat_hist,
    width = COL_WIDTH * 2,
    height = height_plot * 2,
    units = "in", # Units: inches ("in"), centimeters ("cm"), etc.
    dpi = DPI
  )
}


#' Create and save an overview plot of R-hat values for selected models and 
#' population-level or infection-level parameters.
#'
#' @param dtrhat_reshaped A data.table containing Rhat values.
#' @param plot_dir A string specifying the directory path where figures will be 
#' saved.
#' @param model_nos A vector of model numbers to include (refers to stan file 
#' numbers).
#' @param sim_data A logical value indicating whether values are from runs with 
#' simulated data.
#' @param height_plot Height of the plot in inches.
#' @param param_sizes_all_zero A logical value indicating whether to extract 
#' only values from simulation runs where all fixed effects 
#' (infection-group-specific parameters) were set to 0 or only those where 
#' these values were different from 0.
#' @param infection_level A logical value indicating whether Rhat values from 
#' infection-level parameters (as opposed to population-level parameters) 
#' should be plotted.
#' @param sampling_iterations A numeric value specifying the number of sampling 
#' iterations run during model fitting. Show results only for the runs with the 
#' corresponding number of sampling iterations.
#'
#' @return Save a histogram plot of R-hat values, with facets for each model 
#' and main parameter (proliferation time, peak, clearance rate), either for 
#' population-level or infection-level parameters.
make_rhat_overview_plot <- function(dtrhat_reshaped,
                                    plot_dir,
                                    model_nos,
                                    sim_data = TRUE,
                                    height_plot = HEIGHT_SIM_PLOT,
                                    param_sizes_all_zero = TRUE,
                                    infection_level = FALSE,
                                    sampling_iterations = 1000) {
  assert_that(sampling_iterations %in% c(1000, 2000))
  filename_suffix2 <- ifelse(infection_level, "_infection_level", "")
  filename_suffix3 <- ifelse(sampling_iterations == 1000, "", "_iter2000")
  filename_suffix2 <- paste0(filename_suffix2, filename_suffix3)
  for (model_no in model_nos) {
    assertthat::assert_that(model_no %in% c(1, 2, 3))
  }
  if (sim_data) {
    if (param_sizes_all_zero) {
      filename_suffix <- ""
      scales_plot <- "free_y"
    } else {
      filename_suffix <- "GroupDiffs"
      scales_plot <- "fixed"
    }
  } else {
    filename_suffix <- ""
    scales_plot <- "fixed"
  }

  max_val <- max(dtrhat_reshaped$value) + 0.01
  min_val <- min(dtrhat_reshaped$value) - 0.01


  labeller_model <- as_labeller(
    c(
      "Model 1" = "Model~1", "Model 2" = "Model~2",
      "Peak" = "Peak", "Proliferation rate" = "paste('Proliferation rate')",
      "Clearance rate" = "paste('Clearance rate')",
      "Model 1_2" = "Model~1[2]"
    ),
    default = label_parsed
  )

  if (infection_level) {
    min_tests <- c(3, 6)
    for (min_test in min_tests) {
      dt_min_tests <- dtrhat_reshaped[min_tests == eval(min_test)]
      filename_suffix3 <- paste0("_min", min_test, "_tests")
      make_rhat_plot(
        dtrhat_reshaped = dt_min_tests,
        model_nos = model_nos,
        scales_plot = scales_plot,
        labeller_model = labeller_model,
        filename_suffix = filename_suffix,
        filename_suffix2 = paste0(filename_suffix2, filename_suffix3),
        height_plot = height_plot,
        color_by_min_tests = FALSE,
        sampling_iterations = sampling_iterations
      )
    }
  } else {
    make_rhat_plot(
      dtrhat_reshaped = dtrhat_reshaped,
      model_nos = model_nos,
      scales_plot = scales_plot,
      labeller_model = labeller_model,
      filename_suffix = filename_suffix,
      filename_suffix2 = filename_suffix2,
      height_plot = height_plot,
      color_by_min_tests = TRUE,
      sampling_iterations = sampling_iterations
    )
  }
}


#' Create and save separate plots of R-hat values for each selected model 
#' showing values for population-level or infection-level parameters.
#'
#' @param dtrhat_reshaped A data.table containing Rhat values.
#' @param model_nos A vector of model numbers to include (refers to stan file 
#' numbers).
#' @param plot_dir A string specifying the directory path where figures will be 
#' saved.
#' @param infection_level A logical value indicating whether Rhat values from 
#' infection-level parameters (as opposed to population-level parameters) 
#' should be plotted.
#'
#' @return Save separate plots of R-hat values for each selected model showing 
#' values for population-level or infection-level parameters.
make_rhat_model_plots <- function(dtrhat_reshaped,
                                  model_nos,
                                  plot_dir,
                                  infection_level = FALSE,
                                  sampling_iterations = 1000) {
  assert_that(sampling_iterations %in% c(1000, 2000))
  filename_suffix <- ifelse(infection_level, "_inf_level", "")
  filename_suffix2 <- ifelse(sampling_iterations == 1000, "", "_iter2000")
  for (model_no_current in model_nos) {
    assert_that(model_no_current %in% c(1, 2, 3))
    p_rhat_hist <- ggplot(
      data = dtrhat_reshaped[(model_no_stan == model_no_current) &
        (iter_sampling == sampling_iterations)],
      aes(x = value, fill = min_n_pos_tests_factor)
    ) +
      geom_histogram(position = "stack", bins = 30) +
      scale_fill_viridis_d() + # Default: viridis palette
      # coord_cartesian(xlim = c(min_val, max_val)) +
      geom_segment(
        aes(
          x = 1.01,
          xend = 1.01,
          y = 0,
          yend = Inf
        ),
        color = "black",
        linewidth = 0.3
      ) +
      geom_segment(
        aes(
          x = 1.1,
          xend = 1.1,
          y = 0,
          yend = Inf
        ),
        color = "red",
        linewidth = 0.3
      ) +
      facet_wrap(~min_tests,
        ncol = 2,
        scales = "free_y"
      ) +
      labs(
        fill = MINIMUM_POS_PCRS_LABEL,
        x = expression(hat(R)), # hat on top of R
        y = "Count",
      ) +
      ggtitle(MINIMUM_PCRS_LABEL) +
      theme_thesis() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = TEXT_SIZE - 2),
        legend.position = "bottom"
      ) +
      guides(fill = guide_legend(nrow = 1))

    ggsave(
      filename = file.path(plot_dir, paste0(
        "rhatsHistModel", model_no_current,
        filename_suffix, filename_suffix2,
        ".png"
      )),
      plot = p_rhat_hist,
      width = COL_WIDTH * 2,
      height = HEIGHT_RHAT_PLOTS,
      units = "in", # Units: inches ("in"), centimeters ("cm"), etc.
      dpi = DPI
    )
  }
}


#' Create and save plots showing posterior predictions (medians and 94% HPDIs) 
#' of proliferation time, peak and clearance rate for different infection 
#' subgroups (computed from infection-level parameter estimates).
#'
#' @param dt A data.table containing conditional posterior predictions 
#' (computed from infection-level parameter estimates).
#' @param plot_dir A string specifying the directory path where figures will be 
#' saved.
#' @param height_plot A numeric specifying the plot height (in inches).
#' @param diff A logical value indicating whether estimates are given as 
#' differences to the respective reference group:
#' "non-PAMS, non-hospitalized" for disease severity status,
#' "pre-VOC" for SARS-CoV-2 variant,
#' "<30 years" for age group,
#' "no prior infection" for immunization status,
#' "female" for gender.
#' @param sampling_iterations A numeric value specifying the number of sampling 
#' iterations run during model fitting. Show results only for the runs with the 
#' corresponding number of sampling iterations.
#'
#' @return Save plots showing posterior predictions (medians and 94% HPDIs) for
#' different infection subgroups, computed from infection-level parameter 
#' estimates.
plot_group_averages <- function(dt,
                                plot_dir,
                                height_plot,
                                diff = FALSE,
                                sampling_iterations = 1000) {
  # Group averages
  group_cols_avg <- return_group_avg_cols(log_prefix = "", diff = diff)$group_cols_mapping
  group_cols_avg_labels <- return_group_avg_cols(log_prefix = "", diff = diff)$group_cols_mapping_labels
  group_cols_avg_labels_latex <- return_group_avg_cols(
    log_prefix = "", latex = TRUE,
    diff = diff
  )$group_cols_mapping_labels

  # Define patterns and corresponding groups
  patterns <- list(
    "SARS-CoV-2 variant" = "wildtype|alpha|delta|omicron",
    "Disease severity" = "hospitalized|pams",
    "Prior infection" = "prior_infection",
    "Gender" = "female|male",
    "Age" = "age_cat"
  )

  # Function to find the first matching pattern
  find_group <- function(value) {
    for (group in names(patterns)) {
      if (grepl(patterns[[group]], value)) {
        return(group)
      }
    }
    return("unknown")
  }
  labels <- c(PROLIFERATION_LABEL, PEAK_LABEL, CLEARANCE_LABEL)
  counter <- 1
  for (main_param in c("up", "intercept", "down", "slope_up", "slope_down")) {
    cols_current <- unlist(unname(group_cols_avg[[main_param]]))
    cols_hdi_current <- paste0(cols_current, "_hdi")
    cols_hdi50_current <- paste0(cols_current, "_hdi50")
    cols_reshape <- c(
      "model_no_stan", "min_n_tests", cols_current, cols_hdi_current,
      cols_hdi50_current
    )
    name_mapping_current <- group_cols_avg_labels[[main_param]]

    model_levels <- c(1, 2)
    model_labels <- c("Model 1", "Model 2")
    model_nos <- c(1, 2)
    model_n_levels <- rev(c("Model 1 3", "Model 1 6", "Model 2 3", "Model 2 6"))
    model_n_labels <- rev(c(
      "Model 1, 3 PCR days", "Model 1, 6 PCR days", "Model 2, 3 PCR days",
      "Model 2, 6 PCR days"
    ))

    dt_reshaped <- reshape_dt(dt[, ..cols_reshape], hdi_cols = TRUE)
    dt_reshaped[, model := factor(model_no_stan, levels = model_levels, 
                                  labels = model_labels)]
    dt_reshaped[, model_n := factor(paste(model, min_n_tests, sep = " "),
      levels = model_n_levels,
      labels = model_n_labels
    )]
    # Create the new column
    dt_reshaped[, group := factor(sapply(statistic, find_group),
      levels = c(
        "SARS-CoV-2 variant",
        "Disease severity",
        "Gender",
        "Age",
        "Prior infection"
      )
    )]
    dt_reshaped[, statistic := factor(statistic,
      levels = rev(cols_current),
      labels = rev(unname(name_mapping_current))
    )]
    dt_reshaped[, model_n_combination := factor(paste0(model_no_stan, 
                                                       min_n_tests))]

    dodge_height <- 0.7
    p <- ggplot(dt_reshaped, aes(
      x = value, y = statistic,
      color = model_n,
      group = model_n
    ))
    if (diff) {
      p <- p + geom_vline(xintercept = 0, color = "darkgray", 
                          linetype = "solid", linewidth = 0.4)
    }
    p <- p +
      scale_color_manual(
        values = c(
          "Model 1, 3 PCR days" = "#FF6B6B",
          "Model 1, 6 PCR days" = "#C91818",
          "Model 2, 3 PCR days" = "#4ECDC4",
          "Model 2, 6 PCR days" = "#1E9AA8"
        ),
      ) +
      geom_stripped_rows(color = NA) +
      geom_errorbar(aes(xmin = lower, xmax = upper, group = model_n),
        position = ggstance::position_dodgev(height = dodge_height),
        width = 0,
        linewidth = 0.15
      ) +
      geom_errorbar(aes(xmin = lower50, xmax = upper50, group = model_n),
        position = ggstance::position_dodgev(height = dodge_height),
        width = 0,
        linewidth = 0.45
      ) +
      geom_point(
        position = ggstance::position_dodgev(height = dodge_height),
        alpha = ALPHA_SIM_PLOTS,
        size = DATA_POINT_SIZE_SIM_PLOTS - 1.7
      ) +
      labs(x = labels[counter], y = "") +
      theme_thesis() +
      theme(
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),# Remove major horizontal gridlines
        panel.grid.minor.y = element_blank(),# Remove minor horizontal gridlines
        legend.position = "top",
        plot.margin = margin(1, 4, 1, 1, "pt")# Top, right, bottom, left margins
      ) +
      guides(color = guide_legend(nrow = 2, byrow = FALSE, reverse = TRUE))

    diff_suffix <- ifelse(diff, "_diff", "")
    ggsave(
      filename = file.path(plot_dir, paste0(
        "group_avg_",
        main_param,
        diff_suffix,
        "_iter",
        sampling_iterations,
        ".png"
      )),
      plot = p,
      width = COL_WIDTH * 2,
      height = height_plot * 2.5,
      units = "in", # Units: inches ("in"), centimeters ("cm"), etc.
      dpi = DPI
    )
    counter <- counter + 1
  }
}


#' Create and save plots showing posterior predictions (medians and 94% HPDIs) 
#' of proliferation time, peak and clearance rate for different infection 
#' subgroups (computed from population-level parameter estimates).
#'
#' @param plot_dir A string specifying the directory path where figures will be 
#' saved.
#' @param height_plot A numeric specifying the plot height (in inches).
#' @param diff A logical value indicating whether estimates are given as 
#' differences to the respective reference group:
#' "non-PAMS, non-hospitalized" for disease severity status,
#' "pre-VOC" for SARS-CoV-2 variant,
#' "<30 years" for age group,
#' "no prior infection" for immunization status,
#' "female" for gender.
#' @param ref_variant A string indicating the variant used for prediction.
#' @param ref_hosp_status A string indicating the disease severity status 
#' (PAMS or hospitalized) used for prediction.
#' @param model_nos A vector of model numbers to include (refers to stan file 
#' numbers).
#' @param sampling_iterations A numeric value specifying the number of sampling 
#' iterations run during model fitting. Show results only for the runs with the 
#' corresponding number of sampling iterations.
#'
#' @return Save plots showing posterior predictions (medians and 94% HPDIs) for
#' different infection subgroups (computed from population-level parameter 
#' estimates).
plot_group_post_pred <- function(plot_dir, height_plot, diff = FALSE,
                                 ref_variant = "omicron",
                                 ref_hosp_status = "hospitalized",
                                 model_nos = c(1, 2),
                                 sampling_iterations = 1000) {
  assert_that(ref_variant %in% c("wildtype", "omicron"))
  assert_that(ref_hosp_status %in% c("pams", "hospitalized"))
  dt <- get_post_pred_group_dt(
    ref_variant = ref_variant,
    ref_hosp_status = ref_hosp_status,
    model_nos = model_nos,
    sampling_iterations = sampling_iterations
  )
  # Does not average over infection-level parameters to compute main model 
  # parameters within subgroups but rather computes them from group parameters
  variant_suffix <- ifelse(ref_variant == "omicron", "_omi", "_wt")
  hosp_status_suffix <- ifelse(ref_hosp_status == "hospitalized", 
                               "_hosp", "_pams")
  labels <- c(PROLIFERATION_LABEL, PEAK_LABEL, CLEARANCE_LABEL)
  counter <- 1
  for (main_param in c(
    "time2peak", "intercept", "time_from_peak", "slope_up",
    "slope_down"
  )) {
    col_current <- c(main_param)

    colname_suffix <- ifelse(diff, "_diff", "")
    main_param <- paste0(main_param, colname_suffix)
    col_hdi_lower_current <- paste0(col_current, "_hdi_lower", colname_suffix)
    col_hdi_upper_current <- paste0(col_current, "_hdi_upper", colname_suffix)

    col_hdi50_lower_current <- paste0(col_current, "_hdi50_lower", 
                                      colname_suffix)
    col_hdi50_upper_current <- paste0(col_current, "_hdi50_upper", 
                                      colname_suffix)


    model_levels <- c(1, 2)
    model_labels <- c("Model 1", "Model 2")
    model_nos <- c(1, 2)
    model_n_levels <- rev(c("Model 1 3", "Model 1 6", "Model 2 3", "Model 2 6"))
    model_n_labels <- rev(c(
      "Model 1, 3 PCR days", "Model 1, 6 PCR days", "Model 2, 3 PCR days",
      "Model 2, 6 PCR days"
    ))

    if (diff) {
      group_levels <- c(
        unname(unlist(SUBGRP_CATS_LABEL_MAPPING2["variant"]))[2:4],
        unname(unlist(SUBGRP_CATS_LABEL_MAPPING2["PAMS3"]))[2:3],
        unname(unlist(SUBGRP_CATS_LABEL_MAPPING2["gender"]))[2],
        unname(unlist(SUBGRP_CATS_LABEL_MAPPING2["age_category3"]))[2:3],
        unname(unlist(SUBGRP_CATS_LABEL_MAPPING2["prior_infection"]))[2]
      )
    } else {
      group_levels <- unname(unlist(c(
        SUBGRP_CATS_LABEL_MAPPING2["variant"][1:4],
        SUBGRP_CATS_LABEL_MAPPING2["PAMS3"],
        SUBGRP_CATS_LABEL_MAPPING2["gender"],
        SUBGRP_CATS_LABEL_MAPPING2["age_category3"],
        SUBGRP_CATS_LABEL_MAPPING2["prior_infection"]
      )))
    }
    dt[, model := factor(model_no_stan, levels = model_levels, 
                         labels = model_labels)]
    dt[, model_n := factor(paste(model, min_n_tests, sep = " "),
      levels = model_n_levels,
      labels = model_n_labels
    )]
    dt[, model_n_combination := factor(paste0(model_no_stan, min_n_tests))]
    dt <- dt[group != "Unknown"]
    dt[, group := factor(group, levels = rev(group_levels))]
    # Drop rows where 'value' is NA
    dt <- dt[!is.na(group), ]

    dodge_height <- 0.7
    p <- ggplot(dt, aes(x = !!sym(main_param), y = group, color = model_n, 
                        group = model_n))
    if (diff) {
      p <- p + geom_vline(xintercept = 0, color = "darkgray", 
                          linetype = "solid", linewidth = 0.4)
    }
    p <- p +
      scale_color_manual(
        values = c(
          "Model 1, 3 PCR days" = "#FF6B6B",
          "Model 1, 6 PCR days" = "#C91818",
          "Model 2, 3 PCR days" = "#4ECDC4",
          "Model 2, 6 PCR days" = "#1E9AA8"
        ),
      ) +
      geom_stripped_rows(color = NA) +
      geom_errorbar(
        aes(
          xmin = !!sym(col_hdi_lower_current),
          xmax = !!sym(col_hdi_upper_current), group = model_n
        ),
        position = ggstance::position_dodgev(height = dodge_height),
        width = 0,
        linewidth = 0.15
      ) +
      geom_errorbar(
        aes(
          xmin = !!sym(col_hdi50_lower_current),
          xmax = !!sym(col_hdi50_upper_current), group = model_n
        ),
        position = ggstance::position_dodgev(height = dodge_height),
        width = 0,
        linewidth = 0.45
      ) +
      geom_point(
        position = ggstance::position_dodgev(height = dodge_height),
        alpha = ALPHA_SIM_PLOTS,
        size = DATA_POINT_SIZE_SIM_PLOTS - 1.7
      ) +
      labs(x = labels[counter], y = "") +
      theme_thesis() +
      theme(
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),# Remove major horizontal gridlines
        panel.grid.minor.y = element_blank(),# Remove minor horizontal gridlines
        legend.position = "top",
        plot.margin = margin(1, 4, 1, 1, "pt")# Top, right, bottom, left margins
      ) +
      guides(color = guide_legend(nrow = 2, byrow = FALSE, reverse = TRUE))

    if (main_param == "intercept") {
      p <- p + scale_x_continuous(n.breaks = 4) # Show approximately 4 labels
    }

    ggsave(
      filename = file.path(plot_dir, paste0(
        "group_post_pred_",
        main_param,
        variant_suffix,
        hosp_status_suffix,
        "_iter",
        sampling_iterations,
        ".png"
      )),
      plot = p,
      width = COL_WIDTH * 2,
      height = height_plot * 2,
      units = "in", # Units: inches ("in"), centimeters ("cm"), etc.
      dpi = DPI
    )
    counter <- counter + 1
  }
}


#' Create and save plot showing posterior estimates (medians and 94% HPDIs), 
#' e.g. for the proliferation time, for a specific model.
#'
#' @param dt A data.table containing posterior parameter estimates (medians and 
#' 94% HPDIs) for a specific model and different runs, e.g. with different 
#' thresholds for the minimum number of PCR days.
#' @param x Column name for x-axis values (e.g. proliferation time).
#' @param y Column name for y-axis values (e.g. PCR days threshold used - if 
#' several runs with different thresholds were performed).
#' @param color Column name for grouping/coloring points.
#' @param alpha Transparency of points (0 to 1).
#' @param height_dodge Dodge height for points and error bars.
#' @param size_point Size of points.
#' @param xlabel Label for x-axis.
#' @param ylabel Label for y-axis.
#' @param color_label Label for color legend.
#' @param annotation_label Optional annotation tag for the plot.
#' @param left_margin Left margin for plot annotation.
#' @param hjust Horizontal adjustment for annotation.
#' @param as_percent Logical value indicating whether to convert x-axis values 
#' and intervals to percentages.
#'
#' @return Saves plot showing posterior estimates (medians and 94% HPDIs), e.g.
#' for the proliferation time, for a specific model.
plot_main_model_params <- function(dt, x, y, color, alpha,
                                   height_dodge, size_point,
                                   xlabel, ylabel, color_label,
                                   annotation_label = NULL,
                                   left_margin = ANNOTATION_MARGIN_LEFT_SIM_PLOT,
                                   hjust = ANNOTATION_HJUST_SIM_PLOT,
                                   as_percent = FALSE) {
  col_lower <- paste0(x, "_lower")
  col_upper <- paste0(x, "_upper")
  col_lower50 <- paste0(x, "_lower50")
  col_upper50 <- paste0(x, "_upper50")
  dt_new <- copy(dt)
  x_hdi <- paste0(x, "_hdi")
  x_hdi50 <- paste0(x, "_hdi50")

  if (as_percent) {
    dt_new[, eval(x) := get(x) * 100]
    dt_new[, c(col_lower, col_upper) := transpose(get(x_hdi))]
    dt_new[, c(col_lower50, col_upper50) := transpose(get(x_hdi50))]
    dt_new[, eval(col_lower) := get(col_lower) * 100]
    dt_new[, eval(col_upper) := get(col_upper) * 100]
    dt_new[, eval(col_lower50) := get(col_lower50) * 100]
    dt_new[, eval(col_upper50) := get(col_upper50) * 100]
  } else {
    dt_new[, c(col_lower, col_upper) := transpose(get(x_hdi))]
    dt_new[, c(col_lower50, col_upper50) := transpose(get(x_hdi50))]
  }


  p <- ggplot(dt_new, aes(
    x = !!sym(x), y = !!sym(y), color = !!sym(color),
    group = !!sym(color)
  )) +
    scale_color_viridis_d() +
    geom_errorbar(aes(xmin = !!sym(col_lower), xmax = !!sym(col_upper)),
      position = ggstance::position_dodgev(height = height_dodge),
      width = 0,
      linewidth = 0.15
    ) +
    geom_errorbar(aes(xmin = !!sym(col_lower50), xmax = !!sym(col_upper50)),
      position = ggstance::position_dodgev(height = dodge_height),
      width = 0,
      linewidth = 0.45
    ) +
    geom_point(
      position = ggstance::position_dodgev(height = height_dodge), 
      alpha = alpha,
      size = size_point
    ) +
    geom_stripped_rows(color = NA) +
    labs(
      title = "",
      x = xlabel,
      y = ylabel,
      color = paste(color_label, sep = "\n")
    ) +
    theme_thesis() +
    theme(
      panel.grid.major.y = element_blank(), # Remove major horizontal gridlines
      panel.grid.minor.y = element_blank()
    ) # Remove minor horizontal gridlines
  if (!is.null(annotation_label)) {
    p <- plot_annotation(p, letter = annotation_label, hjust = hjust, 
                         left_margin = left_margin)
  }
  return(p)
}
