library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(lubridate)
library(here)

TOP <- here()
source(file.path(TOP, "plot_params.R"))
source(file.path(TOP, "data_params.R"))
source(file.path(TOP, "data_utils.R"))
source(file.path(TOP, "plot_utils.R"))

DODGE_WIDTH <- .7
MY_LEGEND_POSITION <- c(0.3, 0.8)

DIR_FIGURES_SWARMPLOT <- file.path(DIR_FIGURES, "swarmplots")
if (!dir.exists(DIR_FIGURES_SWARMPLOT)) {
  dir.create(DIR_FIGURES_SWARMPLOT, recursive = TRUE)
}


prepare_data_table <- function(dt) {
  # Show first viral loads over time.
  dt <- dt %>%
    .[order(ymd(date))] %>%
    .[, year_month := factor(format(date, "%Y-%m"),
      levels = sort(unique(format(dt$date, "%Y-%m"))),
      ordered = TRUE
    )] %>%
    .[, year_week := factor(format(date, "%Y-%U"),
      levels = sort(unique(format(dt$date, "%Y-%U"))),
      ordered = TRUE
    )] %>%
    .[, first_day_week_date := min(date), by = year_week] %>%
    .[, n_tests := .N, by = year_week] %>%
    .[, mean_vl := mean(log10_load, na.rm = TRUE), by = year_week] %>%
    .[, median_vl := median(log10_load, na.rm = TRUE), by = year_week] %>%
    .[order(first_day_week_date)]

  # Make a new column with every two or four months combined.
  year_month_levels <- levels(dt$year_month)
  year_month1 <- year_month_levels[seq(1, length(year_month_levels), by = 2)]
  year_month2 <- year_month_levels[seq(2, length(year_month_levels), by = 2)]
  year_month3 <- year_month_levels[seq(1, length(year_month_levels), by = 4)]
  year_month4 <- year_month_levels[seq(2, length(year_month_levels), by = 4)]
  year_month5 <- year_month_levels[seq(3, length(year_month_levels), by = 4)]
  year_month6 <- year_month_levels[seq(4, length(year_month_levels), by = 4)]

  year_month_dict <- c(
    set_names(year_month1[-length(year_month1)], year_month2),
    set_names(
      year_month1[length(year_month1)],
      year_month1[length(year_month1)]
    )
  )
  year_month_labels <- c(
    paste(year_month1[-length(year_month1)],
      year_month2,
      sep = "-"
    ),
    year_month1[length(year_month1)]
  )
  year_month_dict2 <- c(
    set_names(year_month3, year_month3),
    set_names(year_month3, year_month4),
    set_names(year_month3, year_month5),
    set_names(year_month3[-length(year_month3)], year_month6),
    set_names(
      year_month3[length(year_month3)],
      year_month3[length(year_month3)]
    )
  )
  year_month_labels2 <- c(
    paste(year_month3[-length(year_month3)],
      year_month6,
      sep = "-"
    ),
    paste0(year_month3[length(year_month3)], "-", "2025-08")
  )
  year_month_labels3 <- c("", "", "", "", "")
  year_month_labels3 <- c(
    "", "", "", "", "", "2021", "", "", "", "", "", "2022",
    "", "", "", "", "", "2023", "", "", "", "", "",
    "2024", "", "", "", "", "", "2025", "", "", "", ""
  )
  year_month_labels4 <- c(
    "", "", "2021", "", "", "2022", "", "", "2023", "", "",
    "2024", "", "", "2025", "", ""
  )

  dt <-
    dt %>%
    .[, year_month2 := as.character(year_month)] %>%
    .[, year_month2 := factor(
      ifelse(year_month2 %in% names(year_month_dict),
        year_month_dict[year_month2], year_month2
      ),
      levels = year_month_dict, ordered = TRUE,
      labels = year_month_labels
    )] %>%
    .[, year_month3 := as.character(year_month)] %>%
    .[, year_month3 := factor(year_month_dict2[year_month3],
      levels = unique(year_month_dict2), ordered = TRUE,
      labels = year_month_labels2
    )]

  return(list(data = dt, label1 = year_month_labels3, label2 = year_month_labels4))
}


make_swarmplots_group <- function(dt, by = "PAMS1", xlabel = "PAMS status") {
  assert_that(by %in% c("PAMS1", "age", "PCR", "variant", "variant2", "gender"))
  if (by == "age") {
    dt[, age_group3 := cut(age,
      breaks = AGE_BREAKS3, ordered = TRUE,
      labels = AGE_LABELS3
    )]
    by <- "age_group3"
  } else if (by == "PAMS1") {
    dt[, PAMS1_new := factor(ifelse(PAMS1 == 0, "non-PAMS", "PAMS"))]
    by <- "PAMS1_new"
  } else if (by == "variant") {
    dt[, variant2_new := factor(variant2,
      levels = VARIANT_LEVELS,
      labels = VARIANT_LEGEND_LABELS
    )]
    by <- "variant2_new"
  } else if (by == "variant2") {
    dt[, variant2_new := factor(variant2,
      levels = VARIANT2_LEVELS,
      labels = VARIANT2_LEGEND_LABELS
    )]
    by <- "variant2_new"
  }

  vl_dist_by_grp <- ggplot(dt, aes(x = get(by), y = log10_load)) +
    geom_quasirandom(
      method = "tukeyDense", size = 0.3, alpha = 0.5, color = "blue",
      stroke = 0
    ) +
    geom_boxplot(
      fill = NA, # No fill
      color = "black", # Outline color
      alpha = 1, # Fully opaque
      width = 0.2 # Adjust width as needed
    ) +
    xlab(xlabel) +
    ylab(VL_LABEL) +
    theme(
      axis.text = element_text(size = TEXT_SIZE),
      axis.title = element_text(size = TEXT_SIZE),
      legend.text = element_text(size = LEGEND_TEXT_SIZE)
    ) +
    theme_minimal()

  if (by == "variant2_new") {
    vl_dist_by_grp <- vl_dist_by_grp + theme(axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ))
  }
  vl_dist_by_grp <- plot_annotation(plot = vl_dist_by_grp, letter = "B")

  return(vl_dist_by_grp)
}


make_swarmplots_by_time <- function(dt, month_var = "year_month2", color_by = "PAMS1",
                                 label = NULL) {
  # assertthat::assert_that(color_by %in% c("PAMS1", "age", "PCR", "variant2"))
  assert_that(color_by %in% c(
    "PAMS1", "age", "PCR", "variant", "variant2",
    "gender"
  ))
  if (month_var %in% c("year_month2", "year_month3")) {
    assert_that(!is.null(label))
  }
  dt_curr <- copy(dt)

  if (color_by == "PAMS1") {
    color_by <- "PAMS1_1"
    dt_curr %>% .[, PAMS1_1 := factor(ifelse(PAMS1 == 0, "non-PAMS", "PAMS"))]
  } else if (color_by == "PCR") {
    color_by <- "PCR_1"
    dt_curr %>% .[, PCR_1 := factor(PCR)]
  }

  if (month_var == "year_month") {
    dateLimit <- "2020-03"
  } else if (month_var == "year_month2") {
    dateLimit <- "2020-02-2020-03"
  } else if (month_var == "year_month3") {
    dateLimit <- "2020-02-2020-05"
  }

  vl_dist_by_month <-
    dt_curr[get(month_var) >= dateLimit] %>%
    ggplot(aes(x = get(month_var), y = log10_load, color = get(color_by))) +
    geom_quasirandom(method = "tukeyDense", size = .3, alpha = 0.7, stroke = 0) +
    geom_boxplot(
      fill = NA, # No fill
      color = "black", # Outline color
      alpha = 1, # Fully opaque
      width = 0.2 # Adjust width as needed
    ) +
    xlab("Time") +
    ylab(VL_LABEL) +
    theme(
      axis.text = element_text(size = TEXT_SIZE),
      axis.title = element_text(size = TEXT_SIZE),
      legend.text = element_text(size = LEGEND_TEXT_SIZE)
    ) +
    theme_minimal()

  if (month_var == "year_month") {
    vl_dist_by_month <- vl_dist_by_month + scale_x_discrete(breaks = levels(dt$year_month)[c(FALSE, TRUE, FALSE)])
  } else if (month_var == "year_month2") {
    vl_dist_by_month <- vl_dist_by_month + scale_x_discrete(labels = label)
  } else if (month_var == "year_month3") {
    vl_dist_by_month <- vl_dist_by_month + scale_x_discrete(labels = label)
  }

  if (color_by == "PAMS1_1") {
    vl_dist_by_month <- vl_dist_by_month + scale_color_manual(values = c("non-PAMS" = "red", "PAMS" = "blue")) +
      # Increase the size of markers in the legend
      labs(color = "PAMS") + guides(color = guide_legend(override.aes = list(size = 3))) +
      theme(legend.title = element_blank()) # getting rid of the legend title here.
  } else if (color_by == "age") {
    vl_dist_by_month <- vl_dist_by_month + labs(color = "Age (years)") + scale_colour_gradientn(colours = rainbow(n = 10)) +
      theme(legend.title = element_text(size = LEGEND_TITLE_SIZE))
  } else if (color_by == "PCR_1") {
    vl_dist_by_month <- vl_dist_by_month + scale_color_manual(values = gg_color_hue(uniqueN(dt_curr$PCR_1))) +
      labs(color = "PCR") +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      theme(legend.title = element_text(size = LEGEND_TITLE_SIZE))
  } else if (color_by == "variant") {
    vl_dist_by_month <- vl_dist_by_month + scale_color_manual(values = VARIANT_COLORS) +
      labs(color = "SARS-CoV-2 variant") +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      theme(legend.title = element_text(size = LEGEND_TITLE_SIZE))
  } else if (color_by == "variant2") {
    vl_dist_by_month <- vl_dist_by_month + scale_color_manual(values = VARIANT2_COLORS) +
      labs(color = "SARS-CoV-2 variant") +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      theme(legend.title = element_text(size = LEGEND_TITLE_SIZE))
  } else if (color_by == "gender") {
    vl_dist_by_month <- vl_dist_by_month + scale_color_manual(values = GENDER_COLORS) +
      labs(color = "SARS-CoV-2 variant") +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      theme(legend.title = element_text(size = LEGEND_TITLE_SIZE))
  }
  return(plot_annotation(vl_dist_by_month, letter = "A"))
}


main <- function() {
  dt <- load_first_pos_vl_data()
  data_n_labels <- prepare_data_table(dt)
  dt <- data_n_labels$data
  month_varVal <- 3
  month_var <- paste0("year_month", month_varVal)
  if (month_varVal == 2) {
    label <- data_n_labels$label1
  } else if (month_varVal == 3) {
    label <- data_n_labels$label2
  }

  vl_dist_by_pams <- make_swarmplots_group(dt, by = "PAMS1", xlabel = PAMS_LABEL)
  vl_dist_by_age <- make_swarmplots_group(dt, by = "age", xlabel = AGE_LABEL)
  vl_dist_by_pcr <- make_swarmplots_group(dt, by = "PCR", xlabel = PCR_LABEL)
  vl_dist_by_variant <- make_swarmplots_group(dt, by = "variant", xlabel = VARIANT_LABEL)
  vl_dist_by_variant2 <- make_swarmplots_group(dt,
    by = "variant2",
    xlabel = VARIANT_LABEL
  )
  vl_dist_by_gender <- make_swarmplots_group(dt, by = "gender", xlabel = GENDER_LABEL)

  ggsave(vl_dist_by_pams,
    file = file.path(DIR_FIGURES_SWARMPLOT, "VL_dist_by_PAMS.png"),
    height = 3, width = COL_WIDTH * 2, units = "in", dpi = DPI
  )
  ggsave(vl_dist_by_age,
    file = file.path(DIR_FIGURES_SWARMPLOT, "VL_dist_by_age.png"),
    height = 3, width = COL_WIDTH * 2, units = "in", dpi = DPI
  )
  ggsave(vl_dist_by_pcr,
    file = file.path(DIR_FIGURES_SWARMPLOT, "VL_dist_by_pcr.png"),
    height = 3, width = COL_WIDTH * 2, units = "in", dpi = DPI
  )
  ggsave(vl_dist_by_variant,
    file = file.path(DIR_FIGURES_SWARMPLOT, "VL_dist_by_variant.png"),
    height = 3, width = COL_WIDTH * 2, units = "in", dpi = DPI
  )
  ggsave(vl_dist_by_variant2,
    file = file.path(DIR_FIGURES_SWARMPLOT, "VL_dist_by_variant2.png"),
    height = 3, width = COL_WIDTH * 2, units = "in", dpi = DPI
  )
  ggsave(vl_dist_by_gender,
    file = file.path(DIR_FIGURES_SWARMPLOT, "VL_dist_by_gender.png"),
    height = 3, width = COL_WIDTH * 2, units = "in", dpi = DPI
  )


  vl_dist_by_month <- make_swarmplots_by_time(dt,
    color_by = "PAMS1",
    month_var = month_var, label = label
  )
  vl_dist_by_month_age <- make_swarmplots_by_time(dt,
    color_by = "age",
    month_var = month_var, label = label
  )
  vl_dist_by_month_pcrs <- make_swarmplots_by_time(dt,
    color_by = "PCR",
    month_var = month_var, label = label
  )
  vl_dist_by_month_variant <- make_swarmplots_by_time(dt,
    color_by = "variant",
    month_var = month_var, label = label
  ) +
    scale_color_manual(values = VARIANT_COLORS, labels = VARIANT_LEGEND_LABELS)
  vl_dist_by_month_variant2 <- make_swarmplots_by_time(dt,
    color_by = "variant2",
    month_var = month_var, label = label
  ) +
    scale_color_manual(
      values = VARIANT2_COLORS,
      labels = VARIANT2_LEGEND_LABELS
    )

  ggsave(vl_dist_by_month,
    file = file.path(DIR_FIGURES_SWARMPLOT, paste0("VL_dist_by_month", month_varVal, ".png")), 
    height = 3, width = COL_WIDTH * 2,
    units = "in", dpi = DPI
  )
  ggsave(vl_dist_by_month_age,
    file = file.path(DIR_FIGURES_SWARMPLOT, paste0("VL_dist_by_month", month_varVal, "_age.png")), 
    height = 3, width = COL_WIDTH * 2,
    units = "in", dpi = DPI
  )
  ggsave(vl_dist_by_month_pcrs,
    file = file.path(DIR_FIGURES_SWARMPLOT, paste0("VL_dist_by_month", month_varVal, "_pcr.png")), 
    height = 3, width = COL_WIDTH * 2,
    units = "in", dpi = DPI
  )
  ggsave(vl_dist_by_month_variant,
    file = file.path(DIR_FIGURES_SWARMPLOT, paste0("VL_dist_by_month", month_varVal, "_variant.png")),
    height = 3, width = COL_WIDTH * 2, units = "in", dpi = DPI
  )
  ggsave(vl_dist_by_month_variant2,
    file = file.path(DIR_FIGURES_SWARMPLOT, paste0("VL_dist_by_month", month_varVal, "_variant2.png")),
    height = 3.55, width = COL_WIDTH * 2, units = "in", dpi = DPI
  )
}

main()
