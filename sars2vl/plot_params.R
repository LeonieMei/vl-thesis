library(here)

TOP <- here()
source(file.path(TOP, "sars2vl", "data_params.R"))

DIR_FIGURES_VL_SIM <- file.path(DIR_FIGURES, "vl_trajectories", "simulations")
DIR_FIGURES_VL <- file.path(DIR_FIGURES, "vl_trajectories", "real")

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

DINA4_WIDTH <- 5.91
# Parameters for plots
COL_WIDTH <- DINA4_WIDTH / 2
DPI <- 300

ADD_SIZE <- 0
GRAY_SHADE <- "gray80"
TITLE_SIZE <- 14 + ADD_SIZE
TEXT_SIZE <- 12 + ADD_SIZE
LEGEND_TITLE_SIZE <- 12 + ADD_SIZE
LEGEND_TEXT_SIZE <- 10 + ADD_SIZE
LEGEND_KEY_SIZE <- 20
PLOT_TITLE_SIZE <- 14 + ADD_SIZE
ANNOTATION_LETTER_SIZE <- 16

# Plot params
ALPHA_SIM_PLOTS <- 0.85
JITTER_WIDTH_SIM_PLOTS <- 0.1
JITTER_HEIGTH_SIM_PLOTS <- 0
DATA_POINT_SIZE_SIM_PLOTS <- 3
HEIGHT_SIM_PLOT <- 2.5
HEIGHT_RHAT_PLOT <- HEIGHT_SIM_PLOT
HEIGHT_RHAT_PLOTS <- HEIGHT_RHAT_PLOT * 1.8
HEIGHT_PRIOR_PRED_PLOT <- HEIGHT_SIM_PLOT * 1.5
ALPHA_SIM_PLOT <- 0.85
WIDTH_JITTER_SIM <- 0.5
HEIGHT_JITTER_SIM <- 0
GEOM_POINT_SIZE_SIM <- 3
ANNOTATION_MARGIN_LEFT_SIM_PLOT <- 20
ANNOTATION_HJUST_SIM_PLOT <- 2


# The dummy "variant" is there so we get a different color for Omicron as otherwise the
# color is too close to the previous one (which is the "Unknown" variant)
VARIANT2_COLORS <- setNames(c("#F8766D", "#7CAE00", "#00BFC4", gg_color_hue(length(VARIANT2_LABELS) + 3)[4:(length(VARIANT2_LABELS) - 1)], "#D3D3D3"), VARIANT2_LABELS)
VARIANT2_COLORS <- VARIANT2_COLORS[VARIANT2_LABELS]
VARIANT2_COLORS2 <- c(VARIANT2_COLORS[1:7], VARIANT2_COLORS[length(VARIANT2_COLORS)])
VARIANT_COLORS <- setNames(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#D3D3D3"), VARIANT_LABELS)
VARIANT_COLORS <- VARIANT_COLORS[VARIANT_LABELS]
AGE_CODE_COLORS <- setNames(gg_color_hue(3), AGE_LABELS3)
AGE_LEGEND_LABELS <- setNames(AGE_LABELS3, AGE_LABELS3)
GENDER_COLORS <- setNames(gg_color_hue(2), c("F", "M"))

VARIANT_LABEL <- "SARS-CoV-2 variant"
VL_LABEL <- expression(Log[10] ~ viral ~ load)
VL_LABEL_LOWER <- expression(log[10] ~ viral ~ load)
VL_LABEL_LATEX <- "$\\textbf{Log}_{\\textbf{10}} \\textbf{viral load}$"
VL_LABEL_LATEX_LOWER <- "$\\textbf{log}_{\\textbf{10}} \\textbf{viral load}$"
COUNT_VL_LABEL_LATEX <- "\\textbf{Counts and} $\\textbf{log}_{\\textbf{10}} \\textbf{viral load}$"
PCR_LABEL <- "PCR target"
PAMS_LABEL <- "PAMS status"
GENDER_LABEL <- "Gender"
AGE_LABEL <- "Age"
DAYS_LABEL <- "Day since peak viral load"

VARIABLE_LABEL_DICT <- list(
  "variant" = VARIANT_LABEL, 
  "PCR" = PCR_LABEL, 
  "PAMS1" = PAMS_LABEL, 
  "age" = AGE_LABEL, 
  "gender" = GENDER_LABEL, 
  "age_group3" = AGE_LABEL
)

VARIABLE_LABEL_DICT_LATEX <- list(
  "variant" = "\\shortstack{SARS-CoV-2\\\\variant}", 
  "PCR" = PCR_LABEL, 
  "PAMS1" = PAMS_LABEL, 
  "age" = AGE_LABEL, 
  "gender" = GENDER_LABEL, 
  "age_group3" = AGE_LABEL
)

VARIABLE_LABEL_DICT_LATEX2 <- list(
  "variant" = "\\shortstack{SARS-CoV-2\\\\variant}", 
  "PAMS3" = "Disease severity", 
  "age" = AGE_LABEL,
  "gender" = GENDER_LABEL, 
  "age_category3" = AGE_LABEL, 
  "age_group3" = AGE_LABEL, 
  "prior_infection" = "Prior infection",
  "material" = "Swabbing site"
)

MINIMUM_PCRS_LABEL <- "Minimum PCR days"
MINIMUM_PCRS_LABEL2 <- "Minimum PCR days"
MINIMUM_POS_PCRS_LABEL <- "Minimum positive PCR days"
MINIMUM_POS_PCRS_LABEL2 <- "Minimum\npositive\nPCR days"


PROLIFERATION_LABEL <- "Proliferation time (days)"
CLEARANCE_LABEL <- "Clearance time (days)"
UP_SLOPE_LABEL <- "Up-slope"
DOWN_SLOPE_LABEL <- "Down-slope"
UP_SLOPE_LABEL2 <- "Proliferation rate"
DOWN_SLOPE_LABEL2 <- "Clearance rate"
PEAK_LABEL <- expression(Peak ~ log[10] ~ viral ~ load)
PROLIFERATION_LABEL2 <- bquote(atop(NA, atop(textstyle("Proliferation time"))))
DOWN_SLOPE_LABEL2_2 <- bquote(atop(NA, atop(textstyle("Clearance rate"))))
CLEARANCE_LABEL2 <- bquote(atop(NA, atop(textstyle("Clearance time"), textstyle("(days)"))))
PEAK_LABEL2 <- bquote("Peak" ~ log[10] ~ "viral load")
