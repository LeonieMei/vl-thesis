library(brms)
library(optparse)
library(here)

TOP <- here()
source(file.path(TOP, "sars2vl", "data_funcs.R"))


option_list <- list(
  make_option(c("--group"),
    default = "variant",
    help = "You can choose from variant, pams, age and pams_variant."
  ),
  make_option(c("--likelihood"),
    default = "student-t",
    help = "You can choose from student-t and skewnormal."
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
assert_that(opt$group %in% c("variant", "pams", "age", "pams_variant"))
assert_that(opt$likelihood %in% c("student-t", "skewnormal"))


priors_studentt <- c(
  set_prior("normal(0, 2)", class = "b"), # Prior for fixed effects (slopes)
  set_prior("normal(7, 2)", class = "Intercept"), # Prior for the intercept
  set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma"),
  set_prior("normal(0, 0.5)", class = "b", dpar = "sigma")
)


main <- function(opt) {
  dt <- load_first_pos_vl_data()
  # We need to make sure that age category is not ordered, else brms does
  # some polynomial fit
  dt$age_category3 <- factor(dt$age_category3, ordered = FALSE)
  dt$age_category3 <- relevel(dt$age_category3, ref = "<30")


  if (opt$likelihood == "student-t") {
    priors <- priors_studentt
  } else {
    priors <- NULL
  }

  # Specify and fit models
  if (opt$group == "variant") {
    if (opt$likelihood == "student-t") {
      formula_variant <- bf(log10_load ~ variant, sigma ~ variant)
    } else {
      formula_variant <- bf(log10_load ~ variant, alpha ~ variant)
    }
    variant_fit <- fit_first_vl_model(formula_variant,
      data = dt[!is.null(variant)],
      fn = "variant_model",
      priors = priors,
      family = opt$likelihood
    )
  } else if (opt$group == "pams") {
    if (opt$likelihood == "student-t") {
      formula_PAMS1 <- bf(log10_load ~ PAMS1, sigma ~ PAMS1)
    } else {
      formula_PAMS1 <- bf(log10_load ~ PAMS1, alpha ~ PAMS1)
    }
    pams_fit <- fit_first_vl_model(formula_PAMS1,
      data = dt[!is.null(PAMS1)],
      fn = "pams_model",
      priors = priors,
      family = opt$likelihood
    )
  } else if (opt$group == "pams_variant") {
    if (opt$likelihood == "student-t") {
      formula_PAMS_variant <- bf(log10_load ~ variant + PAMS1:variant, 
                                 sigma ~ variant + PAMS1:variant)
    } else {
      formula_PAMS_variant <- bf(log10_load ~ variant + PAMS1:variant, 
                                 alpha ~ variant + PAMS1:variant)
    }
    pams_variant_fit <- fit_first_vl_model(formula_PAMS_variant,
      data = dt[!is.null(PAMS1) & !is.null(variant)],
      fn = "pams_variant_model",
      priors = priors,
      family = opt$likelihood
    )
  } else {
    if (opt$likelihood == "student-t") {
      formula_age <- bf(log10_load ~ age_category3, sigma ~ age_category3)
    } else {
      formula_age <- bf(log10_load ~ age_category3, alpha ~ age_category3)
    }
    age_fit <- fit_first_vl_model(formula_age,
      data = dt[!is.null(age_category3)],
      fn = "age_model",
      priors = priors,
      family = opt$likelihood
    )
  }
}

main(opt)
