library(dplyr)
library(data.table)
library(magrittr)
library(dplyr)
library(stringr)
library(assertthat)
library(truncnorm)
library(brms)
library(here)

TOP <- here()
source(file.path(TOP, "sars2vl", "post_analysis.R"))
source(file.path(TOP, "sars2vl", "data_params.R"))


#### Load & prepare data for analysis ####


#' Filter valid factor levels.
#'
#' @param dt A data.table with input data.
#' @param col_name Character. Column name to check for valid levels.
#' @param levels Character vector. Desired levels to filter.
#'
#' @return Character vector. Subset of `levels` present in `dt[, col_name]`.
filter_levels <- function(dt, col_name, levels) {
  categories <- unique(dt[, get(col_name)])
  return(levels[levels %in% categories])
}


#' Identify centre with longest duration between consecutive tests from the same
#' centre.
#'
#' @param centres A vector with test centres.
#' @param day A vector with day of test (first test is day 0).
#'
#' @return The name of the centre with the longest duration between consecutive
#' tests from the same centre.
get_ld_centre <- function(centres, day, max_streaks = 10) {
  streaks <- data.frame(centre = rep("", max_streaks),
                        duration = rep(NA, max_streaks))
  k <- 0
  duration <- 0
  for (d in 2:length(centres)) {
    if (centres[d] == centres[d - 1] && d < length(centres)) {
      duration <- duration + day[d] - day[d - 1]
    } else if (centres[d] != centres[d - 1]) {
      k <- k + 1
      streaks$centre[k] <- centres[d - 1]
      streaks$duration[k] <- duration
      duration <- 0
    } else if (centres[d] == centres[d - 1] && d == length(centres)) {
      duration <- duration + day[d] - day[d - 1]
      k <- k + 1
      streaks$centre[k] <- centres[d - 1]
      streaks$duration[k] <- duration
    }
  }
  streaks <- streaks[complete.cases(streaks), ]
  streaks <- streaks[order(-streaks$duration), ]
  # No centre has a streak longer than a day (I, Leonie, think):
  if (max(streaks$duration) == 0) {
    tbl <- sort(table(streaks$centre), decreasing = TRUE)
    # Or: all tests were taken on the same day, so there is an entry only
    # for tbl[1]
    if (is.na(tbl[2])) {
      return(names(tbl)[1])
    } else {
      return(
        ifelse(tbl[1] != tbl[2],
               names(tbl)[1], "X"
               )
        )
    }
    # There is exactly one centre with the maximum streak duration:
  } else if (nrow(streaks) == 1 || streaks$duration[1] != streaks$duration[2]) {
    return(streaks$centre[streaks$duration == max(streaks$duration)])
    # There is more than one centre with the maximum streak duration:
  } else {
    tmp <- paste0(streaks[streaks$duration == max(streaks$duration), "centre"],
                  collapse = "_"
    )
    if (grepl("ICU", tmp)) {
      return("ICU")
    } else if (grepl("IDW", tmp)) {
      return("IDW")
    } else if ((!grepl("IDW", tmp)) && (!grepl("ICU", tmp)) && grepl("WD", tmp)) {
      return("WD")
    } else if (grepl("C19", tmp) && grepl("ED", tmp)) {
      return("ED")
    } else {
      # Pick one of the centres with the longest duration
      return(streaks$centre[1])
    }
  }
}


#' Filter out unusual infection trajectories.
#'
#' Removes infections with atypical viral load patterns, such as:
#' - Late-stage high viral loads (log10_load ≥ 8 after 20+ days (post first 
#' positive PCR test)).
#' - Infections with durations exceeding 30 days.
#' - Predefined problematic infections (IDs_TO_EXCLUDE).
#'
#' @param data A data.table. Input data with columns: ID, infection_duration, 
#' day, last_pos_day,
#' log10_load.
#'
#' @return data.table. Filtered data with unusual infections removed.
filter_unusual_infections <- function(data) {
  # Filter out infections that have some reactivation pattern, i.e. they have
  # a high viral load followed by a very low viral load, or infections where the 
  # last positive PCR has a very high viral load if it has a sufficient time 
  # distance to the first positive PCR
  IDs_late_high_positives <- data[(infection_duration >= 20) & (last_pos_day == day) & (log10_load >= 8)]$ID
  # Infections with at least two high peaks that are at least 10 days apart with 
  # at least one PCR in between that is much lower.
  # The following only filters out about 92 of 18296 infections so we are 
  # skipping it.
  # IDsMultiplePeaks = find_multiple_spike_infection(data)
  data <- data %>%
     # This is actually already handled when creating the json file in python.
    .[infection_duration <= 30] %>%
    .[!ID %in% IDs_late_high_positives] %>%
    .[!ID %in% IDs_TO_EXCLUDE]
  
  return(data)
}


#' Convert CT values to log10 viral load (NBA data)
#'
#' @param ct Numeric vector. Cycle threshold (CT) values.
#' @return Numeric vector. Log10-transformed viral load estimates.
ct_to_viralload_T1 <- function(ct) {
  mult <- 6.109410646210822e15
  expt <- -0.753893467915248
  return(log10(mult * exp(expt * ct)))
}


#' Load and preprocess data of first positive PCR tests.
#'
#' @return A data.table.
get_first_pos_pcr_data <- function(vl_first_pos_file = VL_FIRST_POS_FILE_TSV) {
  AGE_breaks <- c(seq(0, 25, 5), seq(35, 65, 10), 101)
  AGE_labels <- c(paste0(c(seq(0, 20, 5), seq(25, 55, 10)), "-",
                         c(seq(5, 25, 5), seq(35, 65, 10))), ">65")
  
  fulldt <- read.table(file = vl_first_pos_file, sep = "\t", header = TRUE)
  setDT(fulldt)
  # Turn all None values to NAs.
  cols <- c("infection_hash", "variant", "n_prior_infections",
            "days_to_prior_infection", "is_first_pos_pcr")
  fulldt <- fulldt %>% mutate(across(cols, ~ replace(.x, .x == "None", NA)))
  
  # Remove PCRs with viral load values <= VL_LOWER_LIMIT (except those with
  # viral load levels of -1 which means negative PCR) as it is not clear whether
  # those are false positives.
  fulldt <- fulldt[(log10_load > VL_LOWER_LIMIT) | (log10_load == -1)]
  assert_that(all((fulldt$log10_load > VL_LOWER_LIMIT) | (fulldt$log10_load == -1), 
                   na.rm = FALSE))
  
  fulldt <-
    copy(fulldt) %>%
    # Remove infection_hash from negative PCRs that don't have a later positive 
    # PCR in the data anymore.
    .[, infection_hash := ifelse(max(log10_load) <= VL_LOWER_LIMIT, NA_character_, 
                                infection_hash), by = infection_hash] %>%
    .[, ID := infection_hash] %>%
    .[age <= 0, age := 0.1] %>%
    .[, age_group := cut(age, breaks = AGE_breaks, ordered = TRUE, 
                        labels = AGE_labels)] %>%
    .[, date := as.Date(date, format = "%Y-%m-%d")] %>%
    .[, onset.date := as.Date(onset, format = "%Y-%m-%d")] %>%
    .[, keep_onset := ifelse(is.na(onset.date), FALSE, TRUE)] %>%
    .[is.na(onset.date), onset.date := as.Date(onset, format = "%Y-%m-%d")] %>%
    .[, onset := as.numeric(date - onset.date)] %>%
    .[, month := month(date)] %>%
    .[!is.na(infection_hash) & (log10_load > VL_LOWER_LIMIT), 
      first_pos_test_date := min(date), by = .(infection_hash)] %>%
    .[!is.na(infection_hash), first_pos_test_date := nafill(first_pos_test_date, type = "locf"), by = .(infection_hash)] %>%
    .[!is.na(infection_hash), first_pos_test_date := nafill(first_pos_test_date, type = "nocb"), by = .(infection_hash)] %>%
    .[, day_since_first_pos_test := as.numeric(date - first_pos_test_date)] %>%
    .[, machine := PCR] %>%
    .[PCR %in% c("T1", "T2", "T5", "T6", "T9"), machine := "cobas"] %>%
    .[PCR %in% c("Ta", "Tb", "Tc", "T3", "T4"), machine := "cepheid"] %>%
    .[, machine := factor(machine)] %>%
    .[test_centre_category == "C19" & hospitalized == 0, PAMS1 := 1] %>%
    .[, group := factor(ifelse(PAMS1 == 1, "PAMS", 
                               ifelse(hospitalized == 1, "hospitalized", "Other")),
                        levels = c("Other", "PAMS", "hospitalized")
    )] %>%
    .[, PAMS := ifelse(test_centre_category %in% PAMS_TEST_CENTRE, 1, 0)] %>%
    .[, `Age category` := cut(age,
                              breaks = AGE_BREAKS, ordered_result = TRUE,
                              labels = AGE_LABELS
    )] %>%
    .[, age_category2 := factor(cut(age, breaks = AGE_BREAKS2, 
                                   ordered_result = TRUE, right = FALSE),
                               levels = AGE_LEVELS2, labels = AGE_LABELS2, 
                               ordered = TRUE
    )] %>%
    .[, age_category2_code := as.integer(age_category2)] %>%
    .[, age_category3 := factor(cut(age, breaks = AGE_BREAKS3, 
                                   ordered_result = TRUE, right = FALSE),
                               levels = AGE_LEVELS3, labels = AGE_LABELS3, 
                               ordered = TRUE
    )] %>%
    .[, age_category3_code := as.integer(age_category3)] %>%
    .[, days_to_prior_infection := as.integer(days_to_prior_infection)] %>%
    .[, is_first_pos_pcr := ifelse(is_first_pos_pcr == "True", TRUE, FALSE)]
  
  
  # There should only be one initial negative and first positive PCR test per 
  # infection.
  assert_that(all(fulldt[!is.na(infection_hash)]$day_since_first_pos_test <= 0, 
                  na.rm = TRUE))
  
  # Remove rows with no date.
  fulldt <- fulldt[!is.na(date)] %>%
    setkeyv(c("age", "person_hash"))
  
  fulldt %>%
    .[, month := ordered(months(date), levels = unique(months(sort(fulldt$date))))]
  
  fulldt %>%
    .[, n_age := age + rnorm(1, sd = .5), by = 1:nrow(fulldt)] %>%
    .[, c_week := week(date)] %>%
    .[, centre_week := paste0(test_centre, c_week)] %>%
    .[, variant := stringr::str_to_title(variant)] %>%
    .[(log10_load > VL_LOWER_LIMIT) & is.na(variant), variant := "Unknown"] %>%
    .[, variant := factor(variant, labels = VARIANT_LABELS, 
                          levels = VARIANT_LEVELS)] %>%
    .[, variant2 := variant] %>%
    # Note: there is at least one infection with an Omicron label 
    # (by typing PCR) but it was not considered because it did not happen long 
    # enough after a previous infection, i.e. we cannot be sure that they are 
    # two separate infections. Therefore, there is no infection hash for this 
    # PCR. --> we need to exclude rows without an infection hash
    .[(variant == "Omicron") & !is.na(infection_hash), variant2 := add_omicron_sublineages(date, day_since_first_pos_test, variant, infection_hash), by = .(infection_hash)] %>%
    .[(log10_load > VL_LOWER_LIMIT) & is.na(variant2), variant2 := "Unknown"] %>%
    .[, variant2 := factor(variant2, labels = VARIANT2_LABELS, levels = VARIANT2_LEVELS)] %>%
    .[, wildtype_centre := ifelse(test_centre %in% unique(fulldt[variant == "Wildtype", test_centre]), TRUE, FALSE)] %>%
    .[, alpha_centre := ifelse(test_centre %in% unique(fulldt[variant == "Alpha", test_centre]), TRUE, FALSE)] %>%
    .[, delta_centre := ifelse(test_centre %in% unique(fulldt[variant == "Delta", test_centre]), TRUE, FALSE)] %>%
    .[, omicron_centre := ifelse(test_centre %in% unique(fulldt[variant == "Omicron", test_centre]), TRUE, FALSE)] %>%
    .[, ba1_centre := ifelse(test_centre %in% unique(fulldt[variant2 == "BA1", test_centre]), TRUE, FALSE)] %>%
    .[, ba2_centre := ifelse(test_centre %in% unique(fulldt[variant2 == "BA2", test_centre]), TRUE, FALSE)] %>%
    .[, ba5_centre := ifelse(test_centre %in% unique(fulldt[variant2 == "BA5", test_centre]), TRUE, FALSE)] %>%
    .[, ba2_desc_centre := ifelse(test_centre %in% unique(fulldt[variant2 == "BA2Desc", test_centre]), TRUE, FALSE)]
  
  return(fulldt)
}


#' Load viral load data from first positive PCR tests.
#' @return A data.table with viral load data from first positive PCRs.
load_first_pos_vl_data <- function() {
  dt <- get_first_pos_pcr_data(vl_first_pos_file = VL_FIRST_POS_FILE_TSV) %>%
    .[!is.na(log10_load) & (log10_load > 0)]
  
  return(dt)
}


#' Load and pre-process time course data.
#'
#' @return A data.table with viral load time course data.
get_TC_data <- function(one_infection_per_person = FALSE,
                        use_paper1_data = FALSE,
                        tc_file = TC_FILE_JSON,
                        vl_first_pos_file = VL_FIRST_POS_FILE_TSV) {
  TC.dt <- setDT(jsonlite::stream_in(file(tc_file)))
  summary_cols <- c("n_people", "n_infections", "min_pcr_days", "max_pcr_days", 
                    "max_infection_duration", "include_final_negative", 
                    "max_days_to_negative")
  
  # Extract summary statistics row.
  TC.summary <- TC.dt[1, ..summary_cols]
  # Get rid of the first line (summary statistics).
  TC.dt <- TC.dt[-1]
  # Remove all summary statistic columns.
  TC.dt[, (summary_cols) := NULL]
  # Expand columns with multiple values in one row.
  TC.dt <- setDT(TC.dt %>% tidyr::unnest(cols = c("viral_load", "test_name", 
                                                  "date", "test_centre", 
                                                  "test_centre_category", "age", 
                                                  "material", "ct", "id")))
  # Remove all data that was collected after data used in the first viral load 
  # paper.
  if (use_paper1_data) {
    TC.dt <- TC.dt[(date <= as.Date("2021-04-02")) & (test_name %in% c("T2", "LC480"))]
  }
  
  # Assert that within each infection, rows are ordered by date.
  assert_that(all(TC.dt[, order(date) == seq_len(.N), by = infection_hash]$V1))
  
  TC.dt <- TC.dt %>%
    .[, person_hash := factor(person_hash)] %>%
    .[, person_hash := factor(as.numeric(person_hash))] %>%
    .[, pcr_date := as.Date(date, format = "%Y-%m-%d")] %>%
    .[, day := as.numeric(pcr_date - min(pcr_date)), by = infection_hash] %>%
    .[, date := NULL] %>%
    .[, onset_json := as.Date(onset, format = "%Y-%m-%d")] %>%
    .[, onset := NULL] %>%
    .[, age := age] %>%
    .[, age := min(age), by = infection_hash] %>%
    .[age <= 100] %>%
    .[, age_category2 := factor(cut(age, breaks = AGE_BREAKS2, 
                                   ordered_result = TRUE, right = FALSE),
                               levels = AGE_LEVELS2, labels = AGE_LABELS2, 
                               ordered = TRUE)] %>%
    .[, age_category2_code := as.integer(age_category2)] %>%
    .[, age_category3 := factor(cut(age, breaks = AGE_BREAKS3, 
                                   ordered_result = TRUE, right = FALSE),
                               levels = AGE_LEVELS3, labels = AGE_LABELS3, 
                               ordered = TRUE)] %>%
    .[, age_category3_code := as.integer(age_category3)] %>%
    setnames("viral_load", "log10_load") %>%
    .[, N_tests := .N, by = infection_hash] %>%
    .[, Study := "BER"] %>%
    .[, variant := str_to_title(variant)] %>% # Capitalize first letter
    .[, variant := factor(variant, labels = VARIANT_LABELS, 
                          levels = VARIANT_LEVELS)] %>%
    .[is.na(variant), variant := "Unknown"] %>%
    .[, log10_load_orig := log10_load] %>%
    # At a log10 viral load <=1 the results are unreliable, i.e. some PCRs will 
    # be positive, some not.
    .[log10_load <= VL_LOWER_LIMIT, log10_load := 0] %>%
    # Remove infections where the maximum viral load is smaller than 4.
    .[, .SD[max(log10_load) > 4], by = infection_hash] %>%
    # Note: there is one case 
    # (infection_hash = "4cb54f26d489e512f5da029400d1be1") whose first positive 
    # viral load is just above 2 (and therefore it is included by the python 
    # script used to generate viral-load-with-negatives-paper2.tsv) but here in 
    # the R script, the trailing decimal places are not recorded anymore and 
    # thus this infection will be discarded because the viral load is labelled 
    # as not being above 2 in fulldt. So below, when using fulldt to fill up 
    # column "variant2", variant2 will NA for this infection. 
    # Infection "daca3ce68420a424612e1a29ba0d7a75" has a similar issue, so we 
    # discard it.
    .[!(infection_hash %in% c("4cb54f26d489e512f5da029400d1be14", 
                             "daca3ce68420a424612e1a29ba0d7a75"))]
  
  # Make sure that each infection has a day=0 specified
  assert_that(uniqueN(TC.dt$infection_hash) == uniqueN(TC.dt[day == 0]$infection_hash))
  
  fulldt <- get_first_pos_pcr_data(vl_first_pos_file = vl_first_pos_file)
  
  # If we want to keep only the latest infection per person.
  if (one_infection_per_person) {
    # Identify the date of each infection's first PCR.
    min_date_per_infection <- TC.dt[, .(min_date_infection = min(pcr_date)), 
                                    by = .(person_hash, infection_hash)]
    # Identify the latest infection for each person.
    last_infection <- min_date_per_infection[, .SD[which.max(min_date_infection)], 
                                             by = person_hash]
    # Filter the original data to keep only the rows associated with the last 
    # infection.
    TC.dt <- TC.dt[last_infection, on = .(person_hash, infection_hash), 
                   nomatch = 0]
  }
  
  TC.dt <- TC.dt %>% .[, ID := factor(infection_hash)]
  setkeyv(fulldt, "ID")
  setkeyv(TC.dt, "ID")
  
  TC.dt <- TC.dt %>%
    .[fulldt, `:=`(
      variant = variant,
      variant2 = variant2,
      gender = gender,
      PAMS1 = PAMS1,
      hospitalized = hospitalized,
      symptom_onset = onset.date,
      keep_onset = keep_onset
    )]
  
  # Make sure that all rows have a variant2 value.
  # This is a temporary fix: (2 infections don't fulfill the criterion 
  # currently)
  TC.dt <- TC.dt[!is.na(variant2)]
  assert_that(all(!is.na(TC.dt[["variant2"]])))
  
  # Compute the median of the 1% of the lowest "positive" viral loads.
  TC.dt <- TC.dt %>%
    .[, L_median := {
      non_zero_vl <- log10_load[log10_load > 0]
      sorted_vl <- sort(non_zero_vl)
      median(sorted_vl[1:ceiling(0.01 * length(sorted_vl))])}, by = test_name] %>%
    .[, L_median5 := {
      non_zero_vl <- log10_load[log10_load > 0]
      sorted_vl <- sort(non_zero_vl)
      median(sorted_vl[1:ceiling(0.05 * length(sorted_vl))])}, by = test_name] %>%
    .[, L := min(log10_load[log10_load > 0]), by = test_name] %>%
    .[, U := max(log10_load) + 0.1, by = test_name]
  
  return(TC.dt)
}


#' Filter PCR tests to retain key negative results.
#'
#' Keeps only positive PCR tests (log10_load > VL_LOWER_LIMIT) and the last 
#' leading/first trailing negative tests (log10_load ≤ VL_LOWER_LIMIT) relative 
#' to the first/last positive test.
#'
#' @param tmp data.table. Input data with columns: ID, day, log10_load, 
#' day_to_first_pos_pcr, day_to_last_pos_pcr.
#'
#' @return data.table. Filtered data retaining relevant negative and all 
#' positive PCR tests.
filter_last_leading_first_trailing_neg_pcrs <- function(tmp) {
  tmp <- tmp %>%
    .[, is_last_leading_neg_pcr := day == suppressWarnings(max(day[(log10_load <= VL_LOWER_LIMIT) & (day_to_first_pos_pcr < 0)])), by = ID] %>%
    .[, is_first_trailing_neg_pcr := day == suppressWarnings(min(day[(log10_load <= VL_LOWER_LIMIT) & (day_to_last_pos_pcr > 0)])), by = ID] %>%
    .[(log10_load > VL_LOWER_LIMIT) | is_last_leading_neg_pcr | is_first_trailing_neg_pcr]
  
  return(tmp)
}


#' Pre-process data for viral load time course analysis.
#'
#' @return A data.table with processed viral load time course data.
prep_time_course_data <- function(merged_data,
                                  filter_by_testname,
                                  exclude_lrt_samples,
                                  latest_peak_day = Inf,
                                  selection = 3,
                                  min_n_pos_tests = 2,
                                  min_max_load = 4,
                                  min_duration = 0,
                                  trim_neg_pcrs = FALSE,
                                  days_to_negative = 14,
                                  single_leading_trailing_neg_pcr = FALSE,
                                  simulate_data = FALSE,
                                  grp_var = GRP_VAR) {
  tmp <- merged_data
  
  if (!is.null(filter_by_testname)) {
    assert_that(filter_by_testname %in% unique(tmp$test_name))
    tmp <- tmp[test_name == filter_by_testname]
  }
  if (!is.null(exclude_lrt_samples) && (exclude_lrt_samples == TRUE)) {
    tmp <- tmp[(material != "LRT") | (is.na(material))]
  }
  # Remove all infections for which the maximum viral load is below 2.
  tmp <- tmp[, .SD[!(max(log10_load) <= VL_LOWER_LIMIT)], by = ID]
  tmp[, day := as.numeric(pcr_date - min(pcr_date)), by = ID]
  tmp <- tmp %>%
    .[, day_to_first_pos_pcr := day - min(day[log10_load > VL_LOWER_LIMIT]), by = ID] %>%
    .[, day_to_last_pos_pcr := day - max(day[log10_load > VL_LOWER_LIMIT]), by = ID]
  
  if (trim_neg_pcrs) {
    # At most days_to_negative days from negative to first positive PCR and
    # days_to_negative days from last positive to negative PCR
    tmp <- tmp[(day_to_first_pos_pcr >= -days_to_negative) & (day_to_last_pos_pcr <= days_to_negative)]
  }
  # If we are simulating data, we are first creating the original number of 
  # PCRs and are removing intermediate and leading/trailing negatives later.
  if (single_leading_trailing_neg_pcr & !simulate_data) {
    # of all negative PCRs, keep only the first leading and last trailing 
    # negative PCR test
    tmp <- filter_last_leading_first_trailing_neg_pcrs(tmp)
  }
  
  tmp <- tmp %>%
    .[, gender := as.numeric(factor(gender, levels = GENDER_LEVELS)) - 1] %>%
    .[, hospitalized := as.numeric(sum(test_centre_category %in% HOSPITAL_TEST_CENTRE) > 0), by = ID] %>%
    .[, phosptests := mean(test_centre_category %in% HOSPITAL_TEST_CENTRE), by = ID] %>%
    .[, last_test_negative := ifelse(tail(log10_load, 1) <= VL_LOWER_LIMIT, TRUE, FALSE), by = ID] %>%
    .[, first_test_negative := ifelse(head(log10_load, 1) <= VL_LOWER_LIMIT, TRUE, FALSE), by = ID] %>%
    .[, first_last_test_negative := ifelse(last_test_negative == TRUE & first_test_negative == TRUE, TRUE, FALSE), by = ID] %>%
    # Create new column for prior infection.
    .[, prior_infection := n_prior_infections > 0] %>%
    .[, variant := factor(variant, levels = VARIANT_LEVELS)] %>%
    .[, variant2 := factor(variant2, levels = VARIANT2_LEVELS)] %>%
    .[, multiple_infections := ifelse(length(unique(infection_hash)) > 1, TRUE, FALSE), by = person_hash] %>%
    .[, PAMS3 := ifelse(PAMS1 == 1, 1, ifelse(hospitalized == 1, 2, 3))] %>%
    # Time between earliest and latest positive PCR in an infection.
    .[, is_neg_test := log10_load <= VL_LOWER_LIMIT] %>%
    .[, is_pos_test := log10_load > VL_LOWER_LIMIT] %>%
    .[, N_tests := .N, by = ID] %>%
    .[N_tests >= selection] %>%
    .[, N_pos_tests := sum(log10_load > VL_LOWER_LIMIT), by = ID] %>%
    .[N_pos_tests >= min_n_pos_tests] %>%
    setkeyv(unique(c(grp_var, "ID")))
  
  tmp[, D := as.numeric(log10_load_orig > 0)]
  tmp[, day := as.numeric(pcr_date - min(pcr_date)), by = ID]
  
  tmp <- tmp %>%
    .[order(person_hash, infection_hash, day)] %>%
    .[, highest_vl_day := day[which.max(log10_load)], by = get(grp_var)] %>%
    .[, day2 := day - highest_vl_day] %>%
    .[, starts_with_negative_pcr := log10_load[day == 0] == 0, by = ID] %>%
    # The following is just meant as a check to see if the drop in shift
    # parameter values around day 0 is related to likelihood instability
    # around the peak; here, we check if the drop happens around -5 which
    # would support this hypothesis.
    .[, day3 := day + 5] %>%
    # Label the day of the first positive PCR as day 0.
    .[, day4 := day - min(day[log10_load > VL_LOWER_LIMIT]), by = ID] %>%
    .[, day_to_last_pos_pcr := day - max(day[log10_load > VL_LOWER_LIMIT]), by = ID] %>%
    .[, max_day := max(day), by = ID] %>%
    .[, max_load := max(log10_load), by = ID] %>%
    # Exclude infections where the log10 viral load never exceeds min_max_load.
    .[max_load >= min_max_load] %>%
    .[, shift_sd := 40 / N_pos_tests] %>%
    .[, shift_sd2 := 30 / N_pos_tests] %>%
    .[, shift_sd3 := 30 / N_pos_tests + (9 - max_load)] %>%
    .[, max_load_day := max((log10_load == max_load) * day), by = ID] %>%
    .[max_load_day < latest_peak_day] %>%
    .[, first_pos_day := min(day[log10_load > 0]), by = ID] %>%
    .[, last_pos_day := max(day[log10_load > 0]), by = ID] %>%
    .[, infection_duration := last_pos_day - first_pos_day + 1, by = ID] %>%
    .[, duration := max(day) - min(day) + 1, by = ID] %>%
    .[infection_duration >= min_duration] %>%
    .[, days := max_day - first_pos_day] %>%
    # Calculate differences in log10_load and day between consecutive data points
    .[, vl_diff := c(NA, diff(log10_load)), by = ID] %>%
    .[, day_diff := c(NA, diff(day)), by = ID] %>%
    .[, Test1_negative_malo := ifelse(first_test_negative == TRUE, max_load, 0)] %>%
    .[is.na(recent_prior_infection), recent_prior_infection := FALSE] %>%
    # Label tests where the corresponding person had multiple infections in the 
    # data.
    .[, multiple_infections := ifelse(length(unique(infection_hash)) > 1, TRUE, FALSE), by = person_hash] %>%
    .[, N_test_grp := cut(N_tests, breaks = c(0, 5, 15, 100), ordered_result = TRUE)]
  
  tmp[, onset_day := as.numeric(onset_json - min(pcr_date)), by = .(ID)]
  tmp[onset_day < -14, onset_day := NA]
  
  
  if (selection > 1) {
    # Calculate the maximum rate of change in log10_load per day
    tmp[, diff_load_perday := vl_diff / day_diff]
    tmp[, diff_load_perday_max := max(abs(vl_diff / day_diff), na.rm = TRUE), by = .(ID)]
    tmp[, ld_centre := get_ld_centre(test_centre_category, day), by = ID]
    tmp[is.na(ld_centre), ld_centre := "?"]
    tmp[ld_centre %in% names(which(table(tmp[day == 0, ld_centre]) < 21)), 
        ld_centre := "X"]
    # Commenting the below out because we are reassigning IDs further down.
    ld_centre_lvls <- c("WD", sort(setdiff(unique(tmp$ld_centre), "WD")))
    tmp[, ld_centre := factor(ld_centre, levels = ld_centre_lvls)]
    
    # Make sure the last day of a positive PCR happened after the first day of a
    # positive PCR.
    assert_that(all(tmp[, last_pos_day > first_pos_day]))
  }
  
  # Make sure that each infection has a day=0 specified
  assert_that(uniqueN(tmp$ID) == uniqueN(tmp[day == 0]$ID))
  # Make sure that within each infection, rows are sorted by date.
  assert_that(all(tmp[, order(pcr_date) == seq_len(.N), by = ID]$V1))
  return(filter_unusual_infections(tmp))
}


# Check here for most current variant info:
# https://data.lageso.de/lageso/corona/archiv/berlin-website-2023-09-26.html#abwasser 
# (there is also an archive for older info) --> contains only main sublineages
# https://public.data.rki.de/t/public/views/IGS_Dashboard/DashboardSublineages?%3Aembed=y&%3AisGuestRedirectFromVizportal=y
# and https://public.data.rki.de/t/public/views/IGS_Dashboard/DashboardVariants?%3Aembed=y&%3AisGuestRedirectFromVizportal=y

#' Assign Omicron sublineage based on first PCR date.
#'
#' @param pcr_dates Date vector. Dates of PCR tests for an Omicron infection.
#' @param days Numeric vector. Days relative to first positive PCR (0 = first 
#' positive).
#' @param variants Character vector. Variant names (must be "Omicron").
#' @param infection_hashes Character vector. Unique infection identifiers (for 
#' debugging).
#'
#' @return Character. Assigned Omicron sublineage (e.g., "BA1", "BA2", 
#' "JN1Desc") or "OmicronUnknown".
add_omicron_sublineages <- function(pcr_dates, days, variants, infection_hashes) {
  # Make sure that we are only dealing with Omicron infections.
  assert_that(all(variants == "Omicron", na.rm = TRUE))
  
  # Omicron sublineage time thresholds (>= 50% of samples from that sublineage).
  ba1_end <- as.Date("2022-02-28")
  ba2_start <- as.Date("2022-03-01")
  ba2_end <- as.Date("2022-06-06")
  ba5_start <- as.Date("2022-06-07")
  ba5_end <- as.Date("2023-02-13") # Note: this includes BA.5 and BA.5*
  ba2desc_start <- as.Date("2023-02-14") # This includes the XBBs and EG.5.1
  ba2desc_end <- as.Date("2023-11-20")
  jn1_start <- as.Date("2023-11-21")
  jn1_end <- as.Date("2024-05-20")
  jn1desc_start <- as.Date("2024-05-21")
  jn1desc_end <- as.Date("2024-10-21")
  xec_start <- as.Date("2024-10-22")
  xec_end <- as.Date("2025-04-07")
  other_start <- as.Date("2025-04-08")
  
  
  # It can be the case that several PCRs have been performed on the first day. 
  # --> select 1 entry
  first_pcr_date <- pcr_dates[days == 0][1]
  
  omicronSublineage <- NA_character_
  
  if (first_pcr_date < ba2_start) {
    return("BA1")
  } else if ((ba2_start <= first_pcr_date) & (first_pcr_date < ba5_start)) {
    return("BA2")
  } else if ((ba5_start <= first_pcr_date) & (first_pcr_date < ba2desc_start)) {
    return("BA5")
  } else if ((ba2desc_start <= first_pcr_date) & (first_pcr_date < jn1_start)) {
    return("BA2Desc")
  } else if ((jn1_start <= first_pcr_date) & (first_pcr_date < jn1desc_start)) {
    return("JN1")
  } else if ((jn1desc_start <= first_pcr_date) & (first_pcr_date < xec_start)) {
    return("JN1Desc")
  } else if ((xec_start <= first_pcr_date) & (first_pcr_date < other_start)) {
    return("XEC")
  } else if (other_start <= first_pcr_date) {
    return("OmicronUnknown")
  }
  return(omicronSublineage)
}


#' Create the list passed to the stan model for fitting.
#'
#' @param day_data A data.table with PCR test data.
#' @param grp_var Character. Column name to group by (e.g., infection ID).
#' @param variant_var Character. Column name specifying where variant 
#' information is stored.
#' @param ub_log_slope_up_mu Upper bound for the Stan model parameter 
#' log_slope_up_mu.
#' @param imputation_limit Log10 viral load upper bound for negative PCR tests.
#' @param model_no_stan The number/identifier of the stan file.
#' @param X_PGH A matrix with values of fixed effects (infection-group-specific) 
#' parameters.
#' @param X_ld_centre A matrix with values indicating the test center where most 
#' PCR were sampled (per infection).
#' @param use_true_vls A logical value indicating whether to use the true viral 
#' load values (i.e. also those below the imputation limit) when simulating 
#' viral load data.
#' @param no_shifts A logical value indicating whether to temporally shift 
#' infection data when simulating viral load data.
#' @param sim_data A logical value indicating whether the data were simulated.
#' @param use_paper1_data A logical value indicating whether to use only data 
#' used in https://www.science.org/doi/10.1126/science.abi5273.
#'
#' @return A named list.
make_datalist_day <- function(day_data,
                              grp_var,
                              variant_var,
                              ub_log_slope_up_mu,
                              imputation_limit,
                              model_no_stan,
                              X_PGH = NULL,
                              X_ld_centre = NULL,
                              use_true_vls = FALSE,
                              no_shifts = FALSE,
                              sim_data = FALSE,
                              use_paper1_data = FALSE) {
  vl_lower_limit <- ifelse(sim_data, VL_LOWER_LIMIT_SIM_DATA, VL_LOWER_LIMIT)
  # Precompute some variables so they can be reused.
  idx_GG <- which(day_data[day == 0]$multiple_infections == TRUE)
  idx_omicron <- which(day_data[day == 0]$variant == "Omicron")
  idx_non_omicron <- which(day_data[day == 0]$variant != "Omicron")
  
  day_data <- day_data %>%
    .[, L_infection := median(L), by = ID] %>%
    .[, L_infection_median := median(L_median), by = ID] %>%
    .[, N_pos_tests_infection := sum(log10_load > vl_lower_limit), by = ID]
  
  
  omicron_subvariant_levels <- OMICRON_SUBVARIANT_LEVELS
  
  
  if (variant_var == "variant") {
    variant_levels <- VARIANT_LEVELS
  } else {
    variant_levels <- VARIANT2_LEVELS
  }
  
  # For each infection, compute the number and position of negative tests
  # within that infection.
  neg_test_dt <- day_data %>%
    .[, .(
      indices = list(which(log10_load == 0)),
      count_neg = sum(log10_load == 0)
      ), by = ID]
  neg_indices <- unlist(neg_test_dt[count_neg > 0]$indices)
  cumulative_neg_counts <- cumsum(neg_test_dt$count_neg)
  neg_test_start <- c(1, head(cumulative_neg_counts, -1) + 1)
  neg_test_end <- cumulative_neg_counts
  neg_test_start[neg_test_dt$count_neg == 0] <- 0
  neg_test_end[neg_test_dt$count_neg == 0] <- 0
  # Indicator of whether an infection contains negative PCR tests.
  has_neg_tests <- day_data[day == 0]$ID %in% neg_test_dt[count_neg > 0]$ID
  
  
  if (use_true_vls == TRUE) {
    Y_DAY <- day_data$log10_load_true
  } else {
    Y_DAY <- day_data$log10_load
    Y_DAY_Orig <- day_data$log10_load_orig
  }
  
  if (no_shifts == TRUE) {
    days <- day_data$day_sampled
  } else {
    days <- day_data$day
    # Day with highest viral load in an infection is marked as day 0
    days2 <- day_data$day2
    # For stan model 65: first PCRs are marked as day 5 (i.e. 5 is added to 
    # "day" column)
    days3 <- day_data$day3
    # For stan models 74 and 179: first positive PCR is marked as day 0
    days4 <- day_data$day4
  }
  
  datalist_DAY <- list(
    condition_on_data = 1,
    U2 = day_data$U,
    U = day_data[log10_load > vl_lower_limit, U],
    has_neg_tests = has_neg_tests,
    G2 = cumsum(has_neg_tests),
    # We are assuming that we have 4 cores per chain available, i.e. running
    # 4 threads in parallel.
    grainsize = as.integer(ceiling(uniqueN(day_data[, get(grp_var)])) / 4),
    N_DAY = nrow(day_data),
    N_DAY_pos = sum(day_data$log10_load > vl_lower_limit),
    N_DAY_pos_orig = sum(day_data$log10_load_orig > 0),
    gstart_DAY = c(1, which(diff(as.numeric(day_data[, get(grp_var)])) != 0) + 1),
    gend_DAY = c(
      which(diff(as.numeric(day_data[, get(grp_var)])) != 0),
      nrow(day_data)
    ),
    gstart_DAY_pos = c(1, which(diff(as.numeric(day_data[Y_DAY > vl_lower_limit, 
                                                         get(grp_var)])) != 0) + 1),
    gend_DAY_pos = c(
      which(diff(as.numeric(day_data[Y_DAY > vl_lower_limit, get(grp_var)])) != 0),
      nrow(day_data[Y_DAY > vl_lower_limit])
    ),
    gstart_DAY_pos_orig = c(1, which(diff(as.numeric(day_data[Y_DAY_Orig > 0, 
                                                              get(grp_var)])) != 0) + 1),
    gend_DAY_pos_orig = c(
      which(diff(as.numeric(day_data[Y_DAY_Orig > 0, get(grp_var)])) != 0),
      nrow(day_data[Y_DAY_Orig > 0])
    ),
    gstart_DAY_neg = c(1, which(diff(as.numeric(day_data[Y_DAY == 0, 
                                                         get(grp_var)])) != 0) + 1),
    gend_DAY_neg = c(
      which(diff(as.numeric(day_data[Y_DAY == 0, get(grp_var)])) != 0),
      nrow(day_data[Y_DAY == 0])
    ),
    G = uniqueN(day_data[, get(grp_var)]),
    G_neg = uniqueN(day_data[log10_load <= vl_lower_limit, get(grp_var)]),
    N_subject = uniqueN(day_data[, person_hash]),
    subject = as.numeric(factor(day_data[day == 0, person_hash])),
    GG = sum(day_data[day == 0]$multiple_infections == TRUE),
    idx_GG = idx_GG,
    N_subject_GG = uniqueN(day_data[multiple_infections == TRUE]$person_hash),
    idx_subject_GG = as.numeric(factor(day_data[day == 0][idx_GG]$person_hash)),
    Y_DAY = Y_DAY,
    Y_DAY_Orig = Y_DAY_Orig,
    X_DAY = days,
    Y_DAY_pos = Y_DAY[Y_DAY > vl_lower_limit],
    X_DAY_pos = days[Y_DAY > vl_lower_limit],
    Y_DAY_pos_orig = Y_DAY_Orig[Y_DAY_Orig > 0],
    X_DAY_pos_orig = days[Y_DAY_Orig > 0],
    X_DAY_neg = days[Y_DAY == 0],
    X_DAY2 = days2,
    X_DAY2_pos = days2[Y_DAY > vl_lower_limit],
    X_DAY2_neg = days2[Y_DAY == 0],
    X_DAY3_pos = days3[Y_DAY > vl_lower_limit],
    X_DAY3_neg = days3[Y_DAY == 0],
    X_DAY4 = days4,
    X_DAY4_pos = days4[Y_DAY > vl_lower_limit],
    X_DAY4_neg = days4[Y_DAY == 0],
    starts_with_negative_PCR = day_data[day == 0]$starts_with_negative_pcr,
    G_starts_with_negative_PCR = sum(day_data[day == 0]$starts_with_negative_pcr),
    N_pos_tests_infection = day_data[day == 0]$N_pos_tests_infection,
    idx_omicron = idx_omicron,
    idx_non_omicron = idx_non_omicron,
    N_omicron_infections = length(idx_omicron),
    N_non_omicron_infections = length(idx_non_omicron),
    N_omicron = length(omicron_subvariant_levels),
    N_non_omicron = length(NON_OMICRON_LEVELS),
    omicron_subvariant = factor(day_data[day == 0][idx_omicron]$variant2,
                                levels = filter_levels(
                                  day_data[day == 0],
                                  "variant2",
                                  omicron_subvariant_levels
                                )
    ),
    non_omicron_variant = factor(day_data[day == 0][idx_non_omicron]$variant,
                                 levels = filter_levels(
                                   day_data[day == 0],
                                   "variant",
                                   NON_OMICRON_LEVELS
                                 )
    ),
    N_variant = nlevels(day_data[, get(variant_var)]),
    variant = factor(day_data[day == 0, get(variant_var)],
                     levels = filter_levels(
                       day_data[day == 0],
                       variant_var, variant_levels
                     )
    ),
    N_T1_neg = sum(day_data[day == 0, first_test_negative] == TRUE),
    T1_neg_idx = which(day_data[day == 0, first_test_negative] == TRUE),
    N_T1_pos = sum(day_data[day == 0, first_test_negative] == FALSE),
    T1_pos_idx = which(day_data[day == 0, first_test_negative] == FALSE),
    N_neg_tests = sum(day_data$log10_load == 0),
    N_neg_tests_orig = sum(day_data$log10_load_orig == 0),
    idx_neg_tests = which(day_data$log10_load == 0),
    idx_neg_tests_orig = which(day_data$log10_load_orig == 0),
    N_neg_tests_infection = neg_test_dt$count_neg,
    idx_neg_tests_infection = neg_indices,
    start_neg_tests = neg_test_start,
    end_neg_tests = neg_test_end,
    imputation_limit = imputation_limit,
    PAMS3 = day_data[day == 0, PAMS3],
    prior_infection = as.numeric(factor(day_data[day == 0, prior_infection])),
    max_load = as.numeric(scale(day_data[day == 0, max_load])),
    ub_log_slope_up_mu = ub_log_slope_up_mu,
    # 1 - lambda is the probability of a false negative PCR (model 58)
    lambda = 0.95, 
    # The exponentiated rate for the exponential distribution in the mixture 
    # likelihood for negative PCRs (model 58)
    exprate = 10, 
    # The standard deviation for the shift parameter, depending on the number 
    # of positive PCR tests in the corresponding infection
    shift_sd = day_data[day == 0, shift_sd], 
    shift_sd2 = day_data[day == 0, shift_sd2]
  )
  gender <- day_data[day == 0, gender]
  
  datalist_DAY2 <- list(
    N_centres = uniqueN(day_data$test_centre_category),
    centre = model.matrix(~ 0 + test_centre_category, day_data),
    N_centre1 = uniqueN(day_data[day == 0, test_centre_category]),
    centre1 = as.numeric(factor(day_data[day == 0, test_centre_category])),
    N_ld_centre = uniqueN(day_data[day == 0, ld_centre]),
    ld_centre = as.numeric(factor(day_data[day == 0, ld_centre])),
    K_ld_centre = ifelse(is.null(ncol(X_ld_centre)), 1, ncol(X_ld_centre)),
    X_ld_centre = X_ld_centre,
    # PCR = 1*(day_data$test_name %in% c("T1", "T2")),
    N_PCR = uniqueN(day_data$test_name),
    N_onset = sum(!is.na(day_data[day == 0, onset_day])),
    idx_onset = which(!is.na(day_data[day == 0, onset_day])),
    onset = day_data[day == 0 & !is.na(onset_day), onset_day]
  )
  
  # Change PCR variable (only T2 and LC480) when using data from viral load 
  # paper
  if (use_paper1_data) {
    datalist_DAY2$PCR <- 1 * (day_data$test_name == "T2")
  }
  if (uniqueN(day_data$test_name) > 1) {
    if (!is.null(sim_data) || sim_data == FALSE) {
      test_nameLevels <- filter_levels(day_data, "test_name", 
                                      TEST_NAME_LEVELS_SIM_DATA)
    } else {
      test_nameLevels <- filter_levels(day_data, "test_name", TEST_NAME_LEVELS)
    }
    day_data[, test_name := factor(test_name, levels = test_nameLevels)]
    if (uniqueN(day_data$test_name) > 2) {
      # can then be used for random effects parameterization
      datalist_DAY2$PCR <- model.matrix(~ 0 + test_name, day_data)
      datalist_DAY2$PCR_pos <- as.numeric(day_data[log10_load > vl_lower_limit]$test_name)
      datalist_DAY2$PCR_neg <- as.numeric(day_data[log10_load <= vl_lower_limit]$test_name)
    } else {
      # binary fixed effects with 0 (reference) or 1
      datalist_DAY2$PCR <- model.matrix(~test_name, day_data)[, -1]
      datalist_DAY2$PCR_pos <- as.numeric(day_data[log10_load > vl_lower_limit]$test_name) - 1
      datalist_DAY2$PCR_neg <- as.numeric(day_data[log10_load <= vl_lower_limit]$test_name) - 1
    }
  }
  # Add gender variable.
  genderlist <- list(
    Gender = do.call(c, lapply(gender, function(x) ifelse(is.na(x), 0.5, x)))
  )
  
  assert_that(all(day_data$L_median == day_data$L_median[1]))
  L_list <- list(
    L_median = day_data$L_median[1],
    L_infection_median = day_data[day == 0, L_infection_median],
    L_min = day_data$L[1],
    L_median5 = day_data$L_median5[1]
  )
  
  if (!is.null(X_PGH)) {
    xpgh_list <- list(
      X_PG = X_PGH[, 1:2],
      K_PGH = ncol(X_PGH),
      X_PGH = X_PGH
    )
  } else {
    xpgh_list <- list()
  }
  agelist <- list(
    N_AgeGrp = n_distinct(day_data$age_category3_code),
    AgeGrp = day_data[day == 0, age_category3_code]
  )
  
  if ("day_true" %in% colnames(day_data)) {
    days_true <- day_data$day_true
    day_true_list <- list(
      X_DAY_true = days_true,
      X_DAY_true_pos = days_true[Y_DAY > vl_lower_limit],
      X_DAY_true_neg = days_true[Y_DAY == 0]
    )
    datalist_DAY <- c(
      datalist_DAY, datalist_DAY2, day_true_list, genderlist,
      xpgh_list, agelist, L_list
    )
  } else {
    datalist_DAY <- c(
      datalist_DAY, datalist_DAY2, genderlist, xpgh_list, agelist,
      L_list
    )
  }
  
  return(datalist_DAY)
}


#' Sample infection data.
#'
#' @param samples An integer specifying the number of infections to randomly 
#' select.
#' @param day_data PCR test data.
#' @param model_no_stan The number/identifier of the stan file.
#'
#' @return A data.table containing data from a sub-sample of infections present 
#' in `day_data`.
sample_day_data <- function(samples, day_data, model_no_stan) {
  day_data <- day_data %>%
    .[ID %in% sample(day_data[day == 0, ID], samples, replace = FALSE)] %>%
    .[, N_tests := .N, by = .(ID)] %>%
    .[, min_day := min(day), by = .(ID)] %>%
    .[, day := day - min_day] %>%
    .[, max_load := max(log10_load), by = ID] %>%
    .[, max_load_day := max((log10_load == max_load) * day), by = ID] %>%
    .[, multiple_infections := ifelse(length(unique(infection_hash)) > 1, 
                                      TRUE, FALSE), by = person_hash] %>%
    .[, test_name := factor(test_name)] %>%
    .[, test_centre_category := factor(test_centre_category)] %>%
    .[, test_centre := factor(test_centre)] %>%
    .[, ld_centre := factor(ld_centre)] %>%
    .[, variant := factor(variant)] %>%
    .[, onset_day := as.numeric(onset_json - min(pcr_date)), by = .(ID)] %>%
    .[onset_day < -14, onset_day := NA]
  
  return(day_data)
}


#' Filter out infection courses with viral load increases above a specified 
#' threshold.
#'
#' @param day_data A data.table of daily infection data.
#' @param max_diff_load_perday Maximum allowed daily viral load increase (NULL 
#' to skip filtering).
#'
#' @return Filtered data frame with steep viral load increases removed.
exclude_too_steep_courses <- function(day_data, max_diff_load_perday) {
  if (!is.null(max_diff_load_perday)) {
    exclude_ID <- day_data[diff_load_perday > max_diff_load_perday & day4 == 0 & !is.na(diff_load_perday), ID]
    if (length(exclude_ID) > 0) {
      # remove the datapoint before the first positive PCR
      exclude_idx <- which(day_data$ID %in% exclude_ID & day_data$day4 == 0) - 1 
      day_data <- day_data[-exclude_idx]
      day_data %>%
        .[, N_tests := .N, by = .(ID)] %>%
        .[, min_day := min(day), by = .(ID)] %>%
        .[, day := day - min_day] %>%
        .[, max_load := max(log10_load), by = ID] %>%
        .[, max_load_day := max((log10_load == max_load) * day), by = ID] %>%
        .[, multiple_infections := ifelse(length(unique(infection_hash)) > 1, 
                                          TRUE, FALSE), by = person_hash]
      if (sum(!is.na(day_data$onset_json)) > 0) {
        day_data %>%
          .[, onset_day := as.numeric(onset_json - min(pcr_date)), by = .(ID)] %>%
          .[onset_day < -14, onset_day := NA]
      }
    }
  }
  return(day_data)
}


#' Assert identity of core columns in two data.tables.
#'
#' Checks if `dt1` and `dt2` have identical values for columns: ID, person_hash, 
#' and day. Throws an error if discrepancies are found.
#'
#' @param dt1 data.table. First data.table to compare.
#' @param dt2 data.table. Second data.table to compare.
assert_datatables_identity <- function(dt1, dt2) {
  cols <- c("ID", "person_hash", "day")
  # Step 1: Subset both data.tables to include only common columns
  dt1_common <- dt1[, ..cols, with = FALSE]
  dt2_common <- dt2[, ..cols, with = FALSE]
  
  # order_cols = lapply(common_cols, as.name)
  # Step 2: Sort both data.tables by the common columns
  setorderv(dt1_common, cols = cols)
  setorderv(dt2_common, cols = cols)
  
  # Step 3: Compare the sorted data.tables
  assertthat::assert_that(identical(dt1_common, dt2_common))
}


#' Generate a list with data for the Stan model for time course analysis.
#'
#' @return A list.
make_time_course_standata <- function(selection = 3,
                                      max_n_tests = 20,
                                      min_n_pos_tests = 2,
                                      samples = NULL,
                                      max_diff_load_perday = NULL,
                                      ub_log_slope_up_mu = Inf,
                                      imputation_limit = 3,
                                      latest_peak_day = Inf,
                                      remove_some_onsets = TRUE,
                                      model_no_stan = NULL,
                                      tc_file = TC_FILE_JSON,
                                      vl_first_pos_file = VL_FIRST_POS_FILE_TSV,
                                      simulate_data = FALSE,
                                      simulate_param_sizes_all_zero = FALSE,
                                      simulate_three_slope_model = FALSE,
                                      simulate_infection_wise_alpha = FALSE,
                                      simulate_mixed_courses = FALSE,
                                      error_dist_sim_data = "normal",
                                      use_true_vls = FALSE,
                                      trim_neg_pcrs = FALSE,
                                      days_to_negative = 14,
                                      single_leading_trailing_neg_pcr = FALSE,
                                      add_neg_pcrs = FALSE,
                                      no_shifts = FALSE,
                                      remove_problematic_infections = FALSE,
                                      filter_by_testname = "T2",
                                      only_pos_pcrs_sim_data = FALSE,
                                      infection_wise_vl_error = TRUE,
                                      threading = NULL,
                                      use_paper1_data = FALSE,
                                      exclude_lrt_samples = FALSE,
                                      min_duration = 0,
                                      min_max_load = 4,
                                      one_infection_per_person = TRUE,
                                      variant_var = "variant",
                                      grp_var = GRP_VAR) {
  assert_that(selection >= min_n_pos_tests)
  
  if (!is.null(samples)) {
    assertthat::assert_that(!is.null(model_no_stan))
  }
  
  if (!is.null(filter_by_testname)) {
    assertthat::assert_that(filter_by_testname %in% TEST_NAME_LEVELS)
    if (filter_by_testname == "T1") {
      assert_that(tc_file == TC_FILE_JSON_N_GENE)
      assert_that(vl_first_pos_file == VL_FIRST_POS_FILE_TSV_N_GENE)
    } else {
      if (filter_by_testname == "T2") {
        assert_that(tc_file == TC_FILE_JSON)
        assert_that(vl_first_pos_file == VL_FIRST_POS_FILE_TSV)
      }
    }
  }
  assert_that(!is.null(filter_by_testname))
  
  
  if (use_true_vls == TRUE) {
    assertthat::assert_that(simulate_data == TRUE,
                            msg = "use_true_vls can only be true if working with
                            simulated data"
    )
  }
  TC_data <- get_TC_data(
    one_infection_per_person = one_infection_per_person,
    use_paper1_data = use_paper1_data,
    tc_file = tc_file,
    vl_first_pos_file = vl_first_pos_file
  )
  
  
  # Remove infections for which the time shift parameter did not converge in 
  # previous fits (Rhat > 1.01, model17_17, data from March 2024).
  if (remove_problematic_infections == TRUE) {
    old_length <- nrow(TC_data)
    problematic_infections_dt <- fread(PROBLEMATIC_INFECTIONS_FILE)
    problematic_infections <- problematic_infections_dt$ID
    TC_data <- TC_data[!ID %in% problematic_infections]
    new_length <- nrow(TC_data)
    assert_that(old_length > new_length)
  }
  
  day_data <-
    prep_time_course_data(TC_data,
                          filter_by_testname = filter_by_testname,
                          exclude_lrt_samples = exclude_lrt_samples,
                          latest_peak_day = latest_peak_day,
                          selection = selection,
                          min_n_pos_tests = min_n_pos_tests,
                          min_max_load = min_max_load,
                          min_duration = min_duration,
                          trim_neg_pcrs = trim_neg_pcrs,
                          days_to_negative = days_to_negative,
                          single_leading_trailing_neg_pcr = single_leading_trailing_neg_pcr,
                          simulate_data = simulate_data,
                          grp_var = grp_var
    )
  
  if (!is.null(filter_by_testname)) {
    assert_that(uniqueN(day_data$test_name) == 1)
  }
  
  if (remove_some_onsets == TRUE) {
    day_data[keep_onset == FALSE, onset_day := NA]
  }
  
  day_data <- exclude_too_steep_courses(day_data,
                                        max_diff_load_perday = max_diff_load_perday
  )
  
  
  if (!is.null(samples) & simulate_data == FALSE) {
    day_data <- sample_day_data(samples, day_data, model_no_stan)
  }
  
  centre_matrix <- model.matrix(~ 0 + test_centre_category, day_data[day == 0])
  
  # Adding negative PCRs to the beginning and end of certain infections.
  # Note: we want to keep negative days so that we don't have to adjust the
  # limits of the shift parameter in the stan models.
  if (add_neg_pcrs == TRUE) {
    day_data <- add_neg_PCRs_infection(day_data)
  }
  
  # Make sure the data table is sorted by the variable according to which we
  # want to group.
  day_data <- day_data %>% .[, (grp_var) := factor(get(grp_var))]
  day_data <- day_data[with(day_data, order(get(grp_var), day, person_hash)), ]
  
  if (simulate_data == FALSE) {
    sim_data <- NULL
  } else {
    sim_data <- sim_vl_courses(
      day_data = day_data,
      n = samples,
      sigma_up_slope = 0.4,
      param_sizes_all_zero = simulate_param_sizes_all_zero,
      min_n_tests = selection,
      max_n_tests = max_n_tests,
      min_n_pos_tests = min_n_pos_tests,
      days_to_negative = days_to_negative,
      trim_neg_pcrs = trim_neg_pcrs,
      add_neg_pcrs = add_neg_pcrs,
      only_pos_pcrs = only_pos_pcrs_sim_data,
      infection_wise_vl_error = infection_wise_vl_error,
      error_dist = error_dist_sim_data,
      three_slope_model = simulate_three_slope_model,
      infection_wise_alpha = simulate_infection_wise_alpha,
      mixed_courses = simulate_mixed_courses,
      min_max_load = min_max_load,
      min_duration = min_duration,
      single_leading_trailing_neg_pcr = single_leading_trailing_neg_pcr
    )
    # Assert that the row order is correct.
    expected_order <- order(
      sim_data$dt$infection_hash, sim_data$dt$day,
      sim_data$dt$person_hash
    )
    assert_that(identical(seq_len(nrow(sim_data$dt)), expected_order))
    
    sim_data <- list(
      dt = sim_data$dt, intercept = unique(sim_data$dt$intercept),
      up_slope = unique(sim_data$dt$up_slope),
      down_slope = unique(sim_data$dt$down_slope),
      sigma = unique(sim_data$dt$sigma),
      params = sim_data$params
    )
    assertthat::assert_that(length(sim_data$intercept) == uniqueN(sim_data$dt$ID))
    assertthat::assert_that(length(sim_data$up_slope) == uniqueN(sim_data$dt$ID))
    assertthat::assert_that(length(sim_data$down_slope) == uniqueN(sim_data$dt$ID))
    day_data <- sim_data$dt
  }
  
  # Someone is categorized as PAMS1 if their very first or their first positive
  # PCR was taken at a mild test center and they were never hospitalized
  # Prepare data for model variables.
  X_PGH_data <- day_data[day == 0] %>%
    .[, gender := ifelse(is.na(gender), 0.5, gender)] %>%
    .[, age := scale(age)] %>%
    .[, PAMS2 := ifelse((PAMS1 == 0) & (hospitalized == 0), 0.5, PAMS1)] %>%
    .[, age_cat2 := as.numeric(age_category3_code == 2)] %>%
    .[, age_cat3 := as.numeric(age_category3_code == 3)]
  
  cols.num <- c("PAMS1", "prior_infection", "hospitalized", "age_cat2", "age_cat3")
  X_PGH_data <- X_PGH_data %>%
    .[, (cols.num) := lapply(.SD, as.numeric), .SDcols = cols.num]
  X_PGH <- model.matrix(~ 1 + PAMS1 + gender + hospitalized + prior_infection + age_cat2 + age_cat3 + variant + variant:hospitalized, X_PGH_data)
  
  
  # Remove intercept.
  X_PGH <- X_PGH[, -1]
  # Remove columns with only 0s.
  X_PGH <- X_PGH[, colSums(X_PGH) > 0]
  
  X_ld_centre <- model.matrix(~ld_centre, day_data[day == 0])[, -1]
  
  datalist <- make_datalist_day(day_data,
                                grp_var = grp_var,
                                variant_var = variant_var,
                                ub_log_slope_up_mu = ub_log_slope_up_mu,
                                imputation_limit = imputation_limit,
                                model_no_stan = model_no_stan,
                                X_PGH = X_PGH,
                                X_ld_centre = X_ld_centre,
                                use_true_vls = use_true_vls,
                                no_shifts = no_shifts,
                                sim_data = simulate_data,
                                use_paper1_data = use_paper1_data
  )
  
  assert_that(length(unique(names(datalist))) == length(names(datalist)))
  
  assertthat::assert_that(datalist$G == length(levels(day_data[, ID])))
  return(list(datalist = datalist, day_data = day_data, sim_data = sim_data))
}


#' Generate a table with age group counts.
#'
#' Generates a table of age group counts, optionally filtered by a specific 
#' group. Age groups are defined by `breaks`, and ages are ceiling-rounded 
#' (max 100).
#'
#' @param bdata A data.table. Input data with columns: age, group.
#' @param breaks Numeric vector. Breakpoints for age groups 
#' (e.g., `seq(0, 100, 10)`).
#' @param my_group Character vector. Optional subset of `group` values to 
#' include.
#' @return A data.table. Age group counts with columns: age_group, age, N_age.
make_age_group_N <- function(bdata, breaks = NULL, my_group = NULL) {
  complete_age_table <-
    bdata[, .(age)] %>%
    .[, age := ceiling(age)] %>%
    .[age == 101, age := 100] %>%
    .[, age_group := cut(age, breaks = breaks, right = FALSE)] %>%
    unique()
  
  tmpdata <- bdata[, .(age, group)]
  if (!is.null(my_group)) {
    tmpdata <- tmpdata[group %in% my_group]
  }
  
  age_group_N <-
    tmpdata %>%
    .[, age := ceiling(age)] %>%
    .[, age_group := cut(age, breaks = breaks, right = FALSE)] %>%
    .[, .(N_age = sum(.N)), by = .(age_group, age)]
  
  setkeyv(complete_age_table, names(complete_age_table))
  setkeyv(age_group_N, names(complete_age_table))
  age_group_N <-
    complete_age_table[age_group_N, N_age := N_age] %>%
    .[is.na(N_age), N_age := 0]
  return(age_group_N)
}


#' Initial values for estimation of time course model in Stan.
#'
#' @param datalist A list with data for Stan model.
#' @return A list with initial values.
make_TC_inits <-
  function(datalist) {
    param_list <- list(
      alpha = runif(1, 0.1, 1),
      alpha_mu = runif(1, 1, 4),
      beta_sweight_mu = rnorm(1, 10, .5),
      log_intercept_mu = rnorm(1, log(8), .05),
      log_slope_up_mu = rnorm(1, .6, .025),
      log_slope_down_mu = rnorm(1, -1.75, .05),
      intercept_sigma = runif(1, .0, .05),
      intercept_raw = rnorm(datalist$G, 0, .5),
      intercept_nu = runif(1, 2, 10),
      slope_up_sigma = runif(1, .3, .5),
      slope_up_raw = rnorm(datalist$G, 0, .5),
      slope_up_nu = runif(1, 2, 10),
      slope_down_sigma = runif(1, .0, .05),
      slope_down_raw = rnorm(datalist$G, 0, .5),
      slope_down_nu = runif(1, 2, 10),
      intercept_slope_down_raw = matrix(rnorm(2 * datalist$G, 0, .5),
                                        nrow = 2, ncol = datalist$G
      ),
      intercept_slope_down_sigma = c(
        runif(1, .0, .05),
        runif(1, .0, .05)
      ),
      shift = rep(0, datalist$G),
      b_shift = rep(0, datalist$G),
      shift_neg_infections = rep(0, datalist$G),
      shift_pos_infections = rep(0, datalist$G),
      shift_sigma = runif(1, 1, 5),
      shift_mu = runif(1, -1, 3),
      shift_raw = rnorm(datalist$G, 0, 1),
      theta = rbeta(1, 3, 1),
      lambda = runif(1, 0.4, .7),
      mu = rnorm(1, 0, .5),
      sigma_mu = rnorm(1, 0, .5),
      sigma_sigma = runif(1, .1, .2),
      sigma_raw = rnorm(datalist$G, 0, .5),
      sigma = runif(1, 0.1, 2),
      sigma_qpcr = runif(1, 0.1, 0.5),
      sigma_beta = runif(1, 0.1, 0.5),
      beta_gumbel_mu = rnorm(1, 0, .5),
      beta_gumbel_sigma = runif(1, .1, .2),
      beta_gumbel_raw = rnorm(datalist$G, 0, .5),
      omega_mu = rnorm(1, 0, .5),
      omega_sigma = runif(1, .1, .2),
      omega_raw = rnorm(datalist$G, 0, .5),
      shape_mu = rnorm(1, 0, .5),
      shape_nu = runif(1, 2, 5),
      shape_tau = runif(1, 0.5, 1.5),
      shape_raw = rnorm(datalist$G, 0, .5),
      shape_sigma = runif(1, .1, .2),
      nu = runif(1, 5, 30),
      imp_neg = array(runif(datalist$N_neg_tests, -1, 2),
                      dim = c(1, datalist$N_neg_tests)
      ),
      intercept_sigma_subject = runif(1, 0, .05),
      intercept_raw_subject = rnorm(datalist$N_subject_GG, 0, .5),
      slope_down_sigma_subject = runif(1, 0, .05),
      slope_down_raw_subject = rnorm(datalist$N_subject_GG, 0, .5),
      beta_gender_intercept = rnorm(1, 0, .01),
      beta_gender_slope_up = rnorm(1, 0, .01),
      beta_gender_slope_down = rnorm(1, 0, .01),
      nu_slope_down = rtruncnorm(1, a = 0, b = Inf, mean = 10, sd = 1),
      age_intercept = rnorm(datalist$N_AgeGrp, 0, 0.05),
      age_slope_up = rnorm(datalist$N_AgeGrp, 0, 0.05),
      age_slope_down = rnorm(datalist$N_AgeGrp, 0, 0.05),
      centre_sigma = runif(1, .1, .6),
      centre_raw = rnorm(datalist$N_centres, 0, .5),
      intercept_PCR_raw = rnorm(datalist$N_PCR, 0, .5),
      intercept_PCR_sigma = runif(1, 0, .5),
      int_centre1_sigma = runif(1, .1, .5),
      int_centre1_raw = rnorm(datalist$N_centre1, 0, .5),
      slope_down_ld_centre_sigma = runif(1, .01, .05),
      slope_down_ld_centre_raw = rnorm(datalist$N_ld_centre, 0, .5),
      shift_centre1_sigma = runif(1, .1, .5),
      shift_centre1_raw = array(rnorm(datalist$N_centre1, 0, .5), dim = 1),
      betaPGH_slope_down = rnorm(ncol(datalist$X_PGH), 0, .1),
      betaPGH_intercept = rnorm(ncol(datalist$X_PGH), 0, .05),
      betaPGH_slope_up = rnorm(ncol(datalist$X_PGH), 0, .05)
    )
    assert_that(length(unique(names(param_list))) == length(names(param_list)))
    return(param_list)
  }


#### Simulate data ####

#' Propagate first non-NA value to all rows in a group.
#'
#' Fills NA values in a vector with the first non-NA value encountered.
#' Useful for propagating group-level information to all rows in the group.
#'
#' @param col Vector. Input column, possibly containing NA values.
#' @return Vector. Input vector with NAs replaced by the first non-NA value.
propagate_to_group_rows <- function(col) {
  return(col[!is.na(col)][1])
}


#' Simulate viral load trajectories.
#'
#' Computes mean or stochastic viral load values for given days, using either a 
#' 2-slope or 3-slope model.
#' Optionally adds PCR/centre-specific intercepts and random error (normal, 
#' skew-normal, or student-t).
#'
#' @param days Numeric vector. Days relative to infection.
#' @param beta_sweight_mu Numeric. Slope weight parameter (logistic transition).
#' @param intercept Numeric. Baseline viral load intercept.
#' @param up_slope Numeric. Slope for viral load increase.
#' @param down_slope Numeric. Slope for viral load decrease.
#' @param intercept_pcr Numeric vector. PCR-type-specific intercepts.
#' @param PCRs Character vector. PCR types for each day.
#' @param intercept_centre Numeric vector. Centre-specific intercepts.
#' @param centres Character vector. Testing centres for each day.
#' @param sigma Numeric. Standard deviation for random error.
#' @param param_sizes_all_zero Logical. If `TRUE`, ignores PCR/centre intercepts.
#' @param alpha Numeric. Smoothness parameter for 2/3-slope models.
#' @param three_slope_model Logical. If `TRUE`, uses a 3-slope model.
#' @param error_dist Character. Error distribution: "normal", "skew-normal", or 
#' "student-t".
#' @param alpha_skew_normal Numeric. Skewness parameter for skew-normal 
#' distribution.
#'
#' @return Numeric vector. Simulated viral load values.
compute_viral_loads <- function(days,
                                beta_sweight_mu,
                                intercept,
                                up_slope,
                                down_slope,
                                intercept_pcr,
                                PCRs,
                                intercept_centre,
                                centres,
                                sigma,
                                param_sizes_all_zero = FALSE,
                                alpha = NULL,
                                three_slope_model = FALSE,
                                error_dist = "normal",
                                alpha_skew_normal = -2) {
  if (three_slope_model) {
    assert_that(!is.null(alpha))
    params <- compute_smooth_3_slope_params(intercept, up_slope, down_slope,
                                            L = 2
    )
    vls_mean <- yhat_smooth_3_slopes(
      days, params$a, params$b, params$c1,
      params$c2, params$beta2, alpha
    )
  } else {
    if (is.null(alpha)) {
      weightDown <- plogis(days * beta_sweight_mu)
      vlUp <- c(intercept) + c(up_slope) * days
      vlDown <- c(intercept) + c(down_slope) * days
      vls_mean <- (1 - weightDown) * vlUp + weightDown * vlDown
    } else {
      params <- compute_smooth_2_slope_params(intercept, up_slope, down_slope)
      vls_mean <- yhat_smooth_2_slopes(
        days, params$a, params$b, params$c,
        alpha
      )
    }
  }
  if (param_sizes_all_zero == FALSE) {
    vls_mean <- (vls_mean + intercept_pcr[PCRs] + intercept_centre[centres])
  }
  if (is.null(error_dist)) {
    return(vls_mean)
  } else {
    if (error_dist == "normal") {
      vls <- rnorm(n = length(days), mean = vls_mean, sd = sigma)
      return(vls)
    } else if (error_dist == "skew-normal") {
      vls <- rskew_normal(
        n = length(days), xi = vls_mean, omega = sigma,
        alpha = alpha_skew_normal
      )
      return(vls)
    } else if (error_dist == "student-t") {
      return(rstudent_t(n = length(days), df = 5, mu = vls_mean, sigma = sigma))
    }
  }
}


#' Count valid PCR tests in an infection.
#' Make sure negative PCRs happened at most days_to_negative days before or after
#' the first/last positive PCR test.
#'
#' @param vls A vector of log10 viral load values for an infection.
#' @param days A vector providing the corresponding days of testing (0 - 40
#' for simulated data).
#' @param days_to_negative The maximum number of days between a leading or
#' trailing negative PCR to the closest positive PCR of the same infection.
#'
#' @return An integer specifying the number of PCRs in the infection after
#' trimming.
find_valid_PCRs <- function(vls, days, days_to_negative) {
  # If there is at least one positive PCR
  maks_pos_pcrs <- vls > VL_LOWER_LIMIT_SIM_DATA
  if (sum(maks_pos_pcrs) > 0) {
    first_pos_pcr_idx <- which(maks_pos_pcrs)[1]
    last_pos_pcr_idx <- tail(which(maks_pos_pcrs), 1)
    # Check if surrounding negative PCRs are sufficiently close to positive PCRs
    # Remove all other negative PCRs later
    day_diff <- days[first_pos_pcr_idx] - days
    day_diff2 <- days - days[last_pos_pcr_idx]
    mask_time_pcrs <- (day_diff <= days_to_negative) & (day_diff2 <= days_to_negative)
    mask_all <- maks_pos_pcrs | mask_time_pcrs
  } else {
    mask_all <- maks_pos_pcrs
  }
  return(mask_all)
}


#' Trim and filter PCR test data.
#'
#' Removes negative PCR tests outside a window around the first/last positive 
#' test, shifts days to start at zero, and filters infections by minimum test 
#' count.
#'
#' @param dt A data.table. Input data with columns: ID, days_to_first_pos_pcr, 
#' days_to_last_pos_pcr, day, day_true.
#' @param min_n_tests An integer. Minimum number of tests required per infection.
#' @param days_to_negative An integer. Maximum days from first/last positive 
#' PCR to retain negative tests (default: 14).
#'
#' @return A data.table. Filtered and adjusted data.
trim_negative_PCRs_data <- function(dt, min_n_tests, days_to_negative = 14) {
  # We remove negative tests that happened > 14 days before or after the first
  # or last positive PCR in an infection
  dt <- dt %>%
    .[days_to_first_pos_pcr <= days_to_negative] %>%
    .[days_to_last_pos_pcr >= -days_to_negative] %>%
    .[, day := day - min(day), by = ID] %>%
    .[, day_true := day_true - min(day_true), by = ID] %>%
    .[, n_tests := .N, by = ID] %>%
    .[n_tests >= min_n_tests] # %>%
  
  return(dt)
}


#' Preprocess infection data for negative PCR addition.
#'
#' Identifies and propagates the earliest and latest negative PCR days for each
#' infection, based on viral load thresholds and relative timing to positive
#' tests.
#'
#' @param dt A data.table. Input data with columns: ID, day, log10_load,
#' days_to_first_pos_pcr, days_to_last_pos_pcr.
#'
#' @return A data.table. Input data with added columns: latest_neg_pcr_day,
#' earliest_neg_pcr_day.
add_neg_PCRs_infection_preprocess_dt <- function(dt) {
  dt <- dt %>%
    .[(log10_load <= VL_LOWER_LIMIT) & (days_to_last_pos_pcr < 0), 
      latest_neg_pcr_day := max(day), by = ID] %>%
    .[(log10_load <= VL_LOWER_LIMIT) & (days_to_first_pos_pcr > 0), 
      earliest_neg_pcr_day := min(day), by = ID] %>%
    # Propagate to all rows of the same infection.
    .[, latest_neg_pcr_day := propagate_to_group_rows(latest_neg_pcr_day), by = ID] %>%
    .[, earliest_neg_pcr_day := propagate_to_group_rows(earliest_neg_pcr_day), by = ID]
  
  return(dt)
}


#' Append synthetic negative PCRs to infection end.
#'
#' If the infection ends with a negative PCR, appends up to 5 additional
#' negative PCRs with specified intervals between tests.
#'
#' @param sub_dt A data.table. Infection data for a single infection.
#' @param n_neg_pcrs An integer. Number of negative PCRs to append.
#' @param n_days_between_pcrs An integer vector. Days between each appended
#' negative PCR.
#' @return A data.table. Input data with appended negative PCRs, if applicable.
add_new_neg_pcrs_infection_end <- function(sub_dt, 
                                           n_neg_pcrs, 
                                           n_days_between_pcrs) {
  if (any(!is.na(sub_dt$latest_neg_pcr_day))) {
    row_to_duplicate <- sub_dt[day == latest_neg_pcr_day, ][1]
    new_rows <- row_to_duplicate[rep(1, n_neg_pcrs)]
    new_rows[, day := day + n_days_between_pcrs]
    new_rows[, day_sampled := day_sampled + n_days_between_pcrs]
    return(rbind(sub_dt, new_rows))
  }
  return(sub_dt)
}


#' Prepend synthetic negative PCRs to infection start.
#'
#' If the infection begins with a negative PCR, prepends up to `n_neg_pcrs`
#' additional negative PCRs with specified intervals between tests.
#'
#' @param sub_dt A data.table. Infection data for a single infection.
#' @param n_neg_pcrs An integer. Number of negative PCRs to prepend.
#' @param n_days_between_pcrs An integer vector. Days between each prepended
#' negative PCR.
#'
#' @return A data.table. Input data with prepended negative PCRs, if applicable.
add_new_neg_pcrs_infection_start <- function(sub_dt, 
                                             n_neg_pcrs, 
                                             n_days_between_pcrs) {
  if (any(!is.na(sub_dt$earliest_neg_pcr_day))) {
    row_to_duplicate <- sub_dt[day == earliest_neg_pcr_day, ][1]
    new_rows <- row_to_duplicate[rep(1, n_neg_pcrs)]
    new_rows[, day := day - n_days_between_pcrs]
    new_rows[, day_sampled := day_sampled - n_days_between_pcrs]
    return(rbind(sub_dt, new_rows))
  }
  return(sub_dt)
}


#' Augment infection data with synthetic negative PCRs.
#'
#' Adds simulated negative PCR tests before the first and after the last
#' positive test, using specified counts and intervals. Reorders tests by day
#' and updates test counts.
#'
#' @param dt A data.table. Input data with columns: ID, day.
#' @param n_neg_pcrs_start An integer vector. Number of negative PCRs to add
#' before first positive, per infection.
#' @param n_neg_pcrs_end An integer vector. Number of negative PCRs to add
#' after last positive, per infection.
#' @param n_days_between_pcrs_start An integer vector. Days between negative
#' PCRs added at the start, per infection.
#' @param n_days_between_pcrs_end An integer vector. Days between negative PCRs
#' added at the end, per infection.
#'
#' @return A data.table. Data with added negative PCRs, ordered by day, and
#' updated test counts.
add_neg_PCRs_infection_process_dt <- function(dt, n_neg_pcrs_start,
                                              n_neg_pcrs_end,
                                              n_days_between_pcrs_start,
                                              n_days_between_pcrs_end) {
  # Add negative PCRs to the beginning and end of the infections.
  dt <- dt[, add_new_neg_pcrs_infection_start(.SD,
                                              n_neg_pcrs = n_neg_pcrs_start[.GRP],
                                              n_days_between_pcrs = n_days_between_pcrs_start[[.GRP]]
  ), by = ID]
  dt <- dt[, add_new_neg_pcrs_infection_end(.SD,
                                            n_neg_pcrs = n_neg_pcrs_end[.GRP],
                                            n_days_between_pcrs = n_days_between_pcrs_end[[.GRP]]
  ), by = ID]
  
  dt <- dt %>%
    # We actually don't want to recalculate the day variable here, else we will have
    # to allow for bigger shifts in the stan model which adds unnecessary
    # uncertainty
    .[, n_tests := .N, by = ID] %>%
    .[, .SD[order(day)], by = ID]
  
  return(dt)
}


#' Augment infections with synthetic negative PCRs.
#'
#' Processes infection data and adds 1–5 synthetic negative PCRs at the start 
#' and end of each infection, with random intervals (1–5 days) between tests.
#'
#' @param dt A data.table. Infection data with column: ID.
#' @return A data.table. Input data with added negative PCRs and updated 
#' metadata.
add_neg_PCRs_infection <- function(dt) {
  # Mark which infections qualify for PCR adding.
  dt <- add_neg_PCRs_infection_preprocess_dt(dt)
  
  # Randomly sample 5 negative PCRs with at most 5 days between two consecutive
  # PCRs to add to the beginning of an infection.
  n_infections <- uniqueN(dt$ID)
  n_neg_pcrs_start <- sample(1:5, n_infections, replace = TRUE)
  n_days_between_pcrs_start <- lapply(n_neg_pcrs_start, 
                                      function(n) cumsum(sample(1:5, n, replace = TRUE)))
  
  # Randomly sample 5 negative PCRs with at most 5 days between two consecutive
  # PCRs to add to the end of an infection.
  n_neg_pcrs_end <- sample(1:5, n_infections, replace = TRUE)
  n_days_between_pcrs_end <- lapply(n_neg_pcrs_end, 
                                    function(n) cumsum(sample(1:5, n, replace = TRUE)))
  
  dt <- add_neg_PCRs_infection_process_dt(dt,
                                          n_neg_pcrs_start = n_neg_pcrs_start,
                                          n_neg_pcrs_end = n_neg_pcrs_end,
                                          n_days_between_pcrs_start = n_days_between_pcrs_start,
                                          n_days_between_pcrs_end = n_days_between_pcrs_end
  )
  return(dt)
}


#' Process simulated viral load and associated data.
#'
#' @param dt_all A data.table.
#' @param min_n_tests An integer. Minimum number of tests required per infection.
#'
#' @return A processed data.table
process_sim_data <- function(dt_all, min_n_tests) {
  # Sort by day.
  dt <- dt_all %>%
    .[order(infection_hash, day, person_hash)] %>%
    # All viral load below 2 are considered negative.
    .[log10_load <= VL_LOWER_LIMIT_SIM_DATA, log10_load := 0] %>%
    # If multiple tests were performed on the same day, choose the first one.
    .[, .SD[1], by = .(ID, day)] %>%
    .[log10_load > VL_LOWER_LIMIT_SIM_DATA, first_pos_pcr_day := min(day, na.rm = TRUE), by = ID] %>%
    .[, is_first_pos_pcr := ifelse(day == first_pos_pcr_day, TRUE, FALSE)] %>%
    .[log10_load > VL_LOWER_LIMIT_SIM_DATA, last_pos_pcr_day := max(day, na.rm = TRUE), by = ID] %>%
    # Propagate the result to the negative PCRs in the infection.
    .[, first_pos_pcr_day := propagate_to_group_rows(first_pos_pcr_day), by = ID] %>%
    .[, last_pos_pcr_day := propagate_to_group_rows(last_pos_pcr_day), by = ID] %>%
    .[, days_to_first_pos_pcr := as.numeric(first_pos_pcr_day - day)] %>%
    .[, days_to_last_pos_pcr := as.numeric(last_pos_pcr_day - day)] %>%
    .[, last_test_negative := ifelse(tail(log10_load, 1) <= VL_LOWER_LIMIT_SIM_DATA, TRUE, FALSE), by = ID] %>%
    .[, first_test_negative := ifelse(head(log10_load, 1) <= VL_LOWER_LIMIT_SIM_DATA, TRUE, FALSE), by = ID] %>%
    .[, day := day - min(day), by = ID] %>%
    .[, n_tests := .N, by = ID] %>%
    .[n_tests >= min_n_tests] %>%
    .[, test_name := as.factor(test_name)] %>%
    .[, highest_vl_day := day[which.max(log10_load)], by = ID] %>%
    .[, day2 := day - highest_vl_day] %>%
    .[, starts_with_negative_pcr := log10_load[day == 0] == 0, by = ID] %>%
    # The following is just meant as a check to see if the drop in shift
    # parameter values around day 0 is related to likelihood instability
    # around the peak; here, we check if the drop happens around -5 which
    # would support this hypothesis.
    .[, day3 := day + 5] %>%
    # Label the day of the first positive PCR as day 0.
    .[, day4 := day - min(day[log10_load > VL_LOWER_LIMIT_SIM_DATA]), by = ID] %>%
    # Limit of detection in simulated data.
    .[, L := VL_LOWER_LIMIT_SIM_DATA] %>%
    .[, L_median := VL_LOWER_LIMIT_SIM_DATA] %>%
    .[, L_median5 := VL_LOWER_LIMIT_SIM_DATA] %>%
    .[, infection_duration := max(day[log10_load > 0]) - min(day[log10_load > 0]) + 1, by = ID] %>%
    .[infection_duration <= 30] %>%
    .[, material := "URT"] %>%
    .[, U := max(log10_load) + 0.1, by = test_name] %>%
    .[, multiple_infections := ifelse(length(unique(infection_hash)) > 1, TRUE, FALSE), by = person_hash] %>%
    .[, N_pos_tests := sum(log10_load > VL_LOWER_LIMIT_SIM_DATA), by = ID] %>%
    .[, max_load := max(log10_load), by = ID] %>%
    .[, max_load_day := max((log10_load == max_load) * day), by = ID] %>%
    .[, day4 := day - min(day[log10_load > VL_LOWER_LIMIT]), by = ID] %>%
    .[, shift_sd := 40 / N_pos_tests] %>%
    .[, shift_sd2 := 30 / N_pos_tests] %>%
    .[, shift_sd3 := 30 / N_pos_tests + (9 - max_load)]
  
  return(dt)
}


#' Sample from the discrete truncated normal distribution.
#'
#' @param n An integer. The number of samples.
#' @param min_val An integer. The lower bound of the distribution.
#' @param max_val An integer. The upper bound of the distribution.
#' @param mean The distribution's mean.
#' @param sd The distribution's standard deviation.
#'
#' @return An integer vector with values sampled from the truncated normal 
#' distribution.
discrete_truncnorm_sample <- function(n, min_val, max_val, mean, sd) {
  # Get integer support
  support <- min_val:max_val
  
  # Evaluate truncated normal density at integer points
  probs <- dtruncnorm(support,
                      a = min_val - 0.5, b = max_val + 0.5, mean = mean,
                      sd = sd
  )
  
  # Normalize to get a valid probability distribution
  probs <- probs / sum(probs)
  
  # Sample without replacement
  if (n > length(support)) stop("n is too large for the given range.")
  sample(support, size = n, replace = FALSE, prob = probs)
}


#' Sample viral loads from simulated viral load trajectories
#'
#' @param n number of infections to simulate
#' @param nPCR number of PCR targets (test_name)
#' @param sigma_up_slope standard deviation of up-slope of infection-level random
#' effects
#' @param trim_neg_pcrs Keep only the first preceding and trailing negative PCR
#' before or after a positive PCR in an infection if it is sufficiently close
#' to said positive PCR (at most days_to_negative days apart).
#' @param days_to_negative The maximum number of days between a leading or
#' trailing negative PCR to the closest positive PCR of the same infection.
#' @param add_neg_pcrs A logical value, indicating whether to artificially add
#' negative PCRs to the beginning and/or end of an infection (if that infection
#' starts and/or ends with a negative PCR).
#' @param only_pos_pcrs A logical value indicating whether to omit all negative
#' PCRs from an infection.
#' @param sigma_mu A numerical value specifying the location parameter of the
#' log10 viral load error distribution.
#' @param infection_wise_vl_error A logical value indicating whether each infection
#' should have its own error distribution.
#' @param three_slope_model A logical value indicating whether to simulate data
#' from a three slope instead of a two slope model.
#' @param error_dist A string specifying the log10 viral load error distribution
#' (normal, skew-normal or student-t).
#' @param infection_wise_alpha A logical value indicating whether each infection
#' should have its own smoothness parameter (determining how gradual the
#' infection transitions from proliferation to clearance phase).
#' @param mixed_courses A logical value indicating whether infections from two
#' populations, one with a smoother, one with a sharper transition at peak
#' viral load should be simulated.
#' @param min_max_load A numeric value indicating the minimum viral load that
#' at least one PCR in an infection should display.
#' @param min_duration The minimum number of days an infection should last.
#' @param single_leading_trailing_neg_pcr A logical value indicating whether each
#' infection should only have at most one single leading and/or trailing
#' negative PCR test.
#'
#' @return A named list with a data.table containing simulated viral load data
#' and a named list containing the parameter values used for simulating viral
#' load courses.
sim_vl_courses <- function(day_data,
                           n = NULL,
                           sigma_up_slope = 0.4,
                           param_sizes_all_zero = FALSE,
                           min_n_tests = 3,
                           max_n_tests = 18,
                           min_n_pos_tests = 2,
                           trim_neg_pcrs = TRUE,
                           days_to_negative = 14,
                           add_neg_pcrs = FALSE,
                           only_pos_pcrs = FALSE,
                           sigma_mu = NULL,
                           infection_wise_vl_error = TRUE,
                           three_slope_model = FALSE,
                           error_dist = "normal",
                           infection_wise_alpha = FALSE,
                           mixed_courses = FALSE,
                           min_max_load = 5,
                           min_duration = 0,
                           single_leading_trailing_neg_pcr = FALSE) {
  if (!is.null(error_dist)) {
    assert_that(error_dist %in% c("normal", "skew-normal", "student-t"))
  }
  if (!trim_neg_pcrs) {
    days_to_negative <- 100
  }
  
  assert_that((trim_neg_pcrs + add_neg_pcrs) < 2)
  assert_that((add_neg_pcrs + only_pos_pcrs) < 2)
  assert_that(max(day_data$infection_duration) <= 30)
  
  test_nameLevels <- filter_levels(day_data, "test_name", TEST_NAME_LEVELS)
  day_data[, test_name := factor(test_name, levels = test_nameLevels)]
  
  day_data[["variant_numeric"]] <- as.numeric(day_data$variant)
  day_data[["variant2_numeric"]] <- as.numeric(day_data$variant2)
  day_data[["test_name_numeric"]] <- as.numeric(day_data$test_name)
  
  
  if (!is.null(n)) {
    sampled_infections_unordered <- sample(unique(day_data$ID), n, 
                                         replace = FALSE)
    sampled_dt <- day_data[ID %in% sampled_infections]
    # Order sampled_infections_unordered by factor levels. Don't delete these 
    # lines!
    sampled_infections <- factor(sampled_infections_unordered, 
                                levels = levels(day_data$ID))
    sampled_infections <- sampled_infections[order(sampled_infections)]
  } else {
    sampled_infections <- day_data[day == 0]$ID
    sampled_dt <- day_data
    n <- uniqueN(day_data$ID)
  }
  
  if (mixed_courses) {
    # Mixture of infection courses with a very flat/smooth transition around
    # peak viral load, almost like a plateau phase
    ratio1 <- 0.15
    ratio2 <- 0.85
    assert_that(!infection_wise_alpha)
    sigma_intercept_param <- c(0.05, 0.0001)
    sigma_up_slope_param <- c(sigma_up_slope, 0.05)
    sigma_down_slope_param <- c(0.4, 0.3)
    log_up_slope_mu_param <- 0.9
    log_intercept_mu_param <- c(log(10), log(13.75))
    log_down_slope_mu_param <- c(-1.1, -0.3)
    alpha_mu_param <- c(8, 0.375)
    # Setting parameter sizes.
    # Note: we assume the same mean up-slope for all infections
    log_up_slope_mu <- rep(log_up_slope_mu_param, n)
    log_intercept_mu1 <- rep(log_intercept_mu_param[1], n * ratio1)
    log_down_slope_mu1 <- rep(log_down_slope_mu_param[1], n * ratio1)
    alpha_mu1 <- rep(alpha_mu_param[1], n * ratio1)
    
    # Make a subset of infections with a steeper down-slope
    # and a smaller alpha value
    log_intercept_mu2 <- rep(log_intercept_mu_param[2], n * ratio2)
    log_down_slope_mu2 <- rep(log_down_slope_mu_param[2], n * ratio2)
    alpha_mu2 <- rep(alpha_mu_param[2], n * ratio2)
    
    log_intercept_mu <- c(log_intercept_mu1, log_intercept_mu2)
    log_down_slope_mu <- c(log_down_slope_mu1, log_down_slope_mu2)
    alpha_mu <- c(alpha_mu1, alpha_mu2)
    smooth_transition <- c(rep(FALSE, n * ratio1), rep(TRUE, n * ratio2))
    
    sigma_intercept <- c(rep(sigma_intercept_param[1], n * ratio1), 
                        rep(sigma_intercept_param[2], n * ratio2))
    sigma_up_slope <- c(rep(sigma_up_slope_param[1], n * ratio1), 
                      rep(sigma_up_slope_param[2], n * ratio2))
    sigma_down_slope <- c(rep(sigma_down_slope_param[1], n * ratio1), 
                        rep(sigma_down_slope_param[2], n * ratio2))
  } else {
    ratio1 <- NULL
    ratio2 <- NULL
    sigma_intercept_param <- 0.1
    sigma_up_slope_param <- sigma_up_slope
    sigma_down_slope_param <- 0.4
    log_up_slope_mu_param <- 0.8
    log_intercept_mu_param <- 2.1
    log_down_slope_mu_param <- -1.3
    alpha_mu_param <- 8
    
    log_up_slope_mu <- rep(log_up_slope_mu_param, n)
    log_intercept_mu <- rep(log_intercept_mu_param, n)
    log_down_slope_mu <- rep(log_down_slope_mu_param, n)
    alpha_mu <- rep(alpha_mu_param, n)
    smooth_transition <- FALSE
    
    sigma_intercept <- sigma_intercept_param
    sigma_down_slope <- sigma_down_slope_param
  }
  
  if (infection_wise_alpha) {
    alpha <- abs(alpha_mu + rnorm(n, 0, 5))
  } else {
    alpha <- alpha_mu
  }
  
  # Generate intercept and slopes (categorical parameters should sum to 0).
  log_up_slope_variant <- list(
    "Wildtype" = -0.2, "Alpha" = 0, "Delta" = 0.1,
    "Omicron" = 0.2, "Unknown" = -0.1
  )
  log_intercept_variant <- list(
    "Wildtype" = -.05, "Alpha" = -.05, "Delta" = .05,
    "Omicron" = 0.1, "Unknown" = -.05
  )
  log_down_slope_variant <- list(
    "Wildtype" = -0.2, "Alpha" = -0.1, "Delta" = 0.1,
    "Omicron" = 0.2, "Unknown" = -0.1
  )
  
  log_up_slope_prior_infection <- c(0, -0.05)
  log_intercept_prior_infection <- c(0, -0.05)
  log_down_slope_prior_infection <- c(0, 0.2)
  
  log_up_slope_pams3 <- c(0, 0, 0)
  log_intercept_pams3 <- c(-0.05, 0.05, 0)
  log_down_slope_pams3 <- c(0.1, -0.1, 0)
  
  log_up_slope_infection <- rnorm(n = n, 0, sigma_up_slope)
  log_intercept_infection <- rnorm(n = n, 0, sigma_intercept)
  log_down_slope_infection <- rnorm(n = n, 0, sigma_down_slope)
  
  beta_sweight_mu <- 10
  
  log_up_slope_age <- c(0, 0, 0)
  log_intercept_age <- c(-0.05, 0.0, 0.05)
  log_down_slope_age <- c(0.075, 0.1, -0.075)
  
  centreNames <- c(
    "?", "AIR", "C19", "CP", "ED", "FM", "H", "ICU", "IDW", "L",
    "LW", "OD", "PRI", "RES", "SM", "WD"
  )
  log_intercept_centre1 <- list(
    "?" = 0, "AIR" = 0, "C19" = 0, "CP" = 0,
    "ED" = 0, "FM" = 0, "H" = 0, "ICU" = 0,
    "IDW" = 0, "L" = 0, "LW" = 0, "OD" = 0,
    "PRI" = 0, "RES" = 0, "SM" = 0, "WD" = 0
  )
  log_down_slope_ld_centre <- list(
    "?" = 0, "AIR" = 0, "C19" = 0, "CP" = 0, "ED" = 0,
    "FM" = 0, "H" = 0, "ICU" = 0, "IDW" = 0, "L" = 0,
    "LW" = 0, "OD" = 0, "PRI" = 0, "RES" = 0, "SM" = 0,
    "WD" = 0, "X" = 0
  )
  
  log_intercept_subject <- NULL
  log_up_slope_subject <- NULL
  log_down_slope_subject <- NULL
  
  intercept_pcr <- rnorm(n = uniqueN(day_data$test_name), 0, 0.5)
  intercept_centre <- setNames(rnorm(n = 16, 0, 0.005), centreNames)
  
  if (infection_wise_vl_error) {
    sigma <- numeric(n)
  } else {
    sigma <- 1
  }
  
  alpha_skew_normal <- NULL
  if (error_dist == "skew-normal") {
    alpha_skew_normal <- -2
  }
  
  log_intercept <- numeric(n)
  log_up_slope <- numeric(n)
  log_down_slope <- numeric(n)
  intercept <- numeric(n)
  up_slope <- numeric(n)
  down_slope <- numeric(n)
  days_sampled_list <- list()
  x_days <- list()
  x_days_true <- list()
  x_dates <- list()
  vl_list <- list()
  n_tests_vec <- numeric(n)
  
  if (only_pos_pcrs == TRUE) {
    min_n_pos_tests <- min_n_tests
  }
  
  
  for (i in 1:n) {
    dt_infection <- day_data[ID == sampled_infections[i]]
    row_infection <- dt_infection[day == 0]
    # The following is only used to assert the right order of the resulting 
    # data.table later on.
    n_tests_vec[i] <- row_infection$N_tests
    
    if (infection_wise_vl_error) {
      log_sigma <- rnorm(1, mean = 0, sd = 0.3)
      if (!is.null(sigma_mu)) {
        sigma[i] <- exp(sigma_mu + log_sigma)
      } else {
        sigma[i] <- exp(log_sigma)
      }
    }
    
    if (param_sizes_all_zero) {
      log_intercept[i] <- log_intercept_mu[i] + log_intercept_infection[i]
      log_up_slope[i] <- log_up_slope_mu[i] + log_up_slope_infection[i]
      log_down_slope[i] <- log_down_slope_mu[i] + log_down_slope_infection[i]
    } else {
      log_intercept[i] <- (log_intercept_mu[i]
                          + log_intercept_infection[i]
                          + as.numeric(log_intercept_variant[row_infection$variant_numeric])
                          + log_intercept_prior_infection[row_infection$prior_infection + 1]
                          + log_intercept_pams3[row_infection$PAMS3]
                          + log_intercept_age[row_infection$age_category3_code]
                          + as.numeric(log_intercept_centre1[dt_infection$test_centre_category[1]]))
      log_up_slope[i] <- (log_up_slope_mu[i]
                        + log_up_slope_infection[i]
                        + as.numeric(log_up_slope_variant[row_infection$variant_numeric])
                        + log_up_slope_age[row_infection$age_category3_code]
                        + log_up_slope_prior_infection[row_infection$prior_infection + 1]
                        + log_up_slope_pams3[row_infection$PAMS3])
      log_down_slope[i] <- (log_down_slope_mu[i]
                          + log_down_slope_infection[i]
                          + as.numeric(log_down_slope_variant[row_infection$variant_numeric])
                          + log_down_slope_pams3[row_infection$PAMS3]
                          + log_down_slope_prior_infection[row_infection$prior_infection + 1]
                          + log_down_slope_age[row_infection$age_category3_code]
                          + as.numeric(log_down_slope_ld_centre[row_infection$ld_centre]))
    }
    intercept[i] <- exp(log_intercept[i])
    up_slope[i] <- exp(log_up_slope[i])
    if (param_sizes_all_zero == TRUE) {
      down_slope[i] <- -exp(log_down_slope[i])
    } else {
      ld_centre <- row_infection$ld_centre
      down_slope[i] <- -exp(log_down_slope[i] + as.numeric(log_down_slope_ld_centre[ld_centre]))
    }
    sigma_current <- ifelse(infection_wise_vl_error, sigma[i], sigma)
    
    n_tests <- row_infection$N_tests
    pcrs_numeric <- dt_infection$test_name_numeric
    centres <- dt_infection$test_centre_category
    vls <- c(-1, -1, -1)
    days <- c(0, 0, 0)
    infection_duration <- 0
    vl_last_pos_pcr <- sum(dt_infection$log10_load > VL_LOWER_LIMIT_SIM_DATA)
    vl_last_pos_pcr <- 12
    
    assert_that(n_tests >= min_n_tests)
    n_pcrs_valid <- 0
    
    # We want at least minNPosPcrs positive PCRs 
    # (positive defined as vl > VL_LOWER_LIMIT_SIM_DATA).
    while ((sum(vls > VL_LOWER_LIMIT_SIM_DATA) < min_n_pos_tests) |
           (n_pcrs_valid != n_tests) |
           (max(vls) < min_max_load) |
           (infection_duration > 30) |
           (infection_duration < min_duration) |
           (infection_duration >= 20) & (vl_last_pos_pcr >= 8)) {
      days_sampled <- sort(discrete_truncnorm_sample(
        n = n_tests, min_val = -17,
        max_val = 43, mean = 5, sd = 18
      ))
      
      # Add some jittering but not so much that the day of sampling changes
      days_sampled <- days_sampled + runif(n = n_tests, min = -0.49, max = 0.49)
      # First day should be 0.
      days <- round(days_sampled) - round(min(days_sampled))
      assert_that(uniqueN(days) == length(days))
      days_sampled <- unlist(days_sampled)
      
      vls <- compute_viral_loads(days_sampled,
                                 beta_sweight_mu = beta_sweight_mu,
                                 intercept = intercept[i],
                                 up_slope = up_slope[i],
                                 down_slope = down_slope[i],
                                 intercept_pcr = intercept_pcr,
                                 pcrs_numeric, intercept_centre = intercept_centre,
                                 centres = centres,
                                 sigma = sigma_current,
                                 param_sizes_all_zero = param_sizes_all_zero,
                                 error_dist = error_dist,
                                 alpha = alpha[i],
                                 three_slope_model = three_slope_model,
                                 alpha_skew_normal = alpha_skew_normal
      )
      
      if (trim_neg_pcrs) {
        maks_pcrs_valid <- find_valid_PCRs(
          vls = vls,
          days = days,
          days_to_negative = days_to_negative
        )
        n_pcrs_valid <- sum(maks_pcrs_valid)
      } else {
        n_pcrs_valid <- length(vls)
      }
      # Set a random start date (e.g., within the last year)
      start_date <- as.Date("2024-01-01") + sample(0:365, 1)
      dates <- start_date + days
      days_true <- days_sampled - min(days_sampled)
      
      days_pos <- days[vls > VL_LOWER_LIMIT_SIM_DATA]
      if (length(days_pos)) {
        firstPosDay <- days_pos[1]
        lastPosDay <- tail(days_pos, 1)
        infection_duration <- lastPosDay - firstPosDay + 1
        vl_last_pos_pcr <- tail(vls[vls > VL_LOWER_LIMIT_SIM_DATA], 1)
      } else {
        infection_duration <- 0
        vl_last_pos_pcr <- 12
      }
    }
    x_days[[i]] <- days
    x_days_true[[i]] <- days_true
    x_dates[[i]] <- dates
    vl_list[[i]] <- vls
    days_sampled_list[[i]] <- days_sampled
  }
  
  n_testsPerInfection <- sampled_dt[day == 0]$N_tests
  # The following asserts that we used the right order for creating infections.
  assert_that(all(n_tests_vec == n_testsPerInfection))
  n_testsTotal <- nrow(sampled_dt)
  # Combine all data.
  dt_all <- data.table(
    log10_load_true = unlist(vl_list),
    log10_load = unlist(vl_list),
    log10_load_orig = unlist(vl_list),
    day = unlist(x_days),
    day_true = unlist(x_days_true),
    pcr_date = as.Date(unlist(x_dates)),
    day_sampled = unlist(days_sampled_list),
    test_name = sampled_dt$test_name,
    test_centre = sampled_dt$test_centre,
    test_centre_category = sampled_dt$test_centre_category,
    ld_centre = sampled_dt$ld_centre,
    infection_hash = sampled_dt$infection_hash,
    ID = as.integer(factor(sampled_dt$infection_hash)),
    person_hash = sampled_dt$person_hash,
    prior_infection = sampled_dt$prior_infection,
    hospitalized = sampled_dt$hospitalized,
    PAMS1 = sampled_dt$PAMS1,
    PAMS3 = sampled_dt$PAMS3,
    variant = sampled_dt$variant,
    variant2 = sampled_dt$variant2,
    age = sampled_dt$age,
    age_category3 = sampled_dt$age_category3,
    age_category3_code = sampled_dt$age_category3_code,
    gender = sampled_dt$gender,
    # We only have 1 prior infection at most.
    n_prior_infections = as.integer(sampled_dt$prior_infection == 1),
    recent_prior_infection = rep(NA, n_testsTotal),
    onset_day = rep(NA, n_testsTotal),
    keep_onset = rep(FALSE, n_testsTotal),
    onset_json = rep(NA, n_testsTotal),
    # Expand infection level params.
    intercept = rep(intercept, times = n_testsPerInfection),
    up_slope = rep(up_slope, times = n_testsPerInfection),
    down_slope = rep(down_slope, times = n_testsPerInfection),
    sigma = rep(sigma, times = n_testsPerInfection),
    log_intercept_infection = rep(log_intercept_infection, 
                                times = n_testsPerInfection),
    log_up_slope_infection = rep(log_up_slope_infection, times = n_testsPerInfection),
    log_down_slope_infection = rep(log_down_slope_infection, 
                                times = n_testsPerInfection),
    alpha = rep(alpha, times = n_testsPerInfection),
    smooth_transition = rep(smooth_transition, times = n_testsTotal)
  )
  
  if (only_pos_pcrs == TRUE) {
    dt_all <- dt_all[log10_load > VL_LOWER_LIMIT_SIM_DATA]
  }
  dt <- process_sim_data(dt_all, min_n_tests = min_n_tests)
  
  # Remove negative PCR tests that are too far away in time from the first or
  # last positive PCR test and only use at most one first and one last negative.
  # We also (maybe counterintuitively) do this before adding negative PCRs so
  # that the negative PCRs don't extend too far out.
  if ((trim_neg_pcrs == TRUE) | (add_neg_pcrs == TRUE)) {
    dt <- trim_negative_PCRs_data(dt,
                                  min_n_tests = min_n_tests,
                                  days_to_negative = days_to_negative
    )
  }
  if (single_leading_trailing_neg_pcr) {
    dt[, day_to_first_pos_pcr := day - min(day[log10_load > VL_LOWER_LIMIT]), by = ID]
    dt[, day_to_last_pos_pcr := day - max(day[log10_load > VL_LOWER_LIMIT]), by = ID]
    dt <- filter_last_leading_first_trailing_neg_pcrs(dt)
    dt[, day := day - min(day), by = ID]
  }
  dt <- dt %>% .[, ID := factor(ID)]
  
  if (param_sizes_all_zero == TRUE) {
    params_list <- list(
      sigma_mu = sigma_mu,
      log_up_slope_mu = log_up_slope_mu,
      log_intercept_mu = log_intercept_mu,
      log_down_slope_mu = log_down_slope_mu,
      sigma_up_slope = sigma_up_slope,
      sigma_intercept = sigma_intercept,
      sigma_down_slope = sigma_down_slope,
      min_n_tests = min_n_tests,
      max_n_tests = max_n_tests,
      min_n_pos_tests = min_n_pos_tests,
      only_pos_pcrs = only_pos_pcrs,
      add_neg_pcrs = add_neg_pcrs,
      error_dist = error_dist,
      three_slope_model = three_slope_model,
      alpha_skew_normal = alpha_skew_normal,
      ratio1 = ratio1,
      ratio2 = ratio2,
      param_sizes_all_zero = param_sizes_all_zero,
      mixed_courses = mixed_courses,
      alpha_mu = alpha_mu
    )
  } else {
    params_list <- list(
      sigma_mu = sigma_mu,
      log_up_slope_mu = log_up_slope_mu,
      log_intercept_mu = log_intercept_mu,
      log_down_slope_mu = log_down_slope_mu,
      log_up_slope_variant = log_up_slope_variant,
      log_intercept_variant = log_intercept_variant,
      log_down_slope_variant = log_down_slope_variant,
      log_up_slope_prior_infection = log_up_slope_prior_infection,
      log_intercept_prior_infection = log_intercept_prior_infection,
      log_down_slope_prior_infection = log_down_slope_prior_infection,
      log_up_slope_pams3 = log_up_slope_pams3,
      log_intercept_pams3 = log_intercept_pams3,
      log_down_slope_pams3 = log_down_slope_pams3,
      sigma_up_slope = sigma_up_slope,
      sigma_intercept = sigma_intercept,
      sigma_down_slope = sigma_down_slope,
      beta_sweight_mu = beta_sweight_mu,
      log_up_slope_age = log_up_slope_age,
      log_intercept_age = log_intercept_age,
      log_down_slope_age = log_down_slope_age,
      log_intercept_centre1 = log_intercept_centre1,
      log_down_slope_ld_centre = log_down_slope_ld_centre,
      log_intercept_subject = log_intercept_subject,
      log_up_slope_subject = log_up_slope_subject,
      log_down_slope_subject = log_down_slope_subject,
      intercept_pcr = intercept_pcr,
      intercept_centre = intercept_centre,
      min_n_tests = min_n_tests,
      max_n_tests = max_n_tests,
      min_n_pos_tests = min_n_pos_tests,
      only_pos_pcrs = only_pos_pcrs,
      add_neg_pcrs = add_neg_pcrs,
      error_dist = error_dist,
      three_slope_model = three_slope_model,
      ratio1 = ratio1,
      ratio2 = ratio2,
      param_sizes_all_zero = param_sizes_all_zero,
      mixed_courses = mixed_courses,
      alpha_mu = alpha_mu
    )
  }
  return(list(dt = dt, params = params_list))
}


#### Thesis specific utilities ####

#' Process a data table containing posterior estimates or Rhat values on model 
#' parameters of interest.
#'
#' @param dt A data.table with posterior estimates or Rhat values.
#' @param model_nos A vector of model numbers to include (refers to stan file 
#' numbers).
#' @param includes_sim_data A logical value indicating whether dt includes data 
#' from simulation runs.
#' @param post_pred_groups A logical value indicating whether dt contains 
#' posterior predictions for infection groups.
#' @param rhat A logical value indicating whether dt contains Rhat values (as 
#' opposed to posterior estimates).
#' @param target_accept_value A numeric value indicating the selected target 
#' acceptance rate for MCMC.
#' @param sampling_iterations A numeric value or NULL specifying the number of 
#' sampling iterations run during model fitting. If provided, this value is 
#' used as a filter criterion; if NULL, no filtering is applied.
#'
#' @return A processed data.table.
process_dt <- function(dt,
                       model_nos = NULL,
                       includes_sim_data = TRUE,
                       post_pred_groups = FALSE,
                       rhat = FALSE,
                       target_accept_value = 0.8,
                       sampling_iterations = NULL) {
  assert_that((post_pred_groups + rhat) < 2)
  # Remove duplicates
  if ((post_pred_groups == FALSE) && (rhat == FALSE)) {
    if (includes_sim_data) {
      dt <- unique(dt, by = c("model", "params_true"))
    } else {
      dt <- unique(dt, by = "model")
    }
  }
  if (rhat) {
    dt <- unique(dt, by = "model")
  }
  dt[, target_accept := {
    # Extract the second digit after "adapt_delta_X_", default to 8 if no match
    extracted <- str_extract(model, "(?<=adapt_delta_\\d_)\\d")
    ifelse(is.na(extracted), 0.8, as.numeric(extracted) / 10)
  }]
  dt[target_accept == target_accept_value]
  dt[, simulated := ifelse(grepl("_simvl_", model), TRUE, FALSE)]
  dt <- copy(dt[simulated == FALSE])
  dt[, min_n_tests := as.factor(as.numeric(str_extract(model, "(?<=sel)\\d+")))]
  dt[, min_n_pos_tests := as.numeric(str_extract(model, "(?<=tests)\\d+"))]
  # Extract the number of sampling iterations
  dt[, iter_sampling := as.numeric(str_extract(model, "(?<=n_iter)\\d+"))]
  dt[, iteration := ifelse(grepl("iter2$", model), 2, 1)]
  dt[is.na(min_n_pos_tests), min_n_pos_tests := 2]
  dt[, min_n_pos_tests := as.factor(min_n_pos_tests)]
  
  if (!is.null(model_nos)) {
    assert_that(is.vector(model_nos))
    dt <- dt[model_no_stan %in% model_nos]
  }
  if (!is.null(sampling_iterations)) {
    assert_that(sampling_iterations %in% c(1000, 2000))
    dt <- dt[iter_sampling == sampling_iterations]
  }
  return(dt)
}


#' Process a data table containing posterior estimates or Rhat values on model 
#' parameters of interest, return only estimates/Rhat values from runs with 
#' simulated data.
#'
#' @param dt A data.table with posterior estimates or Rhat values.
#' @param model_no A model number to filter by (refers to the number of the 
#' stan file) (or NULL).
#' @param param_sizes_all_zero A logical value indicating whether to extract 
#' only estimates from simulation runs where all fixed effects 
#' (infection-group-specific parameters) were set to 0 or only those where
#'  these values were different from 0.
#' @param exclude_true_values A logical value indicating whether to exclude 
#' rows that contain true parameter values of used in simulation runs.
#' @param rhat A logical value indicating whether dt contains Rhat values (as 
#' opposed to posterior estimates).
#'
#' @return A processed data table containing only estimates or Rhat values from 
#' runs with simulated data.
process_sim_dt <- function(dt,
                           model_no = NULL,
                           param_sizes_all_zero = TRUE,
                           exclude_true_values = FALSE,
                           rhat = FALSE) {
  if (!is.null(model_no)) {
    assert_that(model_no %in% c(1, 2, 3))
  }
  # We have to do this because model_no is a column name in the dt
  model_no_current <- model_no
  param_sizes_all_zero_current <- param_sizes_all_zero
  # One row per model configuration (somehow the file ended up having 
  # duplicates)
  if (rhat) {
    dt <- unique(dt, by = c("model"))
  } else {
    dt <- unique(dt, by = c("model", "params_true"))
  }
  dt[, iter_sampling := 1000]
  dt[, simulated := ifelse(grepl("_simvl_", model), TRUE, FALSE)]
  dt <- copy(dt[simulated == TRUE])
  dt[, min_n_tests := as.factor(as.numeric(str_extract(model, "(?<=sel)\\d+")))]
  dt[, min_n_pos_tests := as.numeric(str_extract(model, "(?<=tests)\\d+"))]
  dt[, iteration := ifelse(grepl("iter2", model), 2, 1)]
  dt[, param_sizes_all_zero := ifelse(grepl("param_sizes_all_zero", model), 
                                      TRUE, FALSE)]
  dt[is.na(min_n_pos_tests), min_n_pos_tests := 2]
  dt[, min_n_pos_tests := as.factor(min_n_pos_tests)]
  # Find and remove columns with "_avg" in their names
  avg_columns <- names(dt)[grep("_avg", names(dt), fixed = TRUE)]
  dt <- dt[, !avg_columns, with = FALSE]
  if (exclude_true_values) {
    dt <- dt[params_true == 0]
  }
  if (!is.null(model_no_current)) {
    dt <- dt[model_no_stan == model_no_current]
  }
  if (!is.null(param_sizes_all_zero_current)) {
    dt <- dt[param_sizes_all_zero == param_sizes_all_zero_current]
  }
  return(dt)
}

#' Add a column to a data.table with posterior estimates from runs with 
#' simulated data, indicating how estimated and true parameter values compare.
#'
#' @param dt A data.table containing posterior estimates from runs with 
#' simulated data/
#' @param col_names A vector of dt column names (corresponding to model 
#' parameter names) for which to compare estimated and true values.
#'
#' @return A data.table with a new column added indicating how estimated and 
#' true parameter value compare.
create_hdi_combinations <- function(dt, col_names) {
  # Create new column with 4 values:
  # 1) contains true value, excludes 0
  # 2) contains true value, includes 0
  # 3) does not contain true value, includes 0
  # 4) does not contain true value, excludes 0
  process_pair <- function(dt, col_name) {
    true_col <- paste0(col_name, "_hdi_contains_true_value")
    zero_col <- paste0(col_name, "_hdi_contains_zero")
    
    # Create a new column with the four possible combinations
    dt[[paste0(col_name, "_hdi_combination")]] <- with(dt, {
      case_when(
        get(true_col) == TRUE & get(zero_col) == FALSE ~ 1,
        get(true_col) == TRUE & get(zero_col) == TRUE ~ 2,
        get(true_col) == FALSE & get(zero_col) == TRUE ~ 3,
        get(true_col) == FALSE & get(zero_col) == FALSE ~ 4
      )
    })
    return(dt)
  }
  for (col_name in col_names) {
    dt <- process_pair(dt, col_name)
  }
  return(dt)
}

#' Return data tables containing differences between estimated and true 
#' parameter values from runs with simulated data.
#'
#' @param model_no A model number to filter by (refers to the number of the 
#' stan file).
#' @param param_sizes_all_zero A logical value indicating whether to extract 
#' only estimates from simulation runs where all fixed effects 
#' (infection-group-specific parameters) were set to 0 or only those where
#' these values were different from 0.
#'
#' @return A list of data tables, containing absolute and relative differences 
#' between estimated and true parameter values.
get_sim_diff_dts <- function(model_no, param_sizes_all_zero) {
  dt <- process_sim_dt(fread(PARAMS_DIFF_FILE, sep = "\t"),
                       model_no = model_no,
                       param_sizes_all_zero = param_sizes_all_zero
  )
  dtr <- process_sim_dt(fread(PARAMS_REL_DIFF_FILE, sep = "\t"),
                        model_no = model_no,
                        param_sizes_all_zero = param_sizes_all_zero
  )
  dthdi <- process_hdi_cols(process_sim_dt(fread(PARAMS_DIFF_HDI_FILE, sep = "\t"),
                                           model_no = model_no,
                                           param_sizes_all_zero = param_sizes_all_zero
  ))
  dthdi50 <- process_hdi_cols(process_sim_dt(fread(PARAMS_DIFF_HDI50_FILE, sep = "\t"),
                                             model_no = model_no,
                                             param_sizes_all_zero = param_sizes_all_zero
  ))
  dtrhdi <- process_hdi_cols(process_sim_dt(fread(PARAMS_REL_DIFF_HDI_FILE, sep = "\t"),
                                            model_no = model_no,
                                            param_sizes_all_zero = param_sizes_all_zero
  ))
  dtrhdi50 <- process_hdi_cols(process_sim_dt(fread(PARAMS_REL_DIFF_HDI50_FILE, sep = "\t"),
                                              model_no = model_no,
                                              param_sizes_all_zero = param_sizes_all_zero
  ))
  dthdi_orig <- process_hdi_cols(process_sim_dt(fread(PARAMS_HDI_FILE, sep = "\t")[params_true == FALSE],
                                                model_no = model_no,
                                                param_sizes_all_zero = param_sizes_all_zero,
                                                exclude_true_values = TRUE
  ), diffs = FALSE)
  
  contains_zero_value_cols <- c("model", names(dthdi_orig)[grepl("_hdi_contains_zero", names(dthdi_orig))])
  
  dt <- merge(dt, dthdi, by = "model", suffix = c("", "_hdi"))
  dt <- merge(dt, dthdi50, by = "model", suffix = c("", "_hdi50"))
  dt <- merge(dt, dthdi_orig[, ..contains_zero_value_cols], by = "model", suffix = c("", "_hdi_orig"))
  group_cols <- unlist(unname(return_group_cols(log_prefix = "log_")$group_cols_mapping))
  dt <- create_hdi_combinations(dt, group_cols)
  # Find column names containing "_combination"
  combination_cols <- c("model", names(dt)[grepl("_combination", names(dt))])
  
  dtr <- merge(dtr, dtrhdi, by = "model", suffixes = c("", "_hdi"))
  dtr <- merge(dtr, dtrhdi50, by = "model", suffixes = c("", "_hdi50"))
  dtr <- merge(dtr, dt[, ..combination_cols], by = "model", 
               suffixes = c("", "_hdi2col"))
  # Rename combination colums (remove "log" prefix)
  setnames(
    dtr,
    old = grep("^log_.*_hdi_combination$", names(dtr), value = TRUE),
    new = sub("^log_", "", grep("^log_.*_hdi_combination$", names(dtr), 
                                value = TRUE))
  )
  
  return(list(dt = dt, dtr = dtr))
}

#' Return data tables containing estimated and true parameter values from runs 
#' with simulated data.
#'
#' @param model_no A model number to filter by (refers to the number of the 
#' stan file).
#' @param param_sizes_all_zero A logical value indicating whether to extract 
#' only estimates from simulation runs where all fixed effects 
#' (infection-group-specific parameters) were set to 0 or only those where
#' these values were different from 0.
#'
#' @return A list of two data tables, one containing estimated and true 
#' parameter values, the other containing exponentiated true and estimated 
#' parameter values (for those parameters where exponentiation is appropriate,
#' such as infection group specific fixed effects).
get_sim_dts <- function(model_no, param_sizes_all_zero) {
  dt <- process_sim_dt(fread(PARAMS_FILE, sep = "\t"),
                       model_no = model_no,
                       param_sizes_all_zero = param_sizes_all_zero,
                       exclude_true_values = TRUE
  )
  dte <- process_sim_dt(fread(PARAMS_EXP_FILE, sep = "\t"),
                        model_no = model_no,
                        param_sizes_all_zero = param_sizes_all_zero,
                        exclude_true_values = TRUE
  )
  dthdi_diff <- process_hdi_cols(process_sim_dt(fread(PARAMS_DIFF_HDI_FILE, sep = "\t"),
                                                model_no = model_no,
                                                param_sizes_all_zero = param_sizes_all_zero,
                                                exclude_true_values = FALSE
  ))
  dthdi50_diff <- process_hdi_cols(process_sim_dt(fread(PARAMS_DIFF_HDI50_FILE, sep = "\t"),
                                                  model_no = model_no,
                                                  param_sizes_all_zero = param_sizes_all_zero,
                                                  exclude_true_values = FALSE
  ))
  dthdi <- process_hdi_cols(process_sim_dt(fread(PARAMS_HDI_FILE, sep = "\t")[params_true == FALSE],
                                           model_no = model_no,
                                           param_sizes_all_zero = param_sizes_all_zero,
                                           exclude_true_values = TRUE
  ), diffs = FALSE)
  dthdi50 <- process_hdi_cols(process_sim_dt(fread(PARAMS_HDI50_FILE, sep = "\t")[params_true == FALSE],
                                             model_no = model_no,
                                             param_sizes_all_zero = param_sizes_all_zero,
                                             exclude_true_values = TRUE
  ), diffs = FALSE)
  dtehdi <- parse_hdi_cols(process_sim_dt(fread(PARAMS_EXP_HDI_FILE, sep = "\t")[params_true == FALSE],
                                          model_no = model_no,
                                          param_sizes_all_zero = param_sizes_all_zero,
                                          exclude_true_values = TRUE
  ))
  dtehdi50 <- parse_hdi_cols(process_sim_dt(fread(PARAMS_EXP_HDI50_FILE, sep = "\t")[params_true == FALSE],
                                            model_no = model_no,
                                            param_sizes_all_zero = param_sizes_all_zero,
                                            exclude_true_values = TRUE
  ))
  
  contains_true_value_cols <- c("model", names(dthdi_diff)[grepl("_hdi_contains_true_value", names(dthdi_diff))])
  
  dt <- merge(dt, dthdi, by = "model", suffix = c("", "_hdi"))
  dt <- merge(dt, dthdi50, by = "model", suffix = c("", "_hdi50"))
  dt <- merge(dt, dthdi_diff[, ..contains_true_value_cols], by = "model", 
              suffix = c("", "_hdi_diff"))
  group_cols <- unlist(unname(return_group_cols(log_prefix = "log_")$group_cols_mapping))
  dt <- create_hdi_combinations(dt, group_cols)
  # Find column names containing "_combination"
  combination_cols <- c("model", names(dt)[grepl("_combination", names(dt))])
  new_cols <- c(combination_cols, contains_true_value_cols[2:length(contains_true_value_cols)])
  
  dte <- merge(dte, dtehdi, by = "model", suffixes = c("", "_hdi"))
  dte <- merge(dte, dtehdi50, by = "model", suffixes = c("", "_hdi50"))
  dte <- merge(dte, dt[, ..new_cols], by = "model", suffixes = c("", "_hdi2col"))
  # Rename combination colums (remove "log" prefix)
  setnames(
    dte,
    old = grep("^log_.*_hdi_combination$", names(dte), value = TRUE),
    new = sub("^log_", "", grep("^log_.*_hdi_combination$", names(dte), 
                                value = TRUE))
  )
  setnames(
    dte,
    old = grep("^log_.*_contains_true_value$", names(dte), value = TRUE),
    new = sub("^log_", "", grep("^log_.*_hdi_contains_true_value$", names(dte), 
                                value = TRUE))
  )
  
  return(list(dt = dt, dte = dte))
}

#' Return data tables containing conditional posterior predictions, computed 
#' from population-level parameter estimates.
#'
#' @param ref_variant A string indicating the variant used for prediction.
#' @param ref_hosp_status A string indicating the disease severity status (PAMS 
#' or hospitalized) used for prediction.
#' @param model_nos A vector of model numbers to include (refers to stan file 
#' numbers).
#' @param sampling_iterations A numeric value or NULL specifying the number of 
#' sampling iterations run during model fitting. If provided, this value is 
#' used as a filter criterion; if NULL, no filtering is applied.
#'
#' @return A data table with posterior predictions.
get_post_pred_group_dt <- function(ref_variant = "omicron",
                                   ref_hosp_status = "hospitalized",
                                   model_nos = c(1, 2),
                                   sampling_iterations = NULL) {
  assert_that(ref_hosp_status %in% c("hospitalized", "pams"))
  assert_that(ref_variant %in% c("wildtype", "omicron"))
  if ((ref_variant == "omicron") && (ref_hosp_status == "hospitalized")) {
    dt_file <- PARAMS_GROUP_POST_PRED_OMI_HOSP_FILE
  } else if ((ref_variant == "wildtype") && (ref_hosp_status == "hospitalized")) {
    dt_file <- PARAMS_GROUP_POST_PRED_WT_HOSP_FILE
  } else if ((ref_variant == "omicron") && (ref_hosp_status == "pams")) {
    dt_file <- PARAMS_GROUP_POST_PRED_OMI_PAMS_FILE
  } else {
    dt_file <- PARAMS_GROUP_POST_PRED_WT_PAMS_FILE
  }
  dt <- process_dt(fread(dt_file, sep = "\t"),
                   includes_sim_data = FALSE,
                   post_pred_groups = TRUE,
                   model_nos = model_nos,
                   sampling_iterations = sampling_iterations
  )
  return(dt)
}

#' Return data tables containing posterior estimates of model parameter values.
#'
#' @param model_no A model number to filter by (refers to the number of the 
#' stan file).
#' @param sampling_iterations A numeric value or NULL specifying the number of 
#' sampling iterations run during model fitting. If provided, this value is used 
#' as a filter criterion; if NULL, no filtering is applied.
#'
#' @return A list of two data tables, one containing estimated parameter values, 
#' the other containing exponentiated estimated parameter values (for those 
#' parameters where exponentiation is appropriate, such as infection group 
#' specific fixed effects).
get_dts <- function(model_nos, sampling_iterations = NULL) {
  assert_that(is.vector(model_nos))
  dt <- process_dt(fread(PARAMS_FILE, sep = "\t"),
                   model_nos = model_nos,
                   sampling_iterations = sampling_iterations
  )
  dte <- process_dt(fread(PARAMS_EXP_FILE, sep = "\t"),
                    model_nos = model_nos,
                    sampling_iterations = sampling_iterations
  )
  dthdi <- parse_hdi_cols(process_dt(fread(PARAMS_HDI_FILE, sep = "\t"),
                                     model_nos = model_nos,
                                     sampling_iterations = sampling_iterations
  ))
  dtehdi <- parse_hdi_cols(process_dt(fread(PARAMS_EXP_HDI_FILE, sep = "\t"),
                                      model_nos = model_nos,
                                      sampling_iterations = sampling_iterations
  ))
  dthdi50 <- parse_hdi_cols(process_dt(fread(PARAMS_HDI50_FILE, sep = "\t"),
                                       model_nos = model_nos,
                                       sampling_iterations = sampling_iterations
  ))
  dtehdi50 <- parse_hdi_cols(process_dt(fread(PARAMS_EXP_HDI50_FILE, sep = "\t"),
                                        model_nos = model_nos,
                                        sampling_iterations = sampling_iterations
  ))
  
  dt <- merge(dt, dthdi, by = "model", suffix = c("", "_hdi"))
  dte <- merge(dte, dtehdi, by = "model", suffixes = c("", "_hdi"))
  dt <- merge(dt, dthdi50, by = "model", suffix = c("", "_hdi50"))
  dte <- merge(dte, dtehdi50, by = "model", suffixes = c("", "_hdi50"))
  
  return(list(dt = dt, dte = dte))
}

#' Return data tables with conditional posterior predictions computed from 
#' infection-level parameter estimates.
#'
#' @param model_nos A vector of model numbers to include (refers to stan file 
#' numbers).
#' @param sampling_iterations A numeric value or NULL specifying the number of 
#' sampling iterations run during model fitting. If provided, this value is 
#' used as a filter criterion; if NULL, no filtering is applied.
#'
#' @return A data table with posterior predictions.
get_group_avg_dts <- function(model_nos, sampling_iterations = NULL) {
  assert_that(is.vector(model_nos))
  dt <- process_dt(fread(PARAMS_GROUP_AVG_POST_PRED_FILE, sep = "\t"),
                   model_nos = model_nos,
                   sampling_iterations = sampling_iterations,
                   includes_sim_data = FALSE
  )
  dthdi <- parse_hdi_cols(process_dt(fread(PARAMS_GROUP_AVG_POST_PRED_HDI_FILE, 
                                           sep = "\t"),
                                     model_nos = model_nos,
                                     sampling_iterations = sampling_iterations,
                                     includes_sim_data = FALSE
  ))
  dthdi50 <- parse_hdi_cols(process_dt(fread(PARAMS_GROUP_AVG_POST_PRED_HDI50_FILE, 
                                             sep = "\t"),
                                       model_nos = model_nos,
                                       sampling_iterations = sampling_iterations,
                                       includes_sim_data = FALSE
  ))
  
  dt <- merge(dt, dthdi, by = "model", suffix = c("", "_hdi"))
  dt <- merge(dt, dthdi50, by = "model", suffix = c("", "_hdi50"))
  
  return(dt)
}


#' Convert hdi values stored as strings into numeric vectors of length 2.
#'
#' @param dt A data.table containing columns with HPDI values as strings of the 
#' form "c(x, y)".
#' @param col The name of the column containing HPDI values.
#'
#' @return The data.table "dt" with HPDI values in column "col" converted.
parse_hdi_col <- function(dt, col) {
  parsed_col <- lapply(dt[[col]], function(s) {
    # Extract numbers including scientific notation
    nums <- regmatches(s, gregexpr("-?\\d*\\.?\\d+(e[+-]?\\d+)?", s, 
                                   ignore.case = TRUE))[[1]]
    as.numeric(nums)
  })
  dt[, (col) := parsed_col]
}


#' Convert hdi values stored as strings into numeric vectors of length 2.
#'
#' @param dt A data.table containing columns with HPDI values as strings of the 
#' form "c(x, y)".
#'
#' @return The data.table "dt" with HPDI values in all columns (characterized 
#' by the format "c(x, y)") converted to numeric vectors of length 2.
parse_hdi_cols <- function(dt) {
  for (col in colnames(dt)) {
    if (is_vector_column(dt[[col]])) {
      # Step 1: parse the column safely
      parse_hdi_col(dt, col)
    }
  }
  return(dt)
}

# Check if a column contains strings generated from 2-element vectors 
# ("c(x, y)").
is_vector_column <- function(col) {
  return(grepl("c\\(", col[1]))
}


#' Process HPDI columns by adding new logical columns indicating if the HPDI 
#' contains the true value (in the case of simulated data runs) or 0.
#'
#' @param dt A data.table containing columns with HPDI values as strings of the 
#' form "c(x, y)".
#' @param diffs A logical value indicating whether "dt" contains HPDIs of 
#' differences between estimated and true parameters (in the case of simulated 
#' data runs) or of the estimated parameter values.
#'
#' @return The data.table "dt" with a new column for each HPDI column, 
#' indicating whether the HPDI contains the true value (if diffs==TRUE) or 0.
process_hdi_cols <- function(dt, diffs = TRUE) {
  # diffs refers to whether we are processing data.tables containing the
  # differences between estimated and true parameter or the actual estimated
  # parameters
  if (diffs) {
    col_suffix <- "_hdi_contains_true_value"
  } else {
    col_suffix <- "_hdi_contains_zero"
  }
  # For each _hdi column, create a new column with the check
  for (col in colnames(dt)) {
    if (is_vector_column(dt[[col]])) {
      # Step 1: parse the column safely
      parse_hdi_col(dt, col)
      new_col_name <- paste0(col, col_suffix)
      parsed_col <- dt[[col]]
      dt[, (new_col_name) :=
           sapply(parsed_col, function(x) x[1] <= 0 && x[2] >= 0)]
    }
  }
  return(dt)
}


#' Reshape (melt) the data.table such that each row ends up containing the 
#' value (posterior estimate, Rhat, HPDI) for one model parameter for one 
#' specific run.
#'
#' @param dt A data.table with each row containing all posterior parameter 
#' estimates for a specific run.
#' @param id_vars A vector of column names indicating by which columns should 
#' remain as identifiers upon melting.
#' @param contains_true_value_cols A logical value indicating if columns 
#' indicating whether the HPDI contains the true value (in the case of runs 
#' with simulated data) are present in "dt".
#' @param hdi_cols A logical value indicating whether columns containing HPDIs 
#' are present in "dt".
#' @param combination_cols A logical value indicating if columns with HPDI 
#' combination columns (in the case of simulated data runs, see function 
#' create_hdi_combinations) are present in "dt".
#'
#' @return A melted (wide to long format) version of dt.
reshape_dt <- function(dt,
                       id_vars = c("model_no_stan", "min_n_tests"),
                       contains_true_value_cols = FALSE,
                       hdi_cols = FALSE,
                       combination_cols = FALSE) {
  assert_that((combination_cols + contains_true_value_cols + hdi_cols) < 2)
  all_cols <- setdiff(names(dt), id_vars)
  base_stats <- all_cols[!grepl("_hdi", all_cols)]
  
  
  if (contains_true_value_cols) {
    current_cols <- paste0(base_stats, "_hdi_contains_true_value")
    hdi_cols <- paste0(base_stats, "_hdi")
    hdi50_cols <- paste0(base_stats, "_hdi50")
    
    stopifnot(all(current_cols %in% names(dt)))
    stopifnot(all(hdi_cols %in% names(dt)))
    stopifnot(all(hdi50_cols %in% names(dt)))
    
    dt_long <- melt(
      dt,
      id.vars = id_vars,
      measure.vars = list(
        value = base_stats,
        hdi_contains_true_value = current_cols,
        hdi = hdi_cols,
        hdi50 = hdi50_cols
      ),
      variable.name = "statistic"
    )
    # Split HDI cols into lower and upper
    dt_long[, c("lower", "upper") := transpose(hdi)]
    dt_long[, c("lower50", "upper50") := transpose(hdi50)]
  } else if (combination_cols) {
    current_cols <- paste0(base_stats, "_hdi_combination")
    hdi_cols <- paste0(base_stats, "_hdi")
    hdi50_cols <- paste0(base_stats, "_hdi50")
    
    stopifnot(all(current_cols %in% names(dt)))
    stopifnot(all(hdi_cols %in% names(dt)))
    stopifnot(all(hdi50_cols %in% names(dt)))
    
    dt_long <- melt(
      dt,
      id.vars = id_vars,
      measure.vars = list(
        value = base_stats,
        hdi_combination = current_cols,
        hdi = hdi_cols,
        hdi50 = hdi50_cols
      ),
      variable.name = "statistic"
    )
    dt_long[, hdi_combination := factor(hdi_combination)]
    # Split HDI cols into lower and upper
    dt_long[, c("lower", "upper") := transpose(hdi)]
    dt_long[, c("lower50", "upper50") := transpose(hdi50)]
  } else if (hdi_cols) {
    hdi_cols <- paste0(base_stats, "_hdi")
    stopifnot(all(hdi_cols %in% names(dt)))
    hdi50_cols <- paste0(base_stats, "_hdi50")
    
    stopifnot(all(hdi50_cols %in% names(dt)))
    
    dt_long <- melt(
      dt,
      id.vars = id_vars,
      measure.vars = list(
        value = base_stats,
        hdi = hdi_cols,
        hdi50 = hdi50_cols
      ),
      value.name = "value",
      variable.name = "statistic"
    )
    # Split HDI cols into lower and upper
    dt_long[, c("lower", "upper") := transpose(hdi)]
    dt_long[, c("lower50", "upper50") := transpose(hdi50)]
  } else {
    # No HDI columns
    dt_long <- melt(
      dt,
      id.vars = id_vars,
      measure.vars = list(value = base_stats),
      value.name = "value",
      variable.name = "statistic"
    )
  }
  ## Map integer index --> statistic name
  dt_long[, statistic := base_stats[statistic]]
  return(dt_long)
}


#' Return the names of columns that store posterior estimates of the 
#' interaction effect between hospitalization status and SARS-2 variant (for 
#' intercept, up-slope, and down-slope).
#'
#' @param log_prefix A string found as a prefix of column names indicating 
#' whether parameters are in log or natural space.
#' @param hdi A logical value indicating whether the names of columns 
#' containing HPDIs (as opposed to posterior point estimates) are requested.
#'
#' @return A named list of
#' A) group_cols_base: prefixes of columns storing interaction effect estimates,
#' B) group_cols_mapping: A named list with column names storing up-slope, 
#' intercept and down-slope related interaction effects parameters,
#' C) group_cols_mapping_labels: A named list of mappings from column name to 
#' label (to use in tables and figures).
return_hosp_variant_cols <- function(log_prefix, hdi = FALSE) {
  assert_that(log_prefix %in% c("log_", ""))
  
  suffix <- ifelse(hdi, "_hdi", "")
  group_cols_base <- paste0("hospitalized_", c("alpha", "delta", "omicron"))
  
  group_cols_base_labels <- paste0("Hospitalized, ", 
                                   c("Alpha", "Delta", "Omicron"))
  
  group_cols_intercept <- unname(sapply(group_cols_base,
                                        FUN = function(x) {
                                          return(paste0(log_prefix, 
                                                        "intercept_", 
                                                        x, 
                                                        suffix))
                                        }
  ))
  group_cols_slope_up <- unname(sapply(group_cols_base,
                                       FUN = function(x) {
                                         return(paste0(log_prefix, 
                                                       "slope_up_", 
                                                       x, 
                                                       suffix))
                                       }
  ))
  group_cols_slope_down <- unname(sapply(group_cols_base,
                                         FUN = function(x) {
                                           return(paste0(log_prefix, 
                                                         "slope_down_", 
                                                         x,
                                                         suffix))
                                         }
  ))
  
  group_cols_intercept_labels <- paste0("", group_cols_base_labels)
  group_cols_slope_up_labels <- paste0("", group_cols_base_labels)
  group_cols_slope_down_labels <- paste0("", group_cols_base_labels)
  
  group_cols_mapping <- list(
    slope_up = group_cols_slope_up, intercept = group_cols_intercept,
    slope_down = group_cols_slope_down
  )
  group_cols_label_mapping <- list(
    slope_up = setNames(group_cols_slope_up_labels, 
                       as.list(group_cols_slope_up)),
    intercept = setNames(group_cols_intercept_labels, 
                         as.list(group_cols_intercept)),
    slope_down = setNames(group_cols_slope_down_labels, 
                         as.list(group_cols_slope_down))
  )
  
  return(list(
    group_cols_base = group_cols_base, group_cols_mapping = group_cols_mapping,
    group_cols_mapping_labels = group_cols_label_mapping
  ))
}


#' Return the names of columns that store posterior estimates of fixed 
#' effects/infection-group-specific parameters (for intercept, up-slope, and 
#' down-slope).
#'
#' @param log_prefix A string found as a prefix of column names indicating 
#' whether parameters are in log or natural space.
#' @param hdi A logical value indicating whether the names of columns 
#' containing 94% HPDIs (as opposed to posterior point estimates) are requested.
#' @param hdi50 A logical value indicating whether the names of columns 
#' containing 50% HPDIs (as opposed to posterior point estimates) are requested.
#' @param latex A logical value indicating whether parameter labels are 
#' formatted for latex.
#'
#' @return A named list of
#' A) group_cols_base: prefixes of columns storing fixed-effects estimates,
#' B) group_cols_mapping: A named list with column names storing up-slope, 
#' intercept and down-slope related fixed effects parameters,
#' C) group_cols_mapping_labels: A named list of mappings from column name to 
#' label (to use in tables and figures).
return_group_cols <- function(log_prefix, 
                              hdi = FALSE, 
                              hdi50 = FALSE, 
                              latex = FALSE) {
  assert_that((hdi + hdi50) < 2)
  suffix <- ifelse(hdi, "_hdi", ifelse(hdi50, "_hdi50", ""))
  group_cols_base <- c(
    "alpha", "delta", "omicron", "hospitalized", "pams",
    "gender", "age_cat2", "age_cat3",
    "prior_infection", "hospitalized_alpha",
    "hospitalized_delta", "hospitalized_omicron"
  )
  if (latex) {
    group_cols_base_labels <- c(
      "Alpha", "Delta", "Omicron", "Hospitalized",
      "PAMS", "Male", "30-60 years",
      "\\textgreater{60} years", "Prior infection",
      "Hospitalized, Alpha", "Hospitalized, Delta",
      "Hospitalized, Omicron"
    )
  } else {
    group_cols_base_labels <- c(
      "Alpha", "Delta", "Omicron", "Hospitalized",
      "PAMS", "Male", "30-60 years", ">60 years",
      "Prior infection", "Hospitalized, Alpha",
      "Hospitalized, Delta", "Hospitalized, Omicron"
    )
  }
  
  group_cols_intercept <- unname(sapply(group_cols_base,
                                        FUN = function(x) {
                                          return(paste0(log_prefix, 
                                                        "intercept_", 
                                                        x, 
                                                        suffix))
                                        }
  ))
  group_cols_slope_up <- unname(sapply(group_cols_base,
                                       FUN = function(x) {
                                         return(paste0(log_prefix, 
                                                       "slope_up_", 
                                                       x, 
                                                       suffix))
                                       }
  ))
  group_cols_slope_down <- unname(sapply(group_cols_base,
                                         FUN = function(x) {
                                           return(paste0(log_prefix, 
                                                         "slope_down_", 
                                                         x, 
                                                         suffix))
                                         }
  ))
  
  group_cols_intercept_labels <- paste0("", group_cols_base_labels)
  group_cols_slope_up_labels <- paste0("", group_cols_base_labels)
  group_cols_slope_down_labels <- paste0("", group_cols_base_labels)
  
  group_cols_mapping <- list(
    slope_up = group_cols_slope_up, intercept = group_cols_intercept,
    slope_down = group_cols_slope_down
  )
  group_cols_label_mapping <- list(
    slope_up = setNames(group_cols_slope_up_labels, 
                       as.list(group_cols_slope_up)),
    intercept = setNames(group_cols_intercept_labels, 
                         as.list(group_cols_intercept)),
    slope_down = setNames(group_cols_slope_down_labels, 
                         as.list(group_cols_slope_down))
  )
  
  return(list(
    group_cols_base = group_cols_base,
    group_cols_mapping = group_cols_mapping,
    group_cols_mapping_labels = group_cols_label_mapping
  ))
}


#' Return the names of columns that store posterior HPDI-related estimates of 
#' fixed effects/infection-group-specific parameters (for intercept, up-slope, 
#' and down-slope).
#'
#' @param log_prefix A string found as a prefix of column names indicating 
#' whether parameters are in log or natural space.
#' @param include_base_cols A logical value indicating whether to also return 
#' base column names (prefixes) of fixed effects.
#' @param column_suffix A string specifying the suffix to be added to the base 
#' column names to generate the name of the column of interest.
#'
#' @return A named list of HPDI related column names (for up-slope, intercept, 
#' down-slope).
return_group_cols_hdi_mapping <- function(log_prefix, 
                                          include_base_cols = TRUE,
                                          column_suffix = "_hdi_contains_true_value") {
  assert_that(column_suffix %in% c(
    "_hdi_contains_true_value",
    "_hdi_combination", "_hdi", "_hdi50"
  ))
  
  group_cols <- return_group_cols(log_prefix)
  
  group_cols_base_suffix <- paste0(group_cols$group_cols_base, column_suffix)
  if (include_base_cols) {
    group_cols_base2 <- c(group_cols$group_cols_base, group_cols_base_suffix)
  } else {
    group_cols_base2 <- group_cols_base_suffix
  }
  group_cols_intercept2 <- unname(sapply(group_cols_base2,
                                         FUN = function(x) {
                                           return(paste0(log_prefix, 
                                                         "intercept_", 
                                                         x))
                                         }
  ))
  group_cols_slope_up2 <- unname(sapply(group_cols_base2,
                                        FUN = function(x) {
                                          return(paste0(log_prefix, 
                                                        "slope_up_", 
                                                        x))
                                        }
  ))
  group_cols_slope_down2 <- unname(sapply(group_cols_base2,
                                          FUN = function(x) {
                                            return(paste0(log_prefix, 
                                                          "slope_down_", 
                                                          x))
                                          }
  ))
  
  group_cols_mapping2 <- list(
    slope_up = group_cols_slope_up2, intercept = group_cols_intercept2,
    slope_down = group_cols_slope_down2
  )
  
  return(group_cols_mapping2)
}


#' Return the names of columns that store posterior predictions of up-slopes, 
#' peaks and down-slopes for different infection subgroups, computed from 
#' infection-level parameter estimates.
#'
#' @param log_prefix A string found as a prefix of column names indicating 
#' whether parameters are in log or natural space.
#' @param hdi A logical value indicating whether the names of columns 
#' containing 94% HPDIs (as opposed to posterior point estimates) are requested.
#' @param diff A logical value indicating whether estimates are given as 
#' differences to the respective reference group:
#' "non-PAMS, non-hospitalized" for disease severity status,
#' "pre-VOC" for SARS-CoV-2 variant,
#' "<30 years" for age group,
#' "no prior infection" for immunization status,
#' "female" for gender.
#' @param latex A logical value indicating whether labels are formatted for 
#' latex.
#'
#' @return A named list of
#' A) group_cols_base: prefixes of column names storing predictions for the 
#' different infection groups,
#' B) group_cols_mapping: A named list with column names storing predictions of 
#' up-slope, intercept and down-slopes for the different infection subgroups,
#' C) group_cols_mapping_labels: A named list of mappings from column name to 
#' label (to use in tables and figures).
return_group_avg_cols <- function(log_prefix, 
                                  hdi = FALSE, 
                                  latex = FALSE, 
                                  diff = FALSE) {
  if (diff) {
    suffix <- ifelse(hdi, "avg_diff_hdi", "_avg_diff")
  } else {
    suffix <- ifelse(hdi, "_avg_hdi", "_avg")
  }
  if (diff) {
    group_cols_base <- c(
      "alpha", "delta", "omicron", "hospitalized", "pams",
      "gender", "age_cat2", "age_cat3",
      "prior_infection"
    )
  } else {
    group_cols_base <- c(
      "wildtype", "alpha", "delta", "omicron",
      "non_pams_non_hosp", "hospitalized", "pams",
      "female", "male", "age_cat1", "age_cat2", "age_cat3",
      "no_prior_infection", "prior_infection"
    )
  }
  
  if (latex) {
    if (diff) {
      group_cols_base_labels <- c(
        "Alpha", "Delta", "Omicron", "Hospitalized",
        "PAMS", "Male", "30-60 years",
        "\\textgreater{60} years", "Prior infection"
      )
    } else {
      group_cols_base_labels <- c(
        "Pre-VOC", "Alpha", "Delta", "Omicron",
        "Non-PAMS, non-hosp.", "Hospitalized",
        "PAMS", "Female", "Male", "\\textless{30} years",
        "30-60 years", "\\textgreater{60} years",
        "First infection", "Reinfection"
      )
    }
  } else {
    if (diff) {
      group_cols_base_labels <- c(
        "Alpha", "Delta", "Omicron", "Hospitalized",
        "PAMS", "Male", "30-60 years",
        ">60 years", "Prior infection"
      )
    } else {
      group_cols_base_labels <- c(
        "Pre-VOC", "Alpha", "Delta", "Omicron",
        "Non-PAMS, non-hosp.", "Hospitalized",
        "PAMS", "Female", "Male", "<30 years",
        "30-60 years", ">60 years",
        "First infection", "Reinfection"
      )
    }
  }
  
  group_cols_intercept <- unname(sapply(group_cols_base,
                                        FUN = function(x) {
                                          return(paste0(log_prefix, 
                                                        "intercept_", 
                                                        x, 
                                                        suffix))
                                        }
  ))
  group_cols_up <- unname(sapply(group_cols_base,
                                 FUN = function(x) {
                                   return(paste0(log_prefix, 
                                                 "time2peak_", 
                                                 x, 
                                                 suffix))
                                 }
  ))
  group_cols_down <- unname(sapply(group_cols_base,
                                   FUN = function(x) {
                                     return(paste0(log_prefix, 
                                                   "time_from_peak_", 
                                                   x, 
                                                   suffix))
                                   }
  ))
  group_cols_slope_up <- unname(sapply(group_cols_base,
                                       FUN = function(x) {
                                         return(paste0(log_prefix, 
                                                       "slope_up_", 
                                                       x, 
                                                       suffix))
                                       }
  ))
  group_cols_slope_down <- unname(sapply(group_cols_base,
                                         FUN = function(x) {
                                           return(paste0(log_prefix, 
                                                         "slope_down_", 
                                                         x, 
                                                         suffix))
                                         }
  ))
  
  group_cols_intercept_labels <- paste0("", group_cols_base_labels)
  group_cols_up_labels <- paste0("", group_cols_base_labels)
  group_cols_down_labels <- paste0("", group_cols_base_labels)
  group_cols_slope_up_labels <- paste0("", group_cols_base_labels)
  group_cols_slope_down_labels <- paste0("", group_cols_base_labels)
  
  group_cols_mapping <- list(
    up = group_cols_up,
    intercept = group_cols_intercept,
    down = group_cols_down,
    slope_up = group_cols_slope_up,
    slope_down = group_cols_slope_down
  )
  group_cols_label_mapping <- list(
    up = setNames(group_cols_up_labels, as.list(group_cols_up)),
    intercept = setNames(group_cols_intercept_labels, 
                         as.list(group_cols_intercept)),
    down = setNames(group_cols_down_labels, as.list(group_cols_down)),
    slope_up = setNames(group_cols_slope_up_labels, 
                        as.list(group_cols_slope_up)),
    slope_down = setNames(group_cols_slope_down_labels, 
                          as.list(group_cols_slope_down))
  )
  
  return(list(
    group_cols_base = group_cols_base,
    group_cols_mapping = group_cols_mapping,
    group_cols_mapping_labels = group_cols_label_mapping
  ))
}


#' Return the names of columns that store posterior predictions of up-slopes, 
#' peaks and down-slopes for different infection subgroups, computed from 
#' population-level-parameter estimates.
#'
#' @param log_prefix A string found as a prefix of column names indicating 
#' whether parameters are in log or natural space.
#' @param hdi A logical value indicating whether the names of columns 
#' containing 94% HPDIs (as opposed to posterior point estimates) are requested.
#' @param diff A logical value indicating whether estimates are given as 
#' differences to the respective reference group:
#' "non-PAMS, non-hospitalized" for disease severity status,
#' "pre-VOC" for SARS-CoV-2 variant,
#' "<30 years" for age group,
#' "no prior infection" for immunization status,
#' "female" for gender.
#' @param latex A logical value indicating whether parameter labels are 
#' formatted for latex.
#'
#' @return A named list of
#' A) group_names: The column names storing predictions for the different 
#' infection subgroups.
#' B) group_labels: The corresponding labels (to use in tables and figures).
return_group_post_pred_cols <- function(latex = FALSE, diff = FALSE) {
  if (diff) {
    groups <- c(
      "Alpha", "Delta", "Omicron", "Hospitalized", "PAMS", "Male",
      "30-60", ">60", "Prior infection"
    )
  } else {
    groups <- c(
      "Wildtype", "Alpha", "Delta", "Omicron",
      "Non-PAMS, non-hosp.", "Hospitalized", "PAMS",
      "Female", "Male", "<30", "30-60", ">60",
      "First infection", "Prior infection"
    )
  }
  if (latex) {
    if (diff) {
      group_labels <- c(
        "Alpha", "Delta", "Omicron",
        "Hospitalized", "PAMS", "Male", "30-60 years",
        "\\textgreater{60} years", "Reinfection"
      )
    } else {
      group_labels <- c(
        "Pre-VOC", "Alpha", "Delta", "Omicron",
        "Non-PAMS, non-hosp.", "Hospitalized", "PAMS",
        "Female", "Male", "\\textless{30} years",
        "30-60 years", "\\textgreater{60} years",
        "First infection", "Reinfection"
      )
    }
  } else {
    if (diff) {
      group_labels <- c(
        "Alpha", "Delta", "Omicron",
        "Hospitalized", "PAMS", "Male", "30-60 years",
        ">60 years", "Reinfection"
      )
    } else {
      group_labels <- c(
        "Pre-VOC", "Alpha", "Delta", "Omicron",
        "Non-PAMS, non-hosp.", "Hospitalized", "PAMS",
        "Female", "Male", "<30 years",
        "30-60 years", ">60 years",
        "First infection", "Reinfection"
      )
    }
  }
  return(list(group_names = groups, group_labels = group_labels))
}


#' Return the names of columns storing variance parameter values.
#'
#' @param log_prefix A string found as a prefix of column names indicating 
#' whether parameters are in log or natural space.
#'
#' @return A vector of column names storing variance parameter values.
return_sigma_cols <- function(log_prefix = "") {
  assert_that(log_prefix %in% c("log_", ""))
  base_cols <- c("intercept", "slope_up", "slope_down")
  sigma_cols <- paste0(log_prefix, base_cols, "_sigma")
  return(sigma_cols)
}


#' Return a long-format data.table with Rhat values.
#' @param sim_data A logical value indicating whether values are from runs with 
#' simulated data.
#' @param param_sizes_all_zero A logical value indicating whether to extract 
#' only values from simulation runs where all fixed effects 
#' (infection-group-specific parameters) were set to 0 or only those where
#' these values were different from 0.
#' @param model_nos A vector of model numbers to include (refers to stan file 
#' numbers).
#'
#' @return A long-format version of the data.table containing Rhat values.
return_rhat_reshaped <- function(sim_data = FALSE, param_sizes_all_zero = FALSE,
                                 model_nos = c(1, 2, 4, 5)) {
  if (sim_data == FALSE) {
    assert_that(param_sizes_all_zero == FALSE)
  }
  if (sim_data) {
    if (param_sizes_all_zero) {
      assert_that(setequal(model_nos, c(1, 2, 3)))
      model_nos <- c(1, 2, 3)
      model_levels <- c(1, 2, 3)
      model_labels <- c("Model 1", "Model 2", "Model 1_2")
    } else {
      assert_that(setequal(model_nos, c(1, 2)))
      model_nos <- c(1, 2)
      model_levels <- c(1, 2)
      model_labels <- c("Model 1", "Model 2")
    }
    dtrhat <- process_sim_dt(fread(PARAMS_RHAT_FILE, sep = "\t"),
                             param_sizes_all_zero = param_sizes_all_zero,
                             rhat = TRUE
    )
    dtrhat <- dtrhat[(model_no_stan %in% model_nos) &
                       (param_sizes_all_zero == param_sizes_all_zero)]
  } else {
    assert_that((setequal(model_nos, c(1, 2))) |
                  (setequal(model_nos, c(1, 2, 4, 5))))
    if (setequal(model_nos, c(1, 2))) {
      model_nos <- c(1, 2)
      model_levels <- c(1, 2)
      model_labels <- c("Model 1", "Model 2")
    } else {
      model_nos <- c(1, 2, 4, 5)
      model_levels <- c(1, 2, 4, 5)
      model_labels <- c("Model 1", "Model 2", "Model 1", "Model 2")
    }
    filename_suffix <- ""
    scales_plot <- "fixed"
    
    dtrhat <- process_dt(fread(PARAMS_RHAT_FILE, sep = "\t"), rhat = TRUE)
    dtrhat <- dtrhat[model_no_stan %in% model_nos]
  }
  
  dtrhat[, model_no_stan := factor(model_no_stan, levels = model_levels)]
  id_vars <- c("model_no_stan", "min_tests", "min_n_pos_tests", "iter_sampling")
  cols_with_values <- grep("intercept|slope_up|slope_down|sigma|alpha",
                           names(dtrhat),
                           value = TRUE
  )
  cols_to_keep <- c(id_vars, cols_with_values)
  dtrhat_reshaped <- reshape_dt(dtrhat[, ..cols_to_keep], id_vars = id_vars)
  dtrhat_reshaped[, min_tests_factor := factor(min_tests)]
  dtrhat_reshaped[, min_n_pos_tests_factor := factor(min_n_pos_tests)]
  dtrhat_reshaped[, main_param := factor(
    fifelse(
      grepl("intercept", statistic), "Peak",
      fifelse(
        grepl("slope_up", statistic), "Proliferation rate",
        fifelse(
          grepl("slope_down", statistic), "Clearance rate",
          NA_character_
        )
      )
    ),
    levels = c("Proliferation rate", "Peak", "Clearance rate")
  )]
  dtrhat_reshaped[, model_no_plot := factor(model_no_stan, 
                                            levels = model_levels, 
                                            labels = model_labels)]
  
  return(dtrhat_reshaped)
}


#' Return Rhat and ESS values of population-level parameters (either up-slope, 
#' intercept, or down-slope) for a specific run.
#'
#' @param ss A data.table containing summary statistics (including Rhat and ESS 
#' values) from a specific run.
#' @param regex_string A regular expression for selecting rows with estimates 
#' of a specific parameter from ss.
#' @param param_string A string used for removing the column prefix indicating 
#' if a parameter is up-slope, intercept or down-slope-related (only applies to 
#' population-wide location and scale parameters, such as "log_intercept_mu" or 
#' "intercept_sigma").
#' @param min_tests The minimum number of PCR days per infection that were 
#' required for inclusiogn in the run producing estimates stored in ss.
#' @param model_no A model number to filter by (refers to the number of the 
#' stan file).
#'
#' @return A data.table containing Rhat and ESS values for either up-slope,
#' intercept, or down-slope related population-level parameters
process_main_param_rhat <- function(ss, main_param, min_tests,
                                    model_no) {
  assert_that(main_param %in% c("Proliferation rate", "Peak", "Clearance rate",
                                "Other"))
  param_string_mapping <- list(
    "Proliferation rate" = "slope_up_",
    "Peak" = "intercept_",
    "Clearance rate" = "slope_down_",
    "Other" = ""
  )
  param_regex_mapping <- list(
    "Proliferation rate" = "_slope_up|slope_up_sigma",
    "Peak" = "_intercept|intercept_sigma",
    "Clearance rate" = "_slope_down|slope_down_sigma",
    "Other" = "sigma_sigma|sigma_mu|alpha_mu"
  )
  param_str <- param_string_mapping[main_param]
  param_regex <- param_regex_mapping[main_param]
  main_param <- ifelse(main_param == "Other", NA, main_param)

  # Process up-slopes
  param_ss <- ss[grepl(eval(param_regex), variable)]
  # Extract the number, or keep the original value if no match
  param_ss[, number := ifelse(grepl("\\[\\d+\\]", variable),
                              gsub(".*\\[(\\d+)\\].*", "\\1", variable),
                              variable
  )]
  param_ss[, statistic := ifelse(number %in% names(SUBGRP_MAPPING),
                                 SUBGRP_MAPPING[number], number
  )]
  # Remove param_string and log_ prefixes from parameter name
  param_ss[, statistic := gsub(eval(param_str), "", statistic)]
  param_ss[, statistic := gsub("log_", "", statistic)]
  param_ss[, statistic := as.character(SUBGRPS_RHAT_LATEX[statistic])]
  
  dt_param <- data.table(
    param = param_ss$statistic,
    rhat = param_ss$rhat,
    ess_bulk = param_ss$ess_bulk,
    ess_tail = param_ss$ess_tail,
    param_orig = param_ss$variable
  )
  dt_param[, main_param := eval(main_param)]
  dt_param[, min_tests := eval(min_tests)]
  dt_param[, min_n_pos_tests := 2]
  dt_param[, model_no_stan := eval(model_no)]
  
  return(dt_param)
}

#' Return Rhat and ESS values of population-level up-slope, intercept, and 
#' down-slope parameters.
#'
#' @return A data.table containing Rhat and ESS values for up-slope, intercept, 
#' and down-slope related population-level parameters.
return_rhat_reshaped_main_n_groups <- function() {
  model_nos <- c(
    1, 1, 2, 2, 1, 1,
    2, 2, 4, 4, 5, 5
  )
  sampling_iterations <- c(
    1000, 1000, 1000, 1000,
    2000, 2000, 2000, 2000,
    2000, 2000, 2000, 2000
  )
  min_tests <- c(3, 6, 3, 6, 3, 6, 3, 6, 3, 6, 3, 6)
  counter <- 1
  data_tables <- c()
  
  for (ss_file in c(
    SS_MODEL1_3_PCRS_FILE,
    SS_MODEL1_6_PCRS_FILE,
    SS_MODEL2_3_PCRS_FILE,
    SS_MODEL2_6_PCRS_FILE,
    SS_MODEL1_3_PCRS_ITER2000_FILE,
    SS_MODEL1_6_PCRS_ITER2000_FILE,
    SS_MODEL2_3_PCRS_ITER2000_FILE,
    SS_MODEL2_6_PCRS_ITER2000_FILE,
    SS_MODEL4_3_PCRS_ITER2000_FILE,
    SS_MODEL4_6_PCRS_ITER2000_FILE,
    SS_MODEL5_3_PCRS_ITER2000_FILE,
    SS_MODEL5_6_PCRS_ITER2000_FILE
  )) {
    ss <- fread(ss_file, sep = "\t")
    # Process intercepts
    dt_peak <- process_main_param_rhat(
      ss = ss,
      main_param = "Peak",
      min_tests = min_tests[counter],
      model_no = model_nos[counter]
    )
    
    
    # Process up-slopes
    dt_up <- process_main_param_rhat(
      ss = ss,
      main_param = "Proliferation rate",
      min_tests = min_tests[counter],
      model_no = model_nos[counter]
    )
    
    
    # Process down-slopes
    dt_down <- process_main_param_rhat(
      ss = ss,
      main_param = "Clearance rate",
      min_tests = min_tests[counter],
      model_no = model_nos[counter]
    )
    
    dt_other <- process_main_param_rhat(
      ss = ss,
      main_param = "Other",
      min_tests = min_tests[counter],
      model_no = model_nos[counter]
    )
    
    merged_dt <- rbindlist(list(dt_peak, dt_up, dt_down, dt_other), fill = TRUE)
    merged_dt[, iter_sampling := eval(sampling_iterations[counter])]
    
    data_tables <- c(data_tables, list(merged_dt))
    counter <- counter + 1
  }
  dtrhat_reshaped <- rbindlist(data_tables)
  
  model_nos <- c(1, 2, 4, 5)
  model_levels <- c(1, 2, 4, 5)
  model_labels <- c("Model 1", "Model 2", "Model 1", "Model 2")
  filename_suffix <- ""
  scales_plot <- "fixed"
  
  dtrhat_reshaped[, model_no_stan := factor(model_no_stan, 
                                            levels = model_levels)]
  dtrhat_reshaped[, min_tests_factor := factor(min_tests)]
  dtrhat_reshaped[, main_param := factor(main_param, 
                                         levels = c("Proliferation rate", 
                                                    "Peak", 
                                                    "Clearance rate"))]
  dtrhat_reshaped[, model_no_plot := factor(model_no_stan, 
                                            levels = model_levels, 
                                            labels = model_labels)]
  
  return(dtrhat_reshaped)
}


#' Return Rhat and ESS values of infection-level up-slope, intercept, and 
#' down-slope parameters.
#'
#' @return A data.table containing Rhat and ESS values for up-slope, intercept, 
#' and down-slope related infection-level parameters.
return_rhat_reshaped_infection <- function() {
  model_nos <- c(1, 1, 2, 2, 1, 1, 2, 2)
  sampling_iterations <- c(
    1000, 1000, 1000, 1000,
    2000, 2000, 2000, 2000
  )
  min_tests <- c(3, 6, 3, 6, 3, 6, 3, 6)
  counter <- 1
  data_tables <- c()
  for (ss_file in c(
    SS_MODEL1_3_PCRS_FILE,
    SS_MODEL1_6_PCRS_FILE,
    SS_MODEL2_3_PCRS_FILE,
    SS_MODEL2_6_PCRS_FILE,
    SS_MODEL1_3_PCRS_ITER2000_FILE,
    SS_MODEL1_6_PCRS_ITER2000_FILE,
    SS_MODEL2_3_PCRS_ITER2000_FILE,
    SS_MODEL2_6_PCRS_ITER2000_FILE
  )) {
    ss <- fread(ss_file, sep = "\t")
    # Process intercepts
    intercept_ss <- ss[grepl("^intercept\\[", variable)]
    dt_intercept <- data.table(
      param = intercept_ss$variable,
      value = intercept_ss$rhat
    )
    dt_intercept[, main_param := "Peak"]
    dt_intercept[, min_tests := eval(min_tests[counter])]
    dt_intercept[, model_no_stan := eval(model_nos[counter])]
    
    
    # Process up-slopes
    up_ss <- ss[grepl("^slope_up\\[", variable)]
    dt_up <- data.table(
      param = up_ss$variable,
      value = up_ss$rhat
    )
    dt_up[, main_param := "Proliferation rate"]
    dt_up[, min_tests := eval(min_tests[counter])]
    dt_up[, model_no_stan := eval(model_nos[counter])]
    
    
    # Process down-slopes
    down_ss <- ss[grepl("^slope_down\\[", variable)]
    dt_down <- data.table(
      param = down_ss$variable,
      value = down_ss$rhat
    )
    dt_down[, main_param := "Clearance rate"]
    dt_down[, min_tests := eval(min_tests[counter])]
    dt_down[, model_no_stan := eval(model_nos[counter])]
    merged_dt <- rbindlist(list(dt_intercept, dt_up, dt_down), fill = TRUE)
    merged_dt[, iter_sampling := eval(sampling_iterations[counter])]
    
    data_tables <- c(data_tables, list(merged_dt))
    counter <- counter + 1
  }
  
  dtrhat_reshaped <- rbindlist(data_tables)
  
  model_nos <- c(1, 2)
  model_levels <- c(1, 2)
  model_labels <- c("Model 1", "Model 2")
  filename_suffix <- ""
  scales_plot <- "fixed"
  
  dtrhat_reshaped[, model_no_stan := factor(model_no_stan, 
                                            levels = model_levels)]
  dtrhat_reshaped[, min_tests_factor := factor(min_tests)]
  dtrhat_reshaped[, main_param := factor(main_param, 
                                         levels = c("Proliferation rate", 
                                                    "Peak", 
                                                    "Clearance rate"))]
  dtrhat_reshaped[, model_no_plot := factor(model_no_stan, 
                                            levels = model_levels, 
                                            labels = model_labels)]
  
  return(dtrhat_reshaped)
}


#' Return the body of a latex table presenting posterior estimates or Rhat 
#' values across model parameters and runs.
#'
#' @param dt A data.table with a row of values for each parameter and run.
#' @param name_mapping A named list mapping from parameter name in "dt" to the 
#' name used in the latex table.
#' @param rhat A logical value indicating whether Rhat values are being 
#' processed.
#'
#' @return A string representing the body of a latex table with either 
#' posterior estimates or Rhat values across model parameters and runs.
create_latex_body_real_data <- function(dt, name_mapping, rhat = FALSE) {
  param_colname <- ifelse(rhat, "param", "statistic")
  cline_string <- ifelse(rhat, "\\cline{2-8}", "\\cline{2-8}")
  round_n_digits <- ifelse(rhat, 3, 2)
  latex_body <- c()
  for (column in names(name_mapping)) {
    latex_string <- paste0("\\multirow{2}{*}{", name_mapping[column], "}")
    for (model_no in model_nos) {
      model_label <- paste("Model", model_no)
      if (model_label == "Model 1") {
        latex_string <- paste(latex_string, model_no, sep = " & ")
      } else {
        latex_string <- paste(latex_string, paste0("\\\\ ", cline_string, 
                                                   "\n & ", model_no), 
                              sep = "\n")
      }
      for (n_tests in c(3, 6)) {
        row <- dt[(get(param_colname) == column) & (model == model_label) & (min_n_tests == n_tests)]
        if (rhat) {
          new_string <- paste(
            sprintf(paste0("%.", round_n_digits, "f"), row$rhat),
            sprintf("%.0f", row$ess_bulk),
            sprintf("%.0f", row$ess_tail),
            sep = " & "
          )
          latex_string <- paste(latex_string, new_string, sep = " & ")
        } else {
          # Enclose by $ to display minus signs correctly in latex
          new_string <- paste(
            paste0("$", sprintf(paste0("%.", round_n_digits, "f"), row$value), 
                   "$"),
            paste0(
              "$",
              sprintf(paste0("%.", round_n_digits, "f"), row$lower),
              "$",
              " & ",
              "$",
              sprintf(paste0("%.", round_n_digits, "f"), row$upper),
              "$"
            ),
            sep = " & "
          )
          latex_string <- paste(latex_string, new_string, sep = " & ")
        }
      }
    }
    latex_string <- paste0(latex_string, "\n")
    latex_body <- c(latex_body, latex_string)
  }
  latex_body <- paste(latex_body, 
                      collapse = "\\\\ \\specialrule{0.75pt}{0pt}{0pt}\n")
  latex_body <- paste0(latex_body, "\\\\ \\hline")
  
  return(latex_body)
}


#' Summarize viral load (log10_load) by a grouping variable.
#'
#' @param dt A `data.table` containing log10_load and the grouping variable.
#' @param by The column name to group by.
#' @param digits The number of decimal places for rounding (default: 1).
#' @param avg The summary statistic to use (default: "median").
#
#' @return A data table with summary statistics (Count, Median, Q1, Q3, IQR) 
#' for each group.
summary_dt_vl <- function(dt, by, digits = 1, avg = "median") {
  dec_places_formatter <- paste0("%.", digits, "f")
  summary_dt <- dt[!is.na(by) & !is.na(log10_load)] %>%
    group_by(get(by)) %>%
    summarise(
      Count = n(),
      Median = sprintf(dec_places_formatter, 
                       round(median(log10_load, na.rm = TRUE), digits)),
      Q1 = round(quantile(log10_load, 0.25, na.rm = TRUE), digits),
      Q3 = round(quantile(log10_load, 0.75, na.rm = TRUE), digits),
      IQR = paste0(sprintf(dec_places_formatter, Q1), "-", 
                   sprintf(dec_places_formatter, Q3))
    ) %>%
    setnames(old = "get(by)", new = "Variable")
  return(summary_dt)
}


#' Fit brms model for first viral load analysis (skewnormal error distribution).
#'
#' @param model A brm model formula.
#' @param data The data to be used.
#' @param fn The names of the model.
#' @param adapt_delta A parameter for estimation in Stan.
#' @return A brms fit object.
bfit_first_vl_model <- function(formula,
                                data,
                                priors,
                                fn,
                                family,
                                stanvars = NULL,
                                redo = TRUE) {
  assert_that(family %in% c("student-t", "skewnormal"))
  filename_suffix <- paste0("_", family)
  if (family == "skewnormal") {
    family <- skew_normal()
  } else {
    family <- student()
  }
  fn <- file.path(DIR_FIT_FIRST_VL, paste0(fn, filename_suffix, ".Rdata"))
  if (file.exists(fn) && !redo) {
    load(fn)
  } else {
    Bfit <- brm(
      formula = formula,
      family = family,
      data = data,
      stanvars = stanvars,
      chains = 4,
      cores = 4,
      iter = 2000,
      warmup = 1000,
      seed = SEED,
      prior = priors
    )
    
    draws <- as_draws(Bfit$fit)
    sampler_params <-
      nuts_params(Bfit) %>%
      data.table() %>%
      dcast(Chain + Iteration ~ Parameter, value.var = "Value")
    save(Bfit, draws, sampler_params, file = fn)
  }
  return(Bfit)
}
