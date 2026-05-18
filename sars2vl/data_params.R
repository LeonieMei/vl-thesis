library(here)

TOP <- here()

SEED <- 123
TOP_OUTPUT <- file.path(TOP, "output")
DIR_DATA <- file.path(TOP, "data")
DIR_PDATA <- file.path(TOP_OUTPUT, "pdata")
DIR_PDATA_SIM <- DIR_PDATA
DIR_FIGURES <- file.path(TOP_OUTPUT, "figures")
DIR_FIT <- file.path(TOP_OUTPUT, "fits")
DIR_TABLES <- file.path(TOP_OUTPUT, "tables")
DIR_SS_DATA <- file.path(DIR_PDATA, "ss")
DIR_FIT_FIRST_VL <- file.path(DIR_FIT, "first_vl")
DIR_STAN <- file.path(TOP, "stan")

TC_FILE_JSON <- file.path(TOP, "data", "min-3-timeseries-paper2-lines.json")
TC_FILE_JSON_N_GENE <- file.path(TOP_OUTPUT, "data", "min-3-timeseries-paper2-n-gene-lines.json")


VL_FIRST_POS_FILE_TSV <- file.path(TOP, "data", "viral-load-with-negatives-first-positive-paper2.tsv")
VL_FIRST_POS_FILE_TSV_N_GENE <- file.path(TOP_OUTPUT, "data", "viral-load-with-negatives-first-positive-paper2-n-gene.tsv")


PROBLEMATIC_INFECTIONS_FILE <- file.path(DIR_DATA, "problematicInfections.tsv")
HIGH_VL_ERROR_INFECTIONS_FILE <- file.path(TOP, "data", "high_error_infections.tsv")

PARAMS_FILE <- file.path(DIR_PDATA_SIM, "params_stats.tsv")
PARAMS_EXP_FILE <- file.path(DIR_PDATA_SIM, "params_stats_exp.tsv")
PARAMS_HDI_FILE <- file.path(DIR_PDATA_SIM, "params_stats_hdi.tsv")
PARAMS_HDI50_FILE <- file.path(DIR_PDATA_SIM, "params_stats_hdi50.tsv")
PARAMS_EXP_HDI_FILE <- file.path(DIR_PDATA_SIM, "params_stats_exp_hdi.tsv")
PARAMS_EXP_HDI50_FILE <- file.path(DIR_PDATA_SIM, "params_stats_exp_hdi50.tsv")
PARAMS_DIFF_FILE <- file.path(DIR_PDATA_SIM, "params_diff_stats.tsv")
PARAMS_REL_DIFF_FILE <- file.path(DIR_PDATA_SIM, "params_diff_stats_exp.tsv")
PARAMS_DIFF_HDI_FILE <- file.path(DIR_PDATA_SIM, "params_diff_stats_hdi.tsv")
PARAMS_DIFF_HDI50_FILE <- file.path(DIR_PDATA_SIM, "params_diff_stats_hdi50.tsv")
PARAMS_REL_DIFF_HDI_FILE <- file.path(DIR_PDATA_SIM, "params_diff_stats_exp_hdi.tsv")
PARAMS_REL_DIFF_HDI50_FILE <- file.path(DIR_PDATA_SIM, "params_diff_stats_exp_hdi50.tsv")
PARAMS_GROUP_POST_PRED_WT_HOSP_FILE <- file.path(DIR_PDATA_SIM, "params_stats_post_pred_groups_ref_variant_wildtype_hosp.tsv")
PARAMS_GROUP_POST_PRED_OMI_HOSP_FILE <- file.path(DIR_PDATA_SIM, "params_stats_post_pred_groups_ref_variant_omicron_hosp.tsv")
PARAMS_GROUP_POST_PRED_WT_PAMS_FILE <- file.path(DIR_PDATA_SIM, "params_stats_post_pred_groups_ref_variant_wildtype_pams.tsv")
PARAMS_GROUP_POST_PRED_OMI_PAMS_FILE <- file.path(DIR_PDATA_SIM, "params_stats_post_pred_groups_ref_variant_omicron_pams.tsv")
PARAMS_GROUP_AVG_POST_PRED_FILE <- file.path(DIR_PDATA_SIM, "params_stats_post_pred_inf.tsv")
PARAMS_GROUP_AVG_POST_PRED_HDI_FILE <- file.path(DIR_PDATA_SIM, "params_stats_post_pred_inf_hdi.tsv")
PARAMS_GROUP_AVG_POST_PRED_HDI50_FILE <- file.path(DIR_PDATA_SIM, "params_stats_post_pred_inf_hdi50.tsv")

PARAMS_RHAT_FILE <- file.path(DIR_PDATA_SIM, "params_rhat_stats.tsv")

SS_MODEL2_3_PCRS_FILE <- file.path(DIR_SS_DATA, "model2_threading_sel3_chains4_n_iter1000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL2_6_PCRS_FILE <- file.path(DIR_SS_DATA, "model2_threading_sel6_chains4_n_iter1000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL1_3_PCRS_FILE <- file.path(DIR_SS_DATA, "model1_threading_sel3_chains4_n_iter1000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL1_6_PCRS_FILE <- file.path(DIR_SS_DATA, "model1_threading_sel6_chains4_n_iter1000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")

SS_MODEL2_3_PCRS_ITER2000_FILE <- file.path(DIR_SS_DATA, "model2_threading_sel3_chains4_n_iter2000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL2_6_PCRS_ITER2000_FILE <- file.path(DIR_SS_DATA, "model2_threading_sel6_chains4_n_iter2000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL1_3_PCRS_ITER2000_FILE <- file.path(DIR_SS_DATA, "model1_threading_sel3_chains4_n_iter2000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL1_6_PCRS_ITER2000_FILE <- file.path(DIR_SS_DATA, "model1_threading_sel6_chains4_n_iter2000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL4_3_PCRS_ITER2000_FILE <- file.path(DIR_SS_DATA, "model4_threading_sel3_chains4_n_iter2000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL4_6_PCRS_ITER2000_FILE <- file.path(DIR_SS_DATA, "model4_threading_sel6_chains4_n_iter2000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL5_3_PCRS_ITER2000_FILE <- file.path(DIR_SS_DATA, "model5_threading_sel3_chains4_n_iter2000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")
SS_MODEL5_6_PCRS_ITER2000_FILE <- file.path(DIR_SS_DATA, "model5_threading_sel6_chains4_n_iter2000_T2_trim_neg_pcrs_14_days_to_negative_exclude_lrt_samples_min_max_load_5_ss_file.tsv")


GRP_VAR <- "ID"

HOSPITAL_TEST_CENTRE <- c("H", "ICU", "IDW", "WD")
PAMS_TEST_CENTRE <- c("C19")
VARIANT_LEVELS <- c("Wildtype", "Alpha", "Delta", "Omicron", "Unknown")
VARIANT_LABELS <- c("Wildtype", "Alpha", "Delta", "Omicron", "Unknown")
VARIANT_LEGEND_LABELS <- list("Wildtype" = "Pre-VOC", "Alpha" = "Alpha", "Delta" = "Delta", "Omicron" = "Omicron", "Unknown" = "Unknown")

VARIANT2_LEVELS <- c("Wildtype", "Alpha", "Delta", "BA1", "BA2", "BA5", "BA2Desc", "JN1", "JN1Desc", "XEC", "OmicronUnknown", "Unknown")
VARIANT2_LABELS <- c("Wildtype", "Alpha", "Delta", "BA1", "BA2", "BA5", "BA2Desc", "JN1", "JN1Desc", "XEC", "OmicronUnknown", "Unknown")
VARIANT2_LEVELS2 <- c("Wildtype", "Alpha", "Delta", "BA1", "BA2", "BA5", "BA2Desc", "Unknown")
VARIANT2_LABELS2 <- c("Wildtype", "Alpha", "Delta", "BA1", "BA2", "BA5", "BA2Desc", "Unknown")

VARIANT_ALL_LABELS <- c("Wildtype", "Alpha", "Delta", "Omicron", "BA1", "BA2", "BA5", "BA2Desc", "JN1", "JN1Desc", "XEC", "OmicronUnknown", "Unknown")
VARIANT2_LEGEND_LABELS <- c(
  "Wildtype" = "Pre-VOC",
  "Alpha" = "Alpha",
  "Delta" = "Delta",
  "BA1" = "BA.1",
  "BA2" = "BA.2",
  "BA5" = "BA.5",
  "BA2Desc" = "BA.2 descendant",
  "JN1" = "JN.1",
  "JN1Desc" = "JN.1 descendant",
  "XEC" = "XEC",
  "OmicronUnknown" = "Omicron unknown",
  "Unknown" = "Unknown"
)
OMICRON_SUBVARIANT_LEVELS <- c("BA1", "BA2", "BA5", "BA2Desc", "JN1", "JN1Desc", "XEC", "OmicronUnknown")
OMICRON_SUBVARIANT_LEVELS_BA4 <- c("BA1", "BA2", "BA4", "BA5", "BA2Desc")
OMICRON_SUBVARIANT_LEVELS_NO_BA4 <- c("BA1", "BA2", "BA5", "BA2Desc")
NON_OMICRON_LEVELS <- c("Wildtype", "Alpha", "Delta", "Unknown")

MATERIAL_LEVELS <- c("N", "NP", "NPD", "P", "PD", "URT", "Unknown")
MATERIAL_LABELS <- c("Nasal", "Nasopharyngeal", "Nasopharyngeal (deep)", "Pharyngeal", "Pharyngeal (deep)", "URT", "Unknown")


AGE_BREAKS <- c(seq(0, 25, 5), seq(35, 65, 10), 120)
AGE_BREAKS2 <- c(0, 30, 50, 120)
AGE_BREAKS3 <- c(0, 30, 60, 120)
AGE_LEVELS2 <- c("[0,30)", "[30,50)", "[50,120)")
AGE_LEVELS3 <- c("[0,30)", "[30,60)", "[60,120)")
AGE_LABELS <- c("0-5", "5-10", "10-15", "15-20", "20-25", "25-35", "35-45", "45-55", "55-65", ">65")
AGE_LABELS2 <- c("<30", "30-50", ">50")
AGE_LABELS3 <- c("<30", "30-60", ">60")
AGE_SPLINES_MODELS <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 19, 20)

ADD_NBA_MODELS <- c(10, 11)
VL_LOWER_LIMIT_SIM_DATA <- 3.5
VL_LOWER_LIMIT <- 2
VL_LOWER_LIMIT_QUANTIFICATION <- 3


# These infections have log10 viral loads above 9 after day 20 post peak viral load
# in model20_106. I inspected them individually and they either have a very early
# low positive viral load (probably false positive) or apparently more than one peak
# Note: they should be automatically filtered out by filter_unusual_infections
IDs_TO_EXCLUDE1 <- c(
  "0a90d05f102ea057dd7caca8d56ae153", "0f8f7f329afda19e4eb14ff9b3a90a93",
  "13d3a33fc9d96209b90324c1c83ca02b", "1c550d8e27c663db8d2e72700634fe85",
  "2b43a5adc29ef934442d7b5d9ef99496", "3c0089662460668a0ebc7c83951ee223",
  "5b1a49553dbdd9072c19030fd258d038", "5d549fd98b0f855c472011ad84022f9c",
  "75ca756ce3d2eb3ea8aee0734f297198", "f82405211cdcad81065fee39730d5d18"
)

# The following infections aren't filtered out by filter_unusual_infections but
# nevertheless look weird.
IDs_TO_EXCLUDE2 <- c(
  "21f56be8d8128acd4aed0d7adb1e57c3", "5b99e78999236a901d4467b08030ea49",
  "75ca756ce3d2eb3ea8aee0734f297198", "890730c36cbdcf6e432b489b7303b694",
  "bbd91d52a5db4635dcb7e6df9365e9f4", "d7896f402d67d4df32beae3099229f23",
  "e9a167d1d1d7392e140be973b6b43a6d"
)

# These are infections that have at least one data point with a viral load below 4 placed
# between day 0 and 5; they look unusual upon visual inspection.
IDs_TO_EXCLUDE3 <- c(
  "113643bca58e750129b9f5798e2e7acc", "155e02234ccbbf8d04b60df3e988f813",
  "195e3aaad7b66766a7df46644f4c4249", "1b45f936336a220490d7742003a7267b",
  "1b6ab53b1c093f3679c718d5e1f6ddb3", "219617e04faeb6518fd8e3559c13a9fe",
  "2230724f5a5c21760f80e2717a974fd2", "22a6ed5a4e54c11c1a2201458713b7c8",
  "22f15a3b1e6506ef748f740bc3d40cc2", "23700ed1767546b676973421a12644f5",
  "2d87d59a420139d55bf976f5609e065a", "2e4a24ae6dc1fe18fa5a9c3bbdd30a79",
  "2edc1680418be07fbce7484259bc3173", "3ae9ec7edb5fd524ef6c17e2ac93a6a6",
  "3cfad2961b8c6796167c7844362a8d03", "4125a9c4e25c6cac7a105fe11a99b6bd",
  "43f14b30f43dd70f43d2c8e7f8418fca", "478118ea4207676905444c572a001349",
  "49359f93d6a3e739d02dd38f8ef535c7", "4a1f3292e769f013e49eb81ee9363115",
  "5a0a6faeb80a5ca03b46255f0202d29a", "5bec1224cf5904c2227c049f1bb6802f",
  "5c3946e7fc28610cc0dd4b74756fcbfc", "5e799cdd78370af9dc85d3cc782a5c41",
  "5fcbccef109d90e1f8ca96ff4c73f1a7", "65d5875f2abd98bf137c123c0547f369",
  "6c331fe61d6bce142523d46266f29d9d", "701018a428a1892944fb94bff6ee922f",
  "7147a010472001740d1b640048c60a06", "789537806a45e73b90679fe7acf3a136",
  "7c669598c5ff6aa556d07dae41d58d92", "7dc2f7928bcca67693eb5bb10f928613",
  "823a30b95281c7c8ee5009a1530d7cb7", "8e1eea6d4c79a0d799bb1b510785f3da",
  "8ef7fd5cf2a7d9ae06e4ea8f39eacf5f", "9002db26f95c365e8e2cbb8aec01e8e6",
  "a2cd9c99c343e0c8eddf22487c9ebb00", "b89ffec38e872d2d4feaaa190eaa1de0",
  "c7c3b4270f77dab798c5a43e4343c3d4", "cf56d21539cc394ec6fc8ba7f4a61f34",
  "d47d959f8e1c5588c9a943b21aa6a856", "d5de2dcfd6928692deebe85ee902026f",
  "d89e70a5a0570a06cc83681e12b41b73", "e681ec82c62ae51b7000c1c373f3d574",
  "e9d3790cc3261fced03ecfccb6cda215", "eb057a879241af55eb6720f9ab270b2f",
  "f504827db7d842dedc46b0df5eb53f38", "f5aa3cccaf1b3c4c7120c8e4ec105bb2"
)

IDs_TO_EXCLUDE <- c(IDs_TO_EXCLUDE1, IDs_TO_EXCLUDE2, IDs_TO_EXCLUDE3)


# Use T2 as the first test system to have it as the reference level in the statistical model
TEST_NAME_LEVELS <- c("T1", "T2", "LC480", "T3", "T6", "Ta", "T9", "T4", "T5", "Tb")
TEST_NAME_LEVELS_SIM_DATA <- c(1, 2, 3, 4, 5, 6)
GENDER_LEVELS <- c("F", "M")


SUBGRPS <- c(
  "PAMS", "gender", "hospitalized", "prior_infection", "age_cat2", "age_cat3",
  "alpha", "delta", "omicron", "unknown", "hospitalized:alpha",
  "hospitalized:delta", "hospitalized:omicron", "hospitalized:unknown"
)

SUBGRPS_2 <- c(
  "PAMS", "gender", "hospitalized", "prior_infection", "age_cat2", "age_cat3",
  "alpha", "delta", "omicron", "unknown", "hospitalized_alpha",
  "hospitalized_delta", "hospitalized_omicron", "hospitalized_unknown"
)

SUBGRPS_RHAT_LATEX <- list(
  mu = "$\\mu$",
  sigma = "$\\sigma$",
  alpha = "Alpha",
  delta = "Delta",
  omicron = "Omicron",
  PAMS = "PAMS",
  hospitalized = "Hospitalized",
  gender = "Male",
  age_cat2 = "30-60 years",
  age_cat3 = "\\textgreater{60} years",
  prior_infection = "Prior infection",
  hospitalized_alpha = "Hospitalized, Alpha",
  hospitalized_delta = "Hospitalized, Delta",
  hospitalized_omicron = "Hospitalized, Omicron",
  "hospitalized:alpha" = "Hospitalized, Alpha",
  "hospitalized:delta" = "Hospitalized, Delta",
  "hospitalized:omicron" = "Hospitalized, Omicron"
)


SUBGRPS_MAIN <- c(
  "PAMS3", "gender", "prior_infection", "age_category3", "variant"
)


SUBGRP_VARS_MAPPING <- list(
  PAMS1 = c("PAMS"),
  PAMS3 = c("hospitalized", "PAMS"),
  gender = c("gender"),
  hospitalized = c("hospitalized"),
  prior_infection = c("prior_infection"),
  age_category3 = c("age_cat2", "age_cat3"),
  variant = c("alpha", "delta", "omicron", "unknown")
)
SUBGRP_CATS_LABEL_MAPPING <- list(
  PAMS1 = c(FALSE, TRUE),
  PAMS3 = c(3, 2, 1), # Other, Hospitalized, PAMS
  gender = GENDER_LEVELS,
  hospitalized = c(FALSE, TRUE),
  prior_infection = c(FALSE, TRUE),
  age_category3 = AGE_LABELS3,
  variant = VARIANT_LEVELS
)
SUBGRP_CATS_LABEL_MAPPING2 <- list(
  PAMS1 = c("Non-PAMS", "PAMS"),
  PAMS3 = c("Non-PAMS, non-hosp.", "Hospitalized", "PAMS"), # Other, Hospitalized, PAMS
  gender = c("Female", "Male"),
  hospitalized = c("Not hospitalized", "Hospitalized"),
  prior_infection = c("First infection", "Prior infection"),
  age_category3 = AGE_LABELS3,
  variant = VARIANT_LEVELS
)


SUBGRP_MAPPING <- SUBGRPS

names(SUBGRP_MAPPING) <- 1:14
