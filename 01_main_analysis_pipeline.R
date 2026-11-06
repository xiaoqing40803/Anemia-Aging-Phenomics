


## =========================================================
## 01_main_analysis_pipeline.R
##
## Main reproducible analysis pipeline for:
## Anemia-Aging-Phenomics
##
## This script performs cohort construction, laboratory marker
## harmonization, marker-level age association analysis,
## hematometabolic domain score derivation, out-of-fold elastic-net
## modeling, outpatient Hb<120 sensitivity analysis, inpatient
## care-setting contrast analysis, and generation of analysis tables
## and source data.
##
## Figures are intentionally not generated in this main pipeline.
## =========================================================

rm(list = ls())

options(stringsAsFactors = FALSE)
options(na.print = "NA")
set.seed(20261103)

## =========================================================
## 0. Directories and global settings
## =========================================================

DATA_DIR <- Sys.getenv("ANEMIA_PHENOMICS_DATA_DIR", unset = "data_mock")
OUTDIR <- Sys.getenv("ANEMIA_PHENOMICS_OUTDIR", unset = "outputs_anemia_aging_phenomics")

DIR_OUTPATIENT <- file.path(OUTDIR, "01_outpatient_primary_HbIndexed")
DIR_OUTPATIENT_HB120 <- file.path(OUTDIR, "02_outpatient_sensitivity_Hb120")
DIR_INPATIENT <- file.path(OUTDIR, "03_inpatient_contrast_HbIndexed")
DIR_COMPARISON <- file.path(OUTDIR, "04_care_setting_contrast")
DIR_TABLES <- file.path(OUTDIR, "Main_Tables")
DIR_SOURCE <- file.path(OUTDIR, "SourceData")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OUTPATIENT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OUTPATIENT_HB120, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_INPATIENT, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_COMPARISON, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_SOURCE, recursive = TRUE, showWarnings = FALSE)

cat("\n=========================================================\n")
cat("Anemia-Aging-Phenomics main analysis pipeline started.\n")
cat("DATA_DIR:", DATA_DIR, "\n")
cat("OUTDIR:", OUTDIR, "\n")
cat("=========================================================\n")

## =========================================================
## 1. Packages
## =========================================================

required_pkgs <- c(
  "readr", "data.table",
  "dplyr", "tidyr", "purrr", "stringr", "lubridate",
  "ggplot2", "tibble", "forcats", "scales",
  "mgcv", "broom", "glmnet", "splines"
)

new_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]

if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, type = "binary", dependencies = TRUE)
}

invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

## =========================================================
## 2. Helper functions
## =========================================================

clean_id <- function(x) {
  stringr::str_squish(as.character(x))
}

standardize_sex <- function(x) {
  dplyr::case_when(
    x %in% c("男", "M", "Male", "male", "1", 1) ~ 1,
    x %in% c("女", "F", "Female", "female", "0", 0) ~ 0,
    TRUE ~ NA_real_
  )
}

parse_mixed_datetime <- function(x) {
  if (inherits(x, "POSIXct") | inherits(x, "POSIXt")) return(x)
  if (inherits(x, "Date")) return(as.POSIXct(x, tz = "Asia/Shanghai"))
  
  x_chr <- as.character(x)
  
  out1 <- suppressWarnings(lubridate::parse_date_time(
    x_chr,
    orders = c(
      "ymd HMS", "ymd HM", "ymd",
      "Ymd HMS", "Ymd HM", "Ymd",
      "mdy HMS", "mdy HM", "mdy",
      "dmy HMS", "dmy HM", "dmy"
    ),
    tz = "Asia/Shanghai"
  ))
  
  x_num <- suppressWarnings(as.numeric(x_chr))
  out2 <- suppressWarnings(as.POSIXct(
    (x_num - 25569) * 86400,
    origin = "1970-01-01",
    tz = "Asia/Shanghai"
  ))
  
  dplyr::if_else(!is.na(out1), out1, out2)
}

parse_numeric_result <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

format_est_ci <- function(est, low, high, digits = 2) {
  ifelse(
    is.na(est) | is.na(low) | is.na(high),
    NA_character_,
    paste0(
      sprintf(paste0("%.", digits, "f"), est),
      " (",
      sprintf(paste0("%.", digits, "f"), low),
      ", ",
      sprintf(paste0("%.", digits, "f"), high),
      ")"
    )
  )
}

get_col_or_na <- function(data, colname) {
  if (colname %in% names(data)) {
    data[[colname]]
  } else {
    rep(NA_real_, nrow(data))
  }
}

robust_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, constant = 1.4826, na.rm = TRUE)
  if (is.na(mad_val) || mad_val == 0) {
    mad_val <- sd(x, na.rm = TRUE)
  }
  if (is.na(mad_val) || mad_val == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - med) / mad_val
}

cap_z <- function(x, lower = -5, upper = 5) {
  pmin(pmax(x, lower), upper)
}

robust_z_by_sex <- function(x, sex) {
  out <- rep(NA_real_, length(x))
  for (sx in unique(sex[!is.na(sex)])) {
    idx <- which(sex == sx)
    med <- median(x[idx], na.rm = TRUE)
    mad_val <- mad(x[idx], constant = 1.4826, na.rm = TRUE)
    if (is.na(mad_val) || mad_val == 0) {
      mad_val <- sd(x[idx], na.rm = TRUE)
    }
    if (!is.na(mad_val) && mad_val > 0) {
      out[idx] <- (x[idx] - med) / mad_val
    }
  }
  out
}

fmt_n_pct <- function(n, pct, digits = 1) {
  ifelse(
    is.na(n) | is.na(pct),
    NA_character_,
    paste0(format(n, big.mark = ","), " (", sprintf(paste0("%.", digits, "f"), pct), "%)")
  )
}

fmt_median_iqr <- function(med, q1, q3, digits = 1) {
  ifelse(
    is.na(med) | is.na(q1) | is.na(q3),
    NA_character_,
    paste0(
      sprintf(paste0("%.", digits, "f"), med),
      " (",
      sprintf(paste0("%.", digits, "f"), q1),
      "–",
      sprintf(paste0("%.", digits, "f"), q3),
      ")"
    )
  )
}

fmt_int <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    format(round(x), big.mark = ",")
  )
}

## =========================================================
## 3. Laboratory item harmonization
## =========================================================

winsor_lab_range <- function(value, lab_std) {
  dplyr::case_when(
    lab_std == "hb" & (value < 30 | value > 250) ~ NA_real_,
    lab_std == "hct" & (value < 0.05 | value > 0.75) ~ NA_real_,
    lab_std == "rbc" & (value < 0.5 | value > 10) ~ NA_real_,
    lab_std == "mcv" & (value < 40 | value > 150) ~ NA_real_,
    lab_std == "mch" & (value < 10 | value > 60) ~ NA_real_,
    lab_std == "mchc" & (value < 200 | value > 450) ~ NA_real_,
    lab_std == "rdw" & (value < 5 | value > 40) ~ NA_real_,
    lab_std == "wbc" & (value < 0.1 | value > 300) ~ NA_real_,
    lab_std == "neut" & (value < 0 | value > 300) ~ NA_real_,
    lab_std == "lymph" & (value < 0 | value > 300) ~ NA_real_,
    lab_std == "nlr" & (value < 0 | value > 500) ~ NA_real_,
    lab_std == "plt" & (value < 5 | value > 2000) ~ NA_real_,
    lab_std == "ferritin" & (value < 1 | value > 100000) ~ NA_real_,
    lab_std == "tsat" & (value < 0 | value > 100) ~ NA_real_,
    lab_std == "serum_iron" & (value < 0 | value > 200) ~ NA_real_,
    lab_std == "tibc" & (value < 1 | value > 200) ~ NA_real_,
    lab_std == "uibc" & (value < 0 | value > 200) ~ NA_real_,
    lab_std == "transferrin" & (value < 0 | value > 10) ~ NA_real_,
    lab_std == "stfr" & (value < 0 | value > 100) ~ NA_real_,
    lab_std == "retic_count" & (value < 0 | value > 1000) ~ NA_real_,
    lab_std == "retic_pct" & (value < 0 | value > 50) ~ NA_real_,
    lab_std == "hscrp" & (value < 0 | value > 500) ~ NA_real_,
    lab_std == "crp" & (value < 0 | value > 500) ~ NA_real_,
    lab_std == "il6" & (value < 0 | value > 10000) ~ NA_real_,
    lab_std == "saa" & (value < 0 | value > 1000) ~ NA_real_,
    lab_std == "pct_inflam" & (value < 0 | value > 1000) ~ NA_real_,
    lab_std == "esr" & (value < 0 | value > 200) ~ NA_real_,
    lab_std == "alb" & (value < 5 | value > 80) ~ NA_real_,
    lab_std == "prealb" & (value < 10 | value > 1000) ~ NA_real_,
    lab_std == "tp" & (value < 10 | value > 120) ~ NA_real_,
    lab_std == "creatinine" & (value < 10 | value > 3000) ~ NA_real_,
    lab_std == "urea" & (value < 0 | value > 100) ~ NA_real_,
    lab_std == "cysc" & (value < 0 | value > 20) ~ NA_real_,
    lab_std == "egfr" & (value < 0 | value > 300) ~ NA_real_,
    lab_std == "uric_acid" & (value < 0 | value > 2000) ~ NA_real_,
    TRUE ~ value
  )
}

standardize_lab_item <- function(x) {
  x <- stringr::str_squish(as.character(x))
  
  dplyr::case_when(
    x %in% c("血红蛋白", "血色素", "HGB", "Hb", "Hemoglobin", "hemoglobin") ~ "hb",
    x == "红细胞压积" ~ "hct",
    x == "红细胞计数" ~ "rbc",
    x == "平均红细胞体积" ~ "mcv",
    x == "平均红细胞血红蛋白量" ~ "mch",
    x == "平均红细胞血红蛋白浓度" ~ "mchc",
    x == "红细胞分布宽度" ~ "rdw",
    x == "红细胞宽度变异系数" ~ "rdw_cv",
    x == "红细胞宽度标准差" ~ "rdw_sd",
    x == "白细胞计数" ~ "wbc",
    x == "中性粒细胞绝对值" ~ "neut",
    x == "淋巴细胞绝对值" ~ "lymph",
    x == "粒淋比" ~ "nlr",
    x %in% c("血小板计数", "血小板") ~ "plt",
    x %in% c("铁蛋白", "铁蛋白(frp)") ~ "ferritin",
    x == "转铁蛋白饱和度" ~ "tsat",
    x == "铁离子" ~ "serum_iron",
    x == "总铁结合力" ~ "tibc",
    x == "不饱和铁结合力" ~ "uibc",
    x == "转铁蛋白" ~ "transferrin",
    x == "血清转铁蛋白受体(Trf)" ~ "stfr",
    x == "网织红细胞计数" ~ "retic_count",
    x == "网织红细胞百分比" ~ "retic_pct",
    x %in% c("高敏C反应蛋白", "超敏C反应蛋白测定") ~ "hscrp",
    x == "C反应蛋白" ~ "crp",
    x %in% c("白细胞介素-6", "白介素6") ~ "il6",
    x == "淀粉样蛋白A" ~ "saa",
    x %in% c("降钙素原", "降钙素原检测(荧光定量法)") ~ "pct_inflam",
    x == "红细胞沉降率" ~ "esr",
    x == "白蛋白" ~ "alb",
    x == "前白蛋白" ~ "prealb",
    x == "总蛋白" ~ "tp",
    x == "肌酐" ~ "creatinine",
    x == "尿素" ~ "urea",
    x == "胱抑素C" ~ "cysc",
    x == "肾小球滤过率" ~ "egfr",
    x == "尿酸" ~ "uric_acid",
    TRUE ~ NA_character_
  )
}

analysis_labs <- c(
  "hb", "hct", "rbc", "mcv", "mch", "mchc", "rdw", "rdw_cv", "rdw_sd",
  "wbc", "neut", "lymph", "nlr", "plt",
  "ferritin", "tsat", "serum_iron", "tibc", "uibc", "transferrin", "stfr",
  "retic_count", "retic_pct",
  "hscrp", "crp", "il6", "saa", "pct_inflam", "esr",
  "alb", "prealb", "tp",
  "creatinine", "urea", "cysc", "egfr", "uric_acid"
)

baseline_marker_labs <- analysis_labs

marker_labels <- c(
  baseline_hb_final = "Hemoglobin",
  baseline_rbc = "RBC",
  baseline_hct = "Hct",
  baseline_mcv = "MCV",
  baseline_mch = "MCH",
  baseline_mchc = "MCHC",
  baseline_rdw = "RDW",
  baseline_wbc = "WBC",
  baseline_neut = "Neutrophils",
  baseline_lymph = "Lymphocytes",
  baseline_nlr_final = "NLR",
  baseline_plt = "Platelets",
  baseline_ferritin = "Ferritin",
  baseline_tsat = "TSAT",
  baseline_serum_iron = "Serum iron",
  baseline_tibc = "TIBC",
  baseline_transferrin = "Transferrin",
  baseline_crp_hscrp = "CRP/hsCRP",
  baseline_saa = "SAA",
  baseline_esr = "ESR",
  baseline_alb = "Albumin",
  baseline_prealb = "Prealbumin",
  baseline_creatinine = "Creatinine",
  baseline_urea = "Urea",
  baseline_cysc = "Cystatin C",
  baseline_egfr = "eGFR",
  baseline_uric_acid = "Uric acid"
)

log_transform_markers <- c(
  "baseline_ferritin",
  "baseline_crp_hscrp",
  "baseline_saa",
  "baseline_creatinine",
  "baseline_urea",
  "baseline_cysc",
  "baseline_uric_acid",
  "baseline_nlr_final"
)

clock_features_full <- c(
  "baseline_hb_final",
  "baseline_mcv",
  "baseline_mch",
  "baseline_mchc",
  "baseline_rdw",
  "baseline_wbc",
  "baseline_plt",
  "baseline_nlr_final",
  "baseline_ferritin",
  "baseline_alb",
  "baseline_creatinine",
  "baseline_urea",
  "baseline_uric_acid",
  "baseline_cysc"
)

clock_features_nonrenal <- setdiff(
  clock_features_full,
  c("baseline_creatinine", "baseline_urea", "baseline_uric_acid", "baseline_cysc")
)

domain_score_vars <- c(
  "erythroid_disruption_score",
  "renal_metabolic_burden_score",
  "inflammatory_burden_score",
  "nutritional_depletion_score",
  "iron_erythroid_imbalance_score"
)

domain_score_vars_sexstd <- paste0(domain_score_vars, "_sexstd")

## =========================================================
## 4. Marker windows
## =========================================================

make_marker_window <- function(setting = c("outpatient", "inpatient")) {
  setting <- match.arg(setting)
  
  if (setting == "outpatient") {
    tibble(
      lab_std = baseline_marker_labs,
      before_days = case_when(
        lab_std %in% c(
          "ferritin", "tsat", "serum_iron", "tibc", "uibc", "transferrin", "stfr",
          "hscrp", "crp", "il6", "saa", "pct_inflam", "esr",
          "retic_count", "retic_pct"
        ) ~ 90,
        TRUE ~ 30
      ),
      after_days = 7
    )
  } else {
    tibble(
      lab_std = baseline_marker_labs,
      before_days = case_when(
        lab_std %in% c(
          "ferritin", "tsat", "serum_iron", "tibc", "uibc", "transferrin", "stfr",
          "hscrp", "crp", "il6", "saa", "pct_inflam", "esr",
          "retic_count", "retic_pct"
        ) ~ 7,
        TRUE ~ 3
      ),
      after_days = case_when(
        lab_std %in% c(
          "ferritin", "tsat", "serum_iron", "tibc", "uibc", "transferrin", "stfr",
          "hscrp", "crp", "il6", "saa", "pct_inflam", "esr",
          "retic_count", "retic_pct"
        ) ~ 3,
        TRUE ~ 1
      )
    )
  }
}

## =========================================================
## 5. Reading and cohort construction
## =========================================================

read_csv_base <- function(file) {
  read.table(
    file = file,
    sep = ",",
    header = TRUE,
    stringsAsFactors = FALSE,
    skip = 0,
    na.strings = "",
    fileEncoding = "UTF-8"
  )
}

build_lab_long_outpatient <- function(file1, file2) {
  dat1 <- read_csv_base(file1)
  dat2 <- read_csv_base(file2)
  
  bind_rows(dat1, dat2) %>%
    transmute(
      PID = clean_id(pid),
      SEX = standardize_sex(sex),
      AGE = as.numeric(age),
      BARCODE = as.character(条码号),
      lab_item_raw = as.character(检测项目),
      lab_value_raw = as.character(结果),
      lab_unit = as.character(单位),
      lab_range = as.character(检验参考值范围),
      abnormal_flag = as.character(异常提示),
      lab_time = parse_mixed_datetime(报告时间)
    ) %>%
    mutate(
      lab_item = stringr::str_squish(lab_item_raw),
      lab_std = standardize_lab_item(lab_item),
      lab_value_num_raw = parse_numeric_result(lab_value_raw),
      lab_value_num = case_when(
        lab_std == "hb" & lab_unit %in% c("g/dL", "g/dl", "G/dL") ~ lab_value_num_raw * 10,
        lab_std == "hb" & is.na(lab_unit) & lab_value_num_raw <= 25 ~ lab_value_num_raw * 10,
        TRUE ~ lab_value_num_raw
      ),
      lab_value_num = winsor_lab_range(lab_value_num, lab_std)
    ) %>%
    filter(
      !is.na(PID), PID != "",
      !is.na(lab_time),
      lab_std %in% analysis_labs,
      !is.na(lab_value_num)
    ) %>%
    group_by(PID, SEX, AGE, lab_time, lab_std) %>%
    summarise(
      lab_value_num = mean(lab_value_num, na.rm = TRUE),
      lab_unit = dplyr::first(lab_unit),
      .groups = "drop"
    )
}

build_lab_long_inpatient <- function(file1, file2) {
  dat1 <- read_csv_base(file1)
  dat2 <- read_csv_base(file2)
  dat <- bind_rows(dat1, dat2)
  
  admission_col <- if ("ADMINISSION_NO" %in% names(dat)) {
    "ADMINISSION_NO"
  } else if ("ADMISSION_NO" %in% names(dat)) {
    "ADMISSION_NO"
  } else {
    stop("Cannot find ADMINISSION_NO or ADMISSION_NO column in inpatient files.")
  }
  
  dat %>%
    transmute(
      ADMISSION_NO = clean_id(.data[[admission_col]]),
      PID = clean_id(pid),
      PATIENT_NAME = if ("PATIENT_NAME" %in% names(.)) as.character(PATIENT_NAME) else NA_character_,
      SEX = standardize_sex(sex),
      AGE = as.numeric(age),
      WEIGHT = if ("weight" %in% names(.)) as.numeric(weight) else NA_real_,
      HEIGHT = if ("hight" %in% names(.)) as.numeric(hight) else NA_real_,
      BARCODE = as.character(条码号),
      lab_item_raw = as.character(检测项目),
      lab_value_raw = as.character(结果),
      lab_unit = as.character(单位),
      lab_range = as.character(检验参考值范围),
      abnormal_flag = as.character(异常提示),
      lab_time = parse_mixed_datetime(报告时间)
    ) %>%
    mutate(
      lab_item = stringr::str_squish(lab_item_raw),
      lab_std = standardize_lab_item(lab_item),
      lab_value_num_raw = parse_numeric_result(lab_value_raw),
      lab_value_num = case_when(
        lab_std == "hb" & lab_unit %in% c("g/dL", "g/dl", "G/dL") ~ lab_value_num_raw * 10,
        lab_std == "hb" & is.na(lab_unit) & lab_value_num_raw <= 25 ~ lab_value_num_raw * 10,
        TRUE ~ lab_value_num_raw
      ),
      lab_value_num = winsor_lab_range(lab_value_num, lab_std)
    ) %>%
    filter(
      !is.na(ADMISSION_NO), ADMISSION_NO != "",
      !is.na(PID), PID != "",
      !is.na(lab_time),
      lab_std %in% analysis_labs,
      !is.na(lab_value_num)
    ) %>%
    group_by(
      ADMISSION_NO, PID, PATIENT_NAME, SEX, AGE, WEIGHT, HEIGHT,
      lab_time, lab_std
    ) %>%
    summarise(
      lab_value_num = mean(lab_value_num, na.rm = TRUE),
      lab_unit = dplyr::first(lab_unit),
      .groups = "drop"
    )
}

build_hb_index_outpatient <- function(lab_long) {
  hb_true <- lab_long %>%
    filter(lab_std == "hb") %>%
    transmute(
      PID,
      hb_time = lab_time,
      hb_value = lab_value_num,
      AGE,
      SEX
    ) %>%
    arrange(PID, hb_time) %>%
    distinct(PID, hb_time, hb_value, .keep_all = TRUE)
  
  hb_true %>%
    group_by(PID) %>%
    arrange(hb_time, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      PID,
      index_time = hb_time,
      baseline_hb = hb_value,
      baseline_age = AGE,
      baseline_sex = SEX
    ) %>%
    filter(
      !is.na(baseline_age),
      !is.na(baseline_sex),
      !is.na(baseline_hb),
      baseline_age >= 18,
      baseline_age <= 105
    ) %>%
    mutate(
      sex_label = case_when(
        baseline_sex == 1 ~ "Male",
        baseline_sex == 0 ~ "Female",
        TRUE ~ NA_character_
      ),
      age_group_4 = cut(
        baseline_age,
        breaks = c(18, 45, 60, 75, Inf),
        labels = c("18-44", "45-59", "60-74", ">=75"),
        right = FALSE
      ),
      anemia_hb120 = baseline_hb < 120,
      anemia_sex_specific_hb_indexed = case_when(
        baseline_sex == 1 & baseline_hb < 130 ~ TRUE,
        baseline_sex == 0 & baseline_hb < 120 ~ TRUE,
        TRUE ~ FALSE
      )
    )
}

build_hb_index_inpatient <- function(lab_long) {
  hb_true <- lab_long %>%
    filter(lab_std == "hb") %>%
    transmute(
      ADMISSION_NO,
      PID,
      hb_time = lab_time,
      hb_value = lab_value_num,
      AGE,
      SEX
    ) %>%
    arrange(ADMISSION_NO, hb_time) %>%
    distinct(ADMISSION_NO, hb_time, hb_value, .keep_all = TRUE)
  
  hb_true %>%
    group_by(ADMISSION_NO) %>%
    arrange(hb_time, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(
      ADMISSION_NO,
      PID,
      index_time = hb_time,
      baseline_hb = hb_value,
      baseline_age = AGE,
      baseline_sex = SEX
    ) %>%
    filter(
      !is.na(baseline_age),
      !is.na(baseline_sex),
      !is.na(baseline_hb),
      baseline_age >= 18,
      baseline_age <= 105
    ) %>%
    mutate(
      sex_label = case_when(
        baseline_sex == 1 ~ "Male",
        baseline_sex == 0 ~ "Female",
        TRUE ~ NA_character_
      ),
      age_group_4 = cut(
        baseline_age,
        breaks = c(18, 45, 60, 75, Inf),
        labels = c("18-44", "45-59", "60-74", ">=75"),
        right = FALSE
      ),
      anemia_hb120 = baseline_hb < 120,
      anemia_sex_specific_hb_indexed = case_when(
        baseline_sex == 1 & baseline_hb < 130 ~ TRUE,
        baseline_sex == 0 & baseline_hb < 120 ~ TRUE,
        TRUE ~ FALSE
      )
    )
}

build_baseline_wide <- function(lab_long, index_data, setting = c("outpatient", "inpatient")) {
  setting <- match.arg(setting)
  marker_window <- make_marker_window(setting)
  
  if (setting == "outpatient") {
    baseline_long <- lab_long %>%
      filter(lab_std %in% baseline_marker_labs) %>%
      inner_join(index_data %>% select(PID, index_time), by = "PID") %>%
      inner_join(marker_window, by = "lab_std") %>%
      mutate(
        days_from_index = as.numeric(difftime(lab_time, index_time, units = "days")),
        abs_days = abs(days_from_index)
      ) %>%
      filter(days_from_index >= -before_days, days_from_index <= after_days) %>%
      group_by(PID, lab_std) %>%
      slice_min(abs_days, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    value_wide <- baseline_long %>%
      select(PID, lab_std, lab_value_num) %>%
      pivot_wider(
        id_cols = PID,
        names_from = lab_std,
        values_from = lab_value_num,
        names_prefix = "baseline_"
      )
    
    day_wide <- baseline_long %>%
      select(PID, lab_std, days_from_index) %>%
      pivot_wider(
        id_cols = PID,
        names_from = lab_std,
        values_from = days_from_index,
        names_glue = "baseline_{lab_std}_day"
      )
    
    baseline_wide <- value_wide %>%
      left_join(day_wide, by = "PID")
    
    analysis <- index_data %>%
      rename(
        index_hb = baseline_hb,
        index_age = baseline_age,
        index_sex = baseline_sex
      ) %>%
      left_join(baseline_wide, by = "PID") %>%
      mutate(
        baseline_age = index_age,
        baseline_sex = index_sex,
        baseline_hb_final = coalesce(get_col_or_na(., "baseline_hb"), index_hb),
        baseline_crp_hscrp = coalesce(get_col_or_na(., "baseline_hscrp"), get_col_or_na(., "baseline_crp")),
        baseline_nlr_final = case_when(
          !is.na(get_col_or_na(., "baseline_nlr")) ~ get_col_or_na(., "baseline_nlr"),
          !is.na(get_col_or_na(., "baseline_neut")) &
            !is.na(get_col_or_na(., "baseline_lymph")) &
            get_col_or_na(., "baseline_lymph") > 0 ~
            get_col_or_na(., "baseline_neut") / get_col_or_na(., "baseline_lymph"),
          TRUE ~ NA_real_
        )
      )
  } else {
    baseline_long <- lab_long %>%
      filter(lab_std %in% baseline_marker_labs) %>%
      inner_join(index_data %>% select(ADMISSION_NO, PID, index_time), by = c("ADMISSION_NO", "PID")) %>%
      inner_join(marker_window, by = "lab_std") %>%
      mutate(
        days_from_index = as.numeric(difftime(lab_time, index_time, units = "days")),
        abs_days = abs(days_from_index)
      ) %>%
      filter(days_from_index >= -before_days, days_from_index <= after_days) %>%
      group_by(ADMISSION_NO, lab_std) %>%
      slice_min(abs_days, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    baseline_wide <- baseline_long %>%
      filter(lab_std != "hb") %>%
      select(ADMISSION_NO, lab_std, lab_value_num, days_from_index) %>%
      pivot_wider(
        id_cols = ADMISSION_NO,
        names_from = lab_std,
        values_from = c(lab_value_num, days_from_index),
        names_glue = "baseline_{lab_std}{ifelse(.value == 'days_from_index', '_day', '')}"
      )
    
    analysis <- index_data %>%
      left_join(baseline_wide, by = "ADMISSION_NO") %>%
      mutate(
        baseline_hb_final = baseline_hb,
        baseline_crp_hscrp = coalesce(get_col_or_na(., "baseline_hscrp"), get_col_or_na(., "baseline_crp")),
        baseline_nlr_final = case_when(
          !is.na(get_col_or_na(., "baseline_nlr")) ~ get_col_or_na(., "baseline_nlr"),
          !is.na(get_col_or_na(., "baseline_neut")) &
            !is.na(get_col_or_na(., "baseline_lymph")) &
            get_col_or_na(., "baseline_lymph") > 0 ~
            get_col_or_na(., "baseline_neut") / get_col_or_na(., "baseline_lymph"),
          TRUE ~ NA_real_
        )
      )
  }
  
  list(
    baseline_long = baseline_long,
    baseline_wide = baseline_wide,
    analysis = analysis
  )
}

## =========================================================
## 6. Analysis functions
## =========================================================

summarise_cohort <- function(data) {
  data %>%
    summarise(
      n_units = n(),
      n_patients = n_distinct(PID),
      age_mean = mean(baseline_age, na.rm = TRUE),
      age_sd = sd(baseline_age, na.rm = TRUE),
      age_median = median(baseline_age, na.rm = TRUE),
      age_q1 = quantile(baseline_age, 0.25, na.rm = TRUE),
      age_q3 = quantile(baseline_age, 0.75, na.rm = TRUE),
      age_min = min(baseline_age, na.rm = TRUE),
      age_max = max(baseline_age, na.rm = TRUE),
      female_n = sum(baseline_sex == 0, na.rm = TRUE),
      female_pct = mean(baseline_sex == 0, na.rm = TRUE) * 100,
      male_n = sum(baseline_sex == 1, na.rm = TRUE),
      male_pct = mean(baseline_sex == 1, na.rm = TRUE) * 100,
      hb_mean = mean(baseline_hb_final, na.rm = TRUE),
      hb_sd = sd(baseline_hb_final, na.rm = TRUE),
      hb_median = median(baseline_hb_final, na.rm = TRUE),
      hb_q1 = quantile(baseline_hb_final, 0.25, na.rm = TRUE),
      hb_q3 = quantile(baseline_hb_final, 0.75, na.rm = TRUE)
    )
}

summarise_age_sex <- function(data) {
  data %>%
    filter(!is.na(age_group_4), !is.na(sex_label)) %>%
    count(age_group_4, sex_label, name = "n") %>%
    group_by(age_group_4) %>%
    mutate(
      age_group_total = sum(n),
      pct_within_age_group = n / age_group_total * 100
    ) %>%
    ungroup()
}

marker_coverage <- function(data) {
  vars <- c(
    "baseline_hb_final",
    "baseline_rbc", "baseline_hct", "baseline_mcv", "baseline_mch", "baseline_mchc", "baseline_rdw",
    "baseline_wbc", "baseline_neut", "baseline_lymph", "baseline_nlr_final", "baseline_plt",
    "baseline_ferritin", "baseline_tsat", "baseline_serum_iron", "baseline_tibc", "baseline_transferrin",
    "baseline_retic_count", "baseline_retic_pct",
    "baseline_crp_hscrp", "baseline_saa", "baseline_esr",
    "baseline_alb", "baseline_prealb",
    "baseline_creatinine", "baseline_urea", "baseline_cysc", "baseline_egfr", "baseline_uric_acid"
  )
  
  vars <- vars[vars %in% names(data)]
  
  tibble(
    marker = vars,
    marker_label = dplyr::recode(vars, !!!marker_labels, .default = vars),
    n_available = sapply(vars, function(v) sum(!is.na(data[[v]]))),
    pct_available = sapply(vars, function(v) mean(!is.na(data[[v]])) * 100),
    median_value = sapply(vars, function(v) median(data[[v]], na.rm = TRUE)),
    q1 = sapply(vars, function(v) quantile(data[[v]], 0.25, na.rm = TRUE)),
    q3 = sapply(vars, function(v) quantile(data[[v]], 0.75, na.rm = TRUE))
  ) %>%
    arrange(desc(pct_available))
}

fit_marker_gam_table <- function(data) {
  vars <- c(
    "baseline_hb_final",
    "baseline_rbc", "baseline_hct", "baseline_mcv", "baseline_mch", "baseline_mchc", "baseline_rdw",
    "baseline_wbc", "baseline_neut", "baseline_lymph", "baseline_nlr_final", "baseline_plt",
    "baseline_ferritin", "baseline_crp_hscrp", "baseline_saa", "baseline_esr",
    "baseline_alb", "baseline_prealb",
    "baseline_creatinine", "baseline_urea", "baseline_cysc", "baseline_egfr", "baseline_uric_acid"
  )
  
  vars <- vars[vars %in% names(data)]
  
  purrr::map_dfr(vars, function(marker_var) {
    d <- data %>%
      select(
        baseline_age,
        baseline_sex,
        value = all_of(marker_var)
      ) %>%
      filter(
        !is.na(value),
        !is.na(baseline_age),
        !is.na(baseline_sex)
      )
    
    n_available <- nrow(d)
    pct_available <- n_available / nrow(data) * 100
    
    if (n_available < 200 || pct_available < 5) {
      return(tibble(
        marker = marker_var,
        marker_label = dplyr::recode(marker_var, !!!marker_labels, .default = marker_var),
        n = n_available,
        pct_available = pct_available,
        spearman_rho = NA_real_,
        spearman_p = NA_real_,
        gam_p_age = NA_real_,
        direction_18_80 = NA_character_,
        pred_18 = NA_real_,
        pred_80 = NA_real_
      ))
    }
    
    d <- d %>%
      mutate(value_model = if (marker_var %in% log_transform_markers) log1p(value) else value)
    
    sp <- suppressWarnings(cor.test(d$baseline_age, d$value_model, method = "spearman", exact = FALSE))
    
    fit <- mgcv::gam(
      value_model ~ s(baseline_age, k = 5) + baseline_sex,
      data = d,
      method = "REML"
    )
    
    gam_p <- summary(fit)$s.table[1, "p-value"]
    
    age_grid <- tibble(
      baseline_age = seq(
        max(18, floor(min(d$baseline_age, na.rm = TRUE))),
        min(95, ceiling(max(d$baseline_age, na.rm = TRUE))),
        by = 1
      ),
      baseline_sex = median(d$baseline_sex, na.rm = TRUE)
    )
    
    age_grid$pred <- as.numeric(predict(fit, newdata = age_grid, type = "response"))
    
    pred_18 <- age_grid %>%
      filter(abs(baseline_age - 18) == min(abs(baseline_age - 18))) %>%
      slice(1) %>%
      pull(pred)
    
    pred_80 <- age_grid %>%
      filter(abs(baseline_age - 80) == min(abs(baseline_age - 80))) %>%
      slice(1) %>%
      pull(pred)
    
    direction <- case_when(
      is.na(pred_18) | is.na(pred_80) ~ NA_character_,
      pred_80 > pred_18 ~ "Higher at older age",
      pred_80 < pred_18 ~ "Lower at older age",
      TRUE ~ "Flat"
    )
    
    tibble(
      marker = marker_var,
      marker_label = dplyr::recode(marker_var, !!!marker_labels, .default = marker_var),
      n = n_available,
      pct_available = pct_available,
      spearman_rho = unname(sp$estimate),
      spearman_p = sp$p.value,
      gam_p_age = gam_p,
      direction_18_80 = direction,
      pred_18 = pred_18,
      pred_80 = pred_80
    )
  }) %>%
    mutate(
      gam_p_fdr = p.adjust(gam_p_age, method = "BH"),
      spearman_p_fdr = p.adjust(spearman_p, method = "BH"),
      gam_p_fdr_label = format_p(gam_p_fdr),
      spearman_p_fdr_label = format_p(spearman_p_fdr)
    ) %>%
    arrange(gam_p_fdr)
}

add_domain_scores <- function(data) {
  d <- data %>%
    mutate(
      log_ferritin = if ("baseline_ferritin" %in% names(.)) log1p(baseline_ferritin) else NA_real_,
      log_hscrp = if ("baseline_crp_hscrp" %in% names(.)) log1p(baseline_crp_hscrp) else NA_real_,
      log_nlr = if ("baseline_nlr_final" %in% names(.)) log1p(baseline_nlr_final) else NA_real_,
      log_creatinine = if ("baseline_creatinine" %in% names(.)) log1p(baseline_creatinine) else NA_real_,
      log_urea = if ("baseline_urea" %in% names(.)) log1p(baseline_urea) else NA_real_,
      log_cysc = if ("baseline_cysc" %in% names(.)) log1p(baseline_cysc) else NA_real_,
      log_uric_acid = if ("baseline_uric_acid" %in% names(.)) log1p(baseline_uric_acid) else NA_real_
    )
  
  score_vars <- c(
    "baseline_hb_final", "baseline_rbc", "baseline_hct",
    "baseline_mcv", "baseline_rdw",
    "baseline_wbc", "baseline_plt",
    "log_nlr", "log_ferritin", "log_hscrp",
    "baseline_alb", "baseline_prealb",
    "log_creatinine", "log_urea", "log_cysc", "log_uric_acid"
  )
  
  score_vars <- score_vars[score_vars %in% names(d)]
  
  for (v in score_vars) {
    d[[paste0(v, "_z")]] <- cap_z(robust_z(d[[v]]))
    d[[paste0(v, "_sexstd_z")]] <- cap_z(robust_z_by_sex(d[[v]], d$baseline_sex))
  }
  
  d <- d %>%
    mutate(
      erythroid_disruption_score = rowMeans(
        cbind(
          -baseline_hb_final_z,
          -baseline_rbc_z,
          -baseline_hct_z,
          baseline_rdw_z
        ),
        na.rm = TRUE
      ),
      renal_metabolic_burden_score = rowMeans(
        cbind(
          log_creatinine_z,
          log_urea_z,
          log_cysc_z,
          log_uric_acid_z
        ),
        na.rm = TRUE
      ),
      inflammatory_burden_score = rowMeans(
        cbind(
          log_hscrp_z,
          log_nlr_z,
          baseline_wbc_z,
          baseline_plt_z
        ),
        na.rm = TRUE
      ),
      nutritional_depletion_score = rowMeans(
        cbind(
          -baseline_alb_z,
          -baseline_prealb_z
        ),
        na.rm = TRUE
      ),
      iron_erythroid_imbalance_score = rowMeans(
        cbind(
          log_ferritin_z,
          -baseline_mcv_z,
          baseline_rdw_z
        ),
        na.rm = TRUE
      ),
      erythroid_disruption_score_sexstd = rowMeans(
        cbind(
          -baseline_hb_final_sexstd_z,
          -baseline_rbc_sexstd_z,
          -baseline_hct_sexstd_z,
          baseline_rdw_sexstd_z
        ),
        na.rm = TRUE
      ),
      renal_metabolic_burden_score_sexstd = rowMeans(
        cbind(
          log_creatinine_sexstd_z,
          log_urea_sexstd_z,
          log_cysc_sexstd_z,
          log_uric_acid_sexstd_z
        ),
        na.rm = TRUE
      ),
      inflammatory_burden_score_sexstd = rowMeans(
        cbind(
          log_hscrp_sexstd_z,
          log_nlr_sexstd_z,
          baseline_wbc_sexstd_z,
          baseline_plt_sexstd_z
        ),
        na.rm = TRUE
      ),
      nutritional_depletion_score_sexstd = rowMeans(
        cbind(
          -baseline_alb_sexstd_z,
          -baseline_prealb_sexstd_z
        ),
        na.rm = TRUE
      ),
      iron_erythroid_imbalance_score_sexstd = rowMeans(
        cbind(
          log_ferritin_sexstd_z,
          -baseline_mcv_sexstd_z,
          baseline_rdw_sexstd_z
        ),
        na.rm = TRUE
      )
    ) %>%
    mutate(
      across(
        all_of(c(domain_score_vars, domain_score_vars_sexstd)),
        ~ ifelse(is.nan(.x), NA_real_, .x)
      )
    )
  
  d
}

summarise_domain_scores <- function(data) {
  data %>%
    summarise(
      across(
        all_of(c(domain_score_vars, domain_score_vars_sexstd)),
        list(
          n = ~sum(!is.na(.x)),
          mean = ~mean(.x, na.rm = TRUE),
          sd = ~sd(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          q1 = ~quantile(.x, 0.25, na.rm = TRUE),
          q3 = ~quantile(.x, 0.75, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    )
}

domain_by_age_sex <- function(data) {
  data %>%
    select(
      baseline_age,
      age_group_4,
      baseline_sex,
      sex_label,
      all_of(domain_score_vars),
      all_of(domain_score_vars_sexstd)
    ) %>%
    pivot_longer(
      cols = all_of(c(domain_score_vars, domain_score_vars_sexstd)),
      names_to = "domain_score",
      values_to = "score"
    ) %>%
    mutate(
      score_type = if_else(str_detect(domain_score, "_sexstd$"), "Sex-standardized", "Whole-cohort standardized"),
      domain_label = case_when(
        str_detect(domain_score, "erythroid_disruption") ~ "Erythroid disruption",
        str_detect(domain_score, "renal_metabolic_burden") ~ "Renal-metabolic burden",
        str_detect(domain_score, "inflammatory_burden") ~ "Inflammatory burden",
        str_detect(domain_score, "nutritional_depletion") ~ "Nutritional depletion",
        str_detect(domain_score, "iron_erythroid_imbalance") ~ "Iron-erythroid imbalance",
        TRUE ~ domain_score
      )
    ) %>%
    filter(!is.na(score), !is.na(age_group_4), !is.na(sex_label)) %>%
    group_by(score_type, domain_label, age_group_4, sex_label) %>%
    summarise(
      n = n(),
      mean_score = mean(score, na.rm = TRUE),
      sd_score = sd(score, na.rm = TRUE),
      se_score = sd_score / sqrt(n),
      median_score = median(score, na.rm = TRUE),
      q1 = quantile(score, 0.25, na.rm = TRUE),
      q3 = quantile(score, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

fit_oof_model <- function(
    data,
    feature_set,
    model_name,
    id_cols = c("PID"),
    nfolds = 5,
    seed = 20261103
) {
  set.seed(seed)
  
  feature_set <- feature_set[feature_set %in% names(data)]
  
  d <- data %>%
    select(
      any_of(id_cols),
      baseline_age,
      baseline_sex,
      sex_label,
      age_group_4,
      all_of(domain_score_vars_sexstd),
      all_of(feature_set)
    ) %>%
    mutate(
      n_features_available = rowSums(!is.na(select(., all_of(feature_set))))
    ) %>%
    filter(
      !is.na(baseline_age),
      !is.na(baseline_sex),
      baseline_age >= 18,
      baseline_age <= 105,
      n_features_available >= max(4, ceiling(length(feature_set) * 0.5))
    )
  
  cat("\nRunning ", model_name, " with ", nrow(d), " units and ", length(feature_set), " features...\n", sep = "")
  
  if (nrow(d) < 500) {
    warning("Too few units for OOF model: ", model_name)
    return(NULL)
  }
  
  d <- d %>%
    mutate(fold_strata = interaction(age_group_4, baseline_sex, drop = TRUE))
  
  d$fold_id <- NA_integer_
  
  for (s in unique(d$fold_strata)) {
    idx <- which(d$fold_strata == s)
    d$fold_id[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
  }
  
  pred <- rep(NA_real_, nrow(d))
  
  for (fold in seq_len(nfolds)) {
    train_idx <- which(d$fold_id != fold)
    test_idx <- which(d$fold_id == fold)
    
    train_x <- d[train_idx, feature_set, drop = FALSE]
    test_x <- d[test_idx, feature_set, drop = FALSE]
    
    log_vars <- intersect(log_transform_markers, feature_set)
    
    for (v in log_vars) {
      train_x[[v]] <- log1p(train_x[[v]])
      test_x[[v]] <- log1p(test_x[[v]])
    }
    
    train_medians <- sapply(train_x, function(x) median(x, na.rm = TRUE))
    
    for (v in feature_set) {
      train_x[[v]][is.na(train_x[[v]])] <- train_medians[[v]]
      test_x[[v]][is.na(test_x[[v]])] <- train_medians[[v]]
    }
    
    train_means <- sapply(train_x, mean, na.rm = TRUE)
    train_sds <- sapply(train_x, sd, na.rm = TRUE)
    train_sds[is.na(train_sds) | train_sds == 0] <- 1
    
    train_x_scaled <- sweep(as.matrix(train_x), 2, train_means, "-")
    train_x_scaled <- sweep(train_x_scaled, 2, train_sds, "/")
    test_x_scaled <- sweep(as.matrix(test_x), 2, train_means, "-")
    test_x_scaled <- sweep(test_x_scaled, 2, train_sds, "/")
    
    set.seed(seed + fold)
    
    cv_fit <- glmnet::cv.glmnet(
      x = train_x_scaled,
      y = d$baseline_age[train_idx],
      alpha = 0.5,
      family = "gaussian",
      nfolds = 3,
      standardize = FALSE,
      nlambda = 50
    )
    
    pred[test_idx] <- as.numeric(
      predict(cv_fit, newx = test_x_scaled, s = "lambda.min")
    )
  }
  
  d <- d %>%
    mutate(
      model_name = model_name,
      predicted_age_oof = pred,
      raw_age_gap_oof = predicted_age_oof - baseline_age
    )
  
  calibration <- lm(predicted_age_oof ~ baseline_age, data = d)
  
  d <- d %>%
    mutate(age_gap_oof = residuals(calibration))
  
  performance <- tibble(
    model_name = model_name,
    n = nrow(d),
    n_features = length(feature_set),
    correlation_predicted_vs_chronological = cor(d$predicted_age_oof, d$baseline_age, use = "complete.obs"),
    MAE = mean(abs(d$predicted_age_oof - d$baseline_age), na.rm = TRUE),
    RMSE = sqrt(mean((d$predicted_age_oof - d$baseline_age)^2, na.rm = TRUE)),
    age_gap_sd = sd(d$age_gap_oof, na.rm = TRUE)
  )
  
  domain_associations <- purrr::map_dfr(domain_score_vars_sexstd, function(score_var) {
    d_model <- d %>%
      mutate(score_value = .data[[score_var]]) %>%
      filter(
        !is.na(age_gap_oof),
        !is.na(score_value),
        !is.na(baseline_age),
        !is.na(baseline_sex)
      )
    
    if (nrow(d_model) < 200) return(tibble())
    
    fit <- lm(age_gap_oof ~ score_value + baseline_age + baseline_sex, data = d_model)
    
    broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == "score_value") %>%
      mutate(
        model_name = model_name,
        domain_score = score_var,
        domain_label = case_when(
          score_var == "erythroid_disruption_score_sexstd" ~ "Erythroid disruption",
          score_var == "renal_metabolic_burden_score_sexstd" ~ "Renal-metabolic burden",
          score_var == "inflammatory_burden_score_sexstd" ~ "Inflammatory burden",
          score_var == "nutritional_depletion_score_sexstd" ~ "Nutritional depletion",
          score_var == "iron_erythroid_imbalance_score_sexstd" ~ "Iron-erythroid imbalance",
          TRUE ~ score_var
        ),
        n = nrow(d_model),
        p_value_label = format_p(p.value),
        estimate_ci = format_est_ci(estimate, conf.low, conf.high, digits = 2)
      )
  })
  
  list(
    performance = performance,
    domain_associations = domain_associations,
    data = d
  )
}

run_cohort_analysis <- function(data, outdir, prefix, id_cols = c("PID")) {
  cat("\n=========================================================\n")
  cat("Running cohort analysis:", prefix, "\n")
  cat("=========================================================\n")
  
  cohort_summary <- summarise_cohort(data)
  age_sex_summary <- summarise_age_sex(data)
  coverage <- marker_coverage(data)
  gam_table <- fit_marker_gam_table(data)
  
  domain_data <- add_domain_scores(data)
  domain_summary <- summarise_domain_scores(domain_data)
  domain_age_sex <- domain_by_age_sex(domain_data)
  
  oof_results <- list(
    full = fit_oof_model(
      domain_data,
      clock_features_full,
      paste0(prefix, " full high-coverage feature set"),
      id_cols = id_cols
    ),
    nonrenal = fit_oof_model(
      domain_data,
      clock_features_nonrenal,
      paste0(prefix, " non-renal feature set"),
      id_cols = id_cols
    )
  )
  
  oof_results <- oof_results[!sapply(oof_results, is.null)]
  
  if (length(oof_results) > 0) {
    oof_perf <- bind_rows(lapply(oof_results, function(x) x$performance))
    oof_assoc <- bind_rows(lapply(oof_results, function(x) x$domain_associations))
    oof_data <- bind_rows(lapply(oof_results, function(x) x$data))
  } else {
    oof_perf <- tibble()
    oof_assoc <- tibble()
    oof_data <- tibble()
  }
  
  write.csv(cohort_summary, file.path(outdir, paste0(prefix, "_Table_cohort_summary.csv")), row.names = FALSE)
  write.csv(age_sex_summary, file.path(outdir, paste0(prefix, "_Table_age_group_sex_distribution.csv")), row.names = FALSE)
  write.csv(coverage, file.path(outdir, paste0(prefix, "_Table_marker_coverage.csv")), row.names = FALSE)
  write.csv(gam_table, file.path(outdir, paste0(prefix, "_Table_marker_age_GAM_associations.csv")), row.names = FALSE)
  write.csv(domain_summary, file.path(outdir, paste0(prefix, "_Table_domain_score_summary.csv")), row.names = FALSE)
  write.csv(domain_age_sex, file.path(outdir, paste0(prefix, "_Table_domain_scores_by_age_group_sex.csv")), row.names = FALSE)
  write.csv(domain_data, file.path(outdir, paste0(prefix, "_SourceData_domain_scores.csv")), row.names = FALSE)
  write.csv(oof_perf, file.path(outdir, paste0(prefix, "_Table_OOF_model_performance.csv")), row.names = FALSE)
  write.csv(oof_assoc, file.path(outdir, paste0(prefix, "_Table_OOF_domain_associations.csv")), row.names = FALSE)
  write.csv(oof_data, file.path(outdir, paste0(prefix, "_SourceData_OOF_model_data.csv")), row.names = FALSE)
  
  saveRDS(
    list(
      analysis_data = data,
      domain_data = domain_data,
      cohort_summary = cohort_summary,
      age_sex_summary = age_sex_summary,
      coverage = coverage,
      gam_table = gam_table,
      domain_summary = domain_summary,
      domain_age_sex = domain_age_sex,
      oof_performance = oof_perf,
      oof_domain_associations = oof_assoc,
      oof_data = oof_data
    ),
    file.path(outdir, paste0(prefix, "_analysis_object.rds"))
  )
  
  list(
    analysis_data = data,
    domain_data = domain_data,
    cohort_summary = cohort_summary,
    age_sex_summary = age_sex_summary,
    coverage = coverage,
    gam_table = gam_table,
    domain_summary = domain_summary,
    domain_age_sex = domain_age_sex,
    oof_performance = oof_perf,
    oof_domain_associations = oof_assoc,
    oof_data = oof_data
  )
}

## =========================================================
## 7. Read input files
## =========================================================

outpatient_files <- file.path(DATA_DIR, c("门诊检验结果1.csv", "门诊检验结果2.csv"))
inpatient_files <- file.path(DATA_DIR, c("住院检验结果1.csv", "住院检验结果2.csv"))

if (any(!file.exists(outpatient_files))) {
  stop("Missing outpatient files: ", paste(outpatient_files[!file.exists(outpatient_files)], collapse = ", "))
}

if (any(!file.exists(inpatient_files))) {
  stop("Missing inpatient files: ", paste(inpatient_files[!file.exists(inpatient_files)], collapse = ", "))
}

## =========================================================
## 8. Outpatient cohorts
## =========================================================

cat("\nReading outpatient data...\n")

out_lab_long <- build_lab_long_outpatient(outpatient_files[1], outpatient_files[2])

write.csv(out_lab_long, file.path(DIR_SOURCE, "SourceData_outpatient_lab_long_clean.csv"), row.names = FALSE)

out_index_all <- build_hb_index_outpatient(out_lab_long)

write.csv(out_index_all, file.path(DIR_SOURCE, "SourceData_outpatient_Hb_index_all_adults.csv"), row.names = FALSE)

out_index_primary <- out_index_all %>%
  filter(anemia_sex_specific_hb_indexed)

out_baseline_primary <- build_baseline_wide(
  lab_long = out_lab_long,
  index_data = out_index_primary,
  setting = "outpatient"
)

write.csv(out_baseline_primary$baseline_long, file.path(DIR_OUTPATIENT, "Outpatient_primary_SourceData_baseline_marker_long.csv"), row.names = FALSE)
write.csv(out_baseline_primary$analysis, file.path(DIR_OUTPATIENT, "Outpatient_primary_SourceData_analysis_baseline_wide.csv"), row.names = FALSE)

outpatient_primary_results <- run_cohort_analysis(
  out_baseline_primary$analysis,
  outdir = DIR_OUTPATIENT,
  prefix = "Outpatient_primary",
  id_cols = c("PID")
)

out_index_hb120 <- out_index_all %>%
  filter(anemia_hb120)

out_baseline_hb120 <- build_baseline_wide(
  lab_long = out_lab_long,
  index_data = out_index_hb120,
  setting = "outpatient"
)

write.csv(out_baseline_hb120$baseline_long, file.path(DIR_OUTPATIENT_HB120, "Outpatient_Hb120_SourceData_baseline_marker_long.csv"), row.names = FALSE)
write.csv(out_baseline_hb120$analysis, file.path(DIR_OUTPATIENT_HB120, "Outpatient_Hb120_SourceData_analysis_baseline_wide.csv"), row.names = FALSE)

outpatient_hb120_results <- run_cohort_analysis(
  out_baseline_hb120$analysis,
  outdir = DIR_OUTPATIENT_HB120,
  prefix = "Outpatient_Hb120",
  id_cols = c("PID")
)

cohort_definition_summary <- bind_rows(
  out_index_all %>%
    filter(anemia_hb120) %>%
    summarise(
      cohort = "Outpatient Hb <120 g/L sensitivity cohort",
      n = n(),
      age_median = median(baseline_age, na.rm = TRUE),
      age_q1 = quantile(baseline_age, 0.25, na.rm = TRUE),
      age_q3 = quantile(baseline_age, 0.75, na.rm = TRUE),
      female_pct = mean(baseline_sex == 0, na.rm = TRUE) * 100,
      hb_median = median(baseline_hb, na.rm = TRUE)
    ),
  out_index_all %>%
    filter(anemia_sex_specific_hb_indexed) %>%
    summarise(
      cohort = "Outpatient primary sex-specific Hb-indexed anemia cohort",
      n = n(),
      age_median = median(baseline_age, na.rm = TRUE),
      age_q1 = quantile(baseline_age, 0.25, na.rm = TRUE),
      age_q3 = quantile(baseline_age, 0.75, na.rm = TRUE),
      female_pct = mean(baseline_sex == 0, na.rm = TRUE) * 100,
      hb_median = median(baseline_hb, na.rm = TRUE)
    ),
  out_index_all %>%
    summarise(
      cohort = "All outpatient Hb-indexed adults",
      n = n(),
      age_median = median(baseline_age, na.rm = TRUE),
      age_q1 = quantile(baseline_age, 0.25, na.rm = TRUE),
      age_q3 = quantile(baseline_age, 0.75, na.rm = TRUE),
      female_pct = mean(baseline_sex == 0, na.rm = TRUE) * 100,
      hb_median = median(baseline_hb, na.rm = TRUE)
    )
)

write.csv(
  cohort_definition_summary,
  file.path(DIR_SOURCE, "Table_outpatient_cohort_definition_comparison.csv"),
  row.names = FALSE
)

## =========================================================
## 9. Inpatient contrast cohort
## =========================================================

cat("\nReading inpatient data...\n")

ip_lab_long <- build_lab_long_inpatient(inpatient_files[1], inpatient_files[2])

write.csv(ip_lab_long, file.path(DIR_SOURCE, "SourceData_inpatient_lab_long_clean.csv"), row.names = FALSE)

ip_index_all <- build_hb_index_inpatient(ip_lab_long)

write.csv(ip_index_all, file.path(DIR_SOURCE, "SourceData_inpatient_Hb_index_all_adults.csv"), row.names = FALSE)

ip_index_contrast <- ip_index_all %>%
  filter(anemia_sex_specific_hb_indexed)

write.csv(ip_index_contrast, file.path(DIR_INPATIENT, "Inpatient_contrast_SourceData_Hb_index.csv"), row.names = FALSE)

ip_baseline_contrast <- build_baseline_wide(
  lab_long = ip_lab_long,
  index_data = ip_index_contrast,
  setting = "inpatient"
)

write.csv(ip_baseline_contrast$baseline_long, file.path(DIR_INPATIENT, "Inpatient_contrast_SourceData_baseline_marker_long.csv"), row.names = FALSE)
write.csv(ip_baseline_contrast$analysis, file.path(DIR_INPATIENT, "Inpatient_contrast_SourceData_analysis_baseline_wide.csv"), row.names = FALSE)

inpatient_contrast_results <- run_cohort_analysis(
  ip_baseline_contrast$analysis,
  outdir = DIR_INPATIENT,
  prefix = "Inpatient_contrast",
  id_cols = c("ADMISSION_NO", "PID")
)

## =========================================================
## 10. Setting contrast summary
## =========================================================

setting_cohort_summary <- bind_rows(
  outpatient_primary_results$cohort_summary %>%
    mutate(setting = "Outpatient primary", .before = 1),
  outpatient_hb120_results$cohort_summary %>%
    mutate(setting = "Outpatient Hb<120 sensitivity", .before = 1),
  inpatient_contrast_results$cohort_summary %>%
    mutate(setting = "Inpatient contrast", .before = 1)
)

write.csv(
  setting_cohort_summary,
  file.path(DIR_COMPARISON, "Table_setting_cohort_summary.csv"),
  row.names = FALSE
)

setting_oof_summary <- bind_rows(
  outpatient_primary_results$oof_performance %>%
    mutate(setting = "Outpatient primary", .before = 1),
  outpatient_hb120_results$oof_performance %>%
    mutate(setting = "Outpatient Hb<120 sensitivity", .before = 1),
  inpatient_contrast_results$oof_performance %>%
    mutate(setting = "Inpatient contrast", .before = 1)
)

write.csv(
  setting_oof_summary,
  file.path(DIR_COMPARISON, "Table_setting_OOF_model_performance_summary.csv"),
  row.names = FALSE
)

## =========================================================
## 11. Main Table 1 and Table 2
## =========================================================

make_table1_one <- function(result_obj, cohort_label) {
  d <- result_obj$analysis_data
  cov <- result_obj$coverage
  n_units <- nrow(d)
  n_patients <- dplyr::n_distinct(d$PID)
  
  age_dist <- d %>%
    filter(!is.na(age_group_4)) %>%
    count(age_group_4, name = "n") %>%
    mutate(
      pct = n / n_units * 100,
      age_group_4 = as.character(age_group_4)
    )
  
  get_age_row <- function(group_label) {
    tmp <- age_dist %>% filter(age_group_4 == group_label)
    if (nrow(tmp) == 0) return(NA_character_)
    fmt_n_pct(tmp$n[1], tmp$pct[1])
  }
  
  get_cov_row <- function(marker_name) {
    tmp <- cov %>% filter(marker == marker_name)
    if (nrow(tmp) == 0) return(NA_character_)
    fmt_n_pct(tmp$n_available[1], tmp$pct_available[1])
  }
  
  female_n <- sum(d$baseline_sex == 0, na.rm = TRUE)
  male_n <- sum(d$baseline_sex == 1, na.rm = TRUE)
  
  tibble(
    Characteristic = c(
      "Analysis units",
      "Unique patients",
      "Age, years, median (IQR)",
      "Age 18–44 years",
      "Age 45–59 years",
      "Age 60–74 years",
      "Age ≥75 years",
      "Female",
      "Male",
      "Hemoglobin, g/L, median (IQR)",
      "Ferritin available",
      "Albumin available",
      "Prealbumin available",
      "Creatinine available",
      "eGFR available",
      "CRP/hsCRP available",
      "SAA available"
    ),
    value = c(
      fmt_int(n_units),
      fmt_int(n_patients),
      fmt_median_iqr(
        median(d$baseline_age, na.rm = TRUE),
        quantile(d$baseline_age, 0.25, na.rm = TRUE),
        quantile(d$baseline_age, 0.75, na.rm = TRUE),
        digits = 1
      ),
      get_age_row("18-44"),
      get_age_row("45-59"),
      get_age_row("60-74"),
      get_age_row(">=75"),
      fmt_n_pct(female_n, female_n / n_units * 100),
      fmt_n_pct(male_n, male_n / n_units * 100),
      fmt_median_iqr(
        median(d$baseline_hb_final, na.rm = TRUE),
        quantile(d$baseline_hb_final, 0.25, na.rm = TRUE),
        quantile(d$baseline_hb_final, 0.75, na.rm = TRUE),
        digits = 1
      ),
      get_cov_row("baseline_ferritin"),
      get_cov_row("baseline_alb"),
      get_cov_row("baseline_prealb"),
      get_cov_row("baseline_creatinine"),
      get_cov_row("baseline_egfr"),
      get_cov_row("baseline_crp_hscrp"),
      get_cov_row("baseline_saa")
    )
  ) %>%
    rename(!!cohort_label := value)
}

table1 <- make_table1_one(
  outpatient_primary_results,
  "Outpatient primary Hb-indexed anemia cohort"
) %>%
  left_join(
    make_table1_one(
      outpatient_hb120_results,
      "Outpatient Hb <120 g/L sensitivity cohort"
    ),
    by = "Characteristic"
  ) %>%
  left_join(
    make_table1_one(
      inpatient_contrast_results,
      "Inpatient Hb-indexed anemia contrast cohort"
    ),
    by = "Characteristic"
  ) %>%
  mutate(
    Section = case_when(
      Characteristic %in% c("Analysis units", "Unique patients") ~ "Cohort size",
      Characteristic %in% c(
        "Age, years, median (IQR)",
        "Age 18–44 years",
        "Age 45–59 years",
        "Age 60–74 years",
        "Age ≥75 years"
      ) ~ "Age",
      Characteristic %in% c("Female", "Male") ~ "Sex",
      Characteristic == "Hemoglobin, g/L, median (IQR)" ~ "Index hemoglobin",
      TRUE ~ "Baseline marker availability"
    ),
    .before = 1
  )

write.csv(
  table1,
  file.path(DIR_TABLES, "Table1_baseline_characteristics_formatted.csv"),
  row.names = FALSE
)

clean_model_label <- function(x) {
  case_when(
    str_detect(x, "full high-coverage") ~ "Full high-coverage feature set",
    str_detect(x, "non-renal") ~ "Non-renal feature set",
    TRUE ~ x
  )
}

table2_raw <- bind_rows(
  outpatient_primary_results$oof_performance %>%
    mutate(Cohort = "Outpatient primary Hb-indexed anemia cohort"),
  outpatient_hb120_results$oof_performance %>%
    mutate(Cohort = "Outpatient Hb <120 g/L sensitivity cohort"),
  inpatient_contrast_results$oof_performance %>%
    mutate(Cohort = "Inpatient Hb-indexed anemia contrast cohort")
) %>%
  mutate(
    Model = clean_model_label(model_name),
    `Renal-metabolic markers included` = if_else(
      Model == "Full high-coverage feature set",
      "Yes",
      "No"
    ),
    `OOF r` = correlation_predicted_vs_chronological,
    `MAE, years` = MAE,
    `RMSE, years` = RMSE,
    `Age-deviation SD, years` = age_gap_sd
  ) %>%
  select(
    Cohort,
    Model,
    `Renal-metabolic markers included`,
    n,
    n_features,
    `OOF r`,
    `MAE, years`,
    `RMSE, years`,
    `Age-deviation SD, years`
  )

write.csv(
  table2_raw,
  file.path(DIR_TABLES, "Table2_OOF_model_performance_numeric.csv"),
  row.names = FALSE
)

table2_formatted <- table2_raw %>%
  mutate(
    `Analysis units` = fmt_int(n),
    `Number of features` = fmt_int(n_features),
    `OOF r` = sprintf("%.3f", `OOF r`),
    `MAE, years` = sprintf("%.1f", `MAE, years`),
    `RMSE, years` = sprintf("%.1f", `RMSE, years`),
    `Age-deviation SD, years` = sprintf("%.1f", `Age-deviation SD, years`)
  ) %>%
  select(
    Cohort,
    Model,
    `Renal-metabolic markers included`,
    `Analysis units`,
    `Number of features`,
    `OOF r`,
    `MAE, years`,
    `RMSE, years`,
    `Age-deviation SD, years`
  )

write.csv(
  table2_formatted,
  file.path(DIR_TABLES, "Table2_OOF_model_performance_formatted.csv"),
  row.names = FALSE
)

## =========================================================
## 12. Save final objects and session info
## =========================================================

final_objects <- list(
  outpatient_index_all = out_index_all,
  inpatient_index_all = ip_index_all,
  outpatient_primary_results = outpatient_primary_results,
  outpatient_hb120_results = outpatient_hb120_results,
  inpatient_contrast_results = inpatient_contrast_results,
  setting_cohort_summary = setting_cohort_summary,
  setting_oof_summary = setting_oof_summary,
  table1 = table1,
  table2_raw = table2_raw,
  table2_formatted = table2_formatted
)

saveRDS(
  final_objects,
  file.path(OUTDIR, "final_anemia_aging_phenomics_analysis_objects.rds")
)

sink(file.path(OUTDIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\n=========================================================\n")
cat("Main analysis pipeline complete.\n")
cat("Outputs saved to:", OUTDIR, "\n")
cat("Main output folders:\n")
cat("1.", DIR_OUTPATIENT, "\n")
cat("2.", DIR_OUTPATIENT_HB120, "\n")
cat("3.", DIR_INPATIENT, "\n")
cat("4.", DIR_COMPARISON, "\n")
cat("5.", DIR_TABLES, "\n")
cat("6.", DIR_SOURCE, "\n")
cat("=========================================================\n")



