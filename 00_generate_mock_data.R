

## =========================================================
## 00_generate_mock_data.R
##
## Mock EHR laboratory data generator for:
## Anemia-Aging-Phenomics
##
## This script generates fully synthetic outpatient and inpatient
## laboratory data with the same column structure expected by the
## main analysis pipeline. The mock data are intended only to
## demonstrate and test the reproducible workflow. They do not
## contain real patient information and should not be interpreted
## clinically.
## =========================================================

rm(list = ls())

options(stringsAsFactors = FALSE)
set.seed(20261105)

suppressPackageStartupMessages({
  required_pkgs <- c("dplyr", "tibble", "readr")
  new_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
  if (length(new_pkgs) > 0) {
    install.packages(new_pkgs, type = "binary", dependencies = TRUE)
  }
  invisible(lapply(required_pkgs, library, character.only = TRUE))
})

## -----------------------------
## 1. Output directory
## -----------------------------

OUTDIR <- "data_mock"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 2. Helper functions
## -----------------------------

clip <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

fmt_datetime <- function(x) {
  format(x, "%Y-%m-%d %H:%M:%S")
}

random_datetime <- function(n, start = "2021-01-01", end = "2026-06-30") {
  start <- as.POSIXct(start, tz = "Asia/Shanghai")
  end <- as.POSIXct(end, tz = "Asia/Shanghai")
  start + runif(n, 0, as.numeric(difftime(end, start, units = "secs")))
}

sex_to_cn <- function(sex01) {
  ifelse(sex01 == 1, "男", "女")
}

make_barcode <- function(prefix, n) {
  paste0(prefix, sprintf("%09d", seq_len(n)))
}

lab_ref_range <- function(item) {
  dplyr::case_when(
    item == "血红蛋白" ~ "120-160",
    item == "红细胞压积" ~ "0.35-0.50",
    item == "红细胞计数" ~ "3.5-5.8",
    item == "平均红细胞体积" ~ "80-100",
    item == "平均红细胞血红蛋白量" ~ "27-34",
    item == "平均红细胞血红蛋白浓度" ~ "320-360",
    item == "红细胞分布宽度" ~ "11-15",
    item == "白细胞计数" ~ "3.5-9.5",
    item == "中性粒细胞绝对值" ~ "1.8-6.3",
    item == "淋巴细胞绝对值" ~ "1.1-3.2",
    item == "血小板计数" ~ "125-350",
    item == "铁蛋白" ~ "20-300",
    item == "转铁蛋白饱和度" ~ "20-50",
    item == "铁离子" ~ "7-30",
    item == "总铁结合力" ~ "45-75",
    item == "转铁蛋白" ~ "2.0-3.6",
    item == "C反应蛋白" ~ "0-10",
    item == "高敏C反应蛋白" ~ "0-3",
    item == "淀粉样蛋白A" ~ "0-10",
    item == "红细胞沉降率" ~ "0-20",
    item == "白蛋白" ~ "35-55",
    item == "前白蛋白" ~ "180-390",
    item == "总蛋白" ~ "60-85",
    item == "肌酐" ~ "45-110",
    item == "尿素" ~ "2.5-7.5",
    item == "胱抑素C" ~ "0.55-1.15",
    item == "肾小球滤过率" ~ "90-120",
    item == "尿酸" ~ "150-420",
    TRUE ~ ""
  )
}

lab_unit <- function(item) {
  dplyr::case_when(
    item == "血红蛋白" ~ "g/L",
    item == "红细胞压积" ~ "L/L",
    item == "红细胞计数" ~ "10^12/L",
    item == "平均红细胞体积" ~ "fL",
    item == "平均红细胞血红蛋白量" ~ "pg",
    item == "平均红细胞血红蛋白浓度" ~ "g/L",
    item == "红细胞分布宽度" ~ "%",
    item == "白细胞计数" ~ "10^9/L",
    item == "中性粒细胞绝对值" ~ "10^9/L",
    item == "淋巴细胞绝对值" ~ "10^9/L",
    item == "粒淋比" ~ "",
    item == "血小板计数" ~ "10^9/L",
    item == "铁蛋白" ~ "ng/mL",
    item == "转铁蛋白饱和度" ~ "%",
    item == "铁离子" ~ "umol/L",
    item == "总铁结合力" ~ "umol/L",
    item == "转铁蛋白" ~ "g/L",
    item == "C反应蛋白" ~ "mg/L",
    item == "高敏C反应蛋白" ~ "mg/L",
    item == "淀粉样蛋白A" ~ "mg/L",
    item == "红细胞沉降率" ~ "mm/h",
    item == "白蛋白" ~ "g/L",
    item == "前白蛋白" ~ "mg/L",
    item == "总蛋白" ~ "g/L",
    item == "肌酐" ~ "umol/L",
    item == "尿素" ~ "mmol/L",
    item == "胱抑素C" ~ "mg/L",
    item == "肾小球滤过率" ~ "mL/min/1.73m2",
    item == "尿酸" ~ "umol/L",
    TRUE ~ ""
  )
}

abnormal_flag <- function(value, item) {
  rr <- lab_ref_range(item)
  if (is.na(value) || rr == "" || !grepl("-", rr)) return("")
  parts <- strsplit(rr, "-")[[1]]
  lo <- suppressWarnings(as.numeric(parts[1]))
  hi <- suppressWarnings(as.numeric(parts[2]))
  if (is.na(lo) || is.na(hi)) return("")
  if (value < lo) return("L")
  if (value > hi) return("H")
  ""
}

## -----------------------------
## 3. Synthetic patient generator
## -----------------------------

generate_people <- function(n, setting = c("outpatient", "inpatient")) {
  setting <- match.arg(setting)
  
  if (setting == "outpatient") {
    age <- round(clip(rgamma(n, shape = 4.8, scale = 8.5) + 18, 18, 95), 1)
    sex <- rbinom(n, 1, prob = 0.22)  ## 1 male, 0 female
  } else {
    age <- round(clip(rnorm(n, mean = 60, sd = 16), 18, 95), 1)
    sex <- rbinom(n, 1, prob = 0.46)
  }
  
  age_c <- (age - 45) / 10
  
  renal_burden <- scale(0.08 * age_c + 0.35 * (sex == 1) + rnorm(n, 0, 0.8))[, 1]
  inflamm_burden <- scale(0.03 * age_c + rnorm(n, 0, ifelse(setting == "inpatient", 1.2, 0.9)))[, 1]
  nutrition_depletion <- scale(0.04 * age_c + ifelse(setting == "inpatient", 0.5, 0) + rnorm(n, 0, 0.9))[, 1]
  iron_imbalance <- scale(
    0.8 * (sex == 0) * exp(-((age - 42) / 13)^2) +
      0.15 * age_c +
      rnorm(n, 0, 0.8)
  )[, 1]
  
  hb_mean <- ifelse(
    sex == 1,
    118 - 0.12 * pmax(age - 50, 0),
    106 - 0.08 * pmax(age - 50, 0)
  )
  
  if (setting == "inpatient") {
    hb_mean <- hb_mean - 13 - 3 * inflamm_burden - 2 * nutrition_depletion
  }
  
  hb <- clip(
    hb_mean -
      3.5 * renal_burden -
      2.5 * inflamm_burden -
      2.0 * iron_imbalance +
      rnorm(n, 0, ifelse(setting == "inpatient", 13, 10)),
    45, 155
  )
  
  tibble::tibble(
    pid = paste0(ifelse(setting == "outpatient", "MZ", "ZY"), sprintf("%07d", seq_len(n))),
    sex01 = sex,
    sex = sex_to_cn(sex),
    age = age,
    renal_burden = renal_burden,
    inflamm_burden = inflamm_burden,
    nutrition_depletion = nutrition_depletion,
    iron_imbalance = iron_imbalance,
    hb = hb
  )
}

## -----------------------------
## 4. Synthetic lab value generator
## -----------------------------

generate_values <- function(df, setting = c("outpatient", "inpatient")) {
  setting <- match.arg(setting)
  
  n <- nrow(df)
  age <- df$age
  sex <- df$sex01
  hb <- df$hb
  renal <- df$renal_burden
  infl <- df$inflamm_burden
  nutr <- df$nutrition_depletion
  iron <- df$iron_imbalance
  
  out <- tibble::tibble(
    `血红蛋白` = hb,
    `红细胞压积` = clip(hb / 330 + rnorm(n, 0, 0.025), 0.15, 0.55),
    `红细胞计数` = clip(4.6 + 0.012 * (hb - 110) - 0.008 * (age - 45) + rnorm(n, 0, 0.45), 1.8, 6.5),
    `平均红细胞体积` = clip(88 + 0.05 * (age - 45) - 3.0 * iron + rnorm(n, 0, 7), 55, 115),
    `平均红细胞血红蛋白量` = clip(29 + 0.03 * (age - 45) - 1.1 * iron + rnorm(n, 0, 2.2), 18, 38),
    `平均红细胞血红蛋白浓度` = clip(330 - 0.08 * (age - 45) - 2.5 * iron + rnorm(n, 0, 12), 260, 380),
    `红细胞分布宽度` = clip(13.5 + 0.035 * (age - 45) + 0.55 * iron + 0.25 * infl + rnorm(n, 0, 1.6), 10, 25),
    `白细胞计数` = clip(6.0 - 0.012 * (age - 45) + 0.6 * infl + rnorm(n, 0, 1.5), 1.0, 25),
    `中性粒细胞绝对值` = clip(3.5 - 0.006 * (age - 45) + 0.55 * infl + rnorm(n, 0, 1.2), 0.2, 20),
    `淋巴细胞绝对值` = clip(1.8 - 0.010 * (age - 45) - 0.10 * infl + rnorm(n, 0, 0.45), 0.2, 6),
    `血小板计数` = clip(240 - 0.55 * (age - 45) + 8 * infl + rnorm(n, 0, 55), 35, 750),
    `铁蛋白` = clip(exp(log(45) + 0.018 * (age - 45) + 0.35 * infl + 0.45 * iron + rnorm(n, 0, 0.75)), 3, 5000),
    `转铁蛋白饱和度` = clip(23 - 2.0 * iron - 0.8 * infl + rnorm(n, 0, 7), 3, 80),
    `铁离子` = clip(12 - 1.0 * iron - 0.5 * infl + rnorm(n, 0, 4), 2, 40),
    `总铁结合力` = clip(58 + 2.0 * iron - 1.0 * infl + rnorm(n, 0, 8), 20, 95),
    `转铁蛋白` = clip(2.5 + 0.15 * iron - 0.10 * infl + rnorm(n, 0, 0.35), 0.8, 4.5),
    `C反应蛋白` = clip(exp(log(3) + 0.65 * infl + ifelse(setting == "inpatient", 0.8, 0) + rnorm(n, 0, 0.9)), 0.1, 300),
    `高敏C反应蛋白` = clip(exp(log(2) + 0.55 * infl + rnorm(n, 0, 0.8)), 0.05, 200),
    `淀粉样蛋白A` = clip(exp(log(6) + 0.70 * infl + ifelse(setting == "inpatient", 0.9, 0) + rnorm(n, 0, 1.0)), 0.1, 600),
    `红细胞沉降率` = clip(12 + 0.18 * (age - 45) + 6 * infl + rnorm(n, 0, 10), 1, 150),
    `白蛋白` = clip(43 - 0.035 * (age - 45) - 1.8 * nutr - 0.8 * infl + rnorm(n, 0, 3.2), 20, 55),
    `前白蛋白` = clip(260 - 1.4 * (age - 45) - 22 * nutr - 10 * infl + rnorm(n, 0, 45), 50, 450),
    `总蛋白` = clip(70 - 0.02 * (age - 45) - 1.1 * nutr + rnorm(n, 0, 4), 45, 90),
    `肌酐` = clip(exp(log(65) + 0.018 * (age - 45) + 0.28 * renal + 0.12 * sex + rnorm(n, 0, 0.25)), 25, 550),
    `尿素` = clip(exp(log(5.0) + 0.020 * (age - 45) + 0.22 * renal + rnorm(n, 0, 0.28)), 1.5, 45),
    `胱抑素C` = clip(exp(log(0.85) + 0.018 * (age - 45) + 0.25 * renal + rnorm(n, 0, 0.22)), 0.35, 6),
    `肾小球滤过率` = clip(105 - 0.75 * (age - 45) - 10 * renal + rnorm(n, 0, 12), 10, 160),
    `尿酸` = clip(280 + 1.3 * (age - 45) + 28 * renal + 35 * sex + rnorm(n, 0, 65), 80, 850)
  )
  
  out |>
    dplyr::mutate(
      `粒淋比` = clip(`中性粒细胞绝对值` / pmax(`淋巴细胞绝对值`, 0.1), 0.1, 60)
    )
}

## -----------------------------
## 5. Coverage and long table
## -----------------------------

coverage_prob <- function(item, setting = c("outpatient", "inpatient")) {
  setting <- match.arg(setting)
  
  if (item %in% c(
    "血红蛋白", "红细胞压积", "红细胞计数", "平均红细胞体积",
    "平均红细胞血红蛋白量", "平均红细胞血红蛋白浓度",
    "红细胞分布宽度", "白细胞计数", "中性粒细胞绝对值",
    "淋巴细胞绝对值", "粒淋比", "血小板计数"
  )) return(ifelse(setting == "outpatient", 0.95, 0.98))
  
  if (item %in% c("肌酐", "尿素", "胱抑素C", "肾小球滤过率", "尿酸")) {
    return(ifelse(setting == "outpatient", 0.55, 0.85))
  }
  
  if (item %in% c("白蛋白", "前白蛋白", "总蛋白")) {
    return(ifelse(setting == "outpatient", 0.35, 0.88))
  }
  
  if (item %in% c("C反应蛋白", "高敏C反应蛋白", "淀粉样蛋白A", "红细胞沉降率")) {
    return(ifelse(setting == "outpatient", 0.20, 0.65))
  }
  
  if (item %in% c("铁蛋白", "转铁蛋白饱和度", "铁离子", "总铁结合力", "转铁蛋白")) {
    return(ifelse(setting == "outpatient", 0.18, 0.38))
  }
  
  0.10
}

make_lab_long <- function(people, setting = c("outpatient", "inpatient")) {
  setting <- match.arg(setting)
  
  values <- generate_values(people, setting = setting)
  items <- names(values)
  index_time <- random_datetime(nrow(people))
  rows <- vector("list", length = nrow(people))
  
  for (i in seq_len(nrow(people))) {
    keep <- sapply(items, function(item) runif(1) < coverage_prob(item, setting))
    keep[items == "血红蛋白"] <- TRUE
    
    kept_items <- items[keep]
    n_items <- length(kept_items)
    
    if (setting == "outpatient") {
      offset_days <- sapply(kept_items, function(item) {
        if (item == "血红蛋白") {
          0
        } else if (item %in% c(
          "铁蛋白", "转铁蛋白饱和度", "铁离子", "总铁结合力", "转铁蛋白",
          "C反应蛋白", "高敏C反应蛋白", "淀粉样蛋白A", "红细胞沉降率"
        )) {
          sample(-90:7, 1)
        } else {
          sample(-30:7, 1)
        }
      })
    } else {
      offset_days <- sapply(kept_items, function(item) {
        if (item == "血红蛋白") {
          0
        } else if (item %in% c(
          "铁蛋白", "转铁蛋白饱和度", "铁离子", "总铁结合力", "转铁蛋白",
          "C反应蛋白", "高敏C反应蛋白", "淀粉样蛋白A", "红细胞沉降率"
        )) {
          sample(-7:3, 1)
        } else {
          sample(-3:1, 1)
        }
      })
    }
    
    lab_time <- index_time[i] + offset_days * 86400 + sample(0:(12 * 3600), n_items, replace = TRUE)
    
    out <- tibble::tibble(
      pid = people$pid[i],
      sex = people$sex[i],
      age = people$age[i],
      条码号 = make_barcode(ifelse(setting == "outpatient", "MZB", "ZYB"), n_items),
      检测项目 = kept_items,
      结果 = as.character(round(as.numeric(values[i, kept_items]), 2)),
      单位 = sapply(kept_items, lab_unit),
      检验参考值范围 = sapply(kept_items, lab_ref_range),
      异常提示 = mapply(
        abnormal_flag,
        value = as.numeric(values[i, kept_items]),
        item = kept_items,
        USE.NAMES = FALSE
      ),
      报告时间 = fmt_datetime(lab_time)
    )
    
    if (setting == "inpatient") {
      out <- out |>
        dplyr::mutate(
          ADMINISSION_NO = paste0("ADM", sprintf("%08d", i)),
          PATIENT_NAME = paste0("MockPatient", sprintf("%05d", i)),
          weight = round(clip(rnorm(1, 62, 12), 35, 110), 1),
          hight = round(clip(rnorm(1, 165, 8), 140, 190), 1),
          .before = 1
        )
    }
    
    rows[[i]] <- out
  }
  
  dplyr::bind_rows(rows)
}

## -----------------------------
## 6. Generate and save mock files
## -----------------------------

n_outpatient <- 1800
n_inpatient <- 1200

out_people <- generate_people(n_outpatient, setting = "outpatient")
ip_people <- generate_people(n_inpatient, setting = "inpatient")

out_lab <- make_lab_long(out_people, setting = "outpatient") |>
  dplyr::arrange(pid, 报告时间)

ip_lab <- make_lab_long(ip_people, setting = "inpatient") |>
  dplyr::arrange(ADMINISSION_NO, 报告时间)

out_split <- sample(rep(1:2, length.out = nrow(out_lab)))
ip_split <- sample(rep(1:2, length.out = nrow(ip_lab)))

readr::write_csv(out_lab[out_split == 1, ], file.path(OUTDIR, "门诊检验结果1.csv"))
readr::write_csv(out_lab[out_split == 2, ], file.path(OUTDIR, "门诊检验结果2.csv"))
readr::write_csv(ip_lab[ip_split == 1, ], file.path(OUTDIR, "住院检验结果1.csv"))
readr::write_csv(ip_lab[ip_split == 2, ], file.path(OUTDIR, "住院检验结果2.csv"))

cat("\n=========================================================\n")
cat("Mock data generated successfully.\n")
cat("Output directory:", OUTDIR, "\n")
cat("Outpatient lab rows:", nrow(out_lab), "\n")
cat("Inpatient lab rows:", nrow(ip_lab), "\n")
cat("Outpatient synthetic patients:", dplyr::n_distinct(out_lab$pid), "\n")
cat("Inpatient synthetic admissions:", dplyr::n_distinct(ip_lab$ADMINISSION_NO), "\n")
cat("=========================================================\n")



