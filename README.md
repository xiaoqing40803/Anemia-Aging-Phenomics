# Anemia-Aging-Phenomics
This repository contains reproducible R scripts for the manuscript on EHR-based laboratory phenomics of Hb-indexed anemia-care populations.

## Repository structure

- `code/00_generate_mock_data.R`: generates fully synthetic mock EHR laboratory data with the same structure expected by the analysis pipeline.
- `code/01_main_analysis_pipeline.R`: performs cohort construction, laboratory marker harmonization, marker-level GAM analyses, hematometabolic domain score derivation, out-of-fold elastic-net modeling, outpatient sensitivity analysis, inpatient care-setting contrast analysis, and generation of source data and main analysis tables.

## Quick start with mock data

```r
source("code/00_generate_mock_data.R")
source("code/01_main_analysis_pipeline.R")
```

The mock data are written to `data_mock/`, and outputs are written to `outputs_anemia_aging_phenomics/`.

## Running with another data directory

```r
Sys.setenv(ANEMIA_PHENOMICS_DATA_DIR = "path/to/csv_folder")
Sys.setenv(ANEMIA_PHENOMICS_OUTDIR = "outputs_real_analysis")
source("code/01_main_analysis_pipeline.R")
```

The expected input files are:

- `门诊检验结果1.csv`
- `门诊检验结果2.csv`
- `住院检验结果1.csv`
- `住院检验结果2.csv`

The mock data are fully synthetic and are provided only to demonstrate the data structure and reproducible workflow.
