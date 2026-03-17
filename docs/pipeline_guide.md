# Pipeline Guide

## Purpose

Use this guide to reproduce the analysis outputs from the final cohort.

## Run order

### 1) Build or refresh final cohort data
```powershell
python scripts/build_faers_cohort.py --base-dir "." --extracted-root "data/faers_extracted" --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" --acb-csv "data/lookups/acb_lookup.csv" --chunksize 100000 --start-phase 6
```

### 2) Generate Python analysis tables and figures
```powershell
python 06_analysis/build_analysis_outputs.py
```

### 3) Generate complementary R outputs
```powershell
Rscript 06_analysis/build_analysis_outputs.R
```

### 4) Generate poster figures
```powershell
python 06_analysis/build_poster_figures.py
```

## Primary outputs

- `05_final/cohort_analysis.csv`
- `06_analysis/tables/`
- `06_analysis/figures/`
- `07_poster/poster_figures/final/`
