# Pipeline Guide

## Canonical pipeline
- `scripts/build_faers_cohort.py` is the main pipeline script.
- It supports recovery starts to avoid re-running completed phases:
  - `--start-phase 4`: full downstream rebuild from DEMO dedupe onward
  - `--start-phase 6`: rebuild teen-linked DRUG, REAC, and OUTC files from `teen_demo_records.csv`
  - `--start-phase 7`: resume from diphenhydramine identification onward
  - `--start-phase 8`: resume from normalized and aggregated files when DPH-confirmed files already exist

## Legacy pipeline
- `scripts/build_faers_cohort_legacy.py` is kept only as an older reference implementation.

## Recommended run order

### Step 1 — Download FAERS source data (if missing or outdated)
```powershell
python scripts/download_faers_data.py
```

### Step 2 — Build the filtered cohort and final analysis table
```powershell
python scripts/build_faers_cohort.py `
  --base-dir "." `
  --extracted-root "data/faers_extracted" `
  --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" `
  --acb-csv "data/lookups/acb_lookup.csv" `
  --chunksize 100000 `
  --start-phase 6
```

### Step 3 — Generate analysis outputs
```powershell
python 06_analysis/build_analysis_outputs.py
```

## Output path summary
- `03_filtered/teen_demo_records.csv`: teen DEMO cohort (217,507 unique PRIMARYIDs)
- `03_filtered/teen_drug_records.csv`: teen DRUG rows
- `03_filtered/teen_reaction_records.csv`: teen REAC rows
- `03_filtered/teen_outcome_records.csv`: teen OUTC rows
- `03_filtered/dph_confirmed_ids.txt`: confirmed diphenhydramine case IDs (211)
- `03_filtered/dph_drug_records.csv`: confirmed DPH DRUG rows
- `03_filtered/dph_reaction_records.csv`: confirmed DPH REAC rows
- `03_filtered/dph_outcome_records.csv`: confirmed DPH OUTC rows
- `04_processed/dph_drug_normalized.csv`: normalized DPH drug names
- `04_processed/acb_by_case.csv`: case-level ACB summary
- `05_final/cohort_analysis.csv`: final analysis-ready cohort table (211 rows × 14 columns)
- `06_analysis/figures/`: 16 PNG figures
- `06_analysis/tables/`: 7 CSV result tables
- `06_analysis/results_summary.md`: written analysis summary
- `06_analysis/results_metadata.json`: machine-readable run metadata
