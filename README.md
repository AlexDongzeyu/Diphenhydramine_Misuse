# Adolescent Diphenhydramine Misuse and Cardiac Toxicity in FAERS

This repository contains the full data pipeline, cleaned study cohort, and analysis outputs for an adolescent diphenhydramine FAERS project.

## Current submission snapshot
- Active teen DEMO cohort: 217,507 unique PRIMARYIDs in `03_filtered/teen_demo_records.csv`
- Confirmed diphenhydramine cohort: 211 unique PRIMARYIDs in `03_filtered/dph_confirmed_ids.txt`
- Final analysis table: `05_final/cohort_analysis.csv`
- Main pipeline script: `build_faers_cohort.py`
- Main analysis script: `06_analysis/build_analysis_outputs.py`

## Folder map
- `01_raw/`: local copies of raw lookup files used by the pipeline
- `02_combined/`: reserved for optional full-table archival combines
- `03_filtered/`: filtered teen-level and DPH-confirmed FAERS records
- `04_processed/`: normalized drug names and case-level ACB summaries
- `05_final/`: final analysis-ready cohort table
- `06_analysis/`: statistical scripts, tables, figures, and analysis summaries
- `data/`: downloaded FAERS source files and lookup inputs
- `RxNorm/`: RxNorm source release used for generic-name normalization

## Main scripts
- `build_faers_cohort.py`: canonical end-to-end pipeline from raw FAERS tables to the final cohort table
- `build_faers_cohort_legacy.py`: older reference pipeline kept for comparison only
- `download_faers_data.py`: downloads and extracts quarterly FAERS ASCII files
- `repair_teen_demo.py`: targeted DEMO audit and repair utility used to restore the teen DEMO file
- `add_cardiac_acb_flags.py`: standalone helper for joining cardiac outcomes and ACB scores onto a cohort
- `06_analysis/build_analysis_outputs.py`: generates analysis tables, figures, metadata, and the written results summary
- `06_analysis/build_analysis_outputs.R`: R template that mirrors the analysis workflow

## Recommended run order
1. Download or update FAERS source files with `download_faers_data.py` if source data are missing.
2. Build the filtered cohort and final analysis table with `build_faers_cohort.py`.
3. Generate the manuscript tables and figures with `06_analysis/build_analysis_outputs.py`.

## Key outputs
- `03_filtered/teen_demo_records.csv`: canonical adolescent DEMO cohort after dedupe and country/age filtering
- `03_filtered/dph_confirmed_ids.txt`: confirmed diphenhydramine case IDs
- `04_processed/acb_by_case.csv`: case-level ACB totals and co-drug counts
- `05_final/cohort_analysis.csv`: final analysis-ready dataset used by the statistical scripts
- `06_analysis/results_summary.md`: concise analysis summary for the current run
- `06_analysis/results_metadata.json`: run metadata for the current analysis build

Each numbered folder also contains a short `README.md` describing the files in that layer and how they connect to the next step.
