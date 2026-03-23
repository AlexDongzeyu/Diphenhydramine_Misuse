# Diphenhydramine Misuse Project (FAERS)

This repository contains a full, reproducible FAERS pipeline built from raw quarterly files to final analysis figures.

The central analysis question is:

- How co-medication anticholinergic burden is associated with cardiac toxicity and overall severity in confirmed adolescent diphenhydramine overdose reports.

### Step 0 - Raw inputs prepared

We started with quarterly FAERS files in `data/faers_extracted/` (DEMO, DRUG, REAC, OUTC by quarter), lookup tables in `data/lookups/`, and RxNorm files in `RxNorm/RxNorm_full_prescribe_03022026/rrf/`.

These three sources are the required inputs for cohort construction:

- FAERS quarterly case-level records,
- anticholinergic burden lookup values,
- and drug normalization vocabulary (RxNorm).

### Step 1 - Build the cohort foundation (`01_raw`, `02_combined`, `03_filtered`)

The cohort build process is run by `scripts/build_faers_cohort.py`.

What was done:

- FAERS raw quarterly tables were read and harmonized across years.
- The age filter was applied to isolate adolescent/teen records.
- Diphenhydramine exposure was confirmed from drug tables.
- Matched case IDs were used to keep only records tied to the confirmed cohort.

Where outputs were written:

- `01_raw/` stores stable staging lookup content used by the pipeline.
- `02_combined/` is reserved for large intermediate combined tables when needed.
- `03_filtered/` stores the first clean cohort-level outputs.

Key outputs created in `03_filtered/`:

- `teen_demo_records.csv`
- `teen_drug_records.csv`
- `teen_reaction_records.csv`
- `teen_outcome_records.csv`
- `dph_confirmed_ids.txt`
- `dph_drug_records.csv`
- `dph_reaction_records.csv`
- `dph_outcome_records.csv`

### Step 2 - Feature engineering (`04_processed`)

After the confirmed teen diphenhydramine cohort was fixed, drug-level data were converted into case-level analysis features.

What was done:

- Drug names were normalized to reduce spelling/format variation.
- Co-medication counts were computed per case.
- Anticholinergic burden was aggregated per case.

Outputs created in `04_processed/`:

- `dph_drug_normalized.csv` (normalized drug naming table)
- `acb_by_case.csv` (case-level burden and co-medication features)

### Step 3 - Final modeling dataset (`05_final`)

Processed features were merged into one final analysis table.

What was done:

- Demographic, exposure, and outcome fields were aligned by case ID.
- Modeling variables were standardized into one schema.
- Pre/post warning-period indicator and outcome targets were finalized.

Output created in `05_final/`:

- `cohort_analysis.csv`

This file is the single source used by both Python and R analyses.

### Step 4 - Statistical analysis and diagnostics (`06_analysis`)

Two analysis scripts were used to generate the full quantitative evidence package.

Scripts run:

- `06_analysis/build_analysis_outputs.py`
- `06_analysis/build_analysis_outputs.R`

What was done:

- Descriptive summaries were generated for cohort characterization.
- Distribution checks and nonparametric tests were produced.
- Logistic models were fit for cardiac and severity outcomes.
- Cross-validation, ROC/AUC, calibration, and DeLong comparisons were computed.
- Multiple-testing adjustment (BH) outputs were generated.
- Network and subgroup outputs were generated for interpretation support.

Main result folders:

- `06_analysis/tables/`
- `06_analysis/figures/`

### Step 5 - Poster figure build (`07_poster`)

Final presentation figures were built from analysis outputs with `06_analysis/build_poster_figures.py`.

What was done:

- Analysis figures were curated into a publication/presentation set.
- Visual formatting and labeling were standardized for poster use.
- Diagnostic/supporting visuals were kept separate from final poster assets.

Outputs:

- Final poster figures: `07_poster/poster_figures/final/`
- Diagnostics: `07_poster/poster_figures/diagnostics/`
- Figure documentation/support: `07_poster/poster_figures/docs/`

## Run Order

Run these commands in sequence:

```powershell
# Step 1: Build cohort (fresh run — processes all FAERS quarters from scratch)
python scripts/build_faers_cohort.py --base-dir "." --extracted-root "data/faers_extracted" --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" --acb-csv "data/lookups/acb_lookup.csv" --chunksize 100000

# Step 2-4: Generate analysis outputs and poster figures
python 06_analysis/build_analysis_outputs.py
Rscript 06_analysis/build_analysis_outputs.R
python 06_analysis/build_poster_figures.py
```

**Re-run shortcuts** (when cohort files in `03_filtered/` already exist):

```powershell
# Re-run from phase 6 — rebuilds DRUG/REAC/OUTC filters using existing teen_demo_records.csv
python scripts/build_faers_cohort.py ... --start-phase 6

# Re-run from phase 8 — rebuilds feature joins using existing confirmed DPH cohort files
python scripts/build_faers_cohort.py ... --start-phase 8
```

The commands above execute the same process summarized in Steps 1-5.

## Final Deliverables

- `05_final/cohort_analysis.csv`: single analysis-ready cohort dataset.
- `06_analysis/tables/`: inferential and model result tables.
- `06_analysis/figures/`: analysis and diagnostic figure outputs.
- `07_poster/poster_figures/final/`: final curated figures for presentation.

