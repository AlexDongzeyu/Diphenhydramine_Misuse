# Diphenhydramine Misuse Project (FAERS)

This repository contains an end-to-end, reproducible workflow for studying adolescent diphenhydramine misuse in FAERS.

The central analysis question is:

- How anticholinergic burden from co-medications relates to cardiac toxicity risk and case severity in confirmed teen diphenhydramine overdose reports.

The workflow transforms raw/extracted FAERS records into:

- a clean final cohort table,
- publication-style statistical outputs,
- and a poster-ready figure package.

## End-to-End Process (01_raw to Final Figures)

### 01_raw — Pipeline lookup staging

Purpose:

- Store stable lookup copies used during feature engineering.

Key content:

- `acb/acb_scores.csv` (copied from lookup source data).

Why this step matters:

- It fixes lookup inputs for reproducible anticholinergic burden computation.

### 02_combined — Large intermediate combine area

Purpose:

- Hold large intermediate combined FAERS tables when needed.

How it is used now:

- Standard downstream analysis can run without requiring this folder.

### 03_filtered — Cohort filtering and diphenhydramine confirmation

Purpose:

- Build the analysis cohort from extracted FAERS files by applying demographic and case filters.

Main outputs:

- Teen-level filtered tables:
	- `teen_demo_records.csv`
	- `teen_drug_records.csv`
	- `teen_reaction_records.csv`
	- `teen_outcome_records.csv`
- Confirmed diphenhydramine cohort files:
	- `dph_confirmed_ids.txt`
	- `dph_drug_records.csv`
	- `dph_reaction_records.csv`
	- `dph_outcome_records.csv`

Why this step matters:

- It defines the exact population and case set used by all downstream analyses.

### 04_processed — Drug normalization and case-level feature derivation

Purpose:

- Convert filtered drug data into analysis features.

Main outputs:

- `dph_drug_normalized.csv`: normalized/standardized drug names.
- `acb_by_case.csv`: case-level anticholinergic burden and co-medication counts.

Why this step matters:

- It creates interpretable and model-ready predictors rather than raw text drug records.

### 05_final — Final analysis dataset

Purpose:

- Assemble one analysis-ready table for all statistical modeling and figure generation.

Main output:

- `cohort_analysis.csv`

Typical columns include:

- patient/report context (`PRIMARYID`, `AGE`, `SEX`, `YEAR`, `age_group`),
- burden/exposure features (`n_codrugs`, `total_acb_with_dph`, `total_acb_codrugs_only`),
- outcomes (`cardiac_any`, `cardiac_tier1`, `cardiac_tier2`, `max_severity`),
- period indicator (`pre_post_warning`).

Why this step matters:

- It is the single source of truth for both Python and R analysis scripts.

### 06_analysis — Statistical analysis and figure generation

Purpose:

- Produce all formal tables, metrics, and figures from `05_final/cohort_analysis.csv`.

Scripts:

- `build_analysis_outputs.py`: primary analysis engine.
- `build_analysis_outputs.R`: complementary analyses and additional outputs.

Representative table outputs (`06_analysis/tables/`):

- descriptive summaries,
- normality and nonparametric test results,
- logistic model outputs,
- cross-validated model performance,
- DeLong ROC comparison,
- BH-adjusted p-value table,
- network summaries,
- sex-stratified effect-size results.

Representative figure outputs (`06_analysis/figures/`):

- cohort flow,
- age-group ACB comparison,
- severity relationship,
- logistic odds-ratio forest plot,
- ROC and dual ROC,
- calibration curves,
- SHAP explainability beeswarm,
- drug co-occurrence network,
- subgroup/diagnostic plots.

Why this step matters:

- It generates the complete evidence package used for interpretation, reporting, and poster construction.

### 07_poster — Final figure curation for presentation

Purpose:

- Build the exact poster figure set from analysis outputs.

Main output path:

- `07_poster/poster_figures/final/`

Supporting paths:

- `07_poster/poster_figures/diagnostics/`
- `07_poster/poster_figures/docs/`

Why this step matters:

- It separates final presentation figures from diagnostics and keeps poster assets reproducible.

## Command Order

Run these commands in sequence:

```powershell
python scripts/build_faers_cohort.py --base-dir "." --extracted-root "data/faers_extracted" --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" --acb-csv "data/lookups/acb_lookup.csv" --chunksize 100000 --start-phase 6
python 06_analysis/build_analysis_outputs.py
Rscript 06_analysis/build_analysis_outputs.R
python 06_analysis/build_poster_figures.py
```

## Final Deliverables

- `05_final/cohort_analysis.csv`: analysis-ready cohort dataset.
- `06_analysis/tables/`: statistical tables for results reporting.
- `06_analysis/figures/`: full analysis and diagnostic figure set.
- `07_poster/poster_figures/final/`: final curated poster figures.

