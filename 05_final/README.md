# 05_final

Final analysis-ready cohort data.

## Produced by

- `scripts/build_faers_cohort.py`

## Main file

- `cohort_analysis.csv`

## Analysis-critical columns

- `PRIMARYID`: unique FAERS case identifier.
- `AGE`, `SEX`, `EVENT_DT`, `OCCR_COUNTRY`: core demographics/report metadata.
- `n_codrugs`: number of concomitant co-medications.
- `total_acb_with_dph`: total anticholinergic burden including diphenhydramine.
- `total_acb_codrugs_only`: anticholinergic burden from co-medications only.
- `cardiac_any`, `cardiac_tier1`, `cardiac_tier2`: cardiac toxicity indicators.
- `max_severity`: maximum severity indicator for the case.
- `YEAR`, `age_group`, `pre_post_warning`: derived analysis fields.

## Downstream use

- Input to `06_analysis/build_analysis_outputs.py` and `06_analysis/build_analysis_outputs.R`.
