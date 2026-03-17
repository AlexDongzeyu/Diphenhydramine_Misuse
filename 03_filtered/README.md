# 03_filtered

Filtered cohort tables used directly by processing and analysis.

## Produced by

- `scripts/build_faers_cohort.py`

## Core outputs

- Teen cohort: `teen_demo_records.csv`, `teen_drug_records.csv`, `teen_reaction_records.csv`, `teen_outcome_records.csv`
- Confirmed diphenhydramine cohort: `dph_confirmed_ids.txt`, `dph_drug_records.csv`, `dph_reaction_records.csv`, `dph_outcome_records.csv`

## Downstream use

- Feeds `04_processed/` feature generation and `05_final/cohort_analysis.csv` assembly.
