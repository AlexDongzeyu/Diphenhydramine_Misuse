# 01_raw

Purpose: local copies of lookup files used directly by the pipeline.

Files
- acb/acb_scores.csv: pipeline-formatted ACB lookup copied from data/lookups/acb_lookup.csv.

How it connects
- scripts/build_faers_cohort.py writes acb/acb_scores.csv here before case-level ACB scoring.
- 04_processed/acb_by_case.csv is built from this lookup plus normalized DPH drug records.
