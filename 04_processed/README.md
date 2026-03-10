# 04_processed

Purpose: processed case-level features derived from the confirmed diphenhydramine cohort.

Files
- dph_drug_normalized.csv: confirmed DPH DRUG rows with RxNorm-based generic-name normalization.
- acb_by_case.csv: case-level co-drug count and anticholinergic burden summary.

How it connects
- build_faers_cohort.py writes both files after filtering the DPH cohort.
- 05_final/cohort_analysis.csv merges these processed features with DEMO, REAC, and OUTC summaries.
