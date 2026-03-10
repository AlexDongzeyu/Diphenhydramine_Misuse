# 05_final

Purpose: final analysis-ready cohort table.

Files
- cohort_analysis.csv: one row per confirmed adolescent diphenhydramine case with demographics, ACB burden, cardiac flags, severity, and derived analysis fields.

How it connects
- build_faers_cohort.py writes this file from 03_filtered and 04_processed inputs.
- 06_analysis/build_analysis_outputs.py reads this file to create the manuscript tables, figures, and summaries.
