# 06_analysis

Purpose: final statistical outputs for reporting and submission.

Files
- build_analysis_outputs.py: canonical Python script for tables, figures, summary text, and metadata.
- build_analysis_outputs.R: optional R template mirroring the analysis workflow.
- results_summary.md: concise written summary of the current analysis run.
- results_metadata.json: machine-readable run metadata.
- figures/: manuscript and diagnostic figures.
- tables/: statistical result tables used in the manuscript or poster.

How it connects
- Reads 05_final/cohort_analysis.csv.
- Writes publication-ready outputs without changing upstream cohort files.
