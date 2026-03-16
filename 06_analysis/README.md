# 06_analysis

Purpose: final statistical outputs for reporting and submission.

Files
- build_analysis_outputs.py: canonical Python script that generates tables, figures, summary text, and metadata from 05_final/cohort_analysis.csv.
- build_analysis_outputs.R: optional R template mirroring the Python analysis workflow.
- results_summary.md: concise written summary of the current analysis run.
- results_metadata.json: machine-readable run metadata (n=211, 14 columns).
- figures/: 16 manuscript and diagnostic PNG figures (see figures/README.md).
- tables/: 7 statistical result CSV tables (see tables/README.md).

Key statistical results (current run)
- Cohort: 211 confirmed DPH cases from 217,507 adolescent DEMO records
- All continuous variables non-normal by Shapiro-Wilk (p < 0.001)
- ACB vs cardiac outcome: Mann-Whitney U=1098, p=0.991 (ns)
- ACB vs age group: Kruskal-Wallis H=8.034, p=0.018 ★
- ACB vs severity: Spearman ρ=0.137, p=0.046 ★
- Logistic model (cardiac ~ ACB + age + n_codrugs): AUC=0.753, pseudo-R²=0.124, AIC=83.65

How it connects
- Reads 05_final/cohort_analysis.csv.
- Writes all publication-ready outputs without modifying upstream cohort files.
