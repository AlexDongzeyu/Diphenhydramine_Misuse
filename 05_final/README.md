# 05_final

Purpose: final analysis-ready cohort table.

Files
- cohort_analysis.csv: 211 rows × 14 columns — one row per confirmed adolescent diphenhydramine case.

Columns
- PRIMARYID: unique case identifier
- AGE: patient age in years (13–19)
- SEX: reported sex
- EVENT_DT: adverse event date
- OCCR_COUNTRY: reporting country
- n_codrugs: number of co-administered drugs
- total_acb_with_dph: total anticholinergic burden score including diphenhydramine
- total_acb_codrugs_only: total ACB score for co-drugs only (excluding diphenhydramine)
- cardiac_any: flag — any cardiac outcome present (Tier 1 or Tier 2)
- cardiac_tier1: flag — QT-related signal (QT prolonged, QT interval abnormal, Long QT syndrome)
- cardiac_tier2: flag — serious cardiac event (Torsade de pointes, VT, VFib, Cardiac arrest)
- max_severity: maximum outcome severity code for the case
- YEAR: year of adverse event
- age_group: categorical age group (derived)
- pre_post_warning: indicator relative to FDA label update timing

How it connects
- scripts/build_faers_cohort.py writes this file from 03_filtered and 04_processed inputs.
- 06_analysis/build_analysis_outputs.py reads this file to create all manuscript tables, figures, and summaries.
