# Analysis Results Summary

## Cohort
- DEMO teens unique IDs: 217507
- Confirmed DPH cohort: 211
- Final analysis rows: 211

## Normality (Shapiro-Wilk)
- AGE: p=2.92e-09 (non-normal)
- total_acb_with_dph: p=4.85e-17 (non-normal)
- total_acb_codrugs_only: p=1.91e-22 (non-normal)
- n_codrugs: p=6.02e-17 (non-normal)
- max_severity: p=4.07e-13 (non-normal)

## Core tests
- A_mann_whitney_acb_vs_cardiac: statistic=1098, p=0.991 ns
- B_kruskal_acb_vs_age_group: statistic=5.193, p=0.0745 ns
- C_spearman_acb_vs_severity: statistic=0.1373, p=0.0463 *

## Logistic model fit
- full: n=211, pseudo_r2=0.1244, AIC=83.65, AUC=0.753

## CV model performance (leak-free out-of-fold)
- Modeling sample: n_model=211 / n_total=211 (dropped_missing=0)
- logistic_regression_cv: AUC=0.653
- random_forest_cv: AUC=0.820