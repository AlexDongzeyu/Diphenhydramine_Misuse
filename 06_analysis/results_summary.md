# Analysis Results Summary

## Cohort
- DEMO teens unique IDs: 217507
- Confirmed DPH cohort: 3168
- Final analysis rows: 3168

## Normality (Shapiro-Wilk)
- AGE: p=9.23e-37 (non-normal)
- total_acb_with_dph: p=4.85e-63 (non-normal)
- total_acb_codrugs_only: p=1.92e-71 (non-normal)
- n_codrugs: p=2.73e-60 (non-normal)
- max_severity: p=3.59e-47 (non-normal)

## Core tests
- A_mann_whitney_acb_vs_cardiac: statistic=3.758e+05, p=1.3e-06 ***
- B_kruskal_acb_vs_age_group: statistic=19.52, p=5.76e-05 ***
- C_spearman_acb_vs_severity: statistic=0.1424, p=8.06e-16 ***

## Logistic model fit
- full: n=3168, pseudo_r2=0.0294, AIC=1553.90, AUC=0.652

## CV model performance (leak-free out-of-fold)
- Modeling sample: n_model=3168 / n_total=3168 (dropped_missing=0)
- logistic_regression_cv: AUC=0.640
- random_forest_cv: AUC=0.769