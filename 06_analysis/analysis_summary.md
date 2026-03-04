# Advanced Analysis Summary

## Cohort
- DEMO teens unique IDs: 190295
- Confirmed DPH cohort: 86
- Final analysis rows: 86

## Normality (Shapiro-Wilk)
- AGE: p=9.2e-08 (non-normal)
- total_acb_with_dph: p=4.67e-08 (non-normal)
- total_acb_codrugs_only: p=3.97e-14 (non-normal)
- n_codrugs: p=8.81e-07 (non-normal)
- max_severity: p=6.76e-12 (non-normal)

## Core tests
- A_mann_whitney_acb_vs_cardiac: statistic=436.5, p=0.393 ns
- B_kruskal_acb_vs_age_group: statistic=8.713, p=0.0128 *
- C_spearman_acb_vs_severity: statistic=0.06974, p=0.523 ns

## Logistic model fit
- full: n=86, pseudo_r2=nan, AIC=nan
- pre: n=86, pseudo_r2=0.2049, AIC=59.16