# 06_analysis/figures

All generated analysis figures.

## Produced by

- `06_analysis/build_analysis_outputs.py`
- `06_analysis/build_analysis_outputs.R`

## Key figure groups

- Core analysis:
	- `cohort_flow.png`
	- `acb_by_age_group.png`
	- `acb_vs_severity.png`
	- `logistic_odds_ratios.png`
	- `cardiac_risk_roc.png`
	- `overview_dashboard.png`

- Model explainability and comparison:
	- `shap_beeswarm.png`
	- `shap_beeswarm_final.png`
	- `calibration_curve_cv.png`
	- `calibration_curve.png`
	- `roc_dual.png`

- Network and subgroup plots:
	- `drug_cooccurrence_network.png`
	- `acb_by_cardiac_outcome.png`
	- `sex_stratified_acb_violin.png`

- Distribution diagnostics:
	- `age_histogram.png`, `age_qq_plot.png`
	- `acb_with_dph_histogram.png`, `acb_with_dph_qq_plot.png`
	- `acb_codrugs_only_histogram.png`, `acb_codrugs_only_qq_plot.png`
	- `codrug_count_histogram.png`, `codrug_count_qq_plot.png`
	- `severity_histogram.png`, `severity_qq_plot.png`

## Regenerate

```powershell
python 06_analysis/build_analysis_outputs.py
Rscript 06_analysis/build_analysis_outputs.R
```
