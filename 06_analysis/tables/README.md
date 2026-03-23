# 06_analysis/tables

All generated statistical table outputs.

## Produced by

- Python: `06_analysis/build_analysis_outputs.py`
- R: `06_analysis/build_analysis_outputs.R`

## Key table groups

- Descriptive and diagnostics:
	- `descriptive_summary.csv`
	- `descriptive_continuous.csv`
	- `descriptive_categorical.csv`
	- `normality_checks.csv`

- Hypothesis tests:
	- `nonparametric_results.csv`
	- `age_group_posthoc.csv`
	- `p_values_bh.csv` (R output)
	- `p_values_bh_python_core.csv` (Python output)

- Logistic model outputs:
	- `logistic_model_full.csv`
	- `logistic_model_fit.csv`
	- `logistic_model_predictions.csv`
	- `logistic_model_pre_warning.csv`
	- `logistic_model_post_warning.csv`

- Cross-validation and model comparison:
	- `model_predictions_cv.csv`
	- `model_cv_metrics.csv`
	- `modeling_sample_size.csv`
	- `delong_test.csv`

- Network and subgroup analyses:
	- `network_edge_weight_summary.csv`
	- `network_edge_weight_distribution.csv`
	- `network_threshold_metadata.csv`
	- `sex_stratified_results.csv`

## Regenerate

```powershell
python 06_analysis/build_analysis_outputs.py
Rscript 06_analysis/build_analysis_outputs.R
```
