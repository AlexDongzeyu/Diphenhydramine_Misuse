# 06_analysis

This folder contains analysis scripts and generated outputs.

## Scripts

- `build_analysis_outputs.py`: primary analysis builder.
- `build_analysis_outputs.R`: complementary statistical and network analyses.
- `build_poster_figures.py`: generates poster figures.

## Generated outputs

- `figures/`: generated analysis figures.
- `tables/`: generated analysis tables.
- `results_summary.md`: Python analysis summary.
- `results_summary_r.txt`: R analysis summary.
- `results_metadata.json`: run metadata.

## How to regenerate

```powershell
python 06_analysis/build_analysis_outputs.py
Rscript 06_analysis/build_analysis_outputs.R
python 06_analysis/build_poster_figures.py
```

## Dependencies

- Input: `05_final/cohort_analysis.csv`.
- This folder only writes analysis outputs.
