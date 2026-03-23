# scripts

Pipeline scripts that build cohort datasets used by downstream analyses.

## Script list

### Main pipeline

- `download_faers_data.py`: downloads and extracts FAERS quarterly ZIP files into `data/faers_extracted/`.
- `build_faers_cohort.py`: full cohort build pipeline — filters FAERS to teen diphenhydramine cases, normalizes drug names via RxNorm, computes ACB and cardiac flags, and writes outputs through `05_final/`.
- `build_faers_cohort_legacy.py`: earlier version of the cohort builder, retained for reference.

### Utility / recovery scripts

- `repair_teen_demo.py`: repairs incomplete teen DEMO records when a partial pipeline run left gaps.
- `add_cardiac_acb_flags.py`: standalone script to re-append cardiac outcome and ACB features to an existing cohort without re-running the full pipeline.

## Typical usage order

```powershell
# 1. Download raw FAERS data (only needed if faers_extracted/ is empty)
python scripts/download_faers_data.py

# 2. Build the full cohort (fresh run)
python scripts/build_faers_cohort.py --base-dir "." --extracted-root "data/faers_extracted" --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" --acb-csv "data/lookups/acb_lookup.csv" --chunksize 100000
```

Then run the analysis scripts in `06_analysis/`. See `docs/pipeline_guide.md` for re-run shortcuts using `--start-phase`.
