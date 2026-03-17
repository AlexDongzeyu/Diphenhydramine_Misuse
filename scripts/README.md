# scripts

Pipeline scripts that build cohort datasets used by downstream analyses.

## Script list

- `download_faers_data.py`: downloads and extracts FAERS quarterly files.
- `build_faers_cohort.py`: builds the primary cohort table.
- `repair_teen_demo.py`: repairs teen DEMO records when needed.
- `add_cardiac_acb_flags.py`: appends cardiac and ACB flags.

## Typical usage order

```powershell
python scripts/download_faers_data.py
python scripts/build_faers_cohort.py --base-dir "." --extracted-root "data/faers_extracted" --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" --acb-csv "data/lookups/acb_lookup.csv" --chunksize 100000 --start-phase 6
```

Then run the analysis scripts in `06_analysis/`.
