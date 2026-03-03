# Pipeline Notes

## Canonical script
- `run_faers_pipeline_compact.py` is the main pipeline script to use.
- It supports recovery starts:
  - `--start-phase 4` full downstream rebuild (from dedupe onward)
  - `--start-phase 6` rebuild teen-linked DRUG/REAC/OUTC from `DEMO_teens` IDs (recommended for quarter-glitch repair)
  - `--start-phase 7` resume from DPH identification onward
  - `--start-phase 8` resume from normalized/aggregated stage when DPH-confirmed files already exist

## Legacy script
- `run_faers_pipeline.py` is kept as a legacy/non-compact variant.

## Recommended command (current workspace layout)
```powershell
& "C:/Users/dongz/OneDrive/Desktop/Project Code/Diphenhydramine_Misuse/.venv/Scripts/python.exe" "run_faers_pipeline_compact.py" --base-dir "." --extracted-root "data/faers_extracted" --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" --acb-csv "data/lookups/acb_lookup.csv" --chunksize 100000 --start-phase 6
```
