# Anticholinergic Burden-Based Analysis of Cardiac Toxicity in Adolescent Diphenhydramine Overdose Using FAERS

Eureka Research Project (2026)

## Project status
- Pipeline is operational and currently completes through final output using `run_faers_pipeline_compact.py`.
- Core outputs are present:
  - `03_filtered/DEMO_teens.csv`
  - `03_filtered/DRUG_teens.csv`
  - `03_filtered/REAC_teens.csv`
  - `03_filtered/OUTC_teens.csv`
  - `03_filtered/DRUG_dph_confirmed.csv`
  - `04_processed/DRUG_normalized.csv`
  - `04_processed/ACB_per_case.csv`
  - `05_final/analysis_table.csv`
- Current confirmed adolescent DPH cohort size: **86 unique PRIMARYID**.
- Quarter-mismatch/data-loss issue in teen-linked tables was repaired by rebuilding from **Phase 6** and validating ID-subset consistency.

## Current workspace notes
- This workspace uses:
  - FAERS extracted source at `data/faers_extracted/`
  - RxNorm source at `RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF`
  - ACB lookup at `01_raw/acb/acb_scores.csv`
- `run_faers_pipeline_compact.py` is the **main script** for ongoing runs and recovery.
- `run_faers_pipeline.py` is retained as a legacy/non-compact variant.

## What still needs to be done
1. **Finalize structure harmonization** (optional but recommended)
   - Align raw folder layout to one standard (`01_raw/...` vs current mixed `data/...` + `RxNorm/...`) for long-term reproducibility.
2. **Regenerate `02_combined/*_all.csv` if strict archival completeness is required**
   - These files are very large and currently not the active dependency for compact runs.
3. **Run and document analysis/statistics** from `05_final/analysis_table.csv`
   - Define final statistical plan (modeling, subgroup analysis, sensitivity checks).
4. **Prepare research deliverables**
   - Methods write-up, variable definitions, cohort flow diagram, and results tables/figures.

## Minimal run commands
Full downstream rebuild (dedupe onward):
```powershell
& "C:/Users/dongz/OneDrive/Desktop/Project Code/Diphenhydramine_Misuse/.venv/Scripts/python.exe" "run_faers_pipeline_compact.py" --base-dir "." --extracted-root "data/faers_extracted" --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" --acb-csv "data/lookups/acb_lookup.csv" --chunksize 100000 --start-phase 4
```

Targeted teen-table repair (recommended if quarter mismatch reappears):
```powershell
& "C:/Users/dongz/OneDrive/Desktop/Project Code/Diphenhydramine_Misuse/.venv/Scripts/python.exe" "run_faers_pipeline_compact.py" --base-dir "." --extracted-root "data/faers_extracted" --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" --acb-csv "data/lookups/acb_lookup.csv" --chunksize 100000 --start-phase 6
```
