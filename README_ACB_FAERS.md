# FAERS + ACB

This project now supports Anticholinergic Cognitive Burden (ACB) scoring per FAERS case.

## Required input files
Place files here:
- `data/faers_input/cohort.csv` (must include `PRIMARYID`)
- `data/faers_input/reac.csv` (must include `PRIMARYID`, `PT`)
- `data/faers_input/drug.csv` (must include `PRIMARYID` and `DRUGNAME` or `PROD_AI`)
- `data/lookups/acb_lookup.csv` (provided; columns: `generic_name`, `acb_score`)

## Run
```bash
python faers_cardiac_outcomes.py
```

## Output files
- `cohort_with_cardiac_outcomes.csv`:
  - Cardiac flags and tier per case
- `reac_cardiac_cases.csv`:
  - REAC rows matching the hardcoded MedDRA cardiac PTs
- `cohort_with_cardiac_and_acb.csv`:
  - Cardiac + ACB case-level results
  - Includes:
    - `acb_total_score`
    - `acb_max_score`
    - `acb_matched_drug_count`
    - `acb_clinically_relevant` (1 if total score >= 3)
- `drug_acb_matches.csv`:
  - DRUG rows successfully matched to ACB lookup terms

## Notes
- The ACB lookup is in `data/lookups/acb_lookup.csv` and can be edited manually.
- Matching is case-insensitive and normalized for punctuation/spaces.
- No FAERS download is required to keep this project structure ready.
