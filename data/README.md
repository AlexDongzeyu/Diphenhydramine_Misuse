# data

Source FAERS inputs and lookup tables for cohort construction.

## Subfolders

- `faers_raw/`: downloaded FDA FAERS ZIP archives.
- `faers_extracted/`: extracted quarterly text files (`DEMO`, `DRUG`, `REAC`, `OUTC`).
- `lookups/`: lookup tables such as `acb_lookup.csv` and `drug_name_map.csv`.

## How files are created

- `scripts/download_faers_data.py` populates `faers_raw/` and `faers_extracted/`.
- Lookup files are stored in `lookups/`.

## How it is used next

- `scripts/build_faers_cohort.py` reads `faers_extracted/` + `lookups/` to build filtered and processed cohorts.
