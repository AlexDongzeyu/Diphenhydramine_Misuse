# data

Purpose: FAERS source downloads and lookup inputs.

Folders
- faers_raw/: downloaded FAERS ZIP files.
- faers_extracted/: extracted quarterly DEMO, DRUG, REAC, and OUTC text files.
- faers_input/: placeholder inputs for standalone helper scripts.
- lookups/: drug-name and ACB lookup tables used by the pipeline.

How it connects
- scripts/download_faers_data.py fills faers_raw and faers_extracted.
- scripts/build_faers_cohort.py reads faers_extracted and lookups during cohort construction.
