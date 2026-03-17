# RxNorm

## What this folder is

RxNorm source release files used for ingredient-level drug name normalization.

## Required file for current pipeline

- `RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF`

## How it is used

- `scripts/build_faers_cohort.py` reads `RXNCONSO.RRF` to map FAERS drug names to normalized generic/ingredient names in `04_processed/dph_drug_normalized.csv`.
