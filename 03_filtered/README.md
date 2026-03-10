# 03_filtered

Purpose: filtered FAERS records after adolescent screening and DPH cohort identification.

Files
- teen_demo_records.csv: canonical adolescent DEMO cohort after dedupe and country/age filtering.
- teen_drug_records.csv: DRUG rows restricted to the adolescent DEMO cohort.
- teen_reaction_records.csv: REAC rows restricted to the adolescent DEMO cohort.
- teen_outcome_records.csv: OUTC rows restricted to the adolescent DEMO cohort.
- dph_confirmed_ids.txt: unique PRIMARYIDs for the confirmed diphenhydramine cohort.
- dph_drug_records.csv: DRUG rows for the confirmed DPH cohort.
- dph_reaction_records.csv: REAC rows for the confirmed DPH cohort.
- dph_outcome_records.csv: OUTC rows for the confirmed DPH cohort.

How it connects
- build_faers_cohort.py creates this folder from raw FAERS input.
- 04_processed uses dph_drug_records.csv for drug normalization and ACB scoring.
- 05_final/cohort_analysis.csv combines teen_demo_records.csv with processed ACB and outcome summaries.
