# Adolescent Diphenhydramine Misuse and Cardiac Toxicity in FAERS

Pharmacovigilance study examining anticholinergic burden and cardiac toxicity signals among adolescents (ages 13–19) who received diphenhydramine, using FDA Adverse Event Reporting System (FAERS) data. This repository contains the full data pipeline, cleaned study cohort, and publication-ready analysis outputs.

---

## Cohort snapshot

| Metric | Value |
|---|---|
| Adolescent DEMO records (unique PRIMARYIDs) | **217,507** |
| Confirmed diphenhydramine cases | **211** |
| Final analysis table rows | **211** (one per case, 14 columns) |
| FAERS quarterly coverage | 2013 Q1 – present |

---

## Repository structure

```
Diphenhydramine_Misuse/
├── README.md                       # This file
├── docs/
│   └── pipeline_guide.md           # Quick-start command reference
├── scripts/
│   ├── download_faers_data.py      # Step 1 – download & extract FAERS quarterly ZIPs
│   ├── build_faers_cohort.py       # Step 2 – main end-to-end pipeline (Phases 3–12)
│   ├── build_faers_cohort_legacy.py# Legacy reference pipeline (not for active use)
│   ├── repair_teen_demo.py         # Data audit & repair utility for teen DEMO file
│   └── add_cardiac_acb_flags.py    # Standalone helper – adds cardiac flags & ACB scores
│
├── data/
│   ├── faers_raw/                  # Downloaded FAERS ZIP archives
│   ├── faers_extracted/            # Extracted quarterly DEMO / DRUG / REAC / OUTC files
│   ├── faers_input/                # Placeholder inputs for standalone helpers
│   └── lookups/                    # ACB lookup (acb_lookup.csv) & drug name map
│
├── RxNorm/
│   └── RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF   # Generic-name normalization source
│
├── 01_raw/
│   └── acb/acb_scores.csv          # Pipeline-formatted ACB lookup (copied from data/lookups)
│
├── 02_combined/                    # Optional full-table archival (usually empty)
│
├── 03_filtered/                    # Filtered FAERS records
│   ├── teen_demo_records.csv       # 217,507 adolescent DEMO rows
│   ├── teen_drug_records.csv
│   ├── teen_reaction_records.csv
│   ├── teen_outcome_records.csv
│   ├── dph_confirmed_ids.txt       # 211 confirmed DPH PRIMARYIDs
│   ├── dph_drug_records.csv
│   ├── dph_reaction_records.csv
│   └── dph_outcome_records.csv
│
├── 04_processed/                   # Case-level features
│   ├── dph_drug_normalized.csv     # RxNorm-matched generic drug names
│   └── acb_by_case.csv             # Per-case ACB totals and co-drug counts
│
├── 05_final/
│   └── cohort_analysis.csv         # 211 rows × 14 columns – analysis-ready dataset
│
└── 06_analysis/
    ├── build_analysis_outputs.py   # Step 3 – generates all tables, figures, and summaries
    ├── build_analysis_outputs.R    # Optional R template mirroring the Python workflow
    ├── results_summary.md          # Written summary of current analysis run
    ├── results_metadata.json       # Machine-readable run metadata (n=211, 14 cols)
    ├── figures/                    # 16 PNG figures (flow, diagnostic, manuscript)
    └── tables/                     # 7 CSV result tables
```

---

## Scripts

| Script | Purpose |
|---|---|
| `scripts/download_faers_data.py` | Downloads quarterly FAERS ASCII ZIPs from FDA; extracts DEMO, DRUG, REAC, OUTC files; tracks manifest; supports dry-run and year/quarter range arguments |
| `scripts/build_faers_cohort.py` | Canonical pipeline: combines quarterly tables → filters to adolescents (13–19 yrs, US/unknown) → identifies DPH cases → normalizes drug names via RxNorm → calculates ACB per case → aggregates cardiac outcomes and severity → writes `cohort_analysis.csv`. Supports `--start-phase 4/6/7/8` recovery checkpoints |
| `scripts/build_faers_cohort_legacy.py` | Older reference implementation kept for comparison; not recommended for active use |
| `scripts/repair_teen_demo.py` | Audits and repairs the teen DEMO cohort against raw FAERS source; handles 26→23 column schema remapping and missing-record recovery |
| `scripts/add_cardiac_acb_flags.py` | Standalone helper that joins cardiac outcome flags (Tier 1 QT signals, Tier 2 serious events) and ACB scores onto any cohort CSV by PRIMARYID |
| `06_analysis/build_analysis_outputs.py` | Reads `cohort_analysis.csv`; produces descriptive tables, Shapiro-Wilk normality checks, Mann-Whitney / Kruskal-Wallis / Spearman tests, Dunn post-hoc, logistic regression with forest plot and ROC, and the 16-figure diagnostic panel |
| `06_analysis/build_analysis_outputs.R` | Optional R template mirroring the Python analysis workflow |

---

## Recommended run order

### 1. Download / refresh FAERS source data
```powershell
python scripts/download_faers_data.py
```

### 2. Build the filtered cohort and final analysis table
```powershell
python scripts/build_faers_cohort.py `
  --base-dir "." `
  --extracted-root "data/faers_extracted" `
  --rxnorm-rrf "RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF" `
  --acb-csv "data/lookups/acb_lookup.csv" `
  --chunksize 100000 `
  --start-phase 6
```

Use `--start-phase` to skip already-completed phases (see `docs/pipeline_guide.md` for details).

### 3. Generate manuscript tables and figures
```powershell
python 06_analysis/build_analysis_outputs.py
```

---

## Key outputs

| File | Description |
|---|---|
| `03_filtered/teen_demo_records.csv` | Adolescent DEMO cohort (217,507 unique PRIMARYIDs) after dedupe, age filter (13–19), and US/unknown country filter |
| `03_filtered/dph_confirmed_ids.txt` | 211 confirmed diphenhydramine case IDs |
| `04_processed/acb_by_case.csv` | Per-case ACB totals (`total_acb_with_dph`, `total_acb_codrugs_only`) and co-drug count |
| `05_final/cohort_analysis.csv` | Final 211-row analysis dataset with demographics, ACB burden, cardiac flags, severity, and derived fields |
| `06_analysis/results_summary.md` | Concise written summary of current analysis run |
| `06_analysis/results_metadata.json` | Run metadata: n=211, 14 columns |

---

## Statistical results (current run)

All continuous variables are non-normal by Shapiro-Wilk (p < 0.001); non-parametric tests used throughout.

| Test | Result |
|---|---|
| ACB burden vs cardiac outcome (Mann-Whitney U) | U = 1098, p = 0.991 (ns) |
| ACB burden vs age group (Kruskal-Wallis H) | H = 8.034, p = 0.018 ★ |
| ACB burden vs case severity (Spearman ρ) | ρ = 0.137, p = 0.046 ★ |
| Logistic model (cardiac ~ ACB + age + n_codrugs) | AUC = 0.753, pseudo-R² = 0.124, AIC = 83.65 |

Figures and full result tables are in `06_analysis/figures/` and `06_analysis/tables/`.

---

## Cardiac outcome definitions (Tier classification)

- **Tier 1 – QT signals**: QT prolonged, QT interval abnormal, Long QT syndrome
- **Tier 2 – serious cardiac events**: Torsade de pointes, Ventricular tachycardia, Ventricular fibrillation, Cardiac arrest

---

Each numbered folder contains a short `README.md` describing the files at that pipeline layer and how they connect to the next step.
