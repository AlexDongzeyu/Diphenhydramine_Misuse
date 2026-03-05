# Anticholinergic Burden-Based Analysis of Cardiac Toxicity in Adolescent Diphenhydramine Overdose Using FAERS

Youreka Research Project (2026)

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

## What has been done (Alex branch)
- Added an advanced analysis workflow and outputs under `06_analysis/`.
- Implemented reproducible scripts:
   - `06_analysis/run_advanced_analysis.py` (executed in this workspace)
   - `06_analysis/run_advanced_analysis.R` (template for R environment)
- Generated statistical result tables:
   - `06_analysis/tables/table1_descriptive.csv`
   - `06_analysis/tables/normality_shapiro.csv`
   - `06_analysis/tables/nonparametric_tests.csv`
   - `06_analysis/tables/posthoc_dunn_age_group.csv`
   - `06_analysis/tables/table2_logit_full.csv`
   - `06_analysis/tables/table2_logit_pre.csv`
   - `06_analysis/tables/table2_logit_post.csv`
   - `06_analysis/tables/model_fit_stats.csv`
- Generated polished figures for poster/manuscript use:
   - `06_analysis/figures/figure1_flowchart.png`
   - `06_analysis/figures/figure2_acb_cardiac_boxplot.png`
   - `06_analysis/figures/figure3_acb_agegroup_boxplot.png`
   - `06_analysis/figures/figure4_spearman_acb_severity.png`
   - `06_analysis/figures/figure5_forest_logit_full.png`
   - `06_analysis/figures/figure6_dashboard_overview.png`
- Added readability improvements to reduce overlap and clipping (axis headroom, spacing, margins, layout controls).
- Added run metadata and summary:
   - `06_analysis/run_metadata.json`
   - `06_analysis/analysis_summary.md`

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
3. **Finalize inferential reporting for manuscript/presentation**
   - Select primary model specification and finalize interpretation text for ORs/CIs and subgroup findings.
4. **Prepare research deliverables**
   - Methods write-up, variable definitions, cohort flow diagram, and results tables/figures.
