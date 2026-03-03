# FAERS Data Placeholders

This project is set up to **reserve space** for FAERS data without downloading it yet.

Folders:
- `data/faers_raw/` — zipped quarterly FAERS downloads
- `data/faers_extracted/` — extracted quarterly files
- `data/faers_input/` — CSVs used by `faers_cardiac_outcomes.py`

When you're ready, place your files in these folders and run:

```bash
python faers_cardiac_outcomes.py --cohort data/faers_input/cohort.csv --reac data/faers_input/reac.csv
```
