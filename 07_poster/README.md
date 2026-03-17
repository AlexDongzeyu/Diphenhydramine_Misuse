# 07_poster

Poster-facing outputs derived from analysis figures.

## Structure

- `poster_figures/final/`: core figures for the final poster.
- `poster_figures/diagnostics/`: extra diagnostics not intended for the main poster body.
- `poster_figures/docs/`: manifest and layout files.

## How files are created

Run:

```powershell
python 06_analysis/build_poster_figures.py
```

This script copies selected figures from `06_analysis/figures/` and regenerates poster docs.
