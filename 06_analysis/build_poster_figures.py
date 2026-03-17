"""Build a poster-ready figure pack.

Outputs:
- Curated high-value analytical figures copied from 06_analysis/figures
- Conceptual diagram 1: pharmacology mechanism schematic
- Conceptual diagram 2: clinical triage workflow
- Conceptual diagram 3: visual abstract of cohort filtering
- Manifest files that separate poster figures from appendix diagnostics
"""

from __future__ import annotations

from pathlib import Path
import csv
import shutil

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle, Polygon
from matplotlib.gridspec import GridSpec
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "06_analysis" / "figures"
POSTER_DIR = ROOT / "07_poster" / "poster_figures"
APPENDIX_DIR = POSTER_DIR / "diagnostics"
FINAL_DIR = POSTER_DIR / "final"
DOCS_DIR = POSTER_DIR / "docs"

COLORS = {
    "blue": "#1f4e79",
    "teal": "#1f7a8c",
    "orange": "#f59e0b",
    "red": "#b91c1c",
    "slate": "#334155",
    "bg": "#f8fafc",
    "mint": "#a8e6c8",
    "sky": "#c7d7ea",
    "sand": "#f3dd86",
    "gray": "#6b7280",
}


KEEP_FIGURES = [
    "drug_cooccurrence_network.png",
    "shap_beeswarm.png",
    "cardiac_risk_roc.png",
    "logistic_odds_ratios.png",
    "acb_by_age_group.png",
]

APPENDIX_FIGURES = [
    "age_histogram.png",
    "age_qq_plot.png",
    "acb_with_dph_histogram.png",
    "acb_with_dph_qq_plot.png",
    "acb_codrugs_only_histogram.png",
    "acb_codrugs_only_qq_plot.png",
    "codrug_count_histogram.png",
    "codrug_count_qq_plot.png",
    "severity_histogram.png",
    "severity_qq_plot.png",
    "overview_dashboard.png",
]


def setup_dirs() -> None:
    POSTER_DIR.mkdir(parents=True, exist_ok=True)
    APPENDIX_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_DIR.mkdir(parents=True, exist_ok=True)


def clean_legacy_root_files() -> None:
    """Remove previously generated root-level files so folder root stays tidy."""
    legacy_files = [
        *KEEP_FIGURES,
        "mechanism_schematic.png",
        "clinical_triage_workflow.png",
        "visual_abstract.png",
        "drug_cooccurrence_network_poster.png",
        "shap_beeswarm_poster.png",
        "cardiac_risk_roc_poster.png",
        "logistic_odds_ratios_poster.png",
        "poster_layout_blueprint.png",
        "poster_figure_manifest.csv",
        "poster_keep_list.txt",
        "layout_blueprint.png",
        "figure_manifest.csv",
        "selected_files.txt",
    ]
    for name in legacy_files:
        p = POSTER_DIR / name
        if p.exists() and p.is_file():
            p.unlink()

    # Remove legacy docs filenames from prior naming conventions.
    legacy_doc_files = [
        "poster_figure_manifest.csv",
        "poster_keep_list.txt",
        "poster_layout_blueprint.png",
    ]
    for name in legacy_doc_files:
        p = DOCS_DIR / name
        if p.exists() and p.is_file():
            p.unlink()


def copy_curated_figures() -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []

    for name in KEEP_FIGURES:
        src = FIG_DIR / name
        dst = FINAL_DIR / name
        if src.exists():
            shutil.copy2(src, dst)
            rows.append((f"final/{name}", "analysis_core", str(src.relative_to(ROOT)), "Poster all-star (unannotated)"))

    for name in APPENDIX_FIGURES:
        src = FIG_DIR / name
        dst = APPENDIX_DIR / name
        if src.exists():
            shutil.copy2(src, dst)
            rows.append((f"diagnostics/{name}", "appendix_only", str(src.relative_to(ROOT)), "Diagnostic/manuscript appendix"))

    return rows


def trim_background(img: np.ndarray, tolerance: float = 0.03, pad: int = 10) -> np.ndarray:
    """Trim large blank margins using the top-left pixel as background reference."""
    if img.ndim < 2:
        return img
    rgb = img[:, :, :3] if img.ndim == 3 else img
    bg = rgb[0, 0]
    diff = np.max(np.abs(rgb - bg), axis=2)
    mask = diff > tolerance
    if not np.any(mask):
        return img

    ys, xs = np.where(mask)
    y0, y1 = max(0, ys.min() - pad), min(img.shape[0], ys.max() + pad + 1)
    x0, x1 = max(0, xs.min() - pad), min(img.shape[1], xs.max() + pad + 1)
    return img[y0:y1, x0:x1]


def draw_mechanism_schematic(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 6), facecolor=COLORS["bg"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.95, "Pharmacological Mechanism: Anticholinergic Blockade to Cardiac Toxicity", ha="center", va="center", fontsize=16, fontweight="bold", color=COLORS["blue"])

    syn_box = FancyBboxPatch((0.05, 0.20), 0.45, 0.62, boxstyle="round,pad=0.02", facecolor="white", edgecolor=COLORS["sky"], linewidth=2)
    ax.add_patch(syn_box)
    ax.text(0.275, 0.79, "Synapse", ha="center", fontsize=12, color=COLORS["slate"], fontweight="bold")

    pre_neuron = Circle((0.15, 0.52), 0.09, facecolor="#e5eef7", edgecolor=COLORS["blue"], linewidth=2)
    post_neuron = Circle((0.39, 0.52), 0.09, facecolor="#edf7ef", edgecolor=COLORS["teal"], linewidth=2)
    ax.add_patch(pre_neuron)
    ax.add_patch(post_neuron)
    ax.text(0.15, 0.52, "ACh\nrelease", ha="center", va="center", fontsize=10, color=COLORS["slate"])
    ax.text(0.39, 0.52, "M2\nreceptor", ha="center", va="center", fontsize=10, color=COLORS["slate"])

    for x in [0.22, 0.25, 0.28, 0.31]:
        ax.add_patch(Circle((x, 0.52), 0.011, facecolor=COLORS["orange"], edgecolor="none"))

    block = Rectangle((0.305, 0.46), 0.02, 0.12, facecolor=COLORS["red"], edgecolor=COLORS["red"])
    ax.add_patch(block)
    ax.text(0.315, 0.60, "Blocked", ha="center", fontsize=9, color=COLORS["red"], fontweight="bold")

    med_box = FancyBboxPatch((0.11, 0.29), 0.31, 0.10, boxstyle="round,pad=0.02", facecolor="#fff7e6", edgecolor=COLORS["orange"], linewidth=1.5)
    ax.add_patch(med_box)
    ax.text(0.265, 0.34, "Diphenhydramine + co-drugs\n(anticholinergic burden)", ha="center", va="center", fontsize=10, color=COLORS["slate"])

    arrow = FancyArrowPatch((0.52, 0.52), (0.69, 0.52), arrowstyle="-|>", mutation_scale=18, linewidth=2.2, color=COLORS["red"])
    ax.add_patch(arrow)
    ax.text(0.605, 0.56, "Autonomic dysregulation", ha="center", fontsize=10, color=COLORS["red"], fontweight="bold")

    heart_box = FancyBboxPatch((0.72, 0.20), 0.23, 0.62, boxstyle="round,pad=0.02", facecolor="white", edgecolor=COLORS["mint"], linewidth=2)
    ax.add_patch(heart_box)
    ax.text(0.835, 0.79, "Cardiac Effect", ha="center", fontsize=12, color=COLORS["slate"], fontweight="bold")

    heart = Polygon([(0.80, 0.58), (0.76, 0.63), (0.74, 0.58), (0.74, 0.53), (0.80, 0.44), (0.86, 0.53), (0.86, 0.58), (0.84, 0.63)], closed=True, facecolor="#fda4af", edgecolor=COLORS["red"], linewidth=1.8)
    ax.add_patch(heart)

    ax.plot([0.75, 0.78, 0.79, 0.80, 0.81, 0.84, 0.86, 0.90], [0.33, 0.33, 0.37, 0.28, 0.40, 0.33, 0.33, 0.33], color=COLORS["red"], linewidth=2)
    ax.text(0.835, 0.26, "QT prolongation / tachyarrhythmia", ha="center", fontsize=10, color=COLORS["slate"])

    # Data bridge: explicitly connect biology schematic to quantified findings.
    data_box = FancyBboxPatch((0.73, 0.08), 0.22, 0.11, boxstyle="round,pad=0.015", facecolor="white", edgecolor=COLORS["red"], linewidth=1.7)
    ax.add_patch(data_box)
    ax.text(0.84, 0.135, "Our data: 40% increased odds\nof severe cardiac events\nper +1 ACB (p=0.062)", ha="center", va="center", fontsize=9.5, color=COLORS["red"], fontweight="bold")

    fig.tight_layout()
    fig.savefig(out_path, dpi=320, bbox_inches="tight")
    plt.close(fig)


def draw_triage_workflow(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 6), facecolor=COLORS["bg"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.94, "Clinical Triage Workflow for Teen Overdose Cases", ha="center", fontsize=16, fontweight="bold", color=COLORS["blue"])

    step1 = FancyBboxPatch((0.05, 0.40), 0.25, 0.24, boxstyle="round,pad=0.02", facecolor="white", edgecolor=COLORS["sky"], linewidth=2)
    step2 = FancyBboxPatch((0.37, 0.40), 0.25, 0.24, boxstyle="round,pad=0.02", facecolor="white", edgecolor=COLORS["sand"], linewidth=2)
    decision = FancyBboxPatch((0.69, 0.38), 0.26, 0.28, boxstyle="round,pad=0.02", facecolor="white", edgecolor=COLORS["mint"], linewidth=2)
    for p in [step1, step2, decision]:
        ax.add_patch(p)

    ax.text(0.175, 0.60, "Step 1", ha="center", fontsize=11, color=COLORS["blue"], fontweight="bold")
    ax.text(0.175, 0.52, "Teen overdose\nadmission", ha="center", fontsize=11, color=COLORS["slate"])

    ax.text(0.495, 0.60, "Step 2", ha="center", fontsize=11, color=COLORS["orange"], fontweight="bold")
    ax.text(0.495, 0.52, "Calculate cumulative\nco-drug ACB", ha="center", fontsize=11, color=COLORS["slate"])

    ax.text(0.82, 0.60, "Step 3", ha="center", fontsize=11, color=COLORS["teal"], fontweight="bold")
    ax.text(0.82, 0.52, "Risk routing", ha="center", fontsize=11, color=COLORS["slate"], fontweight="bold")
    ax.text(0.82, 0.46, "High ACB -> Immediate ECG\nLow ACB -> Standard observation", ha="center", fontsize=10, color=COLORS["slate"])

    ax.add_patch(FancyArrowPatch((0.30, 0.52), (0.37, 0.52), arrowstyle="-|>", mutation_scale=18, linewidth=2, color=COLORS["slate"]))
    ax.add_patch(FancyArrowPatch((0.62, 0.52), (0.69, 0.52), arrowstyle="-|>", mutation_scale=18, linewidth=2, color=COLORS["slate"]))

    fig.tight_layout()
    fig.savefig(out_path, dpi=320, bbox_inches="tight")
    plt.close(fig)


def draw_visual_abstract(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 6), facecolor=COLORS["bg"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.94, "Visual Abstract: FAERS to Confirmed Teen DPH Cohort", ha="center", fontsize=16, fontweight="bold", color=COLORS["blue"])

    box1 = FancyBboxPatch((0.05, 0.34), 0.24, 0.30, boxstyle="round,pad=0.02", facecolor="white", edgecolor=COLORS["sky"], linewidth=2)
    box2 = FancyBboxPatch((0.38, 0.34), 0.24, 0.30, boxstyle="round,pad=0.02", facecolor="white", edgecolor=COLORS["sand"], linewidth=2)
    box3 = FancyBboxPatch((0.71, 0.34), 0.24, 0.30, boxstyle="round,pad=0.02", facecolor="white", edgecolor=COLORS["mint"], linewidth=2)
    for p in [box1, box2, box3]:
        ax.add_patch(p)

    ax.text(0.17, 0.55, "FAERS database", ha="center", fontsize=12, fontweight="bold", color=COLORS["slate"])
    ax.text(0.17, 0.46, "~190k teen\nrecords", ha="center", fontsize=11, color=COLORS["slate"])
    ax.add_patch(Rectangle((0.13, 0.39), 0.08, 0.05, facecolor="#dbeafe", edgecolor=COLORS["blue"], linewidth=1.5))

    ax.text(0.50, 0.55, "Filter pipeline", ha="center", fontsize=12, fontweight="bold", color=COLORS["slate"])
    ax.text(0.50, 0.46, "US/unknown + age 13-19\nDPH-confirmed cases", ha="center", fontsize=11, color=COLORS["slate"])
    funnel = Polygon([(0.46, 0.43), (0.54, 0.43), (0.52, 0.37), (0.48, 0.37)], closed=True, facecolor="#fef3c7", edgecolor=COLORS["orange"], linewidth=1.5)
    ax.add_patch(funnel)

    ax.text(0.83, 0.55, "Final cohort", ha="center", fontsize=12, fontweight="bold", color=COLORS["slate"])
    ax.text(0.83, 0.46, "N = 211 confirmed\nteen DPH overdoses", ha="center", fontsize=11, color=COLORS["slate"])
    ax.add_patch(Circle((0.80, 0.40), 0.025, facecolor="#bfdbfe", edgecolor=COLORS["blue"], linewidth=1.3))
    ax.add_patch(Rectangle((0.83, 0.375), 0.03, 0.05, facecolor="#fde68a", edgecolor=COLORS["orange"], linewidth=1.3))

    ax.add_patch(FancyArrowPatch((0.29, 0.49), (0.38, 0.49), arrowstyle="-|>", mutation_scale=18, linewidth=2, color=COLORS["slate"]))
    ax.add_patch(FancyArrowPatch((0.62, 0.49), (0.71, 0.49), arrowstyle="-|>", mutation_scale=18, linewidth=2, color=COLORS["slate"]))

    fig.tight_layout()
    fig.savefig(out_path, dpi=320, bbox_inches="tight")
    plt.close(fig)


def write_manifest(rows: list[tuple[str, str, str, str]]) -> None:
    manifest = DOCS_DIR / "figure_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "type", "source", "notes"])
        writer.writerows(rows)

    keep_list = DOCS_DIR / "selected_files.txt"
    keep_lines = [
        "Final poster figures (unannotated):",
        "- final/drug_cooccurrence_network.png",
        "- final/shap_beeswarm.png",
        "- final/cardiac_risk_roc.png",
        "- final/logistic_odds_ratios.png",
        "- final/acb_by_age_group.png",
        "",
        "Move these to appendix only:",
        *[f"- {name}" for name in APPENDIX_FIGURES],
    ]
    keep_list.write_text("\n".join(keep_lines), encoding="utf-8")


def build_layout_blueprint(out_path: Path) -> None:
    """Create a one-page visual blueprint showing exact 3-column poster arrangement."""
    fig = plt.figure(figsize=(18, 10), facecolor="white")
    gs = GridSpec(22, 34, figure=fig, wspace=0.12, hspace=0.14)

    def put_image(cell_slice, file_name: str):
        ax = fig.add_subplot(cell_slice)
        img_path = POSTER_DIR / file_name
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.0)
            spine.set_edgecolor("#cbd5e1")
        if img_path.exists():
            ax.imshow(trim_background(plt.imread(img_path), tolerance=0.02, pad=6))
        else:
            ax.set_facecolor("#f1f5f9")

    # Single-copy layout: each core figure appears exactly once.
    put_image(gs[1:10, 0:11], "final/acb_by_age_group.png")
    put_image(gs[1:10, 11:23], "final/shap_beeswarm.png")
    put_image(gs[1:10, 23:34], "final/logistic_odds_ratios.png")
    put_image(gs[10:21, 0:17], "final/drug_cooccurrence_network.png")
    put_image(gs[10:21, 17:34], "final/cardiac_risk_roc.png")

    fig.savefig(out_path, dpi=280, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    # Stage 1: ensure target directories exist and clear legacy clutter.
    setup_dirs()
    clean_legacy_root_files()
    # Remove old folder names from prior conventions.
    for old_dir_name in ["with_title", "without_title", "conceptual", "appendix_diagnostics", "annotated", "base", "concepts"]:
        old_dir = POSTER_DIR / old_dir_name
        if old_dir.exists() and old_dir.is_dir():
            shutil.rmtree(old_dir)

    # Stage 2: copy finalized figure assets from analysis outputs.
    rows = copy_curated_figures()

    rows.extend([
        ("final/drug_cooccurrence_network.png", "analysis_core", "copied from figures", "Unannotated original"),
        ("final/shap_beeswarm.png", "analysis_core", "copied from figures", "Unannotated original"),
        ("final/cardiac_risk_roc.png", "analysis_core", "copied from figures", "Unannotated original"),
        ("final/logistic_odds_ratios.png", "analysis_core", "copied from figures", "Unannotated original"),
        ("final/acb_by_age_group.png", "analysis_core", "copied from figures", "Unannotated original"),
    ])

    blueprint = DOCS_DIR / "layout_blueprint.png"

    # Stage 3: generate supporting documentation artifacts.
    build_layout_blueprint(blueprint)
    rows.append(("docs/layout_blueprint.png", "layout_guide", "generated", "3-column arrangement guide for final poster assembly"))
    write_manifest(rows)

    print(f"Poster figure pack ready: {POSTER_DIR}")


if __name__ == "__main__":
    main()
