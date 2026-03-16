"""Build publication-ready tables/figures from the final FAERS cohort table.

This script is intentionally end-to-end: it loads the final analysis dataset,
runs descriptive and inferential statistics, fits a logistic model, and writes
all outputs into 06_analysis/tables and 06_analysis/figures.
"""

from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import FancyBboxPatch
from matplotlib import patheffects
from scipy.stats import shapiro, mannwhitneyu, kruskal, spearmanr
import scikit_posthocs as sp
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score, roc_curve


ROOT = Path(__file__).resolve().parents[1]
FINAL_CSV = ROOT / "05_final" / "cohort_analysis.csv"
FILTERED = ROOT / "03_filtered"
OUT_DIR = ROOT / "06_analysis"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"

FILE_NAMES = {
    "demo": "teen_demo_records.csv",
    "drug_teens": "teen_drug_records.csv",
    "confirmed_ids": "dph_confirmed_ids.txt",
    "figure_flow": "cohort_flow.png",
    "figure_age_group": "acb_by_age_group.png",
    "figure_forest": "logistic_odds_ratios.png",
    "figure_roc": "cardiac_risk_roc.png",
    "figure_dashboard": "overview_dashboard.png",
    "figure_severity": "acb_vs_severity.png",
    "table_descriptive": "descriptive_summary.csv",
    "table_normality": "normality_checks.csv",
    "table_nonparametric": "nonparametric_results.csv",
    "table_posthoc": "age_group_posthoc.csv",
    "table_logit": "logistic_model_full.csv",
    "table_fit": "logistic_model_fit.csv",
    "table_predictions": "logistic_model_predictions.csv",
    "summary": "results_summary.md",
    "metadata": "results_metadata.json",
}

HISTOGRAM_NAMES = {
    "AGE": "age_histogram.png",
    "total_acb_with_dph": "acb_with_dph_histogram.png",
    "total_acb_codrugs_only": "acb_codrugs_only_histogram.png",
    "n_codrugs": "codrug_count_histogram.png",
    "max_severity": "severity_histogram.png",
}

QQ_NAMES = {
    "AGE": "age_qq_plot.png",
    "total_acb_with_dph": "acb_with_dph_qq_plot.png",
    "total_acb_codrugs_only": "acb_codrugs_only_qq_plot.png",
    "n_codrugs": "codrug_count_qq_plot.png",
    "max_severity": "severity_qq_plot.png",
}


for path in [FIG_DIR, TAB_DIR]:
    path.mkdir(parents=True, exist_ok=True)


PALETTE = {
    "blue": "#1f4e79",
    "teal": "#1f7a8c",
    "orange": "#f59e0b",
    "red": "#b91c1c",
    "slate": "#334155",
    "bg": "#f8fafc",
}


def p_stars(p: float) -> str:
    """Convert a p-value into a compact significance label."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def flow_counts() -> dict[str, int]:
    """Collect row/ID counts used in the cohort flow diagram."""
    counts: dict[str, int] = {}

    demo = pd.read_csv(FILTERED / FILE_NAMES["demo"], dtype=str, usecols=["PRIMARYID"]) 
    counts["DEMO_teens_unique"] = demo["PRIMARYID"].astype(str).str.strip().nunique()

    drug = pd.read_csv(FILTERED / FILE_NAMES["drug_teens"], dtype=str, usecols=["PRIMARYID"], on_bad_lines="skip", engine="python")
    counts["DRUG_teens_unique"] = drug["PRIMARYID"].astype(str).str.strip().nunique()

    ids_txt = FILTERED / FILE_NAMES["confirmed_ids"]
    ids = [x.strip() for x in ids_txt.read_text(encoding="utf-8", errors="ignore").splitlines() if x.strip()]
    counts["DPH_confirmed_unique"] = len(set(ids))

    final = pd.read_csv(FINAL_CSV, dtype=str)
    counts["Final_rows"] = len(final)
    counts["Final_unique_PRIMARYID"] = final["PRIMARYID"].astype(str).str.strip().nunique()
    return counts


def prepare_df() -> pd.DataFrame:
    """Load the final cohort CSV and coerce analysis columns to numeric types."""
    df = pd.read_csv(FINAL_CSV)

    for col in ["AGE", "n_codrugs", "total_acb_with_dph", "total_acb_codrugs_only", "max_severity", "cardiac_any", "cardiac_tier1", "cardiac_tier2", "pre_post_warning"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df["age_group"] = df["age_group"].astype(str)
    df["cardiac_any"] = df["cardiac_any"].fillna(0).astype(int)
    df["pre_post_warning"] = df["pre_post_warning"].fillna(0).astype(int)
    df["YEAR"] = pd.to_numeric(df.get("YEAR", np.nan), errors="coerce")
    return df


def build_descriptive_table(df: pd.DataFrame) -> pd.DataFrame:
    """Create and save descriptive statistics for continuous and categorical variables."""
    continuous = ["AGE", "total_acb_with_dph", "total_acb_codrugs_only", "n_codrugs", "max_severity"]
    categorical = ["age_group", "pre_post_warning", "cardiac_any"]

    rows = []
    for col in continuous:
        s = pd.to_numeric(df[col], errors="coerce").dropna()
        rows.append({
            "variable": col,
            "type": "continuous",
            "n": int(s.shape[0]),
            "mean": float(s.mean()) if not s.empty else np.nan,
            "sd": float(s.std(ddof=1)) if s.shape[0] > 1 else np.nan,
            "median": float(s.median()) if not s.empty else np.nan,
            "iqr": float(s.quantile(0.75) - s.quantile(0.25)) if not s.empty else np.nan,
        })

    for col in categorical:
        vc = df[col].fillna("Missing").astype(str).value_counts(dropna=False)
        n = vc.sum()
        for level, count in vc.items():
            rows.append({
                "variable": col,
                "type": "categorical",
                "level": str(level),
                "count": int(count),
                "percent": float((count / n) * 100.0) if n else np.nan,
            })

    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / FILE_NAMES["table_descriptive"], index=False)
    return out


def normality_checks(df: pd.DataFrame) -> pd.DataFrame:
    """Run Shapiro-Wilk checks and save histogram/Q-Q plots per continuous variable."""
    continuous = ["AGE", "total_acb_with_dph", "total_acb_codrugs_only", "n_codrugs", "max_severity"]
    rows = []
    for col in continuous:
        s = pd.to_numeric(df[col], errors="coerce").dropna()
        if s.shape[0] >= 3:
            stat, p = shapiro(s)
        else:
            stat, p = np.nan, np.nan
        rows.append({"variable": col, "n": int(s.shape[0]), "shapiro_W": stat, "p_value": p, "non_normal_p_lt_0_05": bool(p < 0.05) if pd.notna(p) else None})

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(s, bins=15, color="#4C78A8", edgecolor="white")
        ax.set_title(f"Histogram: {col}")
        ax.set_xlabel(col)
        ax.set_ylabel("Count")
        fig.tight_layout()
        fig.savefig(FIG_DIR / HISTOGRAM_NAMES[col], dpi=200)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(5, 5))
        # Use scipy's probability plot to visualize normality assumptions.
        from scipy import stats
        stats.probplot(s, dist="norm", plot=ax)
        ax.set_title(f"Q-Q: {col}")
        fig.tight_layout()
        fig.savefig(FIG_DIR / QQ_NAMES[col], dpi=200)
        plt.close(fig)

    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / FILE_NAMES["table_normality"], index=False)
    return out


def analysis_a_b_c(df: pd.DataFrame) -> pd.DataFrame:
    """Run core nonparametric analyses and generate key exploratory figures."""
    rows = []

    g0 = df.loc[df["cardiac_any"] == 0, "total_acb_codrugs_only"].dropna()
    g1 = df.loc[df["cardiac_any"] == 1, "total_acb_codrugs_only"].dropna()
    u_stat, p_mw = mannwhitneyu(g0, g1, alternative="two-sided")
    rows.append({"analysis": "A_mann_whitney_acb_vs_cardiac", "statistic": u_stat, "p_value": p_mw})

    groups = [g.dropna().values for _, g in df.groupby("age_group")["total_acb_with_dph"]]
    k_stat, p_kw = kruskal(*groups)
    rows.append({"analysis": "B_kruskal_acb_vs_age_group", "statistic": k_stat, "p_value": p_kw})

    # Post-hoc Dunn test is used after Kruskal-Wallis for pairwise age-group comparisons.
    dunn_df = sp.posthoc_dunn(df, val_col="total_acb_with_dph", group_col="age_group", p_adjust="bonferroni")
    dunn_df.to_csv(TAB_DIR / FILE_NAMES["table_posthoc"])

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    order = ["13-15", "16-17", "18-19"]
    sns.boxplot(
        data=df,
        x="age_group",
        y="total_acb_with_dph",
        order=order,
        hue="age_group",
        legend=False,
        palette=["#dbeafe", "#a7f3d0", "#fde68a"],
        width=0.55,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="age_group",
        y="total_acb_with_dph",
        order=order,
        color=PALETTE["slate"],
        alpha=0.35,
        size=3.5,
        jitter=0.2,
        ax=ax,
    )
    ymax = np.nanmax(df["total_acb_with_dph"]) if df["total_acb_with_dph"].notna().any() else 1
    ax.set_ylim(top=ymax * 1.30 if ymax > 0 else 1.5)
    ax.text(1, ymax * 1.05 if ymax > 0 else 0.5, f"Kruskal p={p_kw:.3g} {p_stars(p_kw)}", ha="center", color=PALETTE["red"], fontweight="bold")

    level_map = {lvl: idx for idx, lvl in enumerate(order)}
    pair_y = ymax * 1.12 if ymax > 0 else 1
    # Draw significance brackets only for statistically significant pairwise contrasts.
    for a in order:
        for b in order:
            if a >= b:
                continue
            pval = dunn_df.loc[a, b]
            if pd.notna(pval) and pval < 0.05:
                xa, xb = level_map[a], level_map[b]
                ax.plot([xa, xa, xb, xb], [pair_y, pair_y * 1.01, pair_y * 1.01, pair_y], color=PALETTE["slate"], lw=1.3)
                ax.text((xa + xb) / 2, pair_y * 1.015, p_stars(float(pval)), ha="center", va="bottom", color=PALETTE["red"], fontsize=11)
                pair_y *= 1.06

    ax.set_title("Total ACB burden across age groups", color=PALETTE["blue"], fontweight="bold")
    ax.set_xlabel("Age group")
    ax.set_ylabel("Total ACB (with DPH)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / FILE_NAMES["figure_age_group"], dpi=220)
    plt.close(fig)

    rho, p_sp = spearmanr(df["total_acb_codrugs_only"], df["max_severity"], nan_policy="omit")
    rows.append({"analysis": "C_spearman_acb_vs_severity", "statistic": rho, "p_value": p_sp})

    fig, ax = plt.subplots(figsize=(7, 5.2))
    sns.regplot(
        data=df,
        x="total_acb_codrugs_only",
        y="max_severity",
        lowess=True,
        scatter_kws={"alpha": 0.75, "s": 40, "color": PALETTE["teal"]},
        line_kws={"color": PALETTE["orange"], "lw": 2.5},
        ax=ax,
    )
    ax.set_title(f"ACB burden versus case severity (rho={rho:.2f}, p={p_sp:.3g})", color=PALETTE["blue"], fontweight="bold")
    ax.set_xlabel("Total ACB (co-drugs only)")
    ax.set_ylabel("Max severity")
    fig.tight_layout()
    fig.savefig(FIG_DIR / FILE_NAMES["figure_severity"], dpi=220)
    plt.close(fig)

    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / FILE_NAMES["table_nonparametric"], index=False)
    return out


def fit_logit(df: pd.DataFrame):
    """Fit a logistic model and gracefully fall back to regularized fitting if needed."""
    formula = "cardiac_any ~ total_acb_codrugs_only + AGE + n_codrugs"
    model_obj = smf.logit(formula=formula, data=df)
    fit_method = "mle"
    try:
        model = model_obj.fit(disp=False)
        coef = model.params
        conf = model.conf_int()
        pvals = model.pvalues
        pred_prob = model.predict(df)
        auc = roc_auc_score(df["cardiac_any"], pred_prob)
        out = pd.DataFrame({
            "term": coef.index,
            "coef": coef.values,
            "OR": np.exp(coef.values),
            "CI_low": np.exp(conf[0].values),
            "CI_high": np.exp(conf[1].values),
            "p_value": pvals.values,
        })
        fit_stats = {
            "n": int(model.nobs),
            "llf": float(model.llf),
            "aic": float(model.aic),
            "bic": float(model.bic),
            "pseudo_r2_mcfadden": float(model.prsquared),
            "auc": float(auc),
            "fit_method": fit_method,
        }
        return model, out, fit_stats, pd.DataFrame({"PRIMARYID": df["PRIMARYID"], "cardiac_any": df["cardiac_any"], "pred_prob": pred_prob})
    except Exception:
        # Separation/convergence issues can occur in sparse safety datasets.
        # Fallback keeps the pipeline running and still yields probabilities/ORs.
        fit_method = "regularized_l1"
        model = model_obj.fit_regularized(disp=False, alpha=0.1, maxiter=1000)
        coef = model.params
        pred_prob = model.predict(df)
        auc = roc_auc_score(df["cardiac_any"], pred_prob)
        out = pd.DataFrame({
            "term": coef.index,
            "coef": coef.values,
            "OR": np.exp(coef.values),
            "CI_low": np.nan,
            "CI_high": np.nan,
            "p_value": np.nan,
        })
        fit_stats = {
            "n": int(df.shape[0]),
            "llf": np.nan,
            "aic": np.nan,
            "bic": np.nan,
            "pseudo_r2_mcfadden": np.nan,
            "auc": float(auc),
            "fit_method": fit_method,
        }
        return model, out, fit_stats, pd.DataFrame({"PRIMARYID": df["PRIMARYID"], "cardiac_any": df["cardiac_any"], "pred_prob": pred_prob})


def logistic_and_forest(df: pd.DataFrame) -> pd.DataFrame:
    """Run logistic analysis, save model outputs, and render forest/ROC plots."""
    fit_rows = []

    model_full, tbl_full, fit_full, pred_df = fit_logit(df)
    tbl_full.to_csv(TAB_DIR / FILE_NAMES["table_logit"], index=False)
    pred_df.to_csv(TAB_DIR / FILE_NAMES["table_predictions"], index=False)
    fit_rows.append({"model": "full", **fit_full})

    forest_df = tbl_full[~tbl_full["term"].eq("Intercept")].copy()
    forest_df = forest_df.sort_values("OR")

    fig, ax = plt.subplots(figsize=(9.2, max(4.2, 0.55 * len(forest_df))))
    y = np.arange(len(forest_df))
    has_ci = forest_df["CI_low"].notna().all() and forest_df["CI_high"].notna().all()
    if has_ci:
        ax.hlines(y, forest_df["CI_low"], forest_df["CI_high"], color=PALETTE["teal"], lw=2.3)
    ax.scatter(forest_df["OR"], y, color=PALETTE["orange"], s=55, zorder=3, edgecolor="white", linewidth=0.8)
    ax.axvline(1.0, color="black", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(forest_df["term"])
    ax.set_xscale("log")
    ax.set_xlabel("Odds Ratio (log scale)")
    ax.set_title("Adjusted odds ratios for cardiac toxicity", color=PALETTE["blue"], fontweight="bold")
    fig.subplots_adjust(left=0.34, right=0.98, top=0.90, bottom=0.12)
    fig.savefig(FIG_DIR / FILE_NAMES["figure_forest"], dpi=220)
    plt.close(fig)

    roc_df = pred_df.dropna(subset=["pred_prob"]).copy()
    fpr, tpr, _ = roc_curve(roc_df["cardiac_any"], roc_df["pred_prob"])
    auc = fit_full["auc"]
    fig, ax = plt.subplots(figsize=(6.8, 6.0))
    ax.plot(fpr, tpr, color=PALETTE["teal"], lw=2.5, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], color=PALETTE["slate"], ls="--", lw=1.2)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC curve for the cardiac risk model", color=PALETTE["blue"], fontweight="bold")
    ax.legend(loc="lower right", frameon=True)
    fig.tight_layout()
    fig.savefig(FIG_DIR / FILE_NAMES["figure_roc"], dpi=220)
    plt.close(fig)

    pd.DataFrame(fit_rows).to_csv(TAB_DIR / FILE_NAMES["table_fit"], index=False)
    return pd.DataFrame(fit_rows)


def plot_flowchart(counts: dict[str, int]) -> None:
    """Render a simple cohort flowchart from precomputed count checkpoints."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")

    boxes = [
        (0.1, 0.80, 0.8, 0.12, f"DEMO teens (unique PRIMARYID): {counts['DEMO_teens_unique']:,}"),
        (0.1, 0.62, 0.8, 0.12, f"DRUG teens (unique PRIMARYID): {counts['DRUG_teens_unique']:,}"),
        (0.1, 0.44, 0.8, 0.12, f"Confirmed DPH cohort: {counts['DPH_confirmed_unique']:,}"),
        (0.1, 0.26, 0.8, 0.12, f"Final analysis rows: {counts['Final_rows']:,} (unique: {counts['Final_unique_PRIMARYID']:,})"),
    ]

    for (x, y, w, h, text) in boxes:
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.02", edgecolor=PALETTE["blue"], facecolor="#e2e8f0", linewidth=1.4)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10, color=PALETTE["slate"], fontweight="bold")

    for y1, y2 in [(0.80, 0.74), (0.62, 0.56), (0.44, 0.38)]:
        ax.annotate("", xy=(0.5, y2), xytext=(0.5, y1), arrowprops=dict(arrowstyle="->", lw=1.4))

    ax.set_title("Cohort flow summary", fontsize=12, color=PALETTE["blue"], fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / FILE_NAMES["figure_flow"], dpi=220)
    plt.close(fig)


def save_summary(counts: dict[str, int], norm: pd.DataFrame, tests: pd.DataFrame, fit: pd.DataFrame) -> None:
    """Write a concise markdown summary of counts, tests, and model fit statistics."""
    lines = [
        "# Analysis Results Summary",
        "",
        "## Cohort",
        f"- DEMO teens unique IDs: {counts['DEMO_teens_unique']}",
        f"- Confirmed DPH cohort: {counts['DPH_confirmed_unique']}",
        f"- Final analysis rows: {counts['Final_rows']}",
        "",
        "## Normality (Shapiro-Wilk)",
    ]
    for row in norm.itertuples(index=False):
        lines.append(f"- {row.variable}: p={row.p_value:.3g} ({'non-normal' if pd.notna(row.p_value) and row.p_value < 0.05 else 'not significant'})")

    lines.append("")
    lines.append("## Core tests")
    for row in tests.itertuples(index=False):
        lines.append(f"- {row.analysis}: statistic={row.statistic:.4g}, p={row.p_value:.3g} {p_stars(row.p_value)}")

    lines.append("")
    lines.append("## Logistic model fit")
    for row in fit.itertuples(index=False):
        lines.append(f"- {row.model}: n={row.n}, pseudo_r2={row.pseudo_r2_mcfadden:.4f}, AIC={row.aic:.2f}, AUC={row.auc:.3f}")

    (OUT_DIR / FILE_NAMES["summary"]).write_text("\n".join(lines), encoding="utf-8")


def make_dashboard() -> None:
    """Assemble a 2x2 dashboard image from the core generated figures."""
    fig = plt.figure(figsize=(14, 10), facecolor=PALETTE["bg"], constrained_layout=True)
    gs = fig.add_gridspec(2, 2, wspace=0.25, hspace=0.3)

    files = [
        FIG_DIR / FILE_NAMES["figure_flow"],
        FIG_DIR / FILE_NAMES["figure_age_group"],
        FIG_DIR / FILE_NAMES["figure_forest"],
        FIG_DIR / FILE_NAMES["figure_roc"],
    ]
    titles = ["Cohort Flow", "ACB by Age Group", "Adjusted Logistic ORs", "Cardiac Risk ROC"]

    for i, (img_path, title) in enumerate(zip(files, titles)):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        img = plt.imread(img_path)
        ax.imshow(img)
        ax.axis("off")
        ax.set_title(title, fontsize=12, fontweight="bold", color=PALETTE["blue"], pad=10)

    fig.suptitle("Adolescent diphenhydramine FAERS results overview", fontsize=15, fontweight="bold", color=PALETTE["blue"])
    fig.savefig(FIG_DIR / FILE_NAMES["figure_dashboard"], dpi=240)
    plt.close(fig)


if __name__ == "__main__":
    # Set a consistent visual style before any plot is created.
    sns.set_theme(style="whitegrid", rc={"axes.facecolor": "#ffffff", "figure.facecolor": PALETTE["bg"]})

    counts = flow_counts()
    plot_flowchart(counts)

    df = prepare_df()
    build_descriptive_table(df)
    norm = normality_checks(df)
    tests = analysis_a_b_c(df)
    fit = logistic_and_forest(df)

    save_summary(counts, norm, tests, fit)
    make_dashboard()

    # Save lightweight metadata so downstream reports can verify input/version context.
    (OUT_DIR / FILE_NAMES["metadata"]).write_text(
        json.dumps({
            "input": str(FINAL_CSV),
            "n_rows": int(df.shape[0]),
            "n_cols": int(df.shape[1]),
            "cohort_confirmed": int(counts["DPH_confirmed_unique"]),
        }, indent=2),
        encoding="utf-8",
    )

    print("Advanced analysis complete.")
    print(f"Outputs: {OUT_DIR}")
