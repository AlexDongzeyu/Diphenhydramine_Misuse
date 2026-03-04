from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import FancyBboxPatch
from scipy.stats import shapiro, mannwhitneyu, kruskal, spearmanr
import scikit_posthocs as sp
import statsmodels.formula.api as smf


ROOT = Path(__file__).resolve().parents[1]
FINAL_CSV = ROOT / "05_final" / "analysis_table.csv"
FILTERED = ROOT / "03_filtered"
OUT_DIR = ROOT / "06_analysis"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"


for path in [FIG_DIR, TAB_DIR]:
    path.mkdir(parents=True, exist_ok=True)


def p_stars(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def flow_counts() -> dict[str, int]:
    counts: dict[str, int] = {}

    demo = pd.read_csv(FILTERED / "DEMO_teens.csv", dtype=str, usecols=["PRIMARYID"]) 
    counts["DEMO_teens_unique"] = demo["PRIMARYID"].astype(str).str.strip().nunique()

    drug = pd.read_csv(FILTERED / "DRUG_teens.csv", dtype=str, usecols=["PRIMARYID"], on_bad_lines="skip", engine="python")
    counts["DRUG_teens_unique"] = drug["PRIMARYID"].astype(str).str.strip().nunique()

    ids_txt = FILTERED / "teen_dph_ids.txt"
    ids = [x.strip() for x in ids_txt.read_text(encoding="utf-8", errors="ignore").splitlines() if x.strip()]
    counts["DPH_confirmed_unique"] = len(set(ids))

    final = pd.read_csv(FINAL_CSV, dtype=str)
    counts["Final_rows"] = len(final)
    counts["Final_unique_PRIMARYID"] = final["PRIMARYID"].astype(str).str.strip().nunique()
    return counts


def prepare_df() -> pd.DataFrame:
    df = pd.read_csv(FINAL_CSV)

    for col in ["AGE", "n_codrugs", "total_acb_with_dph", "total_acb_codrugs_only", "max_severity", "cardiac_any", "cardiac_tier1", "cardiac_tier2", "pre_post_warning"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df["age_group"] = df["age_group"].astype(str)
    df["cardiac_any"] = df["cardiac_any"].fillna(0).astype(int)
    df["pre_post_warning"] = df["pre_post_warning"].fillna(0).astype(int)
    return df


def build_descriptive_table(df: pd.DataFrame) -> pd.DataFrame:
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
    out.to_csv(TAB_DIR / "table1_descriptive.csv", index=False)
    return out


def normality_checks(df: pd.DataFrame) -> pd.DataFrame:
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
        fig.savefig(FIG_DIR / f"hist_{col}.png", dpi=200)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(5, 5))
        from scipy import stats
        stats.probplot(s, dist="norm", plot=ax)
        ax.set_title(f"Q-Q: {col}")
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"qq_{col}.png", dpi=200)
        plt.close(fig)

    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / "normality_shapiro.csv", index=False)
    return out


def analysis_a_b_c(df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    g0 = df.loc[df["cardiac_any"] == 0, "total_acb_codrugs_only"].dropna()
    g1 = df.loc[df["cardiac_any"] == 1, "total_acb_codrugs_only"].dropna()
    u_stat, p_mw = mannwhitneyu(g0, g1, alternative="two-sided")
    rows.append({"analysis": "A_mann_whitney_acb_vs_cardiac", "statistic": u_stat, "p_value": p_mw})

    fig, ax = plt.subplots(figsize=(6, 5))
    plot_df = df.copy()
    plot_df["cardiac_any"] = plot_df["cardiac_any"].map({0: "No Cardiac", 1: "Cardiac"})
    sns.boxplot(data=plot_df, x="cardiac_any", y="total_acb_codrugs_only", ax=ax)
    ymax = np.nanmax(plot_df["total_acb_codrugs_only"]) if plot_df["total_acb_codrugs_only"].notna().any() else 1
    ax.text(0.5, ymax * 1.05 if ymax > 0 else 0.5, f"p={p_mw:.3g} {p_stars(p_mw)}", ha="center")
    ax.set_title("Figure 2: ACB (co-drugs only) by cardiac outcome")
    ax.set_xlabel("Cardiac outcome")
    ax.set_ylabel("Total ACB (co-drugs only)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure2_acb_cardiac_boxplot.png", dpi=220)
    plt.close(fig)

    groups = [g.dropna().values for _, g in df.groupby("age_group")["total_acb_with_dph"]]
    k_stat, p_kw = kruskal(*groups)
    rows.append({"analysis": "B_kruskal_acb_vs_age_group", "statistic": k_stat, "p_value": p_kw})

    dunn_df = sp.posthoc_dunn(df, val_col="total_acb_with_dph", group_col="age_group", p_adjust="bonferroni")
    dunn_df.to_csv(TAB_DIR / "posthoc_dunn_age_group.csv")

    fig, ax = plt.subplots(figsize=(7, 5))
    order = ["13-15", "16-17", "18-19"]
    sns.boxplot(data=df, x="age_group", y="total_acb_with_dph", order=order, ax=ax)
    ymax = np.nanmax(df["total_acb_with_dph"]) if df["total_acb_with_dph"].notna().any() else 1
    ax.text(1, ymax * 1.05 if ymax > 0 else 0.5, f"Kruskal p={p_kw:.3g} {p_stars(p_kw)}", ha="center")
    ax.set_title("Figure 3: ACB (with DPH) across age groups")
    ax.set_xlabel("Age group")
    ax.set_ylabel("Total ACB (with DPH)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure3_acb_agegroup_boxplot.png", dpi=220)
    plt.close(fig)

    rho, p_sp = spearmanr(df["total_acb_codrugs_only"], df["max_severity"], nan_policy="omit")
    rows.append({"analysis": "C_spearman_acb_vs_severity", "statistic": rho, "p_value": p_sp})

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.regplot(data=df, x="total_acb_codrugs_only", y="max_severity", lowess=True, scatter_kws={"alpha": 0.8}, ax=ax)
    ax.set_title(f"Figure 4: Spearman ACB vs Severity (rho={rho:.2f}, p={p_sp:.3g})")
    ax.set_xlabel("Total ACB (co-drugs only)")
    ax.set_ylabel("Max severity")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure4_spearman_acb_severity.png", dpi=220)
    plt.close(fig)

    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / "nonparametric_tests.csv", index=False)
    return out


def fit_logit(df: pd.DataFrame, include_pre_post: bool = True):
    preds = ["total_acb_codrugs_only", "C(age_group)", "n_codrugs"]
    if include_pre_post:
        preds.append("pre_post_warning")
    formula = "cardiac_any ~ " + " + ".join(preds)
    model_obj = smf.logit(formula=formula, data=df)
    fit_method = "mle"
    try:
        model = model_obj.fit(disp=False)
        coef = model.params
        conf = model.conf_int()
        pvals = model.pvalues
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
            "fit_method": fit_method,
        }
        return model, out, fit_stats
    except Exception:
        fit_method = "regularized_l1"
        model = model_obj.fit_regularized(disp=False, alpha=0.1, maxiter=1000)
        coef = model.params
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
            "fit_method": fit_method,
        }
        return model, out, fit_stats


def logistic_and_forest(df: pd.DataFrame) -> pd.DataFrame:
    fit_rows = []

    model_full, tbl_full, fit_full = fit_logit(df, include_pre_post=True)
    tbl_full.to_csv(TAB_DIR / "table2_logit_full.csv", index=False)
    fit_rows.append({"model": "full", **fit_full})

    pre = df[df["pre_post_warning"] == 0].copy()
    post = df[df["pre_post_warning"] == 1].copy()

    if pre["cardiac_any"].nunique() > 1:
        _, tbl_pre, fit_pre = fit_logit(pre, include_pre_post=False)
        tbl_pre.to_csv(TAB_DIR / "table2_logit_pre.csv", index=False)
        fit_rows.append({"model": "pre", **fit_pre})
    else:
        pd.DataFrame().to_csv(TAB_DIR / "table2_logit_pre.csv", index=False)

    if post["cardiac_any"].nunique() > 1:
        _, tbl_post, fit_post = fit_logit(post, include_pre_post=False)
        tbl_post.to_csv(TAB_DIR / "table2_logit_post.csv", index=False)
        fit_rows.append({"model": "post", **fit_post})
    else:
        pd.DataFrame().to_csv(TAB_DIR / "table2_logit_post.csv", index=False)

    forest_df = tbl_full[~tbl_full["term"].eq("Intercept")].copy()
    forest_df = forest_df.sort_values("OR")

    fig, ax = plt.subplots(figsize=(8, max(4, 0.45 * len(forest_df))))
    y = np.arange(len(forest_df))
    has_ci = forest_df["CI_low"].notna().all() and forest_df["CI_high"].notna().all()
    if has_ci:
        ax.hlines(y, forest_df["CI_low"], forest_df["CI_high"], color="#4C78A8", lw=2)
    ax.scatter(forest_df["OR"], y, color="#F58518", s=40, zorder=3)
    ax.axvline(1.0, color="black", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels(forest_df["term"])
    ax.set_xscale("log")
    ax.set_xlabel("Odds Ratio (log scale)")
    ax.set_title("Figure 5: Logistic regression ORs (95% CI)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure5_forest_logit_full.png", dpi=220)
    plt.close(fig)

    pd.DataFrame(fit_rows).to_csv(TAB_DIR / "model_fit_stats.csv", index=False)
    return pd.DataFrame(fit_rows)


def plot_flowchart(counts: dict[str, int]) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")

    boxes = [
        (0.1, 0.80, 0.8, 0.12, f"DEMO teens (unique PRIMARYID): {counts['DEMO_teens_unique']:,}"),
        (0.1, 0.62, 0.8, 0.12, f"DRUG teens (unique PRIMARYID): {counts['DRUG_teens_unique']:,}"),
        (0.1, 0.44, 0.8, 0.12, f"Confirmed DPH cohort: {counts['DPH_confirmed_unique']:,}"),
        (0.1, 0.26, 0.8, 0.12, f"Final analysis rows: {counts['Final_rows']:,} (unique: {counts['Final_unique_PRIMARYID']:,})"),
    ]

    for (x, y, w, h, text) in boxes:
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02", edgecolor="black", facecolor="#E8EEF7")
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10)

    for y1, y2 in [(0.80, 0.74), (0.62, 0.56), (0.44, 0.38)]:
        ax.annotate("", xy=(0.5, y2), xytext=(0.5, y1), arrowprops=dict(arrowstyle="->", lw=1.4))

    ax.set_title("Figure 1: Cohort flow summary", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "figure1_flowchart.png", dpi=220)
    plt.close(fig)


def save_summary(counts: dict[str, int], norm: pd.DataFrame, tests: pd.DataFrame, fit: pd.DataFrame) -> None:
    lines = [
        "# Advanced Analysis Summary",
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
        lines.append(f"- {row.model}: n={row.n}, pseudo_r2={row.pseudo_r2_mcfadden:.4f}, AIC={row.aic:.2f}")

    (OUT_DIR / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    sns.set_theme(style="whitegrid")

    counts = flow_counts()
    plot_flowchart(counts)

    df = prepare_df()
    build_descriptive_table(df)
    norm = normality_checks(df)
    tests = analysis_a_b_c(df)
    fit = logistic_and_forest(df)

    save_summary(counts, norm, tests, fit)

    (OUT_DIR / "run_metadata.json").write_text(
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
