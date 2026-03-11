from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
from matplotlib.patches import FancyBboxPatch
from matplotlib import patheffects
from scipy import stats
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
    "navy": "#12344d",
    "blue": "#2a6f97",
    "teal": "#2a9d8f",
    "mint": "#8ecae6",
    "gold": "#e9c46a",
    "orange": "#f4a261",
    "rose": "#c1666b",
    "red": "#b23a48",
    "slate": "#52616b",
    "ink": "#1f2933",
    "grid": "#d9e2ec",
    "panel": "#ffffff",
    "bg": "#f4f7fb",
    "shadow": "#cbd5e1",
}

SERIES_COLORS = [PALETTE["blue"], PALETTE["teal"], PALETTE["gold"], PALETTE["orange"], PALETTE["rose"]]

VARIABLE_LABELS = {
    "AGE": "Age (years)",
    "total_acb_with_dph": "Total ACB burden including diphenhydramine",
    "total_acb_codrugs_only": "Total ACB burden from co-drugs only",
    "n_codrugs": "Number of co-drugs",
    "max_severity": "Maximum case severity",
}

TERM_LABELS = {
    "total_acb_codrugs_only": "ACB burden from co-drugs",
    "AGE": "Age (years)",
    "n_codrugs": "Number of co-drugs",
}

CARDIAC_LABELS = {0: "No cardiac toxicity", 1: "Cardiac toxicity"}
CARDIAC_COLORS = {0: PALETTE["blue"], 1: PALETTE["rose"]}
AGE_ORDER = ["13-15", "16-17", "18-19"]
EXPORT_DPI = 320
RNG = np.random.default_rng(42)


def p_stars(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def pretty_label(name: str) -> str:
    return VARIABLE_LABELS.get(name, name.replace("_", " ").title())


def format_p_value(p: float) -> str:
    if pd.isna(p):
        return "p = NA"
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"


def add_shadow(artist) -> None:
    artist.set_path_effects([
        patheffects.withSimplePatchShadow(offset=(2, -2), shadow_rgbFace=PALETTE["shadow"], alpha=0.28),
        patheffects.Normal(),
    ])


def style_axes(ax, grid_axis: str = "y") -> None:
    ax.set_facecolor(PALETTE["panel"])
    ax.tick_params(colors=PALETTE["ink"], labelsize=11)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    for side in ["left", "bottom"]:
        ax.spines[side].set_color(PALETTE["grid"])
        ax.spines[side].set_linewidth(1.2)

    ax.grid(False)
    if grid_axis in {"x", "both"}:
        ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.85, alpha=0.85)
    if grid_axis in {"y", "both"}:
        ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.85, alpha=0.85)
    ax.set_axisbelow(True)


def add_title_block(ax, title: str, subtitle: str = "") -> None:
    ax.set_title(title, loc="left", fontsize=16, fontweight="bold", color=PALETTE["navy"], pad=18)
    if subtitle:
        ax.text(
            0,
            1.02,
            subtitle,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=10.8,
            color=PALETTE["slate"],
        )


def add_badge(
    ax,
    text: str,
    x: float = 0.98,
    y: float = 0.98,
    ha: str = "right",
    va: str = "top",
    facecolor: str = "#ffffff",
) -> None:
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=9.8,
        color=PALETTE["ink"],
        linespacing=1.4,
        bbox={
            "boxstyle": "round,pad=0.45,rounding_size=0.12",
            "facecolor": facecolor,
            "edgecolor": PALETTE["grid"],
            "linewidth": 1.0,
        },
    )


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=EXPORT_DPI, bbox_inches="tight", facecolor=fig.get_facecolor(), edgecolor="none")
    plt.close(fig)


def configure_theme() -> None:
    sns.set_theme(style="whitegrid", font="DejaVu Serif")
    plt.rcParams.update(
        {
            "figure.facecolor": PALETTE["bg"],
            "axes.facecolor": PALETTE["panel"],
            "axes.edgecolor": PALETTE["grid"],
            "axes.labelcolor": PALETTE["ink"],
            "axes.titlecolor": PALETTE["navy"],
            "text.color": PALETTE["ink"],
            "xtick.color": PALETTE["ink"],
            "ytick.color": PALETTE["ink"],
            "grid.color": PALETTE["grid"],
            "font.size": 11,
            "axes.titlesize": 16,
            "axes.labelsize": 12,
            "legend.frameon": False,
            "legend.fontsize": 10,
        }
    )


def flow_counts() -> dict[str, int]:
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
    continuous = ["AGE", "total_acb_with_dph", "total_acb_codrugs_only", "n_codrugs", "max_severity"]
    rows = []
    for idx, col in enumerate(continuous):
        s = pd.to_numeric(df[col], errors="coerce").dropna()
        if s.shape[0] >= 3:
            stat, p = shapiro(s)
        else:
            stat, p = np.nan, np.nan
        rows.append({"variable": col, "n": int(s.shape[0]), "shapiro_W": stat, "p_value": p, "non_normal_p_lt_0_05": bool(p < 0.05) if pd.notna(p) else None})

        color = SERIES_COLORS[idx % len(SERIES_COLORS)]

        fig, ax = plt.subplots(figsize=(9, 5.4), facecolor=PALETTE["bg"])
        sns.histplot(
            s,
            bins="fd" if s.nunique() > 1 else 1,
            kde=True,
            stat="count",
            color=color,
            alpha=0.9,
            edgecolor="white",
            linewidth=1.1,
            line_kws={"linewidth": 2.3, "color": PALETTE["navy"]},
            ax=ax,
        )
        ax.axvline(s.mean(), color=PALETTE["navy"], linewidth=2.2, linestyle="-", label="Mean")
        ax.axvline(s.median(), color=PALETTE["orange"], linewidth=2.2, linestyle="--", label="Median")
        add_title_block(ax, f"Distribution of {pretty_label(col)}", "Histogram with kernel density estimate and central tendency markers")
        style_axes(ax, grid_axis="y")
        ax.set_xlabel(pretty_label(col))
        ax.set_ylabel("Case count")
        ax.legend(loc="upper left")
        add_badge(
            ax,
            "\n".join(
                [
                    f"n = {len(s):,}",
                    f"mean = {s.mean():.2f}",
                    f"median = {s.median():.2f}",
                    f"Shapiro W = {stat:.3f}" if pd.notna(stat) else "Shapiro W = NA",
                    format_p_value(p),
                ]
            ),
            facecolor="#fbfdff",
        )
        save_figure(fig, FIG_DIR / HISTOGRAM_NAMES[col])

        fig, ax = plt.subplots(figsize=(6.6, 6.2), facecolor=PALETTE["bg"])
        (theoretical, ordered), (slope, intercept, r_value) = stats.probplot(s, dist="norm")
        ax.scatter(theoretical, ordered, s=38, color=color, alpha=0.85, edgecolors="white", linewidths=0.55)
        ax.plot(theoretical, slope * np.asarray(theoretical) + intercept, color=PALETTE["navy"], linewidth=2.2)
        add_title_block(ax, f"Q-Q Plot for {pretty_label(col)}", "Departure from the reference line indicates non-normality")
        style_axes(ax, grid_axis="both")
        ax.set_xlabel("Theoretical quantiles")
        ax.set_ylabel("Observed quantiles")
        add_badge(ax, f"Correlation = {r_value:.3f}\n{format_p_value(p)}", facecolor="#fbfdff")
        save_figure(fig, FIG_DIR / QQ_NAMES[col])

    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / FILE_NAMES["table_normality"], index=False)
    return out


def analysis_a_b_c(df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    g0 = df.loc[df["cardiac_any"] == 0, "total_acb_codrugs_only"].dropna()
    g1 = df.loc[df["cardiac_any"] == 1, "total_acb_codrugs_only"].dropna()
    u_stat, p_mw = mannwhitneyu(g0, g1, alternative="two-sided")
    rows.append({"analysis": "A_mann_whitney_acb_vs_cardiac", "statistic": u_stat, "p_value": p_mw})

    groups = [g.dropna().values for _, g in df.groupby("age_group")["total_acb_with_dph"]]
    k_stat, p_kw = kruskal(*groups)
    rows.append({"analysis": "B_kruskal_acb_vs_age_group", "statistic": k_stat, "p_value": p_kw})

    dunn_df = sp.posthoc_dunn(df, val_col="total_acb_with_dph", group_col="age_group", p_adjust="bonferroni")
    dunn_df.to_csv(TAB_DIR / FILE_NAMES["table_posthoc"])

    fig, ax = plt.subplots(figsize=(8.8, 6.4), facecolor=PALETTE["bg"])
    order = AGE_ORDER
    plot_df = df[df["age_group"].isin(order)].copy()
    violin = sns.violinplot(
        data=plot_df,
        x="age_group",
        y="total_acb_with_dph",
        order=order,
        hue="age_group",
        palette=[PALETTE["mint"], "#b7e4c7", "#f8d49d"],
        legend=False,
        inner=None,
        cut=0,
        linewidth=0,
        saturation=1,
        ax=ax,
    )
    for collection in violin.collections:
        collection.set_alpha(0.7)
    sns.boxplot(
        data=plot_df,
        x="age_group",
        y="total_acb_with_dph",
        order=order,
        width=0.22,
        showfliers=False,
        boxprops={"facecolor": "white", "edgecolor": PALETTE["navy"], "linewidth": 1.6, "alpha": 0.96},
        whiskerprops={"color": PALETTE["navy"], "linewidth": 1.4},
        capprops={"color": PALETTE["navy"], "linewidth": 1.4},
        medianprops={"color": PALETTE["red"], "linewidth": 2.0},
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="age_group",
        y="total_acb_with_dph",
        order=order,
        color=PALETTE["ink"],
        alpha=0.26,
        size=3.4,
        jitter=0.2,
        ax=ax,
    )
    ymax = np.nanmax(plot_df["total_acb_with_dph"]) if plot_df["total_acb_with_dph"].notna().any() else 1
    ymin = np.nanmin(plot_df["total_acb_with_dph"]) if plot_df["total_acb_with_dph"].notna().any() else 0
    lower = min(0, ymin - max(0.25, 0.05 * max(ymax, 1)))
    ax.set_ylim(lower, ymax * 1.34 if ymax > 0 else 1.5)

    summary = plot_df.groupby("age_group")["total_acb_with_dph"].agg(["median", "count"]).reindex(order)
    for idx, row in enumerate(summary.itertuples()):
        ax.scatter(idx, row.median, s=65, color=PALETTE["red"], edgecolors="white", linewidths=0.8, zorder=4)
        ax.text(idx, lower + 0.02 * (ax.get_ylim()[1] - lower), f"n = {int(row.count)}", ha="center", va="bottom", fontsize=10, color=PALETTE["slate"])

    level_map = {lvl: idx for idx, lvl in enumerate(order)}
    pair_y = ymax * 1.10 if ymax > 0 else 1
    for i, a in enumerate(order):
        for b in order[i + 1 :]:
            pval = dunn_df.loc[a, b] if a in dunn_df.index and b in dunn_df.columns else np.nan
            if pd.notna(pval) and pval < 0.05:
                xa, xb = level_map[a], level_map[b]
                ax.plot([xa, xa, xb, xb], [pair_y, pair_y * 1.015, pair_y * 1.015, pair_y], color=PALETTE["slate"], lw=1.3)
                ax.text((xa + xb) / 2, pair_y * 1.02, f"{p_stars(float(pval))} ({format_p_value(float(pval))})", ha="center", va="bottom", color=PALETTE["red"], fontsize=10.2)
                pair_y *= 1.06

    add_title_block(ax, "Total ACB Burden Across Age Groups", f"Kruskal-Wallis: H = {k_stat:.2f}, {format_p_value(p_kw)}")
    style_axes(ax, grid_axis="y")
    ax.set_xlabel("Age group")
    ax.set_ylabel("Total ACB (with DPH)")
    add_badge(ax, "Violin = density\nBox = median and IQR\nRed point = group median", x=0.02, ha="left", facecolor="#fbfdff")
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_age_group"])

    rho, p_sp = spearmanr(df["total_acb_codrugs_only"], df["max_severity"], nan_policy="omit")
    rows.append({"analysis": "C_spearman_acb_vs_severity", "statistic": rho, "p_value": p_sp})

    fig, ax = plt.subplots(figsize=(8.4, 6.1), facecolor=PALETTE["bg"])
    severity_df = df[["total_acb_codrugs_only", "max_severity", "cardiac_any"]].dropna().copy()
    severity_df["severity_jittered"] = severity_df["max_severity"] + RNG.normal(0, 0.05, size=len(severity_df))
    for cardiac_value, label in CARDIAC_LABELS.items():
        subset = severity_df[severity_df["cardiac_any"] == cardiac_value]
        ax.scatter(
            subset["total_acb_codrugs_only"],
            subset["severity_jittered"],
            s=34,
            alpha=0.45,
            color=CARDIAC_COLORS[cardiac_value],
            edgecolors="white",
            linewidths=0.45,
            label=label,
        )
    sns.regplot(
        data=severity_df,
        x="total_acb_codrugs_only",
        y="max_severity",
        lowess=True,
        scatter=False,
        line_kws={"color": PALETTE["navy"], "lw": 2.8},
        ax=ax,
    )
    add_title_block(ax, "ACB Burden Versus Maximum Case Severity", f"Spearman correlation: rho = {rho:.2f}, {format_p_value(p_sp)}")
    style_axes(ax, grid_axis="both")
    ax.set_xlabel("Total ACB (co-drugs only)")
    ax.set_ylabel("Max severity")
    ax.set_yticks(sorted(severity_df["max_severity"].dropna().unique()))
    ax.legend(loc="upper left", title="Outcome subgroup")
    add_badge(
        ax,
        "\n".join(
            [
                f"n = {len(severity_df):,}",
                f"Median ACB = {severity_df['total_acb_codrugs_only'].median():.2f}",
                f"Median severity = {severity_df['max_severity'].median():.2f}",
            ]
        ),
        facecolor="#fbfdff",
    )
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_severity"])

    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / FILE_NAMES["table_nonparametric"], index=False)
    return out


def fit_logit(df: pd.DataFrame):
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
    fit_rows = []

    model_full, tbl_full, fit_full, pred_df = fit_logit(df)
    tbl_full.to_csv(TAB_DIR / FILE_NAMES["table_logit"], index=False)
    pred_df.to_csv(TAB_DIR / FILE_NAMES["table_predictions"], index=False)
    fit_rows.append({"model": "full", **fit_full})

    forest_df = tbl_full[~tbl_full["term"].eq("Intercept")].copy()
    forest_df["label"] = forest_df["term"].map(TERM_LABELS).fillna(forest_df["term"])
    forest_df["significant"] = forest_df["p_value"].fillna(1) < 0.05
    forest_df = forest_df.sort_values("OR", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10.2, max(4.4, 0.72 * len(forest_df))), facecolor=PALETTE["bg"])
    y = np.arange(len(forest_df))
    has_ci = forest_df["CI_low"].notna().all() and forest_df["CI_high"].notna().all()
    value_text = []
    for yi, row in forest_df.iterrows():
        point_color = PALETTE["teal"] if row["significant"] else PALETTE["slate"]
        if has_ci:
            ax.plot([row["CI_low"], row["CI_high"]], [yi, yi], color=point_color, lw=3, solid_capstyle="round", alpha=0.95)
            value_text.append(f"{row['OR']:.2f} ({row['CI_low']:.2f}-{row['CI_high']:.2f})")
        else:
            value_text.append(f"{row['OR']:.2f}")
        ax.scatter(row["OR"], yi, color=point_color, s=90, zorder=3, edgecolor="white", linewidth=0.9)

    ax.axvline(1.0, color=PALETTE["red"], ls="--", lw=1.4)
    ax.set_yticks(y)
    ax.set_yticklabels(forest_df["label"])
    ax.set_xscale("log")
    x_low = float(np.nanmin(forest_df["CI_low"] if has_ci else forest_df["OR"]))
    x_high = float(np.nanmax(forest_df["CI_high"] if has_ci else forest_df["OR"]))
    ax.set_xlim(max(0.1, x_low * 0.75), x_high * 2.15)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda value, _: f"{value:g}"))
    style_axes(ax, grid_axis="x")
    ax.set_xlabel("Odds Ratio (log scale)")
    add_title_block(
        ax,
        "Adjusted Odds Ratios for Cardiac Toxicity",
        f"Logistic model adjusted for age and co-drug count; n = {fit_full['n']}, AUC = {fit_full['auc']:.3f}",
    )
    text_x = ax.get_xlim()[1] / 1.02
    for yi, value in enumerate(value_text):
        ax.text(text_x, yi, value, ha="right", va="center", fontsize=10, color=PALETTE["ink"])
    add_badge(ax, "Teal = p < 0.05\nDashed line = null effect", x=0.02, ha="left", facecolor="#fbfdff")
    fig.subplots_adjust(left=0.29, right=0.97, top=0.88, bottom=0.13)
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_forest"])

    roc_df = pred_df.dropna(subset=["pred_prob"]).copy()
    fpr, tpr, _ = roc_curve(roc_df["cardiac_any"], roc_df["pred_prob"])
    auc = fit_full["auc"]
    youden_idx = int(np.argmax(tpr - fpr))
    fig, ax = plt.subplots(figsize=(7.1, 6.4), facecolor=PALETTE["bg"])
    ax.fill_between(fpr, tpr, color=PALETTE["mint"], alpha=0.35)
    ax.plot(fpr, tpr, color=PALETTE["teal"], lw=3.0, label=f"Model ROC (AUC = {auc:.3f})")
    ax.scatter(fpr[youden_idx], tpr[youden_idx], color=PALETTE["orange"], s=85, edgecolors="white", linewidths=0.8, zorder=3)
    ax.plot([0, 1], [0, 1], color=PALETTE["slate"], ls="--", lw=1.2, label="Chance")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    style_axes(ax, grid_axis="both")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    add_title_block(ax, "ROC Curve for the Cardiac Risk Model", "Receiver operating characteristic with highlighted best tradeoff point")
    ax.legend(loc="lower right")
    add_badge(
        ax,
        "\n".join(
            [
                f"Best sensitivity = {tpr[youden_idx]:.2f}",
                f"Best specificity = {1 - fpr[youden_idx]:.2f}",
                f"AUC = {auc:.3f}",
            ]
        ),
        facecolor="#fbfdff",
    )
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_roc"])

    pd.DataFrame(fit_rows).to_csv(TAB_DIR / FILE_NAMES["table_fit"], index=False)
    return pd.DataFrame(fit_rows)


def plot_flowchart(counts: dict[str, int]) -> None:
    fig, ax = plt.subplots(figsize=(9, 7), facecolor=PALETTE["bg"])
    ax.axis("off")

    stages = [
        ("Teen DEMO records", counts["DEMO_teens_unique"], PALETTE["mint"]),
        ("Teen DRUG records", counts["DRUG_teens_unique"], "#d8f3dc"),
        ("Confirmed diphenhydramine cohort", counts["DPH_confirmed_unique"], "#fdecc8"),
        ("Final analytic rows", counts["Final_rows"], "#ffe0db"),
    ]

    y_positions = [0.80, 0.60, 0.40, 0.20]
    previous = None
    for stage_index, ((label, value, fill), y) in enumerate(zip(stages, y_positions), start=1):
        retention = "Baseline" if previous is None else f"{(value / previous) * 100:.1f}% retained"
        rect = FancyBboxPatch(
            (0.10, y),
            0.80,
            0.13,
            boxstyle="round,pad=0.025,rounding_size=0.03",
            edgecolor=PALETTE["navy"],
            facecolor=fill,
            linewidth=1.6,
        )
        add_shadow(rect)
        ax.add_patch(rect)
        ax.text(0.14, y + 0.088, f"Step {stage_index}", fontsize=9.5, color=PALETTE["slate"], fontweight="bold")
        ax.text(0.14, y + 0.054, label, fontsize=13, color=PALETTE["navy"], fontweight="bold")
        ax.text(0.14, y + 0.024, f"n = {value:,}", fontsize=11, color=PALETTE["ink"])
        ax.text(0.86, y + 0.054, retention, ha="right", fontsize=10.2, color=PALETTE["slate"])
        previous = value

    for upper, lower in [(0.80, 0.73), (0.60, 0.53), (0.40, 0.33)]:
        ax.annotate(
            "",
            xy=(0.5, lower),
            xytext=(0.5, upper),
            arrowprops={"arrowstyle": "-|>", "lw": 1.8, "color": PALETTE["navy"]},
        )

    ax.text(0.10, 0.96, "Cohort Flow Summary", fontsize=17, color=PALETTE["navy"], fontweight="bold")
    ax.text(0.10, 0.925, "Progression from all adolescent FAERS records to the analytic diphenhydramine cohort", fontsize=10.8, color=PALETTE["slate"])
    ax.text(0.10, 0.08, f"Final unique PRIMARYID count: {counts['Final_unique_PRIMARYID']:,}", fontsize=10.5, color=PALETTE["ink"])
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_flow"])


def save_summary(counts: dict[str, int], norm: pd.DataFrame, tests: pd.DataFrame, fit: pd.DataFrame) -> None:
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
        panel = FancyBboxPatch((0.01, 0.01), 0.98, 0.98, transform=ax.transAxes, boxstyle="round,pad=0.012,rounding_size=0.03", facecolor="none", edgecolor=PALETTE["grid"], linewidth=1.2)
        ax.add_patch(panel)
        ax.set_title(title, fontsize=12.5, fontweight="bold", color=PALETTE["navy"], pad=10)

    fig.suptitle("Adolescent Diphenhydramine FAERS Results Overview", fontsize=17, fontweight="bold", color=PALETTE["navy"])
    fig.text(0.5, 0.965, "Publication-style summary panel for the main cohort, association, and prediction figures", ha="center", fontsize=10.8, color=PALETTE["slate"])
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_dashboard"])


if __name__ == "__main__":
    configure_theme()

    counts = flow_counts()
    plot_flowchart(counts)

    df = prepare_df()
    build_descriptive_table(df)
    norm = normality_checks(df)
    tests = analysis_a_b_c(df)
    fit = logistic_and_forest(df)

    save_summary(counts, norm, tests, fit)
    make_dashboard()

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
