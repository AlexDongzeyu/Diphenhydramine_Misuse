"""Build publication-ready tables/figures from the final FAERS cohort table.

This script is intentionally end-to-end: it loads the final analysis dataset,
runs descriptive and inferential statistics, fits a logistic model, and writes
all outputs into 06_analysis/tables and 06_analysis/figures.
"""

import json
import shutil
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from matplotlib.patches import FancyBboxPatch
from matplotlib import patheffects
from scipy.stats import shapiro, mannwhitneyu, kruskal, spearmanr
import scikit_posthocs as sp
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibrationDisplay

try:
    import shap
except Exception:  # pragma: no cover - optional dependency for explainability export
    shap = None


ROOT = Path(__file__).resolve().parents[1]
FINAL_CSV = ROOT / "05_final" / "cohort_analysis.csv"
FILTERED = ROOT / "03_filtered"
OUT_DIR = ROOT / "06_analysis"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
FINAL_FIG_DIR = FIG_DIR / "final"
SUBMISSION_FIG_DIR = FIG_DIR / "submission"

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
    "table_cv_predictions": "model_predictions_cv.csv",
    "table_cv_metrics": "model_cv_metrics.csv",
    "table_modeling_n": "modeling_sample_size.csv",
    "table_p_values_bh": "p_values_bh_python_core.csv",
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


for path in [FIG_DIR, TAB_DIR, FINAL_FIG_DIR, SUBMISSION_FIG_DIR]:
    path.mkdir(parents=True, exist_ok=True)


PALETTE = {
    "blue": "#4C78A8",
    "orange": "#D28E3D",
    "gray": "#7F7F7F",
    "light_blue": "#DCE6F2",
    "light_orange": "#F2E1C8",
    "light_gray": "#E6E6E6",
    "black": "#000000",
}

HUMAN_LABELS = {
    "AGE": "Patient age (years)",
    "total_acb_with_dph": "Total anticholinergic burden (including diphenhydramine)",
    "total_acb_codrugs_only": "Total anticholinergic burden from co-medications",
    "n_codrugs": "Number of concomitant co-medications",
    "max_severity": "Maximum reported case severity",
    "cardiac_any": "Cardiac toxicity outcome (yes/no)",
    "pre_post_warning": "Regulatory warning period (pre/post)",
    "age_group": "Age group",
    "Intercept": "Model intercept",
}

# Global journal-style defaults: plain white canvas, black text, no embedded titles.
sns.set_theme(
    style="white",
    rc={
        "axes.facecolor": "#ffffff",
        "figure.facecolor": "#ffffff",
        "font.family": "Arial",
        "text.color": "black",
        "axes.labelcolor": "black",
        "axes.edgecolor": "black",
        "axes.linewidth": 0.8,
        "axes.labelsize": 11,
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "grid.color": "#E6E6E6",
        "grid.linewidth": 0.5,
    },
)


def style_axes(ax: plt.Axes, *, grid: bool = False, grid_axis: str = "both") -> None:
    """Apply the conservative journal axis treatment used across manuscript figures."""
    ax.set_title("")
    ax.set_facecolor("white")
    ax.tick_params(colors="black", width=0.8)
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(0.8)
    if grid:
        ax.grid(True, axis=grid_axis, color="#E6E6E6", linewidth=0.5)
    else:
        ax.grid(False)


def save_figure(fig: plt.Figure, path: Path, dpi: int = 600) -> None:
    """Save figures at journal resolution with a white background."""
    tmp_path = path.with_name(f".{path.stem}.tmp{path.suffix}")
    fig.savefig(tmp_path, dpi=dpi, bbox_inches="tight", pad_inches=0.16, facecolor="white", edgecolor="none")
    for attempt in range(12):
        try:
            tmp_path.replace(path)
            return
        except PermissionError:
            if attempt == 0:
                try:
                    path.unlink(missing_ok=True)
                except PermissionError:
                    pass
            time.sleep(0.5)
    tmp_path.replace(path)


def human_label(name: str) -> str:
    """Convert code-style column names to publication-ready labels."""
    if name in HUMAN_LABELS:
        return HUMAN_LABELS[name]
    return name.replace("_", " ").strip().capitalize()


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
        ax.hist(s, bins=15, color=PALETTE["light_blue"], edgecolor="black", linewidth=0.4)
        ax.set_xlabel(human_label(col))
        ax.set_ylabel("Count")
        style_axes(ax)
        fig.tight_layout()
        save_figure(fig, FIG_DIR / HISTOGRAM_NAMES[col])
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(5, 5))
        # Use scipy's probability plot to visualize normality assumptions.
        from scipy import stats
        stats.probplot(s, dist="norm", plot=ax)
        style_axes(ax)
        fig.tight_layout()
        save_figure(fig, FIG_DIR / QQ_NAMES[col])
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

    groups = [g.dropna().values for _, g in df.groupby("age_group")["total_acb_codrugs_only"]]
    k_stat, p_kw = kruskal(*groups)
    rows.append({"analysis": "B_kruskal_acb_vs_age_group", "statistic": k_stat, "p_value": p_kw})

    # Post-hoc Dunn test is used after Kruskal-Wallis for pairwise age-group comparisons.
    dunn_df = sp.posthoc_dunn(df, val_col="total_acb_codrugs_only", group_col="age_group", p_adjust="bonferroni")
    dunn_df.to_csv(TAB_DIR / FILE_NAMES["table_posthoc"])

    df_age = df.loc[df["age_group"].isin(["13-15", "16-17", "18-19"])].copy()
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    order = ["13-15", "16-17", "18-19"]
    sns.violinplot(
        data=df_age,
        x="age_group",
        y="total_acb_codrugs_only",
        order=order,
        hue="age_group",
        legend=False,
        palette=[PALETTE["light_blue"], PALETTE["light_orange"], PALETTE["light_gray"]],
        inner=None,
        linewidth=0.7,
        cut=0,
        ax=ax,
    )
    sns.boxplot(
        data=df_age,
        x="age_group",
        y="total_acb_codrugs_only",
        order=order,
        hue="age_group",
        legend=False,
        palette=[PALETTE["light_blue"], PALETTE["light_orange"], PALETTE["light_gray"]],
        width=0.28,
        linewidth=0.7,
        showfliers=False,
        ax=ax,
    )
    sns.stripplot(
        data=df_age,
        x="age_group",
        y="total_acb_codrugs_only",
        order=order,
        color=PALETTE["gray"],
        alpha=0.16,
        size=2.6,
        jitter=0.2,
        ax=ax,
    )
    y_cap = df_age["total_acb_codrugs_only"].quantile(0.99)
    if not np.isfinite(y_cap) or y_cap <= 0:
        y_cap = np.nanmax(df_age["total_acb_codrugs_only"]) if df_age["total_acb_codrugs_only"].notna().any() else 1
    ax.set_ylim(0, y_cap * 1.20 if y_cap > 0 else 1.5)

    level_map = {lvl: idx for idx, lvl in enumerate(order)}
    pair_y = y_cap * 1.05 if y_cap > 0 else 1
    # Draw significance brackets only for statistically significant pairwise contrasts.
    for a in order:
        for b in order:
            if a >= b:
                continue
            pval = dunn_df.loc[a, b]
            if pd.notna(pval) and pval < 0.05:
                xa, xb = level_map[a], level_map[b]
                tick = y_cap * 0.025 if y_cap > 0 else 0.04
                ax.plot([xa, xa, xb, xb], [pair_y, pair_y + tick, pair_y + tick, pair_y], color="black", lw=0.8)
                ax.text((xa + xb) / 2, pair_y + tick * 1.15, p_stars(float(pval)), ha="center", va="bottom", color="black", fontsize=10)
                pair_y += y_cap * 0.075 if y_cap > 0 else 0.08

    ax.set_xlabel(human_label("age_group"))
    ax.set_ylabel(human_label("total_acb_codrugs_only"))
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_age_group"])
    plt.close(fig)

    rho, p_sp = spearmanr(df["total_acb_codrugs_only"], df["max_severity"], nan_policy="omit")
    rows.append({"analysis": "C_spearman_acb_vs_severity", "statistic": rho, "p_value": p_sp})

    fig, ax = plt.subplots(figsize=(7, 5.2))
    sns.regplot(
        data=df,
        x="total_acb_codrugs_only",
        y="max_severity",
        lowess=True,
        scatter_kws={"alpha": 0.45, "s": 22, "color": PALETTE["blue"]},
        line_kws={"color": PALETTE["orange"], "lw": 2.5},
        ax=ax,
    )
    ax.set_xlabel(human_label("total_acb_codrugs_only"))
    ax.set_ylabel(human_label("max_severity"))
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_severity"])
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


def run_cv_models(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Build strict out-of-fold probabilities, calibration plots, and SHAP export artifacts."""
    features = ["AGE", "total_acb_codrugs_only", "n_codrugs"]
    model_cols = features + ["cardiac_any", "PRIMARYID"]

    # Explicit missing-data handling avoids sklearn crashes and documents final N.
    df_model = df.loc[:, model_cols].dropna().copy()
    dropped_rows = int(df.shape[0] - df_model.shape[0])
    modeling_n = {
        "n_total": int(df.shape[0]),
        "n_model": int(df_model.shape[0]),
        "n_dropped_missing": dropped_rows,
    }
    pd.DataFrame([modeling_n]).to_csv(TAB_DIR / FILE_NAMES["table_modeling_n"], index=False)

    X_df = df_model.loc[:, features].copy()
    X = X_df.values
    y = df_model["cardiac_any"].astype(int).values

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    rf_model = RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=42)
    lr_model = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(class_weight="balanced", max_iter=1000, random_state=42)),
    ])

    # Strict out-of-fold probabilities prevent data leakage in performance estimates.
    lr_prob = cross_val_predict(lr_model, X, y, cv=cv, method="predict_proba")[:, 1]
    rf_prob = cross_val_predict(rf_model, X, y, cv=cv, method="predict_proba")[:, 1]

    pred_df = pd.DataFrame({
        "PRIMARYID": df_model["PRIMARYID"].astype(str).values,
        "cardiac_any": y,
        "lr_prob": lr_prob,
        "rf_prob": rf_prob,
    })
    pred_df.to_csv(TAB_DIR / FILE_NAMES["table_cv_predictions"], index=False)

    auc_lr = roc_auc_score(y, lr_prob)
    auc_rf = roc_auc_score(y, rf_prob)
    metrics = pd.DataFrame([
        {"model": "logistic_regression_cv", "auc": float(auc_lr), "n": int(len(y))},
        {"model": "random_forest_cv", "auc": float(auc_rf), "n": int(len(y))},
    ])
    metrics.to_csv(TAB_DIR / FILE_NAMES["table_cv_metrics"], index=False)

    fig, ax = plt.subplots(figsize=(7.0, 5.5))
    CalibrationDisplay.from_predictions(y, lr_prob, n_bins=5, name="Logistic (OOF)", ax=ax)
    CalibrationDisplay.from_predictions(y, rf_prob, n_bins=5, name="Random Forest (OOF)", ax=ax)
    style_axes(ax, grid=True, grid_axis="both")
    legend = ax.get_legend()
    if legend is not None:
        legend.set_frame_on(False)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "calibration_curve_cv.png")
    plt.close(fig)

    # This model is strictly for SHAP interpretability. All AUC metrics come from CV above.
    if shap is not None:
        rf_model.fit(X_df, y)
        explainer = shap.TreeExplainer(rf_model)
        shap_values = explainer.shap_values(X_df)
        shap_for_plot = None

        # Explicit class handling avoids plotting the wrong class in binary outputs.
        if isinstance(shap_values, list) and len(shap_values) > 1:
            shap_for_plot = shap_values[1]
        elif isinstance(shap_values, np.ndarray) and shap_values.ndim == 3:
            # Newer SHAP versions may return (n_samples, n_features, n_classes).
            shap_for_plot = shap_values[:, :, 1]
        elif isinstance(shap_values, np.ndarray) and shap_values.ndim == 2:
            if getattr(rf_model, "n_classes_", None) == 2 and shap_values.shape[1] == X_df.shape[1]:
                # For binary TreeExplainer, 2D output corresponds to the positive class.
                shap_for_plot = shap_values
            else:
                print("Unexpected 2D SHAP output shape; skipping SHAP beeswarm export.")
        else:
            print("Unexpected SHAP output format; skipping SHAP beeswarm export.")

        if shap_for_plot is not None:
            # Display publication-ready feature names on SHAP outputs.
            X_display = X_df.rename(columns={c: human_label(c) for c in X_df.columns})
            plt.figure(figsize=(11, 6.5))
            shap.summary_plot(shap_for_plot, X_display, show=False, plot_size=(11, 6.5))
            fig_shap = plt.gcf()
            fig_shap.subplots_adjust(left=0.28, right=0.98, top=0.93, bottom=0.18)
            save_figure(fig_shap, FIG_DIR / "shap_beeswarm.png")
            plt.close()
    else:
        print("SHAP is not installed; skipping SHAP beeswarm export.")

    return pred_df, metrics, modeling_n


def plot_dual_roc(cv_pred_df: pd.DataFrame) -> None:
    """Render the manuscript ROC comparison from strict out-of-fold predictions."""
    roc_df = cv_pred_df.dropna(subset=["cardiac_any", "lr_prob", "rf_prob"]).copy()
    y = roc_df["cardiac_any"].astype(int)
    auc_lr = roc_auc_score(y, roc_df["lr_prob"])
    auc_rf = roc_auc_score(y, roc_df["rf_prob"])
    fpr_lr, tpr_lr, _ = roc_curve(y, roc_df["lr_prob"])
    fpr_rf, tpr_rf, _ = roc_curve(y, roc_df["rf_prob"])

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(fpr_rf, tpr_rf, color=PALETTE["blue"], lw=1.8, label=f"Random forest (AUC={auc_rf:.3f})")
    ax.plot(fpr_lr, tpr_lr, color=PALETTE["orange"], lw=1.8, label=f"Logistic regression (AUC={auc_lr:.3f})")
    ax.plot([0, 1], [0, 1], color=PALETTE["gray"], ls="--", lw=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("1 - Specificity")
    ax.set_ylabel("Sensitivity")
    style_axes(ax)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "roc_dual.png")
    save_figure(fig, FIG_DIR / "roc_dual.pdf")
    plt.close(fig)


def plot_sex_stratified_acb(df: pd.DataFrame) -> None:
    """Render the sex-stratified supplementary burden figure when sex is available."""
    demo_path = FILTERED / FILE_NAMES["demo"]
    if not demo_path.exists():
        print("Sex stratification skipped: teen_demo_records.csv not found.")
        return

    demo = pd.read_csv(demo_path, dtype=str, usecols=["PRIMARYID", "GNDR_COD"])
    sex_df = (
        df.assign(PRIMARYID=df["PRIMARYID"].astype(str))
        .merge(demo.assign(PRIMARYID=demo["PRIMARYID"].astype(str)), on="PRIMARYID", how="left")
        .assign(SEX=lambda d: d["GNDR_COD"].str.upper())
    )
    sex_df = sex_df.loc[sex_df["SEX"].isin(["F", "M"])].copy()
    if sex_df["SEX"].nunique() < 2:
        print("Sex stratification skipped: need both F and M groups.")
        return

    y_cap = sex_df["total_acb_codrugs_only"].quantile(0.99)
    if not np.isfinite(y_cap) or y_cap <= 0:
        y_cap = sex_df["total_acb_codrugs_only"].max()

    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    order = ["F", "M"]
    sns.violinplot(
        data=sex_df,
        x="SEX",
        y="total_acb_codrugs_only",
        order=order,
        hue="SEX",
        legend=False,
        palette=[PALETTE["light_orange"], PALETTE["light_blue"]],
        inner=None,
        linewidth=0.7,
        cut=0,
        ax=ax,
    )
    sns.boxplot(
        data=sex_df,
        x="SEX",
        y="total_acb_codrugs_only",
        order=order,
        hue="SEX",
        legend=False,
        palette=[PALETTE["light_orange"], PALETTE["light_blue"]],
        width=0.22,
        linewidth=0.7,
        showfliers=False,
        ax=ax,
    )
    sns.stripplot(
        data=sex_df,
        x="SEX",
        y="total_acb_codrugs_only",
        order=order,
        color=PALETTE["gray"],
        alpha=0.14,
        size=2.5,
        jitter=0.14,
        ax=ax,
    )
    ax.set_ylim(0, y_cap * 1.08 if y_cap > 0 else 1.5)
    ax.set_xlabel("Sex")
    ax.set_ylabel(human_label("total_acb_codrugs_only"))
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "sex_stratified_acb_violin.png")
    plt.close(fig)


def plot_drug_cooccurrence_network() -> None:
    """Render a conservative drug co-occurrence network without embedded titles."""
    drug_path = ROOT / "04_processed" / "dph_drug_normalized.csv"
    if not drug_path.exists():
        print("Network analysis skipped: dph_drug_normalized.csv not found.")
        return

    drug_df = pd.read_csv(drug_path, dtype=str, usecols=["PRIMARYID", "DRUGNAME", "generic_name"])
    drug_df["drug_name"] = drug_df["generic_name"].fillna("").str.strip().str.lower()
    fallback = drug_df["drug_name"].eq("")
    drug_df.loc[fallback, "drug_name"] = drug_df.loc[fallback, "DRUGNAME"].fillna("").str.strip().str.lower()
    drug_df["drug_name"] = drug_df["drug_name"].replace(
        {
            "benadryl": "diphenhydramine",
            "tylenol": "acetaminophen",
            "advil": "ibuprofen",
            "motrin": "ibuprofen",
            "zofran": "ondansetron",
        }
    )
    drug_df = drug_df.loc[drug_df["PRIMARYID"].notna() & drug_df["drug_name"].ne(""), ["PRIMARYID", "drug_name"]].drop_duplicates()

    pair_counts: dict[tuple[str, str], int] = {}
    for drugs in drug_df.groupby("PRIMARYID")["drug_name"]:
        unique_drugs = sorted(set(drugs[1]))
        for a, b in combinations(unique_drugs, 2):
            pair_counts[(a, b)] = pair_counts.get((a, b), 0) + 1
    if not pair_counts:
        print("Network analysis skipped: no co-occurring drug pairs found.")
        return

    weights = np.array(list(pair_counts.values()), dtype=float)
    threshold = max(2, int(np.floor(np.quantile(weights, 0.75))))
    graph = nx.Graph()
    for (a, b), weight in pair_counts.items():
        if weight >= threshold:
            graph.add_edge(a, b, weight=weight)
    if graph.number_of_edges() == 0:
        print("Network analysis skipped: no edges met the threshold.")
        return

    largest_nodes = max(nx.connected_components(graph), key=len)
    graph = graph.subgraph(largest_nodes).copy()
    strength = dict(graph.degree(weight="weight"))
    keep_nodes = sorted(strength, key=strength.get, reverse=True)[:55]
    graph = graph.subgraph(keep_nodes).copy()
    strength = dict(graph.degree(weight="weight"))

    pos = nx.spring_layout(graph, seed=42, weight="weight", k=0.55)
    edge_weights = np.array([graph.edges[e]["weight"] for e in graph.edges], dtype=float)
    edge_widths = 0.4 + 1.8 * (edge_weights - edge_weights.min()) / (np.ptp(edge_weights) or 1)
    node_sizes = [80 + 520 * strength[n] / max(strength.values()) for n in graph.nodes]

    fig, ax = plt.subplots(figsize=(13, 9))
    ax.set_facecolor("white")
    nx.draw_networkx_edges(graph, pos, ax=ax, width=edge_widths, edge_color="#9A9A9A", alpha=0.45)
    nx.draw_networkx_nodes(graph, pos, ax=ax, node_size=node_sizes, node_color=PALETTE["blue"], edgecolors="black", linewidths=0.4, alpha=0.92)
    nx.draw_networkx_labels(
        graph,
        pos,
        ax=ax,
        font_size=7,
        font_color="black",
        bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "#BDBDBD", "linewidth": 0.3, "alpha": 0.9},
    )
    ax.axis("off")
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "drug_cooccurrence_network.png")
    plt.close(fig)


def save_bh_table(tests: pd.DataFrame) -> pd.DataFrame:
    """Create a BH table for Python core tests only; use R BH table as study-level primary."""

    def bh_adjust(p_values: pd.Series) -> pd.Series:
        vals = p_values.astype(float).values
        n = len(vals)
        order = np.argsort(vals)
        ranked = vals[order]
        adjusted = np.empty(n, dtype=float)
        min_so_far = 1.0
        for i in range(n - 1, -1, -1):
            rank = i + 1
            val = ranked[i] * n / rank
            min_so_far = min(min_so_far, val)
            adjusted[i] = min_so_far
        out = np.empty(n, dtype=float)
        out[order] = np.clip(adjusted, 0.0, 1.0)
        return pd.Series(out, index=p_values.index)

    p_tbl = tests.loc[:, ["analysis", "p_value"]].copy()
    p_tbl = p_tbl.rename(columns={"analysis": "test_name", "p_value": "p_raw"})
    p_tbl["p_bh"] = bh_adjust(p_tbl["p_raw"])
    p_tbl["bh_significant_0_05"] = p_tbl["p_bh"] < 0.05
    p_tbl.to_csv(TAB_DIR / FILE_NAMES["table_p_values_bh"], index=False)
    return p_tbl


def logistic_and_forest(df: pd.DataFrame, cv_pred_df: pd.DataFrame | None = None) -> pd.DataFrame:
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
        ax.hlines(y, forest_df["CI_low"], forest_df["CI_high"], color=PALETTE["blue"], lw=1.6)
    ax.scatter(forest_df["OR"], y, color=PALETTE["orange"], s=48, zorder=3, edgecolor="black", linewidth=0.4)
    ax.axvline(1.0, color="black", ls="--", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([human_label(term) for term in forest_df["term"]])
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio (log scale)")
    ax.set_ylabel("Clinical predictor")
    style_axes(ax, grid=True, grid_axis="x")
    fig.subplots_adjust(left=0.38, right=0.98, top=0.96, bottom=0.18)
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_forest"])
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_forest"].replace(".png", ".pdf"))
    plt.close(fig)

    # Prefer strict out-of-fold probabilities for ROC plotting to avoid in-sample leakage.
    if cv_pred_df is not None and {"cardiac_any", "lr_prob"}.issubset(cv_pred_df.columns):
        roc_df = cv_pred_df.dropna(subset=["cardiac_any", "lr_prob"]).copy()
        score_col = "lr_prob"
        auc = roc_auc_score(roc_df["cardiac_any"], roc_df[score_col])
    else:
        roc_df = pred_df.dropna(subset=["cardiac_any", "pred_prob"]).copy()
        score_col = "pred_prob"
        auc = fit_full["auc"]

    fpr, tpr, _ = roc_curve(roc_df["cardiac_any"], roc_df[score_col])
    fig, ax = plt.subplots(figsize=(6.8, 6.0))
    ax.plot(fpr, tpr, color=PALETTE["blue"], lw=1.8, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], color=PALETTE["gray"], ls="--", lw=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    style_axes(ax)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_roc"])
    plt.close(fig)

    pd.DataFrame(fit_rows).to_csv(TAB_DIR / FILE_NAMES["table_fit"], index=False)
    return pd.DataFrame(fit_rows)


def plot_flowchart(counts: dict[str, int]) -> None:
    """Render a simple cohort flowchart from precomputed count checkpoints."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")

    boxes = [
        (0.1, 0.80, 0.8, 0.12, f"DEMO teen reports (unique case IDs): {counts['DEMO_teens_unique']:,}"),
        (0.1, 0.62, 0.8, 0.12, f"DRUG teen reports (unique case IDs): {counts['DRUG_teens_unique']:,}"),
        (0.1, 0.44, 0.8, 0.12, f"Confirmed DPH cohort: {counts['DPH_confirmed_unique']:,}"),
        (0.1, 0.26, 0.8, 0.12, f"Final analysis records: {counts['Final_rows']:,} (unique case IDs: {counts['Final_unique_PRIMARYID']:,})"),
    ]

    for (x, y, w, h, text) in boxes:
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.02", edgecolor="black", facecolor="white", linewidth=0.8)
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10, color="black")

    for y1, y2 in [(0.80, 0.74), (0.62, 0.56), (0.44, 0.38)]:
        ax.annotate("", xy=(0.5, y2), xytext=(0.5, y1), arrowprops=dict(arrowstyle="->", lw=0.9, color="black"))

    fig.tight_layout()
    save_figure(fig, FIG_DIR / FILE_NAMES["figure_flow"])
    plt.close(fig)


def save_summary(
    counts: dict[str, int],
    norm: pd.DataFrame,
    tests: pd.DataFrame,
    fit: pd.DataFrame,
    cv_metrics: pd.DataFrame,
    modeling_n: dict[str, int],
) -> None:
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

    lines.append("")
    lines.append("## CV model performance (leak-free out-of-fold)")
    lines.append(
        f"- Modeling sample: n_model={modeling_n['n_model']} / n_total={modeling_n['n_total']} "
        f"(dropped_missing={modeling_n['n_dropped_missing']})"
    )
    for row in cv_metrics.itertuples(index=False):
        lines.append(f"- {row.model}: AUC={row.auc:.3f}")

    (OUT_DIR / FILE_NAMES["summary"]).write_text("\n".join(lines), encoding="utf-8")


def make_dashboard() -> None:
    """Assemble a 2x2 dashboard image from the core generated figures."""
    fig = plt.figure(figsize=(14, 10), facecolor="white", constrained_layout=True)
    gs = fig.add_gridspec(2, 2, wspace=0.25, hspace=0.3)

    files = [
        FIG_DIR / FILE_NAMES["figure_flow"],
        FIG_DIR / FILE_NAMES["figure_age_group"],
        FIG_DIR / FILE_NAMES["figure_forest"],
        FIG_DIR / FILE_NAMES["figure_roc"],
    ]
    for i, img_path in enumerate(files):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        img = plt.imread(img_path)
        ax.imshow(img)
        ax.axis("off")

    save_figure(fig, FIG_DIR / FILE_NAMES["figure_dashboard"])
    plt.close(fig)


def export_submission_figures() -> None:
    """Copy manuscript figures into journal submission filenames."""
    final_names = [
        FILE_NAMES["figure_flow"],
        FILE_NAMES["figure_age_group"],
        FILE_NAMES["figure_severity"],
        FILE_NAMES["figure_forest"],
        FILE_NAMES["figure_roc"],
        "roc_dual.png",
        "calibration_curve_cv.png",
        "shap_beeswarm.png",
        "drug_cooccurrence_network.png",
        "sex_stratified_acb_violin.png",
    ]
    for name in final_names:
        src = FIG_DIR / name
        if src.exists():
            shutil.copy2(src, FINAL_FIG_DIR / name)

    mapping = {
        "Figure_1.png": FIG_DIR / FILE_NAMES["figure_flow"],
        "Figure_2.png": FIG_DIR / HISTOGRAM_NAMES["AGE"],
        "Figure_3.png": FIG_DIR / FILE_NAMES["figure_age_group"],
        "Figure_4.png": FIG_DIR / FILE_NAMES["figure_severity"],
        "Figure_5.png": FIG_DIR / "roc_dual.png",
        "Figure_5.pdf": FIG_DIR / "roc_dual.pdf",
        "Figure_6.png": FIG_DIR / FILE_NAMES["figure_forest"],
        "Figure_6.pdf": FIG_DIR / FILE_NAMES["figure_forest"].replace(".png", ".pdf"),
    }
    for filename, src in mapping.items():
        if src.exists():
            shutil.copy2(src, SUBMISSION_FIG_DIR / filename)


if __name__ == "__main__":
    # 1) Build flow/count context and load analysis-ready dataset.
    counts = flow_counts()
    plot_flowchart(counts)

    df = prepare_df()

    # 2) Generate core descriptive, inferential, and predictive outputs.
    build_descriptive_table(df)
    norm = normality_checks(df)
    tests = analysis_a_b_c(df)
    cv_pred_df, cv_metrics, modeling_n = run_cv_models(df)
    plot_dual_roc(cv_pred_df)
    fit = logistic_and_forest(df, cv_pred_df=cv_pred_df)
    save_bh_table(tests)
    plot_sex_stratified_acb(df)
    plot_drug_cooccurrence_network()

    # 3) Write summaries and composite dashboard for reporting.
    save_summary(counts, norm, tests, fit, cv_metrics, modeling_n)
    make_dashboard()
    export_submission_figures()

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
