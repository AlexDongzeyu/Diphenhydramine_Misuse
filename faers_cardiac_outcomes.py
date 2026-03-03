import argparse
from pathlib import Path
import re
import pandas as pd

# MedDRA Preferred Terms (hardcoded; no full MedDRA download required)
CARDIAC_OUTCOMES = [
    "Electrocardiogram QT prolonged",         # PT code: 10014387 — primary outcome
    "Electrocardiogram QT interval abnormal", # PT code: 10063748 — borderline QT
    "Long QT syndrome",                       # PT code: 10024803 — structural QT
    "Torsade de pointes",                     # PT code: 10044066 — malignant arrhythmia
    "Ventricular tachycardia",                # PT code: 10047302 — dangerous rhythm
    "Ventricular fibrillation",               # PT code: 10047281 — pre-terminal
    "Cardiac arrest",                         # PT code: 10007515 — full arrest/death
]

TIER_1_QT_SIGNAL = [
    "Electrocardiogram QT prolonged",
    "Electrocardiogram QT interval abnormal",
    "Long QT syndrome",
]

TIER_2_SERIOUS_CARDIAC_EVENT = [
    "Torsade de pointes",
    "Ventricular tachycardia",
    "Ventricular fibrillation",
    "Cardiac arrest",
]


def normalize_drug_name(name: str) -> str:
    text = str(name).lower().strip()
    text = text.replace("/", " ")
    text = re.sub(r"[^a-z0-9\s]", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def add_cardiac_outcomes(
    cohort_df: pd.DataFrame,
    reac_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Join FAERS REAC cardiac outcomes back to an already-filtered cohort (e.g.,
    diphenhydramine + adolescent cohort) using PRIMARYID.
    """
    if "PRIMARYID" not in cohort_df.columns:
        raise ValueError("cohort_df must contain PRIMARYID")
    if "PRIMARYID" not in reac_df.columns or "PT" not in reac_df.columns:
        raise ValueError("reac_df must contain PRIMARYID and PT")

    cohort_df = cohort_df.copy()
    reac_df = reac_df.copy()

    cohort_df["PRIMARYID"] = cohort_df["PRIMARYID"].astype(str).str.strip()
    reac_df["PRIMARYID"] = reac_df["PRIMARYID"].astype(str).str.strip()
    reac_df["PT"] = reac_df["PT"].astype(str).str.strip()

    # Exact filtering requested
    cardiac_cases = reac_df[reac_df["PT"].isin(CARDIAC_OUTCOMES)].copy()

    def classify_tier(pt: str) -> str:
        if pt in TIER_2_SERIOUS_CARDIAC_EVENT:
            return "Tier 2 — Serious cardiac event"
        if pt in TIER_1_QT_SIGNAL:
            return "Tier 1 — QT signal"
        return "Unknown"

    cardiac_cases["OUTCOME_TIER"] = cardiac_cases["PT"].map(classify_tier)

    # Case-level rollup: if a case has both tiers, Tier 2 dominates
    tier_rank = {
        "Tier 1 — QT signal": 1,
        "Tier 2 — Serious cardiac event": 2,
    }
    case_level = (
        cardiac_cases.assign(TIER_RANK=cardiac_cases["OUTCOME_TIER"].map(tier_rank))
        .sort_values(["PRIMARYID", "TIER_RANK"], ascending=[True, False])
        .drop_duplicates(subset=["PRIMARYID"], keep="first")
        [["PRIMARYID", "OUTCOME_TIER"]]
        .copy()
    )
    case_level["HAS_CARDIAC_OUTCOME"] = 1

    cohort_with_outcomes = cohort_df.merge(case_level, on="PRIMARYID", how="left")
    cohort_with_outcomes["HAS_CARDIAC_OUTCOME"] = (
        cohort_with_outcomes["HAS_CARDIAC_OUTCOME"].fillna(0).astype(int)
    )
    cohort_with_outcomes["OUTCOME_TIER"] = cohort_with_outcomes["OUTCOME_TIER"].fillna("None")

    return cohort_with_outcomes, cardiac_cases


def add_acb_scores(
    cohort_df: pd.DataFrame,
    drug_df: pd.DataFrame,
    acb_lookup_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Join ACB lookup to FAERS DRUG table and aggregate case-level ACB burden.
    Returns:
      1) cohort with ACB features joined on PRIMARYID
      2) matched DRUG rows with assigned ACB score
    """
    if "PRIMARYID" not in cohort_df.columns:
        raise ValueError("cohort_df must contain PRIMARYID")
    if "PRIMARYID" not in drug_df.columns:
        raise ValueError("drug_df must contain PRIMARYID")
    if "generic_name" not in acb_lookup_df.columns or "acb_score" not in acb_lookup_df.columns:
        raise ValueError("acb_lookup_df must contain generic_name and acb_score")

    candidate_name_cols = ["DRUGNAME", "PROD_AI", "drugname", "prod_ai", "generic_name"]
    drug_name_col = next((col for col in candidate_name_cols if col in drug_df.columns), None)
    if drug_name_col is None:
        raise ValueError("drug_df must contain one of: DRUGNAME, PROD_AI, generic_name")

    cohort_df = cohort_df.copy()
    drug_df = drug_df.copy()
    acb_lookup_df = acb_lookup_df.copy()

    cohort_df["PRIMARYID"] = cohort_df["PRIMARYID"].astype(str).str.strip()
    drug_df["PRIMARYID"] = drug_df["PRIMARYID"].astype(str).str.strip()

    drug_df["drug_name_raw"] = drug_df[drug_name_col].astype(str).str.strip()
    drug_df["drug_name_norm"] = drug_df["drug_name_raw"].map(normalize_drug_name)

    acb_lookup_df["generic_name"] = acb_lookup_df["generic_name"].astype(str).str.strip()
    acb_lookup_df["drug_name_norm"] = acb_lookup_df["generic_name"].map(normalize_drug_name)
    acb_lookup_df["acb_score"] = pd.to_numeric(acb_lookup_df["acb_score"], errors="coerce").fillna(0)

    acb_map = acb_lookup_df[["drug_name_norm", "generic_name", "acb_score"]].drop_duplicates(
        subset=["drug_name_norm"], keep="first"
    )

    drug_with_acb = drug_df.merge(acb_map, on="drug_name_norm", how="left")
    matched_drug_rows = drug_with_acb[drug_with_acb["acb_score"].notna()].copy()

    case_scores = (
        matched_drug_rows.groupby("PRIMARYID", as_index=False)
        .agg(
            acb_total_score=("acb_score", "sum"),
            acb_max_score=("acb_score", "max"),
            acb_matched_drug_count=("drug_name_norm", "nunique"),
        )
    )
    case_scores["acb_clinically_relevant"] = (case_scores["acb_total_score"] >= 3).astype(int)

    cohort_with_acb = cohort_df.merge(case_scores, on="PRIMARYID", how="left")
    cohort_with_acb["acb_total_score"] = cohort_with_acb["acb_total_score"].fillna(0)
    cohort_with_acb["acb_max_score"] = cohort_with_acb["acb_max_score"].fillna(0)
    cohort_with_acb["acb_matched_drug_count"] = cohort_with_acb["acb_matched_drug_count"].fillna(0).astype(int)
    cohort_with_acb["acb_clinically_relevant"] = cohort_with_acb["acb_clinically_relevant"].fillna(0).astype(int)

    return cohort_with_acb, matched_drug_rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Flag FAERS cardiac outcomes (hardcoded MedDRA PTs) from REAC and join "
            "back to a filtered cohort by PRIMARYID."
        )
    )
    parser.add_argument(
        "--cohort",
        default="data/faers_input/cohort.csv",
        help="Path to filtered cohort CSV (must include PRIMARYID)",
    )
    parser.add_argument(
        "--reac",
        default="data/faers_input/reac.csv",
        help="Path to REAC CSV (must include PRIMARYID and PT)",
    )
    parser.add_argument(
        "--out-cohort",
        default="cohort_with_cardiac_outcomes.csv",
        help="Output CSV path for cohort joined with cardiac outcome flags",
    )
    parser.add_argument(
        "--out-cardiac-reac",
        default="reac_cardiac_cases.csv",
        help="Output CSV path for REAC rows matching cardiac outcome PTs",
    )
    parser.add_argument(
        "--drug",
        default="data/faers_input/drug.csv",
        help="Path to DRUG CSV (must include PRIMARYID and DRUGNAME or PROD_AI)",
    )
    parser.add_argument(
        "--acb-lookup",
        default="data/lookups/acb_lookup.csv",
        help="Path to ACB lookup CSV with columns generic_name, acb_score",
    )
    parser.add_argument(
        "--out-cohort-acb",
        default="cohort_with_cardiac_and_acb.csv",
        help="Output CSV path for cohort with cardiac outcomes + ACB scores",
    )
    parser.add_argument(
        "--out-drug-acb-matches",
        default="drug_acb_matches.csv",
        help="Output CSV path for DRUG rows matched to ACB lookup",
    )

    args = parser.parse_args()

    cohort_path = Path(args.cohort)
    reac_path = Path(args.reac)
    drug_path = Path(args.drug)
    acb_lookup_path = Path(args.acb_lookup)
    missing = [str(p) for p in (cohort_path, reac_path, drug_path, acb_lookup_path) if not p.exists()]
    if missing:
        print("Missing required input file(s):")
        for path in missing:
            print(f"- {path}")
        print("\nNo FAERS download is required yet.")
        print("Leave your files in these placeholders when ready:")
        print("- data/faers_input/cohort.csv")
        print("- data/faers_input/reac.csv")
        print("- data/faers_input/drug.csv")
        print("- data/lookups/acb_lookup.csv")
        print("\nThen run:")
        print("python faers_cardiac_outcomes.py")
        return

    cohort_df = pd.read_csv(args.cohort, dtype=str, low_memory=False)
    reac_df = pd.read_csv(args.reac, dtype=str, low_memory=False)
    drug_df = pd.read_csv(args.drug, dtype=str, low_memory=False)
    acb_lookup_df = pd.read_csv(args.acb_lookup, dtype=str, low_memory=False)

    cohort_with_outcomes, cardiac_cases = add_cardiac_outcomes(cohort_df, reac_df)
    cohort_with_cardiac_and_acb, drug_acb_matches = add_acb_scores(
        cohort_with_outcomes,
        drug_df,
        acb_lookup_df,
    )

    cohort_with_outcomes.to_csv(args.out_cohort, index=False)
    cardiac_cases.to_csv(args.out_cardiac_reac, index=False)
    cohort_with_cardiac_and_acb.to_csv(args.out_cohort_acb, index=False)
    drug_acb_matches.to_csv(args.out_drug_acb_matches, index=False)

    summary = cohort_with_cardiac_and_acb["OUTCOME_TIER"].value_counts(dropna=False)
    acb_summary = cohort_with_cardiac_and_acb["acb_total_score"].describe()
    print("Saved:")
    print(f"- {args.out_cohort}")
    print(f"- {args.out_cardiac_reac}")
    print(f"- {args.out_cohort_acb}")
    print(f"- {args.out_drug_acb_matches}")
    print("\nOutcome tier counts:")
    print(summary.to_string())
    print("\nACB total score summary:")
    print(acb_summary.to_string())


if __name__ == "__main__":
    main()
