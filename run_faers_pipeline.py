import argparse
import csv
import re
from pathlib import Path

import pandas as pd

TABLES = ("DEMO", "DRUG", "REAC", "OUTC")
DPH_TERMS = {
    "DIPHENHYDRAMINE",
    "DIPHENHYDRAMINE HCL",
    "DIPHENHYDRAMINE HYDROCHLORIDE",
    "DIPHENHYDRAMINE CITRATE",
    "BENADRYL",
    "BENADRYL ALLERGY",
    "NYTOL",
    "SOMINEX",
    "SIMPLY SLEEP",
    "ZZZ QUIL",
    "DPH",
}
TIER1 = {
    "Electrocardiogram QT prolonged",
    "Electrocardiogram QT interval abnormal",
    "Long QT syndrome",
}
TIER2 = {
    "Torsade de pointes",
    "Ventricular tachycardia",
    "Ventricular fibrillation",
    "Cardiac arrest",
}


def log(msg: str) -> None:
    print(f"[PIPELINE] {msg}", flush=True)


def ensure_dirs(base: Path) -> dict[str, Path]:
    paths = {
        "combined": base / "02_combined",
        "filtered": base / "03_filtered",
        "processed": base / "04_processed",
        "final": base / "05_final",
        "raw_acb": base / "01_raw" / "acb",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def collect_quarter_files(extracted_root: Path, table_prefix: str) -> list[Path]:
    candidates: list[Path] = []
    for quarter_dir in sorted([d for d in extracted_root.iterdir() if d.is_dir()]):
        for file_path in sorted(quarter_dir.iterdir()):
            name = file_path.name.upper()
            if name.startswith(table_prefix) and name.endswith(".TXT"):
                candidates.append(file_path)
    return candidates


def clean_columns(columns: list[str]) -> list[str]:
    cleaned: list[str] = []
    for col in columns:
        if col is None:
            continue
        col_name = str(col).strip()
        if not col_name:
            continue
        if col_name.upper().startswith("UNNAMED"):
            continue
        cleaned.append(col_name.upper())
    return cleaned


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(col).strip().upper() for col in out.columns]
    return out


def combine_table(files: list[Path], out_csv: Path, chunksize: int) -> list[str]:
    if not files:
        raise RuntimeError(f"No source files found for {out_csv.name}")

    all_cols: list[str] = []
    seen_cols = set()

    for file_path in files:
        header_df = pd.read_csv(
            file_path,
            sep="$",
            dtype=str,
            encoding="latin-1",
            on_bad_lines="skip",
            engine="python",
            nrows=0,
        )
        cols = clean_columns(header_df.columns.tolist())
        for col in cols:
            if col not in seen_cols:
                seen_cols.add(col)
                all_cols.append(col)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if out_csv.exists():
        out_csv.unlink()

    wrote_header = False
    for idx, file_path in enumerate(files, start=1):
        log(f"Combining {out_csv.name}: {idx}/{len(files)} -> {file_path.name}")
        for chunk in pd.read_csv(
            file_path,
            sep="$",
            dtype=str,
            encoding="latin-1",
            on_bad_lines="skip",
            engine="python",
            chunksize=chunksize,
        ):
            chunk = normalize_columns(chunk)
            drop_cols = [c for c in chunk.columns if c.startswith("UNNAMED") or c.strip() == ""]
            if drop_cols:
                chunk = chunk.drop(columns=drop_cols, errors="ignore")
            chunk = chunk.reindex(columns=all_cols)
            chunk.to_csv(out_csv, mode="a", index=False, header=not wrote_header)
            wrote_header = True

    return all_cols


def to_numeric_id(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(-1).astype("int64")


def dedupe_demo(demo_all: Path, demo_deduped: Path, chunksize: int) -> None:
    log("Phase 4: Deduplicating DEMO by CASEID with max PRIMARYID")
    max_primary_by_case: dict[str, int] = {}

    for idx, chunk in enumerate(pd.read_csv(demo_all, dtype=str, chunksize=chunksize, low_memory=False), start=1):
        chunk = normalize_columns(chunk)
        c = chunk[["CASEID", "PRIMARYID"]].copy()
        c["CASEID"] = c["CASEID"].astype(str).str.strip()
        c["PRIMARYID_NUM"] = to_numeric_id(c["PRIMARYID"])
        grouped = c.groupby("CASEID", as_index=False)["PRIMARYID_NUM"].max()
        for row in grouped.itertuples(index=False):
            caseid = row.CASEID
            primary_num = int(row.PRIMARYID_NUM)
            prev = max_primary_by_case.get(caseid)
            if prev is None or primary_num > prev:
                max_primary_by_case[caseid] = primary_num
        if idx % 20 == 0:
            log(f"DEMO pass1 chunks processed: {idx}")

    if demo_deduped.exists():
        demo_deduped.unlink()

    wrote_header = False
    for idx, chunk in enumerate(pd.read_csv(demo_all, dtype=str, chunksize=chunksize, low_memory=False), start=1):
        chunk = normalize_columns(chunk)
        chunk["CASEID"] = chunk["CASEID"].astype(str).str.strip()
        chunk["PRIMARYID_NUM"] = to_numeric_id(chunk["PRIMARYID"])
        keep_max = chunk["CASEID"].map(max_primary_by_case).fillna(-1).astype("int64")
        deduped = chunk[chunk["PRIMARYID_NUM"] == keep_max].drop(columns=["PRIMARYID_NUM"])
        if not deduped.empty:
            deduped.to_csv(demo_deduped, mode="a", index=False, header=not wrote_header)
            wrote_header = True
        if idx % 20 == 0:
            log(f"DEMO pass2 chunks processed: {idx}")


def filter_demo_teens(demo_deduped: Path, demo_teens: Path, teen_ids_txt: Path, chunksize: int) -> set[str]:
    log("Phase 5: Filtering DEMO to US adolescents with country unknown retained")
    if demo_teens.exists():
        demo_teens.unlink()

    teen_ids: set[str] = set()
    wrote_header = False

    for idx, chunk in enumerate(pd.read_csv(demo_deduped, dtype=str, chunksize=chunksize, low_memory=False), start=1):
        c = normalize_columns(chunk)
        c["OCCR_COUNTRY"] = c.get("OCCR_COUNTRY", "").astype(str).str.strip()
        c["AGE_COD"] = c.get("AGE_COD", "").astype(str).str.strip().str.upper()
        c["AGE_NUM"] = pd.to_numeric(c.get("AGE", ""), errors="coerce")

        country_unknown = c["OCCR_COUNTRY"].eq("") | c["OCCR_COUNTRY"].str.upper().eq("NAN")
        us_or_unknown = c["OCCR_COUNTRY"].str.upper().eq("US") | country_unknown
        teen_age = c["AGE_COD"].eq("YR") & c["AGE_NUM"].between(13, 19, inclusive="both")

        out = c[us_or_unknown & teen_age].copy()
        out["country_unknown"] = country_unknown[us_or_unknown & teen_age].astype(int)
        out = out.drop(columns=["AGE_NUM"])

        if not out.empty:
            out.to_csv(demo_teens, mode="a", index=False, header=not wrote_header)
            wrote_header = True
            teen_ids.update(out["PRIMARYID"].astype(str).str.strip().tolist())
        if idx % 20 == 0:
            log(f"DEMO teens filter chunks processed: {idx}")

    with teen_ids_txt.open("w", encoding="utf-8") as handle:
        for primary_id in sorted(teen_ids):
            handle.write(f"{primary_id}\n")

    return teen_ids


def filter_table_by_ids(input_csv: Path, out_csv: Path, id_set: set[str], chunksize: int) -> None:
    if out_csv.exists():
        out_csv.unlink()

    wrote_header = False
    for idx, chunk in enumerate(pd.read_csv(input_csv, dtype=str, chunksize=chunksize, low_memory=False), start=1):
        chunk = normalize_columns(chunk)
        ids = chunk["PRIMARYID"].astype(str).str.strip()
        filtered = chunk[ids.isin(id_set)]
        if not filtered.empty:
            filtered.to_csv(out_csv, mode="a", index=False, header=not wrote_header)
            wrote_header = True
        if idx % 40 == 0:
            log(f"Filtering {input_csv.name} chunks processed: {idx}")


def identify_dph_ids(drug_teens_csv: Path, chunksize: int) -> set[str]:
    log("Phase 7: Identifying diphenhydramine suspected/suspected-secondary cases")
    confirmed_ids: set[str] = set()

    for chunk in pd.read_csv(drug_teens_csv, dtype=str, chunksize=chunksize, low_memory=False):
        c = normalize_columns(chunk)
        drugname = c.get("DRUGNAME", "").astype(str).str.upper().str.strip()
        prod_ai = c.get("PROD_AI", "").astype(str).str.upper().str.strip()
        role = c.get("ROLE_COD", "").astype(str).str.upper().str.strip()

        pass1 = drugname.isin(DPH_TERMS) | drugname.str.contains("DIPHENHYDRAMIN", na=False)
        pass2 = prod_ai.str.contains("DIPHENHYDRAMIN", na=False)
        role_ok = role.isin({"PS", "SS"})

        matched = c[(pass1 | pass2) & role_ok]
        if not matched.empty:
            confirmed_ids.update(matched["PRIMARYID"].astype(str).str.strip().tolist())

    return confirmed_ids


def normalize_text(value: str) -> str:
    text = str(value).upper().strip()
    text = re.sub(r"[^A-Z0-9]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def build_rxnorm_map(rxnorm_rrf: Path, chunksize: int) -> dict[str, str]:
    log("Phase 8: Building RxNorm IN dictionary")
    rx_map: dict[str, str] = {}
    names = [f"C{i}" for i in range(19)]

    for chunk in pd.read_csv(
        rxnorm_rrf,
        sep="|",
        header=None,
        names=names,
        dtype=str,
        usecols=[0, 11, 12, 14],
        chunksize=chunksize,
        low_memory=False,
    ):
        c = chunk.rename(columns={"C11": "SAB", "C12": "TTY", "C14": "STR"})
        c = c[(c["SAB"] == "RXNORM") & (c["TTY"] == "IN")]
        for row in c.itertuples(index=False):
            term = normalize_text(row.STR)
            if term and term not in rx_map:
                rx_map[term] = str(row.STR).strip().lower()

    return rx_map


def map_generic_name(drugname: str, prod_ai: str, rx_map: dict[str, str]) -> tuple[str | None, int]:
    dn = normalize_text(drugname)
    pa = normalize_text(prod_ai)

    if dn in rx_map:
        return rx_map[dn], 1
    if pa in rx_map:
        return rx_map[pa], 1

    if pa:
        for token in re.split(r"[;,/+]| AND ", pa):
            t = normalize_text(token)
            if t in rx_map:
                return rx_map[t], 1

    return None, 0


def normalize_drug_table(drug_dph_csv: Path, out_csv: Path, rx_map: dict[str, str], chunksize: int) -> None:
    if out_csv.exists():
        out_csv.unlink()

    wrote_header = False
    for chunk in pd.read_csv(drug_dph_csv, dtype=str, chunksize=chunksize, low_memory=False):
        c = normalize_columns(chunk)
        drugnames = c.get("DRUGNAME", "").astype(str)
        prodais = c.get("PROD_AI", "").astype(str)

        mapped = [map_generic_name(dn, pa, rx_map) for dn, pa in zip(drugnames, prodais)]
        c["generic_name"] = [m[0] for m in mapped]
        c["rxnorm_matched"] = [m[1] for m in mapped]
        c.to_csv(out_csv, mode="a", index=False, header=not wrote_header)
        wrote_header = True


def score_acb_per_case(drug_norm_csv: Path, acb_csv: Path, out_csv: Path) -> pd.DataFrame:
    log("Phase 9: Calculating ACB per case")
    d = normalize_columns(pd.read_csv(drug_norm_csv, dtype=str, low_memory=False))
    a = pd.read_csv(acb_csv, dtype=str, low_memory=False)

    a["generic_name"] = a["generic_name"].astype(str).str.strip().str.lower()
    a["acb_score"] = pd.to_numeric(a["acb_score"], errors="coerce").fillna(0)
    acb_map = dict(zip(a["generic_name"], a["acb_score"]))

    d["generic_name"] = d["generic_name"].astype(str).str.strip().str.lower()
    d["acb_score"] = d["generic_name"].map(acb_map).fillna(0)

    d["is_dph"] = d["generic_name"].eq("diphenhydramine")

    d["drug_key"] = d["generic_name"].where(d["generic_name"].ne(""), d["DRUGNAME"].astype(str).str.lower())
    codrugs = d[~d["is_dph"]].copy()
    n_codrugs = codrugs.groupby("PRIMARYID")["drug_key"].nunique().rename("n_codrugs")

    total_with = d.groupby("PRIMARYID")["acb_score"].sum().rename("total_acb_with_dph")
    total_without = codrugs.groupby("PRIMARYID")["acb_score"].sum().rename("total_acb_codrugs_only")

    out = pd.concat([n_codrugs, total_with, total_without], axis=1).fillna(0).reset_index()
    out["n_codrugs"] = out["n_codrugs"].astype(int)
    out.to_csv(out_csv, index=False)
    return out


def aggregate_reac(reac_csv: Path) -> pd.DataFrame:
    log("Phase 10: Aggregating cardiac outcomes from REAC")
    r = normalize_columns(pd.read_csv(reac_csv, dtype=str, low_memory=False))
    r["PT"] = r["PT"].astype(str).str.strip()

    g = r.groupby("PRIMARYID")
    out = pd.DataFrame({"PRIMARYID": sorted(r["PRIMARYID"].astype(str).unique())})
    tier1 = g["PT"].apply(lambda s: int(s.isin(TIER1).any())).rename("cardiac_tier1")
    tier2 = g["PT"].apply(lambda s: int(s.isin(TIER2).any())).rename("cardiac_tier2")

    out = out.merge(tier1.reset_index(), on="PRIMARYID", how="left")
    out = out.merge(tier2.reset_index(), on="PRIMARYID", how="left")
    out[["cardiac_tier1", "cardiac_tier2"]] = out[["cardiac_tier1", "cardiac_tier2"]].fillna(0).astype(int)
    out["cardiac_any"] = ((out["cardiac_tier1"] == 1) | (out["cardiac_tier2"] == 1)).astype(int)
    return out


def aggregate_outc(outc_csv: Path) -> pd.DataFrame:
    log("Phase 11: Aggregating severity from OUTC")
    o = normalize_columns(pd.read_csv(outc_csv, dtype=str, low_memory=False))
    code = o["OUTC_COD"].astype(str).str.strip().str.upper()
    mapping = {"DE": 4, "LT": 3, "HO": 2, "RI": 1, "DS": 1, "OT": 1}
    o["severity"] = code.map(mapping).fillna(0).astype(int)

    agg = o.groupby("PRIMARYID", as_index=False)["severity"].max().rename(columns={"severity": "max_severity"})
    return agg


def build_final_table(
    demo_teens_csv: Path,
    confirmed_ids: set[str],
    acb_case: pd.DataFrame,
    cardiac: pd.DataFrame,
    severity: pd.DataFrame,
    out_csv: Path,
) -> None:
    log("Phase 12: Building final analysis table")
    d = normalize_columns(pd.read_csv(demo_teens_csv, dtype=str, low_memory=False))
    d["PRIMARYID"] = d["PRIMARYID"].astype(str).str.strip()
    d = d[d["PRIMARYID"].isin(confirmed_ids)].copy()

    keep_cols = [col for col in ["PRIMARYID", "AGE", "SEX", "EVENT_DT", "OCCR_COUNTRY"] if col in d.columns]
    final_df = d[keep_cols].drop_duplicates(subset=["PRIMARYID"], keep="first")

    final_df = final_df.merge(acb_case, on="PRIMARYID", how="left")
    final_df = final_df.merge(cardiac, on="PRIMARYID", how="left")
    final_df = final_df.merge(severity, on="PRIMARYID", how="left")

    for col in ["n_codrugs", "total_acb_with_dph", "total_acb_codrugs_only", "cardiac_any", "cardiac_tier1", "cardiac_tier2", "max_severity"]:
        if col not in final_df.columns:
            final_df[col] = 0

    final_df[["n_codrugs", "cardiac_any", "cardiac_tier1", "cardiac_tier2", "max_severity"]] = (
        final_df[["n_codrugs", "cardiac_any", "cardiac_tier1", "cardiac_tier2", "max_severity"]].fillna(0).astype(int)
    )
    final_df[["total_acb_with_dph", "total_acb_codrugs_only"]] = (
        final_df[["total_acb_with_dph", "total_acb_codrugs_only"]].fillna(0)
    )

    final_df["YEAR"] = final_df.get("EVENT_DT", "").astype(str).str.slice(0, 4)
    final_df["MONTH"] = pd.to_numeric(final_df.get("EVENT_DT", "").astype(str).str.slice(4, 6), errors="coerce")
    final_df["YEAR_NUM"] = pd.to_numeric(final_df["YEAR"], errors="coerce")

    age_num = pd.to_numeric(final_df.get("AGE", ""), errors="coerce")
    final_df["age_group"] = pd.cut(age_num, bins=[12, 15, 17, 19], labels=["13-15", "16-17", "18-19"])

    final_df["pre_post_warning"] = (
        (final_df["YEAR_NUM"] > 2020)
        | ((final_df["YEAR_NUM"] == 2020) & (final_df["MONTH"].fillna(0) >= 9))
    ).astype(int)

    final_df = final_df.drop(columns=["MONTH", "YEAR_NUM"])
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(out_csv, index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run full FAERS teen diphenhydramine analysis pipeline (Phases 3-12).")
    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--extracted-root", default="data/faers_extracted")
    parser.add_argument("--rxnorm-rrf", default="RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF")
    parser.add_argument("--acb-csv", default="data/lookups/acb_lookup.csv")
    parser.add_argument("--chunksize", type=int, default=200000)
    parser.add_argument(
        "--rebuild-combined",
        action="store_true",
        help="Rebuild 02_combined DEMO/DRUG/REAC/OUTC even if files already exist",
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    extracted_root = (base_dir / args.extracted_root).resolve()
    rxnorm_rrf = (base_dir / args.rxnorm_rrf).resolve()
    acb_csv = (base_dir / args.acb_csv).resolve()

    if not extracted_root.exists():
        raise FileNotFoundError(f"Extracted FAERS directory not found: {extracted_root}")
    if not rxnorm_rrf.exists():
        raise FileNotFoundError(f"RxNorm RRF not found: {rxnorm_rrf}")
    if not acb_csv.exists():
        raise FileNotFoundError(f"ACB CSV not found: {acb_csv}")

    paths = ensure_dirs(base_dir)

    # copy ACB into requested raw folder location
    acb_out = paths["raw_acb"] / "acb_scores.csv"
    pd.read_csv(acb_csv, dtype=str).loc[:, ["generic_name", "acb_score"]].to_csv(acb_out, index=False)

    combined = paths["combined"]
    filtered = paths["filtered"]
    processed = paths["processed"]

    combined_files = {}
    for table in TABLES:
        table_files = collect_quarter_files(extracted_root, table)
        out_csv = combined / f"{table}_all.csv"
        if out_csv.exists() and not args.rebuild_combined:
            log(f"Reusing existing combined table: {out_csv.name}")
        else:
            combine_table(table_files, out_csv, args.chunksize)
        combined_files[table] = out_csv

    demo_deduped = filtered / "DEMO_deduped.csv"
    demo_teens = filtered / "DEMO_teens.csv"
    teen_ids_txt = filtered / "teen_dph_ids.txt"

    dedupe_demo(combined_files["DEMO"], demo_deduped, args.chunksize)
    teen_ids = filter_demo_teens(demo_deduped, demo_teens, teen_ids_txt, args.chunksize)

    drug_teens = filtered / "DRUG_teens.csv"
    reac_teens = filtered / "REAC_teens.csv"
    outc_teens = filtered / "OUTC_teens.csv"

    log("Phase 6: Filtering DRUG/REAC/OUTC to teen IDs")
    filter_table_by_ids(combined_files["DRUG"], drug_teens, teen_ids, args.chunksize)
    filter_table_by_ids(combined_files["REAC"], reac_teens, teen_ids, args.chunksize)
    filter_table_by_ids(combined_files["OUTC"], outc_teens, teen_ids, args.chunksize)

    confirmed_ids = identify_dph_ids(drug_teens, args.chunksize)
    with teen_ids_txt.open("w", encoding="utf-8") as handle:
        for pid in sorted(confirmed_ids):
            handle.write(f"{pid}\n")

    drug_dph = filtered / "DRUG_dph_confirmed.csv"
    reac_dph = filtered / "REAC_dph_confirmed.csv"
    outc_dph = filtered / "OUTC_dph_confirmed.csv"

    filter_table_by_ids(drug_teens, drug_dph, confirmed_ids, args.chunksize)
    filter_table_by_ids(reac_teens, reac_dph, confirmed_ids, args.chunksize)
    filter_table_by_ids(outc_teens, outc_dph, confirmed_ids, args.chunksize)

    rx_map = build_rxnorm_map(rxnorm_rrf, args.chunksize)
    drug_norm = processed / "DRUG_normalized.csv"
    normalize_drug_table(drug_dph, drug_norm, rx_map, args.chunksize)

    acb_case_csv = processed / "ACB_per_case.csv"
    acb_case = score_acb_per_case(drug_norm, acb_out, acb_case_csv)

    cardiac = aggregate_reac(reac_dph)
    severity = aggregate_outc(outc_dph)

    final_csv = paths["final"] / "analysis_table.csv"
    build_final_table(demo_teens, confirmed_ids, acb_case, cardiac, severity, final_csv)

    log("Pipeline complete.")
    log(f"Final analysis table: {final_csv}")
    log(f"Confirmed teen DPH cohort size (unique PRIMARYID): {len(confirmed_ids)}")


if __name__ == "__main__":
    main()
