import argparse
import re
from pathlib import Path

import pandas as pd

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


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(col).strip().upper() for col in out.columns]
    return out


def to_numeric_id(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(-1).astype("int64")


def get_series_str(df: pd.DataFrame, col: str) -> pd.Series:
    if col in df.columns:
        return df[col].astype(str)
    return pd.Series("", index=df.index, dtype=str)


def collect_quarter_files(extracted_root: Path, table_prefix: str) -> list[Path]:
    candidates: list[Path] = []
    for quarter_dir in sorted([d for d in extracted_root.iterdir() if d.is_dir()]):
        for file_path in sorted(quarter_dir.iterdir()):
            name = file_path.name.upper()
            if name.startswith(table_prefix) and name.endswith(".TXT"):
                candidates.append(file_path)
    return candidates


def read_faers_chunks(file_path: Path, chunksize: int):
    for chunk in pd.read_csv(
        file_path,
        sep="$",
        dtype=str,
        encoding="latin-1",
        on_bad_lines="skip",
        engine="c",
        chunksize=chunksize,
        low_memory=False,
    ):
        chunk = normalize_columns(chunk)
        drop_cols = [c for c in chunk.columns if c.startswith("UNNAMED") or c.strip() == ""]
        if drop_cols:
            chunk = chunk.drop(columns=drop_cols, errors="ignore")
        yield chunk


def filter_by_ids_from_txt_files(
    files: list[Path],
    out_csv: Path,
    id_set: set[str],
    chunksize: int,
) -> None:
    if out_csv.exists():
        out_csv.unlink()
    wrote_header = False
    for file_idx, file_path in enumerate(files, start=1):
        log(f"Filtering {out_csv.name}: {file_idx}/{len(files)} -> {file_path.name}")
        for chunk_idx, chunk in enumerate(read_faers_chunks(file_path, chunksize), start=1):
            ids = chunk["PRIMARYID"].astype(str).str.strip()
            filtered = chunk[ids.isin(id_set)]
            if not filtered.empty:
                filtered.to_csv(out_csv, mode="a", index=False, header=not wrote_header)
                wrote_header = True
            if chunk_idx % 20 == 0:
                log(f"{out_csv.name} chunks in {file_path.name}: {chunk_idx}")


def filter_csv_by_ids(input_csv: Path, out_csv: Path, id_set: set[str], chunksize: int) -> None:
    if out_csv.exists():
        out_csv.unlink()

    wrote_header = False
    for chunk in pd.read_csv(
        input_csv,
        dtype=str,
        chunksize=chunksize,
        on_bad_lines="skip",
        engine="python",
    ):
        c = normalize_columns(chunk)
        ids = get_series_str(c, "PRIMARYID").str.strip()
        filtered = c[ids.isin(id_set)]
        if not filtered.empty:
            filtered.to_csv(out_csv, mode="a", index=False, header=not wrote_header)
            wrote_header = True


def identify_dph_ids(drug_teens_csv: Path, chunksize: int) -> set[str]:
    confirmed_ids: set[str] = set()
    for chunk in pd.read_csv(
        drug_teens_csv,
        dtype=str,
        chunksize=chunksize,
        on_bad_lines="skip",
        engine="python",
    ):
        c = normalize_columns(chunk)
        drugname = get_series_str(c, "DRUGNAME").str.upper().str.strip()
        prod_ai = get_series_str(c, "PROD_AI").str.upper().str.strip()
        role = get_series_str(c, "ROLE_COD").str.upper().str.strip()

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
    log("Building RxNorm IN dictionary")
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
    for chunk in pd.read_csv(
        drug_dph_csv,
        dtype=str,
        chunksize=chunksize,
        on_bad_lines="skip",
        engine="python",
    ):
        c = normalize_columns(chunk)
        drugnames = get_series_str(c, "DRUGNAME")
        prodais = get_series_str(c, "PROD_AI")

        mapped = [map_generic_name(dn, pa, rx_map) for dn, pa in zip(drugnames, prodais)]
        c["generic_name"] = [m[0] for m in mapped]
        c["rxnorm_matched"] = [m[1] for m in mapped]
        c.to_csv(out_csv, mode="a", index=False, header=not wrote_header)
        wrote_header = True


def score_acb_per_case(drug_norm_csv: Path, acb_csv: Path, out_csv: Path) -> pd.DataFrame:
    d = normalize_columns(pd.read_csv(drug_norm_csv, dtype=str, on_bad_lines="skip", engine="python"))
    a = normalize_columns(pd.read_csv(acb_csv, dtype=str, on_bad_lines="skip", engine="python"))

    a["GENERIC_NAME"] = a["GENERIC_NAME"].astype(str).str.strip().str.lower()
    a["ACB_SCORE"] = pd.to_numeric(a["ACB_SCORE"], errors="coerce").fillna(0)
    acb_map = dict(zip(a["GENERIC_NAME"], a["ACB_SCORE"]))

    d["GENERIC_NAME"] = d["GENERIC_NAME"].astype(str).str.strip().str.lower()
    d["acb_score"] = d["GENERIC_NAME"].map(acb_map).fillna(0)

    d["is_dph"] = d["GENERIC_NAME"].eq("diphenhydramine")
    d["drug_key"] = d["GENERIC_NAME"].where(d["GENERIC_NAME"].ne(""), d["DRUGNAME"].astype(str).str.lower())

    codrugs = d[~d["is_dph"]].copy()
    n_codrugs = codrugs.groupby("PRIMARYID")["drug_key"].nunique().rename("n_codrugs")
    total_with = d.groupby("PRIMARYID")["acb_score"].sum().rename("total_acb_with_dph")
    total_without = codrugs.groupby("PRIMARYID")["acb_score"].sum().rename("total_acb_codrugs_only")

    out = pd.concat([n_codrugs, total_with, total_without], axis=1).fillna(0).reset_index()
    out["n_codrugs"] = out["n_codrugs"].astype(int)
    out.to_csv(out_csv, index=False)
    return out


def aggregate_reac(reac_csv: Path) -> pd.DataFrame:
    r = normalize_columns(pd.read_csv(reac_csv, dtype=str, on_bad_lines="skip", engine="python"))
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
    o = normalize_columns(pd.read_csv(outc_csv, dtype=str, on_bad_lines="skip", engine="python"))
    code = o["OUTC_COD"].astype(str).str.strip().str.upper()
    mapping = {"DE": 4, "LT": 3, "HO": 2, "RI": 1, "DS": 1, "OT": 1}
    o["severity"] = code.map(mapping).fillna(0).astype(int)
    return o.groupby("PRIMARYID", as_index=False)["severity"].max().rename(columns={"severity": "max_severity"})


def build_final_table(
    demo_teens_csv: Path,
    confirmed_ids: set[str],
    acb_case: pd.DataFrame,
    cardiac: pd.DataFrame,
    severity: pd.DataFrame,
    out_csv: Path,
) -> None:
    d = normalize_columns(pd.read_csv(demo_teens_csv, dtype=str, on_bad_lines="skip", engine="python"))
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
    parser = argparse.ArgumentParser(description="Compact FAERS pipeline (no huge 02_combined intermediates).")
    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--extracted-root", default="data/faers_extracted")
    parser.add_argument("--rxnorm-rrf", default="RxNorm/RxNorm_full_prescribe_03022026/rrf/RXNCONSO.RRF")
    parser.add_argument("--acb-csv", default="data/lookups/acb_lookup.csv")
    parser.add_argument("--chunksize", type=int, default=100000)
    parser.add_argument(
        "--start-phase",
        type=int,
        default=4,
        choices=[4, 7, 8],
        help="Start at phase 4 (default), 7, or 8 for recovery runs.",
    )
    parser.add_argument(
        "--resume-demo-pass2-from",
        type=int,
        default=1,
        help="Resume DEMO pass2 from this 1-based quarter index (e.g., 41).",
    )
    args = parser.parse_args()

    base = Path(args.base_dir).resolve()
    extracted_root = (base / args.extracted_root).resolve()
    rxnorm_rrf = (base / args.rxnorm_rrf).resolve()
    acb_csv = (base / args.acb_csv).resolve()

    filtered_dir = base / "03_filtered"
    processed_dir = base / "04_processed"
    final_dir = base / "05_final"
    raw_acb_dir = base / "01_raw" / "acb"
    for p in [filtered_dir, processed_dir, final_dir, raw_acb_dir]:
        p.mkdir(parents=True, exist_ok=True)

    if not extracted_root.exists():
        raise FileNotFoundError(f"Extracted FAERS directory not found: {extracted_root}")

    demo_files = collect_quarter_files(extracted_root, "DEMO")
    drug_files = collect_quarter_files(extracted_root, "DRUG")
    reac_files = collect_quarter_files(extracted_root, "REAC")
    outc_files = collect_quarter_files(extracted_root, "OUTC")

    demo_teens = filtered_dir / "DEMO_teens.csv"
    teen_ids_txt = filtered_dir / "teen_dph_ids.txt"
    drug_teens = filtered_dir / "DRUG_teens.csv"
    reac_teens = filtered_dir / "REAC_teens.csv"
    outc_teens = filtered_dir / "OUTC_teens.csv"

    drug_dph = filtered_dir / "DRUG_dph_confirmed.csv"
    reac_dph = filtered_dir / "REAC_dph_confirmed.csv"
    outc_dph = filtered_dir / "OUTC_dph_confirmed.csv"

    if args.start_phase <= 6:
        log("Phase 4+5: DEMO dedupe and teen US/unknown filter")
        max_primary_by_case: dict[str, int] = {}
        for file_idx, file_path in enumerate(demo_files, start=1):
            log(f"DEMO pass1 {file_idx}/{len(demo_files)} -> {file_path.name}")
            for chunk in read_faers_chunks(file_path, args.chunksize):
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

        if args.resume_demo_pass2_from < 1 or args.resume_demo_pass2_from > len(demo_files):
            raise ValueError(f"--resume-demo-pass2-from must be between 1 and {len(demo_files)}")

        teen_ids: set[str] = set()
        if args.resume_demo_pass2_from == 1:
            if demo_teens.exists():
                demo_teens.unlink()
            wrote_header = False
        else:
            if not demo_teens.exists():
                raise FileNotFoundError(
                    f"Cannot resume DEMO pass2 at {args.resume_demo_pass2_from} because DEMO_teens.csv does not exist"
                )
            existing_ids = pd.read_csv(demo_teens, dtype=str, usecols=["PRIMARYID"], low_memory=False)
            teen_ids.update(existing_ids["PRIMARYID"].astype(str).str.strip().tolist())
            wrote_header = demo_teens.stat().st_size > 0
            log(
                f"Resuming DEMO pass2 from {args.resume_demo_pass2_from}/{len(demo_files)} with "
                f"{len(teen_ids)} existing teen IDs"
            )

        for file_idx, file_path in enumerate(demo_files, start=1):
            if file_idx < args.resume_demo_pass2_from:
                continue
            log(f"DEMO pass2 {file_idx}/{len(demo_files)} -> {file_path.name}")
            for chunk in read_faers_chunks(file_path, args.chunksize):
                chunk["CASEID"] = chunk["CASEID"].astype(str).str.strip()
                chunk["PRIMARYID_NUM"] = to_numeric_id(chunk["PRIMARYID"])
                keep_max = chunk["CASEID"].map(max_primary_by_case).fillna(-1).astype("int64")
                deduped = chunk[chunk["PRIMARYID_NUM"] == keep_max].drop(columns=["PRIMARYID_NUM"])

                deduped["OCCR_COUNTRY"] = deduped.get("OCCR_COUNTRY", "").astype(str).str.strip()
                deduped["AGE_COD"] = deduped.get("AGE_COD", "").astype(str).str.strip().str.upper()
                deduped["AGE_NUM"] = pd.to_numeric(deduped.get("AGE", ""), errors="coerce")

                country_unknown = deduped["OCCR_COUNTRY"].eq("") | deduped["OCCR_COUNTRY"].str.upper().eq("NAN")
                us_or_unknown = deduped["OCCR_COUNTRY"].str.upper().eq("US") | country_unknown
                teen_age = deduped["AGE_COD"].eq("YR") & deduped["AGE_NUM"].between(13, 19, inclusive="both")

                out = deduped[us_or_unknown & teen_age].copy()
                out["country_unknown"] = country_unknown[us_or_unknown & teen_age].astype(int)
                out = out.drop(columns=["AGE_NUM"])

                if not out.empty:
                    out.to_csv(demo_teens, mode="a", index=False, header=not wrote_header)
                    wrote_header = True
                    teen_ids.update(out["PRIMARYID"].astype(str).str.strip().tolist())

        with teen_ids_txt.open("w", encoding="utf-8") as handle:
            for pid in sorted(teen_ids):
                handle.write(f"{pid}\n")

        log("Phase 6: filter DRUG/REAC/OUTC to teen IDs")
        filter_by_ids_from_txt_files(drug_files, drug_teens, teen_ids, args.chunksize)
        filter_by_ids_from_txt_files(reac_files, reac_teens, teen_ids, args.chunksize)
        filter_by_ids_from_txt_files(outc_files, outc_teens, teen_ids, args.chunksize)
    else:
        log(f"Skipping phases 4-6 (start-phase={args.start_phase})")

    if args.start_phase <= 7:
        log("Phase 7: identify confirmed diphenhydramine cohort")
        confirmed_ids = identify_dph_ids(drug_teens, args.chunksize)
        with teen_ids_txt.open("w", encoding="utf-8") as handle:
            for pid in sorted(confirmed_ids):
                handle.write(f"{pid}\n")

        filter_csv_by_ids(drug_teens, drug_dph, confirmed_ids, args.chunksize)
        filter_csv_by_ids(reac_teens, reac_dph, confirmed_ids, args.chunksize)
        filter_csv_by_ids(outc_teens, outc_dph, confirmed_ids, args.chunksize)
    else:
        if not teen_ids_txt.exists():
            raise FileNotFoundError("teen_dph_ids.txt is required when --start-phase=8")
        confirmed_ids = {
            line.strip()
            for line in teen_ids_txt.read_text(encoding="utf-8", errors="ignore").splitlines()
            if line.strip()
        }
        if not drug_dph.exists() or not reac_dph.exists() or not outc_dph.exists():
            raise FileNotFoundError("DRUG/REAC/OUTC DPH-confirmed files are required when --start-phase=8")
        log("Skipping phase 7 (start-phase=8); using existing confirmed cohort files")

    log("Phase 8-12: normalize drugs, score ACB, aggregate outcomes, build final table")
    acb_out = raw_acb_dir / "acb_scores.csv"
    normalize_columns(pd.read_csv(acb_csv, dtype=str)).loc[:, ["GENERIC_NAME", "ACB_SCORE"]].rename(
        columns={"GENERIC_NAME": "generic_name", "ACB_SCORE": "acb_score"}
    ).to_csv(acb_out, index=False)

    rx_map = build_rxnorm_map(rxnorm_rrf, args.chunksize)
    drug_norm = processed_dir / "DRUG_normalized.csv"
    normalize_drug_table(drug_dph, drug_norm, rx_map, args.chunksize)

    acb_case_csv = processed_dir / "ACB_per_case.csv"
    acb_case = score_acb_per_case(drug_norm, acb_out, acb_case_csv)

    cardiac = aggregate_reac(reac_dph)
    severity = aggregate_outc(outc_dph)

    final_csv = final_dir / "analysis_table.csv"
    build_final_table(demo_teens, confirmed_ids, acb_case, cardiac, severity, final_csv)

    log(f"Pipeline complete. Final table: {final_csv}")
    log(f"Confirmed cohort (unique PRIMARYID): {len(confirmed_ids)}")


if __name__ == "__main__":
    main()
