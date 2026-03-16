import argparse
import csv
from pathlib import Path

import pandas as pd


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(col).strip().upper() for col in out.columns]
    return out


def to_numeric_id(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(-1).astype("int64")


def collect_quarter_files(extracted_root: Path, table_prefix: str) -> list[Path]:
    candidates: list[Path] = []
    for quarter_dir in sorted(d for d in extracted_root.iterdir() if d.is_dir()):
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


def build_max_primary_by_case(demo_files: list[Path], chunksize: int) -> dict[str, int]:
    max_primary_by_case: dict[str, int] = {}
    for file_path in demo_files:
        for chunk in read_faers_chunks(file_path, chunksize):
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
    return max_primary_by_case


def filtered_demo_rows(file_path: Path, chunksize: int, max_primary_by_case: dict[str, int]):
    for chunk in read_faers_chunks(file_path, chunksize):
        chunk["CASEID"] = chunk["CASEID"].astype(str).str.strip()
        chunk["PRIMARYID"] = chunk["PRIMARYID"].astype(str).str.strip()
        chunk["PRIMARYID_NUM"] = to_numeric_id(chunk["PRIMARYID"])
        keep_max = chunk["CASEID"].apply(lambda caseid: max_primary_by_case.get(caseid, -1)).astype("int64")
        deduped = chunk[chunk["PRIMARYID_NUM"] == keep_max].drop(columns=["PRIMARYID_NUM"])

        deduped["OCCR_COUNTRY"] = deduped.get("OCCR_COUNTRY", "").astype(str).str.strip()
        deduped["AGE_COD"] = deduped.get("AGE_COD", "").astype(str).str.strip().str.upper()
        deduped["AGE_NUM"] = pd.to_numeric(deduped.get("AGE", ""), errors="coerce")

        country_unknown = deduped["OCCR_COUNTRY"].eq("") | deduped["OCCR_COUNTRY"].str.upper().eq("NAN")
        us_or_unknown = deduped["OCCR_COUNTRY"].str.upper().eq("US") | country_unknown
        teen_age = deduped["AGE_COD"].eq("YR") & deduped["AGE_NUM"].between(13, 19, inclusive="both")

        out = deduped[us_or_unknown & teen_age].copy()
        if out.empty:
            continue
        out["country_unknown"] = country_unknown[us_or_unknown & teen_age].astype(int)
        out = out.drop(columns=["AGE_NUM"])
        yield out


def normalize_mixed_demo_csv(input_csv: Path, output_csv: Path, schema_columns: list[str]) -> tuple[int, int, int]:
    if output_csv.exists():
        output_csv.unlink()

    normalized_rows = 0
    passthrough_rows = 0
    remapped_rows = 0

    with input_csv.open("r", encoding="utf-8", errors="ignore", newline="") as source:
        reader = csv.reader(source)
        source_header = next(reader)
        expected_fields = len(schema_columns)
        source_fields = len(source_header)

        with output_csv.open("w", encoding="utf-8", newline="") as target:
            writer = csv.writer(target)
            writer.writerow(schema_columns)

            for row in reader:
                if not row:
                    continue
                if len(row) == expected_fields:
                    writer.writerow(row)
                    passthrough_rows += 1
                    normalized_rows += 1
                    continue

                if len(row) == 26 and expected_fields == 23:
                    mapped = [""] * expected_fields
                    mapped[0:9] = row[0:9]
                    mapped[9] = row[10]
                    mapped[10] = row[11]
                    mapped[11] = row[13]
                    mapped[12] = row[14]
                    mapped[13] = row[16]
                    mapped[14] = row[17]
                    mapped[15] = row[18]
                    mapped[16] = row[19]
                    mapped[17] = row[20]
                    mapped[18] = row[21]
                    mapped[19] = row[22]
                    mapped[20] = row[23]
                    mapped[21] = row[24]
                    mapped[22] = row[25]
                    writer.writerow(mapped)
                    remapped_rows += 1
                    normalized_rows += 1
                    continue

                raise ValueError(
                    f"Unsupported DEMO row width {len(row)} in {input_csv}. "
                    f"Source header has {source_fields} columns and target schema has {expected_fields}."
                )

    return normalized_rows, passthrough_rows, remapped_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit and repair the teen DEMO cohort against raw FAERS DEMO source.")
    parser.add_argument("--base-dir", default=".")
    parser.add_argument("--recovered-csv", default="03_filtered/teen_demo_records.csv")
    parser.add_argument("--extracted-root", default="data/faers_extracted")
    parser.add_argument("--chunksize", type=int, default=100000)
    parser.add_argument("--write-missing-csv", default="")
    parser.add_argument("--write-repaired-csv", default="")
    parser.add_argument("--schema-from-csv", default="")
    parser.add_argument("--normalize-input-csv", default="")
    parser.add_argument("--normalize-output-csv", default="")
    args = parser.parse_args()

    base_dir = Path(args.base_dir).resolve()
    recovered_csv = (base_dir / args.recovered_csv).resolve()
    extracted_root = (base_dir / args.extracted_root).resolve()
    schema_from_csv = (base_dir / args.schema_from_csv).resolve() if args.schema_from_csv else recovered_csv

    schema_columns = list(pd.read_csv(schema_from_csv, dtype=str, nrows=0).columns)

    if args.normalize_input_csv or args.normalize_output_csv:
        if not args.normalize_input_csv or not args.normalize_output_csv:
            raise ValueError("--normalize-input-csv and --normalize-output-csv must be provided together")
        normalize_input_csv = (base_dir / args.normalize_input_csv).resolve()
        normalize_output_csv = (base_dir / args.normalize_output_csv).resolve()
        normalized_rows, passthrough_rows, remapped_rows = normalize_mixed_demo_csv(
            normalize_input_csv,
            normalize_output_csv,
            schema_columns,
        )
        print(f"normalized_csv={normalize_output_csv}")
        print(f"normalized_rows={normalized_rows}")
        print(f"normalized_passthrough_rows={passthrough_rows}")
        print(f"normalized_remapped_rows={remapped_rows}")
        return

    recovered_ids = set(
        pd.read_csv(recovered_csv, dtype=str, usecols=["PRIMARYID"], low_memory=False)["PRIMARYID"]
        .astype(str)
        .str.strip()
    )
    demo_files = collect_quarter_files(extracted_root, "DEMO")
    max_primary_by_case = build_max_primary_by_case(demo_files, args.chunksize)

    expected_total = 0
    recovered_total = 0
    missing_total = 0
    missing_quarters: list[tuple[str, int, int, int]] = []
    missing_rows_written = 0
    missing_writer_header = False
    write_missing_csv = Path(args.write_missing_csv).resolve() if args.write_missing_csv else None
    write_repaired_csv = Path(args.write_repaired_csv).resolve() if args.write_repaired_csv else None

    if write_missing_csv and write_missing_csv.exists():
        write_missing_csv.unlink()

    for file_path in demo_files:
        quarter = file_path.parent.name
        quarter_expected = 0
        quarter_recovered = 0
        quarter_missing = 0

        for out in filtered_demo_rows(file_path, args.chunksize, max_primary_by_case):
            ids = out["PRIMARYID"].astype(str).str.strip()
            in_recovered = ids.isin(recovered_ids)
            missing = out.loc[~in_recovered].copy()

            if not missing.empty:
                for column in schema_columns:
                    if column not in missing.columns:
                        missing[column] = ""
                missing = missing.loc[:, schema_columns]

            quarter_expected += len(out)
            quarter_recovered += int(in_recovered.sum())
            quarter_missing += len(missing)

            if write_missing_csv and not missing.empty:
                missing.to_csv(write_missing_csv, mode="a", index=False, header=not missing_writer_header)
                missing_writer_header = True
                missing_rows_written += len(missing)

        expected_total += quarter_expected
        recovered_total += quarter_recovered
        missing_total += quarter_missing
        missing_quarters.append((quarter, quarter_expected, quarter_recovered, quarter_missing))

    print(f"recovered_rows={len(recovered_ids)}")
    print(f"expected_rows={expected_total}")
    print(f"recovered_matched_rows={recovered_total}")
    print(f"missing_rows={missing_total}")
    if write_missing_csv:
        print(f"missing_csv={write_missing_csv}")
        print(f"missing_rows_written={missing_rows_written}")

    print("quarter_summary")
    for quarter, quarter_expected, quarter_recovered, quarter_missing in missing_quarters:
        if quarter_missing:
            print(
                f"{quarter}\texpected={quarter_expected}\t"
                f"recovered={quarter_recovered}\tmissing={quarter_missing}"
            )

    if write_repaired_csv:
        if not write_missing_csv or not write_missing_csv.exists():
            raise FileNotFoundError("--write-repaired-csv requires a missing CSV produced by --write-missing-csv")

        if write_repaired_csv.exists():
            write_repaired_csv.unlink()

        missing_caseids: set[str] = set()
        missing_raw_lines: list[str] = []
        with write_missing_csv.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
            missing_header = handle.readline()
            for line in handle:
                if not line.strip():
                    continue
                parts = line.split(",", 2)
                if len(parts) < 2:
                    continue
                missing_caseids.add(parts[1].strip())
                missing_raw_lines.append(line)

        removed_rows = 0
        kept_rows = 0
        replaced_caseids: set[str] = set()

        with recovered_csv.open("r", encoding="utf-8", errors="ignore", newline="") as source:
            header = source.readline()
            with write_repaired_csv.open("w", encoding="utf-8", newline="") as repaired:
                repaired.write(header)
                for line in source:
                    if not line.strip():
                        continue
                    parts = line.split(",", 2)
                    if len(parts) < 2:
                        repaired.write(line)
                        kept_rows += 1
                        continue
                    caseid = parts[1].strip()
                    if caseid in missing_caseids:
                        removed_rows += 1
                        replaced_caseids.add(caseid)
                        continue
                    repaired.write(line)
                    kept_rows += 1

                for line in missing_raw_lines:
                    repaired.write(line)

        print(f"repaired_csv={write_repaired_csv}")
        print(f"repair_removed_rows={removed_rows}")
        print(f"repair_kept_rows={kept_rows}")
        print(f"repair_inserted_rows={len(missing_raw_lines)}")
        print(f"repair_replaced_caseids={len(replaced_caseids)}")
        print(f"repair_new_caseids={len(missing_caseids - replaced_caseids)}")


if __name__ == "__main__":
    main()