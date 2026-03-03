import argparse
import csv
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

BASE_URL = "https://fis.fda.gov/content/Exports/faers_ascii_{year}q{quarter}.zip"
TARGET_PREFIXES = ("DEMO", "DRUG", "REAC", "OUTC")


def quarter_iter(start_year: int, end_year: int, end_quarter: int):
    for year in range(start_year, end_year + 1):
        max_quarter = end_quarter if year == end_year else 4
        for quarter in range(1, max_quarter + 1):
            yield year, quarter


def download_file(url: str, output_path: Path, timeout: int, retries: int) -> bool:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                if response.status != 200:
                    raise urllib.error.HTTPError(url, response.status, "HTTP error", response.headers, None)

                with output_path.open("wb") as file_handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        file_handle.write(chunk)
            return True
        except Exception as exc:
            if attempt == retries:
                print(f"  x Download failed after {retries} attempt(s): {exc}")
                return False
            print(f"  ! Download attempt {attempt}/{retries} failed: {exc}; retrying...")
            time.sleep(2)

    return False


def extract_target_files(zip_path: Path, extract_dir: Path) -> int:
    extract_dir.mkdir(parents=True, exist_ok=True)
    extracted_count = 0

    with zipfile.ZipFile(zip_path, "r") as archive:
        members = archive.namelist()
        target_members = []

        for member in members:
            filename = Path(member).name.upper()
            if filename.startswith(TARGET_PREFIXES) and filename.endswith(".TXT"):
                target_members.append(member)

        for member in target_members:
            member_name = Path(member).name
            destination = extract_dir / member_name
            with archive.open(member) as source, destination.open("wb") as target:
                target.write(source.read())
            extracted_count += 1

    return extracted_count


def is_valid_zip(zip_path: Path) -> bool:
    return zip_path.exists() and zipfile.is_zipfile(zip_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download FAERS quarterly ASCII ZIPs and extract DEMO/DRUG/REAC/OUTC tables."
    )
    parser.add_argument("--start-year", type=int, default=2013)
    parser.add_argument("--end-year", type=int, default=2025)
    parser.add_argument("--end-quarter", type=int, default=4, choices=[1, 2, 3, 4])
    parser.add_argument("--raw-dir", default="data/faers_raw")
    parser.add_argument("--extract-dir", default="data/faers_extracted")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--sleep-seconds", type=float, default=1.0)
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions without downloading")
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    extract_root = Path(args.extract_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    extract_root.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    failures = []

    for year, quarter in quarter_iter(args.start_year, args.end_year, args.end_quarter):
        zip_name = f"faers_ascii_{year}q{quarter}.zip"
        zip_path = raw_dir / zip_name
        quarter_dir = extract_root / f"{year}q{quarter}"
        url = BASE_URL.format(year=year, quarter=quarter)

        print(f"\n[{year} Q{quarter}] {zip_name}")
        print(f"  URL: {url}")

        if args.dry_run:
            manifest_rows.append([f"{year}Q{quarter}", zip_name, "dry_run", 0, str(quarter_dir)])
            continue

        if zip_path.exists() and not is_valid_zip(zip_path):
            print("  ! Existing ZIP is invalid/corrupted; deleting and re-downloading")
            zip_path.unlink(missing_ok=True)

        if zip_path.exists():
            print("  - ZIP already exists; skipping download")
        else:
            print("  - Downloading ZIP...")
            success = download_file(url, zip_path, timeout=args.timeout, retries=args.retries)
            if not success:
                failures.append(zip_name)
                manifest_rows.append([f"{year}Q{quarter}", zip_name, "download_failed", 0, str(quarter_dir)])
                continue
            if not is_valid_zip(zip_path):
                print("  x Downloaded file is not a valid ZIP (likely unavailable quarter or HTML response)")
                zip_path.unlink(missing_ok=True)
                failures.append(zip_name)
                manifest_rows.append([f"{year}Q{quarter}", zip_name, "invalid_zip", 0, str(quarter_dir)])
                continue

        try:
            extracted = extract_target_files(zip_path, quarter_dir)
            print(f"  - Extracted {extracted} target file(s) to {quarter_dir}")
            status = "ok" if extracted > 0 else "no_target_files"
            manifest_rows.append([f"{year}Q{quarter}", zip_name, status, extracted, str(quarter_dir)])
        except Exception as exc:
            print(f"  x Extraction failed: {exc}")
            failures.append(zip_name)
            manifest_rows.append([f"{year}Q{quarter}", zip_name, "extract_failed", 0, str(quarter_dir)])
            continue

        time.sleep(args.sleep_seconds)

    manifest_path = extract_root / "faers_download_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as manifest_file:
        writer = csv.writer(manifest_file)
        writer.writerow(["quarter", "zip_file", "status", "extracted_target_files", "extract_dir"])
        writer.writerows(manifest_rows)

    print("\n========== COMPLETE ==========")
    print(f"Manifest: {manifest_path}")
    if failures:
        print(f"Failures ({len(failures)}):")
        for item in failures:
            print(f"- {item}")
    else:
        print("All requested quarters processed.")


if __name__ == "__main__":
    main()
