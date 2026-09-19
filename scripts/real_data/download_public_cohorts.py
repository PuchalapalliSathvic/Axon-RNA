from __future__ import annotations

import argparse
import csv
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REGISTRY = ROOT / "data" / "real" / "dataset_registry.csv"
RAW = ROOT / "data" / "real" / "raw"
RAW.mkdir(parents=True, exist_ok=True)


def geo_suppl_url(accession: str, filename: str) -> str:
    # GEO groups accessions by the first digits + 'nnn', e.g. GSE163877 -> GSE163nnn.
    number = accession.replace("GSE", "")
    family = f"GSE{number[:-3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{family}/{accession}/suppl/{filename}"


def load_registry():
    with REGISTRY.open(newline="") as f:
        return list(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser(description="Download AXON-RNA public GEO cohorts.")
    ap.add_argument("--accessions", nargs="*", default=["GSE53697", "GSE159699", "GSE163877"])
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    rows = {r["accession"]: r for r in load_registry()}
    for accession in args.accessions:
        row = rows.get(accession)
        if not row:
            raise SystemExit(f"Unknown accession: {accession}")
        filename = row["supplementary_file"]
        if filename == "per_sample_tsv":
            print(f"SKIP {accession}: per-sample download is handled in a later integration step.")
            continue
        dest = RAW / filename
        if dest.exists() and not args.force:
            print(f"EXISTS {dest}")
            continue
        url = geo_suppl_url(accession, filename)
        print(f"Downloading {accession}\n  {url}\n  -> {dest}")
        try:
            urllib.request.urlretrieve(url, dest)
        except Exception as exc:
            print(f"FAILED {accession}: {exc}")
            print("You can paste the printed URL into a browser and save the file into data/real/raw/.")


if __name__ == "__main__":
    main()
