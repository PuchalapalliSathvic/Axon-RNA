from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--output", default="data/real/processed/studies.csv")
    args = ap.parse_args()
    frames = [pd.read_csv(f) for f in args.files]
    out = pd.concat(frames, ignore_index=True)
    required = {"study_id", "transcript_id", "log2fc", "pvalue", "padj"}
    missing = required - set(out.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")
    dest = Path(args.output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest, index=False)
    print(out.groupby("study_id").size().to_string())
    print(f"saved={dest}")


if __name__ == "__main__":
    main()
