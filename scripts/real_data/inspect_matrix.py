from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--rows", type=int, default=5)
    args = ap.parse_args()
    path = Path(args.file)
    df = pd.read_csv(path, sep=None, engine="python", compression="infer")
    print(f"shape={df.shape}")
    print("columns:")
    for i, c in enumerate(df.columns):
        print(f"  {i:>3}: {c}")
    print("\nhead:")
    print(df.head(args.rows).to_string(index=False))


if __name__ == "__main__":
    main()
