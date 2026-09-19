from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        running = min(running, ranked[i] * n / rank)
        adjusted[order[i]] = running
    return np.clip(adjusted, 0, 1)


def pick_id_column(df: pd.DataFrame, preferred: str | None):
    if preferred:
        if preferred not in df.columns:
            raise ValueError(f"--id-column {preferred!r} not found")
        return preferred
    candidates = [
        "transcript_id", "gene_symbol", "symbol", "Symbol", "gene", "Gene", "Geneid", "geneID", "EntrezGeneID", "id", "ID"
    ]
    for c in candidates:
        if c in df.columns:
            return c
    return str(df.columns[0])


def select_columns(columns, pattern):
    regex = re.compile(pattern, re.I)
    return [str(c) for c in columns if regex.search(str(c))]


def numeric_matrix(df, cols):
    out = df[cols].apply(pd.to_numeric, errors="coerce")
    return out


def main():
    ap = argparse.ArgumentParser(description="Transparent first-pass DE screen for public AD bulk RNA-seq matrices.")
    ap.add_argument("--input", required=True)
    ap.add_argument("--study-id", required=True)
    ap.add_argument("--ad-pattern", required=True, help="Regex selecting AD sample columns")
    ap.add_argument("--control-pattern", required=True, help="Regex selecting control sample columns")
    ap.add_argument("--id-column")
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-total", type=float, default=10.0)
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep=None, engine="python", compression="infer")
    id_col = pick_id_column(df, args.id_column)
    ad_cols = select_columns(df.columns, args.ad_pattern)
    ctl_cols = select_columns(df.columns, args.control_pattern)
    if not ad_cols or not ctl_cols:
        raise SystemExit(
            "Could not identify both AD and control columns. Run inspect_matrix.py first and adjust --ad-pattern/--control-pattern.\n"
            f"AD matched: {ad_cols}\nControl matched: {ctl_cols}"
        )

    ad = numeric_matrix(df, ad_cols)
    ctl = numeric_matrix(df, ctl_cols)
    values = pd.concat([ad, ctl], axis=1)
    keep = values.fillna(0).sum(axis=1) >= args.min_total
    df, ad, ctl = df.loc[keep].copy(), ad.loc[keep], ctl.loc[keep]

    # Library-size normalize to counts per million, then log2(CPM + 1).
    all_counts = pd.concat([ad, ctl], axis=1).clip(lower=0).fillna(0)
    lib = all_counts.sum(axis=0).replace(0, np.nan)
    logcpm = np.log2(all_counts.div(lib, axis=1) * 1e6 + 1.0)
    ad_log = logcpm[ad_cols]
    ctl_log = logcpm[ctl_cols]

    log2fc = ad_log.mean(axis=1) - ctl_log.mean(axis=1)
    stat, pvalue = ttest_ind(ad_log.to_numpy(), ctl_log.to_numpy(), axis=1, equal_var=False, nan_policy="omit")
    pvalue = np.nan_to_num(pvalue, nan=1.0, posinf=1.0, neginf=1.0)
    padj = bh_adjust(pvalue)

    ids = df[id_col].astype(str)
    # If IDs repeat, keep the strongest row by adjusted p-value.
    out = pd.DataFrame({
        "study_id": args.study_id,
        "transcript_id": ids,
        "log2fc": log2fc.to_numpy(),
        "pvalue": pvalue,
        "padj": padj,
        "n_ad": len(ad_cols),
        "n_control": len(ctl_cols),
    })
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2fc", "pvalue"])
    out = out.sort_values("padj").drop_duplicates("transcript_id", keep="first")
    dest = Path(args.output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest, index=False)

    print(f"study={args.study_id}")
    print(f"id_column={id_col}")
    print(f"AD columns ({len(ad_cols)}): {ad_cols}")
    print(f"Control columns ({len(ctl_cols)}): {ctl_cols}")
    print(f"features tested={len(out)}")
    print(f"FDR<0.05={(out.padj < 0.05).sum()}")
    print(f"saved={dest}")
    print("\nNOTE: this is a transparent screening analysis (logCPM + Welch test), not a replacement for DESeq2/edgeR in publication-grade inference.")


if __name__ == "__main__":
    main()
