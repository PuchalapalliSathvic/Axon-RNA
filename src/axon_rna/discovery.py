from __future__ import annotations
import numpy as np
import pandas as pd

REQUIRED = {"study_id", "transcript_id", "log2fc", "pvalue"}


def build_consensus(df: pd.DataFrame, alpha: float = 0.05, min_studies: int = 2) -> pd.DataFrame:
    """Create transcript-level consensus evidence across AD studies.

    Uses median effect size, sign consistency, Fisher combined p-value,
    and a simple reproducibility score. This is intentionally transparent.
    """
    missing = REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    x = df.copy()
    x["pvalue"] = x["pvalue"].clip(1e-300, 1.0)
    x["sign"] = np.sign(x["log2fc"])

    rows = []
    for tid, g in x.groupby("transcript_id"):
        k = len(g)
        med = float(g["log2fc"].median())
        direction = 1.0 if med >= 0 else -1.0
        consistency = float((g["sign"] == direction).mean())
        fisher_stat = float(-2 * np.log(g["pvalue"]).sum())
        try:
            from scipy.stats import chi2
            combined_p = float(chi2.sf(fisher_stat, 2 * k))
        except Exception:
            combined_p = float(min(1.0, g["pvalue"].min() * k))
        sig_fraction = float((g.get("padj", g["pvalue"]) < alpha).mean())
        evidence = abs(med) * consistency * (0.5 + sig_fraction)
        reproducible = int(k >= min_studies and consistency >= 0.75 and combined_p < alpha)
        rows.append({
            "transcript_id": tid,
            "n_studies": k,
            "median_log2fc": med,
            "direction_consistency": consistency,
            "combined_p": combined_p,
            "significant_fraction": sig_fraction,
            "evidence_score": evidence,
            "ad_dysregulated": reproducible,
        })
    return pd.DataFrame(rows).sort_values(["ad_dysregulated", "evidence_score"], ascending=[False, False]).reset_index(drop=True)
