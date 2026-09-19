"""Fold the first min(1000, transcript length) nucleotides, at 37 C, for every gene."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from axon_rna.real_data import ROOT, parse_rnafold, sequence_features
import pandas as pd


def fold(row, window, version):
    result = dict(gene=row.gene, transcript_id=row.transcript_id, sequence_length=row.sequence_length,
                  structure_status="mapping_failed", window_start=1, window_end=0,
                  window_method=f"5prime_first_min_{window}_nt", rnafold_version=version,
                  temperature_celsius=37, structure_error="")
    if row.mapping_status != "mapped":
        return result
    seq = row.sequence[:window]
    result.update(window_end=len(seq), folded_length=len(seq), was_truncated=len(seq)<len(row.sequence),
                  gc_content=sequence_features(row.sequence)["gc_content"],
                  window_gc_content=sequence_features(seq)["gc_content"])
    cache = ROOT / "data/real/cache/folds"
    cache.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256((seq + version + "--noPS --temp=37").encode()).hexdigest()
    path = cache / f"{key}.txt"
    try:
        if path.exists():
            output = path.read_text()
        else:
            process = subprocess.run(["RNAfold", "--noPS", "--temp=37"], input=seq+"\n", text=True,
                                     capture_output=True, check=True, timeout=120)
            output = process.stdout
            parse_rnafold(output, seq)
            path.write_text(output)
        result.update(parse_rnafold(output, seq), structure_status="folded")
    except Exception as error:
        result.update(structure_status="failed", structure_error=str(error))
    return result


def plot_structures(sequences, structures):
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache/matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Arc
    dest = ROOT / "results/rna_structures/real"
    dest.mkdir(parents=True, exist_ok=True)
    for gene in ["VGF", "RPH3A", "NEUROD6", "FKBP5", "NRN1", "BDNF"]:
        subset = structures[(structures.gene == gene) & (structures.structure_status == "folded")]
        if subset.empty:
            continue
        row = subset.iloc[0]
        seq = sequences.set_index("gene").loc[gene, "sequence"][:int(row.folded_length)]
        fig, ax = plt.subplots(figsize=(11, 4), layout="constrained")
        stack = []
        max_span = 0
        for i, char in enumerate(row.dot_bracket, start=1):
            if char == "(":
                stack.append(i)
            elif char == ")":
                j = stack.pop()
                span = i-j
                max_span = max(max_span, span)
                ax.add_patch(Arc(((i+j)/2, 0), span, span, theta1=0, theta2=180,
                                 linewidth=.6, color="#287c8e", alpha=.65))
        ax.set(xlim=(0,len(seq)+1), ylim=(0,max_span/2+15), xlabel="Representative transcript position (nt)",
               title=f"{gene} • {row.transcript_id}\n5′ region: 1–{len(seq)} nt | RNAfold 37°C | MFE {row.mfe:.2f} kcal/mol")
        ax.set_yticks([])
        ax.spines[["top", "left", "right"]].set_visible(False)
        fig.savefig(dest/f"{gene}_structure.png", dpi=220)
        plt.close(fig)
        subprocess.run(["RNAplot", "--output-format=svg"],
                       input=f">{gene}\n{seq}\n{row.dot_bracket}\n", cwd=dest, text=True,
                       capture_output=True, check=True, timeout=60)
    structures[structures.gene.isin(["VGF", "RPH3A", "NEUROD6", "FKBP5", "NRN1", "BDNF"])].to_csv(dest/"structure_plot_manifest.csv", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--window", type=int, default=1000)
    args = ap.parse_args()
    if not 10 <= args.window <= 2000:
        ap.error("Use a window between 10 and 2000 nt to bound runtime")
    if not shutil.which("RNAfold"):
        raise SystemExit("ViennaRNA/RNAfold is required; no proxy fallback for real results")
    version = subprocess.check_output(["RNAfold", "--version"], text=True).strip()
    path = ROOT / "data/real/sequences/real_model_sequences.csv"
    if not path.exists():
        path = ROOT / "data/real/sequences/real_candidate_sequences.csv"
    sequences = pd.read_csv(path)
    with ThreadPoolExecutor(max_workers=4) as pool:
        structures = pd.DataFrame(pool.map(lambda r: fold(r,args.window,version), sequences.itertuples()))
    dest = ROOT / "data/real/structures"
    dest.mkdir(parents=True, exist_ok=True)
    structures.to_csv(dest/"real_structure_features.csv", index=False)
    plot_structures(sequences, structures)
    print(structures.structure_status.value_counts().to_json())


if __name__ == "__main__":
    main()
