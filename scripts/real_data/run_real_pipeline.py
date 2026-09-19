"""Rebuild new real-data outputs from the preserved consensus and cached resources."""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline",action="store_true",help="Require cached Ensembl/MyGene responses")
    args=parser.parse_args()
    stages=[("fetch_real_sequences.py",["--include-reference"]),
            ("fold_real_sequences.py",[]), ("build_real_dataset.py",[]),
            ("train_real_models.py",[]), ("plot_real_results.py",[]),
            ("extract_rna_embeddings.py",["--check"]), ("validate_real_project.py",[])]
    for script,options in stages:
        if args.offline and script in ["fetch_real_sequences.py","build_real_dataset.py"]:
            options.append("--offline")
        print(f"Running {script}",flush=True)
        subprocess.run([sys.executable,str(ROOT/"scripts/real_data"/script),*options],cwd=ROOT,check=True)


if __name__=="__main__":
    main()
