"""Optional RNA-FM extraction; never substitutes random vectors for embeddings."""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from axon_rna.real_data import ROOT


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check",action="store_true",help="Report optional dependencies without downloads")
    args=ap.parse_args()
    available={module:importlib.util.find_spec(module) is not None for module in ["torch","fm"]}
    if args.check or not all(available.values()):
        status={"completed":False,"dependencies":available,"reason":"RNA-FM extraction not run; requires rna-fm package and pretrained weights",
                "sequence_strategy":"first min(1000, transcript length) nt; mean pool nucleotide representations", "device":"cpu"}
        (ROOT/"results/foundation_model_status.json").write_text(json.dumps(status,indent=2))
        print(json.dumps(status,indent=2))
        if not args.check:
            raise SystemExit("Optional dependencies missing. Classical results are unaffected.")
        return
    import numpy as np
    import pandas as pd
    from axon_rna.features import rnafm_embeddings
    data=pd.read_csv(ROOT/"data/real/modeling/real_model_dataset.csv")
    data=data[data.included_in_models].copy()
    embeddings=rnafm_embeddings(data.sequence.str[:1000],device="cpu")
    dest=ROOT/"data/real/modeling"
    np.save(dest/"rnafm_embeddings.npy",embeddings)
    manifest=data[["gene","transcript_id","sequence_length"]].copy()
    manifest["embedding_window_start"]=1
    manifest["embedding_window_end"]=data.sequence.str.len().clip(upper=1000)
    manifest.to_csv(dest/"rnafm_embedding_manifest.csv",index=False)
    # Comparison remains a separate future experiment until embeddings have been inspected.
    (ROOT/"results/foundation_model_status.json").write_text(json.dumps({"completed":True,"shape":list(embeddings.shape),
        "comparison_completed":False,"device":"cpu","model":"RNA-FM t12","window":"5prime <=1000 nt"},indent=2))


if __name__=="__main__":
    main()
