import argparse, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from axon_rna.discovery import build_consensus
p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--output',required=True); a=p.parse_args()
df=pd.read_csv(a.input); out=build_consensus(df); Path(a.output).parent.mkdir(parents=True,exist_ok=True); out.to_csv(a.output,index=False)
print(out.head(10).to_string(index=False)); print('\nSaved',a.output,'positive=',int(out.ad_dysregulated.sum()),'/',len(out))
