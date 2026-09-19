import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import pandas as pd
from axon_rna.discovery import build_consensus

def test_consensus_reproducible():
    df=pd.DataFrame({
      'study_id':['a','b','a','b'],
      'transcript_id':['x','x','y','y'],
      'log2fc':[1.0,.9,.1,-.1],
      'pvalue':[.001,.002,.8,.7],
      'padj':[.01,.02,.9,.9]})
    out=build_consensus(df).set_index('transcript_id')
    assert out.loc['x','ad_dysregulated']==1
    assert out.loc['y','ad_dysregulated']==0
