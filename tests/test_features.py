import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from axon_rna.features import clean_rna,basic_structure_proxy
import pandas as pd

def test_clean():
    assert clean_rna('acgt-n')=='ACGU'

def test_features_shape():
    x=basic_structure_proxy(pd.Series(['ACGU','GGGG']))
    assert x.shape==(2,4)
