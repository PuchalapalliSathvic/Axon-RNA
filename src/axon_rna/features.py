from __future__ import annotations
import re
import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer


def clean_rna(seq: str) -> str:
    s = re.sub(r"[^ACGTU]", "", str(seq).upper()).replace("T", "U")
    return s


def basic_structure_proxy(seqs: pd.Series) -> pd.DataFrame:
    """Cheap interpretable proxies; replace/augment with ViennaRNA later."""
    out = []
    for raw in seqs:
        s = clean_rna(raw)
        n = max(len(s), 1)
        gc = (s.count("G") + s.count("C")) / n
        au = (s.count("A") + s.count("U")) / n
        ug = sum(1 for i in range(len(s)-1) if s[i:i+2] == "UG") / max(n-1, 1)
        out.append([len(s), gc, au, ug])
    return pd.DataFrame(out, columns=["seq_length", "gc_fraction", "au_fraction", "ug_dinucleotide_fraction"])


def kmer_matrix(train_sequences, test_sequences=None, k: int = 4, max_features: int = 2048):
    analyzer = "char"
    vec = TfidfVectorizer(analyzer=analyzer, ngram_range=(k, k), max_features=max_features, lowercase=False)
    train = [clean_rna(s) for s in train_sequences]
    Xtr = vec.fit_transform(train)
    if test_sequences is None:
        return Xtr, vec
    Xte = vec.transform([clean_rna(s) for s in test_sequences])
    return Xtr, Xte, vec


def rnafm_embeddings(sequences, device: str = "cuda") -> np.ndarray:
    """Optional RNA-FM embedding backend.

    Requires torch + rna-fm. Mean pools layer-12 nucleotide embeddings.
    """
    try:
        import torch
        import fm
    except ImportError as e:
        raise RuntimeError("Install optional dependencies: torch and rna-fm") from e
    model, alphabet = fm.pretrained.rna_fm_t12()
    model = model.to(device).eval()
    converter = alphabet.get_batch_converter()
    pooled = []
    for i, seq in enumerate(sequences):
        data = [(str(i), clean_rna(seq))]
        _, _, tokens = converter(data)
        tokens = tokens.to(device)
        with torch.no_grad():
            rep = model(tokens, repr_layers=[12])["representations"][12]
        # drop BOS/EOS where possible
        arr = rep[0, 1:1+len(clean_rna(seq))].mean(0).detach().cpu().numpy()
        pooled.append(arr)
    return np.vstack(pooled)
