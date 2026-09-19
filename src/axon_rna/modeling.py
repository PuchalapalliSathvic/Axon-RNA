from __future__ import annotations
import json
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from .features import kmer_matrix, basic_structure_proxy


def evaluate_binary(y, prob, threshold=0.5):
    pred = (np.asarray(prob) >= threshold).astype(int)
    metrics = {
        "n": int(len(y)),
        "positive_rate": float(np.mean(y)),
        "accuracy": float(accuracy_score(y, pred)),
        "f1": float(f1_score(y, pred, zero_division=0)),
    }
    if len(np.unique(y)) > 1:
        metrics["roc_auc"] = float(roc_auc_score(y, prob))
        metrics["average_precision"] = float(average_precision_score(y, prob))
    return metrics


def cross_validated_sequence_model(df: pd.DataFrame, n_splits: int = 5, seed: int = 42):
    X, _ = kmer_matrix(df["sequence"], k=4)
    y = df["ad_dysregulated"].astype(int).to_numpy()
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    clf = LogisticRegression(max_iter=3000, class_weight="balanced")
    p = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]
    return evaluate_binary(y, p), p


def cross_validated_knowledge_model(df: pd.DataFrame, knowledge_columns: list[str], n_splits: int = 5, seed: int = 42):
    Xseq, _ = kmer_matrix(df["sequence"], k=4)
    struct = basic_structure_proxy(df["sequence"]).to_numpy(float)
    know = df[knowledge_columns].fillna(0).to_numpy(float)
    dense = np.hstack([struct, know])
    dense = StandardScaler().fit_transform(dense)
    X = sparse.hstack([Xseq, sparse.csr_matrix(dense)]).tocsr()
    y = df["ad_dysregulated"].astype(int).to_numpy()
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    clf = LogisticRegression(max_iter=3000, class_weight="balanced")
    p = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")[:, 1]
    return evaluate_binary(y, p), p


def low_data_curve(df: pd.DataFrame, knowledge_columns: list[str] | None = None, seed: int = 42):
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score
    fractions = [0.1, 0.2, 0.5, 1.0]
    train, test = train_test_split(df, test_size=0.3, stratify=df["ad_dysregulated"], random_state=seed)
    out=[]
    for frac in fractions:
        if frac == 1.0:
            sample = train.copy()
        else:
            parts=[]
            for _, g in train.groupby("ad_dysregulated"):
                n=max(2, int(round(len(g)*frac)))
                n=min(n, len(g))
                parts.append(g.sample(n=n, random_state=seed))
            sample=pd.concat(parts, ignore_index=True)
        Xtr, Xte, _ = kmer_matrix(sample["sequence"], test["sequence"], k=4)
        if knowledge_columns:
            tr_dense=np.hstack([basic_structure_proxy(sample["sequence"]).to_numpy(float), sample[knowledge_columns].fillna(0).to_numpy(float)])
            te_dense=np.hstack([basic_structure_proxy(test["sequence"]).to_numpy(float), test[knowledge_columns].fillna(0).to_numpy(float)])
            scaler=StandardScaler().fit(tr_dense)
            Xtr=sparse.hstack([Xtr, sparse.csr_matrix(scaler.transform(tr_dense))]).tocsr()
            Xte=sparse.hstack([Xte, sparse.csr_matrix(scaler.transform(te_dense))]).tocsr()
        clf=LogisticRegression(max_iter=3000,class_weight="balanced").fit(Xtr, sample["ad_dysregulated"])
        p=clf.predict_proba(Xte)[:,1]
        auc=float(roc_auc_score(test["ad_dysregulated"],p)) if test["ad_dysregulated"].nunique()>1 else float("nan")
        out.append({"train_fraction":frac,"n_train":len(sample),"roc_auc":auc})
    return out
