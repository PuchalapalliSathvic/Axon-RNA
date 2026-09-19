"""Fold-local feature fitting and identical held-out splits for real-data comparisons."""
import hashlib

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.preprocessing import StandardScaler

from .modeling import evaluate_binary
from .real_data import STRUCTURE_FEATURES, KNOWLEDGE_FEATURES, SEED

MODEL_FEATURES = {
    "sequence_only": [],
    "sequence_structure": STRUCTURE_FEATURES,
    "sequence_knowledge": KNOWLEDGE_FEATURES,
    "sequence_structure_knowledge": STRUCTURE_FEATURES + KNOWLEDGE_FEATURES,
}


def make_estimator(numeric_columns):
    allowed = set(STRUCTURE_FEATURES + KNOWLEDGE_FEATURES)
    if set(numeric_columns) - allowed:
        raise ValueError("Only predefined sequence-independent input features are permitted")
    transformers = [("sequence", TfidfVectorizer(analyzer="char", ngram_range=(4,4), lowercase=False), "sequence")]
    if numeric_columns:
        transformers.append(("numeric", make_pipeline(SimpleImputer(strategy="median", keep_empty_features=True), StandardScaler()), numeric_columns))
    return Pipeline([("features", ColumnTransformer(transformers, remainder="drop")),
                     ("classifier", LogisticRegression(C=1.0, max_iter=3000, class_weight="balanced", random_state=SEED))])


def make_splits(frame, n_splits=5):
    if frame.gene.duplicated().any() or frame.transcript_id.duplicated().any():
        raise ValueError("Modeling requires unique genes and representative transcripts")
    groups = frame.sequence.map(lambda s: hashlib.sha256(s.encode()).hexdigest())
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    splits = list(cv.split(frame, frame.label, groups))
    for train, test in splits:
        if set(groups.iloc[train]) & set(groups.iloc[test]):
            raise ValueError("Identical sequences cross folds")
        if frame.iloc[train].label.nunique() != 2 or frame.iloc[test].label.nunique() != 2:
            raise ValueError("Every training and test fold must contain both classes")
    return splits


def metrics(y, probability):
    return {**evaluate_binary(y, probability),
            "confusion_matrix": confusion_matrix(y, np.asarray(probability)>=.5, labels=[0,1]).tolist()}
