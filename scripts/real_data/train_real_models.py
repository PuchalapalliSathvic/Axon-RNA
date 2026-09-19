"""Fixed five-fold comparisons. Labels, DE statistics and literature are never features."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import numpy as np
import pandas as pd
from sklearn.inspection import permutation_importance
from axon_rna.real_data import ROOT, KNOWLEDGE_FEATURES, SEED
from axon_rna.real_modeling import MODEL_FEATURES, make_estimator, make_splits, metrics


def main():
    dataset = pd.read_csv(ROOT/"data/real/modeling/real_model_dataset.csv")
    frame = dataset[dataset.included_in_models].reset_index(drop=True)
    if len(frame)<20 or frame.label.value_counts().min()<5:
        raise SystemExit("Insufficient real mapped/folded genes for five-fold modeling")
    splits = make_splits(frame)
    predictions = frame[["gene", "transcript_id", "label"]].copy()
    comparisons, importance, fold_rows = [], [], []
    for fold, (train, test) in enumerate(splits):
        predictions.loc[test,"fold"] = fold
    for name, columns in MODEL_FEATURES.items():
        probability = np.full(len(frame), np.nan)
        for fold, (train,test) in enumerate(splits):
            estimator = make_estimator(columns)
            estimator.fit(frame.iloc[train], frame.label.iloc[train])
            probability[test] = estimator.predict_proba(frame.iloc[test])[:,1]
            fold_rows.append(dict(model=name, fold=fold, **metrics(frame.label.iloc[test], probability[test])))
            if name == "sequence_knowledge":
                # Permute ONLY generic knowledge columns on held-out genes.
                x = frame.iloc[test][["sequence",*KNOWLEDGE_FEATURES]].copy()
                perm = permutation_importance(estimator, x, frame.label.iloc[test],
                                              scoring="roc_auc", n_repeats=20, random_state=SEED+fold)
                for i,col in enumerate(x.columns):
                    if col in KNOWLEDGE_FEATURES:
                        importance.append(dict(fold=fold, feature=col, auc_decrease=perm.importances_mean[i],
                                               permutation_sd=perm.importances_std[i]))
        if np.isnan(probability).any():
            raise AssertionError("Missing held-out prediction")
        predictions[name] = probability
        score = metrics(frame.label, probability)
        comparisons.append(dict(model=name, **{k:v for k,v in score.items() if k != "confusion_matrix"}))
        report = dict(model=name, **score, cv="StratifiedGroupKFold, 5 splits; exact-sequence groups", seed=SEED,
                      threshold=.5, numeric_features=columns, tfidf="4-mer, fit within training fold", C=1.0,
                      interpretation="Predicts membership in a selected strong gene-level consensus subset; exploratory")
        filename = "sequence_baseline_metrics.json" if name == "sequence_only" else f"{name}_metrics.json"
        (ROOT/"results"/filename).write_text(json.dumps(report, indent=2))
    predictions["fold"] = predictions.fold.astype(int)
    predictions.to_csv(ROOT/"results/model_predictions.csv", index=False)
    predictions[["gene","transcript_id","label","fold","sequence_only"]].rename(columns={"sequence_only":"probability"}).to_csv(ROOT/"results/sequence_baseline_predictions.csv", index=False)
    comparison = pd.DataFrame(comparisons)
    comparison.to_csv(ROOT/"results/model_comparison.csv", index=False)
    pd.DataFrame(fold_rows).to_json(ROOT/"results/rnafold_summary.json", orient="records", indent=2)
    pd.DataFrame(importance).to_csv(ROOT/"results/knowledge_feature_importance.csv", index=False)
    predictions[["gene","transcript_id","label","fold"]].to_csv(ROOT/"data/real/modeling/cv_splits.csv", index=False)
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
