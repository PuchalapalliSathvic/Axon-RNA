"""Verify real artifacts against their source tables and save final status/provenance."""
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from datetime import datetime,timezone
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
import numpy as np
import pandas as pd
from axon_rna.real_data import ROOT
from axon_rna.discovery import build_consensus
from axon_rna.real_modeling import MODEL_FEATURES, metrics


def main():
    consensus=pd.read_csv(ROOT/"results/consensus_genes.csv")
    studies=pd.read_csv(ROOT/"data/real/processed/all_real_studies.csv")
    assert not studies.duplicated(["study_id","transcript_id"]).any()
    rebuilt=build_consensus(studies).sort_values("transcript_id").reset_index(drop=True)
    pd.testing.assert_frame_equal(rebuilt,consensus.sort_values("transcript_id").reset_index(drop=True),check_exact=False,rtol=1e-10,atol=1e-12)
    dataset=pd.read_csv(ROOT/"data/real/modeling/real_model_dataset.csv")
    pred=pd.read_csv(ROOT/"results/model_predictions.csv")
    comparison=pd.read_csv(ROOT/"results/model_comparison.csv")
    assert set(pred.gene)==set(dataset.loc[dataset.included_in_models,"gene"])
    assert not pred.gene.duplicated().any()
    assert set(pred.fold)==set(range(5))
    for name in MODEL_FEATURES:
        assert pred[name].between(0,1).all()
        result=metrics(pred.label,pred[name])
        for metric in ["roc_auc","average_precision","f1","accuracy"]:
            assert np.isclose(result[metric],comparison.set_index("model").loc[name,metric])
    joins=pred.merge(dataset[["gene","sequence"]],on="gene",validate="one_to_one")
    assert joins.groupby("sequence").fold.nunique().max()==1
    structures=pd.read_csv(ROOT/"data/real/structures/real_structure_features.csv")
    for row in structures[structures.structure_status=="folded"].itertuples():
        assert len(row.dot_bracket)==row.folded_length==row.window_end
        assert np.isclose(row.normalized_mfe,row.mfe/row.folded_length)
    figures=["top_consensus_genes","consensus_direction_heatmap","model_comparison","roc_curve",
             "precision_recall_curve","confusion_matrix","structure_feature_comparison",
             "knowledge_feature_importance","project_summary"]
    assert all((ROOT/"results/figures"/f"{name}.png").exists() for name in figures)
    tests=subprocess.run([sys.executable,"-m","pytest","-q"],cwd=ROOT,text=True,capture_output=True,check=True)
    (ROOT/"results/test_output.txt").write_text(tests.stdout+tests.stderr)
    mapping=json.loads((ROOT/"data/real/sequences/mapping_summary.json").read_text())
    foundational=json.loads((ROOT/"results/foundation_model_status.json").read_text())
    paths=[ROOT/"results"/name for name in ["consensus_genes.csv","strong_candidates.csv","model_comparison.csv","model_predictions.csv"]]
    paths+=sorted((ROOT/"data/real/processed").glob("*.csv"))
    paths+=sorted((ROOT/"data/real/sequences").glob("*.csv"))
    paths+=sorted((ROOT/"data/real/structures").glob("*.csv"))
    paths+=sorted((ROOT/"data/real/modeling").glob("*.csv"))
    paths+=[ROOT/"results"/name for name in ["candidate_annotations_enriched.csv", "literature_curation.csv", "literature_screen.csv"]]
    paths+=[ROOT/"data/real/dataset_registry.csv", ROOT/"requirements.txt", ROOT/"tests/test_real_pipeline.py"]
    paths+=sorted((ROOT/"data/real/cache/http").glob("*.json"))
    paths+=sorted((ROOT/"scripts/real_data").glob("*.py"))
    paths+=sorted((ROOT/"src/axon_rna").glob("real*.py"))
    provenance={"validated_utc":datetime.now(timezone.utc).isoformat(),"python":platform.python_version(),
                "packages":{p:importlib.metadata.version(p) for p in ["numpy","pandas","scikit-learn","scipy","matplotlib","requests","pytest","mygene"]},
                "rnafold":subprocess.check_output(["RNAfold","--version"],text=True).strip(),
                "sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (ROOT/"results/run_manifest.json").write_text(json.dumps(provenance,indent=2))
    table="| Model | ROC-AUC | Average precision | F1 | Accuracy |\n|---|---:|---:|---:|---:|\n"
    for r in comparison.itertuples():
        table+=f"| {r.model} | {r.roc_auc:.6f} | {r.average_precision:.6f} | {r.f1:.6f} | {r.accuracy:.6f} |\n"
    base=comparison.set_index("model").loc["sequence_only","roc_auc"]
    knowledge=comparison.set_index("model").loc["sequence_knowledge","roc_auc"]-base
    structure=comparison.set_index("model").loc["sequence_structure","roc_auc"]-base
    text=f"""# Final project status

Validated: {provenance['validated_utc']}. Real results only.

## COMPLETED

- Real datasets: three preserved human brain cohorts; {len(consensus):,} genes in union.
- Differential-expression screening: existing logCPM/Welch outputs retained, not upgraded to covariate-aware inference.
- Consensus: recomputed from the merged study table and matched to preserved results; {int(consensus.ad_dysregulated.sum()):,} consensus-positive, 65 strong candidates.
- Sequence mapping: {mapping['candidates']['mapped']}/65 candidates; {mapping['all_requested']['mapped']}/{mapping['all_requested']['total']} total genes. Candidate MANE Select: {mapping['candidates']['mane_select']}; canonical fallback: {mapping['candidates']['canonical_fallback']}.
- RNA structure: {int((structures.structure_status=='folded').sum())} RNAfold folds; first min(1000, length) nt at 37°C. Six SVG/PNG structure illustrations.
- Biological knowledge: equivalent fresh/cached generic MyGene queries for both classes; {int((dataset.annotation_status=='queried').sum())}/{len(dataset)} queried successfully. Prior annotations, manual evidence and automated literature evidence retained separately.
- Sequence baseline, structure model, knowledge model, and combined model: identical five-fold held-out predictions for {len(pred)} genes ({int(pred.label.sum())} strong, {int((pred.label==0).sum())} reference).
- Visualizations: all nine required real-data figures generated.
- Tests: {tests.stdout.strip().splitlines()[-1]}.
- Offline end-to-end replay is supported by saved HTTP and fold caches; source hashes and package versions are in `results/run_manifest.json`.

## Exact final metrics

Pooled out-of-fold metrics; fixed threshold 0.5. Full floating-point values and confusion matrices are in the metrics JSON files.

{table}
Knowledge ROC-AUC change vs sequence: {knowledge:+.6f}. Structure ROC-AUC change: {structure:+.6f}.
These are descriptive comparisons on one small selected sample, not statistically established improvements. Fold-level scores are retained separately.

## OPTIONAL / NOT COMPLETED

- RNA foundation model: {'extracted; see optional status' if foundational['completed'] else 'not completed; torch is available but RNA-FM (fm) is absent. No model weights downloaded and no embeddings fabricated.'}
- Foundation embeddings versus embeddings + knowledge comparison: not run.

## Remaining limitations

- Labels are gene-level. Selected transcripts are sequence representatives, not measured dysregulated isoforms.
- The 65 positive genes are a strict subset of the 1,255 consensus-positive genes. Reference genes fail the consensus rule; they are not verified biological negatives. Balanced sampling changes prevalence and precision interpretation.
- Screening uses unadjusted Fisher-combined p-values, with no combined-test FDR requirement. Individual-cohort FDR significance is not required by the strong rule.
- Existing DE screens lack covariate/cell-composition adjustment. GSE159699 has region-bearing sample names and possible repeated donors; donor independence needs an explicit audit before inferential claims.
- Five-fold gene splits are not external-cohort validation. Exact duplicate sequences are grouped, but homologous genes can still cross folds.
- Folding a 5′ window omits downstream interactions and in-vivo context. Knowledge term matching is a coarse annotation proxy.
- Candidate-only literature evidence is excluded from prediction to prevent curation coverage revealing the label. Missing/unclear reviews are not encoded as negative evidence.
- Upstream identifier harmonization is retained from supplied files; original API mapping provenance was not present. Current sequence mapping has complete cached response provenance.
- No clinical utility, biomarker, causal, therapeutic-target or novelty claims.
"""
    (ROOT/"VALIDATION.md").write_text(text)
    print(table)
    print(tests.stdout.strip())


if __name__=="__main__":
    main()
