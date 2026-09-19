# Final project status

Validated: 2026-09-19T21:25:42.900367+00:00. Real results only.

## COMPLETED

- Real datasets: three preserved human brain cohorts; 36,907 genes in union.
- Differential-expression screening: existing logCPM/Welch outputs retained, not upgraded to covariate-aware inference.
- Consensus: recomputed from the merged study table and matched to preserved results; 1,255 consensus-positive, 65 strong candidates.
- Sequence mapping: 65/65 candidates; 130/130 total genes. Candidate MANE Select: 64; canonical fallback: 1.
- RNA structure: 130 RNAfold folds; first min(1000, length) nt at 37°C. Six SVG/PNG structure illustrations.
- Biological knowledge: equivalent fresh/cached generic MyGene queries for both classes; 130/130 queried successfully. Prior annotations, manual evidence and automated literature evidence retained separately.
- Sequence baseline, structure model, knowledge model, and combined model: identical five-fold held-out predictions for 130 genes (65 strong, 65 reference).
- Visualizations: all nine required real-data figures generated.
- Tests: 12 passed in 0.84s.
- Offline end-to-end replay is supported by saved HTTP and fold caches; source hashes and package versions are in `results/run_manifest.json`.

## Exact final metrics

Pooled out-of-fold metrics; fixed threshold 0.5. Full floating-point values and confusion matrices are in the metrics JSON files.

| Model | ROC-AUC | Average precision | F1 | Accuracy |
|---|---:|---:|---:|---:|
| sequence_only | 0.507219 | 0.518625 | 0.518519 | 0.500000 |
| sequence_structure | 0.567337 | 0.593838 | 0.558140 | 0.561538 |
| sequence_knowledge | 0.521420 | 0.526670 | 0.461538 | 0.515385 |
| sequence_structure_knowledge | 0.547929 | 0.541395 | 0.504202 | 0.546154 |

Knowledge ROC-AUC change vs sequence: +0.014201. Structure ROC-AUC change: +0.060118.
These are descriptive comparisons on one small selected sample, not statistically established improvements. Fold-level scores are retained separately.

## OPTIONAL / NOT COMPLETED

- RNA foundation model: not completed; torch is available but RNA-FM (fm) is absent. No model weights downloaded and no embeddings fabricated.
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
