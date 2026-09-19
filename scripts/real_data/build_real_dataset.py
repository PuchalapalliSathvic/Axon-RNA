"""Obtain equally queried generic annotations and join the real modeling dataset."""
import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import pandas as pd
from axon_rna.real_data import ROOT, cached_json, annotation_features, encode_manual, assemble_dataset
from annotate_candidates import extract_go_terms, extract_pathways


def annotate(row, offline=False):
    out = {"gene": row.gene, "annotation_status": "failed"}
    try:
        # Same service, query fields and date/cache policy for BOTH classes.
        query = f'ensembl.gene:{row.ensembl_gene_id}' if pd.notna(row.ensembl_gene_id) else f'symbol:"{row.gene}"'
        url = "https://mygene.info/v3/query?" + urlencode(dict(q=query, species="human", fields="symbol,taxid,ensembl.gene,go,pathway", size=10))
        payload = cached_json(url, offline=offline)
        hits = [h for h in payload.get("hits", []) if h.get("taxid") == 9606]
        if len(hits) != 1:
            raise ValueError(f"Expected unique human annotation, found {len(hits)}")
        go = extract_go_terms(hits[0].get("go", {}))
        paths = extract_pathways(hits[0].get("pathway", {}))
        out.update(annotation_status="queried", go_terms=go, pathways=paths, annotation_url=url,
                   **annotation_features(go, paths))
    except Exception as error:
        out.update(annotation_error=str(error), **annotation_features("", "", available=False))
    print(f'{row.gene}: annotation {out["annotation_status"]}', flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    sequences = pd.read_csv(ROOT/"data/real/sequences/real_model_sequences.csv")
    with ThreadPoolExecutor(max_workers=3) as pool:
        knowledge = pd.DataFrame(pool.map(lambda r: annotate(r,args.offline), sequences.itertuples()))
    # Preserve old generic annotation snapshot for audit, not as a class-specific fallback.
    original = pd.read_csv(ROOT/"results/candidate_annotations_enriched.csv")
    original = original[["gene", "go_terms", "pathways"]].rename(columns={"go_terms":"prior_go_terms", "pathways":"prior_pathways"})
    knowledge = knowledge.merge(original, on="gene", how="left", validate="one_to_one")
    manual = pd.read_csv(ROOT/"results/literature_curation.csv")
    evidence_cols = [c for c in manual if c.endswith("evidence") or c == "ad_evidence_level"]
    manual_subset = manual[["gene", *evidence_cols, "pmid", "literature_note"]].copy()
    for col in evidence_cols:
        manual_subset[col+"_encoded"] = manual_subset[col].map(encode_manual)
    manual_subset = manual_subset.rename(columns={c:"manual_"+c for c in manual_subset if c != "gene"})
    knowledge = knowledge.merge(manual_subset, on="gene", how="left", validate="one_to_one")
    knowledge["manual_review_status"] = knowledge.manual_ad_evidence_level.notna().map({True:"reviewed",False:"not_reviewed"})
    auto = pd.read_csv(ROOT/"results/literature_screen.csv")
    auto_cols = [c for c in auto if c.endswith("_hits") or c in ["papers_reviewed_auto", "auto_evidence_level", "pmids_auto"]]
    auto = auto[["gene", *auto_cols]].rename(columns={c:"automated_"+c for c in auto_cols})
    knowledge = knowledge.merge(auto, on="gene", how="left", validate="one_to_one")
    knowledge.to_csv(ROOT/"data/real/modeling/real_knowledge_features.csv", index=False)
    dataset = assemble_dataset(pd.read_csv(ROOT/"data/real/modeling/selected_genes.csv"), sequences,
                               pd.read_csv(ROOT/"data/real/structures/real_structure_features.csv"), knowledge)
    dataset.to_csv(ROOT/"data/real/modeling/real_model_dataset.csv", index=False)
    print(dataset.groupby(["label", "included_in_models"]).size().to_string())


if __name__ == "__main__":
    main()
