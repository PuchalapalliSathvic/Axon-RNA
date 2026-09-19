import pandas as pd
import mygene
from pathlib import Path


INPUT = "results/candidate_annotations.csv"
OUTPUT = "results/candidate_annotations_enriched.csv"


def extract_go_terms(go_block):
    if not isinstance(go_block, dict):
        return ""

    terms = []

    for category in ["BP", "MF", "CC"]:
        values = go_block.get(category, [])

        if isinstance(values, dict):
            values = [values]

        if isinstance(values, list):
            for item in values:
                if isinstance(item, dict):
                    term = item.get("term")
                    if term:
                        terms.append(term)

    return "; ".join(sorted(set(terms)))


def extract_pathways(pathway_block):
    if not isinstance(pathway_block, dict):
        return ""

    pathways = []

    for source, values in pathway_block.items():

        if isinstance(values, dict):
            values = [values]

        if isinstance(values, list):
            for item in values:
                if isinstance(item, dict):
                    name = item.get("name")
                    if name:
                        pathways.append(f"{source}: {name}")

    return "; ".join(sorted(set(pathways)))


def main():

    Path("results").mkdir(exist_ok=True)

    df = pd.read_csv(INPUT)

    mg = mygene.MyGeneInfo()

    genes = df["gene"].dropna().astype(str).tolist()

    print(f"Querying annotations for {len(genes)} genes...")

    results = mg.querymany(
        genes,
        scopes="symbol",
        fields=[
            "symbol",
            "name",
            "type_of_gene",
            "alias",
            "summary",
            "go",
            "pathway"
        ],
        species="human",
        as_dataframe=False
    )

    annotation = {}

    for result in results:

        gene = result.get("query")

        if result.get("notfound"):
            continue

        aliases = result.get("alias", "")

        if isinstance(aliases, list):
            aliases = "; ".join(aliases)

        annotation[gene] = {
            "gene_name": result.get("name", ""),
            "gene_type": result.get("type_of_gene", ""),
            "aliases": aliases,
            "gene_summary": result.get("summary", ""),
            "go_terms": extract_go_terms(result.get("go", {})),
            "pathways": extract_pathways(result.get("pathway", {})),
        }

    annotation_df = pd.DataFrame.from_dict(
        annotation,
        orient="index"
    ).reset_index()

    annotation_df = annotation_df.rename(columns={"index": "gene"})

    enriched = df.merge(
        annotation_df,
        on="gene",
        how="left"
    )

    # Empty fields for the literature-curation phase.
    enriched["known_ad_evidence"] = ""
    enriched["ad_evidence_type"] = ""
    enriched["amyloid_related"] = ""
    enriched["tau_related"] = ""
    enriched["neuroinflammation_related"] = ""
    enriched["synaptic_related"] = ""
    enriched["rna_processing_related"] = ""
    enriched["literature_notes"] = ""

    enriched.to_csv(OUTPUT, index=False)

    print(f"\nSaved: {OUTPUT}")
    print(f"Genes: {len(enriched)}")
    print(
        enriched[
            [
                "gene",
                "median_log2fc",
                "gene_name",
                "gene_type"
            ]
        ].head(20).to_string(index=False)
    )


if __name__ == "__main__":
    main()
