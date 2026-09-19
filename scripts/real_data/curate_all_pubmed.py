import os
import time
import re
import pandas as pd
from Bio import Entrez, Medline

INPUT = "results/literature_curation.csv"
OUTPUT = "results/literature_screen.csv"

EMAIL = os.environ.get("NCBI_EMAIL")

if not EMAIL:
    raise RuntimeError(
        "Set your email first:\n"
        "export NCBI_EMAIL='your_email@example.com'"
    )

Entrez.email = EMAIL


def search_pubmed(gene, max_results=10):
    query = (
        f'("{gene}"[Title/Abstract]) AND '
        f'("Alzheimer Disease"[MeSH Terms] OR Alzheimer*[Title/Abstract])'
    )

    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            sort="relevance"
        )
        result = Entrez.read(handle)
        handle.close()

        return result["IdList"], int(result["Count"])

    except Exception as e:
        print(f"Search error for {gene}: {e}")
        return [], 0


def fetch_pubmed(pmids):
    if not pmids:
        return []

    try:
        handle = Entrez.efetch(
            db="pubmed",
            id=",".join(pmids),
            rettype="medline",
            retmode="text"
        )

        records = list(Medline.parse(handle))
        handle.close()

        return records

    except Exception as e:
        print("Fetch error:", e)
        return []


def contains_any(text, keywords):
    text = text.lower()
    return any(k.lower() in text for k in keywords)


def analyze_records(records):

    human_keywords = [
        "human",
        "patients",
        "patient",
        "postmortem",
        "post-mortem",
        "cerebrospinal",
        "csf",
        "plasma",
        "serum",
        "cohort",
        "brain tissue"
    ]

    mechanistic_keywords = [
        "mouse",
        "mice",
        "cell",
        "cells",
        "knockout",
        "overexpression",
        "inhibition",
        "mechanism",
        "model",
        "transgenic"
    ]

    synaptic_keywords = [
        "synapse",
        "synaptic",
        "neuronal",
        "neuron",
        "plasticity"
    ]

    tau_keywords = [
        "tau",
        "phospho-tau",
        "p-tau",
        "neurofibrillary"
    ]

    amyloid_keywords = [
        "amyloid",
        "amyloid-beta",
        "amyloid beta",
        "aβ",
        "abeta"
    ]

    inflammation_keywords = [
        "microglia",
        "microglial",
        "neuroinflammation",
        "inflammatory",
        "cytokine",
        "astrocyte",
        "astrocytic"
    ]

    human_hits = 0
    mechanism_hits = 0
    synaptic_hits = 0
    tau_hits = 0
    amyloid_hits = 0
    inflammation_hits = 0

    pmids = []
    titles = []

    for rec in records:

        title = rec.get("TI", "")
        abstract = rec.get("AB", "")
        pmid = rec.get("PMID", "")

        text = f"{title} {abstract}"

        if pmid:
            pmids.append(pmid)

        if title:
            titles.append(title)

        if contains_any(text, human_keywords):
            human_hits += 1

        if contains_any(text, mechanistic_keywords):
            mechanism_hits += 1

        if contains_any(text, synaptic_keywords):
            synaptic_hits += 1

        if contains_any(text, tau_keywords):
            tau_hits += 1

        if contains_any(text, amyloid_keywords):
            amyloid_hits += 1

        if contains_any(text, inflammation_keywords):
            inflammation_hits += 1

    return {
        "human_hits": human_hits,
        "mechanistic_hits": mechanism_hits,
        "synaptic_hits": synaptic_hits,
        "tau_hits": tau_hits,
        "amyloid_hits": amyloid_hits,
        "neuroinflammation_hits": inflammation_hits,
        "pmids_auto": ";".join(pmids),
        "top_titles": " || ".join(titles[:5])
    }


def evidence_level(total, human_hits, mechanism_hits):

    # Preliminary automated screening only.
    if total >= 8 and human_hits >= 2 and mechanism_hits >= 2:
        return "strong_candidate"

    if total >= 4 and (human_hits >= 1 or mechanism_hits >= 2):
        return "moderate_candidate"

    if total >= 1:
        return "limited_candidate"

    return "unclear"


def main():

    df = pd.read_csv(INPUT)

    results = []

    total_genes = len(df)

    for i, row in df.iterrows():

        gene = str(row["gene"])

        print(f"[{i+1}/{total_genes}] Searching {gene}...")

        pmids, total_count = search_pubmed(gene)

        time.sleep(0.4)

        records = fetch_pubmed(pmids)

        analysis = analyze_records(records)

        auto_level = evidence_level(
            total_count,
            analysis["human_hits"],
            analysis["mechanistic_hits"]
        )

        results.append({
            "gene": gene,
            "pubmed_total_hits": total_count,
            "papers_reviewed_auto": len(records),
            "auto_evidence_level": auto_level,
            **analysis
        })

        time.sleep(0.4)

    evidence = pd.DataFrame(results)

    merged = df.merge(
        evidence,
        on="gene",
        how="left"
    )

    merged.to_csv(OUTPUT, index=False)

    print("\nSaved:", OUTPUT)

    print("\nAuto evidence summary:")
    print(
        merged["auto_evidence_level"]
        .value_counts(dropna=False)
    )

    print("\nTop genes by PubMed hits:")

    show_cols = [
        "gene",
        "pubmed_total_hits",
        "auto_evidence_level",
        "human_hits",
        "mechanistic_hits",
        "tau_hits",
        "amyloid_hits",
        "synaptic_hits",
        "neuroinflammation_hits"
    ]

    print(
        merged
        .sort_values(
            "pubmed_total_hits",
            ascending=False
        )[show_cols]
        .head(30)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
