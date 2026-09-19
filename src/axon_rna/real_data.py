"""Auditable helpers for gene-level labels and representative RNA sequences."""
from __future__ import annotations

import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
SEED = 42
STRUCTURE_FEATURES = ["normalized_mfe", "paired_fraction"]
KNOWLEDGE_FEATURES = ["go_synaptic", "go_rna_processing", "go_neuroinflammation",
                      "go_neuronal", "go_count", "pathway_count"]


def cached_json(url, cache_dir=None, offline=False):
    cache_dir = Path(cache_dir or ROOT / "data/real/cache/http")
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / (hashlib.sha256(url.encode()).hexdigest() + ".json")
    if path.exists():
        return json.loads(path.read_text())["response"]
    if offline:
        raise RuntimeError(f"Uncached resource in offline mode: {url}")
    for attempt in range(4):
        try:
            response = requests.get(url, headers={"Accept": "application/json"}, timeout=40)
            response.raise_for_status()
            payload = response.json()
            if isinstance(payload, dict) and "error" in payload:
                raise ValueError(payload["error"])
            path.write_text(json.dumps({"url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                                        "response": payload}, indent=2))
            return payload
        except (requests.RequestException, ValueError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def choose_transcript(gene_record):
    transcripts = gene_record.get("Transcript", [])
    mane = [t for t in transcripts if any(m.get("type") == "MANE_Select" for m in t.get("MANE", []))]
    if len(mane) == 1:
        return mane[0], "MANE_Select"
    if len(mane) > 1:
        raise ValueError("Ambiguous MANE Select transcripts")
    canonical_id = str(gene_record.get("canonical_transcript", "")).split(".")[0]
    canonical = [t for t in transcripts if t["id"] == canonical_id]
    if len(canonical) != 1:
        raise ValueError("No unique canonical transcript")
    return canonical[0], "Ensembl_canonical"


def validate_rna(sequence):
    sequence = str(sequence).upper().replace("T", "U")
    if not sequence or re.search("[^ACGU]", sequence):
        raise ValueError("Empty sequence or ambiguous/non-RNA bases; sequence not altered")
    return sequence


def sequence_features(sequence):
    sequence = validate_rna(sequence)
    return {"sequence_length": len(sequence),
            "gc_content": (sequence.count("G") + sequence.count("C")) / len(sequence)}


def select_model_genes(consensus, candidates, seed=SEED):
    c = consensus.rename(columns={"transcript_id": "gene"}).copy()
    if c.gene.duplicated().any():
        raise ValueError("Consensus must contain one row per gene")
    strong = (c.n_studies == 3) & (c.direction_consistency == 1) & (c.combined_p < .01) & (c.median_log2fc.abs() >= .5)
    expected = set(candidates.get("gene", candidates.get("transcript_id")))
    if set(c.loc[strong, "gene"]) != expected:
        raise ValueError("Candidate table does not match documented strong criteria")
    pos = c.loc[strong].copy()
    ref = c.loc[(c.n_studies == 3) & (c.ad_dysregulated == 0)].sort_values("gene")
    ref = ref.sample(n=len(pos), random_state=seed)
    pos["label"], ref["label"] = 1, 0
    return pd.concat([pos, ref], ignore_index=True).sort_values("gene").reset_index(drop=True)


def parse_rnafold(output, sequence):
    match = re.search(r"^([.()]+)\s+\(\s*(-?\d+(?:\.\d+)?)\s*\)", output, re.M)
    if not match or len(match[1]) != len(sequence):
        raise ValueError("RNAfold structure missing or length mismatch")
    bracket, mfe = match[1], float(match[2])
    depth = 0
    for char in bracket:
        depth += (char == "(") - (char == ")")
        if depth < 0:
            raise ValueError("Unbalanced structure")
    if depth:
        raise ValueError("Unbalanced structure")
    return {"dot_bracket": bracket, "mfe": mfe, "normalized_mfe": mfe / len(sequence),
            "paired_fraction": (bracket.count("(") + bracket.count(")")) / len(sequence)}


def annotation_features(go_terms, pathways, available=True):
    if not available:
        return {name: np.nan for name in KNOWLEDGE_FEATURES}
    terms = {s.strip() for s in str(go_terms or "").split("; ") if s.strip()}
    paths = {s.strip() for s in str(pathways or "").split("; ") if s.strip()}
    text = " ".join(terms).lower()
    patterns = {"go_synaptic": r"synap", "go_rna_processing": r"rna (processing|splicing)|spliceosom|rrna|trna processing",
                "go_neuroinflammation": r"inflamm|microglia|immune|cytokine",
                "go_neuronal": r"neuron|axon|dendrit|synap"}
    return {**{k: int(bool(re.search(p, text))) for k, p in patterns.items()},
            "go_count": len(terms), "pathway_count": len(paths)}


def encode_manual(value):
    """Ordinal evidence strength; unknown/review absence stays missing."""
    return {"strong": 3., "moderate": 2., "limited": 1., "yes": 1., "no": 0.}.get(str(value).lower(), np.nan)


def assemble_dataset(selected, sequences, structures, knowledge):
    out = selected.copy()
    for frame, keys in [(sequences, ["gene"]), (structures, ["gene", "transcript_id"]), (knowledge, ["gene"])]:
        if frame.duplicated(keys).any():
            raise ValueError(f"Duplicate mapping keys: {keys}")
        extra = [c for c in frame if c in keys or c not in out]
        out = out.merge(frame[extra], on=keys, how="left", validate="one_to_one")
    out["included_in_models"] = (out.mapping_status == "mapped") & (out.structure_status == "folded")
    return out
