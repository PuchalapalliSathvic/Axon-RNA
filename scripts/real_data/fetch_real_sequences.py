"""Fetch versioned human representative transcripts; retain failures and HTTP provenance."""
import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import pandas as pd
from axon_rna.real_data import ROOT, cached_json, choose_transcript, validate_rna, select_model_genes


def fetch_gene(gene, offline=False):
    row = dict(gene=gene, ensembl_gene_id="", transcript_id="", transcript_selection_method="",
               sequence="", sequence_length=0, mapping_status="failed", mapping_error="")
    try:
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{quote(gene)}?expand=1&mane=1&content-type=application%2Fjson"
        record = cached_json(url, offline=offline)
        if record.get("species") != "homo_sapiens" or record.get("object_type") != "Gene":
            raise ValueError("Not a human gene record")
        transcript, method = choose_transcript(record)
        tid = transcript["id"]
        row.update(ensembl_gene_id=record["id"], transcript_id=f'{tid}.{transcript["version"]}',
                   transcript_selection_method=method, assembly=record.get("assembly_name"),
                   resolved_symbol=record.get("display_name"), gene_biotype=record.get("biotype"), lookup_url=url)
        seq_url = f"https://rest.ensembl.org/sequence/id/{tid}?type=cdna&content-type=application%2Fjson"
        payload = cached_json(seq_url, offline=offline)
        if payload.get("id", "").split(".")[0] != tid:
            raise ValueError("Sequence identifier mismatch")
        seq = validate_rna(payload["seq"])
        expected_length = sum(e["end"] - e["start"] + 1 for e in transcript["Exon"])
        if len(seq) != expected_length:
            raise ValueError("cDNA length differs from selected transcript exons")
        row.update(sequence=seq, sequence_length=len(seq), mapping_status="mapped", sequence_url=seq_url)
    except Exception as error:
        row["mapping_error"] = str(error)
    print(f'{gene}: {row["mapping_status"]} {row.get("transcript_id", "")}', flush=True)
    return row


def save_sequences(frame, stem):
    stem.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(stem.with_suffix(".csv"), index=False)
    with stem.with_suffix(".fasta").open("w") as handle:
        for r in frame.loc[frame.mapping_status == "mapped"].itertuples():
            handle.write(f">{r.gene}|{r.transcript_id}|{r.transcript_selection_method}\n")
            handle.write("\n".join(r.sequence[i:i+80] for i in range(0, len(r.sequence), 80)) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--include-reference", action="store_true")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    candidates = pd.read_csv(ROOT / "results/strong_candidates.csv")
    genes = candidates.transcript_id.tolist()
    if args.include_reference:
        selected = select_model_genes(pd.read_csv(ROOT / "results/consensus_genes.csv"), candidates)
        (ROOT / "data/real/modeling").mkdir(parents=True, exist_ok=True)
        selected.to_csv(ROOT / "data/real/modeling/selected_genes.csv", index=False)
        genes += sorted(set(selected.gene) - set(genes))
    with ThreadPoolExecutor(max_workers=4) as executor:
        frame = pd.DataFrame(executor.map(lambda g: fetch_gene(g, args.offline), genes))
    candidate_frame = frame[frame.gene.isin(candidates.transcript_id)]
    save_sequences(candidate_frame, ROOT / "data/real/sequences/real_candidate_sequences")
    if args.include_reference:
        save_sequences(frame, ROOT / "data/real/sequences/real_model_sequences")
    summary = {}
    for name, df in [("candidates", candidate_frame), ("all_requested", frame)]:
        ok = df[df.mapping_status == "mapped"]
        summary[name] = dict(total=len(df), mapped=len(ok), failed=len(df)-len(ok),
                             mane_select=int((ok.transcript_selection_method == "MANE_Select").sum()),
                             canonical_fallback=int((ok.transcript_selection_method == "Ensembl_canonical").sum()))
    (ROOT / "data/real/sequences/mapping_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
