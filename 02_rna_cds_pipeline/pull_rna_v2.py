import csv
import sys
import json
import time
import urllib.request
import urllib.error
import argparse

csv.field_size_limit(sys.maxsize)

SRC = "/u/project/kappel/asharma/RBP/merged/table_260823.csv"

UNIPROT_ACC_URL = "https://rest.uniprot.org/uniprotkb/accessions"
ENSEMBL_SEQ = "https://rest.ensembl.org/sequence/id?type=cdna;content-type=application/json"


def http_get(url, timeout=60):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def fetch_ensembl_xrefs(accessions, batch_size=100):
    """Returns dict: accession -> list of (ensembl_transcript_id_versioned, isoformId_or_None)."""
    out = {}
    ids = list(accessions)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        url = f"{UNIPROT_ACC_URL}?accessions={','.join(batch)}&fields=accession,xref_ensembl"
        body = None
        for attempt in range(4):
            try:
                body = http_get(url, timeout=90)
                break
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(8)
                    continue
                print(f"  uniprot batch error {e.code} at offset {i}", file=sys.stderr)
                body = ""
                break
            except Exception as e:
                print(f"  uniprot batch exception at offset {i}: {e}", file=sys.stderr)
                time.sleep(3)
                body = ""
        if not body:
            continue
        try:
            d = json.loads(body)
        except json.JSONDecodeError:
            continue
        for entry in d.get("results", []):
            acc = entry.get("primaryAccession")
            xrefs = []
            for xref in entry.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "Ensembl":
                    xrefs.append((xref["id"], xref.get("isoformId")))
            out[acc] = xrefs
        if (i // batch_size) % 10 == 0:
            print(f"  ... UniProt xref lookup {i+len(batch)}/{len(ids)} done", flush=True)
        time.sleep(0.3)
    return out


def fetch_sequences(bare_enst_ids, batch_size=50):
    seqs = {}
    ids = list(bare_enst_ids)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        payload = json.dumps({"ids": batch}).encode()
        req = urllib.request.Request(
            ENSEMBL_SEQ, data=payload,
            headers={"Content-Type": "application/json"},
        )
        body = None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=90) as resp:
                    body = resp.read().decode()
                break
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(8)
                    continue
                body = e.read().decode()
                print(f"  ensembl batch error {e.code}: {body[:200]}", file=sys.stderr)
                body = "[]"
                break
        try:
            results = json.loads(body)
        except json.JSONDecodeError:
            results = []
        if isinstance(results, dict):
            results = [results]
        for r in results:
            if isinstance(r, dict) and "id" in r and "seq" in r:
                seqs[r["id"]] = r["seq"]
        if (i // batch_size) % 20 == 0:
            print(f"  ... Ensembl seq fetch {i+len(batch)}/{len(ids)} done", flush=True)
        time.sleep(0.3)
    return seqs


def pick_transcript(acc, xrefs):
    """Prefer a transcript whose isoformId matches the canonical accession
    (isoformId is None, meaning single-isoform entry, or equals '{acc}-1').
    Falls back to the first available transcript if no canonical-tagged one exists."""
    canonical = [t for t, iso in xrefs if iso is None or iso == f"{acc}-1"]
    if canonical:
        return canonical[0], "canonical_isoform_match"
    if xrefs:
        return xrefs[0][0], "fallback_no_canonical_tag"
    return None, "no_ensembl_xref"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    with open(SRC, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)
            if args.limit and len(rows) >= args.limit:
                break

    uids = [r["uniprot_accession"].strip() for r in rows if r["uniprot_accession"].strip()]
    print(f"Rows: {len(rows)}; distinct uniprot IDs: {len(set(uids))}")

    print("Fetching Ensembl transcript cross-references (with isoform tags) from UniProt...")
    xref_map = fetch_ensembl_xrefs(set(uids))
    print(f"  got xref data for {len(xref_map)} / {len(set(uids))} accessions")

    chosen = {}  # uid -> (versioned_enst, method)
    method_counts = {}
    for uid in set(uids):
        xrefs = xref_map.get(uid, [])
        enst, method = pick_transcript(uid, xrefs)
        chosen[uid] = (enst, method)
        method_counts[method] = method_counts.get(method, 0) + 1

    print("Selection method breakdown:")
    for m, c in method_counts.items():
        print(f"  {m}: {c}")

    bare_ensts = sorted(set(e.split(".")[0] for e, _ in chosen.values() if e))
    print(f"Fetching {len(bare_ensts)} distinct transcript sequences from Ensembl...")
    seqs = fetch_sequences(bare_ensts)
    print(f"  got sequences for {len(seqs)} / {len(bare_ensts)} transcripts")

    out_fields = fieldnames + ["ensembl_transcript_used", "transcript_selection_method", "rna_sequence"]
    n_with_rna = 0
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        for row in rows:
            uid = row["uniprot_accession"].strip()
            enst, method = chosen.get(uid, (None, "no_ensembl_xref"))
            seq = seqs.get(enst.split(".")[0], "") if enst else ""
            row["ensembl_transcript_used"] = enst or ""
            row["transcript_selection_method"] = method
            row["rna_sequence"] = seq
            writer.writerow(row)
            if seq:
                n_with_rna += 1

    print(f"Wrote {len(rows)} rows to {args.out}; {n_with_rna} have an rna_sequence")


if __name__ == "__main__":
    main()
