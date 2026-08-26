import csv
import sys
import json
import time
import urllib.request
import urllib.error

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/id?expand=1"
ENSEMBL_SEQ = "https://rest.ensembl.org/sequence/id?type=cdna;content-type=application/json"
ENSEMBL_SEQ_CDS = "https://rest.ensembl.org/sequence/id?type=cds;content-type=application/json"


def http_post_json(url, payload, timeout=90):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def lookup_genes(bare_ensgs, batch_size=50):
    """Returns dict: ENSG -> chosen transcript ID (translated, preferring canonical) or None"""
    result = {}
    ids = list(bare_ensgs)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        body = None
        for attempt in range(6):
            try:
                body = http_post_json(ENSEMBL_LOOKUP, {"ids": batch})
                break
            except urllib.error.HTTPError as e:
                if e.code in (429, 500, 502, 503, 504):
                    time.sleep(5 * (attempt + 1))
                    continue
                print(f"  lookup batch error {e.code} at offset {i}", file=sys.stderr)
                body = None
                break
            except Exception as e:
                print(f"  lookup batch exception at offset {i}: {e}", file=sys.stderr)
                time.sleep(5)
                continue
        if not body:
            continue
        try:
            d = json.loads(body)
        except json.JSONDecodeError:
            continue
        for gid in batch:
            info = d.get(gid)
            if not info:
                result[gid] = None
                continue
            canonical = info.get("canonical_transcript", "")
            canonical_bare = canonical.split(".")[0] if canonical else None
            chosen = None
            transcripts = info.get("Transcript", [])
            # prefer the canonical transcript if it has a real translation
            for t in transcripts:
                if t.get("id") == canonical_bare and t.get("Translation"):
                    chosen = t["id"]
                    break
            if not chosen:
                for t in transcripts:
                    if t.get("Translation"):
                        chosen = t["id"]
                        break
            result[gid] = chosen
        if (i // batch_size) % 10 == 0:
            print(f"  ... gene lookup {i+len(batch)}/{len(ids)} done", flush=True)
        time.sleep(0.4)
    return result


def fetch_sequences(bare_ids, batch_size=50, url=ENSEMBL_SEQ):
    seqs = {}
    ids = list(bare_ids)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        payload = json.dumps({"ids": batch}).encode()
        body = "[]"
        for attempt in range(6):
            try:
                req = urllib.request.Request(
                    url, data=payload, headers={"Content-Type": "application/json"})
                with urllib.request.urlopen(req, timeout=90) as resp:
                    body = resp.read().decode()
                break
            except urllib.error.HTTPError as e:
                if e.code in (429, 500, 502, 503, 504):
                    time.sleep(5 * (attempt + 1))
                    continue
                body = e.read().decode()
                print(f"  seq batch error {e.code}: {body[:200]}", file=sys.stderr)
                body = "[]"
                break
            except Exception as e:
                print(f"  seq batch exception (attempt {attempt+1}): {e}", file=sys.stderr)
                time.sleep(5)
                continue
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
            print(f"  ... seq fetch {i+len(batch)}/{len(ids)} done", flush=True)
        time.sleep(0.3)
    return seqs


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="only try first N candidates (testing)")
    ap.add_argument("--out", default=TARGET)
    ap.add_argument("--dry-run", action="store_true", help="don't write output, just report")
    args = ap.parse_args()

    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    candidates = {}  # uid -> bare ENSG
    for row in rows:
        if row.get("transcript_selection_method") != "no_ensembl_xref":
            continue
        g = (row.get("ensembl_gene") or "").strip()
        if not g:
            continue
        first_ensg = g.split(";")[0].strip().split(".")[0]
        if first_ensg:
            candidates[row["uniprot_accession"]] = first_ensg
            if args.limit and len(candidates) >= args.limit:
                break

    print(f"Candidates with an ensembl_gene to try: {len(candidates)}")

    unique_ensgs = set(candidates.values())
    print(f"Distinct genes to look up: {len(unique_ensgs)}")

    gene_to_transcript = lookup_genes(unique_ensgs)
    n_found = sum(1 for v in gene_to_transcript.values() if v)
    print(f"Genes with a translated transcript found: {n_found} / {len(unique_ensgs)}")

    transcripts_needed = sorted(set(v for v in gene_to_transcript.values() if v))
    print(f"Fetching mRNA sequences for {len(transcripts_needed)} transcripts...")
    seqs = fetch_sequences(transcripts_needed, url=ENSEMBL_SEQ)
    print(f"Got mRNA sequences for {len(seqs)} / {len(transcripts_needed)}")

    print(f"Fetching CDS for {len(transcripts_needed)} transcripts...")
    cds_seqs = fetch_sequences(transcripts_needed, url=ENSEMBL_SEQ_CDS)
    print(f"Got CDS for {len(cds_seqs)} / {len(transcripts_needed)}")

    out_fields = fieldnames
    n_filled = 0
    for row in rows:
        uid = row["uniprot_accession"]
        if uid not in candidates:
            continue
        ensg = candidates[uid]
        transcript = gene_to_transcript.get(ensg)
        if not transcript:
            continue
        seq = seqs.get(transcript)
        if not seq:
            continue
        row["ensembl_transcript_used"] = transcript
        row["transcript_selection_method"] = "gene_level_fallback"
        row["rna_sequence"] = seq
        row["cds_sequence"] = cds_seqs.get(transcript, "")
        n_filled += 1

    if not args.dry_run:
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=out_fields)
            writer.writeheader()
            writer.writerows(rows)

    print(f"Filled {n_filled} rows via gene-level fallback (dry_run={args.dry_run})")


if __name__ == "__main__":
    main()
