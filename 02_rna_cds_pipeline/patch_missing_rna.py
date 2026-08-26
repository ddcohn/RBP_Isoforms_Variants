import csv
import sys
import json
import time
import urllib.request
import urllib.error

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
ENSEMBL_SEQ = "https://rest.ensembl.org/sequence/id?type=cdna;content-type=application/json"


def fetch_sequences(bare_ids, batch_size=50):
    seqs = {}
    ids = list(bare_ids)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        payload = json.dumps({"ids": batch}).encode()
        body = "[]"
        for attempt in range(8):
            try:
                req = urllib.request.Request(
                    ENSEMBL_SEQ, data=payload,
                    headers={"Content-Type": "application/json"},
                )
                with urllib.request.urlopen(req, timeout=90) as resp:
                    body = resp.read().decode()
                break
            except urllib.error.HTTPError as e:
                if e.code in (429, 503):
                    wait = 5 * (attempt + 1)
                    print(f"  HTTP {e.code}, retrying in {wait}s (attempt {attempt+1})", file=sys.stderr)
                    time.sleep(wait)
                    continue
                body = e.read().decode()
                print(f"  batch error {e.code}: {body[:200]}", file=sys.stderr)
                body = "[]"
                break
            except Exception as e:
                print(f"  batch exception (attempt {attempt+1}): {e}", file=sys.stderr)
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
        print(f"  ... {i+len(batch)}/{len(ids)} done, {len(seqs)} sequences so far", flush=True)
        time.sleep(0.3)
    return seqs


def main():
    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    missing = [r for r in rows if r.get("ensembl_transcript_used") and not r.get("rna_sequence")]
    print(f"Total rows: {len(rows)}; missing rna_sequence despite having a transcript: {len(missing)}")

    bare_ids = sorted(set(r["ensembl_transcript_used"].split(".")[0] for r in missing))
    print(f"Distinct transcripts to re-fetch: {len(bare_ids)}")

    seqs = fetch_sequences(bare_ids)
    print(f"Recovered {len(seqs)} / {len(bare_ids)} sequences")

    n_filled = 0
    for row in rows:
        if row.get("ensembl_transcript_used") and not row.get("rna_sequence"):
            bare = row["ensembl_transcript_used"].split(".")[0]
            seq = seqs.get(bare, "")
            if seq:
                row["rna_sequence"] = seq
                n_filled += 1

    with open(TARGET, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    still_missing = sum(1 for r in rows if r.get("ensembl_transcript_used") and not r.get("rna_sequence"))
    print(f"Filled {n_filled} rows. Still missing (had transcript, no sequence): {still_missing}")
    total_with_rna = sum(1 for r in rows if r.get("rna_sequence"))
    print(f"Total rows with rna_sequence now: {total_with_rna} / {len(rows)}")


if __name__ == "__main__":
    main()
