import csv
import sys
import json
import time
import urllib.request
import urllib.error
import argparse

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
ENSEMBL_SEQ = "https://rest.ensembl.org/sequence/id?type=cds;content-type=application/json"


def fetch_sequences(bare_ids, batch_size=50):
    seqs = {}
    ids = list(bare_ids)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        payload = json.dumps({"ids": batch}).encode()
        body = "[]"
        for attempt in range(6):
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
        if (i // batch_size) % 20 == 0:
            print(f"  ... {i+len(batch)}/{len(ids)} done", flush=True)
        time.sleep(0.3)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=TARGET)
    args = ap.parse_args()

    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)
            if args.limit and len(rows) >= args.limit:
                break

    bare_ids = sorted(set(
        r["ensembl_transcript_used"].split(".")[0]
        for r in rows if r.get("ensembl_transcript_used")
    ))
    print(f"Rows: {len(rows)}; distinct transcripts to fetch CDS for: {len(bare_ids)}")

    seqs = fetch_sequences(bare_ids)
    print(f"Got CDS for {len(seqs)} / {len(bare_ids)} transcripts")

    out_fields = fieldnames if "cds_sequence" in fieldnames else fieldnames + ["cds_sequence"]
    n_with_cds = 0
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        for row in rows:
            enst = row.get("ensembl_transcript_used", "")
            seq = seqs.get(enst.split(".")[0], "") if enst else ""
            row["cds_sequence"] = seq
            writer.writerow(row)
            if seq:
                n_with_cds += 1

    print(f"Wrote {len(rows)} rows; {n_with_cds} have a cds_sequence")


if __name__ == "__main__":
    main()
