import csv
import sys
import json
import time
import urllib.request
import urllib.parse
import argparse

csv.field_size_limit(sys.maxsize)

SRC = "/u/project/kappel/asharma/RBP/merged/table_260823.csv"

UNIPROT_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_STATUS = "https://rest.uniprot.org/idmapping/status/{}"
UNIPROT_STREAM = "https://rest.uniprot.org/idmapping/stream/{}?format=tsv"
ENSEMBL_SEQ = "https://rest.ensembl.org/sequence/id?type=cdna;content-type=application/json"


def http_post(url, data_dict):
    data = urllib.parse.urlencode(data_dict).encode()
    req = urllib.request.Request(url, data=data)
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode()


def http_get(url):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode()


def submit_mapping(ids):
    resp = http_post(UNIPROT_RUN, {
        "ids": ",".join(ids),
        "from": "UniProtKB_AC-ID",
        "to": "Ensembl_Transcript",
    })
    return json.loads(resp)["jobId"]


def wait_for_job(job_id, timeout=600):
    start = time.time()
    while time.time() - start < timeout:
        resp = http_get(UNIPROT_STATUS.format(job_id))
        obj = json.loads(resp)
        if "jobStatus" in obj:
            if obj["jobStatus"] in ("RUNNING", "NEW"):
                time.sleep(3)
                continue
            else:
                raise RuntimeError(f"job failed: {obj}")
        else:
            return
    raise TimeoutError("mapping job did not finish in time")


def get_mapping_results(job_id):
    tsv = http_get(UNIPROT_STREAM.format(job_id))
    lines = tsv.strip().split("\n")
    mapping = {}  # uniprot -> list of ENST (unversioned, deduped, order preserved)
    for line in lines[1:]:
        if not line.strip():
            continue
        frm, to = line.split("\t")
        mapping.setdefault(frm, [])
        if to not in mapping[frm]:
            mapping[frm].append(to)
    return mapping


def fetch_sequences(enst_ids, batch_size=50):
    # Ensembl's sequence endpoint rejects versioned IDs (e.g. ENST00000424496.3);
    # strip the version and map results back onto the bare ID.
    seqs = {}
    enst_ids = list(enst_ids)
    for i in range(0, len(enst_ids), batch_size):
        batch_versioned = enst_ids[i:i+batch_size]
        batch = [e.split(".")[0] for e in batch_versioned]
        payload = json.dumps({"ids": batch}).encode()
        req = urllib.request.Request(
            ENSEMBL_SEQ, data=payload,
            headers={"Content-Type": "application/json"},
        )
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=60) as resp:
                    body = resp.read().decode()
                break
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(5)
                    continue
                body = e.read().decode()
                print(f"  batch error {e.code}: {body[:300]}", file=sys.stderr)
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
        time.sleep(0.3)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="only process first N rows (for testing)")
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
    print(f"Rows to process: {len(rows)}; distinct uniprot IDs: {len(set(uids))}")

    print("Submitting UniProt ID mapping job...")
    job_id = submit_mapping(uids)
    print(f"  job id: {job_id}")
    wait_for_job(job_id)
    print("  job finished, fetching results...")
    mapping = get_mapping_results(job_id)
    print(f"  {len(mapping)} of {len(set(uids))} uniprot IDs mapped to >=1 Ensembl transcript")

    multi = {k: v for k, v in mapping.items() if len(v) > 1}
    print(f"  {len(multi)} uniprot IDs mapped to multiple transcripts (using first returned)")

    chosen_enst = {}  # uid -> enst (unversioned or versioned as returned)
    for uid, ensts in mapping.items():
        chosen_enst[uid] = ensts[0]

    all_ensts = sorted(set(chosen_enst.values()))
    print(f"Fetching {len(all_ensts)} transcript sequences from Ensembl...")
    seqs = fetch_sequences(all_ensts)
    print(f"  got sequences for {len(seqs)} / {len(all_ensts)} transcripts")

    out_fields = fieldnames + ["ensembl_transcript_used", "rna_sequence"]
    n_written = 0
    n_with_rna = 0
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        for row in rows:
            uid = row["uniprot_accession"].strip()
            enst = chosen_enst.get(uid, "")
            seq = seqs.get(enst.split(".")[0], "") if enst else ""
            row["ensembl_transcript_used"] = enst
            row["rna_sequence"] = seq
            writer.writerow(row)
            n_written += 1
            if seq:
                n_with_rna += 1

    print(f"Wrote {n_written} rows to {args.out}; {n_with_rna} have an rna_sequence")


if __name__ == "__main__":
    main()
