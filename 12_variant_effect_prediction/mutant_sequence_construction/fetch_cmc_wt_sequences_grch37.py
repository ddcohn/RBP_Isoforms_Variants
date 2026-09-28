import urllib.request
import json
import time
import os
import sys

# usage: fetch_cmc_wt_sequences_grch37.py <accessions_txt> <out_json> <out_failed_txt>
ACC_FILE = sys.argv[1]
OUT_FILE = sys.argv[2]
FAILED_FILE = sys.argv[3]

BATCH_SIZE = 50
SLEEP_BETWEEN = 0.15
# use the GRCh37 archive exclusively -- CMC's numbering is GRCh37-era,
# and the live/current endpoint can return a structurally different
# (revised) transcript model for the same stable ID.
URL = "https://grch37.rest.ensembl.org/sequence/id?content-type=application/json;type=protein"


def fetch_batch(ids):
    body = json.dumps({"ids": ids}).encode()
    req = urllib.request.Request(URL, data=body, headers={"Content-Type": "application/json"}, method="POST")
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"  batch attempt {attempt+1} failed: {e}")
            time.sleep(2)
    return None


with open(ACC_FILE) as f:
    accessions = [line.strip() for line in f if line.strip()]

print(f"Total accessions to fetch (GRCh37 archive only): {len(accessions)}")

all_results = {}
if os.path.exists(OUT_FILE):
    with open(OUT_FILE) as f:
        all_results = json.load(f)
    print(f"Resuming: {len(all_results)} already fetched")

remaining = [a for a in accessions if a not in all_results]
print(f"Remaining to fetch: {len(remaining)}")

failed = []

for i in range(0, len(remaining), BATCH_SIZE):
    batch = remaining[i:i + BATCH_SIZE]
    print(f"Fetching batch {i // BATCH_SIZE + 1}/{(len(remaining) + BATCH_SIZE - 1) // BATCH_SIZE} ({len(batch)} accessions)...")
    result = fetch_batch(batch)
    if result is None:
        failed.extend(batch)
        continue
    got = set()
    for entry in result:
        if isinstance(entry, dict) and "query" in entry and "seq" in entry:
            all_results[entry["query"]] = {"protein_id": entry.get("id"), "sequence": entry["seq"]}
            got.add(entry["query"])
    for acc in batch:
        if acc not in got:
            failed.append(acc)
    print(f"  got {len(got)}/{len(batch)}")

    with open(OUT_FILE, "w") as f:
        json.dump(all_results, f)

    time.sleep(SLEEP_BETWEEN)

print(f"\nTotal fetched: {len(all_results)}")
print(f"Total failed: {len(failed)}")

with open(FAILED_FILE, "w") as f:
    for acc in failed:
        f.write(acc + "\n")

print(f"Wrote {OUT_FILE}")
print(f"Wrote {FAILED_FILE}")
