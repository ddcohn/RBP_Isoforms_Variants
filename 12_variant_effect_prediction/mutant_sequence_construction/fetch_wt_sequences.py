import urllib.request
import urllib.parse
import time
import re
import json
import os

ACC_FILE = "/u/project/kappel/ddcohn/protein_variant_effects/unique_transcript_accessions.txt"
OUT_FILE = "/u/project/kappel/ddcohn/protein_variant_effects/wt_sequences.json"
FAILED_FILE = "/u/project/kappel/ddcohn/protein_variant_effects/wt_sequences_failed.txt"

BATCH_SIZE = 200
SLEEP_BETWEEN = 0.4  # ~2.5 req/sec, polite without an API key
EMAIL = "daniel.cohn@hhmi.org"
TOOL = "rbp_variant_effect_pipeline"

RECORD_SPLIT_RE = re.compile(r"\n//\n")
ACCVER_RE = re.compile(r"^VERSION\s+(\S+)", re.MULTILINE)
CDS_BLOCK_RE = re.compile(r"^     CDS .*?(?=^     \S|\Z)", re.MULTILINE | re.DOTALL)
PROTEIN_ID_RE = re.compile(r'/protein_id="([^"]+)"')
TRANSLATION_RE = re.compile(r'/translation="([^"]+)"', re.DOTALL)


def fetch_batch(accessions):
    ids = ",".join(accessions)
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": ids,
        "rettype": "gb",
        "retmode": "text",
        "email": EMAIL,
        "tool": TOOL,
    }
    full_url = url + "?" + urllib.parse.urlencode(params)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(full_url, timeout=60) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:
            print(f"  fetch attempt {attempt+1} failed: {e}")
            time.sleep(2)
    return None


def parse_records(text):
    results = {}
    records = RECORD_SPLIT_RE.split(text)
    for rec in records:
        rec = rec.strip()
        if not rec:
            continue
        m = ACCVER_RE.search(rec)
        if not m:
            continue
        accver = m.group(1)
        cds_m = CDS_BLOCK_RE.search(rec)
        if not cds_m:
            continue
        cds_block = cds_m.group(0)
        pid_m = PROTEIN_ID_RE.search(cds_block)
        trans_m = TRANSLATION_RE.search(cds_block)
        if not trans_m:
            continue
        seq = trans_m.group(1)
        seq = re.sub(r"\s+", "", seq)
        protein_id = pid_m.group(1) if pid_m else None
        results[accver] = {"protein_id": protein_id, "sequence": seq}
    return results


with open(ACC_FILE) as f:
    accessions = [line.strip() for line in f if line.strip()]

print(f"Total accessions to fetch: {len(accessions)}")

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
    text = fetch_batch(batch)
    if text is None:
        print("  batch failed entirely, recording as failed")
        failed.extend(batch)
        continue
    parsed = parse_records(text)
    for acc in batch:
        if acc not in parsed:
            failed.append(acc)
    all_results.update(parsed)
    print(f"  got {len(parsed)}/{len(batch)} sequences")

    # checkpoint every batch
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
