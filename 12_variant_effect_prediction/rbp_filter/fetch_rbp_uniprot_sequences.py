import urllib.request
import urllib.parse
import time
import json
import os

D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"
ACC_FILE = f"{D}/rbp_uniprot_ids_to_fetch.txt"
OUT_FILE = f"{D}/rbp_wt_sequences_uniprot_v2.json"

BATCH_SIZE = 100
SLEEP_BETWEEN = 0.3


def fetch_batch(accessions):
    query = " OR ".join(f"accession:{a}" for a in accessions)
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {"query": query, "format": "fasta"}
    )
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:
            print(f"  fetch attempt {attempt+1} failed: {e}")
            time.sleep(2)
    return None


def parse_fasta(text):
    # >sp|P12345|GENE_HUMAN Description OS=... GN=GeneName PE=... SV=...
    results = {}
    acc = None
    gene = None
    seq = []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc is not None:
                results[acc] = {"gene": gene, "sequence": "".join(seq)}
            parts = line[1:].split("|")
            acc = parts[1] if len(parts) >= 2 else line[1:].split()[0]
            gene = None
            for tok in line.split():
                if tok.startswith("GN="):
                    gene = tok[3:]
            seq = []
        else:
            seq.append(line.strip())
    if acc is not None:
        results[acc] = {"gene": gene, "sequence": "".join(seq)}
    return results


with open(ACC_FILE) as f:
    accessions = [line.strip() for line in f if line.strip()]
print(f"Total accessions to fetch: {len(accessions)}")

all_results = {}
if os.path.exists(OUT_FILE):
    all_results = json.load(open(OUT_FILE))
    print(f"Resuming: {len(all_results)} already fetched")

remaining = [a for a in accessions if a not in all_results]
print(f"Remaining to fetch: {len(remaining)}")

no_gn = 0
for i in range(0, len(remaining), BATCH_SIZE):
    batch = remaining[i:i + BATCH_SIZE]
    print(f"Fetching batch {i // BATCH_SIZE + 1}/{(len(remaining) + BATCH_SIZE - 1) // BATCH_SIZE} ({len(batch)} accessions)...")
    text = fetch_batch(batch)
    if text is None:
        print("  batch failed entirely")
        continue
    parsed = parse_fasta(text)
    got = set(parsed.keys()) & set(batch)
    for acc in got:
        all_results[acc] = parsed[acc]
        if parsed[acc]["gene"] is None:
            no_gn += 1
    print(f"  got {len(got)}/{len(batch)} sequences")
    json.dump(all_results, open(OUT_FILE, "w"))
    time.sleep(SLEEP_BETWEEN)

print(f"\nTotal fetched: {len(all_results)}")
print(f"Records missing a GN= gene name: {no_gn}")
print(f"Wrote {OUT_FILE}")
