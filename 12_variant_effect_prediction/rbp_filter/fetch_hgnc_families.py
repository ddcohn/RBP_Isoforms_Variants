import urllib.request
import json
import time
import os
import csv

csv.field_size_limit(2**31 - 1)
D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"
OUT_FILE = f"{D}/hgnc_families.json"

with open(f"{D}/rbp_wt_list_with_go.tsv", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    genes = [row["GeneSymbol"] for row in r]
print(f"Genes to look up: {len(genes)}")

results = {}
if os.path.exists(OUT_FILE):
    results = json.load(open(OUT_FILE))
    print(f"Resuming: {len(results)} already fetched")

remaining = [g for g in genes if g not in results]
print(f"Remaining: {len(remaining)}")


def fetch(symbol):
    url = f"https://rest.genenames.org/fetch/symbol/{symbol}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=20) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:
            time.sleep(1.5)
    return None


n_found = n_no_group = n_no_hit = 0
for i, g in enumerate(remaining):
    data = fetch(g)
    if data is None:
        results[g] = None
        n_no_hit += 1
    else:
        docs = data.get("response", {}).get("docs", [])
        if not docs:
            results[g] = []
            n_no_hit += 1
        else:
            groups = docs[0].get("gene_group", [])
            results[g] = groups
            if groups:
                n_found += 1
            else:
                n_no_group += 1
    if (i + 1) % 100 == 0:
        print(f"  {i+1}/{len(remaining)}: found={n_found} no_group={n_no_group} no_hit={n_no_hit}")
        json.dump(results, open(OUT_FILE, "w"))
    time.sleep(0.2)

json.dump(results, open(OUT_FILE, "w"))
print(f"\nTotal: found_with_family={n_found}, matched_no_family={n_no_group}, no_hgnc_match={n_no_hit}")
print(f"Wrote {OUT_FILE}")
