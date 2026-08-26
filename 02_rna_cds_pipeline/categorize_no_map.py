import json
import time
import urllib.request
import urllib.error
from collections import Counter

with open("/tmp/no_map_ids.txt") as f:
    accs = [l.strip() for l in f if l.strip()]

print(f"Total to categorize: {len(accs)}")

names = {}
batch_size = 100
for i in range(0, len(accs), batch_size):
    batch = accs[i:i+batch_size]
    url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions="
           + ",".join(batch) + "&fields=accession,protein_name,gene_names")
    for attempt in range(5):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=60) as resp:
                d = json.load(resp)
            break
        except Exception as e:
            print(f"  error at {i}: {e}", flush=True)
            time.sleep(5)
            d = {"results": []}
    for e in d.get("results", []):
        acc = e["primaryAccession"]
        pdesc = e.get("proteinDescription", {})
        name = (pdesc.get("recommendedName", {}).get("fullName", {}).get("value")
                or (pdesc.get("submissionNames") or [{}])[0].get("fullName", {}).get("value")
                or "")
        genes = [g.get("geneName", {}).get("value", "") for g in e.get("genes", [])]
        names[acc] = (name, ";".join(genes))
    print(f"  ... {i+len(batch)}/{len(accs)} done", flush=True)
    time.sleep(0.3)

# categorize by keyword heuristics
cats = Counter()
examples = {}
for acc in accs:
    name, genes = names.get(acc, ("", ""))
    combined = (name + " " + genes).lower()
    if any(k in combined for k in ["t cell receptor", "trbv", "trbj", "trbc", "trav", "traj", "trac",
                                     "trgv", "trgj", "trgc", "trdv", "trdj", "trdc"]):
        cat = "TCR gene segment"
    elif any(k in combined for k in ["immunoglobulin", "ig heavy", "ig kappa", "ig lambda",
                                       "igkv", "igkj", "igkc", "iglv", "iglj", "iglc", "ighv", "ighj", "ighd", "ighg", "ighm", "igha"]):
        cat = "Immunoglobulin gene segment"
    elif combined.strip().endswith("p") and "-as1" not in combined and any(g.endswith("P") for g in genes.split(";") if g):
        cat = "Pseudogene (gene symbol ends in P)"
    elif "-as1" in combined or "antisense" in combined or "linc" in combined or "uncharacterized" in combined:
        cat = "Antisense/lncRNA/uncharacterized locus"
    elif "mt-" in genes.lower() or "mitochondri" in combined:
        cat = "Mitochondrial-encoded micropeptide"
    else:
        cat = "Other/unclassified"
    cats[cat] += 1
    if len(examples.setdefault(cat, [])) < 5:
        examples[cat].append((acc, name, genes))

print()
print("=== Category breakdown ===")
for cat, count in cats.most_common():
    print(f"{cat}: {count}")
    for acc, name, genes in examples[cat]:
        print(f"    {acc}  gene={genes}  name={name}")
