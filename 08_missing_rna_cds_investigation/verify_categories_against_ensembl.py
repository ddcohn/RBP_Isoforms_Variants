"""
Verifies the name-keyword-based categorization in categorize_missing_rna.py against
Ensembl's authoritative gene biotype, for every one of the 876 missing-RNA/CDS rows that
has an ensembl_gene. Run locally against RBP_missing_rna_876.csv (the output of
categorize_missing_rna.py).

Finding: the original keyword heuristic was directionally right (most "Other/unclassified"
rows really are pseudogenes/lncRNAs) but individually unreliable -- ~42% of rows changed
category once checked against real Ensembl biotype, including a small number of genuine
mislabels (e.g. FCGR1BP tagged "Immunoglobulin gene segment" by keyword match, but its real
biotype is a plain pseudogene unrelated to Ig).
"""
import csv
import json
import time
import urllib.request
from collections import defaultdict

PATH = "RBP_missing_rna_876.csv"  # output of categorize_missing_rna.py

rows = []
with open(PATH, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        rows.append(row)

with_ensg = [r for r in rows if r["ensembl_gene"].strip()]
no_ensg = [r for r in rows if not r["ensembl_gene"].strip()]
print(f"Total rows: {len(rows)}")
print(f"Rows with ensembl_gene: {len(with_ensg)}")
print(f"Rows with NO ensembl_gene: {len(no_ensg)}")

# batch-query Ensembl for real biotypes
biotypes = {}
ensgs = [r["ensembl_gene"].strip() for r in with_ensg]
BATCH = 100
for i in range(0, len(ensgs), BATCH):
    batch = ensgs[i:i + BATCH]
    url = "https://rest.ensembl.org/lookup/id?content-type=application/json"
    payload = json.dumps({"ids": batch}).encode()
    req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json"}, method="POST")
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                d = json.loads(resp.read().decode())
            break
        except Exception as e:
            print(f"  error batch {i}: {e}", flush=True)
            time.sleep(3)
            d = {}
    for ensg, info in d.items():
        biotypes[ensg] = info.get("biotype", "") if info else "NOT_FOUND"
    print(f"  ... {min(i + BATCH, len(ensgs))}/{len(ensgs)} genes looked up", flush=True)
    time.sleep(0.2)


def true_category(biotype):
    if biotype in ("IG_V_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "IG_D_pseudogene",
                   "IG_V_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene", "IG_pseudogene"):
        return "Immunoglobulin gene segment"
    if biotype in ("TR_V_pseudogene", "TR_J_pseudogene", "TR_D_pseudogene", "TR_C_pseudogene",
                   "TR_V_gene", "TR_J_gene", "TR_D_gene", "TR_C_gene"):
        return "TCR gene segment"
    if biotype in ("Mt_rRNA", "Mt_tRNA"):
        return "Mitochondrial-encoded micropeptide"
    if biotype.endswith("pseudogene"):
        return "Pseudogene"
    if biotype in ("lncRNA", "antisense", "lincRNA", "macro_lncRNA", "sense_intronic",
                   "sense_overlapping", "bidirectional_promoter_lncRNA", "3prime_overlapping_ncRNA",
                   "non_coding", "processed_transcript", "retained_intron", "TEC"):
        return "Antisense/lncRNA/uncharacterized locus"
    if biotype == "protein_coding":
        return "Protein-coding (genuine no-exact-match case)"
    if biotype == "NOT_FOUND":
        return "Ensembl ID not found / retired"
    return f"Other biotype: {biotype}"


def normalize(cat):
    if cat.startswith("Pseudogene"):
        return "Pseudogene"
    if cat.startswith("Other/unclassified"):
        return "Other/unclassified"
    return cat


agree = disagree = 0
resolve_to = defaultdict(int)
real_errors = []
true_dist = defaultdict(int)

for r in with_ensg:
    ensg = r["ensembl_gene"].strip()
    bt = biotypes.get(ensg, "LOOKUP_FAILED")
    tc = true_category(bt)
    true_dist[tc] += 1
    orig_norm = normalize(r["category"])
    tc_norm = normalize(tc)
    if orig_norm == tc_norm:
        agree += 1
    else:
        disagree += 1
        resolve_to[tc_norm] += 1
        if orig_norm != "Other/unclassified":
            real_errors.append((r["uniprot_accession"], r["gene_names"], r["category"], bt, tc))

print(f"\nAgree (normalized): {agree}/{len(with_ensg)} ({100 * agree / len(with_ensg):.1f}%)")
print(f"Disagree: {disagree}/{len(with_ensg)} ({100 * disagree / len(with_ensg):.1f}%)")

print("\nMismatches resolve to:")
for k, v in sorted(resolve_to.items(), key=lambda x: -x[1]):
    print(f"  {k}: {v}")

print(f"\nGenuine heuristic ERRORS (specific wrong label, not just vague 'Other'): {len(real_errors)}")
for x in real_errors:
    print(" ", x)

print("\nTRUE distribution (of rows with ensembl_gene, verified against real Ensembl biotype):")
for cat, count in sorted(true_dist.items(), key=lambda x: -x[1]):
    print(f"  {cat}: {count} ({100 * count / len(with_ensg):.1f}%)")

# how many of the no-ensembl_gene rows were themselves mislabeled by the name heuristic
# instead of being flagged "No ensembl_gene at all"?
from collections import Counter
c = Counter(r["category"] for r in no_ensg)
print(f"\nOf the {len(no_ensg)} rows with NO ensembl_gene, how they were actually labeled:")
for k, v in c.most_common():
    print(f"  labeled '{k}': {v}")
