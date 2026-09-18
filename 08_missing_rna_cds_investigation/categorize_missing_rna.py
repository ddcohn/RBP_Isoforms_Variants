import csv
import sys
import json
import time
import urllib.request
from collections import Counter

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"

rows = []
with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows.append(row)

missing = [r for r in rows if not (r.get("rna_sequence") or "").strip()]
print(f"Rows missing rna_sequence: {len(missing)} / {len(rows)}")

method_counts = Counter((r.get("transcript_selection_method") or "(none)").strip() or "(none)" for r in missing)
print("\nBy transcript_selection_method:")
for m, c in method_counts.most_common():
    print(f"  {m}: {c}")

has_gene = [r for r in missing if (r.get("ensembl_gene") or "").strip()]
print(f"\nOf those missing, {len(has_gene)} DO have an ensembl_gene (searched but no exact match found)")
print(f"and {len(missing) - len(has_gene)} have NO ensembl_gene at all (nothing to search from)")

accs = [r["uniprot_accession"] for r in missing]
has_gene_map = {r["uniprot_accession"]: bool((r.get("ensembl_gene") or "").strip()) for r in missing}

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
    time.sleep(0.3)

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
    elif not has_gene_map.get(acc, False):
        cat = "No ensembl_gene at all"
    else:
        cat = "Other/unclassified (searched all transcripts, no exact match)"
    cats[cat] += 1
    if len(examples.setdefault(cat, [])) < 5:
        examples[cat].append((acc, name, genes))

print()
print("=== Biological category breakdown ===")
for cat, count in cats.most_common():
    print(f"{cat}: {count} ({100*count/len(accs):.1f}%)")
    for acc, name, genes in examples[cat]:
        print(f"    {acc}  gene={genes}  name={name}")

# write full CSV for all 876 missing rows
acc_to_row = {r["uniprot_accession"]: r for r in missing}
acc_to_cat = {}
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
    elif not has_gene_map.get(acc, False):
        cat = "No ensembl_gene at all"
    else:
        cat = "Other/unclassified (searched all transcripts, no exact match)"
    acc_to_cat[acc] = cat

OUT_CSV = "/u/home/d/ddcohn/_claude_missing_rna.csv"
with open(OUT_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["uniprot_accession", "ensembl_gene", "category", "gene_names", "protein_name"])
    for acc in accs:
        row = acc_to_row[acc]
        name, genes = names.get(acc, ("", ""))
        writer.writerow([acc, (row.get("ensembl_gene") or "").strip(), acc_to_cat[acc], genes, name])

print(f"\nWrote {len(accs)} rows to {OUT_CSV}")
