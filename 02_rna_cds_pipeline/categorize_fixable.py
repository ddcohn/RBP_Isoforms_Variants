import csv
import sys
import json
import time
import urllib.request
from collections import Counter

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"


def http_get(url, timeout=60):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def main():
    no_xref_rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["transcript_selection_method"] == "no_ensembl_xref":
                no_xref_rows.append(row)

    with_gene = [r for r in no_xref_rows if (r.get("ensembl_gene") or "").strip()]
    print(f"Total no_ensembl_xref: {len(no_xref_rows)}")
    print(f"With ensembl_gene populated: {len(with_gene)}")

    accs = [r["uniprot_accession"] for r in with_gene]
    names = {}
    batch_size = 100
    for i in range(0, len(accs), batch_size):
        batch = accs[i:i+batch_size]
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions="
               + ",".join(batch) + "&fields=accession,protein_name,gene_names,protein_existence")
        for attempt in range(4):
            try:
                body = http_get(url)
                break
            except Exception as e:
                print(f"  error at {i}: {e}", file=sys.stderr)
                time.sleep(5)
                body = None
        if not body:
            continue
        d = json.loads(body)
        for entry in d.get("results", []):
            acc = entry["primaryAccession"]
            pdesc = entry.get("proteinDescription", {})
            name = (pdesc.get("recommendedName", {}).get("fullName", {}).get("value")
                    or (pdesc.get("submissionNames") or [{}])[0].get("fullName", {}).get("value") or "")
            genes = [g.get("geneName", {}).get("value", "") for g in entry.get("genes", [])]
            existence = entry.get("proteinExistence", "")
            names[acc] = (name, ";".join(genes), existence)
        print(f"  ... {i+len(batch)}/{len(accs)} done", flush=True)
        time.sleep(0.3)

    cats = Counter()
    examples = {}
    for r in with_gene:
        acc = r["uniprot_accession"]
        name, genes, existence = names.get(acc, ("", "", ""))
        combined = (name + " " + genes).lower()
        if any(k in combined for k in ["t cell receptor", "trbv","trbj","trbc","trav","traj","trac",
                                         "trgv","trgj","trgc","trdv","trdj","trdc"]):
            cat = "TCR gene segment"
        elif any(k in combined for k in ["immunoglobulin", "igkv","igkj","igkc","iglv","iglj","iglc",
                                           "ighv","ighj","ighd","ighg","ighm","igha"]):
            cat = "Immunoglobulin gene segment"
        elif genes.upper().rstrip(";").endswith("P") and "-as1" not in combined:
            cat = "Pseudogene (symbol ends in P)"
        elif "-as1" in combined or "antisense" in combined or "linc" in combined or "uncharacterized" in combined:
            cat = "Antisense/lncRNA/uncharacterized locus"
        elif "mt-" in genes.lower() or "mitochondri" in combined:
            cat = "Mitochondrial-encoded micropeptide"
        else:
            cat = "UNCLASSIFIED (potentially fixable)"
        cats[cat] += 1
        if len(examples.setdefault(cat, [])) < 8:
            examples[cat].append((acc, genes, name[:50], existence, r.get("ensembl_gene", "")[:40]))

    print()
    print("=== Categories among the 752 with ensembl_gene populated ===")
    for cat, count in cats.most_common():
        print(f"{cat}: {count}")
        for ex in examples[cat]:
            print(f"    {ex}")


if __name__ == "__main__":
    main()
