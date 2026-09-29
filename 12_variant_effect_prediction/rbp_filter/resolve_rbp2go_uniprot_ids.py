import csv, sys

csv.field_size_limit(sys.maxsize)
BASE = "/u/project/kappel/ddcohn/protein_variant_effects"

with open(f"{BASE}/rbp2go_uniprot_ids.txt") as f:
    rbp2go_ids = {line.strip() for line in f if line.strip()}
print(f"RBP2GO high-confidence human RBP UniProt IDs: {len(rbp2go_ids)}")

# build uniprot accession (primary + secondary) -> gene_symbol from the
# lab's isoform table (68,645 isoforms, has gene_symbol directly)
acc_to_gene = {}
n_rows = 0
with open("/u/project/kappel/RBP/Isoform_Table/Isoform_Post_Merge_PSLab_OpenTargets_Updated_Interim_20260817.csv",
          newline="", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f)
    for row in r:
        n_rows += 1
        gene = row.get("gene_symbol", "").strip()
        if not gene:
            continue
        uid = row.get("uniprot_id", "").strip()
        if uid:
            acc_to_gene.setdefault(uid, gene)
        for sec in row.get("uniprot_secondary_accessions", "").split(";"):
            sec = sec.strip()
            if sec:
                acc_to_gene.setdefault(sec, gene)
        for can in row.get("swissprot_canonical_accessions", "").split(";"):
            can = can.strip()
            if can:
                acc_to_gene.setdefault(can, gene)
print(f"isoform table rows scanned: {n_rows}, unique accession->gene entries: {len(acc_to_gene)}")

resolved = {}
unresolved = []
for uid in rbp2go_ids:
    gene = acc_to_gene.get(uid)
    if gene:
        resolved[uid] = gene
    else:
        unresolved.append(uid)

symbols = set(resolved.values())
print(f"RBP2GO UniProt IDs resolved to a gene symbol: {len(resolved)} ({len(symbols)} unique gene symbols)")
print(f"RBP2GO UniProt IDs NOT resolved: {len(unresolved)}")
if unresolved:
    print("first 20 unresolved:", unresolved[:20])

with open(f"{BASE}/rbp_gene_symbols_rbp2go.txt", "w") as f:
    for s in sorted(symbols):
        f.write(s + "\n")
print(f"Wrote {BASE}/rbp_gene_symbols_rbp2go.txt ({len(symbols)} symbols)")

# compare against the existing domain-based (Has_RBD==1) list
with open(f"{BASE}/rbp_gene_symbols.txt") as f:
    domain_based = {line.strip() for line in f if line.strip()}
print(f"\nExisting domain-based RBP list: {len(domain_based)} genes")
print(f"Overlap: {len(symbols & domain_based)}")
print(f"In domain-based list but NOT in RBP2GO: {len(domain_based - symbols)}")
print(f"In RBP2GO but NOT in domain-based list: {len(symbols - domain_based)}")
