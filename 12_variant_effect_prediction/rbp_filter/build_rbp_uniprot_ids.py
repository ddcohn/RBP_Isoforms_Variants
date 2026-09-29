import csv
import sys

csv.field_size_limit(sys.maxsize)
D = "/u/project/kappel/ddcohn/protein_variant_effects"

with open(f"{D}/rbp_only/rbp_gene_symbols.txt") as f:
    RBP = {l.strip() for l in f if l.strip()}
with open(f"{D}/rbp_only/rbp_gene_symbols_rbp2go.txt") as f:
    rbp2go_symbols = {l.strip() for l in f if l.strip()}
domain_only_symbols = RBP - rbp2go_symbols
print(f"RBP total: {len(RBP)}, from RBP2GO: {len(rbp2go_symbols)}, domain-only (not in RBP2GO): {len(domain_only_symbols)}")

# Source 1: RBP2GO's own UniProt IDs -- direct, no lab table.
rbp2go_uniprot_ids = set()
with open(f"{D}/rbp_only/High_confidence_human_RBPs_rbp2go.txt", encoding="utf-8") as f:
    lines = f.readlines()
hidx = next(i for i, l in enumerate(lines) if l.startswith('"Uniprot_ID"'))
r = csv.DictReader(lines[hidx:], delimiter="\t", quotechar='"')
for row in r:
    if row["RBP_status"] == "RBP":
        rbp2go_uniprot_ids.add(row["Uniprot_ID"])
print(f"UniProt IDs directly from RBP2GO: {len(rbp2go_uniprot_ids)}")

# Source 2: for the domain-only genes, pull uniprot_accession from the
# USER'S OWN table (not the labmate's), matched by the same NCBI GeneIDs
# originally used to build the domain-based list.
symbol_to_geneid = {}
with open("/u/project/kappel/ddcohn/ClinVar_variant_summary_complete.csv", newline="", encoding="utf-8") as f:
    r = csv.DictReader(f)
    for row in r:
        gid = row["GeneID"].strip()
        sym = row["GeneSymbol"].strip()
        if gid and sym and sym not in symbol_to_geneid:
            symbol_to_geneid[sym] = gid

domain_geneids = {symbol_to_geneid[s]: s for s in domain_only_symbols if s in symbol_to_geneid}
print(f"Domain-only symbols resolved to a GeneID: {len(domain_geneids)} of {len(domain_only_symbols)}")

domain_uniprot = {}  # symbol -> uniprot_accession
n_rows = 0
with open("/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv", newline="", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f)
    for row in r:
        n_rows += 1
        for gid in row.get("ncbi_geneid", "").split(";"):
            gid = gid.strip()
            if gid in domain_geneids:
                acc = row.get("uniprot_accession", "").strip()
                if acc:
                    domain_uniprot[domain_geneids[gid]] = acc
print(f"Domain-only symbols with a uniprot_accession from your own table: {len(domain_uniprot)}")

all_ids = rbp2go_uniprot_ids | set(domain_uniprot.values())
print(f"Total unique UniProt IDs to fetch: {len(all_ids)}")

with open(f"{D}/rbp_only/rbp_uniprot_ids_to_fetch.txt", "w") as f:
    for acc in sorted(all_ids):
        f.write(acc + "\n")
print(f"Wrote {D}/rbp_only/rbp_uniprot_ids_to_fetch.txt")
