import csv, sys, gzip, json

csv.field_size_limit(sys.maxsize)
BASE = "/u/project/kappel/ddcohn/protein_variant_effects"

# 1. RBP NCBI gene IDs (Has_RBD == 1) from the user's own table
rbp_geneids = set()
with open("/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv", newline="", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f)
    for row in r:
        if row.get("Has_RBD", "").strip() == "1":
            for gid in row.get("ncbi_geneid", "").split(";"):
                gid = gid.strip()
                if gid:
                    rbp_geneids.add(gid)
print(f"RBP NCBI gene IDs (Has_RBD==1): {len(rbp_geneids)}")

# 2. GeneID -> GeneSymbol map from the raw ClinVar file (covers essentially the whole genome)
geneid_to_symbol = {}
with open("/u/project/kappel/ddcohn/ClinVar_variant_summary_complete.csv", newline="", encoding="utf-8") as f:
    r = csv.DictReader(f)
    for row in r:
        gid = row["GeneID"].strip()
        sym = row["GeneSymbol"].strip()
        if gid and sym and gid not in geneid_to_symbol:
            geneid_to_symbol[gid] = sym
print(f"GeneID->GeneSymbol pairs available from ClinVar: {len(geneid_to_symbol)}")

rbp_symbols = set()
unresolved = []
for gid in rbp_geneids:
    sym = geneid_to_symbol.get(gid)
    if sym:
        rbp_symbols.add(sym)
    else:
        unresolved.append(gid)
print(f"RBP gene IDs resolved to a symbol via ClinVar: {len(rbp_symbols)}")
print(f"RBP gene IDs NOT found in ClinVar (need another source): {len(unresolved)}")
if unresolved:
    print("first 20 unresolved geneids:", unresolved[:20])

with open(f"{BASE}/rbp_gene_symbols.txt", "w") as f:
    for s in sorted(rbp_symbols):
        f.write(s + "\n")
print(f"Wrote {BASE}/rbp_gene_symbols.txt ({len(rbp_symbols)} symbols)")
