BASE = "/u/project/kappel/ddcohn/protein_variant_effects"

with open(f"{BASE}/rbp_gene_symbols.txt") as f:
    domain_based = {l.strip() for l in f if l.strip()}
with open(f"{BASE}/rbp_gene_symbols_rbp2go.txt") as f:
    rbp2go = {l.strip() for l in f if l.strip()}

union = domain_based | rbp2go
print(f"domain-based: {len(domain_based)}")
print(f"rbp2go: {len(rbp2go)}")
print(f"union: {len(union)}")

with open(f"{BASE}/rbp_gene_symbols_domain_based.txt", "w") as f:
    for s in sorted(domain_based):
        f.write(s + "\n")

with open(f"{BASE}/rbp_gene_symbols.txt", "w") as f:
    for s in sorted(union):
        f.write(s + "\n")
print("wrote rbp_gene_symbols_domain_based.txt (533, kept for reference) and overwrote rbp_gene_symbols.txt (union, now the active filter list)")
