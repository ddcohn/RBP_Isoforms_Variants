import json

D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"

data = json.load(open(f"{D}/rbp_wt_sequences_uniprot_v2.json"))
print(f"UniProt records: {len(data)}")

no_gene = [acc for acc, rec in data.items() if not rec["gene"]]
print(f"Records with no GN= at all (excluded): {len(no_gene)} -> {no_gene}")

by_gene = {}
for acc, rec in data.items():
    gene = rec["gene"]
    seq = rec["sequence"]
    if not gene:
        continue
    if gene not in by_gene or len(seq) > len(by_gene[gene][1]):
        by_gene[gene] = (acc, seq)

print(f"Clean, unique gene symbols (from UniProt GN=, no lab table involved): {len(by_gene)}")

old_rbp = {l.strip() for l in open(f"{D}/rbp_gene_symbols.txt") if l.strip()}
old_corrupted = {g for g in old_rbp if "{ECO:" in g}
print(f"\nOld list: {len(old_rbp)} entries, {len(old_corrupted)} corrupted (had evidence-code suffix)")
new_rbp = set(by_gene.keys())
print(f"New clean list: {len(new_rbp)}")
print(f"In old (clean names only, no {{ECO:}}) but not recovered in new: {len((old_rbp - old_corrupted) - new_rbp)}")
print(f"Newly gained vs old clean names: {len(new_rbp - (old_rbp - old_corrupted))}")

total_len = sum(len(seq) for _, seq in by_gene.values())
print(f"\nTotal residues: {total_len:,}")
print(f"Missense-only saturation mutants (19x positions): {total_len*19:,}")

with open(f"{D}/rbp_gene_symbols_clean.txt", "w") as f:
    for g in sorted(new_rbp):
        f.write(g + "\n")
with open(f"{D}/rbp_wt_sequences.fasta", "w") as f:
    for gene, (acc, seq) in sorted(by_gene.items()):
        f.write(f">{gene}|{acc}\n{seq}\n")
print(f"\nWrote {D}/rbp_gene_symbols_clean.txt and rbp_wt_sequences.fasta")
