import json

D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"

new_entries = json.load(open(f"{D}/combine_new_entries.json"))
print(f"New entries to add: {len(new_entries)}")

# existing gene -> (accession, sequence) from the current WT fasta
existing = {}
gene = acc = None
seq = []
with open(f"{D}/rbp_wt_sequences.fasta") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if gene is not None:
                existing[gene] = (acc, "".join(seq))
            gene, acc = line[1:].split("|", 1)
            seq = []
        else:
            seq.append(line)
    if gene is not None:
        existing[gene] = (acc, "".join(seq))
print(f"Existing entries: {len(existing)}")

collisions = set(new_entries) & set(existing)
print(f"Collisions (new gene name already present -- shouldn't happen): {len(collisions)}")
if collisions:
    print("  ", collisions)

combined = dict(existing)
for gene, rec in new_entries.items():
    if gene not in combined:
        combined[gene] = (rec["accession"], rec["sequence"])

print(f"Combined total: {len(combined)}")
total_len = sum(len(seq) for _, seq in combined.values())
print(f"Total residues: {total_len:,}")
print(f"Missense-only saturation mutants (19x positions): {total_len*19:,}")

with open(f"{D}/rbp_gene_symbols.txt", "w") as f:
    for g in sorted(combined):
        f.write(g + "\n")
with open(f"{D}/rbp_wt_sequences.fasta", "w") as f:
    for g, (acc, seq) in sorted(combined.items()):
        f.write(f">{g}|{acc}\n{seq}\n")
with open(f"{D}/rbp_wt_list.tsv", "w") as f:
    f.write("GeneSymbol\tUniProtID\tLength\tSequence\n")
    for g, (acc, seq) in sorted(combined.items()):
        f.write(f"{g}\t{acc}\t{len(seq)}\t{seq}\n")

# also write just the 246 new sequences separately, for the incremental
# DeepLoc/protGPS run (no need to rerun the 2089 already done)
new_only = {g: existing_or_new for g, existing_or_new in combined.items() if g not in existing}
with open(f"{D}/rbp_wt_sequences_new_only.fasta", "w") as f:
    for g, (acc, seq) in sorted(new_only.items()):
        f.write(f">{g}|{acc}\n{seq}\n")
print(f"\nWrote combined rbp_gene_symbols.txt, rbp_wt_sequences.fasta, rbp_wt_list.tsv ({len(combined)} genes)")
print(f"Wrote rbp_wt_sequences_new_only.fasta ({len(new_only)} genes, for incremental prediction run)")
