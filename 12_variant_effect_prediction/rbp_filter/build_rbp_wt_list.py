D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"

records = []
gene = acc = None
seq = []
with open(f"{D}/rbp_wt_sequences.fasta") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if gene is not None:
                records.append((gene, acc, "".join(seq)))
            gene, acc = line[1:].split("|", 1)
            seq = []
        else:
            seq.append(line)
    if gene is not None:
        records.append((gene, acc, "".join(seq)))

records.sort(key=lambda r: r[0])
with open(f"{D}/rbp_wt_list.tsv", "w") as f:
    f.write("GeneSymbol\tUniProtID\tLength\tSequence\n")
    for gene, acc, seq in records:
        f.write(f"{gene}\t{acc}\t{len(seq)}\t{seq}\n")

print(f"Wrote {D}/rbp_wt_list.tsv with {len(records)} rows")
