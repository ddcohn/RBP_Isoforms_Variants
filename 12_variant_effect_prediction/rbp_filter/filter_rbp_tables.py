import csv, sys, gzip

csv.field_size_limit(sys.maxsize)
BASE = "/u/project/kappel/ddcohn/protein_variant_effects"

with open(f"{BASE}/rbp_gene_symbols.txt") as f:
    RBP = {line.strip() for line in f if line.strip()}
print(f"RBP gene symbols: {len(RBP)}")


def filter_by_column(src, dst, gene_col, delim="\t"):
    n_in = n_out = 0
    with open(src, newline="") as f, open(dst, "w", newline="") as out:
        r = csv.DictReader(f, delimiter=delim)
        w = csv.writer(out, delimiter=delim)
        w.writerow(r.fieldnames)
        for row in r:
            n_in += 1
            if row[gene_col] in RBP:
                w.writerow([row[c] for c in r.fieldnames])
                n_out += 1
    print(f"{src} -> {dst}: {n_in} -> {n_out} rows")


# DeepLoc / protGPS deltas already carry a gene column directly
filter_by_column(f"{BASE}/clinvar_deeploc_deltas.tsv", f"{BASE}/clinvar_deeploc_deltas_rbp.tsv", "GeneSymbol")
filter_by_column(f"{BASE}/clinvar_protgps_deltas.tsv", f"{BASE}/clinvar_protgps_deltas_rbp.tsv", "GeneSymbol")
filter_by_column(f"{BASE}/cmc_deeploc_deltas.tsv", f"{BASE}/cmc_deeploc_deltas_rbp.tsv", "GENE_NAME")
filter_by_column(f"{BASE}/cmc_protgps_deltas.tsv", f"{BASE}/cmc_protgps_deltas_rbp.tsv", "GENE_NAME")

# ClinVar metadata table itself is variant -> gene, so filter it too and use
# it to filter clinvar_spliceai_scores.tsv (VariationID has no gene column
# of its own -- SpliceAI needed no protein/gene info to run).
filter_by_column(f"{BASE}/clinvar_metadata_lookup.tsv", f"{BASE}/clinvar_metadata_lookup_rbp.tsv", "GeneSymbol")

clinvar_rbp_vids = set()
with open(f"{BASE}/clinvar_metadata_lookup.tsv", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row["GeneSymbol"] in RBP:
            clinvar_rbp_vids.add(row["VariationID"])
print(f"ClinVar VariationIDs belonging to an RBP gene: {len(clinvar_rbp_vids)}")

n_in = n_out = 0
with open(f"{BASE}/clinvar_spliceai_scores.tsv", newline="") as f, \
     open(f"{BASE}/clinvar_spliceai_scores_rbp.tsv", "w", newline="") as out:
    r = csv.DictReader(f, delimiter="\t")
    w = csv.writer(out, delimiter="\t")
    w.writerow(r.fieldnames)
    for row in r:
        n_in += 1
        if row["VariationID"] in clinvar_rbp_vids:
            w.writerow([row[c] for c in r.fieldnames])
            n_out += 1
print(f"clinvar_spliceai_scores.tsv -> clinvar_spliceai_scores_rbp.tsv: {n_in} -> {n_out} rows")

# CMC has no per-variant gene column in the SpliceAI output either, and
# cmc_protein_changes.tsv only covers the protein-editable categories (not
# synonymous/frameshift/stoploss/noncoding, which SpliceAI also scored) --
# so build the GENOMIC_MUTATION_ID -> GENE_NAME map fresh from the raw file,
# which has a gene name on every row regardless of category.
cmc_id_to_gene = {}
with gzip.open("/u/project/kappel/ddcohn/CancerMutationCensus_AllData_v104_GRCh37.tsv.gz",
               "rt", encoding="utf-8", errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        mid = row.get("GENOMIC_MUTATION_ID", "")
        gene = row.get("GENE_NAME", "")
        if mid and gene and mid not in cmc_id_to_gene:
            cmc_id_to_gene[mid] = gene
print(f"CMC GENOMIC_MUTATION_ID -> GENE_NAME pairs: {len(cmc_id_to_gene)}")

n_in = n_out = 0
with open(f"{BASE}/cmc_spliceai_scores.tsv", newline="") as f, \
     open(f"{BASE}/cmc_spliceai_scores_rbp.tsv", "w", newline="") as out:
    r = csv.DictReader(f, delimiter="\t")
    w = csv.writer(out, delimiter="\t")
    w.writerow(r.fieldnames)
    for row in r:
        n_in += 1
        gene = cmc_id_to_gene.get(row["GENOMIC_MUTATION_ID"])
        if gene in RBP:
            w.writerow([row[c] for c in r.fieldnames])
            n_out += 1
print(f"cmc_spliceai_scores.tsv -> cmc_spliceai_scores_rbp.tsv: {n_in} -> {n_out} rows")
