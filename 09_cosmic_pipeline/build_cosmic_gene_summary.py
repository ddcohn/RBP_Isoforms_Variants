import csv
import sys
import gzip
import time
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

GENES_FILE = "/u/project/kappel/ddcohn/Cosmic_Genes_v104_GRCh38.tsv.gz"
MUTANT_FILE = "/u/project/kappel/ddcohn/Cosmic_GenomeScreensMutant_v104_GRCh38.tsv.gz"
OUT_FILE = "/u/project/kappel/ddcohn/Cosmic_Gene_Summary_v104_GRCh38.csv"

# --- load gene reference table ---
print("Loading Cosmic_Genes...")
genes = {}
with gzip.open(GENES_FILE, mode="rt", newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        cgid = row["COSMIC_GENE_ID"]
        ensg = (row.get("GENE_ACCESSION") or "").strip()
        ensg_bare = ensg.split(".")[0] if ensg else ""
        genes[cgid] = {
            "GENE_SYMBOL": row.get("GENE_SYMBOL", ""),
            "GENE_ACCESSION": ensg,
            "ENSG_bare": ensg_bare,
            "ENTREZ_ID": row.get("ENTREZ_ID", ""),
            "HGNC_ID": row.get("HGNC_ID", ""),
            "IN_CANCER_CENSUS": row.get("IN_CANCER_CENSUS", ""),
            "IS_EXPERT_CURATED": row.get("IS_EXPERT_CURATED", ""),
        }
print(f"  {len(genes):,} genes loaded")

# --- aggregate mutation counts per gene ---
print("Aggregating GenomeScreensMutant...")
total_rows = defaultdict(int)
distinct_mutations = defaultdict(set)
distinct_samples = defaultdict(set)
confirmed_somatic = defaultdict(int)
missense = defaultdict(int)
truncating = defaultdict(int)
synonymous = defaultdict(int)

start = time.time()
n = 0
with gzip.open(MUTANT_FILE, mode="rt", newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        n += 1
        cgid = row.get("COSMIC_GENE_ID", "")
        if not cgid:
            continue
        total_rows[cgid] += 1
        gm = row.get("GENOMIC_MUTATION_ID", "")
        if gm:
            distinct_mutations[cgid].add(gm)
        sid = row.get("COSMIC_SAMPLE_ID", "")
        if sid:
            distinct_samples[cgid].add(sid)
        status = row.get("MUTATION_SOMATIC_STATUS", "")
        if status == "Confirmed somatic variant":
            confirmed_somatic[cgid] += 1
        desc = row.get("MUTATION_DESCRIPTION", "")
        if "missense_variant" in desc:
            missense[cgid] += 1
        if "stop_gained" in desc or "frameshift" in desc:
            truncating[cgid] += 1
        if "synonymous_variant" in desc:
            synonymous[cgid] += 1
        if n % 5_000_000 == 0:
            elapsed = time.time() - start
            print(f"  ... {n:,} rows (elapsed {elapsed/60:.1f}min)", flush=True)

elapsed = time.time() - start
print(f"Done aggregating. Total rows: {n:,} (elapsed {elapsed/60:.1f}min)")

# --- write output: one row per COSMIC gene ---
print(f"Writing {OUT_FILE} ...")
fieldnames = [
    "COSMIC_GENE_ID", "GENE_SYMBOL", "ENSEMBL_GENE_ID", "ENTREZ_ID", "HGNC_ID",
    "IN_CANCER_CENSUS", "IS_EXPERT_CURATED",
    "Total_Mutation_Rows", "Distinct_Genomic_Mutations", "Distinct_Samples",
    "Confirmed_Somatic_Rows", "Missense_Rows", "Truncating_Rows", "Synonymous_Rows",
]

n_written = 0
with open(OUT_FILE, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for cgid, ginfo in genes.items():
        if total_rows.get(cgid, 0) == 0:
            continue  # skip genes with zero mutations in this file
        writer.writerow({
            "COSMIC_GENE_ID": cgid,
            "GENE_SYMBOL": ginfo["GENE_SYMBOL"],
            "ENSEMBL_GENE_ID": ginfo["ENSG_bare"],
            "ENTREZ_ID": ginfo["ENTREZ_ID"],
            "HGNC_ID": ginfo["HGNC_ID"],
            "IN_CANCER_CENSUS": ginfo["IN_CANCER_CENSUS"],
            "IS_EXPERT_CURATED": ginfo["IS_EXPERT_CURATED"],
            "Total_Mutation_Rows": total_rows.get(cgid, 0),
            "Distinct_Genomic_Mutations": len(distinct_mutations.get(cgid, ())),
            "Distinct_Samples": len(distinct_samples.get(cgid, ())),
            "Confirmed_Somatic_Rows": confirmed_somatic.get(cgid, 0),
            "Missense_Rows": missense.get(cgid, 0),
            "Truncating_Rows": truncating.get(cgid, 0),
            "Synonymous_Rows": synonymous.get(cgid, 0),
        })
        n_written += 1

print(f"Wrote {n_written:,} gene rows (of {len(genes):,} total genes in reference file) to {OUT_FILE}")

# genes with mutations but no COSMIC_GENE_ID match in reference file
unmatched = set(total_rows.keys()) - set(genes.keys())
print(f"COSMIC_GENE_IDs seen in mutation file but missing from Cosmic_Genes reference: {len(unmatched)}")
