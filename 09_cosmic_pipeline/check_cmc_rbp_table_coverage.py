"""
Checks the Cancer Mutation Census (CMC) file's coverage against this project's isoform
table, and how reliably its extra columns (ClinVar cross-reference, gnomAD/ExAC frequency,
GERP conservation, SIFT prediction) are actually populated.

Matches by exact Ensembl transcript ID (CMC's ACCESSION_NUMBER vs. the isoform table's
ensembl_transcript_used, both version-stripped) -- a precise but conservative join that
undercounts real gene-level coverage whenever COSMIC and this project's exact-match
pipeline settled on different transcripts of the same gene. See the project README for
why gene-level (ENSG) matching would recover more of the ~31% that doesn't match here.

Requires: rbp_transcripts.csv (uniprot_accession, ensembl_transcript_bare), extracted from
the isoform table's ensembl_transcript_used column with the version suffix stripped, and
the CMC file itself (large, not included in this repo -- see 10_clinvar_repull or the
COSMIC scripted-download flow to obtain it).
"""
import csv
import sys
import gzip

csv.field_size_limit(sys.maxsize)

CMC = "CancerMutationCensus_AllData_v104_GRCh37.tsv.gz"
RBP_TRANSCRIPTS = "rbp_transcripts.csv"

rbp_transcripts = set()
with open(RBP_TRANSCRIPTS, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        t = row["ensembl_transcript_bare"].strip()
        if t:
            rbp_transcripts.add(t)

print(f"Distinct RBP transcripts to match against: {len(rbp_transcripts)}")

n = 0
n_matched_rows = 0
matched_transcripts = set()
matched_genes = set()
all_genes = set()

pop_cols = ["CLINVAR_CLNSIG", "GNOMAD_EXOMES_AF", "GNOMAD_GENOMES_AF", "GERP++_RS",
            "MIN_SIFT_SCORE", "DNDS_DISEASE_QVAL_SIG"]
pop_counts_all = {c: 0 for c in pop_cols}
pop_counts_matched = {c: 0 for c in pop_cols}

with gzip.open(CMC, mode="rt", newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        n += 1
        gene = row.get("GENE_NAME", "")
        all_genes.add(gene)
        acc = (row.get("ACCESSION_NUMBER") or "").strip().split(".")[0]
        is_match = acc in rbp_transcripts
        for c in pop_cols:
            if (row.get(c) or "").strip():
                pop_counts_all[c] += 1
                if is_match:
                    pop_counts_matched[c] += 1
        if is_match:
            n_matched_rows += 1
            matched_transcripts.add(acc)
            matched_genes.add(gene)
        if n % 1_000_000 == 0:
            print(f"  ... {n:,} rows processed", flush=True)

print(f"\nTotal CMC rows: {n:,}")
print(f"Distinct genes in CMC: {len(all_genes):,}")
print(f"Rows matching an RBP transcript: {n_matched_rows:,} ({100 * n_matched_rows / n:.1f}%)")
print(f"Distinct RBP transcripts found in CMC: {len(matched_transcripts):,} / "
      f"{len(rbp_transcripts):,} ({100 * len(matched_transcripts) / len(rbp_transcripts):.1f}%)")
print(f"Distinct RBP genes represented: {len(matched_genes):,}")

print(f"\nColumn population rate (all {n:,} rows):")
for c in pop_cols:
    print(f"  {c}: {pop_counts_all[c]:,} ({100 * pop_counts_all[c] / n:.1f}%)")

print(f"\nColumn population rate (within the {n_matched_rows:,} RBP-matched rows):")
for c in pop_cols:
    pct = 100 * pop_counts_matched[c] / n_matched_rows if n_matched_rows else 0
    print(f"  {c}: {pop_counts_matched[c]:,} ({pct:.1f}%)")
