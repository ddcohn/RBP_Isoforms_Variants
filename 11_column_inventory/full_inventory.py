import csv
import sys
import statistics as st

csv.field_size_limit(sys.maxsize)

RBP_TABLE = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
COSMIC_TABLE = "/u/project/kappel/ddcohn/Cosmic_Gene_Summary_v104_GRCh38.csv"

out_rows = []
out_rows.append(["Table", "Column", "% Populated", "N Populated", "Distinct Values", "Min", "Max", "Mean", "Median", "Notes"])


def add_generic(rows, n, table, col, notes=""):
    vals = [(r.get(col) or "").strip() for r in rows]
    populated = [v for v in vals if v]
    distinct = len(set(populated))
    out_rows.append([table, col, f"{100*len(populated)/n:.1f}%", len(populated), distinct, "", "", "", "", notes])


def add_length_stats(rows, n, table, col, notes=""):
    vals = [(r.get(col) or "").strip() for r in rows]
    populated = [v for v in vals if v]
    lens = [len(v) for v in populated]
    if lens:
        out_rows.append([table, col, f"{100*len(populated)/n:.1f}%", len(populated), "", min(lens), max(lens),
                          f"{st.mean(lens):.1f}", st.median(lens), notes + " (length in chars)"])
    else:
        out_rows.append([table, col, "0.0%", 0, "", "", "", "", "", notes])


def add_numeric_stats(rows, n, table, col, notes="", as_float=False, exclude_zero=True):
    raw = [(r.get(col) or "").strip() for r in rows]
    nums = []
    npop = 0
    for v in raw:
        if not v:
            continue
        try:
            val = float(v) if as_float else int(float(v))
        except ValueError:
            continue
        npop += 1
        if exclude_zero and val == 0:
            continue
        nums.append(val)
    if nums:
        out_rows.append([table, col, f"{100*npop/n:.1f}%", npop, "", min(nums), max(nums),
                          f"{st.mean(nums):.1f}", st.median(nums), notes])
    else:
        out_rows.append([table, col, f"{100*npop/n:.1f}%", npop, "", "", "", "", "", notes])


def add_semicolon_count_stats(rows, n, table, col, notes=""):
    vals = [(r.get(col) or "").strip() for r in rows]
    populated = [v for v in vals if v]
    counts = [len(v.split(";")) for v in populated]
    if counts:
        out_rows.append([table, col, f"{100*len(populated)/n:.1f}%", len(populated), "", min(counts), max(counts),
                          f"{st.mean(counts):.1f}", st.median(counts), notes + " (count = semicolon-separated items)"])
    else:
        out_rows.append([table, col, "0.0%", 0, "", "", "", "", "", notes])


def add_categorical(rows, n, table, col, notes=""):
    vals = [(r.get(col) or "").strip() or "(none)" for r in rows]
    counts = {}
    for v in vals:
        counts[v] = counts.get(v, 0) + 1
    top = sorted(counts.items(), key=lambda x: -x[1])
    top_str = "; ".join(f"{k}={c} ({100*c/n:.1f}%)" for k, c in top[:8])
    out_rows.append([table, col, "100.0%", n, len(counts), "", "", "", "", notes + " | top: " + top_str])


# ============ RBP TABLE ============
print("Reading RBP table...")
rbp_rows = []
with open(RBP_TABLE, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rbp_rows.append(row)
n = len(rbp_rows)
print(f"  {n} rows, {len(fieldnames)} columns")
T = "RBP isoform table (table_260823_with_rna.csv)"

add_generic(rbp_rows, n, T, "uniprot_accession", "primary key")
add_generic(rbp_rows, n, T, "ensembl_gene", "may contain multiple ;-separated IDs")
add_generic(rbp_rows, n, T, "ensembl_protein")
add_generic(rbp_rows, n, T, "ncbi_geneid")
add_length_stats(rbp_rows, n, T, "sequence", "protein sequence (aa)")
add_numeric_stats(rbp_rows, n, T, "IDR_count", "IDRs per protein")
add_numeric_stats(rbp_rows, n, T, "IDR_avg_size", "aa", as_float=True)
add_numeric_stats(rbp_rows, n, T, "IDR_total_size", "aa")
add_generic(rbp_rows, n, T, "IDR_range")
add_generic(rbp_rows, n, T, "IDR_discrete_seq")
add_generic(rbp_rows, n, T, "ensembl_transcript_used")
add_categorical(rbp_rows, n, T, "transcript_selection_method")
add_length_stats(rbp_rows, n, T, "rna_sequence", "nt")
add_length_stats(rbp_rows, n, T, "cds_sequence", "nt")
add_semicolon_count_stats(rbp_rows, n, T, "PPI_UniProt_Partners", "original method")
add_semicolon_count_stats(rbp_rows, n, T, "PPI_UniProt_Partners_in_Dataframe")
add_semicolon_count_stats(rbp_rows, n, T, "PPI_UniProt_Partners_v125_local")
add_semicolon_count_stats(rbp_rows, n, T, "PPI_UniProt_Partners_in_Dataframe_v125_local")
add_semicolon_count_stats(rbp_rows, n, T, "PPI_UniProt_Partners_v125_api")
add_semicolon_count_stats(rbp_rows, n, T, "PPI_UniProt_Partners_in_Dataframe_v125_api")
add_semicolon_count_stats(rbp_rows, n, T, "GO_Cellular_Component", "GO_ID:term pairs")
add_semicolon_count_stats(rbp_rows, n, T, "GO_Biological_Process", "GO_ID:term pairs")
add_semicolon_count_stats(rbp_rows, n, T, "GO_Molecular_Function", "GO_ID:term pairs")
add_semicolon_count_stats(rbp_rows, n, T, "Tissue_Expression_Bulk", "tissue|source|median|unit entries")
add_semicolon_count_stats(rbp_rows, n, T, "OpenTargets_DiseaseId", "disease associations")
add_semicolon_count_stats(rbp_rows, n, T, "OpenTargets_DatatypeId")
add_semicolon_count_stats(rbp_rows, n, T, "OpenTargets_Score")

# ============ COSMIC TABLE ============
print("Reading COSMIC table...")
cosmic_rows = []
with open(COSMIC_TABLE, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    for row in reader:
        cosmic_rows.append(row)
n2 = len(cosmic_rows)
print(f"  {n2} rows")
T2 = "COSMIC gene summary (Cosmic_Gene_Summary_v104_GRCh38.csv)"

add_generic(cosmic_rows, n2, T2, "COSMIC_GENE_ID", "primary key")
add_generic(cosmic_rows, n2, T2, "GENE_SYMBOL")
add_generic(cosmic_rows, n2, T2, "ENSEMBL_GENE_ID", "version-stripped ENSG")
add_generic(cosmic_rows, n2, T2, "ENTREZ_ID")
add_generic(cosmic_rows, n2, T2, "HGNC_ID")
add_categorical(cosmic_rows, n2, T2, "IN_CANCER_CENSUS")
add_categorical(cosmic_rows, n2, T2, "IS_EXPERT_CURATED")
add_numeric_stats(cosmic_rows, n2, T2, "Total_Mutation_Rows", exclude_zero=False)
add_numeric_stats(cosmic_rows, n2, T2, "Distinct_Genomic_Mutations", exclude_zero=False)
add_numeric_stats(cosmic_rows, n2, T2, "Distinct_Samples", exclude_zero=False)
add_numeric_stats(cosmic_rows, n2, T2, "Confirmed_Somatic_Rows", exclude_zero=False)
add_numeric_stats(cosmic_rows, n2, T2, "Missense_Rows", exclude_zero=False)
add_numeric_stats(cosmic_rows, n2, T2, "Truncating_Rows", exclude_zero=False)
add_numeric_stats(cosmic_rows, n2, T2, "Synonymous_Rows", exclude_zero=False)

OUT = "/u/home/d/ddcohn/_claude_full_inventory.tsv"
with open(OUT, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(out_rows)

print(f"\nWrote {len(out_rows)-1} rows to {OUT}")
for r in out_rows:
    print("\t".join(str(x) for x in r))
