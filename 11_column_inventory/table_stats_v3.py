import csv
import sys
import json
import statistics as st

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
OT_RESULT = "/u/home/d/ddcohn/_claude_opentargets_result.json"
OUT = "/u/home/d/ddcohn/_claude_column_stats.tsv"

rows = []
with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)

n = len(rows)
out_rows = []
out_rows.append(["Column", "% Populated", "N Populated", "Distinct Values", "Min", "Max", "Mean", "Median", "Notes"])


def fmt(x):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.2f}"
    return str(x)


def add_generic(col, notes=""):
    vals = [(r.get(col) or "").strip() for r in rows]
    populated = [v for v in vals if v]
    distinct = len(set(populated))
    out_rows.append([col, f"{100*len(populated)/n:.1f}%", len(populated), distinct, "", "", "", "", notes])


def add_length_stats(col, notes=""):
    vals = [(r.get(col) or "").strip() for r in rows]
    populated = [v for v in vals if v]
    lens = [len(v) for v in populated]
    if lens:
        out_rows.append([col, f"{100*len(populated)/n:.1f}%", len(populated), "", min(lens), max(lens),
                          f"{st.mean(lens):.1f}", st.median(lens), notes + " (length in chars)"])
    else:
        out_rows.append([col, "0.0%", 0, "", "", "", "", "", notes])


def add_numeric_stats(col, notes="", as_float=False):
    raw = [(r.get(col) or "").strip() for r in rows]
    nums = []
    for v in raw:
        if not v:
            continue
        try:
            val = float(v) if as_float else int(v)
        except ValueError:
            continue
        if val != 0:
            nums.append(val)
    if nums:
        out_rows.append([col, f"{100*len(nums)/n:.1f}%", len(nums), "", min(nums), max(nums),
                          f"{st.mean(nums):.1f}", st.median(nums), notes])
    else:
        out_rows.append([col, "0.0%", 0, "", "", "", "", "", notes])


def add_semicolon_count_stats(col, notes=""):
    vals = [(r.get(col) or "").strip() for r in rows]
    populated = [v for v in vals if v]
    counts = [len(v.split(";")) for v in populated]
    if counts:
        out_rows.append([col, f"{100*len(populated)/n:.1f}%", len(populated), "", min(counts), max(counts),
                          f"{st.mean(counts):.1f}", st.median(counts), notes + " (count = semicolon-separated partners)"])
    else:
        out_rows.append([col, "0.0%", 0, "", "", "", "", "", notes])


def add_categorical(col, notes=""):
    vals = [(r.get(col) or "").strip() or "(none)" for r in rows]
    counts = {}
    for v in vals:
        counts[v] = counts.get(v, 0) + 1
    top = sorted(counts.items(), key=lambda x: -x[1])
    top_str = "; ".join(f"{k}={c} ({100*c/n:.1f}%)" for k, c in top)
    out_rows.append([col, "100.0%", n, len(counts), "", "", "", "", notes + " | breakdown: " + top_str])


# --- base / ID columns ---
add_generic("uniprot_accession", "primary key, one row per protein")
add_generic("ensembl_gene", "may contain multiple ;-separated IDs")
add_generic("ensembl_protein", "may contain multiple ;-separated IDs")
add_generic("ncbi_geneid", "may contain multiple ;-separated IDs")
add_length_stats("sequence", "input protein sequence (aa)")

# --- IDR ---
add_numeric_stats("IDR_count", "proteins with 0 excluded from stats; count of IDRs per protein")
add_numeric_stats("IDR_avg_size", "avg size of IDRs within a protein (aa)", as_float=True)
add_numeric_stats("IDR_total_size", "total disordered residues in protein (aa)")
add_generic("IDR_range", "list of [start,end] ranges; populated only when IDR_count>0")
add_generic("IDR_discrete_seq", "list of disordered subsequences; populated only when IDR_count>0")

# --- RNA/CDS pipeline ---
add_generic("ensembl_transcript_used", "transcript ID selected as source of rna/cds_sequence")
add_categorical("transcript_selection_method", "method used to select ensembl_transcript_used")
add_length_stats("rna_sequence", "full mRNA (nt)")
add_length_stats("cds_sequence", "coding sequence (nt)")

# --- STRING PPI, 3 methods x 2 columns each ---
add_semicolon_count_stats("PPI_UniProt_Partners", "original method, all STRING partners")
add_semicolon_count_stats("PPI_UniProt_Partners_in_Dataframe", "original method, partners also in this dataset")
add_semicolon_count_stats("PPI_UniProt_Partners_v125_local", "STRING v12.5 local file, all partners")
add_semicolon_count_stats("PPI_UniProt_Partners_in_Dataframe_v125_local", "STRING v12.5 local file, partners also in this dataset")
add_semicolon_count_stats("PPI_UniProt_Partners_v125_api", "STRING v12.5 API, all partners")
add_semicolon_count_stats("PPI_UniProt_Partners_in_Dataframe_v125_api", "STRING v12.5 API, partners also in this dataset")

# --- OpenTargets: not yet a column in the table; computed separately ---
try:
    with open(OT_RESULT) as f:
        ot = json.load(f)
    n_ot_candidates = len(ot)
    counts_per_protein = {uid: len(v.get("diseaseId", [])) for uid, v in ot.items()}
    nonzero = [c for c in counts_per_protein.values() if c > 0]
    if nonzero:
        out_rows.append([
            "OpenTargets_disease_associations (NOT YET IN TABLE - pending merge)",
            f"{100*len(nonzero)/n_ot_candidates:.1f}% of {n_ot_candidates} candidates",
            len(nonzero), "", min(nonzero), max(nonzero),
            f"{st.mean(nonzero):.1f}", st.median(nonzero),
            "count of disease-datatype association rows per protein; stored in _claude_opentargets_result.json"
        ])
    all_diseases = set()
    for v in ot.values():
        all_diseases.update(v.get("diseaseId", []))
    out_rows.append([
        "OpenTargets_distinct_diseases (NOT YET IN TABLE)", "", "", len(all_diseases), "", "", "", "",
        "distinct disease IDs referenced across all proteins"
    ])
except FileNotFoundError:
    pass

with open(OUT, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(out_rows)

print(f"Wrote {len(out_rows)-1} data rows to {OUT}")
print()
for r in out_rows:
    print("\t".join(fmt(x) for x in r))
