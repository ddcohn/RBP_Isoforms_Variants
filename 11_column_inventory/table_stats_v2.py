import csv
import sys
import json
import statistics as st

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
OT_RESULT = "/u/home/d/ddcohn/_claude_opentargets_result.json"

rows = []
with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows.append(row)

n = len(rows)
print(f"Total rows: {n}\n")


def numeric_summary(label, values):
    values = [v for v in values if v is not None]
    if not values:
        print(f"{label}: no data")
        return
    print(f"{label}: n={len(values)}  min={min(values)}  max={max(values)}  "
          f"mean={st.mean(values):.1f}  median={st.median(values)}")


def get_int(row, col):
    v = (row.get(col) or "").strip()
    return int(v) if v.isdigit() else None


def get_float(row, col):
    v = (row.get(col) or "").strip()
    try:
        return float(v)
    except ValueError:
        return None


# --- ID coverage ---
print("=== ID coverage ===")
for col in ["ensembl_gene", "ensembl_protein", "ncbi_geneid"]:
    have = sum(1 for r in rows if (r.get(col) or "").strip())
    print(f"  {col}: {have}/{n} ({100*have/n:.1f}%)")

# --- Protein sequence ---
print("\n=== Protein sequence (input, from base table) ===")
lens = [len(r["sequence"]) for r in rows if (r.get("sequence") or "").strip()]
numeric_summary("  length (aa)", lens)

# --- IDR ---
print("\n=== IDR (intrinsically disordered regions) ===")
idr_rows = [r for r in rows if (r.get("IDR_count") or "0").strip() not in ("", "0")]
print(f"  proteins with >=1 IDR: {len(idr_rows)}/{n} ({100*len(idr_rows)/n:.1f}%)")
numeric_summary("  IDR_count (of proteins with >=1)", [get_int(r, "IDR_count") for r in idr_rows])
numeric_summary("  IDR_avg_size (aa)", [get_float(r, "IDR_avg_size") for r in idr_rows])
numeric_summary("  IDR_total_size (aa)", [get_int(r, "IDR_total_size") for r in idr_rows])
total_idr_aa = sum(get_int(r, "IDR_total_size") or 0 for r in idr_rows)
print(f"  total disordered residues across dataset: {total_idr_aa}")

# --- RNA/CDS ---
print("\n=== RNA / CDS sequences ===")
has_rna = [r for r in rows if (r.get("rna_sequence") or "").strip()]
has_cds = [r for r in rows if (r.get("cds_sequence") or "").strip()]
print(f"  RNA sequence present: {len(has_rna)}/{n} ({100*len(has_rna)/n:.1f}%)")
print(f"  CDS sequence present: {len(has_cds)}/{n} ({100*len(has_cds)/n:.1f}%)")
numeric_summary("  RNA length (nt)", [len(r["rna_sequence"]) for r in has_rna])
numeric_summary("  CDS length (nt)", [len(r["cds_sequence"]) for r in has_cds])
# sanity: CDS should be ~3x protein length + stop codon
ratios = []
for r in has_cds:
    plen = len(r.get("sequence") or "")
    clen = len(r["cds_sequence"])
    if plen:
        ratios.append(clen / plen)
numeric_summary("  CDS_len / protein_len ratio", ratios)

method_counts = {}
for r in rows:
    m = (r.get("transcript_selection_method") or "(none)").strip() or "(none)"
    method_counts[m] = method_counts.get(m, 0) + 1
print("  transcript_selection_method breakdown:")
for m, c in sorted(method_counts.items(), key=lambda x: -x[1]):
    print(f"    {m}: {c} ({100*c/n:.1f}%)")

# --- PPI (3 methods) ---
print("\n=== STRING PPI (3 methods) ===")

def ppi_stats(col_all, col_in_df, label):
    with_partners = [r for r in rows if (r.get(col_all) or "").strip()]
    counts = [len((r.get(col_all) or "").split(";")) for r in with_partners]
    with_df_partners = [r for r in rows if (r.get(col_in_df) or "").strip()]
    total_edges = sum(counts)
    print(f"  [{label}]")
    print(f"    proteins with >=1 partner: {len(with_partners)}/{n} ({100*len(with_partners)/n:.1f}%)")
    numeric_summary("    partners/protein (of those with any)", counts)
    print(f"    total partner edges: {total_edges}")
    print(f"    proteins with >=1 partner ALSO in this dataframe: {len(with_df_partners)}/{n} ({100*len(with_df_partners)/n:.1f}%)")

ppi_stats("PPI_UniProt_Partners", "PPI_UniProt_Partners_in_Dataframe", "original method")
ppi_stats("PPI_UniProt_Partners_v125_local", "PPI_UniProt_Partners_in_Dataframe_v125_local", "v12.5 local file")
ppi_stats("PPI_UniProt_Partners_v125_api", "PPI_UniProt_Partners_in_Dataframe_v125_api", "v12.5 API")

# --- OpenTargets (not yet merged into table; read from separate result file) ---
print("\n=== OpenTargets disease associations (computed, NOT yet merged into table) ===")
try:
    with open(OT_RESULT) as f:
        ot = json.load(f)
    n_ot = len(ot)
    counts_per_protein = [len(v.get("diseaseId", [])) for v in ot.values() if v.get("diseaseId")]
    with_disease = len(counts_per_protein)
    total_assoc = sum(counts_per_protein)
    print(f"  candidate proteins queried: {n_ot}")
    print(f"  proteins with >=1 disease association: {with_disease}/{n_ot} ({100*with_disease/n_ot:.1f}%)")
    print(f"  total disease-datatype association rows: {total_assoc}")
    numeric_summary("  associations/protein (of those with any)", counts_per_protein)

    # datatype breakdown (evidence type: genetic_association, somatic_mutation, etc.)
    dt_counts = {}
    for v in ot.values():
        for dt in v.get("datatypeId", []):
            dt_counts[dt] = dt_counts.get(dt, 0) + 1
    print("  by evidence datatype (association rows):")
    for dt, c in sorted(dt_counts.items(), key=lambda x: -x[1]):
        print(f"    {dt}: {c}")

    # unique diseases
    all_disease_ids = set()
    for v in ot.values():
        all_disease_ids.update(v.get("diseaseId", []))
    print(f"  distinct diseases referenced: {len(all_disease_ids)}")
except FileNotFoundError:
    print("  result file not found.")
