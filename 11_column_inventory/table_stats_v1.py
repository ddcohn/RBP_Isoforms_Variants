import csv
import sys
import json

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

# --- IDR ---
idr_rows = [r for r in rows if (r.get("IDR_count") or "0").strip() not in ("", "0")]
print(f"IDR: {len(idr_rows)}/{n} proteins have >=1 IDR ({100*len(idr_rows)/n:.1f}%)")
total_idr = sum(int(r["IDR_count"]) for r in idr_rows if r["IDR_count"].strip().isdigit())
print(f"  total IDRs across dataset: {total_idr}")

# --- RNA/CDS ---
has_rna = [r for r in rows if (r.get("rna_sequence") or "").strip()]
has_cds = [r for r in rows if (r.get("cds_sequence") or "").strip()]
print(f"\nRNA sequence: {len(has_rna)}/{n} ({100*len(has_rna)/n:.1f}%)")
print(f"CDS sequence: {len(has_cds)}/{n} ({100*len(has_cds)/n:.1f}%)")

method_counts = {}
for r in rows:
    m = (r.get("transcript_selection_method") or "(none)").strip() or "(none)"
    method_counts[m] = method_counts.get(m, 0) + 1
print("  transcript_selection_method breakdown:")
for m, c in sorted(method_counts.items(), key=lambda x: -x[1]):
    print(f"    {m}: {c} ({100*c/n:.1f}%)")

# --- PPI (3 methods) ---
def ppi_stats(col_all, col_in_df, label):
    with_partners = [r for r in rows if (r.get(col_all) or "").strip()]
    counts = [len((r.get(col_all) or "").split(";")) for r in with_partners]
    with_df_partners = [r for r in rows if (r.get(col_in_df) or "").strip()]
    total_edges = sum(counts)
    avg = total_edges / len(with_partners) if with_partners else 0
    print(f"\nPPI [{label}]:")
    print(f"  proteins with >=1 partner: {len(with_partners)}/{n} ({100*len(with_partners)/n:.1f}%)")
    print(f"  total partner edges: {total_edges}, avg partners/protein (of those with any): {avg:.2f}")
    print(f"  proteins with >=1 partner ALSO in this dataframe: {len(with_df_partners)}/{n} ({100*len(with_df_partners)/n:.1f}%)")

ppi_stats("PPI_UniProt_Partners", "PPI_UniProt_Partners_in_Dataframe", "original method")
ppi_stats("PPI_UniProt_Partners_v125_local", "PPI_UniProt_Partners_in_Dataframe_v125_local", "v12.5 local file")
ppi_stats("PPI_UniProt_Partners_v125_api", "PPI_UniProt_Partners_in_Dataframe_v125_api", "v12.5 API")

# --- OpenTargets (not yet merged into table; read from separate result file) ---
try:
    with open(OT_RESULT) as f:
        ot = json.load(f)
    n_ot = len(ot)
    with_disease = [v for v in ot.values() if v.get("diseaseId")]
    total_assoc = sum(len(v.get("diseaseId", [])) for v in ot.values())
    print(f"\nOpenTargets (computed, NOT yet merged into table):")
    print(f"  candidate proteins queried: {n_ot}")
    print(f"  proteins with >=1 disease association: {len(with_disease)}/{n_ot} ({100*len(with_disease)/n_ot:.1f}%)")
    print(f"  total disease-datatype association rows: {total_assoc}")
    avg_assoc = total_assoc / len(with_disease) if with_disease else 0
    print(f"  avg associations/protein (of those with any): {avg_assoc:.1f}")
except FileNotFoundError:
    print("\nOpenTargets result file not found.")
