import csv
import sys
import ast
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

PATH = "/u/project/kappel/RBP/Misc/Classical_RBDs_Annotated.csv"

COLS = ["Domains", "Domains_count", "Domains_avg_size", "Domains_total_size",
        "Domains_range", "Domains_discrete_seq",
        "Has_RBD", "RBD_name_ranges", "RBD_all_ranges", "RBD_names"]

per_protein = defaultdict(list)
n_rows = 0
with open(PATH, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        n_rows += 1
        uid = row.get("uniprot_id", "")
        vals = tuple(row.get(c, "") for c in COLS)
        per_protein[uid].append(vals)
        if n_rows % 5000 == 0:
            print(f"  ... {n_rows} rows read", flush=True)

print(f"\nTotal rows: {n_rows}")
print(f"Distinct uniprot_id: {len(per_protein)}")

n_consistent = n_inconsistent = 0
example_inconsistent = None
for uid, rows in per_protein.items():
    if len(set(rows)) == 1:
        n_consistent += 1
    else:
        n_inconsistent += 1
        if example_inconsistent is None:
            example_inconsistent = (uid, rows[:3])

print(f"Proteins with consistent Domains/RBD data across duplicate rows: {n_consistent}")
print(f"Proteins with INCONSISTENT data across duplicate rows: {n_inconsistent}")
if example_inconsistent:
    print("Example inconsistent:", example_inconsistent)

# --- overlap investigation ---
print("\n=== Domain overlap handling investigation ===")
checked = 0
overlap_examples = []
for uid, rows in per_protein.items():
    vals = rows[0]
    d = dict(zip(COLS, vals))
    try:
        dcount = int(float(d["Domains_count"])) if d["Domains_count"] else 0
    except ValueError:
        dcount = 0
    if dcount < 2:
        continue
    try:
        ranges = ast.literal_eval(d["Domains_range"]) if d["Domains_range"] else []
    except (ValueError, SyntaxError):
        continue
    if not isinstance(ranges, list) or len(ranges) < 2:
        continue
    # check for overlap between consecutive ranges (sorted by start)
    sorted_ranges = sorted(ranges, key=lambda r: r[0])
    has_overlap = any(sorted_ranges[i][1] >= sorted_ranges[i+1][0] for i in range(len(sorted_ranges)-1))
    checked += 1
    if has_overlap:
        try:
            discrete = ast.literal_eval(d["Domains_discrete_seq"]) if d["Domains_discrete_seq"] else []
        except (ValueError, SyntaxError):
            discrete = "PARSE_ERROR"
        overlap_examples.append((uid, d["Domains"], ranges, discrete))
    if len(overlap_examples) >= 5:
        break

print(f"Proteins checked with >=2 domains: {checked}")
print(f"Examples found with OVERLAPPING domain ranges: {len(overlap_examples)}")
for uid, names, ranges, discrete in overlap_examples:
    print(f"\n  {uid}")
    print(f"    Domains: {names}")
    print(f"    Domains_range: {ranges}")
    print(f"    Domains_discrete_seq (n={len(discrete) if isinstance(discrete, list) else 'N/A'}): {str(discrete)[:300]}")
