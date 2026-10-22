import csv
csv.field_size_limit(2**31 - 1)

PATH = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_all_predictions.tsv"
OUT = "/u/home/d/ddcohn/_claude_satmut_column_stats.tsv"

with open(PATH, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    cols = reader.fieldnames
    non_empty = {c: 0 for c in cols}
    examples = {c: "" for c in cols}
    n_total = 0
    for row in reader:
        n_total += 1
        for c in cols:
            v = row[c]
            if v.strip():
                non_empty[c] += 1
                if not examples[c]:
                    examples[c] = v[:80]

print(f"Total rows: {n_total}")

with open(OUT, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["Column", "NonEmpty_rows", "Coverage_pct", "Example"])
    for c in cols:
        pct = 100.0 * non_empty[c] / n_total
        w.writerow([c, non_empty[c], f"{pct:.2f}", examples[c]])

print(f"Wrote {OUT}")
