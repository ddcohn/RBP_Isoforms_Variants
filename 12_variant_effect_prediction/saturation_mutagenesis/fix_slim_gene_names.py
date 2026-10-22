"""
Fixes a gene-symbol mismatch found after the first join attempt:
compute_slim_disruption.py built mutation_id using the merged table's
`gene_symbol` column, but the satmut pipeline (TANGO/protGPS/DeepLoc)
uses asharma's `RBP_name` column from satmut_manifest.tsv -- these
differ for 68 of 2,335 proteins (gene aliases, e.g. O15446 =
CD3EAP/POLR1G), causing 233,985 of 9,297,650 SLiM-derived mutation_ids
to use the "wrong" name and silently fail to join.

Rewrites satmut_slim_disruption.tsv (mutation_id + gene_symbol columns
only) using the manifest's RBP_name as ground truth, keyed by
uniprot_accession, into satmut_slim_disruption_fixed.tsv. The
underlying SLiM match data itself (position, wt_aa, mut_aa,
elm_identifier, disrupts_motif, etc) is untouched and was never wrong
-- only the name used to build the ID string was.

One protein (P0DPA2/VSIG8) has no entry in the manifest at all -- it
was never part of the satmut mutagenesis scan, so its ~2,793 SLiM rows
are left with their original (unmappable) name and will never join to
the TANGO/protGPS/DeepLoc tables. That's expected, not a bug.
"""
import csv
import sys

csv.field_size_limit(sys.maxsize)

BASE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"
MANIFEST = f"{BASE}/satmut_manifest.tsv"
SLIM_IN = f"{BASE}/satmut_slim_disruption.tsv"
SLIM_OUT = f"{BASE}/satmut_slim_disruption_fixed.tsv"

name_map = {}
with open(MANIFEST, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        name_map[row["uniprot_accession"]] = row["RBP_name"]

print(f"Loaded {len(name_map)} accession -> RBP_name mappings from manifest")

n_fixed = 0
n_total = 0
n_unmapped = 0
with open(SLIM_IN, newline="") as f, open(SLIM_OUT, "w", newline="") as out:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for row in reader:
        n_total += 1
        acc = row["uniprot_accession"]
        correct_name = name_map.get(acc)
        if correct_name is None:
            n_unmapped += 1
        elif correct_name != row["gene_symbol"]:
            n_fixed += 1
            row["gene_symbol"] = correct_name
            row["mutation_id"] = f"{acc}|{correct_name}|{row['wt_aa']}{row['position']}{row['mut_aa']}"
        writer.writerow(row)

print(f"Total rows: {n_total}")
print(f"Rows with corrected gene_symbol/mutation_id: {n_fixed}")
print(f"Rows with accession not found in manifest (left as-is): {n_unmapped}")
print(f"Wrote {SLIM_OUT}")
