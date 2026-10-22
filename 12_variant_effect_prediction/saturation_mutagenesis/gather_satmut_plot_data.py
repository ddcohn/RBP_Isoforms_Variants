"""
Single streaming pass over satmut_all_predictions.tsv (31,525,484 rows)
to gather just enough data to plot a quick analysis of TANGO, protGPS,
and DeepLoc results:

  - TANGO: a ~1-in-60 systematic sample of AGG and BETA (for histograms;
    exact histograms over 31.5M rows isn't needed for a "quick" look)
  - protGPS: exact running mean of all 12 compartment scores (cheap,
    one accumulator per column)
  - DeepLoc: exact frequency count of each individual localization
    category (Localizations field is pipe-joined when DeepLoc predicts
    more than one; each individual category is counted once per row it
    appears in)
"""
import csv
import json

csv.field_size_limit(2**31 - 1)

PATH = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_all_predictions.tsv"
OUT = "/u/home/d/ddcohn/_claude_satmut_plot_data.json"

PROTGPS_COLS = [
    "NUCLEAR_SPECKLE_Score", "P-BODY_Score", "PML-BDOY_Score",
    "POST_SYNAPTIC_DENSITY_Score", "STRESS_GRANULE_Score", "CHROMOSOME_Score",
    "NUCLEOLUS_Score", "NUCLEAR_PORE_COMPLEX_Score", "CAJAL_BODY_Score",
    "RNA_GRANULE_Score", "CELL_JUNCTION_Score", "TRANSCRIPTIONAL_Score",
]

SAMPLE_STRIDE = 60

tango_agg_sample = []
tango_beta_sample = []
protgps_sums = {c: 0.0 for c in PROTGPS_COLS}
protgps_valid_counts = {c: 0 for c in PROTGPS_COLS}
protgps_na_counts = {c: 0 for c in PROTGPS_COLS}
deeploc_loc_counts = {}
n_total = 0

with open(PATH, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for i, row in enumerate(reader):
        n_total += 1

        if i % SAMPLE_STRIDE == 0:
            try:
                tango_agg_sample.append(float(row["AGG"]))
                tango_beta_sample.append(float(row["BETA"]))
            except ValueError:
                pass

        for c in PROTGPS_COLS:
            try:
                protgps_sums[c] += float(row[c])
                protgps_valid_counts[c] += 1
            except ValueError:
                protgps_na_counts[c] += 1

        for loc in row["Localizations"].split("|"):
            loc = loc.strip()
            if loc:
                deeploc_loc_counts[loc] = deeploc_loc_counts.get(loc, 0) + 1

print(f"Total rows: {n_total}")

protgps_means = {
    c: (protgps_sums[c] / protgps_valid_counts[c] if protgps_valid_counts[c] else None)
    for c in PROTGPS_COLS
}
print(f"protGPS NA counts per column: {protgps_na_counts}")

out = {
    "n_total": n_total,
    "tango_agg_sample": tango_agg_sample,
    "tango_beta_sample": tango_beta_sample,
    "protgps_means": protgps_means,
    "protgps_na_counts": protgps_na_counts,
    "deeploc_loc_counts": deeploc_loc_counts,
}
with open(OUT, "w") as f:
    json.dump(out, f)

print(f"TANGO sample size: {len(tango_agg_sample)}")
print(f"protGPS means: {protgps_means}")
print(f"DeepLoc localization counts: {deeploc_loc_counts}")
print(f"Wrote {OUT}")
