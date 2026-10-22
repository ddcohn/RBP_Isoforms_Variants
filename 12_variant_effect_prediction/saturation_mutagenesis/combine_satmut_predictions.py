"""
Combines the four satmut mutation-level result files into one wide
table keyed by mutation_id:
  - satmut_tango_all.tsv              (31,525,484 rows -- full mutation universe)
  - satmut_protgps_all.tsv            (31,525,484 rows -- full mutation universe)
  - satmut_deeploc_all.tsv            (31,525,484 rows -- full mutation universe)
  - satmut_slim_disruption_fixed.tsv  (25,718,153 rows -- only mutations that
    fall inside an existing ELM SLiM match; a *subset* of the mutation
    universe, and NOT 1:1 -- a mutation inside >1 overlapping motif has
    multiple rows there)

Approach: memory-safe sorted merge-join, not an in-RAM dict (31.5M-row
dicts across 4 files would be tens of GB). Each input is externally
sorted by mutation_id (GNU sort, LC_ALL=C for byte-order sorting that
matches Python string comparison), then all four sorted streams are
advanced in lockstep as a classic sorted merge-join.

The SLiM file is pre-aggregated to one row per mutation_id first
(semicolon-joined across any overlapping motifs) specifically so the
final join stays a clean 1:1 merge and the output row count is exactly
the tango/protgps/deeploc mutation universe (31,525,484) -- not
inflated by mutations that happen to sit inside more than one
overlapping ELM motif.

Run this on a compute node with real memory reserved (see
submit_combine_satmut.sh, h_data=32G), not the login node -- the
aggregation step builds an in-RAM dict over 25.7M rows and hit a
MemoryError on the shared login node even though the exact same code
had succeeded there once before (login node memory availability is not
reliable for this).

Caches the sorted tango/protgps/deeploc files across reruns (expensive,
unchanging) but always redoes the SLiM aggregation/sort (cheap, and the
input changed once already -- see fix_slim_gene_names.py).
"""
import csv
import os
import subprocess
import sys

csv.field_size_limit(sys.maxsize)

BASE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"
TANGO = f"{BASE}/satmut_tango_all.tsv"
PROTGPS = f"{BASE}/satmut_protgps_all.tsv"
DEEPLOC = f"{BASE}/satmut_deeploc_all.tsv"
SLIM = f"{BASE}/satmut_slim_disruption_fixed.tsv"
SLIM_AGG = f"{BASE}/satmut_slim_disruption_agg.tsv"
OUT_FILE = f"{BASE}/satmut_all_predictions.tsv"

SORT_ENV = {"LC_ALL": "C"}


def sort_by_first_col(infile, outfile):
    """Strips header, sorts body by column 1 (tab-delimited), writes header + sorted body."""
    with open(infile) as f:
        header = f.readline()
    with open(outfile, "w") as out:
        out.write(header)
    cmd = f"tail -n +2 {infile} | LC_ALL=C sort -t $'\\t' -k1,1 -S 4G"
    with open(outfile, "a") as out:
        subprocess.run(["bash", "-c", cmd], stdout=out, check=True, env={**SORT_ENV})
    return header.rstrip("\n").split("\t")


def read_header(path):
    with open(path) as f:
        return f.readline().rstrip("\n").split("\t")


def aggregate_slim(infile, outfile):
    """One row per mutation_id: semicolon-joined elm_identifiers, which ones disrupt, counts."""
    print("Aggregating SLiM disruption file to one row per mutation_id...")
    with open(infile, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        agg = {}
        for row in reader:
            mid = row["mutation_id"]
            entry = agg.setdefault(mid, {"all": [], "disrupted": []})
            entry["all"].append(row["elm_identifier"])
            if row["disrupts_motif"] == "Y":
                entry["disrupted"].append(row["elm_identifier"])

    print(f"  {len(agg)} unique mutation_ids with >=1 overlapping SLiM match")
    n_multi = sum(1 for v in agg.values() if len(v["all"]) > 1)
    print(f"  {n_multi} mutation_ids overlap more than one SLiM match instance")

    with open(outfile, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "mutation_id", "slim_n_motifs_overlapping", "slim_elm_identifiers",
            "slim_any_disrupts_motif", "slim_disrupted_elm_identifiers",
        ])
        for mid, v in agg.items():
            writer.writerow([
                mid,
                len(v["all"]),
                ";".join(v["all"]),
                "Y" if v["disrupted"] else "N",
                ";".join(v["disrupted"]),
            ])
    return len(agg)


def main():
    tango_sorted = TANGO + ".sorted"
    protgps_sorted = PROTGPS + ".sorted"
    deeploc_sorted = DEEPLOC + ".sorted"
    slim_sorted = SLIM_AGG + ".sorted"

    if all(os.path.exists(p) for p in (tango_sorted, protgps_sorted, deeploc_sorted)):
        print("Step 1 already done (sorted files present on disk) -- reusing them.")
        tango_header = read_header(tango_sorted)
        protgps_header = read_header(protgps_sorted)
        deeploc_header = read_header(deeploc_sorted)
    else:
        print("Step 1: sorting tango/protgps/deeploc by mutation_id...")
        tango_header = sort_by_first_col(TANGO, tango_sorted)
        print(f"  sorted TANGO ({len(tango_header)} cols)")
        protgps_header = sort_by_first_col(PROTGPS, protgps_sorted)
        print(f"  sorted protGPS ({len(protgps_header)} cols)")
        deeploc_header = sort_by_first_col(DEEPLOC, deeploc_sorted)
        print(f"  sorted DeepLoc ({len(deeploc_header)} cols)")

    print("\nStep 2: aggregating + sorting SLiM disruption file (always redone -- cheap, and input may have changed)...")
    aggregate_slim(SLIM, SLIM_AGG)
    slim_header = sort_by_first_col(SLIM_AGG, slim_sorted)
    print(f"  sorted SLiM-agg ({len(slim_header)} cols)")

    print("\nStep 3: 4-way sorted merge-join...")
    out_header = (
        tango_header + protgps_header[1:] + deeploc_header[1:] + slim_header[1:]
    )
    slim_blank = [""] * (len(slim_header) - 1)

    with open(tango_sorted, newline="") as ft, \
         open(protgps_sorted, newline="") as fp, \
         open(deeploc_sorted, newline="") as fd, \
         open(slim_sorted, newline="") as fs, \
         open(OUT_FILE, "w", newline="") as out:

        rt = csv.reader(ft, delimiter="\t")
        rp = csv.reader(fp, delimiter="\t")
        rd = csv.reader(fd, delimiter="\t")
        rs = csv.reader(fs, delimiter="\t")
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(out_header)

        next(rt)  # skip header row (sort_by_first_col writes header + sorted body)
        next(rp)
        next(rd)
        next(rs)

        row_p = next(rp)
        row_d = next(rd)
        row_s = next(rs, None)

        n_written = 0
        n_slim_matched = 0

        for row_t in rt:
            mid = row_t[0]

            if row_p[0] != mid:
                raise AssertionError(f"protGPS key mismatch at row {n_written}: {row_p[0]} != {mid}")
            if row_d[0] != mid:
                raise AssertionError(f"DeepLoc key mismatch at row {n_written}: {row_d[0]} != {mid}")

            slim_extra = slim_blank
            if row_s is not None:
                if row_s[0] == mid:
                    slim_extra = row_s[1:]
                    n_slim_matched += 1
                    row_s = next(rs, None)
                elif row_s[0] < mid:
                    while row_s is not None and row_s[0] < mid:
                        row_s = next(rs, None)
                    if row_s is not None and row_s[0] == mid:
                        slim_extra = row_s[1:]
                        n_slim_matched += 1
                        row_s = next(rs, None)

            writer.writerow(row_t + row_p[1:] + row_d[1:] + slim_extra)
            n_written += 1

            row_p = next(rp, None)
            row_d = next(rd, None)
            if row_p is None or row_d is None:
                break

        print(f"\nRows written: {n_written}")
        print(f"Rows with a SLiM match attached: {n_slim_matched}")

    print(f"\nWrote {OUT_FILE}")


if __name__ == "__main__":
    main()
