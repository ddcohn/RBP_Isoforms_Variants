import csv
import sys
import re
from collections import Counter

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
UNIPROT_RE = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$')


def main():
    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    total = len(rows)
    with_partners = 0
    with_in_df = 0
    self_ref = 0
    not_subset = 0
    bad_format = 0
    dup_within_row = 0
    partner_counts = []
    in_df_counts = []

    partners_in_df_by_uid = {}  # uid -> set of partner uids (only in-dataframe list)

    for row in rows:
        uid = row["uniprot_accession"].strip()
        p_str = row.get("PPI_UniProt_Partners", "")
        d_str = row.get("PPI_UniProt_Partners_in_Dataframe", "")
        partners = [x for x in p_str.split(";") if x]
        in_df = [x for x in d_str.split(";") if x]

        if partners:
            with_partners += 1
            partner_counts.append(len(partners))
        if in_df:
            with_in_df += 1
            in_df_counts.append(len(in_df))

        if len(set(partners)) != len(partners):
            dup_within_row += 1

        if uid in partners or uid in in_df:
            self_ref += 1

        if not set(in_df) <= set(partners):
            not_subset += 1

        for p in partners:
            if not UNIPROT_RE.match(p):
                bad_format += 1
                break

        partners_in_df_by_uid[uid] = set(in_df)

    print(f"Total rows: {total}")
    print(f"Rows with >=1 partner: {with_partners}")
    print(f"Rows with >=1 in-dataframe partner: {with_in_df}")
    print()
    print("=== Structural sanity ===")
    print(f"Rows with a self-reference (protein listed as its own partner): {self_ref}")
    print(f"Rows where in_Dataframe is NOT a subset of full partner list: {not_subset}")
    print(f"Rows with duplicate entries within the partner list: {dup_within_row}")
    print(f"Rows containing a malformed-looking UniProt accession: {bad_format}")
    print()
    if partner_counts:
        print(f"Partner count per protein: min={min(partner_counts)}, max={max(partner_counts)}, "
              f"mean={sum(partner_counts)/len(partner_counts):.1f}, median={sorted(partner_counts)[len(partner_counts)//2]}")
    if in_df_counts:
        print(f"In-dataframe partner count: min={min(in_df_counts)}, max={max(in_df_counts)}, "
              f"mean={sum(in_df_counts)/len(in_df_counts):.1f}, median={sorted(in_df_counts)[len(in_df_counts)//2]}")
    print()

    # --- symmetry check on in-dataframe edges ---
    print("=== Symmetry check (A->B in-dataframe implies B->A in-dataframe) ===")
    total_edges = 0
    reciprocated = 0
    examples_non_recip = []
    for uid, partners in partners_in_df_by_uid.items():
        for p in partners:
            total_edges += 1
            if p in partners_in_df_by_uid and uid in partners_in_df_by_uid[p]:
                reciprocated += 1
            else:
                if len(examples_non_recip) < 10:
                    examples_non_recip.append((uid, p))
    print(f"Total directed in-dataframe edges: {total_edges}")
    print(f"Reciprocated (B also lists A): {reciprocated} ({100*reciprocated/total_edges:.1f}%)")
    print(f"Non-reciprocated examples (A, B) where A->B but not B->A:")
    for ex in examples_non_recip:
        print(" ", ex)
    print()

    # --- known biology spot check: look for a couple of well-known RBPs ---
    print("=== Spot check known RBPs ===")
    known = {
        "Q13148": "TDP-43 (TARDBP)",
        "P35637": "FUS",
        "Q15717": "ELAVL1 (HuR)",
        "O43318": "MAP3K7 (TAK1, sanity non-RBP control)",
    }
    by_uid = {r["uniprot_accession"].strip(): r for r in rows}
    for uid, label in known.items():
        if uid in by_uid:
            r = by_uid[uid]
            partners = [x for x in r.get("PPI_UniProt_Partners","").split(";") if x]
            print(f"{uid} ({label}): {len(partners)} partners. First 10: {partners[:10]}")
        else:
            print(f"{uid} ({label}): not present in this dataset")


if __name__ == "__main__":
    main()
