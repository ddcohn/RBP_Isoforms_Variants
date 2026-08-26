import csv
import sys

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"


def jaccard(a, b):
    if not a and not b:
        return None
    u = a | b
    if not u:
        return None
    return len(a & b) / len(u)


def main():
    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in rows_iter(reader):
            rows.append(row)

    # exact vs jaccard comparison: v125_local vs v125_api (same version, two methods)
    exact_match = 0
    both_empty = 0
    one_empty_other_not = 0
    jaccards = []
    low_overlap_examples = []

    for row in rows:
        a = set(x for x in row.get("PPI_UniProt_Partners_v125_local", "").split(";") if x)
        b = set(x for x in row.get("PPI_UniProt_Partners_v125_api", "").split(";") if x)
        if not a and not b:
            both_empty += 1
            continue
        if bool(a) != bool(b):
            one_empty_other_not += 1
            continue
        if a == b:
            exact_match += 1
        j = jaccard(a, b)
        if j is not None:
            jaccards.append(j)
            if j < 0.5 and len(low_overlap_examples) < 8:
                low_overlap_examples.append((row["uniprot_accession"], len(a), len(b), len(a & b), j))

    total_comparable = len(jaccards)
    print(f"Total rows: {len(rows)}")
    print(f"Both v12.5 methods empty: {both_empty}")
    print(f"One method has partners, other doesn't: {one_empty_other_not}")
    print(f"Both non-empty, EXACT same partner set: {exact_match}")
    print(f"Both non-empty, comparable via Jaccard: {total_comparable}")
    if jaccards:
        jaccards.sort()
        print(f"Jaccard similarity (v125_local vs v125_api): "
              f"min={jaccards[0]:.3f}, median={jaccards[len(jaccards)//2]:.3f}, "
              f"mean={sum(jaccards)/len(jaccards):.3f}, max={jaccards[-1]:.3f}")
        n_perfect = sum(1 for j in jaccards if j == 1.0)
        n_high = sum(1 for j in jaccards if j >= 0.9)
        n_low = sum(1 for j in jaccards if j < 0.5)
        print(f"  Jaccard == 1.0 (identical sets): {n_perfect} ({100*n_perfect/total_comparable:.1f}%)")
        print(f"  Jaccard >= 0.9: {n_high} ({100*n_high/total_comparable:.1f}%)")
        print(f"  Jaccard < 0.5: {n_low} ({100*n_low/total_comparable:.1f}%)")
    print()
    print("Examples with low overlap (uid, len_local, len_api, len_intersection, jaccard):")
    for ex in low_overlap_examples:
        print(" ", ex)

    print()
    print("=== v12.0 vs v12.5 (live API) comparison, same access method, different STRING version ===")
    v120_only_new = 0  # protein has partners in v12.0 but not v12.5 api at all
    both_versions = []
    for row in rows:
        old = set(x for x in row.get("PPI_UniProt_Partners", "").split(";") if x)
        new = set(x for x in row.get("PPI_UniProt_Partners_v125_api", "").split(";") if x)
        if old and not new:
            v120_only_new += 1
        if old and new:
            both_versions.append(jaccard(old, new))
    print(f"Rows with v12.0 partners but ZERO v12.5-api partners: {v120_only_new}")
    if both_versions:
        both_versions.sort()
        n = len(both_versions)
        print(f"Jaccard (v12.0 vs v12.5-api) across {n} rows with both: "
              f"min={both_versions[0]:.3f}, median={both_versions[n//2]:.3f}, mean={sum(both_versions)/n:.3f}")


def rows_iter(reader):
    for row in reader:
        yield row


if __name__ == "__main__":
    main()
