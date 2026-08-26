import csv
import sys
from collections import Counter, defaultdict

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

KNOWN_SELENOPROTEINS = {
    "O60613","P07203","P18283","P22352","P36969","P49895","P49908","P55073","P59796",
    "P59797","P62341","P63302","Q16881","Q86VQ6","Q8IZQ5","Q8WWX9","Q92813","Q99611",
    "Q9BQE4","Q9BVL4","Q9C0D9","Q9NNW7","Q9NZV5","Q9Y6D0",
}


def translate(cds):
    prot = []
    for i in range(0, len(cds) - len(cds) % 3, 3):
        codon = cds[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        prot.append(aa)
    return ''.join(prot)


def n_genes(field):
    return len(set(g.split(".")[0] for g in field.split(";") if g.strip()))


def main():
    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    print(f"Total rows in table: {len(rows)}")

    method_counts = Counter(r.get("transcript_selection_method", "") for r in rows)
    print()
    print("=== Rows by transcript_selection_method ===")
    for m, c in method_counts.most_common():
        print(f"  {m!r}: {c}")

    # main sanity pass
    stats = defaultdict(lambda: Counter())
    overall = Counter()
    mismatch_examples = defaultdict(list)

    for row in rows:
        method = row.get("transcript_selection_method", "")
        cds = row.get("cds_sequence", "")
        prot = row.get("sequence", "")
        if not cds or not prot:
            stats[method]["no_cds_or_protein"] += 1
            overall["no_cds_or_protein"] += 1
            continue

        t = translate(cds)
        acc = row.get("uniprot_accession", "")
        is_multigene = n_genes(row.get("ensembl_gene", "")) > 1
        is_seleno = acc in KNOWN_SELENOPROTEINS

        if t == prot:
            cat = "exact_match"
        elif t == "X" + prot:
            cat = "x_prefix_benign"
        elif is_seleno:
            cat = "selenoprotein_explained"
        elif is_multigene:
            cat = "multigene_ambiguous"
        else:
            cat = "unexplained_mismatch"
            if len(mismatch_examples[method]) < 5:
                mismatch_examples[method].append((acc, len(prot), len(t)))

        stats[method][cat] += 1
        overall[cat] += 1

    print()
    print("=== Overall (all methods combined) ===")
    total_checked = sum(overall.values()) - overall["no_cds_or_protein"]
    for cat in ["exact_match", "x_prefix_benign", "selenoprotein_explained",
                "multigene_ambiguous", "unexplained_mismatch"]:
        c = overall[cat]
        pct = 100 * c / total_checked if total_checked else 0
        print(f"  {cat}: {c} ({pct:.1f}% of checked)")
    print(f"  no_cds_or_protein (no RNA available): {overall['no_cds_or_protein']}")
    print(f"  TOTAL checked (had both cds and protein): {total_checked}")

    print()
    print("=== Breakdown by transcript_selection_method ===")
    for method in method_counts:
        s = stats[method]
        checked = sum(s.values()) - s["no_cds_or_protein"]
        print(f"\n  Method: {method!r} (total rows: {method_counts[method]})")
        if checked == 0:
            print(f"    (no rows with both cds and protein)")
            continue
        for cat in ["exact_match", "x_prefix_benign", "selenoprotein_explained",
                    "multigene_ambiguous", "unexplained_mismatch"]:
            c = s[cat]
            pct = 100 * c / checked if checked else 0
            print(f"    {cat}: {c} ({pct:.1f}%)")
        if mismatch_examples[method]:
            print(f"    unexplained_mismatch examples: {mismatch_examples[method]}")

    # rna-contains-cds check across whole table
    print()
    print("=== rna_sequence contains cds_sequence (whole table) ===")
    both = 0
    not_substr = 0
    for row in rows:
        rna = row.get("rna_sequence", "")
        cds = row.get("cds_sequence", "")
        if rna and cds:
            both += 1
            if cds not in rna:
                not_substr += 1
    print(f"Rows with both: {both}; cds NOT found in rna: {not_substr}")


if __name__ == "__main__":
    main()
