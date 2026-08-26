import csv
import sys

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


def translate(cds):
    prot = []
    for i in range(0, len(cds) - len(cds) % 3, 3):
        codon = cds[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        prot.append(aa)
    return ''.join(prot)


def n_genes(ensembl_gene_field):
    genes = set(g.split(".")[0] for g in ensembl_gene_field.split(";") if g.strip())
    return len(genes)


def main():
    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    # baseline: what fraction of ALL rows (with a cds) have multi-gene ensembl_gene?
    baseline_total = 0
    baseline_multi = 0
    for row in rows:
        if not row.get("cds_sequence"):
            continue
        baseline_total += 1
        if n_genes(row.get("ensembl_gene", "")) > 1:
            baseline_multi += 1

    # identify the "other/unexplained" mismatches: exclude the 13 clean X-prefix cases
    x_prefix = 0
    unexplained = []
    other_explained_multigene = 0
    total_mismatch = 0

    for row in rows:
        cds = row.get("cds_sequence", "")
        prot = row.get("sequence", "")
        if not cds or not prot:
            continue
        translated = translate(cds)
        if translated == prot:
            continue
        total_mismatch += 1
        if translated == "X" + prot:
            x_prefix += 1
            continue
        ng = n_genes(row.get("ensembl_gene", ""))
        unexplained.append((row.get("uniprot_accession"), ng, len(prot), len(translated),
                             row.get("ensembl_gene", "")[:80]))

    n_unexplained = len(unexplained)
    n_multigene = sum(1 for u in unexplained if u[1] > 1)
    n_singlegene = n_unexplained - n_multigene

    print(f"Baseline: {baseline_multi} / {baseline_total} rows ({100*baseline_multi/baseline_total:.2f}%) "
          f"have a protein mapped to MULTIPLE genes")
    print()
    print(f"Total mismatches: {total_mismatch}")
    print(f"  X-prefix (explained, Ig/TCR partial-codon quirk): {x_prefix}")
    print(f"  Remaining 'unexplained' mismatches: {n_unexplained}")
    print(f"    of these, protein mapped to MULTIPLE genes: {n_multigene} ({100*n_multigene/n_unexplained:.1f}%)")
    print(f"    of these, protein mapped to a SINGLE gene:   {n_singlegene} ({100*n_singlegene/n_unexplained:.1f}%)")
    print()
    print(f"Enrichment: multi-gene rate in mismatches ({100*n_multigene/n_unexplained:.1f}%) "
          f"vs baseline rate in whole dataset ({100*baseline_multi/baseline_total:.2f}%)")
    print()

    print("Sample of multi-gene mismatches (uid, n_genes, protlen, translen, ensembl_gene):")
    shown = 0
    for u in unexplained:
        if u[1] > 1 and shown < 10:
            print(" ", u)
            shown += 1

    print()
    print("Sample of SINGLE-gene mismatches -- these are still genuinely unexplained (uid, n_genes, protlen, translen, ensembl_gene):")
    shown = 0
    for u in unexplained:
        if u[1] <= 1 and shown < 10:
            print(" ", u)
            shown += 1


if __name__ == "__main__":
    main()
