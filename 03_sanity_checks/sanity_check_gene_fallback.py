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


def main():
    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("transcript_selection_method") == "gene_level_fallback":
                rows.append(row)

    print(f"Total gene-level-fallback rows: {len(rows)}")

    exact = 0
    x_prefix = 0
    other_mismatch = []
    for row in rows:
        cds = row.get("cds_sequence", "")
        prot = row.get("sequence", "")
        if not cds or not prot:
            continue
        t = translate(cds)
        if t == prot:
            exact += 1
        elif t == "X" + prot:
            x_prefix += 1
        else:
            other_mismatch.append((row["uniprot_accession"], row.get("ensembl_transcript_used"),
                                     len(prot), len(t)))

    print(f"Exact match: {exact}")
    print(f"X-prefix pattern (Ig/TCR partial-codon quirk): {x_prefix}")
    print(f"Other mismatches: {len(other_mismatch)}")
    print()
    print("Other mismatch examples (uid, transcript, protlen, translen):")
    for ex in other_mismatch[:20]:
        print(" ", ex)

    # structural sanity
    starts_atg = sum(1 for r in rows if r.get("cds_sequence","").startswith("ATG"))
    div3 = sum(1 for r in rows if len(r.get("cds_sequence","")) % 3 == 0)
    ends_stop = sum(1 for r in rows if r.get("cds_sequence","")[-3:] in ("TAA","TAG","TGA"))
    print()
    print(f"Structural: starts ATG: {starts_atg}/{len(rows)}, div3: {div3}/{len(rows)}, ends stop: {ends_stop}/{len(rows)}")


if __name__ == "__main__":
    main()
