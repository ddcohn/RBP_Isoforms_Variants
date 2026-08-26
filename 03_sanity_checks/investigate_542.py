import csv
import sys
from collections import Counter

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


def n_genes(field):
    return len(set(g.split(".")[0] for g in field.split(";") if g.strip()))


def main():
    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    single_gene_mismatches = []
    for row in rows:
        cds = row.get("cds_sequence", "")
        prot = row.get("sequence", "")
        if not cds or not prot:
            continue
        translated = translate(cds)
        if translated == prot or translated == "X" + prot:
            continue
        if n_genes(row.get("ensembl_gene", "")) > 1:
            continue
        single_gene_mismatches.append((row, translated, prot))

    print(f"Total single-gene mismatches: {len(single_gene_mismatches)}")
    print()

    # length relationship
    shorter = 0
    longer = 0
    same_len_diff = 0
    for row, translated, prot in single_gene_mismatches:
        if len(translated) < len(prot):
            shorter += 1
        elif len(translated) > len(prot):
            longer += 1
        else:
            same_len_diff += 1
    print("=== Length relationship (translated vs stored protein) ===")
    print(f"Translated SHORTER than protein (early stop / truncation): {shorter}")
    print(f"Translated LONGER than protein (stop codon read through / extra residues): {longer}")
    print(f"Same length, different residues (frame/substitution): {same_len_diff}")
    print()

    # gene/description keyword categorization
    cats = Counter()
    examples = {}
    for row, translated, prot in single_gene_mismatches:
        desc = (row.get("gene_description", "") + " " + row.get("Description", "") + " "
                + row.get("gene_symbol", "") + " " + row.get("Name", "")).lower()
        if any(k in desc for k in ["selenoprotein", "selenocysteine"]):
            cat = "Selenoprotein (UGA readthrough)"
        elif any(k in desc for k in ["immunoglobulin", "ig kappa", "ig lambda", "ig heavy",
                                       "t cell receptor", "trbv","trbj","trbc","trav","traj","trac",
                                       "igkv","igkj","igkc","iglv","iglj","iglc","ighv","ighj","ighd","ighg","ighm","igha"]):
            cat = "Immunoglobulin/TCR gene segment"
        elif "mitochondri" in desc or row.get("gene_symbol","").upper().startswith("MT-"):
            cat = "Mitochondrial gene"
        elif "readthrough" in desc:
            cat = "Annotated readthrough gene"
        elif "pseudogene" in desc:
            cat = "Pseudogene"
        else:
            cat = "Other/unclassified"
        cats[cat] += 1
        if len(examples.setdefault(cat, [])) < 6:
            examples[cat].append((row.get("uniprot_accession"), row.get("gene_symbol"),
                                    len(prot), len(translated), row.get("ensembl_transcript_used")))

    print("=== Functional annotation categorization ===")
    for cat, count in cats.most_common():
        print(f"{cat}: {count}")
        for ex in examples[cat]:
            print(f"    {ex}")
    print()

    # for the biggest unexplained bucket, check first few residues for selenocysteine-like early stop
    print("=== Check: does translated sequence end right where a TGA (possible selenocysteine) would be? ===")
    n_tga_stop = 0
    for row, translated, prot in single_gene_mismatches:
        cds = row.get("cds_sequence", "")
        if len(translated) < len(prot):
            stop_pos = len(translated) * 3
            codon_at_stop = cds[stop_pos:stop_pos+3]
            if codon_at_stop == "TGA":
                n_tga_stop += 1
    print(f"Of the 'translated shorter' cases, number where the stopping codon is specifically TGA "
          f"(the selenocysteine-recoded codon): {n_tga_stop} / {shorter}")


if __name__ == "__main__":
    main()
