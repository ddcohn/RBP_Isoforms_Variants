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


def main():
    total = 0
    have_both = 0
    starts_atg = 0
    ends_stop = 0
    div3 = 0
    non_acgt = 0
    exact_match = 0
    mismatch = 0
    mismatch_examples = []
    len_mismatch_only = 0

    rna_missing_cds_substr = 0
    rna_and_cds_present = 0

    valid_bases = set("ACGTN")

    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            total += 1
            cds = row.get("cds_sequence", "")
            rna = row.get("rna_sequence", "")
            prot = row.get("sequence", "")

            if cds:
                if not (set(cds.upper()) <= valid_bases):
                    non_acgt += 1
                if cds.startswith("ATG"):
                    starts_atg += 1
                if len(cds) % 3 == 0:
                    div3 += 1
                if cds[-3:] in ("TAA", "TAG", "TGA"):
                    ends_stop += 1

            if cds and prot:
                have_both += 1
                translated = translate(cds)
                if translated == prot:
                    exact_match += 1
                else:
                    mismatch += 1
                    if len(translated) == len(prot):
                        len_mismatch_only += 1
                    if len(mismatch_examples) < 10:
                        mismatch_examples.append((
                            row.get("uniprot_accession"), row.get("ensembl_transcript_used"),
                            len(prot), len(translated), prot[:30], translated[:30]
                        ))

            if rna and cds:
                rna_and_cds_present += 1
                if cds not in rna:
                    rna_missing_cds_substr += 1

    print(f"Total rows: {total}")
    print()
    print("=== CDS structural sanity ===")
    print(f"rows with non-empty cds_sequence checked above")
    print(f"  starts with ATG: {starts_atg}")
    print(f"  length divisible by 3: {div3}")
    print(f"  ends with stop codon: {ends_stop}")
    print(f"  contains non-ACGTN characters: {non_acgt}")
    print()
    print("=== CDS -> protein translation check ===")
    print(f"rows with both cds_sequence and protein sequence: {have_both}")
    print(f"  translated CDS EXACTLY matches table's protein sequence: {exact_match}")
    print(f"  MISMATCH: {mismatch}  (same length but different residues: {len_mismatch_only})")
    print()
    if mismatch_examples:
        print("Mismatch examples (uniprot, transcript, protlen, translen, prot_start, translated_start):")
        for ex in mismatch_examples:
            print(" ", ex)
    print()
    print("=== rna_sequence contains cds_sequence check ===")
    print(f"rows with both rna_sequence and cds_sequence: {rna_and_cds_present}")
    print(f"  cds_sequence NOT found as substring of rna_sequence: {rna_missing_cds_substr}")


if __name__ == "__main__":
    main()
