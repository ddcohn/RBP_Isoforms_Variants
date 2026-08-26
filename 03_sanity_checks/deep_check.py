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
            rows.append(row)

    # --- 1. dig into A0A087WUL8 specifically ---
    print("=== A0A087WUL8 deep dive ===")
    for row in rows:
        if row.get("uniprot_accession") == "A0A087WUL8":
            print("gene_symbol equivalent fields not in this table; showing raw fields:")
            print("ensembl_gene:", row.get("ensembl_gene"))
            print("ensembl_protein:", row.get("ensembl_protein"))
            print("ensembl_transcript_used:", row.get("ensembl_transcript_used"))
            print("transcript_selection_method:", row.get("transcript_selection_method"))
            print("protein length (sequence col):", len(row.get("sequence", "")))
            print("protein (first 80):", row.get("sequence", "")[:80])
            print("cds_sequence length:", len(row.get("cds_sequence", "")))
            print("cds_sequence (first 80):", row.get("cds_sequence", "")[:80])
            print("translated cds (first 80):", translate(row.get("cds_sequence",""))[:80])
    print()

    # --- 2. categorize all mismatches ---
    print("=== Mismatch categorization ===")
    x_prefix = 0
    x_prefix_atg_fail = 0
    x_prefix_div3_fail = 0
    x_prefix_stop_fail = 0
    other_mismatch = []
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
            if not cds.startswith("ATG"):
                x_prefix_atg_fail += 1
            if len(cds) % 3 != 0:
                x_prefix_div3_fail += 1
            if cds[-3:] not in ("TAA","TAG","TGA"):
                x_prefix_stop_fail += 1
        else:
            other_mismatch.append((row.get("uniprot_accession"), row.get("ensembl_transcript_used"),
                                     len(prot), len(translated)))

    print(f"Total mismatches: {total_mismatch}")
    print(f"  'X' + exact protein prefix pattern: {x_prefix}")
    print(f"    of these, CDS doesn't start ATG: {x_prefix_atg_fail}")
    print(f"    of these, CDS length not div3:   {x_prefix_div3_fail}")
    print(f"    of these, CDS doesn't end stop:  {x_prefix_stop_fail}")
    print(f"  other/unexplained mismatches: {len(other_mismatch)}")
    for ex in other_mismatch[:20]:
        print("   ", ex)
    print()

    # cross-check: how many structural failures are OUTSIDE the x_prefix mismatch group entirely
    struct_fail_rows = []
    for row in rows:
        cds = row.get("cds_sequence", "")
        if not cds:
            continue
        bad = (not cds.startswith("ATG")) or (len(cds) % 3 != 0) or (cds[-3:] not in ("TAA","TAG","TGA"))
        if bad:
            prot = row.get("sequence","")
            translated = translate(cds) if prot else None
            is_x_prefix = (prot and translated == "X"+prot)
            is_exact = (prot and translated == prot)
            struct_fail_rows.append((row.get("uniprot_accession"), is_x_prefix, is_exact))

    print(f"Total rows with a structural CDS issue (bad start/frame/stop): {len(struct_fail_rows)}")
    n_x = sum(1 for _,x,_ in struct_fail_rows if x)
    n_exact = sum(1 for _,_,e in struct_fail_rows if e)
    n_other = len(struct_fail_rows) - n_x - n_exact
    print(f"  of these: X-prefix pattern={n_x}, still exact match={n_exact}, other/unexplained={n_other}")
    print()

    # --- 3. rows where cds not substring of rna ---
    print("=== rna_sequence missing cds_sequence as substring ===")
    for row in rows:
        rna = row.get("rna_sequence","")
        cds = row.get("cds_sequence","")
        if rna and cds and cds not in rna:
            print(row.get("uniprot_accession"), "| transcript:", row.get("ensembl_transcript_used"),
                  "| cds_len:", len(cds), "| rna_len:", len(rna),
                  "| method:", row.get("transcript_selection_method"))


if __name__ == "__main__":
    main()
