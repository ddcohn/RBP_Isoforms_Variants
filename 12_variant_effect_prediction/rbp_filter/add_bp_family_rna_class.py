"""
Extends the existing 9 RNA-class flags (currently MF-only, GO:0003723
descendant terms) by also checking GO_Biological_Process and
HGNC_Gene_Family for proteins that have the general GO:0003723 "RNA
binding" term but didn't match any specific class via MF alone.

All BP GO IDs and HGNC family names below were pulled from this
table's own real data (extract_bp_candidates.py), not guessed from
memory, then manually curated: only terms that imply direct binding/
processing of that RNA type are kept. Explicitly excluded:
  - "positive/negative regulation of X transcription by RNA polymerase
    Y" terms -- these regulate DNA-template transcription of the RNA
    gene, not binding of the mature RNA product.
  - HGNC family "MicroRNA protein coding host genes" (159 hits) -- a
    real false positive caught during curation: this means the gene's
    own intron happens to host a miRNA gene, not that the protein
    product binds miRNA.
"""
import csv
csv.field_size_limit(2**31 - 1)

PATH = "/private/tmp/claude-501/-Users-danielcohn/16426b22-0d42-460f-bc55-5a06e4b1aee2/scratchpad/rbp_wt_list_with_go.tsv"
OUT = "/private/tmp/claude-501/-Users-danielcohn/16426b22-0d42-460f-bc55-5a06e4b1aee2/scratchpad/rbp_wt_list_with_go_v2.tsv"

CLASSES = ["mRNA_binding", "tRNA_binding", "rRNA_binding", "snRNA_binding",
           "snoRNA_binding", "piRNA_binding", "miRNA_binding", "siRNA_binding",
           "lncRNA_binding"]

BP_IDS = {
    "mRNA_binding": [
        "GO:0000398", "GO:0006397", "GO:0000381", "GO:0006406", "GO:0051028",
        "GO:0000184", "GO:0048025", "GO:0006402", "GO:0031124", "GO:0048255",
        "GO:0043488", "GO:0000956", "GO:0061157", "GO:0048024", "GO:0000380",
        "GO:0016556",
    ],
    "tRNA_binding": [
        "GO:0006418", "GO:0008033", "GO:0030488", "GO:0006400", "GO:0001682",
        "GO:0006388", "GO:0042780", "GO:0071038", "GO:0031119", "GO:1990180",
        "GO:0006432", "GO:0002098", "GO:0000379", "GO:0006420",
    ],
    "rRNA_binding": [
        "GO:0006364", "GO:0030490", "GO:0000462", "GO:0000463", "GO:0000460",
        "GO:2000234", "GO:0000470", "GO:0016075", "GO:0031167", "GO:0000467",
        "GO:0000447", "GO:0000448", "GO:0071035", "GO:0000451",
    ],
    "snRNA_binding": [
        "GO:0034475", "GO:0016180", "GO:0034472", "GO:0006408", "GO:0031120",
        "GO:0034473", "GO:0034476", "GO:0016076", "GO:0071037", "GO:0061015",
        "GO:0034477", "GO:0034474", "GO:1990438",
    ],
    "snoRNA_binding": [
        "GO:0071051", "GO:0048254", "GO:0000454", "GO:0071036",
    ],
    "piRNA_binding": [
        "GO:0034587", "GO:0141196", "GO:0140991", "GO:0141006", "GO:0140965",
        "GO:0140990",
    ],
    "miRNA_binding": [
        "GO:0035196", "GO:0035278", "GO:0031054", "GO:0031053", "GO:0035195",
        "GO:0010586", "GO:0010587", "GO:2000637", "GO:0035279", "GO:2000632",
        "GO:0060964", "GO:0060965", "GO:2000633", "GO:1990428",
    ],
    "siRNA_binding": [
        "GO:0030422", "GO:1990432", "GO:0090625", "GO:0070921", "GO:0141194",
    ],
    "lncRNA_binding": [
        "GO:0000512", "GO:0110064", "GO:0180035", "GO:0180034",
    ],
}

FAMILY_KEYWORDS = {
    "mRNA_binding": ["Apolipoprotein B mRNA editing enzyme catalytic subunits",
                      "Regulation of nuclear pre-mRNA domain containing"],
    "tRNA_binding": ["Aminoacyl tRNA synthetases, Class I", "Aminoacyl tRNA synthetases, Class II",
                      "Aminoacyl tRNA synthetase complex interacting multifunctional proteins",
                      "tRNA splicing ligase complex", "tRNA-splicing endonuclease subunits"],
    # miRNA_binding: deliberately NOT adding "MicroRNA protein coding host genes" -- false positive, see docstring.
}


def main():
    with open(PATH, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        orig_fieldnames = reader.fieldnames
        rows = list(reader)

    new_cols = [f"{c}_source" for c in CLASSES]
    fieldnames = orig_fieldnames + [c for c in new_cols if c not in orig_fieldnames]

    n_gap_before = 0
    n_newly_classified = 0
    n_recovered_by_class = {c: 0 for c in CLASSES}

    for row in rows:
        has_general = "GO:0003723" in row["GO_Molecular_Function"]
        any_class_mf = any(row[c].strip() == "Y" for c in CLASSES)
        is_gap = has_general and not any_class_mf

        for c in CLASSES:
            row.setdefault(f"{c}_source", "MF" if row[c].strip() == "Y" else "")

        if not is_gap:
            continue
        n_gap_before += 1

        bp_entries = row["GO_Biological_Process"].split(";")
        bp_ids_present = {e.split(":", 2)[0] + ":" + e.split(":", 2)[1] for e in bp_entries if e.count(":") >= 2}
        fam = row["HGNC_Gene_Family"]

        got_any = False
        for c in CLASSES:
            matched_bp = [gid for gid in BP_IDS.get(c, []) if gid in bp_ids_present]
            matched_fam = [fname for fname in FAMILY_KEYWORDS.get(c, []) if fname in fam]
            if matched_bp or matched_fam:
                row[c] = "Y"
                src = []
                if matched_bp:
                    src.append("BP")
                if matched_fam:
                    src.append("Family")
                row[f"{c}_source"] = "+".join(src)
                n_recovered_by_class[c] += 1
                got_any = True

        if got_any:
            n_newly_classified += 1

    with open(OUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Gap before (general RNA binding, no specific class): {n_gap_before}")
    print(f"Newly classified via BP/Family: {n_newly_classified}")
    print(f"Still unclassified: {n_gap_before - n_newly_classified}")
    print()
    for c in CLASSES:
        print(f"  {c}: +{n_recovered_by_class[c]} newly recovered")
    print()
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
