"""
Joins this project's RNP-machine and RNA-class flags onto a labmate's
(asharma) much larger merged table (PTM predictions, SLiMs, SPURS/RSA
structural data, richer IDR biophysics, HPA localization, condensate
predictions, 137 columns, 2,332 RBPs) -- confirmed to otherwise already
cover everything in rbp_wt_list_with_go.tsv except these 12 base columns
(Ribosomal_protein, Spliceosomal_protein, Other_RNP_machine, and the 9
RNA-class binding flags), each paired with its own `_terms` companion
for 22 new columns total. Writes to the user's own RBP directory, not
asharma's -- the source table is read-only, not modified in place.
"""
import csv

csv.field_size_limit(2**31 - 1)

ASHARMA = "/u/project/kappel/asharma/RBP/merged/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa.csv"
OUR_TABLE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/rbp_wt_list_with_go.tsv"
OUT = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"

ADD_COLS = [
    "Ribosomal_protein", "Spliceosomal_protein", "Other_RNP_machine", "Other_RNP_machine_terms",
    "mRNA_binding", "mRNA_binding_terms",
    "tRNA_binding", "tRNA_binding_terms",
    "rRNA_binding", "rRNA_binding_terms",
    "snRNA_binding", "snRNA_binding_terms",
    "snoRNA_binding", "snoRNA_binding_terms",
    "piRNA_binding", "piRNA_binding_terms",
    "miRNA_binding", "miRNA_binding_terms",
    "siRNA_binding", "siRNA_binding_terms",
    "lncRNA_binding", "lncRNA_binding_terms",
]

our_data = {}
with open(OUR_TABLE, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        our_data[row["UniProtID"]] = {c: row[c] for c in ADD_COLS}

print(f"Our table: {len(our_data)} RBPs loaded")

with open(ASHARMA, newline="") as f:
    r = csv.DictReader(f)
    orig_fieldnames = r.fieldnames
    assert not any(c in orig_fieldnames for c in ADD_COLS), "unexpected pre-existing column collision"
    fieldnames = orig_fieldnames + ADD_COLS
    rows = list(r)

n_matched = 0
for row in rows:
    acc = row["uniprot_accession"]
    extra = our_data.get(acc)
    if extra:
        row.update(extra)
        n_matched += 1
    else:
        for c in ADD_COLS:
            row[c] = ""

print(f"Asharma table: {len(rows)} rows total")
print(f"Matched to our RNP/RNA-class data: {n_matched}")
print(f"Not matched (left blank): {len(rows) - n_matched}")

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

print(f"\nWrote {OUT}")
