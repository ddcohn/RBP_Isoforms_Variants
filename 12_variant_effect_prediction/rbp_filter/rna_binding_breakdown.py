import csv
csv.field_size_limit(2**31 - 1)

PATH = "/private/tmp/claude-501/-Users-danielcohn/16426b22-0d42-460f-bc55-5a06e4b1aee2/scratchpad/rbp_wt_list_with_go.tsv"

CLASSES = ["mRNA_binding", "tRNA_binding", "rRNA_binding", "snRNA_binding",
           "snoRNA_binding", "piRNA_binding", "miRNA_binding", "siRNA_binding",
           "lncRNA_binding"]

n_total = 0
n_general_rna_binding = 0  # GO:0003723 "RNA binding" (general MF term)
class_counts = {c: 0 for c in CLASSES}
n_any_class = 0
n_no_class_but_general = 0

with open(PATH, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        n_total += 1
        mf = row["GO_Molecular_Function"]
        has_general = "GO:0003723" in mf
        if has_general:
            n_general_rna_binding += 1

        any_class = False
        for c in CLASSES:
            if row[c].strip() == "Y":
                class_counts[c] += 1
                any_class = True
        if any_class:
            n_any_class += 1
        if has_general and not any_class:
            n_no_class_but_general += 1

print(f"Total WT RBPs: {n_total}")
print(f"Has general GO:0003723 'RNA binding': {n_general_rna_binding} ({100*n_general_rna_binding/n_total:.1f}%)")
print(f"Has at least one of the 9 specific RNA-class flags: {n_any_class} ({100*n_any_class/n_total:.1f}%)")
print(f"Has GO:0003723 but none of the 9 specific classes matched: {n_no_class_but_general}")
print()
for c in CLASSES:
    print(f"{c}: {class_counts[c]}")
