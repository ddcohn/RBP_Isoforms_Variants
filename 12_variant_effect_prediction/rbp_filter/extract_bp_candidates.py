import csv, re
csv.field_size_limit(2**31 - 1)

PATH = "/private/tmp/claude-501/-Users-danielcohn/16426b22-0d42-460f-bc55-5a06e4b1aee2/scratchpad/rbp_wt_list_with_go.tsv"

CLASSES = ["mRNA_binding", "tRNA_binding", "rRNA_binding", "snRNA_binding",
           "snoRNA_binding", "piRNA_binding", "miRNA_binding", "siRNA_binding",
           "lncRNA_binding"]
KEYWORDS = {
    "mRNA": re.compile(r"\bmRNA\b", re.I),
    "tRNA": re.compile(r"\btRNA\b", re.I),
    "rRNA": re.compile(r"\brRNA\b", re.I),
    "snRNA": re.compile(r"\bsnRNA\b", re.I),
    "snoRNA": re.compile(r"\bsnoRNA\b", re.I),
    "piRNA": re.compile(r"\bpiRNA\b", re.I),
    "miRNA": re.compile(r"\bmiRNA\b|\bmicroRNA\b", re.I),
    "siRNA": re.compile(r"\bsiRNA\b", re.I),
    "lncRNA": re.compile(r"\blncRNA\b|\blong non-?coding RNA\b", re.I),
}

n_total = 0
n_unclassified = []  # proteins with general RNA binding but no specific class

bp_terms_by_keyword = {k: {} for k in KEYWORDS}  # keyword -> {goid:term: count}
family_by_keyword = {k: {} for k in KEYWORDS}

with open(PATH, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        n_total += 1
        has_general = "GO:0003723" in row["GO_Molecular_Function"]
        any_class = any(row[c].strip() == "Y" for c in CLASSES)
        is_gap = has_general and not any_class

        bp = row["GO_Biological_Process"]
        fam = row["HGNC_Gene_Family"]

        for kw, pat in KEYWORDS.items():
            for entry in bp.split(";"):
                if pat.search(entry):
                    bp_terms_by_keyword[kw][entry] = bp_terms_by_keyword[kw].get(entry, 0) + 1
            if fam and pat.search(fam):
                for f_entry in fam.split(";") if ";" in fam else [fam]:
                    if pat.search(f_entry):
                        family_by_keyword[kw][f_entry] = family_by_keyword[kw].get(f_entry, 0) + 1

        if is_gap:
            n_unclassified.append(row)

print(f"Total: {n_total}, gap (general RNA binding, no specific class): {len(n_unclassified)}")
print()

for kw in KEYWORDS:
    print(f"=== BP terms matching '{kw}' (top by frequency) ===")
    for term, cnt in sorted(bp_terms_by_keyword[kw].items(), key=lambda kv: -kv[1])[:15]:
        print(f"  {cnt:5d}  {term}")
    print()

for kw in KEYWORDS:
    if family_by_keyword[kw]:
        print(f"=== HGNC families matching '{kw}' ===")
        for fam, cnt in sorted(family_by_keyword[kw].items(), key=lambda kv: -kv[1]):
            print(f"  {cnt:5d}  {fam}")
        print()
