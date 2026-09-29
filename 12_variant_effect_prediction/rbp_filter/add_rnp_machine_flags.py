import csv

D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"

RIBOSOME_KEYWORDS = ["ribosome", "ribosomal subunit"]
SPLICEOSOME_KEYWORDS = ["spliceosom", "snrnp"]
OTHER_RNP_KEYWORDS = [
    "exosome (rnase complex)",
    "signal recognition particle",
    "rnase p complex",
    "rnase mrp complex",
    "telomerase holoenzyme complex",
    "telomerase complex",
    "vault ribonucleoprotein complex",
    "box c/d snorna",
    "box c/d sno",
    "box h/aca sno",
    "editosome",
    "nuclear pore complex",
    "commitment complex",
    "sm-like protein family complex",
]


def classify(go_cc):
    cc_lower = go_cc.lower()
    is_ribo = any(k in cc_lower for k in RIBOSOME_KEYWORDS)
    is_splice = any(k in cc_lower for k in SPLICEOSOME_KEYWORDS)
    other_hits = [k for k in OTHER_RNP_KEYWORDS if k in cc_lower]
    return is_ribo, is_splice, other_hits


rows_out = []
n_ribo = n_splice = n_other = 0
with open(f"{D}/rbp_wt_list_with_go.tsv", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    fieldnames = r.fieldnames + ["Ribosomal_protein", "Spliceosomal_protein", "Other_RNP_machine", "Other_RNP_machine_terms"]
    for row in r:
        is_ribo, is_splice, other_hits = classify(row.get("GO_Cellular_Component", ""))
        row["Ribosomal_protein"] = "Y" if is_ribo else "N"
        row["Spliceosomal_protein"] = "Y" if is_splice else "N"
        row["Other_RNP_machine"] = "Y" if other_hits else "N"
        row["Other_RNP_machine_terms"] = ";".join(other_hits)
        if is_ribo:
            n_ribo += 1
        if is_splice:
            n_splice += 1
        if other_hits:
            n_other += 1
        rows_out.append(row)

print(f"Total genes: {len(rows_out)}")
print(f"Ribosomal protein: {n_ribo}")
print(f"Spliceosomal protein: {n_splice}")
print(f"Other named large RNP machine: {n_other}")

with open(f"{D}/rbp_wt_list_with_go.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows_out)
print(f"\nUpdated {D}/rbp_wt_list_with_go.tsv in place")
