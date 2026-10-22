"""
Fixes two real false-positive bugs in the original add_rnp_machine_flags.py,
found by inspecting every real GO_Cellular_Component term that matched its
naive substring keywords ("ribosome", "ribosomal subunit", "spliceosom",
"snrnp") against the live data:

1. Spliceosomal_protein wrongly included GO:0005683 "U7 snRNP" (handles
   histone pre-mRNA 3'-end processing, NOT splicing) and GO:0120259
   "7SK snRNP" (regulates RNA Pol II elongation via P-TEFb, NOT
   splicing) -- both matched purely because they contain the substring
   "snRNP". Wrongly flagged 4 proteins (HEXIM1, LARP7, LSM11, MEPCE)
   that have no other spliceosome-related annotation.

2. Ribosomal_protein wrongly included the GO:0030684/0030686/0030687/
   0030688/0034388 "preribosome" family -- these are ribosome assembly
   *intermediates*, not the mature ribosome, and proteins annotated
   there are typically biogenesis factors (chaperones, GTPases, rRNA-
   processing enzymes) that are never part of the final translating
   ribosome. Wrongly flagged 42 proteins (BMS1, BOP1, BYSL, NOP14,
   PES1, etc.) that have no other ribosome-related annotation.

Fix: exact GO ID allowlists instead of substring keywords (matching
the approach used later in this project for the 9 RNA-class flags).
The excluded terms aren't just dropped -- they're relabeled into
Other_RNP_machine under their real identity (U7 snRNP, 7SK snRNP,
ribosome biogenesis/preribosome) so the information isn't lost, just
no longer mislabeled.
"""
import csv

D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"

RIBOSOME_GO_IDS = {
    "GO:0005840",  # ribosome
    "GO:0015934",  # large ribosomal subunit
    "GO:0015935",  # small ribosomal subunit
    "GO:0022625",  # cytosolic large ribosomal subunit
    "GO:0022626",  # cytosolic ribosome
    "GO:0022627",  # cytosolic small ribosomal subunit
    "GO:0005761",  # mitochondrial ribosome
    "GO:0005762",  # mitochondrial large ribosomal subunit
    "GO:0005763",  # mitochondrial small ribosomal subunit
}

SPLICEOSOME_GO_IDS = {
    "GO:0005681",  # spliceosomal complex
    "GO:0005682",  # U5 snRNP
    "GO:0005684",  # U2-type spliceosomal complex
    "GO:0005685",  # U1 snRNP
    "GO:0005686",  # U2 snRNP
    "GO:0005687",  # U4 snRNP
    "GO:0005688",  # U6 snRNP
    "GO:0005689",  # U12-type spliceosomal complex
    "GO:0005690",  # U4atac snRNP
    "GO:0034693",  # U11/U12 snRNP
    "GO:0044530",  # supraspliceosomal complex
    "GO:0046540",  # U4/U6 x U5 tri-snRNP complex
    "GO:0071001",  # U4/U6 snRNP
    "GO:0071004",  # U2-type prespliceosome
    "GO:0071005",  # U2-type precatalytic spliceosome
    "GO:0071006",  # U2-type catalytic step 1 spliceosome
    "GO:0071007",  # U2-type catalytic step 2 spliceosome
    "GO:0071008",  # U2-type post-mRNA release spliceosomal complex
    "GO:0071011",  # precatalytic spliceosome
    "GO:0071013",  # catalytic step 2 spliceosome
    "GO:0071014",  # post-mRNA release spliceosomal complex
    "GO:0071020",  # post-spliceosomal complex
    "GO:0071021",  # U2-type post-spliceosomal complex
    "GO:0071254",  # cytoplasmic U snRNP body
    "GO:0097525",  # spliceosomal snRNP complex
    "GO:0097526",  # spliceosomal tri-snRNP complex
}

# Explicitly NOT spliceosomal despite containing "snRNP" -- confirmed false
# positives in the original keyword-based version.
EXCLUDED_SPLICEOSOME_GO_IDS = {
    "GO:0005683": "U7 snRNP (histone mRNA 3'-end processing, not splicing)",
    "GO:0120259": "7SK snRNP (RNA Pol II elongation control, not splicing)",
}

# Explicitly NOT the mature ribosome -- assembly intermediates / biogenesis.
EXCLUDED_RIBOSOME_GO_IDS = {
    "GO:0030684": "preribosome (biogenesis intermediate)",
    "GO:0030686": "90S preribosome (biogenesis intermediate)",
    "GO:0030687": "preribosome, large subunit precursor (biogenesis intermediate)",
    "GO:0030688": "preribosome, small subunit precursor (biogenesis intermediate)",
    "GO:0034388": "Pwp2p-containing subcomplex of 90S preribosome (biogenesis intermediate)",
}

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

# Confirmed false positive: the SRP RECEPTOR complex (SRPRA/SRPRB, a membrane
# GTPase pair that docks SRP at the ER) is not itself an RNA-binding particle
# -- only the real SRP (GO:0048500 / GO:0005786) contains the 7SL RNA.
# Matched via the "signal recognition particle" substring otherwise.
EXCLUDED_OTHER_RNP_GO_IDS = {
    "GO:0005785": "signal recognition particle receptor complex (not itself RNA-binding; docks SRP at the ER membrane)",
}


def parse_go_entries(go_cc):
    """Returns list of (go_id, full_entry_text) for entries like 'GO:0005840:ribosome'."""
    out = []
    for entry in go_cc.split(";"):
        entry = entry.strip()
        if entry.count(":") >= 2:
            go_id = entry.split(":", 2)[0] + ":" + entry.split(":", 2)[1]
            out.append((go_id, entry))
    return out


RIBOSOME_GO_LABELS = {
    "GO:0005840": "ribosome", "GO:0015934": "large ribosomal subunit",
    "GO:0015935": "small ribosomal subunit", "GO:0022625": "cytosolic large ribosomal subunit",
    "GO:0022626": "cytosolic ribosome", "GO:0022627": "cytosolic small ribosomal subunit",
    "GO:0005761": "mitochondrial ribosome", "GO:0005762": "mitochondrial large ribosomal subunit",
    "GO:0005763": "mitochondrial small ribosomal subunit",
}
SPLICEOSOME_GO_LABELS = {
    "GO:0005681": "spliceosomal complex", "GO:0005682": "U5 snRNP",
    "GO:0005684": "U2-type spliceosomal complex", "GO:0005685": "U1 snRNP",
    "GO:0005686": "U2 snRNP", "GO:0005687": "U4 snRNP", "GO:0005688": "U6 snRNP",
    "GO:0005689": "U12-type spliceosomal complex", "GO:0005690": "U4atac snRNP",
    "GO:0034693": "U11/U12 snRNP", "GO:0044530": "supraspliceosomal complex",
    "GO:0046540": "U4/U6 x U5 tri-snRNP complex", "GO:0071001": "U4/U6 snRNP",
    "GO:0071004": "U2-type prespliceosome", "GO:0071005": "U2-type precatalytic spliceosome",
    "GO:0071006": "U2-type catalytic step 1 spliceosome", "GO:0071007": "U2-type catalytic step 2 spliceosome",
    "GO:0071008": "U2-type post-mRNA release spliceosomal complex", "GO:0071011": "precatalytic spliceosome",
    "GO:0071013": "catalytic step 2 spliceosome", "GO:0071014": "post-mRNA release spliceosomal complex",
    "GO:0071020": "post-spliceosomal complex", "GO:0071021": "U2-type post-spliceosomal complex",
    "GO:0071254": "cytoplasmic U snRNP body", "GO:0097525": "spliceosomal snRNP complex",
    "GO:0097526": "spliceosomal tri-snRNP complex",
}


def classify(go_cc):
    entries = parse_go_entries(go_cc)
    go_ids = {gid for gid, _ in entries}

    ribo_matched_ids = go_ids & RIBOSOME_GO_IDS
    splice_matched_ids = go_ids & SPLICEOSOME_GO_IDS
    is_ribo = bool(ribo_matched_ids)
    is_splice = bool(splice_matched_ids)
    ribo_terms = ";".join(f"{gid}:{RIBOSOME_GO_LABELS[gid]}" for gid in sorted(ribo_matched_ids))
    splice_terms = ";".join(f"{gid}:{SPLICEOSOME_GO_LABELS[gid]}" for gid in sorted(splice_matched_ids))

    other_hits = []
    # per-entry keyword check, skipping entries whose GO ID is a confirmed
    # false positive -- whole-string substring matching would still fire on
    # the excluded entry's own text (e.g. "signal recognition particle
    # receptor complex" contains "signal recognition particle").
    for gid, entry in entries:
        if gid in EXCLUDED_OTHER_RNP_GO_IDS:
            continue
        el = entry.lower()
        for kw in OTHER_RNP_KEYWORDS:
            if kw in el and kw not in other_hits:
                other_hits.append(kw)

    # U7 snRNP / 7SK snRNP / preribosome ARE genuine RNA-protein complexes --
    # just the wrong specific category -- so they're surfaced in
    # Other_RNP_machine with a corrected label, not silently dropped.
    for gid, label in EXCLUDED_SPLICEOSOME_GO_IDS.items():
        if gid in go_ids:
            other_hits.append(label)
    for gid, label in EXCLUDED_RIBOSOME_GO_IDS.items():
        if gid in go_ids:
            other_hits.append(label)
    # The SRP RECEPTOR complex is NOT an RNP at all (no RNA component --
    # it's a pure protein GTPase pair that docks SRP at the ER membrane),
    # so it's excluded entirely rather than relabeled into Other_RNP_machine.

    return is_ribo, is_splice, other_hits, ribo_terms, splice_terms


rows_out = []
n_ribo = n_splice = n_other = 0
n_ribo_fixed = n_splice_fixed = 0

AUDIT_ROWS = []

with open(f"{D}/rbp_wt_list_with_go.tsv", newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    fieldnames = [fn for fn in r.fieldnames if fn not in
                  ("Ribosomal_protein", "Spliceosomal_protein", "Other_RNP_machine", "Other_RNP_machine_terms")]
    fieldnames += ["Ribosomal_protein", "Ribosomal_protein_terms", "Spliceosomal_protein",
                   "Spliceosomal_protein_terms", "Other_RNP_machine", "Other_RNP_machine_terms"]
    for row in r:
        was_ribo = row.get("Ribosomal_protein") == "Y"
        was_splice = row.get("Spliceosomal_protein") == "Y"

        is_ribo, is_splice, other_hits, ribo_terms, splice_terms = classify(row.get("GO_Cellular_Component", ""))
        row["Ribosomal_protein"] = "Y" if is_ribo else "N"
        row["Ribosomal_protein_terms"] = ribo_terms
        row["Spliceosomal_protein"] = "Y" if is_splice else "N"
        row["Spliceosomal_protein_terms"] = splice_terms
        row["Other_RNP_machine"] = "Y" if other_hits else "N"
        row["Other_RNP_machine_terms"] = ";".join(other_hits)

        if is_ribo or is_splice or other_hits:
            AUDIT_ROWS.append([row["GeneSymbol"], row["UniProtID"], row["Ribosomal_protein"],
                                ribo_terms, row["Spliceosomal_protein"], splice_terms,
                                row["Other_RNP_machine"], ";".join(other_hits)])

        if is_ribo:
            n_ribo += 1
        if is_splice:
            n_splice += 1
        if other_hits:
            n_other += 1
        if was_ribo and not is_ribo:
            n_ribo_fixed += 1
        if was_splice and not is_splice:
            n_splice_fixed += 1
        rows_out.append(row)

print(f"Total genes: {len(rows_out)}")
print(f"Ribosomal protein: {n_ribo} ({n_ribo_fixed} corrected from false-positive Y)")
print(f"Spliceosomal protein: {n_splice} ({n_splice_fixed} corrected from false-positive Y)")
print(f"Other named large RNP machine: {n_other}")

with open(f"{D}/rbp_wt_list_with_go.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows_out)
print(f"\nUpdated {D}/rbp_wt_list_with_go.tsv in place")

with open(f"{D}/rnp_machine_flags_audit.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["GeneSymbol", "UniProtID", "Ribosomal_protein", "Ribosomal_protein_terms",
                "Spliceosomal_protein", "Spliceosomal_protein_terms",
                "Other_RNP_machine", "Other_RNP_machine_terms"])
    w.writerows(AUDIT_ROWS)
print(f"Wrote audit file ({len(AUDIT_ROWS)} flagged proteins): {D}/rnp_machine_flags_audit.tsv")
