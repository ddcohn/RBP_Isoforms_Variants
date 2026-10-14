"""
RNA target class flags, derived from GO_Molecular_Function. Supersedes
add_trna_rrna_flags.py's tRNA_binding/rRNA_binding columns with a wider
set of real GO IDs per class (not just the single base term), plus adds
mRNA, snRNA, snoRNA, piRNA, miRNA, siRNA, and lncRNA as explicit flags --
mRNA previously had no flag of its own (implicit default for "neither
tRNA nor rRNA"); it now does, since GO:0003729 (mRNA binding) exists and
was unused.

Every GO ID below was found by grepping this table's own
GO_Molecular_Function column for "RNA binding" and manually sorting the
real results into classes -- not assumed from memory. Left out:
structural/motif terms that don't identify an RNA class (double-stranded,
single-stranded, G-quadruplex, poly(U)/poly(C), regulatory region,
generic "translation factor activity, RNA binding"), and terms for a
single specific ncRNA gene already covered by Other_RNP_machine
(ribonuclease P RNA binding, 7S/SRP RNA binding) or not cleanly
assignable to one class (telomerase RNA binding, telomeric repeat RNA
binding -- TERC is its own thing, not folded into lncRNA here).
"""
import csv

D = "/Users/danielcohn/Desktop/RBP_repo_work/RBP_seq_STRING_OpenTargets/12_variant_effect_prediction/rbp_filter"
TABLE = f"{D}/rbp_wt_list_with_go.tsv"

RNA_CLASS_GO_IDS = {
    "mRNA_binding": ["GO:0003729", "GO:0036002", "GO:1990825"],
    "tRNA_binding": ["GO:0000049", "GO:1904678", "GO:1990856"],
    "rRNA_binding": ["GO:0019843", "GO:0008097", "GO:0070180", "GO:0070181", "GO:1990400"],
    "snRNA_binding": [
        "GO:0017069", "GO:0017070", "GO:0030619", "GO:0030620", "GO:0030621",
        "GO:0030622", "GO:0030623", "GO:0030624", "GO:0030626", "GO:0071209", "GO:0097322",
    ],
    "snoRNA_binding": ["GO:0030515", "GO:0034511", "GO:0034512", "GO:0034513"],
    "piRNA_binding": ["GO:0034584"],
    "miRNA_binding": ["GO:0035198", "GO:0070878", "GO:0070883"],
    "siRNA_binding": ["GO:0035197"],
    "lncRNA_binding": ["GO:0106222", "GO:0002153"],
}

with open(TABLE, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    orig_fieldnames = r.fieldnames
    rows = list(r)

# Drop any columns this script (or its superseded predecessor) already wrote,
# whether from the old flag-only layout or a prior run of this script's
# flag+terms layout, so re-running never produces duplicate columns.
stale = {"tRNA_binding", "rRNA_binding"}
for cls in RNA_CLASS_GO_IDS:
    stale.add(cls)
    stale.add(f"{cls}_terms")
fieldnames = [c for c in orig_fieldnames if c not in stale]
for cls in RNA_CLASS_GO_IDS:
    fieldnames += [cls, f"{cls}_terms"]

counts = {k: 0 for k in RNA_CLASS_GO_IDS}
for row in rows:
    mf = row.get("GO_Molecular_Function", "")
    mf_entries = [e for e in mf.split(";") if e]
    for cls, go_ids in RNA_CLASS_GO_IDS.items():
        matched = [e for e in mf_entries if any(e.startswith(go_id) for go_id in go_ids)]
        row[cls] = "Y" if matched else "N"
        row[f"{cls}_terms"] = ";".join(matched)
        if matched:
            counts[cls] += 1

print(f"Total genes: {len(rows)}")
for cls, n in counts.items():
    print(f"{cls}: {n}")

with open(TABLE, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
print(f"\nUpdated {TABLE} in place")
