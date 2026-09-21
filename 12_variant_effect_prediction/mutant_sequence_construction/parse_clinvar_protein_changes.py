import csv
import sys
import re
import json
from collections import Counter

csv.field_size_limit(sys.maxsize)

SRC = "/u/project/kappel/ddcohn/ClinVar_variant_summary_complete.csv"
OUT_PARSED = "/u/project/kappel/ddcohn/SpliceAI/../protein_variant_effects/clinvar_protein_changes.tsv"
OUT_ACCESSIONS = "/u/project/kappel/ddcohn/protein_variant_effects/unique_transcript_accessions.txt"

import os
os.makedirs("/u/project/kappel/ddcohn/protein_variant_effects", exist_ok=True)
OUT_PARSED = "/u/project/kappel/ddcohn/protein_variant_effects/clinvar_protein_changes.tsv"

AA3 = {
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
    "Ter", "Sec", "Xaa",
}

# transcript accession + protein notation, e.g.
# NM_014630.3(ZNF592):c.3136G>A (p.Gly1046Arg)
NAME_RE = re.compile(r"^(?P<acc>[A-Z]{2}_\d+\.\d+)\([^)]*\):.*?\(p\.(?P<pchange>[^)]+)\)")

SUB_RE = re.compile(r"^([A-Za-z]{3})(\d+)([A-Za-z]{3})$")
FS_RE = re.compile(r"fs")
DELINS_RE = re.compile(r"^([A-Za-z]{3})(\d+)_([A-Za-z]{3})(\d+)delins([A-Za-z]+)$")
INS_TER_RE = re.compile(r"^([A-Za-z]{3})(\d+)_([A-Za-z]{3})(\d+)insTer$")
DEL_RANGE_RE = re.compile(r"^([A-Za-z]{3})(\d+)_([A-Za-z]{3})(\d+)del$")
DEL_SINGLE_RE = re.compile(r"^([A-Za-z]{3})(\d+)del$")
DUP_RANGE_RE = re.compile(r"^([A-Za-z]{3})(\d+)_([A-Za-z]{3})(\d+)dup$")
DUP_SINGLE_RE = re.compile(r"^([A-Za-z]{3})(\d+)dup$")
INS_RANGE_RE = re.compile(r"^([A-Za-z]{3})(\d+)_([A-Za-z]{3})(\d+)ins([A-Za-z]+)$")
DELINS_SINGLE_RE = re.compile(r"^([A-Za-z]{3})(\d+)delins([A-Za-z]+)$")


def classify(pchange):
    if pchange == "=" or pchange.endswith("="):
        return "synonymous_explicit", {}
    if FS_RE.search(pchange):
        return "frameshift", {}
    # order matters: check del/dup/delins/insTer (which contain 3-letter
    # keywords like "del"/"dup"/"ins" that could be mis-parsed by SUB_RE)
    # before the generic single-substitution pattern.
    m = DELINS_RE.match(pchange)
    if m:
        wt1, pos1, wt2, pos2, ins = m.group(1), int(m.group(2)), m.group(3), int(m.group(4)), m.group(5)
        ins_codons = [ins[i:i+3] for i in range(0, len(ins), 3)]
        truncate = bool(ins_codons) and ins_codons[-1] == "Ter"
        if truncate:
            ins_codons = ins_codons[:-1]
        if not all(c in AA3 for c in ins_codons):
            return "unparsed", {}
        return "delins", {"start": pos1, "end": pos2, "wt_start": wt1, "wt_end": wt2, "ins": ins_codons, "truncate": truncate}
    m = DELINS_SINGLE_RE.match(pchange)
    if m:
        wt, pos, ins = m.group(1), int(m.group(2)), m.group(3)
        ins_codons = [ins[i:i+3] for i in range(0, len(ins), 3)]
        truncate = bool(ins_codons) and ins_codons[-1] == "Ter"
        if truncate:
            ins_codons = ins_codons[:-1]
        if not all(c in AA3 for c in ins_codons):
            return "unparsed", {}
        return "delins", {"start": pos, "end": pos, "wt_start": wt, "wt_end": wt, "ins": ins_codons, "truncate": truncate}
    m = INS_TER_RE.match(pchange)
    if m:
        wt1, pos1, wt2, pos2 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
        # e.g. p.Thr33_Trp34insTer: stop inserted immediately after pos1,
        # equivalent to truncating the protein right after residue pos1.
        return "nonsense", {"keep_to": pos1}
    m = INS_RANGE_RE.match(pchange)
    if m:
        wt1, pos1, wt2, pos2, ins = m.group(1), int(m.group(2)), m.group(3), int(m.group(4)), m.group(5)
        ins_codons = [ins[i:i+3] for i in range(0, len(ins), 3)]
        truncate = bool(ins_codons) and ins_codons[-1] == "Ter"
        if truncate:
            ins_codons = ins_codons[:-1]
        if not all(c in AA3 for c in ins_codons):
            return "unparsed", {}
        return "ins", {"after": pos1, "insert": ins_codons, "truncate": truncate}
    m = DEL_RANGE_RE.match(pchange)
    if m:
        wt1, pos1, wt2, pos2 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
        return "del", {"start": pos1, "end": pos2, "wt_start": wt1, "wt_end": wt2}
    m = DEL_SINGLE_RE.match(pchange)
    if m:
        wt, pos = m.group(1), int(m.group(2))
        return "del", {"start": pos, "end": pos, "wt_start": wt, "wt_end": wt}
    m = DUP_RANGE_RE.match(pchange)
    if m:
        wt1, pos1, wt2, pos2 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
        return "dup", {"start": pos1, "end": pos2, "wt_start": wt1, "wt_end": wt2}
    m = DUP_SINGLE_RE.match(pchange)
    if m:
        wt, pos = m.group(1), int(m.group(2))
        return "dup", {"start": pos, "end": pos, "wt_start": wt, "wt_end": wt}
    m = SUB_RE.match(pchange)
    if m:
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        if wt not in AA3 or mut not in AA3:
            return "unparsed", {}
        if wt == mut:
            return "synonymous", {"pos": pos, "wt": wt, "mut": mut}
        if mut == "Ter":
            return "nonsense", {"keep_to": pos - 1}
        if wt == "Ter":
            return "stoploss", {"pos": pos, "mut": mut}
        return "missense", {"pos": pos, "wt": wt, "mut": mut}
    return "unparsed", {}


counts = Counter()
accessions = set()

with open(SRC, newline="", encoding="utf-8") as f, \
     open(OUT_PARSED, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["VariationID", "GeneID", "GeneSymbol", "accession", "category", "details"])
    for row in csv.DictReader(f):
        if row["CoordinateAssembly"] != "GRCh38":
            continue
        name = row["Name"]
        m = NAME_RE.match(name)
        if not m:
            counts["no_protein_notation"] += 1
            continue
        acc = m.group("acc")
        pchange = m.group("pchange")
        category, details = classify(pchange)
        counts[category] += 1
        if category in ("missense", "nonsense", "del", "delins", "dup", "ins"):
            accessions.add(acc)
        writer.writerow([row["VariationID"], row["GeneID"], row["GeneSymbol"], acc, category, json.dumps(details)])

print("Category counts:")
for cat, n in counts.most_common():
    print(f"  {cat}: {n}")
print(f"\nUnique transcript accessions needed (missense/nonsense/del/delins/dup): {len(accessions)}")

with open(OUT_ACCESSIONS, "w") as f:
    for acc in sorted(accessions):
        f.write(acc + "\n")

print(f"\nWrote {OUT_PARSED}")
print(f"Wrote {OUT_ACCESSIONS}")
