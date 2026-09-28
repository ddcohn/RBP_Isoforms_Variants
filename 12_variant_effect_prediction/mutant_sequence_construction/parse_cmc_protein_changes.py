import csv
import sys
import gzip
import re
import json
import os
from collections import Counter

csv.field_size_limit(sys.maxsize)

# usage: parse_cmc_protein_changes.py <raw_tsv_gz> <out_parsed_tsv> <out_accessions_txt>
SRC = sys.argv[1]
OUT_PARSED = sys.argv[2]
OUT_ACCESSIONS = sys.argv[3]
os.makedirs(os.path.dirname(OUT_PARSED), exist_ok=True)

# 1-letter AA codes only (CMC's own convention) + '*' for stop
AA1 = set("ACDEFGHIKLMNPQRSTVWY")

DEFERRED_CATEGORIES = {
    "Deletion - Frameshift", "Insertion - Frameshift", "Complex - frameshift",
    "Frameshift", "Nonstop extension",
}
NO_EFFECT_CATEGORIES = {"Substitution - coding silent"}

SUB_RE = re.compile(r"^p\.([A-Z])(\d+)([A-Z*])$")
DELINS_RANGE_RE = re.compile(r"^p\.([A-Z])(\d+)_([A-Z])(\d+)delins([A-Z]+)$")
DELINS_SINGLE_RE = re.compile(r"^p\.([A-Z])(\d+)delins([A-Z]+)$")
INS_RANGE_RE = re.compile(r"^p\.([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)$")
DEL_RANGE_RE = re.compile(r"^p\.([A-Z])(\d+)_([A-Z])(\d+)del$")
DEL_SINGLE_RE = re.compile(r"^p\.([A-Z])(\d+)del$")
DUP_RANGE_RE = re.compile(r"^p\.([A-Z])(\d+)_([A-Z])(\d+)dup$")
DUP_SINGLE_RE = re.compile(r"^p\.([A-Z])(\d+)dup$")


def classify(pchange):
    m = DELINS_RANGE_RE.match(pchange)
    if m:
        wt1, pos1, wt2, pos2, ins = m.group(1), int(m.group(2)), m.group(3), int(m.group(4)), m.group(5)
        if not all(c in AA1 for c in ins):
            return "unparsed", {}
        return "delins", {"start": pos1, "end": pos2, "ins": list(ins)}
    m = DELINS_SINGLE_RE.match(pchange)
    if m:
        wt, pos, ins = m.group(1), int(m.group(2)), m.group(3)
        if not all(c in AA1 for c in ins):
            return "unparsed", {}
        return "delins", {"start": pos, "end": pos, "ins": list(ins)}
    m = INS_RANGE_RE.match(pchange)
    if m:
        wt1, pos1, wt2, pos2, ins = m.group(1), int(m.group(2)), m.group(3), int(m.group(4)), m.group(5)
        if not all(c in AA1 for c in ins):
            return "unparsed", {}
        return "ins", {"after": pos1, "insert": list(ins)}
    m = DEL_RANGE_RE.match(pchange)
    if m:
        wt1, pos1, wt2, pos2 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
        return "del", {"start": pos1, "end": pos2}
    m = DEL_SINGLE_RE.match(pchange)
    if m:
        wt, pos = m.group(1), int(m.group(2))
        return "del", {"start": pos, "end": pos}
    m = DUP_RANGE_RE.match(pchange)
    if m:
        wt1, pos1, wt2, pos2 = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
        return "dup", {"start": pos1, "end": pos2}
    m = DUP_SINGLE_RE.match(pchange)
    if m:
        wt, pos = m.group(1), int(m.group(2))
        return "dup", {"start": pos, "end": pos}
    m = SUB_RE.match(pchange)
    if m:
        wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
        if wt not in AA1:
            return "unparsed", {}
        if mut == wt:
            return "synonymous", {}
        if mut == "*":
            return "nonsense", {"keep_to": pos - 1}
        if wt == "*":
            return "stoploss", {}
        if mut not in AA1:
            return "unparsed", {}
        return "missense", {"pos": pos, "wt": wt, "mut": mut}
    return "unparsed", {}


counts = Counter()
accessions = set()

with gzip.open(SRC, "rt", encoding="utf-8", errors="replace") as f, \
     open(OUT_PARSED, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["GENOMIC_MUTATION_ID", "GENE_NAME", "accession", "category", "details"])
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        desc_aa = row.get("Mutation Description AA", "")
        pchange = row.get("Mutation AA", "")
        acc_full = row.get("ACCESSION_NUMBER", "")
        acc = acc_full.split(".")[0] if acc_full else ""

        if desc_aa in NO_EFFECT_CATEGORIES:
            counts["synonymous"] += 1
            continue
        if desc_aa in DEFERRED_CATEGORIES:
            counts["frameshift_or_stoploss"] += 1
            continue
        if not pchange or not pchange.startswith("p."):
            counts["no_protein_notation"] += 1
            continue

        category, details = classify(pchange)
        counts[category] += 1
        if category in ("missense", "nonsense", "del", "ins", "delins", "dup") and acc:
            accessions.add(acc)
            writer.writerow([row.get("GENOMIC_MUTATION_ID", ""), row.get("GENE_NAME", ""),
                              acc, category, json.dumps(details)])

print("Category counts:")
for cat, n in counts.most_common():
    print(f"  {cat}: {n}")
print(f"\nUnique transcript accessions needed: {len(accessions)}")

with open(OUT_ACCESSIONS, "w") as f:
    for acc in sorted(accessions):
        f.write(acc + "\n")

print(f"\nWrote {OUT_PARSED}")
print(f"Wrote {OUT_ACCESSIONS}")
