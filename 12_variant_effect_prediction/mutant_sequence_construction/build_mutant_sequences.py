import csv
import json
import sys
from collections import Counter

csv.field_size_limit(sys.maxsize)

# usage: build_mutant_sequences.py <wt_json> <changes_tsv> <out_mutant_fasta> <out_wt_fasta> <out_map_tsv>
WT_FILE = sys.argv[1]
CHANGES_FILE = sys.argv[2]
OUT_FASTA = sys.argv[3]
OUT_WT_FASTA = sys.argv[4]
OUT_MAP = sys.argv[5]

AA3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U",
}

HANDLEABLE = {"missense", "nonsense", "del", "ins", "delins", "dup"}


def build_mutant(seq, category, details):
    """Returns mutant sequence string, or None if inconsistent with WT."""
    if category == "missense":
        pos, wt, mut = details["pos"], details["wt"], details["mut"]
        wt1, mut1 = AA3_TO_1.get(wt), AA3_TO_1.get(mut)
        if wt1 is None or mut1 is None or pos > len(seq) or seq[pos - 1] != wt1:
            return None
        return seq[:pos - 1] + mut1 + seq[pos:]
    if category == "nonsense":
        keep_to = details["keep_to"]
        if keep_to < 0 or keep_to > len(seq):
            return None
        return seq[:keep_to]
    if category == "del":
        start, end = details["start"], details["end"]
        if start < 1 or end > len(seq) or start > end:
            return None
        return seq[:start - 1] + seq[end:]
    if category == "ins":
        after = details["after"]
        insert = "".join(AA3_TO_1.get(c, "") for c in details["insert"])
        if len(insert) != len(details["insert"]) or after > len(seq):
            return None
        new_seq = seq[:after] + insert + seq[after:]
        if details.get("truncate"):
            new_seq = seq[:after] + insert
        return new_seq
    if category == "delins":
        start, end = details["start"], details["end"]
        insert = "".join(AA3_TO_1.get(c, "") for c in details["ins"])
        if len(insert) != len(details["ins"]) or start < 1 or end > len(seq) or start > end:
            return None
        if details.get("truncate"):
            return seq[:start - 1] + insert
        return seq[:start - 1] + insert + seq[end:]
    if category == "dup":
        start, end = details["start"], details["end"]
        if start < 1 or end > len(seq) or start > end:
            return None
        segment = seq[start - 1:end]
        return seq[:end] + segment + seq[end:]
    return None


with open(WT_FILE) as f:
    wt_seqs = json.load(f)

print(f"Loaded {len(wt_seqs)} WT sequences")

# unique WT sequences actually referenced -> write once
counts = Counter()
mismatch_examples = []

written_wt = set()

with open(CHANGES_FILE, newline="") as f, \
     open(OUT_FASTA, "w") as fasta_out, \
     open(OUT_WT_FASTA, "w") as wt_fasta_out, \
     open(OUT_MAP, "w", newline="") as map_out:

    map_writer = csv.writer(map_out, delimiter="\t")
    map_writer.writerow(["VariationID", "GeneID", "GeneSymbol", "accession", "category", "mutant_id"])

    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        category = row["category"]
        if category not in HANDLEABLE:
            counts[f"skip_{category}"] += 1
            continue
        acc = row["accession"]
        entry = wt_seqs.get(acc)
        if entry is None:
            counts["skip_no_wt_seq"] += 1
            continue
        seq = entry["sequence"]

        if acc not in written_wt:
            wt_fasta_out.write(f">{acc}\n{seq}\n")
            written_wt.add(acc)

        details = json.loads(row["details"])
        mut_seq = build_mutant(seq, category, details)
        if mut_seq is None:
            counts["mismatch_or_invalid"] += 1
            if len(mismatch_examples) < 10:
                mismatch_examples.append((row["VariationID"], acc, category, details))
            continue
        if len(mut_seq) == 0:
            counts["empty_mutant"] += 1
            continue

        mutant_id = f"{row['VariationID']}_{acc}"
        fasta_out.write(f">{mutant_id}\n{mut_seq}\n")
        map_writer.writerow([row["VariationID"], row["GeneID"], row["GeneSymbol"], acc, category, mutant_id])
        counts[f"built_{category}"] += 1

print("\nCounts:")
for k, v in counts.most_common():
    print(f"  {k}: {v}")

print(f"\nUnique WT sequences written: {len(written_wt)}")
print(f"\nMismatch examples (first 10):")
for ex in mismatch_examples:
    print(f"  {ex}")

print(f"\nWrote {OUT_FASTA}")
print(f"Wrote {OUT_WT_FASTA}")
print(f"Wrote {OUT_MAP}")
