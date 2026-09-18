import csv
import sys
import ast

csv.field_size_limit(sys.maxsize)

PATH = "/u/project/kappel/RBP/Misc/Classical_RBDs_Annotated.csv"

seen = set()
n_with_domains = 0
n_multi_instance = 0  # proteins with >=2 total domain instances (any family)
n_overlap_within_family = 0
n_overlap_across_family = 0
examples_within = []
examples_across = []

with open(PATH, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    for row in reader:
        uid = row.get("uniprot_id", "")
        if uid in seen:
            continue
        seen.add(uid)
        raw = row.get("Domains_range", "")
        if not raw or raw == "{}":
            continue
        try:
            rangedict = ast.literal_eval(raw)
        except (ValueError, SyntaxError):
            continue
        if not isinstance(rangedict, dict) or not rangedict:
            continue
        n_with_domains += 1

        all_instances = []  # (family, start, end)
        for fam, ranges in rangedict.items():
            for r in ranges:
                if isinstance(r, (list, tuple)) and len(r) == 2:
                    all_instances.append((fam, r[0], r[1]))

        if len(all_instances) < 2:
            continue
        n_multi_instance += 1

        # within-family overlap
        for fam, ranges in rangedict.items():
            if len(ranges) < 2:
                continue
            srt = sorted(ranges, key=lambda r: r[0])
            for i in range(len(srt) - 1):
                if srt[i][1] > srt[i+1][0]:
                    n_overlap_within_family += 1
                    if len(examples_within) < 5:
                        examples_within.append((uid, fam, ranges))
                    break

        # cross-family overlap (any two instances from different families overlapping)
        srt_all = sorted(all_instances, key=lambda x: x[1])
        found_cross = False
        for i in range(len(srt_all) - 1):
            fam1, s1, e1 = srt_all[i]
            fam2, s2, e2 = srt_all[i+1]
            if e1 > s2 and fam1 != fam2:
                found_cross = True
                if len(examples_across) < 5:
                    examples_across.append((uid, rangedict))
                break
        if found_cross:
            n_overlap_across_family += 1

print(f"Distinct proteins: {len(seen)}")
print(f"Proteins with >=1 domain instance: {n_with_domains}")
print(f"Proteins with >=2 total domain instances: {n_multi_instance}")
print(f"Proteins with WITHIN-family overlapping ranges: {n_overlap_within_family}")
print(f"Proteins with CROSS-family overlapping ranges: {n_overlap_across_family}")

print("\n--- within-family overlap examples ---")
for uid, fam, ranges in examples_within:
    print(f"  {uid}  family={fam}  ranges={ranges}")

print("\n--- cross-family overlap examples ---")
for uid, rangedict in examples_across:
    print(f"  {uid}  ranges={rangedict}")
