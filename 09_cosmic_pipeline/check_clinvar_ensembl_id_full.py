import csv
import sys
import time
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

PATH = "/u/project/kappel/tchhabri/clinvar/fullvariantsoct25_annotated_full_fixedRBP_withSASA.csv"

n = 0
n_ensembl = 0
n_geneid = 0
n_symbol = 0
n_all_three = 0

symbol_to_ensg = defaultdict(set)
ensg_to_symbol = defaultdict(set)
symbol_to_geneid = defaultdict(set)

start = time.time()
with open(PATH, newline="", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    for row in reader:
        n += 1
        ens = (row.get("ENSEMBL_ID") or "").strip()
        gid = (row.get("GeneID") or "").strip()
        sym = (row.get("GeneSymbol") or "").strip()
        if ens:
            n_ensembl += 1
        if gid:
            n_geneid += 1
        if sym:
            n_symbol += 1
        if ens and gid and sym:
            n_all_three += 1
        if sym and ens:
            symbol_to_ensg[sym].add(ens)
            ensg_to_symbol[ens].add(sym)
        if sym and gid:
            symbol_to_geneid[sym].add(gid)
        if n % 2_000_000 == 0:
            elapsed = time.time() - start
            print(f"  ... {n:,} rows processed (elapsed {elapsed/60:.1f}min, {elapsed/n*1e6:.2f}us/row)", flush=True)

elapsed = time.time() - start
print(f"\nDone. Total rows: {n:,} (elapsed {elapsed/60:.1f}min)")
print(f"ENSEMBL_ID populated: {n_ensembl:,} ({100*n_ensembl/n:.2f}%)")
print(f"GeneID populated:     {n_geneid:,} ({100*n_geneid/n:.2f}%)")
print(f"GeneSymbol populated: {n_symbol:,} ({100*n_symbol/n:.2f}%)")
print(f"All three populated:  {n_all_three:,} ({100*n_all_three/n:.2f}%)")

inconsistent_symbol_to_ensg = {k: v for k, v in symbol_to_ensg.items() if len(v) > 1}
inconsistent_ensg_to_symbol = {k: v for k, v in ensg_to_symbol.items() if len(v) > 1}
inconsistent_symbol_to_geneid = {k: v for k, v in symbol_to_geneid.items() if len(v) > 1}

print(f"\nDistinct GeneSymbols with ENSEMBL_ID: {len(symbol_to_ensg):,}")
print(f"  ...mapping to >1 distinct ENSEMBL_ID: {len(inconsistent_symbol_to_ensg):,}")
for k, v in list(inconsistent_symbol_to_ensg.items())[:15]:
    print(f"    {k} -> {v}")

print(f"\nDistinct ENSEMBL_IDs: {len(ensg_to_symbol):,}")
print(f"  ...mapping to >1 distinct GeneSymbol: {len(inconsistent_ensg_to_symbol):,}")
for k, v in list(inconsistent_ensg_to_symbol.items())[:15]:
    print(f"    {k} -> {v}")

print(f"\nDistinct GeneSymbols with GeneID: {len(symbol_to_geneid):,}")
print(f"  ...mapping to >1 distinct GeneID: {len(inconsistent_symbol_to_geneid):,}")
for k, v in list(inconsistent_symbol_to_geneid.items())[:15]:
    print(f"    {k} -> {v}")

# genes that have GeneID/GeneSymbol but NEVER have ENSEMBL_ID anywhere (the real gap to fix)
symbols_with_geneid_only = set(symbol_to_geneid.keys()) - set(symbol_to_ensg.keys())
print(f"\nDistinct GeneSymbols that appear with GeneID but NEVER with a populated ENSEMBL_ID: {len(symbols_with_geneid_only):,}")
print("Sample:", list(symbols_with_geneid_only)[:20])
