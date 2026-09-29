import urllib.request
import urllib.parse
import time
import json
import os

D = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"

with open(f"{D}/rbp_gene_symbols.txt") as f:
    ours = {l.strip() for l in f if l.strip()}
with open(f"{D}/protein_ids.txt") as f:
    theirs = {l.strip() for l in f if l.strip()}

new_symbols = sorted(theirs - ours)
print(f"Ours: {len(ours)}, theirs: {len(theirs)}, to resolve: {len(new_symbols)}")

OUT_FILE = f"{D}/combine_resolve.json"
results = {}
if os.path.exists(OUT_FILE):
    results = json.load(open(OUT_FILE))
    print(f"Resuming: {len(results)} already resolved")


def search_gene(symbol):
    query = f"gene:{symbol} AND organism_id:9606 AND reviewed:true"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": query, "format": "fasta", "size": 3}
    )
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                text = r.read().decode("utf-8", errors="replace")
                return text if text.strip() else None
        except Exception as e:
            print(f"  {symbol}: attempt {attempt+1} failed: {e}")
            time.sleep(2)
    return None


def parse_first(text):
    acc = gene = None
    seq = []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc is not None:
                break
            parts = line[1:].split("|")
            acc = parts[1] if len(parts) >= 2 else line[1:].split()[0]
            for tok in line.split():
                if tok.startswith("GN="):
                    gene = tok[3:]
        else:
            seq.append(line.strip())
    return acc, gene, "".join(seq)


remaining = [s for s in new_symbols if s not in results]
print(f"Remaining to search: {len(remaining)}")
for i, sym in enumerate(remaining):
    text = search_gene(sym)
    if text:
        acc, gene, seq = parse_first(text)
        results[sym] = {"accession": acc, "uniprot_gene": gene, "sequence": seq} if seq else None
    else:
        results[sym] = None
    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{len(remaining)} done")
        json.dump(results, open(OUT_FILE, "w"))
    time.sleep(0.25)
json.dump(results, open(OUT_FILE, "w"))

already_covered = 0
genuinely_new = 0
unresolved = 0
new_entries = {}
for sym, rec in results.items():
    if not rec:
        unresolved += 1
        continue
    gn = rec["uniprot_gene"]
    if gn in ours or gn == sym and sym in ours:
        already_covered += 1
    elif gn in ours:
        already_covered += 1
    else:
        genuinely_new += 1
        new_entries[gn or sym] = rec

print(f"\nAlready covered under current gene name (alias resolved to existing entry): {already_covered}")
print(f"Genuinely new genes to add: {genuinely_new}")
print(f"Unresolved (no UniProt hit at all): {unresolved}")

with open(f"{D}/combine_unresolved.txt", "w") as f:
    for sym, rec in results.items():
        if not rec:
            f.write(sym + "\n")
print(f"Wrote {D}/combine_unresolved.txt ({unresolved} symbols)")

with open(f"{D}/combine_new_entries.json", "w") as f:
    json.dump(new_entries, f)
print(f"Wrote {D}/combine_new_entries.json ({len(new_entries)} entries)")
