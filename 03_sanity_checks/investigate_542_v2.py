import csv
import sys
import json
import time
import urllib.request
from collections import Counter

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

# Known human selenoprotein genes (the ~25-member human selenoproteome)
KNOWN_SELENOPROTEINS = {
    "GPX1","GPX2","GPX3","GPX4","GPX6","TXNRD1","TXNRD2","TXNRD3",
    "DIO1","DIO2","DIO3","SELENOP","SELENOK","SELENOS","SELENOW","SEPHS2",
    "SELENOF","SELENOH","SELENOI","SELENOM","SELENON","SELENOO","SELENOT",
    "SELENOV","MSRB1","SELENBP1",
}


def translate(cds):
    prot = []
    for i in range(0, len(cds) - len(cds) % 3, 3):
        codon = cds[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        prot.append(aa)
    return ''.join(prot)


def n_genes(field):
    return len(set(g.split(".")[0] for g in field.split(";") if g.strip()))


def http_get(url, timeout=60):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def main():
    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    single_gene_mismatches = []
    for row in rows:
        cds = row.get("cds_sequence", "")
        prot = row.get("sequence", "")
        if not cds or not prot:
            continue
        translated = translate(cds)
        if translated == prot or translated == "X" + prot:
            continue
        if n_genes(row.get("ensembl_gene", "")) > 1:
            continue
        single_gene_mismatches.append((row.get("uniprot_accession"), cds, prot, translated))

    print(f"Total single-gene mismatches: {len(single_gene_mismatches)}")
    accs = [x[0] for x in single_gene_mismatches]

    # fetch real gene names/descriptions from UniProt in bulk
    gene_info = {}
    batch_size = 100
    for i in range(0, len(accs), batch_size):
        batch = accs[i:i+batch_size]
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions="
               + ",".join(batch) + "&fields=accession,gene_names,protein_name,cc_ptm,ft_non_std")
        for attempt in range(4):
            try:
                body = http_get(url)
                break
            except Exception as e:
                print(f"  batch error at {i}: {e}", file=sys.stderr)
                time.sleep(5)
                body = None
        if not body:
            continue
        d = json.loads(body)
        for entry in d.get("results", []):
            acc = entry["primaryAccession"]
            genes = [g.get("geneName", {}).get("value", "") for g in entry.get("genes", [])]
            pdesc = entry.get("proteinDescription", {})
            name = (pdesc.get("recommendedName", {}).get("fullName", {}).get("value")
                    or (pdesc.get("submissionNames") or [{}])[0].get("fullName", {}).get("value") or "")
            has_selenocysteine_feature = any(
                f.get("description", "").lower().find("selenocysteine") >= 0 or f.get("type") == "Non-standard residue"
                for f in entry.get("features", [])
            )
            gene_info[acc] = (";".join(genes), name, has_selenocysteine_feature)
        print(f"  ... {i+len(batch)}/{len(accs)} done", flush=True)
        time.sleep(0.3)

    cats = Counter()
    examples = {}
    known_seleno_hits = []

    for acc, cds, prot, translated in single_gene_mismatches:
        genes, name, has_sec_feature = gene_info.get(acc, ("", "", False))
        gene_upper = genes.upper()
        desc = (genes + " " + name).lower()

        is_known_seleno = any(g in KNOWN_SELENOPROTEINS for g in gene_upper.split(";"))
        if is_known_seleno or has_sec_feature or "selenoprotein" in desc:
            cat = "Selenoprotein (confirmed via UniProt annotation)"
            known_seleno_hits.append((acc, genes, len(prot), len(translated)))
        elif any(k in desc for k in ["immunoglobulin", "t cell receptor", " ig ", "kappa constant", "lambda constant"]):
            cat = "Immunoglobulin/TCR"
        elif "mitochondri" in desc:
            cat = "Mitochondrial"
        elif "readthrough" in desc:
            cat = "Annotated readthrough gene"
        elif "pseudogene" in desc or "putative" in desc:
            cat = "Pseudogene/putative"
        elif not genes and not name:
            cat = "No UniProt annotation returned"
        else:
            cat = "Other/unclassified"
        cats[cat] += 1
        if len(examples.setdefault(cat, [])) < 8:
            examples[cat].append((acc, genes, name[:50], len(prot), len(translated)))

    print()
    print("=== Functional categorization (real UniProt gene names/descriptions) ===")
    for cat, count in cats.most_common():
        print(f"{cat}: {count} ({100*count/len(single_gene_mismatches):.1f}%)")
        for ex in examples[cat]:
            print(f"    {ex}")
    print()
    print(f"Confirmed selenoproteins found: {len(known_seleno_hits)}")
    for h in known_seleno_hits:
        print(" ", h)


if __name__ == "__main__":
    main()
