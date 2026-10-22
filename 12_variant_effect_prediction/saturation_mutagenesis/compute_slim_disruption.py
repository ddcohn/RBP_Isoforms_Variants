"""
For every ELM SLiM match already annotated on each protein's WT sequence
(asharma's SLiM_range/SLiM_class/SLiM_discrete_seq columns in the merged
table), enumerate every possible point mutation that falls inside that
match's span, and test whether the mutant residue breaks the ELM class's
own regex pattern (full match on the matched substring with one residue
substituted).

Only mutations that fall inside an existing SLiM match are considered --
a mutation outside every annotated match trivially can't disrupt a motif,
so it's left out of the output table rather than writing 31M "N/A" rows.

mutation_id format matches the TANGO/protGPS/DeepLoc satmut tables
exactly (uniprot_accession|gene_symbol|wt_aa+position+mut_aa), so this
can be joined onto those later by mutation_id.

ELM class regexes downloaded from http://elm.eu.org/elms/elms_index.tsv
(note: elm.eu.org serves plain HTTP only -- https connects to nothing on
port 443, confirmed via curl -v; use http:// for any future ELM fetch).
"""
import csv
import json
import re

csv.field_size_limit(2**31 - 1)

MERGED_TABLE = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"
ELM_CLASSES_TSV = "elm_classes.tsv"  # checked into this directory; refresh via: curl http://elm.eu.org/elms/elms_index.tsv
OUT_FILE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_slim_disruption.tsv"

AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")

# --- load ELM class regexes ---
elm_regex = {}
with open(ELM_CLASSES_TSV, newline="") as f:
    lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t", quotechar='"')
    for row in reader:
        ident = row["ELMIdentifier"]
        pattern = row["Regex"]
        try:
            elm_regex[ident] = re.compile(pattern)
        except re.error as e:
            print(f"  skipping uncompilable regex for {ident}: {pattern} ({e})")

print(f"Loaded {len(elm_regex)} ELM class regexes")

# --- load per-protein SLiM annotations (dedupe by uniprot_accession;
#     SLiM columns are protein-level, broadcast across ENCORI duplicate rows) ---
seen_proteins = set()
protein_slims = {}  # uniprot_accession -> (gene_symbol, [(start, end, elm_id, discrete_seq), ...])

with open(MERGED_TABLE, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        acc = row["uniprot_accession"]
        if acc in seen_proteins:
            continue
        seen_proteins.add(acc)

        count = row.get("SLiM_count", "").strip()
        if not count or count == "0":
            continue

        try:
            ranges = json.loads(row["SLiM_range"])
            classes = json.loads(row["SLiM_class"])
            seqs = json.loads(row["SLiM_discrete_seq"])
        except (json.JSONDecodeError, KeyError):
            continue

        if not (len(ranges) == len(classes) == len(seqs)):
            continue

        slims = []
        for (start, end), elm_id, dseq in zip(ranges, classes, seqs):
            slims.append((start, end, elm_id, dseq))
        protein_slims[acc] = (row["gene_symbol"], slims)

print(f"Proteins with >=1 SLiM match: {len(protein_slims)}")
total_matches = sum(len(v[1]) for v in protein_slims.values())
print(f"Total SLiM match instances across all proteins: {total_matches}")

# --- enumerate every point mutation inside each SLiM match, test disruption ---
n_rows = 0
n_skipped_seq_mismatch = 0
n_skipped_unknown_class = 0

with open(OUT_FILE, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "mutation_id", "uniprot_accession", "gene_symbol", "position",
        "wt_aa", "mut_aa", "elm_identifier", "slim_start", "slim_end",
        "wt_discrete_seq", "mut_discrete_seq", "disrupts_motif",
    ])

    for acc, (gene_symbol, slims) in protein_slims.items():
        for start, end, elm_id, dseq in slims:
            regex = elm_regex.get(elm_id)
            if regex is None:
                n_skipped_unknown_class += 1
                continue
            if len(dseq) != (end - start + 1):
                n_skipped_seq_mismatch += 1
                continue
            if not regex.fullmatch(dseq):
                # sanity check: asharma's own match should already satisfy
                # the current regex. If not (e.g. ELM updated a pattern
                # since asharma's run), skip rather than report a bogus call.
                n_skipped_seq_mismatch += 1
                continue

            for offset, wt_aa in enumerate(dseq):
                position = start + offset
                for mut_aa in AA_ALPHABET:
                    if mut_aa == wt_aa:
                        continue
                    mut_dseq = dseq[:offset] + mut_aa + dseq[offset + 1:]
                    disrupts = "N" if regex.fullmatch(mut_dseq) else "Y"
                    mutation_id = f"{acc}|{gene_symbol}|{wt_aa}{position}{mut_aa}"
                    writer.writerow([
                        mutation_id, acc, gene_symbol, position, wt_aa, mut_aa,
                        elm_id, start, end, dseq, mut_dseq, disrupts,
                    ])
                    n_rows += 1

print(f"\nRows written: {n_rows}")
print(f"Skipped SLiM instances (seq/regex mismatch): {n_skipped_seq_mismatch}")
print(f"Skipped SLiM instances (unknown ELM class): {n_skipped_unknown_class}")
print(f"Wrote {OUT_FILE}")
