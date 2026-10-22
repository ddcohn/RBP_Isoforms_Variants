"""
Single streaming pass over satmut_all_predictions.tsv (31,525,484 rows)
that joins each mutation to its protein's WT baseline (TANGO, protGPS,
DeepLoc) and accumulates three views of "how much did this mutation
change the prediction vs WT":

  1. By wt_aa (20 categories): which original residue, when mutated,
     causes the largest average change
  2. By mut_aa (20 categories): which substituted residue causes the
     largest average change
  3. By (wt_aa, mut_aa) pair (20x20): full substitution matrix, for a
     heatmap -- which specific substitutions matter most
  4. By protein (uniprot_accession): mean change across all of a
     protein's mutants, to rank which proteins are most susceptible to
     mutation-driven change

Change metrics (composite, one per tool):
  - tango_abs_delta_agg   = |mutant_AGG - wt_AGG|
  - protgps_total_abs_delta = sum(|mutant_c - wt_c| for c in 12 compartments)
  - deeploc_total_abs_delta = sum(|mutant_c - wt_c| for c in 13 numeric compartment/structure scores)
  - deeploc_loc_changed   = 1 if mutant's Localizations string != WT's, else 0

mutation_id format: uniprot_accession|gene_symbol|wt_aa+position+mut_aa
-- wt_aa/mut_aa are parsed directly from the ID string (single-letter
amino acid codes, position is 1+ digits in between).

Run via submit_analyze_satmut_deltas.sh on a compute node (not login
node -- this file is 15.8GB and a single Python pass over it has
silently died on the login node before).
"""
import csv
import re
import sys

csv.field_size_limit(sys.maxsize)

BASE = "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only"
SATMUT_ALL = f"{BASE}/satmut_all_predictions.tsv"
WT_TANGO = f"{BASE}/wt_tango_scores.tsv"
MERGED_TABLE = "/u/project/kappel/ddcohn/RBP/merged_260903B_ptm_pred_slim_RBPs_domains_spurs_rsa_iso_hpa_with_rnp_rna_class.csv"

OUT_WT_AA = f"{BASE}/satmut_change_by_wt_aa.tsv"
OUT_MUT_AA = f"{BASE}/satmut_change_by_mut_aa.tsv"
OUT_PAIR = f"{BASE}/satmut_change_by_substitution_pair.tsv"
OUT_PROTEIN = f"{BASE}/satmut_change_by_protein.tsv"

PROTGPS_COLS = [
    "NUCLEAR_SPECKLE_Score", "P-BODY_Score", "PML-BDOY_Score",
    "POST_SYNAPTIC_DENSITY_Score", "STRESS_GRANULE_Score", "CHROMOSOME_Score",
    "NUCLEOLUS_Score", "NUCLEAR_PORE_COMPLEX_Score", "CAJAL_BODY_Score",
    "RNA_GRANULE_Score", "CELL_JUNCTION_Score", "TRANSCRIPTIONAL_Score",
]
DEEPLOC_NUMERIC_COLS = [
    "Cytoplasm", "Nucleus", "Extracellular", "Cell membrane", "Mitochondrion",
    "Plastid", "Endoplasmic reticulum", "Lysosome/Vacuole", "Golgi apparatus",
    "Peroxisome", "Peripheral", "Transmembrane", "Lipid anchor", "Soluble",
]

MUT_ID_RE = re.compile(r"^(?P<acc>[^|]+)\|(?P<gene>[^|]+)\|(?P<wt>[A-Z])(?P<pos>\d+)(?P<mut>[A-Z])$")


def load_wt_tango():
    wt = {}
    with open(WT_TANGO, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                wt[row["uniprot_accession"]] = float(row["AGG"])
            except ValueError:
                pass
    print(f"Loaded {len(wt)} WT TANGO baselines")
    return wt


def load_wt_protgps_deeploc():
    protgps_wt = {}
    deeploc_wt = {}
    gene_symbol_of = {}
    seen = set()
    with open(MERGED_TABLE, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row["uniprot_accession"]
            if acc in seen:
                continue
            seen.add(acc)
            gene_symbol_of[acc] = row.get("gene_symbol", "")

            pvals = []
            ok = True
            for c in PROTGPS_COLS:
                v = row.get(f"protgps_{c}", "").strip()
                if not v or v == "NA":
                    ok = False
                    break
                pvals.append(float(v))
            if ok:
                protgps_wt[acc] = pvals

            dvals = []
            ok = True
            for c in DEEPLOC_NUMERIC_COLS:
                v = row.get(f"deeploc_{c}", "").strip()
                if not v or v == "NA":
                    ok = False
                    break
                dvals.append(float(v))
            loc = row.get("deeploc_Localizations", "").strip()
            if ok and loc:
                deeploc_wt[acc] = (dvals, loc)

    print(f"Loaded {len(protgps_wt)} WT protGPS baselines, {len(deeploc_wt)} WT DeepLoc baselines")
    return protgps_wt, deeploc_wt, gene_symbol_of


class Accumulator:
    def __init__(self):
        self.n = 0
        self.sum_tango = 0.0
        self.n_tango = 0
        self.sum_protgps = 0.0
        self.n_protgps = 0
        self.sum_deeploc = 0.0
        self.n_deeploc = 0
        self.n_loc_changed = 0
        self.n_loc_total = 0

    def add(self, tango_d, protgps_d, deeploc_d, loc_changed):
        self.n += 1
        if tango_d is not None:
            self.sum_tango += tango_d
            self.n_tango += 1
        if protgps_d is not None:
            self.sum_protgps += protgps_d
            self.n_protgps += 1
        if deeploc_d is not None:
            self.sum_deeploc += deeploc_d
            self.n_deeploc += 1
        if loc_changed is not None:
            self.n_loc_total += 1
            if loc_changed:
                self.n_loc_changed += 1

    def row(self, key):
        return [
            *key,
            self.n,
            self.sum_tango / self.n_tango if self.n_tango else "",
            self.sum_protgps / self.n_protgps if self.n_protgps else "",
            self.sum_deeploc / self.n_deeploc if self.n_deeploc else "",
            self.n_loc_changed / self.n_loc_total if self.n_loc_total else "",
        ]


def main():
    wt_tango = load_wt_tango()
    wt_protgps, wt_deeploc, gene_symbol_of = load_wt_protgps_deeploc()

    by_wt_aa = {}
    by_mut_aa = {}
    by_pair = {}
    by_protein = {}

    n_total = 0
    n_parsed = 0

    with open(SATMUT_ALL, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            n_total += 1
            m = MUT_ID_RE.match(row["mutation_id"])
            if not m:
                continue
            n_parsed += 1
            acc = m.group("acc")
            wt_aa = m.group("wt")
            mut_aa = m.group("mut")

            tango_d = None
            wt_agg = wt_tango.get(acc)
            if wt_agg is not None:
                try:
                    tango_d = abs(float(row["AGG"]) - wt_agg)
                except ValueError:
                    pass

            protgps_d = None
            wt_p = wt_protgps.get(acc)
            if wt_p is not None:
                try:
                    protgps_d = sum(
                        abs(float(row[c]) - w) for c, w in zip(PROTGPS_COLS, wt_p)
                    )
                except ValueError:
                    pass

            deeploc_d = None
            loc_changed = None
            wt_dl = wt_deeploc.get(acc)
            if wt_dl is not None:
                wt_dvals, wt_loc = wt_dl
                try:
                    deeploc_d = sum(
                        abs(float(row[c]) - w) for c, w in zip(DEEPLOC_NUMERIC_COLS, wt_dvals)
                    )
                except ValueError:
                    pass
                loc_changed = row["Localizations"].strip() != wt_loc

            for d, key in ((by_wt_aa, (wt_aa,)), (by_mut_aa, (mut_aa,)), (by_pair, (wt_aa, mut_aa))):
                acc_obj = d.setdefault(key, Accumulator())
                acc_obj.add(tango_d, protgps_d, deeploc_d, loc_changed)

            prot_acc = by_protein.setdefault(acc, Accumulator())
            prot_acc.add(tango_d, protgps_d, deeploc_d, loc_changed)

    print(f"Total rows: {n_total}")
    print(f"Parsed mutation_id successfully: {n_parsed}")

    header = ["n_mutations", "mean_tango_abs_delta_agg", "mean_protgps_total_abs_delta",
              "mean_deeploc_total_abs_delta", "frac_deeploc_localization_changed"]

    with open(OUT_WT_AA, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["wt_aa"] + header)
        for key, acc_obj in sorted(by_wt_aa.items()):
            w.writerow(acc_obj.row(key))

    with open(OUT_MUT_AA, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["mut_aa"] + header)
        for key, acc_obj in sorted(by_mut_aa.items()):
            w.writerow(acc_obj.row(key))

    with open(OUT_PAIR, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["wt_aa", "mut_aa"] + header)
        for key, acc_obj in sorted(by_pair.items()):
            w.writerow(acc_obj.row(key))

    with open(OUT_PROTEIN, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["uniprot_accession", "gene_symbol"] + header)
        for acc, acc_obj in by_protein.items():
            w.writerow([acc, gene_symbol_of.get(acc, "")] + acc_obj.row(()))

    print(f"Wrote {OUT_WT_AA}")
    print(f"Wrote {OUT_MUT_AA}")
    print(f"Wrote {OUT_PAIR}")
    print(f"Wrote {OUT_PROTEIN}")


if __name__ == "__main__":
    main()
