import csv
import sys
from collections import Counter

path = "/u/project/kappel/RBP/Isoform_Table/Isoform_Post_Merge_PSLab_OpenTargets_Updated_Interim_20260817.csv"
csv.field_size_limit(sys.maxsize)

# Checklist items -> best-guess matching column(s) in the actual CSV.
# None-mapped fields flagged separately as "no matching column found".
checklist_map = {
    "Gene ID": ["gene_symbol", "ncbi_gene_id", "hgnc_ids"],
    "UniProt ID": ["uniprot_id"],
    "Ensembl ID (ENSG)": ["ENSG", "ensembl_gene_ids"],
    "Sequence": ["sequence"],
    "length": ["length_aa"],
    "dominant": ["dominant_isoform"],
    "Domain (predicted)": ["Domains", "InterPro_domains"],
    "Condensate Formation (cd-code)": ["Condensate Name", "UID", "Condensate Type"],
    "GO Cellular Component": ["C_ids", "C_descriptions"],
    "GO Biological Process": ["P_ids", "P_descriptions"],
    "GO Molecular Function": ["F_ids", "F_descriptions"],
    "Protein-Protein Interactions": ["PPI_ENSP_Partners", "PPI_UniProt_Partners", "string_partners_ensp_by_query"],
    "molecular weight": ["molecular_weight"],
    "isoelectric_point_whole": ["isoelectric_point"],
    "IDRs (predicted, metapredict)": ["IDR_count", "IDR_range", "idr_method"],
    "kappa": ["kappa"],
    "delta": ["delta"],
    "deltaMax": ["deltaMax"],
    "IDR_count": ["IDR_count"],
    "IDR_avg_size": ["IDR_avg_size"],
    "IDR_total_size": ["IDR_total_size"],
    "IDR_range": ["IDR_range"],
    "IDR_discrete_seq": ["IDR_discrete_seq"],
    "IDR_isoelectric": ["IDR_isoelectric_point"],
    "fraction_negative": ["fraction_negative"],
    "fraction_positive": ["fraction_positive"],
    "fraction_expanding": ["fraction_expanding"],
    "amino_acid_fractions": ["amino_acid_fractions"],
    "fraction_disorder_promoting": ["fraction_disorder_promoting"],
    "mean_net_charge": ["mean_net_charge"],
    "mean_hydropathy": ["mean_hydropathy"],
    "uversky_hydropathy": ["uversky_hydropathy"],
    "PPII_propensity": ["PPII_propensity"],
    "Domains": ["Domains"],
    "Domains_count": ["Domains_count"],
    "Domains_avg_size": ["Domains_avg_size"],
    "Domains_total_size": ["Domains_total_size"],
    "Domains_range": ["Domains_range"],
    "Domains_discrete_seq": ["Domains_discrete_seq"],
    "canonical RBDs": [],   # no obvious matching column
    "ENSP_clean": ["ENSP_clean"],
    "PPI_ENSP_Partners": ["PPI_ENSP_Partners"],
    "PPI_UniProt_Partners": ["PPI_UniProt_Partners"],
    "PPI_ENSP_Partners_in_Dataframe": ["PPI_ENSP_Partners_in_Dataframe"],
    "PPI_UniProt_Partners_in_Dataframe": ["PPI_UniProt_Partners_in_Dataframe"],
    "Condensate Name": ["Condensate Name"],
    "UID": ["UID"],
    "Condensate Type": ["Condensate Type"],
    "Proteins (condensate)": ["Proteins"],
    "DNA (condensate)": ["DNA"],
    "RNA (condensate)": ["RNA"],
    "C-mods": ["C-mods"],
    "Condensatopathy": ["Condensatopathy"],
    "Confidence Score (condensate)": ["Confidence Score"],
    "diseaseId (OpenTargets, granular)": [],
    "datatypeId (OpenTargets, granular)": [],
    "score (OpenTargets, granular)": [],
    "evidenceCount (OpenTargets, granular)": [],
    "tissues (OpenTargets, granular)": [],
    "Delta G [kT]": ["Delta G [kT]"],
    "Saturation concentration [mg/mL]": ["Saturation concentration [mg/mL]"],
    "Saturation concentration [uM]": ["Saturation concentration [uM]"],
}

all_needed_cols = set()
for cols in checklist_map.values():
    all_needed_cols.update(cols)

nonempty = Counter()
total_rows = 0
fieldnames = None

with open(path, newline='', encoding='utf-8', errors='replace') as f:
    reader = csv.DictReader(f)
    fieldnames = set(reader.fieldnames)
    for row in reader:
        total_rows += 1
        for col in all_needed_cols:
            if col in fieldnames:
                v = (row.get(col) or '').strip()
                if v and v.lower() not in ('nan', 'none', 'null', '[]', '{}'):
                    nonempty[col] += 1

print(f"Total rows: {total_rows}\n")
print(f"{'Checklist item':45s} {'Matched column(s)':55s} {'Status'}")
print("-" * 140)

for item, cols in checklist_map.items():
    present_cols = [c for c in cols if c in fieldnames]
    missing_cols = [c for c in cols if c not in fieldnames]
    if not cols:
        print(f"{item:45s} {'(none)':55s} NO MATCHING COLUMN IN TABLE")
        continue
    if not present_cols:
        print(f"{item:45s} {','.join(cols):55s} COLUMN(S) NOT FOUND: {missing_cols}")
        continue
    parts = []
    for c in present_cols:
        pct = 100.0 * nonempty[c] / total_rows if total_rows else 0
        parts.append(f"{c}={pct:.1f}%")
    status = ", ".join(parts)
    flag = ""
    worst = min(100.0 * nonempty[c] / total_rows for c in present_cols)
    if worst < 5:
        flag = "  <-- ESSENTIALLY EMPTY"
    elif worst < 50:
        flag = "  <-- SPARSE / SUSPECT"
    print(f"{item:45s} {','.join(present_cols):55s} {status}{flag}")
