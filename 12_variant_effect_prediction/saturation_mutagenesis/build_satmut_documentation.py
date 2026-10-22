import openpyxl
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

OUT_PATH = "/Users/danielcohn/Desktop/261016_Satmut_Documentation.xlsx"
TOTAL_ROWS = 31525484

FONT = Font(name="Arial", size=10)
FONT_BOLD = Font(name="Arial", size=10, bold=True)

# (Column, Source Tool, Description, Format/Units, Coverage_rows, Coverage_pct, Example)
COLUMNS = [
    ("mutation_id", "pipeline", "Unique key: uniprot_accession|gene_symbol(RBP_name)|wt_aa+position+mut_aa", "string", TOTAL_ROWS, 100.0, "A0A1B0GTU1|ZC3H11B|A115C"),
    ("AGG", "TANGO", "Whole-sequence aggregation (beta-aggregation) propensity score", "float, sum over sequence", TOTAL_ROWS, 100.0, "799.428"),
    ("AMYLO", "TANGO", "Whole-sequence amyloid-forming propensity score", "float, sum over sequence", TOTAL_ROWS, 100.0, "1954.73"),
    ("TURN", "TANGO", "Whole-sequence turn-region propensity score", "float, sum over sequence", TOTAL_ROWS, 100.0, "439.481"),
    ("HELIX", "TANGO", "Whole-sequence alpha-helix propensity score", "float, sum over sequence", TOTAL_ROWS, 100.0, "799.338"),
    ("HELAGG", "TANGO", "Whole-sequence helical-aggregation propensity score", "float, sum over sequence", TOTAL_ROWS, 100.0, "1.21759e-07"),
    ("BETA", "TANGO", "Whole-sequence beta-sheet propensity score", "float, sum over sequence", TOTAL_ROWS, 100.0, "2264.81"),
    ("NUCLEAR_SPECKLE_Score", "protGPS", "Predicted probability of partitioning into nuclear speckles", "float 0-1", TOTAL_ROWS, 100.0, "0.0007"),
    ("P-BODY_Score", "protGPS", "Predicted probability of partitioning into P-bodies", "float 0-1", TOTAL_ROWS, 100.0, "0.0000"),
    ("PML-BDOY_Score", "protGPS", "Predicted probability of partitioning into PML bodies (sic, protGPS's own column name)", "float 0-1", TOTAL_ROWS, 100.0, "0.0002"),
    ("POST_SYNAPTIC_DENSITY_Score", "protGPS", "Predicted probability of partitioning into the post-synaptic density", "float 0-1", TOTAL_ROWS, 100.0, "0.0000"),
    ("STRESS_GRANULE_Score", "protGPS", "Predicted probability of partitioning into stress granules", "float 0-1", TOTAL_ROWS, 100.0, "0.0001"),
    ("CHROMOSOME_Score", "protGPS", "Predicted probability of partitioning onto chromosomes", "float 0-1", TOTAL_ROWS, 100.0, "0.9960"),
    ("NUCLEOLUS_Score", "protGPS", "Predicted probability of partitioning into the nucleolus", "float 0-1", TOTAL_ROWS, 100.0, "0.0002"),
    ("NUCLEAR_PORE_COMPLEX_Score", "protGPS", "Predicted probability of partitioning into the nuclear pore complex", "float 0-1", TOTAL_ROWS, 100.0, "0.0018"),
    ("CAJAL_BODY_Score", "protGPS", "Predicted probability of partitioning into Cajal bodies", "float 0-1", TOTAL_ROWS, 100.0, "0.0000"),
    ("RNA_GRANULE_Score", "protGPS", "Predicted probability of partitioning into RNA granules", "float 0-1", TOTAL_ROWS, 100.0, "0.0000"),
    ("CELL_JUNCTION_Score", "protGPS", "Predicted probability of partitioning into cell junctions", "float 0-1", TOTAL_ROWS, 100.0, "0.0000"),
    ("TRANSCRIPTIONAL_Score", "protGPS", "Predicted probability of partitioning into transcriptional condensates", "float 0-1", TOTAL_ROWS, 100.0, "0.0000"),
    ("Localizations", "DeepLoc2.1", "Predicted subcellular localization category call(s), pipe-joined if multiple", "string, pipe-joined", TOTAL_ROWS, 100.0, "Cytoplasm|Nucleus"),
    ("Signals", "DeepLoc2.1", "Predicted sorting/targeting signal call(s), pipe-joined if multiple; blank when DeepLoc calls no signal", "string, pipe-joined, blank if none", 28973331, 91.90, "Nuclear localization signal|Nuclear export signal"),
    ("Membrane types", "DeepLoc2.1", "Predicted membrane association type call", "string", TOTAL_ROWS, 100.0, "Soluble"),
    ("Cytoplasm", "DeepLoc2.1", "Predicted probability of cytoplasmic localization", "float 0-1", TOTAL_ROWS, 100.0, "0.5228000283241272"),
    ("Nucleus", "DeepLoc2.1", "Predicted probability of nuclear localization", "float 0-1", TOTAL_ROWS, 100.0, "0.7390999794006348"),
    ("Extracellular", "DeepLoc2.1", "Predicted probability of extracellular localization", "float 0-1", TOTAL_ROWS, 100.0, "0.05480000004172325"),
    ("Cell membrane", "DeepLoc2.1", "Predicted probability of cell membrane localization", "float 0-1", TOTAL_ROWS, 100.0, "0.11640000343322754"),
    ("Mitochondrion", "DeepLoc2.1", "Predicted probability of mitochondrial localization", "float 0-1", TOTAL_ROWS, 100.0, "0.07010000199079514"),
    ("Plastid", "DeepLoc2.1", "Predicted probability of plastid localization", "float 0-1", TOTAL_ROWS, 100.0, "0.0019000000320374966"),
    ("Endoplasmic reticulum", "DeepLoc2.1", "Predicted probability of ER localization", "float 0-1", TOTAL_ROWS, 100.0, "0.03139999881386757"),
    ("Lysosome/Vacuole", "DeepLoc2.1", "Predicted probability of lysosome/vacuole localization", "float 0-1", TOTAL_ROWS, 100.0, "0.06109999865293503"),
    ("Golgi apparatus", "DeepLoc2.1", "Predicted probability of Golgi localization", "float 0-1", TOTAL_ROWS, 100.0, "0.08579999953508377"),
    ("Peroxisome", "DeepLoc2.1", "Predicted probability of peroxisomal localization", "float 0-1", TOTAL_ROWS, 100.0, "0.00839999970048666"),
    ("Peripheral", "DeepLoc2.1", "Predicted probability the protein is a peripheral membrane protein", "float 0-1", TOTAL_ROWS, 100.0, "0.4214000105857849"),
    ("Transmembrane", "DeepLoc2.1", "Predicted probability the protein is transmembrane", "float 0-1", TOTAL_ROWS, 100.0, "0.0803999975323677"),
    ("Lipid anchor", "DeepLoc2.1", "Predicted probability the protein is lipid-anchored", "float 0-1", TOTAL_ROWS, 100.0, "0.23999999463558197"),
    ("Soluble", "DeepLoc2.1", "Predicted probability the protein is soluble", "float 0-1", TOTAL_ROWS, 100.0, "0.8159999847412109"),
    ("slim_n_motifs_overlapping", "ELM (this repo)", "Count of distinct ELM SLiM match instances whose span covers this mutation's position", "int, blank if 0", 9296415, 29.48, "2"),
    ("slim_elm_identifiers", "ELM (this repo)", "Semicolon-joined ELM class identifiers of every overlapping motif", "string, ;-joined, blank if none", 9296415, 29.48, "MOD_SUMO_rev_2;MOD_PIKK_1"),
    ("slim_any_disrupts_motif", "ELM (this repo)", "Y if this mutation breaks at least one overlapping motif's regex, N if all overlapping motifs still match, blank if no overlap", "Y/N, blank if no overlap", 9296415, 29.48, "N"),
    ("slim_disrupted_elm_identifiers", "ELM (this repo)", "Semicolon-joined ELM class identifiers of specifically the motifs this mutation disrupts (subset of slim_elm_identifiers)", "string, ;-joined, blank if none disrupted", 5336455, 16.93, "MOD_SUMO_rev_2"),
]

GROUPS = [
    ("TANGO -- aggregation propensity", "TANGO", 6, TOTAL_ROWS, 100.0,
     "Whole-sequence aggregation/amyloid/structural propensity scores for every point mutant. "
     "Raw mutant values only -- no WT baseline run yet, so no delta column exists."),
    ("protGPS -- condensate partitioning", "protGPS", 12, TOTAL_ROWS, 100.0,
     "Predicted probability of partitioning into each of 12 biomolecular condensate types, per mutant sequence."),
    ("DeepLoc2.1 -- subcellular localization", "DeepLoc2.1", 17, TOTAL_ROWS, 100.0,
     "Predicted subcellular localization/signal/membrane-type calls plus 13 per-compartment probability scores, per mutant sequence."),
    ("ELM SLiM motif disruption", "ELM (this repo)", 4, 9296415, 29.48,
     "For mutations that fall inside an ELM motif match already annotated on the WT sequence (asharma's SLiM_* columns), "
     "whether the mutant residue breaks that motif's own regex. Sparse by construction: only ~29.5% of all mutations "
     "overlap any annotated motif at all; of those, 57.4% are flagged as disrupting it."),
]


def style_header(ws, row, ncols):
    for col in range(1, ncols + 1):
        c = ws.cell(row=row, column=col)
        c.font = FONT_BOLD
        c.alignment = Alignment(wrap_text=False)


def autosize(ws, widths):
    for i, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = w


wb = openpyxl.Workbook()

# --- About ---
ws = wb.active
ws.title = "About"
about_rows = [
    ("satmut_all_predictions.tsv -- documentation", ""),
    ("", ""),
    ("Path (Hoffman2)", "/u/project/kappel/ddcohn/protein_variant_effects/rbp_only/satmut_all_predictions.tsv"),
    ("Total rows", f"{TOTAL_ROWS:,}"),
    ("Total columns", str(len(COLUMNS))),
    ("File size", "15.8 GB"),
    ("Row semantics", "One row per point mutation (full saturation scan: every residue x all 19 other amino acids, across 2,336 RBPs). NOT one row per protein."),
    ("Mutation universe source", "asharma's rbp_satmut_all.csv.gz (33,184,720 rows incl. WT refs; 31,525,484 actual mutant rows used here)"),
    ("", ""),
    ("Build lineage", ""),
    ("1", "build_chunks_from_satmut_source.py -- 226 FASTA chunks of 140,000 mutant sequences each from asharma's source"),
    ("2", "TANGO, protGPS, DeepLoc2.1 each run as a 226-task array job over those chunks (see saturation_mutagenesis/README.md for full scheduling/GPU/wall-time history)"),
    ("3", "Each tool's 226 chunks verified complete via full-array content audit (real output row count == real input FASTA count, every chunk) -- 0 incomplete for all three"),
    ("4", "merge_satmut_tango.py / merge_satmut_protgps.py / merge_satmut_deeploc.py -- concatenate each tool's 226 chunks into one file"),
    ("5", "compute_slim_disruption.py -- for every ELM SLiM match already annotated on each protein's WT sequence (asharma's merged table), enumerate every point mutation inside that match's span and test it against the motif's own ELM regex"),
    ("6", "fix_slim_gene_names.py -- corrected a gene-alias mismatch (68 proteins, e.g. O15446=POLR1G/CD3EAP) between the merged table's gene_symbol and the satmut pipeline's own RBP_name, which had caused 233,985 SLiM-derived mutation_ids to silently fail to join"),
    ("7", "combine_satmut_predictions.py -- external sort + 4-way streaming merge-join (memory-safe; a full in-RAM dict join would need tens of GB) of TANGO + protGPS + DeepLoc + aggregated SLiM disruption calls by mutation_id"),
    ("", ""),
    ("Key caveats", ""),
    ("1", "No WT baseline / delta columns for any of the four sources -- these are raw mutant values only. A WT run would be needed first to compute deltas (see protGPS/DeepLoc WT additions to the separate protein-level table for the analogous pattern)."),
    ("2", "SLiM columns are sparse by design: blank for the ~70.5% of mutations that don't fall inside any ELM motif match already annotated on the WT sequence. A mutation inside >1 overlapping motif has its motifs semicolon-joined, not duplicated as extra rows."),
    ("3", "This table does not re-run ELM's motif scan itself -- it reuses asharma's existing WT-level SLiM_range/SLiM_class/SLiM_discrete_seq annotations (from the protein-level merged table) as the set of candidate motif spans to test mutations against."),
]
for r in about_rows:
    ws.append(r)
ws["A1"].font = Font(name="Arial", size=12, bold=True)
for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
    for cell in row:
        if cell.font.bold:
            continue
        cell.font = FONT
autosize(ws, [20, 110])

# --- Columns ---
ws = wb.create_sheet("Columns")
header = ["#", "Column", "Source Tool", "Description", "Format/Units", "Non-empty rows", "Coverage %", "Example"]
ws.append(header)
for i, (col, tool, desc, fmt, nonempty, pct, example) in enumerate(COLUMNS, start=1):
    ws.append([i, col, tool, desc, fmt, nonempty, f"{pct:.2f}%", example])
style_header(ws, 1, len(header))
for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
    for cell in row:
        cell.font = FONT
    row[5].number_format = "#,##0"
autosize(ws, [4, 30, 16, 65, 22, 15, 11, 45])
ws.freeze_panes = "A2"

# --- Groups ---
ws = wb.create_sheet("Groups")
header = ["Group", "Source Tool", "# Columns", "Rows with data", "Coverage %", "Notes"]
ws.append(header)
for name, tool, ncols, nonempty, pct, notes in GROUPS:
    ws.append([name, tool, ncols, nonempty, f"{pct:.2f}%", notes])
style_header(ws, 1, len(header))
for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
    for cell in row:
        cell.font = FONT
    row[3].number_format = "#,##0"
    row[5].alignment = Alignment(wrap_text=True)
autosize(ws, [32, 16, 11, 16, 11, 90])
ws.freeze_panes = "A2"

wb.save(OUT_PATH)
print(f"Wrote {OUT_PATH}")
