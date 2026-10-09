"""Step 1: dump Ensembl BioMart tables needed for the RBP conservation table.

Writes (skips any file that already exists, so a killed run resumes):
  data/raw/tx_map_<group>.tsv      protein-coding transcript attribute groups
  data/raw/ortho_<species>.tsv     human -> <species> Compara orthology metrics
  data/raw/biomart_release.txt     Ensembl release the dumps came from
Run from the 13_evolutionary_conservation directory.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import biomart as bm  # noqa: E402

RAW = Path("data/raw")
RAW.mkdir(parents=True, exist_ok=True)
HS = "hsapiens_gene_ensembl"
CODING = {"transcript_biotype": "protein_coding"}

# One attribute group per query (BioMart refuses attributes from >1 page).
ID = ["ensembl_gene_id", "ensembl_transcript_id"]
TX_GROUPS = {
    "names_uniprot": ID + ["external_gene_name", "uniprotswissprot", "chromosome_name"],
    "canonical": ID + ["transcript_is_canonical"],
    "mane": ID + ["transcript_mane_select"],
    "cds_length": ID + ["cds_length"],
}

# (BioMart species prefix, label) -- a ladder from close primates to yeast.
SPECIES = [
    "ptroglodytes", "mmulatta", "mmusculus", "rnorvegicus", "ggallus",
    "xtropicalis", "drerio", "dmelanogaster", "celegans", "scerevisiae",
]
# NB: Ensembl 116 BioMart exposes no <species>_homolog_dn / _ds attributes
# (checked against the dataset's attribute list), so dN/dS is not available here.
ORTHO_FIELDS = [
    "homolog_ensembl_gene", "homolog_orthology_type", "homolog_perc_id",
    "homolog_perc_id_r1", "homolog_goc_score",
    "homolog_wga_coverage", "homolog_orthology_confidence",
]


def main():
    (RAW / "biomart_release.txt").write_text(bm.release() + "\n")
    for name, attrs in TX_GROUPS.items():
        out = RAW / f"tx_map_{name}.tsv"
        if out.exists():
            continue
        df = bm.query(HS, attrs, CODING)
        df.to_csv(out, sep="\t", index=False)
        print(f"{out}: {df.shape}", flush=True)

    avail = bm.attribute_names(HS)
    for sp in SPECIES:
        out = RAW / f"ortho_{sp}.tsv"
        if out.exists():
            continue
        # keep only the fields this species actually exposes (distant species
        # lack GOC / WGA-coverage attributes)
        fields = [f"{sp}_{f}" for f in ORTHO_FIELDS if f"{sp}_{f}" in avail]
        missing = [f for f in ORTHO_FIELDS if f"{sp}_{f}" not in avail]
        df = bm.query(HS, ["ensembl_gene_id"] + fields, {"biotype": "protein_coding"})
        df.to_csv(out, sep="\t", index=False)
        print(f"{out}: {df.shape}; unavailable fields: {missing}", flush=True)


if __name__ == "__main__":
    main()
