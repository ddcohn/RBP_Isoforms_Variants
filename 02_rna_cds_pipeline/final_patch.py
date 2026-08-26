import csv
import sys
import json
import time
import urllib.request
import urllib.error

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
ENSEMBL_SEQ = "https://rest.ensembl.org/sequence/id/{}?type=cdna;content-type=application/json"

# alternate transcripts (same isoform tag as original pick) to try, in order, per accession
ALTERNATES = {
    "O95229": ["ENST00000373944", "ENST00000395405", "ENST00000899408"],
    "O95996": ["ENST00000535453", "ENST00000590469"],
    "Q2NL67": ["ENST00000569795", "ENST00000907970", "ENST00000907971"],
    "Q86W28": ["ENST00000590542"],  # isoform -2, best available fallback
    "Q86YA3": ["ENST00000505019"],
    "Q96P50": ["ENST00000354700"],  # isoform -3, best available fallback
    "Q9NT99": ["ENST00000599957", "ENST00000652263"],
    "Q9UPU9": ["ENST00000554335"],
}
# genuinely no replacement exists in current Ensembl (confirmed via /archive/id)
UNRECOVERABLE = {"O60393", "P59047"}


def try_fetch(enst_bare):
    url = ENSEMBL_SEQ.format(enst_bare)
    for attempt in range(5):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=30) as resp:
                body = resp.read().decode()
            d = json.loads(body)
            if "seq" in d:
                return d["seq"]
            return None
        except urllib.error.HTTPError as e:
            if e.code in (429, 503):
                time.sleep(5)
                continue
            return None
        except Exception:
            time.sleep(3)
            continue
    return None


def main():
    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    filled = 0
    for row in rows:
        acc = row["uniprot_accession"].strip()
        if acc in ALTERNATES and not row.get("rna_sequence"):
            for alt in ALTERNATES[acc]:
                seq = try_fetch(alt)
                if seq:
                    row["rna_sequence"] = seq
                    row["ensembl_transcript_used"] = alt
                    row["transcript_selection_method"] = "alternate_after_retired_id"
                    print(f"{acc}: recovered via {alt} ({len(seq)} nt)")
                    filled += 1
                    break
            else:
                print(f"{acc}: no working alternate found among {ALTERNATES[acc]}")
        elif acc in UNRECOVERABLE:
            print(f"{acc}: confirmed unrecoverable (Ensembl transcript retired, no replacement)")

    with open(TARGET, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    total_with_rna = sum(1 for r in rows if r.get("rna_sequence"))
    total_with_transcript = sum(1 for r in rows if r.get("ensembl_transcript_used"))
    print(f"\nFilled {filled} more rows via alternates.")
    print(f"Final: {total_with_rna} / {len(rows)} rows have rna_sequence")
    print(f"({total_with_transcript} rows had a transcript assigned at all; "
          f"{total_with_transcript - total_with_rna} of those still lack a sequence)")


if __name__ == "__main__":
    main()
