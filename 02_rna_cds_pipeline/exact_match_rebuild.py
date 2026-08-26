import csv
import sys
import json
import time
import urllib.request
import urllib.error
import argparse

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/id?expand=1"
ENSEMBL_SEQ_CDS = "https://rest.ensembl.org/sequence/id?type=cds;content-type=application/json"
ENSEMBL_SEQ_CDNA = "https://rest.ensembl.org/sequence/id?type=cdna;content-type=application/json"

GENE_CHECKPOINT = "/u/home/d/ddcohn/_claude_exact_match_gene_checkpoint.json"
CDS_CHECKPOINT = "/u/home/d/ddcohn/_claude_exact_match_cds_checkpoint.json"


def save_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f)
    import os
    os.replace(tmp, path)


def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return None

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


def translate(cds):
    prot = []
    for i in range(0, len(cds) - len(cds) % 3, 3):
        aa = CODON_TABLE.get(cds[i:i+3], 'X')
        if aa == '*':
            break
        prot.append(aa)
    return ''.join(prot)


def http_post_json(url, payload, timeout=90):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def lookup_genes(bare_ensgs, batch_size=100):
    """Returns dict: ENSG -> list of transcript IDs that have a Translation.
    Resumes from checkpoint and saves progress periodically, since this phase
    is slow (large nested JSON per gene) and previously got killed by a
    wall-clock limit partway through."""
    result = load_json(GENE_CHECKPOINT) or {}
    if result:
        print(f"  resuming gene lookup with {len(result)} genes already cached", flush=True)
    ids = [g for g in bare_ensgs if g not in result]
    start = time.time()
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        t0 = time.time()
        body = None
        for attempt in range(6):
            try:
                body = http_post_json(ENSEMBL_LOOKUP, {"ids": batch})
                break
            except urllib.error.HTTPError as e:
                if e.code in (429, 500, 502, 503, 504):
                    time.sleep(3 * (attempt + 1))
                    continue
                print(f"  lookup batch error {e.code} at offset {i}", file=sys.stderr)
                body = None
                break
            except Exception as e:
                print(f"  lookup batch exception at offset {i}: {e}", file=sys.stderr)
                time.sleep(3)
                continue
        if not body:
            for gid in batch:
                result.setdefault(gid, [])
            continue
        try:
            d = json.loads(body)
        except json.JSONDecodeError:
            for gid in batch:
                result.setdefault(gid, [])
            continue
        for gid in batch:
            info = d.get(gid)
            if not info:
                result[gid] = []
                continue
            transcripts = info.get("Transcript", [])
            result[gid] = [t["id"] for t in transcripts if t.get("Translation")]
        elapsed = time.time() - start
        done = i + len(batch)
        print(f"  ... gene lookup {done}/{len(ids)} done (batch took {time.time()-t0:.1f}s, "
              f"avg {elapsed/done:.2f}s/gene, elapsed {elapsed/60:.1f}min)", flush=True)
        if (i // batch_size) % 10 == 0:
            save_json(GENE_CHECKPOINT, result)
    save_json(GENE_CHECKPOINT, result)
    return result


def fetch_sequences(bare_ids, batch_size=50, url=ENSEMBL_SEQ_CDS, label="seq", checkpoint=None):
    seqs = load_json(checkpoint) if checkpoint else None
    seqs = seqs or {}
    if seqs:
        print(f"  resuming {label} fetch with {len(seqs)} already cached", flush=True)
    ids = [x for x in bare_ids if x not in seqs]
    start = time.time()
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        payload = json.dumps({"ids": batch}).encode()
        body = "[]"
        for attempt in range(6):
            try:
                req = urllib.request.Request(url, data=payload, headers={"Content-Type": "application/json"})
                with urllib.request.urlopen(req, timeout=90) as resp:
                    body = resp.read().decode()
                break
            except urllib.error.HTTPError as e:
                if e.code in (429, 500, 502, 503, 504):
                    time.sleep(3 * (attempt + 1))
                    continue
                body = e.read().decode()
                print(f"  {label} batch error {e.code}: {body[:200]}", file=sys.stderr)
                body = "[]"
                break
            except Exception as e:
                print(f"  {label} batch exception (attempt {attempt+1}): {e}", file=sys.stderr)
                time.sleep(3)
                continue
        try:
            results = json.loads(body)
        except json.JSONDecodeError:
            results = []
        if isinstance(results, dict):
            results = [results]
        for r in results:
            if isinstance(r, dict) and "id" in r and "seq" in r:
                seqs[r["id"]] = r["seq"]
        elapsed = time.time() - start
        done = i + len(batch)
        print(f"  ... {label} fetch {done}/{len(ids)} done (elapsed {elapsed/60:.1f}min)", flush=True)
        if checkpoint and (i // batch_size) % 20 == 0:
            save_json(checkpoint, seqs)
    if checkpoint:
        save_json(checkpoint, seqs)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", default=TARGET)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    # candidate rows = every row with an ensembl_gene AND a protein sequence
    candidates = {}  # uid -> bare ENSG
    for row in rows:
        g = (row.get("ensembl_gene") or "").strip()
        prot = (row.get("sequence") or "").strip()
        if not g or not prot:
            continue
        first_ensg = g.split(";")[0].strip().split(".")[0]
        if first_ensg:
            candidates[row["uniprot_accession"]] = first_ensg
            if args.limit and len(candidates) >= args.limit:
                break

    print(f"Candidates (have ensembl_gene + protein sequence): {len(candidates)}")

    unique_ensgs = set(candidates.values())
    print(f"Distinct genes to look up: {len(unique_ensgs)}")

    gene_transcripts = lookup_genes(unique_ensgs)
    total_transcripts = sum(len(v) for v in gene_transcripts.values())
    print(f"Total translated transcripts across all genes: {total_transcripts}")

    all_transcripts = sorted(set(t for v in gene_transcripts.values() for t in v))
    print(f"Fetching CDS for {len(all_transcripts)} distinct transcripts...")
    cds_seqs = fetch_sequences(all_transcripts, url=ENSEMBL_SEQ_CDS, label="cds", checkpoint=CDS_CHECKPOINT)
    print(f"Got CDS for {len(cds_seqs)} / {len(all_transcripts)}")

    # translate all fetched CDS once
    translations = {tid: translate(seq) for tid, seq in cds_seqs.items()}

    # build uid -> protein sequence lookup once (was previously an O(n^2) scan)
    uid_to_protein = {r["uniprot_accession"]: r["sequence"] for r in rows}

    # for each row, find the transcript among its gene's set whose translation exactly matches
    exact_matches = {}   # uid -> transcript_id
    no_exact_match = []  # uid, ensg, protlen, n_transcripts_tried
    for uid, ensg in candidates.items():
        prot = uid_to_protein[uid]
        transcripts = gene_transcripts.get(ensg, [])
        found = None
        for tid in transcripts:
            if translations.get(tid) == prot:
                found = tid
                break
        if found:
            exact_matches[uid] = found
        else:
            no_exact_match.append((uid, ensg, len(prot), len(transcripts)))

    print(f"\nExact-match transcript found: {len(exact_matches)} / {len(candidates)}")
    print(f"No exact match among gene's transcripts: {len(no_exact_match)}")
    print("Sample of no-exact-match rows (uid, ensg, protlen, n_transcripts_tried):")
    for ex in no_exact_match[:15]:
        print(" ", ex)

    # fetch full mRNA (cdna) for the winning transcripts
    winning_transcripts = sorted(set(exact_matches.values()))
    print(f"\nFetching full mRNA for {len(winning_transcripts)} winning transcripts...")
    CDNA_CHECKPOINT = "/u/home/d/ddcohn/_claude_exact_match_cdna_checkpoint.json"
    cdna_seqs = fetch_sequences(winning_transcripts, url=ENSEMBL_SEQ_CDNA, label="cdna", checkpoint=CDNA_CHECKPOINT)
    print(f"Got mRNA for {len(cdna_seqs)} / {len(winning_transcripts)}")

    n_updated = 0
    for row in rows:
        uid = row["uniprot_accession"]
        if uid not in exact_matches:
            continue
        tid = exact_matches[uid]
        row["ensembl_transcript_used"] = tid
        row["transcript_selection_method"] = "exact_translation_match"
        row["cds_sequence"] = cds_seqs.get(tid, "")
        row["rna_sequence"] = cdna_seqs.get(tid, "")
        n_updated += 1

    if not args.dry_run:
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    print(f"\nUpdated {n_updated} rows with exact-match transcripts (dry_run={args.dry_run})")


if __name__ == "__main__":
    main()
