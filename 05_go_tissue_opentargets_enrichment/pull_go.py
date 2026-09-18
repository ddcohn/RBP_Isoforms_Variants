import csv
import sys
import json
import time
import urllib.request
import urllib.error
import argparse

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
CHECKPOINT = "/u/home/d/ddcohn/_claude_go_checkpoint.json"
RESULT = "/u/home/d/ddcohn/_claude_go_result.json"

BATCH_SIZE = 100


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


def fetch_batch(accs):
    """Returns dict: accession -> {"C": [(id,term),...], "P": [...], "F": [...]}"""
    url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions="
           + ",".join(accs) + "&fields=accession,go_p,go_c,go_f&size=" + str(len(accs)))
    for attempt in range(6):
        try:
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=90) as resp:
                d = json.loads(resp.read().decode())
            break
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(3 * (attempt + 1))
                continue
            print(f"  HTTP error {e.code} on batch starting {accs[0]}", file=sys.stderr)
            return {a: {"C": [], "P": [], "F": []} for a in accs}
        except Exception as e:
            print(f"  request exception (attempt {attempt+1}): {e}", file=sys.stderr)
            time.sleep(3)
    else:
        return {a: {"C": [], "P": [], "F": []} for a in accs}

    out = {a: {"C": [], "P": [], "F": []} for a in accs}
    for entry in d.get("results", []):
        acc = entry.get("primaryAccession")
        if acc not in out:
            continue
        for xref in entry.get("uniProtKBCrossReferences", []):
            if xref.get("database") != "GO":
                continue
            go_id = xref.get("id", "")
            term = None
            for prop in xref.get("properties", []):
                if prop.get("key") == "GoTerm":
                    term = prop.get("value", "")
                    break
            if not term or ":" not in term:
                continue
            aspect, name = term.split(":", 1)
            if aspect in out[acc]:
                out[acc][aspect].append((go_id, name))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    accs = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            accs.append(row["uniprot_accession"])
            if args.limit and len(accs) >= args.limit:
                break

    print(f"Total proteins to query: {len(accs)}")

    results = load_json(CHECKPOINT) or {}
    if results:
        print(f"  resuming with {len(results)} already cached", flush=True)

    to_fetch = [a for a in accs if a not in results]
    start = time.time()
    for i in range(0, len(to_fetch), BATCH_SIZE):
        batch = to_fetch[i:i + BATCH_SIZE]
        batch_result = fetch_batch(batch)
        results.update(batch_result)
        elapsed = time.time() - start
        done = i + len(batch)
        print(f"  ... {done}/{len(to_fetch)} done (elapsed {elapsed/60:.1f}min, "
              f"avg {elapsed/done:.3f}s/protein)", flush=True)
        if (i // BATCH_SIZE) % 10 == 0:
            save_json(CHECKPOINT, results)
        time.sleep(0.2)
    save_json(CHECKPOINT, results)

    n_cc = sum(1 for v in results.values() if v["C"])
    n_bp = sum(1 for v in results.values() if v["P"])
    n_mf = sum(1 for v in results.values() if v["F"])
    print(f"\nProteins with >=1 Cellular Component term: {n_cc}/{len(results)}")
    print(f"Proteins with >=1 Biological Process term:  {n_bp}/{len(results)}")
    print(f"Proteins with >=1 Molecular Function term:  {n_mf}/{len(results)}")

    if not args.dry_run:
        save_json(RESULT, results)
        print(f"\nSaved per-protein GO result to {RESULT}")


if __name__ == "__main__":
    main()
