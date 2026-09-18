import csv
import sys
import json
import time
import urllib.request
import urllib.error
import argparse

csv.field_size_limit(sys.maxsize)

TARGET = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

BATCH_SIZE = 5        # genes per GraphQL request
PAGE_SIZE = 1500      # rows per gene per request (covers single-cell noise we'll filter out; well-studied genes top out ~1450-1460 total rows). batch*page above ~10000 trips the API's "too expensive" guard.

CHECKPOINT = "/u/home/d/ddcohn/_claude_tissue_checkpoint.json"
RESULT = "/u/home/d/ddcohn/_claude_tissue_result.json"


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


def graphql_post(query, timeout=120):
    payload = json.dumps({"query": query}).encode()
    req = urllib.request.Request(OT_URL, data=payload, headers={"Content-Type": "application/json"})
    for attempt in range(6):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            body = e.read().decode()
            if e.code in (429, 500, 502, 503, 504):
                time.sleep(3 * (attempt + 1))
                continue
            print(f"  HTTP error {e.code}: {body[:300]}", file=sys.stderr)
            return None
        except Exception as e:
            print(f"  request exception (attempt {attempt+1}): {e}", file=sys.stderr)
            time.sleep(3)
    return None


def bulk_rows_only(rows):
    """Keep only tissue-level rows (drop per-celltype single-cell breakdowns)."""
    out = []
    for r in rows:
        if r.get("celltypeBiosampleFromSource"):
            continue
        out.append({
            "tissue": (r.get("tissueBiosample") or {}).get("biosampleName") or r.get("tissueBiosampleFromSource") or "",
            "datasourceId": r.get("datasourceId"),
            "datatypeId": r.get("datatypeId"),
            "unit": r.get("unit"),
            "median": r.get("median"),
        })
    return out


def fetch_batch(gene_batch):
    """Returns dict: ensg -> list of bulk tissue expression rows."""
    parts = []
    for i, g in enumerate(gene_batch):
        parts.append(
            f't{i}: target(ensemblId: "{g}") {{ baselineExpression(page:{{index:0,size:{PAGE_SIZE}}}) '
            f'{{ count rows {{ datasourceId datatypeId unit median '
            f'celltypeBiosampleFromSource tissueBiosampleFromSource '
            f'tissueBiosample {{ biosampleName }} }} }} }}'
        )
    query = "query{" + " ".join(parts) + "}"
    resp = graphql_post(query)
    result = {}
    if not resp or "data" not in resp:
        return {g: [] for g in gene_batch}
    data = resp["data"]
    for i, g in enumerate(gene_batch):
        v = data.get(f"t{i}")
        if not v:
            result[g] = []
            continue
        be = v["baselineExpression"]
        if be["count"] > len(be["rows"]):
            print(f"  WARNING: {g} truncated ({be['count']} total, only fetched {len(be['rows'])})", file=sys.stderr)
        result[g] = bulk_rows_only(be["rows"])
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    rows = []
    with open(TARGET, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    candidates = {}
    for row in rows:
        g = (row.get("ensembl_gene") or "").strip()
        if not g:
            continue
        first_ensg = g.split(";")[0].strip().split(".")[0]
        if first_ensg:
            candidates[row["uniprot_accession"]] = first_ensg
            if args.limit and len(candidates) >= args.limit:
                break

    print(f"Candidates (have ensembl_gene): {len(candidates)}")
    unique_ensgs = sorted(set(candidates.values()))
    print(f"Distinct genes to query: {len(unique_ensgs)}")

    gene_results = load_json(CHECKPOINT) or {}
    if gene_results:
        print(f"  resuming with {len(gene_results)} genes already cached", flush=True)

    to_fetch = [g for g in unique_ensgs if g not in gene_results]
    start = time.time()
    for i in range(0, len(to_fetch), BATCH_SIZE):
        batch = to_fetch[i:i + BATCH_SIZE]
        batch_result = fetch_batch(batch)
        gene_results.update(batch_result)
        elapsed = time.time() - start
        done = i + len(batch)
        print(f"  ... {done}/{len(to_fetch)} genes done (elapsed {elapsed/60:.1f}min, "
              f"avg {elapsed/done:.3f}s/gene)", flush=True)
        if (i // BATCH_SIZE) % 20 == 0:
            save_json(CHECKPOINT, gene_results)
        time.sleep(0.2)
    save_json(CHECKPOINT, gene_results)

    total_rows = sum(len(v) for v in gene_results.values())
    n_with_data = sum(1 for v in gene_results.values() if v)
    print(f"Genes with >=1 bulk tissue row: {n_with_data}/{len(gene_results)}")
    print(f"Total bulk tissue-expression rows across all genes: {total_rows}")

    result = {}
    for uid, ensg in candidates.items():
        result[uid] = gene_results.get(ensg, [])

    if not args.dry_run:
        save_json(RESULT, result)
        print(f"Saved per-protein tissue result to {RESULT}")


if __name__ == "__main__":
    main()
