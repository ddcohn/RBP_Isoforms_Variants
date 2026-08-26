import csv
import sys
import json
import time
import urllib.request
import urllib.parse
import urllib.error
import argparse

csv.field_size_limit(sys.maxsize)

SRC = "/u/project/kappel/ddcohn/RBP/table_260823_with_rna.csv"

STRING_RESOLVE_URL = "https://string-db.org/api/tsv/get_string_ids"
STRING_URL = "https://string-db.org/api/tsv/interaction_partners"
UNIPROT_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_STATUS = "https://rest.uniprot.org/idmapping/status/{}"
UNIPROT_STREAM = "https://rest.uniprot.org/idmapping/stream/{}?format=tsv"

REQUIRED_SCORE = 700  # STRING "high confidence"
STRING_LIMIT = 100000  # override STRING's default cap of 10 partners/protein

CHECKPOINT = "/u/home/d/ddcohn/_claude_ppi_checkpoint.json"


def save_checkpoint(resolved, partner_map):
    tmp = CHECKPOINT + ".tmp"
    with open(tmp, "w") as f:
        json.dump({
            "resolved": resolved,
            "partner_map": {k: sorted(v) for k, v in partner_map.items()},
        }, f)
    import os
    os.replace(tmp, CHECKPOINT)


def load_checkpoint():
    try:
        with open(CHECKPOINT) as f:
            d = json.load(f)
        return d["resolved"], {k: set(v) for k, v in d["partner_map"].items()}
    except (FileNotFoundError, json.JSONDecodeError, KeyError):
        return None, None


def http_get(url, timeout=60):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def http_post(url, data_dict, timeout=60):
    data = urllib.parse.urlencode(data_dict).encode()
    req = urllib.request.Request(url, data=data)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def resolve_string_ids(bare_ensps, batch_size=100):
    """STRING silently renames/normalizes some input IDs to a current canonical
    STRING ID before returning interaction data, so we must resolve first and
    query using the resolved ID to reliably match results back to the query.
    Returns dict: input_bare_ensp -> resolved_bare_ensp (only for IDs STRING recognizes)."""
    resolved = {}
    ids = list(bare_ensps)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        query = "%0d".join(batch)
        url = f"{STRING_RESOLVE_URL}?identifiers={query}&species=9606"
        body = None
        for attempt in range(4):
            try:
                body = http_get(url, timeout=90)
                break
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(8)
                    continue
                print(f"  resolve batch error {e.code} at offset {i}", file=sys.stderr)
                body = ""
                break
            except Exception as e:
                print(f"  resolve batch exception at offset {i}: {e}", file=sys.stderr)
                time.sleep(3)
                body = ""
        if not body:
            continue
        lines = body.strip().split("\n")
        if not lines or lines[0].startswith("Error"):
            continue
        for line in lines[1:]:
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            try:
                idx = int(cols[0])
            except ValueError:
                continue
            string_id = cols[1]
            bare = string_id.split(".", 1)[1] if "." in string_id else string_id
            if idx < len(batch):
                resolved[batch[idx]] = bare
        if (i // batch_size) % 20 == 0:
            print(f"  ... resolve {i+len(batch)}/{len(ids)} done", flush=True)
        time.sleep(0.4)
    return resolved


def fetch_string_partners(bare_ensps, batch_size=50):
    """Returns dict: query_ensp(bare) -> list of partner_ensp(bare)."""
    partners = {}
    ids = list(bare_ensps)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        query = "%0d".join(batch)
        url = (f"{STRING_URL}?identifiers={query}&species=9606"
               f"&required_score={REQUIRED_SCORE}&limit={STRING_LIMIT}")
        body = None
        for attempt in range(4):
            try:
                body = http_get(url, timeout=90)
                break
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(8)
                    continue
                print(f"  STRING batch error {e.code} at offset {i}", file=sys.stderr)
                body = ""
                break
            except Exception as e:
                print(f"  STRING batch exception at offset {i}: {e}", file=sys.stderr)
                time.sleep(3)
                body = ""
        if not body:
            continue
        lines = body.strip().split("\n")
        if not lines or lines[0].startswith("Error"):
            continue
        for line in lines[1:]:
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            sid_a, sid_b = cols[0], cols[1]
            a = sid_a.split(".", 1)[1] if "." in sid_a else sid_a
            b = sid_b.split(".", 1)[1] if "." in sid_b else sid_b
            partners.setdefault(a, set()).add(b)
        if (i // batch_size) % 20 == 0:
            print(f"  ... STRING queries {i+len(batch)}/{len(ids)} done", flush=True)
        time.sleep(0.5)
    return partners


def submit_mapping(ids, from_db, to_db):
    resp = http_post(UNIPROT_RUN, {
        "ids": ",".join(ids),
        "from": from_db,
        "to": to_db,
    })
    return json.loads(resp)["jobId"]


def wait_for_job(job_id, timeout=1200):
    """Polls job status, tolerating transient HTTP/network errors from UniProt's
    backend (e.g. a momentary 'Connection refused' on their side) instead of
    crashing on the first blip."""
    start = time.time()
    while time.time() - start < timeout:
        try:
            resp = http_get(UNIPROT_STATUS.format(job_id))
        except Exception as e:
            print(f"  status check transient error: {e}, retrying...", file=sys.stderr)
            time.sleep(8)
            continue
        try:
            obj = json.loads(resp)
        except json.JSONDecodeError:
            time.sleep(8)
            continue
        if "jobStatus" in obj:
            if obj["jobStatus"] in ("RUNNING", "NEW"):
                time.sleep(4)
                continue
            else:
                raise RuntimeError(f"job failed: {obj}")
        else:
            return
    raise TimeoutError("mapping job did not finish in time")


def run_mapping_with_retry(ids, from_db, to_db, max_resubmits=5):
    """Submits a UniProt ID mapping job; if the job itself errors out (e.g. a
    transient outage on UniProt's backend), resubmits fresh rather than giving up."""
    for attempt in range(max_resubmits):
        job_id = submit_mapping(ids, from_db, to_db)
        print(f"  job id: {job_id} (attempt {attempt+1}/{max_resubmits})")
        try:
            wait_for_job(job_id)
            return job_id
        except RuntimeError as e:
            print(f"  job errored: {e}; resubmitting...", file=sys.stderr)
            time.sleep(10)
    raise RuntimeError(f"UniProt mapping job failed after {max_resubmits} attempts")


def get_mapping_results(job_id):
    tsv = http_get(UNIPROT_STREAM.format(job_id))
    lines = tsv.strip().split("\n")
    mapping = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        frm, to = line.split("\t")
        mapping.setdefault(frm, [])
        if to not in mapping[frm]:
            mapping[frm].append(to)
    return mapping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="only process first N rows (testing)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = []
    with open(SRC, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)
            if args.limit and len(rows) >= args.limit:
                break

    all_uniprot_ids = set(r["uniprot_accession"].strip() for r in rows if r["uniprot_accession"].strip())
    print(f"Rows: {len(rows)}; distinct uniprot IDs in dataset: {len(all_uniprot_ids)}")

    # query ENSP per row = first listed, version-stripped
    query_ensp_by_uid = {}
    ensp_to_uid = {}  # bare ENSP -> uniprot_accession, built from THIS dataset (for in-dataframe resolution + fallback)
    for r in rows:
        uid = r["uniprot_accession"].strip()
        ensp_field = (r.get("ensembl_protein") or "").strip()
        if not ensp_field:
            continue
        first_ensp = ensp_field.split(";")[0].strip()
        bare = first_ensp.split(".")[0]
        if bare:
            query_ensp_by_uid[uid] = bare
            ensp_to_uid.setdefault(bare, uid)

    print(f"Rows with a queryable Ensembl protein ID: {len(query_ensp_by_uid)}")

    resolved, partner_map = load_checkpoint()
    if resolved is not None:
        print(f"Resuming from checkpoint: {len(resolved)} resolved IDs, {len(partner_map)} query proteins with partner data")
    else:
        print("Resolving input Ensembl protein IDs to current STRING IDs...")
        resolved = resolve_string_ids(set(query_ensp_by_uid.values()))
        print(f"  resolved {len(resolved)} / {len(set(query_ensp_by_uid.values()))} input ENSPs (unresolved = not in STRING, e.g. non-canonical entries)")

        print("Querying STRING for interaction partners (required_score=700, no cap)...")
        partner_map = fetch_string_partners(set(resolved.values()))
        print(f"Got STRING partner data for {len(partner_map)} resolved query proteins")

        save_checkpoint(resolved, partner_map)
        print("  (checkpoint saved)")

    all_partner_ensps = set()
    for plist in partner_map.values():
        all_partner_ensps.update(plist)
    print(f"Distinct partner ENSPs across whole dataset: {len(all_partner_ensps)}")

    # Only need to look up partners not already resolvable via our own table's ENSP->uniprot map
    need_lookup = sorted(e for e in all_partner_ensps if e not in ensp_to_uid)
    print(f"Partner ENSPs needing UniProt lookup via API: {len(need_lookup)}")

    ensp_to_partner_uid = dict(ensp_to_uid)  # start with what we know from our own table
    if need_lookup:
        print("Submitting UniProt ID mapping job for partner ENSPs (Ensembl_Protein -> UniProtKB)...")
        job_id = run_mapping_with_retry(need_lookup, "Ensembl_Protein", "UniProtKB")
        mapping = get_mapping_results(job_id)
        print(f"  mapped {len(mapping)} / {len(need_lookup)} partner ENSPs to >=1 UniProt accession")
        for ensp, uids in mapping.items():
            ensp_to_partner_uid[ensp] = uids[0]  # first returned

    out_fields = fieldnames + ["PPI_UniProt_Partners", "PPI_UniProt_Partners_in_Dataframe"]
    n_with_ppi = 0
    n_with_ppi_in_df = 0
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        for row in rows:
            uid = row["uniprot_accession"].strip()
            q_ensp = query_ensp_by_uid.get(uid)
            q_resolved = resolved.get(q_ensp) if q_ensp else None
            partner_ensps = sorted(partner_map.get(q_resolved, [])) if q_resolved else []

            partner_uids = []
            for pe in partner_ensps:
                pu = ensp_to_partner_uid.get(pe)
                if pu:
                    partner_uids.append(pu)
            partner_uids = sorted(set(partner_uids))

            partner_uids_in_df = sorted(set(pu for pu in partner_uids if pu in all_uniprot_ids))

            row["PPI_UniProt_Partners"] = ";".join(partner_uids)
            row["PPI_UniProt_Partners_in_Dataframe"] = ";".join(partner_uids_in_df)
            writer.writerow(row)

            if partner_uids:
                n_with_ppi += 1
            if partner_uids_in_df:
                n_with_ppi_in_df += 1

    print(f"Wrote {len(rows)} rows to {args.out}")
    print(f"  rows with >=1 PPI_UniProt_Partners: {n_with_ppi}")
    print(f"  rows with >=1 PPI_UniProt_Partners_in_Dataframe: {n_with_ppi_in_df}")


if __name__ == "__main__":
    main()
