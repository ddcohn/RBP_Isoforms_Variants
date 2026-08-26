import sys
import json
import time
import urllib.request
import urllib.parse
import urllib.error

SCRATCH = "/private/tmp/claude-501/-Users-danielcohn/dd8ee4f2-df58-45c2-8006-f25553cf99cc/scratchpad"

STRING_RESOLVE_URL = "https://version-12-5.string-db.org/api/tsv/get_string_ids"
STRING_URL = "https://version-12-5.string-db.org/api/tsv/interaction_partners"
UNIPROT_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_STATUS = "https://rest.uniprot.org/idmapping/status/{}"
UNIPROT_STREAM = "https://rest.uniprot.org/idmapping/stream/{}?format=tsv"

REQUIRED_SCORE = 700
STRING_LIMIT = 100000

CHECKPOINT = f"{SCRATCH}/ppi_api_v125_checkpoint.json"


def http_get(url, timeout=60):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def http_post(url, data_dict, timeout=60):
    data = urllib.parse.urlencode(data_dict).encode()
    req = urllib.request.Request(url, data=data)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def save_checkpoint(resolved, partner_map):
    tmp = CHECKPOINT + ".tmp"
    with open(tmp, "w") as f:
        json.dump({"resolved": resolved, "partner_map": {k: sorted(v) for k, v in partner_map.items()}}, f)
    import os
    os.replace(tmp, CHECKPOINT)


def load_checkpoint():
    try:
        with open(CHECKPOINT) as f:
            d = json.load(f)
        return d["resolved"], {k: set(v) for k, v in d["partner_map"].items()}
    except (FileNotFoundError, json.JSONDecodeError, KeyError):
        return None, None


def resolve_string_ids(bare_ensps, batch_size=100):
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


def wait_for_job(job_id, timeout=1800):
    start = time.time()
    while time.time() - start < timeout:
        try:
            resp = http_get(UNIPROT_STATUS.format(job_id))
            obj = json.loads(resp)
        except Exception as e:
            print(f"  transient status error: {e}", file=sys.stderr)
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
    for attempt in range(max_resubmits):
        resp = http_post(UNIPROT_RUN, {"ids": ",".join(ids), "from": from_db, "to": to_db})
        job_id = json.loads(resp)["jobId"]
        print(f"  job id: {job_id} (attempt {attempt+1}/{max_resubmits})")
        try:
            wait_for_job(job_id)
            return job_id
        except RuntimeError as e:
            print(f"  job errored: {e}; resubmitting...", file=sys.stderr)
            time.sleep(10)
    raise RuntimeError("UniProt mapping job failed after retries")


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
    uid_to_ensp = {}
    ensp_to_uid = {}
    with open(f"{SCRATCH}/query_ensps.txt") as f:
        for line in f:
            uid, ensp_prefixed = line.rstrip("\n").split("\t")
            bare = ensp_prefixed.split(".", 1)[1] if "." in ensp_prefixed else ensp_prefixed
            uid_to_ensp[uid] = bare
            ensp_to_uid.setdefault(bare, uid)

    with open(f"{SCRATCH}/all_uids.txt") as f:
        all_uniprot_ids = set(line.strip() for line in f if line.strip())

    print(f"Query proteins: {len(uid_to_ensp)}; total dataset proteins: {len(all_uniprot_ids)}")

    # reuse the already-computed resolution (shared with the local-file method for
    # a true apples-to-apples comparison) instead of re-resolving here
    with open(f"{SCRATCH}/resolved_v125.json") as f:
        resolved_prefixed = json.load(f)
    # fetch_string_partners/interaction_partners wants bare IDs (no "9606." prefix),
    # unlike get_string_ids/resolve which returns prefixed IDs
    resolved = {k: (v.split(".", 1)[1] if "." in v else v) for k, v in resolved_prefixed.items()}
    print(f"Loaded {len(resolved)} pre-resolved query IDs")

    _, partner_map = load_checkpoint()
    if partner_map is None:
        print("Querying STRING v12.5 for interaction partners...")
        partner_map = fetch_string_partners(set(resolved.values()))
        print(f"Got partner data for {len(partner_map)} resolved query proteins")

        save_checkpoint(resolved, partner_map)
        print("  (checkpoint saved)")
    else:
        print(f"Resuming from checkpoint: {len(partner_map)} with partner data")

    all_partner_ensps = set()
    for s in partner_map.values():
        all_partner_ensps.update(s)
    print(f"Distinct partner ENSPs (API v12.5): {len(all_partner_ensps)}")

    need_lookup = sorted(e for e in all_partner_ensps if e not in ensp_to_uid)
    print(f"Partner ENSPs needing UniProt lookup: {len(need_lookup)}")

    ensp_to_partner_uid = dict(ensp_to_uid)
    if need_lookup:
        print("Submitting UniProt ID mapping job (Ensembl_Protein -> UniProtKB)...")
        job_id = run_mapping_with_retry(need_lookup, "Ensembl_Protein", "UniProtKB")
        mapping = get_mapping_results(job_id)
        print(f"  mapped {len(mapping)} / {len(need_lookup)} partner ENSPs")
        for ensp, uids in mapping.items():
            ensp_to_partner_uid[ensp] = uids[0]

    result = {}
    for uid, q_ensp in uid_to_ensp.items():
        q_resolved = resolved.get(q_ensp)
        partner_ensps = partner_map.get(q_resolved, set()) if q_resolved else set()
        partner_uids = sorted(set(
            ensp_to_partner_uid[pe] for pe in partner_ensps if pe in ensp_to_partner_uid
        ))
        in_df = sorted(set(pu for pu in partner_uids if pu in all_uniprot_ids))
        result[uid] = (partner_uids, in_df)

    with open(f"{SCRATCH}/ppi_api_v125_result.json", "w") as f:
        json.dump(result, f)

    n_with = sum(1 for v in result.values() if v[0])
    n_with_df = sum(1 for v in result.values() if v[1])
    print(f"Done. {n_with} proteins with >=1 partner, {n_with_df} with >=1 in-dataframe partner")


if __name__ == "__main__":
    main()
