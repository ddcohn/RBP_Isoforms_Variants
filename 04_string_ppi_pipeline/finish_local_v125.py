import csv
import sys
import json
import time
import urllib.request
import urllib.parse
import urllib.error

csv.field_size_limit(sys.maxsize)

SCRATCH = "/private/tmp/claude-501/-Users-danielcohn/dd8ee4f2-df58-45c2-8006-f25553cf99cc/scratchpad"
TABLE = "/private/tmp/claude-501/-Users-danielcohn/dd8ee4f2-df58-45c2-8006-f25553cf99cc/scratchpad/table_working.csv"

UNIPROT_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_STATUS = "https://rest.uniprot.org/idmapping/status/{}"
UNIPROT_STREAM = "https://rest.uniprot.org/idmapping/stream/{}?format=tsv"


def http_get(url, timeout=60):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def http_post(url, data_dict, timeout=60):
    data = urllib.parse.urlencode(data_dict).encode()
    req = urllib.request.Request(url, data=data)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


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
    # load uid <-> query bare ensp mapping
    uid_to_ensp = {}
    ensp_to_uid = {}
    with open(f"{SCRATCH}/query_ensps.txt") as f:
        for line in f:
            uid, ensp = line.rstrip("\n").split("\t")
            uid_to_ensp[uid] = ensp
            ensp_to_uid.setdefault(ensp, uid)

    # load query -> partner pairs from local file scan
    partner_map = {}
    with open(f"{SCRATCH}/pairs_v125_dedup.txt") as f:
        for line in f:
            q, p = line.rstrip("\n").split(" ")
            partner_map.setdefault(q, set()).add(p)

    all_partner_ensps = set()
    for s in partner_map.values():
        all_partner_ensps.update(s)
    print(f"Distinct partner ENSPs (local v12.5): {len(all_partner_ensps)}")

    need_lookup = sorted(e for e in all_partner_ensps if e not in ensp_to_uid)
    print(f"Partner ENSPs needing UniProt lookup: {len(need_lookup)}")

    ensp_to_partner_uid = dict(ensp_to_uid)
    if need_lookup:
        # UniProt's Ensembl_Protein mapping wants bare IDs (no "9606." species prefix,
        # which the local STRING file carries) -- strip it before submitting, then
        # reattach the prefix when storing results so keys match partner_map's format.
        need_lookup_bare = [e.split(".", 1)[1] if "." in e else e for e in need_lookup]
        bare_to_prefixed = dict(zip(need_lookup_bare, need_lookup))
        print("Submitting UniProt ID mapping job (Ensembl_Protein -> UniProtKB)...")
        job_id = run_mapping_with_retry(need_lookup_bare, "Ensembl_Protein", "UniProtKB")
        mapping = get_mapping_results(job_id)
        print(f"  mapped {len(mapping)} / {len(need_lookup_bare)} partner ENSPs")
        for bare_ensp, uids in mapping.items():
            prefixed = bare_to_prefixed.get(bare_ensp, bare_ensp)
            ensp_to_partner_uid[prefixed] = uids[0]

    # save intermediate result: uid -> sorted partner uid list
    all_uniprot_ids = set(uid_to_ensp.keys())
    result = {}
    for uid, q_ensp in uid_to_ensp.items():
        partner_ensps = partner_map.get(q_ensp, set())
        partner_uids = sorted(set(
            ensp_to_partner_uid[pe] for pe in partner_ensps if pe in ensp_to_partner_uid
        ))
        in_df = sorted(set(pu for pu in partner_uids if pu in all_uniprot_ids))
        result[uid] = (partner_uids, in_df)

    with open(f"{SCRATCH}/ppi_local_v125_result.json", "w") as f:
        json.dump(result, f)

    n_with = sum(1 for v in result.values() if v[0])
    n_with_df = sum(1 for v in result.values() if v[1])
    print(f"Done. {n_with} proteins with >=1 partner, {n_with_df} with >=1 in-dataframe partner")


if __name__ == "__main__":
    main()
