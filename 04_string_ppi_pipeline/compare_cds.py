import json
import sys
import time
import urllib.request
import urllib.error

JOB_ID = "zv871Jatq2"
UNIPROT_STREAM = f"https://rest.uniprot.org/idmapping/stream/{JOB_ID}?format=tsv"
ENSEMBL_SEQ = "https://rest.ensembl.org/sequence/id?type=cds;content-type=application/json"


def http_get(url):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode()


CHECKPOINT = "/u/home/d/ddcohn/_claude_cds_checkpoint.json"


def load_checkpoint():
    try:
        with open(CHECKPOINT) as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def save_checkpoint(seqs):
    tmp = CHECKPOINT + ".tmp"
    with open(tmp, "w") as f:
        json.dump(seqs, f)
    import os
    os.replace(tmp, CHECKPOINT)


def fetch_cds(enst_bare_ids, batch_size=50):
    seqs = load_checkpoint()
    print(f"  resuming with {len(seqs)} transcripts already cached from checkpoint", flush=True)
    ids = [e for e in enst_bare_ids if e not in seqs]
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        payload = json.dumps({"ids": batch}).encode()
        body = "[]"
        for attempt in range(6):
            try:
                req = urllib.request.Request(
                    ENSEMBL_SEQ, data=payload,
                    headers={"Content-Type": "application/json"},
                )
                with urllib.request.urlopen(req, timeout=60) as resp:
                    body = resp.read().decode()
                break
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(5)
                    continue
                body = e.read().decode()
                print(f"  batch error {e.code}: {body[:200]}", file=sys.stderr)
                body = "[]"
                break
            except Exception as e:
                print(f"  batch exception (attempt {attempt}): {e}", file=sys.stderr)
                time.sleep(5)
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
        if (i // batch_size) % 20 == 0:
            print(f"  ... {i+len(batch)}/{len(ids)} transcripts fetched", flush=True)
        if (i // batch_size) % 100 == 0:
            save_checkpoint(seqs)
        time.sleep(0.3)
    save_checkpoint(seqs)
    return seqs


def main():
    print("Re-fetching UniProt -> Ensembl transcript mapping (all transcripts, not just first)...")
    tsv = http_get(UNIPROT_STREAM)
    lines = tsv.strip().split("\n")
    mapping = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        frm, to = line.split("\t")
        mapping.setdefault(frm, [])
        if to not in mapping[frm]:
            mapping[frm].append(to)

    multi = {k: v for k, v in mapping.items() if len(v) > 1}
    print(f"Proteins with multiple mapped transcripts: {len(multi)}")

    all_bare = set()
    for ensts in multi.values():
        for e in ensts:
            all_bare.add(e.split(".")[0])
    print(f"Distinct transcripts to fetch CDS for: {len(all_bare)}")

    seqs = fetch_cds(all_bare)
    print(f"Got CDS for {len(seqs)} / {len(all_bare)} transcripts")

    all_same = 0
    differ = 0
    incomplete = 0  # couldn't get CDS for all transcripts in the group
    differ_examples = []
    incomplete_examples = []

    for uid, ensts in multi.items():
        bare = [e.split(".")[0] for e in ensts]
        cds_list = [seqs.get(b) for b in bare]
        if any(c is None for c in cds_list):
            incomplete += 1
            if len(incomplete_examples) < 10:
                incomplete_examples.append((uid, ensts))
            continue
        uniq = set(cds_list)
        if len(uniq) == 1:
            all_same += 1
        else:
            differ += 1
            if len(differ_examples) < 10:
                differ_examples.append((uid, ensts, [len(c) for c in cds_list]))

    print()
    print(f"=== Results across {len(multi)} multi-transcript proteins ===")
    print(f"All transcripts share identical CDS: {all_same}")
    print(f"Transcripts have DIFFERENT CDS:      {differ}")
    print(f"Incomplete (missing CDS for >=1):    {incomplete}")
    print()
    print("Examples where CDS differs (uniprot, transcripts, cds lengths):")
    for uid, ensts, lens in differ_examples:
        print(f"  {uid}: {ensts} lengths={lens}")
    print()
    print("Examples with incomplete CDS data:")
    for uid, ensts in incomplete_examples:
        print(f"  {uid}: {ensts}")


if __name__ == "__main__":
    main()
