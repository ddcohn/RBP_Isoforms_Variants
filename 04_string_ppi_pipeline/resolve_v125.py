import sys
import json
import time
import urllib.request
import urllib.error

SCRATCH = "/private/tmp/claude-501/-Users-danielcohn/dd8ee4f2-df58-45c2-8006-f25553cf99cc/scratchpad"
STRING_RESOLVE_URL = "https://version-12-5.string-db.org/api/tsv/get_string_ids"


def http_get(url, timeout=90):
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def resolve_string_ids(bare_ensps, batch_size=100):
    resolved = {}
    ids = list(bare_ensps)
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        query = "%0d".join(batch)
        url = f"{STRING_RESOLVE_URL}?identifiers={query}&species=9606"
        body = None
        for attempt in range(5):
            try:
                body = http_get(url)
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
            string_id = cols[1]  # e.g. "9606.ENSP00000254108"
            if idx < len(batch):
                resolved[batch[idx]] = string_id
        if (i // batch_size) % 20 == 0:
            print(f"  ... resolve {i+len(batch)}/{len(ids)} done", flush=True)
        time.sleep(0.4)
    return resolved


def main():
    uid_to_ensp = {}
    with open(f"{SCRATCH}/query_ensps.txt") as f:
        for line in f:
            uid, ensp_prefixed = line.rstrip("\n").split("\t")
            bare = ensp_prefixed.split(".", 1)[1] if "." in ensp_prefixed else ensp_prefixed
            uid_to_ensp[uid] = bare

    print(f"Resolving {len(set(uid_to_ensp.values()))} distinct query ENSPs against STRING v12.5...")
    resolved = resolve_string_ids(set(uid_to_ensp.values()))  # bare -> "9606.ENSP..." (prefixed, resolved)
    print(f"Resolved {len(resolved)} / {len(set(uid_to_ensp.values()))}")

    with open(f"{SCRATCH}/resolved_v125.json", "w") as f:
        json.dump(resolved, f)

    # also write a plain list of resolved (prefixed) IDs for the awk local-file matching step
    with open(f"{SCRATCH}/resolved_query_ensps_prefixed.txt", "w") as f:
        for uid, bare in uid_to_ensp.items():
            if bare in resolved:
                f.write(uid + "\t" + resolved[bare] + "\n")

    print("Done.")


if __name__ == "__main__":
    main()
