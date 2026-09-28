#!/usr/bin/env python3
"""
Snakemake generic-cluster submit script for Hoffman2 (UGE/SGE).

Snakemake calls this as `<cluster command> <jobscript>` and expects
exactly the numeric job ID printed to stdout. Resource values come from
each rule's `resources:` block (see ../../Snakefile) via the JSON
"properties" comment Snakemake embeds in the generated jobscript --
this is the standard way a generic (non-DRMAA) cluster wrapper reads
per-job resources in Snakemake 7.

Converts Snakemake's `runtime` resource (minutes) to SGE's h_rt
(H:MM:SS) and `mem_mb` to h_data; `sge_extra` is appended verbatim so
individual rules can request e.g. a GPU queue or `highp` without this
script needing to know about every queue.
"""
import json
import os
import re
import subprocess
import sys

jobscript = sys.argv[-1]

with open(jobscript) as f:
    text = f.read()
m = re.search(r"^# properties = (.*)$", text, re.MULTILINE)
props = json.loads(m.group(1)) if m else {}

resources = props.get("resources", {})
rule_name = props.get("rule", "job")
wildcards = props.get("wildcards", {})
tag = "_".join(str(v) for v in wildcards.values()) or "na"

mem_mb = int(resources.get("mem_mb", 8000))
runtime_min = int(resources.get("runtime", 60))
sge_extra = resources.get("sge_extra", "") or ""

hours, minutes = divmod(runtime_min, 60)
h_rt = f"{hours}:{minutes:02d}:00"

log_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "..", "logs")
log_dir = os.path.abspath(log_dir)
os.makedirs(log_dir, exist_ok=True)
log_path = os.path.join(log_dir, f"{rule_name}.{tag}.log")

cmd = [
    # -b y: run the given command directly rather than treating it as an
    # SGE-native script (its own shebang isn't a login shell, so `module`
    # would be undefined otherwise -- Hoffman2/Lmod defines it only in
    # login shells; the fix used throughout this whole pipeline).
    "qsub", "-b", "y", "-cwd", "-j", "y", "-o", log_path,
    "-N", f"smk_{rule_name}",
    "-l", f"h_data={mem_mb}M,h_rt={h_rt}",
]
if sge_extra:
    cmd += sge_extra.split()
cmd += ["/bin/bash", "-l", jobscript]

out = subprocess.check_output(cmd, text=True)
# qsub prints e.g. Your job 12345 ("smk_rule") has been submitted
match = re.search(r"\d+", out)
if not match:
    sys.stderr.write(f"Could not parse job id from qsub output: {out}\n")
    sys.exit(1)
print(match.group(0))
