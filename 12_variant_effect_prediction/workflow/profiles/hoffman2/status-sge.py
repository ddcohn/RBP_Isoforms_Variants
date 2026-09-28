#!/usr/bin/env python3
"""
Snakemake generic-cluster status script for Hoffman2 (UGE/SGE).

Snakemake calls this as `<cluster-status command> <jobid>` and expects
exactly one of "running", "success", "failed" on stdout. This directly
replaces the qstat/qacct polling done by hand throughout this session.

While the job is still in the queue (`qstat`, any state) it's "running".
Once it drops out of qstat, `qacct` has the final exit status -- but
qacct's own accounting file write can lag a couple seconds behind a job
leaving qstat, so a few short retries avoid a false "failed" read in
that gap.
"""
import subprocess
import sys
import time

jobid = sys.argv[-1]


def in_qstat():
    out = subprocess.run(["qstat", "-j", jobid], capture_output=True, text=True)
    return out.returncode == 0


def qacct_exit_status():
    for _ in range(5):
        out = subprocess.run(["qacct", "-j", jobid], capture_output=True, text=True)
        if out.returncode == 0:
            for line in out.stdout.splitlines():
                if line.startswith("exit_status"):
                    return int(line.split()[-1])
        time.sleep(3)
    return None


if in_qstat():
    print("running")
    sys.exit(0)

status = qacct_exit_status()
if status is None:
    # not in qstat and qacct has no record yet -- treat as still running
    # rather than risk a false "failed"; the next poll will resolve it
    print("running")
elif status == 0:
    print("success")
else:
    print("failed")
