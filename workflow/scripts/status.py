#!/usr/bin/env python

import subprocess
import sys
# From Peter Maxwell @ NeSI 11-2023
#
jobid = sys.argv[1]

try:
    output = subprocess.check_output(["squeue", "-j", str(jobid), "--Format=state", "--states=ALL", "--noheader"], 
            universal_newlines=True).strip().split()[0]
except subprocess.CalledProcessError:
    output = subprocess.check_output(["sacct", "-j", str(jobid), "--format=state", "--noheader"], 
            universal_newlines=True).strip().split()[0]

failure_states = {"FAILED", "TIMEOUT", "OUT_OF_MEMORY", "DEADLINE", "CANCELLED"}

if output == "COMPLETED":
    print("success")
elif output in failure_states:
    print("failed")
else:
    print("running")
