#!/usr/bin/env python3
import os
import sys
import subprocess


rootdir = "2024_retest4"
run_kind = "Muon0-Run2024B-PromptReco-v1-AOD"
eos_dir = f"/store/group/lpcmds/ctidmore/test/{rootdir}/{run_kind}"
input_dir = f"root://cmseos.fnal.gov/{eos_dir}"

output_dir = "/uscms/home/ctidmore/nobackup/CMSSW_10_6_30/src/run3_llp_analyzer/Plotting_Scripts/Improve_Eff_Scripts/merged_roots_analyzer"
#output_dir = sys.argv[1]



ls_result = subprocess.run(["xrdfs", "root://cmseos.fnal.gov", "ls", eos_dir], capture_output=True, text=True)
if ls_result.returncode != 0:
    print("ERROR: could not list EOS directory")
    print("===================================================================")
    print(ls_result.stderr)
    print("===================================================================")
    exit(1)


root_files = [f"root://cmseos.fnal.gov/{f.strip()}" for f in ls_result.stdout.splitlines() if f.strip().endswith(".root")]

merged_file = f"{output_dir}/{run_kind}_merged.root"

cmd = ["hadd", "-f", merged_file] + root_files

#result = subprocess.run(cmd, capture_output=True, text=True)
result = subprocess.run(cmd)

if result.returncode != 0:
    print("ERROR: hadd failed")
    print("")
    print("===================================================================")
    print(result.stderr)
    print("===================================================================")
else:
    print(f"Done. {len(root_files)} root files merged successfully.")