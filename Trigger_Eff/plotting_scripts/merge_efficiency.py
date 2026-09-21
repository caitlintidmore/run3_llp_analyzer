#!/usr/bin/env python3
import os
import sys
import subprocess

#script used to merge the chamgers for each dataset into single root
#uses hadd


#option to specify output dir in command line
#output_dir = sys.argv[1]

#hardcoded input dir for now
#variable for finding the unmerged root files
input_dir = "root://cmseos.fnal.gov//store/group/lpcmds/ctidmore/all2024"

#variable to store the newly merged root file
output_dir = "/uscms/home/ctidmore/nobackup/CMSSW_10_6_30/src/run3_llp_analyzer/Plotting_Scripts/Improve_Eff_Scripts/all_merged_roots"

chambers = ["ME13", "ME21", "ME22", "ME31", "ME32", "ME41", "ME42"]

# all 20 dataset base names (without .root)
datasets = [
    "Muon0-Run2024B-PromptReco-v1-AOD_goodLumi",
    "Muon0-Run2024C-PromptReco-v1-AOD_goodLumi",
    "Muon0-Run2024D-PromptReco-v1-AOD_goodLumi",
    "Muon0-Run2024E-PromptReco-v1-AOD_goodLumi",
    "Muon0-Run2024E-PromptReco-v2-AOD_goodLumi",
    "Muon0-Run2024F-PromptReco-v1-AOD_goodLumi",
    "Muon0-Run2024G-PromptReco-v1-AOD_goodLumi",
    "Muon0-Run2024H-PromptReco-v1-AOD_goodLumi",
    "Muon0-Run2024I-PromptReco-v1-AOD_goodLumi",
    "Muon0-Run2024I-PromptReco-v2-AOD_goodLumi",
    "Muon1-Run2024B-PromptReco-v1-AOD_goodLumi",
    "Muon1-Run2024C-PromptReco-v1-AOD_goodLumi",
    "Muon1-Run2024D-PromptReco-v1-AOD_goodLumi",
    "Muon1-Run2024E-PromptReco-v1-AOD_goodLumi",
    "Muon1-Run2024E-PromptReco-v2-AOD_goodLumi",
    "Muon1-Run2024F-PromptReco-v1-AOD_goodLumi",
    "Muon1-Run2024G-PromptReco-v1-AOD_goodLumi",
    "Muon1-Run2024H-PromptReco-v1-AOD_goodLumi",
    "Muon1-Run2024I-PromptReco-v1-AOD_goodLumi",
    "Muon1-Run2024I-PromptReco-v2-AOD_goodLumi",
]

failed = []

for dataset in datasets:
    # build list of the 7 per-chamber input files
    input_files = [
        f"{input_dir}/{dataset}_{chamber}.root"
        for chamber in chambers
    ]

    # # check all 7 exist before trying to merge
    # missing = [f for f in input_files if not os.path.exists(f)]
    # if missing:
    #     print(f"SKIPPING {dataset} — missing files:")
    #     for m in missing:
    #         print(f"  {m}")
    #     failed.append(dataset)
    #     continue

    merged_file = f"{output_dir}/{dataset}.root"

    cmd = ["hadd", "-f", merged_file] + input_files
    print(f"Merging {dataset} ...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: hadd failed for {dataset}")
        print(result.stderr)
        failed.append(dataset)
    else:
        print(f"  -> {merged_file}")
        # optionally remove the per-chamber files now that they're merged
        # comment these lines out if you want to keep them
        # for f in input_files:
        #     os.remove(f)

print("\n========================================")
print(f"Done. {len(datasets) - len(failed)}/{len(datasets)} datasets merged successfully.")
if failed:
    print(f"Failed or skipped: {failed}")
