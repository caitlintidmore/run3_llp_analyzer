#!/usr/bin/env python3

#runs condor jobs by parallelizing over datasets and chambers
#How it works:
#run the below command, but imput a log name
#python3 submit_efficiency.py <log_name> 

import os
import sys
import subprocess

log_name=sys.argv[1]
output_dir="root://cmseos.fnal.gov/store/group/lpcmds/ctidmore/all2024"
#output_dir = sys.argv[2]   # EOS base path for output ROOT files

log_dir    = f"log_{log_name}"
submit_dir = f"submit_{log_name}"
os.makedirs(log_dir,    exist_ok=True)
os.makedirs(submit_dir, exist_ok=True)

#HOME       = os.getenv('HOME')
#CMSSW_BASE = os.getenv('CMSSW_BASE')

HOME="/uscms/home/ctidmore/nobackup/CMSSW_10_6_30/src/run3_llp_analyzer"
CMSSW_BASE="/uscms/home/ctidmore/nobackup/CMSSW_10_6_30"


# #path to the file I made
# list_2024 = ["root://cmseos.fnal.gov//store/group/lpcmds/ctidmore/test/2024_retest4/Muon0-Run2024B-PromptReco-v1-AOD_merged/normalized/Muon0-Run2024B-PromptReco-v1-AOD_merged_goodLumi.root"]

#big list that has all the 2024 files Alex made
list_2024 = [
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024B-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024C-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024D-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024E-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024E-PromptReco-v2-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024F-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024G-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024H-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024I-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon0-Run2024I-PromptReco-v2-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024B-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024C-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024D-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024E-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024E-PromptReco-v2-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024F-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024G-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024H-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024I-PromptReco-v1-AOD_goodLumi.root",
    "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024I-PromptReco-v2-AOD_goodLumi.root",
]

#test chambers
# chambers = ["ME13"]

#Full list of chambers Alex used
chambers = ["ME13", "ME21", "ME22", "ME31", "ME32", "ME41", "ME42"]

#this is where i hope it submits one job per set per chamber
total_jobs = 0

for input_file in list_2024:
    dataset_name = input_file.split("/")[-1].replace(".root", "")

    for chamber in chambers:
        job_name = f"{dataset_name}_{chamber}"
        jdl_path = f"{submit_dir}/{job_name}.jdl"

        with open(jdl_path, "w") as jdl:
            jdl.write("Universe = vanilla\n")
            jdl.write("Executable = run_efficiency.sh\n")
            jdl.write(f"Arguments = {input_file} {chamber} {output_dir} {CMSSW_BASE} {HOME}\n")
            jdl.write(f"Log    = {log_dir}/{job_name}_$(Cluster).$(Process).log\n")
            jdl.write(f"Output = {log_dir}/{job_name}_$(Cluster).$(Process).out\n")
            jdl.write(f"Error  = {log_dir}/{job_name}_$(Cluster).$(Process).err\n")
            jdl.write('+JobFlavour = "tomorrow"\n')   
            jdl.write("RequestMemory = 4096\n")
            jdl.write("RequestCpus = 1\n")
            jdl.write("+RunAsOwner = True\n")
            jdl.write("+InteractiveUser = true\n")
            jdl.write('+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7"\n')         #this is el7; needed for CMSSW_10_6_30
            jdl.write("+SingularityBindCVMFS = True\n")
            jdl.write("run_as_owner = True\n")
            jdl.write("should_transfer_files = YES\n")
            jdl.write("when_to_transfer_output = ON_EXIT_OR_EVICT\n")
            jdl.write("Transfer_Input_Files = run_efficiency.py\n")
            jdl.write("Queue 1\n")

        os.system(f"condor_submit {jdl_path} --batch-name {job_name}")
        total_jobs += 1

print(f"\nSubmitted {total_jobs} jobs")
