#!/usr/bin/env python3

#script to run create histograms for individual chambers/data sets
#changed it to not loop over so we can use condor to parallelize
#can run this on its own, but can also run in tandem with submit_efficiency.py on condor
#Example on its own: python3 run_efficiency.py <input_file> <chamber> <output_directory>

import sys
import os
import numpy as np
from array import array

import uproot
import awkward as ak
import ROOT
import coffea
from coffea import processor

#
input_file       = sys.argv[1]   # full EOS path to one ROOT file
chamber          = sys.argv[2]   # e.g. "ME21"
output_directory = sys.argv[3]   # where to write the output ROOT file


def getLZDF(f, nEvents=-1,version="new"): #lazy dataframe with events that have cluster matched to probe muon
    events_raw = uproot.open(f)['MuonSystem']
    df = coffea.processor.LazyDataFrame(events_raw, entrystop=nEvents)
    start,stop = df._branchargs['entry_start'],df._branchargs['entry_stop']
    events = uproot.lazy(df._tree)
    events = events[start:stop]
    return events


chamber_locations = {
    "ME11": {"MinR": 100, "MaxR": 275, "minZ": 580, "maxZ": 632, "HLT_thresh": 500, "branch_names": ["cscRechitClusterNRechitChamberMinus11", "cscRechitClusterNRechitChamberPlus11"]},
    "ME12": {"MinR": 275, "MaxR": 465, "minZ": 668, "maxZ": 725, "HLT_thresh": 200, "branch_names": ["cscRechitClusterNRechitChamberMinus12", "cscRechitClusterNRechitChamberPlus12"]},
    "ME13": {"MinR": 505, "MaxR": 700, "minZ": 668, "maxZ": 724, "HLT_thresh": 200, "branch_names": ["cscRechitClusterNRechitChamberMinus13", "cscRechitClusterNRechitChamberPlus13"]},
    "ME21": {"MinR": 139, "MaxR": 345, "minZ": 789, "maxZ": 850, "HLT_thresh": 500, "branch_names": ["cscRechitClusterNRechitChamberMinus21", "cscRechitClusterNRechitChamberPlus21"]},
    "ME22": {"MinR": 357, "MaxR": 700, "minZ": 791, "maxZ": 850, "HLT_thresh": 200, "branch_names": ["cscRechitClusterNRechitChamberMinus22", "cscRechitClusterNRechitChamberPlus22"]},
    "ME31": {"MinR": 160, "MaxR": 345, "minZ": 915, "maxZ": 970, "HLT_thresh": 500, "branch_names": ["cscRechitClusterNRechitChamberMinus31", "cscRechitClusterNRechitChamberPlus31"]},
    "ME32": {"MinR": 357, "MaxR": 700, "minZ": 911, "maxZ": 970, "HLT_thresh": 200, "branch_names": ["cscRechitClusterNRechitChamberMinus32", "cscRechitClusterNRechitChamberPlus32"]},
    "ME41": {"MinR": 178, "MaxR": 345, "minZ": 1002, "maxZ": 1063, "HLT_thresh": 500, "branch_names": ["cscRechitClusterNRechitChamberMinus41", "cscRechitClusterNRechitChamberPlus41"]},
    "ME42": {"MinR": 357, "MaxR": 700, "minZ": 1002, "maxZ": 1063, "HLT_thresh": 200, "branch_names": ["cscRechitClusterNRechitChamberMinus42", "cscRechitClusterNRechitChamberPlus42"]}
}

nRechits_branches = ["cscRechitClusterNRechitChamberMinus11", "cscRechitClusterNRechitChamberPlus11", "cscRechitClusterNRechitChamberMinus12", "cscRechitClusterNRechitChamberPlus12",
                    "cscRechitClusterNRechitChamberMinus13", "cscRechitClusterNRechitChamberPlus13", "cscRechitClusterNRechitChamberMinus21", "cscRechitClusterNRechitChamberPlus21",
                    "cscRechitClusterNRechitChamberMinus22", "cscRechitClusterNRechitChamberPlus22", "cscRechitClusterNRechitChamberMinus31", "cscRechitClusterNRechitChamberPlus31",
                    "cscRechitClusterNRechitChamberMinus32", "cscRechitClusterNRechitChamberPlus32", "cscRechitClusterNRechitChamberMinus41", "cscRechitClusterNRechitChamberPlus41",
                    "cscRechitClusterNRechitChamberMinus42", "cscRechitClusterNRechitChamberPlus42"]


def noisy_cluster_events(data):
    #return (np.abs(ak.flatten(data.cscRechitClusterPhi))<0.2) & (ak.flatten(data.cscRechitClusterEta < -2)) & (ak.flatten(data.cscRechitClusterZ < -900)) & (data.runNum >367079)
    #return (np.abs(ak.flatten(data.cscRechitClusterPhi))>3) & (ak.flatten(data.cscRechitClusterEta < -1.9))
    return(((ak.flatten(data.cscRechitClusterPhi)>0.4) & (ak.flatten(data.cscRechitClusterPhi)<0.8))|(np.abs(ak.flatten(data.cscRechitClusterPhi))>2.8))


#code to mask for clusters for which majority of events are in designated chamber
#require greater than 90% of hits to be in given chamber as well
def compute_chamber_mask(data, chamber, threshold=0.9, endcap=None):
    final_masks = []
    for x in range(2):
        correct_chamber_branch    = chamber_locations[chamber]["branch_names"][x]
        incorrect_chamber_branches = []

        for chamber_branch in nRechits_branches:
            if chamber_branch != correct_chamber_branch:
                incorrect_chamber_branches.append(chamber_branch)         

        mask = ak.flatten(data[correct_chamber_branch]) > ak.flatten(data[incorrect_chamber_branches[0]])
        for idx in range(1, len(incorrect_chamber_branches)):
            mask = np.logical_and(mask, ak.flatten(data[correct_chamber_branch])>ak.flatten(data[incorrect_chamber_branches[idx]]))
        mask = np.logical_and(mask, (ak.flatten(data[correct_chamber_branch])/ak.flatten(data["cscRechitClusterSize"])>=threshold))
        #mask = np.logical_and(mask, ak.flatten(data["cscRechitClusterNStation10"])==1)
        #code to mask out negative endcap clusters in problematic region
        if x==0 and (chamber=="ME31" or chamber=="ME41"):
            mask = np.logical_and(mask, np.logical_not(noisy_cluster_events(data)))
        #print(ak.count_nonzero(mask))
        final_masks.append(mask)
    if endcap==None:
        return np.logical_or(final_masks[0], final_masks[1])
    else:
        return final_masks[endcap]


def get_efficiency_hists_ROOT(data, chamber, endcap=None):
    total_mask = compute_chamber_mask(data, chamber, 0.9, endcap)

    #COMMENT OUT WHEN NOT DOING NOISE MASK INVERSION
    #total_mask = np.logical_and(total_mask, noisy_cluster_events(data))

    denom = ak.mask(ak.flatten(data.cscRechitClusterSize), total_mask)
    num = ak.mask(denom, ak.mask(np.logical_or(data.L1_SingleMuShower_Nominal, data.L1_SingleMuShower_Tight),total_mask))
    denom = denom[~ak.is_none(denom)]
    num = num[~ak.is_none(num)]

    denom = ak.to_numpy(denom)
    num = ak.to_numpy(num)

    # print(len(num))
    # print(len(denom))
    inner_bins = np.concatenate([np.arange(0, 500, 100), np.arange(500, 1000, 400)])
    outer_bins = np.concatenate([np.arange(0, 200, 50), np.arange(200, 1000, 100)])


    if chamber[-1:]=='1':
        bins = inner_bins
    else:
        bins = outer_bins

    bins_arr = [binEdge for binEdge in bins]

    num_hist = ROOT.TH1F("num", "num", len(bins_arr)-1, array('f', bins_arr))
    denom_hist = ROOT.TH1F("denom", "denom", len(bins_arr)-1, array('f', bins_arr))

    for val in num:
        num_hist.Fill(val)
    for val in denom:
        denom_hist.Fill(val)

    return num_hist, denom_hist


print(f"Processing: {input_file.split('/')[-1]}  chamber: {chamber}")

# derive output filename: <dataset>_<chamber>.root
dataset_name = input_file.split("/")[-1].replace(".root", "")
os.makedirs(output_directory, exist_ok=True)
output_file = os.path.join(output_directory, f"{dataset_name}_{chamber}.root")


input_events = getLZDF(input_file)
input_events = input_events[(input_events.nCscRechitClusters == 1)]
input_events = ak.mask(input_events, (ak.flatten(input_events.cscRechitCluster_matchToProbeMuon))&(ak.flatten(input_events.cscRechitCluster_PassTimeVeto))&(ak.flatten(input_events.cscRechitClusterNRechitChamberPlus11)+ak.flatten(input_events.cscRechitClusterNRechitChamberPlus12)+ak.flatten(input_events.cscRechitClusterNRechitChamberMinus11)+ak.flatten(input_events.cscRechitClusterNRechitChamberMinus12)==0))
input_events = input_events[~ak.is_none(input_events)]


root_file = ROOT.TFile(output_file, "RECREATE")


num_hist, denom_hist = get_efficiency_hists_ROOT(input_events, chamber)
num_hist.Write(f"{chamber}_num")
denom_hist.Write(f"{chamber}_denom")
eff = ROOT.TEfficiency(num_hist, denom_hist)
eff.SetStatisticOption(ROOT.TEfficiency.kFCP)
eff.SetTitle(f"{chamber};Cluster Size;L1 Efficiency")
eff.Write(chamber)
root_file.Close()

print(f"Written: {output_file}")
