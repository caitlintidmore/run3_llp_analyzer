#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../")
import mplhep as hep
import pickle
import glob
import ROOT
import coffea
import awkward as ak
from coffea import hist, processor
from coffea.nanoevents.methods import candidate
from coffea.nanoevents.methods import vector


# In[2]:





# In[3]:


ak.behavior.update(candidate.behavior)

def getLZDF(f,nEvents=-1,version="new"): #lazy dataframe with events that have cluster matched to probe muon
    events_raw = uproot.open(f)['MuonSystem']
    df = coffea.processor.LazyDataFrame(events_raw,entrystop=nEvents)
    start,stop = df._branchargs['entry_start'],df._branchargs['entry_stop']
    events = uproot.lazy(df._tree)
    events = events[start:stop]
    return events






# In[6]:


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
           "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HMT_L1_Eff_Output/2024_Merged/Muon1-Run2024I-PromptReco-v2-AOD_goodLumi.root"]




# In[7]:


data_events_2024 = [getLZDF(file) for file in list_2024]
events_2024 = ak.concatenate(data_events_2024, axis=0)
#events_2024_probeMuonMatched = events_2024[events_2024.cscRechitCluster_matchToProbeMuon]

events_2024_probeMuonMatched = ak.mask(events_2024, (ak.flatten(events_2024.cscRechitCluster_matchToProbeMuon))&(ak.flatten(events_2024.cscRechitCluster_PassTimeVeto))&(ak.flatten(events_2024.cscRechitClusterNRechitChamberPlus11)+ak.flatten(events_2024.cscRechitClusterNRechitChamberPlus12)+ak.flatten(events_2024.cscRechitClusterNRechitChamberMinus11)+ak.flatten(events_2024.cscRechitClusterNRechitChamberMinus12)==0))
events_2024_probeMuonMatched = events_2024_probeMuonMatched[~ak.is_none(events_2024_probeMuonMatched)]






# In[10]:


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


# In[11]:


nRechits_branches = ["cscRechitClusterNRechitChamberMinus11", "cscRechitClusterNRechitChamberPlus11", "cscRechitClusterNRechitChamberMinus12", "cscRechitClusterNRechitChamberPlus12",
                    "cscRechitClusterNRechitChamberMinus13", "cscRechitClusterNRechitChamberPlus13", "cscRechitClusterNRechitChamberMinus21", "cscRechitClusterNRechitChamberPlus21",
                    "cscRechitClusterNRechitChamberMinus22", "cscRechitClusterNRechitChamberPlus22", "cscRechitClusterNRechitChamberMinus31", "cscRechitClusterNRechitChamberPlus31",
                    "cscRechitClusterNRechitChamberMinus32", "cscRechitClusterNRechitChamberPlus32", "cscRechitClusterNRechitChamberMinus41", "cscRechitClusterNRechitChamberPlus41",
                    "cscRechitClusterNRechitChamberMinus42", "cscRechitClusterNRechitChamberPlus42"]





chambers = ["ME11", "ME12", "ME13", "ME21", "ME22", "ME31", "ME32", "ME41", "ME42"]

chamber_bins = [hist.Bin("cluster_size", "Cluster Size", np.concatenate([np.arange(0, 200, 20), np.arange(200, 1000, 100)])), 
                hist.Bin("cluster_size", "Cluster Size", np.concatenate([np.arange(0, 200, 20), np.arange(200, 1000, 100)])), 
                hist.Bin("cluster_size", "Cluster Size", np.concatenate([np.arange(0, 200, 20), np.arange(200, 1000, 100)])), 
                hist.Bin("cluster_size", "Cluster Size", np.concatenate([np.arange(0, 500, 100), np.arange(500, 1000, 400)])), 
                hist.Bin("cluster_size", "Cluster Size", np.concatenate([np.arange(0, 200, 50), np.arange(200, 1000, 100)])), 
                hist.Bin("cluster_size", "Cluster Size", np.concatenate([np.arange(0, 500, 100), np.arange(500, 1000, 400)])), 
                hist.Bin("cluster_size", "Cluster Size", np.concatenate([np.arange(0, 200, 50), np.arange(200, 1000, 100)])), 
                hist.Bin("cluster_size", "Cluster Size", np.concatenate([np.arange(0, 500, 100), np.arange(500, 1000, 400)])), 
                hist.Bin("cluster_size", "Cluster Size", np.concatenate([np.arange(0, 200, 50), np.arange(200, 1000, 100)]))]


# In[14]:


def noisy_cluster_events(data):
    #return (np.abs(ak.flatten(data.cscRechitClusterPhi))<0.2) & (ak.flatten(data.cscRechitClusterEta < -2)) & (ak.flatten(data.cscRechitClusterZ < -900)) & (data.runNum >367079)
    #return (np.abs(ak.flatten(data.cscRechitClusterPhi))>3) & (ak.flatten(data.cscRechitClusterEta < -1.9))
    return(((ak.flatten(data.cscRechitClusterPhi)>0.4) & (ak.flatten(data.cscRechitClusterPhi)<0.8))|(np.abs(ak.flatten(data.cscRechitClusterPhi))>2.8))


# In[15]:


#code to mask for clusters for which majority of events are in designated chamber
#require greater than 90% of hits to be in given chamber as well
def compute_chamber_mask(data, chamber, chamber_int, threshold, endcap=None):
    nRechits_byChamber = {}
    final_masks = []
    for x in range(2):
        correct_chamber_branch = chamber_locations[chamber]["branch_names"][x]
        incorrect_chamber_branches = []
        for chamber_branch in nRechits_branches:
            if chamber_branch != correct_chamber_branch:
                incorrect_chamber_branches.append(chamber_branch)         
        mask = ak.flatten(data[correct_chamber_branch])>ak.flatten(data[incorrect_chamber_branches[0]])
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


# In[16]:


def get_efficiency_hists_ROOT(data, chamber, bins, endcap=None):
    total_mask = compute_chamber_mask(data, chamber, chamber_int, 0.9, endcap)

    #COMMENT OUT WHEN NOT DOING NOISE MASK INVERSION
    #total_mask = np.logical_and(total_mask, noisy_cluster_events(data))

    denom = ak.mask(ak.flatten(data.cscRechitClusterSize), total_mask)
    num = ak.mask(denom, ak.mask(np.logical_or(data.L1_SingleMuShower_Nominal, data.L1_SingleMuShower_Tight),total_mask))
    denom = denom[~ak.is_none(denom)]
    num = num[~ak.is_none(num)]

    denom = ak.to_numpy(denom)
    num = ak.to_numpy(num)

    print(len(num))
    print(len(denom))
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


# In[17]:


import os
from array import array
chamber_int = 0
#sample = 'all'
output_directory = "2024_trigger_efficiencies_latertest"
os.makedirs(output_directory, exist_ok=True)
for file in list_2024:   
    print(file.split("/")[-1])
    output_file = output_directory+"/"+file.split("/")[-1]
    print(output_file)
    root_file = ROOT.TFile(output_file, "RECREATE")
    input_events = getLZDF(file)
    input_events = input_events[(input_events.nCscRechitClusters==1)]
    input_events = ak.mask(input_events, (ak.flatten(input_events.cscRechitCluster_matchToProbeMuon))&(ak.flatten(input_events.cscRechitCluster_PassTimeVeto))&(ak.flatten(input_events.cscRechitClusterNRechitChamberPlus11)+ak.flatten(input_events.cscRechitClusterNRechitChamberPlus12)+ak.flatten(input_events.cscRechitClusterNRechitChamberMinus11)+ak.flatten(input_events.cscRechitClusterNRechitChamberMinus12)==0))
    input_events = input_events[~ak.is_none(input_events)]
    denom_counts=0


    for chamber_idx, chamber in enumerate(chambers):
        if chamber in ["ME11", "ME12"]:                                 #skip over ME11 and ME12
            continue

        print(chamber)                                                  #print out the chamber name to know where we are in the loop
        
        bins = chamber_bins[chamber_idx]
        num_hist, denom_hist = get_efficiency_hists_ROOT(input_events, chamber, bins)
        #print(type(num_hist))
        denom_counts+=denom_hist.GetEntries()
        num_hist.Write(f"{chamber}_num")
        denom_hist.Write(f"{chamber}_denom")
        eff = ROOT.TEfficiency(num_hist, denom_hist)
        eff.SetStatisticOption(ROOT.TEfficiency.kFCP) 

        #current_plot = dump_chamber_efficiency_plot(chamber, num_hist, denom_hist, bins)
        eff.SetTitle(f"{chamber};Cluster Size;L1 Efficiency")
        eff.Write(chamber)

    root_file.Close()

