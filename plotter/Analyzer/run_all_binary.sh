#!/bin/bash

#root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1600pre3_default", "Phase2_Spring24_1600pre3_default", {"/eos/cms/store/group/phys_muon/sblancof/HLT/PhaseII/MuonHLTPhase2_cmssw1600pre3/20251224/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_default_20251224/251224_085621/0000/*.root"})' >&Phase2_Spring24_1500pre3_default.log&

#root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1600pre3_seedless", "Phase2_Spring24_1600pre3_seedless", {"/eos/cms/store/group/phys_muon/sblancof/HLT/PhaseII/MuonHLTPhase2_cmssw1600pre3/20251229/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_seedless_20251229/251229_114911/0000/*.root"})' >&Phase2_Spring24_1500pre3_seedless.log&

# Updated OI track selection
#root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1600pre3_seedless", "Phase2_Spring24_1600pre3_seedless", {"/eos/cms/store/group/phys_muon/sblancof/HLT/PhaseII/MuonHLTPhase2_cmssw1600pre3_updatedOIFilter/20260112/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_seedless_20260112/260112_094727/0000/*.root"})' >&Phase2_Spring24_1500pre3_seedless.log&

########## Different IO Muon reconstruction strategies

root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1600pre3_generalTracks", "Phase2_Spring24_1600pre3_generalTracks", {"/eos/cms/store/group/phys_muon/sblancof/HLT/PhaseII/MuonHLTPhase2_cmssw1600pre3_updatedOIFilter_generalTracks/20260113/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_generalTracks_20260113_generalTracks/260113_105310/0000/*.root"})' >&Phase2_Spring24_1500pre3_generalTracks.log&

root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1600pre3_offline", "Phase2_Spring24_1600pre3_offline", {"/eos/cms/store/group/phys_muon/sblancof/HLT/PhaseII/MuonHLTPhase2_cmssw1600pre3_updatedOIFilter_offline/20260113/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_offline_20260113_offline/260113_105340/0000/*.root"})' >&Phase2_Spring24_1500pre3_offline.log&

root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1600pre3_singleIterTracking", "Phase2_Spring24_1600pre3_singleIterTracking", {"/eos/cms/store/group/phys_muon/sblancof/HLT/PhaseII/MuonHLTPhase2_cmssw1600pre3_updatedOIFilter_singleIterTracking/20260113/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_singleIterTracking_20260113_singleIterTracking/260113_105608/0000/*.root"})' >&Phase2_Spring24_1500pre3_singleIterTracking.log&

root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1600pre3_singleIterTracking_generalTracks", "Phase2_Spring24_1600pre3_singleIterTracking_generalTracks", {"/eos/cms/store/group/phys_muon/sblancof/HLT/PhaseII/MuonHLTPhase2_cmssw1600pre3_updatedOIFilter_singleIterTracking_generalTracks/20260113/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_singleIterTracking_generalTracks_20260113_singleIterTracking_generalTracks/260113_105637/0000/*.root"})' >&Phase2_Spring24_1500pre3_generalTracks.log&

root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1600pre3_singleIterTracking_offline", "Phase2_Spring24_1600pre3_singleIterTracking_offline", {"/eos/cms/store/group/phys_muon/sblancof/HLT/PhaseII/MuonHLTPhase2_cmssw1600pre3_updatedOIFilter_singleIterTracking_offline/20260113/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_singleIterTracking_offline_20260113_singleIterTracking_offline/260113_105706/0000/*.root"})' >&Phase2_Spring24_1500pre3_singleIterTracking_offline.log&
