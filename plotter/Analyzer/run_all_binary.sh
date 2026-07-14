#!/bin/bash

root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1500pre3_default", "Phase2_Spring24_1500pre3_default", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1500pre3/20250331/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_rerun_1500pre3_default_20250331/250331_040637/0000/*.root"})' >&Phase2_Spring24_1500pre3_default.log&
root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1500pre3_ticl_v5", "Phase2_Spring24_1500pre3_ticl_v5", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1500pre3/20250331/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_rerun_1500pre3_ticl_v5_20250331/250331_040621/0000/*.root"})' >&Phase2_Spring24_1500pre3_ticl_v5.log&


root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1510pre4_OIFirst", "Phase2_Spring24_1510pre4_OIFirst", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1510pre4/20251110/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_OIFirst_dump_20251110/251110_220218/0000/*.root"})' >&Phase2_Spring24_1500pre3_ticl_v5.log&

root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1510pre4_OIFirst", "OIFirst_default", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1510pre4/20251110/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_OIFirst_dump_20251110/251110_220218/0000/*.root"})' >&OIDefault.log&
root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1510pre4_OIFirst_OIFromL1", "OIFirst_FromL1", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1510pre4/20251117/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_OIFirst_OIFromL1_dump_20251117/251117_163637/0000/*.root"})' >&OIFromL1.log&
root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1510pre4_OIFirst_OIFromL1Tk", "OIFirst_FromL1Tk", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1510pre4/20251117/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_OIFirst_OIFromL1Tk_dump_20251117/251117_163624/0000/*.root"})' >&OIFromL1Tk.log&



root -l -b -q 'HLTBDTAnalyzer_binary.C("CMSSW1600pre3_default", "CMSSW1600pre3_default", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1600pre3/20251212/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_20251212/251212_185846/0000/*.root"})' >&log_default.log&
root -l -b -q 'HLTBDTAnalyzer_binary.C("CMSSW1600pre3_newTracking", "CMSSW1600pre3_newTracking", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1600pre3/20251212/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_newTRKbaseline_20251212/251212_185856/0000/*.root"})' >&log_new.log&





