#!/bin/bash

root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1500pre3_default", "Phase2_Spring24_1500pre3_default", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1500pre3/20250331/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_rerun_1500pre3_default_20250331/250331_040637/0000/*.root"})' >&Phase2_Spring24_1500pre3_default.log&
root -l -b -q 'HLTBDTAnalyzer_binary.C("Phase2_Spring24_1500pre3_ticl_v5", "Phase2_Spring24_1500pre3_ticl_v5", {"/eos/cms/store/group/phys_muon/ec/HLT/MuonHLTPhase2_cmssw1500pre3/20250331/DYToLL_M-50_TuneCP5_14TeV-pythia8/crab_DYToLL_M50_Phase2Spring24_hlt_muon_mc_rerun_1500pre3_ticl_v5_20250331/250331_040621/0000/*.root"})' >&Phase2_Spring24_1500pre3_ticl_v5.log&
