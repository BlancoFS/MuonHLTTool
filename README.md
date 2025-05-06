# MuonHLT Ntupler

## Run3 140X Recipe
```
export SCRAM_ARCH=el8_amd64_gcc12
cmsrel CMSSW_14_0_9
cd CMSSW_14_0_9/src
cmsenv
git cms-init

git cms-addpkg HLTrigger/Configuration
git clone -b Run2024 https://github.com/wonpoint4/MuonHLTForRun3.git HLTrigger/Configuration/python/MuonHLTForRun3

## Data (Efficiency) - no L1 emul
hltGetConfiguration /dev/CMSSW_14_0_0/GRun/V148 \
 --process MYHLT \
 --data --globaltag 140X_dataRun3_HLT_for2024TSGStudies_v1 \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_CascadeMu100_v*,\
HLT_HighPtTkMu100_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/data/Run2024C/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/379/774/00000/80a17af7-6739-4ca5-90af-923442f2a321.root \
 --eras Run3 \
 --max-events -1 \
 --full --offline --no-output >hlt_muon_data.py

hltGetConfiguration /users/missirol/test/dev/CMSSW_14_0_0/CMSHLT_3224/Test03/GRun/V2 \
 --process MYHLT \
 --data --globaltag 140X_dataRun3_HLT_for2024TSGStudies_v1 \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_CascadeMu100_v*,\
HLT_HighPtTkMu100_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/data/Run2024C/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/379/774/00000/80a17af7-6739-4ca5-90af-923442f2a321.root \
 --eras Run3 \
 --max-events -1 \
 --full --offline --no-output >hlt_muon_data_CSC.py

# With BDT in Iter3FromL1
hltGetConfiguration /dev/CMSSW_14_0_0/GRun \
 --process MYHLT \
 --data --globaltag 140X_dataRun3_HLT_for2024TSGStudies_v1 \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_CascadeMu100_v*,\
HLT_HighPtTkMu100_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/data/Run2024C/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/379/774/00000/80a17af7-6739-4ca5-90af-923442f2a321.root \
 --eras Run3 --l1-emulator uGT --l1 L1Menu_Collisions2024_v1_1_0_xml \
 --max-events -1 \
 --customise \
HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeIOSeedingPatatrack_withIter3FromL1_wp00,\
HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeIOSeedingPatatrack_withIter3FromL1_wp04 \
 --full --offline --no-output >hlt_muon_data_BDT_wp04.py


## Data (Timing)
### https://twiki.cern.ch/twiki/bin/viewauth/CMS/TriggerStudiesTiming
### Output will be at /eos/cms/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/timing_server_results/wjun/
hltGetConfiguration /dev/CMSSW_14_0_0/GRun/V107 \
 --process MYHLT \
 --data --globaltag 140X_dataRun3_HLT_for2024TSGStudies_v1 \
 --eras Run3 --timing \
 --full --offline --output minimal >hlt_muon_data_Timing.py

hltGetConfiguration /dev/CMSSW_14_0_0/GRun/V107 \
 --process MYHLT \
 --data --globaltag 140X_dataRun3_HLT_for2024TSGStudies_v1 \
 --eras Run3 --timing \
 --customise \
HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeIOSeedingPatatrack_withIter3FromL1_wp00,\
HLTrigger/Configuration/MuonHLTForRun3/customizeMuonHLTForRun3.customizeIOSeedingPatatrack_withIter3FromL1_wp04 \
 --full --offline --output minimal >hlt_muon_data_Timing_BDT_wp04.py


### MC (Efficiency)
hltGetConfiguration /dev/CMSSW_14_0_0/GRun \
 --process MYHLT \
 --mc --globaltag auto:phase1_2024_realistic \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_CascadeMu100_v*,\
HLT_HighPtTkMu100_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/mc/Run3Winter24Reco/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/AODSIM/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50000/1066aeda-aac7-43ef-9599-f531496d32fd.root \
 --eras Run3 \
 --max-events 1000 \
 --full --offline --no-output >hlt_muon_mc.py

# Full L1 emulation : --eras Run3 --l1-emulator FullMC --l1 L1Menu_Collisions2024_v1_1_0_xml \

# Winter23 :  --input /store/mc/Run3Winter23Reco/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/AODSIM/KeepSi_RnD_126X_mcRun3_2023_forPU65_v1-v2/2540000/50afc1f9-64ad-4ae3-8cc3-a43f28776b94.root \
# Summer23BPix :  --input /store/mc/Run3Summer23BPixDRPremix/DYto2L_M-50_TuneCP5_13p6TeV_pythia8/AODSIM/KeepSi_130X_mcRun3_2023_realistic_postBPix_v2-v3/2560000/106765c3-6838-4ef8-bd9c-7487e3ab034e.root \

### 2018Data (Efficiency)
hltGetConfiguration /dev/CMSSW_13_0_0/GRun \
 --process MYHLT \
 --data --globaltag auto:run3_hlt \
 --unprescale \
 --paths \
HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_Mu50_v*,\
HLT_CascadeMu100_v*,\
HLT_HighPtTkMu100_v*,\
HLT_Mu15_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
 --input /store/data/Run2018D/SingleMuon/RAW-RECO/ZMu-12Nov2019_UL2018-v6/40014/F7DED4A7-8B3A-574C-8805-2F4061C87ADA.root \
 --customise HLTrigger/Configuration/customizeHLTforCMSSW.customiseFor2018Input \
 --eras Run3 \
 --max-events 100 \
 --full --offline --no-output >hlt_muon_data_Run2018.py


## Test run
cmsRun hlt_muon_data.py
cmsRun hlt_muon_data_BDT_wp04.py
cmsRun hlt_muon_mc.py
```

## Ntupler
```
git clone -b Run3 https://github.com/wonpoint4/MuonHLTTool.git
scram b -j8

## Test run (Data)
cat <<@EOF >> hlt_muon_data.py
isDIGI = False         # set True (False) for GEN-SIM-DIGI-RAW (GEN-SIM-RAW or Data)

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
process = customizerFuncForMuonHLTNtupler(process, "MYHLT", isDIGI)

#from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
#process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", isDIGI)

# -- L2 seed stat recovery -- #
#process.hltIterL3MuonPixelTracksTrackingRegions.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.hltL3MuonsIterL3IO.L3TrajBuilderParameters.MuonTrackingRegionBuilder.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.HLTIterL3OIAndIOFromL2muonTkCandidateSequence = cms.Sequence(
#    process.HLTIterL3OImuonTkCandidateSequence +
#    process.hltIterL3OIL3MuonsLinksCombination +
#    process.hltIterL3OIL3Muons +
#    process.hltIterL3OIL3MuonCandidates +
#    #process.hltL2SelectorForL3IO +
#    process.HLTIterL3IOmuonTkCandidateSequence +
#    process.hltIterL3MuonsFromL2LinksCombination
#)

process.schedule = cms.Schedule(
     process.HLTriggerFirstPath,
     process.HLT_IsoMu24_v23,
     process.HLT_Mu50_v23,
     process.HLT_CascadeMu100_v11,
     process.HLT_HighPtTkMu100_v10,
     process.HLT_Mu15_v13,
     process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v15,
     process.HLTriggerFinalPath,
     process.mypath,
     process.myendpath,
     #process.myseedpath
)
@EOF
sed -i 's/numberOfThreads = 4/numberOfThreads = 1/g' hlt_muon_data.py
cmsRun hlt_muon_data.py

## Test run (MC)
cat <<@EOF >> hlt_muon_mc.py
isDIGI = True         # set True (False) for GEN-SIM-DIGI-RAW (GEN-SIM-RAW or Data)

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
process = customizerFuncForMuonHLTNtupler(process, "MYHLT", isDIGI)

from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", isDIGI)

# -- L2 seed stat recovery -- #
#process.hltIterL3MuonPixelTracksTrackingRegions.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.hltL3MuonsIterL3IO.L3TrajBuilderParameters.MuonTrackingRegionBuilder.input = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
#process.HLTIterL3OIAndIOFromL2muonTkCandidateSequence = cms.Sequence(
#    process.HLTIterL3OImuonTkCandidateSequence +
#    process.hltIterL3OIL3MuonsLinksCombination +
#    process.hltIterL3OIL3Muons +
#    process.hltIterL3OIL3MuonCandidates +
#    #process.hltL2SelectorForL3IO +
#    process.HLTIterL3IOmuonTkCandidateSequence +
#    process.hltIterL3MuonsFromL2LinksCombination
#)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Run3Winter24Reco/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/AODSIM/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50000/1066aeda-aac7-43ef-9599-f531496d32fd.root ',
    ),
    secondaryFileNames=cms.untracked.vstring(
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/0a3a1e50-ae4e-4105-a88a-04108871b5b3.root',
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/3e745ce9-88f1-4ac0-b945-d714acccfc4b.root',
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/4e877d16-199b-440e-aede-c75e7228ef22.root',
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/8723351d-8657-4cc7-8e47-38645f3f51f8.root',
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/b05a5d60-89ec-497d-8620-ac7df9bcc66c.root',
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/c91c215f-78d5-4e63-b425-bb0699b83186.root',
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/cc63de5e-d7f5-4976-b2be-e64e03977eed.root',
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/d9944467-85bd-4573-a001-72a9946b3199.root',
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/ea271ea5-2f8f-4236-8d4d-81dcfebb5654.root',
        '/store/mc/Run3Winter24Digi/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/KeepSi_133X_mcRun3_2024_realistic_v8-v2/50009/ef0f4597-9c8d-47ca-a5dd-8de787b3d975.root',
    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.schedule = cms.Schedule(
     process.HLTriggerFirstPath,
     process.HLT_IsoMu24_v23,
     process.HLT_Mu50_v23,
     process.HLT_CascadeMu100_v11,
     process.HLT_HighPtTkMu100_v10,
     process.HLT_Mu15_v13,
     process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v15,
     process.HLTriggerFinalPath,
     process.mypath,
     process.myendpath,
     #process.myseedpath
)
@EOF
sed -i 's/numberOfThreads = 4/numberOfThreads = 1/g' hlt_muon_mc.py
cmsRun hlt_muon_mc.py
```
