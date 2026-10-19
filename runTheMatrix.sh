#!/bin/bash

#### ZPrime to mu+ mu-

cmsDriver.py Configuration/GenProduction/python/TSG-Phase2Spring24GS-ZPrimeMM-fragment.py -s GEN,SIM -n 1 --conditions auto:phase2_realistic_T35 --beamspot DBrealisticHLLHC --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry ExtendedRun4D127 --era Phase2C26I13M9 --relval 9000,100 --fileout file:step1.root  --nThreads 8 > step1_ZMM_14TeV+Run4D127.log  2>&1

# in: /afs/cern.ch/work/s/sblancof/private/MuonPOG/hlt_tutorials/PhaseII/TruthGraph/CMSSW_20_1_X_2026-09-28-1100/src going to execute cd 37634.0_TTbar_14TeV+Run4D127
cmsDriver.py step2 -s DIGI:pdigi_valid,L1TrackTrigger,L1,L1P2GT,DIGI2RAW,HLT:@relvalRun4 --conditions auto:phase2_realistic_T35 --datatier GEN-SIM-DIGI-RAW -n 1 --eventcontent FEVTDEBUGHLT --geometry ExtendedRun4D127 --era Phase2C26I13M9 --filein  file:step1.root  --fileout file:step2.root  --nThreads 8 > step2_TTbar_14TeV+Run4D127.log  2>&1
 

# in: /afs/cern.ch/work/s/sblancof/private/MuonPOG/hlt_tutorials/PhaseII/TruthGraph/CMSSW_20_1_X_2026-09-28-1100/src going to execute cd 37634.0_TTbar_14TeV+Run4D127
cmsDriver.py step3  -s RAW2DIGI,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM --conditions auto:phase2_realistic_T35 --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO -n 1 --eventcontent FEVTDEBUGHLT,MINIAODSIM,DQM --geometry ExtendedRun4D127 --era Phase2C26I13M9 --filein  file:step2.root  --fileout file:step3.root  --nThreads 8 > step3_TTbar_14TeV+Run4D127.log  2>&1
 

# in: /afs/cern.ch/work/s/sblancof/private/MuonPOG/hlt_tutorials/PhaseII/TruthGraph/CMSSW_20_1_X_2026-09-28-1100/src going to execute cd 37634.0_TTbar_14TeV+Run4D127
cmsDriver.py step4  -s HARVESTING:@phase2Validation+@phase2+@miniAODValidation+@miniAODDQM --conditions auto:phase2_realistic_T35 --mc  --geometry ExtendedRun4D127 --scenario pp --filetype DQM --era Phase2C26I13M9 -n 1  --filein file:step3_inDQM.root --fileout file:step4.root  > step4_TTbar_14TeV+Run4D127.log  2>&1
 

# in: /afs/cern.ch/work/s/sblancof/private/MuonPOG/hlt_tutorials/PhaseII/TruthGraph/CMSSW_20_1_X_2026-09-28-1100/src going to execute cd 37634.0_TTbar_14TeV+Run4D127
cmsDriver.py step5  -s ALCA:SiPixelCalSingleMuonLoose+SiPixelCalSingleMuonTight+TkAlMuonIsolated+TkAlMinBias+MuAlOverlaps+EcalESAlign+TkAlZMuMu+TkAlDiMuonAndVertex+HcalCalHBHEMuonProducerFilter+TkAlUpsilonMuMu+TkAlJpsiMuMu+TkAlV0s+TkAlJetHT --conditions auto:phase2_realistic_T35 --datatier ALCARECO -n 1 --eventcontent ALCARECO --geometry ExtendedRun4D127 --era Phase2C26I13M9 --filein file:step3.root --fileout file:step5.root  --nThreads 8 > step5_TTbar_14TeV+Run4D127.log  2>&1



##### For ntuplizer

cmsDriver.py Phase2 -s L1,L1TrackTrigger,L1P2GT,HLT:75e33 --processName=MYHLT --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D127 --era Phase2C26I13M9 --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW-MINIAOD --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2FEVTDEBUGHLT.customisePhase2FEVTDEBUGHLT,L1Trigger/Configuration/customisePhase2TTOn110.customisePhase2TTOn110 --filein step3.root --fileout file:muonHLT_ntuple.root --python_filename hlt_muon_phase2.py '--inputCommands=keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' --mc -n -1 --nThreads 1 --no_exec
