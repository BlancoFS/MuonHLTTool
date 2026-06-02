# MuonHLT Ntupler

## HLT Phase2
``` 
cmsrel CMSSW_16_1_0_pre4
cd CMSSW_16_1_0_pre4/src
cmsenv
git cms-init

git cms-addpkg HLTrigger/Configuration

cmsDriver.py Phase2 -s L1,L1TrackTrigger,L1P2GT,HLT:75e33,VALIDATION:@hltValidation \
	     --processName=MYHLT --conditions auto:phase2_realistic_T33 --geometry ExtendedRun4D110 \
	     --era Phase2C17I13M9 --eventcontent FEVTDEBUGHLT,DQMIO \
	     --datatier GEN-SIM-DIGI-RAW-MINIAOD,DQMIO \
	     --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2FEVTDEBUGHLT.customisePhase2FEVTDEBUGHLT,L1Trigger/Configuration/customisePhase2TTOn110.customisePhase2TTOn110	\
	     --filein /store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/67e21bae-f9cd-43f1-8974-e163400220f7.root \
	     --fileout file:output_Phase2_L1T.root \
	     --python_filename hlt_muon_mc_default.py \
	     '--inputCommands=keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' \
	     --mc -n 100 --nThreads 1 --no_exec		
```

## Phase2 MuonHLTNtupler
``` 
git clone git@github.com:BlancoFS/MuonHLTTool.git -b GNN
scram b -j 5
```
## Add lines to configuration
```
cat <<@EOF >> hlt_muon_mc_default.py
       process.source.inputCommands = cms.untracked.vstring(
    'keep *',
    'drop l1tPFJets_*_*_*',
    'drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT',
    'drop *_hlt*_*_HLT',
    'drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT',
    'drop l1tPFCandidates_*_*_RECO'
)

# -- Ntuple, DQMOutput, and EDMOutput -- #
doNtuple = True
if doNtuple:
    from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
    process = cutomizerFuncForHltGeneralTracks(process, "MYHLT")
    process = customizerFuncForMuonHLTNtupler(process, "MYHLT", False)
    process.TFileService.fileName = cms.string("seedNtuple_D110Geo_DYToLL.root")


doDQMOut = False
if doDQMOut:
    process.dqmOutput = cms.OutputModule("DQMRootOutputModule",
        dataset = cms.untracked.PSet(
            dataTier = cms.untracked.string('DQMIO'),
            filterName = cms.untracked.string('')
        ),
        fileName = cms.untracked.string("DQMIO.root"),
        outputCommands = process.DQMEventContent.outputCommands,
        splitLevel = cms.untracked.int32(0)
    )
    process.DQMOutput = cms.EndPath( process.dqmOutput )

doEDMOut = False
if doEDMOut:
    process.writeDataset = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('edmOutput.root'),
        outputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_*_*_MYHLT'
        )
    )
    process.EDMOutput = cms.EndPath(process.writeDataset)
# -- #
process.schedule = cms.Schedule(
    process.L1simulation_step,
    process.L1TrackTrigger_step,
    process.Phase2L1GTProducer,
    process.Phase2L1GTAlgoBlockProducer,
    process.pTripleTkMuon_5_3_0_DoubleTkMuon_5_3_OS_MassTo9,
    process.pTripleTkMuon_5_3p5_2p5_OS_Mass5to17,
    process.pDoubleEGEle37_24,
    process.pDoubleIsoTkPho22_12,
    process.pDoublePuppiJet112_112,
    process.pDoublePuppiJet160_35_mass620,
    process.pDoublePuppiTau52_52,
    process.pDoubleTkEle25_12,
    process.pDoubleTkElePuppiHT_8_8_390,
    process.pDoubleTkMuPuppiHT_3_3_300,
    process.pDoubleTkMuPuppiJetPuppiMet_3_3_60_130,
    process.pDoubleTkMuon15_7,
    process.pDoubleTkMuonTkEle5_5_9,
    process.pDoubleTkMuon_4_4_OS_Dr1p2,
    process.pDoubleTkMuon_4p5_4p5_OS_Er2_Mass7to18,
    process.pDoubleTkMuon_OS_Er1p5_Dr1p4,
    process.pIsoTkEleEGEle22_12,
    process.pNNPuppiTauPuppiMet_55_190,
    process.pPuppiHT400,
    process.pPuppiHT450,
    process.pPuppiMET200,
    process.pPuppiMHT140,
    process.pPuppiTauTkIsoEle45_22,
    process.pPuppiTauTkMuon42_18,
    process.pQuadJet70_55_40_40,
    process.pSingleEGEle51,
    process.pSingleIsoTkEle28,
    process.pSingleIsoTkPho36,
    process.pSinglePuppiJet230,
    process.pSingleTkEle36,
    process.pSingleTkMuon22,
    process.pTkEleIsoPuppiHT_26_190,
    process.pTkElePuppiJet_28_40_MinDR,
    process.pTkEleTkMuon10_20,
    process.pTkMuPuppiJetPuppiMet_3_110_120,
    process.pTkMuTriPuppiJet_12_40_dRMax_DoubleJet_dEtaMax,
    process.pTkMuonDoubleTkEle6_17_17,
    process.pTkMuonPuppiHT6_320,
    process.pTkMuonTkEle7_23,
    process.pTkMuonTkIsoEle7_20,
    process.pTripleTkMuon5_3_3,
    process.HLT_Mu50_FromL1TkMuon,
    process.HLT_IsoMu24_FromL1TkMuon,
    process.HLT_Mu37_Mu27_FromL1TkMuon,
    process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_FromL1TkMuon,
    process.HLT_TriMu_10_5_5_DZ_FromL1TkMuon,
    process.HLTriggerFinalPath,
    process.mypath,
    process.myendpath,
)
```