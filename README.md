# MuonHLT Ntupler

## Setup (15_0_0_pre3, check the P2UG instruction for further updates: https://cmshltupgrade.docs.cern.ch/RunningInstructions/)
```
cmsrel 15_0_0_pre3
cd 15_0_0_pre3/src
cmsenv
git cms-init

git clone -b Phase2_2024 https://github.com/kyHwangs/MuonHLTTool.git
scram b -j8
```

## Configuration
### L1 + HLT re-emulations
```
cmsDriver.py Phase2 -s L1,L1TrackTrigger,L1P2GT,HLT:75e33 --processName=MYHLT \
--conditions auto:phase2_realistic_T33 \
--geometry ExtendedRun4D110 \
--era Phase2C17I13M9 \
--eventcontent FEVTDEBUGHLT \
--datatier GEN-SIM-DIGI-RAW-MINIAOD \
--customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2FEVTDEBUGHLT.customisePhase2FEVTDEBUGHLT,L1Trigger/Configuration/customisePhase2TTOn110.customisePhase2TTOn110 \
--filein /store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/67e21bae-f9cd-43f1-8974-e163400220f7.root \
--fileout file:output_Phase2_L1T.root \
--python_filename hlt_muon_mc.py \
--inputCommands="keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT" \
--mc \
-n 100 --nThreads 1 --no_exec
```

Add some line on configuration file for test run:
```
cat <<@EOF >> hlt_muon_mc.py

# for CMSSW_15_0_0_preX
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
    process = customizerFuncForMuonHLTNtupler(process, "MYHLT", False)

    process.ntupler.offlineMuon                   = cms.untracked.InputTag("slimmedMuons")
    process.ntupler.TkMuonToken                   = cms.InputTag("l1tTkMuonsGmt")
    # process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeeds", "", "MYHLT")
    process.ntupler.doMVA                         = cms.bool(True)
    # Isolation study
    # process.ntupler.trkIsoTags                    = cms.untracked.vstring(   trkIsoTags )
    # process.ntupler.trkIsoLabels                  = cms.untracked.VInputTag( trkIsoLabels )
    # process.ntupler.pfIsoTags                     = cms.untracked.vstring(   pfIsoTags )
    # process.ntupler.pfIsoLabels                   = cms.untracked.VInputTag( pfIsoLabels )

    from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
    process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT", True)

    process.seedNtupler.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "", "MYHLT")
    # process.seedNtupler.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks", "RECO")

    process.TFileService.fileName = cms.string("seedNtuple_D110Geo_DYToLL.root")

    # from HLTrigger.MuonHLTSeedMVAClassifierPhase2.customizerForMuonHLTSeeding import *
    # WPNAME = 'noMVAcut_noSeedMax'
    # doSort = False
    # nSeedMax_B = (-1,)
    # nSeedMax_E = (-1,)
    # mvaCuts_B = (0,)
    # mvaCuts_E = (0,)
    # process = customizerFuncForMuonHLTSeeding(process, "MYHLT", WPNAME, doSort, nSeedMax_B, nSeedMax_E, mvaCuts_B, mvaCuts_E )
    # process.hltIter2Phase2L3FromL1TkMuonPixelSeedsFiltered.L1TkMu = cms.InputTag("l1tTkMuonsGmt", "", "MYHLT")

#process.l1tTkMuonsGmt.applyQualityCuts = cms.bool(False)

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

    # process.L1T_SingleTkMuon_22,
    # process.L1T_DoubleTkMuon_15_7,
    # process.L1T_TripleTkMuon_5_3_3,
    process.HLT_Mu50_FromL1TkMuon,
    process.HLT_IsoMu24_FromL1TkMuon,
    process.HLT_Mu37_Mu27_FromL1TkMuon,
    process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_FromL1TkMuon,
    process.HLT_TriMu_10_5_5_DZ_FromL1TkMuon,
    # process.HLT_IsoStudy,
    process.HLTriggerFinalPath,
    # process.Gen_QCDBCToEFilter,
    # process.Gen_QCDEmEnrichingFilter,
    # process.Gen_QCDEmEnrichingNoBCToEFilter,
    # process.Gen_QCDMuGenFilter,
    # process.Gen_QCDMuNoEmGenFilter,
    # process.Gen_QCDEmNoMuGenFilter,
    # process.myana,
    process.mypath,
    # process.valpath,
    process.myendpath,
    # process.myseedpath
    # process.DQMOutput
    # process.EDMOutput
 )
@EOF

cmsRun hlt_muon_mc.py
```

Before submit the crab job, make sure the # of thread is set to 4.
```
sed -i 's/numberOfThreads = cms.untracked.uint32(1)/numberOfThreads = cms.untracked.uint32(4)/g' hlt_muon_mc.py
```
