# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Phase2 -s L1,L1TrackTrigger,L1P2GT,HLT:75e33 --processName=MYHLT --conditions auto:phase2_realistic_T33 --geometry ExtendedRun4D110 --era Phase2C17I13M9 --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW-MINIAOD --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2FEVTDEBUGHLT.customisePhase2FEVTDEBUGHLT,L1Trigger/Configuration/customisePhase2TTOn110.customisePhase2TTOn110 --filein /store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/67e21bae-f9cd-43f1-8974-e163400220f7.root --fileout file:output_Phase2_L1T.root --python_filename hlt_muon_mc.py --inputCommands=keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT --mc -n 100 --nThreads 1 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('MYHLT',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.SimPhase2L1GlobalTriggerEmulator_cff')
process.load('L1Trigger.Configuration.Phase2GTMenus.SeedDefinitions.step1_2024.l1tGTMenu_cff')
process.load('HLTrigger.Configuration.HLT_75e33_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('/store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/67e21bae-f9cd-43f1-8974-e163400220f7.root'),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop l1tPFJets_*_*_*',
        'drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT',
        'drop *_hlt*_*_HLT',
        'drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Phase2 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-MINIAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:output_Phase2_L1T.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.Phase2L1GTProducer = cms.Path(process.l1tGTProducerSequence)
process.Phase2L1GTAlgoBlockProducer = cms.Path(process.l1tGTAlgoBlockProducerSequence)
process.pDoubleEGEle37_24 = cms.Path(process.DoubleEGEle3724)
process.pDoubleIsoTkPho22_12 = cms.Path(process.DoubleIsoTkPho2212)
process.pDoublePuppiJet112_112 = cms.Path(process.DoublePuppiJet112112)
process.pDoublePuppiJet160_35_mass620 = cms.Path(process.DoublePuppiJet16035Mass620)
process.pDoublePuppiTau52_52 = cms.Path(process.DoublePuppiTau5252)
process.pDoubleTkEle25_12 = cms.Path(process.DoubleTkEle2512)
process.pDoubleTkElePuppiHT_8_8_390 = cms.Path(process.DoubleTkElePuppiHT)
process.pDoubleTkMuPuppiHT_3_3_300 = cms.Path(process.DoubleTkMuPuppiHT)
process.pDoubleTkMuPuppiJetPuppiMet_3_3_60_130 = cms.Path(process.DoubleTkMuPuppiJetPuppiMet)
process.pDoubleTkMuon15_7 = cms.Path(process.DoubleTkMuon157)
process.pDoubleTkMuonTkEle5_5_9 = cms.Path(process.DoubleTkMuonTkEle559)
process.pDoubleTkMuon_4_4_OS_Dr1p2 = cms.Path(process.DoubleTkMuon44OSDr1p2)
process.pDoubleTkMuon_4p5_4p5_OS_Er2_Mass7to18 = cms.Path(process.DoubleTkMuon4p5OSEr2Mass7to18)
process.pDoubleTkMuon_OS_Er1p5_Dr1p4 = cms.Path(process.DoubleTkMuonOSEr1p5Dr1p4)
process.pIsoTkEleEGEle22_12 = cms.Path(process.IsoTkEleEGEle2212)
process.pNNPuppiTauPuppiMet_55_190 = cms.Path(process.NNPuppiTauPuppiMet)
process.pPuppiHT400 = cms.Path(process.PuppiHT400)
process.pPuppiHT450 = cms.Path(process.PuppiHT450)
process.pPuppiMET200 = cms.Path(process.PuppiMET200)
process.pPuppiMHT140 = cms.Path(process.PuppiMHT140)
process.pPuppiTauTkIsoEle45_22 = cms.Path(process.PuppiTauTkIsoEle4522)
process.pPuppiTauTkMuon42_18 = cms.Path(process.PuppiTauTkMuon4218)
process.pQuadJet70_55_40_40 = cms.Path(process.QuadJet70554040)
process.pSingleEGEle51 = cms.Path(process.SingleEGEle51)
process.pSingleIsoTkEle28 = cms.Path(process.SingleIsoTkEle28)
process.pSingleIsoTkPho36 = cms.Path(process.SingleIsoTkPho36)
process.pSinglePuppiJet230 = cms.Path(process.SinglePuppiJet230)
process.pSingleTkEle36 = cms.Path(process.SingleTkEle36)
process.pSingleTkMuon22 = cms.Path(process.SingleTkMuon22)
process.pTkEleIsoPuppiHT_26_190 = cms.Path(process.TkEleIsoPuppiHT)
process.pTkElePuppiJet_28_40_MinDR = cms.Path(process.TkElePuppiJetMinDR)
process.pTkEleTkMuon10_20 = cms.Path(process.TkEleTkMuon1020)
process.pTkMuPuppiJetPuppiMet_3_110_120 = cms.Path(process.TkMuPuppiJetPuppiMet)
process.pTkMuTriPuppiJet_12_40_dRMax_DoubleJet_dEtaMax = cms.Path(process.TkMuTriPuppiJetdRMaxDoubleJetdEtaMax)
process.pTkMuonDoubleTkEle6_17_17 = cms.Path(process.TkMuonDoubleTkEle61717)
process.pTkMuonPuppiHT6_320 = cms.Path(process.TkMuonPuppiHT6320)
process.pTkMuonTkEle7_23 = cms.Path(process.TkMuonTkEle723)
process.pTkMuonTkIsoEle7_20 = cms.Path(process.TkMuonTkIsoEle720)
process.pTripleTkMuon5_3_3 = cms.Path(process.TripleTkMuon533)
process.pTripleTkMuon_5_3_0_DoubleTkMuon_5_3_OS_MassTo9 = cms.Path(process.TripleTkMuon530OSMassMax9)
process.pTripleTkMuon_5_3p5_2p5_OS_Mass5to17 = cms.Path(process.TripleTkMuon53p52p5OSMass5to17)
process.L1T_DoubleNNTau52 = cms.Path(process.HLTL1Sequence+process.hltL1DoubleNNTau52)
process.L1T_SingleNNTau150 = cms.Path(process.HLTL1Sequence+process.hltL1SingleNNTau150)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
# process.schedule imported from cff in HLTrigger.Configuration
process.schedule.insert(0, process.L1simulation_step)
process.schedule.insert(1, process.L1TrackTrigger_step)
process.schedule.insert(2, process.Phase2L1GTProducer)
process.schedule.insert(3, process.Phase2L1GTAlgoBlockProducer)
process.schedule.insert(4, process.pDoubleEGEle37_24)
process.schedule.insert(5, process.pDoubleIsoTkPho22_12)
process.schedule.insert(6, process.pDoublePuppiJet112_112)
process.schedule.insert(7, process.pDoublePuppiJet160_35_mass620)
process.schedule.insert(8, process.pDoublePuppiTau52_52)
process.schedule.insert(9, process.pDoubleTkEle25_12)
process.schedule.insert(10, process.pDoubleTkElePuppiHT_8_8_390)
process.schedule.insert(11, process.pDoubleTkMuPuppiHT_3_3_300)
process.schedule.insert(12, process.pDoubleTkMuPuppiJetPuppiMet_3_3_60_130)
process.schedule.insert(13, process.pDoubleTkMuon15_7)
process.schedule.insert(14, process.pDoubleTkMuonTkEle5_5_9)
process.schedule.insert(15, process.pDoubleTkMuon_4_4_OS_Dr1p2)
process.schedule.insert(16, process.pDoubleTkMuon_4p5_4p5_OS_Er2_Mass7to18)
process.schedule.insert(17, process.pDoubleTkMuon_OS_Er1p5_Dr1p4)
process.schedule.insert(18, process.pIsoTkEleEGEle22_12)
process.schedule.insert(19, process.pNNPuppiTauPuppiMet_55_190)
process.schedule.insert(20, process.pPuppiHT400)
process.schedule.insert(21, process.pPuppiHT450)
process.schedule.insert(22, process.pPuppiMET200)
process.schedule.insert(23, process.pPuppiMHT140)
process.schedule.insert(24, process.pPuppiTauTkIsoEle45_22)
process.schedule.insert(25, process.pPuppiTauTkMuon42_18)
process.schedule.insert(26, process.pQuadJet70_55_40_40)
process.schedule.insert(27, process.pSingleEGEle51)
process.schedule.insert(28, process.pSingleIsoTkEle28)
process.schedule.insert(29, process.pSingleIsoTkPho36)
process.schedule.insert(30, process.pSinglePuppiJet230)
process.schedule.insert(31, process.pSingleTkEle36)
process.schedule.insert(32, process.pSingleTkMuon22)
process.schedule.insert(33, process.pTkEleIsoPuppiHT_26_190)
process.schedule.insert(34, process.pTkElePuppiJet_28_40_MinDR)
process.schedule.insert(35, process.pTkEleTkMuon10_20)
process.schedule.insert(36, process.pTkMuPuppiJetPuppiMet_3_110_120)
process.schedule.insert(37, process.pTkMuTriPuppiJet_12_40_dRMax_DoubleJet_dEtaMax)
process.schedule.insert(38, process.pTkMuonDoubleTkEle6_17_17)
process.schedule.insert(39, process.pTkMuonPuppiHT6_320)
process.schedule.insert(40, process.pTkMuonTkEle7_23)
process.schedule.insert(41, process.pTkMuonTkIsoEle7_20)
process.schedule.insert(42, process.pTripleTkMuon5_3_3)
process.schedule.insert(43, process.pTripleTkMuon_5_3_0_DoubleTkMuon_5_3_OS_MassTo9)
process.schedule.insert(44, process.pTripleTkMuon_5_3p5_2p5_OS_Mass5to17)
process.schedule.extend([process.endjob_step,process.FEVTDEBUGHLToutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000 

#call to customisation function customise_aging_1000 imported from SLHCUpgradeSimulations.Configuration.aging
process = customise_aging_1000(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customisePhase2FEVTDEBUGHLT
from L1Trigger.Configuration.customisePhase2FEVTDEBUGHLT import customisePhase2FEVTDEBUGHLT 

#call to customisation function customisePhase2FEVTDEBUGHLT imported from L1Trigger.Configuration.customisePhase2FEVTDEBUGHLT
process = customisePhase2FEVTDEBUGHLT(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customisePhase2TTOn110
from L1Trigger.Configuration.customisePhase2TTOn110 import customisePhase2TTOn110 

#call to customisation function customisePhase2TTOn110 imported from L1Trigger.Configuration.customisePhase2TTOn110
process = customisePhase2TTOn110(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

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

    process.L1T_SingleTkMuon_22,
    process.L1T_DoubleTkMuon_15_7,
    process.L1T_TripleTkMuon_5_3_3,
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
