import FWCore.ParameterSet.Config as cms
from MuonHLTTool.MuonHLTNtupler.ntupler_cfi import ntuplerBase

def _retagged(tag, newProcessName):
    """Re-point an InputTag at newProcessName, but only if it was
    originally tagged "MYHLT" -- tags with no process name (offline
    collections like "muons", "addPileupInfo", ...) or a different,
    intentional process name are left untouched."""
    if tag.getProcessName() != "MYHLT":
        return tag
    return cms.InputTag(tag.getModuleLabel(), tag.getProductInstanceLabel(), newProcessName)

def customizeToAddNtuplizer(process, newProcessName="MYHLT"):

    if hasattr(process, "DQMOutput"):
        del process.DQMOutput
    
    doNtuple = True
    if doNtuple:
        process.TFileService = cms.Service(
            "TFileService",
            fileName=cms.string("muonHLT_ntuple.root"),
            closeFileFast=cms.untracked.bool(False),
        )        
        process.ntupler = ntuplerBase.clone()

        if newProcessName != "MYHLT":
            for field in (
                    "myTriggerResults", "myTriggerEvent", "lumiScaler",
                    "hltIterL3MuonTrimmedPixelVertices", "hltIterL3FromL1MuonTrimmedPixelVertices",
            ):
                setattr(process.ntupler, field, _retagged(getattr(process.ntupler, field), newProcessName))

            for pset in process.ntupler.muonCollections:
                pset.label = _retagged(pset.label, newProcessName)
        
        process.myendpath = cms.EndPath(process.ntupler)

    doDQMOut = False
    if doDQMOut:
        process.dqmOutput = cms.OutputModule(
            "DQMRootOutputModule",
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
        process.writeDataset = cms.OutputModule(
            "PoolOutputModule",
            fileName = cms.untracked.string('edmOutput.root'),
            outputCommands = cms.untracked.vstring(
                'drop *',
                'keep *_*_*_MYHLT'
            )
        )
        process.EDMOutput = cms.EndPath(process.writeDataset)
    
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
        process.myendpath,
    )
    
