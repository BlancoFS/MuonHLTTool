# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
# process = customizerFuncForMuonHLTNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms

def cutomizerFuncForHltGeneralTracks(process, newProcessName = "MYHLT"):
  
  process.hltInitialStepTracks.TrajectoryInEvent = True
  process.hltInitialStepTrackSelectionHighPurity.copyTrajectories = True
  process.hltGeneralTracks = TrackCollectionMerger.clone(
	      trackProducers   = ["hltInitialStepTrackSelectionHighPurity"],
	      inputClassifiers = ["hltInitialStepTrackCutClassifier"],
        foundHitBonus  = 5.0,
	      lostHitPenalty = 5.0,
        minQuality = cms.string('highPurity'),
        copyExtras = cms.untracked.bool(True),
        copyTrajectories = cms.untracked.bool(True),
  )
  return process

def customizerFuncForMuonHLTNtupler(process, newProcessName = "MYHLT", doDYSkim = False):
    
  from MuonHLTTool.MuonHLTNtupler.ntupler_cfi import ntuplerBase
  process.ntupler = ntuplerBase.clone()

  process.ntupler.associator = cms.untracked.InputTag("hltTrackAssociatorByHits")

  process.TFileService = cms.Service("TFileService",
    fileName = cms.string("seedNtuple_D110Geo_DYToLL.root"),
    closeFileFast = cms.untracked.bool(False),
  )

  #process.mypath    = cms.Path(process.hltTPClusterProducer*process.hltTrackAssociatorByHits)
  process.myendpath = cms.EndPath(process.ntupler)

  return process