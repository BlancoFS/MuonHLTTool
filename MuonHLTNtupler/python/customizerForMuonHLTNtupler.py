# -- customizer for the ntupler, added to the HLT configuration for
#    re-running HLT:
#      from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
#      process = customizerFuncForMuonHLTNtupler(process, "MYHLT")
#
# -----------------------------------------------------------------------
# This used to ALSO be where every process.ntupler.<Field> InputTag got
# set -- ~50 lines that were unconditionally required just to get a
# runnable module, even for the default Phase-2 setup. All of that now
# lives in ntupler_cfi.py as the default configuration. What's left here
# is what's actually customization:
#   1) building the sim-truth track-association producer sequence (a
#      genuine producer-sequence concern, not an ntupler parameter),
#   2) re-tagging process.ntupler's InputTags if you rerun HLT under a
#      process name other than "MYHLT",
#   3) anything else you want to override for a specific study (a
#      different seed-MVA training, a DY skim, non-default DebugMode/
#      SaveAllTracks, etc.) -- add process.ntupler.<Field> = ... lines
#      for just the fields you're changing, same as before.
# -----------------------------------------------------------------------

import FWCore.ParameterSet.Config as cms


def _retagged(tag, newProcessName):
    """Re-point an InputTag at newProcessName, but only if it was
    originally tagged "MYHLT" -- tags with no process name (offline
    collections like "muons", "addPileupInfo", ...) or a different,
    intentional process name are left untouched."""
    if tag.getProcessName() != "MYHLT":
        return tag
    return cms.InputTag(tag.getModuleLabel(), tag.getProductInstanceLabel(), newProcessName)


def customizerFuncForMuonHLTNtupler(process, newProcessName="MYHLT", doDYSkim=False):
    if hasattr(process, "DQMOutput"):
        del process.DQMOutput

    # -----------------------------------------------------------------
    # 1) sim-truth track-association producer sequence. This has to be
    #    built here (not in the cfi) because it creates producers that
    #    must be added to the process and run in the Path -- an ntupler
    #    _parameter_ can't do that. The InputTag *labels* these producers
    #    are given (trackCollectionLabels/associationLabels) already live
    #    in ntupler_cfi.py and must match the names used below.
    # -----------------------------------------------------------------
    import SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi
    from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import tpClusterProducer as _tpClusterProducer

    process.hltTPClusterProducer = _tpClusterProducer.clone(
        phase2OTClusterSrc=cms.InputTag("hltSiPhase2Clusters"),
        pixelClusterSrc=cms.InputTag("hltSiPixelClusters"),
    )
    process.hltTPClusterProducer.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel")

    process.hltTrackAssociatorByHits = SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone()
    process.hltTrackAssociatorByHits.cluster2TPSrc = cms.InputTag("hltTPClusterProducer")
    process.hltTrackAssociatorByHits.UseGrouped = cms.bool(False)
    process.hltTrackAssociatorByHits.UseSplitting = cms.bool(False)
    process.hltTrackAssociatorByHits.ThreeHitTracksAreSpecial = cms.bool(False)

    import SimMuon.MCTruth.MuonTrackProducer_cfi
    process.hltPhase2L3MuonsNoIDTracks = SimMuon.MCTruth.MuonTrackProducer_cfi.muonTrackProducer.clone()
    process.hltPhase2L3MuonsNoIDTracks.muonsTag = cms.InputTag("hltPhase2L3MuonsNoID")
    process.hltPhase2L3MuonsNoIDTracks.selectionTags = ('All',)
    process.hltPhase2L3MuonsNoIDTracks.trackType = "recomuonTrack"
    process.hltPhase2L3MuonsNoIDTracks.ignoreMissingMuonCollection = True
    process.hltPhase2L3MuonsNoIDTracks.inputCSCSegmentCollection = cms.InputTag("hltCscSegments")
    process.hltPhase2L3MuonsNoIDTracks.inputDTRecSegment4DCollection = cms.InputTag("hltDt4DSegments")

    process.hltPhase2L3MuonsTracks = SimMuon.MCTruth.MuonTrackProducer_cfi.muonTrackProducer.clone()
    process.hltPhase2L3MuonsTracks.muonsTag = cms.InputTag("hltPhase2L3Muons")
    process.hltPhase2L3MuonsTracks.selectionTags = ('All',)
    process.hltPhase2L3MuonsTracks.trackType = "recomuonTrack"
    process.hltPhase2L3MuonsTracks.ignoreMissingMuonCollection = True
    process.hltPhase2L3MuonsTracks.inputCSCSegmentCollection = cms.InputTag("hltCscSegments")
    process.hltPhase2L3MuonsTracks.inputDTRecSegment4DCollection = cms.InputTag("hltDt4DSegments")

    from SimMuon.MCTruth.MuonAssociatorByHits_cfi import muonAssociatorByHits as _muonAssociatorByHits
    hltMuonAssociatorByHits = _muonAssociatorByHits.clone()
    hltMuonAssociatorByHits.PurityCut_track = 0.75
    hltMuonAssociatorByHits.PurityCut_muon = 0.75
    hltMuonAssociatorByHits.DTrechitTag = 'hltDt1DRecHits'
    hltMuonAssociatorByHits.ignoreMissingTrackCollection = True
    hltMuonAssociatorByHits.UseTracker = True
    hltMuonAssociatorByHits.UseMuon = True

    process.AhltPhase2L3OIMuonTrackSelectionHighPurity = hltMuonAssociatorByHits.clone(tracksTag='hltPhase2L3OIMuonTrackSelectionHighPurity')
    process.AhltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity = hltMuonAssociatorByHits.clone(tracksTag='hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity')
    process.AhltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity = hltMuonAssociatorByHits.clone(tracksTag='hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity')
    process.AhltIter2Phase2L3FromL1TkMuonMerged = hltMuonAssociatorByHits.clone(tracksTag='hltIter2Phase2L3FromL1TkMuonMerged')
    process.AhltPhase2L3MuonsNoID = hltMuonAssociatorByHits.clone(tracksTag='hltPhase2L3MuonsNoIDTracks')
    process.AhltPhase2L3Muons = hltMuonAssociatorByHits.clone(tracksTag='hltPhase2L3MuonsTracks')

    process.trackAssoSeq = cms.Sequence(
        process.hltPhase2L3MuonsNoIDTracks +
        process.hltPhase2L3MuonsTracks +
        process.AhltPhase2L3OIMuonTrackSelectionHighPurity +
        process.AhltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity +
        process.AhltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity +
        process.AhltIter2Phase2L3FromL1TkMuonMerged +
        process.AhltPhase2L3MuonsNoID +
        process.AhltPhase2L3Muons
    )

    # -----------------------------------------------------------------
    # 2) start from the (now self-sufficient) default configuration
    # -----------------------------------------------------------------
    from MuonHLTTool.MuonHLTNtupler.ntupler_cfi import ntuplerBase
    process.ntupler = ntuplerBase.clone()

    # -----------------------------------------------------------------
    # 3) only re-tag if the caller actually asked for a different HLT
    #    rerun process name than the "MYHLT" default baked into the cfi.
    #    Everything below is a no-op when newProcessName == "MYHLT".
    # -----------------------------------------------------------------
    if newProcessName != "MYHLT":
        for field in (
            "myTriggerResults", "myTriggerEvent", "lumiScaler", "L1Muon",
            "iterL3IOFromL1", "iterL3MuonNoID", "iterL3Muon",
            "hltIterL3MuonTrimmedPixelVertices", "hltIterL3FromL1MuonTrimmedPixelVertices",
            "hltIterL3OISeedsFromL2Muons",
            "hltIter0IterL3MuonPixelSeedsFromPixelTracks", "hltIter2IterL3MuonPixelSeeds",
            "hltIter3IterL3MuonPixelSeeds", "hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks",
            "hltIter2IterL3FromL1MuonPixelSeeds", "hltIter3IterL3FromL1MuonPixelSeeds",
            "hltIterL3OIMuonTrack", "hltIter0IterL3MuonTrack", "hltIter2IterL3MuonTrack",
            "hltIter3IterL3MuonTrack", "hltIter0IterL3FromL1MuonTrack",
            "hltIter2IterL3FromL1MuonTrack", "hltIter3IterL3FromL1MuonTrack",
        ):
            setattr(process.ntupler, field, _retagged(getattr(process.ntupler, field), newProcessName))

        for pset in process.ntupler.muonCollections:
            pset.label = _retagged(pset.label, newProcessName)

    # -----------------------------------------------------------------
    # 4) anything else you want to override for a specific study goes
    #    here, e.g.:
    #      process.ntupler.triggerPaths   = cms.untracked.vstring("HLT_Mu50_v", ...)
    #      process.ntupler.triggerFilters = cms.untracked.vstring("hlt...Filtered50Q", ...)
    #      process.ntupler.doSeed = cms.bool(True)
    #      process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0 = cms.untracked.FileInPath(...)
    #    Only the fields that differ from ntupler_cfi.py's defaults need
    #    to be touched.
    # -----------------------------------------------------------------

    # if doDYSkim:
    #     from MuonHLTTool.MuonHLTNtupler.DYmuSkimmer import DYmuSkimmer
    #     process.Skimmer = DYmuSkimmer.clone()
    #     process.mypath = cms.Path(process.Skimmer * process.hltTPClusterProducer * process.hltTrackAssociatorByHits * process.trackAssoSeq)
    # else:
    #     process.mypath = cms.Path(process.hltTPClusterProducer * process.hltTrackAssociatorByHits * process.trackAssoSeq)

    process.TFileService = cms.Service(
        "TFileService",
        fileName=cms.string("test_ntuple_seedntuple.root"),
        closeFileFast=cms.untracked.bool(False),
    )

    process.mypath = cms.Path(process.hltTPClusterProducer * process.hltTrackAssociatorByHits * process.trackAssoSeq)
    process.myendpath = cms.EndPath(process.ntupler)

    return process
