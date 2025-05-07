# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTSeedNtupler import *
# process = customizerFuncForMuonHLTSeedNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms
# from RecoMuon.TrackerSeedGenerator.mvaScale import *

# -- Std Transform -- #
PU200_Barrel_NThltIter2FromL1_ScaleMean     = [0.00033113700731766336, 1.6825601468762878e-06, 1.790932122524803e-06, 0.010534608406382916, 0.005969459957330139, 0.0009605022254971113, 0.04384189672781466, 7.846741237608237e-05, 0.40725050850004824, 0.41125151617410227, 0.39815551065544846]
PU200_Barrel_NThltIter2FromL1_ScaleStd      = [0.0006042948363798624, 2.445644111872427e-06, 3.454992543447134e-06, 0.09401581628887255, 0.7978806947573766, 0.4932933044535928, 0.04180518265631776, 0.058296511682094855, 0.4071857009373577, 0.41337782307392973, 0.4101160349549534]
PU200_Endcap_NThltIter2FromL1_ScaleMean     = [0.00022658482374555603, 5.358921973784045e-07, 1.010003713549798e-06, 0.0007886873612224615, 0.001197730548842408, -0.0030252353426003594, 0.07151944804171254, -0.0006940626775109026, 0.20535152195939896, 0.2966816533783824, 0.28798220230180455]
PU200_Endcap_NThltIter2FromL1_ScaleStd      = [0.0003857726789049956, 1.4853721474087994e-06, 6.982997036736564e-06, 0.04071340757666084, 0.5897606560095399, 0.33052121398064654, 0.05589386786541949, 0.08806273533388546, 0.3254586902665612, 0.3293354496231377, 0.3179899794578072]


def customizerFuncForMuonHLTSeedNtupler(process, newProcessName = "MYHLT", doDYSkim = False):
    if hasattr(process, "DQMOutput"):
        del process.DQMOutput

    from MuonHLTTool.MuonHLTNtupler.ntupler_seed_cfi import seedNtuplerBase
    import SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi
    from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import tpClusterProducer as _tpClusterProducer

    process.hltTPClusterProducer = _tpClusterProducer.clone(
        phase2OTClusterSrc = cms.InputTag("hltSiPhase2Clusters"),
        pixelClusterSrc = cms.InputTag("hltSiPixelClusters"),
    )
    process.hltTPClusterProducer.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel")
    process.hltTrackAssociatorByHits = SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi.quickTrackAssociatorByHits.clone()
    process.hltTrackAssociatorByHits.cluster2TPSrc            = cms.InputTag("hltTPClusterProducer")
    process.hltTrackAssociatorByHits.UseGrouped               = cms.bool( False )
    process.hltTrackAssociatorByHits.UseSplitting             = cms.bool( False )
    process.hltTrackAssociatorByHits.ThreeHitTracksAreSpecial = cms.bool( False )

    process.seedNtupler = seedNtuplerBase.clone()

    # process.seedNtupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", newProcessName)
    # process.seedNtupler.L1Muon           = cms.untracked.InputTag("gmtStage2Digis",        "Muon", newProcessName)
    # process.seedNtupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", "HLT") #for phaseII w/o emulation
    process.seedNtupler.L1Muon           = cms.untracked.InputTag("simGmtStage2Digis","",newProcessName)  # Phase II sim emul
    process.seedNtupler.L2Muon           = cms.untracked.InputTag("hltL2MuonFromL1TkMuonCandidates", "", "MYHLT")

    process.seedNtupler.L1TkMuon                                          = cms.untracked.InputTag("l1tTkMuonsGmt")
    process.seedNtupler.L1PrimaryVertex                                   = cms.untracked.InputTag("l1tVertexFinderEmulator",            "l1verticesEmulation", newProcessName)

    process.seedNtupler.hltIterL3OISeedsFromL2Muons                       = cms.untracked.InputTag("hltPhase2L3OISeedsFromL2Muons",                         "", newProcessName)
    process.seedNtupler.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.untracked.InputTag("hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks", "", newProcessName)
    process.seedNtupler.hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeeds",                "", newProcessName)

    process.seedNtupler.hltIterL3OIMuonTrack                              = cms.untracked.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity",             "", newProcessName)
    process.seedNtupler.hltIter0IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity",  "", newProcessName)
    process.seedNtupler.hltIter2IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity",  "", newProcessName)

    # process.seedNtupler.associatePixel = cms.bool(True)
    # process.seedNtupler.associateRecoTracks = cms.bool(False)
    # process.seedNtupler.associateStrip = cms.bool(True)
    # process.seedNtupler.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel")
    # process.seedNtupler.stripSimLinkSrc = cms.InputTag("simSiStripDigis")
    # process.seedNtupler.ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof', 'g4SimHitsTrackerHitsPixelBarrelHighTof', 'g4SimHitsTrackerHitsPixelEndcapLowTof', 'g4SimHitsTrackerHitsPixelEndcapHighTof')
    # process.seedNtupler.usePhase2Tracker = cms.bool(True)
    # process.seedNtupler.phase2TrackerSimLinkSrc = cms.InputTag("simSiPixelDigis","Tracker")

    process.seedNtupler.associator = cms.untracked.InputTag("hltTrackAssociatorByHits")
    process.seedNtupler.trackingParticle = cms.untracked.InputTag("mix","MergedTrackTruth")

    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_B_0                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIterL3OI_0.xml")
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_B_1                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIterL3OI_1.xml")
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_B_2                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIterL3OI_2.xml")
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_B_3                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIterL3OI_3.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_0       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter0_0.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_1       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter0_1.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_2       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter0_2.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_B_3       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter0_3.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_0                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter2_0.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_1                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter2_1.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_2                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter2_2.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_B_3                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter2_3.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_0                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter3_0.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_1                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter3_1.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_2                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter3_2.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_B_3                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter3_3.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter0FromL1_0.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter0FromL1_1.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter0FromL1_2.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter0FromL1_3.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Barrel_Binary_NThltIter2FromL1_0.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/xgb_CMSSW_13_0_9_allData_updated_Barrel_NThltIter2FromL1_0.xml")
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Phase2_Iter2FromL1_barrel_v0.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_1                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Barrel_Binary_NThltIter2FromL1_1.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_2                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Barrel_Binary_NThltIter2FromL1_2.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_3                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Barrel_Binary_NThltIter2FromL1_3.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter3FromL1_0.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_1                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter3FromL1_1.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_2                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter3FromL1_2.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_B_3                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Barrel_hltIter3FromL1_3.xml")
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_E_0                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIterL3OI_0.xml")
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_E_1                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIterL3OI_1.xml")
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_E_2                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIterL3OI_2.xml")
    # process.seedNtupler.mvaFileHltIterL3OISeedsFromL2Muons_E_3                       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIterL3OI_3.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_0       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter0_0.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_1       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter0_1.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_2       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter0_2.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3MuonPixelSeedsFromPixelTracks_E_3       = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter0_3.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_0                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter2_0.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_1                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter2_1.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_2                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter2_2.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3MuonPixelSeeds_E_3                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter2_3.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_0                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter3_0.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_1                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter3_1.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_2                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter3_2.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3MuonPixelSeeds_E_3                      = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter3_3.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_0 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter0FromL1_0.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_1 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter0FromL1_1.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_2 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter0FromL1_2.xml")
    # process.seedNtupler.mvaFileHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E_3 = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter0FromL1_3.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Endcap_Binary_NThltIter2FromL1_0.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/xgb_CMSSW_13_0_9_allData_updated_Endcap_NThltIter2FromL1_0.xml")
    process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0                = cms.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Phase2_Iter2FromL1_endcap_v0.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_1                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Endcap_Binary_NThltIter2FromL1_1.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_2                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Endcap_Binary_NThltIter2FromL1_2.xml")
    # process.seedNtupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_3                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/DY_PU200_Endcap_Binary_NThltIter2FromL1_3.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_0                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter3FromL1_0.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_1                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter3FromL1_1.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_2                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter3FromL1_2.xml")
    # process.seedNtupler.mvaFileHltIter3IterL3FromL1MuonPixelSeeds_E_3                = cms.FileInPath("HLTrigger/MuonHLTSeedMVAClassifierPhase2/data/PU180to200Endcap_hltIter3FromL1_3.xml")


    # process.seedNtupler.mvaScaleMeanHltIterL3OISeedsFromL2Muons_B                       = cms.vdouble(PU180to200Barrel_hltIterL3OI_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIterL3OISeedsFromL2Muons_B                        = cms.vdouble(PU180to200Barrel_hltIterL3OI_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_B       = cms.vdouble(PU180to200Barrel_hltIter0_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_B        = cms.vdouble(PU180to200Barrel_hltIter0_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_B                      = cms.vdouble(PU180to200Barrel_hltIter2_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_B                       = cms.vdouble(PU180to200Barrel_hltIter2_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter3IterL3MuonPixelSeeds_B                      = cms.vdouble(PU180to200Barrel_hltIter3_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter3IterL3MuonPixelSeeds_B                       = cms.vdouble(PU180to200Barrel_hltIter3_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B = cms.vdouble(PU180to200Barrel_hltIter0FromL1_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_B  = cms.vdouble(PU180to200Barrel_hltIter0FromL1_ScaleStd)
    process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleMean)
    process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_B                = cms.vdouble(PU180to200Barrel_hltIter3FromL1_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_B                 = cms.vdouble(PU180to200Barrel_hltIter3FromL1_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIterL3OISeedsFromL2Muons_E                       = cms.vdouble(PU180to200Endcap_hltIterL3OI_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIterL3OISeedsFromL2Muons_E                        = cms.vdouble(PU180to200Endcap_hltIterL3OI_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter0IterL3MuonPixelSeedsFromPixelTracks_E       = cms.vdouble(PU180to200Endcap_hltIter0_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter0IterL3MuonPixelSeedsFromPixelTracks_E        = cms.vdouble(PU180to200Endcap_hltIter0_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter2IterL3MuonPixelSeeds_E                      = cms.vdouble(PU180to200Endcap_hltIter2_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter2IterL3MuonPixelSeeds_E                       = cms.vdouble(PU180to200Endcap_hltIter2_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter3IterL3MuonPixelSeeds_E                      = cms.vdouble(PU180to200Endcap_hltIter3_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter3IterL3MuonPixelSeeds_E                       = cms.vdouble(PU180to200Endcap_hltIter3_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E = cms.vdouble(PU180to200Endcap_hltIter0FromL1_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_E  = cms.vdouble(PU180to200Endcap_hltIter0FromL1_ScaleStd)
    process.seedNtupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleMean)
    process.seedNtupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleStd)
    # process.seedNtupler.mvaScaleMeanHltIter3IterL3FromL1MuonPixelSeeds_E                = cms.vdouble(PU180to200Endcap_hltIter3FromL1_ScaleMean)
    # process.seedNtupler.mvaScaleStdHltIter3IterL3FromL1MuonPixelSeeds_E                 = cms.vdouble(PU180to200Endcap_hltIter3FromL1_ScaleStd)


    process.TFileService = cms.Service("TFileService",
      fileName = cms.string("test_ntuple_seedntuple.root"),
      closeFileFast = cms.untracked.bool(False),
    )

    # from MuonHLTTool.MuonHLTNtupler.WmuSkimmer import WmuSkimmer 
    # process.WmuSkimmer = WmuSkimmer.clone()

    process.myseedpath = cms.Path(process.hltTPClusterProducer*process.hltTrackAssociatorByHits*process.seedNtupler)

    return process
