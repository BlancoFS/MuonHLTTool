import FWCore.ParameterSet.Config as cms

# -----------------------------------------------------------------------
# This is now the single source of truth for default running. Everything
# below used to be split between this file (a handful of Run-3-style
# defaults) and customizerForMuonHLTNtupler.py (~50 lines of
# process.ntupler.<Field> = ... overrides that were unconditionally
# required just to get a runnable module). customizerForMuonHLTNtupler.py
# is now only needed when actually changing something -- a different HLT
# process name, a different track-association sequence, a different seed
# MVA training -- see the comments there.
#
# Default assumption: rerunning Phase-2 muon HLT under process name
# "MYHLT". If you rerun under a different name, use the customizer (it
# only needs to re-tag the InputTags below, nothing else).
# -----------------------------------------------------------------------

MYHLT = "MYHLT"

# -- Seed-MVA scale vectors for the one seed stage actually evaluated by
#    default (hltIter2*FromL1*PixelSeeds); PU200, current as of the
#    2026 Phase-2 samples. Swap these in the customizer for a different
#    training/PU scenario -- see customizerForMuonHLTNtupler.py.
PU200_Barrel_NThltIter2FromL1_ScaleMean = [0.00033113700731766336, 1.6825601468762878e-06, 1.790932122524803e-06, 0.010534608406382916, 0.005969459957330139, 0.0009605022254971113, 0.04384189672781466, 7.846741237608237e-05, 0.40725050850004824, 0.41125151617410227, 0.39815551065544846]
PU200_Barrel_NThltIter2FromL1_ScaleStd  = [0.0006042948363798624, 2.445644111872427e-06, 3.454992543447134e-06, 0.09401581628887255, 0.7978806947573766, 0.4932933044535928, 0.04180518265631776, 0.058296511682094855, 0.4071857009373577, 0.41337782307392973, 0.4101160349549534]
PU200_Endcap_NThltIter2FromL1_ScaleMean = [0.00022658482374555603, 5.358921973784045e-07, 1.010003713549798e-06, 0.0007886873612224615, 0.001197730548842408, -0.0030252353426003594, 0.07151944804171254, -0.0006940626775109026, 0.20535152195939896, 0.2966816533783824, 0.28798220230180455]
PU200_Endcap_NThltIter2FromL1_ScaleStd  = [0.0003857726789049956, 1.4853721474087994e-06, 6.982997036736564e-06, 0.04071340757666084, 0.5897606560095399, 0.33052121398064654, 0.05589386786541949, 0.08806273533388546, 0.3254586902665612, 0.3293354496231377, 0.3179899794578072]

# -----------------------------------------------------------------------
# Generic muon-like collections -- adding, removing, or renaming one of
# these (because the HLT sequence changed) is the whole point of this
# refactor: it's a python-only change, no C++ rebuild.
#
# 'type' selects the adapter in GenericMuonCollections.h:
#   RecoChargedCandidate -> reco::RecoChargedCandidateCollection
#   TrackView             -> edm::View<reco::Track> (also covers plain
#                            reco::TrackCollection, e.g. L2Muon)
#   RecoMuon              -> std::vector<reco::Muon>
#   MuonTrackLinks        -> std::vector<reco::MuonTrackLinks> (carries
#                            separate inner/outer/global tracks)
#   L1TkMuon              -> l1t::TrackerMuonCollection
#
# NOTE: fill_trackTemplate()/fill_trackTemplateMva() in the .cc still look
# up the "L3Muon" and "L2Muon" entries BY NAME for an isolation-map size
# check and an MVA gate respectively (see REFACTOR_NOTES.md) -- if you
# rename those two specific entries, update those two lookups in the .cc.
# -----------------------------------------------------------------------
muonCollections = cms.untracked.VPSet(
    cms.PSet(name=cms.string("L2Muon"), type=cms.string("TrackView"),
             label=cms.InputTag("hltL2MuonsFromL1TkMuon", "", MYHLT)),

    cms.PSet(name=cms.string("L3Muon"), type=cms.string("RecoChargedCandidate"),
             label=cms.InputTag("hltPhase2L3MuonCandidates", "", MYHLT)),

    cms.PSet(name=cms.string("iterL3OI"), type=cms.string("TrackView"),
             label=cms.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity", "", MYHLT)),

    cms.PSet(name=cms.string("iterL3IO"), type=cms.string("TrackView"),
             label=cms.InputTag("hltPhase2L3MuonFilter:L3IOTracksFiltered", "", MYHLT)),
)

# -----------------------------------------------------------------------
# Trigger paths / filters for the generic dR-matching in
# MuonHLTTriggerMatching.h. Unversioned path names ("HLT_Mu50_v") are
# resolved against the live menu each run.
#
# These are last-mile filters for the single-L1TkMu22-seeded HLT_IsoMu24/
# HLT_Mu50-style paths in the current dev menu -- update this list if the
# menu's filter names change or you add paths to study.
# -----------------------------------------------------------------------
triggerPaths = cms.untracked.vstring(
    "HLT_IsoMu24_FromL1TkMuon",
    "HLT_Mu50_FromL1TkMuon",
)
triggerFilters = cms.untracked.vstring(
    "hltSingleTkMuon22L1TkMuonFilter",
    "hltL3fL1TkSingleMu22L3Filtered24Q",
    "hltL3fL1TkSingleMu22L3Filtered50Q",
    "hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41",
    "hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40",
    "hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70",
    "hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk",
)
maxDR = cms.untracked.double(0.1)

ntuplerBase = cms.EDAnalyzer(
    "MuonHLTNtupler",

    # -- information stored in the original (non-rerun) edm file
    triggerResults    = cms.untracked.InputTag("TriggerResults::MYHLT"),
    triggerEvent      = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"),
    offlineLumiScaler = cms.untracked.InputTag("scalersRawToDigi"),
    offlineVertex     = cms.untracked.InputTag("offlinePrimaryVertices"),
    offlineMuon       = cms.untracked.InputTag("muons"),

    # -- objects from the HLT rerun (process name = "MYHLT" by default)
    myTriggerResults = cms.untracked.InputTag("TriggerResults",       "", MYHLT),
    myTriggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD", "", MYHLT),
    lumiScaler       = cms.untracked.InputTag("hltScalersRawToDigi",  "", MYHLT),

    L1Muon   = cms.untracked.InputTag("simGmtStage2Digis", "", MYHLT),  # -- Phase-2 sim emulation
    L1TkMuon = cms.untracked.InputTag("l1tTkMuonsGmt"),

    muonCollections = muonCollections,
    triggerPaths    = triggerPaths,
    triggerFilters  = triggerFilters,
    maxDR           = maxDR,

    # -- iterL3Muon/iterL3MuonNoID/iterL3IOFromL1 keep dedicated tokens
    #    (not folded into muonCollections) because Fill_IterL3()/
    #    fill_trackTemplateMva() attach extra structure to them: ID flags,
    #    inner-track pt, and the seed<->muon matching maps
    #    (iterL3IDpassed/iterL3NoIDpassed/MuonIterSeedMap) used by the
    #    seed-MVA study. See REFACTOR_NOTES.md for why these were left out
    #    of the generic migration.
    iterL3IOFromL1 = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonMerged", "", MYHLT),
    iterL3MuonNoID = cms.untracked.InputTag("hltPhase2L3MuonsNoID",               "", MYHLT),
    iterL3Muon     = cms.untracked.InputTag("hltPhase2L3Muons",                   "", MYHLT),

    hltIterL3MuonTrimmedPixelVertices       = cms.untracked.InputTag("hltIterL3MuonTrimmedPixelVertices",       "", MYHLT),
    hltIterL3FromL1MuonTrimmedPixelVertices = cms.untracked.InputTag("hltPhase2L3FromL1TkMuonTrimmedPixelVertices", "", MYHLT),

    doMVA  = cms.bool(False),
    doSeed = cms.bool(False),

    # -- seed collections: only the *FromL1* ones are exercised by default
    #    (doMVA/doSeed are off); the others keep placeholder Run-3-style
    #    tags. Since every getByToken() call for these is wrapped in an
    #    `if (...)` guard in the .cc, a nonexistent product is silently
    #    skipped rather than throwing -- so leaving these unresolved for
    #    Phase-2 running is safe as long as doMVA/doSeed stay False.
    hltIterL3OISeedsFromL2Muons                       = cms.untracked.InputTag("hltPhase2L3OISeedsFromL2Muons",                       "", MYHLT),
    hltIter0IterL3MuonPixelSeedsFromPixelTracks       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks",         "", MYHLT),
    hltIter2IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter2IterL3MuonPixelSeeds",                        "", MYHLT),
    hltIter3IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter3IterL3MuonPixelSeeds",                        "", MYHLT),
    hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.untracked.InputTag("hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks", "", MYHLT),
    hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeeds",              "", MYHLT),
    hltIter3IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter3IterL3FromL1MuonPixelSeeds",                  "", MYHLT),

    hltIterL3OIMuonTrack          = cms.untracked.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity",            "", MYHLT),
    hltIter0IterL3MuonTrack       = cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",           "", MYHLT),
    hltIter2IterL3MuonTrack       = cms.untracked.InputTag("hltIter2IterL3MuonTrackSelectionHighPurity",           "", MYHLT),
    hltIter3IterL3MuonTrack       = cms.untracked.InputTag("hltIter3IterL3MuonTrackSelectionHighPurity",           "", MYHLT),
    hltIter0IterL3FromL1MuonTrack = cms.untracked.InputTag("hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity", "", MYHLT),
    hltIter2IterL3FromL1MuonTrack = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity", "", MYHLT),
    hltIter3IterL3FromL1MuonTrack = cms.untracked.InputTag("hltIter3IterL3FromL1MuonTrackSelectionHighPurity", "", MYHLT),

    # -- generator information
    PUSummaryInfo = cms.untracked.InputTag("addPileupInfo"),
    genEventInfo  = cms.untracked.InputTag("generator"),
    genParticle   = cms.untracked.InputTag("genParticles"),

    # -- L1 tracking + sim-truth association (moved here from the
    #    customizer; the track-association *producers* themselves still
    #    have to be built and put in the Path by the customizer, since
    #    that's genuinely a producer-sequence concern, not an ntupler
    #    parameter)
    L1TrackInputTag = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    l1PrimaryVertex = cms.InputTag("l1tVertexFinderEmulator", "L1VerticesEmulation"),
    associator      = cms.untracked.InputTag("hltTrackAssociatorByHits"),
    trackingParticle = cms.untracked.InputTag("mix", "MergedTrackTruth"),

    trackCollectionNames = cms.untracked.vstring(
        "hltPhase2L3OI",
        "hltIter0Phase2L3FromL1TkMuon",
        "hltIter2Phase2L3FromL1TkMuon",
        "hltPhase2L3IOFromL1",
        "hltPhase2L3MuonsNoID",
        "hltPhase2L3Muons",
    ),
    trackCollectionLabels = cms.untracked.VInputTag(
        cms.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"),
        cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity"),
        cms.InputTag("hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity"),
        cms.InputTag("hltIter2Phase2L3FromL1TkMuonMerged"),
        cms.InputTag("hltPhase2L3MuonsNoIDTracks"),
        cms.InputTag("hltPhase2L3MuonsTracks"),
    ),
    associationLabels = cms.untracked.VInputTag(
        cms.InputTag("AhltPhase2L3OIMuonTrackSelectionHighPurity"),
        cms.InputTag("AhltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity"),
        cms.InputTag("AhltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity"),
        cms.InputTag("AhltIter2Phase2L3FromL1TkMuonMerged"),
        cms.InputTag("AhltPhase2L3MuonsNoID"),
        cms.InputTag("AhltPhase2L3Muons"),
    ),

    trkIsoTags   = cms.untracked.vstring(),
    trkIsoLabels = cms.untracked.VInputTag(),
    pfIsoTags    = cms.untracked.vstring(),
    pfIsoLabels  = cms.untracked.VInputTag(),

    DebugMode     = cms.bool(False),
    SaveAllTracks = cms.bool(True),

    # -- seed-MVA training files/scales for the one seed stage evaluated
    #    by default; swap in the customizer for a different training.
    mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0 = cms.untracked.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Phase2_Iter2FromL1_barrel_v0.xml"),
    mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0 = cms.untracked.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Phase2_Iter2FromL1_endcap_v0.xml"),
    mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B = cms.untracked.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleMean),
    mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B  = cms.untracked.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleStd),
    mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E = cms.untracked.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleMean),
    mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E  = cms.untracked.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleStd),
)

# NOTE: the original cfi/customizer also set a "TkMuonToken" parameter
# (cms.InputTag("L1TkMuons")). It is not read anywhere in MuonHLTNtupler.cc
# -- grep confirms no getParameter/getUntrackedParameter call for that name
# -- so it was already dead configuration before this refactor. Dropped
# here; flagging in case it was meant to feed something that silently
# stopped being wired up at some point.
