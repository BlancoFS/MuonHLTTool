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
# -----------------------------------------------------------------------
muonCollections = cms.untracked.VPSet(
    cms.PSet(name=cms.string("L1TkMuon"), type=cms.string("L1TkMuon"),
            label=cms.InputTag("l1tTkMuonsGmt", "", MYHLT)),

    cms.PSet(name=cms.string("L2Muon"), type=cms.string("TrackView"),
            label=cms.InputTag("hltL2MuonsFromL1TkMuon", "", MYHLT)),

    cms.PSet(name=cms.string("L3MuonNoId"), type=cms.string("RecoMuon"),
             label=cms.InputTag("hltPhase2L3MuonsNoID", "", MYHLT)),

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

    muonCollections = muonCollections,
    triggerPaths    = triggerPaths,
    triggerFilters  = triggerFilters,
    maxDR           = maxDR,

    hltIterL3MuonTrimmedPixelVertices       = cms.untracked.InputTag("hltIterL3MuonTrimmedPixelVertices",       "", MYHLT),
    hltIterL3FromL1MuonTrimmedPixelVertices = cms.untracked.InputTag("hltPhase2L3FromL1TkMuonTrimmedPixelVertices", "", MYHLT),

    # -- generator information
    PUSummaryInfo = cms.untracked.InputTag("addPileupInfo"),
    genEventInfo  = cms.untracked.InputTag("generator"),
    genParticle   = cms.untracked.InputTag("genParticles"),

    DebugMode     = cms.bool(False),
)
