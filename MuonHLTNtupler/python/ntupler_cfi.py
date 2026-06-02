import FWCore.ParameterSet.Config as cms

ntuplerBase = cms.EDAnalyzer("MuonHLTNtupler",
	# -- information stored in edm file
    hltGeneralTracks = cms.InputTag("hltGeneralTracks"),
    trajTrackAssoc = cms.InputTag("hltGeneralTracks"),

	simTrackLabel = cms.untracked.InputTag("g4SimHits"),

	dtSegments = cms.untracked.InputTag('hltDt4DSegments'),
    cscSegments = cms.untracked.InputTag('hltCscSegments'),
    gemSegments = cms.untracked.InputTag('hltGemSegments'),
    me0Segments = cms.untracked.InputTag('hltMe0Segments'),

	dtDigiSimLinkLabel = cms.untracked.InputTag('simMuonDTDigis'),
	cscStripDigiSimLinkLabel = cms.untracked.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
	rpcDigiSimLinkLabel = cms.untracked.InputTag("simMuonRPCDigis"),
    gemDigiSimLinks = cms.untracked.InputTag("simMuonGEMDigis"),
    me0DigiSimLinks = cms.untracked.InputTag("simMuonME0Digis"),
)
