import FWCore.ParameterSet.Config as cms

ntuplerBase = cms.EDAnalyzer("MuonHLTNtupler",
	# -- information stored in edm file
    hltGeneralTracks = cms.InputTag("hltGeneralTracks"),
    trajTrackAssoc = cms.InputTag("hltGeneralTracks"),

    # trackerHitAssociator
    Quality_SimToReco = cms.double(0.5),
    associateRecoTracks = cms.bool(True),
    UseGrouped = cms.bool(True),
    associatePixel = cms.bool(True),
    associateStrip = cms.bool(True),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
    usePhase2Tracker = cms.bool(True),    
    ROUList = cms.VInputTag(
        cms.InputTag("g4SimHits","TrackerHitsTIBLowTof"),
        cms.InputTag("g4SimHits","TrackerHitsTIBHighTof"),
        cms.InputTag("g4SimHits","TrackerHitsTIDLowTof"),
        cms.InputTag("g4SimHits","TrackerHitsTIDHighTof"),
        cms.InputTag("g4SimHits","TrackerHitsTOBLowTof"),
        cms.InputTag("g4SimHits","TrackerHitsTOBHighTof"),
        cms.InputTag("g4SimHits","TrackerHitsTECLowTof"),
        cms.InputTag("g4SimHits","TrackerHitsTECHighTof"),
        cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),
        cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"),
        cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"),
        cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof"),
    ),
    UseSplitting = cms.bool(True),
    UsePixels = cms.bool(True),
    ThreeHitTracksAreSpecial = cms.bool(True),
    AbsoluteNumberOfHits = cms.bool(False),
    Purity_SimToReco = cms.double(0.75),
    Cut_RecoToSim = cms.double(0.75),
    SimToRecoDenominator = cms.string('sim'), ##"reco"
    simHitTpMapTag = cms.InputTag("simHitTPAssocProducer"),
    phase2TrackerSimLinkSrc  = cms.InputTag("simSiPixelDigis","Tracker"),

	simTrackLabel = cms.untracked.InputTag("g4SimHits"),

	dtSegments = cms.InputTag('hltDt4DSegments'),
    cscSegments = cms.InputTag('hltCscSegments'),
    gemSegments = cms.InputTag('hltGemSegments'),
    me0Segments = cms.InputTag('hltMe0Segments'),

	dtDigiSimLinkLabel = cms.untracked.InputTag('simMuonDTDigis'),
	cscStripDigiSimLinkLabel = cms.untracked.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
	rpcDigiSimLinkLabel = cms.untracked.InputTag("simMuonRPCDigis"),
    gemDigiSimLinks = cms.InputTag("simMuonGEMDigis"),
    me0DigiSimLinks = cms.InputTag("simMuonME0Digis"),
)
