import FWCore.ParameterSet.Config as cms

ntuplerBase = cms.EDAnalyzer("MuonHLTNtupler",
                             
                             inputMuonCollection = cms.InputTag( "hltPhase2L3MuonCandidates", "", "MYHLT"),
                             inputMuonFilterCollection = cms.InputTag( "hltL3fL1TkSingleMu22L3Filtered24Q", "", "MYHLT"),

                             #muontrkAssSrc  = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx", "", "MYHLT"),

                             muonbtlMatchChi2 = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:btlMatchChi2"),
                             muonetlMatchChi2 = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:etlMatchChi2"),
                             muonbtlMatchTimeChi2 = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:btlMatchTimeChi2"),
                             muonetlMatchTimeChi2 = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:etlMatchTimeChi2"),
                             muonnpixBarrel = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:npixBarrel"),
                             muonnpixEndcap = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:npixEndcap"),
                             muonTrackOutermostHitPosition = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTrackOutermostHitPosition"),
                             muonTrackp = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTrackp"),
                             muonTrackBeta = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTrackBeta"),
                             muonTrackt0 = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTrackt0"),
                             muonTracksigmat0 = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTracksigmat0"),
                             muonTrackPathLength = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTrackPathLength"),
                             muonTracktmtd = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTracktmtd"),
                             muonTracksigmatmtd = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTracksigmatmtd"),
                             muonTrackmtdpos = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTrackmtdpos"),
                             muonTrackTofMu = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTrackTofMu"),
                             muonTrackSigmaTofMu = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTrackSigmaTofMu"),
                             
                             # PF candidates
                             
                             pfCandidateProducer = cms.InputTag("hltParticleFlowTmp"),
                             drMaxPF = cms.double(0.4),
                             drVetoPF = cms.double(0.01),
                             drVetoPFCh = cms.double(0.0001),
                             minEnergyPF = cms.double(0.0),

                             rhoProducer = cms.InputTag("hltFixedGridRhoFastjetAllCaloForEGamma"),
                             trackBtlMatchChi2 = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:btlMatchChi2"),
                             trackEtlMatchChi2 = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:etlMatchChi2"),
                             trackBtlMatchTimeChi2 = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:btlMatchTimeChi2"),
                             trackEtlMatchTimeChi2 = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:etlMatchTimeChi2"),
                             trackNpixBarrel = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:npixBarrel"),
                             trackNpixEndcap = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:npixEndcap"),
                             trackOutermostHitPosition = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackOutermostHitPosition"),
                             trackp = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackp"),
                             trackBeta = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackBeta"),
                             trackt0 = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackt0"),
                             tracksigmat0 = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracksigmat0"),
                             trackPathLength = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackPathLength"),
                             tracktmtd = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracktmtd"),
                             tracksigmatmtd = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracksigmatmtd"),
                             trackmtdpos = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackmtdpos"),
                             trackTofPi = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackTofPi"),
                             trackSigmaTofPi = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackSigmaTofPi"),
                             trackTofK = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackTofK"),
                             trackSigmaTofK = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackSigmaTofK"),
                             trackTofP = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackTofP"),
                             trackSigmaTofP = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackSigmaTofP"),
                             
	                     # -- information stored in edm file
	                     triggerResults    = cms.untracked.InputTag("TriggerResults::MYHLT"),
	                     triggerEvent      = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"),
	                     offlineLumiScaler = cms.untracked.InputTag("scalersRawToDigi"),
	                     offlineVertex     = cms.untracked.InputTag("offlinePrimaryVertices"),
	                     offlineMuon       = cms.untracked.InputTag("muons"),
                             
	                     # -- newly created objects by HLT rerun
	                     # -- new process name = "MYHLT"
	                     myTriggerResults = cms.untracked.InputTag("TriggerResults",       "", "MYHLT"),
	                     myTriggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "MYHLT"),
	                     lumiScaler       = cms.untracked.InputTag("hltScalersRawToDigi",  "", "MYHLT"),
                             
	                     L1Muon = cms.untracked.InputTag("hltGmtStage2Digis",       "Muon", "MYHLT"), # -- if L1 emulator is used
	                     # L1Muon = cms.untracked.InputTag("gmtStage2Digis",          "Muon", "RECO"), # -- if L1 is not emulated
	                     L2Muon = cms.untracked.InputTag("hltL2MuonCandidates",     "",     "MYHLT"),
	                     L3Muon = cms.untracked.InputTag("hltPhase2L3MuonCandidates", "",     "MYHLT"),
	                     TkMuon = cms.untracked.InputTag("hltHighPtTkMuonCands",    "",     "MYHLT"),
                             
	                     iterL3OI        = cms.untracked.InputTag("hltL3MuonsIterL3OI",                   "", "MYHLT"),
	                     iterL3IOFromL2  = cms.untracked.InputTag("hltL3MuonsIterL3IO",                   "", "MYHLT"),
	                     iterL3FromL2    = cms.untracked.InputTag("hltIterL3MuonsFromL2LinksCombination", "", "MYHLT"),
	                     iterL3IOFromL1  = cms.untracked.InputTag("hltIter3IterL3FromL1MuonMerged",       "", "MYHLT"),
	                     iterL3MuonNoID  = cms.untracked.InputTag("hltIterL3MuonsNoID",                   "", "MYHLT"),
	                     iterL3Muon      = cms.untracked.InputTag("hltIterL3Muons",                       "", "MYHLT"),
                             
	                     hltIterL3MuonTrimmedPixelVertices                 = cms.untracked.InputTag("hltIterL3MuonTrimmedPixelVertices",                   "", "MYHLT"),
	                     hltIterL3FromL1MuonTrimmedPixelVertices           = cms.untracked.InputTag("hltIterL3FromL1MuonTrimmedPixelVertices",             "", "MYHLT"),
                             
	                     doMVA  = cms.bool(False),
	                     doSeed = cms.bool(False),
                             
	                     hltIterL3OISeedsFromL2Muons                       = cms.untracked.InputTag("hltIterL3OISeedsFromL2Muons",                         "", "MYHLT"),
	                     hltIter0IterL3MuonPixelSeedsFromPixelTracks       = cms.untracked.InputTag("hltIter0IterL3MuonPixelSeedsFromPixelTracks",         "", "MYHLT"),
	                     hltIter2IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter2IterL3MuonPixelSeeds",                        "", "MYHLT"),
	                     hltIter3IterL3MuonPixelSeeds                      = cms.untracked.InputTag("hltIter3IterL3MuonPixelSeeds",                        "", "MYHLT"),
	                     hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.untracked.InputTag("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks",   "", "MYHLT"),
	                     hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter2IterL3FromL1MuonPixelSeeds",                  "", "MYHLT"),
	                     hltIter3IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter3IterL3FromL1MuonPixelSeeds",                  "", "MYHLT"),
                             
	                     hltIterL3OIMuonTrack          = cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity",       "", "MYHLT"),
	                     hltIter0IterL3MuonTrack       = cms.untracked.InputTag("hltIter0IterL3MuonTrackSelectionHighPurity",       "", "MYHLT"),
	                     hltIter2IterL3MuonTrack       = cms.untracked.InputTag("hltIter2IterL3MuonTrackSelectionHighPurity",       "", "MYHLT"),
	                     hltIter3IterL3MuonTrack       = cms.untracked.InputTag("hltIter3IterL3MuonTrackSelectionHighPurity",       "", "MYHLT"), 
	                     hltIter0IterL3FromL1MuonTrack = cms.untracked.InputTag("hltIter0IterL3FromL1MuonTrackSelectionHighPurity",       "", "MYHLT"),
	                     hltIter2IterL3FromL1MuonTrack = cms.untracked.InputTag("hltIter2IterL3FromL1MuonTrackSelectionHighPurity",       "", "MYHLT"),
	                     hltIter3IterL3FromL1MuonTrack = cms.untracked.InputTag("hltIter3IterL3FromL1MuonTrackSelectionHighPurity",       "", "MYHLT"),
                             
	                     # -- generator information
	                     PUSummaryInfo = cms.untracked.InputTag("addPileupInfo"),
	                     genEventInfo = cms.untracked.InputTag("generator"),
	                     genParticle = cms.untracked.InputTag("genParticles"),
)
