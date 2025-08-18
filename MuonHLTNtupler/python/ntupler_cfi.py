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
                             
                             pfClusterProducer_ecal = cms.untracked.InputTag( "hltParticleFlowClusterECALUnseeded", "", "MYHLT"),

                             rho_ECAL = cms.untracked.InputTag("hltFixedGridRhoFastjetAllCaloForEGamma", "", "MYHLT"),
                             drMax_ECAL = cms.double( 0.3 ),
                             drVetoBarrel_ECAL = cms.double( 0.05 ),
                             drVetoEndcap_ECAL = cms.double( 0.05 ),
                             etaStripBarrel_ECAL = cms.double( 0.0 ),
                             etaStripEndcap_ECAL = cms.double( 0.0 ),
                             energyBarrel_ECAL = cms.double( 0.0 ),
                             energyEndcap_ECAL = cms.double( 0.0 ),
                             
                             pfClusterProducerHCAL = cms.untracked.InputTag( "hltParticleFlowClusterHCAL", "", "MYHLT"),

                             rho_HCAL = cms.untracked.InputTag("hltFixedGridRhoFastjetAllCaloForEGamma", "", "MYHLT"),
                             drMax_HCAL = cms.double( 0.3 ),
                             drVetoBarrel_HCAL = cms.double( 0.1 ),
                             drVetoEndcap_HCAL = cms.double( 0.1 ),
                             etaStripBarrel_HCAL = cms.double( 0.0 ),
                             etaStripEndcap_HCAL = cms.double( 0.0 ),
                             energyBarrel_HCAL = cms.double( 0.0 ),
                             energyEndcap_HCAL = cms.double( 0.0 ),
                             
                             layerClusterProducer_HGCAL = cms.untracked.InputTag( "hltHgcalMergeLayerClusters", "", "MYHLT"),
                             hgcalLayerClustersTime = cms.untracked.InputTag( "hltHgcalMergeLayerClusters:timeLayerCluster", "", "MYHLT"),

                             drVetoHad_HGCAL = cms.double(0.02),
                             drVetoEM_HGCAL = cms.double(0.0),
                             drMax_HGCAL = cms.double(0.2),

                             inputTrackCollection = cms.untracked.InputTag( "hltPhase2L3MuonGeneralTracks", "", "MYHLT"), # hltPhase2L3MuonGeneralTracksMTDTExtendedVtx
                             #trackAssocSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx"), # hltPhase2L3MuonGeneralTracks?

                             ### From TrackMTDExtender
                             t0Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackt0"),
                             tmtdSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracktmtd"),
                             sigmat0Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracksigmat0"),
                             sigmatmtdSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracksigmatmtd"),
                             tofPiSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackTofPi"),
                             tofKSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackTofK"),
                             tofPSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackTofP"),
                             sigmatofpiSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackSigmaTofPi"),
                             sigmatofkSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackSigmaTofK"),
                             sigmatofpSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackSigmaTofP"),
                             btlMatchChi2Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:btlMatchChi2"),
                             etlMatchChi2Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:etlMatchChi2"),
                             btlMatchTimeChi2Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:btlMatchTimeChi2"),
                             etlMatchTimeChi2Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:etlMatchTimeChi2"),
                             npixBarrelSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:npixBarrel"),
                             npixEndcapSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:npixEndcap"),
                             trackOutermostHitPositionSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackOutermostHitPosition"),
                             trackpSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackp"),
                             trackBetaSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackBeta"),
                             trackPathLengthSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackPathLength"),
                             trackmtdposSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackmtdpos"),
                             
                             # From TofProducer
                             t0PID = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtxTof:t0"),
                             t0SafePID = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtxTof:t0safe"),
                             sigmat0SafePID = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtxTof:sigmat0safe"),
                             probPi = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtxTof:probPi"),
                             probK = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtxTof:probK"),
                             probP = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtxTof:probP"),

                            # From MTDTrackQualityProducer
                                                            
                             trackMVAQual = cms.InputTag("hltPhase2GeneralTracksMTDExtendedVtxMVATrackQuality:mtdQualMVA"),

                             Diff_r = cms.double( 0.1 ),
                             Diff_z = cms.double( 0.2 ),
                             DR_Max = cms.double( 0.3 ),
                             DR_Veto = cms.double( 0.01 ),
                             NHits_Min = cms.uint32( 0 ),
                             Chi2Ndof_Max = cms.double( 1.0E64 ),
                             Chi2Prob_Min = cms.double( -1.0 ),
                             Pt_Min = cms.double( -1.0 ),
                             
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
