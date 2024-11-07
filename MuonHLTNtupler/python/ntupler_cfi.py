import FWCore.ParameterSet.Config as cms


ntuplerBase = cms.EDAnalyzer("MuonHLTNtupler",

                             useHF = cms.bool( False ),
                             drMax_HCAL = cms.double( 0.3 ),

                             drVetoBarrel_HCAL = cms.double( 0.1 ),
                             drVetoEndcap_HCAL = cms.double( 0.1 ),
                             etaStripBarrel_HCAL = cms.double( 0.0 ),
                             etaStripEndcap_HCAL = cms.double( 0.0 ),
                             energyBarrel_HCAL = cms.double( 0.0 ),
                             energyEndcap_HCAL = cms.double( 0.0 ),
                             doRhoCorrection_HCAL = cms.bool( True ),
                             rhoMax_HCAL = cms.double( 9.9999999E7 ),
                             rhoScale_HCAL = cms.double( 1.0 ),
                             effectiveAreas_HCAL = cms.vdouble( 0.227, 0.372 ),
                             absEtaLowEdges_HCAL = cms.vdouble( 0.0, 1.479 ),

                             drMax_ECAL = cms.double( 0.3 ),
                             drVetoBarrel_ECAL = cms.double( 0.05 ),
                             drVetoEndcap_ECAL = cms.double( 0.05 ),
                             etaStripBarrel_ECAL = cms.double( 0.0 ),
                             etaStripEndcap_ECAL = cms.double( 0.0 ),
                             energyBarrel_ECAL = cms.double( 0.0 ),
                             energyEndcap_ECAL = cms.double( 0.0 ),
                             doRhoCorrection_ECAL = cms.bool( True ),
                             rhoMax_ECAL = cms.double( 9.9999999E7 ),
                             rhoScale_ECAL = cms.double( 1.0 ),
                             effectiveAreas_ECAL = cms.vdouble( 0.35, 0.193 ),
                             absEtaLowEdges_ECAL = cms.vdouble( 0.0, 1.479 ),

                             useEt = cms.bool( True ),

                             DepositLabel = cms.untracked.string( "PXLS" ),
                             Diff_r = cms.double( 0.1 ),
                             Diff_z = cms.double( 0.2 ),
                             DR_Max = cms.double( 0.3 ),
                             DR_Veto = cms.double( 0.01 ),
                             NHits_Min = cms.uint32( 0 ),
                             Chi2Ndof_Max = cms.double( 1.0E64 ),
                             Chi2Prob_Min = cms.double( -1.0 ),
                             Pt_Min = cms.double( -1.0 ),
                             
	                     # -- information stored in edm file
	                     triggerResults    = cms.untracked.InputTag("TriggerResults::HLT"),
	                     triggerEvent      = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
	                     offlineLumiScaler = cms.untracked.InputTag("scalersRawToDigi"),
	                     offlineVertex     = cms.untracked.InputTag("offlinePrimaryVertices"),
	                     offlineMuon       = cms.untracked.InputTag("muons"),
	                     beamSpot          = cms.untracked.InputTag("hltOnlineBeamSpot"),

	                     # -- newly created objects by HLT rerun
	                     # -- new process name = "MYHLT"
	                     myTriggerResults = cms.untracked.InputTag("TriggerResults",       "", "MYHLT"),
	                     myTriggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "MYHLT"),
	                     lumiScaler       = cms.untracked.InputTag("hltScalersRawToDigi",  "", "MYHLT"),

                             # -- generator information
                             PUSummaryInfo = cms.untracked.InputTag("addPileupInfo"),
                             genEventInfo = cms.untracked.InputTag("generator"),
                             genParticle = cms.untracked.InputTag("genParticles"),

	                     rho_ECAL = cms.untracked.InputTag("hltFixedGridRhoFastjetECALMFForMuons", "", "MYHLT"),
	                     rho_HCAL = cms.untracked.InputTag("hltFixedGridRhoFastjetHCAL",           "", "MYHLT"),

                             pfClusterProducerHCAL = cms.untracked.InputTag( "hltParticleFlowClusterHCAL", "", "MYHLT"),
                             rhoProducer_HCAL = cms.untracked.InputTag( "hltFixedGridRhoFastjetHCAL", "", "MYHLT"),
                             pfClusterProducerHFEM = cms.untracked.InputTag( "", "", "MYHLT" ),
                             pfClusterProducerHFHAD = cms.untracked.InputTag( "", "", "MYHLT" ),

                             pfClusterProducer = cms.untracked.InputTag( "hltParticleFlowClusterECALForMuonsMF", "", "MYHLT"),
                             inputTrackCollection = cms.untracked.InputTag( "hltIter0L3MuonTrackSelectionHighPurity", "", "MYHLT"),
                             inputMuonCollection = cms.InputTag( "hltL3fL1sSingleMu22L1f0L2f10QL3Filtered24Q", "", "MYHLT"),

                             useSimpleGeometry = cms.bool( True ),
                             useStation2 = cms.bool( True ),
                             fallbackToME1 = cms.bool( False ),
                             cosmicPropagationHypothesis = cms.bool( False ),
                             useMB2InOverlap = cms.bool( False ),
                             useTrack = cms.string( "tracker" ),
                             useState = cms.string( "atVertex" ),
                             propagatorAlong = cms.ESInputTag( "","hltESPSteppingHelixPropagatorAlong" ),
                             propagatorAny = cms.ESInputTag( "","SteppingHelixPropagatorAny" ),
                             propagatorOpposite = cms.ESInputTag( "","hltESPSteppingHelixPropagatorOpposite" ),

)
