# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
# process = customizerFuncForMuonHLTNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms

def customizerForPhase2L3MuonFromGeneralTracksOnly(process, OIFirst = False, newProcessName = "MYHLT", validation = True):

    process.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions = cms.EDProducer("CandidateSeededTrackingRegionsEDProducer",
                                                                                   RegionPSet = cms.PSet(
                                                                                   beamSpot = cms.InputTag("hltOnlineBeamSpot"),
                                                                                   deltaEta = cms.double(0.035),
                                                                                   deltaPhi = cms.double(0.02),
                                                                                   input = cms.InputTag("l1tTkMuonsGmt"),
                                                                                   maxNRegions = cms.int32(10000),
                                                                                   maxNVertices = cms.int32(1),
                                                                                   measurementTrackerName = cms.InputTag(""),
                                                                                   mode = cms.string('BeamSpotSigma'),
                                                                                   nSigmaZBeamSpot = cms.double(4.0),
                                                                                   nSigmaZVertex = cms.double(3.0),
                                                                                   originRadius = cms.double(0.2),
                                                                                   precise = cms.bool(True),
                                                                                   ptMin = cms.double(2.0),
                                                                                   searchOpt = cms.bool(False),
                                                                                   vertexCollection = cms.InputTag("notUsed"),
                                                                                   whereToUseMeasurementTracker = cms.string('Never'),
                                                                                   zErrorBeamSpot = cms.double(24.2),
                                                                                   zErrorVetex = cms.double(0.2)
                                                                               )
    )
    process.hltPhase2L3FromL1TkMuonGeneralTracks = cms.EDProducer("TrackSelectorByRegion",
                                                                  produceMask = cms.bool(False),
                                                                  produceTrackCollection = cms.bool(True),
                                                                  regions = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions"),
                                                                  tracks = cms.InputTag("hltGeneralTracks")
    )

    process.hltPhase2L3FromL1TkMuonGeneralTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
        beamspot = cms.InputTag("hltOnlineBeamSpot"),
        ignoreVertices = cms.bool(False),
        mva = cms.PSet(
            dr_par = cms.PSet(
                d0err = cms.vdouble(0.003, 0.003, 0.003),
                d0err_par = cms.vdouble(0.001, 0.001, 0.001),
                dr_exp = cms.vint32(4, 4, 4),
                dr_par1 = cms.vdouble(0.8, 0.7, 0.6),
                dr_par2 = cms.vdouble(0.6, 0.5, 0.45)
            ),
            dz_par = cms.PSet(
                dz_exp = cms.vint32(4, 4, 4),
                dz_par1 = cms.vdouble(0.9, 0.8, 0.7),
                dz_par2 = cms.vdouble(0.8, 0.7, 0.55)
            ),
            maxChi2 = cms.vdouble(9999.0, 25.0, 16.0),
            maxChi2n = cms.vdouble(2.0, 1.4, 1.2),
            maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
            maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
            maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 15.0),
            maxLostLayers = cms.vint32(3, 2, 2),
            min3DLayers = cms.vint32(3, 3, 3),
            minLayers = cms.vint32(3, 3, 3),
            minNVtxTrk = cms.int32(3),
            minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
            minPixelHits = cms.vint32(0, 0, 3),
            passThroughForAll = cms.bool(False),
            passThroughForDisplaced = cms.bool(False),
            minLayersForDisplaced = cms.int32(4)
        ),
        qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
        src = cms.InputTag("hltPhase2L3FromL1TkMuonGeneralTracks"),
        vertices = cms.InputTag("hltPhase2PixelVertices")
    )

    process.hltPhase2L3FromL1TkMuonGeneralTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
                                                                                    copyExtras = cms.untracked.bool(True),
                                                                                    copyTrajectories = cms.untracked.bool(False),
                                                                                    minQuality = cms.string('highPurity'),
                                                                                    originalMVAVals = cms.InputTag("hltPhase2L3FromL1TkMuonGeneralTrackCutClassifier","MVAValues"),
                                                                                    originalQualVals = cms.InputTag("hltPhase2L3FromL1TkMuonGeneralTrackCutClassifier","QualityMasks"),
                                                                                    originalSource = cms.InputTag("hltPhase2L3FromL1TkMuonGeneralTracks")
    )

    process.hltPhase2L3MuonMerged = cms.EDProducer("TrackListMerger",
                                                   Epsilon = cms.double(-0.001),
                                                   FoundHitBonus = cms.double(5.0),
                                                   LostHitPenalty = cms.double(20.0),
                                                   MaxNormalizedChisq = cms.double(1000.0),
                                                   MinFound = cms.int32(3),
                                                   MinPT = cms.double(0.05),
                                                   ShareFrac = cms.double(0.19),
                                                   TrackProducers = cms.VInputTag(
                                                       "hltPhase2L3OIMuonTrackSelectionHighPurity",
                                                       "hltPhase2L3FromL1TkMuonGeneralTracks",
                                                   ),
                                                   allowFirstHitShare = cms.bool(True),
                                                   copyExtras = cms.untracked.bool(True),
                                                   copyMVA = cms.bool(False),
                                                   hasSelector = cms.vint32(0, 0),
                                                   indivShareFrac = cms.vdouble(1.0, 1.0),
                                                   newQuality = cms.string("confirmed"),
                                                   selectedTrackQuals = cms.VInputTag(
                                                       "hltPhase2L3OIMuonTrackSelectionHighPurity",
                                                       "hltPhase2L3FromL1TkMuonGeneralTrackSelectionHighPurity",
                                                   ),
                                                   setsToMerge = cms.VPSet(cms.PSet(pQual = cms.bool(False), tLists = cms.vint32(0, 1))),
                                                   trackAlgoPriorityOrder = cms.string("hltESPTrackAlgoPriorityOrder"),
                                                   writeOnlyTrkQuals = cms.bool(False),
    )
    process.HLTPhase2L3FromL1TkSequence = cms.Sequence(
        process.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions
        + process.hltPhase2L3FromL1TkMuonGeneralTracks
        + process.hltPhase2L3FromL1TkMuonGeneralTrackCutClassifier
        + process.hltPhase2L3FromL1TkMuonGeneralTrackSelectionHighPurity
    )

    if OIFirst:
        # The alternative HLT Muons sequence (Outside-In first) 
        process.hltPhase2L3MuonFilter = cms.EDProducer("Phase2HLTMuonSelectorForL3",
            l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
            l2MuonsUpdVtx = cms.InputTag("hltL2MuonsFromL1TkMuon:UpdatedAtVtx"),
            l3Tracks = cms.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"),
            IOFirst = cms.bool(False),
            matchingDr = cms.double(0.02),
            applyL3Filters = cms.bool(True),
            MinNhits = cms.int32(1),
            MaxNormalizedChi2 = cms.double(5.0),
            MinNhitsMuons = cms.int32(0),
            MinNhitsPixel = cms.int32(1),
            MinNhitsTracker = cms.int32(6),
            MaxPtDifference = cms.double(999.0),
        )
        process.HLTMuonsSequence = cms.Sequence(
            process.HLTL2MuonsFromL1TkSequence
            + process.HLTPhase2L3OISequence
            + process.hltPhase2L3MuonFilter
            + process.HLTPhase2L3FromL1TkSequence
            + process.HLTPhase2L3MuonsSequence
        )
    else:
        # The default HLT Muons sequence (Inside-Out first)
        process.hltPhase2L3MuonFilter = cms.EDProducer("Phase2HLTMuonSelectorForL3",
            l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
            l2MuonsUpdVtx = cms.InputTag("hltL2MuonsFromL1TkMuon:UpdatedAtVtx"),
            l3Tracks = cms.InputTag("hltPhase2L3FromL1TkMuonGeneralTrackSelectionHighPurity"),
            IOFirst = cms.bool(True),
            matchingDr = cms.double(0.1), # Updated to account for worst resolution. Old value 0.02
            applyL3Filters = cms.bool(True),
            MinNhits = cms.int32(1),
            MaxNormalizedChi2 = cms.double(5.0),
            MinNhitsMuons = cms.int32(0),
            MinNhitsPixel = cms.int32(1),
            MinNhitsTracker = cms.int32(6),
            MaxPtDifference = cms.double(999.0),
            useOfflineSeed = cms.bool(True),
        )
        process.HLTPhase2L3OISequence = cms.Sequence(
            process.hltPhase2L3OISeedsFromL2Muons
            + process.hltPhase2L3OITrackCandidates
            + process.hltPhase2L3OIMuCtfWithMaterialTracks
            + process.hltPhase2L3OIMuonTrackCutClassifier
            + process.hltPhase2L3OIMuonTrackSelectionHighPurity
        )
        process.HLTMuonsSequence = cms.Sequence(
            process.HLTL2MuonsFromL1TkSequence
            + process.HLTPhase2L3FromL1TkSequence
            + process.hltPhase2L3MuonFilter
            + process.HLTPhase2L3OISequence
            + process.HLTPhase2L3MuonsSequence
        )

    process.hltPhase2L3MuonsPFIsodR0p4Dt = cms.EDProducer("MuonHLTPFCandidateIsolationProducer",
        recoCandidateProducer = cms.InputTag("hltPhase2L3MuonCandidates"),
        pfCandidateProducer = cms.InputTag("hltParticleFlowTmp"),
        drMax = cms.double(0.4),
        drVeto = cms.double(0.02),
        drVetoCh = cms.double(0.001),
        minEnergy = cms.double(0.0)
    )

    process.hltPhase2L3MuonsPFChIsoFiltered0p725 = cms.EDFilter("HLTMuonGenericFilter",
        absEtaLowEdges = cms.vdouble(0.0, 1.479),
        candTag = cms.InputTag("hltL3fL1TkSingleMu22L3Filtered24Q"),
        doRhoCorrection = cms.bool(False),
        effectiveAreas = cms.vdouble(0.0, 0.0),
        energyLowEdges = cms.vdouble(0.0),
        l1EGCand = cms.InputTag("hltPhase2L3MuonCandidates"),
        lessThan = cms.bool(True),
        ncandcut = cms.int32(1),
        rhoMax = cms.double(99999999.0),
        rhoScale = cms.double(1.0),
        rhoTag = cms.InputTag(""),
        saveTags = cms.bool(True),
        thrOverE2EB = cms.vdouble(-1.0),
        thrOverE2EE = cms.vdouble(-1.0),
        thrOverEEB = cms.vdouble(0.725),
        thrOverEEE = cms.vdouble(0.725),
        thrRegularEB = cms.vdouble(-1.0),
        thrRegularEE = cms.vdouble(-1.0),
        useEt = cms.bool(True),
        varTag = cms.InputTag("hltPhase2L3MuonsPFIsodR0p4Dt:hltPFRelIsoCh")
    )

    process.hltPhase2L3MuonsPFNhIsoFiltered1p35 = cms.EDFilter("HLTMuonGenericFilter",
        absEtaLowEdges = cms.vdouble(0.0, 1.479),
        candTag = cms.InputTag("hltL3fL1TkSingleMu22L3Filtered24Q"),
        doRhoCorrection = cms.bool(False),
        effectiveAreas = cms.vdouble(0.0, 0.0),
        energyLowEdges = cms.vdouble(0.0),
        l1EGCand = cms.InputTag("hltPhase2L3MuonCandidates"),
        lessThan = cms.bool(True),
        ncandcut = cms.int32(1),
        rhoMax = cms.double(99999999.0),
        rhoScale = cms.double(1.0),
        rhoTag = cms.InputTag(""),
        saveTags = cms.bool(True),
        thrOverE2EB = cms.vdouble(-1.0),
        thrOverE2EE = cms.vdouble(-1.0),
        thrOverEEB = cms.vdouble(1.35),
        thrOverEEE = cms.vdouble(1.35),
        thrRegularEB = cms.vdouble(-1.0),
        thrRegularEE = cms.vdouble(-1.0),
        useEt = cms.bool(True),
        varTag = cms.InputTag("hltPhase2L3MuonsPFIsodR0p4Dt:hltPFRelIsoNh")
    )

    process.hltPhase2L3MuonsPFPhIsoFiltered0p525 = cms.EDFilter("HLTMuonGenericFilter",
        absEtaLowEdges = cms.vdouble(0.0, 1.479),
        candTag = cms.InputTag("hltL3fL1TkSingleMu22L3Filtered24Q"),
        doRhoCorrection = cms.bool(False),
        effectiveAreas = cms.vdouble(0.0, 0.0),
        energyLowEdges = cms.vdouble(0.0),
        l1EGCand = cms.InputTag("hltPhase2L3MuonCandidates"),
        lessThan = cms.bool(True),
        ncandcut = cms.int32(1),
        rhoMax = cms.double(99999999.0),
        rhoScale = cms.double(1.0),
        rhoTag = cms.InputTag(""),
        saveTags = cms.bool(True),
        thrOverE2EB = cms.vdouble(-1.0),
        thrOverE2EE = cms.vdouble(-1.0),
        thrOverEEB = cms.vdouble(0.525),
        thrOverEEE = cms.vdouble(0.525),
        thrRegularEB = cms.vdouble(-1.0),
        thrRegularEE = cms.vdouble(-1.0),
        useEt = cms.bool(True),
        varTag = cms.InputTag("hltPhase2L3MuonsPFIsodR0p4Dt:hltPFRelIsoPh")
    )

    process.HLT_IsoMu24_FromL1TkMuon = cms.Path(
        process.HLTBeginSequence
        + process.hltSingleTkMuon22L1TkMuonFilter
        + process.HLTRawToDigiSequence
        + process.HLTTrackingSequence
        + process.HLTMuonsSequence
        + process.hltL3fL1TkSingleMu22L3Filtered24Q
        + process.HLTLocalrecoSequence
        + process.HLTTICLLocalRecoSequence
        + process.HLTParticleFlowSequence
        + process.hltFixedGridRhoFastjetAllCaloForEGamma
        + process.hltPhase2L3MuonsPFIsodR0p4Dt
        + process.hltPhase2L3MuonsPFChIsoFiltered0p725
        + process.hltPhase2L3MuonsPFNhIsoFiltered1p35
        + process.hltPhase2L3MuonsPFPhIsoFiltered0p525
        + process.HLTEndSequence
    )

    process.HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1 = cms.Path(
        process.HLTBeginSequence
        + process.hltPuppiTauTkMuon4218L1TkFilter
        + process.HLTRawToDigiSequence
        + process.HLTLocalrecoSequence
        + process.HLTTICLLocalRecoSequence
        + process.HLTTrackingSequence
        + process.HLTMuonsSequence
        + process.HLTParticleFlowSequence
        + process.hltParticleFlowRecHitECALUnseeded
        + process.hltParticleFlowClusterECALUncorrectedUnseeded
        + process.hltParticleFlowClusterECALUnseeded
        + process.hltFixedGridRhoFastjetAllCaloForEGamma
        + process.hltPhase2L3MuonsEcalIsodR0p3dRVeto0p000
        + process.hltPhase2L3MuonsHcalIsodR0p3dRVeto0p000
        + process.hltPhase2L3MuonsHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00
        + process.hltL3fL1TkSingleMu18Filtered20
        + process.hltL3crIsoL1TkSingleMu22EcalIso0p41
        + process.hltL3crIsoL1TkSingleMu22HcalIso0p40
        + process.hltL3crIsoL1TkSingleMu22HgcalIso4p70
        + process.HLTPhase2L3MuonGeneralTracksSequence
        + process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p07
        + process.hltL3crIsoL1TkSingleMu22TrkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk
        + process.HLTAK4PFJetsReconstruction
        + process.hltAK4PFJetsForTaus
        + process.HLTPFTauHPS
        + process.HLTHPSDeepTauPFTauSequence
        + process.hltHpsSelectedPFTauLooseTauWPDeepTau
        + process.hltHpsPFTau27LooseTauWPDeepTau
        + process.HLTEndSequence
    )

    process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p4 = cms.EDProducer("L3MuonCombinedRelativeIsolationProducer",
        CaloDepositsLabel = cms.InputTag("notUsed"),
        CaloExtractorPSet = cms.PSet(
            CaloTowerCollectionLabel = cms.InputTag("hltPhase2TowerMakerForAll"),
            ComponentName = cms.string('CaloExtractor'),
            DR_Max = cms.double(0.3),
            DR_Veto_E = cms.double(0.07),
            DR_Veto_H = cms.double(0.1),
            DepositLabel = cms.untracked.string('EcalPlusHcal'),
            Threshold_E = cms.double(0.2),
            Threshold_H = cms.double(0.5),
            Vertex_Constraint_XY = cms.bool(False),
            Vertex_Constraint_Z = cms.bool(False),
            Weight_E = cms.double(1.0),
            Weight_H = cms.double(1.0)
        ),
        CutsPSet = cms.PSet(
            ComponentName = cms.string('SimpleCuts'),
            ConeSizes = cms.vdouble(0.3),
            EtaBounds = cms.vdouble(2.411),
            Thresholds = cms.vdouble(0.4),
            applyCutsORmaxNTracks = cms.bool(False),
            maxNTracks = cms.int32(-1)
        ),
        OutputMuIsoDeposits = cms.bool(True),
        TrackPt_Min = cms.double(-1.0),
        TrkExtractorPSet = cms.PSet(
            BeamSpotLabel = cms.InputTag("hltOnlineBeamSpot"),
            BeamlineOption = cms.string('BeamSpotFromEvent'),
            Chi2Ndof_Max = cms.double(1e+64),
            Chi2Prob_Min = cms.double(-1.0),
            ComponentName = cms.string('PixelTrackExtractor'),
            DR_Max = cms.double(0.3),
            DR_Veto = cms.double(0.005),
            DR_VetoPt = cms.double(0.025),
            DepositLabel = cms.untracked.string('PXLS'),
            Diff_r = cms.double(0.2),
            Diff_z = cms.double(0.25),
            NHits_Min = cms.uint32(0),
            PropagateTracksToRadius = cms.bool(True),
            PtVeto_Min = cms.double(2.0),
            Pt_Min = cms.double(-1.0),
            ReferenceRadius = cms.double(6.0),
            VetoLeadingTrack = cms.bool(True),
            inputTrackCollection = cms.InputTag("hltPhase2PixelTracks")
        ),
        UseCaloIso = cms.bool(False),
        UseRhoCorrectedCaloDeposits = cms.bool(False),
        inputMuonCollection = cms.InputTag("hltPhase2L3MuonCandidates"),
        printDebug = cms.bool(False)
    )

    process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_FromL1TkMuon = cms.Path(
        process.HLTBeginSequence
        + process.hltDoubleTkMuon157L1TkMuonFilter
        + process.hltDoubleMuon7DZ1p0
        + process.HLTRawToDigiSequence
        + process.HLTTrackingSequence
        + process.HLTMuonsSequence
        + process.hltL3fL1DoubleMu155fPreFiltered8
        + process.hltL3fL1DoubleMu155fFiltered17
        + process.HLTPhase2L3MuonGeneralTracksSequence
        + process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p4
        + process.hltDiMuon178RelTrkIsoFiltered0p4
        + process.hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2
        + process.HLTEndSequence
    )

    process.HLT_Mu37_Mu27_FromL1TkMuon = cms.Path(
        process.HLTBeginSequence
        + process.hltDoubleMuon7DZ1p0
        + process.HLTRawToDigiSequence
        + process.HLTTrackingSequence
        + process.HLTMuonsSequence
        + process.hltL3fL1DoubleMu155fPreFiltered27
        + process.hltL3fL1DoubleMu155fFiltered37
        + process.HLTEndSequence
    )

    process.HLT_Mu50_FromL1TkMuon = cms.Path(
        process.HLTBeginSequence
        + process.hltSingleTkMuon22L1TkMuonFilter
        + process.HLTRawToDigiSequence
        + process.HLTTrackingSequence
        + process.HLTMuonsSequence
        + process.hltL3fL1TkSingleMu22L3Filtered50Q
        + process.HLTEndSequence
    )

    process.HLT_TriMu_10_5_5_DZ_FromL1TkMuon = cms.Path(
        process.HLTBeginSequence
        + process.hltTripleMuon3DZ1p0
        + process.hltTripleMuon3DR0
        + process.HLTRawToDigiSequence
        + process.HLTTrackingSequence
        + process.HLTMuonsSequence
        + process.hltL3fL1TkTripleMu533PreFiltered555
        + process.hltL3fL1TkTripleMu533L3Filtered1055
        + process.hltL3fL1TkTripleMu533L31055DZFiltered0p2
        + process.HLTEndSequence
    )
    
    return process

