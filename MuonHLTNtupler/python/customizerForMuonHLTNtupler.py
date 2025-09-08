# -- custoimzer for ntupler that can be added to the HLT configuration for re-running HLT
# -- add two lines in the HLT config.:
# from MuonHLTTool.MuonHLTNtupler.customizerForMuonHLTNtupler import *
# process = customizerFuncForMuonHLTNtupler(process, "MYHLT")

import FWCore.ParameterSet.Config as cms
# from RecoMuon.TrackerSeedGenerator.mvaScale import *

from RecoMTD.TrackExtender.PropagatorWithMaterialForMTD_cfi import *
from RecoMTD.TrackExtender.trackExtenderWithMTDBase_cfi import *
from RecoMTD.TransientTrackingRecHit.MTDTransientTrackingRecHitBuilder_cfi import *
from SimFastTiming.FastTimingCommon.mtdDigitizer_cfi import mtdDigitizer

from TrackingTools.TrackFitters.KFTrajectoryFitter_cfi import *
from TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi import *
from TrackingTools.TrackFitters.KFTrajectorySmoother_cfi import *

from TrackingTools.TrackRefitter.TracksToTrajectories_cff import *
from RecoTracker.FinalTrackSelectors.TrackCollectionMerger_cfi import *
from RecoTracker.FinalTrackSelectors.trackAlgoPriorityOrder_cfi import trackAlgoPriorityOrder

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *
from RecoMuon.TransientTrackingRecHit.MuonTransientTrackingRecHitBuilder_cfi import *

# -- Std Transform -- #
PU200_Barrel_NThltIter2FromL1_ScaleMean     = [0.00033113700731766336, 1.6825601468762878e-06, 1.790932122524803e-06, 0.010534608406382916, 0.005969459957330139, 0.0009605022254971113, 0.04384189672781466, 7.846741237608237e-05, 0.40725050850004824, 0.41125151617410227, 0.39815551065544846]
PU200_Barrel_NThltIter2FromL1_ScaleStd      = [0.0006042948363798624, 2.445644111872427e-06, 3.454992543447134e-06, 0.09401581628887255, 0.7978806947573766, 0.4932933044535928, 0.04180518265631776, 0.058296511682094855, 0.4071857009373577, 0.41337782307392973, 0.4101160349549534]
PU200_Endcap_NThltIter2FromL1_ScaleMean     = [0.00022658482374555603, 5.358921973784045e-07, 1.010003713549798e-06, 0.0007886873612224615, 0.001197730548842408, -0.0030252353426003594, 0.07151944804171254, -0.0006940626775109026, 0.20535152195939896, 0.2966816533783824, 0.28798220230180455]
PU200_Endcap_NThltIter2FromL1_ScaleStd      = [0.0003857726789049956, 1.4853721474087994e-06, 6.982997036736564e-06, 0.04071340757666084, 0.5897606560095399, 0.33052121398064654, 0.05589386786541949, 0.08806273533388546, 0.3254586902665612, 0.3293354496231377, 0.3179899794578072]


def customizerFuncForMuonGeneralTrackerExtender(process, doMuon = True, newProcessName = "MYHLT"):
    
    process.Chi2EstimatorForRefit = Chi2MeasurementEstimator.clone(
        ComponentName = 'Chi2EstimatorForRefit',
        MaxChi2 = 100000.0,
        nSigma = 3.0
    )
    
    process.KFFitterForRefitInsideOut = KFTrajectoryFitter.clone(
        ComponentName = 'KFFitterForRefitInsideOut',
        Propagator = 'SmartPropagatorAnyRK',
        Updator = 'KFUpdator',
        Estimator = 'Chi2EstimatorForRefit',
        minHits = 3
    )

    process.KFSmootherForRefitInsideOut = KFTrajectorySmoother.clone(
        ComponentName = 'KFSmootherForRefitInsideOut',
        Propagator = 'SmartPropagatorAnyRK',
        Updator = 'KFUpdator',
        Estimator = 'Chi2EstimatorForRefit',
        errorRescaling = 100.0,
        minHits = 3
    )
    
    process.hltPhase2L3MuonMtdUncalibratedRecHits = cms.EDProducer(
        "MTDUncalibratedRecHitProducer",
        barrel = cms.PSet(
            algoName = cms.string("BTLUncalibRecHitAlgo"),
            adcNbits = mtdDigitizer.barrelDigitizer.ElectronicsSimulation.adcNbits,
            adcSaturation = mtdDigitizer.barrelDigitizer.ElectronicsSimulation.adcSaturation_MIP,
            toaLSB_ns = mtdDigitizer.barrelDigitizer.ElectronicsSimulation.toaLSB_ns,
            timeResolutionInNs = cms.string("0.308*pow(x,-0.4175)"), # [ns]
            timeCorr_p0 = cms.double( 2.21103),
            timeCorr_p1 = cms.double(-0.933552),
            timeCorr_p2 = cms.double( 0.),
            c_LYSO = cms.double(13.846235)     # in unit cm/ns
        ),
        endcap = cms.PSet(
            algoName      = cms.string("ETLUncalibRecHitAlgo"),
            adcNbits      = mtdDigitizer.endcapDigitizer.ElectronicsSimulation.adcNbits,
            adcSaturation = mtdDigitizer.endcapDigitizer.ElectronicsSimulation.adcSaturation_MIP,
            toaLSB_ns     = mtdDigitizer.endcapDigitizer.ElectronicsSimulation.toaLSB_ns,
            tofDelay      = mtdDigitizer.endcapDigitizer.DeviceSimulation.tofDelay,
            timeResolutionInNs = cms.string("0.0370"), # [ns]
            timeCorr_p0 = cms.double(0.974683),
            timeCorr_p1 = cms.double(-0.237274),
            timeCorr_p2 = cms.double(0.021455),
            timeCorr_p3 = cms.double(-0.000727429)
        ),
        barrelDigis = cms.InputTag('mix:FTLBarrel'),
        endcapDigis = cms.InputTag('mix:FTLEndcap'),
        BarrelHitsName = cms.string('FTLBarrel'),
        EndcapHitsName = cms.string('FTLEndcap')
    )
    
    process.hltPhase2L3MuonMtdRecHits = cms.EDProducer(
        "MTDRecHitProducer", # MTDRecHitProducer
        barrel = cms.PSet(
            algoName = cms.string("MTDRecHitAlgo"),
            thresholdToKeep = cms.double(1.),          # MeV
            calibrationConstant = cms.double(0.03125), # MeV/pC
        ),
        endcap = cms.PSet(
            algoName = cms.string("MTDRecHitAlgo"),
            thresholdToKeep = cms.double(0.0425),    # MeV
            calibrationConstant = cms.double(0.085), # MeV/MIP
        ),
        barrelUncalibratedRecHits = cms.InputTag('hltPhase2L3MuonMtdUncalibratedRecHits:FTLBarrel'),
        endcapUncalibratedRecHits = cms.InputTag('hltPhase2L3MuonMtdUncalibratedRecHits:FTLEndcap'),
        BarrelHitsName = cms.string('FTLBarrel'),
        EndcapHitsName = cms.string('FTLEndcap'),
    )

    process.hltPhase2L3MuonMTDClusters = cms.EDProducer(
        "MTDClusterProducer",
        srcBarrel = cms.InputTag("hltPhase2L3MuonMtdRecHits:FTLBarrel"),
        srcEndcap = cms.InputTag("hltPhase2L3MuonMtdRecHits:FTLEndcap"),        
        BarrelClusterName = cms.string("FTLBarrel"),
        EndcapClusterName = cms.string("FTLEndcap"),
        ClusterMode = cms.string("MTDThresholdClusterizer"),
    )
    
    process.hltPhase2L3MuonMTDRecHits = cms.EDProducer(
        "MTDTrackingRecHitProducer",
        barrelClusters = cms.InputTag("hltPhase2L3MuonMTDClusters:FTLBarrel"),
        endcapClusters = cms.InputTag("hltPhase2L3MuonMTDClusters:FTLEndcap"),
    )    
    
    ### Prepare generalTracks - redefine to access trajectories
    #process.hltPhase2L3MuonInitialStepTracks.TrajectoryInEvent = True
    #process.hltPhase2L3MuonHighPtTripletStepTracks.TrajectoryInEvent = True
    #process.hltPhase2L3MuonInitialStepTracksSelectionHighPurity.copyTrajectories = True
    #process.hltPhase2L3MuonHighPtTripletStepTracksSelectionHighPurity.copyTrajectories = True
    
    #process.hltPhase2L3MuonGeneralTracks = TrackCollectionMerger.clone(
    #    trackProducers   = ["hltPhase2L3MuonInitialStepTracks", "hltPhase2L3MuonHighPtTripletStepTracks"],
    #    inputClassifiers = ["hltPhase2L3MuonInitialStepTrackCutClassifier", "hltPhase2L3MuonHighPtTripletStepTrackCutClassifier"],
    #    foundHitBonus  = 100.0,
    #    lostHitPenalty =   1.0,
    #    minQuality = cms.string('highPurity'),
    #    copyExtras = cms.untracked.bool(True),
    #    copyTrajectories = cms.untracked.bool(True),
    #)

    process.hltInitialStepTracks.TrajectoryInEvent = True
    process.hltHighPtTripletStepTracks.TrajectoryInEvent = True
    process.hltInitialStepTrackSelectionHighPurity.copyTrajectories = True
    process.hltHighPtTripletStepTrackSelectionHighPurity.copyTrajectories = True
    
    process.hltGeneralTracks = TrackCollectionMerger.clone(
        trackProducers   = ["hltInitialStepTracks", "hltHighPtTripletStepTracks"],
        inputClassifiers = ["hltInitialStepTrackCutClassifier", "hltHighPtTripletStepTrackCutClassifier"],
        foundHitBonus  = 5.0,
        lostHitPenalty = 5.0,
        minQuality = cms.string('highPurity'),
        copyExtras = cms.untracked.bool(True),
        copyTrajectories = cms.untracked.bool(True),
    )
    
    process.hltPhase2GeneralTracksMTDTExtendedVtx = cms.EDProducer(
        "TrackExtenderWithMTD",
        #tracksSrc = cms.InputTag("hltPhase2L3MuonGeneralTracks"),
	#trjtrkAssSrc = cms.InputTag("hltPhase2L3MuonGeneralTracks"),
        tracksSrc = cms.InputTag("hltGeneralTracks"),
        trjtrkAssSrc = cms.InputTag("hltGeneralTracks"),
        hitsSrc = cms.InputTag("hltPhase2L3MuonMTDRecHits"),	
        beamSpotSrc = cms.InputTag("hltOnlineBeamSpot"),
	genVtxPositionSrc = cms.InputTag(""),
        genVtxTimeSrc = cms.InputTag(""),
        #vtxSrc = cms.InputTag("hltPhase2L3MuonPixelVertices"),
        vtxSrc = cms.InputTag("hltPhase2PixelVertices"),
        updateTrackTrajectory = cms.bool(True),
        updateTrackExtra = cms.bool(True),
        updateTrackHitPattern = cms.bool(True),
        TransientTrackBuilder = cms.string("TransientTrackBuilder"),
        MTDRecHitBuilder = cms.string("MTDRecHitBuilder"),
        Propagator = cms.string("PropagatorWithMaterialForMTD"),
        TrackTransformer = cms.PSet(
            DoPredictionsOnly = cms.bool(False),
            Fitter = cms.string("KFFitterForRefitInsideOut"),
            Smoother = cms.string("KFSmootherForRefitInsideOut"),
            Propagator = cms.string("PropagatorWithMaterialForMTD"),
            RefitDirection = cms.string("alongMomentum"),
            RefitRPCHits = cms.bool(True),
            TrackerRecHitBuilder = cms.string("WithTrackAngle"),
            MuonRecHitBuilder = cms.string("MuonRecHitBuilder"),
            MTDRecHitBuilder = cms.string("MTDRecHitBuilder"),
        ),
        estimatorMaxChi2 = cms.double(500.),
        estimatorMaxNSigma = cms.double(10.),
        btlChi2Cut = cms.double(50.),
        btlTimeChi2Cut = cms.double(10.),
        etlChi2Cut = cms.double(50.),
        etlTimeChi2Cut = cms.double(10.),
        useVertex = cms.bool(True),
        useSimVertex = cms.bool(False),
        dZCut = cms.double(0.1),
        bsTimeSpread = cms.double(0.2),
    )

    """
    process.hltPhase2GeneralTracksMTDExtendedVtxMVATrackQuality = cms.EDProducer(
        "MTDTrackQualityMVAProducer",
        tracksSrc = cms.InputTag("hltPhase2L3MuonGeneralTracks"),
        btlMatchChi2Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:btlMatchChi2"),
        btlMatchTimeChi2Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:btlMatchTimeChi2"),
        etlMatchChi2Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:etlMatchChi2"),
        etlMatchTimeChi2Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:etlMatchTimeChi2"),
        mtdTimeSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracktmtd"),
        pathLengthSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackPathLength"),
        npixBarrelSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:npixBarrel"),
        npixEndcapSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:npixEndcap"),
        qualityBDT_weights_file = cms.FileInPath("RecoMTD/TimingIDTools/data/clf4D_MTDquality_bo.xml"),
    )
    
    process.hltPhase2GeneralTracksMTDTExtendedVtxTof = cms.EDProducer(
        "TOFPIDProducer",
        tracksSrc = cms.InputTag("hltPhase2L3MuonGeneralTracks"),
        t0Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackt0"),
        tmtdSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracktmtd"),
        sigmat0Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracksigmat0"),
        sigmatmtdSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTracksigmatmtd"),
        tofkSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackTofK"),
        tofpSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackTofP"),
        sigmatofpiSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackSigmaTofPi"),
        sigmatofkSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackSigmaTofK"),
        sigmatofpSrc = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackSigmaTofP"),
        vtxsSrc = cms.InputTag("hltPhase2L3MuonPixelVertices"),
        trackMTDTimeQualityVMapTag = cms.InputTag("hltPhase2GeneralTracksMTDExtendedVtxMVATrackQuality:mtdQualMVA"),
        vtxMaxSigmaT = cms.double(0.025),
        maxDz = cms.double(0.1),
        maxDtSignificance = cms.double(5.0),
        minProbHeavy = cms.double(0.75),
        fixedT0Error = cms.double(0.),
        probPion = cms.double(1.),
        probKaon = cms.double(1.),
        probProton = cms.double(1.),
        minTrackTimeQuality = cms.double(0.8),
        MVASel = cms.bool( False ),
        vertexReassignment = cms.bool( True )
    )
    """
    
    if doMuon:
        
        process.hltPhase2L3OIMuCtfWithMaterialTracks.TrajectoryInEvent = True
        process.hltPhase2L3OIMuonTrackSelectionHighPurity.copyTrajectories = True

        process.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks.TrajectoryInEvent = True
        process.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity.copyTrajectories = True

        process.hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks.TrajectoryInEvent = True
        process.hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity.copyTrajectories = True

        # replace TrackListMerger by TrackCollectionMerger

        #process.hltIter2Phase2L3FromL1TkMuonMerged = TrackCollectionMerger.clone(    # For 2-iteration reconstruction
        #    trackProducers   = ["hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks", "hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks"],
        #    inputClassifiers = ["hltIter0Phase2L3FromL1TkMuonTrackCutClassifier", "hltIter2Phase2L3FromL1TkMuonTrackCutClassifier"],
        #    foundHitBonus  = 100.0,
        #    lostHitPenalty =   1.0,
        #    minQuality = cms.string('highPurity'),
        #    copyExtras = cms.untracked.bool(True),
        #    copyTrajectories = cms.untracked.bool(True),
        #)

        process.hltPhase2L3MuonFilter = cms.EDProducer("Phase2HLTMuonSelectorForL3",
                                                       l1TkMuons = cms.InputTag("l1tTkMuonsGmt"),
                                                       l2MuonsUpdVtx = cms.InputTag("hltL2MuonsFromL1TkMuon:UpdatedAtVtx"),
                                                       #l3Tracks = cms.InputTag("hltIter2Phase2L3FromL1TkMuonMerged"),
                                                       l3Tracks = cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity"), # Alpaka 
                                                       IOFirst = cms.bool(True),
                                                       matchingDr = cms.double(0.02),
                                                       applyL3Filters = cms.bool(True),
                                                       MinNhits = cms.int32(1),
                                                       MaxNormalizedChi2 = cms.double(5.0),
                                                       MinNhitsMuons = cms.int32(0),
                                                       MinNhitsPixel = cms.int32(1),
                                                       MinNhitsTracker = cms.int32(6),
                                                       MaxPtDifference = cms.double(999.0),
                                                       copyTrajectories = cms.bool(True)
        )

        
        process.hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx = cms.EDProducer(
            "TrackExtenderWithMTD",
            tracksSrc = cms.InputTag("hltPhase2L3MuonMerged"),
            trjtrkAssSrc = cms.InputTag("hltPhase2L3MuonMerged"),
            hitsSrc = cms.InputTag("hltPhase2L3MuonMTDRecHits"),
            beamSpotSrc = cms.InputTag("hltOnlineBeamSpot"),
            genVtxPositionSrc = cms.InputTag(""),
            genVtxTimeSrc = cms.InputTag(""),
            #vtxSrc = cms.InputTag("hltPhase2L3MuonPixelVertices"),
            vtxSrc = cms.InputTag("hltPhase2PixelVertices"),
            updateTrackTrajectory = cms.bool(True),
            updateTrackExtra = cms.bool(True),
            updateTrackHitPattern = cms.bool(True),
            TransientTrackBuilder = cms.string("TransientTrackBuilder"),
            MTDRecHitBuilder = cms.string("MTDRecHitBuilder"),
            Propagator = cms.string("PropagatorWithMaterialForMTD"),
            TrackTransformer = cms.PSet(
                DoPredictionsOnly = cms.bool(False),
                Fitter = cms.string("KFFitterForRefitInsideOut"),
                Smoother = cms.string("KFSmootherForRefitInsideOut"),
                Propagator = cms.string("PropagatorWithMaterialForMTD"),
                RefitDirection = cms.string("alongMomentum"),
                RefitRPCHits = cms.bool(True),
                TrackerRecHitBuilder = cms.string("WithTrackAngle"),
                MuonRecHitBuilder = cms.string("MuonRecHitBuilder"),
                MTDRecHitBuilder = cms.string("MTDRecHitBuilder"),
            ),
            estimatorMaxChi2 = cms.double(500.),
            estimatorMaxNSigma = cms.double(10.),
            btlChi2Cut = cms.double(50.),
            btlTimeChi2Cut = cms.double(10.),
            etlChi2Cut = cms.double(50.),
            etlTimeChi2Cut = cms.double(10.),
            useVertex = cms.bool(True),
            useSimVertex = cms.bool(False),
            dZCut = cms.double(0.1),
            bsTimeSpread = cms.double(0.2),
            doMuon = cms.bool(True),
        )
        

    from SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociationDefault_cfi import mtdSimLayerClusterToTPAssociationDefault as _mtdSimLayerClusterToTPAssociationDefault
    mtdSimLayerClusterToTPAssociation = _mtdSimLayerClusterToTPAssociationDefault.clone()
    from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
    premix_stage2.toModify(mtdSimLayerClusterToTPAssociation, mtdSimClustersTag = "mixData:MergedMtdTruthLC")
    premix_stage2.toModify(mtdSimLayerClusterToTPAssociation, trackingParticlesTag = "mixData:MergedTrackTruth")

    from SimFastTiming.MtdAssociatorProducers.mtdRecoClusterToSimLayerClusterAssociationDefault_cfi import mtdRecoClusterToSimLayerClusterAssociationDefault as _mtdRecoClusterToSimLayerClusterAssociationDefault
    mtdRecoClusterToSimLayerClusterAssociation = _mtdRecoClusterToSimLayerClusterAssociationDefault.clone()
    from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
    premix_stage2.toModify(mtdRecoClusterToSimLayerClusterAssociation, mtdSimClustersTag = "mixData:MergedMtdTruthLC")
    
    process.hltPhase2L3MuonsPFIsodR0p4Dt = cms.EDProducer("MuonHLTPFCandidateIsolationWithMTDProducer",
                                                          recoCandidateProducer = cms.InputTag("hltPhase2L3MuonCandidates"),
                                                          pfCandidateProducer = cms.InputTag("hltParticleFlowTmp"),
                                                          #
                                                          # tracksSrc = cms.InputTag("hltGeneralTracks"),
                                                          t0Src = cms.InputTag("hltPhase2GeneralTracksMTDTExtendedVtx:generalTrackt0"),
                                                          #
                                                          candTrackt0 = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTrackt0"),
                                                          candTracksigmat0 = cms.InputTag("hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx:generalTracksigmat0"),
                                                          #
                                                          drMax = cms.double(0.4),
                                                          drVeto = cms.double(0.01),
                                                          drVetoCh = cms.double(0.0001),
                                                          minEnergy = cms.double(0.0)
    )

    process.hltL3crIsoL1TkSingleMu22PFTimeIso0p15 = cms.EDFilter("HLTMuonGenericFilter",
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
                                                             thrOverEEB = cms.vdouble(0.15),
                                                             thrOverEEE = cms.vdouble(0.15),
                                                             thrRegularEB = cms.vdouble(-1.0),
                                                             thrRegularEE = cms.vdouble(-1.0),
                                                             useEt = cms.bool(True),
                                                             varTag = cms.InputTag("hltPhase2L3MuonsPFIsodR0p4Dt")
    )

    """
    if doMuon:
        replaceWith = (
            process.hltPhase2L3MuonMtdUncalibratedRecHits+
            process.hltPhase2L3MuonMtdRecHits+
            process.hltPhase2L3MuonMTDClusters+
            process.hltPhase2L3MuonMTDRecHits+
            process.hltPhase2GeneralTracksMTDTExtendedVtx+
            # Muon
            process.hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx
        )
    else:
        replaceWith = (
            process.hltPhase2L3MuonMtdUncalibratedRecHits+
            process.hltPhase2L3MuonMtdRecHits+
            process.hltPhase2L3MuonMTDClusters+
            process.hltPhase2L3MuonMTDRecHits+
            process.hltPhase2GeneralTracksMTDTExtendedVtx
            #process.hltPhase2GeneralTracksMTDExtendedVtxMVATrackQuality+
            #process.hltPhase2GeneralTracksMTDTExtendedVtxTof
        )

    process.HLTPhase2L3MuonGeneralTracksSequence += replaceWith
    """
    if doMuon:
        process.HLT_IsoMu24_FromL1TkMuon = cms.Path(
            process.HLTBeginSequence
            + process.hltSingleTkMuon22L1TkMuonFilter
            + process.HLTRawToDigiSequence
            + process.HLTItLocalRecoSequence
            + process.HLTOtLocalRecoSequence
            + process.hltPhase2PixelFitterByHelixProjections
            + process.hltPhase2PixelTrackFilterByKinematics
            + process.HLTTrackingSequence
            + process.HLTMuonsSequence
            + process.hltPhase2L3MuonMtdUncalibratedRecHits
            + process.hltPhase2L3MuonMtdRecHits
            + process.hltPhase2L3MuonMTDClusters
            + process.hltPhase2L3MuonMTDRecHits
            + process.hltPhase2GeneralTracksMTDTExtendedVtx
            + process.hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx
            + process.hltL3fL1TkSingleMu22L3Filtered24Q
            + process.HLTLocalrecoSequence
            + process.HLTTICLLocalRecoSequence
            + process.HLTParticleFlowSequence
            + process.hltPhase2L3MuonsPFIsodR0p4Dt
            + process.hltL3crIsoL1TkSingleMu22PFTimeIso0p15
            + process.HLTEndSequence
        )
    else:
        process.HLT_IsoMu24_FromL1TkMuon = cms.Path(
            process.HLTBeginSequence
            + process.hltSingleTkMuon22L1TkMuonFilter
            + process.HLTRawToDigiSequence
            + process.HLTItLocalRecoSequence
            + process.HLTOtLocalRecoSequence
            + process.hltPhase2PixelFitterByHelixProjections
            + process.hltPhase2PixelTrackFilterByKinematics
            + process.HLTTrackingSequence
            + process.HLTMuonsSequence
            + process.hltPhase2L3MuonMtdUncalibratedRecHits
            + process.hltPhase2L3MuonMtdRecHits
            + process.hltPhase2L3MuonMTDClusters
            + process.hltPhase2L3MuonMTDRecHits
            + process.hltPhase2GeneralTracksMTDTExtendedVtx
            + process.hltL3fL1TkSingleMu22L3Filtered24Q
            + process.HLTLocalrecoSequence
            + process.HLTTICLLocalRecoSequence
            + process.HLTParticleFlowSequence
            + process.hltPhase2L3MuonsPFIsodR0p4Dt
            + process.hltL3crIsoL1TkSingleMu22PFTimeIso0p15
            + process.HLTEndSequence
        )

    #process.HLT_IsoMu24_FromL1TkMuon = cms.Path(
    #    process.HLTBeginSequence
    #    + process.hltSingleTkMuon22L1TkMuonFilter
    #    + process.HLTRawToDigiSequence
    #    + process.HLTItLocalRecoSequence
    #    + process.HLTOtLocalRecoSequence
    #    + process.hltPhase2PixelFitterByHelixProjections
    #    + process.hltPhase2PixelTrackFilterByKinematics
    #    + process.HLTTrackingSequence
    #    + process.HLTMuonsSequence
    #    + process.hltL3fL1TkSingleMu22L3Filtered24Q
    #    #
    #    + process.HLTLocalrecoSequence
    #    + process.HLTTICLLocalRecoSequence
    #    + process.HLTParticleFlowSequence
    #    + process.hltPhase2L3MuonMtdUncalibratedRecHits
    #    + process.hltPhase2L3MuonMtdRecHits
    #    + process.hltPhase2L3MuonMTDClusters
    #    + process.hltPhase2L3MuonMTDRecHits
    #    + process.hltPhase2GeneralTracksMTDTExtendedVtx
    #    + process.hltPhase2L3MuonGeneralMuonTrackMTDTExtendedVtx
    #    + process.hltPhase2L3MuonsPFIsodR0p4Dt
    #    + process.hltL3crIsoL1TkSingleMu22PFTimeIso0p15
    #    #
    #    + process.HLTEndSequence
    #)    
    
    return process
      
    
def customizerFuncForMuonHLTNtupler(process, newProcessName = "MYHLT", doDYSkim = False):
    if hasattr(process, "DQMOutput"):
        del process.DQMOutput

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

    # -- track - TP associations
    import SimMuon.MCTruth.MuonTrackProducer_cfi
    process.hltPhase2L3MuonsNoIDTracks = SimMuon.MCTruth.MuonTrackProducer_cfi.muonTrackProducer.clone()
    process.hltPhase2L3MuonsNoIDTracks.muonsTag                      = cms.InputTag("hltPhase2L3MuonsNoID")
    process.hltPhase2L3MuonsNoIDTracks.selectionTags                 = ('All',)
    process.hltPhase2L3MuonsNoIDTracks.trackType                     = "recomuonTrack"
    process.hltPhase2L3MuonsNoIDTracks.ignoreMissingMuonCollection   = True
    process.hltPhase2L3MuonsNoIDTracks.inputCSCSegmentCollection     = cms.InputTag("hltCscSegments")
    process.hltPhase2L3MuonsNoIDTracks.inputDTRecSegment4DCollection = cms.InputTag("hltDt4DSegments")

    process.hltPhase2L3MuonsTracks = SimMuon.MCTruth.MuonTrackProducer_cfi.muonTrackProducer.clone()
    process.hltPhase2L3MuonsTracks.muonsTag                          = cms.InputTag("hltPhase2L3Muons")
    process.hltPhase2L3MuonsTracks.selectionTags                     = ('All',)
    process.hltPhase2L3MuonsTracks.trackType                         = "recomuonTrack"
    process.hltPhase2L3MuonsTracks.ignoreMissingMuonCollection       = True
    process.hltPhase2L3MuonsTracks.inputCSCSegmentCollection         = cms.InputTag("hltCscSegments")
    process.hltPhase2L3MuonsTracks.inputDTRecSegment4DCollection     = cms.InputTag("hltDt4DSegments")

    from SimMuon.MCTruth.MuonAssociatorByHits_cfi import muonAssociatorByHits as _muonAssociatorByHits
    hltMuonAssociatorByHits = _muonAssociatorByHits.clone()
    hltMuonAssociatorByHits.PurityCut_track              = 0.75
    hltMuonAssociatorByHits.PurityCut_muon               = 0.75
    hltMuonAssociatorByHits.DTrechitTag                  = 'hltDt1DRecHits'
    hltMuonAssociatorByHits.ignoreMissingTrackCollection = True
    hltMuonAssociatorByHits.UseTracker                   = True
    hltMuonAssociatorByHits.UseMuon                      = True

    process.AhltPhase2L3OIMuonTrackSelectionHighPurity            = hltMuonAssociatorByHits.clone( tracksTag = 'hltPhase2L3OIMuonTrackSelectionHighPurity' )
    process.AhltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity = hltMuonAssociatorByHits.clone( tracksTag = 'hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity' )
    process.AhltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity = hltMuonAssociatorByHits.clone( tracksTag = 'hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity' )
    process.AhltIter2Phase2L3FromL1TkMuonMerged                   = hltMuonAssociatorByHits.clone( tracksTag = 'hltIter2Phase2L3FromL1TkMuonMerged' )
    process.AhltPhase2L3MuonsNoID                                 = hltMuonAssociatorByHits.clone( tracksTag = 'hltPhase2L3MuonsNoIDTracks' )
    process.AhltPhase2L3Muons                                     = hltMuonAssociatorByHits.clone( tracksTag = 'hltPhase2L3MuonsTracks' )

    trackNames = [
        'hltPhase2L3OI',
        'hltIter0Phase2L3FromL1TkMuon',
        'hltIter2Phase2L3FromL1TkMuon',
        'hltPhase2L3IOFromL1',
        'hltPhase2L3MuonsNoID',
        'hltPhase2L3Muons'
    ]

    trackLabels = [
        'hltPhase2L3OIMuonTrackSelectionHighPurity',
        'hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity',
        'hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity',
        'hltIter2Phase2L3FromL1TkMuonMerged',
        'hltPhase2L3MuonsNoIDTracks',
        'hltPhase2L3MuonsTracks'
    ]

    assoLabels = [
        'AhltPhase2L3OIMuonTrackSelectionHighPurity',
        'AhltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity',
        'AhltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity',
        'AhltIter2Phase2L3FromL1TkMuonMerged',
        'AhltPhase2L3MuonsNoID',
        'AhltPhase2L3Muons'
    ]

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

    # -- Isolations
    trkIsoTags = []
    trkIsoLabels = []
    pfIsoTags = []
    pfIsoLabels = []

    from MuonHLTTool.MuonHLTNtupler.ntupler_cfi import ntuplerBase
    process.ntupler = ntuplerBase.clone()

    process.ntupler.trackCollectionNames  = cms.untracked.vstring(   trackNames )
    process.ntupler.trackCollectionLabels = cms.untracked.VInputTag( trackLabels )
    process.ntupler.associationLabels     = cms.untracked.VInputTag( assoLabels )

    process.ntupler.trkIsoTags   = cms.untracked.vstring(   trkIsoTags )
    process.ntupler.trkIsoLabels = cms.untracked.VInputTag( trkIsoLabels )
    process.ntupler.pfIsoTags    = cms.untracked.vstring(   pfIsoTags )
    process.ntupler.pfIsoLabels  = cms.untracked.VInputTag( pfIsoLabels )

    # -- set to the new process name
    process.ntupler.myTriggerResults = cms.untracked.InputTag("TriggerResults",          "",     newProcessName)
    process.ntupler.myTriggerEvent   = cms.untracked.InputTag("hltTriggerSummaryAOD",    "",     newProcessName)
    process.ntupler.lumiScaler       = cms.untracked.InputTag("hltScalersRawToDigi",     "",     newProcessName)

    # process.ntupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", newProcessName)
    # process.ntupler.L1Muon           = cms.untracked.InputTag("gmtStage2Digis",        "Muon", newProcessName) 
    # process.ntupler.L1Muon           = cms.untracked.InputTag("hltGtStage2Digis",        "Muon", "HLT") #for phaseII w/o emulation
    process.ntupler.L1Muon                        = cms.untracked.InputTag("simGmtStage2Digis",                  "", newProcessName)  # Phase II sim emul
    process.ntupler.L2Muon                        = cms.untracked.InputTag("hltL2MuonFromL1TkMuonCandidates",    "", newProcessName)
    process.ntupler.L3Muon                        = cms.untracked.InputTag("hltPhase2L3MuonCandidates",          "", newProcessName)
    process.ntupler.TkMuon                        = cms.untracked.InputTag("hltHighPtTkMuonCands",               "", newProcessName)

    process.ntupler.iterL3OI                      = cms.untracked.InputTag("hltL3MuonsPhase2L3OI",               "", newProcessName)
    process.ntupler.iterL3IOFromL1                = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonMerged", "", newProcessName)
    process.ntupler.iterL3MuonNoID                = cms.untracked.InputTag("hltPhase2L3MuonsNoID",               "", newProcessName)
    process.ntupler.iterL3Muon                    = cms.untracked.InputTag("hltPhase2L3Muons",                   "", newProcessName)

    process.ntupler.hltIterL3FromL1MuonTrimmedPixelVertices           = cms.untracked.InputTag("hltPhase2L3FromL1TkMuonTrimmedPixelVertices",           "", newProcessName)

    process.ntupler.doMVA  = cms.bool(False)
    process.ntupler.doSeed = cms.bool(False)

    process.ntupler.hltIterL3OISeedsFromL2Muons                       = cms.untracked.InputTag("hltPhase2L3OISeedsFromL2Muons",                         "", newProcessName)
    process.ntupler.hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks = cms.untracked.InputTag("hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks", "", newProcessName)
    process.ntupler.hltIter2IterL3FromL1MuonPixelSeeds                = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeeds",                "", newProcessName)

    process.ntupler.hltIterL3OIMuonTrack                              = cms.untracked.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity",             "", newProcessName)
    process.ntupler.hltIter0IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity",  "", newProcessName)
    process.ntupler.hltIter2IterL3FromL1MuonTrack                     = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity",  "", newProcessName)

    process.ntupler.associator = cms.untracked.InputTag("hltTrackAssociatorByHits")
    process.ntupler.trackingParticle = cms.untracked.InputTag("mix","MergedTrackTruth")

    process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0                = cms.untracked.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Phase2_Iter2FromL1_barrel_v0.xml")
    process.ntupler.mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0                = cms.untracked.FileInPath("RecoMuon/TrackerSeedGenerator/data/xgb_Phase2_Iter2FromL1_endcap_v0.xml")
    process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B                = cms.untracked.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleMean)
    process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B                 = cms.untracked.vdouble(PU200_Barrel_NThltIter2FromL1_ScaleStd)
    process.ntupler.mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E                = cms.untracked.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleMean)
    process.ntupler.mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E                 = cms.untracked.vdouble(PU200_Endcap_NThltIter2FromL1_ScaleStd)

    process.TFileService = cms.Service("TFileService",
      fileName = cms.string("test_ntuple_seedntuple.root"),
      closeFileFast = cms.untracked.bool(False),
    )

    process.ntupler.DebugMode = cms.bool(False)
    process.ntupler.SaveAllTracks = cms.bool(True)
    # process.ntupler.SaveStubs = cms.bool(False)
    process.ntupler.L1TrackInputTag = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks") # TTTrack input
    # process.ntupler.MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks")  ## MCTruth input
    # process.ntupler.L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted")
    process.ntupler.TkMuonToken = cms.InputTag("L1TkMuons")
    process.ntupler.l1PrimaryVertex = cms.InputTag("l1tVertexFinderEmulator", "L1VerticesEmulation")

    # if doDYSkim:
    #     from MuonHLTTool.MuonHLTNtupler.DYmuSkimmer import DYmuSkimmer
    #     process.Skimmer = DYmuSkimmer.clone()
    #     process.mypath = cms.Path(process.Skimmer*process.hltTPClusterProducer*process.hltTrackAssociatorByHits*process.trackAssoSeq*process.ntupler)

    # else:
    #     process.mypath = cms.Path(process.hltTPClusterProducer*process.hltTrackAssociatorByHits*process.trackAssoSeq*process.ntupler)

    process.mypath    = cms.Path(process.hltTPClusterProducer*process.hltTrackAssociatorByHits*process.trackAssoSeq)
    process.myendpath = cms.EndPath(process.ntupler)

    return process
