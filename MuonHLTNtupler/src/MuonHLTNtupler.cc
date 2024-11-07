// -- ntuple maker for Muon HLT study
// -- author: Kyeongpil Lee (Seoul National University, kplee@cern.ch)

#include "MuonHLTTool/MuonHLTNtupler/interface/MuonHLTNtupler.h"

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h" // New
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "CommonTools/Utils/interface/associationMapFilterValues.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractorFactory.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

#include "RecoMuon/MuonIsolation/interface/Range.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace muonisolation;
using reco::isodeposit::Direction;

MuonHLTNtupler::MuonHLTNtupler(const edm::ParameterSet& iConfig):
  //theMuonCollectionLabel(iConfig.getParameter<InputTag>("inputMuonCollection")),  
  
  useHF_                (iConfig.getParameter<bool>("useHF")),
  drMax_HCAL_           (iConfig.getParameter<double>("drMax_HCAL")),
  drVetoBarrel_HCAL_    (iConfig.getParameter<double>("drVetoBarrel_HCAL")),
  drVetoEndcap_HCAL_    (iConfig.getParameter<double>("drVetoEndcap_HCAL")),
  etaStripBarrel_HCAL_  (iConfig.getParameter<double>("etaStripBarrel_HCAL")),
  etaStripEndcap_HCAL_  (iConfig.getParameter<double>("etaStripEndcap_HCAL")),
  energyBarrel_HCAL_    (iConfig.getParameter<double>("energyBarrel_HCAL")),
  energyEndcap_HCAL_    (iConfig.getParameter<double>("energyEndcap_HCAL")),
  doRhoCorrection_HCAL_ (iConfig.getParameter<bool>("doRhoCorrection_HCAL")),
  rhoMax_HCAL_          (iConfig.getParameter<double>("rhoMax_HCAL")),
  rhoScale_HCAL_        (iConfig.getParameter<double>("rhoScale_HCAL")),
  effectiveAreas_HCAL_  (iConfig.getParameter<std::vector<double>>("effectiveAreas_HCAL")),
  absEtaLowEdges_HCAL_  (iConfig.getParameter<std::vector<double>>("absEtaLowEdges_HCAL")),

  drMax_ECAL_           (iConfig.getParameter<double>("drMax_ECAL")),
  drVetoBarrel_ECAL_    (iConfig.getParameter<double>("drVetoBarrel_ECAL")),
  drVetoEndcap_ECAL_    (iConfig.getParameter<double>("drVetoEndcap_ECAL")),
  etaStripBarrel_ECAL_  (iConfig.getParameter<double>("etaStripBarrel_ECAL")),
  etaStripEndcap_ECAL_  (iConfig.getParameter<double>("etaStripEndcap_ECAL")),
  energyBarrel_ECAL_    (iConfig.getParameter<double>("energyBarrel_ECAL")),
  energyEndcap_ECAL_    (iConfig.getParameter<double>("energyEndcap_ECAL")),
  doRhoCorrection_ECAL_ (iConfig.getParameter<bool>("doRhoCorrection_ECAL")),
  rhoMax_ECAL_          (iConfig.getParameter<double>("rhoMax_ECAL")),
  rhoScale_ECAL_        (iConfig.getParameter<double>("rhoScale_ECAL")),
  effectiveAreas_ECAL_  (iConfig.getParameter<std::vector<double>>("effectiveAreas_ECAL")),
  absEtaLowEdges_ECAL_  (iConfig.getParameter<std::vector<double>>("absEtaLowEdges_ECAL")),
  
  useEt_                (iConfig.getParameter<bool>("useEt")),

  theDepositLabel  (iConfig.getUntrackedParameter<string>("DepositLabel")),
  theDiff_r        (iConfig.getParameter<double>("Diff_r")),
  theDiff_z        (iConfig.getParameter<double>("Diff_z")),
  theDR_Max        (iConfig.getParameter<double>("DR_Max")),
  theDR_Veto       (iConfig.getParameter<double>("DR_Veto")),
  theNHits_Min     (iConfig.getParameter<unsigned int>("NHits_Min")),
  theChi2Ndof_Max  (iConfig.getParameter<double>("Chi2Ndof_Max")),
  theChi2Prob_Min  (iConfig.getParameter<double>("Chi2Prob_Min")),
  thePt_Min        (iConfig.getParameter<double>("Pt_Min")),
  
  propSetup_(iConfig, consumesCollector()),
  
  t_beamSpot_          ( consumes< reco::BeamSpot >                         (iConfig.getUntrackedParameter<edm::InputTag>("beamSpot"     )) ),
  t_offlineMuon_       ( consumes< std::vector<reco::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("offlineMuon"       )) ),
  t_offlineVertex_     ( consumes< reco::VertexCollection >                 (iConfig.getUntrackedParameter<edm::InputTag>("offlineVertex"     )) ),
  t_rho_ECAL_          ( consumes< double >                                 (iConfig.getUntrackedParameter<edm::InputTag>("rho_ECAL"          )) ),
  t_rho_HCAL_          ( consumes< double >                                 (iConfig.getUntrackedParameter<edm::InputTag>("rho_HCAL"          )) ),
  
  t_lumiScaler_        ( consumes< LumiScalersCollection >                  (iConfig.getUntrackedParameter<edm::InputTag>("lumiScaler"        )) ),
  t_offlineLumiScaler_ ( consumes< LumiScalersCollection >                  (iConfig.getUntrackedParameter<edm::InputTag>("offlineLumiScaler" )) ),
  t_PUSummaryInfo_     ( consumes< std::vector<PileupSummaryInfo> >         (iConfig.getUntrackedParameter<edm::InputTag>("PUSummaryInfo"     )) ),
  t_genEventInfo_      ( consumes< GenEventInfoProduct >                    (iConfig.getUntrackedParameter<edm::InputTag>("genEventInfo"      )) ),
  t_genParticle_       ( consumes< edm::View<reco::GenParticle> >           (iConfig.getUntrackedParameter<edm::InputTag>("genParticle"       )) ),

  t_triggerResults_    ( consumes< edm::TriggerResults >                    (iConfig.getUntrackedParameter<edm::InputTag>("triggerResults"    )) ),
  t_triggerEvent_      ( consumes< trigger::TriggerEvent >                  (iConfig.getUntrackedParameter<edm::InputTag>("triggerEvent"      )) ),
  t_myTriggerResults_  ( consumes< edm::TriggerResults >                    (iConfig.getUntrackedParameter<edm::InputTag>("myTriggerResults"  )) ),
  t_myTriggerEvent_    ( consumes< trigger::TriggerEvent >                  (iConfig.getUntrackedParameter<edm::InputTag>("myTriggerEvent"    )) ),
  
  pfClusterProducerHCAL_  ( consumes<reco::PFClusterCollection>             (iConfig.getUntrackedParameter<edm::InputTag>("pfClusterProducerHCAL"))  ),
  rhoProducer_HCAL_       ( consumes<double>                                (iConfig.getUntrackedParameter<edm::InputTag>("rhoProducer_HCAL"))       ),
  pfClusterProducerHFEM_  ( consumes<reco::PFClusterCollection>             (iConfig.getUntrackedParameter<edm::InputTag>("pfClusterProducerHFEM"))  ),
  pfClusterProducerHFHAD_ ( consumes<reco::PFClusterCollection>             (iConfig.getUntrackedParameter<edm::InputTag>("pfClusterProducerHFHAD")) ),

  pfClusterProducer_      ( consumes<reco::PFClusterCollection>             (iConfig.getUntrackedParameter<edm::InputTag>("pfClusterProducer"))      ),
  //rhoProducer_ECAL_       ( consumes<double>                                (iConfig.getUntrackedParameter<edm::InputTag>("rhoProducer_ECAL"))       ),

  theTrackCollectionToken ( consumes<TrackCollection>                       (iConfig.getUntrackedParameter<edm::InputTag>("inputTrackCollection"))   ),
  //theMuonCollectionToken  ( consumes<reco::RecoChargedCandidateCollection>  (iConfig.getUntrackedParameter<edm::InputTag>("inputMuonCollection"))    )
  //theMuonCollectionToken  ( consumes<reco::RecoChargedCandidateCollection>  (iConfig.getParameter<edm::InputTag>("inputMuonCollection"))    )
  theMuonCollectionToken (consumes<trigger::TriggerFilterObjectWithRefs>         (iConfig.getParameter<edm::InputTag>("inputMuonCollection")))
{
}  

void MuonHLTNtupler::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  Init();

  // -- basic info.
  isRealData_ = iEvent.isRealData();

  runNum_       = iEvent.id().run();
  lumiBlockNum_ = iEvent.id().luminosityBlock();
  eventNum_     = iEvent.id().event();

  edm::Handle<reco::BeamSpot> h_beamSpot;
  iEvent.getByToken(t_beamSpot_, h_beamSpot);
  bs = h_beamSpot.isValid() ? &*h_beamSpot : 0;
  bs_x0_ = bs->x0();
  bs_y0_ = bs->y0();
  bs_z0_ = bs->z0();
  bs_sigmaZ_ = bs->sigmaZ();
  bs_dxdz_ = bs->dxdz();
  bs_dydz_ = bs->dydz();
  bs_x0Error_ = bs->x0Error();
  bs_y0Error_ = bs->y0Error();
  bs_z0Error_ = bs->z0Error();
  bs_sigmaZ0Error_ = bs->sigmaZ0Error();
  bs_dxdzError_ = bs->dxdzError();
  bs_dydzError_ = bs->dydzError();


  // -- vertex
  edm::Handle<reco::VertexCollection> h_offlineVertex;
  if( iEvent.getByToken(t_offlineVertex_, h_offlineVertex) )
  {
    int nGoodVtx = 0;
    for(reco::VertexCollection::const_iterator it = h_offlineVertex->begin(); it != h_offlineVertex->end(); ++it){
      if( it->isValid() ){
	vertex_x_[nGoodVtx] = it->x();
	vertex_y_[nGoodVtx] = it->y();
	vertex_z_[nGoodVtx] = it->z();
        vertex_t_[nGoodVtx] = it->t();

	nGoodVtx++;
      }
    }
    nVertex_ = nGoodVtx;
  }
  
  if( isRealData_ )
  {
    bunchID_ = iEvent.bunchCrossing();

    // -- lumi scaler @ HLT
    edm::Handle<LumiScalersCollection> h_lumiScaler;
    if( iEvent.getByToken(t_lumiScaler_, h_lumiScaler) && h_lumiScaler->begin() != h_lumiScaler->end() )
    {
      instLumi_  = h_lumiScaler->begin()->instantLumi();
      dataPU_    = h_lumiScaler->begin()->pileup();
      dataPURMS_ = h_lumiScaler->begin()->pileupRMS();
      bunchLumi_ = h_lumiScaler->begin()->bunchLumi();
    }

    // -- lumi scaler @ offline
    edm::Handle<LumiScalersCollection> h_offlineLumiScaler;
    if( iEvent.getByToken(t_offlineLumiScaler_, h_offlineLumiScaler) && h_offlineLumiScaler->begin() != h_offlineLumiScaler->end() )
    {
      offlineInstLumi_  = h_offlineLumiScaler->begin()->instantLumi();
      offlineDataPU_    = h_offlineLumiScaler->begin()->pileup();
      offlineDataPURMS_ = h_offlineLumiScaler->begin()->pileupRMS();
      offlineBunchLumi_ = h_offlineLumiScaler->begin()->bunchLumi();
    }
  }

  // -- True PU info: only for MC -- //
  if( !isRealData_ )
  {
    edm::Handle<std::vector< PileupSummaryInfo > > h_PUSummaryInfo;

    if( iEvent.getByToken(t_PUSummaryInfo_,h_PUSummaryInfo) )
    {
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = h_PUSummaryInfo->begin(); PVI != h_PUSummaryInfo->end(); ++PVI)
      {
        if(PVI->getBunchCrossing()==0)
        {
          truePU_ = PVI->getTrueNumInteractions();
          continue;
        }
      } // -- end of PU iteration -- //
    } // -- end of if ( token exists )
  } // -- end of isMC -- //

  // -- fill each object
  Fill_Muon(iEvent, iSetup);
  Fill_Muon2(iEvent, iSetup);
  Fill_HLT(iEvent, 0);
  Fill_Track(iEvent, iSetup);
  Fill_ECAL(iEvent, iSetup);
  Fill_HCAL(iEvent, iSetup);
  //if( doSeed )  Fill_Seed(iEvent, iSetup);
  if( !isRealData_ ) {
    Fill_GenParticle(iEvent);
  }

  ntuple_->Fill();
}

void MuonHLTNtupler::beginJob()
{
  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TTree>("ntuple","ntuple");

  Make_Branch();
}

void MuonHLTNtupler::Init()
{
  isRealData_ = false;

  runNum_       = -999;
  lumiBlockNum_ = -999;
  eventNum_     = 0;

  bunchID_ = -999;

  bs_x0_ = -999;
  bs_y0_ = -999;
  bs_z0_ = -999;
  bs_sigmaZ_ = -999;
  bs_dxdz_ = -999;
  bs_dydz_ = -999;
  bs_x0Error_ = -999;
  bs_y0Error_ = -999;
  bs_z0Error_ = -999;
  bs_sigmaZ0Error_ = -999;
  bs_dxdzError_ = -999;
  bs_dydzError_ = -999;

  nVertex_ = -999;
  for( int i=0; i<arrSize_; i++)
  {
    vertex_x_[i] = -999;
    vertex_y_[i] = -999;
    vertex_z_[i] = -999;
    vertex_t_[i] = -999;
  }
  
  instLumi_  = -999;
  dataPU_    = -999;
  dataPURMS_ = -999;
  bunchLumi_ = -999;

  offlineInstLumi_  = -999;
  offlineDataPU_    = -999;
  offlineDataPURMS_ = -999;
  offlineBunchLumi_ = -999;

  truePU_ = -999;

  genEventWeight_ = -999;

  nGenParticle_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    genParticle_ID_[i] = -999;
    genParticle_status_[i] = -999;
    genParticle_mother_[i] = -999;

    genParticle_pt_[i]     = -999;
    genParticle_eta_[i]    = -999;
    genParticle_phi_[i]    = -999;
    genParticle_px_[i]     = -999;
    genParticle_py_[i]     = -999;
    genParticle_pz_[i]     = -999;
    genParticle_energy_[i] = -999;
    genParticle_charge_[i] = -999;

    genParticle_isPrompt_[i] = 0;
    genParticle_isPromptFinalState_[i] = 0;
    genParticle_isTauDecayProduct_[i] = 0;
    genParticle_isPromptTauDecayProduct_[i] = 0;
    genParticle_isDirectPromptTauDecayProductFinalState_[i] = 0;
    genParticle_isHardProcess_[i] = 0;
    genParticle_isLastCopy_[i] = 0;
    genParticle_isLastCopyBeforeFSR_[i] = 0;
    genParticle_isPromptDecayed_[i] = 0;
    genParticle_isDecayedLeptonHadron_[i] = 0;
    genParticle_fromHardProcessBeforeFSR_[i] = 0;
    genParticle_fromHardProcessDecayed_[i] = 0;
    genParticle_fromHardProcessFinalState_[i] = 0;
    genParticle_isMostlyLikePythia6Status3_[i] = 0;

  }

  nMuon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    muon_pt_[i] = -999;
    muon_eta_[i] = -999;
    muon_phi_[i] = -999;
    muon_px_[i] = -999;
    muon_py_[i] = -999;
    muon_pz_[i] = -999;
    muon_dB_[i] = -999;
    muon_charge_[i] = -999;
    muon_isGLB_[i] = 0;
    muon_isSTA_[i] = 0;
    muon_isTRK_[i] = 0;
    muon_isPF_[i] = 0;
    muon_isTight_[i] = 0;
    muon_isMedium_[i] = 0;
    muon_isLoose_[i] = 0;
    muon_isHighPt_[i] = 0;
    muon_isHighPtNew_[i] = 0;
    muon_isSoft_[i] = 0;

    muon_iso03_sumPt_[i] = -999;
    muon_iso03_hadEt_[i] = -999;
    muon_iso03_emEt_[i] = -999;

    muon_PFIso03_charged_[i] = -999;
    muon_PFIso03_neutral_[i] = -999;
    muon_PFIso03_photon_[i] = -999;
    muon_PFIso03_sumPU_[i] = -999;

    muon_PFIso04_charged_[i] = -999;
    muon_PFIso04_neutral_[i] = -999;
    muon_PFIso04_photon_[i] = -999;
    muon_PFIso04_sumPU_[i] = -999;

    muon_PFCluster03_ECAL_[i] = -999;
    muon_PFCluster03_HCAL_[i] = -999;

    muon_PFCluster04_ECAL_[i] = -999;
    muon_PFCluster04_HCAL_[i] = -999;

    muon_inner_trkChi2_[i] = -999;
    muon_inner_validFraction_[i] = -999;
    muon_inner_trackerLayers_[i] = -999;
    muon_inner_trackerHits_[i] = -999;
    muon_inner_lostTrackerHits_[i] = -999;
    muon_inner_lostTrackerHitsIn_[i] = -999;
    muon_inner_lostTrackerHitsOut_[i] = -999;
    muon_inner_lostPixelHits_[i] = -999;
    muon_inner_lostPixelBarrelHits_[i] = -999;
    muon_inner_lostPixelEndcapHits_[i] = -999;
    muon_inner_lostStripHits_[i] = -999;
    muon_inner_lostStripTIBHits_[i] = -999;
    muon_inner_lostStripTIDHits_[i] = -999;
    muon_inner_lostStripTOBHits_[i] = -999;
    muon_inner_lostStripTECHits_[i] = -999;
    muon_inner_pixelLayers_[i] = -999;
    muon_inner_pixelHits_[i] = -999;
    muon_global_muonHits_[i] = -999;
    muon_global_trkChi2_[i] = -999;
    muon_global_trackerLayers_[i] = -999;
    muon_global_trackerHits_[i] = -999;
    muon_momentumChi2_[i] = -999;
    muon_positionChi2_[i] = -999;
    muon_glbKink_[i] = -999;
    muon_glbTrackProbability_[i] = -999;
    muon_globalDeltaEtaPhi_[i] = -999;
    muon_localDistance_[i] = -999;
    muon_staRelChi2_[i] = -999;
    muon_tightMatch_[i] = -999;
    muon_trkKink_[i] = -999;
    muon_trkRelChi2_[i] = -999;
    muon_segmentCompatibility_[i] = -999;

    muon_pt_tuneP_[i] = -999;
    muon_ptError_tuneP_[i] = -999;

    muon_dxyVTX_best_[i] = -999;
    muon_dzVTX_best_[i] = -999;

    muon_nMatchedStation_[i] = -999;
    muon_nMatchedRPCLayer_[i] = -999;
    muon_stationMask_[i] = -999;

    muon_dxy_bs_[i] = -999;
    muon_dxyError_bs_[i] = -999;
    muon_dz_bs_[i] = -999;
    muon_dzError_[i] = -999;
    muon_IPSig_[i] = -999;

  }

  
  nMuonCand_ = 0;
  for( int i=0; i<arrSize_; i++)
    {
      muonCand_pt_[i] = -999;
      muonCand_eta_[i] = -999;
      muonCand_phi_[i] = -999;
      muonCand_px_[i] = -999;
      muonCand_py_[i] = -999;
      muonCand_pz_[i] = -999;
      muonCand_vx_[i] = -999;
      muonCand_vy_[i] = -999;
      muonCand_vz_[i] = -999;
      muonCand_dB_[i] = -999;
      muonCand_charge_[i] = -999;
      muonCand_isGLB_[i] = 0;
      muonCand_isSTA_[i] = 0;
      muonCand_isTRK_[i] = 0;
      muonCand_isMuon_[i] = 0;

      muonCand_inner_trkChi2_[i] = -999;
      muonCand_inner_validFraction_[i] = -999;
      muonCand_inner_trackerLayers_[i] = -999;
      muonCand_inner_trackerHits_[i] = -999;
      muonCand_inner_lostTrackerHits_[i] = -999;
      muonCand_inner_lostTrackerHitsIn_[i] = -999;
      muonCand_inner_lostTrackerHitsOut_[i] = -999;
      muonCand_inner_lostPixelHits_[i] = -999;
      muonCand_inner_lostPixelBarrelHits_[i] = -999;
      muonCand_inner_lostPixelEndcapHits_[i] = -999;
      muonCand_inner_lostStripHits_[i] = -999;
      muonCand_inner_lostStripTIBHits_[i] = -999;
      muonCand_inner_lostStripTIDHits_[i] = -999;
      muonCand_inner_lostStripTOBHits_[i] = -999;
      muonCand_inner_lostStripTECHits_[i] = -999;
      muonCand_inner_pixelLayers_[i] = -999;
      muonCand_inner_pixelHits_[i] = -999;
      muonCand_global_muonHits_[i] = -999;
      muonCand_global_trkChi2_[i] = -999;
      muonCand_global_trackerLayers_[i] = -999;
      muonCand_global_trackerHits_[i] = -999;
      
      muonCand_dxy_bs_[i] = -999;
      muonCand_dxyError_bs_[i] = -999;
      muonCand_dz_bs_[i] = -999;
      muonCand_dzError_[i] = -999;
      muonCand_IPSig_[i] = -999;      
    }

  
  nTrack_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    track_pt_[i] = -999;
    track_eta_[i] = -999;
    track_phi_[i] = -999;
    track_charge_[i] = -999;
    track_px_[i] = -999;
    track_py_[i] = -999;
    track_pz_[i] = -999;
    track_vx_[i] = -999;
    track_vy_[i] = -999;
    track_vz_[i] = -999;
    track_dxy_bs_[i] = -999;
    track_dxyError_bs_[i] = -999;
    track_dz_bs_[i] = -999;
    track_dzError_[i] = -999;
    track_trkChi2_[i] = -999;
    track_trackerLayers_[i] = -999;
    track_trackerHits_[i] = -999;
    track_lostTrackerHits_[i] = -999;
    track_lostTrackerHitsIn_[i] = -999;
    track_lostTrackerHitsOut_[i] = -999;
    track_lostPixelHits_[i] = -999;
    track_lostPixelBarrelHits_[i] = -999;
    track_lostPixelEndcapHits_[i] = -999;
    track_lostStripHits_[i] = -999;
    track_lostStripTIBHits_[i] = -999;
    track_lostStripTIDHits_[i] = -999;
    track_lostStripTOBHits_[i] = -999;
    track_lostStripTECHits_[i] = -999;
    track_pixelLayers_[i] = -999;
    track_pixelHits_[i] = -999;
    track_muonHits_[i] = -999;
    track_muonIdx_[i] = -999;
  }

  nECAL_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    ecal_et_[i] = -999;
    ecal_pt_[i] = -999;
    ecal_eta_[i] = -999;
    ecal_phi_[i] = -999;
    ecal_charge_[i] = -999;
    ecal_px_[i] = -999;
    ecal_py_[i] = -999;
    ecal_pz_[i] = -999;
    ecal_vx_[i] = -999;
    ecal_vy_[i] = -999;
    ecal_vz_[i] = -999;
    ecal_time_[i] = -999;
    ecal_depth_[i] = -999;
    ecal_rho_[i] = -999;
    ecal_muonIdx_[i] = -999;
    ecal_nHits_[i] = -999;
  }

  nECALHits_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    ecal_hit_energy_[i] = -999;
    ecal_hit_depth_[i] = -999;
    ecal_hit_time_[i] = -999;
    ecal_hit_pt2_[i] = -999;
    ecal_hit_x_[i] = -999;
    ecal_hit_y_[i] = -999;
    ecal_hit_z_[i] = -999;
    ecal_hit_eta_[i] = -999;
    ecal_hit_phi_[i] = -999;
    ecal_hit_fraction_[i] = -999;
    ecal_hit_idx_[i] = -999;
  }
  
  nHCAL_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    hcal_et_[i] = -999;
    hcal_pt_[i] = -999;
    hcal_eta_[i] = -999;
    hcal_phi_[i] = -999;
    hcal_charge_[i] = -999;
    hcal_px_[i] = -999;
    hcal_py_[i] = -999;
    hcal_pz_[i] = -999;
    hcal_vx_[i] = -999;
    hcal_vy_[i] = -999;
    hcal_vz_[i] = -999;
    hcal_time_[i] = -999;
    hcal_depth_[i] = -999;
    hcal_rho_[i] = -999;
    hcal_muonIdx_[i] = -999;
    hcal_nHits_[i] = -999;
  }

  nHCALHits_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    hcal_hit_energy_[i] = -999;
    hcal_hit_depth_[i] = -999;
    hcal_hit_time_[i] = -999;
    hcal_hit_pt2_[i] = -999;
    hcal_hit_x_[i] = -999;
    hcal_hit_y_[i] = -999;
    hcal_hit_z_[i] = -999;
    hcal_hit_eta_[i] = -999;
    hcal_hit_phi_[i] = -999;
    hcal_hit_fraction_[i] = -999;
    hcal_hit_idx_[i] = -999;
  }
      
  // -- original trigger objects -- //
  vec_firedTrigger_.clear();
  vec_filterName_.clear();
  vec_HLTObj_pt_.clear();
  vec_HLTObj_eta_.clear();
  vec_HLTObj_phi_.clear();

  
}

void MuonHLTNtupler::Make_Branch()
{
  ntuple_->Branch("isRealData", &isRealData_, "isRealData/O"); // -- O: boolean -- //
  ntuple_->Branch("runNum",&runNum_,"runNum/I");
  ntuple_->Branch("lumiBlockNum",&lumiBlockNum_,"lumiBlockNum/I");
  ntuple_->Branch("eventNum",&eventNum_,"eventNum/l"); // -- unsigned long long -- //
  ntuple_->Branch("bs_x0", &bs_x0_, "bs_x0/D");
  ntuple_->Branch("bs_y0", &bs_y0_, "bs_y0/D");
  ntuple_->Branch("bs_z0", &bs_z0_, "bs_z0/D");
  ntuple_->Branch("bs_sigmaZ", &bs_sigmaZ_, "bs_sigmaZ/D");
  ntuple_->Branch("bs_dxdz", &bs_dxdz_, "bs_dxdz/D");
  ntuple_->Branch("bs_dydz", &bs_dydz_, "bs_dydz/D");
  ntuple_->Branch("bs_x0Error", &bs_x0Error_, "bs_x0Error/D");
  ntuple_->Branch("bs_y0Error", &bs_y0Error_, "bs_y0Error/D");
  ntuple_->Branch("bs_z0Error", &bs_z0Error_, "bs_z0Error/D");
  ntuple_->Branch("bs_sigmaZ0Error", &bs_sigmaZ0Error_, "bs_sigmaZ0Error/D");
  ntuple_->Branch("bs_dxdzError", &bs_dxdzError_, "bs_dxdzError/D");
  ntuple_->Branch("bs_dydzError", &bs_dydzError_, "bs_dydzError/D");
  ntuple_->Branch("nVertex", &nVertex_, "nVertex/I");
  ntuple_->Branch("bunchID", &bunchID_, "bunchID/D");
  ntuple_->Branch("instLumi", &instLumi_, "instLumi/D");
  ntuple_->Branch("dataPU", &dataPU_, "dataPU/D");
  ntuple_->Branch("dataPURMS", &dataPURMS_, "dataPURMS/D");
  ntuple_->Branch("bunchLumi", &bunchLumi_, "bunchLumi/D");
  ntuple_->Branch("offlineInstLumi", &offlineInstLumi_, "offlineInstLumi/D");
  ntuple_->Branch("offlineDataPU", &offlineDataPU_, "offlineDataPU/D");
  ntuple_->Branch("offlineDataPURMS", &offlineDataPURMS_, "offlineDataPURMS/D");
  ntuple_->Branch("offlineBunchLumi", &offlineBunchLumi_, "offlineBunchLumi/D");
  ntuple_->Branch("truePU", &truePU_, "truePU/I");

  ntuple_->Branch("vertex_x", &vertex_x_, "vertex_x_[nVertex]/D");
  ntuple_->Branch("vertex_y", &vertex_y_, "vertex_y_[nVertex]/D");
  ntuple_->Branch("vertex_z", &vertex_z_, "vertex_z_[nVertex]/D");
  ntuple_->Branch("vertex_t", &vertex_t_, "vertex_t_[nVertex]/D");
  
  ntuple_->Branch("genEventWeight", &genEventWeight_, "genEventWeight/D");
  ntuple_->Branch("nGenParticle", &nGenParticle_, "nGenParticle/I");
  ntuple_->Branch("genParticle_ID", &genParticle_ID_, "genParticle_ID[nGenParticle]/I");
  ntuple_->Branch("genParticle_status", &genParticle_status_, "genParticle_status[nGenParticle]/I");
  ntuple_->Branch("genParticle_mother", &genParticle_mother_, "genParticle_mother[nGenParticle]/I");
  ntuple_->Branch("genParticle_pt", &genParticle_pt_, "genParticle_pt[nGenParticle]/D");
  ntuple_->Branch("genParticle_eta", &genParticle_eta_, "genParticle_eta[nGenParticle]/D");
  ntuple_->Branch("genParticle_phi", &genParticle_phi_, "genParticle_phi[nGenParticle]/D");
  ntuple_->Branch("genParticle_px", &genParticle_px_, "genParticle_px[nGenParticle]/D");
  ntuple_->Branch("genParticle_py", &genParticle_py_, "genParticle_py[nGenParticle]/D");
  ntuple_->Branch("genParticle_pz", &genParticle_pz_, "genParticle_pz[nGenParticle]/D");
  ntuple_->Branch("genParticle_energy", &genParticle_energy_, "genParticle_energy[nGenParticle]/D");
  ntuple_->Branch("genParticle_charge", &genParticle_charge_, "genParticle_charge[nGenParticle]/D");
  ntuple_->Branch("genParticle_isPrompt", &genParticle_isPrompt_, "genParticle_isPrompt[nGenParticle]/I");
  ntuple_->Branch("genParticle_isPromptFinalState", &genParticle_isPromptFinalState_, "genParticle_isPromptFinalState[nGenParticle]/I");
  ntuple_->Branch("genParticle_isTauDecayProduct", &genParticle_isTauDecayProduct_, "genParticle_isTauDecayProduct[nGenParticle]/I");
  ntuple_->Branch("genParticle_isPromptTauDecayProduct", &genParticle_isPromptTauDecayProduct_, "genParticle_isPromptTauDecayProduct[nGenParticle]/I");
  ntuple_->Branch("genParticle_isDirectPromptTauDecayProductFinalState", &genParticle_isDirectPromptTauDecayProductFinalState_, "genParticle_isDirectPromptTauDecayProductFinalState[nGenParticle]/I");
  ntuple_->Branch("genParticle_isHardProcess", &genParticle_isHardProcess_, "genParticle_isHardProcess[nGenParticle]/I");
  ntuple_->Branch("genParticle_isLastCopy", &genParticle_isLastCopy_, "genParticle_isLastCopy[nGenParticle]/I");
  ntuple_->Branch("genParticle_isLastCopyBeforeFSR", &genParticle_isLastCopyBeforeFSR_, "genParticle_isLastCopyBeforeFSR[nGenParticle]/I");
  ntuple_->Branch("genParticle_isPromptDecayed", &genParticle_isPromptDecayed_, "genParticle_isPromptDecayed[nGenParticle]/I");
  ntuple_->Branch("genParticle_isDecayedLeptonHadron", &genParticle_isDecayedLeptonHadron_, "genParticle_isDecayedLeptonHadron[nGenParticle]/I");
  ntuple_->Branch("genParticle_fromHardProcessBeforeFSR", &genParticle_fromHardProcessBeforeFSR_, "genParticle_fromHardProcessBeforeFSR[nGenParticle]/I");
  ntuple_->Branch("genParticle_fromHardProcessDecayed", &genParticle_fromHardProcessDecayed_, "genParticle_fromHardProcessDecayed[nGenParticle]/I");
  ntuple_->Branch("genParticle_fromHardProcessFinalState", &genParticle_fromHardProcessFinalState_, "genParticle_fromHardProcessFinalState[nGenParticle]/I");
  ntuple_->Branch("genParticle_isMostlyLikePythia6Status3", &genParticle_isMostlyLikePythia6Status3_, "genParticle_isMostlyLikePythia6Status3[nGenParticle]/I");


  ntuple_->Branch("nMuon", &nMuon_, "nMuon/I");

  ntuple_->Branch("muon_pt", &muon_pt_, "muon_pt[nMuon]/D");
  ntuple_->Branch("muon_eta", &muon_eta_, "muon_eta[nMuon]/D");
  ntuple_->Branch("muon_phi", &muon_phi_, "muon_phi[nMuon]/D");
  ntuple_->Branch("muon_px", &muon_px_, "muon_px[nMuon]/D");
  ntuple_->Branch("muon_py", &muon_py_, "muon_py[nMuon]/D");
  ntuple_->Branch("muon_pz", &muon_pz_, "muon_pz[nMuon]/D");
  ntuple_->Branch("muon_dB", &muon_dB_, "muon_dB[nMuon]/D");
  ntuple_->Branch("muon_charge", &muon_charge_, "muon_charge[nMuon]/D");
  ntuple_->Branch("muon_isGLB", &muon_isGLB_, "muon_isGLB[nMuon]/I");
  ntuple_->Branch("muon_isSTA", &muon_isSTA_, "muon_isSTA[nMuon]/I");
  ntuple_->Branch("muon_isTRK", &muon_isTRK_, "muon_isTRK[nMuon]/I");
  ntuple_->Branch("muon_isPF", &muon_isPF_, "muon_isPF[nMuon]/I");
  ntuple_->Branch("muon_isTight", &muon_isTight_, "muon_isTight[nMuon]/I");
  ntuple_->Branch("muon_isMedium", &muon_isMedium_, "muon_isMedium[nMuon]/I");
  ntuple_->Branch("muon_isLoose", &muon_isLoose_, "muon_isLoose[nMuon]/I");
  ntuple_->Branch("muon_isHighPt", &muon_isHighPt_, "muon_isHighPt[nMuon]/I");
  ntuple_->Branch("muon_isHighPtNew", &muon_isHighPtNew_, "muon_isHighPtNew[nMuon]/I");
  ntuple_->Branch("muon_isSoft", &muon_isSoft_, "muon_isSoft[nMuon]/I");

  ntuple_->Branch("muon_iso03_sumPt", &muon_iso03_sumPt_, "muon_iso03_sumPt[nMuon]/D");
  ntuple_->Branch("muon_iso03_hadEt", &muon_iso03_hadEt_, "muon_iso03_hadEt[nMuon]/D");
  ntuple_->Branch("muon_iso03_emEt", &muon_iso03_emEt_, "muon_iso03_emEt[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_charged", &muon_PFIso03_charged_, "muon_PFIso03_charged[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_neutral", &muon_PFIso03_neutral_, "muon_PFIso03_neutral[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_photon", &muon_PFIso03_photon_, "muon_PFIso03_photon[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_sumPU", &muon_PFIso03_sumPU_, "muon_PFIso03_sumPU[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_charged", &muon_PFIso04_charged_, "muon_PFIso04_charged[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_neutral", &muon_PFIso04_neutral_, "muon_PFIso04_neutral[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_photon", &muon_PFIso04_photon_, "muon_PFIso04_photon[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_sumPU", &muon_PFIso04_sumPU_, "muon_PFIso04_sumPU[nMuon]/D");

  ntuple_->Branch("muon_PFCluster03_ECAL", &muon_PFCluster03_ECAL_, "muon_PFCluster03_ECAL[nMuon]/D");
  ntuple_->Branch("muon_PFCluster03_HCAL", &muon_PFCluster03_HCAL_, "muon_PFCluster03_HCAL[nMuon]/D");
  ntuple_->Branch("muon_PFCluster04_ECAL", &muon_PFCluster04_ECAL_, "muon_PFCluster04_ECAL[nMuon]/D");
  ntuple_->Branch("muon_PFCluster04_HCAL", &muon_PFCluster04_HCAL_, "muon_PFCluster04_HCAL[nMuon]/D");

  ntuple_->Branch("muon_inner_trkChi2", &muon_inner_trkChi2_, "muon_inner_trkChi2[nMuon]/D");
  ntuple_->Branch("muon_inner_validFraction", &muon_inner_validFraction_, "muon_inner_validFraction[nMuon]/D");
  ntuple_->Branch("muon_inner_trackerLayers", &muon_inner_trackerLayers_, "muon_inner_trackerLayers[nMuon]/I");
  ntuple_->Branch("muon_inner_trackerHits", &muon_inner_trackerHits_, "muon_inner_trackerHits[nMuon]/I");
  ntuple_->Branch("muon_inner_lostTrackerHits", &muon_inner_lostTrackerHits_, "muon_inner_lostTrackerHits[nMuon]/I");
  ntuple_->Branch("muon_inner_lostTrackerHitsIn", &muon_inner_lostTrackerHitsIn_, "muon_inner_lostTrackerHitsIn[nMuon]/I");
  ntuple_->Branch("muon_inner_lostTrackerHitsOut", &muon_inner_lostTrackerHitsOut_, "muon_inner_lostTrackerHitsOut[nMuon]/I");
  ntuple_->Branch("muon_inner_lostPixelHits", &muon_inner_lostPixelHits_, "muon_inner_lostPixelHits[nMuon]/I");
  ntuple_->Branch("muon_inner_lostPixelBarrelHits", &muon_inner_lostPixelBarrelHits_, "muon_inner_lostPixelBarrelHits[nMuon]/I");
  ntuple_->Branch("muon_inner_lostPixelEndcapHits", &muon_inner_lostPixelEndcapHits_, "muon_inner_lostPixelEndcapHits[nMuon]/I");
  ntuple_->Branch("muon_inner_lostStripHits", &muon_inner_lostStripHits_, "muon_inner_lostStripHits[nMuon]/I");
  ntuple_->Branch("muon_inner_lostStripTIBHits", &muon_inner_lostStripTIBHits_, "muon_inner_lostStripTIBHits[nMuon]/I");
  ntuple_->Branch("muon_inner_lostStripTIDHits", &muon_inner_lostStripTIDHits_, "muon_inner_lostStripTIDHits[nMuon]/I");
  ntuple_->Branch("muon_inner_lostStripTOBHits", &muon_inner_lostStripTOBHits_, "muon_inner_lostStripTOBHits[nMuon]/I");
  ntuple_->Branch("muon_inner_lostStripTECHits", &muon_inner_lostStripTECHits_, "muon_inner_lostStripTECHits[nMuon]/I");
  ntuple_->Branch("muon_inner_pixelLayers", &muon_inner_pixelLayers_, "muon_inner_pixelLayers[nMuon]/I");
  ntuple_->Branch("muon_inner_pixelHits", &muon_inner_pixelHits_, "muon_inner_pixelHits[nMuon]/I");
  ntuple_->Branch("muon_global_muonHits", &muon_global_muonHits_, "muon_global_muonHits[nMuon]/I");
  ntuple_->Branch("muon_global_trkChi2", &muon_global_trkChi2_, "muon_global_trkChi2[nMuon]/D");
  ntuple_->Branch("muon_global_trackerLayers", &muon_global_trackerLayers_, "muon_global_trackerLayers[nMuon]/I");
  ntuple_->Branch("muon_global_trackerHits", &muon_global_trackerHits_, "muon_global_trackerHits[nMuon]/I");
  ntuple_->Branch("muon_momentumChi2", &muon_momentumChi2_, "muon_momentumChi2[nMuon]/D");
  ntuple_->Branch("muon_positionChi2", &muon_positionChi2_, "muon_positionChi2[nMuon]/D");
  ntuple_->Branch("muon_glbKink", &muon_glbKink_, "muon_glbKink[nMuon]/D");
  ntuple_->Branch("muon_glbTrackProbability", &muon_glbTrackProbability_, "muon_glbTrackProbability[nMuon]/D");
  ntuple_->Branch("muon_globalDeltaEtaPhi", &muon_globalDeltaEtaPhi_, "muon_globalDeltaEtaPhi[nMuon]/D");
  ntuple_->Branch("muon_localDistance", &muon_localDistance_, "muon_localDistance[nMuon]/D");
  ntuple_->Branch("muon_staRelChi2", &muon_staRelChi2_, "muon_staRelChi2[nMuon]/D");
  ntuple_->Branch("muon_tightMatch", &muon_tightMatch_, "muon_tightMatch[nMuon]/I");
  ntuple_->Branch("muon_trkKink", &muon_trkKink_, "muon_trkKink[nMuon]/D");
  ntuple_->Branch("muon_trkRelChi2", &muon_trkRelChi2_, "muon_trkRelChi2[nMuon]/D");
  ntuple_->Branch("muon_segmentCompatibility", &muon_segmentCompatibility_, "muon_segmentCompatibility[nMuon]/D");

  ntuple_->Branch("muon_pt_tuneP", &muon_pt_tuneP_, "muon_pt_tuneP[nMuon]/D");
  ntuple_->Branch("muon_ptError_tuneP", &muon_ptError_tuneP_, "muon_ptError_tuneP[nMuon]/D");
  ntuple_->Branch("muon_dxyVTX_best", &muon_dxyVTX_best_, "muon_dxyVTX_best[nMuon]/D");
  ntuple_->Branch("muon_dzVTX_best", &muon_dzVTX_best_, "muon_dzVTX_best[nMuon]/D");
  ntuple_->Branch("muon_nMatchedStation", &muon_nMatchedStation_, "muon_nMatchedStation[nMuon]/I");
  ntuple_->Branch("muon_nMatchedRPCLayer", &muon_nMatchedRPCLayer_, "muon_nMatchedRPCLayer[nMuon]/I");
  ntuple_->Branch("muon_stationMask", &muon_stationMask_, "muon_stationMask[nMuon]/I");
  ntuple_->Branch("muon_dxy_bs", &muon_dxy_bs_, "muon_dxy_bs[nMuon]/D");
  ntuple_->Branch("muon_dxyError_bs", &muon_dxyError_bs_, "muon_dxyError_bs[nMuon]/D");
  ntuple_->Branch("muon_dz_bs", &muon_dz_bs_, "muon_dz_bs[nMuon]/D");
  ntuple_->Branch("muon_dzError", &muon_dzError_, "muon_dzError[nMuon]/D");
  ntuple_->Branch("muon_IPSig", &muon_IPSig_, "muon_IPSig[nMuon]/D");


  ntuple_->Branch("nMuonCand", &nMuonCand_, "nMuonCand/I");

  ntuple_->Branch("muonCand_pt", &muonCand_pt_, "muonCand_pt[nMuonCand]/D");
  ntuple_->Branch("muonCand_eta", &muonCand_eta_, "muonCand_eta[nMuonCand]/D");
  ntuple_->Branch("muonCand_phi", &muonCand_phi_, "muonCand_phi[nMuonCand]/D");
  ntuple_->Branch("muonCand_px", &muonCand_px_, "muonCand_px[nMuonCand]/D");
  ntuple_->Branch("muonCand_py", &muonCand_py_, "muonCand_py[nMuonCand]/D");
  ntuple_->Branch("muonCand_pz", &muonCand_pz_, "muonCand_pz[nMuonCand]/D");
  ntuple_->Branch("muonCand_vx", &muonCand_vx_, "muonCand_vx[nMuonCand]/D");
  ntuple_->Branch("muonCand_vy", &muonCand_vy_, "muonCand_vy[nMuonCand]/D");
  ntuple_->Branch("muonCand_vz", &muonCand_vz_, "muonCand_vz[nMuonCand]/D");
  ntuple_->Branch("muonCand_dB", &muonCand_dB_, "muonCand_dB[nMuonCand]/D");
  ntuple_->Branch("muonCand_charge", &muonCand_charge_, "muonCand_charge[nMuonCand]/D");
  ntuple_->Branch("muonCand_isGLB", &muonCand_isGLB_, "muonCand_isGLB[nMuonCand]/I");
  ntuple_->Branch("muonCand_isSTA", &muonCand_isSTA_, "muonCand_isSTA[nMuonCand]/I");
  ntuple_->Branch("muonCand_isTRK", &muonCand_isTRK_, "muonCand_isTRK[nMuonCand]/I");
  ntuple_->Branch("muonCand_isMuon", &muonCand_isMuon_, "muonCand_isMuon[nMuonCand]/I");

  ntuple_->Branch("muonCand_inner_trkChi2", &muonCand_inner_trkChi2_, "muonCand_inner_trkChi2[nMuonCand]/D");
  ntuple_->Branch("muonCand_inner_validFraction", &muonCand_inner_validFraction_, "muonCand_inner_validFraction[nMuonCand]/D");
  ntuple_->Branch("muonCand_inner_trackerLayers", &muonCand_inner_trackerLayers_, "muonCand_inner_trackerLayers[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_trackerHits", &muonCand_inner_trackerHits_, "muonCand_inner_trackerHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostTrackerHits", &muonCand_inner_lostTrackerHits_, "muonCand_inner_lostTrackerHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostTrackerHitsIn", &muonCand_inner_lostTrackerHitsIn_, "muonCand_inner_lostTrackerHitsIn[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostTrackerHitsOut", &muonCand_inner_lostTrackerHitsOut_, "muonCand_inner_lostTrackerHitsOut[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostPixelHits", &muonCand_inner_lostPixelHits_, "muonCand_inner_lostPixelHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostPixelBarrelHits", &muonCand_inner_lostPixelBarrelHits_, "muonCand_inner_lostPixelBarrelHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostPixelEndcapHits", &muonCand_inner_lostPixelEndcapHits_, "muonCand_inner_lostPixelEndcapHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripHits", &muonCand_inner_lostStripHits_, "muonCand_inner_lostStripHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripTIBHits", &muonCand_inner_lostStripTIBHits_, "muonCand_inner_lostStripTIBHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripTIDHits", &muonCand_inner_lostStripTIDHits_, "muonCand_inner_lostStripTIDHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripTOBHits", &muonCand_inner_lostStripTOBHits_, "muonCand_inner_lostStripTOBHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripTECHits", &muonCand_inner_lostStripTECHits_, "muonCand_inner_lostStripTECHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_pixelLayers", &muonCand_inner_pixelLayers_, "muonCand_inner_pixelLayers[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_pixelHits", &muonCand_inner_pixelHits_, "muonCand_inner_pixelHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_global_muonHits", &muonCand_global_muonHits_, "muonCand_global_muonHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_global_trkChi2", &muonCand_global_trkChi2_, "muonCand_global_trkChi2[nMuonCand]/D");
  ntuple_->Branch("muonCand_global_trackerLayers", &muonCand_global_trackerLayers_, "muonCand_global_trackerLayers[nMuonCand]/I");
  ntuple_->Branch("muonCand_global_trackerHits", &muonCand_global_trackerHits_, "muonCand_global_trackerHits[nMuonCand]/I");

  ntuple_->Branch("muonCand_dxy_bs", &muonCand_dxy_bs_, "muonCand_dxy_bs[nMuonCand]/D");
  ntuple_->Branch("muonCand_dxyError_bs", &muonCand_dxyError_bs_, "muonCand_dxyError_bs[nMuonCand]/D");
  ntuple_->Branch("muonCand_dz_bs", &muonCand_dz_bs_, "muonCand_dz_bs[nMuonCand]/D");
  ntuple_->Branch("muonCand_dzError", &muonCand_dzError_, "muonCand_dzError[nMuonCand]/D");
  ntuple_->Branch("muonCand_IPSig", &muonCand_IPSig_, "muonCand_IPSig[nMuonCand]/D");
  
  
  ntuple_->Branch("nTrack", &nTrack_, "nTrack/I");
  ntuple_->Branch("track_pt", &track_pt_, "track_pt[nTrack]/D");
  ntuple_->Branch("track_eta", &track_eta_, "track_eta[nTrack]/D");
  ntuple_->Branch("track_phi", &track_phi_, "track_phi[nTrack]/D");
  ntuple_->Branch("track_charge", &track_charge_, "track_charge[nTrack]/D");
  ntuple_->Branch("track_px", &track_px_, "track_px[nTrack]/D");
  ntuple_->Branch("track_py", &track_py_, "track_py[nTrack]/D");
  ntuple_->Branch("track_pz", &track_pz_, "track_pz[nTrack]/D");
  ntuple_->Branch("track_vx", &track_vx_, "track_vx[nTrack]/D");
  ntuple_->Branch("track_vy", &track_vy_, "track_vy[nTrack]/D");
  ntuple_->Branch("track_vz", &track_vz_, "track_vz[nTrack]/D");
  ntuple_->Branch("track_dxy_bs", &track_dxy_bs_, "track_dxy_bs[nTrack]/D");
  ntuple_->Branch("track_dxyError_bs", &track_dxyError_bs_, "track_dxyError_bs[nTrack]/D");
  ntuple_->Branch("track_dz_bs", &track_dz_bs_, "track_dz_bs[nTrack]/D");
  ntuple_->Branch("track_dzError", &track_dzError_, "track_dzError[nTrack]/D");
  ntuple_->Branch("track_trkChi2", &track_trkChi2_, "track_trkChi2[nTrack]/D");
  ntuple_->Branch("track_trackerLayers", &track_trackerLayers_, "track_trackerLayers[nTrack]/I");
  ntuple_->Branch("track_trackerHits", &track_trackerHits_, "track_trackerHits[nTrack]/I");
  ntuple_->Branch("track_lostTrackerHits", &track_lostTrackerHits_, "track_lostTrackerHits[nTrack]/I");
  ntuple_->Branch("track_lostTrackerHitsIn", &track_lostTrackerHitsIn_, "track_lostTrackerHitsIn[nTrack]/I");
  ntuple_->Branch("track_lostTrackerHitsOut", &track_lostTrackerHitsOut_, "track_lostTrackerHitsOut[nTrack]/I");
  ntuple_->Branch("track_lostPixelHits", &track_lostPixelHits_, "track_lostPixelHits[nTrack]/I");
  ntuple_->Branch("track_lostPixelBarrelHits", &track_lostPixelBarrelHits_, "track_lostPixelBarrelHits[nTrack]/I");
  ntuple_->Branch("track_lostPixelEndcapHits", &track_lostPixelEndcapHits_, "track_lostPixelEndcapHits[nTrack]/I");
  ntuple_->Branch("track_lostStripHits", &track_lostStripHits_, "track_lostStripHits[nTrack]/I");
  ntuple_->Branch("track_lostStripTIBHits", &track_lostStripTIBHits_, "track_lostStripTIBHits[nTrack]/I");
  ntuple_->Branch("track_lostStripTIDHits", &track_lostStripTIDHits_, "track_lostStripTIDHits[nTrack]/I");
  ntuple_->Branch("track_lostStripTOBHits", &track_lostStripTOBHits_, "track_lostStripTOBHits[nTrack]/I");
  ntuple_->Branch("track_lostStripTECHits", &track_lostStripTECHits_, "track_lostStripTECHits[nTrack]/I");
  ntuple_->Branch("track_pixelLayers", &track_pixelLayers_, "track_pixelLayers[nTrack]/I");
  ntuple_->Branch("track_pixelHits", &track_pixelHits_, "track_pixelHits[nTrack]/I");
  ntuple_->Branch("track_muonHits", &track_muonHits_, "track_muonHits[nTrack]/I");
  ntuple_->Branch("track_muonIdx", &track_muonIdx_, "track_muonIdx[nTrack]/I");

  ntuple_->Branch("nECAL", &nECAL_, "nECAL/I");
  ntuple_->Branch("ecal_et", &ecal_et_, "ecal_et_[nECAL]/D");
  ntuple_->Branch("ecal_pt", &ecal_pt_, "ecal_pt_[nECAL]/D");
  ntuple_->Branch("ecal_eta", &ecal_eta_, "ecal_eta_[nECAL]/D");
  ntuple_->Branch("ecal_phi", &ecal_phi_, "ecal_phi_[nECAL]/D");
  ntuple_->Branch("ecal_charge", &ecal_charge_, "ecal_charge_[nECAL]/D");
  ntuple_->Branch("ecal_px", &ecal_px_, "ecal_px_[nECAL]/D");
  ntuple_->Branch("ecal_py", &ecal_py_, "ecal_py_[nECAL]/D");
  ntuple_->Branch("ecal_pz", &ecal_pz_, "ecal_pz_[nECAL]/D");
  ntuple_->Branch("ecal_vx", &ecal_vx_, "ecal_vx_[nECAL]/D");
  ntuple_->Branch("ecal_vy", &ecal_vy_, "ecal_vy_[nECAL]/D");
  ntuple_->Branch("ecal_vz", &ecal_vz_, "ecal_vz_[nECAL]/D");
  ntuple_->Branch("ecal_time", &ecal_time_, "ecal_time_[nECAL]/D");
  ntuple_->Branch("ecal_depth", &ecal_depth_, "ecal_depth_[nECAL]/D");
  ntuple_->Branch("ecal_rho", &ecal_rho_, "ecal_rho_[nECAL]/D");
  ntuple_->Branch("ecal_muonIdx", &ecal_muonIdx_, "ecal_muonIdx_[nECAL]/I");
  ntuple_->Branch("ecal_nHits", &ecal_nHits_, "ecal_nHits_[nECAL]/I");
  
  ntuple_->Branch("nECALHits", &nECALHits_, "nECALHits/I");
  ntuple_->Branch("ecal_hit_energy", &ecal_hit_energy_, "ecal_hit_energy_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_depth", &ecal_hit_depth_, "ecal_hit_depth_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_time", &ecal_hit_time_, "ecal_hit_time_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_pt2", &ecal_hit_pt2_, "ecal_hit_pt2_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_x", &ecal_hit_x_, "ecal_hit_x_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_y", &ecal_hit_y_, "ecal_hit_y_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_z", &ecal_hit_z_, "ecal_hit_z_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_eta", &ecal_hit_eta_, "ecal_hit_eta_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_phi", &ecal_hit_phi_, "ecal_hit_phi_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_fraction", &ecal_hit_fraction_, "ecal_hit_fraction_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_idx", &ecal_hit_idx_, "ecal_hit_idx_[nECALHits]/D");

  ntuple_->Branch("nHCAL", &nHCAL_, "nHCAL/I");
  ntuple_->Branch("hcal_et", &hcal_et_, "hcal_et_[nHCAL]/D");
  ntuple_->Branch("hcal_pt", &hcal_pt_, "hcal_pt_[nHCAL]/D");
  ntuple_->Branch("hcal_eta", &hcal_eta_, "hcal_eta_[nHCAL]/D");
  ntuple_->Branch("hcal_phi", &hcal_phi_, "hcal_phi_[nHCAL]/D");
  ntuple_->Branch("hcal_charge", &hcal_charge_, "hcal_charge_[nHCAL]/D");
  ntuple_->Branch("hcal_px", &hcal_px_, "hcal_px_[nHCAL]/D");
  ntuple_->Branch("hcal_py", &hcal_py_, "hcal_py_[nHCAL]/D");
  ntuple_->Branch("hcal_pz", &hcal_pz_, "hcal_pz_[nHCAL]/D");
  ntuple_->Branch("hcal_vx", &hcal_vx_, "hcal_vx_[nHCAL]/D");
  ntuple_->Branch("hcal_vy", &hcal_vy_, "hcal_vy_[nHCAL]/D");
  ntuple_->Branch("hcal_vz", &hcal_vz_, "hcal_vz_[nHCAL]/D");
  ntuple_->Branch("hcal_time", &hcal_time_, "hcal_time_[nHCAL]/D");
  ntuple_->Branch("hcal_depth", &hcal_depth_, "hcal_depth_[nHCAL]/D");
  ntuple_->Branch("hcal_rho", &hcal_rho_, "hcal_rho_[nHCAL]/D");
  ntuple_->Branch("hcal_muonIdx", &hcal_muonIdx_, "hcal_muonIdx_[nHCAL]/I");
  ntuple_->Branch("hcal_nHits", &hcal_nHits_, "hcal_nHits_[nHCAL]/I");

  ntuple_->Branch("nHCALHits", &nHCALHits_, "nHCALHits/I");
  ntuple_->Branch("hcal_hit_energy", &hcal_hit_energy_, "hcal_hit_energy_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_depth", &hcal_hit_depth_, "hcal_hit_depth_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_time", &hcal_hit_time_, "hcal_hit_time_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_pt2", &hcal_hit_pt2_, "hcal_hit_pt2_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_x", &hcal_hit_x_, "hcal_hit_x_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_y", &hcal_hit_y_, "hcal_hit_y_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_z", &hcal_hit_z_, "hcal_hit_z_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_eta", &hcal_hit_eta_, "hcal_hit_eta_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_phi", &hcal_hit_phi_, "hcal_hit_phi_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_fraction", &hcal_hit_fraction_, "hcal_hit_fraction_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_idx", &hcal_hit_idx_, "hcal_hit_idx_[nHCALHits]/D");
  
  ntuple_->Branch("vec_firedTrigger", &vec_firedTrigger_);
  ntuple_->Branch("vec_filterName", &vec_filterName_);
  ntuple_->Branch("vec_HLTObj_pt", &vec_HLTObj_pt_);
  ntuple_->Branch("vec_HLTObj_eta", &vec_HLTObj_eta_);
  ntuple_->Branch("vec_HLTObj_phi", &vec_HLTObj_phi_);

}

void MuonHLTNtupler::Fill_Muon(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  auto const prop = propSetup_.init(iSetup);

  edm::Handle<std::vector<reco::Muon> > h_offlineMuon;  
  if( iEvent.getByToken(t_offlineMuon_, h_offlineMuon) ) // -- only when the dataset has offline muon collection (e.g. AOD) -- //
  {
    edm::Handle<reco::VertexCollection> h_offlineVertex;
    iEvent.getByToken(t_offlineVertex_, h_offlineVertex);
    const reco::Vertex & pv = h_offlineVertex->at(0);

    int _nMuon = 0;
    for(std::vector<reco::Muon>::const_iterator mu=h_offlineMuon->begin(); mu!=h_offlineMuon->end(); ++mu)
    {
      muon_pt_[_nMuon]  = mu->pt();
      muon_eta_[_nMuon] = mu->eta();
      muon_phi_[_nMuon] = mu->phi();
      muon_px_[_nMuon]  = mu->px();
      muon_py_[_nMuon]  = mu->py();
      muon_pz_[_nMuon]  = mu->pz();
      muon_charge_[_nMuon] = mu->charge();

      if( mu->isGlobalMuon() ) muon_isGLB_[_nMuon] = 1;
      if( mu->isStandAloneMuon() ) muon_isSTA_[_nMuon] = 1;
      if( mu->isTrackerMuon() ) muon_isTRK_[_nMuon] = 1;
      if( mu->isPFMuon() ) muon_isPF_[_nMuon] = 1;

      // -- defintion of ID functions: http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_4_0/doc/html/da/d18/namespacemuon.html#ac122b2516e5711ce206256d7945473d2 -- //
      if( muon::isTightMuon( (*mu), pv ) )  muon_isTight_[_nMuon] = 1;
      if( muon::isMediumMuon( (*mu) ) )     muon_isMedium_[_nMuon] = 1;
      if( muon::isLooseMuon( (*mu) ) )      muon_isLoose_[_nMuon] = 1;
      if( muon::isHighPtMuon( (*mu), pv ) ) muon_isHighPt_[_nMuon] = 1;
      if( isNewHighPtMuon( (*mu), pv ) )    muon_isHighPtNew_[_nMuon] = 1;

      // -- bool muon::isSoftMuon(const reco::Muon& muon, const reco::Vertex& vtx, bool run2016_hip_mitigation)
      // -- it is different under CMSSW_8_0_29: bool muon::isSoftMuon(const reco::Muon& muon, const reco::Vertex& vtx)
      // -- Remove this part to avoid compile error (and soft muon would not be used for now) - need to be fixed at some point
      // if( muon::isSoftMuon( (*mu), pv, 0) ) muon_isSoft_[_nMuon] = 1;

      muon_iso03_sumPt_[_nMuon] = mu->isolationR03().sumPt;
      muon_iso03_hadEt_[_nMuon] = mu->isolationR03().hadEt;
      muon_iso03_emEt_[_nMuon]  = mu->isolationR03().emEt;

      muon_PFIso03_charged_[_nMuon] = mu->pfIsolationR03().sumChargedHadronPt;
      muon_PFIso03_neutral_[_nMuon] = mu->pfIsolationR03().sumNeutralHadronEt;
      muon_PFIso03_photon_[_nMuon]  = mu->pfIsolationR03().sumPhotonEt;
      muon_PFIso03_sumPU_[_nMuon]   = mu->pfIsolationR03().sumPUPt;

      muon_PFIso04_charged_[_nMuon] = mu->pfIsolationR04().sumChargedHadronPt;
      muon_PFIso04_neutral_[_nMuon] = mu->pfIsolationR04().sumNeutralHadronEt;
      muon_PFIso04_photon_[_nMuon]  = mu->pfIsolationR04().sumPhotonEt;
      muon_PFIso04_sumPU_[_nMuon]   = mu->pfIsolationR04().sumPUPt;

      reco::MuonRef muRef = reco::MuonRef(h_offlineMuon, _nMuon);

      reco::TrackRef innerTrk = mu->innerTrack();
      if( innerTrk.isNonnull() )
        {
          muon_dxy_bs_[_nMuon] = innerTrk->dxy(bs->position());
          muon_dxyError_bs_[_nMuon] = innerTrk->dxyError(*bs);
          muon_dz_bs_[_nMuon] = innerTrk->dz(bs->position());
          muon_dzError_[_nMuon] = innerTrk->dzError();
          if (innerTrk->dxyError(*bs) > 0.) {
            muon_IPSig_[_nMuon] = abs(innerTrk->dxy(bs->position()) / innerTrk->dxyError(*bs));
          }
          muon_inner_trkChi2_[_nMuon]             = innerTrk->normalizedChi2();
          muon_inner_validFraction_[_nMuon]       = innerTrk->validFraction();
          muon_inner_trackerLayers_[_nMuon]       = innerTrk->hitPattern().trackerLayersWithMeasurement();
          muon_inner_trackerHits_[_nMuon]         = innerTrk->hitPattern().numberOfValidTrackerHits();
          muon_inner_lostTrackerHits_[_nMuon]     = innerTrk->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS);
          muon_inner_lostTrackerHitsIn_[_nMuon]   = innerTrk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
          muon_inner_lostTrackerHitsOut_[_nMuon]  = innerTrk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS);
          muon_inner_lostPixelHits_[_nMuon]       = innerTrk->hitPattern().numberOfLostPixelHits(HitPattern::TRACK_HITS);
          muon_inner_lostPixelBarrelHits_[_nMuon] = innerTrk->hitPattern().numberOfLostPixelBarrelHits(HitPattern::TRACK_HITS);
          muon_inner_lostPixelEndcapHits_[_nMuon] = innerTrk->hitPattern().numberOfLostPixelEndcapHits(HitPattern::TRACK_HITS);
          muon_inner_lostStripHits_[_nMuon]       = innerTrk->hitPattern().numberOfLostStripHits(HitPattern::TRACK_HITS);
          muon_inner_lostStripTIBHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTIBHits(HitPattern::TRACK_HITS);
          muon_inner_lostStripTIDHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTIDHits(HitPattern::TRACK_HITS);
          muon_inner_lostStripTOBHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTOBHits(HitPattern::TRACK_HITS);
          muon_inner_lostStripTECHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTECHits(HitPattern::TRACK_HITS);
          muon_inner_pixelLayers_[_nMuon]         = innerTrk->hitPattern().pixelLayersWithMeasurement();
          muon_inner_pixelHits_[_nMuon]           = innerTrk->hitPattern().numberOfValidPixelHits();
        }

      reco::TrackRef globalTrk = mu->globalTrack();
      if( globalTrk.isNonnull() )
      {
        muon_global_muonHits_[_nMuon]           = globalTrk->hitPattern().numberOfValidMuonHits();
        muon_global_trkChi2_[_nMuon]            = globalTrk->normalizedChi2();
        muon_global_trackerLayers_[_nMuon]      = globalTrk->hitPattern().trackerLayersWithMeasurement();
        muon_global_trackerHits_[_nMuon]        = globalTrk->hitPattern().numberOfValidTrackerHits();
      }
      muon_momentumChi2_[_nMuon]         = mu->combinedQuality().chi2LocalMomentum;
      muon_positionChi2_[_nMuon]         = mu->combinedQuality().chi2LocalPosition;
      muon_glbKink_[_nMuon]              = mu->combinedQuality().glbKink;
      muon_glbTrackProbability_[_nMuon]  = mu->combinedQuality().glbTrackProbability;
      muon_globalDeltaEtaPhi_[_nMuon]    = mu->combinedQuality().globalDeltaEtaPhi;
      muon_localDistance_[_nMuon]        = mu->combinedQuality().localDistance;
      muon_staRelChi2_[_nMuon]           = mu->combinedQuality().staRelChi2;
      muon_tightMatch_[_nMuon]           = mu->combinedQuality().tightMatch;
      muon_trkKink_[_nMuon]              = mu->combinedQuality().trkKink;
      muon_trkRelChi2_[_nMuon]           = mu->combinedQuality().trkRelChi2;
      muon_segmentCompatibility_[_nMuon] = muon::segmentCompatibility(*mu);

      reco::TrackRef tunePTrk = mu->tunePMuonBestTrack();
      if( tunePTrk.isNonnull() )
      {
        muon_pt_tuneP_[_nMuon]      = tunePTrk->pt();
        muon_ptError_tuneP_[_nMuon] = tunePTrk->ptError();
      }

      muon_dxyVTX_best_[_nMuon] = mu->muonBestTrack()->dxy( pv.position() );
      muon_dzVTX_best_[_nMuon]  = mu->muonBestTrack()->dz( pv.position() );

      muon_nMatchedStation_[_nMuon] = mu->numberOfMatchedStations();
      muon_nMatchedRPCLayer_[_nMuon] = mu->numberOfMatchedRPCLayers();
      muon_stationMask_[_nMuon] = mu->stationMask();
      
      _nMuon++;
    }

    nMuon_ = _nMuon;
  }
}



void MuonHLTNtupler::Fill_HLT(const edm::Event &iEvent, bool isMYHLT)
{
  edm::Handle<edm::TriggerResults>  h_triggerResults;
  edm::Handle<trigger::TriggerEvent> h_triggerEvent;

  if( isMYHLT )
  {
    iEvent.getByToken(t_myTriggerResults_, h_triggerResults);
    iEvent.getByToken(t_myTriggerEvent_,   h_triggerEvent);
  }
  else
  {
    iEvent.getByToken(t_triggerResults_, h_triggerResults);
    iEvent.getByToken(t_triggerEvent_,   h_triggerEvent);
  }

  edm::TriggerNames triggerNames = iEvent.triggerNames(*h_triggerResults);

  for(unsigned int itrig=0; itrig<triggerNames.size(); ++itrig)
  {
    LogDebug("triggers") << triggerNames.triggerName(itrig);

    if( h_triggerResults->accept(itrig) )
    {
      std::string pathName = triggerNames.triggerName(itrig);
      if( SavedTriggerCondition(pathName) || isMYHLT )
      {
        if( isMYHLT ) vec_myFiredTrigger_.push_back( pathName );
        else          vec_firedTrigger_.push_back( pathName );
      }
    } // -- end of if fired -- //

  } // -- end of iteration over all trigger names -- //

  const trigger::size_type nFilter(h_triggerEvent->sizeFilters());
  for( trigger::size_type i_filter=0; i_filter<nFilter; i_filter++)
  {
    std::string filterName = h_triggerEvent->filterTag(i_filter).encode();

    if( SavedFilterCondition(filterName) || isMYHLT )
    {
      trigger::Keys objectKeys = h_triggerEvent->filterKeys(i_filter);
      const trigger::TriggerObjectCollection& triggerObjects(h_triggerEvent->getObjects());

      for( trigger::size_type i_key=0; i_key<objectKeys.size(); i_key++)
      {
        trigger::size_type objKey = objectKeys.at(i_key);
        const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);

        if( isMYHLT )
        {
          vec_myFilterName_.push_back( filterName );
          vec_myHLTObj_pt_.push_back( triggerObj.pt() );
          vec_myHLTObj_eta_.push_back( triggerObj.eta() );
          vec_myHLTObj_phi_.push_back( triggerObj.phi() );
        }
        else
        {
          vec_filterName_.push_back( filterName );
          vec_HLTObj_pt_.push_back( triggerObj.pt() );
          vec_HLTObj_eta_.push_back( triggerObj.eta() );
          vec_HLTObj_phi_.push_back( triggerObj.phi() );
        }
      }
    } // -- end of if( muon filters )-- //
  } // -- end of filter iteration -- //
}

bool MuonHLTNtupler::SavedTriggerCondition( std::string& pathName )
{
  bool flag = false;

  // -- muon triggers
  if( pathName.find("Mu")           != std::string::npos ||
      pathName.find("HLT_IsoMu")    != std::string::npos ||
      pathName.find("HLT_Mu")       != std::string::npos ||
      pathName.find("HLT_OldMu")    != std::string::npos ||
      pathName.find("HLT_TkMu")     != std::string::npos ||
      pathName.find("HLT_IsoTkMu")  != std::string::npos ||
      pathName.find("HLT_DoubleMu") != std::string::npos ||
      pathName.find("HLT_Mu8_T")    != std::string::npos ) flag = true;

  return flag;
}

bool MuonHLTNtupler::SavedFilterCondition( std::string& filterName )
{
  bool flag = false;

  // -- muon filters
  if( (
        filterName.find("sMu") != std::string::npos ||
        filterName.find("SingleMu") != std::string::npos ||
        filterName.find("TkMu") != std::string::npos
      ) &&
       filterName.find("Tau")      == std::string::npos &&
       filterName.find("EG")       == std::string::npos &&
       filterName.find("MultiFit") == std::string::npos ) flag = true;

  return flag;
}

void MuonHLTNtupler::Fill_GenParticle(const edm::Event &iEvent)
{
  // -- Gen-weight info -- //
  edm::Handle<GenEventInfoProduct> h_genEventInfo;
  iEvent.getByToken(t_genEventInfo_, h_genEventInfo);
  genEventWeight_ = h_genEventInfo->weight();

  // -- Gen-particle info -- //
  edm::Handle<edm::View<reco::GenParticle>> h_genParticle;
  iEvent.getByToken(t_genParticle_, h_genParticle);

  int _nGenParticle = 0;
  for( size_t i=0; i< h_genParticle->size(); ++i)
  {
    const auto &parCand = (*h_genParticle)[i];
    auto genRef = h_genParticle->refAt(i);

    genParticle_ID_[_nGenParticle]     = parCand.pdgId();
    genParticle_status_[_nGenParticle] = parCand.status();
    genParticle_mother_[_nGenParticle] = parCand.mother(0)? parCand.mother(0)->pdgId(): -999;
    
    genParticle_pt_[_nGenParticle]  = parCand.pt();
    genParticle_eta_[_nGenParticle] = parCand.eta();
    genParticle_phi_[_nGenParticle] = parCand.phi();
    genParticle_px_[_nGenParticle]  = parCand.px();
    genParticle_py_[_nGenParticle]  = parCand.py();
    genParticle_pz_[_nGenParticle]  = parCand.pz();
    genParticle_energy_[_nGenParticle] = parCand.energy();
    genParticle_charge_[_nGenParticle] = parCand.charge();
    
    if( parCand.statusFlags().isPrompt() )                genParticle_isPrompt_[_nGenParticle] = 1;
    if( parCand.statusFlags().isTauDecayProduct() )       genParticle_isTauDecayProduct_[_nGenParticle] = 1;
    if( parCand.statusFlags().isPromptTauDecayProduct() ) genParticle_isPromptTauDecayProduct_[_nGenParticle] = 1;
    if( parCand.statusFlags().isDecayedLeptonHadron() )   genParticle_isDecayedLeptonHadron_[_nGenParticle] = 1;
    
    if( parCand.isPromptFinalState() ) genParticle_isPromptFinalState_[_nGenParticle] = 1;
    if( parCand.isDirectPromptTauDecayProductFinalState() ) genParticle_isDirectPromptTauDecayProductFinalState_[_nGenParticle] = 1;
    if( parCand.isHardProcess() ) genParticle_isHardProcess_[_nGenParticle] = 1;
    if( parCand.isLastCopy() ) genParticle_isLastCopy_[_nGenParticle] = 1;
    if( parCand.isLastCopyBeforeFSR() ) genParticle_isLastCopyBeforeFSR_[_nGenParticle] = 1;
    
    if( parCand.isPromptDecayed() )           genParticle_isPromptDecayed_[_nGenParticle] = 1;
    if( parCand.fromHardProcessBeforeFSR() )  genParticle_fromHardProcessBeforeFSR_[_nGenParticle] = 1;
    if( parCand.fromHardProcessDecayed() )    genParticle_fromHardProcessDecayed_[_nGenParticle] = 1;
    if( parCand.fromHardProcessFinalState() ) genParticle_fromHardProcessFinalState_[_nGenParticle] = 1;
    // if( parCand.isMostlyLikePythia6Status3() ) this->genParticle_isMostlyLikePythia6Status3[_nGenParticle] = 1;
    
    _nGenParticle++;
    
  }
  nGenParticle_ = _nGenParticle;
}

/////
void MuonHLTNtupler::Fill_Muon2(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //auto const prop = propSetup_.init(iSetup);

  //theMuonCollectionToken = consumes<RecoChargedCandidateCollection>(theMuonCollectionLabel);
  ///// Muons
  //Handle<RecoChargedCandidateCollection> muons;
  //iEvent.getByToken(theMuonCollectionToken, muons)
  //if( iEvent.getByToken(theMuonCollectionToken, muons) ){ 

  using namespace trigger;
  
  edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;
  if ( iEvent.getByToken(theMuonCollectionToken, PrevFilterOutput) ){

    std::vector<reco::RecoChargedCandidateRef> muons;
    PrevFilterOutput->getObjects(TriggerMuon, muons);
  
    int _nMuon = 0;
    //for(std::vector<reco::RecoChargedCandidateRef>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu)
    //  {
    
    for (unsigned int iMu = 0; iMu < muons.size(); iMu++)
      {

	reco::RecoChargedCandidateRef mu = muons[iMu]; 
	
	muonCand_pt_[_nMuon]  = mu->pt();
	muonCand_eta_[_nMuon] = mu->eta();
	muonCand_phi_[_nMuon] = mu->phi();
	muonCand_px_[_nMuon]  = mu->px();
	muonCand_py_[_nMuon]  = mu->py();
	muonCand_pz_[_nMuon]  = mu->pz();
	muonCand_charge_[_nMuon] = mu->charge();
	
	if( mu->isGlobalMuon() ) muonCand_isGLB_[_nMuon] = 1;
	if( mu->isStandAloneMuon() ) muonCand_isSTA_[_nMuon] = 1;
	if( mu->isTrackerMuon() ) muonCand_isTRK_[_nMuon] = 1;
	if( mu->isMuon() ) muonCand_isMuon_[_nMuon] = 1;
	
	const reco::Track* innerTrk = mu->bestTrack();

        muonCand_vx_[_nMuon] = innerTrk->vx();
        muonCand_vy_[_nMuon] = innerTrk->vy();
	muonCand_vz_[_nMuon] = innerTrk->vz();
	  
	muonCand_dxy_bs_[_nMuon] = innerTrk->dxy(bs->position());
	muonCand_dxyError_bs_[_nMuon] = innerTrk->dxyError(*bs);
	muonCand_dz_bs_[_nMuon] = innerTrk->dz(bs->position());
	muonCand_dzError_[_nMuon] = innerTrk->dzError();
	if (innerTrk->dxyError(*bs) > 0.) {
	  muonCand_IPSig_[_nMuon] = abs(innerTrk->dxy(bs->position()) / innerTrk->dxyError(*bs));
	}
	muonCand_inner_trkChi2_[_nMuon]             = innerTrk->normalizedChi2();
	muonCand_inner_validFraction_[_nMuon]       = innerTrk->validFraction();
	muonCand_inner_trackerLayers_[_nMuon]       = innerTrk->hitPattern().trackerLayersWithMeasurement();
	muonCand_inner_trackerHits_[_nMuon]         = innerTrk->hitPattern().numberOfValidTrackerHits();
	muonCand_inner_lostTrackerHits_[_nMuon]     = innerTrk->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostTrackerHitsIn_[_nMuon]   = innerTrk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
	muonCand_inner_lostTrackerHitsOut_[_nMuon]  = innerTrk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS);
	muonCand_inner_lostPixelHits_[_nMuon]       = innerTrk->hitPattern().numberOfLostPixelHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostPixelBarrelHits_[_nMuon] = innerTrk->hitPattern().numberOfLostPixelBarrelHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostPixelEndcapHits_[_nMuon] = innerTrk->hitPattern().numberOfLostPixelEndcapHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripHits_[_nMuon]       = innerTrk->hitPattern().numberOfLostStripHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripTIBHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTIBHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripTIDHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTIDHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripTOBHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTOBHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripTECHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTECHits(HitPattern::TRACK_HITS);
	muonCand_inner_pixelLayers_[_nMuon]         = innerTrk->hitPattern().pixelLayersWithMeasurement();
	muonCand_inner_pixelHits_[_nMuon]           = innerTrk->hitPattern().numberOfValidPixelHits();
	
	muonCand_global_muonHits_[_nMuon]           = innerTrk->hitPattern().numberOfValidMuonHits();
	muonCand_global_trkChi2_[_nMuon]            = innerTrk->normalizedChi2();
	muonCand_global_trackerLayers_[_nMuon]      = innerTrk->hitPattern().trackerLayersWithMeasurement();
	muonCand_global_trackerHits_[_nMuon]        = innerTrk->hitPattern().numberOfValidTrackerHits();
	
	_nMuon++;
      }
    
    nMuonCand_ = _nMuon;
  }
}
/////
void MuonHLTNtupler::Fill_Track(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  auto const prop = propSetup_.init(iSetup);

  //theMuonCollectionToken = consumes<RecoChargedCandidateCollection>(theMuonCollectionLabel);
  ///// Muons
  //Handle<RecoChargedCandidateCollection> muons;
  // iEvent.getByToken(theMuonCollectionToken, muons);
  //if( iEvent.getByToken(theMuonCollectionToken, muons) ){

  using namespace trigger;
  
  edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;
  if ( iEvent.getByToken(theMuonCollectionToken, PrevFilterOutput) ){

    std::vector<reco::RecoChargedCandidateRef> muons;
    PrevFilterOutput->getObjects(TriggerMuon, muons);
    
    //// Tracks
    Handle<TrackCollection> tracksH;
    //iEvent.getByToken(theTrackCollectionToken, tracksH);
    if ( iEvent.getByToken(theTrackCollectionToken, tracksH) ){
      const TrackCollection tracks = *(tracksH.product());
      
      int _nTrack = 0;  
      unsigned int nMuons = muons.size();
      for (unsigned int iMu = 0; iMu < nMuons; ++iMu) {  
	
	//RecoChargedCandidateRef candref(muons, iMu);
	RecoChargedCandidateRef candref = muons[iMu];
	TrackRef muRef = candref->track();
	const Track& muon = *muRef;
	
	reco::isodeposit::Direction muonDir(muon.eta(), muon.phi());
	
	double vtx_z = muon.vz();
	reco::TrackBase::Point beamPoint(0, 0, 0);
	beamPoint = bs->position();
	
	muonisolation::Range<float> zRange(vtx_z - theDiff_z, vtx_z + theDiff_z);
	muonisolation::Range<float> rRange(0, theDiff_r);    
	for (auto const& tk : tracks) {
	  
	  /// Piece of code from RecoMuon/MuonIsolation/plugin/TrackSelector.cc
	  float tZ = tk.vz();
	  float tPt = tk.pt();
	  float tD0Cor = fabs(tk.dxy(beamPoint));
	  float tChi2Ndof = tk.normalizedChi2();
	  
	  if (!zRange.inside(tZ))
	    continue;
	  if (tPt < thePt_Min)
	    continue;
	  if (!rRange.inside(tD0Cor))
	    continue;
	  if (tChi2Ndof > theChi2Ndof_Max)
	    continue;
	  
	  float tEta = tk.eta();
	  float tPhi = tk.phi();
	  if (muonDir.deltaR2(reco::isodeposit::Direction(tEta, tPhi)) > theDR_Max)
	    continue;

	  if (muonDir.deltaR2(reco::isodeposit::Direction(tEta, tPhi)) < theDR_Veto)
            continue;
	  
	  if (theNHits_Min > 0) {
	    unsigned int tHits = tk.numberOfValidHits();
	    if (tHits < theNHits_Min)
	      continue;
	  }
	  
	  if (theChi2Prob_Min > 0) {
	    float tChi2Prob = ChiSquaredProbability(tk.chi2(), tk.ndof());
	    if (tChi2Prob < theChi2Prob_Min)
	      continue;
	  }
	  
	  ////---------
	  
	  track_pt_[_nTrack] = tk.pt();
	  track_eta_[_nTrack] = tk.eta();
	  track_phi_[_nTrack] = tk.phi();
	  track_charge_[_nTrack] = tk.charge();
	  track_px_[_nTrack] = tk.px();
	  track_py_[_nTrack] = tk.py();
	  track_pz_[_nTrack] = tk.pz();
	  track_vx_[_nTrack] = tk.vx();
	  track_vy_[_nTrack] = tk.vy();
	  track_vz_[_nTrack] = tk.vz();
	  track_dxy_bs_[_nTrack] = tk.dxy(bs->position());
	  track_dxyError_bs_[_nTrack] = tk.dxyError(*bs);
	  track_dz_bs_[_nTrack] = tk.dz(bs->position());
	  track_dzError_[_nTrack] = tk.dzError();
	  track_trkChi2_[_nTrack] = tk.normalizedChi2();
	  track_trackerLayers_[_nTrack] = tk.hitPattern().trackerLayersWithMeasurement();
	  track_trackerHits_[_nTrack] = tk.hitPattern().numberOfValidTrackerHits();
	  track_lostTrackerHits_[_nTrack] = tk.hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS);
	  track_lostTrackerHitsIn_[_nTrack] = tk.hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
	  track_lostTrackerHitsOut_[_nTrack] = tk.hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS);
	  track_lostPixelHits_[_nTrack] = tk.hitPattern().numberOfLostPixelHits(HitPattern::TRACK_HITS);
	  track_lostPixelBarrelHits_[_nTrack] = tk.hitPattern().numberOfLostPixelBarrelHits(HitPattern::TRACK_HITS);
	  track_lostPixelEndcapHits_[_nTrack] = tk.hitPattern().numberOfLostPixelEndcapHits(HitPattern::TRACK_HITS);
	  track_lostStripHits_[_nTrack] = tk.hitPattern().numberOfLostStripHits(HitPattern::TRACK_HITS);
	  track_lostStripTIBHits_[_nTrack] = tk.hitPattern().numberOfLostStripTIBHits(HitPattern::TRACK_HITS);
	  track_lostStripTIDHits_[_nTrack] = tk.hitPattern().numberOfLostStripTIDHits(HitPattern::TRACK_HITS);
	  track_lostStripTOBHits_[_nTrack] = tk.hitPattern().numberOfLostStripTOBHits(HitPattern::TRACK_HITS);
	  track_lostStripTECHits_[_nTrack] = tk.hitPattern().numberOfLostStripTECHits(HitPattern::TRACK_HITS);
	  track_pixelLayers_[_nTrack] = tk.hitPattern().pixelLayersWithMeasurement();
	  track_pixelHits_[_nTrack] = tk.hitPattern().numberOfValidPixelHits();
	  track_muonIdx_[_nTrack] = iMu;
	  _nTrack++;
	}
	
      }
      nTrack_ = _nTrack;  
    }
  }
}

bool MuonHLTNtupler::computedRVeto(RecoChargedCandidateRef candRef, reco::PFClusterRef pfclu, double drMAX, double drVeto2_) {
  float dR2 = deltaR2(candRef->eta(), candRef->phi(), pfclu->eta(), pfclu->phi());
  if (dR2 > (drMAX * drMAX) || dR2 < drVeto2_)
    return false;
  else
    return true;
}


void MuonHLTNtupler::Fill_ECAL(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  auto const prop = propSetup_.init(iSetup);

  //theMuonCollectionToken = consumes<RecoChargedCandidateCollection>(theMuonCollectionLabel);
  ///// Muons
  //Handle<RecoChargedCandidateCollection> muons;
  //iEvent.getByToken(theMuonCollectionToken, muons);
  //if ( iEvent.getByToken(theMuonCollectionToken, muons) ){

  using namespace trigger;
  
  edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;
  if ( iEvent.getByToken(theMuonCollectionToken, PrevFilterOutput) ){
    
    std::vector<reco::RecoChargedCandidateRef> muons;
    PrevFilterOutput->getObjects(TriggerMuon, muons);
    
    if (doRhoCorrection_ECAL_) {
      if (absEtaLowEdges_ECAL_.size() != effectiveAreas_ECAL_.size())
	throw cms::Exception("IncompatibleVects") << "absEtaLowEdges and effectiveAreas should be of the same size. \n";
      
      if (absEtaLowEdges_ECAL_.at(0) != 0.0)
	throw cms::Exception("IncompleteCoverage") << "absEtaLowEdges should start from 0. \n";
      
      for (unsigned int aIt = 0; aIt < absEtaLowEdges_ECAL_.size() - 1; aIt++) {
	if (!(absEtaLowEdges_ECAL_.at(aIt) < absEtaLowEdges_ECAL_.at(aIt + 1)))
	  throw cms::Exception("ImproperBinning") << "absEtaLowEdges entries should be in increasing order. \n";
      }
    }

    edm::Handle<double> rhoHandle_ECAL;
    if (doRhoCorrection_ECAL_) {
      iEvent.getByToken(t_rho_ECAL_, rhoHandle_ECAL);
      rho_ECAL = *(rhoHandle_ECAL.product());
    }
    if (rho_ECAL > rhoMax_ECAL_)
      rho_ECAL = rhoMax_ECAL_;
    
    rho_ECAL = rho_ECAL * rhoScale_ECAL_;
    
    edm::Handle<reco::PFClusterCollection> clusterHandle;
    //iEvent.getByToken(pfClusterProducer_, clusterHandle);
    if ( iEvent.getByToken(pfClusterProducer_, clusterHandle) ){
      
      int _nECAL = 0;
      int _nHits = 0;
      
      unsigned int nMuons = muons.size();
      for (unsigned int iMu = 0; iMu < nMuons; ++iMu) {
	
	//RecoChargedCandidateRef candref(muons, iMu);
	RecoChargedCandidateRef candref = muons[iMu];
	TrackRef muRef = candref->track();
	const Track& muon = *muRef;
	
	reco::isodeposit::Direction muonDir(muon.eta(), muon.phi());
	
	/// Piece of code from cmssw/RecoEgamma/EgammaIsolationAlgos/src/EcalPFClusterIsolation.cc
	double drVeto2_ = -1;
	float etaStrip = -1;
	
	if (std::abs(muon.eta()) < 1.479) {
	  drVeto2_ = drVetoBarrel_ECAL_ * drVetoBarrel_ECAL_;
	  etaStrip = etaStripBarrel_ECAL_;
	} else {
	  drVeto2_ = drVetoEndcap_ECAL_ * drVetoEndcap_ECAL_;
	  etaStrip = etaStripEndcap_ECAL_;
	}
	
	for (size_t i = 0; i < clusterHandle->size(); i++) {
	  reco::PFClusterRef pfcluRef(clusterHandle, i);
	  const reco::PFCluster& pfclu = *pfcluRef;
	  
	  if (std::abs(muon.eta()) < 1.479) {
	    if (std::abs(pfclu.pt()) < energyBarrel_ECAL_)
	      continue;
	  } else {
	    if (std::abs(pfclu.energy()) < energyEndcap_ECAL_)
	      continue;
	  }
	  
	  float dEta = std::abs(muon.eta() - pfclu.eta());
	  if (dEta < etaStrip)
	    continue;
	  if (not computedRVeto(candref, pfcluRef, drMax_ECAL_, drVeto2_))
	    continue;
	  
	  ecal_et_[_nECAL] = pfclu.energy();
	  ecal_pt_[_nECAL] = pfclu.pt();
	  ecal_eta_[_nECAL] = pfclu.eta();
	  ecal_phi_[_nECAL] = pfclu.phi();
	  ecal_charge_[_nECAL] = pfclu.charge();      
	  ecal_px_[_nECAL] = pfclu.x();
	  ecal_py_[_nECAL] = pfclu.y();
	  ecal_pz_[_nECAL] = pfclu.z();
	  ecal_vx_[_nECAL] = pfclu.vx();
	  ecal_vy_[_nECAL] = pfclu.vy();
	  ecal_vz_[_nECAL] = pfclu.vz();      
	  ecal_time_[_nECAL] = (double)pfclu.time();
	  ecal_depth_[_nECAL] = pfclu.depth();
	  ecal_rho_[_nECAL] = rho_ECAL;
	  ecal_muonIdx_[_nECAL] = iMu;
	  
	  const std::vector<reco::PFRecHitFraction>& ecal_recHitFractions = pfclu.recHitFractions();
	  
	  unsigned int nHits = ecal_recHitFractions.size();
	  for (unsigned int iH = 0; iH < nHits; ++iH) {
	    reco::PFRecHitFraction ecal_hit = ecal_recHitFractions[iH];
	    const PFRecHitRef& ecal_hitRef = ecal_hit.recHitRef();
	    const reco::PFRecHit ecal_rechit = *ecal_hitRef;
	    
	    ecal_hit_energy_[_nHits] = ecal_rechit.energy();
	    ecal_hit_depth_[_nHits] = ecal_rechit.depth();
	    ecal_hit_time_[_nHits] = ecal_rechit.time();
	    ecal_hit_pt2_[_nHits] = ecal_rechit.pt2();
	    ecal_hit_x_[_nHits] = ecal_rechit.position().x();
	    ecal_hit_y_[_nHits] = ecal_rechit.position().y();
	    ecal_hit_z_[_nHits] = ecal_rechit.position().z();
	    
	    ecal_hit_eta_[_nHits] = ecal_rechit.positionREP().eta();
	    ecal_hit_phi_[_nHits] = ecal_rechit.positionREP().phi();
	    
	    ecal_hit_fraction_[_nHits] = ecal_hit.fraction();
	    ecal_hit_idx_[_nHits] = _nECAL;
	    
	    _nHits++;
	  }
	  
	  ecal_nHits_[_nECAL] = nHits;
	  
	  _nECAL++;
	}
      }
      nECAL_ = _nECAL;
      nECALHits_ = _nHits;
    }
  }
}


void MuonHLTNtupler::Fill_HCAL(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  auto const prop = propSetup_.init(iSetup);

  //theMuonCollectionToken = consumes<RecoChargedCandidateCollection>(theMuonCollectionLabel);
  ///// Muons
  //Handle<RecoChargedCandidateCollection> muons;
  //
  //iEvent.getByToken(theMuonCollectionToken, muons);   
  //if ( iEvent.getByToken(theMuonCollectionToken, muons) ){

  using namespace trigger;
  
  edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;  
  if ( iEvent.getByToken(theMuonCollectionToken, PrevFilterOutput) ){

    std::vector<reco::RecoChargedCandidateRef> muons;
    PrevFilterOutput->getObjects(TriggerMuon, muons);
    
    if (doRhoCorrection_HCAL_) {
      if (absEtaLowEdges_HCAL_.size() != effectiveAreas_HCAL_.size())
	throw cms::Exception("IncompatibleVects") << "absEtaLowEdges and effectiveAreas should be of the same size. \n";
      
      if (absEtaLowEdges_HCAL_.at(0) != 0.0)
	throw cms::Exception("IncompleteCoverage") << "absEtaLowEdges should start from 0. \n";
      
      for (unsigned int aIt = 0; aIt < absEtaLowEdges_HCAL_.size() - 1; aIt++) {
	if (!(absEtaLowEdges_HCAL_.at(aIt) < absEtaLowEdges_HCAL_.at(aIt + 1)))
	  throw cms::Exception("ImproperBinning") << "absEtaLowEdges entries should be in increasing order. \n";
      }
    }
    edm::Handle<double> rhoHandle_HCAL;
    if (doRhoCorrection_HCAL_) {
      iEvent.getByToken(t_rho_HCAL_, rhoHandle_HCAL);
      rho_HCAL = *(rhoHandle_HCAL.product());
    }
    
    if (rho_HCAL > rhoMax_HCAL_)
      rho_HCAL = rhoMax_HCAL_;
    
    rho_HCAL = rho_HCAL * rhoScale_HCAL_;
    
    std::vector<edm::Handle<reco::PFClusterCollection>> clusterHandles;
    edm::Handle<reco::PFClusterCollection> clusterHcalHandle;
    edm::Handle<reco::PFClusterCollection> clusterHfemHandle;
    edm::Handle<reco::PFClusterCollection> clusterHfhadHandle;  
  
    //iEvent.getByToken(pfClusterProducerHCAL_, clusterHcalHandle);
    if ( iEvent.getByToken(pfClusterProducerHCAL_, clusterHcalHandle) ){
      
      clusterHandles.push_back(clusterHcalHandle);
      if (useHF_) {
	iEvent.getByToken(pfClusterProducerHFEM_, clusterHfemHandle);
	clusterHandles.push_back(clusterHfemHandle);
	iEvent.getByToken(pfClusterProducerHFHAD_, clusterHfhadHandle);
	clusterHandles.push_back(clusterHfhadHandle);
      }
      
      int _nHCAL = 0;
      int _nHits = 0;
      
      unsigned int nMuons = muons.size();
      for (unsigned int iMu = 0; iMu < nMuons; ++iMu) {
	
	RecoChargedCandidateRef candref = muons[iMu];
	TrackRef muRef = candref->track();
	const Track& muon = *muRef;
	
	reco::isodeposit::Direction muonDir(muon.eta(), muon.phi());
	
	double candAbsEta = std::abs(muon.eta());
	float etaStrip = 0;
	float dRVeto = 0;
	if (candAbsEta < 1.479) {
	  dRVeto = drVetoBarrel_HCAL_;
	  etaStrip = etaStripBarrel_HCAL_;
	} else {
	  dRVeto = drVetoEndcap_HCAL_;
	  etaStrip = etaStripEndcap_HCAL_;
	}
	
	for (unsigned int nHandle = 0; nHandle < clusterHandles.size(); nHandle++) {
	  for (unsigned i = 0; i < clusterHandles[nHandle]->size(); i++) {
	    const reco::PFClusterRef pfcluRef(clusterHandles[nHandle], i);
	    const reco::PFCluster& pfclu = *pfcluRef;
	    
	    if (candAbsEta < 1.479) {
	      if (std::abs(pfclu.pt()) < energyBarrel_HCAL_)
		continue;
	    } else {
	      if (std::abs(pfclu.energy()) < energyEndcap_HCAL_)
		continue;
	    }
	    
	    float dEta = std::abs(muon.eta() - pfclu.eta());
	    if (dEta < etaStrip)
	      continue;
	    
	    float dR2 = deltaR2(muon.eta(), muon.phi(), pfclu.eta(), pfclu.phi());
	    if (dR2 > (drMax_HCAL_ * drMax_HCAL_) || dR2 < (dRVeto * dRVeto))
	      continue;
	    
	    
	    hcal_et_[_nHCAL] = pfclu.energy();
	    hcal_pt_[_nHCAL] = pfclu.pt();
	    hcal_eta_[_nHCAL] = pfclu.eta();
	    hcal_phi_[_nHCAL] = pfclu.phi();
	    hcal_charge_[_nHCAL] = pfclu.charge();
	    hcal_px_[_nHCAL] = pfclu.x();
	    hcal_py_[_nHCAL] = pfclu.y();
	    hcal_pz_[_nHCAL] = pfclu.z();
	    hcal_vx_[_nHCAL] = pfclu.vx();
	    hcal_vy_[_nHCAL] = pfclu.vy();
	    hcal_vz_[_nHCAL] = pfclu.vz();
	    hcal_time_[_nHCAL] = (double)pfclu.time();
	    hcal_depth_[_nHCAL] = pfclu.depth();
	    hcal_rho_[_nHCAL] = rho_HCAL;
	    hcal_muonIdx_[_nHCAL] = iMu;
	    
	    const std::vector<reco::PFRecHitFraction>& hcal_recHitFractions = pfclu.recHitFractions();
	    
	    unsigned int nHits = hcal_recHitFractions.size();
	    for (unsigned int iH = 0; iH < nHits; ++iH) {
	      reco::PFRecHitFraction hcal_hit = hcal_recHitFractions[iH];
	      const PFRecHitRef& hcal_hitRef = hcal_hit.recHitRef();
	      const reco::PFRecHit hcal_rechit = *hcal_hitRef;
	      
	      hcal_hit_energy_[_nHits] = hcal_rechit.energy();
	      hcal_hit_depth_[_nHits] = hcal_rechit.depth();
	      hcal_hit_time_[_nHits] = hcal_rechit.time();
	      hcal_hit_pt2_[_nHits] = hcal_rechit.pt2();
	      hcal_hit_x_[_nHits] = hcal_rechit.position().x();
	      hcal_hit_y_[_nHits] = hcal_rechit.position().y();
	      hcal_hit_z_[_nHits] = hcal_rechit.position().z();
	      
	      hcal_hit_eta_[_nHits] = hcal_rechit.positionREP().eta();
	      hcal_hit_phi_[_nHits] = hcal_rechit.positionREP().phi();
	      
	      hcal_hit_fraction_[_nHits] = hcal_hit.fraction();
	      hcal_hit_idx_[_nHits] = _nHCAL;
	      
	      _nHits++;
	    }
	    
	    hcal_nHits_[_nHCAL] = nHits;	
	    _nHCAL++;
	    
	  }
	}
      }
      nHCAL_ = _nHCAL;
      nHCALHits_ = _nHits;
    }
  }
}

// -- reference: https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonReco/src/MuonSelectors.cc#L910-L938
bool MuonHLTNtupler::isNewHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx){
  if(!muon.isGlobalMuon()) return false;

  bool muValHits = ( muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0 ||
                     muon.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits()>0 );

  bool muMatchedSt = muon.numberOfMatchedStations()>1;
  if(!muMatchedSt) {
    if( muon.isTrackerMuon() && muon.numberOfMatchedStations()==1 ) {
      if( muon.expectedNnumberOfMatchedStations()<2 ||
          !(muon.stationMask()==1 || muon.stationMask()==16) ||
          muon.numberOfMatchedRPCLayers()>2
        )
        muMatchedSt = true;
    }
  }

  bool muID = muValHits && muMatchedSt;

  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

  bool momQuality = muon.tunePMuonBestTrack()->ptError()/muon.tunePMuonBestTrack()->pt() < 0.3;

  bool ip = fabs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.innerTrack()->dz(vtx.position())) < 0.5;

  return muID && hits && momQuality && ip;
}

void MuonHLTNtupler::endJob() {
  //for( int i=0; i<4; ++i ) {
  // for( int i=0; i<1; ++i ) {
  //   delete mvaHltIter2IterL3MuonPixelSeeds_.at(i).first;
  //   delete mvaHltIter2IterL3MuonPixelSeeds_.at(i).second;
  //   delete mvaHltIter2IterL3FromL1MuonPixelSeeds_.at(i).first;
  //   delete mvaHltIter2IterL3FromL1MuonPixelSeeds_.at(i).second;
  // }

  // for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
  //   delete trkTemplates_.at(i);
  //   delete tpTemplates_.at(i);
  // }
}

// void MuonHLTNtupler::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}
// void MuonHLTNtupler::endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}

DEFINE_FWK_MODULE(MuonHLTNtupler);
