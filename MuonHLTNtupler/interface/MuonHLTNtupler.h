// -- ntuple maker for Muon HLT study
// -- author: Kyeongpil Lee (Seoul National University, kplee@cern.ch)

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
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
// #include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
// #include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
// #include "DataFormats/L1TrackTrigger/interface/TTStub.h"
// #include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

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
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

//--- for SimHit association
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
////////////////////////////
// DETECTOR GEOMETRY HEADERS
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

// #include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
// #include "Geometry/CommonTopologies/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
// #include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
// #include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
// #include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
// #include "DataFormats/L1TCorrelator/interface/TkMuon.h"
// #include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"
// #include "DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h"

#include "RecoMuon/TrackerSeedGenerator/interface/SeedMvaEstimator.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

// #include "MuonHLTTool/MuonHLTNtupler/interface/MuonHLTobjCorrelator.h"

#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractorFactory.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"


#include "TTree.h"
#include "TString.h"

using namespace std;
using namespace reco;
using namespace edm;

class MuonHLTNtupler : public edm::one::EDAnalyzer<>
{
public:
  explicit MuonHLTNtupler(const edm::ParameterSet &iConfig);
  virtual ~MuonHLTNtupler() {};

  virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  virtual void beginJob();
  virtual void endJob();

  // virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);
  // virtual void endRun(const edm::Run &iRun, const edm::EventSetup &iSetup);

private:
  void Init();
  void Make_Branch();

  //// Fill functions -------------------------------------------------------
  void Fill_Muon(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  void Fill_HLT(const edm::Event &iEvent, bool isMYHLT);

  bool SavedTriggerCondition( std::string& pathName );
  bool SavedFilterCondition( std::string& filterName );

  void Fill_GenParticle(const edm::Event &iEvent);
  
  void Fill_Muon2(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  void Fill_Track(const edm::Event &iEvent, const edm::EventSetup &iSetup);

  bool computedRVeto(RecoChargedCandidateRef candRef, reco::PFClusterRef pfclu, double drMAX, double drVeto2_);

  void Fill_ECAL(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  void Fill_HCAL(const edm::Event &iEvent, const edm::EventSetup &iSetup);

  bool isNewHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx);

  //const edm::InputTag theMuonCollectionLabel;

  // HCAL

  const bool useHF_;
  const double drMax_HCAL_;
  const double drVetoBarrel_HCAL_;
  const double drVetoEndcap_HCAL_;
  const double etaStripBarrel_HCAL_;
  const double etaStripEndcap_HCAL_;
  const double energyBarrel_HCAL_;
  const double energyEndcap_HCAL_;

  const bool doRhoCorrection_HCAL_;
  const double rhoMax_HCAL_;
  const double rhoScale_HCAL_;
  const std::vector<double> effectiveAreas_HCAL_;
  const std::vector<double> absEtaLowEdges_HCAL_;

  // ECAL
  
  const double drMax_ECAL_;
  const double drVetoBarrel_ECAL_;
  const double drVetoEndcap_ECAL_;
  const double etaStripBarrel_ECAL_;
  const double etaStripEndcap_ECAL_;
  const double energyBarrel_ECAL_;
  const double energyEndcap_ECAL_;

  const bool doRhoCorrection_ECAL_;
  const double rhoMax_ECAL_;
  const double rhoScale_ECAL_;
  const std::vector<double> effectiveAreas_ECAL_;
  const std::vector<double> absEtaLowEdges_ECAL_;

  // General CALO
  
  const bool useEt_;                          // use E or Et in relative isolation cuts

  // Tracker isolation
  
  const std::string theDepositLabel;                                      //! name for deposit
  const double theDiff_r;                                                 //! transverse distance to vertex 
  const double theDiff_z;                                                 //! z distance to vertex
  const double theDR_Max;                                                 //! Maximum cone angle for deposits 
  const double theDR_Veto;                                                //! Veto cone angle
  const unsigned int theNHits_Min;                                        //! trk.numberOfValidHits >= theNHits_Min
  const double theChi2Ndof_Max;                                           //! trk.normalizedChi2 < theChi2Ndof_Max
  const double theChi2Prob_Min;  //! ChiSquaredProbability(trk.chi2,trk.ndof) > theChi2Prob_Min
  const double thePt_Min;        //! min track pt to include into iso deposit     

  // Muon propagator

  const PropagateToMuonSetup propSetup_;

  
  //// General parameters ------------------------------------------------------------------------------------------------
  //edm::ParameterSet theConfig;
  //// Inputs ----------------------------------------------------------------------------------------------------------

  const edm::EDGetTokenT< reco::BeamSpot >                         t_beamSpot_;
  const edm::EDGetTokenT< std::vector<reco::Muon> >                t_offlineMuon_;
  const edm::EDGetTokenT< reco::VertexCollection >                 t_offlineVertex_;
  const edm::EDGetTokenT< double >                                 t_rho_ECAL_;
  const edm::EDGetTokenT< double >                                 t_rho_HCAL_;

  const edm::EDGetTokenT< LumiScalersCollection >                  t_lumiScaler_;
  const edm::EDGetTokenT< LumiScalersCollection >                  t_offlineLumiScaler_;
  const edm::EDGetTokenT< std::vector<PileupSummaryInfo> >         t_PUSummaryInfo_;
  const edm::EDGetTokenT< GenEventInfoProduct >                    t_genEventInfo_;
  const edm::EDGetTokenT< edm::View<reco::GenParticle> >           t_genParticle_;

  const edm::EDGetTokenT< edm::TriggerResults >                    t_triggerResults_;
  const edm::EDGetTokenT< trigger::TriggerEvent >                  t_triggerEvent_;
  const edm::EDGetTokenT< edm::TriggerResults >                    t_myTriggerResults_;
  const edm::EDGetTokenT< trigger::TriggerEvent >                  t_myTriggerEvent_;
  
  // HCAL PF Candidates
  const edm::EDGetTokenT<reco::PFClusterCollection> pfClusterProducerHCAL_;
  const edm::EDGetTokenT<double> rhoProducer_HCAL_;
  const edm::EDGetTokenT<reco::PFClusterCollection> pfClusterProducerHFEM_;
  const edm::EDGetTokenT<reco::PFClusterCollection> pfClusterProducerHFHAD_;
  
  // ECAL PF Candidates
  const edm::EDGetTokenT<reco::PFClusterCollection> pfClusterProducer_;
  const edm::EDGetTokenT<double> rhoProducer_ECAL_;

  // Track
  const edm::EDGetTokenT<reco::TrackCollection> theTrackCollectionToken;
  
  // Muon track Collection Label  
  //const edm::EDGetTokenT<reco::RecoChargedCandidateCollection> theMuonCollectionToken;
  const edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> theMuonCollectionToken;
  
  // ECAL Isolation
  
  std::vector<double> energyLowEdges_;  // lower bin edges for energy-dependent cuts
  bool lessThan_;                       // the cut is "<" or ">" ?
  bool useAbs_;                         // use the standard abs of the variable (before any rho corr)
  
  std::vector<double> thrRegularEB_ECAL_;    // threshold for regular cut (x < thr) - ECAL barrel
  std::vector<double> thrRegularEE_ECAL_;    // threshold for regular cut (x < thr) - ECAL endcap
  std::vector<double> thrOverEEB_ECAL_;      // threshold for x/E < thr cut (isolations) - ECAL barrel
  std::vector<double> thrOverEEE_ECAL_;      // threshold for x/E < thr cut (isolations) - ECAL endcap
  std::vector<double> thrOverE2EB_ECAL_;     // threshold for x/E^2 < thr cut (isolations) - ECAL barrel
  std::vector<double> thrOverE2EE_ECAL_;     // threshold for x/E^2 < thr cut (isolations) - ECAL endcap

  std::vector<double> thrRegularEB_HCAL_;
  std::vector<double> thrRegularEE_HCAL_;
  std::vector<double> thrOverEEB_HCAL_;
  std::vector<double> thrOverEEE_HCAL_;
  std::vector<double> thrOverE2EB_HCAL_; 
  std::vector<double> thrOverE2EE_HCAL_;
  
  int ncandcut_;                        // number of candidates required

  edm::InputTag l1EGTag_;

  edm::InputTag rhoTag_;  // input tag identifying rho producer
  edm::EDGetTokenT<double> rhoToken_;
  
  
  // Tracker parameters 
  
  //! pt cut to consider track in sumPt after extracting iso deposit
  //! better split this off into a filter
  double theTrackPt_Min;

  //! max number of tracks to allow in the sum
  //! count <= maxN
  int theMaxNTracks;

  //! apply or not the maxN cut on top of the sumPt (or nominall eff) < cuts
  bool theApplyCutsORmaxNTracks;

  //// Ntuples branches --------------------------------------------------------------------------------------------------

  const reco::BeamSpot* bs;

  TTree *ntuple_;
  static const int arrSize_ = 5000;

  // -- general event information
  bool isRealData_;
  int runNum_;
  int lumiBlockNum_;
  unsigned long long eventNum_;

  double bs_x0_;
  double bs_y0_;
  double bs_z0_;
  double bs_sigmaZ_;
  double bs_dxdz_;
  double bs_dydz_;
  double bs_x0Error_;
  double bs_y0Error_;
  double bs_z0Error_;
  double bs_sigmaZ0Error_;
  double bs_dxdzError_;
  double bs_dydzError_;

  int nVertex_;

  double vertex_x_[arrSize_];
  double vertex_y_[arrSize_];
  double vertex_z_[arrSize_];
  double vertex_t_[arrSize_];

  double bunchID_;
  double instLumi_;
  double dataPU_;
  double dataPURMS_;
  double bunchLumi_;
  double offlineInstLumi_;
  double offlineDataPU_;
  double offlineDataPURMS_;
  double offlineBunchLumi_;
  int truePU_;
  double genEventWeight_;

  double rho_ECAL;
  double rho_HCAL;

  // -- generator level particles (only MC)
  int nGenParticle_;
  int genParticle_ID_[arrSize_];
  int genParticle_status_[arrSize_];
  int genParticle_mother_[arrSize_];

  double genParticle_pt_[arrSize_];
  double genParticle_eta_[arrSize_];
  double genParticle_phi_[arrSize_];
  double genParticle_px_[arrSize_];
  double genParticle_py_[arrSize_];
  double genParticle_pz_[arrSize_];
  double genParticle_energy_[arrSize_];
  double genParticle_charge_[arrSize_];

  int genParticle_isPrompt_[arrSize_];
  int genParticle_isPromptFinalState_[arrSize_];
  int genParticle_isTauDecayProduct_[arrSize_];
  int genParticle_isPromptTauDecayProduct_[arrSize_];
  int genParticle_isDirectPromptTauDecayProductFinalState_[arrSize_];
  int genParticle_isHardProcess_[arrSize_];
  int genParticle_isLastCopy_[arrSize_];
  int genParticle_isLastCopyBeforeFSR_[arrSize_];
  int genParticle_isPromptDecayed_[arrSize_];
  int genParticle_isDecayedLeptonHadron_[arrSize_];
  int genParticle_fromHardProcessBeforeFSR_[arrSize_];
  int genParticle_fromHardProcessDecayed_[arrSize_];
  int genParticle_fromHardProcessFinalState_[arrSize_];
  int genParticle_isMostlyLikePythia6Status3_[arrSize_];

  // -- offline muon
  int nMuon_;

  double muon_pt_[arrSize_];
  double muon_eta_[arrSize_];
  double muon_phi_[arrSize_];
  double muon_px_[arrSize_];
  double muon_py_[arrSize_];
  double muon_pz_[arrSize_];
  double muon_dB_[arrSize_];
  double muon_charge_[arrSize_];
  int muon_isGLB_[arrSize_];
  int muon_isSTA_[arrSize_];
  int muon_isTRK_[arrSize_];
  int muon_isPF_[arrSize_];
  int muon_isTight_[arrSize_];
  int muon_isMedium_[arrSize_];
  int muon_isLoose_[arrSize_];
  int muon_isHighPt_[arrSize_];
  int muon_isHighPtNew_[arrSize_];
  int muon_isSoft_[arrSize_];

  double muon_iso03_sumPt_[arrSize_];
  double muon_iso03_hadEt_[arrSize_];
  double muon_iso03_emEt_[arrSize_];

  double muon_PFIso03_charged_[arrSize_];
  double muon_PFIso03_neutral_[arrSize_];
  double muon_PFIso03_photon_[arrSize_];
  double muon_PFIso03_sumPU_[arrSize_];

  double muon_PFIso04_charged_[arrSize_];
  double muon_PFIso04_neutral_[arrSize_];
  double muon_PFIso04_photon_[arrSize_];
  double muon_PFIso04_sumPU_[arrSize_];

  double muon_PFCluster03_ECAL_[arrSize_];
  double muon_PFCluster03_HCAL_[arrSize_];

  double muon_PFCluster04_ECAL_[arrSize_];
  double muon_PFCluster04_HCAL_[arrSize_];

  double muon_inner_trkChi2_[arrSize_];
  double muon_inner_validFraction_[arrSize_];
  int    muon_inner_trackerLayers_[arrSize_];
  int    muon_inner_trackerHits_[arrSize_];
  int    muon_inner_lostTrackerHits_[arrSize_];
  int    muon_inner_lostTrackerHitsIn_[arrSize_];
  int    muon_inner_lostTrackerHitsOut_[arrSize_];
  int    muon_inner_lostPixelHits_[arrSize_];
  int    muon_inner_lostPixelBarrelHits_[arrSize_];
  int    muon_inner_lostPixelEndcapHits_[arrSize_];
  int    muon_inner_lostStripHits_[arrSize_];
  int    muon_inner_lostStripTIBHits_[arrSize_];
  int    muon_inner_lostStripTIDHits_[arrSize_];
  int    muon_inner_lostStripTOBHits_[arrSize_];
  int    muon_inner_lostStripTECHits_[arrSize_];
  int    muon_inner_pixelLayers_[arrSize_];
  int    muon_inner_pixelHits_[arrSize_];
  int    muon_global_muonHits_[arrSize_];
  double muon_global_trkChi2_[arrSize_];
  int    muon_global_trackerLayers_[arrSize_];
  int    muon_global_trackerHits_[arrSize_];
  double muon_momentumChi2_[arrSize_];
  double muon_positionChi2_[arrSize_];
  double muon_glbKink_[arrSize_];
  double muon_glbTrackProbability_[arrSize_];
  double muon_globalDeltaEtaPhi_[arrSize_];
  double muon_localDistance_[arrSize_];
  double muon_staRelChi2_[arrSize_];
  int    muon_tightMatch_[arrSize_];
  double muon_trkKink_[arrSize_];
  double muon_trkRelChi2_[arrSize_];
  double muon_segmentCompatibility_[arrSize_];

  double muon_pt_tuneP_[arrSize_];
  double muon_ptError_tuneP_[arrSize_];

  double muon_dxyVTX_best_[arrSize_];
  double muon_dzVTX_best_[arrSize_];

  int muon_nMatchedStation_[arrSize_];
  int muon_nMatchedRPCLayer_[arrSize_];
  int muon_stationMask_[arrSize_];

  double muon_dxy_bs_[arrSize_];
  double muon_dxyError_bs_[arrSize_];
  double muon_dz_bs_[arrSize_];
  double muon_dzError_[arrSize_];
  double muon_IPSig_[arrSize_];


  int nMuonCand_;

  double muonCand_pt_[arrSize_];
  double muonCand_eta_[arrSize_];
  double muonCand_phi_[arrSize_];
  double muonCand_px_[arrSize_];
  double muonCand_py_[arrSize_];
  double muonCand_pz_[arrSize_];
  double muonCand_vx_[arrSize_];
  double muonCand_vy_[arrSize_];
  double muonCand_vz_[arrSize_];
  double muonCand_dB_[arrSize_];
  double muonCand_charge_[arrSize_];
  int muonCand_isGLB_[arrSize_];
  int muonCand_isSTA_[arrSize_];
  int muonCand_isTRK_[arrSize_];
  int muonCand_isMuon_[arrSize_];

  double muonCand_inner_trkChi2_[arrSize_];
  double muonCand_inner_validFraction_[arrSize_];
  int    muonCand_inner_trackerLayers_[arrSize_];
  int    muonCand_inner_trackerHits_[arrSize_];
  int    muonCand_inner_lostTrackerHits_[arrSize_];
  int    muonCand_inner_lostTrackerHitsIn_[arrSize_];
  int    muonCand_inner_lostTrackerHitsOut_[arrSize_];
  int    muonCand_inner_lostPixelHits_[arrSize_];
  int    muonCand_inner_lostPixelBarrelHits_[arrSize_];
  int    muonCand_inner_lostPixelEndcapHits_[arrSize_];
  int    muonCand_inner_lostStripHits_[arrSize_];
  int    muonCand_inner_lostStripTIBHits_[arrSize_];
  int    muonCand_inner_lostStripTIDHits_[arrSize_];
  int    muonCand_inner_lostStripTOBHits_[arrSize_];
  int    muonCand_inner_lostStripTECHits_[arrSize_];
  int    muonCand_inner_pixelLayers_[arrSize_];
  int    muonCand_inner_pixelHits_[arrSize_];
  int    muonCand_global_muonHits_[arrSize_];
  double muonCand_global_trkChi2_[arrSize_];
  int    muonCand_global_trackerLayers_[arrSize_];
  int    muonCand_global_trackerHits_[arrSize_];

  double muonCand_dxy_bs_[arrSize_];
  double muonCand_dxyError_bs_[arrSize_];
  double muonCand_dz_bs_[arrSize_];
  double muonCand_dzError_[arrSize_];
  double muonCand_IPSig_[arrSize_];
  
  int nTrack_;
  
  double track_pt_[arrSize_];
  double track_eta_[arrSize_];
  double track_phi_[arrSize_];
  double track_charge_[arrSize_];
  double track_px_[arrSize_];
  double track_py_[arrSize_];
  double track_pz_[arrSize_];
  double track_vx_[arrSize_];
  double track_vy_[arrSize_];
  double track_vz_[arrSize_];
  double track_dxy_bs_[arrSize_];
  double track_dxyError_bs_[arrSize_];
  double track_dz_bs_[arrSize_];
  double track_dzError_[arrSize_];
  double track_trkChi2_[arrSize_];
  int    track_trackerLayers_[arrSize_];
  int    track_trackerHits_[arrSize_];
  int    track_lostTrackerHits_[arrSize_];
  int    track_lostTrackerHitsIn_[arrSize_];
  int    track_lostTrackerHitsOut_[arrSize_];
  int    track_lostPixelHits_[arrSize_];
  int    track_lostPixelBarrelHits_[arrSize_];
  int    track_lostPixelEndcapHits_[arrSize_];
  int    track_lostStripHits_[arrSize_];
  int    track_lostStripTIBHits_[arrSize_];
  int    track_lostStripTIDHits_[arrSize_];
  int    track_lostStripTOBHits_[arrSize_];
  int    track_lostStripTECHits_[arrSize_];
  int    track_pixelLayers_[arrSize_];
  int    track_pixelHits_[arrSize_];
  int    track_muonHits_[arrSize_];
  int    track_muonIdx_[arrSize_];

  int nECAL_;

  double ecal_et_[arrSize_];
  double ecal_pt_[arrSize_];
  double ecal_eta_[arrSize_];
  double ecal_phi_[arrSize_];
  double ecal_charge_[arrSize_];
  double ecal_px_[arrSize_];
  double ecal_py_[arrSize_];
  double ecal_pz_[arrSize_];
  double ecal_vx_[arrSize_];
  double ecal_vy_[arrSize_];
  double ecal_vz_[arrSize_];
  double ecal_time_[arrSize_];
  double ecal_depth_[arrSize_];
  double ecal_rho_[arrSize_];
  int ecal_muonIdx_[arrSize_];
  int ecal_nHits_[arrSize_];

  int nECALHits_;
  
  double ecal_hit_energy_[arrSize_];
  double ecal_hit_depth_[arrSize_];
  double ecal_hit_time_[arrSize_];
  double ecal_hit_pt2_[arrSize_];
  double ecal_hit_x_[arrSize_];
  double ecal_hit_y_[arrSize_];
  double ecal_hit_z_[arrSize_];
  double ecal_hit_eta_[arrSize_];
  double ecal_hit_phi_[arrSize_];
  double ecal_hit_fraction_[arrSize_];
  double ecal_hit_idx_[arrSize_];
  
  int nHCAL_;

  double hcal_et_[arrSize_];
  double hcal_pt_[arrSize_];
  double hcal_eta_[arrSize_];
  double hcal_phi_[arrSize_];
  double hcal_charge_[arrSize_];
  double hcal_px_[arrSize_];
  double hcal_py_[arrSize_];
  double hcal_pz_[arrSize_];
  double hcal_vx_[arrSize_];
  double hcal_vy_[arrSize_];
  double hcal_vz_[arrSize_];
  double hcal_time_[arrSize_];
  double hcal_depth_[arrSize_];
  double hcal_rho_[arrSize_];
  int hcal_muonIdx_[arrSize_];
  int hcal_nHits_[arrSize_];

  int nHCALHits_;

  double hcal_hit_energy_[arrSize_];
  double hcal_hit_depth_[arrSize_];
  double hcal_hit_time_[arrSize_];
  double hcal_hit_pt2_[arrSize_];
  double hcal_hit_x_[arrSize_];
  double hcal_hit_y_[arrSize_];
  double hcal_hit_z_[arrSize_];
  double hcal_hit_eta_[arrSize_];
  double hcal_hit_phi_[arrSize_];
  double hcal_hit_fraction_[arrSize_];
  double hcal_hit_idx_[arrSize_];
  
  // -- trigger info.
  vector< std::string > vec_firedTrigger_;
  vector< std::string > vec_filterName_;
  vector< double > vec_HLTObj_pt_;
  vector< double > vec_HLTObj_eta_;
  vector< double > vec_HLTObj_phi_;

  vector< std::string > vec_myFiredTrigger_;
  vector< std::string > vec_myFilterName_;
  vector< double > vec_myHLTObj_pt_;
  vector< double > vec_myHLTObj_eta_;
  vector< double > vec_myHLTObj_phi_;
  
};

