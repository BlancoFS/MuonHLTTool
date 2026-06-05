// -- ntuple maker for GNN Track - Muon HLT study
// -- author: Sergio Blanco (Institute of Physics of Cantabria)

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

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
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
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

#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Associations/interface/TTClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/TTStubAssociationMap.h"
#include "SimDataFormats/Associations/interface/TTTrackAssociationMap.h"
#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"

#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/ME0DigiCollection.h"

#include "DataFormats/DTRecHit/interface/DTRangeMapAccessor.h"
#include "SimDataFormats/DigiSimLinks/interface/DTDigiSimLinkCollection.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/RPCDigiSimLink/interface/RPCDigiSimLink.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"
#include "SimDataFormats/GEMDigiSimLink/interface/ME0DigiSimLink.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"


#include "TTree.h"
#include "TString.h"

using namespace std;
using namespace reco;
using namespace edm;

class MuonHLTNtupler : public edm::one::EDAnalyzer<>
{
public:
  MuonHLTNtupler(const edm::ParameterSet &iConfig);
  virtual ~MuonHLTNtupler() {};

  virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup);
  virtual void endRun(const edm::Run &iRun, const edm::EventSetup &iSetup);

private:
  void Init();
  void Make_Branch();
  void Fill_Tracks_And_Segments(const edm::Event &iEvent, const edm::EventSetup &iSetup);

  // Tracks
  edm::EDGetTokenT<reco::TrackCollection>           t_hltGeneralTracks_;
  edm::EDGetTokenT<TrajTrackAssociationCollection>  t_trajTrackAssoc_;
  TrackerHitAssociator::Config                      trackerHitAssociatorConfig_;

  edm::EDGetTokenT< edm::SimTrackContainer >        simTrackToken_;

  // Segments
  edm::EDGetTokenT<DTRecSegment4DCollection>        t_dtSegments_;
  edm::EDGetTokenT<CSCSegmentCollection>            t_cscSegments_;
  edm::EDGetTokenT<GEMSegmentCollection>            t_gemSegments_;
  edm::EDGetTokenT<ME0SegmentCollection>            t_me0Segments_;

  edm::EDGetTokenT< DTDigiSimLinkCollection >                dtdigiLinkToken_;
  edm::EDGetTokenT< edm::DetSetVector<StripDigiSimLink> >    cscStripDigiSimLinksToken_;
  edm::EDGetTokenT< edm::DetSetVector<RPCDigiSimLink> >      rpcDigiSimLinksToken_;
  edm::EDGetTokenT<edm::DetSetVector<GEMDigiSimLink>>        gemdigiSimLinkToken_;
  edm::EDGetTokenT<edm::DetSetVector<ME0DigiSimLink>>        me0digiSimLinkToken_;

  // Geometry tokens
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomToken_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  edm::ESGetToken<ME0Geometry, MuonGeometryRecord> me0GeomToken_;


  TTree *ntuple_;
  static const int arrSize_ = 5000;

  static const int MAX_TRACKS       = 5000;
  static const int MAX_DT_SEG       = 500;
  static const int MAX_CSC_SEG      = 1000;
  static const int MAX_GEM_SEG      = 500;
  static const int MAX_ME0_SEG      = 500;

  int   nTracks_;
  float track_pt_[MAX_TRACKS];
  float track_eta_[MAX_TRACKS]; 
  float track_phi_[MAX_TRACKS];
  int   track_charge_[MAX_TRACKS];
  float track_chi2_[MAX_TRACKS];
  int   track_ndof_[MAX_TRACKS];
  int   track_nValidHits_[MAX_TRACKS]; 
  int   track_nLostHits_[MAX_TRACKS];
  float track_vx_[MAX_TRACKS];
  float track_vy_[MAX_TRACKS];  
  float track_vz_[MAX_TRACKS];
  float track_dxy_[MAX_TRACKS];
  float track_dz_[MAX_TRACKS];
  float track_outerX_[MAX_TRACKS];  
  float track_outerY_[MAX_TRACKS];  
  float track_outerZ_[MAX_TRACKS];
  float track_outerPx_[MAX_TRACKS]; 
  float track_outerPy_[MAX_TRACKS]; 
  float track_outerPz_[MAX_TRACKS];
  float track_lastTSOS_x_[MAX_TRACKS];
  float track_lastTSOS_y_[MAX_TRACKS];  
  float track_lastTSOS_z_[MAX_TRACKS];
  float track_lastTSOS_px_[MAX_TRACKS];
  float track_lastTSOS_py_[MAX_TRACKS]; 
  float track_lastTSOS_pz_[MAX_TRACKS];
  float track_lastTSOS_qoverp_[MAX_TRACKS];
  float track_lastTSOS_lambda_[MAX_TRACKS];
  int   track_lastTSOS_isValid_[MAX_TRACKS];
  int   track_simTrackId_[MAX_TRACKS];
  int   track_simTrack_pdgId_[MAX_TRACKS];
  int   track_bestVotes_[MAX_TRACKS];
  float track_simTrack_pt_[MAX_TRACKS]; 
  float track_simTrack_eta_[MAX_TRACKS]; 
  float track_simTrack_phi_[MAX_TRACKS];

  int   nDTSegments_;
  float dtSeg_x_[MAX_DT_SEG], dtSeg_y_[MAX_DT_SEG], dtSeg_z_[MAX_DT_SEG];
  float dtSeg_dx_[MAX_DT_SEG], dtSeg_dy_[MAX_DT_SEG], dtSeg_dz_[MAX_DT_SEG];
  float dtSeg_eta_[MAX_DT_SEG], dtSeg_phi_[MAX_DT_SEG];
  float dtSeg_chi2_[MAX_DT_SEG];
  int   dtSeg_ndof_[MAX_DT_SEG];
  int   dtSeg_station_[MAX_DT_SEG], dtSeg_wheel_[MAX_DT_SEG], dtSeg_sector_[MAX_DT_SEG];
  int   dtSeg_nHits_phi_[MAX_DT_SEG], dtSeg_nHits_z_[MAX_DT_SEG];
  int   dtSeg_hasPhi_[MAX_DT_SEG],  dtSeg_hasZed_[MAX_DT_SEG];
  int   dtSeg_simTrackId_[MAX_DT_SEG], dtSeg_simTrack_pdgId_[MAX_DT_SEG];

  int   nCSCSegments_;
  float cscSeg_x_[MAX_CSC_SEG], cscSeg_y_[MAX_CSC_SEG], cscSeg_z_[MAX_CSC_SEG];
  float cscSeg_dx_[MAX_CSC_SEG], cscSeg_dy_[MAX_CSC_SEG], cscSeg_dz_[MAX_CSC_SEG];
  float cscSeg_eta_[MAX_CSC_SEG], cscSeg_phi_[MAX_CSC_SEG];
  float cscSeg_chi2_[MAX_CSC_SEG];
  int   cscSeg_ndof_[MAX_CSC_SEG], cscSeg_nHits_[MAX_CSC_SEG];
  int   cscSeg_endcap_[MAX_CSC_SEG], cscSeg_station_[MAX_CSC_SEG];
  int   cscSeg_ring_[MAX_CSC_SEG],   cscSeg_chamber_[MAX_CSC_SEG];
  int   cscSeg_simTrackId_[MAX_CSC_SEG], cscSeg_simTrack_pdgId_[MAX_CSC_SEG];

  int   nGEMSegments_;
  float gemSeg_x_[MAX_GEM_SEG], gemSeg_y_[MAX_GEM_SEG], gemSeg_z_[MAX_GEM_SEG];
  float gemSeg_dx_[MAX_GEM_SEG], gemSeg_dy_[MAX_GEM_SEG], gemSeg_dz_[MAX_GEM_SEG];
  float gemSeg_eta_[MAX_GEM_SEG], gemSeg_phi_[MAX_GEM_SEG];
  float gemSeg_chi2_[MAX_GEM_SEG];
  int   gemSeg_ndof_[MAX_GEM_SEG], gemSeg_nHits_[MAX_GEM_SEG];
  int   gemSeg_region_[MAX_GEM_SEG], gemSeg_station_[MAX_GEM_SEG];
  int   gemSeg_ring_[MAX_GEM_SEG],   gemSeg_chamber_[MAX_GEM_SEG];
  int   gemSeg_simTrackId_[MAX_GEM_SEG], gemSeg_simTrack_pdgId_[MAX_GEM_SEG];

  int   nME0Segments_;
  float me0Seg_x_[MAX_ME0_SEG], me0Seg_y_[MAX_ME0_SEG], me0Seg_z_[MAX_ME0_SEG];
  float me0Seg_dx_[MAX_ME0_SEG], me0Seg_dy_[MAX_ME0_SEG], me0Seg_dz_[MAX_ME0_SEG];
  float me0Seg_eta_[MAX_ME0_SEG], me0Seg_phi_[MAX_ME0_SEG];
  float me0Seg_chi2_[MAX_ME0_SEG];
  int   me0Seg_ndof_[MAX_ME0_SEG], me0Seg_nHits_[MAX_ME0_SEG];
  int   me0Seg_region_[MAX_ME0_SEG], me0Seg_station_[MAX_ME0_SEG], me0Seg_chamber_[MAX_ME0_SEG];
  int   me0Seg_simTrackId_[MAX_ME0_SEG], me0Seg_simTrack_pdgId_[MAX_ME0_SEG];

};