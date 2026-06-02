// -- ntuple maker for GNN Track - Muon HLT study
// -- author: Sergio Blanco (Institute of Physics of Cantabria)

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



#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"

using namespace std;
using namespace reco;
using namespace edm;


MuonHLTNtupler::MuonHLTNtupler(const edm::ParameterSet& iConfig):
t_hltGeneralTracks_     (consumes<reco::TrackCollection>          (iConfig.getParameter<edm::InputTag>("hltGeneralTracks"))),
t_trajTrackAssoc_       (consumes<TrajTrackAssociationCollection> (iConfig.getParameter<edm::InputTag>("trajTrackAssoc"))),
trackerHitAssociatorConfig_(iConfig, consumesCollector()),

simTrackToken_ (consumes<edm::SimTrackContainer> (iConfig.getUntrackedParameter<edm::InputTag>("simTrackLabel")) ),

t_dtSegments_           (consumes<DTRecSegment4DCollection>       (iConfig.getParameter<edm::InputTag>("dtSegments"))),
t_cscSegments_          (consumes<CSCSegmentCollection>           (iConfig.getParameter<edm::InputTag>("cscSegments"))),
t_gemSegments_          (consumes<GEMSegmentCollection>           (iConfig.getParameter<edm::InputTag>("gemSegments"))),
t_me0Segments_          (consumes<ME0SegmentCollection>           (iConfig.getParameter<edm::InputTag>("me0Segments"))),

dtdigiLinkToken_ (consumes<DTDigiSimLinkCollection> (iConfig.getUntrackedParameter<edm::InputTag>("dtDigiSimLinkLabel")) ),
cscStripDigiSimLinksToken_ (consumes<edm::DetSetVector<StripDigiSimLink>> (iConfig.getUntrackedParameter<edm::InputTag>("cscStripDigiSimLinkLabel")) ),
rpcDigiSimLinksToken_ (consumes<edm::DetSetVector<RPCDigiSimLink>> (iConfig.getUntrackedParameter<edm::InputTag>("rpcDigiSimLinkLabel")) ),
gemdigiSimLinkToken_    (consumes<edm::DetSetVector<GEMDigiSimLink>>  (iConfig.getParameter<edm::InputTag>("gemDigiSimLinks"))),
me0digiSimLinkToken_    (consumes<edm::DetSetVector<ME0DigiSimLink>>  (iConfig.getParameter<edm::InputTag>("me0DigiSimLinks")))
{
  dtGeomToken_ = esConsumes();
  cscGeomToken_ = esConsumes();
  gemGeomToken_ = esConsumes();
  me0GeomToken_ = esConsumes();
}

void MuonHLTNtupler::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  Init(); 

  // -- fill each object
  Fill_Tracks_And_Segments(iEvent, iSetup);

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
  
  nTracks_ = 0;
  for ( int i=0; i<MAX_TRACKS; i++)
  {
    track_pt_[i] = -999;
    track_eta_[i] = -999;
    track_phi_[i] = -999;
    track_charge_[i] = -999;
    track_chi2_[i] = -999;
    track_ndof_[i] = -999;
    track_nValidHits_[i] = -999;
    track_nLostHits_[i] = -999;
    track_vx_[i] = -999;
    track_vy_[i] = -999;
    track_vz_[i] = -999;
    track_dxy_[i] = -999;
    track_dz_[i] = -999;
    track_outerX_[i] = -999;
    track_outerY_[i] = -999;
    track_outerZ_[i] = -999;
    track_outerPx_[i] = -999;
    track_outerPy_[i] = -999;
    track_outerPz_[i] = -999;
    track_lastTSOS_x_[i] = -999;
    track_lastTSOS_y_[i] = -999;
    track_lastTSOS_z_[i] = -999;
    track_lastTSOS_px_[i] = -999;
    track_lastTSOS_py_[i] = -999;
    track_lastTSOS_pz_[i] = -999;
    track_lastTSOS_qoverp_[i] = -999;
    track_lastTSOS_lambda_[i] = -999;
    track_lastTSOS_isValid_[i] = 0;
    track_simTrackId_[i] = -1;
    track_simTrack_pdgId_[i] = -999;
    track_simTrack_pt_[i] = -999;
    track_simTrack_eta_[i] = -999;
    track_simTrack_phi_[i] = -999;
  }

  nDTSegments_ = 0;
  for ( int i=0; i<MAX_DT_SEG; i++)
  {
    dtSeg_x_[i] = -999;
    dtSeg_y_[i] = -999;
    dtSeg_z_[i] = -999;
    dtSeg_dx_[i] = -999;
    dtSeg_dy_[i] = -999;
    dtSeg_dz_[i] = -999;
    dtSeg_eta_[i] = -999;
    dtSeg_phi_[i] = -999;
    dtSeg_chi2_[i] = -999;
    dtSeg_ndof_[i] = -999;
    dtSeg_station_[i] = -1;
    dtSeg_wheel_[i] = -1;
    dtSeg_sector_[i] = -1;
    dtSeg_nHits_phi_[i] = -1;
    dtSeg_nHits_z_[i] = -1;
    dtSeg_hasPhi_[i] = 0;
    dtSeg_hasZed_[i] = 0;
    dtSeg_simTrackId_[i] = -1;
    dtSeg_simTrack_pdgId_[i] = -999;
  }

  nCSCSegments_ = 0;
  for ( int i=0; i<MAX_CSC_SEG; i++)
  {
    cscSeg_x_[i] = -999;
    cscSeg_y_[i] = -999;
    cscSeg_z_[i] = -999;
    cscSeg_dx_[i] = -999;
    cscSeg_dy_[i] = -999;
    cscSeg_dz_[i] = -999;
    cscSeg_eta_[i] = -999;
    cscSeg_phi_[i] = -999;
    cscSeg_chi2_[i] = -999;
    cscSeg_ndof_[i] = -999;
    cscSeg_nHits_[i] = -1;
    cscSeg_endcap_[i] = -1;
    cscSeg_station_[i] = -1;
    cscSeg_ring_[i] = -1;
    cscSeg_chamber_[i] = -1;
    cscSeg_simTrackId_[i] = -1;
    cscSeg_simTrack_pdgId_[i] = -999;
  }

  nGEMSegments_ = 0;
  for ( int i=0; i<MAX_GEM_SEG; i++)
  {
    gemSeg_x_[i] = -999;
    gemSeg_y_[i] = -999;
    gemSeg_z_[i] = -999;
    gemSeg_dx_[i] = -999;
    gemSeg_dy_[i] = -999;
    gemSeg_dz_[i] = -999;
    gemSeg_eta_[i] = -999;
    gemSeg_phi_[i] = -999;
    gemSeg_chi2_[i] = -999;
    gemSeg_ndof_[i] = -999;
    gemSeg_nHits_[i] = -1;
    gemSeg_region_[i] = -1;
    gemSeg_station_[i] = -1;
    gemSeg_ring_[i] = -1;
    gemSeg_chamber_[i] = -1;
    gemSeg_simTrackId_[i] = -1;
    gemSeg_simTrack_pdgId_[i] = -999;
  }

  nME0Segments_ = 0;
  for ( int i=0; i<MAX_ME0_SEG; i++)
  {
    me0Seg_x_[i] = -999;
    me0Seg_y_[i] = -999;
    me0Seg_z_[i] = -999;
    me0Seg_dx_[i] = -999;
    me0Seg_dy_[i] = -999;
    me0Seg_dz_[i] = -999;
    me0Seg_eta_[i] = -999;
    me0Seg_phi_[i] = -999;
    me0Seg_chi2_[i] = -999;
    me0Seg_ndof_[i] = -999;
    me0Seg_nHits_[i] = -1;
    me0Seg_region_[i] = -1;
    me0Seg_station_[i] = -1;
    me0Seg_chamber_[i] = -1;
    me0Seg_simTrackId_[i] = -1;
    me0Seg_simTrack_pdgId_[i] = -999;
  }

}

void MuonHLTNtupler::Make_Branch()
{

  ntuple_->Branch("nTracks", &nTracks_, "nTracks/I");
  ntuple_->Branch("track_pt", &track_pt_, "track_pt[nTracks]/F");
  ntuple_->Branch("track_eta", &track_eta_, "track_eta[nTracks]/F");
  ntuple_->Branch("track_phi", &track_phi_, "track_phi[nTracks]/F");
  ntuple_->Branch("track_charge", &track_charge_, "track_charge[nTracks]/I");
  ntuple_->Branch("track_chi2", &track_chi2_, "track_chi2[nTracks]/F");
  ntuple_->Branch("track_ndof", &track_ndof_, "track_ndof[nTracks]/I");
  ntuple_->Branch("track_nValidHits", &track_nValidHits_, "track_nValidHits[nTracks]/I");
  ntuple_->Branch("track_nLostHits", &track_nLostHits_, "track_nLostHits[nTracks]/I");
  ntuple_->Branch("track_vx", &track_vx_, "track_vx[nTracks]/F");
  ntuple_->Branch("track_vy", &track_vy_, "track_vy[nTracks]/F");
  ntuple_->Branch("track_vz", &track_vz_, "track_vz[nTracks]/F");
  ntuple_->Branch("track_dxy", &track_dxy_, "track_dxy[nTracks]/F");
  ntuple_->Branch("track_dz", &track_dz_, "track_dz[nTracks]/F");
  ntuple_->Branch("track_outerX", &track_outerX_, "track_outerX[nTracks]/F");
  ntuple_->Branch("track_outerY", &track_outerY_, "track_outerY[nTracks]/F");
  ntuple_->Branch("track_outerZ", &track_outerZ_, "track_outerZ[nTracks]/F");
  ntuple_->Branch("track_outerPx", &track_outerPx_, "track_outerPx[nTracks]/F");
  ntuple_->Branch("track_outerPy", &track_outerPy_, "track_outerPy[nTracks]/F");
  ntuple_->Branch("track_outerPz", &track_outerPz_, "track_outerPz[nTracks]/F");
  ntuple_->Branch("track_lastTSOS_x", &track_lastTSOS_x_, "track_lastTSOS_x[nTracks]/F");
  ntuple_->Branch("track_lastTSOS_y", &track_lastTSOS_y_, "track_lastTSOS_y[nTracks]/F");
  ntuple_->Branch("track_lastTSOS_z", &track_lastTSOS_z_, "track_lastTSOS_z[nTracks]/F");
  ntuple_->Branch("track_lastTSOS_px", &track_lastTSOS_px_, "track_lastTSOS_px[nTracks]/F");
  ntuple_->Branch("track_lastTSOS_py", &track_lastTSOS_py_, "track_lastTSOS_py[nTracks]/F");
  ntuple_->Branch("track_lastTSOS_pz", &track_lastTSOS_pz_, "track_lastTSOS_pz[nTracks]/F");
  ntuple_->Branch("track_lastTSOS_qoverp", &track_lastTSOS_qoverp_, "track_lastTSOS_qoverp[nTracks]/F");
  ntuple_->Branch("track_lastTSOS_lambda", &track_lastTSOS_lambda_, "track_lastTSOS_lambda[nTracks]/F");
  ntuple_->Branch("track_lastTSOS_isValid", &track_lastTSOS_isValid_, "track_lastTSOS_isValid[nTracks]/I");
  ntuple_->Branch("track_simTrackId", &track_simTrackId_, "track_simTrackId[nTracks]/I");
  ntuple_->Branch("track_simTrack_pdgId", &track_simTrack_pdgId_, "track_simTrack_pdgId[nTracks]/I");
  ntuple_->Branch("track_simTrack_pt", &track_simTrack_pt_, "track_simTrack_pt[nTracks]/F");
  ntuple_->Branch("track_simTrack_eta", &track_simTrack_eta_, "track_simTrack_eta[nTracks]/F");
  ntuple_->Branch("track_simTrack_phi", &track_simTrack_phi_, "track_simTrack_phi[nTracks]/F");

  ntuple_->Branch("nDTSegments", &nDTSegments_, "nDTSegments/I");
  ntuple_->Branch("dtSeg_x", &dtSeg_x_, "dtSeg_x[nDTSegments]/F");
  ntuple_->Branch("dtSeg_y", &dtSeg_y_, "dtSeg_y[nDTSegments]/F");
  ntuple_->Branch("dtSeg_z", &dtSeg_z_, "dtSeg_z[nDTSegments]/F");
  ntuple_->Branch("dtSeg_dx", &dtSeg_dx_, "dtSeg_dx[nDTSegments]/F");
  ntuple_->Branch("dtSeg_dy", &dtSeg_dy_, "dtSeg_dy[nDTSegments]/F");
  ntuple_->Branch("dtSeg_dz", &dtSeg_dz_, "dtSeg_dz[nDTSegments]/F");
  ntuple_->Branch("dtSeg_eta", &dtSeg_eta_, "dtSeg_eta[nDTSegments]/F");
  ntuple_->Branch("dtSeg_phi", &dtSeg_phi_, "dtSeg_phi[nDTSegments]/F");
  ntuple_->Branch("dtSeg_chi2", &dtSeg_chi2_, "dtSeg_chi2[nDTSegments]/F");
  ntuple_->Branch("dtSeg_ndof", &dtSeg_ndof_, "dtSeg_ndof[nDTSegments]/I");
  ntuple_->Branch("dtSeg_station", &dtSeg_station_, "dtSeg_station[nDTSegments]/I");
  ntuple_->Branch("dtSeg_wheel", &dtSeg_wheel_, "dtSeg_wheel[nDTSegments]/I");
  ntuple_->Branch("dtSeg_sector", &dtSeg_sector_, "dtSeg_sector[nDTSegments]/I");
  ntuple_->Branch("dtSeg_nHits_phi", &dtSeg_nHits_phi_, "dtSeg_nHits_phi[nDTSegments]/I");
  ntuple_->Branch("dtSeg_nHits_z", &dtSeg_nHits_z_, "dtSeg_nHits_z[nDTSegments]/I");
  ntuple_->Branch("dtSeg_hasPhi", &dtSeg_hasPhi_, "dtSeg_hasPhi[nDTSegments]/I");
  ntuple_->Branch("dtSeg_hasZed", &dtSeg_hasZed_, "dtSeg_hasZed[nDTSegments]/I");
  ntuple_->Branch("dtSeg_simTrackId", &dtSeg_simTrackId_, "dtSeg_simTrackId[nDTSegments]/I");
  ntuple_->Branch("dtSeg_simTrack_pdgId", &dtSeg_simTrack_pdgId_, "dtSeg_simTrack_pdgId[nDTSegments]/I");

  ntuple_->Branch("nCSCSegments", &nCSCSegments_, "nCSCSegments/I");
  ntuple_->Branch("cscSeg_x", &cscSeg_x_, "cscSeg_x[nCSCSegments]/F");
  ntuple_->Branch("cscSeg_y", &cscSeg_y_, "cscSeg_y[nCSCSegments]/F");
  ntuple_->Branch("cscSeg_z", &cscSeg_z_, "cscSeg_z[nCSCSegments]/F");
  ntuple_->Branch("cscSeg_dx", &cscSeg_dx_, "cscSeg_dx[nCSCSegments]/F");
  ntuple_->Branch("cscSeg_dy", &cscSeg_dy_, "cscSeg_dy[nCSCSegments]/F");
  ntuple_->Branch("cscSeg_dz", &cscSeg_dz_, "cscSeg_dz[nCSCSegments]/F");
  ntuple_->Branch("cscSeg_eta", &cscSeg_eta_, "cscSeg_eta[nCSCSegments]/F");
  ntuple_->Branch("cscSeg_phi", &cscSeg_phi_, "cscSeg_phi[nCSCSegments]/F");
  ntuple_->Branch("cscSeg_chi2", &cscSeg_chi2_, "cscSeg_chi2[nCSCSegments]/F");
  ntuple_->Branch("cscSeg_ndof", &cscSeg_ndof_, "cscSeg_ndof[nCSCSegments]/I");
  ntuple_->Branch("cscSeg_nHits", &cscSeg_nHits_, "cscSeg_nHits[nCSCSegments]/I");
  ntuple_->Branch("cscSeg_endcap", &cscSeg_endcap_, "cscSeg_endcap[nCSCSegments]/I");
  ntuple_->Branch("cscSeg_station", &cscSeg_station_, "cscSeg_station[nCSCSegments]/I");
  ntuple_->Branch("cscSeg_ring", &cscSeg_ring_, "cscSeg_ring[nCSCSegments]/I");
  ntuple_->Branch("cscSeg_chamber", &cscSeg_chamber_, "cscSeg_chamber[nCSCSegments]/I");
  ntuple_->Branch("cscSeg_simTrackId", &cscSeg_simTrackId_, "cscSeg_simTrackId[nCSCSegments]/I");
  ntuple_->Branch("cscSeg_simTrack_pdgId", &cscSeg_simTrack_pdgId_, "cscSeg_simTrack_pdgId[nCSCSegments]/I");

  ntuple_->Branch("nGEMSegments", &nGEMSegments_, "nGEMSegments/I");
  ntuple_->Branch("gemSeg_x", &gemSeg_x_, "gemSeg_x[nGEMSegments]/F");
  ntuple_->Branch("gemSeg_y", &gemSeg_y_, "gemSeg_y[nGEMSegments]/F");
  ntuple_->Branch("gemSeg_z", &gemSeg_z_, "gemSeg_z[nGEMSegments]/F");
  ntuple_->Branch("gemSeg_dx", &gemSeg_dx_, "gemSeg_dx[nGEMSegments]/F");
  ntuple_->Branch("gemSeg_dy", &gemSeg_dy_, "gemSeg_dy[nGEMSegments]/F");
  ntuple_->Branch("gemSeg_dz", &gemSeg_dz_, "gemSeg_dz[nGEMSegments]/F");
  ntuple_->Branch("gemSeg_eta", &gemSeg_eta_, "gemSeg_eta[nGEMSegments]/F");
  ntuple_->Branch("gemSeg_phi", &gemSeg_phi_, "gemSeg_phi[nGEMSegments]/F");
  ntuple_->Branch("gemSeg_chi2", &gemSeg_chi2_, "gemSeg_chi2[nGEMSegments]/F");
  ntuple_->Branch("gemSeg_ndof", &gemSeg_ndof_, "gemSeg_ndof[nGEMSegments]/I");
  ntuple_->Branch("gemSeg_nHits", &gemSeg_nHits_, "gemSeg_nHits[nGEMSegments]/I");
  ntuple_->Branch("gemSeg_region", &gemSeg_region_, "gemSeg_region[nGEMSegments]/I");
  ntuple_->Branch("gemSeg_station", &gemSeg_station_, "gemSeg_station[nGEMSegments]/I");
  ntuple_->Branch("gemSeg_ring", &gemSeg_ring_, "gemSeg_ring[nGEMSegments]/I");
  ntuple_->Branch("gemSeg_chamber", &gemSeg_chamber_, "gemSeg_chamber[nGEMSegments]/I");
  ntuple_->Branch("gemSeg_simTrackId", &gemSeg_simTrackId_, "gemSeg_simTrackId[nGEMSegments]/I");
  ntuple_->Branch("gemSeg_simTrack_pdgId", &gemSeg_simTrack_pdgId_, "gemSeg_simTrack_pdgId[nGEMSegments]/I");

  ntuple_->Branch("nME0Segments", &nME0Segments_, "nME0Segments/I");
  ntuple_->Branch("me0Seg_x", &me0Seg_x_, "me0Seg_x[nME0Segments]/F");
  ntuple_->Branch("me0Seg_y", &me0Seg_y_, "me0Seg_y[nME0Segments]/F");
  ntuple_->Branch("me0Seg_z", &me0Seg_z_, "me0Seg_z[nME0Segments]/F");
  ntuple_->Branch("me0Seg_dx", &me0Seg_dx_, "me0Seg_dx[nME0Segments]/F");
  ntuple_->Branch("me0Seg_dy", &me0Seg_dy_, "me0Seg_dy[nME0Segments]/F");
  ntuple_->Branch("me0Seg_dz", &me0Seg_dz_, "me0Seg_dz[nME0Segments]/F");
  ntuple_->Branch("me0Seg_eta", &me0Seg_eta_, "me0Seg_eta[nME0Segments]/F");
  ntuple_->Branch("me0Seg_phi", &me0Seg_phi_, "me0Seg_phi[nME0Segments]/F");
  ntuple_->Branch("me0Seg_chi2", &me0Seg_chi2_, "me0Seg_chi2[nME0Segments]/F");
  ntuple_->Branch("me0Seg_ndof", &me0Seg_ndof_, "me0Seg_ndof[nME0Segments]/I");
  ntuple_->Branch("me0Seg_nHits", &me0Seg_nHits_, "me0Seg_nHits[nME0Segments]/I");
  ntuple_->Branch("me0Seg_region", &me0Seg_region_, "me0Seg_region[nME0Segments]/I");
  ntuple_->Branch("me0Seg_station", &me0Seg_station_, "me0Seg_station[nME0Segments]/I");
  ntuple_->Branch("me0Seg_chamber", &me0Seg_chamber_, "me0Seg_chamber[nME0Segments]/I");
  ntuple_->Branch("me0Seg_simTrackId", &me0Seg_simTrackId_, "me0Seg_simTrackId[nME0Segments]/I");
  ntuple_->Branch("me0Seg_simTrack_pdgId", &me0Seg_simTrack_pdgId_, "me0Seg_simTrack_pdgId[nME0Segments]/I");

}

// =============================================================================
// Fill_Tracks_And_Segments — MuonHLTNtupler
// =============================================================================
//
// This function does three things:
//   1. Reads hltGeneralTracks, matches each to a SimMuon via TrackerHitAssociator
//      majority voting, and stores kinematics, vertex, outer/inner states, and
//      the last TrajectoryStateOnSurface (TSOS) from the associated Trajectory.
//   2. Reads DT / CSC / GEM / ME0 segments, stores their geometry and basic
//      quality variables, and assigns a SimMuon track ID to each via majority
//      voting over the component RecHit → digi → DigiSimLink chain.
//
// =============================================================================


void MuonHLTNtupler::Fill_Tracks_And_Segments(const edm::Event& iEvent,
                                               const edm::EventSetup& iSetup) {

  bool debug = true;

  if (debug) {
    std::cout << "\n\n";
    std::cout << "-------------------------------------------------------------------\n";
    std::cout << "-------------- Fill_Tracks_And_Segments: New Event ----------------\n";
    std::cout << "-------------------------------------------------------------------\n\n";
  }

  // ===========================================================================
  // 0.  Retrieve shared handles
  // ===========================================================================

  // SimTracks
  edm::Handle<edm::SimTrackContainer> simTracks;
  iEvent.getByToken(simTrackToken_, simTracks);

  // DT/CSC digi–sim links
  edm::Handle<DTDigiSimLinkCollection> dtdigiLinks;
  iEvent.getByToken(dtdigiLinkToken_, dtdigiLinks);

  edm::Handle<edm::DetSetVector<StripDigiSimLink>> cscStripDigiSimLink;
  iEvent.getByToken(cscStripDigiSimLinksToken_, cscStripDigiSimLink);

  // GEM / ME0 digi–sim links
  edm::Handle<edm::DetSetVector<GEMDigiSimLink>> gemDigiSimLinks;
  iEvent.getByToken(gemdigiSimLinkToken_, gemDigiSimLinks);

  edm::Handle<edm::DetSetVector<ME0DigiSimLink>> me0DigiSimLinks;
  iEvent.getByToken(me0digiSimLinkToken_, me0DigiSimLinks);

  // Geometry
  ESHandle<DTGeometry> dtGeom = iSetup.getHandle(dtGeomToken_);
  ESHandle<CSCGeometry> cscGeom = iSetup.getHandle(cscGeomToken_);
  ESHandle<GEMGeometry> gemGeom = iSetup.getHandle(gemGeomToken_);
  ESHandle<ME0Geometry> me0Geom = iSetup.getHandle(me0GeomToken_);

  // const DTGeometry&  dtGeom  = iSetup.getData(dtGeomToken_);
  // const CSCGeometry& cscGeom = iSetup.getData(cscGeomToken_);
  // const GEMGeometry& gemGeom = iSetup.getData(gemGeomToken_);
  // const ME0Geometry& me0Geom = iSetup.getData(me0GeomToken_);

  // --------------------------------------------------------------------------
  // Build a fast lookup map: simTrackId → PDG id  (primary muons only)
  // --------------------------------------------------------------------------
  std::map<unsigned int, int> simTrackToPdgId;
  std::map<unsigned int, const SimTrack*> simTrackById;
  for (const auto& st : *simTracks) {
    if (std::abs(st.type()) != 13) continue;
    if (!st.isPrimary())           continue;
    simTrackToPdgId[st.trackId()] = st.type();
    simTrackById  [st.trackId()] = &st;
  }

  // ===========================================================================
  // 1.  HLT General Tracks
  // ===========================================================================

  edm::Handle<reco::TrackCollection> hltTracks;
  iEvent.getByToken(t_hltGeneralTracks_, hltTracks);

  edm::Handle<TrajTrackAssociationCollection> trajTrackAssoc;
  iEvent.getByToken(t_trajTrackAssoc_, trajTrackAssoc);

  // Build a quick map  TrackRef → TrajectoryRef  from the association map
  std::map<reco::TrackRef, edm::Ref<std::vector<Trajectory>>> trackToTraj;
  if (trajTrackAssoc.isValid()) {
    for (const auto& assoc : *trajTrackAssoc)
      trackToTraj[assoc.val] = assoc.key;
  }

  // TrackerHitAssociator: maps each tracker RecHit to its contributing SimTrack(s)
  TrackerHitAssociator hitAssociator(iEvent, trackerHitAssociatorConfig_);

  int _nTracks = 0;

  if (hltTracks.isValid()) {
    for (size_t iTrack = 0; iTrack < hltTracks->size(); ++iTrack) {

      const reco::Track& track = (*hltTracks)[iTrack];
      reco::TrackRef    tRef(hltTracks, iTrack);

      // -- Basic kinematics and fit quality --
      track_pt_         [_nTracks] = track.pt();
      track_eta_        [_nTracks] = track.eta();
      track_phi_        [_nTracks] = track.phi();
      track_charge_     [_nTracks] = track.charge();
      track_chi2_       [_nTracks] = track.normalizedChi2();
      track_ndof_       [_nTracks] = track.ndof();
      track_nValidHits_ [_nTracks] = track.numberOfValidHits();
      track_nLostHits_  [_nTracks] = track.numberOfLostHits();

      // -- Vertex / impact parameters --
      track_vx_  [_nTracks] = track.vx();
      track_vy_  [_nTracks] = track.vy();
      track_vz_  [_nTracks] = track.vz();
      track_dxy_ [_nTracks] = track.dxy();
      track_dz_  [_nTracks] = track.dz();

      // -- Outer measurement (last hit in the tracker, closest to muon system) --
      // track.outerX/Y/Z() and outerPx/Py/Pz() are filled by the track builder
      // directly from the outermost valid RecHit's state.
      track_outerX_  [_nTracks] = track.outerX();
      track_outerY_  [_nTracks] = track.outerY();
      track_outerZ_  [_nTracks] = track.outerZ();
      track_outerPx_ [_nTracks] = track.outerPx();
      track_outerPy_ [_nTracks] = track.outerPy();
      track_outerPz_ [_nTracks] = track.outerPz();

      // -- Last TSOS from the Trajectory object --
      // The Trajectory carries the full Kalman-filter state at every measurement.
      // lastMeasurement() gives the outermost one; we prefer the updated state
      // (post-hit update) and fall back to the forward prediction otherwise.
      track_lastTSOS_x_      [_nTracks] = -999.f;
      track_lastTSOS_y_      [_nTracks] = -999.f;
      track_lastTSOS_z_      [_nTracks] = -999.f;
      track_lastTSOS_px_     [_nTracks] = -999.f;
      track_lastTSOS_py_     [_nTracks] = -999.f;
      track_lastTSOS_pz_     [_nTracks] = -999.f;
      track_lastTSOS_qoverp_ [_nTracks] = -999.f;
      track_lastTSOS_lambda_ [_nTracks] = -999.f;
      track_lastTSOS_isValid_[_nTracks] = 0;

      auto trajIt = trackToTraj.find(tRef);
      if (trajIt != trackToTraj.end() && trajIt->second.isNonnull()) {
        const Trajectory& traj = *(trajIt->second);
        if (!traj.empty()) {
          const TrajectoryMeasurement& lastMeas = traj.lastMeasurement();

          // Prefer the updated (smoothed) TSOS; fall back to forward prediction
          const TrajectoryStateOnSurface& lastTSOS =
            lastMeas.updatedState().isValid() ? lastMeas.updatedState()
                                              : lastMeas.forwardPredictedState();

          if (lastTSOS.isValid()) {
            const GlobalPoint&  pos = lastTSOS.globalPosition();
            const GlobalVector& mom = lastTSOS.globalMomentum();

            track_lastTSOS_x_  [_nTracks] = pos.x();
            track_lastTSOS_y_  [_nTracks] = pos.y();
            track_lastTSOS_z_  [_nTracks] = pos.z();
            track_lastTSOS_px_ [_nTracks] = mom.x();
            track_lastTSOS_py_ [_nTracks] = mom.y();
            track_lastTSOS_pz_ [_nTracks] = mom.z();

            // q/p from local parameters (same sign convention as Kalman)
            track_lastTSOS_qoverp_[_nTracks] =
              lastTSOS.localParameters().qbp() * mom.mag(); // = q / |p|

            // lambda = pi/2 - polar angle (dip angle, signed)
            track_lastTSOS_lambda_[_nTracks] =
              static_cast<float>(M_PI / 2.0) - mom.theta();

            track_lastTSOS_isValid_[_nTracks] = 1;

            if (debug) {
              std::cout << "  Track " << _nTracks
                        << ": last TSOS at (" << pos.x() << ", "
                        << pos.y() << ", " << pos.z() << ")"
                        << ", |p|=" << mom.mag()
                        << ", lambda=" << track_lastTSOS_lambda_[_nTracks]
                        << std::endl;
            }
          }
        }
      }

      // -- SimTrack matching: majority-vote over tracker RecHits --
      // For each valid hit on the track, TrackerHitAssociator returns the list
      // of (SimTrackId, eventId) pairs that deposited charge there.  The SimTrack
      // that appears most often across all hits "owns" the reconstructed track.
      track_simTrackId_      [_nTracks] = -1;
      track_simTrack_pdgId_  [_nTracks] = -999;
      track_simTrack_pt_     [_nTracks] = -999.f;
      track_simTrack_eta_    [_nTracks] = -999.f;
      track_simTrack_phi_    [_nTracks] = -999.f;

      std::map<std::pair<unsigned int, unsigned int>, int> simHitVotes;
      for (auto hitIt = track.recHitsBegin(); hitIt != track.recHitsEnd(); ++hitIt) {
        
        const TrackingRecHit* hit = *hitIt;
        // const TrackingRecHit* hit = hitIt->get();

        if (!hit->isValid()) continue;
        for (const SimHitIdpr& id : hitAssociator.associateHitId(*hit)){
          std::pair<unsigned int, unsigned int> tmp_id(id.first, id.second.event());
          simHitVotes[tmp_id]++;
        }
      }

      unsigned int bestId    = 0;
      int          bestVotes = 0;
      for (const auto& kv : simHitVotes) {
        if (kv.second > bestVotes) {
          bestVotes = kv.second;
          bestId    = kv.first.first; // SimTrack ID
        }
      }

      if (bestVotes > 0) {
        auto stIt = simTrackById.find(bestId);
        if (stIt != simTrackById.end()) {
          const SimTrack* st = stIt->second;
          track_simTrackId_    [_nTracks] = static_cast<int>(bestId);
          track_simTrack_pdgId_[_nTracks] = st->type();
          track_simTrack_pt_   [_nTracks] = st->momentum().pt();
          track_simTrack_eta_  [_nTracks] = st->momentum().eta();
          track_simTrack_phi_  [_nTracks] = st->momentum().phi();
        }
      }

      if (debug) {
        std::cout << "Track " << _nTracks
                  << ": pT=" << track.pt()
                  << " eta=" << track.eta()
                  << " phi=" << track.phi()
                  << " nHits=" << track.numberOfValidHits()
                  << " → SimTrack " << track_simTrackId_[_nTracks]
                  << " (pdgId=" << track_simTrack_pdgId_[_nTracks] << ")"
                  << std::endl;
      }

      ++_nTracks;
      if (_nTracks >= MAX_TRACKS) {
        edm::LogWarning("MuonHLTNtupler") << "MAX_TRACKS reached, truncating.";
        break;
      }
    }
  }
  nTracks_ = _nTracks;

  if (debug)
    std::cout << "Total tracks stored: " << nTracks_ << "\n";

  // ===========================================================================
  // Helper lambda: majority-vote SimMuon ID from an arbitrary vote map
  // ===========================================================================
  auto bestSimMuon = [&](const std::map<unsigned int, int>& votes)
                         -> std::pair<int,int> // {simTrackId, pdgId}
  {
    unsigned int winId    = 0;
    int          winVotes = 0;
    for (const auto& kv : votes) {
      if (kv.second > winVotes) {
        winVotes = kv.second;
        winId    = kv.first;
      }
    }
    if (winVotes == 0) return {-1, -999};
    auto it = simTrackToPdgId.find(winId);
    int pdg = (it != simTrackToPdgId.end()) ? it->second : -999;
    return {static_cast<int>(winId), pdg};
  };

  // ===========================================================================
  // 2.  DT Segments
  // ===========================================================================

  edm::Handle<DTRecSegment4DCollection> dtSegments;
  iEvent.getByToken(t_dtSegments_, dtSegments);

  int _nDTSeg = 0;

  if (dtSegments.isValid()) {
    for (auto segIt = dtSegments->begin(); segIt != dtSegments->end(); ++segIt) {

      const DTRecSegment4D& seg = *segIt;
      DTChamberId cid           = seg.chamberId();
      const GeomDet* gdet       = dtGeom->idToDet(cid);

      GlobalPoint  gp = gdet->toGlobal(seg.localPosition());
      GlobalVector gv = gdet->toGlobal(seg.localDirection());

      // Geometry
      dtSeg_x_  [_nDTSeg] = gp.x();
      dtSeg_y_  [_nDTSeg] = gp.y();
      dtSeg_z_  [_nDTSeg] = gp.z();
      dtSeg_dx_ [_nDTSeg] = gv.x();
      dtSeg_dy_ [_nDTSeg] = gv.y();
      dtSeg_dz_ [_nDTSeg] = gv.z();
      dtSeg_eta_[_nDTSeg] = gp.eta();
      dtSeg_phi_[_nDTSeg] = gp.phi();

      // DetId
      dtSeg_station_[_nDTSeg] = cid.station();
      dtSeg_wheel_  [_nDTSeg] = cid.wheel();
      dtSeg_sector_ [_nDTSeg] = cid.sector();

      // Quality
      dtSeg_chi2_    [_nDTSeg] = static_cast<float>(seg.chi2());
      dtSeg_ndof_    [_nDTSeg] = seg.degreesOfFreedom();
      dtSeg_hasPhi_  [_nDTSeg] = seg.hasPhi() ? 1 : 0;
      dtSeg_hasZed_  [_nDTSeg] = seg.hasZed()  ? 1 : 0;

      // Collect component 1D hits (phi and Z super-layers)
      std::vector<const DTRecHit1D*> hits1D;
      dtSeg_nHits_phi_[_nDTSeg] = 0;
      dtSeg_nHits_z_  [_nDTSeg] = 0;

      if (seg.hasPhi()) {
        const DTChamberRecSegment2D* phiSeg = seg.phiSegment();
        if (phiSeg) {
          const auto& phiHits = phiSeg->specificRecHits();
          dtSeg_nHits_phi_[_nDTSeg] = static_cast<int>(phiHits.size());
          for (const DTRecHit1D& h : phiHits)
            hits1D.push_back(&h);
        }
      }
      if (seg.hasZed()) {
        const DTSLRecSegment2D* zSeg = seg.zSegment();
        if (zSeg) {
          const auto& zHits = zSeg->specificRecHits();
          dtSeg_nHits_z_[_nDTSeg] = static_cast<int>(zHits.size());
          for (const DTRecHit1D& h : zHits)
            hits1D.push_back(&h);
        }
      }

      // SimMuon matching:
      //   For each component 1D hit → find its layer in DTDigiSimLinkCollection
      //   → look for a link on the same wire → check if the linked SimTrack is a muon
      std::map<unsigned int, int> votes;
      for (const DTRecHit1D* hit1D : hits1D) {
        DTLayerId layId = hit1D->wireId().layerId();
        int       wire  = hit1D->wireId().wire();

        for (const auto& linkUnit : *dtdigiLinks) {
          if (linkUnit.first.rawId() != layId.rawId()) continue;

          const DTDigiSimLinkCollection::Range& lRange = linkUnit.second;
          for (auto lnk = lRange.first; lnk != lRange.second; ++lnk) {
            if (lnk->wire() != wire) continue;
            unsigned int stId = lnk->SimTrackId();
            if (!simTrackToPdgId.count(stId)) continue; // not a primary muon
            votes[stId]++;
          }
        }
      }

      auto [winId, winPdg]      = bestSimMuon(votes);
      dtSeg_simTrackId_   [_nDTSeg] = winId;
      dtSeg_simTrack_pdgId_[_nDTSeg] = winPdg;

      if (debug) {
        std::cout << "DT Segment " << _nDTSeg
                  << ": MB" << cid.station() << "/W" << cid.wheel() << "/S" << cid.sector()
                  << "  gp=(" << gp.x() << "," << gp.y() << "," << gp.z() << ")"
                  << "  nPhi=" << dtSeg_nHits_phi_[_nDTSeg]
                  << "  nZ="   << dtSeg_nHits_z_  [_nDTSeg]
                  << "  chi2/ndof=" << dtSeg_chi2_[_nDTSeg] << "/" << dtSeg_ndof_[_nDTSeg]
                  << "  → SimTrack " << winId
                  << std::endl;
      }

      ++_nDTSeg;
      if (_nDTSeg >= MAX_DT_SEG) {
        edm::LogWarning("MuonHLTNtupler") << "MAX_DT_SEG reached, truncating.";
        break;
      }
    }
  }
  nDTSegments_ = _nDTSeg;
  if (debug) std::cout << "DT segments stored: " << nDTSegments_ << "\n";

  // ===========================================================================
  // 3.  CSC Segments
  // ===========================================================================

  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByToken(t_cscSegments_, cscSegments);

  int _nCSCSeg = 0;

  if (cscSegments.isValid()) {
    for (auto segIt = cscSegments->begin(); segIt != cscSegments->end(); ++segIt) {

      const CSCSegment& seg = *segIt;
      CSCDetId          cid = seg.cscDetId();
      const GeomDet*   gdet = cscGeom->idToDet(cid);

      GlobalPoint  gp = gdet->toGlobal(seg.localPosition());
      GlobalVector gv = gdet->toGlobal(seg.localDirection());

      // Geometry
      cscSeg_x_  [_nCSCSeg] = gp.x();
      cscSeg_y_  [_nCSCSeg] = gp.y();
      cscSeg_z_  [_nCSCSeg] = gp.z();
      cscSeg_dx_ [_nCSCSeg] = gv.x();
      cscSeg_dy_ [_nCSCSeg] = gv.y();
      cscSeg_dz_ [_nCSCSeg] = gv.z();
      cscSeg_eta_[_nCSCSeg] = gp.eta();
      cscSeg_phi_[_nCSCSeg] = gp.phi();

      // DetId
      cscSeg_endcap_ [_nCSCSeg] = cid.endcap();
      cscSeg_station_[_nCSCSeg] = cid.station();
      cscSeg_ring_   [_nCSCSeg] = cid.ring();
      cscSeg_chamber_[_nCSCSeg] = cid.chamber();

      // Quality
      cscSeg_nHits_[_nCSCSeg] = seg.nRecHits();
      cscSeg_chi2_ [_nCSCSeg] = static_cast<float>(seg.chi2());
      cscSeg_ndof_ [_nCSCSeg] = seg.degreesOfFreedom();

      // SimMuon matching via CSC strip digi–sim links
      // Each CSCRecHit2D carries the list of strip channels it was built from.
      // StripDigiSimLink.channel() gives the strip number; we check overlap.
      std::map<unsigned int, int> votes;

      if (cscStripDigiSimLink.isValid()) {
        for (const CSCRecHit2D& rh : seg.specificRecHits()) {
          CSCDetId rhId = rh.cscDetId();

          for (const auto& linkSet : *cscStripDigiSimLink) {
            if (linkSet.detId() != rhId.rawId()) continue;
            // if (linkSet.id.rawId() != rhId.rawId()) continue;

            for (const StripDigiSimLink& lnk : linkSet.data) {
              int linkStrip = static_cast<int>(lnk.channel());

              // Check whether this strip was part of the RecHit cluster
              bool stripMatched = false;
              for (unsigned int is = 0; is < rh.nStrips(); ++is) {
                if (rh.channels(is) == linkStrip) { stripMatched = true; break; }
              }
              if (!stripMatched) continue;

              unsigned int stId = lnk.SimTrackId();
              if (!simTrackToPdgId.count(stId)) continue;
              votes[stId]++;
            }
          }
        }
      }

      auto [winId, winPdg]       = bestSimMuon(votes);
      cscSeg_simTrackId_   [_nCSCSeg] = winId;
      cscSeg_simTrack_pdgId_[_nCSCSeg] = winPdg;

      if (debug) {
        std::cout << "CSC Segment " << _nCSCSeg
                  << ": ME" << cid.station() << "/" << cid.ring()
                  << "  endcap=" << cid.endcap()
                  << "  chamber=" << cid.chamber()
                  << "  gp=(" << gp.x() << "," << gp.y() << "," << gp.z() << ")"
                  << "  nHits=" << seg.nRecHits()
                  << "  chi2/ndof=" << cscSeg_chi2_[_nCSCSeg] << "/" << cscSeg_ndof_[_nCSCSeg]
                  << "  → SimTrack " << winId
                  << std::endl;
      }

      ++_nCSCSeg;
      if (_nCSCSeg >= MAX_CSC_SEG) {
        edm::LogWarning("MuonHLTNtupler") << "MAX_CSC_SEG reached, truncating.";
        break;
      }
    }
  }
  nCSCSegments_ = _nCSCSeg;
  if (debug) std::cout << "CSC segments stored: " << nCSCSegments_ << "\n";

  // ===========================================================================
  // 4.  GEM Segments
  // ===========================================================================

  edm::Handle<GEMSegmentCollection> gemSegments;
  iEvent.getByToken(t_gemSegments_, gemSegments);

  int _nGEMSeg = 0;

  if (gemSegments.isValid()) {
    for (auto segIt = gemSegments->begin(); segIt != gemSegments->end(); ++segIt) {

      const GEMSegment& seg = *segIt;
      GEMDetId          gid = seg.gemDetId();
      const GeomDet*   gdet = gemGeom->idToDet(gid);

      GlobalPoint  gp = gdet->toGlobal(seg.localPosition());
      GlobalVector gv = gdet->toGlobal(seg.localDirection());

      // Geometry
      gemSeg_x_  [_nGEMSeg] = gp.x();
      gemSeg_y_  [_nGEMSeg] = gp.y();
      gemSeg_z_  [_nGEMSeg] = gp.z();
      gemSeg_dx_ [_nGEMSeg] = gv.x();
      gemSeg_dy_ [_nGEMSeg] = gv.y();
      gemSeg_dz_ [_nGEMSeg] = gv.z();
      gemSeg_eta_[_nGEMSeg] = gp.eta();
      gemSeg_phi_[_nGEMSeg] = gp.phi();

      // DetId
      gemSeg_region_ [_nGEMSeg] = gid.region();
      gemSeg_station_[_nGEMSeg] = gid.station();
      gemSeg_ring_   [_nGEMSeg] = gid.ring();
      gemSeg_chamber_[_nGEMSeg] = gid.chamber();

      // Quality
      gemSeg_nHits_[_nGEMSeg] = seg.nRecHits();
      gemSeg_chi2_ [_nGEMSeg] = static_cast<float>(seg.chi2());
      gemSeg_ndof_ [_nGEMSeg] = seg.degreesOfFreedom();

      // SimMuon matching via GEM digi–sim links
      // GEMRecHit stores first strip of the cluster and cluster size;
      // any strip in [firstClusterStrip, firstClusterStrip+clusterSize) is a match.
      std::map<unsigned int, int> votes;

      if (gemDigiSimLinks.isValid()) {
        for (const GEMRecHit& rh : seg.specificRecHits()) {
          GEMDetId rhId      = rh.gemId();
          int firstStrip     = rh.firstClusterStrip();
          int clusterSize    = rh.clusterSize();

          for (const auto& linkSet : *gemDigiSimLinks) {
            if (linkSet.detId() != rhId.rawId()) continue;
            // if (linkSet.id.rawId() != rhId.rawId()) continue;

            for (const GEMDigiSimLink& lnk : linkSet.data) {
              int strip = static_cast<int>(lnk.getStrip());
              if (strip < firstStrip || strip >= firstStrip + clusterSize) continue;

              unsigned int stId = lnk.getTrackId();
              if (!simTrackToPdgId.count(stId)) continue;
              votes[stId]++;
            }
          }
        }
      }

      auto [winId, winPdg]       = bestSimMuon(votes);
      gemSeg_simTrackId_   [_nGEMSeg] = winId;
      gemSeg_simTrack_pdgId_[_nGEMSeg] = winPdg;

      if (debug) {
        std::cout << "GEM Segment " << _nGEMSeg
                  << ": region=" << gid.region()
                  << "  station=" << gid.station()
                  << "  ring=" << gid.ring()
                  << "  chamber=" << gid.chamber()
                  << "  gp=(" << gp.x() << "," << gp.y() << "," << gp.z() << ")"
                  << "  nHits=" << seg.nRecHits()
                  << "  → SimTrack " << winId
                  << std::endl;
      }

      ++_nGEMSeg;
      if (_nGEMSeg >= MAX_GEM_SEG) {
        edm::LogWarning("MuonHLTNtupler") << "MAX_GEM_SEG reached, truncating.";
        break;
      }
    }
  }
  nGEMSegments_ = _nGEMSeg;
  if (debug) std::cout << "GEM segments stored: " << nGEMSegments_ << "\n";

  // ===========================================================================
  // 5.  ME0 Segments
  // ===========================================================================

  edm::Handle<ME0SegmentCollection> me0Segments;
  iEvent.getByToken(t_me0Segments_, me0Segments);

  int _nME0Seg = 0;

  if (me0Segments.isValid()) {
    for (auto segIt = me0Segments->begin(); segIt != me0Segments->end(); ++segIt) {

      const ME0Segment& seg = *segIt;
      ME0DetId          mid = seg.me0DetId();
      const GeomDet*   gdet = me0Geom->idToDet(mid);

      GlobalPoint  gp = gdet->toGlobal(seg.localPosition());
      GlobalVector gv = gdet->toGlobal(seg.localDirection());

      // Geometry
      me0Seg_x_  [_nME0Seg] = gp.x();
      me0Seg_y_  [_nME0Seg] = gp.y();
      me0Seg_z_  [_nME0Seg] = gp.z();
      me0Seg_dx_ [_nME0Seg] = gv.x();
      me0Seg_dy_ [_nME0Seg] = gv.y();
      me0Seg_dz_ [_nME0Seg] = gv.z();
      me0Seg_eta_[_nME0Seg] = gp.eta();
      me0Seg_phi_[_nME0Seg] = gp.phi();

      // DetId
      me0Seg_region_ [_nME0Seg] = mid.region();
      me0Seg_station_[_nME0Seg] = mid.station();
      me0Seg_chamber_[_nME0Seg] = mid.chamber();

      // Quality
      me0Seg_nHits_[_nME0Seg] = seg.nRecHits();
      me0Seg_chi2_ [_nME0Seg] = static_cast<float>(seg.chi2());
      me0Seg_ndof_ [_nME0Seg] = seg.degreesOfFreedom();

      // SimMuon matching via ME0 digi–sim links
      // ME0DigiSimLink has strip() and SimTrackId(), same pattern as GEM.
      std::map<unsigned int, int> votes;

      if (me0DigiSimLinks.isValid()) {
        for (const ME0RecHit& rh : seg.specificRecHits()) {
          ME0DetId rhId   = rh.me0Id();

          for (const auto& linkSet : *me0DigiSimLinks) {
            if (linkSet.detId() != rhId.rawId()) continue;

            for (const ME0DigiSimLink& lnk : linkSet.data) {

              unsigned int stId = lnk.getTrackId();
              if (!simTrackToPdgId.count(stId)) continue;
              votes[stId]++;
            }
          }
        }
      }

      auto [winId, winPdg]        = bestSimMuon(votes);
      me0Seg_simTrackId_   [_nME0Seg] = winId;
      me0Seg_simTrack_pdgId_[_nME0Seg] = winPdg;

      if (debug) {
        std::cout << "ME0 Segment " << _nME0Seg
                  << ": region=" << mid.region()
                  << "  station=" << mid.station()
                  << "  chamber=" << mid.chamber()
                  << "  gp=(" << gp.x() << "," << gp.y() << "," << gp.z() << ")"
                  << "  nHits=" << seg.nRecHits()
                  << "  → SimTrack " << winId
                  << std::endl;
      }

      ++_nME0Seg;
      if (_nME0Seg >= MAX_ME0_SEG) {
        edm::LogWarning("MuonHLTNtupler") << "MAX_ME0_SEG reached, truncating.";
        break;
      }
    }
  }
  nME0Segments_ = _nME0Seg;
  if (debug) std::cout << "ME0 segments stored: " << nME0Segments_ << "\n";

} // end Fill_Tracks_And_Segments


void MuonHLTNtupler::endJob() {}
void MuonHLTNtupler::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}
void MuonHLTNtupler::endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}

DEFINE_FWK_MODULE(MuonHLTNtupler);
