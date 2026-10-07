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
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/ESInputTag.h"


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
#include "Geometry/CommonTopologies/interface/GeomDetType.h"
#include "Geometry/CommonTopologies/interface/GeomDet.h"

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
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"
#include "DataFormats/Math/interface/angle_units.h"

#include "RecoMuon/TrackerSeedGenerator/interface/SeedMvaEstimatorPhase2.h"

// -- for L1TkMu propagation
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTrackerBuilder.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "MuonHLTTool/MuonHLTNtupler/interface/MuonHLTobjCorrelator.h"
#include "MuonHLTTool/MuonHLTNtupler/interface/GenericMuonCollections.h"
#include "MuonHLTTool/MuonHLTNtupler/interface/MuonHLTTriggerMatching.h"

#include "TTree.h"
#include "TString.h"

using namespace std;
using namespace reco;
using namespace edm;

class MuonHLTNtupler : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
public:
  MuonHLTNtupler(const edm::ParameterSet &iConfig);
  virtual ~MuonHLTNtupler() {};

  virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) override;
  virtual void endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) override;

private:
  void Init();
  void Make_Branch();
  void Fill_HLT(const edm::Event &iEvent, bool isMYHLT);
  void Fill_Muon(const edm::Event &iEvent);
  void Fill_GenericMuonCollections(const edm::Event &iEvent);
  void Fill_GenParticle(const edm::Event &iEvent);
  
  bool SavedTriggerCondition( std::string& pathName );
  bool SavedFilterCondition( std::string& filterName );

  bool isNewHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx);

  bool isMuonCollectionValid(const std::string& name) const;
  
  bool DebugMode;
  
  // edm::EDGetTokenT< std::vector<reco::Muon> >                t_offlineMuon_;
  edm::EDGetTokenT< edm::View<reco::Muon> >                  t_offlineMuon_;
  edm::EDGetTokenT< reco::VertexCollection >                 t_offlineVertex_;
  edm::EDGetTokenT< edm::TriggerResults >                    t_triggerResults_;
  edm::EDGetTokenT< trigger::TriggerEvent >                  t_triggerEvent_;
  edm::EDGetTokenT< edm::TriggerResults >                    t_myTriggerResults_;
  edm::EDGetTokenT< trigger::TriggerEvent >                  t_myTriggerEvent_;
 
  // -- L2Muon, L3Muon, TkMuon, iterL3OI, iterL3IOFromL2, iterL3FromL2 used
  //    to each have a dedicated token + dedicated Fill_* code here. They are
  //    now handled generically -- see muonCollectionCfgs_/muonCollectionAdapters_
  //    below and GenericMuonCollections.h. Adding, removing, or renaming one
  //    of these (or any future HLT muon collection) is now a python-only
  //    change; L1Muon/L1TkMuon/iterL3Muon/iterL3MuonNoID/iterL3IOFromL1 stay
  //    as dedicated members above because they carry additional structure
  //    (hardware words, ID flags, seed-matching side effects) that a plain
  //    pt/eta/phi/charge record doesn't capture -- see REFACTOR_NOTES.md.
  std::vector<MuonCollectionCfg> muonCollectionCfgs_;
  std::vector<std::unique_ptr<IMuonCollectionAdapter>> muonCollectionAdapters_;
  // -- names of muonCollectionCfgs_, in order; written to the tree as its
  //    own branch (see Make_Branch()) so downstream analysis code can
  //    discover which collections exist in a given file instead of
  //    hardcoding a list, the same way configuredTriggerFilters/Paths
  //    already let it discover filters/paths.
  std::vector<std::string> muonCollectionNamesForBranch_;
  // -- one entry per muonCollectionCfgs_; set after the first time a given
  //    collection is found invalid, so the LogWarning in
  //    Fill_GenericMuonCollections() fires once per job (per collection)
  //    instead of once per event.
  mutable std::vector<bool> warnedCollectionInvalid_;

  std::vector<std::string> triggerPathsCfg_;       // unversioned, from python
  std::vector<std::string> triggerFiltersCfg_;      // from python
  std::vector<std::string> resolvedTriggerPaths_;   // versioned, filled in beginRun
  std::vector<std::string> resolvedPathLastFilter_; // one per resolvedTriggerPaths_ entry,
                                                     // filled in beginRun -- see lastFilterOfPath()
  double maxDR_;                                    // trigger-matching cone size
  HLTConfigProvider hltConfig_;

  // -- one entry per name in muonCollectionCfgs_; std::map node addresses
  //    are stable across insertion, so &genPt_[name] taken once in
  //    Make_Branch() stays valid after later insertions into the map.
  std::map<std::string, std::vector<float>> muonTrack_Pt_, muonTrack_Eta_, muonTrack_Phi_, muonTrack_Charge_, muonTrack_TrackPt_;
  // -- only populated (non -99) for collections of type "MuonTrackLinks"
  //    (iterL3OI, iterL3IOFromL2, iterL3FromL2 by default); -99 elsewhere.
  //    Preserves the old inner/outer/global branch triplets without a
  //    dedicated array per collection.
  std::map<std::string, std::vector<float>> muonTrack_InnerPt_, muonTrack_InnerEta_, muonTrack_InnerPhi_, muonTrack_InnerCharge_;
  std::map<std::string, std::vector<float>> muonTrack_OuterPt_, muonTrack_OuterEta_, muonTrack_OuterPhi_, muonTrack_OuterCharge_;
  std::map<std::string, std::vector<float>> muonTrack_GlobalPt_, muonTrack_GlobalEta_, muonTrack_GlobalPhi_, muonTrack_GlobalCharge_;
  // -- one INNER vector per muon (same indexing as genPt_[name], etc.):
  //    genTrigMatchedFilterIdx_[name][i] lists the indices (into
  //    triggerFiltersCfg_/"configuredTriggerFilters") of every configured
  //    filter that muon i matched within maxDR_; genTrigDR_[name][i] is
  //    the parallel list of dR values. Usually 0 or 1 entries per muon for
  //    an exclusive single-object trigger, but can be more.
  std::map<std::string, std::vector<std::vector<int>>>   muonTrack_TrigMatchedFilterIdx_;
  std::map<std::string, std::vector<std::vector<float>>> muonTrack_TrigDR_;

  // -- same idea, but per configured PATH (matched against that path's
  //    last filter -- see resolvedPathLastFilter_) rather than per
  //    individually-configured filter. Indices refer to
  //    triggerPathsCfg_/"configuredTriggerPaths", same ordering.
  std::map<std::string, std::vector<std::vector<int>>>   muonTrack_PathMatchedIdx_;
  std::map<std::string, std::vector<std::vector<float>>> muonTrack_PathDR_;

  edm::EDGetTokenT< reco::VertexCollection >                 t_hltIterL3MuonTrimmedPixelVertices_;
  edm::EDGetTokenT< reco::VertexCollection >                 t_hltIterL3FromL1MuonTrimmedPixelVertices_;

  edm::EDGetTokenT< LumiScalersCollection >                  t_lumiScaler_;
  edm::EDGetTokenT< LumiScalersCollection >                  t_offlineLumiScaler_;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> >         t_PUSummaryInfo_;
  edm::EDGetTokenT< GenEventInfoProduct >                    t_genEventInfo_;
  edm::EDGetTokenT< reco::GenParticleCollection >            t_genParticle_;

  // ------ For CMSSW_12 -------
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyESToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryESToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldESToken_;
  const edm::ESGetToken<GeometricDet, IdealGeometryRecord> geomDetESToken_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorESToken_;

  const TrackerTopology* _trackerTopology;
  const TrackerGeometry* _trackerGeometry;
  // const GeometricDet* _geometryDetector;
  // const MagneticField* _magneticField;
  // --------------------------------------

  TTree *ntuple_;
  static const int arrSize_ = 5000;

  // -- general event information
  bool isRealData_;
  int runNum_;
  int lumiBlockNum_;
  unsigned long long eventNum_;

  int nVertex_;

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
  double qScale_;

  vector<float> PU_pT_hats_;

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

  int muon_isLooseTriggerMuon_[arrSize_];
  int muon_isME0Muon_[arrSize_];
  int muon_isGEMMuon_[arrSize_];
  int muon_isRPCMuon_[arrSize_];
  int muon_isGoodMuon_TMOneStationTight_[arrSize_];

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

  double muon_normChi2_global_[arrSize_];
  int muon_nTrackerHit_global_[arrSize_];
  int muon_nTrackerLayer_global_[arrSize_];
  int muon_nPixelHit_global_[arrSize_];
  int muon_nMuonHit_global_[arrSize_];

  double muon_normChi2_inner_[arrSize_];
  int muon_nTrackerHit_inner_[arrSize_];
  int muon_nTrackerLayer_inner_[arrSize_];
  int muon_nPixelHit_inner_[arrSize_];

  double muon_pt_tuneP_[arrSize_];
  double muon_ptError_tuneP_[arrSize_];

  double muon_dxyVTX_best_[arrSize_];
  double muon_dzVTX_best_[arrSize_];

  int muon_nMatchedStation_[arrSize_];
  int muon_nMatchedRPCLayer_[arrSize_];
  int muon_stationMask_[arrSize_];
  int muon_expectedNnumberOfMatchedStations_[arrSize_];

  class vtxTemplate {
  private:
    int nVtxs;
    std::vector<int> isValid;
    std::vector<double> chi2;
    std::vector<double> ndof;
    std::vector<double> nTracks;
    std::vector<double> vtxX;
    std::vector<double> vtxXerr;
    std::vector<double> vtxY;
    std::vector<double> vtxYerr;
    std::vector<double> vtxZ;
    std::vector<double> vtxZerr;
  public:
    void clear() {
      nVtxs = 0;
      isValid.clear();
      chi2.clear();
      ndof.clear();
      nTracks.clear();
      vtxX.clear();
      vtxXerr.clear();
      vtxY.clear();
      vtxYerr.clear();
      vtxZ.clear();
      vtxZerr.clear();

      return;
    }

    void setBranch(TTree* tmpntpl, TString name) {
      tmpntpl->Branch("n"+name, &nVtxs);
      tmpntpl->Branch(name+"_isValid", &isValid);
      tmpntpl->Branch(name+"_chi2", &chi2);
      tmpntpl->Branch(name+"_ndof", &ndof);
      tmpntpl->Branch(name+"_nTracks", &nTracks);
      tmpntpl->Branch(name+"_x", &vtxX);
      tmpntpl->Branch(name+"_xerr", &vtxXerr);
      tmpntpl->Branch(name+"_y", &vtxY);
      tmpntpl->Branch(name+"_yerr", &vtxYerr);
      tmpntpl->Branch(name+"_z", &vtxZ);
      tmpntpl->Branch(name+"_zerr", &vtxZerr);

      return;
    }

    void fill(const reco::Vertex vtx) {
      isValid.push_back(vtx.isValid());
      chi2.push_back(vtx.chi2());
      ndof.push_back(vtx.ndof());
      nTracks.push_back(vtx.nTracks());
      vtxX.push_back(vtx.x());
      vtxXerr.push_back(vtx.xError());
      vtxY.push_back(vtx.y());
      vtxYerr.push_back(vtx.yError());
      vtxZ.push_back(vtx.z());
      vtxZerr.push_back(vtx.zError());
      nVtxs++;

      return;
    }
  };

  vtxTemplate* VThltIterL3MuonTrimmedPixelVertices = new vtxTemplate();
  vtxTemplate* VThltIterL3FromL1MuonTrimmedPixelVertices = new vtxTemplate();

};
