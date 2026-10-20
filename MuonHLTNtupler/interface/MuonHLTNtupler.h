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

#include "DataFormats/DTRecHit/interface/DTRecHit1D.h"
#include "DataFormats/DTRecHit/interface/DTRecHit1DPair.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <span>

#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"
#include "PhysicsTools/TruthInfo/interface/BranchHitAssociator.h"

// -- muon segments (DT / CSC / GEM / ME0) and their geometry
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
// -- muon rechits (DT / CSC / RPC / GEM / ME0)
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0RecHitCollection.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "TTree.h"
#include "TString.h"

using namespace std;
using namespace reco;
using namespace edm;

namespace {
  // Ordinary hadron containing quark `flavor` (5 = b, 4 = c).
  inline bool hadronHasQuark(int pdgId, int flavor) {
    const int id = std::abs(pdgId);
    if (id < 100 || id >= 1000000000)
      return false;
    return (id / 1000) % 10 == flavor || (id / 100) % 10 == flavor || (id / 10) % 10 == flavor;
  }
 
  // Bring a muon-system reco hit into the DetId space of the stored PSimHit
  // detUnitId(), so the shared-DetId match works.
  //   DT : PSimHits are stored per WIRE (DTWireId); DTRecHit1D lives on the
  //        layer but remembers its wire -> use wireId().
  //   CSC: sim hits of the whole ME1/1 chamber are ring 1, while ME1/1a reco
  //        hits are ring 4 -> fold ring 4 into ring 1.
  //   RPC/GEM: roll / eta-partition on both sides already.
  inline std::optional<uint32_t> muonHitDetId(TrackingRecHit const& hit) {
    const DetId id = hit.geographicalId();
    if (id.det() != DetId::Muon)
      return std::nullopt;
    switch (id.subdetId()) {
      case MuonSubdetId::DT:
        if (auto const* h = dynamic_cast<DTRecHit1D const*>(&hit))
          return h->wireId().rawId();
        if (auto const* h = dynamic_cast<DTRecHit1DPair const*>(&hit))
          return h->wireId().rawId();
        return std::nullopt;  // DT segments (4D/2D) have no single wire: skip
      case MuonSubdetId::CSC: {
        CSCDetId cscId(id);
        if (cscId.station() == 1 && cscId.ring() == 4)
          cscId = CSCDetId(cscId.endcap(), 1, 1, cscId.chamber(), cscId.layer());
        return cscId.rawId();
      }
      default:
        return id.rawId();
    }
  }
  inline void sortUnique(std::vector<truth::RecoHit>& hits) {
    std::sort(hits.begin(), hits.end(), [](auto const& a, auto const& b) { return a.detId < b.detId; });
    hits.erase(std::unique(hits.begin(), hits.end(), [](auto const& a, auto const& b) { return a.detId == b.detId; }),
               hits.end());
  }
  // Sorted, unique tracker-channel DetIds of each candidate root's subgraph.
  inline std::vector<std::vector<uint32_t>> rootDetIds(truth::LogicalGraphHitIndex const& idx,
                                                       std::vector<uint32_t> const& roots,
                                                       truth::HitChannel ch) {
    std::vector<std::vector<uint32_t>> out(roots.size());
    std::vector<truth::LogicalGraphHitIndex::Hit> tmp;
    for (std::size_t r = 0; r < roots.size(); ++r) {
      tmp.clear();
      idx.appendSubgraphHits(ch, roots[r], tmp);
      auto& ids = out[r];
      ids.reserve(tmp.size());
      for (auto const& h : tmp) ids.push_back(h.detId);
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    }
    return out;
  }

  // Best root by number of shared DetIds; recoHits must be sorted+unique.
  // Same convention as BranchHitAssociator's SharedHits metric.
  inline std::optional<truth::BranchMatch> matchByDetIds(std::vector<truth::RecoHit> const& recoHits,
                                                         std::vector<uint32_t> const& roots,
                                                         std::vector<std::vector<uint32_t>> const& rootIds,
                                                         double maxScore) {
    if (recoHits.empty()) return std::nullopt;
    std::optional<truth::BranchMatch> best;
    for (std::size_t r = 0; r < roots.size(); ++r) {
      auto const& ids = rootIds[r];
      uint32_t shared = 0;
      for (auto const& h : recoHits)
        shared += std::binary_search(ids.begin(), ids.end(), h.detId);
      if (shared == 0) continue;
      const float score = 1.f - float(shared) / float(recoHits.size());
      if (!best || score < best->score || (score == best->score && shared > best->sharedEnergy)) {
        truth::BranchMatch m;
        m.rootParticleId = roots[r];
        m.sharedEnergy = float(shared);
        m.score = score;
        m.reverseScore = ids.empty() ? 1.f : 1.f - float(shared) / float(ids.size());
        best = m;
      }
    }
    if (best && best->score > maxScore) return std::nullopt;
    return best;
  }
}  // namespace

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
  void Fill_Muon(const edm::Event &iEvent, bool isData);
  void Fill_GenericMuonCollections(const edm::Event &iEvent);
  void Fill_GenParticle(const edm::Event &iEvent);
  void Fill_MuonSystem(const edm::Event &iEvent, const edm::EventSetup &iSetup, bool isData);

  std::optional<truth::BranchMatch> bestMatch(truth::BranchHitAssociator const&, std::vector<truth::RecoHit>&) const;
  void fillTruth(truth::Graph const&, std::optional<truth::BranchMatch> const&, float nRecoHits, int _nMuon);
  
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

  // Truth-Graph
  edm::EDGetTokenT<truth::Graph>                             graphToken_;
  edm::EDGetTokenT<truth::LogicalGraphHitIndex>              hitIndexToken_;
  const double maxScore_;

  // Muon segments
  edm::EDGetTokenT<DTRecSegment4DCollection>                 t_dtSegments_;
  edm::EDGetTokenT<CSCSegmentCollection>                     t_cscSegments_;
  edm::EDGetTokenT<GEMSegmentCollection>                     t_gemSegments_;
  edm::EDGetTokenT<ME0SegmentCollection>                     t_me0Segments_;
  const bool runME0_;  // legacy ME0 path; in recent Phase-2 geometries ME0 = GE0 (GEM station 0)
  const edm::ESGetToken<DTGeometry,  MuonGeometryRecord>     dtGeomToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord>     cscGeomToken_;
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord>     gemGeomToken_;
  edm::ESGetToken<ME0Geometry, MuonGeometryRecord>           me0GeomToken_;  // only consumed if runME0_

  // Muon rechits
  edm::EDGetTokenT<DTRecHitCollection>                       t_dtRecHits_;
  edm::EDGetTokenT<CSCRecHit2DCollection>                    t_cscRecHits_;
  edm::EDGetTokenT<RPCRecHitCollection>                      t_rpcRecHits_;
  edm::EDGetTokenT<GEMRecHitCollection>                      t_gemRecHits_;
  edm::EDGetTokenT<ME0RecHitCollection>                      t_me0RecHits_;   // only read if runME0_
  const edm::ESGetToken<RPCGeometry, MuonGeometryRecord>     rpcGeomToken_;
 
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

  std::map<std::string, std::vector<int>>   muonTrack_OfflineMatchedIdx_;
  std::map<std::string, std::vector<float>> muonTrack_OfflineDR_;

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

  int muon_matchChannel_[arrSize_];
  int muon_muSysSimId_[arrSize_];
  double muon_muSysHitPurity_[arrSize_];
  int muon_matched_[arrSize_];
  int muon_simId_[arrSize_];
  int muon_simPdgId_[arrSize_];
  int muon_simHasGen_[arrSize_];
  int muon_simMotherPdgId_[arrSize_];
  double muon_simPt_[arrSize_];
  double muon_simEta_[arrSize_];
  double muon_simPhi_[arrSize_];
  int muon_sharedHits_[arrSize_];
  double muon_hitPurity_[arrSize_];
  int muon_fromW_[arrSize_]; 
  int muon_fromZ_[arrSize_];
  int muon_fromTau_[arrSize_];
  int muon_fromLightHadron_[arrSize_];
  int muon_fromB_[arrSize_];
  int muon_fromC_[arrSize_];
  int muon_fromBSM_[arrSize_];
  int muon_isPrompt_[arrSize_];
  int muon_truthHead_[arrSize_];

public:
  // ---------------------------------------------------------------------------
  // Truth result for one segment (filled by classifySegment in the .cc)
  //   origin:  0 no truth match          1 the muon itself
  //            2 muon delta ray (ionisation)  3 muon bremsstrahlung
  //            4 other EM process of the muon (pair production, ...)
  //            5 muon decay in flight      6 muon nuclear interaction
  //            7 not from a muon (punch-through, neutrons, ...)
  //   originReason: truth::VertexReason of the process behind `origin`
  //            (for origin 1: production of the muon chain head, e.g. Decay
  //            for a pi/K decay-in-flight muon; for 2-6: production of the
  //            first-generation secondary of the muon; for 7: production of
  //            the depositing particle).
  // ---------------------------------------------------------------------------
  struct SegTruth {
    int   truthId = -1, truthPdgId = 0, nTruthHits = 0;
    float purity = 0.f;
    int   origin = 0, originReason = 0;
    int   muonHeadId = -1, genMuonId = -1, offlineMuonIdx = -1;
    int   nCandidates = 0;  // distinct truth particles sharing the matched cells (>1 = ambiguous)
    float genMuonPt = -999.f, genMuonEta = -999.f, genMuonPhi = -999.f;
  };

private:
  // -- one instance per segment type; branches are <prefix>_<field>, counter n<Prefix>
  class segTemplate {
  public:
    std::vector<float> x, y, z, dx, dy, dz, eta, phi, chi2;
    std::vector<int>   ndof, nHits, region, station, ring, chamber;
    std::vector<int>   truthId, truthPdgId, nTruthHits, origin, originReason;
    std::vector<int>   muonHeadId, genMuonId, offlineMuonIdx, nCandidates;
    std::vector<float> purity, genMuonPt, genMuonEta, genMuonPhi;
    int n = 0;
 
    void clear() {
      n = 0;
      for (auto* v : {&x, &y, &z, &dx, &dy, &dz, &eta, &phi, &chi2, &purity, &genMuonPt, &genMuonEta, &genMuonPhi})
        v->clear();
      for (auto* v : {&ndof, &nHits, &region, &station, &ring, &chamber, &truthId, &truthPdgId, &nTruthHits,
                      &origin, &originReason, &muonHeadId, &genMuonId, &offlineMuonIdx, &nCandidates})
        v->clear();
    }
 
    void setBranch(TTree* t, const TString& p, const TString& counter) {
      t->Branch(counter, &n);
      t->Branch(p + "_nCandidates", &nCandidates);
      t->Branch(p + "_x", &x);           t->Branch(p + "_y", &y);           t->Branch(p + "_z", &z);
      t->Branch(p + "_dx", &dx);         t->Branch(p + "_dy", &dy);         t->Branch(p + "_dz", &dz);
      t->Branch(p + "_eta", &eta);       t->Branch(p + "_phi", &phi);
      t->Branch(p + "_chi2", &chi2);     t->Branch(p + "_ndof", &ndof);     t->Branch(p + "_nHits", &nHits);
      t->Branch(p + "_region", &region); t->Branch(p + "_station", &station);
      t->Branch(p + "_ring", &ring);     t->Branch(p + "_chamber", &chamber);
      t->Branch(p + "_truthId", &truthId);           t->Branch(p + "_truthPdgId", &truthPdgId);
      t->Branch(p + "_nTruthHits", &nTruthHits);     t->Branch(p + "_purity", &purity);
      t->Branch(p + "_origin", &origin);             t->Branch(p + "_originReason", &originReason);
      t->Branch(p + "_muonHeadId", &muonHeadId);     t->Branch(p + "_genMuonId", &genMuonId);
      t->Branch(p + "_genMuonPt", &genMuonPt);       t->Branch(p + "_genMuonEta", &genMuonEta);
      t->Branch(p + "_genMuonPhi", &genMuonPhi);     t->Branch(p + "_offlineMuonIdx", &offlineMuonIdx);
    }

    // region/station/ring/chamber are the detector's own indices:
    //   DT : region = 0,       ring = wheel,  chamber = sector
    //   CSC: region = endcap (1 = +z, 2 = -z)
    //   GEM/ME0: region = +-1
    void fill(const GlobalPoint& gp, const GlobalVector& gv, float c2, int nd, int nh,
              int reg, int st, int rg, int ch, const SegTruth& t) {
      x.push_back(gp.x());   y.push_back(gp.y());   z.push_back(gp.z());
      dx.push_back(gv.x());  dy.push_back(gv.y());  dz.push_back(gv.z());
      eta.push_back(gp.eta()); phi.push_back(gp.phi());
      chi2.push_back(c2); ndof.push_back(nd); nHits.push_back(nh);
      region.push_back(reg); station.push_back(st); ring.push_back(rg); chamber.push_back(ch);
      truthId.push_back(t.truthId);       truthPdgId.push_back(t.truthPdgId);
      nTruthHits.push_back(t.nTruthHits); purity.push_back(t.purity);
      origin.push_back(t.origin);         originReason.push_back(t.originReason);
      muonHeadId.push_back(t.muonHeadId); genMuonId.push_back(t.genMuonId);
      genMuonPt.push_back(t.genMuonPt);   genMuonEta.push_back(t.genMuonEta); genMuonPhi.push_back(t.genMuonPhi);
      offlineMuonIdx.push_back(t.offlineMuonIdx);
      nCandidates.push_back(t.nCandidates);
      ++n;
    }
  };

  segTemplate dtSegs_, cscSegs_, gemSegs_, me0Segs_;
  std::vector<int> dtSeg_nHitsPhi_, dtSeg_nHitsZ_;  // DT-specific

  // -- one instance per rechit type; branches are <prefix>_<field>, counter n<Prefix>
  //    region/station/ring/sector/layer/roll are the detector's own indices:
  //      DT : region 0, ring = wheel, sector, layer = 4*(superlayer-1)+layer (1..12), roll = superlayer
  //      CSC: region = endcap, ring, sector = chamber, layer, roll = 0
  //      RPC: region (0 barrel, +-1 endcap), ring (= wheel in barrel), station, sector, layer, roll
  //      GEM/ME0: region, ring, sector = chamber, layer, roll = eta partition
  //    channel/clusterSize: DT wire / 1, CSC first strip / nStrips (wireGroup in `wire`),
  //                         RPC & GEM first strip / cluster size, ME0 -1 / -1
  //    time: DT digi time, CSC tpeak, RPC time, GEM bunch crossing, ME0 tof
  //    segIdx: index into the matching <det>Seg_* arrays of the segment that
  //            uses this hit (-1 if none; RPC has no segments)
  class rhTemplate {
  public:
    std::vector<float> x, y, z, eta, phi, time;
    std::vector<int>   region, station, ring, sector, layer, roll, channel, clusterSize, wire, segIdx;
    std::vector<int>   truthId, truthPdgId, nCandidates, origin, originReason;
    std::vector<int>   muonHeadId, genMuonId, offlineMuonIdx;
    int n = 0;
 
    void clear() {
      n = 0;
      for (auto* v : {&x, &y, &z, &eta, &phi, &time})
        v->clear();
      for (auto* v : {&region, &station, &ring, &sector, &layer, &roll, &channel, &clusterSize, &wire, &segIdx,
                      &truthId, &truthPdgId, &nCandidates, &origin, &originReason, &muonHeadId, &genMuonId,
                      &offlineMuonIdx})
        v->clear();
    }
 
    void setBranch(TTree* t, const TString& p, const TString& counter) {
      t->Branch(counter, &n);
      t->Branch(p + "_x", &x);     t->Branch(p + "_y", &y);       t->Branch(p + "_z", &z);
      t->Branch(p + "_eta", &eta); t->Branch(p + "_phi", &phi);   t->Branch(p + "_time", &time);
      t->Branch(p + "_region", &region);   t->Branch(p + "_station", &station); t->Branch(p + "_ring", &ring);
      t->Branch(p + "_sector", &sector);   t->Branch(p + "_layer", &layer);     t->Branch(p + "_roll", &roll);
      t->Branch(p + "_channel", &channel); t->Branch(p + "_clusterSize", &clusterSize);
      t->Branch(p + "_wire", &wire);       t->Branch(p + "_segIdx", &segIdx);
      t->Branch(p + "_truthId", &truthId);         t->Branch(p + "_truthPdgId", &truthPdgId);
      t->Branch(p + "_nCandidates", &nCandidates); t->Branch(p + "_origin", &origin);
      t->Branch(p + "_originReason", &originReason);
      t->Branch(p + "_muonHeadId", &muonHeadId);   t->Branch(p + "_genMuonId", &genMuonId);
      t->Branch(p + "_offlineMuonIdx", &offlineMuonIdx);
    }
 
    struct Ids {
      int region = 0, station = 0, ring = 0, sector = 0, layer = 0, roll = 0;
      int channel = -1, clusterSize = -1, wire = -1;
    };

    void fill(const GlobalPoint& gp, float t, const Ids& id, int seg, const SegTruth& tr) {
      x.push_back(gp.x()); y.push_back(gp.y()); z.push_back(gp.z());
      eta.push_back(gp.eta()); phi.push_back(gp.phi()); time.push_back(t);
      region.push_back(id.region); station.push_back(id.station); ring.push_back(id.ring);
      sector.push_back(id.sector); layer.push_back(id.layer);     roll.push_back(id.roll);
      channel.push_back(id.channel); clusterSize.push_back(id.clusterSize); wire.push_back(id.wire);
      segIdx.push_back(seg);
      truthId.push_back(tr.truthId);         truthPdgId.push_back(tr.truthPdgId);
      nCandidates.push_back(tr.nCandidates); origin.push_back(tr.origin);
      originReason.push_back(tr.originReason);
      muonHeadId.push_back(tr.muonHeadId);   genMuonId.push_back(tr.genMuonId);
      offlineMuonIdx.push_back(tr.offlineMuonIdx);
      ++n;
    }
  };

  rhTemplate dtRecHits_, cscRecHits_, rpcRecHits_, gemRecHits_, me0RecHits_;
  std::vector<float> dtRh_xLeft_, dtRh_xRight_;  // DT-specific: local x of the two L/R solutions (layer frame)



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
