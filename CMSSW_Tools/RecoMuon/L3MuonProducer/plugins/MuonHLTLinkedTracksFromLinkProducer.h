#ifndef MuonHLTSegmentMatcher_MuonHLTSegmentMatcherProducer_MuonHLTMuonFromLinksProducer_h
#define MuonHLTSegmentMatcher_MuonHLTSegmentMatcherProducer_MuonHLTMuonFromLinksProducer_h
/*
 * MuonHLTMuonFromLinksProducer.h
 * ================================
 * Builds TrackerMuons, GlobalMuons and StandAloneMuons from the two link
 * collections written by MuonHLTSegmentMatcherProducer (a GNN that has
 * already solved track<->segment and segment<->segment association).
 *
 * - TrackerMuon : tracker track decorated with chamber/segment matches
 *                 (no independent muon-system fit).
 * - GlobalMuon  : tracker track refit together with the GNN-selected muon
 *                 segments in a single Kalman filter (GlobalMuonRefitter),
 *                 exactly as GlobalTrajectoryBuilderBase::build() does for
 *                 offline/HLT L3 muons.
 * - StandAloneMuon : segment clusters that the GNN linked to EACH OTHER
 *                 (segment-segment links) but NOT to any tracker track are
 *                 fit on their own, seeded from a simple two-point circle
 *                 estimate through the segment positions (see
 *                 estimateChargeAndMomentum below), then refined by the
 *                 same generic Kalman fitter/smoother used offline
 *                 (KFFitterSmootherSTA).
 *
 * NOTE ON VERIFICATION: everything routed through MuonServiceProxy,
 * GlobalMuonRefitter and MuonTrackLoader is exercised against the actual
 * cms-sw/cmssw headers/source (checked against `master` while writing
 * this). The standalone seed estimate (estimateChargeAndMomentum /
 * buildSeedState) is NEW code written for this producer, not lifted from
 * an existing CMSSW class -- treat it as a first cut to tune, not as
 * something copied from a validated CMS algorithm.
 */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
//#include "TrackingTools/PatternTools/interface/TrajectoryFwd.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/TrackingTools/interface/MuonTrackLoader.h"
#include "RecoMuon/TrackingTools/interface/MuonCandidate.h"
#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonRefitter.h"

#include "RecoMuon/L3MuonProducer/interface/MuonSegmentLink.h"

#include <vector>
#include <map>
#include <memory>

class Trajectory;

class MuonHLTLinkedTracksFromLinkProducer : public edm::stream::EDProducer<> {
public:
  explicit MuonHLTLinkedTracksFromLinkProducer(const edm::ParameterSet&);
  ~MuonHLTLinkedTracksFromLinkProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  // ---- internal record type (defined in the .cc, forward-declared here
  //      via an opaque struct is avoided for simplicity: see .cc) ----
  struct SegRecord;

  // helper: build the tracker-side TransientTrackingRecHits for a track,
  // mirroring GlobalTrajectoryBuilderBase::getTransientRecHits()
  TransientTrackingRecHit::ConstRecHitContainer trackerTransientHits(const reco::Track& trk) const;

  // helper: crude two-point-plus-origin circle estimate of |q/pt| and a
  // seed state for the standalone-only fit; returns false if the two
  // positions are (numerically) collinear with the origin.
  bool estimateSeedPt(const GlobalPoint& pInner, const GlobalPoint& pOuter, double& ptGuess) const;

  // helper: attempt the standalone-only fit for one segment cluster,
  // trying both charge hypotheses; returns the winning Trajectory or
  // nullptr. `fitter` is set up (cloned + hit-cloner attached) once per
  // event by the caller and reused across all clusters.
  std::unique_ptr<Trajectory> fitStandaloneCluster(
      const TransientTrackingRecHit::ConstRecHitContainer& orderedHits, TrajectoryFitter& fitter) const;

  // ---- inputs ----
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegsToken_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegsToken_;
  edm::EDGetTokenT<GEMSegmentCollection> gemSegsToken_;
  edm::EDGetTokenT<std::vector<MuonSegmentLink>> tsLinksToken_;
  edm::EDGetTokenT<std::vector<MuonSegmentLink>> ssLinksToken_;

  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken_;

  // rechit builders (same names GlobalTrajectoryBuilderBase uses)
  const std::string TrackerRecHitConfig_,MuonRecHitConfig_;
  
  edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> trackerRecHitBuilderToken_;
  edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> muonRecHitBuilderToken_;

  const TransientTrackingRecHitBuilder* trackerRecHitBuilder_ = nullptr;
  const TransientTrackingRecHitBuilder* muonRecHitBuilder_ = nullptr;

  //edm::ESHandle<TransientTrackingRecHitBuilder> trackerRecHitBuilder_;
  //edm::ESHandle<TransientTrackingRecHitBuilder> muonRecHitBuilder_;

  // generic Kalman fitter/smoother used ONLY for the standalone-only fit
  // (KFFitterSmootherSTA in the standard config -- see fillDescriptions)
  edm::ESGetToken<TrajectoryFitter, TrajectoryFitter::Record> staFitterToken_;

  // ---- thresholds / behaviour ----
  double minTSScore_;
  double minSSScore_;
  bool expandBySS_;
  unsigned int minSegmentsPerMuon_;  // TrackerMuon chamber-match threshold
  unsigned int minStationsForSTA_;   // min distinct stations to attempt an orphan STA fit
  double glbPtCut_;
  double glbPCut_;
  int glbMuonHitsOption_;  // passed straight to GlobalMuonRefitter::refit
  std::string trackerPropagatorName_;
  double staSeedErrorRescale_;

  // ---- reused CMSSW reconstruction machinery ----
  std::unique_ptr<MuonServiceProxy> service_;
  std::unique_ptr<GlobalMuonRefitter> glbRefitter_;
  std::unique_ptr<MuonTrackLoader> trackLoader_;

  bool debug_ = true;
};

#endif
