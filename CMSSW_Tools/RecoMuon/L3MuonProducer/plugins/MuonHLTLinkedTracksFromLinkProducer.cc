/*
 * MuonHLTL3FromLinkProducer.cc
 * ================================
 * Produces ONLY what MuonIdProducer needs to build muons itself:
 *   - "inner"     reco::TrackCollection : GNN-filtered subset of `tracks`
 *                 (only tracks the graph linked to >= minSegmentsPerMuon_
 *                 segments -- MuonIdProducer never even sees the rest)
 *   - (unnamed)   reco::MuonTrackLinksCollection : tracker+standalone+global
 *                 triplets, for candidates with both a successful combined
 *                 fit AND a companion standalone-only fit
 *   - "standAlone" reco::TrackCollection : every segment cluster (claimed by
 *                 a track or not) that fits standalone on its own
 *   - "global"    reco::TrackCollection : backing store for the global leg
 *                 of each MuonTrackLinks entry (not itself one of
 *                 MuonIdProducer's three inputCollectionLabels, but required
 *                 for the links' globalTrack() refs to be valid)
 *
 * No reco::Muon objects are built here -- MuonIdProducer does that itself,
 * including its own chamber-to-track propagation and matching, from these
 * raw collections. Building MuonChamberMatch/MuonSegmentMatch objects here
 * would just be redundant work MuonIdProducer throws away.
 *
 *   REUSED, verified against cms-sw/cmssw:
 *     - MuonServiceProxy            (RecoMuon/TrackingTools)
 *     - GlobalMuonRefitter          (RecoMuon/GlobalTrackingTools)  -- tracker+muon combined KF fit
 *     - MuonTrackLoader             (RecoMuon/TrackingTools)        -- Trajectory -> persisted reco::Track
 *     - trackerTransientHits()      -- mirrors GlobalTrajectoryBuilderBase::getTransientRecHits()
 *
 *   NEW, written here (tune/replace as needed):
 *     - estimateSeedPt()            -- 3-point (origin + 2 hits) circle estimate of pt
 *     - fitStandaloneCluster()      -- seeds a bare Kalman fit through segments with no tracker track,
 *                                      trying both charge hypotheses
 *
 * Efficiency notes:
 *   - segments, union-find clusters, and track<->segment grouping are each built ONCE in a
 *     single pass, and reused across the standalone and global fits alike.
 *   - each raw segment is wrapped into a MuonTransientTrackingRecHit exactly once (SegRecord::hit)
 *     and that same hit pointer is reused for both fits.
 *   - MuonTrackLoader::loadTracks() is called ONCE for all standalone candidates and ONCE for all
 *     global candidates (batched), not once per candidate.
 */

#include "RecoMuon/L3MuonProducer/plugins/MuonHLTLinkedTracksFromLinkProducer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <unordered_map>

#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"

#include "Geometry/CommonTopologies/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/GeomDetType.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

// -------------------------------------------------------------------------
// Internal per-segment record. Built once; reused by both fits. No local-
// frame chamber-match fields here anymore -- MuonIdProducer builds its own
// MuonChamberMatch/MuonSegmentMatch objects via its own track propagation,
// so keeping a parallel copy of that here would just be redundant work.
// -------------------------------------------------------------------------
struct MuonHLTLinkedTracksFromLinkProducer::SegRecord {
  enum Subdet { DT = 0, CSC = 1, GEM = 2 };

  Subdet subdet;
  int station;

  // built once via MuonRecHitBuilder; reused by both the global fit and the
  // standalone-only fit
  TransientTrackingRecHit::ConstRecHitPointer hit;
};

// -------------------------------------------------------------------------
// Constructor
// -------------------------------------------------------------------------
MuonHLTLinkedTracksFromLinkProducer::MuonHLTLinkedTracksFromLinkProducer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      dtSegsToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("dtSegments"))),
      cscSegsToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegments"))),
      gemSegsToken_(consumes<GEMSegmentCollection>(iConfig.getParameter<edm::InputTag>("gemSegments"))),
      tsLinksToken_(
          consumes<std::vector<MuonSegmentLink>>(iConfig.getParameter<edm::InputTag>("trackSegmentLinks"))),
      ssLinksToken_(
          consumes<std::vector<MuonSegmentLink>>(iConfig.getParameter<edm::InputTag>("segmentSegmentLinks"))),
      topoToken_(esConsumes()),
      fieldToken_(esConsumes()),
      TrackerRecHitConfig_(iConfig.getParameter<std::string>("TrackerRecHitBuilder")),
      MuonRecHitConfig_(iConfig.getParameter<std::string>("MuonRecHitBuilder")),
      staFitterToken_(esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("StandaloneFitter")))),
      minTSScore_(iConfig.getParameter<double>("minTSScore")),
      minSSScore_(iConfig.getParameter<double>("minSSScore")),
      expandBySS_(iConfig.getParameter<bool>("expandBySS")),
      minSegmentsPerMuon_(iConfig.getParameter<unsigned int>("minSegmentsPerMuon")),
      minStationsForSTA_(iConfig.getParameter<unsigned int>("minStationsForSTA")),
      glbPtCut_(iConfig.getParameter<double>("glbPtCut")),
      glbPCut_(iConfig.getParameter<double>("glbPCut")),
      glbMuonHitsOption_(iConfig.getParameter<int>("glbMuonHitsOption")),
      trackerPropagatorName_(iConfig.getParameter<std::string>("TrackerPropagator")),
      staSeedErrorRescale_(iConfig.getParameter<double>("staSeedErrorRescale")) {
  edm::ConsumesCollector iC = consumesCollector();

  service_ = std::make_unique<MuonServiceProxy>(iConfig.getParameter<edm::ParameterSet>("ServiceParameters"),
                                                consumesCollector());

  edm::ParameterSet glbRefitterPSet = iConfig.getParameter<edm::ParameterSet>("GlbRefitterParameters");
  glbRefitter_ = std::make_unique<GlobalMuonRefitter>(glbRefitterPSet, service_.get(), iC);

  edm::ParameterSet trackLoaderPSet = iConfig.getParameter<edm::ParameterSet>("TrackLoaderParameters");
  trackLoader_ = std::make_unique<MuonTrackLoader>(trackLoaderPSet, iC, service_.get());

  trackerRecHitBuilderToken_ = esConsumes<TransientTrackingRecHitBuilder, TransientRecHitRecord>(edm::ESInputTag("", TrackerRecHitConfig_));
  muonRecHitBuilderToken_ = esConsumes<TransientTrackingRecHitBuilder, TransientRecHitRecord>(edm::ESInputTag("", MuonRecHitConfig_));
  
  // GNN-filtered subset of `tracks` (segSet.size() >= minSegmentsPerMuon_)
  // -- "inner tracks" input for MuonIdProducer.
  produces<reco::TrackCollection>("inner").setBranchAlias("filteredInnerTracks");

  produces<reco::MuonTrackLinksCollection>();  // "links" input for MuonIdProducer

  produces<reco::TrackCollection>("global").setBranchAlias("globalMuonTracks");
  produces<reco::TrackExtraCollection>("global");
  produces<TrackingRecHitCollection>("global");
  produces<std::vector<Trajectory>>("global");
  produces<TrajTrackAssociationCollection>("global");

  produces<reco::TrackCollection>("standAlone").setBranchAlias("standAloneMuonTracks");  // "outer tracks"
  produces<reco::TrackExtraCollection>("standAlone");
  produces<TrackingRecHitCollection>("standAlone");
  produces<std::vector<Trajectory>>("standAlone");
  produces<TrajTrackAssociationCollection>("standAlone");
}

MuonHLTLinkedTracksFromLinkProducer::~MuonHLTLinkedTracksFromLinkProducer() {}

// -------------------------------------------------------------------------
// trackerTransientHits: mirrors GlobalTrajectoryBuilderBase::getTransientRecHits(),
// restricted to the Tracker-det branch (we add muon-system hits separately,
// from SegRecord::hit, since they come from the GNN's own selection rather
// than from track.recHits()).
// -------------------------------------------------------------------------
TransientTrackingRecHit::ConstRecHitContainer MuonHLTLinkedTracksFromLinkProducer::trackerTransientHits(
    const reco::Track& track) const {
  TransientTrackingRecHit::ConstRecHitContainer result;

  TrajectoryStateOnSurface currTsos =
      trajectoryStateTransform::innerStateOnSurface(track, *service_->trackingGeometry(), &*service_->magneticField());

  auto tkbuilder = static_cast<TkTransientTrackingRecHitBuilder const*>(trackerRecHitBuilder_);
  auto hitCloner = tkbuilder->cloner();

  for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
    if (!(*hit)->isValid())
      continue;
    if ((*hit)->geographicalId().det() != DetId::Tracker)
      continue;

    if (!(*hit)->hasPositionAndError()) {
      TrajectoryStateOnSurface predTsos =
          service_->propagator(trackerPropagatorName_)
              ->propagate(currTsos, service_->trackingGeometry()->idToDet((*hit)->geographicalId())->surface());
      if (!predTsos.isValid()) {
        edm::LogWarning("MuonHLTL3FromLinkProducer") << "tracker hit propagation failed; skipping hit";
        continue;
      }
      currTsos = predTsos;
      auto h = (**hit).cloneForFit(*tkbuilder->geometry()->idToDet((**hit).geographicalId()));
      result.emplace_back(hitCloner.makeShared(h, predTsos));
    } else {
      result.push_back((*hit)->cloneSH());
    }
  }
  return result;
}

// -------------------------------------------------------------------------
// estimateSeedPt: crude circle-through-the-origin estimate.
// Treats the beamline (0,0) as a third point on the trajectory in the
// transverse plane -- reasonable for a prompt muon, less so for a muon
// from a displaced vertex, but it only needs to get the fit's Kalman
// filter into roughly the right basin; the fit corrects it hit-by-hit.
// -------------------------------------------------------------------------
bool MuonHLTLinkedTracksFromLinkProducer::estimateSeedPt(const GlobalPoint& pInner,
                                                const GlobalPoint& pOuter,
                                                double& ptGuess) const {
  const double x1 = pInner.x(), y1 = pInner.y();
  const double x2 = pOuter.x(), y2 = pOuter.y();

  // 2 * area of the triangle (origin, p1, p2)
  const double area2 = std::abs(x1 * y2 - x2 * y1);
  if (area2 < 1e-3)
    return false;  // (nearly) collinear with the origin: no curvature handle

  const double a = std::hypot(x1, y1);
  const double b = std::hypot(x2 - x1, y2 - y1);
  const double c = std::hypot(x2, y2);
  const double R_cm = (a * b * c) / (2.0 * area2);  // circumradius, cm

  const GlobalPoint mid((x1 + x2) / 2.0, (y1 + y2) / 2.0, (pInner.z() + pOuter.z()) / 2.0);
  const double Bz = service_->magneticField()->inTesla(mid).z();

  ptGuess = 0.3 * std::abs(Bz) * (R_cm / 100.0);  // GeV
  ptGuess = std::clamp(ptGuess, 2.0, 200.0);
  return true;
}

// -------------------------------------------------------------------------
// fitStandaloneCluster: bare Kalman fit through GNN-clustered segments with
// no tracker seed. Tries both charge hypotheses, keeps the lower chi2/ndof.
// -------------------------------------------------------------------------
std::unique_ptr<Trajectory> MuonHLTLinkedTracksFromLinkProducer::fitStandaloneCluster(
    const TransientTrackingRecHit::ConstRecHitContainer& orderedHits, TrajectoryFitter& fitter) const {
  if (orderedHits.size() < 2)
    return nullptr;

  const GlobalPoint pIn = orderedHits.front()->globalPosition();
  const GlobalPoint pOut = orderedHits.back()->globalPosition();

  double ptGuess = 5.0;  // fallback if the circle estimate is degenerate
  estimateSeedPt(pIn, pOut, ptGuess);

  const GlobalVector dir = (pOut - pIn).unit();
  const double perp = std::max(0.05, static_cast<double>(dir.perp()));  // guard very forward tracks
  const GlobalVector mom = dir * (ptGuess / perp);

  AlgebraicSymMatrix55 seedErr;
  for (int i = 0; i < 5; ++i)
    seedErr(i, i) = staSeedErrorRescale_ * staSeedErrorRescale_;

  std::unique_ptr<Trajectory> best;
  double bestChi2Ndof = std::numeric_limits<double>::max();

  for (int charge : {1, -1}) {
    FreeTrajectoryState fts(GlobalTrajectoryParameters(pIn, mom, charge, &*service_->magneticField()),
                            CurvilinearTrajectoryError(seedErr));

    TrajectoryStateOnSurface firstTSOS(fts, orderedHits.front()->det()->surface());
    if (!firstTSOS.isValid())
      continue;

    PTrajectoryStateOnDet garbage1;
    edm::OwnVector<TrackingRecHit> garbage2;
    TrajectorySeed seed(garbage1, garbage2, alongMomentum);

    std::vector<Trajectory> result = fitter.fit(seed, orderedHits, firstTSOS);
    if (result.empty() || !result.front().isValid())
      continue;

    const double ndof = result.front().ndof();
    const double chi2Ndof = (ndof > 0) ? result.front().chiSquared() / ndof : std::numeric_limits<double>::max();
    if (chi2Ndof < bestChi2Ndof) {
      bestChi2Ndof = chi2Ndof;
      best = std::make_unique<Trajectory>(result.front());
    }
  }
  return best;
}

// -------------------------------------------------------------------------
// produce
// -------------------------------------------------------------------------
void MuonHLTLinkedTracksFromLinkProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  service_->update(iSetup);

  const TrackerTopology& tTopo = iSetup.getData(topoToken_);

  trackerRecHitBuilder_ = &iSetup.getData(trackerRecHitBuilderToken_);
  muonRecHitBuilder_ = &iSetup.getData(muonRecHitBuilderToken_);

  //trackerRecHitBuilder_ = iSetup.getHandle(trackerRecHitBuilderToken_);
  //muonRecHitBuilder_ = iSetup.getHandle(muonRecHitBuilderToken_);

  glbRefitter_->setEvent(iEvent);
  glbRefitter_->setServices(iSetup);

  // fitter for the standalone-only path: clone + attach a hit cloner once
  // per event, mirroring GlobalMuonRefitter::setServices().
  std::unique_ptr<TrajectoryFitter> staFitter = iSetup.getData(staFitterToken_).clone();
  {
    auto tkbuilder = static_cast<TkTransientTrackingRecHitBuilder const*>(trackerRecHitBuilder_);
    TkClonerImpl cloner = tkbuilder->cloner();
    staFitter->setHitCloner(&cloner);
  }

  auto tracksHandle = iEvent.getHandle(tracksToken_);
  const auto& tracks = *tracksHandle;

  const auto& dtSegs = iEvent.get(dtSegsToken_);
  const auto& cscSegs = iEvent.get(cscSegsToken_);
  const auto& gemSegs = iEvent.get(gemSegsToken_);
  const auto& tsLinks = iEvent.get(tsLinksToken_);
  const auto& ssLinks = iEvent.get(ssLinksToken_);

  const int nTracks = static_cast<int>(tracks.size());
  const int nDT = static_cast<int>(dtSegs.size());
  const int nCSC = static_cast<int>(cscSegs.size());
  const int nGEM = static_cast<int>(gemSegs.size());
  const int nSegs = nDT + nCSC + nGEM;

  // ---- pre-collect segment records, same flat order as the matcher: DT, CSC, GEM ----
  std::vector<SegRecord> segRecs;
  segRecs.reserve(nSegs);

  for (const auto& seg : dtSegs) {
    SegRecord r;
    r.subdet = SegRecord::DT;
    r.station = seg.chamberId().station();
    r.hit = muonRecHitBuilder_->build(&seg);
    segRecs.push_back(std::move(r));
  }
  for (const auto& seg : cscSegs) {
    SegRecord r;
    r.subdet = SegRecord::CSC;
    r.station = seg.cscDetId().station();
    r.hit = muonRecHitBuilder_->build(&seg);
    segRecs.push_back(std::move(r));
  }
  for (const auto& seg : gemSegs) {
    SegRecord r;
    r.subdet = SegRecord::GEM;
    r.station = seg.gemDetId().station();
    r.hit = muonRecHitBuilder_->build(&seg);
    segRecs.push_back(std::move(r));
  }

  if (nTracks == 0 && nSegs == 0) {
    iEvent.put(std::make_unique<reco::TrackCollection>(), "inner");
    iEvent.put(std::make_unique<reco::MuonTrackLinksCollection>());
    MuonCandidate::TrajectoryContainer emptyTraj;
    std::vector<bool> emptyOk;
    trackLoader_->loadTracks(emptyTraj, iEvent, emptyOk, tTopo, "global", true);
    trackLoader_->loadTracks(emptyTraj, iEvent, emptyOk, tTopo, "standAlone", true);
    return;
  }

  // ---- union-find over SS links: segment clusters, computed unconditionally ----
  std::vector<int> parent(nSegs);
  std::iota(parent.begin(), parent.end(), 0);
  auto find = [&](int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  auto unite = [&](int a, int b) {
    int ra = find(a);
    int rb = find(b);
    if (ra != rb)
      parent[ra] = rb;
  };

  // The format is the same as for tsLinks so that's why it calls trackIdx
  // However, it's not the trackIdx but the mother segment
  for (const auto& link : ssLinks) {
    if (link.score < minSSScore_)
      continue;
    if (link.trackIdx < 0 || link.trackIdx >= nSegs)
      continue;
    if (link.segmentIdx < 0 || link.segmentIdx >= nSegs)
      continue;
    unite(link.trackIdx, link.segmentIdx);
  }

  std::unordered_map<int, std::vector<int>> byRoot;
  byRoot.reserve(nSegs);
  for (int s = 0; s < nSegs; ++s)
    byRoot[find(s)].push_back(s);

  // ---- group TS links by track, optionally expand via clusters ----
  std::vector<std::set<int>> trackToSegs(nTracks);
  for (const auto& link : tsLinks) {
    if (link.score < minTSScore_)
      continue;
    if (link.trackIdx < 0 || link.trackIdx >= nTracks)
      continue;
    if (link.segmentIdx < 0 || link.segmentIdx >= nSegs)
      continue;
    trackToSegs[link.trackIdx].insert(link.segmentIdx);
  }

  if (expandBySS_) {
    for (auto& segSet : trackToSegs) {
      if (segSet.empty())
        continue;
      std::set<int> roots;
      for (int s : segSet)
        roots.insert(find(s));
      for (int root : roots) {
        const auto it = byRoot.find(root);
        if (it == byRoot.end())
          continue;
        for (int s : it->second)
          segSet.insert(s);
      }
    }
  }

  // =========================================================================
  // PASS 1: run the standalone-only fit for every union-find cluster that
  //         spans enough stations. Every fit that succeeds is persisted --
  //         this collection IS MuonIdProducer's "outer tracks", and also
  //         supplies the STA leg for any "link" built in PASS 2/3 below.
  // =========================================================================
  MuonCandidate::TrajectoryContainer staTrajs;
  std::vector<int> staTrajRoot;  // parallel to staTrajs: source cluster root

  std::set<int> processedRoots;
  for (int s = 0; s < nSegs; ++s) {
    const int root = find(s);
    if (!processedRoots.insert(root).second)
      continue;

    const auto& members = byRoot[root];

    std::set<int> stations;
    for (int m : members)
      stations.insert(segRecs[m].station * 10 + static_cast<int>(segRecs[m].subdet));
    if (stations.size() < minStationsForSTA_)
      continue;

    std::vector<int> ordered(members.begin(), members.end());
    std::sort(ordered.begin(), ordered.end(), [&](int a, int b) {
      return segRecs[a].hit->globalPosition().mag() < segRecs[b].hit->globalPosition().mag();
    });

    TransientTrackingRecHit::ConstRecHitContainer hits;
    hits.reserve(ordered.size());
    for (int idx : ordered)
      if (segRecs[idx].hit)
        hits.push_back(segRecs[idx].hit);
    if (hits.size() < 2)
      continue;

    std::unique_ptr<Trajectory> traj = fitStandaloneCluster(hits, *staFitter);
    if (!traj)
      continue;

    staTrajs.push_back(std::move(traj));
    staTrajRoot.push_back(root);
  }

  // batch-persist ALL standalone tracks now ("outer tracks"), and build a
  // root -> TrackRef map so PASS 2 can find each candidate's companion STA leg.
  std::unordered_map<int, reco::TrackRef> rootToStaRef;
  {
    std::vector<bool> ok(staTrajs.size(), false);
    edm::OrphanHandle<reco::TrackCollection> staHandle =
        trackLoader_->loadTracks(staTrajs, iEvent, ok, tTopo, "standAlone", /*reallyDoSmoothing=*/true);

    unsigned int outIdx = 0;
    for (size_t i = 0; i < ok.size(); ++i) {
      if (!ok[i])
        continue;
      rootToStaRef[staTrajRoot[i]] = reco::TrackRef(staHandle, outIdx++);
    }
  }

  // =========================================================================
  // PASS 2: for every track the GNN linked to >= minSegmentsPerMuon_
  //         segments, attempt the combined tracker+muon-system fit via
  //         GlobalMuonRefitter.
  // =========================================================================
  MuonCandidate::TrajectoryContainer globalTrajs;
  std::vector<int> globalTrajTrackIdx;
  std::vector<reco::TrackRef> globalTrajStaRef;  // parallel to globalTrajs: companion STA ref, or null
  std::vector<int> filteredInnerTrackIdx;        // GNN-filtered "inner tracks"

  for (int iTrack = 0; iTrack < nTracks; ++iTrack) {
    const auto& segSet = trackToSegs[iTrack];
    if (segSet.size() < minSegmentsPerMuon_)
      continue;

    filteredInnerTrackIdx.push_back(iTrack);

    const reco::Track& trk = tracks[iTrack];
    const reco::TrackRef trkRef(tracksHandle, iTrack);

    // ---- pT threshold for global-muon refit ----
    if (trk.p() < glbPCut_ || trk.pt() < glbPtCut_)
      continue;

    TransientTrackingRecHit::ConstRecHitContainer allHits = trackerTransientHits(trk);
    allHits.reserve(allHits.size() + segSet.size());
    for (int s : segSet)
      if (segRecs[s].hit)
        allHits.push_back(segRecs[s].hit);

    reco::TransientTrack tTT(trkRef, &*service_->magneticField(), service_->trackingGeometry());
    std::vector<Trajectory> refit = glbRefitter_->refit(trk, tTT, allHits, glbMuonHitsOption_, &tTopo);
    if (refit.empty())
      continue;

    // which cluster (union-find root) does this track's segment set draw
    // most of its segments from? that cluster's persisted STA track (if
    // any, from PASS 1) becomes the "muon track" leg of this candidate's
    // MuonTrackLinks entry.
    std::unordered_map<int, int> rootCounts;
    for (int s : segSet)
      rootCounts[find(s)]++;
    int domRoot = -1, domCount = -1;
    for (const auto& kv : rootCounts) {
      if (kv.second > domCount) {
        domCount = kv.second;
        domRoot = kv.first;
      }
    }

    reco::TrackRef staRef;  // left null if no companion STA fit exists
    const auto itSta = rootToStaRef.find(domRoot);
    if (itSta != rootToStaRef.end())
      staRef = itSta->second;

    globalTrajs.push_back(std::make_unique<Trajectory>(refit.front()));
    globalTrajTrackIdx.push_back(iTrack);
    globalTrajStaRef.push_back(staRef);
  }

  // GNN-filtered "inner tracks": plain copies of the surviving reco::Track
  // objects. Their TrackExtraRef / RecHit refs still point into the ORIGINAL
  // `tracks` product's extra/rechit collections (refs are position-based,
  // not memory-based, so this stays valid), so no need to also clone
  // TrackExtra/RecHits ourselves -- MuonIdProducer's track->extra() calls
  // resolve exactly as they would reading `tracks` directly.
  {
    auto filteredTracks = std::make_unique<reco::TrackCollection>();
    filteredTracks->reserve(filteredInnerTrackIdx.size());
    for (int idx : filteredInnerTrackIdx)
      filteredTracks->push_back(tracks[idx]);
    iEvent.put(std::move(filteredTracks), "inner");
  }

  // =========================================================================
  // PASS 3: batch-persist the global trajectories ("global", the backing
  //         store for the links' globalTrack() refs), then build
  //         MuonTrackLinksCollection ("links") for every candidate that has
  //         BOTH a persisted global track AND a companion STA track from
  //         PASS 1.
  //
  //         Together with "outer tracks" (the "standAlone" instance from
  //         PASS 1) and "inner tracks" (the "inner" instance emitted just
  //         above), this is exactly what MuonIdProducer wants for
  //         inputCollectionTypes = ['inner tracks', 'links', 'outer tracks'].
  // =========================================================================
  auto linksCollection = std::make_unique<reco::MuonTrackLinksCollection>();

  {
    std::vector<bool> ok(globalTrajs.size(), false);
    edm::OrphanHandle<reco::TrackCollection> glbHandle =
        trackLoader_->loadTracks(globalTrajs, iEvent, ok, tTopo, "global", /*reallyDoSmoothing=*/true);

    unsigned int outIdx = 0;
    for (size_t i = 0; i < ok.size(); ++i) {
      if (!ok[i])
        continue;
      const reco::TrackRef gref(glbHandle, outIdx++);
      const reco::TrackRef& staRef = globalTrajStaRef[i];
      if (staRef.isNull())
        continue;  // global fit succeeded but no companion STA fit -- no link;
                   // MuonIdProducer will still see this track via "inner tracks"

      const reco::TrackRef trkRef(tracksHandle, globalTrajTrackIdx[i]);
      reco::MuonTrackLinks links;
      links.setTrackerTrack(trkRef);
      links.setStandAloneTrack(staRef);
      links.setGlobalTrack(gref);
      linksCollection->push_back(links);
    }
  }

  if (debug_)
    std::cout << "Event: nTracks=" << nTracks << " nSegs=" << nSegs << " -> innerTracks="
              << filteredInnerTrackIdx.size() << " (global attempts=" << globalTrajs.size()
              << ", standalone attempts=" << staTrajs.size() << ", links=" << linksCollection->size() << ")"
              << std::endl;

  iEvent.put(std::move(linksCollection));
}

// -------------------------------------------------------------------------
// fillDescriptions
// -------------------------------------------------------------------------
void MuonHLTLinkedTracksFromLinkProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("hltGeneralTracks"));
  desc.add<edm::InputTag>("dtSegments", edm::InputTag("hltDt4DSegments"));
  desc.add<edm::InputTag>("cscSegments", edm::InputTag("hltCscSegments"));
  desc.add<edm::InputTag>("gemSegments", edm::InputTag("hltGemSegments"));
  desc.add<edm::InputTag>("trackSegmentLinks", edm::InputTag("muonHLTSegmentMatcher", "trackSegmentLinks"));
  desc.add<edm::InputTag>("segmentSegmentLinks", edm::InputTag("muonHLTSegmentMatcher", "segmentSegmentLinks"));

  desc.add<double>("minTSScore", 0.5);
  desc.add<double>("minSSScore", 0.5);
  desc.add<bool>("expandBySS", true);
  desc.add<unsigned int>("minSegmentsPerMuon", 1);
  desc.add<unsigned int>("minStationsForSTA", 2);

  desc.add<double>("glbPtCut", 1.0);
  desc.add<double>("glbPCut", 2.5);
  desc.add<int>("glbMuonHitsOption", 1);  // see GlobalMuonRefitter::refit: 1 = use all supplied muon hits

  desc.add<std::string>("TrackerRecHitBuilder", "WithAngleAndTemplate");
  desc.add<std::string>("MuonRecHitBuilder", "MuonRecHitBuilder");
  desc.add<std::string>("TrackerPropagator", "SteppingHelixPropagatorAny");
  desc.add<std::string>("StandaloneFitter", "KFFitterSmootherSTA");
  desc.add<double>("staSeedErrorRescale", 100.0);

  // Pass-through PSets: copy these from the real cff files rather than
  // hand-declaring every nested key here (their schemas are owned by the
  // classes that consume them):
  //   ServiceParameters    <- RecoMuon.TrackingTools.MuonServiceProxy_cff.MuonServiceProxy
  //   GlbRefitterParameters <- RecoMuon.GlobalTrackingTools.GlobalMuonRefitter_cff.GlobalMuonRefitter
  //   TrackLoaderParameters <- RecoMuon.TrackingTools.MuonTrackLoader_cff.MuonTrackLoaderForGLB
  //     (set TrackLoaderParameters.DoSmoothing = True; VertexConstraint = False
  //      keeps the product menu declared in the constructor above sufficient)
  edm::ParameterSetDescription servicePSet;
  servicePSet.setAllowAnything();
  desc.add<edm::ParameterSetDescription>("ServiceParameters", servicePSet);

  edm::ParameterSetDescription glbRefitterPSet;
  glbRefitterPSet.setAllowAnything();
  desc.add<edm::ParameterSetDescription>("GlbRefitterParameters", glbRefitterPSet);

  edm::ParameterSetDescription trackLoaderPSet;
  trackLoaderPSet.setAllowAnything();
  desc.add<edm::ParameterSetDescription>("TrackLoaderParameters", trackLoaderPSet);

  descriptions.addWithDefaultLabel(desc);
}
#include <FWCore/Framework/interface/MakerMacros.h>

DEFINE_FWK_MODULE(MuonHLTLinkedTracksFromLinkProducer);
