/*
 * MuonHLTSegmentMatcherProducer.cc
 * =================================
 * CMSSW EDProducer that runs the exported TorchScript GNN to score
 * track-segment and segment-segment links.
 *
 * Workflow per event:
 *  1. Read tracks, DT/CSC/GEM segments from the event.
 *  2. Build node-feature tensors (same order & normalisation as training).
 *  3. Build message-passing edge_index tensors using a dR heuristic
 *     (no sim info available at HLT) -- this MUST match the geometric
 *     heuristic used to build message-passing edges in training
 *     (transformData.build_mp_edges_ts/ss), not sim truth.
 *  4. Build query edge_label_index for (track, segment) and
 *     (segment, segment) pairs, restricted to a generous dR window
 *     (not an exhaustive all-pairs scan). Geometric edge features
 *     (deta, dphi, dR) are computed inside the exported TorchScript
 *     model itself (export_model.py's InferenceWrapper), not here.
 *  5. Run the TorchScript model and creates sigmoid scores.
 *  6. Store links above threshold as the output data product.
 */

#include "RecoMuon/L3MuonProducer/plugins/MuonMLSegmentMatcherProducer.h"

#include <cmath>
#include <numeric>

// -------------------------------------------------------------------------
// Constructor
// -------------------------------------------------------------------------

MuonMLSegmentMatcherProducer::MuonMLSegmentMatcherProducer(
    const edm::ParameterSet& iConfig)
    : tracksToken_  (consumes<reco::TrackCollection>    (iConfig.getParameter<edm::InputTag>("tracks")))
    , dtSegsToken_  (consumes<DTRecSegment4DCollection> (iConfig.getParameter<edm::InputTag>("dtSegments")))
    , cscSegsToken_ (consumes<CSCSegmentCollection>     (iConfig.getParameter<edm::InputTag>("cscSegments")))
    , gemSegsToken_ (consumes<GEMSegmentCollection>     (iConfig.getParameter<edm::InputTag>("gemSegments")))
    , dtGeomToken_  (esConsumes<DTGeometry,  MuonGeometryRecord>())
    , cscGeomToken_ (esConsumes<CSCGeometry, MuonGeometryRecord>())
    , gemGeomToken_ (esConsumes<GEMGeometry, MuonGeometryRecord>())
    , modelPath_    (iConfig.getParameter<edm::FileInPath>("modelPath")) //, modelPath_    (iConfig.getParameter<std::string>("modelPath"))
    , scoreCutTS_   (iConfig.getParameter<double>("scoreCutTS"))
    , scoreCutSS_   (iConfig.getParameter<double>("scoreCutSS"))
    , mpDRCutTS_    (iConfig.getParameter<double>("mpDRCutTS"))
    , mpDRCutSS_    (iConfig.getParameter<double>("mpDRCutSS"))
{
  // Load the TorchScript model once at construction time
  try {
    model_ = torch::jit::load(modelPath_.fullPath());
    model_.eval();
    edm::LogInfo("MuonHLTSegmentMatcher") << "Model loaded on CPU: " << modelPath_;
  } catch (const c10::Error& e) {
    throw cms::Exception("MuonHLTSegmentMatcher")
      << "Failed to load TorchScript model from " << modelPath_
      << "\n" << e.what();
  }
  // Move to GPU if available
  // if (torch::cuda::is_available()) {
  //    model_.to(torch::kCUDA);
  //    edm::LogInfo("MuonHLTSegmentMatcher") << "Model loaded on GPU: " << modelPath_;
  //  } else {
  //    edm::LogInfo("MuonHLTSegmentMatcher") << "Model loaded on CPU: " << modelPath_;
  //  }
  //} catch (const c10::Error& e) {
  //  throw cms::Exception("MuonHLTSegmentMatcher")
  //      << "Failed to load TorchScript model from " << modelPath_
  //      << "\n" << e.what();
  //}

  produces<std::vector<MuonSegmentLink>>("trackSegmentLinks");
  produces<std::vector<MuonSegmentLink>>("segmentSegmentLinks");
}

// -------------------------------------------------------------------------
// Geometric helpers (moved up from below fillDescriptions so `produce()`
// can also use them to filter query candidate pairs, not just the
// message-passing edges)
// -------------------------------------------------------------------------

namespace {
  inline float deltaPhi(float a, float b) {
    float d = a - b;
    while (d >  M_PI) d -= 2 * M_PI;
    while (d < -M_PI) d += 2 * M_PI;
    return d;
  }
  inline float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = deltaPhi(phi1, phi2);
    return std::sqrt(deta * deta + dphi * dphi);
  }
}

// -------------------------------------------------------------------------
// produce
// -------------------------------------------------------------------------

void MuonMLSegmentMatcherProducer::produce(
    edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // ---- geometry ----
  const auto& dtGeom  = iSetup.getData(dtGeomToken_);
  const auto& cscGeom = iSetup.getData(cscGeomToken_);
  const auto& gemGeom = iSetup.getData(gemGeomToken_);

  // ---- collections ----
  const auto& tracks  = iEvent.get(tracksToken_);
  const auto& dtSegs  = iEvent.get(dtSegsToken_);
  const auto& cscSegs = iEvent.get(cscSegsToken_);
  const auto& gemSegs = iEvent.get(gemSegsToken_);

  const int nTracks = static_cast<int>(tracks.size());
  const int nDT     = static_cast<int>(dtSegs.size());
  const int nCSC    = static_cast<int>(cscSegs.size());
  const int nGEM    = static_cast<int>(gemSegs.size());
  const int nSegs   = nDT + nCSC + nGEM;

  auto tsLinks = std::make_unique<std::vector<MuonSegmentLink>>();
  auto ssLinks = std::make_unique<std::vector<MuonSegmentLink>>();

  // Skip empty events
  if (nTracks == 0 || nSegs == 0) {
    iEvent.put(std::move(tsLinks), "trackSegmentLinks");
    iEvent.put(std::move(ssLinks), "segmentSegmentLinks");
    return;
  }

  // ---- build node features ----
  // Track features: (nTracks, 13)
  std::vector<std::vector<float>> trackFeats;
  trackFeats.reserve(nTracks);
  for (const auto& trk : tracks)
    trackFeats.push_back(buildTrackFeatures(trk));

  // Segment features: (nSegs, 8) in order DT | CSC | GEM
  std::vector<std::vector<float>> segFeats;
  segFeats.reserve(nSegs);
  for (const auto& seg : dtSegs)
    segFeats.push_back(buildDTSegmentFeatures(seg, dtGeom));
  for (const auto& seg : cscSegs)
    segFeats.push_back(buildCSCSegmentFeatures(seg, cscGeom));
  for (const auto& seg : gemSegs)
    segFeats.push_back(buildGEMSegmentFeatures(seg, gemGeom));

  // ---- flatten into tensors ----
  const torch::Device device = torch::kCPU;
  //torch::Device device = torch::cuda::is_available()
  //                       ? torch::Device(torch::kCUDA)
  //                       : torch::Device(torch::kCPU);

  // track_x : (nTracks, N_TRACK_FEAT)
  auto track_x = torch::zeros({nTracks, N_TRACK_FEAT});
  for (int i = 0; i < nTracks; ++i)
    for (int f = 0; f < N_TRACK_FEAT; ++f)
      track_x[i][f] = trackFeats[i][f];

  // segment_x : (nSegs, N_SEG_FEAT)
  auto segment_x = torch::zeros({nSegs, N_SEG_FEAT});
  for (int i = 0; i < nSegs; ++i)
    for (int f = 0; f < N_SEG_FEAT; ++f)
      segment_x[i][f] = segFeats[i][f];

  // ---- message-passing edges (geometric heuristic) ----
  auto mpEdgesTS = buildMPEdgesTS(trackFeats, segFeats);
  auto mpEdgesSS = buildMPEdgesSS(segFeats);

  // ei_matched / ei_rev_matched : (2, E_mp_ts)
  const int E_mp_ts = static_cast<int>(mpEdgesTS.size());
  auto ei_matched     = torch::zeros({2, E_mp_ts}, torch::kLong);
  auto ei_rev_matched = torch::zeros({2, E_mp_ts}, torch::kLong);
  for (int e = 0; e < E_mp_ts; ++e) {
    ei_matched[0][e]     = mpEdgesTS[e].first;   // track index
    ei_matched[1][e]     = mpEdgesTS[e].second;  // segment index
    ei_rev_matched[0][e] = mpEdgesTS[e].second;
    ei_rev_matched[1][e] = mpEdgesTS[e].first;
  }

  // ei_co_seg : (2, E_mp_ss)
  const int E_mp_ss = static_cast<int>(mpEdgesSS.size());
  auto ei_co_seg = torch::zeros({2, E_mp_ss}, torch::kLong);
  for (int e = 0; e < E_mp_ss; ++e) {
    ei_co_seg[0][e] = mpEdgesSS[e].first;
    ei_co_seg[1][e] = mpEdgesSS[e].second;
  }

  // ---- query edges: geometric-window candidates, NOT exhaustive all-pairs ----
  // Scoring literally every track x segment / segment x segment combination
  // means the vast majority of candidates the model ever sees are wildly
  // uncorrelated pairs on opposite sides of the detector -- combinatorics
  // the training procedure never meaningfully taught the network to reject,
  // and which dominate both the "too many spurious links" symptom and the
  // per-event runtime. We restrict candidates to a window well beyond the
  // message-passing dR cut (mpDRCutTS_ / mpDRCutSS_) so no genuine
  // wide-angle match is silently dropped -- this is a generous combinatorial
  // pre-filter, not a physics cut; the model score + scoreCutTS_/scoreCutSS_
  // below still makes the final accept/reject decision.
  //
  // NOTE: we do NOT compute (deta, dphi, dR) edge features here. The
  // exported TorchScript model (see export_model.py's InferenceWrapper)
  // derives them internally from track_x/segment_x + eli_ts/eli_ss, using
  // the exact same code path as training. Duplicating that arithmetic here
  // in C++ would reintroduce a second, independent implementation of the
  // same computation -- precisely the kind of train/inference split that
  // caused the original truth-edges-vs-dR-edges bug.
  const double queryDRCutTS = 3.0 * mpDRCutTS_;
  const double queryDRCutSS = 3.0 * mpDRCutSS_;

  std::vector<int> tsTrackIdx, tsSegIdx;
  for (int t = 0; t < nTracks; ++t) {
    float t_eta = trackFeats[t][1];
    float t_phi = trackFeats[t][2];
    for (int s = 0; s < nSegs; ++s) {
      float s_eta = segFeats[s][3];
      float s_phi = segFeats[s][4];
      if (deltaR(t_eta, t_phi, s_eta, s_phi) < queryDRCutTS) {
        tsTrackIdx.push_back(t);
        tsSegIdx.push_back(s);
      }
    }
  }
  const int E_ts = static_cast<int>(tsTrackIdx.size());
  auto eli_ts = torch::zeros({2, E_ts}, torch::kLong);
  for (int e = 0; e < E_ts; ++e) {
    eli_ts[0][e] = tsTrackIdx[e];
    eli_ts[1][e] = tsSegIdx[e];
  }

  std::vector<int> ssI, ssJ;
  for (int i = 0; i < nSegs; ++i) {
    float eta_i = segFeats[i][3];
    float phi_i = segFeats[i][4];
    for (int j = 0; j < nSegs; ++j) {
      if (i == j) continue;
      float eta_j = segFeats[j][3];
      float phi_j = segFeats[j][4];
      if (deltaR(eta_i, phi_i, eta_j, phi_j) < queryDRCutSS) {
        ssI.push_back(i);
        ssJ.push_back(j);
      }
    }
  }
  const int E_ss = static_cast<int>(ssI.size());
  auto eli_ss = torch::zeros({2, E_ss}, torch::kLong);
  for (int e = 0; e < E_ss; ++e) {
    eli_ss[0][e] = ssI[e];
    eli_ss[1][e] = ssJ[e];
  }

  // ---- move all tensors to device ----
  track_x         = track_x.to(device);
  segment_x       = segment_x.to(device);
  ei_matched      = ei_matched.to(device);
  ei_rev_matched  = ei_rev_matched.to(device);
  ei_co_seg       = ei_co_seg.to(device);
  eli_ts          = eli_ts.to(device);
  eli_ss          = eli_ss.to(device);

  // ---- run inference ----
  // Same 7-tensor signature InferenceWrapper.forward() expects in
  // export_model.py -- track/segment features plus message-passing and
  // query edge indices. Edge attributes are computed inside the traced
  // model, not passed in from here.
  std::vector<torch::jit::IValue> inputs = {
      track_x, segment_x,
      ei_matched, ei_rev_matched, ei_co_seg,
      eli_ts, eli_ss
  };

  torch::jit::IValue output;
  try {
    torch::NoGradGuard no_grad;
    output = model_.forward(inputs);
  } catch (const c10::Error& e) {
    throw cms::Exception("MuonHLTSegmentMatcher")
        << "Model inference failed: " << e.what();
  }

  // ---- unpack outputs ----
  auto out_tuple = output.toTuple();
  auto scores_ts = torch::sigmoid(out_tuple->elements()[0].toTensor()).to(torch::kCPU);
  auto scores_ss = torch::sigmoid(out_tuple->elements()[1].toTensor()).to(torch::kCPU);

  // ---- store links above threshold ----
  auto scores_ts_acc = scores_ts.accessor<float, 1>();
  for (int e = 0; e < E_ts; ++e) {
    float score = scores_ts_acc[e];
    if (score >= scoreCutTS_) {
      MuonSegmentLink link;
      link.trackIdx   = static_cast<int>(eli_ts[0][e].item<int64_t>());
      link.segmentIdx = static_cast<int>(eli_ts[1][e].item<int64_t>());
      link.score      = score;
      tsLinks->push_back(link);
    }
  }

  auto scores_ss_acc = scores_ss.accessor<float, 1>();
  for (int e = 0; e < E_ss; ++e) {
    float score = scores_ss_acc[e];
    if (score >= scoreCutSS_) {
      MuonSegmentLink link;
      link.trackIdx   = static_cast<int>(eli_ss[0][e].item<int64_t>());
      link.segmentIdx = static_cast<int>(eli_ss[1][e].item<int64_t>());
      link.score      = score;
      ssLinks->push_back(link);
    }
  }

  bool debug = true;
  //edm::LogDebug("MuonMLSegmentMatcherProducer")
  if (debug)
    std::cout << "Event: nTracks=" << nTracks << " nSegs=" << nSegs
	      << " -> tsLinks=" << tsLinks->size()
	      << " ssLinks=" << ssLinks->size() << std::endl;

  iEvent.put(std::move(tsLinks), "trackSegmentLinks");
  iEvent.put(std::move(ssLinks), "segmentSegmentLinks");
}

// -------------------------------------------------------------------------
// Feature builders  (must mirror train_sage_muon.py exactly)
// -------------------------------------------------------------------------

std::vector<float>
MuonMLSegmentMatcherProducer::buildTrackFeatures(const reco::Track& trk) const
{
  float chi2_ndof = (trk.ndof() > 0) ? trk.chi2() / trk.ndof() : 0.f;

  return {
    static_cast<float>(trk.pt()),
    static_cast<float>(trk.eta()),
    static_cast<float>(trk.phi()),
    static_cast<float>(trk.charge()),
    chi2_ndof,
    static_cast<float>(trk.numberOfValidHits()),
    static_cast<float>(trk.numberOfLostHits()),
    static_cast<float>(trk.outerPosition().x()),
    static_cast<float>(trk.outerPosition().y()),
    static_cast<float>(trk.outerPosition().z()),
    static_cast<float>(trk.outerMomentum().x()),
    static_cast<float>(trk.outerMomentum().y()),
    static_cast<float>(trk.outerMomentum().z()),
  };
}

std::vector<float>
MuonMLSegmentMatcherProducer::buildDTSegmentFeatures(
    const DTRecSegment4D& seg, const DTGeometry& dtGeom) const
{
  const auto* chamber = dtGeom.chamber(seg.chamberId());
  const auto  gp      = chamber->toGlobal(seg.localPosition());
  const float eta     = gp.eta();
  const float phi     = gp.phi();
  const float chi2    = seg.hasPhi() ? seg.phiSegment()->chi2() : 0.f;

  return {
    static_cast<float>(gp.x()),
    static_cast<float>(gp.y()),
    static_cast<float>(gp.z()),
    eta,
    phi,
    chi2,
    static_cast<float>(seg.chamberId().station()),
    1,
  };
}

std::vector<float>
MuonMLSegmentMatcherProducer::buildCSCSegmentFeatures(
    const CSCSegment& seg, const CSCGeometry& cscGeom) const
{
  const auto* chamber = cscGeom.chamber(seg.cscDetId());
  const auto  gp      = chamber->toGlobal(seg.localPosition());

  return {
    static_cast<float>(gp.x()),
    static_cast<float>(gp.y()),
    static_cast<float>(gp.z()),
    static_cast<float>(gp.eta()),
    static_cast<float>(gp.phi()),
    static_cast<float>(seg.chi2()),
    static_cast<float>(seg.cscDetId().station()),
    2,
  };
}

std::vector<float>
MuonMLSegmentMatcherProducer::buildGEMSegmentFeatures(
    const GEMSegment& seg, const GEMGeometry& gemGeom) const
{
  const auto* superChamber = gemGeom.superChamber(seg.gemDetId());
  const auto  gp           = superChamber->toGlobal(seg.localPosition());

  return {
    static_cast<float>(gp.x()),
    static_cast<float>(gp.y()),
    static_cast<float>(gp.z()),
    static_cast<float>(gp.eta()),
    static_cast<float>(gp.phi()),
    static_cast<float>(seg.chi2()),
    static_cast<float>(seg.gemDetId().station()),
    3,
  };
}

// -------------------------------------------------------------------------
// Message-passing edge builders (dR heuristic)
// Feature layout: [x, y, z, eta(3), phi(4), chi2, station]
//                  0  1  2   3        4       5     6
// -------------------------------------------------------------------------

std::vector<std::pair<int,int>>
MuonMLSegmentMatcherProducer::buildMPEdgesTS(
    const std::vector<std::vector<float>>& trackFeats,
    const std::vector<std::vector<float>>& segFeats) const
{
  // Track features: eta=1, phi=2  (from buildTrackFeatures)
  // Segment features: eta=3, phi=4 (from buildSegmentFeatures)
  std::vector<std::pair<int,int>> edges;
  for (int t = 0; t < static_cast<int>(trackFeats.size()); ++t) {
    float t_eta = trackFeats[t][1];
    float t_phi = trackFeats[t][2];
    for (int s = 0; s < static_cast<int>(segFeats.size()); ++s) {
      float s_eta = segFeats[s][3];
      float s_phi = segFeats[s][4];
      if (deltaR(t_eta, t_phi, s_eta, s_phi) < mpDRCutTS_)
        edges.emplace_back(t, s);
    }
  }
  return edges;
}

std::vector<std::pair<int,int>>
MuonMLSegmentMatcherProducer::buildMPEdgesSS(
    const std::vector<std::vector<float>>& segFeats) const
{
  std::vector<std::pair<int,int>> edges;
  const int n = static_cast<int>(segFeats.size());
  for (int i = 0; i < n; ++i) {
    float eta_i = segFeats[i][3];
    float phi_i = segFeats[i][4];
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      float eta_j = segFeats[j][3];
      float phi_j = segFeats[j][4];
      if (deltaR(eta_i, phi_i, eta_j, phi_j) < mpDRCutSS_)
        edges.emplace_back(i, j);
    }
  }
  return edges;
}

// -------------------------------------------------------------------------
// fillDescriptions
// -------------------------------------------------------------------------

void MuonMLSegmentMatcherProducer::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks",      edm::InputTag("hltGeneralTracks"));
  desc.add<edm::InputTag>("dtSegments",  edm::InputTag("hltDt4DSegments"));
  desc.add<edm::InputTag>("cscSegments", edm::InputTag("hltCscSegments"));
  desc.add<edm::InputTag>("gemSegments", edm::InputTag("hltGemSegments"));
  desc.add<edm::FileInPath>("modelPath", edm::FileInPath("RecoMuon/L3MuonProducer/data/MuonSegmentMatcher_Jul26_CMSSW16X_v1.pt"));
  desc.add<double>       ("scoreCutTS",  0.5);
  desc.add<double>       ("scoreCutSS",  0.5);
  desc.add<double>       ("mpDRCutTS",   0.8);
  desc.add<double>       ("mpDRCutSS",   0.8);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(MuonMLSegmentMatcherProducer);
