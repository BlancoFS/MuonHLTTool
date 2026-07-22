#ifndef MuonMLSegmentMatcher_MuonMLSegmentMatcherProducer_h
#define MuonMLSegmentMatcher_MuonMLSegmentMatcherProducer_h

// Standard
#include <string>
#include <vector>

// LibTorch
#include <torch/script.h>

// CMSSW Framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Muon reco objects
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"

// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

// -------------------------------------------------------------------------
// Output data product
// -------------------------------------------------------------------------

#include "RecoMuon/L3MuonProducer/interface/MuonSegmentLink.h"

// -------------------------------------------------------------------------
// EDProducer declaration
// -------------------------------------------------------------------------

class MuonMLSegmentMatcherProducer : public edm::stream::EDProducer<> {
public:
  explicit MuonMLSegmentMatcherProducer(const edm::ParameterSet&);
  ~MuonMLSegmentMatcherProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  // --- main produce method ---
  void produce(edm::Event&, const edm::EventSetup&) override;

  // --- feature builders (mirror training code exactly) ---
  std::vector<float> buildTrackFeatures(const reco::Track& trk) const;

  std::vector<float> buildDTSegmentFeatures(
      const DTRecSegment4D& seg,
      const DTGeometry& dtGeom) const;

  std::vector<float> buildCSCSegmentFeatures(
      const CSCSegment& seg,
      const CSCGeometry& cscGeom) const;

  std::vector<float> buildGEMSegmentFeatures(
      const GEMSegment& seg,
      const GEMGeometry& gemGeom) const;

  // --- edge builders for message-passing graph ---
  // At inference we have no sim info, so we connect nodes that are
  // geometrically compatible (dR or same hemisphere).
  std::vector<std::pair<int,int>> buildMPEdgesTS(
      const std::vector<std::vector<float>>& trackFeats,
      const std::vector<std::vector<float>>& segFeats) const;

  std::vector<std::pair<int,int>> buildMPEdgesSS(
      const std::vector<std::vector<float>>& segFeats) const;

  // --- tokens ---
  edm::EDGetTokenT<reco::TrackCollection>     tracksToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection>  dtSegsToken_;
  edm::EDGetTokenT<CSCSegmentCollection>      cscSegsToken_;
  edm::EDGetTokenT<GEMSegmentCollection>      gemSegsToken_;

  // --- geometry tokens ---
  edm::ESGetToken<DTGeometry,  MuonGeometryRecord> dtGeomToken_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;

  // --- model ---
  torch::jit::script::Module model_;

  // --- config ---
  edm::FileInPath modelPath_;
  const float       scoreCutTS_;       // minimum score to store a track-segment link
  const float       scoreCutSS_;       // minimum score to store a segment-segment link
  const float       mpDRCutTS_;        // dR cone for message-passing track–segment edges
  const float       mpDRCutSS_;        // dR cone for message-passing segment–segment edges

  // Number of input features (must match training)
  static constexpr int N_TRACK_FEAT = 13;
  static constexpr int N_SEG_FEAT   =  8;
};

#endif
