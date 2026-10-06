#ifndef GenericMuonCollections_h
#define GenericMuonCollections_h

// -----------------------------------------------------------------------
// Purpose
// -----------------------------------------------------------------------
// Today MuonHLTNtupler has one hardcoded edm::EDGetTokenT<T> + one hardcoded
// set of branch arrays PER named collection (L2Muon, L3Muon, TkMuon,
// iterL3OI, iterL3IOFromL2, ...). Every time the muon HLT sequence is
// restructured, someone has to edit the .h, the .cc, and the .cfi/customizer
// in three places just to follow a rename.
//
// This file replaces that with a small type-erasure layer: the *type* of
// EDM product (RecoChargedCandidate, Track view, reco::Muon, L1TkMuon, ...)
// is still known at compile time (CMSSW requires that), but *which* input
// tag maps to *which* branch name is entirely data-driven from a python
// cms.VPSet. Adding, removing, or renaming a collection becomes a one-line
// python change; no recompilation of ntupler logic is needed. Supporting a
// genuinely NEW EDM type (not just a new instance name) is the only case
// that still touches this file, and it's a single well-isolated
// specialization added at the bottom.
// -----------------------------------------------------------------------

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include <memory>
#include <string>
#include <vector>

// -- One flat, name-agnostic record per reconstructed "muon-like" object,
//    no matter which underlying EDM collection it came from.
//
//    reco::MuonTrackLinks-based collections (iterL3OI, iterL3IOFromL2,
//    iterL3FromL2 in the original code) carry THREE tracks per link
//    (tracker/inner, standalone/outer, global) rather than one -- the
//    inner_*/outer_*/global_* fields below cover that case and stay at
//    -99 for every other collection type.
struct GenericMuon {
  float pt = -99.f;
  float eta = -99.f;
  float phi = -99.f;
  float charge = 0.f;
  float trackPt = -99.f;  // inner/tracker-track pt when available, else == pt

  bool hasInner = false, hasOuter = false, hasGlobal = false;
  float inner_pt = -99.f, inner_eta = -99.f, inner_phi = -99.f, inner_charge = -99.f;
  float outer_pt = -99.f, outer_eta = -99.f, outer_phi = -99.f, outer_charge = -99.f;
  float global_pt = -99.f, global_eta = -99.f, global_phi = -99.f, global_charge = -99.f;
};

// -- One entry of the "muonCollections" VPSet coming from python.
struct MuonCollectionCfg {
  std::string name;  // branch-name prefix, e.g. "hltL3Muon"
  std::string type;  // selects which adapter specialization below to use
  edm::InputTag tag;
};

// -- Interface every per-type adapter implements. The ntupler only ever
//    talks to this interface, never to the concrete EDM type.
class IMuonCollectionAdapter {
public:
  virtual ~IMuonCollectionAdapter() = default;
  virtual void getByEvent(const edm::Event& iEvent) = 0;
  virtual bool isValid() const = 0;
  virtual std::vector<GenericMuon> extract() const = 0;
};

// -- One adapter template per underlying EDM type. `extract()` is
//    specialized per T below; the token plumbing is shared.
template <typename T>
class MuonCollectionAdapter : public IMuonCollectionAdapter {
public:
  MuonCollectionAdapter(const edm::InputTag& tag, edm::ConsumesCollector&& cc)
      : token_(cc.consumes<T>(tag)) {}

  void getByEvent(const edm::Event& iEvent) override { handle_ = iEvent.getHandle(token_); }
  bool isValid() const override { return handle_.isValid(); }
  std::vector<GenericMuon> extract() const override;  // specialized below

private:
  edm::EDGetTokenT<T> token_;
  edm::Handle<T> handle_;
};

// ---------------------------------------------------------------------
// Specializations for the EDM types actually seen in Phase-2 muon HLT.
// Add a new one here ONLY when a genuinely new product type shows up
// (e.g. a new L1 format) -- never for a rename of an existing collection.
// ---------------------------------------------------------------------

template <>
inline std::vector<GenericMuon> MuonCollectionAdapter<reco::RecoChargedCandidateCollection>::extract() const {
  std::vector<GenericMuon> out;
  out.reserve(handle_->size());
  for (const auto& c : *handle_) {
    GenericMuon m;
    m.pt = c.pt();
    m.eta = c.eta();
    m.phi = c.phi();
    m.charge = c.charge();
    m.trackPt = c.track().isNonnull() ? c.track()->pt() : c.pt();
    out.push_back(m);
  }
  return out;
}

template <>
inline std::vector<GenericMuon> MuonCollectionAdapter<edm::View<reco::Track>>::extract() const {
  std::vector<GenericMuon> out;
  out.reserve(handle_->size());
  for (const auto& t : *handle_) {
    GenericMuon m;
    m.pt = t.pt();
    m.eta = t.eta();
    m.phi = t.phi();
    m.charge = t.charge();
    m.trackPt = t.pt();
    out.push_back(m);
  }
  return out;
}

template <>
inline std::vector<GenericMuon> MuonCollectionAdapter<std::vector<reco::Muon>>::extract() const {
  std::vector<GenericMuon> out;
  out.reserve(handle_->size());
  for (const auto& mu : *handle_) {
    GenericMuon m;
    m.pt = mu.pt();
    m.eta = mu.eta();
    m.phi = mu.phi();
    m.charge = mu.charge();
    m.trackPt = mu.innerTrack().isNonnull() ? mu.innerTrack()->pt() : mu.pt();
    out.push_back(m);
  }
  return out;
}

// -- reco::MuonTrackLinks: one entry can carry an inner (tracker) track, an
//    outer (standalone) track, and a global track, any of which may be
//    null. This is the type behind iterL3OI / iterL3IOFromL2 / iterL3FromL2
//    in the original code, which is why those three had "_inner_",
//    "_outer_", "_global_" branch triplets instead of a single pt/eta/phi.
//    `pt`/`eta`/`phi`/`charge` below default to the global track when
//    present (falling back to inner, then outer) so this collection type
//    can still be treated as "just another GenericMuon" by generic code
//    that only cares about a single 4-vector.
template <>
inline std::vector<GenericMuon> MuonCollectionAdapter<std::vector<reco::MuonTrackLinks>>::extract() const {
  std::vector<GenericMuon> out;
  out.reserve(handle_->size());
  for (const auto& link : *handle_) {
    GenericMuon m;

    if (link.trackerTrack().isNonnull()) {
      const auto& t = *link.trackerTrack();
      m.hasInner = true;
      m.inner_pt = t.pt(); m.inner_eta = t.eta(); m.inner_phi = t.phi(); m.inner_charge = t.charge();
    }
    if (link.standAloneTrack().isNonnull()) {
      const auto& t = *link.standAloneTrack();
      m.hasOuter = true;
      m.outer_pt = t.pt(); m.outer_eta = t.eta(); m.outer_phi = t.phi(); m.outer_charge = t.charge();
    }
    if (link.globalTrack().isNonnull()) {
      const auto& t = *link.globalTrack();
      m.hasGlobal = true;
      m.global_pt = t.pt(); m.global_eta = t.eta(); m.global_phi = t.phi(); m.global_charge = t.charge();
    }

    // -- best-available single 4-vector, preferring global > inner > outer
    if (m.hasGlobal)      { m.pt = m.global_pt; m.eta = m.global_eta; m.phi = m.global_phi; m.charge = m.global_charge; }
    else if (m.hasInner)  { m.pt = m.inner_pt;  m.eta = m.inner_eta;  m.phi = m.inner_phi;  m.charge = m.inner_charge; }
    else if (m.hasOuter)  { m.pt = m.outer_pt;  m.eta = m.outer_eta;  m.phi = m.outer_phi;  m.charge = m.outer_charge; }
    m.trackPt = m.inner_pt;

    out.push_back(m);
  }
  return out;
}

template <>
inline std::vector<GenericMuon> MuonCollectionAdapter<l1t::TrackerMuonCollection>::extract() const {
  std::vector<GenericMuon> out;
  out.reserve(handle_->size());
  for (const auto& mu : *handle_) {
    GenericMuon m;
    m.pt = mu.phPt();
    m.eta = mu.phEta();
    m.phi = mu.phPhi();
    m.charge = mu.phCharge();
    out.push_back(m);
  }
  return out;
}

// -- Factory: turns the python-configured "type" string into the right
//    adapter. This is the ONLY switch statement in the whole mechanism.
inline std::unique_ptr<IMuonCollectionAdapter> makeMuonCollectionAdapter(const std::string& type,
                                                                          const edm::InputTag& tag,
                                                                          edm::ConsumesCollector&& cc) {
  if (type == "RecoChargedCandidate")
    return std::make_unique<MuonCollectionAdapter<reco::RecoChargedCandidateCollection>>(tag, std::move(cc));
  if (type == "TrackView")
    // edm::View<reco::Track> transparently wraps both edm::View-producing
    // modules AND plain std::vector<reco::Track>-based products (e.g. the
    // reco::TrackCollection produced for L2 muons), so this one adapter
    // covers both without needing a separate "Track" type.
    return std::make_unique<MuonCollectionAdapter<edm::View<reco::Track>>>(tag, std::move(cc));
  if (type == "RecoMuon")
    return std::make_unique<MuonCollectionAdapter<std::vector<reco::Muon>>>(tag, std::move(cc));
  if (type == "L1TkMuon")
    return std::make_unique<MuonCollectionAdapter<l1t::TrackerMuonCollection>>(tag, std::move(cc));
  if (type == "MuonTrackLinks")
    return std::make_unique<MuonCollectionAdapter<std::vector<reco::MuonTrackLinks>>>(tag, std::move(cc));

  throw cms::Exception("ConfigurationError")
      << "MuonHLTNtupler: unknown muonCollections 'type' = '" << type << "' for tag " << tag.encode() << ".\n"
      << "If this is a genuinely new EDM product (not just a renamed instance of an existing one), "
      << "add a MuonCollectionAdapter<> specialization + factory branch in GenericMuonCollections.h.";
}

#endif
