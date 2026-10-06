#ifndef MuonHLTTriggerMatching_h
#define MuonHLTTriggerMatching_h

// -----------------------------------------------------------------------
// Purpose
// -----------------------------------------------------------------------
// Today, Fill_HLT() in MuonHLTNtupler.cc dumps EVERY fired path and EVERY
// filter's trigger objects into flat vectors (the SavedTriggerCondition /
// SavedFilterCondition gates are both hardcoded to `if(true)`), with no
// per-candidate dR matching. Your Run-2 MuonMiniAODAnalyzer had the right
// idea (HLTaccept + embedTriggerMatching: configurable path/filter list,
// dR-matched), but it's wired specifically to tag/probe pat::Muon objects.
//
// This header generalizes that idea: it takes the *configured* list of
// path names / filter labels from python, resolves versioned path names
// once per run, and exposes a plain dR-match function that works against
// ANY candidate exposing pt()/eta()/phi() -- an L1 muon, an L2 muon, an
// offline muon, whatever collection GenericMuonCollections.h handed you.
// There is no tag/probe distinction baked in here; the caller decides what
// "tag" or "probe" means, if anything, for its own analysis.
// -----------------------------------------------------------------------

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TString.h"

#include <string>
#include <vector>

namespace MuonHLTTriggerMatching {

  // -- Resolve versioned path names ("HLT_Mu50_v12") from unversioned
  //    python config entries ("HLT_Mu50_v"). Call once in beginRun() after
  //    hltConfig.init(...). This means the python trigger list never needs
  //    to be touched when the menu bumps a trigger version.
  inline std::vector<std::string> resolvePathNames(const HLTConfigProvider& hltConfig,
                                                     const std::vector<std::string>& configuredPaths) {
    std::vector<std::string> resolved;
    resolved.reserve(configuredPaths.size());
    for (const auto& configured : configuredPaths) {
      TString configuredT(configured);
      bool found = false;
      for (const auto& actual : hltConfig.triggerNames()) {
        if (TString(actual).Contains(configuredT)) {
          resolved.push_back(actual);
          found = true;
          break;
        }
      }
      // keep the configured name even if unresolved -- it simply won't be
      // found in TriggerResults this run, rather than crashing the job.
      if (!found)
        resolved.push_back(configured);
    }
    return resolved;
  }

  struct FilterObject {
    float pt, eta, phi;
  };

  // -- Pull, for ONE configured filter label, all trigger objects it fired
  //    on in this event. Computed once per filter per event and reused
  //    against every candidate collection -- O(nFilters) TriggerEvent
  //    scans per event instead of O(nFilters * nCandidates).
  inline std::vector<FilterObject> getFilterObjects(const trigger::TriggerEvent& triggerEvent,
                                                      const std::string& filterLabel) {
    std::vector<FilterObject> objs;
    const trigger::TriggerObjectCollection& all = triggerEvent.getObjects();
    for (trigger::size_type iFilter = 0; iFilter < triggerEvent.sizeFilters(); ++iFilter) {
      if (triggerEvent.filterTag(iFilter).label() != filterLabel)
        continue;
      for (auto key : triggerEvent.filterKeys(iFilter)) {
        objs.push_back({all[key].pt(), all[key].eta(), all[key].phi()});
      }
    }
    return objs;
  }

  // -- Generic dR match of ANY candidate (eta/phi) against one filter's
  //    objects. Works identically for L1, L2, L3, TkMuon, offline muons --
  //    whatever GenericMuon you hand it. Returns false (bestDR = 99) if
  //    there is no object within maxDR.
  inline bool matchToFilter(float eta,
                             float phi,
                             const std::vector<FilterObject>& filterObjs,
                             float maxDR,
                             float& bestDR,
                             float& bestPt) {
    bestDR = 99.f;
    bestPt = -99.f;
    for (const auto& obj : filterObjs) {
      float dr = reco::deltaR(eta, phi, obj.eta, obj.phi);
      if (dr < bestDR) {
        bestDR = dr;
        bestPt = obj.pt;
      }
    }
    return bestDR < maxDR;
  }

  // -- Find the last EDFilter-type module in a path's (resolved, versioned)
  //    module sequence. An object is conventionally said to have "fired"
  //    a path if it matches this filter's trigger objects -- passing the
  //    final filter of a leg implies everything upstream of it in that leg
  //    was already satisfied. This is the same convention used for HLT
  //    trigger-object-to-path embedding in miniAOD.
  inline std::string lastFilterOfPath(const HLTConfigProvider& hltConfig, const std::string& resolvedPathName) {
    const auto& modules = hltConfig.moduleLabels(resolvedPathName);
    for (auto it = modules.rbegin()+1; it != modules.rend(); ++it) {      
      if (hltConfig.moduleEDMType(*it) == "EDFilter") {
        return *it;
      }
    }
    return "";  // path not found or has no EDFilter module -- callers treat as "never matches"
  }

}  // namespace MuonHLTTriggerMatching

#endif
