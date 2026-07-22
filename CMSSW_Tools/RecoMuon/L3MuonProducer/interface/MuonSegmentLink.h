#ifndef MuonHLTSegmentMatcher_MuonSegmentLink_h
#define MuonHLTSegmentMatcher_MuonSegmentLink_h

// -------------------------------------------------------------------------
// Small standalone header shared by:
//   - MuonHLTSegmentMatcherProducer  (writes the links)
//   - MuonHLTMuonFromLinksProducer   (reads the links)
//
// Kept torch-free so consumers of the links don't need LibTorch.
//
// Semantics
// ---------
// For a track-segment link (produced under label "trackSegmentLinks"):
//   trackIdx    = index into the input track collection
//   segmentIdx  = index into the concatenated [DT | CSC | GEM] segment list
//
// For a segment-segment link (produced under label "segmentSegmentLinks"):
//   trackIdx    = index of the FIRST  segment in the concatenated list
//   segmentIdx  = index of the SECOND segment in the concatenated list
//
// The fields are named after the more common track-segment case; the
// segment-segment producer overloads their meaning to keep a single POD.
// -------------------------------------------------------------------------

struct MuonSegmentLink {
  int   trackIdx;      // track index  OR  first segment index (SS links)
  int   segmentIdx;    // segment index OR second segment index (SS links)
  float score;         // sigmoid(logit), higher = more likely same muon
};

#endif
