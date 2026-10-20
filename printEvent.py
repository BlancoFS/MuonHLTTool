import uproot
import awkward as ak

fname = "muonHLT_ntuple.root"

branches = [
    "muon_pt","muon_eta","muon_phi",
    "L3Muon_pt","L3Muon_eta","L3Muon_phi",
    "L1TkMuon_pt","L1TkMuon_eta","L1TkMuon_phi",
    "muon_matchChannel",
    "muon_muSysSimId",
    "muon_muSysHitPurity",
    "muon_matched",
    "muon_simId",
    "muon_simPdgId",
    "muon_simHasGen",    
    "muon_simPt",
    "muon_simEta", 
    "muon_simPhi", 
    "muon_sharedHit",
    "muon_hitPurity",
    "muon_fromW_", 
    "muon_fromZ",
    "muon_fromTau",
    "muon_fromLightHadron",
    "muon_fromB",
    "muon_fromC",
    "muon_fromBSM",
    "muon_simMotherPdgId",
    "muon_isGenerator",    
    "muon_isPrompt",
]

with uproot.open(fname) as f:
    tree = f["ntupler/ntuple"]
    missing = [b for b in branches if b not in tree.keys()]
    if missing:
        print("Missing branches:", missing)
    events = tree.arrays([b for b in branches if b not in missing], library="ak")

print(f"{len(events)} events")
for i, ev in enumerate(events):
    print(f"--- event {i} ---")
    for b in events.fields:
        print(f"  {b:14s}: {ev[b].tolist()}")
