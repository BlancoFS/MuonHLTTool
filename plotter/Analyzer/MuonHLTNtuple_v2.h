//////////////////////////////////////////////////////////////////////////
// MuonHLTNtuple_v2.h
//
// Generic reader for the "configuredMuonCollections" / "configuredTriggerFilters"
// / "configuredTriggerPaths" ntuple schema produced by the refactored
// MuonHLTNtupler. Unlike the old MuonHLTNtuple_binary.h -- which hardcoded
// one get_XXX() method and one set of branch pointers per named collection,
// and required re-editing this file whenever the HLT sequence changed --
// this reader discovers which collections/filters/paths exist by reading
// the three "configured*" branches from the first entry of the chain, then
// binds branches dynamically. Adding, removing, or renaming a collection,
// filter, or path in ntupler_cfi.py requires ZERO changes here.
//
// ASSUMPTION: every file in the TChain was produced with the SAME
// muonCollections/triggerFilters/triggerPaths configuration (true for a
// normal single-campaign/single-JobId processing, as in run_all_binary.sh).
// Only the first file's lists are read; mixing incompatible configurations
// in one chain will silently misbehave (branches from later files simply
// won't match what was bound from the first file).
//////////////////////////////////////////////////////////////////////////

#ifndef MuonHLTNtuple_v2_h
#define MuonHLTNtuple_v2_h

#define ArrSize 50000

#include <TROOT.h>
#include <TChain.h>
#include <TString.h>

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

//////////////////////////////////////////////////////////////////////////
// Object: same interface as the original binary-era wrapper's Object
// (pt/eta/phi, addVar/get, deltaR/matched helpers) so existing matching
// logic ports over unchanged. Two additions for objects coming from a
// muonCollections entry: matchedFilterIdx/matchedPathIdx (+ parallel dR
// vectors), and passedFilter()/passedPath() to query them. Empty/unused
// for GenParticle objects.
//////////////////////////////////////////////////////////////////////////
class Object
{
private:
    map<TString, double> vars;

public:
    double pt, eta, phi;

    vector<int>   matchedFilterIdx;
    vector<float> matchedFilterDR;
    vector<int>   matchedPathIdx;
    vector<float> matchedPathDR;

    Object() { pt = -99999.; eta = -99999.; phi = -99999.; }

    Object(double _pt, double _eta, double _phi) : Object()
    {
        pt = _pt; eta = _eta; phi = _phi;
    }

    void addVar(TString name, double value) { vars[name] = value; }

    bool has(TString key) { return (vars.find(key) != vars.end()); }

    double get(TString key)
    {
        if (vars.find(key) == vars.end()) {
            cout << key << " does not exist -> return -99999" << endl;
            return -99999.0;
        }
        return vars[key];
    }

    // -- did this object match configured filter/path index idx (as
    //    produced by MuonHLTNtupler's genTrigMatchedFilterIdx_/
    //    genPathMatchedIdx_)? idx comes from
    //    MuonHLTNtuple_v2::filterIndex()/pathIndex().
    bool passedFilter(int idx) const
    {
        if (idx < 0) return false;
        for (int m : matchedFilterIdx) if (m == idx) return true;
        return false;
    }
    bool passedPath(int idx) const
    {
        if (idx < 0) return false;
        for (int m : matchedPathIdx) if (m == idx) return true;
        return false;
    }

    Object clone()
    {
        Object out;
        out.pt = pt; out.eta = eta; out.phi = phi;
        out.vars = vars;
        out.matchedFilterIdx = matchedFilterIdx;
        out.matchedFilterDR  = matchedFilterDR;
        out.matchedPathIdx   = matchedPathIdx;
        out.matchedPathDR    = matchedPathDR;
        return out;
    }

    static double reduceRange(double x)
    {
        double o2pi = 1. / (2. * M_PI);
        if (std::abs(x) <= double(M_PI)) return x;
        double n = std::round(x * o2pi);
        return x - n * double(2. * M_PI);
    }

    double deltaR(double _eta, double _phi) const
    {
        double dphi = reduceRange(this->phi - _phi);
        return sqrt((this->eta - _eta) * (this->eta - _eta) + dphi * dphi);
    }
    double deltaR(const Object& other) const { return deltaR(other.eta, other.phi); }

    bool matched(const Object& other, double dR_match = 0.1, double dpt_match = 1.e9) const
    {
        double dR  = deltaR(other);
        double dpt = fabs(this->pt - other.pt) / this->pt;
        return (dR < dR_match) && (dpt < dpt_match);
    }

    // -- match against a collection, marking off whichever object is used
    //    (via `used`) so two probes can't both claim the same reco object.
    //    Returns the matched index, or -1 if none.
    int matched(const vector<Object>& objects, vector<int>& used, double dR_match = 0.1, double dpt_match = 1.e9) const
    {
        double bestDR = dR_match;
        int bestI = -1;
        for (size_t i = 0; i < objects.size(); ++i) {
            if (used[i] > 0) continue;
            double dR  = deltaR(objects[i]);
            double dpt = fabs(this->pt - objects[i].pt) / this->pt;
            if (dR < bestDR && dpt < dpt_match) { bestDR = dR; bestI = (int)i; }
        }
        if (bestI >= 0) used[bestI] = 1;
        return bestI;
    }

    bool matched(const vector<Object>& objects, double dR_match = 0.1, double dpt_match = 1.e9) const
    {
        for (auto& o : objects)
            if (deltaR(o) < dR_match && fabs(this->pt - o.pt) / this->pt < dpt_match) return true;
        return false;
    }

    friend ostream& operator<<(ostream& os, const Object& obj)
    {
        os << "(" << obj.pt << ", " << obj.eta << ", " << obj.phi << ")";
        return os;
    }
};

//////////////////////////////////////////////////////////////////////////
// MuonHLTNtuple_v2
//////////////////////////////////////////////////////////////////////////
class MuonHLTNtuple_v2 {
public:
    TChain  *fChain;
    Int_t    fCurrent;

    MuonHLTNtuple_v2(TChain *tree);
    virtual ~MuonHLTNtuple_v2();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual Bool_t   Notify();

    // -- discovered from the file itself; no hardcoded collection/filter/
    //    path list to maintain.
    const vector<TString>& muonCollectionNames() const { return v_muonCollectionNames; }
    const vector<TString>& triggerFilterNames()  const { return v_triggerFilterNames; }
    const vector<TString>& triggerPathNames()    const { return v_triggerPathNames; }

    int filterIndex(TString name) const;  // -1 if not configured this run
    int pathIndex(TString name) const;

    vector<Object> get_GenParticles();
    vector<Object> get_Collection(TString name);  // works for ANY discovered muonCollections entry

    // -- event-level scalars
    Int_t    truePU;
    Double_t genEventWeight;

private:
    void bindGenAndScalarBranches();
    void bindCollectionBranches();

    vector<TString> v_muonCollectionNames, v_triggerFilterNames, v_triggerPathNames;

    vector<string> *p_muonCollectionNames = nullptr;
    vector<string> *p_triggerFilterNames  = nullptr;
    vector<string> *p_triggerPathNames    = nullptr;

    // -- one entry per discovered collection name (string key: TString's
    //    ordering operators are best avoided as a map key across ROOT
    //    versions, so we key on std::string internally)
    map<string, vector<float>*>         m_pt, m_eta, m_phi, m_charge, m_trkPt;
    map<string, vector<vector<int>>*>   m_trigIdx, m_pathIdx;
    map<string, vector<vector<float>>*> m_trigDR,  m_pathDR;

    // -- GenParticle: fixed-size legacy arrays, untouched by the
    //    muonCollections refactor, trimmed to only the fields the
    //    efficiency analyzer actually uses (status/ID/hard-process flags
    //    + kinematics). Extend here if you need more.
    Int_t    nGenParticle;
    Int_t    genParticle_ID[ArrSize];
    Int_t    genParticle_status[ArrSize];
    Int_t    genParticle_isHardProcess[ArrSize];
    Int_t    genParticle_fromHardProcessFinalState[ArrSize];
    Double_t genParticle_pt[ArrSize];
    Double_t genParticle_eta[ArrSize];
    Double_t genParticle_phi[ArrSize];
    Double_t genParticle_charge[ArrSize];
};

MuonHLTNtuple_v2::MuonHLTNtuple_v2(TChain *tree) : fChain(0), fCurrent(-1)
{
    fChain = tree;
    fChain->SetBranchStatus("*", 0);

    // -- Phase 1: bootstrap. Discover collection/filter/path names from
    //    the first entry before binding anything else.
    fChain->SetBranchStatus("configuredMuonCollections", 1);
    fChain->SetBranchStatus("configuredTriggerFilters",  1);
    fChain->SetBranchStatus("configuredTriggerPaths",    1);
    fChain->SetBranchAddress("configuredMuonCollections", &p_muonCollectionNames);
    fChain->SetBranchAddress("configuredTriggerFilters",  &p_triggerFilterNames);
    fChain->SetBranchAddress("configuredTriggerPaths",    &p_triggerPathNames);

    if (fChain->GetEntries() == 0) {
        cerr << "MuonHLTNtuple_v2: chain has zero entries, cannot discover configuration." << endl;
        return;
    }
    fChain->GetEntry(0);

    if (p_muonCollectionNames) for (auto& s : *p_muonCollectionNames) v_muonCollectionNames.push_back(s);
    if (p_triggerFilterNames)  for (auto& s : *p_triggerFilterNames)  v_triggerFilterNames.push_back(s);
    if (p_triggerPathNames)    for (auto& s : *p_triggerPathNames)    v_triggerPathNames.push_back(s);

    cout << "MuonHLTNtuple_v2: discovered " << v_muonCollectionNames.size() << " muon collections, "
         << v_triggerFilterNames.size() << " trigger filters, "
         << v_triggerPathNames.size() << " trigger paths." << endl;
    cout << "  Collections: ";
    for (auto& n : v_muonCollectionNames) cout << n << " ";
    cout << endl;

    // -- Phase 2: bind everything else now that we know the names.
    bindGenAndScalarBranches();
    bindCollectionBranches();
}

MuonHLTNtuple_v2::~MuonHLTNtuple_v2()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t MuonHLTNtuple_v2::GetEntry(Long64_t entry)
{
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t MuonHLTNtuple_v2::LoadTree(Long64_t entry)
{
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry >= 0 && fChain->GetTreeNumber() != fCurrent)
        fCurrent = fChain->GetTreeNumber();
    return centry;
}

Bool_t MuonHLTNtuple_v2::Notify()
{
    // See the file-level ASSUMPTION note: this does not re-validate that a
    // newly-loaded file in the chain shares the same configured
    // collections/filters/paths as the first one.
    return kTRUE;
}

void MuonHLTNtuple_v2::bindGenAndScalarBranches()
{
    fChain->SetBranchStatus("truePU", 1);
    fChain->SetBranchStatus("genEventWeight", 1);
    fChain->SetBranchStatus("nGenParticle", 1);
    fChain->SetBranchStatus("genParticle_ID", 1);
    fChain->SetBranchStatus("genParticle_status", 1);
    fChain->SetBranchStatus("genParticle_isHardProcess", 1);
    fChain->SetBranchStatus("genParticle_fromHardProcessFinalState", 1);
    fChain->SetBranchStatus("genParticle_pt", 1);
    fChain->SetBranchStatus("genParticle_eta", 1);
    fChain->SetBranchStatus("genParticle_phi", 1);
    fChain->SetBranchStatus("genParticle_charge", 1);

    fChain->SetBranchAddress("truePU", &truePU);
    fChain->SetBranchAddress("genEventWeight", &genEventWeight);
    fChain->SetBranchAddress("nGenParticle", &nGenParticle);
    fChain->SetBranchAddress("genParticle_ID", genParticle_ID);
    fChain->SetBranchAddress("genParticle_status", genParticle_status);
    fChain->SetBranchAddress("genParticle_isHardProcess", genParticle_isHardProcess);
    fChain->SetBranchAddress("genParticle_fromHardProcessFinalState", genParticle_fromHardProcessFinalState);
    fChain->SetBranchAddress("genParticle_pt", genParticle_pt);
    fChain->SetBranchAddress("genParticle_eta", genParticle_eta);
    fChain->SetBranchAddress("genParticle_phi", genParticle_phi);
    fChain->SetBranchAddress("genParticle_charge", genParticle_charge);
}

void MuonHLTNtuple_v2::bindCollectionBranches()
{
    for (auto& tname : v_muonCollectionNames) {
        string name = tname.Data();

        auto bind = [&](const char* suffix, auto*& ptr) {
            TString bname = tname + suffix;
            if (!fChain->GetBranch(bname)) {
                cerr << "MuonHLTNtuple_v2: expected branch '" << bname << "' not found -- skipping." << endl;
                return;
            }
            fChain->SetBranchStatus(bname, 1);
            fChain->SetBranchAddress(bname, &ptr);
        };

        bind("_pt",     m_pt[name]);
        bind("_eta",    m_eta[name]);
        bind("_phi",    m_phi[name]);
        bind("_charge", m_charge[name]);
        bind("_trkPt",  m_trkPt[name]);
        bind("_trigMatchedFilterIdx", m_trigIdx[name]);
        bind("_trigDR",               m_trigDR[name]);
        bind("_pathMatchedIdx",       m_pathIdx[name]);
        bind("_pathDR",               m_pathDR[name]);
    }
}

int MuonHLTNtuple_v2::filterIndex(TString name) const
{
    for (size_t i = 0; i < v_triggerFilterNames.size(); ++i)
        if (v_triggerFilterNames[i] == name) return (int)i;
    return -1;
}

int MuonHLTNtuple_v2::pathIndex(TString name) const
{
    for (size_t i = 0; i < v_triggerPathNames.size(); ++i)
        if (v_triggerPathNames[i] == name) return (int)i;
    return -1;
}

vector<Object> MuonHLTNtuple_v2::get_GenParticles()
{
    vector<Object> out;
    for (int i = 0; i < nGenParticle; ++i) {
        if (genParticle_status[i] != 1 && genParticle_isHardProcess[i] != 1)
            continue;

        Object obj(genParticle_pt[i], genParticle_eta[i], genParticle_phi[i]);
        obj.addVar("ID", genParticle_ID[i]);
        obj.addVar("status", genParticle_status[i]);
        obj.addVar("isHardProcess", genParticle_isHardProcess[i]);
        obj.addVar("fromHardProcessFinalState", genParticle_fromHardProcessFinalState[i]);
        obj.addVar("charge", genParticle_charge[i]);
        out.push_back(obj);
    }
    return out;
}

vector<Object> MuonHLTNtuple_v2::get_Collection(TString tname)
{
    vector<Object> out;
    string name = tname.Data();

    if (m_pt.find(name) == m_pt.end() || m_pt[name] == nullptr) {
        cerr << "MuonHLTNtuple_v2::get_Collection: unknown or unbound collection '" << tname << "'" << endl;
        return out;
    }

    auto* pt      = m_pt[name];
    auto* eta     = m_eta[name];
    auto* phi     = m_phi[name];
    auto* charge  = m_charge[name];
    auto* trkPt   = m_trkPt[name];
    auto* trigIdx = m_trigIdx[name];
    auto* trigDR  = m_trigDR[name];
    auto* pathIdx = m_pathIdx[name];
    auto* pathDR  = m_pathDR[name];

    for (size_t i = 0; i < pt->size(); ++i) {
        Object obj(pt->at(i), eta->at(i), phi->at(i));
        if (charge) obj.addVar("charge", charge->at(i));
        if (trkPt)  obj.addVar("trkPt",  trkPt->at(i));
        if (trigIdx && i < trigIdx->size()) obj.matchedFilterIdx = trigIdx->at(i);
        if (trigDR  && i < trigDR->size())  obj.matchedFilterDR  = trigDR->at(i);
        if (pathIdx && i < pathIdx->size()) obj.matchedPathIdx   = pathIdx->at(i);
        if (pathDR  && i < pathDR->size())  obj.matchedPathDR    = pathDR->at(i);
        out.push_back(obj);
    }
    return out;
}

#endif
