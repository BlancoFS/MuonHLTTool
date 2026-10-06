//////////////////////////////////////////////////////////////////////////
// HLTEfficiency.C
//
// Adapted from HLTBDTAnalyzer_binary.C for the new, self-describing ntuple
// schema (MuonHLTNtuple_v2.h). Main simplifications vs the original:
//
//  - No more one hardcoded get_XXX()/branch-set per named collection: the
//    reader discovers muonCollections/triggerFilters/triggerPaths from the
//    file itself, so this macro loops over whatever is actually configured
//    rather than a fixed L3types list that had to be hand-edited whenever
//    the HLT sequence changed.
//  - Trigger filter/path efficiency no longer needs its own dR-matching
//    pass per filter (the old get_HLTObjects(filterName) + manual match):
//    the ntupler already computed and stored per-object filter/path
//    matches, so this macro just gen-matches ONE reference reco collection
//    and then reads off whether that matched object passed each filter/path.
//  - The core Eff/Purity/Res machinery (HistContainer/ResContainer) is
//    unchanged so output files stay compatible with existing plotting
//    scripts. The "+L1" conditional variants (hc_Eff_L1Tk etc.) are kept,
//    generalized to whatever l1Collection/l1PtCut you configure, so you
//    can still factorize HLT-only vs HLT+L1 efficiency.
//
// Usage (unchanged from the original, same positional args your
// run_all_binary.sh already uses):
//   root -l -b -q 'HLTEfficiency.C("v00", "TEST", {"/path/to/ntuple*.root"})' >&log&
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include <TH1.h>
#include <TF1.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>

#include "MuonHLTNtuple_v2.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
// small utilities (unchanged from the original)
//////////////////////////////////////////////////////////////////////////
void printRunTime(TStopwatch timer_)
{
    cout << endl << "************************************************" << endl;
    cout << "Total real time: " << timer_.RealTime() << " (seconds)" << endl;
    cout << "Total CPU time:  " << timer_.CpuTime()  << " (seconds)" << endl;
    cout << "************************************************" << endl;
}

static inline void loadBar(int x, int n, int r, int w)
{
    if (x == n) cout << endl;
    if (x % (n / r + 1) != 0) return;
    float ratio = x / (float)n;
    int c = ratio * w;
    printf("%3d%% [", (int)(ratio * 100));
    for (int i = 0; i < c; i++) cout << "=";
    for (int i = c; i < w; i++) cout << " ";
    cout << "]\r" << flush;
}

bool acceptance(const Object& obj) { return fabs(obj.eta) < 2.4; }

// -- collections built from L1/L2 objects need a looser dR window (coarser
//    position resolution) and no relative-pt cut; everything downstream of
//    L3 tracking gets the tight window. Same heuristic as the original
//    code's substring checks, generalized to work on discovered names.
bool isLooseCollection(const TString& name)
{
    return name.Contains("L1") || name.Contains("L2");
}

//////////////////////////////////////////////////////////////////////////
// ResContainer / HistContainer -- unchanged from the original (already
// generic; kept identical so downstream plotting scripts still work).
//////////////////////////////////////////////////////////////////////////
class ResContainer
{
public:
    ResContainer(TString _Tag) : Tag(_Tag) { Init(); }
    void Init()
    {
        TH1::SetDefaultSumw2(kTRUE);
        TH1::AddDirectory(kFALSE);
        v_res = new TH1D(Tag, "", 500, -2.5, 2.5);
    }
    void Fill(double reco_pt, double gen_pt, double weight = 1.0)
    {
        v_res->Fill((reco_pt - gen_pt) / gen_pt, weight);
    }
    void Save(TDirectory* dir) { dir->cd(); v_res->SetDirectory(dir); v_res->Write(); }
    ~ResContainer() { delete v_res; }

    TH1D* v_res;
    TString Tag;
};

class HistContainer
{
public:
    HistContainer(
        TString _Tag,
        vector<TString> _variables = {"pt", "eta", "phi", "pu"},
        vector<vector<double>> _ranges = {
            {1000, 0, 1000}, {48, -2.4, 2.4}, {60, -TMath::Pi(), TMath::Pi()}, {250, 0, 250}
        }
    ) : Tag(_Tag), variables(_variables), ranges(_ranges)
    {
        if (variables.size() != ranges.size()) {
            cout << "HistContainer: variables.size() != ranges.size()" << endl;
            exit(1);
        }
        nVar = variables.size();
        Init();
    }

    void fill_den(Object& obj, double PU, double weight = 1.0)
    {
        for (int k = 0; k < nVar; ++k) {
            if (variables[k] == "pu") v_den[k]->Fill(PU, weight);
            else if (obj.has(variables[k])) v_den[k]->Fill(obj.get(variables[k]), weight);
        }
    }
    void fill_num(Object& obj, double PU, double weight = 1.0)
    {
        for (int k = 0; k < nVar; ++k) {
            if (variables[k] == "pu") v_num[k]->Fill(PU, weight);
            else if (obj.has(variables[k])) v_num[k]->Fill(obj.get(variables[k]), weight);
        }
    }
    void Save(TDirectory* dir)
    {
        dir->cd();
        for (int k = 0; k < nVar; ++k) {
            v_den[k]->SetDirectory(dir); v_num[k]->SetDirectory(dir);
            v_den[k]->Write(); v_num[k]->Write();
            delete v_den[k]; delete v_num[k];
        }
    }
    ~HistContainer() {}

private:
    TString Tag;
    int nVar;
    vector<TString> variables;
    vector<vector<double>> ranges;
    vector<TH1D*> v_den, v_num;

    void Init()
    {
        TH1::SetDefaultSumw2(kTRUE);
        TH1::AddDirectory(kFALSE);
        TString tag_ = (Tag == "") ? "" : "_" + Tag;
        for (int k = 0; k < nVar; ++k) {
            TString name = TString::Format("%s_%s", tag_.Data(), variables[k].Data());
            v_den.push_back(new TH1D("den" + name, "", ranges[k][0], ranges[k][1], ranges[k][2]));
            v_num.push_back(new TH1D("num" + name, "", ranges[k][0], ranges[k][1], ranges[k][2]));
        }
    }
};

// -- gen muon passes the fixed kinematic/status selection used throughout
//    ("probe" muon in the tag-and-probe/gen-matched sense: status==1 stable
//    muon, in acceptance, and -- if doDimuon -- actually from the hard
//    process DY decay rather than e.g. a heavy-flavor decay muon)
bool isGoodGenMuon(Object& genmu, vector<Object>& genMuonsFromHardProcess, bool doDimuon)
{
    if (fabs(genmu.get("ID")) != 13) return false;
    if (genmu.get("status") != 1) return false;
    if (!acceptance(genmu)) return false;
    if (doDimuon && !genmu.matched(genMuonsFromHardProcess, 0.001)) return false;
    return true;
}

//////////////////////////////////////////////////////////////////////////
void HLTEfficiency(
    TString ver = "v00", TString tag = "TEST",
    vector<TString> vec_Dataset = {}, TString JobId = "",
    TString refCollection = "L3Muon",   // reference reco collection for filter/path efficiency
    const bool doDimuon = false,
    TString l1Collection = "L1TkMuon",  // collection used for the "+L1" conditional histograms
    double l1PtCut = 22.0,              // pt cut defining the L1 seed used to factorize HLT-only vs HLT+L1
    double PU_min = -1, double PU_max = 1e6,
    Int_t maxEv = -1, bool doBar = true
) {
    TH1::SetDefaultSumw2(kTRUE);
    TH1::AddDirectory(kFALSE);

    TStopwatch timer_total; timer_total.Start();

    vector<TString> paths = vec_Dataset;
    if (tag == "TEST") paths = {"./ntuple*.root"};

    TString fileName = TString::Format("hist-%s-%s", ver.Data(), tag.Data());
    if (PU_min >= 0.) fileName += TString::Format("-PU%.0fto%.0f", PU_min, PU_max);
    if (JobId != "")  fileName += TString::Format("--%s", JobId.Data());
    TFile* f_output = TFile::Open("./Phase2_Muon_" + fileName + ".root", "RECREATE");

    TChain* _chain_Ev = new TChain("ntupler/ntuple");
    for (auto& p : paths) { _chain_Ev->Add(p); cout << "Adding path: " << p << endl; }

    unsigned nEvent = _chain_Ev->GetEntries();
    if (maxEv >= 0) nEvent = maxEv;
    cout << "\t nEvent: " << nEvent << endl;

    // -- no branch_tags whitelist needed: MuonHLTNtuple_v2 only binds the
    //    branches it actually reads.
    unique_ptr<MuonHLTNtuple_v2> nt(new MuonHLTNtuple_v2(_chain_Ev));

    const vector<TString>& collNames   = nt->muonCollectionNames();
    const vector<TString>& filterNames = nt->triggerFilterNames();
    const vector<TString>& pathNames   = nt->triggerPathNames();

    int refIdx = -1;
    for (size_t i = 0; i < collNames.size(); ++i) if (collNames[i] == refCollection) refIdx = (int)i;
    if (refIdx < 0) {
        cout << "HLTEfficiency: refCollection '" << refCollection << "' not found among configured "
             << "muonCollections -- filter/path efficiency will be skipped. Check the name against "
             << "what's printed above." << endl;
    }

    // -- basic gen-level bookkeeping histograms (unchanged from the original)
    TH1D* h_nEvents = new TH1D("h_nEvents", "", 3, -1, 2);
    TH1D* h_gen_pt  = new TH1D("h_gen_pt",  "", 1000, 0, 1000);
    TH1D* h_gen_eta = new TH1D("h_gen_eta", "", 60, -3, 3);
    TH1D* h_gen_acc_pt  = new TH1D("h_gen_acc_pt",  "", 1000, 0, 1000);
    TH1D* h_gen_acc_eta = new TH1D("h_gen_acc_eta", "", 60, -3, 3);

    // -- efficiency binning thresholds (same defaults as the original)
    vector<double> Eff_genpt_mins  = {0, 26, 53};
    vector<double> Eff_L3pt_mins   = {0, 8, 22, 24, 50};
    vector<double> Purity_L3pt_mins = {0, 26, 53};

    // -- per-collection reconstruction efficiency: Eff[collection][gen pt min]
    //    "_L1Tk" variants condition on the gen muon also matching an
    //    l1Collection object above l1PtCut, factorizing HLT-only efficiency
    //    (hc_Eff) from HLT+L1 efficiency (hc_Eff_L1Tk).
    vector<vector<HistContainer*>> hc_Eff, hc_EffTO, hc_Eff_L1Tk, hc_EffTO_L1Tk;
    for (auto& name : collNames) {
        vector<HistContainer*> row0, row1, row0L1, row1L1;
        for (auto& ptmin : Eff_genpt_mins) {
            row0.push_back(new HistContainer(TString::Format("Eff_%s_genpt%.0f", name.Data(), ptmin)));
            row0L1.push_back(new HistContainer(TString::Format("Eff_L1Tk_%s_genpt%.0f", name.Data(), ptmin)));
        }
        for (auto& ptmin : Eff_L3pt_mins) {
            row1.push_back(new HistContainer(TString::Format("EffTO_%s_L3pt%.0f", name.Data(), ptmin)));
            row1L1.push_back(new HistContainer(TString::Format("EffTO_L1Tk_%s_L3pt%.0f", name.Data(), ptmin)));
        }
        hc_Eff.push_back(row0);       hc_Eff_L1Tk.push_back(row0L1);
        hc_EffTO.push_back(row1);     hc_EffTO_L1Tk.push_back(row1L1);
    }

    // -- per-filter / per-path efficiency wrt refCollection: Eff[filter][gen pt min]
    vector<vector<HistContainer*>> hc_FilterEff, hc_FilterEff_L1Tk, hc_PathEff, hc_PathEff_L1Tk;
    for (auto& name : filterNames) {
        vector<HistContainer*> row, rowL1;
        for (auto& ptmin : Eff_genpt_mins) {
            row.push_back(new HistContainer(TString::Format("Eff_filter_%s_genpt%.0f", name.Data(), ptmin)));
            rowL1.push_back(new HistContainer(TString::Format("Eff_L1Tk_filter_%s_genpt%.0f", name.Data(), ptmin)));
        }
        hc_FilterEff.push_back(row);
        hc_FilterEff_L1Tk.push_back(rowL1);
    }
    for (auto& name : pathNames) {
        vector<HistContainer*> row, rowL1;
        for (auto& ptmin : Eff_genpt_mins) {
            row.push_back(new HistContainer(TString::Format("Eff_path_%s_genpt%.0f", name.Data(), ptmin)));
            rowL1.push_back(new HistContainer(TString::Format("Eff_L1Tk_path_%s_genpt%.0f", name.Data(), ptmin)));
        }
        hc_PathEff.push_back(row);
        hc_PathEff_L1Tk.push_back(rowL1);
    }

    // -- per-collection purity + pt resolution. "_L1Tk" purity conditions
    //    on the RECO object itself (not the gen muon) matching l1Collection
    //    above l1PtCut -- same convention as the original.
    vector<vector<HistContainer*>> hc_Purity_Sig, hc_Purity_Sig_L1Tk;
    vector<vector<ResContainer*>>  hc_Res;
    for (auto& name : collNames) {
        vector<HistContainer*> row0, row0L1;
        vector<ResContainer*>  row1;
        for (auto& ptmin : Purity_L3pt_mins) {
            row0.push_back(new HistContainer(TString::Format("Purity_%s_L3pt%.0f", name.Data(), ptmin)));
            row0L1.push_back(new HistContainer(TString::Format("Purity_L1Tk_%s_L3pt%.0f", name.Data(), ptmin)));
            row1.push_back(new ResContainer(TString::Format("Res_%s_L3pt%.0f", name.Data(), ptmin)));
        }
        hc_Purity_Sig.push_back(row0);
        hc_Purity_Sig_L1Tk.push_back(row0L1);
        hc_Res.push_back(row1);
    }

    int l1Idx = -1;
    for (size_t i = 0; i < collNames.size(); ++i) if (collNames[i] == l1Collection) l1Idx = (int)i;
    if (l1Idx < 0) {
        cout << "HLTEfficiency: l1Collection '" << l1Collection << "' not found among configured "
             << "muonCollections -- all '_L1Tk' histograms will be filled with zero entries." << endl;
    }


    cout << "Evt loop start" << endl;
    for (unsigned i_ev = 0; i_ev < nEvent; ++i_ev) {
        if (doBar) loadBar(i_ev + 1, nEvent, 100, 100);

        nt->GetEntry(i_ev);
        double genWeight = nt->genEventWeight > 0.0 ? 1.0 : -1.0;
        h_nEvents->Fill(genWeight);

        vector<Object> GenParticles = nt->get_GenParticles();

        vector<Object> genMuonsFromHardProcess;
        bool found0 = false, found1 = false;
        for (auto& g : GenParticles) {
            if (fabs(g.get("ID")) == 13 && g.get("fromHardProcessFinalState") == 1)
                genMuonsFromHardProcess.push_back(g);
            if (g.get("ID") == 13  && g.get("isHardProcess") == 1) found0 = true;
            if (g.get("ID") == -13 && g.get("isHardProcess") == 1) found1 = true;

            if (fabs(g.get("ID")) == 13 && g.get("status") == 1) {
                h_gen_pt->Fill(g.pt, genWeight);
                h_gen_eta->Fill(g.eta, genWeight);
                if (acceptance(g)) { h_gen_acc_pt->Fill(g.pt, genWeight); h_gen_acc_eta->Fill(g.eta, genWeight); }
            }
        }
        bool isDimuon = found0 && found1;
        if (doDimuon && !isDimuon) continue;

        // -- fetch every discovered collection once per event
        map<TString, vector<Object>> collections;
        for (auto& name : collNames) collections[name] = nt->get_Collection(name);

        // -- L1 seed muons (l1Collection above l1PtCut) used to factorize
        //    HLT-only vs HLT+L1 efficiency in every "_L1Tk" histogram below.
        vector<Object> l1SeedMuons;
        if (l1Idx >= 0)
            for (auto& l1mu : collections[collNames[l1Idx]])
                if (l1mu.pt > l1PtCut) l1SeedMuons.push_back(l1mu);

        // ---------------------------------------------------------------
        // Reconstruction efficiency, per collection
        // ---------------------------------------------------------------
        for (size_t i = 0; i < collNames.size(); ++i) {
            vector<Object>& coll = collections[collNames[i]];
            bool loose = isLooseCollection(collNames[i]);
            vector<int> used(coll.size(), 0);

            for (auto& genmu : GenParticles) {
                if (!isGoodGenMuon(genmu, genMuonsFromHardProcess, doDimuon)) continue;

                int matchedIdx = loose ? genmu.matched(coll, used, 0.3)
                                       : genmu.matched(coll, used, 0.1, 0.3);
                bool matchedL1 = genmu.matched(l1SeedMuons, 0.3);

                for (size_t j = 0; j < Eff_genpt_mins.size(); ++j) {
                    if (genmu.pt <= Eff_genpt_mins[j]) continue;
                    hc_Eff[i][j]->fill_den(genmu, nt->truePU, genWeight);
                    if (matchedIdx >= 0) hc_Eff[i][j]->fill_num(genmu, nt->truePU, genWeight);

                    if (matchedL1) {
                        hc_Eff_L1Tk[i][j]->fill_den(genmu, nt->truePU, genWeight);
                        if (matchedIdx >= 0) hc_Eff_L1Tk[i][j]->fill_num(genmu, nt->truePU, genWeight);
                    }
                }
                for (size_t j = 0; j < Eff_L3pt_mins.size(); ++j) {
                    hc_EffTO[i][j]->fill_den(genmu, nt->truePU, genWeight);
                    if (matchedIdx >= 0 && coll[matchedIdx].pt > Eff_L3pt_mins[j])
                        hc_EffTO[i][j]->fill_num(genmu, nt->truePU, genWeight);

                    if (matchedL1) {
                        hc_EffTO_L1Tk[i][j]->fill_den(genmu, nt->truePU, genWeight);
                        if (matchedIdx >= 0 && coll[matchedIdx].pt > Eff_L3pt_mins[j])
                            hc_EffTO_L1Tk[i][j]->fill_num(genmu, nt->truePU, genWeight);
                    }
                }
            }
        }

        // ---------------------------------------------------------------
        // Filter / path efficiency wrt refCollection -- this is the part
        // that's essentially free now: no per-filter dR matching needed,
        // just gen-match ONE reference collection and read off whether the
        // matched object already has the filter/path index in its list.
        // ---------------------------------------------------------------
        if (refIdx >= 0) {
            vector<Object>& refColl = collections[collNames[refIdx]];
            bool loose = isLooseCollection(collNames[refIdx]);
            vector<int> used(refColl.size(), 0);

            for (auto& genmu : GenParticles) {
                if (!isGoodGenMuon(genmu, genMuonsFromHardProcess, doDimuon)) continue;

                int matchedIdx = loose ? genmu.matched(refColl, used, 0.3)
                                       : genmu.matched(refColl, used, 0.1, 0.3);
                bool matched = (matchedIdx >= 0);
                Object* matchedObj = matched ? &refColl[matchedIdx] : nullptr;
                bool matchedL1 = genmu.matched(l1SeedMuons, 0.3);

                for (size_t j = 0; j < Eff_genpt_mins.size(); ++j) {
                    if (genmu.pt <= Eff_genpt_mins[j]) continue;

                    for (size_t k = 0; k < filterNames.size(); ++k) {
                        bool passNum = matched && matchedObj->passedFilter((int)k);
                        hc_FilterEff[k][j]->fill_den(genmu, nt->truePU, genWeight);
                        if (passNum) hc_FilterEff[k][j]->fill_num(genmu, nt->truePU, genWeight);

                        if (matchedL1) {
                            hc_FilterEff_L1Tk[k][j]->fill_den(genmu, nt->truePU, genWeight);
                            if (passNum) hc_FilterEff_L1Tk[k][j]->fill_num(genmu, nt->truePU, genWeight);
                        }
                    }
                    for (size_t k = 0; k < pathNames.size(); ++k) {
                        bool passNum = matched && matchedObj->passedPath((int)k);
                        hc_PathEff[k][j]->fill_den(genmu, nt->truePU, genWeight);
                        if (passNum) hc_PathEff[k][j]->fill_num(genmu, nt->truePU, genWeight);

                        if (matchedL1) {
                            hc_PathEff_L1Tk[k][j]->fill_den(genmu, nt->truePU, genWeight);
                            if (passNum) hc_PathEff_L1Tk[k][j]->fill_num(genmu, nt->truePU, genWeight);
                        }
                    }
                }
            }
        }

        // ---------------------------------------------------------------
        // Purity + resolution, per collection. "_L1Tk" purity conditions on
        // the RECO object mu itself matching an l1SeedMuon (not the gen
        // muon) -- same convention the original used.
        // ---------------------------------------------------------------
        for (size_t i = 0; i < collNames.size(); ++i) {
            vector<Object>& coll = collections[collNames[i]];
            bool loose = isLooseCollection(collNames[i]);

            for (auto& mu : coll) {
                bool muMatchedL1 = mu.matched(l1SeedMuons, 0.3);

                for (size_t j = 0; j < Purity_L3pt_mins.size(); ++j) {
                    if (mu.pt > Purity_L3pt_mins[j]) {
                        hc_Purity_Sig[i][j]->fill_den(mu, nt->truePU, genWeight);
                        if (muMatchedL1) hc_Purity_Sig_L1Tk[i][j]->fill_den(mu, nt->truePU, genWeight);
                    }
                }

                for (auto& genmu : GenParticles) {
                    if (!isGoodGenMuon(genmu, genMuonsFromHardProcess, doDimuon)) continue;
                    bool matched = loose ? mu.matched(genmu, 0.3) : mu.matched(genmu, 0.1, 0.3);
                    if (!matched) continue;

                    for (size_t j = 0; j < Purity_L3pt_mins.size(); ++j) {
                        if (mu.pt > Purity_L3pt_mins[j]) {
                            hc_Purity_Sig[i][j]->fill_num(mu, nt->truePU, genWeight);
                            hc_Res[i][j]->Fill(mu.pt, genmu.pt, genWeight);
                            if (muMatchedL1) hc_Purity_Sig_L1Tk[i][j]->fill_num(mu, nt->truePU, genWeight);
                        }
                    }
                    break;  // -- one-to-one: stop at the first gen match
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // Save
    // -----------------------------------------------------------------
    f_output->cd();
    h_nEvents->Write();
    h_gen_pt->Write(); h_gen_eta->Write();
    h_gen_acc_pt->Write(); h_gen_acc_eta->Write();

    TDirectory* dirEff        = f_output->mkdir("Eff");
    TDirectory* dirEffFilters = dirEff->mkdir("Filters");
    TDirectory* dirEffPaths   = dirEff->mkdir("Paths");
    TDirectory* dirPur        = f_output->mkdir("Pur");
    TDirectory* dirRes        = f_output->mkdir("Res");

    for (size_t i = 0; i < collNames.size(); ++i) {
        for (size_t j = 0; j < Eff_genpt_mins.size(); ++j) {
            hc_Eff[i][j]->Save(dirEff);      delete hc_Eff[i][j];
            hc_Eff_L1Tk[i][j]->Save(dirEff); delete hc_Eff_L1Tk[i][j];
        }
        for (size_t j = 0; j < Eff_L3pt_mins.size(); ++j) {
            hc_EffTO[i][j]->Save(dirEff);      delete hc_EffTO[i][j];
            hc_EffTO_L1Tk[i][j]->Save(dirEff); delete hc_EffTO_L1Tk[i][j];
        }
        for (size_t j = 0; j < Purity_L3pt_mins.size(); ++j) {
            hc_Purity_Sig[i][j]->Save(dirPur);      delete hc_Purity_Sig[i][j];
            hc_Purity_Sig_L1Tk[i][j]->Save(dirPur); delete hc_Purity_Sig_L1Tk[i][j];
            hc_Res[i][j]->Save(dirRes);             delete hc_Res[i][j];
        }
    }
    for (size_t k = 0; k < filterNames.size(); ++k)
        for (size_t j = 0; j < Eff_genpt_mins.size(); ++j) {
            hc_FilterEff[k][j]->Save(dirEffFilters);      delete hc_FilterEff[k][j];
            hc_FilterEff_L1Tk[k][j]->Save(dirEffFilters); delete hc_FilterEff_L1Tk[k][j];
        }
    for (size_t k = 0; k < pathNames.size(); ++k)
        for (size_t j = 0; j < Eff_genpt_mins.size(); ++j) {
            hc_PathEff[k][j]->Save(dirEffPaths);      delete hc_PathEff[k][j];
            hc_PathEff_L1Tk[k][j]->Save(dirEffPaths); delete hc_PathEff_L1Tk[k][j];
        }

    f_output->Close();
    delete f_output;

    printRunTime(timer_total);
}
