#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TSystem.h>

#include "MuonHLTNtuple.h"

using namespace std;


void DEBUG(int i, TString str = "")
{
    cout << "DEBUG: " << i << "\t" << str << endl;
}

void printRunTime(TStopwatch timer_)
{
    Double_t cpuTime = timer_.CpuTime();
    Double_t realTime = timer_.RealTime();

    cout << endl;
    cout << "************************************************" << endl;
    cout << "Total real time: " << realTime << " (seconds)" << endl;
    cout << "Total CPU time:  " << cpuTime << " (seconds)" << endl;
    cout << "  CPU time / real time = " << cpuTime / realTime << endl;
    cout << "************************************************" << endl;
}

void printRunTimeShort(TStopwatch timer_)
{
    // Double_t cpuTime = timer_.CpuTime();
    Double_t realTime = timer_.RealTime();

    cout << "Total real time: " << realTime << " (seconds)" << endl;
    // cout << "Total CPU time:  " << cpuTime << " (seconds)" << endl;
}

static inline void printMemory( TString tab = "" )
{
    ifstream proc_status("/proc/self/status");
    string buffer;
    while (proc_status.peek() != EOF) {
        getline(proc_status, buffer);
        TString str = buffer;
        if(str.Contains("RSS")) {
            cout << tab << str << endl;
            break;
        }
    }
}

static inline void loadBar(int x, int n, int r, int w)
{
    // Only update r times.
    if( x == n )
        cout << endl;

    if ( x % (n/r +1) != 0 ) return;

    // Calculuate the ratio of complete-to-incomplete.
    float ratio = x/(float)n;
    int   c     = ratio * w;

    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );

    // Show the load bar.
    for (int x=0; x<c; x++) cout << "=";

    for (int x=c; x<w; x++) cout << " ";

    // ANSI Control codes to go back to the
    // previous line and clear it.
    cout << "]\r" << flush;
}




const static int n_eta_bins = 15;
double eta_bins[n_eta_bins] = {
  -2.4, -2.1, -1.6, -1.2, -0.9,
  -0.3, -0.2,  0.0,  0.2,  0.3,
   0.9,  1.2,  1.6,  2.1,  2.4
};

// const static int n_pt_bins = 26;
// double pt_bins[n_pt_bins] = {
//    0, 10, 15, 16, 17,
//   18, 19, 20, 21, 22,
//   23, 24, 25, 26, 27,
//   30, 35, 40, 50, 60,
//   100, 150, 200, 300, 500,
//   1000
// };

const static int n_pt_bins = 16;
double pt_bins[n_pt_bins] = {
     0, 10, 12, 15, 19,
    25, 31, 39, 50, 63,
    79, 100, 140, 200, 500,
    1000
};

bool acceptance( Object obj )
{
    return ( fabs(obj.eta) < 2.4 );
}

struct sort_by_pt
{
    inline bool operator() (const Object& a, const Object& b)
    {
        return (a.pt > b.pt);
    }
};

// echo 'gROOT->LoadMacro("HLTResAnalyzer.C+"); gSystem->Exit(0);' | root -b -l
// root -l -b -q 'HLTResAnalyzer.C("v00", "TEST")' >&log&

void HLTResAnalyzer(
    TString ver = "v00", TString tag = "TEST",
    vector<TString> vec_Dataset = {}, TString JobId = "",
    const bool doDimuon = false, const bool doGenMatchForIso = false,
    double PU_min = -1, Double_t PU_max = 1e6,
    Int_t maxEv = -1, bool doMem = false, int nMem = 10001, bool doBar = true  // HERE
) {
    bool xxx = (doDimuon && doGenMatchForIso);  xxx = true;

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);

    TStopwatch timer_total;
    timer_total.Start();

    // -- Input
        vector<TString> paths = vec_Dataset;

        TString base_path = "";

        if(tag == "TEST") {
            base_path = "../";
            paths = { "./ntuple*.root" };
        }

    // -- Output
        TString fileName = TString::Format( "hist-%s-%s", ver.Data(), tag.Data() );
        if(PU_min >= 0.) fileName = fileName + TString::Format("-PU%.0fto%.0f", PU_min, PU_max);
        if(JobId != "")  fileName = fileName + TString::Format("--%s", JobId.Data());
        TFile *f_output = TFile::Open(fileName+"-Res.root", "RECREATE");

    // -- Event chain
        TChain *_chain_Ev          = new TChain("ntupler/ntuple");
        for(size_t f = 0; f<paths.size(); ++f) {
            _chain_Ev->Add(paths[f]);
            cout << "Adding path: " << paths[f] << endl;
        }
        cout << endl;

        unsigned nEvent      = _chain_Ev->GetEntries();
        if(maxEv >= 0)  nEvent = maxEv;
        cout << "\t nEvent: " << nEvent << endl;

        bool isIterL3 = false;
        isIterL3 = false;
        if(_chain_Ev->GetListOfBranches()->FindObject( "hltIterL3Muons_pt" ))
            isIterL3 = true;

        vector<TString> branch_tags = {
            "genParticle",
            "vec_my",
            "L1TkMu",
            "L2Muon",
            "hltPhase2L3OI",
            "hltPhase2L3IOFromL1",
            // "hltPhase2L3MuonsNoID",
            // "hltPhase2L3Muons",
            "L3MuonsNoId",
            "L3Muons"
            // "tpTo_hltPhase2L3OI",
            // "tpTo_hltPhase2L3IOFromL1",
            // "tpTo_hltPhase2L3MuonsNoID",
            // "tpTo_hltPhase2L3Muons"
        };

        unique_ptr<MuonHLTNtuple>  nt( new MuonHLTNtuple( _chain_Ev, branch_tags ) );

    // -- Histograms
        TH1D *h_nEvents = new TH1D("h_nEvents",  "", 3, -1, 2);

        TH1D *h_gen_pt  = new TH1D("h_gen_pt",  "", 1000, 0, 1000);
        TH1D *h_gen_eta = new TH1D("h_gen_eta", "", 60, -3, 3);
        TH1D *h_gen_phi = new TH1D("h_gen_phi", "", 64, -3.2, 3.2);

        TH1D *h_gen_acc_pt  = new TH1D("h_gen_acc_pt",  "", 1000, 0, 1000);
        TH1D *h_gen_acc_eta = new TH1D("h_gen_acc_eta", "", 60, -3, 3);
        TH1D *h_gen_acc_phi = new TH1D("h_gen_acc_phi", "", 64, -3.2, 3.2);

        TH1D *h_gen_hard_pt  = new TH1D("h_gen_hard_pt",  "", 1000, 0, 1000);
        TH1D *h_gen_hard_eta = new TH1D("h_gen_hard_eta", "", 60, -3, 3);
        TH1D *h_gen_hard_phi = new TH1D("h_gen_hard_phi", "", 64, -3.2, 3.2);

        TH1D *h_gen_hard_acc_pt  = new TH1D("h_gen_hard_acc_pt",  "", 1000, 0, 1000);
        TH1D *h_gen_hard_acc_eta = new TH1D("h_gen_hard_acc_eta", "", 60, -3, 3);
        TH1D *h_gen_hard_acc_phi = new TH1D("h_gen_hard_acc_phi", "", 64, -3.2, 3.2);

        double L1_pt_cut = 22.;
        TString thefilter = "hltL3fL1TkSingleMu22L3Filtered24Q";

        vector<TString> L3types = {
            "L1TkMuon",
            "L2Muon",
            "L3OI",
            "L3IO",
            "L3MuonNoId",
            "L3MuonNoIdInner",
            "L3Muon",
            "L3MuonInner",
            "L3Filter"
        };

        vector<vector<TH1D *>> vh_L3_qbpt_pt  = {};
        vector<vector<TH1D *>> vh_L3_qbpt_eta = {};
        vector<vector<TH1D *>> vh_L3_pt_pt    = {};
        vector<vector<TH1D *>> vh_L3_pt_eta   = {};

        for(unsigned itype=0; itype < L3types.size(); ++itype) {

            vh_L3_qbpt_pt.push_back( {} );
            vh_L3_qbpt_eta.push_back( {} );
            vh_L3_pt_pt.push_back( {} );
            vh_L3_pt_eta.push_back( {} );

            for(int ipt=0; ipt<n_pt_bins-1; ++ipt) {
                TString name = TString::Format("h_%s_qbpt_pt_%d_L1pt%.0f", L3types.at(itype).Data(), ipt, L1_pt_cut);
                TH1D *h_L3_qbpt  = new TH1D(name,  "", 4000, -2, 2);
                vh_L3_qbpt_pt.at(itype).push_back( h_L3_qbpt );

                name = TString::Format("h_%s_pt_pt_%d_L1pt%.0f", L3types.at(itype).Data(), ipt, L1_pt_cut);
                TH1D *h_L3_pt  = new TH1D(name,  "", 4000, -2, 2);
                vh_L3_pt_pt.at(itype).push_back( h_L3_pt );
            }

            for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {
                TString name = TString::Format("h_%s_qbpt_eta_%d_L1pt%.0f", L3types.at(itype).Data(), ieta, L1_pt_cut);
                TH1D *h_L3_qbpt  = new TH1D(name,  "", 4000, -2, 2);
                vh_L3_qbpt_eta.at(itype).push_back( h_L3_qbpt );

                name = TString::Format("h_%s_pt_eta_%d_L1pt%.0f", L3types.at(itype).Data(), ieta, L1_pt_cut);
                TH1D *h_L3_pt  = new TH1D(name,  "", 4000, -2, 2);
                vh_L3_pt_eta.at(itype).push_back( h_L3_pt );
            }
        }

    // -- Event loop
    for(unsigned i_ev=0; i_ev<nEvent; i_ev++) {
        if(doBar)
            loadBar(i_ev+1, nEvent, 100, 100);
        else if( doMem && i_ev !=0 && i_ev % nMem == 0 )
            printMemory("\t");
        else
            printRunTimeShort(timer_total);

        nt->GetEntry( i_ev );

        double genWeight = nt->genEventWeight > 0.0 ? 1.0 : -1.0;
        h_nEvents->Fill( genWeight );

        vector<Object> GenParticles                   = nt->get_GenParticles();

        bool isDimuon = false;
        vector<Object> GenMuonsFromHardProcess = {};
        // -- Dimuon skim for DY
            bool found0 = false;
            bool found1 = false;
            for(auto& genP: GenParticles) {

                if( fabs(genP.get("ID")) == 13 && genP.get("fromHardProcessFinalState") == 1 ) {
                    GenMuonsFromHardProcess.push_back( genP );

                    h_gen_hard_pt->Fill( genP.pt, genWeight );
                    h_gen_hard_eta->Fill( genP.eta, genWeight );
                    h_gen_hard_phi->Fill( genP.phi, genWeight );

                    if( acceptance( genP ) ) {
                        h_gen_hard_acc_pt->Fill( genP.pt, genWeight );
                        h_gen_hard_acc_eta->Fill( genP.eta, genWeight );
                        h_gen_hard_acc_phi->Fill( genP.phi, genWeight );
                    }
                }

                if( fabs(genP.get("ID")) == 13 && genP.get("status") == 1 ) {
                    h_gen_pt->Fill( genP.pt, genWeight );
                    h_gen_eta->Fill( genP.eta, genWeight );
                    h_gen_phi->Fill( genP.phi, genWeight );

                    if( acceptance( genP ) ) {
                        h_gen_acc_pt->Fill( genP.pt, genWeight );
                        h_gen_acc_eta->Fill( genP.eta, genWeight );
                        h_gen_acc_phi->Fill( genP.phi, genWeight );
                    }
                }

                if( genP.get("ID") == 13 && genP.get("isHardProcess") == 1 )
                    found0 = true;
                if( genP.get("ID") == -13 && genP.get("isHardProcess") == 1 )
                    found1 = true;
            }
            isDimuon = (found0 && found1);

        // if( (doDimuon && !isDimuon) )
        //     continue;

        // -- Get object collections
            vector<Object> HLTObjects                = nt->get_HLTObjects( thefilter );
            vector<Object> L1TkMuons                 = nt->get_L1TkMuons();
            vector<Object> L2Muons                   = nt->get_L2Muons();
            vector<Object> hltPhase2L3OI             = nt->get_hltPhase2L3OI();
            vector<Object> hltPhase2L3IOFromL1       = nt->get_hltPhase2L3IOFromL1();
            // vector<Object> hltPhase2L3MuonsNoID      = nt->get_hltPhase2L3MuonsNoID();
            // vector<Object> hltPhase2L3Muons          = nt->get_hltPhase2L3Muons();
            vector<Object> L3MuonsNoId                   = nt->get_L3MuonsNoId();
            vector<Object> L3Muons                   = nt->get_L3Muons();
            // vector<Object> tpTo_hltPhase2L3OI        = nt->get_tpTo_hltPhase2L3OI();
            // vector<Object> tpTo_hltPhase2L3IOFromL1  = nt->get_tpTo_hltPhase2L3IOFromL1();
            // vector<Object> tpTo_hltPhase2L3MuonsNoID = nt->get_tpTo_hltPhase2L3MuonsNoID();
            // vector<Object> tpTo_hltPhase2L3Muons     = nt->get_tpTo_hltPhase2L3Muons();

            // std::sort(GenParticles.begin(), GenParticles.end(), sort_by_pt());
            // std::sort(HLTObjects.begin(), HLTObjects.end(), sort_by_pt());
            // std::sort(L1TkMuons.begin(), L1TkMuons.end(), sort_by_pt());
            // std::sort(L2Muons.begin(), L2Muons.end(), sort_by_pt());
            // std::sort(hltPhase2L3OI.begin(), hltPhase2L3OI.end(), sort_by_pt());
            // std::sort(hltPhase2L3IOFromL1.begin(), hltPhase2L3IOFromL1.end(), sort_by_pt());
            // std::sort(hltPhase2L3MuonsNoID.begin(), hltPhase2L3MuonsNoID.end(), sort_by_pt());
            // std::sort(hltPhase2L3Muons.begin(), hltPhase2L3Muons.end(), sort_by_pt());
            // std::sort(L3Muons.begin(), L3Muons.end(), sort_by_pt());
            // std::sort(tpTo_hltPhase2L3OI.begin(), tpTo_hltPhase2L3OI.end(), sort_by_pt());
            // std::sort(tpTo_hltPhase2L3IOFromL1.begin(), tpTo_hltPhase2L3IOFromL1.end(), sort_by_pt());
            // std::sort(tpTo_hltPhase2L3MuonsNoID.begin(), tpTo_hltPhase2L3MuonsNoID.end(), sort_by_pt());
            // std::sort(tpTo_hltPhase2L3Muons.begin(), tpTo_hltPhase2L3Muons.end(), sort_by_pt());

            vector<vector<Object>*> L3MuonColls {
                &L1TkMuons,
                &L2Muons,
                &hltPhase2L3OI,
                &hltPhase2L3IOFromL1,
                &L3MuonsNoId,
                &L3MuonsNoId, // inner
                &L3Muons,
                &L3Muons, // inner
                &HLTObjects
            };

            // vector<vector<Object>*> tpTo_L3MuonColls {
            //     &tpTo_hltPhase2L3OI,
            //     &tpTo_hltPhase2L3IOFromL1,
            //     &tpTo_hltPhase2L3MuonsNoID,
            //     &tpTo_hltPhase2L3Muons
            // };

        // -- using Gen
            for(auto& genmu: GenMuonsFromHardProcess) {
            // for(auto& genmu: GenParticles) {

                if( fabs(genmu.get("ID")) != 13 )
                        continue;

                if( !acceptance( genmu ) )
                    continue;

                if( genmu.get("status") != 1 )
                    continue;

                if( !genmu.matched( L1TkMuons, 0.3 ) )
                    continue;

                for(unsigned itype=0; itype<L3types.size(); ++itype) {
                    vector<Object>* L3Coll = L3MuonColls.at(itype);

                    bool looseMatch = (L3types.at(itype).Contains("L1TkMuon") || L3types.at(itype).Contains("L2Muon"));

                    vector<int> L3map(L3Coll->size(), -1);
                    int matched_idx = looseMatch ? genmu.matched( *L3Coll, L3map, 0.3 ) : genmu.matched( *L3Coll, L3map, 0.1 );

                    if( matched_idx > -1 ) {

                        bool nocharge = (L3types.at(itype).Contains("L3Filter"));

                        if(L3types.at(itype).Contains("L1TkMuon") && L3Coll->at(matched_idx).has("trackCurvature")) {
                            double L1TkMuon_charge = L3Coll->at(matched_idx).get("trackCurvature") > 0. ? 1. : -1;
                            L3Coll->at(matched_idx).addVar("charge", L1TkMuon_charge);
                        }

                        double genpt   = genmu.pt;
                        double genqbpt = nocharge ? 1. / genmu.pt : genmu.get("charge") / genmu.pt;
                        double L3pt    = L3Coll->at(matched_idx).pt;
                        double L3qbpt  = nocharge ? 1. / L3Coll->at(matched_idx).pt : L3Coll->at(matched_idx).get("charge") / L3Coll->at(matched_idx).pt;

                        if(L3types.at(itype).Contains("Inner")) {
                            L3pt    = L3Coll->at(matched_idx).get("inner_pt");
                            L3qbpt  = L3Coll->at(matched_idx).get("charge") / L3Coll->at(matched_idx).get("inner_pt");
                        }

                        double res_qbpt = (genqbpt - L3qbpt) / genqbpt;
                        double res_pt   = (L3pt - genpt) / genpt;

                        for(int ipt=0; ipt<n_pt_bins-1; ++ipt) {
                            if( genmu.pt >= pt_bins[ipt] && genmu.pt < pt_bins[ipt+1] ) {
                                vh_L3_qbpt_pt.at(itype)[ipt]->Fill( res_qbpt, genWeight );
                                vh_L3_pt_pt.at(itype)[ipt]->Fill(   res_pt,   genWeight );
                                break;
                            }
                        }
                        for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {
                            if( genmu.eta >= eta_bins[ieta] && genmu.eta < eta_bins[ieta+1] && genmu.pt > 26.0 ) {
                                vh_L3_qbpt_eta.at(itype)[ieta]->Fill( res_qbpt, genWeight );
                                vh_L3_pt_eta.at(itype)[ieta]->Fill(   res_pt,   genWeight );
                                break;
                            }
                        }
                    }
                }
            }

        // -- using TP
            // for(unsigned i=0; i<L3types.size(); ++i) {
            //     vector<Object>* tpColl = tpTo_L3MuonColls.at(i);
            //     // vector<Object>* L3Coll = L3MuonColls.at(i);

            //     int imu = 0;
            //     for(auto& mu: *tpColl) {

            //         if( !mu.matched( GenMuonsFromHardProcess, 0.1 ) )
            //             continue;

            //         if( !mu.matched( L1TkMuons, 0.3 ) )
            //             continue;

            //         if( mu.get("bestMatchTrk_pt") > 0 ) {

            //             double genpt   = mu.pt;
            //             double genqbpt = mu.get("charge") / mu.pt;
            //             double L3pt    = mu.get("bestMatchTrk_pt");
            //             double L3qbpt  = mu.get("bestMatchTrk_charge") / mu.get("bestMatchTrk_pt");

            //             double res_qbpt = (genqbpt - L3qbpt) / genqbpt;
            //             double res_pt   = (L3pt - genpt) / genpt;

            //             for(int ipt=0; ipt<n_pt_bins-1; ++ipt) {
            //                 if( mu.pt >= pt_bins[ipt] && mu.pt < pt_bins[ipt+1] ) {
            //                     vh_L3_qbpt_pt[ipt]->Fill( res_qbpt, genWeight );
            //                     vh_L3_pt_pt[ipt]->Fill(   res_pt,   genWeight );
            //                     break;
            //                 }
            //             }
            //             for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {
            //                 if( mu.eta >= eta_bins[ieta] && mu.eta < eta_bins[ieta+1] && mu.pt > 26.0 ) {
            //                     vh_L3_qbpt_eta[ieta]->Fill( res_qbpt, genWeight );
            //                     vh_L3_pt_eta[ieta]->Fill(   res_pt,   genWeight );
            //                     break;
            //                 }
            //             }
            //         }

            //         imu += 1;
            //     }
            // }
    }

    // -- Save output and Clear memory
        // delete _chain_Ev;

        f_output->cd();

        h_nEvents->Write();

        h_gen_pt->Write();
        h_gen_eta->Write();
        h_gen_phi->Write();

        h_gen_acc_pt->Write();
        h_gen_acc_eta->Write();
        h_gen_acc_phi->Write();

        h_gen_hard_pt->Write();
        h_gen_hard_eta->Write();
        h_gen_hard_phi->Write();

        h_gen_hard_acc_pt->Write();
        h_gen_hard_acc_eta->Write();
        h_gen_hard_acc_phi->Write();

        for(unsigned itype=0; itype<L3types.size(); ++itype) {
            for(int ipt=0; ipt<n_pt_bins-1; ++ipt) {
                vh_L3_qbpt_pt.at(itype)[ipt]->Write();
                vh_L3_pt_pt.at(itype)[ipt]->Write();
            }
            for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {
                vh_L3_qbpt_eta.at(itype)[ieta]->Write();
                vh_L3_pt_eta.at(itype)[ieta]->Write();
            }
        }

        f_output->Close();

    delete f_output;

    printRunTime(timer_total);
}


