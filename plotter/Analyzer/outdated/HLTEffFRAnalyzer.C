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

class HistContainer
{
public:
    HistContainer(
      TString _Tag,
      vector<TString> _variables = { "pt", "eta", "phi", "pu" },
      vector<vector<double>> _ranges = {
        { 1000, 0, 1000 },
        { 48, -2.4, 2.4 },
        { 60, -TMath::Pi(), TMath::Pi() },
        { 250, 0, 250}
      }
    ) {

      if(_variables.size() != _ranges.size()) {
        cout << "HistContainer: _variables.size() != _ranges.size()" << endl;
        exit(1);
      }

      this->Tag = _Tag;
      this->variables = _variables;
      this->ranges = _ranges;
      this->nVar = variables.size();

      this->Init();
    }

    void fill_den( Object obj, double PU, double weight = 1.0 ) {
      for(int k=0; k < nVar; ++k) {
        if( variables[k] == "pu" ) {
          v_den[k]->Fill( PU, weight );
        }
        else if(obj.has(variables[k])) {
          v_den[k]->Fill( obj.get(variables[k]), weight );
        }
      }
    }

    void fill_num( Object obj, double PU, double weight = 1.0 ) {
      for(int k=0; k < nVar; ++k) {
        if( variables[k] == "pu" ) {
          v_num[k]->Fill( PU, weight );
        }
        else if(obj.has(variables[k])) {
          v_num[k]->Fill( obj.get(variables[k]), weight );
        }
      }
    }

    void fill_den( TString varname, int ivar, double value, double weight = 1.0 ) {
      if( variables[ivar] == varname ) {
        v_den[ivar]->Fill( value, weight );
      }
    }

    void fill_num( TString varname, int ivar, double value, double weight = 1.0 ) {
      if( variables[ivar] == varname ) {
        v_num[ivar]->Fill( value, weight );
      }
    }

    void Save( TFile *f_output )
    {
      f_output->cd();

      for(int k=0; k < nVar; ++k) {
        v_den[k]->Write();
        v_num[k]->Write();
        delete v_den[k];
        delete v_num[k];
      }
    }

    void Save( TDirectory *dir )
    {
      dir->cd();

      for(int k=0; k < nVar; ++k) {
        v_den[k]->SetDirectory(dir);
        v_num[k]->SetDirectory(dir);
        v_den[k]->Write();
        v_num[k]->Write();
        delete v_den[k];
        delete v_num[k];
      }
    }

    ~HistContainer() {}

private:
    TString Tag;
    int nVar;
    vector<TString> variables;
    vector<vector<double>> ranges;

    vector<TH1D*> v_den;
    vector<TH1D*> v_num;

    void Init()
    {
      TH1::SetDefaultSumw2(kTRUE);
      TH2::SetDefaultSumw2(kTRUE);
      TH1::AddDirectory(kFALSE);
      TH2::AddDirectory(kFALSE);

      TString Tag_tmp = this->Tag == "" ? "" : "_"+this->Tag;

      for(int k=0; k < nVar; ++k) {
        TString name = TString::Format("%s_%s", Tag_tmp.Data(), variables[k].Data() );

        TH1D *den = new TH1D("den"+name, "", ranges[k][0], ranges[k][1], ranges[k][2]);
        TH1D *num = new TH1D("num"+name, "", ranges[k][0], ranges[k][1], ranges[k][2]);

        v_den.push_back( den );
        v_num.push_back( num );
      }
    }
};

struct sort_by_pt
{
    inline bool operator() (const Object& a, const Object& b)
    {
        return (a.pt > b.pt);
    }
};

bool acceptance( Object obj )
{
    return ( fabs(obj.eta) < 2.4 );
}

// echo 'gROOT->LoadMacro("HLTEffFRAnalyzer.C+"); gSystem->Exit(0);' | root -b -l
// root -l -b -q 'HLTEffFRAnalyzer.C("v00", "TEST")' >&log&

void HLTEffFRAnalyzer(
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
        TFile *f_output = TFile::Open(fileName+"-EffFR.root", "RECREATE");

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
            "hltPhase2L3MuonsNoID",
            "hltPhase2L3Muons",
            "L3Muons" // inner
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

        TString thefilter = "hltL3fL1TkSingleMu22L3Filtered24Q";

        vector<TString> L3types = {
            "L1TkMuon",
            "L2Muon",
            "L3OI",
            "L3IO",
            "L3MuonNoId",
            "L3Muon",
            "L3MuonInner",
            "L3Filter"
        };

        unsigned L3_type_purity_init = 2;
        unsigned L3_type_purity_fini = 5;
        unsigned L3_type_iso = 5;

        // -- Efficiency and Purity (FR)
            vector<double> Eff_genpt_mins = {
                0,
                26
            };
            vector<double> Eff_L3pt_mins = {
                24,
                50
            };
            vector<double> FR_L3pt_mins = {
                0,
                24
            };

            vector<vector<HistContainer*>> hc_Eff        = {};  // Eff[L3 type][gen pt min]
            vector<vector<HistContainer*>> hc_Eff_L1Tk   = {};  // Eff[L3 type][gen pt min]

            vector<vector<HistContainer*>> hc_EffTO      = {};  // Eff[L3 type][L3 pt min]
            vector<vector<HistContainer*>> hc_EffTO_L1Tk = {};  // Eff[L3 type][L3 pt min]

            vector<vector<HistContainer*>> hc_FR         = {};  // FR[L3 type][pt threshold]

            int iL3type = 0;
            for(auto& L3type: L3types) {

                hc_Eff.push_back( {} );
                hc_Eff_L1Tk.push_back( {} );
                hc_EffTO.push_back( {} );
                hc_EffTO_L1Tk.push_back( {} );
                hc_FR.push_back( {} );

                for(auto& Eff_genpt_min: Eff_genpt_mins) {
                    HistContainer* hc_tmp0 = new HistContainer( TString::Format("Eff_%s_genpt%.0f",      L3type.Data(), Eff_genpt_min) );
                    HistContainer* hc_tmp1 = new HistContainer( TString::Format("Eff_L1Tk_%s_genpt%.0f", L3type.Data(), Eff_genpt_min) );
                    hc_Eff.at(iL3type).push_back( hc_tmp0 );
                    hc_Eff_L1Tk.at(iL3type).push_back( hc_tmp1 );
                }

                for(auto& Eff_L3pt_min: Eff_L3pt_mins) {
                    HistContainer* hc_tmp0 = new HistContainer( TString::Format("Eff_%s_L3pt%.0f",      L3type.Data(), Eff_L3pt_min) );
                    HistContainer* hc_tmp1 = new HistContainer( TString::Format("Eff_L1Tk_%s_L3pt%.0f", L3type.Data(), Eff_L3pt_min) );
                    hc_EffTO.at(iL3type).push_back( hc_tmp0 );
                    hc_EffTO_L1Tk.at(iL3type).push_back( hc_tmp1 );
                }

                for(auto& FR_L3pt_min: FR_L3pt_mins) {
                    HistContainer* hc_tmp = new HistContainer( TString::Format("FR_%s_L3pt%.0f", L3type.Data(), FR_L3pt_min) );
                    hc_FR.at(iL3type).push_back( hc_tmp );
                }

                iL3type += 1;
            }

        // -- Isolations
            double iso_pt_min = 24.0;

            vector<HistContainer*> hc_IsoRel = {};

            int ivar_inner_pt = 3;
            vector<TString> variables_IsoRel = { "pt", "eta", "pu", "inner_pt" };
            vector<vector<double>> ranges_IsoRel = {
                { 1000, 0, 1000 },
                { 48, -2.4, 2.4 },
                { 250, 0, 250},
                { 1000, 0, 1000 }
            };

            vector<TString> pfIsoEcalTags = {
                "pfEcalIsodR0p3dRVeto0p000"
            };

            vector<vector<double>> pfIsoEcalCuts = {
                {0.32, 0.41, 1e9}
            };

            vector<TString> pfIsoHcalTags = {
                "pfHcalIsodR0p3dRVeto0p000"
            };

            vector<vector<double>> pfIsoHcalCuts = {
                {0.40, 0.50, 1e9}
            };

            vector<TString> pfIsoHgcalTags = {
                "pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00"
            };

            vector<vector<double>> pfIsoHgcalCuts = {
                {4.70, 5.52, 1e9}
            };

            vector<TString> trkIsoTags = {
                   "trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0",
                "trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0",
                       "trkIsoFulldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0",
                    "trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0"
            };

            vector<vector<double>> trkIsoCuts = {
                {0.06, 0.09, 1e9},
                {0.07, 0.11, 1e9},
                {0.08, 0.13, 1e9},
                {0.11, 0.17, 1e9}
            };

            int itag = 0;
            for(unsigned ipfIsoEcalTag=0; ipfIsoEcalTag<pfIsoEcalTags.size(); ++ipfIsoEcalTag) {
                TString pfIsoEcalTag = pfIsoEcalTags.at(ipfIsoEcalTag);

                for(unsigned ipfIsoEcalCut=0; ipfIsoEcalCut<pfIsoEcalCuts.at(ipfIsoEcalTag).size(); ++ipfIsoEcalCut) {
                    double pfIsoEcalCut = pfIsoEcalCuts.at(ipfIsoEcalTag).at(ipfIsoEcalCut);


                    for(unsigned ipfIsoHcalTag=0; ipfIsoHcalTag<pfIsoHcalTags.size(); ++ipfIsoHcalTag) {
                        TString pfIsoHcalTag = pfIsoHcalTags.at(ipfIsoHcalTag);

                        for(unsigned ipfIsoHcalCut=0; ipfIsoHcalCut<pfIsoHcalCuts.at(ipfIsoHcalTag).size(); ++ipfIsoHcalCut) {
                            double pfIsoHcalCut = pfIsoHcalCuts.at(ipfIsoHcalTag).at(ipfIsoHcalCut);


                            for(unsigned ipfIsoHgcalTag=0; ipfIsoHgcalTag<pfIsoHgcalTags.size(); ++ipfIsoHgcalTag) {
                                TString pfIsoHgcalTag = pfIsoHgcalTags.at(ipfIsoHgcalTag);

                                for(unsigned ipfIsoHgcalCut=0; ipfIsoHgcalCut<pfIsoHgcalCuts.at(ipfIsoHgcalTag).size(); ++ipfIsoHgcalCut) {
                                    double pfIsoHgcalCut = pfIsoHgcalCuts.at(ipfIsoHgcalTag).at(ipfIsoHgcalCut);


                                    for(unsigned itrkIsoTag=0; itrkIsoTag<trkIsoTags.size(); ++itrkIsoTag) {
                                        TString trkIsoTag = trkIsoTags.at(itrkIsoTag);

                                        for(unsigned itrkIsoCut=0; itrkIsoCut<trkIsoCuts.at(itrkIsoTag).size(); ++itrkIsoCut) {
                                            double trkIsoCut = trkIsoCuts.at(itrkIsoTag).at(itrkIsoCut);

                                            TString name_tmp = TString::Format(
                                                                    "Iso_%s_RelCut%.2f_%s_RelCut%.2f_%s_RelCut%.2f_%s_RelCut%.2f",
                                                                    pfIsoEcalTag.Data(),  pfIsoEcalCut,
                                                                    pfIsoHcalTag.Data(),  pfIsoHcalCut,
                                                                    pfIsoHgcalTag.Data(), pfIsoHgcalCut,
                                                                    trkIsoTag.Data(),     trkIsoCut );
                                            TString inf_tmp = TString::Format("RelCut%.2f", 1e9);
                                            name_tmp = name_tmp.ReplaceAll(inf_tmp, "RelCutInf");
                                            name_tmp = name_tmp.ReplaceAll(".", "p");
                                            HistContainer* hc_tmp = new HistContainer(
                                                name_tmp,
                                                variables_IsoRel,
                                                ranges_IsoRel
                                            );

                                            hc_IsoRel.push_back(hc_tmp);

                                            itag += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
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

        // -- Get object collections
            vector<Object> HLTObjects                = nt->get_HLTObjects( thefilter );
            vector<Object> L1TkMuons                 = nt->get_L1TkMuons();
            vector<Object> L2Muons                   = nt->get_L2Muons();
            vector<Object> hltPhase2L3OI             = nt->get_hltPhase2L3OI();
            vector<Object> hltPhase2L3IOFromL1       = nt->get_hltPhase2L3IOFromL1();
            vector<Object> hltPhase2L3MuonsNoID      = nt->get_hltPhase2L3MuonsNoID();
            vector<Object> hltPhase2L3Muons          = nt->get_hltPhase2L3Muons();
            vector<Object> L3Muons                   = nt->get_L3Muons();
            vector<Object> hltL1TkSingleMuFiltered22 = nt->get_HLTObjects( "hltL1TkSingleMuFiltered22::MYHLT" );

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
            // std::sort(tpTo_hltPhase2L3OI.begin(), tpTo_hltPhase2L3OI.end(), sort_by_pt());
            // std::sort(tpTo_hltPhase2L3IOFromL1.begin(), tpTo_hltPhase2L3IOFromL1.end(), sort_by_pt());
            // std::sort(tpTo_hltPhase2L3MuonsNoID.begin(), tpTo_hltPhase2L3MuonsNoID.end(), sort_by_pt());
            // std::sort(tpTo_hltPhase2L3Muons.begin(), tpTo_hltPhase2L3Muons.end(), sort_by_pt());

            vector<vector<Object>*> L3MuonColls {
                &L1TkMuons,
                &L2Muons,
                &hltPhase2L3OI,
                &hltPhase2L3IOFromL1,
                &hltPhase2L3MuonsNoID,
                &hltPhase2L3Muons,
                &L3Muons,
                &HLTObjects
            };

            // vector<vector<Object>*> tpTo_L3MuonColls {
            //     &tpTo_hltPhase2L3OI,
            //     &tpTo_hltPhase2L3IOFromL1,
            //     &tpTo_hltPhase2L3MuonsNoID,
            //     &tpTo_hltPhase2L3Muons
            // };

        for(unsigned i=0; i<L3types.size(); ++i) {
            // vector<Object>* tpColl = tpTo_L3MuonColls.at(i);
            vector<Object>* L3Coll = L3MuonColls.at(i);

            bool looseMatch = (L3types.at(i).Contains("L1TkMuon") || L3types.at(i).Contains("L2Muon"));

            // -- Efficiency
            if( !doDimuon || (doDimuon && isDimuon) ) {

                for(auto& genmu: GenParticles) {

                    if( fabs(genmu.get("ID")) != 13 )
                        continue;

                    if( !acceptance( genmu ) )
                        continue;

                    if( genmu.get("status") != 1 )
                        continue;

                    bool fromHardProcess = genmu.matched( GenMuonsFromHardProcess, 0.001 );
                    if( doDimuon && !fromHardProcess )
                        continue;

                    bool matched_L1TkMu = genmu.matched( L1TkMuons, 0.3 );

                    vector<int> L3map(L3Coll->size(), -1);
                    int matched_idx = looseMatch ? genmu.matched( *L3Coll, L3map, 0.3 ) : genmu.matched( *L3Coll, L3map, 0.1, 0.3 );

                    // --  Efficiency / Gen or L1Tk
                    for(unsigned j=0; j<Eff_genpt_mins.size(); ++j) {
                        if( genmu.pt > Eff_genpt_mins.at(j) ) {
                            hc_Eff.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if( matched_L1TkMu ) {
                                hc_Eff_L1Tk.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                                if( matched_idx > -1 ) {
                                    hc_Eff.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                    hc_Eff_L1Tk.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                        }
                    }

                    // --  Efficiency turn-on / Gen or L1Tk
                    for(unsigned j=0; j<Eff_L3pt_mins.size(); ++j) {
                        hc_EffTO.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                        if( matched_L1TkMu ) {
                            hc_EffTO_L1Tk.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
                                hc_EffTO.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                hc_EffTO_L1Tk.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                            }
                        }
                    }
                }
            }

            // -- using TP
                // if( !doDimuon || (doDimuon && isDimuon) ) {
                //     for(auto& mu: *tpColl) {

                //         if( fabs(mu.get("pdgId")) != 13 )
                //             continue;

                //         if( !acceptance( mu ) )
                //             continue;

                //         if( mu.get("status") != 1 )
                //             continue;

                //         bool fromHardProcess = mu.matched( GenMuonsFromHardProcess, 0.1 );
                //         if( doDimuon && !fromHardProcess )
                //             continue;

                //         bool matched_L1TkMu = mu.matched( L1TkMuons, 0.3 );

                //         // --  Efficiency / Gen or L1Tk
                //         for(unsigned j=0; j<Eff_genpt_mins.size(); ++j) {
                //             if( mu.pt > Eff_genpt_mins.at(j) ) {
                //                 hc_Eff.at(i).at(j)->fill_den( mu, nt->truePU, genWeight );

                //                 if( matched_L1TkMu ) {
                //                     hc_Eff_L1Tk.at(i).at(j)->fill_den( mu, nt->truePU, genWeight );

                //                     if( mu.get("bestMatchTrk_pt") > 0 ) {
                //                         hc_Eff.at(i).at(j)->fill_num( mu, nt->truePU, genWeight );
                //                         hc_Eff_L1Tk.at(i).at(j)->fill_num( mu, nt->truePU, genWeight );
                //                     }
                //                 }
                //             }
                //         }

                //         // --  Efficiency turn-on / Gen or L1Tk
                //         for(unsigned j=0; j<Eff_L3pt_mins.size(); ++j) {
                //             hc_EffTO.at(i).at(j)->fill_den( mu, nt->truePU, genWeight );

                //             if( matched_L1TkMu ) {
                //                 hc_EffTO_L1Tk.at(i).at(j)->fill_den( mu, nt->truePU, genWeight );

                //                 if( mu.get("bestMatchTrk_pt") > Eff_L3pt_mins.at(j) ) {
                //                     hc_EffTO.at(i).at(j)->fill_num( mu, nt->truePU, genWeight );
                //                     hc_EffTO_L1Tk.at(i).at(j)->fill_num( mu, nt->truePU, genWeight );
                //                 }
                //             }
                //         }
                //     }
                // }

            // -- Purity and Isolations
            int imu = 0;
            for(auto& mu: *L3Coll) {

                // -- Purity
                if( i >= L3_type_purity_init && i <= L3_type_purity_fini && hltL1TkSingleMuFiltered22.size() > 0 ) {
                    for(unsigned j=0; j<FR_L3pt_mins.size(); ++j) {
                        if( mu.pt > FR_L3pt_mins.at(j) ) {
                            hc_FR.at(i).at(j)->fill_den( mu, nt->truePU, genWeight );
                            if( mu.get("matchedTPsize") > 0 && fabs(mu.get("bestMatchTP_pdgId")) == 13 )
                                hc_FR.at(i).at(j)->fill_num( mu, nt->truePU, genWeight );
                        }
                    }
                }

                // -- Isolations
                if( i == L3_type_iso ) {

                    Object mu2 = L3Muons.at(imu).clone();
                    if( mu.deltaR(mu2) > 0.1 ) {
                        cout << "WARNING: dR(mu, mu2) > 0.1 !!!" << endl;
                        cout << "\t" << mu << "\t" << mu2 << "\t" << mu.deltaR(mu2) << endl;
                    }

                    double thept = mu2.get("inner_pt");

                    if(!doGenMatchForIso || mu2.matched( GenMuonsFromHardProcess, 0.1 )) {
                        if( thept > iso_pt_min ) {

                            int itag2 = 0;
                            for(unsigned ipfIsoEcalTag=0; ipfIsoEcalTag<pfIsoEcalTags.size(); ++ipfIsoEcalTag) {
                                TString pfIsoEcalTag = pfIsoEcalTags.at(ipfIsoEcalTag);

                                for(unsigned ipfIsoEcalCut=0; ipfIsoEcalCut<pfIsoEcalCuts.at(ipfIsoEcalTag).size(); ++ipfIsoEcalCut) {
                                    double pfIsoEcalCut = pfIsoEcalCuts.at(ipfIsoEcalTag).at(ipfIsoEcalCut);


                                    for(unsigned ipfIsoHcalTag=0; ipfIsoHcalTag<pfIsoHcalTags.size(); ++ipfIsoHcalTag) {
                                        TString pfIsoHcalTag = pfIsoHcalTags.at(ipfIsoHcalTag);

                                        for(unsigned ipfIsoHcalCut=0; ipfIsoHcalCut<pfIsoHcalCuts.at(ipfIsoHcalTag).size(); ++ipfIsoHcalCut) {
                                            double pfIsoHcalCut = pfIsoHcalCuts.at(ipfIsoHcalTag).at(ipfIsoHcalCut);


                                            for(unsigned ipfIsoHgcalTag=0; ipfIsoHgcalTag<pfIsoHgcalTags.size(); ++ipfIsoHgcalTag) {
                                                TString pfIsoHgcalTag = pfIsoHgcalTags.at(ipfIsoHgcalTag);

                                                for(unsigned ipfIsoHgcalCut=0; ipfIsoHgcalCut<pfIsoHgcalCuts.at(ipfIsoHgcalTag).size(); ++ipfIsoHgcalCut) {
                                                    double pfIsoHgcalCut = pfIsoHgcalCuts.at(ipfIsoHgcalTag).at(ipfIsoHgcalCut);


                                                    for(unsigned itrkIsoTag=0; itrkIsoTag<trkIsoTags.size(); ++itrkIsoTag) {
                                                        TString trkIsoTag = trkIsoTags.at(itrkIsoTag);

                                                        for(unsigned itrkIsoCut=0; itrkIsoCut<trkIsoCuts.at(itrkIsoTag).size(); ++itrkIsoCut) {
                                                            double trkIsoCut = trkIsoCuts.at(itrkIsoTag).at(itrkIsoCut);

                                                            double rel_Ecal_iso  = mu.get(pfIsoEcalTag) / thept;
                                                            double rel_Hcal_iso  = mu.get(pfIsoHcalTag) / thept;
                                                            double rel_Hgcal_iso = mu.get(pfIsoHgcalTag) / thept;
                                                            double rel_trk_iso   = mu.get(trkIsoTag) / thept;

                                                            hc_IsoRel.at(itag2)->fill_den( mu2, nt->truePU, genWeight );
                                                            // hc_IsoRel.at(itag2)->fill_den( "inner_pt", ivar_inner_pt, thept, genWeight );
                                                            if(
                                                                rel_Ecal_iso  < pfIsoEcalCut && 
                                                                rel_Hcal_iso  < pfIsoHcalCut && 
                                                                rel_Hgcal_iso < pfIsoHgcalCut && 
                                                                rel_trk_iso   < trkIsoCut
                                                            ) {
                                                                hc_IsoRel.at(itag2)->fill_num( mu2, nt->truePU, genWeight );
                                                                // hc_IsoRel.at(itag2)->fill_num( "inner_pt", ivar_inner_pt, thept, genWeight );
                                                            }

                                                            itag2 += 1;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            // -- Old Iso
                                // for(unsigned j=0; j<isoTags.size(); ++j) {
                                //     // -- Abs Isolations
                                //     for(unsigned k=0; k<cuts_IsoAbs.size(); ++k) {
                                //         double abs_iso = mu.get(isoTags.at(j));
                                //         hc_IsoAbs.at(j).at(k)->fill_den( mu, nt->truePU, genWeight );
                                //         hc_IsoAbs.at(j).at(k)->fill_den( "absiso", ivar_iso, abs_iso, genWeight );
                                //         if(abs_iso < cuts_IsoAbs.at(k)) {
                                //             hc_IsoAbs.at(j).at(k)->fill_num( mu, nt->truePU, genWeight );
                                //             hc_IsoAbs.at(j).at(k)->fill_num( "absiso", ivar_iso, abs_iso, genWeight );
                                //         }
                                //     }

                                //     // -- Rel Isolations
                                //     for(unsigned k=0; k<trkIsoCuts.size(); ++k) {
                                //         double rel_iso = mu.get(isoTags.at(j)) / mu.pt;
                                //         hc_IsoRel.at(j).at(k)->fill_den( mu, nt->truePU, genWeight );
                                //         hc_IsoRel.at(j).at(k)->fill_den( "reliso", ivar_iso, rel_iso, genWeight );
                                //         if(rel_iso < trkIsoCuts.at(k)) {
                                //             hc_IsoRel.at(j).at(k)->fill_num( mu, nt->truePU, genWeight );
                                //             hc_IsoRel.at(j).at(k)->fill_num( "reliso", ivar_iso, rel_iso, genWeight );
                                //         }
                                //     }
                                // }
                        }
                    }
                }

                imu += 1;
            }
        }
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

        TDirectory* dir0 = f_output->mkdir("EffFR");
        dir0->cd();

        for(unsigned i=0; i<L3types.size(); ++i) {
            for(unsigned j=0; j<Eff_genpt_mins.size(); ++j) {
                hc_Eff.at(i).at(j)->Save( dir0 );
                hc_Eff_L1Tk.at(i).at(j)->Save( dir0 );
                delete hc_Eff.at(i).at(j);
                delete hc_Eff_L1Tk.at(i).at(j);
            }

            for(unsigned j=0; j<Eff_L3pt_mins.size(); ++j) {
                hc_EffTO.at(i).at(j)->Save( dir0 );
                hc_EffTO_L1Tk.at(i).at(j)->Save( dir0 );
                delete hc_EffTO.at(i).at(j);
                delete hc_EffTO_L1Tk.at(i).at(j);
            }

            for(unsigned j=0; j<FR_L3pt_mins.size(); ++j) {
                hc_FR.at(i).at(j)->Save( dir0 );
                delete hc_FR.at(i).at(j);
            }
        }


        f_output->cd();
        TDirectory* dir2 = f_output->mkdir("Iso");
        dir2->cd();

        for(unsigned i=0; i<hc_IsoRel.size(); ++i) {
            hc_IsoRel.at(i)->Save( dir2 );
            delete hc_IsoRel.at(i);
        }


        // -- Old Iso
            // f_output->cd();
            // TDirectory* dir2 = f_output->mkdir("Iso");
            // dir2->cd();

            // for(unsigned j=0; j<isoTags.size(); ++j) {

            //     TDirectory* dirtmp = f_output->mkdir( TString::Format("Iso/%s", isoTags.at(j).Data()) );

            //     for(unsigned k=0; k<cuts_IsoAbs.size(); ++k) {
            //         hc_IsoAbs.at(j).at(k)->Save( f_output->GetDirectory( TString::Format("Iso/%s", isoTags.at(j).Data()) ) );
            //         delete hc_IsoAbs.at(j).at(k);
            //     }

            //     for(unsigned k=0; k<trkIsoCuts.size(); ++k) {
            //         hc_IsoRel.at(j).at(k)->Save( f_output->GetDirectory( TString::Format("Iso/%s", isoTags.at(j).Data()) ) );
            //         delete hc_IsoRel.at(j).at(k);
            //     }
            // }

        f_output->cd();
        f_output->Close();

    delete f_output;

    printRunTime(timer_total);
}


