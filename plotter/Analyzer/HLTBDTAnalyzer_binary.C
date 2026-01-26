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

#include "MuonHLTNtuple_binary.h"

#define getName(var) #var

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

class ResContainer
{
public:
ResContainer(TString _Tag) {
    this->Tag = _Tag;

    this->Init();
}

void Init() {
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);
    TH1::AddDirectory(kFALSE);
    TH2::AddDirectory(kFALSE);

    v_res = new TH1D(this->Tag, "", 500, -2.5, 2.5);
}

void Fill(double hlt_pT, double gen_pT, double weight = 1.0) {
    v_res->Fill( (hlt_pT - gen_pT) / gen_pT , weight );
}

void Save( TDirectory *dir )
{
  dir->cd();

  v_res->SetDirectory(dir);
  v_res->Write();
//   delete v_res[0];
}

~ResContainer() {
    delete v_res;
}

TH1D* v_res;
TString Tag;

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

bool isDoubleMuonFilter( TString filter )
{
    vector<TString> filters = {
        "hltL1TkDoubleMuFiltered7",
        "hltL1TkSingleMuFiltered15",
        "hltDoubleMuon7DZ1p0",
        "hltL3fL1DoubleMu155fPreFiltered27",
        "hltL3fL1DoubleMu155fFiltered37",
        "hltL3fL1DoubleMu155fPreFiltered8",
        "hltL3fL1DoubleMu155fFiltered17",
        "hltDiMuon178RelTrkIsoFiltered0p4",
        "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2",
	    "hltL1TkSingleMuFiltered22",
        "hltL3fL1TkSingleMu22L3Filtered50Q",
        "hltL3fL1TkSingleMu22L3Filtered24Q",
        "hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41",
        "hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40",
        "hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70",
        "hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk"
    };

    for(unsigned i=0; i<filters.size(); ++i) {
        if( filters.at(i) == filter )
            return true;
    }

    return false;
}


// echo 'gROOT->LoadMacro("HLTBDTAnalyzer_binary.C+"); gSystem->Exit(0);' | root -b -l
// root -l -b -q 'HLTBDTAnalyzer_binary.C("v00", "TEST")' >&log&

void HLTBDTAnalyzer_binary(
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

    TFile *f_output = TFile::Open("./Phase2_Muon_" + fileName + ".root", "RECREATE");

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

    vector<TString> branch_tags = {
        "genParticle",
        "vec_my",
        "L1Muon",
        "L1TkMuon",
        "L1TkMu",
        "L2Muon",
        "hltPhase2L3OI",
        "hltPhase2L3IOFromL1",
	"hltPhase2L3MuonMerged",
        "hltPhase2L3MuonsNoID",
        "hltPhase2L3Muons",
        "L3Muons", // inner
        "hltIter0Phase2L3FromL1TkMuon",
        "hltIter2Phase2L3FromL1TkMuon",
	
    };

    unique_ptr<MuonHLTNtuple_binary>  nt( new MuonHLTNtuple_binary( _chain_Ev, branch_tags ) );

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

    vector<TString> L3types = {
        "L1Muon",
        "L2Muon",
        "L1TkMuon",
        "L3OI",
        "L3Iter0FromL1",
        "L3Iter2FromL1",
        
        "L3IOFromL1",
	"L3MuonMerged",
        "L3MuonNoId",
        "L3Muon",

        "L3MuonInner",
        "hltL1TkSingleMuFiltered22",
        "hltL3fL1TkSingleMu22L3Filtered50Q",
        "hltL3fL1TkSingleMu22L3Filtered24Q",
        "hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41",
        "hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40",
        "hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70",
        "hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk"
    };

    // -- Efficiency
    vector<double> Eff_genpt_mins = {
        0,
        26,
        53
    };
    vector<double> Eff_L3pt_mins = {
        0,
        8,
        22,
        24,
        50
    };

    vector<vector<HistContainer*>> hc_Eff        = {};  // Eff[L3 type][gen pt min]
    vector<vector<HistContainer*>> hc_Eff_L1Tk   = {};  // Eff[L3 type][gen pt min]

    vector<vector<HistContainer*>> hc_EffTO      = {};  // Eff[L3 type][L3 pt min]
    vector<vector<HistContainer*>> hc_EffTO_L1Tk = {};  // Eff[L3 type][L3 pt min]

    // std::cout << "Start" << std::endl;
    int iL3type = 0;
    for(auto& L3type: L3types) {
        // std::cout << L3type << std::endl;
        hc_Eff.push_back( {} );
        hc_Eff_L1Tk.push_back( {} );
        hc_EffTO.push_back( {} );
        hc_EffTO_L1Tk.push_back( {} );

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

        iL3type += 1;
    }

    // -- Purity
    vector<double> Purity_L3pt_mins = {
        0,
        26,
        53
    };

    vector<vector<HistContainer*>> hc_Purity_Sig           = {};
    vector<vector<HistContainer*>> hc_Purity_Sig_L1Tk      = {};
    vector<vector<HistContainer*>> hc_Purity_Sig_L3pT      = {};
    vector<vector<HistContainer*>> hc_Purity_Sig_L1Tk_L3pT = {};
    vector<vector<ResContainer*>> hc_Res                  = {};


    iL3type = 0;
    for(auto& L3type: L3types) {
        // std::cout << L3type << std::endl;
        hc_Purity_Sig.push_back( {} );
        hc_Purity_Sig_L1Tk.push_back( {} );
        hc_Purity_Sig_L3pT.push_back( {} );
        hc_Purity_Sig_L1Tk_L3pT.push_back( {} );
        hc_Res.push_back( {} );

        for(auto& Purity_L3pt_min: Purity_L3pt_mins) {
            HistContainer* hc_tmp0 = new HistContainer( TString::Format("Purity_sig_%s_genpt%.0f", L3type.Data(), Purity_L3pt_min) );
            HistContainer* hc_tmp1 = new HistContainer( TString::Format("Purity_L1Tk_sig_%s_genpt%.0f", L3type.Data(), Purity_L3pt_min) );
            HistContainer* hc_tmp2 = new HistContainer( TString::Format("Purity_sig_%s_L3pt%.0f", L3type.Data(), Purity_L3pt_min) );
            HistContainer* hc_tmp3 = new HistContainer( TString::Format("Purity_L1Tk_sig_%s_L3pt%.0f", L3type.Data(), Purity_L3pt_min) );

            ResContainer* hc_tmp4 = new ResContainer( TString::Format("Res_%s_genpt%.0f", L3type.Data(), Purity_L3pt_min) );
            
            hc_Purity_Sig.at(iL3type).push_back( hc_tmp0 );
            hc_Purity_Sig_L1Tk.at(iL3type).push_back( hc_tmp1 );
            hc_Purity_Sig_L3pT.at(iL3type).push_back( hc_tmp2 );
            hc_Purity_Sig_L1Tk_L3pT.at(iL3type).push_back( hc_tmp3 );
            hc_Res.at(iL3type).push_back( hc_tmp4 );
        }

        iL3type += 1;
    }


    std::cout << "Evt loop start" << std::endl;
    // -- Event loop
    for(unsigned i_ev=0; i_ev<nEvent; i_ev++) {
        if(doBar)
            loadBar(i_ev+1, nEvent, 100, 100);
        else if( doMem && i_ev !=0 && i_ev % nMem == 0 )
            printMemory("\t");
        else
            printRunTimeShort(timer_total);

        nt->GetEntry( i_ev );
        // std::cout << "Got entry : " << i_ev << std::endl;

        double genWeight = nt->genEventWeight > 0.0 ? 1.0 : -1.0;
        h_nEvents->Fill( genWeight );

        vector<Object> GenParticles                   = nt->get_GenParticles();
        // std::cout << "Got gen ptc" << std::endl;
        bool isDimuon = false;
        vector<Object> GenMuonsFromHardProcess = {};
        vector<Object> GenMuonsInAcc = {};
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
        // std::cout << "Is Dimuon : " << isDimuon << std::endl;

        // -- Get object collections
        vector<Object> L1Muons              = nt->get_L1Muons();
        vector<Object> L1TkMuons            = nt->get_L1TkMuons();
        vector<Object> L2Muons              = nt->get_L2Muons();
        vector<Object> L3Muons              = nt->get_L3Muons();

        vector<Object> hltPhase2L3OI          = nt->get_hltPhase2L3OI();
        vector<Object> hltPhase2L3IOFromL1    = nt->get_hltPhase2L3IOFromL1();
	vector<Object> hltPhase2L3MuonMerged = nt->get_hltPhase2L3MuonMerged();
        vector<Object> hltPhase2L3MuonsNoID   = nt->get_hltPhase2L3MuonsNoID();
        vector<Object> hltPhase2L3Muons       = nt->get_hltPhase2L3Muons();

        vector<Object> iterL3OI             = nt->get_iterL3OI();
        vector<Object> iterL3IOFromL2       = nt->get_iterL3IOFromL2();
        vector<Object> iterL3IOFromL1       = nt->get_iterL3IOFromL1();
	vector<Object> iterL3MuonMerged     = nt->get_iterL3MuonMerged();
	vector<Object> iterL3MuonNoID       = nt->get_iterL3MuonNoID();	
        vector<Object> iterL3Muon           = nt->get_iterL3Muon();

        vector<Object> hltL1TkSingleMuFiltered22 = nt->get_HLTObjects( "hltL1TkSingleMuFiltered22::MYHLT" );
        vector<Object> hltL3fL1TkSingleMu22L3Filtered50Q = nt->get_HLTObjects( "hltL3fL1TkSingleMu22L3Filtered50Q::MYHLT" );
        vector<Object> hltL3fL1TkSingleMu22L3Filtered24Q = nt->get_HLTObjects( "hltL3fL1TkSingleMu22L3Filtered24Q::MYHLT" );
        vector<Object> hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41 = nt->get_HLTObjects( "hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41::MYHLT" );
        vector<Object> hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40 = nt->get_HLTObjects( "hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40::MYHLT" );
        vector<Object> hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70 = nt->get_HLTObjects( "hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70::MYHLT" );
        vector<Object> hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk = nt->get_HLTObjects( "hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk::MYHLT" );
        vector<Object> hltIter0Phase2L3FromL1TkMuon = nt->get_hltIter0Phase2L3FromL1TkMuon();
        vector<Object> hltIter2Phase2L3FromL1TkMuon = nt->get_hltIter2Phase2L3FromL1TkMuon();

        // std::cout << "hltIter0Phase2L3FromL1TkMuon " << hltIter0Phase2L3FromL1TkMuon.size() << std::endl;
        // for (int i = 0; i < hltIter0Phase2L3FromL1TkMuon.size(); i++) {

        //     std::cout << i << " " << hltIter0Phase2L3FromL1TkMuon.at(i) << std::endl;
        // }
        
        // std::cout << "hltIter2Phase2L3FromL1TkMuon " << hltIter2Phase2L3FromL1TkMuon.size() << std::endl;
        // for (int i = 0; i < hltIter2Phase2L3FromL1TkMuon.size(); i++) {

        //     std::cout << i << " " << hltIter2Phase2L3FromL1TkMuon.at(i) << std::endl;
        // }

        vector<Object> theL1Muons_pt22 = {};
        for(auto& l1mu: L1TkMuons) {
            if(l1mu.pt > 22.0) {
                theL1Muons_pt22.push_back( l1mu.clone() );
            }
        }

        vector<vector<Object>*> L3MuonColls {
            &theL1Muons_pt22,
            &L2Muons,
            &L1TkMuons,
            &hltPhase2L3OI,
            &hltIter0Phase2L3FromL1TkMuon,
            &hltIter2Phase2L3FromL1TkMuon,
            &hltPhase2L3IOFromL1,
	    &hltPhase2L3MuonMerged,
            &hltPhase2L3MuonsNoID,
            &hltPhase2L3Muons,
	    &L3Muons,
            &hltL1TkSingleMuFiltered22,
            &hltL3fL1TkSingleMu22L3Filtered50Q,
            &hltL3fL1TkSingleMu22L3Filtered24Q,
            &hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41,
            &hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40,
            &hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70,
            &hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk
        };

        // // DEBUG >> Print out objects per evt loop //
        // for (auto & obj : hltIter2IterL3FromL1MuonTrack) {
        //     std::cout << "Evt num     : " << i_ev << std::endl;
        //     std::cout << "Obj content : " << std ::endl;
        //     obj.print();
        // }

        for(unsigned i=0; i<L3types.size(); ++i) {
            vector<Object>* L3Coll = L3MuonColls.at(i);

            bool looseMatch = (L3types.at(i).Contains("L1Muon") ||
                               L3types.at(i).Contains("L2Muon") ||
                               L3types.at(i).Contains("hltL1TkSingleMuFiltered22") ||
                               L3types.at(i).Contains("hltL1TkDoubleMuFiltered7") );

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

                    bool matched_L1Muon = false;
                    int matched_idx = -1e6;

                    matched_L1Muon = genmu.matched( theL1Muons_pt22, 0.3 );
                    vector<int> L3map(L3Coll->size(), -1);
                    matched_idx = looseMatch ? genmu.matched( *L3Coll, L3map, 0.3 ) : genmu.matched( *L3Coll, L3map, 0.1, 0.3 );

                    // HERE temporary...
                    if( L3types.at(i).Contains("L1TkSingleMu22") ) {
                        matched_L1Muon = genmu.matched( hltL1TkSingleMuFiltered22, 0.3 );
                    }

                    // HERE temporary...
                    if( matched_idx > -1 && L3types.at(i).Contains("IsoL1TkSingleMu22") ) {
                        vector<Object>* L3Coll_prev = L3MuonColls.at(i-1);
                        if( !L3Coll->at(matched_idx).matched( *L3Coll_prev, 0.01, 0.01 ) )
                            matched_idx = -123123;
                    }

                    // HERE temporary...
                    if( matched_idx == -123123 ) {
                        cout << "\t\t not matched prev filter --> return   " << L3types.at(i) << endl;
                        return;
                    }


                    // --  Efficiency / Gen or L1Tk
                    for(unsigned j=0; j<Eff_genpt_mins.size(); ++j) {
                        if( genmu.pt > Eff_genpt_mins.at(j) ) {
                            hc_Eff.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

			    if( matched_idx > -1 ) {
			      hc_Eff.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
			    }
                            if( matched_L1Muon ) {
                                hc_Eff_L1Tk.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );
                                if( matched_idx > -1 ) {
                                    hc_Eff_L1Tk.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                }
                            }
                        }
                    }

                    // --  Efficiency turn-on / Gen or L1Tk
                    for(unsigned j=0; j<Eff_L3pt_mins.size(); ++j) {
                        hc_EffTO.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );
			if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
			  hc_EffTO.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
			}
                        if( matched_L1Muon ) {
                            hc_EffTO_L1Tk.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
                                hc_EffTO_L1Tk.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                            }
                        }
                    }
                }
            }
        }

        for(unsigned i=0; i<L3types.size(); ++i) {
            vector<Object>* L3Coll = L3MuonColls.at(i);

            bool looseMatch = (L3types.at(i).Contains("L1Muon") ||
                               L3types.at(i).Contains("L1TkMuon") ||
                               L3types.at(i).Contains("L2Muon") ||
                               L3types.at(i).Contains("hltL1TkSingleMuFiltered22") ||
                               L3types.at(i).Contains("hltL1TkDoubleMuFiltered7") );

            // bool matched( const Object& other, double dR_match = 0.1, double dpt_match = 1.e9 ) {
            //   double dR  = deltaR( other );
            //   double dpt = fabs( this->pt - other.pt ) / this->pt;
            //   return ( (dR < dR_match) && (dpt < dpt_match) );
            // }

            // hc_Purity_Sig.push_back( {} ); i,j: i - L3 type, j - minpT
            // hc_Purity_Sig_L1Tk.push_back( {} );
            // hc_Purity_Bkg.push_back( {} );
            // hc_Purity_Bkg_L1Tk.push_back( {} );

            for(auto& mu: *L3Coll) {

                bool matched_L1Muon = false;
                matched_L1Muon = mu.matched( theL1Muons_pt22, 0.3 );

                for(unsigned j=0; j<Purity_L3pt_mins.size(); ++j) {
                    hc_Purity_Sig.at(i).at(j)->fill_den( mu, nt->truePU, genWeight );
                    if (mu.pt > Purity_L3pt_mins.at(j))
                        hc_Purity_Sig_L3pT.at(i).at(j)->fill_den( mu, nt->truePU, genWeight );
		    
                    if( matched_L1Muon ) {
                        hc_Purity_Sig_L1Tk.at(i).at(j)->fill_den( mu, nt->truePU, genWeight );
                        if (mu.pt > Purity_L3pt_mins.at(j))
                            hc_Purity_Sig_L1Tk_L3pT.at(i).at(j)->fill_den( mu, nt->truePU, genWeight );
                    }
                }

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

                    bool matched_idx_res = false;
                    matched_idx_res = looseMatch ? mu.matched( genmu, 0.3 ) : mu.matched( genmu, 0.3);
                    
                    for(unsigned j=0; j<Purity_L3pt_mins.size(); ++j)
                        if (matched_idx_res && mu.pt > Purity_L3pt_mins.at(j))
                            hc_Res.at(i).at(j)->Fill( mu.pt, genmu.pt, genWeight );
                    
                    if (matched_idx_res)
                        break;
                }


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

                    bool matched_bool = false;
                    matched_bool = looseMatch ? mu.matched( genmu, 0.3 ) : mu.matched( genmu, 0.1, 0.3 );

                    if (!matched_bool) 
                        continue;

                    bool gen_matched_L1 = false;
                    gen_matched_L1 = genmu.matched( theL1Muons_pt22, 0.3 );

                    for(unsigned j=0; j<Purity_L3pt_mins.size(); ++j) {

                        if (matched_bool) { 
                            hc_Purity_Sig.at(i).at(j)->fill_num( mu, nt->truePU, genWeight );
                            if (mu.pt > Purity_L3pt_mins.at(j)) {
                                hc_Purity_Sig_L3pT.at(i).at(j)->fill_num( mu, nt->truePU, genWeight );
                            }
                            
                            if( gen_matched_L1 ) {
                                hc_Purity_Sig_L1Tk.at(i).at(j)->fill_num( mu, nt->truePU, genWeight );
                                if (mu.pt > Purity_L3pt_mins.at(j)) {
                                    hc_Purity_Sig_L1Tk_L3pT.at(i).at(j)->fill_num( mu, nt->truePU, genWeight );
                                }
                            }

                        }

                    }

                    if (matched_bool)
                        break;
                }
            }
        }

        // -- Purity
        // for(auto& mu: hltIter2IterL3FromL1MuonTrack) {
        //     // double mva = (
        //     //     exp(mu.get("mva3")) / (
        //     //         exp(mu.get("mva0")) + 
        //     //         exp(mu.get("mva1")) + 
        //     //         exp(mu.get("mva2")) + 
        //     //         exp(mu.get("mva3"))
        //     //     )
        //     double mva = (
        //         1 / (1 + exp(-mu.get("mva0")))
        //     );
        //     mu.addVar("mva", mva);

        //     for(unsigned j=0; j<Purity_L3pt_mins.size(); ++j) {
        //         if( mu.pt > Purity_L3pt_mins.at(j) ) {
        //             hc_Purity_Sig.at(j)->fill_den( mu, nt->truePU, genWeight );
        //             hc_Purity_Bkg.at(j)->fill_den( mu, nt->truePU, genWeight );
        //             if( mu.get("matchedTPsize") > 0 && fabs(mu.get("bestMatchTP_pdgId")) == 13 )
        //                 hc_Purity_Sig.at(j)->fill_num( mu, nt->truePU, genWeight );
        //             else
        //                 hc_Purity_Bkg.at(j)->fill_num( mu, nt->truePU, genWeight );
        //         }
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

    TDirectory* dir0 = f_output->mkdir("Eff");
    TDirectory* dir1 = f_output->mkdir("Pur");
    TDirectory* dir2 = f_output->mkdir("Res");
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
    }

    dir1->cd();
    for(unsigned i=0; i<L3types.size(); ++i) {
        for(unsigned j=0; j<Purity_L3pt_mins.size(); ++j) {
            hc_Purity_Sig.at(i).at(j)->Save(dir1);
            hc_Purity_Sig_L1Tk.at(i).at(j)->Save(dir1);
            hc_Purity_Sig_L3pT.at(i).at(j)->Save(dir1);
            hc_Purity_Sig_L1Tk_L3pT.at(i).at(j)->Save(dir1);
            delete hc_Purity_Sig.at(i).at(j);
            delete hc_Purity_Sig_L1Tk.at(i).at(j);
            delete hc_Purity_Sig_L3pT.at(i).at(j);
            delete hc_Purity_Sig_L1Tk_L3pT.at(i).at(j);
        }
    }

    dir2->cd();
    for (unsigned i=0; i<L3types.size(); ++i) {
        for (unsigned j=0; j<Purity_L3pt_mins.size(); ++j) {
            hc_Res.at(i).at(j)->Save(dir2);
            delete hc_Res.at(i).at(j);
        }
    }


    f_output->Close();

    delete f_output;

    printRunTime(timer_total);
}


