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
        "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2"
    };

    for(unsigned i=0; i<filters.size(); ++i) {
        if( filters.at(i) == filter )
            return true;
    }

    return false;
}


// echo 'gROOT->LoadMacro("HLTEffAnalyzer.C+"); gSystem->Exit(0);' | root -b -l
// root -l -b -q 'HLTEffAnalyzer.C("v00", "TEST")' >&log&

void HLTEffAnalyzer(
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
        TString outputDir = TString::Format("../Outputs_%s/", ver.Data());
        if (gSystem->mkdir(outputDir, kTRUE) != -1) gSystem->mkdir(outputDir,kTRUE);
        TFile *f_output = TFile::Open(fileName+"-Eff.root", "RECREATE");
        
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
        //if(_chain_Ev->GetListOfBranches()->FindObject( "hltIterL3Muons_pt" ))
        //    isIterL3 = true;

        vector<TString> branch_tags_L3FromL1TkMuon = {
            "genParticle",
            "vec_my",
            "L1TkMu",
            "L2Muon",
            "hltPhase2L3OI",
            "hltPhase2L3IOFromL1",
            "hltPhase2L3MuonsNoID",
            "hltPhase2L3Muons",
            "L3Muons", // inner
            "hltIter2IterL3FromL1MuonTrack"
        };

        vector<TString> branch_tags_IterL3 = {
            "genParticle",
            "vec_my",
            "L1Muon",
            "L2Muon",
            "iterL3OI",
            "iterL3IOFromL2",
            "iterL3FromL2",
            "iterL3IOFromL1",
            "iterL3MuonNoID",
            "iterL3Muon",
            "L3Muons" // inner
        };

        vector<TString> branch_tags = isIterL3 ? branch_tags_IterL3 : branch_tags_L3FromL1TkMuon;

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

        vector<TString> L3types_L3FromL1TkMuon = {
            "L1Muon",
            "L2Muon",
            "L3OI",
            "L3IOFromL2",
            "L3IOFromL1",
            "L3MuonNoId",
            "L3Muon",
            "L3MuonInner",
            "hltL1TkSingleMuFiltered22",
            "hltL3fL1TkSingleMu22L3Filtered50Q",
            "hltL3fL1TkSingleMu22L3Filtered24Q",
            "hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41",
            "hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40",
            "hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70",
            "hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk",
            "hltL1TkDoubleMuFiltered7",
            "hltL3fL1DoubleMu155fPreFiltered8",
            "hltDiMuon178RelTrkIsoFiltered0p4",
            "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2"
        };

        vector<TString> L3types_IterL3 = {
            "L1Muon",
            "L2Muon",
            "L3OI",
            "L3IOFromL2",
            "L3IOFromL1",
            "L3MuonNoId",
            "L3Muon",
            "L3MuonInner",
            "hltL1fL1sMu22or25L1Filtered0",
            "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"
        };

        vector<TString> L3types = isIterL3 ? L3types_IterL3 : L3types_L3FromL1TkMuon;

        unsigned L3_type_iso = 6;

        // -- Efficiency
            vector<double> Eff_genpt_mins = {
                0,
                10,
                26,
                53
            };
            vector<double> Eff_L3pt_mins = {
                7,
                8,
                22,
                24,
                50
            };

            vector<vector<HistContainer*>> hc_Eff        = {};  // Eff[L3 type][gen pt min]
            vector<vector<HistContainer*>> hc_Eff_L1Tk   = {};  // Eff[L3 type][gen pt min]

            vector<vector<HistContainer*>> hc_EffTO      = {};  // Eff[L3 type][L3 pt min]
            vector<vector<HistContainer*>> hc_EffTO_L1Tk = {};  // Eff[L3 type][L3 pt min]

            int iL3type = 0;
            for(auto& L3type: L3types) {

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

        // -- Isolations
            double iso_L3pt_min = 24.0;
            double iso_genpt_min = 26.0;

            vector<HistContainer*> hc_IsoRel = {};
            vector<HistContainer*> hc_IsoRelTO = {};

            int ivar_inner_pt = 3;
            vector<TString> variables_IsoRel = { "pt", "eta", "pu"};  // , "inner_pt" };
            vector<vector<double>> ranges_IsoRel = {
                { 1000, 0, 1000 },
                { 48, -2.4, 2.4 },
                { 250, 0, 250}
                // { 1000, 0, 1000 }
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

                                            TString name_tmp2 = name_tmp;
                                            name_tmp2 = name_tmp2.ReplaceAll("Iso_", "IsoTO_");
                                            HistContainer* hc_tmp2 = new HistContainer(
                                                name_tmp2,
                                                variables_IsoRel,
                                                ranges_IsoRel
                                            );

                                            hc_IsoRel.push_back(hc_tmp);
                                            hc_IsoRelTO.push_back(hc_tmp2);

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

        // -- Get object collections
            vector<Object> L1Muons              = nt->get_L1Muons();
            vector<Object> L1TkMuons            = nt->get_L1TkMuons();
            vector<Object> L2Muons              = nt->get_L2Muons();
            vector<Object> L3Muons              = nt->get_L3Muons();

            vector<Object> hltPhase2L3OI        = nt->get_hltPhase2L3OI();
            vector<Object> hltPhase2L3IOFromL2  = {};  // dummy
            vector<Object> hltPhase2L3IOFromL1  = nt->get_hltPhase2L3IOFromL1();
            vector<Object> hltPhase2L3MuonsNoID = nt->get_hltPhase2L3MuonsNoID();
            vector<Object> hltPhase2L3Muons     = nt->get_hltPhase2L3Muons();

            vector<Object> iterL3OI             = nt->get_iterL3OI();
            vector<Object> iterL3IOFromL2       = nt->get_iterL3IOFromL2();
            vector<Object> iterL3IOFromL1       = nt->get_iterL3IOFromL1();
            vector<Object> iterL3MuonNoID       = nt->get_iterL3MuonNoID();
            vector<Object> iterL3Muon           = nt->get_iterL3Muon();
            vector<Object> hltIter2IterL3FromL1MuonTrack = nt->get_hltIter2IterL3FromL1MuonTrack();


            vector<Object> hltL1TkSingleMuFiltered22 = nt->get_HLTObjects( "hltL1TkSingleMuFiltered22::MYHLT" );
            vector<Object> hltL3fL1TkSingleMu22L3Filtered50Q = nt->get_HLTObjects( "hltL3fL1TkSingleMu22L3Filtered50Q::MYHLT" );
            vector<Object> hltL3fL1TkSingleMu22L3Filtered24Q = nt->get_HLTObjects( "hltL3fL1TkSingleMu22L3Filtered24Q::MYHLT" );
            vector<Object> hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41 = nt->get_HLTObjects( "hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41::MYHLT" );
            vector<Object> hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40 = nt->get_HLTObjects( "hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40::MYHLT" );
            vector<Object> hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70 = nt->get_HLTObjects( "hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70::MYHLT" );
            vector<Object> hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk = nt->get_HLTObjects( "hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk::MYHLT" );

            vector<Object> hltL1TkDoubleMuFiltered7 = nt->get_HLTObjects( "hltL1TkDoubleMuFiltered7::MYHLT" );
            vector<Object> hltL1TkSingleMuFiltered15 = nt->get_HLTObjects( "hltL1TkSingleMuFiltered15::MYHLT" );
            vector<Object> hltDoubleMuon7DZ1p0 = nt->get_HLTObjects( "hltDoubleMuon7DZ1p0::MYHLT" );
            // vector<Object> hltL3fL1DoubleMu155fPreFiltered27 = nt->get_HLTObjects( "hltL3fL1DoubleMu155fPreFiltered27::MYHLT" );
            // vector<Object> hltL3fL1DoubleMu155fFiltered37 = nt->get_HLTObjects( "hltL3fL1DoubleMu155fFiltered37::MYHLT" );
            vector<Object> hltL3fL1DoubleMu155fPreFiltered8 = nt->get_HLTObjects( "hltL3fL1DoubleMu155fPreFiltered8::MYHLT" );
            // vector<Object> hltL3fL1DoubleMu155fFiltered17 = nt->get_HLTObjects( "hltL3fL1DoubleMu155fFiltered17::MYHLT" );
            vector<Object> hltDiMuon178RelTrkIsoFiltered0p4 = nt->get_HLTObjects( "hltDiMuon178RelTrkIsoFiltered0p4::MYHLT" );
            vector<Object> hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2 = nt->get_HLTObjects( "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2::MYHLT" );

            vector<Object> hltL1fL1sMu22or25L1Filtered0 = nt->get_HLTObjects( "hltL1fL1sMu22or25L1Filtered0::MYHLT" );
            vector<Object> hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q = nt->get_HLTObjects( "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q::MYHLT" );

            vector<Object> theL1Muons = isIterL3 ? L1Muons : L1TkMuons;
            vector<Object> theL1Muons_pt22 = {};
            for(auto& l1mu: theL1Muons) {
                if(l1mu.pt > 22.0) {
                    if(isIterL3) {
                        if(l1mu.get("quality") > 11) {
                            theL1Muons_pt22.push_back( l1mu.clone() );
                        }
                    }
                    else {
                        theL1Muons_pt22.push_back( l1mu.clone() );
                    }
                }
            }

            vector<vector<Object>*> L3MuonColls_L3FromL1TkMuon {
                &theL1Muons_pt22,
                &L2Muons,
                &hltPhase2L3OI,
                &hltIter2IterL3FromL1MuonTrack,
                &hltPhase2L3IOFromL2,
                &hltPhase2L3IOFromL1,
                &hltPhase2L3MuonsNoID,
                &hltPhase2L3Muons,
                &L3Muons,
                &hltL1TkSingleMuFiltered22,
                &hltL3fL1TkSingleMu22L3Filtered50Q,
                &hltL3fL1TkSingleMu22L3Filtered24Q,
                &hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41,
                &hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40,
                &hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70,
                &hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk,
                &hltL1TkDoubleMuFiltered7,
                &hltL3fL1DoubleMu155fPreFiltered8,
                &hltDiMuon178RelTrkIsoFiltered0p4,
                &hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2
            };

            vector<vector<Object>*> L3MuonColls_IterL3 {
                &theL1Muons_pt22,
                &L2Muons,
                &iterL3OI,
                &iterL3IOFromL2,
                &iterL3IOFromL1,
                &iterL3MuonNoID,
                &iterL3Muon,
                &L3Muons,
                &hltL1fL1sMu22or25L1Filtered0,
                &hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q
            };
            
            // DEBUG >> Print out objects per evt loop //
            for (auto & obj : hltIter2IterL3FromL1MuonTrack) {
                std::cout << "Evt num     : " << i_ev << std::endl;
                std::cout << "Obj content : " << std ::endl;
                obj.print();
            }

            vector<vector<Object>*> L3MuonColls = isIterL3 ? L3MuonColls_IterL3 : L3MuonColls_L3FromL1TkMuon;

        bool L1_DoubleTkMuon_15_7 = false;
        for(auto& path: *nt->vec_myFiredTrigger) {
            if( path == "L1_DoubleTkMuon_15_7" ) {
                L1_DoubleTkMuon_15_7 = true;
                break;
            }
        }

        for(unsigned i=0; i<L3types.size(); ++i) {
            vector<Object>* L3Coll = L3MuonColls.at(i);

            bool looseMatch = (L3types.at(i).Contains("L1Muon") ||
                               L3types.at(i).Contains("L2Muon") ||
                               L3types.at(i).Contains("hltL1TkSingleMuFiltered22") ||
                               L3types.at(i).Contains("hltL1TkDoubleMuFiltered7") );
            bool looseMatchIterL3 = (isIterL3 &&
                                        (L3types.at(i).Contains("L1Muon") ||
                                         L3types.at(i).Contains("hltL1fL1sMu22or25L1Filtered0") )
                                    );

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

                    // L1 seed with pT > 22 GeV
                    if( !isDoubleMuonFilter( L3types.at(i) ) ) {

                        matched_L1Muon = genmu.matched( theL1Muons_pt22, 0.3 );
                        vector<int> L3map(L3Coll->size(), -1);
                        matched_idx = looseMatch ? genmu.matched( *L3Coll, L3map, 0.3 ) : genmu.matched( *L3Coll, L3map, 0.1, 0.3 );
                        if(looseMatchIterL3) {
                            matched_idx = genmu.matched( *L3Coll, L3map, 0.5 );
                        }

                        // HERE temporary...
                        if(isIterL3) {
                            matched_L1Muon = genmu.matched( theL1Muons_pt22, 0.5 );
                            if(L3types.at(i).Contains("L1Muon")) {
                                vector<int> L3map2(L3Coll->size(), -1);
                                matched_idx = genmu.matched( *L3Coll, L3map2, 0.5 );
                            }
                        }

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
                    }

                    // HERE temporary...
                    // Double muon efficiency...
                    else {
                        matched_L1Muon = ( L1_DoubleTkMuon_15_7 &&
                                           genmu.matched( hltL1TkDoubleMuFiltered7, 0.3 ) &&
                                           genmu.matched( hltDoubleMuon7DZ1p0, 0.3 ) &&
                                           hltL1TkDoubleMuFiltered7.size() > 1 &&
                                           hltL1TkSingleMuFiltered15.size() > 0 &&
                                           hltDoubleMuon7DZ1p0.size() > 1
                                         );

                        vector<int> L3map(L3Coll->size(), -1);
                        if(L3types.at(i) == "hltL1TkDoubleMuFiltered7") {
                            matched_idx = genmu.matched( *L3Coll, L3map, 0.3 );
                        }
                        else if(L3types.at(i) == "hltL3fL1DoubleMu155fPreFiltered8") {
                            matched_idx = genmu.matched( *L3Coll, L3map, 0.1, 0.3 );
                        }
                        else if(L3types.at(i) == "hltDiMuon178RelTrkIsoFiltered0p4") {
                            matched_idx = genmu.matched( *L3Coll, L3map, 0.1, 0.3 );
                            if( matched_idx > -1 && !L3Coll->at(matched_idx).matched( hltL3fL1DoubleMu155fPreFiltered8, 0.01, 0.01 ) )
                                matched_idx = -123123;
                        }
                        else if(L3types.at(i) == "hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2") {
                            matched_idx = genmu.matched( *L3Coll, L3map, 0.1, 0.3 );
                            if( matched_idx > -1 && !L3Coll->at(matched_idx).matched( hltDiMuon178RelTrkIsoFiltered0p4, 0.01, 0.01 ) )
                                matched_idx = -123123;
                        }

                        // HERE
                        // if(matched_L1Muon && genmu.pt < 20) {
                        //     cout << endl;
                        //     cout << L3types.at(i) << "  " << hltL1TkDoubleMuFiltered7.size()
                        //                           << "  " << hltL1TkSingleMuFiltered15.size()
                        //                           << "  " << hltDoubleMuon7DZ1p0.size()
                        //                           << "  " << hltL3fL1DoubleMu155fPreFiltered8.size()
                        //                           << "  " << hltDiMuon178RelTrkIsoFiltered0p4.size()
                        //                           << "  " << hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2.size()
                        //                           << endl;
                        //     cout << "\tgenmu: " << genmu << endl;
                        //     if(matched_idx > -1)
                        //         cout << "\tL3mu:  " << L3Coll->at(matched_idx) << endl;
                        // }
                    }

                    // cout << endl;
                    // cout << L3types.at(i) << "  " << matched_L1Muon << "  " << matched_idx << endl;
                    // cout << "\t gen" << genmu << endl;
                    // if(matched_idx > -1)
                    //     cout << "\t trg" << L3Coll->at(matched_idx) << endl;

                    // HERE temporary...
                    if( matched_idx == -123123 ) {
                        cout << "\t\t not matched prev filter --> return   " << L3types.at(i) << endl;
                        return;
                    }


                    // --  Efficiency / Gen or L1Tk
                    for(unsigned j=0; j<Eff_genpt_mins.size(); ++j) {
                        if( genmu.pt > Eff_genpt_mins.at(j) ) {
                            hc_Eff.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if( matched_L1Muon ) {
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

                        if( matched_L1Muon ) {
                            hc_EffTO_L1Tk.at(i).at(j)->fill_den( genmu, nt->truePU, genWeight );

                            if( matched_idx > -1 && L3Coll->at(matched_idx).pt > Eff_L3pt_mins.at(j) ) {
                                hc_EffTO.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                                hc_EffTO_L1Tk.at(i).at(j)->fill_num( genmu, nt->truePU, genWeight );
                            }
                        }
                    }

                    // -- Iso / L3
                    if( i == L3_type_iso && matched_L1Muon && matched_idx > -1 && !isIterL3 ) {

                        Object mu  = L3Coll->at(matched_idx).clone();
                        Object mu2 = L3Muons.at(matched_idx).clone();
                        if( mu.deltaR(mu2) > 0.1 ) {
                            cout << "WARNING: dR(mu, mu2) > 0.1 !!!" << endl;
                            cout << "\t" << mu << "\t" << mu2 << "\t" << mu.deltaR(mu2) << endl;
                        }

                        double thept = mu2.get("inner_pt");

                        if( thept > iso_L3pt_min ) {

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

                                                            hc_IsoRelTO.at(itag2)->fill_den( genmu, nt->truePU, genWeight );

                                                            if(genmu.pt > iso_genpt_min)
                                                                hc_IsoRel.at(itag2)->fill_den( genmu, nt->truePU, genWeight );

                                                            if(
                                                                rel_Ecal_iso  < pfIsoEcalCut && 
                                                                rel_Hcal_iso  < pfIsoHcalCut && 
                                                                rel_Hgcal_iso < pfIsoHgcalCut && 
                                                                rel_trk_iso   < trkIsoCut
                                                            ) {
                                                                hc_IsoRelTO.at(itag2)->fill_num( genmu, nt->truePU, genWeight );

                                                                if(genmu.pt > iso_genpt_min)
                                                                    hc_IsoRel.at(itag2)->fill_num( genmu, nt->truePU, genWeight );
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
                        }
                    }
                }
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

        TDirectory* dir0 = f_output->mkdir("Eff");
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

        if(!isIterL3) {
            f_output->cd();
            TDirectory* dir2 = f_output->mkdir("Iso");
            dir2->cd();

            for(unsigned i=0; i<hc_IsoRel.size(); ++i) {
                hc_IsoRel.at(i)->Save( dir2 );
                delete hc_IsoRel.at(i);
            }

            for(unsigned i=0; i<hc_IsoRelTO.size(); ++i) {
                hc_IsoRelTO.at(i)->Save( dir2 );
                delete hc_IsoRelTO.at(i);
            }
        }

        f_output->cd();
        f_output->Close();

    delete f_output;

    printRunTime(timer_total);
}


