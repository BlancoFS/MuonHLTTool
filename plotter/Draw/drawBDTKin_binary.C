#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TTree.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TSystem.h>

#define DEBUG (0)
#include "PlotTools.h"

using namespace std;

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

static inline void loadBar(int x, int num, int r, int w)
{
  // Only update r times.
  if( x == num )
    cout << endl;

  if ( x % (num/r +1) != 0 ) return;

  // Calculuate the ratio of complete-to-incomplete.
  float ratio = x/(float)num;
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

TH1F* getCumulative(TH1F* h)
{
  const Int_t nbinsx = h->GetNbinsX();
  TH1F* hintegrated = (TH1F*)h->Clone( (TString)h->GetName() + "_cum");
  hintegrated->Reset();

  Double_t sum   = 0.;
  Double_t sumw2 = 0.;
  for (Int_t binx = nbinsx+1; binx >= 1; --binx) {
    sum   += h->GetBinContent(binx);
    sumw2 += h->GetBinError(binx) * h->GetBinError(binx);
    hintegrated->SetBinContent(binx, sum);
    hintegrated->SetBinError(binx, sqrt(sumw2));
  }

  return hintegrated;
}

int n_pt_bins = 13-1;
double pt_bins[13] = {
  0, 10, 20, 22, 23, 24, 26, 30, 40, 60, 100, 200, 1000
};

int n_pt_bins_L122 = 25-1;
double pt_bins_L122[25] = {
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 18, 20, 21, 22, 23, 26, 30, 40, 60, 100, 200, 1000
};

int n_pt_bins_50 = 14-1;
double pt_bins_50[14] = {
  0, 40, 45, 47, 48, 49, 50, 51, 52, 55, 60, 100, 200, 1000
};

int n_eta_bins = 15-1;
double eta_bins[15] = {
  -2.4, -2.1, -1.6, -1.2, -0.9,
  -0.3, -0.2,  0.0,  0.2,  0.3,
   0.9,  1.2,  1.6,  2.1,  2.4
};

int n_eta_bins_more = 21-1;
double eta_bins_more[21] = {
  -2.4, -2.1, -1.9, -1.7, -1.6, -1.5, -1.2, -0.9,
  -0.3, -0.2,  0.0,  0.2,  0.3,
   0.9,  1.2,  1.5, 1.6, 1.7, 1.9, 2.1,  2.4
};


// echo 'gROOT->LoadMacro("drawBDTEff_binary.C+"); gSystem->Exit(0);' | root -b -l
// rootbq 'drawBDTEff_binary.C("v30", "DY PU200", "PU200-DYToLL_M50", "L1Tk")'
// rootbq 'drawBDTEff_binary.C("v30", "DY PU200", "PU200-DYToLL_M50", "")'

void drawBDTKin_binary(
  TString ver = "v00", TString SAMPLE = "DY PU200", TString tag = "PU200-DYToLL_M50",
  TString eff_tag = "L3IOFromL1", bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString Dir = "../plot_Binary/plots_"+ver+"/plots_BDTEff_"+ver+"/"+tag+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);


  TString L3_pt_min_str = "p_{T}^{HLT} > 24 GeV";
  TString gen_pt_min_str = "p_{T}^{gen} > 26 GeV";

  vector<Color_t> v_color = {
    kBlack,
    static_cast<short>(TColor::GetColor("#5790fc")),
    static_cast<short>(TColor::GetColor("#f89c20")),
    static_cast<short>(TColor::GetColor("#e42536")),
    static_cast<short>(TColor::GetColor("#964a8b")),
    static_cast<short>(TColor::GetColor("#9c9ca1")),
    static_cast<short>(TColor::GetColor("#7a21dd")),

    kMagenta,

    kGreen+2,
    kYellow+2,
  };

  vector<int> v_marker = {
    20,
    21,
    22,
    20,
    21,
    22,
    20,
    21,
    22,
    20,
    20,
    20,
    20,
    22
  };

  vector<TString> v_var = {"pt", "eta", "pu"};
  vector< vector<double> > range = {
    {1, 10, 150},  // pt
    {1, -2.4, 2.4},  // eta
    {1, 200, 201}  // PU
  };

  if(tag.Contains("PU140")) {
    range.at(2).at(1) = 140;
    range.at(2).at(2) = 141;
  }

  vector<TString> types_file = {
    // "../Analyzer/OI_default.root",
    // "../Analyzer/OI_default.root",
    // "../Analyzer/OI_default.root",

    "../Analyzer/OI_default.root",
    "../Analyzer/OI_FromL1.root",
    "../Analyzer/OI_FromL1Tk.root",
  };

  vector<TString> types = {
    // "Eff/den_Eff_L1Muon_genpt26",
    // "Eff/den_Eff_L1TkMuon_genpt26",
    // "Eff/den_Eff_L2Muon_genpt26",

    // "Eff/den_Eff_L3OI_L3pt24",
    // "Eff/den_Eff_L3OI_L3pt24",
    // "Eff/den_Eff_L3OI_L3pt24",

    "Eff/num_Eff_L3OI_L3pt24",
    "Eff/num_Eff_L3OI_L3pt24",
    "Eff/num_Eff_L3OI_L3pt24",

    // "Res/Res_L3OI_genpt26",
    // "Res/Res_L3OI_genpt26",
    // "Res/Res_L3OI_genpt26",

    // "Res/Res_sig_L3OI_genpt26",
    // "Res/Res_sig_L3OI_genpt26",
    // "Res/Resk_sig_L3OI_genpt26",

    // "Eff/num_Eff_L1Muon_genpt26",
    // "Eff/num_Eff_L1TkMuon_genpt26",
    // "Eff/num_Eff_L2Muon_genpt26",
    // "Eff/num_Eff_L3OI_genpt26",
    // "Eff/num_Eff_L3OI_genpt26",
    // "Eff/num_Eff_L3OI_genpt26",
  };

  vector<TString> types_str = {
    // "L1 Muon",
    // "L1 TkMuon",
    // "L2 Muon (From L1 TkMuon)",

    "OI FromL2 (default)",
    "OI FromL1",
    "OI FromL1Tk",

    // "Unlimited",
    // "Maximum # of seeds = 100",
    // "Maximum # of seeds = 50",
    // "Maximum # of seeds = 10",
    // "Maximum # of seeds = 0",
  };

  for(int ivar=0; ivar<v_var.size(); ++ivar) {

    double xmin = -0.5;
    double xmax = 0.5;

    TString canvasName = TString::Format("Eff_%s_%s_KIN", tag.Data(), eff_tag.Data());
    canvasName.ReplaceAll(".","p").ReplaceAll("-","_");
    TCanvas *c;
    SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
    c->cd();

    TLegend *legend;
    SetLegend( legend, 0.15, 0.70, 0.90, 0.87, -1);

    bool isFirst = true;
    for(int i = 0; i<(int)types_file.size(); ++i) {

      TString titleX = GetTitleX(v_var.at(ivar)+"_offline");
      TString titleY = "nMuons";

      TString den_name = types.at(i) + "_" + v_var.at(ivar);

      if(v_var.at(ivar) == "pt") {
        den_name = den_name.ReplaceAll("L3pt24", "L3pt0");
      }

  
      TH1F* den_clone = Get_Hist( types_file.at(i), den_name );
      // auto den_clone = (TH1F*)den->Clone(den_name+"_clone");
      // den_clone->Reset("ICES");
      
      // if(v_var.at(ivar) == "pt") {
      //   den_clone = (TH1F*)den->Rebin(n_pt_bins, den_name+"_rb", pt_bins);
      // }
      // else if(v_var.at(ivar) == "eta") {
      //   den_clone = (TH1F*)den->Rebin(n_eta_bins, den_name+"_rb", eta_bins);
      // }
      // else{
      //   den_clone = (TH1F*)den->Rebin(range.at(ivar)[0]);
      // }
      



      c->cd();

      den_clone->SetTitle("");
      den_clone->SetStats(0);
      // den_clone->SetMarkerSize(markersize);
      // den_clone->SetMarkerStyle(v_marker.at(i));
      // den_clone->SetMarkerColor(v_color.at(i));
      den_clone->SetLineColor(  v_color.at(i));
      den_clone->SetLineWidth(2);

      std::cout << " " << den_clone->GetMaximum() << " " << den_clone->GetEntries() << std::endl;

      // SetAxis_SinglePad( den->GetXaxis(), den->GetYaxis(), titleX, titleY );

      if(isFirst) {
        TH1D* frame = new TH1D("frame", "", 1, range.at(ivar)[1], range.at(ivar)[2]);
        frame->SetStats(0);
        frame->GetYaxis()->SetRangeUser(0, den_clone->GetMaximum()*1.5);
        SetAxis_SinglePad( frame->GetXaxis(), frame->GetYaxis(), titleX, titleY );
        frame->Draw("HIST");
        isFirst = false;
      }

      if (i == 1) den_clone->Scale(984000./953000.);

      den_clone->Draw("HIST same");



      legend->AddEntry( den_clone, types_str.at(i), "lep" );

    }

    legend->Draw();

    TLatex latex;
    Latex_Simulation_14TeV( latex );
    latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
    if(!(v_var.at(ivar) == "pt")) latex.DrawLatexNDC(0.16, 0.89, "#font[42]{#scale[0.8]{"+L3_pt_min_str+"}}");

    c->Modified();  c->Update();  c->RedrawAxis();
    gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
    c->SaveAs(Dir+canvasName+"_"+v_var.at(ivar)+"_genMatching.pdf","pdf");
    // c->SaveAs(Dir+canvasName+logy_tag+".png","png");
    gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

    c->Close();
  }


  printRunTime(timer_total);
}
