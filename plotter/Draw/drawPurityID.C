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
#include <PlotTools.h>
#include "tdrstyle.C"
#include "CMS_lumi.C"

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

// echo 'gROOT->LoadMacro("drawPurityID.C+"); gSystem->Exit(0);' | root -b -l
// rootbq 'drawPurityID.C("v30", "DY PU200", "DYToLL_M50", "L1Tk")'

void drawPurityID(
  TString ver = "v30", TString SAMPLE = "DY PU200", TString tag = "PU200-DYToLL_M50",
  TString eff_tag = "", bool isLogy = false
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  //-- Input file

  vector<TString> pt_ranges = {
    "L3pt0",
    "L3pt24"
  };

  vector<TString> pt_ranges_s = {
    "p_{T}^{HLT} > 0 GeV",
    "p_{T}^{HLT} > 24 GeV"
  };

  TString fileName = TString::Format("../Outputs_%s/hist-%s-%s-EffFR.root", ver.Data(), ver.Data(), tag.Data());

  vector<Color_t> v_color = {
    // kGray+2,
    // kBlack,
    // kBlue,
    kRed,
    // kMagenta,
    kGreen+2,
    kYellow
  };

  vector<int> v_marker = {
    // 26,
    // 20,
    // 22,
    23,
    21,
    // 24,
    // 25
  };

  vector<TString> types = {
    // "FR_L3OI",
    // "FR_L3IO",
    "FR_L3MuonNoId",
    "FR_L3Muon"
  };

  vector<TString> types_str = {
    // "Gen muon / L3 OI",
    // "Gen muon / L3 IO",
    // "Gen muon / L3 Muon",
    // "Gen muon / L3 Muon passing ID"
    "L3 muons",
    "L3 muons + muon ID"
  };

  vector<TString> v_var = {"pt", "eta"};  // , "pu"};
  vector< vector<double> > range = {
    {1, 5, 150},  // pt
    {1, -2.4, 2.4},  // eta
    {1, 200, 201}  // PU
  };

  int n_pt_bins = 28-1;
  double pt_bins[28] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 27, 30, 35, 40, 50, 60, 90, 150, 250, 1000
  };

  int n_eta_bins = 15-1;
  double eta_bins[15] = {
    -2.4, -2.1, -1.6, -1.2, -0.9,
    -0.3, -0.2,  0.0,  0.2,  0.3,
     0.9,  1.2,  1.6,  2.1,  2.4
  };

  TString Dir = "./plots_PurityID_"+ver+"/"+tag+"/";
    if (gSystem->mkdir(Dir,kTRUE) != -1)
      gSystem->mkdir(Dir,kTRUE);

  for(int ivar=0; ivar<(int)v_var.size(); ++ivar) {

    double xmin = range[ivar][1];
    double xmax = range[ivar][2];
    double ymin = 0.0;
    double ymax = 1.5;
    if(eff_tag.Contains("L1Tk")) {
      ymin = 0.8;
      ymax = 1.15;
    }

    for(int ipt=0; ipt<(int)pt_ranges.size(); ++ipt) {
      if(v_var[ivar]=="pt" && pt_ranges[ipt]=="L3pt24")
        continue;
      if(v_var[ivar]!="pt" && pt_ranges[ipt]=="L3pt0")
        continue;

      if(pt_ranges.at(ipt) == "L3pt24") {
        ymin = 0.5;
        ymax = 1.3;
      }

      TString canvasName = TString::Format("PurityID%s_%s_%s_%s", eff_tag.Data(), tag.Data(), v_var[ivar].Data(), pt_ranges[ipt].Data() );
      canvasName.ReplaceAll(".","p").ReplaceAll("-","_");
      TCanvas *c;
      SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
      c->cd();
      if(isLogy)
        c->SetLogy();

      TLegend *legend;
      SetLegend( legend, 0.15, 0.69, 0.80, 0.84, 0.035);

      bool isFirst = true;
      for(int i = 0; i<(int)types.size(); ++i) {

        TString the_type = types[i];
        TString the_type_str = types_str[i];
        // if(eff_tag.Contains("L1Tk")) {
        //   the_type.ReplaceAll("Eff_", "Eff_L1Tk_");
        //   the_type_str.ReplaceAll("Gen muon", "L1TkMuon");
        // }

        TString titleX = the_type.Contains("Eff_") ? GetTitleX(v_var[ivar]+"_gen") : GetTitleX(v_var[ivar]+"_gen").ReplaceAll("{gen}", "{HLT}");
        TString titleY = the_type.Contains("Eff_") ? "Efficiency" : "Purity";

        TString den_name = TString::Format("EffFR/den_%s_%s_%s", the_type.Data(), pt_ranges[ipt].Data(), v_var[ivar].Data() );
        TString num_name = TString::Format("EffFR/num_%s_%s_%s", the_type.Data(), pt_ranges[ipt].Data(), v_var[ivar].Data() );

        TH1F* den = Get_Hist( fileName, den_name );
        TH1F* num = Get_Hist( fileName, num_name );

        if(v_var[ivar] == "pt") {
          den = (TH1F*)den->Rebin(n_pt_bins, den_name+"_rb", pt_bins);
          num = (TH1F*)num->Rebin(n_pt_bins, num_name+"_rb", pt_bins);
        }
        else if(v_var[ivar] == "eta") {
          den = (TH1F*)den->Rebin(n_eta_bins, den_name+"_rb", eta_bins);
          num = (TH1F*)num->Rebin(n_eta_bins, num_name+"_rb", eta_bins);
        }
        else{
          den = (TH1F*)den->Rebin(range[ivar][0]);
          num = (TH1F*)num->Rebin(range[ivar][0]);
        }

        int nbins = den->GetNbinsX();

        c->cd();

        TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins);
        g->Divide(num, den, "n e0");

        for(int ip=0; ip<nbins; ++ip) {
          if(g->GetPointY(ip) == 0.)  g->SetPointEYhigh(ip, 0.0);
        }

        g->SetTitle("");
        g->SetMarkerSize(1.5);
        g->SetMarkerStyle(v_marker[i]);
        g->SetMarkerColor(v_color[i]);
        g->SetLineColor(  v_color[i]);
        g->SetLineWidth(1);

        g->GetXaxis()->SetLimits( xmin, xmax );
        g->GetXaxis()->SetRangeUser( xmin, xmax );
        g->GetYaxis()->SetRangeUser( ymin, ymax );

        SetAxis_SinglePad( g->GetXaxis(), g->GetYaxis(), titleX, titleY );

        if(isFirst) {
          g->Draw("APE");
          isFirst = false;
        }
        else {
          g->Draw("PE same");
        }

        if( v_var[ivar] == "pu" ) {
          legend->AddEntry( g, TString::Format("%s  %.3f", the_type_str.Data(), g->GetPointY(200)), "lep" );
        }
        else
          legend->AddEntry( g, TString::Format("%s", the_type_str.Data()), "lep" );
      }

      legend->Draw();

      TLatex latex;
      // Latex_Simulation_14TeV( latex );
      // latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
      if(v_var[ivar] == "eta" || pt_ranges[ipt] != "L3pt0" )
        latex.DrawLatexNDC(0.65, 0.87, "#font[42]{#scale[0.8]{"+pt_ranges_s[ipt]+"}}");

      TString logy_tag = isLogy ? "_log" : "";
      CMS_lumi(c, 98, 11);
      c->Modified();  c->Update();  c->RedrawAxis();
      gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
      c->SaveAs(Dir+canvasName+logy_tag+".pdf","pdf");
      c->SaveAs(Dir+canvasName+logy_tag+".C","C");
      c->SaveAs(Dir+canvasName+logy_tag+".root","root");
      gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

      c->Close();
    }
  }

  printRunTime(timer_total);
}
