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


// echo 'gROOT->LoadMacro("drawBDTEffFull.C+"); gSystem->Exit(0);' | root -b -l
// rootbq 'drawBDTEffFull.C("v30", "DY PU200", "PU200-DYToLL_M50", "L1Tk")'
// rootbq 'drawBDTEffFull.C("v30", "DY PU200", "PU200-DYToLL_M50", "")'

void drawBDTEffFull(
  TString ver = "v30", TString SAMPLE = "DY PU200", TString tag = "PU200-DYToLL_M50",
  TString eff_tag = "L3Iter2FromL1", bool isLogy = false  // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString Dir = "../plots_"+ver+"/plots_BDTEffFull_"+ver+"/"+tag+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);


  TString L3_pt_min_str = "p_{T}^{HLT} > 24 GeV";
  TString gen_pt_min_str = "p_{T}^{gen} > 26 GeV";

  vector<Color_t> v_color = {
    kBlack,
    kRed,
    kBlue,
    kMagenta,

    kGreen+2,
    kYellow
  };

  vector<int> v_marker = {
    20,
    21,
    22,
    23,
    24,
    25
  };

  vector<TString> v_var = {"pt", "eta", "pu"};
  vector< vector<double> > range = {
    {1, 10, 200},  // pt
    {1, -2.4, 2.4},  // eta
    {1, 200, 201}  // PU
  };

  if(tag.Contains("PU140")) {
    range.at(2).at(1) = 140;
    range.at(2).at(2) = 141;
  }

  vector<TString> types = {
    "Eff/num_Eff_"+eff_tag+"_genpt26",
    //"Eff/num_Eff_"+eff_tag+"_genpt26",
    //"Eff/num_Eff_"+eff_tag+"_genpt26",
    //"Eff/num_Eff_"+eff_tag+"_genpt26",
    //"Eff/num_Eff_"+eff_tag+"_genpt26",
  };

  vector<TString> types_den = {
    "Eff/den_Eff_"+eff_tag+"_genpt26",
    //"Eff/den_Eff_"+eff_tag+"_genpt26",
    //"Eff/den_Eff_"+eff_tag+"_genpt26",
    //"Eff/den_Eff_"+eff_tag+"_genpt26",
    //"Eff/den_Eff_"+eff_tag+"_genpt26",
  };

  vector<TString> types_file = {
    //TString::Format("../Outputs_%s/hist-%s-%s_No-BDT.root", ver.Data(), ver.Data(), tag.Data()),
    TString::Format("../Outputs_%s/hist-%s-%s-BDT.root", ver.Data(), ver.Data(), tag.Data()),
    //TString::Format("../Outputs_%s/hist-%s-%s_sort100-BDT.root", ver.Data(), ver.Data(), tag.Data()),
    //TString::Format("../Outputs_%s/hist-%s-%s_sort50-BDT.root", ver.Data(), ver.Data(), tag.Data()),
    //TString::Format("../Outputs_%s/hist-%s-%s_sort10-BDT.root", ver.Data(), ver.Data(), tag.Data()),
    //TString::Format("../Outputs_%s/hist-%s-%s_sort0-BDT.root", ver.Data(), ver.Data(), tag.Data()),
  };

  vector<TString> types_str = {
    "Unlimited",
    //"Maximum # of seeds = 100",
    //"Maximum # of seeds = 50",
    //"Maximum # of seeds = 10",
    //"Maximum # of seeds = 0",
  };

  for(int ivar=0; ivar<(int)v_var.size(); ++ivar) {

    double xmin = range.at(ivar).at(1);
    double xmax = range.at(ivar).at(2);
    double ymin = 0.0;
    double ymax = 1.5;

    if(!v_var.at(ivar).Contains("pt")) {
      ymin = 0.8;
      ymax = 1.15;
    }

    if(eff_tag=="L3Iter2FromL1") {
      ymin = 0.0;
      ymax = 0.8;
    }

    TString canvasName = TString::Format("EffFull_%s_%s_%s", tag.Data(), eff_tag.Data(), v_var.at(ivar).Data() );
    canvasName.ReplaceAll(".","p").ReplaceAll("-","_");
    TCanvas *c;
    SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
    c->cd();
    if(isLogy)
      c->SetLogy();

    TLegend *legend;
    SetLegend( legend, 0.15, 0.70, 0.90, 0.87, -1);

    bool isFirst = true;
    for(int i = 0; i<(int)types.size(); ++i) {
    // for(int i = (int)types.size()-1; i>-1; --i) {

      TString the_type_num = types.at(i);
      TString the_type_den = types_den.at(i);
      TString the_type_str = types_str.at(i);
      TString fileName = types_file.at(i);

      TString titleX = GetTitleX(v_var.at(ivar)+"_gen");
      TString titleY = "L3 reconstruction efficiency";

      TString den_name = TString::Format("%s_%s", the_type_den.Data(), v_var.at(ivar).Data() );
      TString num_name = TString::Format("%s_%s", the_type_num.Data(), v_var.at(ivar).Data() );

      if(v_var.at(ivar) == "pt") {
        den_name = den_name.ReplaceAll("genpt26", "L3pt24");
        num_name = num_name.ReplaceAll("genpt26", "L3pt24");
      }

      TH1F* den = Get_Hist( fileName, den_name );
      TH1F* num = Get_Hist( fileName, num_name );

      if(v_var.at(ivar) == "pt") {
        den = (TH1F*)den->Rebin(n_pt_bins, den_name+"_rb", pt_bins);
        num = (TH1F*)num->Rebin(n_pt_bins, num_name+"_rb", pt_bins);
      }
      else if(v_var.at(ivar) == "eta") {
        den = (TH1F*)den->Rebin(n_eta_bins, den_name+"_rb", eta_bins);
        num = (TH1F*)num->Rebin(n_eta_bins, num_name+"_rb", eta_bins);
      }
      else{
        den = (TH1F*)den->Rebin(range.at(ivar)[0]);
        num = (TH1F*)num->Rebin(range.at(ivar)[0]);
      }

      int nbins = den->GetNbinsX();

      c->cd();

      TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins);
      g->Divide(num, den, "n e0");

      for(int ip=0; ip<nbins; ++ip) {
        if(g->GetPointY(ip) == 0.)  g->SetPointEYhigh(ip, 0.0);
      }

      double markersize = 2.0 - 0.5*i;
      if(markersize < 1.0)
        markersize = 1.0;

      g->SetTitle("");
      g->SetMarkerSize(markersize);
      g->SetMarkerStyle(v_marker.at(i));
      g->SetMarkerColor(v_color.at(i));
      g->SetLineColor(  v_color.at(i));
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

      if( v_var.at(ivar) == "pu" ) {
        int thep = tag.Contains("PU200") ? 200 : 140;
        legend->AddEntry( g, TString::Format("%s  %.5f", the_type_str.Data(), g->GetPointY(thep)), "lep" );
      }
      else
        legend->AddEntry( g, TString::Format("%s", the_type_str.Data()), "lep" );
    }

    legend->Draw();

    TLatex latex;
    Latex_Simulation_14TeV( latex );
    latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
    if(v_var.at(ivar) != "pt" )
      latex.DrawLatexNDC(0.16, 0.89, "#font[42]{#scale[0.8]{"+gen_pt_min_str+"}}");

    TString logy_tag = isLogy ? "_log" : "";
    c->Modified();  c->Update();  c->RedrawAxis();
    gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
    c->SaveAs(Dir+canvasName+logy_tag+".pdf","pdf");
    c->SaveAs(Dir+canvasName+logy_tag+".png","png");
    gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

    c->Close();
  }


  printRunTime(timer_total);
}
