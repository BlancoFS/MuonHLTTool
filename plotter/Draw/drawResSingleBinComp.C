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

// echo 'gROOT->LoadMacro("drawResSingleBinComp.C+"); gSystem->Exit(0);' | root -b -l
// root -l -b -q 'drawResSingleBinComp.C("qbpt")'
// root -l -b -q 'drawResSingleBinComp.C("pt")'
void drawResSingleBinComp(
  TString ver = "v12",
  TString res_tag = "pt",
  bool isLogy = false,
  TString SAMPLE = "DY PU 200",
  TString tag = "PU200-DYToLL_M50" // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  vector<TString> L3types = {
    "L1TkMuon",
    "L2Muon",
    // "L3OI",
    // "L3IO",
    // "L3MuonNoId",
    "L3MuonNoIdInner",
    // "L3Muon",
    "L3MuonInner"
    // "L3Filter"
  };

  vector<TString> L3typestr = {
    "L1TkMuon",
    "L2 muon",
    // "Outside-in",
    // "Inside-out",
    "L3 muon",
    // "L3 muon (inner)",
    "L3 muon + ID"
    // "L3 muon + ID (inner)"
  };

  vector<Color_t> v_color = {
    kMagenta,
    kBlue,
    kRed,
    // kRed+2,
    // kCyan,
    kBlack
  };

  vector<int> v_marker = {
    20,
    22,
    // 21,
    21,
    // 33,
    33
  };

  vector<int> v_line = {
    8,
    7,
    // 2,
    2,
    // 9,
    9 // 1
  };

  vector<TString> v_var = {"mean_pt", "sigma_pt", "mean_eta", "sigma_eta"};
  vector< vector<double> > range = {
    {1, 25, 200, -0.5, 0.5},  // pt
    {1, 25, 200, 0, 0.3},  // pt
    {1, -2.4, 2.4, -0.5, 0.5},  // eta
    {1, -2.4, 2.4, 0, 0.3}  // eta
  };

  const static int n_eta_bins = 15;
  double eta_bins[n_eta_bins] = {
    -2.4, -2.1, -1.6, -1.2, -0.9,
    -0.3, -0.2,  0.0,  0.2,  0.3,
     0.9,  1.2,  1.6,  2.1,  2.4
  };

  const static int n_pt_bins = 26;
  double pt_bins[n_pt_bins] = {
     0, 10, 15, 16, 17,
    18, 19, 20, 21, 22,
    23, 24, 25, 26, 27,
    30, 35, 40, 50, 60,
    100, 150, 200, 300, 500,
    1000
  };

  TString Dir = "./plots_Res_"+ver+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);

  int reb = 10;
  double xmin = -0.5;
  double xmax = 0.5;
  double ymin = 0.0;
  double ymax = 0.6;
  if(isLogy) {
    ymin = 3e-6;
    ymax = 1e3;
  }
  TString logy_tag = isLogy ? "_log" : "";

  TString titleX = "(p_{T}^{HLT} - p_{T}^{gen}) / p_{T}^{gen}";
  if(res_tag == "qbpt")  titleX = "(q^{gen}/p_{T}^{gen} - q^{HLT}/p_{T}^{HLT}) / q^{gen}/p_{T}^{gen}";

  TString canvasName = TString::Format("ResDist_%s", res_tag.Data() );
  TCanvas *c;
  SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
  c->cd();
  if(isLogy)
    c->SetLogy();

  TLegend *legend;
  SetLegend( legend, 0.20, 0.65, 0.90, 0.83, 0.035);

  const static int nL3types = (int)L3types.size();
  TH1D* h[nL3types];
  double norm = -1;
  for(int iL3type=0; iL3type<(int)L3types.size(); ++iL3type) {

    for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {
      TString L1_pt_cut = "L1pt22";
      TString fileName = TString::Format("../Outputs_%s/hist-%s-%s-Res.root", ver.Data(), ver.Data(), tag.Data() );
      TString name = TString::Format("h_%s_%s_eta_%d_%s", L3types[iL3type].Data(), res_tag.Data(), ieta, L1_pt_cut.Data());
      TH1D* h_tmp = Get_Hist_1D( fileName, name );
      if(ieta==0) h[iL3type] = (TH1D*)h_tmp->Clone(TString::Format("h_%s", tag.Data()));
      else        h[iL3type]->Add(h_tmp);

      /*
      if(L3types[iL3type] == "L3MuonNoId") {
        TString name2 = TString::Format("h_%s_%s_eta_%d_%s", "L3Muon", res_tag.Data(), ieta, L1_pt_cut.Data());
        TH1D* h_tmp2 = Get_Hist_1D( fileName, name2 );
        for(int ib=1; ib<h_tmp->GetNbinsX()+1; ++ib) {
          double x0 = h_tmp->GetBinLowEdge(ib);
          double x1 = h_tmp->GetBinLowEdge(ib+1);

          double a_0 = h_tmp->GetBinContent(ib-2);
          double a_1 = h_tmp->GetBinContent(ib-1);
          double a_2 = h_tmp->GetBinContent(ib);
          double a_3 = h_tmp->GetBinContent(ib+1);
          double a_4 = h_tmp->GetBinContent(ib+1);

          double b_0 = h_tmp2->GetBinContent(ib-2);
          double b_1 = h_tmp2->GetBinContent(ib-1);
          double b_2 = h_tmp2->GetBinContent(ib);
          double b_3 = h_tmp2->GetBinContent(ib+1);
          double b_4 = h_tmp2->GetBinContent(ib+1);

          if(a_2 < b_2) {
            cout << endl;
            cout << name << " [" << x0 << "-" << x1 << "]" << endl;
            cout << "\t L3MuonNoId: " << a_0 << ", " << a_1 << ", " << a_2 << ", " << a_3 << ", " << a_4 << endl;
            cout << "\t     L3Muon: " << b_0 << ", " << b_1 << ", " << b_2 << ", " << b_3 << ", " << b_4 << endl;
          }
        }
      }
      */
    }

    h[iL3type]->SetTitle("");
    h[iL3type]->SetStats(0);
    h[iL3type]->SetMarkerSize(1.0);
    h[iL3type]->SetMarkerStyle(v_marker[iL3type]);
    h[iL3type]->SetMarkerColor(v_color[iL3type]);
    h[iL3type]->SetLineColor(  v_color[iL3type]);
    h[iL3type]->SetLineStyle(  v_line[iL3type]);
    h[iL3type]->SetLineWidth(2);

    double int_full = h[iL3type]->Integral(0,-1);
    int ib0 = h[iL3type]->FindBin(-0.1+0.00001);
    int ib1 = h[iL3type]->FindBin( 0.1-0.00001);
    double int_core = h[iL3type]->Integral(ib0,ib1);
    cout << L3types[iL3type] << ": "
         << int_full << ", "
         << int_core << ", "
         << (int_full - int_core) << ", "
         << ((int_full - int_core)/int_full) << endl;

    if(L3types.at(iL3type) == "L2Muon") {
      h[iL3type]->Rebin(2*reb);
    }
    else
      h[iL3type]->Rebin(reb);

    // h[iL3type]->Scale(1./h[iL3type]->Integral(0,-1));
    if(iL3type == 0)
      norm = h[iL3type]->Integral(0,-1);
    h[iL3type]->Scale(1./norm);

    h[iL3type]->GetXaxis()->SetRangeUser( xmin, xmax );
    h[iL3type]->GetYaxis()->SetRangeUser( ymin, ymax );

    SetAxis_SinglePad( h[iL3type]->GetXaxis(), h[iL3type]->GetYaxis(), titleX, "a.u." );

    double mean  = h[iL3type]->GetMean();
    double sigma = h[iL3type]->GetStdDev();
    double rms   = h[iL3type]->GetRMS();

    // TFitResultPtr fr = h[iL3type]->Fit("gaus", "R S", "SAME", mean-2*sigma, mean+2*sigma);
    // mean  = fr->Parameter(1);
    // sigma = fr->Parameter(2);

    // fr = h[iL3type]->Fit("gaus", "R S", "", mean-2*sigma, mean+2*sigma);
    // mean  = fr->Parameter(1);
    // sigma = fr->Parameter(2);

    // double uff = h[iL3type]->GetBinContent(0);
    // double off = h[iL3type]->GetBinContent(h[iL3type]->GetNbinsX()+1);

    if(iL3type==0)  h[iL3type]->Draw("HIST");
    else            h[iL3type]->Draw("HIST SAME");

    legend->AddEntry( h[iL3type], TString::Format("%s #scale[0.8]{(RMS=%.3f)}", L3typestr[iL3type].Data(), rms), "l" );
    // legend->AddEntry( h[iL3type], TString::Format("%s (mean=%.3f, sigma=%.3f, uf=%.3f, of=%.3f)", L3typestr[iL3type].ReplaceAll("-L1-", "").Data(), mean, sigma, uff, off), "lep" );
  }

  legend->Draw();

  TLatex latex;
  // Latex_Simulation_14TeV( latex );
  // latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
  latex.DrawLatexNDC(0.65, 0.87, "#font[42]{#scale[0.8]{p_{T}^{gen} > 26 GeV}}");

  CMS_lumi(c, 98, 11);
  c->Modified();  c->Update();  c->RedrawAxis();
  gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
  c->SaveAs(Dir+canvasName+logy_tag+".pdf","pdf");
  c->SaveAs(Dir+canvasName+logy_tag+".C","C");
  c->SaveAs(Dir+canvasName+logy_tag+".root","root");
  gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

  c->Close();

  printRunTime(timer_total);
}
