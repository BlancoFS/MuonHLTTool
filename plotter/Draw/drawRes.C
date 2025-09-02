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

double DSCB(double *x, double *par) {
  double N       = par[0];
  double mean    = par[1];
  double sigma   = par[2];
  double alpha_l = par[3];
  double alpha_h = par[4];
  double n_l     = par[5];
  double n_h     = par[6];

  float t = (x[0]-mean)/sigma;
  double result = -1e-9;
  double cons_l = alpha_l/n_l;
  double fact_l = (n_l/alpha_l) - alpha_l -t;
  double cons_h = alpha_h/n_h;
  double fact_h = (n_h/alpha_h) - alpha_h +t;

  if (-alpha_l <= t && alpha_h >= t) {
    result = exp(-0.5*t*t);
  }
  else if (t < -alpha_l) {
    result = exp(-0.5*alpha_l*alpha_l)*pow(cons_l*fact_l, -n_l);
  }
  else if (t > alpha_h) {
    result = exp(-0.5*alpha_h*alpha_h)*pow(cons_h*fact_h, -n_h);
  }
  return N*result;
}

// echo 'gROOT->LoadMacro("drawRes.C+"); gSystem->Exit(0);' | root -b -l

void drawRes(
  TString ver = "v30",
  TString res_tag = "pt",  // "pt"  "qbpt"
  bool isLogy = true,
  TString SAMPLE = "DY PU 200",
  TString tag = "PU200-DYToLL_M50" // HERE
) {
  TStopwatch timer_total;
  timer_total.Start();

  gStyle->SetPalette(kRainBow);
  TH1::SetDefaultSumw2(kTRUE);

  TString L3typestr = "";
  // TString L3typestr = "L3 Muon passing ID";
  double scaleL2 = 0.2;
  unsigned fitIter = 20;
  TString fitOpt = "R S Q";
  bool doGaus = true;

  vector<TString> types = {
    "L1TkMuon",
    "L2Muon",
    // "L3OI",
    // "L3IO",
    // "L3MuonNoId",
    "L3MuonNoIdInner",
    // "L3Muon"
    "L3MuonInner"
    // "L3Filter"
  };

  vector<TString> types_str = {
    "L1TkMuon",
    "L2 muon #times 0.2",
    // "Outside-in",
    // "Inside-out",
    // "L3 muon",
    // "L3 muon inner",
    // "L3 muon + ID"
    "L3 muon",
    "L3 muon + ID"
    // "L3 Filter"
  };

  vector<Color_t> v_color = {
    kMagenta,
    kBlue,
    kRed,
    kBlack,

    kGreen+2,
    kYellow
  };

  vector<int> v_marker = {
    20,
    22,
    21,
    33,
    23,
    24,
    25
  };

  vector<TString> v_var = {"mean_pt", "sigma_pt", "mean_eta", "sigma_eta"};
  vector< vector<double> > range = {
    {1, 10, 200, -0.05, 0.05},  // pt
    {1, 10, 200, 0, 0.06},  // pt
    {1, -2.4, 2.4, -0.05, 0.05},  // eta
    {1, -2.4, 2.4, 0, 0.08}  // eta
  };

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

  TString Dir = "./plots_Res_"+ver+"/";
  if (gSystem->mkdir(Dir,kTRUE) != -1)
    gSystem->mkdir(Dir,kTRUE);
  // TString canvasDir_pt = Dir+"Canvas/"+res_tag+"/pt/";
  // if (gSystem->mkdir(canvasDir_pt,kTRUE) != -1)
  //   gSystem->mkdir(canvasDir_pt,kTRUE);
  // TString canvasDir_eta = Dir+"Canvas/"+res_tag+"/eta/";
  // if (gSystem->mkdir(canvasDir_eta,kTRUE) != -1)
  //   gSystem->mkdir(canvasDir_eta,kTRUE);

  vector<vector<double>> means_pt    = { {}, {}, {}, {} };
  vector<vector<double>> means_pt_e  = { {}, {}, {}, {} };
  vector<vector<double>> sigma_pt    = { {}, {}, {}, {} };
  vector<vector<double>> sigma_pt_e  = { {}, {}, {}, {} };
  vector<vector<double>> means_eta   = { {}, {}, {}, {} };
  vector<vector<double>> means_eta_e = { {}, {}, {}, {} };
  vector<vector<double>> sigma_eta   = { {}, {}, {}, {} };
  vector<vector<double>> sigma_eta_e = { {}, {}, {}, {} };

  vector<double> pt_bins_e  = {};
  vector<double> eta_bins_e = {};

  int reb = 5;
  double xmin = -1.0;
  double xmax = 1.0;
  double ymin = 0.0;
  double ymax = 0.7;
  if(isLogy) {
    ymin = 1e-5;
    ymax = 1e2;
  }
  TString logy_tag = isLogy ? "_log" : "";

  TString titleX = "(p_{T}^{HLT} - p_{T}^{gen}) / p_{T}^{gen}";
  if(res_tag == "qbpt")  titleX = "(q^{gen}/p_{T}^{gen} - q^{HLT}/p_{T}^{HLT}) / q^{gen}/p_{T}^{gen}";

  for(int ipt=0; ipt<n_pt_bins-1; ++ipt) {
    if(pt_bins[ipt] < 10 || pt_bins[ipt] > 190)  continue;

    pt_bins_e.push_back(pt_bins[ipt]);

    TString binstr = TString::Format("%.0f < p_{T} < %.0f GeV", pt_bins[ipt], pt_bins[ipt+1]);

    TString canvasName = TString::Format("Res_%s_pt_%d_%.0f_%.0f", res_tag.Data(), ipt, pt_bins[ipt], pt_bins[ipt+1] );
    canvasName.ReplaceAll(".","p").ReplaceAll("-","_");
    TCanvas *c;
    SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
    c->cd();
    if(isLogy)
      c->SetLogy();

    TLegend *legend;
    SetLegend( legend, 0.15, 0.70, 0.90, 0.87, -1);

    for(int itype=0; itype<(int)types.size(); ++itype) {
      TString L1_pt_cut = "L1pt22";
      TString fileName = TString::Format("../Outputs_%s/hist-%s-%s-Res.root", ver.Data(), ver.Data(), tag.Data() );
      TString name = TString::Format("h_%s_%s_pt_%d_%s", types[itype].Data(), res_tag.Data(), ipt, L1_pt_cut.Data());

      TH1D* h = Get_Hist_1D( fileName, name );
      h->SetTitle("");
      h->SetStats(0);
      h->SetMarkerSize(1.0);
      h->SetMarkerStyle(v_marker[itype]);
      h->SetMarkerColor(v_color[itype]);
      h->SetLineColor(  v_color[itype]);
      h->SetLineWidth(1);

      if(types.at(itype) == "L1TkMuon") {
        if(pt_bins[ipt] > 70) {
          h->Rebin(5*reb);
        }
        else {
          h->Rebin(2*reb);
        }
      }
      else if(types.at(itype) == "L2Muon") {
        if(pt_bins[ipt] > 70) {
          h->Rebin(10*reb);
        }
        else {
          h->Rebin(5*reb);
        }
      }
      else {
        if(pt_bins[ipt] > 100) {
          h->Rebin(2*reb);
        }
        else {
          h->Rebin(reb);
        }
      }
      h->Scale(1./h->Integral(0,-1));

      h->GetXaxis()->SetRangeUser( xmin, xmax );
      h->GetYaxis()->SetRangeUser( ymin, ymax );

      SetAxis_SinglePad( h->GetXaxis(), h->GetYaxis(), titleX, "a.u." );

      double mean  = h->GetMean();
      double sigma = h->GetStdDev();

      TFitResultPtr fr = nullptr;
      // -- Gaussian
      if(doGaus) {
        TF1 *fitGaus = new TF1("fitGaus", "gaus", -0.5, 0.5);
        fr = h->Fit(fitGaus, fitOpt, "SAME", mean-2.0*sigma, mean+2.0*sigma);
        for(unsigned iter=0; iter<fitIter; ++iter) {
          fitGaus->SetParameters(
            fr->Parameter(0),
            fr->Parameter(1),
            fr->Parameter(2)
          );
          double fmin = fr->Parameter(1)-2.2*fr->Parameter(2) > -0.5 ? fr->Parameter(1)-2.2*fr->Parameter(2) : -0.5;
          double fmax = fr->Parameter(1)+2.2*fr->Parameter(2) <  0.5 ? fr->Parameter(1)+2.2*fr->Parameter(2) :  0.5;
          fr = h->Fit(fitGaus, fitOpt, "SAME", fmin, fmax);
        }
      }
      // -- DSCB
      else {
        TF1 *fitDSCB = new TF1("fitDSCB", DSCB, -0.5, 0.5, 7);
        fitDSCB->SetParNames(
          "Norm",
          "mean",
          "sigma",
          "alpha_l",
          "alpha_h",
          "n_l",
          "n_h"
        );
        fitDSCB->SetParLimits(3, 0.5, 5);
        fitDSCB->SetParLimits(4, 0.5, 5);
        fitDSCB->SetParLimits(5, 1e-1, 1e1);
        fitDSCB->SetParLimits(6, 1e-1, 1e1);
        fitDSCB->SetParameters(
          0.1,
          mean,
          sigma,
          1,
          1,
          2,
          2
        );

        fr = h->Fit(fitDSCB, fitOpt, "SAME", mean-4.0*sigma, mean+4.0*sigma);
        for(unsigned iter=0; iter<fitIter; ++iter) {
          fitDSCB->SetParameters(
            fr->Parameter(0),
            fr->Parameter(1),
            fr->Parameter(2),
            fr->Parameter(3),
            fr->Parameter(4),
            fr->Parameter(5),
            fr->Parameter(6)
          );
          double fmin = fr->Parameter(3) > 1. ? fr->Parameter(3)*fr->Parameter(2) : fr->Parameter(2);
          double fmax = fr->Parameter(4) > 1. ? fr->Parameter(4)*fr->Parameter(2) : fr->Parameter(2);
          fmin = fr->Parameter(1)-3.0*fmin > -0.5 ? fr->Parameter(1)-3.0*fmin : -0.5;
          fmax = fr->Parameter(1)+3.0*fmax <  0.5 ? fr->Parameter(1)+3.0*fmax :  0.5;
          fr = h->Fit(fitDSCB, fitOpt, "SAME", fmin, fmax);
        }
      }

      mean  = fr->Parameter(1);
      sigma = fr->Parameter(2);
      double mean_e = fr->ParError(1);
      double sigma_e = fr->ParError(2);

      cout << endl;
      cout << name << endl;
      fr->Print();

      means_pt[itype].push_back( mean );
      means_pt_e[itype].push_back( mean_e );
      sigma_pt[itype].push_back( sigma );
      sigma_pt_e[itype].push_back( sigma_e );

      double uff = h->GetBinContent(0);
      double off = h->GetBinContent(h->GetNbinsX()+1);

      // HERE
      h->GetXaxis()->SetRangeUser( -0.3, 0.3 );

      if(itype==0)  h->Draw("PE");
      else          h->Draw("PE SAME");

      TString legstr = types_str[itype];
      legstr.ReplaceAll("-L1-", "");
      legend->AddEntry( h, TString::Format("%s (mean=%.3f, sigma=%.3f, uf=%.3f, of=%.3f)", legstr.Data(), mean, sigma, uff, off), "lep" );
    }

    legend->Draw();

    TLatex latex;
    Latex_Simulation_14TeV( latex );
    latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
    latex.DrawLatexNDC(0.65, 0.87, "#font[42]{#scale[0.8]{"+binstr+"}}");

    c->Modified();  c->Update();  c->RedrawAxis();
    // gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
    // c->SaveAs(canvasDir_pt+canvasName+logy_tag+".pdf","pdf");
    // gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

    c->Close();
  }
  pt_bins_e.push_back(200.0);

  for(int ieta=0; ieta<n_eta_bins-1; ++ieta) {

    eta_bins_e.push_back(eta_bins[ieta]);

    TString binstr = TString::Format("%.1f < #eta < %.1f GeV", eta_bins[ieta], eta_bins[ieta+1]);

    TString canvasName = TString::Format("Res_%s_eta_%d_%.1f_%.1f", res_tag.Data(), ieta, eta_bins[ieta], eta_bins[ieta+1] );
    canvasName.ReplaceAll(".","p").ReplaceAll("-","_");
    TCanvas *c;
    SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
    c->cd();
    if(isLogy)
      c->SetLogy();

    TLegend *legend;
    SetLegend( legend, 0.15, 0.70, 0.90, 0.87, -1);

    for(int itype=0; itype<(int)types.size(); ++itype) {
      TString L1_pt_cut = "L1pt22";
      TString fileName = TString::Format("../Outputs_%s/hist-%s-%s-Res.root", ver.Data(), ver.Data(), tag.Data() );
      TString name = TString::Format("h_%s_%s_eta_%d_%s", types[itype].Data(), res_tag.Data(), ieta, L1_pt_cut.Data());

      TH1D* h = Get_Hist_1D( fileName, name );
      h->SetTitle("");
      h->SetStats(0);
      h->SetMarkerSize(1.0);
      h->SetMarkerStyle(v_marker[itype]);
      h->SetMarkerColor(v_color[itype]);
      h->SetLineColor(  v_color[itype]);
      h->SetLineWidth(1);


      if(types.at(itype) == "L1TkMuon") {
        h->Rebin(2*reb);
      }
      else if(types.at(itype) == "L2Muon") {
        h->Rebin(5*reb);
      }
      else {
        h->Rebin(reb);
      }
      h->Scale(1./h->Integral(0,-1));

      h->GetXaxis()->SetRangeUser( xmin, xmax );
      h->GetYaxis()->SetRangeUser( ymin, ymax );

      SetAxis_SinglePad( h->GetXaxis(), h->GetYaxis(), titleX, "a.u." );

      double mean  = h->GetMean();
      double sigma = h->GetStdDev();

      TFitResultPtr fr = nullptr;
      // -- Gaussian
      if(doGaus) {
        TF1 *fitGaus = new TF1("fitGaus", "gaus", -0.5, 0.5);
        fr = h->Fit(fitGaus, fitOpt, "SAME", mean-2.0*sigma, mean+2.0*sigma);
        for(unsigned iter=0; iter<fitIter; ++iter) {
          fitGaus->SetParameters(
            fr->Parameter(0),
            fr->Parameter(1),
            fr->Parameter(2)
          );
          double fmin = fr->Parameter(1)-2.2*fr->Parameter(2) > -0.5 ? fr->Parameter(1)-2.2*fr->Parameter(2) : -0.5;
          double fmax = fr->Parameter(1)+2.2*fr->Parameter(2) <  0.5 ? fr->Parameter(1)+2.2*fr->Parameter(2) :  0.5;
          fr = h->Fit(fitGaus, fitOpt, "SAME", fmin, fmax);
        }
      }
      // -- DSCB
      else {
        TF1 *fitDSCB = new TF1("fitDSCB", DSCB, -0.5, 0.5, 7);
        fitDSCB->SetParNames(
          "Norm",
          "mean",
          "sigma",
          "alpha_l",
          "alpha_h",
          "n_l",
          "n_h"
        );
        fitDSCB->SetParLimits(3, 0.5, 5);
        fitDSCB->SetParLimits(4, 0.5, 5);
        fitDSCB->SetParLimits(5, 1e-1, 1e1);
        fitDSCB->SetParLimits(6, 1e-1, 1e1);
        fitDSCB->SetParameters(
          0.1,
          mean,
          sigma,
          1,
          1,
          2,
          2
        );

        fr = h->Fit(fitDSCB, fitOpt, "SAME", mean-4.0*sigma, mean+4.0*sigma);
        for(unsigned iter=0; iter<fitIter; ++iter) {
          fitDSCB->SetParameters(
            fr->Parameter(0),
            fr->Parameter(1),
            fr->Parameter(2),
            fr->Parameter(3),
            fr->Parameter(4),
            fr->Parameter(5),
            fr->Parameter(6)
          );
          double fmin = fr->Parameter(3) > 1. ? fr->Parameter(3)*fr->Parameter(2) : fr->Parameter(2);
          double fmax = fr->Parameter(4) > 1. ? fr->Parameter(4)*fr->Parameter(2) : fr->Parameter(2);
          fmin = fr->Parameter(1)-3.0*fmin > -0.5 ? fr->Parameter(1)-3.0*fmin : -0.5;
          fmax = fr->Parameter(1)+3.0*fmax <  0.5 ? fr->Parameter(1)+3.0*fmax :  0.5;
          fr = h->Fit(fitDSCB, fitOpt, "SAME", fmin, fmax);
        }
      }

      mean  = fr->Parameter(1);
      sigma = fr->Parameter(2);
      double mean_e = fr->ParError(1);
      double sigma_e = fr->ParError(2);

      cout << endl;
      cout << name << endl;
      fr->Print();

      means_eta[itype].push_back( mean );
      means_eta_e[itype].push_back( mean_e );
      sigma_eta[itype].push_back( sigma );
      sigma_eta_e[itype].push_back( sigma_e );

      double uff = h->GetBinContent(0);
      double off = h->GetBinContent(h->GetNbinsX()+1);

      // HERE
      h->GetXaxis()->SetRangeUser( -0.3, 0.3 );

      if(itype==0)  h->Draw("PE");
      else          h->Draw("PE SAME");

      TString legstr = types_str[itype];
      legstr.ReplaceAll("-L1-", "");
      legend->AddEntry( h, TString::Format("%s (mean=%.3f, sigma=%.3f, uf=%.3f, of=%.3f)", legstr.Data(), mean, sigma, uff, off), "lep" );
    }

    legend->Draw();

    TLatex latex;
    Latex_Simulation_14TeV( latex );
    latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
    latex.DrawLatexNDC(0.65, 0.87, "#font[42]{#scale[0.8]{"+binstr+"}}");

    c->Modified();  c->Update();  c->RedrawAxis();
    // gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
    // c->SaveAs(canvasDir_eta+canvasName+logy_tag+".pdf","pdf");
    // gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

    c->Close();
  }
  eta_bins_e.push_back(2.4);

  vector< vector< vector<double> > > v_points = {
    means_pt,
    sigma_pt,
    means_eta,
    sigma_eta
  };

  vector< vector< vector<double> > > v_points_e = {
    means_pt_e,
    sigma_pt_e,
    means_eta_e,
    sigma_eta_e
  };


  for(int ivar=0; ivar<(int)v_var.size(); ++ivar) {

    double xmin = range[ivar][1];
    double xmax = range[ivar][2];
    double ymin = range[ivar][3];
    double ymax = range[ivar][4];

    TString titleX = v_var[ivar].Contains("pt") ? GetTitleX("pt_gen") : GetTitleX("eta_gen");
    TString titleY = res_tag.Contains("qbpt") ? "q/p_{T} residual" : "p_{T} residual";
    if( v_var[ivar].Contains("mean") )  titleY = "Mean of "+titleY;
    else                                titleY = "Sigma of "+titleY;

    vector<double> bin_e = v_var[ivar].Contains("pt") ? pt_bins_e : eta_bins_e;

    TString canvasName = TString::Format("Res_%s_%s", res_tag.Data(), v_var[ivar].Data() );
    canvasName.ReplaceAll(".","p").ReplaceAll("-","_");
    TCanvas *c;
    SetCanvas_Square( c, canvasName, kFALSE, kFALSE, 900, 900 );
    c->cd();
    c->SetLogy(false);
    if(v_var[ivar].Contains("pt"))
      c->SetLogx(true);

    TLegend *legend;
    if(v_var[ivar].Contains("pt"))
      SetLegend( legend, 0.20, 0.65, 0.90, 0.83, 0.035);
    else
      SetLegend( legend, 0.20, 0.65, 0.90, 0.83, 0.035);

    bool isFirst = true;
    for(int itype=0; itype<(int)types.size(); ++itype) {
      vector<double> points   = v_points[ivar][itype];
      vector<double> points_e = v_points_e[ivar][itype];

      double scale = (types.at(itype) == "L2Muon") ? scaleL2 : 1.0;
      int nbins = (int)points.size();
      TGraphAsymmErrors* g = new TGraphAsymmErrors(nbins);
      for(int ip=0; ip<nbins; ++ip) {
        double binc = (bin_e[ip] + bin_e[ip+1]) / 2.0;
        double binw = (bin_e[ip+1] - bin_e[ip]) / 2.0;
        g->SetPoint(ip, binc, scale*points[ip]);
        g->SetPointError(ip, binw, binw, scale*points_e[ip], scale*points_e[ip]);
      }

      g->SetTitle("");
      g->SetMarkerSize(1.5);
      g->SetMarkerStyle(v_marker[itype]);
      g->SetMarkerColor(v_color[itype]);
      g->SetLineColor(  v_color[itype]);
      g->SetLineWidth(1);

      g->GetXaxis()->SetLimits( xmin, xmax );
      g->GetXaxis()->SetRangeUser( xmin, xmax );
      g->GetYaxis()->SetRangeUser( ymin, ymax );

      SetAxis_SinglePad( g->GetXaxis(), g->GetYaxis(), titleX, titleY );

      c->cd();
      if(isFirst) {
        g->Draw("APE");
        isFirst = false;
      }
      else {
        g->Draw("PE same");
      }

      TString legstr = types_str[itype];
      legstr.ReplaceAll("-L1-", "");
      legend->AddEntry( g, TString::Format("%s", legstr.Data()), "lep" );
    }

    legend->Draw();

    TLatex latex;
    // Latex_Simulation_14TeV( latex );
    // latex.DrawLatexNDC( 0.45,0.96, "#scale[0.8]{#font[42]{"+SAMPLE+"}}");
    // latex.DrawLatexNDC(0.65, 0.87, "#font[42]{#scale[0.8]{"+L3typestr+"}}");
    if(v_var[ivar].Contains("eta"))
      latex.DrawLatexNDC(0.65, 0.87, "#font[42]{#scale[0.8]{p_{T}^{gen} > 26 GeV}}");

    if(v_var[ivar].Contains("mean"))
      continue;

    // TString logy_tag = isLogy ? "_log" : "";
    CMS_lumi(c, 98, 11);
    c->Modified();  c->Update();  c->RedrawAxis();
    gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
    c->SaveAs(Dir+canvasName+".pdf","pdf");
    c->SaveAs(Dir+canvasName+".C","C");
    c->SaveAs(Dir+canvasName+".root","root");
    gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

    c->Close();
  }

  printRunTime(timer_total);
}
