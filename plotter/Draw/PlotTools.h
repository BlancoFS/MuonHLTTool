#pragma once

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TProfile.h>
#include <TColor.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLatex.h>
#include <TFile.h>
#include <TPad.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include "Math/SpecFunc.h"
#include "Math/DistFunc.h"

#include <vector>

using namespace std;


// -- Basic -- //
  TString makeDir( TString path = "./plots/")
  {
    if(gSystem->mkdir(path,kTRUE) != -1)
      gSystem->mkdir(path,kTRUE);

    return path;
  }

  // http://igotit.tistory.com/entry/C-함수-인자로-포인터-전달하고-함수내에서-동적-메모리-할당-받기-2가지-방식
  TH1F* Get_Hist(TString FileName, TString HistName, TString HistName_New = "" )
  { 
    std::cout << "Get_Hist : FileName : " << FileName << std::endl;
    std::cout << "Get_Hist : HistName : " << HistName << std::endl;
    TH1::AddDirectory(kFALSE);

    TFile *f_input = TFile::Open( FileName );
    TH1F* h_temp = (TH1F*)f_input->Get(HistName)->Clone();
    if( HistName_New != "" )
      h_temp->SetName( HistName_New );

    f_input->Close();
    // delete f_input;

    return h_temp;
  }

  TH1D* Get_Hist_1D(TString FileName, TString HistName, TString HistName_New = "" )
  {
    TH1::AddDirectory(kFALSE);

    TFile *f_input = TFile::Open( FileName );
    TH1D* h_temp = (TH1D*)f_input->Get(HistName)->Clone();
    if( HistName_New != "" )
      h_temp->SetName( HistName_New );

    f_input->Close();
    // delete f_input;

    return h_temp;
  }

  TH2F* Get_Hist_2D(TString FileName, TString HistName, TString HistName_New = "" )
  {
    TH2::AddDirectory(kFALSE);

    TFile *f_input = TFile::Open( FileName );
    TH2F* h_temp = (TH2F*)f_input->Get(HistName)->Clone();
    if( HistName_New != "" )
      h_temp->SetName( HistName_New );

    f_input->Close();
    // delete f_input;

    return h_temp;
  }

  TH3D* Get_Hist_3D(TString FileName, TString HistName, TString HistName_New = "" )
  {
    TH3::AddDirectory(kFALSE);

    TFile *f_input = TFile::Open( FileName );
    TH3D* h_temp = (TH3D*)f_input->Get(HistName)->Clone();
    if( HistName_New != "" )
      h_temp->SetName( HistName_New );

    f_input->Close();
    // delete f_input;

    return h_temp;
  }

  TProfile* Get_Profile(TString FileName, TString HistName, TString HistName_New = "" )
  {
    TH1::AddDirectory(kFALSE);

    TFile *f_input = TFile::Open( FileName );
    TProfile* h_temp = (TProfile*)f_input->Get(HistName)->Clone();
    if( HistName_New != "" )
      h_temp->SetName( HistName_New );

    f_input->Close();

    return h_temp;
  }

  TH1F *Get_Projection( TH2F* _h2, TString axisProj,
                        Double_t minP = -1e99, Double_t maxP = 1e99, Bool_t isAbs = kFALSE )
  {
    TString name = _h2->GetName();
    TH2F *h2 = (TH2F*)_h2->Clone();
    TH1F *h;
    TH1F *hm;

    TH1F *h_p;

    Int_t minB = 0;
    Int_t maxB = -1;
    Int_t minBm = 0;
    Int_t maxBm = -1;

    if(axisProj=="X" || axisProj=="x") {
      h_p = (TH1F*)h2->ProjectionY();
      for(Int_t i=0; i<=h_p->GetNbinsX()+1; ++i) {
        if( h_p->GetBinCenter(i-1) < minP && h_p->GetBinCenter(i) > minP )
          minB = i;
        if( h_p->GetBinCenter(i) < maxP && h_p->GetBinCenter(i+1) > maxP )
          maxB = i;
        if( h_p->GetBinCenter(i-1) < -1*maxP && h_p->GetBinCenter(i) > -1*maxP )
          minBm = i;
        if( h_p->GetBinCenter(i) < -1*minP && h_p->GetBinCenter(i+1) > -1*minP )
          maxBm = i;
      }

      h = (TH1F*)h2->ProjectionX(name+"_px_",minB,maxB,"e");
      if(isAbs) {
        hm = (TH1F*)h2->ProjectionX(name+"_px__",minBm,maxBm,"e");
        h->Add(hm);
      }
    }

    else if(axisProj=="Y" || axisProj=="y") {
      h_p = (TH1F*)h2->ProjectionX();
      for(Int_t i=0; i<=h_p->GetNbinsX()+1; ++i) {
        if( h_p->GetBinCenter(i-1) < minP && h_p->GetBinCenter(i) > minP )
          minB = i;
        if( h_p->GetBinCenter(i) < maxP && h_p->GetBinCenter(i+1) > maxP )
          maxB = i;
        if( h_p->GetBinCenter(i-1) < -1*maxP && h_p->GetBinCenter(i) > -1*maxP )
          minBm = i;
        if( h_p->GetBinCenter(i) < -1*minP && h_p->GetBinCenter(i+1) > -1*minP )
          maxBm = i;
      }

      h = (TH1F*)h2->ProjectionY(name+"_py_",minB,maxB,"e");
      if(isAbs) {
        hm = (TH1F*)h2->ProjectionY(name+"_py__",minBm,maxBm,"e");
        h->Add(hm);
      }
    }

    else {
      h = new TH1F("failed", "", 1, 0, 1);
      cout << "WARNING : Get_Projection : Wrong projection axis" << endl;
    }

    //cout << minB << "\t" << maxB << "\t" << minBm << "\t" << maxBm << endl;

    return h;
  }

  TGraphAsymmErrors* Get_Graph(TString FileName, TString GraphName, TString GraphName_New = "" )
  {
    TFile *f_input = TFile::Open( FileName );
    TGraphAsymmErrors* g_temp = (TGraphAsymmErrors*)f_input->Get(GraphName)->Clone();
    if( GraphName_New != "" )
      g_temp->SetName( GraphName_New );

    f_input->Close();

    return g_temp;
  }

  void CopyVectorToArray( std::vector<double>& vec, double*& arr ) {
    int nBin = (int)vec.size()-1; // -- # bins = # bin edges-1 -- //
    arr = new Double_t[nBin+1]; // -- dynamic allocation -- //
    for(int i=0; i<nBin+1; i++)
      arr[i] = vec[i];
  }

  TString GetTitleX( TString  varName, Bool_t noUnit = kFALSE )
  {
    TString titleX = "";
    if( varName.Contains("IPSig") )            titleX = "dxy/#sigma(dxy)";
    if( varName.Contains("pu") )               titleX = "#scale[1.0]{True number of interactions}";  // tmp
    if( varName.Contains("vtx") )              titleX = "PU";  // tmp
    if( varName.Contains("pt_trk") )           titleX = "p_{T}(trk) [GeV]";
    if( varName.Contains("eta_trk") )          titleX = "#eta(trk)";
    if( varName.Contains("pt_sim") )           titleX = "p_{T}(sim) [GeV]";
    if( varName.Contains("eta_sim") )          titleX = "#eta(sim)";
    if( varName.Contains("mass_gen") )         titleX = "m_{#mu^{+}#mu^{-}, GEN} [GeV]";
    if( varName.Contains("pt_gen") )           titleX = "p_{T}^{gen} [GeV]";
    
    if( varName.Contains("pt_offline") )           titleX = "p_{T}^{offline} [GeV]";
    if( varName.Contains("eta_offline") )          titleX = "#eta^{offline}";
    if( varName.Contains("pu_offline") )          titleX = "PU";

    if( varName.Contains("eta_gen") )          titleX = "#eta^{gen}";
    if( varName.Contains("phi_gen") )          titleX = "#phi^{gen}";
    if( varName.Contains("Pt") )               titleX = "p_{T} [GeV]";
    if( varName.Contains("AbsP") )             titleX = "p [GeV]";
    if( varName == "P" )                       titleX = "p [GeV]";
    if( varName=="LeptonP" )                   titleX = "p [GeV]";
    if( varName.Contains("Eta") )              titleX = "#eta";
    if( varName.Contains("Phi") )              titleX = "#phi";
    if( varName.Contains("Vtx") )              titleX = "# VTX";
    if( varName.Contains("Vertices") )         titleX = "# VTX";
    if( varName.Contains("Mass") )             titleX = "m_{#mu^{+}#mu^{-}} [GeV]";
    if( varName.Contains("NHits") )            titleX = "# hits";
    if( varName.Contains("NSegs") )            titleX = "# extra segments";
    if( varName.Contains("NShower") )          titleX = "# showers";
    if( varName.Contains("Shower") )           titleX = "# showers";
    if( varName.Contains("RunNum") )           titleX = "run no.";
    if( varName.Contains("RelIsoSumPt") )      titleX = "Rel Trk iso";
    if( varName.Contains("RelTkIso") )         titleX = "Rel Trk iso";

    if( varName.Contains("res_mass") )         titleX = "m (#mu^{+}#mu^{-}, GEN) [GeV]";
    if( varName.Contains("res_pt") )           titleX = "p_{T #mu^{+}#mu^{-}}(GEN) [GeV]";

    if( varName.Contains("Zpt") )              titleX = "p_{T} (#mu^{+}#mu^{-}, GEN) [GeV]";
    if( varName.Contains("Zpz") )              titleX = "p_{z} (#mu^{+}#mu^{-}, GEN) [GeV]";
    if( varName.Contains("Zy") )               titleX = "y (#mu^{+}#mu^{-}, GEN)";
    if( varName.Contains("Zphi") )             titleX = "#phi (#mu^{+}#mu^{-}, GEN)";
    if( varName == "l_pt" )                    titleX = "p_{T} (leading #mu, GEN) [GeV]";
    if( varName.Contains("l_eta") )            titleX = "#eta (leading #mu, GEN)";
    if( varName.Contains("l_phi") )            titleX = "#phi (leading #mu, GEN)";
    if( varName == "s_pt" )                    titleX = "p_{T} (sub-leading #mu, GEN) [GeV]";
    if( varName.Contains("s_eta") )            titleX = "#eta (sub-leading #mu, GEN)";
    if( varName.Contains("s_phi") )            titleX = "#phi (sub-leading #mu, GEN)";
    if( varName == "p_pt" )                    titleX = "p_{T} (#mu^{+}, GEN) [GeV]";
    if( varName.Contains("p_eta") )            titleX = "#eta (#mu^{+}, GEN)";
    if( varName.Contains("p_phi") )            titleX = "#phi (#mu^{+}, GEN)";
    if( varName == "m_pt" )                    titleX = "p_{T} (#mu^{-}, GEN) [GeV]";
    if( varName.Contains("m_eta") )            titleX = "#eta (#mu^{-}, GEN)";
    if( varName.Contains("m_phi") )            titleX = "#phi (#mu^{-}, GEN)";
    if( varName == "mu_pt" )                   titleX = "p_{T} (#mu, GEN) [GeV]";
    if( varName.Contains("mu_eta") )           titleX = "#eta (#mu, GEN)";

    if( varName.Contains("DptPt") )            titleX = "#sigma(p_{T})/p_{T}";

    if( varName.Contains("tsos_eta") )         titleX = "#eta(seed)";
    if( varName.Contains("pt_l1") )            titleX = "p_{T}(L1 muon) [GeV]";
    if( varName.Contains("eta_cut_l1") )       titleX = "#eta(L1 muon)";
    if( varName.Contains("pt_l1tk") )          titleX = "p_{T}(L1 Tk muon) [GeV]";
    if( varName.Contains("eta_cut_l1tk") )     titleX = "#eta(L1 Tk muon)";

    if( varName.Contains("pt_jet") )           titleX = "p_{T}^{jet} [GeV]";
    if( varName.Contains("eta_jet") )          titleX = "#eta^{jet}";

    if(varName.Contains("abseta"))
      titleX = "|"+titleX+"|";

    if(noUnit)
      titleX.ReplaceAll(" [GeV]", "").ReplaceAll(" (#mu)","");

    return titleX;
  }

  TString GetTitleXOld( TString  varName, Bool_t noUnit = kFALSE )
  {
    TString titleX = "";
    if( varName.Contains("Pt") )               titleX = "p_{T} (#mu) [GeV]";
    if( varName.Contains("AbsP") )             titleX = "|p| (#mu) [GeV]";
    if( varName == "P" )                       titleX = "|p| (#mu) [GeV]";
    if( varName=="LeptonP" )                   titleX = "|p| (#mu) [GeV]";
    if( varName.Contains("Eta") )              titleX = "#eta (#mu)";
    if( varName.Contains("Phi") )              titleX = "#phi (#mu)";
    if( varName.Contains("Vtx") )              titleX = "# VTX";
    if( varName.Contains("Vertices") )         titleX = "# VTX";
    if( varName.Contains("Mass") )             titleX = "m_{#mu^{+}#mu^{-}} [GeV]";
    if( varName.Contains("NHits") )            titleX = "# hits";
    if( varName.Contains("NSegs") )            titleX = "# extra segments";
    if( varName.Contains("NShower") )          titleX = "# showers";
    if( varName.Contains("Shower") )           titleX = "# showers";
    if( varName.Contains("RunNum") )           titleX = "run no.";
    if( varName.Contains("RelIsoSumPt") )      titleX = "Rel Trk iso";
    if( varName.Contains("RelTkIso") )         titleX = "Rel Trk iso";

    if( varName.Contains("res_mass") )         titleX = "m (#mu^{+}#mu^{-}, GEN) [GeV]";
    if( varName.Contains("res_pt") )           titleX = "p_{T #mu^{+}#mu^{-}}(GEN) [GeV]";

    if( varName.Contains("Zpt") )              titleX = "p_{T} (#mu^{+}#mu^{-}, GEN) [GeV]";
    if( varName.Contains("Zpz") )              titleX = "p_{z} (#mu^{+}#mu^{-}, GEN) [GeV]";
    if( varName.Contains("Zy") )               titleX = "y (#mu^{+}#mu^{-}, GEN)";
    if( varName.Contains("Zphi") )             titleX = "#phi (#mu^{+}#mu^{-}, GEN)";
    if( varName == "l_pt" )                    titleX = "p_{T} (leading #mu, GEN) [GeV]";
    if( varName.Contains("l_eta") )            titleX = "#eta (leading #mu, GEN)";
    if( varName.Contains("l_phi") )            titleX = "#phi (leading #mu, GEN)";
    if( varName == "s_pt" )                    titleX = "p_{T} (sub-leading #mu, GEN) [GeV]";
    if( varName.Contains("s_eta") )            titleX = "#eta (sub-leading #mu, GEN)";
    if( varName.Contains("s_phi") )            titleX = "#phi (sub-leading #mu, GEN)";
    if( varName == "p_pt" )                    titleX = "p_{T} (#mu^{+}, GEN) [GeV]";
    if( varName.Contains("p_eta") )            titleX = "#eta (#mu^{+}, GEN)";
    if( varName.Contains("p_phi") )            titleX = "#phi (#mu^{+}, GEN)";
    if( varName == "m_pt" )                    titleX = "p_{T} (#mu^{-}, GEN) [GeV]";
    if( varName.Contains("m_eta") )            titleX = "#eta (#mu^{-}, GEN)";
    if( varName.Contains("m_phi") )            titleX = "#phi (#mu^{-}, GEN)";
    if( varName == "mu_pt" )                   titleX = "p_{T} (#mu, GEN) [GeV]";
    if( varName.Contains("mu_eta") )           titleX = "#eta (#mu, GEN)";

    if( varName.Contains("DptPt") )            titleX = "#sigma(p_{T})/p_{T}";

    if(noUnit)
      titleX.ReplaceAll(" [GeV]", "").ReplaceAll(" (#mu)","");

    return titleX;
  }

  TH1F *Get_Cumulative( TH1F* h )
  {
    TH1F *h_cumul = (TH1F*)h->Clone("h_cumul");
    Int_t nBins = h->GetNbinsX();

    Double_t sum = 0.;
    Double_t sumW2 = 0.;
    for(Int_t i = nBins+1; i >= 0; --i) {
      sum   += h->GetBinContent(i);
      sumW2 += h->GetBinError(i) * h->GetBinError(i);
      h_cumul->SetBinContent(i, sum);
      h_cumul->SetBinError(i, sqrt(sumW2));
    }

    return h_cumul;
  }

  void RemoveUnderOverFlow( TH1F* h )
  {
    //cout << h->GetBinContent(0) << endl;
    //cout << h->GetBinError(0) << endl;
    //cout << h->GetBinContent(h->GetNbinsX()+1) << endl;
    //cout << h->GetBinError(h->GetNbinsX()+1) << endl;

    h->SetBinContent(0, 0);
    h->SetBinError(0, 0);
    h->SetBinContent( h->GetNbinsX()+1, 0);
    h->SetBinError( h->GetNbinsX()+1, 0);
  }

  void BinSizeNorm( TH1F *h, Double_t scale = 1. )
  {
    Int_t nBin = h->GetNbinsX();
    for (Int_t i=1; i<=nBin; ++i){
      h->SetBinContent( i, scale*(h->GetBinContent(i)/h->GetXaxis()->GetBinWidth(i)) );
      h->SetBinError( i, scale*(h->GetBinError(i)/h->GetXaxis()->GetBinWidth(i)) );
    }
  }

  Double_t Get_ClosestX( TH1F *h, Double_t minX, Double_t maxX, Double_t theY)
  {
    Double_t var = -99999.;
    Double_t delta = 1e9;

    Int_t nBins = h->GetNbinsX();
    for(Int_t i=1; i<=nBins; ++i) {

      Double_t x_temp = h->GetBinCenter(i);
      if( x_temp < minX || x_temp > maxX )
        continue;

      Double_t y_temp = h->GetBinContent(i);
      Double_t d_temp = fabs(theY - y_temp);
      if(d_temp < delta) {
        delta = d_temp;
        var = x_temp;
      }

    }

    return var;
  }

  Double_t Get_NonZeroMin( TH1F *h, Double_t minX = -999., Double_t maxX = -999. )
  {
    Double_t value = 1e99;

    Bool_t scan_x = kTRUE;
    if(minX==maxX)
      scan_x = kFALSE;

    Int_t nBin = h->GetNbinsX();
    for (Int_t i=1; i<=nBin; ++i){

      Double_t temp_x = h->GetBinCenter(i);
      if(scan_x && ( temp_x < minX || temp_x > maxX ) )
        continue;

      Double_t temp = h->GetBinContent(i);
      if(temp > 0 && temp < value)
        value = temp;
    }

    if( value==1e99 )
      value = 0;

    return value;
  }

  Double_t Get_NonZeroMin( TGraphAsymmErrors *g, Double_t minX = -999., Double_t maxX = -999., Double_t minDecimal = -1. )
  {
    Double_t value = 1e99;

    Bool_t scan_x = kTRUE;
    if(minX==maxX)
      scan_x = kFALSE;

    Int_t nPoint = g->GetN();
    for(Int_t i=0; i<nPoint; ++i) {
      Double_t x, y;
      g->GetPoint(i, x, y);

      if(scan_x && ( x < minX || x > maxX ) )
        continue;

      if(y > 0 && y < value)
        value = y;

    }

    if( value==1e99 )
      value = 0;

    if( minDecimal > 0 ) {
      if( ( log10(minDecimal) - roundf(log10(minDecimal)) ) == 0 ) {  // 0.1, 0.01, ...
        value = ((int)( (value/minDecimal) + 0.) * minDecimal);
      }
      else if( ( log10(minDecimal*2.) - roundf(log10(minDecimal*2.)) ) == 0 ) {  // 0.5, 0.05, ...
        Double_t temp1 = ((int)( (value/(minDecimal*2.)) + 0.)      * (minDecimal*2.));
        Double_t temp2 = ((int)( ((value*2.)/(minDecimal*2.)) + 0.) * (minDecimal));
        value = temp1 > temp2 ? temp1 : temp2;
      }
      else {
        cout << "WARNING: Get_NonZeroMin: wrong minDecimal" << endl;
      }
    }

    return value;
  }

  Double_t Get_MaxInRange( TH1F *h, Double_t minX, Double_t maxX, Double_t defualtMax = 1.5 )
  {
    Double_t value = defualtMax;
    Int_t nBin = h->GetNbinsX();
    for (Int_t i=1; i<=nBin; ++i){
      if( h->GetBinCenter(i)<minX || h->GetBinCenter(i)>maxX )
        continue;
      Double_t temp = h->GetBinContent(i);
      if(temp > value)
        value = temp;
    }

    return value;
  }

  Double_t Get_MaxInRange( TGraphAsymmErrors *g, Double_t minX, Double_t maxX, Double_t defualtMax = 1.5, Double_t minDecimal = -1. )
  {
    Double_t value = defualtMax;
    Int_t nPoint = g->GetN();
    for(Int_t i=0; i<nPoint; ++i) {
      Double_t x, y;
      g->GetPoint(i, x, y);

      if( x < minX || x > maxX )
        continue;

      if(y > value)
        value = y;
    }

    if( minDecimal > 0 ) {
      if( ( log10(minDecimal) - roundf(log10(minDecimal)) ) == 0 ) {  // 0.1, 0.01, ...
        value = ((int)( (value/minDecimal) + 1.) * minDecimal);
      }
      else if( ( log10(minDecimal*2.) - roundf(log10(minDecimal*2.)) ) == 0 ) {  // 0.5, 0.05, ...
        Double_t temp1 = ((int)( (value/(minDecimal*2.)) + 1.)      * (minDecimal*2.));
        Double_t temp2 = ((int)( ((value*2.)/(minDecimal*2.)) + 1.) * (minDecimal));
        value = temp1 < temp2 ? temp1 : temp2;
      }
      else {
        cout << "WARNING: Get_MaxInRange: wrong minDecimal" << endl;
      }
    }

    return value;
  }


  void DrawLine( TF1*& f_line, Int_t color = kRed )
  {
    f_line = new TF1("f_line", "1", -10000, 10000);
    f_line->SetLineColor(color);
    f_line->SetLineWidth(1);
    f_line->Draw("PSAME");
  }

// -- Non-trivial -- //
  TH2F *Get_ReverseAxis( TH2F *_h, TString ex = "" )
  {

    TH1F *hX = (TH1F*)_h->ProjectionX();
    Int_t nBinX = hX->GetNbinsX();
    Double_t minX = hX->GetBinLowEdge(1);
    Double_t maxX = hX->GetBinLowEdge(nBinX+1);

    TH1F *hY = (TH1F*)_h->ProjectionY();
    Int_t nBinY = hY->GetNbinsX();
    Double_t minY = hY->GetBinLowEdge(1);
    Double_t maxY = hY->GetBinLowEdge(nBinY+1);

    TH2F *h = new TH2F( (TString)_h->GetName()+"_reverse"+ex,"", nBinY, minY, maxY, nBinX, minX, maxX);
    for(Int_t ix=0; ix<=nBinX+1; ++ix) {
      for(Int_t iy=0; iy<=nBinY+1; ++iy) {
        Double_t temp_val = _h->GetBinContent(ix, iy);
        Double_t temp_err = _h->GetBinError(ix, iy);
        h->SetBinContent(iy, ix, temp_val);
        h->SetBinError(iy, ix, temp_err);
      }
    }

    return h;
  }

  TH1F* Get_FailFromTotaAndPass(TH1F* _hPass, TH1F* _hTota, TString HistName)
  {

    TH1F* hPass = (TH1F*)_hPass->Clone("hPass");
    TH1F* hTota = (TH1F*)_hTota->Clone("hTota");

    Int_t nBin = (Int_t)hPass->GetNbinsX();
    if( nBin != hTota->GetNbinsX() ) {
      cout << "Get_FailFromTotaAndPass : WARNING : hPass->GetNbinsX() != hTota->GetNbinsX()" << endl;
      return NULL;
    }

    TH1F* hFail = (TH1F*)hPass->Clone(HistName);

    for(Int_t i=0; i<=nBin+1; ++i) {
      hFail->SetBinContent(i, 0);
      hFail->SetBinError(i, 0);

      Double_t nP = hPass->GetBinContent(i);
      Double_t eP = hPass->GetBinError(i);
      Double_t nT = hTota->GetBinContent(i);
      Double_t eT = hTota->GetBinError(i);
      if(nP > nT){
        nT = nP;
        eP = eT;
        //return NULL;
      }

      Double_t nF = nT - nP;
      Double_t eF = sqrt(eT*eT - eP*eP);

      hFail->SetBinContent(i, nF);
      hFail->SetBinError(i, eF);
    }

    return hFail;
  }

  TGraphAsymmErrors* Convert_EffHistToGraph( TH1F* _h )
  {

    TH1F* h = (TH1F*)_h->Clone();
    for(Int_t i=1; i<=h->GetNbinsX(); ++i) {
      if( h->GetBinContent(i) < 0 || h->GetBinContent(i) > 1 ) {
        cout << "PlotTools||Convert_EffHistToGraph : input histogram is not efficiency" << endl;
        return NULL;
      }
    }

    TGraphAsymmErrors* g = new TGraphAsymmErrors( h );
    const Int_t nBin = g->GetN();
    for(Int_t i=0; i<nBin; i++) {

      Double_t x, y;
      g->GetPoint(i, x, y);

      Double_t ErrX_Low = g->GetErrorXlow(i);
      Double_t ErrX_High = g->GetErrorXhigh(i);
      Double_t ErrY_Low = g->GetErrorYlow(i);
      Double_t ErrY_High = g->GetErrorYhigh(i);

      Double_t ErrY_Low_New = (y-ErrY_Low) < 0 ? y : ErrY_Low;
      Double_t ErrY_High_New = (y+ErrY_High) > 1 ? (1-y) : ErrY_High;

      g->SetPointEYlow(i,ErrY_Low_New);
      g->SetPointEYhigh(i,ErrY_High_New);
    }

    return g;
  }

// -- Axis -- //
  void SetAxis_SinglePad( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle )
  {
    X_axis->SetTitle( XTitle );
    X_axis->SetLabelSize(0.04);
    X_axis->SetTitleOffset(1.1);
    X_axis->SetTitleSize(0.05);
    X_axis->SetNoExponent();
    X_axis->SetMoreLogLabels();

    Y_axis->SetTitle( YTitle );
    Y_axis->SetTitleSize(0.05);
    Y_axis->SetTitleOffset(1.25);
    Y_axis->SetLabelSize(0.04);
  }

  void SetAxis_SinglePad_Long( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle )
  {
    X_axis->SetTitle( XTitle );
    X_axis->SetLabelSize(0.04);
    X_axis->SetTitleOffset(1);
    X_axis->SetTitleSize(0.05);
    X_axis->SetNoExponent();
    X_axis->SetMoreLogLabels();

    Y_axis->SetTitle( YTitle );
    Y_axis->SetTitleSize(0.05);
    Y_axis->SetTitleOffset(0.8);
    Y_axis->SetLabelSize(0.04);
  }

  void SetAxis_TopPad( TAxis *X_axis, TAxis *Y_axis, TString YTitle )
  {
    X_axis->SetLabelFont(42);
    X_axis->SetLabelSize(0.000);
    X_axis->SetTitleSize(0.000);

    Y_axis->SetTitleFont(42);
    Y_axis->SetTitle( YTitle );
    Y_axis->SetTitleSize(0.05);
    Y_axis->SetTitleFont(42);
    Y_axis->SetTitleOffset(1.35);  // 1.25
    Y_axis->SetLabelFont(42);
    Y_axis->SetLabelSize(0.038);
  }

  void SetAxis_TopPad_Axis( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle )
  {
    X_axis->SetTitleFont(42);
    X_axis->SetTitle( XTitle );
    X_axis->SetTitleOffset(4.05);  // 1
    X_axis->SetTitleSize(0.05);  // 0.05
    X_axis->SetLabelFont(42);
    X_axis->SetLabelOffset(0.25);
    X_axis->SetLabelSize(0.04);  // 0.04
    X_axis->SetNoExponent();
    X_axis->SetMoreLogLabels();

    Y_axis->SetTitleFont(42);
    Y_axis->SetTitle( YTitle );
    Y_axis->SetTitleSize(0.05);
    Y_axis->SetTitleFont(42);
    Y_axis->SetTitleOffset(1.35);  // 1.25
    Y_axis->SetLabelFont(42);
    Y_axis->SetLabelSize(0.04);
  }

  void SetAxis_TopPad_forNM1( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle )
  {
    X_axis->SetTitleFont(42);
    X_axis->SetTitle( XTitle );
    X_axis->SetTitleOffset(4.05);  // 1
    X_axis->SetTitleSize(0.05);  // 0.05
    X_axis->SetLabelFont(42);
    X_axis->SetLabelOffset(0.172);
    X_axis->SetLabelSize(0.047);  // 0.04
    X_axis->SetNoExponent();
    X_axis->SetMoreLogLabels();

    Y_axis->SetTitleFont(42);
    Y_axis->SetTitle( YTitle );
    Y_axis->SetTitleSize(0.05);
    Y_axis->SetTitleFont(42);
    Y_axis->SetTitleOffset(1.35);  // 1.25
    Y_axis->SetLabelFont(42);
    Y_axis->SetLabelSize(0.035);
  }

  void SetAxis_BottomPad( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle, Double_t yMin = 0.51, Double_t yMax = 1.49 )
  {
    X_axis->SetMoreLogLabels();
    X_axis->SetNoExponent();
    X_axis->SetTitle( XTitle );
    X_axis->SetTitleFont(42);
    X_axis->SetTitleOffset( 0.85 );
    X_axis->SetTitleSize( 0.2 );
    X_axis->SetLabelColor(1);
    X_axis->SetLabelFont(42);
    X_axis->SetLabelOffset(0.01);
    X_axis->SetLabelSize(0.13);

    Y_axis->SetTitle( YTitle );
    Y_axis->SetTitleOffset( 0.55 );
    Y_axis->SetTitleSize( 0.12);
    Y_axis->SetTitleFont(42);
    Y_axis->SetRangeUser( yMin, yMax );
    Y_axis->SetLabelSize( 0.1 );
    Y_axis->SetLabelFont(42);
    Y_axis->SetNdivisions(505);

    // HERE tmp for NM1 single value
    // X_axis->SetLabelSize(0.11);
    // Y_axis->SetLabelSize( 0.08 );
  }

  void SetAxis_BottomPad46( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle, Double_t yMin = 0.51, Double_t yMax = 1.49 )
  {
    X_axis->SetMoreLogLabels();
    X_axis->SetNoExponent();
    X_axis->SetTitle( XTitle );
    X_axis->SetTitleFont(42);
    X_axis->SetTitleOffset( 1.2 );
    X_axis->SetTitleSize( 0.1 );
    X_axis->SetLabelColor(1);
    X_axis->SetLabelFont(42);
    X_axis->SetLabelOffset(0.02);
    X_axis->SetLabelSize(0.11);

    Y_axis->SetTitle( YTitle );
    Y_axis->SetTitleOffset( 0.65 );
    Y_axis->SetTitleSize( 0.1);  // 0.12
    Y_axis->SetTitleFont(42);
    Y_axis->SetRangeUser( yMin, yMax );
    Y_axis->SetLabelSize( 0.08 );  // 0.1
    Y_axis->SetLabelFont(42);
    Y_axis->SetNdivisions(505);
  }

  void SetAxis_BottomPad_NoLabel( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle, Double_t yMin = 0.51, Double_t yMax = 1.49, Int_t nDiv = 505 )
  {
    X_axis->SetTitle(XTitle);
    X_axis->SetTitle("");
    X_axis->SetLabelFont(42);
    X_axis->SetLabelSize(0.000);
    X_axis->SetTitleSize(0.000);

    Y_axis->SetTitle( YTitle );
    Y_axis->SetTitleOffset( 0.261 );  // 0.55
    Y_axis->SetTitleSize( 0.23 );   // 0.12
    Y_axis->SetTitleFont(42);
    Y_axis->SetRangeUser( yMin, yMax );
    Y_axis->SetLabelSize( 0.22 );  // 0.10
    Y_axis->SetLabelFont(42);
    Y_axis->SetNdivisions(nDiv);  // 505
  }

  void SetAxis_2D( TAxis *X_axis, TAxis *Y_axis, TString XTitle, TString YTitle )
  {
    X_axis->SetTitle( XTitle );
    X_axis->SetTitleSize(0.04);
    X_axis->SetTitleOffset(1.35);
    X_axis->SetLabelSize(0.04);
    // X_axis->SetNoExponent();
    // X_axis->SetMoreLogLabels();

    Y_axis->SetTitle( YTitle );
    Y_axis->SetTitleSize(0.04);
    Y_axis->SetTitleOffset(1.5);
    Y_axis->SetLabelSize(0.04);
  }

// -- Canvas -- //
  void SetCanvas_Square( TCanvas*& c, TString CanvasName, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE, Double_t SizeX = 800, Double_t SizeY = 800 )
  {
    c = new TCanvas(CanvasName, "", SizeX, SizeY);
    c->cd();

    c->SetTicky();
    c->SetTickx();
    c->SetTopMargin(0.05);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.045);
    c->SetBottomMargin(0.13);

    if( isLogx )
      c->SetLogx();
    if( isLogy )
      c->SetLogy();
  }

  void SetCanvas_Square_Long( TCanvas*& c, TString CanvasName, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE, Double_t SizeX = 800, Double_t SizeY = 800 )
  {
    c = new TCanvas(CanvasName, "", SizeX, SizeY);
    c->cd();

    c->SetTicky();
    c->SetTickx();
    c->SetTopMargin(0.05);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.045);
    c->SetBottomMargin(0.13);

    if( isLogx )
      c->SetLogx();
    if( isLogy )
      c->SetLogy();
  }

  void SetCanvas_Square2D( TCanvas*& c, TString CanvasName, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE, Double_t SizeX = 800, Double_t SizeY = 800 )
  {
    SetCanvas_Square( c, CanvasName, isLogx, isLogy, SizeX, SizeY );
    c->SetRightMargin(0.12);
  }

  void SetCanvas_Ratio( TCanvas*& c, TString CanvasName, TPad*& TopPad, TPad*& BottomPad, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE )
  {
    c = new TCanvas(CanvasName, "", 800, 800);
    c->cd();

    TopPad = new TPad("TopPad","TopPad", 0.01, 0.01, 0.99, 0.99 );
    TopPad->Draw();
    TopPad->cd();
    TopPad->SetFillStyle(0);
    TopPad->SetTicky();
    TopPad->SetTickx();
    TopPad->SetTopMargin(0.05);
    TopPad->SetLeftMargin(0.13);
    TopPad->SetRightMargin(0.045);
    TopPad->SetBottomMargin(0.3);

    if( isLogx )
      TopPad->SetLogx();
    if( isLogy )
      TopPad->SetLogy();

    c->cd();
    BottomPad = new TPad( "BottomPad", "BottomPad", 0.01, 0.01, 0.99, 0.29 );
    BottomPad->Draw();
    BottomPad->cd();
    BottomPad->SetFillStyle(0);
    BottomPad->SetTicky();
    BottomPad->SetTickx();
    BottomPad->SetGridx();
    BottomPad->SetGridy();
    BottomPad->SetTopMargin(0.05);
    BottomPad->SetBottomMargin(0.4);
    BottomPad->SetRightMargin(0.045);
    BottomPad->SetLeftMargin(0.13);

    if( isLogx )
      BottomPad->SetLogx();
  }

  void SetCanvas_Ratio( TPad*& c, TString CanvasName, TPad*& TopPad, TPad*& BottomPad, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE )
  {
    //c = new TPad(CanvasName, "", 800, 800);
    c->SetName(CanvasName);
    c->cd();

    TopPad = new TPad("TopPad","TopPad", 0.01, 0.01, 0.99, 0.99 );
    TopPad->Draw();
    TopPad->cd();
    TopPad->SetFillStyle(0);
    TopPad->SetTicky();
    TopPad->SetTickx();
    TopPad->SetTopMargin(0.05);
    TopPad->SetLeftMargin(0.13);
    TopPad->SetRightMargin(0.045);
    TopPad->SetBottomMargin(0.3);

    if( isLogx )
      TopPad->SetLogx();
    if( isLogy )
      TopPad->SetLogy();

    c->cd();
    BottomPad = new TPad( "BottomPad", "BottomPad", 0.01, 0.01, 0.99, 0.29 );
    BottomPad->Draw();
    BottomPad->cd();
    BottomPad->SetFillStyle(0);
    BottomPad->SetTicky();
    BottomPad->SetTickx();
    BottomPad->SetGridx();
    BottomPad->SetGridy();
    BottomPad->SetTopMargin(0.05);
    BottomPad->SetBottomMargin(0.4);
    BottomPad->SetRightMargin(0.045);
    BottomPad->SetLeftMargin(0.13);

    if( isLogx )
      BottomPad->SetLogx();
  }

  void SetCanvas_Ratio46( TCanvas*& c, TString CanvasName, TPad*& TopPad, TPad*& BottomPad, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE, Double_t leftMargin = 0.13 )
  {
    c = new TCanvas(CanvasName, "", 800, 800);
    c->cd();

    TopPad = new TPad("TopPad","TopPad", 0.01, 0.01, 0.99, 0.99 );
    TopPad->Draw();
    TopPad->cd();
    TopPad->SetFillStyle(0);
    TopPad->SetTicky();
    TopPad->SetTickx();
    TopPad->SetTopMargin(0.05);
    TopPad->SetLeftMargin(leftMargin);
    TopPad->SetRightMargin(0.045);
    TopPad->SetBottomMargin(0.4);

    if( isLogx )
      TopPad->SetLogx();
    if( isLogy )
      TopPad->SetLogy();

    c->cd();
    BottomPad = new TPad( "BottomPad", "BottomPad", 0.01, 0.01, 0.99, 0.39 );
    BottomPad->Draw();
    BottomPad->cd();
    BottomPad->SetFillStyle(0);
    BottomPad->SetTicky();
    BottomPad->SetTickx();
    BottomPad->SetGridx();
    BottomPad->SetGridy();
    BottomPad->SetTopMargin(0.05);
    BottomPad->SetBottomMargin(0.25);
    BottomPad->SetRightMargin(0.045);
    BottomPad->SetLeftMargin(leftMargin);

    if( isLogx )
      BottomPad->SetLogx();
  }

  void SetCanvas_RatioLong( TCanvas*& c, TString CanvasName, TPad*& TopPad, TPad*& BottomPad, Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE )
  {
    c = new TCanvas(CanvasName, "", 800, 850);
    c->cd();

    TopPad = new TPad("TopPad","TopPad", 0.01, 0.01, 0.99, 0.99 );
    TopPad->Draw();
    TopPad->cd();
    TopPad->SetFillStyle(0);
    TopPad->SetTicky();
    TopPad->SetTickx();
    TopPad->SetTopMargin(0.05);
    TopPad->SetLeftMargin(0.13);
    TopPad->SetRightMargin(0.045);
    TopPad->SetBottomMargin(0.4);

    if( isLogx )
      TopPad->SetLogx();
    if( isLogy )
      TopPad->SetLogy();

    c->cd();
    BottomPad = new TPad( "BottomPad", "BottomPad", 0.01, 0.01, 0.99, 0.39 );
    BottomPad->Draw();
    BottomPad->cd();
    BottomPad->SetFillStyle(0);
    BottomPad->SetTicky();
    BottomPad->SetTickx();
    BottomPad->SetGridx();
    BottomPad->SetGridy();
    BottomPad->SetTopMargin(0.05);
    BottomPad->SetBottomMargin(0.5);
    BottomPad->SetRightMargin(0.045);
    BottomPad->SetLeftMargin(0.13);

    if( isLogx )
      BottomPad->SetLogx();
  }

  void SetCanvas_DoubleRatioLong( TCanvas*& c, TString CanvasName, TPad*& TopPad, TPad*& BottomPad1, TPad*& BottomPad2,
                                  Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE )
  {
    c = new TCanvas(CanvasName, "", 800, 1000);  // 950
    c->cd();

    // Double_t r_top = 0.4;
    // Double_t r_bot = 0.15;
    Double_t r_top = 0.42;
    Double_t r_bot = 0.17;
    Double_t r_mid = (r_top + r_bot) / 2.;

    TopPad = new TPad("TopPad","TopPad", 0.01, 0.01, 0.99, 0.99 );
    TopPad->Draw();
    TopPad->cd();
    TopPad->SetFillStyle(0);
    TopPad->SetTicky();
    TopPad->SetTickx();
    TopPad->SetTopMargin(0.05);
    TopPad->SetBottomMargin(r_top+0.01);
    TopPad->SetLeftMargin(0.13);
    TopPad->SetRightMargin(0.045);

    if( isLogx )
      TopPad->SetLogx();
    if( isLogy )
      TopPad->SetLogy();

    c->cd();
    BottomPad1 = new TPad( "BottomPad1", "BottomPad1", 0.01, r_mid, 0.99, r_top );
    BottomPad1->Draw();
    BottomPad1->cd();
    BottomPad1->SetFillStyle(0);
    BottomPad1->SetTicky();
    BottomPad1->SetTickx();
    BottomPad1->SetGridx();
    BottomPad1->SetGridy();
    BottomPad1->SetTopMargin(0.08);  // 0.05
    BottomPad1->SetBottomMargin(0.08);  // 0.05
    BottomPad1->SetRightMargin(0.045);
    BottomPad1->SetLeftMargin(0.13);

    if( isLogx )
      BottomPad1->SetLogx();

    c->cd();
    BottomPad2 = new TPad( "BottomPad2", "BottomPad2", 0.01, r_bot, 0.99, r_mid );
    BottomPad2->Draw();
    BottomPad2->cd();
    BottomPad2->SetFillStyle(0);
    BottomPad2->SetTicky();
    BottomPad2->SetTickx();
    BottomPad2->SetGridx();
    BottomPad2->SetGridy();
    BottomPad2->SetTopMargin(0.08); // 0.05
    BottomPad2->SetBottomMargin(0.08);   // 0.05
    BottomPad2->SetRightMargin(0.045);
    BottomPad2->SetLeftMargin(0.13);

    if( isLogx )
      BottomPad2->SetLogx();
  }

  void SetCanvas_DoubleRatio( TCanvas*& c, TString CanvasName, TPad*& TopPad, TPad*& BottomPad1, TPad*& BottomPad2,
                              Bool_t isLogx = kFALSE, Bool_t isLogy = kFALSE )
  {
    c = new TCanvas(CanvasName, "", 800, 900);
    c->cd();

    Double_t r_top = 0.34;
    Double_t r_bot = 0.1;
    Double_t r_mid = (r_top + r_bot) / 2.;

    TopPad = new TPad("TopPad","TopPad", 0.01, 0.01, 0.99, 0.99 );
    TopPad->Draw();
    TopPad->cd();
    TopPad->SetFillStyle(0);
    TopPad->SetTicky();
    TopPad->SetTickx();
    TopPad->SetTopMargin(0.05);
    TopPad->SetBottomMargin(0.35);
    TopPad->SetLeftMargin(0.13);
    TopPad->SetRightMargin(0.045);

    if( isLogx )
      TopPad->SetLogx();
    if( isLogy )
      TopPad->SetLogy();

    c->cd();
    BottomPad1 = new TPad( "BottomPad1", "BottomPad1", 0.01, r_mid, 0.99, r_top );
    BottomPad1->Draw();
    BottomPad1->cd();
    BottomPad1->SetFillStyle(0);
    BottomPad1->SetTicky();
    BottomPad1->SetTickx();
    BottomPad1->SetGridx();
    BottomPad1->SetGridy();
    BottomPad1->SetTopMargin(0.08);  // 0.05
    BottomPad1->SetBottomMargin(0.08);  // 0.05
    BottomPad1->SetRightMargin(0.045);
    BottomPad1->SetLeftMargin(0.13);

    if( isLogx )
      BottomPad1->SetLogx();

    c->cd();
    BottomPad2 = new TPad( "BottomPad2", "BottomPad2", 0.01, r_bot, 0.99, r_mid );
    BottomPad2->Draw();
    BottomPad2->cd();
    BottomPad2->SetFillStyle(0);
    BottomPad2->SetTicky();
    BottomPad2->SetTickx();
    BottomPad2->SetGridx();
    BottomPad2->SetGridy();
    BottomPad2->SetTopMargin(0.08); // 0.05
    BottomPad2->SetBottomMargin(0.08);   // 0.05
    BottomPad2->SetRightMargin(0.045);
    BottomPad2->SetLeftMargin(0.13);

    if( isLogx )
      BottomPad2->SetLogx();
  }

// -- Legend -- //
  void SetLegend( TLegend *& legend, Double_t xMin = 0.75, Double_t yMin = 0.75, Double_t xMax = 0.95, Double_t yMax = 0.95, double textsize = 0.045 )
  {
    legend = new TLegend( xMin, yMin, xMax, yMax );
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont( 42 );
    legend->SetMargin(0.15);
    if(textsize > 0.)
      legend->SetTextSize(textsize);
  }

// -- Print -- //
  void Print_Histogram( TH1F* h, Bool_t NegativeCheck = kFALSE )
  {
    h->Print();

    // -- underflow -- //
    Double_t value_uf = h->GetBinContent(0);
    Double_t errorAbs_uf = h->GetBinError(0);
    Double_t errorRel_uf = value_uf == 0 ? 0 : errorAbs_uf / value_uf;

    printf( "Underflow: (value, error) = (%lf, %lf (%7.3lf %%))\n",
           value_uf, errorAbs_uf, errorRel_uf*100 );

    if( NegativeCheck && value_uf < 0 )
      printf("################## NEGATIVE BIN ##################");

    Int_t nBin = h->GetNbinsX();
    for(Int_t i=0; i<nBin; i++)
    {
      Int_t i_bin = i+1;
      Double_t LowerEdge = h->GetBinLowEdge(i_bin);
      Double_t UpperEdge = h->GetBinLowEdge(i_bin+1);

      Double_t value = h->GetBinContent(i_bin);
      Double_t errorAbs = h->GetBinError(i_bin);
      Double_t errorRel;
      if( value != 0 )
        errorRel = errorAbs / value;
      else
        errorRel = 0;

      printf( "%02d bin: [%6.1lf, %6.1lf] (value, error) = (%lf, %lf (%7.3lf %%))\n",
             i_bin, LowerEdge, UpperEdge, value, errorAbs, errorRel*100 );

      if( NegativeCheck && value < 0 )
        printf("################## NEGATIVE BIN ##################");
    }

    // -- overflow -- //
    Double_t value_of = h->GetBinContent(nBin+1);
    Double_t errorAbs_of = h->GetBinError(nBin+1);
    Double_t errorRel_of = value_of == 0 ? 0 : errorAbs_of / value_of;

    printf( "Overflow: (value, error) = (%lf, %lf (%7.3lf %%))\n",
           value_of, errorAbs_of, errorRel_of*100 );

    if( NegativeCheck && value_of < 0 )
      printf("################## NEGATIVE BIN ##################");

    printf("\n\n");
  }

  void Print_Graph( TGraphAsymmErrors* g )
  {
    TString GraphName = g->GetName();
    printf("[GraphName: %s]\n", GraphName.Data());
    Int_t nPoint = g->GetN();
    for(Int_t i=0; i<nPoint; i++)
    {
      Double_t x, y;
      g->GetPoint(i, x, y);

      Double_t xErrLow = g->GetErrorXlow(i);
      Double_t xErrHigh = g->GetErrorXhigh(i);
      Double_t LowerEdge = x - xErrLow;
      Double_t UpperEdge = x + xErrHigh;

      Double_t yErrLow = g->GetErrorYlow(i);
      Double_t yRelErrLow = yErrLow / y;
      Double_t yErrHigh = g->GetErrorYhigh(i);
      Double_t yRelErrHigh = yErrHigh / y;

      printf( "%02d point: [%6.1lf, %6.1lf] (value, errorLow, errorHigh) = (%lf, %lf (%7.3lf %%), %lf (%7.3lf %%))\n",
             i, LowerEdge, UpperEdge, y, yErrLow, yRelErrLow*100, yErrHigh, yRelErrHigh*100 );
    }
    printf("\n\n");
  }

  void Print_TH2ValLowerUpper( TH2F* h_val, TH2F* h_errlo, TH2F* h_errup, Double_t minX, Double_t maxX, Double_t minY, Double_t maxY )
  {
    Int_t nBinsX = h_val->GetNbinsX();
    Int_t nBinsY = h_val->GetNbinsY();

    TH1F *hX = (TH1F*)h_val->ProjectionX();
    TH1F *hY = (TH1F*)h_val->ProjectionY();

    vector<Double_t> vec_Xrange, vec_Yrange;

    Bool_t isFirst = kTRUE;
    for(Int_t iy=1; iy<=nBinsY; ++iy) {
      Double_t temp_y   = hY->GetBinCenter(iy);
      Double_t temp_ylo = hY->GetBinLowEdge(iy);
      Double_t temp_yup = hY->GetBinLowEdge(iy+1);
      if(temp_y<minY || temp_y>maxY)
        continue;

      vec_Yrange.push_back(temp_ylo);

      for(Int_t ix=1; ix<=nBinsX; ++ix) {
        Double_t temp_x   = hX->GetBinCenter(ix);
        Double_t temp_xlo = hX->GetBinLowEdge(ix);
        Double_t temp_xup = hX->GetBinLowEdge(ix+1);
        if(temp_x<minX || temp_x>maxX)
          continue;

        if(isFirst) {
          vec_Xrange.push_back(temp_xlo);
        }

        //-- print
        Double_t temp_val = h_val->GetBinContent(ix, iy);
        Double_t temp_lo  = h_errlo->GetBinContent(ix, iy);
        Double_t temp_up  = h_errup->GetBinContent(ix, iy);
        cout << TString::Format("%.3f+%.3f-%.3f", temp_val, temp_up, temp_lo) << "\t";

      }
      cout << "\n";
      isFirst = kFALSE;
    }
    cout << "\n";
    cout << "X range: " << endl;
    for(Int_t i=0; i<(int)vec_Xrange.size(); ++i)
      cout << "\t" << vec_Xrange[i];
    cout << "\n";

    cout << "Y range: " << endl;
    for(Int_t i=0; i<(int)vec_Yrange.size(); ++i)
      cout << "\n" << vec_Yrange[i];
    cout << "\n\n";

  }

// -- Latex -- //
  void Latex_NoPreliminary( TLatex &latex )
  {
    latex.DrawLatexNDC(0.13, 0.96, "#font[62]{CMS}");
  }

  void Latex_Preliminary_TopRight( TLatex &latex )
  {
    latex.DrawLatexNDC(0.13, 0.96, "#font[62]{CMS}");
    latex.DrawLatexNDC(0.75, 0.89, "#font[62]{#bf{#it{#scale[0.7]{ Preliminary}}}}");
  }

  void Latex_Preliminary_NoDataInfo( TLatex &latex )
  {
    latex.DrawLatexNDC(0.13, 0.96, "#font[62]{CMS}#font[42]{#it{#scale[0.8]{ Preliminary}}}");
  }

  void Latex_Preliminary_13TeV( TLatex &latex )
  {
    Latex_Preliminary_NoDataInfo( latex );
    latex.DrawLatexNDC(0.82, 0.96, "#font[42]{#scale[0.8]{13 TeV}}");
  }

  void Latex_Preliminary( TLatex &latex, Double_t lumi  )
  {
    Latex_Preliminary_NoDataInfo( latex );
    latex.DrawLatexNDC(0.70, 0.96, "#font[42]{#scale[0.8]{"+TString::Format("%.2lf", lumi)+" fb^{-1} (13 TeV)}}");
  }

  void Latex_Preliminary( TLatex &latex, Double_t lumi, Int_t E_CM  )
  {
    Latex_Preliminary_NoDataInfo( latex );
    latex.DrawLatexNDC(0.69, 0.96, "#font[42]{#scale[0.8]{"+TString::Format("%.1lf fb^{-1} (%d TeV)", lumi, E_CM)+"}}");
  }

  void Latex_Simulation( TLatex &latex, TString offset = "" )
  {
    latex.DrawLatexNDC(0.82, 0.96, offset+"#font[42]{#scale[0.8]{13 TeV}}");
    latex.DrawLatexNDC(0.13, 0.96, offset+"#font[62]{CMS}#font[42]{#it{#scale[0.8]{ Simulation}}}");
  }

  void Latex_Simulation_14TeV( TLatex &latex, TString offset = "" )
  {
    latex.DrawLatexNDC(0.82, 0.96, offset+"#font[42]{#scale[0.8]{14 TeV}}");
    latex.DrawLatexNDC(0.13, 0.96, offset+"#font[62]{CMS}#font[42]{#it{#scale[0.8]{ Simulation}}}");
  }

  void Latex_Preliminary_EffPlotStyle( TLatex &latex, Int_t year, Int_t E_CM = 13 )
  {
    latex.DrawLatexNDC(0.70, 0.96, "#font[42]{#scale[0.8]{"+TString::Format("%d, %d TeV", year, E_CM)+"}}");
    latex.DrawLatexNDC(0.13, 0.96, "#font[62]{CMS}");
    latex.DrawLatexNDC(0.24, 0.96, "#font[42]{#it{#scale[0.8]{Preliminary}}}");
  }

  void Latex_Geant4( TLatex &latex, TString target, TString particle, Double_t energy, TString unit )
  {
    TString str_energy = TString::Format("%d", (int)energy);
    latex.DrawLatexNDC(0.82, 0.96, "#font[42]{#scale[0.8]{#font[62]{"+particle+"},#font[42]{ "+str_energy+" "+unit+"}}}");
    latex.DrawLatexNDC(0.17, 0.96, "#font[62]{Geant4}#font[42]{#it{#scale[0.8]{ Simulation}}}");
    latex.DrawLatexNDC(0.50, 0.97, "#font[62]{#scale[0.7]{"+target+"}}");
  }



// -- Error -- //
  pair<Double_t, Double_t> GetClopperPearsonErrorLU( Double_t nPass, Double_t nTota, Double_t alpha = 1-0.6827, bool equal_tailed = true )
  {
    Double_t alpha_min = -999.;
    if( equal_tailed )
      alpha_min = alpha/2;
    else
      alpha_min = alpha;

    Double_t Eff = 0;
    Double_t lower = 0;
    Double_t upper = 1;
    if( nPass > 0 ) {
      Eff = ( nPass / nTota );
      lower = ROOT::Math::beta_quantile(alpha_min, nPass, nTota - nPass + 1);
    }
    if( nTota - nPass > 0 ) {
      Eff = ( nPass / nTota );
      upper = ROOT::Math::beta_quantile_c(alpha_min, nPass + 1, nTota - nPass);
    }
    else if( nTota == 0 ) {
      Eff = 0;
      lower = 0;
      upper = 0;
    }

    return make_pair(Eff-lower, upper-Eff);
  }

  pair<Double_t, Double_t> GetWeightedBinomialErrorLU( Double_t nPass, Double_t ePass, Double_t nFail, Double_t eFail )
  {
    Double_t nTota = nPass + nFail;
    Double_t Eff = nPass / nTota;
    Double_t Err = 1./(nTota*nTota) * sqrt( nPass*nPass* eFail*eFail + nFail*nFail * ePass*ePass );
    
    Double_t ErrL = Eff - Err < 0 ? Eff : Err;
    Double_t ErrU = Eff + Err > 1 ? 1 - Eff : Err;
    return make_pair(ErrL, ErrU);
  }

  Double_t GetUncorrelatedError(Double_t A, Double_t sigma_A, Double_t B, Double_t sigma_B)
  {
    Double_t ratio_A, ratio_B;
    if(sigma_A<=0)
      ratio_A = 0;
    else
      ratio_A = (sigma_A) / A;

    if(sigma_B<=0)
      ratio_B = 0;
    else
      ratio_B = (sigma_B) / B;

    Double_t errorSquare = ratio_A * ratio_A + ratio_B * ratio_B;

    return (A/B) * sqrt(errorSquare);
  }

/*
  TH1F* QuadSum_NoError( TH1F* h1, TH1F* h2 )
  {
    TH1F* h_QuadSum = (TH1F*)h1->Clone( "h_QuadSum" );
    Int_t nBin = h1->GetNbinsX();
    for(Int_t i=0; i<nBin; i++)
    {
      Int_t i_bin = i+1;

      Double_t value1 = h1->GetBinContent(i_bin);
      Double_t value2 = h2->GetBinContent(i_bin);

      Double_t QuadSum = sqrt( value1*value1 + value2*value2 );

      h_QuadSum->SetBinContent( i_bin, QuadSum );
      h_QuadSum->SetBinError( i_bin, 0 );
    }
    return h_QuadSum;
  }

  void AssignErrors( TH1F* h_cv, TH1F* h_RelUnc, Bool_t isPercent = kFALSE )
  {
    if( h_cv->GetNbinsX() != h_RelUnc->GetNbinsX() )
    {
      printf("# bins for central value and relative uncertainty histograms are not same! .. %d != %d\n",
        h_cv->GetNbinsX(), h_RelUnc->GetNbinsX() );

      return;
    }

    Int_t nBin = h_cv->GetNbinsX();
    for(Int_t i=0; i<nBin; i++)
    {
      Int_t i_bin = i+1;

      Double_t value = h_cv->GetBinContent(i_bin);
      Double_t RelUnc = h_RelUnc->GetBinContent(i_bin);
      if( isPercent ) RelUnc = RelUnc / 100.0;

      Double_t AbsUnc = value * RelUnc;

      h_cv->SetBinError(i_bin, AbsUnc);
    }
  }

  void AssignErrors_Graph( TGraphAsymmErrors* g, TH1F* h_RelUnc, Bool_t isPercent = kFALSE )
  {
    Int_t nPoint = g->GetN();
    Int_t nBin = h_RelUnc->GetNbinsX();
    if( nPoint != nBin )
    {
      printf("[nPoint != nBin: %d, %d]\n", nPoint, nBin);
      return;
    }

    for(Int_t i=0; i<nPoint; i++)
    {
      Double_t x = 0;
      Double_t y = 0;
      g->GetPoint( i, x, y );

      Int_t i_bin = i+1;
      Double_t RelUnc = h_RelUnc->GetBinContent(i_bin);
      if( isPercent ) RelUnc = RelUnc / 100.0;

      Double_t AbsUnc = y * RelUnc;

      g->SetPointEYlow(i, AbsUnc );
      g->SetPointEYhigh(i, AbsUnc );
    }
  }

  TH1F* Convert_GraphToHist( TGraphAsymmErrors *g )
  {
    const Int_t nBin = g->GetN();
    Double_t *BinEdges = new Double_t[nBin+1];
    Double_t *value = new Double_t[nBin];
    Double_t *error = new Double_t[nBin];

    for(Int_t i=0; i<nBin; i++)
    {
      Double_t x, y;
      g->GetPoint(i, x, y);

      // -- make BinEdges array -- //
      Double_t ErrX_Low = g->GetErrorXlow(i);
      Double_t ErrX_High = g->GetErrorXhigh(i);

      if( i == nBin-1 )
      {
        BinEdges[i] = x - ErrX_Low;
        BinEdges[i+1] = x + ErrX_High;
      }
      else
        BinEdges[i] = x - ErrX_Low;


      // -- store graph information -- //
      value[i] = y;

      Double_t ErrY_Low = g->GetErrorYlow(i);
      Double_t ErrY_High = g->GetErrorYhigh(i);

      // -- take the larger one -- //
      error[i] = ErrY_Low > ErrY_High ? ErrY_Low : ErrY_High;
    }

    TString GraphName = g->GetName();
    TH1F* h_temp = new TH1F( "h_"+GraphName, "", nBin, BinEdges );

    // -- fill this histogram using graph information -- //
    for(Int_t i=0; i<nBin; i++)
    {
      Int_t i_bin = i+1;
      h_temp->SetBinContent( i_bin, value[i] );
      h_temp->SetBinError( i_bin, error[i] );
    }

    return h_temp;
  }

  TH1F* Extract_RelUnc( TH1F* h, TString HistName = "", Bool_t ConvertToPercent = kFALSE )
  {
    TH1F* h_RelUnc = (TH1F*)h->Clone();
    if( HistName != "" )
      h_RelUnc->SetName(HistName);

    Int_t nBin = h->GetNbinsX();
    for(Int_t i=0; i<nBin; i++)
    {
      Int_t i_bin = i+1;

      Double_t value = h->GetBinContent(i_bin);
      Double_t error = h->GetBinError(i_bin);

      Double_t RelUnc = error / value;
      if( ConvertToPercent )
        RelUnc = RelUnc * 100;

      h_RelUnc->SetBinContent(i_bin, RelUnc );
      h_RelUnc->SetBinError(i_bin, 0);
    }

    return h_RelUnc;
  }

  TH1F* ConvertHist_AbsToRel( TH1F* h_CenV, TH1F* h_AbsUnc, Bool_t ConvertToPercent = kFALSE )
  {
    TH1F* h_RelUnc = (TH1F*)h_AbsUnc->Clone();

    Int_t nBin = h_CenV->GetNbinsX();
    for(Int_t i=0; i<nBin; i++)
    {
      Int_t i_bin = i+1;

      Double_t CenV = h_CenV->GetBinContent(i_bin);
      Double_t AbsUnc = h_AbsUnc->GetBinContent(i_bin);
      Double_t RelUnc = AbsUnc / CenV;
      if( ConvertToPercent )
        RelUnc = RelUnc * 100;

      h_RelUnc->SetBinContent(i_bin, RelUnc );
      h_RelUnc->SetBinError(i_bin, 0);
    }

    return h_RelUnc;
  }
*/

// -- Old -- //
/*
  TH1F* DivideEachBin_ByBinWidth( TH1F* h, TString HistName = "" )
  {
    TH1F* h_return = (TH1F*)h->Clone();
    if( HistName != "" )
      h_return->SetName(HistName);

    Int_t nBin = h->GetNbinsX();
    for(Int_t i=0; i<nBin; i++)
    {
      Int_t i_bin = i+1;
      Double_t Entry_before = h->GetBinContent(i_bin);
      Double_t Error_before = h->GetBinError(i_bin);
      Double_t BinWidth = h->GetBinWidth(i_bin);

      Double_t Entry_after = Entry_before / BinWidth;
      Double_t Error_after = Error_before / BinWidth;

      h_return->SetBinContent(i_bin, Entry_after);
      h_return->SetBinError(i_bin, Error_after);
    }

    return h_return;
  }

  TH1F* MultiplyEachBin_byBinWidth( TH1F* h, TString HistName = "" )
  {
    TH1F* h_return = (TH1F*)h->Clone();
    if( HistName != "" )
      h_return->SetName(HistName);

    Int_t nBin = h->GetNbinsX();
    for(Int_t i=0; i<nBin; i++)
    {
      Int_t i_bin = i+1;
      Double_t Entry_before = h->GetBinContent(i_bin);
      Double_t Error_before = h->GetBinError(i_bin);
      Double_t BinWidth = h->GetBinWidth(i_bin);

      Double_t Entry_after = Entry_before * BinWidth;
      Double_t Error_after = Error_before * BinWidth;

      h_return->SetBinContent(i_bin, Entry_after);
      h_return->SetBinError(i_bin, Error_after);
    }

    return h_return;
  }

  void SaveAsHist_OneContent( Double_t content, TString HistName, TFile *f_output )
  {
    TH1F *h = new TH1F(HistName, "", 1, 0, 1);
    h->SetBinContent(1, content);
    h->SetBinError(1, 0);
    f_output->cd();
    h->Write();
  }

  void SaveAsHist_OneContent_WithError( Double_t content, Double_t error, TString HistName, TFile *f_output )
  {
    TH1F *h = new TH1F(HistName, "", 1, 0, 1);
    h->SetBinContent(1, content);
    h->SetBinError(1, error);
    f_output->cd();
    h->Write();
  }

  Double_t GetContent_OneBinHist( TString FileName, TString HistName )
  {
    TH1F* h_temp = Get_Hist( FileName, HistName );
    if( h_temp->GetNbinsX() != 1 )
    {
      cout << "This histogram has more than 1 bin! ... please check. Return -999" << endl;
      return -999;
    }

    return h_temp->GetBinContent(1);
  }
*/



// -- SampleInfo -- //
  class SampleInfo
  {
  public:
    TString name; // -- short name for the convenience in the code: ex> DYMuMu_M50 -- //
    TString fullName; // -- for legend: ex> DY #rightarrow Z/#gamma* -- //

    Bool_t isRealData;
    TString fileName;
    Int_t color;
    Int_t markerStyle;

    Double_t xSec;
    Double_t sumWeight; // -- same with # events for the samples without negative weights -- //
    Double_t normFactor;

    SampleInfo()
    {
      this->Init();
    }

    SampleInfo( TString _fullName, Int_t _color, Int_t _markerStyle )
    {
      this->fullName = _fullName;
      this->SetColor( _color );
      this->SetMarkerStyle( _markerStyle );
    }

    SampleInfo( Bool_t _isRealData, TString _name, TString _fullName )
    {
      this->isRealData = _isRealData;
      this->SetName( _name, _fullName );
    }

    void SetName( TString _name, TString _fullName )
    {
      this->name = _name;
      this->fullName = _fullName;
    }

    void SetFileName( TString _name )
    {
      this->fileName = _name;
    }

    void SetColor( Int_t _color )
    {
      this->color = _color;
    }

    void SetMarkerStyle( Int_t _markerStyle )
    {
      this->markerStyle = _markerStyle;
    }

    void SetNormFactor( Double_t lumi, Double_t _xSec, Double_t _sumWeight )
    {
      if( this->isRealData )
      {
        cout << "[SetNormFactor] This is real data! ... something goes wrong" << endl;
        return;
      }

      this->xSec = _xSec;
      this->sumWeight = _sumWeight;
      this->normFactor = (lumi * this->xSec) / this->sumWeight;
      printf("[SetNormFactor] Sample: %s\n", name.Data() );
      printf("\tNormalization factor = (%.3lf * %.3e) / (%.1lf) = %.3e\n", lumi, this->xSec, this->sumWeight, this->normFactor);
    }

  private:
    void Init()
    {
      this->name = "";
      this->fullName = "";

      this->isRealData = kFALSE;
      this->fileName = "";
      this->color = 0;
      this->markerStyle = 20;

      this->xSec = -999;
      this->sumWeight = -999;
      this->normFactor = -999;
    }
  };

// -- HistInfo -- //
  class HistInfo
  {
  public:
    TString name;
    TString titleX;
    TString titleY;

    Bool_t hasXRange;
    Double_t minX;
    Double_t maxX;

    Bool_t hasYRange;
    Double_t minY;
    Double_t maxY;

    Bool_t hasZRange;
    Double_t minZ;
    Double_t maxZ;

    Double_t hasRebinX;
    Double_t nRebinX;

    Double_t hasRebinY;
    Double_t nRebinY;

    Bool_t isFilled;

    // -- log axis: property of canvases, not histograms -- //
    // Bool_t isLogX;
    // Bool_t isLogY;

    HistInfo()
    {
      this->Init();
    }

    HistInfo( TString _name ): HistInfo()
    {
      this->name = _name;
    }

    HistInfo( TString _name, TString _titleX, TString _titleY ): HistInfo()
    {
      this->name = _name;
      this->SetTitle( _titleX, _titleY );
    }

    HistInfo( TString _titleX, TString _titleY ): HistInfo()
    {
      this->SetTitle( _titleX, _titleY );
    }

    HistInfo( TString _titleX, TString _titleY, Double_t minX, Double_t maxX, Double_t minY, Double_t maxY ): HistInfo()
    {
      this->SetTitle( _titleX, _titleY );
      this->SetXRange( minX, maxX );
      this->SetYRange( minY, maxY );
    }

    void SetTitle( TString X, TString Y )
    {
      this->titleX = X;
      this->titleY = Y;
    }

    // void SetLogAxis( Bool_t X, Bool_t Y )
    // {
    //  this->isLogX = X;
    //  this->isLogY = Y;
    // }

    void SetXRange( Double_t min, Double_t max )
    {
      this->hasXRange = kTRUE;
      this->minX = min;
      this->maxX = max;
    }

    void SetYRange( Double_t min, Double_t max )
    {
      this->hasYRange = kTRUE;
      this->minY = min;
      this->maxY = max;
    }

    void SetZRange( Double_t min, Double_t max )
    {
      this->hasZRange = kTRUE;
      this->minZ = min;
      this->maxZ = max;
    }

    void SetRebinX( Int_t _nRebin )
    {
      this->hasRebinX = kTRUE;
      this->nRebinX = _nRebin;
    }

    void SetRebinY( Int_t _nRebin )
    {
      this->hasRebinY = kTRUE;
      this->nRebinY = _nRebin;
    }

    void IsFilled( Bool_t _isFilled = kTRUE )
    {
      this->isFilled = _isFilled;
    }

  private:
    void Init()
    {
      this->name = "";
      this->titleX = "";
      this->titleY = "";

      this->hasXRange = kFALSE;
      this->minX = -999;
      this->maxX = -999;

      this->hasYRange = kFALSE;
      this->minY = -999;
      this->maxY = -999;

      this->hasZRange = kFALSE;
      this->minZ = -999;
      this->maxZ = -999;

      this->hasRebinX = kFALSE;
      this->nRebinX = 1;

      this->hasRebinY = kFALSE;
      this->nRebinY = 1;

      this->isFilled = kFALSE;

      // this->isLogX = kFALSE;
      // this->isLogY = kFALSE;
    }

  };

// -- TH1 extension -- //
  class TH1Ext
  {
  public:
    TH1F* h;

    Bool_t hasRatio;
    TH1F* h_ratio;

    SampleInfo* sampleInfo;
    HistInfo* histInfo;

    TH1Ext()
    {
      TH1::AddDirectory( kFALSE );

      this->h = NULL;

      this->hasRatio = kFALSE;
      this->h_ratio = NULL;
    }

    TH1Ext( SampleInfo* _sampleInfo, HistInfo* _histInfo, TString histName = "" ): TH1Ext()
    {
      this->sampleInfo = _sampleInfo;
      this->histInfo = _histInfo;

      if( histName == "" )
        this->h = Get_Hist( this->sampleInfo->fileName, this->histInfo->name );
      else
        this->h = Get_Hist( this->sampleInfo->fileName, histName );

      // -- it should be done earlier: to be consistent with the ratio calculation -- //
      if( this->histInfo->hasRebinX )
        h->Rebin( this->histInfo->nRebinX );
    }

    TH1Ext( SampleInfo* _sampleInfo, HistInfo* _histInfo, TH1F* _h ): TH1Ext()
    {
      this->sampleInfo = _sampleInfo;
      this->histInfo = _histInfo;

      this->h = (TH1F*)_h->Clone();

      // -- it should be done earlier: to be consistent with the ratio calculation -- //
      if( this->histInfo->hasRebinX )
        h->Rebin( this->histInfo->nRebinX );
    }

    void DrawAndSet( TString drawOp )
    {
      this->h->Draw( drawOp );
      this->SetAttributes(); // -- setting after drawing: to be consistent with TGraphExt case -- //
    }

    void DrawRatioAndSet( TString DrawOp, TString ratioTitle, Double_t minRatio = 0.5, Double_t maxRatio = 1.5 )
    {
      this->h_ratio->Draw( DrawOp );
      this->SetAttributesRatio(ratioTitle, minRatio, maxRatio);
    }

    void AddToLegend( TLegend *legend, Bool_t addBold = kFALSE, TString Opt = "lep" )
    {
      if(addBold)
        legend->AddEntry( this->h, "#bf{"+this->sampleInfo->fullName+"}", Opt );
      else
        legend->AddEntry( this->h, this->sampleInfo->fullName, Opt );
    }

    void CalcRatio_DEN( TH1F* h_DEN )
    {
      this->hasRatio = kTRUE;

      if( h == NULL )
      {
        cout << "Histogram is not assigned yet!" << endl;
        return;
      }

      gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
      h->Sumw2();
      h_DEN->Sumw2();
      gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

      this->h_ratio = (TH1F*)this->h->Clone();
      h_ratio->Divide( this->h, h_DEN );
    }

    void CalcRatio_NUM( TH1F* h_NUM )
    {
      this->hasRatio = kTRUE;

      if( h == NULL )
      {
        cout << "Histogram is not assigned yet!" << endl;
        return;
      }

      gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
      h->Sumw2();
      h_NUM->Sumw2();
      gROOT->ProcessLine( "gErrorIgnoreLevel = kPrint;");

      this->h_ratio = (TH1F*)this->h->Clone();
      h_ratio->Divide( h_NUM, this->h );
    }

    void SetAttributes()
    {
      this->h->SetTitle("");
      this->h->SetStats(kFALSE);

      this->h->SetLineColor( this->sampleInfo->color );
      this->h->SetFillColorAlpha( kWhite, 0 );
      this->h->SetMarkerStyle( this->sampleInfo->markerStyle );
      this->h->SetMarkerSize( 1.3 );
      this->h->SetMarkerColor( this->sampleInfo->color );
      if( this->histInfo->isFilled )
      {
        this->h->SetLineColorAlpha( kWhite, 0 );
        this->h->SetMarkerColorAlpha( kWhite, 0 );
        this->h->SetFillColorAlpha( this->sampleInfo->color, 1 );
      }

      if( this->histInfo->hasXRange )
        h->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );

      if( this->histInfo->hasYRange )
        h->GetYaxis()->SetRangeUser( this->histInfo->minY, this->histInfo->maxY );

      if( this->histInfo->hasZRange )
        h->GetZaxis()->SetRangeUser( this->histInfo->minZ, this->histInfo->maxZ );

      if( this->hasRatio )
        SetAxis_TopPad( this->h->GetXaxis(), this->h->GetYaxis(), this->histInfo->titleY );
      else
        SetAxis_SinglePad( this->h->GetXaxis(), this->h->GetYaxis(), this->histInfo->titleX, this->histInfo->titleY );
    }

    void SetAttributesRatio(TString ratioTitle, Double_t minRatio, Double_t maxRatio)
    {
      if( this->h_ratio == NULL ) return;

      this->h_ratio->SetTitle("");
      this->h_ratio->SetStats(kFALSE);

      this->h_ratio->SetLineColor( this->sampleInfo->color );
      this->h_ratio->SetFillColorAlpha( kWhite, 0 );
      if( this->histInfo->isFilled )
      {
        this->h_ratio->SetLineColorAlpha( kWhite, 0 );
        this->h_ratio->SetMarkerColorAlpha( kWhite, 0 );
        this->h_ratio->SetFillColorAlpha( this->sampleInfo->color, 1 );
      }
      this->h_ratio->SetMarkerStyle( this->sampleInfo->markerStyle );
      this->h_ratio->SetMarkerColor( this->sampleInfo->color );

      if( this->histInfo->hasXRange ) {
        h_ratio->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );
      }

      SetAxis_BottomPad( this->h_ratio->GetXaxis(), this->h_ratio->GetYaxis(), this->histInfo->titleX, ratioTitle, minRatio, maxRatio );
    }
  };

// -- TH2 extension -- //
  class TH2Ext
  {
  public:
    TH2F* h;

    SampleInfo* sampleInfo;
    HistInfo* histInfo;

    TH2Ext()
    {
      TH1::AddDirectory( kFALSE );

      this->h = NULL;
    }

    TH2Ext( SampleInfo* _sampleInfo, HistInfo* _histInfo, TString histName = "" ): TH2Ext()
    {
      this->sampleInfo = _sampleInfo;
      this->histInfo = _histInfo;

      if( histName == "" )
        this->h = Get_Hist_2D( this->sampleInfo->fileName, this->histInfo->name );
      else
        this->h = Get_Hist_2D( this->sampleInfo->fileName, histName );
    }

    // -- when the histogram already exists -- //
    TH2Ext( SampleInfo* _sampleInfo, HistInfo* _histInfo, TH2F* _h ): TH2Ext()
    {
      this->sampleInfo = _sampleInfo;
      this->histInfo = _histInfo;
      this->h = (TH2F*)_h->Clone();
    }

    void DrawAndSet( TString drawOp )
    {
      this->h->Draw( drawOp );
      this->SetAttributes( drawOp ); // -- setting after drawing: to be consistent with TGraphExt case -- //
    }

    void AddToLegend( TLegend *legend )
    {
      legend->AddEntry( this->h, this->sampleInfo->fullName );
    }

  protected:
    void SetAttributes( TString drawOp )
    {
      this->h->SetTitle("");
      this->h->SetStats(kFALSE);

      if( drawOp.Contains("SCAT") )
      {
        this->h->SetMarkerStyle( this->sampleInfo->markerStyle );
        this->h->SetLineColorAlpha( kWhite, 0 );
        this->h->SetMarkerColor( this->sampleInfo->color );
      }

      if( this->histInfo->hasXRange )
        h->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );

      if( this->histInfo->hasYRange )
        h->GetYaxis()->SetRangeUser( this->histInfo->minY, this->histInfo->maxY );

      if( this->histInfo->hasZRange )
        h->GetZaxis()->SetRangeUser( this->histInfo->minZ, this->histInfo->maxZ );

      if( this->histInfo->hasRebinX )
        h->RebinX( this->histInfo->nRebinX );

      if( this->histInfo->hasRebinY )
        h->RebinY( this->histInfo->nRebinY );

      SetAxis_2D( this->h->GetXaxis(), this->h->GetYaxis(), this->histInfo->titleX, this->histInfo->titleY );
    }
  };

// -- TGraph extension -- //
  class TGraphExt
  {
  public:
    TGraphAsymmErrors* g;
    TH1F *hFrame;
    TH1F *hFrame_ratio;

    Bool_t hasRatio;
    TGraphAsymmErrors* g_ratio;

    SampleInfo* sampleInfo;
    HistInfo* histInfo;

    TGraphExt()
    {
      this->g = NULL;

      this->hasRatio = kFALSE;
      this->g_ratio = NULL;
    }

    TGraphExt( SampleInfo* _sampleInfo, HistInfo* _histInfo, TString graphName = "" ): TGraphExt()
    {
      this->sampleInfo = _sampleInfo;
      this->histInfo = _histInfo;

      if( graphName == "" )
        this->g = Get_Graph( this->sampleInfo->fileName, this->histInfo->name );
      else
        this->g = Get_Graph( this->sampleInfo->fileName, graphName );

      if(!this->histInfo->hasXRange)
        this->hFrame = (TH1F*)((this->g)->GetHistogram())->Clone("hFrame");
      else
        this->hFrame = new TH1F("hFrame", "", 1, this->histInfo->minX, this->histInfo->maxX );
      this->hFrame->SetLineWidth(0);
      this->hFrame->SetStats(0);
      this->hFrame->SetTitle(0);

      if(!this->histInfo->hasXRange)
        this->hFrame_ratio = (TH1F*)((this->g)->GetHistogram())->Clone("hFrame_ratio");
      else
        this->hFrame_ratio = new TH1F("hFrame_ratio", "", 1, this->histInfo->minX, this->histInfo->maxX );
      this->hFrame_ratio->SetLineWidth(0);
      this->hFrame_ratio->SetStats(0);
      this->hFrame_ratio->SetTitle(0);
    }

    TGraphExt( SampleInfo* _sampleInfo, HistInfo* _histInfo, TGraphAsymmErrors* _g ): TGraphExt()
    {
      this->sampleInfo = _sampleInfo;
      this->histInfo = _histInfo;
      this->g = (TGraphAsymmErrors*)_g->Clone();

      if(!this->histInfo->hasXRange)
        this->hFrame = (TH1F*)((this->g)->GetHistogram())->Clone("hFrame");
      else
        this->hFrame = new TH1F("hFrame", "", 1, this->histInfo->minX, this->histInfo->maxX );
      this->hFrame->SetLineWidth(0);
      this->hFrame->SetStats(0);
      this->hFrame->SetTitle(0);

      if(!this->histInfo->hasXRange)
        this->hFrame_ratio = (TH1F*)((this->g)->GetHistogram())->Clone("hFrame_ratio");
      else
        this->hFrame_ratio = new TH1F("hFrame_ratio", "", 1, this->histInfo->minX, this->histInfo->maxX );
      this->hFrame_ratio->SetLineWidth(0);
      this->hFrame_ratio->SetStats(0);
      this->hFrame_ratio->SetTitle(0);
    }

    void DrawAndSet( TString drawOp, Bool_t isFrame = kFALSE, Bool_t drawAxis = kFALSE )
    {
      if(isFrame) this->hFrame->Draw();
      this->g->Draw( drawOp );
      this->SetAttributes( drawAxis );
    }

    void CalcRatio_DEN( TGraphAsymmErrors* g_DEN, TString errorPropagation = "AoverB" )
    {
      this->g_ratio = this->MakeRatioGraph( g, g_DEN, errorPropagation );
      //this->g_ratio->SetHistogram( this->hFrame_ratio );
      this->hasRatio = kTRUE;
    }

    void CalcRatio_NUM( TGraphAsymmErrors* g_NUM, TString errorPropagation = "AoverB" )
    {
      this->g_ratio = this->MakeRatioGraph( g_NUM, g, errorPropagation );
      //this->g_ratio->SetHistogram( this->hFrame_ratio );
      this->hasRatio = kTRUE;
    }

    void CalcDynamicRatio_DEN( TGraphAsymmErrors* g_DEN, TString errorPropagation = "AoverB" )
    {
      this->g_ratio = this->MakeDynamicRatioGraph( g, g_DEN, errorPropagation );
      //this->g_ratio->SetHistogram( this->hFrame_ratio );
      this->hasRatio = kTRUE;
    }

    void CalcDynamicRatio_NUM( TGraphAsymmErrors* g_NUM, TString errorPropagation = "AoverB" )
    {
      this->g_ratio = this->MakeDynamicRatioGraph( g_NUM, g, errorPropagation );
      //this->g_ratio->SetHistogram( this->hFrame_ratio );
      this->hasRatio = kTRUE;
    }

    void DrawRatioAndSet( TString DrawOp, TString ratioTitle, Double_t minRatio = 0.5, Double_t maxRatio = 1.5,
                          Bool_t isFrame = kFALSE, Bool_t isNoLabel = kFALSE, Int_t nDiv = 505, Bool_t is46 = kFALSE )
    {
      if(isFrame) this->hFrame_ratio->Draw();
      this->g_ratio->Draw( DrawOp );
      this->SetAttributesRatio(ratioTitle, minRatio, maxRatio, isNoLabel, nDiv, is46 );
    }

    void AddToLegend( TLegend *legend, Bool_t addBold = kFALSE, TString Opt = "lep" )
    {
      if(addBold)
        legend->AddEntry( this->g, "#bf{"+this->sampleInfo->fullName+"}", Opt );
      else
        legend->AddEntry( this->g, this->sampleInfo->fullName, Opt );
    }

  protected:
    void SetAttributes( Bool_t drawAxis = kFALSE )
    {
      this->g->SetTitle("");
      this->g->SetLineColor( this->sampleInfo->color );
      this->g->SetLineWidth( 1 );
      this->g->SetMarkerStyle( this->sampleInfo->markerStyle );
      this->g->SetMarkerSize( 1.5 );
      this->g->SetMarkerColor( this->sampleInfo->color );
      this->g->SetFillColorAlpha( kWhite, 0 );
      if( this->histInfo->isFilled )
      {
        this->g->SetLineColorAlpha( kWhite, 0 );
        this->g->SetMarkerColorAlpha( kWhite, 0 );
        this->g->SetFillColorAlpha( this->sampleInfo->color, 1 );
      }

      if( this->histInfo->hasXRange ) {
        this->g->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );
        this->hFrame->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );
      }

      if( this->histInfo->hasYRange ) {
        this->g->GetYaxis()->SetRangeUser( this->histInfo->minY, this->histInfo->maxY );
        this->hFrame->GetYaxis()->SetRangeUser( this->histInfo->minY, this->histInfo->maxY );
      }

      if( this->hasRatio ) {
        SetAxis_TopPad( this->g->GetXaxis(), this->g->GetYaxis(), this->histInfo->titleY );
        SetAxis_TopPad( this->hFrame->GetXaxis(), this->hFrame->GetYaxis(), this->histInfo->titleY );
        if( drawAxis ) {
          SetAxis_TopPad_Axis( this->g->GetXaxis(), this->g->GetYaxis(), this->histInfo->titleX, this->histInfo->titleY );
          SetAxis_TopPad_Axis( this->hFrame->GetXaxis(), this->hFrame->GetYaxis(), this->histInfo->titleX, this->histInfo->titleY );
        }
      }
      else {
        SetAxis_SinglePad( this->g->GetXaxis(), this->g->GetYaxis(), this->histInfo->titleX, this->histInfo->titleY );
        SetAxis_SinglePad( this->hFrame->GetXaxis(), this->hFrame->GetYaxis(), this->histInfo->titleX, this->histInfo->titleY );
      }

      //-- for showers
      if( (this->histInfo->titleX).Contains("shower") ) {
        this->g->GetXaxis()->SetNdivisions(005);
        this->hFrame->GetXaxis()->SetNdivisions(005);
      }
    }

    void SetAttributesRatio(TString ratioTitle, Double_t minRatio, Double_t maxRatio, Bool_t isNoLabel = kFALSE, Int_t nDiv = 505, Bool_t is46 = kFALSE )
    {
      if( this->g_ratio == NULL ) return;

      this->g_ratio->SetTitle( "" );

      // this->g_ratio->GetYaxis()->SetNdivisions(505);
      this->g_ratio->SetLineWidth( 1 );
      this->g_ratio->SetMarkerSize( 1.3 );

      if(!isNoLabel) { // temp
        this->g_ratio->SetLineColor( this->sampleInfo->color );
        this->g_ratio->SetMarkerStyle( this->sampleInfo->markerStyle );
        this->g_ratio->SetMarkerColor( this->sampleInfo->color );
      }

      this->g_ratio->SetFillColorAlpha( kWhite, 0 );
      if( this->histInfo->isFilled )
      {
        this->g_ratio->SetLineColorAlpha( kWhite, 0 );
        this->g_ratio->SetMarkerColorAlpha( kWhite, 0 );
        this->g_ratio->SetFillColorAlpha( this->sampleInfo->color, 1 );
      }

      if( this->histInfo->hasXRange ) {
        this->g_ratio->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );
        this->hFrame_ratio->GetXaxis()->SetRangeUser( this->histInfo->minX, this->histInfo->maxX );
      }

      if(isNoLabel) {
        SetAxis_BottomPad_NoLabel( this->g_ratio->GetXaxis(), this->g_ratio->GetYaxis(), this->histInfo->titleX, ratioTitle, minRatio, maxRatio, nDiv );
        SetAxis_BottomPad_NoLabel( this->hFrame_ratio->GetXaxis(), this->hFrame_ratio->GetYaxis(), this->histInfo->titleX, ratioTitle, minRatio, maxRatio, nDiv );
      }
      else {
        if(is46) {
          SetAxis_BottomPad46( this->g_ratio->GetXaxis(), this->g_ratio->GetYaxis(), this->histInfo->titleX, ratioTitle, minRatio, maxRatio );
          SetAxis_BottomPad46( this->hFrame_ratio->GetXaxis(), this->hFrame_ratio->GetYaxis(), this->histInfo->titleX, ratioTitle, minRatio, maxRatio );
        }
        else {
          SetAxis_BottomPad( this->g_ratio->GetXaxis(), this->g_ratio->GetYaxis(), this->histInfo->titleX, ratioTitle, minRatio, maxRatio );
          SetAxis_BottomPad( this->hFrame_ratio->GetXaxis(), this->hFrame_ratio->GetYaxis(), this->histInfo->titleX, ratioTitle, minRatio, maxRatio );
        }
      }

      //-- for showers
      if( (this->histInfo->titleX).Contains("shower") ) {
        this->g_ratio->GetXaxis()->SetNdivisions(005);
        this->hFrame_ratio->GetXaxis()->SetNdivisions(005);
      }
    }

    TGraphAsymmErrors* MakeRatioGraph(TGraphAsymmErrors *_g_Type1, TGraphAsymmErrors *_g_Type2, TString errorPropagation = "AoverB")
    {
      TGraphAsymmErrors* g_Type1 = (TGraphAsymmErrors*)_g_Type1->Clone();
      TGraphAsymmErrors* g_Type2 = (TGraphAsymmErrors*)_g_Type2->Clone();

      g_ratio = (TGraphAsymmErrors*)g_Type1->Clone();
      g_ratio->Set(0); // --Remove all points (reset) -- //

      Int_t nPoint = g_Type1->GetN();
      Int_t nPoint_2 = g_Type2->GetN();
      if( nPoint != nPoint_2 ) {
        printf("# points is different bewteen two graph... return NULL\n");
        cout << "\tnPoint Num = " << nPoint << endl;
        cout << "\tnPoint Den = " << nPoint_2 << endl;
        return NULL;
      }

      for(Int_t i_p=0; i_p<nPoint; i_p++)
      {
        // cout << i_p << "th Point" << endl;
        //Get Type1 point
        Double_t x_Type1, y_Type1;
        g_Type1->GetPoint(i_p, x_Type1, y_Type1);
        Double_t error_Type1 = this->ReturnLargerValue( g_Type1->GetErrorYhigh(i_p), g_Type1->GetErrorYlow(i_p) );
        // cout << "x_Type1: " << x_Type1 << " y_Type1: " << y_Type1 << " error_Type1: " << error_Type1 << " g_Type1->GetErrorYhigh: " << g_Type1->GetErrorYhigh(i_p) << " g_Type1->GetErrorYlow: " << g_Type1->GetErrorYlow(i_p) << endl;

        //Get Type2 point
        Double_t x_Type2, y_Type2;
        g_Type2->GetPoint(i_p, x_Type2, y_Type2);
        Double_t error_Type2 = this->ReturnLargerValue( g_Type2->GetErrorYhigh(i_p), g_Type2->GetErrorYlow(i_p) );
        // cout << "x_Type2: " << x_Type2 << " y_Type2: " << y_Type2 << " error_Type2: " << error_Type2 << " g_Type2->GetErrorYhigh: " << g_Type2->GetErrorYhigh(i_p) << " g_Type2->GetErrorYlow: " << g_Type2->GetErrorYlow(i_p) << endl;

        Double_t ratio;
        Double_t ratio_error = 999.;
        if(y_Type2 != 0)
        {
          ratio = y_Type1 / y_Type2;
          if(errorPropagation == "AoverB")
            ratio_error = GetUncorrelatedError(y_Type1, error_Type1, y_Type2, error_Type2);
          // else if(errorPropagation == "WeightedBinomial") {
          //   pair<Double_t, Double_t> ErrLU = GetWeightedBinomialErrorLU( y_Type1, error_Type1, (y_Type2-y_Type1), sqrt(y_Type2*y_Type2-error_Type1*error_Type1) );
          //   ratio_error = ErrLU.first;
          // }
          //calculate Scale Factor(Type1/Type2) & error

          // cout << "ratio: " << ratio << " ratio_error: " << ratio_error << endl;
        }
        else
        {
          //if(DEBUG)  cout << "Denominator is 0! ... ratio and its error are set as 0" << endl;
          ratio = 0;
          ratio_error = 0;
        }

        //Set Central value
        g_ratio->SetPoint(i_p, x_Type1, ratio);

        //Set the error
        Double_t error_XLow = g_Type1->GetErrorXlow(i_p);
        Double_t error_Xhigh = g_Type1->GetErrorXhigh(i_p);
        g_ratio->SetPointError(i_p, error_XLow, error_Xhigh, ratio_error, ratio_error);

        // cout << endl;
      }

      return g_ratio;
    }

    TGraphAsymmErrors* MakeDynamicRatioGraph(TGraphAsymmErrors *_g_1, TGraphAsymmErrors *_g_2, TString errorPropagation = "AoverB")
    {
      if(errorPropagation!="AoverB") {
        cout << "WARNING::MakeDynamicRatioGraph errorPropagation!=AoverB" << endl;
        return 0;
      }

      TGraphAsymmErrors* g_1 = (TGraphAsymmErrors*)_g_1->Clone();
      TGraphAsymmErrors* g_2 = (TGraphAsymmErrors*)_g_2->Clone();

      Int_t nPoint_1 = g_1->GetN();
      Int_t nPoint_2 = g_2->GetN();
      Int_t nBin = nPoint_1 > nPoint_2 ? nPoint_1 : nPoint_2;

      TGraphAsymmErrors* g_r = new TGraphAsymmErrors(nBin);
      g_r->Set(0); // --Remove all points (reset) -- //

      Int_t i_fill = 0;

      for(Int_t i_1=0; i_1<nPoint_1; ++i_1) {

        vector<Double_t> vec_x_temp = {};
        vector<Double_t> vec_e_XLo_temp = {};
        vector<Double_t> vec_e_Xhi_temp = {};
        vector<Double_t> vec_r_temp = {};
        vector<Double_t> vec_e_temp = {};

        Double_t x_1, y_1;
        g_1->GetPoint(i_1, x_1, y_1);
        Double_t e_1 = this->ReturnLargerValue( g_1->GetErrorYhigh(i_1), g_1->GetErrorYlow(i_1) );
        Double_t e_XLo_1 = g_1->GetErrorXlow(i_1);
        Double_t e_Xhi_1 = g_1->GetErrorXhigh(i_1);

        // cout << "\n\tType_1: " << i_1 << endl;
        // cout << "\t\t   x_1: " << x_1 << ", " << e_XLo_1 << ", " << e_Xhi_1 << endl;

        for(Int_t i_2=0; i_2<nPoint_2; ++i_2) {

          Double_t x_2, y_2;
          g_2->GetPoint(i_2, x_2, y_2);
          Double_t e_2 = this->ReturnLargerValue( g_2->GetErrorYhigh(i_2), g_2->GetErrorYlow(i_2) );
          Double_t e_XLo_2 = g_2->GetErrorXlow(i_2);
          Double_t e_Xhi_2 = g_2->GetErrorXhigh(i_2);

          Double_t x_temp, e_XLo_temp, e_Xhi_temp, r_temp, e_temp;

          r_temp     = y_2 > 0 ? y_1 / y_2 : 0;
          e_temp     = GetUncorrelatedError( y_1, e_1, y_2, e_2 );

          Int_t iCase = -1;
          //-- case 0
          if( (x_1-e_XLo_1)==(x_2-e_XLo_2) &&
              (x_1+e_Xhi_1)==(x_2+e_Xhi_2) ) {
            x_temp     = x_1;
            e_XLo_temp = e_XLo_1;
            e_Xhi_temp = e_Xhi_1;
            iCase = 0;
          }

          //-- case 1
          else if( (x_1-e_XLo_1)==(x_2-e_XLo_2) &&
                   (x_1+e_Xhi_1)>(x_2+e_Xhi_2) ) {
            x_temp     = x_2;
            e_XLo_temp = e_XLo_2;
            e_Xhi_temp = e_Xhi_2;
            iCase = 1;
          }

          //-- case 2
          else if( (x_1-e_XLo_1)==(x_2-e_XLo_2) &&
                   (x_1+e_Xhi_1)<(x_2+e_Xhi_2) ) {
            x_temp     = x_1;
            e_XLo_temp = e_XLo_1;
            e_Xhi_temp = e_Xhi_1;
            iCase = 2;
          }

          //-- case 3
          else if( (x_1-e_XLo_1)<(x_2-e_XLo_2) &&
                   (x_1+e_Xhi_1)==(x_2+e_Xhi_2) ) {
            x_temp     = x_2;
            e_XLo_temp = e_XLo_2;
            e_Xhi_temp = e_Xhi_2;
            iCase = 3;
          }

          //-- case 4
          else if( (x_1-e_XLo_1)>(x_2-e_XLo_2) &&
                   (x_1+e_Xhi_1)==(x_2+e_Xhi_2) ) {
            x_temp     = x_1;
            e_XLo_temp = e_XLo_1;
            e_Xhi_temp = e_Xhi_1;
            iCase = 4;
          }

          else
            continue;

          if(i_1==nPoint_1-1 && iCase!=0)
            continue;

          // cout << "\t\tType_2: " << i_2 << endl;
          // cout << "\t\t\tiCase: " << iCase << endl;
          // cout << "\t\t\tx_temp: " << x_temp << endl;
          // cout << "\t\t\te_XLo_temp: " << e_XLo_temp << endl;
          // cout << "\t\t\te_Xhi_temp: " << e_Xhi_temp << endl;

          vec_x_temp.push_back(x_temp);
          vec_e_XLo_temp.push_back(e_XLo_temp);
          vec_e_Xhi_temp.push_back(e_Xhi_temp);
          vec_r_temp.push_back(r_temp);
          vec_e_temp.push_back(e_temp);

          if(iCase!=1)
            continue;
        }

        if(vec_r_temp.size() > 0) {
          for(Int_t i_t=0; i_t<(int)vec_r_temp.size(); ++ i_t) {
            g_r->SetPoint(i_fill, vec_x_temp[i_t], vec_r_temp[i_t]);
            g_r->SetPointError(i_fill, vec_e_XLo_temp[i_t], vec_e_Xhi_temp[i_t], vec_e_temp[i_t], vec_e_temp[i_t]);
            i_fill += 1;
          }
        }

      }

      return g_r;
    }

    Double_t ReturnLargerValue(Double_t a, Double_t b)
    {
      if( a > b )
        return a;
      else
        return b;
    }
  };
