//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 23 15:41:01 2020 by ROOT version 6.22/02
// from TChain ntuple/ntuple
//////////////////////////////////////////////////////////

#ifndef MuonHLTNtuple_h
#define MuonHLTNtuple_h

#define ArrSize 50000

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include "vector"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

float TkMuonOfflineEt(double Et, double Eta) {
    vector<float> barrelScalings_ = {0.820128, 1.04124, 0.0};
    vector<float> overlapScalings_ = {0.920897, 1.03712, 0.0};
    vector<float> endcapScalings_ = {0.864715, 1.03215, 0.0};

    if (std::abs(Eta) < 0.9)
        return (barrelScalings_.at(0) + Et * barrelScalings_.at(1) + Et * Et * barrelScalings_.at(2));
    else if (std::abs(Eta) < 1.2)
        return (overlapScalings_.at(0) + Et * overlapScalings_.at(1) + Et * Et * overlapScalings_.at(2));
    else
        return (endcapScalings_.at(0) + Et * endcapScalings_.at(1) + Et * Et * endcapScalings_.at(2));
}

class Object
{
private:
    map<TString, double> vars;
    map<TString, TString> strvars;

public:
    double pt;
    double eta;
    double phi;

    Object()
    {
        pt  = -99999.;
        eta = -99999.;
        phi = -99999.;
    }

    Object( double _pt, double _eta, double _phi ): Object()
    {
      pt            = _pt;
      eta           = _eta;
      phi           = _phi;
    }

    void addVar( TString name, double value )
    {
        vars[name] = value;
    }

    void addStrVar( TString name, TString value )
    {
        strvars[name] = value;
    }

    bool has( TString key )
    {
        return (vars.find(key) != vars.end());
    }

    bool hasstr( TString key )
    {
        return (strvars.find(key) != strvars.end());
    }

    double get( TString key )
    {
        if( vars.find(key) == vars.end() ) {
            cout << key << " does not exist -> return -99999" << endl;
            this->print();
            return -99999.0;
        }
        else {
            return vars[key];
        }
    }

    TString getstr( TString key )
    {
        if( strvars.find(key) == strvars.end() ) {
            cout << key << " does not exist -> return XXXXX" << endl;
            return "XXXXX";
        }
        else {
            return strvars[key];
        }
    }

    void print()
    {
        cout << endl;
        cout << (*this) << endl;
        for( auto& it: vars ) {
            cout << "\t" << it.first << ": " << it.second << endl;
        }
        for( auto& it: strvars ) {
            cout << "\t" << it.first << ": " << it.second << endl;
        }
    }

    Object clone()
    {
      Object out = Object();

      out.pt      = this->pt;
      out.eta     = this->eta;
      out.phi     = this->phi;
      out.vars    = this->vars;
      out.strvars = this->strvars;

      return out;
    }

    double reduceRange(double x) {
        double o2pi = 1. / (2. * M_PI);
        if (std::abs(x) <= double(M_PI))
            return x;
        double n = std::round(x * o2pi);
        return x - n * double(2. * M_PI);
    }

    double deltaPhi(double phi)
    {
        return reduceRange(this->phi - phi);
    }

    double deltaPhi(const Object& other)
    {
        return reduceRange(this->phi - other.phi);
    }

    double deltaR( double eta, double phi )
    {
      double dR = sqrt( (this->eta - eta)*(this->eta - eta)
                      + reduceRange(this->phi - phi)*reduceRange(this->phi - phi) );
      return dR;
    }

    double deltaR( const Object& other )
    {
      double dR = sqrt( (this->eta - other.eta)*(this->eta - other.eta)
                      + reduceRange(this->phi - other.phi)*reduceRange(this->phi - other.phi) );
      return dR;
    }

    bool matched( const Object& other, double dR_match = 0.1, double dpt_match = 1.e9 ) {
      double dR  = deltaR( other );
      double dpt = fabs( this->pt - other.pt ) / this->pt;
      return ( (dR < dR_match) && (dpt < dpt_match) );
    }

    int matched( vector<Object>& objects, vector<int>& map, double dR_match = 0.1, double dpt_match = 1.e9 ) {
      bool found     = false;
      double the_dR  = dR_match;
      int the_i = -1e9;

      unsigned n = objects.size();
      for(unsigned i=0; i<n; ++i) {
        if( map[i] > 0 )  continue;

        double dR  = deltaR( objects.at(i) );
        double dpt = fabs( this->pt - objects.at(i).pt ) / this->pt;
        if( (dR < the_dR) && (dpt < dpt_match) ) {
          found = true;
          the_dR = dR;
          the_i = i;
        }
      }

      if(found) {
        map[the_i] = 1;
      }

      return the_i;
    }

    bool matched( vector<Object>& objects, double dR_match = 0.1, double dpt_match = 1.e9 ) {
      bool found     = false;

      unsigned n = objects.size();
      for(unsigned i=0; i<n; ++i) {
        double dR  = deltaR( objects.at(i) );
        double dpt = fabs( this->pt - objects.at(i).pt ) / this->pt;
        if( (dR < dR_match) && (dpt < dpt_match) ) {
          found = true;
          break;
        }
      }

      return found;
    }

    friend ostream& operator<<(ostream& os, const Object& obj) {
        os << "(" << obj.pt << ", " << obj.eta << ", " << obj.phi << ")";
        return os;
    }

    bool operator==(const Object& other) const {
      return (
        (fabs(this->pt - other.pt) / this->pt < 1.e-3) &&
        (fabs(this->eta - other.eta) < 1.e-3) &&
        (fabs(this->phi - other.phi) < 1.e-3)
      );
    }
};

class MuonHLTNtuple {
public:
    TChain          *fChain;   //!pointer to the analyzed TChain or TChain
    Int_t           fCurrent; //!current Tree number in a TChain
    bool            isIterL3;

    MuonHLTNtuple(TChain *tree=0, vector<TString> branches = {});
    virtual ~MuonHLTNtuple();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TChain *tree);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);

    virtual void     PrintTemplate( TString head );

    vector<TString> getIsoTags();

    vector<Object> get_GenParticles();
    vector<Object> get_L1TkMuons();
    vector<Object> get_L1TTs();
    vector<Object> get_L1Muons();
    vector<Object> get_L2Muons();
    vector<Object> get_L3MuonsNoId();
    vector<Object> get_L3Muons();
    vector<Object> get_iterL3OI();
    vector<Object> get_iterL3IOFromL2();
    vector<Object> get_iterL3FromL2();
    vector<Object> get_iterL3IOFromL1();
    vector<Object> get_iterL3MuonNoID();
    vector<Object> get_iterL3Muon();
    vector<Object> get_hltPhase2L3OI();
    vector<Object> get_hltIter0Phase2L3FromL1TkMuon();
    vector<Object> get_hltIter2Phase2L3FromL1TkMuon();
    vector<Object> get_hltPhase2L3IOFromL1();
    vector<Object> get_hltPhase2L3MuonsNoID();
    vector<Object> get_hltPhase2L3Muons();
    vector<Object> get_tpTo_hltPhase2L3OI();
    vector<Object> get_tpTo_hltIter0Phase2L3FromL1TkMuon();
    vector<Object> get_tpTo_hltIter2Phase2L3FromL1TkMuon();
    vector<Object> get_tpTo_hltPhase2L3IOFromL1();
    vector<Object> get_tpTo_hltPhase2L3MuonsNoID();
    vector<Object> get_tpTo_hltPhase2L3Muons();
    vector<Object> get_hltIterL3OI();
    vector<Object> get_hltIter0IterL3();
    vector<Object> get_hltIter2IterL3();
    vector<Object> get_hltIter0IterL3FromL1Muon();
    vector<Object> get_hltIter2IterL3FromL1Muon();

    //
    vector<Object> get_hltIter2IterL3FromL1MuonTrack();
    //

    vector<Object> get_hltIterL3IOFromL1();
    vector<Object> get_hltIterL3MuonsNoID();
    vector<Object> get_hltIterL3Muons();
    vector<Object> get_tpTo_hltIterL3OI();
    vector<Object> get_tpTo_hltIter0IterL3();
    vector<Object> get_tpTo_hltIter2IterL3();
    vector<Object> get_tpTo_hltIter0IterL3FromL1Muon();
    vector<Object> get_tpTo_hltIter2IterL3FromL1Muon();
    vector<Object> get_tpTo_hltIterL3IOFromL1();
    vector<Object> get_tpTo_hltIterL3MuonsNoID();
    vector<Object> get_tpTo_hltIterL3Muons();
    vector<Object> get_HLTObjects( TString filter = "hltL3fL1TkSingleMu22L3Filtered24Q" );

    vector<TString> branch_names;

    // -- Leaf types -- //
        vector<float>   *trk_pt;
        vector<float>   *trk_eta;
        vector<float>   *trk_phi;
        vector<float>   *trk_d0;
        vector<float>   *trk_z0;
        vector<float>   *trk_rInv;
        vector<float>   *trk_tanL;
        vector<float>   *trk_MVA1;
        vector<float>   *trk_MVA2;
        vector<float>   *trk_MVA3;
        vector<float>   *trk_chi2;
        vector<float>   *trk_bendchi2;
        vector<int>     *trk_nstub;
        vector<int>     *trk_lhits;
        vector<int>     *trk_dhits;
        vector<int>     *trk_seed;
        vector<unsigned int> *trk_phiSector;
        vector<int>     *trk_genuine;
        vector<int>     *trk_loose;
        vector<int>     *trk_unknown;
        vector<int>     *trk_combinatoric;
        vector<int>     *trk_fake;
        vector<int>     *trk_matchtp_pdgid;
        vector<float>   *trk_matchtp_pt;
        vector<float>   *trk_matchtp_eta;
        vector<float>   *trk_matchtp_phi;
        vector<float>   *trk_matchtp_z0;
        vector<float>   *trk_matchtp_dxy;
        vector<float>   *stub_x;
        vector<float>   *stub_y;
        vector<float>   *stub_z;
        vector<int>     *stub_isBarrel;
        vector<int>     *stub_layer;
        vector<float>   *L1TkMu_pt;
        vector<float>   *L1TkMu_eta;
        vector<float>   *L1TkMu_phi;
        vector<float>   *L1TkMu_trkIsol;
        vector<float>   *L1TkMu_trkzVtx;
        vector<float>   *L1TkMu_dR;
        vector<int>     *L1TkMu_nTracksMatched;
        vector<float>   *L1TkMu_trackCurvature;
        vector<unsigned int> *L1TkMu_quality;
        vector<unsigned int> *L1TkMu_pattern;
        vector<unsigned int> *L1TkMu_muonDetector;
        vector<int>     *L1TkMu_TTTpointer;
        vector<float>   *L1TkMu_muRefHwPt;
        vector<int>     *L1TkMu_muRefHwDXY;
        vector<float>   *L1TkMu_muRefHwEta;
        vector<float>   *L1TkMu_muRefHwPhi;
        vector<int>     *L1TkMu_muRefHwSign;
        vector<int>     *L1TkMu_muRefHwSignValid;
        vector<int>     *L1TkMu_muRefHwQual;
        Bool_t          isRealData;
        Int_t           runNum;
        Int_t           lumiBlockNum;
        ULong64_t       eventNum;
        Int_t           nVertex;
        Double_t        bunchID;
        Double_t        instLumi;
        Double_t        dataPU;
        Double_t        dataPURMS;
        Double_t        bunchLumi;
        Double_t        offlineInstLumi;
        Double_t        offlineDataPU;
        Double_t        offlineDataPURMS;
        Double_t        offlineBunchLumi;
        Int_t           truePU;
        Double_t        qScale;
        Double_t        genEventWeight;
        vector<float>   *PU_pT_hats;
        Int_t           nGenParticle;
        Int_t           genParticle_ID[ArrSize];   //[nGenParticle]
        Int_t           genParticle_status[ArrSize];   //[nGenParticle]
        Int_t           genParticle_mother[ArrSize];   //[nGenParticle]
        Double_t        genParticle_pt[ArrSize];   //[nGenParticle]
        Double_t        genParticle_eta[ArrSize];   //[nGenParticle]
        Double_t        genParticle_phi[ArrSize];   //[nGenParticle]
        Double_t        genParticle_px[ArrSize];   //[nGenParticle]
        Double_t        genParticle_py[ArrSize];   //[nGenParticle]
        Double_t        genParticle_pz[ArrSize];   //[nGenParticle]
        Double_t        genParticle_energy[ArrSize];   //[nGenParticle]
        Double_t        genParticle_charge[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isPrompt[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isPromptFinalState[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isTauDecayProduct[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isPromptTauDecayProduct[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isDirectPromptTauDecayProductFinalState[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isHardProcess[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isLastCopy[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isLastCopyBeforeFSR[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isPromptDecayed[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isDecayedLeptonHadron[ArrSize];   //[nGenParticle]
        Int_t           genParticle_fromHardProcessBeforeFSR[ArrSize];   //[nGenParticle]
        Int_t           genParticle_fromHardProcessDecayed[ArrSize];   //[nGenParticle]
        Int_t           genParticle_fromHardProcessFinalState[ArrSize];   //[nGenParticle]
        Int_t           genParticle_isMostlyLikePythia6Status3[ArrSize];   //[nGenParticle]
        vector<string>  *vec_firedTrigger;
        vector<string>  *vec_filterName;
        vector<double>  *vec_HLTObj_pt;
        vector<double>  *vec_HLTObj_eta;
        vector<double>  *vec_HLTObj_phi;
        vector<string>  *vec_myFiredTrigger;
        vector<string>  *vec_myFilterName;
        vector<double>  *vec_myHLTObj_pt;
        vector<double>  *vec_myHLTObj_eta;
        vector<double>  *vec_myHLTObj_phi;
        Int_t           nMuon;
        Double_t        muon_pt[ArrSize];   //[nMuon]
        Double_t        muon_eta[ArrSize];   //[nMuon]
        Double_t        muon_phi[ArrSize];   //[nMuon]
        Double_t        muon_px[ArrSize];   //[nMuon]
        Double_t        muon_py[ArrSize];   //[nMuon]
        Double_t        muon_pz[ArrSize];   //[nMuon]
        Double_t        muon_dB[ArrSize];   //[nMuon]
        Double_t        muon_charge[ArrSize];   //[nMuon]
        Int_t           muon_isGLB[ArrSize];   //[nMuon]
        Int_t           muon_isSTA[ArrSize];   //[nMuon]
        Int_t           muon_isTRK[ArrSize];   //[nMuon]
        Int_t           muon_isPF[ArrSize];   //[nMuon]
        Int_t           muon_isTight[ArrSize];   //[nMuon]
        Int_t           muon_isMedium[ArrSize];   //[nMuon]
        Int_t           muon_isLoose[ArrSize];   //[nMuon]
        Int_t           muon_isHighPt[ArrSize];   //[nMuon]
        Int_t           muon_isHighPtNew[ArrSize];   //[nMuon]
        Int_t           muon_isSoft[ArrSize];   //[nMuon]
        Int_t           muon_isLooseTriggerMuon[ArrSize];   //[nMuon]
        Int_t           muon_isME0Muon[ArrSize];   //[nMuon]
        Int_t           muon_isGEMMuon[ArrSize];   //[nMuon]
        Int_t           muon_isRPCMuon[ArrSize];   //[nMuon]
        Int_t           muon_isGoodMuon_TMOneStationTight[ArrSize];   //[nMuon]
        Double_t        muon_iso03_sumPt[ArrSize];   //[nMuon]
        Double_t        muon_iso03_hadEt[ArrSize];   //[nMuon]
        Double_t        muon_iso03_emEt[ArrSize];   //[nMuon]
        Double_t        muon_PFIso03_charged[ArrSize];   //[nMuon]
        Double_t        muon_PFIso03_neutral[ArrSize];   //[nMuon]
        Double_t        muon_PFIso03_photon[ArrSize];   //[nMuon]
        Double_t        muon_PFIso03_sumPU[ArrSize];   //[nMuon]
        Double_t        muon_PFIso04_charged[ArrSize];   //[nMuon]
        Double_t        muon_PFIso04_neutral[ArrSize];   //[nMuon]
        Double_t        muon_PFIso04_photon[ArrSize];   //[nMuon]
        Double_t        muon_PFIso04_sumPU[ArrSize];   //[nMuon]
        Double_t        muon_PFCluster03_ECAL[ArrSize];   //[nMuon]
        Double_t        muon_PFCluster03_HCAL[ArrSize];   //[nMuon]
        Double_t        muon_PFCluster04_ECAL[ArrSize];   //[nMuon]
        Double_t        muon_PFCluster04_HCAL[ArrSize];   //[nMuon]
        Double_t        muon_normChi2_global[ArrSize];   //[nMuon]
        Int_t           muon_nTrackerHit_global[ArrSize];   //[nMuon]
        Int_t           muon_nTrackerLayer_global[ArrSize];   //[nMuon]
        Int_t           muon_nPixelHit_global[ArrSize];   //[nMuon]
        Int_t           muon_nMuonHit_global[ArrSize];   //[nMuon]
        Double_t        muon_normChi2_inner[ArrSize];   //[nMuon]
        Int_t           muon_nTrackerHit_inner[ArrSize];   //[nMuon]
        Int_t           muon_nTrackerLayer_inner[ArrSize];   //[nMuon]
        Int_t           muon_nPixelHit_inner[ArrSize];   //[nMuon]
        Double_t        muon_pt_tuneP[ArrSize];   //[nMuon]
        Double_t        muon_ptError_tuneP[ArrSize];   //[nMuon]
        Double_t        muon_dxyVTX_best[ArrSize];   //[nMuon]
        Double_t        muon_dzVTX_best[ArrSize];   //[nMuon]
        Int_t           muon_nMatchedStation[ArrSize];   //[nMuon]
        Int_t           muon_nMatchedRPCLayer[ArrSize];   //[nMuon]
        Int_t           muon_stationMask[ArrSize];   //[nMuon]
        Int_t           muon_expectedNnumberOfMatchedStations[ArrSize];   //[nMuon]
        Int_t           nL3Muon;
        Double_t        L3Muon_pt[ArrSize];   //[nL3Muon]
        Double_t        L3Muon_eta[ArrSize];   //[nL3Muon]
        Double_t        L3Muon_phi[ArrSize];   //[nL3Muon]
        Double_t        L3Muon_charge[ArrSize];   //[nL3Muon]
        Double_t        L3Muon_trkPt[ArrSize];   //[nL3Muon]
        Int_t           nL2Muon;
        Double_t        L2Muon_pt[ArrSize];   //[nL2Muon]
        Double_t        L2Muon_eta[ArrSize];   //[nL2Muon]
        Double_t        L2Muon_phi[ArrSize];   //[nL2Muon]
        Double_t        L2Muon_charge[ArrSize];   //[nL2Muon]
        Double_t        L2Muon_trkPt[ArrSize];   //[nL2Muon]
        Int_t           nTkMuon;
        Double_t        TkMuon_pt[ArrSize];   //[nTkMuon]
        Double_t        TkMuon_eta[ArrSize];   //[nTkMuon]
        Double_t        TkMuon_phi[ArrSize];   //[nTkMuon]
        Double_t        TkMuon_charge[ArrSize];   //[nTkMuon]
        Double_t        TkMuon_trkPt[ArrSize];   //[nTkMuon]
        Int_t           nL1Muon;
        Double_t        L1Muon_pt[ArrSize];   //[nL1Muon]
        Double_t        L1Muon_eta[ArrSize];   //[nL1Muon]
        Double_t        L1Muon_phi[ArrSize];   //[nL1Muon]
        Double_t        L1Muon_charge[ArrSize];   //[nL1Muon]
        Double_t        L1Muon_quality[ArrSize];   //[nL1Muon]
        Double_t        L1Muon_etaAtVtx[ArrSize];   //[nL1Muon]
        Double_t        L1Muon_phiAtVtx[ArrSize];   //[nL1Muon]
        Int_t           nIterL3OI;
        Double_t        iterL3OI_inner_pt[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_inner_eta[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_inner_phi[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_inner_charge[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_outer_pt[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_outer_eta[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_outer_phi[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_outer_charge[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_global_pt[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_global_eta[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_global_phi[ArrSize];   //[nIterL3OI]
        Double_t        iterL3OI_global_charge[ArrSize];   //[nIterL3OI]
        Int_t           nIterL3IOFromL2;
        Double_t        iterL3IOFromL2_inner_pt[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_inner_eta[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_inner_phi[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_inner_charge[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_outer_pt[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_outer_eta[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_outer_phi[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_outer_charge[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_global_pt[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_global_eta[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_global_phi[ArrSize];   //[nIterL3IOFromL2]
        Double_t        iterL3IOFromL2_global_charge[ArrSize];   //[nIterL3IOFromL2]
        Int_t           nIterL3IOFromL1;
        Double_t        iterL3IOFromL1_pt[ArrSize];   //[nIterL3IOFromL1]
        Double_t        iterL3IOFromL1_eta[ArrSize];   //[nIterL3IOFromL1]
        Double_t        iterL3IOFromL1_phi[ArrSize];   //[nIterL3IOFromL1]
        Double_t        iterL3IOFromL1_charge[ArrSize];   //[nIterL3IOFromL1]
        Int_t           nIterL3FromL2;
        Double_t        iterL3FromL2_inner_pt[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_inner_eta[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_inner_phi[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_inner_charge[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_outer_pt[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_outer_eta[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_outer_phi[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_outer_charge[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_global_pt[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_global_eta[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_global_phi[ArrSize];   //[nIterL3FromL2]
        Double_t        iterL3FromL2_global_charge[ArrSize];   //[nIterL3FromL2]
        Int_t           nIterL3MuonNoID;
        Double_t        iterL3MuonNoID_pt[ArrSize];   //[nIterL3MuonNoID]
        Double_t        iterL3MuonNoID_innerPt[ArrSize];   //[nIterL3MuonNoID]
        Double_t        iterL3MuonNoID_eta[ArrSize];   //[nIterL3MuonNoID]
        Double_t        iterL3MuonNoID_phi[ArrSize];   //[nIterL3MuonNoID]
        Double_t        iterL3MuonNoID_charge[ArrSize];   //[nIterL3MuonNoID]
        Int_t           iterL3MuonNoID_isGLB[ArrSize];   //[nIterL3MuonNoID]
        Int_t           iterL3MuonNoID_isSTA[ArrSize];   //[nIterL3MuonNoID]
        Int_t           iterL3MuonNoID_isTRK[ArrSize];   //[nIterL3MuonNoID]
        Int_t           nIterL3Muon;
        Double_t        iterL3Muon_pt[ArrSize];   //[nIterL3Muon]
        Double_t        iterL3Muon_innerPt[ArrSize];   //[nIterL3Muon]
        Double_t        iterL3Muon_eta[ArrSize];   //[nIterL3Muon]
        Double_t        iterL3Muon_phi[ArrSize];   //[nIterL3Muon]
        Double_t        iterL3Muon_charge[ArrSize];   //[nIterL3Muon]
        Int_t           iterL3Muon_isGLB[ArrSize];   //[nIterL3Muon]
        Int_t           iterL3Muon_isSTA[ArrSize];   //[nIterL3Muon]
        Int_t           iterL3Muon_isTRK[ArrSize];   //[nIterL3Muon]
        Int_t           nhltIterL3OISeedsFromL2Muons;
        vector<int>     *hltIterL3OISeedsFromL2Muons_dir;
        vector<unsigned int> *hltIterL3OISeedsFromL2Muons_tsos_detId;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_pt;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_pt_val;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_eta;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_phi;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_glob_x;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_glob_y;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_glob_z;
        vector<int>     *hltIterL3OISeedsFromL2Muons_tsos_hasErr;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err0;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err1;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err2;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err3;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err4;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err5;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err6;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err7;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err8;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err9;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err10;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err11;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err12;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err13;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_err14;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_x;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_y;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_dxdz;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_dydz;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_px;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_py;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_pz;
        vector<float>   *hltIterL3OISeedsFromL2Muons_tsos_qbp;
        vector<int>     *hltIterL3OISeedsFromL2Muons_tsos_charge;
        vector<int>     *hltIterL3OISeedsFromL2Muons_iterL3Matched;
        vector<int>     *hltIterL3OISeedsFromL2Muons_iterL3Ref;
        vector<int>     *hltIterL3OISeedsFromL2Muons_tmpL3Ref;
        Int_t           nhltIter0IterL3MuonPixelSeedsFromPixelTracks;
        vector<int>     *hltIter0IterL3MuonPixelSeedsFromPixelTracks_dir;
        vector<unsigned int> *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_detId;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt_val;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_eta;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_phi;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_x;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_y;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_z;
        vector<int>     *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_hasErr;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err0;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err1;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err2;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err3;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err4;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err5;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err6;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err7;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err8;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err9;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err10;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err11;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err12;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err13;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err14;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_x;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_y;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dxdz;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dydz;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_px;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_py;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pz;
        vector<float>   *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_qbp;
        vector<int>     *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_charge;
        vector<int>     *hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Matched;
        vector<int>     *hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Ref;
        vector<int>     *hltIter0IterL3MuonPixelSeedsFromPixelTracks_tmpL3Ref;
        Int_t           nhltIter2IterL3MuonPixelSeeds;
        vector<int>     *hltIter2IterL3MuonPixelSeeds_dir;
        vector<unsigned int> *hltIter2IterL3MuonPixelSeeds_tsos_detId;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_pt;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_pt_val;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_eta;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_phi;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_glob_x;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_glob_y;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_glob_z;
        vector<int>     *hltIter2IterL3MuonPixelSeeds_tsos_hasErr;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err0;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err1;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err2;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err3;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err4;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err5;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err6;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err7;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err8;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err9;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err10;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err11;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err12;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err13;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_err14;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_x;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_y;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_dxdz;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_dydz;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_px;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_py;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_pz;
        vector<float>   *hltIter2IterL3MuonPixelSeeds_tsos_qbp;
        vector<int>     *hltIter2IterL3MuonPixelSeeds_tsos_charge;
        vector<int>     *hltIter2IterL3MuonPixelSeeds_iterL3Matched;
        vector<int>     *hltIter2IterL3MuonPixelSeeds_iterL3Ref;
        vector<int>     *hltIter2IterL3MuonPixelSeeds_tmpL3Ref;
        Int_t           nhltIter3IterL3MuonPixelSeeds;
        vector<int>     *hltIter3IterL3MuonPixelSeeds_dir;
        vector<unsigned int> *hltIter3IterL3MuonPixelSeeds_tsos_detId;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_pt;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_pt_val;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_eta;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_phi;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_glob_x;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_glob_y;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_glob_z;
        vector<int>     *hltIter3IterL3MuonPixelSeeds_tsos_hasErr;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err0;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err1;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err2;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err3;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err4;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err5;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err6;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err7;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err8;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err9;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err10;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err11;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err12;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err13;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_err14;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_x;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_y;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_dxdz;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_dydz;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_px;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_py;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_pz;
        vector<float>   *hltIter3IterL3MuonPixelSeeds_tsos_qbp;
        vector<int>     *hltIter3IterL3MuonPixelSeeds_tsos_charge;
        vector<int>     *hltIter3IterL3MuonPixelSeeds_iterL3Matched;
        vector<int>     *hltIter3IterL3MuonPixelSeeds_iterL3Ref;
        vector<int>     *hltIter3IterL3MuonPixelSeeds_tmpL3Ref;
        Int_t           nhltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks;
        vector<int>     *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_dir;
        vector<unsigned int> *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_detId;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt_val;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_eta;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_phi;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_x;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_y;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_z;
        vector<int>     *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_hasErr;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err0;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err1;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err2;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err3;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err4;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err5;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err6;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err7;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err8;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err9;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err10;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err11;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err12;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err13;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err14;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_x;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_y;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dxdz;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dydz;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_px;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_py;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pz;
        vector<float>   *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_qbp;
        vector<int>     *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_charge;
        vector<int>     *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Matched;
        vector<int>     *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Ref;
        vector<int>     *hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tmpL3Ref;
        Int_t           nhltIter2IterL3FromL1MuonPixelSeeds;
        vector<int>     *hltIter2IterL3FromL1MuonPixelSeeds_dir;
        vector<unsigned int> *hltIter2IterL3FromL1MuonPixelSeeds_tsos_detId;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt_val;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_eta;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_phi;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_x;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_y;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_z;
        vector<int>     *hltIter2IterL3FromL1MuonPixelSeeds_tsos_hasErr;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err0;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err1;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err2;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err3;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err4;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err5;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err6;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err7;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err8;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err9;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err10;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err11;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err12;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err13;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_err14;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_x;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_y;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_dxdz;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_dydz;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_px;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_py;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_pz;
        vector<float>   *hltIter2IterL3FromL1MuonPixelSeeds_tsos_qbp;
        vector<int>     *hltIter2IterL3FromL1MuonPixelSeeds_tsos_charge;
        vector<int>     *hltIter2IterL3FromL1MuonPixelSeeds_iterL3Matched;
        vector<int>     *hltIter2IterL3FromL1MuonPixelSeeds_iterL3Ref;
        vector<int>     *hltIter2IterL3FromL1MuonPixelSeeds_tmpL3Ref;
        Int_t           nhltIter3IterL3FromL1MuonPixelSeeds;
        vector<int>     *hltIter3IterL3FromL1MuonPixelSeeds_dir;
        vector<unsigned int> *hltIter3IterL3FromL1MuonPixelSeeds_tsos_detId;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt_val;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_eta;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_phi;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_x;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_y;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_z;
        vector<int>     *hltIter3IterL3FromL1MuonPixelSeeds_tsos_hasErr;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err0;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err1;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err2;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err3;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err4;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err5;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err6;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err7;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err8;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err9;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err10;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err11;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err12;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err13;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_err14;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_x;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_y;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_dxdz;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_dydz;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_px;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_py;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_pz;
        vector<float>   *hltIter3IterL3FromL1MuonPixelSeeds_tsos_qbp;
        vector<int>     *hltIter3IterL3FromL1MuonPixelSeeds_tsos_charge;
        vector<int>     *hltIter3IterL3FromL1MuonPixelSeeds_iterL3Matched;
        vector<int>     *hltIter3IterL3FromL1MuonPixelSeeds_iterL3Ref;
        vector<int>     *hltIter3IterL3FromL1MuonPixelSeeds_tmpL3Ref;
        Int_t           nhltIterL3OIMuonTrack;
        vector<double>  *hltIterL3OIMuonTrack_pt;
        vector<double>  *hltIterL3OIMuonTrack_ptError;
        vector<double>  *hltIterL3OIMuonTrack_eta;
        vector<double>  *hltIterL3OIMuonTrack_phi;
        vector<int>     *hltIterL3OIMuonTrack_charge;
        vector<int>     *hltIterL3OIMuonTrack_matchedL3;
        vector<int>     *hltIterL3OIMuonTrack_matchedL3NoId;
        vector<float>   *hltIterL3OIMuonTrack_bestMatchTP_charge;
        vector<int>     *hltIterL3OIMuonTrack_bestMatchTP_pdgId;
        vector<double>  *hltIterL3OIMuonTrack_bestMatchTP_energy;
        vector<double>  *hltIterL3OIMuonTrack_bestMatchTP_pt;
        vector<double>  *hltIterL3OIMuonTrack_bestMatchTP_eta;
        vector<double>  *hltIterL3OIMuonTrack_bestMatchTP_phi;
        vector<double>  *hltIterL3OIMuonTrack_bestMatchTP_parentVx;
        vector<double>  *hltIterL3OIMuonTrack_bestMatchTP_parentVy;
        vector<double>  *hltIterL3OIMuonTrack_bestMatchTP_parentVz;
        vector<int>     *hltIterL3OIMuonTrack_bestMatchTP_status;
        vector<int>     *hltIterL3OIMuonTrack_bestMatchTP_numberOfHits;
        vector<int>     *hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIterL3OIMuonTrack_bestMatchTP_sharedFraction;
        vector<int>     *hltIterL3OIMuonTrack_matchedTPsize;
        vector<float>   *hltIterL3OIMuonTrack_mva0;
        vector<float>   *hltIterL3OIMuonTrack_mva1;
        vector<float>   *hltIterL3OIMuonTrack_mva2;
        vector<float>   *hltIterL3OIMuonTrack_mva3;
        Int_t           nhltIter0IterL3MuonTrack;
        vector<double>  *hltIter0IterL3MuonTrack_pt;
        vector<double>  *hltIter0IterL3MuonTrack_ptError;
        vector<double>  *hltIter0IterL3MuonTrack_eta;
        vector<double>  *hltIter0IterL3MuonTrack_phi;
        vector<int>     *hltIter0IterL3MuonTrack_charge;
        vector<int>     *hltIter0IterL3MuonTrack_matchedL3;
        vector<int>     *hltIter0IterL3MuonTrack_matchedL3NoId;
        vector<float>   *hltIter0IterL3MuonTrack_bestMatchTP_charge;
        vector<int>     *hltIter0IterL3MuonTrack_bestMatchTP_pdgId;
        vector<double>  *hltIter0IterL3MuonTrack_bestMatchTP_energy;
        vector<double>  *hltIter0IterL3MuonTrack_bestMatchTP_pt;
        vector<double>  *hltIter0IterL3MuonTrack_bestMatchTP_eta;
        vector<double>  *hltIter0IterL3MuonTrack_bestMatchTP_phi;
        vector<double>  *hltIter0IterL3MuonTrack_bestMatchTP_parentVx;
        vector<double>  *hltIter0IterL3MuonTrack_bestMatchTP_parentVy;
        vector<double>  *hltIter0IterL3MuonTrack_bestMatchTP_parentVz;
        vector<int>     *hltIter0IterL3MuonTrack_bestMatchTP_status;
        vector<int>     *hltIter0IterL3MuonTrack_bestMatchTP_numberOfHits;
        vector<int>     *hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter0IterL3MuonTrack_bestMatchTP_sharedFraction;
        vector<int>     *hltIter0IterL3MuonTrack_matchedTPsize;
        vector<float>   *hltIter0IterL3MuonTrack_mva0;
        vector<float>   *hltIter0IterL3MuonTrack_mva1;
        vector<float>   *hltIter0IterL3MuonTrack_mva2;
        vector<float>   *hltIter0IterL3MuonTrack_mva3;
        Int_t           nhltIter2IterL3MuonTrack;
        vector<double>  *hltIter2IterL3MuonTrack_pt;
        vector<double>  *hltIter2IterL3MuonTrack_ptError;
        vector<double>  *hltIter2IterL3MuonTrack_eta;
        vector<double>  *hltIter2IterL3MuonTrack_phi;
        vector<int>     *hltIter2IterL3MuonTrack_charge;
        vector<int>     *hltIter2IterL3MuonTrack_matchedL3;
        vector<int>     *hltIter2IterL3MuonTrack_matchedL3NoId;
        vector<float>   *hltIter2IterL3MuonTrack_bestMatchTP_charge;
        vector<int>     *hltIter2IterL3MuonTrack_bestMatchTP_pdgId;
        vector<double>  *hltIter2IterL3MuonTrack_bestMatchTP_energy;
        vector<double>  *hltIter2IterL3MuonTrack_bestMatchTP_pt;
        vector<double>  *hltIter2IterL3MuonTrack_bestMatchTP_eta;
        vector<double>  *hltIter2IterL3MuonTrack_bestMatchTP_phi;
        vector<double>  *hltIter2IterL3MuonTrack_bestMatchTP_parentVx;
        vector<double>  *hltIter2IterL3MuonTrack_bestMatchTP_parentVy;
        vector<double>  *hltIter2IterL3MuonTrack_bestMatchTP_parentVz;
        vector<int>     *hltIter2IterL3MuonTrack_bestMatchTP_status;
        vector<int>     *hltIter2IterL3MuonTrack_bestMatchTP_numberOfHits;
        vector<int>     *hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter2IterL3MuonTrack_bestMatchTP_sharedFraction;
        vector<int>     *hltIter2IterL3MuonTrack_matchedTPsize;
        vector<float>   *hltIter2IterL3MuonTrack_mva0;
        vector<float>   *hltIter2IterL3MuonTrack_mva1;
        vector<float>   *hltIter2IterL3MuonTrack_mva2;
        vector<float>   *hltIter2IterL3MuonTrack_mva3;
        Int_t           nhltIter3IterL3MuonTrack;
        vector<double>  *hltIter3IterL3MuonTrack_pt;
        vector<double>  *hltIter3IterL3MuonTrack_ptError;
        vector<double>  *hltIter3IterL3MuonTrack_eta;
        vector<double>  *hltIter3IterL3MuonTrack_phi;
        vector<int>     *hltIter3IterL3MuonTrack_charge;
        vector<int>     *hltIter3IterL3MuonTrack_matchedL3;
        vector<int>     *hltIter3IterL3MuonTrack_matchedL3NoId;
        vector<float>   *hltIter3IterL3MuonTrack_bestMatchTP_charge;
        vector<int>     *hltIter3IterL3MuonTrack_bestMatchTP_pdgId;
        vector<double>  *hltIter3IterL3MuonTrack_bestMatchTP_energy;
        vector<double>  *hltIter3IterL3MuonTrack_bestMatchTP_pt;
        vector<double>  *hltIter3IterL3MuonTrack_bestMatchTP_eta;
        vector<double>  *hltIter3IterL3MuonTrack_bestMatchTP_phi;
        vector<double>  *hltIter3IterL3MuonTrack_bestMatchTP_parentVx;
        vector<double>  *hltIter3IterL3MuonTrack_bestMatchTP_parentVy;
        vector<double>  *hltIter3IterL3MuonTrack_bestMatchTP_parentVz;
        vector<int>     *hltIter3IterL3MuonTrack_bestMatchTP_status;
        vector<int>     *hltIter3IterL3MuonTrack_bestMatchTP_numberOfHits;
        vector<int>     *hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter3IterL3MuonTrack_bestMatchTP_sharedFraction;
        vector<int>     *hltIter3IterL3MuonTrack_matchedTPsize;
        vector<float>   *hltIter3IterL3MuonTrack_mva0;
        vector<float>   *hltIter3IterL3MuonTrack_mva1;
        vector<float>   *hltIter3IterL3MuonTrack_mva2;
        vector<float>   *hltIter3IterL3MuonTrack_mva3;
        Int_t           nhltIter0IterL3FromL1MuonTrack;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_pt;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_ptError;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_eta;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_phi;
        vector<int>     *hltIter0IterL3FromL1MuonTrack_charge;
        vector<int>     *hltIter0IterL3FromL1MuonTrack_matchedL3;
        vector<int>     *hltIter0IterL3FromL1MuonTrack_matchedL3NoId;
        vector<float>   *hltIter0IterL3FromL1MuonTrack_bestMatchTP_charge;
        vector<int>     *hltIter0IterL3FromL1MuonTrack_bestMatchTP_pdgId;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_bestMatchTP_energy;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_bestMatchTP_pt;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_bestMatchTP_eta;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_bestMatchTP_phi;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVx;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVy;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVz;
        vector<int>     *hltIter0IterL3FromL1MuonTrack_bestMatchTP_status;
        vector<int>     *hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfHits;
        vector<int>     *hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter0IterL3FromL1MuonTrack_bestMatchTP_sharedFraction;
        vector<int>     *hltIter0IterL3FromL1MuonTrack_matchedTPsize;
        vector<float>   *hltIter0IterL3FromL1MuonTrack_mva0;
        vector<float>   *hltIter0IterL3FromL1MuonTrack_mva1;
        vector<float>   *hltIter0IterL3FromL1MuonTrack_mva2;
        vector<float>   *hltIter0IterL3FromL1MuonTrack_mva3;
        Int_t           nhltIter2IterL3FromL1MuonTrack;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_pt;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_ptError;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_eta;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_phi;
        vector<int>     *hltIter2IterL3FromL1MuonTrack_charge;
        vector<int>     *hltIter2IterL3FromL1MuonTrack_matchedL3;
        vector<int>     *hltIter2IterL3FromL1MuonTrack_matchedL3NoId;
        vector<float>   *hltIter2IterL3FromL1MuonTrack_bestMatchTP_charge;
        vector<int>     *hltIter2IterL3FromL1MuonTrack_bestMatchTP_pdgId;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_bestMatchTP_energy;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_bestMatchTP_pt;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_bestMatchTP_eta;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_bestMatchTP_phi;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVx;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVy;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVz;
        vector<int>     *hltIter2IterL3FromL1MuonTrack_bestMatchTP_status;
        vector<int>     *hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfHits;
        vector<int>     *hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter2IterL3FromL1MuonTrack_bestMatchTP_sharedFraction;
        vector<int>     *hltIter2IterL3FromL1MuonTrack_matchedTPsize;
        vector<float>   *hltIter2IterL3FromL1MuonTrack_mva0;
        vector<float>   *hltIter2IterL3FromL1MuonTrack_mva1;
        vector<float>   *hltIter2IterL3FromL1MuonTrack_mva2;
        vector<float>   *hltIter2IterL3FromL1MuonTrack_mva3;
        Int_t           nhltIter3IterL3FromL1MuonTrack;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_pt;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_ptError;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_eta;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_phi;
        vector<int>     *hltIter3IterL3FromL1MuonTrack_charge;
        vector<int>     *hltIter3IterL3FromL1MuonTrack_matchedL3;
        vector<int>     *hltIter3IterL3FromL1MuonTrack_matchedL3NoId;
        vector<float>   *hltIter3IterL3FromL1MuonTrack_bestMatchTP_charge;
        vector<int>     *hltIter3IterL3FromL1MuonTrack_bestMatchTP_pdgId;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_bestMatchTP_energy;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_bestMatchTP_pt;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_bestMatchTP_eta;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_bestMatchTP_phi;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVx;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVy;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVz;
        vector<int>     *hltIter3IterL3FromL1MuonTrack_bestMatchTP_status;
        vector<int>     *hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfHits;
        vector<int>     *hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter3IterL3FromL1MuonTrack_bestMatchTP_sharedFraction;
        vector<int>     *hltIter3IterL3FromL1MuonTrack_matchedTPsize;
        vector<float>   *hltIter3IterL3FromL1MuonTrack_mva0;
        vector<float>   *hltIter3IterL3FromL1MuonTrack_mva1;
        vector<float>   *hltIter3IterL3FromL1MuonTrack_mva2;
        vector<float>   *hltIter3IterL3FromL1MuonTrack_mva3;

        Int_t           nL3MuonsNoId;
        vector<double>  *L3MuonsNoId_pt;
        vector<double>  *L3MuonsNoId_inner_pt;
        vector<double>  *L3MuonsNoId_inner_ptError;
        vector<double>  *L3MuonsNoId_eta;
        vector<double>  *L3MuonsNoId_phi;
        vector<double>  *L3MuonsNoId_charge;
        vector<double>  *L3MuonsNoId_isGlobalMuon;
        vector<double>  *L3MuonsNoId_isStandAloneMuon;
        vector<double>  *L3MuonsNoId_isTrackerMuon;
        vector<double>  *L3MuonsNoId_isLooseTriggerMuon;
        vector<double>  *L3MuonsNoId_isME0Muon;
        vector<double>  *L3MuonsNoId_isGEMMuon;
        vector<double>  *L3MuonsNoId_isRPCMuon;
        vector<double>  *L3MuonsNoId_isGoodMuon_TMOneStationTight;
        vector<double>  *L3MuonsNoId_numberOfMatchedStations;
        vector<double>  *L3MuonsNoId_numberOfMatchedRPCLayers;
        vector<double>  *L3MuonsNoId_expectedNnumberOfMatchedStations;
        vector<double>  *L3MuonsNoId_inner_normalizedChi2;
        vector<double>  *L3MuonsNoId_inner_numberOfValidTrackerHits;
        vector<double>  *L3MuonsNoId_inner_trackerLayersWithMeasurement;
        vector<double>  *L3MuonsNoId_inner_numberOfValidPixelHits;
        vector<double>  *L3MuonsNoId_inner_dz_l1vtx;
        Int_t           nL3Muons;
        vector<double>  *L3Muons_pt;
        vector<double>  *L3Muons_inner_pt;
        vector<double>  *L3Muons_inner_ptError;
        vector<double>  *L3Muons_eta;
        vector<double>  *L3Muons_phi;
        vector<double>  *L3Muons_charge;
        vector<double>  *L3Muons_isGlobalMuon;
        vector<double>  *L3Muons_isStandAloneMuon;
        vector<double>  *L3Muons_isTrackerMuon;
        vector<double>  *L3Muons_isLooseTriggerMuon;
        vector<double>  *L3Muons_isME0Muon;
        vector<double>  *L3Muons_isGEMMuon;
        vector<double>  *L3Muons_isRPCMuon;
        vector<double>  *L3Muons_isGoodMuon_TMOneStationTight;
        vector<double>  *L3Muons_numberOfMatchedStations;
        vector<double>  *L3Muons_numberOfMatchedRPCLayers;
        vector<double>  *L3Muons_expectedNnumberOfMatchedStations;
        vector<double>  *L3Muons_inner_normalizedChi2;
        vector<double>  *L3Muons_inner_numberOfValidTrackerHits;
        vector<double>  *L3Muons_inner_trackerLayersWithMeasurement;
        vector<double>  *L3Muons_inner_numberOfValidPixelHits;
        vector<double>  *L3Muons_inner_dz_l1vtx;

        Int_t           nTP;
        vector<float>   *TP_charge;
        vector<int>     *TP_pdgId;
        vector<double>  *TP_energy;
        vector<double>  *TP_pt;
        vector<double>  *TP_eta;
        vector<double>  *TP_phi;
        vector<double>  *TP_parentVx;
        vector<double>  *TP_parentVy;
        vector<double>  *TP_parentVz;
        vector<int>     *TP_status;
        vector<int>     *TP_numberOfHits;
        vector<int>     *TP_numberOfTrackerHits;
        vector<int>     *TP_numberOfTrackerLayers;
        vector<float>   *TP_gen_charge;
        vector<int>     *TP_gen_pdgId;
        vector<double>  *TP_gen_pt;
        vector<double>  *TP_gen_eta;
        vector<double>  *TP_gen_phi;
        vector<double>  *TP_bestMatchTrk_pt;
        vector<double>  *TP_bestMatchTrk_eta;
        vector<double>  *TP_bestMatchTrk_phi;
        vector<int>     *TP_bestMatchTrk_charge;
        vector<double>  *TP_bestMatchTrk_quality;
        vector<int>     *TP_bestMatchTrk_NValidHits;
        Int_t           nhltIterL3MuonTrimmedPixelVertices;
        vector<int>     *hltIterL3MuonTrimmedPixelVertices_isValid;
        vector<double>  *hltIterL3MuonTrimmedPixelVertices_chi2;
        vector<double>  *hltIterL3MuonTrimmedPixelVertices_ndof;
        vector<double>  *hltIterL3MuonTrimmedPixelVertices_nTracks;
        vector<double>  *hltIterL3MuonTrimmedPixelVertices_x;
        vector<double>  *hltIterL3MuonTrimmedPixelVertices_xerr;
        vector<double>  *hltIterL3MuonTrimmedPixelVertices_y;
        vector<double>  *hltIterL3MuonTrimmedPixelVertices_yerr;
        vector<double>  *hltIterL3MuonTrimmedPixelVertices_z;
        vector<double>  *hltIterL3MuonTrimmedPixelVertices_zerr;
        Int_t           nhltIterL3FromL1MuonTrimmedPixelVertices;
        vector<int>     *hltIterL3FromL1MuonTrimmedPixelVertices_isValid;
        vector<double>  *hltIterL3FromL1MuonTrimmedPixelVertices_chi2;
        vector<double>  *hltIterL3FromL1MuonTrimmedPixelVertices_ndof;
        vector<double>  *hltIterL3FromL1MuonTrimmedPixelVertices_nTracks;
        vector<double>  *hltIterL3FromL1MuonTrimmedPixelVertices_x;
        vector<double>  *hltIterL3FromL1MuonTrimmedPixelVertices_xerr;
        vector<double>  *hltIterL3FromL1MuonTrimmedPixelVertices_y;
        vector<double>  *hltIterL3FromL1MuonTrimmedPixelVertices_yerr;
        vector<double>  *hltIterL3FromL1MuonTrimmedPixelVertices_z;
        vector<double>  *hltIterL3FromL1MuonTrimmedPixelVertices_zerr;

        Int_t           nhltPhase2L3OI;
        vector<double>  *hltPhase2L3OI_pt;
        vector<double>  *hltPhase2L3OI_ptError;
        vector<double>  *hltPhase2L3OI_eta;
        vector<double>  *hltPhase2L3OI_phi;
        vector<int>     *hltPhase2L3OI_charge;
        vector<int>     *hltPhase2L3OI_matchedL3;
        vector<int>     *hltPhase2L3OI_matchedL3NoId;
        vector<float>   *hltPhase2L3OI_bestMatchTP_charge;
        vector<int>     *hltPhase2L3OI_bestMatchTP_pdgId;
        vector<double>  *hltPhase2L3OI_bestMatchTP_energy;
        vector<double>  *hltPhase2L3OI_bestMatchTP_pt;
        vector<double>  *hltPhase2L3OI_bestMatchTP_eta;
        vector<double>  *hltPhase2L3OI_bestMatchTP_phi;
        vector<double>  *hltPhase2L3OI_bestMatchTP_parentVx;
        vector<double>  *hltPhase2L3OI_bestMatchTP_parentVy;
        vector<double>  *hltPhase2L3OI_bestMatchTP_parentVz;
        vector<int>     *hltPhase2L3OI_bestMatchTP_status;
        vector<int>     *hltPhase2L3OI_bestMatchTP_numberOfHits;
        vector<int>     *hltPhase2L3OI_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltPhase2L3OI_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltPhase2L3OI_bestMatchTP_sharedFraction;
        vector<int>     *hltPhase2L3OI_matchedTPsize;
        vector<float>   *hltPhase2L3OI_mva0;
        vector<float>   *hltPhase2L3OI_mva1;
        vector<float>   *hltPhase2L3OI_mva2;
        vector<float>   *hltPhase2L3OI_mva3;
        Int_t           ntpTo_hltPhase2L3OI;
        vector<float>   *tpTo_hltPhase2L3OI_charge;
        vector<int>     *tpTo_hltPhase2L3OI_pdgId;
        vector<double>  *tpTo_hltPhase2L3OI_energy;
        vector<double>  *tpTo_hltPhase2L3OI_pt;
        vector<double>  *tpTo_hltPhase2L3OI_eta;
        vector<double>  *tpTo_hltPhase2L3OI_phi;
        vector<double>  *tpTo_hltPhase2L3OI_parentVx;
        vector<double>  *tpTo_hltPhase2L3OI_parentVy;
        vector<double>  *tpTo_hltPhase2L3OI_parentVz;
        vector<int>     *tpTo_hltPhase2L3OI_status;
        vector<int>     *tpTo_hltPhase2L3OI_numberOfHits;
        vector<int>     *tpTo_hltPhase2L3OI_numberOfTrackerHits;
        vector<int>     *tpTo_hltPhase2L3OI_numberOfTrackerLayers;
        vector<float>   *tpTo_hltPhase2L3OI_gen_charge;
        vector<int>     *tpTo_hltPhase2L3OI_gen_pdgId;
        vector<double>  *tpTo_hltPhase2L3OI_gen_pt;
        vector<double>  *tpTo_hltPhase2L3OI_gen_eta;
        vector<double>  *tpTo_hltPhase2L3OI_gen_phi;
        vector<double>  *tpTo_hltPhase2L3OI_bestMatchTrk_pt;
        vector<double>  *tpTo_hltPhase2L3OI_bestMatchTrk_eta;
        vector<double>  *tpTo_hltPhase2L3OI_bestMatchTrk_phi;
        vector<int>     *tpTo_hltPhase2L3OI_bestMatchTrk_charge;
        vector<double>  *tpTo_hltPhase2L3OI_bestMatchTrk_quality;
        vector<int>     *tpTo_hltPhase2L3OI_bestMatchTrk_NValidHits;

        Int_t           nhltIter0Phase2L3FromL1TkMuon;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_pt;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_ptError;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_eta;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_phi;
        vector<int>     *hltIter0Phase2L3FromL1TkMuon_charge;
        vector<int>     *hltIter0Phase2L3FromL1TkMuon_matchedL3;
        vector<int>     *hltIter0Phase2L3FromL1TkMuon_matchedL3NoId;
        vector<float>   *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_charge;
        vector<int>     *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pdgId;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_energy;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pt;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_eta;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_phi;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVx;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVy;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVz;
        vector<int>     *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_status;
        vector<int>     *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits;
        vector<int>     *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter0Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction;
        vector<int>     *hltIter0Phase2L3FromL1TkMuon_matchedTPsize;
        vector<float>   *hltIter0Phase2L3FromL1TkMuon_mva0;
        vector<float>   *hltIter0Phase2L3FromL1TkMuon_mva1;
        vector<float>   *hltIter0Phase2L3FromL1TkMuon_mva2;
        vector<float>   *hltIter0Phase2L3FromL1TkMuon_mva3;
        Int_t           ntpTo_hltIter0Phase2L3FromL1TkMuon;
        vector<float>   *tpTo_hltIter0Phase2L3FromL1TkMuon_charge;
        vector<int>     *tpTo_hltIter0Phase2L3FromL1TkMuon_pdgId;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_energy;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_pt;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_eta;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_phi;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_parentVx;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_parentVy;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_parentVz;
        vector<int>     *tpTo_hltIter0Phase2L3FromL1TkMuon_status;
        vector<int>     *tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfHits;
        vector<int>     *tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerHits;
        vector<int>     *tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIter0Phase2L3FromL1TkMuon_gen_charge;
        vector<int>     *tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pdgId;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pt;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_gen_eta;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_gen_phi;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits;
        Int_t           nhltIter2Phase2L3FromL1TkMuon;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_pt;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_ptError;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_eta;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_phi;
        vector<int>     *hltIter2Phase2L3FromL1TkMuon_charge;
        vector<int>     *hltIter2Phase2L3FromL1TkMuon_matchedL3;
        vector<int>     *hltIter2Phase2L3FromL1TkMuon_matchedL3NoId;
        vector<float>   *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_charge;
        vector<int>     *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pdgId;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_energy;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pt;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_eta;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_phi;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVx;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVy;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVz;
        vector<int>     *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_status;
        vector<int>     *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits;
        vector<int>     *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter2Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction;
        vector<int>     *hltIter2Phase2L3FromL1TkMuon_matchedTPsize;
        vector<float>   *hltIter2Phase2L3FromL1TkMuon_mva0;
        vector<float>   *hltIter2Phase2L3FromL1TkMuon_mva1;
        vector<float>   *hltIter2Phase2L3FromL1TkMuon_mva2;
        vector<float>   *hltIter2Phase2L3FromL1TkMuon_mva3;
        Int_t           ntpTo_hltIter2Phase2L3FromL1TkMuon;
        vector<float>   *tpTo_hltIter2Phase2L3FromL1TkMuon_charge;
        vector<int>     *tpTo_hltIter2Phase2L3FromL1TkMuon_pdgId;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_energy;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_pt;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_eta;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_phi;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_parentVx;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_parentVy;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_parentVz;
        vector<int>     *tpTo_hltIter2Phase2L3FromL1TkMuon_status;
        vector<int>     *tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfHits;
        vector<int>     *tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerHits;
        vector<int>     *tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIter2Phase2L3FromL1TkMuon_gen_charge;
        vector<int>     *tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pdgId;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pt;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_gen_eta;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_gen_phi;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits;
        Int_t           nhltPhase2L3IOFromL1;
        vector<double>  *hltPhase2L3IOFromL1_pt;
        vector<double>  *hltPhase2L3IOFromL1_ptError;
        vector<double>  *hltPhase2L3IOFromL1_eta;
        vector<double>  *hltPhase2L3IOFromL1_phi;
        vector<int>     *hltPhase2L3IOFromL1_charge;
        vector<int>     *hltPhase2L3IOFromL1_matchedL3;
        vector<int>     *hltPhase2L3IOFromL1_matchedL3NoId;
        vector<float>   *hltPhase2L3IOFromL1_bestMatchTP_charge;
        vector<int>     *hltPhase2L3IOFromL1_bestMatchTP_pdgId;
        vector<double>  *hltPhase2L3IOFromL1_bestMatchTP_energy;
        vector<double>  *hltPhase2L3IOFromL1_bestMatchTP_pt;
        vector<double>  *hltPhase2L3IOFromL1_bestMatchTP_eta;
        vector<double>  *hltPhase2L3IOFromL1_bestMatchTP_phi;
        vector<double>  *hltPhase2L3IOFromL1_bestMatchTP_parentVx;
        vector<double>  *hltPhase2L3IOFromL1_bestMatchTP_parentVy;
        vector<double>  *hltPhase2L3IOFromL1_bestMatchTP_parentVz;
        vector<int>     *hltPhase2L3IOFromL1_bestMatchTP_status;
        vector<int>     *hltPhase2L3IOFromL1_bestMatchTP_numberOfHits;
        vector<int>     *hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltPhase2L3IOFromL1_bestMatchTP_sharedFraction;
        vector<int>     *hltPhase2L3IOFromL1_matchedTPsize;
        vector<float>   *hltPhase2L3IOFromL1_mva0;
        vector<float>   *hltPhase2L3IOFromL1_mva1;
        vector<float>   *hltPhase2L3IOFromL1_mva2;
        vector<float>   *hltPhase2L3IOFromL1_mva3;
        Int_t           ntpTo_hltPhase2L3IOFromL1;
        vector<float>   *tpTo_hltPhase2L3IOFromL1_charge;
        vector<int>     *tpTo_hltPhase2L3IOFromL1_pdgId;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_energy;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_pt;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_eta;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_phi;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_parentVx;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_parentVy;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_parentVz;
        vector<int>     *tpTo_hltPhase2L3IOFromL1_status;
        vector<int>     *tpTo_hltPhase2L3IOFromL1_numberOfHits;
        vector<int>     *tpTo_hltPhase2L3IOFromL1_numberOfTrackerHits;
        vector<int>     *tpTo_hltPhase2L3IOFromL1_numberOfTrackerLayers;
        vector<float>   *tpTo_hltPhase2L3IOFromL1_gen_charge;
        vector<int>     *tpTo_hltPhase2L3IOFromL1_gen_pdgId;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_gen_pt;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_gen_eta;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_gen_phi;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_bestMatchTrk_pt;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_bestMatchTrk_eta;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_bestMatchTrk_phi;
        vector<int>     *tpTo_hltPhase2L3IOFromL1_bestMatchTrk_charge;
        vector<double>  *tpTo_hltPhase2L3IOFromL1_bestMatchTrk_quality;
        vector<int>     *tpTo_hltPhase2L3IOFromL1_bestMatchTrk_NValidHits;
        Int_t           nhltPhase2L3MuonsNoID;
        vector<double>  *hltPhase2L3MuonsNoID_pt;
        vector<double>  *hltPhase2L3MuonsNoID_ptError;
        vector<double>  *hltPhase2L3MuonsNoID_eta;
        vector<double>  *hltPhase2L3MuonsNoID_phi;
        vector<int>     *hltPhase2L3MuonsNoID_charge;
        vector<int>     *hltPhase2L3MuonsNoID_matchedL3;
        vector<int>     *hltPhase2L3MuonsNoID_matchedL3NoId;
        vector<float>   *hltPhase2L3MuonsNoID_bestMatchTP_charge;
        vector<int>     *hltPhase2L3MuonsNoID_bestMatchTP_pdgId;
        vector<double>  *hltPhase2L3MuonsNoID_bestMatchTP_energy;
        vector<double>  *hltPhase2L3MuonsNoID_bestMatchTP_pt;
        vector<double>  *hltPhase2L3MuonsNoID_bestMatchTP_eta;
        vector<double>  *hltPhase2L3MuonsNoID_bestMatchTP_phi;
        vector<double>  *hltPhase2L3MuonsNoID_bestMatchTP_parentVx;
        vector<double>  *hltPhase2L3MuonsNoID_bestMatchTP_parentVy;
        vector<double>  *hltPhase2L3MuonsNoID_bestMatchTP_parentVz;
        vector<int>     *hltPhase2L3MuonsNoID_bestMatchTP_status;
        vector<int>     *hltPhase2L3MuonsNoID_bestMatchTP_numberOfHits;
        vector<int>     *hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltPhase2L3MuonsNoID_bestMatchTP_sharedFraction;
        vector<int>     *hltPhase2L3MuonsNoID_matchedTPsize;
        vector<float>   *hltPhase2L3MuonsNoID_mva0;
        vector<float>   *hltPhase2L3MuonsNoID_mva1;
        vector<float>   *hltPhase2L3MuonsNoID_mva2;
        vector<float>   *hltPhase2L3MuonsNoID_mva3;
        Int_t           ntpTo_hltPhase2L3MuonsNoID;
        vector<float>   *tpTo_hltPhase2L3MuonsNoID_charge;
        vector<int>     *tpTo_hltPhase2L3MuonsNoID_pdgId;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_energy;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_pt;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_eta;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_phi;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_parentVx;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_parentVy;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_parentVz;
        vector<int>     *tpTo_hltPhase2L3MuonsNoID_status;
        vector<int>     *tpTo_hltPhase2L3MuonsNoID_numberOfHits;
        vector<int>     *tpTo_hltPhase2L3MuonsNoID_numberOfTrackerHits;
        vector<int>     *tpTo_hltPhase2L3MuonsNoID_numberOfTrackerLayers;
        vector<float>   *tpTo_hltPhase2L3MuonsNoID_gen_charge;
        vector<int>     *tpTo_hltPhase2L3MuonsNoID_gen_pdgId;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_gen_pt;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_gen_eta;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_gen_phi;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_pt;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_eta;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_phi;
        vector<int>     *tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_charge;
        vector<double>  *tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_quality;
        vector<int>     *tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_NValidHits;
        Int_t           nhltPhase2L3Muons;
        vector<double>  *hltPhase2L3Muons_pt;
        vector<double>  *hltPhase2L3Muons_ptError;
        vector<double>  *hltPhase2L3Muons_eta;
        vector<double>  *hltPhase2L3Muons_phi;
        vector<int>     *hltPhase2L3Muons_charge;
        vector<int>     *hltPhase2L3Muons_matchedL3;
        vector<int>     *hltPhase2L3Muons_matchedL3NoId;
        vector<float>   *hltPhase2L3Muons_bestMatchTP_charge;
        vector<int>     *hltPhase2L3Muons_bestMatchTP_pdgId;
        vector<double>  *hltPhase2L3Muons_bestMatchTP_energy;
        vector<double>  *hltPhase2L3Muons_bestMatchTP_pt;
        vector<double>  *hltPhase2L3Muons_bestMatchTP_eta;
        vector<double>  *hltPhase2L3Muons_bestMatchTP_phi;
        vector<double>  *hltPhase2L3Muons_bestMatchTP_parentVx;
        vector<double>  *hltPhase2L3Muons_bestMatchTP_parentVy;
        vector<double>  *hltPhase2L3Muons_bestMatchTP_parentVz;
        vector<int>     *hltPhase2L3Muons_bestMatchTP_status;
        vector<int>     *hltPhase2L3Muons_bestMatchTP_numberOfHits;
        vector<int>     *hltPhase2L3Muons_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltPhase2L3Muons_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltPhase2L3Muons_bestMatchTP_sharedFraction;
        vector<int>     *hltPhase2L3Muons_matchedTPsize;
        vector<float>   *hltPhase2L3Muons_mva0;
        vector<float>   *hltPhase2L3Muons_mva1;
        vector<float>   *hltPhase2L3Muons_mva2;
        vector<float>   *hltPhase2L3Muons_mva3;

        // -- Isolation
        vector<float>   *hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0;
        vector<float>   *hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p000;
        vector<float>   *hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p000;
        vector<float>   *hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p030;
        vector<float>   *hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p030;
        vector<float>   *hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p050;
        vector<float>   *hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p050;
        vector<float>   *hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p00minEEM0p00minEHad0p00;
        vector<float>   *hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00;
        vector<float>   *hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p04minEEM0p00minEHad0p00;
        vector<float>   *hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p00minEEM0p00minEHad0p00;
        vector<float>   *hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p02minEEM0p00minEHad0p00;
        vector<float>   *hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p04minEEM0p00minEHad0p00;
        vector<float>   *hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p00minEEM0p00minEHad0p00;
        vector<float>   *hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p02minEEM0p00minEHad0p00;
        vector<float>   *hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p04minEEM0p00minEHad0p00;

        Int_t           ntpTo_hltPhase2L3Muons;
        vector<float>   *tpTo_hltPhase2L3Muons_charge;
        vector<int>     *tpTo_hltPhase2L3Muons_pdgId;
        vector<double>  *tpTo_hltPhase2L3Muons_energy;
        vector<double>  *tpTo_hltPhase2L3Muons_pt;
        vector<double>  *tpTo_hltPhase2L3Muons_eta;
        vector<double>  *tpTo_hltPhase2L3Muons_phi;
        vector<double>  *tpTo_hltPhase2L3Muons_parentVx;
        vector<double>  *tpTo_hltPhase2L3Muons_parentVy;
        vector<double>  *tpTo_hltPhase2L3Muons_parentVz;
        vector<int>     *tpTo_hltPhase2L3Muons_status;
        vector<int>     *tpTo_hltPhase2L3Muons_numberOfHits;
        vector<int>     *tpTo_hltPhase2L3Muons_numberOfTrackerHits;
        vector<int>     *tpTo_hltPhase2L3Muons_numberOfTrackerLayers;
        vector<float>   *tpTo_hltPhase2L3Muons_gen_charge;
        vector<int>     *tpTo_hltPhase2L3Muons_gen_pdgId;
        vector<double>  *tpTo_hltPhase2L3Muons_gen_pt;
        vector<double>  *tpTo_hltPhase2L3Muons_gen_eta;
        vector<double>  *tpTo_hltPhase2L3Muons_gen_phi;
        vector<double>  *tpTo_hltPhase2L3Muons_bestMatchTrk_pt;
        vector<double>  *tpTo_hltPhase2L3Muons_bestMatchTrk_eta;
        vector<double>  *tpTo_hltPhase2L3Muons_bestMatchTrk_phi;
        vector<int>     *tpTo_hltPhase2L3Muons_bestMatchTrk_charge;
        vector<double>  *tpTo_hltPhase2L3Muons_bestMatchTrk_quality;
        vector<int>     *tpTo_hltPhase2L3Muons_bestMatchTrk_NValidHits;

        Int_t           nhltIterL3OI;
        vector<double>  *hltIterL3OI_pt;
        vector<double>  *hltIterL3OI_ptError;
        vector<double>  *hltIterL3OI_eta;
        vector<double>  *hltIterL3OI_phi;
        vector<int>     *hltIterL3OI_charge;
        vector<int>     *hltIterL3OI_matchedL3;
        vector<int>     *hltIterL3OI_matchedL3NoId;
        vector<float>   *hltIterL3OI_bestMatchTP_charge;
        vector<int>     *hltIterL3OI_bestMatchTP_pdgId;
        vector<double>  *hltIterL3OI_bestMatchTP_energy;
        vector<double>  *hltIterL3OI_bestMatchTP_pt;
        vector<double>  *hltIterL3OI_bestMatchTP_eta;
        vector<double>  *hltIterL3OI_bestMatchTP_phi;
        vector<double>  *hltIterL3OI_bestMatchTP_parentVx;
        vector<double>  *hltIterL3OI_bestMatchTP_parentVy;
        vector<double>  *hltIterL3OI_bestMatchTP_parentVz;
        vector<int>     *hltIterL3OI_bestMatchTP_status;
        vector<int>     *hltIterL3OI_bestMatchTP_numberOfHits;
        vector<int>     *hltIterL3OI_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIterL3OI_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIterL3OI_bestMatchTP_sharedFraction;
        vector<int>     *hltIterL3OI_matchedTPsize;
        vector<float>   *hltIterL3OI_mva0;
        vector<float>   *hltIterL3OI_mva1;
        vector<float>   *hltIterL3OI_mva2;
        vector<float>   *hltIterL3OI_mva3;
        Int_t           ntpTo_hltIterL3OI;
        vector<float>   *tpTo_hltIterL3OI_charge;
        vector<int>     *tpTo_hltIterL3OI_pdgId;
        vector<double>  *tpTo_hltIterL3OI_energy;
        vector<double>  *tpTo_hltIterL3OI_pt;
        vector<double>  *tpTo_hltIterL3OI_eta;
        vector<double>  *tpTo_hltIterL3OI_phi;
        vector<double>  *tpTo_hltIterL3OI_parentVx;
        vector<double>  *tpTo_hltIterL3OI_parentVy;
        vector<double>  *tpTo_hltIterL3OI_parentVz;
        vector<int>     *tpTo_hltIterL3OI_status;
        vector<int>     *tpTo_hltIterL3OI_numberOfHits;
        vector<int>     *tpTo_hltIterL3OI_numberOfTrackerHits;
        vector<int>     *tpTo_hltIterL3OI_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIterL3OI_gen_charge;
        vector<int>     *tpTo_hltIterL3OI_gen_pdgId;
        vector<double>  *tpTo_hltIterL3OI_gen_pt;
        vector<double>  *tpTo_hltIterL3OI_gen_eta;
        vector<double>  *tpTo_hltIterL3OI_gen_phi;
        vector<double>  *tpTo_hltIterL3OI_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIterL3OI_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIterL3OI_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIterL3OI_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIterL3OI_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIterL3OI_bestMatchTrk_NValidHits;
        Int_t           nhltIter0IterL3;
        vector<double>  *hltIter0IterL3_pt;
        vector<double>  *hltIter0IterL3_ptError;
        vector<double>  *hltIter0IterL3_eta;
        vector<double>  *hltIter0IterL3_phi;
        vector<int>     *hltIter0IterL3_charge;
        vector<int>     *hltIter0IterL3_matchedL3;
        vector<int>     *hltIter0IterL3_matchedL3NoId;
        vector<float>   *hltIter0IterL3_bestMatchTP_charge;
        vector<int>     *hltIter0IterL3_bestMatchTP_pdgId;
        vector<double>  *hltIter0IterL3_bestMatchTP_energy;
        vector<double>  *hltIter0IterL3_bestMatchTP_pt;
        vector<double>  *hltIter0IterL3_bestMatchTP_eta;
        vector<double>  *hltIter0IterL3_bestMatchTP_phi;
        vector<double>  *hltIter0IterL3_bestMatchTP_parentVx;
        vector<double>  *hltIter0IterL3_bestMatchTP_parentVy;
        vector<double>  *hltIter0IterL3_bestMatchTP_parentVz;
        vector<int>     *hltIter0IterL3_bestMatchTP_status;
        vector<int>     *hltIter0IterL3_bestMatchTP_numberOfHits;
        vector<int>     *hltIter0IterL3_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter0IterL3_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter0IterL3_bestMatchTP_sharedFraction;
        vector<int>     *hltIter0IterL3_matchedTPsize;
        vector<float>   *hltIter0IterL3_mva0;
        vector<float>   *hltIter0IterL3_mva1;
        vector<float>   *hltIter0IterL3_mva2;
        vector<float>   *hltIter0IterL3_mva3;
        Int_t           ntpTo_hltIter0IterL3;
        vector<float>   *tpTo_hltIter0IterL3_charge;
        vector<int>     *tpTo_hltIter0IterL3_pdgId;
        vector<double>  *tpTo_hltIter0IterL3_energy;
        vector<double>  *tpTo_hltIter0IterL3_pt;
        vector<double>  *tpTo_hltIter0IterL3_eta;
        vector<double>  *tpTo_hltIter0IterL3_phi;
        vector<double>  *tpTo_hltIter0IterL3_parentVx;
        vector<double>  *tpTo_hltIter0IterL3_parentVy;
        vector<double>  *tpTo_hltIter0IterL3_parentVz;
        vector<int>     *tpTo_hltIter0IterL3_status;
        vector<int>     *tpTo_hltIter0IterL3_numberOfHits;
        vector<int>     *tpTo_hltIter0IterL3_numberOfTrackerHits;
        vector<int>     *tpTo_hltIter0IterL3_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIter0IterL3_gen_charge;
        vector<int>     *tpTo_hltIter0IterL3_gen_pdgId;
        vector<double>  *tpTo_hltIter0IterL3_gen_pt;
        vector<double>  *tpTo_hltIter0IterL3_gen_eta;
        vector<double>  *tpTo_hltIter0IterL3_gen_phi;
        vector<double>  *tpTo_hltIter0IterL3_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIter0IterL3_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIter0IterL3_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIter0IterL3_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIter0IterL3_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIter0IterL3_bestMatchTrk_NValidHits;
        Int_t           nhltIter2IterL3;
        vector<double>  *hltIter2IterL3_pt;
        vector<double>  *hltIter2IterL3_ptError;
        vector<double>  *hltIter2IterL3_eta;
        vector<double>  *hltIter2IterL3_phi;
        vector<int>     *hltIter2IterL3_charge;
        vector<int>     *hltIter2IterL3_matchedL3;
        vector<int>     *hltIter2IterL3_matchedL3NoId;
        vector<float>   *hltIter2IterL3_bestMatchTP_charge;
        vector<int>     *hltIter2IterL3_bestMatchTP_pdgId;
        vector<double>  *hltIter2IterL3_bestMatchTP_energy;
        vector<double>  *hltIter2IterL3_bestMatchTP_pt;
        vector<double>  *hltIter2IterL3_bestMatchTP_eta;
        vector<double>  *hltIter2IterL3_bestMatchTP_phi;
        vector<double>  *hltIter2IterL3_bestMatchTP_parentVx;
        vector<double>  *hltIter2IterL3_bestMatchTP_parentVy;
        vector<double>  *hltIter2IterL3_bestMatchTP_parentVz;
        vector<int>     *hltIter2IterL3_bestMatchTP_status;
        vector<int>     *hltIter2IterL3_bestMatchTP_numberOfHits;
        vector<int>     *hltIter2IterL3_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter2IterL3_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter2IterL3_bestMatchTP_sharedFraction;
        vector<int>     *hltIter2IterL3_matchedTPsize;
        vector<float>   *hltIter2IterL3_mva0;
        vector<float>   *hltIter2IterL3_mva1;
        vector<float>   *hltIter2IterL3_mva2;
        vector<float>   *hltIter2IterL3_mva3;
        Int_t           ntpTo_hltIter2IterL3;
        vector<float>   *tpTo_hltIter2IterL3_charge;
        vector<int>     *tpTo_hltIter2IterL3_pdgId;
        vector<double>  *tpTo_hltIter2IterL3_energy;
        vector<double>  *tpTo_hltIter2IterL3_pt;
        vector<double>  *tpTo_hltIter2IterL3_eta;
        vector<double>  *tpTo_hltIter2IterL3_phi;
        vector<double>  *tpTo_hltIter2IterL3_parentVx;
        vector<double>  *tpTo_hltIter2IterL3_parentVy;
        vector<double>  *tpTo_hltIter2IterL3_parentVz;
        vector<int>     *tpTo_hltIter2IterL3_status;
        vector<int>     *tpTo_hltIter2IterL3_numberOfHits;
        vector<int>     *tpTo_hltIter2IterL3_numberOfTrackerHits;
        vector<int>     *tpTo_hltIter2IterL3_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIter2IterL3_gen_charge;
        vector<int>     *tpTo_hltIter2IterL3_gen_pdgId;
        vector<double>  *tpTo_hltIter2IterL3_gen_pt;
        vector<double>  *tpTo_hltIter2IterL3_gen_eta;
        vector<double>  *tpTo_hltIter2IterL3_gen_phi;
        vector<double>  *tpTo_hltIter2IterL3_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIter2IterL3_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIter2IterL3_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIter2IterL3_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIter2IterL3_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIter2IterL3_bestMatchTrk_NValidHits;
        Int_t           nhltIter0IterL3FromL1Muon;
        vector<double>  *hltIter0IterL3FromL1Muon_pt;
        vector<double>  *hltIter0IterL3FromL1Muon_ptError;
        vector<double>  *hltIter0IterL3FromL1Muon_eta;
        vector<double>  *hltIter0IterL3FromL1Muon_phi;
        vector<int>     *hltIter0IterL3FromL1Muon_charge;
        vector<int>     *hltIter0IterL3FromL1Muon_matchedL3;
        vector<int>     *hltIter0IterL3FromL1Muon_matchedL3NoId;
        vector<float>   *hltIter0IterL3FromL1Muon_bestMatchTP_charge;
        vector<int>     *hltIter0IterL3FromL1Muon_bestMatchTP_pdgId;
        vector<double>  *hltIter0IterL3FromL1Muon_bestMatchTP_energy;
        vector<double>  *hltIter0IterL3FromL1Muon_bestMatchTP_pt;
        vector<double>  *hltIter0IterL3FromL1Muon_bestMatchTP_eta;
        vector<double>  *hltIter0IterL3FromL1Muon_bestMatchTP_phi;
        vector<double>  *hltIter0IterL3FromL1Muon_bestMatchTP_parentVx;
        vector<double>  *hltIter0IterL3FromL1Muon_bestMatchTP_parentVy;
        vector<double>  *hltIter0IterL3FromL1Muon_bestMatchTP_parentVz;
        vector<int>     *hltIter0IterL3FromL1Muon_bestMatchTP_status;
        vector<int>     *hltIter0IterL3FromL1Muon_bestMatchTP_numberOfHits;
        vector<int>     *hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter0IterL3FromL1Muon_bestMatchTP_sharedFraction;
        vector<int>     *hltIter0IterL3FromL1Muon_matchedTPsize;
        vector<float>   *hltIter0IterL3FromL1Muon_mva0;
        vector<float>   *hltIter0IterL3FromL1Muon_mva1;
        vector<float>   *hltIter0IterL3FromL1Muon_mva2;
        vector<float>   *hltIter0IterL3FromL1Muon_mva3;
        Int_t           ntpTo_hltIter0IterL3FromL1Muon;
        vector<float>   *tpTo_hltIter0IterL3FromL1Muon_charge;
        vector<int>     *tpTo_hltIter0IterL3FromL1Muon_pdgId;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_energy;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_pt;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_eta;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_phi;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_parentVx;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_parentVy;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_parentVz;
        vector<int>     *tpTo_hltIter0IterL3FromL1Muon_status;
        vector<int>     *tpTo_hltIter0IterL3FromL1Muon_numberOfHits;
        vector<int>     *tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerHits;
        vector<int>     *tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIter0IterL3FromL1Muon_gen_charge;
        vector<int>     *tpTo_hltIter0IterL3FromL1Muon_gen_pdgId;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_gen_pt;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_gen_eta;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_gen_phi;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_NValidHits;
        Int_t           nhltIter2IterL3FromL1Muon;
        vector<double>  *hltIter2IterL3FromL1Muon_pt;
        vector<double>  *hltIter2IterL3FromL1Muon_ptError;
        vector<double>  *hltIter2IterL3FromL1Muon_eta;
        vector<double>  *hltIter2IterL3FromL1Muon_phi;
        vector<int>     *hltIter2IterL3FromL1Muon_charge;
        vector<int>     *hltIter2IterL3FromL1Muon_matchedL3;
        vector<int>     *hltIter2IterL3FromL1Muon_matchedL3NoId;
        vector<float>   *hltIter2IterL3FromL1Muon_bestMatchTP_charge;
        vector<int>     *hltIter2IterL3FromL1Muon_bestMatchTP_pdgId;
        vector<double>  *hltIter2IterL3FromL1Muon_bestMatchTP_energy;
        vector<double>  *hltIter2IterL3FromL1Muon_bestMatchTP_pt;
        vector<double>  *hltIter2IterL3FromL1Muon_bestMatchTP_eta;
        vector<double>  *hltIter2IterL3FromL1Muon_bestMatchTP_phi;
        vector<double>  *hltIter2IterL3FromL1Muon_bestMatchTP_parentVx;
        vector<double>  *hltIter2IterL3FromL1Muon_bestMatchTP_parentVy;
        vector<double>  *hltIter2IterL3FromL1Muon_bestMatchTP_parentVz;
        vector<int>     *hltIter2IterL3FromL1Muon_bestMatchTP_status;
        vector<int>     *hltIter2IterL3FromL1Muon_bestMatchTP_numberOfHits;
        vector<int>     *hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIter2IterL3FromL1Muon_bestMatchTP_sharedFraction;
        vector<int>     *hltIter2IterL3FromL1Muon_matchedTPsize;
        vector<float>   *hltIter2IterL3FromL1Muon_mva0;
        vector<float>   *hltIter2IterL3FromL1Muon_mva1;
        vector<float>   *hltIter2IterL3FromL1Muon_mva2;
        vector<float>   *hltIter2IterL3FromL1Muon_mva3;
        Int_t           ntpTo_hltIter2IterL3FromL1Muon;
        vector<float>   *tpTo_hltIter2IterL3FromL1Muon_charge;
        vector<int>     *tpTo_hltIter2IterL3FromL1Muon_pdgId;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_energy;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_pt;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_eta;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_phi;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_parentVx;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_parentVy;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_parentVz;
        vector<int>     *tpTo_hltIter2IterL3FromL1Muon_status;
        vector<int>     *tpTo_hltIter2IterL3FromL1Muon_numberOfHits;
        vector<int>     *tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerHits;
        vector<int>     *tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIter2IterL3FromL1Muon_gen_charge;
        vector<int>     *tpTo_hltIter2IterL3FromL1Muon_gen_pdgId;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_gen_pt;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_gen_eta;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_gen_phi;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_NValidHits;
        Int_t           nhltIterL3IOFromL1;
        vector<double>  *hltIterL3IOFromL1_pt;
        vector<double>  *hltIterL3IOFromL1_ptError;
        vector<double>  *hltIterL3IOFromL1_eta;
        vector<double>  *hltIterL3IOFromL1_phi;
        vector<int>     *hltIterL3IOFromL1_charge;
        vector<int>     *hltIterL3IOFromL1_matchedL3;
        vector<int>     *hltIterL3IOFromL1_matchedL3NoId;
        vector<float>   *hltIterL3IOFromL1_bestMatchTP_charge;
        vector<int>     *hltIterL3IOFromL1_bestMatchTP_pdgId;
        vector<double>  *hltIterL3IOFromL1_bestMatchTP_energy;
        vector<double>  *hltIterL3IOFromL1_bestMatchTP_pt;
        vector<double>  *hltIterL3IOFromL1_bestMatchTP_eta;
        vector<double>  *hltIterL3IOFromL1_bestMatchTP_phi;
        vector<double>  *hltIterL3IOFromL1_bestMatchTP_parentVx;
        vector<double>  *hltIterL3IOFromL1_bestMatchTP_parentVy;
        vector<double>  *hltIterL3IOFromL1_bestMatchTP_parentVz;
        vector<int>     *hltIterL3IOFromL1_bestMatchTP_status;
        vector<int>     *hltIterL3IOFromL1_bestMatchTP_numberOfHits;
        vector<int>     *hltIterL3IOFromL1_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIterL3IOFromL1_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIterL3IOFromL1_bestMatchTP_sharedFraction;
        vector<int>     *hltIterL3IOFromL1_matchedTPsize;
        vector<float>   *hltIterL3IOFromL1_mva0;
        vector<float>   *hltIterL3IOFromL1_mva1;
        vector<float>   *hltIterL3IOFromL1_mva2;
        vector<float>   *hltIterL3IOFromL1_mva3;
        Int_t           ntpTo_hltIterL3IOFromL1;
        vector<float>   *tpTo_hltIterL3IOFromL1_charge;
        vector<int>     *tpTo_hltIterL3IOFromL1_pdgId;
        vector<double>  *tpTo_hltIterL3IOFromL1_energy;
        vector<double>  *tpTo_hltIterL3IOFromL1_pt;
        vector<double>  *tpTo_hltIterL3IOFromL1_eta;
        vector<double>  *tpTo_hltIterL3IOFromL1_phi;
        vector<double>  *tpTo_hltIterL3IOFromL1_parentVx;
        vector<double>  *tpTo_hltIterL3IOFromL1_parentVy;
        vector<double>  *tpTo_hltIterL3IOFromL1_parentVz;
        vector<int>     *tpTo_hltIterL3IOFromL1_status;
        vector<int>     *tpTo_hltIterL3IOFromL1_numberOfHits;
        vector<int>     *tpTo_hltIterL3IOFromL1_numberOfTrackerHits;
        vector<int>     *tpTo_hltIterL3IOFromL1_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIterL3IOFromL1_gen_charge;
        vector<int>     *tpTo_hltIterL3IOFromL1_gen_pdgId;
        vector<double>  *tpTo_hltIterL3IOFromL1_gen_pt;
        vector<double>  *tpTo_hltIterL3IOFromL1_gen_eta;
        vector<double>  *tpTo_hltIterL3IOFromL1_gen_phi;
        vector<double>  *tpTo_hltIterL3IOFromL1_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIterL3IOFromL1_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIterL3IOFromL1_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIterL3IOFromL1_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIterL3IOFromL1_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIterL3IOFromL1_bestMatchTrk_NValidHits;
        Int_t           nhltIterL3MuonsNoID;
        vector<double>  *hltIterL3MuonsNoID_pt;
        vector<double>  *hltIterL3MuonsNoID_ptError;
        vector<double>  *hltIterL3MuonsNoID_eta;
        vector<double>  *hltIterL3MuonsNoID_phi;
        vector<int>     *hltIterL3MuonsNoID_charge;
        vector<int>     *hltIterL3MuonsNoID_matchedL3;
        vector<int>     *hltIterL3MuonsNoID_matchedL3NoId;
        vector<float>   *hltIterL3MuonsNoID_bestMatchTP_charge;
        vector<int>     *hltIterL3MuonsNoID_bestMatchTP_pdgId;
        vector<double>  *hltIterL3MuonsNoID_bestMatchTP_energy;
        vector<double>  *hltIterL3MuonsNoID_bestMatchTP_pt;
        vector<double>  *hltIterL3MuonsNoID_bestMatchTP_eta;
        vector<double>  *hltIterL3MuonsNoID_bestMatchTP_phi;
        vector<double>  *hltIterL3MuonsNoID_bestMatchTP_parentVx;
        vector<double>  *hltIterL3MuonsNoID_bestMatchTP_parentVy;
        vector<double>  *hltIterL3MuonsNoID_bestMatchTP_parentVz;
        vector<int>     *hltIterL3MuonsNoID_bestMatchTP_status;
        vector<int>     *hltIterL3MuonsNoID_bestMatchTP_numberOfHits;
        vector<int>     *hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIterL3MuonsNoID_bestMatchTP_sharedFraction;
        vector<int>     *hltIterL3MuonsNoID_matchedTPsize;
        vector<float>   *hltIterL3MuonsNoID_mva0;
        vector<float>   *hltIterL3MuonsNoID_mva1;
        vector<float>   *hltIterL3MuonsNoID_mva2;
        vector<float>   *hltIterL3MuonsNoID_mva3;
        Int_t           ntpTo_hltIterL3MuonsNoID;
        vector<float>   *tpTo_hltIterL3MuonsNoID_charge;
        vector<int>     *tpTo_hltIterL3MuonsNoID_pdgId;
        vector<double>  *tpTo_hltIterL3MuonsNoID_energy;
        vector<double>  *tpTo_hltIterL3MuonsNoID_pt;
        vector<double>  *tpTo_hltIterL3MuonsNoID_eta;
        vector<double>  *tpTo_hltIterL3MuonsNoID_phi;
        vector<double>  *tpTo_hltIterL3MuonsNoID_parentVx;
        vector<double>  *tpTo_hltIterL3MuonsNoID_parentVy;
        vector<double>  *tpTo_hltIterL3MuonsNoID_parentVz;
        vector<int>     *tpTo_hltIterL3MuonsNoID_status;
        vector<int>     *tpTo_hltIterL3MuonsNoID_numberOfHits;
        vector<int>     *tpTo_hltIterL3MuonsNoID_numberOfTrackerHits;
        vector<int>     *tpTo_hltIterL3MuonsNoID_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIterL3MuonsNoID_gen_charge;
        vector<int>     *tpTo_hltIterL3MuonsNoID_gen_pdgId;
        vector<double>  *tpTo_hltIterL3MuonsNoID_gen_pt;
        vector<double>  *tpTo_hltIterL3MuonsNoID_gen_eta;
        vector<double>  *tpTo_hltIterL3MuonsNoID_gen_phi;
        vector<double>  *tpTo_hltIterL3MuonsNoID_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIterL3MuonsNoID_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIterL3MuonsNoID_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIterL3MuonsNoID_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIterL3MuonsNoID_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIterL3MuonsNoID_bestMatchTrk_NValidHits;
        Int_t           nhltIterL3Muons;
        vector<double>  *hltIterL3Muons_pt;
        vector<double>  *hltIterL3Muons_ptError;
        vector<double>  *hltIterL3Muons_eta;
        vector<double>  *hltIterL3Muons_phi;
        vector<int>     *hltIterL3Muons_charge;
        vector<int>     *hltIterL3Muons_matchedL3;
        vector<int>     *hltIterL3Muons_matchedL3NoId;
        vector<float>   *hltIterL3Muons_bestMatchTP_charge;
        vector<int>     *hltIterL3Muons_bestMatchTP_pdgId;
        vector<double>  *hltIterL3Muons_bestMatchTP_energy;
        vector<double>  *hltIterL3Muons_bestMatchTP_pt;
        vector<double>  *hltIterL3Muons_bestMatchTP_eta;
        vector<double>  *hltIterL3Muons_bestMatchTP_phi;
        vector<double>  *hltIterL3Muons_bestMatchTP_parentVx;
        vector<double>  *hltIterL3Muons_bestMatchTP_parentVy;
        vector<double>  *hltIterL3Muons_bestMatchTP_parentVz;
        vector<int>     *hltIterL3Muons_bestMatchTP_status;
        vector<int>     *hltIterL3Muons_bestMatchTP_numberOfHits;
        vector<int>     *hltIterL3Muons_bestMatchTP_numberOfTrackerHits;
        vector<int>     *hltIterL3Muons_bestMatchTP_numberOfTrackerLayers;
        vector<double>  *hltIterL3Muons_bestMatchTP_sharedFraction;
        vector<int>     *hltIterL3Muons_matchedTPsize;
        vector<float>   *hltIterL3Muons_mva0;
        vector<float>   *hltIterL3Muons_mva1;
        vector<float>   *hltIterL3Muons_mva2;
        vector<float>   *hltIterL3Muons_mva3;
        Int_t           ntpTo_hltIterL3Muons;
        vector<float>   *tpTo_hltIterL3Muons_charge;
        vector<int>     *tpTo_hltIterL3Muons_pdgId;
        vector<double>  *tpTo_hltIterL3Muons_energy;
        vector<double>  *tpTo_hltIterL3Muons_pt;
        vector<double>  *tpTo_hltIterL3Muons_eta;
        vector<double>  *tpTo_hltIterL3Muons_phi;
        vector<double>  *tpTo_hltIterL3Muons_parentVx;
        vector<double>  *tpTo_hltIterL3Muons_parentVy;
        vector<double>  *tpTo_hltIterL3Muons_parentVz;
        vector<int>     *tpTo_hltIterL3Muons_status;
        vector<int>     *tpTo_hltIterL3Muons_numberOfHits;
        vector<int>     *tpTo_hltIterL3Muons_numberOfTrackerHits;
        vector<int>     *tpTo_hltIterL3Muons_numberOfTrackerLayers;
        vector<float>   *tpTo_hltIterL3Muons_gen_charge;
        vector<int>     *tpTo_hltIterL3Muons_gen_pdgId;
        vector<double>  *tpTo_hltIterL3Muons_gen_pt;
        vector<double>  *tpTo_hltIterL3Muons_gen_eta;
        vector<double>  *tpTo_hltIterL3Muons_gen_phi;
        vector<double>  *tpTo_hltIterL3Muons_bestMatchTrk_pt;
        vector<double>  *tpTo_hltIterL3Muons_bestMatchTrk_eta;
        vector<double>  *tpTo_hltIterL3Muons_bestMatchTrk_phi;
        vector<int>     *tpTo_hltIterL3Muons_bestMatchTrk_charge;
        vector<double>  *tpTo_hltIterL3Muons_bestMatchTrk_quality;
        vector<int>     *tpTo_hltIterL3Muons_bestMatchTrk_NValidHits;



    // -- Branches -- //
        TBranch        *b_trk_pt;   //!
        TBranch        *b_trk_eta;   //!
        TBranch        *b_trk_phi;   //!
        TBranch        *b_trk_d0;   //!
        TBranch        *b_trk_z0;   //!
        TBranch        *b_trk_rInv;   //!
        TBranch        *b_trk_tanL;   //!
        TBranch        *b_trk_MVA1;   //!
        TBranch        *b_trk_MVA2;   //!
        TBranch        *b_trk_MVA3;   //!
        TBranch        *b_trk_chi2;   //!
        TBranch        *b_trk_bendchi2;   //!
        TBranch        *b_trk_nstub;   //!
        TBranch        *b_trk_lhits;   //!
        TBranch        *b_trk_dhits;   //!
        TBranch        *b_trk_seed;   //!
        TBranch        *b_trk_phiSector;   //!
        TBranch        *b_trk_genuine;   //!
        TBranch        *b_trk_loose;   //!
        TBranch        *b_trk_unknown;   //!
        TBranch        *b_trk_combinatoric;   //!
        TBranch        *b_trk_fake;   //!
        TBranch        *b_trk_matchtp_pdgid;   //!
        TBranch        *b_trk_matchtp_pt;   //!
        TBranch        *b_trk_matchtp_eta;   //!
        TBranch        *b_trk_matchtp_phi;   //!
        TBranch        *b_trk_matchtp_z0;   //!
        TBranch        *b_trk_matchtp_dxy;   //!
        TBranch        *b_stub_x;   //!
        TBranch        *b_stub_y;   //!
        TBranch        *b_stub_z;   //!
        TBranch        *b_stub_isBarrel;   //!
        TBranch        *b_stub_layer;   //!
        TBranch        *b_L1TkMu_pt;   //!
        TBranch        *b_L1TkMu_eta;   //!
        TBranch        *b_L1TkMu_phi;   //!
        TBranch        *b_L1TkMu_trkIsol;   //!
        TBranch        *b_L1TkMu_trkzVtx;   //!
        TBranch        *b_L1TkMu_dR;   //!
        TBranch        *b_L1TkMu_nTracksMatched;   //!
        TBranch        *b_L1TkMu_trackCurvature;   //!
        TBranch        *b_L1TkMu_quality;   //!
        TBranch        *b_L1TkMu_pattern;   //!
        TBranch        *b_L1TkMu_muonDetector;   //!
        TBranch        *b_L1TkMu_TTTpointer;   //!
        TBranch        *b_L1TkMu_muRefHwPt;   //!
        TBranch        *b_L1TkMu_muRefHwDXY;   //!
        TBranch        *b_L1TkMu_muRefHwEta;   //!
        TBranch        *b_L1TkMu_muRefHwPhi;   //!
        TBranch        *b_L1TkMu_muRefHwSign;   //!
        TBranch        *b_L1TkMu_muRefHwSignValid;   //!
        TBranch        *b_L1TkMu_muRefHwQual;   //!
        TBranch        *b_isRealData;   //!
        TBranch        *b_runNum;   //!
        TBranch        *b_lumiBlockNum;   //!
        TBranch        *b_eventNum;   //!
        TBranch        *b_nVertex;   //!
        TBranch        *b_bunchID;   //!
        TBranch        *b_instLumi;   //!
        TBranch        *b_dataPU;   //!
        TBranch        *b_dataPURMS;   //!
        TBranch        *b_bunchLumi;   //!
        TBranch        *b_offlineInstLumi;   //!
        TBranch        *b_offlineDataPU;   //!
        TBranch        *b_offlineDataPURMS;   //!
        TBranch        *b_offlineBunchLumi;   //!
        TBranch        *b_truePU;   //!
        TBranch        *b_qScale;   //!
        TBranch        *b_genEventWeight;   //!
        TBranch        *b_PU_pT_hats;   //!
        TBranch        *b_nGenParticle;   //!
        TBranch        *b_genParticle_ID;   //!
        TBranch        *b_genParticle_status;   //!
        TBranch        *b_genParticle_mother;   //!
        TBranch        *b_genParticle_pt;   //!
        TBranch        *b_genParticle_eta;   //!
        TBranch        *b_genParticle_phi;   //!
        TBranch        *b_genParticle_px;   //!
        TBranch        *b_genParticle_py;   //!
        TBranch        *b_genParticle_pz;   //!
        TBranch        *b_genParticle_energy;   //!
        TBranch        *b_genParticle_charge;   //!
        TBranch        *b_genParticle_isPrompt;   //!
        TBranch        *b_genParticle_isPromptFinalState;   //!
        TBranch        *b_genParticle_isTauDecayProduct;   //!
        TBranch        *b_genParticle_isPromptTauDecayProduct;   //!
        TBranch        *b_genParticle_isDirectPromptTauDecayProductFinalState;   //!
        TBranch        *b_genParticle_isHardProcess;   //!
        TBranch        *b_genParticle_isLastCopy;   //!
        TBranch        *b_genParticle_isLastCopyBeforeFSR;   //!
        TBranch        *b_genParticle_isPromptDecayed;   //!
        TBranch        *b_genParticle_isDecayedLeptonHadron;   //!
        TBranch        *b_genParticle_fromHardProcessBeforeFSR;   //!
        TBranch        *b_genParticle_fromHardProcessDecayed;   //!
        TBranch        *b_genParticle_fromHardProcessFinalState;   //!
        TBranch        *b_genParticle_isMostlyLikePythia6Status3;   //!
        TBranch        *b_vec_firedTrigger;   //!
        TBranch        *b_vec_filterName;   //!
        TBranch        *b_vec_HLTObj_pt;   //!
        TBranch        *b_vec_HLTObj_eta;   //!
        TBranch        *b_vec_HLTObj_phi;   //!
        TBranch        *b_vec_myFiredTrigger;   //!
        TBranch        *b_vec_myFilterName;   //!
        TBranch        *b_vec_myHLTObj_pt;   //!
        TBranch        *b_vec_myHLTObj_eta;   //!
        TBranch        *b_vec_myHLTObj_phi;   //!
        TBranch        *b_nMuon;   //!
        TBranch        *b_muon_pt;   //!
        TBranch        *b_muon_eta;   //!
        TBranch        *b_muon_phi;   //!
        TBranch        *b_muon_px;   //!
        TBranch        *b_muon_py;   //!
        TBranch        *b_muon_pz;   //!
        TBranch        *b_muon_dB;   //!
        TBranch        *b_muon_charge;   //!
        TBranch        *b_muon_isGLB;   //!
        TBranch        *b_muon_isSTA;   //!
        TBranch        *b_muon_isTRK;   //!
        TBranch        *b_muon_isPF;   //!
        TBranch        *b_muon_isTight;   //!
        TBranch        *b_muon_isMedium;   //!
        TBranch        *b_muon_isLoose;   //!
        TBranch        *b_muon_isHighPt;   //!
        TBranch        *b_muon_isHighPtNew;   //!
        TBranch        *b_muon_isSoft;   //!
        TBranch        *b_muon_isLooseTriggerMuon;   //!
        TBranch        *b_muon_isME0Muon;   //!
        TBranch        *b_muon_isGEMMuon;   //!
        TBranch        *b_muon_isRPCMuon;   //!
        TBranch        *b_muon_isGoodMuon_TMOneStationTight;   //!
        TBranch        *b_muon_iso03_sumPt;   //!
        TBranch        *b_muon_iso03_hadEt;   //!
        TBranch        *b_muon_iso03_emEt;   //!
        TBranch        *b_muon_PFIso03_charged;   //!
        TBranch        *b_muon_PFIso03_neutral;   //!
        TBranch        *b_muon_PFIso03_photon;   //!
        TBranch        *b_muon_PFIso03_sumPU;   //!
        TBranch        *b_muon_PFIso04_charged;   //!
        TBranch        *b_muon_PFIso04_neutral;   //!
        TBranch        *b_muon_PFIso04_photon;   //!
        TBranch        *b_muon_PFIso04_sumPU;   //!
        TBranch        *b_muon_PFCluster03_ECAL;   //!
        TBranch        *b_muon_PFCluster03_HCAL;   //!
        TBranch        *b_muon_PFCluster04_ECAL;   //!
        TBranch        *b_muon_PFCluster04_HCAL;   //!
        TBranch        *b_muon_normChi2_global;   //!
        TBranch        *b_muon_nTrackerHit_global;   //!
        TBranch        *b_muon_nTrackerLayer_global;   //!
        TBranch        *b_muon_nPixelHit_global;   //!
        TBranch        *b_muon_nMuonHit_global;   //!
        TBranch        *b_muon_normChi2_inner;   //!
        TBranch        *b_muon_nTrackerHit_inner;   //!
        TBranch        *b_muon_nTrackerLayer_inner;   //!
        TBranch        *b_muon_nPixelHit_inner;   //!
        TBranch        *b_muon_pt_tuneP;   //!
        TBranch        *b_muon_ptError_tuneP;   //!
        TBranch        *b_muon_dxyVTX_best;   //!
        TBranch        *b_muon_dzVTX_best;   //!
        TBranch        *b_muon_nMatchedStation;   //!
        TBranch        *b_muon_nMatchedRPCLayer;   //!
        TBranch        *b_muon_stationMask;   //!
        TBranch        *b_muon_expectedNnumberOfMatchedStations;   //!
        TBranch        *b_nL3Muon;   //!
        TBranch        *b_L3Muon_pt;   //!
        TBranch        *b_L3Muon_eta;   //!
        TBranch        *b_L3Muon_phi;   //!
        TBranch        *b_L3Muon_charge;   //!
        TBranch        *b_L3Muon_trkPt;   //!
        TBranch        *b_nL2Muon;   //!
        TBranch        *b_L2Muon_pt;   //!
        TBranch        *b_L2Muon_eta;   //!
        TBranch        *b_L2Muon_phi;   //!
        TBranch        *b_L2Muon_charge;   //!
        TBranch        *b_L2Muon_trkPt;   //!
        TBranch        *b_nTkMuon;   //!
        TBranch        *b_TkMuon_pt;   //!
        TBranch        *b_TkMuon_eta;   //!
        TBranch        *b_TkMuon_phi;   //!
        TBranch        *b_TkMuon_charge;   //!
        TBranch        *b_TkMuon_trkPt;   //!
        TBranch        *b_nL1Muon;   //!
        TBranch        *b_L1Muon_pt;   //!
        TBranch        *b_L1Muon_eta;   //!
        TBranch        *b_L1Muon_phi;   //!
        TBranch        *b_L1Muon_charge;   //!
        TBranch        *b_L1Muon_quality;   //!
        TBranch        *b_L1Muon_etaAtVtx;   //!
        TBranch        *b_L1Muon_phiAtVtx;   //!
        TBranch        *b_nIterL3OI;   //!
        TBranch        *b_iterL3OI_inner_pt;   //!
        TBranch        *b_iterL3OI_inner_eta;   //!
        TBranch        *b_iterL3OI_inner_phi;   //!
        TBranch        *b_iterL3OI_inner_charge;   //!
        TBranch        *b_iterL3OI_outer_pt;   //!
        TBranch        *b_iterL3OI_outer_eta;   //!
        TBranch        *b_iterL3OI_outer_phi;   //!
        TBranch        *b_iterL3OI_outer_charge;   //!
        TBranch        *b_iterL3OI_global_pt;   //!
        TBranch        *b_iterL3OI_global_eta;   //!
        TBranch        *b_iterL3OI_global_phi;   //!
        TBranch        *b_iterL3OI_global_charge;   //!
        TBranch        *b_nIterL3IOFromL2;   //!
        TBranch        *b_iterL3IOFromL2_inner_pt;   //!
        TBranch        *b_iterL3IOFromL2_inner_eta;   //!
        TBranch        *b_iterL3IOFromL2_inner_phi;   //!
        TBranch        *b_iterL3IOFromL2_inner_charge;   //!
        TBranch        *b_iterL3IOFromL2_outer_pt;   //!
        TBranch        *b_iterL3IOFromL2_outer_eta;   //!
        TBranch        *b_iterL3IOFromL2_outer_phi;   //!
        TBranch        *b_iterL3IOFromL2_outer_charge;   //!
        TBranch        *b_iterL3IOFromL2_global_pt;   //!
        TBranch        *b_iterL3IOFromL2_global_eta;   //!
        TBranch        *b_iterL3IOFromL2_global_phi;   //!
        TBranch        *b_iterL3IOFromL2_global_charge;   //!
        TBranch        *b_nIterL3IOFromL1;   //!
        TBranch        *b_iterL3IOFromL1_pt;   //!
        TBranch        *b_iterL3IOFromL1_eta;   //!
        TBranch        *b_iterL3IOFromL1_phi;   //!
        TBranch        *b_iterL3IOFromL1_charge;   //!
        TBranch        *b_nIterL3FromL2;   //!
        TBranch        *b_iterL3FromL2_inner_pt;   //!
        TBranch        *b_iterL3FromL2_inner_eta;   //!
        TBranch        *b_iterL3FromL2_inner_phi;   //!
        TBranch        *b_iterL3FromL2_inner_charge;   //!
        TBranch        *b_iterL3FromL2_outer_pt;   //!
        TBranch        *b_iterL3FromL2_outer_eta;   //!
        TBranch        *b_iterL3FromL2_outer_phi;   //!
        TBranch        *b_iterL3FromL2_outer_charge;   //!
        TBranch        *b_iterL3FromL2_global_pt;   //!
        TBranch        *b_iterL3FromL2_global_eta;   //!
        TBranch        *b_iterL3FromL2_global_phi;   //!
        TBranch        *b_iterL3FromL2_global_charge;   //!
        TBranch        *b_nIterL3MuonNoID;   //!
        TBranch        *b_iterL3MuonNoID_pt;   //!
        TBranch        *b_iterL3MuonNoID_innerPt;   //!
        TBranch        *b_iterL3MuonNoID_eta;   //!
        TBranch        *b_iterL3MuonNoID_phi;   //!
        TBranch        *b_iterL3MuonNoID_charge;   //!
        TBranch        *b_iterL3MuonNoID_isGLB;   //!
        TBranch        *b_iterL3MuonNoID_isSTA;   //!
        TBranch        *b_iterL3MuonNoID_isTRK;   //!
        TBranch        *b_nIterL3Muon;   //!
        TBranch        *b_iterL3Muon_pt;   //!
        TBranch        *b_iterL3Muon_innerPt;   //!
        TBranch        *b_iterL3Muon_eta;   //!
        TBranch        *b_iterL3Muon_phi;   //!
        TBranch        *b_iterL3Muon_charge;   //!
        TBranch        *b_iterL3Muon_isGLB;   //!
        TBranch        *b_iterL3Muon_isSTA;   //!
        TBranch        *b_iterL3Muon_isTRK;   //!
        TBranch        *b_nhltIterL3OISeedsFromL2Muons;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_dir;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_detId;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_pt;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_pt_val;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_eta;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_phi;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_glob_x;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_glob_y;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_glob_z;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_hasErr;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err0;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err1;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err2;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err3;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err4;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err5;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err6;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err7;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err8;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err9;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err10;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err11;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err12;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err13;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_err14;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_x;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_y;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_dxdz;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_dydz;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_px;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_py;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_pz;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_qbp;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tsos_charge;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_iterL3Matched;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_iterL3Ref;   //!
        TBranch        *b_hltIterL3OISeedsFromL2Muons_tmpL3Ref;   //!
        TBranch        *b_nhltIter0IterL3MuonPixelSeedsFromPixelTracks;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_dir;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_detId;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt_val;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_eta;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_phi;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_x;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_y;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_z;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_hasErr;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err0;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err1;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err2;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err3;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err4;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err5;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err6;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err7;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err8;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err9;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err10;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err11;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err12;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err13;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err14;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_x;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_y;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dxdz;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dydz;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_px;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_py;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pz;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_qbp;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_charge;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Matched;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Ref;   //!
        TBranch        *b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tmpL3Ref;   //!
        TBranch        *b_nhltIter2IterL3MuonPixelSeeds;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_dir;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_detId;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_pt;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_pt_val;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_eta;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_phi;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_glob_x;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_glob_y;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_glob_z;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_hasErr;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err0;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err1;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err2;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err3;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err4;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err5;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err6;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err7;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err8;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err9;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err10;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err11;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err12;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err13;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_err14;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_x;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_y;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_dxdz;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_dydz;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_px;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_py;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_pz;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_qbp;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tsos_charge;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_iterL3Matched;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_iterL3Ref;   //!
        TBranch        *b_hltIter2IterL3MuonPixelSeeds_tmpL3Ref;   //!
        TBranch        *b_nhltIter3IterL3MuonPixelSeeds;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_dir;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_detId;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_pt;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_pt_val;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_eta;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_phi;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_glob_x;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_glob_y;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_glob_z;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_hasErr;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err0;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err1;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err2;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err3;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err4;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err5;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err6;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err7;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err8;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err9;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err10;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err11;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err12;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err13;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_err14;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_x;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_y;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_dxdz;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_dydz;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_px;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_py;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_pz;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_qbp;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tsos_charge;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_iterL3Matched;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_iterL3Ref;   //!
        TBranch        *b_hltIter3IterL3MuonPixelSeeds_tmpL3Ref;   //!
        TBranch        *b_nhltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_dir;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_detId;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt_val;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_eta;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_phi;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_x;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_y;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_z;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_hasErr;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err0;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err1;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err2;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err3;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err4;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err5;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err6;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err7;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err8;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err9;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err10;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err11;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err12;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err13;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err14;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_x;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_y;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dxdz;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dydz;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_px;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_py;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pz;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_qbp;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_charge;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Matched;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Ref;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tmpL3Ref;   //!
        TBranch        *b_nhltIter2IterL3FromL1MuonPixelSeeds;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_dir;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_detId;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt_val;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_eta;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_phi;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_x;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_y;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_z;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_hasErr;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err0;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err1;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err2;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err3;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err4;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err5;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err6;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err7;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err8;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err9;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err10;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err11;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err12;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err13;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err14;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_x;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_y;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_dxdz;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_dydz;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_px;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_py;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_pz;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_qbp;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_charge;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_iterL3Matched;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_iterL3Ref;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonPixelSeeds_tmpL3Ref;   //!
        TBranch        *b_nhltIter3IterL3FromL1MuonPixelSeeds;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_dir;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_detId;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt_val;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_eta;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_phi;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_x;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_y;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_z;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_hasErr;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err0;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err1;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err2;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err3;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err4;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err5;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err6;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err7;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err8;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err9;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err10;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err11;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err12;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err13;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err14;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_x;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_y;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_dxdz;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_dydz;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_px;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_py;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_pz;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_qbp;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_charge;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_iterL3Matched;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_iterL3Ref;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonPixelSeeds_tmpL3Ref;   //!
        TBranch        *b_nhltIterL3OIMuonTrack;   //!
        TBranch        *b_hltIterL3OIMuonTrack_pt;   //!
        TBranch        *b_hltIterL3OIMuonTrack_ptError;   //!
        TBranch        *b_hltIterL3OIMuonTrack_eta;   //!
        TBranch        *b_hltIterL3OIMuonTrack_phi;   //!
        TBranch        *b_hltIterL3OIMuonTrack_charge;   //!
        TBranch        *b_hltIterL3OIMuonTrack_matchedL3;   //!
        TBranch        *b_hltIterL3OIMuonTrack_matchedL3NoId;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_charge;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_energy;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_pt;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_eta;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_phi;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_status;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIterL3OIMuonTrack_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIterL3OIMuonTrack_matchedTPsize;   //!
        TBranch        *b_hltIterL3OIMuonTrack_mva0;   //!
        TBranch        *b_hltIterL3OIMuonTrack_mva1;   //!
        TBranch        *b_hltIterL3OIMuonTrack_mva2;   //!
        TBranch        *b_hltIterL3OIMuonTrack_mva3;   //!
        TBranch        *b_nhltIter0IterL3MuonTrack;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_pt;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_ptError;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_eta;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_phi;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_charge;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_matchedL3;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_matchedL3NoId;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_charge;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_energy;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_pt;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_eta;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_phi;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_status;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_matchedTPsize;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_mva0;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_mva1;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_mva2;   //!
        TBranch        *b_hltIter0IterL3MuonTrack_mva3;   //!
        TBranch        *b_nhltIter2IterL3MuonTrack;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_pt;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_ptError;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_eta;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_phi;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_charge;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_matchedL3;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_matchedL3NoId;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_charge;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_energy;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_pt;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_eta;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_phi;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_status;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_matchedTPsize;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_mva0;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_mva1;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_mva2;   //!
        TBranch        *b_hltIter2IterL3MuonTrack_mva3;   //!
        TBranch        *b_nhltIter3IterL3MuonTrack;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_pt;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_ptError;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_eta;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_phi;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_charge;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_matchedL3;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_matchedL3NoId;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_charge;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_energy;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_pt;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_eta;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_phi;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_status;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_matchedTPsize;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_mva0;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_mva1;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_mva2;   //!
        TBranch        *b_hltIter3IterL3MuonTrack_mva3;   //!
        TBranch        *b_nhltIter0IterL3FromL1MuonTrack;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_pt;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_ptError;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_eta;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_phi;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_charge;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_matchedL3;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_matchedL3NoId;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_charge;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_energy;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_pt;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_eta;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_phi;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_status;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_matchedTPsize;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_mva0;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_mva1;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_mva2;   //!
        TBranch        *b_hltIter0IterL3FromL1MuonTrack_mva3;   //!
        TBranch        *b_nhltIter2IterL3FromL1MuonTrack;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_pt;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_ptError;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_eta;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_phi;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_charge;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_matchedL3;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_matchedL3NoId;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_charge;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_energy;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_pt;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_eta;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_phi;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_status;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_matchedTPsize;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_mva0;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_mva1;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_mva2;   //!
        TBranch        *b_hltIter2IterL3FromL1MuonTrack_mva3;   //!
        TBranch        *b_nhltIter3IterL3FromL1MuonTrack;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_pt;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_ptError;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_eta;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_phi;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_charge;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_matchedL3;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_matchedL3NoId;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_charge;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_energy;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_pt;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_eta;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_phi;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_status;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_matchedTPsize;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_mva0;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_mva1;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_mva2;   //!
        TBranch        *b_hltIter3IterL3FromL1MuonTrack_mva3;   //!

        TBranch        *b_nL3MuonsNoId;   //!
        TBranch        *b_L3MuonsNoId_pt;   //!
        TBranch        *b_L3MuonsNoId_inner_pt;   //!
        TBranch        *b_L3MuonsNoId_inner_ptError;   //!
        TBranch        *b_L3MuonsNoId_eta;   //!
        TBranch        *b_L3MuonsNoId_phi;   //!
        TBranch        *b_L3MuonsNoId_charge;   //!
        TBranch        *b_L3MuonsNoId_isGlobalMuon;   //!
        TBranch        *b_L3MuonsNoId_isStandAloneMuon;   //!
        TBranch        *b_L3MuonsNoId_isTrackerMuon;   //!
        TBranch        *b_L3MuonsNoId_isLooseTriggerMuon;   //!
        TBranch        *b_L3MuonsNoId_isME0Muon;   //!
        TBranch        *b_L3MuonsNoId_isGEMMuon;   //!
        TBranch        *b_L3MuonsNoId_isRPCMuon;   //!
        TBranch        *b_L3MuonsNoId_isGoodMuon_TMOneStationTight;   //!
        TBranch        *b_L3MuonsNoId_numberOfMatchedStations;   //!
        TBranch        *b_L3MuonsNoId_numberOfMatchedRPCLayers;   //!
        TBranch        *b_L3MuonsNoId_expectedNnumberOfMatchedStations;   //!
        TBranch        *b_L3MuonsNoId_inner_normalizedChi2;   //!
        TBranch        *b_L3MuonsNoId_inner_numberOfValidTrackerHits;   //!
        TBranch        *b_L3MuonsNoId_inner_trackerLayersWithMeasurement;   //!
        TBranch        *b_L3MuonsNoId_inner_numberOfValidPixelHits;   //!
        TBranch        *b_L3MuonsNoId_inner_dz_l1vtx;   //!
        TBranch        *b_nL3Muons;   //!
        TBranch        *b_L3Muons_pt;   //!
        TBranch        *b_L3Muons_inner_pt;   //!
        TBranch        *b_L3Muons_inner_ptError;   //!
        TBranch        *b_L3Muons_eta;   //!
        TBranch        *b_L3Muons_phi;   //!
        TBranch        *b_L3Muons_charge;   //!
        TBranch        *b_L3Muons_isGlobalMuon;   //!
        TBranch        *b_L3Muons_isStandAloneMuon;   //!
        TBranch        *b_L3Muons_isTrackerMuon;   //!
        TBranch        *b_L3Muons_isLooseTriggerMuon;   //!
        TBranch        *b_L3Muons_isME0Muon;   //!
        TBranch        *b_L3Muons_isGEMMuon;   //!
        TBranch        *b_L3Muons_isRPCMuon;   //!
        TBranch        *b_L3Muons_isGoodMuon_TMOneStationTight;   //!
        TBranch        *b_L3Muons_numberOfMatchedStations;   //!
        TBranch        *b_L3Muons_numberOfMatchedRPCLayers;   //!
        TBranch        *b_L3Muons_expectedNnumberOfMatchedStations;   //!
        TBranch        *b_L3Muons_inner_normalizedChi2;   //!
        TBranch        *b_L3Muons_inner_numberOfValidTrackerHits;   //!
        TBranch        *b_L3Muons_inner_trackerLayersWithMeasurement;   //!
        TBranch        *b_L3Muons_inner_numberOfValidPixelHits;   //!
        TBranch        *b_L3Muons_inner_dz_l1vtx;   //!

        TBranch        *b_nTP;   //!
        TBranch        *b_TP_charge;   //!
        TBranch        *b_TP_pdgId;   //!
        TBranch        *b_TP_energy;   //!
        TBranch        *b_TP_pt;   //!
        TBranch        *b_TP_eta;   //!
        TBranch        *b_TP_phi;   //!
        TBranch        *b_TP_parentVx;   //!
        TBranch        *b_TP_parentVy;   //!
        TBranch        *b_TP_parentVz;   //!
        TBranch        *b_TP_status;   //!
        TBranch        *b_TP_numberOfHits;   //!
        TBranch        *b_TP_numberOfTrackerHits;   //!
        TBranch        *b_TP_numberOfTrackerLayers;   //!
        TBranch        *b_TP_gen_charge;   //!
        TBranch        *b_TP_gen_pdgId;   //!
        TBranch        *b_TP_gen_pt;   //!
        TBranch        *b_TP_gen_eta;   //!
        TBranch        *b_TP_gen_phi;   //!
        TBranch        *b_TP_bestMatchTrk_pt;   //!
        TBranch        *b_TP_bestMatchTrk_eta;   //!
        TBranch        *b_TP_bestMatchTrk_phi;   //!
        TBranch        *b_TP_bestMatchTrk_charge;   //!
        TBranch        *b_TP_bestMatchTrk_quality;   //!
        TBranch        *b_TP_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltIterL3MuonTrimmedPixelVertices;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_isValid;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_chi2;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_ndof;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_nTracks;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_x;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_xerr;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_y;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_yerr;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_z;   //!
        TBranch        *b_hltIterL3MuonTrimmedPixelVertices_zerr;   //!
        TBranch        *b_nhltIterL3FromL1MuonTrimmedPixelVertices;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_isValid;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_chi2;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_ndof;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_nTracks;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_x;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_xerr;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_y;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_yerr;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_z;   //!
        TBranch        *b_hltIterL3FromL1MuonTrimmedPixelVertices_zerr;   //!

        TBranch        *b_nhltPhase2L3OI;   //!
        TBranch        *b_hltPhase2L3OI_pt;   //!
        TBranch        *b_hltPhase2L3OI_ptError;   //!
        TBranch        *b_hltPhase2L3OI_eta;   //!
        TBranch        *b_hltPhase2L3OI_phi;   //!
        TBranch        *b_hltPhase2L3OI_charge;   //!
        TBranch        *b_hltPhase2L3OI_matchedL3;   //!
        TBranch        *b_hltPhase2L3OI_matchedL3NoId;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_charge;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_pdgId;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_energy;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_pt;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_eta;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_phi;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_parentVx;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_parentVy;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_parentVz;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_status;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltPhase2L3OI_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltPhase2L3OI_matchedTPsize;   //!
        TBranch        *b_hltPhase2L3OI_mva0;   //!
        TBranch        *b_hltPhase2L3OI_mva1;   //!
        TBranch        *b_hltPhase2L3OI_mva2;   //!
        TBranch        *b_hltPhase2L3OI_mva3;   //!
        TBranch        *b_ntpTo_hltPhase2L3OI;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_pdgId;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_energy;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_parentVx;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_parentVy;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_parentVz;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_status;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_numberOfHits;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_gen_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_gen_pdgId;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_gen_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_gen_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_gen_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltPhase2L3OI_bestMatchTrk_NValidHits;   //!

        TBranch        *b_nhltIter0Phase2L3FromL1TkMuon;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_pt;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_ptError;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_eta;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_phi;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_charge;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_matchedL3;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_matchedL3NoId;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_charge;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_energy;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pt;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_eta;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_phi;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_status;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_matchedTPsize;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_mva0;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_mva1;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_mva2;   //!
        TBranch        *b_hltIter0Phase2L3FromL1TkMuon_mva3;   //!
        TBranch        *b_ntpTo_hltIter0Phase2L3FromL1TkMuon;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_charge;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_pdgId;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_energy;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_pt;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_eta;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_phi;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_parentVx;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_parentVy;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_parentVz;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_status;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfHits;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_charge;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pt;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_eta;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_phi;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltIter2Phase2L3FromL1TkMuon;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_pt;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_ptError;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_eta;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_phi;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_charge;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_matchedL3;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_matchedL3NoId;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_charge;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_energy;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pt;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_eta;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_phi;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_status;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_matchedTPsize;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_mva0;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_mva1;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_mva2;   //!
        TBranch        *b_hltIter2Phase2L3FromL1TkMuon_mva3;   //!
        TBranch        *b_ntpTo_hltIter2Phase2L3FromL1TkMuon;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_charge;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_pdgId;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_energy;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_pt;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_eta;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_phi;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_parentVx;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_parentVy;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_parentVz;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_status;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfHits;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_charge;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pt;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_eta;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_phi;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltPhase2L3IOFromL1;   //!
        TBranch        *b_hltPhase2L3IOFromL1_pt;   //!
        TBranch        *b_hltPhase2L3IOFromL1_ptError;   //!
        TBranch        *b_hltPhase2L3IOFromL1_eta;   //!
        TBranch        *b_hltPhase2L3IOFromL1_phi;   //!
        TBranch        *b_hltPhase2L3IOFromL1_charge;   //!
        TBranch        *b_hltPhase2L3IOFromL1_matchedL3;   //!
        TBranch        *b_hltPhase2L3IOFromL1_matchedL3NoId;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_charge;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_pdgId;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_energy;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_pt;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_eta;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_phi;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_parentVx;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_parentVy;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_parentVz;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_status;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltPhase2L3IOFromL1_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltPhase2L3IOFromL1_matchedTPsize;   //!
        TBranch        *b_hltPhase2L3IOFromL1_mva0;   //!
        TBranch        *b_hltPhase2L3IOFromL1_mva1;   //!
        TBranch        *b_hltPhase2L3IOFromL1_mva2;   //!
        TBranch        *b_hltPhase2L3IOFromL1_mva3;   //!
        TBranch        *b_ntpTo_hltPhase2L3IOFromL1;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_pdgId;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_energy;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_parentVx;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_parentVy;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_parentVz;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_status;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_numberOfHits;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_gen_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_gen_pdgId;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_gen_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_gen_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_gen_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltPhase2L3MuonsNoID;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_pt;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_ptError;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_eta;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_phi;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_charge;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_matchedL3;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_matchedL3NoId;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_charge;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_pdgId;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_energy;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_pt;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_eta;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_phi;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_parentVx;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_parentVy;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_parentVz;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_status;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_matchedTPsize;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_mva0;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_mva1;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_mva2;   //!
        TBranch        *b_hltPhase2L3MuonsNoID_mva3;   //!
        TBranch        *b_ntpTo_hltPhase2L3MuonsNoID;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_pdgId;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_energy;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_parentVx;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_parentVy;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_parentVz;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_status;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_numberOfHits;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_gen_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_gen_pdgId;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_gen_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_gen_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_gen_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltPhase2L3Muons;   //!
        TBranch        *b_hltPhase2L3Muons_pt;   //!
        TBranch        *b_hltPhase2L3Muons_ptError;   //!
        TBranch        *b_hltPhase2L3Muons_eta;   //!
        TBranch        *b_hltPhase2L3Muons_phi;   //!
        TBranch        *b_hltPhase2L3Muons_charge;   //!
        TBranch        *b_hltPhase2L3Muons_matchedL3;   //!
        TBranch        *b_hltPhase2L3Muons_matchedL3NoId;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_charge;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_pdgId;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_energy;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_pt;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_eta;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_phi;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_parentVx;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_parentVy;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_parentVz;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_status;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltPhase2L3Muons_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltPhase2L3Muons_matchedTPsize;   //!
        TBranch        *b_hltPhase2L3Muons_mva0;   //!
        TBranch        *b_hltPhase2L3Muons_mva1;   //!
        TBranch        *b_hltPhase2L3Muons_mva2;   //!
        TBranch        *b_hltPhase2L3Muons_mva3;   //!

        // -- Isolation
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0;   //!
        TBranch        *b_hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p000;   //!
        TBranch        *b_hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p000;   //!
        TBranch        *b_hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p030;   //!
        TBranch        *b_hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p030;   //!
        TBranch        *b_hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p050;   //!
        TBranch        *b_hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p050;   //!
        TBranch        *b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p00minEEM0p00minEHad0p00;   //!
        TBranch        *b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00;   //!
        TBranch        *b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p04minEEM0p00minEHad0p00;   //!
        TBranch        *b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p00minEEM0p00minEHad0p00;   //!
        TBranch        *b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p02minEEM0p00minEHad0p00;   //!
        TBranch        *b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p04minEEM0p00minEHad0p00;   //!
        TBranch        *b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p00minEEM0p00minEHad0p00;   //!
        TBranch        *b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p02minEEM0p00minEHad0p00;   //!
        TBranch        *b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p04minEEM0p00minEHad0p00;   //!

        TBranch        *b_ntpTo_hltPhase2L3Muons;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_pdgId;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_energy;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_parentVx;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_parentVy;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_parentVz;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_status;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_numberOfHits;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_gen_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_gen_pdgId;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_gen_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_gen_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_gen_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltPhase2L3Muons_bestMatchTrk_NValidHits;   //!

        TBranch        *b_nhltIterL3OI;   //!
        TBranch        *b_hltIterL3OI_pt;   //!
        TBranch        *b_hltIterL3OI_ptError;   //!
        TBranch        *b_hltIterL3OI_eta;   //!
        TBranch        *b_hltIterL3OI_phi;   //!
        TBranch        *b_hltIterL3OI_charge;   //!
        TBranch        *b_hltIterL3OI_matchedL3;   //!
        TBranch        *b_hltIterL3OI_matchedL3NoId;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_charge;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_energy;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_pt;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_eta;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_phi;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_status;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIterL3OI_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIterL3OI_matchedTPsize;   //!
        TBranch        *b_hltIterL3OI_mva0;   //!
        TBranch        *b_hltIterL3OI_mva1;   //!
        TBranch        *b_hltIterL3OI_mva2;   //!
        TBranch        *b_hltIterL3OI_mva3;   //!
        TBranch        *b_ntpTo_hltIterL3OI;   //!
        TBranch        *b_tpTo_hltIterL3OI_charge;   //!
        TBranch        *b_tpTo_hltIterL3OI_pdgId;   //!
        TBranch        *b_tpTo_hltIterL3OI_energy;   //!
        TBranch        *b_tpTo_hltIterL3OI_pt;   //!
        TBranch        *b_tpTo_hltIterL3OI_eta;   //!
        TBranch        *b_tpTo_hltIterL3OI_phi;   //!
        TBranch        *b_tpTo_hltIterL3OI_parentVx;   //!
        TBranch        *b_tpTo_hltIterL3OI_parentVy;   //!
        TBranch        *b_tpTo_hltIterL3OI_parentVz;   //!
        TBranch        *b_tpTo_hltIterL3OI_status;   //!
        TBranch        *b_tpTo_hltIterL3OI_numberOfHits;   //!
        TBranch        *b_tpTo_hltIterL3OI_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIterL3OI_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIterL3OI_gen_charge;   //!
        TBranch        *b_tpTo_hltIterL3OI_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIterL3OI_gen_pt;   //!
        TBranch        *b_tpTo_hltIterL3OI_gen_eta;   //!
        TBranch        *b_tpTo_hltIterL3OI_gen_phi;   //!
        TBranch        *b_tpTo_hltIterL3OI_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIterL3OI_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIterL3OI_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIterL3OI_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIterL3OI_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIterL3OI_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltIter0IterL3;   //!
        TBranch        *b_hltIter0IterL3_pt;   //!
        TBranch        *b_hltIter0IterL3_ptError;   //!
        TBranch        *b_hltIter0IterL3_eta;   //!
        TBranch        *b_hltIter0IterL3_phi;   //!
        TBranch        *b_hltIter0IterL3_charge;   //!
        TBranch        *b_hltIter0IterL3_matchedL3;   //!
        TBranch        *b_hltIter0IterL3_matchedL3NoId;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_charge;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_energy;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_pt;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_eta;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_phi;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_status;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter0IterL3_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter0IterL3_matchedTPsize;   //!
        TBranch        *b_hltIter0IterL3_mva0;   //!
        TBranch        *b_hltIter0IterL3_mva1;   //!
        TBranch        *b_hltIter0IterL3_mva2;   //!
        TBranch        *b_hltIter0IterL3_mva3;   //!
        TBranch        *b_ntpTo_hltIter0IterL3;   //!
        TBranch        *b_tpTo_hltIter0IterL3_charge;   //!
        TBranch        *b_tpTo_hltIter0IterL3_pdgId;   //!
        TBranch        *b_tpTo_hltIter0IterL3_energy;   //!
        TBranch        *b_tpTo_hltIter0IterL3_pt;   //!
        TBranch        *b_tpTo_hltIter0IterL3_eta;   //!
        TBranch        *b_tpTo_hltIter0IterL3_phi;   //!
        TBranch        *b_tpTo_hltIter0IterL3_parentVx;   //!
        TBranch        *b_tpTo_hltIter0IterL3_parentVy;   //!
        TBranch        *b_tpTo_hltIter0IterL3_parentVz;   //!
        TBranch        *b_tpTo_hltIter0IterL3_status;   //!
        TBranch        *b_tpTo_hltIter0IterL3_numberOfHits;   //!
        TBranch        *b_tpTo_hltIter0IterL3_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIter0IterL3_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIter0IterL3_gen_charge;   //!
        TBranch        *b_tpTo_hltIter0IterL3_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIter0IterL3_gen_pt;   //!
        TBranch        *b_tpTo_hltIter0IterL3_gen_eta;   //!
        TBranch        *b_tpTo_hltIter0IterL3_gen_phi;   //!
        TBranch        *b_tpTo_hltIter0IterL3_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIter0IterL3_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIter0IterL3_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIter0IterL3_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIter0IterL3_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIter0IterL3_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltIter2IterL3;   //!
        TBranch        *b_hltIter2IterL3_pt;   //!
        TBranch        *b_hltIter2IterL3_ptError;   //!
        TBranch        *b_hltIter2IterL3_eta;   //!
        TBranch        *b_hltIter2IterL3_phi;   //!
        TBranch        *b_hltIter2IterL3_charge;   //!
        TBranch        *b_hltIter2IterL3_matchedL3;   //!
        TBranch        *b_hltIter2IterL3_matchedL3NoId;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_charge;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_energy;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_pt;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_eta;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_phi;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_status;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter2IterL3_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter2IterL3_matchedTPsize;   //!
        TBranch        *b_hltIter2IterL3_mva0;   //!
        TBranch        *b_hltIter2IterL3_mva1;   //!
        TBranch        *b_hltIter2IterL3_mva2;   //!
        TBranch        *b_hltIter2IterL3_mva3;   //!
        TBranch        *b_ntpTo_hltIter2IterL3;   //!
        TBranch        *b_tpTo_hltIter2IterL3_charge;   //!
        TBranch        *b_tpTo_hltIter2IterL3_pdgId;   //!
        TBranch        *b_tpTo_hltIter2IterL3_energy;   //!
        TBranch        *b_tpTo_hltIter2IterL3_pt;   //!
        TBranch        *b_tpTo_hltIter2IterL3_eta;   //!
        TBranch        *b_tpTo_hltIter2IterL3_phi;   //!
        TBranch        *b_tpTo_hltIter2IterL3_parentVx;   //!
        TBranch        *b_tpTo_hltIter2IterL3_parentVy;   //!
        TBranch        *b_tpTo_hltIter2IterL3_parentVz;   //!
        TBranch        *b_tpTo_hltIter2IterL3_status;   //!
        TBranch        *b_tpTo_hltIter2IterL3_numberOfHits;   //!
        TBranch        *b_tpTo_hltIter2IterL3_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIter2IterL3_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIter2IterL3_gen_charge;   //!
        TBranch        *b_tpTo_hltIter2IterL3_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIter2IterL3_gen_pt;   //!
        TBranch        *b_tpTo_hltIter2IterL3_gen_eta;   //!
        TBranch        *b_tpTo_hltIter2IterL3_gen_phi;   //!
        TBranch        *b_tpTo_hltIter2IterL3_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIter2IterL3_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIter2IterL3_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIter2IterL3_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIter2IterL3_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIter2IterL3_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltIter0IterL3FromL1Muon;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_pt;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_ptError;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_eta;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_phi;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_charge;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_matchedL3;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_matchedL3NoId;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_charge;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_energy;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_pt;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_eta;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_phi;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_status;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_matchedTPsize;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_mva0;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_mva1;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_mva2;   //!
        TBranch        *b_hltIter0IterL3FromL1Muon_mva3;   //!
        TBranch        *b_ntpTo_hltIter0IterL3FromL1Muon;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_charge;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_pdgId;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_energy;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_pt;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_eta;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_phi;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_parentVx;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_parentVy;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_parentVz;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_status;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_numberOfHits;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_gen_charge;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_gen_pt;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_gen_eta;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_gen_phi;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltIter2IterL3FromL1Muon;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_pt;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_ptError;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_eta;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_phi;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_charge;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_matchedL3;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_matchedL3NoId;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_charge;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_energy;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_pt;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_eta;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_phi;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_status;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_matchedTPsize;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_mva0;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_mva1;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_mva2;   //!
        TBranch        *b_hltIter2IterL3FromL1Muon_mva3;   //!
        TBranch        *b_ntpTo_hltIter2IterL3FromL1Muon;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_charge;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_pdgId;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_energy;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_pt;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_eta;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_phi;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_parentVx;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_parentVy;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_parentVz;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_status;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_numberOfHits;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_gen_charge;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_gen_pt;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_gen_eta;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_gen_phi;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltIterL3IOFromL1;   //!
        TBranch        *b_hltIterL3IOFromL1_pt;   //!
        TBranch        *b_hltIterL3IOFromL1_ptError;   //!
        TBranch        *b_hltIterL3IOFromL1_eta;   //!
        TBranch        *b_hltIterL3IOFromL1_phi;   //!
        TBranch        *b_hltIterL3IOFromL1_charge;   //!
        TBranch        *b_hltIterL3IOFromL1_matchedL3;   //!
        TBranch        *b_hltIterL3IOFromL1_matchedL3NoId;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_charge;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_energy;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_pt;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_eta;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_phi;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_status;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIterL3IOFromL1_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIterL3IOFromL1_matchedTPsize;   //!
        TBranch        *b_hltIterL3IOFromL1_mva0;   //!
        TBranch        *b_hltIterL3IOFromL1_mva1;   //!
        TBranch        *b_hltIterL3IOFromL1_mva2;   //!
        TBranch        *b_hltIterL3IOFromL1_mva3;   //!
        TBranch        *b_ntpTo_hltIterL3IOFromL1;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_charge;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_pdgId;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_energy;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_pt;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_eta;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_phi;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_parentVx;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_parentVy;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_parentVz;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_status;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_numberOfHits;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_gen_charge;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_gen_pt;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_gen_eta;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_gen_phi;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIterL3IOFromL1_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltIterL3MuonsNoID;   //!
        TBranch        *b_hltIterL3MuonsNoID_pt;   //!
        TBranch        *b_hltIterL3MuonsNoID_ptError;   //!
        TBranch        *b_hltIterL3MuonsNoID_eta;   //!
        TBranch        *b_hltIterL3MuonsNoID_phi;   //!
        TBranch        *b_hltIterL3MuonsNoID_charge;   //!
        TBranch        *b_hltIterL3MuonsNoID_matchedL3;   //!
        TBranch        *b_hltIterL3MuonsNoID_matchedL3NoId;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_charge;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_energy;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_pt;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_eta;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_phi;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_status;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIterL3MuonsNoID_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIterL3MuonsNoID_matchedTPsize;   //!
        TBranch        *b_hltIterL3MuonsNoID_mva0;   //!
        TBranch        *b_hltIterL3MuonsNoID_mva1;   //!
        TBranch        *b_hltIterL3MuonsNoID_mva2;   //!
        TBranch        *b_hltIterL3MuonsNoID_mva3;   //!
        TBranch        *b_ntpTo_hltIterL3MuonsNoID;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_charge;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_pdgId;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_energy;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_pt;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_eta;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_phi;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_parentVx;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_parentVy;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_parentVz;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_status;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_numberOfHits;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_gen_charge;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_gen_pt;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_gen_eta;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_gen_phi;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_NValidHits;   //!
        TBranch        *b_nhltIterL3Muons;   //!
        TBranch        *b_hltIterL3Muons_pt;   //!
        TBranch        *b_hltIterL3Muons_ptError;   //!
        TBranch        *b_hltIterL3Muons_eta;   //!
        TBranch        *b_hltIterL3Muons_phi;   //!
        TBranch        *b_hltIterL3Muons_charge;   //!
        TBranch        *b_hltIterL3Muons_matchedL3;   //!
        TBranch        *b_hltIterL3Muons_matchedL3NoId;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_charge;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_pdgId;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_energy;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_pt;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_eta;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_phi;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_parentVx;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_parentVy;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_parentVz;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_status;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_numberOfHits;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_numberOfTrackerHits;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_numberOfTrackerLayers;   //!
        TBranch        *b_hltIterL3Muons_bestMatchTP_sharedFraction;   //!
        TBranch        *b_hltIterL3Muons_matchedTPsize;   //!
        TBranch        *b_hltIterL3Muons_mva0;   //!
        TBranch        *b_hltIterL3Muons_mva1;   //!
        TBranch        *b_hltIterL3Muons_mva2;   //!
        TBranch        *b_hltIterL3Muons_mva3;   //!
        TBranch        *b_ntpTo_hltIterL3Muons;   //!
        TBranch        *b_tpTo_hltIterL3Muons_charge;   //!
        TBranch        *b_tpTo_hltIterL3Muons_pdgId;   //!
        TBranch        *b_tpTo_hltIterL3Muons_energy;   //!
        TBranch        *b_tpTo_hltIterL3Muons_pt;   //!
        TBranch        *b_tpTo_hltIterL3Muons_eta;   //!
        TBranch        *b_tpTo_hltIterL3Muons_phi;   //!
        TBranch        *b_tpTo_hltIterL3Muons_parentVx;   //!
        TBranch        *b_tpTo_hltIterL3Muons_parentVy;   //!
        TBranch        *b_tpTo_hltIterL3Muons_parentVz;   //!
        TBranch        *b_tpTo_hltIterL3Muons_status;   //!
        TBranch        *b_tpTo_hltIterL3Muons_numberOfHits;   //!
        TBranch        *b_tpTo_hltIterL3Muons_numberOfTrackerHits;   //!
        TBranch        *b_tpTo_hltIterL3Muons_numberOfTrackerLayers;   //!
        TBranch        *b_tpTo_hltIterL3Muons_gen_charge;   //!
        TBranch        *b_tpTo_hltIterL3Muons_gen_pdgId;   //!
        TBranch        *b_tpTo_hltIterL3Muons_gen_pt;   //!
        TBranch        *b_tpTo_hltIterL3Muons_gen_eta;   //!
        TBranch        *b_tpTo_hltIterL3Muons_gen_phi;   //!
        TBranch        *b_tpTo_hltIterL3Muons_bestMatchTrk_pt;   //!
        TBranch        *b_tpTo_hltIterL3Muons_bestMatchTrk_eta;   //!
        TBranch        *b_tpTo_hltIterL3Muons_bestMatchTrk_phi;   //!
        TBranch        *b_tpTo_hltIterL3Muons_bestMatchTrk_charge;   //!
        TBranch        *b_tpTo_hltIterL3Muons_bestMatchTrk_quality;   //!
        TBranch        *b_tpTo_hltIterL3Muons_bestMatchTrk_NValidHits;   //!
};

MuonHLTNtuple::MuonHLTNtuple(TChain *tree, vector<TString> branch_tags ) : fChain(0) 
{
    Init(tree);

    fChain->SetBranchStatus( "*", 0 );
    for( auto& tag: branch_tags) {
        fChain->SetBranchStatus( tag+"*", 1 );
    }

    fChain->SetBranchStatus( "n*", 1 );
    fChain->SetBranchStatus( "isRealData", 1 );
    fChain->SetBranchStatus( "runNum", 1 );
    fChain->SetBranchStatus( "lumiBlockNum", 1 );
    fChain->SetBranchStatus( "eventNum", 1 );
    fChain->SetBranchStatus( "nVertex", 1 );
    fChain->SetBranchStatus( "truePU", 1 );
    fChain->SetBranchStatus( "qScale", 1 );
    fChain->SetBranchStatus( "genEventWeight", 1 );
    fChain->SetBranchStatus( "PU_pT_hats", 1 );

    branch_names = {};
    TObjArray* listofb = fChain->GetListOfBranches();
    for(const auto b: *listofb) {
        TString name = b->GetName();
        if(fChain->GetBranchStatus(name.Data())) {
            branch_names.push_back(name);
        }
    }
}

MuonHLTNtuple::~MuonHLTNtuple()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

vector<TString> MuonHLTNtuple::getIsoTags()
{

    vector<TString> tags = {};

    vector<double> trkIsodRs     = { 0.3 };  // , 0.4 };
    vector<double> trkIsodRVetos = { 0.005 };  // , 0.01 };  // , 0.02 };
    vector<double> trkIsodzs     = { 0.1, 0.2, 0.25 };  // { 0.1, 0.15, 0.2, 0.25 };
    vector<double> trkIsodrs     = { 0.1, 0.2 };
    vector<double> trkIsoChi2s   = { 1.0E64 };
    vector<double> trkIsoPtMins  = { -1.0 };
    vector<TString> trkIsoTypes   = { "Regional", "Full", "RegionalNew", "Offline" };

    for(auto& trkIsodR: trkIsodRs) {
        for(auto& trkIsodRVeto: trkIsodRVetos) {
            for(auto& trkIsodz: trkIsodzs) {
                for(auto& trkIsodr: trkIsodrs) {
                    for(auto& trkIsoChi2: trkIsoChi2s) {
                        for(auto& trkIsoPtMin: trkIsoPtMins) {
                            for(auto& trkIsoType: trkIsoTypes) {
                                TString trkIsoTag = TString::Format("trkIso%sdR%.1fdRVeto%.3fdz%.2fdr%.2fChisq%.0fPtMin%.1f", trkIsoType.Data(), trkIsodR, trkIsodRVeto, trkIsodz, trkIsodr, trkIsoChi2, trkIsoPtMin);
                                trkIsoTag = trkIsoTag.ReplaceAll(".", "p");
                                trkIsoTag = trkIsoTag.ReplaceAll(TString::Format("Chisq%.0f", 1.0E64), "ChisqInf");
                                trkIsoTag = trkIsoTag.ReplaceAll("PtMin-1p0", "PtMin0p0");

                                tags.push_back( trkIsoTag );
                            }
                        }
                    }
                }
            }
        }
    }


    vector<double> pfIsodRs     = { 0.3 };  // , 0.4 };  // , 0.5 };
    vector<double> pfIsodRVetos = { 0.00, 0.03, 0.05};  // , 0.07, 0.1 };

    for(auto& pfIsodR: pfIsodRs) {
        for(auto& pfIsodRVeto: pfIsodRVetos) {
            TString pfIsoTag = TString::Format("IsodR%.1fdRVeto%.3f", pfIsodR, pfIsodRVeto);
            pfIsoTag = pfIsoTag.ReplaceAll(".", "p");

            tags.push_back( "pfEcal"+pfIsoTag );
            tags.push_back( "pfHcal"+pfIsoTag );
            // tags.push_back( "pfHgcal"+pfIsoTag );
        }
    }

    vector<double> lcIsodRs        = { 0.2 };  // , 0.3, 0.4 };
    vector<double> lcIsodRVetoEMs  = { 0.00, 0.02, 0.04 };
    vector<double> lcIsodRVetoHads = { 0.00, 0.02, 0.04 };
    vector<double> lcIsoMinEEMs    = { 0.00 };  // , 0.02, 0.05, 0.07 };
    vector<double> lcIsoMinEHads   = { 0.00 };  // , 0.02, 0.05, 0.07 };

    for(auto& lcIsodR: lcIsodRs) {
        for(auto& lcIsodRVetoEM: lcIsodRVetoEMs) {
            for(auto& lcIsodRVetoHad: lcIsodRVetoHads) {
                for(auto& lcIsoMinEEM: lcIsoMinEEMs) {
                    for(auto& lcIsoMinEHad: lcIsoMinEHads) {
                        TString HgcalLCIsoTag = TString::Format("HgcalLCIsodR%.1fdRVetoEM%.2fdRVetoHad%.2fminEEM%.2fminEHad%.2f", lcIsodR, lcIsodRVetoEM, lcIsodRVetoHad, lcIsoMinEEM, lcIsoMinEHad);
                        HgcalLCIsoTag = HgcalLCIsoTag.ReplaceAll(".", "p");

                        tags.push_back( "pf"+HgcalLCIsoTag );
                    }
                }
            }
        }
    }

    return tags;
}

void MuonHLTNtuple::PrintTemplate( TString head )
{
    cout << endl;
    cout << TString::Format("vector<Object> MuonHLTNtuple::get_%s()", head.Data()) << endl;
    cout << TString::Format("{") << endl;
    cout << TString::Format("    vector<Object> out = {};" ) << endl;
    cout << TString::Format("    if(%s_pt == 0 || %s_pt == nullptr)", head.Data(), head.Data()) << endl;
    cout << TString::Format("        return out;" ) << endl;
    cout << endl;
    cout << TString::Format("    for(unsigned i=0; i<%s_pt->size(); ++i) {", head.Data()) << endl;
    cout << TString::Format("        Object obj = Object( %s_pt->at(i), %s_eta->at(i), %s_phi->at(i) );", head.Data(), head.Data(), head.Data()) << endl;
    cout << endl;

    TObjArray* listofb = fChain->GetListOfBranches();
    for(const auto b: *listofb) {
        TString brname = b->GetName();

        if(brname.BeginsWith(head+"_")) {
            TString varname = brname;
            varname = varname.ReplaceAll(head+"_", "");
            // if( varname == "pt" || varname == "eta" || varname == "phi" )
            //     continue;
            cout << TString::Format("        obj.addVar( \"%s\", %s->at(i) );", varname.Data(), brname.Data() ) << endl;
        }
    }

    cout << endl;
    cout << "        out.push_back(obj);" << endl;
    cout << "    }" << endl;
    cout << endl;
    cout << "    return out;" << endl;
    cout << "}" << endl;
}

vector<Object> MuonHLTNtuple::get_GenParticles()
{
    vector<Object> out = {};

    for(int i=0; i<nGenParticle; ++i) {

        if( genParticle_status[i] != 1 && genParticle_isHardProcess[i] != 1 )
            continue;

        Object obj = Object( genParticle_pt[i], genParticle_eta[i], genParticle_phi[i] );

        obj.addVar( "ID", genParticle_ID[i] );
        obj.addVar( "status", genParticle_status[i] );
        obj.addVar( "mother", genParticle_mother[i] );
        obj.addVar( "pt", genParticle_pt[i] );
        obj.addVar( "eta", genParticle_eta[i] );
        obj.addVar( "phi", genParticle_phi[i] );
        obj.addVar( "px", genParticle_px[i] );
        obj.addVar( "py", genParticle_py[i] );
        obj.addVar( "pz", genParticle_pz[i] );
        obj.addVar( "energy", genParticle_energy[i] );
        obj.addVar( "charge", genParticle_charge[i] );
        obj.addVar( "isPrompt", genParticle_isPrompt[i] );
        obj.addVar( "isPromptFinalState", genParticle_isPromptFinalState[i] );
        obj.addVar( "isTauDecayProduct", genParticle_isTauDecayProduct[i] );
        obj.addVar( "isPromptTauDecayProduct", genParticle_isPromptTauDecayProduct[i] );
        obj.addVar( "isDirectPromptTauDecayProductFinalState", genParticle_isDirectPromptTauDecayProductFinalState[i] );
        obj.addVar( "isHardProcess", genParticle_isHardProcess[i] );
        obj.addVar( "isLastCopy", genParticle_isLastCopy[i] );
        obj.addVar( "isLastCopyBeforeFSR", genParticle_isLastCopyBeforeFSR[i] );
        obj.addVar( "isPromptDecayed", genParticle_isPromptDecayed[i] );
        obj.addVar( "isDecayedLeptonHadron", genParticle_isDecayedLeptonHadron[i] );
        obj.addVar( "fromHardProcessBeforeFSR", genParticle_fromHardProcessBeforeFSR[i] );
        obj.addVar( "fromHardProcessDecayed", genParticle_fromHardProcessDecayed[i] );
        obj.addVar( "fromHardProcessFinalState", genParticle_fromHardProcessFinalState[i] );
        obj.addVar( "isMostlyLikePythia6Status3", genParticle_isMostlyLikePythia6Status3[i] );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_L1TkMuons()
{
    vector<Object> out = {};
    if(L1TkMu_pt == 0 || L1TkMu_pt == nullptr)
        return out;

    for(unsigned i=0; i<L1TkMu_pt->size(); ++i) {
        Object obj = Object( L1TkMu_pt->at(i), L1TkMu_eta->at(i), L1TkMu_phi->at(i) );

        obj.addVar( "pt", L1TkMu_pt->at(i) );
        obj.addVar( "eta", L1TkMu_eta->at(i) );
        obj.addVar( "phi", L1TkMu_phi->at(i) );
        obj.addVar( "trkIsol", L1TkMu_trkIsol->at(i) );
        obj.addVar( "trkzVtx", L1TkMu_trkzVtx->at(i) );
        obj.addVar( "dR", L1TkMu_dR->at(i) );
        obj.addVar( "nTracksMatched", L1TkMu_nTracksMatched->at(i) );
        obj.addVar( "trackCurvature", L1TkMu_trackCurvature->at(i) );
        obj.addVar( "quality", L1TkMu_quality->at(i) );
        obj.addVar( "pattern", L1TkMu_pattern->at(i) );
        obj.addVar( "muonDetector", L1TkMu_muonDetector->at(i) );
        obj.addVar( "TTTpointer", L1TkMu_TTTpointer->at(i) );
        obj.addVar( "muRefHwPt", L1TkMu_muRefHwPt->at(i) );
        obj.addVar( "muRefHwDXY", L1TkMu_muRefHwDXY->at(i) );
        obj.addVar( "muRefHwEta", L1TkMu_muRefHwEta->at(i) );
        obj.addVar( "muRefHwPhi", L1TkMu_muRefHwPhi->at(i) );
        obj.addVar( "muRefHwSign", L1TkMu_muRefHwSign->at(i) );
        obj.addVar( "muRefHwSignValid", L1TkMu_muRefHwSignValid->at(i) );
        obj.addVar( "muRefHwQual", L1TkMu_muRefHwQual->at(i) );
        obj.addVar( "offlinePt", TkMuonOfflineEt(L1TkMu_pt->at(i), L1TkMu_eta->at(i)) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_L1TTs()
{
    vector<Object> out = {};
    if(trk_pt == 0 || trk_pt == nullptr)
        return out;

    for(unsigned i=0; i<trk_pt->size(); ++i) {

        if( trk_pt->at(i) < 5.0 )
            continue;

        Object obj = Object( trk_pt->at(i), trk_eta->at(i), trk_phi->at(i) );

        obj.addVar( "pt", trk_pt->at(i) );
        obj.addVar( "eta", trk_eta->at(i) );
        obj.addVar( "phi", trk_phi->at(i) );
        obj.addVar( "d0", trk_d0->at(i) );
        obj.addVar( "z0", trk_z0->at(i) );
        obj.addVar( "rInv", trk_rInv->at(i) );
        obj.addVar( "tanL", trk_tanL->at(i) );
        // obj.addVar( "MVA1", trk_MVA1->at(i) );
        // obj.addVar( "MVA2", trk_MVA2->at(i) );
        // obj.addVar( "MVA3", trk_MVA3->at(i) );
        // obj.addVar( "chi2", trk_chi2->at(i) );
        // obj.addVar( "bendchi2", trk_bendchi2->at(i) );
        // obj.addVar( "nstub", trk_nstub->at(i) );
        // obj.addVar( "lhits", trk_lhits->at(i) );
        // obj.addVar( "dhits", trk_dhits->at(i) );
        // obj.addVar( "seed", trk_seed->at(i) );
        // obj.addVar( "phiSector", trk_phiSector->at(i) );
        // obj.addVar( "genuine", trk_genuine->at(i) );
        // obj.addVar( "loose", trk_loose->at(i) );
        // obj.addVar( "unknown", trk_unknown->at(i) );
        // obj.addVar( "combinatoric", trk_combinatoric->at(i) );
        // obj.addVar( "fake", trk_fake->at(i) );
        // obj.addVar( "matchtp_pdgid", trk_matchtp_pdgid->at(i) );
        // obj.addVar( "matchtp_pt", trk_matchtp_pt->at(i) );
        // obj.addVar( "matchtp_eta", trk_matchtp_eta->at(i) );
        // obj.addVar( "matchtp_phi", trk_matchtp_phi->at(i) );
        // obj.addVar( "matchtp_z0", trk_matchtp_z0->at(i) );
        // obj.addVar( "matchtp_dxy", trk_matchtp_dxy->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_L1Muons()
{
    vector<Object> out = {};
    for(int i=0; i<nL1Muon; ++i) {

        if(L1Muon_quality[i] < 12)
            continue;

        Object obj = Object( L1Muon_pt[i], L1Muon_eta[i], L1Muon_phi[i] );

        obj.addVar( "pt", L1Muon_pt[i] );
        obj.addVar( "eta", L1Muon_eta[i] );
        obj.addVar( "phi", L1Muon_phi[i] );
        obj.addVar( "charge", L1Muon_charge[i] );
        obj.addVar( "quality", L1Muon_quality[i] );
        obj.addVar( "etaAtVtx", L1Muon_etaAtVtx[i] );
        obj.addVar( "phiAtVtx", L1Muon_phiAtVtx[i] );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_L2Muons()
{
    vector<Object> out = {};

    for(int i=0; i<nL2Muon; ++i) {

        Object obj = Object( L2Muon_pt[i], L2Muon_eta[i], L2Muon_phi[i] );

        obj.addVar( "pt", L2Muon_pt[i] );
        obj.addVar( "eta", L2Muon_eta[i] );
        obj.addVar( "phi", L2Muon_phi[i] );
        obj.addVar( "charge", L2Muon_charge[i] );
        obj.addVar( "trkPt", L2Muon_trkPt[i] );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_L3MuonsNoId()
{
    vector<Object> out = {};
    if(L3MuonsNoId_pt == 0 || L3MuonsNoId_pt == nullptr)
        return out;

    for(unsigned i=0; i<L3MuonsNoId_pt->size(); ++i) {
        Object obj = Object( L3MuonsNoId_pt->at(i), L3MuonsNoId_eta->at(i), L3MuonsNoId_phi->at(i) );

        obj.addVar( "pt", L3MuonsNoId_pt->at(i) );
        obj.addVar( "inner_pt", L3MuonsNoId_inner_pt->at(i) );
        obj.addVar( "inner_ptError", L3MuonsNoId_inner_ptError->at(i) );
        obj.addVar( "eta", L3MuonsNoId_eta->at(i) );
        obj.addVar( "phi", L3MuonsNoId_phi->at(i) );
        obj.addVar( "charge", L3MuonsNoId_charge->at(i) );
        obj.addVar( "isGlobalMuon", L3MuonsNoId_isGlobalMuon->at(i) );
        obj.addVar( "isStandAloneMuon", L3MuonsNoId_isStandAloneMuon->at(i) );
        obj.addVar( "isTrackerMuon", L3MuonsNoId_isTrackerMuon->at(i) );
        obj.addVar( "isLooseTriggerMuon", L3MuonsNoId_isLooseTriggerMuon->at(i) );
        obj.addVar( "isME0Muon", L3MuonsNoId_isME0Muon->at(i) );
        obj.addVar( "isGEMMuon", L3MuonsNoId_isGEMMuon->at(i) );
        obj.addVar( "isRPCMuon", L3MuonsNoId_isRPCMuon->at(i) );
        obj.addVar( "isGoodMuon_TMOneStationTight", L3MuonsNoId_isGoodMuon_TMOneStationTight->at(i) );
        obj.addVar( "numberOfMatchedStations", L3MuonsNoId_numberOfMatchedStations->at(i) );
        obj.addVar( "numberOfMatchedRPCLayers", L3MuonsNoId_numberOfMatchedRPCLayers->at(i) );
        obj.addVar( "expectedNnumberOfMatchedStations", L3MuonsNoId_expectedNnumberOfMatchedStations->at(i) );
        obj.addVar( "inner_normalizedChi2", L3MuonsNoId_inner_normalizedChi2->at(i) );
        obj.addVar( "inner_numberOfValidTrackerHits", L3MuonsNoId_inner_numberOfValidTrackerHits->at(i) );
        obj.addVar( "inner_trackerLayersWithMeasurement", L3MuonsNoId_inner_trackerLayersWithMeasurement->at(i) );
        obj.addVar( "inner_numberOfValidPixelHits", L3MuonsNoId_inner_numberOfValidPixelHits->at(i) );
        obj.addVar( "inner_dz_l1vtx", L3MuonsNoId_inner_dz_l1vtx->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_L3Muons()
{
    vector<Object> out = {};
    if(L3Muons_pt == 0 || L3Muons_pt == nullptr)
        return out;

    for(unsigned i=0; i<L3Muons_pt->size(); ++i) {
        Object obj = Object( L3Muons_pt->at(i), L3Muons_eta->at(i), L3Muons_phi->at(i) );

        obj.addVar( "pt", L3Muons_pt->at(i) );
        obj.addVar( "inner_pt", L3Muons_inner_pt->at(i) );
        obj.addVar( "inner_ptError", L3Muons_inner_ptError->at(i) );
        obj.addVar( "eta", L3Muons_eta->at(i) );
        obj.addVar( "phi", L3Muons_phi->at(i) );
        obj.addVar( "charge", L3Muons_charge->at(i) );
        obj.addVar( "isGlobalMuon", L3Muons_isGlobalMuon->at(i) );
        obj.addVar( "isStandAloneMuon", L3Muons_isStandAloneMuon->at(i) );
        obj.addVar( "isTrackerMuon", L3Muons_isTrackerMuon->at(i) );
        obj.addVar( "isLooseTriggerMuon", L3Muons_isLooseTriggerMuon->at(i) );
        obj.addVar( "isME0Muon", L3Muons_isME0Muon->at(i) );
        obj.addVar( "isGEMMuon", L3Muons_isGEMMuon->at(i) );
        obj.addVar( "isRPCMuon", L3Muons_isRPCMuon->at(i) );
        obj.addVar( "isGoodMuon_TMOneStationTight", L3Muons_isGoodMuon_TMOneStationTight->at(i) );
        obj.addVar( "numberOfMatchedStations", L3Muons_numberOfMatchedStations->at(i) );
        obj.addVar( "numberOfMatchedRPCLayers", L3Muons_numberOfMatchedRPCLayers->at(i) );
        obj.addVar( "expectedNnumberOfMatchedStations", L3Muons_expectedNnumberOfMatchedStations->at(i) );
        obj.addVar( "inner_normalizedChi2", L3Muons_inner_normalizedChi2->at(i) );
        obj.addVar( "inner_numberOfValidTrackerHits", L3Muons_inner_numberOfValidTrackerHits->at(i) );
        obj.addVar( "inner_trackerLayersWithMeasurement", L3Muons_inner_trackerLayersWithMeasurement->at(i) );
        obj.addVar( "inner_numberOfValidPixelHits", L3Muons_inner_numberOfValidPixelHits->at(i) );
        obj.addVar( "inner_dz_l1vtx", L3Muons_inner_dz_l1vtx->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_iterL3OI()
{
    vector<Object> out = {};
    for(int i=0; i<nIterL3OI; ++i) {
        Object obj = Object( iterL3OI_inner_pt[i], iterL3OI_inner_eta[i], iterL3OI_inner_phi[i] );

        obj.addVar( "inner_pt", iterL3OI_inner_pt[i] );
        obj.addVar( "inner_eta", iterL3OI_inner_eta[i] );
        obj.addVar( "inner_phi", iterL3OI_inner_phi[i] );
        obj.addVar( "inner_charge", iterL3OI_inner_charge[i] );
        obj.addVar( "outer_pt", iterL3OI_outer_pt[i] );
        obj.addVar( "outer_eta", iterL3OI_outer_eta[i] );
        obj.addVar( "outer_phi", iterL3OI_outer_phi[i] );
        obj.addVar( "outer_charge", iterL3OI_outer_charge[i] );
        obj.addVar( "global_pt", iterL3OI_global_pt[i] );
        obj.addVar( "global_eta", iterL3OI_global_eta[i] );
        obj.addVar( "global_phi", iterL3OI_global_phi[i] );
        obj.addVar( "global_charge", iterL3OI_global_charge[i] );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_iterL3IOFromL2()
{
    vector<Object> out = {};
    for(int i=0; i<nIterL3IOFromL2; ++i) {
        Object obj = Object( iterL3IOFromL2_inner_pt[i], iterL3IOFromL2_inner_eta[i], iterL3IOFromL2_inner_phi[i] );

        obj.addVar( "inner_pt", iterL3IOFromL2_inner_pt[i] );
        obj.addVar( "inner_eta", iterL3IOFromL2_inner_eta[i] );
        obj.addVar( "inner_phi", iterL3IOFromL2_inner_phi[i] );
        obj.addVar( "inner_charge", iterL3IOFromL2_inner_charge[i] );
        obj.addVar( "outer_pt", iterL3IOFromL2_outer_pt[i] );
        obj.addVar( "outer_eta", iterL3IOFromL2_outer_eta[i] );
        obj.addVar( "outer_phi", iterL3IOFromL2_outer_phi[i] );
        obj.addVar( "outer_charge", iterL3IOFromL2_outer_charge[i] );
        obj.addVar( "global_pt", iterL3IOFromL2_global_pt[i] );
        obj.addVar( "global_eta", iterL3IOFromL2_global_eta[i] );
        obj.addVar( "global_phi", iterL3IOFromL2_global_phi[i] );
        obj.addVar( "global_charge", iterL3IOFromL2_global_charge[i] );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_iterL3FromL2()
{
    vector<Object> out = {};
    for(int i=0; i<nIterL3FromL2; ++i) {
        Object obj = Object( iterL3FromL2_inner_pt[i], iterL3FromL2_inner_eta[i], iterL3FromL2_inner_phi[i] );

        obj.addVar( "inner_pt", iterL3FromL2_inner_pt[i] );
        obj.addVar( "inner_eta", iterL3FromL2_inner_eta[i] );
        obj.addVar( "inner_phi", iterL3FromL2_inner_phi[i] );
        obj.addVar( "inner_charge", iterL3FromL2_inner_charge[i] );
        obj.addVar( "outer_pt", iterL3FromL2_outer_pt[i] );
        obj.addVar( "outer_eta", iterL3FromL2_outer_eta[i] );
        obj.addVar( "outer_phi", iterL3FromL2_outer_phi[i] );
        obj.addVar( "outer_charge", iterL3FromL2_outer_charge[i] );
        obj.addVar( "global_pt", iterL3FromL2_global_pt[i] );
        obj.addVar( "global_eta", iterL3FromL2_global_eta[i] );
        obj.addVar( "global_phi", iterL3FromL2_global_phi[i] );
        obj.addVar( "global_charge", iterL3FromL2_global_charge[i] );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_iterL3IOFromL1()
{
    vector<Object> out = {};
    for(int i=0; i<nIterL3IOFromL1; ++i) {
        Object obj = Object( iterL3IOFromL1_pt[i], iterL3IOFromL1_eta[i], iterL3IOFromL1_phi[i] );

        obj.addVar( "pt", iterL3IOFromL1_pt[i] );
        obj.addVar( "eta", iterL3IOFromL1_eta[i] );
        obj.addVar( "phi", iterL3IOFromL1_phi[i] );
        obj.addVar( "charge", iterL3IOFromL1_charge[i] );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_iterL3MuonNoID()
{
    vector<Object> out = {};
    for(int i=0; i<nIterL3MuonNoID; ++i) {
        Object obj = Object( iterL3MuonNoID_pt[i], iterL3MuonNoID_eta[i], iterL3MuonNoID_phi[i] );

        obj.addVar( "pt", iterL3MuonNoID_pt[i] );
        obj.addVar( "innerPt", iterL3MuonNoID_innerPt[i] );
        obj.addVar( "eta", iterL3MuonNoID_eta[i] );
        obj.addVar( "phi", iterL3MuonNoID_phi[i] );
        obj.addVar( "charge", iterL3MuonNoID_charge[i] );
        obj.addVar( "isGLB", iterL3MuonNoID_isGLB[i] );
        obj.addVar( "isSTA", iterL3MuonNoID_isSTA[i] );
        obj.addVar( "isTRK", iterL3MuonNoID_isTRK[i] );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_iterL3Muon()
{
    vector<Object> out = {};
    for(int i=0; i<nIterL3Muon; ++i) {
        Object obj = Object( iterL3Muon_pt[i], iterL3Muon_eta[i], iterL3Muon_phi[i] );

        obj.addVar( "pt", iterL3Muon_pt[i] );
        obj.addVar( "innerPt", iterL3Muon_innerPt[i] );
        obj.addVar( "eta", iterL3Muon_eta[i] );
        obj.addVar( "phi", iterL3Muon_phi[i] );
        obj.addVar( "charge", iterL3Muon_charge[i] );
        obj.addVar( "isGLB", iterL3Muon_isGLB[i] );
        obj.addVar( "isSTA", iterL3Muon_isSTA[i] );
        obj.addVar( "isTRK", iterL3Muon_isTRK[i] );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_HLTObjects( TString filter )
{
    vector<Object> out = {};
    if(vec_myFilterName == 0 || vec_myFilterName == nullptr)
        return out;

    for(unsigned i=0; i<vec_myFilterName->size(); ++i) {

        TString ifilter = TString(vec_myFilterName->at(i));
        if( !ifilter.Contains(filter) )
            continue;

        Object obj = Object( vec_myHLTObj_pt->at(i), vec_myHLTObj_eta->at(i), vec_myHLTObj_phi->at(i) );

        obj.addStrVar( "filter", ifilter );
        obj.addVar( "pt", vec_myHLTObj_pt->at(i) );
        obj.addVar( "eta", vec_myHLTObj_eta->at(i) );
        obj.addVar( "phi", vec_myHLTObj_phi->at(i) );

        if(ifilter.Contains("hltL1TkSingleMuFiltered22") ||
           ifilter.Contains("hltL1TkDoubleMuFiltered7") ||
           ifilter.Contains("hltL1TkSingleMuFiltered15") ||
           ifilter.Contains("hltDoubleMuon7DZ1p0") ||
           ifilter.Contains("hltL1TripleMuFiltered3") ||
           ifilter.Contains("hltTripleMuon3DZ1p0") ||
           ifilter.Contains("hltTripleMuon3DR0")) {
            obj.addVar( "offlinePt", TkMuonOfflineEt(vec_myHLTObj_pt->at(i), vec_myHLTObj_eta->at(i)) );
        }

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltPhase2L3OI()
{
    vector<Object> out = {};
    if(hltPhase2L3OI_pt == 0 || hltPhase2L3OI_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltPhase2L3OI_pt->size(); ++i) {
        Object obj = Object( hltPhase2L3OI_pt->at(i), hltPhase2L3OI_eta->at(i), hltPhase2L3OI_phi->at(i) );

        obj.addVar( "pt", hltPhase2L3OI_pt->at(i) );
        obj.addVar( "ptError", hltPhase2L3OI_ptError->at(i) );
        obj.addVar( "eta", hltPhase2L3OI_eta->at(i) );
        obj.addVar( "phi", hltPhase2L3OI_phi->at(i) );
        obj.addVar( "charge", hltPhase2L3OI_charge->at(i) );
        obj.addVar( "matchedL3", hltPhase2L3OI_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltPhase2L3OI_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltPhase2L3OI_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltPhase2L3OI_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltPhase2L3OI_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltPhase2L3OI_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltPhase2L3OI_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltPhase2L3OI_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltPhase2L3OI_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltPhase2L3OI_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltPhase2L3OI_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltPhase2L3OI_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltPhase2L3OI_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltPhase2L3OI_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltPhase2L3OI_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltPhase2L3OI_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltPhase2L3OI_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltPhase2L3OI_mva0->at(i) );
        obj.addVar( "mva1", hltPhase2L3OI_mva1->at(i) );
        obj.addVar( "mva2", hltPhase2L3OI_mva2->at(i) );
        obj.addVar( "mva3", hltPhase2L3OI_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIter0Phase2L3FromL1TkMuon()
{
    vector<Object> out = {};
    if(hltIter0Phase2L3FromL1TkMuon_pt == 0 || hltIter0Phase2L3FromL1TkMuon_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIter0Phase2L3FromL1TkMuon_pt->size(); ++i) {
        Object obj = Object( hltIter0Phase2L3FromL1TkMuon_pt->at(i), hltIter0Phase2L3FromL1TkMuon_eta->at(i), hltIter0Phase2L3FromL1TkMuon_phi->at(i) );

        obj.addVar( "pt", hltIter0Phase2L3FromL1TkMuon_pt->at(i) );
        obj.addVar( "ptError", hltIter0Phase2L3FromL1TkMuon_ptError->at(i) );
        obj.addVar( "eta", hltIter0Phase2L3FromL1TkMuon_eta->at(i) );
        obj.addVar( "phi", hltIter0Phase2L3FromL1TkMuon_phi->at(i) );
        obj.addVar( "charge", hltIter0Phase2L3FromL1TkMuon_charge->at(i) );
        obj.addVar( "matchedL3", hltIter0Phase2L3FromL1TkMuon_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIter0Phase2L3FromL1TkMuon_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIter0Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIter0Phase2L3FromL1TkMuon_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIter0Phase2L3FromL1TkMuon_mva0->at(i) );
        obj.addVar( "mva1", hltIter0Phase2L3FromL1TkMuon_mva1->at(i) );
        obj.addVar( "mva2", hltIter0Phase2L3FromL1TkMuon_mva2->at(i) );
        obj.addVar( "mva3", hltIter0Phase2L3FromL1TkMuon_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIter2Phase2L3FromL1TkMuon()
{
    vector<Object> out = {};
    if(hltIter2Phase2L3FromL1TkMuon_pt == 0 || hltIter2Phase2L3FromL1TkMuon_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIter2Phase2L3FromL1TkMuon_pt->size(); ++i) {
        Object obj = Object( hltIter2Phase2L3FromL1TkMuon_pt->at(i), hltIter2Phase2L3FromL1TkMuon_eta->at(i), hltIter2Phase2L3FromL1TkMuon_phi->at(i) );

        obj.addVar( "pt", hltIter2Phase2L3FromL1TkMuon_pt->at(i) );
        obj.addVar( "ptError", hltIter2Phase2L3FromL1TkMuon_ptError->at(i) );
        obj.addVar( "eta", hltIter2Phase2L3FromL1TkMuon_eta->at(i) );
        obj.addVar( "phi", hltIter2Phase2L3FromL1TkMuon_phi->at(i) );
        obj.addVar( "charge", hltIter2Phase2L3FromL1TkMuon_charge->at(i) );
        obj.addVar( "matchedL3", hltIter2Phase2L3FromL1TkMuon_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIter2Phase2L3FromL1TkMuon_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIter2Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIter2Phase2L3FromL1TkMuon_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIter2Phase2L3FromL1TkMuon_mva0->at(i) );
        obj.addVar( "mva1", hltIter2Phase2L3FromL1TkMuon_mva1->at(i) );
        obj.addVar( "mva2", hltIter2Phase2L3FromL1TkMuon_mva2->at(i) );
        obj.addVar( "mva3", hltIter2Phase2L3FromL1TkMuon_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltPhase2L3IOFromL1()
{
    vector<Object> out = {};
    if(hltPhase2L3IOFromL1_pt == 0 || hltPhase2L3IOFromL1_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltPhase2L3IOFromL1_pt->size(); ++i) {
        Object obj = Object( hltPhase2L3IOFromL1_pt->at(i), hltPhase2L3IOFromL1_eta->at(i), hltPhase2L3IOFromL1_phi->at(i) );

        obj.addVar( "pt", hltPhase2L3IOFromL1_pt->at(i) );
        obj.addVar( "ptError", hltPhase2L3IOFromL1_ptError->at(i) );
        obj.addVar( "eta", hltPhase2L3IOFromL1_eta->at(i) );
        obj.addVar( "phi", hltPhase2L3IOFromL1_phi->at(i) );
        obj.addVar( "charge", hltPhase2L3IOFromL1_charge->at(i) );
        obj.addVar( "matchedL3", hltPhase2L3IOFromL1_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltPhase2L3IOFromL1_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltPhase2L3IOFromL1_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltPhase2L3IOFromL1_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltPhase2L3IOFromL1_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltPhase2L3IOFromL1_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltPhase2L3IOFromL1_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltPhase2L3IOFromL1_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltPhase2L3IOFromL1_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltPhase2L3IOFromL1_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltPhase2L3IOFromL1_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltPhase2L3IOFromL1_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltPhase2L3IOFromL1_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltPhase2L3IOFromL1_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltPhase2L3IOFromL1_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltPhase2L3IOFromL1_mva0->at(i) );
        obj.addVar( "mva1", hltPhase2L3IOFromL1_mva1->at(i) );
        obj.addVar( "mva2", hltPhase2L3IOFromL1_mva2->at(i) );
        obj.addVar( "mva3", hltPhase2L3IOFromL1_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltPhase2L3MuonsNoID()
{
    vector<Object> out = {};
    if(hltPhase2L3MuonsNoID_pt == 0 || hltPhase2L3MuonsNoID_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltPhase2L3MuonsNoID_pt->size(); ++i) {
        Object obj = Object( hltPhase2L3MuonsNoID_pt->at(i), hltPhase2L3MuonsNoID_eta->at(i), hltPhase2L3MuonsNoID_phi->at(i) );

        obj.addVar( "pt", hltPhase2L3MuonsNoID_pt->at(i) );
        obj.addVar( "ptError", hltPhase2L3MuonsNoID_ptError->at(i) );
        obj.addVar( "eta", hltPhase2L3MuonsNoID_eta->at(i) );
        obj.addVar( "phi", hltPhase2L3MuonsNoID_phi->at(i) );
        obj.addVar( "charge", hltPhase2L3MuonsNoID_charge->at(i) );
        obj.addVar( "matchedL3", hltPhase2L3MuonsNoID_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltPhase2L3MuonsNoID_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltPhase2L3MuonsNoID_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltPhase2L3MuonsNoID_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltPhase2L3MuonsNoID_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltPhase2L3MuonsNoID_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltPhase2L3MuonsNoID_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltPhase2L3MuonsNoID_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltPhase2L3MuonsNoID_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltPhase2L3MuonsNoID_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltPhase2L3MuonsNoID_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltPhase2L3MuonsNoID_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltPhase2L3MuonsNoID_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltPhase2L3MuonsNoID_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltPhase2L3MuonsNoID_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltPhase2L3MuonsNoID_mva0->at(i) );
        obj.addVar( "mva1", hltPhase2L3MuonsNoID_mva1->at(i) );
        obj.addVar( "mva2", hltPhase2L3MuonsNoID_mva2->at(i) );
        obj.addVar( "mva3", hltPhase2L3MuonsNoID_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltPhase2L3Muons()
{
    vector<Object> out = {};
    if(hltPhase2L3Muons_pt == 0 || hltPhase2L3Muons_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltPhase2L3Muons_pt->size(); ++i) {
        Object obj = Object( hltPhase2L3Muons_pt->at(i), hltPhase2L3Muons_eta->at(i), hltPhase2L3Muons_phi->at(i) );

        obj.addVar( "pt", hltPhase2L3Muons_pt->at(i) );
        obj.addVar( "ptError", hltPhase2L3Muons_ptError->at(i) );
        obj.addVar( "eta", hltPhase2L3Muons_eta->at(i) );
        obj.addVar( "phi", hltPhase2L3Muons_phi->at(i) );
        obj.addVar( "charge", hltPhase2L3Muons_charge->at(i) );
        obj.addVar( "matchedL3", hltPhase2L3Muons_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltPhase2L3Muons_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltPhase2L3Muons_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltPhase2L3Muons_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltPhase2L3Muons_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltPhase2L3Muons_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltPhase2L3Muons_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltPhase2L3Muons_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltPhase2L3Muons_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltPhase2L3Muons_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltPhase2L3Muons_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltPhase2L3Muons_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltPhase2L3Muons_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltPhase2L3Muons_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltPhase2L3Muons_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltPhase2L3Muons_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltPhase2L3Muons_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltPhase2L3Muons_mva0->at(i) );
        obj.addVar( "mva1", hltPhase2L3Muons_mva1->at(i) );
        obj.addVar( "mva2", hltPhase2L3Muons_mva2->at(i) );
        obj.addVar( "mva3", hltPhase2L3Muons_mva3->at(i) );

        // -- Isolation
        //obj.addVar( "trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoFulldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoFulldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoFulldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoFulldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoFulldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoFulldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0", hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0->at(i) );
        //obj.addVar( "pfEcalIsodR0p3dRVeto0p000", hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p000->at(i) );
        //obj.addVar( "pfHcalIsodR0p3dRVeto0p000", hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p000->at(i) );
        //obj.addVar( "pfEcalIsodR0p3dRVeto0p030", hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p030->at(i) );
        //obj.addVar( "pfHcalIsodR0p3dRVeto0p030", hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p030->at(i) );
        //obj.addVar( "pfEcalIsodR0p3dRVeto0p050", hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p050->at(i) );
        //obj.addVar( "pfHcalIsodR0p3dRVeto0p050", hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p050->at(i) );
        //obj.addVar( "pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p00minEEM0p00minEHad0p00", hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p00minEEM0p00minEHad0p00->at(i) );
        //obj.addVar( "pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00", hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00->at(i) );
        //obj.addVar( "pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p04minEEM0p00minEHad0p00", hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p04minEEM0p00minEHad0p00->at(i) );
        //obj.addVar( "pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p00minEEM0p00minEHad0p00", hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p00minEEM0p00minEHad0p00->at(i) );
        //obj.addVar( "pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p02minEEM0p00minEHad0p00", hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p02minEEM0p00minEHad0p00->at(i) );
        //obj.addVar( "pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p04minEEM0p00minEHad0p00", hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p04minEEM0p00minEHad0p00->at(i) );
        //obj.addVar( "pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p00minEEM0p00minEHad0p00", hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p00minEEM0p00minEHad0p00->at(i) );
        //obj.addVar( "pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p02minEEM0p00minEHad0p00", hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p02minEEM0p00minEHad0p00->at(i) );
        //obj.addVar( "pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p04minEEM0p00minEHad0p00", hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p04minEEM0p00minEHad0p00->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltPhase2L3OI()
{
    vector<Object> out = {};
    if(tpTo_hltPhase2L3OI_pt == 0 || tpTo_hltPhase2L3OI_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltPhase2L3OI_pt->size(); ++i) {
        Object obj = Object( tpTo_hltPhase2L3OI_pt->at(i), tpTo_hltPhase2L3OI_eta->at(i), tpTo_hltPhase2L3OI_phi->at(i) );

        obj.addVar( "charge", tpTo_hltPhase2L3OI_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltPhase2L3OI_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltPhase2L3OI_energy->at(i) );
        obj.addVar( "pt", tpTo_hltPhase2L3OI_pt->at(i) );
        obj.addVar( "eta", tpTo_hltPhase2L3OI_eta->at(i) );
        obj.addVar( "phi", tpTo_hltPhase2L3OI_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltPhase2L3OI_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltPhase2L3OI_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltPhase2L3OI_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltPhase2L3OI_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltPhase2L3OI_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltPhase2L3OI_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltPhase2L3OI_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltPhase2L3OI_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltPhase2L3OI_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltPhase2L3OI_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltPhase2L3OI_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltPhase2L3OI_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltPhase2L3OI_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltPhase2L3OI_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltPhase2L3OI_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltPhase2L3OI_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltPhase2L3OI_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltPhase2L3OI_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIter0Phase2L3FromL1TkMuon()
{
    vector<Object> out = {};
    if(tpTo_hltIter0Phase2L3FromL1TkMuon_pt == 0 || tpTo_hltIter0Phase2L3FromL1TkMuon_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIter0Phase2L3FromL1TkMuon_pt->size(); ++i) {

        if( tpTo_hltIter0Phase2L3FromL1TkMuon_pt->at(i) < 5.0 || tpTo_hltIter0Phase2L3FromL1TkMuon_status->at(i) != 1 || tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pt->at(i) < 0.0 )
            continue;

        Object obj = Object( tpTo_hltIter0Phase2L3FromL1TkMuon_pt->at(i), tpTo_hltIter0Phase2L3FromL1TkMuon_eta->at(i), tpTo_hltIter0Phase2L3FromL1TkMuon_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIter0Phase2L3FromL1TkMuon_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIter0Phase2L3FromL1TkMuon_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIter0Phase2L3FromL1TkMuon_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIter0Phase2L3FromL1TkMuon_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIter0Phase2L3FromL1TkMuon_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIter0Phase2L3FromL1TkMuon_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIter0Phase2L3FromL1TkMuon_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIter0Phase2L3FromL1TkMuon_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIter0Phase2L3FromL1TkMuon_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIter0Phase2L3FromL1TkMuon_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIter0Phase2L3FromL1TkMuon_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIter0Phase2L3FromL1TkMuon_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIter0Phase2L3FromL1TkMuon_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIter2Phase2L3FromL1TkMuon()
{
    vector<Object> out = {};
    if(tpTo_hltIter2Phase2L3FromL1TkMuon_pt == 0 || tpTo_hltIter2Phase2L3FromL1TkMuon_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIter2Phase2L3FromL1TkMuon_pt->size(); ++i) {

        if( tpTo_hltIter2Phase2L3FromL1TkMuon_pt->at(i) < 5.0 || tpTo_hltIter2Phase2L3FromL1TkMuon_status->at(i) != 1 || tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pt->at(i) < 0.0 )
            continue;

        Object obj = Object( tpTo_hltIter2Phase2L3FromL1TkMuon_pt->at(i), tpTo_hltIter2Phase2L3FromL1TkMuon_eta->at(i), tpTo_hltIter2Phase2L3FromL1TkMuon_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIter2Phase2L3FromL1TkMuon_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIter2Phase2L3FromL1TkMuon_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIter2Phase2L3FromL1TkMuon_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIter2Phase2L3FromL1TkMuon_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIter2Phase2L3FromL1TkMuon_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIter2Phase2L3FromL1TkMuon_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIter2Phase2L3FromL1TkMuon_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIter2Phase2L3FromL1TkMuon_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIter2Phase2L3FromL1TkMuon_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIter2Phase2L3FromL1TkMuon_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIter2Phase2L3FromL1TkMuon_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIter2Phase2L3FromL1TkMuon_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIter2Phase2L3FromL1TkMuon_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltPhase2L3IOFromL1()
{
    vector<Object> out = {};
    if(tpTo_hltPhase2L3IOFromL1_pt == 0 || tpTo_hltPhase2L3IOFromL1_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltPhase2L3IOFromL1_pt->size(); ++i) {

        if( tpTo_hltPhase2L3IOFromL1_pt->at(i) < 5.0 || tpTo_hltPhase2L3IOFromL1_status->at(i) != 1 || tpTo_hltPhase2L3IOFromL1_gen_pt->at(i) < 0.0 )
            continue;

        Object obj = Object( tpTo_hltPhase2L3IOFromL1_pt->at(i), tpTo_hltPhase2L3IOFromL1_eta->at(i), tpTo_hltPhase2L3IOFromL1_phi->at(i) );

        obj.addVar( "charge", tpTo_hltPhase2L3IOFromL1_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltPhase2L3IOFromL1_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltPhase2L3IOFromL1_energy->at(i) );
        obj.addVar( "pt", tpTo_hltPhase2L3IOFromL1_pt->at(i) );
        obj.addVar( "eta", tpTo_hltPhase2L3IOFromL1_eta->at(i) );
        obj.addVar( "phi", tpTo_hltPhase2L3IOFromL1_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltPhase2L3IOFromL1_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltPhase2L3IOFromL1_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltPhase2L3IOFromL1_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltPhase2L3IOFromL1_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltPhase2L3IOFromL1_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltPhase2L3IOFromL1_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltPhase2L3IOFromL1_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltPhase2L3IOFromL1_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltPhase2L3IOFromL1_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltPhase2L3IOFromL1_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltPhase2L3IOFromL1_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltPhase2L3IOFromL1_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltPhase2L3IOFromL1_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltPhase2L3IOFromL1_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltPhase2L3IOFromL1_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltPhase2L3IOFromL1_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltPhase2L3IOFromL1_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltPhase2L3IOFromL1_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltPhase2L3MuonsNoID()
{
    vector<Object> out = {};
    if(tpTo_hltPhase2L3MuonsNoID_pt == 0 || tpTo_hltPhase2L3MuonsNoID_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltPhase2L3MuonsNoID_pt->size(); ++i) {

        if( tpTo_hltPhase2L3MuonsNoID_pt->at(i) < 5.0 || tpTo_hltPhase2L3MuonsNoID_status->at(i) != 1 || tpTo_hltPhase2L3MuonsNoID_gen_pt->at(i) < 0.0 )
            continue;

        Object obj = Object( tpTo_hltPhase2L3MuonsNoID_pt->at(i), tpTo_hltPhase2L3MuonsNoID_eta->at(i), tpTo_hltPhase2L3MuonsNoID_phi->at(i) );

        obj.addVar( "charge", tpTo_hltPhase2L3MuonsNoID_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltPhase2L3MuonsNoID_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltPhase2L3MuonsNoID_energy->at(i) );
        obj.addVar( "pt", tpTo_hltPhase2L3MuonsNoID_pt->at(i) );
        obj.addVar( "eta", tpTo_hltPhase2L3MuonsNoID_eta->at(i) );
        obj.addVar( "phi", tpTo_hltPhase2L3MuonsNoID_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltPhase2L3MuonsNoID_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltPhase2L3MuonsNoID_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltPhase2L3MuonsNoID_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltPhase2L3MuonsNoID_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltPhase2L3MuonsNoID_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltPhase2L3MuonsNoID_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltPhase2L3MuonsNoID_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltPhase2L3MuonsNoID_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltPhase2L3MuonsNoID_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltPhase2L3MuonsNoID_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltPhase2L3MuonsNoID_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltPhase2L3MuonsNoID_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltPhase2L3Muons()
{
    vector<Object> out = {};
    if(tpTo_hltPhase2L3Muons_pt == 0 || tpTo_hltPhase2L3Muons_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltPhase2L3Muons_pt->size(); ++i) {

        if( tpTo_hltPhase2L3Muons_pt->at(i) < 5.0 || tpTo_hltPhase2L3Muons_status->at(i) != 1 || tpTo_hltPhase2L3Muons_gen_pt->at(i) < 0.0 )
            continue;

        Object obj = Object( tpTo_hltPhase2L3Muons_pt->at(i), tpTo_hltPhase2L3Muons_eta->at(i), tpTo_hltPhase2L3Muons_phi->at(i) );

        obj.addVar( "charge", tpTo_hltPhase2L3Muons_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltPhase2L3Muons_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltPhase2L3Muons_energy->at(i) );
        obj.addVar( "pt", tpTo_hltPhase2L3Muons_pt->at(i) );
        obj.addVar( "eta", tpTo_hltPhase2L3Muons_eta->at(i) );
        obj.addVar( "phi", tpTo_hltPhase2L3Muons_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltPhase2L3Muons_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltPhase2L3Muons_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltPhase2L3Muons_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltPhase2L3Muons_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltPhase2L3Muons_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltPhase2L3Muons_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltPhase2L3Muons_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltPhase2L3Muons_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltPhase2L3Muons_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltPhase2L3Muons_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltPhase2L3Muons_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltPhase2L3Muons_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltPhase2L3Muons_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltPhase2L3Muons_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltPhase2L3Muons_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltPhase2L3Muons_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltPhase2L3Muons_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltPhase2L3Muons_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIterL3OI()
{
    vector<Object> out = {};
    if(hltIterL3OI_pt == 0 || hltIterL3OI_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIterL3OI_pt->size(); ++i) {
        Object obj = Object( hltIterL3OI_pt->at(i), hltIterL3OI_eta->at(i), hltIterL3OI_phi->at(i) );

        obj.addVar( "pt", hltIterL3OI_pt->at(i) );
        obj.addVar( "ptError", hltIterL3OI_ptError->at(i) );
        obj.addVar( "eta", hltIterL3OI_eta->at(i) );
        obj.addVar( "phi", hltIterL3OI_phi->at(i) );
        obj.addVar( "charge", hltIterL3OI_charge->at(i) );
        obj.addVar( "matchedL3", hltIterL3OI_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIterL3OI_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIterL3OI_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIterL3OI_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIterL3OI_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIterL3OI_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIterL3OI_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIterL3OI_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIterL3OI_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIterL3OI_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIterL3OI_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIterL3OI_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIterL3OI_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIterL3OI_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIterL3OI_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIterL3OI_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIterL3OI_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIterL3OI_mva0->at(i) );
        obj.addVar( "mva1", hltIterL3OI_mva1->at(i) );
        obj.addVar( "mva2", hltIterL3OI_mva2->at(i) );
        obj.addVar( "mva3", hltIterL3OI_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIter0IterL3()
{
    vector<Object> out = {};
    if(hltIter0IterL3_pt == 0 || hltIter0IterL3_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIter0IterL3_pt->size(); ++i) {
        Object obj = Object( hltIter0IterL3_pt->at(i), hltIter0IterL3_eta->at(i), hltIter0IterL3_phi->at(i) );

        obj.addVar( "pt", hltIter0IterL3_pt->at(i) );
        obj.addVar( "ptError", hltIter0IterL3_ptError->at(i) );
        obj.addVar( "eta", hltIter0IterL3_eta->at(i) );
        obj.addVar( "phi", hltIter0IterL3_phi->at(i) );
        obj.addVar( "charge", hltIter0IterL3_charge->at(i) );
        obj.addVar( "matchedL3", hltIter0IterL3_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIter0IterL3_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIter0IterL3_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIter0IterL3_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIter0IterL3_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIter0IterL3_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIter0IterL3_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIter0IterL3_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIter0IterL3_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIter0IterL3_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIter0IterL3_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIter0IterL3_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIter0IterL3_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIter0IterL3_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIter0IterL3_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIter0IterL3_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIter0IterL3_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIter0IterL3_mva0->at(i) );
        obj.addVar( "mva1", hltIter0IterL3_mva1->at(i) );
        obj.addVar( "mva2", hltIter0IterL3_mva2->at(i) );
        obj.addVar( "mva3", hltIter0IterL3_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIter2IterL3()
{
    vector<Object> out = {};
    if(hltIter2IterL3_pt == 0 || hltIter2IterL3_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIter2IterL3_pt->size(); ++i) {
        Object obj = Object( hltIter2IterL3_pt->at(i), hltIter2IterL3_eta->at(i), hltIter2IterL3_phi->at(i) );

        obj.addVar( "pt", hltIter2IterL3_pt->at(i) );
        obj.addVar( "ptError", hltIter2IterL3_ptError->at(i) );
        obj.addVar( "eta", hltIter2IterL3_eta->at(i) );
        obj.addVar( "phi", hltIter2IterL3_phi->at(i) );
        obj.addVar( "charge", hltIter2IterL3_charge->at(i) );
        obj.addVar( "matchedL3", hltIter2IterL3_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIter2IterL3_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIter2IterL3_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIter2IterL3_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIter2IterL3_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIter2IterL3_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIter2IterL3_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIter2IterL3_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIter2IterL3_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIter2IterL3_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIter2IterL3_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIter2IterL3_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIter2IterL3_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIter2IterL3_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIter2IterL3_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIter2IterL3_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIter2IterL3_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIter2IterL3_mva0->at(i) );
        obj.addVar( "mva1", hltIter2IterL3_mva1->at(i) );
        obj.addVar( "mva2", hltIter2IterL3_mva2->at(i) );
        obj.addVar( "mva3", hltIter2IterL3_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIter0IterL3FromL1Muon()
{
    vector<Object> out = {};
    if(hltIter0IterL3FromL1Muon_pt == 0 || hltIter0IterL3FromL1Muon_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIter0IterL3FromL1Muon_pt->size(); ++i) {
        Object obj = Object( hltIter0IterL3FromL1Muon_pt->at(i), hltIter0IterL3FromL1Muon_eta->at(i), hltIter0IterL3FromL1Muon_phi->at(i) );

        obj.addVar( "pt", hltIter0IterL3FromL1Muon_pt->at(i) );
        obj.addVar( "ptError", hltIter0IterL3FromL1Muon_ptError->at(i) );
        obj.addVar( "eta", hltIter0IterL3FromL1Muon_eta->at(i) );
        obj.addVar( "phi", hltIter0IterL3FromL1Muon_phi->at(i) );
        obj.addVar( "charge", hltIter0IterL3FromL1Muon_charge->at(i) );
        obj.addVar( "matchedL3", hltIter0IterL3FromL1Muon_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIter0IterL3FromL1Muon_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIter0IterL3FromL1Muon_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIter0IterL3FromL1Muon_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIter0IterL3FromL1Muon_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIter0IterL3FromL1Muon_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIter0IterL3FromL1Muon_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIter0IterL3FromL1Muon_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIter0IterL3FromL1Muon_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIter0IterL3FromL1Muon_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIter0IterL3FromL1Muon_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIter0IterL3FromL1Muon_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIter0IterL3FromL1Muon_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIter0IterL3FromL1Muon_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIter0IterL3FromL1Muon_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIter0IterL3FromL1Muon_mva0->at(i) );
        obj.addVar( "mva1", hltIter0IterL3FromL1Muon_mva1->at(i) );
        obj.addVar( "mva2", hltIter0IterL3FromL1Muon_mva2->at(i) );
        obj.addVar( "mva3", hltIter0IterL3FromL1Muon_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIter2IterL3FromL1Muon()
{
    vector<Object> out = {};
    if(hltIter2IterL3FromL1Muon_pt == 0 || hltIter2IterL3FromL1Muon_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIter2IterL3FromL1Muon_pt->size(); ++i) {
        Object obj = Object( hltIter2IterL3FromL1Muon_pt->at(i), hltIter2IterL3FromL1Muon_eta->at(i), hltIter2IterL3FromL1Muon_phi->at(i) );

        obj.addVar( "pt", hltIter2IterL3FromL1Muon_pt->at(i) );
        obj.addVar( "ptError", hltIter2IterL3FromL1Muon_ptError->at(i) );
        obj.addVar( "eta", hltIter2IterL3FromL1Muon_eta->at(i) );
        obj.addVar( "phi", hltIter2IterL3FromL1Muon_phi->at(i) );
        obj.addVar( "charge", hltIter2IterL3FromL1Muon_charge->at(i) );
        obj.addVar( "matchedL3", hltIter2IterL3FromL1Muon_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIter2IterL3FromL1Muon_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIter2IterL3FromL1Muon_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIter2IterL3FromL1Muon_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIter2IterL3FromL1Muon_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIter2IterL3FromL1Muon_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIter2IterL3FromL1Muon_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIter2IterL3FromL1Muon_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIter2IterL3FromL1Muon_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIter2IterL3FromL1Muon_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIter2IterL3FromL1Muon_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIter2IterL3FromL1Muon_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIter2IterL3FromL1Muon_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIter2IterL3FromL1Muon_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIter2IterL3FromL1Muon_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIter2IterL3FromL1Muon_mva0->at(i) );
        obj.addVar( "mva1", hltIter2IterL3FromL1Muon_mva1->at(i) );
        obj.addVar( "mva2", hltIter2IterL3FromL1Muon_mva2->at(i) );
        obj.addVar( "mva3", hltIter2IterL3FromL1Muon_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIter2IterL3FromL1MuonTrack()
{
    vector<Object> out = {};
    if(hltIter2IterL3FromL1MuonTrack_pt == 0 || hltIter2IterL3FromL1MuonTrack_pt == nullptr)
        return out;
    
    for(unsigned i=0; i<hltIter2IterL3FromL1MuonTrack_pt->size(); ++i) {
        Object obj = Object( hltIter2IterL3FromL1MuonTrack_pt->at(i), hltIter2IterL3FromL1MuonTrack_eta->at(i), hltIter2IterL3FromL1MuonTrack_phi->at(i) );

        obj.addVar( "pt", hltIter2IterL3FromL1MuonTrack_pt->at(i) );
        obj.addVar( "ptError", hltIter2IterL3FromL1MuonTrack_ptError->at(i) );
        obj.addVar( "eta", hltIter2IterL3FromL1MuonTrack_eta->at(i) );
        obj.addVar( "phi", hltIter2IterL3FromL1MuonTrack_phi->at(i) );
        obj.addVar( "charge", hltIter2IterL3FromL1MuonTrack_charge->at(i) );
        obj.addVar( "matchedL3", hltIter2IterL3FromL1MuonTrack_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIter2IterL3FromL1MuonTrack_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIter2IterL3FromL1MuonTrack_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIter2IterL3FromL1MuonTrack_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIter2IterL3FromL1MuonTrack_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIter2IterL3FromL1MuonTrack_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIter2IterL3FromL1MuonTrack_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIter2IterL3FromL1MuonTrack_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIter2IterL3FromL1MuonTrack_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIter2IterL3FromL1MuonTrack_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIter2IterL3FromL1MuonTrack_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIter2IterL3FromL1MuonTrack_mva0->at(i) );
        obj.addVar( "mva1", hltIter2IterL3FromL1MuonTrack_mva1->at(i) );
        obj.addVar( "mva2", hltIter2IterL3FromL1MuonTrack_mva2->at(i) );
        obj.addVar( "mva3", hltIter2IterL3FromL1MuonTrack_mva3->at(i) );

        out.push_back( obj );
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIterL3IOFromL1()
{
    vector<Object> out = {};
    if(hltIterL3IOFromL1_pt == 0 || hltIterL3IOFromL1_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIterL3IOFromL1_pt->size(); ++i) {
        Object obj = Object( hltIterL3IOFromL1_pt->at(i), hltIterL3IOFromL1_eta->at(i), hltIterL3IOFromL1_phi->at(i) );

        obj.addVar( "pt", hltIterL3IOFromL1_pt->at(i) );
        obj.addVar( "ptError", hltIterL3IOFromL1_ptError->at(i) );
        obj.addVar( "eta", hltIterL3IOFromL1_eta->at(i) );
        obj.addVar( "phi", hltIterL3IOFromL1_phi->at(i) );
        obj.addVar( "charge", hltIterL3IOFromL1_charge->at(i) );
        obj.addVar( "matchedL3", hltIterL3IOFromL1_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIterL3IOFromL1_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIterL3IOFromL1_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIterL3IOFromL1_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIterL3IOFromL1_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIterL3IOFromL1_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIterL3IOFromL1_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIterL3IOFromL1_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIterL3IOFromL1_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIterL3IOFromL1_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIterL3IOFromL1_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIterL3IOFromL1_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIterL3IOFromL1_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIterL3IOFromL1_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIterL3IOFromL1_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIterL3IOFromL1_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIterL3IOFromL1_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIterL3IOFromL1_mva0->at(i) );
        obj.addVar( "mva1", hltIterL3IOFromL1_mva1->at(i) );
        obj.addVar( "mva2", hltIterL3IOFromL1_mva2->at(i) );
        obj.addVar( "mva3", hltIterL3IOFromL1_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIterL3MuonsNoID()
{
    vector<Object> out = {};
    if(hltIterL3MuonsNoID_pt == 0 || hltIterL3MuonsNoID_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIterL3MuonsNoID_pt->size(); ++i) {
        Object obj = Object( hltIterL3MuonsNoID_pt->at(i), hltIterL3MuonsNoID_eta->at(i), hltIterL3MuonsNoID_phi->at(i) );

        obj.addVar( "pt", hltIterL3MuonsNoID_pt->at(i) );
        obj.addVar( "ptError", hltIterL3MuonsNoID_ptError->at(i) );
        obj.addVar( "eta", hltIterL3MuonsNoID_eta->at(i) );
        obj.addVar( "phi", hltIterL3MuonsNoID_phi->at(i) );
        obj.addVar( "charge", hltIterL3MuonsNoID_charge->at(i) );
        obj.addVar( "matchedL3", hltIterL3MuonsNoID_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIterL3MuonsNoID_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIterL3MuonsNoID_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIterL3MuonsNoID_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIterL3MuonsNoID_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIterL3MuonsNoID_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIterL3MuonsNoID_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIterL3MuonsNoID_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIterL3MuonsNoID_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIterL3MuonsNoID_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIterL3MuonsNoID_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIterL3MuonsNoID_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIterL3MuonsNoID_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIterL3MuonsNoID_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIterL3MuonsNoID_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIterL3MuonsNoID_mva0->at(i) );
        obj.addVar( "mva1", hltIterL3MuonsNoID_mva1->at(i) );
        obj.addVar( "mva2", hltIterL3MuonsNoID_mva2->at(i) );
        obj.addVar( "mva3", hltIterL3MuonsNoID_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_hltIterL3Muons()
{
    vector<Object> out = {};
    if(hltIterL3Muons_pt == 0 || hltIterL3Muons_pt == nullptr)
        return out;

    for(unsigned i=0; i<hltIterL3Muons_pt->size(); ++i) {
        Object obj = Object( hltIterL3Muons_pt->at(i), hltIterL3Muons_eta->at(i), hltIterL3Muons_phi->at(i) );

        obj.addVar( "pt", hltIterL3Muons_pt->at(i) );
        obj.addVar( "ptError", hltIterL3Muons_ptError->at(i) );
        obj.addVar( "eta", hltIterL3Muons_eta->at(i) );
        obj.addVar( "phi", hltIterL3Muons_phi->at(i) );
        obj.addVar( "charge", hltIterL3Muons_charge->at(i) );
        obj.addVar( "matchedL3", hltIterL3Muons_matchedL3->at(i) );
        obj.addVar( "matchedL3NoId", hltIterL3Muons_matchedL3NoId->at(i) );
        obj.addVar( "bestMatchTP_charge", hltIterL3Muons_bestMatchTP_charge->at(i) );
        obj.addVar( "bestMatchTP_pdgId", hltIterL3Muons_bestMatchTP_pdgId->at(i) );
        obj.addVar( "bestMatchTP_energy", hltIterL3Muons_bestMatchTP_energy->at(i) );
        obj.addVar( "bestMatchTP_pt", hltIterL3Muons_bestMatchTP_pt->at(i) );
        obj.addVar( "bestMatchTP_eta", hltIterL3Muons_bestMatchTP_eta->at(i) );
        obj.addVar( "bestMatchTP_phi", hltIterL3Muons_bestMatchTP_phi->at(i) );
        obj.addVar( "bestMatchTP_parentVx", hltIterL3Muons_bestMatchTP_parentVx->at(i) );
        obj.addVar( "bestMatchTP_parentVy", hltIterL3Muons_bestMatchTP_parentVy->at(i) );
        obj.addVar( "bestMatchTP_parentVz", hltIterL3Muons_bestMatchTP_parentVz->at(i) );
        obj.addVar( "bestMatchTP_status", hltIterL3Muons_bestMatchTP_status->at(i) );
        obj.addVar( "bestMatchTP_numberOfHits", hltIterL3Muons_bestMatchTP_numberOfHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerHits", hltIterL3Muons_bestMatchTP_numberOfTrackerHits->at(i) );
        obj.addVar( "bestMatchTP_numberOfTrackerLayers", hltIterL3Muons_bestMatchTP_numberOfTrackerLayers->at(i) );
        obj.addVar( "bestMatchTP_sharedFraction", hltIterL3Muons_bestMatchTP_sharedFraction->at(i) );
        obj.addVar( "matchedTPsize", hltIterL3Muons_matchedTPsize->at(i) );
        obj.addVar( "mva0", hltIterL3Muons_mva0->at(i) );
        obj.addVar( "mva1", hltIterL3Muons_mva1->at(i) );
        obj.addVar( "mva2", hltIterL3Muons_mva2->at(i) );
        obj.addVar( "mva3", hltIterL3Muons_mva3->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIterL3OI()
{
    vector<Object> out = {};
    if(tpTo_hltIterL3OI_pt == 0 || tpTo_hltIterL3OI_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIterL3OI_pt->size(); ++i) {
        Object obj = Object( tpTo_hltIterL3OI_pt->at(i), tpTo_hltIterL3OI_eta->at(i), tpTo_hltIterL3OI_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIterL3OI_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIterL3OI_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIterL3OI_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIterL3OI_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIterL3OI_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIterL3OI_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIterL3OI_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIterL3OI_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIterL3OI_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIterL3OI_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIterL3OI_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIterL3OI_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIterL3OI_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIterL3OI_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIterL3OI_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIterL3OI_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIterL3OI_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIterL3OI_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIterL3OI_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIterL3OI_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIterL3OI_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIterL3OI_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIterL3OI_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIterL3OI_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIter0IterL3()
{
    vector<Object> out = {};
    if(tpTo_hltIter0IterL3_pt == 0 || tpTo_hltIter0IterL3_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIter0IterL3_pt->size(); ++i) {
        Object obj = Object( tpTo_hltIter0IterL3_pt->at(i), tpTo_hltIter0IterL3_eta->at(i), tpTo_hltIter0IterL3_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIter0IterL3_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIter0IterL3_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIter0IterL3_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIter0IterL3_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIter0IterL3_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIter0IterL3_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIter0IterL3_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIter0IterL3_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIter0IterL3_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIter0IterL3_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIter0IterL3_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIter0IterL3_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIter0IterL3_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIter0IterL3_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIter0IterL3_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIter0IterL3_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIter0IterL3_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIter0IterL3_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIter0IterL3_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIter0IterL3_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIter0IterL3_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIter0IterL3_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIter0IterL3_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIter0IterL3_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIter2IterL3()
{
    vector<Object> out = {};
    if(tpTo_hltIter2IterL3_pt == 0 || tpTo_hltIter2IterL3_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIter2IterL3_pt->size(); ++i) {
        Object obj = Object( tpTo_hltIter2IterL3_pt->at(i), tpTo_hltIter2IterL3_eta->at(i), tpTo_hltIter2IterL3_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIter2IterL3_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIter2IterL3_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIter2IterL3_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIter2IterL3_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIter2IterL3_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIter2IterL3_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIter2IterL3_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIter2IterL3_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIter2IterL3_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIter2IterL3_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIter2IterL3_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIter2IterL3_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIter2IterL3_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIter2IterL3_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIter2IterL3_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIter2IterL3_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIter2IterL3_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIter2IterL3_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIter2IterL3_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIter2IterL3_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIter2IterL3_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIter2IterL3_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIter2IterL3_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIter2IterL3_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIter0IterL3FromL1Muon()
{
    vector<Object> out = {};
    if(tpTo_hltIter0IterL3FromL1Muon_pt == 0 || tpTo_hltIter0IterL3FromL1Muon_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIter0IterL3FromL1Muon_pt->size(); ++i) {
        Object obj = Object( tpTo_hltIter0IterL3FromL1Muon_pt->at(i), tpTo_hltIter0IterL3FromL1Muon_eta->at(i), tpTo_hltIter0IterL3FromL1Muon_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIter0IterL3FromL1Muon_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIter0IterL3FromL1Muon_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIter0IterL3FromL1Muon_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIter0IterL3FromL1Muon_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIter0IterL3FromL1Muon_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIter0IterL3FromL1Muon_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIter0IterL3FromL1Muon_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIter0IterL3FromL1Muon_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIter0IterL3FromL1Muon_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIter0IterL3FromL1Muon_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIter0IterL3FromL1Muon_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIter0IterL3FromL1Muon_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIter0IterL3FromL1Muon_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIter0IterL3FromL1Muon_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIter0IterL3FromL1Muon_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIter0IterL3FromL1Muon_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIter2IterL3FromL1Muon()
{
    vector<Object> out = {};
    if(tpTo_hltIter2IterL3FromL1Muon_pt == 0 || tpTo_hltIter2IterL3FromL1Muon_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIter2IterL3FromL1Muon_pt->size(); ++i) {
        Object obj = Object( tpTo_hltIter2IterL3FromL1Muon_pt->at(i), tpTo_hltIter2IterL3FromL1Muon_eta->at(i), tpTo_hltIter2IterL3FromL1Muon_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIter2IterL3FromL1Muon_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIter2IterL3FromL1Muon_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIter2IterL3FromL1Muon_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIter2IterL3FromL1Muon_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIter2IterL3FromL1Muon_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIter2IterL3FromL1Muon_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIter2IterL3FromL1Muon_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIter2IterL3FromL1Muon_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIter2IterL3FromL1Muon_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIter2IterL3FromL1Muon_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIter2IterL3FromL1Muon_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIter2IterL3FromL1Muon_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIter2IterL3FromL1Muon_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIter2IterL3FromL1Muon_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIter2IterL3FromL1Muon_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIter2IterL3FromL1Muon_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIterL3IOFromL1()
{
    vector<Object> out = {};
    if(tpTo_hltIterL3IOFromL1_pt == 0 || tpTo_hltIterL3IOFromL1_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIterL3IOFromL1_pt->size(); ++i) {
        Object obj = Object( tpTo_hltIterL3IOFromL1_pt->at(i), tpTo_hltIterL3IOFromL1_eta->at(i), tpTo_hltIterL3IOFromL1_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIterL3IOFromL1_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIterL3IOFromL1_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIterL3IOFromL1_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIterL3IOFromL1_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIterL3IOFromL1_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIterL3IOFromL1_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIterL3IOFromL1_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIterL3IOFromL1_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIterL3IOFromL1_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIterL3IOFromL1_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIterL3IOFromL1_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIterL3IOFromL1_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIterL3IOFromL1_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIterL3IOFromL1_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIterL3IOFromL1_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIterL3IOFromL1_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIterL3IOFromL1_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIterL3IOFromL1_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIterL3IOFromL1_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIterL3IOFromL1_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIterL3IOFromL1_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIterL3IOFromL1_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIterL3IOFromL1_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIterL3IOFromL1_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIterL3MuonsNoID()
{
    vector<Object> out = {};
    if(tpTo_hltIterL3MuonsNoID_pt == 0 || tpTo_hltIterL3MuonsNoID_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIterL3MuonsNoID_pt->size(); ++i) {
        Object obj = Object( tpTo_hltIterL3MuonsNoID_pt->at(i), tpTo_hltIterL3MuonsNoID_eta->at(i), tpTo_hltIterL3MuonsNoID_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIterL3MuonsNoID_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIterL3MuonsNoID_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIterL3MuonsNoID_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIterL3MuonsNoID_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIterL3MuonsNoID_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIterL3MuonsNoID_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIterL3MuonsNoID_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIterL3MuonsNoID_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIterL3MuonsNoID_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIterL3MuonsNoID_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIterL3MuonsNoID_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIterL3MuonsNoID_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIterL3MuonsNoID_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIterL3MuonsNoID_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIterL3MuonsNoID_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIterL3MuonsNoID_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIterL3MuonsNoID_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIterL3MuonsNoID_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIterL3MuonsNoID_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIterL3MuonsNoID_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIterL3MuonsNoID_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIterL3MuonsNoID_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIterL3MuonsNoID_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIterL3MuonsNoID_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}

vector<Object> MuonHLTNtuple::get_tpTo_hltIterL3Muons()
{
    vector<Object> out = {};
    if(tpTo_hltIterL3Muons_pt == 0 || tpTo_hltIterL3Muons_pt == nullptr)
        return out;

    for(unsigned i=0; i<tpTo_hltIterL3Muons_pt->size(); ++i) {
        Object obj = Object( tpTo_hltIterL3Muons_pt->at(i), tpTo_hltIterL3Muons_eta->at(i), tpTo_hltIterL3Muons_phi->at(i) );

        obj.addVar( "charge", tpTo_hltIterL3Muons_charge->at(i) );
        obj.addVar( "pdgId", tpTo_hltIterL3Muons_pdgId->at(i) );
        obj.addVar( "energy", tpTo_hltIterL3Muons_energy->at(i) );
        obj.addVar( "pt", tpTo_hltIterL3Muons_pt->at(i) );
        obj.addVar( "eta", tpTo_hltIterL3Muons_eta->at(i) );
        obj.addVar( "phi", tpTo_hltIterL3Muons_phi->at(i) );
        obj.addVar( "parentVx", tpTo_hltIterL3Muons_parentVx->at(i) );
        obj.addVar( "parentVy", tpTo_hltIterL3Muons_parentVy->at(i) );
        obj.addVar( "parentVz", tpTo_hltIterL3Muons_parentVz->at(i) );
        obj.addVar( "status", tpTo_hltIterL3Muons_status->at(i) );
        obj.addVar( "numberOfHits", tpTo_hltIterL3Muons_numberOfHits->at(i) );
        obj.addVar( "numberOfTrackerHits", tpTo_hltIterL3Muons_numberOfTrackerHits->at(i) );
        obj.addVar( "numberOfTrackerLayers", tpTo_hltIterL3Muons_numberOfTrackerLayers->at(i) );
        obj.addVar( "gen_charge", tpTo_hltIterL3Muons_gen_charge->at(i) );
        obj.addVar( "gen_pdgId", tpTo_hltIterL3Muons_gen_pdgId->at(i) );
        obj.addVar( "gen_pt", tpTo_hltIterL3Muons_gen_pt->at(i) );
        obj.addVar( "gen_eta", tpTo_hltIterL3Muons_gen_eta->at(i) );
        obj.addVar( "gen_phi", tpTo_hltIterL3Muons_gen_phi->at(i) );
        obj.addVar( "bestMatchTrk_pt", tpTo_hltIterL3Muons_bestMatchTrk_pt->at(i) );
        obj.addVar( "bestMatchTrk_eta", tpTo_hltIterL3Muons_bestMatchTrk_eta->at(i) );
        obj.addVar( "bestMatchTrk_phi", tpTo_hltIterL3Muons_bestMatchTrk_phi->at(i) );
        obj.addVar( "bestMatchTrk_charge", tpTo_hltIterL3Muons_bestMatchTrk_charge->at(i) );
        obj.addVar( "bestMatchTrk_quality", tpTo_hltIterL3Muons_bestMatchTrk_quality->at(i) );
        obj.addVar( "bestMatchTrk_NValidHits", tpTo_hltIterL3Muons_bestMatchTrk_NValidHits->at(i) );

        out.push_back(obj);
    }

    return out;
}



Int_t MuonHLTNtuple::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t MuonHLTNtuple::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
       fCurrent = fChain->GetTreeNumber();
       Notify();
    }
    return centry;
}

void MuonHLTNtuple::Init(TChain *tree)
{

    // -- Initialize -- #
        trk_pt = 0;
        trk_eta = 0;
        trk_phi = 0;
        trk_d0 = 0;
        trk_z0 = 0;
        trk_rInv = 0;
        trk_tanL = 0;
        trk_MVA1 = 0;
        trk_MVA2 = 0;
        trk_MVA3 = 0;
        trk_chi2 = 0;
        trk_bendchi2 = 0;
        trk_nstub = 0;
        trk_lhits = 0;
        trk_dhits = 0;
        trk_seed = 0;
        trk_phiSector = 0;
        trk_genuine = 0;
        trk_loose = 0;
        trk_unknown = 0;
        trk_combinatoric = 0;
        trk_fake = 0;
        trk_matchtp_pdgid = 0;
        trk_matchtp_pt = 0;
        trk_matchtp_eta = 0;
        trk_matchtp_phi = 0;
        trk_matchtp_z0 = 0;
        trk_matchtp_dxy = 0;
        stub_x = 0;
        stub_y = 0;
        stub_z = 0;
        stub_isBarrel = 0;
        stub_layer = 0;
        L1TkMu_pt = 0;
        L1TkMu_eta = 0;
        L1TkMu_phi = 0;
        L1TkMu_trkIsol = 0;
        L1TkMu_trkzVtx = 0;
        L1TkMu_dR = 0;
        L1TkMu_nTracksMatched = 0;
        L1TkMu_trackCurvature = 0;
        L1TkMu_quality = 0;
        L1TkMu_pattern = 0;
        L1TkMu_muonDetector = 0;
        L1TkMu_TTTpointer = 0;
        L1TkMu_muRefHwPt = 0;
        L1TkMu_muRefHwDXY = 0;
        L1TkMu_muRefHwEta = 0;
        L1TkMu_muRefHwPhi = 0;
        L1TkMu_muRefHwSign = 0;
        L1TkMu_muRefHwSignValid = 0;
        L1TkMu_muRefHwQual = 0;
        vec_firedTrigger = 0;
        vec_filterName = 0;
        vec_HLTObj_pt = 0;
        vec_HLTObj_eta = 0;
        vec_HLTObj_phi = 0;
        vec_myFiredTrigger = 0;
        vec_myFilterName = 0;
        vec_myHLTObj_pt = 0;
        vec_myHLTObj_eta = 0;
        vec_myHLTObj_phi = 0;
        hltIterL3OISeedsFromL2Muons_dir = 0;
        hltIterL3OISeedsFromL2Muons_tsos_detId = 0;
        hltIterL3OISeedsFromL2Muons_tsos_pt = 0;
        hltIterL3OISeedsFromL2Muons_tsos_pt_val = 0;
        hltIterL3OISeedsFromL2Muons_tsos_eta = 0;
        hltIterL3OISeedsFromL2Muons_tsos_phi = 0;
        hltIterL3OISeedsFromL2Muons_tsos_glob_x = 0;
        hltIterL3OISeedsFromL2Muons_tsos_glob_y = 0;
        hltIterL3OISeedsFromL2Muons_tsos_glob_z = 0;
        hltIterL3OISeedsFromL2Muons_tsos_hasErr = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err0 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err1 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err2 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err3 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err4 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err5 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err6 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err7 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err8 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err9 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err10 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err11 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err12 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err13 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_err14 = 0;
        hltIterL3OISeedsFromL2Muons_tsos_x = 0;
        hltIterL3OISeedsFromL2Muons_tsos_y = 0;
        hltIterL3OISeedsFromL2Muons_tsos_dxdz = 0;
        hltIterL3OISeedsFromL2Muons_tsos_dydz = 0;
        hltIterL3OISeedsFromL2Muons_tsos_px = 0;
        hltIterL3OISeedsFromL2Muons_tsos_py = 0;
        hltIterL3OISeedsFromL2Muons_tsos_pz = 0;
        hltIterL3OISeedsFromL2Muons_tsos_qbp = 0;
        hltIterL3OISeedsFromL2Muons_tsos_charge = 0;
        hltIterL3OISeedsFromL2Muons_iterL3Matched = 0;
        hltIterL3OISeedsFromL2Muons_iterL3Ref = 0;
        hltIterL3OISeedsFromL2Muons_tmpL3Ref = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_dir = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_detId = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt_val = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_eta = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_phi = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_x = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_y = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_z = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_hasErr = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err0 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err1 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err2 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err3 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err4 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err5 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err6 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err7 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err8 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err9 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err10 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err11 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err12 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err13 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err14 = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_x = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_y = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dxdz = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dydz = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_px = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_py = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pz = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_qbp = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_charge = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Matched = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Ref = 0;
        hltIter0IterL3MuonPixelSeedsFromPixelTracks_tmpL3Ref = 0;
        hltIter2IterL3MuonPixelSeeds_dir = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_detId = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_pt = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_pt_val = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_eta = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_phi = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_glob_x = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_glob_y = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_glob_z = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_hasErr = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err0 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err1 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err2 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err3 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err4 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err5 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err6 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err7 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err8 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err9 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err10 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err11 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err12 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err13 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_err14 = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_x = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_y = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_dxdz = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_dydz = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_px = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_py = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_pz = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_qbp = 0;
        hltIter2IterL3MuonPixelSeeds_tsos_charge = 0;
        hltIter2IterL3MuonPixelSeeds_iterL3Matched = 0;
        hltIter2IterL3MuonPixelSeeds_iterL3Ref = 0;
        hltIter2IterL3MuonPixelSeeds_tmpL3Ref = 0;
        hltIter3IterL3MuonPixelSeeds_dir = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_detId = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_pt = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_pt_val = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_eta = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_phi = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_glob_x = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_glob_y = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_glob_z = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_hasErr = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err0 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err1 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err2 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err3 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err4 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err5 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err6 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err7 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err8 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err9 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err10 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err11 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err12 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err13 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_err14 = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_x = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_y = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_dxdz = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_dydz = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_px = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_py = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_pz = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_qbp = 0;
        hltIter3IterL3MuonPixelSeeds_tsos_charge = 0;
        hltIter3IterL3MuonPixelSeeds_iterL3Matched = 0;
        hltIter3IterL3MuonPixelSeeds_iterL3Ref = 0;
        hltIter3IterL3MuonPixelSeeds_tmpL3Ref = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_dir = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_detId = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt_val = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_eta = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_phi = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_x = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_y = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_z = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_hasErr = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err0 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err1 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err2 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err3 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err4 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err5 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err6 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err7 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err8 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err9 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err10 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err11 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err12 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err13 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err14 = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_x = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_y = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dxdz = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dydz = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_px = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_py = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pz = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_qbp = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_charge = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Matched = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Ref = 0;
        hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tmpL3Ref = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_dir = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_detId = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt_val = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_eta = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_phi = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_x = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_y = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_z = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_hasErr = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err0 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err1 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err2 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err3 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err4 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err5 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err6 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err7 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err8 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err9 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err10 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err11 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err12 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err13 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_err14 = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_x = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_y = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_dxdz = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_dydz = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_px = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_py = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_pz = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_qbp = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tsos_charge = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_iterL3Matched = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_iterL3Ref = 0;
        hltIter2IterL3FromL1MuonPixelSeeds_tmpL3Ref = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_dir = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_detId = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt_val = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_eta = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_phi = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_x = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_y = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_z = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_hasErr = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err0 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err1 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err2 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err3 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err4 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err5 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err6 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err7 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err8 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err9 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err10 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err11 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err12 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err13 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_err14 = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_x = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_y = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_dxdz = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_dydz = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_px = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_py = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_pz = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_qbp = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tsos_charge = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_iterL3Matched = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_iterL3Ref = 0;
        hltIter3IterL3FromL1MuonPixelSeeds_tmpL3Ref = 0;
        hltIterL3OIMuonTrack_pt = 0;
        hltIterL3OIMuonTrack_ptError = 0;
        hltIterL3OIMuonTrack_eta = 0;
        hltIterL3OIMuonTrack_phi = 0;
        hltIterL3OIMuonTrack_charge = 0;
        hltIterL3OIMuonTrack_matchedL3 = 0;
        hltIterL3OIMuonTrack_matchedL3NoId = 0;
        hltIterL3OIMuonTrack_bestMatchTP_charge = 0;
        hltIterL3OIMuonTrack_bestMatchTP_pdgId = 0;
        hltIterL3OIMuonTrack_bestMatchTP_energy = 0;
        hltIterL3OIMuonTrack_bestMatchTP_pt = 0;
        hltIterL3OIMuonTrack_bestMatchTP_eta = 0;
        hltIterL3OIMuonTrack_bestMatchTP_phi = 0;
        hltIterL3OIMuonTrack_bestMatchTP_parentVx = 0;
        hltIterL3OIMuonTrack_bestMatchTP_parentVy = 0;
        hltIterL3OIMuonTrack_bestMatchTP_parentVz = 0;
        hltIterL3OIMuonTrack_bestMatchTP_status = 0;
        hltIterL3OIMuonTrack_bestMatchTP_numberOfHits = 0;
        hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerHits = 0;
        hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerLayers = 0;
        hltIterL3OIMuonTrack_bestMatchTP_sharedFraction = 0;
        hltIterL3OIMuonTrack_matchedTPsize = 0;
        hltIterL3OIMuonTrack_mva0 = 0;
        hltIterL3OIMuonTrack_mva1 = 0;
        hltIterL3OIMuonTrack_mva2 = 0;
        hltIterL3OIMuonTrack_mva3 = 0;
        hltIter0IterL3MuonTrack_pt = 0;
        hltIter0IterL3MuonTrack_ptError = 0;
        hltIter0IterL3MuonTrack_eta = 0;
        hltIter0IterL3MuonTrack_phi = 0;
        hltIter0IterL3MuonTrack_charge = 0;
        hltIter0IterL3MuonTrack_matchedL3 = 0;
        hltIter0IterL3MuonTrack_matchedL3NoId = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_charge = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_pdgId = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_energy = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_pt = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_eta = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_phi = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_parentVx = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_parentVy = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_parentVz = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_status = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_numberOfHits = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerHits = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter0IterL3MuonTrack_bestMatchTP_sharedFraction = 0;
        hltIter0IterL3MuonTrack_matchedTPsize = 0;
        hltIter0IterL3MuonTrack_mva0 = 0;
        hltIter0IterL3MuonTrack_mva1 = 0;
        hltIter0IterL3MuonTrack_mva2 = 0;
        hltIter0IterL3MuonTrack_mva3 = 0;
        hltIter2IterL3MuonTrack_pt = 0;
        hltIter2IterL3MuonTrack_ptError = 0;
        hltIter2IterL3MuonTrack_eta = 0;
        hltIter2IterL3MuonTrack_phi = 0;
        hltIter2IterL3MuonTrack_charge = 0;
        hltIter2IterL3MuonTrack_matchedL3 = 0;
        hltIter2IterL3MuonTrack_matchedL3NoId = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_charge = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_pdgId = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_energy = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_pt = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_eta = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_phi = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_parentVx = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_parentVy = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_parentVz = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_status = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_numberOfHits = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerHits = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter2IterL3MuonTrack_bestMatchTP_sharedFraction = 0;
        hltIter2IterL3MuonTrack_matchedTPsize = 0;
        hltIter2IterL3MuonTrack_mva0 = 0;
        hltIter2IterL3MuonTrack_mva1 = 0;
        hltIter2IterL3MuonTrack_mva2 = 0;
        hltIter2IterL3MuonTrack_mva3 = 0;
        hltIter3IterL3MuonTrack_pt = 0;
        hltIter3IterL3MuonTrack_ptError = 0;
        hltIter3IterL3MuonTrack_eta = 0;
        hltIter3IterL3MuonTrack_phi = 0;
        hltIter3IterL3MuonTrack_charge = 0;
        hltIter3IterL3MuonTrack_matchedL3 = 0;
        hltIter3IterL3MuonTrack_matchedL3NoId = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_charge = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_pdgId = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_energy = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_pt = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_eta = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_phi = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_parentVx = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_parentVy = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_parentVz = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_status = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_numberOfHits = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerHits = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter3IterL3MuonTrack_bestMatchTP_sharedFraction = 0;
        hltIter3IterL3MuonTrack_matchedTPsize = 0;
        hltIter3IterL3MuonTrack_mva0 = 0;
        hltIter3IterL3MuonTrack_mva1 = 0;
        hltIter3IterL3MuonTrack_mva2 = 0;
        hltIter3IterL3MuonTrack_mva3 = 0;
        hltIter0IterL3FromL1MuonTrack_pt = 0;
        hltIter0IterL3FromL1MuonTrack_ptError = 0;
        hltIter0IterL3FromL1MuonTrack_eta = 0;
        hltIter0IterL3FromL1MuonTrack_phi = 0;
        hltIter0IterL3FromL1MuonTrack_charge = 0;
        hltIter0IterL3FromL1MuonTrack_matchedL3 = 0;
        hltIter0IterL3FromL1MuonTrack_matchedL3NoId = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_charge = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_pdgId = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_energy = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_pt = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_eta = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_phi = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVx = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVy = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVz = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_status = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfHits = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter0IterL3FromL1MuonTrack_bestMatchTP_sharedFraction = 0;
        hltIter0IterL3FromL1MuonTrack_matchedTPsize = 0;
        hltIter0IterL3FromL1MuonTrack_mva0 = 0;
        hltIter0IterL3FromL1MuonTrack_mva1 = 0;
        hltIter0IterL3FromL1MuonTrack_mva2 = 0;
        hltIter0IterL3FromL1MuonTrack_mva3 = 0;
        hltIter2IterL3FromL1MuonTrack_pt = 0;
        hltIter2IterL3FromL1MuonTrack_ptError = 0;
        hltIter2IterL3FromL1MuonTrack_eta = 0;
        hltIter2IterL3FromL1MuonTrack_phi = 0;
        hltIter2IterL3FromL1MuonTrack_charge = 0;
        hltIter2IterL3FromL1MuonTrack_matchedL3 = 0;
        hltIter2IterL3FromL1MuonTrack_matchedL3NoId = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_charge = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_pdgId = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_energy = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_pt = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_eta = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_phi = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVx = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVy = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVz = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_status = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfHits = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter2IterL3FromL1MuonTrack_bestMatchTP_sharedFraction = 0;
        hltIter2IterL3FromL1MuonTrack_matchedTPsize = 0;
        hltIter2IterL3FromL1MuonTrack_mva0 = 0;
        hltIter2IterL3FromL1MuonTrack_mva1 = 0;
        hltIter2IterL3FromL1MuonTrack_mva2 = 0;
        hltIter2IterL3FromL1MuonTrack_mva3 = 0;
        hltIter3IterL3FromL1MuonTrack_pt = 0;
        hltIter3IterL3FromL1MuonTrack_ptError = 0;
        hltIter3IterL3FromL1MuonTrack_eta = 0;
        hltIter3IterL3FromL1MuonTrack_phi = 0;
        hltIter3IterL3FromL1MuonTrack_charge = 0;
        hltIter3IterL3FromL1MuonTrack_matchedL3 = 0;
        hltIter3IterL3FromL1MuonTrack_matchedL3NoId = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_charge = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_pdgId = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_energy = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_pt = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_eta = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_phi = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVx = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVy = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVz = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_status = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfHits = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter3IterL3FromL1MuonTrack_bestMatchTP_sharedFraction = 0;
        hltIter3IterL3FromL1MuonTrack_matchedTPsize = 0;
        hltIter3IterL3FromL1MuonTrack_mva0 = 0;
        hltIter3IterL3FromL1MuonTrack_mva1 = 0;
        hltIter3IterL3FromL1MuonTrack_mva2 = 0;
        hltIter3IterL3FromL1MuonTrack_mva3 = 0;

        L3MuonsNoId_pt = 0;
        L3MuonsNoId_inner_pt = 0;
        L3MuonsNoId_inner_ptError = 0;
        L3MuonsNoId_eta = 0;
        L3MuonsNoId_phi = 0;
        L3MuonsNoId_charge = 0;
        L3MuonsNoId_isGlobalMuon = 0;
        L3MuonsNoId_isStandAloneMuon = 0;
        L3MuonsNoId_isTrackerMuon = 0;
        L3MuonsNoId_isLooseTriggerMuon = 0;
        L3MuonsNoId_isME0Muon = 0;
        L3MuonsNoId_isGEMMuon = 0;
        L3MuonsNoId_isRPCMuon = 0;
        L3MuonsNoId_isGoodMuon_TMOneStationTight = 0;
        L3MuonsNoId_numberOfMatchedStations = 0;
        L3MuonsNoId_numberOfMatchedRPCLayers = 0;
        L3MuonsNoId_expectedNnumberOfMatchedStations = 0;
        L3MuonsNoId_inner_normalizedChi2 = 0;
        L3MuonsNoId_inner_numberOfValidTrackerHits = 0;
        L3MuonsNoId_inner_trackerLayersWithMeasurement = 0;
        L3MuonsNoId_inner_numberOfValidPixelHits = 0;
        L3MuonsNoId_inner_dz_l1vtx = 0;
        L3Muons_pt = 0;
        L3Muons_inner_pt = 0;
        L3Muons_inner_ptError = 0;
        L3Muons_eta = 0;
        L3Muons_phi = 0;
        L3Muons_charge = 0;
        L3Muons_isGlobalMuon = 0;
        L3Muons_isStandAloneMuon = 0;
        L3Muons_isTrackerMuon = 0;
        L3Muons_isLooseTriggerMuon = 0;
        L3Muons_isME0Muon = 0;
        L3Muons_isGEMMuon = 0;
        L3Muons_isRPCMuon = 0;
        L3Muons_isGoodMuon_TMOneStationTight = 0;
        L3Muons_numberOfMatchedStations = 0;
        L3Muons_numberOfMatchedRPCLayers = 0;
        L3Muons_expectedNnumberOfMatchedStations = 0;
        L3Muons_inner_normalizedChi2 = 0;
        L3Muons_inner_numberOfValidTrackerHits = 0;
        L3Muons_inner_trackerLayersWithMeasurement = 0;
        L3Muons_inner_numberOfValidPixelHits = 0;
        L3Muons_inner_dz_l1vtx = 0;

        TP_charge = 0;
        TP_pdgId = 0;
        TP_energy = 0;
        TP_pt = 0;
        TP_eta = 0;
        TP_phi = 0;
        TP_parentVx = 0;
        TP_parentVy = 0;
        TP_parentVz = 0;
        TP_status = 0;
        TP_numberOfHits = 0;
        TP_numberOfTrackerHits = 0;
        TP_numberOfTrackerLayers = 0;
        TP_gen_charge = 0;
        TP_gen_pdgId = 0;
        TP_gen_pt = 0;
        TP_gen_eta = 0;
        TP_gen_phi = 0;
        TP_bestMatchTrk_pt = 0;
        TP_bestMatchTrk_eta = 0;
        TP_bestMatchTrk_phi = 0;
        TP_bestMatchTrk_charge = 0;
        TP_bestMatchTrk_quality = 0;
        TP_bestMatchTrk_NValidHits = 0;
        hltIterL3MuonTrimmedPixelVertices_isValid = 0;
        hltIterL3MuonTrimmedPixelVertices_chi2 = 0;
        hltIterL3MuonTrimmedPixelVertices_ndof = 0;
        hltIterL3MuonTrimmedPixelVertices_nTracks = 0;
        hltIterL3MuonTrimmedPixelVertices_x = 0;
        hltIterL3MuonTrimmedPixelVertices_xerr = 0;
        hltIterL3MuonTrimmedPixelVertices_y = 0;
        hltIterL3MuonTrimmedPixelVertices_yerr = 0;
        hltIterL3MuonTrimmedPixelVertices_z = 0;
        hltIterL3MuonTrimmedPixelVertices_zerr = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_isValid = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_chi2 = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_ndof = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_nTracks = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_x = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_xerr = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_y = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_yerr = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_z = 0;
        hltIterL3FromL1MuonTrimmedPixelVertices_zerr = 0;

        hltPhase2L3OI_pt = 0;
        hltPhase2L3OI_ptError = 0;
        hltPhase2L3OI_eta = 0;
        hltPhase2L3OI_phi = 0;
        hltPhase2L3OI_charge = 0;
        hltPhase2L3OI_matchedL3 = 0;
        hltPhase2L3OI_matchedL3NoId = 0;
        hltPhase2L3OI_bestMatchTP_charge = 0;
        hltPhase2L3OI_bestMatchTP_pdgId = 0;
        hltPhase2L3OI_bestMatchTP_energy = 0;
        hltPhase2L3OI_bestMatchTP_pt = 0;
        hltPhase2L3OI_bestMatchTP_eta = 0;
        hltPhase2L3OI_bestMatchTP_phi = 0;
        hltPhase2L3OI_bestMatchTP_parentVx = 0;
        hltPhase2L3OI_bestMatchTP_parentVy = 0;
        hltPhase2L3OI_bestMatchTP_parentVz = 0;
        hltPhase2L3OI_bestMatchTP_status = 0;
        hltPhase2L3OI_bestMatchTP_numberOfHits = 0;
        hltPhase2L3OI_bestMatchTP_numberOfTrackerHits = 0;
        hltPhase2L3OI_bestMatchTP_numberOfTrackerLayers = 0;
        hltPhase2L3OI_bestMatchTP_sharedFraction = 0;
        hltPhase2L3OI_matchedTPsize = 0;
        hltPhase2L3OI_mva0 = 0;
        hltPhase2L3OI_mva1 = 0;
        hltPhase2L3OI_mva2 = 0;
        hltPhase2L3OI_mva3 = 0;
        tpTo_hltPhase2L3OI_charge = 0;
        tpTo_hltPhase2L3OI_pdgId = 0;
        tpTo_hltPhase2L3OI_energy = 0;
        tpTo_hltPhase2L3OI_pt = 0;
        tpTo_hltPhase2L3OI_eta = 0;
        tpTo_hltPhase2L3OI_phi = 0;
        tpTo_hltPhase2L3OI_parentVx = 0;
        tpTo_hltPhase2L3OI_parentVy = 0;
        tpTo_hltPhase2L3OI_parentVz = 0;
        tpTo_hltPhase2L3OI_status = 0;
        tpTo_hltPhase2L3OI_numberOfHits = 0;
        tpTo_hltPhase2L3OI_numberOfTrackerHits = 0;
        tpTo_hltPhase2L3OI_numberOfTrackerLayers = 0;
        tpTo_hltPhase2L3OI_gen_charge = 0;
        tpTo_hltPhase2L3OI_gen_pdgId = 0;
        tpTo_hltPhase2L3OI_gen_pt = 0;
        tpTo_hltPhase2L3OI_gen_eta = 0;
        tpTo_hltPhase2L3OI_gen_phi = 0;
        tpTo_hltPhase2L3OI_bestMatchTrk_pt = 0;
        tpTo_hltPhase2L3OI_bestMatchTrk_eta = 0;
        tpTo_hltPhase2L3OI_bestMatchTrk_phi = 0;
        tpTo_hltPhase2L3OI_bestMatchTrk_charge = 0;
        tpTo_hltPhase2L3OI_bestMatchTrk_quality = 0;
        tpTo_hltPhase2L3OI_bestMatchTrk_NValidHits = 0;

        hltIter0Phase2L3FromL1TkMuon_pt = 0;
        hltIter0Phase2L3FromL1TkMuon_ptError = 0;
        hltIter0Phase2L3FromL1TkMuon_eta = 0;
        hltIter0Phase2L3FromL1TkMuon_phi = 0;
        hltIter0Phase2L3FromL1TkMuon_charge = 0;
        hltIter0Phase2L3FromL1TkMuon_matchedL3 = 0;
        hltIter0Phase2L3FromL1TkMuon_matchedL3NoId = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_charge = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pdgId = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_energy = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pt = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_eta = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_phi = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVx = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVy = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVz = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_status = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter0Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction = 0;
        hltIter0Phase2L3FromL1TkMuon_matchedTPsize = 0;
        hltIter0Phase2L3FromL1TkMuon_mva0 = 0;
        hltIter0Phase2L3FromL1TkMuon_mva1 = 0;
        hltIter0Phase2L3FromL1TkMuon_mva2 = 0;
        hltIter0Phase2L3FromL1TkMuon_mva3 = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_charge = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_pdgId = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_energy = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_pt = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_eta = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_phi = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_parentVx = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_parentVy = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_parentVz = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_status = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfHits = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerHits = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerLayers = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_gen_charge = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pdgId = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pt = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_gen_eta = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_gen_phi = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_pt = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_eta = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_phi = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_charge = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_quality = 0;
        tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits = 0;
        hltIter2Phase2L3FromL1TkMuon_pt = 0;
        hltIter2Phase2L3FromL1TkMuon_ptError = 0;
        hltIter2Phase2L3FromL1TkMuon_eta = 0;
        hltIter2Phase2L3FromL1TkMuon_phi = 0;
        hltIter2Phase2L3FromL1TkMuon_charge = 0;
        hltIter2Phase2L3FromL1TkMuon_matchedL3 = 0;
        hltIter2Phase2L3FromL1TkMuon_matchedL3NoId = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_charge = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pdgId = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_energy = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pt = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_eta = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_phi = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVx = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVy = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVz = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_status = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter2Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction = 0;
        hltIter2Phase2L3FromL1TkMuon_matchedTPsize = 0;
        hltIter2Phase2L3FromL1TkMuon_mva0 = 0;
        hltIter2Phase2L3FromL1TkMuon_mva1 = 0;
        hltIter2Phase2L3FromL1TkMuon_mva2 = 0;
        hltIter2Phase2L3FromL1TkMuon_mva3 = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_charge = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_pdgId = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_energy = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_pt = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_eta = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_phi = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_parentVx = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_parentVy = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_parentVz = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_status = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfHits = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerHits = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerLayers = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_gen_charge = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pdgId = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pt = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_gen_eta = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_gen_phi = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_pt = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_eta = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_phi = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_charge = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_quality = 0;
        tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits = 0;
        hltPhase2L3IOFromL1_pt = 0;
        hltPhase2L3IOFromL1_ptError = 0;
        hltPhase2L3IOFromL1_eta = 0;
        hltPhase2L3IOFromL1_phi = 0;
        hltPhase2L3IOFromL1_charge = 0;
        hltPhase2L3IOFromL1_matchedL3 = 0;
        hltPhase2L3IOFromL1_matchedL3NoId = 0;
        hltPhase2L3IOFromL1_bestMatchTP_charge = 0;
        hltPhase2L3IOFromL1_bestMatchTP_pdgId = 0;
        hltPhase2L3IOFromL1_bestMatchTP_energy = 0;
        hltPhase2L3IOFromL1_bestMatchTP_pt = 0;
        hltPhase2L3IOFromL1_bestMatchTP_eta = 0;
        hltPhase2L3IOFromL1_bestMatchTP_phi = 0;
        hltPhase2L3IOFromL1_bestMatchTP_parentVx = 0;
        hltPhase2L3IOFromL1_bestMatchTP_parentVy = 0;
        hltPhase2L3IOFromL1_bestMatchTP_parentVz = 0;
        hltPhase2L3IOFromL1_bestMatchTP_status = 0;
        hltPhase2L3IOFromL1_bestMatchTP_numberOfHits = 0;
        hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerHits = 0;
        hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerLayers = 0;
        hltPhase2L3IOFromL1_bestMatchTP_sharedFraction = 0;
        hltPhase2L3IOFromL1_matchedTPsize = 0;
        hltPhase2L3IOFromL1_mva0 = 0;
        hltPhase2L3IOFromL1_mva1 = 0;
        hltPhase2L3IOFromL1_mva2 = 0;
        hltPhase2L3IOFromL1_mva3 = 0;
        tpTo_hltPhase2L3IOFromL1_charge = 0;
        tpTo_hltPhase2L3IOFromL1_pdgId = 0;
        tpTo_hltPhase2L3IOFromL1_energy = 0;
        tpTo_hltPhase2L3IOFromL1_pt = 0;
        tpTo_hltPhase2L3IOFromL1_eta = 0;
        tpTo_hltPhase2L3IOFromL1_phi = 0;
        tpTo_hltPhase2L3IOFromL1_parentVx = 0;
        tpTo_hltPhase2L3IOFromL1_parentVy = 0;
        tpTo_hltPhase2L3IOFromL1_parentVz = 0;
        tpTo_hltPhase2L3IOFromL1_status = 0;
        tpTo_hltPhase2L3IOFromL1_numberOfHits = 0;
        tpTo_hltPhase2L3IOFromL1_numberOfTrackerHits = 0;
        tpTo_hltPhase2L3IOFromL1_numberOfTrackerLayers = 0;
        tpTo_hltPhase2L3IOFromL1_gen_charge = 0;
        tpTo_hltPhase2L3IOFromL1_gen_pdgId = 0;
        tpTo_hltPhase2L3IOFromL1_gen_pt = 0;
        tpTo_hltPhase2L3IOFromL1_gen_eta = 0;
        tpTo_hltPhase2L3IOFromL1_gen_phi = 0;
        tpTo_hltPhase2L3IOFromL1_bestMatchTrk_pt = 0;
        tpTo_hltPhase2L3IOFromL1_bestMatchTrk_eta = 0;
        tpTo_hltPhase2L3IOFromL1_bestMatchTrk_phi = 0;
        tpTo_hltPhase2L3IOFromL1_bestMatchTrk_charge = 0;
        tpTo_hltPhase2L3IOFromL1_bestMatchTrk_quality = 0;
        tpTo_hltPhase2L3IOFromL1_bestMatchTrk_NValidHits = 0;
        hltPhase2L3MuonsNoID_pt = 0;
        hltPhase2L3MuonsNoID_ptError = 0;
        hltPhase2L3MuonsNoID_eta = 0;
        hltPhase2L3MuonsNoID_phi = 0;
        hltPhase2L3MuonsNoID_charge = 0;
        hltPhase2L3MuonsNoID_matchedL3 = 0;
        hltPhase2L3MuonsNoID_matchedL3NoId = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_charge = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_pdgId = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_energy = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_pt = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_eta = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_phi = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_parentVx = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_parentVy = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_parentVz = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_status = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_numberOfHits = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerHits = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerLayers = 0;
        hltPhase2L3MuonsNoID_bestMatchTP_sharedFraction = 0;
        hltPhase2L3MuonsNoID_matchedTPsize = 0;
        hltPhase2L3MuonsNoID_mva0 = 0;
        hltPhase2L3MuonsNoID_mva1 = 0;
        hltPhase2L3MuonsNoID_mva2 = 0;
        hltPhase2L3MuonsNoID_mva3 = 0;
        tpTo_hltPhase2L3MuonsNoID_charge = 0;
        tpTo_hltPhase2L3MuonsNoID_pdgId = 0;
        tpTo_hltPhase2L3MuonsNoID_energy = 0;
        tpTo_hltPhase2L3MuonsNoID_pt = 0;
        tpTo_hltPhase2L3MuonsNoID_eta = 0;
        tpTo_hltPhase2L3MuonsNoID_phi = 0;
        tpTo_hltPhase2L3MuonsNoID_parentVx = 0;
        tpTo_hltPhase2L3MuonsNoID_parentVy = 0;
        tpTo_hltPhase2L3MuonsNoID_parentVz = 0;
        tpTo_hltPhase2L3MuonsNoID_status = 0;
        tpTo_hltPhase2L3MuonsNoID_numberOfHits = 0;
        tpTo_hltPhase2L3MuonsNoID_numberOfTrackerHits = 0;
        tpTo_hltPhase2L3MuonsNoID_numberOfTrackerLayers = 0;
        tpTo_hltPhase2L3MuonsNoID_gen_charge = 0;
        tpTo_hltPhase2L3MuonsNoID_gen_pdgId = 0;
        tpTo_hltPhase2L3MuonsNoID_gen_pt = 0;
        tpTo_hltPhase2L3MuonsNoID_gen_eta = 0;
        tpTo_hltPhase2L3MuonsNoID_gen_phi = 0;
        tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_pt = 0;
        tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_eta = 0;
        tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_phi = 0;
        tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_charge = 0;
        tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_quality = 0;
        tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_NValidHits = 0;
        hltPhase2L3Muons_pt = 0;
        hltPhase2L3Muons_ptError = 0;
        hltPhase2L3Muons_eta = 0;
        hltPhase2L3Muons_phi = 0;
        hltPhase2L3Muons_charge = 0;
        hltPhase2L3Muons_matchedL3 = 0;
        hltPhase2L3Muons_matchedL3NoId = 0;
        hltPhase2L3Muons_bestMatchTP_charge = 0;
        hltPhase2L3Muons_bestMatchTP_pdgId = 0;
        hltPhase2L3Muons_bestMatchTP_energy = 0;
        hltPhase2L3Muons_bestMatchTP_pt = 0;
        hltPhase2L3Muons_bestMatchTP_eta = 0;
        hltPhase2L3Muons_bestMatchTP_phi = 0;
        hltPhase2L3Muons_bestMatchTP_parentVx = 0;
        hltPhase2L3Muons_bestMatchTP_parentVy = 0;
        hltPhase2L3Muons_bestMatchTP_parentVz = 0;
        hltPhase2L3Muons_bestMatchTP_status = 0;
        hltPhase2L3Muons_bestMatchTP_numberOfHits = 0;
        hltPhase2L3Muons_bestMatchTP_numberOfTrackerHits = 0;
        hltPhase2L3Muons_bestMatchTP_numberOfTrackerLayers = 0;
        hltPhase2L3Muons_bestMatchTP_sharedFraction = 0;
        hltPhase2L3Muons_matchedTPsize = 0;
        hltPhase2L3Muons_mva0 = 0;
        hltPhase2L3Muons_mva1 = 0;
        hltPhase2L3Muons_mva2 = 0;
        hltPhase2L3Muons_mva3 = 0;

        // -- Isolation
        hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0 = 0;
        hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p000 = 0;
        hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p000 = 0;
        hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p030 = 0;
        hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p030 = 0;
        hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p050 = 0;
        hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p050 = 0;
        hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p00minEEM0p00minEHad0p00 = 0;
        hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00 = 0;
        hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p04minEEM0p00minEHad0p00 = 0;
        hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p00minEEM0p00minEHad0p00 = 0;
        hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p02minEEM0p00minEHad0p00 = 0;
        hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p04minEEM0p00minEHad0p00 = 0;
        hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p00minEEM0p00minEHad0p00 = 0;
        hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p02minEEM0p00minEHad0p00 = 0;
        hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p04minEEM0p00minEHad0p00 = 0;

        tpTo_hltPhase2L3Muons_charge = 0;
        tpTo_hltPhase2L3Muons_pdgId = 0;
        tpTo_hltPhase2L3Muons_energy = 0;
        tpTo_hltPhase2L3Muons_pt = 0;
        tpTo_hltPhase2L3Muons_eta = 0;
        tpTo_hltPhase2L3Muons_phi = 0;
        tpTo_hltPhase2L3Muons_parentVx = 0;
        tpTo_hltPhase2L3Muons_parentVy = 0;
        tpTo_hltPhase2L3Muons_parentVz = 0;
        tpTo_hltPhase2L3Muons_status = 0;
        tpTo_hltPhase2L3Muons_numberOfHits = 0;
        tpTo_hltPhase2L3Muons_numberOfTrackerHits = 0;
        tpTo_hltPhase2L3Muons_numberOfTrackerLayers = 0;
        tpTo_hltPhase2L3Muons_gen_charge = 0;
        tpTo_hltPhase2L3Muons_gen_pdgId = 0;
        tpTo_hltPhase2L3Muons_gen_pt = 0;
        tpTo_hltPhase2L3Muons_gen_eta = 0;
        tpTo_hltPhase2L3Muons_gen_phi = 0;
        tpTo_hltPhase2L3Muons_bestMatchTrk_pt = 0;
        tpTo_hltPhase2L3Muons_bestMatchTrk_eta = 0;
        tpTo_hltPhase2L3Muons_bestMatchTrk_phi = 0;
        tpTo_hltPhase2L3Muons_bestMatchTrk_charge = 0;
        tpTo_hltPhase2L3Muons_bestMatchTrk_quality = 0;
        tpTo_hltPhase2L3Muons_bestMatchTrk_NValidHits = 0;

        hltIterL3OI_pt = 0;
        hltIterL3OI_ptError = 0;
        hltIterL3OI_eta = 0;
        hltIterL3OI_phi = 0;
        hltIterL3OI_charge = 0;
        hltIterL3OI_matchedL3 = 0;
        hltIterL3OI_matchedL3NoId = 0;
        hltIterL3OI_bestMatchTP_charge = 0;
        hltIterL3OI_bestMatchTP_pdgId = 0;
        hltIterL3OI_bestMatchTP_energy = 0;
        hltIterL3OI_bestMatchTP_pt = 0;
        hltIterL3OI_bestMatchTP_eta = 0;
        hltIterL3OI_bestMatchTP_phi = 0;
        hltIterL3OI_bestMatchTP_parentVx = 0;
        hltIterL3OI_bestMatchTP_parentVy = 0;
        hltIterL3OI_bestMatchTP_parentVz = 0;
        hltIterL3OI_bestMatchTP_status = 0;
        hltIterL3OI_bestMatchTP_numberOfHits = 0;
        hltIterL3OI_bestMatchTP_numberOfTrackerHits = 0;
        hltIterL3OI_bestMatchTP_numberOfTrackerLayers = 0;
        hltIterL3OI_bestMatchTP_sharedFraction = 0;
        hltIterL3OI_matchedTPsize = 0;
        hltIterL3OI_mva0 = 0;
        hltIterL3OI_mva1 = 0;
        hltIterL3OI_mva2 = 0;
        hltIterL3OI_mva3 = 0;
        tpTo_hltIterL3OI_charge = 0;
        tpTo_hltIterL3OI_pdgId = 0;
        tpTo_hltIterL3OI_energy = 0;
        tpTo_hltIterL3OI_pt = 0;
        tpTo_hltIterL3OI_eta = 0;
        tpTo_hltIterL3OI_phi = 0;
        tpTo_hltIterL3OI_parentVx = 0;
        tpTo_hltIterL3OI_parentVy = 0;
        tpTo_hltIterL3OI_parentVz = 0;
        tpTo_hltIterL3OI_status = 0;
        tpTo_hltIterL3OI_numberOfHits = 0;
        tpTo_hltIterL3OI_numberOfTrackerHits = 0;
        tpTo_hltIterL3OI_numberOfTrackerLayers = 0;
        tpTo_hltIterL3OI_gen_charge = 0;
        tpTo_hltIterL3OI_gen_pdgId = 0;
        tpTo_hltIterL3OI_gen_pt = 0;
        tpTo_hltIterL3OI_gen_eta = 0;
        tpTo_hltIterL3OI_gen_phi = 0;
        tpTo_hltIterL3OI_bestMatchTrk_pt = 0;
        tpTo_hltIterL3OI_bestMatchTrk_eta = 0;
        tpTo_hltIterL3OI_bestMatchTrk_phi = 0;
        tpTo_hltIterL3OI_bestMatchTrk_charge = 0;
        tpTo_hltIterL3OI_bestMatchTrk_quality = 0;
        tpTo_hltIterL3OI_bestMatchTrk_NValidHits = 0;
        hltIter0IterL3_pt = 0;
        hltIter0IterL3_ptError = 0;
        hltIter0IterL3_eta = 0;
        hltIter0IterL3_phi = 0;
        hltIter0IterL3_charge = 0;
        hltIter0IterL3_matchedL3 = 0;
        hltIter0IterL3_matchedL3NoId = 0;
        hltIter0IterL3_bestMatchTP_charge = 0;
        hltIter0IterL3_bestMatchTP_pdgId = 0;
        hltIter0IterL3_bestMatchTP_energy = 0;
        hltIter0IterL3_bestMatchTP_pt = 0;
        hltIter0IterL3_bestMatchTP_eta = 0;
        hltIter0IterL3_bestMatchTP_phi = 0;
        hltIter0IterL3_bestMatchTP_parentVx = 0;
        hltIter0IterL3_bestMatchTP_parentVy = 0;
        hltIter0IterL3_bestMatchTP_parentVz = 0;
        hltIter0IterL3_bestMatchTP_status = 0;
        hltIter0IterL3_bestMatchTP_numberOfHits = 0;
        hltIter0IterL3_bestMatchTP_numberOfTrackerHits = 0;
        hltIter0IterL3_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter0IterL3_bestMatchTP_sharedFraction = 0;
        hltIter0IterL3_matchedTPsize = 0;
        hltIter0IterL3_mva0 = 0;
        hltIter0IterL3_mva1 = 0;
        hltIter0IterL3_mva2 = 0;
        hltIter0IterL3_mva3 = 0;
        tpTo_hltIter0IterL3_charge = 0;
        tpTo_hltIter0IterL3_pdgId = 0;
        tpTo_hltIter0IterL3_energy = 0;
        tpTo_hltIter0IterL3_pt = 0;
        tpTo_hltIter0IterL3_eta = 0;
        tpTo_hltIter0IterL3_phi = 0;
        tpTo_hltIter0IterL3_parentVx = 0;
        tpTo_hltIter0IterL3_parentVy = 0;
        tpTo_hltIter0IterL3_parentVz = 0;
        tpTo_hltIter0IterL3_status = 0;
        tpTo_hltIter0IterL3_numberOfHits = 0;
        tpTo_hltIter0IterL3_numberOfTrackerHits = 0;
        tpTo_hltIter0IterL3_numberOfTrackerLayers = 0;
        tpTo_hltIter0IterL3_gen_charge = 0;
        tpTo_hltIter0IterL3_gen_pdgId = 0;
        tpTo_hltIter0IterL3_gen_pt = 0;
        tpTo_hltIter0IterL3_gen_eta = 0;
        tpTo_hltIter0IterL3_gen_phi = 0;
        tpTo_hltIter0IterL3_bestMatchTrk_pt = 0;
        tpTo_hltIter0IterL3_bestMatchTrk_eta = 0;
        tpTo_hltIter0IterL3_bestMatchTrk_phi = 0;
        tpTo_hltIter0IterL3_bestMatchTrk_charge = 0;
        tpTo_hltIter0IterL3_bestMatchTrk_quality = 0;
        tpTo_hltIter0IterL3_bestMatchTrk_NValidHits = 0;
        hltIter2IterL3_pt = 0;
        hltIter2IterL3_ptError = 0;
        hltIter2IterL3_eta = 0;
        hltIter2IterL3_phi = 0;
        hltIter2IterL3_charge = 0;
        hltIter2IterL3_matchedL3 = 0;
        hltIter2IterL3_matchedL3NoId = 0;
        hltIter2IterL3_bestMatchTP_charge = 0;
        hltIter2IterL3_bestMatchTP_pdgId = 0;
        hltIter2IterL3_bestMatchTP_energy = 0;
        hltIter2IterL3_bestMatchTP_pt = 0;
        hltIter2IterL3_bestMatchTP_eta = 0;
        hltIter2IterL3_bestMatchTP_phi = 0;
        hltIter2IterL3_bestMatchTP_parentVx = 0;
        hltIter2IterL3_bestMatchTP_parentVy = 0;
        hltIter2IterL3_bestMatchTP_parentVz = 0;
        hltIter2IterL3_bestMatchTP_status = 0;
        hltIter2IterL3_bestMatchTP_numberOfHits = 0;
        hltIter2IterL3_bestMatchTP_numberOfTrackerHits = 0;
        hltIter2IterL3_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter2IterL3_bestMatchTP_sharedFraction = 0;
        hltIter2IterL3_matchedTPsize = 0;
        hltIter2IterL3_mva0 = 0;
        hltIter2IterL3_mva1 = 0;
        hltIter2IterL3_mva2 = 0;
        hltIter2IterL3_mva3 = 0;
        tpTo_hltIter2IterL3_charge = 0;
        tpTo_hltIter2IterL3_pdgId = 0;
        tpTo_hltIter2IterL3_energy = 0;
        tpTo_hltIter2IterL3_pt = 0;
        tpTo_hltIter2IterL3_eta = 0;
        tpTo_hltIter2IterL3_phi = 0;
        tpTo_hltIter2IterL3_parentVx = 0;
        tpTo_hltIter2IterL3_parentVy = 0;
        tpTo_hltIter2IterL3_parentVz = 0;
        tpTo_hltIter2IterL3_status = 0;
        tpTo_hltIter2IterL3_numberOfHits = 0;
        tpTo_hltIter2IterL3_numberOfTrackerHits = 0;
        tpTo_hltIter2IterL3_numberOfTrackerLayers = 0;
        tpTo_hltIter2IterL3_gen_charge = 0;
        tpTo_hltIter2IterL3_gen_pdgId = 0;
        tpTo_hltIter2IterL3_gen_pt = 0;
        tpTo_hltIter2IterL3_gen_eta = 0;
        tpTo_hltIter2IterL3_gen_phi = 0;
        tpTo_hltIter2IterL3_bestMatchTrk_pt = 0;
        tpTo_hltIter2IterL3_bestMatchTrk_eta = 0;
        tpTo_hltIter2IterL3_bestMatchTrk_phi = 0;
        tpTo_hltIter2IterL3_bestMatchTrk_charge = 0;
        tpTo_hltIter2IterL3_bestMatchTrk_quality = 0;
        tpTo_hltIter2IterL3_bestMatchTrk_NValidHits = 0;
        hltIter0IterL3FromL1Muon_pt = 0;
        hltIter0IterL3FromL1Muon_ptError = 0;
        hltIter0IterL3FromL1Muon_eta = 0;
        hltIter0IterL3FromL1Muon_phi = 0;
        hltIter0IterL3FromL1Muon_charge = 0;
        hltIter0IterL3FromL1Muon_matchedL3 = 0;
        hltIter0IterL3FromL1Muon_matchedL3NoId = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_charge = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_pdgId = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_energy = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_pt = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_eta = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_phi = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_parentVx = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_parentVy = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_parentVz = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_status = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_numberOfHits = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter0IterL3FromL1Muon_bestMatchTP_sharedFraction = 0;
        hltIter0IterL3FromL1Muon_matchedTPsize = 0;
        hltIter0IterL3FromL1Muon_mva0 = 0;
        hltIter0IterL3FromL1Muon_mva1 = 0;
        hltIter0IterL3FromL1Muon_mva2 = 0;
        hltIter0IterL3FromL1Muon_mva3 = 0;
        tpTo_hltIter0IterL3FromL1Muon_charge = 0;
        tpTo_hltIter0IterL3FromL1Muon_pdgId = 0;
        tpTo_hltIter0IterL3FromL1Muon_energy = 0;
        tpTo_hltIter0IterL3FromL1Muon_pt = 0;
        tpTo_hltIter0IterL3FromL1Muon_eta = 0;
        tpTo_hltIter0IterL3FromL1Muon_phi = 0;
        tpTo_hltIter0IterL3FromL1Muon_parentVx = 0;
        tpTo_hltIter0IterL3FromL1Muon_parentVy = 0;
        tpTo_hltIter0IterL3FromL1Muon_parentVz = 0;
        tpTo_hltIter0IterL3FromL1Muon_status = 0;
        tpTo_hltIter0IterL3FromL1Muon_numberOfHits = 0;
        tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerHits = 0;
        tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerLayers = 0;
        tpTo_hltIter0IterL3FromL1Muon_gen_charge = 0;
        tpTo_hltIter0IterL3FromL1Muon_gen_pdgId = 0;
        tpTo_hltIter0IterL3FromL1Muon_gen_pt = 0;
        tpTo_hltIter0IterL3FromL1Muon_gen_eta = 0;
        tpTo_hltIter0IterL3FromL1Muon_gen_phi = 0;
        tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_pt = 0;
        tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_eta = 0;
        tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_phi = 0;
        tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_charge = 0;
        tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_quality = 0;
        tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_NValidHits = 0;
        hltIter2IterL3FromL1Muon_pt = 0;
        hltIter2IterL3FromL1Muon_ptError = 0;
        hltIter2IterL3FromL1Muon_eta = 0;
        hltIter2IterL3FromL1Muon_phi = 0;
        hltIter2IterL3FromL1Muon_charge = 0;
        hltIter2IterL3FromL1Muon_matchedL3 = 0;
        hltIter2IterL3FromL1Muon_matchedL3NoId = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_charge = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_pdgId = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_energy = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_pt = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_eta = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_phi = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_parentVx = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_parentVy = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_parentVz = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_status = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_numberOfHits = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers = 0;
        hltIter2IterL3FromL1Muon_bestMatchTP_sharedFraction = 0;
        hltIter2IterL3FromL1Muon_matchedTPsize = 0;
        hltIter2IterL3FromL1Muon_mva0 = 0;
        hltIter2IterL3FromL1Muon_mva1 = 0;
        hltIter2IterL3FromL1Muon_mva2 = 0;
        hltIter2IterL3FromL1Muon_mva3 = 0;
        tpTo_hltIter2IterL3FromL1Muon_charge = 0;
        tpTo_hltIter2IterL3FromL1Muon_pdgId = 0;
        tpTo_hltIter2IterL3FromL1Muon_energy = 0;
        tpTo_hltIter2IterL3FromL1Muon_pt = 0;
        tpTo_hltIter2IterL3FromL1Muon_eta = 0;
        tpTo_hltIter2IterL3FromL1Muon_phi = 0;
        tpTo_hltIter2IterL3FromL1Muon_parentVx = 0;
        tpTo_hltIter2IterL3FromL1Muon_parentVy = 0;
        tpTo_hltIter2IterL3FromL1Muon_parentVz = 0;
        tpTo_hltIter2IterL3FromL1Muon_status = 0;
        tpTo_hltIter2IterL3FromL1Muon_numberOfHits = 0;
        tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerHits = 0;
        tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerLayers = 0;
        tpTo_hltIter2IterL3FromL1Muon_gen_charge = 0;
        tpTo_hltIter2IterL3FromL1Muon_gen_pdgId = 0;
        tpTo_hltIter2IterL3FromL1Muon_gen_pt = 0;
        tpTo_hltIter2IterL3FromL1Muon_gen_eta = 0;
        tpTo_hltIter2IterL3FromL1Muon_gen_phi = 0;
        tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_pt = 0;
        tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_eta = 0;
        tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_phi = 0;
        tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_charge = 0;
        tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_quality = 0;
        tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_NValidHits = 0;
        hltIterL3IOFromL1_pt = 0;
        hltIterL3IOFromL1_ptError = 0;
        hltIterL3IOFromL1_eta = 0;
        hltIterL3IOFromL1_phi = 0;
        hltIterL3IOFromL1_charge = 0;
        hltIterL3IOFromL1_matchedL3 = 0;
        hltIterL3IOFromL1_matchedL3NoId = 0;
        hltIterL3IOFromL1_bestMatchTP_charge = 0;
        hltIterL3IOFromL1_bestMatchTP_pdgId = 0;
        hltIterL3IOFromL1_bestMatchTP_energy = 0;
        hltIterL3IOFromL1_bestMatchTP_pt = 0;
        hltIterL3IOFromL1_bestMatchTP_eta = 0;
        hltIterL3IOFromL1_bestMatchTP_phi = 0;
        hltIterL3IOFromL1_bestMatchTP_parentVx = 0;
        hltIterL3IOFromL1_bestMatchTP_parentVy = 0;
        hltIterL3IOFromL1_bestMatchTP_parentVz = 0;
        hltIterL3IOFromL1_bestMatchTP_status = 0;
        hltIterL3IOFromL1_bestMatchTP_numberOfHits = 0;
        hltIterL3IOFromL1_bestMatchTP_numberOfTrackerHits = 0;
        hltIterL3IOFromL1_bestMatchTP_numberOfTrackerLayers = 0;
        hltIterL3IOFromL1_bestMatchTP_sharedFraction = 0;
        hltIterL3IOFromL1_matchedTPsize = 0;
        hltIterL3IOFromL1_mva0 = 0;
        hltIterL3IOFromL1_mva1 = 0;
        hltIterL3IOFromL1_mva2 = 0;
        hltIterL3IOFromL1_mva3 = 0;
        tpTo_hltIterL3IOFromL1_charge = 0;
        tpTo_hltIterL3IOFromL1_pdgId = 0;
        tpTo_hltIterL3IOFromL1_energy = 0;
        tpTo_hltIterL3IOFromL1_pt = 0;
        tpTo_hltIterL3IOFromL1_eta = 0;
        tpTo_hltIterL3IOFromL1_phi = 0;
        tpTo_hltIterL3IOFromL1_parentVx = 0;
        tpTo_hltIterL3IOFromL1_parentVy = 0;
        tpTo_hltIterL3IOFromL1_parentVz = 0;
        tpTo_hltIterL3IOFromL1_status = 0;
        tpTo_hltIterL3IOFromL1_numberOfHits = 0;
        tpTo_hltIterL3IOFromL1_numberOfTrackerHits = 0;
        tpTo_hltIterL3IOFromL1_numberOfTrackerLayers = 0;
        tpTo_hltIterL3IOFromL1_gen_charge = 0;
        tpTo_hltIterL3IOFromL1_gen_pdgId = 0;
        tpTo_hltIterL3IOFromL1_gen_pt = 0;
        tpTo_hltIterL3IOFromL1_gen_eta = 0;
        tpTo_hltIterL3IOFromL1_gen_phi = 0;
        tpTo_hltIterL3IOFromL1_bestMatchTrk_pt = 0;
        tpTo_hltIterL3IOFromL1_bestMatchTrk_eta = 0;
        tpTo_hltIterL3IOFromL1_bestMatchTrk_phi = 0;
        tpTo_hltIterL3IOFromL1_bestMatchTrk_charge = 0;
        tpTo_hltIterL3IOFromL1_bestMatchTrk_quality = 0;
        tpTo_hltIterL3IOFromL1_bestMatchTrk_NValidHits = 0;
        hltIterL3MuonsNoID_pt = 0;
        hltIterL3MuonsNoID_ptError = 0;
        hltIterL3MuonsNoID_eta = 0;
        hltIterL3MuonsNoID_phi = 0;
        hltIterL3MuonsNoID_charge = 0;
        hltIterL3MuonsNoID_matchedL3 = 0;
        hltIterL3MuonsNoID_matchedL3NoId = 0;
        hltIterL3MuonsNoID_bestMatchTP_charge = 0;
        hltIterL3MuonsNoID_bestMatchTP_pdgId = 0;
        hltIterL3MuonsNoID_bestMatchTP_energy = 0;
        hltIterL3MuonsNoID_bestMatchTP_pt = 0;
        hltIterL3MuonsNoID_bestMatchTP_eta = 0;
        hltIterL3MuonsNoID_bestMatchTP_phi = 0;
        hltIterL3MuonsNoID_bestMatchTP_parentVx = 0;
        hltIterL3MuonsNoID_bestMatchTP_parentVy = 0;
        hltIterL3MuonsNoID_bestMatchTP_parentVz = 0;
        hltIterL3MuonsNoID_bestMatchTP_status = 0;
        hltIterL3MuonsNoID_bestMatchTP_numberOfHits = 0;
        hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerHits = 0;
        hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerLayers = 0;
        hltIterL3MuonsNoID_bestMatchTP_sharedFraction = 0;
        hltIterL3MuonsNoID_matchedTPsize = 0;
        hltIterL3MuonsNoID_mva0 = 0;
        hltIterL3MuonsNoID_mva1 = 0;
        hltIterL3MuonsNoID_mva2 = 0;
        hltIterL3MuonsNoID_mva3 = 0;
        tpTo_hltIterL3MuonsNoID_charge = 0;
        tpTo_hltIterL3MuonsNoID_pdgId = 0;
        tpTo_hltIterL3MuonsNoID_energy = 0;
        tpTo_hltIterL3MuonsNoID_pt = 0;
        tpTo_hltIterL3MuonsNoID_eta = 0;
        tpTo_hltIterL3MuonsNoID_phi = 0;
        tpTo_hltIterL3MuonsNoID_parentVx = 0;
        tpTo_hltIterL3MuonsNoID_parentVy = 0;
        tpTo_hltIterL3MuonsNoID_parentVz = 0;
        tpTo_hltIterL3MuonsNoID_status = 0;
        tpTo_hltIterL3MuonsNoID_numberOfHits = 0;
        tpTo_hltIterL3MuonsNoID_numberOfTrackerHits = 0;
        tpTo_hltIterL3MuonsNoID_numberOfTrackerLayers = 0;
        tpTo_hltIterL3MuonsNoID_gen_charge = 0;
        tpTo_hltIterL3MuonsNoID_gen_pdgId = 0;
        tpTo_hltIterL3MuonsNoID_gen_pt = 0;
        tpTo_hltIterL3MuonsNoID_gen_eta = 0;
        tpTo_hltIterL3MuonsNoID_gen_phi = 0;
        tpTo_hltIterL3MuonsNoID_bestMatchTrk_pt = 0;
        tpTo_hltIterL3MuonsNoID_bestMatchTrk_eta = 0;
        tpTo_hltIterL3MuonsNoID_bestMatchTrk_phi = 0;
        tpTo_hltIterL3MuonsNoID_bestMatchTrk_charge = 0;
        tpTo_hltIterL3MuonsNoID_bestMatchTrk_quality = 0;
        tpTo_hltIterL3MuonsNoID_bestMatchTrk_NValidHits = 0;
        hltIterL3Muons_pt = 0;
        hltIterL3Muons_ptError = 0;
        hltIterL3Muons_eta = 0;
        hltIterL3Muons_phi = 0;
        hltIterL3Muons_charge = 0;
        hltIterL3Muons_matchedL3 = 0;
        hltIterL3Muons_matchedL3NoId = 0;
        hltIterL3Muons_bestMatchTP_charge = 0;
        hltIterL3Muons_bestMatchTP_pdgId = 0;
        hltIterL3Muons_bestMatchTP_energy = 0;
        hltIterL3Muons_bestMatchTP_pt = 0;
        hltIterL3Muons_bestMatchTP_eta = 0;
        hltIterL3Muons_bestMatchTP_phi = 0;
        hltIterL3Muons_bestMatchTP_parentVx = 0;
        hltIterL3Muons_bestMatchTP_parentVy = 0;
        hltIterL3Muons_bestMatchTP_parentVz = 0;
        hltIterL3Muons_bestMatchTP_status = 0;
        hltIterL3Muons_bestMatchTP_numberOfHits = 0;
        hltIterL3Muons_bestMatchTP_numberOfTrackerHits = 0;
        hltIterL3Muons_bestMatchTP_numberOfTrackerLayers = 0;
        hltIterL3Muons_bestMatchTP_sharedFraction = 0;
        hltIterL3Muons_matchedTPsize = 0;
        hltIterL3Muons_mva0 = 0;
        hltIterL3Muons_mva1 = 0;
        hltIterL3Muons_mva2 = 0;
        hltIterL3Muons_mva3 = 0;
        tpTo_hltIterL3Muons_charge = 0;
        tpTo_hltIterL3Muons_pdgId = 0;
        tpTo_hltIterL3Muons_energy = 0;
        tpTo_hltIterL3Muons_pt = 0;
        tpTo_hltIterL3Muons_eta = 0;
        tpTo_hltIterL3Muons_phi = 0;
        tpTo_hltIterL3Muons_parentVx = 0;
        tpTo_hltIterL3Muons_parentVy = 0;
        tpTo_hltIterL3Muons_parentVz = 0;
        tpTo_hltIterL3Muons_status = 0;
        tpTo_hltIterL3Muons_numberOfHits = 0;
        tpTo_hltIterL3Muons_numberOfTrackerHits = 0;
        tpTo_hltIterL3Muons_numberOfTrackerLayers = 0;
        tpTo_hltIterL3Muons_gen_charge = 0;
        tpTo_hltIterL3Muons_gen_pdgId = 0;
        tpTo_hltIterL3Muons_gen_pt = 0;
        tpTo_hltIterL3Muons_gen_eta = 0;
        tpTo_hltIterL3Muons_gen_phi = 0;
        tpTo_hltIterL3Muons_bestMatchTrk_pt = 0;
        tpTo_hltIterL3Muons_bestMatchTrk_eta = 0;
        tpTo_hltIterL3Muons_bestMatchTrk_phi = 0;
        tpTo_hltIterL3Muons_bestMatchTrk_charge = 0;
        tpTo_hltIterL3Muons_bestMatchTrk_quality = 0;
        tpTo_hltIterL3Muons_bestMatchTrk_NValidHits = 0;

    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    isIterL3 = false;
    if(fChain->GetListOfBranches()->FindObject( "hltIterL3Muons_pt" ))
        isIterL3 = true;

    // cout << "MuonHLTNtuple::Init   isIterL3 = " << isIterL3 << endl;

    // -- Set branch addresses -- //
        fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
        fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
        fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
        fChain->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
        fChain->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
        fChain->SetBranchAddress("trk_rInv", &trk_rInv, &b_trk_rInv);
        fChain->SetBranchAddress("trk_tanL", &trk_tanL, &b_trk_tanL);
        fChain->SetBranchAddress("trk_MVA1", &trk_MVA1, &b_trk_MVA1);
        fChain->SetBranchAddress("trk_MVA2", &trk_MVA2, &b_trk_MVA2);
        fChain->SetBranchAddress("trk_MVA3", &trk_MVA3, &b_trk_MVA3);
        fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
        fChain->SetBranchAddress("trk_bendchi2", &trk_bendchi2, &b_trk_bendchi2);
        fChain->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
        fChain->SetBranchAddress("trk_lhits", &trk_lhits, &b_trk_lhits);
        fChain->SetBranchAddress("trk_dhits", &trk_dhits, &b_trk_dhits);
        fChain->SetBranchAddress("trk_seed", &trk_seed, &b_trk_seed);
        fChain->SetBranchAddress("trk_phiSector", &trk_phiSector, &b_trk_phiSector);
        fChain->SetBranchAddress("trk_genuine", &trk_genuine, &b_trk_genuine);
        fChain->SetBranchAddress("trk_loose", &trk_loose, &b_trk_loose);
        fChain->SetBranchAddress("trk_unknown", &trk_unknown, &b_trk_unknown);
        fChain->SetBranchAddress("trk_combinatoric", &trk_combinatoric, &b_trk_combinatoric);
        fChain->SetBranchAddress("trk_fake", &trk_fake, &b_trk_fake);
        fChain->SetBranchAddress("trk_matchtp_pdgid", &trk_matchtp_pdgid, &b_trk_matchtp_pdgid);
        fChain->SetBranchAddress("trk_matchtp_pt", &trk_matchtp_pt, &b_trk_matchtp_pt);
        fChain->SetBranchAddress("trk_matchtp_eta", &trk_matchtp_eta, &b_trk_matchtp_eta);
        fChain->SetBranchAddress("trk_matchtp_phi", &trk_matchtp_phi, &b_trk_matchtp_phi);
        fChain->SetBranchAddress("trk_matchtp_z0", &trk_matchtp_z0, &b_trk_matchtp_z0);
        fChain->SetBranchAddress("trk_matchtp_dxy", &trk_matchtp_dxy, &b_trk_matchtp_dxy);
        fChain->SetBranchAddress("stub_x", &stub_x, &b_stub_x);
        fChain->SetBranchAddress("stub_y", &stub_y, &b_stub_y);
        fChain->SetBranchAddress("stub_z", &stub_z, &b_stub_z);
        fChain->SetBranchAddress("stub_isBarrel", &stub_isBarrel, &b_stub_isBarrel);
        fChain->SetBranchAddress("stub_layer", &stub_layer, &b_stub_layer);
        fChain->SetBranchAddress("L1TkMu_pt", &L1TkMu_pt, &b_L1TkMu_pt);
        fChain->SetBranchAddress("L1TkMu_eta", &L1TkMu_eta, &b_L1TkMu_eta);
        fChain->SetBranchAddress("L1TkMu_phi", &L1TkMu_phi, &b_L1TkMu_phi);
        fChain->SetBranchAddress("L1TkMu_trkIsol", &L1TkMu_trkIsol, &b_L1TkMu_trkIsol);
        fChain->SetBranchAddress("L1TkMu_trkzVtx", &L1TkMu_trkzVtx, &b_L1TkMu_trkzVtx);
        fChain->SetBranchAddress("L1TkMu_dR", &L1TkMu_dR, &b_L1TkMu_dR);
        fChain->SetBranchAddress("L1TkMu_nTracksMatched", &L1TkMu_nTracksMatched, &b_L1TkMu_nTracksMatched);
        fChain->SetBranchAddress("L1TkMu_trackCurvature", &L1TkMu_trackCurvature, &b_L1TkMu_trackCurvature);
        fChain->SetBranchAddress("L1TkMu_quality", &L1TkMu_quality, &b_L1TkMu_quality);
        fChain->SetBranchAddress("L1TkMu_pattern", &L1TkMu_pattern, &b_L1TkMu_pattern);
        fChain->SetBranchAddress("L1TkMu_muonDetector", &L1TkMu_muonDetector, &b_L1TkMu_muonDetector);
        fChain->SetBranchAddress("L1TkMu_TTTpointer", &L1TkMu_TTTpointer, &b_L1TkMu_TTTpointer);
        fChain->SetBranchAddress("L1TkMu_muRefHwPt", &L1TkMu_muRefHwPt, &b_L1TkMu_muRefHwPt);
        fChain->SetBranchAddress("L1TkMu_muRefHwDXY", &L1TkMu_muRefHwDXY, &b_L1TkMu_muRefHwDXY);
        fChain->SetBranchAddress("L1TkMu_muRefHwEta", &L1TkMu_muRefHwEta, &b_L1TkMu_muRefHwEta);
        fChain->SetBranchAddress("L1TkMu_muRefHwPhi", &L1TkMu_muRefHwPhi, &b_L1TkMu_muRefHwPhi);
        fChain->SetBranchAddress("L1TkMu_muRefHwSign", &L1TkMu_muRefHwSign, &b_L1TkMu_muRefHwSign);
        fChain->SetBranchAddress("L1TkMu_muRefHwSignValid", &L1TkMu_muRefHwSignValid, &b_L1TkMu_muRefHwSignValid);
        fChain->SetBranchAddress("L1TkMu_muRefHwQual", &L1TkMu_muRefHwQual, &b_L1TkMu_muRefHwQual);
        fChain->SetBranchAddress("isRealData", &isRealData, &b_isRealData);
        fChain->SetBranchAddress("runNum", &runNum, &b_runNum);
        fChain->SetBranchAddress("lumiBlockNum", &lumiBlockNum, &b_lumiBlockNum);
        fChain->SetBranchAddress("eventNum", &eventNum, &b_eventNum);
        fChain->SetBranchAddress("nVertex", &nVertex, &b_nVertex);
        fChain->SetBranchAddress("bunchID", &bunchID, &b_bunchID);
        fChain->SetBranchAddress("instLumi", &instLumi, &b_instLumi);
        fChain->SetBranchAddress("dataPU", &dataPU, &b_dataPU);
        fChain->SetBranchAddress("dataPURMS", &dataPURMS, &b_dataPURMS);
        fChain->SetBranchAddress("bunchLumi", &bunchLumi, &b_bunchLumi);
        fChain->SetBranchAddress("offlineInstLumi", &offlineInstLumi, &b_offlineInstLumi);
        fChain->SetBranchAddress("offlineDataPU", &offlineDataPU, &b_offlineDataPU);
        fChain->SetBranchAddress("offlineDataPURMS", &offlineDataPURMS, &b_offlineDataPURMS);
        fChain->SetBranchAddress("offlineBunchLumi", &offlineBunchLumi, &b_offlineBunchLumi);
        fChain->SetBranchAddress("truePU", &truePU, &b_truePU);
        fChain->SetBranchAddress("qScale", &qScale, &b_qScale);
        fChain->SetBranchAddress("genEventWeight", &genEventWeight, &b_genEventWeight);
        fChain->SetBranchAddress("PU_pT_hats", &PU_pT_hats, &b_PU_pT_hats);
        fChain->SetBranchAddress("nGenParticle", &nGenParticle, &b_nGenParticle);
        fChain->SetBranchAddress("genParticle_ID", genParticle_ID, &b_genParticle_ID);
        fChain->SetBranchAddress("genParticle_status", genParticle_status, &b_genParticle_status);
        fChain->SetBranchAddress("genParticle_mother", genParticle_mother, &b_genParticle_mother);
        fChain->SetBranchAddress("genParticle_pt", genParticle_pt, &b_genParticle_pt);
        fChain->SetBranchAddress("genParticle_eta", genParticle_eta, &b_genParticle_eta);
        fChain->SetBranchAddress("genParticle_phi", genParticle_phi, &b_genParticle_phi);
        fChain->SetBranchAddress("genParticle_px", genParticle_px, &b_genParticle_px);
        fChain->SetBranchAddress("genParticle_py", genParticle_py, &b_genParticle_py);
        fChain->SetBranchAddress("genParticle_pz", genParticle_pz, &b_genParticle_pz);
        fChain->SetBranchAddress("genParticle_energy", genParticle_energy, &b_genParticle_energy);
        fChain->SetBranchAddress("genParticle_charge", genParticle_charge, &b_genParticle_charge);
        fChain->SetBranchAddress("genParticle_isPrompt", genParticle_isPrompt, &b_genParticle_isPrompt);
        fChain->SetBranchAddress("genParticle_isPromptFinalState", genParticle_isPromptFinalState, &b_genParticle_isPromptFinalState);
        fChain->SetBranchAddress("genParticle_isTauDecayProduct", genParticle_isTauDecayProduct, &b_genParticle_isTauDecayProduct);
        fChain->SetBranchAddress("genParticle_isPromptTauDecayProduct", genParticle_isPromptTauDecayProduct, &b_genParticle_isPromptTauDecayProduct);
        fChain->SetBranchAddress("genParticle_isDirectPromptTauDecayProductFinalState", genParticle_isDirectPromptTauDecayProductFinalState, &b_genParticle_isDirectPromptTauDecayProductFinalState);
        fChain->SetBranchAddress("genParticle_isHardProcess", genParticle_isHardProcess, &b_genParticle_isHardProcess);
        fChain->SetBranchAddress("genParticle_isLastCopy", genParticle_isLastCopy, &b_genParticle_isLastCopy);
        fChain->SetBranchAddress("genParticle_isLastCopyBeforeFSR", genParticle_isLastCopyBeforeFSR, &b_genParticle_isLastCopyBeforeFSR);
        fChain->SetBranchAddress("genParticle_isPromptDecayed", genParticle_isPromptDecayed, &b_genParticle_isPromptDecayed);
        fChain->SetBranchAddress("genParticle_isDecayedLeptonHadron", genParticle_isDecayedLeptonHadron, &b_genParticle_isDecayedLeptonHadron);
        fChain->SetBranchAddress("genParticle_fromHardProcessBeforeFSR", genParticle_fromHardProcessBeforeFSR, &b_genParticle_fromHardProcessBeforeFSR);
        fChain->SetBranchAddress("genParticle_fromHardProcessDecayed", genParticle_fromHardProcessDecayed, &b_genParticle_fromHardProcessDecayed);
        fChain->SetBranchAddress("genParticle_fromHardProcessFinalState", genParticle_fromHardProcessFinalState, &b_genParticle_fromHardProcessFinalState);
        fChain->SetBranchAddress("genParticle_isMostlyLikePythia6Status3", genParticle_isMostlyLikePythia6Status3, &b_genParticle_isMostlyLikePythia6Status3);
        fChain->SetBranchAddress("vec_firedTrigger", &vec_firedTrigger, &b_vec_firedTrigger);
        fChain->SetBranchAddress("vec_filterName", &vec_filterName, &b_vec_filterName);
        fChain->SetBranchAddress("vec_HLTObj_pt", &vec_HLTObj_pt, &b_vec_HLTObj_pt);
        fChain->SetBranchAddress("vec_HLTObj_eta", &vec_HLTObj_eta, &b_vec_HLTObj_eta);
        fChain->SetBranchAddress("vec_HLTObj_phi", &vec_HLTObj_phi, &b_vec_HLTObj_phi);
        fChain->SetBranchAddress("vec_myFiredTrigger", &vec_myFiredTrigger, &b_vec_myFiredTrigger);
        fChain->SetBranchAddress("vec_myFilterName", &vec_myFilterName, &b_vec_myFilterName);
        fChain->SetBranchAddress("vec_myHLTObj_pt", &vec_myHLTObj_pt, &b_vec_myHLTObj_pt);
        fChain->SetBranchAddress("vec_myHLTObj_eta", &vec_myHLTObj_eta, &b_vec_myHLTObj_eta);
        fChain->SetBranchAddress("vec_myHLTObj_phi", &vec_myHLTObj_phi, &b_vec_myHLTObj_phi);
        fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
        fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
        fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
        fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
        fChain->SetBranchAddress("muon_px", muon_px, &b_muon_px);
        fChain->SetBranchAddress("muon_py", muon_py, &b_muon_py);
        fChain->SetBranchAddress("muon_pz", muon_pz, &b_muon_pz);
        fChain->SetBranchAddress("muon_dB", muon_dB, &b_muon_dB);
        fChain->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
        fChain->SetBranchAddress("muon_isGLB", muon_isGLB, &b_muon_isGLB);
        fChain->SetBranchAddress("muon_isSTA", muon_isSTA, &b_muon_isSTA);
        fChain->SetBranchAddress("muon_isTRK", muon_isTRK, &b_muon_isTRK);
        fChain->SetBranchAddress("muon_isPF", muon_isPF, &b_muon_isPF);
        fChain->SetBranchAddress("muon_isTight", muon_isTight, &b_muon_isTight);
        fChain->SetBranchAddress("muon_isMedium", muon_isMedium, &b_muon_isMedium);
        fChain->SetBranchAddress("muon_isLoose", muon_isLoose, &b_muon_isLoose);
        fChain->SetBranchAddress("muon_isHighPt", muon_isHighPt, &b_muon_isHighPt);
        fChain->SetBranchAddress("muon_isHighPtNew", muon_isHighPtNew, &b_muon_isHighPtNew);
        fChain->SetBranchAddress("muon_isSoft", muon_isSoft, &b_muon_isSoft);
        fChain->SetBranchAddress("muon_isLooseTriggerMuon", muon_isLooseTriggerMuon, &b_muon_isLooseTriggerMuon);
        fChain->SetBranchAddress("muon_isME0Muon", muon_isME0Muon, &b_muon_isME0Muon);
        fChain->SetBranchAddress("muon_isGEMMuon", muon_isGEMMuon, &b_muon_isGEMMuon);
        fChain->SetBranchAddress("muon_isRPCMuon", muon_isRPCMuon, &b_muon_isRPCMuon);
        fChain->SetBranchAddress("muon_isGoodMuon_TMOneStationTight", muon_isGoodMuon_TMOneStationTight, &b_muon_isGoodMuon_TMOneStationTight);
        fChain->SetBranchAddress("muon_iso03_sumPt", muon_iso03_sumPt, &b_muon_iso03_sumPt);
        fChain->SetBranchAddress("muon_iso03_hadEt", muon_iso03_hadEt, &b_muon_iso03_hadEt);
        fChain->SetBranchAddress("muon_iso03_emEt", muon_iso03_emEt, &b_muon_iso03_emEt);
        fChain->SetBranchAddress("muon_PFIso03_charged", muon_PFIso03_charged, &b_muon_PFIso03_charged);
        fChain->SetBranchAddress("muon_PFIso03_neutral", muon_PFIso03_neutral, &b_muon_PFIso03_neutral);
        fChain->SetBranchAddress("muon_PFIso03_photon", muon_PFIso03_photon, &b_muon_PFIso03_photon);
        fChain->SetBranchAddress("muon_PFIso03_sumPU", muon_PFIso03_sumPU, &b_muon_PFIso03_sumPU);
        fChain->SetBranchAddress("muon_PFIso04_charged", muon_PFIso04_charged, &b_muon_PFIso04_charged);
        fChain->SetBranchAddress("muon_PFIso04_neutral", muon_PFIso04_neutral, &b_muon_PFIso04_neutral);
        fChain->SetBranchAddress("muon_PFIso04_photon", muon_PFIso04_photon, &b_muon_PFIso04_photon);
        fChain->SetBranchAddress("muon_PFIso04_sumPU", muon_PFIso04_sumPU, &b_muon_PFIso04_sumPU);
        fChain->SetBranchAddress("muon_PFCluster03_ECAL", muon_PFCluster03_ECAL, &b_muon_PFCluster03_ECAL);
        fChain->SetBranchAddress("muon_PFCluster03_HCAL", muon_PFCluster03_HCAL, &b_muon_PFCluster03_HCAL);
        fChain->SetBranchAddress("muon_PFCluster04_ECAL", muon_PFCluster04_ECAL, &b_muon_PFCluster04_ECAL);
        fChain->SetBranchAddress("muon_PFCluster04_HCAL", muon_PFCluster04_HCAL, &b_muon_PFCluster04_HCAL);
        fChain->SetBranchAddress("muon_normChi2_global", muon_normChi2_global, &b_muon_normChi2_global);
        fChain->SetBranchAddress("muon_nTrackerHit_global", muon_nTrackerHit_global, &b_muon_nTrackerHit_global);
        fChain->SetBranchAddress("muon_nTrackerLayer_global", muon_nTrackerLayer_global, &b_muon_nTrackerLayer_global);
        fChain->SetBranchAddress("muon_nPixelHit_global", muon_nPixelHit_global, &b_muon_nPixelHit_global);
        fChain->SetBranchAddress("muon_nMuonHit_global", muon_nMuonHit_global, &b_muon_nMuonHit_global);
        fChain->SetBranchAddress("muon_normChi2_inner", muon_normChi2_inner, &b_muon_normChi2_inner);
        fChain->SetBranchAddress("muon_nTrackerHit_inner", muon_nTrackerHit_inner, &b_muon_nTrackerHit_inner);
        fChain->SetBranchAddress("muon_nTrackerLayer_inner", muon_nTrackerLayer_inner, &b_muon_nTrackerLayer_inner);
        fChain->SetBranchAddress("muon_nPixelHit_inner", muon_nPixelHit_inner, &b_muon_nPixelHit_inner);
        fChain->SetBranchAddress("muon_pt_tuneP", muon_pt_tuneP, &b_muon_pt_tuneP);
        fChain->SetBranchAddress("muon_ptError_tuneP", muon_ptError_tuneP, &b_muon_ptError_tuneP);
        fChain->SetBranchAddress("muon_dxyVTX_best", muon_dxyVTX_best, &b_muon_dxyVTX_best);
        fChain->SetBranchAddress("muon_dzVTX_best", muon_dzVTX_best, &b_muon_dzVTX_best);
        fChain->SetBranchAddress("muon_nMatchedStation", muon_nMatchedStation, &b_muon_nMatchedStation);
        fChain->SetBranchAddress("muon_nMatchedRPCLayer", muon_nMatchedRPCLayer, &b_muon_nMatchedRPCLayer);
        fChain->SetBranchAddress("muon_stationMask", muon_stationMask, &b_muon_stationMask);
        fChain->SetBranchAddress("muon_expectedNnumberOfMatchedStations", muon_expectedNnumberOfMatchedStations, &b_muon_expectedNnumberOfMatchedStations);
        fChain->SetBranchAddress("nL3Muon", &nL3Muon, &b_nL3Muon);
        fChain->SetBranchAddress("L3Muon_pt", L3Muon_pt, &b_L3Muon_pt);
        fChain->SetBranchAddress("L3Muon_eta", L3Muon_eta, &b_L3Muon_eta);
        fChain->SetBranchAddress("L3Muon_phi", L3Muon_phi, &b_L3Muon_phi);
        fChain->SetBranchAddress("L3Muon_charge", L3Muon_charge, &b_L3Muon_charge);
        fChain->SetBranchAddress("L3Muon_trkPt", L3Muon_trkPt, &b_L3Muon_trkPt);
        fChain->SetBranchAddress("nL2Muon", &nL2Muon, &b_nL2Muon);
        fChain->SetBranchAddress("L2Muon_pt", &L2Muon_pt, &b_L2Muon_pt);
        fChain->SetBranchAddress("L2Muon_eta", &L2Muon_eta, &b_L2Muon_eta);
        fChain->SetBranchAddress("L2Muon_phi", &L2Muon_phi, &b_L2Muon_phi);
        fChain->SetBranchAddress("L2Muon_charge", &L2Muon_charge, &b_L2Muon_charge);
        fChain->SetBranchAddress("L2Muon_trkPt", &L2Muon_trkPt, &b_L2Muon_trkPt);
        fChain->SetBranchAddress("nTkMuon", &nTkMuon, &b_nTkMuon);
        fChain->SetBranchAddress("TkMuon_pt", &TkMuon_pt, &b_TkMuon_pt);
        fChain->SetBranchAddress("TkMuon_eta", &TkMuon_eta, &b_TkMuon_eta);
        fChain->SetBranchAddress("TkMuon_phi", &TkMuon_phi, &b_TkMuon_phi);
        fChain->SetBranchAddress("TkMuon_charge", &TkMuon_charge, &b_TkMuon_charge);
        fChain->SetBranchAddress("TkMuon_trkPt", &TkMuon_trkPt, &b_TkMuon_trkPt);
        fChain->SetBranchAddress("nL1Muon", &nL1Muon, &b_nL1Muon);
        fChain->SetBranchAddress("L1Muon_pt", L1Muon_pt, &b_L1Muon_pt);
        fChain->SetBranchAddress("L1Muon_eta", L1Muon_eta, &b_L1Muon_eta);
        fChain->SetBranchAddress("L1Muon_phi", L1Muon_phi, &b_L1Muon_phi);
        fChain->SetBranchAddress("L1Muon_charge", L1Muon_charge, &b_L1Muon_charge);
        fChain->SetBranchAddress("L1Muon_quality", L1Muon_quality, &b_L1Muon_quality);
        fChain->SetBranchAddress("L1Muon_etaAtVtx", L1Muon_etaAtVtx, &b_L1Muon_etaAtVtx);
        fChain->SetBranchAddress("L1Muon_phiAtVtx", L1Muon_phiAtVtx, &b_L1Muon_phiAtVtx);
        fChain->SetBranchAddress("nIterL3OI", &nIterL3OI, &b_nIterL3OI);
        fChain->SetBranchAddress("iterL3OI_inner_pt", &iterL3OI_inner_pt, &b_iterL3OI_inner_pt);
        fChain->SetBranchAddress("iterL3OI_inner_eta", &iterL3OI_inner_eta, &b_iterL3OI_inner_eta);
        fChain->SetBranchAddress("iterL3OI_inner_phi", &iterL3OI_inner_phi, &b_iterL3OI_inner_phi);
        fChain->SetBranchAddress("iterL3OI_inner_charge", &iterL3OI_inner_charge, &b_iterL3OI_inner_charge);
        fChain->SetBranchAddress("iterL3OI_outer_pt", &iterL3OI_outer_pt, &b_iterL3OI_outer_pt);
        fChain->SetBranchAddress("iterL3OI_outer_eta", &iterL3OI_outer_eta, &b_iterL3OI_outer_eta);
        fChain->SetBranchAddress("iterL3OI_outer_phi", &iterL3OI_outer_phi, &b_iterL3OI_outer_phi);
        fChain->SetBranchAddress("iterL3OI_outer_charge", &iterL3OI_outer_charge, &b_iterL3OI_outer_charge);
        fChain->SetBranchAddress("iterL3OI_global_pt", &iterL3OI_global_pt, &b_iterL3OI_global_pt);
        fChain->SetBranchAddress("iterL3OI_global_eta", &iterL3OI_global_eta, &b_iterL3OI_global_eta);
        fChain->SetBranchAddress("iterL3OI_global_phi", &iterL3OI_global_phi, &b_iterL3OI_global_phi);
        fChain->SetBranchAddress("iterL3OI_global_charge", &iterL3OI_global_charge, &b_iterL3OI_global_charge);
        fChain->SetBranchAddress("nIterL3IOFromL2", &nIterL3IOFromL2, &b_nIterL3IOFromL2);
        fChain->SetBranchAddress("iterL3IOFromL2_inner_pt", &iterL3IOFromL2_inner_pt, &b_iterL3IOFromL2_inner_pt);
        fChain->SetBranchAddress("iterL3IOFromL2_inner_eta", &iterL3IOFromL2_inner_eta, &b_iterL3IOFromL2_inner_eta);
        fChain->SetBranchAddress("iterL3IOFromL2_inner_phi", &iterL3IOFromL2_inner_phi, &b_iterL3IOFromL2_inner_phi);
        fChain->SetBranchAddress("iterL3IOFromL2_inner_charge", &iterL3IOFromL2_inner_charge, &b_iterL3IOFromL2_inner_charge);
        fChain->SetBranchAddress("iterL3IOFromL2_outer_pt", &iterL3IOFromL2_outer_pt, &b_iterL3IOFromL2_outer_pt);
        fChain->SetBranchAddress("iterL3IOFromL2_outer_eta", &iterL3IOFromL2_outer_eta, &b_iterL3IOFromL2_outer_eta);
        fChain->SetBranchAddress("iterL3IOFromL2_outer_phi", &iterL3IOFromL2_outer_phi, &b_iterL3IOFromL2_outer_phi);
        fChain->SetBranchAddress("iterL3IOFromL2_outer_charge", &iterL3IOFromL2_outer_charge, &b_iterL3IOFromL2_outer_charge);
        fChain->SetBranchAddress("iterL3IOFromL2_global_pt", &iterL3IOFromL2_global_pt, &b_iterL3IOFromL2_global_pt);
        fChain->SetBranchAddress("iterL3IOFromL2_global_eta", &iterL3IOFromL2_global_eta, &b_iterL3IOFromL2_global_eta);
        fChain->SetBranchAddress("iterL3IOFromL2_global_phi", &iterL3IOFromL2_global_phi, &b_iterL3IOFromL2_global_phi);
        fChain->SetBranchAddress("iterL3IOFromL2_global_charge", &iterL3IOFromL2_global_charge, &b_iterL3IOFromL2_global_charge);
        fChain->SetBranchAddress("nIterL3IOFromL1", &nIterL3IOFromL1, &b_nIterL3IOFromL1);
        fChain->SetBranchAddress("iterL3IOFromL1_pt", iterL3IOFromL1_pt, &b_iterL3IOFromL1_pt);
        fChain->SetBranchAddress("iterL3IOFromL1_eta", iterL3IOFromL1_eta, &b_iterL3IOFromL1_eta);
        fChain->SetBranchAddress("iterL3IOFromL1_phi", iterL3IOFromL1_phi, &b_iterL3IOFromL1_phi);
        fChain->SetBranchAddress("iterL3IOFromL1_charge", iterL3IOFromL1_charge, &b_iterL3IOFromL1_charge);
        fChain->SetBranchAddress("nIterL3FromL2", &nIterL3FromL2, &b_nIterL3FromL2);
        fChain->SetBranchAddress("iterL3FromL2_inner_pt", &iterL3FromL2_inner_pt, &b_iterL3FromL2_inner_pt);
        fChain->SetBranchAddress("iterL3FromL2_inner_eta", &iterL3FromL2_inner_eta, &b_iterL3FromL2_inner_eta);
        fChain->SetBranchAddress("iterL3FromL2_inner_phi", &iterL3FromL2_inner_phi, &b_iterL3FromL2_inner_phi);
        fChain->SetBranchAddress("iterL3FromL2_inner_charge", &iterL3FromL2_inner_charge, &b_iterL3FromL2_inner_charge);
        fChain->SetBranchAddress("iterL3FromL2_outer_pt", &iterL3FromL2_outer_pt, &b_iterL3FromL2_outer_pt);
        fChain->SetBranchAddress("iterL3FromL2_outer_eta", &iterL3FromL2_outer_eta, &b_iterL3FromL2_outer_eta);
        fChain->SetBranchAddress("iterL3FromL2_outer_phi", &iterL3FromL2_outer_phi, &b_iterL3FromL2_outer_phi);
        fChain->SetBranchAddress("iterL3FromL2_outer_charge", &iterL3FromL2_outer_charge, &b_iterL3FromL2_outer_charge);
        fChain->SetBranchAddress("iterL3FromL2_global_pt", &iterL3FromL2_global_pt, &b_iterL3FromL2_global_pt);
        fChain->SetBranchAddress("iterL3FromL2_global_eta", &iterL3FromL2_global_eta, &b_iterL3FromL2_global_eta);
        fChain->SetBranchAddress("iterL3FromL2_global_phi", &iterL3FromL2_global_phi, &b_iterL3FromL2_global_phi);
        fChain->SetBranchAddress("iterL3FromL2_global_charge", &iterL3FromL2_global_charge, &b_iterL3FromL2_global_charge);
        fChain->SetBranchAddress("nIterL3MuonNoID", &nIterL3MuonNoID, &b_nIterL3MuonNoID);
        fChain->SetBranchAddress("iterL3MuonNoID_pt", iterL3MuonNoID_pt, &b_iterL3MuonNoID_pt);
        fChain->SetBranchAddress("iterL3MuonNoID_innerPt", iterL3MuonNoID_innerPt, &b_iterL3MuonNoID_innerPt);
        fChain->SetBranchAddress("iterL3MuonNoID_eta", iterL3MuonNoID_eta, &b_iterL3MuonNoID_eta);
        fChain->SetBranchAddress("iterL3MuonNoID_phi", iterL3MuonNoID_phi, &b_iterL3MuonNoID_phi);
        fChain->SetBranchAddress("iterL3MuonNoID_charge", iterL3MuonNoID_charge, &b_iterL3MuonNoID_charge);
        fChain->SetBranchAddress("iterL3MuonNoID_isGLB", iterL3MuonNoID_isGLB, &b_iterL3MuonNoID_isGLB);
        fChain->SetBranchAddress("iterL3MuonNoID_isSTA", iterL3MuonNoID_isSTA, &b_iterL3MuonNoID_isSTA);
        fChain->SetBranchAddress("iterL3MuonNoID_isTRK", iterL3MuonNoID_isTRK, &b_iterL3MuonNoID_isTRK);
        fChain->SetBranchAddress("nIterL3Muon", &nIterL3Muon, &b_nIterL3Muon);
        fChain->SetBranchAddress("iterL3Muon_pt", iterL3Muon_pt, &b_iterL3Muon_pt);
        fChain->SetBranchAddress("iterL3Muon_innerPt", iterL3Muon_innerPt, &b_iterL3Muon_innerPt);
        fChain->SetBranchAddress("iterL3Muon_eta", iterL3Muon_eta, &b_iterL3Muon_eta);
        fChain->SetBranchAddress("iterL3Muon_phi", iterL3Muon_phi, &b_iterL3Muon_phi);
        fChain->SetBranchAddress("iterL3Muon_charge", iterL3Muon_charge, &b_iterL3Muon_charge);
        fChain->SetBranchAddress("iterL3Muon_isGLB", iterL3Muon_isGLB, &b_iterL3Muon_isGLB);
        fChain->SetBranchAddress("iterL3Muon_isSTA", iterL3Muon_isSTA, &b_iterL3Muon_isSTA);
        fChain->SetBranchAddress("iterL3Muon_isTRK", iterL3Muon_isTRK, &b_iterL3Muon_isTRK);
        fChain->SetBranchAddress("nhltIterL3OISeedsFromL2Muons", &nhltIterL3OISeedsFromL2Muons, &b_nhltIterL3OISeedsFromL2Muons);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_dir", &hltIterL3OISeedsFromL2Muons_dir, &b_hltIterL3OISeedsFromL2Muons_dir);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_detId", &hltIterL3OISeedsFromL2Muons_tsos_detId, &b_hltIterL3OISeedsFromL2Muons_tsos_detId);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_pt", &hltIterL3OISeedsFromL2Muons_tsos_pt, &b_hltIterL3OISeedsFromL2Muons_tsos_pt);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_pt_val", &hltIterL3OISeedsFromL2Muons_tsos_pt_val, &b_hltIterL3OISeedsFromL2Muons_tsos_pt_val);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_eta", &hltIterL3OISeedsFromL2Muons_tsos_eta, &b_hltIterL3OISeedsFromL2Muons_tsos_eta);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_phi", &hltIterL3OISeedsFromL2Muons_tsos_phi, &b_hltIterL3OISeedsFromL2Muons_tsos_phi);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_glob_x", &hltIterL3OISeedsFromL2Muons_tsos_glob_x, &b_hltIterL3OISeedsFromL2Muons_tsos_glob_x);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_glob_y", &hltIterL3OISeedsFromL2Muons_tsos_glob_y, &b_hltIterL3OISeedsFromL2Muons_tsos_glob_y);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_glob_z", &hltIterL3OISeedsFromL2Muons_tsos_glob_z, &b_hltIterL3OISeedsFromL2Muons_tsos_glob_z);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_hasErr", &hltIterL3OISeedsFromL2Muons_tsos_hasErr, &b_hltIterL3OISeedsFromL2Muons_tsos_hasErr);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err0", &hltIterL3OISeedsFromL2Muons_tsos_err0, &b_hltIterL3OISeedsFromL2Muons_tsos_err0);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err1", &hltIterL3OISeedsFromL2Muons_tsos_err1, &b_hltIterL3OISeedsFromL2Muons_tsos_err1);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err2", &hltIterL3OISeedsFromL2Muons_tsos_err2, &b_hltIterL3OISeedsFromL2Muons_tsos_err2);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err3", &hltIterL3OISeedsFromL2Muons_tsos_err3, &b_hltIterL3OISeedsFromL2Muons_tsos_err3);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err4", &hltIterL3OISeedsFromL2Muons_tsos_err4, &b_hltIterL3OISeedsFromL2Muons_tsos_err4);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err5", &hltIterL3OISeedsFromL2Muons_tsos_err5, &b_hltIterL3OISeedsFromL2Muons_tsos_err5);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err6", &hltIterL3OISeedsFromL2Muons_tsos_err6, &b_hltIterL3OISeedsFromL2Muons_tsos_err6);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err7", &hltIterL3OISeedsFromL2Muons_tsos_err7, &b_hltIterL3OISeedsFromL2Muons_tsos_err7);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err8", &hltIterL3OISeedsFromL2Muons_tsos_err8, &b_hltIterL3OISeedsFromL2Muons_tsos_err8);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err9", &hltIterL3OISeedsFromL2Muons_tsos_err9, &b_hltIterL3OISeedsFromL2Muons_tsos_err9);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err10", &hltIterL3OISeedsFromL2Muons_tsos_err10, &b_hltIterL3OISeedsFromL2Muons_tsos_err10);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err11", &hltIterL3OISeedsFromL2Muons_tsos_err11, &b_hltIterL3OISeedsFromL2Muons_tsos_err11);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err12", &hltIterL3OISeedsFromL2Muons_tsos_err12, &b_hltIterL3OISeedsFromL2Muons_tsos_err12);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err13", &hltIterL3OISeedsFromL2Muons_tsos_err13, &b_hltIterL3OISeedsFromL2Muons_tsos_err13);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_err14", &hltIterL3OISeedsFromL2Muons_tsos_err14, &b_hltIterL3OISeedsFromL2Muons_tsos_err14);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_x", &hltIterL3OISeedsFromL2Muons_tsos_x, &b_hltIterL3OISeedsFromL2Muons_tsos_x);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_y", &hltIterL3OISeedsFromL2Muons_tsos_y, &b_hltIterL3OISeedsFromL2Muons_tsos_y);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_dxdz", &hltIterL3OISeedsFromL2Muons_tsos_dxdz, &b_hltIterL3OISeedsFromL2Muons_tsos_dxdz);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_dydz", &hltIterL3OISeedsFromL2Muons_tsos_dydz, &b_hltIterL3OISeedsFromL2Muons_tsos_dydz);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_px", &hltIterL3OISeedsFromL2Muons_tsos_px, &b_hltIterL3OISeedsFromL2Muons_tsos_px);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_py", &hltIterL3OISeedsFromL2Muons_tsos_py, &b_hltIterL3OISeedsFromL2Muons_tsos_py);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_pz", &hltIterL3OISeedsFromL2Muons_tsos_pz, &b_hltIterL3OISeedsFromL2Muons_tsos_pz);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_qbp", &hltIterL3OISeedsFromL2Muons_tsos_qbp, &b_hltIterL3OISeedsFromL2Muons_tsos_qbp);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tsos_charge", &hltIterL3OISeedsFromL2Muons_tsos_charge, &b_hltIterL3OISeedsFromL2Muons_tsos_charge);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_iterL3Matched", &hltIterL3OISeedsFromL2Muons_iterL3Matched, &b_hltIterL3OISeedsFromL2Muons_iterL3Matched);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_iterL3Ref", &hltIterL3OISeedsFromL2Muons_iterL3Ref, &b_hltIterL3OISeedsFromL2Muons_iterL3Ref);
        fChain->SetBranchAddress("hltIterL3OISeedsFromL2Muons_tmpL3Ref", &hltIterL3OISeedsFromL2Muons_tmpL3Ref, &b_hltIterL3OISeedsFromL2Muons_tmpL3Ref);
        fChain->SetBranchAddress("nhltIter0IterL3MuonPixelSeedsFromPixelTracks", &nhltIter0IterL3MuonPixelSeedsFromPixelTracks, &b_nhltIter0IterL3MuonPixelSeedsFromPixelTracks);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_dir", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_dir, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_dir);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_detId", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_detId, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_detId);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt_val", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt_val, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pt_val);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_eta", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_eta, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_eta);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_phi", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_phi, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_phi);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_x", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_x, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_x);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_y", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_y, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_y);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_z", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_z, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_glob_z);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_hasErr", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_hasErr, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_hasErr);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err0", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err0, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err0);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err1", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err1, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err1);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err2", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err2, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err2);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err3", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err3, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err3);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err4", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err4, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err4);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err5", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err5, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err5);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err6", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err6, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err6);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err7", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err7, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err7);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err8", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err8, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err8);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err9", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err9, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err9);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err10", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err10, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err10);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err11", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err11, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err11);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err12", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err12, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err12);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err13", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err13, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err13);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err14", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err14, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_err14);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_x", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_x, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_x);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_y", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_y, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_y);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dxdz", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dxdz, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dxdz);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dydz", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dydz, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_dydz);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_px", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_px, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_px);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_py", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_py, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_py);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pz", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pz, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_pz);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_qbp", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_qbp, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_qbp);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_charge", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_charge, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tsos_charge);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Matched", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Matched, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Matched);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Ref", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Ref, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_iterL3Ref);
        fChain->SetBranchAddress("hltIter0IterL3MuonPixelSeedsFromPixelTracks_tmpL3Ref", &hltIter0IterL3MuonPixelSeedsFromPixelTracks_tmpL3Ref, &b_hltIter0IterL3MuonPixelSeedsFromPixelTracks_tmpL3Ref);
        fChain->SetBranchAddress("nhltIter2IterL3MuonPixelSeeds", &nhltIter2IterL3MuonPixelSeeds, &b_nhltIter2IterL3MuonPixelSeeds);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_dir", &hltIter2IterL3MuonPixelSeeds_dir, &b_hltIter2IterL3MuonPixelSeeds_dir);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_detId", &hltIter2IterL3MuonPixelSeeds_tsos_detId, &b_hltIter2IterL3MuonPixelSeeds_tsos_detId);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_pt", &hltIter2IterL3MuonPixelSeeds_tsos_pt, &b_hltIter2IterL3MuonPixelSeeds_tsos_pt);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_pt_val", &hltIter2IterL3MuonPixelSeeds_tsos_pt_val, &b_hltIter2IterL3MuonPixelSeeds_tsos_pt_val);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_eta", &hltIter2IterL3MuonPixelSeeds_tsos_eta, &b_hltIter2IterL3MuonPixelSeeds_tsos_eta);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_phi", &hltIter2IterL3MuonPixelSeeds_tsos_phi, &b_hltIter2IterL3MuonPixelSeeds_tsos_phi);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_glob_x", &hltIter2IterL3MuonPixelSeeds_tsos_glob_x, &b_hltIter2IterL3MuonPixelSeeds_tsos_glob_x);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_glob_y", &hltIter2IterL3MuonPixelSeeds_tsos_glob_y, &b_hltIter2IterL3MuonPixelSeeds_tsos_glob_y);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_glob_z", &hltIter2IterL3MuonPixelSeeds_tsos_glob_z, &b_hltIter2IterL3MuonPixelSeeds_tsos_glob_z);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_hasErr", &hltIter2IterL3MuonPixelSeeds_tsos_hasErr, &b_hltIter2IterL3MuonPixelSeeds_tsos_hasErr);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err0", &hltIter2IterL3MuonPixelSeeds_tsos_err0, &b_hltIter2IterL3MuonPixelSeeds_tsos_err0);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err1", &hltIter2IterL3MuonPixelSeeds_tsos_err1, &b_hltIter2IterL3MuonPixelSeeds_tsos_err1);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err2", &hltIter2IterL3MuonPixelSeeds_tsos_err2, &b_hltIter2IterL3MuonPixelSeeds_tsos_err2);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err3", &hltIter2IterL3MuonPixelSeeds_tsos_err3, &b_hltIter2IterL3MuonPixelSeeds_tsos_err3);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err4", &hltIter2IterL3MuonPixelSeeds_tsos_err4, &b_hltIter2IterL3MuonPixelSeeds_tsos_err4);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err5", &hltIter2IterL3MuonPixelSeeds_tsos_err5, &b_hltIter2IterL3MuonPixelSeeds_tsos_err5);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err6", &hltIter2IterL3MuonPixelSeeds_tsos_err6, &b_hltIter2IterL3MuonPixelSeeds_tsos_err6);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err7", &hltIter2IterL3MuonPixelSeeds_tsos_err7, &b_hltIter2IterL3MuonPixelSeeds_tsos_err7);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err8", &hltIter2IterL3MuonPixelSeeds_tsos_err8, &b_hltIter2IterL3MuonPixelSeeds_tsos_err8);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err9", &hltIter2IterL3MuonPixelSeeds_tsos_err9, &b_hltIter2IterL3MuonPixelSeeds_tsos_err9);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err10", &hltIter2IterL3MuonPixelSeeds_tsos_err10, &b_hltIter2IterL3MuonPixelSeeds_tsos_err10);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err11", &hltIter2IterL3MuonPixelSeeds_tsos_err11, &b_hltIter2IterL3MuonPixelSeeds_tsos_err11);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err12", &hltIter2IterL3MuonPixelSeeds_tsos_err12, &b_hltIter2IterL3MuonPixelSeeds_tsos_err12);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err13", &hltIter2IterL3MuonPixelSeeds_tsos_err13, &b_hltIter2IterL3MuonPixelSeeds_tsos_err13);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_err14", &hltIter2IterL3MuonPixelSeeds_tsos_err14, &b_hltIter2IterL3MuonPixelSeeds_tsos_err14);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_x", &hltIter2IterL3MuonPixelSeeds_tsos_x, &b_hltIter2IterL3MuonPixelSeeds_tsos_x);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_y", &hltIter2IterL3MuonPixelSeeds_tsos_y, &b_hltIter2IterL3MuonPixelSeeds_tsos_y);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_dxdz", &hltIter2IterL3MuonPixelSeeds_tsos_dxdz, &b_hltIter2IterL3MuonPixelSeeds_tsos_dxdz);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_dydz", &hltIter2IterL3MuonPixelSeeds_tsos_dydz, &b_hltIter2IterL3MuonPixelSeeds_tsos_dydz);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_px", &hltIter2IterL3MuonPixelSeeds_tsos_px, &b_hltIter2IterL3MuonPixelSeeds_tsos_px);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_py", &hltIter2IterL3MuonPixelSeeds_tsos_py, &b_hltIter2IterL3MuonPixelSeeds_tsos_py);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_pz", &hltIter2IterL3MuonPixelSeeds_tsos_pz, &b_hltIter2IterL3MuonPixelSeeds_tsos_pz);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_qbp", &hltIter2IterL3MuonPixelSeeds_tsos_qbp, &b_hltIter2IterL3MuonPixelSeeds_tsos_qbp);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tsos_charge", &hltIter2IterL3MuonPixelSeeds_tsos_charge, &b_hltIter2IterL3MuonPixelSeeds_tsos_charge);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_iterL3Matched", &hltIter2IterL3MuonPixelSeeds_iterL3Matched, &b_hltIter2IterL3MuonPixelSeeds_iterL3Matched);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_iterL3Ref", &hltIter2IterL3MuonPixelSeeds_iterL3Ref, &b_hltIter2IterL3MuonPixelSeeds_iterL3Ref);
        fChain->SetBranchAddress("hltIter2IterL3MuonPixelSeeds_tmpL3Ref", &hltIter2IterL3MuonPixelSeeds_tmpL3Ref, &b_hltIter2IterL3MuonPixelSeeds_tmpL3Ref);
        fChain->SetBranchAddress("nhltIter3IterL3MuonPixelSeeds", &nhltIter3IterL3MuonPixelSeeds, &b_nhltIter3IterL3MuonPixelSeeds);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_dir", &hltIter3IterL3MuonPixelSeeds_dir, &b_hltIter3IterL3MuonPixelSeeds_dir);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_detId", &hltIter3IterL3MuonPixelSeeds_tsos_detId, &b_hltIter3IterL3MuonPixelSeeds_tsos_detId);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_pt", &hltIter3IterL3MuonPixelSeeds_tsos_pt, &b_hltIter3IterL3MuonPixelSeeds_tsos_pt);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_pt_val", &hltIter3IterL3MuonPixelSeeds_tsos_pt_val, &b_hltIter3IterL3MuonPixelSeeds_tsos_pt_val);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_eta", &hltIter3IterL3MuonPixelSeeds_tsos_eta, &b_hltIter3IterL3MuonPixelSeeds_tsos_eta);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_phi", &hltIter3IterL3MuonPixelSeeds_tsos_phi, &b_hltIter3IterL3MuonPixelSeeds_tsos_phi);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_glob_x", &hltIter3IterL3MuonPixelSeeds_tsos_glob_x, &b_hltIter3IterL3MuonPixelSeeds_tsos_glob_x);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_glob_y", &hltIter3IterL3MuonPixelSeeds_tsos_glob_y, &b_hltIter3IterL3MuonPixelSeeds_tsos_glob_y);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_glob_z", &hltIter3IterL3MuonPixelSeeds_tsos_glob_z, &b_hltIter3IterL3MuonPixelSeeds_tsos_glob_z);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_hasErr", &hltIter3IterL3MuonPixelSeeds_tsos_hasErr, &b_hltIter3IterL3MuonPixelSeeds_tsos_hasErr);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err0", &hltIter3IterL3MuonPixelSeeds_tsos_err0, &b_hltIter3IterL3MuonPixelSeeds_tsos_err0);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err1", &hltIter3IterL3MuonPixelSeeds_tsos_err1, &b_hltIter3IterL3MuonPixelSeeds_tsos_err1);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err2", &hltIter3IterL3MuonPixelSeeds_tsos_err2, &b_hltIter3IterL3MuonPixelSeeds_tsos_err2);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err3", &hltIter3IterL3MuonPixelSeeds_tsos_err3, &b_hltIter3IterL3MuonPixelSeeds_tsos_err3);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err4", &hltIter3IterL3MuonPixelSeeds_tsos_err4, &b_hltIter3IterL3MuonPixelSeeds_tsos_err4);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err5", &hltIter3IterL3MuonPixelSeeds_tsos_err5, &b_hltIter3IterL3MuonPixelSeeds_tsos_err5);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err6", &hltIter3IterL3MuonPixelSeeds_tsos_err6, &b_hltIter3IterL3MuonPixelSeeds_tsos_err6);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err7", &hltIter3IterL3MuonPixelSeeds_tsos_err7, &b_hltIter3IterL3MuonPixelSeeds_tsos_err7);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err8", &hltIter3IterL3MuonPixelSeeds_tsos_err8, &b_hltIter3IterL3MuonPixelSeeds_tsos_err8);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err9", &hltIter3IterL3MuonPixelSeeds_tsos_err9, &b_hltIter3IterL3MuonPixelSeeds_tsos_err9);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err10", &hltIter3IterL3MuonPixelSeeds_tsos_err10, &b_hltIter3IterL3MuonPixelSeeds_tsos_err10);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err11", &hltIter3IterL3MuonPixelSeeds_tsos_err11, &b_hltIter3IterL3MuonPixelSeeds_tsos_err11);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err12", &hltIter3IterL3MuonPixelSeeds_tsos_err12, &b_hltIter3IterL3MuonPixelSeeds_tsos_err12);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err13", &hltIter3IterL3MuonPixelSeeds_tsos_err13, &b_hltIter3IterL3MuonPixelSeeds_tsos_err13);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_err14", &hltIter3IterL3MuonPixelSeeds_tsos_err14, &b_hltIter3IterL3MuonPixelSeeds_tsos_err14);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_x", &hltIter3IterL3MuonPixelSeeds_tsos_x, &b_hltIter3IterL3MuonPixelSeeds_tsos_x);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_y", &hltIter3IterL3MuonPixelSeeds_tsos_y, &b_hltIter3IterL3MuonPixelSeeds_tsos_y);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_dxdz", &hltIter3IterL3MuonPixelSeeds_tsos_dxdz, &b_hltIter3IterL3MuonPixelSeeds_tsos_dxdz);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_dydz", &hltIter3IterL3MuonPixelSeeds_tsos_dydz, &b_hltIter3IterL3MuonPixelSeeds_tsos_dydz);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_px", &hltIter3IterL3MuonPixelSeeds_tsos_px, &b_hltIter3IterL3MuonPixelSeeds_tsos_px);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_py", &hltIter3IterL3MuonPixelSeeds_tsos_py, &b_hltIter3IterL3MuonPixelSeeds_tsos_py);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_pz", &hltIter3IterL3MuonPixelSeeds_tsos_pz, &b_hltIter3IterL3MuonPixelSeeds_tsos_pz);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_qbp", &hltIter3IterL3MuonPixelSeeds_tsos_qbp, &b_hltIter3IterL3MuonPixelSeeds_tsos_qbp);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tsos_charge", &hltIter3IterL3MuonPixelSeeds_tsos_charge, &b_hltIter3IterL3MuonPixelSeeds_tsos_charge);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_iterL3Matched", &hltIter3IterL3MuonPixelSeeds_iterL3Matched, &b_hltIter3IterL3MuonPixelSeeds_iterL3Matched);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_iterL3Ref", &hltIter3IterL3MuonPixelSeeds_iterL3Ref, &b_hltIter3IterL3MuonPixelSeeds_iterL3Ref);
        fChain->SetBranchAddress("hltIter3IterL3MuonPixelSeeds_tmpL3Ref", &hltIter3IterL3MuonPixelSeeds_tmpL3Ref, &b_hltIter3IterL3MuonPixelSeeds_tmpL3Ref);
        fChain->SetBranchAddress("nhltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks", &nhltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks, &b_nhltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_dir", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_dir, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_dir);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_detId", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_detId, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_detId);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt_val", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt_val, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pt_val);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_eta", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_eta, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_eta);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_phi", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_phi, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_phi);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_x", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_x, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_x);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_y", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_y, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_y);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_z", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_z, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_glob_z);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_hasErr", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_hasErr, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_hasErr);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err0", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err0, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err0);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err1", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err1, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err1);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err2", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err2, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err2);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err3", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err3, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err3);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err4", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err4, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err4);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err5", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err5, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err5);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err6", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err6, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err6);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err7", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err7, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err7);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err8", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err8, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err8);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err9", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err9, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err9);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err10", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err10, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err10);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err11", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err11, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err11);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err12", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err12, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err12);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err13", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err13, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err13);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err14", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err14, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_err14);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_x", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_x, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_x);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_y", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_y, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_y);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dxdz", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dxdz, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dxdz);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dydz", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dydz, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_dydz);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_px", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_px, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_px);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_py", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_py, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_py);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pz", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pz, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_pz);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_qbp", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_qbp, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_qbp);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_charge", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_charge, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tsos_charge);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Matched", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Matched, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Matched);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Ref", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Ref, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_iterL3Ref);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tmpL3Ref", &hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tmpL3Ref, &b_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_tmpL3Ref);
        fChain->SetBranchAddress("nhltIter2IterL3FromL1MuonPixelSeeds", &nhltIter2IterL3FromL1MuonPixelSeeds, &b_nhltIter2IterL3FromL1MuonPixelSeeds);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_dir", &hltIter2IterL3FromL1MuonPixelSeeds_dir, &b_hltIter2IterL3FromL1MuonPixelSeeds_dir);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_detId", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_detId, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_detId);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt_val", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt_val, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_pt_val);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_eta", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_eta, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_eta);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_phi", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_phi, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_phi);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_x", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_x, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_x);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_y", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_y, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_y);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_z", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_z, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_glob_z);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_hasErr", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_hasErr, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_hasErr);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err0", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err0, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err0);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err1", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err1, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err1);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err2", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err2, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err2);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err3", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err3, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err3);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err4", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err4, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err4);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err5", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err5, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err5);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err6", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err6, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err6);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err7", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err7, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err7);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err8", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err8, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err8);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err9", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err9, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err9);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err10", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err10, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err10);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err11", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err11, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err11);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err12", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err12, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err12);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err13", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err13, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err13);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_err14", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_err14, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_err14);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_x", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_x, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_x);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_y", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_y, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_y);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_dxdz", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_dxdz, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_dxdz);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_dydz", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_dydz, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_dydz);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_px", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_px, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_px);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_py", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_py, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_py);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_pz", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_pz, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_pz);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_qbp", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_qbp, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_qbp);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tsos_charge", &hltIter2IterL3FromL1MuonPixelSeeds_tsos_charge, &b_hltIter2IterL3FromL1MuonPixelSeeds_tsos_charge);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_iterL3Matched", &hltIter2IterL3FromL1MuonPixelSeeds_iterL3Matched, &b_hltIter2IterL3FromL1MuonPixelSeeds_iterL3Matched);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_iterL3Ref", &hltIter2IterL3FromL1MuonPixelSeeds_iterL3Ref, &b_hltIter2IterL3FromL1MuonPixelSeeds_iterL3Ref);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonPixelSeeds_tmpL3Ref", &hltIter2IterL3FromL1MuonPixelSeeds_tmpL3Ref, &b_hltIter2IterL3FromL1MuonPixelSeeds_tmpL3Ref);
        fChain->SetBranchAddress("nhltIter3IterL3FromL1MuonPixelSeeds", &nhltIter3IterL3FromL1MuonPixelSeeds, &b_nhltIter3IterL3FromL1MuonPixelSeeds);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_dir", &hltIter3IterL3FromL1MuonPixelSeeds_dir, &b_hltIter3IterL3FromL1MuonPixelSeeds_dir);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_detId", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_detId, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_detId);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt_val", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt_val, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_pt_val);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_eta", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_eta, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_eta);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_phi", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_phi, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_phi);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_x", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_x, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_x);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_y", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_y, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_y);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_z", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_z, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_glob_z);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_hasErr", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_hasErr, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_hasErr);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err0", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err0, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err0);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err1", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err1, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err1);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err2", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err2, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err2);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err3", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err3, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err3);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err4", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err4, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err4);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err5", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err5, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err5);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err6", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err6, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err6);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err7", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err7, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err7);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err8", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err8, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err8);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err9", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err9, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err9);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err10", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err10, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err10);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err11", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err11, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err11);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err12", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err12, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err12);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err13", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err13, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err13);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_err14", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_err14, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_err14);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_x", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_x, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_x);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_y", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_y, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_y);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_dxdz", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_dxdz, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_dxdz);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_dydz", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_dydz, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_dydz);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_px", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_px, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_px);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_py", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_py, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_py);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_pz", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_pz, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_pz);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_qbp", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_qbp, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_qbp);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tsos_charge", &hltIter3IterL3FromL1MuonPixelSeeds_tsos_charge, &b_hltIter3IterL3FromL1MuonPixelSeeds_tsos_charge);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_iterL3Matched", &hltIter3IterL3FromL1MuonPixelSeeds_iterL3Matched, &b_hltIter3IterL3FromL1MuonPixelSeeds_iterL3Matched);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_iterL3Ref", &hltIter3IterL3FromL1MuonPixelSeeds_iterL3Ref, &b_hltIter3IterL3FromL1MuonPixelSeeds_iterL3Ref);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonPixelSeeds_tmpL3Ref", &hltIter3IterL3FromL1MuonPixelSeeds_tmpL3Ref, &b_hltIter3IterL3FromL1MuonPixelSeeds_tmpL3Ref);
        fChain->SetBranchAddress("nhltIterL3OIMuonTrack", &nhltIterL3OIMuonTrack, &b_nhltIterL3OIMuonTrack);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_pt", &hltIterL3OIMuonTrack_pt, &b_hltIterL3OIMuonTrack_pt);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_ptError", &hltIterL3OIMuonTrack_ptError, &b_hltIterL3OIMuonTrack_ptError);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_eta", &hltIterL3OIMuonTrack_eta, &b_hltIterL3OIMuonTrack_eta);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_phi", &hltIterL3OIMuonTrack_phi, &b_hltIterL3OIMuonTrack_phi);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_charge", &hltIterL3OIMuonTrack_charge, &b_hltIterL3OIMuonTrack_charge);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_matchedL3", &hltIterL3OIMuonTrack_matchedL3, &b_hltIterL3OIMuonTrack_matchedL3);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_matchedL3NoId", &hltIterL3OIMuonTrack_matchedL3NoId, &b_hltIterL3OIMuonTrack_matchedL3NoId);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_charge", &hltIterL3OIMuonTrack_bestMatchTP_charge, &b_hltIterL3OIMuonTrack_bestMatchTP_charge);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_pdgId", &hltIterL3OIMuonTrack_bestMatchTP_pdgId, &b_hltIterL3OIMuonTrack_bestMatchTP_pdgId);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_energy", &hltIterL3OIMuonTrack_bestMatchTP_energy, &b_hltIterL3OIMuonTrack_bestMatchTP_energy);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_pt", &hltIterL3OIMuonTrack_bestMatchTP_pt, &b_hltIterL3OIMuonTrack_bestMatchTP_pt);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_eta", &hltIterL3OIMuonTrack_bestMatchTP_eta, &b_hltIterL3OIMuonTrack_bestMatchTP_eta);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_phi", &hltIterL3OIMuonTrack_bestMatchTP_phi, &b_hltIterL3OIMuonTrack_bestMatchTP_phi);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_parentVx", &hltIterL3OIMuonTrack_bestMatchTP_parentVx, &b_hltIterL3OIMuonTrack_bestMatchTP_parentVx);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_parentVy", &hltIterL3OIMuonTrack_bestMatchTP_parentVy, &b_hltIterL3OIMuonTrack_bestMatchTP_parentVy);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_parentVz", &hltIterL3OIMuonTrack_bestMatchTP_parentVz, &b_hltIterL3OIMuonTrack_bestMatchTP_parentVz);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_status", &hltIterL3OIMuonTrack_bestMatchTP_status, &b_hltIterL3OIMuonTrack_bestMatchTP_status);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_numberOfHits", &hltIterL3OIMuonTrack_bestMatchTP_numberOfHits, &b_hltIterL3OIMuonTrack_bestMatchTP_numberOfHits);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerHits", &hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerHits, &b_hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerHits);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerLayers", &hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerLayers, &b_hltIterL3OIMuonTrack_bestMatchTP_numberOfTrackerLayers);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_bestMatchTP_sharedFraction", &hltIterL3OIMuonTrack_bestMatchTP_sharedFraction, &b_hltIterL3OIMuonTrack_bestMatchTP_sharedFraction);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_matchedTPsize", &hltIterL3OIMuonTrack_matchedTPsize, &b_hltIterL3OIMuonTrack_matchedTPsize);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_mva0", &hltIterL3OIMuonTrack_mva0, &b_hltIterL3OIMuonTrack_mva0);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_mva1", &hltIterL3OIMuonTrack_mva1, &b_hltIterL3OIMuonTrack_mva1);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_mva2", &hltIterL3OIMuonTrack_mva2, &b_hltIterL3OIMuonTrack_mva2);
        fChain->SetBranchAddress("hltIterL3OIMuonTrack_mva3", &hltIterL3OIMuonTrack_mva3, &b_hltIterL3OIMuonTrack_mva3);
        fChain->SetBranchAddress("nhltIter0IterL3MuonTrack", &nhltIter0IterL3MuonTrack, &b_nhltIter0IterL3MuonTrack);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_pt", &hltIter0IterL3MuonTrack_pt, &b_hltIter0IterL3MuonTrack_pt);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_ptError", &hltIter0IterL3MuonTrack_ptError, &b_hltIter0IterL3MuonTrack_ptError);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_eta", &hltIter0IterL3MuonTrack_eta, &b_hltIter0IterL3MuonTrack_eta);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_phi", &hltIter0IterL3MuonTrack_phi, &b_hltIter0IterL3MuonTrack_phi);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_charge", &hltIter0IterL3MuonTrack_charge, &b_hltIter0IterL3MuonTrack_charge);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_matchedL3", &hltIter0IterL3MuonTrack_matchedL3, &b_hltIter0IterL3MuonTrack_matchedL3);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_matchedL3NoId", &hltIter0IterL3MuonTrack_matchedL3NoId, &b_hltIter0IterL3MuonTrack_matchedL3NoId);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_charge", &hltIter0IterL3MuonTrack_bestMatchTP_charge, &b_hltIter0IterL3MuonTrack_bestMatchTP_charge);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_pdgId", &hltIter0IterL3MuonTrack_bestMatchTP_pdgId, &b_hltIter0IterL3MuonTrack_bestMatchTP_pdgId);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_energy", &hltIter0IterL3MuonTrack_bestMatchTP_energy, &b_hltIter0IterL3MuonTrack_bestMatchTP_energy);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_pt", &hltIter0IterL3MuonTrack_bestMatchTP_pt, &b_hltIter0IterL3MuonTrack_bestMatchTP_pt);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_eta", &hltIter0IterL3MuonTrack_bestMatchTP_eta, &b_hltIter0IterL3MuonTrack_bestMatchTP_eta);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_phi", &hltIter0IterL3MuonTrack_bestMatchTP_phi, &b_hltIter0IterL3MuonTrack_bestMatchTP_phi);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_parentVx", &hltIter0IterL3MuonTrack_bestMatchTP_parentVx, &b_hltIter0IterL3MuonTrack_bestMatchTP_parentVx);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_parentVy", &hltIter0IterL3MuonTrack_bestMatchTP_parentVy, &b_hltIter0IterL3MuonTrack_bestMatchTP_parentVy);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_parentVz", &hltIter0IterL3MuonTrack_bestMatchTP_parentVz, &b_hltIter0IterL3MuonTrack_bestMatchTP_parentVz);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_status", &hltIter0IterL3MuonTrack_bestMatchTP_status, &b_hltIter0IterL3MuonTrack_bestMatchTP_status);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_numberOfHits", &hltIter0IterL3MuonTrack_bestMatchTP_numberOfHits, &b_hltIter0IterL3MuonTrack_bestMatchTP_numberOfHits);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerHits", &hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerHits, &b_hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerHits);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers", &hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers, &b_hltIter0IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_bestMatchTP_sharedFraction", &hltIter0IterL3MuonTrack_bestMatchTP_sharedFraction, &b_hltIter0IterL3MuonTrack_bestMatchTP_sharedFraction);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_matchedTPsize", &hltIter0IterL3MuonTrack_matchedTPsize, &b_hltIter0IterL3MuonTrack_matchedTPsize);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_mva0", &hltIter0IterL3MuonTrack_mva0, &b_hltIter0IterL3MuonTrack_mva0);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_mva1", &hltIter0IterL3MuonTrack_mva1, &b_hltIter0IterL3MuonTrack_mva1);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_mva2", &hltIter0IterL3MuonTrack_mva2, &b_hltIter0IterL3MuonTrack_mva2);
        fChain->SetBranchAddress("hltIter0IterL3MuonTrack_mva3", &hltIter0IterL3MuonTrack_mva3, &b_hltIter0IterL3MuonTrack_mva3);
        fChain->SetBranchAddress("nhltIter2IterL3MuonTrack", &nhltIter2IterL3MuonTrack, &b_nhltIter2IterL3MuonTrack);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_pt", &hltIter2IterL3MuonTrack_pt, &b_hltIter2IterL3MuonTrack_pt);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_ptError", &hltIter2IterL3MuonTrack_ptError, &b_hltIter2IterL3MuonTrack_ptError);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_eta", &hltIter2IterL3MuonTrack_eta, &b_hltIter2IterL3MuonTrack_eta);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_phi", &hltIter2IterL3MuonTrack_phi, &b_hltIter2IterL3MuonTrack_phi);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_charge", &hltIter2IterL3MuonTrack_charge, &b_hltIter2IterL3MuonTrack_charge);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_matchedL3", &hltIter2IterL3MuonTrack_matchedL3, &b_hltIter2IterL3MuonTrack_matchedL3);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_matchedL3NoId", &hltIter2IterL3MuonTrack_matchedL3NoId, &b_hltIter2IterL3MuonTrack_matchedL3NoId);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_charge", &hltIter2IterL3MuonTrack_bestMatchTP_charge, &b_hltIter2IterL3MuonTrack_bestMatchTP_charge);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_pdgId", &hltIter2IterL3MuonTrack_bestMatchTP_pdgId, &b_hltIter2IterL3MuonTrack_bestMatchTP_pdgId);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_energy", &hltIter2IterL3MuonTrack_bestMatchTP_energy, &b_hltIter2IterL3MuonTrack_bestMatchTP_energy);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_pt", &hltIter2IterL3MuonTrack_bestMatchTP_pt, &b_hltIter2IterL3MuonTrack_bestMatchTP_pt);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_eta", &hltIter2IterL3MuonTrack_bestMatchTP_eta, &b_hltIter2IterL3MuonTrack_bestMatchTP_eta);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_phi", &hltIter2IterL3MuonTrack_bestMatchTP_phi, &b_hltIter2IterL3MuonTrack_bestMatchTP_phi);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_parentVx", &hltIter2IterL3MuonTrack_bestMatchTP_parentVx, &b_hltIter2IterL3MuonTrack_bestMatchTP_parentVx);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_parentVy", &hltIter2IterL3MuonTrack_bestMatchTP_parentVy, &b_hltIter2IterL3MuonTrack_bestMatchTP_parentVy);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_parentVz", &hltIter2IterL3MuonTrack_bestMatchTP_parentVz, &b_hltIter2IterL3MuonTrack_bestMatchTP_parentVz);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_status", &hltIter2IterL3MuonTrack_bestMatchTP_status, &b_hltIter2IterL3MuonTrack_bestMatchTP_status);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_numberOfHits", &hltIter2IterL3MuonTrack_bestMatchTP_numberOfHits, &b_hltIter2IterL3MuonTrack_bestMatchTP_numberOfHits);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerHits", &hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerHits, &b_hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerHits);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers", &hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers, &b_hltIter2IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_bestMatchTP_sharedFraction", &hltIter2IterL3MuonTrack_bestMatchTP_sharedFraction, &b_hltIter2IterL3MuonTrack_bestMatchTP_sharedFraction);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_matchedTPsize", &hltIter2IterL3MuonTrack_matchedTPsize, &b_hltIter2IterL3MuonTrack_matchedTPsize);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_mva0", &hltIter2IterL3MuonTrack_mva0, &b_hltIter2IterL3MuonTrack_mva0);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_mva1", &hltIter2IterL3MuonTrack_mva1, &b_hltIter2IterL3MuonTrack_mva1);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_mva2", &hltIter2IterL3MuonTrack_mva2, &b_hltIter2IterL3MuonTrack_mva2);
        fChain->SetBranchAddress("hltIter2IterL3MuonTrack_mva3", &hltIter2IterL3MuonTrack_mva3, &b_hltIter2IterL3MuonTrack_mva3);
        fChain->SetBranchAddress("nhltIter3IterL3MuonTrack", &nhltIter3IterL3MuonTrack, &b_nhltIter3IterL3MuonTrack);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_pt", &hltIter3IterL3MuonTrack_pt, &b_hltIter3IterL3MuonTrack_pt);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_ptError", &hltIter3IterL3MuonTrack_ptError, &b_hltIter3IterL3MuonTrack_ptError);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_eta", &hltIter3IterL3MuonTrack_eta, &b_hltIter3IterL3MuonTrack_eta);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_phi", &hltIter3IterL3MuonTrack_phi, &b_hltIter3IterL3MuonTrack_phi);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_charge", &hltIter3IterL3MuonTrack_charge, &b_hltIter3IterL3MuonTrack_charge);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_matchedL3", &hltIter3IterL3MuonTrack_matchedL3, &b_hltIter3IterL3MuonTrack_matchedL3);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_matchedL3NoId", &hltIter3IterL3MuonTrack_matchedL3NoId, &b_hltIter3IterL3MuonTrack_matchedL3NoId);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_charge", &hltIter3IterL3MuonTrack_bestMatchTP_charge, &b_hltIter3IterL3MuonTrack_bestMatchTP_charge);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_pdgId", &hltIter3IterL3MuonTrack_bestMatchTP_pdgId, &b_hltIter3IterL3MuonTrack_bestMatchTP_pdgId);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_energy", &hltIter3IterL3MuonTrack_bestMatchTP_energy, &b_hltIter3IterL3MuonTrack_bestMatchTP_energy);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_pt", &hltIter3IterL3MuonTrack_bestMatchTP_pt, &b_hltIter3IterL3MuonTrack_bestMatchTP_pt);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_eta", &hltIter3IterL3MuonTrack_bestMatchTP_eta, &b_hltIter3IterL3MuonTrack_bestMatchTP_eta);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_phi", &hltIter3IterL3MuonTrack_bestMatchTP_phi, &b_hltIter3IterL3MuonTrack_bestMatchTP_phi);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_parentVx", &hltIter3IterL3MuonTrack_bestMatchTP_parentVx, &b_hltIter3IterL3MuonTrack_bestMatchTP_parentVx);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_parentVy", &hltIter3IterL3MuonTrack_bestMatchTP_parentVy, &b_hltIter3IterL3MuonTrack_bestMatchTP_parentVy);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_parentVz", &hltIter3IterL3MuonTrack_bestMatchTP_parentVz, &b_hltIter3IterL3MuonTrack_bestMatchTP_parentVz);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_status", &hltIter3IterL3MuonTrack_bestMatchTP_status, &b_hltIter3IterL3MuonTrack_bestMatchTP_status);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_numberOfHits", &hltIter3IterL3MuonTrack_bestMatchTP_numberOfHits, &b_hltIter3IterL3MuonTrack_bestMatchTP_numberOfHits);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerHits", &hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerHits, &b_hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerHits);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers", &hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers, &b_hltIter3IterL3MuonTrack_bestMatchTP_numberOfTrackerLayers);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_bestMatchTP_sharedFraction", &hltIter3IterL3MuonTrack_bestMatchTP_sharedFraction, &b_hltIter3IterL3MuonTrack_bestMatchTP_sharedFraction);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_matchedTPsize", &hltIter3IterL3MuonTrack_matchedTPsize, &b_hltIter3IterL3MuonTrack_matchedTPsize);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_mva0", &hltIter3IterL3MuonTrack_mva0, &b_hltIter3IterL3MuonTrack_mva0);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_mva1", &hltIter3IterL3MuonTrack_mva1, &b_hltIter3IterL3MuonTrack_mva1);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_mva2", &hltIter3IterL3MuonTrack_mva2, &b_hltIter3IterL3MuonTrack_mva2);
        fChain->SetBranchAddress("hltIter3IterL3MuonTrack_mva3", &hltIter3IterL3MuonTrack_mva3, &b_hltIter3IterL3MuonTrack_mva3);
        fChain->SetBranchAddress("nhltIter0IterL3FromL1MuonTrack", &nhltIter0IterL3FromL1MuonTrack, &b_nhltIter0IterL3FromL1MuonTrack);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_pt", &hltIter0IterL3FromL1MuonTrack_pt, &b_hltIter0IterL3FromL1MuonTrack_pt);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_ptError", &hltIter0IterL3FromL1MuonTrack_ptError, &b_hltIter0IterL3FromL1MuonTrack_ptError);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_eta", &hltIter0IterL3FromL1MuonTrack_eta, &b_hltIter0IterL3FromL1MuonTrack_eta);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_phi", &hltIter0IterL3FromL1MuonTrack_phi, &b_hltIter0IterL3FromL1MuonTrack_phi);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_charge", &hltIter0IterL3FromL1MuonTrack_charge, &b_hltIter0IterL3FromL1MuonTrack_charge);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_matchedL3", &hltIter0IterL3FromL1MuonTrack_matchedL3, &b_hltIter0IterL3FromL1MuonTrack_matchedL3);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_matchedL3NoId", &hltIter0IterL3FromL1MuonTrack_matchedL3NoId, &b_hltIter0IterL3FromL1MuonTrack_matchedL3NoId);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_charge", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_charge, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_charge);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_pdgId", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_pdgId, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_pdgId);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_energy", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_energy, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_energy);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_pt", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_pt, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_pt);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_eta", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_eta, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_eta);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_phi", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_phi, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_phi);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVx", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVx, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVx);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVy", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVy, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVy);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVz", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVz, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_parentVz);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_status", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_status, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_status);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfHits", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfHits, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfHits);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_bestMatchTP_sharedFraction", &hltIter0IterL3FromL1MuonTrack_bestMatchTP_sharedFraction, &b_hltIter0IterL3FromL1MuonTrack_bestMatchTP_sharedFraction);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_matchedTPsize", &hltIter0IterL3FromL1MuonTrack_matchedTPsize, &b_hltIter0IterL3FromL1MuonTrack_matchedTPsize);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_mva0", &hltIter0IterL3FromL1MuonTrack_mva0, &b_hltIter0IterL3FromL1MuonTrack_mva0);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_mva1", &hltIter0IterL3FromL1MuonTrack_mva1, &b_hltIter0IterL3FromL1MuonTrack_mva1);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_mva2", &hltIter0IterL3FromL1MuonTrack_mva2, &b_hltIter0IterL3FromL1MuonTrack_mva2);
        fChain->SetBranchAddress("hltIter0IterL3FromL1MuonTrack_mva3", &hltIter0IterL3FromL1MuonTrack_mva3, &b_hltIter0IterL3FromL1MuonTrack_mva3);
        fChain->SetBranchAddress("nhltIter2IterL3FromL1MuonTrack", &nhltIter2IterL3FromL1MuonTrack, &b_nhltIter2IterL3FromL1MuonTrack);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_pt", &hltIter2IterL3FromL1MuonTrack_pt, &b_hltIter2IterL3FromL1MuonTrack_pt);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_ptError", &hltIter2IterL3FromL1MuonTrack_ptError, &b_hltIter2IterL3FromL1MuonTrack_ptError);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_eta", &hltIter2IterL3FromL1MuonTrack_eta, &b_hltIter2IterL3FromL1MuonTrack_eta);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_phi", &hltIter2IterL3FromL1MuonTrack_phi, &b_hltIter2IterL3FromL1MuonTrack_phi);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_charge", &hltIter2IterL3FromL1MuonTrack_charge, &b_hltIter2IterL3FromL1MuonTrack_charge);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_matchedL3", &hltIter2IterL3FromL1MuonTrack_matchedL3, &b_hltIter2IterL3FromL1MuonTrack_matchedL3);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_matchedL3NoId", &hltIter2IterL3FromL1MuonTrack_matchedL3NoId, &b_hltIter2IterL3FromL1MuonTrack_matchedL3NoId);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_charge", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_charge, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_charge);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_pdgId", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_pdgId, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_pdgId);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_energy", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_energy, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_energy);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_pt", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_pt, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_pt);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_eta", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_eta, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_eta);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_phi", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_phi, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_phi);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVx", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVx, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVx);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVy", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVy, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVy);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVz", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVz, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_parentVz);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_status", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_status, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_status);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfHits", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfHits, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfHits);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_bestMatchTP_sharedFraction", &hltIter2IterL3FromL1MuonTrack_bestMatchTP_sharedFraction, &b_hltIter2IterL3FromL1MuonTrack_bestMatchTP_sharedFraction);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_matchedTPsize", &hltIter2IterL3FromL1MuonTrack_matchedTPsize, &b_hltIter2IterL3FromL1MuonTrack_matchedTPsize);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_mva0", &hltIter2IterL3FromL1MuonTrack_mva0, &b_hltIter2IterL3FromL1MuonTrack_mva0);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_mva1", &hltIter2IterL3FromL1MuonTrack_mva1, &b_hltIter2IterL3FromL1MuonTrack_mva1);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_mva2", &hltIter2IterL3FromL1MuonTrack_mva2, &b_hltIter2IterL3FromL1MuonTrack_mva2);
        fChain->SetBranchAddress("hltIter2IterL3FromL1MuonTrack_mva3", &hltIter2IterL3FromL1MuonTrack_mva3, &b_hltIter2IterL3FromL1MuonTrack_mva3);
        fChain->SetBranchAddress("nhltIter3IterL3FromL1MuonTrack", &nhltIter3IterL3FromL1MuonTrack, &b_nhltIter3IterL3FromL1MuonTrack);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_pt", &hltIter3IterL3FromL1MuonTrack_pt, &b_hltIter3IterL3FromL1MuonTrack_pt);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_ptError", &hltIter3IterL3FromL1MuonTrack_ptError, &b_hltIter3IterL3FromL1MuonTrack_ptError);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_eta", &hltIter3IterL3FromL1MuonTrack_eta, &b_hltIter3IterL3FromL1MuonTrack_eta);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_phi", &hltIter3IterL3FromL1MuonTrack_phi, &b_hltIter3IterL3FromL1MuonTrack_phi);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_charge", &hltIter3IterL3FromL1MuonTrack_charge, &b_hltIter3IterL3FromL1MuonTrack_charge);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_matchedL3", &hltIter3IterL3FromL1MuonTrack_matchedL3, &b_hltIter3IterL3FromL1MuonTrack_matchedL3);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_matchedL3NoId", &hltIter3IterL3FromL1MuonTrack_matchedL3NoId, &b_hltIter3IterL3FromL1MuonTrack_matchedL3NoId);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_charge", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_charge, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_charge);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_pdgId", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_pdgId, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_pdgId);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_energy", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_energy, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_energy);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_pt", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_pt, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_pt);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_eta", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_eta, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_eta);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_phi", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_phi, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_phi);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVx", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVx, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVx);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVy", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVy, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVy);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVz", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVz, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_parentVz);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_status", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_status, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_status);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfHits", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfHits, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfHits);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerHits);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_numberOfTrackerLayers);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_bestMatchTP_sharedFraction", &hltIter3IterL3FromL1MuonTrack_bestMatchTP_sharedFraction, &b_hltIter3IterL3FromL1MuonTrack_bestMatchTP_sharedFraction);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_matchedTPsize", &hltIter3IterL3FromL1MuonTrack_matchedTPsize, &b_hltIter3IterL3FromL1MuonTrack_matchedTPsize);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_mva0", &hltIter3IterL3FromL1MuonTrack_mva0, &b_hltIter3IterL3FromL1MuonTrack_mva0);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_mva1", &hltIter3IterL3FromL1MuonTrack_mva1, &b_hltIter3IterL3FromL1MuonTrack_mva1);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_mva2", &hltIter3IterL3FromL1MuonTrack_mva2, &b_hltIter3IterL3FromL1MuonTrack_mva2);
        fChain->SetBranchAddress("hltIter3IterL3FromL1MuonTrack_mva3", &hltIter3IterL3FromL1MuonTrack_mva3, &b_hltIter3IterL3FromL1MuonTrack_mva3);

        fChain->SetBranchAddress("nL3MuonsNoId", &nL3MuonsNoId, &b_nL3MuonsNoId);
        fChain->SetBranchAddress("L3MuonsNoId_pt", &L3MuonsNoId_pt, &b_L3MuonsNoId_pt);
        fChain->SetBranchAddress("L3MuonsNoId_inner_pt", &L3MuonsNoId_inner_pt, &b_L3MuonsNoId_inner_pt);
        fChain->SetBranchAddress("L3MuonsNoId_inner_ptError", &L3MuonsNoId_inner_ptError, &b_L3MuonsNoId_inner_ptError);
        fChain->SetBranchAddress("L3MuonsNoId_eta", &L3MuonsNoId_eta, &b_L3MuonsNoId_eta);
        fChain->SetBranchAddress("L3MuonsNoId_phi", &L3MuonsNoId_phi, &b_L3MuonsNoId_phi);
        fChain->SetBranchAddress("L3MuonsNoId_charge", &L3MuonsNoId_charge, &b_L3MuonsNoId_charge);
        fChain->SetBranchAddress("L3MuonsNoId_isGlobalMuon", &L3MuonsNoId_isGlobalMuon, &b_L3MuonsNoId_isGlobalMuon);
        fChain->SetBranchAddress("L3MuonsNoId_isStandAloneMuon", &L3MuonsNoId_isStandAloneMuon, &b_L3MuonsNoId_isStandAloneMuon);
        fChain->SetBranchAddress("L3MuonsNoId_isTrackerMuon", &L3MuonsNoId_isTrackerMuon, &b_L3MuonsNoId_isTrackerMuon);
        fChain->SetBranchAddress("L3MuonsNoId_isLooseTriggerMuon", &L3MuonsNoId_isLooseTriggerMuon, &b_L3MuonsNoId_isLooseTriggerMuon);
        fChain->SetBranchAddress("L3MuonsNoId_isME0Muon", &L3MuonsNoId_isME0Muon, &b_L3MuonsNoId_isME0Muon);
        fChain->SetBranchAddress("L3MuonsNoId_isGEMMuon", &L3MuonsNoId_isGEMMuon, &b_L3MuonsNoId_isGEMMuon);
        fChain->SetBranchAddress("L3MuonsNoId_isRPCMuon", &L3MuonsNoId_isRPCMuon, &b_L3MuonsNoId_isRPCMuon);
        fChain->SetBranchAddress("L3MuonsNoId_isGoodMuon_TMOneStationTight", &L3MuonsNoId_isGoodMuon_TMOneStationTight, &b_L3MuonsNoId_isGoodMuon_TMOneStationTight);
        fChain->SetBranchAddress("L3MuonsNoId_numberOfMatchedStations", &L3MuonsNoId_numberOfMatchedStations, &b_L3MuonsNoId_numberOfMatchedStations);
        fChain->SetBranchAddress("L3MuonsNoId_numberOfMatchedRPCLayers", &L3MuonsNoId_numberOfMatchedRPCLayers, &b_L3MuonsNoId_numberOfMatchedRPCLayers);
        fChain->SetBranchAddress("L3MuonsNoId_expectedNnumberOfMatchedStations", &L3MuonsNoId_expectedNnumberOfMatchedStations, &b_L3MuonsNoId_expectedNnumberOfMatchedStations);
        fChain->SetBranchAddress("L3MuonsNoId_inner_normalizedChi2", &L3MuonsNoId_inner_normalizedChi2, &b_L3MuonsNoId_inner_normalizedChi2);
        fChain->SetBranchAddress("L3MuonsNoId_inner_numberOfValidTrackerHits", &L3MuonsNoId_inner_numberOfValidTrackerHits, &b_L3MuonsNoId_inner_numberOfValidTrackerHits);
        fChain->SetBranchAddress("L3MuonsNoId_inner_trackerLayersWithMeasurement", &L3MuonsNoId_inner_trackerLayersWithMeasurement, &b_L3MuonsNoId_inner_trackerLayersWithMeasurement);
        fChain->SetBranchAddress("L3MuonsNoId_inner_numberOfValidPixelHits", &L3MuonsNoId_inner_numberOfValidPixelHits, &b_L3MuonsNoId_inner_numberOfValidPixelHits);
        fChain->SetBranchAddress("L3MuonsNoId_inner_dz_l1vtx", &L3MuonsNoId_inner_dz_l1vtx, &b_L3MuonsNoId_inner_dz_l1vtx);
        fChain->SetBranchAddress("nL3Muons", &nL3Muons, &b_nL3Muons);
        fChain->SetBranchAddress("L3Muons_pt", &L3Muons_pt, &b_L3Muons_pt);
        fChain->SetBranchAddress("L3Muons_inner_pt", &L3Muons_inner_pt, &b_L3Muons_inner_pt);
        fChain->SetBranchAddress("L3Muons_inner_ptError", &L3Muons_inner_ptError, &b_L3Muons_inner_ptError);
        fChain->SetBranchAddress("L3Muons_eta", &L3Muons_eta, &b_L3Muons_eta);
        fChain->SetBranchAddress("L3Muons_phi", &L3Muons_phi, &b_L3Muons_phi);
        fChain->SetBranchAddress("L3Muons_charge", &L3Muons_charge, &b_L3Muons_charge);
        fChain->SetBranchAddress("L3Muons_isGlobalMuon", &L3Muons_isGlobalMuon, &b_L3Muons_isGlobalMuon);
        fChain->SetBranchAddress("L3Muons_isStandAloneMuon", &L3Muons_isStandAloneMuon, &b_L3Muons_isStandAloneMuon);
        fChain->SetBranchAddress("L3Muons_isTrackerMuon", &L3Muons_isTrackerMuon, &b_L3Muons_isTrackerMuon);
        fChain->SetBranchAddress("L3Muons_isLooseTriggerMuon", &L3Muons_isLooseTriggerMuon, &b_L3Muons_isLooseTriggerMuon);
        fChain->SetBranchAddress("L3Muons_isME0Muon", &L3Muons_isME0Muon, &b_L3Muons_isME0Muon);
        fChain->SetBranchAddress("L3Muons_isGEMMuon", &L3Muons_isGEMMuon, &b_L3Muons_isGEMMuon);
        fChain->SetBranchAddress("L3Muons_isRPCMuon", &L3Muons_isRPCMuon, &b_L3Muons_isRPCMuon);
        fChain->SetBranchAddress("L3Muons_isGoodMuon_TMOneStationTight", &L3Muons_isGoodMuon_TMOneStationTight, &b_L3Muons_isGoodMuon_TMOneStationTight);
        fChain->SetBranchAddress("L3Muons_numberOfMatchedStations", &L3Muons_numberOfMatchedStations, &b_L3Muons_numberOfMatchedStations);
        fChain->SetBranchAddress("L3Muons_numberOfMatchedRPCLayers", &L3Muons_numberOfMatchedRPCLayers, &b_L3Muons_numberOfMatchedRPCLayers);
        fChain->SetBranchAddress("L3Muons_expectedNnumberOfMatchedStations", &L3Muons_expectedNnumberOfMatchedStations, &b_L3Muons_expectedNnumberOfMatchedStations);
        fChain->SetBranchAddress("L3Muons_inner_normalizedChi2", &L3Muons_inner_normalizedChi2, &b_L3Muons_inner_normalizedChi2);
        fChain->SetBranchAddress("L3Muons_inner_numberOfValidTrackerHits", &L3Muons_inner_numberOfValidTrackerHits, &b_L3Muons_inner_numberOfValidTrackerHits);
        fChain->SetBranchAddress("L3Muons_inner_trackerLayersWithMeasurement", &L3Muons_inner_trackerLayersWithMeasurement, &b_L3Muons_inner_trackerLayersWithMeasurement);
        fChain->SetBranchAddress("L3Muons_inner_numberOfValidPixelHits", &L3Muons_inner_numberOfValidPixelHits, &b_L3Muons_inner_numberOfValidPixelHits);
        fChain->SetBranchAddress("L3Muons_inner_dz_l1vtx", &L3Muons_inner_dz_l1vtx, &b_L3Muons_inner_dz_l1vtx);

        fChain->SetBranchAddress("nTP", &nTP, &b_nTP);
        fChain->SetBranchAddress("TP_charge", &TP_charge, &b_TP_charge);
        fChain->SetBranchAddress("TP_pdgId", &TP_pdgId, &b_TP_pdgId);
        fChain->SetBranchAddress("TP_energy", &TP_energy, &b_TP_energy);
        fChain->SetBranchAddress("TP_pt", &TP_pt, &b_TP_pt);
        fChain->SetBranchAddress("TP_eta", &TP_eta, &b_TP_eta);
        fChain->SetBranchAddress("TP_phi", &TP_phi, &b_TP_phi);
        fChain->SetBranchAddress("TP_parentVx", &TP_parentVx, &b_TP_parentVx);
        fChain->SetBranchAddress("TP_parentVy", &TP_parentVy, &b_TP_parentVy);
        fChain->SetBranchAddress("TP_parentVz", &TP_parentVz, &b_TP_parentVz);
        fChain->SetBranchAddress("TP_status", &TP_status, &b_TP_status);
        fChain->SetBranchAddress("TP_numberOfHits", &TP_numberOfHits, &b_TP_numberOfHits);
        fChain->SetBranchAddress("TP_numberOfTrackerHits", &TP_numberOfTrackerHits, &b_TP_numberOfTrackerHits);
        fChain->SetBranchAddress("TP_numberOfTrackerLayers", &TP_numberOfTrackerLayers, &b_TP_numberOfTrackerLayers);
        fChain->SetBranchAddress("TP_gen_charge", &TP_gen_charge, &b_TP_gen_charge);
        fChain->SetBranchAddress("TP_gen_pdgId", &TP_gen_pdgId, &b_TP_gen_pdgId);
        fChain->SetBranchAddress("TP_gen_pt", &TP_gen_pt, &b_TP_gen_pt);
        fChain->SetBranchAddress("TP_gen_eta", &TP_gen_eta, &b_TP_gen_eta);
        fChain->SetBranchAddress("TP_gen_phi", &TP_gen_phi, &b_TP_gen_phi);
        fChain->SetBranchAddress("TP_bestMatchTrk_pt", &TP_bestMatchTrk_pt, &b_TP_bestMatchTrk_pt);
        fChain->SetBranchAddress("TP_bestMatchTrk_eta", &TP_bestMatchTrk_eta, &b_TP_bestMatchTrk_eta);
        fChain->SetBranchAddress("TP_bestMatchTrk_phi", &TP_bestMatchTrk_phi, &b_TP_bestMatchTrk_phi);
        fChain->SetBranchAddress("TP_bestMatchTrk_charge", &TP_bestMatchTrk_charge, &b_TP_bestMatchTrk_charge);
        fChain->SetBranchAddress("TP_bestMatchTrk_quality", &TP_bestMatchTrk_quality, &b_TP_bestMatchTrk_quality);
        fChain->SetBranchAddress("TP_bestMatchTrk_NValidHits", &TP_bestMatchTrk_NValidHits, &b_TP_bestMatchTrk_NValidHits);
        fChain->SetBranchAddress("nhltIterL3MuonTrimmedPixelVertices", &nhltIterL3MuonTrimmedPixelVertices, &b_nhltIterL3MuonTrimmedPixelVertices);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_isValid", &hltIterL3MuonTrimmedPixelVertices_isValid, &b_hltIterL3MuonTrimmedPixelVertices_isValid);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_chi2", &hltIterL3MuonTrimmedPixelVertices_chi2, &b_hltIterL3MuonTrimmedPixelVertices_chi2);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_ndof", &hltIterL3MuonTrimmedPixelVertices_ndof, &b_hltIterL3MuonTrimmedPixelVertices_ndof);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_nTracks", &hltIterL3MuonTrimmedPixelVertices_nTracks, &b_hltIterL3MuonTrimmedPixelVertices_nTracks);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_x", &hltIterL3MuonTrimmedPixelVertices_x, &b_hltIterL3MuonTrimmedPixelVertices_x);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_xerr", &hltIterL3MuonTrimmedPixelVertices_xerr, &b_hltIterL3MuonTrimmedPixelVertices_xerr);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_y", &hltIterL3MuonTrimmedPixelVertices_y, &b_hltIterL3MuonTrimmedPixelVertices_y);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_yerr", &hltIterL3MuonTrimmedPixelVertices_yerr, &b_hltIterL3MuonTrimmedPixelVertices_yerr);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_z", &hltIterL3MuonTrimmedPixelVertices_z, &b_hltIterL3MuonTrimmedPixelVertices_z);
        fChain->SetBranchAddress("hltIterL3MuonTrimmedPixelVertices_zerr", &hltIterL3MuonTrimmedPixelVertices_zerr, &b_hltIterL3MuonTrimmedPixelVertices_zerr);
        fChain->SetBranchAddress("nhltIterL3FromL1MuonTrimmedPixelVertices", &nhltIterL3FromL1MuonTrimmedPixelVertices, &b_nhltIterL3FromL1MuonTrimmedPixelVertices);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_isValid", &hltIterL3FromL1MuonTrimmedPixelVertices_isValid, &b_hltIterL3FromL1MuonTrimmedPixelVertices_isValid);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_chi2", &hltIterL3FromL1MuonTrimmedPixelVertices_chi2, &b_hltIterL3FromL1MuonTrimmedPixelVertices_chi2);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_ndof", &hltIterL3FromL1MuonTrimmedPixelVertices_ndof, &b_hltIterL3FromL1MuonTrimmedPixelVertices_ndof);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_nTracks", &hltIterL3FromL1MuonTrimmedPixelVertices_nTracks, &b_hltIterL3FromL1MuonTrimmedPixelVertices_nTracks);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_x", &hltIterL3FromL1MuonTrimmedPixelVertices_x, &b_hltIterL3FromL1MuonTrimmedPixelVertices_x);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_xerr", &hltIterL3FromL1MuonTrimmedPixelVertices_xerr, &b_hltIterL3FromL1MuonTrimmedPixelVertices_xerr);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_y", &hltIterL3FromL1MuonTrimmedPixelVertices_y, &b_hltIterL3FromL1MuonTrimmedPixelVertices_y);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_yerr", &hltIterL3FromL1MuonTrimmedPixelVertices_yerr, &b_hltIterL3FromL1MuonTrimmedPixelVertices_yerr);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_z", &hltIterL3FromL1MuonTrimmedPixelVertices_z, &b_hltIterL3FromL1MuonTrimmedPixelVertices_z);
        fChain->SetBranchAddress("hltIterL3FromL1MuonTrimmedPixelVertices_zerr", &hltIterL3FromL1MuonTrimmedPixelVertices_zerr, &b_hltIterL3FromL1MuonTrimmedPixelVertices_zerr);

        if(!isIterL3) {
            fChain->SetBranchAddress("nhltPhase2L3OI", &nhltPhase2L3OI, &b_nhltPhase2L3OI);
            fChain->SetBranchAddress("hltPhase2L3OI_pt", &hltPhase2L3OI_pt, &b_hltPhase2L3OI_pt);
            fChain->SetBranchAddress("hltPhase2L3OI_ptError", &hltPhase2L3OI_ptError, &b_hltPhase2L3OI_ptError);
            fChain->SetBranchAddress("hltPhase2L3OI_eta", &hltPhase2L3OI_eta, &b_hltPhase2L3OI_eta);
            fChain->SetBranchAddress("hltPhase2L3OI_phi", &hltPhase2L3OI_phi, &b_hltPhase2L3OI_phi);
            fChain->SetBranchAddress("hltPhase2L3OI_charge", &hltPhase2L3OI_charge, &b_hltPhase2L3OI_charge);
            fChain->SetBranchAddress("hltPhase2L3OI_matchedL3", &hltPhase2L3OI_matchedL3, &b_hltPhase2L3OI_matchedL3);
            fChain->SetBranchAddress("hltPhase2L3OI_matchedL3NoId", &hltPhase2L3OI_matchedL3NoId, &b_hltPhase2L3OI_matchedL3NoId);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_charge", &hltPhase2L3OI_bestMatchTP_charge, &b_hltPhase2L3OI_bestMatchTP_charge);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_pdgId", &hltPhase2L3OI_bestMatchTP_pdgId, &b_hltPhase2L3OI_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_energy", &hltPhase2L3OI_bestMatchTP_energy, &b_hltPhase2L3OI_bestMatchTP_energy);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_pt", &hltPhase2L3OI_bestMatchTP_pt, &b_hltPhase2L3OI_bestMatchTP_pt);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_eta", &hltPhase2L3OI_bestMatchTP_eta, &b_hltPhase2L3OI_bestMatchTP_eta);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_phi", &hltPhase2L3OI_bestMatchTP_phi, &b_hltPhase2L3OI_bestMatchTP_phi);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_parentVx", &hltPhase2L3OI_bestMatchTP_parentVx, &b_hltPhase2L3OI_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_parentVy", &hltPhase2L3OI_bestMatchTP_parentVy, &b_hltPhase2L3OI_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_parentVz", &hltPhase2L3OI_bestMatchTP_parentVz, &b_hltPhase2L3OI_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_status", &hltPhase2L3OI_bestMatchTP_status, &b_hltPhase2L3OI_bestMatchTP_status);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_numberOfHits", &hltPhase2L3OI_bestMatchTP_numberOfHits, &b_hltPhase2L3OI_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_numberOfTrackerHits", &hltPhase2L3OI_bestMatchTP_numberOfTrackerHits, &b_hltPhase2L3OI_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_numberOfTrackerLayers", &hltPhase2L3OI_bestMatchTP_numberOfTrackerLayers, &b_hltPhase2L3OI_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltPhase2L3OI_bestMatchTP_sharedFraction", &hltPhase2L3OI_bestMatchTP_sharedFraction, &b_hltPhase2L3OI_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltPhase2L3OI_matchedTPsize", &hltPhase2L3OI_matchedTPsize, &b_hltPhase2L3OI_matchedTPsize);
            fChain->SetBranchAddress("hltPhase2L3OI_mva0", &hltPhase2L3OI_mva0, &b_hltPhase2L3OI_mva0);
            fChain->SetBranchAddress("hltPhase2L3OI_mva1", &hltPhase2L3OI_mva1, &b_hltPhase2L3OI_mva1);
            fChain->SetBranchAddress("hltPhase2L3OI_mva2", &hltPhase2L3OI_mva2, &b_hltPhase2L3OI_mva2);
            fChain->SetBranchAddress("hltPhase2L3OI_mva3", &hltPhase2L3OI_mva3, &b_hltPhase2L3OI_mva3);
            fChain->SetBranchAddress("ntpTo_hltPhase2L3OI", &ntpTo_hltPhase2L3OI, &b_ntpTo_hltPhase2L3OI);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_charge", &tpTo_hltPhase2L3OI_charge, &b_tpTo_hltPhase2L3OI_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_pdgId", &tpTo_hltPhase2L3OI_pdgId, &b_tpTo_hltPhase2L3OI_pdgId);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_energy", &tpTo_hltPhase2L3OI_energy, &b_tpTo_hltPhase2L3OI_energy);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_pt", &tpTo_hltPhase2L3OI_pt, &b_tpTo_hltPhase2L3OI_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_eta", &tpTo_hltPhase2L3OI_eta, &b_tpTo_hltPhase2L3OI_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_phi", &tpTo_hltPhase2L3OI_phi, &b_tpTo_hltPhase2L3OI_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_parentVx", &tpTo_hltPhase2L3OI_parentVx, &b_tpTo_hltPhase2L3OI_parentVx);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_parentVy", &tpTo_hltPhase2L3OI_parentVy, &b_tpTo_hltPhase2L3OI_parentVy);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_parentVz", &tpTo_hltPhase2L3OI_parentVz, &b_tpTo_hltPhase2L3OI_parentVz);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_status", &tpTo_hltPhase2L3OI_status, &b_tpTo_hltPhase2L3OI_status);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_numberOfHits", &tpTo_hltPhase2L3OI_numberOfHits, &b_tpTo_hltPhase2L3OI_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_numberOfTrackerHits", &tpTo_hltPhase2L3OI_numberOfTrackerHits, &b_tpTo_hltPhase2L3OI_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_numberOfTrackerLayers", &tpTo_hltPhase2L3OI_numberOfTrackerLayers, &b_tpTo_hltPhase2L3OI_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_gen_charge", &tpTo_hltPhase2L3OI_gen_charge, &b_tpTo_hltPhase2L3OI_gen_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_gen_pdgId", &tpTo_hltPhase2L3OI_gen_pdgId, &b_tpTo_hltPhase2L3OI_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_gen_pt", &tpTo_hltPhase2L3OI_gen_pt, &b_tpTo_hltPhase2L3OI_gen_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_gen_eta", &tpTo_hltPhase2L3OI_gen_eta, &b_tpTo_hltPhase2L3OI_gen_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_gen_phi", &tpTo_hltPhase2L3OI_gen_phi, &b_tpTo_hltPhase2L3OI_gen_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_bestMatchTrk_pt", &tpTo_hltPhase2L3OI_bestMatchTrk_pt, &b_tpTo_hltPhase2L3OI_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_bestMatchTrk_eta", &tpTo_hltPhase2L3OI_bestMatchTrk_eta, &b_tpTo_hltPhase2L3OI_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_bestMatchTrk_phi", &tpTo_hltPhase2L3OI_bestMatchTrk_phi, &b_tpTo_hltPhase2L3OI_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_bestMatchTrk_charge", &tpTo_hltPhase2L3OI_bestMatchTrk_charge, &b_tpTo_hltPhase2L3OI_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_bestMatchTrk_quality", &tpTo_hltPhase2L3OI_bestMatchTrk_quality, &b_tpTo_hltPhase2L3OI_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltPhase2L3OI_bestMatchTrk_NValidHits", &tpTo_hltPhase2L3OI_bestMatchTrk_NValidHits, &b_tpTo_hltPhase2L3OI_bestMatchTrk_NValidHits);

            fChain->SetBranchAddress("nhltIter0Phase2L3FromL1TkMuon", &nhltIter0Phase2L3FromL1TkMuon, &b_nhltIter0Phase2L3FromL1TkMuon);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_pt", &hltIter0Phase2L3FromL1TkMuon_pt, &b_hltIter0Phase2L3FromL1TkMuon_pt);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_ptError", &hltIter0Phase2L3FromL1TkMuon_ptError, &b_hltIter0Phase2L3FromL1TkMuon_ptError);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_eta", &hltIter0Phase2L3FromL1TkMuon_eta, &b_hltIter0Phase2L3FromL1TkMuon_eta);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_phi", &hltIter0Phase2L3FromL1TkMuon_phi, &b_hltIter0Phase2L3FromL1TkMuon_phi);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_charge", &hltIter0Phase2L3FromL1TkMuon_charge, &b_hltIter0Phase2L3FromL1TkMuon_charge);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_matchedL3", &hltIter0Phase2L3FromL1TkMuon_matchedL3, &b_hltIter0Phase2L3FromL1TkMuon_matchedL3);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_matchedL3NoId", &hltIter0Phase2L3FromL1TkMuon_matchedL3NoId, &b_hltIter0Phase2L3FromL1TkMuon_matchedL3NoId);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_charge", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_charge, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pdgId", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pdgId, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_energy", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_energy, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pt", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pt, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_eta", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_eta, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_phi", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_phi, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVx", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVx, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVy", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVy, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVz", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVz, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_status", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_status, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_status);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction", &hltIter0Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction, &b_hltIter0Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_matchedTPsize", &hltIter0Phase2L3FromL1TkMuon_matchedTPsize, &b_hltIter0Phase2L3FromL1TkMuon_matchedTPsize);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_mva0", &hltIter0Phase2L3FromL1TkMuon_mva0, &b_hltIter0Phase2L3FromL1TkMuon_mva0);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_mva1", &hltIter0Phase2L3FromL1TkMuon_mva1, &b_hltIter0Phase2L3FromL1TkMuon_mva1);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_mva2", &hltIter0Phase2L3FromL1TkMuon_mva2, &b_hltIter0Phase2L3FromL1TkMuon_mva2);
            fChain->SetBranchAddress("hltIter0Phase2L3FromL1TkMuon_mva3", &hltIter0Phase2L3FromL1TkMuon_mva3, &b_hltIter0Phase2L3FromL1TkMuon_mva3);
            fChain->SetBranchAddress("ntpTo_hltIter0Phase2L3FromL1TkMuon", &ntpTo_hltIter0Phase2L3FromL1TkMuon, &b_ntpTo_hltIter0Phase2L3FromL1TkMuon);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_charge", &tpTo_hltIter0Phase2L3FromL1TkMuon_charge, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_charge);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_pdgId", &tpTo_hltIter0Phase2L3FromL1TkMuon_pdgId, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_energy", &tpTo_hltIter0Phase2L3FromL1TkMuon_energy, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_energy);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_pt", &tpTo_hltIter0Phase2L3FromL1TkMuon_pt, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_pt);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_eta", &tpTo_hltIter0Phase2L3FromL1TkMuon_eta, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_eta);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_phi", &tpTo_hltIter0Phase2L3FromL1TkMuon_phi, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_phi);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_parentVx", &tpTo_hltIter0Phase2L3FromL1TkMuon_parentVx, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_parentVx);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_parentVy", &tpTo_hltIter0Phase2L3FromL1TkMuon_parentVy, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_parentVy);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_parentVz", &tpTo_hltIter0Phase2L3FromL1TkMuon_parentVz, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_parentVz);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_status", &tpTo_hltIter0Phase2L3FromL1TkMuon_status, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_status);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfHits", &tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfHits, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerHits", &tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerHits, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerLayers", &tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerLayers, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_gen_charge", &tpTo_hltIter0Phase2L3FromL1TkMuon_gen_charge, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pdgId", &tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pdgId, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pt", &tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pt, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_gen_eta", &tpTo_hltIter0Phase2L3FromL1TkMuon_gen_eta, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_gen_phi", &tpTo_hltIter0Phase2L3FromL1TkMuon_gen_phi, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_pt", &tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_pt, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_eta", &tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_eta, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_phi", &tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_phi, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_charge", &tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_charge, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_quality", &tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_quality, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits", &tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits, &b_tpTo_hltIter0Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltIter2Phase2L3FromL1TkMuon", &nhltIter2Phase2L3FromL1TkMuon, &b_nhltIter2Phase2L3FromL1TkMuon);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_pt", &hltIter2Phase2L3FromL1TkMuon_pt, &b_hltIter2Phase2L3FromL1TkMuon_pt);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_ptError", &hltIter2Phase2L3FromL1TkMuon_ptError, &b_hltIter2Phase2L3FromL1TkMuon_ptError);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_eta", &hltIter2Phase2L3FromL1TkMuon_eta, &b_hltIter2Phase2L3FromL1TkMuon_eta);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_phi", &hltIter2Phase2L3FromL1TkMuon_phi, &b_hltIter2Phase2L3FromL1TkMuon_phi);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_charge", &hltIter2Phase2L3FromL1TkMuon_charge, &b_hltIter2Phase2L3FromL1TkMuon_charge);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_matchedL3", &hltIter2Phase2L3FromL1TkMuon_matchedL3, &b_hltIter2Phase2L3FromL1TkMuon_matchedL3);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_matchedL3NoId", &hltIter2Phase2L3FromL1TkMuon_matchedL3NoId, &b_hltIter2Phase2L3FromL1TkMuon_matchedL3NoId);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_charge", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_charge, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pdgId", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pdgId, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_energy", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_energy, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pt", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pt, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_eta", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_eta, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_phi", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_phi, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVx", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVx, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVy", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVy, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVz", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVz, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_status", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_status, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_status);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction", &hltIter2Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction, &b_hltIter2Phase2L3FromL1TkMuon_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_matchedTPsize", &hltIter2Phase2L3FromL1TkMuon_matchedTPsize, &b_hltIter2Phase2L3FromL1TkMuon_matchedTPsize);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_mva0", &hltIter2Phase2L3FromL1TkMuon_mva0, &b_hltIter2Phase2L3FromL1TkMuon_mva0);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_mva1", &hltIter2Phase2L3FromL1TkMuon_mva1, &b_hltIter2Phase2L3FromL1TkMuon_mva1);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_mva2", &hltIter2Phase2L3FromL1TkMuon_mva2, &b_hltIter2Phase2L3FromL1TkMuon_mva2);
            fChain->SetBranchAddress("hltIter2Phase2L3FromL1TkMuon_mva3", &hltIter2Phase2L3FromL1TkMuon_mva3, &b_hltIter2Phase2L3FromL1TkMuon_mva3);
            fChain->SetBranchAddress("ntpTo_hltIter2Phase2L3FromL1TkMuon", &ntpTo_hltIter2Phase2L3FromL1TkMuon, &b_ntpTo_hltIter2Phase2L3FromL1TkMuon);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_charge", &tpTo_hltIter2Phase2L3FromL1TkMuon_charge, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_charge);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_pdgId", &tpTo_hltIter2Phase2L3FromL1TkMuon_pdgId, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_energy", &tpTo_hltIter2Phase2L3FromL1TkMuon_energy, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_energy);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_pt", &tpTo_hltIter2Phase2L3FromL1TkMuon_pt, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_pt);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_eta", &tpTo_hltIter2Phase2L3FromL1TkMuon_eta, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_eta);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_phi", &tpTo_hltIter2Phase2L3FromL1TkMuon_phi, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_phi);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_parentVx", &tpTo_hltIter2Phase2L3FromL1TkMuon_parentVx, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_parentVx);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_parentVy", &tpTo_hltIter2Phase2L3FromL1TkMuon_parentVy, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_parentVy);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_parentVz", &tpTo_hltIter2Phase2L3FromL1TkMuon_parentVz, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_parentVz);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_status", &tpTo_hltIter2Phase2L3FromL1TkMuon_status, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_status);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfHits", &tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfHits, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerHits", &tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerHits, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerLayers", &tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerLayers, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_gen_charge", &tpTo_hltIter2Phase2L3FromL1TkMuon_gen_charge, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pdgId", &tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pdgId, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pt", &tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pt, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_gen_eta", &tpTo_hltIter2Phase2L3FromL1TkMuon_gen_eta, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_gen_phi", &tpTo_hltIter2Phase2L3FromL1TkMuon_gen_phi, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_pt", &tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_pt, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_eta", &tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_eta, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_phi", &tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_phi, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_charge", &tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_charge, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_quality", &tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_quality, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits", &tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits, &b_tpTo_hltIter2Phase2L3FromL1TkMuon_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltPhase2L3IOFromL1", &nhltPhase2L3IOFromL1, &b_nhltPhase2L3IOFromL1);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_pt", &hltPhase2L3IOFromL1_pt, &b_hltPhase2L3IOFromL1_pt);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_ptError", &hltPhase2L3IOFromL1_ptError, &b_hltPhase2L3IOFromL1_ptError);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_eta", &hltPhase2L3IOFromL1_eta, &b_hltPhase2L3IOFromL1_eta);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_phi", &hltPhase2L3IOFromL1_phi, &b_hltPhase2L3IOFromL1_phi);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_charge", &hltPhase2L3IOFromL1_charge, &b_hltPhase2L3IOFromL1_charge);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_matchedL3", &hltPhase2L3IOFromL1_matchedL3, &b_hltPhase2L3IOFromL1_matchedL3);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_matchedL3NoId", &hltPhase2L3IOFromL1_matchedL3NoId, &b_hltPhase2L3IOFromL1_matchedL3NoId);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_charge", &hltPhase2L3IOFromL1_bestMatchTP_charge, &b_hltPhase2L3IOFromL1_bestMatchTP_charge);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_pdgId", &hltPhase2L3IOFromL1_bestMatchTP_pdgId, &b_hltPhase2L3IOFromL1_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_energy", &hltPhase2L3IOFromL1_bestMatchTP_energy, &b_hltPhase2L3IOFromL1_bestMatchTP_energy);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_pt", &hltPhase2L3IOFromL1_bestMatchTP_pt, &b_hltPhase2L3IOFromL1_bestMatchTP_pt);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_eta", &hltPhase2L3IOFromL1_bestMatchTP_eta, &b_hltPhase2L3IOFromL1_bestMatchTP_eta);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_phi", &hltPhase2L3IOFromL1_bestMatchTP_phi, &b_hltPhase2L3IOFromL1_bestMatchTP_phi);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_parentVx", &hltPhase2L3IOFromL1_bestMatchTP_parentVx, &b_hltPhase2L3IOFromL1_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_parentVy", &hltPhase2L3IOFromL1_bestMatchTP_parentVy, &b_hltPhase2L3IOFromL1_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_parentVz", &hltPhase2L3IOFromL1_bestMatchTP_parentVz, &b_hltPhase2L3IOFromL1_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_status", &hltPhase2L3IOFromL1_bestMatchTP_status, &b_hltPhase2L3IOFromL1_bestMatchTP_status);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_numberOfHits", &hltPhase2L3IOFromL1_bestMatchTP_numberOfHits, &b_hltPhase2L3IOFromL1_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerHits", &hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerHits, &b_hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerLayers", &hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerLayers, &b_hltPhase2L3IOFromL1_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_bestMatchTP_sharedFraction", &hltPhase2L3IOFromL1_bestMatchTP_sharedFraction, &b_hltPhase2L3IOFromL1_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_matchedTPsize", &hltPhase2L3IOFromL1_matchedTPsize, &b_hltPhase2L3IOFromL1_matchedTPsize);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_mva0", &hltPhase2L3IOFromL1_mva0, &b_hltPhase2L3IOFromL1_mva0);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_mva1", &hltPhase2L3IOFromL1_mva1, &b_hltPhase2L3IOFromL1_mva1);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_mva2", &hltPhase2L3IOFromL1_mva2, &b_hltPhase2L3IOFromL1_mva2);
            fChain->SetBranchAddress("hltPhase2L3IOFromL1_mva3", &hltPhase2L3IOFromL1_mva3, &b_hltPhase2L3IOFromL1_mva3);
            fChain->SetBranchAddress("ntpTo_hltPhase2L3IOFromL1", &ntpTo_hltPhase2L3IOFromL1, &b_ntpTo_hltPhase2L3IOFromL1);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_charge", &tpTo_hltPhase2L3IOFromL1_charge, &b_tpTo_hltPhase2L3IOFromL1_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_pdgId", &tpTo_hltPhase2L3IOFromL1_pdgId, &b_tpTo_hltPhase2L3IOFromL1_pdgId);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_energy", &tpTo_hltPhase2L3IOFromL1_energy, &b_tpTo_hltPhase2L3IOFromL1_energy);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_pt", &tpTo_hltPhase2L3IOFromL1_pt, &b_tpTo_hltPhase2L3IOFromL1_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_eta", &tpTo_hltPhase2L3IOFromL1_eta, &b_tpTo_hltPhase2L3IOFromL1_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_phi", &tpTo_hltPhase2L3IOFromL1_phi, &b_tpTo_hltPhase2L3IOFromL1_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_parentVx", &tpTo_hltPhase2L3IOFromL1_parentVx, &b_tpTo_hltPhase2L3IOFromL1_parentVx);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_parentVy", &tpTo_hltPhase2L3IOFromL1_parentVy, &b_tpTo_hltPhase2L3IOFromL1_parentVy);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_parentVz", &tpTo_hltPhase2L3IOFromL1_parentVz, &b_tpTo_hltPhase2L3IOFromL1_parentVz);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_status", &tpTo_hltPhase2L3IOFromL1_status, &b_tpTo_hltPhase2L3IOFromL1_status);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_numberOfHits", &tpTo_hltPhase2L3IOFromL1_numberOfHits, &b_tpTo_hltPhase2L3IOFromL1_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_numberOfTrackerHits", &tpTo_hltPhase2L3IOFromL1_numberOfTrackerHits, &b_tpTo_hltPhase2L3IOFromL1_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_numberOfTrackerLayers", &tpTo_hltPhase2L3IOFromL1_numberOfTrackerLayers, &b_tpTo_hltPhase2L3IOFromL1_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_gen_charge", &tpTo_hltPhase2L3IOFromL1_gen_charge, &b_tpTo_hltPhase2L3IOFromL1_gen_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_gen_pdgId", &tpTo_hltPhase2L3IOFromL1_gen_pdgId, &b_tpTo_hltPhase2L3IOFromL1_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_gen_pt", &tpTo_hltPhase2L3IOFromL1_gen_pt, &b_tpTo_hltPhase2L3IOFromL1_gen_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_gen_eta", &tpTo_hltPhase2L3IOFromL1_gen_eta, &b_tpTo_hltPhase2L3IOFromL1_gen_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_gen_phi", &tpTo_hltPhase2L3IOFromL1_gen_phi, &b_tpTo_hltPhase2L3IOFromL1_gen_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_bestMatchTrk_pt", &tpTo_hltPhase2L3IOFromL1_bestMatchTrk_pt, &b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_bestMatchTrk_eta", &tpTo_hltPhase2L3IOFromL1_bestMatchTrk_eta, &b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_bestMatchTrk_phi", &tpTo_hltPhase2L3IOFromL1_bestMatchTrk_phi, &b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_bestMatchTrk_charge", &tpTo_hltPhase2L3IOFromL1_bestMatchTrk_charge, &b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_bestMatchTrk_quality", &tpTo_hltPhase2L3IOFromL1_bestMatchTrk_quality, &b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltPhase2L3IOFromL1_bestMatchTrk_NValidHits", &tpTo_hltPhase2L3IOFromL1_bestMatchTrk_NValidHits, &b_tpTo_hltPhase2L3IOFromL1_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltPhase2L3MuonsNoID", &nhltPhase2L3MuonsNoID, &b_nhltPhase2L3MuonsNoID);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_pt", &hltPhase2L3MuonsNoID_pt, &b_hltPhase2L3MuonsNoID_pt);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_ptError", &hltPhase2L3MuonsNoID_ptError, &b_hltPhase2L3MuonsNoID_ptError);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_eta", &hltPhase2L3MuonsNoID_eta, &b_hltPhase2L3MuonsNoID_eta);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_phi", &hltPhase2L3MuonsNoID_phi, &b_hltPhase2L3MuonsNoID_phi);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_charge", &hltPhase2L3MuonsNoID_charge, &b_hltPhase2L3MuonsNoID_charge);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_matchedL3", &hltPhase2L3MuonsNoID_matchedL3, &b_hltPhase2L3MuonsNoID_matchedL3);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_matchedL3NoId", &hltPhase2L3MuonsNoID_matchedL3NoId, &b_hltPhase2L3MuonsNoID_matchedL3NoId);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_charge", &hltPhase2L3MuonsNoID_bestMatchTP_charge, &b_hltPhase2L3MuonsNoID_bestMatchTP_charge);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_pdgId", &hltPhase2L3MuonsNoID_bestMatchTP_pdgId, &b_hltPhase2L3MuonsNoID_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_energy", &hltPhase2L3MuonsNoID_bestMatchTP_energy, &b_hltPhase2L3MuonsNoID_bestMatchTP_energy);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_pt", &hltPhase2L3MuonsNoID_bestMatchTP_pt, &b_hltPhase2L3MuonsNoID_bestMatchTP_pt);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_eta", &hltPhase2L3MuonsNoID_bestMatchTP_eta, &b_hltPhase2L3MuonsNoID_bestMatchTP_eta);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_phi", &hltPhase2L3MuonsNoID_bestMatchTP_phi, &b_hltPhase2L3MuonsNoID_bestMatchTP_phi);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_parentVx", &hltPhase2L3MuonsNoID_bestMatchTP_parentVx, &b_hltPhase2L3MuonsNoID_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_parentVy", &hltPhase2L3MuonsNoID_bestMatchTP_parentVy, &b_hltPhase2L3MuonsNoID_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_parentVz", &hltPhase2L3MuonsNoID_bestMatchTP_parentVz, &b_hltPhase2L3MuonsNoID_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_status", &hltPhase2L3MuonsNoID_bestMatchTP_status, &b_hltPhase2L3MuonsNoID_bestMatchTP_status);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_numberOfHits", &hltPhase2L3MuonsNoID_bestMatchTP_numberOfHits, &b_hltPhase2L3MuonsNoID_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerHits", &hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerHits, &b_hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerLayers", &hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerLayers, &b_hltPhase2L3MuonsNoID_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_bestMatchTP_sharedFraction", &hltPhase2L3MuonsNoID_bestMatchTP_sharedFraction, &b_hltPhase2L3MuonsNoID_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_matchedTPsize", &hltPhase2L3MuonsNoID_matchedTPsize, &b_hltPhase2L3MuonsNoID_matchedTPsize);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_mva0", &hltPhase2L3MuonsNoID_mva0, &b_hltPhase2L3MuonsNoID_mva0);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_mva1", &hltPhase2L3MuonsNoID_mva1, &b_hltPhase2L3MuonsNoID_mva1);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_mva2", &hltPhase2L3MuonsNoID_mva2, &b_hltPhase2L3MuonsNoID_mva2);
            fChain->SetBranchAddress("hltPhase2L3MuonsNoID_mva3", &hltPhase2L3MuonsNoID_mva3, &b_hltPhase2L3MuonsNoID_mva3);
            fChain->SetBranchAddress("ntpTo_hltPhase2L3MuonsNoID", &ntpTo_hltPhase2L3MuonsNoID, &b_ntpTo_hltPhase2L3MuonsNoID);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_charge", &tpTo_hltPhase2L3MuonsNoID_charge, &b_tpTo_hltPhase2L3MuonsNoID_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_pdgId", &tpTo_hltPhase2L3MuonsNoID_pdgId, &b_tpTo_hltPhase2L3MuonsNoID_pdgId);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_energy", &tpTo_hltPhase2L3MuonsNoID_energy, &b_tpTo_hltPhase2L3MuonsNoID_energy);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_pt", &tpTo_hltPhase2L3MuonsNoID_pt, &b_tpTo_hltPhase2L3MuonsNoID_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_eta", &tpTo_hltPhase2L3MuonsNoID_eta, &b_tpTo_hltPhase2L3MuonsNoID_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_phi", &tpTo_hltPhase2L3MuonsNoID_phi, &b_tpTo_hltPhase2L3MuonsNoID_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_parentVx", &tpTo_hltPhase2L3MuonsNoID_parentVx, &b_tpTo_hltPhase2L3MuonsNoID_parentVx);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_parentVy", &tpTo_hltPhase2L3MuonsNoID_parentVy, &b_tpTo_hltPhase2L3MuonsNoID_parentVy);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_parentVz", &tpTo_hltPhase2L3MuonsNoID_parentVz, &b_tpTo_hltPhase2L3MuonsNoID_parentVz);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_status", &tpTo_hltPhase2L3MuonsNoID_status, &b_tpTo_hltPhase2L3MuonsNoID_status);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_numberOfHits", &tpTo_hltPhase2L3MuonsNoID_numberOfHits, &b_tpTo_hltPhase2L3MuonsNoID_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_numberOfTrackerHits", &tpTo_hltPhase2L3MuonsNoID_numberOfTrackerHits, &b_tpTo_hltPhase2L3MuonsNoID_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_numberOfTrackerLayers", &tpTo_hltPhase2L3MuonsNoID_numberOfTrackerLayers, &b_tpTo_hltPhase2L3MuonsNoID_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_gen_charge", &tpTo_hltPhase2L3MuonsNoID_gen_charge, &b_tpTo_hltPhase2L3MuonsNoID_gen_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_gen_pdgId", &tpTo_hltPhase2L3MuonsNoID_gen_pdgId, &b_tpTo_hltPhase2L3MuonsNoID_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_gen_pt", &tpTo_hltPhase2L3MuonsNoID_gen_pt, &b_tpTo_hltPhase2L3MuonsNoID_gen_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_gen_eta", &tpTo_hltPhase2L3MuonsNoID_gen_eta, &b_tpTo_hltPhase2L3MuonsNoID_gen_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_gen_phi", &tpTo_hltPhase2L3MuonsNoID_gen_phi, &b_tpTo_hltPhase2L3MuonsNoID_gen_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_pt", &tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_pt, &b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_eta", &tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_eta, &b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_phi", &tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_phi, &b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_charge", &tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_charge, &b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_quality", &tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_quality, &b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_NValidHits", &tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_NValidHits, &b_tpTo_hltPhase2L3MuonsNoID_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltPhase2L3Muons", &nhltPhase2L3Muons, &b_nhltPhase2L3Muons);
            fChain->SetBranchAddress("hltPhase2L3Muons_pt", &hltPhase2L3Muons_pt, &b_hltPhase2L3Muons_pt);
            fChain->SetBranchAddress("hltPhase2L3Muons_ptError", &hltPhase2L3Muons_ptError, &b_hltPhase2L3Muons_ptError);
            fChain->SetBranchAddress("hltPhase2L3Muons_eta", &hltPhase2L3Muons_eta, &b_hltPhase2L3Muons_eta);
            fChain->SetBranchAddress("hltPhase2L3Muons_phi", &hltPhase2L3Muons_phi, &b_hltPhase2L3Muons_phi);
            fChain->SetBranchAddress("hltPhase2L3Muons_charge", &hltPhase2L3Muons_charge, &b_hltPhase2L3Muons_charge);
            fChain->SetBranchAddress("hltPhase2L3Muons_matchedL3", &hltPhase2L3Muons_matchedL3, &b_hltPhase2L3Muons_matchedL3);
            fChain->SetBranchAddress("hltPhase2L3Muons_matchedL3NoId", &hltPhase2L3Muons_matchedL3NoId, &b_hltPhase2L3Muons_matchedL3NoId);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_charge", &hltPhase2L3Muons_bestMatchTP_charge, &b_hltPhase2L3Muons_bestMatchTP_charge);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_pdgId", &hltPhase2L3Muons_bestMatchTP_pdgId, &b_hltPhase2L3Muons_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_energy", &hltPhase2L3Muons_bestMatchTP_energy, &b_hltPhase2L3Muons_bestMatchTP_energy);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_pt", &hltPhase2L3Muons_bestMatchTP_pt, &b_hltPhase2L3Muons_bestMatchTP_pt);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_eta", &hltPhase2L3Muons_bestMatchTP_eta, &b_hltPhase2L3Muons_bestMatchTP_eta);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_phi", &hltPhase2L3Muons_bestMatchTP_phi, &b_hltPhase2L3Muons_bestMatchTP_phi);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_parentVx", &hltPhase2L3Muons_bestMatchTP_parentVx, &b_hltPhase2L3Muons_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_parentVy", &hltPhase2L3Muons_bestMatchTP_parentVy, &b_hltPhase2L3Muons_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_parentVz", &hltPhase2L3Muons_bestMatchTP_parentVz, &b_hltPhase2L3Muons_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_status", &hltPhase2L3Muons_bestMatchTP_status, &b_hltPhase2L3Muons_bestMatchTP_status);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_numberOfHits", &hltPhase2L3Muons_bestMatchTP_numberOfHits, &b_hltPhase2L3Muons_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_numberOfTrackerHits", &hltPhase2L3Muons_bestMatchTP_numberOfTrackerHits, &b_hltPhase2L3Muons_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_numberOfTrackerLayers", &hltPhase2L3Muons_bestMatchTP_numberOfTrackerLayers, &b_hltPhase2L3Muons_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltPhase2L3Muons_bestMatchTP_sharedFraction", &hltPhase2L3Muons_bestMatchTP_sharedFraction, &b_hltPhase2L3Muons_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltPhase2L3Muons_matchedTPsize", &hltPhase2L3Muons_matchedTPsize, &b_hltPhase2L3Muons_matchedTPsize);
            fChain->SetBranchAddress("hltPhase2L3Muons_mva0", &hltPhase2L3Muons_mva0, &b_hltPhase2L3Muons_mva0);
            fChain->SetBranchAddress("hltPhase2L3Muons_mva1", &hltPhase2L3Muons_mva1, &b_hltPhase2L3Muons_mva1);
            fChain->SetBranchAddress("hltPhase2L3Muons_mva2", &hltPhase2L3Muons_mva2, &b_hltPhase2L3Muons_mva2);
            fChain->SetBranchAddress("hltPhase2L3Muons_mva3", &hltPhase2L3Muons_mva3, &b_hltPhase2L3Muons_mva3);

            // -- Isolation
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p10dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p20dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p10ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionaldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoFulldR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0", &hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0, &b_hltPhase2L3Muons_trkIsoOfflinedR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p000", &hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p000, &b_hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p000);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p000", &hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p000, &b_hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p000);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p030", &hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p030, &b_hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p030);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p030", &hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p030, &b_hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p030);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p050", &hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p050, &b_hltPhase2L3Muons_pfEcalIsodR0p3dRVeto0p050);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p050", &hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p050, &b_hltPhase2L3Muons_pfHcalIsodR0p3dRVeto0p050);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p00minEEM0p00minEHad0p00", &hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p00minEEM0p00minEHad0p00, &b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p00minEEM0p00minEHad0p00);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00", &hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00, &b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p04minEEM0p00minEHad0p00", &hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p04minEEM0p00minEHad0p00, &b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p04minEEM0p00minEHad0p00);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p00minEEM0p00minEHad0p00", &hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p00minEEM0p00minEHad0p00, &b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p00minEEM0p00minEHad0p00);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p02minEEM0p00minEHad0p00", &hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p02minEEM0p00minEHad0p00, &b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p02minEEM0p00minEHad0p00);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p04minEEM0p00minEHad0p00", &hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p04minEEM0p00minEHad0p00, &b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p02dRVetoHad0p04minEEM0p00minEHad0p00);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p00minEEM0p00minEHad0p00", &hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p00minEEM0p00minEHad0p00, &b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p00minEEM0p00minEHad0p00);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p02minEEM0p00minEHad0p00", &hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p02minEEM0p00minEHad0p00, &b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p02minEEM0p00minEHad0p00);
            //fChain->SetBranchAddress("hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p04minEEM0p00minEHad0p00", &hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p04minEEM0p00minEHad0p00, &b_hltPhase2L3Muons_pfHgcalLCIsodR0p2dRVetoEM0p04dRVetoHad0p04minEEM0p00minEHad0p00);

            fChain->SetBranchAddress("ntpTo_hltPhase2L3Muons", &ntpTo_hltPhase2L3Muons, &b_ntpTo_hltPhase2L3Muons);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_charge", &tpTo_hltPhase2L3Muons_charge, &b_tpTo_hltPhase2L3Muons_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_pdgId", &tpTo_hltPhase2L3Muons_pdgId, &b_tpTo_hltPhase2L3Muons_pdgId);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_energy", &tpTo_hltPhase2L3Muons_energy, &b_tpTo_hltPhase2L3Muons_energy);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_pt", &tpTo_hltPhase2L3Muons_pt, &b_tpTo_hltPhase2L3Muons_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_eta", &tpTo_hltPhase2L3Muons_eta, &b_tpTo_hltPhase2L3Muons_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_phi", &tpTo_hltPhase2L3Muons_phi, &b_tpTo_hltPhase2L3Muons_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_parentVx", &tpTo_hltPhase2L3Muons_parentVx, &b_tpTo_hltPhase2L3Muons_parentVx);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_parentVy", &tpTo_hltPhase2L3Muons_parentVy, &b_tpTo_hltPhase2L3Muons_parentVy);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_parentVz", &tpTo_hltPhase2L3Muons_parentVz, &b_tpTo_hltPhase2L3Muons_parentVz);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_status", &tpTo_hltPhase2L3Muons_status, &b_tpTo_hltPhase2L3Muons_status);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_numberOfHits", &tpTo_hltPhase2L3Muons_numberOfHits, &b_tpTo_hltPhase2L3Muons_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_numberOfTrackerHits", &tpTo_hltPhase2L3Muons_numberOfTrackerHits, &b_tpTo_hltPhase2L3Muons_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_numberOfTrackerLayers", &tpTo_hltPhase2L3Muons_numberOfTrackerLayers, &b_tpTo_hltPhase2L3Muons_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_gen_charge", &tpTo_hltPhase2L3Muons_gen_charge, &b_tpTo_hltPhase2L3Muons_gen_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_gen_pdgId", &tpTo_hltPhase2L3Muons_gen_pdgId, &b_tpTo_hltPhase2L3Muons_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_gen_pt", &tpTo_hltPhase2L3Muons_gen_pt, &b_tpTo_hltPhase2L3Muons_gen_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_gen_eta", &tpTo_hltPhase2L3Muons_gen_eta, &b_tpTo_hltPhase2L3Muons_gen_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_gen_phi", &tpTo_hltPhase2L3Muons_gen_phi, &b_tpTo_hltPhase2L3Muons_gen_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_bestMatchTrk_pt", &tpTo_hltPhase2L3Muons_bestMatchTrk_pt, &b_tpTo_hltPhase2L3Muons_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_bestMatchTrk_eta", &tpTo_hltPhase2L3Muons_bestMatchTrk_eta, &b_tpTo_hltPhase2L3Muons_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_bestMatchTrk_phi", &tpTo_hltPhase2L3Muons_bestMatchTrk_phi, &b_tpTo_hltPhase2L3Muons_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_bestMatchTrk_charge", &tpTo_hltPhase2L3Muons_bestMatchTrk_charge, &b_tpTo_hltPhase2L3Muons_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_bestMatchTrk_quality", &tpTo_hltPhase2L3Muons_bestMatchTrk_quality, &b_tpTo_hltPhase2L3Muons_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltPhase2L3Muons_bestMatchTrk_NValidHits", &tpTo_hltPhase2L3Muons_bestMatchTrk_NValidHits, &b_tpTo_hltPhase2L3Muons_bestMatchTrk_NValidHits);
        }

        else {
            fChain->SetBranchAddress("nhltIterL3OI", &nhltIterL3OI, &b_nhltIterL3OI);
            fChain->SetBranchAddress("hltIterL3OI_pt", &hltIterL3OI_pt, &b_hltIterL3OI_pt);
            fChain->SetBranchAddress("hltIterL3OI_ptError", &hltIterL3OI_ptError, &b_hltIterL3OI_ptError);
            fChain->SetBranchAddress("hltIterL3OI_eta", &hltIterL3OI_eta, &b_hltIterL3OI_eta);
            fChain->SetBranchAddress("hltIterL3OI_phi", &hltIterL3OI_phi, &b_hltIterL3OI_phi);
            fChain->SetBranchAddress("hltIterL3OI_charge", &hltIterL3OI_charge, &b_hltIterL3OI_charge);
            fChain->SetBranchAddress("hltIterL3OI_matchedL3", &hltIterL3OI_matchedL3, &b_hltIterL3OI_matchedL3);
            fChain->SetBranchAddress("hltIterL3OI_matchedL3NoId", &hltIterL3OI_matchedL3NoId, &b_hltIterL3OI_matchedL3NoId);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_charge", &hltIterL3OI_bestMatchTP_charge, &b_hltIterL3OI_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_pdgId", &hltIterL3OI_bestMatchTP_pdgId, &b_hltIterL3OI_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_energy", &hltIterL3OI_bestMatchTP_energy, &b_hltIterL3OI_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_pt", &hltIterL3OI_bestMatchTP_pt, &b_hltIterL3OI_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_eta", &hltIterL3OI_bestMatchTP_eta, &b_hltIterL3OI_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_phi", &hltIterL3OI_bestMatchTP_phi, &b_hltIterL3OI_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_parentVx", &hltIterL3OI_bestMatchTP_parentVx, &b_hltIterL3OI_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_parentVy", &hltIterL3OI_bestMatchTP_parentVy, &b_hltIterL3OI_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_parentVz", &hltIterL3OI_bestMatchTP_parentVz, &b_hltIterL3OI_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_status", &hltIterL3OI_bestMatchTP_status, &b_hltIterL3OI_bestMatchTP_status);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_numberOfHits", &hltIterL3OI_bestMatchTP_numberOfHits, &b_hltIterL3OI_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_numberOfTrackerHits", &hltIterL3OI_bestMatchTP_numberOfTrackerHits, &b_hltIterL3OI_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_numberOfTrackerLayers", &hltIterL3OI_bestMatchTP_numberOfTrackerLayers, &b_hltIterL3OI_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIterL3OI_bestMatchTP_sharedFraction", &hltIterL3OI_bestMatchTP_sharedFraction, &b_hltIterL3OI_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIterL3OI_matchedTPsize", &hltIterL3OI_matchedTPsize, &b_hltIterL3OI_matchedTPsize);
            fChain->SetBranchAddress("hltIterL3OI_mva0", &hltIterL3OI_mva0, &b_hltIterL3OI_mva0);
            fChain->SetBranchAddress("hltIterL3OI_mva1", &hltIterL3OI_mva1, &b_hltIterL3OI_mva1);
            fChain->SetBranchAddress("hltIterL3OI_mva2", &hltIterL3OI_mva2, &b_hltIterL3OI_mva2);
            fChain->SetBranchAddress("hltIterL3OI_mva3", &hltIterL3OI_mva3, &b_hltIterL3OI_mva3);
            fChain->SetBranchAddress("ntpTo_hltIterL3OI", &ntpTo_hltIterL3OI, &b_ntpTo_hltIterL3OI);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_charge", &tpTo_hltIterL3OI_charge, &b_tpTo_hltIterL3OI_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_pdgId", &tpTo_hltIterL3OI_pdgId, &b_tpTo_hltIterL3OI_pdgId);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_energy", &tpTo_hltIterL3OI_energy, &b_tpTo_hltIterL3OI_energy);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_pt", &tpTo_hltIterL3OI_pt, &b_tpTo_hltIterL3OI_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_eta", &tpTo_hltIterL3OI_eta, &b_tpTo_hltIterL3OI_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_phi", &tpTo_hltIterL3OI_phi, &b_tpTo_hltIterL3OI_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_parentVx", &tpTo_hltIterL3OI_parentVx, &b_tpTo_hltIterL3OI_parentVx);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_parentVy", &tpTo_hltIterL3OI_parentVy, &b_tpTo_hltIterL3OI_parentVy);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_parentVz", &tpTo_hltIterL3OI_parentVz, &b_tpTo_hltIterL3OI_parentVz);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_status", &tpTo_hltIterL3OI_status, &b_tpTo_hltIterL3OI_status);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_numberOfHits", &tpTo_hltIterL3OI_numberOfHits, &b_tpTo_hltIterL3OI_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_numberOfTrackerHits", &tpTo_hltIterL3OI_numberOfTrackerHits, &b_tpTo_hltIterL3OI_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_numberOfTrackerLayers", &tpTo_hltIterL3OI_numberOfTrackerLayers, &b_tpTo_hltIterL3OI_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_gen_charge", &tpTo_hltIterL3OI_gen_charge, &b_tpTo_hltIterL3OI_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_gen_pdgId", &tpTo_hltIterL3OI_gen_pdgId, &b_tpTo_hltIterL3OI_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_gen_pt", &tpTo_hltIterL3OI_gen_pt, &b_tpTo_hltIterL3OI_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_gen_eta", &tpTo_hltIterL3OI_gen_eta, &b_tpTo_hltIterL3OI_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_gen_phi", &tpTo_hltIterL3OI_gen_phi, &b_tpTo_hltIterL3OI_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_bestMatchTrk_pt", &tpTo_hltIterL3OI_bestMatchTrk_pt, &b_tpTo_hltIterL3OI_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_bestMatchTrk_eta", &tpTo_hltIterL3OI_bestMatchTrk_eta, &b_tpTo_hltIterL3OI_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_bestMatchTrk_phi", &tpTo_hltIterL3OI_bestMatchTrk_phi, &b_tpTo_hltIterL3OI_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_bestMatchTrk_charge", &tpTo_hltIterL3OI_bestMatchTrk_charge, &b_tpTo_hltIterL3OI_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_bestMatchTrk_quality", &tpTo_hltIterL3OI_bestMatchTrk_quality, &b_tpTo_hltIterL3OI_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIterL3OI_bestMatchTrk_NValidHits", &tpTo_hltIterL3OI_bestMatchTrk_NValidHits, &b_tpTo_hltIterL3OI_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltIter0IterL3", &nhltIter0IterL3, &b_nhltIter0IterL3);
            fChain->SetBranchAddress("hltIter0IterL3_pt", &hltIter0IterL3_pt, &b_hltIter0IterL3_pt);
            fChain->SetBranchAddress("hltIter0IterL3_ptError", &hltIter0IterL3_ptError, &b_hltIter0IterL3_ptError);
            fChain->SetBranchAddress("hltIter0IterL3_eta", &hltIter0IterL3_eta, &b_hltIter0IterL3_eta);
            fChain->SetBranchAddress("hltIter0IterL3_phi", &hltIter0IterL3_phi, &b_hltIter0IterL3_phi);
            fChain->SetBranchAddress("hltIter0IterL3_charge", &hltIter0IterL3_charge, &b_hltIter0IterL3_charge);
            fChain->SetBranchAddress("hltIter0IterL3_matchedL3", &hltIter0IterL3_matchedL3, &b_hltIter0IterL3_matchedL3);
            fChain->SetBranchAddress("hltIter0IterL3_matchedL3NoId", &hltIter0IterL3_matchedL3NoId, &b_hltIter0IterL3_matchedL3NoId);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_charge", &hltIter0IterL3_bestMatchTP_charge, &b_hltIter0IterL3_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_pdgId", &hltIter0IterL3_bestMatchTP_pdgId, &b_hltIter0IterL3_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_energy", &hltIter0IterL3_bestMatchTP_energy, &b_hltIter0IterL3_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_pt", &hltIter0IterL3_bestMatchTP_pt, &b_hltIter0IterL3_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_eta", &hltIter0IterL3_bestMatchTP_eta, &b_hltIter0IterL3_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_phi", &hltIter0IterL3_bestMatchTP_phi, &b_hltIter0IterL3_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_parentVx", &hltIter0IterL3_bestMatchTP_parentVx, &b_hltIter0IterL3_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_parentVy", &hltIter0IterL3_bestMatchTP_parentVy, &b_hltIter0IterL3_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_parentVz", &hltIter0IterL3_bestMatchTP_parentVz, &b_hltIter0IterL3_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_status", &hltIter0IterL3_bestMatchTP_status, &b_hltIter0IterL3_bestMatchTP_status);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_numberOfHits", &hltIter0IterL3_bestMatchTP_numberOfHits, &b_hltIter0IterL3_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_numberOfTrackerHits", &hltIter0IterL3_bestMatchTP_numberOfTrackerHits, &b_hltIter0IterL3_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_numberOfTrackerLayers", &hltIter0IterL3_bestMatchTP_numberOfTrackerLayers, &b_hltIter0IterL3_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIter0IterL3_bestMatchTP_sharedFraction", &hltIter0IterL3_bestMatchTP_sharedFraction, &b_hltIter0IterL3_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIter0IterL3_matchedTPsize", &hltIter0IterL3_matchedTPsize, &b_hltIter0IterL3_matchedTPsize);
            fChain->SetBranchAddress("hltIter0IterL3_mva0", &hltIter0IterL3_mva0, &b_hltIter0IterL3_mva0);
            fChain->SetBranchAddress("hltIter0IterL3_mva1", &hltIter0IterL3_mva1, &b_hltIter0IterL3_mva1);
            fChain->SetBranchAddress("hltIter0IterL3_mva2", &hltIter0IterL3_mva2, &b_hltIter0IterL3_mva2);
            fChain->SetBranchAddress("hltIter0IterL3_mva3", &hltIter0IterL3_mva3, &b_hltIter0IterL3_mva3);
            fChain->SetBranchAddress("ntpTo_hltIter0IterL3", &ntpTo_hltIter0IterL3, &b_ntpTo_hltIter0IterL3);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_charge", &tpTo_hltIter0IterL3_charge, &b_tpTo_hltIter0IterL3_charge);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_pdgId", &tpTo_hltIter0IterL3_pdgId, &b_tpTo_hltIter0IterL3_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_energy", &tpTo_hltIter0IterL3_energy, &b_tpTo_hltIter0IterL3_energy);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_pt", &tpTo_hltIter0IterL3_pt, &b_tpTo_hltIter0IterL3_pt);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_eta", &tpTo_hltIter0IterL3_eta, &b_tpTo_hltIter0IterL3_eta);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_phi", &tpTo_hltIter0IterL3_phi, &b_tpTo_hltIter0IterL3_phi);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_parentVx", &tpTo_hltIter0IterL3_parentVx, &b_tpTo_hltIter0IterL3_parentVx);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_parentVy", &tpTo_hltIter0IterL3_parentVy, &b_tpTo_hltIter0IterL3_parentVy);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_parentVz", &tpTo_hltIter0IterL3_parentVz, &b_tpTo_hltIter0IterL3_parentVz);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_status", &tpTo_hltIter0IterL3_status, &b_tpTo_hltIter0IterL3_status);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_numberOfHits", &tpTo_hltIter0IterL3_numberOfHits, &b_tpTo_hltIter0IterL3_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_numberOfTrackerHits", &tpTo_hltIter0IterL3_numberOfTrackerHits, &b_tpTo_hltIter0IterL3_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_numberOfTrackerLayers", &tpTo_hltIter0IterL3_numberOfTrackerLayers, &b_tpTo_hltIter0IterL3_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_gen_charge", &tpTo_hltIter0IterL3_gen_charge, &b_tpTo_hltIter0IterL3_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_gen_pdgId", &tpTo_hltIter0IterL3_gen_pdgId, &b_tpTo_hltIter0IterL3_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_gen_pt", &tpTo_hltIter0IterL3_gen_pt, &b_tpTo_hltIter0IterL3_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_gen_eta", &tpTo_hltIter0IterL3_gen_eta, &b_tpTo_hltIter0IterL3_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_gen_phi", &tpTo_hltIter0IterL3_gen_phi, &b_tpTo_hltIter0IterL3_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_bestMatchTrk_pt", &tpTo_hltIter0IterL3_bestMatchTrk_pt, &b_tpTo_hltIter0IterL3_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_bestMatchTrk_eta", &tpTo_hltIter0IterL3_bestMatchTrk_eta, &b_tpTo_hltIter0IterL3_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_bestMatchTrk_phi", &tpTo_hltIter0IterL3_bestMatchTrk_phi, &b_tpTo_hltIter0IterL3_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_bestMatchTrk_charge", &tpTo_hltIter0IterL3_bestMatchTrk_charge, &b_tpTo_hltIter0IterL3_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_bestMatchTrk_quality", &tpTo_hltIter0IterL3_bestMatchTrk_quality, &b_tpTo_hltIter0IterL3_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3_bestMatchTrk_NValidHits", &tpTo_hltIter0IterL3_bestMatchTrk_NValidHits, &b_tpTo_hltIter0IterL3_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltIter2IterL3", &nhltIter2IterL3, &b_nhltIter2IterL3);
            fChain->SetBranchAddress("hltIter2IterL3_pt", &hltIter2IterL3_pt, &b_hltIter2IterL3_pt);
            fChain->SetBranchAddress("hltIter2IterL3_ptError", &hltIter2IterL3_ptError, &b_hltIter2IterL3_ptError);
            fChain->SetBranchAddress("hltIter2IterL3_eta", &hltIter2IterL3_eta, &b_hltIter2IterL3_eta);
            fChain->SetBranchAddress("hltIter2IterL3_phi", &hltIter2IterL3_phi, &b_hltIter2IterL3_phi);
            fChain->SetBranchAddress("hltIter2IterL3_charge", &hltIter2IterL3_charge, &b_hltIter2IterL3_charge);
            fChain->SetBranchAddress("hltIter2IterL3_matchedL3", &hltIter2IterL3_matchedL3, &b_hltIter2IterL3_matchedL3);
            fChain->SetBranchAddress("hltIter2IterL3_matchedL3NoId", &hltIter2IterL3_matchedL3NoId, &b_hltIter2IterL3_matchedL3NoId);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_charge", &hltIter2IterL3_bestMatchTP_charge, &b_hltIter2IterL3_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_pdgId", &hltIter2IterL3_bestMatchTP_pdgId, &b_hltIter2IterL3_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_energy", &hltIter2IterL3_bestMatchTP_energy, &b_hltIter2IterL3_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_pt", &hltIter2IterL3_bestMatchTP_pt, &b_hltIter2IterL3_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_eta", &hltIter2IterL3_bestMatchTP_eta, &b_hltIter2IterL3_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_phi", &hltIter2IterL3_bestMatchTP_phi, &b_hltIter2IterL3_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_parentVx", &hltIter2IterL3_bestMatchTP_parentVx, &b_hltIter2IterL3_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_parentVy", &hltIter2IterL3_bestMatchTP_parentVy, &b_hltIter2IterL3_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_parentVz", &hltIter2IterL3_bestMatchTP_parentVz, &b_hltIter2IterL3_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_status", &hltIter2IterL3_bestMatchTP_status, &b_hltIter2IterL3_bestMatchTP_status);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_numberOfHits", &hltIter2IterL3_bestMatchTP_numberOfHits, &b_hltIter2IterL3_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_numberOfTrackerHits", &hltIter2IterL3_bestMatchTP_numberOfTrackerHits, &b_hltIter2IterL3_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_numberOfTrackerLayers", &hltIter2IterL3_bestMatchTP_numberOfTrackerLayers, &b_hltIter2IterL3_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIter2IterL3_bestMatchTP_sharedFraction", &hltIter2IterL3_bestMatchTP_sharedFraction, &b_hltIter2IterL3_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIter2IterL3_matchedTPsize", &hltIter2IterL3_matchedTPsize, &b_hltIter2IterL3_matchedTPsize);
            fChain->SetBranchAddress("hltIter2IterL3_mva0", &hltIter2IterL3_mva0, &b_hltIter2IterL3_mva0);
            fChain->SetBranchAddress("hltIter2IterL3_mva1", &hltIter2IterL3_mva1, &b_hltIter2IterL3_mva1);
            fChain->SetBranchAddress("hltIter2IterL3_mva2", &hltIter2IterL3_mva2, &b_hltIter2IterL3_mva2);
            fChain->SetBranchAddress("hltIter2IterL3_mva3", &hltIter2IterL3_mva3, &b_hltIter2IterL3_mva3);
            fChain->SetBranchAddress("ntpTo_hltIter2IterL3", &ntpTo_hltIter2IterL3, &b_ntpTo_hltIter2IterL3);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_charge", &tpTo_hltIter2IterL3_charge, &b_tpTo_hltIter2IterL3_charge);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_pdgId", &tpTo_hltIter2IterL3_pdgId, &b_tpTo_hltIter2IterL3_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_energy", &tpTo_hltIter2IterL3_energy, &b_tpTo_hltIter2IterL3_energy);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_pt", &tpTo_hltIter2IterL3_pt, &b_tpTo_hltIter2IterL3_pt);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_eta", &tpTo_hltIter2IterL3_eta, &b_tpTo_hltIter2IterL3_eta);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_phi", &tpTo_hltIter2IterL3_phi, &b_tpTo_hltIter2IterL3_phi);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_parentVx", &tpTo_hltIter2IterL3_parentVx, &b_tpTo_hltIter2IterL3_parentVx);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_parentVy", &tpTo_hltIter2IterL3_parentVy, &b_tpTo_hltIter2IterL3_parentVy);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_parentVz", &tpTo_hltIter2IterL3_parentVz, &b_tpTo_hltIter2IterL3_parentVz);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_status", &tpTo_hltIter2IterL3_status, &b_tpTo_hltIter2IterL3_status);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_numberOfHits", &tpTo_hltIter2IterL3_numberOfHits, &b_tpTo_hltIter2IterL3_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_numberOfTrackerHits", &tpTo_hltIter2IterL3_numberOfTrackerHits, &b_tpTo_hltIter2IterL3_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_numberOfTrackerLayers", &tpTo_hltIter2IterL3_numberOfTrackerLayers, &b_tpTo_hltIter2IterL3_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_gen_charge", &tpTo_hltIter2IterL3_gen_charge, &b_tpTo_hltIter2IterL3_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_gen_pdgId", &tpTo_hltIter2IterL3_gen_pdgId, &b_tpTo_hltIter2IterL3_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_gen_pt", &tpTo_hltIter2IterL3_gen_pt, &b_tpTo_hltIter2IterL3_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_gen_eta", &tpTo_hltIter2IterL3_gen_eta, &b_tpTo_hltIter2IterL3_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_gen_phi", &tpTo_hltIter2IterL3_gen_phi, &b_tpTo_hltIter2IterL3_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_bestMatchTrk_pt", &tpTo_hltIter2IterL3_bestMatchTrk_pt, &b_tpTo_hltIter2IterL3_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_bestMatchTrk_eta", &tpTo_hltIter2IterL3_bestMatchTrk_eta, &b_tpTo_hltIter2IterL3_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_bestMatchTrk_phi", &tpTo_hltIter2IterL3_bestMatchTrk_phi, &b_tpTo_hltIter2IterL3_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_bestMatchTrk_charge", &tpTo_hltIter2IterL3_bestMatchTrk_charge, &b_tpTo_hltIter2IterL3_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_bestMatchTrk_quality", &tpTo_hltIter2IterL3_bestMatchTrk_quality, &b_tpTo_hltIter2IterL3_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3_bestMatchTrk_NValidHits", &tpTo_hltIter2IterL3_bestMatchTrk_NValidHits, &b_tpTo_hltIter2IterL3_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltIter0IterL3FromL1Muon", &nhltIter0IterL3FromL1Muon, &b_nhltIter0IterL3FromL1Muon);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_pt", &hltIter0IterL3FromL1Muon_pt, &b_hltIter0IterL3FromL1Muon_pt);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_ptError", &hltIter0IterL3FromL1Muon_ptError, &b_hltIter0IterL3FromL1Muon_ptError);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_eta", &hltIter0IterL3FromL1Muon_eta, &b_hltIter0IterL3FromL1Muon_eta);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_phi", &hltIter0IterL3FromL1Muon_phi, &b_hltIter0IterL3FromL1Muon_phi);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_charge", &hltIter0IterL3FromL1Muon_charge, &b_hltIter0IterL3FromL1Muon_charge);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_matchedL3", &hltIter0IterL3FromL1Muon_matchedL3, &b_hltIter0IterL3FromL1Muon_matchedL3);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_matchedL3NoId", &hltIter0IterL3FromL1Muon_matchedL3NoId, &b_hltIter0IterL3FromL1Muon_matchedL3NoId);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_charge", &hltIter0IterL3FromL1Muon_bestMatchTP_charge, &b_hltIter0IterL3FromL1Muon_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_pdgId", &hltIter0IterL3FromL1Muon_bestMatchTP_pdgId, &b_hltIter0IterL3FromL1Muon_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_energy", &hltIter0IterL3FromL1Muon_bestMatchTP_energy, &b_hltIter0IterL3FromL1Muon_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_pt", &hltIter0IterL3FromL1Muon_bestMatchTP_pt, &b_hltIter0IterL3FromL1Muon_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_eta", &hltIter0IterL3FromL1Muon_bestMatchTP_eta, &b_hltIter0IterL3FromL1Muon_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_phi", &hltIter0IterL3FromL1Muon_bestMatchTP_phi, &b_hltIter0IterL3FromL1Muon_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_parentVx", &hltIter0IterL3FromL1Muon_bestMatchTP_parentVx, &b_hltIter0IterL3FromL1Muon_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_parentVy", &hltIter0IterL3FromL1Muon_bestMatchTP_parentVy, &b_hltIter0IterL3FromL1Muon_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_parentVz", &hltIter0IterL3FromL1Muon_bestMatchTP_parentVz, &b_hltIter0IterL3FromL1Muon_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_status", &hltIter0IterL3FromL1Muon_bestMatchTP_status, &b_hltIter0IterL3FromL1Muon_bestMatchTP_status);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_numberOfHits", &hltIter0IterL3FromL1Muon_bestMatchTP_numberOfHits, &b_hltIter0IterL3FromL1Muon_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits", &hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits, &b_hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers", &hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers, &b_hltIter0IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_bestMatchTP_sharedFraction", &hltIter0IterL3FromL1Muon_bestMatchTP_sharedFraction, &b_hltIter0IterL3FromL1Muon_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_matchedTPsize", &hltIter0IterL3FromL1Muon_matchedTPsize, &b_hltIter0IterL3FromL1Muon_matchedTPsize);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_mva0", &hltIter0IterL3FromL1Muon_mva0, &b_hltIter0IterL3FromL1Muon_mva0);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_mva1", &hltIter0IterL3FromL1Muon_mva1, &b_hltIter0IterL3FromL1Muon_mva1);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_mva2", &hltIter0IterL3FromL1Muon_mva2, &b_hltIter0IterL3FromL1Muon_mva2);
            fChain->SetBranchAddress("hltIter0IterL3FromL1Muon_mva3", &hltIter0IterL3FromL1Muon_mva3, &b_hltIter0IterL3FromL1Muon_mva3);
            fChain->SetBranchAddress("ntpTo_hltIter0IterL3FromL1Muon", &ntpTo_hltIter0IterL3FromL1Muon, &b_ntpTo_hltIter0IterL3FromL1Muon);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_charge", &tpTo_hltIter0IterL3FromL1Muon_charge, &b_tpTo_hltIter0IterL3FromL1Muon_charge);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_pdgId", &tpTo_hltIter0IterL3FromL1Muon_pdgId, &b_tpTo_hltIter0IterL3FromL1Muon_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_energy", &tpTo_hltIter0IterL3FromL1Muon_energy, &b_tpTo_hltIter0IterL3FromL1Muon_energy);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_pt", &tpTo_hltIter0IterL3FromL1Muon_pt, &b_tpTo_hltIter0IterL3FromL1Muon_pt);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_eta", &tpTo_hltIter0IterL3FromL1Muon_eta, &b_tpTo_hltIter0IterL3FromL1Muon_eta);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_phi", &tpTo_hltIter0IterL3FromL1Muon_phi, &b_tpTo_hltIter0IterL3FromL1Muon_phi);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_parentVx", &tpTo_hltIter0IterL3FromL1Muon_parentVx, &b_tpTo_hltIter0IterL3FromL1Muon_parentVx);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_parentVy", &tpTo_hltIter0IterL3FromL1Muon_parentVy, &b_tpTo_hltIter0IterL3FromL1Muon_parentVy);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_parentVz", &tpTo_hltIter0IterL3FromL1Muon_parentVz, &b_tpTo_hltIter0IterL3FromL1Muon_parentVz);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_status", &tpTo_hltIter0IterL3FromL1Muon_status, &b_tpTo_hltIter0IterL3FromL1Muon_status);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_numberOfHits", &tpTo_hltIter0IterL3FromL1Muon_numberOfHits, &b_tpTo_hltIter0IterL3FromL1Muon_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerHits", &tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerHits, &b_tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerLayers", &tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerLayers, &b_tpTo_hltIter0IterL3FromL1Muon_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_gen_charge", &tpTo_hltIter0IterL3FromL1Muon_gen_charge, &b_tpTo_hltIter0IterL3FromL1Muon_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_gen_pdgId", &tpTo_hltIter0IterL3FromL1Muon_gen_pdgId, &b_tpTo_hltIter0IterL3FromL1Muon_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_gen_pt", &tpTo_hltIter0IterL3FromL1Muon_gen_pt, &b_tpTo_hltIter0IterL3FromL1Muon_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_gen_eta", &tpTo_hltIter0IterL3FromL1Muon_gen_eta, &b_tpTo_hltIter0IterL3FromL1Muon_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_gen_phi", &tpTo_hltIter0IterL3FromL1Muon_gen_phi, &b_tpTo_hltIter0IterL3FromL1Muon_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_pt", &tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_pt, &b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_eta", &tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_eta, &b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_phi", &tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_phi, &b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_charge", &tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_charge, &b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_quality", &tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_quality, &b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_NValidHits", &tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_NValidHits, &b_tpTo_hltIter0IterL3FromL1Muon_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltIter2IterL3FromL1Muon", &nhltIter2IterL3FromL1Muon, &b_nhltIter2IterL3FromL1Muon);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_pt", &hltIter2IterL3FromL1Muon_pt, &b_hltIter2IterL3FromL1Muon_pt);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_ptError", &hltIter2IterL3FromL1Muon_ptError, &b_hltIter2IterL3FromL1Muon_ptError);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_eta", &hltIter2IterL3FromL1Muon_eta, &b_hltIter2IterL3FromL1Muon_eta);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_phi", &hltIter2IterL3FromL1Muon_phi, &b_hltIter2IterL3FromL1Muon_phi);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_charge", &hltIter2IterL3FromL1Muon_charge, &b_hltIter2IterL3FromL1Muon_charge);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_matchedL3", &hltIter2IterL3FromL1Muon_matchedL3, &b_hltIter2IterL3FromL1Muon_matchedL3);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_matchedL3NoId", &hltIter2IterL3FromL1Muon_matchedL3NoId, &b_hltIter2IterL3FromL1Muon_matchedL3NoId);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_charge", &hltIter2IterL3FromL1Muon_bestMatchTP_charge, &b_hltIter2IterL3FromL1Muon_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_pdgId", &hltIter2IterL3FromL1Muon_bestMatchTP_pdgId, &b_hltIter2IterL3FromL1Muon_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_energy", &hltIter2IterL3FromL1Muon_bestMatchTP_energy, &b_hltIter2IterL3FromL1Muon_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_pt", &hltIter2IterL3FromL1Muon_bestMatchTP_pt, &b_hltIter2IterL3FromL1Muon_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_eta", &hltIter2IterL3FromL1Muon_bestMatchTP_eta, &b_hltIter2IterL3FromL1Muon_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_phi", &hltIter2IterL3FromL1Muon_bestMatchTP_phi, &b_hltIter2IterL3FromL1Muon_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_parentVx", &hltIter2IterL3FromL1Muon_bestMatchTP_parentVx, &b_hltIter2IterL3FromL1Muon_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_parentVy", &hltIter2IterL3FromL1Muon_bestMatchTP_parentVy, &b_hltIter2IterL3FromL1Muon_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_parentVz", &hltIter2IterL3FromL1Muon_bestMatchTP_parentVz, &b_hltIter2IterL3FromL1Muon_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_status", &hltIter2IterL3FromL1Muon_bestMatchTP_status, &b_hltIter2IterL3FromL1Muon_bestMatchTP_status);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_numberOfHits", &hltIter2IterL3FromL1Muon_bestMatchTP_numberOfHits, &b_hltIter2IterL3FromL1Muon_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits", &hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits, &b_hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers", &hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers, &b_hltIter2IterL3FromL1Muon_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_bestMatchTP_sharedFraction", &hltIter2IterL3FromL1Muon_bestMatchTP_sharedFraction, &b_hltIter2IterL3FromL1Muon_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_matchedTPsize", &hltIter2IterL3FromL1Muon_matchedTPsize, &b_hltIter2IterL3FromL1Muon_matchedTPsize);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_mva0", &hltIter2IterL3FromL1Muon_mva0, &b_hltIter2IterL3FromL1Muon_mva0);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_mva1", &hltIter2IterL3FromL1Muon_mva1, &b_hltIter2IterL3FromL1Muon_mva1);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_mva2", &hltIter2IterL3FromL1Muon_mva2, &b_hltIter2IterL3FromL1Muon_mva2);
            fChain->SetBranchAddress("hltIter2IterL3FromL1Muon_mva3", &hltIter2IterL3FromL1Muon_mva3, &b_hltIter2IterL3FromL1Muon_mva3);
            fChain->SetBranchAddress("ntpTo_hltIter2IterL3FromL1Muon", &ntpTo_hltIter2IterL3FromL1Muon, &b_ntpTo_hltIter2IterL3FromL1Muon);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_charge", &tpTo_hltIter2IterL3FromL1Muon_charge, &b_tpTo_hltIter2IterL3FromL1Muon_charge);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_pdgId", &tpTo_hltIter2IterL3FromL1Muon_pdgId, &b_tpTo_hltIter2IterL3FromL1Muon_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_energy", &tpTo_hltIter2IterL3FromL1Muon_energy, &b_tpTo_hltIter2IterL3FromL1Muon_energy);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_pt", &tpTo_hltIter2IterL3FromL1Muon_pt, &b_tpTo_hltIter2IterL3FromL1Muon_pt);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_eta", &tpTo_hltIter2IterL3FromL1Muon_eta, &b_tpTo_hltIter2IterL3FromL1Muon_eta);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_phi", &tpTo_hltIter2IterL3FromL1Muon_phi, &b_tpTo_hltIter2IterL3FromL1Muon_phi);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_parentVx", &tpTo_hltIter2IterL3FromL1Muon_parentVx, &b_tpTo_hltIter2IterL3FromL1Muon_parentVx);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_parentVy", &tpTo_hltIter2IterL3FromL1Muon_parentVy, &b_tpTo_hltIter2IterL3FromL1Muon_parentVy);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_parentVz", &tpTo_hltIter2IterL3FromL1Muon_parentVz, &b_tpTo_hltIter2IterL3FromL1Muon_parentVz);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_status", &tpTo_hltIter2IterL3FromL1Muon_status, &b_tpTo_hltIter2IterL3FromL1Muon_status);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_numberOfHits", &tpTo_hltIter2IterL3FromL1Muon_numberOfHits, &b_tpTo_hltIter2IterL3FromL1Muon_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerHits", &tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerHits, &b_tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerLayers", &tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerLayers, &b_tpTo_hltIter2IterL3FromL1Muon_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_gen_charge", &tpTo_hltIter2IterL3FromL1Muon_gen_charge, &b_tpTo_hltIter2IterL3FromL1Muon_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_gen_pdgId", &tpTo_hltIter2IterL3FromL1Muon_gen_pdgId, &b_tpTo_hltIter2IterL3FromL1Muon_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_gen_pt", &tpTo_hltIter2IterL3FromL1Muon_gen_pt, &b_tpTo_hltIter2IterL3FromL1Muon_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_gen_eta", &tpTo_hltIter2IterL3FromL1Muon_gen_eta, &b_tpTo_hltIter2IterL3FromL1Muon_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_gen_phi", &tpTo_hltIter2IterL3FromL1Muon_gen_phi, &b_tpTo_hltIter2IterL3FromL1Muon_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_pt", &tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_pt, &b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_eta", &tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_eta, &b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_phi", &tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_phi, &b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_charge", &tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_charge, &b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_quality", &tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_quality, &b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_NValidHits", &tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_NValidHits, &b_tpTo_hltIter2IterL3FromL1Muon_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltIterL3IOFromL1", &nhltIterL3IOFromL1, &b_nhltIterL3IOFromL1);
            fChain->SetBranchAddress("hltIterL3IOFromL1_pt", &hltIterL3IOFromL1_pt, &b_hltIterL3IOFromL1_pt);
            fChain->SetBranchAddress("hltIterL3IOFromL1_ptError", &hltIterL3IOFromL1_ptError, &b_hltIterL3IOFromL1_ptError);
            fChain->SetBranchAddress("hltIterL3IOFromL1_eta", &hltIterL3IOFromL1_eta, &b_hltIterL3IOFromL1_eta);
            fChain->SetBranchAddress("hltIterL3IOFromL1_phi", &hltIterL3IOFromL1_phi, &b_hltIterL3IOFromL1_phi);
            fChain->SetBranchAddress("hltIterL3IOFromL1_charge", &hltIterL3IOFromL1_charge, &b_hltIterL3IOFromL1_charge);
            fChain->SetBranchAddress("hltIterL3IOFromL1_matchedL3", &hltIterL3IOFromL1_matchedL3, &b_hltIterL3IOFromL1_matchedL3);
            fChain->SetBranchAddress("hltIterL3IOFromL1_matchedL3NoId", &hltIterL3IOFromL1_matchedL3NoId, &b_hltIterL3IOFromL1_matchedL3NoId);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_charge", &hltIterL3IOFromL1_bestMatchTP_charge, &b_hltIterL3IOFromL1_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_pdgId", &hltIterL3IOFromL1_bestMatchTP_pdgId, &b_hltIterL3IOFromL1_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_energy", &hltIterL3IOFromL1_bestMatchTP_energy, &b_hltIterL3IOFromL1_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_pt", &hltIterL3IOFromL1_bestMatchTP_pt, &b_hltIterL3IOFromL1_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_eta", &hltIterL3IOFromL1_bestMatchTP_eta, &b_hltIterL3IOFromL1_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_phi", &hltIterL3IOFromL1_bestMatchTP_phi, &b_hltIterL3IOFromL1_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_parentVx", &hltIterL3IOFromL1_bestMatchTP_parentVx, &b_hltIterL3IOFromL1_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_parentVy", &hltIterL3IOFromL1_bestMatchTP_parentVy, &b_hltIterL3IOFromL1_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_parentVz", &hltIterL3IOFromL1_bestMatchTP_parentVz, &b_hltIterL3IOFromL1_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_status", &hltIterL3IOFromL1_bestMatchTP_status, &b_hltIterL3IOFromL1_bestMatchTP_status);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_numberOfHits", &hltIterL3IOFromL1_bestMatchTP_numberOfHits, &b_hltIterL3IOFromL1_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_numberOfTrackerHits", &hltIterL3IOFromL1_bestMatchTP_numberOfTrackerHits, &b_hltIterL3IOFromL1_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_numberOfTrackerLayers", &hltIterL3IOFromL1_bestMatchTP_numberOfTrackerLayers, &b_hltIterL3IOFromL1_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIterL3IOFromL1_bestMatchTP_sharedFraction", &hltIterL3IOFromL1_bestMatchTP_sharedFraction, &b_hltIterL3IOFromL1_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIterL3IOFromL1_matchedTPsize", &hltIterL3IOFromL1_matchedTPsize, &b_hltIterL3IOFromL1_matchedTPsize);
            fChain->SetBranchAddress("hltIterL3IOFromL1_mva0", &hltIterL3IOFromL1_mva0, &b_hltIterL3IOFromL1_mva0);
            fChain->SetBranchAddress("hltIterL3IOFromL1_mva1", &hltIterL3IOFromL1_mva1, &b_hltIterL3IOFromL1_mva1);
            fChain->SetBranchAddress("hltIterL3IOFromL1_mva2", &hltIterL3IOFromL1_mva2, &b_hltIterL3IOFromL1_mva2);
            fChain->SetBranchAddress("hltIterL3IOFromL1_mva3", &hltIterL3IOFromL1_mva3, &b_hltIterL3IOFromL1_mva3);
            fChain->SetBranchAddress("ntpTo_hltIterL3IOFromL1", &ntpTo_hltIterL3IOFromL1, &b_ntpTo_hltIterL3IOFromL1);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_charge", &tpTo_hltIterL3IOFromL1_charge, &b_tpTo_hltIterL3IOFromL1_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_pdgId", &tpTo_hltIterL3IOFromL1_pdgId, &b_tpTo_hltIterL3IOFromL1_pdgId);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_energy", &tpTo_hltIterL3IOFromL1_energy, &b_tpTo_hltIterL3IOFromL1_energy);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_pt", &tpTo_hltIterL3IOFromL1_pt, &b_tpTo_hltIterL3IOFromL1_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_eta", &tpTo_hltIterL3IOFromL1_eta, &b_tpTo_hltIterL3IOFromL1_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_phi", &tpTo_hltIterL3IOFromL1_phi, &b_tpTo_hltIterL3IOFromL1_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_parentVx", &tpTo_hltIterL3IOFromL1_parentVx, &b_tpTo_hltIterL3IOFromL1_parentVx);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_parentVy", &tpTo_hltIterL3IOFromL1_parentVy, &b_tpTo_hltIterL3IOFromL1_parentVy);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_parentVz", &tpTo_hltIterL3IOFromL1_parentVz, &b_tpTo_hltIterL3IOFromL1_parentVz);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_status", &tpTo_hltIterL3IOFromL1_status, &b_tpTo_hltIterL3IOFromL1_status);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_numberOfHits", &tpTo_hltIterL3IOFromL1_numberOfHits, &b_tpTo_hltIterL3IOFromL1_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_numberOfTrackerHits", &tpTo_hltIterL3IOFromL1_numberOfTrackerHits, &b_tpTo_hltIterL3IOFromL1_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_numberOfTrackerLayers", &tpTo_hltIterL3IOFromL1_numberOfTrackerLayers, &b_tpTo_hltIterL3IOFromL1_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_gen_charge", &tpTo_hltIterL3IOFromL1_gen_charge, &b_tpTo_hltIterL3IOFromL1_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_gen_pdgId", &tpTo_hltIterL3IOFromL1_gen_pdgId, &b_tpTo_hltIterL3IOFromL1_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_gen_pt", &tpTo_hltIterL3IOFromL1_gen_pt, &b_tpTo_hltIterL3IOFromL1_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_gen_eta", &tpTo_hltIterL3IOFromL1_gen_eta, &b_tpTo_hltIterL3IOFromL1_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_gen_phi", &tpTo_hltIterL3IOFromL1_gen_phi, &b_tpTo_hltIterL3IOFromL1_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_bestMatchTrk_pt", &tpTo_hltIterL3IOFromL1_bestMatchTrk_pt, &b_tpTo_hltIterL3IOFromL1_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_bestMatchTrk_eta", &tpTo_hltIterL3IOFromL1_bestMatchTrk_eta, &b_tpTo_hltIterL3IOFromL1_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_bestMatchTrk_phi", &tpTo_hltIterL3IOFromL1_bestMatchTrk_phi, &b_tpTo_hltIterL3IOFromL1_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_bestMatchTrk_charge", &tpTo_hltIterL3IOFromL1_bestMatchTrk_charge, &b_tpTo_hltIterL3IOFromL1_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_bestMatchTrk_quality", &tpTo_hltIterL3IOFromL1_bestMatchTrk_quality, &b_tpTo_hltIterL3IOFromL1_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIterL3IOFromL1_bestMatchTrk_NValidHits", &tpTo_hltIterL3IOFromL1_bestMatchTrk_NValidHits, &b_tpTo_hltIterL3IOFromL1_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltIterL3MuonsNoID", &nhltIterL3MuonsNoID, &b_nhltIterL3MuonsNoID);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_pt", &hltIterL3MuonsNoID_pt, &b_hltIterL3MuonsNoID_pt);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_ptError", &hltIterL3MuonsNoID_ptError, &b_hltIterL3MuonsNoID_ptError);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_eta", &hltIterL3MuonsNoID_eta, &b_hltIterL3MuonsNoID_eta);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_phi", &hltIterL3MuonsNoID_phi, &b_hltIterL3MuonsNoID_phi);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_charge", &hltIterL3MuonsNoID_charge, &b_hltIterL3MuonsNoID_charge);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_matchedL3", &hltIterL3MuonsNoID_matchedL3, &b_hltIterL3MuonsNoID_matchedL3);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_matchedL3NoId", &hltIterL3MuonsNoID_matchedL3NoId, &b_hltIterL3MuonsNoID_matchedL3NoId);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_charge", &hltIterL3MuonsNoID_bestMatchTP_charge, &b_hltIterL3MuonsNoID_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_pdgId", &hltIterL3MuonsNoID_bestMatchTP_pdgId, &b_hltIterL3MuonsNoID_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_energy", &hltIterL3MuonsNoID_bestMatchTP_energy, &b_hltIterL3MuonsNoID_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_pt", &hltIterL3MuonsNoID_bestMatchTP_pt, &b_hltIterL3MuonsNoID_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_eta", &hltIterL3MuonsNoID_bestMatchTP_eta, &b_hltIterL3MuonsNoID_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_phi", &hltIterL3MuonsNoID_bestMatchTP_phi, &b_hltIterL3MuonsNoID_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_parentVx", &hltIterL3MuonsNoID_bestMatchTP_parentVx, &b_hltIterL3MuonsNoID_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_parentVy", &hltIterL3MuonsNoID_bestMatchTP_parentVy, &b_hltIterL3MuonsNoID_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_parentVz", &hltIterL3MuonsNoID_bestMatchTP_parentVz, &b_hltIterL3MuonsNoID_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_status", &hltIterL3MuonsNoID_bestMatchTP_status, &b_hltIterL3MuonsNoID_bestMatchTP_status);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_numberOfHits", &hltIterL3MuonsNoID_bestMatchTP_numberOfHits, &b_hltIterL3MuonsNoID_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerHits", &hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerHits, &b_hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerLayers", &hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerLayers, &b_hltIterL3MuonsNoID_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_bestMatchTP_sharedFraction", &hltIterL3MuonsNoID_bestMatchTP_sharedFraction, &b_hltIterL3MuonsNoID_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_matchedTPsize", &hltIterL3MuonsNoID_matchedTPsize, &b_hltIterL3MuonsNoID_matchedTPsize);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_mva0", &hltIterL3MuonsNoID_mva0, &b_hltIterL3MuonsNoID_mva0);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_mva1", &hltIterL3MuonsNoID_mva1, &b_hltIterL3MuonsNoID_mva1);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_mva2", &hltIterL3MuonsNoID_mva2, &b_hltIterL3MuonsNoID_mva2);
            fChain->SetBranchAddress("hltIterL3MuonsNoID_mva3", &hltIterL3MuonsNoID_mva3, &b_hltIterL3MuonsNoID_mva3);
            fChain->SetBranchAddress("ntpTo_hltIterL3MuonsNoID", &ntpTo_hltIterL3MuonsNoID, &b_ntpTo_hltIterL3MuonsNoID);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_charge", &tpTo_hltIterL3MuonsNoID_charge, &b_tpTo_hltIterL3MuonsNoID_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_pdgId", &tpTo_hltIterL3MuonsNoID_pdgId, &b_tpTo_hltIterL3MuonsNoID_pdgId);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_energy", &tpTo_hltIterL3MuonsNoID_energy, &b_tpTo_hltIterL3MuonsNoID_energy);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_pt", &tpTo_hltIterL3MuonsNoID_pt, &b_tpTo_hltIterL3MuonsNoID_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_eta", &tpTo_hltIterL3MuonsNoID_eta, &b_tpTo_hltIterL3MuonsNoID_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_phi", &tpTo_hltIterL3MuonsNoID_phi, &b_tpTo_hltIterL3MuonsNoID_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_parentVx", &tpTo_hltIterL3MuonsNoID_parentVx, &b_tpTo_hltIterL3MuonsNoID_parentVx);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_parentVy", &tpTo_hltIterL3MuonsNoID_parentVy, &b_tpTo_hltIterL3MuonsNoID_parentVy);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_parentVz", &tpTo_hltIterL3MuonsNoID_parentVz, &b_tpTo_hltIterL3MuonsNoID_parentVz);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_status", &tpTo_hltIterL3MuonsNoID_status, &b_tpTo_hltIterL3MuonsNoID_status);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_numberOfHits", &tpTo_hltIterL3MuonsNoID_numberOfHits, &b_tpTo_hltIterL3MuonsNoID_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_numberOfTrackerHits", &tpTo_hltIterL3MuonsNoID_numberOfTrackerHits, &b_tpTo_hltIterL3MuonsNoID_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_numberOfTrackerLayers", &tpTo_hltIterL3MuonsNoID_numberOfTrackerLayers, &b_tpTo_hltIterL3MuonsNoID_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_gen_charge", &tpTo_hltIterL3MuonsNoID_gen_charge, &b_tpTo_hltIterL3MuonsNoID_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_gen_pdgId", &tpTo_hltIterL3MuonsNoID_gen_pdgId, &b_tpTo_hltIterL3MuonsNoID_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_gen_pt", &tpTo_hltIterL3MuonsNoID_gen_pt, &b_tpTo_hltIterL3MuonsNoID_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_gen_eta", &tpTo_hltIterL3MuonsNoID_gen_eta, &b_tpTo_hltIterL3MuonsNoID_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_gen_phi", &tpTo_hltIterL3MuonsNoID_gen_phi, &b_tpTo_hltIterL3MuonsNoID_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_bestMatchTrk_pt", &tpTo_hltIterL3MuonsNoID_bestMatchTrk_pt, &b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_bestMatchTrk_eta", &tpTo_hltIterL3MuonsNoID_bestMatchTrk_eta, &b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_bestMatchTrk_phi", &tpTo_hltIterL3MuonsNoID_bestMatchTrk_phi, &b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_bestMatchTrk_charge", &tpTo_hltIterL3MuonsNoID_bestMatchTrk_charge, &b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_bestMatchTrk_quality", &tpTo_hltIterL3MuonsNoID_bestMatchTrk_quality, &b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIterL3MuonsNoID_bestMatchTrk_NValidHits", &tpTo_hltIterL3MuonsNoID_bestMatchTrk_NValidHits, &b_tpTo_hltIterL3MuonsNoID_bestMatchTrk_NValidHits);
            fChain->SetBranchAddress("nhltIterL3Muons", &nhltIterL3Muons, &b_nhltIterL3Muons);
            fChain->SetBranchAddress("hltIterL3Muons_pt", &hltIterL3Muons_pt, &b_hltIterL3Muons_pt);
            fChain->SetBranchAddress("hltIterL3Muons_ptError", &hltIterL3Muons_ptError, &b_hltIterL3Muons_ptError);
            fChain->SetBranchAddress("hltIterL3Muons_eta", &hltIterL3Muons_eta, &b_hltIterL3Muons_eta);
            fChain->SetBranchAddress("hltIterL3Muons_phi", &hltIterL3Muons_phi, &b_hltIterL3Muons_phi);
            fChain->SetBranchAddress("hltIterL3Muons_charge", &hltIterL3Muons_charge, &b_hltIterL3Muons_charge);
            fChain->SetBranchAddress("hltIterL3Muons_matchedL3", &hltIterL3Muons_matchedL3, &b_hltIterL3Muons_matchedL3);
            fChain->SetBranchAddress("hltIterL3Muons_matchedL3NoId", &hltIterL3Muons_matchedL3NoId, &b_hltIterL3Muons_matchedL3NoId);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_charge", &hltIterL3Muons_bestMatchTP_charge, &b_hltIterL3Muons_bestMatchTP_charge);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_pdgId", &hltIterL3Muons_bestMatchTP_pdgId, &b_hltIterL3Muons_bestMatchTP_pdgId);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_energy", &hltIterL3Muons_bestMatchTP_energy, &b_hltIterL3Muons_bestMatchTP_energy);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_pt", &hltIterL3Muons_bestMatchTP_pt, &b_hltIterL3Muons_bestMatchTP_pt);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_eta", &hltIterL3Muons_bestMatchTP_eta, &b_hltIterL3Muons_bestMatchTP_eta);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_phi", &hltIterL3Muons_bestMatchTP_phi, &b_hltIterL3Muons_bestMatchTP_phi);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_parentVx", &hltIterL3Muons_bestMatchTP_parentVx, &b_hltIterL3Muons_bestMatchTP_parentVx);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_parentVy", &hltIterL3Muons_bestMatchTP_parentVy, &b_hltIterL3Muons_bestMatchTP_parentVy);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_parentVz", &hltIterL3Muons_bestMatchTP_parentVz, &b_hltIterL3Muons_bestMatchTP_parentVz);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_status", &hltIterL3Muons_bestMatchTP_status, &b_hltIterL3Muons_bestMatchTP_status);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_numberOfHits", &hltIterL3Muons_bestMatchTP_numberOfHits, &b_hltIterL3Muons_bestMatchTP_numberOfHits);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_numberOfTrackerHits", &hltIterL3Muons_bestMatchTP_numberOfTrackerHits, &b_hltIterL3Muons_bestMatchTP_numberOfTrackerHits);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_numberOfTrackerLayers", &hltIterL3Muons_bestMatchTP_numberOfTrackerLayers, &b_hltIterL3Muons_bestMatchTP_numberOfTrackerLayers);
            fChain->SetBranchAddress("hltIterL3Muons_bestMatchTP_sharedFraction", &hltIterL3Muons_bestMatchTP_sharedFraction, &b_hltIterL3Muons_bestMatchTP_sharedFraction);
            fChain->SetBranchAddress("hltIterL3Muons_matchedTPsize", &hltIterL3Muons_matchedTPsize, &b_hltIterL3Muons_matchedTPsize);
            fChain->SetBranchAddress("hltIterL3Muons_mva0", &hltIterL3Muons_mva0, &b_hltIterL3Muons_mva0);
            fChain->SetBranchAddress("hltIterL3Muons_mva1", &hltIterL3Muons_mva1, &b_hltIterL3Muons_mva1);
            fChain->SetBranchAddress("hltIterL3Muons_mva2", &hltIterL3Muons_mva2, &b_hltIterL3Muons_mva2);
            fChain->SetBranchAddress("hltIterL3Muons_mva3", &hltIterL3Muons_mva3, &b_hltIterL3Muons_mva3);
            fChain->SetBranchAddress("ntpTo_hltIterL3Muons", &ntpTo_hltIterL3Muons, &b_ntpTo_hltIterL3Muons);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_charge", &tpTo_hltIterL3Muons_charge, &b_tpTo_hltIterL3Muons_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_pdgId", &tpTo_hltIterL3Muons_pdgId, &b_tpTo_hltIterL3Muons_pdgId);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_energy", &tpTo_hltIterL3Muons_energy, &b_tpTo_hltIterL3Muons_energy);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_pt", &tpTo_hltIterL3Muons_pt, &b_tpTo_hltIterL3Muons_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_eta", &tpTo_hltIterL3Muons_eta, &b_tpTo_hltIterL3Muons_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_phi", &tpTo_hltIterL3Muons_phi, &b_tpTo_hltIterL3Muons_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_parentVx", &tpTo_hltIterL3Muons_parentVx, &b_tpTo_hltIterL3Muons_parentVx);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_parentVy", &tpTo_hltIterL3Muons_parentVy, &b_tpTo_hltIterL3Muons_parentVy);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_parentVz", &tpTo_hltIterL3Muons_parentVz, &b_tpTo_hltIterL3Muons_parentVz);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_status", &tpTo_hltIterL3Muons_status, &b_tpTo_hltIterL3Muons_status);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_numberOfHits", &tpTo_hltIterL3Muons_numberOfHits, &b_tpTo_hltIterL3Muons_numberOfHits);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_numberOfTrackerHits", &tpTo_hltIterL3Muons_numberOfTrackerHits, &b_tpTo_hltIterL3Muons_numberOfTrackerHits);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_numberOfTrackerLayers", &tpTo_hltIterL3Muons_numberOfTrackerLayers, &b_tpTo_hltIterL3Muons_numberOfTrackerLayers);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_gen_charge", &tpTo_hltIterL3Muons_gen_charge, &b_tpTo_hltIterL3Muons_gen_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_gen_pdgId", &tpTo_hltIterL3Muons_gen_pdgId, &b_tpTo_hltIterL3Muons_gen_pdgId);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_gen_pt", &tpTo_hltIterL3Muons_gen_pt, &b_tpTo_hltIterL3Muons_gen_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_gen_eta", &tpTo_hltIterL3Muons_gen_eta, &b_tpTo_hltIterL3Muons_gen_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_gen_phi", &tpTo_hltIterL3Muons_gen_phi, &b_tpTo_hltIterL3Muons_gen_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_bestMatchTrk_pt", &tpTo_hltIterL3Muons_bestMatchTrk_pt, &b_tpTo_hltIterL3Muons_bestMatchTrk_pt);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_bestMatchTrk_eta", &tpTo_hltIterL3Muons_bestMatchTrk_eta, &b_tpTo_hltIterL3Muons_bestMatchTrk_eta);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_bestMatchTrk_phi", &tpTo_hltIterL3Muons_bestMatchTrk_phi, &b_tpTo_hltIterL3Muons_bestMatchTrk_phi);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_bestMatchTrk_charge", &tpTo_hltIterL3Muons_bestMatchTrk_charge, &b_tpTo_hltIterL3Muons_bestMatchTrk_charge);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_bestMatchTrk_quality", &tpTo_hltIterL3Muons_bestMatchTrk_quality, &b_tpTo_hltIterL3Muons_bestMatchTrk_quality);
            fChain->SetBranchAddress("tpTo_hltIterL3Muons_bestMatchTrk_NValidHits", &tpTo_hltIterL3Muons_bestMatchTrk_NValidHits, &b_tpTo_hltIterL3Muons_bestMatchTrk_NValidHits);
        }

    Notify();
}

Bool_t MuonHLTNtuple::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TChain in a TChain or when when a new TChain
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void MuonHLTNtuple::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

#endif
