// -- ntuple maker for Muon HLT study
// -- author: Kyeongpil Lee (Seoul National University, kplee@cern.ch)

#include "MuonHLTTool/MuonHLTNtupler/interface/MuonHLTNtupler.h"

// -------- For CMSSW_12 -----------
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
// ---------------------------------

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "CommonTools/Utils/interface/associationMapFilterValues.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"

using namespace std;
using namespace reco;
using namespace edm;


MuonHLTNtupler::MuonHLTNtupler(const edm::ParameterSet& iConfig):
DebugMode(iConfig.getParameter<bool>("DebugMode")),

// t_offlineMuon_       ( consumes< std::vector<reco::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("offlineMuon"       )) ),
t_offlineMuon_       ( consumes< edm::View<reco::Muon> >                  (iConfig.getUntrackedParameter<edm::InputTag>("offlineMuon"       )) ),
t_offlineVertex_     ( consumes< reco::VertexCollection >                 (iConfig.getUntrackedParameter<edm::InputTag>("offlineVertex"     )) ),
t_triggerResults_    ( consumes< edm::TriggerResults >                    (iConfig.getUntrackedParameter<edm::InputTag>("triggerResults"    )) ),
t_triggerEvent_      ( consumes< trigger::TriggerEvent >                  (iConfig.getUntrackedParameter<edm::InputTag>("triggerEvent"      )) ),
t_myTriggerResults_  ( consumes< edm::TriggerResults >                    (iConfig.getUntrackedParameter<edm::InputTag>("myTriggerResults"  )) ),
t_myTriggerEvent_    ( consumes< trigger::TriggerEvent >                  (iConfig.getUntrackedParameter<edm::InputTag>("myTriggerEvent"    )) ),

t_hltIterL3MuonTrimmedPixelVertices_       ( consumes< reco::VertexCollection >(iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3MuonTrimmedPixelVertices")) ),
t_hltIterL3FromL1MuonTrimmedPixelVertices_ ( consumes< reco::VertexCollection >(iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3FromL1MuonTrimmedPixelVertices")) ),

t_lumiScaler_        ( consumes< LumiScalersCollection >                  (iConfig.getUntrackedParameter<edm::InputTag>("lumiScaler"        )) ),
t_offlineLumiScaler_ ( consumes< LumiScalersCollection >                  (iConfig.getUntrackedParameter<edm::InputTag>("offlineLumiScaler" )) ),
t_PUSummaryInfo_     ( consumes< std::vector<PileupSummaryInfo> >         (iConfig.getUntrackedParameter<edm::InputTag>("PUSummaryInfo"     )) ),
t_genEventInfo_      ( consumes< GenEventInfoProduct >                    (iConfig.getUntrackedParameter<edm::InputTag>("genEventInfo"      )) ),
t_genParticle_       ( consumes< reco::GenParticleCollection >            (iConfig.getUntrackedParameter<edm::InputTag>("genParticle"       )) ),

trackerTopologyESToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>()),
trackerGeometryESToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
magFieldESToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
geomDetESToken_(esConsumes<GeometricDet, IdealGeometryRecord>()),
propagatorESToken_(esConsumes<Propagator, TrackingComponentsRecord>(edm::ESInputTag("", "PropagatorWithMaterialParabolicMf")))
{
  // -- generic muon-like collections (L2Muon, L3Muon, TkMuon, iterL3OI,
  //    iterL3IOFromL2, iterL3FromL2 by default -- see ntupler_cfi.py).
  //    Adding/renaming a collection is now a python-only change: extend the
  //    "muonCollections" VPSet, no C++ edit needed here.
  for (const auto& pset : iConfig.getUntrackedParameter<std::vector<edm::ParameterSet>>("muonCollections")) {
    MuonCollectionCfg cfg;
    cfg.name = pset.getParameter<std::string>("name");
    cfg.type = pset.getParameter<std::string>("type");
    cfg.tag  = pset.getParameter<edm::InputTag>("label");
    muonCollectionCfgs_.push_back(cfg);
    muonCollectionAdapters_.push_back(makeMuonCollectionAdapter(cfg.type, cfg.tag, consumesCollector()));
  }

  // -- generic trigger path/filter matching config -- see MuonHLTTriggerMatching.h
  triggerPathsCfg_   = iConfig.getUntrackedParameter<std::vector<std::string>>("triggerPaths");
  triggerFiltersCfg_ = iConfig.getUntrackedParameter<std::vector<std::string>>("triggerFilters");
  maxDR_             = iConfig.getUntrackedParameter<double>("maxDR");

  warnedCollectionInvalid_.assign(muonCollectionCfgs_.size(), false);
  for (const auto& cfg : muonCollectionCfgs_) {
    muonCollectionNamesForBranch_.push_back(cfg.name);
  }
}

void MuonHLTNtupler::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  Init();

  _trackerTopology = &iSetup.getData(trackerTopologyESToken_);
  _trackerGeometry = &iSetup.getData(trackerGeometryESToken_);
  //_geometryDetector = &iSetup.getData(geomDetESToken_);
  // _magneticField = &iSetup.getData(magFieldESToken_);

  // -- basic info.
  isRealData_ = iEvent.isRealData();

  runNum_       = iEvent.id().run();
  lumiBlockNum_ = iEvent.id().luminosityBlock();
  eventNum_     = iEvent.id().event();

  // -- vertex
  edm::Handle<reco::VertexCollection> h_offlineVertex;
  if( iEvent.getByToken(t_offlineVertex_, h_offlineVertex) )
  {
    int nGoodVtx = 0;
    for(reco::VertexCollection::const_iterator it = h_offlineVertex->begin(); it != h_offlineVertex->end(); ++it)
      if( it->isValid() ) nGoodVtx++;

    nVertex_ = nGoodVtx;
  }

  // -- hltIterL3MuonTrimmedPixelVertices - needed?
  edm::Handle<reco::VertexCollection> h_hltIterL3MuonTrimmedPixelVertices;
  if( iEvent.getByToken(t_hltIterL3MuonTrimmedPixelVertices_, h_hltIterL3MuonTrimmedPixelVertices) )
  {
    for(reco::VertexCollection::const_iterator it = h_hltIterL3MuonTrimmedPixelVertices->begin(); it != h_hltIterL3MuonTrimmedPixelVertices->end(); ++it)
      VThltIterL3MuonTrimmedPixelVertices->fill(*it);
  }
  
  // -- hltIterL3FromL1MuonTrimmedPixelVertices
  edm::Handle<reco::VertexCollection> h_hltIterL3FromL1MuonTrimmedPixelVertices;
  if( iEvent.getByToken(t_hltIterL3FromL1MuonTrimmedPixelVertices_, h_hltIterL3FromL1MuonTrimmedPixelVertices) )
  {
    for(reco::VertexCollection::const_iterator it = h_hltIterL3FromL1MuonTrimmedPixelVertices->begin(); it != h_hltIterL3FromL1MuonTrimmedPixelVertices->end(); ++it)
      VThltIterL3FromL1MuonTrimmedPixelVertices->fill(*it);
  }

  if( isRealData_ )
  {
    bunchID_ = iEvent.bunchCrossing();

    // -- lumi scaler @ HLT
    edm::Handle<LumiScalersCollection> h_lumiScaler;
    if( iEvent.getByToken(t_lumiScaler_, h_lumiScaler) && h_lumiScaler->begin() != h_lumiScaler->end() )
    {
      instLumi_  = h_lumiScaler->begin()->instantLumi();
      dataPU_    = h_lumiScaler->begin()->pileup();
      dataPURMS_ = h_lumiScaler->begin()->pileupRMS();
      bunchLumi_ = h_lumiScaler->begin()->bunchLumi();
    }

    // -- lumi scaler @ offline
    edm::Handle<LumiScalersCollection> h_offlineLumiScaler;
    if( iEvent.getByToken(t_offlineLumiScaler_, h_offlineLumiScaler) && h_offlineLumiScaler->begin() != h_offlineLumiScaler->end() )
    {
      offlineInstLumi_  = h_offlineLumiScaler->begin()->instantLumi();
      offlineDataPU_    = h_offlineLumiScaler->begin()->pileup();
      offlineDataPURMS_ = h_offlineLumiScaler->begin()->pileupRMS();
      offlineBunchLumi_ = h_offlineLumiScaler->begin()->bunchLumi();
    }
  }

  // -- True PU info: only for MC -- //
  if( !isRealData_ )
  {
    edm::Handle<std::vector< PileupSummaryInfo > > h_PUSummaryInfo;

    if( iEvent.getByToken(t_PUSummaryInfo_,h_PUSummaryInfo) )
    {
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = h_PUSummaryInfo->begin(); PVI != h_PUSummaryInfo->end(); ++PVI)
      {
        if(PVI->getBunchCrossing()==0)
        {
          truePU_ = PVI->getTrueNumInteractions();
          PU_pT_hats_ = PVI->getPU_pT_hats();
          continue;
        }
      } // -- end of PU iteration -- //
    } // -- end of if ( token exists )
  } // -- end of isMC -- //

  // -- fill each object -------------
  Fill_HLT(iEvent, false);
  Fill_Muon(iEvent);
  Fill_GenericMuonCollections(iEvent);
  if( !isRealData_ ) {
    Fill_GenParticle(iEvent);
  }

  ntuple_->Fill();
}

void MuonHLTNtupler::beginJob()
{
  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TTree>("ntuple","ntuple");

  Make_Branch();
}

void MuonHLTNtupler::Init()
{
  isRealData_ = false;

  runNum_       = -999;
  lumiBlockNum_ = -999;
  eventNum_     = 0;

  bunchID_ = -999;

  nVertex_ = -999;

  instLumi_  = -999;
  dataPU_    = -999;
  dataPURMS_ = -999;
  bunchLumi_ = -999;

  offlineInstLumi_  = -999;
  offlineDataPU_    = -999;
  offlineDataPURMS_ = -999;
  offlineBunchLumi_ = -999;

  truePU_ = -999;

  genEventWeight_ = -999;
  qScale_ = -999;

  PU_pT_hats_.clear();

  nGenParticle_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    genParticle_ID_[i] = -999;
    genParticle_status_[i] = -999;
    genParticle_mother_[i] = -999;

    genParticle_pt_[i]     = -999;
    genParticle_eta_[i]    = -999;
    genParticle_phi_[i]    = -999;
    genParticle_px_[i]     = -999;
    genParticle_py_[i]     = -999;
    genParticle_pz_[i]     = -999;
    genParticle_energy_[i] = -999;
    genParticle_charge_[i] = -999;

    genParticle_isPrompt_[i] = 0;
    genParticle_isPromptFinalState_[i] = 0;
    genParticle_isTauDecayProduct_[i] = 0;
    genParticle_isPromptTauDecayProduct_[i] = 0;
    genParticle_isDirectPromptTauDecayProductFinalState_[i] = 0;
    genParticle_isHardProcess_[i] = 0;
    genParticle_isLastCopy_[i] = 0;
    genParticle_isLastCopyBeforeFSR_[i] = 0;
    genParticle_isPromptDecayed_[i] = 0;
    genParticle_isDecayedLeptonHadron_[i] = 0;
    genParticle_fromHardProcessBeforeFSR_[i] = 0;
    genParticle_fromHardProcessDecayed_[i] = 0;
    genParticle_fromHardProcessFinalState_[i] = 0;
    genParticle_isMostlyLikePythia6Status3_[i] = 0;
  }

  // -- original trigger objects- - //
  vec_firedTrigger_.clear();
  vec_filterName_.clear();
  vec_HLTObj_pt_.clear();
  vec_HLTObj_eta_.clear();
  vec_HLTObj_phi_.clear();

  // -- HLT rerun objects -- //
  vec_myFiredTrigger_.clear();
  vec_myFilterName_.clear();
  vec_myHLTObj_pt_.clear();
  vec_myHLTObj_eta_.clear();
  vec_myHLTObj_phi_.clear();

  nMuon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    muon_pt_[i] = -999;
    muon_eta_[i] = -999;
    muon_phi_[i] = -999;
    muon_px_[i] = -999;
    muon_py_[i] = -999;
    muon_pz_[i] = -999;
    muon_dB_[i] = -999;
    muon_charge_[i] = -999;
    muon_isGLB_[i] = 0;
    muon_isSTA_[i] = 0;
    muon_isTRK_[i] = 0;
    muon_isPF_[i] = 0;
    muon_isTight_[i] = 0;
    muon_isMedium_[i] = 0;
    muon_isLoose_[i] = 0;
    muon_isHighPt_[i] = 0;
    muon_isHighPtNew_[i] = 0;
    muon_isSoft_[i] = 0;

    muon_isLooseTriggerMuon_[i] = 0;
    muon_isME0Muon_[i] = 0;
    muon_isGEMMuon_[i] = 0;
    muon_isRPCMuon_[i] = 0;
    muon_isGoodMuon_TMOneStationTight_[i] = 0;

    muon_iso03_sumPt_[i] = -999;
    muon_iso03_hadEt_[i] = -999;
    muon_iso03_emEt_[i] = -999;

    muon_PFIso03_charged_[i] = -999;
    muon_PFIso03_neutral_[i] = -999;
    muon_PFIso03_photon_[i] = -999;
    muon_PFIso03_sumPU_[i] = -999;

    muon_PFIso04_charged_[i] = -999;
    muon_PFIso04_neutral_[i] = -999;
    muon_PFIso04_photon_[i] = -999;
    muon_PFIso04_sumPU_[i] = -999;

    muon_PFCluster03_ECAL_[i] = -999;
    muon_PFCluster03_HCAL_[i] = -999;

    muon_PFCluster04_ECAL_[i] = -999;
    muon_PFCluster04_HCAL_[i] = -999;

    muon_normChi2_global_[i] = -999;
    muon_nTrackerHit_global_[i] = -999;
    muon_nTrackerLayer_global_[i] = -999;
    muon_nPixelHit_global_[i] = -999;
    muon_nMuonHit_global_[i] = -999;

    muon_normChi2_inner_[i] = -999;
    muon_nTrackerHit_inner_[i] = -999;
    muon_nTrackerLayer_inner_[i] = -999;
    muon_nPixelHit_inner_[i] = -999;

    muon_pt_tuneP_[i] = -999;
    muon_ptError_tuneP_[i] = -999;

    muon_dxyVTX_best_[i] = -999;
    muon_dzVTX_best_[i] = -999;

    muon_nMatchedStation_[i] = -999;
    muon_nMatchedRPCLayer_[i] = -999;
    muon_stationMask_[i] = -999;
    muon_expectedNnumberOfMatchedStations_[i] = -999;
  }

  // -- generic muon-like collections (replaces the old per-collection
  //    clearing blocks for L3Muon/L2Muon/TkMuon/iterL3OI/iterL3IOFromL2/
  //    iterL3FromL2 above)
  for (const auto& cfg : muonCollectionCfgs_) {
    muonTrack_Pt_[cfg.name].clear();
    muonTrack_Eta_[cfg.name].clear();
    muonTrack_Phi_[cfg.name].clear();
    muonTrack_Charge_[cfg.name].clear();
    muonTrack_TrackPt_[cfg.name].clear();
    muonTrack_InnerPt_[cfg.name].clear();  muonTrack_InnerEta_[cfg.name].clear();
    muonTrack_InnerPhi_[cfg.name].clear(); muonTrack_InnerCharge_[cfg.name].clear();
    muonTrack_OuterPt_[cfg.name].clear();  muonTrack_OuterEta_[cfg.name].clear();
    muonTrack_OuterPhi_[cfg.name].clear(); muonTrack_OuterCharge_[cfg.name].clear();
    muonTrack_GlobalPt_[cfg.name].clear(); muonTrack_GlobalEta_[cfg.name].clear();
    muonTrack_GlobalPhi_[cfg.name].clear();muonTrack_GlobalCharge_[cfg.name].clear();
    muonTrack_TrigMatchedFilterIdx_[cfg.name].clear();
    muonTrack_TrigDR_[cfg.name].clear();
    muonTrack_PathMatchedIdx_[cfg.name].clear();
    muonTrack_PathDR_[cfg.name].clear();
  }

  VThltIterL3MuonTrimmedPixelVertices->clear();
  VThltIterL3FromL1MuonTrimmedPixelVertices->clear();
}

void MuonHLTNtupler::Make_Branch()
{
  ntuple_->Branch("isRealData", &isRealData_, "isRealData/O"); // -- O: boolean -- //
  ntuple_->Branch("runNum",&runNum_,"runNum/I");
  ntuple_->Branch("lumiBlockNum",&lumiBlockNum_,"lumiBlockNum/I");
  ntuple_->Branch("eventNum",&eventNum_,"eventNum/l"); // -- unsigned long long -- //
  ntuple_->Branch("nVertex", &nVertex_, "nVertex/I");
  ntuple_->Branch("bunchID", &bunchID_, "bunchID/D");
  ntuple_->Branch("instLumi", &instLumi_, "instLumi/D");
  ntuple_->Branch("dataPU", &dataPU_, "dataPU/D");
  ntuple_->Branch("dataPURMS", &dataPURMS_, "dataPURMS/D");
  ntuple_->Branch("bunchLumi", &bunchLumi_, "bunchLumi/D");
  ntuple_->Branch("offlineInstLumi", &offlineInstLumi_, "offlineInstLumi/D");
  ntuple_->Branch("offlineDataPU", &offlineDataPU_, "offlineDataPU/D");
  ntuple_->Branch("offlineDataPURMS", &offlineDataPURMS_, "offlineDataPURMS/D");
  ntuple_->Branch("offlineBunchLumi", &offlineBunchLumi_, "offlineBunchLumi/D");
  ntuple_->Branch("truePU", &truePU_, "truePU/I");

  ntuple_->Branch("PU_pT_hats", &PU_pT_hats_);

  ntuple_->Branch("genEventWeight", &genEventWeight_, "genEventWeight/D");
  ntuple_->Branch("qScale", &qScale_, "qScale/D");

  ntuple_->Branch("nGenParticle", &nGenParticle_, "nGenParticle/I");
  ntuple_->Branch("genParticle_ID", &genParticle_ID_, "genParticle_ID[nGenParticle]/I");
  ntuple_->Branch("genParticle_status", &genParticle_status_, "genParticle_status[nGenParticle]/I");
  ntuple_->Branch("genParticle_mother", &genParticle_mother_, "genParticle_mother[nGenParticle]/I");
  ntuple_->Branch("genParticle_pt", &genParticle_pt_, "genParticle_pt[nGenParticle]/D");
  ntuple_->Branch("genParticle_eta", &genParticle_eta_, "genParticle_eta[nGenParticle]/D");
  ntuple_->Branch("genParticle_phi", &genParticle_phi_, "genParticle_phi[nGenParticle]/D");
  ntuple_->Branch("genParticle_px", &genParticle_px_, "genParticle_px[nGenParticle]/D");
  ntuple_->Branch("genParticle_py", &genParticle_py_, "genParticle_py[nGenParticle]/D");
  ntuple_->Branch("genParticle_pz", &genParticle_pz_, "genParticle_pz[nGenParticle]/D");
  ntuple_->Branch("genParticle_energy", &genParticle_energy_, "genParticle_energy[nGenParticle]/D");
  ntuple_->Branch("genParticle_charge", &genParticle_charge_, "genParticle_charge[nGenParticle]/D");
  ntuple_->Branch("genParticle_isPrompt", &genParticle_isPrompt_, "genParticle_isPrompt[nGenParticle]/I");
  ntuple_->Branch("genParticle_isPromptFinalState", &genParticle_isPromptFinalState_, "genParticle_isPromptFinalState[nGenParticle]/I");
  ntuple_->Branch("genParticle_isTauDecayProduct", &genParticle_isTauDecayProduct_, "genParticle_isTauDecayProduct[nGenParticle]/I");
  ntuple_->Branch("genParticle_isPromptTauDecayProduct", &genParticle_isPromptTauDecayProduct_, "genParticle_isPromptTauDecayProduct[nGenParticle]/I");
  ntuple_->Branch("genParticle_isDirectPromptTauDecayProductFinalState", &genParticle_isDirectPromptTauDecayProductFinalState_, "genParticle_isDirectPromptTauDecayProductFinalState[nGenParticle]/I");
  ntuple_->Branch("genParticle_isHardProcess", &genParticle_isHardProcess_, "genParticle_isHardProcess[nGenParticle]/I");
  ntuple_->Branch("genParticle_isLastCopy", &genParticle_isLastCopy_, "genParticle_isLastCopy[nGenParticle]/I");
  ntuple_->Branch("genParticle_isLastCopyBeforeFSR", &genParticle_isLastCopyBeforeFSR_, "genParticle_isLastCopyBeforeFSR[nGenParticle]/I");
  ntuple_->Branch("genParticle_isPromptDecayed", &genParticle_isPromptDecayed_, "genParticle_isPromptDecayed[nGenParticle]/I");
  ntuple_->Branch("genParticle_isDecayedLeptonHadron", &genParticle_isDecayedLeptonHadron_, "genParticle_isDecayedLeptonHadron[nGenParticle]/I");
  ntuple_->Branch("genParticle_fromHardProcessBeforeFSR", &genParticle_fromHardProcessBeforeFSR_, "genParticle_fromHardProcessBeforeFSR[nGenParticle]/I");
  ntuple_->Branch("genParticle_fromHardProcessDecayed", &genParticle_fromHardProcessDecayed_, "genParticle_fromHardProcessDecayed[nGenParticle]/I");
  ntuple_->Branch("genParticle_fromHardProcessFinalState", &genParticle_fromHardProcessFinalState_, "genParticle_fromHardProcessFinalState[nGenParticle]/I");
  ntuple_->Branch("genParticle_isMostlyLikePythia6Status3", &genParticle_isMostlyLikePythia6Status3_, "genParticle_isMostlyLikePythia6Status3[nGenParticle]/I");

  ntuple_->Branch("vec_firedTrigger", &vec_firedTrigger_);
  ntuple_->Branch("vec_filterName", &vec_filterName_);
  ntuple_->Branch("vec_HLTObj_pt", &vec_HLTObj_pt_);
  ntuple_->Branch("vec_HLTObj_eta", &vec_HLTObj_eta_);
  ntuple_->Branch("vec_HLTObj_phi", &vec_HLTObj_phi_);

  ntuple_->Branch("vec_myFiredTrigger", &vec_myFiredTrigger_);
  ntuple_->Branch("vec_myFilterName", &vec_myFilterName_);
  ntuple_->Branch("vec_myHLTObj_pt", &vec_myHLTObj_pt_);
  ntuple_->Branch("vec_myHLTObj_eta", &vec_myHLTObj_eta_);
  ntuple_->Branch("vec_myHLTObj_phi", &vec_myHLTObj_phi_);

  ntuple_->Branch("nMuon", &nMuon_, "nMuon/I");

  ntuple_->Branch("muon_pt", &muon_pt_, "muon_pt[nMuon]/D");
  ntuple_->Branch("muon_eta", &muon_eta_, "muon_eta[nMuon]/D");
  ntuple_->Branch("muon_phi", &muon_phi_, "muon_phi[nMuon]/D");
  ntuple_->Branch("muon_px", &muon_px_, "muon_px[nMuon]/D");
  ntuple_->Branch("muon_py", &muon_py_, "muon_py[nMuon]/D");
  ntuple_->Branch("muon_pz", &muon_pz_, "muon_pz[nMuon]/D");
  ntuple_->Branch("muon_dB", &muon_dB_, "muon_dB[nMuon]/D");
  ntuple_->Branch("muon_charge", &muon_charge_, "muon_charge[nMuon]/D");
  ntuple_->Branch("muon_isGLB", &muon_isGLB_, "muon_isGLB[nMuon]/I");
  ntuple_->Branch("muon_isSTA", &muon_isSTA_, "muon_isSTA[nMuon]/I");
  ntuple_->Branch("muon_isTRK", &muon_isTRK_, "muon_isTRK[nMuon]/I");
  ntuple_->Branch("muon_isPF", &muon_isPF_, "muon_isPF[nMuon]/I");
  ntuple_->Branch("muon_isTight", &muon_isTight_, "muon_isTight[nMuon]/I");
  ntuple_->Branch("muon_isMedium", &muon_isMedium_, "muon_isMedium[nMuon]/I");
  ntuple_->Branch("muon_isLoose", &muon_isLoose_, "muon_isLoose[nMuon]/I");
  ntuple_->Branch("muon_isHighPt", &muon_isHighPt_, "muon_isHighPt[nMuon]/I");
  ntuple_->Branch("muon_isHighPtNew", &muon_isHighPtNew_, "muon_isHighPtNew[nMuon]/I");
  ntuple_->Branch("muon_isSoft", &muon_isSoft_, "muon_isSoft[nMuon]/I");

  ntuple_->Branch("muon_isLooseTriggerMuon", &muon_isLooseTriggerMuon_, "muon_isLooseTriggerMuon[nMuon]/I");
  ntuple_->Branch("muon_isME0Muon", &muon_isME0Muon_, "muon_isME0Muon[nMuon]/I");
  ntuple_->Branch("muon_isGEMMuon", &muon_isGEMMuon_, "muon_isGEMMuon[nMuon]/I");
  ntuple_->Branch("muon_isRPCMuon", &muon_isRPCMuon_, "muon_isRPCMuon[nMuon]/I");
  ntuple_->Branch("muon_isGoodMuon_TMOneStationTight", &muon_isGoodMuon_TMOneStationTight_, "muon_isGoodMuon_TMOneStationTight[nMuon]/I");

  ntuple_->Branch("muon_iso03_sumPt", &muon_iso03_sumPt_, "muon_iso03_sumPt[nMuon]/D");
  ntuple_->Branch("muon_iso03_hadEt", &muon_iso03_hadEt_, "muon_iso03_hadEt[nMuon]/D");
  ntuple_->Branch("muon_iso03_emEt", &muon_iso03_emEt_, "muon_iso03_emEt[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_charged", &muon_PFIso03_charged_, "muon_PFIso03_charged[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_neutral", &muon_PFIso03_neutral_, "muon_PFIso03_neutral[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_photon", &muon_PFIso03_photon_, "muon_PFIso03_photon[nMuon]/D");
  ntuple_->Branch("muon_PFIso03_sumPU", &muon_PFIso03_sumPU_, "muon_PFIso03_sumPU[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_charged", &muon_PFIso04_charged_, "muon_PFIso04_charged[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_neutral", &muon_PFIso04_neutral_, "muon_PFIso04_neutral[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_photon", &muon_PFIso04_photon_, "muon_PFIso04_photon[nMuon]/D");
  ntuple_->Branch("muon_PFIso04_sumPU", &muon_PFIso04_sumPU_, "muon_PFIso04_sumPU[nMuon]/D");

  ntuple_->Branch("muon_PFCluster03_ECAL", &muon_PFCluster03_ECAL_, "muon_PFCluster03_ECAL[nMuon]/D");
  ntuple_->Branch("muon_PFCluster03_HCAL", &muon_PFCluster03_HCAL_, "muon_PFCluster03_HCAL[nMuon]/D");
  ntuple_->Branch("muon_PFCluster04_ECAL", &muon_PFCluster04_ECAL_, "muon_PFCluster04_ECAL[nMuon]/D");
  ntuple_->Branch("muon_PFCluster04_HCAL", &muon_PFCluster04_HCAL_, "muon_PFCluster04_HCAL[nMuon]/D");
  ntuple_->Branch("muon_normChi2_global", &muon_normChi2_global_, "muon_normChi2_global[nMuon]/D");
  ntuple_->Branch("muon_nTrackerHit_global", &muon_nTrackerHit_global_, "muon_nTrackerHit_global[nMuon]/I");
  ntuple_->Branch("muon_nTrackerLayer_global", &muon_nTrackerLayer_global_, "muon_nTrackerLayer_global[nMuon]/I");
  ntuple_->Branch("muon_nPixelHit_global", &muon_nPixelHit_global_, "muon_nPixelHit_global[nMuon]/I");
  ntuple_->Branch("muon_nMuonHit_global", &muon_nMuonHit_global_, "muon_nMuonHit_global[nMuon]/I");
  ntuple_->Branch("muon_normChi2_inner", &muon_normChi2_inner_, "muon_normChi2_inner[nMuon]/D");
  ntuple_->Branch("muon_nTrackerHit_inner", &muon_nTrackerHit_inner_, "muon_nTrackerHit_inner[nMuon]/I");
  ntuple_->Branch("muon_nTrackerLayer_inner", &muon_nTrackerLayer_inner_, "muon_nTrackerLayer_inner[nMuon]/I");
  ntuple_->Branch("muon_nPixelHit_inner", &muon_nPixelHit_inner_, "muon_nPixelHit_inner[nMuon]/I");
  ntuple_->Branch("muon_pt_tuneP", &muon_pt_tuneP_, "muon_pt_tuneP[nMuon]/D");
  ntuple_->Branch("muon_ptError_tuneP", &muon_ptError_tuneP_, "muon_ptError_tuneP[nMuon]/D");
  ntuple_->Branch("muon_dxyVTX_best", &muon_dxyVTX_best_, "muon_dxyVTX_best[nMuon]/D");
  ntuple_->Branch("muon_dzVTX_best", &muon_dzVTX_best_, "muon_dzVTX_best[nMuon]/D");
  ntuple_->Branch("muon_nMatchedStation", &muon_nMatchedStation_, "muon_nMatchedStation[nMuon]/I");
  ntuple_->Branch("muon_nMatchedRPCLayer", &muon_nMatchedRPCLayer_, "muon_nMatchedRPCLayer[nMuon]/I");
  ntuple_->Branch("muon_stationMask", &muon_stationMask_, "muon_stationMask[nMuon]/I");
  ntuple_->Branch("muon_expectedNnumberOfMatchedStations", &muon_expectedNnumberOfMatchedStations_, "muon_expectedNnumberOfMatchedStations[nMuon]/I");

  // -- generic muon-like collections (default: L3Muon, L2Muon, TkMuon,
  //    iterL3OI, iterL3IOFromL2, iterL3FromL2 -- see ntupler_cfi.py). One
  //    set of branches per configured collection name; the branch *schema*
  //    stays fixed even as triggerFilters_ changes, since matches are
  //    stored as filter *indices* (decode via the "configuredTriggerFilters"
  //    branch written once per event below) rather than one branch per
  //    filter.
  for (const auto& cfg : muonCollectionCfgs_) {
    ntuple_->Branch((cfg.name + "_pt").c_str(),     &muonTrack_Pt_[cfg.name]);
    ntuple_->Branch((cfg.name + "_eta").c_str(),    &muonTrack_Eta_[cfg.name]);
    ntuple_->Branch((cfg.name + "_phi").c_str(),    &muonTrack_Phi_[cfg.name]);
    ntuple_->Branch((cfg.name + "_charge").c_str(), &muonTrack_Charge_[cfg.name]);
    ntuple_->Branch((cfg.name + "_trkPt").c_str(),  &muonTrack_TrackPt_[cfg.name]);

    // -- only meaningfully filled for type == "MuonTrackLinks"; -99 for
    //    every other collection, preserving the old inner/outer/global
    //    branch triplets that iterL3OI/iterL3IOFromL2/iterL3FromL2 used to
    //    have as dedicated arrays.
    ntuple_->Branch((cfg.name + "_inner_pt").c_str(),     &muonTrack_InnerPt_[cfg.name]);
    ntuple_->Branch((cfg.name + "_inner_eta").c_str(),    &muonTrack_InnerEta_[cfg.name]);
    ntuple_->Branch((cfg.name + "_inner_phi").c_str(),    &muonTrack_InnerPhi_[cfg.name]);
    ntuple_->Branch((cfg.name + "_inner_charge").c_str(), &muonTrack_InnerCharge_[cfg.name]);
    ntuple_->Branch((cfg.name + "_outer_pt").c_str(),     &muonTrack_OuterPt_[cfg.name]);
    ntuple_->Branch((cfg.name + "_outer_eta").c_str(),    &muonTrack_OuterEta_[cfg.name]);
    ntuple_->Branch((cfg.name + "_outer_phi").c_str(),    &muonTrack_OuterPhi_[cfg.name]);
    ntuple_->Branch((cfg.name + "_outer_charge").c_str(), &muonTrack_OuterCharge_[cfg.name]);
    ntuple_->Branch((cfg.name + "_global_pt").c_str(),     &muonTrack_GlobalPt_[cfg.name]);
    ntuple_->Branch((cfg.name + "_global_eta").c_str(),    &muonTrack_GlobalEta_[cfg.name]);
    ntuple_->Branch((cfg.name + "_global_phi").c_str(),    &muonTrack_GlobalPhi_[cfg.name]);
    ntuple_->Branch((cfg.name + "_global_charge").c_str(), &muonTrack_GlobalCharge_[cfg.name]);

    ntuple_->Branch((cfg.name + "_trigMatchedFilterIdx").c_str(), &muonTrack_TrigMatchedFilterIdx_[cfg.name]);
    ntuple_->Branch((cfg.name + "_trigDR").c_str(),               &muonTrack_TrigDR_[cfg.name]);

    // -- "did this object fire path P" (matched against P's last filter,
    //    see MuonHLTTriggerMatching::lastFilterOfPath), as opposed to
    //    "_trigMatchedFilterIdx" above which is per individually-configured
    //    filter. Indices refer to "configuredTriggerPaths" below.
    ntuple_->Branch((cfg.name + "_pathMatchedIdx").c_str(), &muonTrack_PathMatchedIdx_[cfg.name]);
    ntuple_->Branch((cfg.name + "_pathDR").c_str(),         &muonTrack_PathDR_[cfg.name]);
  }
  // -- lookup table for the filter/path indices above; written once per
  //    event (not once per file) so it survives merging ntuples across menus.
  ntuple_->Branch("configuredTriggerFilters", &triggerFiltersCfg_);
  ntuple_->Branch("configuredTriggerPaths",   &resolvedTriggerPaths_);
  ntuple_->Branch("configuredMuonCollections", &muonCollectionNamesForBranch_);

  VThltIterL3MuonTrimmedPixelVertices->setBranch(ntuple_,"hltIterL3MuonTrimmedPixelVertices");
  VThltIterL3FromL1MuonTrimmedPixelVertices->setBranch(ntuple_,"hltIterL3FromL1MuonTrimmedPixelVertices");

}

void MuonHLTNtupler::Fill_Muon(const edm::Event &iEvent)
{
  // edm::Handle<std::vector<reco::Muon> > h_offlineMuon;
  edm::Handle< edm::View<reco::Muon> > h_offlineMuon;
  if( iEvent.getByToken(t_offlineMuon_, h_offlineMuon) ) // -- only when the dataset has offline muon collection (e.g. AOD) -- //
  {
    edm::Handle<reco::VertexCollection> h_offlineVertex;
    bool isVertex = iEvent.getByToken(t_offlineVertex_, h_offlineVertex);
    // const reco::Vertex & pv = h_offlineVertex->at(0);

    int _nMuon = 0;
    for(auto mu=h_offlineMuon->begin(); mu!=h_offlineMuon->end(); ++mu)
    {
      muon_pt_[_nMuon]  = mu->pt();
      muon_eta_[_nMuon] = mu->eta();
      muon_phi_[_nMuon] = mu->phi();
      muon_px_[_nMuon]  = mu->px();
      muon_py_[_nMuon]  = mu->py();
      muon_pz_[_nMuon]  = mu->pz();
      // muon_dB_[_nMuon] = mu->dB(); // -- dB is only availabe in pat::Muon -- //
      muon_charge_[_nMuon] = mu->charge();

      if( mu->isGlobalMuon() ) muon_isGLB_[_nMuon] = 1;
      if( mu->isStandAloneMuon() ) muon_isSTA_[_nMuon] = 1;
      if( mu->isTrackerMuon() ) muon_isTRK_[_nMuon] = 1;
      if( mu->isPFMuon() ) muon_isPF_[_nMuon] = 1;

      // -- defintion of ID functions: http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_4_0/doc/html/da/d18/namespacemuon.html#ac122b2516e5711ce206256d7945473d2 -- //
      if( muon::isMediumMuon( (*mu) ) )     muon_isMedium_[_nMuon] = 1;
      if( muon::isLooseMuon( (*mu) ) )      muon_isLoose_[_nMuon] = 1;
      if(isVertex) {
        const reco::Vertex & pv = h_offlineVertex->at(0);
        if( muon::isTightMuon( (*mu), pv ) )  muon_isTight_[_nMuon] = 1;
        if( muon::isHighPtMuon( (*mu), pv ) ) muon_isHighPt_[_nMuon] = 1;
        if( isNewHighPtMuon( (*mu), pv ) )    muon_isHighPtNew_[_nMuon] = 1;
      }

      if( muon::isLooseTriggerMuon( (*mu) ) )                   muon_isLooseTriggerMuon_[_nMuon] = 1;
      if( mu->isME0Muon() )                                     muon_isME0Muon_[_nMuon] = 1;
      if( mu->isGEMMuon() )                                     muon_isGEMMuon_[_nMuon] = 1;
      if( mu->isRPCMuon() )                                     muon_isRPCMuon_[_nMuon] = 1;
      if( muon::isGoodMuon( (*mu), muon::TMOneStationTight ) )  muon_isGoodMuon_TMOneStationTight_[_nMuon] = 1;

      // -- bool muon::isSoftMuon(const reco::Muon& muon, const reco::Vertex& vtx, bool run2016_hip_mitigation)
      // -- it is different under CMSSW_8_0_29: bool muon::isSoftMuon(const reco::Muon& muon, const reco::Vertex& vtx)
      // -- Remove this part to avoid compile error (and soft muon would not be used for now) - need to be fixed at some point
      // if( muon::isSoftMuon( (*mu), pv, 0) ) muon_isSoft_[_nMuon] = 1;

      muon_iso03_sumPt_[_nMuon] = mu->isolationR03().sumPt;
      muon_iso03_hadEt_[_nMuon] = mu->isolationR03().hadEt;
      muon_iso03_emEt_[_nMuon]  = mu->isolationR03().emEt;

      muon_PFIso03_charged_[_nMuon] = mu->pfIsolationR03().sumChargedHadronPt;
      muon_PFIso03_neutral_[_nMuon] = mu->pfIsolationR03().sumNeutralHadronEt;
      muon_PFIso03_photon_[_nMuon]  = mu->pfIsolationR03().sumPhotonEt;
      muon_PFIso03_sumPU_[_nMuon]   = mu->pfIsolationR03().sumPUPt;

      muon_PFIso04_charged_[_nMuon] = mu->pfIsolationR04().sumChargedHadronPt;
      muon_PFIso04_neutral_[_nMuon] = mu->pfIsolationR04().sumNeutralHadronEt;
      muon_PFIso04_photon_[_nMuon]  = mu->pfIsolationR04().sumPhotonEt;
      muon_PFIso04_sumPU_[_nMuon]   = mu->pfIsolationR04().sumPUPt;

      // reco::MuonRef muRef = reco::MuonRef(h_offlineMuon, _nMuon);

      reco::TrackRef globalTrk = mu->globalTrack();
      if( globalTrk.isNonnull() )
      {
        muon_normChi2_global_[_nMuon] = globalTrk->normalizedChi2();

        const reco::HitPattern & globalTrkHit = globalTrk->hitPattern();
        muon_nTrackerHit_global_[_nMuon]   = globalTrkHit.numberOfValidTrackerHits();
        muon_nTrackerLayer_global_[_nMuon] = globalTrkHit.trackerLayersWithMeasurement();
        muon_nPixelHit_global_[_nMuon]     = globalTrkHit.numberOfValidPixelHits();
        muon_nMuonHit_global_[_nMuon]      = globalTrkHit.numberOfValidMuonHits();
      }

      reco::TrackRef innerTrk = mu->innerTrack();
      if( innerTrk.isNonnull() )
      {
        muon_normChi2_inner_[_nMuon] = innerTrk->normalizedChi2();

        const reco::HitPattern & innerTrkHit = innerTrk->hitPattern();
        muon_nTrackerHit_inner_[_nMuon]   = innerTrkHit.numberOfValidTrackerHits();
        muon_nTrackerLayer_inner_[_nMuon] = innerTrkHit.trackerLayersWithMeasurement();
        muon_nPixelHit_inner_[_nMuon]     = innerTrkHit.numberOfValidPixelHits();
      }

      reco::TrackRef tunePTrk = mu->tunePMuonBestTrack();
      if( tunePTrk.isNonnull() )
      {
        muon_pt_tuneP_[_nMuon]      = tunePTrk->pt();
        muon_ptError_tuneP_[_nMuon] = tunePTrk->ptError();
      }

      if(isVertex) {
        const reco::Vertex & pv = h_offlineVertex->at(0);
        muon_dxyVTX_best_[_nMuon] = mu->muonBestTrack()->dxy( pv.position() );
        muon_dzVTX_best_[_nMuon]  = mu->muonBestTrack()->dz( pv.position() );
      }

      muon_nMatchedStation_[_nMuon] = mu->numberOfMatchedStations();
      muon_nMatchedRPCLayer_[_nMuon] = mu->numberOfMatchedRPCLayers();
      muon_stationMask_[_nMuon] = mu->stationMask();
      muon_expectedNnumberOfMatchedStations_[_nMuon] = mu->expectedNnumberOfMatchedStations();

      _nMuon++;
    }

    nMuon_ = _nMuon;
  }
}

void MuonHLTNtupler::Fill_HLT(const edm::Event &iEvent, bool isMYHLT)
{
  edm::Handle<edm::TriggerResults>  h_triggerResults;
  edm::Handle<trigger::TriggerEvent> h_triggerEvent;

  if( isMYHLT )
  {
    iEvent.getByToken(t_myTriggerResults_, h_triggerResults);
    iEvent.getByToken(t_myTriggerEvent_,   h_triggerEvent);
  }
  else
  {
    iEvent.getByToken(t_triggerResults_, h_triggerResults);
    iEvent.getByToken(t_triggerEvent_,   h_triggerEvent);
  }

  edm::TriggerNames triggerNames = iEvent.triggerNames(*h_triggerResults);

  for(unsigned int itrig=0; itrig<triggerNames.size(); ++itrig)
  {
    LogDebug("triggers") << triggerNames.triggerName(itrig);

    if( h_triggerResults->accept(itrig) )
    {
      std::string pathName = triggerNames.triggerName(itrig);
      // if( SavedTriggerCondition(pathName) )
      if( true )
      {
        if( isMYHLT ) vec_myFiredTrigger_.push_back( pathName );
        else          vec_firedTrigger_.push_back( pathName );
      }
    } // -- end of if fired -- //

  } // -- end of iteration over all trigger names -- //

  const trigger::size_type nFilter(h_triggerEvent->sizeFilters());
  for( trigger::size_type i_filter=0; i_filter<nFilter; i_filter++)
  {
    std::string filterName = h_triggerEvent->filterTag(i_filter).encode();

    // if( SavedFilterCondition(filterName) )
    if( true )
    {
      trigger::Keys objectKeys = h_triggerEvent->filterKeys(i_filter);
      const trigger::TriggerObjectCollection& triggerObjects(h_triggerEvent->getObjects());

      for( trigger::size_type i_key=0; i_key<objectKeys.size(); i_key++)
      {
        trigger::size_type objKey = objectKeys.at(i_key);
        const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);

        if( isMYHLT )
        {
          vec_myFilterName_.push_back( filterName );
          vec_myHLTObj_pt_.push_back( triggerObj.pt() );
          vec_myHLTObj_eta_.push_back( triggerObj.eta() );
          vec_myHLTObj_phi_.push_back( triggerObj.phi() );
        }
        else
        {
          vec_filterName_.push_back( filterName );
          vec_HLTObj_pt_.push_back( triggerObj.pt() );
          vec_HLTObj_eta_.push_back( triggerObj.eta() );
          vec_HLTObj_phi_.push_back( triggerObj.phi() );
        }
      }
    } // -- end of if( muon filters )-- //
  } // -- end of filter iteration -- //
}

bool MuonHLTNtupler::SavedTriggerCondition( std::string& pathName )
{
  bool flag = false;

  // -- muon triggers
  if( pathName.find("HLT_IsoMu")    != std::string::npos ||
      pathName.find("HLT_Mu")       != std::string::npos ||
      pathName.find("HLT_OldMu")    != std::string::npos ||
      pathName.find("HLT_TkMu")     != std::string::npos ||
      pathName.find("HLT_IsoTkMu")  != std::string::npos ||
      pathName.find("HLT_DoubleMu") != std::string::npos ||
      pathName.find("HLT_Mu8_T")    != std::string::npos ) flag = true;

  return flag;
}

bool MuonHLTNtupler::SavedFilterCondition( std::string& filterName )
{
  bool flag = false;

  // -- muon filters
  if( (filterName.find("sMu") != std::string::npos || filterName.find("SingleMu") != std::string::npos) &&
       filterName.find("Tau")      == std::string::npos &&
       filterName.find("EG")       == std::string::npos &&
       filterName.find("MultiFit") == std::string::npos ) flag = true;

  return flag;
}

bool MuonHLTNtupler::isMuonCollectionValid(const std::string& name) const
{
  for (size_t i = 0; i < muonCollectionCfgs_.size(); ++i) {
    if (muonCollectionCfgs_[i].name == name) return muonCollectionAdapters_[i]->isValid();
  }
  return false;
}

void MuonHLTNtupler::Fill_GenericMuonCollections(const edm::Event &iEvent)
{
  // -- Every collection configured in the python
  //    "muonCollections" VPSet is filled the same way here, and every one
  //    of them additionally gets real per-object trigger matching against
  //    the configured "triggerFilters" -- both were previously handled by
  //    hardcoded, per-collection blocks with no dR matching at all.

  edm::Handle<trigger::TriggerEvent> h_triggerEvent;
  bool haveTriggerEvent = iEvent.getByToken(t_myTriggerEvent_, h_triggerEvent);

  // Get filter and path objects to later match them to the muon collection
  std::vector<std::vector<MuonHLTTriggerMatching::FilterObject>> filterObjs;
  std::vector<std::vector<MuonHLTTriggerMatching::FilterObject>> pathFilterObjs;
  if (haveTriggerEvent) {
    filterObjs.reserve(triggerFiltersCfg_.size());
    for (const auto& filt : triggerFiltersCfg_) {
      filterObjs.push_back(MuonHLTTriggerMatching::getFilterObjects(*h_triggerEvent, filt));
    }
    pathFilterObjs.reserve(resolvedPathLastFilter_.size());
    for (const auto& lastFilt : resolvedPathLastFilter_) {
      pathFilterObjs.push_back(lastFilt.empty()
                                    ? std::vector<MuonHLTTriggerMatching::FilterObject>()
                                    : MuonHLTTriggerMatching::getFilterObjects(*h_triggerEvent, lastFilt));
    }
  }

  // Loop over the collections initialized in the configuration
  // and store information for each muon track type/step
  //
  // Allows: [RecoChargedCandidate, TrackView, RecoMuon, MuonTrackLinks, L1TkMuon]
  //
  for (size_t iColl = 0; iColl < muonCollectionAdapters_.size(); ++iColl) {
    const auto& cfg = muonCollectionCfgs_[iColl];
    auto& adapter = muonCollectionAdapters_[iColl];

    adapter->getByEvent(iEvent);
    if (!adapter->isValid()) {
      if (!warnedCollectionInvalid_[iColl]) {
        edm::LogWarning("MuonHLTNtupler")
            << "muonCollections entry '" << cfg.name << "' (type=\"" << cfg.type
            << "\", tag=" << cfg.tag.encode() << ") was not found/valid in this event.\n"
            << "If this collection is genuinely absent in some events/runs, that's fine and "
            << "this warning is harmless. But if it should always be present, the most likely "
            << "cause is that 'type' doesn't match this InputTag's actual EDM product type -- "
            << "a mismatch makes getByToken silently return an invalid handle instead of "
            << "throwing, every event, forever. Double check the producer that makes \""
            << cfg.tag.label() << "\" and confirm its output type matches the adapter you asked "
            << "for (RecoChargedCandidate/TrackView/RecoMuon/MuonTrackLinks/L1TkMuon in "
            << "GenericMuonCollections.h). This warning will not repeat for this collection.";
        warnedCollectionInvalid_[iColl] = true;
      }
      continue;
    }

    for (const auto& mu : adapter->extract()) {
      muonTrack_Pt_[cfg.name].push_back(mu.pt);
      muonTrack_Eta_[cfg.name].push_back(mu.eta);
      muonTrack_Phi_[cfg.name].push_back(mu.phi);
      muonTrack_Charge_[cfg.name].push_back(mu.charge);
      muonTrack_TrackPt_[cfg.name].push_back(mu.trackPt);

      muonTrack_InnerPt_[cfg.name].push_back(mu.inner_pt);
      muonTrack_InnerEta_[cfg.name].push_back(mu.inner_eta);
      muonTrack_InnerPhi_[cfg.name].push_back(mu.inner_phi);
      muonTrack_InnerCharge_[cfg.name].push_back(mu.inner_charge);
      muonTrack_OuterPt_[cfg.name].push_back(mu.outer_pt);
      muonTrack_OuterEta_[cfg.name].push_back(mu.outer_eta);
      muonTrack_OuterPhi_[cfg.name].push_back(mu.outer_phi);
      muonTrack_OuterCharge_[cfg.name].push_back(mu.outer_charge);
      muonTrack_GlobalPt_[cfg.name].push_back(mu.global_pt);
      muonTrack_GlobalEta_[cfg.name].push_back(mu.global_eta);
      muonTrack_GlobalPhi_[cfg.name].push_back(mu.global_phi);
      muonTrack_GlobalCharge_[cfg.name].push_back(mu.global_charge);

      // -- real per-object dR matching against every configured filter,
      //    replacing the old "dump everything, match nothing" behavior.
      //    Build one inner vector for THIS muon (index = current size of
      //    muonTrack_Pt_[cfg.name], i.e. this muon's own index in the collection)
      //    listing every configured filter it matched, decoded via
      //    "configuredTriggerFilters". Pushed even when empty, so indexing
      //    genTrigMatchedFilterIdx_[name][i] always lines up with
      //    muonTrack_Pt_[name][i]/genEta_[name][i]/etc.
      std::vector<int> matchedFilterIdx;
      std::vector<float> matchedDR;
      if (haveTriggerEvent) {
        for (size_t iFilt = 0; iFilt < triggerFiltersCfg_.size(); ++iFilt) {
          float dr, ptAtFilter;
          if (MuonHLTTriggerMatching::matchToFilter(mu.eta, mu.phi, filterObjs[iFilt], maxDR_, dr, ptAtFilter)) {
            matchedFilterIdx.push_back(static_cast<int>(iFilt));
            matchedDR.push_back(dr);
          }
        }
      }
      muonTrack_TrigMatchedFilterIdx_[cfg.name].push_back(matchedFilterIdx);
      muonTrack_TrigDR_[cfg.name].push_back(matchedDR);

      // -- same thing, but "did this object fire path P" (matched against
      //    P's last filter) rather than per individually-configured
      //    filter -- see resolvedPathLastFilter_/lastFilterOfPath().
      std::vector<int> matchedPathIdx;
      std::vector<float> matchedPathDR;
      if (haveTriggerEvent) {
        for (size_t iPath = 0; iPath < resolvedPathLastFilter_.size(); ++iPath) {
          float dr, ptAtFilter;
          if (MuonHLTTriggerMatching::matchToFilter(mu.eta, mu.phi, pathFilterObjs[iPath], maxDR_, dr, ptAtFilter)) {
            matchedPathIdx.push_back(static_cast<int>(iPath));
            matchedPathDR.push_back(dr);
          }
        }
      }
      muonTrack_PathMatchedIdx_[cfg.name].push_back(matchedPathIdx);
      muonTrack_PathDR_[cfg.name].push_back(matchedPathDR);
    }
  }
}

void MuonHLTNtupler::Fill_GenParticle(const edm::Event &iEvent)
{
  // -- Gen-weight info -- //
  edm::Handle<GenEventInfoProduct> h_genEventInfo;
  iEvent.getByToken(t_genEventInfo_, h_genEventInfo);
  genEventWeight_ = h_genEventInfo->weight();
  qScale_ = h_genEventInfo->qScale();

  // -- Gen-particle info -- //
  edm::Handle<reco::GenParticleCollection> h_genParticle;
  iEvent.getByToken(t_genParticle_, h_genParticle);

  int _nGenParticle = 0;
  for( size_t i=0; i< h_genParticle->size(); ++i)
  {
    const reco::GenParticle &parCand = (*h_genParticle)[i];

    if( abs(parCand.pdgId()) == 13 || parCand.isHardProcess() ) // -- only muons -- //
    {
      genParticle_ID_[_nGenParticle]     = parCand.pdgId();
      genParticle_status_[_nGenParticle] = parCand.status();
      // genParticle_mother_[_nGenParticle] = parCand.mother(0)->pdgId();

      genParticle_pt_[_nGenParticle]  = parCand.pt();
      genParticle_eta_[_nGenParticle] = parCand.eta();
      genParticle_phi_[_nGenParticle] = parCand.phi();
      genParticle_px_[_nGenParticle]  = parCand.px();
      genParticle_py_[_nGenParticle]  = parCand.py();
      genParticle_pz_[_nGenParticle]  = parCand.pz();
      genParticle_energy_[_nGenParticle] = parCand.energy();
      genParticle_charge_[_nGenParticle] = parCand.charge();

      if( parCand.statusFlags().isPrompt() )                genParticle_isPrompt_[_nGenParticle] = 1;
      if( parCand.statusFlags().isTauDecayProduct() )       genParticle_isTauDecayProduct_[_nGenParticle] = 1;
      if( parCand.statusFlags().isPromptTauDecayProduct() ) genParticle_isPromptTauDecayProduct_[_nGenParticle] = 1;
      if( parCand.statusFlags().isDecayedLeptonHadron() )   genParticle_isDecayedLeptonHadron_[_nGenParticle] = 1;

      if( parCand.isPromptFinalState() ) genParticle_isPromptFinalState_[_nGenParticle] = 1;
      if( parCand.isDirectPromptTauDecayProductFinalState() ) genParticle_isDirectPromptTauDecayProductFinalState_[_nGenParticle] = 1;
      if( parCand.isHardProcess() ) genParticle_isHardProcess_[_nGenParticle] = 1;
      if( parCand.isLastCopy() ) genParticle_isLastCopy_[_nGenParticle] = 1;
      if( parCand.isLastCopyBeforeFSR() ) genParticle_isLastCopyBeforeFSR_[_nGenParticle] = 1;

      if( parCand.isPromptDecayed() )           genParticle_isPromptDecayed_[_nGenParticle] = 1;
      if( parCand.fromHardProcessBeforeFSR() )  genParticle_fromHardProcessBeforeFSR_[_nGenParticle] = 1;
      if( parCand.fromHardProcessDecayed() )    genParticle_fromHardProcessDecayed_[_nGenParticle] = 1;
      if( parCand.fromHardProcessFinalState() ) genParticle_fromHardProcessFinalState_[_nGenParticle] = 1;
      // if( parCand.isMostlyLikePythia6Status3() ) this->genParticle_isMostlyLikePythia6Status3[_nGenParticle] = 1;

      _nGenParticle++;
    }
  }
  nGenParticle_ = _nGenParticle;
}


// -- reference: https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonReco/src/MuonSelectors.cc#L910-L938
bool MuonHLTNtupler::isNewHighPtMuon(const reco::Muon& muon, const reco::Vertex& vtx){
  if(!muon.isGlobalMuon()) return false;

  bool muValHits = ( muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0 ||
                     muon.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits()>0 );

  bool muMatchedSt = muon.numberOfMatchedStations()>1;
  if(!muMatchedSt) {
    if( muon.isTrackerMuon() && muon.numberOfMatchedStations()==1 ) {
      if( muon.expectedNnumberOfMatchedStations()<2 ||
          !(muon.stationMask()==1 || muon.stationMask()==16) ||
          muon.numberOfMatchedRPCLayers()>2
        )
        muMatchedSt = true;
    }
  }

  bool muID = muValHits && muMatchedSt;

  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

  bool momQuality = muon.tunePMuonBestTrack()->ptError()/muon.tunePMuonBestTrack()->pt() < 0.3;

  bool ip = fabs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.innerTrack()->dz(vtx.position())) < 0.5;

  return muID && hits && momQuality && ip;
}

void MuonHLTNtupler::endJob() {

  delete VThltIterL3MuonTrimmedPixelVertices;
  delete VThltIterL3FromL1MuonTrimmedPixelVertices;
}

void MuonHLTNtupler::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {
  // -- resolve unversioned python path names ("HLT_Mu50_v") against this
  //    run's actual menu ("HLT_Mu50_v12") once per run, so triggerPathsCfg_
  //    never needs to track trigger versions. See MuonHLTTriggerMatching.h.

  std::cout << "Start beginRun " << std::endl;
  
  bool changed = true;
  resolvedPathLastFilter_.assign(triggerPathsCfg_.size(), "");
  if (hltConfig_.init(iRun, iSetup, "MYHLT", changed)) {
    resolvedTriggerPaths_ = MuonHLTTriggerMatching::resolvePathNames(hltConfig_, triggerPathsCfg_);
    
    // -- an object is said to have "fired" a path if it matches that
    //    path's LAST filter (see lastFilterOfPath()); resolve that once
    //    per run here, same as the path names themselves.
    for (size_t i = 0; i < resolvedTriggerPaths_.size(); ++i) {
      bool pathExistsInMenu = false;
      for (const auto& actual : hltConfig_.triggerNames()) {
        if (actual == resolvedTriggerPaths_[i]) { pathExistsInMenu = true; break; }
      }
      if (!pathExistsInMenu) {
        edm::LogWarning("MuonHLTNtupler")
            << "Configured path '" << triggerPathsCfg_[i] << "' could not be resolved against the "
            << "'MYHLT' menu this run -- no trigger name in this menu contains that string. Check the "
            << "exact path name (e.g. via hltConfig.dump(\"Triggers\") or edmConfigDump on your rerun "
            << "config). Path-level matching for this entry ('" << triggerPathsCfg_[i] << "') will "
            << "never match anything this run.";
        continue;
      }
      resolvedPathLastFilter_[i] = MuonHLTTriggerMatching::lastFilterOfPath(hltConfig_, resolvedTriggerPaths_[i]);

      std::cout << "[MuonHLTTriggerMatching] Found trigger filter " << resolvedPathLastFilter_[i] << " for path: " << resolvedTriggerPaths_[i] << std::endl;
      
      if (resolvedPathLastFilter_[i].empty()) {
        edm::LogWarning("MuonHLTNtupler")
            << "Path '" << resolvedTriggerPaths_[i] << "' resolved fine but no EDFilter-type module was "
            << "found in its module sequence -- path-level matching for this entry will never match "
            << "anything this run.";
      }
    }
  } else {
    std::cout << "not found" << std::endl;
    edm::LogWarning("MuonHLTNtupler") << "HLTConfigProvider::init failed for process name 'MYHLT' -- "
                                       << "trigger-path resolution will fall back to the unversioned "
                                       << "names configured in python, and path-level matching will "
                                       << "not work this run.";
    resolvedTriggerPaths_ = triggerPathsCfg_;
  }
}
void MuonHLTNtupler::endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}

DEFINE_FWK_MODULE(MuonHLTNtupler);
