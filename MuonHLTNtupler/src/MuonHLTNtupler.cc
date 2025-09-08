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
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h" // New
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
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociatorBaseImpl.h"
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociationMap.h"
#include "CommonTools/Utils/interface/associationMapFilterValues.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "RecoMuon/MuonIsolation/interface/Range.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoMuon/TrackerSeedGenerator/interface/SeedMvaEstimator.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"

#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractorFactory.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace muonisolation;

MuonHLTNtupler::MuonHLTNtupler(const edm::ParameterSet& iConfig):
doMVA(iConfig.getParameter<bool>("doMVA")),
doSeed(iConfig.getParameter<bool>("doSeed")),
DebugMode(iConfig.getParameter<bool>("DebugMode")),
SaveAllTracks(iConfig.getParameter<bool>("SaveAllTracks")),
// SaveStubs(iConfig.getParameter<bool>("SaveStubs")),
ttTrackToken_        ( consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_> > >(iConfig.getParameter<edm::InputTag>("L1TrackInputTag"     )) ),
// ttTrackMCTruthToken_ ( consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >(iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag"))),
// ttStubToken_         ( consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(iConfig.getParameter<edm::InputTag>("L1StubInputTag"))),
TkMuonToken_         ( consumes<l1t::TrackerMuonCollection>     (iConfig.getParameter<edm::InputTag>("TkMuonToken"))),
l1PrimaryVertexToken_( consumes<l1t::VertexWordCollection>      (iConfig.getParameter<edm::InputTag>("l1PrimaryVertex"))),

// trackerHitAssociatorConfig_(iConfig, consumesCollector()),
associatorToken(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getUntrackedParameter<edm::InputTag>("associator"))),
trackingParticleToken(consumes<TrackingParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackingParticle"))),
// t_offlineMuon_       ( consumes< std::vector<reco::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("offlineMuon"       )) ),
//propSetup_(iConfig, consumesCollector()),
//t_beamSpot_          ( consumes< reco::BeamSpot >                         (iConfig.getUntrackedParameter<edm::InputTag>("beamSpot"     )) ),

t_offlineMuon_       ( consumes< edm::View<reco::Muon> >                  (iConfig.getUntrackedParameter<edm::InputTag>("offlineMuon"       )) ),
t_offlineVertex_     ( consumes< reco::VertexCollection >                 (iConfig.getUntrackedParameter<edm::InputTag>("offlineVertex"     )) ),
t_triggerResults_    ( consumes< edm::TriggerResults >                    (iConfig.getUntrackedParameter<edm::InputTag>("triggerResults"    )) ),
t_triggerEvent_      ( consumes< trigger::TriggerEvent >                  (iConfig.getUntrackedParameter<edm::InputTag>("triggerEvent"      )) ),
t_myTriggerResults_  ( consumes< edm::TriggerResults >                    (iConfig.getUntrackedParameter<edm::InputTag>("myTriggerResults"  )) ),
t_myTriggerEvent_    ( consumes< trigger::TriggerEvent >                  (iConfig.getUntrackedParameter<edm::InputTag>("myTriggerEvent"    )) ),
t_L3Muon_            ( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("L3Muon"            )) ),
t_L2Muon_            ( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("L2Muon"            )) ),
t_L1Muon_            ( consumes< l1t::MuonBxCollection  >                 (iConfig.getUntrackedParameter<edm::InputTag>("L1Muon"            )) ),
t_TkMuon_            ( consumes< reco::RecoChargedCandidateCollection >   (iConfig.getUntrackedParameter<edm::InputTag>("TkMuon"            )) ),

t_iterL3OI_          ( consumes< std::vector<reco::MuonTrackLinks> >      (iConfig.getUntrackedParameter<edm::InputTag>("iterL3OI"          )) ),
t_iterL3IOFromL2_    ( consumes< std::vector<reco::MuonTrackLinks> >      (iConfig.getUntrackedParameter<edm::InputTag>("iterL3IOFromL2"    )) ),
t_iterL3FromL2_      ( consumes< std::vector<reco::MuonTrackLinks> >      (iConfig.getUntrackedParameter<edm::InputTag>("iterL3FromL2"      )) ),
t_iterL3IOFromL1_    ( consumes< std::vector<reco::Track> >               (iConfig.getUntrackedParameter<edm::InputTag>("iterL3IOFromL1"    )) ),
t_iterL3MuonNoID_    ( consumes< std::vector<reco::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("iterL3MuonNoID"    )) ),
t_iterL3Muon_        ( consumes< std::vector<reco::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("iterL3Muon"        )) ),

t_hltIterL3MuonTrimmedPixelVertices_       ( consumes< reco::VertexCollection > (iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3MuonTrimmedPixelVertices"     )) ),
t_hltIterL3FromL1MuonTrimmedPixelVertices_ ( consumes< reco::VertexCollection > (iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3FromL1MuonTrimmedPixelVertices"     )) ),

t_hltIterL3OISeedsFromL2Muons_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3OISeedsFromL2Muons")) ),
t_hltIter0IterL3MuonPixelSeedsFromPixelTracks_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter0IterL3MuonPixelSeedsFromPixelTracks")) ),
t_hltIter2IterL3MuonPixelSeeds_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter2IterL3MuonPixelSeeds")) ),
t_hltIter3IterL3MuonPixelSeeds_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter3IterL3MuonPixelSeeds")) ),
t_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks")) ),
t_hltIter2IterL3FromL1MuonPixelSeeds_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter2IterL3FromL1MuonPixelSeeds")) ),
t_hltIter3IterL3FromL1MuonPixelSeeds_ ( consumes< TrajectorySeedCollection >     (iConfig.getUntrackedParameter<edm::InputTag>("hltIter3IterL3FromL1MuonPixelSeeds")) ),

t_hltIterL3OIMuonTrack_    ( consumes< edm::View<reco::Track> >                  (iConfig.getUntrackedParameter<edm::InputTag>("hltIterL3OIMuonTrack"    )) ),
t_hltIter0IterL3MuonTrack_    ( consumes< edm::View<reco::Track> >               (iConfig.getUntrackedParameter<edm::InputTag>("hltIter0IterL3MuonTrack"    )) ),
t_hltIter2IterL3MuonTrack_    ( consumes< edm::View<reco::Track> >               (iConfig.getUntrackedParameter<edm::InputTag>("hltIter2IterL3MuonTrack"    )) ),
t_hltIter3IterL3MuonTrack_    ( consumes< edm::View<reco::Track> >               (iConfig.getUntrackedParameter<edm::InputTag>("hltIter3IterL3MuonTrack"    )) ),
t_hltIter0IterL3FromL1MuonTrack_    ( consumes< edm::View<reco::Track> >         (iConfig.getUntrackedParameter<edm::InputTag>("hltIter0IterL3FromL1MuonTrack"    )) ),
t_hltIter2IterL3FromL1MuonTrack_    ( consumes< edm::View<reco::Track> >         (iConfig.getUntrackedParameter<edm::InputTag>("hltIter2IterL3FromL1MuonTrack"    )) ),
t_hltIter3IterL3FromL1MuonTrack_    ( consumes< edm::View<reco::Track> >         (iConfig.getUntrackedParameter<edm::InputTag>("hltIter3IterL3FromL1MuonTrack"    )) ),

// Muon timing

theMuonCollectionToken_ ( consumes<reco::RecoChargedCandidateCollection>                            (iConfig.getParameter<edm::InputTag>("inputMuonCollection"))            ),
theMuonFilteredCollectionToken_ ( consumes<trigger::TriggerFilterObjectWithRefs>             (iConfig.getParameter<edm::InputTag>("inputMuonFilterCollection"))            ),
//muonAssocToken_    ( consumes<edm::ValueMap<int>>                            (iConfig.getParameter<edm::InputTag>("muontrkAssSrc"))            ),

muonbtlMatchChi2Token_                  ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonbtlMatchChi2"))                    ),
muonetlMatchChi2Token_                  ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonetlMatchChi2"))                    ),
muonbtlMatchTimeChi2Token_              ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonbtlMatchTimeChi2"))                    ),
muonetlMatchTimeChi2Token_              ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonetlMatchTimeChi2"))                    ),
muonnpixBarrelToken_                    ( consumes<edm::ValueMap<int>>                          (iConfig.getParameter<edm::InputTag>("muonnpixBarrel"))                    ),
muonnpixEndcapToken_                    ( consumes<edm::ValueMap<int>>                          (iConfig.getParameter<edm::InputTag>("muonnpixEndcap"))                    ),
muonTrackOutermostHitPositionToken_ ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTrackOutermostHitPosition"))                    ),
muonTrackpToken_                    ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTrackp"))                    ),
muonTrackBetaToken_                 ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTrackBeta"))                    ),
muonTrackt0Token_                   ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTrackt0"))                    ),
muonTracksigmat0Token_              ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTracksigmat0"))                    ),
muonTrackPathLengthToken_           ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTrackPathLength"))                    ),
muonTracktmtdToken_                 ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTracktmtd"))                    ),
muonTracksigmatmtdToken_            ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTracksigmatmtd"))                    ),
muonTrackmtdposToken_               ( consumes<edm::ValueMap<GlobalPoint>>                    (iConfig.getParameter<edm::InputTag>("muonTrackmtdpos"))                    ),
muonTrackTofMuToken_                ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTrackTofMu"))                    ),
muonTrackSigmaTofMuToken_           ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("muonTrackSigmaTofMu"))                    ),

// Isolation PF Candidates

pfCandidateProducer_ (consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidateProducer"))),
drMaxPf_     (iConfig.getParameter<double>("drMaxPF")),
drVetoPf_    (iConfig.getParameter<double>("drVetoPF")),
drVetoPfCh_  (iConfig.getParameter<double>("drVetoPFCh")),
minEnergyPf_ (iConfig.getParameter<double>("minEnergyPF")),

rhoProducer_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoProducer"))),
trackBtlMatchChi2Token_                  ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("btlMatchChi2"))                    ),
trackEtlMatchChi2Token_                  ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("etlMatchChi2"))                    ),
trackBtlMatchTimeChi2Token_              ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("btlMatchTimeChi2"))                    ),
trackEtlMatchTimeChi2Token_              ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("etlMatchTimeChi2"))                    ),
trackNpixBarrelToken_                    ( consumes<edm::ValueMap<int>>                          (iConfig.getParameter<edm::InputTag>("npixBarrel"))                    ),
trackNpixEndcapToken_                    ( consumes<edm::ValueMap<int>>                          (iConfig.getParameter<edm::InputTag>("npixEndcap"))                    ),
trackOutermostHitPositionToken_ ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackOutermostHitPosition"))                    ),
trackpToken_                    ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackp"))                    ),
trackBetaToken_                 ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackBeta"))                    ),
trackt0Token_                   ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackt0"))                    ),
tracksigmat0Token_              ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("tracksigmat0"))                    ),
trackPathLengthToken_           ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackPathLength"))                    ),
tracktmtdToken_                 ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("tracktmtd"))                    ),
tracksigmatmtdToken_            ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("tracksigmatmtd"))                    ),
trackmtdposToken_               ( consumes<edm::ValueMap<GlobalPoint>>                    (iConfig.getParameter<edm::InputTag>("trackmtdpos"))                    ),
trackTofPiToken_                ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackTofPi"))                    ),
trackSigmaTofPiToken_           ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackSigmaTofPi"))                    ),
trackTofKToken_                ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackTofK"))                    ),
trackSigmaTofKToken_           ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackSigmaTofK"))                    ),
trackTofPToken_                ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackTofP"))                    ),
trackSigmaTofPToken_           ( consumes<edm::ValueMap<float>>                          (iConfig.getParameter<edm::InputTag>("trackSigmaTofP"))                    ),


// End of isolation inputs

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

  trackCollectionNames_ = iConfig.getUntrackedParameter<std::vector<std::string>>("trackCollectionNames");
  trackCollectionLabels_ = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("trackCollectionLabels");
  associationLabels_      = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("associationLabels");
  if( trackCollectionNames_.size() != trackCollectionLabels_.size() || trackCollectionLabels_.size() != associationLabels_.size()) {
    throw cms::Exception("ConfigurationError")
      << "Number of track collection names is different from number of track collection names or association labels";
  }
  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    trackCollectionTokens_.push_back(     consumes< edm::View<reco::Track> >(  trackCollectionLabels_[i]) );
    simToRecoCollectionTokens_.push_back( consumes<reco::SimToRecoCollection>( associationLabels_[i]) );
    recoToSimCollectionTokens_.push_back( consumes<reco::RecoToSimCollection>( associationLabels_[i]) );
    trkTemplates_.push_back( new trkTemplate() );
    tpTemplates_.push_back(  new tpTemplate()  );
  }

  trkIsoTags_ = iConfig.getUntrackedParameter<std::vector<std::string>>("trkIsoTags");
  trkIsoLabels_ = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("trkIsoLabels");
  if( trkIsoTags_.size() != trkIsoLabels_.size() ) {
    throw cms::Exception("ConfigurationError")
      << "trkIsoTags_.size() != trkIsoLabels_.size()";
  }
  for( unsigned int i = 0; i < trkIsoTags_.size(); ++i) {
    trkIsoTokens_.push_back( consumes<reco::IsoDepositMap>( trkIsoLabels_.at(i) ) );
  }

  pfIsoTags_ = iConfig.getUntrackedParameter<std::vector<std::string>>("pfIsoTags");
  pfIsoLabels_ = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("pfIsoLabels");
  if( pfIsoTags_.size() != pfIsoLabels_.size() ) {
    throw cms::Exception("ConfigurationError")
      << "pfIsoTags_.size() != pfIsoLabels_.size()";
  }
  for( unsigned int i = 0; i < pfIsoTags_.size(); ++i) {
    pfIsoTokens_.push_back( consumes<reco::RecoChargedCandidateIsolationMap>( pfIsoLabels_.at(i) ) );
  }
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0_                = iConfig.getUntrackedParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0");
  mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0_                = iConfig.getUntrackedParameter<edm::FileInPath>("mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0");
  mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B_ =                iConfig.getUntrackedParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B");
  mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B_ =                 iConfig.getUntrackedParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B");
  mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E_ =                iConfig.getUntrackedParameter<std::vector<double>>("mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E");
  mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E_ =                 iConfig.getUntrackedParameter<std::vector<double>>("mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E");
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

  // -- hltIterL3MuonTrimmedPixelVertices
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

  // -- fill each object

  std::cout << "------- MUON ANALYZER --------" << std::endl;
  
  Fill_L1Track(iEvent, iSetup);
  Fill_Muon(iEvent);
  Fill_Muon2(iEvent);
  Fill_HLT(iEvent, 0); // -- original HLT objects saved in data taking
  Fill_HLT(iEvent, 1); // -- rerun objects
  Fill_HLTMuon(iEvent);
  Fill_L1Muon(iEvent);
  Fill_IterL3(iEvent, iSetup);
  if( doSeed )  Fill_Seed(iEvent, iSetup);
  if( !isRealData_ ) {
    Fill_GenParticle(iEvent);
    // Fill_TP(iEvent, TrkParticle);
  }
  Fill_PFCand(iEvent, iSetup);

  /**
  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    bool doIso = (i == trackCollectionNames_.size()-1);
    fill_trackTemplate( iEvent, trackCollectionTokens_.at(i), recoToSimCollectionTokens_.at(i), trkTemplates_.at(i), doIso );
    fill_tpTemplate(    iEvent,                               simToRecoCollectionTokens_.at(i), tpTemplates_.at(i) );
  }

  // -- Fill L3 Muon
  edm::Handle<std::vector<l1t::VertexWord> > l1PrimaryVertex;
  iEvent.getByToken(l1PrimaryVertexToken_,l1PrimaryVertex);
  double l1vtx_z = l1PrimaryVertex->size() > 0 ? l1PrimaryVertex->at(0).z0() : -1e9;

  edm::Handle< std::vector<reco::Muon> > h_iterL3Muon;
  if( iEvent.getByToken( t_iterL3Muon_, h_iterL3Muon) ) {
    for( auto i=0U; i<h_iterL3Muon->size(); ++i )
    {
      const auto& muon(h_iterL3Muon->at(i));
      MTL3Muons->fill(muon, l1vtx_z);
    }
  }

  edm::Handle< std::vector<reco::Muon> > h_iterL3MuonNoID;
  if( iEvent.getByToken( t_iterL3MuonNoID_, h_iterL3MuonNoID) )
  {
    for( auto i=0U; i<h_iterL3MuonNoID->size(); ++i )
    {
      const auto& muon(h_iterL3MuonNoID->at(i));
      MTL3MuonsNoId->fill(muon, l1vtx_z);
    }
  }
  **/
  
  ntuple_->Fill();
}

void MuonHLTNtupler::beginJob()
{
  edm::Service<TFileService> fs;
  ntuple_ = fs->make<TTree>("ntuple","ntuple");

  Make_Branch();

  if(doMVA) {
    mvaPhase2HltIter2IterL3FromL1MuonPixelSeeds_ = std::make_pair(
      std::make_unique<SeedMvaEstimatorPhase2>(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_B_0_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_B_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_B_),
      std::make_unique<SeedMvaEstimatorPhase2>(mvaFileHltIter2IterL3FromL1MuonPixelSeeds_E_0_, mvaScaleMeanHltIter2IterL3FromL1MuonPixelSeeds_E_, mvaScaleStdHltIter2IterL3FromL1MuonPixelSeeds_E_) );
  }
}

void MuonHLTNtupler::Init()
{
  if (SaveAllTracks){
  m_trk_pt.clear();
  m_trk_eta.clear();
  m_trk_phi.clear();
  m_trk_d0.clear();
  m_trk_z0.clear();
  m_trk_rInv.clear();
  m_trk_tanL.clear();
  m_trk_MVA1.clear();
  m_trk_MVA2.clear();
  m_trk_MVA3.clear();
  m_trk_chi2.clear();
  m_trk_bendchi2.clear();
  m_trk_nstub.clear();
  m_trk_lhits.clear();
  m_trk_dhits.clear();
  m_trk_seed.clear();
  m_trk_phiSector.clear();
  m_trk_genuine.clear();
  m_trk_loose.clear();
  m_trk_unknown.clear();
  m_trk_combinatoric.clear();
  m_trk_fake.clear();
  m_trk_matchtp_pdgid.clear();
  m_trk_matchtp_pt.clear();
  m_trk_matchtp_eta.clear();
  m_trk_matchtp_phi.clear();
  m_trk_matchtp_z0.clear();
  m_trk_matchtp_dxy.clear();

  m_stub_x.clear();
  m_stub_y.clear();
  m_stub_z.clear();
  m_stub_isBarrel.clear();
  m_stub_layer.clear();

  // l1TkMuon
  mL1TkMu_pt.clear();
  mL1TkMu_eta.clear();
  mL1TkMu_phi.clear();
  mL1TkMu_trkIsol.clear();
  mL1TkMu_trkzVtx.clear();
  mL1TkMu_dR.clear();

  mL1TkMu_nTracksMatched.clear();
  mL1TkMu_trackCurvature.clear();

  mL1TkMu_quality.clear();
  mL1TkMu_pattern.clear();
  mL1TkMu_muonDetector.clear();

  mL1TkMu_TTTpointer.clear();

  mL1TkMu_muRefHwPt.clear();
  mL1TkMu_muRefHwDXY.clear();
  mL1TkMu_muRefHwEta.clear();
  mL1TkMu_muRefHwPhi.clear();
  mL1TkMu_muRefHwSign.clear();
  mL1TkMu_muRefHwSignValid.clear();
  mL1TkMu_muRefHwQual.clear();

  mTTTrackMap.clear();
  }


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

  MuonIterSeedMap.clear();
  // MuonIterNoIdSeedMap.clear();
  hltIterL3OIMuonTrackMap.clear();
  hltIter0IterL3MuonTrackMap.clear();
  hltIter2IterL3MuonTrackMap.clear();
  hltIter3IterL3MuonTrackMap.clear();
  hltIter0IterL3FromL1MuonTrackMap.clear();
  hltIter2IterL3FromL1MuonTrackMap.clear();
  hltIter3IterL3FromL1MuonTrackMap.clear();
  iterL3IDpassed.clear();
  iterL3NoIDpassed.clear();

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

  nL3Muon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    L3Muon_pt_[i] = -999;
    L3Muon_eta_[i] = -999;
    L3Muon_phi_[i] = -999;
    L3Muon_charge_[i] = -999;
    L3Muon_trkPt_[i] = -999;
  }

  nL2Muon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    L2Muon_pt_[i] = -999;
    L2Muon_eta_[i] = -999;
    L2Muon_phi_[i] = -999;
    L2Muon_charge_[i] = -999;
    L2Muon_trkPt_[i] = -999;
  }

  nTkMuon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    TkMuon_pt_[i] = -999;
    TkMuon_eta_[i] = -999;
    TkMuon_phi_[i] = -999;
    TkMuon_charge_[i] = -999;
    TkMuon_trkPt_[i] = -999;
  }

  nL1Muon_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    L1Muon_pt_[i] = -999;
    L1Muon_eta_[i] = -999;
    L1Muon_phi_[i] = -999;
    L1Muon_charge_[i] = -999;
    L1Muon_quality_[i] = -999;
    L1Muon_etaAtVtx_[i] = -999;
    L1Muon_phiAtVtx_[i] = -999;
  }

  nIterL3OI_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    iterL3OI_inner_pt_[i] = -999;
    iterL3OI_inner_eta_[i] = -999;
    iterL3OI_inner_phi_[i] = -999;
    iterL3OI_inner_charge_[i] = -999;
    iterL3OI_outer_pt_[i] = -999;
    iterL3OI_outer_eta_[i] = -999;
    iterL3OI_outer_phi_[i] = -999;
    iterL3OI_outer_charge_[i] = -999;
    iterL3OI_global_pt_[i] = -999;
    iterL3OI_global_eta_[i] = -999;
    iterL3OI_global_phi_[i] = -999;
    iterL3OI_global_charge_[i] = -999;
  }

  nIterL3IOFromL2_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    iterL3IOFromL2_inner_pt_[i] = -999;
    iterL3IOFromL2_inner_eta_[i] = -999;
    iterL3IOFromL2_inner_phi_[i] = -999;
    iterL3IOFromL2_inner_charge_[i] = -999;
    iterL3IOFromL2_outer_pt_[i] = -999;
    iterL3IOFromL2_outer_eta_[i] = -999;
    iterL3IOFromL2_outer_phi_[i] = -999;
    iterL3IOFromL2_outer_charge_[i] = -999;
    iterL3IOFromL2_global_pt_[i] = -999;
    iterL3IOFromL2_global_eta_[i] = -999;
    iterL3IOFromL2_global_phi_[i] = -999;
    iterL3IOFromL2_global_charge_[i] = -999;
  }

  nIterL3FromL2_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    iterL3FromL2_inner_pt_[i] = -999;
    iterL3FromL2_inner_eta_[i] = -999;
    iterL3FromL2_inner_phi_[i] = -999;
    iterL3FromL2_inner_charge_[i] = -999;
    iterL3FromL2_outer_pt_[i] = -999;
    iterL3FromL2_outer_eta_[i] = -999;
    iterL3FromL2_outer_phi_[i] = -999;
    iterL3FromL2_outer_charge_[i] = -999;
    iterL3FromL2_global_pt_[i] = -999;
    iterL3FromL2_global_eta_[i] = -999;
    iterL3FromL2_global_phi_[i] = -999;
    iterL3FromL2_global_charge_[i] = -999;
  }

  nIterL3IOFromL1_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    iterL3IOFromL1_pt_[i] = -999;
    iterL3IOFromL1_eta_[i] = -999;
    iterL3IOFromL1_phi_[i] = -999;
    iterL3IOFromL1_charge_[i] = -999;
  }


  nIterL3MuonNoID_ = 0;
  for (int i=0; i<arrSize_; ++i)
  {
    iterL3MuonNoID_pt_[i] = -999;
    iterL3MuonNoID_innerPt_[i] = -999;
    iterL3MuonNoID_eta_[i] = -999;
    iterL3MuonNoID_phi_[i] = -999;
    iterL3MuonNoID_charge_[i] = -999;

    iterL3MuonNoID_isGLB_[i] = 0;
    iterL3MuonNoID_isSTA_[i] = 0;
    iterL3MuonNoID_isTRK_[i] = 0;
  }

  nIterL3Muon_ = 0;
  for (int i=0; i<arrSize_; ++i)
  {
    iterL3Muon_pt_[i] = -999;
    iterL3Muon_innerPt_[i] = -999;
    iterL3Muon_eta_[i] = -999;
    iterL3Muon_phi_[i] = -999;
    iterL3Muon_charge_[i] = -999;

    iterL3Muon_isGLB_[i] = 0;
    iterL3Muon_isSTA_[i] = 0;
    iterL3Muon_isTRK_[i] = 0;
  }

  nMuonCand_ = 0;
  for( int i=0; i<arrSize_; i++)
    {
      muonCand_pt_[i] = -999;
      muonCand_eta_[i] = -999;
      muonCand_phi_[i] = -999;
      muonCand_px_[i] = -999;
      muonCand_py_[i] = -999;
      muonCand_pz_[i] = -999;
      muonCand_vx_[i] = -999;
      muonCand_vy_[i] = -999;
      muonCand_vz_[i] = -999;
      muonCand_dB_[i] = -999;
      muonCand_charge_[i] = -999;
      muonCand_isGLB_[i] = 0;
      muonCand_isSTA_[i] = 0;
      muonCand_isTRK_[i] = 0;
      muonCand_isMuon_[i] = 0;
      muonCand_inner_trkChi2_[i] = -999;
      muonCand_inner_validFraction_[i] = -999;
      muonCand_inner_trackerLayers_[i] = -999;
      muonCand_inner_trackerHits_[i] = -999;
      muonCand_inner_lostTrackerHits_[i] = -999;
      muonCand_inner_lostTrackerHitsIn_[i] = -999;
      muonCand_inner_lostTrackerHitsOut_[i] = -999;
      muonCand_inner_lostPixelHits_[i] = -999;
      muonCand_inner_lostPixelBarrelHits_[i] = -999;
      muonCand_inner_lostPixelEndcapHits_[i] = -999;
      muonCand_inner_lostStripHits_[i] = -999;
      muonCand_inner_lostStripTIBHits_[i] = -999;
      muonCand_inner_lostStripTIDHits_[i] = -999;
      muonCand_inner_lostStripTOBHits_[i] = -999;
      muonCand_inner_lostStripTECHits_[i] = -999;
      muonCand_inner_pixelLayers_[i] = -999;
      muonCand_inner_pixelHits_[i] = -999;
      muonCand_global_muonHits_[i] = -999;
      muonCand_global_trkChi2_[i] = -999;
      muonCand_global_trackerLayers_[i] = -999;
      muonCand_global_trackerHits_[i] = -999;
      
      muonCand_dxy_[i] = -999;
      muonCand_dxyError_bs_[i] = -999;
      muonCand_dz_[i] = -999;
      muonCand_dzError_[i] = -999;
      muonCand_IPSig_[i] = -999;

      muonCand_btlMatchChi2_[i] = -999;
      muonCand_etlMatchChi2_[i] = -999;
      muonCand_btlMatchTimeChi2_[i] = -999;
      muonCand_etlMatchTimeChi2_[i] = -999;
      muonCand_npixBarrel_[i] = -999;
      muonCand_npixEndcap_[i] = -999;
      muonCand_outermostHitPosition_[i] = -999;
      muonCand_p_[i] = -999;
      muonCand_beta_[i] = -999;
      muonCand_t0_[i] = -999;
      muonCand_sigmat0_[i] = -999;
      muonCand_pathLength_[i] = -999;
      muonCand_tmtd_[i] = -999;
      muonCand_sigmatmtd_[i] = -999;
      muonCand_tofMu_[i] = -999;
      muonCand_sigmaTofMu_[i] = -999;
      muonCand_mtdpos_x_[i] = -999;
      muonCand_mtdpos_y_[i] = -999;
      muonCand_mtdpos_z_[i] = -999;
    }

  nECAL_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    ecal_et_[i] = -999;
    ecal_pt_[i] = -999;
    ecal_eta_[i] = -999;
    ecal_phi_[i] = -999;
    ecal_charge_[i] = -999;
    ecal_px_[i] = -999;
    ecal_py_[i] = -999;
    ecal_pz_[i] = -999;
    ecal_vx_[i] = -999;
    ecal_vy_[i] = -999;
    ecal_vz_[i] = -999;
    ecal_time_[i] = -999;
    ecal_timeErr_[i] = -999;
    ecal_depth_[i] = -999;
    ecal_rho_[i] = -999;
    ecal_muonIdx_[i] = -999;
    ecal_nHits_[i] = -999;
  }

  nECALHits_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    ecal_hit_energy_[i] = -999;
    ecal_hit_depth_[i] = -999;
    ecal_hit_time_[i] = -999;
    ecal_hit_timeErr_[i] = -999;
    ecal_hit_pt2_[i] = -999;
    ecal_hit_x_[i] = -999;
    ecal_hit_y_[i] = -999;
    ecal_hit_z_[i] = -999;
    ecal_hit_eta_[i] = -999;
    ecal_hit_phi_[i] = -999;
    ecal_hit_fraction_[i] = -999;
    ecal_hit_idx_[i] = -999;
  }

  nHCAL_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    hcal_et_[i] = -999;
    hcal_pt_[i] = -999;
    hcal_eta_[i] = -999;
    hcal_phi_[i] = -999;
    hcal_charge_[i] = -999;
    hcal_px_[i] = -999;
    hcal_py_[i] = -999;
    hcal_pz_[i] = -999;
    hcal_vx_[i] = -999;
    hcal_vy_[i] = -999;
    hcal_vz_[i] = -999;
    hcal_time_[i] = -999;
    hcal_timeErr_[i] = -999;
    hcal_depth_[i] = -999;
    hcal_rho_[i] = -999;
    hcal_muonIdx_[i] = -999;
    hcal_nHits_[i] = -999;
  }

  nHCALHits_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    hcal_hit_energy_[i] = -999;
    hcal_hit_depth_[i] = -999;
    hcal_hit_time_[i] = -999;
    hcal_hit_timeErr_[i] = -999;
    hcal_hit_pt2_[i] = -999;
    hcal_hit_x_[i] = -999;
    hcal_hit_y_[i] = -999;
    hcal_hit_z_[i] = -999;
    hcal_hit_eta_[i] = -999;
    hcal_hit_phi_[i] = -999;
    hcal_hit_fraction_[i] = -999;
    hcal_hit_idx_[i] = -999;
  }
  
  nHGCAL_em_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    hgcal_em_et_[i] = -999;
    hgcal_em_pt_[i] = -999;
    hgcal_em_eta_[i] = -999;
    hgcal_em_phi_[i] = -999;
    hgcal_em_charge_[i] = -999;
    hgcal_em_px_[i] = -999;
    hgcal_em_py_[i] = -999;
    hgcal_em_pz_[i] = -999;
    hgcal_em_vx_[i] = -999;
    hgcal_em_vy_[i] = -999;
    hgcal_em_vz_[i] = -999;
    hgcal_em_time_[i] = -999;
    hgcal_em_timeErr_[i] = -999;
    hgcal_em_depth_[i] = -999;
    hgcal_em_algoID_[i] = -999;
    hgcal_em_muonIdx_[i] = -999;
  }
  nHGCAL_had_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    hgcal_had_et_[i] = -999;
    hgcal_had_pt_[i] = -999;
    hgcal_had_eta_[i] = -999;
    hgcal_had_phi_[i] = -999;
    hgcal_had_charge_[i] = -999;
    hgcal_had_px_[i] = -999;
    hgcal_had_py_[i] = -999;
    hgcal_had_pz_[i] = -999;
    hgcal_had_vx_[i] = -999;
    hgcal_had_vy_[i] = -999;
    hgcal_had_vz_[i] = -999;
    hgcal_had_time_[i] = -999;
    hgcal_had_timeErr_[i] = -999;
    hgcal_had_depth_[i] = -999;
    hgcal_had_algoID_[i] = -999;
    hgcal_had_muonIdx_[i] = -999;
  }

  nTrack_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    track_pt_[i] = -999;
    track_eta_[i] = -999;
    track_phi_[i] = -999;
    track_charge_[i] = -999;
    track_px_[i] = -999;
    track_py_[i] = -999;
    track_pz_[i] = -999;
    track_vx_[i] = -999;
    track_vy_[i] = -999;
    track_vz_[i] = -999;
    track_dxy_bs_[i] = -999;
    track_dxyError_bs_[i] = -999;
    track_dz_bs_[i] = -999;
    track_dzError_[i] = -999;
    track_trkChi2_[i] = -999;
    track_trackerLayers_[i] = -999;
    track_trackerHits_[i] = -999;
    track_lostTrackerHits_[i] = -999;
    track_lostTrackerHitsIn_[i] = -999;
    track_lostTrackerHitsOut_[i] = -999;
    track_lostPixelHits_[i] = -999;
    track_lostPixelBarrelHits_[i] = -999;
    track_lostPixelEndcapHits_[i] = -999;
    track_lostStripHits_[i] = -999;
    track_lostStripTIBHits_[i] = -999;
    track_lostStripTIDHits_[i] = -999;
    track_lostStripTOBHits_[i] = -999;
    track_lostStripTECHits_[i] = -999;
    track_pixelLayers_[i] = -999;
    track_pixelHits_[i] = -999;
    track_muonHits_[i] = -999;
    track_t0Src_[i] = -999;
    track_Sigmat0Src_[i] = -999;
    track_t0Pid_[i] = -999;
    track_t0Safe_[i] = -999;
    track_sigmat0Safe_[i] = -999;
    track_mtdQualMVA_[i] = -999;
    track_tMtd_[i] = -999;
    track_tofPi_[i] = -999;
    track_tofK_[i] = -999;
    track_tofP_[i] = -999;
    track_probPi_[i] = -999;
    track_probK_[i] = -999;
    track_probP_[i] = -999;
    track_sigmatofpi_[i] = -999;
    track_sigmatofk_[i] = -999;
    track_sigmatofp_[i] = -999;
    track_btlMatchChi2_[i] = -999;
    track_btlMatchTimeChi2_[i] = -999;
    track_etlMatchChi2_[i] = -999;
    track_etlMatchTimeChi2_[i] = -999;
    track_npixBarrel_[i] = -999;
    track_npixEndcap_[i] = -999;
    track_muonIdx_[i] = -999;
    track_outermostHitPosition_[i] = -999;
    track_p_[i] = -999;
    track_beta_[i] = -999;
    track_pathLength_[i] = -999;
    track_mtdpos_x_[i] = -999;
    track_mtdpos_y_[i] = -999;
    track_mtdpos_z_[i] = -999;
    track_TPcharge_[i] = -999;
    track_TPpdgId_[i] = -999;
    track_TPenergy_[i] = -999;
    track_TPpt_[i] = -999;
    track_TPeta_[i] = -999;
    track_TPphi_[i] = -999;
    track_TPparentVx_[i] = -999;
    track_TPparentVy_[i] = -999;
    track_TPparentVz_[i] = -999;
    track_TPstatus_[i] = -999;
    track_TPnumberOfHits_[i] = -999;
    track_TPnumberOfTrackerHits_[i] = -999;
    track_TPnumberOfTrackerLayers_[i] = -999;    
  }

  nPFCand_ = 0;
  for( int i=0; i<arrSize_; i++)
  {
    pfcand_pt_[i] = -999;
    pfcand_eta_[i] = -999;
    pfcand_phi_[i] = -999;
    pfcand_charge_[i] = -999;
    pfcand_pdgId_[i] = -999;
    pfcand_px_[i] = -999;
    pfcand_py_[i] = -999;
    pfcand_pz_[i] = -999;
    pfcand_vx_[i] = -999;
    pfcand_vy_[i] = -999;
    pfcand_vz_[i] = -999;
    pfcand_time_[i] = -999;
    pfcand_timeErr_[i] = -999;
    pfcand_dxy_[i] = -999;
    pfcand_dz_[i] = -999;
    pfcand_dxyErr_[i] = -999;
    pfcand_dzErr_[i] = -999;
    pfcand_vChi2NoF_[i] = -999;
    pfcand_muonIdx_[i] = -999;

    pfcand_rho_[i] = -999;
    pfcand_btlMatchChi2_[i] = -999;
    pfcand_etlMatchChi2_[i] = -999;
    pfcand_btlMatchTimeChi2_[i] = -999;
    pfcand_etlMatchTimeChi2_[i] = -999;
    pfcand_npixBarrel_[i] = -999;
    pfcand_npixEndcap_[i] = -999;
    pfcand_outermostHitPosition_[i] = -999;
    pfcand_p_[i] = -999;
    pfcand_beta_[i] = -999;
    pfcand_t0_[i] = -999;
    pfcand_sigmat0_[i] = -999;
    pfcand_pathLength_[i] = -999;
    pfcand_tmtd_[i] = -999;
    pfcand_sigmatmtd_[i] = -999;
    pfcand_tofPi_[i] = -999;
    pfcand_sigmaTofPi_[i] = -999;
    pfcand_tofK_[i] = -999;
    pfcand_sigmaTofK_[i] = -999;
    pfcand_tofP_[i] = -999;
    pfcand_sigmaTofP_[i] = -999;
    pfcand_mtdpos_x_[i] = -999;
    pfcand_mtdpos_y_[i] = -999;
    pfcand_mtdpos_z_[i] = -999;
  }
  
  SThltIterL3OISeedsFromL2Muons->clear();
  SThltIter0IterL3MuonPixelSeedsFromPixelTracks->clear();
  SThltIter2IterL3MuonPixelSeeds->clear();
  SThltIter3IterL3MuonPixelSeeds->clear();
  SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->clear();
  SThltIter2IterL3FromL1MuonPixelSeeds->clear();
  SThltIter3IterL3FromL1MuonPixelSeeds->clear();

  TThltIterL3OIMuonTrack->clear();
  TThltIter0IterL3MuonTrack->clear();
  TThltIter2IterL3MuonTrack->clear();
  TThltIter3IterL3MuonTrack->clear();
  TThltIter0IterL3FromL1MuonTrack->clear();
  TThltIter2IterL3FromL1MuonTrack->clear();
  TThltIter3IterL3FromL1MuonTrack->clear();

  MTL3MuonsNoId->clear();
  MTL3Muons->clear();

  TrkParticle->clear();

  VThltIterL3MuonTrimmedPixelVertices->clear();
  VThltIterL3FromL1MuonTrimmedPixelVertices->clear();

  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    trkTemplates_.at(i)->clear();
    tpTemplates_.at(i)->clear();
  }
}

void MuonHLTNtupler::Make_Branch()
{
  if (SaveAllTracks){
  ntuple_->Branch("trk_pt",    &m_trk_pt);
  ntuple_->Branch("trk_eta",   &m_trk_eta);
  ntuple_->Branch("trk_phi",   &m_trk_phi);
  ntuple_->Branch("trk_d0",    &m_trk_d0);
  ntuple_->Branch("trk_z0",    &m_trk_z0);
  ntuple_->Branch("trk_rInv", &m_trk_rInv);
  ntuple_->Branch("trk_tanL", &m_trk_tanL);
  ntuple_->Branch("trk_MVA1", &m_trk_MVA1);
  ntuple_->Branch("trk_MVA2", &m_trk_MVA2);
  ntuple_->Branch("trk_MVA3", &m_trk_MVA3);

  ntuple_->Branch("trk_chi2",  &m_trk_chi2);
  ntuple_->Branch("trk_bendchi2",  &m_trk_bendchi2);
  ntuple_->Branch("trk_nstub", &m_trk_nstub);
  ntuple_->Branch("trk_lhits", &m_trk_lhits);
  ntuple_->Branch("trk_dhits", &m_trk_dhits);
  ntuple_->Branch("trk_seed",    &m_trk_seed);
  ntuple_->Branch("trk_phiSector", &m_trk_phiSector);
  ntuple_->Branch("trk_genuine",      &m_trk_genuine);
  ntuple_->Branch("trk_loose",        &m_trk_loose);
  ntuple_->Branch("trk_unknown",      &m_trk_unknown);
  ntuple_->Branch("trk_combinatoric", &m_trk_combinatoric);
  ntuple_->Branch("trk_fake",         &m_trk_fake);
  ntuple_->Branch("trk_matchtp_pdgid",&m_trk_matchtp_pdgid);
  ntuple_->Branch("trk_matchtp_pt",   &m_trk_matchtp_pt);
  ntuple_->Branch("trk_matchtp_eta",  &m_trk_matchtp_eta);
  ntuple_->Branch("trk_matchtp_phi",  &m_trk_matchtp_phi);
  ntuple_->Branch("trk_matchtp_z0",   &m_trk_matchtp_z0);
  ntuple_->Branch("trk_matchtp_dxy",  &m_trk_matchtp_dxy);

  ntuple_->Branch("stub_x", &m_stub_x);
  ntuple_->Branch("stub_y", &m_stub_y);
  ntuple_->Branch("stub_z", &m_stub_z);
  ntuple_->Branch("stub_isBarrel",   &m_stub_isBarrel);
  ntuple_->Branch("stub_layer",      &m_stub_layer);

  ntuple_->Branch("L1TkMu_pt", &mL1TkMu_pt);
  ntuple_->Branch("L1TkMu_eta", &mL1TkMu_eta);
  ntuple_->Branch("L1TkMu_phi", &mL1TkMu_phi);

  ntuple_->Branch("L1TkMu_trkIsol", &mL1TkMu_trkIsol);
  ntuple_->Branch("L1TkMu_trkzVtx", &mL1TkMu_trkzVtx);
  ntuple_->Branch("L1TkMu_dR", &mL1TkMu_dR);

  ntuple_->Branch("L1TkMu_nTracksMatched", &mL1TkMu_nTracksMatched);
  ntuple_->Branch("L1TkMu_trackCurvature", &mL1TkMu_trackCurvature);

  ntuple_->Branch("L1TkMu_quality", &mL1TkMu_quality);
  ntuple_->Branch("L1TkMu_pattern", &mL1TkMu_pattern);
  ntuple_->Branch("L1TkMu_muonDetector", &mL1TkMu_muonDetector);

  ntuple_->Branch("L1TkMu_TTTpointer", &mL1TkMu_TTTpointer);

  ntuple_->Branch("L1TkMu_muRefHwPt", &mL1TkMu_muRefHwPt);
  ntuple_->Branch("L1TkMu_muRefHwDXY", &mL1TkMu_muRefHwDXY);
  ntuple_->Branch("L1TkMu_muRefHwEta", &mL1TkMu_muRefHwEta);
  ntuple_->Branch("L1TkMu_muRefHwPhi", &mL1TkMu_muRefHwPhi);
  ntuple_->Branch("L1TkMu_muRefHwSign", &mL1TkMu_muRefHwSign);
  ntuple_->Branch("L1TkMu_muRefHwSignValid", &mL1TkMu_muRefHwSignValid);
  ntuple_->Branch("L1TkMu_muRefHwQual", &mL1TkMu_muRefHwQual);
  }

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

  ntuple_->Branch("nL3Muon", &nL3Muon_, "nL3Muon/I");
  ntuple_->Branch("L3Muon_pt", &L3Muon_pt_, "L3Muon_pt[nL3Muon]/D");
  ntuple_->Branch("L3Muon_eta", &L3Muon_eta_, "L3Muon_eta[nL3Muon]/D");
  ntuple_->Branch("L3Muon_phi", &L3Muon_phi_, "L3Muon_phi[nL3Muon]/D");
  ntuple_->Branch("L3Muon_charge", &L3Muon_charge_, "L3Muon_charge[nL3Muon]/D");
  ntuple_->Branch("L3Muon_trkPt", &L3Muon_trkPt_, "L3Muon_trkPt[nL3Muon]/D");

  ntuple_->Branch("nL2Muon", &nL2Muon_, "nL2Muon/I");
  ntuple_->Branch("L2Muon_pt", &L2Muon_pt_, "L2Muon_pt[nL2Muon]/D");
  ntuple_->Branch("L2Muon_eta", &L2Muon_eta_, "L2Muon_eta[nL2Muon]/D");
  ntuple_->Branch("L2Muon_phi", &L2Muon_phi_, "L2Muon_phi[nL2Muon]/D");
  ntuple_->Branch("L2Muon_charge", &L2Muon_charge_, "L2Muon_charge[nL2Muon]/D");
  ntuple_->Branch("L2Muon_trkPt", &L2Muon_trkPt_, "L2Muon_trkPt[nL2Muon]/D");

  ntuple_->Branch("nTkMuon", &nTkMuon_, "nTkMuon/I");
  ntuple_->Branch("TkMuon_pt", &TkMuon_pt_, "TkMuon_pt[nTkMuon]/D");
  ntuple_->Branch("TkMuon_eta", &TkMuon_eta_, "TkMuon_eta[nTkMuon]/D");
  ntuple_->Branch("TkMuon_phi", &TkMuon_phi_, "TkMuon_phi[nTkMuon]/D");
  ntuple_->Branch("TkMuon_charge", &TkMuon_charge_, "TkMuon_charge[nTkMuon]/D");
  ntuple_->Branch("TkMuon_trkPt", &TkMuon_trkPt_, "TkMuon_trkPt[nTkMuon]/D");

  ntuple_->Branch("nL1Muon", &nL1Muon_, "nL1Muon/I");
  ntuple_->Branch("L1Muon_pt", &L1Muon_pt_, "L1Muon_pt[nL1Muon]/D");
  ntuple_->Branch("L1Muon_eta", &L1Muon_eta_, "L1Muon_eta[nL1Muon]/D");
  ntuple_->Branch("L1Muon_phi", &L1Muon_phi_, "L1Muon_phi[nL1Muon]/D");
  ntuple_->Branch("L1Muon_charge", &L1Muon_charge_, "L1Muon_charge[nL1Muon]/D");
  ntuple_->Branch("L1Muon_quality", &L1Muon_quality_, "L1Muon_quality[nL1Muon]/D");
  ntuple_->Branch("L1Muon_etaAtVtx", &L1Muon_etaAtVtx_, "L1Muon_etaAtVtx[nL1Muon]/D");
  ntuple_->Branch("L1Muon_phiAtVtx", &L1Muon_phiAtVtx_, "L1Muon_phiAtVtx[nL1Muon]/D");

  ntuple_->Branch("nIterL3OI", &nIterL3OI_, "nIterL3OI/I");
  ntuple_->Branch("iterL3OI_inner_pt", &iterL3OI_inner_pt_, "iterL3OI_inner_pt[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_inner_eta", &iterL3OI_inner_eta_, "iterL3OI_inner_eta[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_inner_phi", &iterL3OI_inner_phi_, "iterL3OI_inner_phi[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_inner_charge", &iterL3OI_inner_charge_, "iterL3OI_inner_charge[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_outer_pt", &iterL3OI_outer_pt_, "iterL3OI_outer_pt[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_outer_eta", &iterL3OI_outer_eta_, "iterL3OI_outer_eta[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_outer_phi", &iterL3OI_outer_phi_, "iterL3OI_outer_phi[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_outer_charge", &iterL3OI_outer_charge_, "iterL3OI_outer_charge[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_global_pt", &iterL3OI_global_pt_, "iterL3OI_global_pt[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_global_eta", &iterL3OI_global_eta_, "iterL3OI_global_eta[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_global_phi", &iterL3OI_global_phi_, "iterL3OI_global_phi[nIterL3OI]/D");
  ntuple_->Branch("iterL3OI_global_charge", &iterL3OI_global_charge_, "iterL3OI_global_charge[nIterL3OI]/D");

  ntuple_->Branch("nIterL3IOFromL2", &nIterL3IOFromL2_, "nIterL3IOFromL2/I");
  ntuple_->Branch("iterL3IOFromL2_inner_pt", &iterL3IOFromL2_inner_pt_, "iterL3IOFromL2_inner_pt[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_inner_eta", &iterL3IOFromL2_inner_eta_, "iterL3IOFromL2_inner_eta[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_inner_phi", &iterL3IOFromL2_inner_phi_, "iterL3IOFromL2_inner_phi[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_inner_charge", &iterL3IOFromL2_inner_charge_, "iterL3IOFromL2_inner_charge[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_outer_pt", &iterL3IOFromL2_outer_pt_, "iterL3IOFromL2_outer_pt[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_outer_eta", &iterL3IOFromL2_outer_eta_, "iterL3IOFromL2_outer_eta[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_outer_phi", &iterL3IOFromL2_outer_phi_, "iterL3IOFromL2_outer_phi[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_outer_charge", &iterL3IOFromL2_outer_charge_, "iterL3IOFromL2_outer_charge[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_global_pt", &iterL3IOFromL2_global_pt_, "iterL3IOFromL2_global_pt[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_global_eta", &iterL3IOFromL2_global_eta_, "iterL3IOFromL2_global_eta[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_global_phi", &iterL3IOFromL2_global_phi_, "iterL3IOFromL2_global_phi[nIterL3IOFromL2]/D");
  ntuple_->Branch("iterL3IOFromL2_global_charge", &iterL3IOFromL2_global_charge_, "iterL3IOFromL2_global_charge[nIterL3IOFromL2]/D");

  ntuple_->Branch("nIterL3IOFromL1", &nIterL3IOFromL1_, "nIterL3IOFromL1/I");
  ntuple_->Branch("iterL3IOFromL1_pt", &iterL3IOFromL1_pt_, "iterL3IOFromL1_pt[nIterL3IOFromL1]/D");
  ntuple_->Branch("iterL3IOFromL1_eta", &iterL3IOFromL1_eta_, "iterL3IOFromL1_eta[nIterL3IOFromL1]/D");
  ntuple_->Branch("iterL3IOFromL1_phi", &iterL3IOFromL1_phi_, "iterL3IOFromL1_phi[nIterL3IOFromL1]/D");
  ntuple_->Branch("iterL3IOFromL1_charge", &iterL3IOFromL1_charge_, "iterL3IOFromL1_charge[nIterL3IOFromL1]/D");

  ntuple_->Branch("nIterL3FromL2", &nIterL3FromL2_, "nIterL3FromL2/I");
  ntuple_->Branch("iterL3FromL2_inner_pt", &iterL3FromL2_inner_pt_, "iterL3FromL2_inner_pt[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_inner_eta", &iterL3FromL2_inner_eta_, "iterL3FromL2_inner_eta[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_inner_phi", &iterL3FromL2_inner_phi_, "iterL3FromL2_inner_phi[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_inner_charge", &iterL3FromL2_inner_charge_, "iterL3FromL2_inner_charge[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_outer_pt", &iterL3FromL2_outer_pt_, "iterL3FromL2_outer_pt[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_outer_eta", &iterL3FromL2_outer_eta_, "iterL3FromL2_outer_eta[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_outer_phi", &iterL3FromL2_outer_phi_, "iterL3FromL2_outer_phi[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_outer_charge", &iterL3FromL2_outer_charge_, "iterL3FromL2_outer_charge[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_global_pt", &iterL3FromL2_global_pt_, "iterL3FromL2_global_pt[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_global_eta", &iterL3FromL2_global_eta_, "iterL3FromL2_global_eta[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_global_phi", &iterL3FromL2_global_phi_, "iterL3FromL2_global_phi[nIterL3FromL2]/D");
  ntuple_->Branch("iterL3FromL2_global_charge", &iterL3FromL2_global_charge_, "iterL3FromL2_global_charge[nIterL3FromL2]/D");

  ntuple_->Branch("nIterL3MuonNoID",       &nIterL3MuonNoID_,       "nIterL3MuonNoID/I");
  ntuple_->Branch("iterL3MuonNoID_pt",     &iterL3MuonNoID_pt_,     "iterL3MuonNoID_pt[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_innerPt",     &iterL3MuonNoID_innerPt_,     "iterL3MuonNoID_innerPt[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_eta",    &iterL3MuonNoID_eta_,    "iterL3MuonNoID_eta[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_phi",    &iterL3MuonNoID_phi_,    "iterL3MuonNoID_phi[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_charge", &iterL3MuonNoID_charge_, "iterL3MuonNoID_charge[nIterL3MuonNoID]/D");
  ntuple_->Branch("iterL3MuonNoID_isGLB",  &iterL3MuonNoID_isGLB_,  "iterL3MuonNoID_isGLB[nIterL3MuonNoID]/I");
  ntuple_->Branch("iterL3MuonNoID_isSTA",  &iterL3MuonNoID_isSTA_,  "iterL3MuonNoID_isSTA[nIterL3MuonNoID]/I");
  ntuple_->Branch("iterL3MuonNoID_isTRK",  &iterL3MuonNoID_isTRK_,  "iterL3MuonNoID_isTRK[nIterL3MuonNoID]/I");

  ntuple_->Branch("nIterL3Muon",       &nIterL3Muon_,       "nIterL3Muon/I");
  ntuple_->Branch("iterL3Muon_pt",     &iterL3Muon_pt_,     "iterL3Muon_pt[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_innerPt", &iterL3Muon_innerPt_, "iterL3Muon_innerPt[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_eta",    &iterL3Muon_eta_,    "iterL3Muon_eta[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_phi",    &iterL3Muon_phi_,    "iterL3Muon_phi[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_charge", &iterL3Muon_charge_, "iterL3Muon_charge[nIterL3Muon]/D");
  ntuple_->Branch("iterL3Muon_isGLB",  &iterL3Muon_isGLB_,  "iterL3Muon_isGLB[nIterL3Muon]/I");
  ntuple_->Branch("iterL3Muon_isSTA",  &iterL3Muon_isSTA_,  "iterL3Muon_isSTA[nIterL3Muon]/I");
  ntuple_->Branch("iterL3Muon_isTRK",  &iterL3Muon_isTRK_,  "iterL3Muon_isTRK[nIterL3Muon]/I");

  ntuple_->Branch("nMuonCand", &nMuonCand_, "nMuonCand/I");

  ntuple_->Branch("muonCand_pt", &muonCand_pt_, "muonCand_pt[nMuonCand]/D");
  ntuple_->Branch("muonCand_eta", &muonCand_eta_, "muonCand_eta[nMuonCand]/D");
  ntuple_->Branch("muonCand_phi", &muonCand_phi_, "muonCand_phi[nMuonCand]/D");
  ntuple_->Branch("muonCand_px", &muonCand_px_, "muonCand_px[nMuonCand]/D");
  ntuple_->Branch("muonCand_py", &muonCand_py_, "muonCand_py[nMuonCand]/D");
  ntuple_->Branch("muonCand_pz", &muonCand_pz_, "muonCand_pz[nMuonCand]/D");
  ntuple_->Branch("muonCand_vx", &muonCand_vx_, "muonCand_vx[nMuonCand]/D");
  ntuple_->Branch("muonCand_vy", &muonCand_vy_, "muonCand_vy[nMuonCand]/D");
  ntuple_->Branch("muonCand_vz", &muonCand_vz_, "muonCand_vz[nMuonCand]/D");
  ntuple_->Branch("muonCand_dz", &muonCand_dz_, "muonCand_dz[nMuonCand]/D");
  ntuple_->Branch("muonCand_dxy", &muonCand_dxy_, "muonCand_dxy[nMuonCand]/D");
  ntuple_->Branch("muonCand_dB", &muonCand_dB_, "muonCand_dB[nMuonCand]/D");
  ntuple_->Branch("muonCand_charge", &muonCand_charge_, "muonCand_charge[nMuonCand]/D");
  ntuple_->Branch("muonCand_isGLB", &muonCand_isGLB_, "muonCand_isGLB[nMuonCand]/I");
  ntuple_->Branch("muonCand_isSTA", &muonCand_isSTA_, "muonCand_isSTA[nMuonCand]/I");
  ntuple_->Branch("muonCand_isTRK", &muonCand_isTRK_, "muonCand_isTRK[nMuonCand]/I");
  ntuple_->Branch("muonCand_isMuon", &muonCand_isMuon_, "muonCand_isMuon[nMuonCand]/I");

  ntuple_->Branch("muonCand_inner_trkChi2", &muonCand_inner_trkChi2_, "muonCand_inner_trkChi2[nMuonCand]/D");
  ntuple_->Branch("muonCand_inner_validFraction", &muonCand_inner_validFraction_, "muonCand_inner_validFraction[nMuonCand]/D");
  ntuple_->Branch("muonCand_inner_trackerLayers", &muonCand_inner_trackerLayers_, "muonCand_inner_trackerLayers[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_trackerHits", &muonCand_inner_trackerHits_, "muonCand_inner_trackerHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostTrackerHits", &muonCand_inner_lostTrackerHits_, "muonCand_inner_lostTrackerHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostTrackerHitsIn", &muonCand_inner_lostTrackerHitsIn_, "muonCand_inner_lostTrackerHitsIn[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostTrackerHitsOut", &muonCand_inner_lostTrackerHitsOut_, "muonCand_inner_lostTrackerHitsOut[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostPixelHits", &muonCand_inner_lostPixelHits_, "muonCand_inner_lostPixelHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostPixelBarrelHits", &muonCand_inner_lostPixelBarrelHits_, "muonCand_inner_lostPixelBarrelHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostPixelEndcapHits", &muonCand_inner_lostPixelEndcapHits_, "muonCand_inner_lostPixelEndcapHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripHits", &muonCand_inner_lostStripHits_, "muonCand_inner_lostStripHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripTIBHits", &muonCand_inner_lostStripTIBHits_, "muonCand_inner_lostStripTIBHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripTIDHits", &muonCand_inner_lostStripTIDHits_, "muonCand_inner_lostStripTIDHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripTOBHits", &muonCand_inner_lostStripTOBHits_, "muonCand_inner_lostStripTOBHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_lostStripTECHits", &muonCand_inner_lostStripTECHits_, "muonCand_inner_lostStripTECHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_pixelLayers", &muonCand_inner_pixelLayers_, "muonCand_inner_pixelLayers[nMuonCand]/I");
  ntuple_->Branch("muonCand_inner_pixelHits", &muonCand_inner_pixelHits_, "muonCand_inner_pixelHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_global_muonHits", &muonCand_global_muonHits_, "muonCand_global_muonHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_global_trkChi2", &muonCand_global_trkChi2_, "muonCand_global_trkChi2[nMuonCand]/D");
  ntuple_->Branch("muonCand_global_trackerLayers", &muonCand_global_trackerLayers_, "muonCand_global_trackerLayers[nMuonCand]/I");
  ntuple_->Branch("muonCand_global_trackerHits", &muonCand_global_trackerHits_, "muonCand_global_trackerHits[nMuonCand]/I");
  ntuple_->Branch("muonCand_btlMatchChi2", &muonCand_btlMatchChi2_, "muonCand_btlMatchChi2[nMuonCand]/D");
  ntuple_->Branch("muonCand_etlMatchChi2", &muonCand_etlMatchChi2_, "muonCand_etlMatchChi2[nMuonCand]/D");
  ntuple_->Branch("muonCand_btlMatchTimeChi2", &muonCand_btlMatchTimeChi2_, "muonCand_btlMatchTimeChi2[nMuonCand]/D");
  ntuple_->Branch("muonCand_etlMatchTimeChi2", &muonCand_etlMatchTimeChi2_, "muonCand_etlMatchTimeChi2[nMuonCand]/D");
  ntuple_->Branch("muonCand_npixBarrel", &muonCand_npixBarrel_, "muonCand_npixBarrel[nMuonCand]/I");
  ntuple_->Branch("muonCand_npixEndcap", &muonCand_npixEndcap_, "muonCand_npixEndcap[nMuonCand]/I");
  ntuple_->Branch("muonCand_outermostHitPosition", &muonCand_outermostHitPosition_, "muonCand_outermostHitPosition[nMuonCand]/D");
  ntuple_->Branch("muonCand_p", &muonCand_p_, "muonCand_p[nMuonCand]/D");
  ntuple_->Branch("muonCand_beta", &muonCand_beta_, "muonCand_beta[nMuonCand]/D");
  ntuple_->Branch("muonCand_t0", &muonCand_t0_, "muonCand_t0[nMuonCand]/D");
  ntuple_->Branch("muonCand_sigmat0", &muonCand_sigmat0_, "muonCand_sigmat0[nMuonCand]/D");
  ntuple_->Branch("muonCand_pathLength", &muonCand_pathLength_, "muonCand_pathLength[nMuonCand]/D");
  ntuple_->Branch("muonCand_tmtd", &muonCand_tmtd_, "muonCand_tmtd[nMuonCand]/D");
  ntuple_->Branch("muonCand_sigmatmtd", &muonCand_sigmatmtd_, "muonCand_sigmatmtd[nMuonCand]/D");
  ntuple_->Branch("muonCand_tofMu", &muonCand_tofMu_, "muonCand_tofMu[nMuonCand]/D");
  ntuple_->Branch("muonCand_sigmaTofMu", &muonCand_sigmaTofMu_, "muonCand_sigmaTofMu[nMuonCand]/D");
  ntuple_->Branch("muonCand_mtdpos_x", &muonCand_mtdpos_x_, "muonCand_mtdpos_x[nMuonCand]/D");
  ntuple_->Branch("muonCand_mtdpos_y", &muonCand_mtdpos_y_, "muonCand_mtdpos_y[nMuonCand]/D");
  ntuple_->Branch("muonCand_mtdpos_z", &muonCand_mtdpos_z_, "muonCand_mtdpos_z[nMuonCand]/D");
  
  // ECAL Isolation ----------
  ntuple_->Branch("nECAL", &nECAL_, "nECAL/I");
  ntuple_->Branch("ecal_et", &ecal_et_, "ecal_et_[nECAL]/D");
  ntuple_->Branch("ecal_pt", &ecal_pt_, "ecal_pt_[nECAL]/D");
  ntuple_->Branch("ecal_eta", &ecal_eta_, "ecal_eta_[nECAL]/D");
  ntuple_->Branch("ecal_phi", &ecal_phi_, "ecal_phi_[nECAL]/D");
  ntuple_->Branch("ecal_charge", &ecal_charge_, "ecal_charge_[nECAL]/D");
  ntuple_->Branch("ecal_px", &ecal_px_, "ecal_px_[nECAL]/D");
  ntuple_->Branch("ecal_py", &ecal_py_, "ecal_py_[nECAL]/D");
  ntuple_->Branch("ecal_pz", &ecal_pz_, "ecal_pz_[nECAL]/D");
  ntuple_->Branch("ecal_vx", &ecal_vx_, "ecal_vx_[nECAL]/D");
  ntuple_->Branch("ecal_vy", &ecal_vy_, "ecal_vy_[nECAL]/D");
  ntuple_->Branch("ecal_vz", &ecal_vz_, "ecal_vz_[nECAL]/D");
  ntuple_->Branch("ecal_time", &ecal_time_, "ecal_time_[nECAL]/D");
  ntuple_->Branch("ecal_timeErr", &ecal_timeErr_, "ecal_timeErr_[nECAL]/D");
  ntuple_->Branch("ecal_depth", &ecal_depth_, "ecal_depth_[nECAL]/D");
  ntuple_->Branch("ecal_rho", &ecal_rho_, "ecal_rho_[nECAL]/D");
  ntuple_->Branch("ecal_muonIdx", &ecal_muonIdx_, "ecal_muonIdx_[nECAL]/I");
  ntuple_->Branch("ecal_nHits", &ecal_nHits_, "ecal_nHits_[nECAL]/I");
  
  ntuple_->Branch("nECALHits", &nECALHits_, "nECALHits/I");
  ntuple_->Branch("ecal_hit_energy", &ecal_hit_energy_, "ecal_hit_energy_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_depth", &ecal_hit_depth_, "ecal_hit_depth_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_time", &ecal_hit_time_, "ecal_hit_time_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_timeErr", &ecal_hit_timeErr_, "ecal_hit_timeErr_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_pt2", &ecal_hit_pt2_, "ecal_hit_pt2_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_x", &ecal_hit_x_, "ecal_hit_x_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_y", &ecal_hit_y_, "ecal_hit_y_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_z", &ecal_hit_z_, "ecal_hit_z_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_eta", &ecal_hit_eta_, "ecal_hit_eta_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_phi", &ecal_hit_phi_, "ecal_hit_phi_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_fraction", &ecal_hit_fraction_, "ecal_hit_fraction_[nECALHits]/D");
  ntuple_->Branch("ecal_hit_idx", &ecal_hit_idx_, "ecal_hit_idx_[nECALHits]/D");

  ntuple_->Branch("nHCAL", &nHCAL_, "nHCAL/I");
  ntuple_->Branch("hcal_et", &hcal_et_, "hcal_et_[nHCAL]/D");
  ntuple_->Branch("hcal_pt", &hcal_pt_, "hcal_pt_[nHCAL]/D");
  ntuple_->Branch("hcal_eta", &hcal_eta_, "hcal_eta_[nHCAL]/D");
  ntuple_->Branch("hcal_phi", &hcal_phi_, "hcal_phi_[nHCAL]/D");
  ntuple_->Branch("hcal_charge", &hcal_charge_, "hcal_charge_[nHCAL]/D");
  ntuple_->Branch("hcal_px", &hcal_px_, "hcal_px_[nHCAL]/D");
  ntuple_->Branch("hcal_py", &hcal_py_, "hcal_py_[nHCAL]/D");
  ntuple_->Branch("hcal_pz", &hcal_pz_, "hcal_pz_[nHCAL]/D");
  ntuple_->Branch("hcal_vx", &hcal_vx_, "hcal_vx_[nHCAL]/D");
  ntuple_->Branch("hcal_vy", &hcal_vy_, "hcal_vy_[nHCAL]/D");
  ntuple_->Branch("hcal_vz", &hcal_vz_, "hcal_vz_[nHCAL]/D");
  ntuple_->Branch("hcal_time", &hcal_time_, "hcal_time_[nHCAL]/D");
  ntuple_->Branch("hcal_timeErr", &hcal_timeErr_, "hcal_timeErr_[nHCAL]/D");
  ntuple_->Branch("hcal_depth", &hcal_depth_, "hcal_depth_[nHCAL]/D");
  ntuple_->Branch("hcal_rho", &hcal_rho_, "hcal_rho_[nHCAL]/D");
  ntuple_->Branch("hcal_muonIdx", &hcal_muonIdx_, "hcal_muonIdx_[nHCAL]/I");
  ntuple_->Branch("hcal_nHits", &hcal_nHits_, "hcal_nHits_[nHCAL]/I");

  ntuple_->Branch("nHCALHits", &nHCALHits_, "nHCALHits/I");
  ntuple_->Branch("hcal_hit_energy", &hcal_hit_energy_, "hcal_hit_energy_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_depth", &hcal_hit_depth_, "hcal_hit_depth_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_time", &hcal_hit_time_, "hcal_hit_time_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_timeErr", &hcal_hit_timeErr_, "hcal_hit_timeErr_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_pt2", &hcal_hit_pt2_, "hcal_hit_pt2_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_x", &hcal_hit_x_, "hcal_hit_x_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_y", &hcal_hit_y_, "hcal_hit_y_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_z", &hcal_hit_z_, "hcal_hit_z_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_eta", &hcal_hit_eta_, "hcal_hit_eta_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_phi", &hcal_hit_phi_, "hcal_hit_phi_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_fraction", &hcal_hit_fraction_, "hcal_hit_fraction_[nHCALHits]/D");
  ntuple_->Branch("hcal_hit_idx", &hcal_hit_idx_, "hcal_hit_idx_[nHCALHits]/D");

  ntuple_->Branch("nHGCAL_em", &nHGCAL_em_, "nHGCAL_em/I");
  ntuple_->Branch("hgcal_em_et", &hgcal_em_et_, "hgcal_em_et/D");
  ntuple_->Branch("hgcal_em_pt", &hgcal_em_pt_, "hgcal_em_pt/D");
  ntuple_->Branch("hgcal_em_eta", &hgcal_em_eta_, "hgcal_em_eta/D");
  ntuple_->Branch("hgcal_em_phi", &hgcal_em_phi_, "hgcal_em_phi/D");
  ntuple_->Branch("hgcal_em_px", &hgcal_em_px_, "hgcal_em_px/D");
  ntuple_->Branch("hgcal_em_py", &hgcal_em_py_, "hgcal_em_py/D");
  ntuple_->Branch("hgcal_em_pz", &hgcal_em_pz_, "hgcal_em_pz/D");
  ntuple_->Branch("hgcal_em_vx", &hgcal_em_vx_, "hgcal_em_vx/D");
  ntuple_->Branch("hgcal_em_vy", &hgcal_em_vy_, "hgcal_em_vy/D");
  ntuple_->Branch("hgcal_em_vz", &hgcal_em_vz_, "hgcal_em_vz/D");
  ntuple_->Branch("hgcal_em_time", &hgcal_em_time_, "hgcal_em_time/D");
  ntuple_->Branch("hgcal_em_timeErr", &hgcal_em_timeErr_, "hgcal_em_timeErr/D");
  ntuple_->Branch("hgcal_em_depth", &hgcal_em_depth_, "hgcal_em_depth/D");
  ntuple_->Branch("hgcal_em_algoID", &hgcal_em_algoID_, "hgcal_em_algoID/D");
  ntuple_->Branch("hgcal_em_muonIdx", &hgcal_em_muonIdx_, "hgcal_em_muonIdx/I");
  
  ntuple_->Branch("nHGCAL_had", &nHGCAL_had_, "nHGCAL_had/I");
  ntuple_->Branch("hgcal_had_et", &hgcal_had_et_, "hgcal_had_et/D");
  ntuple_->Branch("hgcal_had_pt", &hgcal_had_pt_, "hgcal_had_pt/D");
  ntuple_->Branch("hgcal_had_eta", &hgcal_had_eta_, "hgcal_had_eta/D");
  ntuple_->Branch("hgcal_had_phi", &hgcal_had_phi_, "hgcal_had_phi/D");
  ntuple_->Branch("hgcal_had_px", &hgcal_had_px_, "hgcal_had_px/D");
  ntuple_->Branch("hgcal_had_py", &hgcal_had_py_, "hgcal_had_py/D");
  ntuple_->Branch("hgcal_had_pz", &hgcal_had_pz_, "hgcal_had_pz/D");
  ntuple_->Branch("hgcal_had_vx", &hgcal_had_vx_, "hgcal_had_vx/D");
  ntuple_->Branch("hgcal_had_vy", &hgcal_had_vy_, "hgcal_had_vy/D");
  ntuple_->Branch("hgcal_had_vz", &hgcal_had_vz_, "hgcal_had_vz/D");
  ntuple_->Branch("hgcal_had_time", &hgcal_had_time_, "hgcal_had_time/D");
  ntuple_->Branch("hgcal_had_timeErr", &hgcal_had_timeErr_, "hgcal_had_timeErr/D");
  ntuple_->Branch("hgcal_had_depth", &hgcal_had_depth_, "hgcal_had_depth/D");
  ntuple_->Branch("hgcal_had_algoID", &hgcal_had_algoID_, "hgcal_had_algoID/D");
  ntuple_->Branch("hgcal_had_muonIdx", &hgcal_had_muonIdx_, "hgcal_had_muonIdx/I");

  ntuple_->Branch("nTrack", &nTrack_, "nTrack/I");
  ntuple_->Branch("track_pt", &track_pt_, "track_pt[nTrack]/D");
  ntuple_->Branch("track_eta", &track_eta_, "track_eta[nTrack]/D");
  ntuple_->Branch("track_phi", &track_phi_, "track_phi[nTrack]/D");
  ntuple_->Branch("track_charge", &track_charge_, "track_charge[nTrack]/D");
  ntuple_->Branch("track_px", &track_px_, "track_px[nTrack]/D");
  ntuple_->Branch("track_py", &track_py_, "track_py[nTrack]/D");
  ntuple_->Branch("track_pz", &track_pz_, "track_pz[nTrack]/D");
  ntuple_->Branch("track_vx", &track_vx_, "track_vx[nTrack]/D");
  ntuple_->Branch("track_vy", &track_vy_, "track_vy[nTrack]/D");
  ntuple_->Branch("track_vz", &track_vz_, "track_vz[nTrack]/D");
  ntuple_->Branch("track_dxy_bs", &track_dxy_bs_, "track_dxy_bs[nTrack]/D");
  ntuple_->Branch("track_dxyError_bs", &track_dxyError_bs_, "track_dxyError_bs[nTrack]/D");
  ntuple_->Branch("track_dz_bs", &track_dz_bs_, "track_dz_bs[nTrack]/D");
  ntuple_->Branch("track_dzError", &track_dzError_, "track_dzError[nTrack]/D");
  ntuple_->Branch("track_trkChi2", &track_trkChi2_, "track_trkChi2[nTrack]/D");
  ntuple_->Branch("track_trackerLayers", &track_trackerLayers_, "track_trackerLayers[nTrack]/I");
  ntuple_->Branch("track_trackerHits", &track_trackerHits_, "track_trackerHits[nTrack]/I");
  ntuple_->Branch("track_lostTrackerHits", &track_lostTrackerHits_, "track_lostTrackerHits[nTrack]/I");
  ntuple_->Branch("track_lostTrackerHitsIn", &track_lostTrackerHitsIn_, "track_lostTrackerHitsIn[nTrack]/I");
  ntuple_->Branch("track_lostTrackerHitsOut", &track_lostTrackerHitsOut_, "track_lostTrackerHitsOut[nTrack]/I");
  ntuple_->Branch("track_lostPixelHits", &track_lostPixelHits_, "track_lostPixelHits[nTrack]/I");
  ntuple_->Branch("track_lostPixelBarrelHits", &track_lostPixelBarrelHits_, "track_lostPixelBarrelHits[nTrack]/I");
  ntuple_->Branch("track_lostPixelEndcapHits", &track_lostPixelEndcapHits_, "track_lostPixelEndcapHits[nTrack]/I");
  ntuple_->Branch("track_lostStripHits", &track_lostStripHits_, "track_lostStripHits[nTrack]/I");
  ntuple_->Branch("track_lostStripTIBHits", &track_lostStripTIBHits_, "track_lostStripTIBHits[nTrack]/I");
  ntuple_->Branch("track_lostStripTIDHits", &track_lostStripTIDHits_, "track_lostStripTIDHits[nTrack]/I");
  ntuple_->Branch("track_lostStripTOBHits", &track_lostStripTOBHits_, "track_lostStripTOBHits[nTrack]/I");
  ntuple_->Branch("track_lostStripTECHits", &track_lostStripTECHits_, "track_lostStripTECHits[nTrack]/I");
  ntuple_->Branch("track_pixelLayers", &track_pixelLayers_, "track_pixelLayers[nTrack]/I");
  ntuple_->Branch("track_pixelHits", &track_pixelHits_, "track_pixelHits[nTrack]/I");
  ntuple_->Branch("track_muonHits", &track_muonHits_, "track_muonHits[nTrack]/I");
  ntuple_->Branch("track_muonIdx", &track_muonIdx_, "track_muonIdx[nTrack]/I");
  ntuple_->Branch("track_t0Src", &track_t0Src_, "track_t0Src[nTrack]/D");
  ntuple_->Branch("track_Sigmat0Src", &track_Sigmat0Src_, "track_Sigmat0Src[nTrack]/D");
  ntuple_->Branch("track_t0Pid", &track_t0Pid_, "track_t0Pid[nTrack]/D");
  ntuple_->Branch("track_t0Safe", &track_t0Safe_, "track_t0Safe[nTrack]/D");
  ntuple_->Branch("track_sigmat0Safe", &track_sigmat0Safe_, "track_sigmat0Safe[nTrack]/D");
  ntuple_->Branch("track_mtdQualMVA", &track_mtdQualMVA_, "track_mtdQualMVA[nTrack]/D");
  ntuple_->Branch("track_tMtd", &track_tMtd_, "track_tMtd[nTrack]/D");
  ntuple_->Branch("track_tofPi", &track_tofPi_, "track_tofPi[nTrack]/D");
  ntuple_->Branch("track_tofK", &track_tofK_, "track_tofK[nTrack]/D");
  ntuple_->Branch("track_tofP", &track_tofP_, "track_tofP[nTrack]/D");
  ntuple_->Branch("track_probPi", &track_probPi_, "track_probPi[nTrack]/D");
  ntuple_->Branch("track_probK", &track_probK_, "track_probK[nTrack]/D");
  ntuple_->Branch("track_probP", &track_probP_, "track_probP[nTrack]/D");
  ntuple_->Branch("track_sigmatofpi", &track_sigmatofpi_, "track_sigmatofpi[nTrack]/D");
  ntuple_->Branch("track_sigmatofk", &track_sigmatofk_, "track_sigmatofk[nTrack]/D");
  ntuple_->Branch("track_sigmatofp", &track_sigmatofp_, "track_sigmatofp[nTrack]/D");
  ntuple_->Branch("track_btlMatchChi2", &track_btlMatchChi2_, "track_btlMatchChi2[nTrack]/D");
  ntuple_->Branch("track_btlMatchTimeChi2", &track_btlMatchTimeChi2_, "track_btlMatchTimeChi2[nTrack]/D");
  ntuple_->Branch("track_etlMatchChi2", &track_etlMatchChi2_, "track_etlMatchChi2[nTrack]/D");
  ntuple_->Branch("track_etlMatchTimeChi2", &track_etlMatchTimeChi2_, "track_etlMatchTimeChi2[nTrack]/D");
  ntuple_->Branch("track_npixBarrel", &track_npixBarrel_, "track_npixBarrel[nTrack]/I");
  ntuple_->Branch("track_npixEndcap", &track_npixEndcap_, "track_npixEndcap[nTrack]/I");
  ntuple_->Branch("track_outermostHitPosition", &track_outermostHitPosition_, "track_outermostHitPosition[nTrack]/D");
  ntuple_->Branch("track_p", &track_p_, "track_p[nTrack]/D");
  ntuple_->Branch("track_beta", &track_beta_, "track_beta[nTrack]/D");
  ntuple_->Branch("track_pathLength", &track_pathLength_, "track_pathLength[nTrack]/D");
  ntuple_->Branch("track_mtdpos_x", &track_mtdpos_x_, "track_mtdpos_x[nTrack]/D");
  ntuple_->Branch("track_mtdpos_y", &track_mtdpos_y_, "track_mtdpos_y[nTrack]/D");
  ntuple_->Branch("track_mtdpos_z", &track_mtdpos_z_, "track_mtdpos_z[nTrack]/D");
  ntuple_->Branch("track_TPcharge", &track_TPcharge_, "track_TPcharge[nTrack]/I");
  ntuple_->Branch("track_TPpdgId", &track_TPpdgId_, "track_TPpdgId[nTrack]/I");
  ntuple_->Branch("track_TPenergy", &track_TPenergy_, "track_TPenergy[nTrack]/D");
  ntuple_->Branch("track_TPpt", &track_TPpt_, "track_TPpt[nTrack]/D");
  ntuple_->Branch("track_TPeta", &track_TPeta_, "track_TPeta[nTrack]/D");
  ntuple_->Branch("track_TPphi", &track_TPphi_, "track_TPphi[nTrack]/D");
  ntuple_->Branch("track_TPparentVx", &track_TPparentVx_, "track_TPparentVx[nTrack]/D");
  ntuple_->Branch("track_TPparentVy", &track_TPparentVy_, "track_TPparentVy[nTrack]/D");
  ntuple_->Branch("track_TPparentVz", &track_TPparentVz_, "track_TPparentVz[nTrack]/D");
  ntuple_->Branch("track_TPstatus", &track_TPstatus_, "track_TPstatus[nTrack]/I");
  ntuple_->Branch("track_TPnumberOfHits", &track_TPnumberOfHits_, "track_TPnumberOfHits[nTrack]/I");
  ntuple_->Branch("track_TPnumberOfTrackerHits", &track_TPnumberOfTrackerHits_, "track_TPnumberOfTrackerHits[nTrack]/I");
  ntuple_->Branch("track_TPnumberOfTrackerLayers", &track_TPnumberOfTrackerLayers_, "track_TPnumberOfTrackerLayers[nTrack]/I");

  ntuple_->Branch("nPFCand", &nPFCand_, "nPFCand/I");
  ntuple_->Branch("pfcand_pt", &pfcand_pt_, "pfcand_pt[nPFCand]/D");
  ntuple_->Branch("pfcand_eta", &pfcand_eta_, "pfcand_eta[nPFCand]/D");
  ntuple_->Branch("pfcand_phi", &pfcand_phi_, "pfcand_phi[nPFCand]/D");
  ntuple_->Branch("pfcand_charge", &pfcand_charge_, "pfcand_charge[nPFCand]/I");
  ntuple_->Branch("pfcand_pdgId", &pfcand_pdgId_, "pfcand_pdgId[nPFCand]/I");
  ntuple_->Branch("pfcand_px", &pfcand_px_, "pfcand_px[nPFCand]/D");
  ntuple_->Branch("pfcand_py", &pfcand_py_, "pfcand_py[nPFCand]/D");
  ntuple_->Branch("pfcand_pz", &pfcand_pz_, "pfcand_pz[nPFCand]/D");
  ntuple_->Branch("pfcand_vx", &pfcand_vx_, "pfcand_vx[nPFCand]/D");
  ntuple_->Branch("pfcand_vy", &pfcand_vy_, "pfcand_vy[nPFCand]/D");
  ntuple_->Branch("pfcand_vz", &pfcand_vz_, "pfcand_vz[nPFCand]/D");
  ntuple_->Branch("pfcand_time", &pfcand_time_, "pfcand_time[nPFCand]/D");
  ntuple_->Branch("pfcand_timeErr", &pfcand_timeErr_, "pfcand_timeErr[nPFCand]/D");
  ntuple_->Branch("pfcand_dxy", &pfcand_dxy_, "pfcand_dxy[nPFCand]/D");
  ntuple_->Branch("pfcand_dz", &pfcand_dz_, "pfcand_dz[nPFCand]/D");
  ntuple_->Branch("pfcand_dxyErr", &pfcand_dxyErr_, "pfcand_dxyErr[nPFCand]/D");
  ntuple_->Branch("pfcand_dzErr", &pfcand_dzErr_, "pfcand_dzErr[nPFCand]/D");
  ntuple_->Branch("pfcand_vChi2NoF", &pfcand_vChi2NoF_, "pfcand_vChi2NoF[nPFCand]/D");
  ntuple_->Branch("pfcand_muonIdx", &pfcand_muonIdx_, "pfcand_muonIdx[nPFCand]/I");

  ntuple_->Branch("pfcand_rho", &pfcand_rho_, "pfcand_rho/D");
  ntuple_->Branch("pfcand_btlMatchChi2", &pfcand_btlMatchChi2_, "pfcand_btlMatchChi2[nPFCand]/D");
  ntuple_->Branch("pfcand_etlMatchChi2", &pfcand_etlMatchChi2_, "pfcand_etlMatchChi2[nPFCand]/D");
  ntuple_->Branch("pfcand_btlMatchTimeChi2", &pfcand_btlMatchTimeChi2_, "pfcand_btlMatchTimeChi2[nPFCand]/D");
  ntuple_->Branch("pfcand_etlMatchTimeChi2", &pfcand_etlMatchTimeChi2_, "pfcand_etlMatchTimeChi2[nPFCand]/D");
  ntuple_->Branch("pfcand_npixBarrel", &pfcand_npixBarrel_, "pfcand_npixBarrel[nPFCand]/I");
  ntuple_->Branch("pfcand_npixEndcap", &pfcand_npixEndcap_, "pfcand_npixEndcap[nPFCand]/I");
  ntuple_->Branch("pfcand_outermostHitPosition", &pfcand_outermostHitPosition_, "pfcand_outermostHitPosition[nPFCand]/D");
  ntuple_->Branch("pfcand_p", &pfcand_p_, "pfcand_p[nPFCand]/D");
  ntuple_->Branch("pfcand_beta", &pfcand_beta_, "pfcand_beta[nPFCand]/D");
  ntuple_->Branch("pfcand_t0", &pfcand_t0_, "pfcand_t0[nPFCand]/D");
  ntuple_->Branch("pfcand_sigmat0", &pfcand_sigmat0_, "pfcand_sigmat0[nPFCand]/D");
  ntuple_->Branch("pfcand_pathLength", &pfcand_pathLength_, "pfcand_pathLength[nPFCand]/D");
  ntuple_->Branch("pfcand_tmtd", &pfcand_tmtd_, "pfcand_tmtd[nPFCand]/D");
  ntuple_->Branch("pfcand_sigmatmtd", &pfcand_sigmatmtd_, "pfcand_sigmatmtd[nPFCand]/D");
  ntuple_->Branch("pfcand_tofPi", &pfcand_tofPi_, "pfcand_tofPi[nPFCand]/D");
  ntuple_->Branch("pfcand_sigmaTofPi", &pfcand_sigmaTofPi_, "pfcand_sigmaTofPi[nPFCand]/D");
  ntuple_->Branch("pfcand_tofK", &pfcand_tofK_, "pfcand_tofK[nPFCand]/D");
  ntuple_->Branch("pfcand_sigmaTofK", &pfcand_sigmaTofK_, "pfcand_sigmaTofK[nPFCand]/D");
  ntuple_->Branch("pfcand_tofP", &pfcand_tofP_, "pfcand_tofP[nPFCand]/D");
  ntuple_->Branch("pfcand_sigmaTofP", &pfcand_sigmaTofP_, "pfcand_sigmaTofP[nPFCand]/D");
  ntuple_->Branch("pfcand_mtdpos_x", &pfcand_mtdpos_x_, "pfcand_mtdpos_x[nPFCand]/D");
  ntuple_->Branch("pfcand_mtdpos_y", &pfcand_mtdpos_y_, "pfcand_mtdpos_y[nPFCand]/D");
  ntuple_->Branch("pfcand_mtdpos_z", &pfcand_mtdpos_z_, "pfcand_mtdpos_z[nPFCand]/D");
  
  /// End of isolation information
 
  SThltIterL3OISeedsFromL2Muons->setBranch(ntuple_,"hltIterL3OISeedsFromL2Muons");
  SThltIter0IterL3MuonPixelSeedsFromPixelTracks->setBranch(ntuple_,"hltIter0IterL3MuonPixelSeedsFromPixelTracks");
  SThltIter2IterL3MuonPixelSeeds->setBranch(ntuple_,"hltIter2IterL3MuonPixelSeeds");
  SThltIter3IterL3MuonPixelSeeds->setBranch(ntuple_,"hltIter3IterL3MuonPixelSeeds");
  SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->setBranch(ntuple_,"hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks");
  SThltIter2IterL3FromL1MuonPixelSeeds->setBranch(ntuple_,"hltIter2IterL3FromL1MuonPixelSeeds");
  SThltIter3IterL3FromL1MuonPixelSeeds->setBranch(ntuple_,"hltIter3IterL3FromL1MuonPixelSeeds");

  TThltIterL3OIMuonTrack->setBranch(ntuple_,"hltIterL3OIMuonTrack");
  TThltIter0IterL3MuonTrack->setBranch(ntuple_,"hltIter0IterL3MuonTrack");
  TThltIter2IterL3MuonTrack->setBranch(ntuple_,"hltIter2IterL3MuonTrack");
  TThltIter3IterL3MuonTrack->setBranch(ntuple_,"hltIter3IterL3MuonTrack");
  TThltIter0IterL3FromL1MuonTrack->setBranch(ntuple_,"hltIter0IterL3FromL1MuonTrack");
  TThltIter2IterL3FromL1MuonTrack->setBranch(ntuple_,"hltIter2IterL3FromL1MuonTrack");
  TThltIter3IterL3FromL1MuonTrack->setBranch(ntuple_,"hltIter3IterL3FromL1MuonTrack");

  MTL3MuonsNoId->setBranch(ntuple_,"L3MuonsNoId");
  MTL3Muons->setBranch(ntuple_,"L3Muons");

  TrkParticle->setBranch(ntuple_,"TP");

  VThltIterL3MuonTrimmedPixelVertices->setBranch(ntuple_,"hltIterL3MuonTrimmedPixelVertices");
  VThltIterL3FromL1MuonTrimmedPixelVertices->setBranch(ntuple_,"hltIterL3FromL1MuonTrimmedPixelVertices");

  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    TString trkName = TString(trackCollectionNames_.at(i));
    TString tpName  = "tpTo_" + TString(trackCollectionNames_.at(i));

    bool doIso = (i == trackCollectionNames_.size()-1);

    if(doIso)
      trkTemplates_.at(i)->setIsoTags( trkIsoTags_, pfIsoTags_ );

    trkTemplates_.at(i)->setBranch(ntuple_, trkName, doIso );
    tpTemplates_.at(i)->setBranch(ntuple_,  tpName );
  }
}

void MuonHLTNtupler::Fill_L1Track(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > TTTrackHandle;
  iEvent.getByToken(ttTrackToken_, TTTrackHandle);

  if (SaveAllTracks){
  if (DebugMode) {
    cout << endl << "Loop over L1 tracks!" << endl;
    // cout << endl << "Looking at " << L1Tk_nPar << "-parameter tracks!" << endl;
  }

  int this_l1track = 0;
  std::vector< TTTrack< Ref_Phase2TrackerDigi_ > >::const_iterator iterL1Track;
  for ( iterL1Track = TTTrackHandle->begin(); iterL1Track != TTTrackHandle->end(); iterL1Track++ ) {
    edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > l1track_ptr(TTTrackHandle, this_l1track);

    float tmp_trk_pt   = iterL1Track->momentum().perp();
    float tmp_trk_eta  = iterL1Track->momentum().eta();
    float tmp_trk_phi  = iterL1Track->momentum().phi();
    float tmp_trk_z0   = iterL1Track->z0(); //cm

    float tmp_trk_rInv = static_cast<float>(iterL1Track->rInv());
    float tmp_trk_tanL = static_cast<float>(iterL1Track->tanL());
    float tmp_trk_MVA1 = static_cast<float>(iterL1Track->trkMVA1());
    float tmp_trk_MVA2 = static_cast<float>(iterL1Track->trkMVA2());
    float tmp_trk_MVA3 = static_cast<float>(iterL1Track->trkMVA3());

    float tmp_trk_d0 = static_cast<float>(iterL1Track->d0());

    // if (L1Tk_nPar == 5) {
    //   float tmp_trk_x0   = iterL1Track->POCA().x();
    //   float tmp_trk_y0   = iterL1Track->POCA().y();
    //   tmp_trk_d0 = -tmp_trk_x0*sin(tmp_trk_phi) + tmp_trk_y0*cos(tmp_trk_phi);
    //   }

    float tmp_trk_chi2 = iterL1Track->chi2();
    float tmp_trk_bendchi2 = iterL1Track->stubPtConsistency();

    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > stubRefs = iterL1Track->getStubRefs();
    int tmp_trk_nstub  = (int) stubRefs.size();

    int tmp_trk_seed = 0;
    // if (SaveTracklet) tmp_trk_seed = (int) iterL1Track->trackSeedType();

    unsigned int tmp_trk_phiSector = iterL1Track->phiSector();

      /*
      int tmp_trk_nPSstub = 0;
      if (SaveTracklet) {
  for (int is=0; is<tmp_trk_nstub; is++) {

    DetId detIdStub = theTrackerGeom->idToDet( (stubRefs.at(is)->clusterRef(0))->getDetId() )->geographicalId();
    DetId stackDetid = tTopo->stack(detIdStub);

    bool isPS = (theTrackerGeom->getDetectorType(stackDetid)==TrackerGeometry::ModuleType::Ph2PSP);
    if (isPS) tmp_trk_nPSstub++;
  }
      }
      */

    // ----------------------------------------------------------------------------------------------
    // loop over stubs on tracks

    //float tmp_trk_bend_chi2 = 0;
    int tmp_trk_dhits=0;
    int tmp_trk_lhits=0;

    // loop over stubs
    for (int is=0; is<tmp_trk_nstub; is++) {

      //detID of stub
      // ---------- For CMSSW_12 --------------
      DetId detIdStub = _trackerGeometry->idToDet( (stubRefs.at(is)->clusterRef(0))->getDetId() )->geographicalId();
      MeasurementPoint coords = stubRefs.at(is)->clusterRef(0)->findAverageLocalCoordinatesCentered();
      const GeomDet* theGeomDet = _trackerGeometry->idToDet(detIdStub);
      Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(coords) );

      double x=posStub.x();
      double y=posStub.y();
      double z=posStub.z();

      int isBarrel = 0;
      int layer=-999999;
      if ( detIdStub.subdetId()==StripSubdetector::TOB ) {
        isBarrel = 1;
        // ---------- For CMSSW_12 ----------------
        layer  = static_cast<int>(_trackerTopology->layer(detIdStub));
        if (DebugMode) cout << "   stub in layer " << layer << " at position x y z = " << x << " " << y << " " << z << endl;
        tmp_trk_lhits+=pow(10,layer-1);
      }
      else if ( detIdStub.subdetId()==StripSubdetector::TID ) {
        isBarrel = 0;
        // ---------- For CMSSW_12 ----------------
        layer  = static_cast<int>(_trackerTopology->layer(detIdStub));
        if (DebugMode) cout << "   stub in disk " << layer << " at position x y z = " << x << " " << y << " " << z << endl;
        tmp_trk_dhits+=pow(10,layer-1);
      }

      m_stub_x.push_back(x);
      m_stub_y.push_back(y);
      m_stub_z.push_back(z);
      m_stub_isBarrel.push_back(isBarrel);
      m_stub_layer.push_back(layer);
    }//end loop over stubs

    int tmp_trk_genuine = 0;
    int tmp_trk_loose = 0;
    int tmp_trk_unknown = 0;
    int tmp_trk_combinatoric = 0;
    // if (MCTruthTTTrackHandle->isLooselyGenuine(l1track_ptr)) tmp_trk_loose = 1;
    // if (MCTruthTTTrackHandle->isGenuine(l1track_ptr)) tmp_trk_genuine = 1;
    // if (MCTruthTTTrackHandle->isUnknown(l1track_ptr)) tmp_trk_unknown = 1;
    // if (MCTruthTTTrackHandle->isCombinatoric(l1track_ptr)) tmp_trk_combinatoric = 1;

    if (DebugMode) {
    cout << "L1 track, pt: " << tmp_trk_pt << " eta: " << tmp_trk_eta << " phi: " << tmp_trk_phi
       << " z0: " << tmp_trk_z0 << " chi2: " << tmp_trk_chi2 << " nstub: " << tmp_trk_nstub;
    if (tmp_trk_genuine) cout << " (is genuine)" << endl;
    if (tmp_trk_unknown) cout << " (is unknown)" << endl;
    if (tmp_trk_combinatoric) cout << " (is combinatoric)" << endl;
    }

  m_trk_pt.push_back(tmp_trk_pt);
  m_trk_eta.push_back(tmp_trk_eta);
  m_trk_phi.push_back(tmp_trk_phi);
  m_trk_z0.push_back(tmp_trk_z0);
  m_trk_d0.push_back(tmp_trk_d0);
  m_trk_rInv.push_back(tmp_trk_rInv);
  m_trk_tanL.push_back(tmp_trk_tanL);
  m_trk_MVA1.push_back(tmp_trk_MVA1);
  m_trk_MVA2.push_back(tmp_trk_MVA2);
  m_trk_MVA3.push_back(tmp_trk_MVA3);
  m_trk_chi2.push_back(tmp_trk_chi2);
  m_trk_bendchi2.push_back(tmp_trk_bendchi2);
  m_trk_nstub.push_back(tmp_trk_nstub);
  // m_trk_dhits.push_back(tmp_trk_dhits);
  // m_trk_lhits.push_back(tmp_trk_lhits);
  m_trk_seed.push_back(tmp_trk_seed);
  m_trk_phiSector.push_back(tmp_trk_phiSector);
  m_trk_genuine.push_back(tmp_trk_genuine);
  m_trk_loose.push_back(tmp_trk_loose);
  m_trk_unknown.push_back(tmp_trk_unknown);
  m_trk_combinatoric.push_back(tmp_trk_combinatoric);

  MuonHLTobjCorrelator::L1TTTrack theTrack( *iterL1Track );

  mTTTrackMap.insert( make_pair(theTrack, static_cast<unsigned int>(this_l1track)) );

  // ----------------------------------------------------------------------------------------------
  // for studying the fake rate
  // ----------------------------------------------------------------------------------------------
  /*
  edm::Ptr< TrackingParticle > my_tp = MCTruthTTTrackHandle->findTrackingParticlePtr(l1track_ptr);

  int myFake = 0;

  int myTP_pdgid = -999;
  float myTP_pt = -999;
  float myTP_eta = -999;
  float myTP_phi = -999;
  float myTP_z0 = -999;
  float myTP_dxy = -999;

  if (my_tp.isNull()) myFake = 0;
  else {
  int tmp_eventid = my_tp->eventId().event();

  if (tmp_eventid > 0) myFake = 2;
  else myFake = 1;

  myTP_pdgid = my_tp->pdgId();
  myTP_pt = my_tp->p4().pt();
  myTP_eta = my_tp->p4().eta();
  myTP_phi = my_tp->p4().phi();
  myTP_z0 = my_tp->vertex().z();

  float myTP_x0 = my_tp->vertex().x();
  float myTP_y0 = my_tp->vertex().y();
  myTP_dxy = sqrt(myTP_x0*myTP_x0 + myTP_y0*myTP_y0);

  if (DebugMode) {
    cout << "TP matched to track has pt = " << my_tp->p4().pt() << " eta = " << my_tp->momentum().eta()
         << " phi = " << my_tp->momentum().phi() << " z0 = " << my_tp->vertex().z()
         << " pdgid = " <<  my_tp->pdgId() << " dxy = " << myTP_dxy << endl;
        }
      }

  m_trk_fake->push_back(myFake);

  m_trk_matchtp_pdgid->push_back(myTP_pdgid);
  m_trk_matchtp_pt->push_back(myTP_pt);
  m_trk_matchtp_eta->push_back(myTP_eta);
  m_trk_matchtp_phi->push_back(myTP_phi);
  m_trk_matchtp_z0->push_back(myTP_z0);
  m_trk_matchtp_dxy->push_back(myTP_dxy);
  */

  this_l1track++;
    }//l1 track loop
  }//save all tracks

  // Now(May2023), TkMuon.h (in HLT-TDR) --> TrackerMuon.h (See comments in https://github.com/cms-sw/cmssw/blob/a7e908cdc4a22f35aaf0775d9f84b068d9bd2f7d/HLTrigger/HLTfilters/plugins/L1TTkMuonFilter.cc#L130-L133)
  edm::Handle<l1t::TrackerMuonCollection> TkMuon;
  iEvent.getByToken(TkMuonToken_,TkMuon);
  for(auto Tkmu=TkMuon->begin(); Tkmu!=TkMuon->end(); ++Tkmu)
  {
    // https://github.com/cms-sw/cmssw/blob/09b17fcfb3900782ab78ad6e0c76e1957c94ff71/DataFormats/L1TCorrelator/interface/TkMuon.h
    // https://github.com/cms-sw/cmssw/blob/2ee90040994dde0e6a4fe686de614b09916bd80a/DataFormats/L1TrackTrigger/interface/TTTrack.h
    // cout<< "phi"<<Tkmu->trkPtr()->phi()<<endl;
    // cout<< "eta"<<Tkmu->trkPtr()->eta()<<endl;
    // cout<< "z0"<<Tkmu->trkPtr()->z0()<<endl;

    // l1TkMuon
    mL1TkMu_pt.push_back( Tkmu->pt() );
    mL1TkMu_eta.push_back( Tkmu->eta() );
    mL1TkMu_phi.push_back( Tkmu->phi() );

    mL1TkMu_trkIsol.push_back( Tkmu->hwIso() );
    mL1TkMu_trkzVtx.push_back( Tkmu->hwZ0() );
    //mL1TkMu_dR.push_back( Tkmu->dR() );

    //mL1TkMu_nTracksMatched.push_back( Tkmu->nTracksMatched() );
    //mL1TkMu_trackCurvature.push_back( Tkmu->trackCurvature() );

    mL1TkMu_quality.push_back( Tkmu->hwQual() );
    //mL1TkMu_pattern.push_back( Tkmu->pattern() );
    //mL1TkMu_muonDetector.push_back( Tkmu->muonDetector() );

    // if (Tkmu->muRef().isNonnull()) {
    //if (Tkmu->muonDetector() != 3 && Tkmu->muRef().isNonnull()) {
    /*
    if (Tkmu->muonRef().isNonnull()) {
      auto regionalCandidate = Tkmu->muonRef().get();
      mL1TkMu_muRefHwPt.push_back( static_cast<float>(regionalCandidate->hwPt())*0.5 );
      mL1TkMu_muRefHwDXY.push_back( regionalCandidate->hwDXY() ); // 4 bit information, don't know how to decode
      mL1TkMu_muRefHwEta.push_back( static_cast<float>(regionalCandidate->hwEta())*0.010875 );
      // mL1TkMu_muRefHwPhi.push_back( static_cast<float>(regionalCandidate->hwPhi())*2.*M_PI/576. );
      float muRefHwPhi = static_cast<float>(
        l1t::MicroGMTConfiguration::calcGlobalPhi(
          regionalCandidate->hwPhi(),
          regionalCandidate->trackFinderType(),
          regionalCandidate->processor()
        )
      );
      mL1TkMu_muRefHwPhi.push_back( reco::reduceRange( muRefHwPhi*2.*M_PI/576. ) );
      mL1TkMu_muRefHwSign.push_back( ( regionalCandidate->hwSign()%2==0 ) ? 1 : -1 );
      mL1TkMu_muRefHwSignValid.push_back( regionalCandidate->hwSignValid() );
      mL1TkMu_muRefHwQual.push_back( regionalCandidate->hwQual() );
    }
    else if( Tkmu->emtfTrk().isNonnull() ) {
      mL1TkMu_muRefHwPt.push_back( Tkmu->emtfTrk()->Pt() );
      mL1TkMu_muRefHwDXY.push_back( -99999 );
      mL1TkMu_muRefHwEta.push_back( Tkmu->emtfTrk()->Eta() );
      mL1TkMu_muRefHwPhi.push_back( reco::reduceRange( angle_units::operators::convertDegToRad(Tkmu->emtfTrk()->Phi_glob()) ) );
      mL1TkMu_muRefHwSign.push_back( Tkmu->emtfTrk()->Charge() );
      mL1TkMu_muRefHwSignValid.push_back( -99999 );
      mL1TkMu_muRefHwQual.push_back( -99999 );
    }
    else {
      // this should never happen
      mL1TkMu_muRefHwPt.push_back( -99999. );
      mL1TkMu_muRefHwDXY.push_back( -99999 );
      mL1TkMu_muRefHwEta.push_back( -99999. );
      mL1TkMu_muRefHwPhi.push_back( -99999. );
      mL1TkMu_muRefHwSign.push_back( -99999 );
      mL1TkMu_muRefHwSignValid.push_back( -99999 );
      mL1TkMu_muRefHwQual.push_back( -99999 );
    }
    */

    auto theTTTrack = Tkmu->trkPtr();
    MuonHLTobjCorrelator::L1TTTrack theTTTrackTmp( *theTTTrack );

    std::map<MuonHLTobjCorrelator::L1TTTrack,unsigned int>::const_iterator where = mTTTrackMap.find(theTTTrackTmp);
    int linkNo = (where==mTTTrackMap.end()) ? -1 : static_cast<int>(mTTTrackMap[theTTTrackTmp]);
    mL1TkMu_TTTpointer.push_back( linkNo );
  }
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

void MuonHLTNtupler::Fill_HLTMuon(const edm::Event &iEvent)
{
  ///////////////////
  // -- L3 Muon -- //
  ///////////////////
  edm::Handle<reco::RecoChargedCandidateCollection> h_L3Muon;
  if( iEvent.getByToken( t_L3Muon_, h_L3Muon ) )
  {
    int _nL3Muon = 0;
    for(unsigned int i_L3=0; i_L3<h_L3Muon->size(); i_L3++)
    {
      reco::RecoChargedCandidateRef ref_L3Mu(h_L3Muon, _nL3Muon);

      L3Muon_pt_[_nL3Muon]     = ref_L3Mu->pt();
      L3Muon_eta_[_nL3Muon]    = ref_L3Mu->eta();
      L3Muon_phi_[_nL3Muon]    = ref_L3Mu->phi();
      L3Muon_charge_[_nL3Muon] = ref_L3Mu->charge();

      reco::TrackRef trackRef = ref_L3Mu->track();
      L3Muon_trkPt_[_nL3Muon] = trackRef->pt();

      _nL3Muon++;
    }
    nL3Muon_ = _nL3Muon;
  } // -- if( L3 handle is valid ) -- //


  ///////////////////
  // -- L2 Muon -- //
  ///////////////////
  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
  if( iEvent.getByToken( t_L2Muon_, h_L2Muon ) )
  {
    int _nL2Muon = 0;
    for( unsigned int i_L2=0; i_L2<h_L2Muon->size(); i_L2++)
    {
      reco::RecoChargedCandidateRef ref_L2Mu(h_L2Muon, _nL2Muon);

      L2Muon_pt_[_nL2Muon]     = ref_L2Mu->pt();
      L2Muon_eta_[_nL2Muon]    = ref_L2Mu->eta();
      L2Muon_phi_[_nL2Muon]    = ref_L2Mu->phi();
      L2Muon_charge_[_nL2Muon] = ref_L2Mu->charge();

      reco::TrackRef trackRef = ref_L2Mu->track();
      L2Muon_trkPt_[_nL2Muon] = trackRef->pt();

      _nL2Muon++;
    }
    nL2Muon_ = _nL2Muon;
  }

  ///////////////////
  // -- Tk Muon -- //
  ///////////////////
  edm::Handle<reco::RecoChargedCandidateCollection> h_TkMuon;
  if( iEvent.getByToken( t_TkMuon_, h_TkMuon ) )
  {
    int _nTkMuon = 0;
    for( unsigned int i_Tk=0; i_Tk<h_TkMuon->size(); i_Tk++)
    {
      reco::RecoChargedCandidateRef ref_TkMu(h_TkMuon, _nTkMuon);

      TkMuon_pt_[_nTkMuon]     = ref_TkMu->pt();
      TkMuon_eta_[_nTkMuon]    = ref_TkMu->eta();
      TkMuon_phi_[_nTkMuon]    = ref_TkMu->phi();
      TkMuon_charge_[_nTkMuon] = ref_TkMu->charge();

      reco::TrackRef trackRef = ref_TkMu->track();
      TkMuon_trkPt_[_nTkMuon] = trackRef->pt();

      _nTkMuon++;
    }
    nTkMuon_ = _nTkMuon;
  }
}

void MuonHLTNtupler::Fill_L1Muon(const edm::Event &iEvent)
{
  edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  if( iEvent.getByToken(t_L1Muon_, h_L1Muon) )
  {
    int _nL1Muon = 0;
    for(int ibx = h_L1Muon->getFirstBX(); ibx<=h_L1Muon->getLastBX(); ++ibx)
    {
      if(ibx != 0) continue; // -- only take when ibx == 0 -- //
      for(auto it=h_L1Muon->begin(ibx); it!=h_L1Muon->end(ibx); it++)
      {
        l1t::MuonRef ref_L1Mu(h_L1Muon, distance(h_L1Muon->begin(h_L1Muon->getFirstBX()), it) );

        L1Muon_pt_[_nL1Muon]      = ref_L1Mu->pt();
        L1Muon_eta_[_nL1Muon]     = ref_L1Mu->eta();
        L1Muon_phi_[_nL1Muon]     = ref_L1Mu->phi();
        L1Muon_charge_[_nL1Muon]  = ref_L1Mu->charge();
        L1Muon_quality_[_nL1Muon] = ref_L1Mu->hwQual();

        L1Muon_etaAtVtx_[_nL1Muon] = ref_L1Mu->etaAtVtx();
        L1Muon_phiAtVtx_[_nL1Muon] = ref_L1Mu->phiAtVtx();

        _nL1Muon++;
      }
    }
    nL1Muon_ = _nL1Muon;
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

    //if( abs(parCand.pdgId()) == 13 || parCand.isHardProcess() ) // -- only muons -- //
    //{
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
  nGenParticle_ = _nGenParticle;
}

void MuonHLTNtupler::Fill_IterL3(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //////////////////////////
  // -- IterL3Muon -- //
  //////////////////////////
  edm::Handle< std::vector<reco::Muon> > h_iterL3Muon;
  edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
  iEvent.getByToken(associatorToken, theAssociator);
  edm::Handle<TrackingParticleCollection> TPCollection;
  iEvent.getByToken(trackingParticleToken, TPCollection);

  // edm::ESHandle<MagneticField> magfieldH;
  // iSetup.get<IdealMagneticFieldRecord>().get(magfieldH);
  edm::ESHandle<MagneticField> magfieldH = iSetup.getHandle(magFieldESToken_);

  // edm::ESHandle<GeometricDet> geomDet;
  // iSetup.get<IdealGeometryRecord>().get(geomDet);
  edm::ESHandle<GeometricDet> geomDet = iSetup.getHandle(geomDetESToken_);

  // edm::ESHandle<TrackerTopology> trkTopo;
  // iSetup.get<TrackerTopologyRcd>().get(trkTopo);
  edm::ESHandle<TrackerTopology> trkTopo = iSetup.getHandle(trackerTopologyESToken_);

  // edm::ESHandle<TrackerGeometry> trkGeom;
  // iSetup.get<TrackerDigiGeometryRecord>().get(trkGeom);
  edm::ESHandle<TrackerGeometry> trkGeom = iSetup.getHandle(trackerGeometryESToken_);

  GeometricSearchTrackerBuilder builder;
  GeometricSearchTracker geomTracker = *(builder.build(&(*geomDet), &(*trkGeom), &(*trkTopo)));

  if( iEvent.getByToken( t_iterL3Muon_, h_iterL3Muon) )
  {
    int _nIterL3Muon = 0;
    for( auto i=0U; i<h_iterL3Muon->size(); ++i )
    {
      const auto& muon(h_iterL3Muon->at(i));

      iterL3Muon_pt_[_nIterL3Muon]     = muon.pt();
      iterL3Muon_eta_[_nIterL3Muon]    = muon.eta();
      iterL3Muon_phi_[_nIterL3Muon]    = muon.phi();
      iterL3Muon_charge_[_nIterL3Muon] = muon.charge();

      if( muon.isGlobalMuon() )     iterL3Muon_isGLB_[_nIterL3Muon] = 1;
      if( muon.isStandAloneMuon() ) iterL3Muon_isSTA_[_nIterL3Muon] = 1;
      if( muon.isTrackerMuon() )    iterL3Muon_isTRK_[_nIterL3Muon] = 1;

      reco::TrackRef innerTrk = muon.innerTrack();
      if( innerTrk.isNonnull() ) {
        iterL3Muon_innerPt_[_nIterL3Muon] = innerTrk->pt();

        tmpTrk trkTmp(innerTrk);
        iterL3IDpassed.push_back(trkTmp);

        const PTrajectoryStateOnDet tmpseed = innerTrk->seedRef()->startingState();
        tmpTSOD tsod(tmpseed);
        MuonIterSeedMap.insert(make_pair(tsod,i));
      }
      else {
        cout << "IterL3Muon: innerTrk.isNonnull(): this should never happen" << endl;
        tmpTrk trkTmp( -99999. );  // dummy tmpTrk
        iterL3IDpassed.push_back(trkTmp);
      }

      _nIterL3Muon++;
    } // -- end of muon iteration

    nIterL3Muon_ = _nIterL3Muon;
  } // -- if getByToken is valid

  ////////////////////
  // -- IterL3OI -- //
  ////////////////////
  edm::Handle< std::vector<reco::MuonTrackLinks> > h_iterL3OI;
  if( iEvent.getByToken( t_iterL3OI_, h_iterL3OI ) )
  {
    int _nIterL3OI = 0;
    for( unsigned int i=0; i<h_iterL3OI->size(); i++)
    {
      if( h_iterL3OI->at(i).trackerTrack().isNonnull() )
      {
        iterL3OI_inner_pt_[_nIterL3OI]     = h_iterL3OI->at(i).trackerTrack()->pt();
        iterL3OI_inner_eta_[_nIterL3OI]    = h_iterL3OI->at(i).trackerTrack()->eta();
        iterL3OI_inner_phi_[_nIterL3OI]    = h_iterL3OI->at(i).trackerTrack()->phi();
        iterL3OI_inner_charge_[_nIterL3OI] = h_iterL3OI->at(i).trackerTrack()->charge();
      }
      if( h_iterL3OI->at(i).standAloneTrack().isNonnull() )
      {
        iterL3OI_outer_pt_[_nIterL3OI]     = h_iterL3OI->at(i).standAloneTrack()->pt();
        iterL3OI_outer_eta_[_nIterL3OI]    = h_iterL3OI->at(i).standAloneTrack()->eta();
        iterL3OI_outer_phi_[_nIterL3OI]    = h_iterL3OI->at(i).standAloneTrack()->phi();
        iterL3OI_outer_charge_[_nIterL3OI] = h_iterL3OI->at(i).standAloneTrack()->charge();
      }
      if( h_iterL3OI->at(i).globalTrack().isNonnull() )
      {
        iterL3OI_global_pt_[_nIterL3OI]     = h_iterL3OI->at(i).globalTrack()->pt();
        iterL3OI_global_eta_[_nIterL3OI]    = h_iterL3OI->at(i).globalTrack()->eta();
        iterL3OI_global_phi_[_nIterL3OI]    = h_iterL3OI->at(i).globalTrack()->phi();
        iterL3OI_global_charge_[_nIterL3OI] = h_iterL3OI->at(i).globalTrack()->charge();
      }
      _nIterL3OI++;
    }
    nIterL3OI_ = _nIterL3OI;
  }

  //////////////////////////
  // -- IterL3IOFromL2 -- //
  //////////////////////////
  edm::Handle< std::vector<reco::MuonTrackLinks> > h_iterL3IOFromL2;
  if( iEvent.getByToken( t_iterL3IOFromL2_, h_iterL3IOFromL2 ) )
  {
    int _nIterL3IOFromL2 = 0;
    for( unsigned int i=0; i<h_iterL3IOFromL2->size(); i++)
    {
      if( h_iterL3IOFromL2->at(i).trackerTrack().isNonnull() )
      {
        iterL3IOFromL2_inner_pt_[_nIterL3IOFromL2]     = h_iterL3IOFromL2->at(i).trackerTrack()->pt();
        iterL3IOFromL2_inner_eta_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).trackerTrack()->eta();
        iterL3IOFromL2_inner_phi_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).trackerTrack()->phi();
        iterL3IOFromL2_inner_charge_[_nIterL3IOFromL2] = h_iterL3IOFromL2->at(i).trackerTrack()->charge();
      }
      if( h_iterL3IOFromL2->at(i).standAloneTrack().isNonnull() )
      {
        iterL3IOFromL2_outer_pt_[_nIterL3IOFromL2]     = h_iterL3IOFromL2->at(i).standAloneTrack()->pt();
        iterL3IOFromL2_outer_eta_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).standAloneTrack()->eta();
        iterL3IOFromL2_outer_phi_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).standAloneTrack()->phi();
        iterL3IOFromL2_outer_charge_[_nIterL3IOFromL2] = h_iterL3IOFromL2->at(i).standAloneTrack()->charge();
      }
      if( h_iterL3IOFromL2->at(i).globalTrack().isNonnull() )
      {
        iterL3IOFromL2_global_pt_[_nIterL3IOFromL2]     = h_iterL3IOFromL2->at(i).globalTrack()->pt();
        iterL3IOFromL2_global_eta_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).globalTrack()->eta();
        iterL3IOFromL2_global_phi_[_nIterL3IOFromL2]    = h_iterL3IOFromL2->at(i).globalTrack()->phi();
        iterL3IOFromL2_global_charge_[_nIterL3IOFromL2] = h_iterL3IOFromL2->at(i).globalTrack()->charge();
      }
      _nIterL3IOFromL2++;
    }
    nIterL3IOFromL2_ = _nIterL3IOFromL2;
  }

  ////////////////////////////////
  // -- IterL3FromL2 (OI+IO) -- //
  ////////////////////////////////
  edm::Handle< std::vector<reco::MuonTrackLinks> > h_iterL3FromL2;
  if( iEvent.getByToken( t_iterL3FromL2_, h_iterL3FromL2 ) )
  {
    int _nIterL3FromL2 = 0;
    for( unsigned int i=0; i<h_iterL3FromL2->size(); i++)
    {
      if( h_iterL3FromL2->at(i).trackerTrack().isNonnull() )
      {
        iterL3FromL2_inner_pt_[_nIterL3FromL2]     = h_iterL3FromL2->at(i).trackerTrack()->pt();
        iterL3FromL2_inner_eta_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).trackerTrack()->eta();
        iterL3FromL2_inner_phi_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).trackerTrack()->phi();
        iterL3FromL2_inner_charge_[_nIterL3FromL2] = h_iterL3FromL2->at(i).trackerTrack()->charge();
      }
      if( h_iterL3FromL2->at(i).standAloneTrack().isNonnull() )
      {
        iterL3FromL2_outer_pt_[_nIterL3FromL2]     = h_iterL3FromL2->at(i).standAloneTrack()->pt();
        iterL3FromL2_outer_eta_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).standAloneTrack()->eta();
        iterL3FromL2_outer_phi_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).standAloneTrack()->phi();
        iterL3FromL2_outer_charge_[_nIterL3FromL2] = h_iterL3FromL2->at(i).standAloneTrack()->charge();
      }
      if( h_iterL3FromL2->at(i).globalTrack().isNonnull() )
      {
        iterL3FromL2_global_pt_[_nIterL3FromL2]     = h_iterL3FromL2->at(i).globalTrack()->pt();
        iterL3FromL2_global_eta_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).globalTrack()->eta();
        iterL3FromL2_global_phi_[_nIterL3FromL2]    = h_iterL3FromL2->at(i).globalTrack()->phi();
        iterL3FromL2_global_charge_[_nIterL3FromL2] = h_iterL3FromL2->at(i).globalTrack()->charge();
      }
      _nIterL3FromL2++;
    }
    nIterL3FromL2_ = _nIterL3FromL2;
  }

  //////////////////////////
  // -- IterL3IOFromL1 -- //
  //////////////////////////
  edm::Handle< std::vector<reco::Track> > h_iterL3IOFromL1;
  if( iEvent.getByToken( t_iterL3IOFromL1_, h_iterL3IOFromL1 ) )
  {
    int _nIterL3IOFromL1 = 0;
    for( unsigned int i=0; i<h_iterL3IOFromL1->size(); i++)
    {
      iterL3IOFromL1_pt_[_nIterL3IOFromL1]     = h_iterL3IOFromL1->at(i).pt();
      iterL3IOFromL1_eta_[_nIterL3IOFromL1]    = h_iterL3IOFromL1->at(i).eta();
      iterL3IOFromL1_phi_[_nIterL3IOFromL1]    = h_iterL3IOFromL1->at(i).phi();
      iterL3IOFromL1_charge_[_nIterL3IOFromL1] = h_iterL3IOFromL1->at(i).charge();
      _nIterL3IOFromL1++;
    }
    nIterL3IOFromL1_ = _nIterL3IOFromL1;
  }

  //////////////////////////
  // -- IterL3MuonNoID -- //
  //////////////////////////
  edm::Handle< std::vector<reco::Muon> > h_iterL3MuonNoID;
  if( iEvent.getByToken( t_iterL3MuonNoID_, h_iterL3MuonNoID) )
  {
    int _nIterL3MuonNoID = 0;
    for( auto i=0U; i<h_iterL3MuonNoID->size(); ++i )
    {
      const auto& muon(h_iterL3MuonNoID->at(i));

      iterL3MuonNoID_pt_[_nIterL3MuonNoID]     = muon.pt();
      iterL3MuonNoID_eta_[_nIterL3MuonNoID]    = muon.eta();
      iterL3MuonNoID_phi_[_nIterL3MuonNoID]    = muon.phi();
      iterL3MuonNoID_charge_[_nIterL3MuonNoID] = muon.charge();

      if( muon.isGlobalMuon() )     iterL3MuonNoID_isGLB_[_nIterL3MuonNoID] = 1;
      if( muon.isStandAloneMuon() ) iterL3MuonNoID_isSTA_[_nIterL3MuonNoID] = 1;
      if( muon.isTrackerMuon() )    iterL3MuonNoID_isTRK_[_nIterL3MuonNoID] = 1;

      reco::TrackRef innerTrk = muon.innerTrack();
      if( innerTrk.isNonnull() ) {
        iterL3MuonNoID_innerPt_[_nIterL3MuonNoID] = innerTrk->pt();

        tmpTrk trkTmp(innerTrk);
        iterL3NoIDpassed.push_back(trkTmp);

        // const PTrajectoryStateOnDet tmpseed = innerTrk->seedRef()->startingState();
        // tmpTSOD tsod(tmpseed);
        // MuonIterNoIdSeedMap.insert(make_pair(tsod,i));
      }
      else {
        // cout << "IterL3MuonNoID: innerTrk.isNonnull(): this should never happen" << endl;
        tmpTrk trkTmp( -99999. );  // dummy tmpTrk
        iterL3NoIDpassed.push_back(trkTmp);
      }

      _nIterL3MuonNoID++;
    } // -- end of muon iteration

    nIterL3MuonNoID_ = _nIterL3MuonNoID;
  } // -- if getByToken is valid

  //////////////////////////
  // -- Tracks from each algo -- //
  //////////////////////////
  // edm::ESHandle<TrackerGeometry> tracker;
  // iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
  edm::ESHandle<TrackerGeometry> tracker = iSetup.getHandle(trackerGeometryESToken_);

  if(doMVA) {
    fill_trackTemplateMva(iEvent, t_hltIter2IterL3FromL1MuonTrack_, theAssociator, TPCollection, tracker, mvaPhase2HltIter2IterL3FromL1MuonPixelSeeds_, hltIter2IterL3FromL1MuonTrackMap,TThltIter2IterL3FromL1MuonTrack, magfieldH, iSetup, geomTracker);
  }
}

void MuonHLTNtupler::Fill_Seed(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  // TrackerHitAssociator associate(iEvent, trackerHitAssociatorConfig_);
  edm::ESHandle<TrackerGeometry> tracker = iSetup.getHandle(trackerGeometryESToken_);

  //////////////////////////
  // -- hltIterL3OISeedsFromL2Muons -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIterL3OISeedsFromL2Muons;
  if( iEvent.getByToken( t_hltIterL3OISeedsFromL2Muons_, h_hltIterL3OISeedsFromL2Muons) )
  {
    for( auto i=0U; i<h_hltIterL3OISeedsFromL2Muons->size(); ++i )
    {
      const auto& seed(h_hltIterL3OISeedsFromL2Muons->at(i));

      SThltIterL3OISeedsFromL2Muons->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIterL3OISeedsFromL2Muons->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIterL3OIMuonTrackMap.find(seedTsod);
      SThltIterL3OISeedsFromL2Muons->linktmpL3((where2==hltIterL3OIMuonTrackMap.end()) ? -1 : hltIterL3OIMuonTrackMap[seedTsod]);
      SThltIterL3OISeedsFromL2Muons->matchingL3((where2==hltIterL3OIMuonTrackMap.end()) ? -1 : TThltIterL3OIMuonTrack->matchedIDpassedL3(hltIterL3OIMuonTrackMap[seedTsod]));

      // std::cout << "OI RecHit loop start" << std::endl;
      // for (auto rechit = seed.recHits().first; rechit != seed.recHits().second; ++rechit) {
      //   // std::cout << "OI RecHit is valid : " << rechit->isValid() << std::endl;
      //   auto matched = associate.associateHit(*rechit);
      //   // std::cout << "OI Matched RecHit size = " << matched.size() << std::endl;
      // }
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter0IterL3MuonPixelSeedsFromPixelTracks -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter0IterL3MuonPixelSeedsFromPixelTracks;
  if( iEvent.getByToken( t_hltIter0IterL3MuonPixelSeedsFromPixelTracks_, h_hltIter0IterL3MuonPixelSeedsFromPixelTracks) )
  {
    for( auto i=0U; i<h_hltIter0IterL3MuonPixelSeedsFromPixelTracks->size(); ++i )
    {
      const auto& seed(h_hltIter0IterL3MuonPixelSeedsFromPixelTracks->at(i));

      SThltIter0IterL3MuonPixelSeedsFromPixelTracks->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter0IterL3MuonPixelSeedsFromPixelTracks->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter0IterL3MuonTrackMap.find(seedTsod);
      SThltIter0IterL3MuonPixelSeedsFromPixelTracks->linktmpL3((where2==hltIter0IterL3MuonTrackMap.end()) ? -1 : hltIter0IterL3MuonTrackMap[seedTsod]);
      SThltIter0IterL3MuonPixelSeedsFromPixelTracks->matchingL3((where2==hltIter0IterL3MuonTrackMap.end()) ? -1 : TThltIter0IterL3MuonTrack->matchedIDpassedL3(hltIter0IterL3MuonTrackMap[seedTsod]));

      // std::cout << "IO RecHit loop start" << std::endl;
      // for (auto rechit = seed.recHits().first; rechit != seed.recHits().second; ++rechit) {
      //   // std::cout << "IO RecHit is valid : " << rechit->isValid() << std::endl;
      //   auto matched = associate.associateHit(*rechit);
      //   // std::cout << "IO Matched RecHit size = " << matched.size() << std::endl;
      // }
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter2IterL3MuonPixelSeeds -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter2IterL3MuonPixelSeeds;
  if( iEvent.getByToken( t_hltIter2IterL3MuonPixelSeeds_, h_hltIter2IterL3MuonPixelSeeds) )
  {
    for( auto i=0U; i<h_hltIter2IterL3MuonPixelSeeds->size(); ++i )
    {
      const auto& seed(h_hltIter2IterL3MuonPixelSeeds->at(i));

      SThltIter2IterL3MuonPixelSeeds->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter2IterL3MuonPixelSeeds->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter2IterL3MuonTrackMap.find(seedTsod);
      SThltIter2IterL3MuonPixelSeeds->linktmpL3((where2==hltIter2IterL3MuonTrackMap.end()) ? -1 : hltIter2IterL3MuonTrackMap[seedTsod]);
      SThltIter2IterL3MuonPixelSeeds->matchingL3((where2==hltIter2IterL3MuonTrackMap.end()) ? -1 : TThltIter2IterL3MuonTrack->matchedIDpassedL3(hltIter2IterL3MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter3IterL3MuonPixelSeeds -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter3IterL3MuonPixelSeeds;
  if( iEvent.getByToken( t_hltIter3IterL3MuonPixelSeeds_, h_hltIter3IterL3MuonPixelSeeds) )
  {
    for( auto i=0U; i<h_hltIter3IterL3MuonPixelSeeds->size(); ++i )
    {
      const auto& seed(h_hltIter3IterL3MuonPixelSeeds->at(i));

      SThltIter3IterL3MuonPixelSeeds->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter3IterL3MuonPixelSeeds->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter3IterL3MuonTrackMap.find(seedTsod);
      SThltIter3IterL3MuonPixelSeeds->linktmpL3((where2==hltIter3IterL3MuonTrackMap.end()) ? -1 : hltIter3IterL3MuonTrackMap[seedTsod]);
      SThltIter3IterL3MuonPixelSeeds->matchingL3((where2==hltIter3IterL3MuonTrackMap.end()) ? -1 : TThltIter3IterL3MuonTrack->matchedIDpassedL3(hltIter3IterL3MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks;
  if( iEvent.getByToken( t_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks_, h_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks) )
  {
    for( auto i=0U; i<h_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->size(); ++i )
    {
      const auto& seed(h_hltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->at(i));

      SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter0IterL3FromL1MuonTrackMap.find(seedTsod);
      SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->linktmpL3((where2==hltIter0IterL3FromL1MuonTrackMap.end()) ? -1 : hltIter0IterL3FromL1MuonTrackMap[seedTsod]);
      SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks->matchingL3((where2==hltIter0IterL3FromL1MuonTrackMap.end()) ? -1 : TThltIter0IterL3FromL1MuonTrack->matchedIDpassedL3(hltIter0IterL3FromL1MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter2IterL3FromL1MuonPixelSeeds -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter2IterL3FromL1MuonPixelSeeds;
  if( iEvent.getByToken( t_hltIter2IterL3FromL1MuonPixelSeeds_, h_hltIter2IterL3FromL1MuonPixelSeeds) )
  {
    for( auto i=0U; i<h_hltIter2IterL3FromL1MuonPixelSeeds->size(); ++i )
    {
      const auto& seed(h_hltIter2IterL3FromL1MuonPixelSeeds->at(i));

      SThltIter2IterL3FromL1MuonPixelSeeds->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter2IterL3FromL1MuonPixelSeeds->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter2IterL3FromL1MuonTrackMap.find(seedTsod);
      SThltIter2IterL3FromL1MuonPixelSeeds->linktmpL3((where2==hltIter2IterL3FromL1MuonTrackMap.end()) ? -1 : hltIter2IterL3FromL1MuonTrackMap[seedTsod]);
      SThltIter2IterL3FromL1MuonPixelSeeds->matchingL3((where2==hltIter2IterL3FromL1MuonTrackMap.end()) ? -1 : TThltIter2IterL3FromL1MuonTrack->matchedIDpassedL3(hltIter2IterL3FromL1MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid

  //////////////////////////
  // -- hltIter3IterL3FromL1MuonPixelSeeds -- //
  //////////////////////////
  edm::Handle< TrajectorySeedCollection > h_hltIter3IterL3FromL1MuonPixelSeeds;
  if( iEvent.getByToken( t_hltIter3IterL3FromL1MuonPixelSeeds_, h_hltIter3IterL3FromL1MuonPixelSeeds) )
  {
    for( auto i=0U; i<h_hltIter3IterL3FromL1MuonPixelSeeds->size(); ++i )
    {
      const auto& seed(h_hltIter3IterL3FromL1MuonPixelSeeds->at(i));

      SThltIter3IterL3FromL1MuonPixelSeeds->fill(seed, tracker);

      tmpTSOD seedTsod(seed.startingState());
      std::map<tmpTSOD,unsigned int>::const_iterator where = MuonIterSeedMap.find(seedTsod);
      SThltIter3IterL3FromL1MuonPixelSeeds->linkIterL3((where==MuonIterSeedMap.end()) ? -1 : MuonIterSeedMap[seedTsod]);

      std::map<tmpTSOD,unsigned int>::const_iterator where2 = hltIter3IterL3FromL1MuonTrackMap.find(seedTsod);
      SThltIter3IterL3FromL1MuonPixelSeeds->linktmpL3((where2==hltIter3IterL3FromL1MuonTrackMap.end()) ? -1 : hltIter3IterL3FromL1MuonTrackMap[seedTsod]);
      SThltIter3IterL3FromL1MuonPixelSeeds->matchingL3((where2==hltIter3IterL3FromL1MuonTrackMap.end()) ? -1 : TThltIter3IterL3FromL1MuonTrack->matchedIDpassedL3(hltIter3IterL3FromL1MuonTrackMap[seedTsod]));
    } // -- end of seed iteration
  } // -- if getByToken is valid
}

void MuonHLTNtupler::Fill_TP( const edm::Event &iEvent, tpTemplate* tpTmp )
{
  edm::Handle<TrackingParticleCollection> TPCollection;
  if( iEvent.getByToken(trackingParticleToken, TPCollection) ) {
    for( auto i=0U; i< TPCollection->size(); ++i) {
      if( abs(TPCollection->at(i).pdgId()) == 13 ) {
        tpTmp->fill( TPCollection->at(i) );
      }
    }
  }
}


void MuonHLTNtupler::Fill_Muon2(const edm::Event &iEvent)
{

  using namespace trigger;
  
  edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;
  if ( iEvent.getByToken(theMuonFilteredCollectionToken_, PrevFilterOutput) ){

    std::vector<reco::RecoChargedCandidateRef> muons;
    PrevFilterOutput->getObjects(trigger::TriggerMuon, muons);

    Handle<RecoChargedCandidateCollection> mucands;
    iEvent.getByToken(theMuonCollectionToken_, mucands);
    
    //const auto& muonAssoc = iEvent.get(muonAssocToken_);

    edm::Handle<edm::ValueMap<float>> btlMatchChi2Handle;
    edm::ValueMap<float> btlMatchChi2;
    edm::ValueMap<float> etlMatchChi2;
    edm::ValueMap<float> btlMatchTimeChi2;
    edm::ValueMap<float> etlMatchTimeChi2;
    edm::ValueMap<int> npixBarrel;
    edm::ValueMap<int> npixEndcap;
    edm::ValueMap<float> muonTrackOutermostHitPosition;
    edm::ValueMap<float> muonTrackp;
    edm::ValueMap<float> muonTrackBeta;
    edm::ValueMap<float> muonTrackt0;
    edm::ValueMap<float> muonTracksigmat0;
    edm::ValueMap<float> muonTrackPathLength;
    edm::ValueMap<float> muonTracktmtd;
    edm::ValueMap<float> muonTracksigmatmtd;
    edm::ValueMap<GlobalPoint> muonTrackmtdpos;
    edm::ValueMap<float> muonTrackTofMu;
    edm::ValueMap<float> muonTrackSigmaTofMu;

    bool skipMTD = false;
    if (iEvent.getByToken(muonbtlMatchChi2Token_, btlMatchChi2Handle)){
      btlMatchChi2 = *btlMatchChi2Handle;
      //btlMatchChi2 = iEvent.get(muonbtlMatchChi2Token_);
      etlMatchChi2 = iEvent.get(muonetlMatchChi2Token_);
      btlMatchTimeChi2 = iEvent.get(muonbtlMatchTimeChi2Token_);
      etlMatchTimeChi2 = iEvent.get(muonetlMatchTimeChi2Token_);
      npixBarrel = iEvent.get(muonnpixBarrelToken_);
      npixEndcap = iEvent.get(muonnpixEndcapToken_);
      muonTrackOutermostHitPosition = iEvent.get(muonTrackOutermostHitPositionToken_);
      muonTrackp = iEvent.get(muonTrackpToken_);
      muonTrackBeta = iEvent.get(muonTrackBetaToken_);
      muonTrackt0 = iEvent.get(muonTrackt0Token_);
      muonTracksigmat0 = iEvent.get(muonTracksigmat0Token_);
      muonTrackPathLength = iEvent.get(muonTrackPathLengthToken_);
      muonTracktmtd = iEvent.get(muonTracktmtdToken_);
      muonTracksigmatmtd = iEvent.get(muonTracksigmatmtdToken_);
      muonTrackmtdpos = iEvent.get(muonTrackmtdposToken_);
      muonTrackTofMu = iEvent.get(muonTrackTofMuToken_);
      muonTrackSigmaTofMu = iEvent.get(muonTrackSigmaTofMuToken_);
    }else{
      skipMTD = true;
    }
    
    int _nMuon = 0;
    //unsigned int index = 0;
    for (unsigned int iMu = 0; iMu < mucands->size(); iMu++)
      {

	RecoChargedCandidateRef cand(mucands, iMu);
	if (!triggerdByPreviousLevel(cand, muons))
	  continue;

	reco::RecoChargedCandidateRef mu(mucands, iMu);

	bool hasMTDInfo = false;
	//const reco::TrackRef muontrackref(iEvent.getHandle(theMuonCollectionToken_), index);
	reco::TrackRef muontrackref = mu->track();
	//index++;

	if (skipMTD) {
	  std::cout << "Event without MTD information" << std::endl;
	  hasMTDInfo = false;
	}else if (muonTrackPathLength[muontrackref] == -1) {
	  //if (muonAssoc[muontrackref] == -1) {
	  std::cout << "MuonExtenderWithMTD: track not associated" << std::endl;
	  hasMTDInfo = false;
	}else{
	  std::cout << "MuonExtenderWithMTD: found a valid track" << std::endl;
	  hasMTDInfo = true;
	}
	
	muonCand_pt_[_nMuon]  = mu->pt();
	muonCand_eta_[_nMuon] = mu->eta();
	muonCand_phi_[_nMuon] = mu->phi();
	muonCand_px_[_nMuon]  = mu->px();
	muonCand_py_[_nMuon]  = mu->py();
	muonCand_pz_[_nMuon]  = mu->pz();
	muonCand_charge_[_nMuon] = mu->charge();
	
	if( mu->isGlobalMuon() ) muonCand_isGLB_[_nMuon] = 1;
	if( mu->isStandAloneMuon() ) muonCand_isSTA_[_nMuon] = 1;
	if( mu->isTrackerMuon() ) muonCand_isTRK_[_nMuon] = 1;
	if( mu->isMuon() ) muonCand_isMuon_[_nMuon] = 1;
	
	const reco::Track* innerTrk = mu->bestTrack();

        muonCand_vx_[_nMuon] = innerTrk->vx();
        muonCand_vy_[_nMuon] = innerTrk->vy();
	muonCand_vz_[_nMuon] = innerTrk->vz();
	muonCand_dxy_[_nMuon] = mu->bestTrack().dxy();
	//muonCand_dxyError_bs_[_nMuon] = innerTrk->dxyError(*bs);
	muonCand_dz_[_nMuon] = mu->bestTrack().dz();
	muonCand_dzError_[_nMuon] = innerTrk->dzError();
	//if (innerTrk->dxyError(*bs) > 0.) {
	//  muonCand_IPSig_[_nMuon] = abs(innerTrk->dxy(bs->position()) / innerTrk->dxyError(*bs));
	//}
	muonCand_inner_trkChi2_[_nMuon]             = innerTrk->normalizedChi2();
	muonCand_inner_validFraction_[_nMuon]       = innerTrk->validFraction();
	muonCand_inner_trackerLayers_[_nMuon]       = innerTrk->hitPattern().trackerLayersWithMeasurement();
	muonCand_inner_trackerHits_[_nMuon]         = innerTrk->hitPattern().numberOfValidTrackerHits();
	muonCand_inner_lostTrackerHits_[_nMuon]     = innerTrk->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostTrackerHitsIn_[_nMuon]   = innerTrk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
	muonCand_inner_lostTrackerHitsOut_[_nMuon]  = innerTrk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS);
	muonCand_inner_lostPixelHits_[_nMuon]       = innerTrk->hitPattern().numberOfLostPixelHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostPixelBarrelHits_[_nMuon] = innerTrk->hitPattern().numberOfLostPixelBarrelHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostPixelEndcapHits_[_nMuon] = innerTrk->hitPattern().numberOfLostPixelEndcapHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripHits_[_nMuon]       = innerTrk->hitPattern().numberOfLostStripHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripTIBHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTIBHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripTIDHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTIDHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripTOBHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTOBHits(HitPattern::TRACK_HITS);
	muonCand_inner_lostStripTECHits_[_nMuon]    = innerTrk->hitPattern().numberOfLostStripTECHits(HitPattern::TRACK_HITS);
	muonCand_inner_pixelLayers_[_nMuon]         = innerTrk->hitPattern().pixelLayersWithMeasurement();
	muonCand_inner_pixelHits_[_nMuon]           = innerTrk->hitPattern().numberOfValidPixelHits();
	
	muonCand_global_muonHits_[_nMuon]           = innerTrk->hitPattern().numberOfValidMuonHits();
	muonCand_global_trkChi2_[_nMuon]            = innerTrk->normalizedChi2();
	muonCand_global_trackerLayers_[_nMuon]      = innerTrk->hitPattern().trackerLayersWithMeasurement();
	muonCand_global_trackerHits_[_nMuon]        = innerTrk->hitPattern().numberOfValidTrackerHits();

	if (hasMTDInfo){
	  muonCand_btlMatchChi2_[_nMuon] = btlMatchChi2[muontrackref];
	  muonCand_etlMatchChi2_[_nMuon] =  etlMatchChi2[muontrackref];
	  muonCand_btlMatchTimeChi2_[_nMuon] =  btlMatchTimeChi2[muontrackref];
	  muonCand_etlMatchTimeChi2_[_nMuon] =  etlMatchTimeChi2[muontrackref];
	  muonCand_npixBarrel_[_nMuon] =  npixBarrel[muontrackref];
	  muonCand_npixEndcap_[_nMuon] =  npixEndcap[muontrackref];
	  muonCand_outermostHitPosition_[_nMuon] =  muonTrackOutermostHitPosition[muontrackref];
	  muonCand_p_[_nMuon] =  muonTrackp[muontrackref];
	  muonCand_beta_[_nMuon] =  muonTrackBeta[muontrackref];
	  muonCand_t0_[_nMuon] =  muonTrackt0[muontrackref];
	  muonCand_sigmat0_[_nMuon] =  muonTracksigmat0[muontrackref];
	  muonCand_pathLength_[_nMuon] =  muonTrackPathLength[muontrackref];
	  muonCand_tmtd_[_nMuon] =  muonTracktmtd[muontrackref];
	  muonCand_sigmatmtd_[_nMuon] =  muonTracksigmatmtd[muontrackref];
	  muonCand_tofMu_[_nMuon] =  muonTrackTofMu[muontrackref];
	  muonCand_sigmaTofMu_[_nMuon] =  muonTrackSigmaTofMu[muontrackref];
	  muonCand_mtdpos_x_[_nMuon] =  muonTrackmtdpos[muontrackref].x();
	  muonCand_mtdpos_y_[_nMuon] =  muonTrackmtdpos[muontrackref].y();
	  muonCand_mtdpos_z_[_nMuon] =  muonTrackmtdpos[muontrackref].z();
	}
	
	_nMuon++;
      }
    
    nMuonCand_ = _nMuon;
  }
}

bool MuonHLTNtupler::triggerdByPreviousLevel(const reco::RecoChargedCandidateRef& candref,
					     const std::vector<reco::RecoChargedCandidateRef>& vcands) {
  unsigned int i = 0;
  unsigned int i_max = vcands.size();
  for (; i != i_max; ++i) {
    if (candref == vcands[i])
      return true;
  }

  return false;
}


void MuonHLTNtupler::Fill_PFCand(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  
  edm::Handle<trigger::TriggerFilterObjectWithRefs> PrevFilterOutput;
  if ( iEvent.getByToken(theMuonFilteredCollectionToken_, PrevFilterOutput) ){
    
    std::vector<reco::RecoChargedCandidateRef> muons;
    PrevFilterOutput->getObjects(trigger::TriggerMuon, muons);
    
    Handle<RecoChargedCandidateCollection> mucands;
    iEvent.getByToken(theMuonCollectionToken_, mucands);
    
    edm::Handle<reco::PFCandidateCollection> pfCandidateHandle;
    if ( iEvent.getByToken(pfCandidateProducer_, pfCandidateHandle) ){
      
      int _nPFCand = 0;

      edm::Handle<edm::ValueMap<float>> trackBtlMatchChi2Handle;
      edm::ValueMap<float> trackEtlMatchChi2;
      edm::ValueMap<float> trackEtlMatchChi2;
      edm::ValueMap<float> trackEtlMatchTimeChi2;
      edm::ValueMap<float> trackEtlMatchTimeChi2;
      edm::ValueMap<int> trackNpixBarrel;
      edm::ValueMap<int> trackNpixEndcap;
      edm::ValueMap<float> trackOutermostHitPosition;
      edm::ValueMap<float> trackp;
      edm::ValueMap<float> trackBeta;
      edm::ValueMap<float> trackt0;
      edm::ValueMap<float> tracksigmat0;
      edm::ValueMap<float> trackPathLength;
      edm::ValueMap<float> tracktmtd;
      edm::ValueMap<float> tracksigmatmtd;
      edm::ValueMap<GlobalPoint> trackmtdpos;
      edm::ValueMap<float> trackTofPi;
      edm::ValueMap<float> trackSigmaTofPi;
      edm::ValueMap<float> trackTofK;
      edm::ValueMap<float> trackSigmaTofK;
      edm::ValueMap<float> trackTofP;
      edm::ValueMap<float> trackSigmaTofP;

      bool skipMTD = false;
      if (iEvent.getByToken(trackBtlMatchChi2Token_, trackBtlMatchChi2Handle)){
	trackBtlMatchChi2 = *trackBtlMatchChi2Handle;
	trackEtlMatchChi2 = iEvent.get(trackEtlMatchChi2Token_);
	trackBtlMatchTimeChi2 = iEvent.get(trackBtlMatchTimeChi2Token_);
	trackEtlMatchTimeChi2 = iEvent.get(trackEtlMatchTimeChi2Token_);
	trackNpixBarrel = iEvent.get(trackNpixBarrelToken_);
	trackNpixEndcap = iEvent.get(trackNpixEndcapToken_);
	trackOutermostHitPosition = iEvent.get(trackOutermostHitPositionToken_);
	trackp = iEvent.get(trackpToken_);
	trackBeta = iEvent.get(trackBetaToken_);
	trackt0 = iEvent.get(trackt0Token_);
	tracksigmat0 = iEvent.get(tracksigmat0Token_);
	trackPathLength = iEvent.get(trackPathLengthToken_);
	tracktmtd = iEvent.get(tracktmtdToken_);
	tracksigmatmtd = iEvent.get(tracksigmatmtdToken_);
	trackmtdpos = iEvent.get(trackmtdposToken_);
	trackTofPi = iEvent.get(trackTofPiToken_);
	trackSigmaTofPi = iEvent.get(trackSigmaTofPiToken_);
	trackTofK = iEvent.get(trackTofKToken_);
        trackSigmaTofK = iEvent.get(trackSigmaTofKToken_);
	trackTofP = iEvent.get(trackTofPToken_);
        trackSigmaTofP = iEvent.get(trackSigmaTofPToken_);
      }else{
	skipMTD = true;
      }

      edm::Handle<double> rhoHandle;
      iEvent.getByToken(rhoProducer_, rhoHandle);
      pfcand_rho_ = *(rhoHandle.product());
      
      unsigned int nMuons = mucands->size();
      for (unsigned int iMu = 0; iMu < nMuons; iMu++){
	
        reco::RecoChargedCandidateRef muRef(mucands, iMu);
        reco::TrackRef trackRef = muRef->track();
        const Track& muon = *trackRef;
        if (!triggerdByPreviousLevel(muRef, muons))
          continue;

	for (unsigned int iPF = 0; iPF < pfCandidateHandle->size(); iPF++) {
	  reco::PFCandidateRef pc(pfCandidateHandle, iPF);
	  
	  float dr2 = reco::deltaR2(muon.eta(), muon.phi(), pc->eta(), pc->phi());
	  if (dr2 > drMaxPf_ * drMaxPf_)
	    continue;
	  if (pc->pt() < minEnergyPf_)
	    continue;
	  if (fabs(pc->charge())){
	    if (dr2 < drVetoPfCh_)
	      continue;
	  }else{
	    if (dr2 < drVetoPf_)
	      continue;
	  }

	  bool trackHasMTDInfo = false;
	  const reco::TrackRef trackref = pc->trackRef();
	  if (trackref.isNonnull() && !skipMTD){
	    try {
	      if (t0Src[trackref] == -1) {
		trackHasMTDInfo = false;
	      }else{
		trackHasMTDInfo = true;
	      }
	    }catch (...) {
	      trackHasMTDInfo = false;
	    }
	  }

	  pfcand_pt_[_nPFCand] = pc->pt();
	  pfcand_eta_[_nPFCand] = pc->eta();
	  pfcand_phi_[_nPFCand] = pc->phi();
	  pfcand_charge_[_nPFCand] = pc->charge();
	  pfcand_pdgId_[_nPFCand] = pc->pdgId();
	  pfcand_px_[_nPFCand] = pc->px();
	  pfcand_py_[_nPFCand] = pc->py();
	  pfcand_pz_[_nPFCand] = pc->pz();
	  pfcand_vx_[_nPFCand] = pc->vx();
	  pfcand_vy_[_nPFCand] = pc->vy();
	  pfcand_vz_[_nPFCand] = pc->vz();
	  pfcand_time_[_nPFCand] = (double)pc->time();
	  pfcand_timeErr_[_nPFCand] = (double)pc->timeError();
	  pfcand_dxyErr_[_nPFCand] = pc->dxyError();
	  pfcand_dzErr_[_nPFCand] = pc->dzError();
	  pfcand_vChi2NoF_[_nPFCand] = pc->vertexNdof();
	  if (pc->trackRef().isNonnull()){
	    pfcand_dxy_[_nPFCand] = pc->bestTrack()->dxy();
	    pfcand_dz_[_nPFCand] = pc->bestTrack()->dz();
	  }

	  if (trackHasMTDInfo){
	    pfcand_btlMatchChi2_[_nPFCand] = trackBtlMatchChi2[trackref];
	    pfcand_etlMatchChi2_[_nPFCand] =  trackEtlMatchChi2[trackref];
	    pfcand_btlMatchTimeChi2_[_nPFCand] =  trackBtlMatchTimeChi2[trackref];
	    pfcand_etlMatchTimeChi2_[_nPFCand] =  trackEtlMatchTimeChi2[trackref];
	    pfcand_npixBarrel_[_nPFCand] =  trackNpixBarrel[trackref];
	    pfcand_npixEndcap_[_nPFCand] =  trackNpixEndcap[trackref];
	    pfcand_outermostHitPosition_[_nPFCand] =  trackOutermostHitPosition[trackref];
	    pfcand_p_[_nPFCand] =  trackp[trackref];
	    pfcand_beta_[_nPFCand] =  trackBeta[trackref];
	    pfcand_t0_[_nPFCand] =  trackt0[trackref];
	    pfcand_sigmat0_[_nPFCand] =  tracksigmat0[trackref];
	    pfcand_pathLength_[_nPFCand] =  trackPathLength[trackref];
	    pfcand_tmtd_[_nPFCand] =  tracktmtd[trackref];
	    pfcand_sigmatmtd_[_nPFCand] =  tracksigmatmtd[trackref];
	    pfcand_tofPi_[_nPFCand] =  trackTofPi[trackref];
	    pfcand_sigmaTofPi_[_nPFCand] =  trackSigmaTofPi[trackref];
	    pfcand_tofK_[_nPFCand] =  trackTofK[trackref];
	    pfcand_sigmaTofK_[_nPFCand] =  trackSigmaTofK[trackref];
	    pfcand_tofP_[_nPFCand] =  trackTofP[trackref];
	    pfcand_sigmaTofP_[_nPFCand] =  trackSigmaTofP[trackref];
	    pfcand_mtdpos_x_[_nPFCand] =  trackmtdpos[trackref].x();
	    pfcand_mtdpos_y_[_nPFCand] =  trackmtdpos[trackref].y();
	    pfcand_mtdpos_z_[_nPFCand] =  trackmtdpos[trackref].z();
	  }

	  pfcand_muonIdx_[_nPFCand] = iMu;
	  
	  _nPFCand++;	  
	}
      }
      nPFCand_ = _nPFCand;
    }
  }
}

bool MuonHLTNtupler::computedRVeto(RecoChargedCandidateRef candRef, reco::PFClusterRef pfclu, double drMAX, double drVeto2_) {
  float dR2 = deltaR2(candRef->eta(), candRef->phi(), pfclu->eta(), pfclu->phi());
  if (dR2 > (drMAX * drMAX) || dR2 < drVeto2_)
    return false;
  else
    return true;
}


void MuonHLTNtupler::fill_trackTemplate(
  const edm::Event &iEvent,
  edm::EDGetTokenT<edm::View<reco::Track>>& trkToken,
  edm::EDGetTokenT<reco::RecoToSimCollection>& assoToken,
  trkTemplate* TTtrack,
  bool doIso = false
) {

  edm::Handle<edm::View<reco::Track>> trkHandle;
  if( iEvent.getByToken( trkToken, trkHandle ) ) {

    edm::Handle<reco::RecoChargedCandidateCollection> h_L3Muon;
    if(doIso) {
      bool hasL3Muon = iEvent.getByToken( t_L3Muon_, h_L3Muon );
      if(!hasL3Muon) {
        throw cms::Exception("ConfigurationError")
              << "getByToken failed for L3 Muon Candidates, it is needed for Iso maps";
      }
      if( h_L3Muon->size() != trkHandle->size() ) {
        throw cms::Exception("ConfigurationError")
            << "h_L3Muon->size() != trkHandle->size()";
      }
    }

    edm::Handle<reco::RecoToSimCollection> assoHandle;
    if( iEvent.getByToken( assoToken, assoHandle ) ) {
      auto recSimColl = *assoHandle.product();

      for( unsigned int i = 0; i < trkHandle->size(); i++ ) {
        TTtrack->fill(trkHandle->at(i));

        auto track = trkHandle->refAt(i);
        auto TPfound = recSimColl.find(track);
        if (TPfound != recSimColl.end()) {
          const auto& TPmatch = TPfound->val;
          TTtrack->fillBestTP(TPmatch[0].first);
          TTtrack->fillBestTPsharedFrac(TPmatch[0].second);
          TTtrack->fillmatchedTPsize(TPmatch.size());
        } else {  // sync vector size
          TTtrack->fillDummyTP();
          TTtrack->fillBestTPsharedFrac(-99999.);
          TTtrack->fillmatchedTPsize(0);
        }

        // -- fill dummy
        TTtrack->linkIterL3(-1);
        TTtrack->linkIterL3NoId(-1);
        TTtrack->fillMva( -99999. );

        if(doIso) {
          vector<float> trkIsolations = {};
          vector<float> pfIsolations = {};

          reco::RecoChargedCandidateRef muRef(h_L3Muon, i);
          // if( fabs(trkHandle->at(i).pt() - muRef->pt()) / muRef->pt() > 0.001 ||
          //     fabs(trkHandle->at(i).eta() - muRef->eta()) > 0.001 ||
          //     fabs(trkHandle->at(i).phi() - muRef->phi()) > 0.001
          // ) {
          //   throw cms::Exception("ConfigurationError")
          //       << "L3 muon candidate != corresponding track";
          // }

          for( unsigned int ii = 0; ii < trkIsoTags_.size(); ++ii) {

            edm::Handle<reco::IsoDepositMap> trkIsoMap;
            if( iEvent.getByToken(trkIsoTokens_.at(ii), trkIsoMap) ) {

              reco::IsoDeposit trkIso = (*trkIsoMap)[muRef];
              TString Tag_tstr = TString(trkIsoTags_.at(ii));
              float dR = 0.0;
              if(Tag_tstr.Contains("dR0p1"))       dR = 0.1;
              else if(Tag_tstr.Contains("dR0p2"))  dR = 0.2;
              else if(Tag_tstr.Contains("dR0p3"))  dR = 0.3;
              else if(Tag_tstr.Contains("dR0p4"))  dR = 0.4;
              else if(Tag_tstr.Contains("dR0p5"))  dR = 0.5;
              else if(Tag_tstr.Contains("dR0p6"))  dR = 0.6;
              else if(Tag_tstr.Contains("dR0p7"))  dR = 0.7;
              else if(Tag_tstr.Contains("dR0p8"))  dR = 0.8;
              else if(Tag_tstr.Contains("dR0p9"))  dR = 0.9;
              else if(Tag_tstr.Contains("dR1p0"))  dR = 1.0;
              else                                 dR = 0.0;

              trkIsolations.push_back( trkIso.depositWithin(dR) );
            }
            else {
              trkIsolations.push_back( -99999. );
            }
          }

          for( unsigned int ii = 0; ii < pfIsoTags_.size(); ++ii) {

            edm::Handle<reco::RecoChargedCandidateIsolationMap> pfIsoMap;
            if( iEvent.getByToken(pfIsoTokens_.at(ii), pfIsoMap) ) {

              reco::RecoChargedCandidateIsolationMap::const_iterator pfIso = (*pfIsoMap).find( muRef );

              pfIsolations.push_back( pfIso->val );
            }
            else {
              pfIsolations.push_back( -99999. );
            }
          }

          TTtrack->fillIso( trkIsolations, pfIsolations );
        }
      }
    }
  }
}

void MuonHLTNtupler::fill_tpTemplate(
  const edm::Event &iEvent,
  edm::EDGetTokenT<reco::SimToRecoCollection>& assoToken,
  tpTemplate* TTtp
) {

  edm::Handle<TrackingParticleCollection> TPCollection;
  if( iEvent.getByToken(trackingParticleToken, TPCollection) ) {

    edm::Handle<reco::SimToRecoCollection> assoHandle;
    if( iEvent.getByToken( assoToken, assoHandle ) ) {
      auto simRecColl = *assoHandle.product();

      for( unsigned int i = 0; i < TPCollection->size(); i++ ) {

        auto tp = TPCollection->at(i);

        if( !(tp.eventId().bunchCrossing() == 0 && tp.eventId().event() == 0) )
          continue;

        if( abs( tp.pdgId() ) != 13 )
          continue;

        bool isStable = true;
        for (TrackingParticle::genp_iterator j = tp.genParticle_begin(); j != tp.genParticle_end(); ++j) {
          if (j->get() == nullptr || j->get()->status() != 1) {
            isStable = false;
          }
        }
        if( tp.status() == -99 && (std::abs(tp.pdgId()) != 11 && std::abs(tp.pdgId()) != 13 && std::abs(tp.pdgId()) != 211 &&
                                   std::abs(tp.pdgId()) != 321 && std::abs(tp.pdgId()) != 2212 && std::abs(tp.pdgId()) != 3112 &&
                                   std::abs(tp.pdgId()) != 3222 && std::abs(tp.pdgId()) != 3312 && std::abs(tp.pdgId()) != 3334)
        ) {
          isStable = false;
        }

        if( !isStable )
          continue;

        TTtp->fill( tp );

        TrackingParticleRef tpref(TPCollection, i);
        auto TrkFound = simRecColl.find( tpref );
        if( TrkFound != simRecColl.end() ) {
          const auto& trkMatch = TrkFound->val;
          TTtp->fill_matchedTrk(
            trkMatch[0].first->pt(),
            trkMatch[0].first->eta(),
            trkMatch[0].first->phi(),
            trkMatch[0].first->charge(),
            trkMatch[0].second,
            trkMatch[0].first->numberOfValidHits()
          );
        }
        else {
          TTtp->fill_matchedTrk(
            -99999.,
            -99999.,
            -99999.,
            -99999,
            -99999.,
            -99999
          );
        }
      }
    }
  }
}


void MuonHLTNtupler::fill_trackTemplateMva(
  const edm::Event &iEvent,
  edm::EDGetTokenT<edm::View<reco::Track>>& theToken,
  edm::Handle<reco::TrackToTrackingParticleAssociator>& theAssociator_,
  edm::Handle<TrackingParticleCollection>& TPCollection_,
  edm::ESHandle<TrackerGeometry>& tracker,
  const pairSeedMvaEstimatorPhase2& pairSeedMvaEstimatorPhase2,
  std::map<tmpTSOD,unsigned int>& trkMap,
  trkTemplate* TTtrack,
  edm::ESHandle<MagneticField> magfieldH,
  const edm::EventSetup &iSetup,
  const GeometricSearchTracker& geomTracker
) {

  // edm::Handle<l1t::MuonBxCollection> h_L1Muon;
  // bool hasL1 = iEvent.getByToken( t_L1Muon_, h_L1Muon);

  edm::Handle<reco::RecoChargedCandidateCollection> h_L2Muon;
  bool hasL2 = iEvent.getByToken( t_L2Muon_, h_L2Muon );

  edm::Handle<l1t::TrackerMuonCollection> h_L1TkMu;
  bool hasL1TkMu = iEvent.getByToken( TkMuonToken_, h_L1TkMu);

  edm::Handle<edm::View<reco::Track>> trkHandle;
  if( iEvent.getByToken( theToken, trkHandle ) )
  {
    auto recSimColl = theAssociator_->associateRecoToSim(trkHandle,TPCollection_);

    for( unsigned int i = 0; i < trkHandle->size(); i++ )
    {
      TTtrack->fill(trkHandle->at(i));

      int linkNo = -1;
      for (unsigned int idxL3passed = 0; idxL3passed < iterL3IDpassed.size(); idxL3passed++) {
        if ( iterL3IDpassed.at(idxL3passed).isMatched(trkHandle->at(i)) ) linkNo = idxL3passed;
      }
      TTtrack->linkIterL3(linkNo);

      int linkNoNoId = -1;
      for (unsigned int idxL3passed = 0; idxL3passed < iterL3NoIDpassed.size(); idxL3passed++) {
        if ( iterL3NoIDpassed.at(idxL3passed).isMatched(trkHandle->at(i)) ) linkNoNoId = idxL3passed;
      }
      TTtrack->linkIterL3NoId(linkNoNoId);

      const PTrajectoryStateOnDet tmpseed = trkHandle->at(i).seedRef()->startingState();
      tmpTSOD tsod(tmpseed);
      trkMap.insert(make_pair(tsod,i));

      auto track = trkHandle->refAt(i);
      auto TPfound = recSimColl.find(track);
      if (TPfound != recSimColl.end()) {
        const auto& TPmatch = TPfound->val;
        TTtrack->fillBestTP(TPmatch[0].first);
        TTtrack->fillBestTPsharedFrac(TPmatch[0].second);
        TTtrack->fillmatchedTPsize(TPmatch.size());
      } else {  // sync vector size
        TTtrack->fillDummyTP();
        TTtrack->fillBestTPsharedFrac(-99999.);
        TTtrack->fillmatchedTPsize(0);
      }

      // if( doMVA && hasL1 && hasL2 && hasL1TkMu ) {
      if( doMVA && hasL2 && hasL1TkMu ) {
        const TrajectorySeed seed = *(trkHandle->at(i).seedRef());
        GlobalVector global_p = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().momentum());
        GlobalPoint  global_x = tracker->idToDet(seed.startingState().detId())->surface().toGlobal(seed.startingState().parameters().position());

        double mva = getSeedMva(
          pairSeedMvaEstimatorPhase2,
          seed,
          global_p,
          global_x,
          h_L1TkMu,
          magfieldH,
          iSetup,
          geomTracker
        );
        TTtrack->fillMva(mva);
      }
      else {
        TTtrack->fillMva( -99999. );
      }
    }
  }
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
  delete SThltIterL3OISeedsFromL2Muons;
  delete SThltIter0IterL3MuonPixelSeedsFromPixelTracks;
  delete SThltIter2IterL3MuonPixelSeeds;
  delete SThltIter3IterL3MuonPixelSeeds;
  delete SThltIter0IterL3FromL1MuonPixelSeedsFromPixelTracks;
  delete SThltIter2IterL3FromL1MuonPixelSeeds;
  delete SThltIter3IterL3FromL1MuonPixelSeeds;

  delete TThltIterL3OIMuonTrack;
  delete TThltIter0IterL3MuonTrack;
  delete TThltIter2IterL3MuonTrack;
  delete TThltIter3IterL3MuonTrack;
  delete TThltIter0IterL3FromL1MuonTrack;
  delete TThltIter2IterL3FromL1MuonTrack;
  delete TThltIter3IterL3FromL1MuonTrack;

  delete MTL3MuonsNoId;
  delete MTL3Muons;

  delete TrkParticle;

  delete VThltIterL3MuonTrimmedPixelVertices;
  delete VThltIterL3FromL1MuonTrimmedPixelVertices;

  for( unsigned int i = 0; i < trackCollectionNames_.size(); ++i) {
    delete trkTemplates_.at(i);
    delete tpTemplates_.at(i);
  }

  // if(doMVA) {
  //   delete mvaPhase2HltIter2IterL3FromL1MuonPixelSeeds_.first;
  //   delete mvaPhase2HltIter2IterL3FromL1MuonPixelSeeds_.second;
  // }
}

void MuonHLTNtupler::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}
void MuonHLTNtupler::endRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {}

DEFINE_FWK_MODULE(MuonHLTNtupler);
