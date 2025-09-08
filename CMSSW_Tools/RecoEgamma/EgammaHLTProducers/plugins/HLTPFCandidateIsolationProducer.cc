#include <iostream>
#include <vector>
#include <memory>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"

template <typename T1>
class HLTPFCandidateIsolationProducer : public edm::stream::EDProducer<> {
  typedef std::vector<T1> T1Collection;
  typedef edm::Ref<T1Collection> T1Ref;
  typedef edm::AssociationMap<edm::OneToValue<std::vector<T1>, float>> T1IsolationMap;

public:
  explicit HLTPFCandidateIsolationProducer(const edm::ParameterSet&);
  ~HLTPFCandidateIsolationProducer() override = default;

  void produce(edm::Event&, const edm::EventSetup&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  //// To change
  edm::EDGetTokenT<T1Collection> recoCandidateProducer_;
  const edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidateProducer_;
  //const edm::EDGetTokenT<double> rhoProducer_;

  const double drMax_;
  const double drVeto_;
  const double drVetoCh_;
  const double minEnergy_;
  //const bool doRhoCorrection_;
  //const double rhoMax_;
  //const double rhoScale_;
  //const std::vector<double> effectiveAreas_;
};

template <typename T1>
HLTPFCandidateIsolationProducer<T1>::HLTPFCandidateIsolationProducer(const edm::ParameterSet& config)
    : pfCandidateProducer_(consumes<reco::PFCandidateCollection>(config.getParameter<edm::InputTag>("pfCandidateProducer"))),
      //rhoProducer_(consumes<double>(config.getParameter<edm::InputTag>("rhoProducer"))),
      drMax_(config.getParameter<double>("drMax")),
      drVeto_(config.getParameter<double>("drVeto")),
      drVetoCh_(config.getParameter<double>("drVetoCh")),
      minEnergy_(config.getParameter<double>("minEnergy"))
      //doRhoCorrection_(config.getParameter<bool>("doRhoCorrection")),
      //rhoMax_(config.getParameter<double>("rhoMax")),
      //rhoScale_(config.getParameter<double>("rhoScale")),
      //effectiveAreas_(config.getParameter<std::vector<double>>("effectiveAreas"))
{
  //  if (doRhoCorrection_) {
  //  if (effectiveAreas_.size() != 2)
  //    throw cms::Exception("IncompatibleVects")
  //        << "effectiveAreas should have two elements for em and had components. \n";

  std::string recoCandidateProducerName = "recoCandidateProducer";
  if ((typeid(HLTPFCandidateIsolationProducer<T1>) ==
       typeid(HLTPFCandidateIsolationProducer<reco::RecoEcalCandidate>)))
    recoCandidateProducerName = "recoEcalCandidateProducer";

  recoCandidateProducer_ = consumes<T1Collection>(config.getParameter<edm::InputTag>(recoCandidateProducerName));
  produces<T1IsolationMap>();
}

template <typename T1>
void HLTPFCandidateIsolationProducer<T1>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  std::string recoCandidateProducerName = "recoCandidateProducer";
  if ((typeid(HLTPFCandidateIsolationProducer<T1>) ==
       typeid(HLTPFCandidateIsolationProducer<reco::RecoEcalCandidate>)))
    recoCandidateProducerName = "recoEcalCandidateProducer";

  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>(recoCandidateProducerName, edm::InputTag("hltL1SeededRecoEcalCandidatePF"));
  desc.add<edm::InputTag>("pfCandidateProducer", edm::InputTag("hltParticleFlowTmp"));
  // desc.add<edm::InputTag>("rhoProducer", edm::InputTag("fixedGridRhoFastjetAllCalo"));
  // desc.add<bool>("doRhoCorrection", false);
  // desc.add<double>("rhoMax", 9.9999999E7);
  // desc.add<double>("rhoScale", 1.0);
  desc.add<double>("drMax", 0.4);
  desc.add<double>("drVeto", 0.01);
  desc.add<double>("drVetoCh", 0.0001);
  desc.add<double>("minEnergy", 0.0);
  //desc.add<std::vector<double>>("effectiveAreas", {0.0, 0.0});  // for em and had components
  descriptions.add(defaultModuleLabel<HLTPFCandidateIsolationProducer<T1>>(), desc);
}

template <typename T1>
void HLTPFCandidateIsolationProducer<T1>::produce(edm::Event& iEvent, const edm::EventSetup&) {

  // edm::Handle<double> rhoHandle;
  // double rho = 0.0;
  // if (doRhoCorrection_) {
  //   iEvent.getByToken(rhoProducer_, rhoHandle);
  //   rho = *(rhoHandle.product());
  // }
  //
  // rho = std::min(rho, rhoMax_);
  // rho = rho * rhoScale_;

  edm::Handle<T1Collection> recoCandHandle;
  edm::Handle<reco::PFCandidateCollection> pfCandidateHandle;

  iEvent.getByToken(recoCandidateProducer_, recoCandHandle);
  iEvent.getByToken(pfCandidateProducer_, pfCandidateHandle);

  const std::vector<reco::PFCandidate> pfCandidates = *(pfCandidateHandle.product());

  T1IsolationMap recoCandMap(recoCandHandle);
  for (unsigned int iReco = 0; iReco < recoCandHandle->size(); iReco++) {
    T1Ref candRef(recoCandHandle, iReco);

    float sum = 0.0;
    for (unsigned int iPF = 0; iPF < pfCandidateHandle->size(); iPF++) {
      reco::PFCandidateRef pc(pfCandidateHandle, iPF);
      
      float dr2 = reco::deltaR2(candRef->eta(), candRef->phi(), pc->eta(), pc->phi());
      if (dr2 > drMax_ * drMax_)
	continue;
      if (fabs(pc->charge())){
	if (dr2 > drVetoCh_ && pc->pt() > minEnergy_)
	  sum += pc->pt();
      }else{
	if (dr2 > drVeto_ && pc->pt() > minEnergy_)
          sum += pc->pt();
      }
    }
    recoCandMap.insert(candRef, sum);
  }

  iEvent.put(std::make_unique<T1IsolationMap>(recoCandMap));
}

typedef HLTPFCandidateIsolationProducer<reco::RecoEcalCandidate> EgammaHLTPFCandidateIsolationProducer;
typedef HLTPFCandidateIsolationProducer<reco::RecoChargedCandidate> MuonHLTPFCandidateIsolationProducer;

DEFINE_FWK_MODULE(EgammaHLTPFCandidateIsolationProducer);
DEFINE_FWK_MODULE(MuonHLTPFCandidateIsolationProducer);
