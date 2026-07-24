#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

namespace l1t {
  class HGC3DClusterDummySelector : public edm::stream::EDProducer<> {
  public:
    explicit HGC3DClusterDummySelector(const edm::ParameterSet &);
    ~HGC3DClusterDummySelector() override {}

  private:
    edm::EDGetTokenT<l1t::HGCalMulticlusterBxCollection> src_;

    void produce(edm::Event &, const edm::EventSetup &) override;

  };  // class
}  // namespace l1t

l1t::HGC3DClusterDummySelector::HGC3DClusterDummySelector(const edm::ParameterSet &iConfig)
    : src_(consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getParameter<edm::InputTag>("src"))){
  produces<l1t::HGCalMulticlusterBxCollection>();
}

void l1t::HGC3DClusterDummySelector::produce(edm::Event &iEvent, const edm::EventSetup &) {
  auto out = std::make_unique<l1t::HGCalMulticlusterBxCollection>();

  edm::Handle<l1t::HGCalMulticlusterBxCollection> multiclusters;
  iEvent.getByToken(src_, multiclusters);

  for (int bx = multiclusters->getFirstBX(); bx <= multiclusters->getLastBX(); ++bx) {
    for (auto it = multiclusters->begin(bx), ed = multiclusters->end(bx); it != ed; ++it){
      const auto &multicluster = *it;
      out->push_back(bx, multicluster);
    }
  }

  iEvent.put(std::move(out));
}
using l1t::HGC3DClusterDummySelector;
DEFINE_FWK_MODULE(HGC3DClusterDummySelector);
