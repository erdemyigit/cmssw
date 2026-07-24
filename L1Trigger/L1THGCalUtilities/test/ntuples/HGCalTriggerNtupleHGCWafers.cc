
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToMany.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCalUtilities/interface/HGCalTriggerNtupleBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"

#include "L1Trigger/L1THGCal/interface/concentrator/AEinputUtil.h"

class HGCalTriggerNtupleHGCWafers : public HGCalTriggerNtupleBase {
public:
  HGCalTriggerNtupleHGCWafers(const edm::ParameterSet& conf);
  ~HGCalTriggerNtupleHGCWafers() override{};
  void initialize(TTree&, const edm::ParameterSet&, edm::ConsumesCollector&&) final;
  void fill(const edm::Event& e, const HGCalTriggerNtupleEventSetup& es) final;

private:
  double calibrate(double, int, unsigned);
  void simhits(const edm::Event& e,
               std::unordered_map<uint32_t, double>& simhits_ee,
               std::unordered_map<uint32_t, double>& simhits_fh,
               std::unordered_map<uint32_t, double>& simhits_bh);
  float getsimenergy(const unsigned& i, 
        const std::unordered_map<uint32_t, double>& simhits_ee, 
        const std::unordered_map<uint32_t, double>& simhits_fh,
        const std::unordered_map<uint32_t, double>& simhits_bh);
  void clear() final;

  HGCalTriggerTools triggerTools_;

  edm::EDGetToken trigger_cells_token_;
  edm::EDGetToken simhits_ee_token_, simhits_fh_token_, simhits_bh_token_;

  unsigned bitsPerADC_;
  unsigned bitsPerNorm_;
  unsigned bitsPerCALQ_;
  unsigned bitsPerInput_;

  bool useModuleFactor_;
  bool bitShiftNormalize_;
  bool useTransverseADC_;
  bool normByMax_;

  AEinputUtil aeInputUtil_;

  double keV2fC_;
  std::vector<double> fcPerMip_;
  std::vector<double> layerWeights_;
  std::vector<double> thicknessCorrections_;

  int wafer_n_;
  std::vector<uint32_t> wafer_id_;
  std::vector<int> wafer_side_;
  std::vector<int> wafer_layer_;
  std::vector<int> wafer_subdet_;
  std::vector<int> wafer_waferu_;
  std::vector<int> wafer_waferv_;
  std::vector<int> wafer_wafertype_;

  std::vector<float> wafer_mipPt_;
  std::vector<float> wafer_energy_;
  std::vector<float> wafer_simenergy_;

  std::vector<float> wafer_eta_;
  std::vector<float> wafer_phi_;
  std::vector<float> wafer_x_;
  std::vector<float> wafer_y_;
  std::vector<float> wafer_z_;

  std::vector<uint64_t> wafer_sumCALQ_;

  std::vector<std::vector<uint32_t>> wafer_norm_;
  std::vector<std::vector<uint32_t>> wafer_ADC_;
  std::vector<std::vector<uint32_t>> wafer_CALQ_;
  std::vector<std::vector<float>> wafer_AEin_;

  static constexpr int nInputs_ = 64;
};

DEFINE_EDM_PLUGIN(HGCalTriggerNtupleFactory, HGCalTriggerNtupleHGCWafers, "HGCalTriggerNtupleHGCWafers");

float HGCalTriggerNtupleHGCWafers::getsimenergy(const unsigned& id, 
        const std::unordered_map<uint32_t, double>& simhits_ee, 
        const std::unordered_map<uint32_t, double>& simhits_fh,
        const std::unordered_map<uint32_t, double>& simhits_bh) {
    double energy = 0;
    unsigned layer = triggerTools_.layerWithOffset(id);
    // search for simhit for all the cells inside the trigger cell
    const auto& cells = triggerTools_.getTriggerGeometry()->getCellsFromTriggerCell(id);
    for (uint32_t c_id : cells) {
      int thickness = triggerTools_.thicknessIndex(c_id);
      if (triggerTools_.isEm(id)) {
        auto itr = simhits_ee.find(c_id);
        if (itr != simhits_ee.end()){
          energy += calibrate(itr->second, thickness, layer);
        }
      } else if (triggerTools_.isSilicon(id)) {
        auto itr = simhits_fh.find(c_id);
        if (itr != simhits_fh.end()){
          energy += calibrate(itr->second, thickness, layer);
        }
      } else {
        auto itr = simhits_bh.find(c_id);
        if (itr != simhits_bh.end()){
          energy += itr->second;
        }
      }
    }
    return energy;
}

HGCalTriggerNtupleHGCWafers::HGCalTriggerNtupleHGCWafers(const edm::ParameterSet& conf)
    : HGCalTriggerNtupleBase(conf),
      bitsPerADC_(conf.getParameter<unsigned>("bitsPerADC")),
      bitsPerNorm_(conf.getParameter<unsigned>("bitsPerNorm")),
      bitsPerCALQ_(conf.getParameter<unsigned>("bitsPerCALQ")),
      bitsPerInput_(conf.getParameter<unsigned>("bitsPerInput")),
      useModuleFactor_(conf.getParameter<bool>("useModuleFactor")),
      bitShiftNormalize_(conf.getParameter<bool>("bitShiftNormalize")),
      useTransverseADC_(conf.getParameter<bool>("useTransverseADC")),
      normByMax_(conf.getParameter<bool>("normByMax")),
      aeInputUtil_(bitsPerADC_, bitsPerNorm_, bitsPerCALQ_, bitsPerInput_, useModuleFactor_, bitShiftNormalize_, useTransverseADC_, normByMax_),
      keV2fC_(conf.getParameter<double>("keV2fC")),
      fcPerMip_(conf.getParameter<std::vector<double>>("fcPerMip")),
      layerWeights_(conf.getParameter<std::vector<double>>("layerWeights")),
      thicknessCorrections_(conf.getParameter<std::vector<double>>("thicknessCorrections")){
  accessEventSetup_ = false;
}

void HGCalTriggerNtupleHGCWafers::initialize(TTree& tree,
                                                   const edm::ParameterSet& conf,
                                                   edm::ConsumesCollector&& collector) {
  trigger_cells_token_ =
      collector.consumes<l1t::HGCalTriggerCellBxCollection>(conf.getParameter<edm::InputTag>("TriggerCells"));

  simhits_ee_token_ = collector.consumes<edm::PCaloHitContainer>(conf.getParameter<edm::InputTag>("eeSimHits"));
  simhits_fh_token_ = collector.consumes<edm::PCaloHitContainer>(conf.getParameter<edm::InputTag>("fhSimHits"));
  simhits_bh_token_ = collector.consumes<edm::PCaloHitContainer>(conf.getParameter<edm::InputTag>("bhSimHits"));

  std::string prefix(conf.getUntrackedParameter<std::string>("Prefix", "wafer"));

  std::string bname;
  auto withPrefix([&prefix, &bname](char const* vname) -> char const* {
    bname = prefix + "_" + vname;
    return bname.c_str();
  });

  tree.Branch(withPrefix("n"), &wafer_n_, (prefix + "_n/I").c_str());
  tree.Branch(withPrefix("id"), &wafer_id_);
  tree.Branch(withPrefix("zside"), &wafer_side_);
  tree.Branch(withPrefix("layer"), &wafer_layer_);
  tree.Branch(withPrefix("subdet"), &wafer_subdet_);
  tree.Branch(withPrefix("waferu"), &wafer_waferu_);
  tree.Branch(withPrefix("waferv"), &wafer_waferv_);
  tree.Branch(withPrefix("wafertype"), &wafer_wafertype_);
  tree.Branch(withPrefix("mipPt"), &wafer_mipPt_);
  tree.Branch(withPrefix("energy"), &wafer_energy_);
  tree.Branch(withPrefix("simenergy"), &wafer_simenergy_);
  tree.Branch(withPrefix("eta"), &wafer_eta_);
  tree.Branch(withPrefix("phi"), &wafer_phi_);
  tree.Branch(withPrefix("x"), &wafer_x_);
  tree.Branch(withPrefix("y"), &wafer_y_);
  tree.Branch(withPrefix("z"), &wafer_z_);
  tree.Branch(withPrefix("sumCALQ"), &wafer_sumCALQ_);

  wafer_CALQ_.resize(nInputs_);
  wafer_AEin_.resize(nInputs_);
  wafer_ADC_.resize(nInputs_);
  wafer_norm_.resize(nInputs_);
  for(unsigned i=0; i<nInputs_; ++i){
    tree.Branch(withPrefix(("CALQ"+std::to_string(i)).c_str()), &wafer_CALQ_[i]);
    tree.Branch(withPrefix(("AEin"+std::to_string(i)).c_str()), &wafer_AEin_[i]);
    tree.Branch(withPrefix(("ADC"+std::to_string(i)).c_str()), &wafer_ADC_[i]);
    tree.Branch(withPrefix(("norm"+std::to_string(i)).c_str()), &wafer_norm_[i]);
  }
}

void HGCalTriggerNtupleHGCWafers::fill(const edm::Event& e, const HGCalTriggerNtupleEventSetup& es) {
  // retrieve trigger cells
  edm::Handle<l1t::HGCalTriggerCellBxCollection> trigger_cells_h;
  e.getByToken(trigger_cells_token_, trigger_cells_h);
  const l1t::HGCalTriggerCellBxCollection& trigger_cells = *trigger_cells_h;

  // sim hit association
  std::unordered_map<uint32_t, double> simhits_ee;
  std::unordered_map<uint32_t, double> simhits_fh;
  std::unordered_map<uint32_t, double> simhits_bh;
  std::unordered_map<uint32_t, bool> simhit_has_rechit;
  simhits(e, simhits_ee, simhits_fh, simhits_bh);

  triggerTools_.setGeometry(es.geometry.product());
  aeInputUtil_.setGeometry(es.geometry.product());

  clear();

  std::unordered_map<uint32_t, std::vector<l1t::HGCalTriggerCell>> wafers;
  for(const auto& tc : trigger_cells){
    if (triggerTools_.isScintillator(tc.detId())) continue;
    uint32_t wafer_id = triggerTools_.getTriggerGeometry()->getModuleFromTriggerCell(tc.detId());
    wafers[wafer_id].push_back(tc);
  }

  for(const auto& wafer : wafers){
      ++wafer_n_;

      DetId id0(wafer.second[0].detId());

      wafer_id_.emplace_back(wafer.first);
      wafer_side_.emplace_back(triggerTools_.zside(id0));
      wafer_layer_.emplace_back(triggerTools_.layerWithOffset(id0));
      if(id0.det() == DetId::HGCalTrigger){
          HGCalTriggerDetId idtrg0(id0);
          wafer_subdet_.emplace_back(idtrg0.subdet());
          wafer_waferu_.emplace_back(idtrg0.waferU());
          wafer_waferv_.emplace_back(idtrg0.waferV());
          wafer_wafertype_.emplace_back(idtrg0.type());
      //} else if(id0.det() == DetId::HGCalHSc){
      //    HGCScintillatorDetId idsc0(id0);
      //    wafer_subdet_.emplace_back(idsc0.subdet());
      //    wafer_waferu_.emplace_back(-999);
      //    wafer_waferv_.emplace_back(-999);
      //    wafer_wafertype_.emplace_back(idsc0.type());
      } else {
          throw cms::Exception("HGCalTriggerNtupleWafers") << "Unknown detid " << id0.det() << std::endl;
      }

      float mipPt=0, energy=0, simenergy=0;
      for(const auto& tc : wafer.second){
          mipPt += tc.mipPt();
          energy += tc.energy();
          simenergy += getsimenergy(tc.detId(), simhits_ee, simhits_fh, simhits_bh);
      }
      wafer_mipPt_.emplace_back(mipPt);
      wafer_energy_.emplace_back(energy);
      wafer_simenergy_.emplace_back(simenergy);

      auto pos = triggerTools_.getTriggerGeometry()->getModulePosition(wafer.first);
      wafer_eta_.emplace_back(pos.eta());
      wafer_phi_.emplace_back(pos.phi());
      wafer_x_.emplace_back(pos.x());
      wafer_y_.emplace_back(pos.y());
      wafer_z_.emplace_back(pos.z());

      aeInputUtil_.run(wafer.second);

      wafer_sumCALQ_.emplace_back(aeInputUtil_.getModSum());

      for(unsigned i=0; i<nInputs_; ++i){
          wafer_ADC_[i].emplace_back(aeInputUtil_.getADC(i));
          wafer_norm_[i].emplace_back(aeInputUtil_.getNorm(i));
          wafer_CALQ_[i].emplace_back(aeInputUtil_.getCALQ(i));
          wafer_AEin_[i].emplace_back(
                  aeInputUtil_.getInput(i)/aeInputUtil_.getInputNorm());
      }

      /*
      printf("\n------------------------------------\n");
      printf("ADCs\n");
      unsigned i=0;
      for(unsigned u=0; u<8; ++u){
          for(unsigned v=0; v<8; ++v){
            printf("%03u ", aeInputUtil_.getADC(i));
            ++i;
          }
          printf("\n");
      }
      printf("\n");

      printf("norms\n");
      i=0;
      for(unsigned u=0; u<8; ++u){
          for(unsigned v=0; v<8; ++v){
              printf("%03u ", aeInputUtil_.getNorm(i));
              ++i;
          }
          printf("\n");
      }
      printf("\n");

      printf("CALQ\n");
      i=0;
      for(unsigned u=0; u<8; ++u){
          for(unsigned v=0; v<8; ++v){
              printf("%03u ", aeInputUtil_.getCALQ(i));
              ++i;
          }
          printf("\n");
      }
      printf("\n");

      printf("AEin\n");
      i=0;
      for(unsigned u=0; u<8; ++u){
          for(unsigned v=0; v<8; ++v){
              printf("%0.3f ", aeInputUtil_.getInput(i)/aeInputUtil_.getInputNorm());
              ++i;
          }
          printf("\n");
      }
      printf("\n------------------------------------\n");*/
  }
}

double HGCalTriggerNtupleHGCWafers::calibrate(double energy, int thickness, unsigned layer) {
  double fcPerMip = fcPerMip_[thickness];
  double thicknessCorrection = thicknessCorrections_[thickness];
  double layerWeight = layerWeights_[layer];
  double TeV2GeV = 1.e3;
  return energy * keV2fC_ / fcPerMip * layerWeight * TeV2GeV / thicknessCorrection;
}

void HGCalTriggerNtupleHGCWafers::simhits(const edm::Event& e,
                                                std::unordered_map<uint32_t, double>& simhits_ee,
                                                std::unordered_map<uint32_t, double>& simhits_fh,
                                                std::unordered_map<uint32_t, double>& simhits_bh) {
  edm::Handle<edm::PCaloHitContainer> ee_simhits_h;
  e.getByToken(simhits_ee_token_, ee_simhits_h);
  const edm::PCaloHitContainer& ee_simhits = *ee_simhits_h;
  edm::Handle<edm::PCaloHitContainer> fh_simhits_h;
  e.getByToken(simhits_fh_token_, fh_simhits_h);
  const edm::PCaloHitContainer& fh_simhits = *fh_simhits_h;
  edm::Handle<edm::PCaloHitContainer> bh_simhits_h;
  e.getByToken(simhits_bh_token_, bh_simhits_h);
  const edm::PCaloHitContainer& bh_simhits = *bh_simhits_h;

  // EE
  for (const auto& simhit : ee_simhits) {
    DetId id = simhit.id();
    if (id.rawId() == 0)
      continue;
    auto itr_insert = simhits_ee.emplace(id, 0.);
    itr_insert.first->second += simhit.energy();
  }
  //  FH
  for (const auto& simhit : fh_simhits) {
    DetId id = simhit.id();
    if (id.rawId() == 0)
      continue;
    auto itr_insert = simhits_fh.emplace(id, 0.);
    itr_insert.first->second += simhit.energy();
  }
  //  BH
  for (const auto& simhit : bh_simhits) {
    DetId id = simhit.id();
    if (id.rawId() == 0)
      continue;
    auto itr_insert = simhits_bh.emplace(id, 0.);
    itr_insert.first->second += simhit.energy();
  }
}

void HGCalTriggerNtupleHGCWafers::clear() {
  wafer_n_ = 0;

  wafer_id_.clear();
  wafer_side_.clear();
  wafer_layer_.clear();
  wafer_subdet_.clear();
  wafer_waferu_.clear();
  wafer_waferv_.clear();
  wafer_wafertype_.clear();

  wafer_mipPt_.clear();
  wafer_energy_.clear();
  wafer_simenergy_.clear();

  wafer_eta_.clear();
  wafer_phi_.clear();
  wafer_x_.clear();
  wafer_y_.clear();
  wafer_z_.clear();

  for(unsigned i=0; i<nInputs_; ++i){
      wafer_CALQ_[i].clear();
      wafer_AEin_[i].clear();
      wafer_ADC_[i].clear();
      wafer_norm_[i].clear();
  }
}
