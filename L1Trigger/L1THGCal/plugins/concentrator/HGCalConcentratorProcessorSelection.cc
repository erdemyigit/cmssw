#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorProcessorSelection.h"
#include <limits>
#include <set>

#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

DEFINE_EDM_PLUGIN(HGCalConcentratorFactory, HGCalConcentratorProcessorSelection, "HGCalConcentratorProcessorSelection");

HGCalConcentratorProcessorSelection::HGCalConcentratorProcessorSelection(const edm::ParameterSet& conf)
    : HGCalConcentratorProcessorBase(conf),
      fixedDataSizePerHGCROC_(conf.getParameter<bool>("fixedDataSizePerHGCROC")),
      allTrigCellsInTrigSums_(conf.getParameter<bool>("allTrigCellsInTrigSums")),
      coarsenTriggerCells_(conf.getParameter<std::vector<unsigned>>("coarsenTriggerCells")),
      selectionType_(kNSubDetectors_){
  std::vector<std::string> selectionType(conf.getParameter<std::vector<std::string>>("Method"));
  if (selectionType.size() != kNSubDetectors_ || coarsenTriggerCells_.size() != kNSubDetectors_) {
    throw cms::Exception("HGCTriggerParameterError")
        << "Inconsistent number of sub-detectors (should be " << kNSubDetectors_ << ")";
  }

  for (int subdet = 0; subdet < kNSubDetectors_; subdet++) {
    if (selectionType[subdet] == "thresholdSelect") {
      selectionType_[subdet] = thresholdSelect;
      if (!thresholdImpl_)
        thresholdImpl_ = std::make_unique<HGCalConcentratorThresholdImpl>(conf);
      if (!trigSumImpl_)
        trigSumImpl_ = std::make_unique<HGCalConcentratorTrigSumImpl>(conf);
    } else if (selectionType[subdet] == "bestChoiceSelect") {
      selectionType_[subdet] = bestChoiceSelect;
      if (!bestChoiceImpl_)
        bestChoiceImpl_ = std::make_unique<HGCalConcentratorBestChoiceImpl>(conf);
      if (!trigSumImpl_)
        trigSumImpl_ = std::make_unique<HGCalConcentratorTrigSumImpl>(conf);
    } else if (selectionType[subdet] == "superTriggerCellSelect") {
      selectionType_[subdet] = superTriggerCellSelect;
      if (!superTriggerCellImpl_)
        superTriggerCellImpl_ = std::make_unique<HGCalConcentratorSuperTriggerCellImpl>(conf);
      if (conf.existsAs<bool>("stcLinkAware") && conf.getParameter<bool>("stcLinkAware") && !stcFormatImpl_[0]) {
        stcLinkAware_ = true;
        stcLinkAwareDebug_ = conf.getUntrackedParameter<bool>("stcLinkAwareDebug", false);
        const size_t nsize = conf.getParameter<std::vector<unsigned>>("stcSize").size();
        // {STC size, exponent bits, mantissa bits, energy division}; no LSB dropped for STC/CTC values (spec sec. 13.4)
        const std::array<std::tuple<unsigned, uint32_t, uint32_t, std::string>, kNStcFormats> formats = {
            {{4, 5, 4, "superTriggerCell"}, {4, 4, 3, "superTriggerCell"}, {4, 4, 3, "equalShare"},
             {16, 5, 4, "superTriggerCell"}}};
        for (unsigned f = 0; f < kNStcFormats; f++) {
          edm::ParameterSet pf(conf);
          pf.addParameter<std::vector<unsigned>>("stcSize", std::vector<unsigned>(nsize, std::get<0>(formats[f])));
          pf.addParameter<std::string>("type_energy_division", std::get<3>(formats[f]));
          edm::ParameterSet comp(conf.getParameterSet("superTCCompression"));
          comp.addParameter<uint32_t>("exponentBits", std::get<1>(formats[f]));
          comp.addParameter<uint32_t>("mantissaBits", std::get<2>(formats[f]));
          comp.addParameter<uint32_t>("truncationBits", 0);
          pf.addParameter<edm::ParameterSet>("superTCCompression", comp);
          stcFormatImpl_[f] = std::make_unique<HGCalConcentratorSuperTriggerCellImpl>(pf);
        }
      }
    } else if (selectionType[subdet] == "autoEncoder") {
      selectionType_[subdet] = autoEncoderSelect;
      if (!autoEncoderImpl_)
        autoEncoderImpl_ = std::make_unique<HGCalConcentratorAutoEncoderImpl>(conf);
    } else if (selectionType[subdet] == "noSelection") {
      selectionType_[subdet] = noSelection;
    } else {
      throw cms::Exception("HGCTriggerParameterError")
          << "Unknown type of concentrator selection '" << selectionType[subdet] << "'";
    }
  }

  if (std::find(coarsenTriggerCells_.begin(), coarsenTriggerCells_.end(), true) != coarsenTriggerCells_.end() ||
      fixedDataSizePerHGCROC_) {
    coarsenerImpl_ = std::make_unique<HGCalConcentratorCoarsenerImpl>(conf);
  }
}

void HGCalConcentratorProcessorSelection::run(const edm::Handle<l1t::HGCalTriggerCellBxCollection>& triggerCellCollInput,
                                              std::tuple<l1t::HGCalTriggerCellBxCollection,
                                                         l1t::HGCalTriggerSumsBxCollection,
                                                         l1t::HGCalConcentratorDataBxCollection>& triggerCollOutput) {
  if (thresholdImpl_)
    thresholdImpl_->setGeometry(geometry());
  if (bestChoiceImpl_)
    bestChoiceImpl_->setGeometry(geometry());
  if (superTriggerCellImpl_)
    superTriggerCellImpl_->setGeometry(geometry());
  for (auto& impl : stcFormatImpl_)
    if (impl)
      impl->setGeometry(geometry());
  if (autoEncoderImpl_)
    autoEncoderImpl_->setGeometry(geometry());
  if (coarsenerImpl_)
    coarsenerImpl_->setGeometry(geometry());
  if (trigSumImpl_)
    trigSumImpl_->setGeometry(geometry());
  
  triggerTools_.setGeometry(geometry());

  auto& triggerCellCollOutput = std::get<0>(triggerCollOutput);
  auto& triggerSumCollOutput = std::get<1>(triggerCollOutput);
  auto& autoEncoderCollOutput = std::get<2>(triggerCollOutput);

  const l1t::HGCalTriggerCellBxCollection& collInput = *triggerCellCollInput;

  std::unordered_map<uint32_t, std::vector<l1t::HGCalTriggerCell>> tc_modules;
  for (const auto& trigCell : collInput) {
    uint32_t module = geometry()->getModuleFromTriggerCell(trigCell.detId());
    tc_modules[module].push_back(trigCell);
  }

  for (const auto& module_trigcell : tc_modules) {
    std::vector<l1t::HGCalTriggerCell> trigCellVecOutput;
    std::vector<l1t::HGCalTriggerCell> trigCellVecCoarsened;
    std::vector<l1t::HGCalTriggerCell> trigCellVecNotSelected;
    std::vector<l1t::HGCalTriggerSums> trigSumsVecOutput;
    std::vector<l1t::HGCalConcentratorData> ae_EncodedLayerOutput;

    int thickness = triggerTools_.thicknessIndex(module_trigcell.second.at(0).detId());

    HGCalTriggerTools::SubDetectorType subdet = triggerTools_.getSubDetectorType(module_trigcell.second.at(0).detId());

    if (coarsenTriggerCells_[subdet] || (fixedDataSizePerHGCROC_ && thickness > kHighDensityThickness_)) {
      coarsenerImpl_->coarsen(module_trigcell.second, trigCellVecCoarsened);

      switch (selectionType_[subdet]) {
        case thresholdSelect:
          thresholdImpl_->select(trigCellVecCoarsened, trigCellVecOutput, trigCellVecNotSelected);
          break;
        case bestChoiceSelect:
          if (triggerTools_.isEm(module_trigcell.first)) {
            bestChoiceImpl_->select(geometry()->getLinksInModule(module_trigcell.first),
                                    geometry()->getModuleSize(module_trigcell.first),
                                    module_trigcell.second,
                                    trigCellVecOutput,
                                    trigCellVecNotSelected);
          } else {
            bestChoiceImpl_->select(geometry()->getLinksInModule(module_trigcell.first),
                                    geometry()->getModuleSize(module_trigcell.first),
                                    trigCellVecCoarsened,
                                    trigCellVecOutput,
                                    trigCellVecNotSelected);
          }
          break;
        case superTriggerCellSelect:
          if (stcLinkAware_)
            selectSuperTriggerCellsLinkAware(module_trigcell.first, trigCellVecCoarsened, trigCellVecOutput);
          else
            superTriggerCellImpl_->select(trigCellVecCoarsened, trigCellVecOutput);
          break;
        case autoEncoderSelect:
          autoEncoderImpl_->select(geometry()->getLinksInModule(module_trigcell.first),
                                   trigCellVecCoarsened,
                                   trigCellVecOutput,
                                   ae_EncodedLayerOutput);
          break;
        case noSelection:
          trigCellVecOutput = trigCellVecCoarsened;
          break;
        default:
          // Should not happen, selection type checked in constructor
          break;
      }

    } else {
      switch (selectionType_[subdet]) {
        case thresholdSelect:
          thresholdImpl_->select(module_trigcell.second, trigCellVecOutput, trigCellVecNotSelected);
          break;
        case bestChoiceSelect:
          bestChoiceImpl_->select(geometry()->getLinksInModule(module_trigcell.first),
                                  geometry()->getModuleSize(module_trigcell.first),
                                  module_trigcell.second,
                                  trigCellVecOutput,
                                  trigCellVecNotSelected);
          break;
        case superTriggerCellSelect:
          if (stcLinkAware_)
            selectSuperTriggerCellsLinkAware(module_trigcell.first, module_trigcell.second, trigCellVecOutput);
          else
            superTriggerCellImpl_->select(module_trigcell.second, trigCellVecOutput);
          break;
        case autoEncoderSelect:
          autoEncoderImpl_->select(geometry()->getLinksInModule(module_trigcell.first),
                                   module_trigcell.second,
                                   trigCellVecOutput,
                                   ae_EncodedLayerOutput);
          break;
        case noSelection:
          trigCellVecOutput = module_trigcell.second;
          break;
        default:
          // Should not happen, selection type checked in constructor
          break;
      }
    }

    // trigger sum
    if (trigSumImpl_) {
      if (allTrigCellsInTrigSums_) {  // using all TCs
        trigSumImpl_->doSum(module_trigcell.first, module_trigcell.second, trigSumsVecOutput);
      } else {  // using only unselected TCs
        trigSumImpl_->doSum(module_trigcell.first, trigCellVecNotSelected, trigSumsVecOutput);
      }
    }

    for (const auto& trigCell : trigCellVecOutput) {
      triggerCellCollOutput.push_back(0, trigCell);
    }
    for (const auto& trigSums : trigSumsVecOutput) {
      triggerSumCollOutput.push_back(0, trigSums);
    }
    for (const auto& aeVal : ae_EncodedLayerOutput) {
      autoEncoderCollOutput.push_back(0, aeVal);
    }
  }
}

void HGCalConcentratorProcessorSelection::selectSuperTriggerCellsLinkAware(
    uint32_t module,
    const std::vector<l1t::HGCalTriggerCell>& trigCellVecInput,
    std::vector<l1t::HGCalTriggerCell>& trigCellVecOutput) {
  if (trigCellVecInput.empty())
    return;
  const unsigned nLinks = geometry()->getLinksInModule(module);
  if (nLinks == 0)  // no active eTx: nothing is transmitted
    return;
  const unsigned etx = std::min(nLinks, kMaxEtx_);
  const uint32_t tc0 = trigCellVecInput.at(0).detId();
  const bool highDensity =
      triggerTools_.isSilicon(tc0) && triggerTools_.thicknessIndex(tc0) == kHighDensityThickness_;
  // high density: finest format that fits every sum of the module; low density and scintillator: STC16 (nominal)
  static const std::vector<unsigned> candidatesHD = {kSTC4_5E4M, kSTC4_4E3M, kCTC4_4E3M, kSTC16_5E4M};
  static const std::vector<unsigned> candidatesLD = {kSTC16_5E4M};
  const auto& candidates = highDensity ? candidatesHD : candidatesLD;
  const auto moduleTCs = geometry()->getTriggerCellsFromModule(module);

  unsigned chosen = candidates.back();
  std::vector<uint32_t> ids;
  for (unsigned f : candidates) {
    std::set<uint32_t> uniq;
    for (uint32_t tc : moduleTCs)
      uniq.insert(stcFormatImpl_[f]->stcId(tc));
    ids.assign(uniq.begin(), uniq.end());
    chosen = f;
    if (kMaxSTC_[f][etx - 1] >= ids.size())
      break;  // fits (if none does, the last candidate is used with truncation)
  }
  const size_t nBefore = trigCellVecOutput.size();
  stcFormatImpl_[chosen]->select(trigCellVecInput, trigCellVecOutput, kMaxSTC_[chosen][etx - 1], ids);
  if (stcLinkAwareDebug_) {
    std::set<uint32_t> sent;
    for (size_t i = nBefore; i < trigCellVecOutput.size(); i++)
      sent.insert(stcFormatImpl_[chosen]->stcId(trigCellVecOutput[i].detId()));
    edm::LogVerbatim("STCLINK") << "STCLINK module=" << module << " nlinks=" << nLinks << " hd=" << highDensity
                                << " format=" << chosen << " units=" << ids.size()
                                << " cap=" << kMaxSTC_[chosen][etx - 1] << " sums_sent=" << sent.size()
                                << " tcs_out=" << trigCellVecOutput.size() - nBefore;
  }
}
