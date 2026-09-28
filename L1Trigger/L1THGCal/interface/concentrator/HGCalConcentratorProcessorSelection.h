#ifndef __L1Trigger_L1THGCal_HGCalConcentratorProcessorSelection_h__
#define __L1Trigger_L1THGCal_HGCalConcentratorProcessorSelection_h__

#include "L1Trigger/L1THGCal/interface/HGCalProcessorBase.h"
#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorThresholdImpl.h"
#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorBestChoiceImpl.h"
#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorSuperTriggerCellImpl.h"
#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorCoarsenerImpl.h"
#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorTrigSumImpl.h"
#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorAutoEncoderImpl.h"

#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerSums.h"
#include "DataFormats/L1THGCal/interface/HGCalConcentratorData.h"

#include <array>
#include <utility>
#include <tuple>

class HGCalConcentratorProcessorSelection : public HGCalConcentratorProcessorBase {
private:
  enum SelectionType { thresholdSelect, bestChoiceSelect, superTriggerCellSelect, autoEncoderSelect, noSelection };

public:
  HGCalConcentratorProcessorSelection(const edm::ParameterSet& conf);

  void run(const edm::Handle<l1t::HGCalTriggerCellBxCollection>& triggerCellCollInput,
           std::tuple<l1t::HGCalTriggerCellBxCollection,
                      l1t::HGCalTriggerSumsBxCollection,
                      l1t::HGCalConcentratorDataBxCollection>& triggerCollOutput) override;

private:
  bool fixedDataSizePerHGCROC_;
  bool allTrigCellsInTrigSums_;
  std::vector<unsigned> coarsenTriggerCells_;
  static constexpr int kHighDensityThickness_ = 0;
  static constexpr int kNSubDetectors_ = 3;

  std::vector<SelectionType> selectionType_;

  std::unique_ptr<HGCalConcentratorThresholdImpl> thresholdImpl_;
  std::unique_ptr<HGCalConcentratorBestChoiceImpl> bestChoiceImpl_;
  std::unique_ptr<HGCalConcentratorSuperTriggerCellImpl> superTriggerCellImpl_;
  std::unique_ptr<HGCalConcentratorCoarsenerImpl> coarsenerImpl_;
  std::unique_ptr<HGCalConcentratorTrigSumImpl> trigSumImpl_;
  std::unique_ptr<HGCalConcentratorAutoEncoderImpl> autoEncoderImpl_;

  // eLink-aware super trigger cells (ECON-T spec v10, sec. 14.2, Fig. 13): per module, the STC format is chosen from the
  // module's number of eLinks (active eTx) so that the transmitted data fit; enabled by the optional bool stcLinkAware.
  enum StcFormat { kSTC4_5E4M = 0, kSTC4_4E3M, kCTC4_4E3M, kSTC16_5E4M, kNStcFormats };
  static constexpr unsigned kMaxEtx_ = 5;
  static constexpr unsigned kMaxSTC_[kNStcFormats][kMaxEtx_] = {
      {2, 5, 8, 11, 14},   // STC4 5E+4M: 9-bit sum + 2-bit address
      {3, 6, 10, 13, 17},  // STC4 4E+3M: 7-bit sum + 2-bit address
      {4, 8, 13, 17, 22},  // CTC4 4E+3M: 7-bit sum, no address
      {2, 4, 7, 9, 12}};   // STC16 5E+4M: 9-bit sum + 4-bit address
  bool stcLinkAware_ = false;
  std::array<std::unique_ptr<HGCalConcentratorSuperTriggerCellImpl>, kNStcFormats> stcFormatImpl_;
  void selectSuperTriggerCellsLinkAware(uint32_t module,
                                        const std::vector<l1t::HGCalTriggerCell>& trigCellVecInput,
                                        std::vector<l1t::HGCalTriggerCell>& trigCellVecOutput);

  HGCalTriggerTools triggerTools_;
};

#endif
