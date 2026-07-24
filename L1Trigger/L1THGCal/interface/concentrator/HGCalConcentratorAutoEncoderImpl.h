#ifndef __L1Trigger_L1THGCal_HGCalConcentratorAutoEncoderImpl_h__
#define __L1Trigger_L1THGCal_HGCalConcentratorAutoEncoderImpl_h__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalConcentratorData.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
#include <vector>

#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#include "L1Trigger/L1THGCal/interface/concentrator/AEinputUtil.h"

class HGCalConcentratorAutoEncoderImpl {
public:
  HGCalConcentratorAutoEncoderImpl(const edm::ParameterSet& conf);

  void select(unsigned nLinks,
              const std::vector<l1t::HGCalTriggerCell>& trigCellVecInput,
              std::vector<l1t::HGCalTriggerCell>& trigCellVecOutput,
              std::vector<l1t::HGCalConcentratorData>& ae_EncodedOutput);

  void setGeometry(const HGCalTriggerGeometryBase* const geom) { 
      triggerTools_.setGeometry(geom); 
      aeInputUtil_.setGeometry(geom);
  }

  inline const HGCalTriggerGeometryBase* geometry() const {
      return triggerTools_.getTriggerGeometry();
  }

private:
  static constexpr int nEncodedLayerNodes_ = 16;

  static constexpr unsigned int maxNumberOfLinks_ = 13;

  static constexpr int encoderTensorDims_ = 4;

  static constexpr int decoderTensorDims_ = 2;

  unsigned int nInputs_;
  std::vector<uint> encoderShape_;
  std::vector<uint> decoderShape_;

  unsigned bitsPerADC_;
  unsigned bitsPerNorm_;
  unsigned bitsPerCALQ_;
  unsigned bitsPerInput_;
  int maxBitsPerOutput_;

  std::vector<int> outputBitsPerLink_;

  std::vector<edm::ParameterSet> modelFilePaths_;

  std::string inputTensorName_encoder_;
  std::string outputTensorName_encoder_;
  std::string inputCondTensorName_encoder_;

  std::unique_ptr<tensorflow::GraphDef> graphDef_encoder_;
  std::vector<std::unique_ptr<tensorflow::Session>> session_encoder_;

  std::string inputTensorName_decoder_;
  std::string outputTensorName_decoder_;
  std::unique_ptr<tensorflow::GraphDef> graphDef_decoder_;
  std::vector<std::unique_ptr<tensorflow::Session>> session_decoder_;

  std::vector<unsigned int> linkToGraphMap_;

  double zeroSuppresionThreshold_;

  bool useModuleFactor_;
  bool bitShiftNormalization_;
  bool useTransverseADC_;
  bool normByMax_;

  bool skipAE_;

  bool saveEncodedValues_;
  bool preserveModuleSum_;

  std::array<double, nEncodedLayerNodes_> ae_encodedLayer_;

  HGCalTriggerTools triggerTools_;
  AEinputUtil aeInputUtil_;

  int verbose_;
  
};

#endif
