#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorAutoEncoderImpl.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iomanip>
#include <cmath> 
#include <sstream>

// Following example of implementing graphloading from here:
// https://gitlab.cern.ch/mrieger/CMSSW-TensorFlowExamples/-/blob/master/GraphLoading/

HGCalConcentratorAutoEncoderImpl::HGCalConcentratorAutoEncoderImpl(const edm::ParameterSet& conf)
    : encoderShape_(conf.getParameter<std::vector<uint>>("encoderShape")),
      decoderShape_(conf.getParameter<std::vector<uint>>("decoderShape")),
      bitsPerADC_(conf.getParameter<unsigned>("bitsPerADC")),
      bitsPerNorm_(conf.getParameter<unsigned>("bitsPerNorm")),
      bitsPerCALQ_(conf.getParameter<unsigned>("bitsPerCALQ")),
      bitsPerInput_(conf.getParameter<unsigned>("nBitsPerInput")),
      maxBitsPerOutput_(conf.getParameter<int>("maxBitsPerOutput")),
      outputBitsPerLink_(conf.getParameter<std::vector<int>>("bitsPerLink")),
      modelFilePaths_(conf.getParameter<std::vector<edm::ParameterSet>>("modelFiles")),
      linkToGraphMap_(conf.getParameter<std::vector<unsigned int>>("linkToGraphMap")),
      zeroSuppresionThreshold_(conf.getParameter<double>("zeroSuppresionThreshold")),
      useModuleFactor_(conf.getParameter<bool>("useModuleFactor")),
      bitShiftNormalization_(conf.getParameter<bool>("bitShiftNormalization")),
      useTransverseADC_(conf.getParameter<bool>("useTransverseADC")),
      normByMax_(conf.getParameter<bool>("normByMax")),
      skipAE_(conf.getParameter<bool>("skipAE")),
      saveEncodedValues_(conf.getParameter<bool>("saveEncodedValues")),
      preserveModuleSum_(conf.getParameter<bool>("preserveModuleSum")),
      aeInputUtil_(bitsPerADC_, bitsPerNorm_, bitsPerCALQ_, bitsPerInput_,
                  useModuleFactor_, bitShiftNormalization_, useTransverseADC_, 
                  normByMax_),
      verbose_(conf.getParameter<int>("verbose")),
      condInputDebug_(conf.getUntrackedParameter<std::string>("condInputDebug", "")),
      condInputDebugFull_(conf.getUntrackedParameter<bool>("condInputDebugFull", false)){
  // find total size of the expected input shape
  nInputs_ = 1;
  for (const auto& i : encoderShape_) {
//     printf("%d\n", i);
    nInputs_ *= i;
  }
  if (verbose_) {
    printf("nInputs %d\n", nInputs_);
  }
  
  // check the size of the inputs shapes
  if (encoderShape_.size() != encoderTensorDims_) {
    throw cms::Exception("BadInitialization")
        << "Encoder input shapes are currently expected to be " << encoderTensorDims_ << " values";
  }

  if (decoderShape_.size() != decoderTensorDims_) {
    throw cms::Exception("BadInitialization")
        << "Encoder input shapes are currently expected to be " << decoderTensorDims_ << " values long";
  }


  tensorflow::setLogging("0");

  for (const auto& modelFilePset : modelFilePaths_) {
    std::string encoderPath = modelFilePset.getParameter<edm::FileInPath>("encoderModelFile").fullPath();
    std::string decoderPath = modelFilePset.getParameter<edm::FileInPath>("decoderModelFile").fullPath();

    // optional condition preprocessing of this graph: both vectors or neither, one value per condition input
    const bool hasCondOffset = modelFilePset.existsAs<std::vector<double>>("condOffset");
    const bool hasCondScale = modelFilePset.existsAs<std::vector<double>>("condScale");
    if (hasCondOffset != hasCondScale) {
      throw cms::Exception("BadInitialization") << "condOffset and condScale must be given together (graph "
                                                << condOffset_.size() << ")";
    }
    std::vector<double> condOffset, condScale;
    if (hasCondOffset) {
      condOffset = modelFilePset.getParameter<std::vector<double>>("condOffset");
      condScale = modelFilePset.getParameter<std::vector<double>>("condScale");
      if (decoderShape_[1] <= nEncodedLayerNodes_) {
        throw cms::Exception("BadInitialization") << "condOffset/condScale given for a decoder without condition inputs";
      }
      const unsigned nCond = decoderShape_[1] - nEncodedLayerNodes_;
      if (condOffset.size() != nCond || condScale.size() != nCond) {
        throw cms::Exception("BadInitialization") << "condOffset/condScale need " << nCond << " values each, got "
                                                  << condOffset.size() << "/" << condScale.size();
      }
      for (unsigned i = 0; i < nCond; ++i) {
        if (!std::isfinite(condOffset[i]) || !std::isfinite(condScale[i]) || condScale[i] == 0.) {
          throw cms::Exception("BadInitialization") << "condOffset/condScale value " << i << " not usable";
        }
      }
    }
    if (!condInputDebug_.empty()) {
      std::ostringstream os;
      os << std::hexfloat << "CONDAFFINE " << condInputDebug_ << " graph=" << condOffset_.size() << " offset=";
      for (unsigned i = 0; i < condOffset.size(); ++i)
        os << (i ? "," : "") << condOffset[i];
      os << " scale=";
      for (unsigned i = 0; i < condScale.size(); ++i)
        os << (i ? "," : "") << condScale[i];
      edm::LogVerbatim("CONDIN") << os.str();
    }
    condOffset_.push_back(condOffset);
    condScale_.push_back(condScale);

    graphDef_encoder_ = std::unique_ptr<tensorflow::GraphDef>{tensorflow::loadGraphDef(encoderPath)};

    // create a new session and add the graphDef
    session_encoder_.push_back(
        std::unique_ptr<tensorflow::Session>{tensorflow::createSession(graphDef_encoder_.get())});

    graphDef_decoder_ = std::unique_ptr<tensorflow::GraphDef>{tensorflow::loadGraphDef(decoderPath)};

    // create a new session and add the graphDef
    session_decoder_.push_back(
        std::unique_ptr<tensorflow::Session>{tensorflow::createSession(graphDef_decoder_.get())});

    //extract encoder tenser names from first graph, check that rest of the names are consistent
    if (modelFilePset == modelFilePaths_.front()) {
      inputTensorName_encoder_ = graphDef_encoder_.get()->node(0).name();
      //Could be it (print this out)
      inputCondTensorName_encoder_ = graphDef_encoder_.get()->node(1).name();
      
      outputTensorName_encoder_ = graphDef_encoder_.get()->node(graphDef_encoder_.get()->node_size() - 1).name();
        
      inputTensorName_decoder_ = graphDef_decoder_.get()->node(0).name();
      outputTensorName_decoder_ = graphDef_decoder_.get()->node(graphDef_decoder_.get()->node_size() - 1).name();
    } else {
      if (inputTensorName_encoder_ != graphDef_encoder_.get()->node(0).name()) {
        throw cms::Exception("BadInitialization") << "provided list of encoder graphs have different input nodes";
      }
      if (outputTensorName_encoder_ != graphDef_encoder_.get()->node(graphDef_encoder_.get()->node_size() - 1).name()) {
        throw cms::Exception("BadInitialization") << "provided list of encoder graphs have different output nodes";
      }
      if (inputTensorName_decoder_ != graphDef_decoder_.get()->node(0).name()) {
        throw cms::Exception("BadInitialization") << "provided list of decoder graphs have different input nodes";
      }
      if (outputTensorName_decoder_ != graphDef_decoder_.get()->node(graphDef_decoder_.get()->node_size() - 1).name()) {
        throw cms::Exception("BadInitialization") << "provided list of decoder graphs have different output nodes";
      }
    }
  }

  // check that the appropriate number of links have been specified
  if (linkToGraphMap_.size() <= maxNumberOfLinks_) {
    throw cms::Exception("BadInitialization")
        << "Autoencoder graph number must be specified for all link allocation possibilities. Only "
        << linkToGraphMap_.size() << " values specified while " << maxNumberOfLinks_ << "links are possible";
  }

  // check that all graph indices specified exist in the model file lists
  for (const auto& graphNumber : linkToGraphMap_) {
    if (graphNumber >= modelFilePaths_.size()) {
      throw cms::Exception("BadInitialization")
          << "Autoencoder graph number  " << graphNumber << " is larger than the size of the provided list of graphs "
          << modelFilePaths_.size();
    }
  }
}

void HGCalConcentratorAutoEncoderImpl::select(
                    unsigned nLinks,
                    const std::vector<l1t::HGCalTriggerCell>& trigCellVecInput,
                    std::vector<l1t::HGCalTriggerCell>& trigCellVecOutput,
                    std::vector<l1t::HGCalConcentratorData>& ae_encodedLayer_Output) {
    if(verbose_){
      printf("\n----------------------------------------------------------\n");
    }
  if(trigCellVecInput.empty()){
      return;
  }
  if(triggerTools_.isScintillator(trigCellVecInput[0].detId())){
      return;
  }

  std::vector<double> uncompressedCharge(nInputs_, 0);
  std::vector<double> compressedCharge(nInputs_, 0);

  std::vector<double> ae_inputArray(nInputs_, 0);
  std::vector<double> ae_outputArray(nInputs_, 0);
  
  //Nate Added
  //std::vector<double> ae_condArray(5, 0);


 

  //reset inputs to 0 to account for zero suppressed trigger cells
  double modSum = 0;

  int bitsPerOutput = outputBitsPerLink_.at(nLinks);
  int nIntegerBits = 1;
  int nDecimalBits = bitsPerOutput - nIntegerBits;
  double outputSaturationValue = (1 << nIntegerBits) - 1./(1 << nDecimalBits);

  // largest expected input and output values, used for bit truncation
  // values of -1 for the number of bits used to keep full precision, 
  // in which case the MaxIntSize variables are not used
  // NB from Simon: I haven't touched this
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-but-set-variable"

  double outputMaxIntSize = 1;

  #pragma GCC diagnostic pop
   
  if (bitsPerOutput > 0)
    outputMaxIntSize = 1 << nDecimalBits;
  double outputMaxIntSizeGlobal = 1;
  if (maxBitsPerOutput_ > 0)
   outputMaxIntSizeGlobal = 1 << (maxBitsPerOutput_ - nIntegerBits);

  aeInputUtil_.run(trigCellVecInput);
  if(aeInputUtil_.getModSum() <=0){
      return;
  }

  
  
    
  
  for (const auto& trigCell : trigCellVecInput) {
    HGCalTriggerDetId id(trigCell.detId());
    uint cellu = id.triggerCellU();
    uint cellv = id.triggerCellV();
    
    
    int inputIndex = aeInputUtil_.getAEIndex(cellu, cellv);
    if (inputIndex < 0) {
      throw cms::Exception("BadInitialization")
          << "Invalid index provided for trigger cell u=" << cellu << " v=" << cellv;
    }
    uncompressedCharge[inputIndex] = trigCell.uncompressedCharge();
    compressedCharge[inputIndex] = trigCell.compressedCharge();

    ae_inputArray[inputIndex] = aeInputUtil_.getInput(inputIndex)/
        aeInputUtil_.getInputNorm();
    if(verbose_){
        printf("tc (%u, %u) has ADC %u\n", cellu, cellv, aeInputUtil_.getADC(inputIndex));
    }
  
  }
  
  modSum = aeInputUtil_.getModSum();

  double originalADCsum = 0;
  double originalCALQsum = 0;
  double originalINPUTsum = 0;

  for(unsigned u=0; u<8; ++u){
    for(unsigned v=0; v<8; ++v){
      originalADCsum += aeInputUtil_.getADC(u,v);
      // Skip CALQ at index (0,0) to match training preprocessing
      // (process_data.py uses range(1,64) for CALQ, skipping CALQ0)
      if (!(u == 0 && v == 0)) {
        originalCALQsum += aeInputUtil_.getCALQ(u,v);
      }
      originalINPUTsum += aeInputUtil_.getInput(u,v)/aeInputUtil_.getInputNorm();
    }
  }
//   ae_condArray[4] = originalCALQsum;
//   ae_condArray[0] = trigCellVecInput[0].eta();
  
  tensorflow::Tensor encoder_input(tensorflow::DT_FLOAT,
                                   {encoderShape_[0], encoderShape_[1],
                                   encoderShape_[2], encoderShape_[3]});
    
    
//   std::cout << "Value of inputTensorName_encoder_: " << inputTensorName_encoder_ << std::endl;
//   std::cout << "Value of inputCondTensorName_encoder_: " << inputCondTensorName_encoder_ << std::endl;
//   std::cout << "Value of outputTensorName_encoder_: " << outputTensorName_encoder_ << std::endl;
  if (verbose_) {
    printf("\n originalCALQsum \n");
    printf("%0.3f ", originalCALQsum);
    printf("\n originalADCsum \n");
    printf("%0.3f ", originalADCsum);
    printf("\n originalINPUTsum \n");
    printf("%0.3f ", originalINPUTsum);
    printf("\n modSum \n");
    printf("%0.3f ", modSum);
  }
      
//   printf("Cond Data\n");
    
//   HGCalTriggerDetId id(trigCellVecInput.at(0).detId());
//   int type = id.type();
//   int layer = id.layer();
//   int waferU = id.waferU();
//   int waferV = id.waferV();
  
//   // Try not filling with anything
//   ae_condArray[1] = waferU;
//   ae_condArray[2] = waferV;
//   ae_condArray[3] = type;
//   tensorflow::Tensor encoder_cond(tensorflow::DT_FLOAT,
//                                    {encoderShape_[0], encoderShape_[1],
//                                    encoderShape_[2], encoderShape_[3]});
  

  for (unsigned i = 0; i < nInputs_; ++i) {
    encoder_input.flat<float>().data()[i] = ae_inputArray[i];
  }
  if(verbose_){
      printf("INPUT\n");
      for (unsigned i = 0; i < nInputs_; ++i) {
        printf("%0.3f ", ae_inputArray[i]);
        if((i+1)%8==0){
            printf("\n");
        }
      }
  }
    if(verbose_){
      printf("CALQ INPUT\n");
      for(unsigned u=0; u<8; ++u){
        for(unsigned v=0; v<8; ++v){
          printf("%d ", aeInputUtil_.getCALQ(u,v));
        }
      printf("\n");
      }
    }
  

  
                                   
  if(!skipAE_){//run AE
      int graphIndex = linkToGraphMap_.at(nLinks);
      
      //Nate Modified
      std::vector<tensorflow::Tensor> encoder_outputs;
      

      

      tensorflow::run(session_encoder_.at(graphIndex).get(),
                      {{inputTensorName_encoder_, encoder_input}
                       },
                      {outputTensorName_encoder_},
                      &encoder_outputs);
      

      if (encoder_outputs.empty()) {
        throw cms::Exception("BadInitialization") << "Autoencoder graph returning empty output vector";
      }
//       if(verbose_){
//           printf("LATENT SPACE\n");
//           for (int i = 0; i < encoder_outputs[0].NumElements(); i++) {
//             printf("%f\n", encoder_outputs[0].flat<float>().data()[i]);
//           }
//       }
      if(verbose_){
          
          printf("bitsPerOutput\n");
          printf("%d\n",bitsPerOutput);
          printf("outputMaxIntSize\n");
          printf("%f\n",outputMaxIntSize);
          printf("outputSaturationValue\n");
          printf("%f\n",outputSaturationValue);
          printf("nDecimalBits\n");
          printf("%d\n",nDecimalBits);
          printf("nIntegerBits\n");
          printf("%d\n",nIntegerBits);
           

          
          }
      
      for (int i = 0; i < encoder_outputs[0].NumElements(); i++) {
        ae_encodedLayer_[i] = encoder_outputs[0].flat<float>().data()[i];
        //truncate the encoded layer bits
      
        if (bitsPerOutput > 0 && maxBitsPerOutput_ > 0) {
          ae_encodedLayer_[i] = std::min(std::floor(ae_encodedLayer_[i] * outputMaxIntSize) / outputMaxIntSize, outputSaturationValue);
        }
      }
      
//       for (int i = 0; i < encoder_outputs[0].NumElements(); i++) {
//         ae_encodedLayer_[i] = encoder_outputs[0].flat<float>().data()[i];
//         //truncate the encoded layer bits
//         if (bitsPerOutput > 0 && maxBitsPerOutput_ > 0) {
//           ae_encodedLayer_[i] = ae_encodedLayer_[i];
//         }
//       }

      tensorflow::Tensor decoder_input(tensorflow::DT_FLOAT, 
              {decoderShape_[0], decoderShape_[1]});
      std::fill_n(decoder_input.flat<float>().data(), decoderShape_[1], 1.0f);
      for (int i = 0; i < nEncodedLayerNodes_; i++) {
        decoder_input.flat<float>().data()[i] = ae_encodedLayer_[i];
      }
      std::string condDebugLine;  // CONDIN printout, emitted after the decoder ran (condInputDebug only)

if (decoderShape_[1] > 16){
          if (verbose_) {
            printf("\nConditional Info\n");
            fflush(stdout);
          }
          HGCalTriggerDetId id(trigCellVecInput.at(0).detId());
          
          // 1. Eta: use geometric center of the silicon module (matches wafer.eta in training)
          unsigned moduleId = triggerTools_.getTriggerGeometry()->getModuleFromTriggerCell(
              trigCellVecInput[0].detId());
          GlobalPoint modulePos = triggerTools_.getTriggerGeometry()->getModulePosition(moduleId);
          decoder_input.flat<float>().data()[16] = modulePos.eta() / 3.1;
          
          // V and U
          decoder_input.flat<float>().data()[17] = (double)id.waferV() / 12.0;
          decoder_input.flat<float>().data()[18] = (double)id.waferU() / 12.0;
          
          // 2. One-hot encode the wafertype
          int wafertype0 = (id.type() == 0) ? 1 : 0;
          int wafertype1 = (id.type() == 1) ? 1 : 0;
          int wafertype2 = (id.type() == 2) ? 1 : 0;

          decoder_input.flat<float>().data()[19] = wafertype0;
          decoder_input.flat<float>().data()[20] = wafertype1;
          decoder_input.flat<float>().data()[21] = wafertype2;
          
          // 3. log(1 + CALQ sum over all 64 cells): training (preprocessing v2) sums CALQ0..CALQ63, so CALQ(0,0),
          //    which originalCALQsum skips, is added back here. originalCALQsum itself feeds only the verbose printout.
          decoder_input.flat<float>().data()[22] = std::log(originalCALQsum + aeInputUtil_.getCALQ(0, 0) + 1.0);
          
          // 4. (layer - 1)/46 with the CE-E offset for CE-H (id.layer() restarts at 1 in CE-H), as the training
          //    ntuple's wafer layer (HGCalTriggerNtupleHGCWafers: layerWithOffset)
          decoder_input.flat<float>().data()[23] = ((double)triggerTools_.layerWithOffset(id) - 1.0) / 46.0;

          // 5. Optional per-graph condition preprocessing, the arithmetic of the training's cond_norm.apply: the stored
          //    float32 value -> double, (x - offset) / scale in double, -> float32. Without constants: unchanged.
          const unsigned nCond = decoderShape_[1] - nEncodedLayerNodes_;
          float* condIn = decoder_input.flat<float>().data() + nEncodedLayerNodes_;
          const std::vector<float> rawCond(condIn, condIn + nCond);
          const std::vector<double>& condOffset = condOffset_.at(graphIndex);
          if (!condOffset.empty()) {
            const std::vector<double>& condScale = condScale_.at(graphIndex);
            for (unsigned i = 0; i < nCond; ++i)
              condIn[i] = static_cast<float>((static_cast<double>(condIn[i]) - condOffset[i]) / condScale[i]);
          }

          if (verbose_) {
            printf("INPUT\n");
            for (unsigned i = 0; i < decoderShape_[1]; ++i) {
              printf("%0.3f \n", decoder_input.flat<float>().data()[i]);
            }
            printf("END OF Conditionals \n");
            fflush(stdout);
          }
          if (!condInputDebug_.empty()) {
            std::ostringstream os;
            os << std::setprecision(17) << "CONDIN " << condInputDebug_ << " module=" << moduleId
               << " layer=" << id.layer() << " layerWithOffset=" << triggerTools_.layerWithOffset(id)
               << " nLinks=" << nLinks << " modSum=" << modSum << " calq00=" << aeInputUtil_.getCALQ(0, 0)
               << " in=";
            for (unsigned i = 0; i < decoderShape_[1]; ++i)
              os << (i ? "," : "") << decoder_input.flat<float>().data()[i];
            os << " graph=" << graphIndex << " raw=";
            for (unsigned i = 0; i < nCond; ++i)
              os << (i ? "," : "") << rawCond[i];
            if (condInputDebugFull_) {
              os << " enc=";
              for (unsigned i = 0; i < nInputs_; ++i)
                os << (i ? "," : "") << encoder_input.flat<float>().data()[i];
            }
            condDebugLine = os.str();
          }
      }


      std::vector<tensorflow::Tensor> decoder_outputs;
      tensorflow::run(session_decoder_.at(graphIndex).get(),
                      {{inputTensorName_decoder_, decoder_input}},
                      {outputTensorName_decoder_},
                      &decoder_outputs);
      if (!condDebugLine.empty()) {
        if (condInputDebugFull_) {
          std::ostringstream os;
          os << std::setprecision(17) << " out=";
          for (unsigned i = 0; i < nInputs_; ++i)
            os << (i ? "," : "") << decoder_outputs[0].flat<float>().data()[i];
          condDebugLine += os.str();
        }
        edm::LogVerbatim("CONDIN") << condDebugLine;
      }
      if (verbose_) {
        printf("Past Decoder\n");
        fflush(stdout);
      }
      for (uint i = 0; i < nInputs_; i++) {
        ae_outputArray[i] = decoder_outputs[0].flat<float>().data()[i];
      }
  } else { //skipAE
      for (uint i = 0; i < nInputs_; i++) {
        ae_outputArray[i] = ae_inputArray[i];
      }
  }//endif skipAE

  if (verbose_){
      printf("\nOUTPUTS\n");
      for(unsigned i=0; i<nInputs_; ++i){
          printf("%0.3f ", ae_outputArray[i]);
          if((i+1)%8==0){
              printf("\n");
          }
      }
  }

  // Add data back into trigger cells
  if (modSum >= 0) {
    //get detID for everything but cell, take first entry detID and subtract off cellU and cellV contribution
    HGCalTriggerDetId id(trigCellVecInput.at(0).detId());
    int subdet = id.subdet();
    int zp = id.zside();
    int type = id.type();
    int layer = id.layer();
    int waferU = id.waferU();
    int waferV = id.waferV();
    int cellU = id.triggerCellU();
    int cellV = id.triggerCellV();



    //use first TC to find mipPt conversions to Et and ADC
    float mipPtToEt_conv = trigCellVecInput[0].et() / trigCellVecInput[0].mipPt();
    float mipToADC_conv = trigCellVecInput[0].hwPt() / (trigCellVecInput[0].mipPt() * cosh(trigCellVecInput[0].eta()));
    double outputSum = 0;
    fflush(stdout);  
     
        
    for (unsigned i = 0; i < nInputs_; i++) {
        fflush(stdout);
        cellU = aeInputUtil_.getUtc(i);
        cellV = aeInputUtil_.getVtc(i);
    
        fflush(stdout);
        if(cellU<0 || cellV<0){
            continue;
        }
        HGCalTriggerDetId id(subdet, zp, type, layer, waferU, waferV, cellU, cellV);
        if(triggerTools_.getTriggerGeometry()->validTriggerCell(id)){
            outputSum += ae_outputArray[i];
        }
    }
    double renormalizationFactor = 1.;
    if (preserveModuleSum_ && outputSum > 0) {
      renormalizationFactor = modSum / outputSum;
    }

    double finalADCsum = 0;
    double finalCALQsum = 0;
    double finalINPUTsum = 0;
    for (unsigned i = 0; i < nInputs_; i++) {
      if (ae_outputArray[i] > 0) {
        cellU = aeInputUtil_.getUtc(i);
        cellV = aeInputUtil_.getVtc(i);
        if(cellU<0 || cellV<0){
            continue;
        }
          
        HGCalTriggerDetId id(subdet, zp, type, layer, waferU, waferV, cellU, cellV);
        GlobalPoint point = triggerTools_.getTCPosition(id);
          
        double CALQ = ae_outputArray[i] * renormalizationFactor;
        double adc = aeInputUtil_.CALQtoADC(CALQ, i);

        double mipPt = adc / mipToADC_conv / cosh(point.eta());
        double et = mipPt * mipPtToEt_conv;
        finalINPUTsum += ae_outputArray[i];
        finalCALQsum += CALQ;
        finalADCsum += adc;
        if (adc < zeroSuppresionThreshold_)
          continue;

        if (!triggerTools_.getTriggerGeometry()->validTriggerCell(id))
          continue;

        if(verbose_){
            printf("tc (%u, %u) has ADC %f\n", cellU, cellV, adc);
        }
        l1t::HGCalTriggerCell triggerCell(reco::LeafCandidate::LorentzVector(), adc, 0, 0, 0, id);
          
        //Keep the pre-autoencoder charge for this cell
        triggerCell.setUncompressedCharge(uncompressedCharge[i]);
        triggerCell.setCompressedCharge(compressedCharge[i]);
        triggerCell.setMipPt(mipPt);

        math::PtEtaPhiMLorentzVector p4(et, point.eta(), point.phi(), 0.);

        triggerCell.setP4(p4);
        triggerCell.setPosition(point);

        trigCellVecOutput.push_back(triggerCell);
      }
    }
    

    if (saveEncodedValues_) {
      id = HGCalTriggerDetId(subdet, zp, type, layer, waferU, waferV, 0, 0);
      for (int i = 0; i < nEncodedLayerNodes_; i++) {
        l1t::HGCalConcentratorData encodedLayerData(ae_encodedLayer_[i] * outputMaxIntSizeGlobal, i, id);
        ae_encodedLayer_Output.push_back(encodedLayerData);
          }
    }
    if(verbose_){
      printf("------------------------------------------------------------\n");
    }
  }
}
