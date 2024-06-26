#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorAutoEncoderImpl.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include <iomanip>
#include <cmath> 

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
      verbose_(conf.getParameter<int>("verbose")){
  // find total size of the expected input shape
  nInputs_ = 1;
  for (const auto& i : encoderShape_) {
//     printf("%d\n", i);
    nInputs_ *= i;
  }
  printf("%d\n", nInputs_);
  
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
  std::vector<double> ae_condArray(5, 0);

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
      originalCALQsum += aeInputUtil_.getCALQ(u,v);
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
  printf("\n originalCALQsum \n");
  printf("%0.3f ", originalCALQsum);
  printf("\n originalADCsum \n");
  printf("%0.3f ", originalADCsum);
  printf("\n originalINPUTsum \n");  
  printf("%0.3f ", originalINPUTsum);
      
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
  

  if(verbose_){
      printf("INPUT\n");
      for (unsigned i = 0; i < nInputs_; ++i) {
        encoder_input.flat<float>().data()[i] = ae_inputArray[i];
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
      
      if (decoderShape_[1] > 16){
          printf("\nConditional Info\n");
          fflush(stdout);
          HGCalTriggerDetId id(trigCellVecInput.at(0).detId());
          decoder_input.flat<float>().data()[16] = trigCellVecInput[0].eta()/3.1;
          decoder_input.flat<float>().data()[17] = id.waferV()/12;
          decoder_input.flat<float>().data()[18] = id.waferU()/12;
          int wafertype0 = 0;
          int wafertype1 = 0;
          int wafertype2 = 0;
          if (id.type() == 0) {
                wafertype0 = 1;
            } else if (id.type() == 1) {
                wafertype1 = 1;
            } else if (id.type() == 2) {
                wafertype2 = 1;
            }
          

          decoder_input.flat<float>().data()[19] = wafertype0;
          decoder_input.flat<float>().data()[20] = wafertype1;
          decoder_input.flat<float>().data()[21] = wafertype2;
          decoder_input.flat<float>().data()[22] = log(originalCALQsum+1);
          decoder_input.flat<float>().data()[23] = (id.layer()-1)/(47-1);
          printf("INPUT\n");
          for (unsigned i = 0; i < decoderShape_[1]; ++i) {
            
            printf("%0.3f \n", decoder_input.flat<float>().data()[i]);
            
          }
          fflush(stdout);
          printf("END OF Conditionals \n");
      }
      fflush(stdout);
      
         
      std::vector<tensorflow::Tensor> decoder_outputs;
      tensorflow::run(session_decoder_.at(graphIndex).get(),
                      {{inputTensorName_decoder_, decoder_input}},
                      {outputTensorName_decoder_},
                      &decoder_outputs);
      printf("Past Decoder\n");
      fflush(stdout);
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
