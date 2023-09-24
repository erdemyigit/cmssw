#include "L1Trigger/L1THGCal/interface/concentrator/HGCalConcentratorAutoEncoderImpl.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include <iomanip>


// Following example of implementing graphloading from here:
// https://gitlab.cern.ch/mrieger/CMSSW-TensorFlowExamples/-/blob/master/GraphLoading/

HGCalConcentratorAutoEncoderImpl::HGCalConcentratorAutoEncoderImpl(const edm::ParameterSet& conf)
    : cellRemap_(conf.getParameter<std::vector<int>>("cellRemap")),
      cellRemapNoDuplicates_(conf.getParameter<std::vector<int>>("cellRemapNoDuplicates")),
      encoderShape_(conf.getParameter<std::vector<uint>>("encoderShape")),
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
      aeInputUtil_(bitsPerADC_, bitsPerNorm_, bitsPerCALQ_, bitsPerInput_, useModuleFactor_, bitShiftNormalization_, useTransverseADC_, normByMax_){
  // find total size of the expected input shape
  // used for checking the maximum size used in cell Remap
  nInputs_ = 1;
  for (const auto& i : encoderShape_) {
    nInputs_ *= i;
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

  if (cellRemap_.size() != nInputs_) {
    throw cms::Exception("BadInitialization")
        << "Size of cellRemap (" << cellRemap_.size()
        << ") does not agree with the total size specified for the encoder inputs based on the encoderShape variable ("
        << nInputs_ << ")";
  }

  if (cellRemap_.size() != cellRemapNoDuplicates_.size()) {
    throw cms::Exception("BadInitialization")
        << "Size of cellRemap (" << cellRemap_.size() << ") does not agree with size of cellRemapNoDuplicates ("
        << cellRemapNoDuplicates_.size() << ")";
  }

  for (unsigned i = 0; i < cellRemap_.size(); i++) {
    if (cellRemap_[i] > nTriggerCells_ - 1) {
      throw cms::Exception("BadInitialization")
          << "cellRemap value " << cellRemap_[i] << " is larger than the number of trigger cells " << nTriggerCells_;
    }
    if (cellRemapNoDuplicates_[i] > nTriggerCells_ - 1) {
      throw cms::Exception("BadInitialization") << "cellRemapNoDuplicates value " << cellRemapNoDuplicates_[i]
                                                << " is larger than the number of trigger cells " << nTriggerCells_;
    }
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

void HGCalConcentratorAutoEncoderImpl::select(unsigned nLinks,
                                              const std::vector<l1t::HGCalTriggerCell>& trigCellVecInput,
                                              std::vector<l1t::HGCalTriggerCell>& trigCellVecOutput,
                                              std::vector<l1t::HGCalConcentratorData>& ae_encodedLayer_Output) {
  if(trigCellVecInput.empty()){
      return;
  }
  if(triggerTools_.isScintillator(trigCellVecInput[0].detId())){
      return;
  }
  std::cout << "1" << std::endl;
    fflush(stdout);
  std::array<double, nTriggerCells_> uncompressedCharge;
  std::array<double, nTriggerCells_> compressedCharge;

  std::array<double, maxAEInputSize_> ae_inputArray;
  std::array<double, nTriggerCells_> ae_outputArray;

  //reset inputs to 0 to account for zero suppressed trigger cells
  uncompressedCharge.fill(0);
  compressedCharge.fill(0);
  ae_inputArray.fill(0);
  ae_outputArray.fill(0);

  double modSum = 0;

  int bitsPerOutput = outputBitsPerLink_.at(nLinks);
  int nIntegerBits = 1;
  int nDecimalBits = bitsPerOutput - nIntegerBits;
  double outputSaturationValue = (1 << nIntegerBits) - 1./(1 << nDecimalBits);

  // largest expected input and output values, used for bit truncation
  // values of -1 for the number of bits used to keep full precision, in which case the MaxIntSize variables are not used
  double outputMaxIntSize = 1;
  if (bitsPerOutput > 0)
    outputMaxIntSize = 1 << nDecimalBits;
  double outputMaxIntSizeGlobal = 1;
  if (maxBitsPerOutput_ > 0)
    outputMaxIntSizeGlobal = 1 << (maxBitsPerOutput_ - nIntegerBits);

  aeInputUtil_.run(trigCellVecInput);
  if(aeInputUtil_.getModSum() <=0){
      return;
  }

  //legacy loop for uncompressed and compressed charges (don't really know what these are....)
  //unsigned aeInput2d[8][8] = {};
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

    ae_inputArray[inputIndex] = aeInputUtil_.getInput(inputIndex)/aeInputUtil_.getInputNorm();
    //aeInput2d[cellu][cellv] = aeInputUtil_.getInput(inputIndex);
  }
  modSum = aeInputUtil_.getModSum();

  double originalADCsum = 0;
  float originalCALQsum = 0;
  double originalINPUTsum = 0;

  for(unsigned u=0; u<8; ++u){
    for(unsigned v=0; v<8; ++v){
      originalADCsum += aeInputUtil_.getADC(u,v);
      originalCALQsum += aeInputUtil_.getCALQ(u,v);
      originalINPUTsum += aeInputUtil_.getInput(u,v)/aeInputUtil_.getInputNorm();
    }
  }
  std::cout << "2" << std::endl;
    fflush(stdout);
  tensorflow::Tensor encoder_input(tensorflow::DT_FLOAT,
                                   {encoderShape_[0], encoderShape_[1], encoderShape_[2], encoderShape_[3]});
  
  int n = encoderShape_[0]; // Number of 8x8 matrices
  // 1 x 4 x 4 x 3

  
  int numRows = encoderShape_[1]; // Number of rows in each 8x8 matrix
  int numCols = encoderShape_[2]; // Number of columns in each 8x8 matrix

    // Initialize a vector to store the sum results

    // Loop over each 8x8 matrix and calculate the sum using std::accumulate
  for (int i = 0; i < n; i++) {
      float* d = encoder_input.flat<float>().data() + i * numRows * numCols;
  }
  
  float* d = encoder_input.flat<float>().data();
  for (uint i = 0; i < nInputs_; i++, d++) {
    *d = ae_inputArray[i];
  }
  // Building conditional tensors
  tensorflow::Tensor inputTensor2(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, 2}));
    float* tensorData2 = inputTensor2.flat<float>().data();
    tensorData2[0] = static_cast<float>(originalCALQsum);
    tensorData2[1] = static_cast<float>(trigCellVecInput[0].eta());
  int num_elements = encoder_input.flat<float>().size();

  int graphIndex = linkToGraphMap_.at(nLinks);
//   for (int i = 0; i < num_elements; i++) {
//         printf("%f\n", d[i]);
//   }
   std::cout << "3" << std::endl;

//   std::cout << inputTensorName_encoder_ << std::endl;
//   std::cout << inputCondTensorName_encoder_ << std::endl;
  
  
//   inputMap[inputTensorName_encoder_] = encoder_input;
//   inputMap[inputCondTensorName_encoder_] = inputTensor2;  
  
  
//   tensorflow::NamedTensorList inputs = {
//     {inputTensorName_encoder_, encoder_input},
//     {inputCondTensorName_encoder_, inputTensor2} // Assuming "cond" is the name of your second input tensor
//   };
//   std::cout << "6" << std::endl;
//   fflush(stdout);
//   std::cout << inputTensorName_encoder_ << std::endl;
//   fflush(stdout);
//   std::cout << inputCondTensorName_encoder_ << std::endl;
//   fflush(stdout);
    
    
 // Assign to the second element

  
  std::vector<tensorflow::Tensor> encoder_outputs;
// std::map<std::string, tensorflow::Tensor> inputMap;
// std::cout << encoder_input << std::endl;
// std::cout << inputTensor2 << std::endl;
// inputMap[inputTensorName_encoder_] = encoder_input;
// inputMap[inputCondTensorName_encoder_] = inputTensor2;

// Construct input_args dynamically
// std::vector<std::pair<std::string, tensorflow::Tensor>> input_args;
// for (const auto& entry : inputMap) {
//     input_args.push_back(entry);
// }

// tensorflow::run(
//     session_encoder_.at(graphIndex).get(),
//     input_args,
//     {outputTensorName_encoder_},
//     &encoder_outputs
// );
// Maybe need to try putting the map as a class object?
// std::vector<std::pair<std::string, tensorflow::Tensor>> inputVector;
// inputVector.emplace_back(inputTensorName_encoder_, encoder_input);
// inputVector.emplace_back(inputCondTensorName_encoder_, inputTensor2);
//    vector<pair<string, Tensor>> feed_dict = {
//     {inputTensorName_encoder_, encoder_input},
//     {inputCondTensorName_encoder_, inputTensor2}
// };
// Construct input_args dynamically
// std::vector<std::pair<std::string, tensorflow::Tensor>> input_args = inputVector;

// tensorflow::run(
//     session_encoder_.at(graphIndex).get(),
//     feed_dict,
//     {inputTensorName_encoder_,inputCondTensorName_encoder_}
//     {outputTensorName_encoder_},
//     &encoder_outputs
// );
  tensorflow::run(session_encoder_.at(graphIndex).get(),
                  {{inputTensorName_encoder_,encoder_input},
                   {inputCondTensorName_encoder_,inputTensor2}},
                  {outputTensorName_encoder_},
                  &encoder_outputs);
  std::cout << "7" << std::endl;
  
//   std::vector<tensorflow::Tensor> encoder_outputs;
//   tensorflow::run(session_encoder_.at(graphIndex).get(),
//                   {{inputTensorName_encoder_, encoder_input}},
//                   {outputTensorName_encoder_},
//                   &encoder_outputs);

  if (encoder_outputs.empty()) {
    throw cms::Exception("BadInitialization") << "Autoencoder graph returning empty output vector";
  }
 
  //Here we can add conditioning to CMSSW
  
    
  std::cout << "8" << std::endl;
    fflush(stdout);
  d = encoder_outputs[0].flat<float>().data();
  for (int i = 0; i < encoder_outputs[0].NumElements(); i++, d++) {
    ae_encodedLayer_[i] = *d;
    //truncate the encoded layer bits
    if (bitsPerOutput > 0 && maxBitsPerOutput_ > 0) {
      ae_encodedLayer_[i] = std::min(std::floor(ae_encodedLayer_[i] * outputMaxIntSize) / outputMaxIntSize, outputSaturationValue);
    }
  }
    
  std::cout << "9" << std::endl;
    fflush(stdout);
  tensorflow::Tensor decoder_input(tensorflow::DT_FLOAT, {decoderShape_[0], decoderShape_[1]});
  d = decoder_input.flat<float>().data();
  for (int i = 0; i < nEncodedLayerNodes_; i++, d++) {
    *d = ae_encodedLayer_[i];
  }
  
  // Get the shape of the 'decoder_input' tensor
tensorflow::TensorShape decoder_input_shape = decoder_input.shape();

// Access and print the dimensions of the shape
// for (int i = 0; i < decoder_input_shape.dims(); i++) {
//     std::cout << "Dimension " << i << ": " << decoder_input_shape.dim_size(i) << std::endl;
// }
  std::cout << "10" << std::endl;
    fflush(stdout);
  std::vector<tensorflow::Tensor> decoder_outputs;
  tensorflow::run(session_decoder_.at(graphIndex).get(),
                  {{inputTensorName_decoder_, decoder_input}},
                  {outputTensorName_decoder_},
                  &decoder_outputs);


  d = decoder_outputs[0].flat<float>().data();
  for (uint i = 0; i < nInputs_; i++, d++) {
    if(skipAE_){
        ae_outputArray[i] = ae_inputArray[i];
    } else {
        ae_outputArray[i] = *d;
    }
  }
  std::cout << "11" << std::endl;
    fflush(stdout);
  // Add data back into trigger cells
  if (modSum > 0) {
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
    std::cout << "12" << std::endl;
      fflush(stdout);

    //use first TC to find mipPt conversions to Et and ADC
    float mipPtToEt_conv = trigCellVecInput[0].et() / trigCellVecInput[0].mipPt();
    float mipToADC_conv = trigCellVecInput[0].hwPt() / (trigCellVecInput[0].mipPt() * cosh(trigCellVecInput[0].eta()));

    double outputSum = 0;
    std::cout << "13" << std::endl;
      fflush(stdout);
    for (int i = 0; i < nTriggerCells_; i++) {
        cellU = aeInputUtil_.getU(i);
        cellV = aeInputUtil_.getV(i);
        
        HGCalTriggerDetId id(subdet, zp, type, layer, waferU, waferV, cellU, cellV);

        if(triggerTools_.getTriggerGeometry()->validTriggerCell(id)){
            outputSum += ae_outputArray[i];
        }
    }
    std::cout << "14" << std::endl;
      fflush(stdout);
    double renormalizationFactor = 1.;
    if (preserveModuleSum_ && outputSum > 0) {
      renormalizationFactor = modSum / outputSum;
    }
    std::cout << "15" << std::endl;
      fflush(stdout);
    //unsigned aeOutput2d[8][8] = {};
    //unsigned aeOutputCALQ2d[8][8] = {};
    //unsigned aeOutputADC2d[8][8] = {};
    double finalADCsum = 0;
    double finalCALQsum = 0;
    double finalINPUTsum = 0;
    //unsigned anID;
    std::cout << "16" << std::endl;
      fflush(stdout);
    for (int i = 0; i < nTriggerCells_; i++) {
      if (ae_outputArray[i] > 0) {
        cellU = aeInputUtil_.getU(i);
        cellV = aeInputUtil_.getV(i);

        HGCalTriggerDetId id(subdet, zp, type, layer, waferU, waferV, cellU, cellV);
        //anID = id;

        GlobalPoint point = triggerTools_.getTCPosition(id);

        double CALQ = ae_outputArray[i] * renormalizationFactor;
        double adc = aeInputUtil_.CALQtoADC(CALQ, cellU, cellV);

        double mipPt = adc / mipToADC_conv / cosh(point.eta());
        double et = mipPt * mipPtToEt_conv;

        finalINPUTsum += ae_outputArray[i];
        finalCALQsum += CALQ;
        finalADCsum += adc;

        //if (adc < zeroSuppresionThreshold_)
        //  continue;

        if (!triggerTools_.getTriggerGeometry()->validTriggerCell(id))
          continue;

        //aeOutput2d[cellU][cellV] = ae_outputArray[i] * aeInputUtil_.getInputNorm();

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
    std::cout << "17" << std::endl;
      fflush(stdout);
    //printf("%d %d %d %d\n", useModuleFactor_, bitShiftNormalization_, useTransverseADC_, normByMax_);
    //printf("%d\n", triggerTools_.getTriggerGeometry()->getModuleFromTriggerCell(anID));
    //printf("'input' sum: %0.3f -> %0.3f\n", originalINPUTsum, finalINPUTsum);
    //printf("'CALQ' sum: %0.3f -> %0.3f\n", originalCALQsum, finalCALQsum);
    //printf("'ADC' sum: %0.3f -> %0.3f\n", originalADCsum, finalADCsum);
    //printf("\n");
    //printf("ae output [raw]\n");
    //aeInputUtil_.print2d(aeOutput2d);
    //printf("AE output [renormalized]\n");
    //aeInputUtil_.print2d(aeOutputCALQ2d);
    //printf("AE output [ADC]\n");
    //aeInputUtil_.print2d(aeOutputADC2d);
    //printf("AE input\n");
    //aeInputUtil_.print2d(aeInput2d);
    //printf("AE output\n");
    //aeInputUtil_.print2d(aeOutput2d);
    
    if (saveEncodedValues_) {
      id = HGCalTriggerDetId(subdet, zp, type, layer, waferU, waferV, 0, 0);
      for (int i = 0; i < nEncodedLayerNodes_; i++) {
        l1t::HGCalConcentratorData encodedLayerData(ae_encodedLayer_[i] * outputMaxIntSizeGlobal, i, id);
        ae_encodedLayer_Output.push_back(encodedLayerData);
      }
    }
    std::cout << "18" << std::endl;
      fflush(stdout);
  }
}
