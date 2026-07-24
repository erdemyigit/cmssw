import FWCore.ParameterSet.Config as cms

model = 'dummy'
AEProducer = cms.EDProducer("ECONTritonProducer",
    Client = cms.PSet(
        mode = cms.string("Async"),
        modelName = cms.string(model),
        modelConfigPath = cms.FileInPath("L1Trigger/L1THGCal/data/models/%s/config.pbtxt"%model),
        allowedTries = cms.untracked.uint32(1),
        timeout = cms.untracked.uint32(10),
        useSharedMemory=cms.untracked.bool(False)
    ),
    TriggerCells = cms.InputTag("Floatingpoint","HGCalVFEProcessorSums"),
    cut = cms.string(""),
    verbose = cms.uint32(0),
    bitsPerADC = cms.uint32(22),
    bitsPerNorm = cms.uint32(12),
    bitsPerCALQ = cms.uint32(23),
    nBitsPerInput = cms.uint32(8),
    bitShiftNormalization = cms.bool(True),
    useTransverseADC = cms.bool(True),
    useModuleFactor=cms.bool(False),
    normByMax=cms.bool(False),
)
