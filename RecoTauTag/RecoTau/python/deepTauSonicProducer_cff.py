import FWCore.ParameterSet.Config as cms

deepTauSonicProducer = cms.EDProducer("DeepTauIdSonicProducer",
    Client = cms.PSet(
        timeout = cms.untracked.uint32(300),
        mode = cms.string("Async"),
        modelName = cms.string("deeptau_ensemble"),
        modelConfigPath = cms.FileInPath("HeterogeneousCore/SonicTriton/data/models/deeptau_ensemble/config.pbtxt"),
        modelVersion = cms.string(""),
        verbose = cms.untracked.bool(False),
        allowedTries = cms.untracked.uint32(0),
        useSharedMemory = cms.untracked.bool(True),
        compression = cms.untracked.string(""),
    ),
    electrons = cms.InputTag('slimmedElectrons'),
    muons = cms.InputTag('slimmedMuons'),
    taus = cms.InputTag('slimmedTaus'),
    pfcands = cms.InputTag('packedPFCandidates'),
    vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
    rho = cms.InputTag('fixedGridRhoAll'),
    disable_dxy_pca = cms.bool(True)
)
