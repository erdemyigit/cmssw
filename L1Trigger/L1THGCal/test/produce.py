import FWCore.ParameterSet.Config as cms 
from Configuration.ProcessModifiers.enableSonicTriton_cff import enableSonicTriton

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('DIGI',Phase2C9, enableSonicTriton)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3)
)

# Input source
process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring('file:/uscms/home/cmantill/nobackup/hgcal/CMSSW_11_0_0_patch1/src/EGM-Phase2HLTTDRWinter20DIGI-00009.root'),
       fileNames = cms.untracked.vstring('file:/home/submit/srothman/cmsdata/hgcal/Egun_200PU/MINIAOD/003B8BCB-93B0-4040-854A-04C77E4BD066.root'),
       inputCommands=cms.untracked.vstring(
           'keep *',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
           'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016s_simEmtfDigis__HLT',
           )
       )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleElectronPt10_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("ntuple_tfile.root")
)

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:nano.root'),
    outputCommands = cms.untracked.vstring(
      'drop *',
      'keep nanoaodFlatTable_*_*_*',
    )
)

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# generation


# path and endpath definitions
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
process.hgcl1tpg_step = cms.Path(process.L1THGCalTriggerPrimitivesTask)
process.endjob_step = cms.EndPath(process.endOfProcess)

# load HGCAL TPG simulation
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
process.load('L1Trigger.L1THGCalUtilities.HGC3DClusterGenMatchSelector_cff')
process.load('L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cff')
from L1Trigger.L1THGCalUtilities.hgcalTriggerChains import HGCalTriggerChains
import L1Trigger.L1THGCalUtilities.vfe as vfe
import L1Trigger.L1THGCalUtilities.concentrator as concentrator
import L1Trigger.L1THGCalUtilities.clustering2d as clustering2d
import L1Trigger.L1THGCalUtilities.clustering3d as clustering3d
import L1Trigger.L1THGCalUtilities.selectors as selectors
import L1Trigger.L1THGCalUtilities.customNtuples as ntuple

# fill cluster layer info
process.ntuple_multiclusters.FillLayerInfo = True
process.ntuple_triggercells.FillSimEnergy = True

chains = HGCalTriggerChains()
# Register algorithms
## VFE
chains.register_vfe("Floatingpoint", vfe.CreateVfe())
## ECON
#chains.register_concentrator("Supertriggercell", concentrator.CreateSuperTriggerCell())
#chains.register_concentrator("Threshold", concentrator.CreateThreshold())
chains.register_concentrator("Threshold0", concentrator.CreateThreshold(threshold_silicon=cms.double(-1.), threshold_scintillator=cms.double(-1.)))
#chains.register_concentrator("Bestchoice", concentrator.CreateBestChoice())

# AE models
triggerCellRemap = [28,29,30,31,0,4,8,12,
                    24,25,26,27,1,5,9,13,
                    20,21,22,23,2,6,10,14,
                    16,17,18,19,3,7,11,15,
                    47,43,39,35,-1,-1,-1,-1,
                    46,42,38,34,-1,-1,-1,-1,
                    45,41,37,33,-1,-1,-1,-1,
                    44,40,36,32,-1,-1,-1,-1]

AE_8x8_stride_telescope = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/QAEModels/encoder_lphe.pb'),
                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/QAEModels/decoder_lphe.pb'))
chains.register_concentrator("AE", 
                             concentrator.CreateAutoencoder(modelFiles = cms.VPSet([AE_8x8_stride_telescope]),
                                                            linkToGraphMap = cms.vuint32([0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
                                                            encoderShape=cms.vuint32([1,8,8,1]),
                                                            cellRemap = cms.vint32(triggerCellRemap),
                                                            saveEncodedValues = cms.bool(True),
                                                            cellRemapNoDuplicates = cms.vint32(triggerCellRemap),
                                                            #threshold_silicon = cms.double(-1.),
                                                            #threshold_scintillator = cms.double(-1.)
                                                            zeroSuppresionThreshold = cms.double(-1.)
                             ))


## BE1                                                                                                                                                                            
chains.register_backend1("Dummy", clustering2d.CreateDummy())
## BE2                                                                                                                                                                             
chains.register_backend2("Histomax", clustering3d.CreateHistoMax())
# Register selector                                                                                                                                                               
chains.register_selector("Genmatch", selectors.CreateGenMatch(distance=10000))

# Register ntuples
ntuple_list = ['event', 'gen', 'econdata', 'triggersums', 'triggercells', 'multiclusters']
chains.register_ntuple("Genclustersntuple", ntuple.CreateNtuple(ntuple_list))

# Register trigger chains
concentrator_algos = [
    #'Threshold','Threshold0','Bestchoice',
    'AE', 'Threshold0'
]
backend_algos = ['Histomax']

## Make cross product fo ECON and BE algos
import itertools
for cc,be in itertools.product(concentrator_algos,backend_algos):
    chains.register_chain('Floatingpoint', cc, 'Dummy', be, 'Genmatch', 'Genclustersntuple')

process = chains.create_sequences(process)

# Remove towers from sequence
#process.L1THGCalTriggerPrimitivesTask.remove(process.hgcalTowerMap)
#process.L1THGCalTriggerPrimitivesTask.remove(process.hgcalTower)

process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives)
process.selector_step = cms.Path(process.hgcalTriggerSelector)
process.ntuple_step = cms.Path(process.hgcalTriggerNtuples)

process.TCTable = cms.EDProducer("TCTableProducer",
    TriggerCells = cms.InputTag("FloatingpointAE", "HGCalConcentratorProcessorSelection"),
    #TriggerCells = cms.InputTag("Floatingpoint","HGCalVFEProcessorSums"),
    MultiClusters = cms.InputTag("FloatingpointAEDummyHistomaxGenmatch"),
    eeSimHits = cms.InputTag("g4SimHits","HGCHitsEE"),
    fhSimHits = cms.InputTag("g4SimHits","HGCHitsHEfront"),
    bhSimHits = cms.InputTag("g4SimHits","HGCHitsHEback"),
    keV2fC = cms.double(0.044259),
    fcPerMip = cms.vdouble(2.06, 3.43, 5.15),
    layerWeights = cms.vdouble(
        0, 8.894541, 0, 21.875814, 0,
        21.875814, 0, 21.875814, 0, 21.875814,
        0, 21.865764, 0, 21.875814, 0,
        21.876338, 0, 21.876338, 0, 21.876338,
        0, 21.876338, 0, 21.876338, 0,
        21.876338, 0, 21.876338, 0, 83.906398,
        51.444192, 51.444192, 51.444192, 51.444192, 51.444192,
        51.444192, 51.444192, 51.444192, 51.444192, 51.444192,
        69.513118, 87.582044, 87.582044, 87.582044, 87.582044,
        87.582044, 87.214571, 86.888309, 86.92952, 86.92952,
        86.92952
    ),
    thicknessCorrections = cms.vdouble(
        0.77, 0.77, 0.77, 0.84, 0.84,
        0.84
    )
)

process.table_step = cms.Path(process.TCTable) 
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

model='encode_decode'
process.AEProducer = cms.EDProducer("ECONTritonProducer",
    Client = cms.PSet(
        mode = cms.string("Async"),
        modelName = cms.string(model),
        modelConfigPath = cms.FileInPath("L1Trigger/L1THGCal/data/models/%s/config.pbtxt"%model),
        allowedTries = cms.untracked.uint32(1),
        timeout = cms.untracked.uint32(1)
    ),
    TriggerCells = cms.InputTag("Floatingpoint","HGCalVFEProcessorSums"),
)
process.AE_step = cms.Path(process.AEProducer)

# Schedule definition
process.schedule = cms.Schedule(process.hgcl1tpg_step, process.selector_step, process.ntuple_step, process.table_step, process.AE_step, process.NANOAODSIMoutput_step)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


