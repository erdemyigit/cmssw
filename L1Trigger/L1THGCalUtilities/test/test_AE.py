import FWCore.ParameterSet.Config as cms 

import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.ProcessModifiers.enableSonicTriton_cff import enableSonicTriton
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('DIGI',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:DoubleElectron_FlatPt-1To100-gun_noPU.root'),
    fileNames = cms.untracked.vstring('file:/home/submit/srothman/cmsdata/ECON_datasets/DoubleElectron_FlatPt-1To100_PU200/MINIAOD/008bdb86-d408-48f5-9ee7-ab93cda56cc1.root'),

    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tTkPrimaryVertexs_L1TkPrimaryVertex__RECO',
    )
)

process.options = cms.untracked.PSet()
process.options.numberOfThreads = cms.untracked.uint32(4)
process.options.numberOfStreams = cms.untracked.uint32(4)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleElectronPt10_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition
process.TFileService = cms.Service(
    "TFileService",
    #fileName = cms.string("/home/submit/srothman/cmsdata/hgcal/myntuples/ntuple_noAE6.root")
    #fileName = cms.string("ntuple.root")
    fileName = cms.string('ntuple.root')
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

# load HGCAL TPG simulation
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
process.load('L1Trigger.L1THGCalUtilities.HGC3DClusterSelectors_cff')
process.load('L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cff')
from L1Trigger.L1THGCalUtilities.hgcalTriggerChains import HGCalTriggerChains
import L1Trigger.L1THGCalUtilities.vfe as vfe
import L1Trigger.L1THGCalUtilities.concentrator as concentrator
import L1Trigger.L1THGCalUtilities.clustering2d as clustering2d
import L1Trigger.L1THGCalUtilities.clustering3d as clustering3d
import L1Trigger.L1THGCalUtilities.selectors as selectors
import L1Trigger.L1THGCalUtilities.customNtuples as ntuple
process.ntuple_triggercells.FillSimEnergy=True


chains = HGCalTriggerChains()
# Register algorithms
## VFE
chains.register_vfe("Floatingpoint", vfe.CreateVfe())
## ECON
ntuple_list = ['event', 'gen', 'multiclusters', 'triggercells']
chains.register_ntuple("nTuple", ntuple.CreateNtuple(ntuple_list))

chains.register_concentrator("Threshold0", concentrator.CreateThreshold(
  threshold_scintillator=cms.double(-1),
  threshold_silicon=cms.double(-1)
))
chains.register_concentrator("Threshold135", concentrator.CreateThreshold())
chains.register_concentrator("Bestchoice", concentrator.CreateBestChoice())
chains.register_concentrator("Supertriggercell", concentrator.CreateSuperTriggerCell())
chains.register_concentrator("Badae", concentrator.CreateAutoencoder(
    useTransverseADC=True,
    skipAE=False,
    useModuleFactor=False,
    bitShiftNormalization=True,
    normByMax=False,
))

#test models
autoEncoder_training_2eLinks = cms.PSet(encoderModelFile = cms.FileInPath('/work/submit/nswood/HGCAL/CMSSW_12_5_2_patch1/src/L1Trigger/L1THGCal/data/models/dummy_8_8/encoder_dummy_8_8.pb'),
                                        decoderModelFile = cms.FileInPath('/work/submit/nswood/HGCAL/CMSSW_12_5_2_patch1/src/L1Trigger/L1THGCal/data/models/dummy_8_8/decoder_dummy_8_8.pb'))

# autoEncoder_training_2eLinks = cms.PSet(encoderModelFile = cms.FileInPath('/work/submit/nswood/HGCAL/CMSSW_12_5_2_patch1/src/L1Trigger/L1THGCal/data/models/dummy_AE_4_4_3/encoder_dummy_4_4_3.pb'),
#                                         decoderModelFile = cms.FileInPath('/work/submit/nswood/HGCAL/CMSSW_12_5_2_patch1/src/L1Trigger/L1THGCal/data/models/dummy_AE_4_4_3/decoder_dummy_4_4_3.pb'))



linkToGraphMapping = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
remap_8x8 = [4, 12, 20, 28,  5, 13, 21, 29,  6, 14, 22, 30,  7, 15, 23, 31, 
             24, 25, 26, 27, 16, 17, 18, 19,  8,  9, 10, 11,  0,  1,  2,  3, 
             59, 51, 43, 35, 58, 50, 42, 34, 57, 49, 41, 33, 56, 48, 40, 32]
#https://github.com/cms-sw/cmssw/blob/master/L1Trigger/L1THGCal/python/l1tHGCalConcentratorProducer_cfi.py#L200

modelFiles = cms.VPSet([autoEncoder_training_2eLinks])

#modelFiles = cms.VPSet([autoEncoder_training_2eLinks, autoEncoder_training_2eLinks, autoEncoder_training_2eLinks, autoEncoder_training_2eLinks])



#Main problem is dealing with the conditioning in CMSSW
#Conditioning on sum is easy enough, but not sure how to condition on eta

#may have to reformat input as 58 inputs and reshape everything in keras model

#what inputs are given in CMSSW to the AE? all 58

#For inputs, just make sure it's the exact same as Rohan's and it should run

#https://github.com/ssrothman/cmssw/blob/ECON_12_5_2_patch1/L1Trigger/L1THGCal/src/concentrator/HGCalConcentratorAutoEncoderImpl.cc
triggerCellRemap = [28,29,30,31,0,4,8,12,
                    24,25,26,27,1,5,9,13,
                    20,21,22,23,2,6,10,14,
                    16,17,18,19,3,7,11,15,
                    47,43,39,35,-1,-1,-1,-1,
                    46,42,38,34,-1,-1,-1,-1,
                    45,41,37,33,-1,-1,-1,-1,
                    44,40,36,32,-1,-1,-1,-1]

chains.register_concentrator("NateAE", concentrator.CreateAutoencoder(
    useTransverseADC=True,
    skipAE=False,
    modelFiles = modelFiles,
    useModuleFactor=False,
    bitShiftNormalization=True,
    normByMax=False,
    linkToGraphMap = linkToGraphMapping,
    encoderShape=cms.vuint32([1,8,8,1]),
    cellRemap = cms.vint32(triggerCellRemap),
    cellRemapNoDuplicates = cms.vint32(triggerCellRemap)
))


# AE models                                                                                                                                                                                                                       

# AE_8x8_s2_tele_elec = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCalUtilities/test/AEmodels/8x8_c8_stride2_tele_ele_0PU_7_19/encoder_8x8_c8_S2_tele.pb'),
#                                decoderModelFile = cms.FileInPath('L1Trigger/L1THGCalUtilities/test/AEmodels/8x8_c8_stride2_tele_ele_0PU_7_19/decoder_8x8_c8_S2_tele.pb'))

# chains.register_concentrator("AutoEncoderStrideTelescopeEle",
#                              lambda p, i : concentrator.create_autoencoder(p, i,
#                            modelFiles = cms.VPSet([AE_8x8_s2_tele_elec]),
#                            linkToGraphMap = cms.vuint32([0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
#                            encoderShape=cms.vuint32([1,8,8,1]),
#                            cellRemap = cms.vint32(triggerCellRemap),
#                            cellRemapNoDuplicates = cms.vint32(triggerCellRemap)))


## BE1
chains.register_backend1("Dummy", clustering2d.CreateDummy())
## BE2
chains.register_backend2("Histomax", clustering3d.CreateHistoMax())
# Register selector
chains.register_selector("Dummy", selectors.CreateDummy())


# Register trigger chains
standard_concentrators = ['Threshold0', 'Threshold135', 'Bestchoice', 'Supertriggercell', 'NateAE']
for cc in standard_concentrators:
    chains.register_chain('Floatingpoint', cc, 'Dummy', 'Histomax', 'Dummy', 'nTuple')

process = chains.create_sequences(process)

# Remove towers from sequence
process.L1THGCalTriggerPrimitives.remove(process.L1THGCalTowerMap)
process.L1THGCalTriggerPrimitives.remove(process.L1THGCalTower)

from CommonTools.CandAlgos.genParticleCustomSelector_cfi import genParticleCustomSelector
process.filter = genParticleCustomSelector.clone(
    minRapidity = -1.444,
    maxRapidity = 1.444,
    invertRapidityCut = True,
    filter = cms.bool(True)
)

process.hgcl1tpg_step = cms.Path(process.L1THGCalTriggerPrimitives)
process.selector_step = cms.Path(process.L1THGCalTriggerSelector)
process.ntuple_step = cms.Path(process.L1THGCalTriggerNtuples)

process.fullpath = cms.Path(process.filter + process.L1THGCalTriggerPrimitives + process.L1THGCalTriggerSelector + process.L1THGCalTriggerNtuples)

# Schedule definition
process.schedule = cms.Schedule(process.fullpath)
#process.schedule = cms.Schedule(process.hgcl1tpg_step, process.selector_step, process.ntuple_step)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

