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
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:DoubleElectron_FlatPt-1To100-gun_noPU.root'),
    #fileNames = cms.untracked.vstring('/store/mc/Phase2Fall22DRMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/noPU_125X_mcRun4_realistic_v2-v1/2550000/066944a3-a061-42ac-ba45-9faadb46407a.root'),
    fileNames = cms.untracked.vstring('/store/mc/Phase2Fall22DRMiniAOD/DoubleElectron_FlatPt-1To100-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/65ce4640-c197-4c07-9fa4-cb505ab72738.root'),                       
    #fileNames = cms.untracked.vstring('file:/uscms/home/eertorer/nobackup/CMSSW_12_5_2_patch1/src/L1Trigger/L1THGCalUtilities/test/data_for_local_test/65ce4640-c197-4c07-9fa4-cb505ab72738.root'), #For local tests
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tTkPrimaryVertexs_L1TkPrimaryVertex__RECO',
    )
)

process.options = cms.untracked.PSet()
process.options.numberOfThreads = cms.untracked.uint32(4)
process.options.numberOfStreams = cms.untracked.uint32(4)
# process.options.numberOfThreads = cms.untracked.uint32(1)
# process.options.numberOfStreams = cms.untracked.uint32(1)

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



## Dev
# eLinks [1,2,3,4,5,6,7,8,9,10,11]

# eLinkCAE_1 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_1_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_1_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_2 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_2_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_2_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_3 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_3_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_3_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_4 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_4_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_4_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_5 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_5_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_5_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_6 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_6_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_6_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_7 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_7_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_7_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_8 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_8_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_8_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_9 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_9_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_9_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_10 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_10_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_10_eLinks/decoder_vanilla_AE.pb'))

# eLinkCAE_11 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_11_eLinks/encoder_vanilla_AE.pb'),
#                                   decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/CAE_6_20_CMSSW/model_11_eLinks/decoder_vanilla_AE.pb'))



# chains.register_concentrator("Newelink", concentrator.CreateAutoencoder(
#     useTransverseADC=True,
#     skipAE=False,
#     modelFiles = [eLinkCAE_1,eLinkCAE_2,eLinkCAE_3,eLinkCAE_4,eLinkCAE_5,eLinkCAE_6,eLinkCAE_7,
#                  eLinkCAE_8,eLinkCAE_9,eLinkCAE_10,eLinkCAE_11],
#     useModuleFactor=False,
#     bitShiftNormalization=True,
#     normByMax=False,
#     verbose =True, 
#     linkToGraphMap = cms.vuint32([0,1,2,3,4,5,6,7,8,9,10,10,10,10]),
#     encoderShape=cms.vuint32([1,8,8,1]),
#     decoderShape=cms.vuint32([1,24]),
# ))

#eLinkCAE_3 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/Tele_CAE_biased_90_CMSSW/model_3_bits/encoder_vanilla_AE.pb'),
#                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/Tele_CAE_biased_90_CMSSW/model_3_bits/decoder_vanilla_AE.pb'))
#
#eLinkCAE_5 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/Tele_CAE_biased_90_CMSSW/model_5_bits/encoder_vanilla_AE.pb'),
#                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/Tele_CAE_biased_90_CMSSW/model_5_bits/decoder_vanilla_AE.pb'))
#
#eLinkCAE_7 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/Tele_CAE_biased_90_CMSSW/model_7_bits/encoder_vanilla_AE.pb'),
#                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/Tele_CAE_biased_90_CMSSW/model_7_bits/decoder_vanilla_AE.pb'))
#
#eLinkCAE_9 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/Tele_CAE_biased_90_CMSSW/model_9_bits/encoder_vanilla_AE.pb'),
#                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/Tele_CAE_biased_90_CMSSW/model_9_bits/decoder_vanilla_AE.pb'))



#chains.register_concentrator("Papermodel", concentrator.CreateAutoencoder(
#    useTransverseADC=True,
#    skipAE=False,
#    modelFiles = [eLinkCAE_3,eLinkCAE_5,eLinkCAE_7,eLinkCAE_9],
#    useModuleFactor=False,
#    bitShiftNormalization=True,
#    normByMax=False,
#    verbose =True, 
#    linkToGraphMap = cms.vuint32([0,0,0,1,2,3,3,3,3,3,3,3,3,3]), #Edited this to match the number of bits
#    encoderShape=cms.vuint32([1,8,8,1]),
#    decoderShape=cms.vuint32([1,24]), #If you want to change the model you have to change the decoder shape but not the encoder shape (encoder shape is already written on the ASIC
#))


#

#L1Trigger/L1THGCal/data/models/elink_90_20_files_500_epoch/CMSSW_models/model_2_eLinks

#redo_tele_90_latest/CMSSW_models/model_2_eLinks/encoder_vanilla_AE.pb

# CAE Concentrators
#eLinkCAE_2 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/HyperbandNobias_May262025/best_model_eLink_2_post_seed_variation_larger_dataset_for_CMSSW/encoder_search.pb'),                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/HyperbandNobias_May262025/best_model_eLink_2_post_seed_variation_larger_dataset_for_CMSSW/decoder_search.pb')) 
#eLinkCAE_3 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/HyperbandNobias_May262025/best_model_eLink_3_post_seed_variation_larger_dataset_for_CMSSW/encoder_search.pb'),                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/HyperbandNobias_May262025/best_model_eLink_3_post_seed_variation_larger_dataset_for_CMSSW/decoder_search.pb'))
#eLinkCAE_4 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/HyperbandNobias_May262025/best_model_eLink_4_post_seed_variation_larger_dataset_for_CMSSW/encoder_search.pb'),                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/HyperbandNobias_May262025/best_model_eLink_4_post_seed_variation_larger_dataset_for_CMSSW/decoder_search.pb'))
#eLinkCAE_5 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/HyperbandNobias_May262025/best_model_eLink_5_post_seed_variation_larger_dataset_for_CMSSW/encoder_search.pb'),                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/HyperbandNobias_May262025/best_model_eLink_5_post_seed_variation_larger_dataset_for_CMSSW/decoder_search.pb'))

# AE Concentrators
#eLinkAE_2 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/AE_hyperband_8_15/best_model_eLink_2_post_seed_variation_larger_dataset_for_CMSSW/encoder_search.pb'),                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/AE_hyperband_8_15/best_model_eLink_2_post_seed_variation_larger_dataset_for_CMSSW/decoder_search.pb'))

#eLinkAE_3 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/AE_hyperband_8_15/best_model_eLink_3_post_seed_variation_larger_dataset_for_CMSSW/encoder_search.pb'),                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/AE_hyperband_8_15/best_model_eLink_3_post_seed_variation_larger_dataset_for_CMSSW/decoder_search.pb'))

#eLinkAE_4 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/AE_hyperband_8_15/best_model_eLink_4_post_seed_variation_larger_dataset_for_CMSSW/encoder_search.pb'),                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/AE_hyperband_8_15/best_model_eLink_4_post_seed_variation_larger_dataset_for_CMSSW/decoder_search.pb'))

#eLinkAE_5 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/AE_hyperband_8_15/best_model_eLink_5_post_seed_variation_larger_dataset_for_CMSSW/encoder_search.pb'),                                  decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/AE_hyperband_8_15/best_model_eLink_5_post_seed_variation_larger_dataset_for_CMSSW/decoder_search.pb'))

# =============================================================================
# July-2026 production block (v3-trained results_march30_2026 models) -- kept
# commented for the record.  It overrode linkToGraphMap without touching
# bitsPerLink; see ECON_V4_BRANCH_NOTES.md.  Superseded by the v4 block below.
# =============================================================================
# # CAE Concentrators
# eLinkCAE_2 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/CAE/eLink2/encoder_CAE_model.pb'), decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/CAE/eLink2/decoder_CAE_model.pb'))
# eLinkCAE_3 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/CAE/eLink3/encoder_CAE_model.pb'), decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/CAE/eLink3/decoder_CAE_model.pb'))
# eLinkCAE_4 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/CAE/eLink4/encoder_CAE_model.pb'), decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/CAE/eLink4/decoder_CAE_model.pb'))
# eLinkCAE_5 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/CAE/eLink5/encoder_CAE_model.pb'), decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/CAE/eLink5/decoder_CAE_model.pb'))
#
# # AE Concentrators
# eLinkAE_2 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/AE/eLink2/encoder_AE_model.pb'), decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/AE/eLink2/decoder_AE_model.pb'))
# eLinkAE_3 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/AE/eLink3/encoder_AE_model.pb'), decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/AE/eLink3/decoder_AE_model.pb'))
# eLinkAE_4 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/AE/eLink4/encoder_AE_model.pb'), decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/AE/eLink4/decoder_AE_model.pb'))
# eLinkAE_5 = cms.PSet(encoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/AE/eLink5/encoder_AE_model.pb'), decoderModelFile = cms.FileInPath('L1Trigger/L1THGCal/data/models/results_march30_2026/AE/eLink5/decoder_AE_model.pb'))
#
#
# # NOTE: Decoder dimensions for CAE (1,24) and AE (1,16) vary due to the concetanation of 8D Conditional Variables Vector 
# chains.register_concentrator("CAE", concentrator.CreateAutoencoder(
#     useTransverseADC=True,
#     skipAE=False,
#     modelFiles = [eLinkCAE_2,eLinkCAE_3,eLinkCAE_4,eLinkCAE_5],
#     useModuleFactor=False,
#     bitShiftNormalization=True,
#     normByMax=False,
#     verbose =True, 
#     linkToGraphMap = cms.vuint32([0,0,0,1,1,2,2,3,3,3,3,3,3,3]),
#     encoderShape=cms.vuint32([1,8,8,1]),
#     decoderShape=cms.vuint32([1,24]),
# ))
#
#
# chains.register_concentrator("AE", concentrator.CreateAutoencoder(
#     useTransverseADC=True,
#     skipAE=False,
#     modelFiles = [eLinkAE_2,eLinkAE_3,eLinkAE_4,eLinkAE_5],
#     useModuleFactor=False,
#     bitShiftNormalization=True,
#     normByMax=False,
#     verbose =True, 
#     linkToGraphMap = cms.vuint32([0,0,0,1,1,2,2,3,3,3,3,3,3,3]),
#     encoderShape=cms.vuint32([1,8,8,1]),
#     decoderShape=cms.vuint32([1,16]),
# ))

# =============================================================================
# ECON v4 models (branch econ-v4-production)
# Full derivation: L1Trigger/L1THGCal/test/ECON_V4_BRANCH_NOTES.md
#
# (i) ROUTING nLinks -> graph.  NO linkToGraphMap override: the upstream default
#     from L1Trigger/L1THGCal/python/l1tHGCalConcentratorProducer_cfi.py applies:
#         linkToGraphMapping = [0,0,0,1,2,3,3,3,3,3,3,3,3,3,3]
#         index (= nLinks):     0 1 2 3 4 5 6 7 8 9 . . . . 14
#     nLinks 2 -> graph 0 (eLink-2 model), 3 -> graph 1 (eLink-3 model),
#     4 -> graph 2 (eLink-4 model), >=5 -> graph 3 (eLink-5 model).
#     The v4 training classes (class_from_nlinks / Table 1a) are exactly this
#     table, so the July-2026 override [0,0,0,1,1,2,2,3,...] must NOT be used.
#
# (ii) bitsPerLink is indexed by nLinks INDEPENDENTLY of linkToGraphMap
#     (HGCalConcentratorAutoEncoderImpl.cc: bitsPerOutput = bitsPerLink.at(nLinks);
#     graphIndex = linkToGraphMap.at(nLinks)).  Each graph must therefore be fed
#     the latent depth it was trained at: eLink-2/3/4/5 models -> 3/5/7/9 bits.
#                          nLinks:  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
econ_v4_bitsPerLink = cms.vint32([0, 1, 3, 5, 7, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9])
#                    graph index:  0  0  0  1  2  3  3  3  3  3  3  3  3  3  3
#            graph's trained bits: 3  3  3  5  7  9  9  9  9  9  9  9  9  9  9
#     Indices 2..14 match their graph.  Index 0 is unreachable (nLinks = 0 never
#     occurs in the CMSSW link table).  This vector equals the upstream default
#     autoEncoder_bitsPerOutputLink; it is written out so the pairing is explicit.
#
# (iii) ######## DECISION REQUIRED: nLinks == 1 (index 1) ########
#     Under CMSSW's own link table (hgcal_trigger_link_mapping_120links_v1.json)
#     485,731 wafers = 2.51 % have nLinks = 1, spread over planes 7-47.
#     linkToGraphMap[1] = 0 routes them to the eLink-2 model (trained at 3 bits)
#     while bitsPerLink[1] = 1 truncates their latent to 1 bit.
#       Option A (CURRENT, upstream value): bitsPerLink[1] = 1.  Honest 1-eLink
#                bandwidth, but the 3-bit graph decodes a 1-bit latent it never saw.
#       Option B: bitsPerLink[1] = 3.  The graph sees its trained depth; the module
#                is credited 2 bits per latent value more than one eLink carries.
#       (A 5th, 1-bit-trained graph would be Option C; out of scope here.)
#     NOT decided on this branch -- index 1 is deliberately left at 1.
#     ################################################################
#
# (iv) decoderShape: AE decoder input is [1,16]; CAE and FiLM decoders take
#     [1,24] (16 latent ++ 8 conditions; the FiLM graph is exported single-input
#     and slices internally).  encoderShape stays [1,8,8,1] (fixed by the ASIC).
#
# Model files: training/final_run/export_cmssw.py in the ECON_CAE repo writes
#     <out>/<ARM>/eLink<N>/encoder_<ARM>_model.pb and decoder_<ARM>_model.pb
# for ARM in {AE, CAE, FILM}, N in {2,3,4,5}.  Copy them under
# L1Trigger/L1THGCal/data/models/<tag>/ and replace the PLACEHOLDER root below.
# =============================================================================
ECON_V4_MODEL_ROOT = '/PATH/TO/econ_v4_models'   # PLACEHOLDER -- e.g. 'L1Trigger/L1THGCal/data/models/econ_v4'

def _econ_v4_model(arm, nLinks):
    d = '%s/%s/eLink%d' % (ECON_V4_MODEL_ROOT, arm, nLinks)
    return cms.PSet(encoderModelFile = cms.FileInPath('%s/encoder_%s_model.pb' % (d, arm)),
                    decoderModelFile = cms.FileInPath('%s/decoder_%s_model.pb' % (d, arm)))

# ---- AE arm (plain autoencoder, decoderShape [1,16]) ------------------------
eLinkAE_2, eLinkAE_3, eLinkAE_4, eLinkAE_5 = [_econ_v4_model('AE', n) for n in (2, 3, 4, 5)]
chains.register_concentrator("AE", concentrator.CreateAutoencoder(
    useTransverseADC=True,
    skipAE=False,
    modelFiles = [eLinkAE_2, eLinkAE_3, eLinkAE_4, eLinkAE_5],   # graph 0,1,2,3
    useModuleFactor=False,
    bitShiftNormalization=True,
    normByMax=False,
    verbose =True,
    bitsPerLink = econ_v4_bitsPerLink,
    # linkToGraphMap: upstream default, see (i) -- do not override
    encoderShape=cms.vuint32([1,8,8,1]),
    decoderShape=cms.vuint32([1,16]),
    preserveModuleSum=True,
))

# ---- CAE arm (concat-conditioned, decoderShape [1,24]) ----------------------
eLinkCAE_2, eLinkCAE_3, eLinkCAE_4, eLinkCAE_5 = [_econ_v4_model('CAE', n) for n in (2, 3, 4, 5)]
chains.register_concentrator("CAE", concentrator.CreateAutoencoder(
    useTransverseADC=True,
    skipAE=False,
    modelFiles = [eLinkCAE_2, eLinkCAE_3, eLinkCAE_4, eLinkCAE_5],   # graph 0,1,2,3
    useModuleFactor=False,
    bitShiftNormalization=True,
    normByMax=False,
    verbose =True,
    bitsPerLink = econ_v4_bitsPerLink,
    # linkToGraphMap: upstream default, see (i) -- do not override
    encoderShape=cms.vuint32([1,8,8,1]),
    decoderShape=cms.vuint32([1,24]),
    preserveModuleSum=True,
))

# ---- FiLM arm (film-conditioned decoder, single-input graph, decoderShape [1,24])
# Registered but not in standard_concentrators by default; append 'FILM' to run it.
eLinkFILM_2, eLinkFILM_3, eLinkFILM_4, eLinkFILM_5 = [_econ_v4_model('FILM', n) for n in (2, 3, 4, 5)]
chains.register_concentrator("FILM", concentrator.CreateAutoencoder(
    useTransverseADC=True,
    skipAE=False,
    modelFiles = [eLinkFILM_2, eLinkFILM_3, eLinkFILM_4, eLinkFILM_5],   # graph 0,1,2,3
    useModuleFactor=False,
    bitShiftNormalization=True,
    normByMax=False,
    verbose =True,
    bitsPerLink = econ_v4_bitsPerLink,
    # linkToGraphMap: upstream default, see (i) -- do not override
    encoderShape=cms.vuint32([1,8,8,1]),
    decoderShape=cms.vuint32([1,24]),
    preserveModuleSum=True,
))


## BE1
chains.register_backend1("Dummy", clustering2d.CreateDummy())
## BE2
chains.register_backend2("Histomax", clustering3d.CreateHistoMax())
# Register selector
chains.register_selector("Dummy", selectors.CreateDummy())


# Register trigger chains
# standard_concentrators = ['Threshold0', 'Threshold135', 'Bestchoice', 'Supertriggercell', 'NateAE']
# standard_concentrators = ['Threshold0', 'Threshold135', 'adam500','SimonAdam500']SimonAdam500_Testing

# standard_concentrators = ['eLinkTele']

standard_concentrators = ['CAE','AE','Threshold0']


# standard_concentrators = ['eLinkemd']


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

