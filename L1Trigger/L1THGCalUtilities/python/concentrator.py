
import FWCore.ParameterSet.Config as cms
import SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi as digiparam
from L1Trigger.L1THGCal.l1tHGCalConcentratorProducer_cfi import threshold_conc_proc, best_conc_proc, supertc_conc_proc, coarsetc_onebitfraction_proc, custom_conc_proc, autoEncoder_conc_proc 


class CreateSuperTriggerCell(object):
    def __init__(self,
            stcSize=supertc_conc_proc.stcSize,
            type_energy_division=supertc_conc_proc.type_energy_division,
            fixedDataSizePerHGCROC=supertc_conc_proc.fixedDataSizePerHGCROC,
            coarsenTriggerCells=supertc_conc_proc.coarsenTriggerCells,
            ctcSize=supertc_conc_proc.ctcSize,
            ):
        self.processor = supertc_conc_proc.clone(
                stcSize = stcSize,
                type_energy_division = type_energy_division,
                fixedDataSizePerHGCROC = fixedDataSizePerHGCROC,
                coarsenTriggerCells = coarsenTriggerCells,
                ctcSize = ctcSize,
                )

    def __call__(self, process, inputs):
        producer = process.l1tHGCalConcentratorProducer.clone(
                InputTriggerCells = cms.InputTag(inputs),
                InputTriggerSums = cms.InputTag(inputs),
                ProcessorParameters = self.processor
                )
        return producer


class CreateThreshold(object):
    def __init__(self,
            threshold_silicon=threshold_conc_proc.threshold_silicon,  # in mipT
            threshold_scintillator=threshold_conc_proc.threshold_scintillator  # in mipT
            ):
        self.processor = threshold_conc_proc.clone(
                threshold_silicon = threshold_silicon,  # MipT
                threshold_scintillator = threshold_scintillator  # MipT
                )

    def __call__(self, process, inputs):
        producer = process.l1tHGCalConcentratorProducer.clone(
                InputTriggerCells = cms.InputTag(inputs),
                InputTriggerSums = cms.InputTag(inputs),
                ProcessorParameters = self.processor
                )
        return producer


class CreateBestChoice(object):
    def __init__(self,
            triggercells=best_conc_proc.NData,
            coarsenTriggerCells=best_conc_proc.coarsenTriggerCells,
            ctcSize=best_conc_proc.ctcSize,
            ):
        self.processor = best_conc_proc.clone(
                NData = triggercells,
                coarsenTriggerCells = coarsenTriggerCells,
                ctcSize=ctcSize,
                )

    def __call__(self, process, inputs):
        producer = process.l1tHGCalConcentratorProducer.clone(
                InputTriggerCells = cms.InputTag(inputs),
                InputTriggerSums = cms.InputTag(inputs),
                ProcessorParameters = self.processor
                )
        return producer

from L1Trigger.L1THGCal.ECONTritonProducer_cfi import AEProducer
from L1Trigger.L1THGCal.MapMergeProducer_cfi import MapMergeProducer

class CreateAutoencoder(object):
    def __init__(self,
            nBitsPerInput = autoEncoder_conc_proc.nBitsPerInput,
            maxBitsPerOutput = autoEncoder_conc_proc.maxBitsPerOutput,
            bitsPerLink = autoEncoder_conc_proc.bitsPerLink,
            modelFiles = autoEncoder_conc_proc.modelFiles,
            linkToGraphMap = autoEncoder_conc_proc.linkToGraphMap,
            encoderShape = autoEncoder_conc_proc.encoderShape,
            decoderShape = autoEncoder_conc_proc.decoderShape,
            zeroSuppresionThreshold = autoEncoder_conc_proc.zeroSuppresionThreshold,
            saveEncodedValues = autoEncoder_conc_proc.saveEncodedValues,
            preserveModuleSum = autoEncoder_conc_proc.preserveModuleSum,
            threshold_scintillator = autoEncoder_conc_proc.threshold_scintillator,
            threshold_silicon = autoEncoder_conc_proc.threshold_silicon,
            useTransverseADC = autoEncoder_conc_proc.useTransverseADC,
            clipTransverseADC = autoEncoder_conc_proc.clipTransverseADC,
            skipAE = autoEncoder_conc_proc.skipAE,
            bitsPerADC = autoEncoder_conc_proc.bitsPerADC,
            bitsPerNorm = autoEncoder_conc_proc.bitsPerNorm,
            bitsPerCALQ = autoEncoder_conc_proc.bitsPerCALQ,
            useModuleFactor = autoEncoder_conc_proc.useModuleFactor, 
            bitShiftNormalization = autoEncoder_conc_proc.bitShiftNormalization,
            normByMax = autoEncoder_conc_proc.normByMax,
            scintillatorMethod = 'thresholdSelect',
            verbose = 0):
        self.processor = autoEncoder_conc_proc.clone(
            encoderShape = encoderShape,
            decoderShape = decoderShape,
            bitsPerADC = bitsPerADC,
            bitsPerNorm = bitsPerNorm,
            bitsPerCALQ = bitsPerCALQ,
            nBitsPerInput = nBitsPerInput,
            maxBitsPerOutput = maxBitsPerOutput,
            bitsPerLink = bitsPerLink,
            modelFiles = modelFiles,
            linkToGraphMap = linkToGraphMap,
            zeroSuppresionThreshold = zeroSuppresionThreshold,
            useModuleFactor = useModuleFactor,
            bitShiftNormalization = bitShiftNormalization,
            useTransverseADC = useTransverseADC,
            normByMax = normByMax,
            skipAE = skipAE,
            saveEncodedValues = saveEncodedValues,
            preserveModuleSum = preserveModuleSum,
            threshold_scintillator = threshold_scintillator,
            threshold_silicon = threshold_silicon,
            Method = cms.vstring(['autoEncoder','autoEncoder', scintillatorMethod]),
            verbose = cms.int32(verbose)
        )

    def __call__(self, process, inputs):
        producer = process.l1tHGCalConcentratorProducer.clone(
                InputTriggerCells = cms.InputTag(inputs),
                InputTriggerSums = cms.InputTag(inputs),
                ProcessorParameters = self.processor
        )
        return producer


class CreateOneBitFraction(object):
    def __init__(self,
            stcSize=coarsetc_onebitfraction_proc.stcSize,
            fixedDataSizePerHGCROC=coarsetc_onebitfraction_proc.fixedDataSizePerHGCROC
            ):
        self.processor = coarsetc_onebitfraction_proc.clone(
                stcSize = stcSize,
                fixedDataSizePerHGCROC = fixedDataSizePerHGCROC
                )

    def __call__(self, process, inputs):
        producer = process.l1tHGCalConcentratorProducer.clone(
                InputTriggerCells = cms.InputTag(inputs),
                InputTriggerSums = cms.InputTag(inputs),
                ProcessorParameters = self.processor
                )
        return producer


class CreateMixedFeOptions(object):
    def __init__(self,
            stcSize=custom_conc_proc.stcSize,
            type_energy_division=custom_conc_proc.type_energy_division,
            fixedDataSizePerHGCROC=custom_conc_proc.fixedDataSizePerHGCROC,
            triggercells=custom_conc_proc.NData
            ):
        self.processor = custom_conc_proc.clone(
                stcSize = stcSize,
                type_energy_division = type_energy_division,
                fixedDataSizePerHGCROC = fixedDataSizePerHGCROC,
                NData = triggercells,
                Method = cms.vstring('bestChoiceSelect','superTriggerCellSelect','superTriggerCellSelect'),
                )

    def __call__(self, process, inputs):
        producer = process.l1tHGCalConcentratorProducer.clone(
                InputTriggerCells = cms.InputTag(inputs),
                InputTriggerSums = cms.InputTag(inputs),
                ProcessorParameters = self.processor
                )
        return producer


class CreateCustom(object):
    def __init__(self,
            stcSize=custom_conc_proc.stcSize,
            type_energy_division=custom_conc_proc.type_energy_division,
            fixedDataSizePerHGCROC=custom_conc_proc.fixedDataSizePerHGCROC,
            triggercells=custom_conc_proc.NData,
            threshold_silicon=custom_conc_proc.threshold_silicon,  # in mipT
            threshold_scintillator=custom_conc_proc.threshold_scintillator,  # in mipT
            Method = custom_conc_proc.Method,
            coarsenTriggerCells=custom_conc_proc.coarsenTriggerCells,
            ctcSize=custom_conc_proc.ctcSize,
            ):
        self.processor = custom_conc_proc.clone(
                stcSize = stcSize,
                type_energy_division = type_energy_division,
                fixedDataSizePerHGCROC = fixedDataSizePerHGCROC,
                NData = triggercells,
                threshold_silicon = threshold_silicon,  # MipT
                threshold_scintillator = threshold_scintillator,  # MipT
                Method = Method,
                coarsenTriggerCells=coarsenTriggerCells,
                ctcSize = ctcSize,
                )

    def __call__(self, process, inputs):
        producer = process.l1tHGCalConcentratorProducer.clone(
                InputTriggerCells = cms.InputTag(inputs),
                InputTriggerSums = cms.InputTag(inputs),
                ProcessorParameters = self.processor
                )
        return producer

