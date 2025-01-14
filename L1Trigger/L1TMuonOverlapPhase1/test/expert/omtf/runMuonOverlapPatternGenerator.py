# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys
import commands
import re
from os import listdir
from os.path import isfile, join

process.load("FWCore.MessageLogger.MessageLogger_cfi")

verbose = True

if verbose: 
    process.MessageLogger = cms.Service("MessageLogger",
       #suppressInfo       = cms.untracked.vstring('AfterSource', 'PostModule'),
       destinations   = cms.untracked.vstring(
                                               #'detailedInfo',
                                               #'critical',
                                               #'cout',
                                               #'cerr',
                                               'omtfEventPrint'
                    ),
       categories        = cms.untracked.vstring('l1tOmtfEventPrint', 'OMTFReconstruction'),
       omtfEventPrint = cms.untracked.PSet(    
                         filename  = cms.untracked.string('log_Patterns_0x00011_oldSample_3_30Files_test'),
                         extension = cms.untracked.string('.txt'),                
                         threshold = cms.untracked.string('DEBUG'),
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tOmtfEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) ),
                         OMTFReconstruction = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) )
                       ),
       debugModules = cms.untracked.vstring('simOmtfDigis') 
       #debugModules = cms.untracked.vstring('*')
    )

    #process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
if not verbose:
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(-1)
    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), 
                                         #SkipEvent = cms.untracked.vstring('ProductNotFound') 
                                     )
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D41_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

path = '/eos/user/k/kbunkow/cms_data/SingleMuFullEta/721_FullEta_v4/' #old sample, but very big
#path = '/eos/user/a/akalinow/Data/SingleMu/9_3_14_FullEta_v2/' #new sample, but small and more noisy
#path = '/eos/user/a/akalinow/Data/SingleMu/9_3_14_FullEta_v1/'

#path = '/afs/cern.ch/work/a/akalinow/public/MuCorrelator/Data/SingleMu/9_3_14_FullEta_v1/'
#path = '/afs/cern.ch/work/k/kbunkow/public/data/SingleMuFullEta/721_FullEta_v4/'

onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
#print onlyfiles

filesNameLike = sys.argv[2]
#chosenFiles = ['file://' + path + f for f in onlyfiles if (('_p_10_' in f) or ('_m_10_' in f))]
#chosenFiles = ['file://' + path + f for f in onlyfiles if (('_10_p_10_' in f))]
#chosenFiles = ['file://' + path + f for f in onlyfiles if (re.match('.*_._p_10.*', f))]
#chosenFiles = ['file://' + path + f for f in onlyfiles if ((filesNameLike in f))]

#print onlyfiles

chosenFiles = []

filesPerPtBin = 10 #TODO max is 200 for the 721_FullEta_v4 and 100 for 9_3_14_FullEta_v2

if filesNameLike == 'allPt' :
    for ptCode in range(31, 4, -1) : #the rigt bound of range is not included 
        if ptCode == 5 : #5 is 3-4 GeV (maybe 3-3.5 GeV). 4 is 2-3GeV (maybe 2.5-3 GeV), very small fraction makes candidates, and even less reaches the second station
            filesPerPtBin = 30
        elif ptCode == 6 : #5 is 3-4 GeV (maybe 3-3.5 GeV). 4 is 2-3GeV (maybe 2.5-3 GeV), very small fraction makes candidates, and even less reaches the second station
            filesPerPtBin = 20    
        elif ptCode <= 7 : 
            filesPerPtBin = 10
        elif ptCode <= 12 :
            filesPerPtBin = 5
        else :    
            filesPerPtBin = 3
            
        for sign in ['_m', '_p'] : #, m
            selFilesPerPtBin = 0
            for i in range(1, 201, 1): #TODO
                for f in onlyfiles:
                   if (( '_' + str(ptCode) + sign + '_' + str(i) + '_') in f): #TODO for 721_FullEta_v4/
                   #if (( '_' + str(ptCode) + sign + '_' + str(i) + ".") in f):  #TODO for 9_3_14_FullEta_v2
                        #print f
                        chosenFiles.append('file://' + path + f) 
                        selFilesPerPtBin += 1
                if(selFilesPerPtBin >= filesPerPtBin):
                    break
                        
else :
    for i in range(1, 2, 1):
        for f in onlyfiles:
            if (( filesNameLike + '_' + str(i) + '_') in f):  #TODO for 721_FullEta_v4/
            #if (( filesNameLike + '_' + str(i) + '.') in f): #TODO for 9_3_14_FullEta_v2
                print f
                chosenFiles.append('file://' + path + f) 
         

print "chosenFiles"
for chFile in chosenFiles:
    print chFile

if len(chosenFiles) == 0 :
    print "no files selected!!!!!!!!!!!!!!!"
    exit

firstEv = 0#40000
#nEvents = 1000

# input files (up to 255 files accepted)
process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring( 
    #'file:/eos/user/k/kbunkow/cms_data/SingleMuFullEta/721_FullEta_v4/SingleMu_16_p_1_1_xTE.root',
    #'file:/afs/cern.ch/user/k/kpijanow/Neutrino_Pt-2to20_gun_50.root',
    list(chosenFiles), ),
    skipEvents =  cms.untracked.uint32(0),
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
        'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016s_simEmtfDigis__HLT')
)

	                    
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))


####Event Setup Producer
process.load('L1Trigger.L1TMuonOverlapPhase1.fakeOmtfParams_cff')
#process.omtfParams.configXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/hwToLogicLayer_0x0008_patGen.xml"),

process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
   toGet = cms.VPSet(
      cms.PSet(record = cms.string('L1TMuonOverlapParamsRcd'),
               data = cms.vstring('L1TMuonOverlapParams'))
                   ),
   verbose = cms.untracked.bool(False)
)

#process.TFileService = cms.Service("TFileService", fileName = cms.string('omtfAnalysis1_1.root'), closeFileFast = cms.untracked.bool(True) )
                                   
####OMTF Emulator
process.load('L1Trigger.L1TMuonOverlapPhase1.simOmtfDigis_cfi')

process.simOmtfDigis.bxMin = cms.int32(0)
process.simOmtfDigis.bxMax = cms.int32(0)

process.simOmtfDigis.dumpResultToXML = cms.bool(False)
process.simOmtfDigis.dumpResultToROOT = cms.bool(False)
process.simOmtfDigis.eventCaptureDebug = cms.bool(False)

process.simOmtfDigis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/Patterns_template.xml")
#process.simOmtfDigis.patternsXMLFiles = cms.VPSet(cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/GPs_parametrised_plus_v1.xml")),
#                                                       cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/GPs_parametrised_minus_v1.xml"))
#)

process.simOmtfDigis.patternGenerator = cms.string("patternGen")

process.simOmtfDigis.patternType = cms.string("GoldenPatternWithStat")
process.simOmtfDigis.generatePatterns = cms.bool(True)
process.simOmtfDigis.optimisedPatsXmlFile = cms.string("Patterns_0x00011_oldSample_3_30Files_layerStat_test.xml")

process.simOmtfDigis.rpcMaxClusterSize = cms.int32(3)
process.simOmtfDigis.rpcMaxClusterCnt = cms.int32(2)
process.simOmtfDigis.rpcDropAllClustersIfMoreThanMax = cms.bool(True)

process.simOmtfDigis.goldenPatternResultFinalizeFunction = cms.int32(3) ## is needed here , becasue it just counts the number of layers with a stub
process.simOmtfDigis.lctCentralBx = cms.int32(6);#<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!TODO this was changed in CMSSW 10(?) to 8. if the data were generated with the previous CMSSW then you have to use 6


#process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")
#process.dumpES = cms.EDAnalyzer("PrintEventSetupContent")

process.L1TMuonSeq = cms.Sequence( process.esProd          
                                   + process.simOmtfDigis 
                                   #+ process.dumpED
                                   #+ process.dumpES
)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)

#process.out = cms.OutputModule("PoolOutputModule", 
#   fileName = cms.untracked.string("l1tomtf_superprimitives1.root")
#)

#process.output_step = cms.EndPath(process.out)
#process.schedule = cms.Schedule(process.L1TMuonPath)
#process.schedule.extend([process.output_step])
