import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var,P3Vars

hgcEESimHitsTable = cms.EDProducer("SimplePCaloHitFlatTableProducer",
    src = cms.InputTag("g4SimHits:HGCHitsEE"),
    cut = cms.string(""), 
    name = cms.string("SimHitHGCEE"),
    doc  = cms.string("Geant4 SimHits in HGCAL Electromagnetic endcap"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('id', 'int', precision=-1, doc='detId'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        trackId = Var('geantTrackId', 'int', precision=-1, doc='Geant4 track ID'),
    )
)

hgcEEHitsToSimClusterTable = cms.EDProducer("CaloHitToSimClusterIndexTableProducer",
    cut = hgcEESimHitsTable.cut,
    src = hgcEESimHitsTable.src,
    objName = hgcEESimHitsTable.name,
    branchName = cms.string("SimCluster"),
    objMap = cms.InputTag("mix:simHitHGCEEToSimCluster"),
    docString = cms.string("SimCluster containing SimHit")
)

hgcHEfrontSimHitsTable = hgcEESimHitsTable.clone()
hgcHEfrontSimHitsTable.src = "g4SimHits:HGCHitsHEfront"
hgcHEfrontSimHitsTable.name = "SimHitHGCHEF"

hgcHEfrontHitsToSimClusterTable = hgcEEHitsToSimClusterTable.clone()
hgcHEfrontHitsToSimClusterTable.src = hgcHEfrontSimHitsTable.src
hgcHEfrontHitsToSimClusterTable.objName = hgcHEfrontSimHitsTable.name
hgcHEfrontHitsToSimClusterTable.objMap = "mix:simHitHGCHEfrontToSimCluster"

hgcHEbackSimHitsTable = hgcEESimHitsTable.clone()
hgcHEbackSimHitsTable.src = "g4SimHits:HGCHitsHEback"
hgcHEbackSimHitsTable.name = "SimHitHGCHEB"

hgcHEbackHitsToSimClusterTable = hgcEEHitsToSimClusterTable.clone()
hgcHEbackHitsToSimClusterTable.src = hgcHEbackSimHitsTable.src
hgcHEbackHitsToSimClusterTable.objName = hgcHEbackSimHitsTable.name
hgcHEbackHitsToSimClusterTable.objMap = "mix:simHitHGCHEbackToSimCluster"

HcalSimHitsTable = hgcEESimHitsTable.clone()
HcalSimHitsTable.src = "g4SimHits:HcalHits"
HcalSimHitsTable.name = "SimHitHCAL"

HcalHitsToSimClusterTable = hgcEEHitsToSimClusterTable.clone()
HcalHitsToSimClusterTable.src = HcalSimHitsTable.src
HcalHitsToSimClusterTable.objName = HcalSimHitsTable.name
HcalHitsToSimClusterTable.objMap = "mix:simHitHcalToSimCluster"

EcalEBSimHitsTable = hgcEESimHitsTable.clone()
EcalEBSimHitsTable.src = "g4SimHits:EcalHitsEB"
EcalEBSimHitsTable.name = "SimHitEcalEB"

EcalEBHitsToSimClusterTable = hgcEEHitsToSimClusterTable.clone()
EcalEBHitsToSimClusterTable.src = EcalEBSimHitsTable.src
EcalEBHitsToSimClusterTable.objName = EcalEBSimHitsTable.name
EcalEBHitsToSimClusterTable.objMap = "mix:simHitEcalToSimCluster"

hgcEESimHitsPositionTable = cms.EDProducer("HGCalSimHitPositionTableProducer",
    src = hgcEESimHitsTable.src,
    cut = hgcEESimHitsTable.cut, 
    name = hgcEESimHitsTable.name,
    doc  = hgcEESimHitsTable.doc,
)

hgcHEfrontSimHitsPositionTable = hgcEESimHitsPositionTable.clone()
hgcHEfrontSimHitsPositionTable.name = hgcHEfrontSimHitsTable.name
hgcHEfrontSimHitsPositionTable.src = hgcHEfrontSimHitsTable.src

hgcHEbackSimHitsPositionTable = hgcEESimHitsPositionTable.clone()
hgcHEbackSimHitsPositionTable.name = hgcHEbackSimHitsTable.name
hgcHEbackSimHitsPositionTable.src = hgcHEbackSimHitsTable.src

HcalSimHitsPositionTable = hgcEESimHitsPositionTable.clone()
HcalSimHitsPositionTable.name = HcalSimHitsTable.name
HcalSimHitsPositionTable.src = HcalSimHitsTable.src

EcalEBSimHitsPositionTable = hgcEESimHitsPositionTable.clone()
EcalEBSimHitsPositionTable.name = EcalEBSimHitsTable.name
EcalEBSimHitsPositionTable.src = EcalEBSimHitsTable.src

hgcSimHitsSequence = cms.Sequence(
                hgcEESimHitsTable
                +hgcHEbackSimHitsTable
                +hgcHEfrontSimHitsTable
                +HcalSimHitsTable
                +EcalEBSimHitsTable
                
                +hgcEESimHitsPositionTable
                +hgcHEfrontSimHitsPositionTable
                +hgcHEbackSimHitsPositionTable
                +HcalSimHitsPositionTable
                +EcalEBSimHitsPositionTable

                #+hgcEEHitsToSimClusterTable
                #+hgcHEfrontHitsToSimClusterTable
                #+hgcHEbackHitsToSimClusterTable
                +hgcHEfrontSimHitsTable
                +hgcHEbackSimHitsTable)

