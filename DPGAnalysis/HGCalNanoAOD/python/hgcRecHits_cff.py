import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var,P3Vars

hgcRecHits = cms.EDProducer("HGCRecHitCollectionMerger",
	src = cms.VInputTag("HGCalRecHit:HGCEERecHits",
        "HGCalRecHit:HGCHEFRecHits", "HGCalRecHit:HGCHEBRecHits",
	)	
)

hgcRecHitsTable = cms.EDProducer("SimpleCaloRecHitFlatTableProducer",
    src = cms.InputTag("hgcRecHits"),
    cut = cms.string(""), 
    name = cms.string("RecHitHGC"),
    doc  = cms.string("HGCAL RecHits"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(
        detId = Var('detid().rawId()', 'int', precision=-1, doc='detId'),
        energy = Var('energy', 'float', precision=14, doc='energy'),
        time = Var('time', 'float', precision=14, doc='hit time'),
    )
)

hgcRecHitsToLayerClusters = cms.EDProducer("RecHitToLayerClusterAssociationProducer",
    caloRecHits = cms.VInputTag("hgcRecHits"),
    layerClusters = cms.InputTag("hgcalLayerClusters"),
)

hgcRecHitsToLayerClusterTable = cms.EDProducer("HGCRecHitToLayerClusterIndexTableProducer",
    cut = hgcRecHitsTable.cut,
    src = hgcRecHitsTable.src,
    objName = hgcRecHitsTable.name,
    branchName = cms.string("LayerCluster"),
    objMap = cms.InputTag("hgcRecHitsToLayerClusters:hgcRecHitsToLayerCluster"),
    docString = cms.string("LayerCluster assigned largest RecHit fraction")
)

hgcRecHitsPositionTable = cms.EDProducer("HGCALRecHitPositionTableProducer",
    src = hgcRecHitsTable.src,
    cut = hgcRecHitsTable.cut, 
    name = hgcRecHitsTable.name,
    doc  = hgcRecHitsTable.doc,
)

hgcRecHitsSequence = cms.Sequence(hgcRecHits
				+hgcRecHitsTable
                +hgcRecHitsToLayerClusters 
                +hgcRecHitsToLayerClusterTable
                +hgcRecHitsPositionTable
)
