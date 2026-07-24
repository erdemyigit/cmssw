import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCalUtilities.HGC3DClusterSelectors_cfi import *

L1THGCalTriggerSelector = cms.Sequence(l1tHGCal3DClusterGenMatchSelector+ l1tHGCal3DClusterDummySelector)
