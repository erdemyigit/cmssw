import FWCore.ParameterSet.Config as cms

from DQMOffline.EGamma.electronOfflineClient_cfi import *

dqmElectronClientAllElectrons = dqmElectronOfflineClient.clone() ;
dqmElectronClientAllElectrons.InputFolderName = cms.string("Egamma/Electrons/Ele2_All") ;
dqmElectronClientAllElectrons.OutputFolderName = cms.string("Egamma/Electrons/Ele2_All") ;

dqmElectronClientAllElectronsHGC = dqmElectronOfflineClient.clone() ;
dqmElectronClientAllElectronsHGC.InputFolderName = cms.string("Egamma/Electrons/Ele2HGC_All") ;
dqmElectronClientAllElectronsHGC.OutputFolderName = cms.string("Egamma/Electrons/Ele2HGC_All") ;

dqmElectronClientSelectionEt = dqmElectronOfflineClient.clone() ;
dqmElectronClientSelectionEt.InputFolderName = cms.string("Egamma/Electrons/Ele3_Et10") ;
dqmElectronClientSelectionEt.OutputFolderName = cms.string("Egamma/Electrons/Ele3_Et10") ;

dqmElectronClientSelectionEtIso = dqmElectronOfflineClient.clone() ;
dqmElectronClientSelectionEtIso.InputFolderName = cms.string("Egamma/Electrons/Ele4_Et10TkIso1") ;
dqmElectronClientSelectionEtIso.OutputFolderName = cms.string("Egamma/Electrons/Ele4_Et10TkIso1") ;

#dqmElectronClientSelectionEtIsoElID = dqmElectronOfflineClient.clone() ;
#dqmElectronClientSelectionEtIsoElID.InputFolderName = cms.string("Egamma/Electrons/Ele4_Et10TkIso1ElID") ;
#dqmElectronClientSelectionEtIsoElID.OutputFolderName = cms.string("Egamma/Electrons/Ele4_Et10TkIso1ElID") ;

dqmElectronClientTagAndProbe = dqmElectronOfflineClient.clone() ;
dqmElectronClientTagAndProbe.InputFolderName = cms.string("Egamma/Electrons/Ele5_TagAndProbe") ;
dqmElectronClientTagAndProbe.OutputFolderName = cms.string("Egamma/Electrons/Ele5_TagAndProbe") ;
dqmElectronClientTagAndProbe.EffHistoTitle = cms.string("")

electronOfflineClientSequence = cms.Sequence(
   dqmElectronClientAllElectrons
 * dqmElectronClientSelectionEt
 * dqmElectronClientSelectionEtIso
# * dqmElectronClientSelectionEtIsoElID
 * dqmElectronClientTagAndProbe
)
_electronOfflineClientSequenceHGC = electronOfflineClientSequence.copy()
_electronOfflineClientSequenceHGC += dqmElectronClientAllElectronsHGC

from Configuration.Eras.Modifier_phase2_hgcal_cff import phase2_hgcal
phase2_hgcal.toReplaceWith(
  electronOfflineClientSequence, _electronOfflineClientSequenceHGC
)

