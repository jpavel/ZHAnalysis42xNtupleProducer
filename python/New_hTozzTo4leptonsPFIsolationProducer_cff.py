# Electron PF isolation
from CommonTools.ParticleFlow.PFBRECO_cff import *
from CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff import *
elPFIsoDepositCharged.src    = cms.InputTag("selectedPatElectrons")
elPFIsoDepositChargedAll.src = cms.InputTag("selectedPatElectrons")
elPFIsoDepositNeutral.src    = cms.InputTag("selectedPatElectrons")
elPFIsoDepositGamma.src      = cms.InputTag("selectedPatElectrons")
elPFIsoDepositPU.src         = cms.InputTag("selectedPatElectrons")


# Muon PF isolation
from CommonTools.ParticleFlow.PFBRECO_cff import *
from CommonTools.ParticleFlow.Isolation.pfMuonIsolation_cff import *
muPFIsoDepositCharged.src    = cms.InputTag("selectedPatMuons")
muPFIsoDepositChargedAll.src = cms.InputTag("selectedPatMuons")
muPFIsoDepositNeutral.src    = cms.InputTag("selectedPatMuons")
muPFIsoDepositGamma.src      = cms.InputTag("selectedPatMuons")
muPFIsoDepositPU.src         = cms.InputTag("selectedPatMuons")

