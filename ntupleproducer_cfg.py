import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleProducer")

process.out = cms.OutputModule("PoolOutputModule",
                               fileName=cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
                               # save only events passing the full path
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                               outputCommands=cms.untracked.vstring('drop *')
                               )

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("RecoEcal.EgammaClusterProducers.ecalClusteringSequence_cff")
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
#process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi")

#trigger
#process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
#process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerMatchEmbedder_cfi")
#process.load("TauAnalysis.RecoTools.recoVertexSelection_cff")
process.load("Analysis.NtupleProducer.New_hTozzTo4leptonsPFIsolationProducer_cff")
#process.load("Analysis.NtupleProducer.JES_Uncertainty_FR_cff")
process.MessageLogger = cms.Service("MessageLogger")
process.load("RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi")
process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.load('CommonTools.ParticleFlow.pfNoPileUp_cff')
# load the PU JetID sequence
process.load("CMGTools.External.pujetidsequence_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load('JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_42X_cff') #MVA MET
#################################################   Samples and GlobalTag   ############################
#Source File
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
                            #                            'file:/tmp/abdollah/WToTauNu_TuneZ2_AODSIM_E7TeV_FlatDist10_2011EarlyData_50ns_START311_V1G1.root')
                            #                            'file:/tmp/abdollah/WToTauNu_TuneZ2_7TeV_Spring11_AODSIM_PU_S1_START311_V1G1.root')data2011F66C812D-E28C-E111-B27E-001D09F28F25.root
#                            'file:/afs/cern.ch/work/a/abdollah/data2011F66C812D-E28C-E111-B27E-001D09F28F25.root')
#'file:/tmp/abdollah/5CC144A1-42F0-E011-B6C1-00215E21DC90.root')
    'file:/afs/cern.ch/work/j/jez/ntuples/analysis/mc/Summer11/WH_ZH_TTH_HToTauTau_M-130_7TeV-pythia6-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0EC8AE34-128E-E011-A1F3-90E6BA19A25E.root'
   # 'file:/afs/cern.ch/work/a/abdollah/Sample/WH_ZH_TTH_HToTauTau_M-130_7TeV_PU_S6_START42_V14B-v1.root'
    )
 
#                                                                                 'file:/tmp/abdollah/E23A3C95-74BA-E011-8A75-00237DA435CA.root')
                            #                            'file:/tmp/evan_zz_fall11_1k_numEvent1000.root')
                            #                            'file:/afs/cern.ch/user/f/friis/public/DoubleElectron-Run2011A-05Aug2011-v1-AOD-uw_eemt_l1_anti_iso.root',
                            #                                                        'file:/tmp/abdollah/skim_WtoTauNu_data_RunData42Xv4.root')
                            #'rfio:/castor/cern.ch/user/a/abdollah/HLT/TrigEfficiency/Wtaunu/TauHLTOutput_Fall10_981_2_F9K.root')
                            )
process.maxEvents = cms.untracked.PSet(
                                       input=cms.untracked.int32(100)
                                       )

isMC = True      #comment it for Data
#isMC = False     #comment it for MC

if isMC:
    HLTProcessName = 'HLT'
    process.GlobalTag.globaltag = cms.string('START42_V14B::All') # Redigi of summer 11 samples with more pileup 
#    process.GlobalTag.globaltag = cms.string('START41_V0::All')
#    process.GlobalTag.globaltag = cms.string('START311_V2::All')

else:
    HLTProcessName = 'HLT'
    process.GlobalTag.globaltag = cms.string('FT42_V24_AN1::All') #Jan 16 processing
#    process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')
#    process.GlobalTag.globaltag = cms.string('GR_R_311_V2::All')



#################################################   EDANALYZER   ##################################
process.myanalysis = cms.EDAnalyzer("NtupleProducer",

                                    #OutPut Filles (NTuples)
                                    HistOutFile=cms.untracked.string('output_Ntuples.root'),

                                    # Include or Exclude the objects
                                    Include_HPSTau=cms.bool(True),
                                    Include_Muon=cms.bool(True),
                                    Include_Electron=cms.bool(True),
                                    Include_Jet=cms.bool(True),
                                    Include_JetCorrection=cms.bool(False),
                                    Include_Jet_ESU=cms.bool(False),
                                    Include_MET=cms.bool(True),
                                    Include_GenPartiles=cms.bool(True),
                                    Include_HLT=cms.bool(True),
                                    Include_Vertex=cms.bool(True),
                                    Is_MC=cms.bool(isMC),
                                    # storing only certain trigger strings
                                    filterTriggerResults=cms.bool(True),
                                    el_trigger_name = cms.string("Ele"),
                                    mu_trigger_name = cms.string("Mu"),
                                    #vetrex and Track
                                    vertices=cms.InputTag('offlinePrimaryVertices'),
                                    tracks=cms.InputTag("generalTracks"),

                                    #MET
                                    met=cms.InputTag("met"),
                                    PFmet=cms.InputTag("pfMet"),
                                    tcmet=cms.InputTag("tcMet"),
                                    Type1CorMET=cms.InputTag("patType1CorrectedPFMet"),
                                    MVAmet=cms.InputTag("pfMEtMVA"),
                                    #Jets
                                    bjets=cms.InputTag("trackCountingHighPurBJetTags"),
                                    PFAK5=cms.InputTag("selectedPatJets"),
                                    rhoJetsLabel=cms.InputTag("kt6PFJets", "rho"),

                                    #Leptons
                                    preselectedelectrons=cms.InputTag("selectedPatElectrons"),
                                    preselectedmuons=cms.InputTag("selectedPatMuons"),
                                    preselectedHPSTaus=cms.InputTag("selectedPatTaus"),
                                    tauPtCut=cms.double(15.0),

#                                    # CIC electron Identification
#                                     eleID_VeryLooseTag=cms.InputTag("eidVeryLoose"),
#                                     eleID_LooseTag=cms.InputTag("eidLoose"),
#                                     eleID_MediumTag=cms.InputTag("eidMedium"),
#                                     eleID_TightTag=cms.InputTag("eidTight"),

                                    #Trigger


                                    srcTriggerResults=cms.InputTag("TriggerResults", "", HLTProcessName),

                                    # MC Information
                                    PileUpInfo=cms.InputTag("addPileupInfo"),
                                    genParticlesInfo=cms.InputTag("genParticles"),

                                                                      #Trigger and TriggerMatching
#                                    trigger=cms.InputTag("patTrigger"),
                                    triggerEvent=cms.InputTag("patTriggerEvent"),
#                                    taus=cms.InputTag("selectedPatTaus"),
                                    tauMatch_Loose=cms.string('tauTriggerMatchHLTTausLoose'),
                                    tauMatch_Medium=cms.string('tauTriggerMatchHLTTausMedium'),
                                    electronMatch_Loose=cms.string('electronTriggerMatchHLTElectronsLoose'),
                                    muonMatch_Loose=cms.string('muonTriggerMatchHLTMuonsLoose'),
                                    puJetIdFlag=cms.InputTag("puJetMva","fullId"),
                                    #rho
                                    rhoProducer = cms.InputTag('kt6PFJetsForRhoComputationVoronoi','rho')

                                    )
#################################################   PFIsolation  ################################
from CommonTools.ParticleFlow.pfParticleSelection_cff import *
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
process.eleIsoSequence = setupPFElectronIso(process, 'selectedPatElectrons')

# cone vetos as used in Htautau
process.elPFIsoValueChargedAll04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueGamma04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')
process.elPFIsoValueGamma04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')


#################################################   scrapingVeto  ################################
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter=cms.untracked.bool(True),
                                    debugOn=cms.untracked.bool(False),
                                    numtrack=cms.untracked.uint32(10),
                                    thresh=cms.untracked.double(0.25)
                                    )
if not isMC:
    process.dataFilter = cms.Sequence( process.scrapingVeto )
else:
    process.dataFilter = cms.Sequence()


#################################################   Good Primary Vertex ################################
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodVertices = cms.EDFilter(
                                    "PrimaryVertexObjectFilter",
                                    filterParams=pvSelector.clone(minNdof=cms.double(7.0), maxZ=cms.double(24.0)),
                                    src=cms.InputTag('offlinePrimaryVertices')
                                    )

#################################################   PAT APPLICATIONS   ##################################

#Removing MC Matching
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'])

#Switch to ak5PFJets (L2 and L3 Corrections are included)
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,
                    cms.InputTag('ak5PFJets'),
                    doJTA=True,
                    doBTagging=True,
                    jetCorrLabel=('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                    #                    jetCorrLabel=('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                    #                    jetCorrLabel=('AK5PF', cms.vstring(['L2Relative', 'L3Absolute'])),
                    doType1MET=True,
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID=True
                    )
if not isMC:
    process.patJetCorrFactors.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual')

### Jet Corrections ##########################################################
# See https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookJetEnergyCorrections
# note: this runs the L1Fast-Jet corrections for PF jets. not applied on Calo
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#..process.load('RecoJets.Configuration.RecoPFJets_cff')
# Turn-on the FastJet density calculation -----------------------
#### FastJet corrections
from RecoJets.JetProducers.kt4PFJets_cfi import *
import RecoJets.JetProducers.kt4PFJets_cfi
process.kt6PFJets = RecoJets.JetProducers.kt4PFJets_cfi.kt4PFJets.clone()
process.kt6PFJets.rParam = cms.double(0.6)
process.kt6PFJets.doRhoFastjet = cms.bool(True)
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

from RecoJets.Configuration.RecoPFJets_cff import *
process.kt6PFJetsCentral = kt6PFJets.clone(
    Ghost_EtaMax = cms.double(2.5),
    Rho_EtaMax = cms.double(2.5),
)
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
#process.load('RecoJets.Configuration.RecoPFJets_cff')
#process.load('RecoJets.JetProducers.ak5PFJets_cfi')
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
#process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets = process.kt4PFJets.clone()
process.kt6PFJets.rParam = 0.6
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)
process.kt6PFJets.Ghost_EtaMax = cms.double(2.5)
process.kt6PFJets.doAreaFastjet = True

#CV: need to rerun 'ak5PFJets' module with jet area computation enabled,
#    since it has not been enabled per default in CMSSW_4_2_x
#   (if the jet area of 'ak5PFJets' is zero, the L1FastjetCorrector::correction function always returns 1.0)
#process.load("RecoJets.JetProducers.ak5PFJets_cfi")
process.ak5PFJets.doAreaFastjet = cms.bool(True)



##############changed
# switch to hadron-plus-strip(s) (HPS) PFTau collection (HPS is defuat PATTau now)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
#switchToPFTauHPS(process,
#                 cms.InputTag('shrinkingConePFTauProducer'),
#                 cms.InputTag('hpsPFTauProducer'),
#                 patTauLabel="",
#                 postfix="")
#process.cleanPatTaus.preselection = cms.string('pt > 15 & abs(eta) < 2.3 & tauID("decayModeFinding") > 0.5 & tauID("byLooseIsolation") > 0.5 & tauID("againstMuonTight") > 0.5 & tauID("againstElectronLoose") > 0.5')

## Add CiC electron ID's and WP eleID's
#process.load("cicEleIdSequence_42_cff")
#process.load("hzzEleIdSequence_cff")
#process.load("simpleEleIdSequence_cff")

# Electron Id
#from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi import *
#process.CiC_Ele_Id = cms.Sequence(process.eidVeryLoose + process.eidLoose + process.eidMedium + process.eidTight + process.eidSuperTight)

# Electron Id
#from RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi import *
#process.CiC_Ele_Id_HZZ = cms.Sequence(process.eidHZZVeryLoose + process.eidHZZLoose + process.eidHZZMedium + process.eidHZZTight + process.eidHZZSuperTight)
# Adding different WP for electron

process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectrons.electronIDSources = cms.PSet(
    mvaTrigV0=cms.InputTag("mvaTrigV0"),
     mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    )

# process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
# process.makePatElectrons = cms.Sequence(process.patElectronIDs * process.patElectronIsolation * process.patElectrons)
# process.patElectrons.addElectronID = cms.bool(True)
# process.patElectrons.electronIDSources = cms.PSet(
#                                                   simpleEleId95relIso=cms.InputTag("simpleEleId95relIso"),
#                                                   simpleEleId90relIso=cms.InputTag("simpleEleId90relIso"),
#                                                   simpleEleId85relIso=cms.InputTag("simpleEleId85relIso"),
#                                                   simpleEleId80relIso=cms.InputTag("simpleEleId80relIso"),
#                                                   simpleEleId70relIso=cms.InputTag("simpleEleId70relIso"),
#                                                   simpleEleId60relIso=cms.InputTag("simpleEleId60relIso"),

#                                                   simpleEleId95cIso=cms.InputTag("simpleEleId95cIso"),
#                                                   simpleEleId90cIso=cms.InputTag("simpleEleId90cIso"),
#                                                   simpleEleId85cIso=cms.InputTag("simpleEleId85cIso"),
#                                                   simpleEleId80cIso=cms.InputTag("simpleEleId80cIso"),
#                                                   simpleEleId70cIso=cms.InputTag("simpleEleId70cIso"),
#                                                   simpleEleId60cIso=cms.InputTag("simpleEleId60cIso"),

#                                                   eidVeryLoose=cms.InputTag("eidVeryLoose"),
#                                                   eidLoose=cms.InputTag("eidLoose"),
#                                                   eidMedium=cms.InputTag("eidMedium"),
#                                                   eidTight=cms.InputTag("eidTight"),
#                                                   eidSuperTight=cms.InputTag("eidSuperTight"),

#                                                   eidHZZVeryLoose=cms.InputTag("eidHZZVeryLoose"),
#                                                   eidHZZLoose=cms.InputTag("eidHZZLoose"),
#                                                   eidHZZMedium=cms.InputTag("eidHZZMedium"),
#                                                   eidHZZTight=cms.InputTag("eidHZZTight"),
#                                                   eidHZZSuperTight=cms.InputTag("eidHZZSuperTight"),

#                                                   )



#Trigger matching "HLT_SingleIsoTau20_Trk15_MET25_v*"
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.tauTriggerMatchHLTTausLoose = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
                                                src=cms.InputTag("selectedPatTaus"),
                                                matched=cms.InputTag("patTrigger"),
                                                matchedCuts=cms.string('path( "HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v*" )'),
                                                maxDPtRel=cms.double(0.5),
                                                maxDeltaR=cms.double(0.5),
                                                resolveAmbiguities=cms.bool(True),
                                                resolveByMatchQuality=cms.bool(True)
                                                )
process.tauTriggerMatchHLTTausMedium = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
                                                src=cms.InputTag("selectedPatTaus"),
                                                matched=cms.InputTag("patTrigger"),
                                                matchedCuts=cms.string('path( "HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v*" )'),
                                                maxDPtRel=cms.double(0.5),
                                                maxDeltaR=cms.double(0.5),
                                                resolveAmbiguities=cms.bool(True),
                                                resolveByMatchQuality=cms.bool(True)
                                                )

process.electronTriggerMatchHLTElectronsLoose = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
                                                               src=cms.InputTag("selectedPatElectrons"),
                                                               matched=cms.InputTag("patTrigger"),
                                                               matchedCuts=cms.string('path( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*" )'),
                                                               maxDPtRel=cms.double(0.5),
                                                               maxDeltaR=cms.double(0.5),
                                                               resolveAmbiguities=cms.bool(True),
                                                               resolveByMatchQuality=cms.bool(True)
                                                               )

process.muonTriggerMatchHLTMuonsLoose = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
                                                       src=cms.InputTag("selectedPatMuons"),
                                                       matched=cms.InputTag("patTrigger"),
                                                       matchedCuts=cms.string('path( "HLT_Mu17_Mu8_v*" ) || path( "HLT_Mu17_TkMu8_v*" )'),
                                                       maxDPtRel=cms.double(0.5),
                                                       maxDeltaR=cms.double(0.5),
                                                       resolveAmbiguities=cms.bool(True),
                                                       resolveByMatchQuality=cms.bool(True)
                                                       )


process.pfPileUp.Enable = cms.bool(True)

process.AddPUInfo = process.myanalysis.clone()
process.AddGenInfo = process.myanalysis.clone()
if not isMC:
    process.AddPUInfo = process.myanalysis.clone(PileUpInfo=cms.InputTag(""))
    process.AddGenInfo = process.myanalysis.clone(genParticlesInfo=cms.InputTag(""))


##From MOnica
process.pfPileUp.PFCandidates = 'particleFlow'
process.pfNoPileUp.bottomCollection = 'particleFlow'
process.pfPileUpIso.PFCandidates = 'particleFlow'
process.pfNoPileUpIso.bottomCollection = 'particleFlow'

process.pfAllPhotons.src=cms.InputTag("pfNoPileUp")
process.pfAllChargedHadrons.src=cms.InputTag("pfNoPileUp")
process.pfAllNeutralHadrons.src=cms.InputTag("pfNoPileUp")
process.pfAllChargedParticles.src=cms.InputTag("pfNoPileUp")

## MVA MET configuration
from JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_cfi import isotaus

process.isotaus.discriminators = cms.VPSet(
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),       selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"), selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"), selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),     selectionCut=cms.double(0.5))
    )

#################################################   PATH of CONFIG FILE ##################################
#print process.dumpPython()
process.p = cms.Path (
                      process.dataFilter *
                      process.goodVertices *
                      process.PFTau * #prescribed by TAU POG
                      process.mvaID   *
                      process.pfNoPileUpSequence *
#                      process.makePatElectrons *
#                      process.ak5PFJets *
#                      process.kt6PFJets  *
#                      process.producePFMETCorrections*
                      #                      process.kt6PFJetsCentral *                     
                      process.patDefaultSequence +
#                      process.pfPileUpIso*
                      process.pfParticleSelectionSequence +                
                      process.eleIsoSequence*
#                      process.pfElectronIsolationSequence                 *
                      process.pfMuonIsolationSequence                     *
#                      process.pfPostSequence *
#                      process.produceJEC_Uncertainty_FR *
#                      process.electronPFIsoDeposits *
#                      process.electronPFIsoValues *
#                      process.muPFIsoDeposits *
#                      process.muPFIsoValues *
                      process.puJetIdSqeuence * 
                      process.pfMEtMVAsequence *
                      process.myanalysis *
                      process.AddPUInfo *
                      process.AddGenInfo
                      )

#################################################   NEEDED FOR TRIGGER MATCHING   #######################



from PhysicsTools.PatAlgos.tools.trigTools import *

switchOnTrigger( process ) # This is optional and can be omitted.

switchOnTriggerMatching(process, ['tauTriggerMatchHLTTausLoose','tauTriggerMatchHLTTausMedium','electronTriggerMatchHLTElectronsLoose','muonTriggerMatchHLTMuonsLoose'])
#switchOnTriggerMatching(process, ['tauTriggerMatchHLTTausLoose'])
#switchOnTriggerMatching(process, ['tauTriggerMatchHLTTausMedium'])

# Switch to selected PAT objects in the trigger matching
removeCleaningFromTriggerMatching( process ) #,outputModule='out') #process,'out')

#switchOnTrigger(process)
process.patTrigger.processName = HLTProcessName
process.patTriggerEvent.processName = HLTProcessName
