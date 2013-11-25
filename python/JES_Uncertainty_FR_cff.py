import FWCore.ParameterSet.Config as cms
#--------------------------------------------------------------------------------------------
    # produce collection of jets shifted up/down in energy
    #--------------------------------------------------------------------------------------------

    # in case of "raw" (uncorrected) MET,
    # add residual jet energy corrections in quadrature to jet energy uncertainties:
    # cf. https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
#from PhysicsTools.PatUtils.tools.metUncertaintyTools import *
#load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
jetsEnUpSubTotalDataMC = cms.EDProducer("ShiftedPATJetProducer",
                                        src=cms.InputTag("selectedPatJets"),
                                        #jetCorrPayloadName = cms.string(jetCorrPayloadName),
                                        #jetCorrUncertaintyTag = cms.string('Uncertainty'),
                                        jetCorrInputFileName=cms.FileInPath('PhysicsTools/PatUtils/data/JEC11_V12_AK5PF_UncertaintySources.txt'),
                                        jetCorrUncertaintyTag=cms.string("SubTotalDataMC"),
                                        addResidualJES=cms.bool(True),
                                        jetCorrLabelUpToL3=cms.string("ak5PFL1FastL2L3"),
                                        jetCorrLabelUpToL3Res=cms.string("ak5PFL1FastL2L3Residual"),
                                        shiftBy=cms.double(+ 1. * 1)
                                        )

jetsEnDownSubTotalDataMC = jetsEnUpSubTotalDataMC.clone(
                                                        shiftBy=cms.double(-1. * 1)
                                                        )

jetsEnUpTotal = jetsEnUpSubTotalDataMC.clone(
                                             jetCorrUncertaintyTag=cms.string("Total"),
                                             )


jetsEnDownTotal = jetsEnUpSubTotalDataMC.clone(
                                               jetCorrUncertaintyTag=cms.string("Total"),
                                               shiftBy=cms.double(-1. * 1)
                                               )



produceJEC_Uncertainty_FR = cms.Sequence(
                                         jetsEnUpSubTotalDataMC
                                         * jetsEnDownSubTotalDataMC
                                         * jetsEnUpTotal
                                         * jetsEnDownTotal
                                         )