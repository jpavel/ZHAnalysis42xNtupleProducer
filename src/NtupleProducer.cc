// -*- C++ -*-
//
// Package:    NtupleProducer
// Class:      NtupleProducer
// 
/**\class NtupleProducer NtupleProducer.cc Analysis/NtupleProducer/src/NtupleProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Abdollah Mohammadi,22 R-013,+41227672649,
//         Created:  Wed Apr 20 14:56:20 CEST 2011
// $Id: NtupleProducer.cc,v 1.4 2013/02/22 00:52:37 jez Exp $
//
//



#include "Analysis/NtupleProducer/interface/NtupleProducer.h"

// ------------ method called once each job just before starting event loop  ------------

void
NtupleProducer::beginJob() {
    hOutputFile = new TFile(fOutputFileName.c_str(), "RECREATE");
    t = new TTree("t", "tree");
    m = new myevent;
    t->Branch("myevent", "myevent", &m, 256000, 1);

}

// ------------ method called once each job just after ending the event loop  ------------

void
NtupleProducer::endJob() {

    hOutputFile->Write();
    hOutputFile->Close();
}

NtupleProducer::NtupleProducer(const edm::ParameterSet& iConfig) {

    fOutputFileName = iConfig.getUntrackedParameter<string > ("HistOutFile");

    Include_Vertex = iConfig.getParameter<bool>("Include_Vertex");
    Include_Electron = iConfig.getParameter<bool>("Include_Electron");
    Include_Muon = iConfig.getParameter<bool>("Include_Muon");
    Include_HPSTau = iConfig.getParameter<bool>("Include_HPSTau");
    Include_JetCorrection = iConfig.getParameter<bool>("Include_JetCorrection");
    Include_GenPartiles = iConfig.getParameter<bool>("Include_GenPartiles");
    Include_Jet = iConfig.getParameter<bool>("Include_Jet");
    Include_Jet_ESU = iConfig.getParameter<bool>("Include_Jet_ESU");
    Include_MET = iConfig.getParameter<bool>("Include_MET");
    Include_HLT = iConfig.getParameter<bool>("Include_HLT");
    IsMC = iConfig.getParameter<bool>("Is_MC");

    TrackCollection_ = iConfig.getParameter<edm::InputTag > ("tracks");
    VertexCollection_ = iConfig.getParameter<edm::InputTag > ("vertices");
    rhoCorrection_ = iConfig.getParameter<edm::InputTag > ("rhoJetsLabel");

    JetCollection_ = iConfig.getParameter<edm::InputTag > ("PFAK5");
    BJetCollection_ = iConfig.getParameter<edm::InputTag > ("bjets");

    MetCollection_ = iConfig.getParameter<edm::InputTag > ("met");
    PFMetCollection_ = iConfig.getParameter<edm::InputTag > ("PFmet");
    tcMetCollection_ = iConfig.getParameter<edm::InputTag > ("tcmet");
    //    Type1CorMETCollection_ = iConfig.getParameter<edm::InputTag > ("Type1CorMET");
    MVAMetCollection_ = iConfig.getParameter<edm::InputTag > ("MVAmet");

    PreSelectedElectronCollection_ = iConfig.getParameter<edm::InputTag > ("preselectedelectrons");
    PreSelectedMuonCollection_ = iConfig.getParameter<edm::InputTag > ("preselectedmuons");
    PreSelectedhpsCollection_ = iConfig.getParameter<edm::InputTag > ("preselectedHPSTaus");


//     EleID_VeryLooseTag_ = iConfig.getParameter<edm::InputTag > ("eleID_VeryLooseTag");
//     EleID_LooseTag_ = iConfig.getParameter<edm::InputTag > ("eleID_LooseTag");
//     EleID_MediumTag_ = iConfig.getParameter<edm::InputTag > ("eleID_MediumTag");
//     EleID_TightTag_ = iConfig.getParameter<edm::InputTag > ("eleID_TightTag");
    srcTriggerResults_ = iConfig.getParameter<edm::InputTag > ("srcTriggerResults");

    PileUpInfo_ = iConfig.getParameter<edm::InputTag > ("PileUpInfo");
    GenParticlesInfo_ = iConfig.getParameter<edm::InputTag > ("genParticlesInfo");

    //needed for trigger matching
    triggerEvent_ = (iConfig.getParameter< edm::InputTag > ("triggerEvent"));
    tauMatch_Loose_ = (iConfig.getParameter< std::string > ("tauMatch_Loose"));
    tauMatch_Medium_ = (iConfig.getParameter< std::string > ("tauMatch_Medium"));
    electronMatch_Loose_ = (iConfig.getParameter< std::string > ("electronMatch_Loose"));
    muonMatch_Loose_ = (iConfig.getParameter< std::string > ("muonMatch_Loose"));

    filterTriggerResults = iConfig.exists("filterTriggerResults") ? iConfig.getParameter<bool>("filterTriggerResults") : 0 ;
    puJetIdFlag_ = iConfig.getParameter<edm::InputTag>("puJetIdFlag");
    rhoProducer_ = iConfig.getParameter<edm::InputTag>("rhoProducer");

    el_trigger_name = iConfig.exists("el_trigger_name") ? iConfig.getParameter<std::string> ("el_trigger_name") : "Ele" ;
    mu_trigger_name = iConfig.exists("mu_trigger_name") ? iConfig.getParameter<std::string> ("mu_trigger_name") : "Mu" ;

    tauPtcut_ = iConfig.exists("tauPtcut_") ? iConfig.getParameter<double> ("tauPtCut") : 15.0;
}

NtupleProducer::~NtupleProducer() {

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}

void
NtupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


    if (Include_Vertex) DoVertexAnalysis(iEvent);
    if (Include_Electron) DoElectronAnalysis(iEvent, iSetup);
    if (Include_Muon) DoMuonAnalysis(iEvent, iSetup);
    if (Include_HPSTau) DoHPSTauAnalysis(iEvent, iSetup);
    if (Include_GenPartiles) DoGenParticlesAnalysis(iEvent);
    if (Include_Jet) DoJetAnalysis(iEvent);
    if (Include_JetCorrection) DoJetEC_Uncer_Analysis(iEvent);
    if (Include_MET) DoMetAnalysis(iEvent);
    if (Include_HLT) DoHLTAnalysis(iEvent);

    t->Fill();

}//analyze



//define this as a plug-in
//DEFINE_FWK_MODULE(NtupleProducer);


//DEFINE_SEAL_MODULE();
DEFINE_FWK_MODULE(NtupleProducer);
