#include "Analysis/NtupleProducer/interface/NtupleProducer.h"

void NtupleProducer::DoJetEC_Uncer_Analysis(const edm::Event& iEvent) {

    (m->RecPFJetsAK5_Up_SubTotal).clear();
    (m->RecPFJetsAK5_Down_SubTotal).clear();
    (m->RecPFJetsAK5_Up_Total).clear();
    (m->RecPFJetsAK5_Down_Total).clear();

    //=================== PF Jet AK5 =================================
    using pat::Jet;
    using pat::JetCollection;

    Handle<pat::JetCollection> PFjetsAK5_Up_SubTotal;
    iEvent.getByLabel("jetsEnUpSubTotalDataMC", PFjetsAK5_Up_SubTotal);

    Handle<pat::JetCollection> PFjetsAK5_Down_SubTotal;
    iEvent.getByLabel("jetsEnDownSubTotalDataMC", PFjetsAK5_Down_SubTotal);

    Handle<pat::JetCollection> PFjetsAK5_Up_Total;
    iEvent.getByLabel("jetsEnUpTotal", PFjetsAK5_Up_Total);

    Handle<pat::JetCollection> PFjetsAK5_Down_Total;
    iEvent.getByLabel("jetsEnDownTotal", PFjetsAK5_Down_Total);



    //=================== PF Jet AK5 =================================
    for (JetCollection::const_iterator jet = PFjetsAK5_Up_SubTotal->begin(); jet != PFjetsAK5_Up_SubTotal->end(); jet++) {
        myobject myjet;
        myjet.pt = jet->pt();
        myjet.eta = jet->eta();
        myjet.phi = jet->phi();
        myjet.et = jet->et();


        if (jet->pt() > 10) {
            (m->RecPFJetsAK5_Up_SubTotal).push_back(myjet);
        }
    }//loop over jets

    for (JetCollection::const_iterator jet = PFjetsAK5_Down_SubTotal->begin(); jet != PFjetsAK5_Down_SubTotal->end(); jet++) {
        myobject myjet;
        myjet.pt = jet->pt();
        myjet.eta = jet->eta();
        myjet.phi = jet->phi();
        myjet.et = jet->et();


        if (jet->pt() > 10) {
            (m->RecPFJetsAK5_Down_SubTotal).push_back(myjet);
        }
    }//loop over jets
    for (JetCollection::const_iterator jet = PFjetsAK5_Up_Total->begin(); jet != PFjetsAK5_Up_Total->end(); jet++) {
        myobject myjet;
        myjet.pt = jet->pt();
        myjet.eta = jet->eta();
        myjet.phi = jet->phi();
        myjet.et = jet->et();


        if (jet->pt() > 10) {
            (m->RecPFJetsAK5_Up_Total).push_back(myjet);
        }
    }//loop over jets

    for (JetCollection::const_iterator jet = PFjetsAK5_Down_Total->begin(); jet != PFjetsAK5_Down_Total->end(); jet++) {
        myobject myjet;
        myjet.pt = jet->pt();
        myjet.eta = jet->eta();
        myjet.phi = jet->phi();
        myjet.et = jet->et();


        if (jet->pt() > 10) {
            (m->RecPFJetsAK5_Down_Total).push_back(myjet);
        }
    }//loop over jets



}

