#include <FWCore/Framework/interface/Event.h>

#include "Analysis/NtupleProducer/interface/NtupleProducer.h"

void NtupleProducer::DoMetAnalysis(const edm::Event& iEvent) {
    //    (m->RecMet).clear();
    (m->RecPFMet).clear();
    (m->RectcMet).clear();
    (m->RecPFMetCor).clear();
    (m->RecMVAMet).clear();

    //    Handle<CaloMETCollection> met;
    //    iEvent.getByLabel(MetCollection_, met);

    Handle<METCollection> tcmet;
    iEvent.getByLabel(tcMetCollection_, tcmet);

    Handle<PFMETCollection> PFmet;
    iEvent.getByLabel(PFMetCollection_, PFmet);

    Handle<PFMETCollection> MVAmet;
    iEvent.getByLabel(MVAMetCollection_, MVAmet); 

    Handle< edm::View<reco::PFMET> > mvaMEThandle;
    iEvent.getByLabel(MVAMetCollection_, mvaMEThandle);


    Handle< edm::View<reco::PFMET> > pfMEThandle;
    iEvent.getByLabel(PFMetCollection_, pfMEThandle);


    //    for (reco::CaloMETCollection::const_iterator met_iter = met->begin();
    //            met_iter != met->end(); met_iter++) {
    //
    //        myobject mymet;
    //        mymet.pt = met_iter->pt();
    ////        mymet.px = met_iter->px();
    ////        mymet.py = met_iter->py();
    //        mymet.phi = met_iter->phi();
    //        mymet.et = met_iter->et();
    //        (m->RecMet).push_back(mymet);
    //
    //    }

    //============== PF Met ==================================

    for (reco::PFMETCollection::const_iterator met_iter = PFmet->begin();
            met_iter != PFmet->end(); met_iter++) {

        myobject mymet;
        mymet.pt = met_iter->pt();
        mymet.px = met_iter->px();
        mymet.py = met_iter->py();
        mymet.phi = met_iter->phi();
        mymet.et = met_iter->et();
        (m->RecPFMet).push_back(mymet);

    }

    // MET significance for SV fit calculation

    TMatrixD covMET(2, 2);
    covMET = (pfMEThandle->front()).getSignificanceMatrix();
    m->MET_sigMatrix_00 = covMET[0][0];
    m->MET_sigMatrix_01 = covMET[0][1];
    m->MET_sigMatrix_10 = covMET[1][0];
    m->MET_sigMatrix_11 = covMET[1][1];


    // MVA MET

    for (reco::PFMETCollection::const_iterator met_iter = MVAmet->begin();
	 met_iter != MVAmet->end(); met_iter++) {
      myobject mymet;
      mymet.pt = met_iter->pt();
      mymet.px = met_iter->px();
      mymet.py = met_iter->py();
      mymet.phi = met_iter->phi();
      mymet.et = met_iter->et();
      (m->RecMVAMet).push_back(mymet);

    }

    TMatrixD covMVAMet(2, 2);
    covMVAMet = (mvaMEThandle->front()).getSignificanceMatrix();
    m->MVAMet_sigMatrix_00 = covMVAMet[0][0];
    m->MVAMet_sigMatrix_01 = covMVAMet[0][1];
    m->MVAMet_sigMatrix_10 = covMVAMet[1][0];
    m->MVAMet_sigMatrix_11 = covMVAMet[1][1];


    //============== tc Met ==================================

    for (reco::METCollection::const_iterator met_iter = tcmet->begin();
            met_iter != tcmet->end(); met_iter++) {

        myobject mymet;
        mymet.pt = met_iter->pt();
        mymet.px = met_iter->px();
        mymet.py = met_iter->py();
        mymet.phi = met_iter->phi();
        mymet.et = met_iter->et();
        (m->RectcMet).push_back(mymet);

    }
    //============== Type1Corrected MET ==================================
//     typedef edm::View<reco::MET> METView;

//     Handle<METView> Type1CorMEThandle;
//     iEvent.getByLabel(Type1CorMETCollection_, Type1CorMEThandle);
//     for (METView::const_iterator met_iter = Type1CorMEThandle->begin();
//             met_iter != Type1CorMEThandle->end(); met_iter++) {
//         myobject mymet;
//         mymet.pt = met_iter->pt();
//         mymet.px = met_iter->px();
//         mymet.py = met_iter->py();
//         mymet.phi = met_iter->phi();
//         mymet.et = met_iter->et();
//         (m->RecPFMetCor).push_back(mymet);
//     }

}
