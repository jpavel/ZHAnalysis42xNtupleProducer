#include <FWCore/Framework/interface/Event.h>
#include "Analysis/NtupleProducer/interface/NtupleProducer.h"

void NtupleProducer::DoHLTAnalysis(const edm::Event& iEvent) {
    (m->HLT).clear();
    (m->HLT_DiElectron) = 0;
    (m->HLT_DiMuon) = 0;

    edm::Handle<TriggerResults> triggerResults;
    iEvent.getByLabel(srcTriggerResults_, triggerResults);
    
    
    string eleTrigger (el_trigger_name);
    string muTrigger (mu_trigger_name);


    if (triggerResults.isValid()) {
        int ntrigs = triggerResults->size();
        TriggerNames const &triggerNames = iEvent.triggerNames(*triggerResults);

        for (int itrig = 0; itrig < ntrigs; itrig++) {
            string name = triggerNames.triggerName(itrig);
            bool result = triggerResults->accept(itrig);
	    size_t foundEl=name.find(eleTrigger);
	    size_t foundMu=name.find(muTrigger);
	    if(filterTriggerResults && (foundEl!=string::npos || foundMu!=string::npos))
             (m->HLT)[name] = result;
	    else if(!filterTriggerResults)
	      (m->HLT)[name] = result;


            if (name == "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") {
                m->HLT_DiElectron = result;
                //                cout << "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2" << endl;
            }
            if (name == "HLT_DoubleMu7_v1") {
                m->HLT_DiMuon = result;
                //                cout << "HLT_DoubleMu7_v1" << endl;
            }


        }//for itrig
    }//if triggerResults valid



//
//
//edm::Handle<pat::TriggerPathCollection> paths;
//    if(iEvent.getByLabel(src_, paths)) {
//
//      //get the names of the triggers
//      for(unsigned int i=0;i<paths_.size();++i) {
//	bool found=false;
//		bool fired_t=false;
//	for(unsigned int j=0;j<paths->size() && fired_t==false;++j) {
//          size_t trigPath = paths->at(j).name().find(paths_.at(i));
//          if ( trigPath == 0) {
//	    found=true;
//	    fired[i]=paths->at(j).wasAccept();
//	    wasRun[i]=paths->at(j).wasRun();
//	    prescale[i]=paths->at(j).prescale();
//	    error[i]=paths->at(j).wasError();
//	    if(paths->at(j).wasRun()&&paths->at(j).wasAccept()&&paths->at(j).prescale()==1){
//	      any=1; fired_t=true;
//		}
//	    break;
//	  }
//	}
//	if(!found) {
//	  fired[i]=0;
//	  wasRun[i]=0;
//	  prescale[i]=0;
//	  error[i]=0;
//
//	}
//      }
//
//    }






}
