#include "Analysis/NtupleProducer/interface/NtupleProducer.h"

void NtupleProducer::DoGenParticlesAnalysis(const edm::Event& iEvent) {

    (m->RecGenParticle).clear();
    (m->RecGenMet).clear();
    (m->RecGenJet).clear();
    (m->RecGenTauVisible).clear();
    if (IsMC) {

        //=================== PF gen AK5 =================================
        using reco::GenParticle;
        using reco::GenMET;
        using reco::GenJet;

        Handle<reco::GenParticleCollection> GenPar;
        iEvent.getByLabel(GenParticlesInfo_, GenPar);
	int index = 0;
        for (reco::GenParticleCollection::const_iterator gen = GenPar->begin(); gen != GenPar->end(); gen++) {

            myGenobject mygen;
	    
            mygen.pdgId = gen->pdgId();
            mygen.status = gen->status();
            mygen.pt = gen->pt();
            mygen.eta = gen->eta();
            mygen.phi = gen->phi();
            mygen.charge = gen->charge();
            mygen.z = gen->vz();
	    mygen.mass = gen->mass();

            mygen.mod_pdgId = (gen->mother() != NULL ? gen->mother()->pdgId() : -1000);
            mygen.mod_status = (gen->mother() != NULL ? gen->mother()->status() : -1000);
            mygen.mod_pt = (gen->mother() != NULL ? gen->mother()->pt() : -1000);
            mygen.mod_eta = (gen->mother() != NULL ? gen->mother()->eta() : -1000);
            mygen.mod_phi = (gen->mother() != NULL ? gen->mother()->phi() : -1000);
            mygen.mod_charge = (gen->mother() != NULL ? gen->mother()->charge() : -1000);
            mygen.mod_z = (gen->mother() != NULL ? gen->mother()->vz() : -1000);
	    mygen.mod_mass = (gen->mother() != NULL ? gen->mother()->mass() : -1000);

            mygen.Gmod_pdgId = ((gen->mother() != NULL && gen->mother()->mother() != NULL) ? gen->mother()->mother()->pdgId() : -1000);
            mygen.Gmod_status = ((gen->mother() != NULL && gen->mother()->mother() != NULL) ? gen->mother()->mother()->status() : -1000);
            mygen.Gmod_pt = ((gen->mother() != NULL && gen->mother()->mother() != NULL) ? gen->mother()->mother()->pt() : -1000);
            mygen.Gmod_eta = ((gen->mother() != NULL && gen->mother()->mother() != NULL) ? gen->mother()->mother()->eta() : -1000);
            mygen.Gmod_phi = ((gen->mother() != NULL && gen->mother()->mother() != NULL) ? gen->mother()->mother()->phi() : -1000);
            mygen.Gmod_charge = ((gen->mother() != NULL && gen->mother()->mother() != NULL) ? gen->mother()->mother()->charge() : -1000);
            mygen.Gmod_z = ((gen->mother() != NULL && gen->mother()->mother() != NULL) ? gen->mother()->mother()->vz() : -1000);
	    mygen.Gmod_mass = ((gen->mother() != NULL && gen->mother()->mother() != NULL) ? gen->mother()->mother()->mass() : -1000);

            (m->RecGenParticle).push_back(mygen);
	    // storing visible tau 4-vectors
	    if(abs(gen->pdgId())==15)
	      {
		const reco::GenParticleRefVector& mRefs = gen->daughterRefVector();
		reco::Particle::LorentzVector invisibleP4( 0.0, 0.0, 0.0, 0.0 );
		int decayMode = -1; // 0= hadronic, 1=electron, 2=muon
		unsigned int nNu = 0;
		for(reco::GenParticleRefVector::const_iterator imr = mRefs.begin(); imr!= mRefs.end(); ++imr) {
		  if(abs((*imr)->pdgId())==11) decayMode = 1;
		  if(abs((*imr)->pdgId())==13) decayMode = 2;
		  if(abs((*imr)->pdgId())==16 || abs((*imr)->pdgId())==14 || abs((*imr)->pdgId())==12)
		    {
		      invisibleP4+=(*imr)->p4();
		      nNu++;
		    }
		}
		if(nNu==1) decayMode = 0;
		reco::Particle::LorentzVector visibleP4 = gen->p4()-invisibleP4;
		myGenobject myVisTau;
		myVisTau.pdgId=gen->pdgId();
		myVisTau.status = gen->status();
		myVisTau.gen_index = index;
		myVisTau.pt = visibleP4.pt();
		myVisTau.eta = visibleP4.eta();
		myVisTau.phi = visibleP4.phi();
		myVisTau.charge = gen->charge();
		myVisTau.mass = visibleP4.M();
		myVisTau.decay_mode = decayMode;
		(m->RecGenTauVisible).push_back(myVisTau);
	      }
	    index++;
        }//loop over gens





        Handle<reco::GenMETCollection> genMet_;
        iEvent.getByLabel("genMetTrue", genMet_);
        for (reco::GenMETCollection::const_iterator genmet = genMet_->begin(); genmet != genMet_->end(); genmet++) {


            myGenobject myGenMet;
            myGenMet.pt = genmet->pt();
            myGenMet.phi = genmet->phi();
            myGenMet.et = genmet->et();
            (m->RecGenMet).push_back(myGenMet);


        }


        Handle<reco::GenJetCollection> genJet_;
        iEvent.getByLabel("ak5GenJets", genJet_);
        for (reco::GenJetCollection::const_iterator genjet = genJet_->begin(); genjet != genJet_->end(); genjet++) {

            myGenobject myGenJet;
            myGenJet.pt = genjet->pt();
            myGenJet.eta = genjet->eta();
            myGenJet.phi = genjet->phi();
            myGenJet.et = genjet->et();
            (m->RecGenJet).push_back(myGenJet);

        }
    }
}

