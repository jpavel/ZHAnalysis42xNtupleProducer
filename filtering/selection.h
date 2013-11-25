#ifndef selection_H
#define selection_H

#include "myevent.h"
#include "helper.h"

#include "TLorentzVector.h"

bool PassFilter(myevent* m, bool debug,
				std::string doubEle,
				std::string doubEle2,
				
			    std::string doubMu,
			    std::string doubMu2,
			    std::string doubMu3,
			    
			    bool is2011,
			    int N_ALL,
			    int N_ELMU,
			    
			    double cut_mu_pt,
			    double cut_mu_eta,    
			    double cut_el_pt,
			    double cut_el_eta,
			    double cut_tau_pt,
			    double cut_tau_eta,
			    
			    double cut_dR,
			    double lepton_mass_min,
			    double lepton_mass_max)
{
	
	bool printout = debug;
	if(printout) std::cout << "Here!" << std::endl;
	
	bool trigPass = Trg_MC_12(m,doubEle, doubEle2, doubMu, doubMu2, doubMu3, is2011);
	if(printout) std::cout << " Trigger decision " << trigPass << std::endl;
	if(!trigPass)
	{
		if(printout) std::cout << "Trigger!" << std::endl;
		return false;
	}
	if(printout) std::cout << "Here2!" << std::endl;
	
	std::vector<myobject> vertex = m->Vertex;
	bool goodVertex = false;
	
	for (uint i = 0; i < vertex.size() && !goodVertex; i++) {
		if (vertex[i].isValid && vertex[i].normalizedChi2 > 0 && vertex[i].ndof > 4 && fabs(vertex[i].z) < 24)
		goodVertex=true;
	}
	
	if(!goodVertex){		
		if(printout) std::cout << "vertex!" << std::endl;
		return false;
	}
	
	
	std::vector<myobject> muon = m->PreSelectedMuons;
	if(printout) std::cout << " There are " << muon.size() << " preselected muons " << std::endl;

	std::vector<myobject> goodMuon;
	goodMuon.clear();

	for (uint i = 0; i < muon.size(); i++) {

		double muPt = muon[i].pt;
		double muEta = muon[i].eta;	
		if (muPt > cut_mu_pt && fabs(muEta) < cut_mu_eta)
		{
				goodMuon.push_back(muon[i]);
		}

	}
	
	uint N_MU=goodMuon.size();
	
	std::vector<myobject> electron = m->PreSelectedElectrons;
	if(printout) std::cout << " There are " << electron.size() << " preselected electrons " << std::endl;
	
	std::vector<myobject> goodElectron;
	goodElectron.clear();

	for (uint i = 0; i < electron.size(); i++) {

		double elPt = electron[i].pt;
		double elEta = electron[i].eta;
		
		if (elPt > cut_el_pt && fabs(elEta) < cut_el_eta )
		{
			goodElectron.push_back(electron[i]);
		}
	}
	
	uint N_EL=goodElectron.size();
	
	std::vector<myobject> tau = m->PreSelectedHPSTaus;

	if(printout)  std::cout << " There are " << tau.size() << " preselected taus " << std::endl;

	std::vector<myobject> goodTau;
	goodTau.clear();

	for (uint i = 0; i < tau.size(); i++) {

		double tauPt = tau[i].pt;
		double tauEta = tau[i].eta;
		bool DecayMode = (tau[i].discriminationByDecayModeFinding > 0.5);


		if (tauPt > cut_tau_pt && fabs(tauEta) < cut_tau_eta && DecayMode)
			goodTau.push_back(tau[i]);
	}
	
	uint N_TAU=goodTau.size();
	if(printout)  std::cout << " There are " << N_MU << " muons " << N_EL << " electrons and " << N_TAU << " taus." << std::endl;

	if(N_ALL > 0 && (N_MU+N_EL+N_TAU) < N_ALL){		if(printout) std::cout << "N_ALL!" << std::endl; return false;}
	if(N_ELMU > 0 && (N_MU+N_EL) < N_ELMU){		if(printout) std::cout << "N_ELMU!" << std::endl; return false;}
	
	if(printout) std::cout << "Here3!" << std::endl;
	
	bool failDR = false;
	if(cut_dR > 0){
		if(N_MU > 0)
		{
			for(uint iMu=0; iMu < N_MU && !failDR; iMu++)
			{
				for(uint jMu =iMu+1; jMu < N_MU && !failDR; jMu++)
				{
					double dR=deltaR(goodMuon[iMu],goodMuon[jMu]);
					if(printout) std::cout << " MuMu: " << iMu << " " << jMu << " :" << dR << std::endl;	
					if(dR < cut_dR) failDR=true;
				}
				
				for(uint iEl =0; iEl < N_EL && !failDR; iEl++)
				{
					double dR=deltaR(goodMuon[iMu],goodElectron[iEl]);
					if(printout) std::cout << " MuEl: " << iMu << " " << iEl << " :" << dR << std::endl;	
					if(dR < cut_dR) failDR=true;
				}
				
				for(uint iTau =0; iTau < N_TAU && !failDR; iTau++)
				{
					double dR=deltaR(goodMuon[iMu],goodTau[iTau]);
					if(printout) std::cout << " MuTau: " << iMu << " " << iTau << " :" << dR << std::endl;	
					if(dR < cut_dR) failDR=true;
				}
			}
		}else if(N_EL > 0)
		{
			for(uint iEl=0; iEl < N_EL && !failDR; iEl++)
			{
				for(uint jEl =iEl+1; jEl < N_EL && !failDR; jEl++)
				{
					double dR=deltaR(goodElectron[iEl],goodElectron[jEl]);
					if(printout) std::cout << " ElEl: " << iEl << " " << jEl << " :" << dR << std::endl;	
					if(dR < cut_dR) failDR=true;
				}
				
				for(uint iTau =0; iTau < N_TAU && !failDR; iTau++)
				{
					double dR=deltaR(goodElectron[iEl],goodTau[iTau]);
					if(printout) std::cout << " ElTau: " << iEl << " " << iTau << " :" << dR << std::endl;	
					if(dR < cut_dR) failDR=true;
				}
			}
		
		}else if(N_TAU > 0)
		{
			for(uint iTau=0; iTau < N_TAU && !failDR; iTau++)
			{
				for(uint jTau =iTau+1; jTau < N_TAU && !failDR; jTau++)
				{
					double dR=deltaR(goodTau[iTau],goodTau[jTau]);
					if(printout) std::cout << " TauTau: " << iTau << " " << jTau << " :" << dR << std::endl;	
					if(dR < cut_dR) failDR=true;
				}
			}
		}
	}
	
	if(failDR){		if(printout) std::cout << "dR" << std::endl; return false;}
	if(printout) std::cout << "Here4!" << std::endl;
	
	
	// dMass cut
	
	bool passMass = false;
	
	for (uint i = 0; i < goodMuon.size() && !passMass; i++) {
		for (uint j = i+1; j < goodMuon.size() && !passMass; j++) {
			TLorentzVector cand;
			cand.SetPxPyPzE(goodMuon[i].px+goodMuon[j].px,
					goodMuon[i].py+goodMuon[j].py,
					goodMuon[i].pz+goodMuon[j].pz,
					goodMuon[i].E+goodMuon[j].E);
			double mass = cand.M();
			if(printout) std::cout << mass << std::endl;
			if(mass > lepton_mass_min && mass < lepton_mass_max ) passMass=true;
			if(printout) {
				if(passMass) std::cout << " Passed mass cut with value " << mass << std::endl;
				else std::cout << " Failed mass cut with value " << mass << std::endl;
			}
		}
	}
	
	for (uint i = 0; i < goodElectron.size() && !passMass; i++) {
		for (uint j = i+1; j < goodElectron.size() && !passMass; j++) {
			TLorentzVector cand;
			cand.SetPxPyPzE(goodElectron[i].px+goodElectron[j].px,
					goodElectron[i].py+goodElectron[j].py,
					goodElectron[i].pz+goodElectron[j].pz,
					goodElectron[i].E+goodElectron[j].E);
			double mass = cand.M();
			if(printout) std::cout << mass << std::endl;
			if(mass > lepton_mass_min && mass < lepton_mass_max ) passMass=true;
			if(printout) {
				if(passMass) std::cout << " Passed mass cut with value " << mass << std::endl;
				else std::cout << " Failed mass cut with value " << mass << std::endl;
			}
		}
	}
	
	if(!passMass){		if(printout) std::cout << "Mass" << std::endl; return false;}
	if(printout) std::cout << "Here5!" << std::endl;
	
	if(printout) std::cout << " Passed!" << std::endl;
	
	
	return true;
}

#endif // selection_H

