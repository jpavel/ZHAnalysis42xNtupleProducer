#ifndef helper_H
#define helper_H

#include "myevent.h"

#include "TMath.h"

#include <string>
#include <vector>
#include <cmath>

double deltaR(double eta1, double phi1, double eta2, double phi2){
	double dR, dPhi, dEta;
	dR=dPhi=dEta=0.0;
	dPhi = phi1-phi2;
	if (dPhi < -TMath::Pi()) dPhi += TMath::TwoPi();
	if (dPhi > +TMath::Pi()) dPhi -= TMath::TwoPi();
	dEta = eta1-eta2;
	dR = sqrt(dPhi*dPhi+dEta*dEta);

	return dR;
}

double deltaR(myobject o1, myobject o2)
{
	return deltaR(o1.eta,o1.phi,o2.eta,o2.phi);
}


bool Trg_MC_12(myevent* m,
	std::string doubEle,
	std::string doubEle2,
    std::string doubMu,
    std::string doubMu2,
    std::string doubMu3,
    bool is2011 )
{
	map<string, int> myHLT = m->HLT;
	bool Trigger = false;
	bool TriggerEle = false;
	bool TriggerMu = false;


	for (map<string, int> ::iterator ihlt = myHLT.begin(); ihlt != myHLT.end() && !TriggerEle && !TriggerMu; ihlt++) {
		//	std::cout << ihlt->first << std::endl; 
		size_t foundEl=(ihlt->first).find(doubEle);
		size_t foundEl2=(ihlt->first).find(doubEle2);
		if(!is2011) foundEl2=foundEl;

		size_t foundMu1=(ihlt->first).find(doubMu);
		size_t foundMu2=(ihlt->first).find(doubMu2);
		size_t foundMu3=(ihlt->first).find(doubMu3);
		if(!is2011) foundMu3=foundMu2;

		if (foundEl!=string::npos || foundEl2!=string::npos)
			TriggerEle = ihlt->second;
		if (foundMu1!=string::npos || foundMu2!=string::npos || foundMu3!=string::npos)
			TriggerMu = ihlt->second;
	}
	Trigger = TriggerEle || TriggerMu;
	return Trigger;    
}

#endif //helper_H
