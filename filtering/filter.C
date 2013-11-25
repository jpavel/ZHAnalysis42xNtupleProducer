#include "TChain.h"
#include "TFile.h"
#include "myevent.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "selection.h"

#include "LinkDef.h"


int filter(TString inputDir, long n_evt=100,bool is11=false,bool isFR=false)
{
	long TotalEvents=n_evt;
	
	//options
	bool is2011=is11;
    int N_ALL=4;
    int N_ELMU=-3;
    if(isFR){
		N_ALL=-4;
		N_ELMU=3;
	}
           
    double cut_mu_pt=9.0;
    double cut_mu_eta=2.4; 
    double cut_el_pt=9.0;
    double cut_el_eta=2.5;
    double cut_tau_pt=19.0;
    double cut_tau_eta=2.3;
    
    if(isFR) cut_tau_pt=0.0;
    
    double cut_dR=-0.09;
    double lepton_mass_min=50.0;
    double lepton_mass_max=150.0;
    
	std::string doubEle="HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
	std::string doubEle2="HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v";
	
    std::string doubMu="HLT_Mu17_Mu8_v";
    std::string doubMu2="HLT_Mu17_TkMu8_v"; // "HLT_Mu13_Mu8_v" in 2011
    if(is2011) doubMu2="HLT_Mu13_Mu8_v";
    std::string doubMu3="HLT_DoubleMu7_v";
    
    bool debug = false;
    

	
	TChain *chain = new TChain("t");
	chain->Add(inputDir+"*root*");
	
	
	
	//////////////////////// end of modular part ////////////////////////////////////////
	/////////////////////////COMMON part //////////////////////////////////////////
	myevent* in = new myevent();
	chain->SetBranchAddress("myevent",&in);
	
	long MaxEvents = chain->GetEntries();
	if(TotalEvents > MaxEvents || TotalEvents < 0) TotalEvents=MaxEvents;
	
	std::cout << "There are " << MaxEvents << " entries." << std::endl;
	
	
	//bookkeeping
	ofstream lumi;
	long current_run=-1;
	long current_lumi=-1;
    long m_allEvents=0;
	lumi.open("lumi.csv");
	
	//output definition
	
	 TFile *newFile = new TFile("output.root","recreate");
	TChain *newchain = (TChain*)chain->CloneTree(0);
	 TTree *tree = newchain->GetTree();
	
	
	for(uint i =0; i < TotalEvents; i++)
	{
		chain->GetEntry(i);
		m_allEvents++;
		
		if(debug) std::cout << "Processing entry " << i << " event: " << in->eventNumber << std::endl;
		if(i%1000==0)std::cout << "Analyzing entry no. " << i << std::endl;
		
		// storing lumis
		if(in->runNumber!=current_run || in->lumiNumber!=current_lumi){
			lumi << in->runNumber << " " << in->lumiNumber << std::endl;
			current_run=in->runNumber;
			current_lumi=in->lumiNumber;
		}
		
		if(PassFilter(in, debug, doubEle, doubEle2, doubMu, doubMu2, doubMu3,
			         is2011, N_ALL, N_ELMU, 
			         cut_mu_pt, cut_mu_eta, cut_el_pt, cut_el_eta, cut_tau_pt, cut_tau_eta,		    
			         cut_dR, lepton_mass_min, lepton_mass_max)
		   ) tree->Fill();
	}
	
	newFile->cd();
	newFile->Write();
	newFile->Close();
	
	 ofstream log1;       
     log1.open("total.txt");
     log1 << m_allEvents << std::endl;
     log1.close();
     
	return 0;
}
