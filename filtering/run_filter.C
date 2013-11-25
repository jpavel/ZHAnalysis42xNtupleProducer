#include <iostream>
void run_filter()
{
	long n_evt=2000;
	TString input="/home/jpavel/ntuples/CMS/analysis/ZHtautau/HCP_sync/";
	bool is2011 = false;
	bool isFR = false;
	
	cout << "events: " << n_evt << endl;
	cout << "input directory: " << input << endl;
	if(is2011) cout << "Processing 2011 sample." << endl;
	else cout << "Processing 2012 sample." << endl;
	
	if(isFR) cout << "Filtering for fake rate estimation." << endl;
	else cout << "Filtering for analysis." << endl;
	
	gROOT->ProcessLine(" .L filter.C++ ");
	filter(input, n_evt, is2011, isFR);
	return;
}
