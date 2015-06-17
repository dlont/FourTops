
#define _USE_MATH_DEFINES

#include "TTree.h"
#include "TNtuple.h"
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <errno.h>
#include <iostream>
#include <cstdlib>
#include "TChain.h"

using namespace std;


int main(){

	cout<<"start"<<endl;
	TChain ch("EventTree");
	ch.Add("/pnfs/iihe/cms/store/user/heilman/FOURTOP_NLO_13TeV_aMCatNLO_LHE_pythia8/crab_TOPTREE_742_FOURTOP_1M_v3/150605_093838/0000/*");
	cout<<"added files"<<endl;
	int nEntries = ch.GetEntries();
	cout<<nEntries<<endl;

	return 0;
}

