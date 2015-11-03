#include "../interface/Trigger.h"

Trigger::Trigger(bool isMuon, bool isElectron):
muon(false), electron(false), trigged(false), redotrigmap(false), triggerList(), currentRun(0), previousRun(-1), currentFilename(""),
 previousFilename(""), iFile(-1), triggermap(), runInfos2(new TRootRun()), previousDatasetName("")
{
	if(isMuon){
		muon = true;
	}
	else if(isElectron){
		electron = true;
	}
	else{
		cout<<"neither lepton selection"<<endl;
	}
}

Trigger::~Trigger(){

}

void Trigger::bookTriggers(){
	if(muon){
	    //triggerList.push_back("HLT_IsoMu20_eta2p1_v2");
	    triggerList.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v*");
	    triggerList.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v*");	
	    triggerList.push_back("HLT_IsoMu20_eta2p1_v*");
	    triggerList.push_back("HLT_IsoMu18_v*");
	    triggerList.push_back("HLT_IsoMu18_TriCentralPFJet50_40_30_v*");		    
	    triggerList.push_back("HLT_TkIsoMu20_eta2p1_v*");
   	}

    if (electron){
	    triggerList.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v*");
	    triggerList.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v*");
	    triggerList.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v*");
	    triggerList.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v*");
	    triggerList.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v*");
	    triggerList.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v*");    	
	    triggerList.push_back("HLT_Ele23_WPLoose_Gsf_v*"); 
	    triggerList.push_back("HLT_Ele23_WPLoose_Gsf_TriCentralPFJet50_40_30_v*");    	


    }

    for(UInt_t itrig=0; itrig<triggerList.size(); itrig++){
        triggermap[triggerList[itrig]]=std::pair<int,bool>(-999,false);
    }
 
}

void Trigger::checkAvail(int currentRun, vector < Dataset* > datasets, unsigned int d, TTreeLoader *treeLoader, TRootEvent* event){
	redotrigmap=false;
	currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	// currentRun = event->runId();
	if(previousFilename != currentFilename){
	    previousFilename = currentFilename;
	    iFile++;
	    redotrigmap=true;
	    cout<<"File changed!!! => iFile = "<<iFile << " new file is " << datasets[d]->eventTree()->GetFile()->GetName() << " in sample " << datasets[d]->Name() << endl;
	}
	if(previousRun != currentRun){
	    previousRun = currentRun;
	    cout<<"*****!!!!new run!!!! => new run = "<<previousRun<<" *****"<<endl;
	    redotrigmap=true;
	}
	// if(datasets[d]->Name() != previousDatasetName){
	//     datasets[d]->runTree()->SetBranchStatus("runInfos*",1);
	//     datasets[d]->runTree()->SetBranchAddress("runInfos",&runInfos2);		
	// }
	// datasets[d]->runTree()->GetEntry(iFile);

	if(redotrigmap){
		cout<<"redo map::print names list"<<endl;
		// runInfos2->getHLTinfo(currentRun).gethltNameList();
		//treeLoader->ListTriggers(currentRun, iFile);
	}


	// get trigger info:
	for(std::map<std::string,std::pair<int,bool> >::iterator iter = triggermap.begin(); iter != triggermap.end(); iter++){
	    if(redotrigmap){
			//  cout << "Bytes read out by GetEntry: " << rBytes << endl;
			//cout << "Getting HLT Path Info" << endl;
			// Int_t loc;
			// if(runInfos2 == 0){
			// 	loc = -9999;
			// }
			// else {
			// 	loc = runInfos2->getHLTinfo(currentRun).hltPath(iter->first);
			// }

	        Int_t loc = treeLoader->iTrigger(iter->first, currentRun, iFile);
	        string trigname = iter->first;
	        cout<<"trigname: "<<trigname<<"  location: "<<loc<<endl;

	        iter->second.first=loc;
	    }
	    // and check if it exists and if it fired:
	    if(iter->second.first>=0 && iter->second.first!=9999) // trigger exists
	        iter->second.second=event->trigHLT(iter->second.first);
	    else
	        iter->second.second=false;
	}   
}


int Trigger::checkIfFired(){
	// now check if the appropriate triggers fired for each analysis:
	trigged =0;

	for(UInt_t itrig=0; itrig<triggerList.size() && trigged==0; itrig++){
	    if(triggermap[triggerList[itrig]].second)   trigged=1;
	}

	return trigged;
}

