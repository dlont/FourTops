#include "../interface/Trigger.h"

Trigger::Trigger(bool isMuon, bool isElectron):
muon(false), electron(false), trigged(false), redotrigmap(false), triggerListDataC(), triggerListDataD(), triggerListMC(), triggerList(), currentRun(0), previousRun(-1), currentFilename(""),
 previousFilename(""), iFile(-1), triggermapDataC(), triggermapDataD(), triggermapMC(), triggermap(), runInfos2(new TRootRun()), previousDatasetName("")
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
	    // triggerListDataC.push_back("HLT_IsoMu20_eta2p1_v*");	    
	    // triggerListDataC.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v*");
	    // triggerListDataC.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v*");	

	    // triggerListDataD.push_back("HLT_IsoMu18_v*");
	    // triggerListDataD.push_back("HLT_IsoMu18_TriCentralPFJet50_40_30_v*");	

	    triggerListDataC.push_back("HLT_IsoTkMu20_v*");
   	    triggerListDataC.push_back("HLT_IsoMu20_v*");
	    triggerListDataD.push_back("HLT_IsoTkMu20_v*");
   	    triggerListDataD.push_back("HLT_IsoMu20_v*");

	    triggerListMC.push_back("HLT_IsoTkMu20_v*");
   	    triggerListMC.push_back("HLT_IsoMu20_v*");
	    // triggerListMC.push_back("HLT_IsoMu17_eta2p1_v*");
   	    // triggerListMC.push_back("HLT_IsoMu20_eta2p1_v*");
   	    // triggerListMC.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v*");
   	    // triggerListMC.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v*");
	    // triggerList.push_back("HLT_TkIsoMu20_eta2p1_v*");
	    // triggerList.push_back("HLT_IsoMu20_eta2p1_v2");
	    // triggerList.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v2");
	    // triggerList.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v2");	
	    // triggerList.push_back("HLT_IsoMu20_eta2p1_v1");
	    // triggerList.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet30_v1");
	    // triggerList.push_back("HLT_IsoMu20_eta2p1_TriCentralPFJet50_40_30_v1");		    
	    // triggerList.push_back("HLT_TkIsoMu20_eta2p1_v1");	    

	}

    if (electron){
	    triggerListDataC.push_back("HLT_Ele23_WPLoose_Gsf_v*");    	
	    triggerListDataC.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v*");
	    triggerListDataC.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v*");
	    triggerListDataC.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v*");

	    triggerListDataD.push_back("HLT_Ele23_WPLoose_Gsf_v*"); 
	    triggerListDataD.push_back("HLT_Ele23_WPLoose_Gsf_TriCentralPFJet50_40_30_v*");    	


	    triggerListMC.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v*");
	    triggerListMC.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v*");
	    triggerListMC.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v*");    

	    // triggerList.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v1");
	    // triggerList.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet30_v1");
	    // triggerList.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_TriCentralPFJet50_40_30_v1");
	    // triggerList.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v1");
	    // triggerList.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet30_v1");
	    // triggerList.push_back("HLT_Ele27_eta2p1_WP75_Gsf_TriCentralPFJet50_40_30_v1");  



    }

    for(UInt_t itrig=0; itrig<triggerListMC.size(); itrig++){
        // triggermapDataC[triggerListDataC[itrig]]=std::pair<int,bool>(-999,false);
        // triggermapDataD[triggerListDataC[itrig]]=std::pair<int,bool>(-999,false);
        triggermapMC[triggerListMC[itrig]]=std::pair<int,bool>(-999,false);

        // triggermap[triggerList[itrig]]=std::pair<int,bool>(-999,false);
    }

        for(UInt_t itrig=0; itrig<triggerListDataC.size(); itrig++){
        // triggermapDataC[triggerListDataC[itrig]]=std::pair<int,bool>(-999,false);
        // triggermapDataD[triggerListDataC[itrig]]=std::pair<int,bool>(-999,false);
        triggermapDataC[triggerListDataC[itrig]]=std::pair<int,bool>(-999,false);

        // triggermap[triggerList[itrig]]=std::pair<int,bool>(-999,false);
    }

    for(UInt_t itrig=0; itrig<triggerListDataD.size(); itrig++){
        // triggermapDataC[triggerListDataC[itrig]]=std::pair<int,bool>(-999,false);
        // triggermapDataD[triggerListDataC[itrig]]=std::pair<int,bool>(-999,false);
        triggermapDataD[triggerListDataD[itrig]]=std::pair<int,bool>(-999,false);

        // triggermap[triggerList[itrig]]=std::pair<int,bool>(-999,false);
    }
 
}

void Trigger::checkAvail(int currentRun, vector < Dataset* > datasets, unsigned int d, TTreeLoader *treeLoader, TRootEvent* event, int treeNumber){
	redotrigmap=false;
	currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
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

	iFile = treeNumber;
	// if(redotrigmap){
	// 	treeLoader->ListTriggers(currentRun, iFile);
	// }


	// get trigger info:
	string datasetName = datasets[d]->Name();
	// bool bdizzle = datasetName.find("Data")==string::npos;
	if(datasetName.find("Data")==string::npos){
	
		for(std::map<std::string,std::pair<int,bool> >::iterator iter = triggermapMC.begin(); iter != triggermapMC.end(); iter++){
		    if(redotrigmap){
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
	else if(currentRun>253658 && currentRun<256465){
		for(std::map<std::string,std::pair<int,bool> >::iterator iter = triggermapDataC.begin(); iter != triggermapDataC.end(); iter++){
		    if(redotrigmap){
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
	else if(currentRun>256629){
		for(std::map<std::string,std::pair<int,bool> >::iterator iter = triggermapDataD.begin(); iter != triggermapDataD.end(); iter++){
		    if(redotrigmap){
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

}


int Trigger::checkIfFired(int currentRun, vector < Dataset* > datasets, unsigned int d){
	// now check if the appropriate triggers fired for each analysis:
	trigged =0;
	string datasetName = datasets[d]->Name();

	if(datasetName.find("Data")==string::npos){
		for(UInt_t itrig=0; itrig<triggerListMC.size() && trigged==0; itrig++){
			// cout<<"fired: "<<triggermapMC[triggerListMC[itrig]].second<<endl;
		    if(triggermapMC[triggerListMC[itrig]].second)   trigged=1;
		}
	}
	else if(currentRun>253658 && currentRun<256465){
		for(UInt_t itrig=0; itrig<triggerListDataC.size() && trigged==0; itrig++){
			// cout<<"fired: "<<triggermapDataC[triggerListDataC[itrig]].second<<endl;
		    if(triggermapDataC[triggerListDataC[itrig]].second)   trigged=1;
		}		
	}
	else if(currentRun>256629){
		for(UInt_t itrig=0; itrig<triggerListDataD.size() && trigged==0; itrig++){
			// cout<<"fired: "<<triggermapDataD[triggerListDataD[itrig]].second<<endl;
		    if(triggermapDataD[triggerListDataD[itrig]].second)   trigged=1;
		}		
	}
	else{
		cout<<"currentRun::: "<<currentRun<<endl;
	}
	return trigged;
}

