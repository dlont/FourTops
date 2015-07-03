//////////////////////////////////////////////
//          Make cutflow table              //
//////////////////////////////////////////////

#include "../interface/CutsTable.h"

CutsTable::CutsTable(bool isMuon, bool isElectron){
	if (Muon){
		leptonChoice = "Muon";
	}
	else if (Electron){
		leptonChoice = "Electron";
	}
}

CutsTable::~CutsTable(){

}

void CutsTable::AddSelections(){
	CutsselecTable.push_back(string("initial"));
	CutsselecTable.push_back(string("Event cleaning and Trigger"));
	if(leptonChoice = "Muon"){
	    CutsselecTable.push_back(string("Exactly 1 Tight Isolated Muon"));
	    CutsselecTable.push_back(string("Exactly 1 Loose Isolated Muon"));
	    CutsselecTable.push_back(string("Exactly zero electrons"));

	}
	else if(leptonChoice = "Electron"){
	    CutsselecTable.push_back(string("Exactly 1 Tight Electron"));
	    CutsselecTable.push_back(string("Exactly 1 Loose Electron"));
	    CutsselecTable.push_back(string("Exactly zero muons"));
	}
	CutsselecTable.push_back(string("At least 6 Jets"));
	CutsselecTable.push_back(string("At least 1 CSVM Jet"));
	CutsselecTable.push_back(string("At least 2 CSVM Jets"));
	//CutsselecTable.push_back(string("HT $\\geq 100 GeV$"));
	//CutsselecTable.push_back(string("HT $\\geq 200 GeV$"));
	//CutsselecTable.push_back(string("HT $\\geq 300 GeV$"));
	//CutsselecTable.push_back(string("HT $\\geq 400 GeV$"));
}

void CutsTable::CreateTable(vector < Dataset* > datasets){

    SelectionTable selecTable(CutsselecTable, datasets);
    selecTable.SetLuminosity(Luminosity);
    selecTable.SetPrecision(1);

}

void CutsTable::FillTable(){
	//selectTable -> Fill..

}


void CutsTable::Calc_Write(){


}




