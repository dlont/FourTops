#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"
#include "TNtuple.h"

//user code
#include "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/TopTreeProducer/interface/TRootRun.h"
#include "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/TopTreeProducer/interface/TRootEvent.h"
#include "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
//#include "../macros/Style.C"

using namespace std;
using namespace TopTree;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,MultiSamplePlot*> MSPlot;

void SystematicsAnalyser(int nBins, string lepton, string leptoAbbr, bool Normalise, TFile *shapefile, TFile *errorfile);
void DatasetPlotter(int nBins, string lepton, string leptoAbbr, TFile *shapefile, TFile *errorfile);


int main()
{
	int NumberOfBins = 10;	//fixed width nBins
	TFile *shapefile = new TFile("shapefile.root", "RECREATE");
	TFile *errorfile = new TFile("caleFilesMu_light/Error.root","RECREATE");
	for(int num=0; num<1; num++)
	{
		string lepton;
		string leptoAbbr;
		if(num==0)
		{
			lepton = "Muon";
			leptoAbbr = "Mu";
		}
		else
		{
			lepton = "Electron";
			leptoAbbr = "El";
		}

		SystematicsAnalyser(NumberOfBins, lepton, leptoAbbr, false, shapefile, errorfile);
		DatasetPlotter(NumberOfBins, lepton, leptoAbbr, shapefile, errorfile);
	}
	delete shapefile;
	delete errorfile;
}


void DatasetPlotter(int nBins, string lepton, string leptoAbbr, TFile *shapefile, TFile *errorfile)
{
	string xmlFileName = "";  xmlFileName = "config/Run2SingleLepton_samples.xml";		//xmlFileName += lepton;		xmlFileName += "Full.xml";
	const char *xmlfile = xmlFileName.c_str();	cout << "used config file: " << xmlfile << endl; 

	string pathPNG = "FourTop_Light"; 	pathPNG += leptoAbbr;//name for output directory and root plots
	//TFile *fout = new TFile ((pathPNG+".root").c_str(), "RECREATE");   
	//TFile *fnom = new TFile("LimitSettingLight/NominalShapes_Mu_light.root","RECREATE"); //file for writing

	pathPNG += "_MSPlots/"; 	//add MSplot name to directory
	mkdir(pathPNG.c_str(),0777);  	cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

	///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
    //cout<<"loading...."<<endl;
	TTreeLoader treeLoader;
	vector < Dataset* > datasets; cout<<"vector filled"<<endl;
	treeLoader.LoadDatasets (datasets, xmlfile);	cout<<"datasets loaded"<<endl;

	//***************************************************CREATING PLOTS****************************************************
	////////////////////for i bins////////////////
	//string ni = static_cast<ostringstream*>( &(ostringstream() << nBins) )->str();
	string plainplotname = "HT";   ///// Non Jet Split plot
	//string nk;
	MSPlot[plainplotname.c_str()] = new MultiSamplePlot(datasets, plainplotname.c_str(), nBins, 0, 1000, "BDT Discriminator");

	//***********************************************OPEN FILES & GET NTUPLES**********************************************
	string dataSetName, filename, filepath;
	int nEntries; 
	float MVA, ScaleFactor, NormFactor, Luminosity, nJets, Ntags, HTb, HTX, HTH, HT;
	//cout<<"plotname  "<<plainplotname<<endl;
	cout<<"dataset size: "<<datasets.size()<<endl;

	for (int d = 0; d < datasets.size(); d++){ //Loop through datasets
		dataSetName = datasets[d]->Name();	cout<<"Dataset:  :"<<dataSetName<<endl;
		
		if(lepton == "Muon")//muon
		{
			filepath = "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/FourTop/Craneens_Mu/Craneens5_2_2015/Craneen_"+dataSetName + "_Run2_TopTree_Study_"+dataSetName+".root";
		}
		else if (lepton == "Electron") 
		{
			//filepath = "/Users/lb8075/JamesCode/TopBrussels/LightTrees_El/140206/"+ dataSetName + "_El.root";
		}
		else{cout<<"lepton error"<<endl;}
		filename = dataSetName;
		cout<<"filepath: "<<filepath<<endl; 
		cout<<"filename"<<filename<<endl;
		FileObj[filename.c_str()] = new TFile((filepath).c_str());
		string nTuplename = "Craneen_"+ dataSetName;
		nTuple[filename.c_str()] = (TNtuple*)FileObj[filename.c_str()]->Get(nTuplename.c_str());
		nEntries = (int)nTuple[filename.c_str()]->GetEntries();		cout<<"                 nEntries: "<<nEntries<<endl;

		nTuple[filename.c_str()]->SetBranchAddress("HT",&HT);
		nTuple[filename.c_str()]->SetBranchAddress("nJets",&nJets);
		nTuple[filename.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
		nTuple[filename.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
		nTuple[filename.c_str()]->SetBranchAddress("Luminosity",&Luminosity);			//cout<<"got ntuple vars"<<endl;

		//for fixed bin width
		histo1D[filename.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, 0, 1000);

		/////*****loop through entries and fill plots*****
		for (int j = 0; j<nEntries; j++){
			nTuple[filename.c_str()]->GetEntry(j);

			MSPlot[plainplotname.c_str()]->Fill(HT, datasets[d], true, ScaleFactor*Luminosity);
			histo1D[filename.c_str()]->Fill(HT,NormFactor*ScaleFactor*Luminosity);
		}

		if(dataSetName == "TTJets")
		{
			errorfile->cd();
			errorfile->cd("MultiSamplePlot_HT");

			histo1D[filename.c_str()]->Write("Nominal");
									cout<<"write ttjets error"<<endl;
			errorfile->Write();
												cout<<"write ttjets error"<<endl;
			errorfile->Close();
			//delete errorfile;

		}

		shapefile->cd();
		TCanvas *canv = new TCanvas();

		histo1D[filename.c_str()]->Draw();
		histo1D[filename.c_str()]->Write((dataSetName).c_str());

		canv->SaveAs((pathPNG+dataSetName+".pdf").c_str());
	}

		treeLoader.UnLoadDataset();
	
		MSPlot[plainplotname.c_str()]->setErrorBandFile("ScaleFilesMu_light/Error.root");

		for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){
			string name = it->first;
			MultiSamplePlot *temp = it->second;
			//temp->Draw("HT", 0, false, true, false, 100);
			//temp->Draw_wSysUnc(false,"ScaleFilesMu_light", name, true, true, false, false, false,100,true, false, false, true); // merge TT/QCD/W/Z/ST/
			//temp->Draw(false, ("CMSPlot"), true, false, true, true, false, 100, false, false, false);
			//temp->Write(shapefile, name, true, pathPNG, "pdf");
		}

	//delete fout;
};


void SystematicsAnalyser(int nBins, string lepton, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile)
{
	string xmlFileName = ""; xmlFileName = "config/Run2SingleLepton_samples_Sys.xml"; //xmlFileName += lepton; xmlFileName += "Sys.xml";
	const char *xmlfile = xmlFileName.c_str();	cout << "used config file: " << xmlfile << endl; 

	shapefile->cd();

	string pathPNG = "FourTop_SysPlots_" + leptoAbbr; 	//add MSplot name to directory
	mkdir(pathPNG.c_str(),0777);  	cout <<"Making directory :"<< pathPNG  <<endl;		//make directory


	///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
    
    cout<<"loading...."<<endl;
	TTreeLoader treeLoader;
	vector < Dataset* > datasets; cout<<"vector filled"<<endl;
	treeLoader.LoadDatasets (datasets, xmlfile);	cout<<"datasets loaded"<<endl;
	///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
	string dataSetName, filepath;
	int nEntries; 
	float varofInterest, ScaleFactor,NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH;
	//int nBinnings[] = {1, 2, 3, 5, 6, 10, 15, 25, 30, 50, 75, 150};
	cout<<"datsetsize sys: "<<datasets.size()<<endl;
	for (int d = 0; d < datasets.size(); d++){ //Loop through datasets
		dataSetName = datasets[d]->Name();		cout<<"Dataset:  :"<<dataSetName<<endl;
		//filename = "sys_"+dataSetName+"_";      filename += leptoAbbr;
		if (lepton =="Muon")
		{
			filepath = "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/FourTop/Craneens_Mu/Craneens6_2_2015/Craneen_"+dataSetName + "_Run2_TopTree_Study_"+dataSetName+".root";
		}
		else if(lepton == "Electron")
		{
			//filepath = "/Users/lb8075/JamesCode/TopBrussels/LightTrees_El/140206/Sys/"+ filename + ".root";
		}
		else{cout<<"lepton error"<<endl;}
		//cout<<"file in use: "<<filepath<<endl;

		FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());		cout<<"initialised file"<<endl;
		string nTupleName = "Craneen_"+dataSetName;

		//cout<<"nTuple name = "<<nTupleName<<endl;
		nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
		nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();		cout<<"                 nEntries: "<<nEntries<<endl;

		nTuple[dataSetName.c_str()]->SetBranchAddress("HT",&varofInterest);
		//nTuple[dataSetName.c_str()]->SetBranchAddress("nJets",&Njets);
		nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
		nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
		nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);

		////****************************************************Define plots************************************************////
		string plotnamePlain = "MVA_"+dataSetName;
		histo1D[plotnamePlain.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, 0, 1000);

		for (int i = 0; i<nEntries; i++){
			nTuple[dataSetName.c_str()]->GetEntry(i);
			histo1D[plotnamePlain.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
		}

		////****************************************************Fill plots************************************************////
		shapefile->cd();

		TCanvas *canv2 = new TCanvas();
		if(Normalise)
		{
			double dIntegral = histo1D[plotnamePlain.c_str()]->Integral();
			histo1D[plotnamePlain.c_str()]->Scale(1./dIntegral);
		}
		histo1D[plotnamePlain.c_str()]->Draw();
		histo1D[plotnamePlain.c_str()]->Write((plotnamePlain).c_str());
		canv2->SaveAs(("Sys_"+plotnamePlain+".pdf").c_str());

		if(dataSetName == "TTScaledown")
		{
			errorfile->cd();
			errorfile->mkdir("MultiSamplePlot_HT");
			errorfile->cd("MultiSamplePlot_HT");
			histo1D[plotnamePlain.c_str()]->Write("Minus");
			errorfile->Write();
			errorfile->Close();
			//delete errorfile;
		}

		if(dataSetName == "TTScaleup")
		{
			cout<<"    SCALE UP    "<<endl;
			string filenameError = "ScaleFilesMu_light/Error_MVA.root";
			TFile* errorfile = new TFile(filenameError.c_str(),"UPDATE");
			errorfile->cd();
			errorfile->cd("MultiSamplePlot_HT");
			histo1D[plotnamePlain.c_str()]->Write("Plus");
			errorfile->Write();
			errorfile->Close();
			//delete errorfile;
		}
	}
};

