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

void SystematicsAnalyser(int nBins, string leptoAbbr, bool Normalise, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlSys, string CraneenPath);
void DatasetPlotter(int nBins, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath);


int main()
{
	int NumberOfBins = 10;	//fixed width nBins

	//------- Set Channel --------//
	bool DileptonMuEl = false;
	bool SingleMu = true;
	bool SingleEl = false;

	string VoI = "HT"; //variable of interest for plotting

	string leptoAbbr;
	string channel;
	string xmlFileName;
	string xmlFileNameSys;
	string CraneenPath;

	if(SingleMu)
	{
		leptoAbbr = "Mu";
		channel = "ttttmu__";
		xmlFileName = "config/Run2SingleLepton_samples.xml";	
		xmlFileNameSys = "config/Run2SingleLepton_samples_Sys.xml";	
		CraneenPath = "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/FourTop/Craneens_Mu/Craneens5_2_2015/Craneen_";
	}
	else if(SingleEl)
	{
		leptoAbbr = "El";
		channel = "ttttel__";
		xmlFileName = "config/Run2SingleLepton_samples.xml";
		xmlFileNameSys = "config/Run2SingleLepton_samples_Sys.xml";
		CraneenPath = "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/FourTop/Craneens_Mu/Craneens5_2_2015/Craneen_";
	}
	else if(DileptonMuEl)
	{
		leptoAbbr = "MuEl";
		channel = "ttttmuel__";
		xmlFileName = "config/Run2_Samples.xml";
		xmlFileNameSys = "config/Run2SingleLepton_samples_Sys.xml";
		CraneenPath = "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/FourTop/Craneens_Mu/Craneens5_2_2015/Craneen_";

	}

	TFile *shapefile = new TFile(("shapefile"+leptoAbbr+".root").c_str(), "RECREATE");
	TFile *errorfile = new TFile("ScaleFilesMu_light/Error.root","RECREATE");

	SystematicsAnalyser(NumberOfBins, leptoAbbr, false, shapefile, errorfile, channel, VoI, xmlFileNameSys, CraneenPath);
	DatasetPlotter(NumberOfBins, leptoAbbr, shapefile, errorfile, channel, VoI, xmlFileName, CraneenPath);

	errorfile->Close();
	shapefile->Close();
	delete shapefile;
	delete errorfile;
}


void DatasetPlotter(int nBins, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath)
{	
	cout<<""<<endl; cout<<"RUNNING NOMINAL DATASETS"<<endl; cout<<""<<endl;
	shapefile->cd();

	const char *xmlfile = xmlNom.c_str();	cout << "used config file: " << xmlfile << endl; 

	string pathPNG = "FourTop_Light"; 	pathPNG += leptoAbbr;	pathPNG += "_MSPlots/";
	mkdir(pathPNG.c_str(),0777);  	cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

	///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
	TTreeLoader treeLoader;
	vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
	treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

	//***************************************************CREATING PLOTS****************************************************
	string plotname = sVarofinterest;   ///// Non Jet Split plot
	MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, 0, 1000, "BDT Discriminator");

	//***********************************************OPEN FILES & GET NTUPLES**********************************************
	string dataSetName, filepath;		int nEntries; 
	float ScaleFactor, NormFactor, Luminosity, varofInterest;

	for (int d = 0; d < datasets.size(); d++){ //Loop through datasets
		dataSetName = datasets[d]->Name();		cout<<"Dataset:  :"<<dataSetName<<endl;
		
		filepath = CraneenPath+dataSetName + "_Run2_TopTree_Study_"+dataSetName+".root";
		//cout<<"filepath: "<<filepath<<endl; 

		FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset
		string nTuplename = "Craneen_"+ dataSetName;
		nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
		nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();			cout<<"                 nEntries: "<<nEntries<<endl;

		nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
		nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
		nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
		nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);			

		//for fixed bin width
		histo1D[dataSetName.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, 0, 1000);
		/////*****loop through entries and fill plots*****
		for (int j = 0; j<nEntries; j++){
			nTuple[dataSetName.c_str()]->GetEntry(j);

			MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, ScaleFactor*Luminosity);
			histo1D[dataSetName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
		}
		if(dataSetName == "TTJets")  //to put nominal histo into error file
		{
			errorfile->cd();
			errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
			histo1D[dataSetName.c_str()]->Write("Nominal");
			errorfile->Write();
		}

		shapefile->cd();
		TCanvas *canv = new TCanvas();

		histo1D[dataSetName.c_str()]->Draw();
		string writename = ""; writename = channel + dataSetName +"__nominal";	cout<<"writename  :"<<writename<<endl;
		histo1D[dataSetName.c_str()]->Write((writename).c_str());

		canv->SaveAs((pathPNG+dataSetName+".pdf").c_str());
	}


	treeLoader.UnLoadDataset();

	MSPlot[plotname.c_str()]->setErrorBandFile("ScaleFilesMu_light/Error.root"); //set error file for uncertainty bands on multisample plot

	for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++){
		string name = it->first;
		MultiSamplePlot *temp = it->second;
		temp->Draw(sVarofinterest.c_str(), 0, false, true, false, 100);
		temp->Write(shapefile, name, true, pathPNG, "pdf");
	}
};


void SystematicsAnalyser(int nBins, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlSys, string CraneenPath)
{
	cout<<""<<endl; cout<<"RUNNING SYS"<<endl; cout<<""<<endl;
	const char *xmlfile = xmlSys.c_str();	cout << "used config file: " << xmlfile << endl; 

	string pathPNG = "FourTop_SysPlots_" + leptoAbbr; 	//add MSplot name to directory
	mkdir(pathPNG.c_str(),0777);  	cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

	///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
	TTreeLoader treeLoader;
	vector < Dataset* > datasets; //cout<<"vector filled"<<endl; 
	treeLoader.LoadDatasets (datasets, xmlfile);	cout<<"datasets loaded"<<endl;

	///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
	string dataSetName, filepath;	int nEntries;	float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH;

	for (int d = 0; d < datasets.size(); d++){ //Loop through datasets
		dataSetName = datasets[d]->Name();		cout<<"	Dataset:  :"<<dataSetName<<endl;
		filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study_"+dataSetName+".root";

		cout<<""<<endl;	cout<<"file in use: "<<filepath<<endl;

		FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());		cout<<"initialised file"<<endl;
		string nTupleName = "Craneen_"+dataSetName;

		nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
		nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();		cout<<"                 nEntries: "<<nEntries<<endl;

		//sVarofinterest = "HT";
		nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
		nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
		nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
		nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);

		////****************************************************Define plots***********************************************
		string plotname = sVarofinterest+"_"+dataSetName;
		histo1D[plotname.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, 0, 1000);

		for (int i = 0; i<nEntries; i++){  //Fill histo with variable of interest
			nTuple[dataSetName.c_str()]->GetEntry(i);
			histo1D[plotname.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
		}

		////****************************************************Fill plots***********************************************
		shapefile->cd();

		TCanvas *canv2 = new TCanvas();
		if(Normalise)
		{
			double dIntegral = histo1D[plotname.c_str()]->Integral();
			histo1D[plotname.c_str()]->Scale(1./dIntegral);
		}
		histo1D[plotname.c_str()]->Draw();
		string writename = ""; writename = channel + "TTJets__" + dataSetName;	cout<<"writename  :"<<writename<<endl;

		histo1D[plotname.c_str()]->Write((writename).c_str());
		canv2->SaveAs(("Sys_"+plotname+".pdf").c_str());

		if(dataSetName == "TTScaledown")
		{
			errorfile->cd();
			errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
			errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
			histo1D[plotname.c_str()]->Write("Minus");
			//errorfile->Write();
		}

		if(dataSetName == "TTScaleup")
		{
			errorfile->cd();
			errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
			histo1D[plotname.c_str()]->Write("Plus");
			//errorfile->Write();
		}
	}
};

