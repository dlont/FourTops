#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include "TRandom3.h"
#include "TNtuple.h"

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
//#include "../macros/Style.C"

#include <sstream>

using namespace std;
using namespace TopTree;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TFile*> FileObj;
map<string,TNtuple*> nTuple;
map<string,MultiSamplePlot*> MSPlot;

std::string intToStr (int number);

void SystematicsAnalyser(int nBins, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlSys, string CraneenPath);
void DatasetPlotter(int nBins, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath);

void SplitDatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlNom, string CraneenPath);
void SplitSystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlSys, string CraneenPath);



int main()
{
    int NumberOfBins = 30;	//fixed width nBins
    int lumiScale = 50;  //Amount of luminosity to scale to in fb^-1

    //------- Set Channel --------//
    bool DileptonMuEl = false;
    bool DileptonMuMu = false;
    bool DileptonElEl = true;
    bool SingleMu = false;
    bool SingleEl = false;
    bool jetSplit = true;

    string VoI = "BDT"; //variable of interest for plotting
    float uBound = 1.0;
    float lBound = -0.5;
    vector<string> vars;
    vars.push_back("HT");
    vars.push_back("LeadingMuonPt");
    vars.push_back("LeadingElectronPt");
    vars.push_back("LeadingBJetPt");
    vars.push_back("HT2M");
    vars.push_back("MVAvals1");

    string leptoAbbr;
    string channel;
    string xmlFileName;
    string xmlFileNameSys;
    string CraneenPath;
    string splitVar;

    float bSplit, tSplit, wSplit;  //the bottom, top, and width of the splitting

    if(SingleMu)
    {
        leptoAbbr = "Mu";
        channel = "ttttmu";
        xmlFileName = "config/Run2SingleLepton_samples.xml";
        xmlFileNameSys = "config/Run2SingleLepton_samples_Sys.xml";
        CraneenPath = "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/FourTop/Craneens_Mu/Craneens5_2_2015/Craneen_";
    }
    else if(SingleEl)
    {
        leptoAbbr = "El";
        channel = "ttttel";
        xmlFileName = "config/Run2SingleLepton_samples.xml";
        xmlFileNameSys = "config/Run2SingleLepton_samples_Sys.xml";
        CraneenPath = "/user/lbeck/ThirteenTeV/CMSSW_7_2_1_patch1/src/TopBrussels/FourTop/Craneens_Mu/Craneens5_2_2015/Craneen_";
    }
    else if(DileptonMuEl)
    {
        leptoAbbr = "MuEl";
        channel = "ttttmuel";
        xmlFileName = "config/Run2DiLepton_Craneens_Nom.xml";
        xmlFileNameSys = "config/Run2DiLepton_Craneens_Sys.xml";
        CraneenPath = "/user/heilman/CMSSW_7_2_1_patch1/src/TopBrussels/FourTopsLight/Craneens_MuEl/Craneens18_3_2015/merge/Craneen_";

    }
    else if(DileptonMuMu)
    {
        leptoAbbr = "MuMu";
        channel = "ttttmumu";
        xmlFileName = "config/Run2DiLepton_Craneens_Nom.xml";
        xmlFileNameSys = "config/Run2DiLepton_Craneens_Sys.xml";
        CraneenPath = "/user/heilman/CMSSW_7_2_1_patch1/src/TopBrussels/FourTopsLight/Craneens_MuMu/Craneens17_3_2015/merge/Craneen_";

    }
    else if(DileptonElEl)
    {
        leptoAbbr = "ElEl";
        channel = "ttttelel";
        xmlFileName = "config/Run2DiLepton_Craneens_Nom.xml";
        xmlFileNameSys = "config/Run2DiLepton_Craneens_Sys.xml";
        CraneenPath = "/user/heilman/CMSSW_7_2_1_patch1/src/TopBrussels/FourTopsLight/Craneens_ElEl/Craneens18_3_2015/merge/Craneen_";

    }

    std::string slumiScale = intToStr(lumiScale);

    TFile *shapefile = new TFile(("shapefile"+leptoAbbr+"_"+slumiScale+"_"+VoI+".root").c_str(), "RECREATE");
    TFile *errorfile = new TFile(("ScaleFiles"+leptoAbbr+"_light/Error.root").c_str(),"RECREATE");

    if(jetSplit)
    {
        //Control variables for splitting in nJets
        splitVar = "nJets";
        bSplit = 4; //Lower bound of jetSplit bins
        tSplit = 8; //First bin no longer bound by bin width.  This bin contains all information up to infinity in the splitVar
        wSplit = 1; //width of the bins
        SplitSystematicsAnalyser(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, false, shapefile, errorfile, channel, VoI, splitVar, bSplit, tSplit, wSplit, xmlFileNameSys, CraneenPath);
        SplitDatasetPlotter(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, shapefile, errorfile, channel, VoI, splitVar, bSplit, tSplit, wSplit, xmlFileName, CraneenPath);
//        for(int k=0; k<vars.size(); k++)
//        {
//            string varchannel = channel + vars[k];
//            SplitSystematicsAnalyser(NumberOfBins, leptoAbbr, false, shapefile, errorfile, varchannel, vars[k], splitVar, bSplit, tSplit, wSplit, xmlFileNameSys, CraneenPath);
//            SplitDatasetPlotter(NumberOfBins, leptoAbbr, shapefile, errorfile, varchannel, vars[k], splitVar, bSplit, tSplit, wSplit, xmlFileName, CraneenPath);
//        }
    }
    else
    {
        SystematicsAnalyser(NumberOfBins, lBound, uBound, leptoAbbr, false, shapefile, errorfile, channel, VoI, xmlFileNameSys, CraneenPath);
        DatasetPlotter(NumberOfBins, lBound, uBound, leptoAbbr, shapefile, errorfile, channel, VoI, xmlFileName, CraneenPath);
    }


    errorfile->Close();
    shapefile->Close();
    delete shapefile;
    delete errorfile;
}


void DatasetPlotter(int nBins, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING NOMINAL DATASETS"<<endl;
    cout<<""<<endl;
    shapefile->cd();

    const char *xmlfile = xmlNom.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_Light";
    pathPNG += leptoAbbr;
    pathPNG += "_MSPlots/";
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

    //***************************************************CREATING PLOTS****************************************************
    string plotname = sVarofinterest;   ///// Non Jet Split plot
    MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, "HT2M");

    //***********************************************OPEN FILES & GET NTUPLES**********************************************
    string dataSetName, filepath;
    int nEntries;
    float ScaleFactor, NormFactor, Luminosity, varofInterest;

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"Dataset:  :"<<dataSetName<<endl;

        filepath = CraneenPath+dataSetName + "_Run2_TopTree_Study.root";
        //cout<<"filepath: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset
        string nTuplename = "Craneen__"+ leptoAbbr;
        nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
        nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

        nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
        nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);



        //for fixed bin width
        histo1D[dataSetName.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, plotLow, plotHigh);
        /////*****loop through entries and fill plots*****
        for (int j = 0; j<nEntries; j++)
        {
            nTuple[dataSetName.c_str()]->GetEntry(j);
            //artificial Lumi
            Luminosity = 15000;

            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, NormFactor*ScaleFactor*Luminosity);
                histo1D[dataSetName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }
            else
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, ScaleFactor*Luminosity);
                histo1D[dataSetName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }

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
        string writename = "";
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
        {
            writename = channel + "__data_obs__nominal";
        }
        else
        {
            writename = channel + "__" + dataSetName +"__nominal";
        }
        cout<<"writename  :"<<writename<<endl;
        histo1D[dataSetName.c_str()]->Write((writename).c_str());

        canv->SaveAs((pathPNG+dataSetName+".pdf").c_str());
    }


    treeLoader.UnLoadDataset();

    string scaleFileDir = "ScaleFiles" + leptoAbbr + "_light";
    mkdir(scaleFileDir.c_str(),0777);
    string scaleFileName = scaleFileDir + "/Error.root";
    MSPlot[plotname.c_str()]->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot

    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->Draw(sVarofinterest.c_str(), 0, false, true, true, 1);
        temp->Write(shapefile, name, true, pathPNG, "pdf");
    }
};


void SystematicsAnalyser(int nBins, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlSys, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING SYS"<<endl;
    cout<<""<<endl;
    const char *xmlfile = xmlSys.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_SysPlots_" + leptoAbbr; 	//add MSplot name to directory
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; //cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);
    cout<<"datasets loaded"<<endl;

    ///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
    string dataSetName, filepath;
    int nEntries;
    float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH;

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"	Dataset:  :"<<dataSetName<<endl;
        filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study"+".root";

        cout<<""<<endl;
        cout<<"file in use: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());
        cout<<"initialised file"<<endl;
        string nTupleName = "Craneen__"+leptoAbbr;

        nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
        nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

        //sVarofinterest = "HT";
        nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
        nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);



        ////****************************************************Define plots***********************************************
        string plotname = sVarofinterest+"_"+dataSetName;
        histo1D[plotname.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, plotLow, plotHigh);

        for (int i = 0; i<nEntries; i++)   //Fill histo with variable of interest
        {
            nTuple[dataSetName.c_str()]->GetEntry(i);
            //artificial Lumi
            Luminosity = 15000;
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
        string writename = "";
        writename = channel + "__TTJets__" + dataSetName;
        cout<<"writename  :"<<writename<<endl;

        histo1D[plotname.c_str()]->Write((writename).c_str());
        canv2->SaveAs(("Sys_"+plotname+".pdf").c_str());

        if(dataSetName == "TTScaleDown")
        {
            errorfile->cd();
            errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
            histo1D[plotname.c_str()]->Write("Minus");
            //errorfile->Write();
        }

        if(dataSetName == "TTScaleUp")
        {
            errorfile->cd();
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
            histo1D[plotname.c_str()]->Write("Plus");
            //errorfile->Write();
        }
    }
};

void SplitDatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlNom, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING NOMINAL DATASETS"<<endl;
    cout<<""<<endl;
    shapefile->cd();

    const char *xmlfile = xmlNom.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_Light";
    pathPNG += leptoAbbr;
    pathPNG += "_MSPlots/";
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

    //***************************************************CREATING PLOTS****************************************************
    string plotname;   ///// Jet Split plot
    string numStr;
    string sbSplit = static_cast<ostringstream*>( &(ostringstream() << fbSplit) )->str();
    string stSplit = static_cast<ostringstream*>( &(ostringstream() << ftSplit) )->str();
    string swSplit = static_cast<ostringstream*>( &(ostringstream() << fwSplit) )->str();
    for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
    {
        numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
        plotname = sVarofinterest + numStr + sSplitVar;
        MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
    }
    plotname = "";

    //***********************************************OPEN FILES & GET NTUPLES**********************************************
    string dataSetName, filepath, histoName;
    int nEntries;
    float ScaleFactor, NormFactor, Luminosity, varofInterest, splitVar;

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"Dataset:  :"<<dataSetName<<endl;

        filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study.root";
        //cout<<"filepath: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset
        string nTuplename = "Craneen__"+ leptoAbbr;
        nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
        nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

        nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
        nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
        nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar.c_str(), &splitVar);



        //for fixed bin width
        for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
        {
            numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            histoName = dataSetName + numStr + sSplitVar;
            histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);
        }
        /////*****loop through entries and fill plots*****
        for (int j = 0; j<nEntries; j++)
        {
            nTuple[dataSetName.c_str()]->GetEntry(j);
            //artificial Lumi
            if(lScale > 0 )
            {
                Luminosity = 1000*lScale;
            }


            if(splitVar >= ftSplit) //Check if this entry belongs in the last bin.  Done here to optimize number of checks
            {
                plotname = sVarofinterest + stSplit + sSplitVar;
                histoName = dataSetName + stSplit + sSplitVar;
            }
            else //If it doesn't belong in the last bin, find out which it belongs in
            {
                for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
                {
                    if(splitVar>=s && splitVar<(s+fwSplit)) //splitVar falls inside one of the bins
                    {
                        numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
                        plotname = sVarofinterest + numStr + sSplitVar;
                        histoName = dataSetName + numStr + sSplitVar;
                        break;
                    }
                }
            }

            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, NormFactor*ScaleFactor*Luminosity);
                histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }
            else
            {
                MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, ScaleFactor*Luminosity);
                histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
            }

        }
        if(dataSetName == "TTJets")  //to put nominal histo into error file
        {
            errorfile->cd();
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
            histo1D[histoName.c_str()]->Write("Nominal");
            errorfile->Write();
        }

        shapefile->cd();
        TCanvas *canv = new TCanvas();
        for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
        {
            numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            plotname = sVarofinterest + numStr + sSplitVar;
            histoName = dataSetName + numStr + sSplitVar;
            histo1D[histoName.c_str()]->Draw();
            string writename = "";
            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                writename = channel + numStr + sSplitVar + "__data_obs__nominal";
            }
            else
            {
                writename = channel + numStr + sSplitVar + "__" + dataSetName +"__nominal";
            }
            cout<<"writename  :"<<writename<<endl;
            histo1D[histoName.c_str()]->Write((writename).c_str());

            canv->SaveAs((pathPNG+histoName+".pdf").c_str());
        }
    }


    treeLoader.UnLoadDataset();

    string scaleFileDir = "ScaleFiles" + leptoAbbr + "_light";
    mkdir(scaleFileDir.c_str(),0777);
    string scaleFileName = scaleFileDir + "/Error.root";
    //MSPlot[plotname.c_str()]->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot

    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot
        temp->Draw(sVarofinterest.c_str(), 0, false, true, true, 1);
        temp->Write(shapefile, name, true, pathPNG, "pdf");
    }
};


void SplitSystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlSys, string CraneenPath)
{
    cout<<""<<endl;
    cout<<"RUNNING SYS"<<endl;
    cout<<""<<endl;
    const char *xmlfile = xmlSys.c_str();
    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_SysPlots_" + leptoAbbr; 	//add MSplot name to directory
    mkdir(pathPNG.c_str(),0777);
    cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; //cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);
    cout<<"datasets loaded"<<endl;

    ///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
    string dataSetName, filepath, numStr, histoName;
    int nEntries;
    float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH, splitVar;
    string sbSplit = static_cast<ostringstream*>( &(ostringstream() << fbSplit) )->str();
    string stSplit = static_cast<ostringstream*>( &(ostringstream() << ftSplit) )->str();
    string swSplit = static_cast<ostringstream*>( &(ostringstream() << fwSplit) )->str();

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<"	Dataset:  :"<<dataSetName<<endl;
        filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study"+".root";

        cout<<""<<endl;
        cout<<"file in use: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());
        cout<<"initialised file"<<endl;
        string nTupleName = "Craneen__"+leptoAbbr;

        nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
        nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();
        cout<<"                 nEntries: "<<nEntries<<endl;

        //sVarofinterest = "HT";
        nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
        nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
        nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar.c_str(), &splitVar);



        ////****************************************************Define plots***********************************************
        for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
        {
            numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            histoName = dataSetName + numStr + sSplitVar;
            histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);
        }

        for (int i = 0; i<nEntries; i++)   //Fill histo with variable of interest
        {
            nTuple[dataSetName.c_str()]->GetEntry(i);
            //artificial Lumi
            if(lScale > 0 )
            {
                Luminosity = 1000*lScale;
            }

            if(splitVar >= ftSplit) //Check if this entry belongs in the last bin.  Done here to optimize number of checks
            {
                histoName = dataSetName + stSplit + sSplitVar;
            }
            else //If it doesn't belong in the last bin, find out which it belongs in
            {
                for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
                {
                    if(splitVar>=s && splitVar<(s+fwSplit)) //splitVar falls inside one of the bins
                    {
                        numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
                        histoName = dataSetName + numStr + sSplitVar;
                        break;
                    }
                }
            }


            histo1D[histoName.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
        }

        ////****************************************************Fill plots***********************************************


        TCanvas *canv2 = new TCanvas();

        for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
        {
            shapefile->cd();

            numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
            histoName = dataSetName + numStr + sSplitVar;
            if(Normalise)
            {
                double dIntegral = histo1D[histoName.c_str()]->Integral();
                histo1D[histoName.c_str()]->Scale(1./dIntegral);
            }
            histo1D[histoName.c_str()]->Draw();
            string writename = "";

//            writename = channel + numStr + sSplitVar + "__" + dataSetName +"__nominal";
//
//            cout<<"writename  :"<<writename<<endl;
//            histo1D[histoName.c_str()]->Write((writename).c_str());
//
//            canv->SaveAs(("Sys_"+histoName+".pdf").c_str());

            if(dataSetName == "TTScaleDown")
            {
                writename = channel + numStr + sSplitVar + "__TTJets__scaleDown";
                cout<<"writename  :"<<writename<<endl;
                histo1D[histoName.c_str()]->Write((writename).c_str());
                canv2->SaveAs(("Sys_"+histoName+".pdf").c_str());
                errorfile->cd();
                errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
                errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
                histo1D[histoName.c_str()]->Write("Minus");
                //errorfile->Write();
            }

            if(dataSetName == "TTScaleUp")
            {
                writename = channel + numStr + sSplitVar + "__TTJets__scaleUp";
                cout<<"writename  :"<<writename<<endl;
                histo1D[histoName.c_str()]->Write((writename).c_str());
                canv2->SaveAs(("Sys_"+histoName+".pdf").c_str());
                errorfile->cd();
                errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
                histo1D[histoName.c_str()]->Write("Plus");
                //errorfile->Write();
            }
        }


//        histo1D[plotname.c_str()]->Draw();
//        string writename = "";
//        writename = channel + "__TTJets__" + dataSetName;
//        cout<<"writename  :"<<writename<<endl;
//
//        histo1D[plotname.c_str()]->Write((writename).c_str());
//        canv2->SaveAs(("Sys_"+plotname+".pdf").c_str());
//
//        if(dataSetName == "TTScaleDown")
//        {
//            errorfile->cd();
//            errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
//            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//            histo1D[plotname.c_str()]->Write("Minus");
//            //errorfile->Write();
//        }
//
//        if(dataSetName == "TTScaleUp")
//        {
//            errorfile->cd();
//            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//            histo1D[plotname.c_str()]->Write("Plus");
//            //errorfile->Write();
//        }
    }
};

std::string intToStr (int number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

