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
#include "TopTreeAnalysisBase/tinyxml/tinyxml.h"
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
map<std::string,std::string> MessageMap;

std::string intToStr (int number);
int howManyElements(TiXmlElement* child, string variableName);

// void SystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlSys, string CraneenPath);
// void DatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath);

// void SplitDatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlNom, string CraneenPath);
// void SplitSystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlSys, string CraneenPath);

void SplitPlusDatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, string sSplitVar2, float ltSplit, float gtSplit, string xmlNom, string CraneenPath, vector<float> splits);
void SplitPlusSystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, string sSplitVar2, float ltSplit, float gtSplit, string xmlSys, string CraneenPath, vector<float> splits);

// void Split2DatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom, string CraneenPath);
// void Split2SystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlSys, string CraneenPath);

// void DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, string xmlNom, float lScale);
// void Split_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string xmlNom, float lScale);
void SplitPlus_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, string sSplitVar1, string sSplitVar2, float ltSplit, float gtSplit, string xmlNom, float lScale, vector<float> splits);
// void Split2_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, bool jetTagsplit, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom, float lScale);


int main(int argc, char** argv)
{
    int NumberOfBins;	//fixed width nBins
    int lumiScale = 50;  //Amount of luminosity to scale to in fb^-1

    bool jetSplit = false; 
    bool jetTagsplit = false;
    bool jetPlusTagSplit = false;

    string DatacardVar = "BDT"; //variable of interest for plotting
    //upper and lower bound of variable in plot
    float lBound, uBound, bSplit, tSplit, wSplit, bSplit1, tSplit1, wSplit1, bSplit2, tSplit2, wSplit2, ltSplit, gtSplit;  // + the bottom, top, and width of the splitting for 1 & 2 variables

    string leptoAbbr, channel, chan, xmlFileName, xmlFileNameSys, CraneenPath, splitVar, splitVar1, splitVar2, VoI;
    string splitting = "inc";

    if(argc>1){
        int iarg=1;
        if(argc>1){
            iarg=1;
            string val = argv[iarg];
            if(val.find("--chan")!=std::string::npos){
                iarg++;
                chan=argv[iarg];
                cout << "now selecting channel : " << chan << endl;
            }
        }
        while(iarg<argc){
            string val = argv[iarg];
            iarg++;
            if(val.find("--JTS")!=std::string::npos){
                jetTagsplit=true;
                splitting = "JTS";
                cout << "!!!!!!! doing Jet tag split !!!" << endl;
            }
            else if(val.find("--JS")!=std::string::npos){
                jetSplit=true;
                splitting = "JS";
                cout << "!!!!!!! doing Jet split !!!" << endl;
            }
            else if(val.find("--JPTS2")!=std::string::npos){
                jetPlusTagSplit=true;
                splitting = "JPTS2";
                cout << "!!!!!!! doing Jet + tag split !!!" << splitting<<endl;
            }
        }
    }
    else throw std::invalid_argument("*** INPUT FLAGS REQUIRED!!! ***");


    string xmlstring = "config/Vars.xml";
    const char * xmlchar = xmlstring.c_str();
    TiXmlDocument doc(xmlchar);
    bool loadOkay = doc.LoadFile();
    if (loadOkay)
    {
        printf("\n%s:\n", xmlchar);
    }
    else
    {
        printf("Failed to load file \"%s\"\n", xmlchar);
        return 0;
    }


    std::string slumiScale = intToStr(lumiScale);

    TiXmlHandle hDoc(&doc);
    TiXmlNode *node = 0;
    node = hDoc.Node ();
    TiXmlHandle hRoot(0);
    cout<<"defined"<<endl;


    {
        TiXmlElement* child = hDoc.FirstChild("analyses").FirstChild("channel").Element();  //get name of channel and check it matches input before continuing
        for (child; child; child=child->NextSiblingElement()){
            string pName = child->Attribute("name");
            if (pName == chan){
                cout<<pName<<endl;
                break;
            }
        }
        //set paths
        leptoAbbr = child->FirstChild( "leptoAbbr")->ToElement()->GetText();
        channel = child->FirstChild( "chan")->ToElement()->GetText();
        xmlFileName = child->FirstChild( "fileName")->ToElement()->GetText();
        xmlFileNameSys = child->FirstChild( "fileNameSys")->ToElement()->GetText();
        CraneenPath = child->FirstChild( "craneenPath")->ToElement()->GetText();
        child->FirstChild("lumiScale")->ToElement()->QueryIntAttribute("lumi", &lumiScale);

        cout<<"leptoAbbr: "<<leptoAbbr<<"  channel: "<<channel<<"  lumiScale: : "<<lumiScale<<"xmlFileName: "<<xmlFileName<<"  xmlFileNameSys: "<<xmlFileNameSys<<"  CraneenPath: "<<CraneenPath<<endl;

        //Get splittings from xml depending on JS or JTS
        TiXmlElement* child2 = child->FirstChild( "splitting" )->ToElement();
        int numberino = howManyElements(child2, "splitting");
        cout<<numberino<<"!!!!!"<<endl;
        if (jetSplit || jetTagsplit || jetPlusTagSplit){
            for (child2; child2; child2=child2->NextSiblingElement("splitting")){
                const char *p3Key=child2->Value();
                const char *p3Text=child2->GetText();
                if (p3Key && p3Text && p3Text==splitting) 
                {
                    cout<<"   splitting: "<<splitting<<endl;
                    break;
                }
            }
        }
        if(jetPlusTagSplit){
            splitVar1 = child2->Attribute("splitVar1");
            splitVar2 = child2->Attribute("splitVar2");
            cout<<"splitVar1: "<<splitVar1<<"   splitVar2: "<<splitVar2<<endl;
            // child2->QueryFloatAttribute("ltSplit", &ltSplit);
            // child2->QueryFloatAttribute("gtSplit", &gtSplit);
        }

        vector<float> splits;
        TiXmlElement* child33 = child2->FirstChild("split")->ToElement();
        for (child33; child33; child33=child33->NextSiblingElement("split")){
            float numberHolder;
            child33->QueryFloatAttribute("n", &numberHolder);
            splits.push_back(numberHolder);
            cout<<"NumberHolder"<< numberHolder<<endl;
        }

        TiXmlElement* child3 = child->FirstChild( "var" )->ToElement();
        for (child3; child3; child3=child3->NextSiblingElement("var")){
            const char *ppKey=child3->Value();
            const char *ppText=child3->GetText();
            if (ppKey && ppText) 
            {
                VoI = ppText;
                string shapefileName = "";
                shapefileName = "shapefile"+leptoAbbr+"_"+slumiScale+"_"+VoI+"_"+splitting+".root";
                cout<<shapefileName<<endl;
                TFile *shapefile = new TFile((shapefileName).c_str(), "RECREATE");
                TFile *errorfile = new TFile(("ScaleFiles"+leptoAbbr+"_light/Error"+VoI+".root").c_str(),"RECREATE");

                child3->QueryFloatAttribute("lBound", &lBound);
                child3->QueryFloatAttribute("uBound", &uBound);
                child3->QueryIntAttribute("nBins", &NumberOfBins);   cout<<"Variable : "<<ppText<<"  lBound : "<<lBound<<"   uBound : "<<uBound<<"  nBins: "<<NumberOfBins<<endl;
                if(jetPlusTagSplit){
                    SplitPlusSystematicsAnalyser(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, false, shapefile, errorfile, channel, VoI, splitVar1, splitVar2, ltSplit, gtSplit, xmlFileNameSys, CraneenPath, splits);
                    SplitPlusDatasetPlotter(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, shapefile, errorfile, channel, VoI, splitVar1, splitVar2, ltSplit, gtSplit, xmlFileName, CraneenPath, splits);
                    if(ppText == DatacardVar){
                        SplitPlus_DataCardProducer(shapefile, shapefileName ,channel, leptoAbbr, jetSplit, splitVar1, splitVar2, ltSplit, gtSplit, xmlFileName, lumiScale, splits);
                    }                    
                }

                errorfile->Close();
                shapefile->Close();
                delete shapefile;
                delete errorfile;
                cout<<""<<endl;                cout<<"end var"<<endl;
            }
        }
    } 
    cout<<" DONE !! "<<endl;
}

void SplitPlusDatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, string sSplitVar2, float ltSplit, float gtSplit, string xmlNom, string CraneenPath, vector<float> splits)
{
    cout<<""<<endl;    cout<<"RUNNING NOMINAL DATASETS"<<endl;    cout<<""<<endl;
    shapefile->cd();

    const char *xmlfile = xmlNom.c_str();    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_Light";
    pathPNG += leptoAbbr;
    pathPNG += "_MSPlots/";
    mkdir(pathPNG.c_str(),0777);    cout <<"Making directory :"<< pathPNG  <<endl;      //make directory

    ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;

    //***************************************************CREATING PLOTS****************************************************
    string plotname, plotname1, plotname0, numStr, dataSetName, filepath, histoName, histoName1;   ///// Jet Split plot

    // plotname0 = "channel0";
    // MSPlot[plotname0.c_str()] = new MultiSamplePlot(datasets, plotname0.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());

    for ( int s = 0; s<splits.size(); s++){
        plotname1 = ("channel"+intToStr(ceil(splits[s])));
        MSPlot[plotname1] = new MultiSamplePlot(datasets, plotname1.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
    }

    plotname = "";

    //***********************************************OPEN FILES & GET NTUPLES**********************************************
    int nEntries;
    float ScaleFactor, NormFactor, Luminosity, varofInterest, splitVar1, splitVar2;

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();    cout<<"Dataset:  :"<<dataSetName<<endl;

        filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study.root";    //cout<<"filepath: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset
        string nTuplename = "Craneen__"+ leptoAbbr;
        nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
        nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();    cout<<"                 nEntries: "<<nEntries<<endl;

        nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
        nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
        nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar1.c_str(), &splitVar1);
        nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar2.c_str(), &splitVar2);

        // string histoName0 = leptoAbbr+"channel0";
        // histo1D[histoName0.c_str()] = new TH1F(histoName1.c_str(),histoName0.c_str(), nBins, plotLow, plotHigh);

        for ( int s = 0; s<splits.size(); s++){
            histoName1 = (channel+intToStr(ceil(splits[s])));
            histo1D[histoName1] = new TH1F(histoName1.c_str(),histoName1.c_str(), nBins, plotLow, plotHigh);
        }

        /////*****loop through entries and fill plots*****
        for (int j = 0; j<nEntries; j++)
        {
            nTuple[dataSetName.c_str()]->GetEntry(j);
            //artificial Lumi
            if(lScale > 0 )  Luminosity = 1000*lScale;

            int splitTotal = splitVar1+splitVar2;
            //cout<<"splitTotal:  "<<splitTotal<<endl;
            // if(splitTotal<splits[0]){
            //     plotname=plotname0;
            //     histoName = histoName0;
            // }
            if (splitTotal>splits[ splits.size()-1 ]){
                plotname = "channel"+intToStr(ceil(splits[ splits.size()-1 ]));
                histoName = channel+intToStr(ceil(splits[ splits.size()-1 ]));
            }
            else{
                for(int j=0; j < (splits.size()-1); j++){
                    //cout<<"j = "<<j<<endl;
                    if(splits[j]<splitTotal && splitTotal<splits[j+1]){
                        plotname = "channel"+intToStr(ceil(splits[ j ]));
                        histoName = channel+intToStr(ceil(splits[ j ]));
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
        /*
        if(dataSetName == "TTJets")  //to put nominal histo into error file
        {
            errorfile->cd();
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest+"/chan1").c_str());
            histo1D[dataSetName + "channel1"]->Write("Nominal");
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest+"/chan2").c_str());
            histo1D[dataSetName + "channel2"]->Write("Nominal");
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest+"/chan3").c_str());
            histo1D[dataSetName + "channel3"]->Write("Nominal");
            errorfile->Write();
        }
        */
        shapefile->cd();
        TCanvas *canv = new TCanvas();

        for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
        {
            string name = it->first;   
            TH1F *temp = it->second;
            histo1D[name]->Draw();
            string writename = "";
            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                writename = name  + "__data_obs__nominal";
            }
            else
            {
                writename = name + "__" + dataSetName +"__nominal";
            }
            histo1D[name]->Write((writename).c_str());
            histo1D.erase(name);
        }    
    }


    treeLoader.UnLoadDataset();
/*
    string scaleFileDir = "ScaleFiles" + leptoAbbr + "_light";
    mkdir(scaleFileDir.c_str(),0777);
    string scaleFileName = scaleFileDir + "/Error"+sVarofinterest+".root";
    //MSPlot[plotname.c_str()]->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot
*/
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;  cout<<name<< " ** "<<sVarofinterest<<endl;
        MultiSamplePlot *temp = it->second;
        //temp->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot
        temp->Draw(sVarofinterest.c_str(), 0, false, false, false, 1);
        temp->Write(shapefile, name, true, pathPNG, "pdf");
        MSPlot.erase(name);

    }

};


void SplitPlusSystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, string sSplitVar2, float ltSplit, float gtSplit, string xmlSys, string CraneenPath, vector<float> splits)
{
    cout<<""<<endl;    cout<<"RUNNING SYS"<<endl;    cout<<""<<endl;
    const char *xmlfile = xmlSys.c_str();    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_SysPlots_" + leptoAbbr;   //add MSplot name to directory
    mkdir(pathPNG.c_str(),0777);    cout <<"Making directory :"<< pathPNG  <<endl;      //make directory

    ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; //cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);    cout<<"datasets loaded"<<endl;

    ///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
    string dataSetName, filepath, numStr, histoName1, histoName2, histoName3, histoName;
    int nEntries;
    float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH, splitVar1, splitVar2;
    string sltSplit = static_cast<ostringstream*>( &(ostringstream() << ltSplit) )->str();
    string sgtSplit = static_cast<ostringstream*>( &(ostringstream() << gtSplit) )->str();
    errorfile->cd();
/*
    errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest+"/chan1").c_str());
    errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest+"/chan2").c_str());
    errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest+"/chan3").c_str());
*/
    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();
        cout<<" Dataset:  :"<<dataSetName<<endl;
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
        nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar1.c_str(), &splitVar1);
        nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar2.c_str(), &splitVar2);



        ////****************************************************Define plots***********************************************
        // string histoName0 = leptoAbbr+"channel0";
        // histo1D[histoName0.c_str()] = new TH1F(histoName1.c_str(),histoName0.c_str(), nBins, plotLow, plotHigh);

        for ( int s = 0; s<splits.size(); s++){
            histoName1 = (channel+intToStr(ceil(splits[s])));
            histo1D[histoName1] = new TH1F(histoName1.c_str(),histoName1.c_str(), nBins, plotLow, plotHigh);
        }

        
        for (int i = 0; i<nEntries; i++)   //Fill histo with variable of interest
        {
            nTuple[dataSetName.c_str()]->GetEntry(i);
            //artificial Lumi

            if(lScale > 0 )
            {
                Luminosity = 1000*lScale;
            }

            int splitTotal = splitVar1+splitVar2;           
            // if(splitTotal<splits[0]){
            //     histoName = histoName0;
            // }
            if (splitTotal>splits[ splits.size()-1 ]){
                histoName = channel+intToStr(ceil(splits[ splits.size()-1 ]));
            }
            else{
                for(int j=0; j < (splits.size()-1); j++){
                    if(splits[j]<splitTotal && splitTotal<splits[j+1]){
                        histoName = channel+intToStr(ceil(splits[ j ]));
                    }
                }
            }

            histo1D[histoName.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
        }

        ////****************************************************Fill plots***********************************************
            

        TCanvas *canv2 = new TCanvas();

        for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
        {
            shapefile->cd();
            string name = it->first;    cout<<name<<endl;
            TH1F *temp = it->second;
            if(Normalise)
            {
                double dIntegral = temp->Integral();
                temp->Scale(1./dIntegral);
            }
            temp->Draw();

            if (dataSetName == "TTScaleDown"){
            histo1D[name]->Write((name+"__TTJets__scaleDown").c_str());    
            /*errorfile->cd();
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest+"/chan1").c_str());
            histo1D[histoName1]->Write("Minus");
            */
            }
            else if (dataSetName == "TTScaleUp"){
            histo1D[name]->Write((name+"__TTJets__scaleUp").c_str());    
            /*errorfile->cd();
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest+"/chan1").c_str());
            histo1D[histoName1]->Write("Minus");
            */
            }
            histo1D.erase(name);
        }

    }
};


void SplitPlus_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, string sSplitVar1, string sSplitVar2, float ltSplit, float gtSplit, string xmlNom, float lScale, vector<float> splits){
    TTreeLoader treeLoader;
    vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
    const char *xmlfile = xmlNom.c_str();
    treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;
    int nDatasets = datasets.size();
    TH1F *tempHisto;
    float tempEntries;
    int nChannels = 0;

    cout<<""<<endl;
    cout<<"PRODUCING DATACARD"<<endl;

    string numStr1, binname, histoName, dataSetName;
    ofstream card;
    std::string slScale = intToStr(lScale);
    string datacardname = "datacard" + leptoAbbr + "_" + slScale + "_JPTS2" + ".txt";
    card.open (datacardname.c_str());
    cout<<"datacard name "<<datacardname<<endl;
    for(int s = 0; s < splits.size(); s++)
    {
        nChannels += 1;
    } 
    card << "imax " + static_cast<ostringstream*>( &(ostringstream() << nChannels) )->str() + "\n"; 
    card << "jmax " + static_cast<ostringstream*>( &(ostringstream() << nDatasets-2) )->str() + "\n";  //number of background (so minus data and signal)
    card << "kmax *\n";
    card << "---------------\n";
    card << "shapes * * "+shapefileName+"  $CHANNEL__$PROCESS__nominal  $CHANNEL__$PROCESS__$SYSTEMATIC\n";
    card << "---------------\n";
    card << "bin                               ";
    for(int s = 0; s < splits.size(); s++)
    {
        numStr1 = static_cast<ostringstream*>( &( ostringstream() << ceil(splits[s]) ) )->str();
        binname = channel + numStr1;
        card << binname + " ";
    }    
    card << "\n";
    card << "observation                               ";

    tempEntries = 0;
    /*    for (int j = 0; j<nDatasets; j++){
        dataSetName=datasets[j]->Name();
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
        {
            histoName = "Muchannel0" + "__data_obs__nominal";
            tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
            tempEntries = tempHisto->GetSumOfWeights();  
            card<<static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str()+"         ";                
        }
        else{
            continue;
        }
    }*/

    for(int s = 0; s < splits.size(); s++)
    {
        int floor_s = ceil(splits[s]);
        numStr1 = static_cast<ostringstream*>( &(ostringstream() << floor_s) )->str();
        //binname = channel + numStr1 + sSplitVar1;
        tempEntries = 0;
        for (int j = 0; j<nDatasets; j++){
            dataSetName=datasets[j]->Name();
            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                histoName = channel + numStr1 + "__data_obs__nominal";  cout<<histoName<<endl;
                tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
                tempEntries = tempHisto->GetSumOfWeights();  
                card<<static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str()+"         ";                
            }
            else{
                continue;
            }
        }
    } 
    card << "\n";

    card << "---------------------------\n";
    card << "bin                               ";
    // for(int i = 0; i<nDatasets-1; i++)    card << "Muchannel0" + " ";

    for(int s = 0; s < splits.size(); s++)
    {
        numStr1 = static_cast<ostringstream*>( &(ostringstream() << ceil(splits[s])) )->str();
        binname = channel + numStr1;
        for(int i = 0; i<nDatasets-1; i++)    card << binname + " ";
    }
    card << "\n";
    card << "process                             ";
    for(int s = 0; s < splits.size(); s++)
    {
        binname = channel + numStr1;
        for (int j = 0; j<nDatasets; j++){
            dataSetName=datasets[j]->Name();
            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                continue;
            }
            else{
                card<<dataSetName+"           ";
            }
        }
    }    

    card << "\n";
    card << "process                                ";
    for(int s = 0; s < splits.size(); s++)
    {
        for(int i=0; i<datasets.size()-1; i++){
            card << static_cast<ostringstream*>( &(ostringstream() << i) )->str() + "                    ";
        }
    }  
    card << "\n";
    card << "rate                                ";
    for(int s = 0; s < splits.size(); s++)
    {
        numStr1 = static_cast<ostringstream*>( &(ostringstream() << ceil(splits[s])) )->str();

        // binname = channel + numStr1;
        tempEntries = 0;
        for (int j = 0; j<nDatasets; j++){
            dataSetName=datasets[j]->Name();
            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
            {
                continue;
            }
            else{
                histoName = channel + numStr1 + "__" + dataSetName +"__nominal"; cout<<histoName<<endl;
                tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
                tempEntries = tempHisto->GetSumOfWeights();  
                card << static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str() + "               ";
            }
        }
    }     
    card << "\n";
    card << "---------------------------\n";
    card << "lumi                  lnN           ";
    
    for (int i=0; i<nChannels*(nDatasets-1); i++){
        card << "1.026                ";
    }
    card << "\n";
    for (int d = 0; d<nDatasets; d++){
        dataSetName = datasets[d]->Name();
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
        {
            continue;
        }
        else {
            if(dataSetName.find("tttt")!=string::npos){
                card << "tttt_norm              lnN          ";
                for (int j = 0; j<nChannels; j++){
                    card << "1.1                 ";
                    for(int k = 0; k<nDatasets-2; k++){
                        card<<"-                    ";
                    }
                }
                card << "\n";
            }
            else{
                card << dataSetName + "_norm      lnN           ";
                for(int k = 0; k<nChannels; k++){
                    for (int dash1 = 0; dash1<d-1; dash1++){
                        card << "-                  ";
                    }
                    card << "1.04                  ";                   
                    for (int dash2 = 5; dash2>d; dash2-- ){
                        card << "-                  ";
                    }
                }
                card<<"\n";
            }   
        }        
    }

    card << "scale                shape           ";
    for (int d = 0; d<nDatasets; d++){
        dataSetName = datasets[d]->Name();
        if(dataSetName.find("TTJets")!=string::npos)
        {
            for(int k = 0; k<nChannels; k++){
                for (int dash1 = 0; dash1<d-1; dash1++){
                    card << "-                  ";
                }
                card << "1                      ";                   
                for (int dash2 = 5; dash2>d; dash2-- ){
                    card << "-                  ";
                }
            }
            card<<"\n";
        }
        else {
            continue;
        }        
    }    

    card.close();
    
};


std::string intToStr (int number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

int howManyElements(TiXmlElement* child, string variableName){
    int counter = 0;
    for (child; child; child=child->NextSiblingElement((variableName).c_str())) {
        counter++;
    }
    return counter;
}
