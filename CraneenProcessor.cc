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
void dump_to_stdout(const char* pFilename);

void SystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlSys, string CraneenPath);
void DatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath, string shapefileName);

void SplitDatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlNom, string CraneenPath);
void SplitSystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlSys, string CraneenPath);

void Split2DatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom, string CraneenPath);
void Split2SystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlSys, string CraneenPath);

void DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, string xmlNom, float lScale);
void Split_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string xmlNom, float lScale);
void Split2_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, bool jetTagsplit, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom, float lScale);

void GetScaleEnvelope(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath, string shapefileName);
void GetScaleEnvelopeSplit(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlNom, string CraneenPath, string shapefileName);

string DatacardVar = "BDT"; //variable of interest for plotting //global


int main(int argc, char** argv)
{
    int NumberOfBins;	//fixed width nBins
    float lumiScale = -1;  //Amount of luminosity to scale to in fb^-1
    cout<<"!!!"<<endl;
    cout<<"!!!"<<endl;
    cout<<"!!!"<<endl;
    cout<<"!!!"<<endl;
    cout<<"!!!"<<endl;

    bool jetSplit = false; 
    bool jetTagsplit = false;
    int isplit_ttbar = 0;
    bool split_ttbar = false;

    //upper and lower bound of variable in plot
    float lBound, uBound, bSplit, tSplit, wSplit, bSplit1, tSplit1, wSplit1, bSplit2, tSplit2, wSplit2;  // + the bottom, top, and width of the splitting for 1 & 2 variables

    string leptoAbbr, channel, chan, xmlFileName, xmlFileNameSys, CraneenPath, splitVar, splitVar1, splitVar2, VoI;
    string splitting = "inc";

    if(argc>0){
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
        else{
            cout<<"no arguments supplied"<<endl;
            return 0;
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
        }

    }



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


    //std::string slumiScale = intToStr(lumiScale);

    TiXmlHandle hDoc(&doc);
    TiXmlElement* pElem;
    TiXmlElement* qElem;
    TiXmlNode *node = 0;
    node = hDoc.Node ();
    TiXmlHandle hRoot(0);
    cout<<"defined"<<endl;


    {
        pElem=hDoc.FirstChildElement().Element();
        hRoot=TiXmlHandle(pElem);
        //TiXmlElement* child = hRoot.FirstChild("analyses").FirstChild().Element();
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

        //split ttbar
        //TiXmlElement* childttbarsplit = child->FirstChild( "split_ttbar" )->ToElement();
        // split_ttbar =  childttbarsplit->Attribute("split"); 
        // childttbarsplit->QueryIntAttribute("split", &isplit_ttbar);
        // split_ttbar = isplit_ttbar;
        // cout<<"TTBAR SPLIT "<<split_ttbar<<endl;
        cout<<"leptoAbbr: "<<leptoAbbr<<"  channel: "<<channel<<"  xmlFileName: "<<xmlFileName<<"  xmlFileNameSys: "<<xmlFileNameSys<<"  CraneenPath: "<<CraneenPath<<endl;

        //Get splittings from xml depending on JS or JTS
        TiXmlElement* child3 = child->FirstChild( "splitting" )->ToElement();
        if (jetSplit || jetTagsplit){
            for (child3; child3; child3=child3->NextSiblingElement()){
                const char *p3Key=child3->Value();
                const char *p3Text=child3->GetText();
                if (p3Key && p3Text && p3Text==splitting) 
                {
                    cout<<"splitting: "<<splitting<<endl;
                    break;
                }
            }
        }
        if(jetSplit){
            splitVar = child3->Attribute("splitVar");
            child3->QueryFloatAttribute("bSplit", &bSplit);
            child3->QueryFloatAttribute("tSplit", &tSplit);
            child3->QueryFloatAttribute("wSplit", &wSplit);     cout<<"splitVar: "<<splitVar<<"  b: "<<bSplit<<"  t: "<<tSplit<<"  w: "<<wSplit<<endl;
        }
        else if(jetTagsplit){
            splitVar1 = child3->Attribute("splitVar1");
            splitVar2 = child3->Attribute("splitVar2");
            child3->QueryFloatAttribute("bSplit1", &bSplit1);
            child3->QueryFloatAttribute("tSplit1", &tSplit1);
            child3->QueryFloatAttribute("wSplit1", &wSplit1);
            child3->QueryFloatAttribute("bSplit2", &bSplit2);
            child3->QueryFloatAttribute("tSplit2", &tSplit2);
            child3->QueryFloatAttribute("wSplit2", &wSplit2);   cout<<"splitVar1: "<<splitVar1<<"splitVar2: "<<splitVar2<<"  b1: "<<bSplit1<<"  t1: "<<tSplit1<<"  w1: "<<wSplit1<<"  b2: "<<bSplit2<<"  t2: "<<tSplit2<<"  w2: "<<wSplit2<<endl;                 
        }

        TiXmlElement* child2 = child->FirstChild( "var" )->ToElement();
        for (child2; child2; child2=child2->NextSiblingElement()){
            const char *ppKey=child2->Value();
            const char *ppText=child2->GetText();
            if (ppKey && ppText) 
            {
                VoI = ppText;
                string shapefileName = "";
                shapefileName = "shapefile"+leptoAbbr+"_"+"DATA"+"_"+VoI+"_"+splitting+".root";
                cout<<shapefileName<<endl;
                TFile *shapefile = new TFile((shapefileName).c_str(), "RECREATE");
                TFile *errorfile = new TFile(("ScaleFiles"+leptoAbbr+"_light/Error"+VoI+".root").c_str(),"RECREATE");

                child2->QueryFloatAttribute("lBound", &lBound);
                child2->QueryFloatAttribute("uBound", &uBound);
                child2->QueryIntAttribute("nBins", &NumberOfBins);   cout<<"Variable : "<<ppText<<"  lBound : "<<lBound<<"   uBound : "<<uBound<<"  nBins: "<<NumberOfBins<<endl;
                if(jetSplit)
                {
                    // SplitSystematicsAnalyser(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, false, shapefile, errorfile, channel, VoI, splitVar, bSplit, tSplit, wSplit, xmlFileNameSys, CraneenPath);
                    // SplitDatasetPlotter(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, shapefile, errorfile, channel, VoI, splitVar, bSplit, tSplit, wSplit, xmlFileName, CraneenPath);
                    // if(ppText == DatacardVar){
                    //     Split_DataCardProducer(shapefile, shapefileName ,channel, leptoAbbr, jetSplit, splitVar, bSplit, tSplit, wSplit, xmlFileName, lumiScale);
                    // }
                }
                else if(jetTagsplit){

                    // Split2SystematicsAnalyser(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, false, shapefile, errorfile, channel, VoI, splitVar1, bSplit1, tSplit1, wSplit1,splitVar2, bSplit2, tSplit2, wSplit2, xmlFileNameSys, CraneenPath);
                    // Split2DatasetPlotter(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, shapefile, errorfile, channel, VoI, splitVar1, bSplit1, tSplit1, wSplit1, splitVar2, bSplit2, tSplit2, wSplit2, xmlFileName, CraneenPath);
                    // if(ppText == DatacardVar){     
                    //     Split2_DataCardProducer(shapefile, shapefileName ,channel, leptoAbbr, jetSplit, jetTagsplit, splitVar1, bSplit1, tSplit1, wSplit1,splitVar2, bSplit2, tSplit2, wSplit2, xmlFileName, lumiScale);
                    // }
                }
                else
                {
                    //SystematicsAnalyser(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, false, shapefile, errorfile, channel, VoI, xmlFileNameSys, CraneenPath);
                    DatasetPlotter(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, shapefile, errorfile, channel, VoI, xmlFileName, CraneenPath, shapefileName);
                    if(ppText == DatacardVar){
                        DataCardProducer(shapefile, shapefileName ,channel, leptoAbbr, xmlFileName, lumiScale);
                    }
                }

                shapefile->Close();
                // GetScaleEnvelope(NumberOfBins, lumiScale, lBound, uBound, leptoAbbr, shapefile, errorfile, channel, VoI, xmlFileName, CraneenPath, shapefileName);

                // errorfile->Close();


                delete shapefile;
                delete errorfile;
                cout<<""<<endl;                cout<<"end var"<<endl;
            }
        }
    } 
    cout<<" DONE !! "<<endl;
}

void GetScaleEnvelope(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath, string shapefileName){

    const char *xmlfile = xmlNom.c_str();    cout << "used config file: " << xmlfile << endl;

    ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
    TTreeLoader treeLoader;
    vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;

    for (int d = 0; d < datasets.size(); d++){  //Loop through datasets
        string dataSetName = datasets[d]->Name();
        if(dataSetName.find("TTJets")!=string::npos){
            cout<<"Producing envelope for scale uncertainty"<<endl;
            histo1D["weightPlus"] = new TH1F("Plus","Plus", nBins, plotLow, plotHigh);
            histo1D["weightMinus"] = new TH1F("Minus","Minus", nBins, plotLow, plotHigh);

            TFile * reopenSF = new TFile(shapefileName.c_str());
            histo1D["weight0"] = (TH1F*)reopenSF->Get("Genweight_tt");
            histo1D["weight1"] = (TH1F*)reopenSF->Get("weight1");
            histo1D["weight2"] = (TH1F*)reopenSF->Get("weight2");
            histo1D["weight3"] = (TH1F*)reopenSF->Get("weight3");
            histo1D["weight4"] = (TH1F*)reopenSF->Get("weight4");
            histo1D["weight5"] = (TH1F*)reopenSF->Get("weight5");
            histo1D["weight6"] = (TH1F*)reopenSF->Get("weight6");
            histo1D["weight7"] = (TH1F*)reopenSF->Get("weight7");
            histo1D["weight8"] = (TH1F*)reopenSF->Get("weight8");

            // float weighty = weight0->Integral(); cout<<weighty<<endl;

            cout<<"Got weights"<<endl;


            for(int binN = 1; binN <nBins+1; ++binN){
            float binContentMax = -1, binContentMin =  1000000000000000;  

                for(int weightIt = 1; weightIt<9; ++weightIt){
                    if(weightIt == 6 || weightIt == 8){
                        continue; //extreme weights with 1/2u and 2u
                    }

                    string weightStr = static_cast<ostringstream*>( &(ostringstream() << weightIt) )->str();
                    string weightHistoName = ("weight"+weightStr).c_str();
                    cout<<"weight histo name: "<<weightHistoName<<endl;
                    cout<<histo1D[weightHistoName.c_str()]->GetBinContent(binN)<<endl;
                    if(binContentMin > histo1D[weightHistoName]->GetBinContent(binN) ) binContentMin = histo1D[weightHistoName]->GetBinContent(binN);
                    if(binContentMax < histo1D[weightHistoName]->GetBinContent(binN) ) binContentMax = histo1D[weightHistoName]->GetBinContent(binN);
                    cout<<"binContentMax: "<<binContentMax<<"binContentMin: "<<binContentMin<<endl;
                    // cout<<"bin number : "<<binN<<"   bincontent: "<<histo1D["weight0"]->GetBinContent(binN)<<endl;
                }
                histo1D["weightPlus"]->SetBinContent(binN, binContentMax);
                histo1D["weightMinus"]->SetBinContent(binN, binContentMin);
            }  


            errorfile->cd();
            // errorfile->SetWritable();
            // errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
            errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
            histo1D["weightMinus"]->Write("Minus");
            histo1D["weightPlus"]->Write("Plus");
            histo1D["weight0"]->Write("Nominal");
            shapefile->cd();
            string scalesysname = channel+"__TTJets_MLM__scale";
            histo1D["weightMinus"]->Write((scalesysname+"Down").c_str());
            histo1D["weightPlus"]->Write((scalesysname+"Up").c_str());            
        }    
    }
}

void GetScaleEnvelopeSplit(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlNom, string CraneenPath, string shapefileName){

    // const char *xmlfile = xmlNom.c_str();    cout << "used config file: " << xmlfile << endl;

    // ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
    // TTreeLoader treeLoader;
    // vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
    // treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;

    // string plotname;   ///// Jet Split plot
    // string numStr;
    // string sbSplit = static_cast<ostringstream*>( &(ostringstream() << fbSplit) )->str();
    // string stSplit = static_cast<ostringstream*>( &(ostringstream() << ftSplit) )->str();
    // string swSplit = static_cast<ostringstream*>( &(ostringstream() << fwSplit) )->str();
    // for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
    // {
    //     numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
    //     plotname = sVarofinterest + numStr + sSplitVar;
    //     MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
    // }

    // for (int d = 0; d < datasets.size(); d++){  //Loop through datasets
    //     string dataSetName = datasets[d]->Name();
    //     if(dataSetName.find("TTJets")!=string::npos){
    //         cout<<"Producing envelope for scale uncertainty"<<endl;
    //         for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
    //         {
    //             numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
    //             histo1D[("weightPlus"+numStr).c_str()] = new TH1F("Plus","Plus", nBins, plotLow, plotHigh);
    //             histo1D[("weightMinus"+numStr).c_str()] = new TH1F("Minus","Minus", nBins, plotLow, plotHigh);

    //             TFile * reopenSF = new TFile(shapefileName.c_str());
    //             histo1D[("Genweight_tt"+numStr).c_str()] = (TH1F*)reopenSF->Get(("Genweight_tt"+numStr).c_str());
    //             histo1D[("weight1"+numStr).c_str()] = (TH1F*)reopenSF->Get(("weight1"+numStr).c_str());
    //             histo1D[("weight2"+numStr).c_str()] = (TH1F*)reopenSF->Get(("weight2"+numStr).c_str());
    //             histo1D[("weight3"+numStr).c_str()] = (TH1F*)reopenSF->Get(("weight3"+numStr).c_str());
    //             histo1D[("weight4"+numStr).c_str()] = (TH1F*)reopenSF->Get(("weight4"+numStr).c_str());
    //             histo1D[("weight5"+numStr).c_str()] = (TH1F*)reopenSF->Get(("weight5"+numStr).c_str());
    //             histo1D[("weight6"+numStr).c_str()] = (TH1F*)reopenSF->Get(("weight6"+numStr).c_str());
    //             histo1D[("weight7"+numStr).c_str()] = (TH1F*)reopenSF->Get(("weight7"+numStr).c_str());
    //             histo1D[("weight8"+numStr).c_str()] = (TH1F*)reopenSF->Get(("weight8"+numStr).c_str());

    //             cout<<"Got weights"<<endl;


    //             for(int binN = 1; binN <nBins+1; ++binN){
    //             float binContentMax = -1, binContentMin =  1000000000000000;  

    //                 for(int weightIt = 1; weightIt<9; ++weightIt){
    //                     if(weightIt == 6 || weightIt == 8){
    //                         continue; //extreme weights with 1/2u and 2u
    //                     }

    //                     string weightStr = static_cast<ostringstream*>( &(ostringstream() << weightIt) )->str();
    //                     string weightHistoName = ("weight"+weightStr+numStr).c_str();
    //                     if(binContentMin > histo1D[weightHistoName]->GetBinContent(binN) ) binContentMin = histo1D[weightHistoName]->GetBinContent(binN);
    //                     if(binContentMax < histo1D[weightHistoName]->GetBinContent(binN) ) binContentMax = histo1D[weightHistoName]->GetBinContent(binN);
    //                     cout<<"binContentMax: "<<binContentMax<<"binContentMin: "<<binContentMin<<endl;
    //                     // cout<<"bin number : "<<binN<<"   bincontent: "<<histo1D["weight0"]->GetBinContent(binN)<<endl;
    //                 }
    //                 histo1D[("weightPlus"+numStr).c_str()]->SetBinContent(binN, binContentMax);
    //                 histo1D[("weightMinus"+numStr).c_str()]->SetBinContent(binN, binContentMin);
    //             }  


    //             errorfile->cd();
    //             errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
    //             errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
    //             histo1D[("weightMinus"+numStr).c_str()]->Write("Minus");
    //             histo1D[("weightPlus"+numStr).c_str()]->Write("Plus");
    //         }
    //     }    
    // }
}

void DatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string xmlNom, string CraneenPath, string shapefileName)
{
    cout<<""<<endl;    cout<<"RUNNING NOMINAL DATASETS"<<endl;    cout<<""<<endl;
    shapefile->cd();
    errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());

    const char *xmlfile = xmlNom.c_str();    cout << "used config file: " << xmlfile << endl;

    string pathPNG = "FourTop_Light";    pathPNG += leptoAbbr;    pathPNG += "_MSPlots/";
    mkdir(pathPNG.c_str(),0777);    //cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

    ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////
    TTreeLoader treeLoader;
    vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
    treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

    //***************************************************CREATING PLOTS****************************************************
    string plotname = sVarofinterest;   ///// Non Jet Split plot

    //***********************************************OPEN FILES & GET NTUPLES**********************************************
    string dataSetName, filepath;
    int nEntries;
    float ScaleFactor, NormFactor, Luminosity, varofInterest, GenWeight, weight1, weight2, weight3, weight4, weight5, weight6, weight7, weight8, ttbar_flav, SFlepton, SFbtag, SFbtagUp, SFbtagDown, SFPU;
    Dataset* ttbar_ll;
    Dataset* ttbar_cc;
    Dataset* ttbar_bb;
    bool split_ttbar = false;
    if (split_ttbar){
        int ndatasets = datasets.size();
        cout << " - splitting TTBar dataset ..." << ndatasets   << endl;
        vector<string> ttbar_filenames = datasets[ndatasets-1]->Filenames();
        cout <<"ttbar filenames =  "<< ttbar_filenames[0] <<endl;

        ttbar_ll = new Dataset("TTJets_ll", "tt + ll" , true, 633, 2, 2, 1, 213.4, ttbar_filenames );
        ttbar_cc = new Dataset("TTJets_cc", "tt + cc" , true, 633, 2, 2, 1, 6.9,   ttbar_filenames );
        ttbar_bb = new Dataset("TTJets_bb", "tt + bb" , true, 633, 2, 2, 1, 4.8,   ttbar_filenames );

        ///heavy flav re-weight -> scaling ttbb up and ttjj down so that ttbb/ttjj matches CMS measurement.
        double ll_rw = 0.976;
        double bb_rw = 3.;
        int newlumi = datasets[ndatasets-1]->EquivalentLumi();
        ttbar_ll->SetEquivalentLuminosity(newlumi/ll_rw);
        ttbar_cc->SetEquivalentLuminosity(newlumi/ll_rw);
        ttbar_bb->SetEquivalentLuminosity(newlumi/bb_rw);

        ttbar_ll->SetColor(kRed);
        ttbar_cc->SetColor(kRed-3);
        ttbar_bb->SetColor(kRed+2);


        //datasets.pop_back();
        datasets.push_back(ttbar_bb);
        datasets.push_back(ttbar_cc);
        datasets.push_back(ttbar_ll);     
    }     
    MSPlot[plotname] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());

    for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
    {
        dataSetName = datasets[d]->Name();        cout<<"Dataset:  :"<<dataSetName<<endl;
        if (dataSetName.find("TTJets_ll")!=string::npos ||dataSetName.find("TTJets_cc")!=string::npos ||dataSetName.find("TTJets_bb")!=string::npos){
            cout<<"skip ttbar split datasets!"<<endl;
            continue;
        }
        filepath = CraneenPath+dataSetName + "_Run2_TopTree_Study.root";        //cout<<"filepath: "<<filepath<<endl;

        FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset
        string nTuplename = "Craneen__"+ leptoAbbr;
        nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
        nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();        cout<<"                 nEntries: "<<nEntries<<endl;

        nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
        nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
        nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
        nTuple[dataSetName.c_str()]->SetBranchAddress("GenWeight",&GenWeight);
        nTuple[dataSetName.c_str()]->SetBranchAddress("weight1",&weight1);
        nTuple[dataSetName.c_str()]->SetBranchAddress("weight2",&weight2);
        nTuple[dataSetName.c_str()]->SetBranchAddress("weight3",&weight3);
        nTuple[dataSetName.c_str()]->SetBranchAddress("weight4",&weight4);
        nTuple[dataSetName.c_str()]->SetBranchAddress("weight5",&weight5);
        nTuple[dataSetName.c_str()]->SetBranchAddress("weight6",&weight6);
        nTuple[dataSetName.c_str()]->SetBranchAddress("weight7",&weight7);
        nTuple[dataSetName.c_str()]->SetBranchAddress("weight8",&weight8);
        nTuple[dataSetName.c_str()]->SetBranchAddress("SFlepton",&SFlepton);
        nTuple[dataSetName.c_str()]->SetBranchAddress("SFbtag",&SFbtag);
        nTuple[dataSetName.c_str()]->SetBranchAddress("SFbtagUp",&SFbtagUp);
        nTuple[dataSetName.c_str()]->SetBranchAddress("SFbtagDown",&SFbtagDown);
        nTuple[dataSetName.c_str()]->SetBranchAddress("SFPU",&SFPU);       
        if(split_ttbar){nTuple[dataSetName.c_str()]->SetBranchAddress("ttbar_flav",&ttbar_flav);}


        float eqlumi = 1./datasets[d]->EquivalentLumi();
        cout<<"eqlumi: "<<eqlumi<<endl;
        // cout<<"got variables"<<endl;
        //for fixed bin width
        histo1D[dataSetName.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, plotLow, plotHigh);
        if(dataSetName.find("TTJets")!=string::npos && dataSetName.find("JES")==string::npos && dataSetName.find("JER")==string::npos){
            histo1D["Genweight_tt"] = new TH1F("Genweight","Genweight", nBins, plotLow, plotHigh);
            histo1D["weight1_tt"] = new TH1F("weight1","weight1", nBins, plotLow, plotHigh);
            histo1D["weight2_tt"] = new TH1F("weight2","weight2", nBins, plotLow, plotHigh);
            histo1D["weight3_tt"] = new TH1F("weight3","weight3", nBins, plotLow, plotHigh);
            histo1D["weight4_tt"] = new TH1F("weight4","weight4", nBins, plotLow, plotHigh);
            histo1D["weight5_tt"] = new TH1F("weight5","weight5", nBins, plotLow, plotHigh);
            histo1D["weight6_tt"] = new TH1F("weight6","weight6", nBins, plotLow, plotHigh);
            histo1D["weight7_tt"] = new TH1F("weight7","weight7", nBins, plotLow, plotHigh);
            histo1D["weight8_tt"] = new TH1F("weight8","weight8", nBins, plotLow, plotHigh);      
            histo1D["PU_Up"] = new TH1F("PU_Up","PU_Up", nBins, plotLow, plotHigh);
            histo1D["PU_Down"] = new TH1F("PU_Down","PU_Down", nBins, plotLow, plotHigh);
            histo1D["btag_Up"] = new TH1F("btag_Up","btag_Up", nBins, plotLow, plotHigh);
            histo1D["btag_Down"] = new TH1F("btag_Down","btag_Down", nBins, plotLow, plotHigh);    
        }     
        else if(dataSetName.find("tttt")!=string::npos){
            histo1D["Genweight_tttt"] = new TH1F("Genweight_tttt","Genweight_tttt", nBins, plotLow, plotHigh);
        }
        else if(dataSetName.find("JESUp")!=string::npos){
            histo1D["JES_Up"] = new TH1F("JES_Up","JES_Up", nBins, plotLow, plotHigh);
        }
        else if(dataSetName.find("JESDown")!=string::npos){
            histo1D["JES_Down"] = new TH1F("JES_Down","JES_Down", nBins, plotLow, plotHigh);
        }
        else if(dataSetName.find("JERUp")!=string::npos){
            histo1D["JER_Up"] = new TH1F("JER_Up","JER_Up", nBins, plotLow, plotHigh);
        }
        else if(dataSetName.find("JERDown")!=string::npos){
            histo1D["JER_Down"] = new TH1F("JER_Down","JER_Down", nBins, plotLow, plotHigh);
        }

        // cout<<"made th1fs"<<endl;
        /////*****loop through entries and fill plots*****
        for (int j = 0; j<nEntries; j++)
        {
            nTuple[dataSetName.c_str()]->GetEntry(j);
            float OverallSF = SFlepton*SFbtag*SFPU;

            //cout<<"OverallSF: "<<OverallSF<<"   SF: "<<ScaleFactor<<" btag: "<<SFbtag<<" lepton: "<<SFlepton<<" PU: "<<SFPU<<endl;

            // cout<<"genweight: "<<GenWeight<<endl;
            // cout<<"weight1: "<<weight1<<endl;
            if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos || dataSetName.find("NP_overlay_Data")!=string::npos)
            {
                MSPlot[plotname]->Fill(varofInterest, datasets[d], true,  1.);
                histo1D[dataSetName.c_str()]->Fill(varofInterest,eqlumi);
            }
            else
            {
                if(lScale > 0 ) //artificial Lumi
                {
                    Luminosity = 1000*lScale;
                }               
                Luminosity = 2460.38; 
                NormFactor =1;
                if (dataSetName.find("tttt")!=string::npos) {
                    NormFactor = 0.5635;                
                }
                else if (dataSetName.find("WJets")!=string::npos){
                    NormFactor = 0.522796;
                    // cout<<"voi: "<<varofInterest<<"  NormFactor: "<<NormFactor<<"  ScaleFactor: "<<ScaleFactor<<"  lumi: "<<Luminosity<<endl; 
                }
                else if (dataSetName.find("DYJets")!=string::npos){
                    NormFactor = 0.566;
                    // cout<<"voi: "<<varofInterest<<"  NormFactor: "<<NormFactor<<"  ScaleFactor: "<<ScaleFactor<<"  lumi: "<<Luminosity<<endl; 
                }

                // ScaleFactor =1;
                if (dataSetName.find("WJets")!=string::npos || dataSetName.find("DYJets")!=string::npos) cout<<"btagsf: "<<SFbtag<<"  lepton: "<<SFlepton<<"  SF PU: "<<SFPU<<" totalSF:"<<ScaleFactor<<" weight sign: "<<GenWeight<<"  total weight: "<<NormFactor*ScaleFactor*Luminosity*GenWeight*eqlumi<<endl;
                histo1D[dataSetName]->Fill(varofInterest,NormFactor*ScaleFactor*GenWeight*Luminosity*eqlumi);

                if (dataSetName.find("TTJets")!=string::npos && dataSetName.find("JES")==string::npos &&dataSetName.find("JER")==string::npos){

                    histo1D["Genweight_tt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight*eqlumi);
                    histo1D["weight1_tt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight1*eqlumi);
                    histo1D["weight2_tt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight2*eqlumi);
                    histo1D["weight3_tt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight3*eqlumi);
                    histo1D["weight4_tt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight4*eqlumi);
                    histo1D["weight5_tt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight5*eqlumi);
                    histo1D["weight6_tt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight6*eqlumi);
                    histo1D["weight7_tt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight7*eqlumi);
                    histo1D["weight8_tt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight8*eqlumi);   

                    histo1D["PU_Up"]->Fill(varofInterest,NormFactor*SFlepton*SFbtag*Luminosity*GenWeight*eqlumi);   
                    histo1D["PU_Down"]->Fill(varofInterest,NormFactor*SFlepton*SFbtag*SFPU*SFPU*Luminosity*GenWeight*eqlumi);   
                    histo1D["btag_Up"]->Fill(varofInterest,NormFactor*SFlepton*SFbtagUp*SFPU*Luminosity*GenWeight*eqlumi);   
                    histo1D["btag_Down"]->Fill(varofInterest,NormFactor*SFlepton*SFbtagDown*SFPU*Luminosity*GenWeight*eqlumi);             
                     // cout<<"njets: "<<nJets<<endl;
                    if(split_ttbar){
                        // if (ttbar_flav>1.5) cout<<"ttbar_flav:  "<<ttbar_flav<<endl;
                        if(ttbar_flav<0.5){//light
                            
                            MSPlot[plotname]->Fill(varofInterest, ttbar_ll, true, NormFactor*GenWeight*ScaleFactor*Luminosity);                
                        }     
                        else if(ttbar_flav>0.5 && ttbar_flav<1.5){//cc
                            MSPlot[plotname]->Fill(varofInterest, ttbar_cc, true, NormFactor*GenWeight*ScaleFactor*Luminosity);                
                        }   
                        else{
                            MSPlot[plotname]->Fill(varofInterest, ttbar_bb, true, NormFactor*GenWeight*ScaleFactor*Luminosity);                
                        }          
                    }
                    else{
                        MSPlot[plotname]->Fill(varofInterest, datasets[d], true, NormFactor*GenWeight*ScaleFactor*Luminosity);   
                    }             
                }
                else if (dataSetName.find("tttt")!=string::npos){
                    // histo1D["Genweight_tttt"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight*eqlumi);
                    MSPlot[plotname]->Fill(varofInterest, datasets[d], true, NormFactor*GenWeight*ScaleFactor*Luminosity*GenWeight);
                }
                else if(dataSetName.find("JESUp")!=string::npos){
                    histo1D["JES_Up"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight*eqlumi);
                }
                else if(dataSetName.find("JESDown")!=string::npos){
                    histo1D["JES_Down"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight*eqlumi);
                }
                else if(dataSetName.find("JERUp")!=string::npos){
                    histo1D["JER_Up"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight*eqlumi);
                }
                else if(dataSetName.find("JERDown")!=string::npos){
                    histo1D["JER_Down"]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight*eqlumi);
                }            
                else{
                    MSPlot[plotname]->Fill(varofInterest, datasets[d], true, NormFactor*GenWeight*ScaleFactor*Luminosity*GenWeight);
                }
            }
        }//end of event loop
        // cout<<"after for loop entries"<<endl;
        //if scaleenv
        // if(dataSetName.find("TTJets")!=string::npos)  //to put nominal histo into error file
        // {
        //     errorfile->cd();
        //     errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
        //     histo1D["Genweight_tt"]->Write("Nominal");
        //     errorfile->Write();
        // }

        shapefile->cd();
        TCanvas *canv = new TCanvas();
        histo1D[dataSetName.c_str()]->Draw();
        string writename = "";
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos){
            writename = channel + "__data_obs__nominal";
        }
        else if(dataSetName.find("JESUp")!=string::npos){
            writename = channel + "__TTJets_MLM__JESUp";
        }
        else if(dataSetName.find("JESDown")!=string::npos){
            writename = channel + "__TTJets_MLM__JESDown";
        }   
        else{
            writename = channel + "__" + dataSetName +"__nominal";
        }
        // cout<<"writename  :"<<writename<<endl;
        histo1D[dataSetName.c_str()]->Write((writename).c_str());
        canv->SaveAs((pathPNG+dataSetName+".png").c_str());
        // cout<<"datasetname:   "<<dataSetName<<endl;
        /////
        string scalesysname = channel+"__TTJets_MLM__";

        if ( dataSetName.find("TTJets")!= string::npos && dataSetName.find("JES")==string::npos &&dataSetName.find("JER")==string::npos){
            //cout<<"  making weights histos"<<endl;
       
            TCanvas *canv0 = new TCanvas();
            histo1D["Genweight_tt"]->Draw();
            histo1D["Genweight_tt"]->Write("Genweight_tt");                  
            canv0->SaveAs("Genweight_tt.png");

            for(int ii=1; ii<9; ii++){
                TCanvas *canv1 = new TCanvas();
                string weightstring = static_cast<ostringstream*>( &(ostringstream() << ii) )->str();
                string weighthisto = "weight"+weightstring+"_tt";
                //cout<<"weight histo:  "<<weighthisto<<endl;
                histo1D[weighthisto]->Draw();
                histo1D[weighthisto]->Write(("weight"+weightstring).c_str());
                canv1->SaveAs(("weight"+weightstring+".png").c_str());
            }    

            histo1D["PU_Up"]->Write((scalesysname+"PUUP").c_str());   
            histo1D["PU_Down"]->Write((scalesysname+"PUDown").c_str());   
            histo1D["btag_Up"]->Write((scalesysname+"btagUp").c_str());   
            histo1D["btag_Down"]->Write((scalesysname+"btagDown").c_str());  
        }
        else if(dataSetName.find("JESUp")!=string::npos){
            histo1D["JES_Up"]->Write((scalesysname+"JESUp").c_str());
        }
        else if(dataSetName.find("JESDown")!=string::npos){
            histo1D["JES_Down"]->Write((scalesysname+"JESDown").c_str());
        }
        else if(dataSetName.find("JERUp")!=string::npos){
            histo1D["JER_Up"]->Write((scalesysname+"JERUp").c_str());
        }
        else if(dataSetName.find("JERDown")!=string::npos){
            histo1D["JER_Down"]->Write((scalesysname+"JERDown").c_str());
        }  
        else if ( sVarofinterest==DatacardVar && dataSetName.find("tttt")!=string::npos ){
            TCanvas *canv2 = new TCanvas();
            histo1D["Genweight_tttt"]->Draw();
            histo1D["Genweight_tttt"]->Write("Genweight_tttt");                  
            canv2->SaveAs("Genweight_tttt.png");            
        }
        // cout<<"saved histos"<<endl;
    } //end of dataset loop
    //cout<<"unloading datasets"<<endl;

    GetScaleEnvelope(nBins, lScale, plotLow, plotHigh, leptoAbbr, shapefile, errorfile, channel, sVarofinterest, xmlNom, CraneenPath, shapefileName);
    errorfile->Close();

    // treeLoader.UnLoadDataset();
    // cout<<"make scale dirs"<<endl;
    string scaleFileDir = "ScaleFiles" + leptoAbbr + "_light";
    // mkdir(scaleFileDir.c_str(),0777);
    string scaleFileName = scaleFileDir + "/Error"+sVarofinterest+".root";
    // cout<<"plotname: "<<plotname<<endl;
    // MSPlot[plotname.c_str()]->Draw();
    // MSPlot[plotname.c_str()]->setErrorBandFile("ScaleFilesMu_light/ErrorBDTSave2.root"); //set error file for uncertainty bands on multisample plot
    MSPlot[plotname.c_str()]->setErrorBandFile((scaleFileName).c_str()); //set error file for uncertainty bands on multisample plot
    MSPlot[plotname.c_str()]->setDataLumi(2500);
    MSPlot[plotname.c_str()]->setMaxY(1000000);
    cout<<scaleFileName<<endl;
    // cout<<"running delete loop"<<endl;
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;  
        MultiSamplePlot *temp = it->second;
        temp->Draw(sVarofinterest.c_str(), 1, true, true, true, 100);
        temp->Write(shapefile, name, true, pathPNG, "png");
    }
    MSPlot.erase(plotname);

    for(map<string,TH1F*>::const_iterator it = histo1D.begin(); it!=histo1D.end();it++){
        histo1D.erase(it->first);
    }
};


// void SystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile* errorfile, string channel, string sVarofinterest, string xmlSys, string CraneenPath)
// {
//     // cout<<""<<endl;
//     // cout<<"RUNNING SYS"<<endl;
//     // cout<<""<<endl;
//     // const char *xmlfile = xmlSys.c_str();
//     // cout << "used config file: " << xmlfile << endl;

//     // string pathPNG = "FourTop_SysPlots_" + leptoAbbr; 	//add MSplot name to directory
//     // mkdir(pathPNG.c_str(),0777);
//     // cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

//     // ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
//     // TTreeLoader treeLoader;
//     // vector < Dataset* > datasets; //cout<<"vector filled"<<endl;
//     // treeLoader.LoadDatasets (datasets, xmlfile);
//     // cout<<"datasets loaded"<<endl;

//     // ///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
//     // string dataSetName, filepath;
//     // int nEntries;
//     // float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH;

//     // for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
//     // {
//     //     dataSetName = datasets[d]->Name();
//     //     cout<<"	Dataset:  :"<<dataSetName<<endl;
//     //     filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study"+".root";

//     //     cout<<""<<endl;
//     //     cout<<"file in use: "<<filepath<<endl;

//     //     FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());
//     //     cout<<"initialised file"<<endl;
//     //     string nTupleName = "Craneen__"+leptoAbbr;

//     //     nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
//     //     nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();
//     //     cout<<"                 nEntries: "<<nEntries<<endl;

//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);



//     //     ////****************************************************Define plots***********************************************
//     //     string plotname = sVarofinterest+"_"+dataSetName;
//     //     histo1D[plotname.c_str()] = new TH1F(dataSetName.c_str(),dataSetName.c_str(), nBins, plotLow, plotHigh);

//     //     for (int i = 0; i<nEntries; i++)   //Fill histo with variable of interest
//     //     {
//     //         nTuple[dataSetName.c_str()]->GetEntry(i);
//     //         //artificial Lumi
//     //         if(lScale > 0 )
//     //         {
//     //             Luminosity = 1000*lScale;
//     //         }
//     //         histo1D[plotname.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
//     //     }

//     //     ////****************************************************Fill plots***********************************************
//     //     shapefile->cd();

//     //     TCanvas *canv2 = new TCanvas();
//     //     if(Normalise)
//     //     {
//     //         double dIntegral = histo1D[plotname.c_str()]->Integral();
//     //         histo1D[plotname.c_str()]->Scale(1./dIntegral);
//     //     }
//     //     histo1D[plotname.c_str()]->Draw();
//     //     string writename = "";
//     //     writename = channel + "__TTJets__" + dataSetName;
//     //     //cout<<"writename  :"<<writename<<endl;

//     //     histo1D[plotname.c_str()]->Write((writename).c_str());
//     //     canv2->SaveAs(("Sys_"+plotname+".png").c_str());

//     //     if(dataSetName == "TTScaleDown")
//     //     {
//     //         errorfile->cd();
//     //         errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //         errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //         histo1D[plotname.c_str()]->Write("Minus");
//     //         //errorfile->Write();
//     //     }

//     //     if(dataSetName == "TTScaleUp")
//     //     {
//     //         errorfile->cd();
//     //         errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //         histo1D[plotname.c_str()]->Write("Plus");
//     //         //errorfile->Write();
//     //     }
//     // }
// };

// void SplitDatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlNom, string CraneenPath)
// {
//     // cout<<""<<endl;
//     // cout<<"RUNNING NOMINAL DATASETS"<<endl;
//     // cout<<""<<endl;
//     // shapefile->cd();

//     // const char *xmlfile = xmlNom.c_str();
//     // cout << "used config file: " << xmlfile << endl;

//     // string pathPNG = "FourTop_Light";
//     // pathPNG += leptoAbbr;
//     // pathPNG += "_MSPlots/";
//     // mkdir(pathPNG.c_str(),0777);
//     // cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

//     // ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
//     // TTreeLoader treeLoader;
//     // vector < Dataset* > datasets; 					//cout<<"vector filled"<<endl;
//     // treeLoader.LoadDatasets (datasets, xmlfile);	//cout<<"datasets loaded"<<endl;

//     // //***************************************************CREATING PLOTS****************************************************
//     // string plotname;   ///// Jet Split plot
//     // string numStr;
//     // string sbSplit = static_cast<ostringstream*>( &(ostringstream() << fbSplit) )->str();
//     // string stSplit = static_cast<ostringstream*>( &(ostringstream() << ftSplit) )->str();
//     // string swSplit = static_cast<ostringstream*>( &(ostringstream() << fwSplit) )->str();
//     // for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
//     // {
//     //     numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     plotname = sVarofinterest + numStr + sSplitVar;
//     //     MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
//     // }
//     // plotname = "";

//     // //***********************************************OPEN FILES & GET NTUPLES**********************************************
//     // string dataSetName, filepath, histoName;
//     // int nEntries;
//     // float ScaleFactor, NormFactor, Luminosity, varofInterest, splitVar, GenWeight, weight1, weight2, weight3, weight4, weight5, weight6, weight7, weight8;

//     // for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
//     // {
//     //     dataSetName = datasets[d]->Name();
//     //     cout<<"Dataset:  :"<<dataSetName<<endl;

//     //     filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study.root";
//     //     //cout<<"filepath: "<<filepath<<endl;

//     //     FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset
//     //     string nTuplename = "Craneen__"+ leptoAbbr;
//     //     nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
//     //     nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();
//     //     cout<<"                 nEntries: "<<nEntries<<endl;

//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar.c_str(), &splitVar);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("GenWeight",&GenWeight);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("weight1",&weight1);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("weight2",&weight2);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("weight3",&weight3);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("weight4",&weight4);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("weight5",&weight5);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("weight6",&weight6);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("weight7",&weight7);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("weight8",&weight8);


//     //     //for fixed bin width
//     //     for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
//     //     {
//     //         numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //         histoName = dataSetName + numStr + sSplitVar;
//     //         histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);

//     //         if(dataSetName.find("TTJets")!=string::npos){
//     //             histo1D[("Genweight_tt"+numStr).c_str()] = new TH1F(("Genweight_tt"+numStr).c_str(),"Genweight", nBins, plotLow, plotHigh);
//     //             histo1D[("weight1_tt"+numStr).c_str()] = new TH1F(("weight1_tt"+numStr).c_str(),"weight1", nBins, plotLow, plotHigh);
//     //             histo1D[("weight2_tt"+numStr).c_str()] = new TH1F(("weight2_tt"+numStr).c_str(),"weight2", nBins, plotLow, plotHigh);
//     //             histo1D[("weight3_tt"+numStr).c_str()] = new TH1F(("weight3_tt"+numStr).c_str(),"weight3", nBins, plotLow, plotHigh);
//     //             histo1D[("weight4_tt"+numStr).c_str()] = new TH1F(("weight4_tt"+numStr).c_str(),"weight4", nBins, plotLow, plotHigh);
//     //             histo1D[("weight5_tt"+numStr).c_str()] = new TH1F(("weight5_tt"+numStr).c_str(),"weight5", nBins, plotLow, plotHigh);
//     //             histo1D[("weight6_tt"+numStr).c_str()] = new TH1F(("weight6_tt"+numStr).c_str(),"weight6", nBins, plotLow, plotHigh);
//     //             histo1D[("weight7_tt"+numStr).c_str()] = new TH1F(("weight7_tt"+numStr).c_str(),"weight7", nBins, plotLow, plotHigh);
//     //             histo1D[("weight8_tt"+numStr).c_str()] = new TH1F(("weight8_tt"+numStr).c_str(),"weight8", nBins, plotLow, plotHigh);            
//     //         }
//     //         else if(dataSetName.find("tttt")!=string::npos){
//     //             histo1D[("Genweight_tttt"+numStr).c_str()] = new TH1F(("Genweight_tttt"+numStr).c_str(),"Genweight_tttt", nBins, plotLow, plotHigh);
//     //         }

//     //     }
//     //     /////*****loop through entries and fill plots*****
//     //     for (int j = 0; j<nEntries; j++)
//     //     {
//     //         nTuple[dataSetName.c_str()]->GetEntry(j);
//     //         //artificial Lumi

//     //         if(lScale > 0 )
//     //         {
//     //             Luminosity = 1000*lScale;
//     //         }

            
//     //         if(splitVar >= ftSplit) //Check if this entry belongs in the last bin.  Done here to optimize number of checks
//     //         {
//     //             plotname = sVarofinterest + stSplit + sSplitVar;
//     //             histoName = dataSetName + stSplit + sSplitVar;
//     //             //outsiderange = false;
//     //             if (dataSetName == "TTJets"){
//     //                 histo1D[("Genweight_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight);
//     //                 //cout<<"genweight: "<<GenWeight<<endl;
//     //                 histo1D[("weight1_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight1);
//     //                 histo1D[("weight2_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight2);
//     //                 histo1D[("weight3_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight3);
//     //                 histo1D[("weight4_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight4);
//     //                 histo1D[("weight5_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight5);
//     //                 histo1D[("weight6_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight6);
//     //                 histo1D[("weight7_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight7);
//     //                 histo1D[("weight8_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight8);                    
//     //             }
//     //             else if (dataSetName.find("tttt")!=string::npos){
//     //                 histo1D[("Genweight_tttt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight);
//     //             }

//     //         }
//     //         else //If it doesn't belong in the last bin, find out which it belongs in
//     //         {
//     //             bool outsiderange = true;
//     //             for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
//     //             {
//     //                 if(splitVar>=s && splitVar<(s+fwSplit)) //splitVar falls inside one of the bins
//     //                 {
//     //                     numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //                     plotname = sVarofinterest + numStr + sSplitVar;
//     //                     histoName = dataSetName + numStr + sSplitVar;
//     //                     if (dataSetName == "TTJets"){
//     //                         histo1D[("Genweight_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight);
//     //                         //cout<<"genweight: "<<GenWeight<<endl;
//     //                         histo1D[("weight1_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight1);
//     //                         histo1D[("weight2_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight2);
//     //                         histo1D[("weight3_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight3);
//     //                         histo1D[("weight4_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight4);
//     //                         histo1D[("weight5_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight5);
//     //                         histo1D[("weight6_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight6);
//     //                         histo1D[("weight7_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight7);
//     //                         histo1D[("weight8_tt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*weight8);                    
//     //                     }
//     //                     else if (dataSetName.find("tttt")!=string::npos){
//     //                         histo1D[("Genweight_tttt"+numStr).c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity*GenWeight);
//     //                     }

//     //                     outsiderange = false;
//     //                     break;
//     //                 }
//     //             }
//     //             if (outsiderange == true){
//     //                 continue;
//     //             }
//     //         }


//     //         if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //         {
//     //             MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, NormFactor*ScaleFactor*Luminosity);
//     //             histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
//     //         }
//     //         else
//     //         {
//     //             MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, ScaleFactor*Luminosity);
//     //             histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
//     //         }

//     //     }
//     //     if(dataSetName == "TTJets")  //to put nominal histo into error file
//     //     {
//     //         errorfile->cd();
//     //         errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //         histo1D[histoName.c_str()]->Write("Nominal");
//     //         errorfile->Write();
//     //     }

//     //     shapefile->cd();
//     //     TCanvas *canv = new TCanvas();
//     //     for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
//     //     {
//     //         numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //         plotname = sVarofinterest + numStr + sSplitVar;
//     //         histoName = dataSetName + numStr + sSplitVar;
//     //         histo1D[histoName.c_str()]->Draw();
//     //         string writename = "";
//     //         if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //         {
//     //             writename = channel + numStr + sSplitVar + "__data_obs__nominal";
//     //         }
//     //         else
//     //         {
//     //             writename = channel + numStr + sSplitVar + "__" + dataSetName +"__nominal";
//     //         }
//     //         //cout<<"writename  :"<<writename<<endl;
//     //         histo1D[histoName.c_str()]->Write((writename).c_str());

//     //         canv->SaveAs((pathPNG+histoName+".png").c_str());
//     //         if ( sVarofinterest==DatacardVar && dataSetName.find("TTJets")!= string::npos ){
//     //             cout<<"  making weights histos"<<endl;
           
//     //             TCanvas *canv0 = new TCanvas();
//     //             histo1D[("Genweight_tt"+numStr).c_str()]->Draw();
//     //             histo1D[("Genweight_tt"+numStr).c_str()]->Write(("Genweight_tt"+numStr).c_str());                  
//     //             canv0->SaveAs(("Genweight_tt"+numStr+".png").c_str());

//     //             for(int ii=1; ii<9; ii++){
//     //                 TCanvas *canv1 = new TCanvas();
//     //                 string weightstring = static_cast<ostringstream*>( &(ostringstream() << ii) )->str();
//     //                 string weighthisto = "weight"+weightstring+"_tt";
//     //                 histo1D[weighthisto]->Draw();
//     //                 histo1D[weighthisto]->Write(("weight"+weightstring).c_str());
//     //                 canv1->SaveAs(("weight"+weightstring+".png").c_str());

//     //             }    
       
//     //         }
//     //         else if ( sVarofinterest==DatacardVar && dataSetName.find("tttt")!=string::npos ){
//     //             TCanvas *canv2 = new TCanvas();
//     //             histo1D[("Genweight_tttt"+numStr).c_str()]->Draw();
//     //             histo1D[("Genweight_tttt"+numStr).c_str()]->Write(("Genweight_tttt"+numStr).c_str());                  
//     //             canv2->SaveAs(("Genweight_tttt"+numStr+".png").c_str());            
//     //         }
//     //     }
//     // }


//     // treeLoader.UnLoadDataset();

//     // string scaleFileDir = "ScaleFiles" + leptoAbbr + "_light";
//     // mkdir(scaleFileDir.c_str(),0777);
//     // string scaleFileName = scaleFileDir + "/Error"+sVarofinterest+".root";
//     // //MSPlot[plotname.c_str()]->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot

//     // for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
//     // {
//     //     string name = it->first;  cout<<name<< " ** "<<sVarofinterest<<endl;
//     //     MultiSamplePlot *temp = it->second;
//     //     temp->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot
//     //     temp->Draw(sVarofinterest.c_str(), 0, false, false, false, 100);
//     //     temp->Write(shapefile, name, true, pathPNG, "png");
//     // }

//     // for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
//     // {
//     //     numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     plotname = sVarofinterest + numStr + sSplitVar;
//     //     MSPlot.erase(plotname);
//     // }

// };


// void SplitSystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar, float fbSplit, float ftSplit, float fwSplit, string xmlSys, string CraneenPath)
// {
//     // cout<<""<<endl;
//     // cout<<"RUNNING SYS"<<endl;
//     // cout<<""<<endl;
//     // const char *xmlfile = xmlSys.c_str();
//     // cout << "used config file: " << xmlfile << endl;

//     // string pathPNG = "FourTop_SysPlots_" + leptoAbbr; 	//add MSplot name to directory
//     // mkdir(pathPNG.c_str(),0777);
//     // cout <<"Making directory :"<< pathPNG  <<endl;		//make directory

//     // ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
//     // TTreeLoader treeLoader;
//     // vector < Dataset* > datasets; //cout<<"vector filled"<<endl;
//     // treeLoader.LoadDatasets (datasets, xmlfile);
//     // cout<<"datasets loaded"<<endl;

//     // ///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
//     // string dataSetName, filepath, numStr, histoName;
//     // int nEntries;
//     // float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH, splitVar;
//     // string sbSplit = static_cast<ostringstream*>( &(ostringstream() << fbSplit) )->str();
//     // string stSplit = static_cast<ostringstream*>( &(ostringstream() << ftSplit) )->str();
//     // string swSplit = static_cast<ostringstream*>( &(ostringstream() << fwSplit) )->str();
//     // errorfile->cd();
//     // errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
//     // for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
//     // {
//     //     dataSetName = datasets[d]->Name();
//     //     cout<<"	Dataset:  :"<<dataSetName<<endl;
//     //     filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study"+".root";

//     //     cout<<""<<endl;
//     //     cout<<"file in use: "<<filepath<<endl;

//     //     FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());
//     //     cout<<"initialised file"<<endl;
//     //     string nTupleName = "Craneen__"+leptoAbbr;

//     //     nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
//     //     nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();
//     //     cout<<"                 nEntries: "<<nEntries<<endl;

//     //     //sVarofinterest = "HT";
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar.c_str(), &splitVar);



//     //     ////****************************************************Define plots***********************************************
//     //     for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
//     //     {
//     //         numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //         histoName = dataSetName + numStr + sSplitVar;
//     //         histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);
//     //     }

//     //     for (int i = 0; i<nEntries; i++)   //Fill histo with variable of interest
//     //     {
//     //         nTuple[dataSetName.c_str()]->GetEntry(i);
//     //         //artificial Lumi

//     //         if(lScale > 0 )
//     //         {
//     //             Luminosity = 1000*lScale;
//     //         }
            
//     //         //bool outsiderange = true;
//     //         if(splitVar >= ftSplit) //Check if this entry belongs in the last bin.  Done here to optimize number of checks
//     //         {
//     //             histoName = dataSetName + stSplit + sSplitVar;
//     //             //outsiderange = false;
//     //         }
//     //         else //If it doesn't belong in the last bin, find out which it belongs in
//     //         {
//     //             bool outsiderange = true;
//     //             for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
//     //             {
//     //                 if(splitVar>=s && splitVar<(s+fwSplit)) //splitVar falls inside one of the bins
//     //                 {
//     //                     numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //                     histoName = dataSetName + numStr + sSplitVar;
//     //                     outsiderange = false;
//     //                     break;
//     //                 }
//     //             }
//     //             if (outsiderange == true){
//     //                 continue;
//     //             }
//     //         }


//     //         histo1D[histoName.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
//     //     }

//     //     ////****************************************************Fill plots***********************************************


//     //     TCanvas *canv2 = new TCanvas();

//     //     for(int s = fbSplit; s <= ftSplit; s+=fwSplit)
//     //     {
//     //         shapefile->cd();

//     //         numStr = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //         histoName = dataSetName + numStr + sSplitVar;
//     //         if(Normalise)
//     //         {
//     //             double dIntegral = histo1D[histoName.c_str()]->Integral();
//     //             histo1D[histoName.c_str()]->Scale(1./dIntegral);
//     //         }
//     //         histo1D[histoName.c_str()]->Draw();
//     //         string writename = "";

//     //         if(dataSetName == "TTScaleDown")
//     //         {
//     //             writename = channel + numStr + sSplitVar + "__TTJets__scaleDown";
//     //             //cout<<"writename  :"<<writename<<endl;
//     //             histo1D[histoName.c_str()]->Write((writename).c_str());
//     //             canv2->SaveAs(("Sys_"+histoName+".png").c_str());
//     //             errorfile->cd();
//     //             //errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //             errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //             histo1D[histoName.c_str()]->Write("Minus");
//     //             //errorfile->Write();
//     //         }

//     //         if(dataSetName == "TTScaleUp")
//     //         {
//     //             writename = channel + numStr + sSplitVar + "__TTJets__scaleUp";
//     //             //cout<<"writename  :"<<writename<<endl;
//     //             histo1D[histoName.c_str()]->Write((writename).c_str());
//     //             canv2->SaveAs(("Sys_"+histoName+".png").c_str());
//     //             errorfile->cd();
//     //             errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //             histo1D[histoName.c_str()]->Write("Plus");
//     //             //errorfile->Write();
//     //         }
//     //     }
//     // }
// };

// void Split2DatasetPlotter(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, TFile *shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom, string CraneenPath)
// {
//     // cout<<""<<endl;
//     // cout<<"RUNNING NOMINAL DATASETS"<<endl;
//     // cout<<""<<endl;
//     // shapefile->cd();

//     // const char *xmlfile = xmlNom.c_str();
//     // cout << "used config file: " << xmlfile << endl;

//     // string pathPNG = "FourTop_Light";
//     // pathPNG += leptoAbbr;
//     // pathPNG += "_MSPlots/";
//     // mkdir(pathPNG.c_str(),0777);
//     // cout <<"Making directory :"<< pathPNG  <<endl;      //make directory

//     // ///////////////////////////////////////////////////////////// Load Datasets //////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
//     // TTreeLoader treeLoader;
//     // vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
//     // treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;

//     // //***************************************************CREATING PLOTS****************************************************
//     // string plotname;   ///// Jet Split plot
//     // string numStr1;
//     // string numStr2;
//     // string sbSplit1 = static_cast<ostringstream*>( &(ostringstream() << fbSplit1) )->str();
//     // string stSplit1 = static_cast<ostringstream*>( &(ostringstream() << ftSplit1) )->str();
//     // string swSplit1 = static_cast<ostringstream*>( &(ostringstream() << fwSplit1) )->str();
//     // string sbSplit2 = static_cast<ostringstream*>( &(ostringstream() << fbSplit2) )->str();
//     // string stSplit2 = static_cast<ostringstream*>( &(ostringstream() << ftSplit2) )->str();
//     // string swSplit2 = static_cast<ostringstream*>( &(ostringstream() << fwSplit2) )->str();
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //         numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //         plotname = sVarofinterest + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //         MSPlot[plotname.c_str()] = new MultiSamplePlot(datasets, plotname.c_str(), nBins, plotLow, plotHigh, sVarofinterest.c_str());
//     //     }  
//     // }
//     // plotname = "";

//     // //***********************************************OPEN FILES & GET NTUPLES**********************************************
//     // string dataSetName, filepath, histoName;
//     // int nEntries;
//     // float ScaleFactor, NormFactor, Luminosity, varofInterest, splitVar1, splitVar2;

//     // for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
//     // {
//     //     dataSetName = datasets[d]->Name();
//     //     cout<<"Dataset:  :"<<dataSetName<<endl;

//     //     filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study.root";
//     //     //cout<<"filepath: "<<filepath<<endl;

//     //     FileObj[dataSetName.c_str()] = new TFile((filepath).c_str()); //create TFile for each dataset
//     //     string nTuplename = "Craneen__"+ leptoAbbr;
//     //     nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTuplename.c_str()); //get ntuple for each dataset
//     //     nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();
//     //     cout<<"                 nEntries: "<<nEntries<<endl;

//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar1.c_str(), &splitVar1);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar2.c_str(), &splitVar2);

//     //     //for fixed bin width
//     //     for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     //     {
//     //         numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //         for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2){
//     //             numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();   
//     //             histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //             histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);
//     //         }
//     //     }
//     //     /////*****loop through entries and fill plots*****
//     //     for (int j = 0; j<nEntries; j++)
//     //     {
//     //         nTuple[dataSetName.c_str()]->GetEntry(j);
//     //         //artificial Lumi
//     //         if(lScale > 0 )
//     //         {
//     //             Luminosity = 1000*lScale;
//     //         }
//     //         bool outsiderange = true;
//     //         if(splitVar1 >= ftSplit1) //Check if this entry belongs in the last bin in var1.  Done here to optimize number of checks
//     //         {   

//     //             if(splitVar2 >= ftSplit2){ //Check if this entry belongs in the last bin in var2.
//     //                 plotname = sVarofinterest + stSplit1 + sSplitVar1 +stSplit2 + sSplitVar2;
//     //                 histoName = dataSetName + stSplit1 + sSplitVar1 +stSplit2 + sSplitVar2;
//     //                 outsiderange = false;
//     //                 //                                cout<<"splitvar2: "<<splitVar2<<endl;
//     //             }    
//     //             else //If it doesn't belong in the last bin in var2, find out which it belongs in
//     //             {

//     //                 for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2)
//     //                 {
//     //                     if(splitVar2>=t2 && splitVar2<(t2+fwSplit2)) //splitVar falls inside one of the bins
//     //                     {
//     //                         numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //                         plotname = sVarofinterest + stSplit1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //                         histoName = dataSetName + stSplit1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //                         outsiderange = false;
//     //                         break;
//     //                     }
//     //                 }
//     //             }            
//     //         }
//     //         else //If it doesn't belong in the last bin, find out which it belongs in
//     //         {
//     //             for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     //             {
//     //                 if(splitVar1>=s && splitVar1<(s+fwSplit1)) //splitVar falls inside one of the bins
//     //                 {
//     //                     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //                     if(splitVar2 >= ftSplit2){ //Check if this entry belongs in the last bin in var2.
//     //                         plotname = sVarofinterest + numStr1 + sSplitVar1 +stSplit2 + sSplitVar2;
//     //                         histoName = dataSetName + numStr1 + sSplitVar1 +stSplit2 + sSplitVar2;
//     //                         outsiderange = false;
//     //                     }    
//     //                     else //If it doesn't belong in the last bin, find out which it belongs in
//     //                     {
//     //                         for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2)
//     //                         {
//     //                             if(splitVar2>=t2 && splitVar2<(t2+fwSplit2)) //splitVar falls inside one of the bins
//     //                             {
//     //                                 numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //                                 plotname = sVarofinterest + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //                                 histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //                                 outsiderange = false;
//     //                                 break;
//     //                             }
//     //                         }
//     //                     } 
//     //                 }
//     //             }
//     //         }
//     //         if (outsiderange == true){
//     //             continue;
//     //         }
//     //         if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //         {
//     //             MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, NormFactor*ScaleFactor*Luminosity);
//     //             histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
//     //         }
//     //         else
//     //         {
//     //             MSPlot[plotname.c_str()]->Fill(varofInterest, datasets[d], true, ScaleFactor*Luminosity);
//     //             histo1D[histoName.c_str()]->Fill(varofInterest,NormFactor*ScaleFactor*Luminosity);
//     //         }
//     //     }
//     //     if(dataSetName == "TTJets")  //to put nominal histo into error file
//     //     {
//     //         errorfile->cd();
//     //         errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //         histo1D[histoName.c_str()]->Write("Nominal");
//     //         errorfile->Write();
//     //     }

//     //     shapefile->cd();
//     //     TCanvas *canv = new TCanvas();
//     //     for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     //     {
//     //         numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //         for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2){
//     //             numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //             plotname = sVarofinterest + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //             histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //             histo1D[histoName.c_str()]->Draw();
//     //             string writename = "";
//     //             if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //             {
//     //                 writename = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__data_obs__nominal";
//     //             }
//     //             else
//     //             {
//     //                 writename = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__" + dataSetName +"__nominal";
//     //             }
//     //             //cout<<"writename  :"<<writename<<endl;
//     //             histo1D[histoName.c_str()]->Write((writename).c_str());

//     //             canv->SaveAs((pathPNG+histoName+".png").c_str());
//     //         }
//     //     }
//     // }


//     // treeLoader.UnLoadDataset();

//     // string scaleFileDir = "ScaleFiles" + leptoAbbr + "_light";
//     // mkdir(scaleFileDir.c_str(),0777);
//     // string scaleFileName = scaleFileDir + "/Error"+sVarofinterest+".root";
//     // //MSPlot[plotname.c_str()]->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot

//     // for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
//     // {
//     //     string name = it->first;
//     //     MultiSamplePlot *temp = it->second;
//     //     temp->setErrorBandFile(scaleFileName.c_str()); //set error file for uncertainty bands on multisample plot
//     //     temp->Draw(sVarofinterest.c_str(), 0, false, false, false, 1);
//     //     temp->Write(shapefile, name, true, pathPNG, "png");
//     // }
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //         numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //         plotname = sVarofinterest + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //         MSPlot.erase(plotname);
//     //     }  
//     // }
// };


// void Split2SystematicsAnalyser(int nBins, float lScale, float plotLow, float plotHigh, string leptoAbbr, bool Normalise, TFile* shapefile, TFile *errorfile, string channel, string sVarofinterest, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlSys, string CraneenPath)
// {
//     // cout<<""<<endl;
//     // cout<<"RUNNING SYS"<<endl;
//     // cout<<""<<endl;
//     // const char *xmlfile = xmlSys.c_str();
//     // cout << "used config file: " << xmlfile << endl;

//     // string pathPNG = "FourTop_SysPlots_" + leptoAbbr;   //add MSplot name to directory
//     // mkdir(pathPNG.c_str(),0777);
//     // cout <<"Making directory :"<< pathPNG  <<endl;      //make directory

//     // ///////////////////////////////////////////////////////////// Load Datasets ////////////////////////////////////////////////////////////////////cout<<"loading...."<<endl;
//     // TTreeLoader treeLoader;
//     // vector < Dataset* > datasets; //cout<<"vector filled"<<endl;
//     // treeLoader.LoadDatasets (datasets, xmlfile);
//     // cout<<"datasets loaded"<<endl;

//     // ///////////////////////////Open files and get ntuples//////////////////////////////////////////////////
//     // string dataSetName, filepath, numStr1, numStr2, histoName;
//     // int nEntries;
//     // float varofInterest, ScaleFactor, NormFactor, Luminosity, Njets, Ntags, HTb, HTX, HTH, splitVar1, splitVar2;
//     // string sbSplit1 = static_cast<ostringstream*>( &(ostringstream() << fbSplit1) )->str();
//     // string stSplit1 = static_cast<ostringstream*>( &(ostringstream() << ftSplit1) )->str();
//     // string swSplit1 = static_cast<ostringstream*>( &(ostringstream() << fwSplit1) )->str();
//     // string sbSplit2 = static_cast<ostringstream*>( &(ostringstream() << fbSplit2) )->str();
//     // string stSplit2 = static_cast<ostringstream*>( &(ostringstream() << ftSplit2) )->str();
//     // string swSplit2 = static_cast<ostringstream*>( &(ostringstream() << fwSplit2) )->str();

//     // errorfile->cd();
//     // errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
//     // for (int d = 0; d < datasets.size(); d++)  //Loop through datasets
//     // {
//     //     dataSetName = datasets[d]->Name();
//     //     cout<<" Dataset:  :"<<dataSetName<<endl;
//     //     filepath = CraneenPath + dataSetName + "_Run2_TopTree_Study"+".root";

//     //     cout<<""<<endl;
//     //     cout<<"file in use: "<<filepath<<endl;

//     //     FileObj[dataSetName.c_str()] = new TFile((filepath).c_str());
//     //     cout<<"initialised file"<<endl;
//     //     string nTupleName = "Craneen__"+leptoAbbr;

//     //     nTuple[dataSetName.c_str()] = (TNtuple*)FileObj[dataSetName.c_str()]->Get(nTupleName.c_str());
//     //     nEntries = (int)nTuple[dataSetName.c_str()]->GetEntries();
//     //     cout<<"                 nEntries: "<<nEntries<<endl;

//     //     //sVarofinterest = "HT";
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sVarofinterest.c_str(),&varofInterest);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("ScaleFactor",&ScaleFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("NormFactor",&NormFactor);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress("Luminosity",&Luminosity);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar1.c_str(), &splitVar1);
//     //     nTuple[dataSetName.c_str()]->SetBranchAddress(sSplitVar2.c_str(), &splitVar2);

//     //     ////****************************************************Define plots***********************************************
//     //     for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     //     {
//     //         numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //         for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //             numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();            
//     //             histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //             //cout<<"histoName: "<< histoName<<endl;
//     //             histo1D[histoName.c_str()] = new TH1F(histoName.c_str(),histoName.c_str(), nBins, plotLow, plotHigh);
//     //         }
//     //     }
//     //     for (int i = 0; i<nEntries; i++)   //Fill histo with variable of interest
//     //     {
//     //         nTuple[dataSetName.c_str()]->GetEntry(i);
//     //         //artificial Lumi
//     //         if(lScale > 0 )
//     //         {
//     //             Luminosity = 1000*lScale;
//     //         }
//     //         bool outsiderange = true;
//     //         if(splitVar1 >= ftSplit1) //Check if this entry belongs in the last bin.  Done here to optimize number of checks
//     //         {
//     //             if(splitVar2 >= ftSplit2){ //Check if this entry belongs in the last bin in var2.
//     //                 histoName = dataSetName + stSplit1 + sSplitVar1 + stSplit2 + sSplitVar2;
//     //                 outsiderange = false;
//     //             } 
//     //             else //If it doesn't belong in the last bin, find out which it belongs in
//     //             {
//     //                 for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2)
//     //                 {
//     //                     if(splitVar2>=t2 && splitVar2<(t2+fwSplit2)) //splitVar falls inside one of the bins
//     //                     {
//     //                         numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //                         histoName = dataSetName + stSplit1+ sSplitVar1 + numStr2 + sSplitVar2;
//     //                         outsiderange = false;
//     //                         break;
//     //                     }
//     //                 }
//     //             } 
//     //         }
//     //         else //If it doesn't belong in the last bin for var1, find out which it belongs in
//     //         {
//     //             for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     //             {
//     //                 numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //                 if(splitVar1>=s && splitVar1<(s+fwSplit1)) //splitVar falls inside one of the bins
//     //                 {

//     //                     if(splitVar2 >= ftSplit2){ //Check if this entry belongs in the last bin in var2.
//     //                         histoName = dataSetName + numStr1 + sSplitVar1 + stSplit2 + sSplitVar2;
//     //                         outsiderange = false;
//     //                     }    
//     //                     else //If it doesn't belong in the last bin, find out which it belongs in
//     //                     {
//     //                         for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2)
//     //                         {
//     //                             if(splitVar2>=t2 && splitVar2<(t2+fwSplit2)) //splitVar falls inside one of the bins
//     //                             {
//     //                                 numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //                                 histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //                                 outsiderange =false;
//     //                                 break;
//     //                             }
//     //                         }
//     //                     } 
//     //                     //break;
//     //                 }
//     //             }
//     //         }
//     //         if (outsiderange == true){
//     //             continue;
//     //         }
//     //         histo1D[histoName.c_str()]->Fill(varofInterest,ScaleFactor*NormFactor*Luminosity);
//     //     }

//     //     ////****************************************************Fill plots***********************************************

//     //     TCanvas *canv2 = new TCanvas();

//     //     for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     //     {
//     //         numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //         for(int t2 = fbSplit2; t2 <= ftSplit2; t2+=fwSplit2){
//     //             shapefile->cd();
//     //             numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //             histoName = dataSetName + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //             if(Normalise)
//     //             {
//     //                 double dIntegral = histo1D[histoName.c_str()]->Integral();
//     //                 histo1D[histoName.c_str()]->Scale(1./dIntegral);
//     //             }
//     //             histo1D[histoName.c_str()]->Draw();
//     //             string writename = "";

//     //             if(dataSetName == "TTScaleDown")
//     //             {
//     //                 writename = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__TTJets__scaleDown";
//     //                 //cout<<"writename  :"<<writename<<endl;
//     //                 histo1D[histoName.c_str()]->Write((writename).c_str());
//     //                 canv2->SaveAs(("Sys_"+histoName+".png").c_str());
//     //                 errorfile->cd();
//     //                 //errorfile->mkdir(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //                 errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //                 histo1D[histoName.c_str()]->Write("Minus");
//     //                 //errorfile->Write();
//     //             }

//     //             if(dataSetName == "TTScaleUp")
//     //             {
//     //                 writename = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__TTJets__scaleUp";
//     //                 //cout<<"writename  :"<<writename<<endl;
//     //                 histo1D[histoName.c_str()]->Write((writename).c_str());
//     //                 canv2->SaveAs(("Sys_"+histoName+".png").c_str());
//     //                 errorfile->cd();
//     //                 errorfile->cd(("MultiSamplePlot_"+sVarofinterest).c_str());
//     //                 histo1D[histoName.c_str()]->Write("Plus");
//     //                 //errorfile->Write();
//     //             }
//     //         }
//     //     }
//     // }
// };

void DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, string xmlNom, float lScale){
    TTreeLoader treeLoader;
    vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
    const char *xmlfile = xmlNom.c_str();
    treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;
    int nDatasets = datasets.size();
    TH1F *tempHisto;
    float tempEntries;
    int nChannels = 1;
    cout<<""<<endl;

    cout<<"PRODUCING DATACARD"<<endl;

    string binname, histoName, dataSetName, datacardname;
    ofstream card;
    string slScale = intToStr(lScale);

    if(lScale>0){
        datacardname = "datacard" + leptoAbbr + "_" + slScale + "_inc.txt";   
    }
    else{
        datacardname = "datacard" + leptoAbbr + "_inc.txt";          
    }    
    card.open (datacardname.c_str());

    card << "imax " + static_cast<ostringstream*>( &(ostringstream() << nChannels) )->str() + "\n";
    card << "jmax " + static_cast<ostringstream*>( &(ostringstream() << nDatasets-2) )->str() + "\n";
    card << "kmax *\n";
    card << "---------------\n";
    card << "shapes * * "+shapefileName+"  $CHANNEL__$PROCESS__nominal  $CHANNEL__$PROCESS__$SYSTEMATIC\n";
    card << "---------------\n";
    card << "bin                               ";

    binname = channel;
    card << binname + " ";
  
    card << "\n";
    card << "observation                               ";

    tempEntries = 0;
    for (int j = 0; j<nDatasets; j++){
        dataSetName=datasets[j]->Name();
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
        {
            histoName = channel + "__data_obs__nominal";
            tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
            tempEntries = tempHisto->GetSumOfWeights();  
            card<<static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str()+"         ";
        }
        else{
            continue;
        }
    }

    card << "\n";

    card << "---------------------------\n";
    card << "bin                               ";

    for(int i = 0; i<nDatasets-1; i++)    card << binname + " ";

    card << "\n";
    card << "process                             ";

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
 

    card << "\n";
    card << "process                                ";

    for(int i=0; i<datasets.size()-1; i++){
        card << static_cast<ostringstream*>( &(ostringstream() << i) )->str() + "                    ";
    }
  
    card << "\n";
    card << "rate                                ";

            tempEntries = 0;
    for (int j = 0; j<nDatasets; j++){
        dataSetName=datasets[j]->Name();
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
        {
            continue;
        }
        else{
            histoName = channel + "__" + dataSetName +"__nominal";
            tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
            tempEntries = tempHisto->GetSumOfWeights();  
            card << static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str() + "               ";
        }
    }
   
    card << "\n";
    card << "---------------------------\n";
    card << "lumi                  lnN           ";
    for (int i=0; i<nChannels*(nDatasets-1); i++){
        card << "1.046                ";
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
                    card << "1.1                 -                    -                    -                  -                    ";
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

    card << "btag                shape           ";
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

    card << "PU                shape           ";
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

    card << "JES                shape           ";
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

    card << "JER                shape           ";
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

// void Split_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string xmlNom, float lScale){
//     // TTreeLoader treeLoader;
//     // vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
//     // const char *xmlfile = xmlNom.c_str();
//     // treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;
//     // int nDatasets = datasets.size();
//     // TH1F *tempHisto;
//     // float tempEntries;
//     // int nChannels = 0;

//     // cout<<""<<endl;
//     // cout<<"PRODUCING DATACARD"<<endl;

//     // string numStr1, binname, histoName, dataSetName;
//     // ofstream card;
//     // std::string slScale = intToStr(lScale);
//     // string datacardname = "datacard" + leptoAbbr + "_" + slScale + "_JS" + ".txt";
//     // card.open (datacardname.c_str());
//     // cout<<"datacard name "<<datacardname<<endl;
//     // cout<< fbSplit1 <<"  "<<ftSplit1<<"  "<<fwSplit1<<endl;
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     nChannels += 1;
//     // } 
//     // card << "imax " + static_cast<ostringstream*>( &(ostringstream() << nChannels) )->str() + "\n"; 
//     // card << "jmax " + static_cast<ostringstream*>( &(ostringstream() << nDatasets-2) )->str() + "\n";  //number of background (so minus data and signal)
//     // card << "kmax *\n";
//     // card << "---------------\n";
//     // card << "shapes * * "+shapefileName+"  $CHANNEL__$PROCESS__nominal  $CHANNEL__$PROCESS__$SYSTEMATIC\n";
//     // card << "---------------\n";
//     // card << "bin                               ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     binname = channel + numStr1 + sSplitVar1;
//     //     card << binname + " ";
//     // }    
//     // card << "\n";
//     // card << "observation                               ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     binname = channel + numStr1 + sSplitVar1;
//     //     tempEntries = 0;
//     //     for (int j = 0; j<nDatasets; j++){
//     //         dataSetName=datasets[j]->Name();
//     //         if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //         {
//     //             histoName = channel + numStr1 + sSplitVar1 + "__data_obs__nominal";
//     //             tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
//     //             tempEntries = tempHisto->GetSumOfWeights();  
//     //             card<<static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str()+"         ";                
//     //         }
//     //         else{
//     //             continue;
//     //         }
//     //     }
//     // } 
//     // card << "\n";

//     // card << "---------------------------\n";
//     // card << "bin                               ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     binname = channel + numStr1 + sSplitVar1;
//     //     for(int i = 0; i<nDatasets-1; i++)    card << binname + " ";
//     // }
//     // card << "\n";
//     // card << "process                             ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     binname = channel + numStr1 + sSplitVar1;
//     //     for (int j = 0; j<nDatasets; j++){
//     //         dataSetName=datasets[j]->Name();
//     //         if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //         {
//     //             continue;
//     //         }
//     //         else{
//     //             card<<dataSetName+"           ";
//     //         }
//     //     }
//     // }    

//     // card << "\n";
//     // card << "process                                ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     for(int i=0; i<datasets.size()-1; i++){
//     //         card << static_cast<ostringstream*>( &(ostringstream() << i) )->str() + "                    ";
//     //     }
//     // }  
//     // card << "\n";
//     // card << "rate                                ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();

//     //     binname = channel + numStr1 + sSplitVar1;
//     //     tempEntries = 0;
//     //     for (int j = 0; j<nDatasets; j++){
//     //         dataSetName=datasets[j]->Name();
//     //         if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //         {
//     //             continue;
//     //         }
//     //         else{
//     //             histoName = channel + numStr1 + sSplitVar1 + "__" + dataSetName +"__nominal";
//     //             tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
//     //             tempEntries = tempHisto->GetSumOfWeights();  
//     //             card << static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str() + "               ";
//     //         }
//     //     }
//     // }     
//     // card << "\n";
//     // card << "---------------------------\n";
//     // card << "lumi                  lnN           ";
//     // for (int i=0; i<nChannels*(nDatasets-1); i++){
//     //     card << "1.026                ";
//     // }
//     // card << "\n";
//     // for (int d = 0; d<nDatasets; d++){
//     //     dataSetName = datasets[d]->Name();
//     //     if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //     {
//     //         continue;
//     //     }
//     //     else {
//     //         if(dataSetName.find("tttt")!=string::npos){
//     //             card << "tttt_norm              lnN          ";
//     //             for (int j = 0; j<nChannels; j++){
//     //                 card << "1.1                 -                    -                    -                  -                    ";
//     //             }
//     //             card << "\n";
//     //         }
//     //         else{
//     //             card << dataSetName + "_norm      lnN           ";
//     //             for(int k = 0; k<nChannels; k++){
//     //                 for (int dash1 = 0; dash1<d-1; dash1++){
//     //                     card << "-                  ";
//     //                 }
//     //                 card << "1.04                  ";                   
//     //                 for (int dash2 = 5; dash2>d; dash2-- ){
//     //                     card << "-                  ";
//     //                 }
//     //             }
//     //             card<<"\n";
//     //         }   
//     //     }        
//     // }

//     // card << "scale                shape           ";
//     // for (int d = 0; d<nDatasets; d++){
//     //     dataSetName = datasets[d]->Name();
//     //     if(dataSetName.find("TTJets")!=string::npos)
//     //     {
//     //         for(int k = 0; k<nChannels; k++){
//     //             for (int dash1 = 0; dash1<d-1; dash1++){
//     //                 card << "-                  ";
//     //             }
//     //             card << "1                      ";                   
//     //             for (int dash2 = 5; dash2>d; dash2-- ){
//     //                 card << "-                  ";
//     //             }
//     //         }
//     //         card<<"\n";
//     //     }
//     //     else {
//     //         continue;
//     //     }        
//     // }    

//     // card.close();
// };

// void Split2_DataCardProducer(TFile *shapefile, string shapefileName, string channel, string leptoAbbr, bool jetSplit, bool jetTagsplit, string sSplitVar1, float fbSplit1, float ftSplit1, float fwSplit1, string sSplitVar2, float fbSplit2, float ftSplit2, float fwSplit2, string xmlNom, float lScale){
//     // TTreeLoader treeLoader;
//     // vector < Dataset* > datasets;                   //cout<<"vector filled"<<endl;
//     // const char *xmlfile = xmlNom.c_str();
//     // treeLoader.LoadDatasets (datasets, xmlfile);    //cout<<"datasets loaded"<<endl;
//     // int nDatasets = datasets.size();
//     // TH1F *tempHisto;
//     // float tempEntries;
//     // int nChannels = 0;
//     // cout<<""<<endl;
//     // cout<<"PRODUCING DATACARD"<<endl;

//     // string numStr1, numStr2, binname, histoName, dataSetName;
//     // ofstream card;
//     // std::string slScale = intToStr(lScale);
//     // string datacardname = "datacard" + leptoAbbr + "_" + slScale + "_JTS" + ".txt";
//     // card.open (datacardname.c_str());

//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //         //for(int i=0; i<datasets.size()-1; i++){
//     //             nChannels += 1;
//     //         //}
//     //     }  
//     // } 

//     // card << "imax " + static_cast<ostringstream*>( &(ostringstream() << nChannels) )->str() + "\n";
//     // card << "jmax " + static_cast<ostringstream*>( &(ostringstream() << nDatasets-2) )->str() + "\n";
//     // card << "kmax *\n";
//     // card << "---------------\n";
//     // card << "shapes * * "+shapefileName+"  $CHANNEL__$PROCESS__nominal  $CHANNEL__$PROCESS__$SYSTEMATIC\n";
//     // card << "---------------\n";
//     // card << "bin                               ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //         numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //         binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //         card << binname + " ";
//     //     }  
//     // }    
//     // card << "\n";
//     // card << "observation                               ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //         numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //         binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //         tempEntries = 0;
//     //         for (int j = 0; j<nDatasets; j++){
//     //             dataSetName=datasets[j]->Name();
//     //             if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //             {
//     //                 histoName = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__data_obs__nominal";
//     //                 tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
//     //                 tempEntries = tempHisto->GetSumOfWeights();  
//     //                 card<<static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str()+"         ";                }
//     //             else{
//     //                 continue;
//     //             }
//     //         }
//     //     }  
//     // } 
//     // card << "\n";

//     // card << "---------------------------\n";
//     // card << "bin                               ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //         numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //         binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //         for(int i = 0; i<nDatasets-1; i++)    card << binname + " ";
//     //     }  
//     // }
//     // card << "\n";
//     // card << "process                             ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //         binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //         for (int j = 0; j<nDatasets; j++){
//     //             dataSetName=datasets[j]->Name();
//     //             if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //             {
//     //                 continue;
//     //             }
//     //             else{
//     //                 card<<dataSetName+"           ";
//     //             }
//     //         }
//     //     }  
//     // }    

//     // card << "\n";
//     // card << "process                                ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //         for(int i=0; i<datasets.size()-1; i++){
//     //             card << static_cast<ostringstream*>( &(ostringstream() << i) )->str() + "                    ";
//     //         }
//     //     }  
//     // }  
//     // card << "\n";
//     // card << "rate                                ";
//     // for(int s = fbSplit1; s <= ftSplit1; s+=fwSplit1)
//     // {
//     //     numStr1 = static_cast<ostringstream*>( &(ostringstream() << s) )->str();
//     //     for(int t2 = fbSplit2; t2<= ftSplit2; t2+=fwSplit2){
//     //         numStr2 = static_cast<ostringstream*>( &(ostringstream() << t2) )->str();
//     //         binname = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2;
//     //         tempEntries = 0;
//     //         for (int j = 0; j<nDatasets; j++){
//     //             dataSetName=datasets[j]->Name();
//     //             if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //             {
//     //                 continue;
//     //             }
//     //             else{
//     //                 histoName = channel + numStr1 + sSplitVar1 + numStr2 + sSplitVar2 + "__" + dataSetName +"__nominal";
//     //                 tempHisto = (TH1F*)shapefile->Get(histoName.c_str());
//     //                 tempEntries = tempHisto->GetSumOfWeights();  
//     //                 card << static_cast<ostringstream*>( &(ostringstream() << tempEntries) )->str() + "               ";
//     //             }
//     //         }
//     //     }  
//     // }     
//     // card << "\n";
//     // card << "---------------------------\n";
//     // card << "lumi                  lnN           ";
//     // for (int i=0; i<nChannels*(nDatasets-1); i++){
//     //     card << "1.026                ";
//     // }
//     // card << "\n";
//     // for (int d = 0; d<nDatasets; d++){
//     //     dataSetName = datasets[d]->Name();
//     //     if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
//     //     {
//     //         continue;
//     //     }
//     //     else {
//     //         if(dataSetName.find("tttt")!=string::npos){
//     //             card << "tttt_norm              lnN          ";
//     //             for (int j = 0; j<nChannels; j++){
//     //                 card << "1.1                 -                    -                    -                  -                    ";
//     //             }
//     //             card << "\n";
//     //         }
//     //         else{
//     //             card << dataSetName + "_norm      lnN           ";
//     //             for(int k = 0; k<nChannels; k++){
//     //                 for (int dash1 = 0; dash1<d-1; dash1++){
//     //                     card << "-                  ";
//     //                 }
//     //                 card << "1.04                  ";                   
//     //                 for (int dash2 = 5; dash2>d; dash2-- ){
//     //                     card << "-                  ";
//     //                 }
//     //             }
//     //             card<<"\n";
//     //         }   
//     //     }        
//     // }

//     // card << "scale                shape           ";
//     // for (int d = 0; d<nDatasets; d++){
//     //     dataSetName = datasets[d]->Name();
//     //     if(dataSetName.find("TTJets")!=string::npos)
//     //     {
//     //         for(int k = 0; k<nChannels; k++){
//     //             for (int dash1 = 0; dash1<d-1; dash1++){
//     //                 card << "-                  ";
//     //             }
//     //             card << "1                      ";                   
//     //             for (int dash2 = 5; dash2>d; dash2-- ){
//     //                 card << "-                  ";
//     //             }
//     //         }
//     //         card<<"\n";
//     //     }
//     //     else {
//     //         continue;
//     //     }        
//     // }    

//     // card.close();
// };




std::string intToStr (int number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

