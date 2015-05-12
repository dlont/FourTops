//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.                  ////
////////////////////////////////////////////////////////////////////////////////////

// ttbar @ NLO 13 TeV:
//all-had ->679 * .46 = 312.34
//semi-lep ->679 *.45 = 305.55
//di-lep-> 679* .09 = 61.11

//ttbar @ NNLO 8 TeV:
//all-had -> 245.8 * .46 = 113.068
//semi-lep-> 245.8 * .45 = 110.61
//di-lep ->  245.8 * .09 = 22.122

#define _USE_MATH_DEFINES
#include "TStyle.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TNtuple.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <ctime>

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <errno.h>
#include "TRandom3.h"
#include "TRandom.h"
#include "TProfile.h"
#include <iostream>
#include <map>
#include <cstdlib>

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Selection/interface/Run2Selection.h"
//#include "TopTreeAnalysisBase/Selection/interface/FourTopSelectionTable.h"

#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"

#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h
//#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"


#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

using namespace std;
using namespace TopTree;
using namespace reweight;

bool split_ttbar = false;
bool debug = false;
float topness;

pair<float, vector<unsigned int> > MVAvals1;
pair<float, vector<unsigned int> > MVAvals2;
pair<float, vector<unsigned int> > MVAvals2ndPass;
pair<float, vector<unsigned int> > MVAvals3rdPass;

int nMVASuccesses=0;
int nMatchedEvents=0;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TProfile*> histoProfile;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;

struct HighestTCHEBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_trackCountingHighEffBJetTags() > j2->btag_trackCountingHighEffBJetTags();
    }
};
struct HighestCVSBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_combinedInclusiveSecondaryVertexV2BJetTags() > j2->btag_combinedInclusiveSecondaryVertexV2BJetTags();
    }
};

bool match;

//To cout the Px, Py, Pz, E and Pt of objects
int Factorial(int N);
float Sphericity(vector<TLorentzVector> parts );
float Centrality(vector<TLorentzVector> parts);

int main (int argc, char *argv[])
{

    //Checking Passed Arguments to ensure proper execution of MACRO
    if(argc < 14)
    {
        std::cerr << "TOO FEW INPUTs FROM XMLFILE.  CHECK XML INPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
        return 1;
    }

    //Placing arguments in properly typed variables for Dataset creation

    const string dName              = argv[1];
    const string dTitle             = argv[2];
    const int color                 = strtol(argv[4], NULL, 10);
    const int ls                    = strtol(argv[5], NULL, 10);
    const int lw                    = strtol(argv[6], NULL, 10);
    const float normf               = strtod(argv[7], NULL);
    const float EqLumi              = strtod(argv[8], NULL);
    const float xSect               = strtod(argv[9], NULL);
    const float PreselEff           = strtod(argv[10], NULL);
    string fileName                 = argv[11];
    const int startEvent            = strtol(argv[argc-2], NULL, 10);
    const int endEvent              = strtol(argv[argc-1], NULL, 10);
    vector<string> vecfileNames;
    for(int args = 11; args < argc-2; args++)
    {
        vecfileNames.push_back(argv[args]);
    }




    cout << "---Dataset accepted from command line---" << endl;
    cout << "Dataset Name: " << dName << endl;
    cout << "Dataset Title: " << dTitle << endl;
    cout << "Dataset color: " << color << endl;
    cout << "Dataset ls: " << ls << endl;
    cout << "Dataset lw: " << lw << endl;
    cout << "Dataset normf: " << normf << endl;
    cout << "Dataset EqLumi: " << EqLumi << endl;
    cout << "Dataset xSect: " << xSect << endl;
    cout << "Dataset File Name: " << vecfileNames[0] << endl;
    cout << "Beginning Event: " << startEvent << endl;
    cout << "Ending Event: " << endEvent << endl;
    cout << "----------------------------------------" << endl;
//    cin.get();



    ofstream eventlist;
    eventlist.open ("interesting_events_mu.txt");

    int passed = 0;
    int ndefs =0;
    int negWeights = 0;
    float weightCount = 0.0;
    int eventCount = 0;

    string btagger = "CSVL";
    float scalefactorbtageff, mistagfactor;
    float workingpointvalue = 0.679; //working points updated to 2012 BTV-POG recommendations.
    bool bx25 = false;

    if(btagger == "CSVL")
        workingpointvalue = .244;
    else if(btagger == "CSVM")
        workingpointvalue = .679;
    else if(btagger == "CSVT")
        workingpointvalue = .898;

    clock_t start = clock();

    BTagWeightTools * bTool = new BTagWeightTools("SFb-pt_NOttbar_payload_EPS13.txt", "CSVM") ;

    int doJESShift = 0; // 0: off 1: minus 2: plus
    cout << "doJESShift: " << doJESShift << endl;

    int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
    cout << "doJERShift: " << doJERShift << endl;

    int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
    cout << "dobTagEffShift: " << dobTagEffShift << endl;

    int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
    cout << "domisTagEffShift: " << domisTagEffShift << endl;

    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the FourTop search ! "           << endl;
    cout << "*************************************************************" << endl;


    string postfix = "_Run2_TopTree_Study_" + dName; // to relabel the names of the output file

    if (doJESShift == 1)
        postfix= postfix+"_JESMinus";
    if (doJESShift == 2)
        postfix= postfix+"_JESPlus";
    if (doJERShift == 1)
        postfix= postfix+"_JERMinus";
    if (doJERShift == 2)
        postfix= postfix+"_JERPlus";
    if (dobTagEffShift == -1)
        postfix= postfix+"_bTagMinus";
    if (dobTagEffShift == 1)
        postfix= postfix+"_bTagPlus";
    if(domisTagEffShift == -1)
        postfix= postfix+"_misTagMinus";
    if(domisTagEffShift == 1)
        postfix= postfix+"_misTagPlus";

    ///////////////////////////////////////
    // Configuration
    ///////////////////////////////////////

    string channelpostfix = "";
    string xmlFileName = "";

    //Setting Lepton Channels (Setting both flags true will select Muon-Electron Channel when dilepton is also true)
    bool dilepton = true;
    bool Muon = true;
    bool Electron = false;

    if(Muon && Electron && dilepton)
    {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_MuEl";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else if(Muon && !Electron && dilepton)
    {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_MuMu";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else if(!Muon && Electron && dilepton)
    {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_ElEl";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else
    {
        cerr<<"Correct Di-lepton Channel not selected."<<endl;
        exit(1);
    }

    bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
    bool trainEventMVA = false; // If false, the previously trained MVA will be used to calculate stuff
    bool computeEventMVA = false;


    const char *xmlfile = xmlFileName.c_str();
    cout << "used config file: " << xmlfile << endl;

    /////////////////////////////
    //  Set up AnalysisEnvironment
    /////////////////////////////

    AnalysisEnvironment anaEnv;
    cout<<" - Creating environment ..."<<endl;
//    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    anaEnv.PrimaryVertexCollection = "PrimaryVertex";
    anaEnv.JetCollection = "PFJets_slimmedJets";
    anaEnv.FatJetCollection = "FatJets_slimmedJetsAK8";
    anaEnv.METCollection = "PFMET_slimmedMETs";
    anaEnv.MuonCollection = "Muons_slimmedMuons";
    anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
    anaEnv.GenJetCollection   = "GenJets_slimmedGenJets";
    anaEnv.TrackMETCollection = "";
    anaEnv.GenEventCollection = "GenEvent";
    anaEnv.NPGenEventCollection = "NPGenEvent";
    anaEnv.MCParticlesCollection = "MCParticles";
    anaEnv.loadFatJetCollection = true;
    anaEnv.loadGenJetCollection = false;
    anaEnv.loadGenEventCollection = false;
    anaEnv.loadNPGenEventCollection = false;
    anaEnv.loadMCParticles = true;
    anaEnv.loadTrackMETCollection = false;
    anaEnv.JetType = 2;
    anaEnv.METType = 2;
    int verbose = 2;//anaEnv.Verbose;



    ////////////////////////////////
    //  Load datasets
    ////////////////////////////////

    TTreeLoader treeLoader;
    vector < Dataset* > datasets;
    cout << " - Creating Dataset ..." << endl;
    Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
    theDataset->SetEquivalentLuminosity(EqLumi);
    datasets.push_back(theDataset);
    float Luminosity = 5000.0; //pb^-1??
    vector<string> MVAvars;

    MVAvars.push_back("topness");
    MVAvars.push_back("muonpt");
    MVAvars.push_back("muoneta");
    MVAvars.push_back("HTH");
    MVAvars.push_back("HTRat");
    MVAvars.push_back("HTb");
    MVAvars.push_back("nLtags");
    MVAvars.push_back("nMtags");
    MVAvars.push_back("nTtags");
    MVAvars.push_back("nJets");
    MVAvars.push_back("Jet3Pt");
    MVAvars.push_back("Jet4Pt");

    MVAComputer* Eventcomputer_ = new MVAComputer("BDT","MasterMVA_Mu_10thMarch.root","MasterMVA_Mu_10March",MVAvars, "_10thMarch2015");

    cout << " Initialized Eventcomputer_" << endl;



    string dataSetName;

    string MVAmethod = "BDT"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)

    cout <<"Instantiating jet combiner..."<<endl;

    JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, false);
    cout <<"Instantiated jet combiner..."<<endl;


    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////

    cout <<"found sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    dataSetName = theDataset->Name();
    if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
    {
        Luminosity = theDataset->EquivalentLumi();
        cout <<"found DATA sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    }

    cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;
    int ndatasets = datasets.size() - 1 ;

    double currentLumi;
    double newlumi;

    //Output ROOT file
    string outputDirectory("MACRO_Output"+channelpostfix);
    mkdir(outputDirectory.c_str(),0777);
    string rootFileName (outputDirectory+"/FourTop"+postfix+channelpostfix+".root");
    TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

    //vector of objects
    cout << " - Variable declaration ..." << endl;
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
    vector < TRootJet* >      init_fatjets;
    vector < TRootMET* >      mets;

    //Global variable
    TRootEvent* event = 0;

    ////////////////////////////////////////////////////////////////////
    ////////////////// MultiSample plots  //////////////////////////////
    ////////////////////////////////////////////////////////////////////

    MSPlot["NbOfVertices"]                                  = new MultiSamplePlot(datasets, "NbOfVertices", 60, 0, 60, "Nb. of vertices");
    //Muons
    MSPlot["MuonPt"]                                        = new MultiSamplePlot(datasets, "MuonPt", 30, 0, 300, "PT_{#mu}");
    //Electrons
    MSPlot["ElectronRelIsolation"]                          = new MultiSamplePlot(datasets, "ElectronRelIsolation", 10, 0, .25, "RelIso");
    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]                             = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV", 20, 0, 1, "CSV b-disc.");
    //Jets
    MSPlot["JetEta"]                                        = new MultiSamplePlot(datasets, "JetEta", 40,-4, 4, "Jet #eta");
    MSPlot["HT_SelectedJets"]                               = new MultiSamplePlot(datasets, "HT_SelectedJets", 30, 0, 1500, "HT");
    MSPlot["HTExcess2M"]                                    = new MultiSamplePlot(datasets, "HTExcess2M", 30, 0, 1500, "HT_{Excess 2 M b-tags}");
    MSPlot["HTH"]                                           = new MultiSamplePlot(datasets, "HTH", 20, 0, 1, "HTH");
    //MET
    MSPlot["MET"]                                           = new MultiSamplePlot(datasets, "MET", 70, 0, 700, "MET");
    MSPlot["METCutEff"]                                     = new MultiSamplePlot(datasets, "METCutEff", 30, 0, 150, "MET cut pre-jet selection");
    //Topology Plots
    MSPlot["TotalSphericity"]                               = new MultiSamplePlot(datasets, "TotalSphericity", 20, 0, 1, "Sphericity_{all}");
    MSPlot["TotalCentrality"]                               = new MultiSamplePlot(datasets, "TotalCentrality", 20, 0, 1, "Centrality_{all}");
    MSPlot["DiLepSphericity"]                               = new MultiSamplePlot(datasets, "DiLepSphericity", 20, 0, 1, "Sphericity_{ll}");
    MSPlot["DiLepCentrality"]                               = new MultiSamplePlot(datasets, "DiLepCentrality", 20, 0, 1, "Centrality_{ll}");
    MSPlot["TopDiLepSphericity"]                            = new MultiSamplePlot(datasets, "TopDiLepSphericity", 20, 0, 1, "Sphericity_{tll}");
    MSPlot["TopDiLepCentrality"]                            = new MultiSamplePlot(datasets, "TopDiLepCentrality", 20, 0, 1, "Centrality_{tll}");
    //MVA Top Roconstruction Plots
    MSPlot["MVA1TriJet"]                                    = new MultiSamplePlot(datasets, "MVA1TriJet", 30, -1.0, 0.2, "MVA1TriJet");
    MSPlot["MVA1TriJetMass"]                                = new MultiSamplePlot(datasets, "MVA1TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1DiJetMass"]                                 = new MultiSamplePlot(datasets, "MVA1DiJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1PtRat"]                                     = new MultiSamplePlot(datasets, "MVA1PtRat", 25, 0, 2, "P_{t}^{Rat}");
    MSPlot["MVA1BTag"]                                      = new MultiSamplePlot(datasets, "MVA1BTag", 35, 0, 1, "BTag");
    MSPlot["MVA1AnThBh"]                                    = new MultiSamplePlot(datasets, "MVA1AnThBh", 35, 0, 3.14, "AnThBh");
    MSPlot["MVA1AnThWh"]                                    = new MultiSamplePlot(datasets, "MVA1AnThWh", 35, 0, 3.14, "AnThWh");
    MSPlot["MVA1dPhiThDiLep"]                               = new MultiSamplePlot(datasets, "MVA1dPhiThDiLep", 35, 0, 3.14, "dPhiThDiLep");
    MSPlot["MVA1dRThDiLep"]                                 = new MultiSamplePlot(datasets, "MVA1dRThDiLep", 35, 0, 3.14, "dRThDiLep");

    //ZMass window plots
    MSPlot["ZMassWindowWidthEff"]                           = new MultiSamplePlot(datasets, "ZMassWindowWidthEff", 20, 0, 100, "Z mass window width");
    MSPlot["DiLepMass"]                                     = new MultiSamplePlot(datasets, "DiLepMass", 30, 0, 150, "m_{ll}");


    ///////////////////
    // 1D histograms //
    ///////////////////

    ///////////////////
    // 2D histograms //
    ///////////////////
    histo2D["HTLepSep"] = new TH2F("HTLepSep","dR_{ll}:HT",50,0,1000, 20, 0,4);

    //Plots
    string pathPNG = "MSPlots_FourTop"+postfix+channelpostfix;
    pathPNG += "_MSPlots/";
    //pathPNG = pathPNG +"/";
    mkdir(pathPNG.c_str(),0777);

    cout <<"Making directory :"<< pathPNG  <<endl;
    vector<string> CutsselecTable;
    if(dilepton)
    {
        /////////////////////////////////
        // Selection table: Dilepton + jets
        /////////////////////////////////
        if(Muon && Electron)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("Exactly 1 Loose Isolated Muon"));
            CutsselecTable.push_back(string("Exactly 1 Loose Electron"));
            CutsselecTable.push_back(string("At least 4 Jets"));
            CutsselecTable.push_back(string("At least 1 CSVM Jet"));
            CutsselecTable.push_back(string("At least 2 CSVM Jets"));
            CutsselecTable.push_back(string("Exactly 5 Jets"));
            CutsselecTable.push_back(string("Exactly 6 Jets"));
            CutsselecTable.push_back(string("Exactly 7 jets"));
            CutsselecTable.push_back(string("At Least 8 Jets"));
        }
        if(Muon && !Electron)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("Exactly 2 Loose Isolated Muon"));
            CutsselecTable.push_back(string("Z Mass Veto"));
            CutsselecTable.push_back(string("At least 4 Jets"));
            CutsselecTable.push_back(string("At least 1 CSVM Jet"));
            CutsselecTable.push_back(string("At least 2 CSVM Jets"));
            CutsselecTable.push_back(string("Exactly 5 Jets"));
            CutsselecTable.push_back(string("Exactly 6 Jets"));
            CutsselecTable.push_back(string("Exactly 7 jets"));
            CutsselecTable.push_back(string("At Least 8 Jets"));
        }
        if(!Muon && Electron)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("Exactly 2 Loose Electron"));
            CutsselecTable.push_back(string("Z Mass Veto"));
            CutsselecTable.push_back(string("At least 4 Jets"));
            CutsselecTable.push_back(string("At least 1 CSVM Jet"));
            CutsselecTable.push_back(string("At least 2 CSVM Jets"));
            CutsselecTable.push_back(string("Exactly 5 Jets"));
            CutsselecTable.push_back(string("Exactly 6 Jets"));
            CutsselecTable.push_back(string("Exactly 7 jets"));
            CutsselecTable.push_back(string("At Least 8 Jets"));
        }
    }


    SelectionTable selecTable(CutsselecTable, datasets);
    selecTable.SetLuminosity(Luminosity);
    selecTable.SetPrecision(1);

    /////////////////////////////////
    // Loop on datasets
    /////////////////////////////////

    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        string previousFilename = "";
        int iFile = -1;
        bool nlo = false;
        dataSetName = datasets[d]->Name();
        if(dataSetName.find("bx50") != std::string::npos) bx25 = false;
        else bx25 = true;

        if(dataSetName.find("NLO") != std::string::npos || dataSetName.find("nlo") !=std::string::npos) nlo = true;
        else nlo = false;

        if(bx25) cout << "Dataset with 25ns Bunch Spacing!" <<endl;
        else cout << "Dataset with 50ns Bunch Spacing!" <<endl;
        if(nlo) cout << "NLO Dataset!" <<endl;
        else cout << "LO Dataset!" << endl;


        //////////////////////////////////////////////
        // Setup Date string and nTuple for output  //
        //////////////////////////////////////////////

        time_t t = time(0);   // get time now
        struct tm * now = localtime( & t );

        int year = now->tm_year + 1900;
        int month =  now->tm_mon + 1;
        int day = now->tm_mday;
        int hour = now->tm_hour;
        int min = now->tm_min;
        int sec = now->tm_sec;

        string year_str;
        string month_str;
        string day_str;
        string hour_str;
        string min_str;
        string sec_str;

        ostringstream convert;   // stream used for the conversion
        convert << year;      // insert the textual representation of 'Number' in the characters in the stream
        year_str = convert.str();
        convert.str("");
        convert.clear();
        convert << month;      // insert the textual representation of 'Number' in the characters in the stream
        month_str = convert.str();
        convert.str("");
        convert.clear();
        convert << day;      // insert the textual representation of 'Number' in the characters in the stream
        day_str = convert.str();
        convert.str("");
        convert.clear();
        convert << hour;      // insert the textual representation of 'Number' in the characters in the stream
        hour_str = convert.str();
        convert.str("");
        convert.clear();
        convert << min;      // insert the textual representation of 'Number' in the characters in the stream
        min_str = convert.str();
        convert.str("");
        convert.clear();
        convert << day;      // insert the textual representation of 'Number' in the characters in the stream
        sec_str = convert.str();
        convert.str("");
        convert.clear();


        string date_str = day_str + "_" + month_str + "_" + year_str;

        cout <<"DATE STRING   "<<date_str << endl;

        //string dataSetName = datasets[d]->Name();
        string channel_dir = "Craneens"+channelpostfix;
        string date_dir = channel_dir+"/Craneens" + date_str +"/";
        int mkdirstatus = mkdir(channel_dir.c_str(),0777);
        mkdirstatus = mkdir(date_dir.c_str(),0777);



        //     string Ntupname = "Craneens/Craneen_" + dataSetName +postfix + "_" + date_str+  ".root";

        string Ntupname = "Craneens"+channelpostfix+"/Craneens"+ date_str  +"/Craneen_" + dataSetName +postfix + ".root";
        string Ntuptitle = "Craneen_" + channelpostfix;

        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");

        // TNtuple * tup = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(),"nJets:nLtags:nMtags:nTtags:HT:LeadingMuonPt:LeadingMuonEta:LeadingElectronPt:LeadingBJetPt:HT2M:HTb:HTH:HTRat:topness:ScaleFactor:PU:NormFactor:Luminosity:GenWeight");

        TNtuple * tup = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(),"BDT:nJets:nFatJets:nWTags:nTopTags:nLtags:nMtags:nTtags:HT:LeadingMuonPt:LeadingMuonEta:LeadingElectronPt:LeadingBJetPt:HT2L:HTb:HTH:HTRat:topness:ScaleFactor:PU:NormFactor:Luminosity:GenWeight");


        //////////////////////////////////////////////////
        /// Initialize JEC factors ///////////////////////
        //////////////////////////////////////////////////

        vector<JetCorrectorParameters> vCorrParam;

        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) // Data!
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
        }
        else
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L1FastJet_AK5PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L2Relative_AK5PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L3Absolute_AK5PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
        }
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_Uncertainty_AK5PFchs.txt");
//    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "SubTotalMC")));
//    JetCorrectionUncertainty *jecUncTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "Total")));

        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);

        //////////////////////////////////////////////////
        // Loop on events
        /////////////////////////////////////////////////

        int itrigger = -1, previousRun = -1;

        int start = 0;
        cout << "teh bugz!" << endl;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();

        cout <<"Number of events in total dataset = "<<  ending  <<endl;

        int event_start = startEvent;
        if (verbose > 1) cout << " - Loop over events " << endl;

        float BDTScore, MHT, MHTSig, STJet,muoneta, muonpt,electronpt,bjetpt, EventMass, EventMassX , SumJetMass, SumJetMassX,H,HX ,HTHi,HTRat, HT, HTX,HTH,HTXHX, sumpx_X, sumpy_X, sumpz_X, sume_X, sumpx, sumpy, sumpz, sume, jetpt,PTBalTopEventX,PTBalTopSumJetX , PTBalTopMuMet;

        double currentfrac =0.;
        double end_d;
        if(endEvent > ending)
            end_d = ending;
        else
            end_d = endEvent;

        //end_d = 10000; //artifical ending for debug
        int nEvents = end_d - event_start;
        cout <<"Will run over "<<  (end_d - event_start) << " events..."<<endl;
        cout <<"Starting event = = = = "<< event_start  << endl;

        //define object containers
        vector<TRootElectron*> selectedElectrons;
        vector<TRootPFJet*>    selectedJets;
        vector<TRootSubstructureJet*>    selectedFatJets;
        vector<TRootPFJet*>    MVASelJets1;
        vector<TRootMuon*>     selectedMuons;
        vector<TRootElectron*> selectedExtraElectrons;
        vector<TRootMuon*>     selectedExtraMuons;
        selectedElectrons.reserve(10);
        selectedMuons.reserve(10);

        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////

        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
            BDTScore= -99999.0, MHT = 0.,MHTSig = 0.,muoneta = 0., muonpt =0., electronpt=0., bjetpt =0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0., HTHi =0., HTRat = 0;
            H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0.;

            double ievt_d = ievt;
            currentfrac = ievt_d/end_d;
            if (debug)cout <<"event loop 1"<<endl;

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, init_fatjets,  mets, debug);  //load event
            if (debug)cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;
            float weight_0 = event->weight0();
            if (debug)cout <<"Weight0: " << weight_0 <<endl;
            if(nlo)
            {
                if(weight_0 < 0.0)
                {
                    scaleFactor = -1.0;  //Taking into account negative weights in NLO Monte Carlo
                    negWeights++;
                }
            }

            float rho = event->fixedGridRhoFastjetAll();
            if (debug)cout <<"Rho: " << rho <<endl;
            string graphName;

            //////////////////
            //Loading Gen jets
            //////////////////

            vector<TRootGenJet*> genjets;
            if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
            {
                // loading GenJets as I need them for JER
                genjets = treeLoader.LoadGenJet(ievt);
            }

            ///////////////////////
            // JER smearing
            //////////////////////

//            if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
//            {
//                //JER
//                doJERShift == 0;
//                if(doJERShift == 1)
//                    jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
//                else if(doJERShift == 2)
//                    jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
//                else
//                    jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
//
//                //     coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");
//
//                // JES sysematic!
//                if (doJESShift == 1)
//                    jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
//                else if (doJESShift == 2)
//                    jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
//
//                //            coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");
//
//            }

            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////

            // Apply trigger selection
//            trigged = treeLoader.EventTrigged (itrigger);
            bool trigged = true;  // Disabling the HLT requirement
            if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;
            if(!trigged)		   continue;  //If an HLT condition is not present, skip this event in the loop.
            // Declare selection instance
            Run2Selection selection(init_jets, init_fatjets, init_muons, init_electrons, mets);

            // Define object selection cuts
            if (Muon && Electron && dilepton)
            {
                if (debug)cout<<"Getting Jets"<<endl;
                selectedJets                                        = selection.GetSelectedJets(); // Relying solely on cuts defined in setPFJetCuts()
                selectedFatJets                                        = selection.GetSelectedFatJets(); // Relying solely on cuts defined in setPFJetCuts()

                if (debug)cout<<"Getting Tight Muons"<<endl;
                selectedMuons                                       = selection.GetSelectedMuons();
                if (debug)cout<<"Getting Loose Electrons"<<endl;
                selectedElectrons                                   = selection.GetSelectedElectrons("Loose","PHYS14",true); // VBTF ID
            }
            if (Muon && !Electron && dilepton)
            {
                if (debug)cout<<"Getting Jets"<<endl;
                selectedJets                                        = selection.GetSelectedJets(); // Relying solely on cuts defined in setPFJetCuts()
                if (debug)cout<<"Getting loose Muons"<<endl;
                selectedMuons                                       = selection.GetSelectedDiMuons();
                if (debug)cout<<"Getting Loose Electrons"<<endl;
                selectedElectrons                                   = selection.GetSelectedElectrons("Loose","PHYS14",true); // VBTF ID
            }
            if (!Muon && Electron && dilepton)
            {
                if (debug)cout<<"Getting Jets"<<endl;
                selectedJets                                        = selection.GetSelectedJets(); // Relying solely on cuts defined in setPFJetCuts()
                if (debug)cout<<"Getting Tight Muons"<<endl;
                selectedMuons                                       = selection.GetSelectedMuons();
                if (debug)cout<<"Getting Loose Electrons"<<endl;
                selectedElectrons                                   = selection.GetSelectedElectrons("Loose","PHYS14",true); // VBTF ID
            }


            vector<TRootJet*>      selectedLBJets;
            vector<TRootJet*>      selectedMBJets;
            vector<TRootJet*>      selectedTBJets;
            vector<TRootJet*>      selectedLightJets;

            int JetCut =0;
            int nMu, nEl, nLooseIsoMu;

            nMu = selectedMuons.size(); //Number of Muons in Event
            nEl = selectedElectrons.size(); //Number of Electrons in Event


            bool isTagged =false;
            vector<TLorentzVector> selectedMuonsTLV_JC, selectedElectronsTLV_JC, selectedLooseIsoMuonsTLV;
            vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV, mcMuonsTLV, mcPartonsTLV;
            vector<TRootMCParticle*> mcParticlesMatching_;
            vector<int> mcMuonIndex, mcPartonIndex;
            JetPartonMatching muonMatching, jetMatching;

            //////////////////////////////////
            // Preselection Lepton Operations //
            //////////////////////////////////

            float diLepMass = 0, diMuMass = 0;
            bool ZVeto = false;
            float ZMass = 91, ZMassWindow = 15;
            TLorentzVector lep1, lep2, diLep;

            for(int selmu = 0; selmu < selectedMuons.size(); selmu++)
            {
                selectedMuonsTLV_JC.push_back(*selectedMuons[selmu]);
            }

            for(int selel = 0; selel < selectedElectrons.size(); selel++)
            {
                selectedElectronsTLV_JC.push_back(*selectedElectrons[selel]);
            }

            if(nMu ==2 && nEl == 0 && Muon && !Electron)
            {
                lep1 = selectedMuonsTLV_JC[0];
                lep2 = selectedMuonsTLV_JC[1];
            }

            if(nEl == 2 && nMu == 0 && Electron && !Muon)
            {
                lep1 = selectedElectronsTLV_JC[0];
                lep2 = selectedElectronsTLV_JC[1];
            }
            if(nEl == 1 && nMu == 1 && Electron && Muon)
            {
                lep1 = selectedMuonsTLV_JC[0];
                lep2 = selectedElectronsTLV_JC[0];
            }
            diLep = lep1 + lep2;
            diLepMass = diLep.M();
            MSPlot["DiLepMass"]->Fill(diLepMass, datasets[d], true, Luminosity*scaleFactor );
            for(int mass = 0; mass < 20; mass++)
            {
                float windowRes = 2.5;
                if((diLepMass < (ZMass-(mass*windowRes))) || (diLepMass > (ZMass+(mass*windowRes)))) MSPlot["ZMassWindowWidthEff"]->Fill((2*mass*windowRes), datasets[d], true, Luminosity*scaleFactor );
            }

            ///////////////////////////////////////////////////////////////////////////////////
            // Preselection looping over Jet Collection                                      //
            // Summing HT and calculating leading, lagging, and ratio for Selected and BJets //
            ///////////////////////////////////////////////////////////////////////////////////
            float temp_HT = 0., HTb = 0.;

            double p_tags_tagged_mc = 1.;
            double p_tags_untagged_mc = 1.;
            double p_tags_tagged_data = 1.;
            double p_tags_untagged_data = 1.;
            double p_mc = 1., p_data = 1.;
            int jet_flavor;
            float eff=1   ;
            float scaled_eff=1 ;
            float a_eff = 1;
            float sf_a_eff = 1;
            double LightJeteff;
            double JetPt, JetEta;
            double SF_tag =1.;
            double event_weight = 1.;
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
                selectedJetsTLV.push_back(*selectedJets[seljet]);
                jet_flavor = selectedJets[seljet]->partonFlavour();
                JetPt = selectedJets[seljet]->Pt() ;
                JetEta = selectedJets[seljet]->Eta() ;
                if (JetPt > 800.) JetPt = 800;
                if (JetEta > 2.4)
                {
                    JetEta = 2.4;
                }
                else if (JetEta < -2.4)
                {
                    JetEta = -2.4;
                }


                if(fabs(jet_flavor) == 5 || fabs(jet_flavor) == 4  )
                {
                    SF_tag =  bTool->getSF(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,dobTagEffShift );
                    //  cout <<" "<<endl;
                    ////cout <<"jet SF nom "<< bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,0 )    <<endl;
                    //cout <<"jet SF minus "<< bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,-1 )    <<endl;
                    //cout <<"jet SF plus "<< bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,1 )    <<endl;
                }
                else
                {
                    //  cout <<" light jet... "<<endl;
                    SF_tag =  bTool->getSF(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,domisTagEffShift);
                }
                if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.244   )
                {

                    selectedLBJets.push_back(selectedJets[seljet]);
                    if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.679)
                    {
                        selectedMBJets.push_back(selectedJets[seljet]);

                        if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.898)
                        {
                            selectedTBJets.push_back(selectedJets[seljet]);
                        }

                    }
                }
                else
                {
                    selectedLightJets.push_back(selectedJets[seljet]);
                }
            }

            float nJets = selectedJets.size(); //Number of Jets in Event
            float nMtags = selectedMBJets.size(); //Number of CSVM tags in Event
            float nLtags = selectedLBJets.size(); //Number of CSVL tags in Event (includes jets that pass CSVM)
            float nTtags = selectedTBJets.size(); //Number of CSVL tags in Event (includes jets that pass CSVM)
            float nFatJets = selectedFatJets.size();

            //            cout <<" med tags ...   "<< nMtags   <<endl;

            float nTopTags = 0;
            float nWTags = 0;
            ////
            //// W/Top tagging
            ////

            for (int fatjet = 0; fatjet < nFatJets; fatjet++)
            {

                float tau1 = selectedFatJets[fatjet]->Tau1();
                float tau2 = selectedFatJets[fatjet]->Tau2();
                float prunedmass = selectedFatJets[fatjet]->PrunedMass();
                float nsubjets =  selectedFatJets[fatjet]->CmsTopTagNsubjets();
                float minmass =  selectedFatJets[fatjet]->CmsTopTagMinMass();
                float topmass =  selectedFatJets[fatjet]->CmsTopTagMass();

                // cout <<"llop in fat jet "<< " tau1 "   <<  tau1  << " tau2 "<< tau2  << " prunedmass " << prunedmass  << " nsubjets " << nsubjets  << " minmass "  <<minmass<< " topmass "  << topmass <<endl;



                //W-tagging
                if(  (tau2/tau1)  > 0.6 && prunedmass > 50.0 )
                {

                    nWTags++;
                    //cout <<"W-TAG!"<<endl;
                }

                //Top-tagging
                if(  nsubjets  > 2 && minmass > 50.0 &&  topmass > 150.0 )
                {
                    //cout <<"TOP-TAG!"<<endl;
                    nTopTags++;

                }


            }


            //////////////////////
            // Sync'ing cutflow //
            //////////////////////

            if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;
            // Apply primary vertex selection
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV << " TriggerBit: " << trigged <<endl;
            if (debug) cin.get();
            selecTable.Fill(d,0,scaleFactor);
            weightCount += scaleFactor;
            eventCount++;
            if(Muon && Electron && dilepton)   //Muon-Electron Selection Table
            {
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,scaleFactor);
                    if (nMu==1)
                    {
                        selecTable.Fill(d,2,scaleFactor);
                        if(nEl==1)
                        {
                            selecTable.Fill(d,3,scaleFactor);
                            if(nJets>=4)
                            {
                                selecTable.Fill(d,4,scaleFactor);
                                if(nMtags>=1)
                                {
                                    selecTable.Fill(d,5,scaleFactor);
                                    if(nMtags>=2)
                                    {
                                        selecTable.Fill(d,6,scaleFactor);
                                        if(nJets==5)
                                        {
                                            selecTable.Fill(d,7,scaleFactor);
                                            if(nJets==6)
                                            {
                                                selecTable.Fill(d,8,scaleFactor);
                                                if(nJets==7)
                                                {
                                                    selecTable.Fill(d,9,scaleFactor);
                                                    if(nJets==8)
                                                    {
                                                        selecTable.Fill(d,10,scaleFactor);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(Muon && !Electron && dilepton)   //Muon-Electron Selection Table
            {
                if(diLepMass < 20 || (diLepMass > (ZMass-ZMassWindow) && diLepMass < (ZMass+ZMassWindow))) ZVeto = true;
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,scaleFactor);
                    if (nMu == 2 && nEl == 0)
                    {
                        selecTable.Fill(d,2,scaleFactor);
                        if(!ZVeto)
                        {
                            selecTable.Fill(d,3,scaleFactor);
                            if(nJets>=4)
                            {
                                selecTable.Fill(d,4,scaleFactor);
                                if(nMtags>=1)
                                {
                                    selecTable.Fill(d,5,scaleFactor);
                                    if(nMtags>=2)
                                    {
                                        selecTable.Fill(d,6,scaleFactor);
                                        if(nJets==5)
                                        {
                                            selecTable.Fill(d,7,scaleFactor);
                                            if(nJets==6)
                                            {
                                                selecTable.Fill(d,8,scaleFactor);
                                                if(nJets==7)
                                                {
                                                    selecTable.Fill(d,9,scaleFactor);
                                                    if(nJets==8)
                                                    {
                                                        selecTable.Fill(d,10,scaleFactor);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(!Muon && Electron && dilepton)   //Muon-Electron Selection Table
            {
                if(diLepMass < 20 || (diLepMass > (ZMass-ZMassWindow) && diLepMass < (ZMass-ZMassWindow))) ZVeto = true;
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,scaleFactor);
                    if (nEl==2 && nMu ==0)
                    {
                        selecTable.Fill(d,2,scaleFactor);
                        if(!ZVeto)
                        {
                            selecTable.Fill(d,3,scaleFactor);
                            if(nJets>=4)
                            {
                                selecTable.Fill(d,4,scaleFactor);
                                if(nMtags>=1)
                                {
                                    selecTable.Fill(d,5,scaleFactor);
                                    if(nMtags>=2)
                                    {
                                        selecTable.Fill(d,6,scaleFactor);
                                        if(nJets==5)
                                        {
                                            selecTable.Fill(d,7,scaleFactor);
                                            if(nJets==6)
                                            {
                                                selecTable.Fill(d,8,scaleFactor);
                                                if(nJets==7)
                                                {
                                                    selecTable.Fill(d,9,scaleFactor);
                                                    if(nJets==8)
                                                    {
                                                        selecTable.Fill(d,10,scaleFactor);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            /////////////////////////////////
            // Applying baseline selection //
            /////////////////////////////////

            //Filling Histogram of the number of vertices before Event Selection

//            if (!trigged) continue;  // Redunant check that an HLT was triggered
            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
////            if (!(selectedJets.size() >= 6)) continue; //Selection of a minimum of 6 Jets in Event
//
//          cout <<"Number of Muons, Electrons, Jets, BJets, JetCut, MuonChannel, ElectronChannel ===>  "<< nMu <<"  "  <<nEl<<" "<< selectedJets.size()   <<"  " <<  nLtags   <<"  "<<JetCut  <<"  "<<Muon<<" "<<Electron<<endl;


            if (debug)	cout <<" applying baseline event selection..."<<endl;
            //Apply the lepton, btag and HT selections
            if (Muon && Electron && dilepton)
            {
                if  (  !( nMu == 1 && nEl == 1 )) continue; // Muon-Electron Channel Selection
            }
            else if (Muon && !Electron && dilepton)
            {
                if  (  !( nMu == 2 && nEl == 0 && !ZVeto)) continue; // Muon-Electron Channel Selection
            }
            else if (!Muon && Electron && dilepton)
            {
                if  (  !( nMu == 0 && nEl == 2 && !ZVeto)) continue; // Muon-Electron Channel Selection
            }
            else
            {
                cerr<<"Correct Channel not selected."<<endl;
                exit(1);
            }
            sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());

            //Scan for best MET Cut before Jet requirements
            for(int metCut = 0; metCut < 100; metCut = metCut + 5)
            {
                if(mets[0]->Et() > metCut) MSPlot["METCutEff"]->Fill(metCut, datasets[d], true, Luminosity*scaleFactor );
            }



            if (dilepton && Muon && Electron)
            {
                if (!(nJets>=4 && nMtags >=2 )) continue; //Jet Tag Event Selection Requirements for Mu-El dilepton channel
//                if (!(temp_HT >= 400)) continue; //Jet Tag Event Selection Requirements for Mu-El dilepton channel
            }
            else if (dilepton && Muon && !Electron)
            {
                if (!(nJets>=4 && nMtags >=2 )) continue; //Jet Tag Event Selection Requirements for Mu-El dilepton channel
//                if (!(temp_HT >= 400)) continue; //Jet Tag Event Selection Requirements for Mu-El dilepton channel
            }
            else if (dilepton && !Muon && Electron)
            {
                if (!(nJets>=4 && nMtags >=2 )) continue; //Jet Tag Event Selection Requirements for Mu-El dilepton channel
//                if (!(temp_HT >= 400)) continue; //Jet Tag Event Selection Requirements for Mu-El dilepton channel
            }
            if(debug)
            {
                cout<<"Selection Passed."<<endl;
                cin.get();
            }
            passed++;

            vector<TLorentzVector*> selectedMuonTLV_JC;
            selectedMuonTLV_JC.push_back(selectedMuons[0]);

            ///////////////////////
            // Getting Gen Event //
            ///////////////////////

            TRootGenEvent* genEvt = 0;

            if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "Data")
            {
                vector<TRootMCParticle*> mcParticles;
                vector<TRootMCParticle*> mcTops;
                mcParticlesMatching_.clear();
                mcParticlesTLV.clear();
                //selectedJetsTLV.clear();
                mcParticles.clear();
                mcTops.clear();

                int leptonPDG, muonPDG = 13, electronPDG = 11;
                leptonPDG = muonPDG;

                genEvt = treeLoader.LoadGenEvent(ievt,false);
                treeLoader.LoadMCEvent(ievt, genEvt, 0, mcParticlesMatching_,false);
                if (debug) cout <<"size   "<< mcParticlesMatching_.size()<<endl;
            }

            //////////////////////////////////////
            // MVA Hadronic Top Reconstructions //
            //////////////////////////////////////

            jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticlesMatching_, selectedJets, selectedMuonTLV_JC[0], genEvt, scaleFactor);
            double TriJetMass, DiJetMass;
            vector<TRootPFJet*> MVASelJets1;
            TLorentzVector Wh, Bh, Th;
            int wj1;
                int wj2;
                int bj1;

            if(!TrainMVA)
            {
                MVAvals1 = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value
                MSPlot["MVA1TriJet"]->Fill(MVAvals1.first, datasets[d], true, Luminosity*scaleFactor );
                topness = MVAvals1.first;
                for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
                {
                    if (seljet1 == MVAvals1.second[0] || seljet1 == MVAvals1.second[1] || seljet1 == MVAvals1.second[2])
                    {
                        MVASelJets1.push_back(selectedJets[seljet1]);
                    }

                }

                //check data-mc agreement of kin. reco. variables.
                float mindeltaR =100.;
                float mindeltaR_temp =100.;


                //define the jets from W as the jet pair with smallest deltaR
                for (int m=0; m<MVASelJets1.size(); m++)
                {
                    for (int n=0; n<MVASelJets1.size(); n++)
                    {
                        if(n==m) continue;
                        TLorentzVector lj1  = *MVASelJets1[m];
                        TLorentzVector lj2  = *MVASelJets1[n];
                        mindeltaR_temp  = lj1.DeltaR(lj2);
                        if (mindeltaR_temp < mindeltaR)
                        {
                            mindeltaR = mindeltaR_temp;
                            wj1 = m;
                            wj2 = n;
                        }
                    }
                }
                // find the index of the jet not chosen as a W-jet
                for (unsigned int p=0; p<MVASelJets1.size(); p++)
                {
                    if(p!=wj1 && p!=wj2) bj1 = p;
                }

                if (debug) cout <<"Processing event with jetcombiner : 3 "<< endl;

                //now that putative b and W jets are chosen, calculate the six kin. variables.
                Wh = *MVASelJets1[wj1]+*MVASelJets1[wj2];
                Bh = *MVASelJets1[bj1];
                Th = Wh+Bh;

                TriJetMass = Th.M();

                DiJetMass = Wh.M();
                //DeltaR
                float AngleThWh = fabs(Th.DeltaPhi(Wh));
                float AngleThBh = fabs(Th.DeltaPhi(Bh));

                float btag = MVASelJets1[bj1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();

                double PtRat = (( *MVASelJets1[0] + *MVASelJets1[1] + *MVASelJets1[2] ).Pt())/( MVASelJets1[0]->Pt() + MVASelJets1[1]->Pt() + MVASelJets1[2]->Pt() );
                double diLepThdR = fabs(Th.DeltaR(diLep));
                double diLepThdPhi = fabs(Th.DeltaPhi(diLep));
                if (debug) cout <<"Processing event with jetcombiner : 4 "<< endl;

                MSPlot["MVA1TriJetMass"]->Fill(TriJetMass,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1DiJetMass"]->Fill(DiJetMass,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1BTag"]->Fill(btag,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1PtRat"]->Fill(PtRat,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1AnThWh"]->Fill(AngleThWh,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1AnThBh"]->Fill(AngleThBh,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1dRThDiLep"]->Fill(diLepThdR,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1dPhiThDiLep"]->Fill(diLepThdPhi,  datasets[d], true, Luminosity*scaleFactor );


                if (debug) cout <<"Processing event with jetcombiner : 8 "<< endl;


            }



            ///////////////////////////////////
            // Filling histograms / plotting //
            ///////////////////////////////////

            MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);




            //////////////////////
            // Muon Based Plots //
            //////////////////////
            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ )
            {
                MSPlot["MuonPt"]->Fill(selectedMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            }

            //////////////////////////
            // Electron Based Plots //
            //////////////////////////

            for (Int_t selel =0; selel < selectedElectrons.size(); selel++ )
            {
                float reliso = selectedElectrons[selel]->relPfIso(3, 0.5);
                MSPlot["ElectronRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
            }

            //////////////////////
            // Jets Based Plots //
            //////////////////////

            HT = 0;
            float HT1M2L=0, H1M2L=0, HTbjets=0, HT2M=0, H2M=0;


            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
                if(nMtags>=2 && seljet1>=2)
                {
                    jetpt = selectedJets[seljet1]->Pt();
                    HT2M = HT2M + jetpt;
                    H2M = H2M + selectedJets[seljet1]->P();
                }
                if(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.244)
                {
                    HTb += selectedJets[seljet1]->Pt();
                }
                MSPlot["BdiscBJetCand_CSV"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(),datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                //Event-level variables
                jetpt = selectedJets[seljet1]->Pt();
                HT = HT + jetpt;
                H = H +  selectedJets[seljet1]->P();
                if (seljet1 > 2  )  HTHi +=  selectedJets[seljet1]->Pt();
            }

            HTRat = HTHi/HT;
            HTH = HT/H;

            MSPlot["HTExcess2M"]->Fill(HT2M, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["HTH"]->Fill(HTH, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
            histo2D["HTLepSep"]->Fill(HT, lep1.DeltaR(lep2));
            sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order Jets wrt Pt for tuple output

            if(dilepton && Muon && Electron)
            {
                muonpt = selectedMuons[0]->Pt();
                muoneta = selectedMuons[0]->Eta();
                electronpt  = selectedElectrons[0]->Pt();
            }
            if(dilepton && Muon && !Electron)
            {
                muonpt = selectedMuons[0]->Pt();
                muoneta = selectedMuons[0]->Eta();
            }
            if(dilepton && !Muon && Electron)
            {
                muonpt = selectedElectrons[0]->Pt();
                muoneta = selectedElectrons[0]->Eta();
            }


            bjetpt= selectedLBJets[0]->Pt();

            Eventcomputer_->FillVar("topness",topness);
            Eventcomputer_->FillVar("muonpt",muonpt);
            Eventcomputer_->FillVar("muoneta",muoneta);
            Eventcomputer_->FillVar("HTH", HTH);
            Eventcomputer_->FillVar("HTRat",HTRat);
            Eventcomputer_->FillVar("HTb", HTb);
            Eventcomputer_->FillVar("nLtags",nLtags );
            Eventcomputer_->FillVar("nMtags",nMtags );
            Eventcomputer_->FillVar("nTtags",nTtags );
            Eventcomputer_->FillVar("nJets", selectedJets.size() );
            Eventcomputer_->FillVar("Jet3Pt", selectedJets[2]->Pt() );
            Eventcomputer_->FillVar("Jet4Pt", selectedJets[3]->Pt() );

            std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();

            for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it)
            {

                //  cout <<"MVA Method : "<< it->first    <<" Score : "<< it->second <<endl;
                BDTScore = it->second;

            }


            float nvertices = vertex.size();
            float normfactor = datasets[d]->NormFactor();

            ///////////////////
            //MET Based Plots//
            ///////////////////

            MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);

            //////////////////
            //Topology Plots//
            /////////////////

            vector<TLorentzVector> selectedParticlesTLV, diLepSystemTLV, topDiLepSystemTLV;
            //collection Total Event TLVs
            selectedParticlesTLV.insert(selectedParticlesTLV.end(), selectedElectronsTLV_JC.begin(), selectedElectronsTLV_JC.end());
            selectedParticlesTLV.insert(selectedParticlesTLV.end(), selectedMuonsTLV_JC.begin(), selectedMuonsTLV_JC.end());
            selectedParticlesTLV.insert(selectedParticlesTLV.end(), selectedJetsTLV.begin(), selectedJetsTLV.end());
            selectedParticlesTLV.push_back(*mets[0]);
            //collecting diLep TLVs
            diLepSystemTLV.push_back(lep1);
            diLepSystemTLV.push_back(lep2);
            diLepSystemTLV.push_back(*mets[0]);
            //collecting topDiLep TLVs
            topDiLepSystemTLV.insert(topDiLepSystemTLV.end(), diLepSystemTLV.begin(), diLepSystemTLV.end());
            topDiLepSystemTLV.push_back(*MVASelJets1[wj1]);
            topDiLepSystemTLV.push_back(*MVASelJets1[wj2]);
            topDiLepSystemTLV.push_back(*MVASelJets1[bj1]);

            MSPlot["TotalSphericity"]->Fill(Sphericity(selectedParticlesTLV), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["TotalCentrality"]->Fill(Centrality(selectedParticlesTLV), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["DiLepSphericity"]->Fill(Sphericity(diLepSystemTLV), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["DiLepCentrality"]->Fill(Centrality(diLepSystemTLV), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["TopDiLepSphericity"]->Fill(Sphericity(topDiLepSystemTLV), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["TopDiLepCentrality"]->Fill(Centrality(topDiLepSystemTLV), datasets[d], true, Luminosity*scaleFactor);



            //////////////////
            //Filling nTuple//
            //////////////////

            //	  tup->Fill(nJets,nLtags,nMtags,nTtags,HT,muonpt,muoneta,electronpt,bjetpt,HT2M,HTb,HTH,HTRat,topness,scaleFactor,nvertices,normfactor,Luminosity,weight_0);



            //                "BDT:nJets:nFatJets:nWTags:nTopTags:nLtags:nMtags:nTtags:HT:LeadingMuonPt:LeadingMuonEta:LeadingElectronPt:LeadingBJetPt:HT2L:HTb:HTH:HTRat:topness:ScaleFactor:PU:NormFactor:Luminosity:GenWeight");

            float vals[23] = {BDTScore,nJets,nFatJets,nWTags,nTopTags,nLtags,nMtags,nTtags,HT,muonpt,muoneta,electronpt,bjetpt,HT2M,HTb,HTH,HTRat,topness,scaleFactor,nvertices,normfactor,Luminosity,weight_0};

            tup->Fill(vals);


        } //End Loop on Events

        tup->Write();
        tupfile->Close();
        cout <<"n events passed  =  "<<passed <<endl;
        cout <<"n events with negative weights = "<<negWeights << endl;
        cout << "Event Count: " << eventCount << endl;
        cout << "Weight Count: " << weightCount << endl;
        //important: free memory
        treeLoader.UnLoadDataset();
    } //End Loop on Datasets

    eventlist.close();

    /////////////
    // Writing //
    /////////////

    cout << " - Writing outputs to the files ..." << endl;

    //////////////////////
    // Selection tables //
    //////////////////////

    //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
    selecTable.TableCalculator(  true, true, true, true, true);

    //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
    selecTable.Write(  outputDirectory+"/FourTop"+postfix+"_Table"+channelpostfix+".tex",    false,true,true,true,false,false,true);

    fout->cd();
    TFile *foutmva = new TFile ("foutMVA.root","RECREATE");
    cout <<" after cd .."<<endl;

    string pathPNGJetCombi = pathPNG + "JetCombination/";
    mkdir(pathPNGJetCombi.c_str(),0777);
//    if(TrainMVA)jetCombiner->Write(foutmva, true, pathPNGJetCombi.c_str());

//Output ROOT file
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin();
            it != MSPlot.end();
            it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->Write(fout, name, true, pathPNG, "pdf");
    }


    TDirectory* th2dir = fout->mkdir("Histos2D");
    th2dir->cd();


    for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
    {

        TH2F *temp = it->second;
        temp->Write();
        //TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
        //tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
    }
    delete fout;
    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

    return 0;
}

int Factorial(int N = 1)
{
    int fact = 1;
    for( int i=1; i<=N; i++ )
        fact = fact * i;  // OR fact *= i;
    return fact;
}

float Sphericity(vector<TLorentzVector> parts )
{
    if(parts.size()>0)
    {
        double spTensor[3*3] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
        int counter = 0;
        float tensorNorm = 0, y1 = 0, y2 = 0, y3 = 0;

        for(int tenx = 0; tenx < 3; tenx++)
        {
            for(int teny = 0; teny < 3; teny++)
            {
                for(int selpart = 0; selpart < parts.size(); selpart++)
                {

                    spTensor[counter] += ((parts[selpart][tenx])*(parts[selpart][teny]));
//                    if((tenx == 0 && teny == 2) || (tenx == 2 && teny == 1))
//                    {
//                    cout << "nan debug term " << counter+1 << ": " << (parts[selpart][tenx])*(parts[selpart][teny]) << endl;
//                    cout << "Tensor Building Term " << counter+1 << ": " << spTensor[counter] << endl;
//                    }
                    if(tenx ==0 && teny == 0)
                    {
                        tensorNorm += parts[selpart].Vect().Mag2();
                    }
                }
                if((tenx == 0 && teny == 2) || (tenx == 2 && teny == 1))
                {
//                    cout << "Tensor term pre-norm " << counter+1 << ": " << spTensor[counter] << endl;
                }
                spTensor[counter] /= tensorNorm;
//                cout << "Tensor Term " << counter+1 << ": " << spTensor[counter] << endl;
                counter++;
            }
        }
        TMatrixDSym m(3, spTensor);
        //m.Print();
        TMatrixDSymEigen me(m);
        TVectorD eigenval = me.GetEigenValues();
        vector<float> eigenVals;
        eigenVals.push_back(eigenval[0]);
        eigenVals.push_back(eigenval[1]);
        eigenVals.push_back(eigenval[2]);
        sort(eigenVals.begin(), eigenVals.end());
        //cout << "EigenVals: "<< eigenVals[0] << ", " << eigenVals[1] << ", " << eigenVals[2] << ", " << endl;
        float sp = 3.0*(eigenVals[0] + eigenVals[1])/2.0;
        //cout << "Sphericity: " << sp << endl;
        return sp;
    }
    else
    {
        return 0;
    }
}
float Centrality(vector<TLorentzVector> parts)
{
    float E = 0, ET = 0;
    for(int selpart = 0; selpart < parts.size(); selpart++)
    {
        E += parts[selpart].E();
        ET += parts[selpart].Et();
    }
    return ET/E;
}
