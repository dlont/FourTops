//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.            ////
//////////////////////////////////////////////////////////////////////////////

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
//#include "TopTreeAnalysisBase/Selection/interface/FourTopSelectionTable.h"
#include "TopTreeAnalysisBase/Selection/interface/Run2Selection.h"

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
#include "TopTreeAnalysisBase/Tools/interface/SourceDate.h"

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


pair<float, vector<unsigned int> > MVAvals1;
pair<float, vector<unsigned int> > MVAvals2;
pair<float, vector<unsigned int> > MVAvals2ndPass;
pair<float, vector<unsigned int> > MVAvals3rdPass;

int nMVASuccesses=0;
int nMatchedEvents=0;

/// Normal Plots (TH1F* and TH2F*)

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;
/*
struct HighestTCHEBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_trackCountingHighEffBJetTags() > j2->btag_trackCountingHighEffBJetTags();
    }
};
*/
struct HighestCVSBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_combinedInclusiveSecondaryVertexV2BJetTags() > j2->btag_combinedInclusiveSecondaryVertexV2BJetTags();
    }
};


bool match;

int main (int argc, char *argv[])
{

    //Checking Passed Arguments to ensure proper execution of MACRO
    if(argc != 14)
    {
        std::cerr << "INVALID INPUT FROM XMLFILE.  CHECK XML IMPUT FROM SCRIPT.  " << argc << " ARGUMENTS HAVE BEEN PASSED." << std::endl;
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
    const int startEvent            = strtol(argv[12], NULL, 10);
    const int endEvent              = strtol(argv[13], NULL, 10);
    vector<string> vecfileNames;
    vecfileNames.push_back(fileName);

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
    cout << "EDITS" <<endl;
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

    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the FourTop search ! "           << endl;
    cout << "*************************************************************" << endl;

    string postfix = "_Run2_TopTree_Study"; // to relabel the names of the output file

    ///////////////////////////////////////
    // Configuration
    ///////////////////////////////////////

    string channelpostfix = "";
    string xmlFileName = "";

    //Setting Lepton Channels (Setting both flags true will select Muon-Electron Channel when dilepton is also true)
    bool dilepton = false;
    bool SingleLepton = true;

    bool Muon = true;
    bool Electron = false;


    if(Muon && Electron && dilepton)
    {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_MuEl";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else if(Muon && SingleLepton){
        cout<<" ***** USING SINGLE MUON CHANNEL  ******"<<endl;
        channelpostfix = "_Mu";
        xmlFileName = "config/Run2_SingleLeptonsamples.xml";
    }
    else if(Electron && SingleLepton){
        cout<<" ***** Using SINGLE ELECTRON CHANNEL *****"<<endl;
        channelpostfix = "_El";
        xmlFileName = "config/Run2_SingleLeptonsamples.xml";
    }
    else
    {
        cerr<<"Correct Channel not selected."<<endl;
        exit(1);
    }

    const char *xmlfile = xmlFileName.c_str();
    cout << "used config file: " << xmlfile << endl;

    bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff


    /////////////////////////////
    //  Set up AnalysisEnvironment
    /////////////////////////////

    AnalysisEnvironment anaEnv;
    cout<<" - Creating environment ..."<<endl;
//    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    anaEnv.PrimaryVertexCollection = "PrimaryVertex";
    anaEnv.JetCollection = "PFJets_slimmedJets";
    anaEnv.METCollection = "PFMET_slimmedMETs";
    anaEnv.MuonCollection = "Muons_slimmedMuons";
    anaEnv.ElectronCollection = "Electrons_slimmedElectrons";
    anaEnv.GenJetCollection   = "GenJets_slimmedGenJets";
    anaEnv.TrackMETCollection = "";
    anaEnv.GenEventCollection = "GenEvent";
    anaEnv.NPGenEventCollection = "NPGenEvent";
    anaEnv.MCParticlesCollection = "MCParticles";
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

    ////////BDT stuff/////////

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

    /////////end BDT stuff ////////////

    string dataSetName;
    
    string MVAmethod = "BDT"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)

    JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, false);
    cout <<"Instantiated jet combiner..."<<endl;
    double bestTopMass1 =0.;
    double bestTopMass2 = 0.;
    double bestTopMass2ndPass = 0.;
    double bestTopPt =0.;
    float topness;


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
    string rootFileName ("FourTop"+postfix+"_"+dName+channelpostfix+".root");
    TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

    //vector of objects
    cout << " - Variable declaration ..." << endl;
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
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
    MSPlot["NbOfSelectedBJets"]                             = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of tags");
    MSPlot["HTb_SelectedJets"]                              = new MultiSamplePlot(datasets, "HTb_SelectedJets", 50, 0, 1500, "HTb");    
    //Jets
    MSPlot["JetEta"]                                        = new MultiSamplePlot(datasets, "JetEta", 40,-4, 4, "Jet #eta");
    MSPlot["HT_SelectedJets"]                               = new MultiSamplePlot(datasets, "HT_SelectedJets", 30, 0, 1500, "HT");
    MSPlot["HTRat"]                                         = new MultiSamplePlot(datasets, "HTRat", 50, 0, 20, "HTRat");
    MSPlot["5thJetPt"]                                      = new MultiSamplePlot(datasets, "5thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["6thJetPt"]                                      = new MultiSamplePlot(datasets, "6thJetPt", 60, 0, 400, "PT_{jet}");
    //MSPlot["7thJetPt"]                                      = new MultiSamplePlot(datasets, "7thJetPt", 60, 0, 400, "PT_{jet}");
    //MET
    MSPlot["MET"]                                           = new MultiSamplePlot(datasets, "MET", 70, 0, 700, "MET");

    MSPlot["NbOfBadTrijets"]                  = new MultiSamplePlot(datasets, "NbOfBadTriJets", 150, 0, 150, "Nb. of Bad Combs");
    
    MSPlot["MVA1TriJetMass"] = new MultiSamplePlot(datasets, "MVA1TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1DiJetMass"] = new MultiSamplePlot(datasets, "MVA1DiJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1PtRat"] = new MultiSamplePlot(datasets, "MVA1PtRat", 25, 0, 2, "P_{t}^{Rat}");
    MSPlot["MVA1BTag"] = new MultiSamplePlot(datasets, "MVA1BTag", 35, 0, 1, "BTag");
    MSPlot["MVA1AnThBh"] = new MultiSamplePlot(datasets, "MVA1AnThBh", 35, 0, 3.14, "AnThBh");
    MSPlot["MVA1AnThWh"] = new MultiSamplePlot(datasets, "MVA1AnThWh", 35, 0, 3.14, "AnThWh");

    MSPlot["TriJetMass_Matched"] = new MultiSamplePlot(datasets, "TriJetMassMatched", 100, 0, 1000, "m_{bjj}");
    MSPlot["TriJetMass_UnMatched"] = new MultiSamplePlot(datasets, "TriJetMassUnMatched", 100, 0, 1000, "m_{bjj}");
    MSPlot["Chi2_Matched"] = new MultiSamplePlot(datasets, "Chi2Matched", 100, 0, 100, "#chi^{2}");
    MSPlot["Chi2_UnMatched"] = new MultiSamplePlot(datasets, "Chi2UnMatched", 100, 0, 100, "#chi^{2}");
 
    MSPlot["MVA1TriJet"] = new MultiSamplePlot(datasets, "MVA1TriJet", 35, -1., 0.5, "BDT Discriminator");

    MSPlot["MVA2TriJetMass"] = new MultiSamplePlot(datasets, "MVA2TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA2ndPassTriJetMass"] = new MultiSamplePlot(datasets, "MVA2ndPassTriJetMass", 30, 0, 1000, "m_{bjj}");
    MSPlot["MultiTopness"] = new MultiSamplePlot(datasets, "MultiTopness", 35, -1., 0.5, "MultiTopness");

    MSPlot["MVA1TriJetMassMatched"] = new MultiSamplePlot(datasets, "MVA1TriJetMassMatched", 75, 0, 500, "m_{bjj}");

    ///////////////////
    // 1D histograms
    ///////////////////

    //Plots
    string pathPNG = "FourTop"+postfix+channelpostfix;
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
            CutsselecTable.push_back(string("At least 1 Loose Isolated Muon"));
            CutsselecTable.push_back(string("At least 1 Loose Electron"));
            CutsselecTable.push_back(string("At least 4 Jets"));
            CutsselecTable.push_back(string("At least 1 CSVM Jet"));
            CutsselecTable.push_back(string("At least 2 CSVM Jets"));
            CutsselecTable.push_back(string("HT $\\geq 100 GeV$"));
            CutsselecTable.push_back(string("HT $\\geq 200 GeV$"));
            CutsselecTable.push_back(string("HT $\\geq 300 GeV$"));
            CutsselecTable.push_back(string("HT $\\geq 400 GeV$"));
        }
    }
    else if(SingleLepton)
    {
        /////////////////////////////////
        // Selection table: SingleLepton + jets
        /////////////////////////////////
        if(Muon)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("At least 1 Loose Isolated Muon"));
            CutsselecTable.push_back(string("At least 6 Jets"));
            CutsselecTable.push_back(string("At least 1 CSVL Jet"));
            CutsselecTable.push_back(string("At least 2 CSVL Jets"));
            CutsselecTable.push_back(string("HT $\\geq 100 GeV$"));
            CutsselecTable.push_back(string("HT $\\geq 200 GeV$"));
            CutsselecTable.push_back(string("HT $\\geq 300 GeV$"));
            CutsselecTable.push_back(string("HT $\\geq 400 GeV$"));
        }

        else if(Electron)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("At least 1 Loose Electron"));
            CutsselecTable.push_back(string("At least 6 Jets"));
            CutsselecTable.push_back(string("At least 1 CSVL Jet"));
            CutsselecTable.push_back(string("At least 2 CSVL Jets"));
            CutsselecTable.push_back(string("HT $\\geq 100 GeV$"));
            CutsselecTable.push_back(string("HT $\\geq 200 GeV$"));
            CutsselecTable.push_back(string("HT $\\geq 300 GeV$"));
            CutsselecTable.push_back(string("HT $\\geq 400 GeV$"));
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

        SourceDate *strdate = new SourceDate();
        string date_str = strdate->ReturnDateStr();

        string dataSetName = datasets[d]->Name();
        string channel_dir = "Craneens"+channelpostfix;
        string date_dir = channel_dir+"/Craneens" + date_str +"/";
        int mkdirstatus = mkdir(channel_dir.c_str(),0777);
        mkdirstatus = mkdir(date_dir.c_str(),0777);

        //     string Ntupname = "Craneens/Craneen_" + dataSetName +postfix + "_" + date_str+  ".root";

        string Ntupname = "Craneens"+channelpostfix+"/Craneens"+ date_str  +"/Craneen_" + dataSetName +postfix + ".root";
        string Ntuptitle = "Craneen_" + channelpostfix;

        TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");

        //TNtuple * tup = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(),"HT:nJets:nTags:nLights:5thJetPt:6thJetPt:HTRat:HTb:HT2L2J:Multitopness:ScaleFactor:NormFactor:Luminosity");
        TNtuple * tup = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(),"BDT:nJets:nLtags:nMtags:nTtags:HT:LeadingMuonPt:LeadingMuonEta:LeadingBJetPt:HT2M:HTb:HTH:HTRat:topness:ScaleFactor:PU:NormFactor:Luminosity:GenWeight");


        //////////////////////////////////////////////////
        // Loop on events
        /////////////////////////////////////////////////

        int itrigger = -1, previousRun = -1;

        int start = 0;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();

        cout <<"Number of events = "<<  ending  <<endl;

        //int event_start = 0;
        int event_start = startEvent;

        if (dataSetName == "Data")TrainMVA=false;

        if (verbose > 1) cout << " - Loop over events " << endl;

        //double MHT, MHTSig, STJet, EventMass, EventMassX , SumJetMass, SumJetMassX,H,HX , HT, HTX,HTH,HTXHX, sumpx_X, sumpy_X, sumpz_X, sume_X, sumpx, sumpy, sumpz, sume, jetpt,PTBalTopEventX,PTBalTopSumJetX , PTBalTopMuMet;
        float BDTScore, MHT, MHTSig, STJet,muoneta, muonpt,electronpt,bjetpt, EventMass, EventMassX , SumJetMass, SumJetMassX,H,HX ,HTHi,HTRat, HT, HTX,HTH,HTXHX, sumpx_X, sumpy_X, sumpz_X, sume_X, sumpx, sumpy, sumpz, sume, jetpt,PTBalTopEventX,PTBalTopSumJetX , PTBalTopMuMet;

        double currentfrac =0.;
        double end_d = ending;

        if(endEvent > ending)
            end_d = ending;
        else
            end_d = endEvent;

        cout <<"Will run over "<<  end_d<< " events..."<<endl;
        cout <<"Starting event = = = = "<< event_start  << endl;

        //define object containers
        vector<TRootElectron*> selectedElectrons;
        vector<TRootPFJet*>    selectedJets;
        vector<TRootPFJet*>    selectedJets2; //after removal of 2 highest CSVL btags
        vector<TRootMuon*>     selectedMuons;
        vector<TRootElectron*> selectedExtraElectrons;
        vector<TRootMuon*>     selectedExtraMuons;
        selectedElectrons.reserve(10);
        selectedMuons.reserve(10);
        vector<TRootPFJet*>      selectedJets2ndPass;
        vector<TRootPFJet*>      selectedJets3rdPass;
        vector<TRootPFJet*>      MVASelJets1;
        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////
        //end_d = 10;
        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
            //MHT = 0.,MHTSig = 0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0.  ,H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0.;
            BDTScore= -99999.0, MHT = 0.,MHTSig = 0.,muoneta = 0., muonpt =0., electronpt=0., bjetpt =0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0., HTHi =0., HTRat = 0;  H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0.;

            double ievt_d = ievt;
            currentfrac = ievt_d/end_d;
            if (debug)cout <<"event loop 1"<<endl;

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);  //load event
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


            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////

            // Apply trigger selection
            // trigged = treeLoader.EventTrigged (itrigger);
            bool trigged = true;  // Disabling the HLT requirement
            if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;
            if(!trigged)		   continue;  //If an HLT condition is not present, skip this event in the loop.
            // Declare selection instance
            Run2Selection r2selection(init_jets, init_muons, init_electrons, mets);
            // Define object selection cuts
            if (Muon && Electron && dilepton)
            {


                if (debug)cout<<"Getting Jets"<<endl;
                selectedJets                                        = r2selection.GetSelectedJets(); // ApplyJetId
                if (debug)cout<<"Getting Tight Muons"<<endl;
                selectedMuons                                       = r2selection.GetSelectedMuons();
                if (debug)cout<<"Getting Loose Electrons"<<endl;
                selectedElectrons                                   = r2selection.GetSelectedElectrons("Tight", "PHYS14", true); // VBTF ID
            }

            else if(SingleLepton && Electron){

                if (debug)cout<<"Getting Jets"<<endl;
                selectedJets                                        = r2selection.GetSelectedJets(); // ApplyJetId
                if (debug)cout<<"Getting Tight Muons"<<endl;
                selectedMuons                                       = r2selection.GetSelectedMuons();
                if (debug)cout<<"Getting Tight Electrons"<<endl;
                selectedElectrons                                   = r2selection.GetSelectedElectrons("Tight", "PHYS14", true); // VBTF ID                       
                if (debug)cout<<"Getting Loose Electrons"<<endl;
                selectedExtraElectrons                              = r2selection.GetSelectedElectrons("Loose", "PHYS14", true);
            }
            else if(SingleLepton && Muon){

                if (debug)cout<<"Getting Jets"<<endl;
                selectedJets                                        = r2selection.GetSelectedJets(); // ApplyJetId
                if (debug)cout<<"Getting Tight Muons"<<endl;
                selectedMuons                                       = r2selection.GetSelectedMuons();
                if (debug)cout<<"Getting Tight Electrons"<<endl;
                selectedElectrons                                   = r2selection.GetSelectedElectrons("Tight", "PHYS14",true); // VBTF ID    
                if (debug)cout<<"Getting Loose Muons"<<endl;
                selectedExtraMuons                                  = r2selection.GetSelectedMuons(20, 2.4, 0.20);                                   
            }

            vector<TRootJet*>      selectedMBJets;
            vector<TRootJet*>      selectedLBJets;
            vector<TRootJet*>      selectedLightJets;
            vector<TRootJet*>      selectedTBJets;


            int JetCut =0;
            int nMu, nEl, nLooseMu, nLooseEl;
            if(dilepton && Muon && Electron)
            {
                nMu = selectedMuons.size(); //Number of Muons in Event
                nEl = selectedElectrons.size(); //Number of Electrons in Event
            }
            else if(SingleLepton && Muon){
                nMu = selectedMuons.size(); //Number of Muons in Event
                nEl = selectedElectrons.size(); //Number of Electrons in Event   
                nLooseMu = selectedExtraMuons.size();   //Number of loose muons      
            }
            else if(SingleLepton && Electron){
                nMu = selectedMuons.size(); //Number of Muons in Event
                nEl = selectedElectrons.size(); //Number of Electrons in Event
                nLooseEl = selectedExtraElectrons.size(); //Number of loose muons
            }
            bool isTagged =false;

            ///////////////////////////////////////////////////////////////////////////////////
            // Preselection looping over Jet Collection                                      //
            // Summing HT and calculating leading, lagging, and ratio for Selected and BJets //
            ///////////////////////////////////////////////////////////////////////////////////
            float HTb = 0.;
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
                if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.244   )
                {
                    selectedLBJets.push_back(selectedJets[seljet]);
                    if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.679)
                    {
                        selectedMBJets.push_back(selectedJets[seljet]);
                        HTb += selectedJets[seljet]->Pt();

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

            float temp_HT = 0.;
            double HT_leading = 0.;
            double HT_lagging = 0.;
            double HTRat = 0;
            for (Int_t seljet0 =0; seljet0 < selectedJets.size(); seljet0++ ){
              temp_HT += selectedJets[seljet0]->Pt();
              if (seljet0 < 4){
                HT_leading += selectedJets[seljet0]->Pt();
            }else{

                HT_lagging += selectedJets[seljet0]->Pt();
            }
        }

        HTRat = HT_leading/HT_lagging;

            //////////////////////
            // Sync'ing cutflow //
            //////////////////////

        if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;
            // Apply primary vertex selection
        bool isGoodPV = r2selection.isPVSelected(vertex, 4, 24., 2);
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
                    if (nMu>=1)
                    {
                        selecTable.Fill(d,2,scaleFactor);
                        if(nEl>=1)
                        {
                            selecTable.Fill(d,3,scaleFactor);
                            if(nJets>=4)
                            {
                                selecTable.Fill(d,4,scaleFactor);
                                if(nLtags>=1)
                                {
                                    selecTable.Fill(d,5,scaleFactor);
                                    if(nLtags>=2)
                                    {
                                        selecTable.Fill(d,6,scaleFactor);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if(Muon && SingleLepton)   //Muon-Electron Selection Table
            {
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,scaleFactor);
                    if (nMu>=1)
                    {
                        selecTable.Fill(d,2,scaleFactor);
                        if(nEl>=1)
                        {
                            selecTable.Fill(d,3,scaleFactor);
                            if(nJets>=4)
                            {
                                selecTable.Fill(d,4,scaleFactor);
                                if(nLtags>=1)
                                {
                                    selecTable.Fill(d,5,scaleFactor);
                                    if(nLtags>=2)
                                    {
                                        selecTable.Fill(d,6,scaleFactor);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if(Electron && SingleLepton)   //Muon-Electron Selection Table
            {
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,scaleFactor);
                    if (nMu>=1)
                    {
                        selecTable.Fill(d,2,scaleFactor);
                        if(nEl>=1)
                        {
                            selecTable.Fill(d,3,scaleFactor);
                            if(nJets>=4)
                            {
                                selecTable.Fill(d,4,scaleFactor);
                                if(nLtags>=1)
                                {
                                    selecTable.Fill(d,5,scaleFactor);
                                    if(nLtags>=2)
                                    {
                                        selecTable.Fill(d,6,scaleFactor);
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

            // if (!trigged) continue;  // Redunant check that an HLT was triggered
            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
            ////  if (!(selectedJets.size() >= 6)) continue; //Selection of a minimum of 6 Jets in Event

            if (debug) cout <<"Number of Muons, Electrons, Jets, BJets, JetCut, MuonChannel, ElectronChannel ===>  "<< nMu <<"  "  <<nEl<<" "<< selectedJets.size()   <<"  " <<  nLtags   <<"  "<<JetCut  <<"  "<<Muon<<" "<<Electron<<endl;

            if (debug)	cout <<" applying baseline event selection..."<<endl;
            //Apply the lepton, btag and HT selections
            if (Muon && Electron && dilepton)
            {
                if  (  !( nMu >= 1 && nEl >= 1 )) continue; // Muon-Electron Channel Selection
            }
            else if (SingleLepton && Muon)
            {
                if  (  !( nMu == 1 && nEl == 0 && nLooseMu == 1)) continue; // Muon Channel Selection
            }
            else if(SingleLepton && Electron){
                if  (  !( nMu == 0 && nEl == 1 && nLooseEl == 1)) continue; // Electron Channel Selection
            }
            else{
                cerr<<"Correct Channel not selected."<<endl;
                exit(1);
            }

            sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());

            if (dilepton && Muon && Electron)
            {
                if (!(nJets>=4 && nLtags >=2 )) continue; //Jet Tag Event Selection Requirements for Mu-El dilepton channel
                //if (!(temp_HT >= 400)) continue; //Jet Tag Event Selection Requirements for Mu-El dilepton channel
            }
            else if (SingleLepton && Muon)
            {
                if (!(nJets>=6 && nLtags >=2 )) continue; //Jet Tag Event Selection Requirements for Mu dilepton channel
                //if (!(temp_HT >= 400)) continue; //Jet Tag Event Selection Requirements for Mu dilepton channel
            }
            else if (SingleLepton && Electron)
            {
                if (!(nJets>=6 && nLtags >=2 )) continue; //Jet Tag Event Selection Requirements for El dilepton channel
                //if (!(temp_HT >= 400)) continue; //Jet Tag Event Selection Requirements for El dilepton channel
            }
            if(debug)
            {
                cout<<"Selection Passed."<<endl;
                cin.get();
            }
            passed++;

            /////////////////////////////////
            /// Find indices of jets from Tops
            ////////////////////////////////


            vector<TLorentzVector*> selectedMuonTLV_JC;
            selectedMuonTLV_JC.push_back(selectedMuons[0]);        

             
            vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV, mcMuonsTLV;
            vector<TRootMCParticle*> mcParticlesMatching_;
            bool muPlusFromTop = false, muMinusFromTop = false;
            bool elPlusFromTop = false, elMinusFromTop = false;
            pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
            leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
            vector<TRootMCParticle*> mcParticles_flav;
            TRootGenEvent* genEvt_flav = 0;

            for(unsigned int i=0; i<selectedJets.size(); i++) selectedJetsTLV.push_back(*selectedJets[i]);
                JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);

                vector< pair<unsigned int, unsigned int> > JetPartonPair;
            for(unsigned int i=0; i<mcParticlesTLV.size(); i++) //loop through mc particles and find matched jets
            {
            int matchedJetNumber = matching.getMatchForParton(i, 0);
                if(matchedJetNumber != -1)
                    JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );// Jet index, MC Particle index
            }
        
            if (debug) cout <<"n sel jets  "<<selectedJets.size()  << "   n mc particles tlv : "<< mcParticlesTLV.size() << " jet parton pari size :   "<< JetPartonPair.size()<<"  "<< muPlusFromTop<<muMinusFromTop<<endl;
        
            for(unsigned int i=0; i<JetPartonPair.size(); i++)//looping through matched jet-parton pairs
            {
                unsigned int j = JetPartonPair[i].second;     //get index of matched mc particle

                if( fabs(mcParticlesMatching_[j]->type()) < 5 )
                {
                    if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -6 )
                    || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == 6 ) )
                    {
                        if(hadronicWJet1_.first == 9999)
                        {
                            hadronicWJet1_ = JetPartonPair[i];
                            // MCPermutation[0] = JetPartonPair[i].first;
                        }
                        else if(hadronicWJet2_.first == 9999)
                        {
                            hadronicWJet2_ = JetPartonPair[i];
                            //MCPermutation[1] = JetPartonPair[i].first;
                        }
                        //else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
                    }
                }
                else if( fabs(mcParticlesMatching_[j]->type()) == 5 )
                {

                    if(  ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -6) || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 6 ) )
                    {
                        hadronicBJet_ = JetPartonPair[i];
                        //MCPermutation[2] = JetPartonPair[i].first;
                    }
                    else if((muPlusFromTop && mcParticlesMatching_[j]->motherType() == 6) ||  ( muMinusFromTop &&mcParticlesMatching_[j]->motherType() == -6) )
                    {
                        leptonicBJet_ = JetPartonPair[i];
                        //MCPermutation[3] = JetPartonPair[i].first;
                    }
                }
            }

            //  cout <<"  "<<endl;
            if (debug) cout <<"Indices of matched jets are :  "<< hadronicBJet_.first<<"  "<< hadronicWJet1_.first  <<" " << hadronicWJet2_.first <<endl;
            
            /////////////////////////////////
            /// TMVA for mass Reconstruction
            ////////////////////////////////

            jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticles_flav, selectedJets, selectedMuonTLV_JC[0], genEvt_flav, scaleFactor);
            double TriJetMass, DiJetMass;
            if (debug) cout <<"Processing event with jetcombiner :  "<< endl;

            float MultiTopness;

            if(!TrainMVA){
                MVAvals1 = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value

                MSPlot["MVA1TriJet"]->Fill(MVAvals1.first   ,  datasets[d], true, Luminosity*scaleFactor );

                selectedJets2ndPass.clear();
                selectedJets3rdPass.clear();
                MVASelJets1.clear();    

                if (debug) cout <<"Processing event with jetcombiner : 1 "<< endl;

                //make vector of jets excluding thise selected by 1st pass of mass reco
                for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ){
                    if (seljet1 == MVAvals1.second[0] || seljet1 == MVAvals1.second[1] || seljet1 == MVAvals1.second[2]){ 
                        MVASelJets1.push_back(selectedJets[seljet1]);
                        continue;
                    }
                    selectedJets2ndPass.push_back(selectedJets[seljet1]);
                }

                jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticles_flav, selectedJets2ndPass, selectedMuonTLV_JC[0], genEvt_flav, scaleFactor);
                MVAvals2ndPass = jetCombiner->getMVAValue(MVAmethod, 1);

                MultiTopness = MVAvals2ndPass.first;
                topness = MultiTopness;

                MSPlot["MultiTopness"]->Fill(MultiTopness,  datasets[d], true, Luminosity*scaleFactor );

                if (debug) cout <<"Processing event with jetcombiner : 2 "<< endl;

                //check data-mc agreement of kin. reco. variables.
                float mindeltaR =100.;
                float mindeltaR_temp =100.;
                int wj1;
                int wj2;
                int bj1;

                //define the jets from W as the jet pair with smallest deltaR
                for (int m=0; m<MVASelJets1.size(); m++) {
                    for (int n=0; n<MVASelJets1.size(); n++) {
                        if(n==m) continue;
                        TLorentzVector lj1  = *MVASelJets1[m];
                        TLorentzVector lj2  = *MVASelJets1[n];
                        mindeltaR_temp  = lj1.DeltaR(lj2);
                        if (mindeltaR_temp < mindeltaR){
                            mindeltaR = mindeltaR_temp;
                            wj1 = m;
                            wj2 = n;
                        }
                    }
                }
                // find the index of the jet not chosen as a W-jet
                for (unsigned int p=0; p<MVASelJets1.size(); p++) {
                    if(p!=wj1 && p!=wj2) bj1 = p;
                }

                if (debug) cout <<"Processing event with jetcombiner : 3 "<< endl;

                //now that putative b and W jets are chosen, calculate the six kin. variables.
                TLorentzVector Wh = *MVASelJets1[wj1]+*MVASelJets1[wj2];
                TLorentzVector Bh = *MVASelJets1[bj1];
                TLorentzVector Th = Wh+Bh;

                TriJetMass = Th.M();

                DiJetMass = Wh.M();
                    //DeltaR
                float AngleThWh = fabs(Th.DeltaPhi(Wh));
                float AngleThBh = fabs(Th.DeltaPhi(Bh));

                float btag = MVASelJets1[bj1]->btag_combinedSecondaryVertexBJetTags();

                double PtRat = (( *MVASelJets1[0] + *MVASelJets1[1] + *MVASelJets1[2] ).Pt())/( MVASelJets1[0]->Pt() + MVASelJets1[1]->Pt() + MVASelJets1[2]->Pt() );
                if (debug) cout <<"Processing event with jetcombiner : 4 "<< endl;    

                MSPlot["MVA1TriJetMass"]->Fill(TriJetMass,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1DiJetMass"]->Fill(DiJetMass,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1BTag"]->Fill(btag,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1PtRat"]->Fill(PtRat,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1AnThWh"]->Fill(AngleThWh,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA1AnThBh"]->Fill(AngleThBh,  datasets[d], true, Luminosity*scaleFactor );


                if (debug) cout <<"Processing event with jetcombiner : 8 "<< endl;    

                bestTopMass1 =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]).M();
                bestTopMass2ndPass =( *selectedJets[MVAvals2ndPass.second[0]] + *selectedJets[MVAvals2ndPass.second[1]] + *selectedJets[MVAvals2ndPass.second[2]]).M();
                bestTopPt =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]).Pt();

                // cout <<"Indices of best MVA jets are :  "<< MVAvals1.second[0] <<"  "<< MVAvals1.second[1]  <<" " << MVAvals1.second[2]<<endl;

                //   cout <<"MVA Mass 1 = "<< bestTopMass1 << " MVA Mass 2 = "<< bestTopMass2ndPass << endl;

                // cout <<"   "<<endl;

                if (debug) cout <<"Processing event with jetcombiner : 9 "<< endl;    






                MSPlot["MVA1TriJetMass"]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MVA2ndPassTriJetMass"]->Fill(bestTopMass2ndPass,  datasets[d], true, Luminosity*scaleFactor );

                if (debug)  cout <<"MVA Mass 1 = "<< bestTopMass1 << " MVA Mass 2 = "<< bestTopMass2 << endl;



                ////////////////////////////////////////////////////////////////////////////////////////////////////
                /////  Calculating how well the MVA jet selection is doing: Fraction of ttbar events            ////
                ////    where the jets selected by the TMVA massReco match the true jets from the hadronic top) ////
                ////////////////////////////////////////////////////////////////////////////////////////////////////

                if(   ( hadronicBJet_.first == MVAvals1.second[0] || hadronicBJet_.first == MVAvals1.second[1] || hadronicBJet_.first == MVAvals1.second[2]   )  && ( hadronicWJet1_.first == MVAvals1.second[0] || hadronicWJet1_.first == MVAvals1.second[1] || hadronicWJet1_.first == MVAvals1.second[2]   )    && ( hadronicWJet2_.first == MVAvals1.second[0] || hadronicWJet2_.first == MVAvals1.second[1] || hadronicWJet2_.first == MVAvals1.second[2]   )      ){
                    nMVASuccesses++;
                }

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
            //double HT1M2L=0, H1M2L=0, HTbjets=0, HT2M=0, H2M=0, HT2L2J=0,nLights;
            float HT1M2L=0, H1M2L=0, HTbjets=0, HT2M=0, H2M=0, HT2L2J=0, nLights;


            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
                MSPlot["BdiscBJetCand_CSV"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(),datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                            //Event-level variables
                jetpt = selectedJets[seljet1]->Pt();
                HT = HT + jetpt;
                H = H +  selectedJets[seljet1]->P();
                if (seljet1 > 2  )  HTHi +=  selectedJets[seljet1]->Pt();
            }

            //HTH = HT/H;
            HTH = HTHi/HT;

            vector<double> ptList;
            double jetptTemp;
            for(Int_t csvJets = 2; csvJets<selectedJets.size(); csvJets++){
                jetptTemp = (double)selectedJets[csvJets]->Pt();
                ptList.push_back(jetptTemp);
                //cout<<csvJets<<"   ptlist "<<ptList[csvJets-2]<<endl;
                //selectedJets2.push_back(selectedJets[csvJets]);  //created array of selected jets without 2 highest CSVL btags
                //cout<<csvJets<<" jet pt "<<selectedJets2[csvJets-2]->Pt()<<"   "<<selectedJets[csvJets]->Pt()<<endl;
            }

            sort(selectedJets2.begin(),selectedJets2.end(),HighestPt()); //order Jets wrt Pt for tuple output
            nLights=selectedLightJets.size();

            HT2L2J = HT - selectedJets[0]->Pt() - selectedJets[1]->Pt() - ptList[0] - ptList[1];    
            //cout<<"HT:  "<<HT<<"  "<<selectedJets[0]->Pt()<<"  "<<selectedJets[1]->Pt()<<"  "<<ptList[0]<<"  "<<ptList[1]<<"  HT2l2J"<<HT2L2J<<endl;    

            //HT - (2 highest CSVL btags) and (2 highest pt jets from remaining jets)

            MSPlot["HTb_SelectedJets"]->Fill(HTb, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["HTRat"]->Fill(HTRat, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["NbOfSelectedBJets"]->Fill(selectedLBJets.size(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);


            MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            //MSPlot["7thJetPt"]->Fill(selectedJets[6]->Pt(), datasets[d], true, Luminosity*scaleFactor);

            MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
            sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order Jets wrt Pt for tuple output

            muonpt = selectedMuons[0]->Pt();
            muoneta = selectedMuons[0]->Eta();
            //electronpt  = selectedMuons[0]->Pt();
            //electroneta  = selectedMuons[0]->Eta();

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
            
            for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it){

                //  cout <<"MVA Method : "<< it->first    <<" Score : "<< it->second <<endl;
                BDTScore = it->second;

            }

                float nvertices = vertex.size();
                float normfactor = datasets[d]->NormFactor();

            float vals[20] = {BDTScore,nJets,nLtags,nMtags,nTtags,HT,muonpt,muoneta,bjetpt,HT2M,HTb,HTH,HTRat,topness,scaleFactor,nvertices,normfactor,Luminosity,weight_0};

            tup->Fill(vals);

            //tup->Fill(HT,nJets,nLtags,nLights,selectedJets[4]->Pt(),selectedJets[5]->Pt(),HTRat,HTb,HT2L2J,MultiTopness,scaleFactor,datasets[d]->NormFactor(),Luminosity);

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
    selecTable.Write(  "FourTop"+postfix+"_Table"+channelpostfix+".tex",    false,true,true,true,false,false,true);

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
    delete fout;

    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

    return 0;
}
