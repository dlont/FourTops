//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.            ////
//////////////////////////////////////////////////////////////////////////////

// ttbar @ NLO 13 TeV:                              //ttbar @ NNLO 8 TeV:
//all-had ->679 * .46 = 312.34                      //all-had -> 245.8 * .46 = 113.068
//semi-lep ->679 *.45 = 305.55                      //semi-lep-> 245.8 * .45 = 110.61
//di-lep-> 679* .09 = 61.113                        //di-lep ->  245.8 * .09 = 22.122
#define _USE_MATH_DEFINES
#include "TStyle.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TNtuple.h"
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
#include "TopTreeAnalysisBase/Tools/interface/BTagCalibrationStandalone.h"

#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
// #include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

#include "TopTreeAnalysisBase/../FourTops/SingleLepAnalysis/interface/CutsTable.h"
#include "TopTreeAnalysisBase/../FourTops/SingleLepAnalysis/interface/HadronicTopReco.h"
#include "TopTreeAnalysisBase/../FourTops/SingleLepAnalysis/interface/EventBDT.h"
#include "TopTreeAnalysisBase/../FourTops/SingleLepAnalysis/interface/Zpeak.h"
#include "TopTreeAnalysisBase/../FourTops/SingleLepAnalysis/interface/Trigger.h"

using namespace std;
using namespace TopTree;
using namespace reweight;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

struct HighestCVSBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_combinedInclusiveSecondaryVertexV2BJetTags() > j2->btag_combinedInclusiveSecondaryVertexV2BJetTags();
    }
};


int main (int argc, char *argv[])
{

    //Checking Passed Arguments to ensure proper execution of MACRO
    if(argc < 12)
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
    const int startEvent            = strtol(argv[argc-2], NULL, 10);
    const int endEvent              = strtol(argv[argc-1], NULL, 10);

    vector<string> vecfileNames;
    cout<<"arg 10: "<<argv[10]<<endl;
    for(int args = 11; args < argc-2; args++) //argc-2 starting at 11
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
    for(int vecfiles=0; vecfiles<vecfileNames.size(); vecfiles++){
        cout<<"vecfile names: "<<vecfiles<<" : "<<vecfileNames[vecfiles]<<endl;
    }
    ofstream eventlist;
    eventlist.open ("interesting_events_mu2.txt");

    int passed = 0;
    int preTrig = 0;
    int postTrig = 0;
    int ndefs =0;
    int negWeights = 0;
    float weightCount = 0.0;
    int eventCount = 0;
    float scalefactorbtageff, mistagfactor;
    string dataSetName = "";
    string channelpostfix = "";
    string postfix = "_Run2_TopTree_Study"; // to relabel the names of the output file

    clock_t start = clock();

    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the FourTop search ! "           << endl;
    cout << "*************************************************************" << endl;


    ///////////////////////////////////////
    //      Configuration                //
    ///////////////////////////////////////

    bool SingleLepton  = true;
    bool Muon          = true;
    bool Electron      = false;
    bool HadTopOn      = true;
    bool EventBDTOn    = true;
    bool TrainMVA      = false; // If false, the previously trained MVA will be used to calculate stuff
    bool bx25          = true;
    bool bTagReweight  = true;
    bool bLeptonSF     = true;
    bool debug         = false;
    bool applyJER      = false;
    bool applyJEC      = false;
    bool JERNom        = true;
    bool JERUp         = false;
    bool JERDown       = false;
    bool JESUp         = false;
    bool JESDown       = false;
    bool fillingbTagHistos = true;
    string MVAmethod   = "BDT"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
    float Luminosity   = 2460.0 ; //pb^-1 shown is C+D, D only is 2094.08809124; silverJson
    //bool split_ttbar = false;

    if(JERDown){
        postfix += "_JER_Down";
    }
    else if(JERUp){
        postfix += "_JER_Up";
    }
    if(JESDown){
        postfix += "_JES_Down";
    }
    else if(JESUp){
        postfix += "_JES_Up";
    }    

    if(Muon && SingleLepton){
        cout<<" ***** USING SINGLE MUON CHANNEL  ******"<<endl;
        channelpostfix = "_Mu";
    }
    else if(Electron && SingleLepton){
        cout<<" ***** Using SINGLE ELECTRON CHANNEL *****"<<endl;
        channelpostfix = "_El";
    }
    else    {
        cerr<<"Correct Channel not selected."<<endl;
        exit(1);
    }

    /////////////////////////////////
    //  Set up AnalysisEnvironment //
    /////////////////////////////////

    AnalysisEnvironment anaEnv;
    cout<<" - Creating environment ..."<<endl;
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

    ///////////////////////////////////////////
    //            Load datasets              //
    ///////////////////////////////////////////

    TTreeLoader treeLoader;
    vector < Dataset* > datasets;    cout << " - Creating Dataset ..." << endl;
    Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
    theDataset->SetEquivalentLuminosity(EqLumi);
    datasets.push_back(theDataset);
    dataSetName = theDataset->Name();

    //////////////////////////////////////////////////////
    //     bTag calibration reader and weight tools     //
    //////////////////////////////////////////////////////

    BTagCalibration * bTagCalib;   
    BTagCalibrationReader * bTagReader;
    BTagCalibrationReader * bTagReaderUp;
    BTagCalibrationReader * bTagReaderDown;
    BTagWeightTools *btwt;
    BTagWeightTools *btwtUp;
    BTagWeightTools *btwtDown;

    if(bTagReweight && dataSetName.find("Data")==string::npos){
        //Btag documentation : http://mon.iihe.ac.be/~smoortga/TopTrees/BTagSF/BTaggingSF_inTopTrees.pdf //v2 or _v2
        bTagCalib = new BTagCalibration("CSVv2","../TopTreeAnalysisBase/Calibrations/BTagging/CSVv2_13TeV_25ns_combToMujets.csv");
        bTagReader = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"mujets","central"); //mujets
        bTagReaderUp = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"mujets","up"); //mujets
        bTagReaderDown = new BTagCalibrationReader(bTagCalib,BTagEntry::OP_MEDIUM,"mujets","down"); //mujets

        if(fillingbTagHistos) {
            if (Muon) {
                btwt = new BTagWeightTools(bTagReader,"HistosPtEta_"+dataSetName+"_Mu.root",false,30,500,2.4);
                btwtUp = new BTagWeightTools(bTagReader,"HistosPtEta_"+dataSetName+"_Mu_Up.root",false,30,500,2.4);
                btwtDown = new BTagWeightTools(bTagReader,"HistosPtEta_"+dataSetName+"_Mu_Down.root",false,30,500,2.4);
            }
            else if (Electron) {
                btwt = new BTagWeightTools(bTagReader,"HistosPtEta_"+dataSetName+"_El.root",false,30,500,2.4);
                btwtUp = new BTagWeightTools(bTagReader,"HistosPtEta_"+dataSetName+"_El_Up.root",false,30,500,2.4);
                btwtDown = new BTagWeightTools(bTagReader,"HistosPtEta_"+dataSetName+"_El_Down.root",false,30,500,2.4);
            }
        }    
        else {
            if(Muon){
                btwt = new BTagWeightTools(bTagReader,"btaghistos/HistosPtEta_TTJets_MLM_Mu_4J_2.root",false,30,500,2.4); 
                btwtUp = new BTagWeightTools(bTagReader,"btaghistos/HistosPtEta_TTJets_MLM_Mu_4J_2_Up.root",false,30,500,2.4); 
                btwtDown = new BTagWeightTools(bTagReader,"btaghistos/HistosPtEta_TTJets_MLM_Mu_4J_2_Down.root",false,30,500,2.4); 
            }
            else if (Electron){
                btwt = new BTagWeightTools(bTagReader,"btaghistos/HistosPtEta_TTJets_MLM_El_4j_2.root",false,30,500,2.4); 
                btwtUp = new BTagWeightTools(bTagReader,"btaghistos/HistosPtEta_TTJets_MLM_El_4j_2_Up.root",false,30,500,2.4); 
                btwtDown = new BTagWeightTools(bTagReader,"btaghistos/HistosPtEta_TTJets_MLM_El_4j_2_Down.root",false,30,500,2.4); 
            }
        }
    }

    /////////////////////////////////////////////////
    //                   Lepton SF                 //
    /////////////////////////////////////////////////
    MuonSFWeight* muonSFWeightID_T;   
    MuonSFWeight *muonSFWeightIso_TT;
    MuonSFWeight* muonSFWeightTrigC_TT;
    MuonSFWeight* muonSFWeightTrigD1_TT;
    MuonSFWeight* muonSFWeightTrigD2_TT;


    ElectronSFWeight* electronSFWeight; 
    if(bLeptonSF){
        if(Muon){ 
            // muonSFWeight = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/Muon_SF_TopEA.root","SF_totErr",false,false);  OLD SF WEIGHT
            muonSFWeightID_T = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonID_Z_RunD_Reco74X_Nov20.root", "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);
            muonSFWeightIso_TT = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/MuonIso_Z_RunD_Reco74X_Nov20.root", "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio", true, false, false);  // Tight RelIso, Tight ID
            muonSFWeightTrigC_TT = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root", "runCreRECO_IsoMu20_OR_IsoTkMu20_PtEtaBins/abseta_pt_ratio", true, false, false);
            muonSFWeightTrigD1_TT = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins", true, false, false);
            muonSFWeightTrigD2_TT = new MuonSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root", "runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins", true, false, false);

        }
        else if(Electron){
            electronSFWeight = new ElectronSFWeight("../TopTreeAnalysisBase/Calibrations/LeptonSF/Elec_SF_TopEA.root","GlobalSF",false,false);    
        }
    }
    /////////////////////////////////////////////////
    ////                Event BDT                 ///
    /////////////////////////////////////////////////

    EventBDT* eventBDT;
    if (EventBDTOn){
        HadTopOn = true;
        eventBDT = new EventBDT();
        eventBDT->initaliseEventComp();
    }
    cout << " Initialized Eventcomputer_ for event_level BDT" << endl;

    /////////////////////////////////////////////////
    ////                 Trigger                  ///
    /////////////////////////////////////////////////
    Trigger* trigger = new Trigger(Muon, Electron);
    trigger->bookTriggers();

    /////////////////////////////////////////////////
    //             Get Luminosity for data         //
    /////////////////////////////////////////////////

    cout <<"found sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    if(dataSetName.find("Data") != string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos)
    {
        Luminosity = theDataset->EquivalentLumi();      cout <<"found DATA sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    }
    cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;

    /////////////////////////////////////////////////
    //               Output ROOT file              //
    /////////////////////////////////////////////////
    string rootFileName ("FourTop"+postfix+"_"+dName+channelpostfix+".root"); //eg. FourTop_Run2_TopTree_Study_Data_Mu.root
    TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");


    /////////////////////////////////////////////////
    //                Top Reco MVA                 //
    /////////////////////////////////////////////////
    HadronicTopReco *hadronicTopReco;
    if(HadTopOn){
        hadronicTopReco = new HadronicTopReco(fout, Muon, Electron, TrainMVA, datasets, MVAmethod, debug, Luminosity);
    }
    /////////////////////////////////////////////////
    //            vectors of objects               //
    /////////////////////////////////////////////////
    cout << " - Variable declaration ..." << endl;
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
    vector < TRootMET* >      mets;
    vector < TRootGenJet* > genjets;

    /////////////////////////////////////////////////
    //              Global variable                //
    /////////////////////////////////////////////////
    TRootEvent* event = 0;
    TRootRun *runInfos = new TRootRun();

    ///////////////////////////////////////////////////////
    //               MultiSample plots                   //
    ///////////////////////////////////////////////////////

    MSPlot["NbOfVertices"]          = new MultiSamplePlot(datasets, "NbOfVertices", 60, 0, 60, "Nb. of vertices");
    //Muons
    MSPlot["MuonPt"]                = new MultiSamplePlot(datasets, "MuonPt", 30, 0, 300, "PT_{#mu}");
    MSPlot["leptonIso"]             = new MultiSamplePlot(datasets, "LeptonIso", 10, 0, 0.25, "RelIso");    
    //Electrons
    MSPlot["ElectronRelIsolation"]  = new MultiSamplePlot(datasets, "ElectronRelIsolation", 10, 0, .25, "RelIso");
    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]     = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV", 20, 0, 1, "CSV b-disc.");
    MSPlot["NbOfSelectedBJets"]     = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of tags");
    MSPlot["HTb_SelectedJets"]      = new MultiSamplePlot(datasets, "HTb_SelectedJets", 50, 0, 1500, "HTb");    
    //Jets
    MSPlot["JetEta"]                = new MultiSamplePlot(datasets, "JetEta", 40,-4, 4, "Jet #eta");
    MSPlot["JetpT"]                 = new MultiSamplePlot(datasets, "JetpT", 100, 0, 400, "Jet p_{T}");
    MSPlot["HT_SelectedJets"]       = new MultiSamplePlot(datasets, "HT_SelectedJets", 30, 0, 1500, "HT");
    MSPlot["HTRat"]                 = new MultiSamplePlot(datasets, "HTRat", 50, 0, 20, "HTRat");
    MSPlot["5thJetPt"]              = new MultiSamplePlot(datasets, "5thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["6thJetPt"]              = new MultiSamplePlot(datasets, "6thJetPt", 60, 0, 400, "PT_{jet}");
    
    MSPlot["MET"]                   = new MultiSamplePlot(datasets, "MET", 70, 0, 700, "MET");
    MSPlot["NbOfBadTrijets"]        = new MultiSamplePlot(datasets, "NbOfBadTriJets", 150, 0, 150, "Nb. of Bad Combs");

    MSPlot["TriJetMass_Matched"]    = new MultiSamplePlot(datasets, "TriJetMassMatched", 100, 0, 1000, "m_{bjj}");
    MSPlot["TriJetMass_UnMatched"]  = new MultiSamplePlot(datasets, "TriJetMassUnMatched", 100, 0, 1000, "m_{bjj}");

    MSPlot["MVA2TriJetMass"]        = new MultiSamplePlot(datasets, "MVA2TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1TriJetMassMatched"] = new MultiSamplePlot(datasets, "MVA1TriJetMassMatched", 75, 0, 500, "m_{bjj}");

    /////////////////////////////////////////////////
    //                  Plots path                 //
    /////////////////////////////////////////////////
    string pathPNG = "FourTop"+postfix+channelpostfix+"_MSPlots/";
    mkdir(pathPNG.c_str(),0777);    cout <<"Making directory :"<< pathPNG  <<endl;

    /////////////////////////////////////////////////
    //                 Cuts table                  //
    /////////////////////////////////////////////////
    CutsTable *cutsTable = new CutsTable(Muon, Electron);
    cutsTable->AddSelections();
    cutsTable->CreateTable(datasets, Luminosity);

    /////////////////////////////////////////////////
    //                Z peak maker                 //
    /////////////////////////////////////////////////
    Zpeak *zPeakMaker;
    zPeakMaker = new Zpeak(datasets);

    /////////////////////////////////////////////////
    //               Pu reweighting                //
    /////////////////////////////////////////////////
    LumiReWeighting LumiWeights;
    LumiWeights = LumiReWeighting("/user/lbeck/CMSSW_7_6_0/src/TopBrussels/TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_MC_RunIISpring15DR74-Asympt25ns.root", "/user/lbeck/CMSSW_7_6_0/src/TopBrussels/TopTreeAnalysisBase/Calibrations/PileUpReweighting/pileup_2015Data74X_25ns-Run246908-260627Cert_Silver.root", "pileup50", "pileup");    

    ///////////////////////////////////////////
    ///  Initialise Jet Energy Corrections  ///
    ///////////////////////////////////////////
    
    vector<JetCorrectorParameters> vCorrParam;
    string pathCalJEC = "../TopTreeAnalysisBase/Calibrations/JECFiles/";

    if(dataSetName.find("Data")!=string::npos)
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
    }
    else
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters(pathCalJEC+"Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(pathCalJEC+"Summer15_25nsV6_MC_Uncertainty_AK4PFchs.txt");
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true); //true means redo also L1
    /////////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                             //
    //                                      Loop on datasets                                       //
    //                                                                                             //
    /////////////////////////////////////////////////////////////////////////////////////////////////
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;    
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        string previousFilename2 = "";
        string currentfilename2 = "";

        /////////////////////////////////////////////////
        //                 nlo or bx25?                //
        /////////////////////////////////////////////////
        bool nlo = true;
        dataSetName = datasets[d]->Name();
        if(dataSetName.find("bx50") != std::string::npos) bx25 = false;
        else bx25 = true;
        if(bx25) cout << "Dataset with 25ns Bunch Spacing!" <<endl;
        else cout << "Dataset with 50ns Bunch Spacing!" <<endl;

        if(dataSetName.find("NLO") != std::string::npos || dataSetName.find("nlo") !=std::string::npos) nlo = true;
        else nlo = false;
        if(nlo) cout << "NLO Dataset!" <<endl;
        else cout << "LO Dataset!" << endl;

        ///////////////////////////////////////////////////////
        //      Setup Date string and nTuple for output      //
        ///////////////////////////////////////////////////////

        SourceDate *strdate = new SourceDate();
        string date_str = strdate->ReturnDateStr();
        if(debug)cout<<"date print"<<endl;

        /////////////////////////////////////////////////
        //               Craneen setup                 //
        /////////////////////////////////////////////////
        string channel_dir = "Craneens"+channelpostfix;
        string date_dir = channel_dir+"/Craneens" + date_str +"/";
        int mkdirstatus = mkdir(channel_dir.c_str(),0777);
        mkdirstatus = mkdir(date_dir.c_str(),0777);
        if(debug)cout<<"created dirs"<<endl;
        string Ntuptitle   = "Craneen_" + channelpostfix;
        
        string Ntupname    = "Craneens" + channelpostfix + "/Craneens" + date_str + "/Craneen_" + dataSetName + postfix + ".root";     
        TFile * tupfile    = new TFile(Ntupname.c_str(),"RECREATE");
        TNtuple * tup      = new TNtuple(Ntuptitle.c_str(), Ntuptitle.c_str(), "BDT:nJets:NOrigJets:nLtags:nMtags:nTtags:HT:LeptonPt:LeptonEta:LeadingBJetPt:HT2M:HTb:HTH:HTRat:HTX:SumJetMassX:multitopness:nbb:ncc:nll:ttbar_flav:ScaleFactor:SFlepton:SFbtag:SFbtagUp:SFbtagDown§:SFPU:PU:NormFactor:Luminosity:GenWeight:weight1:weight2:weight3:weight4:weight5:weight6:weight7:weight8:met:angletop1top2:angletoplep:1stjetpt:2ndjetpt:leptonIso:leptonphi:chargedHIso:neutralHIso:photonIso:PUIso:5thjetpt:6thjetpt:jet5and6pt");
       
        // string Ntup4j0bname    = "Craneens" + channelpostfix + "/Craneens" + date_str + "/Craneen_4j0b_" + dataSetName + postfix + ".root";     
        // TFile * tup4j0bfile    = new TFile(Ntupname.c_str(),"RECREATE");
        // TNtuple * tup4j0b      = new TNtuple(Ntuptitle.c_str(), Ntuptitle.c_str(), "BDT:nJets:NOrigJets:nLtags:nMtags:nTtags:HT:LeptonPt:LeptonEta:LeadingBJetPt:HT2M:HTb:HTH:HTRat:multitopness:nbb:ncc:nll:ttbar_flav:ScaleFactor:SFlepton:SFbtag:SFPU:PU:NormFactor:Luminosity:GenWeight:weight1:weight2:weight3:weight4:weight5:weight6:weight7:weight8:met:angletop1top2:angletoplep:1stjetpt:2ndjetpt:leptonIso:leptonphi:chargedHIso:neutralHIso:photonIso:PUIso");
        

        string Ntupjetname = "Craneens" + channelpostfix + "/Craneens" + date_str + "/CraneenJets_" + dataSetName + postfix + ".root";
        TFile * tupjetfile = new TFile(Ntupjetname.c_str(),"RECREATE");
        TNtuple * tupjet   = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(), "jetpT:csvDisc:jeteta:jetphi:jetLeptDR:ScaleFactor:NormFactor:Luminosity");
        
        string NtupZname   = "Craneens" + channelpostfix + "/Craneens" + date_str + "/CraneenZ_" + dataSetName + postfix + ".root";
        TFile * tupZfile   = new TFile(NtupZname.c_str(),"RECREATE");
        TNtuple * tupZ     = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(), "invMassll:ScaleFactor:NormFactor:Luminosity");
        if(debug)cout<<"created craneens"<<endl;


        ////////////////////////////////////////////////////////////
        //       Define object containers and initalisations      //
        ////////////////////////////////////////////////////////////

        float BDTScore, MHT, MHTSig, STJet,leptoneta, leptonpt, leptonphi, electronpt, electroneta, bjetpt, EventMass, EventMassX, SumJetMass, SumJetMassX, H, HX;
        float HTHi, HTRat, HT, HTX, HTH, HTXHX, sumpx_X, sumpy_X, sumpz_X, sume_X, sumpx, sumpy, sumpz, sume, jetpt, PTBalTopEventX, PTBalTopSumJetX, PTBalTopMuMet;     
        int itrigger = -1, previousRun = -1;
        int currentRun;           
        vector<TRootElectron*> selectedElectrons;
        vector<TRootPFJet*>    selectedOrigJets; //all original jets before jet lepton cleaning
        vector<TRootPFJet*>    selectedJets; //all jets after jet lepton cleaning
        vector<TRootPFJet*>    selectedJets2; //after removal of 2 highest CSVL btags
        vector<TRootMuon*>     selectedMuons;
        vector<TRootElectron*> selectedExtraElectrons;
        vector<TRootElectron*> selectedOrigElectrons;
        vector<TRootElectron*> selectedOrigExtraElectrons;

        vector<TRootMuon*>     selectedExtraMuons;
        selectedElectrons.reserve(10);
        selectedMuons.reserve(10);
        vector<TRootPFJet*>      selectedJets2ndPass;
        vector<TRootPFJet*>      selectedJets3rdPass;
        vector<TRootPFJet*>      MVASelJets1;
        vector<TRootJet*>      selectedLBJets; //CSVL btags
        vector<TRootJet*>      selectedMBJets; //CSVM btags
        vector<TRootJet*>      selectedTBJets; //CSVT btags
        vector<TRootJet*>      selectedLightJets;
        int iFile2 = -1;
        datasets[d]->runTree()->SetBranchStatus("runInfos*",1);
        datasets[d]->runTree()->SetBranchAddress("runInfos",&runInfos);
        if (dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos || dataSetName.find("DATA")!=string::npos) TrainMVA=false;

        ///////////////////////////////////////////////////////////
        //             Get # of events to run over               //
        ///////////////////////////////////////////////////////////
//
        int start = 0;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();    cout <<"Number of events = "<<  ending  <<endl;
        int event_start = startEvent;
        double end_d = ending;
        if(dataSetName.find("Data")!=string::npos || dataSetName.find("data")!=string::npos  || dataSetName.find("DATA")!=string::npos ){
            // end_d=ending;
            if(endEvent > ending)            end_d = ending;
            else            end_d = endEvent;  
        }
        else{  //this only works if using parallel launcher to split into jobs
            if(endEvent > ending)            end_d = ending;
            else            end_d = endEvent;            
        }

        cout <<"Will run over "<<  end_d<< " events..."<<endl;    cout <<"Starting event = = = = "<< event_start  << endl;

        ////////////////////////////////////////////////////////////////////////////////
        //                                 Loop on events                             //
        ////////////////////////////////////////////////////////////////////////////////

        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
            if(debug) cout<<"START OF EVENT LOOP"<<endl;
            BDTScore= -99999.0, MHT = 0., MHTSig = 0.,leptoneta = 0., leptonpt =0., electronpt=0., electroneta=0., bjetpt =0., STJet = 0.;
            EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0., HTHi =0., HTRat = 0;  H = 0., HX =0., HT = 0., HTX = 0.;
            HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0.;
            PTBalTopEventX = 0., PTBalTopSumJetX =0.;

            double ievt_d = ievt;

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC 
                << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            bool trigged = false;  // Disabling the HLT requirement
            if(debug) cout<<"before tree load"<<endl;
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);      if(debug)cout<<"after tree load"<<endl; //load event
            if (dataSetName.find("Data")==string::npos) {
                genjets = treeLoader.LoadGenJet(ievt,false);
            //sort(genjets.begin(),genjets.end(),HighestPt()); // HighestPt() is included from the Selection class
            }
            float rho = event->fixedGridRhoFastjetAll();        if (debug)cout <<"Rho: " << rho <<endl;


            ///////////////////////////////////////////
            //      Set up for miniAOD weights       //
            ///////////////////////////////////////////

            currentRun = event->runId(); if(debug) cout<<"got run ID"<<endl;
            datasets[d]->eventTree()->LoadTree(ievt); if(debug) cout<<"load tree"<<endl;
            int treenumber = datasets[d]->eventTree()->GetTreeNumber(); 
            // cout<<"treenumber: "<< treenumber<<endl;
            currentfilename2 = datasets[d]->eventTree()->GetFile()->GetName();
            if(previousFilename2 != currentfilename2){
                previousFilename2 = currentfilename2;
                iFile2++;
                cout<<"got tree number: "<<treenumber<<endl;
                cout<<"File changed!!!! => iFile2 = "<<iFile2 << " new file is " << datasets[d]->eventTree()->GetFile()->GetName() << " in sample " << datasets[d]->Name() << endl;
            }

            //int rBytes = datasets[d]->runTree()->GetEntry(treenumber);  if(debug) cout<<"get entry with treenumber"<<endl;

            ///////////////////////////////////////////////////////////
            //           Object definitions for selection            //
            ///////////////////////////////////////////////////////////

            Run2Selection r2selection(init_jets, init_muons, init_electrons, mets);
            int nMu = 0, nEl = 0, nLooseMu = 0, nLooseEl = 0; //number of (loose) muons/electrons

            if(Electron){
                                                                                                                                           
                selectedOrigJets                                    = r2selection.GetSelectedJets(); if (debug)cout<<"Getting Jets"<<endl; // ApplyJetId
                                                                                                                                            
                selectedMuons                                       = r2selection.GetSelectedMuons(10, 2.5, 0.25, "Loose", "Spring15"); if (debug)cout<<"Getting Loose Muons"<<endl;
                nMu = selectedMuons.size(); //Number of Muons in Event
                                                                                                                                            
                selectedOrigElectrons                                   = r2selection.GetSelectedElectrons(30, 2.1, "Tight", "Spring15_25ns", true); if (debug)cout<<"Getting Tight Electrons"<<endl; // VBTF ID       
                                                                                                                                            
                selectedOrigExtraElectrons                              = r2selection.GetSelectedElectrons(15, 2.5, "Veto", "Spring15_25ns", true); if (debug)cout<<"Getting Loose Electrons"<<endl;

            }
            else if(Muon){
                                                                                                                                           
                selectedOrigJets                                    = r2selection.GetSelectedJets();  if (debug)cout<<"Getting Jets"<<endl; // ApplyJetId
                                                                                                                                            
                selectedMuons                                       = r2selection.GetSelectedMuons(26, 2.1, 0.15, "Tight", "Spring15"); if (debug)cout<<"Getting Tight Muons"<<endl;
                nMu = selectedMuons.size(); //Number of Muons in Event
                                                                                                                                            
                selectedOrigElectrons                                   = r2selection.GetSelectedElectrons(15, 2.5, "Veto", "Spring15_25ns", true); if (debug)cout<<"Getting Loose Electrons"<<endl; // VBTF ID    
                                                                                                                                            
                selectedExtraMuons                                  = r2selection.GetSelectedMuons(10, 2.5, 0.25, "Loose", "Spring15"); if (debug)cout<<"Getting Loose Muons"<<endl;
                nLooseMu = selectedExtraMuons.size();   //Number of loose muons      
            }
            //if(nEl>0) cout<<"nEl: "<<nEl<<endl;

            //remove electrons between 1.4442 and 1.5660

            for(int e_iter=0; e_iter<selectedOrigElectrons.size();e_iter++){
                if(selectedOrigElectrons[e_iter]->Eta()<=1.4442 || selectedOrigElectrons[e_iter]->Eta()>=1.5660){
                    selectedElectrons.push_back(selectedOrigElectrons[e_iter]);
                }
            }
            nEl = selectedElectrons.size(); //Number of Electrons in Event   

            if(Electron){
                for(int e_iter=0; e_iter<selectedOrigExtraElectrons.size();e_iter++){
                    if(selectedOrigExtraElectrons[e_iter]->Eta()<=1.4442 || selectedOrigExtraElectrons[e_iter]->Eta()>=1.5660){
                        selectedExtraElectrons.push_back(selectedOrigExtraElectrons[e_iter]);
                    }
                }
                nLooseEl = selectedExtraElectrons.size(); //Number of loose muons
             cout<<"nel: "<<nEl<<"  nLooseEl: "<<nLooseEl<<"  origel: "<<selectedOrigElectrons.size()<<"  origextra el: "<<selectedOrigExtraElectrons.size()<<endl;
            }

            ////////////////////////////////////////////////////////////
            //   Z peak before nMu==1 or nEl==1 for jetlep cleaning   //
            ////////////////////////////////////////////////////////////
            zPeakMaker->invariantMass(r2selection);
            float invMassll = zPeakMaker->returnInvMass();

            /////////////////////////////////////////////////
            //            Jet lepton cleaning              //
            /////////////////////////////////////////////////
            selectedJets.clear();
            //cout<<nMu<<"<--nmu  nEl-->"<<nEl<<endl;
            if(Muon && nMu>0){
                for (int origJets=0; origJets<selectedOrigJets.size(); origJets++){
                    //cout<<"DR: "<< selectedOrigJets[origJets]->DeltaR(*selectedMuons[0])<<endl;
                    if(selectedOrigJets[origJets]->DeltaR(*selectedMuons[0])>0.4){
                        selectedJets.push_back(selectedOrigJets[origJets]);
                    }                    
                }
            }
            else if(Electron && nEl>0){
                for (int origJets=0; origJets<selectedOrigJets.size(); origJets++){
                    //cout<<"DR: "<< selectedOrigJets[origJets]->DeltaR(*selectedMuons[0])<<endl;
                    if(selectedOrigJets[origJets]->DeltaR(*selectedElectrons[0])>0.4){
                        selectedJets.push_back(selectedOrigJets[origJets]);
                    }                       
                }
            }
            else selectedJets = selectedOrigJets;

            ///////////////////////////////////////////////////////////////////////////////////
            // Preselection looping over Jet Collection                                      //
            // Summing HT and calculating leading, lagging, and ratio for Selected and BJets //
            ///////////////////////////////////////////////////////////////////////////////////

            selectedLBJets.clear();
            selectedMBJets.clear();
            selectedTBJets.clear();
            selectedLightJets.clear();

            float HTb = 0.;  //calculate assigning loose, medium and tight tags
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
                if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.605   )
                {
                    selectedLBJets.push_back(selectedJets[seljet]);
                    if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.890)
                    {
                        HTb += selectedJets[seljet]->Pt();
                        selectedMBJets.push_back(selectedJets[seljet]);
                        if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.970)
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
            float nLtags = selectedLBJets.size(); //Number of CSVL tags in Event (includes jets that pass CSVM & CSVT)
            float nMtags = selectedMBJets.size(); //Number of CSVM tags in Event (includes jets that pass CSVT)
            float nTtags = selectedTBJets.size(); //Number of CSVT tags in Event 
            float nLights=selectedLightJets.size();


            //htrat leading lagging calculation
            float temp_HT = 0.;
            float HT_leading = 0.;
            float HT_lagging = 0.;
            float HTRat = 0;
            for (Int_t seljet0 =0; seljet0 < selectedJets.size(); seljet0++ ){
                temp_HT += selectedJets[seljet0]->Pt();
                if (seljet0 < 4){
                    HT_leading += selectedJets[seljet0]->Pt();
                }else{
                    HT_lagging += selectedJets[seljet0]->Pt();
                }
            }

            //HTRat = HT_leading/HT_lagging;

            ///////////////////////////////////////////
            //     Apply primary vertex selection    //
            ///////////////////////////////////////////

            bool isGoodPV = r2selection.isPVSelected(vertex, 4, 24., 2);
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV << " TriggerBit: " << trigged <<endl;
            if (debug) cin.get();


            /////////////////////////////////
            //        Trigger              //
            /////////////////////////////////

            // cout<<"!!CHECK AVAIL!!"<<endl;
            trigger->checkAvail(currentRun, datasets, d, &treeLoader, event, treenumber);
            // cout<<"!!CHECK FIRED!!"<<endl;

            trigged = trigger->checkIfFired(currentRun, datasets, d);
            // trigged=true;
            cutsTable->FillTable(d, isGoodPV, trigged, scaleFactor, nMu, nLooseMu, nEl, nLooseEl, nJets, nLtags, nMtags, nTtags);   if(debug) cout<<"cuts table filled"<<endl;
 
            /////////////////////////////////
            //       Primary vertex        //
            /////////////////////////////////
            //Filling Histogram of the number of vertices before Event Selection
            MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
            if (!isGoodPV) continue; // Check that there is a good Primary Vertex

            /////////////////////////////////
            //        Trigger              //
            /////////////////////////////////            
            preTrig++;
            if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;
            //if(dataSetName.find("Data") != string::npos || dataSetName.find("data") != string::npos || dataSetName.find("DATA") != string::npos){
            if (!trigged)          continue;  //If an HLT condition is not present, skip this event in the loop.       
            //}
            postTrig++; 


            //////////////////////////////////////
            ///  Jet Energy Scale Corrections  ///
            //////////////////////////////////////

            if (applyJER && dataSetName.find("Data")==string::npos)
            {
                if(JERNom) jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal", false);
                else if(JERDown) jetTools->correctJetJER(init_jets, genjets, mets[0], "minus", false);
                else if (JERUp) jetTools->correctJetJER(init_jets, genjets, mets[0], "plus", false);

                /// Example how to apply JES systematics
                if(JESDown) jetTools->correctJetJESUnc(init_jets, "minus", 1);
                else if(JESUp) jetTools->correctJetJESUnc(init_jets, "plus", 1);
                //cout << "JER smeared!!! " << endl;
            }

            bool isData = false;
            if(dataSetName.find("Data")==string::npos){
                isData = true;
            }
            if (applyJEC)   ///should this have  && dataSetName.find("Data")==string::npos
            {
                jetTools->correctJets(init_jets, event->fixedGridRhoFastjetAll(), isData);
            }



            /////////////////////////////////////////////////
            //               Pu reweighting                //
            /////////////////////////////////////////////////

            float lumiWeight;
            if(dataSetName.find("Data") !=string::npos || dataSetName.find("data") != string::npos || dataSetName.find("DATA") != string::npos)
            {
                lumiWeight=1;
            }
            else{
                // lumiWeight = LumiWeights.ITweight( vertex.size() ); 
                lumiWeight = LumiWeights.ITweight( (int)event->nTruePU()); 
                
            }
            // if(lumiWeight<0.2)     cout<<"PU:  "<<(int)event->nTruePU()<<"    LUMI WEIGHT   :   "<<lumiWeight<<" ! "<<endl;
            scaleFactor = scaleFactor * lumiWeight;


            /////////////////////////////////////////////////
            //                    bTag SF                  //
            /////////////////////////////////////////////////

            float bTagEff(-1);
            float bTagEffUp(-1);
            float bTagEffDown(-1);
            if(fillingbTagHistos){
                if(bTagReweight && dataSetName.find("Data")==string::npos){
                //get btag weight info
                    for(int jetbtag = 0; jetbtag<selectedJets.size(); jetbtag++){
                        float jetpt = selectedJets[jetbtag]->Pt();
                        float jeteta = selectedJets[jetbtag]->Eta();
                        float jetdisc = selectedJets[jetbtag]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                        BTagEntry::JetFlavor jflav;
                        int jetpartonflav = std::abs(selectedJets[jetbtag]->partonFlavour());
                        if(debug) cout<<"parton flavour: "<<jetpartonflav<<"  jet eta: "<<jeteta<<" jet pt: "<<jetpt<<"  jet disc: "<<jetdisc<<endl;
                        if(jetpartonflav == 5){
                            jflav = BTagEntry::FLAV_B;
                        }
                        else if(jetpartonflav == 4){
                            jflav = BTagEntry::FLAV_C;
                        }
                        else{
                            jflav = BTagEntry::FLAV_UDSG;
                        }
                        bTagEff = bTagReader->eval(jflav, jeteta, jetpt, jetdisc);    
                        bTagEffUp = bTagReaderUp->eval(jflav, jeteta, jetpt, jetdisc);                     
                        bTagEffDown = bTagReaderDown->eval(jflav, jeteta, jetpt, jetdisc);                     
                 
                        if(debug)cout<<"btag efficiency = "<<bTagEff<<endl;       
                    }      
                    btwt->FillMCEfficiencyHistos(selectedJets); 
                    btwtUp->FillMCEfficiencyHistos(selectedJets); 
                    btwtDown->FillMCEfficiencyHistos(selectedJets); 

                }
            }
            if (debug) cout<<"getMCEventWeight for btag"<<endl;
            float btagWeight = 1;
            float btagWeightUp = 1;
            float btagWeightDown = 1;
            if(bTagReweight && dataSetName.find("Data")==string::npos){
                if(!fillingbTagHistos){
                    btagWeight =  btwt->getMCEventWeight(selectedJets, false);
                    btagWeightUp =  btwtUp->getMCEventWeight(selectedJets, false);
                    btagWeightDown =  btwtDown->getMCEventWeight(selectedJets, false);
                }
                
                if(debug) cout<<"btag weight "<<btagWeight<<"  btag weight Up "<<btagWeightUp<<"   btag weight Down "<<btagWeightDown<<endl;
            }
            scaleFactor *= btagWeight;
            if(ievt<1000)
            {
                cout<<"btagweight: "<<btagWeight<<endl;
            }

            /////////////////////////////////////////////////
            //                   Lepton SF                 //
            /////////////////////////////////////////////////

            float fleptonSF = 1;
            if(bLeptonSF){
                if(Muon && nMu>0){
                    fleptonSF = muonSFWeightID_T->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0) * muonSFWeightIso_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                }
                else if(Electron && nEl>0){
                    fleptonSF = electronSFWeight->at(selectedElectrons[0]->Eta(),selectedElectrons[0]->Pt(),0);
                }
            }

            float trigSFC = 1;
            float trigSFD1 = 1;
            float trigSFD2 = 1;
            float trigSFTot = 1;
            if(bLeptonSF){
                if(dataSetName.find("Data")==string::npos){
                    trigSFC = muonSFWeightTrigC_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                    trigSFD1 = muonSFWeightTrigD1_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);
                    trigSFD2 = muonSFWeightTrigD2_TT->at(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), 0);       
                    trigSFTot =( (trigSFC*17.2) + (trigSFD1*923.88) + (trigSFD2*1639.20) )/Luminosity;  
                }
                fleptonSF*=trigSFTot;
            }

            if(debug) cout<<"lepton SF:  "<<fleptonSF<<endl;
            if(dataSetName.find("Data")==string::npos)   scaleFactor *= fleptonSF;


            /////////////////////////////////
            //        Gen weights          //
            /////////////////////////////////

            float weight_0 = 1; //nominal
            float weight_1 = 1, weight_2 = 1, weight_3 = 1, weight_4 = 1, weight_5 = 1, weight_6 = 1, weight_7 = 1, weight_8 = 1;

            // if(dataSetName.find("tttt") != string::npos){
            //     weight_0 = ( event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("scale_variation 1"))/(abs(event->originalXWGTUP()))); //nominal                   
            // }
            // else if(dataSetName.find("TTJets") != string::npos){
            //     weight_0 = ( event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 1"))/(abs(event->originalXWGTUP()))); //nominal               
            // }
            if(dataSetName.find("Data")==string::npos){
                if(event->getWeight(1)!= -9999){
                    weight_0 = (event->getWeight(1))/(abs(event->originalXWGTUP()));  
                    weight_1 = (event->getWeight(2))/(abs(event->originalXWGTUP()));                
                    weight_2 = (event->getWeight(3))/(abs(event->originalXWGTUP()));                
                    weight_3 = (event->getWeight(4))/(abs(event->originalXWGTUP()));                
                    weight_4 = (event->getWeight(5))/(abs(event->originalXWGTUP()));                
                    weight_5 = (event->getWeight(6))/(abs(event->originalXWGTUP()));                
                    weight_6 = (event->getWeight(7))/(abs(event->originalXWGTUP()));                
                    weight_7 = (event->getWeight(8))/(abs(event->originalXWGTUP()));                
                    weight_8 = (event->getWeight(9))/(abs(event->originalXWGTUP()));                    
                }
                else if (event->getWeight(1001)!= -9999){
                    weight_0 = (event->getWeight(1001))/(abs(event->originalXWGTUP()));  
                    weight_1 = (event->getWeight(1002))/(abs(event->originalXWGTUP()));                
                    weight_2 = (event->getWeight(1003))/(abs(event->originalXWGTUP()));                
                    weight_3 = (event->getWeight(1004))/(abs(event->originalXWGTUP()));                
                    weight_4 = (event->getWeight(1005))/(abs(event->originalXWGTUP()));                
                    weight_5 = (event->getWeight(1006))/(abs(event->originalXWGTUP()));                
                    weight_6 = (event->getWeight(1007))/(abs(event->originalXWGTUP()));                
                    weight_7 = (event->getWeight(1008))/(abs(event->originalXWGTUP()));                
                    weight_8 = (event->getWeight(1009))/(abs(event->originalXWGTUP()));                    
                }
                else {
                    cout<<endl;
                    cout<<"no weights found"<<endl;
                    cout<<endl;
                }
            }
            // if(dataSetName.find("tttt")!=string::npos){
            //     weight_0 = (event->getWeight(1))/ ( abs(event->originalXWGTUP()) );
            // }
            // if(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 1")>=0){
            //     weight_0 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 1")))/(abs(event->originalXWGTUP())); 
            //     cout<<"central"<<endl;  
            // }
            // else if(runInfos->getWeightInfo(currentRun).weightIndex("scale_variation 1")>=0){
            //     weight_0 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("scale_variation 1")))/(abs(event->originalXWGTUP())); 
            //                     cout<<"scale"<<endl;  
 
            // }

            // // cout<<"weight: "<<weight_0<<endl;

            // if(debug) cout<<"weight0 obtained"<<endl;

            // if(dataSetName.find("TTJets")!=string::npos){
            //     cout<<"get weight info for event "<<ievt<<endl;
            //     // weight_1 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 2")))/(abs(event->originalXWGTUP()));                
            //     // weight_2 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 3")))/(abs(event->originalXWGTUP()));                
            //     // weight_3 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 4")))/(abs(event->originalXWGTUP()));                
            //     // weight_4 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 5")))/(abs(event->originalXWGTUP()));                
            //     // weight_5 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 6")))/(abs(event->originalXWGTUP()));                
            //     // weight_6 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 7")))/(abs(event->originalXWGTUP()));                
            //     // weight_7 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 8")))/(abs(event->originalXWGTUP()));                
            //     // weight_8 = (event->getWeight(runInfos->getWeightInfo(currentRun).weightIndex("Central scale variation 9")))/(abs(event->originalXWGTUP()));   
            //     weight_0 = (event->getWeight(1))/(abs(event->originalXWGTUP()));  
            //     weight_1 = (event->getWeight(2))/(abs(event->originalXWGTUP()));                
            //     weight_2 = (event->getWeight(3))/(abs(event->originalXWGTUP()));                
            //     weight_3 = (event->getWeight(4))/(abs(event->originalXWGTUP()));                
            //     weight_4 = (event->getWeight(5))/(abs(event->originalXWGTUP()));                
            //     weight_5 = (event->getWeight(6))/(abs(event->originalXWGTUP()));                
            //     weight_6 = (event->getWeight(7))/(abs(event->originalXWGTUP()));                
            //     weight_7 = (event->getWeight(8))/(abs(event->originalXWGTUP()));                
            //     weight_8 = (event->getWeight(9))/(abs(event->originalXWGTUP()));                    
            // }
            // cout<<"weight0: "<<weight_0<<" weight1: "<<weight_1<<" weight2: "<<weight_2<<" weight3: "<<weight_3<<" weight4: "<<weight_4<<" weight5: "<<weight_5<<" weight6: "<<weight_6<<" weight7: "<<weight_7<<" weight8: "<<weight_8<<endl;
            // cout<<"event->getWeight(1001) :"<<event->getWeight(1)<<endl;
            // cout<<"abs(event->originalXWGTUP()) : "<<abs(event->originalXWGTUP())<<endl;


            /////////////////////////////////////////////////
            //            neg weights counter              //
            /////////////////////////////////////////////////


            if(weight_0 < 0.0)
            {
                //scaleFactor *= -1.0;  //Taking into account negative weights in NLO Monte Carlo
                negWeights++;
            }



            /////////////////////////////////////
            //    Fill cuts table + z peak     //
            /////////////////////////////////////
            zPeakMaker->fillPlot(datasets, d, Luminosity, scaleFactor);

            /////////////////////////////////////////////////
            //                Z peak maker                 //
            /////////////////////////////////////////////////

            float normfactor = datasets[d]->NormFactor();
            float vals2[4] = {invMassll,scaleFactor,normfactor,Luminosity};
            bool isTwoLeptons=zPeakMaker->requireTwoLeptons();     if(debug) cout<<"isTwoLeptons  "<<isTwoLeptons<<endl;

            if (isTwoLeptons){
                tupZfile->cd();        
                tupZ->Fill(vals2);
            }
            if(debug) cout<<"Z peak filled"<<endl;

            //Apply the lepton, jet, btag and HT & MET selections

            if (debug)  cout<<"Number of Muons = "<< nMu <<"    electrons =  "  <<nEl<<"     Jets = "<< selectedJets.size()   <<" loose BJets = "<<  nLtags   <<
                "  MuonChannel = "<<Muon<<" Electron Channel"<<Electron<<endl;

            ///////////////////////////////////////////////////////////
            //               Baseline Event selection                //
            ///////////////////////////////////////////////////////////

            if (Muon)
            {   
                if  (  (!( nMu == 1 /*&& nEl == 0 */&& nLooseMu == 1 && nJets>=6 && nMtags >=2)) )continue; // Muon Channel Selection
            }
            else if(Electron){
                if  (  !( nMu == 0 && nEl == 1 && nLooseEl == 1 && nJets>=6 && nMtags >=2)) continue; // Electron Channel Selection
            }
            else{
                cerr<<"Correct Channel not selected."<<endl;
                exit(1);
            }
            if(debug) cout<<"after baseline"<<endl;
            weightCount += scaleFactor;
            eventCount++;
            if(debug)
            {
                cout<<"Selection Passed."<<endl;
                cin.get();
            }
            passed++;
            /////////////////////////////////////////////////
            //            ttbb reweighting                 //
            /////////////////////////////////////////////////
            float numOfbb = 0;
            float numOfcc = 0;
            float numOfll = 0;
            float ttbar_flav = -1;
            vector<TRootMCParticle*> mcParticles_flav;
            TRootGenEvent* genEvt_flav = 0;
            if(dataSetName.find("TTJets")!=string::npos){
                genEvt_flav = treeLoader.LoadGenEvent(ievt,false);
                treeLoader.LoadMCEvent(ievt, genEvt_flav, 0, mcParticles_flav,false);
                for(unsigned int p=0; p<mcParticles_flav.size(); p++) {
                    //cout<<"status: "<<mcParticles_flav[p]->status()<<"  id: "<<mcParticles_flav[p]->type()<<" mother: "<<mcParticles_flav[p]->motherType()<<endl;
                    if(mcParticles_flav[p]->status()<30 && mcParticles_flav[p]->status()>20 && abs(mcParticles_flav[p]->motherType())!=6){

                        if (abs(mcParticles_flav[p]->type())==5)
                        {
                            // ttbar_flav=2;
                            numOfbb++;  
                        }
                        
                        else if (abs(mcParticles_flav[p]->type())==4 && abs(mcParticles_flav[p]->motherType())!=5 && abs(mcParticles_flav[p]->motherType())!=24)
                        {
                            // ttbar_flav=1;
                            numOfcc++; 
                        }
                        
                        else if (abs(mcParticles_flav[p]->type())<4){
                            // ttbar_flav=1;
                            numOfll++; 
                        }
                    }

                }
            } 

            if(numOfbb>=2){
                ttbar_flav = 2;
            }
            else if(numOfcc>=2){
                ttbar_flav = 1;
            }
            else{
                ttbar_flav = 0;
            }


            //////////////////////////////////////////
            //     TMVA for mass Reconstruction     //
            //////////////////////////////////////////
            float diTopness = 0;
            if (debug) cout<<"TMVA mass reco"<<endl;
            sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());

            if (HadTopOn){
                hadronicTopReco->SetCollections(selectedJets, selectedMuons, selectedElectrons, scaleFactor);
            }

            if(HadTopOn){
                if(!TrainMVA){ //if not training, but computing 
                    hadronicTopReco->Compute1st(d, selectedJets, datasets);
                    hadronicTopReco->Compute2nd(d, selectedJets, datasets);
                    diTopness = hadronicTopReco->ReturnDiTopness();
                    SumJetMassX = hadronicTopReco->ReturnSumJetMassX();
                    HTX = hadronicTopReco->ReturnHTX();// cout<<"HTX: "<<HTX<<endl;
                }
                hadronicTopReco->FillDiagnosticPlots(fout, d, selectedJets, datasets);
            }

            //cout<<"SumJetMassX: "<<SumJetMassX<<endl;
            ///////////////////////////////////
            // Filling histograms / plotting //
            ///////////////////////////////////
            if (debug) cout<<"Plots"<<endl;

            //////////////////////
            // Muon Based Plots //
            //////////////////////
            float leptonIso=0;
            float chargedHIso = 0;
            float neutralHIso = 0;
            float photonIso = 0;
            float PUIso = 0;
            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ )
            {
                float relisomu = (selectedMuons[0]->chargedHadronIso(4) + max( 0.0, selectedMuons[0]->neutralHadronIso(4) + selectedMuons[0]->photonIso(4) - 0.) ) / selectedMuons[0]->Pt();
                chargedHIso = selectedMuons[0]->chargedHadronIso(4);
                neutralHIso = selectedMuons[0]->neutralHadronIso(4);
                photonIso = selectedMuons[0]->photonIso(4);
                PUIso = selectedMuons[0]->puChargedHadronIso(4);

                
                MSPlot["leptonIso"]->Fill(relisomu, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonPt"]->Fill(selectedMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                leptonIso = relisomu;
            }

            //////////////////////////
            // Electron Based Plots //
            //////////////////////////

            for (Int_t selel =0; selel < selectedElectrons.size(); selel++ )
            {
                float reliso = selectedElectrons[selel]->relPfIso(3, 0.5);
                MSPlot["ElectronRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                leptonIso = reliso;
            }

            ////////////////////////////////////////////
            //       calculating HT rat and HTH       //
            ////////////////////////////////////////////
            if (debug) cout<<"HT rat and HTH"<<endl;
            HT = 0;
            float HT1M2L=0, H1M2L=0, HTbjets=0, HT2M=0, H2M=0, HT2L2J=0;
            sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order Jets wrt Pt for tuple output

            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
                MSPlot["BdiscBJetCand_CSV"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(),datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                //Event-level variables
                jetpt = selectedJets[seljet1]->Pt();
                HT = HT + jetpt;
                H = H +  selectedJets[seljet1]->P();
                if (seljet1 > 4  )  HTHi +=  selectedJets[seljet1]->Pt();
            }

            HTH = HT/H;
            HTRat = HTHi/HT;

            //////////////////////
            // Jets Based Plots //
            //////////////////////
            MSPlot["HTb_SelectedJets"]->Fill(HTb, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["HTRat"]->Fill(HTRat, datasets[d], true, Luminosity*scaleFactor);
            MSPlot["NbOfSelectedBJets"]->Fill(selectedMBJets.size(), datasets[d], true, Luminosity*scaleFactor);
            float met = mets[0]->Et();
            MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);

            if(nJets>5){
                MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            }

            MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);


            if(debug) cout<<"lepton vars"<<endl;
            ////////////////////////////////////////////
            //              Get lepton pt             //
            ////////////////////////////////////////////
            float selectedLeptonPt = 0 ;
            if(Muon){
                selectedLeptonPt = selectedMuons[0]->Pt();
                leptoneta = selectedMuons[0]->Eta();
                leptonphi = selectedMuons[0]->Phi();                
            }
            else if(Electron){
                selectedLeptonPt = selectedElectrons[0]->Pt();
                leptoneta = selectedElectrons[0]->Eta();
                leptonphi = selectedElectrons[0]->Phi();
            }
            ////////////////////////////////////////////
            //    Output special events to .txt       //
            ////////////////////////////////////////////
            if(debug) cout<<"special event output"<<endl;
            if(nJets > 7 && (dataSetName.find("Data") || dataSetName.find("data") || dataSetName.find("DATA")) ){
                //cout<<event->runId()  << " " << event->lumiBlockId() <<" " <<event->eventId() << "  jets "  << nJets <<"  nmtags "<<nMtags<<" muon pt "<<selectedMuons[0]->Pt()<<" 1stjetpt "<<selectedJets[0]->Pt()<<"  2ndjet pt "<<selectedJets[1]->Pt()<<endl;        

                eventlist <<event->runId()  << " " << event->lumiBlockId() <<" " <<event->eventId() << "  jets "  << nJets <<" nmtags "<<nMtags<<" muon pt "<<selectedLeptonPt<<" 1stjetpt "<<selectedJets[0]->Pt()<<"  2ndjet pt "<<selectedJets[1]->Pt()<<endl;        
                for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
                {
                    eventlist<<"  jet pt  "<<selectedJets[seljet1]->Pt()<<"   btag csv "<<selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags()<<endl;
                }
            }

            ////////////////////////////////////////////
            //    Jet variables + jet craneen         //
            ////////////////////////////////////////////
            tupjetfile->cd();
            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
                float jeteta = selectedJets[seljet1]->Eta();
                float jetphi = selectedJets[seljet1]->Phi();
                float csvDisc = selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
                float jetpT = selectedJets[seljet1]->Pt();
                float jetLepDR = 0;
                if (Muon){
                    jetLepDR = selectedJets[seljet1]->DeltaR(*selectedMuons[0]);
                }
                else if (Electron){
                    jetLepDR = selectedJets[seljet1]->DeltaR(*selectedElectrons[0]);
                }
                float jetvals[8] = {jetpT,csvDisc,jeteta,jetphi,jetLepDR,scaleFactor,normfactor,Luminosity};
                tupjet->Fill(jetvals);
            }

            // if(Muon){
            //     // muonpt  = selectedMuons[0]->Pt();
            //     // muoneta = selectedMuons[0]->Eta();
            //     leptonphi = selectedMuons[0]->Phi();
            // }
            // else if(Electron){
            //     // muonpt  = selectedElectrons[0]->Pt();
            //     // muoneta = selectedElectrons[0]->Eta();
            //     leptonphi = selectedElectrons[0]->Phi();

            // }

            ////////////////////////////////////////////
            //       Fill BDT & compute score         //
            ////////////////////////////////////////////
            if(debug) cout<<"event BDT computer"<<endl;
            float jet5Pt = 0;
            float jet6Pt = 0;
            if (EventBDTOn){
                if (nJets>5){
                    jet5Pt =  selectedJets[4]->Pt();
                    jet6Pt = selectedJets[5]->Pt();
                    //cout<<"5thjetpt "<<jet5Pt<<"  jet6pt: "<<jet6Pt<<endl;
                }
                eventBDT->fillVariables(diTopness, selectedLeptonPt, leptoneta, HTH, HTRat, HTb, nLtags, nMtags, nTtags, nJets, jet5Pt, jet6Pt);
            }
           
            BDTScore = 0 ;
            if(EventBDTOn){
                eventBDT->computeBDTScore();
                BDTScore = eventBDT->returnBDTScore();
            }

            ////////////////////////////////////////////
            //       Return variables for ntup        //
            ////////////////////////////////////////////
            if (selectedMBJets.size()>0){
                bjetpt = selectedMBJets[0]->Pt();
            }
            float firstjetpt = selectedJets[0]->Pt();
            float secondjetpt = selectedJets[1]->Pt();
            float nvertices = vertex.size();
            float angletoplep = 0;
            float angletop1top2 = 0;
            if(HadTopOn){
                angletop1top2 = hadronicTopReco->ReturnAnglet1t2();
                angletoplep = hadronicTopReco->ReturnAngletoplep();                
            }
            float nOrigJets = (float)selectedOrigJets.size();
            float jet5and6Pt = jet5Pt+jet6Pt;
            float vals[53] = {BDTScore,nJets,nOrigJets,nLtags,nMtags,nTtags,HT,selectedLeptonPt,leptoneta,bjetpt,HT2M,HTb,HTH,HTRat,HTX,SumJetMassX,diTopness,numOfbb,numOfcc,numOfll,ttbar_flav,scaleFactor,fleptonSF,btagWeight,btagWeightUp,btagWeightDown,lumiWeight,nvertices,normfactor,Luminosity,weight_0,weight_1,weight_2,weight_3,weight_4,weight_5,weight_6,weight_7,weight_8,met,angletop1top2,angletoplep,firstjetpt,secondjetpt,leptonIso,leptonphi,chargedHIso,neutralHIso,photonIso,PUIso,jet5Pt,jet6Pt,jet5and6Pt};
            tupfile->cd();
            tup->Fill(vals);
        } //End Loop on Events
        cout<<"Write files"<<endl;
        tupfile->cd();
        tup->Write();
        tupfile->Close();
        
        tupjetfile->cd();
        tupjet->Write();
        tupjetfile->Close();

        tupZfile->cd();
        tupZ->Write();
        tupZfile->Close();
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

    zPeakMaker->writeErase(fout, pathPNG);
    delete zPeakMaker;

    cutsTable->Calc_Write(postfix, dName, channelpostfix);
    delete cutsTable;

    if(fillingbTagHistos && bTagReweight && dataSetName.find("Data")==string::npos)    delete btwt;

    cout<<"TRIGGGG"<<endl;

    cout<<"preTrig: "<<preTrig<<"   postTrig: "<<postTrig<<endl;
    cout<<"********"<<endl;
    if(postTrig>0){
        cout<<"negative weight NormFactor: "<< ( (postTrig - (2*negWeights))/postTrig )<<endl;
    }

    delete trigger;

    fout->cd();
    TFile *foutmva = new TFile ("foutMVA.root","RECREATE");
    cout <<" after cd .."<<endl;

    string pathPNGJetCombi = pathPNG + "JetCombination/";
    mkdir(pathPNGJetCombi.c_str(),0777);
    if (HadTopOn){
        //if(TrainMVA)jetCombiner->Write(foutmva, true, pathPNGJetCombi.c_str());
        hadronicTopReco->WriteDiagnosticPlots(fout, pathPNG);
        delete hadronicTopReco;        
    }
    if(EventBDTOn){
        delete eventBDT;
    }

    if(jetTools) delete jetTools;


    //Output ROOT file
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
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


            // sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());

            // vector<double> ptList;
            // double jetptTemp;
            // for(Int_t csvJets = 2; csvJets<selectedJets.size(); csvJets++){
            //     jetptTemp = (double)selectedJets[csvJets]->Pt();
            //     ptList.push_back(jetptTemp);
            //     //cout<<csvJets<<"   ptlist "<<ptList[csvJets-2]<<endl;
            //     //selectedJets2.push_back(selectedJets[csvJets]);  //created array of selected jets without 2 highest CSVL btags
            //     //cout<<csvJets<<" jet pt "<<selectedJets2[csvJets-2]->Pt()<<"   "<<selectedJets[csvJets]->Pt()<<endl;
            // }

            // sort(selectedJets2.begin(),selectedJets2.end(),HighestPt()); //order Jets wrt Pt for tuple output

            // HT2L2J = HT - selectedJets[0]->Pt() - selectedJets[1]->Pt() - ptList[0] - ptList[1];    
            //cout<<"HT:  "<<HT<<"  "<<selectedJets[0]->Pt()<<"  "<<selectedJets[1]->Pt()<<"  "<<ptList[0]<<"  "<<ptList[1]<<"  HT2l2J"<<HT2L2J<<endl;    

            //HT - (2 highest CSVL btags) and (2 highest pt jets from remaining jets)

