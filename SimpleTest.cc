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
#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

#include "TopTreeAnalysisBase/Analysis/interface/CutsTable.h"
#include "TopTreeAnalysisBase/Analysis/interface/HadronicTopReco.h"
#include "TopTreeAnalysisBase/Analysis/interface/EventBDT.h"
#include "TopTreeAnalysisBase/Analysis/interface/Zpeak.h"

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
    if(argc < 14)
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

    ofstream eventlist;
    eventlist.open ("interesting_events_mu3.txt");

    int passed = 0;
    int mupassed = 0.;
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
    // Configuration
    ///////////////////////////////////////

    //Setting Lepton Channels 
    bool SingleLepton = true;
    bool Muon = true;
    bool Electron = false;
    bool HadTopOn = false;
    bool EventBDTOn = false;
    bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
    bool bx25 = false;
    bool split_ttbar = false;
    bool debug = false;
    string MVAmethod = "BDT"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
    float Luminosity = 7.342; //pb^-1

    if(Muon){
        cout<<" ***** USING SINGLE MUON CHANNEL  ******"<<endl;
        channelpostfix = "_Mu";
    }
    else if(Electron){
        cout<<" ***** Using SINGLE ELECTRON CHANNEL *****"<<endl;
        channelpostfix = "_El";
    }
    else
    {
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

    ////////////////////////////////
    //  Load datasets
    ////////////////////////////////

    TTreeLoader treeLoader;
    vector < Dataset* > datasets;    cout << " - Creating Dataset ..." << endl;
    Dataset* theDataset = new Dataset(dName, dTitle, true, color, ls, lw, normf, xSect, vecfileNames);
    theDataset->SetEquivalentLuminosity(EqLumi);
    datasets.push_back(theDataset);


    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////

    cout <<"found sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    dataSetName = theDataset->Name();
    if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
    {
        Luminosity = theDataset->EquivalentLumi();      cout <<"found DATA sample with equivalent lumi "<<  theDataset->EquivalentLumi() <<endl;
    }

    cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;

    //vector of objects
    cout << " - Variable declaration ..." << endl;
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
    vector < TRootMET* >      mets;

    CutsTable *cutsTableMu = new CutsTable(Muon, Electron);
    cutsTableMu->AddSelectionsMuons();
    cutsTableMu->CreateTable(datasets, Luminosity);    

    //Global variable
    TRootEvent* event = 0;

    ////////////////////////////////////////////////////////////////////
    ////////////////// MultiSample plots  //////////////////////////////
    ////////////////////////////////////////////////////////////////////

    MSPlot["NbOfVertices"]          = new MultiSamplePlot(datasets, "NbOfVertices", 60, 0, 60, "Nb. of vertices");

    /////////////////////////////////
    //       Loop on datasets      //
    /////////////////////////////////
    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;    treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
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

      
        //////////////////////////////////////////////////
        // Loop on events
        /////////////////////////////////////////////////

        int itrigger = -1, previousRun = -1, start = 0;

        unsigned int ending = datasets[d]->NofEvtsToRunOver();    cout <<"Number of events = "<<  ending  <<endl;

        int event_start = startEvent;

        if (dataSetName == "Data") TrainMVA=false;
        if (debug) cout << " - Loop over events " << endl;

        float BDTScore, MHT, MHTSig, STJet,muoneta, muonpt, leptonphi, electronpt, electroneta, bjetpt, EventMass, EventMassX, SumJetMass, SumJetMassX, H, HX;
        float HTHi, HTRat, HT, HTX, HTH, HTXHX, sumpx_X, sumpy_X, sumpz_X, sume_X, sumpx, sumpy, sumpz, sume, jetpt, PTBalTopEventX, PTBalTopSumJetX, PTBalTopMuMet;

        double end_d = ending;
        if(endEvent > ending)
            end_d = ending;
        else
            end_d = endEvent;

        cout <<"Will run over "<<  end_d<< " events..."<<endl;    cout <<"Starting event = = = = "<< event_start  << endl;

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
        vector<TRootJet*>      selectedLBJets; //CSVL btags
        vector<TRootJet*>      selectedMBJets; //CSVM btags
        vector<TRootJet*>      selectedTBJets; //CSVT btags
        vector<TRootJet*>      selectedLightJets;


        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////
        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
            double ievt_d = ievt;

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC 
                << " ("<<100*(ievt-start)/(ending-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);  //load event
            float weight_0 = event->weight0();

            if(nlo)
            {
                if(weight_0 < 0.0)
                {
                    scaleFactor = -1.0;  //Taking into account negative weights in NLO Monte Carlo
                    negWeights++;
                }
            }

            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////

            // Apply trigger selection
            // trigged = treeLoader.EventTrigged (itrigger);
            bool trigged = true;  // Disabling the HLT requirement

            cutsTableMu->FillTableMuons(d, scaleFactor, init_muons);

            // Declare selection instance
            Run2Selection r2selection(init_jets, init_muons, init_electrons, mets);

            // Define object selection cuts

            int nMu, nEl, nLooseMu, nLooseEl; //number of (loose) muons/electrons
            nMu = 0, nEl = 0, nLooseMu=0, nLooseEl=0;

            selectedJets                                        = r2selection.GetSelectedJets(); // ApplyJetId
            selectedMuons                                       = r2selection.GetSelectedMuons();
            selectedExtraMuons                                  = r2selection.GetSelectedMuons(20, 2.4, 0.20);     
            selectedElectrons                                   = r2selection.GetSelectedElectrons("Tight", "PHYS14", true); // VBTF ID       
            selectedExtraElectrons                              = r2selection.GetSelectedElectrons("Loose", "PHYS14", true);

            nMu = selectedMuons.size(); //Number of Muons in Event
            nEl = selectedElectrons.size(); //Number of Electrons in Event
            nLooseMu = selectedExtraMuons.size();   //Number of loose muons      
            nLooseEl = selectedExtraElectrons.size(); //Number of loose muons
            int nJets = selectedJets.size();

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
                        HTb += selectedJets[seljet]->Pt();
                        selectedMBJets.push_back(selectedJets[seljet]);
                        if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.898)
                        {
                            selectedTBJets.push_back(selectedJets[seljet]);
                        }
                    }
                }
            }
            float nLtags = selectedLBJets.size(); //Number of CSVL tags in Event (includes jets that pass CSVM & CSVT)
            float nMtags = selectedMBJets.size(); //Number of CSVM tags in Event (includes jets that pass CSVT)
            float nTtags = selectedTBJets.size(); //Number of CSVT tags in Event 

            selectedLBJets.clear();
            selectedMBJets.clear();
            selectedTBJets.clear();
            selectedLightJets.clear();

            //////////////////////
            // Sync'ing cutflow //
            //////////////////////

            if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;

            // Apply primary vertex selection
            bool isGoodPV = r2selection.isPVSelected(vertex, 4, 24., 2);
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV << " TriggerBit: " << trigged <<endl;
            if (debug) cin.get();

            weightCount += scaleFactor;
            eventCount++;



            /////////////////////////////////
            // Applying baseline selection //
            /////////////////////////////////

            //Filling Histogram of the number of vertices before Event Selection
            MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

            if (!isGoodPV) continue; // Check that there is a good Primary Vertex

            if(Muon && nMu==1){
                mupassed++;
            }

            //Apply the lepton, jet, btag and HT & MET selections
            if (Muon)
            {
                if  (  !( nMu == 1 && nEl == 0 && nLooseMu == 1 && nJets>=4 && nMtags >=2)) continue; // Muon Channel Selection

            }
            else if(Electron){
                if  (  !( nMu == 0 && nEl == 1 && nLooseEl == 1 && nJets>=4 && nMtags >=2)) continue; // Electron Channel Selection
            }
            else{
                cerr<<"Correct Channel not selected."<<endl;
                exit(1);
            }

            passed++;
        }
        //important: free memory
        treeLoader.UnLoadDataset();
    } //End Loop on Datasets

    cutsTableMu->Calc_Write("Mu", dName, "MuonCuts");
    eventlist.close();

    /////////////
    // Writing //
    /////////////



    cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
    cout << "********************************************" << endl;
    cout << "           End of the program !!            " << endl;
    cout << "********************************************" << endl;

    cout <<"n events passed  =  "<<passed <<endl;
    cout <<"n events mupassed  =  "<<mupassed <<endl;

    cout <<"n events with negative weights = "<<negWeights << endl;    
    
    return 0;
}