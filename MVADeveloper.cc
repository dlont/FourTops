//////////////////////////////////////////////////////////////////////////////
////     Stand-alone code to train the event-level BDT
////////////////////////////////////////////////////////////////////////////////////
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

struct HighestCVSBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_combinedInclusiveSecondaryVertexV2BJetTags() > j2->btag_combinedInclusiveSecondaryVertexV2BJetTags();
    }
};

float Sphericity(vector<TLorentzVector> parts );
float Centrality(vector<TLorentzVector> parts);

int main ()
{


    clock_t start = clock();
    string xmlFileName = "config/Run2DiLepton_TOPTREES_BDTTrain.xml";

    const char *xmlfile = xmlFileName.c_str();
    cout << "used config file: " << xmlfile << endl;

/////////////////////////////
/// AnalysisEnvironment
/////////////////////////////

    AnalysisEnvironment anaEnv;
    cout<<" - Loading environment ..."<<endl;
    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    int verbose = 2;//anaEnv.Verbose;

////////////////////////////////
//  Load datasets
////////////////////////////////

    TTreeLoader treeLoader;
    vector < Dataset* > datasets;
    cout << " - Loading datasets ..." << endl;
    treeLoader.LoadDatasets (datasets, xmlfile);
    float Luminosity = 2000.0; //pb^-1??
    vector<string> MVAvars;

    //A few bools to steer the MassReco and Event MVAs
    string MVAmethod = "BDT";

    int nTags ;
    double H, HTH, HT, HTHi,HTb,HTRat,HT2M;
    int jet_flavor;
    double JetPt, JetEta;
    double SF_tag =1.;
    float workingpointvalue = 0.679; //working points updated to 2012 BTV-POG recommendations.

    cout <<"Instantiating jet combiner..."<<endl;

    JetCombiner* jetCombiner = new JetCombiner(false, 1000.0, datasets, MVAmethod, false);
    cout <<"Instantiated jet combiner..."<<endl;


    //vector of objects
    cout << " - Variable declaration ..." << endl;
    vector < TRootVertex* >   vertex;
    vector < TRootMuon* >     init_muons;
    vector < TRootElectron* > init_electrons;
    vector < TRootJet* >      init_jets;
    vector < TRootMET* >      mets;

    bool debug = false;
    bool singlelep = false;
    bool dilep = true;
    bool Muon = true;
    bool Electron = true;
    int nPassed = 0;

    BTagWeightTools * bTool = new BTagWeightTools("SFb-pt_NOttbar_payload_EPS13.txt", "CSVM") ;

    int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
    cout << "dobTagEffShift: " << dobTagEffShift << endl;

    int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
    cout << "domisTagEffShift: " << domisTagEffShift << endl;


    //Global variable
    TRootEvent* event = 0;

    MVATrainer* Eventtrainer_;

    if (dilep)
    {
        Eventtrainer_ = new MVATrainer("BDT","MasterMVA_MuEl_25thAugust", "MVA/MasterMVA_MuEl_25thAugust.root");
    }
    else if (singlelep)
    {
        Eventtrainer_ = new MVATrainer("BDT","MasterMVA_SingleMuon_24thMarch", "MVA/MasterMVA_SingleMuon_24thMarch.root");
    }


    //define object containers
    vector<TRootJet*>      selectedLBJets;
    vector<TRootJet*>      selectedMBJets;
    vector<TRootJet*>      selectedTBJets;
    vector<TRootElectron*> selectedElectrons;
    vector<TRootPFJet*>    selectedJets;
    vector<TRootJet*>      selectedLightJets;
    vector<TRootJet*>      selectedBJets;
    vector<TRootMuon*>     selectedMuons;
    vector<TRootElectron*> selectedExtraElectrons;
    vector<TRootMuon*>     selectedExtraMuons;
    vector<TRootMCParticle*> mcParticlesMatching_;
    vector<TRootPFJet*>  selectedJets2ndPass;


    if (dilep)
    {
        Eventtrainer_->bookWeight("Weight");
        Eventtrainer_->bookInputVar("topness");
        Eventtrainer_->bookInputVar("muonpt");
        Eventtrainer_->bookInputVar("muoneta");
        Eventtrainer_->bookInputVar("HTH");
        Eventtrainer_->bookInputVar("HTRat");
        Eventtrainer_->bookInputVar("HTb");
        Eventtrainer_->bookInputVar("nLtags");
        Eventtrainer_->bookInputVar("nMtags");
        Eventtrainer_->bookInputVar("nTtags");
        Eventtrainer_->bookInputVar("nJets");
        Eventtrainer_->bookInputVar("Jet3Pt");
        Eventtrainer_->bookInputVar("Jet4Pt");
        Eventtrainer_->bookInputVar("HT2M");
        Eventtrainer_->bookInputVar("EventSph");
//        Eventtrainer_->bookInputVar("EventCen");
        Eventtrainer_->bookInputVar("DiLepSph");
//        Eventtrainer_->bookInputVar("DiLepCen");
//        Eventtrainer_->bookInputVar("TopDiLepSph");
//        Eventtrainer_->bookInputVar("TopDiLepCen");
    }
    else if (singlelep)
    {

        // Eventtrainer_->bookWeight("Weight");
        Eventtrainer_->bookInputVar("multitopness");
        Eventtrainer_->bookInputVar("muonpt");
        Eventtrainer_->bookInputVar("muoneta");
        Eventtrainer_->bookInputVar("HTH");
        Eventtrainer_->bookInputVar("HTRat");
        Eventtrainer_->bookInputVar("HTb");
        Eventtrainer_->bookInputVar("nLtags");
        Eventtrainer_->bookInputVar("nMtags");
        Eventtrainer_->bookInputVar("nTtags");
        Eventtrainer_->bookInputVar("nJets");
        Eventtrainer_->bookInputVar("Jet5Pt");
        Eventtrainer_->bookInputVar("Jet6Pt");
    }


    TRootGenEvent* genEvt = 0;
    int event_start = 0;

    cout <<"  N datsets = "<< datasets.size () << endl;


    for (unsigned int d = 0; d < datasets.size (); d++)
    {
        cout <<"Beginning dataset loop: " << datasets[d]->Name() << " with " << datasets[d]->NofEvtsToRunOver() << " events." << endl;
        unsigned int ending = datasets[d]->NofEvtsToRunOver();
        int start = 0;
        cout <<"Number of events in total dataset = "<< ending <<endl;
        event_start = 0;
        nPassed = 0;

        cout << "Dataset loop..." << endl;
        cout <<" dataset name "  << datasets[d]->Name() <<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset

        for (unsigned int ievt = 0; ievt < ending; ievt++)
        {

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(100000-event_start)<<"%)"<<flush<<"\r"<<endl;
            }
            if(nPassed >= 7000) continue;

            mcParticlesMatching_.clear();
            genEvt = treeLoader.LoadGenEvent(ievt,false);
            treeLoader.LoadMCEvent(ievt, genEvt, 0, mcParticlesMatching_,false);
            HT =0;
            H =0;
            HTH=0;
            HTHi =0;
            nTags =0;
            HTb = 0;
            HT2M = 0;
            selectedLBJets.clear();
            selectedMBJets.clear();
            selectedTBJets.clear();

            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, false);  //load event

// Declare selection instance
            Run2Selection selection(init_jets, init_muons, init_electrons, mets);

            if (debug)cout<<"Getting Jets"<<endl;
            selectedJets                                        = selection.GetSelectedJets(); // Relying solely on cuts defined in setPFJetCuts()
            if (debug)cout<<"Getting Tight Muons"<<endl;
            selectedMuons                                       = selection.GetSelectedDiMuons();
            selectedExtraMuons = selection.GetSelectedMuons(20, 2.4, 0.20);
            if (debug)cout<<"Getting Loose Electrons"<<endl;
            selectedElectrons                                   = selection.GetSelectedElectrons("Loose","PHYS14",true); // VBTF ID

            int nJets = selectedJets.size(); //Number of Jets in Event
            int nMu = selectedMuons.size(); //Number of Muons in Event
            int nEl = selectedElectrons.size(); //Number of Electrons in Event

            vector<TLorentzVector> selectedMuonsTLV_JC, selectedElectronsTLV_JC, selectedJetsTLV;

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

            else if(nEl == 2 && nMu == 0 && Electron && !Muon)
            {
                lep1 = selectedElectronsTLV_JC[0];
                lep2 = selectedElectronsTLV_JC[1];
            }
            else if(nEl == 1 && nMu == 1 && Electron && Muon)
            {
                lep1 = selectedMuonsTLV_JC[0];
                lep2 = selectedElectronsTLV_JC[0];
            }
            sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());

            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
                selectedJetsTLV.push_back(*selectedJets[seljet1]);
                HT +=  selectedJets[seljet1]->Pt();
                H +=  selectedJets[seljet1]->P();


                if (dilep)
                {
                    if(seljet1>=2)
                    {
                        HT2M = HT2M + selectedJets[seljet1]->Pt();
                    }
                    if (seljet1 > 2 )     HTHi +=  selectedJets[seljet1]->Pt();
                }
                else if (singlelep)
                {

                    if (seljet1 > 4 )     HTHi +=  selectedJets[seljet1]->Pt();
                }

            }

            HTRat = HTHi/HT;
            HTH = HT/H;

            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
                if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue)
                {
                    selectedMBJets.push_back(selectedJets[seljet]);
                    HTb += selectedJets[seljet]->Pt();
                }

                if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.244) selectedLBJets.push_back(selectedJets[seljet]);
                if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.898) selectedTBJets.push_back(selectedJets[seljet]);
            }



            int nLtags = selectedLBJets.size(); //Number of CSVM tags in Event
            int nMtags = selectedMBJets.size(); //Number of CSVM tags in Event
            int nTtags = selectedTBJets.size(); //Number of CSVM tags in Event
            int nLooseMu = selectedExtraMuons.size();
            bool ZVeto = false;

            if (dilep && Muon && Electron)
            {
                if  (  !( nMu == 1 && nEl == 1 )) continue; // Muon-Electron Channel Selection
                if  (  !( selectedMuons[0]->Pt() <= 1000 && selectedElectrons[0]->Pt() <= 1000 )) continue; // Muon-Electron Channel Selection
                if  (  !( nJets >= 4 )) continue; // Muon-Electron Channel Selection
                if  (  !( nMtags >= 2 )) continue; // Muon-Electron Channel Selection
                if  (  !( HT >= 500 )) continue; // Muon-Electron Channel Selection
                if  (  !( mets[0]->Et() <= 1000 )) continue; // Muon-Electron Channel Selection
            }
            else if (dilep && Muon && !Electron)
            {
                if  (  !( nMu == 2 && nEl == 0 )) continue; // Muon-Electron Channel Selection
                if  (  !( selectedMuons[0]->Pt() <= 1000 && selectedMuons[1]->Pt() <= 1000 )) continue; // Muon-Electron Channel Selection
                TLorentzVector diMu = (TLorentzVector)(*selectedMuons[0]) + (TLorentzVector)(*selectedMuons[1]);
                float diMuMass = diMu.M();
                if(diMuMass < 20 || (diMuMass > 76 && diMuMass < 106)) ZVeto = true;
                if(ZVeto) continue;
                if  (  !( nJets >= 4 )) continue; // Muon-Electron Channel Selection
                if  (  !( nMtags >= 2 )) continue; // Muon-Electron Channel Selection
                if  (  !( HT >= 500 )) continue; // Muon-Electron Channel Selection
                if  (  !( mets[0]->Et() <= 1000 )) continue; // Muon-Electron Channel Selection
            }
            else if (dilep && !Muon && Electron)
            {
                if  (  !( nMu == 0 && nEl == 2 )) continue; // Muon-Electron Channel Selection
                if  (  !( selectedElectrons[0]->Pt() <= 1000 && selectedElectrons[1]->Pt() <= 1000 )) continue; // Muon-Electron Channel Selection
                TLorentzVector diEl = (TLorentzVector)(*selectedElectrons[0]) + (TLorentzVector)(*selectedElectrons[1]);
                float diElMass = diEl.M();
                if(diElMass < 20 || (diElMass > 76 && diElMass < 106)) ZVeto = true;
                if(ZVeto) continue;
                if  (  !( nJets >= 4 )) continue; // Muon-Electron Channel Selection
                if  (  !( nMtags >= 2 )) continue; // Muon-Electron Channel Selection
                if  (  !( HT >= 500 )) continue; // Muon-Electron Channel Selection
                if  (  !( mets[0]->Et() <= 1000 )) continue; // Muon-Electron Channel Selection
            }
            else if (singlelep)
            {
                if  (  !( nMu == 1 && nEl == 0 && nLooseMu ==1  )) continue; // Muon + jets Channel Selection
                if  (  !( nJets >= 6 )) continue; // Muon + jets Channel Selection
                if  (  !( nMtags >= 2 )) continue; // Muon + jets Channel Selection
            }

            if(debug) cout <<"event passed - "<< datasets[d]->Name() <<" njets  "  << nJets  << " ntags "<< nMtags  <<  endl;
            nPassed++;
            if(nPassed%1000 == 0)
            {
                std::cout<<nPassed<<" events have passed baseline selection"<<endl;
            }

            if (dilep && !Muon && Electron) {
                jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticlesMatching_, selectedJets, &selectedElectronsTLV_JC[0], genEvt, 1.);
            }
            else {
                jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticlesMatching_, selectedJets, &selectedMuonsTLV_JC[0], genEvt, 1.);
            }

            pair<float, vector<unsigned int> > MVAvals1;
            MVAvals1 = jetCombiner->getMVAValue(MVAmethod, 1);

            double topness = MVAvals1.first;

            selectedJets2ndPass.clear();

            pair<float, vector<unsigned int> > MVAvals2ndPass;

//make vector of jets excluding thise selected by 1st pass of mass reco
            for (Int_t seljet3 =0; seljet3 < selectedJets.size(); seljet3++ )
            {
                if (seljet3 == MVAvals1.second[0] || seljet3 == MVAvals1.second[1] || seljet3 == MVAvals1.second[2])
                {
                    //   MVASelJets1.push_back(selectedJets[seljet3]);
                    continue;
                }
                selectedJets2ndPass.push_back(selectedJets[seljet3]);
            }

            double multitopness = 0;
            double muonpt, muoneta;
            if(nJets >=6)
            {
                if (dilep && !Muon && Electron) {
                jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticlesMatching_, selectedJets, &selectedElectronsTLV_JC[0], genEvt, 1.);
                muonpt = selectedElectronsTLV_JC[0].Pt();
                muoneta = selectedElectronsTLV_JC[0].Eta();
            }
            else {
                jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticlesMatching_, selectedJets, &selectedMuonsTLV_JC[0], genEvt, 1.);
                muonpt = selectedMuons[0]->Pt();
                muoneta = selectedMuons[0]->Eta();
            }
                MVAvals2ndPass = jetCombiner->getMVAValue(MVAmethod, 1);

                multitopness = MVAvals2ndPass.first;
            }

            //  cout <<" # leftover jets " <<  selectedJets2ndPass.size()     << " multitopness  =  "<<  multitopness   <<endl;
            float weight_0 = event->weight0();
            float scaleFactor;
            if(weight_0 < 0) scaleFactor = -1.0;
            else scaleFactor = 1.0;


            if (dilep)
            {
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
                topDiLepSystemTLV.push_back(*selectedJets[MVAvals1.second[0]]);
                topDiLepSystemTLV.push_back(*selectedJets[MVAvals1.second[1]]);
                topDiLepSystemTLV.push_back(*selectedJets[MVAvals1.second[2]]);

                float tSph = Sphericity(selectedParticlesTLV), tCen = Centrality(selectedParticlesTLV);
                float dSph = Sphericity(diLepSystemTLV), dCen = Centrality(diLepSystemTLV);
                float tdSph = Sphericity(topDiLepSystemTLV), tdCen = Centrality(topDiLepSystemTLV);

                if(datasets[d]->Name()=="ttttNLO_SUS_MuElSkim")
                {
                    Eventtrainer_->Fill("S","Weight",scaleFactor);
                    Eventtrainer_->Fill("S","topness",topness );
                    Eventtrainer_->Fill("S","muonpt",muonpt);
                    Eventtrainer_->Fill("S","muoneta",muoneta);
                    Eventtrainer_->Fill("S","HTb", HTb);
                    Eventtrainer_->Fill("S","HTH", HTH);
                    Eventtrainer_->Fill("S","HTRat", HTRat);
                    Eventtrainer_->Fill("S","nLtags",nLtags  );
                    Eventtrainer_->Fill("S","nMtags",nMtags  );
                    Eventtrainer_->Fill("S","nTtags",nTtags  );
                    Eventtrainer_->Fill("S","nJets",selectedJets.size()  );
                    Eventtrainer_->Fill("S","Jet3Pt",selectedJets[2]->Pt());
                    Eventtrainer_->Fill("S","Jet4Pt",selectedJets[3]->Pt());
                    Eventtrainer_->Fill("S","HT2M",HT2M  );
                    Eventtrainer_->Fill("S","EventSph",tSph  );
//                    Eventtrainer_->Fill("S","EventCen",tCen  );
                    Eventtrainer_->Fill("S","DiLepSph",dSph  );
//                    Eventtrainer_->Fill("S","DiLepCen",dCen  );
//                    Eventtrainer_->Fill("S","TopDiLepSph",tdSph  );
//                    Eventtrainer_->Fill("S","TopDiLepCen",tdCen  );
                }

                if(datasets[d]->Name()!="ttttNLO_SUS_MuElSkim")
                {
                    Eventtrainer_->Fill("B","Weight", scaleFactor);
                    Eventtrainer_->Fill("B","topness",topness);
                    Eventtrainer_->Fill("B","muonpt",muonpt);
                    Eventtrainer_->Fill("B","muoneta",muoneta);
                    Eventtrainer_->Fill("B","HTb", HTb);
                    Eventtrainer_->Fill("B","HTRat", HTRat);
                    Eventtrainer_->Fill("B","HTH", HTH);
                    Eventtrainer_->Fill("B","nLtags", nLtags );
                    Eventtrainer_->Fill("B","nMtags", nMtags );
                    Eventtrainer_->Fill("B","nTtags", nTtags );
                    Eventtrainer_->Fill("B","nJets", selectedJets.size() );
                    Eventtrainer_->Fill("B","Jet3Pt", selectedJets[2]->Pt() );
                    Eventtrainer_->Fill("B","Jet4Pt", selectedJets[3]->Pt() );
                    Eventtrainer_->Fill("B","HT2M",HT2M  );
                    Eventtrainer_->Fill("B","EventSph",tSph  );
//                    Eventtrainer_->Fill("B","EventCen",tCen  );
                    Eventtrainer_->Fill("B","DiLepSph",dSph  );
//                    Eventtrainer_->Fill("B","DiLepCen",dCen  );
//                    Eventtrainer_->Fill("B","TopDiLepSph",tdSph  );
//                    Eventtrainer_->Fill("B","TopDiLepCen",tdCen  );
                }

            }
            else if (singlelep)
            {

                if(datasets[d]->Name()=="NP_overlay_TTTT")
                {
                    //        Eventtrainer_->FillWeight("S","Weight",scaleFactor);
                    Eventtrainer_->Fill("S","multitopness",multitopness );
                    Eventtrainer_->Fill("S","muonpt",muonpt);
                    Eventtrainer_->Fill("S","muoneta",muoneta);
                    Eventtrainer_->Fill("S","HTb", HTb);
                    Eventtrainer_->Fill("S","HTH", HTH);
                    Eventtrainer_->Fill("S","HTRat", HTRat);
                    Eventtrainer_->Fill("S","nLtags",nLtags  );
                    Eventtrainer_->Fill("S","nMtags",nMtags  );
                    Eventtrainer_->Fill("S","nTtags",nTtags  );
                    Eventtrainer_->Fill("S","nJets",selectedJets.size()  );
                    Eventtrainer_->Fill("S","Jet5Pt",selectedJets[4]->Pt());
                    Eventtrainer_->Fill("S","Jet6Pt",selectedJets[5]->Pt());
                }

                if(datasets[d]->Name()=="TTJets")
                {

                    cout <<" fill ttjets tree"<<endl;

                    //	   Eventtrainer_->FillWeight("B","Weight", scaleFactor);
                    Eventtrainer_->Fill("B","multitopness",multitopness);
                    Eventtrainer_->Fill("B","muonpt",muonpt);
                    Eventtrainer_->Fill("B","muoneta",muoneta);
                    Eventtrainer_->Fill("B","HTb", HTb);
                    Eventtrainer_->Fill("B","HTRat", HTRat);
                    Eventtrainer_->Fill("B","HTH", HTH);
                    Eventtrainer_->Fill("B","nLtags", nLtags );
                    Eventtrainer_->Fill("B","nMtags", nMtags );
                    Eventtrainer_->Fill("B","nTtags", nTtags );
                    Eventtrainer_->Fill("B","nJets", selectedJets.size() );
                    Eventtrainer_->Fill("B","Jet5Pt", selectedJets[4]->Pt() );
                    Eventtrainer_->Fill("B","Jet6Pt", selectedJets[5]->Pt() );
                }

            }


        }

    }
    Eventtrainer_->TrainMVA("Random","",0,0,"",0,0,"_MuElAugust25th2015", true);


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
