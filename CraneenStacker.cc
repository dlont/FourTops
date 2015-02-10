///////////////////////////////////////////////////////////////////////
//// This is a macro that produces publication-level plots         ////
////   from the ultra-light "craneen" ntuples produced by          ////
////   the main analysis macros.                                   ////
////                                                               ////
////   This macro will favour speed over flexibility as it is meant////
////   to be the final macro ran O(100) times in order to perfect  ////
////   the plots etc, NOT implement changes that will affect the   ////
////   results like SFs etc...                                     ////
////                                                               ////
////    James                                                      ////
/////////////////////////////////////////////////////////////////////// 
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "CraneenStacker.h"
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

using namespace std;

int main ()
{

  TFile * f1;
  bool debug = true;

  string xmlFileName = "Craneens.xml"; 
  //  string xmlFileName = "config/Run2sgluon_Samples.xml"; 
  const char *xmlfile = xmlFileName.c_str();

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);

  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  float Luminosity = 5000.0; //pb^-1??
 
  vector<string> samp_name;
  vector<plot_params> plots;
 
  //read in vector of craneens
  string craneen_dir = "Craneens_MuEl/Craneens6_2_2015/merged/";  
  //first read in list of samples 
  ifstream samplist;
  samplist.open ("samples_emu.txt");
  string line;

  vector < string > samps;

  if (samplist.is_open()){ 
    while ( getline (samplist,line) ) {
      samps.push_back(line);
    }
}
    samplist.close();
 int nsamps = samps.size();

 //read in list of variables for plotting
  ifstream varlist;
  varlist.open ("variables_emu.txt");
  vector < string > vars;
  if (varlist.is_open()){ 
    while ( getline (varlist,line) ){

      stringstream   linestream(line);
      plot_params param; 
      string var1;
      int var2;
      int var3;
      int var4;

      getline(linestream, var1, '\t'); 
      linestream >> var2 >>  var3 >> var4;

      param.name = var1;
      param.nbins = var2;
      param.xlo = var3;
      param.xhi = var4;

      cout << "var1 " << var1 <<  "  var2 " << var2 << " var3 " << var3 << " var4 " << var4 << endl;

      MSPlot[var1] = new MultiSamplePlot(datasets, var1, var2, var3, var4, var1);

    vars.push_back(line);
   
    }
  }

 varlist.close();

  //  for (int i=0;  i< vars.size(); i++  )  cout << vars[i] << '\n';

 //to be replaced with input from text file
 // MSPlot["HT"] = new MultiSamplePlot(datasets, "HT", 20, 0, 2500, "HT");

  for (unsigned int d = 0; d < datasets.size(); d++){
  samp_name =  datasets[d]->Filenames();
  TFile *f = new TFile(samp_name[0].c_str());
  TNtuple * tup  = (TNtuple*)f->Get("Craneen__MuEl");
  int nEvents = tup->GetEntries();
  Float_t HT, NormFactor, GenWeight, HT2M, nMtags,LeadingMuonPt, LeadingElectronPt, LeadingBJetPt, PU, ScaleFactor, nJets;

    tup->SetBranchAddress("NormFactor",&NormFactor);
    tup->SetBranchAddress("GenWeight",&GenWeight);
    tup->SetBranchAddress("HT",&HT);
    tup->SetBranchAddress("HT2M",&HT2M);
    tup->SetBranchAddress("nMtags",&nMtags);
    tup->SetBranchAddress("nJets",&nJets);
    tup->SetBranchAddress("LeadingMuonPt",&LeadingMuonPt);
    tup->SetBranchAddress("LeadingElectronPt",&LeadingElectronPt);
    tup->SetBranchAddress("LeadingBJetPt",&LeadingBJetPt);
    tup->SetBranchAddress("PU",&PU);

  for (int ev = 0; ev < nEvents; ev++){
    tup->GetEntry(ev);
    MSPlot["PU"]->Fill(PU, datasets[d], false, NormFactor*Luminosity );
    MSPlot["HT"]->Fill(HT, datasets[d], false, NormFactor*Luminosity );
    MSPlot["HT2M"]->Fill(HT2M, datasets[d], false, NormFactor*Luminosity );
    MSPlot["nMtags"]->Fill(nMtags, datasets[d], false, NormFactor*Luminosity );
    MSPlot["nJets"]->Fill(nJets, datasets[d], false, NormFactor*Luminosity );
    MSPlot["LeadingMuonPt"]->Fill(LeadingMuonPt, datasets[d], false, NormFactor*Luminosity );
    MSPlot["LeadingElectronPt"]->Fill(LeadingElectronPt, datasets[d], false, NormFactor*Luminosity );
    MSPlot["LeadingBJetPt"]->Fill(LeadingBJetPt, datasets[d], false, NormFactor*Luminosity );
 
}
  f1 = new TFile("StackCraneens.root", "RECREATE");

}

 for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
{

string name = it->first;
MultiSamplePlot *temp = it->second;

temp->Draw("CMSPlot", 1, true, true, true, 1);
void Write(TFile* file, string label = string(""), bool savePNG = false, string pathPNG = string(""), string ext = string("png"));
temp->Write(f1, name, true, ".", "pdf");


}

}
