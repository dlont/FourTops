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
 if (varlist.is_open())  while ( getline (varlist,line) ) vars.push_back(line);
    varlist.close();
 for (int i=0;  i< vars.size(); i++  )  cout << vars[i] << '\n';

 
//to be replaced with input from text file
MSPlot["HT"] = new MultiSamplePlot(datasets, "HT", 20, 0, 2500, "HT");

  for (unsigned int d = 0; d < datasets.size(); d++){
  samp_name =  datasets[d]->Filenames();
  TFile *f = new TFile(samp_name[0].c_str());
  TNtuple * tup  = (TNtuple*)f->Get("Craneen__MuEl");
  int nEvents = tup->GetEntries();
  Float_t HT, NormFactor, GenWeight;
  tup->SetBranchAddress("HT",&HT);
  tup->SetBranchAddress("NormFactor",&NormFactor);
  tup->SetBranchAddress("GenWeight",&GenWeight);


  for (int ev = 0; ev < nEvents; ev++){
    tup->GetEntry(ev);
    MSPlot["HT"]->Fill(HT, datasets[d], false, NormFactor*Luminosity );
    //  if (debug)cout <<" here 4 looping events "<< HT << endl;

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
