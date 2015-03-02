////////////////////////////////////////////////////////////////////////
////   This is a macro that produces publication-level plots        ////
////   from the ultra-light "craneen" ntuples produced by           ////
////   the main analysis macros.                                    ////
////                                                                ////
////   This macro will favour speed over flexibility as it is meant ////
////   to be the final macro ran O(100) times in order to perfect   ////
////   the plots etc, NOT implement changes that will affect the    ////
////   results like SFs etc...                                      ////
////                                                                ////
////    James                                                       ////
////////////////////////////////////////////////////////////////////////
 
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

/// MultiSamplePlots
map<string,MultiSamplePlot*> MSPlot;
map<string,TH1F*> histos1D_Up;
map<string,TH1F*> histos1D_Down;
map<string,TH1F*> histos1D;

using namespace std;

int main ()
{

  string name;
  string sys_name;
  string dir_name;
  TFile * f1;
  TFile * f2;
  double nom_norm;
  bool debug = true;
  double scalefac;
  bool write_sys = false;
  TH1F *hnom;
  // *hup, *hdown;

  string xmlFileName;

  if (write_sys){
  xmlFileName = "config/Craneens_sgluons_sys.xml"; 
}else{
  xmlFileName = "config/Craneens_sgluons.xml"; 
}


  const char *xmlfile = xmlFileName.c_str();

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);

  TTreeLoader treeLoader; 
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  float Luminosity = 5000.0; //pb^-1??
 
  string dataset_name;
  vector<string> samp_name;
  // vector<plot_params> plots;
 
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
      //  plot_params param; 
      string var1;
      int var2;
      int var3;
      int var4;

      getline(linestream, var1, '\t'); 
      linestream >> var2 >>  var3 >> var4;
      MSPlot[var1] = new MultiSamplePlot(datasets, var1, var2, var3, var4, var1);
      histos1D_Up[var1] = new TH1F(var1.c_str(), var1.c_str(), var2, var3, var4);
      histos1D_Down[var1] = new TH1F(var1.c_str(), var1.c_str(), var2, var3, var4);
      histos1D[var1] = new TH1F(var1.c_str(), var1.c_str(), var2, var3, var4);
      vars.push_back(line);
 
    }
  }

 varlist.close();

  for (unsigned int d = 0; d < datasets.size(); d++){
  samp_name =  datasets[d]->Filenames();

  dataset_name =  datasets[d]->Name();

  cout <<"datset name "<< dataset_name << endl; 

  TFile *f = new TFile(samp_name[0].c_str());
  TNtuple * tup  = (TNtuple*)f->Get("Craneen__MuEl");
  int nEvents = tup->GetEntries();
  Float_t HT, NormFactor, GenWeight, HT2M, nMtags,LeadingMuonPt, LeadingElectronPt, LeadingBJetPt, PU, ScaleFactor, nJets;

    tup->SetBranchAddress("NormFactor",&NormFactor);
    tup->SetBranchAddress("GenWeight",&GenWeight);
    tup->SetBranchAddress("ScaleFactor",&ScaleFactor);
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
    double final_weight = Luminosity;
    // double final_weight = 100.;

    if ( write_sys ){

 if (dataset_name == "TTJets"){
    histos1D["PU"]->Fill(PU);
    histos1D["HT"]->Fill(HT);
    histos1D["HT2M"]->Fill(HT2M);
    histos1D["nMtags"]->Fill(nMtags);
    histos1D["nJets"]->Fill(nJets);
    histos1D["LeadingMuonPt"]->Fill(LeadingMuonPt);
    histos1D["LeadingElectronPt"]->Fill(LeadingElectronPt);
    histos1D["LeadingBJetPt"]->Fill(LeadingBJetPt);
  
    cout <<"filling ttjets... "<<  endl; 
 }

else    if (dataset_name == "Plus"){

    histos1D_Up["PU"]->Fill(PU);
    histos1D_Up["HT"]->Fill(HT);
    histos1D_Up["HT2M"]->Fill(HT2M);
    histos1D_Up["nMtags"]->Fill(nMtags);
    histos1D_Up["nJets"]->Fill(nJets);
    histos1D_Up["LeadingMuonPt"]->Fill(LeadingMuonPt);
    histos1D_Up["LeadingElectronPt"]->Fill(LeadingElectronPt);
    histos1D_Up["LeadingBJetPt"]->Fill(LeadingBJetPt);

    }else if(dataset_name == "Minus"){

    histos1D_Down["PU"]->Fill(PU );
    histos1D_Down["HT"]->Fill(HT );
    histos1D_Down["HT2M"]->Fill(HT2M);
    histos1D_Down["nMtags"]->Fill(nMtags );
    histos1D_Down["nJets"]->Fill(nJets );
    histos1D_Down["LeadingMuonPt"]->Fill(LeadingMuonPt);
    histos1D_Down["LeadingElectronPt"]->Fill(LeadingElectronPt );
    histos1D_Down["LeadingBJetPt"]->Fill(LeadingBJetPt);

    } 

    }
    if(dataset_name != "Plus" && dataset_name != "Minus" ){
    MSPlot["PU"]->Fill(PU, datasets[d], true, final_weight );
    MSPlot["HT"]->Fill(HT, datasets[d], true, final_weight );
    MSPlot["HT2M"]->Fill(HT2M, datasets[d], true, final_weight );
    MSPlot["nMtags"]->Fill(nMtags, datasets[d], true, final_weight );
    MSPlot["nJets"]->Fill(nJets, datasets[d], true, final_weight );
    MSPlot["LeadingMuonPt"]->Fill(LeadingMuonPt, datasets[d], true, final_weight );
    MSPlot["LeadingElectronPt"]->Fill(LeadingElectronPt, datasets[d], true, final_weight );
    MSPlot["LeadingBJetPt"]->Fill(LeadingBJetPt, datasets[d], true, final_weight );
    
}
 
}
  }

 f1 = new TFile("StackCraneens.root", "RECREATE");

 cout <<" made file... "<<  endl; 


 if(!write_sys)
 {
 for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
{
 name = it->first;
 MultiSamplePlot *temp = it->second;
 cout <<" before draw ... "<<  endl; 
 temp->Draw("CMSPlot", 2, true, true, true, 1);

 cout <<" after draw ... "<<  endl; 

 void Write(TFile* file, string label = string(""), bool savePNG = false, string pathPNG = string(""), string ext = string("png"));
 temp->Write(f1, name, true, ".", "pdf");
}

 }
 cout <<" msplost written ... "<<  endl; 


 if(write_sys){
 f2 = new TFile("systematics.root", "RECREATE");

for(map<string,TH1F*>::const_iterator i = histos1D.begin(); i !=  histos1D.end(); i++)
{
 cout <<" syshis loop ... "<<  endl; 

 name = i->first;
    cout <<" writing histos... "<<  endl; 
 dir_name = "MultiSamplePlot_" + name;
 f2->mkdir(dir_name.c_str());
 f2->cd(dir_name.c_str());
 hnom = i->second;
 nom_norm = 20366.0;
 double nom_int = hnom->Integral();
 hnom->Scale(nom_norm/nom_int);

 //TH1F* hup = hnom->Clone();
 //TH1F* hdown = hnom->Clone();
 // hup->Scale(1.5);
 //hdown->Scale(0.5);
  hnom->SetName("Nominal");
 hnom->Write("Nominal");
 //hdown->SetName("Minus");
 //hdown->Write("Minus");
 //hup->SetName("Plus");
 // hup->Write("Plus");
 }



 for(map<string,TH1F*>::const_iterator iu = histos1D_Up.begin(); iu !=  histos1D_Up.end(); iu++)
{

 name = iu->first;
 dir_name = "MultiSamplePlot_" + name;
 f2->cd(dir_name.c_str());
 // TH1F *hup = iu->second;
 TH1F* hup = hnom->Clone();
 scalefac = nom_norm/hup->Integral();
 hup->Scale(1.5* scalefac);
 hup->SetName("Plus");
 hup->Write("Plus");

}

for(map<string,TH1F*>::const_iterator id = histos1D_Down.begin(); id !=  histos1D_Down.end(); id++)
{
 name = id->first;
 dir_name = "MultiSamplePlot_" + name;
 f2->cd(dir_name.c_str());
 // TH1F *hdown = id->second;
 TH1F* hdown = hnom->Clone();
 scalefac = nom_norm/hdown->Integral();
 hdown->Scale(0.5*scalefac);
 hdown->SetName("Minus");
 hdown->Write("Minus");
 }



 }


}
