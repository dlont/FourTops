#include "../interface/Zpeak.h"


Zpeak::Zpeak(vector < Dataset* > datasets):
MSPlot(),
looseMuons_(),
tightMuons_(),
looseElectrons_(),
tightElectrons_(),
diLeptonCollection_(),
diLeptonInvMass_(0),
isTwoLeptons_(false)
{
	MSPlot["invarMassLeptons"] = new MultiSamplePlot(datasets, "invarMassLeptons", 40, 0, 200, "Invm_{ll}");
}

Zpeak::~Zpeak(){

}

void Zpeak::invariantMass(Run2Selection r2selection){
	isTwoLeptons_=false;
	looseMuons_ = r2selection.GetSelectedMuons(26, 2.1, 0.12, "Loose", "Spring15");
	tightMuons_ = r2selection.GetSelectedMuons(26, 2.1, 0.12, "Tight", "Spring15");
	looseElectrons_ = r2selection.GetSelectedElectrons("Loose", "Spring15_50ns", true);
	tightElectrons_ = r2selection.GetSelectedElectrons("Tight", "Spring15_50ns", true);
	if (tightMuons_.size() == 2){
		diLeptonCollection_ = *tightMuons_[0] + *tightMuons_[1];
		diLeptonInvMass_ = diLeptonCollection_.M();
		isTwoLeptons_=true;
	}
	else if (tightElectrons_.size() == 2){
		diLeptonCollection_ = *tightElectrons_[0] + *tightElectrons_[1];
		diLeptonInvMass_ = diLeptonCollection_.M();
		isTwoLeptons_=true;
	}
}

void Zpeak::fillPlot(vector < Dataset* > datasets, unsigned int d, float Luminosity, float scaleFactor){
	if (tightMuons_.size() == 2 || tightElectrons_.size() == 2){
		MSPlot["invarMassLeptons"]->Fill(diLeptonInvMass_, datasets[d], true, Luminosity*scaleFactor);
	}
}

void Zpeak::writeErase(TFile *fout, string pathPNG){
	for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        MultiSamplePlot *temp = it->second;
        temp->Write(fout, name, true, pathPNG, "pdf");
		MSPlot.erase(name);
    }
}

float Zpeak::returnInvMass(){
	return diLeptonInvMass_;
}

bool Zpeak::requireTwoLeptons(){
	return isTwoLeptons_;
}