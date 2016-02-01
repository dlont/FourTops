{
	gStyle->SetOptStat(0);
	TFile *_file0 = TFile::Open("shapefileMu_DATA_BDT_inc.root");
	TH1F* nominal = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__nominal");

	TH1F* PUUp = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__PUUP");
	TH1F* PUDown = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__PUDown");
	TH1F* PUUpOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__PUUP");
	TH1F* PUDownOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__PUDown");

	TH1F* btagDown = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__btagDown");
	TH1F* btagUp = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__btagUp");
	TH1F* btagDownOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__btagDown");
	TH1F* btagUpOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__btagUp");

	TH1F* scaleDown = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__scaleDown");
	TH1F* scaleUp = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__scaleUp");
	TH1F* scaleDownOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__scaleDown");
	TH1F* scaleUpOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__scaleUp");

	TH1F* JESUp = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__JESUp");
	TH1F* JESDown = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__JESDown");
	TH1F* JESUpOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__JESUp");
	TH1F* JESDownOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__JESDown");

	/*
	TH1F* JERUp = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__JERUp");
	TH1F* JERDown = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__JERDown");
	TH1F* JERUpOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__JERUp");
	TH1F* JERDownOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__JERDown"); 				

	TH1F* MatchingUp = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__MatchingUp");
	TH1F* MatchingDown = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__MatchingDown");
	TH1F* MatchingUpOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__MatchingUp");
	TH1F* MatchingDownOrig = (TH1F*)_file0->Get("ttttmu__TTJets_MLM__MatchingDown"); 	
	*/

 	//btag

	TCanvas *c1 = new TCanvas();
	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	c1->cd();
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();

	btagUp->SetLineColor(kRed);
	btagDown->SetLineColor(kCyan);

	btagUpOrig->SetLineColor(kRed);
	btagDownOrig->SetLineColor(kCyan);
	btagUpOrig->SetTitle("TTbar BDT disciminator");
	btagUpOrig->Draw();
	btagDownOrig->Draw("same");
	nominal->Draw("same");

	TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
	leg->AddEntry(btagUp,"btag Up","l");
	leg->AddEntry(nominal,"nominal ttbar","l");
	leg->AddEntry(btagDown,"btag Down","l");
	leg->Draw();
	c1->cd();

	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);

	pad2->SetTopMargin(0);
	pad2->SetGridy();
	pad2->Draw();

	btagUp->SetTitle("");

	pad2->cd();
	// btagUp->Sumw2();
	btagUp->SetStats(0);
	btagUp->Divide(nominal);
	btagUp->SetMarkerStyle(21);
	btagUp->Draw("ep");
	// btagDown>Sumw2();
	btagDown->SetStats(0);
	btagDown->Divide(nominal);
	btagDown->SetMarkerStyle(21);
	btagDown->Draw("epsame");
	c1->cd();
	c1->SaveAs("btagsys.pdf");


 	// PU
	TCanvas *cPU1 = new TCanvas();
	TPad *padPU1 = new TPad("pad1","pad1",0,0.3,1,1);
	cPU1->cd();
	padPU1->SetBottomMargin(0);
	padPU1->Draw();
	padPU1->cd();

	PUUp->SetLineColor(kRed);
	PUDown->SetLineColor(kCyan);

	PUUpOrig->SetLineColor(kRed);
	PUDownOrig->SetLineColor(kCyan);
	PUDownOrig->SetTitle("TTbar BDT disciminator");
	PUDownOrig->Draw();
	PUUpOrig->Draw("same");
	nominal->Draw("same");

	TLegend *legPU = new TLegend(0.5,0.7,0.9,0.9);
	legPU->AddEntry(PUDown,"PU Up","l");
	legPU->AddEntry(nominal,"nominal ttbar","l");
	legPU->AddEntry(PUUp,"PU Down","l");

	legPU->Draw();
	cPU1->cd();

	TPad *padPU2 = new TPad("pad2","pad2",0,0,1,0.3);

	padPU2->SetTopMargin(0);
	padPU2->SetGridy();
	padPU2->Draw();

	PUUp->SetTitle("");

	padPU2->cd();
	// PUUp->Sumw2();
	PUUp->SetStats(0);
	PUUp->Divide(nominal);
	PUUp->SetMarkerStyle(21);
	PUUp->Draw("ep");
	// PUDown>Sumw2();
	PUDown->SetStats(0);
	PUDown->Divide(nominal);
	PUDown->SetMarkerStyle(21);
	PUDown->Draw("epsame");
	cPU1->cd();
	cPU1->SaveAs("PUsys.pdf");

	//scale

	TCanvas *cscale1 = new TCanvas();
	TPad *padscale1 = new TPad("pad1","pad1",0,0.3,1,1);
	cscale1->cd();
	padscale1->SetBottomMargin(0);
	padscale1->Draw();
	padscale1->cd();

	scaleUp->SetLineColor(kRed);
	scaleDown->SetLineColor(kCyan);

	scaleUpOrig->SetLineColor(kRed);
	scaleDownOrig->SetLineColor(kCyan);
	scaleUpOrig->SetTitle("TTbar BDT disciminator");
	scaleUpOrig->Draw();
	scaleDownOrig->Draw("same");
	nominal->Draw("same");

	TLegend *legscale = new TLegend(0.5,0.7,0.9,0.9);
	legscale->AddEntry(scaleUp,"scale Up","l");
	legscale->AddEntry(nominal,"nominal ttbar","l");
	legscale->AddEntry(scaleDown,"scale Down","l");
	legscale->Draw();
	cscale1->cd();

	TPad *pad2scale = new TPad("pad2","pad2",0,0,1,0.3);

	pad2scale->SetTopMargin(0);
	pad2scale->SetGridy();
	pad2scale->Draw();

	scaleUp->SetTitle("");

	pad2scale->cd();
	// scaleUp->Sumw2();
	scaleUp->SetStats(0);
	scaleUp->Divide(nominal);
	scaleUp->SetMarkerStyle(21);
	scaleUp->Draw("ep");
	// scaleDown>Sumw2();
	scaleDown->SetStats(0);
	scaleDown->Divide(nominal);
	scaleDown->SetMarkerStyle(21);
	scaleDown->Draw("epsame");
	cscale1->cd();
	cscale1->SaveAs("scalesys.pdf");

	//JES
	TCanvas *cJES1 = new TCanvas();
	TPad *padJES1 = new TPad("pad1","pad1",0,0.3,1,1);
	cJES1->cd();
	padJES1->SetBottomMargin(0);
	padJES1->Draw();
	padJES1->cd();

	JESUp->SetLineColor(kRed);
	JESDown->SetLineColor(kCyan);

	JESUpOrig->SetLineColor(kRed);
	JESDownOrig->SetLineColor(kCyan);
	JESDownOrig->SetTitle("TTbar BDT disciminator");
	JESDownOrig->Draw();
	JESUpOrig->Draw("same");
	nominal->Draw("same");

	TLegend *legJES = new TLegend(0.5,0.7,0.9,0.9);
	legJES->AddEntry(JESDown,"JES Up","l");
	legJES->AddEntry(nominal,"nominal ttbar","l");
	legJES->AddEntry(JESUp,"JES Down","l");

	legJES->Draw();
	cJES1->cd();

	TPad *padJES2 = new TPad("pad2","pad2",0,0,1,0.3);

	padJES2->SetTopMargin(0);
	padJES2->SetGridy();
	padJES2->Draw();

	JESUp->SetTitle("");

	padJES2->cd();
	// JESUp->Sumw2();
	JESUp->SetStats(0);
	JESUp->Divide(nominal);
	JESUp->SetMarkerStyle(21);
	JESUp->Draw("ep");
	// JESDown>Sumw2();
	JESDown->SetStats(0);
	JESDown->Divide(nominal);
	JESDown->SetMarkerStyle(21);
	JESDown->Draw("epsame");
	cJES1->cd();
	cJES1->SaveAs("JESsys.pdf");

}
