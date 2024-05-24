void BeautifyPlots(TH1F *hist, int FillColor)
{
 hist->SetFillColorAlpha(FillColor, 0.8); 
}

void MakePlots(TString Filename, TString Treename, TString outFileName, bool Writeleg = false, styleMacro = "")
{

	gROOT->ProcessLine(".x " + );
	gStyle->SetBarWidth(0.45);
	gStyle->SetOptFit(111);
	TFile* _file0 = TFile::Open(Filename); 
	TTree* OmegaPi0Eff_FlatTree = (TTree*)_file0->Get(Treename);	
	
	TLegend *leg = new TLegend(0.75, 0.7, 0.9, 0.95);
	leg->SetName("CutNames");
	
	TFile *outFile = TFile::Open(outFileName.Data(), "RECREATE");

	TCanvas *c1 = new TCanvas("RecoilOmega", "RecoilOmega", 1200, 800);
	c1->cd();
	//c1->SetLogy();
	
	// Graph of recoil omega
	
	//No cuts
	OmegaPi0Eff_FlatTree->Draw("omega_14_mass>>Omega14Mass(1000, 0.2, 1.7)", "", "HIST");
	TH1F *Omega14Mass = (TH1F*)gDirectory->Get("Omega14Mass");
	Omega14Mass->SetTitle("Recoil Omega Mass No Cuts;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
	Omega14Mass->GetXaxis()->SetTitleColor(kBlack);
	Omega14Mass->SetFillColorAlpha(kGray, 0.9);
	leg->AddEntry(Omega14Mass, "No Cuts", "f");
	Omega14Mass->Write();
	
	//pi0 cuts
	OmegaPi0Eff_FlatTree->Draw("omega_14_mass>>Omega14MassCut1(1000, 0.2, 1.7)", "(pi0_14_passed) & (pi0_23_passed)", "SAME HIST");
	TH1F *Omega14MassCut1 = (TH1F*)gDirectory->Get("Omega14MassCut1");
	Omega14MassCut1->SetTitle("Recoil Omega Mass No Cuts;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
	Omega14MassCut1->GetXaxis()->SetTitleColor(kBlack);
	Omega14MassCut1->SetFillColorAlpha(kViolet, 0.8);
	leg->AddEntry(Omega14MassCut1, "M(#gamma_{1/2}#gamma_{4/3}) cut", "f");
	Omega14MassCut1->Write();
	// antipi0 12 cut
	OmegaPi0Eff_FlatTree->Draw("omega_14_mass>>Omega14MassCut3(1000, 0.2, 1.7)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed)", "SAME HIST");
	TH1F *Omega14MassCut3 = (TH1F*)gDirectory->Get("Omega14MassCut3");
	Omega14MassCut3->SetTitle("Recoil Omega Mass No Cuts;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
	Omega14MassCut3->GetXaxis()->SetTitleColor(kBlack);
	Omega14MassCut3->SetFillColorAlpha(kBlue, 0.8);
	leg->AddEntry(Omega14MassCut3, "Anti M(#gamma_{1}#gamma_{2}) cut", "f");	
	Omega14MassCut3->Write();
	leg->Draw("SAME");
	
	// antipi0 13 cut
	OmegaPi0Eff_FlatTree->Draw("omega_14_mass>>Omega14MassCut4(1000, 0.2, 1.7)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed) & (antipi0_13_passed)", "SAME HIST");
	TH1F *Omega14MassCut4 = (TH1F*)gDirectory->Get("Omega14MassCut4");
	Omega14MassCut4->SetTitle("Recoil Omega Mass No Cuts;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
	Omega14MassCut4->GetXaxis()->SetTitleColor(kBlack);
	Omega14MassCut4->SetFillColorAlpha(kYellow, 0.8);
	leg->AddEntry(Omega14MassCut4, "Anti M(#gamma_{1}#gamma_{3}) cut", "f");
	Omega14MassCut4->Write();
	// antipi0 24 cut
	OmegaPi0Eff_FlatTree->Draw("omega_14_mass>>Omega14MassCut5(1000, 0.2, 1.7)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed) & (antipi0_13_passed) & (antipi0_24_passed)", "SAME HIST");
	TH1F *Omega14MassCut5 = (TH1F*)gDirectory->Get("Omega14MassCut5");
	Omega14MassCut5->SetTitle("Recoil Omega Mass No Cuts;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
	Omega14MassCut5->GetXaxis()->SetTitleColor(kBlack);
	Omega14MassCut5->SetFillColorAlpha(kCyan, 0.8);
	leg->AddEntry(Omega14MassCut5, "Anti M(#gamma_{2}#gamma_{4}) cut", "f");	
	Omega14MassCut5->Write();
	// antipi0 34 cut
	OmegaPi0Eff_FlatTree->Draw("omega_14_mass>>Omega14MassCut6(1000, 0.2, 1.7)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed) & (antipi0_13_passed) & (antipi0_24_passed) & (antipi0_34_passed)", "SAME HIST");
	TH1F *Omega14MassCut6 = (TH1F*)gDirectory->Get("Omega14MassCut6");
	Omega14MassCut6->SetTitle("Recoil Omega Mass No Cuts;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
	Omega14MassCut6->GetXaxis()->SetTitleColor(kBlack);
	Omega14MassCut6->SetFillColorAlpha(kMagenta, 0.8);
	leg->AddEntry(Omega14MassCut6, "Anti M(#gamma_{3}#gamma_{4}) cut", "f");
	Omega14MassCut6->Write();
	// antiomega 23 cut
	OmegaPi0Eff_FlatTree->Draw("omega_14_mass>>Omega14MassCut2(1000, 0.2, 1.7)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_13_passed) & (antipi0_24_passed) & (antipi0_34_passed)  & (antipi0_12_passed) & (antiomega_23_passed)", "SAME HIST");
	TH1F *Omega14MassCut2 = (TH1F*)gDirectory->Get("Omega14MassCut2");
	Omega14MassCut2->SetTitle("Recoil Omega Mass No Cuts;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
	Omega14MassCut2->GetXaxis()->SetTitleColor(kBlack);
	Omega14MassCut2->SetFillColorAlpha(kGreen, 0.8);
	leg->AddEntry(Omega14MassCut2, "Anti #omega(#gamma_{2}#gamma_{3}) cut", "f");
	Omega14MassCut2->Write();
	
	// antiomega 23 cut
	OmegaPi0Eff_FlatTree->Draw("omega_14_mass>>Omega14MassCut7(1000, 0.2, 1.7)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_13_passed) & (antipi0_24_passed) & (antipi0_34_passed)  & (antipi0_12_passed) & (antiomega_23_passed) & (P4_g4.Theta()*TMath::RadToDeg() > 10.)", "SAME HIST");
	TH1F *Omega14MassCut7 = (TH1F*)gDirectory->Get("Omega14MassCut7");
	Omega14MassCut7->SetTitle("Recoil Omega Mass No Cuts;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
	Omega14MassCut7->GetXaxis()->SetTitleColor(kBlack);
	Omega14MassCut7->SetFillColorAlpha(kBlack, 0.8);
	leg->AddEntry(Omega14MassCut7, "#gamma_{4} in BCAL cut", "f");
	Omega14MassCut7->Write();
	leg->Draw("SAME");
	if(Writeleg) leg->Write();


	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4Theta(180, 0, 180.)", "", "HIST");
	TH1F *g4Theta = (TH1F*)gDirectory->Get("g4Theta");
	g4Theta->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4Theta->GetXaxis()->SetTitleColor(kBlack);
	g4Theta->SetFillColorAlpha(kGray, 0.9);
	g4Theta->Write();
	
	//pi0 cuts
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut1(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed)", "SAME HIST");
	TH1F *g4ThetaCut1 = (TH1F*)gDirectory->Get("g4ThetaCut1");
	g4ThetaCut1->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut1->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut1->SetFillColorAlpha(kViolet, 0.8);
	g4ThetaCut1->Write();
	
	// antipi0 12 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut3(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed)", "SAME HIST");
	TH1F *g4ThetaCut3 = (TH1F*)gDirectory->Get("g4ThetaCut3");
	g4ThetaCut3->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut3->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut3->SetFillColorAlpha(kBlue, 0.8);
	g4ThetaCut3->Write();

	
	// antipi0 13 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut4(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed) & (antipi0_13_passed)", "SAME HIST");
	TH1F *g4ThetaCut4 = (TH1F*)gDirectory->Get("g4ThetaCut4");
	g4ThetaCut4->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut4->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut4->SetFillColorAlpha(kYellow, 0.8);
	g4ThetaCut4->Write();
	
	// antipi0 24 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut5(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed) & (antipi0_13_passed) & (antipi0_24_passed)", "SAME HIST");
	TH1F *g4ThetaCut5 = (TH1F*)gDirectory->Get("g4ThetaCut5");
	g4ThetaCut5->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut5->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut5->SetFillColorAlpha(kCyan, 0.8);
	
	
	// antipi0 34 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut6(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed) & (antipi0_13_passed) & (antipi0_24_passed) & (antipi0_34_passed)", "SAME HIST");
	TH1F *g4ThetaCut6 = (TH1F*)gDirectory->Get("g4ThetaCut6");
	g4ThetaCut6->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut6->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut6->SetFillColorAlpha(kMagenta, 0.8);
	g4ThetaCut6->Write();
	
	// antiomega 23 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut2(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_13_passed) & (antipi0_24_passed) & (antipi0_34_passed)  & (antipi0_12_passed) & (antiomega_23_passed)", "SAME HIST");
	TH1F *g4ThetaCut2 = (TH1F*)gDirectory->Get("g4ThetaCut2");
	g4ThetaCut2->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut2->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut2->SetFillColorAlpha(kGreen, 0.8);
	g4ThetaCut2->Write();

	//c1->Print("Cuts.pdf");
	TF1* fitfunc = new TF1("doublegaus", "gaus(0) + gaus(3) + gaus(6)", 0.6, 1.2);
	TF1* gaus1 = new TF1("gaus1", "gaus", 0.6, 1.2);
	gaus1->SetLineColor(kOrange);
	TF1* gaus2 = new TF1("gaus2", "gaus", 0.6, 1.2);
	gaus2->SetLineColor(kCyan);
	TF1* gaus3 = new TF1("gaus3", "gaus", 0.6, 1.2);
	gaus3->SetLineColor(kYellow);
	fitfunc->SetLineWidth(3);
	fitfunc->SetLineColor(kGreen);
	fitfunc->SetParameters(Omega14MassCut7->GetMaximum()*0.8, 0.782, 0.008, Omega14MassCut7->GetMaximum()*0.1, 0.782, 0.10, Omega14MassCut7->GetMaximum()*0.1, 0.782, 0.2);
	fitfunc->SetParLimits(2, 0., 0.1);
	fitfunc->SetParLimits(5, 0., 1.);
	fitfunc->SetParLimits(8, 0., 1.);
	TCanvas *C3 = new TCanvas("FitFunc", "FitFunc", 1200, 800);
	C3->cd();
	//Omega14MassCut7->SetFillColorAlpha(kWhite, 1.);
	Omega14MassCut7->Fit(fitfunc, "R");
	gaus1->SetParameters(fitfunc->GetParameter(0), fitfunc->GetParameter(1), fitfunc->GetParameter(2));
	gaus2->SetParameters(fitfunc->GetParameter(3), fitfunc->GetParameter(4), fitfunc->GetParameter(5));
	gaus3->SetParameters(fitfunc->GetParameter(6), fitfunc->GetParameter(7), fitfunc->GetParameter(8));
	Omega14MassCut7->Draw("HIST");
	fitfunc->Draw("SAME");
	gaus1->Draw("SAME");
	gaus2->Draw("SAME");
	gaus3->Draw("SAME");
	
	cout << fitfunc->GetChisquare() << " " << fitfunc->GetNDF() << endl;
	
	//C3->Print("Fits.pdf");
	
	const int NPoints = 7;
	TH1F *dummy = new TH1F("dummy", "", NPoints, 0, NPoints);
	
	TCanvas *c2 = new TCanvas("StatReduction", "StatReduction", 1200, 800);
	c2->cd();
	//c2->SetLogy();
	
	TString Labels[NPoints] = {"No Cuts", "M(#gamma_{1/2}#gamma_{4/3}) cut",
								"Anti M(#gamma_{1}#gamma_{2}) cut", "Anti M(#gamma_{1}#gamma_{3}) cut",
								"Anti M(#gamma_{2}#gamma_{4}) cut", "Anti M(#gamma_{3}#gamma_{4}) cut",
								"Anti #omega(#gamma_{2}#gamma_{3}) cut"};
	
	for(unsigned int i = 0; i < NPoints; i++)
	{
		dummy->Fill(i + 0.5, Omega14Mass->GetEntries());
		dummy->GetXaxis()->SetBinLabel(i+1, Labels[i]);
	}
	
	dummy->SetLineColor(kWhite);
	dummy->SetLineWidth(0);
	dummy->SetMarkerSize(0);
	dummy->SetFillColorAlpha(kWhite, 0);
	dummy->SetMaximum(Omega14Mass->GetEntries()*1.0);
	dummy->SetMinimum(0);
	dummy->SetBarWidth(0.45);
	
	TH1F *gr = new TH1F("NoCut", "NoCut", NPoints, 0, NPoints);
	gr->Fill(0.5, Omega14Mass->GetEntries());
	gr->SetFillColor(kGray);
	gr->SetLineColor(kBlack);
	gr->SetBarWidth(0.45);
	
	TH1F *gr1 = new TH1F("Cut1", "Cut1", NPoints, 0, NPoints);
	gr1->Fill(1.5, Omega14MassCut1->GetEntries());
	gr1->SetFillColor(kViolet);
	gr1->SetLineColor(kBlack);
	gr1->SetBarWidth(0.45);
	
	TH1F *gr2 = new TH1F("Cut3", "Cut3", NPoints, 0, NPoints);
	gr2->Fill(2.5, Omega14MassCut3->GetEntries());
	gr2->SetFillColor(kBlue);
	gr2->SetLineColor(kBlack);
	gr2->SetBarWidth(0.45);

	TH1F *gr3 = new TH1F("Cut4", "Cut4", NPoints, 0, NPoints);
	gr3->Fill(3.5, Omega14MassCut4->GetEntries());
	gr3->SetFillColor(kYellow);
	gr3->SetLineColor(kBlack);
	gr3->SetBarWidth(0.45);
	
	TH1F *gr4 = new TH1F("Cut5", "Cut5", NPoints, 0, NPoints);
	gr4->Fill(4.5, Omega14MassCut5->GetEntries());
	gr4->SetFillColor(kCyan);
	gr4->SetLineColor(kBlack);
	gr4->SetBarWidth(0.45);
	
	TH1F *gr5 = new TH1F("Cut6", "Cut6", NPoints, 0, NPoints);
	gr5->Fill(5.5, Omega14MassCut6->GetEntries());
	gr5->SetFillColor(kMagenta);
	gr5->SetLineColor(kBlack);
	gr5->SetBarWidth(0.45);
	
	TH1F *gr6 = new TH1F("Cut2", "Cut2", NPoints, 0, NPoints);
	gr6->Fill(6.5, Omega14MassCut2->GetEntries());
	gr6->SetFillColor(kGreen);
	gr6->SetLineColor(kBlack);
	gr6->SetBarWidth(0.45);
	
	
	dummy->Draw("HIST");
	gr->Draw("SAME HIST");
	gr1->Draw("SAME HIST");
	gr2->Draw("SAME HIST");
	gr3->Draw("SAME HIST");
	gr4->Draw("SAME HIST");
	gr5->Draw("SAME HIST");
	gr6->Draw("SAME HIST");
	//c2->Print("CUTS.pdf");

	dummy->Write();
	gr->Write();
	gr2->Write();
	gr3->Write();
	gr4->Write();
	gr5->Write();
	gr6->Write();
	outFile->Close();
	
	
	
}
