void BeautifyPlots(TH1F *hist, int FillColor)
{
 hist->SetFillColorAlpha(FillColor, 0.8); 
}

void MakeThetaPlots(TString Filename, TString Treename)
{

	gROOT->ProcessLine(".x $ROOTSYS/macros/gluex_style.C");
	gStyle->SetBarWidth(0.45);
	TFile* _file0 = TFile::Open(Filename); 
	TTree* OmegaPi0Eff_FlatTree = (TTree*)_file0->Get(Treename);	
	
	TLegend *leg = new TLegend(0.75, 0.7, 0.9, 0.95);
	
	TCanvas *c1 = new TCanvas("RecoilOmega", "RecoilOmega", 1200, 800);
	c1->cd();
	c1->SetLogy();
	
	// Graph of recoil omega
	
	//No cuts
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4Theta(180, 0, 180.)", "", "HIST");
	TH1F *g4Theta = (TH1F*)gDirectory->Get("g4Theta");
	g4Theta->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4Theta->GetXaxis()->SetTitleColor(kBlack);
	g4Theta->SetFillColorAlpha(kGray, 0.9);
	leg->AddEntry(g4Theta, "No Cuts", "f");
	
	//pi0 cuts
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut1(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed)", "SAME HIST");
	TH1F *g4ThetaCut1 = (TH1F*)gDirectory->Get("g4ThetaCut1");
	g4ThetaCut1->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut1->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut1->SetFillColorAlpha(kViolet, 0.8);
	leg->AddEntry(g4ThetaCut1, "M(#gamma_{1/2}#gamma_{4/3}) cut", "f");
	
	// antipi0 12 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut3(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed)", "SAME HIST");
	TH1F *g4ThetaCut3 = (TH1F*)gDirectory->Get("g4ThetaCut3");
	g4ThetaCut3->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut3->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut3->SetFillColorAlpha(kBlue, 0.8);
	leg->AddEntry(g4ThetaCut3, "Anti M(#gamma_{1}#gamma_{2}) cut", "f");	
	
	leg->Draw("SAME");
	
	// antipi0 13 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut4(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed) & (antipi0_13_passed)", "SAME HIST");
	TH1F *g4ThetaCut4 = (TH1F*)gDirectory->Get("g4ThetaCut4");
	g4ThetaCut4->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut4->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut4->SetFillColorAlpha(kYellow, 0.8);
	leg->AddEntry(g4ThetaCut4, "Anti M(#gamma_{1}#gamma_{3}) cut", "f");
	
	// antipi0 24 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut5(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed) & (antipi0_13_passed) & (antipi0_24_passed)", "SAME HIST");
	TH1F *g4ThetaCut5 = (TH1F*)gDirectory->Get("g4ThetaCut5");
	g4ThetaCut5->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut5->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut5->SetFillColorAlpha(kCyan, 0.8);
	leg->AddEntry(g4ThetaCut5, "Anti M(#gamma_{2}#gamma_{4}) cut", "f");	
	
	// antipi0 34 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut6(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_12_passed) & (antipi0_13_passed) & (antipi0_24_passed) & (antipi0_34_passed)", "SAME HIST");
	TH1F *g4ThetaCut6 = (TH1F*)gDirectory->Get("g4ThetaCut6");
	g4ThetaCut6->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut6->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut6->SetFillColorAlpha(kMagenta, 0.8);
	leg->AddEntry(g4ThetaCut6, "Anti M(#gamma_{3}#gamma_{4}) cut", "f");
	
	// antiomega 23 cut
	OmegaPi0Eff_FlatTree->Draw("P4_g4.Theta()*TMath::RadToDeg()>>g4ThetaCut2(180, 0, 180.)", "(pi0_14_passed) & (pi0_23_passed) & (antipi0_13_passed) & (antipi0_24_passed) & (antipi0_34_passed)  & (antipi0_12_passed) & (antiomega_23_passed)", "SAME HIST");
	TH1F *g4ThetaCut2 = (TH1F*)gDirectory->Get("g4ThetaCut2");
	g4ThetaCut2->SetTitle(";Missing #gamma_{4} Theta (#theta) [deg] [GeV];Counts / deg");
	g4ThetaCut2->GetXaxis()->SetTitleColor(kBlack);
	g4ThetaCut2->SetFillColorAlpha(kGreen, 0.8);
	leg->AddEntry(g4ThetaCut2, "Anti #omega(#gamma_{2}#gamma_{3}) cut", "f");
	
	leg->Draw("SAME");
	
	
	
}
