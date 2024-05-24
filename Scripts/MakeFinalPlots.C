void MakeFinalPlots(TString FileName)
{
    gROOT->ProcessLine(".x /work/Thesis/Scripts/gluex_style.C");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    TFile *inFile = TFile::Open(FileName.Data());

    TLegend *leg = new TLegend(0.2, 0.5, 0.4, 0.95);
    leg->SetName("CutNames");

    TH1F *Omega14Mass = (TH1F*)inFile->Get("Omega14Mass_0_0");
    Omega14Mass->SetFillColorAlpha(kGray, 0.9);
    leg->AddEntry(Omega14Mass, "No Cuts", "f");

    TH1F *Omega14MassCut1 = (TH1F*)inFile->Get("Omega14MassCut1_0_0");
 	Omega14MassCut1->SetFillColorAlpha(kViolet, 0.8);
    leg->AddEntry(Omega14MassCut1, "M(#gamma_{1/2}#gamma_{4/3}) cut", "f");
    
    TH1F *Omega14MassCut3 = (TH1F*)inFile->Get("Omega14MassCut3_0_0");
	Omega14MassCut3->SetFillColorAlpha(kBlue, 0.8); 
    leg->AddEntry(Omega14MassCut3, "Anti M(#gamma_{1}#gamma_{2}) cut", "f");
    
    TH1F *Omega14MassCut4 = (TH1F*)inFile->Get("Omega14MassCut4_0_0");
 	Omega14MassCut4->SetFillColorAlpha(kYellow, 0.8);   
    leg->AddEntry(Omega14MassCut4, "Anti M(#gamma_{1}#gamma_{3}) cut", "f");
    
    TH1F *Omega14MassCut5 = (TH1F*)inFile->Get("Omega14MassCut5_0_0");
	Omega14MassCut5->SetFillColorAlpha(kCyan, 0.8);	
    leg->AddEntry(Omega14MassCut5, "Anti M(#gamma_{2}#gamma_{4}) cut", "f");
    
    TH1F *Omega14MassCut6 = (TH1F*)inFile->Get("Omega14MassCut6_0_0");
	Omega14MassCut6->SetFillColorAlpha(kMagenta, 0.8);
	leg->AddEntry(Omega14MassCut6, "Anti M(#gamma_{3}#gamma_{4}) cut", "f");
    
    TH1F *Omega14MassCut2 = (TH1F*)inFile->Get("Omega14MassCut2_0_0");
	Omega14MassCut2->SetFillColorAlpha(kGreen, 0.8);	
	leg->AddEntry(Omega14MassCut2, "Anti #omega(#gamma_{2}#gamma_{3}) cut", "f");
    
    TH1F *Omega14MassCut7 = (TH1F*)inFile->Get("Omega14MassCut7_0_0");
	Omega14MassCut7->SetFillColorAlpha(kBlack, 0.8);
	leg->AddEntry(Omega14MassCut7, "#gamma_{4} in BCAL cut", "f");

    TH1F *Omega14InvMassCut7 = (TH1F*)inFile->Get("Omega14InvariantMassCut7_0_0");
	Omega14InvMassCut7->SetFillColorAlpha(kBlack, 0.8);
	//leg->AddEntry(Omega14InvMassCut7, "#gamma_{4} in BCAL cut", "f");

    TCanvas *C1 = new TCanvas("RecoilOmega", "RecoilOmega", 1200, 800);
    C1->cd();
    Omega14Mass->Draw("HIST");
    Omega14MassCut1->Draw("SAME HIST");
    Omega14MassCut3->Draw("SAME HIST");
    Omega14MassCut4->Draw("SAME HIST");
    Omega14MassCut5->Draw("SAME HIST");
    Omega14MassCut6->Draw("SAME HIST");
    Omega14MassCut2->Draw("SAME HIST");
    Omega14MassCut7->Draw("SAME HIST");
    leg->Draw("SAME");
    C1->Print("Recoid3PiMass_Data_Part.pdf");


    TF1* fitfunc = new TF1("doublegaus", "gaus(0) + gaus(3) + pol2(6)", 0.6, 1.2);
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
	TCanvas *C2 = new TCanvas("FitFunc", "FitFunc", 2400, 800);
    C2->Divide(2, 1);
	C2->cd(1);
	//Omega14MassCut7->SetFillColorAlpha(kWhite, 1.);
	Omega14MassCut7->Fit(fitfunc, "R0");
	gaus1->SetParameters(fitfunc->GetParameter(0), fitfunc->GetParameter(1), fitfunc->GetParameter(2));
	gaus2->SetParameters(fitfunc->GetParameter(3), fitfunc->GetParameter(4), fitfunc->GetParameter(5));
	gaus3->SetParameters(fitfunc->GetParameter(6), fitfunc->GetParameter(7), fitfunc->GetParameter(8));
	Omega14MassCut7->Draw("HIST");
    Omega14MassCut7->GetXaxis()->SetTitleSize(0.05);
    C2->cd(2);
    Omega14InvMassCut7->Draw("HIST");
    Omega14InvMassCut7->GetXaxis()->SetTitleSize(0.05);
	//fitfunc->Draw("SAME");
	//gaus1->Draw("SAME");
	//gaus2->Draw("SAME");
	//gaus3->Draw("SAME");
	
	cout << fitfunc->GetChisquare() << " " << fitfunc->GetNDF() << endl;
	
	//C2->Print("Recoil3PiMassFit.pdf");
    C2->Print("3PiMassFit.pdf");
    /*
    TLegend *leg1 = new TLegend(0.75, 0.7, 0.9, 0.95);
    leg1->SetName("CutNames");
    
    TH1F *g4Theta = (TH1F*)inFile->Get("g4Theta");
    g4Theta->SetFillColorAlpha(kGray, 0.9);
    leg1->AddEntry(g4Theta, "No Cuts", "f");

	TH1F *g4ThetaCut1 = (TH1F*)inFile->Get("g4ThetaCut1");
	g4ThetaCut1->SetFillColorAlpha(kViolet, 0.8);
    leg1->AddEntry(g4ThetaCut1, "M(#gamma_{1/2}#gamma_{4/3}) cut", "f");

	TH1F *g4ThetaCut3 = (TH1F*)inFile->Get("g4ThetaCut3");
	g4ThetaCut3->SetFillColorAlpha(kBlue, 0.8);
    leg1->AddEntry(g4ThetaCut3, "Anti M(#gamma_{1}#gamma_{2}) cut", "f");

	TH1F *g4ThetaCut4 = (TH1F*)inFile->Get("g4ThetaCut4");
	g4ThetaCut4->SetFillColorAlpha(kYellow, 0.8);
    leg1->AddEntry(g4ThetaCut4, "Anti M(#gamma_{1}#gamma_{3}) cut", "f");

	
    TH1F *g4ThetaCut5 = (TH1F*)inFile->Get("g4ThetaCut5");
	g4ThetaCut5->SetFillColorAlpha(kCyan, 0.8);
    leg1->AddEntry(g4ThetaCut5, "Anti M(#gamma_{2}#gamma_{4}) cut", "f");
    

	TH1F *g4ThetaCut6 = (TH1F*)inFile->Get("g4ThetaCut6");
	g4ThetaCut6->SetFillColorAlpha(kMagenta, 0.8);
    leg1->AddEntry(g4ThetaCut6, "Anti M(#gamma_{3}#gamma_{4}) cut", "f");

	TH1F *g4ThetaCut2 = (TH1F*)inFile->Get("g4ThetaCut2");
	g4ThetaCut2->SetFillColorAlpha(kGreen, 0.8);
    leg1->AddEntry(g4ThetaCut2, "Anti #omega(#gamma_{2}#gamma_{3}) cut", "f");


    TCanvas *C3 = new TCanvas("thetaCanvas", "thetaCanvas", 1200, 800);
    C3->cd();
    C3->SetLogy();

    g4Theta->Draw("HIST");
    g4ThetaCut1->Draw("SAME HIST");
    g4ThetaCut3->Draw("SAME HIST");
    g4ThetaCut4->Draw("SAME HIST");
    //g4ThetaCut5->Draw("SAME HIST");
    g4ThetaCut6->Draw("SAME HIST");
    g4ThetaCut2->Draw("SAME HIST");
    leg1->Draw("SAME");

    C3->Print("ThetaDistribution_Data.pdf");
    */
    

}