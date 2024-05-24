void BeautifyPlots(TH1F *hist, int FillColor)
{
 hist->SetFillColorAlpha(FillColor, 0.8); 
}

void MakePlots(TString Filename, TString Treename, 
				TString outFileName, bool Writeleg = false, 
				TString styleMacro = "/work/Thesis/Scripts/gluex_style.C"
			)
{

	gROOT->ProcessLine(".x " + styleMacro);
	gStyle->SetBarWidth(0.45);
	gStyle->SetOptFit(111);
	TFile* _file0 = TFile::Open(Filename); 
	TTree* OmegaPi0Eff_FlatTree = (TTree*)_file0->Get(Treename);	

	Double_t omega_14_mass = 0., omega_23_mass = 0., pi0_14_mass = 0., pi0_23_mass = 0., pi0_12_mass = 0., pi0_13_mass = 0., pi0_24_mass = 0., pi0_34_mass = 0.;
	bool pi0_14_passed = false, pi0_23_passed = false;
	bool antipi0_12_passed = false, antipi0_13_passed = false, antipi0_24_passed = false, antipi0_34_passed = false;
	bool antiomega_23_passed = false;
	Double_t P4_g4_Theta = 0.;
	Double_t accidweight = 0.;
	TLorentzVector *P4_g4 = nullptr;

	OmegaPi0Eff_FlatTree->SetBranchAddress("accidweight", &accidweight);
	OmegaPi0Eff_FlatTree->SetBranchAddress("omega_14_mass", &omega_14_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("omega_23_mass", &omega_23_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_14_mass", &pi0_14_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_23_mass", &pi0_23_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_12_mass", &pi0_12_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_13_mass", &pi0_13_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_24_mass", &pi0_24_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_34_mass", &pi0_34_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_14_passed", &pi0_14_passed);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_23_passed", &pi0_23_passed);
	OmegaPi0Eff_FlatTree->SetBranchAddress("antipi0_12_passed", &antipi0_12_passed);
	OmegaPi0Eff_FlatTree->SetBranchAddress("antipi0_13_passed", &antipi0_13_passed);
	OmegaPi0Eff_FlatTree->SetBranchAddress("antipi0_24_passed", &antipi0_24_passed);
	OmegaPi0Eff_FlatTree->SetBranchAddress("antipi0_34_passed", &antipi0_34_passed);
	OmegaPi0Eff_FlatTree->SetBranchAddress("antiomega_23_passed", &antiomega_23_passed);
	//OmegaPi0Eff_FlatTree->SetBranchAddress("P4_g4_Theta", &P4_g4_Theta);
	OmegaPi0Eff_FlatTree->SetBranchAddress("P4_g4", &P4_g4);

	int nEntries = OmegaPi0Eff_FlatTree->GetEntries();
	cout << "Number of Entries: " << nEntries << endl;
	const int nThetas = (int)((180. - 0.)/5.);
	const int nMomentum = (int)((3.0 - 0.)/0.50);
	TH1F *Omega14Mass[nMomentum][nThetas];
	TH1F *Omega14MassCut1[nMomentum][nThetas];
	TH1F *Omega14MassCut2[nMomentum][nThetas];
	TH1F *Omega14MassCut3[nMomentum][nThetas];
	TH1F *Omega14MassCut4[nMomentum][nThetas];
	TH1F *Omega14MassCut5[nMomentum][nThetas];
	TH1F *Omega14MassCut6[nMomentum][nThetas];
	TH1F *Omega14MassCut7[nMomentum][nThetas];

	for (int p = 0; p < nMomentum; p++)
	{
		for (int th = 0; th < nThetas; th++)
		{
			Omega14Mass[p][th] = new TH1F(Form("Omega14Mass_%i_%i", p, th), "", 1000, 0.2, 1.7);
			Omega14Mass[p][th]->SetTitle("Recoil Omega Mass No Cuts;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
			
			Omega14MassCut1[p][th] = new TH1F(Form("Omega14MassCut1_%i_%i", p, th), "", 1000, 0.2, 1.7);
			Omega14MassCut1[p][th]->SetTitle("M(#gamma_{1/2}#gamma_{4/3}) cut;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");
			
			Omega14MassCut2[p][th] = new TH1F(Form("Omega14MassCut2_%i_%i", p, th), "", 1000, 0.2, 1.7);
			Omega14MassCut2[p][th]->SetTitle("Anti #omega(#gamma_{2}#gamma_{3}) cut;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");

			Omega14MassCut3[p][th] = new TH1F(Form("Omega14MassCut3_%i_%i", p, th), "", 1000, 0.2, 1.7);
			Omega14MassCut3[p][th]->SetTitle("Anti M(#gamma_{1}#gamma_{2}) cut;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");

			Omega14MassCut4[p][th] = new TH1F(Form("Omega14MassCut4_%i_%i", p, th), "", 1000, 0.2, 1.7);
			Omega14MassCut4[p][th]->SetTitle("Anti M(#gamma_{1}#gamma_{3}) cut;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");

			Omega14MassCut5[p][th] = new TH1F(Form("Omega14MassCut5_%i_%i", p, th), "", 1000, 0.2, 1.7);
			Omega14MassCut5[p][th]->SetTitle("Anti M(#gamma_{2}#gamma_{4}) cut;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");

			Omega14MassCut6[p][th] = new TH1F(Form("Omega14MassCut6_%i_%i", p, th), "", 1000, 0.2, 1.7);
			Omega14MassCut6[p][th]->SetTitle("Anti M(#gamma_{3}#gamma_{4}) cut;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");

			Omega14MassCut7[p][th] = new TH1F(Form("Omega14MassCut7_%i_%i", p, th), "", 1000, 0.2, 1.7);
			Omega14MassCut7[p][th]->SetTitle("#gamma_{4} in BCAL cut;Recoil #omega(#pi^{+}#pi^{-}[#gamma_{1}(#gamma_{4})]) [GeV];Counts / 1.5MeV");

		}
	}

	for (unsigned int i = 0; i < nEntries; i++)
	{
		OmegaPi0Eff_FlatTree->GetEntry(i);
		if (i % 1000000 == 0) cout << "Event: " << i << " of " << nEntries << " processed..." <<endl;
		
		float theta = P4_g4->Theta()*TMath::RadToDeg();
		float momentum = P4_g4->Energy();
		int thetaBin = (int)(theta/5.);
		int momentumBin = (int)((momentum)/0.50);
		
		if (thetaBin >= nThetas) continue;
		if (momentumBin >= nMomentum) continue;

		// Make the cuts
		if (momentum < 0.10) continue;
		if (momentum > 3.0) continue;
		if (theta < 0.) continue;
		if (theta > 180.) continue;

		// Fill in uncut histogram
		Omega14Mass[momentumBin][thetaBin]->Fill(omega_14_mass, accidweight);
		if (pi0_14_passed && pi0_23_passed) 
		{
			Omega14MassCut1[momentumBin][thetaBin]->Fill(omega_14_mass, accidweight);
		}
		if (pi0_14_passed && pi0_23_passed && antipi0_12_passed) 
		{
			Omega14MassCut3[momentumBin][thetaBin]->Fill(omega_14_mass, accidweight);
		}
		if (pi0_14_passed && pi0_23_passed && antipi0_12_passed && antipi0_13_passed) 
		{
			Omega14MassCut4[momentumBin][thetaBin]->Fill(omega_14_mass, accidweight);
		}
		if (pi0_14_passed && pi0_23_passed && antipi0_12_passed && antipi0_13_passed && antipi0_24_passed) 
		{
			Omega14MassCut5[momentumBin][thetaBin]->Fill(omega_14_mass, accidweight);
		}
		if (pi0_14_passed && pi0_23_passed && antipi0_12_passed && antipi0_13_passed && antipi0_24_passed && antipi0_34_passed) 
		{
			Omega14MassCut6[momentumBin][thetaBin]->Fill(omega_14_mass, accidweight);
		}
		if (pi0_14_passed && pi0_23_passed && antipi0_12_passed && antipi0_13_passed && antipi0_24_passed && antipi0_34_passed && antiomega_23_passed) 
		{
			Omega14MassCut2[momentumBin][thetaBin]->Fill(omega_14_mass, accidweight);
		}
		if (pi0_14_passed && pi0_23_passed && antipi0_12_passed && antipi0_13_passed && antipi0_24_passed && antipi0_34_passed && antiomega_23_passed && (theta > 10.)) 
		{
			Omega14MassCut7[momentumBin][thetaBin]->Fill(omega_14_mass, accidweight);
		}

		
	}
	
	TFile *outFile = new TFile(outFileName, "RECREATE");
	for (int p = 0; p < nMomentum; p++)
	{
		for (int th = 0; th < nThetas; th++)
		{
			Omega14Mass[p][th]->Write();
			Omega14MassCut1[p][th]->Write();
			Omega14MassCut2[p][th]->Write();
			Omega14MassCut3[p][th]->Write();
			Omega14MassCut4[p][th]->Write();
			Omega14MassCut5[p][th]->Write();
			Omega14MassCut6[p][th]->Write();
			Omega14MassCut7[p][th]->Write();
		}
	}
	outFile->Close();
}
