void BeautifyPlots(TH1F *hist, int FillColor)
{
 hist->SetFillColorAlpha(FillColor, 0.8); 
}

void MakeAntiPlots(TString Filename, TString Treename)
{

	gROOT->ProcessLine(".x $ROOTSYS/macros/gluex_style.C");
	gStyle->SetBarWidth(0.45);
	TFile* _file0 = TFile::Open(Filename); 
	TTree* OmegaPi0Eff_FlatTree = (TTree*)_file0->Get(Treename);
	
	double pi0_12_mass, pi0_13_mass, pi0_14_mass, pi0_23_mass, pi0_24_mass, pi0_34_mass;
	double omega_12_mass, omega_13_mass, omega_14_mass, omega_23_mass, omega_24_mass, omega_34_mass;
	TLorentzVector *P4_g1 = NULL, *P4_g2 = NULL, *P4_g3 = NULL, *P4_g4 = NULL, *P4_pip = NULL, *P4_pim = NULL;
	double kinfit_CL, ExtraShowerE;
	
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_12_mass", &pi0_12_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_13_mass", &pi0_13_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_14_mass", &pi0_14_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_23_mass", &pi0_23_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_24_mass", &pi0_24_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("pi0_34_mass", &pi0_34_mass);
	
	OmegaPi0Eff_FlatTree->SetBranchAddress("omega_12_mass", &omega_12_mass);	
	OmegaPi0Eff_FlatTree->SetBranchAddress("omega_13_mass", &omega_13_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("omega_14_mass", &omega_14_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("omega_23_mass", &omega_23_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("omega_24_mass", &omega_24_mass);
	OmegaPi0Eff_FlatTree->SetBranchAddress("omega_34_mass", &omega_34_mass);

	OmegaPi0Eff_FlatTree->SetBranchAddress("P4_g1", &P4_g1);
	OmegaPi0Eff_FlatTree->SetBranchAddress("P4_g2", &P4_g2);
	OmegaPi0Eff_FlatTree->SetBranchAddress("P4_g3", &P4_g3);
	OmegaPi0Eff_FlatTree->SetBranchAddress("P4_g4", &P4_g4);
	OmegaPi0Eff_FlatTree->SetBranchAddress("P4_pip", &P4_pip);
	OmegaPi0Eff_FlatTree->SetBranchAddress("P4_pim", &P4_pim);
	
	OmegaPi0Eff_FlatTree->SetBranchAddress("kinfit_CL", &kinfit_CL);
	OmegaPi0Eff_FlatTree->SetBranchAddress("ExtraShowerE", &ExtraShowerE);

	const int nCuts = 5;
	const double pi0M = 0.135;
	const double pi0S = 0.0045936;
	const double omegaM = 0.782;
	const double omegaS = 0.00849;
	const double omegaT = omegaS*3;
	TH1F *omegaHistsPiCuts[nCuts];
	double Cuts[nCuts] = {pi0S*3, pi0S*5, pi0S*10, pi0S*15, pi0S*20};
	TString Labels[nCuts] = {"3#sigma #pi^{0} cut", "5#sigma #pi^{0} cut", "10#sigma #pi^{0} cut", "15#sigma #pi^{0} cut", "20#sigma #pi^{0} cut"};
	const int Colors[nCuts] = {kRed, kBlue, kGreen, kCyan, kYellow};
	
	TLegend *leg = new TLegend(0.75, 0.75, 0.9, 0.9);
	
	const int Entries = OmegaPi0Eff_FlatTree->GetEntries();
	
	for (unsigned int i = 0; i < nCuts; i++)
	{
		omegaHistsPiCuts[i] = new TH1F(Form("PiCutHist_%i", i),";M(#pi^{+}#pi^{-}#gamma_{1}#gamma_{4}) [GeV];Counts / 1.5 MeV", 1000, 0.2, 1.7);
		omegaHistsPiCuts[i]->GetXaxis()->SetTitleColor(kBlack);
		omegaHistsPiCuts[i]->SetFillColorAlpha(Colors[i], 0.8);
		omegaHistsPiCuts[i]->SetLineColor(kBlack);
		omegaHistsPiCuts[i]->SetMinimum(1);
		leg->AddEntry(omegaHistsPiCuts[i], Labels[i], "f");
	}
	
	for (unsigned int i = 0; i < Entries; i++)
	{
		if(i % 100000 == 0) { cout << "Done analyzing " << i << " combos" << endl;}
		
		OmegaPi0Eff_FlatTree->GetEntry(i);
		if(fabs(omega_12_mass - omegaM) > omegaT && fabs(omega_13_mass - omegaM) > omegaT && 
			fabs(omega_23_mass - omegaM) > omegaT && fabs(omega_24_mass - omegaM) > omegaT &&
			fabs(omega_34_mass - omegaM) > omegaT
		)
		{
			for (unsigned j = 0; j < nCuts; j++)
			{
				if(fabs(pi0_12_mass - pi0M) > Cuts[j] && fabs(pi0_13_mass - pi0M > Cuts[j]) &&
					fabs(pi0_24_mass - pi0M) > Cuts[j] && fabs(pi0_34_mass - pi0M) > Cuts[j] &&
					fabs(pi0_14_mass - pi0M) < Cuts[j] && fabs(pi0_23_mass - pi0M) < Cuts[j]
				)
				{
					omegaHistsPiCuts[j]->Fill(omega_14_mass);
				}
			}
		}
	}
	
	TCanvas *C1 = new TCanvas("RecoilMass", "RecoilMass", 1200, 800);
	C1->cd();
	C1->SetLogy();
	omegaHistsPiCuts[2]->Draw("HIST");
	cout << omegaHistsPiCuts[0]->GetEntries() <<endl;
	for(unsigned int i = 0; i < nCuts; i++)
	{
		if (i!=2) omegaHistsPiCuts[i]->Draw("SAME HIST");
		cout << omegaHistsPiCuts[i]->GetEntries() <<endl;
	}
	leg->Draw("SAME");
	//C1->Print("Testing.pdf");
}
