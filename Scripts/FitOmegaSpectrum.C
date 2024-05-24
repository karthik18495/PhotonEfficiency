

//***************************************//
//**************************************//
//Defining all Signal functions here...//
//************************************//
//***********************************//

Double_t gaus(Double_t *x, Double_t *par)
{
    return par[0]*TMath::Exp(-0.5*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]));
}
Double_t doublegaus(Double_t *x, Double_t *par)
{
    return gaus (x, par) + gaus (x, &par[3]);
}
Double_t jz_3gaus(Double_t* x, Double_t* parms) {
	Double_t amplitude1 = parms[0];
	Double_t mean1      = parms[1];
	Double_t sigma1     = parms[2];
	Double_t amplitude2 = parms[3];
	Double_t mean2      = parms[4];
	Double_t sigma2     = parms[5];
	Double_t amplitude3 = parms[6];
	Double_t mean3      = parms[7]; 
	Double_t sigma3     = parms[8]; 

	Double_t gaus_return1 = parms[0]* (1/parms[2]) * (1/sqrt(2*3.14159) ) *exp(-(x[0]-parms[1])*(x[0]-parms[1])/(2*parms[2]*parms[2]) );
	Double_t gaus_return2 = parms[0+3]* (1/parms[2+3]) * (1/sqrt(2*3.14159) ) *exp(-(x[0]-parms[1+3])*(x[0]-parms[1+3])/(2*parms[2+3]*parms[2+3]) );
	Double_t gaus_return3 = parms[0+6]* (1/parms[2+6]) * (1/sqrt(2*3.14159) ) *exp(-(x[0]-parms[1+6])*(x[0]-parms[1+6])/(2*parms[2+6]*parms[2+6]) );
	
	return gaus_return1+gaus_return2+gaus_return3;
}

Double_t jz_2gaus(Double_t* x, Double_t* parms) {
	
	Double_t sig_yield    = parms[0];
	Double_t mean1        = parms[1];
	Double_t sigma1       = parms[2];
	Double_t rel_fraction = parms[3];
	Double_t mean2        = parms[4]; 
	Double_t sigma2       = parms[5]; 

	Double_t gaus_return1 = parms[0]* (1/parms[2]) * (1/sqrt(2*3.14159) ) *exp(-(x[0]-parms[1])*(x[0]-parms[1])/(2*parms[2]*parms[2]) );
	Double_t gaus_return2 = parms[0+3]* (1/parms[2+3]) * (1/sqrt(2*3.14159) ) *exp(-(x[0]-parms[1+3])*(x[0]-parms[1+3])/(2*parms[2+3]*parms[2+3]) );
	
	return sig_yield*(rel_fraction*gaus_return1+(1-rel_fraction)*gaus_return2); 
}

//********************************************//
//*******************************************//
//**Defining all background functions here**//
//*****************************************//
//****************************************//

Double_t poly1(Double_t* x, Double_t* parms) {
    Double_t slope     = parms[0];
    Double_t offset    = parms[1];

    return slope*x[0]+offset;
}

Double_t poly2(Double_t*x, Double_t* parms) {
    Double_t p0 = parms[0];
    Double_t p1     = parms[1];
    Double_t p2    = parms[2];

    return p0*x[0]*x[0]+p1*x[0]+p2;
}

Double_t expo(Double_t* x, Double_t* parms) {
    Double_t amplitude = parms[0];
    Double_t slope     = parms[1];
    
    return amplitude*exp(slope*x[0]);
}

Double_t expo2(Double_t* x, Double_t* parms) {
    Double_t amplitude = parms[0];
    Double_t p0     = parms[1];
    Double_t p1    = parms[2];
    
    return amplitude*exp(p0*x[0]*x[0]+p1*x[0]);
}

Double_t jz_gaus2_pol1(Double_t* x, Double_t* parms) {
    return poly1(x, parms) + jz_2gaus(x, &parms[2]);
}

Double_t jz_gaus2_pol2(Double_t* x, Double_t* parms) {
    return poly2(x, parms) + jz_2gaus(x, &parms[3]);
}

Double_t jz_gaus2_expo(Double_t* x, Double_t* parms) {
    return expo(x, parms) + jz_2gaus(x, &parms[2]);
}

Double_t jz_gaus2_expo2(Double_t* x, Double_t* parms) {
    return expo2(x, parms) + jz_2gaus(x, &parms[3]);
}

Double_t jz_gaus3_pol1(Double_t* x, Double_t* parms) {
    return poly1(x, parms) + jz_3gaus(x, &parms[2]);
}

Double_t jz_gaus3_pol2(Double_t* x, Double_t* parms) {
    return poly2(x, parms) + jz_3gaus(x, &parms[3]);
}

Double_t jz_gaus3_expo(Double_t* x, Double_t* parms) {
    return expo(x, parms) + jz_3gaus(x, &parms[2]);
}

Double_t jz_gaus3_expo2(Double_t* x, Double_t* parms) {
    return expo2(x, parms) + jz_3gaus(x, &parms[3]);
}

// FitFunc4
Double_t FITFUNC4(Double_t *x, Double_t *par)
{
    return gaus (x, par) + gaus (x, &par[3]) + expo2 (x, &par[6]);
}

struct FitUtils {
    TF1* RecoilFitFunc[3][20];
    TF1* InvMassFitFunc[3][20];
    TF1* RecoilSigFunc[3][20];
    TF1* InvMassSigFunc[3][20];
    TF1* RecoilBkgFunc[3][20];
    TF1* InvMassBkgFunc[3][20];

    double NDF[3][20];
    double Chi2[3][20];

    double YieldRecoil[3][20];
    double YieldRecoilErr[3][20];
    double BkgRecoil[3][20];
    double BkgRecoilErr[3][20];

    double YieldInv[3][20];
    double YieldInvErr[3][20];
    double BkgInv[3][20];
    double BkgInvErr[3][20];

    double Effic[3][20];
    double EfficErr[3][20];

    double PurityRecoil[3][20];
    double PurityRecoilErr[3][20];

    double PurityInv[3][20];
    double PurityInvErr[3][20];

};

void FitOmegaSpectrum(TString inFileName)
{
    gStyle->SetOptFit(111);
    TFile *_file = TFile::Open(inFileName);
    const int nMomentum = 3;
    const int nTheta = 20;
    const int NFITS = 20;
    const double lowMomBins[nMomentum] = {0.10, 0.50, 1.00};
    const double highMomBins[nMomentum] = {0.50, 1.00, 3.00};
    const double lowthetaBins[nTheta] =  {10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 25., 28., 31., 35., 40., 45., 50., 60.};
    const double highthetaBins[nTheta] = {11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 25., 28., 31., 35., 40., 45., 50., 60., 80.};
    
    const double lowMassRange = 0.65;
    const double highMassRange = 0.95;
    const double omegaMass = 0.782;
    const double lowMassWidth = 0.008;
    const double highMassWidth = 0.10;

    TH1F *histOmegaRecoilAfterCuts[nMomentum][nTheta];
    TH1F *histOmegaInvMassAfterCuts[nMomentum][nTheta];

    FitUtils JZ_GAUSSIAN_2_POL1;
    FitUtils JZ_GAUSSIAN_2_POL2;
    FitUtils JZ_GAUSSIAN_2_EXPO;
    FitUtils JZ_GAUSSIAN_2_EXPO2;

    FitUtils JZ_GAUSSIAN_3_POL1;
    FitUtils JZ_GAUSSIAN_3_POL2;
    FitUtils JZ_GAUSSIAN_3_EXPO;
    FitUtils JZ_GAUSSIAN_3_EXPO2;

    TCanvas *C1[nMomentum][nTheta];
    TCanvas *C2[nMomentum][nTheta];
    TF1 *temp_fit, *temp_fit_inv;
    for (int i = 0; i < nMomentum; i++)
    {
        for (int j = 0; j < nTheta; j++)
        {
            cout << "i = " << i << " j = " << j << endl;
            C1[i][j] = new TCanvas(Form("C1_%i_%i", i, j), Form("C1_%i_%i", i, j), 1200, 800);
            C1[i][j]->Divide(2, 2);

            C2[i][j] = new TCanvas(Form("C2_%i_%i", i, j), Form("C2_%i_%i", i, j), 1200, 800);
            C2[i][j]->Divide(2, 2);

            histOmegaRecoilAfterCuts[i][j] = (TH1F *)_file->Get(Form("Omega14MassCut7_%i_%i", i, j));
            histOmegaRecoilAfterCuts[i][j]->Rebin(4);
            int binl = histOmegaRecoilAfterCuts[i][j]->GetXaxis()->FindBin(lowMassRange);
            int binh = histOmegaRecoilAfterCuts[i][j]->GetXaxis()->FindBin(highMassRange);
            double Max = histOmegaRecoilAfterCuts[i][j]->GetMaximum();

            histOmegaInvMassAfterCuts[i][j] = (TH1F *)_file->Get(Form("Omega14InvariantMassCut7_%i_%i", i, j));
            histOmegaInvMassAfterCuts[i][j]->Rebin(4);
            int binlInv = histOmegaInvMassAfterCuts[i][j]->GetXaxis()->FindBin(lowMassRange);
            int binhInv = histOmegaInvMassAfterCuts[i][j]->GetXaxis()->FindBin(highMassRange);
            double MaxInv = histOmegaInvMassAfterCuts[i][j]->GetMaximum();

            int lowMass = histOmegaRecoilAfterCuts[i][j]->FindFirstBinAbove(lowMassRange, 1);
            int highMass = histOmegaRecoilAfterCuts[i][j]->FindFirstBinAbove(highMassRange, 1);
            cout << histOmegaRecoilAfterCuts[i][j]->GetNbinsX() << " ; " << binl << " ; " << binh << endl;
            cout << "p = " << lowMomBins[i]/2 + highMomBins[i]/2 << " th = " << lowthetaBins[j]/2. + highthetaBins[j]/2. << " , Efficiency = " << 
            histOmegaInvMassAfterCuts[i][j]->Integral(binl, binh)/histOmegaRecoilAfterCuts[i][j]->Integral(binlInv, binhInv) << endl;
            //continue;
            // Fit function of jz_2gaus + pol1
            temp_fit = new TF1("temp_fit", "gaus(0) + gaus(3)", 0.65, 0.85);
            temp_fit->SetParameters(0.6*Max, omegaMass, lowMassWidth, 0.3*Max, omegaMass, lowMassWidth);
            temp_fit->SetParLimits(0, 0., Max);
            temp_fit->SetParLimits(1, lowMassRange, highMassRange);
            temp_fit->SetParLimits(2, lowMassWidth, highMassWidth);
            temp_fit->SetParLimits(3, 0., Max);
            temp_fit->SetParLimits(4, lowMassRange, lowMassRange);
            temp_fit->SetParLimits(5, lowMassWidth, highMassWidth);
            histOmegaRecoilAfterCuts[i][j]->Fit(temp_fit, "R0");
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL1_RecoilFitFunc_%i_%i", i, j), jz_gaus2_pol1, lowMassRange, highMassRange, 8);
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetParameters(1., 0., temp_fit->GetParameter(0), 
                                                                            temp_fit->GetParameter(1), temp_fit->GetParameter(2), 
                                                                            temp_fit->GetParameter(3), temp_fit->GetParameter(4), 
                                                                            temp_fit->GetParameter(5)
                                                                );
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetParLimits(2, 0.5, Max);
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetParLimits(3, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetParLimits(4, 0.003, temp_fit->GetParameter(2));
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetParLimits(5, 0., 0.5);
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetParLimits(6, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetParLimits(7, 0.003, temp_fit->GetParameter(5));
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetLineColor(kRed);
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->SetLineStyle(2);
            
            histOmegaRecoilAfterCuts[i][j]->Fit(JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j], "R");

            histOmegaRecoilAfterCuts[i][j]->Draw("HIST");
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->Draw("SAME");
            
            JZ_GAUSSIAN_2_POL1.RecoilSigFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL1_RecoilSigFunc_%i_%i", i, j), jz_2gaus, lowMassRange, highMassRange, 6);
            JZ_GAUSSIAN_2_POL1.RecoilSigFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameter(2), 
                                                                    JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameter(3), 
                                                                    JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameter(4), 
                                                                    JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameter(5),
                                                                    JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameter(6),
                                                                    JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameter(7)
                                                                    );
            JZ_GAUSSIAN_2_POL1.RecoilSigFunc[i][j]->SetLineColor(kBlue); // Signal is blue
            JZ_GAUSSIAN_2_POL1.RecoilSigFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_POL1.RecoilSigFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_POL1.RecoilBkgFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL1_RecoilBkgFunc_%i_%i", i, j), poly1, lowMassRange, highMassRange, 2);
            JZ_GAUSSIAN_2_POL1.RecoilBkgFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameter(0), 
                                                                    JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameter(1)
                                                                    );
            JZ_GAUSSIAN_2_POL1.RecoilBkgFunc[i][j]->SetLineColor(kBlack); // Background is Black
            JZ_GAUSSIAN_2_POL1.RecoilBkgFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_POL1.RecoilBkgFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_POL1.YieldRecoil[i][j] = JZ_GAUSSIAN_2_POL1.RecoilSigFunc[i][j]->Integral(lowMassRange, highMassRange);

            JZ_GAUSSIAN_2_POL1.YieldRecoilErr[i][j] = JZ_GAUSSIAN_2_POL1.RecoilSigFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_POL1.BkgRecoil[i][j] = JZ_GAUSSIAN_2_POL1.RecoilBkgFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_POL1.BkgRecoilErr[i][j] = JZ_GAUSSIAN_2_POL1.RecoilBkgFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_POL1.PurityRecoil[i][j] = (JZ_GAUSSIAN_2_POL1.YieldRecoil[i][j] - JZ_GAUSSIAN_2_POL1.BkgRecoil[i][j])/JZ_GAUSSIAN_2_POL1.YieldRecoil[i][j];
            JZ_GAUSSIAN_2_POL1.PurityRecoilErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_POL1.YieldRecoilErr[i][j]/JZ_GAUSSIAN_2_POL1.YieldRecoil[i][j], 2) + pow(JZ_GAUSSIAN_2_POL1.BkgRecoilErr[i][j]/JZ_GAUSSIAN_2_POL1.BkgRecoil[i][j], 2))*JZ_GAUSSIAN_2_POL1.PurityRecoil[i][j];


            temp_fit_inv = new TF1("temp_fit_inv", "gaus(0) + gaus(3)", lowMassRange, highMassRange);
            temp_fit_inv->SetParameters(0.6*MaxInv, omegaMass, lowMassWidth, 0.3*MaxInv, omegaMass, lowMassWidth);
            temp_fit_inv->SetParLimits(0, 0., MaxInv);
            temp_fit_inv->SetParLimits(1, lowMassRange, highMassRange);
            temp_fit_inv->SetParLimits(2, lowMassWidth, highMassWidth);
            temp_fit_inv->SetParLimits(3, 0., MaxInv);
            temp_fit_inv->SetParLimits(4, lowMassRange, lowMassRange);
            temp_fit_inv->SetParLimits(5, lowMassWidth, highMassWidth);

            histOmegaInvMassAfterCuts[i][j]->Fit(temp_fit_inv, "R0");
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL1_InvMassFitFunc_%i_%i", i, j), jz_gaus2_pol1, lowMassRange, highMassRange, 8);
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->SetParameters(0., 0., temp_fit_inv->GetParameter(0), 
                                                                            temp_fit_inv->GetParameter(1), temp_fit_inv->GetParameter(2), 
                                                                            temp_fit_inv->GetParameter(3), temp_fit_inv->GetParameter(4), 
                                                                            temp_fit_inv->GetParameter(5)
                                                                );
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->SetParLimits(2, 0., MaxInv);
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->SetParLimits(3, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->SetParLimits(4, 0.003, temp_fit_inv->GetParameter(2));
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->SetParLimits(5, 0., 0.5);
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->SetParLimits(6, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->SetParLimits(7, 0.003, temp_fit_inv->GetParameter(5));
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->SetLineColor(kRed);
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->SetLineWidth(2);


            histOmegaInvMassAfterCuts[i][j]->Fit(JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j], "R0");

            JZ_GAUSSIAN_2_POL1.InvMassSigFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL1_InvMassSigFunc_%i_%i", i, j), jz_2gaus, lowMassRange, highMassRange, 6);
            JZ_GAUSSIAN_2_POL1.InvMassSigFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameter(2), 
                                                                    JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameter(3), 
                                                                    JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameter(4), 
                                                                    JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameter(5),
                                                                    JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameter(6),
                                                                    JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameter(7)
                                                                    );
            JZ_GAUSSIAN_2_POL1.InvMassSigFunc[i][j]->SetLineColor(kBlue); // Signal is blue
            JZ_GAUSSIAN_2_POL1.InvMassSigFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_POL1.InvMassSigFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_POL1.InvMassBkgFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL1_InvMassBkgFunc_%i_%i", i, j), poly1, lowMassRange, highMassRange, 2);
            JZ_GAUSSIAN_2_POL1.InvMassBkgFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameter(0), 
                                                                    JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameter(1)
                                                                    );
            JZ_GAUSSIAN_2_POL1.InvMassBkgFunc[i][j]->SetLineColor(kBlack); // Background is Black
            JZ_GAUSSIAN_2_POL1.InvMassBkgFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_POL1.InvMassBkgFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_POL1.YieldInv[i][j] = JZ_GAUSSIAN_2_POL1.InvMassSigFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_POL1.YieldInvErr[i][j] = JZ_GAUSSIAN_2_POL1.InvMassSigFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_POL1.BkgInv[i][j] = JZ_GAUSSIAN_2_POL1.InvMassBkgFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_POL1.BkgInvErr[i][j] = JZ_GAUSSIAN_2_POL1.InvMassBkgFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_POL1.PurityInv[i][j] = (JZ_GAUSSIAN_2_POL1.YieldInv[i][j] - JZ_GAUSSIAN_2_POL1.BkgInv[i][j])/JZ_GAUSSIAN_2_POL1.YieldInv[i][j];
            JZ_GAUSSIAN_2_POL1.PurityInvErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_POL1.YieldInvErr[i][j]/JZ_GAUSSIAN_2_POL1.YieldInv[i][j], 2) + pow(JZ_GAUSSIAN_2_POL1.BkgInvErr[i][j]/JZ_GAUSSIAN_2_POL1.BkgInv[i][j], 2))*JZ_GAUSSIAN_2_POL1.PurityInv[i][j];

            JZ_GAUSSIAN_2_POL1.Effic[i][j] = JZ_GAUSSIAN_2_POL1.YieldInv[i][j]/JZ_GAUSSIAN_2_POL1.YieldRecoil[i][j];
            JZ_GAUSSIAN_2_POL1.EfficErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_POL1.YieldInvErr[i][j]/JZ_GAUSSIAN_2_POL1.YieldInv[i][j], 2) + pow(JZ_GAUSSIAN_2_POL1.YieldRecoilErr[i][j]/JZ_GAUSSIAN_2_POL1.YieldRecoil[i][j], 2))*JZ_GAUSSIAN_2_POL1.Effic[i][j];

            // Fit function of jz_2gaus + pol2

            temp_fit = new TF1("temp_fit", "gaus(0) + gaus(3)", lowMassRange, highMassRange);
            temp_fit->SetParameters(0.6*Max, omegaMass, lowMassWidth, 0.3*Max, omegaMass, lowMassWidth);
            temp_fit->SetParLimits(0, 0., Max);
            temp_fit->SetParLimits(1, lowMassRange, highMassRange);
            temp_fit->SetParLimits(2, lowMassWidth, highMassWidth);
            temp_fit->SetParLimits(3, 0., Max);
            temp_fit->SetParLimits(4, lowMassRange, lowMassRange);
            temp_fit->SetParLimits(5, lowMassWidth, highMassWidth);
            histOmegaRecoilAfterCuts[i][j]->Fit(temp_fit, "R0");
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL2_RecoilFitFunc_%i_%i", i, j), jz_gaus2_pol2, lowMassRange, highMassRange, 9);
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->SetParameters(0., 0., 0., temp_fit->GetParameter(0), 
                                                                            temp_fit->GetParameter(1), temp_fit->GetParameter(2), 
                                                                            temp_fit->GetParameter(3), temp_fit->GetParameter(4), 
                                                                            temp_fit->GetParameter(5)
                                                                );
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->SetParLimits(3, 0., Max);
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->SetParLimits(4, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->SetParLimits(5, 0.003, temp_fit->GetParameter(2));
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->SetParLimits(6, 0., 0.5);
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->SetParLimits(7, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->SetParLimits(8, 0.003, temp_fit->GetParameter(5));
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->SetLineColor(kRed);
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->SetLineWidth(2);


            histOmegaRecoilAfterCuts[i][j]->Fit(JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j], "R0");
            JZ_GAUSSIAN_2_POL2.RecoilSigFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL2_RecoilSigFunc_%i_%i", i, j), jz_2gaus, lowMassRange, highMassRange, 6);
            JZ_GAUSSIAN_2_POL2.RecoilSigFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameter(3), 
                                                                    JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameter(4), 
                                                                    JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameter(5), 
                                                                    JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameter(6),
                                                                    JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameter(7),
                                                                    JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameter(8)
                                                                    );
            JZ_GAUSSIAN_2_POL2.RecoilSigFunc[i][j]->SetLineColor(kBlue); // Signal is blue
            JZ_GAUSSIAN_2_POL2.RecoilSigFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_POL2.RecoilSigFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_POL2.RecoilBkgFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL2_RecoilBkgFunc_%i_%i", i, j), poly2, lowMassRange, highMassRange, 3);
            JZ_GAUSSIAN_2_POL2.RecoilBkgFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameter(0), 
                                                                    JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameter(1),
                                                                    JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameter(2)
                                                                    );
            JZ_GAUSSIAN_2_POL2.RecoilBkgFunc[i][j]->SetLineColor(kBlack); // Background is Black
            JZ_GAUSSIAN_2_POL2.RecoilBkgFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_POL2.RecoilBkgFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_POL2.YieldRecoil[i][j] = JZ_GAUSSIAN_2_POL2.RecoilSigFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_POL2.YieldRecoilErr[i][j] = JZ_GAUSSIAN_2_POL2.RecoilSigFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_POL2.BkgRecoil[i][j] = JZ_GAUSSIAN_2_POL2.RecoilBkgFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_POL2.BkgRecoilErr[i][j] = JZ_GAUSSIAN_2_POL2.RecoilBkgFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_POL2.PurityRecoil[i][j] = (JZ_GAUSSIAN_2_POL2.YieldRecoil[i][j] - JZ_GAUSSIAN_2_POL2.BkgRecoil[i][j])/JZ_GAUSSIAN_2_POL2.YieldRecoil[i][j];
            JZ_GAUSSIAN_2_POL2.PurityRecoilErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_POL2.YieldRecoilErr[i][j]/JZ_GAUSSIAN_2_POL2.YieldRecoil[i][j], 2) + pow(JZ_GAUSSIAN_2_POL2.BkgRecoilErr[i][j]/JZ_GAUSSIAN_2_POL2.BkgRecoil[i][j], 2))*JZ_GAUSSIAN_2_POL2.PurityRecoil[i][j];


            temp_fit_inv = new TF1("temp_fit_inv", "gaus(0) + gaus(3)", lowMassRange, highMassRange);
            temp_fit_inv->SetParameters(0.6*MaxInv, omegaMass, lowMassWidth, 0.3*MaxInv, omegaMass, lowMassWidth);
            temp_fit_inv->SetParLimits(0, 0., MaxInv);
            temp_fit_inv->SetParLimits(1, lowMassRange, highMassRange);
            temp_fit_inv->SetParLimits(2, lowMassWidth, highMassWidth);
            temp_fit_inv->SetParLimits(3, 0., MaxInv);
            temp_fit_inv->SetParLimits(4, lowMassRange, lowMassRange);
            temp_fit_inv->SetParLimits(5, lowMassWidth, highMassWidth);

            histOmegaInvMassAfterCuts[i][j]->Fit(temp_fit_inv, "R0");
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL2_InvMassFitFunc_%i_%i", i, j), jz_gaus2_pol2, lowMassRange, highMassRange, 9);
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->SetParameters(0., 0., 0., temp_fit_inv->GetParameter(0), 
                                                                            temp_fit_inv->GetParameter(1), temp_fit_inv->GetParameter(2), 
                                                                            temp_fit_inv->GetParameter(3), temp_fit_inv->GetParameter(4), 
                                                                            temp_fit_inv->GetParameter(5)
                                                                );
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->SetParLimits(3, 0., MaxInv);
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->SetParLimits(4, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->SetParLimits(5, 0.003, temp_fit_inv->GetParameter(2));
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->SetParLimits(6, 0., 0.5);
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->SetParLimits(7, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->SetParLimits(8, 0.003, temp_fit_inv->GetParameter(5));
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->SetLineColor(kRed);
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->SetLineWidth(2);

            histOmegaInvMassAfterCuts[i][j]->Fit(JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j], "R0");

            JZ_GAUSSIAN_2_POL2.InvMassSigFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL2_InvMassSigFunc_%i_%i", i, j), jz_2gaus, lowMassRange, highMassRange, 6);
            JZ_GAUSSIAN_2_POL2.InvMassSigFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameter(3), 
                                                                    JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameter(4), 
                                                                    JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameter(5), 
                                                                    JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameter(6),
                                                                    JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameter(7),
                                                                    JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameter(8)
                                                                    );
            JZ_GAUSSIAN_2_POL2.InvMassSigFunc[i][j]->SetLineColor(kBlue); // Signal is blue
            JZ_GAUSSIAN_2_POL2.InvMassSigFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_POL2.InvMassSigFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_POL2.InvMassBkgFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_POL2_InvMassBkgFunc_%i_%i", i, j), poly2, lowMassRange, highMassRange, 3);
            JZ_GAUSSIAN_2_POL2.InvMassBkgFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameter(0), 
                                                                    JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameter(1),
                                                                    JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameter(2)
                                                                    );
            JZ_GAUSSIAN_2_POL2.InvMassBkgFunc[i][j]->SetLineColor(kBlack); // Background is Black
            JZ_GAUSSIAN_2_POL2.InvMassBkgFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_POL2.InvMassBkgFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_POL2.YieldInv[i][j] = JZ_GAUSSIAN_2_POL2.InvMassSigFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_POL2.YieldInvErr[i][j] = JZ_GAUSSIAN_2_POL2.InvMassSigFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_POL2.BkgInv[i][j] = JZ_GAUSSIAN_2_POL2.InvMassBkgFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_POL2.BkgInvErr[i][j] = JZ_GAUSSIAN_2_POL2.InvMassBkgFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_POL2.PurityInv[i][j] = (JZ_GAUSSIAN_2_POL2.YieldInv[i][j] - JZ_GAUSSIAN_2_POL2.BkgInv[i][j])/JZ_GAUSSIAN_2_POL2.YieldInv[i][j];
            JZ_GAUSSIAN_2_POL2.PurityInvErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_POL2.YieldInvErr[i][j]/JZ_GAUSSIAN_2_POL2.YieldInv[i][j], 2) + pow(JZ_GAUSSIAN_2_POL2.BkgInvErr[i][j]/JZ_GAUSSIAN_2_POL2.BkgInv[i][j], 2))*JZ_GAUSSIAN_2_POL2.PurityInv[i][j];

            JZ_GAUSSIAN_2_POL2.Effic[i][j] = JZ_GAUSSIAN_2_POL2.YieldInv[i][j]/JZ_GAUSSIAN_2_POL2.YieldRecoil[i][j];
            JZ_GAUSSIAN_2_POL2.EfficErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_POL2.YieldInvErr[i][j]/JZ_GAUSSIAN_2_POL2.YieldInv[i][j], 2) + pow(JZ_GAUSSIAN_2_POL2.YieldRecoilErr[i][j]/JZ_GAUSSIAN_2_POL2.YieldRecoil[i][j], 2))*JZ_GAUSSIAN_2_POL2.Effic[i][j];


            // Fit function of jz_2gaus + expo

            temp_fit = new TF1("temp_fit", "gaus(0) + gaus(3)", lowMassRange, highMassRange);
            temp_fit->SetParameters(0.6*Max, omegaMass, lowMassWidth, 0.3*Max, omegaMass, lowMassWidth);
            temp_fit->SetParLimits(0, 0., Max);
            temp_fit->SetParLimits(1, lowMassRange, highMassRange);
            temp_fit->SetParLimits(2, lowMassWidth, highMassWidth);
            temp_fit->SetParLimits(3, 0., Max);
            temp_fit->SetParLimits(4, lowMassRange, lowMassRange);
            temp_fit->SetParLimits(5, lowMassWidth, highMassWidth);
            histOmegaRecoilAfterCuts[i][j]->Fit(temp_fit, "R0");
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO_RecoilFitFunc_%i_%i", i, j), jz_gaus2_expo, lowMassRange, highMassRange, 8);
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->SetParameters(0., 0., temp_fit->GetParameter(0), 
                                                                            temp_fit->GetParameter(1), temp_fit->GetParameter(2), 
                                                                            temp_fit->GetParameter(3), temp_fit->GetParameter(4), 
                                                                            temp_fit->GetParameter(5)
                                                                );
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->SetParLimits(2, 0., Max);
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->SetParLimits(3, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->SetParLimits(4, 0.003, temp_fit->GetParameter(2));
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->SetParLimits(5, 0., 0.5);
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->SetParLimits(6, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->SetParLimits(7, 0.003, temp_fit->GetParameter(5));
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->SetLineColor(kRed);
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->SetLineWidth(2);


            histOmegaRecoilAfterCuts[i][j]->Fit(JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j], "R0");
            JZ_GAUSSIAN_2_EXPO.RecoilSigFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO_RecoilSigFunc_%i_%i", i, j), jz_2gaus, lowMassRange, highMassRange, 6);
            JZ_GAUSSIAN_2_EXPO.RecoilSigFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameter(2), 
                                                                    JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameter(3), 
                                                                    JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameter(4), 
                                                                    JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameter(5),
                                                                    JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameter(6),
                                                                    JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameter(7)
                                                                    );
            JZ_GAUSSIAN_2_EXPO.RecoilSigFunc[i][j]->SetLineColor(kBlue); // Signal is blue
            JZ_GAUSSIAN_2_EXPO.RecoilSigFunc[i][j]->SetLineWidth(2);

            JZ_GAUSSIAN_2_EXPO.RecoilBkgFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO_RecoilBkgFunc_%i_%i", i, j), expo, lowMassRange, highMassRange, 2);
            JZ_GAUSSIAN_2_EXPO.RecoilBkgFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameter(0), 
                                                                    JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameter(1)
                                                                    );
            JZ_GAUSSIAN_2_EXPO.RecoilBkgFunc[i][j]->SetLineColor(kBlack); // Background is Black
            JZ_GAUSSIAN_2_EXPO.RecoilBkgFunc[i][j]->SetLineWidth(2);

            JZ_GAUSSIAN_2_EXPO.YieldRecoil[i][j] = JZ_GAUSSIAN_2_EXPO.RecoilSigFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_EXPO.YieldRecoilErr[i][j] = JZ_GAUSSIAN_2_EXPO.RecoilSigFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_EXPO.BkgRecoil[i][j] = JZ_GAUSSIAN_2_EXPO.RecoilBkgFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_EXPO.BkgRecoilErr[i][j] = JZ_GAUSSIAN_2_EXPO.RecoilBkgFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_EXPO.PurityRecoil[i][j] = (JZ_GAUSSIAN_2_EXPO.YieldRecoil[i][j] - JZ_GAUSSIAN_2_EXPO.BkgRecoil[i][j])/JZ_GAUSSIAN_2_EXPO.YieldRecoil[i][j];
            JZ_GAUSSIAN_2_EXPO.PurityRecoilErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_EXPO.YieldRecoilErr[i][j]/JZ_GAUSSIAN_2_EXPO.YieldRecoil[i][j], 2) + pow(JZ_GAUSSIAN_2_EXPO.BkgRecoilErr[i][j]/JZ_GAUSSIAN_2_EXPO.BkgRecoil[i][j], 2))*JZ_GAUSSIAN_2_EXPO.PurityRecoil[i][j];


            temp_fit_inv = new TF1("temp_fit_inv", "gaus(0) + gaus(3)", lowMassRange, highMassRange);
            temp_fit_inv->SetParameters(0.6*MaxInv, omegaMass, lowMassWidth, 0.3*MaxInv, omegaMass, lowMassWidth);
            temp_fit_inv->SetParLimits(0, 0., MaxInv);
            temp_fit_inv->SetParLimits(1, lowMassRange, highMassRange);
            temp_fit_inv->SetParLimits(2, lowMassWidth, highMassWidth);
            temp_fit_inv->SetParLimits(3, 0., MaxInv);
            temp_fit_inv->SetParLimits(4, lowMassRange, lowMassRange);
            temp_fit_inv->SetParLimits(5, lowMassWidth, highMassWidth);

            histOmegaInvMassAfterCuts[i][j]->Fit(temp_fit_inv, "R0");
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO_InvMassFitFunc_%i_%i", i, j), jz_gaus2_expo, lowMassRange, highMassRange, 8);
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->SetParameters(0., 0., temp_fit_inv->GetParameter(0), 
                                                                            temp_fit_inv->GetParameter(1), temp_fit_inv->GetParameter(2), 
                                                                            temp_fit_inv->GetParameter(3), temp_fit_inv->GetParameter(4), 
                                                                            temp_fit_inv->GetParameter(5)
                                                                );
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->SetParLimits(3, 0., MaxInv);
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->SetParLimits(4, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->SetParLimits(5, 0.003, temp_fit_inv->GetParameter(2));
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->SetParLimits(6, 0., 0.5);
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->SetParLimits(7, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->SetParLimits(8, 0.003, temp_fit_inv->GetParameter(5));
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->SetLineColor(kRed);
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->SetLineWidth(2);

            histOmegaInvMassAfterCuts[i][j]->Fit(JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j], "R0");

            JZ_GAUSSIAN_2_EXPO.InvMassSigFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO_InvMassSigFunc_%i_%i", i, j), jz_2gaus, lowMassRange, highMassRange, 6);
            JZ_GAUSSIAN_2_EXPO.InvMassSigFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameter(2), 
                                                                    JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameter(3), 
                                                                    JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameter(4), 
                                                                    JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameter(5),
                                                                    JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameter(6),
                                                                    JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameter(7)
                                                                    );
            JZ_GAUSSIAN_2_EXPO.InvMassSigFunc[i][j]->SetLineColor(kBlue); // Signal is blue
            JZ_GAUSSIAN_2_EXPO.InvMassSigFunc[i][j]->SetLineWidth(2);

            JZ_GAUSSIAN_2_EXPO.InvMassBkgFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO_InvMassBkgFunc_%i_%i", i, j), expo, lowMassRange, highMassRange, 2);
            JZ_GAUSSIAN_2_EXPO.InvMassBkgFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameter(0), 
                                                                    JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameter(1)
                                                                    );
            JZ_GAUSSIAN_2_EXPO.InvMassBkgFunc[i][j]->SetLineColor(kBlack); // Background is Black
            JZ_GAUSSIAN_2_EXPO.InvMassBkgFunc[i][j]->SetLineWidth(2);

            JZ_GAUSSIAN_2_EXPO.YieldInv[i][j] = JZ_GAUSSIAN_2_EXPO.InvMassSigFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_EXPO.YieldInvErr[i][j] = JZ_GAUSSIAN_2_EXPO.InvMassSigFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_EXPO.BkgInv[i][j] = JZ_GAUSSIAN_2_EXPO.InvMassBkgFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_EXPO.BkgInvErr[i][j] = JZ_GAUSSIAN_2_EXPO.InvMassBkgFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_EXPO.PurityInv[i][j] = (JZ_GAUSSIAN_2_EXPO.YieldInv[i][j] - JZ_GAUSSIAN_2_EXPO.BkgInv[i][j])/JZ_GAUSSIAN_2_EXPO.YieldInv[i][j];
            JZ_GAUSSIAN_2_EXPO.PurityInvErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_EXPO.YieldInvErr[i][j]/JZ_GAUSSIAN_2_EXPO.YieldInv[i][j], 2) + pow(JZ_GAUSSIAN_2_EXPO.BkgInvErr[i][j]/JZ_GAUSSIAN_2_EXPO.BkgInv[i][j], 2))*JZ_GAUSSIAN_2_EXPO.PurityInv[i][j];

            JZ_GAUSSIAN_2_EXPO.Effic[i][j] = JZ_GAUSSIAN_2_EXPO.YieldInv[i][j]/JZ_GAUSSIAN_2_EXPO.YieldRecoil[i][j];
            JZ_GAUSSIAN_2_EXPO.EfficErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_EXPO.YieldInvErr[i][j]/JZ_GAUSSIAN_2_EXPO.YieldInv[i][j], 2) + pow(JZ_GAUSSIAN_2_EXPO.YieldRecoilErr[i][j]/JZ_GAUSSIAN_2_EXPO.YieldRecoil[i][j], 2))*JZ_GAUSSIAN_2_EXPO.Effic[i][j];


            // Fit function of jz_2gaus + expo2

            temp_fit = new TF1("temp_fit", "gaus(0) + gaus(3)", lowMassRange, highMassRange);
            temp_fit->SetParameters(0.6*Max, omegaMass, lowMassWidth, 0.3*Max, omegaMass, lowMassWidth);
            temp_fit->SetParLimits(0, 0., Max);
            temp_fit->SetParLimits(1, lowMassRange, highMassRange);
            temp_fit->SetParLimits(2, lowMassWidth, highMassWidth);
            temp_fit->SetParLimits(3, 0., Max);
            temp_fit->SetParLimits(4, lowMassRange, lowMassRange);
            temp_fit->SetParLimits(5, lowMassWidth, highMassWidth);

            histOmegaRecoilAfterCuts[i][j]->Fit(temp_fit, "R0");

            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO2_RecoilFitFunc_%i_%i", i, j), jz_gaus2_expo2, lowMassRange, highMassRange, 9);
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->SetParameters(0., 0., 0., temp_fit->GetParameter(0), 
                                                                            temp_fit->GetParameter(1), temp_fit->GetParameter(2), 
                                                                            temp_fit->GetParameter(3), temp_fit->GetParameter(4), 
                                                                            temp_fit->GetParameter(5)
                                                                );
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->SetParLimits(3, 0., Max);
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->SetParLimits(4, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->SetParLimits(5, 0.003, temp_fit->GetParameter(2));
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->SetParLimits(6, 0., 0.5);
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->SetParLimits(7, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->SetParLimits(8, 0.003, temp_fit->GetParameter(5));
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->SetLineColor(kRed);
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->SetLineWidth(2);

            histOmegaRecoilAfterCuts[i][j]->Fit(JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j], "R0");
            JZ_GAUSSIAN_2_EXPO2.RecoilSigFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO2_RecoilSigFunc_%i_%i", i, j), jz_2gaus, lowMassRange, highMassRange, 6);
            JZ_GAUSSIAN_2_EXPO2.RecoilSigFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameter(3), 
                                                                    JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameter(4), 
                                                                    JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameter(5), 
                                                                    JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameter(6),
                                                                    JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameter(7),
                                                                    JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameter(8)
                                                                    );
            JZ_GAUSSIAN_2_EXPO2.RecoilSigFunc[i][j]->SetLineColor(kBlue); // Signal is blue
            JZ_GAUSSIAN_2_EXPO2.RecoilSigFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_EXPO2.RecoilSigFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_EXPO2.RecoilBkgFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO2_RecoilBkgFunc_%i_%i", i, j), expo2, lowMassRange, highMassRange, 3);
            JZ_GAUSSIAN_2_EXPO2.RecoilBkgFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameter(0), 
                                                                    JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameter(1),
                                                                    JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameter(2)
                                                                    );
            JZ_GAUSSIAN_2_EXPO2.RecoilBkgFunc[i][j]->SetLineColor(kBlack); // Background is Black
            JZ_GAUSSIAN_2_EXPO2.RecoilBkgFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_EXPO2.RecoilBkgFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_EXPO2.YieldRecoil[i][j] = JZ_GAUSSIAN_2_EXPO2.RecoilSigFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_EXPO2.YieldRecoilErr[i][j] = JZ_GAUSSIAN_2_EXPO2.RecoilSigFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_EXPO2.BkgRecoil[i][j] = JZ_GAUSSIAN_2_EXPO2.RecoilBkgFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_EXPO2.BkgRecoilErr[i][j] = JZ_GAUSSIAN_2_EXPO2.RecoilBkgFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_EXPO2.PurityRecoil[i][j] = (JZ_GAUSSIAN_2_EXPO2.YieldRecoil[i][j] - JZ_GAUSSIAN_2_EXPO2.BkgRecoil[i][j])/JZ_GAUSSIAN_2_EXPO2.YieldRecoil[i][j];
            JZ_GAUSSIAN_2_EXPO2.PurityRecoilErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_EXPO2.YieldRecoilErr[i][j]/JZ_GAUSSIAN_2_EXPO2.YieldRecoil[i][j], 2) + pow(JZ_GAUSSIAN_2_EXPO2.BkgRecoilErr[i][j]/JZ_GAUSSIAN_2_EXPO2.BkgRecoil[i][j], 2))*JZ_GAUSSIAN_2_EXPO2.PurityRecoil[i][j];


            temp_fit_inv = new TF1("temp_fit_inv", "gaus(0) + gaus(3)", lowMassRange, highMassRange); // Change to single gaussian.
            temp_fit_inv->SetParameters(0.6*MaxInv, omegaMass, lowMassWidth, 0.3*MaxInv, omegaMass, lowMassWidth);
            temp_fit_inv->SetParLimits(0, 0., MaxInv);
            temp_fit_inv->SetParLimits(1, lowMassRange, highMassRange);
            temp_fit_inv->SetParLimits(2, lowMassWidth, highMassWidth);
            temp_fit_inv->SetParLimits(3, 0., MaxInv);
            temp_fit_inv->SetParLimits(4, lowMassRange, lowMassRange);
            temp_fit_inv->SetParLimits(5, lowMassWidth, highMassWidth);

            histOmegaInvMassAfterCuts[i][j]->Fit(temp_fit_inv, "R0");
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO2_InvMassFitFunc_%i_%i", i, j), jz_gaus2_expo2, lowMassRange, highMassRange, 9);
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->SetParameters(0., 0., 0., temp_fit_inv->GetParameter(0), 
                                                                            temp_fit_inv->GetParameter(1), temp_fit_inv->GetParameter(2), 
                                                                            temp_fit_inv->GetParameter(3), temp_fit_inv->GetParameter(4), 
                                                                            temp_fit_inv->GetParameter(5)
                                                                );
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->SetParLimits(3, 0., MaxInv);
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->SetParLimits(4, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->SetParLimits(5, 0.003, temp_fit_inv->GetParameter(2));
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->SetParLimits(6, 0., 0.5);
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->SetParLimits(7, lowMassRange, highMassRange); //
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->SetParLimits(8, 0.003, temp_fit_inv->GetParameter(5));
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->SetLineColor(kRed);
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->SetLineWidth(2);

            histOmegaInvMassAfterCuts[i][j]->Fit(JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j], "R0");

            JZ_GAUSSIAN_2_EXPO2.InvMassSigFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO2_InvMassSigFunc_%i_%i", i, j), jz_2gaus, lowMassRange, highMassRange, 6);
            JZ_GAUSSIAN_2_EXPO2.InvMassSigFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameter(3), 
                                                                    JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameter(4), 
                                                                    JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameter(5), 
                                                                    JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameter(6),
                                                                    JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameter(7),
                                                                    JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameter(8)
                                                                    );
            JZ_GAUSSIAN_2_EXPO2.InvMassSigFunc[i][j]->SetLineColor(kBlue); // Signal is blue
            JZ_GAUSSIAN_2_EXPO2.InvMassSigFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_EXPO2.InvMassSigFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_EXPO2.InvMassBkgFunc[i][j] = new TF1(Form("JZ_GAUSSIAN_2_EXPO2_InvMassBkgFunc_%i_%i", i, j), expo2, lowMassRange, highMassRange, 3);
            JZ_GAUSSIAN_2_EXPO2.InvMassBkgFunc[i][j]->SetParameters(JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameter(0), 
                                                                    JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameter(1),
                                                                    JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameter(2)
                                                                    );
            JZ_GAUSSIAN_2_EXPO2.InvMassBkgFunc[i][j]->SetLineColor(kBlack); // Background is Black
            JZ_GAUSSIAN_2_EXPO2.InvMassBkgFunc[i][j]->SetLineWidth(2);
            JZ_GAUSSIAN_2_EXPO2.InvMassBkgFunc[i][j]->SetLineStyle(2);

            JZ_GAUSSIAN_2_EXPO2.YieldInv[i][j] = JZ_GAUSSIAN_2_EXPO2.InvMassSigFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_EXPO2.YieldInvErr[i][j] = JZ_GAUSSIAN_2_EXPO2.InvMassSigFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_EXPO2.BkgInv[i][j] = JZ_GAUSSIAN_2_EXPO2.InvMassBkgFunc[i][j]->Integral(lowMassRange, highMassRange);
            JZ_GAUSSIAN_2_EXPO2.BkgInvErr[i][j] = JZ_GAUSSIAN_2_EXPO2.InvMassBkgFunc[i][j]->IntegralError(lowMassRange, highMassRange, JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->GetParameters());

            JZ_GAUSSIAN_2_EXPO2.PurityInv[i][j] = (JZ_GAUSSIAN_2_EXPO2.YieldInv[i][j] - JZ_GAUSSIAN_2_EXPO2.BkgInv[i][j])/JZ_GAUSSIAN_2_EXPO2.YieldInv[i][j];
            JZ_GAUSSIAN_2_EXPO2.PurityInvErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_EXPO2.YieldInvErr[i][j]/JZ_GAUSSIAN_2_EXPO2.YieldInv[i][j], 2) + pow(JZ_GAUSSIAN_2_EXPO2.BkgInvErr[i][j]/JZ_GAUSSIAN_2_EXPO2.BkgInv[i][j], 2))*JZ_GAUSSIAN_2_EXPO2.PurityInv[i][j];

            JZ_GAUSSIAN_2_EXPO2.Effic[i][j] = JZ_GAUSSIAN_2_EXPO2.YieldInv[i][j]/JZ_GAUSSIAN_2_EXPO2.YieldRecoil[i][j];
            JZ_GAUSSIAN_2_EXPO2.EfficErr[i][j] = sqrt(pow(JZ_GAUSSIAN_2_EXPO2.YieldInvErr[i][j]/JZ_GAUSSIAN_2_EXPO2.YieldInv[i][j], 2) + pow(JZ_GAUSSIAN_2_EXPO2.YieldRecoilErr[i][j]/JZ_GAUSSIAN_2_EXPO2.YieldRecoil[i][j], 2))*JZ_GAUSSIAN_2_EXPO2.Effic[i][j];


            // Fit function of jz_2gaus + pol1 draw
            C1[i][j]->cd(1);
            TLatex *lat1 = new TLatex(0.3, 0.6*Max, "pol1 + jz_gaus2");
            histOmegaRecoilAfterCuts[i][j]->Draw("HIST");
            JZ_GAUSSIAN_2_POL1.RecoilFitFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_POL1.RecoilSigFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_POL1.RecoilBkgFunc[i][j]->Draw("SAME");
            C1[i][j]->cd(1)->Update();
            //TPaveStats *stats1 = (TPaveStats*)histOmegaRecoilAfterCuts[i][j]->FindObject("stats")->Clone();
            //stats1->SetStats(0);
            //stats1->Draw("SAME");
            lat1->Draw("SAME");

            // Fit function of jz_2gaus + pol2 draw
            C1[i][j]->cd(2);
            TLatex *lat2 = new TLatex(0.3, 0.6*Max, "pol2 + jz_gaus2");
            histOmegaRecoilAfterCuts[i][j]->Draw("HIST");
            JZ_GAUSSIAN_2_POL2.RecoilFitFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_POL2.RecoilSigFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_POL2.RecoilBkgFunc[i][j]->Draw("SAME");
            C1[i][j]->cd(2)->Update();
            //TPaveStats *stats2 = (TPaveStats*)histOmegaInvMassAfterCuts[i][j]->FindObject("stats")->Clone();
            //stats2->SetStats(0);
            //stats2->Draw("SAME");
            lat2->Draw("SAME");

            // Fit function of jz_2gaus + expo draw
            C1[i][j]->cd(3);
            TLatex *lat3 = new TLatex(0.3, 0.6*Max, "expo + jz_gaus2");
            histOmegaRecoilAfterCuts[i][j]->Draw("HIST");
            JZ_GAUSSIAN_2_EXPO.RecoilFitFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_EXPO.RecoilSigFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_EXPO.RecoilBkgFunc[i][j]->Draw("SAME");
            C1[i][j]->cd(3)->Update();
            //TPaveStats *stats3 = (TPaveStats*)histOmegaRecoilAfterCuts[i][j]->FindObject("stats")->Clone();
            //stats3->SetStats(0);
            //stats3->Draw("SAME");
            lat3->Draw("SAME");

            // Fit function of jz_2gaus + expo2 draw
            C1[i][j]->cd(4);
            TLatex *lat4 = new TLatex(0.3, 0.6*Max, "expo2 + jz_gaus2");
            histOmegaRecoilAfterCuts[i][j]->Draw("HIST");
            JZ_GAUSSIAN_2_EXPO2.RecoilFitFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_EXPO2.RecoilSigFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_EXPO2.RecoilBkgFunc[i][j]->Draw("SAME");
            C1[i][j]->cd(4)->Update();
            //TPaveStats *stats4 = (TPaveStats*)histOmegaRecoilAfterCuts[i][j]->FindObject("stats")->Clone();
            //stats4->SetStats(0);
            //stats4->Draw("SAME");
            lat4->Draw("SAME");

            // Fit function of jz_2gaus + pol1 inv draw
            C2[i][j]->cd(1);
            TLatex *latInv = new TLatex(0.3, 0.6*MaxInv, "pol1 + jz_gaus2");
            histOmegaInvMassAfterCuts[i][j]->Draw("HIST");
            JZ_GAUSSIAN_2_POL1.InvMassFitFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_POL1.InvMassSigFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_POL1.InvMassBkgFunc[i][j]->Draw("SAME");
            C2[i][j]->cd(1)->Update();
            //TPaveStats *statsInv = (TPaveStats*)histOmegaInvMassAfterCuts[i][j]->FindObject("stats")->Clone();
            //statsInv->SetStats(0);
            //statsInv->Draw("SAME");
            latInv->Draw("SAME");

            // Fit function of jz_2gaus + pol2 inv draw
            C2[i][j]->cd(2);
            TLatex *latInv2 = new TLatex(0.3, 0.6*MaxInv, "pol2 + jz_gaus2");
            histOmegaInvMassAfterCuts[i][j]->Draw("HIST");
            JZ_GAUSSIAN_2_POL2.InvMassFitFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_POL2.InvMassSigFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_POL2.InvMassBkgFunc[i][j]->Draw("SAME");
            C2[i][j]->cd(2)->Update();
            TPaveStats *statsInv2 = (TPaveStats*)histOmegaInvMassAfterCuts[i][j]->FindObject("stats")->Clone();
            //statsInv2->SetStats(0);
            //statsInv2->Draw("SAME");
            latInv2->Draw("SAME");

            // Fit function of jz_2gaus + expo inv draw
            C2[i][j]->cd(3);
            TLatex *latInv3 = new TLatex(0.3, 0.6*MaxInv, "expo + jz_gaus2");
            histOmegaInvMassAfterCuts[i][j]->Draw("HIST");
            JZ_GAUSSIAN_2_EXPO.InvMassFitFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_EXPO.InvMassSigFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_EXPO.InvMassBkgFunc[i][j]->Draw("SAME");
            C2[i][j]->cd(3)->Update();
            TPaveStats *statsInv3 = (TPaveStats*)histOmegaInvMassAfterCuts[i][j]->FindObject("stats")->Clone();
            //statsInv3->SetStats(0);
            //statsInv3->Draw("SAME");
            latInv3->Draw("SAME");

            // Fit function of jz_2gaus + expo2 inv draw
            C2[i][j]->cd(4);
            TLatex *latInv4 = new TLatex(0.3, 0.6*MaxInv, "expo2 + jz_gaus2");
            histOmegaInvMassAfterCuts[i][j]->Draw("HIST");
            JZ_GAUSSIAN_2_EXPO2.InvMassFitFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_EXPO2.InvMassSigFunc[i][j]->Draw("SAME");
            JZ_GAUSSIAN_2_EXPO2.InvMassBkgFunc[i][j]->Draw("SAME");
            C2[i][j]->cd(4)->Update();
            TPaveStats *statsInv4 = (TPaveStats*)histOmegaInvMassAfterCuts[i][j]->FindObject("stats")->Clone();
            //statsInv4->SetStats(0);
            //statsInv4->Draw("SAME");
            latInv4->Draw("SAME");


            // Print the canvas
            
            if(i == 0 && j == 0) 
            {
                C1[i][j]->Print("RecoilOmegaMassFit.pdf(");
                
            }
            else if (i == nMomentum - 1 && j == nTheta - 1) 
            {
                C1[i][j]->Print("RecoilOmegaMassFit.pdf)");
                
            }
            else 
            {
                C1[i][j]->Print("RecoilOmegaMassFit.pdf");
            }

            if (i == 0 && j == 0)
            {
                C2[i][j]->Print("InvMassOmegaMassFit.pdf(");
            }
            else if (i == nMomentum - 1 && j == nTheta - 1)
            {
                C2[i][j]->Print("InvMassOmegaMassFit.pdf)");
                cout << "Done!" << endl;
            }
            else
            {
                C2[i][j]->Print("InvMassOmegaMassFit.pdf");
            }
        }
    }
    
    // Efficiency and Purity as a function of theta in 3 momentum
    TGraphErrors *grEff_JZ_GAUSSIAN_2_POL1[nMomentum], *grEff_JZ_GAUSSIAN_2_POL2[nMomentum], *grEff_JZ_GAUSSIAN_2_EXPO[nMomentum], *grEff_JZ_GAUSSIAN_2_EXPO2[nMomentum];
    TGraphErrors *grRecoilPur_JZ_GAUSSIAN_2_POL1[nMomentum], *grRecoilPur_JZ_GAUSSIAN_2_POL2[nMomentum], *grRecoilPur_JZ_GAUSSIAN_2_EXPO[nMomentum], *grRecoilPur_JZ_GAUSSIAN_2_EXPO2[nMomentum];
    TGraphErrors *grInvPur_JZ_GAUSSIAN_2_POL1[nMomentum], *grInvPur_JZ_GAUSSIAN_2_POL2[nMomentum], *grInvPur_JZ_GAUSSIAN_2_EXPO[nMomentum], *grInvPur_JZ_GAUSSIAN_2_EXPO2[nMomentum];

    const int nColors[nMomentum] = {kRed, kBlue, kMagenta};
    const int nMStyles[4] = {20, 21, 22, 23};
    TLegend *leg = new TLegend(0.2, 0.1, 0.6, 0.4);
    TLegend *leg2 = new TLegend(0.2, 0.1, 0.6, 0.4);
    TCanvas *C3[nMomentum], *C4[nMomentum], *C5[nMomentum];
    for (int p = 0; p < nMomentum; p++)
        {
            grEff_JZ_GAUSSIAN_2_POL1[p] = new TGraphErrors(nTheta);
            grEff_JZ_GAUSSIAN_2_POL1[p]->SetName(Form("Efficiency_JZ_gaus_2_pol1_p_%i", p));
            grEff_JZ_GAUSSIAN_2_POL1[p]->SetMarkerStyle(20);
            grEff_JZ_GAUSSIAN_2_POL1[p]->SetMarkerSize(1.5);
            grEff_JZ_GAUSSIAN_2_POL1[p]->SetMarkerColor(nColors[p]);
            grEff_JZ_GAUSSIAN_2_POL1[p]->SetLineColor(kBlack);
            grEff_JZ_GAUSSIAN_2_POL1[p]->SetLineWidth(2);
            grEff_JZ_GAUSSIAN_2_POL1[p]->SetTitle("Efficiency vs. #theta; #theta [deg]; Efficiency");
            leg->AddEntry(grEff_JZ_GAUSSIAN_2_POL1[p], Form("pol1 + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grEff_JZ_GAUSSIAN_2_POL1[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL1.Effic[p][t]);
                grEff_JZ_GAUSSIAN_2_POL1[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL1.EfficErr[p][t]);
            }
            grEff_JZ_GAUSSIAN_2_POL1[p]->SetMinimum(0.);
            grEff_JZ_GAUSSIAN_2_POL1[p]->SetMaximum(1.);

            grEff_JZ_GAUSSIAN_2_POL2[p] = new TGraphErrors(nTheta);
            grEff_JZ_GAUSSIAN_2_POL2[p]->SetName(Form("Efficiency_JZ_gaus_2_pol2_p_%i", p));
            grEff_JZ_GAUSSIAN_2_POL2[p]->SetMarkerStyle(22);
            grEff_JZ_GAUSSIAN_2_POL2[p]->SetMarkerSize(1.5);
            grEff_JZ_GAUSSIAN_2_POL2[p]->SetMarkerColor(nColors[p]);
            grEff_JZ_GAUSSIAN_2_POL2[p]->SetLineColor(kBlack);
            grEff_JZ_GAUSSIAN_2_POL2[p]->SetLineWidth(2);
            grEff_JZ_GAUSSIAN_2_POL2[p]->SetTitle("Efficiency vs. #theta; #theta [deg]; Efficiency");
            leg->AddEntry(grEff_JZ_GAUSSIAN_2_POL2[p], Form("pol2 + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grEff_JZ_GAUSSIAN_2_POL2[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL2.Effic[p][t]);
                grEff_JZ_GAUSSIAN_2_POL2[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL2.EfficErr[p][t]);
            }
            grEff_JZ_GAUSSIAN_2_POL2[p]->SetMinimum(0.);
            grEff_JZ_GAUSSIAN_2_POL2[p]->SetMaximum(1.);

            grEff_JZ_GAUSSIAN_2_EXPO[p] = new TGraphErrors(nTheta);
            grEff_JZ_GAUSSIAN_2_EXPO[p]->SetName(Form("Efficiency_JZ_gaus_2_expo_p_%i", p));
            grEff_JZ_GAUSSIAN_2_EXPO[p]->SetMarkerStyle(24);
            grEff_JZ_GAUSSIAN_2_EXPO[p]->SetMarkerSize(1.5);
            grEff_JZ_GAUSSIAN_2_EXPO[p]->SetMarkerColor(nColors[p]);
            grEff_JZ_GAUSSIAN_2_EXPO[p]->SetLineColor(kBlack);
            grEff_JZ_GAUSSIAN_2_EXPO[p]->SetLineWidth(2);
            grEff_JZ_GAUSSIAN_2_EXPO[p]->SetTitle("Efficiency vs. #theta; #theta [deg]; Efficiency");
            leg->AddEntry(grEff_JZ_GAUSSIAN_2_EXPO[p], Form("expo + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grEff_JZ_GAUSSIAN_2_EXPO[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO.Effic[p][t]);
                grEff_JZ_GAUSSIAN_2_EXPO[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO.EfficErr[p][t]);
            }
            grEff_JZ_GAUSSIAN_2_EXPO[p]->SetMinimum(0.);
            grEff_JZ_GAUSSIAN_2_EXPO[p]->SetMaximum(1.);

            grEff_JZ_GAUSSIAN_2_EXPO2[p] = new TGraphErrors(nTheta);
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetName(Form("Efficiency_JZ_gaus_2_expo2_p_%i", p));
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetMarkerStyle(26);
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetMarkerSize(1.5);
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetMarkerColor(nColors[p]);
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetLineColor(kBlack);
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetLineWidth(2);
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetTitle(Form("Efficiency vs. #theta, %f < p < %f [GeV]; #theta [deg]; Efficiency", lowMomBins[p], highMomBins[p]));
            leg->AddEntry(grEff_JZ_GAUSSIAN_2_EXPO2[p], Form("expo2 + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                int binl = histOmegaRecoilAfterCuts[p][t]->GetXaxis()->FindBin(lowMassRange);
                int binh = histOmegaRecoilAfterCuts[p][t]->GetXaxis()->FindBin(highMassRange);
                int binlInv = histOmegaInvMassAfterCuts[p][t]->GetXaxis()->FindBin(lowMassRange);
                int binhInv = histOmegaInvMassAfterCuts[p][t]->GetXaxis()->FindBin(highMassRange);
                double eff = histOmegaInvMassAfterCuts[p][t]->Integral(binl, binh)/histOmegaRecoilAfterCuts[p][t]->Integral(binl, binh);
                grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., eff);
                grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO2.EfficErr[p][t]);
            }
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetMinimum(0.);
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->SetMaximum(1.);

            C3[p] = new TCanvas (Form("Efficiency_%i", p), Form("Efficiency_%i", p), 1200, 800);
            C3[p]->cd();
            //grEff_JZ_GAUSSIAN_2_POL1[p]->Draw("AP");
            //grEff_JZ_GAUSSIAN_2_POL2[p]->Draw("P SAME");
            //grEff_JZ_GAUSSIAN_2_EXPO[p]->Draw("P SAME");
            grEff_JZ_GAUSSIAN_2_EXPO2[p]->Draw("AP");
            leg->Draw("SAME");
            if(p == 0) {
                C3[p]->Print(Form("Efficiency_JZ_GAUS2.pdf("));
            }
            else if (p == nMomentum - 1) {
                C3[p]->Print(Form("Efficiency_JZ_GAUS2.pdf)"));
            }
            else {
                C3[p]->Print(Form("Efficiency_JZ_GAUS2.pdf"));
            }

            grRecoilPur_JZ_GAUSSIAN_2_POL1[p] = new TGraphErrors(nTheta);
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetName(Form("RecoilPurity_JZ_gaus_2_pol1_p_%i", p));
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetMarkerStyle(20);
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetMarkerSize(1.5);
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetMarkerColor(nColors[p]);
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetLineColor(kBlack);
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetLineWidth(2);
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetTitle("Recoil Purity vs. #theta; #theta [deg]; Purity");
            leg->AddEntry(grRecoilPur_JZ_GAUSSIAN_2_POL1[p], Form("pol1 + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL1.PurityRecoil[p][t]);
                grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL1.PurityRecoilErr[p][t]);
            }
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetMinimum(0.);
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->SetMaximum(1.);

            grRecoilPur_JZ_GAUSSIAN_2_POL2[p] = new TGraphErrors(nTheta);
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetName(Form("RecoilPurity_JZ_gaus_2_pol2_p_%i", p));
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetMarkerStyle(20);
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetMarkerSize(1.5);
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetMarkerColor(nColors[p]);
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetLineColor(kBlack);
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetLineWidth(2);
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetTitle("Recoil Purity vs. #theta; #theta [deg]; Purity");
            leg->AddEntry(grRecoilPur_JZ_GAUSSIAN_2_POL2[p], Form("pol2 + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL2.PurityRecoil[p][t]);
                grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL2.PurityRecoilErr[p][t]);
            }
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetMinimum(0.);
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->SetMaximum(1.);

            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p] = new TGraphErrors(nTheta);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetName(Form("RecoilPurity_JZ_gaus_2_expo_p_%i", p));
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetMarkerStyle(20);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetMarkerSize(1.5);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetMarkerColor(nColors[p]);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetLineColor(kBlack);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetLineWidth(2);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetTitle("Recoil Purity vs. #theta; #theta [deg]; Purity");
            leg->AddEntry(grRecoilPur_JZ_GAUSSIAN_2_EXPO[p], Form("expo + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO.PurityRecoil[p][t]);
                grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO.PurityRecoilErr[p][t]);
            }
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetMinimum(0.);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->SetMaximum(1.);

            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p] = new TGraphErrors(nTheta);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetName(Form("RecoilPurity_JZ_gaus_2_expo2_p_%i", p));
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMarkerStyle(20);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMarkerSize(1.5);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMarkerColor(nColors[p]);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetLineColor(kBlack);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetLineWidth(2);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetTitle("Recoil Purity vs. #theta; #theta [deg]; Purity");
            leg->AddEntry(grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p], Form("expo2 + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO2.PurityRecoil[p][t]);
                grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO2.PurityRecoilErr[p][t]);
            }
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMinimum(0.);
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMaximum(1.);

            C4[p] = new TCanvas (Form("RecoilPurity_%i", p), Form("RecoilPurity_%i", p), 1200, 800);
            C4[p]->cd();
            grRecoilPur_JZ_GAUSSIAN_2_POL1[p]->Draw("AP");
            grRecoilPur_JZ_GAUSSIAN_2_POL2[p]->Draw("P SAME");
            grRecoilPur_JZ_GAUSSIAN_2_EXPO[p]->Draw("P SAME");
            grRecoilPur_JZ_GAUSSIAN_2_EXPO2[p]->Draw("P SAME");
            leg->Draw("SAME");
            if(p == 0) {
                C4[p]->Print(Form("RecoilPurity_JZ_GAUS2.pdf("));
            }
            else if (p == nMomentum - 1) {
                C4[p]->Print(Form("RecoilPurity_JZ_GAUS2.pdf)"));
            }
            else {
                C4[p]->Print(Form("RecoilPurity_JZ_GAUS2.pdf"));
            }

            grInvPur_JZ_GAUSSIAN_2_POL1[p] = new TGraphErrors(nTheta);
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetName(Form("InvPurity_JZ_gaus_2_pol1_p_%i", p));
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetMarkerStyle(20);
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetMarkerSize(1.5);
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetMarkerColor(nColors[p]);
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetLineColor(kBlack);
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetLineWidth(2);
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetTitle("InvMass Purity vs. #theta; #theta [deg]; Purity");
            leg->AddEntry(grInvPur_JZ_GAUSSIAN_2_POL1[p], Form("pol1 + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL1.PurityInv[p][t]);
                grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL1.PurityInvErr[p][t]);
            }
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetMinimum(0.);
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->SetMaximum(1.);

            grInvPur_JZ_GAUSSIAN_2_POL2[p] = new TGraphErrors(nTheta);
            grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetName(Form("InvPurity_JZ_gaus_2_pol2_p_%i", p));
            grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetMarkerStyle(20);
            grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetMarkerSize(1.5);
            grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetMarkerColor(nColors[p]);

            grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetLineColor(kBlack);
            grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetLineWidth(2);
            grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetTitle("InvMass Purity vs. #theta; #theta [deg]; Purity");
            leg->AddEntry(grInvPur_JZ_GAUSSIAN_2_POL2[p], Form("pol2 + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL2.PurityInv[p][t]);
                grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_POL2.PurityInvErr[p][t]);
            }
            grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetMinimum(0.);
            grInvPur_JZ_GAUSSIAN_2_POL2[p]->SetMaximum(1.);

            grInvPur_JZ_GAUSSIAN_2_EXPO[p] = new TGraphErrors(nTheta);
            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetName(Form("InvPurity_JZ_gaus_2_expo_p_%i", p));
            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetMarkerStyle(20);
            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetMarkerSize(1.5);
            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetMarkerColor(nColors[p]);

            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetLineColor(kBlack);
            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetLineWidth(2);
            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetTitle("InvMass Purity vs. #theta; #theta [deg]; Purity");
            leg->AddEntry(grInvPur_JZ_GAUSSIAN_2_EXPO[p], Form("expo + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO.PurityInv[p][t]);
                grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO.PurityInvErr[p][t]);
            }
            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetMinimum(0.);
            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->SetMaximum(1.);

            grInvPur_JZ_GAUSSIAN_2_EXPO2[p] = new TGraphErrors(nTheta);
            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetName(Form("InvPurity_JZ_gaus_2_expo2_p_%i", p));
            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMarkerStyle(20);
            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMarkerSize(1.5);
            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMarkerColor(nColors[p]);

            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetLineColor(kBlack);
            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetLineWidth(2);
            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetTitle("InvMass Purity vs. #theta; #theta [deg]; Purity");
            leg->AddEntry(grInvPur_JZ_GAUSSIAN_2_EXPO2[p], Form("expo2 + jz_gaus2"), "lep");
            for (int t = 0; t < nTheta; t++)
            {
                grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO2.PurityInv[p][t]);
                grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., JZ_GAUSSIAN_2_EXPO2.PurityInvErr[p][t]);
            }
            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMinimum(0.);
            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->SetMaximum(1.);

            C5[p] = new TCanvas (Form("InvPurity_%i", p), Form("InvPurity_%i", p), 1200, 800);
            C5[p]->cd();
            grInvPur_JZ_GAUSSIAN_2_POL1[p]->Draw("AP");
            grInvPur_JZ_GAUSSIAN_2_POL2[p]->Draw("P SAME");
            grInvPur_JZ_GAUSSIAN_2_EXPO[p]->Draw("P SAME");
            grInvPur_JZ_GAUSSIAN_2_EXPO2[p]->Draw("P SAME");
            leg->Draw("SAME");
            if(p == 0) {
                C5[p]->Print(Form("InvPurity_JZ_GAUS2.pdf("));
            }
            else if (p == nMomentum - 1) {
                C5[p]->Print(Form("InvPurity_JZ_GAUS2.pdf)"));
            }
            else {
                C5[p]->Print(Form("InvPurity_JZ_GAUS2.pdf"));
            }

        }
    

}
