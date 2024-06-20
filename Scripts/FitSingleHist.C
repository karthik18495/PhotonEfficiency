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
	Double_t amplitude1 = parms[0]; // Signal Yield
	Double_t mean1      = parms[1]; // Mean of the first Gaussian (omega mass)
	Double_t sigma1     = parms[2]; // Width of the first Gaussian (omega width)
	Double_t amplitude2 = parms[3]; // Second gaussian yield
	Double_t mean2      = parms[4]; // Mean of the second Gaussian (larger omega mass)
	Double_t sigma2     = parms[5]; // Width of the second Gaussian (larger omega width)
	Double_t amplitude3 = parms[6]; // Third gaussian yield
	Double_t mean3      = parms[7]; // Mean of the third Gaussian (larger omega mass)
	Double_t sigma3     = parms[8]; // Width of the third Gaussian (larger omega width)

	Double_t gaus_return1 = parms[0]* (1/parms[2]) * (1/sqrt(2*3.14159) ) *exp(-(x[0]-parms[1])*(x[0]-parms[1])/(2*parms[2]*parms[2]) );
	Double_t gaus_return2 = parms[0+3]* (1/parms[2+3]) * (1/sqrt(2*3.14159) ) *exp(-(x[0]-parms[1+3])*(x[0]-parms[1+3])/(2*parms[2+3]*parms[2+3]) );
	Double_t gaus_return3 = parms[0+6]* (1/parms[2+6]) * (1/sqrt(2*3.14159) ) *exp(-(x[0]-parms[1+6])*(x[0]-parms[1+6])/(2*parms[2+6]*parms[2+6]) );
	
	return gaus_return1+gaus_return2+gaus_return3;
}

Double_t jz_2gaus(Double_t* x, Double_t* parms) {
	
	Double_t sig_yield    = parms[0]; // Signal Yield
	Double_t mean1        = parms[1]; // Mean of first Gaussian, Close to omega peak
	Double_t sigma1       = parms[2]; // Sigma of first Gaussian, Close to omega width
	Double_t rel_fraction = parms[3]; // Relative fraction of second Gaussian, Close to 0
	Double_t mean2        = parms[4]; // Mean of second Gaussian, usually, the same as the first Gaussian
	Double_t sigma2       = parms[5]; // Sigma of second Gaussian, usually, very large compared to the first Gaussian

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

struct GetNParams
{
    // signal functions
    int jz_3gaus = 9;
    int jz_2gaus = 6;
    int gaus = 3;
    int doublegaus = 6;
    // background functions
    int expo2 = 3;
    int poly2 = 3;
    int poly1 = 2;
    int expo = 2;
};


void FitSingleHist(TString inFileName, TString histName, double min_x = -1, double max_x = -1)
// inFileName: input file name
// histName: histogram name
// min_x: minimum x value of your fit
// max_x: maximum x value of your fit
// usage example: root 'FitSingleHist.C("path/to/Data/Fall2018_PhotonEff.root","Omega14MassCut7_0_9", 0.4, 1.0)'
{
    gStyle->SetOptFit(111);
    TFile *_file = TFile::Open(inFileName);
    if (! _file->Get(histName))
    {
        cout << "Error: " << histName << " not found in file " << inFileName << endl;
        cout << "Available keys:" << endl;
        for (int i = 0; i < _file->GetListOfKeys()->GetEntries(); i++)
        cout << _file->GetListOfKeys()->At(i)->GetName() << endl;
        return;
    }

    TH1F *h = (TH1F*)_file->Get(histName);
    cout << "Working with " << h->GetName() << "..." << endl;

    // Define some useful functions here
    GetNParams nParams;
    int total_params = nParams.jz_2gaus + nParams.poly2;
    min_x = max(h->GetXaxis()->GetXmin(), min_x);
    max_x = max(h->GetXaxis()->GetXmax(), max_x);
    cout << "min_x = " << min_x << ", max_x = " << max_x << endl;
    // We are defining the fit functions in here
    TF1 *fitfunc = new TF1("fitfunc", jz_gaus2_pol2, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), 7);
    // Define functional start points. This will be different for different fit functions.
    // Refer to the respective definitions for the parameters of the fit function.
    fitfunc->SetParameters(0, 0, 0, 0, 0, 0, 0, 0, 0);
    // Set limits on fit parameters. Again this depends on the fit function one uses.
    fitfunc->SetParLimits(0, 0, h->GetMaximum()); 
    fitfunc->SetParLimits(1, -h->GetMaximum(), h->GetMaximum());
    fitfunc->SetParLimits(2, 0, h->GetMaximum());
    fitfunc->SetParLimits(3, -h->GetMaximum(), h->GetMaximum());
    fitfunc->SetParLimits(4, 0, h->GetMaximum());
    fitfunc->SetParLimits(5, -h->GetMaximum(), h->GetMaximum());

    h->Fit(fitfunc, "R");
    h->Fit(fitfunc, "R");

    cout << "Fit results:" << endl;
    cout << "chi2 = " << fitfunc->GetChisquare() << endl;
    cout << "ndf = " << fitfunc->GetNDF() << endl;
    for (int i = 0; i < fitfunc->GetNpar(); i++)
    cout << "par[" << i << "] " << fitfunc->GetParName(i) << " = " << fitfunc->GetParameter(i) << ", err = " << fitfunc->GetParError(i) << endl;
    
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    h->Draw();
    fitfunc->Draw("SAME");
    c->SaveAs("FitSingleHist.pdf");

}