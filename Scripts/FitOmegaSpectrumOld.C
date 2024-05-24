Double_t gaus(Double_t *x, Double_t *par)
{
    return par[0]*TMath::Exp(-0.5*(x[0]-par[1])*(x[0]-par[1])/(par[2]*par[2]));
}
Double_t expo2(Double_t *x, Double_t *par)
{
    return TMath::Exp(par[0]*x[0]*x[0] + par[1]*x[0] + par[2]);
}

Double_t doublegaus(Double_t *x, Double_t *par)
{
    return gaus (x, par) + gaus (x, &par[3]);
}

Double_t FITFUNC4(Double_t *x, Double_t *par)
{
    return gaus (x, par) + gaus (x, &par[3]) + expo2 (x, &par[6]);
}

void FitOmegaSpectrumOld(TString inFileName)
{
    gStyle->SetOptFit(111);
    TFile *_file = TFile::Open(inFileName);
    const int nMomentum = 3;
    const int nTheta = 20;
    const double lowMomBins[nMomentum] = {0.10, 0.50, 1.00};
    const double highMomBins[nMomentum] = {0.50, 1.00, 3.00};
    const double lowthetaBins[nTheta] =  {10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 25., 28., 31., 35., 40., 45., 50., 60.};
    const double highthetaBins[nTheta] = {11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 25., 28., 31., 35., 40., 45., 50., 60., 80.};
    
    double YieldRecoil[nMomentum][nTheta];
    double YieldRecoilErr[nMomentum][nTheta];
    double BkgRecoil[nMomentum][nTheta];
    double BkgRecoilErr[nMomentum][nTheta];

    double YieldInv[nMomentum][nTheta];
    double YieldInvErr[nMomentum][nTheta];
    double BkgInv[nMomentum][nTheta];
    double BkgInvErr[nMomentum][nTheta];

    double Effic[nMomentum][nTheta];
    double EfficErr[nMomentum][nTheta];
 
    TH1F *histOmegaRecoilAfterCuts[nMomentum][nTheta];
    TH1F *histOmegaInvMassAfterCuts[nMomentum][nTheta];
    // Fit Function gaus + gaus + expo
    TF1 *fitFunc1[nMomentum][nTheta];
    TF1 *gaus1_1[nMomentum][nTheta];
    TF1 *gaus2_1[nMomentum][nTheta];
    TF1 *sig1[nMomentum][nTheta];
    TF1 *expo_1[nMomentum][nTheta];

    TF1 *fitFuncInv1[nMomentum][nTheta];
    TF1 *gausInv1_1[nMomentum][nTheta];
    TF1 *gausInv2_1[nMomentum][nTheta];
    TF1 *sigInv1[nMomentum][nTheta];
    TF1 *expoInv_1[nMomentum][nTheta];

    // Fit Function gaus + gaus + pol2
    TF1 *fitFunc2[nMomentum][nTheta];
    TF1 *gaus1_2[nMomentum][nTheta];
    TF1 *gaus2_2[nMomentum][nTheta];
    TF1 *sig2[nMomentum][nTheta];
    TF1 *pol2_2[nMomentum][nTheta];

    TF1 *fitFuncInv2[nMomentum][nTheta];
    TF1 *gausInv1_2[nMomentum][nTheta];
    TF1 *gausInv2_2[nMomentum][nTheta];
    TF1 *sigInv2[nMomentum][nTheta];
    TF1 *polInv2_2[nMomentum][nTheta];

    // Fit Function gaus + gaus + pol1
    TF1 *fitFunc3[nMomentum][nTheta];
    TF1 *gaus1_3[nMomentum][nTheta];
    TF1 *gaus2_3[nMomentum][nTheta];
    TF1 *sig3[nMomentum][nTheta];
    TF1 *pol1_3[nMomentum][nTheta];

    TF1 *fitFuncInv3[nMomentum][nTheta];
    TF1 *gausInv1_3[nMomentum][nTheta];
    TF1 *gausInv2_3[nMomentum][nTheta];
    TF1 *sigInv3[nMomentum][nTheta];
    TF1 *polInv1_3[nMomentum][nTheta];

    // Fit Function gaus + gaus + exp(x**2)
    TF1 *fitFunc4[nMomentum][nTheta];
    TF1 *gaus1_4[nMomentum][nTheta];
    TF1 *gaus2_4[nMomentum][nTheta];
    TF1 *sig4[nMomentum][nTheta];
    TF1 *expo_4[nMomentum][nTheta];

    TF1 *fitFuncInv4[nMomentum][nTheta];
    TF1 *gausInv1_4[nMomentum][nTheta];
    TF1 *gausInv2_4[nMomentum][nTheta];
    TF1 *sigInv4[nMomentum][nTheta];
    TF1 *expoInv_4[nMomentum][nTheta];

    TCanvas *C1[nMomentum][nTheta];
    TCanvas *C2[nMomentum][nTheta];
 

    for (int i = 0; i < nMomentum; i++)
    {
        for (int j = 0; j < nTheta; j++)
        {
            C1[i][j] = new TCanvas(Form("C1_%i_%i", i, j), Form("C1_%i_%i", i, j), 1200, 800);
            C1[i][j]->Divide(2, 2);

            C2[i][j] = new TCanvas(Form("C2_%i_%i", i, j), Form("C2_%i_%i", i, j), 1200, 800);
            C2[i][j]->Divide(2, 2);

            histOmegaRecoilAfterCuts[i][j] = (TH1F *)_file->Get(Form("Omega14MassCut7_%i_%i", i, j));
            histOmegaRecoilAfterCuts[i][j]->Rebin(4);
            double Max = histOmegaRecoilAfterCuts[i][j]->GetMaximum();

            histOmegaInvMassAfterCuts[i][j] = (TH1F *)_file->Get(Form("Omega14InvariantMassCut7_%i_%i", i, j));
            histOmegaInvMassAfterCuts[i][j]->Rebin(4);
            double MaxInv = histOmegaInvMassAfterCuts[i][j]->GetMaximum();
            
            // Fit Function gaus + gaus + expo
            fitFunc1[i][j] = new TF1(Form("fitFunc_%i_%i", i, j), "gaus(0) + gaus(3) + expo(6)", 0.65, 1.);
            fitFunc1[i][j]->SetParameters(0.6*Max, 0.782, 0.008, 0.3*Max, 0.782, 0.10, 0., -1.);
            fitFunc1[i][j]->SetParLimits(0, 0., Max);
            fitFunc1[i][j]->SetParLimits(1, 0.7, 0.85);
            fitFunc1[i][j]->SetParLimits(2, 0.005, 0.5);
            fitFunc1[i][j]->SetParLimits(3, 0., Max);
            fitFunc1[i][j]->SetParLimits(4, 0.7, 0.85);
            fitFunc1[i][j]->SetParLimits(5, 0.005, 0.5);
            histOmegaRecoilAfterCuts[i][j]->Fit(fitFunc1[i][j], "R0");

            gaus1_1[i][j] = new TF1(Form("gaus1_%i_%i", i, j), "gaus", 0.65, 1.);
            gaus1_1[i][j]->SetLineColor(kGreen);
            gaus1_1[i][j]->SetParameters(fitFunc1[i][j]->GetParameter(0), fitFunc1[i][j]->GetParameter(1), fitFunc1[i][j]->GetParameter(2));
            gaus2_1[i][j] = new TF1(Form("gaus2_%i_%i", i, j), "gaus", 0.65, 1.);
            gaus2_1[i][j]->SetLineColor(kBlue);
            gaus2_1[i][j]->SetParameters(fitFunc1[i][j]->GetParameter(3), fitFunc1[i][j]->GetParameter(4), fitFunc1[i][j]->GetParameter(5));
            expo_1[i][j] = new TF1(Form("expo_%i_%i", i, j), "expo", 0.65, 1.);
            expo_1[i][j]->SetLineColor(kBlack);
            expo_1[i][j]->SetParameters(fitFunc1[i][j]->GetParameter(6), fitFunc1[i][j]->GetParameter(7));
            
            fitFuncInv1[i][j] = new TF1(Form("fitFuncInv_%i_%i", i, j), "gaus(0) + gaus(3) + expo(6)", 0.65, 1.);
            fitFuncInv1[i][j]->SetParameters(0.6*MaxInv, 0.782, 0.008, 0.3*MaxInv, 0.782, 0.10, 0., -1.);
            fitFuncInv1[i][j]->SetParLimits(0, 0., MaxInv);
            fitFuncInv1[i][j]->SetParLimits(1, 0.7, 0.85);
            fitFuncInv1[i][j]->SetParLimits(2, 0.005, 0.5);
            fitFuncInv1[i][j]->SetParLimits(3, 0., MaxInv);
            fitFuncInv1[i][j]->SetParLimits(4, 0.7, 0.85);
            fitFuncInv1[i][j]->SetParLimits(5, 0.005, 0.5);
            histOmegaInvMassAfterCuts[i][j]->Fit(fitFuncInv1[i][j], "R0");
            gausInv1_1[i][j] = new TF1(Form("gausInv1_%i_%i", i, j), "gaus", 0.65, 1.);
            gausInv1_1[i][j]->SetLineColor(kGreen);
            gausInv1_1[i][j]->SetParameters(fitFuncInv1[i][j]->GetParameter(0), fitFuncInv1[i][j]->GetParameter(1), fitFuncInv1[i][j]->GetParameter(2));
            gausInv2_1[i][j] = new TF1(Form("gausInv2_%i_%i", i, j), "gaus", 0.65, 1.);
            gausInv2_1[i][j]->SetLineColor(kBlue);
            gausInv2_1[i][j]->SetParameters(fitFuncInv1[i][j]->GetParameter(3), fitFuncInv1[i][j]->GetParameter(4), fitFuncInv1[i][j]->GetParameter(5));
            expoInv_1[i][j] = new TF1(Form("expoInv_%i_%i", i, j), "expo", 0.65, 1.);
            expoInv_1[i][j]->SetLineColor(kBlack);
            expoInv_1[i][j]->SetParameters(fitFuncInv1[i][j]->GetParameter(6), fitFuncInv1[i][j]->GetParameter(7));

            C1[i][j]->cd(1);
            TLatex *lat1 = new TLatex(0.3, 0.6*Max, "gaus + gaus + expo");
            histOmegaRecoilAfterCuts[i][j]->Draw("HIST");
            fitFunc1[i][j]->Draw("SAME");
            gaus1_1[i][j]->Draw("SAME");
            gaus2_1[i][j]->Draw("SAME");
            expo_1[i][j]->Draw("SAME");
            C1[i][j]->cd(1)->Update();
            TPaveStats *stats1 = (TPaveStats*)histOmegaRecoilAfterCuts[i][j]->FindObject("stats")->Clone();
            //stats1->SetStats(0);
            stats1->Draw("SAME");
            lat1->Draw("SAME");

            C2[i][j]->cd(1);
            TLatex *latInv1 = new TLatex(0.3, 0.6*MaxInv, "gaus + gaus + expo");
            histOmegaInvMassAfterCuts[i][j]->Draw("HIST");
            fitFuncInv1[i][j]->Draw("SAME");
            gausInv1_1[i][j]->Draw("SAME");
            gausInv2_1[i][j]->Draw("SAME");
            expoInv_1[i][j]->Draw("SAME");
            C2[i][j]->cd(1)->Update();
            TPaveStats *statsInv1 = (TPaveStats*)histOmegaInvMassAfterCuts[i][j]->FindObject("stats")->Clone();
            //statsInv1->SetStats(0);
            statsInv1->Draw("SAME");
            latInv1->Draw("SAME");

            // Fit Function gaus + gaus + pol2
            fitFunc2[i][j] = new TF1(Form("fitFunc_%i_%i", i, j), "gaus(0) + gaus(3) + pol2(6)", 0.65, 1.);
            fitFunc2[i][j]->SetParameters(0.6*Max, 0.782, 0.008, 0.3*Max, 0.782, 0.10, 0., 0., 0.);
            fitFunc2[i][j]->SetParLimits(0, 0., Max);
            fitFunc2[i][j]->SetParLimits(1, 0.7, 0.85);
            fitFunc2[i][j]->SetParLimits(2, 0.005, 0.5);
            fitFunc2[i][j]->SetParLimits(3, 0., Max);
            fitFunc2[i][j]->SetParLimits(4, 0.7, 0.85);
            fitFunc2[i][j]->SetParLimits(5, 0.005, 0.5);
            histOmegaRecoilAfterCuts[i][j]->Fit(fitFunc2[i][j], "R0");
            gaus1_2[i][j] = new TF1(Form("gaus1_%i_%i", i, j), "gaus", 0.65, 1.);
            gaus1_2[i][j]->SetLineColor(kGreen);
            gaus1_2[i][j]->SetParameters(fitFunc2[i][j]->GetParameter(0), fitFunc2[i][j]->GetParameter(1), fitFunc2[i][j]->GetParameter(2));
            gaus2_2[i][j] = new TF1(Form("gaus2_%i_%i", i, j), "gaus", 0.65, 1.);
            gaus2_2[i][j]->SetLineColor(kBlue);
            gaus2_2[i][j]->SetParameters(fitFunc2[i][j]->GetParameter(3), fitFunc2[i][j]->GetParameter(4), fitFunc2[i][j]->GetParameter(5));
            pol2_2[i][j] = new TF1(Form("pol2_%i_%i", i, j), "pol2", 0.65, 1.);
            pol2_2[i][j]->SetLineColor(kBlack);
            pol2_2[i][j]->SetParameters(fitFunc2[i][j]->GetParameter(6), fitFunc2[i][j]->GetParameter(7), fitFunc2[i][j]->GetParameter(8));

            fitFuncInv2[i][j] = new TF1(Form("fitFuncInv_%i_%i", i, j), "gaus(0) + gaus(3) + pol2(6)", 0.65, 1.);
            fitFuncInv2[i][j]->SetParameters(0.6*MaxInv, 0.782, 0.008, 0.3*MaxInv, 0.782, 0.10, 0., 0., 0.);
            fitFuncInv2[i][j]->SetParLimits(0, 0., MaxInv);
            fitFuncInv2[i][j]->SetParLimits(1, 0.7, 0.85);
            fitFuncInv2[i][j]->SetParLimits(2, 0.005, 0.5);
            fitFuncInv2[i][j]->SetParLimits(3, 0., MaxInv);
            fitFuncInv2[i][j]->SetParLimits(4, 0.7, 0.85);
            fitFuncInv2[i][j]->SetParLimits(5, 0.005, 0.5);
            histOmegaInvMassAfterCuts[i][j]->Fit(fitFuncInv2[i][j], "R0");
            gausInv1_2[i][j] = new TF1(Form("gausInv1_%i_%i", i, j), "gaus", 0.65, 1.);
            gausInv1_2[i][j]->SetLineColor(kGreen);
            gausInv1_2[i][j]->SetParameters(fitFuncInv2[i][j]->GetParameter(0), fitFuncInv2[i][j]->GetParameter(1), fitFuncInv2[i][j]->GetParameter(2));
            gausInv2_2[i][j] = new TF1(Form("gausInv2_%i_%i", i, j), "gaus", 0.65, 1.);
            gausInv2_2[i][j]->SetLineColor(kBlue);
            gausInv2_2[i][j]->SetParameters(fitFuncInv2[i][j]->GetParameter(3), fitFuncInv2[i][j]->GetParameter(4), fitFuncInv2[i][j]->GetParameter(5));
            polInv2_2[i][j] = new TF1(Form("polInv2_%i_%i", i, j), "pol2", 0.65, 1.);
            polInv2_2[i][j]->SetLineColor(kBlack);
            polInv2_2[i][j]->SetParameters(fitFuncInv2[i][j]->GetParameter(6), fitFuncInv2[i][j]->GetParameter(7), fitFuncInv2[i][j]->GetParameter(8));
            
            C1[i][j]->cd(2);
            TLatex *lat2 = new TLatex(0.3, 0.6*Max, "gaus + gaus + pol2");
            histOmegaRecoilAfterCuts[i][j]->Draw("HIST");
            fitFunc2[i][j]->Draw("SAME");
            gaus1_2[i][j]->Draw("SAME");
            gaus2_2[i][j]->Draw("SAME");
            pol2_2[i][j]->Draw("SAME");
            C1[i][j]->cd(2)->Update();
            TPaveStats *stats2 = (TPaveStats*)histOmegaRecoilAfterCuts[i][j]->FindObject("stats")->Clone();
            //stats2->SetName("gaus + gaus + pol2");
            stats2->Draw("SAME");
            lat2->Draw("SAME");

            C2[i][j]->cd(2);
            TLatex *latInv2 = new TLatex(0.3, 0.6*MaxInv, "gaus + gaus + pol2");
            histOmegaInvMassAfterCuts[i][j]->Draw("HIST");
            fitFuncInv2[i][j]->Draw("SAME");
            gausInv1_2[i][j]->Draw("SAME");
            gausInv2_2[i][j]->Draw("SAME");
            polInv2_2[i][j]->Draw("SAME");
            C2[i][j]->cd(2)->Update();
            TPaveStats *statsInv2 = (TPaveStats*)histOmegaInvMassAfterCuts[i][j]->FindObject("stats")->Clone();
            //statsInv2->SetName("gaus + gaus + pol2");
            statsInv2->Draw("SAME");
            latInv2->Draw("SAME");

            // Fit Function gaus + gaus + pol1
            fitFunc3[i][j] = new TF1(Form("fitFunc_%i_%i", i, j), "gaus(0) + gaus(3) + pol1(6)", 0.65, 1.);
            fitFunc3[i][j]->SetParameters(0.6*Max, 0.782, 0.008, 0.3*Max, 0.782, 0.10, 0., 0.);
            fitFunc3[i][j]->SetParLimits(0, 0., Max);
            fitFunc3[i][j]->SetParLimits(1, 0.7, 0.85);
            fitFunc3[i][j]->SetParLimits(2, 0.005, 0.5);
            fitFunc3[i][j]->SetParLimits(3, 0., Max);
            fitFunc3[i][j]->SetParLimits(4, 0.7, 0.85);
            fitFunc3[i][j]->SetParLimits(5, 0.005, 0.5);
            histOmegaRecoilAfterCuts[i][j]->Fit(fitFunc3[i][j], "R0");
            gaus1_3[i][j] = new TF1(Form("gaus1_%i_%i", i, j), "gaus", 0.65, 1.);
            gaus1_3[i][j]->SetLineColor(kGreen);
            gaus1_3[i][j]->SetParameters(fitFunc3[i][j]->GetParameter(0), fitFunc3[i][j]->GetParameter(1), fitFunc3[i][j]->GetParameter(2));
            gaus2_3[i][j] = new TF1(Form("gaus2_%i_%i", i, j), "gaus", 0.65, 1.);
            gaus2_3[i][j]->SetLineColor(kBlue);
            gaus2_3[i][j]->SetParameters(fitFunc3[i][j]->GetParameter(3), fitFunc3[i][j]->GetParameter(4), fitFunc3[i][j]->GetParameter(5));
            pol1_3[i][j] = new TF1(Form("pol1_%i_%i", i, j), "pol1", 0.65, 1.);
            pol1_3[i][j]->SetLineColor(kBlack);
            pol1_3[i][j]->SetParameters(fitFunc3[i][j]->GetParameter(6), fitFunc3[i][j]->GetParameter(7));

            fitFuncInv3[i][j] = new TF1(Form("fitFuncInv_%i_%i", i, j), "gaus(0) + gaus(3) + pol1(6)", 0.65, 1.);
            fitFuncInv3[i][j]->SetParameters(0.6*MaxInv, 0.782, 0.008, 0.3*MaxInv, 0.782, 0.10, 0., 0.);
            fitFuncInv3[i][j]->SetParLimits(0, 0., MaxInv);
            fitFuncInv3[i][j]->SetParLimits(1, 0.7, 0.85);
            fitFuncInv3[i][j]->SetParLimits(2, 0.005, 0.5);
            fitFuncInv3[i][j]->SetParLimits(3, 0., MaxInv);
            fitFuncInv3[i][j]->SetParLimits(4, 0.7, 0.85);
            fitFuncInv3[i][j]->SetParLimits(5, 0.005, 0.5);
            histOmegaInvMassAfterCuts[i][j]->Fit(fitFuncInv3[i][j], "R0");
            gausInv1_3[i][j] = new TF1(Form("gausInv1_%i_%i", i, j), "gaus", 0.65, 1.);
            gausInv1_3[i][j]->SetLineColor(kGreen);
            gausInv1_3[i][j]->SetParameters(fitFuncInv3[i][j]->GetParameter(0), fitFuncInv3[i][j]->GetParameter(1), fitFuncInv3[i][j]->GetParameter(2));
            gausInv2_3[i][j] = new TF1(Form("gausInv2_%i_%i", i, j), "gaus", 0.65, 1.);
            gausInv2_3[i][j]->SetLineColor(kBlue);
            gausInv2_3[i][j]->SetParameters(fitFuncInv3[i][j]->GetParameter(3), fitFuncInv3[i][j]->GetParameter(4), fitFuncInv3[i][j]->GetParameter(5));
            polInv1_3[i][j] = new TF1(Form("polInv1_%i_%i", i, j), "pol1", 0.65, 1.);
            polInv1_3[i][j]->SetLineColor(kBlack);
            polInv1_3[i][j]->SetParameters(fitFuncInv3[i][j]->GetParameter(6), fitFuncInv3[i][j]->GetParameter(7));

            C1[i][j]->cd(3);
            TLatex *lat3 = new TLatex(0.3, 0.6*Max, "gaus + gaus + pol1");
            histOmegaRecoilAfterCuts[i][j]->Draw("HIST");
            fitFunc3[i][j]->Draw("SAME");
            gaus1_3[i][j]->Draw("SAME");
            gaus2_3[i][j]->Draw("SAME");
            pol1_3[i][j]->Draw("SAME");
            C1[i][j]->cd(3)->Update();
            TPaveStats *stats3 = (TPaveStats*)histOmegaRecoilAfterCuts[i][j]->FindObject("stats")->Clone();
            //stats3->SetName("gaus + gaus + pol1");
            stats3->Draw("SAME");
            lat3->Draw("SAME");

            C2[i][j]->cd(3);
            TLatex *latInv3 = new TLatex(0.3, 0.6*MaxInv, "gaus + gaus + pol1");
            histOmegaInvMassAfterCuts[i][j]->Draw("HIST");
            fitFuncInv3[i][j]->Draw("SAME");
            gausInv1_3[i][j]->Draw("SAME");
            gausInv2_3[i][j]->Draw("SAME");
            polInv1_3[i][j]->Draw("SAME");
            C2[i][j]->cd(3)->Update();
            TPaveStats *statsInv3 = (TPaveStats*)histOmegaInvMassAfterCuts[i][j]->FindObject("stats")->Clone();
            //statsInv3->SetName("gaus + gaus + pol1");
            statsInv3->Draw("SAME");
            latInv3->Draw("SAME");

            // Fit Function gaus + gaus + exp(x**2)
            fitFunc4[i][j] = new TF1(Form("fitFunc_%i_%i", i, j), FITFUNC4, 0.65, 1., 9);
            fitFunc4[i][j]->SetParameters(0.6*Max, 0.782, 0.008, 0.3*Max, 0.782, 0.10, 0., 0., 0.);
            fitFunc4[i][j]->SetParLimits(0, 0., Max);
            fitFunc4[i][j]->SetParLimits(1, 0.7, 0.85);
            fitFunc4[i][j]->SetParLimits(2, 0.005, 0.5);
            fitFunc4[i][j]->SetParLimits(3, 0., Max);
            fitFunc4[i][j]->SetParLimits(4, 0.7, 0.85);
            fitFunc4[i][j]->SetParLimits(5, 0.005, 1.0);
            histOmegaRecoilAfterCuts[i][j]->Fit(fitFunc4[i][j], "R0");
            gaus1_4[i][j] = new TF1(Form("gaus1_%i_%i", i, j), gaus, 0.65, 1., 3);
            gaus1_4[i][j]->SetLineColor(kGreen);
            gaus1_4[i][j]->SetParameters(fitFunc4[i][j]->GetParameter(0), fitFunc4[i][j]->GetParameter(1), fitFunc4[i][j]->GetParameter(2));
            gaus2_4[i][j] = new TF1(Form("gaus2_%i_%i", i, j), gaus, 0.65, 1., 3);
            gaus2_4[i][j]->SetLineColor(kBlue);
            gaus2_4[i][j]->SetParameters(fitFunc4[i][j]->GetParameter(3), fitFunc4[i][j]->GetParameter(4), fitFunc4[i][j]->GetParameter(5));
            expo_4[i][j] = new TF1(Form("expo_%i_%i", i, j), expo2, 0.65, 1., 3);
            expo_4[i][j]->SetLineColor(kBlack);
            expo_4[i][j]->SetParameters(fitFunc4[i][j]->GetParameter(6), fitFunc4[i][j]->GetParameter(7), fitFunc4[i][j]->GetParameter(8));

            /*
            sig4[i][j] = new TF1(Form("sig_%i_%i", i, j), doublegaus, 0.65, 1., 6);
            sig4[i][j]->SetParameters(fitFunc4[i][j]->GetParameter(0), fitFunc4[i][j]->GetParameter(1), 
                                        fitFunc4[i][j]->GetParameter(2), fitFunc4[i][j]->GetParameter(3), 
                                        fitFunc4[i][j]->GetParameter(4), fitFunc4[i][j]->GetParameter(5)
                                    );
            */
            sig4[i][j] = new TF1(Form("sig_%i_%i", i, j), "gaus", 0.65, 1.);
            sig4[i][j]->SetParameters(fitFunc4[i][j]->GetParameter(0), fitFunc4[i][j]->GetParameter(1), fitFunc4[i][j]->GetParameter(2));
            YieldRecoil[i][j] = sig4[i][j]->Integral(0.65, 1.);
            YieldRecoilErr[i][j] = (fitFunc4[i][j]->GetChisquare()/fitFunc4[i][j]->GetNDF())*sig4[i][j]->IntegralError(0.65, 1.);

            BkgRecoil[i][j] = expo_4[i][j]->Integral(0.65, 1.);
            BkgRecoilErr[i][j] = expo_4[i][j]->IntegralError(0.65, 1.);

            fitFuncInv4[i][j] = new TF1(Form("fitFuncInv_%i_%i", i, j), FITFUNC4, 0.65, 1., 9);
            fitFuncInv4[i][j]->SetParameters(0.6*MaxInv, 0.782, 0.008, 0.3*MaxInv, 0.782, 0.10, 0., 0., 0.);
            fitFuncInv4[i][j]->SetParLimits(0, 0., MaxInv);
            fitFuncInv4[i][j]->SetParLimits(1, 0.7, 0.85);
            fitFuncInv4[i][j]->SetParLimits(2, 0.005, 0.5);
            fitFuncInv4[i][j]->SetParLimits(3, 0., MaxInv);
            fitFuncInv4[i][j]->SetParLimits(4, 0.7, 0.85);
            fitFuncInv4[i][j]->SetParLimits(5, 0.005, 1.00);
            histOmegaInvMassAfterCuts[i][j]->Fit(fitFuncInv4[i][j], "R0");
            gausInv1_4[i][j] = new TF1(Form("gausInv1_%i_%i", i, j), gaus, 0.65, 1., 3);
            gausInv1_4[i][j]->SetLineColor(kGreen);
            gausInv1_4[i][j]->SetParameters(fitFuncInv4[i][j]->GetParameter(0), fitFuncInv4[i][j]->GetParameter(1), fitFuncInv4[i][j]->GetParameter(2));
            gausInv2_4[i][j] = new TF1(Form("gausInv2_%i_%i", i, j), gaus, 0.65, 1., 3);
            gausInv2_4[i][j]->SetLineColor(kBlue);
            gausInv2_4[i][j]->SetParameters(fitFuncInv4[i][j]->GetParameter(3), fitFuncInv4[i][j]->GetParameter(4), fitFuncInv4[i][j]->GetParameter(5));
            expoInv_4[i][j] = new TF1(Form("expoInv_%i_%i", i, j), expo2, 0.65, 1., 3);
            expoInv_4[i][j]->SetLineColor(kBlack);
            expoInv_4[i][j]->SetParameters(fitFuncInv4[i][j]->GetParameter(6), fitFuncInv4[i][j]->GetParameter(7), fitFuncInv4[i][j]->GetParameter(8));

            /*
            sigInv4[i][j] = new TF1(Form("sigInv_%i_%i", i, j), doublegaus, 0.65, 1., 6);
            sigInv4[i][j]->SetParameters(fitFuncInv4[i][j]->GetParameter(0), fitFuncInv4[i][j]->GetParameter(1), 
                                        fitFuncInv4[i][j]->GetParameter(2), fitFuncInv4[i][j]->GetParameter(3), 
                                        fitFuncInv4[i][j]->GetParameter(4), fitFuncInv4[i][j]->GetParameter(5)
                                    );
            */
            sigInv4[i][j] = new TF1(Form("sigInv_%i_%i", i, j), "gaus", 0.65, 1.);
            sigInv4[i][j]->SetParameters(fitFuncInv4[i][j]->GetParameter(0), fitFuncInv4[i][j]->GetParameter(1), fitFuncInv4[i][j]->GetParameter(2));
            YieldInv[i][j] = sigInv4[i][j]->Integral(0.65, 1.);
            YieldInvErr[i][j] = (fitFuncInv4[i][j]->GetChisquare()/fitFuncInv4[i][j]->GetNDF())*sigInv4[i][j]->IntegralError(0.65, 1.);

            BkgInv[i][j] = expoInv_4[i][j]->Integral(0.65, 1.);
            BkgInvErr[i][j] = expoInv_4[i][j]->IntegralError(0.65, 1.);

            Effic[i][j] = YieldInv[i][j]/(YieldInv[i][j] + YieldRecoil[i][j]);
            EfficErr[i][j] = Effic[i][j]*sqrt(pow(YieldInvErr[i][j]/YieldInv[i][j], 2) + pow(YieldRecoilErr[i][j]/YieldRecoil[i][j], 2));
            
            C1[i][j]->cd(4);
            TLatex *lat4 = new TLatex(0.3, 0.6*Max, "gaus + gaus + exp(ax**2 + bx + c)");
            histOmegaRecoilAfterCuts[i][j]->Draw("HIST");
            fitFunc4[i][j]->Draw("SAME");
            gaus1_4[i][j]->Draw("SAME");
            gaus2_4[i][j]->Draw("SAME");
            expo_4[i][j]->Draw("SAME");
            C1[i][j]->cd(4)->Update();
            TPaveStats *stats4 = (TPaveStats*)histOmegaRecoilAfterCuts[i][j]->FindObject("stats")->Clone();
            //stats4->SetName("gaus + gaus + exp(ax**2 + bx + c)");
            stats4->Draw("SAME");
            lat4->Draw("SAME");

            C2[i][j]->cd(4);
            TLatex *latInv4 = new TLatex(0.3, 0.6*MaxInv, "gaus + gaus + exp(ax**2 + bx + c)");
            histOmegaInvMassAfterCuts[i][j]->Draw("HIST");
            fitFuncInv4[i][j]->Draw("SAME");
            gausInv1_4[i][j]->Draw("SAME");
            gausInv2_4[i][j]->Draw("SAME");
            expoInv_4[i][j]->Draw("SAME");
            C2[i][j]->cd(4)->Update();
            TPaveStats *statsInv4 = (TPaveStats*)histOmegaInvMassAfterCuts[i][j]->FindObject("stats")->Clone();
            //statsInv4->SetName("gaus + gaus + exp(ax**2 + bx + c)");
            statsInv4->Draw("SAME");
            latInv4->Draw("SAME");
            
            if(i == 0 && j == 0) C1[i][j]->Print("RecoilOmegaMassFit.pdf(");
            else if (i == nMomentum - 1 && j == nTheta - 1) C1[i][j]->Print("RecoilOmegaMassFit.pdf)");
            else C1[i][j]->Print("RecoilOmegaMassFit.pdf");

            if(i == 0 && j == 0) C2[i][j]->Print("InvMassOmegaMassFit.pdf(");
            else if (i == nMomentum - 1 && j == nTheta - 1) C2[i][j]->Print("InvMassOmegaMassFit.pdf)");
            else C2[i][j]->Print("InvMassOmegaMassFit.pdf");

            TGraphErrors *grErr[nMomentum];
            const int nColors[nMomentum] = {kRed, kBlue, kMagenta};
            TLegend *leg = new TLegend(0.2, 0.1, 0.6, 0.4);
            for (int p = 0; p < nMomentum; p++)
            {
                grErr[p] = new TGraphErrors(nTheta);
                grErr[p]->SetName(Form("Efficiency_p_%i", p));
                grErr[p]->SetMarkerStyle(20);
                grErr[p]->SetMarkerSize(1.5);
                grErr[p]->SetMarkerColor(nColors[p]);
                grErr[p]->SetLineColor(kBlack);
                grErr[p]->SetLineWidth(2);
                grErr[p]->SetTitle("Efficiency vs. #theta; #theta [deg]; Efficiency");
                leg->AddEntry(grErr[p], Form("%.2f #leq ~ p < %.2f [GeV/c]", lowMomBins[p], highMomBins[p]), "lep");
                for (int t = 0; t < nTheta; t++)
                {
                    grErr[p]->SetPoint(t, lowthetaBins[t] + (highthetaBins[t] - lowthetaBins[t])/2., Effic[p][t]);
                    grErr[p]->SetPointError(t, (highthetaBins[t] - lowthetaBins[t])/2., EfficErr[p][t]);
                }
                grErr[p]->SetMinimum(0.);
            }
            TCanvas *C3 = new TCanvas ("Efficiency", "Efficiency", 1200, 800);
            C3->cd();
            for (int p = 0; p < nMomentum; p++)
            {
                if (p == 0) grErr[p]->Draw("AP");
                //else grErr[p]->Draw("P SAME");
            }
            //leg->Draw("SAME");
            C3->Print("Efficiency.pdf");
            

        }
    }

}
