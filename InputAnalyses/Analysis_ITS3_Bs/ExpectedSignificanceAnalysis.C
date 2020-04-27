void InitializeProbCuts();

void fit_signal(TH1F* hsig, int j, int i, double &xmin, double &xmax, bool draw);
void calculate_background(TH1F* hbkg, TF1* f1, int j, int i, double xmin, double xmax, double &bkgcent, TString filenamefits, bool draw, bool finalscan, double &bkglow, double &bkghigh);
void extract_fonll(TString filnam, int j, double &fonllmin, double &fonllcent, double &fonllmax, bool draw);
void extract_TAMU(TString filn, int j, double &tamucent, bool draw);
void calculate_efficiency(TH1F* heff, int j, int i, double &effcent, double &erreffcent, bool draw);

void plot_expected_significance(TString filenamemass, TString filenameeff, TString filenamefits);

const Int_t nptbins = 7;
Int_t ptbins[nptbins+1] = {0, 2, 4, 6, 8, 12, 16, 24};
Float_t ptbinsfl[nptbins+1] = {0, 2, 4, 6, 8, 12, 16, 24};

const Int_t trials = 5000;
Float_t probcuts[trials+1];
Int_t selbin[nptbins] = {4988, 4983, 4995, 4978, 4945, 4908, 4956};

Double_t fitsplit = 0.0001;
Int_t nfits = 8;

Int_t rebin[nptbins] = {16, 16, 16, 16, 16, 16, 16};
Bool_t bincountBkg[nptbins] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
Int_t sgnfmax = 10;

Double_t nEv335 = 852644; //756420
TString nEv335String = "8.53 * 10^{5}"; //"7.56 * 10^{5}";
Double_t nEvExpected = 8000000000;
TString nEvExpectedString = "8 * 10^{9}";
TString nLumiExpectedString = "10 nb^{-1}";

TString filenameBkgCorr = "theory/BkgCorrFactor_Bs_1DataFile_25MCFile.root";
TString filenameTAMU = "theory/input_RAA_TAMU_Bs.txt";
TString filnameFONLL = "theory/inputFONLL.txt";
Double_t fbtoB = 0.407; //http://pdg.lbl.gov/2019/reviews/rpp2018-rev-b-meson-prod-decay.pdf, table 85.1
Double_t fbtoBUnc = 0.007;
Double_t fLHCbBBs = 2*0.122; //https://journals.aps.org/prd/pdf/10.1103/PhysRevD.100.031102, factor 2 because B0 + B+
Double_t fLHCbBBsUnc = 2*0.006;
Double_t relSystFF = fbtoBUnc/fbtoB;
Double_t relSystLHCb = fLHCbBBsUnc/fLHCbBBs;

Double_t BRBs = 0.00304;
Double_t errBRBs = 0.00023;
Double_t BRDs = 0.0227;
Double_t errBRDs = 0.0008;
Double_t TAA = 23.07 * 1e-3; //mb^-1 -> mub^-1 in which we put FONLL
Double_t gauss3sigmafactor = 0.9973;

void ExpectedSignificanceAnalysis(TString filenamemass, TString filenameeff, TString filenamefits){
  InitializeProbCuts();
  plot_expected_significance(filenamemass, filenameeff, filenamefits);
}

void plot_expected_significance(TString filenamemass, TString filenameeff, TString filenamefits){
  
  TGaxis::SetMaxDigits(3);
  
  TFile* fmass = new TFile(filenamemass.Data());
  TFile* feff = new TFile(filenameeff.Data());
  TFile* fBkgCorr = new TFile(filenameBkgCorr.Data());

  TString nameeff = "eff_";

  TH1F *heff[trials+1];

  TString namemasssig = "hmass_sigpt_cand";
  TString namemassbkg = "hmass_bkgpt_cand";

  TH1F *hsig[nptbins][trials+1];
  TH1F *hbkg[nptbins][trials+1];
  TF1  *fbkg[nptbins][trials+1];

  for(int i = 0; i <= trials; i++){
      
    if(i != 0 &&
       i != selbin[0] &&
       i != selbin[1] &&
       i != selbin[2] &&
       i != selbin[3] &&
       i != selbin[4] &&
       i != selbin[5] &&
       i != selbin[6] ) continue;
      
    heff[i] = (TH1F*)feff->Get(Form("%s%d",nameeff.Data(),i));
    heff[i]->GetBinContent(1);

    for(int j = 0; j < nptbins; j++){
      hsig[j][i] = (TH1F*)fmass->Get(Form("%s%d_%d_%d",namemasssig.Data(),ptbins[j],ptbins[j+1],i));
      hbkg[j][i] = (TH1F*)fmass->Get(Form("%s%d_%d_%d",namemassbkg.Data(),ptbins[j],ptbins[j+1],i));
      fbkg[j][i] = new TF1(Form("f_%d_%d",j,i), "expo",  5.07, 5.65);
    }
  }
  
  TH1F* hBkgCorr = (TH1F*)fBkgCorr->Get("hCorrFacBs");

  double xmin[nptbins];
  double xmax[nptbins];
  double fonllmin[nptbins];
  double fonllcent[nptbins];
  double fonllmax[nptbins];
  double tamucent[nptbins];
  double effcent[nptbins];
  double erreffcent[nptbins];
  double expectedsignalcent[nptbins];
  double expectedsignalfonllmin[nptbins];
  double expectedsignalfonllmax[nptbins];
  double expectedsignalerreffmin[nptbins];
  double expectedsignalerreffmax[nptbins];
  double expectedbkg[nptbins];
  double expectedbkglow[nptbins];
  double expectedbkghigh[nptbins];
  double expectedsgnfcent[nptbins];
  double expectedsgnffonllmin[nptbins];
  double expectedsgnffonllmax[nptbins];
  double expectedsgnferreffmin[nptbins];
  double expectedsgnferreffmax[nptbins];
  double expectedsgnfbkgmin[nptbins];
  double expectedsgnfbkgmax[nptbins];
  TCanvas *canv[nptbins];
  
  for(int j = 0; j < nptbins; j++){
    for(int i = 0; i <= trials; i++){
      if(j == 0 && i != selbin[0]) continue;
      if(j == 1 && i != selbin[1]) continue;
      if(j == 2 && i != selbin[2]) continue;
      if(j == 3 && i != selbin[3]) continue;
      if(j == 4 && i != selbin[4]) continue;
      if(j == 5 && i != selbin[5]) continue;
      if(j == 6 && i != selbin[6]) continue;

      canv[j] = new TCanvas(Form("canv%d_%d_%d",ptbins[j],ptbins[j+1],i),Form("canv%d_%d_%d",ptbins[j],ptbins[j+1],i),1050,600);
      canv[j]->Divide(3,2);

      canv[j]->cd(1);
      cout << "Using signal histogram = " << hsig[j][i]->GetName() << endl;
      fit_signal(hsig[j][i], j, i, xmin[j], xmax[j], kTRUE);

      canv[j]->cd(2);
      cout << "Using background histogram = " << hbkg[j][i]->GetName() << endl;
      calculate_background(hbkg[j][i], fbkg[j][i], j, i, xmin[j], xmax[j], expectedbkg[j], filenamefits, kTRUE, kTRUE, expectedbkglow[j], expectedbkghigh[j]);

      canv[j]->cd(3);
      cout << "Using input FONLL file = " << filnameFONLL << endl;
      extract_fonll(filnameFONLL, j, fonllmin[j], fonllcent[j], fonllmax[j], kTRUE);

      canv[j]->cd(4);
      cout << "Using input TAMU file = " << filenameTAMU << endl;
      extract_TAMU(filenameTAMU, j, tamucent[j], kTRUE);

      canv[j]->cd(5);
      cout << "Using input efficiency histograms = " << heff[i]->GetName() << endl;
      calculate_efficiency(heff[i], j, i, effcent[j], erreffcent[j], kTRUE);

      expectedsignalcent[j] = 2 * (ptbins[j+1] - ptbins[j]) * 1 * (BRBs * BRDs) * nEvExpected * effcent[j] * TAA * fonllcent[j] * tamucent[j];
      expectedsignalfonllmin[j] = 2 * (ptbins[j+1] - ptbins[j]) * 1 * (BRBs * BRDs) * nEvExpected * effcent[j] * TAA * fonllmin[j] * tamucent[j];
      expectedsignalfonllmax[j] = 2 * (ptbins[j+1] - ptbins[j]) * 1 * (BRBs * BRDs) * nEvExpected * effcent[j] * TAA * fonllmax[j] * tamucent[j];
      expectedsignalerreffmin[j] = 2 * (ptbins[j+1] - ptbins[j]) * 1 * (BRBs * BRDs) * nEvExpected * (effcent[j] - erreffcent[j]) * TAA * fonllcent[j] * tamucent[j];
      expectedsignalerreffmax[j] = 2 * (ptbins[j+1] - ptbins[j]) * 1 * (BRBs * BRDs) * nEvExpected * (effcent[j] + erreffcent[j]) * TAA * fonllcent[j] * tamucent[j];

      if(hBkgCorr->FindBin(ptbins[j]+0.01) != hBkgCorr->FindBin(ptbins[j+1]-0.01)){
        cout << "Warning! Different pT binning HIJING correction factor" << endl;
      }
      expectedbkg[j] = hBkgCorr->GetBinContent(hBkgCorr->FindBin(ptbins[j]+0.01)) * nEvExpected * expectedbkg[j] / nEv335;
      expectedbkglow[j] = hBkgCorr->GetBinContent(hBkgCorr->FindBin(ptbins[j]+0.01)) * nEvExpected * expectedbkglow[j] / nEv335;
      expectedbkghigh[j] = hBkgCorr->GetBinContent(hBkgCorr->FindBin(ptbins[j]+0.01)) * nEvExpected * expectedbkghigh[j] / nEv335;
      
      expectedsignalcent[j] *= gauss3sigmafactor;
      expectedsignalfonllmin[j] *= gauss3sigmafactor;
      expectedsignalfonllmax[j] *= gauss3sigmafactor;
      expectedsignalerreffmin[j] *= gauss3sigmafactor;
      expectedsignalerreffmax[j] *= gauss3sigmafactor;

      expectedsgnfcent[j] = expectedsignalcent[j] / TMath::Sqrt(expectedsignalcent[j] + expectedbkg[j]);
      expectedsgnffonllmin[j] = expectedsignalfonllmin[j] / TMath::Sqrt(expectedsignalfonllmin[j] + expectedbkg[j]);
      expectedsgnffonllmax[j] = expectedsignalfonllmax[j] / TMath::Sqrt(expectedsignalfonllmax[j] + expectedbkg[j]);
      expectedsgnferreffmin[j] = expectedsignalerreffmin[j] / TMath::Sqrt(expectedsignalerreffmin[j] + expectedbkg[j]);
      expectedsgnferreffmax[j] = expectedsignalerreffmax[j] / TMath::Sqrt(expectedsignalerreffmax[j] + expectedbkg[j]);
      expectedsgnfbkgmin[j] = expectedsignalerreffmin[j] / TMath::Sqrt(expectedsignalcent[j] + expectedbkghigh[j]);
      expectedsgnfbkgmax[j] = expectedsignalerreffmax[j] / TMath::Sqrt(expectedsignalcent[j] + expectedbkglow[j]);

      cout << probcuts[i] << " " << expectedsignalcent[j] << " " << expectedbkg[j] << " " << expectedsgnfcent[j] << endl;
    }
  }
  
  //Plot final distributions
  TH1F* hefftot = (TH1F*)heff[0]->Clone("hefftot");
  hefftot->SetTitle("Acceptance-times-efficiency;#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon) #times 2#it{y}_{fid}");
  hefftot->SetMarkerStyle(20);
  hefftot->SetMarkerColor(kAzure+1);
  hefftot->SetLineColor(kAzure+1);
  hefftot->SetLineWidth(2);
  hefftot->SetStats(0);

  TH1F* hSgnftot = new TH1F("hSngftot",Form("Expected Significance (%s);#it{p}_{T} (GeV/#it{c});Expected Significance",nLumiExpectedString.Data()),nptbins,ptbinsfl);
  hSgnftot->SetMarkerStyle(20);
  hSgnftot->SetLineColor(kBlack);
  hSgnftot->SetMarkerColor(kBlack);
  hSgnftot->SetLineWidth(2);
  hSgnftot->SetStats(0);
  
  TGraphAsymmErrors* grSgnftot = new TGraphAsymmErrors(0);
  grSgnftot->SetName("grSgnftot");
  grSgnftot->SetFillColor(kAzure+1);
  grSgnftot->SetLineColor(kAzure+1);
  grSgnftot->SetLineWidth(2);
  grSgnftot->SetFillStyle(3002);

  TGraphAsymmErrors* grSgnfBkgUnc = new TGraphAsymmErrors(0);
  grSgnfBkgUnc->SetName("grSgnfBkgUnc");
  grSgnfBkgUnc->SetFillColor(kBlue);
  grSgnfBkgUnc->SetLineColor(kBlue);
  grSgnfBkgUnc->SetLineWidth(2);
  grSgnfBkgUnc->SetFillStyle(0);
  
  for(int j = 0; j < nptbins; j++){
    for(int i = 0; i <= trials; i++){
      if(j == 0 && i != selbin[0]) continue;
      if(j == 1 && i != selbin[1]) continue;
      if(j == 2 && i != selbin[2]) continue;
      if(j == 3 && i != selbin[3]) continue;
      if(j == 4 && i != selbin[4]) continue;
      if(j == 5 && i != selbin[5]) continue;
      if(j == 6 && i != selbin[6]) continue;

      hefftot->SetBinContent(j+1, effcent[j]);
      hefftot->SetBinError(j+1, erreffcent[j]);
      
      hSgnftot->SetBinContent(j+1,expectedsgnfcent[j]);
      hSgnftot->SetBinError(j+1,expectedsgnfcent[j] - expectedsgnferreffmin[j]);
      
      grSgnftot->SetPoint(j, ptbins[j] + 0.5*(ptbins[j+1]-ptbins[j]), expectedsgnfcent[j]);
      grSgnftot->SetPointError(j, 0.25*(ptbins[j+1]-ptbins[j]), 0.25*(ptbins[j+1]-ptbins[j]), expectedsgnfcent[j] - expectedsgnffonllmin[j], expectedsgnffonllmax[j] - expectedsgnfcent[j]);
      
      grSgnfBkgUnc->SetPoint(j, ptbins[j] + 0.5*(ptbins[j+1]-ptbins[j]), expectedsgnfcent[j]);
      grSgnfBkgUnc->SetPointError(j,0, 0, expectedsgnfcent[j] - expectedsgnfbkgmin[j], expectedsgnfbkgmax[j] - expectedsgnfcent[j]);
    }
  }
  
  for(int j = 0; j < nptbins; j++){
    for(int i = 0; i <= trials; i++){
      if(j == 0 && i != selbin[0]) continue;
      if(j == 1 && i != selbin[1]) continue;
      if(j == 2 && i != selbin[2]) continue;
      if(j == 3 && i != selbin[3]) continue;
      if(j == 4 && i != selbin[4]) continue;
      if(j == 5 && i != selbin[5]) continue;
      if(j == 6 && i != selbin[6]) continue;

      canv[j]->cd(5);
      hefftot->Draw("ep");
      gPad->SetTickx();
      gPad->SetTicky();
      hefftot->GetYaxis()->SetRangeUser(0.001,1);
      gPad->SetLogy();

      TH1F* heffcent = new TH1F(Form("heffcent_set%d_%d_v2",i,j),"",1,ptbins[j],ptbins[j+1]);
      heffcent->SetBinContent(1, hefftot->GetBinContent(hefftot->FindBin(ptbins[j] + 0.5*(ptbins[j+1]-ptbins[j]))));
      heffcent->SetBinError(1, hefftot->GetBinError(hefftot->FindBin(ptbins[j] + 0.5*(ptbins[j+1]-ptbins[j]))));
      heffcent->SetMarkerStyle(20);
      heffcent->SetMarkerColor(kRed);
      heffcent->SetLineColor(kRed);
      heffcent->SetLineWidth(2);
      heffcent->Draw("same ep");

      canv[j]->cd(6);
      hSgnftot->Draw("ep");
      gPad->SetTickx();
      gPad->SetTicky();
      hSgnftot->GetYaxis()->SetRangeUser(0.,sgnfmax);
      grSgnftot->Draw("same e2");
      hSgnftot->Draw("same ep");
      grSgnfBkgUnc->Draw("same []");
      TLegend* leg = new TLegend(0.2, 0.6, 0.68, 0.8, 0, "NDC");
      leg->SetTextFont(43); leg->SetTextSize(14); leg->SetFillColor(0); leg->SetLineColor(0);
      leg->AddEntry(hSgnftot,"Stat. unc. from efficiency", "lpm");
      leg->AddEntry(grSgnftot,"FONLL and FF unc.","f");
      leg->AddEntry(grSgnfBkgUnc,"Background estimation unc.","e2");
      leg->Draw();
    }
  }
  
  for(int j = 0; j < nptbins; j++){
    TFile* ffits = new TFile(filenamefits.Data());
    TCanvas* cfits = (TCanvas*)ffits->Get(Form("BkgFits_parametrised_%d",j));
    //gStyle->SetOptFit(0);
    //gStyle->SetOptStat(0);
    cfits->Draw();
  }
  /*
  TFile* fSgnf = new TFile("ExpectedSignificance_ITS23_100320.root","UPDATE");
  fSgnf->cd();
  hSgnftot->Write(Form("hSgnftot_%d",selbin[0]));
  grSgnftot->Write(Form("grSgnftot_%d",selbin[0]));
  grSgnfBkgUnc->Write(Form("grSgnfBkgUnc_%d",selbin[0]));
  fSgnf->Close();
  */
}

void fit_signal(TH1F* hsig, int j, int i, double &xmin, double &xmax, bool draw){

  TF1* f1 = new TF1("f1", "gaus",  5.32, 5.42);
  hsig->Fit("f1", "R");

  if(draw){
    hsig->SetTitle(Form("%d < #it{p}_{T} < %d GeV/#it{c};#it{M}_{ inv} B_{s} (GeV/#it{c}^{2});Counts",ptbins[j], ptbins[j+1]));
    hsig->SetName(Form("hsig_set%d_pt%d", i, j+1));
    hsig->SetLineColor(kBlack);
    
    gStyle->SetOptFit(1111);
    hsig->Draw("hist");
    gPad->SetTickx();
    gPad->SetTicky();
    f1->Draw("same");
  }
  
  double mean = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);

  xmin = mean - 3 * sigma;
  xmax = mean + 3 * sigma;
}

void calculate_background(TH1F* hbkg, TF1* f1, int j, int i, double xmin, double xmax, double &bkgcent, TString filenamefits, bool draw, bool finalscan, double &bkglow, double &bkghigh){

  hbkg->Rebin(rebin[j]);
  
  if(draw){
    hbkg->SetTitle(Form("%d < #it{p}_{T} < %d GeV/#it{c};#it{M}_{ inv} B_{s} (GeV/#it{c}^{2});Counts",ptbins[j], ptbins[j+1]));
    hbkg->SetName(Form("hbkg_set%d_pt%d", i, j+1));
    hbkg->SetStats(0);
    hbkg->SetLineColor(kBlack);
    hbkg->SetMarkerColor(kBlack);
    hbkg->SetMarkerStyle(20);
    
    gStyle->SetOptFit(1111);
    hbkg->Draw("ep");
  }
  
  //=Chi2 fit, add "L" for likelyhood
  hbkg->Fit(Form("f_%d_%d",j,i), "R,E,+");

  TF1* fcent;
  TF1* flow;
  TF1* fup;
  if(draw || finalscan){
    TFile* ffits = new TFile(filenamefits.Data());
    cout << filenamefits << endl;
    TGraph* grpar0cent = (TGraph*)ffits->Get(Form("gpar0cent_%d",j));
    TGraph* grpar0low = (TGraph*)ffits->Get(Form("gpar0low_%d",j));
    TGraph* grpar0up = (TGraph*)ffits->Get(Form("gpar0high_%d",j));
    TF1* fpar1cent = (TF1*)ffits->Get(Form("fpar1cent_%d",j));
    TF1* fpar1low = (TF1*)ffits->Get(Form("fpar1low_%d",j));
    TF1* fpar1up = (TF1*)ffits->Get(Form("fpar1high_%d",j));

    fcent = (TF1*)f1->Clone(Form("fcent_%d",j));
    fcent->SetParameter(0,grpar0cent->Eval(probcuts[i]));
    double dpar1cent = fpar1cent->Eval(probcuts[i]);
    if(dpar1cent > 0) dpar1cent = 0;
    fcent->SetParameter(1,dpar1cent);
    fcent->SetLineColor(kBlue);
    flow = (TF1*)f1->Clone(Form("flow_%d",j));
    flow->SetParameter(0,grpar0low->Eval(probcuts[i]));
    double dpar1low = fpar1low->Eval(probcuts[i]);
    if(dpar1low > 0) dpar1low = 0;
    flow->SetParameter(1,dpar1low);
    flow->SetLineStyle(2);
    flow->SetLineColor(kBlue);
    fup = (TF1*)f1->Clone(Form("fup_%d",j));
    fup->SetParameter(0,grpar0up->Eval(probcuts[i]));
    double dpar1up = fpar1up->Eval(probcuts[i]);
    if(dpar1up > 0) dpar1up = 0;
    fup->SetParameter(1,dpar1up);
    fup->SetLineStyle(2);
    fup->SetLineColor(kBlue);
  }
  
  Double_t bkgcentbc = hbkg->Integral(hbkg->FindBin(xmin),hbkg->FindBin(xmax));
  Double_t bkglowbc = hbkg->Integral(hbkg->FindBin(xmin),hbkg->FindBin(xmax)) - TMath::Sqrt(hbkg->Integral(hbkg->FindBin(xmin),hbkg->FindBin(xmax)));
  Double_t bkghighbc = hbkg->Integral(hbkg->FindBin(xmin),hbkg->FindBin(xmax)) + TMath::Sqrt(hbkg->Integral(hbkg->FindBin(xmin),hbkg->FindBin(xmax)));
  
  Double_t bkgcentfit = f1->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  Double_t bkglowfit = f1->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1); //No error
  Double_t bkghighfit = f1->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1); //No error

  if(draw || finalscan){
    cout << "Background entries parametrised fit: " << fup->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1) << " " << fcent->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1) << " " << flow->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1) << endl;
    
    bkgcentfit = fcent->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
    bkglowfit = flow->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
    bkghighfit = fup->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  }
  if(bincountBkg[j]){
    bkgcent = bkgcentbc;
    bkglow = bkglowbc;
    bkghigh = bkghighbc;
  } else {
    bkgcent = bkgcentfit;
    bkglow = bkglowfit;
    bkghigh = bkghighfit;
  }
  
  if(draw){
    fcent->Draw("same");
    fup->Draw("same");
    flow->Draw("same");
    
    gPad->SetTickx();
    gPad->SetTicky();
    
    TPaveText *pinfos=new TPaveText(0.12,0.72,0.47,0.85,"NDC");
    pinfos->SetBorderSize(0);
    pinfos->SetFillStyle(0);
    pinfos->AddText(Form("Events = %s",nEv335String.Data()));
    pinfos->AddText(Form("B (3#sigma) = %.2f + %.2f - %.2f (BC)",bkgcentbc, bkghighbc - bkgcentbc, bkgcentbc - bkglowbc));
    pinfos->AddText(Form("B (3#sigma) = %.2f + %.2f - %.2f (fit)",bkgcentfit, bkghighfit - bkgcentfit, bkgcentfit - bkglowfit));
    pinfos->Draw();
    
    Double_t maxplotax = hbkg->GetMaximum();
    hbkg->GetYaxis()->SetRangeUser(0., 1.3 * maxplotax);
    
    Double_t maxplot = hbkg->GetBinContent(hbkg->FindBin(5.366));
    if(maxplot < 10) maxplot = 2;
    Double_t xgr[1] = {xmin + 0.5*(xmax - xmin)};
    Double_t ygr[1] = {0.5 * maxplot};
    Double_t xerrgr[1] = {0.5*(xmax - xmin)};
    Double_t yerrgr[1] = {0.5 * maxplot};
    TGraphErrors* grrange = new TGraphErrors(1, xgr, ygr, xerrgr, yerrgr);
    grrange->SetFillColor(kRed); grrange->SetLineWidth(2); grrange->SetLineColor(kGray+1); grrange->SetFillStyle(3445);
    grrange->Draw("same2");
    
    hbkg->Draw("same ep");
  }
}

void extract_fonll(TString filnam, int j, double &fonllmin, double &fonllcent, double &fonllmax, bool draw){

  if(filnam=="") return 0x0;
  FILE* infil=fopen(filnam.Data(),"r");
  Char_t line[200];

  for(Int_t il=0; il<18; il++){
    fgets(line,200,infil);
    if(strstr(line,"central")) break;
  }
  Float_t ptmin,ptmax,csc,csmin,csmax,dum;
  fscanf(infil,"%f",&ptmin);

  //--Draw--//
  Int_t iPt=0;
  Int_t iPt2=0;
  TGraphAsymmErrors* gf=new TGraphAsymmErrors(0);
  TGraphAsymmErrors* gfpt=new TGraphAsymmErrors(0);
  //--Draw--//
  
  while(!feof(infil)){
    fscanf(infil,"%f %f %f",&csc,&csmin,&csmax);
    for(Int_t i=0; i<12;i++) fscanf(infil,"%f",&dum);
    if(feof(infil)) break;
    fscanf(infil,"%f",&ptmax);
    Double_t ptmed=0.5*(ptmin+ptmax);
    Double_t dpt=(ptmax-ptmin);
    
    //Multiplying by F(b->B) and f(B/Bs) factor
    Double_t normFact=fLHCbBBs*fbtoB*1e-6/dpt; //from pb to ub/GeV/c.
    csc*=normFact;
    csmin*=normFact;
    csmax*=normFact;
    
    Double_t systFF=relSystFF*csc;
    Double_t systLHCb=relSystLHCb*csc;
    Double_t errup=csmax-csc;
    Double_t errdw=csc-csmin;
    Double_t errtotup=TMath::Sqrt(errup*errup+systFF*systFF+systLHCb*systLHCb);
    Double_t errtotdw=TMath::Sqrt(errdw*errdw+systFF*systFF+systLHCb*systLHCb);
    
    if(draw){
      gf->SetPoint(iPt,ptmed,csc);
      gf->SetPointError(iPt,0.5*dpt,0.5*dpt,errtotdw,errtotup);
      iPt++;
    }
    ptmin=ptmax;

    if(ptmed > ptbins[j] && ptmed < ptbins[j+1]){
      fonllmin = csc - errtotdw;
      fonllcent = csc;
      fonllmax = csc + errtotup;

      if(draw){
        gfpt->SetPoint(iPt2,ptmed,csc);
        gfpt->SetPointError(iPt2,0.5*dpt,0.5*dpt,errtotdw,errtotup);
        iPt2++;
      }
    }
  }

  if(draw){
    gf->SetFillColor(kAzure+1);
    gf->SetLineColor(kAzure+1);
    gf->SetLineWidth(2);
    gf->SetFillStyle(3002);
    gfpt->SetFillColor(2);
    gfpt->SetLineColor(2);
    gfpt->SetLineWidth(2);
    gfpt->SetFillStyle(3002);
    
    //printf("----FONLL-----\n");
    //gf->Print();
    TGraphAsymmErrors* gfcent=(TGraphAsymmErrors*)gf->Clone("gFONLLcent");
    for(Int_t jp=0; jp<gfcent->GetN(); jp++){
      gfcent->SetPointEYlow(jp,0.);
      gfcent->SetPointEYhigh(jp,0.);
      gfcent->SetLineWidth(3);
    }
    TGraphAsymmErrors* gfptcent=(TGraphAsymmErrors*)gfpt->Clone("gFONLLptcent");
    for(Int_t jp=0; jp<gfptcent->GetN(); jp++){
      gfptcent->SetPointEYlow(jp,0.);
      gfptcent->SetPointEYhigh(jp,0.);
      gfptcent->SetLineWidth(3);
    }
    
    gf->SetTitle("FONLL b #rightarrow B (5.02 TeV)");
    gf->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    gf->GetYaxis()->SetTitle("d^{2}#sigma/(d#it{p}_{T}d#it{y}) (#mub GeV^{-1} #it{c})");
    gf->Draw("ae5s");
    gPad->SetLogy();
    gPad->SetTickx();
    gPad->SetTicky();
    
    gf->GetYaxis()->SetRangeUser(0.001,10);
    gfcent->Draw("ezsame");
    
    gfpt->Draw("same e5s");
    gfptcent->Draw("ezsame");
    
    TPaveText *pinfos=new TPaveText(0.22,0.72,0.87,0.85,"NDC");
    pinfos->SetBorderSize(0);
    pinfos->SetFillStyle(0);
    pinfos->AddText(Form("Multiplied by F(b->B) = %.3f #pm %.3f",fbtoB, fbtoBUnc));
    pinfos->AddText(Form("Multiplied by F(B/B_{s})^{LHCb} = %.3f #pm %.3f",fLHCbBBs, fLHCbBBsUnc));
    pinfos->Draw();
  }
  
  fclose(infil);
}

void extract_TAMU(TString filn, int j, double &tamucent, bool draw){
  
  FILE* f=fopen(filn.Data(),"r");

  //--Draw--//
  TGraph* g=new TGraph(0);
  TGraphErrors* gpt=new TGraphErrors(0);
  //--Draw--//
  
  Float_t pt, raa;
  Int_t iPt=0;
  Double_t meanRAA = 0.;
  Int_t ncount = 0;
  while(!feof(f)){
    fscanf(f,"%f %f\n",&pt,&raa);
    if(draw) g->SetPoint(iPt++,pt,raa);

    if(pt > ptbins[j] && pt <= ptbins[j+1]){
      meanRAA += raa;
      ncount++;
    }
  }
  tamucent = meanRAA / ((double)ncount);
  
  if(draw){
    gpt->SetPoint(0, ptbins[j] + 0.5*(ptbins[j+1] - ptbins[j]), tamucent);
    gpt->SetPointError(0,0.5*(ptbins[j+1] - ptbins[j]),0);
    
    g->SetLineColor(kAzure+1);
    g->SetLineWidth(3);
    g->SetLineStyle(1);
    
    gpt->SetLineColor(kRed);
    gpt->SetLineWidth(2);
    gpt->SetLineStyle(1);
    
    g->SetTitle("TAMU B_{s};#it{p}_{T} (GeV/#it{c});#it{R}_{AA}");
    g->Draw("al");
    gPad->SetTickx();
    gPad->SetTicky();
    g->GetXaxis()->SetRangeUser(0.,25.);
    g->GetYaxis()->SetRangeUser(0.,2.5);
    gpt->Draw("same");
    
    TLine* l = new TLine(0,1,25,1);
    l->SetLineStyle(2);
    l->Draw();
    
    TPaveText *pinfos=new TPaveText(0.52,0.72,0.87,0.85,"NDC");
    pinfos->SetBorderSize(0);
    pinfos->SetFillStyle(0);
    pinfos->AddText("Pb-Pb, 5.02 TeV, 0-20%");
    pinfos->AddText("at #it{T}_{kin}");
    pinfos->Draw();
  }
  
  fclose(f);
}

void calculate_efficiency(TH1F* heff, int j, int i, double &effcent, double &erreffcent, bool draw){

  if(draw){
    heff->SetTitle("Acceptance-times-efficiency;#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon) #times 2#it{y}_{fid}");
    heff->SetName(Form("heff_set%d_pt%d", i, j+1));
    heff->SetMarkerStyle(20);
    heff->SetMarkerColor(kAzure+1);
    heff->SetLineColor(kAzure+1);
    heff->SetLineWidth(2);
    heff->SetStats(0);
    
    heff->Draw("ep");
    gPad->SetTickx();
    gPad->SetTicky();
    heff->GetYaxis()->SetRangeUser(0.001,1);
    gPad->SetLogy();
  }
  TH1F* heffcent = new TH1F(Form("heffcent_set%d_%d",i,j),"",1,ptbins[j],ptbins[j+1]);
  heffcent->SetBinContent(1, heff->GetBinContent(heff->FindBin(ptbins[j] + 0.5*(ptbins[j+1]-ptbins[j]))));
  heffcent->SetBinError(1, heff->GetBinError(heff->FindBin(ptbins[j] + 0.5*(ptbins[j+1]-ptbins[j]))));
  
  if(draw){
    heffcent->SetMarkerStyle(20);
    heffcent->SetMarkerColor(kRed);
    heffcent->SetLineColor(kRed);
    heffcent->SetLineWidth(2);
    //heffcent->Draw("same ep");
  }
  
  effcent = heffcent->GetBinContent(1);
  erreffcent = heffcent->GetBinError(1);
}

void InitializeProbCuts(){
  
  Float_t probcuts_temp[trials+1] = {0.95, 0.95001, 0.95002, 0.95003, 0.95004, 0.95005, 0.95006, 0.95007, 0.95008, 0.95009, 0.9501, 0.95011, 0.95012, 0.95013, 0.95014, 0.95015, 0.95016, 0.95017, 0.95018, 0.95019, 0.9502, 0.95021, 0.95022, 0.95023, 0.95024, 0.95025, 0.95026, 0.95027, 0.95028, 0.95029, 0.95030, 0.95031, 0.95032, 0.95033, 0.95034, 0.95035, 0.95036, 0.95037, 0.95038, 0.95039, 0.9504, 0.95041, 0.95042, 0.95043, 0.95044, 0.95045, 0.95046, 0.95047, 0.95048, 0.95049, 0.95050, 0.95051, 0.95052, 0.95053, 0.95054, 0.95055, 0.95056, 0.95057, 0.95058, 0.95059, 0.9506, 0.95061, 0.95062, 0.95063, 0.95064, 0.95065, 0.95066, 0.95067, 0.95068, 0.95069, 0.9507, 0.95071, 0.95072, 0.95073, 0.95074, 0.95075, 0.95076, 0.95077, 0.95078, 0.95079, 0.9508, 0.95081, 0.95082, 0.95083, 0.95084, 0.95085, 0.95086, 0.95087, 0.95088, 0.95089, 0.9509, 0.95091, 0.95092, 0.95093, 0.95094, 0.95095, 0.95096, 0.95097, 0.95098, 0.95099, 0.951, 0.95101, 0.95102, 0.95103, 0.95104, 0.95105, 0.95106, 0.95107, 0.95108, 0.95109, 0.9511, 0.95111, 0.95112, 0.95113, 0.95114, 0.95115, 0.95116, 0.95117, 0.95118, 0.95119, 0.9512, 0.95121, 0.95122, 0.95123, 0.95124, 0.95125, 0.95126, 0.95127, 0.95128, 0.95129, 0.95130, 0.95131, 0.95132, 0.95133, 0.95134, 0.95135, 0.95136, 0.95137, 0.95138, 0.95139, 0.9514, 0.95141, 0.95142, 0.95143, 0.95144, 0.95145, 0.95146, 0.95147, 0.95148, 0.95149, 0.9515, 0.95151, 0.95152, 0.95153, 0.95154, 0.95155, 0.95156, 0.95157, 0.95158, 0.95159, 0.9516, 0.95161, 0.95162, 0.95163, 0.95164, 0.95165, 0.95166, 0.95167, 0.95168, 0.95169, 0.9517, 0.95171, 0.95172, 0.95173, 0.95174, 0.95175, 0.95176, 0.95177, 0.95178, 0.95179, 0.9518, 0.95181, 0.95182, 0.95183, 0.95184, 0.95185, 0.95186, 0.95187, 0.95188, 0.95189, 0.9519, 0.95191, 0.95192, 0.95193, 0.95194, 0.95195, 0.95196, 0.95197, 0.95198, 0.95199, 0.952, 0.95201, 0.95202, 0.95203, 0.95204, 0.95205, 0.95206, 0.95207, 0.95208, 0.95209, 0.9521, 0.95211, 0.95212, 0.95213, 0.95214, 0.95215, 0.95216, 0.95217, 0.95218, 0.95219, 0.9522, 0.95221, 0.95222, 0.95223, 0.95224, 0.95225, 0.95226, 0.95227, 0.95228, 0.95229, 0.95230, 0.95231, 0.95232, 0.95233, 0.95234, 0.95235, 0.95236, 0.95237, 0.95238, 0.95239, 0.9524, 0.95241, 0.95242, 0.95243, 0.95244, 0.95245, 0.95246, 0.95247, 0.95248, 0.95249, 0.95250, 0.95251, 0.95252, 0.95253, 0.95254, 0.95255, 0.95256, 0.95257, 0.95258, 0.95259, 0.9526, 0.95261, 0.95262, 0.95263, 0.95264, 0.95265, 0.95266, 0.95267, 0.95268, 0.95269, 0.9527, 0.95271, 0.95272, 0.95273, 0.95274, 0.95275, 0.95276, 0.95277, 0.95278, 0.95279, 0.9528, 0.95281, 0.95282, 0.95283, 0.95284, 0.95285, 0.95286, 0.95287, 0.95288, 0.95289, 0.9529, 0.95291, 0.95292, 0.95293, 0.95294, 0.95295, 0.95296, 0.95297, 0.95298, 0.95299, 0.953, 0.95301, 0.95302, 0.95303, 0.95304, 0.95305, 0.95306, 0.95307, 0.95308, 0.95309, 0.9531, 0.95311, 0.95312, 0.95313, 0.95314, 0.95315, 0.95316, 0.95317, 0.95318, 0.95319, 0.9532, 0.95321, 0.95322, 0.95323, 0.95324, 0.95325, 0.95326, 0.95327, 0.95328, 0.95329, 0.95330, 0.95331, 0.95332, 0.95333, 0.95334, 0.95335, 0.95336, 0.95337, 0.95338, 0.95339, 0.9534, 0.95341, 0.95342, 0.95343, 0.95344, 0.95345, 0.95346, 0.95347, 0.95348, 0.95349, 0.9535, 0.95351, 0.95352, 0.95353, 0.95354, 0.95355, 0.95356, 0.95357, 0.95358, 0.95359, 0.9536, 0.95361, 0.95362, 0.95363, 0.95364, 0.95365, 0.95366, 0.95367, 0.95368, 0.95369, 0.9537, 0.95371, 0.95372, 0.95373, 0.95374, 0.95375, 0.95376, 0.95377, 0.95378, 0.95379, 0.9538, 0.95381, 0.95382, 0.95383, 0.95384, 0.95385, 0.95386, 0.95387, 0.95388, 0.95389, 0.9539, 0.95391, 0.95392, 0.95393, 0.95394, 0.95395, 0.95396, 0.95397, 0.95398, 0.95399, 0.954, 0.95401, 0.95402, 0.95403, 0.95404, 0.95405, 0.95406, 0.95407, 0.95408, 0.95409, 0.9541, 0.95411, 0.95412, 0.95413, 0.95414, 0.95415, 0.95416, 0.95417, 0.95418, 0.95419, 0.9542, 0.95421, 0.95422, 0.95423, 0.95424, 0.95425, 0.95426, 0.95427, 0.95428, 0.95429, 0.95430, 0.95431, 0.95432, 0.95433, 0.95434, 0.95435, 0.95436, 0.95437, 0.95438, 0.95439, 0.9544, 0.95441, 0.95442, 0.95443, 0.95444, 0.95445, 0.95446, 0.95447, 0.95448, 0.95449, 0.95450, 0.95451, 0.95452, 0.95453, 0.95454, 0.95455, 0.95456, 0.95457, 0.95458, 0.95459, 0.9546, 0.95461, 0.95462, 0.95463, 0.95464, 0.95465, 0.95466, 0.95467, 0.95468, 0.95469, 0.9547, 0.95471, 0.95472, 0.95473, 0.95474, 0.95475, 0.95476, 0.95477, 0.95478, 0.95479, 0.9548, 0.95481, 0.95482, 0.95483, 0.95484, 0.95485, 0.95486, 0.95487, 0.95488, 0.95489, 0.9549, 0.95491, 0.95492, 0.95493, 0.95494, 0.95495, 0.95496, 0.95497, 0.95498, 0.95499, 0.95500, 0.95501, 0.95502, 0.95503, 0.95504, 0.95505, 0.95506, 0.95507, 0.95508, 0.95509, 0.9551, 0.95511, 0.95512, 0.95513, 0.95514, 0.95515, 0.95516, 0.95517, 0.95518, 0.95519, 0.9552, 0.95521, 0.95522, 0.95523, 0.95524, 0.95525, 0.95526, 0.95527, 0.95528, 0.95529, 0.95530, 0.95531, 0.95532, 0.95533, 0.95534, 0.95535, 0.95536, 0.95537, 0.95538, 0.95539, 0.9554, 0.95541, 0.95542, 0.95543, 0.95544, 0.95545, 0.95546, 0.95547, 0.95548, 0.95549, 0.9555, 0.95551, 0.95552, 0.95553, 0.95554, 0.95555, 0.95556, 0.95557, 0.95558, 0.95559, 0.9556, 0.95561, 0.95562, 0.95563, 0.95564, 0.95565, 0.95566, 0.95567, 0.95568, 0.95569, 0.9557, 0.95571, 0.95572, 0.95573, 0.95574, 0.95575, 0.95576, 0.95577, 0.95578, 0.95579, 0.9558, 0.95581, 0.95582, 0.95583, 0.95584, 0.95585, 0.95586, 0.95587, 0.95588, 0.95589, 0.95590, 0.95591, 0.95592, 0.95593, 0.95594, 0.95595, 0.95596, 0.95597, 0.95598, 0.95599, 0.956, 0.95601, 0.95602, 0.95603, 0.95604, 0.95605, 0.95606, 0.95607, 0.95608, 0.95609, 0.9561, 0.95611, 0.95612, 0.95613, 0.95614, 0.95615, 0.95616, 0.95617, 0.95618, 0.95619, 0.9562, 0.95621, 0.95622, 0.95623, 0.95624, 0.95625, 0.95626, 0.95627, 0.95628, 0.95629, 0.95630, 0.95631, 0.95632, 0.95633, 0.95634, 0.95635, 0.95636, 0.95637, 0.95638, 0.95639, 0.9564, 0.95641, 0.95642, 0.95643, 0.95644, 0.95645, 0.95646, 0.95647, 0.95648, 0.95649, 0.95650, 0.95651, 0.95652, 0.95653, 0.95654, 0.95655, 0.95656, 0.95657, 0.95658, 0.95659, 0.9566, 0.95661, 0.95662, 0.95663, 0.95664, 0.95665, 0.95666, 0.95667, 0.95668, 0.95669, 0.9567, 0.95671, 0.95672, 0.95673, 0.95674, 0.95675, 0.95676, 0.95677, 0.95678, 0.95679, 0.9568, 0.95681, 0.95682, 0.95683, 0.95684, 0.95685, 0.95686, 0.95687, 0.95688, 0.95689, 0.9569, 0.95691, 0.95692, 0.95693, 0.95694, 0.95695, 0.95696, 0.95697, 0.95698, 0.95699, 0.957, 0.95701, 0.95702, 0.95703, 0.95704, 0.95705, 0.95706, 0.95707, 0.95708, 0.95709, 0.9571, 0.95711, 0.95712, 0.95713, 0.95714, 0.95715, 0.95716, 0.95717, 0.95718, 0.95719, 0.9572, 0.95721, 0.95722, 0.95723, 0.95724, 0.95725, 0.95726, 0.95727, 0.95728, 0.95729, 0.95730, 0.95731, 0.95732, 0.95733, 0.95734, 0.95735, 0.95736, 0.95737, 0.95738, 0.95739, 0.9574, 0.95741, 0.95742, 0.95743, 0.95744, 0.95745, 0.95746, 0.95747, 0.95748, 0.95749, 0.9575, 0.95751, 0.95752, 0.95753, 0.95754, 0.95755, 0.95756, 0.95757, 0.95758, 0.95759, 0.9576, 0.95761, 0.95762, 0.95763, 0.95764, 0.95765, 0.95766, 0.95767, 0.95768, 0.95769, 0.9577, 0.95771, 0.95772, 0.95773, 0.95774, 0.95775, 0.95776, 0.95777, 0.95778, 0.95779, 0.9578, 0.95781, 0.95782, 0.95783, 0.95784, 0.95785, 0.95786, 0.95787, 0.95788, 0.95789, 0.9579, 0.95791, 0.95792, 0.95793, 0.95794, 0.95795, 0.95796, 0.95797, 0.95798, 0.95799, 0.958, 0.95801, 0.95802, 0.95803, 0.95804, 0.95805, 0.95806, 0.95807, 0.95808, 0.95809, 0.9581, 0.95811, 0.95812, 0.95813, 0.95814, 0.95815, 0.95816, 0.95817, 0.95818, 0.95819, 0.9582, 0.95821, 0.95822, 0.95823, 0.95824, 0.95825, 0.95826, 0.95827, 0.95828, 0.95829, 0.95830, 0.95831, 0.95832, 0.95833, 0.95834, 0.95835, 0.95836, 0.95837, 0.95838, 0.95839, 0.9584, 0.95841, 0.95842, 0.95843, 0.95844, 0.95845, 0.95846, 0.95847, 0.95848, 0.95849, 0.95850, 0.95851, 0.95852, 0.95853, 0.95854, 0.95855, 0.95856, 0.95857, 0.95858, 0.95859, 0.9586, 0.95861, 0.95862, 0.95863, 0.95864, 0.95865, 0.95866, 0.95867, 0.95868, 0.95869, 0.9587, 0.95871, 0.95872, 0.95873, 0.95874, 0.95875, 0.95876, 0.95877, 0.95878, 0.95879, 0.9588, 0.95881, 0.95882, 0.95883, 0.95884, 0.95885, 0.95886, 0.95887, 0.95888, 0.95889, 0.9589, 0.95891, 0.95892, 0.95893, 0.95894, 0.95895, 0.95896, 0.95897, 0.95898, 0.95899, 0.959, 0.95901, 0.95902, 0.95903, 0.95904, 0.95905, 0.95906, 0.95907, 0.95908, 0.95909, 0.9591, 0.95911, 0.95912, 0.95913, 0.95914, 0.95915, 0.95916, 0.95917, 0.95918, 0.95919, 0.9592, 0.95921, 0.95922, 0.95923, 0.95924, 0.95925, 0.95926, 0.95927, 0.95928, 0.95929, 0.95930, 0.95931, 0.95932, 0.95933, 0.95934, 0.95935, 0.95936, 0.95937, 0.95938, 0.95939, 0.9594, 0.95941, 0.95942, 0.95943, 0.95944, 0.95945, 0.95946, 0.95947, 0.95948, 0.95949, 0.9595, 0.95951, 0.95952, 0.95953, 0.95954, 0.95955, 0.95956, 0.95957, 0.95958, 0.95959, 0.9596, 0.95961, 0.95962, 0.95963, 0.95964, 0.95965, 0.95966, 0.95967, 0.95968, 0.95969, 0.9597, 0.95971, 0.95972, 0.95973, 0.95974, 0.95975, 0.95976, 0.95977, 0.95978, 0.95979, 0.9598, 0.95981, 0.95982, 0.95983, 0.95984, 0.95985, 0.95986, 0.95987, 0.95988, 0.95989, 0.9599, 0.95991, 0.95992, 0.95993, 0.95994, 0.95995, 0.95996, 0.95997, 0.95998, 0.95999, 0.96, 0.96001, 0.96002, 0.96003, 0.96004, 0.96005, 0.96006, 0.96007, 0.96008, 0.96009, 0.9601, 0.96011, 0.96012, 0.96013, 0.96014, 0.96015, 0.96016, 0.96017, 0.96018, 0.96019, 0.9602, 0.96021, 0.96022, 0.96023, 0.96024, 0.96025, 0.96026, 0.96027, 0.96028, 0.96029, 0.96030, 0.96031, 0.96032, 0.96033, 0.96034, 0.96035, 0.96036, 0.96037, 0.96038, 0.96039, 0.9604, 0.96041, 0.96042, 0.96043, 0.96044, 0.96045, 0.96046, 0.96047, 0.96048, 0.96049, 0.96050, 0.96051, 0.96052, 0.96053, 0.96054, 0.96055, 0.96056, 0.96057, 0.96058, 0.96059, 0.9606, 0.96061, 0.96062, 0.96063, 0.96064, 0.96065, 0.96066, 0.96067, 0.96068, 0.96069, 0.9607, 0.96071, 0.96072, 0.96073, 0.96074, 0.96075, 0.96076, 0.96077, 0.96078, 0.96079, 0.9608, 0.96081, 0.96082, 0.96083, 0.96084, 0.96085, 0.96086, 0.96087, 0.96088, 0.96089, 0.9609, 0.96091, 0.96092, 0.96093, 0.96094, 0.96095, 0.96096, 0.96097, 0.96098, 0.96099, 0.961, 0.96101, 0.96102, 0.96103, 0.96104, 0.96105, 0.96106, 0.96107, 0.96108, 0.96109, 0.9611, 0.96111, 0.96112, 0.96113, 0.96114, 0.96115, 0.96116, 0.96117, 0.96118, 0.96119, 0.9612, 0.96121, 0.96122, 0.96123, 0.96124, 0.96125, 0.96126, 0.96127, 0.96128, 0.96129, 0.96130, 0.96131, 0.96132, 0.96133, 0.96134, 0.96135, 0.96136, 0.96137, 0.96138, 0.96139, 0.9614, 0.96141, 0.96142, 0.96143, 0.96144, 0.96145, 0.96146, 0.96147, 0.96148, 0.96149, 0.9615, 0.96151, 0.96152, 0.96153, 0.96154, 0.96155, 0.96156, 0.96157, 0.96158, 0.96159, 0.9616, 0.96161, 0.96162, 0.96163, 0.96164, 0.96165, 0.96166, 0.96167, 0.96168, 0.96169, 0.9617, 0.96171, 0.96172, 0.96173, 0.96174, 0.96175, 0.96176, 0.96177, 0.96178, 0.96179, 0.9618, 0.96181, 0.96182, 0.96183, 0.96184, 0.96185, 0.96186, 0.96187, 0.96188, 0.96189, 0.9619, 0.96191, 0.96192, 0.96193, 0.96194, 0.96195, 0.96196, 0.96197, 0.96198, 0.96199, 0.962, 0.96201, 0.96202, 0.96203, 0.96204, 0.96205, 0.96206, 0.96207, 0.96208, 0.96209, 0.9621, 0.96211, 0.96212, 0.96213, 0.96214, 0.96215, 0.96216, 0.96217, 0.96218, 0.96219, 0.9622, 0.96221, 0.96222, 0.96223, 0.96224, 0.96225, 0.96226, 0.96227, 0.96228, 0.96229, 0.96230, 0.96231, 0.96232, 0.96233, 0.96234, 0.96235, 0.96236, 0.96237, 0.96238, 0.96239, 0.9624, 0.96241, 0.96242, 0.96243, 0.96244, 0.96245, 0.96246, 0.96247, 0.96248, 0.96249, 0.96250, 0.96251, 0.96252, 0.96253, 0.96254, 0.96255, 0.96256, 0.96257, 0.96258, 0.96259, 0.9626, 0.96261, 0.96262, 0.96263, 0.96264, 0.96265, 0.96266, 0.96267, 0.96268, 0.96269, 0.9627, 0.96271, 0.96272, 0.96273, 0.96274, 0.96275, 0.96276, 0.96277, 0.96278, 0.96279, 0.9628, 0.96281, 0.96282, 0.96283, 0.96284, 0.96285, 0.96286, 0.96287, 0.96288, 0.96289, 0.9629, 0.96291, 0.96292, 0.96293, 0.96294, 0.96295, 0.96296, 0.96297, 0.96298, 0.96299, 0.963, 0.96301, 0.96302, 0.96303, 0.96304, 0.96305, 0.96306, 0.96307, 0.96308, 0.96309, 0.9631, 0.96311, 0.96312, 0.96313, 0.96314, 0.96315, 0.96316, 0.96317, 0.96318, 0.96319, 0.9632, 0.96321, 0.96322, 0.96323, 0.96324, 0.96325, 0.96326, 0.96327, 0.96328, 0.96329, 0.96330, 0.96331, 0.96332, 0.96333, 0.96334, 0.96335, 0.96336, 0.96337, 0.96338, 0.96339, 0.9634, 0.96341, 0.96342, 0.96343, 0.96344, 0.96345, 0.96346, 0.96347, 0.96348, 0.96349, 0.9635, 0.96351, 0.96352, 0.96353, 0.96354, 0.96355, 0.96356, 0.96357, 0.96358, 0.96359, 0.9636, 0.96361, 0.96362, 0.96363, 0.96364, 0.96365, 0.96366, 0.96367, 0.96368, 0.96369, 0.9637, 0.96371, 0.96372, 0.96373, 0.96374, 0.96375, 0.96376, 0.96377, 0.96378, 0.96379, 0.9638, 0.96381, 0.96382, 0.96383, 0.96384, 0.96385, 0.96386, 0.96387, 0.96388, 0.96389, 0.9639, 0.96391, 0.96392, 0.96393, 0.96394, 0.96395, 0.96396, 0.96397, 0.96398, 0.96399, 0.964, 0.96401, 0.96402, 0.96403, 0.96404, 0.96405, 0.96406, 0.96407, 0.96408, 0.96409, 0.9641, 0.96411, 0.96412, 0.96413, 0.96414, 0.96415, 0.96416, 0.96417, 0.96418, 0.96419, 0.9642, 0.96421, 0.96422, 0.96423, 0.96424, 0.96425, 0.96426, 0.96427, 0.96428, 0.96429, 0.96430, 0.96431, 0.96432, 0.96433, 0.96434, 0.96435, 0.96436, 0.96437, 0.96438, 0.96439, 0.9644, 0.96441, 0.96442, 0.96443, 0.96444, 0.96445, 0.96446, 0.96447, 0.96448, 0.96449, 0.96450, 0.96451, 0.96452, 0.96453, 0.96454, 0.96455, 0.96456, 0.96457, 0.96458, 0.96459, 0.9646, 0.96461, 0.96462, 0.96463, 0.96464, 0.96465, 0.96466, 0.96467, 0.96468, 0.96469, 0.9647, 0.96471, 0.96472, 0.96473, 0.96474, 0.96475, 0.96476, 0.96477, 0.96478, 0.96479, 0.9648, 0.96481, 0.96482, 0.96483, 0.96484, 0.96485, 0.96486, 0.96487, 0.96488, 0.96489, 0.9649, 0.96491, 0.96492, 0.96493, 0.96494, 0.96495, 0.96496, 0.96497, 0.96498, 0.96499, 0.965, 0.96501, 0.96502, 0.96503, 0.96504, 0.96505, 0.96506, 0.96507, 0.96508, 0.96509, 0.9651, 0.96511, 0.96512, 0.96513, 0.96514, 0.96515, 0.96516, 0.96517, 0.96518, 0.96519, 0.9652, 0.96521, 0.96522, 0.96523, 0.96524, 0.96525, 0.96526, 0.96527, 0.96528, 0.96529, 0.96530, 0.96531, 0.96532, 0.96533, 0.96534, 0.96535, 0.96536, 0.96537, 0.96538, 0.96539, 0.9654, 0.96541, 0.96542, 0.96543, 0.96544, 0.96545, 0.96546, 0.96547, 0.96548, 0.96549, 0.9655, 0.96551, 0.96552, 0.96553, 0.96554, 0.96555, 0.96556, 0.96557, 0.96558, 0.96559, 0.9656, 0.96561, 0.96562, 0.96563, 0.96564, 0.96565, 0.96566, 0.96567, 0.96568, 0.96569, 0.9657, 0.96571, 0.96572, 0.96573, 0.96574, 0.96575, 0.96576, 0.96577, 0.96578, 0.96579, 0.9658, 0.96581, 0.96582, 0.96583, 0.96584, 0.96585, 0.96586, 0.96587, 0.96588, 0.96589, 0.9659, 0.96591, 0.96592, 0.96593, 0.96594, 0.96595, 0.96596, 0.96597, 0.96598, 0.96599, 0.966, 0.96601, 0.96602, 0.96603, 0.96604, 0.96605, 0.96606, 0.96607, 0.96608, 0.96609, 0.9661, 0.96611, 0.96612, 0.96613, 0.96614, 0.96615, 0.96616, 0.96617, 0.96618, 0.96619, 0.9662, 0.96621, 0.96622, 0.96623, 0.96624, 0.96625, 0.96626, 0.96627, 0.96628, 0.96629, 0.96630, 0.96631, 0.96632, 0.96633, 0.96634, 0.96635, 0.96636, 0.96637, 0.96638, 0.96639, 0.9664, 0.96641, 0.96642, 0.96643, 0.96644, 0.96645, 0.96646, 0.96647, 0.96648, 0.96649, 0.96650, 0.96651, 0.96652, 0.96653, 0.96654, 0.96655, 0.96656, 0.96657, 0.96658, 0.96659, 0.9666, 0.96661, 0.96662, 0.96663, 0.96664, 0.96665, 0.96666, 0.96667, 0.96668, 0.96669, 0.9667, 0.96671, 0.96672, 0.96673, 0.96674, 0.96675, 0.96676, 0.96677, 0.96678, 0.96679, 0.9668, 0.96681, 0.96682, 0.96683, 0.96684, 0.96685, 0.96686, 0.96687, 0.96688, 0.96689, 0.9669, 0.96691, 0.96692, 0.96693, 0.96694, 0.96695, 0.96696, 0.96697, 0.96698, 0.96699, 0.967, 0.96701, 0.96702, 0.96703, 0.96704, 0.96705, 0.96706, 0.96707, 0.96708, 0.96709, 0.9671, 0.96711, 0.96712, 0.96713, 0.96714, 0.96715, 0.96716, 0.96717, 0.96718, 0.96719, 0.9672, 0.96721, 0.96722, 0.96723, 0.96724, 0.96725, 0.96726, 0.96727, 0.96728, 0.96729, 0.96730, 0.96731, 0.96732, 0.96733, 0.96734, 0.96735, 0.96736, 0.96737, 0.96738, 0.96739, 0.9674, 0.96741, 0.96742, 0.96743, 0.96744, 0.96745, 0.96746, 0.96747, 0.96748, 0.96749, 0.9675, 0.96751, 0.96752, 0.96753, 0.96754, 0.96755, 0.96756, 0.96757, 0.96758, 0.96759, 0.9676, 0.96761, 0.96762, 0.96763, 0.96764, 0.96765, 0.96766, 0.96767, 0.96768, 0.96769, 0.9677, 0.96771, 0.96772, 0.96773, 0.96774, 0.96775, 0.96776, 0.96777, 0.96778, 0.96779, 0.9678, 0.96781, 0.96782, 0.96783, 0.96784, 0.96785, 0.96786, 0.96787, 0.96788, 0.96789, 0.9679, 0.96791, 0.96792, 0.96793, 0.96794, 0.96795, 0.96796, 0.96797, 0.96798, 0.96799, 0.968, 0.96801, 0.96802, 0.96803, 0.96804, 0.96805, 0.96806, 0.96807, 0.96808, 0.96809, 0.9681, 0.96811, 0.96812, 0.96813, 0.96814, 0.96815, 0.96816, 0.96817, 0.96818, 0.96819, 0.9682, 0.96821, 0.96822, 0.96823, 0.96824, 0.96825, 0.96826, 0.96827, 0.96828, 0.96829, 0.96830, 0.96831, 0.96832, 0.96833, 0.96834, 0.96835, 0.96836, 0.96837, 0.96838, 0.96839, 0.9684, 0.96841, 0.96842, 0.96843, 0.96844, 0.96845, 0.96846, 0.96847, 0.96848, 0.96849, 0.96850, 0.96851, 0.96852, 0.96853, 0.96854, 0.96855, 0.96856, 0.96857, 0.96858, 0.96859, 0.9686, 0.96861, 0.96862, 0.96863, 0.96864, 0.96865, 0.96866, 0.96867, 0.96868, 0.96869, 0.9687, 0.96871, 0.96872, 0.96873, 0.96874, 0.96875, 0.96876, 0.96877, 0.96878, 0.96879, 0.9688, 0.96881, 0.96882, 0.96883, 0.96884, 0.96885, 0.96886, 0.96887, 0.96888, 0.96889, 0.9689, 0.96891, 0.96892, 0.96893, 0.96894, 0.96895, 0.96896, 0.96897, 0.96898, 0.96899, 0.969, 0.96901, 0.96902, 0.96903, 0.96904, 0.96905, 0.96906, 0.96907, 0.96908, 0.96909, 0.9691, 0.96911, 0.96912, 0.96913, 0.96914, 0.96915, 0.96916, 0.96917, 0.96918, 0.96919, 0.9692, 0.96921, 0.96922, 0.96923, 0.96924, 0.96925, 0.96926, 0.96927, 0.96928, 0.96929, 0.96930, 0.96931, 0.96932, 0.96933, 0.96934, 0.96935, 0.96936, 0.96937, 0.96938, 0.96939, 0.9694, 0.96941, 0.96942, 0.96943, 0.96944, 0.96945, 0.96946, 0.96947, 0.96948, 0.96949, 0.9695, 0.96951, 0.96952, 0.96953, 0.96954, 0.96955, 0.96956, 0.96957, 0.96958, 0.96959, 0.9696, 0.96961, 0.96962, 0.96963, 0.96964, 0.96965, 0.96966, 0.96967, 0.96968, 0.96969, 0.9697, 0.96971, 0.96972, 0.96973, 0.96974, 0.96975, 0.96976, 0.96977, 0.96978, 0.96979, 0.9698, 0.96981, 0.96982, 0.96983, 0.96984, 0.96985, 0.96986, 0.96987, 0.96988, 0.96989, 0.9699, 0.96991, 0.96992, 0.96993, 0.96994, 0.96995, 0.96996, 0.96997, 0.96998, 0.96999, 0.97, 0.97001, 0.97002, 0.97003, 0.97004, 0.97005, 0.97006, 0.97007, 0.97008, 0.97009, 0.9701, 0.97011, 0.97012, 0.97013, 0.97014, 0.97015, 0.97016, 0.97017, 0.97018, 0.97019, 0.9702, 0.97021, 0.97022, 0.97023, 0.97024, 0.97025, 0.97026, 0.97027, 0.97028, 0.97029, 0.97030, 0.97031, 0.97032, 0.97033, 0.97034, 0.97035, 0.97036, 0.97037, 0.97038, 0.97039, 0.9704, 0.97041, 0.97042, 0.97043, 0.97044, 0.97045, 0.97046, 0.97047, 0.97048, 0.97049, 0.97050, 0.97051, 0.97052, 0.97053, 0.97054, 0.97055, 0.97056, 0.97057, 0.97058, 0.97059, 0.9706, 0.97061, 0.97062, 0.97063, 0.97064, 0.97065, 0.97066, 0.97067, 0.97068, 0.97069, 0.9707, 0.97071, 0.97072, 0.97073, 0.97074, 0.97075, 0.97076, 0.97077, 0.97078, 0.97079, 0.9708, 0.97081, 0.97082, 0.97083, 0.97084, 0.97085, 0.97086, 0.97087, 0.97088, 0.97089, 0.9709, 0.97091, 0.97092, 0.97093, 0.97094, 0.97095, 0.97096, 0.97097, 0.97098, 0.97099, 0.971, 0.97101, 0.97102, 0.97103, 0.97104, 0.97105, 0.97106, 0.97107, 0.97108, 0.97109, 0.9711, 0.97111, 0.97112, 0.97113, 0.97114, 0.97115, 0.97116, 0.97117, 0.97118, 0.97119, 0.9712, 0.97121, 0.97122, 0.97123, 0.97124, 0.97125, 0.97126, 0.97127, 0.97128, 0.97129, 0.97130, 0.97131, 0.97132, 0.97133, 0.97134, 0.97135, 0.97136, 0.97137, 0.97138, 0.97139, 0.9714, 0.97141, 0.97142, 0.97143, 0.97144, 0.97145, 0.97146, 0.97147, 0.97148, 0.97149, 0.9715, 0.97151, 0.97152, 0.97153, 0.97154, 0.97155, 0.97156, 0.97157, 0.97158, 0.97159, 0.9716, 0.97161, 0.97162, 0.97163, 0.97164, 0.97165, 0.97166, 0.97167, 0.97168, 0.97169, 0.9717, 0.97171, 0.97172, 0.97173, 0.97174, 0.97175, 0.97176, 0.97177, 0.97178, 0.97179, 0.9718, 0.97181, 0.97182, 0.97183, 0.97184, 0.97185, 0.97186, 0.97187, 0.97188, 0.97189, 0.9719, 0.97191, 0.97192, 0.97193, 0.97194, 0.97195, 0.97196, 0.97197, 0.97198, 0.97199, 0.972, 0.97201, 0.97202, 0.97203, 0.97204, 0.97205, 0.97206, 0.97207, 0.97208, 0.97209, 0.9721, 0.97211, 0.97212, 0.97213, 0.97214, 0.97215, 0.97216, 0.97217, 0.97218, 0.97219, 0.9722, 0.97221, 0.97222, 0.97223, 0.97224, 0.97225, 0.97226, 0.97227, 0.97228, 0.97229, 0.97230, 0.97231, 0.97232, 0.97233, 0.97234, 0.97235, 0.97236, 0.97237, 0.97238, 0.97239, 0.9724, 0.97241, 0.97242, 0.97243, 0.97244, 0.97245, 0.97246, 0.97247, 0.97248, 0.97249, 0.97250, 0.97251, 0.97252, 0.97253, 0.97254, 0.97255, 0.97256, 0.97257, 0.97258, 0.97259, 0.9726, 0.97261, 0.97262, 0.97263, 0.97264, 0.97265, 0.97266, 0.97267, 0.97268, 0.97269, 0.9727, 0.97271, 0.97272, 0.97273, 0.97274, 0.97275, 0.97276, 0.97277, 0.97278, 0.97279, 0.9728, 0.97281, 0.97282, 0.97283, 0.97284, 0.97285, 0.97286, 0.97287, 0.97288, 0.97289, 0.9729, 0.97291, 0.97292, 0.97293, 0.97294, 0.97295, 0.97296, 0.97297, 0.97298, 0.97299, 0.973, 0.97301, 0.97302, 0.97303, 0.97304, 0.97305, 0.97306, 0.97307, 0.97308, 0.97309, 0.9731, 0.97311, 0.97312, 0.97313, 0.97314, 0.97315, 0.97316, 0.97317, 0.97318, 0.97319, 0.9732, 0.97321, 0.97322, 0.97323, 0.97324, 0.97325, 0.97326, 0.97327, 0.97328, 0.97329, 0.97330, 0.97331, 0.97332, 0.97333, 0.97334, 0.97335, 0.97336, 0.97337, 0.97338, 0.97339, 0.9734, 0.97341, 0.97342, 0.97343, 0.97344, 0.97345, 0.97346, 0.97347, 0.97348, 0.97349, 0.9735, 0.97351, 0.97352, 0.97353, 0.97354, 0.97355, 0.97356, 0.97357, 0.97358, 0.97359, 0.9736, 0.97361, 0.97362, 0.97363, 0.97364, 0.97365, 0.97366, 0.97367, 0.97368, 0.97369, 0.9737, 0.97371, 0.97372, 0.97373, 0.97374, 0.97375, 0.97376, 0.97377, 0.97378, 0.97379, 0.9738, 0.97381, 0.97382, 0.97383, 0.97384, 0.97385, 0.97386, 0.97387, 0.97388, 0.97389, 0.9739, 0.97391, 0.97392, 0.97393, 0.97394, 0.97395, 0.97396, 0.97397, 0.97398, 0.97399, 0.974, 0.97401, 0.97402, 0.97403, 0.97404, 0.97405, 0.97406, 0.97407, 0.97408, 0.97409, 0.9741, 0.97411, 0.97412, 0.97413, 0.97414, 0.97415, 0.97416, 0.97417, 0.97418, 0.97419, 0.9742, 0.97421, 0.97422, 0.97423, 0.97424, 0.97425, 0.97426, 0.97427, 0.97428, 0.97429, 0.97430, 0.97431, 0.97432, 0.97433, 0.97434, 0.97435, 0.97436, 0.97437, 0.97438, 0.97439, 0.9744, 0.97441, 0.97442, 0.97443, 0.97444, 0.97445, 0.97446, 0.97447, 0.97448, 0.97449, 0.97450, 0.97451, 0.97452, 0.97453, 0.97454, 0.97455, 0.97456, 0.97457, 0.97458, 0.97459, 0.9746, 0.97461, 0.97462, 0.97463, 0.97464, 0.97465, 0.97466, 0.97467, 0.97468, 0.97469, 0.9747, 0.97471, 0.97472, 0.97473, 0.97474, 0.97475, 0.97476, 0.97477, 0.97478, 0.97479, 0.9748, 0.97481, 0.97482, 0.97483, 0.97484, 0.97485, 0.97486, 0.97487, 0.97488, 0.97489, 0.9749, 0.97491, 0.97492, 0.97493, 0.97494, 0.97495, 0.97496, 0.97497, 0.97498, 0.97499, 0.975, 0.97501, 0.97502, 0.97503, 0.97504, 0.97505, 0.97506, 0.97507, 0.97508, 0.97509, 0.9751, 0.97511, 0.97512, 0.97513, 0.97514, 0.97515, 0.97516, 0.97517, 0.97518, 0.97519, 0.9752, 0.97521, 0.97522, 0.97523, 0.97524, 0.97525, 0.97526, 0.97527, 0.97528, 0.97529, 0.9753, 0.97531, 0.97532, 0.97533, 0.97534, 0.97535, 0.97536, 0.97537, 0.97538, 0.97539, 0.9754, 0.97541, 0.97542, 0.97543, 0.97544, 0.97545, 0.97546, 0.97547, 0.97548, 0.97549, 0.9755, 0.97551, 0.97552, 0.97553, 0.97554, 0.97555, 0.97556, 0.97557, 0.97558, 0.97559, 0.9756, 0.97561, 0.97562, 0.97563, 0.97564, 0.97565, 0.97566, 0.97567, 0.97568, 0.97569, 0.9757, 0.97571, 0.97572, 0.97573, 0.97574, 0.97575, 0.97576, 0.97577, 0.97578, 0.97579, 0.9758, 0.97581, 0.97582, 0.97583, 0.97584, 0.97585, 0.97586, 0.97587, 0.97588, 0.97589, 0.9759, 0.97591, 0.97592, 0.97593, 0.97594, 0.97595, 0.97596, 0.97597, 0.97598, 0.97599, 0.976, 0.97601, 0.97602, 0.97603, 0.97604, 0.97605, 0.97606, 0.97607, 0.97608, 0.97609, 0.9761, 0.97611, 0.97612, 0.97613, 0.97614, 0.97615, 0.97616, 0.97617, 0.97618, 0.97619, 0.9762, 0.97621, 0.97622, 0.97623, 0.97624, 0.97625, 0.97626, 0.97627, 0.97628, 0.97629, 0.9763, 0.97631, 0.97632, 0.97633, 0.97634, 0.97635, 0.97636, 0.97637, 0.97638, 0.97639, 0.9764, 0.97641, 0.97642, 0.97643, 0.97644, 0.97645, 0.97646, 0.97647, 0.97648, 0.97649, 0.97650, 0.97651, 0.97652, 0.97653, 0.97654, 0.97655, 0.97656, 0.97657, 0.97658, 0.97659, 0.9766, 0.97661, 0.97662, 0.97663, 0.97664, 0.97665, 0.97666, 0.97667, 0.97668, 0.97669, 0.9767, 0.97671, 0.97672, 0.97673, 0.97674, 0.97675, 0.97676, 0.97677, 0.97678, 0.97679, 0.9768, 0.97681, 0.97682, 0.97683, 0.97684, 0.97685, 0.97686, 0.97687, 0.97688, 0.97689, 0.9769, 0.97691, 0.97692, 0.97693, 0.97694, 0.97695, 0.97696, 0.97697, 0.97698, 0.97699, 0.977, 0.97701, 0.97702, 0.97703, 0.97704, 0.97705, 0.97706, 0.97707, 0.97708, 0.97709, 0.9771, 0.97711, 0.97712, 0.97713, 0.97714, 0.97715, 0.97716, 0.97717, 0.97718, 0.97719, 0.9772, 0.97721, 0.97722, 0.97723, 0.97724, 0.97725, 0.97726, 0.97727, 0.97728, 0.97729, 0.9773, 0.97731, 0.97732, 0.97733, 0.97734, 0.97735, 0.97736, 0.97737, 0.97738, 0.97739, 0.9774, 0.97741, 0.97742, 0.97743, 0.97744, 0.97745, 0.97746, 0.97747, 0.97748, 0.97749, 0.9775, 0.97751, 0.97752, 0.97753, 0.97754, 0.97755, 0.97756, 0.97757, 0.97758, 0.97759, 0.9776, 0.97761, 0.97762, 0.97763, 0.97764, 0.97765, 0.97766, 0.97767, 0.97768, 0.97769, 0.9777, 0.97771, 0.97772, 0.97773, 0.97774, 0.97775, 0.97776, 0.97777, 0.97778, 0.97779, 0.9778, 0.97781, 0.97782, 0.97783, 0.97784, 0.97785, 0.97786, 0.97787, 0.97788, 0.97789, 0.9779, 0.97791, 0.97792, 0.97793, 0.97794, 0.97795, 0.97796, 0.97797, 0.97798, 0.97799, 0.978, 0.97801, 0.97802, 0.97803, 0.97804, 0.97805, 0.97806, 0.97807, 0.97808, 0.97809, 0.9781, 0.97811, 0.97812, 0.97813, 0.97814, 0.97815, 0.97816, 0.97817, 0.97818, 0.97819, 0.9782, 0.97821, 0.97822, 0.97823, 0.97824, 0.97825, 0.97826, 0.97827, 0.97828, 0.97829, 0.9783, 0.97831, 0.97832, 0.97833, 0.97834, 0.97835, 0.97836, 0.97837, 0.97838, 0.97839, 0.9784, 0.97841, 0.97842, 0.97843, 0.97844, 0.97845, 0.97846, 0.97847, 0.97848, 0.97849, 0.97850, 0.97851, 0.97852, 0.97853, 0.97854, 0.97855, 0.97856, 0.97857, 0.97858, 0.97859, 0.9786, 0.97861, 0.97862, 0.97863, 0.97864, 0.97865, 0.97866, 0.97867, 0.97868, 0.97869, 0.9787, 0.97871, 0.97872, 0.97873, 0.97874, 0.97875, 0.97876, 0.97877, 0.97878, 0.97879, 0.9788, 0.97881, 0.97882, 0.97883, 0.97884, 0.97885, 0.97886, 0.97887, 0.97888, 0.97889, 0.9789, 0.97891, 0.97892, 0.97893, 0.97894, 0.97895, 0.97896, 0.97897, 0.97898, 0.97899, 0.979, 0.97901, 0.97902, 0.97903, 0.97904, 0.97905, 0.97906, 0.97907, 0.97908, 0.97909, 0.9791, 0.97911, 0.97912, 0.97913, 0.97914, 0.97915, 0.97916, 0.97917, 0.97918, 0.97919, 0.9792, 0.97921, 0.97922, 0.97923, 0.97924, 0.97925, 0.97926, 0.97927, 0.97928, 0.97929, 0.9793, 0.97931, 0.97932, 0.97933, 0.97934, 0.97935, 0.97936, 0.97937, 0.97938, 0.97939, 0.9794, 0.97941, 0.97942, 0.97943, 0.97944, 0.97945, 0.97946, 0.97947, 0.97948, 0.97949, 0.9795, 0.97951, 0.97952, 0.97953, 0.97954, 0.97955, 0.97956, 0.97957, 0.97958, 0.97959, 0.9796, 0.97961, 0.97962, 0.97963, 0.97964, 0.97965, 0.97966, 0.97967, 0.97968, 0.97969, 0.9797, 0.97971, 0.97972, 0.97973, 0.97974, 0.97975, 0.97976, 0.97977, 0.97978, 0.97979, 0.9798, 0.97981, 0.97982, 0.97983, 0.97984, 0.97985, 0.97986, 0.97987, 0.97988, 0.97989, 0.9799, 0.97991, 0.97992, 0.97993, 0.97994, 0.97995, 0.97996, 0.97997, 0.97998, 0.97999, 0.98, 0.98001, 0.98002, 0.98003, 0.98004, 0.98005, 0.98006, 0.98007, 0.98008, 0.98009, 0.9801, 0.98011, 0.98012, 0.98013, 0.98014, 0.98015, 0.98016, 0.98017, 0.98018, 0.98019, 0.9802, 0.98021, 0.98022, 0.98023, 0.98024, 0.98025, 0.98026, 0.98027, 0.98028, 0.98029, 0.9803, 0.98031, 0.98032, 0.98033, 0.98034, 0.98035, 0.98036, 0.98037, 0.98038, 0.98039, 0.9804, 0.98041, 0.98042, 0.98043, 0.98044, 0.98045, 0.98046, 0.98047, 0.98048, 0.98049, 0.98050, 0.98051, 0.98052, 0.98053, 0.98054, 0.98055, 0.98056, 0.98057, 0.98058, 0.98059, 0.9806, 0.98061, 0.98062, 0.98063, 0.98064, 0.98065, 0.98066, 0.98067, 0.98068, 0.98069, 0.9807, 0.98071, 0.98072, 0.98073, 0.98074, 0.98075, 0.98076, 0.98077, 0.98078, 0.98079, 0.9808, 0.98081, 0.98082, 0.98083, 0.98084, 0.98085, 0.98086, 0.98087, 0.98088, 0.98089, 0.9809, 0.98091, 0.98092, 0.98093, 0.98094, 0.98095, 0.98096, 0.98097, 0.98098, 0.98099, 0.981, 0.98101, 0.98102, 0.98103, 0.98104, 0.98105, 0.98106, 0.98107, 0.98108, 0.98109, 0.9811, 0.98111, 0.98112, 0.98113, 0.98114, 0.98115, 0.98116, 0.98117, 0.98118, 0.98119, 0.9812, 0.98121, 0.98122, 0.98123, 0.98124, 0.98125, 0.98126, 0.98127, 0.98128, 0.98129, 0.9813, 0.98131, 0.98132, 0.98133, 0.98134, 0.98135, 0.98136, 0.98137, 0.98138, 0.98139, 0.9814, 0.98141, 0.98142, 0.98143, 0.98144, 0.98145, 0.98146, 0.98147, 0.98148, 0.98149, 0.9815, 0.98151, 0.98152, 0.98153, 0.98154, 0.98155, 0.98156, 0.98157, 0.98158, 0.98159, 0.9816, 0.98161, 0.98162, 0.98163, 0.98164, 0.98165, 0.98166, 0.98167, 0.98168, 0.98169, 0.9817, 0.98171, 0.98172, 0.98173, 0.98174, 0.98175, 0.98176, 0.98177, 0.98178, 0.98179, 0.9818, 0.98181, 0.98182, 0.98183, 0.98184, 0.98185, 0.98186, 0.98187, 0.98188, 0.98189, 0.9819, 0.98191, 0.98192, 0.98193, 0.98194, 0.98195, 0.98196, 0.98197, 0.98198, 0.98199, 0.982, 0.98201, 0.98202, 0.98203, 0.98204, 0.98205, 0.98206, 0.98207, 0.98208, 0.98209, 0.9821, 0.98211, 0.98212, 0.98213, 0.98214, 0.98215, 0.98216, 0.98217, 0.98218, 0.98219, 0.9822, 0.98221, 0.98222, 0.98223, 0.98224, 0.98225, 0.98226, 0.98227, 0.98228, 0.98229, 0.9823, 0.98231, 0.98232, 0.98233, 0.98234, 0.98235, 0.98236, 0.98237, 0.98238, 0.98239, 0.9824, 0.98241, 0.98242, 0.98243, 0.98244, 0.98245, 0.98246, 0.98247, 0.98248, 0.98249, 0.98250, 0.98251, 0.98252, 0.98253, 0.98254, 0.98255, 0.98256, 0.98257, 0.98258, 0.98259, 0.9826, 0.98261, 0.98262, 0.98263, 0.98264, 0.98265, 0.98266, 0.98267, 0.98268, 0.98269, 0.9827, 0.98271, 0.98272, 0.98273, 0.98274, 0.98275, 0.98276, 0.98277, 0.98278, 0.98279, 0.9828, 0.98281, 0.98282, 0.98283, 0.98284, 0.98285, 0.98286, 0.98287, 0.98288, 0.98289, 0.9829, 0.98291, 0.98292, 0.98293, 0.98294, 0.98295, 0.98296, 0.98297, 0.98298, 0.98299, 0.983, 0.98301, 0.98302, 0.98303, 0.98304, 0.98305, 0.98306, 0.98307, 0.98308, 0.98309, 0.9831, 0.98311, 0.98312, 0.98313, 0.98314, 0.98315, 0.98316, 0.98317, 0.98318, 0.98319, 0.9832, 0.98321, 0.98322, 0.98323, 0.98324, 0.98325, 0.98326, 0.98327, 0.98328, 0.98329, 0.9833, 0.98331, 0.98332, 0.98333, 0.98334, 0.98335, 0.98336, 0.98337, 0.98338, 0.98339, 0.9834, 0.98341, 0.98342, 0.98343, 0.98344, 0.98345, 0.98346, 0.98347, 0.98348, 0.98349, 0.9835, 0.98351, 0.98352, 0.98353, 0.98354, 0.98355, 0.98356, 0.98357, 0.98358, 0.98359, 0.9836, 0.98361, 0.98362, 0.98363, 0.98364, 0.98365, 0.98366, 0.98367, 0.98368, 0.98369, 0.9837, 0.98371, 0.98372, 0.98373, 0.98374, 0.98375, 0.98376, 0.98377, 0.98378, 0.98379, 0.9838, 0.98381, 0.98382, 0.98383, 0.98384, 0.98385, 0.98386, 0.98387, 0.98388, 0.98389, 0.9839, 0.98391, 0.98392, 0.98393, 0.98394, 0.98395, 0.98396, 0.98397, 0.98398, 0.98399, 0.984, 0.98401, 0.98402, 0.98403, 0.98404, 0.98405, 0.98406, 0.98407, 0.98408, 0.98409, 0.9841, 0.98411, 0.98412, 0.98413, 0.98414, 0.98415, 0.98416, 0.98417, 0.98418, 0.98419, 0.9842, 0.98421, 0.98422, 0.98423, 0.98424, 0.98425, 0.98426, 0.98427, 0.98428, 0.98429, 0.9843, 0.98431, 0.98432, 0.98433, 0.98434, 0.98435, 0.98436, 0.98437, 0.98438, 0.98439, 0.9844, 0.98441, 0.98442, 0.98443, 0.98444, 0.98445, 0.98446, 0.98447, 0.98448, 0.98449, 0.98450, 0.98451, 0.98452, 0.98453, 0.98454, 0.98455, 0.98456, 0.98457, 0.98458, 0.98459, 0.9846, 0.98461, 0.98462, 0.98463, 0.98464, 0.98465, 0.98466, 0.98467, 0.98468, 0.98469, 0.9847, 0.98471, 0.98472, 0.98473, 0.98474, 0.98475, 0.98476, 0.98477, 0.98478, 0.98479, 0.9848, 0.98481, 0.98482, 0.98483, 0.98484, 0.98485, 0.98486, 0.98487, 0.98488, 0.98489, 0.9849, 0.98491, 0.98492, 0.98493, 0.98494, 0.98495, 0.98496, 0.98497, 0.98498, 0.98499, 0.985, 0.98501, 0.98502, 0.98503, 0.98504, 0.98505, 0.98506, 0.98507, 0.98508, 0.98509, 0.9851, 0.98511, 0.98512, 0.98513, 0.98514, 0.98515, 0.98516, 0.98517, 0.98518, 0.98519, 0.9852, 0.98521, 0.98522, 0.98523, 0.98524, 0.98525, 0.98526, 0.98527, 0.98528, 0.98529, 0.9853, 0.98531, 0.98532, 0.98533, 0.98534, 0.98535, 0.98536, 0.98537, 0.98538, 0.98539, 0.9854, 0.98541, 0.98542, 0.98543, 0.98544, 0.98545, 0.98546, 0.98547, 0.98548, 0.98549, 0.9855, 0.98551, 0.98552, 0.98553, 0.98554, 0.98555, 0.98556, 0.98557, 0.98558, 0.98559, 0.9856, 0.98561, 0.98562, 0.98563, 0.98564, 0.98565, 0.98566, 0.98567, 0.98568, 0.98569, 0.9857, 0.98571, 0.98572, 0.98573, 0.98574, 0.98575, 0.98576, 0.98577, 0.98578, 0.98579, 0.9858, 0.98581, 0.98582, 0.98583, 0.98584, 0.98585, 0.98586, 0.98587, 0.98588, 0.98589, 0.9859, 0.98591, 0.98592, 0.98593, 0.98594, 0.98595, 0.98596, 0.98597, 0.98598, 0.98599, 0.986, 0.98601, 0.98602, 0.98603, 0.98604, 0.98605, 0.98606, 0.98607, 0.98608, 0.98609, 0.9861, 0.98611, 0.98612, 0.98613, 0.98614, 0.98615, 0.98616, 0.98617, 0.98618, 0.98619, 0.9862, 0.98621, 0.98622, 0.98623, 0.98624, 0.98625, 0.98626, 0.98627, 0.98628, 0.98629, 0.9863, 0.98631, 0.98632, 0.98633, 0.98634, 0.98635, 0.98636, 0.98637, 0.98638, 0.98639, 0.9864, 0.98641, 0.98642, 0.98643, 0.98644, 0.98645, 0.98646, 0.98647, 0.98648, 0.98649, 0.98650, 0.98651, 0.98652, 0.98653, 0.98654, 0.98655, 0.98656, 0.98657, 0.98658, 0.98659, 0.9866, 0.98661, 0.98662, 0.98663, 0.98664, 0.98665, 0.98666, 0.98667, 0.98668, 0.98669, 0.9867, 0.98671, 0.98672, 0.98673, 0.98674, 0.98675, 0.98676, 0.98677, 0.98678, 0.98679, 0.9868, 0.98681, 0.98682, 0.98683, 0.98684, 0.98685, 0.98686, 0.98687, 0.98688, 0.98689, 0.9869, 0.98691, 0.98692, 0.98693, 0.98694, 0.98695, 0.98696, 0.98697, 0.98698, 0.98699, 0.987, 0.98701, 0.98702, 0.98703, 0.98704, 0.98705, 0.98706, 0.98707, 0.98708, 0.98709, 0.9871, 0.98711, 0.98712, 0.98713, 0.98714, 0.98715, 0.98716, 0.98717, 0.98718, 0.98719, 0.9872, 0.98721, 0.98722, 0.98723, 0.98724, 0.98725, 0.98726, 0.98727, 0.98728, 0.98729, 0.9873, 0.98731, 0.98732, 0.98733, 0.98734, 0.98735, 0.98736, 0.98737, 0.98738, 0.98739, 0.9874, 0.98741, 0.98742, 0.98743, 0.98744, 0.98745, 0.98746, 0.98747, 0.98748, 0.98749, 0.9875, 0.98751, 0.98752, 0.98753, 0.98754, 0.98755, 0.98756, 0.98757, 0.98758, 0.98759, 0.9876, 0.98761, 0.98762, 0.98763, 0.98764, 0.98765, 0.98766, 0.98767, 0.98768, 0.98769, 0.9877, 0.98771, 0.98772, 0.98773, 0.98774, 0.98775, 0.98776, 0.98777, 0.98778, 0.98779, 0.9878, 0.98781, 0.98782, 0.98783, 0.98784, 0.98785, 0.98786, 0.98787, 0.98788, 0.98789, 0.9879, 0.98791, 0.98792, 0.98793, 0.98794, 0.98795, 0.98796, 0.98797, 0.98798, 0.98799, 0.988, 0.98801, 0.98802, 0.98803, 0.98804, 0.98805, 0.98806, 0.98807, 0.98808, 0.98809, 0.9881, 0.98811, 0.98812, 0.98813, 0.98814, 0.98815, 0.98816, 0.98817, 0.98818, 0.98819, 0.9882, 0.98821, 0.98822, 0.98823, 0.98824, 0.98825, 0.98826, 0.98827, 0.98828, 0.98829, 0.9883, 0.98831, 0.98832, 0.98833, 0.98834, 0.98835, 0.98836, 0.98837, 0.98838, 0.98839, 0.9884, 0.98841, 0.98842, 0.98843, 0.98844, 0.98845, 0.98846, 0.98847, 0.98848, 0.98849, 0.98850, 0.98851, 0.98852, 0.98853, 0.98854, 0.98855, 0.98856, 0.98857, 0.98858, 0.98859, 0.9886, 0.98861, 0.98862, 0.98863, 0.98864, 0.98865, 0.98866, 0.98867, 0.98868, 0.98869, 0.9887, 0.98871, 0.98872, 0.98873, 0.98874, 0.98875, 0.98876, 0.98877, 0.98878, 0.98879, 0.9888, 0.98881, 0.98882, 0.98883, 0.98884, 0.98885, 0.98886, 0.98887, 0.98888, 0.98889, 0.9889, 0.98891, 0.98892, 0.98893, 0.98894, 0.98895, 0.98896, 0.98897, 0.98898, 0.98899, 0.989, 0.98901, 0.98902, 0.98903, 0.98904, 0.98905, 0.98906, 0.98907, 0.98908, 0.98909, 0.9891, 0.98911, 0.98912, 0.98913, 0.98914, 0.98915, 0.98916, 0.98917, 0.98918, 0.98919, 0.9892, 0.98921, 0.98922, 0.98923, 0.98924, 0.98925, 0.98926, 0.98927, 0.98928, 0.98929, 0.9893, 0.98931, 0.98932, 0.98933, 0.98934, 0.98935, 0.98936, 0.98937, 0.98938, 0.98939, 0.9894, 0.98941, 0.98942, 0.98943, 0.98944, 0.98945, 0.98946, 0.98947, 0.98948, 0.98949, 0.9895, 0.98951, 0.98952, 0.98953, 0.98954, 0.98955, 0.98956, 0.98957, 0.98958, 0.98959, 0.9896, 0.98961, 0.98962, 0.98963, 0.98964, 0.98965, 0.98966, 0.98967, 0.98968, 0.98969, 0.9897, 0.98971, 0.98972, 0.98973, 0.98974, 0.98975, 0.98976, 0.98977, 0.98978, 0.98979, 0.9898, 0.98981, 0.98982, 0.98983, 0.98984, 0.98985, 0.98986, 0.98987, 0.98988, 0.98989, 0.9899, 0.98991, 0.98992, 0.98993, 0.98994, 0.98995, 0.98996, 0.98997, 0.98998, 0.98999, 0.99, 0.99001, 0.99002, 0.99003, 0.99004, 0.99005, 0.99006, 0.99007, 0.99008, 0.99009, 0.9901, 0.99011, 0.99012, 0.99013, 0.99014, 0.99015, 0.99016, 0.99017, 0.99018, 0.99019, 0.9902, 0.99021, 0.99022, 0.99023, 0.99024, 0.99025, 0.99026, 0.99027, 0.99028, 0.99029, 0.9903, 0.99031, 0.99032, 0.99033, 0.99034, 0.99035, 0.99036, 0.99037, 0.99038, 0.99039, 0.9904, 0.99041, 0.99042, 0.99043, 0.99044, 0.99045, 0.99046, 0.99047, 0.99048, 0.99049, 0.99050, 0.99051, 0.99052, 0.99053, 0.99054, 0.99055, 0.99056, 0.99057, 0.99058, 0.99059, 0.9906, 0.99061, 0.99062, 0.99063, 0.99064, 0.99065, 0.99066, 0.99067, 0.99068, 0.99069, 0.9907, 0.99071, 0.99072, 0.99073, 0.99074, 0.99075, 0.99076, 0.99077, 0.99078, 0.99079, 0.9908, 0.99081, 0.99082, 0.99083, 0.99084, 0.99085, 0.99086, 0.99087, 0.99088, 0.99089, 0.9909, 0.99091, 0.99092, 0.99093, 0.99094, 0.99095, 0.99096, 0.99097, 0.99098, 0.99099, 0.991, 0.99101, 0.99102, 0.99103, 0.99104, 0.99105, 0.99106, 0.99107, 0.99108, 0.99109, 0.9911, 0.99111, 0.99112, 0.99113, 0.99114, 0.99115, 0.99116, 0.99117, 0.99118, 0.99119, 0.9912, 0.99121, 0.99122, 0.99123, 0.99124, 0.99125, 0.99126, 0.99127, 0.99128, 0.99129, 0.9913, 0.99131, 0.99132, 0.99133, 0.99134, 0.99135, 0.99136, 0.99137, 0.99138, 0.99139, 0.9914, 0.99141, 0.99142, 0.99143, 0.99144, 0.99145, 0.99146, 0.99147, 0.99148, 0.99149, 0.9915, 0.99151, 0.99152, 0.99153, 0.99154, 0.99155, 0.99156, 0.99157, 0.99158, 0.99159, 0.9916, 0.99161, 0.99162, 0.99163, 0.99164, 0.99165, 0.99166, 0.99167, 0.99168, 0.99169, 0.9917, 0.99171, 0.99172, 0.99173, 0.99174, 0.99175, 0.99176, 0.99177, 0.99178, 0.99179, 0.9918, 0.99181, 0.99182, 0.99183, 0.99184, 0.99185, 0.99186, 0.99187, 0.99188, 0.99189, 0.9919, 0.99191, 0.99192, 0.99193, 0.99194, 0.99195, 0.99196, 0.99197, 0.99198, 0.99199, 0.992, 0.99201, 0.99202, 0.99203, 0.99204, 0.99205, 0.99206, 0.99207, 0.99208, 0.99209, 0.9921, 0.99211, 0.99212, 0.99213, 0.99214, 0.99215, 0.99216, 0.99217, 0.99218, 0.99219, 0.9922, 0.99221, 0.99222, 0.99223, 0.99224, 0.99225, 0.99226, 0.99227, 0.99228, 0.99229, 0.9923, 0.99231, 0.99232, 0.99233, 0.99234, 0.99235, 0.99236, 0.99237, 0.99238, 0.99239, 0.9924, 0.99241, 0.99242, 0.99243, 0.99244, 0.99245, 0.99246, 0.99247, 0.99248, 0.99249, 0.99250, 0.99251, 0.99252, 0.99253, 0.99254, 0.99255, 0.99256, 0.99257, 0.99258, 0.99259, 0.9926, 0.99261, 0.99262, 0.99263, 0.99264, 0.99265, 0.99266, 0.99267, 0.99268, 0.99269, 0.9927, 0.99271, 0.99272, 0.99273, 0.99274, 0.99275, 0.99276, 0.99277, 0.99278, 0.99279, 0.9928, 0.99281, 0.99282, 0.99283, 0.99284, 0.99285, 0.99286, 0.99287, 0.99288, 0.99289, 0.9929, 0.99291, 0.99292, 0.99293, 0.99294, 0.99295, 0.99296, 0.99297, 0.99298, 0.99299, 0.993, 0.99301, 0.99302, 0.99303, 0.99304, 0.99305, 0.99306, 0.99307, 0.99308, 0.99309, 0.9931, 0.99311, 0.99312, 0.99313, 0.99314, 0.99315, 0.99316, 0.99317, 0.99318, 0.99319, 0.9932, 0.99321, 0.99322, 0.99323, 0.99324, 0.99325, 0.99326, 0.99327, 0.99328, 0.99329, 0.9933, 0.99331, 0.99332, 0.99333, 0.99334, 0.99335, 0.99336, 0.99337, 0.99338, 0.99339, 0.9934, 0.99341, 0.99342, 0.99343, 0.99344, 0.99345, 0.99346, 0.99347, 0.99348, 0.99349, 0.9935, 0.99351, 0.99352, 0.99353, 0.99354, 0.99355, 0.99356, 0.99357, 0.99358, 0.99359, 0.9936, 0.99361, 0.99362, 0.99363, 0.99364, 0.99365, 0.99366, 0.99367, 0.99368, 0.99369, 0.9937, 0.99371, 0.99372, 0.99373, 0.99374, 0.99375, 0.99376, 0.99377, 0.99378, 0.99379, 0.9938, 0.99381, 0.99382, 0.99383, 0.99384, 0.99385, 0.99386, 0.99387, 0.99388, 0.99389, 0.9939, 0.99391, 0.99392, 0.99393, 0.99394, 0.99395, 0.99396, 0.99397, 0.99398, 0.99399, 0.994, 0.99401, 0.99402, 0.99403, 0.99404, 0.99405, 0.99406, 0.99407, 0.99408, 0.99409, 0.9941, 0.99411, 0.99412, 0.99413, 0.99414, 0.99415, 0.99416, 0.99417, 0.99418, 0.99419, 0.9942, 0.99421, 0.99422, 0.99423, 0.99424, 0.99425, 0.99426, 0.99427, 0.99428, 0.99429, 0.9943, 0.99431, 0.99432, 0.99433, 0.99434, 0.99435, 0.99436, 0.99437, 0.99438, 0.99439, 0.9944, 0.99441, 0.99442, 0.99443, 0.99444, 0.99445, 0.99446, 0.99447, 0.99448, 0.99449, 0.99450, 0.99451, 0.99452, 0.99453, 0.99454, 0.99455, 0.99456, 0.99457, 0.99458, 0.99459, 0.9946, 0.99461, 0.99462, 0.99463, 0.99464, 0.99465, 0.99466, 0.99467, 0.99468, 0.99469, 0.9947, 0.99471, 0.99472, 0.99473, 0.99474, 0.99475, 0.99476, 0.99477, 0.99478, 0.99479, 0.9948, 0.99481, 0.99482, 0.99483, 0.99484, 0.99485, 0.99486, 0.99487, 0.99488, 0.99489, 0.9949, 0.99491, 0.99492, 0.99493, 0.99494, 0.99495, 0.99496, 0.99497, 0.99498, 0.99499, 0.995, 0.99501, 0.99502, 0.99503, 0.99504, 0.99505, 0.99506, 0.99507, 0.99508, 0.99509, 0.9951, 0.99511, 0.99512, 0.99513, 0.99514, 0.99515, 0.99516, 0.99517, 0.99518, 0.99519, 0.9952, 0.99521, 0.99522, 0.99523, 0.99524, 0.99525, 0.99526, 0.99527, 0.99528, 0.99529, 0.9953, 0.99531, 0.99532, 0.99533, 0.99534, 0.99535, 0.99536, 0.99537, 0.99538, 0.99539, 0.9954, 0.99541, 0.99542, 0.99543, 0.99544, 0.99545, 0.99546, 0.99547, 0.99548, 0.99549, 0.9955, 0.99551, 0.99552, 0.99553, 0.99554, 0.99555, 0.99556, 0.99557, 0.99558, 0.99559, 0.9956, 0.99561, 0.99562, 0.99563, 0.99564, 0.99565, 0.99566, 0.99567, 0.99568, 0.99569, 0.9957, 0.99571, 0.99572, 0.99573, 0.99574, 0.99575, 0.99576, 0.99577, 0.99578, 0.99579, 0.9958, 0.99581, 0.99582, 0.99583, 0.99584, 0.99585, 0.99586, 0.99587, 0.99588, 0.99589, 0.9959, 0.99591, 0.99592, 0.99593, 0.99594, 0.99595, 0.99596, 0.99597, 0.99598, 0.99599, 0.996, 0.99601, 0.99602, 0.99603, 0.99604, 0.99605, 0.99606, 0.99607, 0.99608, 0.99609, 0.9961, 0.99611, 0.99612, 0.99613, 0.99614, 0.99615, 0.99616, 0.99617, 0.99618, 0.99619, 0.9962, 0.99621, 0.99622, 0.99623, 0.99624, 0.99625, 0.99626, 0.99627, 0.99628, 0.99629, 0.9963, 0.99631, 0.99632, 0.99633, 0.99634, 0.99635, 0.99636, 0.99637, 0.99638, 0.99639, 0.9964, 0.99641, 0.99642, 0.99643, 0.99644, 0.99645, 0.99646, 0.99647, 0.99648, 0.99649, 0.99650, 0.99651, 0.99652, 0.99653, 0.99654, 0.99655, 0.99656, 0.99657, 0.99658, 0.99659, 0.9966, 0.99661, 0.99662, 0.99663, 0.99664, 0.99665, 0.99666, 0.99667, 0.99668, 0.99669, 0.9967, 0.99671, 0.99672, 0.99673, 0.99674, 0.99675, 0.99676, 0.99677, 0.99678, 0.99679, 0.9968, 0.99681, 0.99682, 0.99683, 0.99684, 0.99685, 0.99686, 0.99687, 0.99688, 0.99689, 0.9969, 0.99691, 0.99692, 0.99693, 0.99694, 0.99695, 0.99696, 0.99697, 0.99698, 0.99699, 0.997, 0.99701, 0.99702, 0.99703, 0.99704, 0.99705, 0.99706, 0.99707, 0.99708, 0.99709, 0.9971, 0.99711, 0.99712, 0.99713, 0.99714, 0.99715, 0.99716, 0.99717, 0.99718, 0.99719, 0.9972, 0.99721, 0.99722, 0.99723, 0.99724, 0.99725, 0.99726, 0.99727, 0.99728, 0.99729, 0.9973, 0.99731, 0.99732, 0.99733, 0.99734, 0.99735, 0.99736, 0.99737, 0.99738, 0.99739, 0.9974, 0.99741, 0.99742, 0.99743, 0.99744, 0.99745, 0.99746, 0.99747, 0.99748, 0.99749, 0.9975, 0.99751, 0.99752, 0.99753, 0.99754, 0.99755, 0.99756, 0.99757, 0.99758, 0.99759, 0.9976, 0.99761, 0.99762, 0.99763, 0.99764, 0.99765, 0.99766, 0.99767, 0.99768, 0.99769, 0.9977, 0.99771, 0.99772, 0.99773, 0.99774, 0.99775, 0.99776, 0.99777, 0.99778, 0.99779, 0.9978, 0.99781, 0.99782, 0.99783, 0.99784, 0.99785, 0.99786, 0.99787, 0.99788, 0.99789, 0.9979, 0.99791, 0.99792, 0.99793, 0.99794, 0.99795, 0.99796, 0.99797, 0.99798, 0.99799, 0.998, 0.99801, 0.99802, 0.99803, 0.99804, 0.99805, 0.99806, 0.99807, 0.99808, 0.99809, 0.9981, 0.99811, 0.99812, 0.99813, 0.99814, 0.99815, 0.99816, 0.99817, 0.99818, 0.99819, 0.9982, 0.99821, 0.99822, 0.99823, 0.99824, 0.99825, 0.99826, 0.99827, 0.99828, 0.99829, 0.9983, 0.99831, 0.99832, 0.99833, 0.99834, 0.99835, 0.99836, 0.99837, 0.99838, 0.99839, 0.9984, 0.99841, 0.99842, 0.99843, 0.99844, 0.99845, 0.99846, 0.99847, 0.99848, 0.99849, 0.99850, 0.99851, 0.99852, 0.99853, 0.99854, 0.99855, 0.99856, 0.99857, 0.99858, 0.99859, 0.9986, 0.99861, 0.99862, 0.99863, 0.99864, 0.99865, 0.99866, 0.99867, 0.99868, 0.99869, 0.9987, 0.99871, 0.99872, 0.99873, 0.99874, 0.99875, 0.99876, 0.99877, 0.99878, 0.99879, 0.9988, 0.99881, 0.99882, 0.99883, 0.99884, 0.99885, 0.99886, 0.99887, 0.99888, 0.99889, 0.9989, 0.99891, 0.99892, 0.99893, 0.99894, 0.99895, 0.99896, 0.99897, 0.99898, 0.99899, 0.999, 0.99901, 0.99902, 0.99903, 0.99904, 0.99905, 0.99906, 0.99907, 0.99908, 0.99909, 0.9991, 0.99911, 0.99912, 0.99913, 0.99914, 0.99915, 0.99916, 0.99917, 0.99918, 0.99919, 0.9992, 0.99921, 0.99922, 0.99923, 0.99924, 0.99925, 0.99926, 0.99927, 0.99928, 0.99929, 0.9993, 0.99931, 0.99932, 0.99933, 0.99934, 0.99935, 0.99936, 0.99937, 0.99938, 0.99939, 0.9994, 0.99941, 0.99942, 0.99943, 0.99944, 0.99945, 0.99946, 0.99947, 0.99948, 0.99949, 0.9995, 0.99951, 0.99952, 0.99953, 0.99954, 0.99955, 0.99956, 0.99957, 0.99958, 0.99959, 0.9996, 0.99961, 0.99962, 0.99963, 0.99964, 0.99965, 0.99966, 0.99967, 0.99968, 0.99969, 0.9997, 0.99971, 0.99972, 0.99973, 0.99974, 0.99975, 0.99976, 0.99977, 0.99978, 0.99979, 0.9998, 0.99981, 0.99982, 0.99983, 0.99984, 0.99985, 0.99986, 0.99987, 0.99988, 0.99989, 0.9999, 0.99991, 0.99992, 0.99993, 0.99994, 0.99995, 0.99996, 0.99997, 0.99998, 0.99999, 1.0};
  for(int k = 0; k < trials+1; k++){
    probcuts[k] = probcuts_temp[k];
  }
}

