void InitializeProbCuts();

void calculate_Dsbackground_MLindep();
void calculate_Dsbackground_MLindep_corr();
void calculate_Dsbackground_MLdep(Int_t ptbin, Int_t iscan);
void calculate_Signal_CombBackground_MLdep(Int_t ptbin, Int_t iscan, Float_t expBssig, Float_t &expBsbkg, Float_t &expBssgnf);
void SetStyleHisto(TH1D *h);

void fit_signal(TH1F* hsig, int j, int i, double &xmin, double &xmax, bool draw);
void calculate_background(TH1F* hbkg, TF1* f1, int j, int i, double xmin, double xmax, double &bkgcent, TString filenamefits, bool draw, bool finalscan, double &bkglow, double &bkghigh);
void extract_fonll(TString filnam, int j, double &fonllmin, double &fonllcent, double &fonllmax, bool draw);
void extract_TAMU(TString filn, int j, double &tamucent, bool draw);
void calculate_efficiency(TH1F* heff, int j, int i, double &effcent, double &erreffcent, bool draw);

void plot_expected_significance(TString filenamemass, TString filenameeff, TString filenamefits);

//------------------------------------------------------

//Files related to Bs
TFile* fmass;
TFile* feff;
TFile* ffits;
TFile* fBkgCorr;
TString filenameTAMU;
TString filnameFONLL;

//Files related to Ds + pi background
TFile* fbkgshape;
TFile* fExpDsSignal;
TFile* fExpDsSignalScaled;
TFile* ftamuDspr;
TFile* ftamuDsfd;
TFile* ffonllDs;
TFile* feffDs;
TFile* feffDsmatchBs;

Bool_t drawall = kTRUE;

//TGraphs for saving final result
TGraphAsymmErrors* grBkgHIJING;
TGraphAsymmErrors* grBkgFONLL;
TGraphAsymmErrors* grBkgHIJINGrel;
TGraphAsymmErrors* grBkgFONLLrel;

TGraphAsymmErrors* grSigPythia;
TGraphAsymmErrors* grSigFONLL;
TGraphAsymmErrors* grSigPythiarel;
TGraphAsymmErrors* grSigFONLLrel;

TH1F* hEffPythia;

//------------------------------------------------------

Double_t nEv335 = 852644;
TString nEv335String = "8.53e5";
Double_t nEvExpected = 8000000000;
TString nLumiExpectedString = "10 nb^{-1}";

Double_t TAA = 23.07 * 1e-3; //mb^-1 -> mub^-1 in which we put FONLL (for Bs)
Double_t TAA2 = 23.07 * 1e-9; //mb^-1 -> pb^-1 in which we put FONLL (for Ds)
Double_t gauss3sigmafactor = 0.9973;

//Fragmentation / BR things for Bs
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

//Fragmentation + BR things for Ds
//https://indico.cern.ch/event/873138/contributions/3823930/attachments/2019548/3376384/D2H_FONLL_predictions_feeddown.pdf
Double_t ff_bplus = 0.09 * 0.412 / (0.117 * 0.412 + 0.09 * 0.412 + 0.93 * 0.088 + 0.011 * 0.089);
Double_t ff_bzero = 0.117 * 0.412 / (0.117 * 0.412 + 0.09 * 0.412 + 0.93 * 0.088 + 0.011 * 0.089);
Double_t ff_lambdab = 0.011 * 0.089 / (0.117 * 0.412 + 0.09 * 0.412 + 0.93 * 0.088 + 0.011 * 0.089);
Double_t ff_bs = 0.93 * 0.088 / (0.117 * 0.412 + 0.09 * 0.412 + 0.93 * 0.088 + 0.011 * 0.089);

//------------------------------------------------------

const Int_t nptbins = 5;
Int_t ptbins[nptbins+1] = {0, 4, 8, 12, 16, 24};
Float_t ptbinsfl[nptbins+1] = {0, 4, 8, 12, 16, 24};
Double_t ptbinsdb[nptbins+1] = {0, 4, 8, 12, 16, 24};
Int_t rebin[nptbins] = {16, 16, 16, 16, 16};
Int_t sgnfmax = 10;

/*
//ITS2 05/05 (double training 5000).
//--FINAL VALUES--
//   0-4: 4980; 4-8: 4750-4850; 8-12: 4400; 12-16: 4050; 16-24: 4800
const Int_t trials = 5000;
Float_t probcuts[trials+1];
Int_t selbin[nptbins] = {4980, 4800, 4400, 4050, 4800};//{0, 0, 0, 0, 0}; //
Float_t preselML[nptbins] = {0.6, 0.7 ,0.5, 0.3, 0.2};
TString bkgfitoption[nptbins] = {"expo", "expo", "pol1", "pol1", "pol1"};
Bool_t useRealFit[nptbins] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
Bool_t bincountBkg[nptbins] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
Int_t opt = 3;
Int_t selbinDs[nptbins] = {selbin[0]/10 + 1, selbin[1]/10 + 1, selbin[2]/10 + 1, selbin[3]/10 + 1, selbin[4]/10 + 1}; //Need +1, because injected Ds starts with ML=0
*/


//ITS3 05/05 (double training 5000)
//--FINAL VALUES--
//   0-4: 4600 ; 4-8: 4800 ; 8-12: >4850 ; 12-16: 4800 ; 16-24: 4550
const Int_t trials = 5000;
Float_t probcuts[trials+1];
Int_t selbin[nptbins] = {4600, 4800, 4850, 4800, 4550};//{0, 0, 0, 0, 0}; //
Float_t preselML[nptbins] = {0.4,0.5,0.4,0.3,0.3};
TString bkgfitoption[nptbins] = {"expo", "expo", "pol1", "pol1", "pol1"};
Bool_t useRealFit[nptbins] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
Bool_t bincountBkg[nptbins] = {kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
Int_t opt = 3;
Int_t selbinDs[nptbins] = {selbin[0]/10 + 1, selbin[1]/10 + 1, selbin[2]/10 + 1, selbin[3]/10 + 1, selbin[4]/10 + 1}; //Need +1, because injected Ds starts with ML=0


void ExpectedSignificanceAnalysis(TString path = "", Int_t opt = 0, Bool_t drawessentials = kFALSE){

  if(drawessentials) drawall = kFALSE;

  TGaxis::SetMaxDigits(3);

  TString filenamemass = path + "masshisto.root";
  TString filenameeff = path + "effhisto.root";
  TString filenamefits = path + "parametrised_bkgpars.root";
  TString filenameBkgCorr = "theory/BkgCorrFactor_Bs_1DataFile_25MCFile_CoarsepTBinning.root";
  TString filenamebkgshape = path + "masshisto_bkgshape.root";
  TString filenameexpDssignal = path + "expected_Ds_signal_perevent_correct.root";
  TString filenameexpDssignalscaled = path + "expected_Ds_signal_fullyscaled.root";
  TString filenametamuDspr = "theory/input_RAA_TAMU_prDs.root";
  TString filenametamuDsfd = "theory/input_RAA_TAMU_fdDs.root";
  TString filenamefonllDs = "fonll/DmesonLcPredictions_sqrt5500_y05_pythia8_FFee_BRPDG.root";
  TString filenameeffDs = path + "efficiency_Ds_050_50MeVbins.root";
  TString filenameeffDsmatchBs = path + "BspT_for_50MeVDsbins_Tree.root";

  //Files related to Bs
  fmass = new TFile(filenamemass.Data());
  feff = new TFile(filenameeff.Data());
  ffits = new TFile(filenamefits.Data());
  fBkgCorr = new TFile(filenameBkgCorr.Data());
  filenameTAMU = "theory/input_RAA_TAMU_Bs.txt";
  filnameFONLL = "fonll/FONLL-Bhadron-ds-sqrts5500-CoarseBinning.txt";

  //Files related to Ds + pi background
  fbkgshape = new TFile(filenamebkgshape.Data());
  if(opt < 1) fExpDsSignal = new TFile(filenameexpDssignal.Data(), "RECREATE");
  else        fExpDsSignal = new TFile(filenameexpDssignal.Data());
  if(opt < 2) fExpDsSignalScaled = new TFile(filenameexpDssignalscaled.Data(), "RECREATE");
  else        fExpDsSignalScaled = new TFile(filenameexpDssignalscaled.Data());
  ftamuDspr = new TFile(filenametamuDspr.Data());
  ftamuDsfd = new TFile(filenametamuDsfd.Data());
  ffonllDs = new TFile(filenamefonllDs.Data());
  feffDs = new TFile(filenameeffDs.Data());
  feffDsmatchBs = new TFile(filenameeffDsmatchBs.Data());

  InitializeProbCuts();

/*
  for(int p = 0; p < 100; p++){
    selbin[0] = p*50;
    selbin[1] = p*50;
    selbin[2] = p*50;
    selbin[3] = p*50;
    selbin[4] = p*50;

    selbinDs[0] = selbin[0] / 10 + 1;
    selbinDs[1] = selbin[1] / 10 + 1;
    selbinDs[2] = selbin[2] / 10 + 1;
    selbinDs[3] = selbin[3] / 10 + 1;
    selbinDs[4] = selbin[4] / 10 + 1;

    //if(opt < 1) calculate_Dsbackground_MLindep();
    if(opt < 1 && p == 0) calculate_Dsbackground_MLindep_corr();
*/

  if(opt < 1) calculate_Dsbackground_MLindep_corr();

  if(opt < 2){
    //Need +1, because injected Ds starts with ML=0
    for(int iscan = 0; iscan <= trials/10; iscan++){
      for(int ptbin = 0; ptbin < nptbins; ptbin++){
        if(iscan != selbinDs[ptbin]) continue;
        calculate_Dsbackground_MLdep(ptbin, iscan);
      }
    }
  }

  Float_t expectedsignalcent[nptbins];
  Float_t expectedsignalfonllmin[nptbins];
  Float_t expectedsignalfonllmax[nptbins];
  Float_t expectedsignalerreffmin[nptbins];
  Float_t expectedsignalerreffmax[nptbins];
  Float_t expectedbkg[nptbins];
  Float_t expectedbkglow[nptbins];
  Float_t expectedbkghigh[nptbins];
  Float_t expectedsgnfcent[nptbins];
  Float_t expectedsgnffonllmin[nptbins];
  Float_t expectedsgnffonllmax[nptbins];
  Float_t expectedsgnferreffmin[nptbins];
  Float_t expectedsgnferreffmax[nptbins];
  Float_t expectedsgnfbkgmin[nptbins];
  Float_t expectedsgnfbkgmax[nptbins];

  grBkgHIJING = new TGraphAsymmErrors(nptbins);
  grBkgFONLL = new TGraphAsymmErrors(nptbins);
  grBkgHIJINGrel = new TGraphAsymmErrors(nptbins);
  grBkgFONLLrel = new TGraphAsymmErrors(nptbins);

  grSigPythia = new TGraphAsymmErrors(nptbins);
  grSigFONLL = new TGraphAsymmErrors(nptbins);
  grSigPythiarel = new TGraphAsymmErrors(nptbins);
  grSigFONLLrel = new TGraphAsymmErrors(nptbins);

  hEffPythia = new TH1F("hEffPythia", ";#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon)",nptbins,ptbinsfl);

  for(int iscan = 0; iscan <= trials; iscan++){
    for(int ptbin = 0; ptbin < nptbins; ptbin++){
      if(iscan != selbin[ptbin]) continue;

      cout << "Using input FONLL file = " << filnameFONLL << endl;
      Double_t fonllmin, fonllcent, fonllmax;
      extract_fonll(filnameFONLL, ptbin, fonllmin, fonllcent, fonllmax, kFALSE);
      cout << " Central FONLL [" << ptbins[ptbin] << "," << ptbins[ptbin+1] << "]: " << fonllcent << endl;

      cout << "Using input TAMU file = " << filenameTAMU << endl;
      Double_t tamucent;
      extract_TAMU(filenameTAMU, ptbin, tamucent, kFALSE);
      cout << " Central TAMU [" << ptbins[ptbin] << "," << ptbins[ptbin+1] << "]: " << tamucent << endl;

      TH1F* heff = (TH1F*)feff->Get(Form("eff_%d",iscan));
      cout << "Using input efficiency histograms = " << heff->GetName() << endl;
      Double_t effcent, erreffcent;
      calculate_efficiency(heff, ptbin, iscan, effcent, erreffcent, kFALSE);

      expectedsignalcent[ptbin] = 2 * (ptbins[ptbin+1] - ptbins[ptbin]) * 1 * (BRBs * BRDs) * nEvExpected * effcent * TAA * fonllcent * tamucent;
      expectedsignalcent[ptbin] *= gauss3sigmafactor;

      calculate_Signal_CombBackground_MLdep(ptbin, iscan, expectedsignalcent[ptbin], expectedbkg[ptbin], expectedsgnfcent[ptbin]);

      expectedsignalfonllmin[ptbin] = 2 * (ptbins[ptbin+1] - ptbins[ptbin]) * 1 * (BRBs * BRDs) * nEvExpected * effcent * TAA * fonllmin * tamucent;
      expectedsignalfonllmax[ptbin] = 2 * (ptbins[ptbin+1] - ptbins[ptbin]) * 1 * (BRBs * BRDs) * nEvExpected * effcent * TAA * fonllmax * tamucent;
      expectedsignalerreffmin[ptbin] = 2 * (ptbins[ptbin+1] - ptbins[ptbin]) * 1 * (BRBs * BRDs) * nEvExpected * (effcent - erreffcent) * TAA * fonllcent * tamucent;
      expectedsignalerreffmax[ptbin] = 2 * (ptbins[ptbin+1] - ptbins[ptbin]) * 1 * (BRBs * BRDs) * nEvExpected * (effcent + erreffcent) * TAA * fonllcent * tamucent;

      expectedsignalfonllmin[ptbin] *= gauss3sigmafactor;
      expectedsignalfonllmax[ptbin] *= gauss3sigmafactor;
      expectedsignalerreffmin[ptbin] *= gauss3sigmafactor;
      expectedsignalerreffmax[ptbin] *= gauss3sigmafactor;

      Double_t ptmin = ptbins[ptbin];
      Double_t ptmax = ptbins[ptbin+1];
      Double_t ptcent = 0.5*(ptmax + ptmin);
      Double_t pterr = 0.5*(ptmax - ptmin);
      grSigPythia->SetPoint(ptbin, ptcent, expectedsignalcent[ptbin]);
      grSigPythia->SetPointError(ptbin, pterr, pterr, expectedsignalcent[ptbin]-expectedsignalerreffmin[ptbin], expectedsignalerreffmax[ptbin]-expectedsignalcent[ptbin]);
      grSigFONLL->SetPoint(ptbin, ptcent, expectedsignalcent[ptbin]);
      grSigFONLL->SetPointError(ptbin, pterr, pterr, expectedsignalcent[ptbin]-expectedsignalfonllmin[ptbin], expectedsignalfonllmax[ptbin]-expectedsignalcent[ptbin]);

      grSigPythiarel->SetPoint(ptbin, ptcent, 0);
      grSigPythiarel->SetPointError(ptbin, pterr, pterr, (expectedsignalcent[ptbin]-expectedsignalerreffmin[ptbin])/expectedsignalcent[ptbin], (expectedsignalerreffmax[ptbin]-expectedsignalcent[ptbin])/expectedsignalcent[ptbin]);
      grSigFONLLrel->SetPoint(ptbin, ptcent, 0);
      grSigFONLLrel->SetPointError(ptbin, pterr, pterr, (expectedsignalcent[ptbin]-expectedsignalfonllmin[ptbin])/expectedsignalcent[ptbin], (expectedsignalfonllmax[ptbin]-expectedsignalcent[ptbin])/expectedsignalcent[ptbin]);

      hEffPythia->SetBinContent(hEffPythia->FindBin(ptcent),effcent);
      hEffPythia->SetBinError(hEffPythia->FindBin(ptcent),erreffcent);
    }
  }
  for(int ptbin = 0; ptbin < nptbins; ptbin++){
    cout << "[" << ptbins[ptbin] << "-" << ptbins[ptbin+1] << "], sig=" << expectedsignalcent[ptbin] << ", bkg=" << expectedbkg[ptbin] << ", sgnf=" << expectedsgnfcent[ptbin] << endl;
  }
  /*
  TCanvas* ctest = new TCanvas("ctest","ctest", 800,400);
  ctest->Divide(2);
  ctest->cd(1);
  grBkgHIJING->Draw("alp");
  grBkgFONLL->SetLineColor(kRed);
  grBkgFONLL->Draw("same");
  ctest->cd(2);
  grBkgHIJINGrel->Draw("alp");
  grBkgFONLLrel->SetLineColor(kRed);
  grBkgFONLLrel->Draw("same");

  TCanvas* ctest2 = new TCanvas("ctest2","ctest2", 800,400);
  ctest2->Divide(2);
  ctest2->cd(1);
  grSigFONLL->SetLineColor(kRed);
  grSigFONLL->Draw("alp");
  grSigPythia->Draw("same");
  ctest2->cd(2);
  grSigFONLLrel->SetLineColor(kRed);
  grSigFONLLrel->Draw("alp");
  grSigPythiarel->Draw("same");

  TCanvas* ctest3 = new TCanvas();
  ctest3->cd();
  hEffPythia->Draw("ep");

  TString foutname_sbe = path + "signal_background_efficiency.root";
  TFile* fout_sbe = new TFile(foutname_sbe.Data(), "RECREATE");
  fout_sbe->cd();
  grSigFONLL->Write("grSigFONLL");
  grSigPythia->Write("grSigPythia");
  grSigFONLLrel->Write("grSigFONLLrel");
  grSigPythiarel->Write("grSigPythiarel");

  grBkgFONLL->Write("grBkgFONLL");
  grBkgHIJING->Write("grBkgHIJING");
  grBkgFONLLrel->Write("grBkgFONLLrel");
  grBkgFONLLrel->Write("grBkgFONLLrel");

  hEffPythia->Write("hEffPythia");
  cout << "Wrote output in: " << foutname_sbe.Data() << endl;
  fout_sbe->Close();
  */
  plot_expected_significance(filenamemass, filenameeff, filenamefits);
//  }
}

void calculate_Dsbackground_MLindep_corr(){

  TH1F* htamu_pr = (TH1F*)ftamuDspr->Get("hTAMUcent");
  TH1F* htamu_fd = (TH1F*)ftamuDsfd->Get("hTAMUcent");

  TH1F* hfonll_pr = (TH1F*)ffonllDs->Get("hDsPhipitoKkpipred_central");
  TH1F* hfonll_fd = (TH1F*)ffonllDs->Get("hDsPhipitoKkpifromBpred_central_corr");
  TH1F* hfonllmin_pr = (TH1F*)ffonllDs->Get("hDsPhipitoKkpipred_min");
  TH1F* hfonllmin_fd = (TH1F*)ffonllDs->Get("hDsPhipitoKkpifromBpred_min_corr");
  TH1F* hfonllmax_pr = (TH1F*)ffonllDs->Get("hDsPhipitoKkpipred_max");
  TH1F* hfonllmax_fd = (TH1F*)ffonllDs->Get("hDsPhipitoKkpifromBpred_max_corr");

  TH1F* h_gen_pr = (TH1F*)feffDs->Get("h_gen_pr");
  TH1F* h_gen_fd = (TH1F*)feffDs->Get("h_gen_fd");
  TH1F* h_presel_pr = (TH1F*)feffDs->Get("h_presel_pr");
  TH1F* h_presel_fd = (TH1F*)feffDs->Get("h_presel_fd");

  //50 MeV bins, prediction from 0-50 GeV/c
  TH1F* h_match_pr[1000];
  TH1F* h_match_fd[1000];
  TH1F* h_match_pr_min[1000];
  TH1F* h_match_fd_min[1000];
  TH1F* h_match_pr_max[1000];
  TH1F* h_match_fd_max[1000];

  for(int i = 0; i < 1000; i++){
    h_match_pr[i] = (TH1F*)feffDsmatchBs->Get(Form("h_bs_pr_%d",i));
    h_match_pr_min[i] = (TH1F*)h_match_pr[i]->Clone(Form("h_bs_pr_min_%d",i));
    h_match_pr_max[i] = (TH1F*)h_match_pr[i]->Clone(Form("h_bs_pr_max_%d",i));
    h_match_fd[i] = (TH1F*)feffDsmatchBs->Get(Form("h_bs_fd_%d",i));
    h_match_fd_min[i] = (TH1F*)h_match_fd[i]->Clone(Form("h_bs_fd_min_%d",i));
    h_match_fd_max[i] = (TH1F*)h_match_fd[i]->Clone(Form("h_bs_fd_max_%d",i));
  }

  Double_t testNPDs12;

  for(Int_t i = 1; i <= hfonll_pr->GetNbinsX(); i++){
    Double_t pt = hfonll_pr->GetBinCenter(i);
    if(pt > 50) continue;
    Double_t dpt = hfonll_pr->GetBinWidth(i);
    Double_t yfonll = hfonll_pr->GetBinContent(i);
    Double_t yfonllmin = hfonllmin_pr->GetBinContent(i);
    Double_t yfonllmax = hfonllmax_pr->GetBinContent(i);
    Double_t ytamu = htamu_pr->GetBinContent(htamu_pr->FindBin(pt));
    Double_t dy = 1;
    Double_t acceffDs;
    if(h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt)) > 0){
      acceffDs = h_presel_pr->GetBinContent(h_presel_pr->FindBin(pt)) / h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt));
    } else {
      acceffDs = 0;
    }

    Double_t signal = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonll * ytamu;
    Double_t signalmin = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonllmin * ytamu;
    Double_t signalmax = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonllmax * ytamu;
    h_match_pr[i-1]->Scale(signal);
    h_match_pr_min[i-1]->Scale(signalmin);
    h_match_pr_max[i-1]->Scale(signalmax);

    if(pt >= 1 && pt < 2){
      //Test to non-prompt Ds ITS3 analysis, bin [1-2]
      //https://indico.cern.ch/event/893570/contributions/3769353/attachments/1999626/3337235/Ds_ITS3_estimates_2020_03_06.pdf#page=7
      //  Acc x Eff Bs TTree: 0.00565959
      //  Acc x Eff NP Ds: 0.002
      //    Exp Ev Fabrizio: "Indeed I always write 8, but I use 7.7"
      testNPDs12 += (0.002/0.00565959) * signal*7700000000.;
    }
  }

  for(Int_t i = 1; i <= hfonll_fd->GetNbinsX(); i++){
    Double_t pt = hfonll_fd->GetBinCenter(i);
    if(pt > 50) continue;
    Double_t dpt = hfonll_fd->GetBinWidth(i);
    Double_t yfonll = hfonll_fd->GetBinContent(i);
    Double_t yfonllmin = hfonllmin_fd->GetBinContent(i);
    Double_t yfonllmax = hfonllmax_fd->GetBinContent(i);
    Double_t ytamu = htamu_fd->GetBinContent(htamu_fd->FindBin(pt));
    Double_t dy = 1;
    Double_t acceffDs;
    if(h_gen_fd->GetBinContent(h_gen_fd->FindBin(pt)) > 0){
      acceffDs = h_presel_fd->GetBinContent(h_presel_fd->FindBin(pt)) / h_gen_fd->GetBinContent(h_gen_fd->FindBin(pt));
    } else {
      acceffDs = 0;
    }

    Double_t signal = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonll * ytamu;
    Double_t signalmin = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonllmin * ytamu;
    Double_t signalmax = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonllmax * ytamu;
    h_match_fd[i-1]->Scale(signal);
    h_match_fd_min[i-1]->Scale(signalmin);
    h_match_fd_max[i-1]->Scale(signalmax);

    if(pt >= 1 && pt < 2){
      //Test to non-prompt Ds ITS3 analysis, bin [1-2]
      //https://indico.cern.ch/event/893570/contributions/3769353/attachments/1999626/3337235/Ds_ITS3_estimates_2020_03_06.pdf#page=7
      //  Acc x Eff Bs TTree: 0.00523262
      //  Acc x Eff NP Ds: 0.006
      //    Exp Ev Fabrizio: "Indeed I always write 8, but I use 7.7"
      testNPDs12 += (0.006/0.00523262) * signal*7700000000.;
    }
  }

  cout << "\nCrossCheck expected Ds signal for pT [1-2] with NP Ds analysis: " << endl;
  cout << "  NP: S = f(Sgnf, S/B) = ((70*70)*0.9 + 70*70)/0.9 = 10344.444" << endl;
  cout << "  Bs: S (scaled to match efficiency NP analysis) = " << testNPDs12 << endl;
  cout << "    Should be very close (03/05/20)!" << endl << endl;

  for(int i = 1; i < 1000; i++){
    h_match_pr[0]->Add(h_match_pr[i]);
    h_match_pr_min[0]->Add(h_match_pr_min[i]);
    h_match_pr_max[0]->Add(h_match_pr_max[i]);
    h_match_fd[0]->Add(h_match_fd[i]);
    h_match_fd_min[0]->Add(h_match_fd_min[i]);
    h_match_fd_max[0]->Add(h_match_fd_max[i]);
  }

  TH1F* h_expsig_pr = (TH1F*)h_match_pr[0]->Rebin(nptbins, "h_expsig_pr", ptbinsdb);
  TH1F* h_expsigmin_pr = (TH1F*)h_match_pr_min[0]->Rebin(nptbins, "h_expsigmin_pr", ptbinsdb);
  TH1F* h_expsigmax_pr = (TH1F*)h_match_pr_max[0]->Rebin(nptbins, "h_expsigmax_pr", ptbinsdb);
  h_expsig_pr->SetTitle(";#it{p}_{T}(B_{s}) (GeV/it{c});Expected pr. signal per event");
  h_expsigmin_pr->SetTitle(";#it{p}_{T}(B_{s}) (GeV/it{c});Min. Expected pr. signal per event");
  h_expsigmax_pr->SetTitle(";#it{p}_{T}(B_{s}) (GeV/it{c});Max. Expected pr. signal per event");

  TH1F* h_expsig_fd = (TH1F*)h_match_fd[0]->Rebin(nptbins, "h_expsig_fd", ptbinsdb);
  TH1F* h_expsigmin_fd = (TH1F*)h_match_fd_min[0]->Rebin(nptbins, "h_expsigmin_fd", ptbinsdb);
  TH1F* h_expsigmax_fd = (TH1F*)h_match_fd_max[0]->Rebin(nptbins, "h_expsigmax_fd", ptbinsdb);
  h_expsig_fd->SetTitle(";#it{p}_{T}(B_{s}) (GeV/it{c});Expected fd. signal per event");
  h_expsigmin_fd->SetTitle(";#it{p}_{T}(B_{s}) (GeV/it{c});Min. Expected fd. signal per event");
  h_expsigmax_fd->SetTitle(";#it{p}_{T}(B_{s}) (GeV/it{c});Max. Expected fd. signal per event");

  fExpDsSignal->cd();
  h_match_pr[0]->Write();
  h_match_pr_min[0]->Write();
  h_match_pr_max[0]->Write();
  h_match_fd[0]->Write();
  h_match_fd_min[0]->Write();
  h_match_fd_max[0]->Write();
  h_expsig_pr->Write();
  h_expsigmin_pr->Write();
  h_expsigmax_pr->Write();
  h_expsig_fd->Write();
  h_expsigmin_fd->Write();
  h_expsigmax_fd->Write();

}

void calculate_Dsbackground_MLdep(Int_t ptbin, Int_t iscan){

  Int_t ptmin = ptbins[ptbin];
  Int_t ptmax = ptbins[ptbin+1];
  cout << "\nAnalysing: " << ptmin << " - " << ptmax << endl;
  TString suffixpt = Form("pt_cand%d_%d", ptmin, ptmax);
  TString suffix = Form("pt_cand%d_%d_%d", ptmin, ptmax, iscan);
  TString suffix0 = Form("pt_cand%d_%d_0", ptmin, ptmax);

  TH1F* h_expsig_pr = (TH1F*)fExpDsSignal->Get("h_expsig_pr");
  TH1F* h_expsigmin_pr = (TH1F*)fExpDsSignal->Get("h_expsigmin_pr");
  TH1F* h_expsigmax_pr = (TH1F*)fExpDsSignal->Get("h_expsigmax_pr");
  TH1F* h_expsig_fd = (TH1F*)fExpDsSignal->Get("h_expsig_fd");
  TH1F* h_expsigmin_fd = (TH1F*)fExpDsSignal->Get("h_expsigmin_fd");
  TH1F* h_expsigmax_fd = (TH1F*)fExpDsSignal->Get("h_expsigmax_fd");

  TH1F* hUniqueDs = (TH1F*)fbkgshape->Get(Form("hNormUniqueDs%s",suffixpt.Data()));
  TH1F* hdspr0 = (TH1F*)fbkgshape->Get(Form("hmass_DsPr%s", suffix0.Data()));
  TH1F* hdsfdbplus0 = (TH1F*)fbkgshape->Get(Form("hmass_DsFDBplus%s", suffix0.Data()));
  TH1F* hdsfdbzero0 = (TH1F*)fbkgshape->Get(Form("hmass_DsFDBzero%s", suffix0.Data()));
  TH1F* hdsfdlambdab0 = (TH1F*)fbkgshape->Get(Form("hmass_DsFDLambdab%s", suffix0.Data()));
  TH1F* hdsfdbs0 = (TH1F*)fbkgshape->Get(Form("hmass_DsFDBs%s", suffix0.Data()));

  TH1F* hdspr = (TH1F*)fbkgshape->Get(Form("hmass_DsPr%s", suffix.Data()));
  TH1F* hdsfdbplus = (TH1F*)fbkgshape->Get(Form("hmass_DsFDBplus%s", suffix.Data()));
  TH1F* hdsfdbzero = (TH1F*)fbkgshape->Get(Form("hmass_DsFDBzero%s", suffix.Data()));
  TH1F* hdsfdlambdab = (TH1F*)fbkgshape->Get(Form("hmass_DsFDLambdab%s", suffix.Data()));
  TH1F* hdsfdbs = (TH1F*)fbkgshape->Get(Form("hmass_DsFDBs%s", suffix.Data()));

  Double_t exp_pr_cent_perev = h_expsig_pr->GetBinContent(h_expsig_pr->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_pr_min_perev = h_expsigmin_pr->GetBinContent(h_expsigmin_pr->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_pr_max_perev = h_expsigmax_pr->GetBinContent(h_expsigmax_pr->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_fd_cent_perev = h_expsig_fd->GetBinContent(h_expsig_fd->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_fd_min_perev = h_expsigmin_fd->GetBinContent(h_expsigmin_fd->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_fd_max_perev = h_expsigmax_fd->GetBinContent(h_expsigmax_fd->FindBin(0.5 * (ptmin + ptmax)));

  //Not fully correct, now I'm using #N_{Ds}(unique in Bs TTree), but I think it should be N_{Ds}(selected Ds TTree)
  //   Difference are the Ds candidates that only form Bs cands that are *all* rejected by the prefilter cuts in TTreeCreator
  //   This would reduce the efficiency. I should match the Ds
  //Double_t effBs_dspr = hdspr->GetEntries()/hUniqueDs->GetBinContent(1);
  //Double_t effBs_dsfdbplus = hdsfdbplus->GetEntries()/hUniqueDs->GetBinContent(2) * ff_bplus;
  //Double_t effBs_dsfdbzero = hdsfdbzero->GetEntries()/hUniqueDs->GetBinContent(3) * ff_bzero;
  //Double_t effBs_dsfdlambdab = hdsfdlambdab->GetEntries()/hUniqueDs->GetBinContent(4) * ff_lambdab;
  //Double_t effBs_dsfdbs = hdsfdbs->GetEntries()/hUniqueDs->GetBinContent(5) * ff_bs;
  //cout << "  Efficiency Ds -> Bs for isc = " << iscan*10 << " is: " << effBs_dspr << " " << effBs_dsfdbplus << " " << effBs_dsfdbzero << " " << effBs_dsfdlambdab << " " << effBs_dsfdbs << endl;
  //
  //This should be fixed now (04/05/20)... Look's ok

  Double_t effBs_dspr = hdspr->GetEntries()/hdspr0->GetEntries();
  Double_t effBs_dsfdbplus = hdsfdbplus->GetEntries()/hdsfdbplus0->GetEntries() * ff_bplus;
  Double_t effBs_dsfdbzero = hdsfdbzero->GetEntries()/hdsfdbzero0->GetEntries() * ff_bzero;
  Double_t effBs_dsfdlambdab = hdsfdlambdab->GetEntries()/hdsfdlambdab0->GetEntries() * ff_lambdab;
  Double_t effBs_dsfdbs = hdsfdbs->GetEntries()/hdsfdbs0->GetEntries() * ff_bs;
  cout << "  Efficiency Ds -> Bs (+FF factor for beauty) for isc = " << iscan*10 << " is: " << effBs_dspr << " " << effBs_dsfdbplus << " " << effBs_dsfdbzero << " " << effBs_dsfdlambdab << " " << effBs_dsfdbs << endl;

  Double_t exp_Bspr_cent_perev = effBs_dspr * exp_pr_cent_perev;
  Double_t exp_Bspr_min_perev =  effBs_dspr * exp_pr_min_perev;
  Double_t exp_Bspr_max_perev =  effBs_dspr * exp_pr_max_perev;
  Double_t exp_Bsfdbplus_cent_perev = effBs_dsfdbplus * exp_fd_cent_perev;
  Double_t exp_Bsfdbplus_min_perev =  effBs_dsfdbplus * exp_fd_min_perev;
  Double_t exp_Bsfdbplus_max_perev =  effBs_dsfdbplus * exp_fd_max_perev;
  Double_t exp_Bsfdbzero_cent_perev = effBs_dsfdbzero * exp_fd_cent_perev;
  Double_t exp_Bsfdbzero_min_perev =  effBs_dsfdbzero * exp_fd_min_perev;
  Double_t exp_Bsfdbzero_max_perev =  effBs_dsfdbzero * exp_fd_max_perev;
  Double_t exp_Bsfdlambdab_cent_perev = effBs_dsfdlambdab * exp_fd_cent_perev;
  Double_t exp_Bsfdlambdab_min_perev =  effBs_dsfdlambdab * exp_fd_min_perev;
  Double_t exp_Bsfdlambdab_max_perev =  effBs_dsfdlambdab * exp_fd_max_perev;
  Double_t exp_Bsfdbs_cent_perev = effBs_dsfdbs * exp_fd_cent_perev;
  Double_t exp_Bsfdbs_min_perev =  effBs_dsfdbs * exp_fd_min_perev;
  Double_t exp_Bsfdbs_max_perev =  effBs_dsfdbs * exp_fd_max_perev;
  cout << "  Expected (pr. Ds) + pi -> Bs cand: " << exp_Bspr_cent_perev << " + " << exp_Bspr_max_perev - exp_Bspr_cent_perev << " - " << exp_Bspr_cent_perev - exp_Bspr_min_perev << endl;
  cout << "  Expected (B+ fd. Ds) + pi -> Bs cand: " << exp_Bsfdbplus_cent_perev << " + " << exp_Bsfdbplus_max_perev - exp_Bsfdbplus_cent_perev << " - " << exp_Bsfdbplus_cent_perev - exp_Bsfdbplus_min_perev << endl;
  cout << "  Expected (B0 fd. Ds) + pi -> Bs cand: " << exp_Bsfdbzero_cent_perev << " + " << exp_Bsfdbzero_max_perev - exp_Bsfdbzero_cent_perev << " - " << exp_Bsfdbzero_cent_perev - exp_Bsfdbzero_min_perev << endl;
  cout << "  Expected (Lb fd. Ds) + pi -> Bs cand: " << exp_Bsfdlambdab_cent_perev << " + " << exp_Bsfdlambdab_max_perev - exp_Bsfdlambdab_cent_perev << " - " << exp_Bsfdlambdab_cent_perev - exp_Bsfdlambdab_min_perev << endl;
  cout << "  Expected (Bs fd. Ds) + pi -> Bs cand: " << exp_Bsfdbs_cent_perev << " + " << exp_Bsfdbs_max_perev - exp_Bsfdbs_cent_perev << " - " << exp_Bsfdbs_cent_perev - exp_Bsfdbs_min_perev << endl;

  TF1* fdspr2 = (TF1*)fbkgshape->Get(Form("fbkg_dspr%s", suffix.Data()));
  TF1* fdsfdbzero2 = (TF1*)fbkgshape->Get(Form("fbkg_dsfdbzero%s", suffix.Data()));
  TF1* fdsfdbplus2 = (TF1*)fbkgshape->Get(Form("fbkg_dsfdbplus%s", suffix.Data()));
  TF1* fdsfdlambdab2 = (TF1*)fbkgshape->Get(Form("fbkg_dsfdlambdab%s", suffix.Data()));
  TF1* fdsfdbs2 = (TF1*)fbkgshape->Get(Form("fbkg_dsfdbs%s", suffix.Data()));
  fdspr2->SetLineWidth(2);
  fdsfdbzero2->SetLineWidth(2);
  fdsfdbplus2->SetLineWidth(2);
  fdsfdlambdab2->SetLineWidth(2);
  fdsfdbs2->SetLineWidth(2);

  TF1* fdspr1 = (TF1*)fbkgshape->Get(Form("fbkg_dspr1%s", suffix.Data()));
  TF1* fdsfdbzero1 = (TF1*)fbkgshape->Get(Form("fbkg_dsfdbzero1%s", suffix.Data()));
  TF1* fdsfdbplus1 = (TF1*)fbkgshape->Get(Form("fbkg_dsfdbplus1%s", suffix.Data()));
  TF1* fdsfdlambdab1 = (TF1*)fbkgshape->Get(Form("fbkg_dsfdlambdab1%s", suffix.Data()));
  TF1* fdsfdbs1 = (TF1*)fbkgshape->Get(Form("fbkg_dsfdbs1%s", suffix.Data()));
  fdspr1->SetLineWidth(2);
  fdsfdbzero1->SetLineWidth(2);
  fdsfdbplus1->SetLineWidth(2);
  fdsfdlambdab1->SetLineWidth(2);
  fdsfdbs1->SetLineWidth(2);

  TH1D *href = (TH1D*)fbkgshape->Get(Form("h_norm_dspr%s", suffix.Data()));
  href->Reset("ICEMS");
  href->SetTitle(Form("%.1f < #it{p}_{T} < %.1f (GeV/#it{c})", (float)ptmin, (float)ptmax));
  href->SetMinimum(0.);
  href->SetMaximum(0.05);
  SetStyleHisto(href);

  TCanvas* cpol1 = new TCanvas(Form("cpol1_%d_%d_%d",ptmin,ptmax,iscan*10), Form("cpol1_%d_%d_%d",ptmin,ptmax,iscan*10), 450, 400);
  cpol1->cd();
  href->Draw();
  gPad->SetTickx();
  gPad->SetTicky();
  fdspr1->Draw("same");
  fdsfdbzero1->Draw("same");
  fdsfdbplus1->Draw("same");
  fdsfdbs1->Draw("same");
  fdsfdlambdab1->Draw("same");

  TLegend* leg = new TLegend(0.5, 0.63, 0.92, 0.88, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(16); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->SetHeader("Injected D_{s}^{+} with HIJING #pi^{-}");
  leg->AddEntry(fdspr1, "(Prompt)   #scale[0.6]{ }D_{s} + #pi", "l");
  leg->AddEntry(fdsfdbzero1, "(B^{0} #rightarrow X +) D_{s} + #pi", "l");
  leg->AddEntry(fdsfdbplus1, "(B^{+} #rightarrow X +) D_{s} + #pi", "l");
  leg->AddEntry(fdsfdbs1, "(B_{s} #rightarrow X +) D_{s} + #pi", "l");
  leg->AddEntry(fdsfdlambdab1, "(#Lambda_{b} #rightarrow X +) D_{s} + #pi", "l");
  leg->Draw();

  TLatex info; info.SetNDC(); info.SetTextFont(43); info.SetTextSize(23);
  if(iscan == 0){
    info.DrawLatex(0.19, 0.14, "No ML selection");
  } else {
    info.DrawLatex(0.19, 0.14, Form("ML_{1} > %.3f + ML_{2} > %.3f", preselML[ptbin], probcuts[(iscan-1)*10]));
  }

  TF1* fdspr1norm = new TF1(Form("fdspr1norm%d_%d",ptbin,iscan),"[0]",fdspr1->GetXmin(), fdspr1->GetXmax());
  TF1* fdspr1normmin = (TF1*)fdspr1norm->Clone(Form("fdspr1normmin%d_%d",ptbin,iscan));
  TF1* fdspr1normmax = (TF1*)fdspr1norm->Clone(Form("fdspr1normmax%d_%d",ptbin,iscan));
  fdspr1norm->SetParameter(0,exp_Bspr_cent_perev * nEvExpected);
  fdspr1normmin->SetParameter(0,exp_Bspr_min_perev * nEvExpected);
  fdspr1normmax->SetParameter(0,exp_Bspr_max_perev * nEvExpected);
  cout << "   Setting par0 for fdspr1norm = " << exp_Bspr_cent_perev * nEvExpected << endl;
  TF1* fdspr1func = new TF1(Form("fdspr1func%d_%d",ptbin,iscan),Form("%f+%f*x",fdspr1->GetParameter(0),fdspr1->GetParameter(1)),fdspr1->GetXmin(), fdspr1->GetXmax());
  TF1* fdspr1scaled = new TF1(Form("fdspr1scaled%d_%d",ptbin,iscan),Form("%s*%s",fdspr1norm->GetName(), fdspr1func->GetName()), fdspr1->GetXmin(), fdspr1->GetXmax());
  TF1* fdspr1scaledmin = new TF1(Form("fdspr1scaledmin%d_%d",ptbin,iscan),Form("%s*%s",fdspr1normmin->GetName(), fdspr1func->GetName()), fdspr1->GetXmin(), fdspr1->GetXmax());
  TF1* fdspr1scaledmax = new TF1(Form("fdspr1scaledmax%d_%d",ptbin,iscan),Form("%s*%s",fdspr1normmax->GetName(), fdspr1func->GetName()), fdspr1->GetXmin(), fdspr1->GetXmax());
  fdspr1scaled->SetLineColor(fdspr1->GetLineColor());
  fdspr1scaled->SetLineWidth(fdspr1->GetLineWidth());
  fdspr1scaledmin->SetLineColor(fdspr1->GetLineColor());
  fdspr1scaledmin->SetLineWidth(fdspr1->GetLineWidth());
  fdspr1scaledmin->SetLineStyle(2);
  fdspr1scaledmax->SetLineColor(fdspr1->GetLineColor());
  fdspr1scaledmax->SetLineWidth(fdspr1->GetLineWidth());
  fdspr1scaledmax->SetLineStyle(2);

  //37 bins in 0.6
  //Double_t prefac = 37/0.6;
  //cout << "   Setting par0 for fdspr1norm = " << exp_Bspr_cent_perev * nEvExpected << endl;
  //cout << "   Setting par0 for fdspr1norm = " << exp_Bspr_min_perev * nEvExpected << endl;
  //cout << "   Setting par0 for fdspr1norm = " << exp_Bspr_max_perev * nEvExpected << endl;
  //cout << "Check integral: " << prefac*fdspr1func->Integral(fdspr1->GetXmin(), fdspr1->GetXmax()) << " ";
  //cout << "  " << prefac*fdspr1scaled->Integral(fdspr1->GetXmin(), fdspr1->GetXmax()) << " ";
  //cout << "  " << prefac*fdspr1scaledmin->Integral(fdspr1->GetXmin(), fdspr1->GetXmax()) << " ";
  //cout << "  " << prefac*fdspr1scaledmax->Integral(fdspr1->GetXmin(), fdspr1->GetXmax()) << endl;
  //==> Correct, although the unc functions look a bit weird in one plot, but this is as expected:
    //TF1* f1 = new TF1("f1",Form("%f - %f * x",10.,1.5),0,5)
    //TH1F* h100 = new TH1F("h100","",1000,0,5)
    //TH1F* h200 = new TH1F("h200","",1000,0,5)
    //h100->FillRandom("f1",100000)
    //h200->FillRandom("f1",200000)
    //h100->Draw("hist")
    //h200->Draw("same")

  TF1* fdsfdbzero1norm = new TF1(Form("fdsfdbzero1norm%d_%d",ptbin,iscan),"[0]",fdsfdbzero1->GetXmin(), fdsfdbzero1->GetXmax());
  TF1* fdsfdbzero1normmin = (TF1*)fdsfdbzero1norm->Clone(Form("fdsfdbzero1normmin%d_%d",ptbin,iscan));
  TF1* fdsfdbzero1normmax = (TF1*)fdsfdbzero1norm->Clone(Form("fdsfdbzero1normmax%d_%d",ptbin,iscan));
  fdsfdbzero1norm->SetParameter(0,exp_Bsfdbzero_cent_perev * nEvExpected);
  fdsfdbzero1normmin->SetParameter(0,exp_Bsfdbzero_min_perev * nEvExpected);
  fdsfdbzero1normmax->SetParameter(0,exp_Bsfdbzero_max_perev * nEvExpected);
  cout << "   Setting par0 for fdsfdbzero1norm = " << exp_Bsfdbzero_cent_perev * nEvExpected << endl;
  TF1* fdsfdbzero1func = new TF1(Form("fdsfdbzero1func%d_%d",ptbin,iscan),Form("%f+%f*x",fdsfdbzero1->GetParameter(0),fdsfdbzero1->GetParameter(1)),fdsfdbzero1->GetXmin(), fdsfdbzero1->GetXmax());
  TF1* fdsfdbzero1scaled = new TF1(Form("fdsfdbzero1scaled%d_%d",ptbin,iscan),Form("%s*%s",fdsfdbzero1norm->GetName(), fdsfdbzero1func->GetName()), fdsfdbzero1->GetXmin(), fdsfdbzero1->GetXmax());
  TF1* fdsfdbzero1scaledmin = new TF1(Form("fdsfdbzero1scaledmin%d_%d",ptbin,iscan),Form("%s*%s",fdsfdbzero1normmin->GetName(), fdsfdbzero1func->GetName()), fdsfdbzero1->GetXmin(), fdsfdbzero1->GetXmax());
  TF1* fdsfdbzero1scaledmax = new TF1(Form("fdsfdbzero1scaledmax%d_%d",ptbin,iscan),Form("%s*%s",fdsfdbzero1normmax->GetName(), fdsfdbzero1func->GetName()), fdsfdbzero1->GetXmin(), fdsfdbzero1->GetXmax());
  fdsfdbzero1scaled->SetLineColor(fdsfdbzero1->GetLineColor());
  fdsfdbzero1scaled->SetLineWidth(fdsfdbzero1->GetLineWidth());
  fdsfdbzero1scaledmin->SetLineColor(fdsfdbzero1->GetLineColor());
  fdsfdbzero1scaledmin->SetLineWidth(fdsfdbzero1->GetLineWidth());
  fdsfdbzero1scaledmin->SetLineStyle(2);
  fdsfdbzero1scaledmax->SetLineColor(fdsfdbzero1->GetLineColor());
  fdsfdbzero1scaledmax->SetLineWidth(fdsfdbzero1->GetLineWidth());
  fdsfdbzero1scaledmax->SetLineStyle(2);

  TF1* fdsfdbplus1norm = new TF1(Form("fdsfdbplus1norm%d_%d",ptbin,iscan),"[0]",fdsfdbplus1->GetXmin(), fdsfdbplus1->GetXmax());
  TF1* fdsfdbplus1normmin = (TF1*)fdsfdbplus1norm->Clone(Form("fdsfdbplus1normmin%d_%d",ptbin,iscan));
  TF1* fdsfdbplus1normmax = (TF1*)fdsfdbplus1norm->Clone(Form("fdsfdbplus1normmax%d_%d",ptbin,iscan));
  fdsfdbplus1norm->SetParameter(0,exp_Bsfdbplus_cent_perev * nEvExpected);
  fdsfdbplus1normmin->SetParameter(0,exp_Bsfdbplus_min_perev * nEvExpected);
  fdsfdbplus1normmax->SetParameter(0,exp_Bsfdbplus_max_perev * nEvExpected);
  cout << "   Setting par0 for fdsfdbplus1norm = " << exp_Bsfdbplus_cent_perev * nEvExpected << endl;
  TF1* fdsfdbplus1func = new TF1(Form("fdsfdbplus1func%d_%d",ptbin,iscan),Form("%f+%f*x",fdsfdbplus1->GetParameter(0),fdsfdbplus1->GetParameter(1)),fdsfdbplus1->GetXmin(), fdsfdbplus1->GetXmax());
  TF1* fdsfdbplus1scaled = new TF1(Form("fdsfdbplus1scaled%d_%d",ptbin,iscan),Form("%s*%s",fdsfdbplus1norm->GetName(), fdsfdbplus1func->GetName()), fdsfdbplus1->GetXmin(), fdsfdbplus1->GetXmax());
  TF1* fdsfdbplus1scaledmin = new TF1(Form("fdsfdbplus1scaledmin%d_%d",ptbin,iscan),Form("%s*%s",fdsfdbplus1normmin->GetName(), fdsfdbplus1func->GetName()), fdsfdbplus1->GetXmin(), fdsfdbplus1->GetXmax());
  TF1* fdsfdbplus1scaledmax = new TF1(Form("fdsfdbplus1scaledmax%d_%d",ptbin,iscan),Form("%s*%s",fdsfdbplus1normmax->GetName(), fdsfdbplus1func->GetName()), fdsfdbplus1->GetXmin(), fdsfdbplus1->GetXmax());
  fdsfdbplus1scaled->SetLineColor(fdsfdbplus1->GetLineColor());
  fdsfdbplus1scaled->SetLineWidth(fdsfdbplus1->GetLineWidth());
  fdsfdbplus1scaledmin->SetLineColor(fdsfdbplus1->GetLineColor());
  fdsfdbplus1scaledmin->SetLineWidth(fdsfdbplus1->GetLineWidth());
  fdsfdbplus1scaledmin->SetLineStyle(2);
  fdsfdbplus1scaledmax->SetLineColor(fdsfdbplus1->GetLineColor());
  fdsfdbplus1scaledmax->SetLineWidth(fdsfdbplus1->GetLineWidth());
  fdsfdbplus1scaledmax->SetLineStyle(2);

  TF1* fdsfdlambdab1norm = new TF1(Form("fdsfdlambdab1norm%d_%d",ptbin,iscan),"[0]",fdsfdlambdab1->GetXmin(), fdsfdlambdab1->GetXmax());
  TF1* fdsfdlambdab1normmin = (TF1*)fdsfdlambdab1norm->Clone(Form("fdsfdlambdab1normmin%d_%d",ptbin,iscan));
  TF1* fdsfdlambdab1normmax = (TF1*)fdsfdlambdab1norm->Clone(Form("fdsfdlambdab1normmax%d_%d",ptbin,iscan));
  fdsfdlambdab1norm->SetParameter(0,exp_Bsfdlambdab_cent_perev * nEvExpected);
  fdsfdlambdab1normmin->SetParameter(0,exp_Bsfdlambdab_min_perev * nEvExpected);
  fdsfdlambdab1normmax->SetParameter(0,exp_Bsfdlambdab_max_perev * nEvExpected);
  cout << "   Setting par0 for fdsfdlambdab1norm = " << exp_Bsfdlambdab_cent_perev * nEvExpected << endl;
  TF1* fdsfdlambdab1func = new TF1(Form("fdsfdlambdab1func%d_%d",ptbin,iscan),Form("%f+%f*x",fdsfdlambdab1->GetParameter(0),fdsfdlambdab1->GetParameter(1)),fdsfdlambdab1->GetXmin(), fdsfdlambdab1->GetXmax());
  TF1* fdsfdlambdab1scaled = new TF1(Form("fdsfdlambdab1scaled%d_%d",ptbin,iscan),Form("%s*%s",fdsfdlambdab1norm->GetName(), fdsfdlambdab1func->GetName()), fdsfdlambdab1->GetXmin(), fdsfdlambdab1->GetXmax());
  TF1* fdsfdlambdab1scaledmin = new TF1(Form("fdsfdlambdab1scaledmin%d_%d",ptbin,iscan),Form("%s*%s",fdsfdlambdab1normmin->GetName(), fdsfdlambdab1func->GetName()), fdsfdlambdab1->GetXmin(), fdsfdlambdab1->GetXmax());
  TF1* fdsfdlambdab1scaledmax = new TF1(Form("fdsfdlambdab1scaledmax%d_%d",ptbin,iscan),Form("%s*%s",fdsfdlambdab1normmax->GetName(), fdsfdlambdab1func->GetName()), fdsfdlambdab1->GetXmin(), fdsfdlambdab1->GetXmax());
  fdsfdlambdab1scaled->SetLineColor(fdsfdlambdab1->GetLineColor());
  fdsfdlambdab1scaled->SetLineWidth(fdsfdlambdab1->GetLineWidth());
  fdsfdlambdab1scaledmin->SetLineColor(fdsfdlambdab1->GetLineColor());
  fdsfdlambdab1scaledmin->SetLineWidth(fdsfdlambdab1->GetLineWidth());
  fdsfdlambdab1scaledmin->SetLineStyle(2);
  fdsfdlambdab1scaledmax->SetLineColor(fdsfdlambdab1->GetLineColor());
  fdsfdlambdab1scaledmax->SetLineWidth(fdsfdlambdab1->GetLineWidth());
  fdsfdlambdab1scaledmax->SetLineStyle(2);

  TF1* fdsfdbs1norm = new TF1(Form("fdsfdbs1norm%d_%d",ptbin,iscan),"[0]",fdsfdbs1->GetXmin(), fdsfdbs1->GetXmax());
  TF1* fdsfdbs1normmin = (TF1*)fdsfdbs1norm->Clone(Form("fdsfdbs1normmin%d_%d",ptbin,iscan));
  TF1* fdsfdbs1normmax = (TF1*)fdsfdbs1norm->Clone(Form("fdsfdbs1normmax%d_%d",ptbin,iscan));
  fdsfdbs1norm->SetParameter(0,exp_Bsfdbs_cent_perev * nEvExpected);
  fdsfdbs1normmin->SetParameter(0,exp_Bsfdbs_min_perev * nEvExpected);
  fdsfdbs1normmax->SetParameter(0,exp_Bsfdbs_max_perev * nEvExpected);
  cout << "   Setting par0 for fdsfdbs1norm = " << exp_Bsfdbs_cent_perev * nEvExpected << endl;
  TF1* fdsfdbs1func = new TF1(Form("fdsfdbs1func%d_%d",ptbin,iscan),Form("%f+%f*x",fdsfdbs1->GetParameter(0),fdsfdbs1->GetParameter(1)),fdsfdbs1->GetXmin(), fdsfdbs1->GetXmax());
  TF1* fdsfdbs1scaled = new TF1(Form("fdsfdbs1scaled%d_%d",ptbin,iscan),Form("%s*%s",fdsfdbs1norm->GetName(), fdsfdbs1func->GetName()), fdsfdbs1->GetXmin(), fdsfdbs1->GetXmax());
  TF1* fdsfdbs1scaledmin = new TF1(Form("fdsfdbs1scaledmin%d_%d",ptbin,iscan),Form("%s*%s",fdsfdbs1normmin->GetName(), fdsfdbs1func->GetName()), fdsfdbs1->GetXmin(), fdsfdbs1->GetXmax());
  TF1* fdsfdbs1scaledmax = new TF1(Form("fdsfdbs1scaledmax%d_%d",ptbin,iscan),Form("%s*%s",fdsfdbs1normmax->GetName(), fdsfdbs1func->GetName()), fdsfdbs1->GetXmin(), fdsfdbs1->GetXmax());
  fdsfdbs1scaled->SetLineColor(fdsfdbs1->GetLineColor());
  fdsfdbs1scaled->SetLineWidth(fdsfdbs1->GetLineWidth());
  fdsfdbs1scaledmin->SetLineColor(fdsfdbs1->GetLineColor());
  fdsfdbs1scaledmin->SetLineWidth(fdsfdbs1->GetLineWidth());
  fdsfdbs1scaledmin->SetLineStyle(2);
  fdsfdbs1scaledmax->SetLineColor(fdsfdbs1->GetLineColor());
  fdsfdbs1scaledmax->SetLineWidth(fdsfdbs1->GetLineWidth());
  fdsfdbs1scaledmax->SetLineStyle(2);

  TH1D *href2 = (TH1D*)fbkgshape->Get(Form("h_norm_dsfdbzero%s", suffix.Data()));
  href2->Reset("ICEMS");
  href2->SetTitle(Form("%.1f < #it{p}_{T} < %.1f (GeV/#it{c})", (float)ptmin, (float)ptmax));
  href2->GetYaxis()->SetTitle(href->GetYaxis()->GetTitle());
  href2->GetXaxis()->SetTitle(href->GetXaxis()->GetTitle());
  href2->SetMinimum(0.01);
  href2->SetMaximum(100000);
  SetStyleHisto(href2);

  TCanvas* cpol1scaled = new TCanvas(Form("cpol1scaled_%d_%d_%d",ptmin,ptmax,iscan*10), Form("cpol1scaled_%d_%d_%d",ptmin,ptmax,iscan*10), 450, 400);
  cpol1scaled->cd();
  href2->Draw();
  gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  fdspr1scaled->Draw("same");
  fdspr1scaledmin->Draw("same");
  fdspr1scaledmax->Draw("same");
  fdsfdbzero1scaled->Draw("same");
  fdsfdbzero1scaledmin->Draw("same");
  fdsfdbzero1scaledmax->Draw("same");
  fdsfdbplus1scaled->Draw("same");
  fdsfdbplus1scaledmin->Draw("same");
  fdsfdbplus1scaledmax->Draw("same");
  fdsfdbs1scaled->Draw("same");
  fdsfdbs1scaledmin->Draw("same");
  fdsfdbs1scaledmax->Draw("same");
  fdsfdlambdab1scaled->Draw("same");
  fdsfdlambdab1scaledmin->Draw("same");
  fdsfdlambdab1scaledmax->Draw("same");

  TLegend* legscaled = new TLegend(0.5, 0.63, 0.92, 0.88, 0, "NDC");
  legscaled->SetTextFont(43); legscaled->SetTextSize(16); legscaled->SetFillColor(0); legscaled->SetFillStyle(0); legscaled->SetLineColor(0);
  legscaled->SetHeader("Injected D_{s}^{+} with HIJING #pi^{-}");
  legscaled->AddEntry(fdspr1scaled, "(Prompt)   #scale[0.6]{ }D_{s} + #pi", "l");
  legscaled->AddEntry(fdsfdbzero1scaled, "(B^{0} #rightarrow X +) D_{s} + #pi", "l");
  legscaled->AddEntry(fdsfdbplus1scaled, "(B^{+} #rightarrow X +) D_{s} + #pi", "l");
  legscaled->AddEntry(fdsfdbs1scaled, "(B_{s} #rightarrow X +) D_{s} + #pi", "l");
  legscaled->AddEntry(fdsfdlambdab1scaled, "(#Lambda_{b} #rightarrow X +) D_{s} + #pi", "l");
  legscaled->Draw();

  TLatex infoscaled; infoscaled.SetNDC(); infoscaled.SetTextFont(43); infoscaled.SetTextSize(23);
  if(iscan == 0){
    infoscaled.DrawLatex(0.19, 0.14, "No ML selection");
  } else {
    infoscaled.DrawLatex(0.19, 0.14, Form("ML_{1} > %.3f + ML_{2} > %.3f", preselML[ptbin], probcuts[(iscan - 1)*10]));
  }

  fExpDsSignalScaled->cd();
  cpol1->Write();
  cpol1scaled->Write();

  fdspr1scaled->Write();
  fdsfdbzero1scaled->Write();
  fdsfdbplus1scaled->Write();
  fdsfdbs1scaled->Write();
  fdsfdlambdab1scaled->Write();

  fdspr1scaledmin->Write();
  fdspr1scaledmax->Write();
  fdsfdbzero1scaledmin->Write();
  fdsfdbzero1scaledmax->Write();
  fdsfdbplus1scaledmin->Write();
  fdsfdbplus1scaledmax->Write();
  fdsfdbs1scaledmin->Write();
  fdsfdbs1scaledmax->Write();
  fdsfdlambdab1scaledmin->Write();
  fdsfdlambdab1scaledmax->Write();

  fdspr1func->Write();
  fdsfdbzero1func->Write();
  fdsfdbplus1func->Write();
  fdsfdbs1func->Write();
  fdsfdlambdab1func->Write();
  fdspr1norm->Write();
  fdspr1normmin->Write();
  fdspr1normmax->Write();
  fdsfdbzero1norm->Write();
  fdsfdbzero1normmin->Write();
  fdsfdbzero1normmax->Write();
  fdsfdbplus1norm->Write();
  fdsfdbplus1normmin->Write();
  fdsfdbplus1normmax->Write();
  fdsfdbs1norm->Write();
  fdsfdbs1normmin->Write();
  fdsfdbs1normmax->Write();
  fdsfdlambdab1norm->Write();
  fdsfdlambdab1normmin->Write();
  fdsfdlambdab1normmax->Write();

  if(!drawall) cpol1->Close();
  if(!drawall) cpol1scaled->Close();
}

void calculate_Signal_CombBackground_MLdep(Int_t ptbin, Int_t iscan, Float_t expBssig, Float_t &expBsbkg, Float_t &expBssgnf){

  TH1D* hsig = (TH1D*)fmass->Get(Form("hmass_sigpt_cand%d_%d_%d",ptbins[ptbin],ptbins[ptbin+1],iscan));
  TH1D* hbkg = (TH1D*)fmass->Get(Form("hmass_bkgpt_cand%d_%d_%d",ptbins[ptbin],ptbins[ptbin+1],iscan));

  TH1F* hBkgCorr = (TH1F*)fBkgCorr->Get("hCorrFacBs");

  //FIT SIGNAL AS INITIALISATION
  TCanvas* csignal = new TCanvas(Form("csignal_%d_%d",ptbin, iscan), Form("csignal_%d_%d",ptbin, iscan), 450, 400);
  csignal->cd();

  cout << "Using signal histogram = " << hsig->GetName() << endl;
  hsig->Rebin(rebin[ptbin]);
  hsig->Scale(expBssig / hsig->GetEntries());

  TF1* fsig = new TF1(Form("fsig_%d_%d",ptbin,iscan), "gaus",  5.32, 5.42);
  TF1* fsigfull = new TF1(Form("fsigfull_%d_%d",ptbin,iscan), "gaus", 5.07, 5.65);
  hsig->Fit(Form("fsig_%d_%d",ptbin,iscan), "R");
  for(int ip = 0; ip < fsig->GetNpar(); ip++) fsigfull->SetParameter(ip, fsig->GetParameter(ip));

  hsig->SetTitle(Form("%d < #it{p}_{T} < %d GeV/#it{c};#it{M}(KK#pi#pi) (GeV/#it{c}^{2});Counts per %.3f MeV/#it{c}^{2} (norm.)",ptbins[ptbin], ptbins[ptbin+1], hsig->GetBinWidth(1)));
  hsig->SetName(Form("signal_%d_%d", ptbin, iscan));
  hsig->SetLineColor(kBlack);
  SetStyleHisto(hsig);
  hsig->SetStats(1);

  gStyle->SetOptFit(1111);
  hsig->Draw("hist");
  gPad->SetTickx();
  gPad->SetTicky();
  fsigfull->Draw("same");

  Double_t xmin = fsig->GetParameter(1) - 3 * fsig->GetParameter(2);
  Double_t xmax = fsig->GetParameter(1) + 3 * fsig->GetParameter(2);
  if(!drawall) csignal->Close();

  //EXTRACT COMBINATORIAL BACKGROUND
  TCanvas* ccombbkg = new TCanvas(Form("ccombbkg_%d_%d",ptbin, iscan), Form("ccombbkg_%d_%d",ptbin, iscan), 450, 400);
  ccombbkg->cd();

  cout << "Using background histogram = " << hbkg->GetName() << endl;
  hbkg->Rebin(rebin[ptbin]);

  hbkg->SetTitle(Form("%d < #it{p}_{T} < %d GeV/#it{c};#it{M}(KK#pi#pi) (GeV/#it{c}^{2});Counts per %.2f MeV/#it{c}^{2}",ptbins[ptbin], ptbins[ptbin+1], hbkg->GetBinWidth(1)));
  hbkg->SetName(Form("combbkg_%d_%d", ptbin, iscan));
  hbkg->SetLineColor(kBlack);
  hbkg->SetMarkerColor(kBlack);
  hbkg->SetMarkerStyle(20);
  SetStyleHisto(hbkg);

  //=Chi2 fit, add "L" for likelyhood
  TF1* fbkg = new TF1(Form("f_%d_%d",ptbin,iscan), bkgfitoption[ptbin], 5.07, 5.65);
  hbkg->Fit(Form("f_%d_%d",ptbin,iscan), "R,E,+");

  TGraph* grpar0cent = (TGraph*)ffits->Get(Form("gpar0cent_%d",ptbin));
  TGraph* grpar0low = (TGraph*)ffits->Get(Form("gpar0low_%d",ptbin));
  TGraph* grpar0up = (TGraph*)ffits->Get(Form("gpar0high_%d",ptbin));
  TF1* fpar1cent = (TF1*)ffits->Get(Form("fpar1cent_%d",ptbin));
  TF1* fpar1low = (TF1*)ffits->Get(Form("fpar1low_%d",ptbin));
  TF1* fpar1up = (TF1*)ffits->Get(Form("fpar1high_%d",ptbin));
  TGraph* grpar1cent = (TGraph*)ffits->Get(Form("gpar1cent_%d",ptbin));
  TGraph* grpar1low = (TGraph*)ffits->Get(Form("gpar1low_%d",ptbin));
  TGraph* grpar1up = (TGraph*)ffits->Get(Form("gpar1high_%d",ptbin));

  TF1* fcent = (TF1*)fbkg->Clone(Form("fcent_%d_%d",ptbin,iscan));
  fcent->SetParameter(0,grpar0cent->Eval(probcuts[iscan]));
  double dpar1cent;
  if(grpar1cent) dpar1cent = grpar1cent->Eval(probcuts[iscan]);
  else           dpar1cent = fpar1cent->Eval(probcuts[iscan]);
  fcent->SetParameter(1,dpar1cent);
  fcent->SetLineColor(kBlue);

  TF1* flow = (TF1*)fbkg->Clone(Form("flow_%d_%d",ptbin,iscan));
  flow->SetParameter(0,grpar0low->Eval(probcuts[iscan]));
  double dpar1low;
  if(grpar1low) dpar1low = grpar1low->Eval(probcuts[iscan]);
  else          dpar1low = fpar1low->Eval(probcuts[iscan]);
  flow->SetParameter(1,dpar1low);
  flow->SetLineStyle(2);
  flow->SetLineColor(kBlue);

  TF1* fup = (TF1*)fbkg->Clone(Form("fup_%d_%d",ptbin,iscan));
  fup->SetParameter(0,grpar0up->Eval(probcuts[iscan]));
  double dpar1up;
  if(grpar1up) dpar1up = grpar1up->Eval(probcuts[iscan]);
  else         dpar1up = fpar1up->Eval(probcuts[iscan]);
  fup->SetParameter(1,dpar1up);
  fup->SetLineStyle(2);
  fup->SetLineColor(kBlue);

  hbkg->Draw("ep");
  gPad->SetTickx();
  gPad->SetTicky();
  fbkg->Draw("same");
  fcent->Draw("same");
  fup->Draw("same");
  flow->Draw("same");

  Double_t maxplotax = hbkg->GetMaximum();
  Double_t maxplotaxnew = 1.3 * maxplotax;
  if(maxplotax < 10) maxplotaxnew = 15;
  if(maxplotax > 10 &&  maxplotax < 15) maxplotaxnew = 1.5 * maxplotax;
  hbkg->GetYaxis()->SetRangeUser(0., maxplotaxnew);

  Double_t maxplot = hbkg->GetBinContent(hbkg->FindBin(5.366));
  Double_t xgr[1] = {xmin + 0.5*(xmax - xmin)};
  Double_t ygr[1] = {0.2 * maxplot};
  Double_t xerrgr[1] = {0.5*(xmax - xmin)};
  Double_t yerrgr[1] = {0.2 * maxplot};
  TGraphErrors* grrange = new TGraphErrors(1, xgr, ygr, xerrgr, yerrgr);
  grrange->SetFillColor(kRed); grrange->SetLineWidth(2); grrange->SetLineColor(kRed); grrange->SetFillStyle(3445);
  grrange->Draw("same2");

  hbkg->Draw("same ep");

  TLegend* leg = new TLegend(0.4, 0.60, 0.92, 0.84, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(16); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->SetHeader(Form("Combinatorial;  %s event",nEv335String.Data()));
  leg->AddEntry(hbkg, "HIJING", "lpm");
  leg->AddEntry(fbkg, Form("%s fit",bkgfitoption[ptbin].Data()), "l");
  leg->AddEntry(fcent, Form("%s fit (parametrised)",bkgfitoption[ptbin].Data()), "l");
  leg->AddEntry(grrange, "3#sigma signal region", "f");
  leg->Draw();
  /*if(!drawall) ccombbkg->Close();*/

  TCanvas* cfits = (TCanvas*)ffits->Get(Form("BkgFits_parametrised_%d",ptbin));
  cfits->Draw();
  if(!drawall) cfits->Close();


  //EXTRACT SCALED COMBINATORIAL BACKGROUND
  TCanvas* ccombbkgscaled = new TCanvas(Form("ccombbkgscaled_%d_%d",ptbin, iscan), Form("ccombbkgscaled_%d_%d",ptbin, iscan), 450, 400);
  ccombbkgscaled->cd();

  Double_t bkgscalingHIJING = hBkgCorr->GetBinContent(hBkgCorr->FindBin(ptbins[ptbin]+0.01)); //HIJING correction
  Double_t bkgscalingEv = nEvExpected / nEv335; //Number of event correction
  cout << "\n Background HIJING scaling [" << ptbins[ptbin] << "," << ptbins[ptbin+1] << "]: " << bkgscalingHIJING << endl;
  cout << " Background HIJING scaling [" << ptbins[ptbin] << "," << ptbins[ptbin+1] << "]: " << bkgscalingEv << endl;

  TF1* fbkgnorm = new TF1(Form("fbkgnorm%d_%d",ptbin,iscan),"[0]",fbkg->GetXmin(), fbkg->GetXmax());
  fbkgnorm->SetParameter(0, bkgscalingHIJING * bkgscalingEv);
  cout << "   Setting par0 for fbkgnorm = " << bkgscalingHIJING * bkgscalingEv << endl;

  TF1* fbkgscaled;
  if(bkgfitoption[ptbin] == "pol1"){
    fbkgscaled = new TF1(Form("fbkgscaled%d_%d",ptbin,iscan),"pol1", fbkg->GetXmin(), fbkg->GetXmax());
    fbkgscaled->SetParameter(0, fbkgnorm->GetParameter(0) * fbkg->GetParameter(0));
    fbkgscaled->SetParameter(1, fbkgnorm->GetParameter(0) * fbkg->GetParameter(1));
  } else fbkgscaled = new TF1(Form("fbkgscaled%d_%d",ptbin,iscan),Form("%s*%s",fbkgnorm->GetName(), fbkg->GetName()), fbkg->GetXmin(), fbkg->GetXmax());
  fbkgscaled->SetLineColor(kMagenta+2);
  fbkgscaled->SetLineWidth(fbkg->GetLineWidth());
  fbkgscaled->SetLineStyle(fbkg->GetLineStyle());

  TF1* fcentscaled;
  if(bkgfitoption[ptbin] == "pol1"){
    fcentscaled = new TF1(Form("fcentscaled%d_%d",ptbin,iscan),"pol1", fcent->GetXmin(), fcent->GetXmax());
    fcentscaled->SetParameter(0, fbkgnorm->GetParameter(0) * fcent->GetParameter(0));
    fcentscaled->SetParameter(1, fbkgnorm->GetParameter(0) * fcent->GetParameter(1));
  } else fcentscaled = new TF1(Form("fcentscaled%d_%d",ptbin,iscan),Form("%s*%s",fbkgnorm->GetName(), fcent->GetName()), fcent->GetXmin(), fcent->GetXmax());
  fcentscaled->SetLineColor(kMagenta+2);
  fcentscaled->SetLineWidth(fcent->GetLineWidth());
  fcentscaled->SetLineStyle(fcent->GetLineStyle());
  TF1* fupscaled;
  if(bkgfitoption[ptbin] == "pol1"){
    fupscaled = new TF1(Form("fupscaled%d_%d",ptbin,iscan),"pol1", fup->GetXmin(), fup->GetXmax());
    fupscaled->SetParameter(0, fbkgnorm->GetParameter(0) * fup->GetParameter(0));
    fupscaled->SetParameter(1, fbkgnorm->GetParameter(0) * fup->GetParameter(1));
  } else fupscaled = new TF1(Form("fupscaled%d_%d",ptbin,iscan),Form("%s*%s",fbkgnorm->GetName(), fup->GetName()), fup->GetXmin(), fup->GetXmax());
  fupscaled->SetLineColor(kMagenta+2);
  fupscaled->SetLineWidth(fup->GetLineWidth());
  fupscaled->SetLineStyle(fup->GetLineStyle());
  TF1* flowscaled;
  if(bkgfitoption[ptbin] == "pol1"){
    flowscaled = new TF1(Form("flowscaled%d_%d",ptbin,iscan),"pol1", flow->GetXmin(), flow->GetXmax());
    flowscaled->SetParameter(0, fbkgnorm->GetParameter(0) * flow->GetParameter(0));
    flowscaled->SetParameter(1, fbkgnorm->GetParameter(0) * flow->GetParameter(1));
  } else flowscaled = new TF1(Form("flowscaled%d_%d",ptbin,iscan),Form("%s*%s",fbkgnorm->GetName(), flow->GetName()), flow->GetXmin(), flow->GetXmax());
  flowscaled->SetLineColor(kMagenta+2);
  flowscaled->SetLineWidth(flow->GetLineWidth());
  flowscaled->SetLineStyle(flow->GetLineStyle());

  //+1 because injected Ds + HIJING pion starts with ML=0, where HIJING bkg doesn't have this
  TF1* fdspr1scaled = (TF1*)fExpDsSignalScaled->Get(Form("fdspr1scaled%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdbzero1scaled = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdbzero1scaled%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdbplus1scaled = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdbplus1scaled%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdbs1scaled = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdbs1scaled%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdlambdab1scaled = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdlambdab1scaled%d_%d",ptbin,iscan/10+1));

  TF1* fdspr1scaledmin = (TF1*)fExpDsSignalScaled->Get(Form("fdspr1scaledmin%d_%d",ptbin,iscan/10+1));
  TF1* fdspr1scaledmax = (TF1*)fExpDsSignalScaled->Get(Form("fdspr1scaledmax%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdbzero1scaledmin = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdbzero1scaledmin%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdbzero1scaledmax = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdbzero1scaledmax%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdbplus1scaledmin = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdbplus1scaledmin%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdbplus1scaledmax = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdbplus1scaledmax%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdbs1scaledmin = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdbs1scaledmin%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdbs1scaledmax = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdbs1scaledmax%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdlambdab1scaledmin = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdlambdab1scaledmin%d_%d",ptbin,iscan/10+1));
  TF1* fdsfdlambdab1scaledmax = (TF1*)fExpDsSignalScaled->Get(Form("fdsfdlambdab1scaledmax%d_%d",ptbin,iscan/10+1));

  TH1D *href = (TH1D*)hbkg->Clone(Form("href%d_%d",ptbin,iscan));
  href->Reset("ICEMS");
  href->SetMinimum(0.01);
  href->SetMaximum(10e7);
  SetStyleHisto(href);

  href->Draw();
  gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  if(useRealFit[ptbin]) fbkgscaled->Draw("same");
  else fcentscaled->Draw("same");
  //fupscaled->Draw("same");
  //flowscaled->Draw("same");

  fdspr1scaled->Draw("same");
  fdspr1scaledmax->Draw("same");
  fdspr1scaledmin->Draw("same");
  fdsfdbzero1scaled->Draw("same");
  fdsfdbzero1scaledmax->Draw("same");
  fdsfdbzero1scaledmin->Draw("same");
  fdsfdbplus1scaled->Draw("same");
  fdsfdbplus1scaledmax->Draw("same");
  fdsfdbplus1scaledmin->Draw("same");
  fdsfdbs1scaled->Draw("same");
  fdsfdbs1scaledmax->Draw("same");
  fdsfdbs1scaledmin->Draw("same");
  fdsfdlambdab1scaled->Draw("same");
  fdsfdlambdab1scaledmax->Draw("same");
  fdsfdlambdab1scaledmin->Draw("same");
  fsigfull->Draw("same");

  TLegend* leg2 = new TLegend(0.4, 0.50, 0.92, 0.84, 0, "NDC");
  leg2->SetTextFont(43); leg2->SetTextSize(16); leg2->SetFillColor(0); leg2->SetFillStyle(0); leg2->SetLineColor(0);
  leg2->SetHeader("Background;  10 nb^{-1}");
  leg2->AddEntry(fcentscaled, "Combinatorial", "l");
  leg2->AddEntry(fdspr1scaled, "(Prompt)   #scale[0.6]{ }D_{s} + #pi", "l");
  leg2->AddEntry(fdsfdbzero1scaled, "(B^{0} #rightarrow X +) D_{s} + #pi", "l");
  leg2->AddEntry(fdsfdbplus1scaled, "(B^{+} #rightarrow X +) D_{s} + #pi", "l");
  leg2->AddEntry(fdsfdbs1scaled, "(B_{s} #rightarrow X +) D_{s} + #pi", "l");
  leg2->AddEntry(fdsfdlambdab1scaled, "(#Lambda_{b} #rightarrow X +) D_{s} + #pi", "l");
  leg2->Draw();

  TLatex info; info.SetNDC(); info.SetTextFont(43); info.SetTextSize(23);
  info.DrawLatex(0.19, 0.14, Form("ML_{1} > %.3f + ML_{2} > %.3f", preselML[ptbin], probcuts[iscan]));
  if(!drawall) ccombbkgscaled->Close();

  Double_t bkgcentfitcomb = 0;
  if(useRealFit[ptbin]) bkgcentfitcomb += fbkgscaled->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  else                  bkgcentfitcomb += fcentscaled->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  cout << "    Bkg combinatorial: " << bkgcentfitcomb << endl;

  Double_t bkgcentfonll = 0;
  bkgcentfonll += fdspr1scaled->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgcentfonll += fdsfdbs1scaled->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgcentfonll += fdsfdbzero1scaled->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgcentfonll += fdsfdbplus1scaled->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgcentfonll += fdsfdlambdab1scaled->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  expBsbkg = bkgcentfonll + bkgcentfitcomb;
  cout << "    Bkg combinatorial+injected: " << expBsbkg << " " << expBsbkg/bkgcentfitcomb << endl;

  Double_t bkgfitmax = 0;
  if(useRealFit[ptbin]) bkgfitmax += 0;
  else                  bkgfitmax += fupscaled->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfitmax += bkgcentfonll;

  Double_t bkgfonllmax = 0;
  bkgfonllmax += fdspr1scaledmax->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmax += fdsfdbs1scaledmax->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmax += fdsfdbzero1scaledmax->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmax += fdsfdbplus1scaledmax->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmax += fdsfdlambdab1scaledmax->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmax += bkgcentfitcomb;

  Double_t bkgfitmin = 0;
  if(useRealFit[ptbin]) bkgfitmin += 0;
  else                  bkgfitmin += flowscaled->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfitmin += bkgcentfonll;

  Double_t bkgfonllmin = 0;
  bkgfonllmin += fdspr1scaledmin->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmin += fdsfdbs1scaledmin->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmin += fdsfdbzero1scaledmin->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmin += fdsfdbplus1scaledmin->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmin += fdsfdlambdab1scaledmin->Integral(xmin,xmax)/(Double_t)hbkg->GetBinWidth(1);
  bkgfonllmin += bkgcentfitcomb;

  Double_t ptmin = ptbins[ptbin];
  Double_t ptmax = ptbins[ptbin+1];
  Double_t ptcent = 0.5*(ptmax + ptmin);
  Double_t pterr = 0.5*(ptmax - ptmin);

  grBkgHIJING->SetPoint(ptbin, ptcent, expBsbkg);
  grBkgHIJING->SetPointError(ptbin, pterr, pterr, expBsbkg-bkgfitmin, bkgfitmax-expBsbkg);
  grBkgFONLL->SetPoint(ptbin, ptcent, expBsbkg);
  grBkgFONLL->SetPointError(ptbin, pterr, pterr, expBsbkg-bkgfonllmin, bkgfonllmax-expBsbkg);

  grBkgHIJINGrel->SetPoint(ptbin, ptcent, 0);
  grBkgHIJINGrel->SetPointError(ptbin, pterr, pterr, (expBsbkg-bkgfitmin)/expBsbkg, (bkgfitmax-expBsbkg)/expBsbkg);
  grBkgFONLLrel->SetPoint(ptbin, ptcent, 0);
  grBkgFONLLrel->SetPointError(ptbin, pterr, pterr, (expBsbkg-bkgfonllmin)/expBsbkg, (bkgfonllmax-expBsbkg)/expBsbkg);

  expBssgnf = expBssig / TMath::Sqrt(expBssig + expBsbkg);
  cout << "     Sngf = " << expBssgnf << " (old = " << expBssig / TMath::Sqrt(expBssig + bkgcentfitcomb) << ")\n\n" << endl;

  //ccombbkg->cd();
  //TLatex info2; info2.SetNDC(); info2.SetTextFont(43); info2.SetTextSize(23);
  //info2.DrawLatex(0.19, 0.14, Form("Sgnf_comb = %.3f, Sgnf_all = %.3f", expBssgnf, expBssig / TMath::Sqrt(expBssig + bkgcentfitcomb)));
  //ccombbkg->SaveAs(Form("tempfig_ITS2/%s.eps",ccombbkg->GetName()));

}

void plot_expected_significance(TString filenamemass, TString filenameeff, TString filenamefits){
  
  TGaxis::SetMaxDigits(3);
  
  //TFile* fmass = new TFile(filenamemass.Data());
  //TFile* feff = new TFile(filenameeff.Data());
  //TFile* fBkgCorr = new TFile(filenameBkgCorr.Data());

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
       i != selbin[4] ) continue;
      
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

expectedbkg[j] = grBkgHIJING->GetY()[j];
expectedbkglow[j] = grBkgHIJING->GetY()[j] - grBkgHIJING->GetEYlow()[j];
expectedbkghigh[j] = grBkgHIJING->GetY()[j] + grBkgHIJING->GetEYhigh()[j];

  //if(j == 0) expectedbkghigh[j] = grBkgHIJING->GetY()[j] + 3*grBkgHIJING->GetEYlow()[j];
  //if(j == 1) expectedbkghigh[j] = grBkgHIJING->GetY()[j] + 2*grBkgHIJING->GetEYlow()[j];
  //if(j == 2) expectedbkghigh[j] = grBkgHIJING->GetY()[j] + grBkgHIJING->GetEYlow()[j];
  //grBkgHIJING2->SetPointEYhigh(0, 3*grBkgHIJING2->GetErrorYlow(0));
  //grBkgHIJING2->SetPointEYhigh(1, 2*grBkgHIJING2->GetErrorYlow(1));
  //grBkgHIJING2->SetPointEYhigh(2, grBkgHIJING2->GetErrorYlow(2));

  if(j == 0) expectedbkghigh[j] = grBkgHIJING->GetY()[j] + 2*grBkgHIJING->GetEYlow()[j];
  if(j == 1) expectedbkghigh[j] = grBkgHIJING->GetY()[j] + grBkgHIJING->GetEYlow()[j];
  if(j == 2) expectedbkglow[j] = grBkgHIJING->GetY()[j] - 0.75 * grBkgHIJING->GetEYhigh()[j];
  if(j == 2) expectedbkghigh[j] = grBkgHIJING->GetY()[j] + 0.75 * grBkgHIJING->GetEYhigh()[j];\
  if(j == 3) expectedbkghigh[j] = grBkgHIJING->GetY()[j] + grBkgHIJING->GetEYlow()[j];
  //grBkgHIJING3->SetPointEYhigh(0, 2*grBkgHIJING3->GetErrorYlow(0));
  //grBkgHIJING3->SetPointEYhigh(1, grBkgHIJING3->GetErrorYlow(1));
  //grBkgHIJING3->SetPointEYlow(2, 0.75*grBkgHIJING3->GetErrorYhigh(2));
  //grBkgHIJING3->SetPointEYhigh(2, 0.75*grBkgHIJING3->GetErrorYhigh(2));
  //grBkgHIJING3->SetPointEYhigh(3, grBkgHIJING3->GetErrorYlow(3));

      expectedsignalcent[j] *= gauss3sigmafactor;
      expectedsignalfonllmin[j] *= gauss3sigmafactor;
      expectedsignalfonllmax[j] *= gauss3sigmafactor;
      expectedsignalerreffmin[j] *= gauss3sigmafactor;
      expectedsignalerreffmax[j] *= gauss3sigmafactor;

      expectedsgnfcent[j] = expectedsignalcent[j] / TMath::Sqrt(expectedsignalcent[j] + expectedbkg[j]);

      Double_t significance = expectedsgnfcent[j];
      Double_t signal = expectedsignalcent[j];
      Double_t background = expectedbkg[j];
      Double_t sigPlusBkg = signal + background;

      Double_t errBkgSq = 0;
      Double_t errSigSq = (expectedsignalcent[j] - expectedsignalfonllmin[j])*(expectedsignalcent[j] - expectedsignalfonllmin[j]);
      expectedsgnffonllmin[j] = significance - significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
      errSigSq = (expectedsignalfonllmax[j] - expectedsignalcent[j])*(expectedsignalfonllmax[j] - expectedsignalcent[j]);
      expectedsgnffonllmax[j] = significance + significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

      errSigSq = (expectedsignalcent[j] - expectedsignalerreffmin[j])*(expectedsignalcent[j] - expectedsignalerreffmin[j]);
      expectedsgnferreffmin[j] = significance - significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
      errSigSq = (expectedsignalerreffmax[j] - expectedsignalcent[j])*(expectedsignalerreffmax[j] - expectedsignalcent[j]);
      expectedsgnferreffmax[j] = significance + significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

      errSigSq = 0;
      errBkgSq = (expectedbkghigh[j] - expectedbkg[j])*(expectedbkghigh[j] - expectedbkg[j]);
      expectedsgnfbkgmin[j] = significance - significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
      errBkgSq = (expectedbkg[j] - expectedbkglow[j])*(expectedbkg[j] - expectedbkglow[j]);
      expectedsgnfbkgmax[j] = significance + significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

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
      canv[j]->SaveAs(Form("tempfig/canv%d_%d_%d.eps",ptbins[j],ptbins[j+1],i));
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
    cout << filenamefits << endl;
    TGraph* grpar0cent = (TGraph*)ffits->Get(Form("gpar0cent_%d",j));
    TGraph* grpar0low = (TGraph*)ffits->Get(Form("gpar0low_%d",j));
    TGraph* grpar0up = (TGraph*)ffits->Get(Form("gpar0high_%d",j));
    //TF1* fpar1cent = (TF1*)ffits->Get(Form("fpar1cent_%d",j));
    //TF1* fpar1low = (TF1*)ffits->Get(Form("fpar1low_%d",j));
    //TF1* fpar1up = (TF1*)ffits->Get(Form("fpar1high_%d",j));
    TGraph* grpar1cent = (TGraph*)ffits->Get(Form("gpar1cent_%d",j));
    TGraph* grpar1low = (TGraph*)ffits->Get(Form("gpar1low_%d",j));
    TGraph* grpar1up = (TGraph*)ffits->Get(Form("gpar1high_%d",j));

    fcent = (TF1*)f1->Clone(Form("fcent_%d",j));
    fcent->SetParameter(0,grpar0cent->Eval(probcuts[i]));
    //double dpar1cent = fpar1cent->Eval(probcuts[i]);
    double dpar1cent = grpar1cent->Eval(probcuts[i]);
    fcent->SetParameter(1,dpar1cent);
    fcent->SetLineColor(kBlue);
    flow = (TF1*)f1->Clone(Form("flow_%d",j));
    flow->SetParameter(0,grpar0low->Eval(probcuts[i]));
    //double dpar1low = fpar1low->Eval(probcuts[i]);
    double dpar1low = grpar1low->Eval(probcuts[i]);
    flow->SetParameter(1,dpar1low);
    flow->SetLineStyle(2);
    flow->SetLineColor(kBlue);
    fup = (TF1*)f1->Clone(Form("fup_%d",j));
    fup->SetParameter(0,grpar0up->Eval(probcuts[i]));
    //double dpar1up = fpar1up->Eval(probcuts[i]);
    double dpar1up = grpar1up->Eval(probcuts[i]);
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
    for(Int_t i=0; i<10;i++) fscanf(infil,"%f",&dum);
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

void SetStyleHisto(TH1D *h){

  h->SetStats(0);
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetLabelSize(0.04);
  //h->GetYaxis()->SetDecimals(kTRUE);
  //h->GetYaxis()->SetNdivisions(507);
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetNdivisions(505);

}

void InitializeProbCuts(){

  if(opt == 2){
    Float_t probcuts_temp[trials+1] = {0.7, 0.7001, 0.7001999999999999, 0.7002999999999999, 0.7003999999999999, 0.7004999999999999, 0.7006, 0.7007, 0.7008, 0.7009, 0.701, 0.7011, 0.7011999999999999, 0.7012999999999999, 0.7013999999999999, 0.7014999999999999, 0.7016, 0.7017, 0.7018, 0.7019, 0.702, 0.7021, 0.7021999999999999, 0.7022999999999999, 0.7023999999999999, 0.7024999999999999, 0.7026, 0.7027, 0.7028, 0.7029, 0.703, 0.7031, 0.7031999999999999, 0.7032999999999999, 0.7033999999999999, 0.7034999999999999, 0.7036, 0.7037, 0.7038, 0.7039, 0.704, 0.7041, 0.7041999999999999, 0.7042999999999999, 0.7043999999999999, 0.7044999999999999, 0.7046, 0.7047, 0.7048, 0.7049, 0.705, 0.7051, 0.7051999999999999, 0.7052999999999999, 0.7053999999999999, 0.7054999999999999, 0.7056, 0.7057, 0.7058, 0.7059, 0.706, 0.7061, 0.7061999999999999, 0.7062999999999999, 0.7063999999999999, 0.7064999999999999, 0.7066, 0.7067, 0.7068, 0.7069, 0.707, 0.7071, 0.7071999999999999, 0.7072999999999999, 0.7073999999999999, 0.7074999999999999, 0.7076, 0.7077, 0.7078, 0.7079, 0.708, 0.7081, 0.7081999999999999, 0.7082999999999999, 0.7083999999999999, 0.7084999999999999, 0.7086, 0.7087, 0.7088, 0.7089, 0.709, 0.7091, 0.7091999999999999, 0.7092999999999999, 0.7093999999999999, 0.7094999999999999, 0.7096, 0.7097, 0.7098, 0.7099, 0.71, 0.7101, 0.7101999999999999, 0.7102999999999999, 0.7103999999999999, 0.7104999999999999, 0.7106, 0.7107, 0.7108, 0.7109, 0.711, 0.7111, 0.7111999999999999, 0.7112999999999999, 0.7113999999999999, 0.7114999999999999, 0.7116, 0.7117, 0.7118, 0.7119, 0.712, 0.7121, 0.7121999999999999, 0.7122999999999999, 0.7123999999999999, 0.7124999999999999, 0.7126, 0.7127, 0.7128, 0.7129, 0.713, 0.7131, 0.7132, 0.7132999999999999, 0.7133999999999999, 0.7134999999999999, 0.7136, 0.7137, 0.7138, 0.7139, 0.714, 0.7141, 0.7142, 0.7142999999999999, 0.7143999999999999, 0.7144999999999999, 0.7145999999999999, 0.7147, 0.7148, 0.7149, 0.715, 0.7151, 0.7152, 0.7152999999999999, 0.7153999999999999, 0.7154999999999999, 0.7156, 0.7157, 0.7158, 0.7159, 0.716, 0.7161, 0.7162, 0.7162999999999999, 0.7163999999999999, 0.7164999999999999, 0.7165999999999999, 0.7167, 0.7168, 0.7169, 0.717, 0.7171, 0.7172, 0.7172999999999999, 0.7173999999999999, 0.7174999999999999, 0.7176, 0.7177, 0.7178, 0.7179, 0.718, 0.7181, 0.7182, 0.7182999999999999, 0.7183999999999999, 0.7184999999999999, 0.7185999999999999, 0.7187, 0.7188, 0.7189, 0.719, 0.7191, 0.7192, 0.7192999999999999, 0.7193999999999999, 0.7194999999999999, 0.7195999999999999, 0.7197, 0.7198, 0.7199, 0.72, 0.7201, 0.7202, 0.7202999999999999, 0.7203999999999999, 0.7204999999999999, 0.7205999999999999, 0.7207, 0.7208, 0.7209, 0.721, 0.7211, 0.7212, 0.7212999999999999, 0.7213999999999999, 0.7214999999999999, 0.7215999999999999, 0.7217, 0.7218, 0.7219, 0.722, 0.7221, 0.7222, 0.7222999999999999, 0.7223999999999999, 0.7224999999999999, 0.7225999999999999, 0.7227, 0.7228, 0.7229, 0.723, 0.7231, 0.7232, 0.7232999999999999, 0.7233999999999999, 0.7234999999999999, 0.7235999999999999, 0.7237, 0.7238, 0.7239, 0.724, 0.7241, 0.7242, 0.7242999999999999, 0.7243999999999999, 0.7244999999999999, 0.7245999999999999, 0.7247, 0.7248, 0.7249, 0.725, 0.7251, 0.7252, 0.7253, 0.7253999999999999, 0.7254999999999999, 0.7255999999999999, 0.7257, 0.7258, 0.7259, 0.726, 0.7261, 0.7262, 0.7263, 0.7263999999999999, 0.7264999999999999, 0.7265999999999999, 0.7267, 0.7268, 0.7269, 0.727, 0.7271, 0.7272, 0.7273, 0.7273999999999999, 0.7274999999999999, 0.7275999999999999, 0.7277, 0.7278, 0.7279, 0.728, 0.7281, 0.7282, 0.7283, 0.7283999999999999, 0.7284999999999999, 0.7285999999999999, 0.7287, 0.7288, 0.7289, 0.729, 0.7291, 0.7292, 0.7293, 0.7293999999999999, 0.7294999999999999, 0.7295999999999999, 0.7297, 0.7298, 0.7299, 0.73, 0.7301, 0.7302, 0.7303, 0.7303999999999999, 0.7304999999999999, 0.7305999999999999, 0.7306999999999999, 0.7308, 0.7309, 0.731, 0.7311, 0.7312, 0.7313, 0.7313999999999999, 0.7314999999999999, 0.7315999999999999, 0.7317, 0.7318, 0.7319, 0.732, 0.7321, 0.7322, 0.7323, 0.7323999999999999, 0.7324999999999999, 0.7325999999999999, 0.7326999999999999, 0.7328, 0.7329, 0.733, 0.7331, 0.7332, 0.7333, 0.7333999999999999, 0.7334999999999999, 0.7335999999999999, 0.7337, 0.7338, 0.7339, 0.734, 0.7341, 0.7342, 0.7343, 0.7343999999999999, 0.7344999999999999, 0.7345999999999999, 0.7346999999999999, 0.7348, 0.7349, 0.735, 0.7351, 0.7352, 0.7353, 0.7353999999999999, 0.7354999999999999, 0.7355999999999999, 0.7357, 0.7358, 0.7359, 0.736, 0.7361, 0.7362, 0.7363, 0.7363999999999999, 0.7364999999999999, 0.7365999999999999, 0.7366999999999999, 0.7368, 0.7369, 0.737, 0.7371, 0.7372, 0.7373, 0.7373999999999999, 0.7374999999999999, 0.7375999999999999, 0.7376999999999999, 0.7378, 0.7379, 0.738, 0.7381, 0.7382, 0.7383, 0.7384, 0.7384999999999999, 0.7385999999999999, 0.7386999999999999, 0.7388, 0.7389, 0.739, 0.7391, 0.7392, 0.7393, 0.7394, 0.7394999999999999, 0.7395999999999999, 0.7396999999999999, 0.7398, 0.7399, 0.74, 0.7401, 0.7402, 0.7403, 0.7404, 0.7404999999999999, 0.7405999999999999, 0.7406999999999999, 0.7408, 0.7409, 0.741, 0.7411, 0.7412, 0.7413, 0.7414, 0.7414999999999999, 0.7415999999999999, 0.7416999999999999, 0.7418, 0.7419, 0.742, 0.7421, 0.7422, 0.7423, 0.7424, 0.7424999999999999, 0.7425999999999999, 0.7426999999999999, 0.7427999999999999, 0.7429, 0.743, 0.7431, 0.7432, 0.7433, 0.7434, 0.7434999999999999, 0.7435999999999999, 0.7436999999999999, 0.7438, 0.7439, 0.744, 0.7441, 0.7442, 0.7443, 0.7444, 0.7444999999999999, 0.7445999999999999, 0.7446999999999999, 0.7447999999999999, 0.7449, 0.745, 0.7451, 0.7452, 0.7453, 0.7454, 0.7454999999999999, 0.7455999999999999, 0.7456999999999999, 0.7458, 0.7459, 0.746, 0.7461, 0.7462, 0.7463, 0.7464, 0.7464999999999999, 0.7465999999999999, 0.7466999999999999, 0.7467999999999999, 0.7469, 0.747, 0.7471, 0.7472, 0.7473, 0.7474, 0.7474999999999999, 0.7475999999999999, 0.7476999999999999, 0.7478, 0.7479, 0.748, 0.7481, 0.7482, 0.7483, 0.7484, 0.7484999999999999, 0.7485999999999999, 0.7486999999999999, 0.7487999999999999, 0.7489, 0.749, 0.7491, 0.7492, 0.7493, 0.7494, 0.7494999999999999, 0.7495999999999999, 0.7496999999999999, 0.7498, 0.7499, 0.75, 0.7501, 0.7502, 0.7503, 0.7504, 0.7505, 0.7505999999999999, 0.7506999999999999, 0.7507999999999999, 0.7509, 0.751, 0.7511, 0.7512, 0.7513, 0.7514, 0.7515, 0.7515999999999999, 0.7516999999999999, 0.7518, 0.7519, 0.752, 0.7521, 0.7522, 0.7523, 0.7524, 0.7525, 0.7525999999999999, 0.7526999999999999, 0.7527999999999999, 0.7529, 0.753, 0.7531, 0.7532, 0.7533, 0.7534, 0.7535, 0.7535999999999999, 0.7536999999999999, 0.7537999999999999, 0.7539, 0.754, 0.7541, 0.7542, 0.7543, 0.7544, 0.7545, 0.7545999999999999, 0.7546999999999999, 0.7547999999999999, 0.7549, 0.755, 0.7551, 0.7552, 0.7553, 0.7554, 0.7555, 0.7555999999999999, 0.7556999999999999, 0.7557999999999999, 0.7559, 0.756, 0.7561, 0.7562, 0.7563, 0.7564, 0.7565, 0.7565999999999999, 0.7566999999999999, 0.7567999999999999, 0.7569, 0.757, 0.7571, 0.7572, 0.7573, 0.7574, 0.7575, 0.7575999999999999, 0.7576999999999999, 0.7577999999999999, 0.7579, 0.758, 0.7581, 0.7582, 0.7583, 0.7584, 0.7585, 0.7585999999999999, 0.7586999999999999, 0.7587999999999999, 0.7588999999999999, 0.759, 0.7591, 0.7592, 0.7593, 0.7594, 0.7595, 0.7595999999999999, 0.7596999999999999, 0.7597999999999999, 0.7599, 0.76, 0.7601, 0.7602, 0.7603, 0.7604, 0.7605, 0.7605999999999999, 0.7606999999999999, 0.7607999999999999, 0.7608999999999999, 0.761, 0.7611, 0.7612, 0.7613, 0.7614, 0.7615, 0.7615999999999999, 0.7616999999999999, 0.7617999999999999, 0.7619, 0.762, 0.7621, 0.7622, 0.7623, 0.7624, 0.7625, 0.7626, 0.7626999999999999, 0.7627999999999999, 0.7628999999999999, 0.763, 0.7631, 0.7632, 0.7633, 0.7634, 0.7635, 0.7636, 0.7636999999999999, 0.7637999999999999, 0.7639, 0.764, 0.7641, 0.7642, 0.7643, 0.7644, 0.7645, 0.7646, 0.7646999999999999, 0.7647999999999999, 0.7648999999999999, 0.765, 0.7651, 0.7652, 0.7653, 0.7654, 0.7655, 0.7656, 0.7656999999999999, 0.7657999999999999, 0.7659, 0.766, 0.7661, 0.7662, 0.7663, 0.7664, 0.7665, 0.7666, 0.7666999999999999, 0.7667999999999999, 0.7668999999999999, 0.767, 0.7671, 0.7672, 0.7673, 0.7674, 0.7675, 0.7676, 0.7676999999999999, 0.7677999999999999, 0.7679, 0.768, 0.7681, 0.7682, 0.7683, 0.7684, 0.7685, 0.7686, 0.7686999999999999, 0.7687999999999999, 0.7688999999999999, 0.7689999999999999, 0.7691, 0.7692, 0.7693, 0.7694, 0.7695, 0.7696, 0.7696999999999999, 0.7697999999999999, 0.7699, 0.77, 0.7701, 0.7702, 0.7703, 0.7704, 0.7705, 0.7706, 0.7706999999999999, 0.7707999999999999, 0.7708999999999999, 0.7709999999999999, 0.7711, 0.7712, 0.7713, 0.7714, 0.7715, 0.7716, 0.7716999999999999, 0.7717999999999999, 0.7719, 0.772, 0.7721, 0.7722, 0.7723, 0.7724, 0.7725, 0.7726, 0.7726999999999999, 0.7727999999999999, 0.7728999999999999, 0.7729999999999999, 0.7731, 0.7732, 0.7733, 0.7734, 0.7735, 0.7736, 0.7736999999999999, 0.7737999999999999, 0.7738999999999999, 0.774, 0.7741, 0.7742, 0.7743, 0.7744, 0.7745, 0.7746, 0.7746999999999999, 0.7747999999999999, 0.7748999999999999, 0.7749999999999999, 0.7751, 0.7752, 0.7753, 0.7754, 0.7755, 0.7756, 0.7757, 0.7757999999999999, 0.7758999999999999, 0.776, 0.7761, 0.7762, 0.7763, 0.7764, 0.7765, 0.7766, 0.7767, 0.7767999999999999, 0.7768999999999999, 0.7769999999999999, 0.7771, 0.7772, 0.7773, 0.7774, 0.7775, 0.7776, 0.7777, 0.7777999999999999, 0.7778999999999999, 0.778, 0.7781, 0.7782, 0.7783, 0.7784, 0.7785, 0.7786, 0.7787, 0.7787999999999999, 0.7788999999999999, 0.7789999999999999, 0.7791, 0.7792, 0.7793, 0.7794, 0.7795, 0.7796, 0.7797, 0.7797999999999999, 0.7798999999999999, 0.78, 0.7801, 0.7802, 0.7803, 0.7804, 0.7805, 0.7806, 0.7807, 0.7807999999999999, 0.7808999999999999, 0.7809999999999999, 0.7811, 0.7812, 0.7813, 0.7814, 0.7815, 0.7816, 0.7817, 0.7817999999999999, 0.7818999999999999, 0.782, 0.7821, 0.7822, 0.7823, 0.7824, 0.7825, 0.7826, 0.7827, 0.7827999999999999, 0.7828999999999999, 0.7829999999999999, 0.7831, 0.7832, 0.7833, 0.7834, 0.7835, 0.7836, 0.7837, 0.7837999999999999, 0.7838999999999999, 0.784, 0.7841, 0.7842, 0.7843, 0.7844, 0.7845, 0.7846, 0.7847, 0.7847999999999999, 0.7848999999999999, 0.7849999999999999, 0.7850999999999999, 0.7852, 0.7853, 0.7854, 0.7855, 0.7856, 0.7857, 0.7857999999999999, 0.7858999999999999, 0.786, 0.7861, 0.7862, 0.7863, 0.7864, 0.7865, 0.7866, 0.7867, 0.7867999999999999, 0.7868999999999999, 0.7869999999999999, 0.7870999999999999, 0.7872, 0.7873, 0.7874, 0.7875, 0.7876, 0.7877, 0.7878, 0.7878999999999999, 0.788, 0.7881, 0.7882, 0.7883, 0.7884, 0.7885, 0.7886, 0.7887, 0.7888, 0.7888999999999999, 0.7889999999999999, 0.7890999999999999, 0.7892, 0.7893, 0.7894, 0.7895, 0.7896, 0.7897, 0.7898, 0.7898999999999999, 0.7899999999999999, 0.7901, 0.7902, 0.7903, 0.7904, 0.7905, 0.7906, 0.7907, 0.7908, 0.7908999999999999, 0.7909999999999999, 0.7910999999999999, 0.7912, 0.7913, 0.7914, 0.7915, 0.7916, 0.7917, 0.7918, 0.7918999999999999, 0.7919999999999999, 0.7921, 0.7922, 0.7923, 0.7924, 0.7925, 0.7926, 0.7927, 0.7928, 0.7928999999999999, 0.7929999999999999, 0.7930999999999999, 0.7932, 0.7933, 0.7934, 0.7935, 0.7936, 0.7937, 0.7938, 0.7938999999999999, 0.7939999999999999, 0.7941, 0.7942, 0.7943, 0.7944, 0.7945, 0.7946, 0.7947, 0.7948, 0.7948999999999999, 0.7949999999999999, 0.7950999999999999, 0.7952, 0.7953, 0.7954, 0.7955, 0.7956, 0.7957, 0.7958, 0.7958999999999999, 0.7959999999999999, 0.7961, 0.7962, 0.7963, 0.7964, 0.7965, 0.7966, 0.7967, 0.7968, 0.7968999999999999, 0.7969999999999999, 0.7970999999999999, 0.7972, 0.7973, 0.7974, 0.7975, 0.7976, 0.7977, 0.7978, 0.7978999999999999, 0.7979999999999999, 0.7981, 0.7982, 0.7983, 0.7984, 0.7985, 0.7986, 0.7987, 0.7988, 0.7988999999999999, 0.7989999999999999, 0.7990999999999999, 0.7992, 0.7993, 0.7994, 0.7995, 0.7996, 0.7997, 0.7998, 0.7998999999999999, 0.7999999999999999, 0.8001, 0.8002, 0.8003, 0.8004, 0.8005, 0.8006, 0.8007, 0.8008, 0.8009, 0.8009999999999999, 0.8010999999999999, 0.8011999999999999, 0.8013, 0.8014, 0.8015, 0.8016, 0.8017, 0.8018, 0.8019, 0.8019999999999999, 0.8021, 0.8022, 0.8023, 0.8024, 0.8025, 0.8026, 0.8027, 0.8028, 0.8029, 0.8029999999999999, 0.8030999999999999, 0.8031999999999999, 0.8033, 0.8034, 0.8035, 0.8036, 0.8037, 0.8038, 0.8039, 0.8039999999999999, 0.8041, 0.8042, 0.8043, 0.8044, 0.8045, 0.8046, 0.8047, 0.8048, 0.8049, 0.8049999999999999, 0.8050999999999999, 0.8051999999999999, 0.8053, 0.8054, 0.8055, 0.8056, 0.8057, 0.8058, 0.8059, 0.8059999999999999, 0.8060999999999999, 0.8062, 0.8063, 0.8064, 0.8065, 0.8066, 0.8067, 0.8068, 0.8069, 0.8069999999999999, 0.8070999999999999, 0.8071999999999999, 0.8073, 0.8074, 0.8075, 0.8076, 0.8077, 0.8078, 0.8079, 0.8079999999999999, 0.8080999999999999, 0.8082, 0.8083, 0.8084, 0.8085, 0.8086, 0.8087, 0.8088, 0.8089, 0.8089999999999999, 0.8090999999999999, 0.8091999999999999, 0.8093, 0.8094, 0.8095, 0.8096, 0.8097, 0.8098, 0.8099, 0.8099999999999999, 0.8100999999999999, 0.8102, 0.8103, 0.8104, 0.8105, 0.8106, 0.8107, 0.8108, 0.8109, 0.8109999999999999, 0.8110999999999999, 0.8111999999999999, 0.8113, 0.8114, 0.8115, 0.8116, 0.8117, 0.8118, 0.8119, 0.8119999999999999, 0.8120999999999999, 0.8122, 0.8123, 0.8124, 0.8125, 0.8126, 0.8127, 0.8128, 0.8129, 0.813, 0.8130999999999999, 0.8131999999999999, 0.8133, 0.8134, 0.8135, 0.8136, 0.8137, 0.8138, 0.8139, 0.814, 0.8140999999999999, 0.8142, 0.8143, 0.8144, 0.8145, 0.8146, 0.8147, 0.8148, 0.8149, 0.815, 0.8150999999999999, 0.8151999999999999, 0.8153, 0.8154, 0.8155, 0.8156, 0.8157, 0.8158, 0.8159, 0.816, 0.8160999999999999, 0.8162, 0.8163, 0.8164, 0.8165, 0.8166, 0.8167, 0.8168, 0.8169, 0.817, 0.8170999999999999, 0.8171999999999999, 0.8172999999999999, 0.8174, 0.8175, 0.8176, 0.8177, 0.8178, 0.8179, 0.818, 0.8180999999999999, 0.8182, 0.8183, 0.8184, 0.8185, 0.8186, 0.8187, 0.8188, 0.8189, 0.819, 0.8190999999999999, 0.8191999999999999, 0.8192999999999999, 0.8194, 0.8195, 0.8196, 0.8197, 0.8198, 0.8199, 0.82, 0.8200999999999999, 0.8201999999999999, 0.8203, 0.8204, 0.8205, 0.8206, 0.8207, 0.8208, 0.8209, 0.821, 0.8210999999999999, 0.8211999999999999, 0.8212999999999999, 0.8214, 0.8215, 0.8216, 0.8217, 0.8218, 0.8219, 0.822, 0.8220999999999999, 0.8221999999999999, 0.8223, 0.8224, 0.8225, 0.8226, 0.8227, 0.8228, 0.8229, 0.823, 0.8230999999999999, 0.8231999999999999, 0.8232999999999999, 0.8234, 0.8235, 0.8236, 0.8237, 0.8238, 0.8239, 0.824, 0.8240999999999999, 0.8241999999999999, 0.8243, 0.8244, 0.8245, 0.8246, 0.8247, 0.8248, 0.8249, 0.825, 0.8251, 0.8251999999999999, 0.8252999999999999, 0.8253999999999999, 0.8255, 0.8256, 0.8257, 0.8258, 0.8259, 0.826, 0.8261, 0.8262, 0.8263, 0.8264, 0.8265, 0.8266, 0.8267, 0.8268, 0.8269, 0.827, 0.8271, 0.8271999999999999, 0.8272999999999999, 0.8273999999999999, 0.8275, 0.8276, 0.8277, 0.8278, 0.8279, 0.828, 0.8281, 0.8282, 0.8283, 0.8284, 0.8285, 0.8286, 0.8287, 0.8288, 0.8289, 0.829, 0.8291, 0.8291999999999999, 0.8292999999999999, 0.8293999999999999, 0.8295, 0.8296, 0.8297, 0.8298, 0.8299, 0.83, 0.8301, 0.8301999999999999, 0.8303, 0.8304, 0.8305, 0.8306, 0.8307, 0.8308, 0.8309, 0.831, 0.8311, 0.8311999999999999, 0.8312999999999999, 0.8313999999999999, 0.8315, 0.8316, 0.8317, 0.8318, 0.8319, 0.832, 0.8321, 0.8321999999999999, 0.8323, 0.8324, 0.8325, 0.8326, 0.8327, 0.8328, 0.8329, 0.833, 0.8331, 0.8331999999999999, 0.8332999999999999, 0.8333999999999999, 0.8335, 0.8336, 0.8337, 0.8338, 0.8339, 0.834, 0.8341, 0.8341999999999999, 0.8343, 0.8344, 0.8345, 0.8346, 0.8347, 0.8348, 0.8349, 0.835, 0.8351, 0.8351999999999999, 0.8352999999999999, 0.8353999999999999, 0.8355, 0.8356, 0.8357, 0.8358, 0.8359, 0.836, 0.8361, 0.8361999999999999, 0.8363, 0.8364, 0.8365, 0.8366, 0.8367, 0.8368, 0.8369, 0.837, 0.8371, 0.8371999999999999, 0.8372999999999999, 0.8373999999999999, 0.8375, 0.8376, 0.8377, 0.8378, 0.8379, 0.838, 0.8381, 0.8382, 0.8383, 0.8384, 0.8385, 0.8386, 0.8387, 0.8388, 0.8389, 0.839, 0.8391, 0.8392, 0.8392999999999999, 0.8393999999999999, 0.8394999999999999, 0.8396, 0.8397, 0.8398, 0.8399, 0.84, 0.8401, 0.8402, 0.8403, 0.8404, 0.8405, 0.8406, 0.8407, 0.8408, 0.8409, 0.841, 0.8411, 0.8412, 0.8412999999999999, 0.8413999999999999, 0.8414999999999999, 0.8416, 0.8417, 0.8418, 0.8419, 0.842, 0.8421, 0.8422, 0.8423, 0.8424, 0.8425, 0.8426, 0.8427, 0.8428, 0.8429, 0.843, 0.8431, 0.8432, 0.8432999999999999, 0.8433999999999999, 0.8434999999999999, 0.8436, 0.8437, 0.8438, 0.8439, 0.844, 0.8441, 0.8442, 0.8443, 0.8444, 0.8445, 0.8446, 0.8447, 0.8448, 0.8449, 0.845, 0.8451, 0.8452, 0.8452999999999999, 0.8453999999999999, 0.8454999999999999, 0.8456, 0.8457, 0.8458, 0.8459, 0.846, 0.8461, 0.8462, 0.8462999999999999, 0.8464, 0.8465, 0.8466, 0.8467, 0.8468, 0.8469, 0.847, 0.8471, 0.8472, 0.8472999999999999, 0.8473999999999999, 0.8474999999999999, 0.8476, 0.8477, 0.8478, 0.8479, 0.848, 0.8481, 0.8482, 0.8482999999999999, 0.8484, 0.8485, 0.8486, 0.8487, 0.8488, 0.8489, 0.849, 0.8491, 0.8492, 0.8492999999999999, 0.8493999999999999, 0.8494999999999999, 0.8496, 0.8497, 0.8498, 0.8499, 0.85, 0.8501, 0.8502, 0.8503, 0.8504, 0.8505, 0.8506, 0.8507, 0.8508, 0.8509, 0.851, 0.8511, 0.8512, 0.8513, 0.8513999999999999, 0.8514999999999999, 0.8516, 0.8517, 0.8518, 0.8519, 0.852, 0.8521, 0.8522, 0.8523, 0.8524, 0.8525, 0.8526, 0.8527, 0.8528, 0.8529, 0.853, 0.8531, 0.8532, 0.8533, 0.8533999999999999, 0.8534999999999999, 0.8535999999999999, 0.8537, 0.8538, 0.8539, 0.854, 0.8541, 0.8542, 0.8543, 0.8544, 0.8545, 0.8546, 0.8547, 0.8548, 0.8549, 0.855, 0.8551, 0.8552, 0.8553, 0.8553999999999999, 0.8554999999999999, 0.8555999999999999, 0.8557, 0.8558, 0.8559, 0.856, 0.8561, 0.8562, 0.8563, 0.8564, 0.8565, 0.8566, 0.8567, 0.8568, 0.8569, 0.857, 0.8571, 0.8572, 0.8573, 0.8573999999999999, 0.8574999999999999, 0.8575999999999999, 0.8577, 0.8578, 0.8579, 0.858, 0.8581, 0.8582, 0.8583, 0.8584, 0.8585, 0.8586, 0.8587, 0.8588, 0.8589, 0.859, 0.8591, 0.8592, 0.8593, 0.8593999999999999, 0.8594999999999999, 0.8595999999999999, 0.8597, 0.8598, 0.8599, 0.86, 0.8601, 0.8602, 0.8603, 0.8604, 0.8605, 0.8606, 0.8607, 0.8608, 0.8609, 0.861, 0.8611, 0.8612, 0.8613, 0.8613999999999999, 0.8614999999999999, 0.8615999999999999, 0.8617, 0.8618, 0.8619, 0.862, 0.8621, 0.8622, 0.8623, 0.8623999999999999, 0.8625, 0.8626, 0.8627, 0.8628, 0.8629, 0.863, 0.8631, 0.8632, 0.8633, 0.8634, 0.8634999999999999, 0.8635999999999999, 0.8637, 0.8638, 0.8639, 0.864, 0.8641, 0.8642, 0.8643, 0.8644, 0.8645, 0.8646, 0.8647, 0.8648, 0.8649, 0.865, 0.8651, 0.8652, 0.8653, 0.8654, 0.8654999999999999, 0.8655999999999999, 0.8657, 0.8658, 0.8659, 0.866, 0.8661, 0.8662, 0.8663, 0.8664, 0.8665, 0.8666, 0.8667, 0.8668, 0.8669, 0.867, 0.8671, 0.8672, 0.8673, 0.8674, 0.8674999999999999, 0.8675999999999999, 0.8677, 0.8678, 0.8679, 0.868, 0.8681, 0.8682, 0.8683, 0.8684, 0.8685, 0.8686, 0.8687, 0.8688, 0.8689, 0.869, 0.8691, 0.8692, 0.8693, 0.8694, 0.8694999999999999, 0.8695999999999999, 0.8696999999999999, 0.8698, 0.8699, 0.87, 0.8701, 0.8702, 0.8703, 0.8704, 0.8705, 0.8706, 0.8707, 0.8708, 0.8709, 0.871, 0.8711, 0.8712, 0.8713, 0.8714, 0.8714999999999999, 0.8715999999999999, 0.8716999999999999, 0.8718, 0.8719, 0.872, 0.8721, 0.8722, 0.8723, 0.8724, 0.8725, 0.8726, 0.8727, 0.8728, 0.8729, 0.873, 0.8731, 0.8732, 0.8733, 0.8734, 0.8734999999999999, 0.8735999999999999, 0.8736999999999999, 0.8738, 0.8739, 0.874, 0.8741, 0.8742, 0.8743, 0.8744, 0.8745, 0.8746, 0.8747, 0.8748, 0.8749, 0.875, 0.8751, 0.8752, 0.8753, 0.8754, 0.8755, 0.8755999999999999, 0.8756999999999999, 0.8758, 0.8759, 0.876, 0.8761, 0.8762, 0.8763, 0.8764, 0.8765000000000001, 0.8766, 0.8767, 0.8768, 0.8769, 0.877, 0.8771, 0.8772, 0.8773, 0.8774, 0.8775, 0.8775999999999999, 0.8776999999999999, 0.8778, 0.8779, 0.878, 0.8781, 0.8782, 0.8783, 0.8784, 0.8785, 0.8786, 0.8787, 0.8788, 0.8789, 0.879, 0.8791, 0.8792, 0.8793, 0.8794, 0.8795, 0.8795999999999999, 0.8796999999999999, 0.8798, 0.8799, 0.88, 0.8801, 0.8802, 0.8803, 0.8804, 0.8805, 0.8806, 0.8807, 0.8808, 0.8809, 0.881, 0.8811, 0.8812, 0.8813, 0.8814, 0.8815, 0.8815999999999999, 0.8816999999999999, 0.8818, 0.8819, 0.882, 0.8821, 0.8822, 0.8823, 0.8824, 0.8825, 0.8826, 0.8827, 0.8828, 0.8829, 0.883, 0.8831, 0.8832, 0.8833, 0.8834, 0.8835, 0.8835999999999999, 0.8836999999999999, 0.8838, 0.8839, 0.884, 0.8841, 0.8842, 0.8843, 0.8844, 0.8845, 0.8846, 0.8847, 0.8848, 0.8849, 0.885, 0.8851, 0.8852, 0.8853, 0.8854, 0.8855, 0.8855999999999999, 0.8856999999999999, 0.8857999999999999, 0.8859, 0.886, 0.8861, 0.8862, 0.8863, 0.8864, 0.8865, 0.8866, 0.8867, 0.8868, 0.8869, 0.887, 0.8871, 0.8872, 0.8873, 0.8874, 0.8875, 0.8876, 0.8876999999999999, 0.8877999999999999, 0.8879, 0.888, 0.8881, 0.8882, 0.8883, 0.8884, 0.8885, 0.8886000000000001, 0.8887, 0.8888, 0.8889, 0.889, 0.8891, 0.8892, 0.8893, 0.8894, 0.8895, 0.8896, 0.8896999999999999, 0.8897999999999999, 0.8899, 0.89, 0.8901, 0.8902, 0.8903, 0.8904, 0.8905, 0.8906000000000001, 0.8907, 0.8908, 0.8909, 0.891, 0.8911, 0.8912, 0.8913, 0.8914, 0.8915, 0.8916, 0.8916999999999999, 0.8917999999999999, 0.8919, 0.892, 0.8921, 0.8922, 0.8923, 0.8924, 0.8925, 0.8926000000000001, 0.8927, 0.8928, 0.8929, 0.893, 0.8931, 0.8932, 0.8933, 0.8934, 0.8935, 0.8936, 0.8936999999999999, 0.8937999999999999, 0.8939, 0.894, 0.8941, 0.8942, 0.8943, 0.8944, 0.8945, 0.8946, 0.8947, 0.8948, 0.8949, 0.895, 0.8951, 0.8952, 0.8953, 0.8954, 0.8955, 0.8956, 0.8956999999999999, 0.8957999999999999, 0.8959, 0.896, 0.8961, 0.8962, 0.8963, 0.8964, 0.8965, 0.8966, 0.8967, 0.8968, 0.8969, 0.897, 0.8971, 0.8972, 0.8973, 0.8974, 0.8975, 0.8976, 0.8976999999999999, 0.8977999999999999, 0.8979, 0.898, 0.8981, 0.8982, 0.8983, 0.8984, 0.8985, 0.8986, 0.8987, 0.8988, 0.8989, 0.899, 0.8991, 0.8992, 0.8993, 0.8994, 0.8995, 0.8996, 0.8996999999999999, 0.8997999999999999, 0.8999, 0.9, 0.9001, 0.9002, 0.9003, 0.9004, 0.9005, 0.9006, 0.9007000000000001, 0.9008, 0.9009, 0.901, 0.9011, 0.9012, 0.9013, 0.9014, 0.9015, 0.9016, 0.9017, 0.9017999999999999, 0.9018999999999999, 0.902, 0.9021, 0.9022, 0.9023, 0.9024, 0.9025, 0.9026, 0.9027000000000001, 0.9028, 0.9029, 0.903, 0.9031, 0.9032, 0.9033, 0.9034, 0.9035, 0.9036, 0.9037, 0.9037999999999999, 0.9038999999999999, 0.904, 0.9041, 0.9042, 0.9043, 0.9044, 0.9045, 0.9046, 0.9047000000000001, 0.9048, 0.9049, 0.905, 0.9051, 0.9052, 0.9053, 0.9054, 0.9055, 0.9056, 0.9057, 0.9057999999999999, 0.9058999999999999, 0.906, 0.9061, 0.9062, 0.9063, 0.9064, 0.9065, 0.9066, 0.9067000000000001, 0.9068, 0.9069, 0.907, 0.9071, 0.9072, 0.9073, 0.9074, 0.9075, 0.9076, 0.9077, 0.9077999999999999, 0.9078999999999999, 0.908, 0.9081, 0.9082, 0.9083, 0.9084, 0.9085, 0.9086, 0.9087000000000001, 0.9088, 0.9089, 0.909, 0.9091, 0.9092, 0.9093, 0.9094, 0.9095, 0.9096, 0.9097, 0.9097999999999999, 0.9098999999999999, 0.91, 0.9101, 0.9102, 0.9103, 0.9104, 0.9105, 0.9106, 0.9107, 0.9108, 0.9109, 0.911, 0.9111, 0.9112, 0.9113, 0.9114, 0.9115, 0.9116, 0.9117, 0.9117999999999999, 0.9118999999999999, 0.912, 0.9121, 0.9122, 0.9123, 0.9124, 0.9125, 0.9126, 0.9127, 0.9128000000000001, 0.9129, 0.913, 0.9131, 0.9132, 0.9133, 0.9134, 0.9135, 0.9136, 0.9137, 0.9138, 0.9138999999999999, 0.914, 0.9141, 0.9142, 0.9143, 0.9144, 0.9145, 0.9146, 0.9147, 0.9148000000000001, 0.9149, 0.915, 0.9151, 0.9152, 0.9153, 0.9154, 0.9155, 0.9156, 0.9157, 0.9158, 0.9158999999999999, 0.916, 0.9161, 0.9162, 0.9163, 0.9164, 0.9165, 0.9166, 0.9167, 0.9168000000000001, 0.9169, 0.917, 0.9171, 0.9172, 0.9173, 0.9174, 0.9175, 0.9176, 0.9177, 0.9178, 0.9178999999999999, 0.9179999999999999, 0.9181, 0.9182, 0.9183, 0.9184, 0.9185, 0.9186, 0.9187, 0.9188000000000001, 0.9189, 0.919, 0.9191, 0.9192, 0.9193, 0.9194, 0.9195, 0.9196, 0.9197, 0.9198, 0.9198999999999999, 0.9199999999999999, 0.9201, 0.9202, 0.9203, 0.9204, 0.9205, 0.9206, 0.9207, 0.9208000000000001, 0.9209, 0.921, 0.9211, 0.9212, 0.9213, 0.9214, 0.9215, 0.9216, 0.9217, 0.9218, 0.9218999999999999, 0.9219999999999999, 0.9221, 0.9222, 0.9223, 0.9224, 0.9225, 0.9226, 0.9227, 0.9228000000000001, 0.9229, 0.923, 0.9231, 0.9232, 0.9233, 0.9234, 0.9235, 0.9236, 0.9237, 0.9238, 0.9238999999999999, 0.9239999999999999, 0.9241, 0.9242, 0.9243, 0.9244, 0.9245, 0.9246, 0.9247, 0.9248, 0.9249, 0.925, 0.9251, 0.9252, 0.9253, 0.9254, 0.9255, 0.9256, 0.9257, 0.9258, 0.9259, 0.9259999999999999, 0.9261, 0.9262, 0.9263, 0.9264, 0.9265, 0.9266, 0.9267, 0.9268, 0.9269000000000001, 0.927, 0.9271, 0.9272, 0.9273, 0.9274, 0.9275, 0.9276, 0.9277, 0.9278, 0.9279, 0.9279999999999999, 0.9281, 0.9282, 0.9283, 0.9284, 0.9285, 0.9286, 0.9287, 0.9288, 0.9289000000000001, 0.929, 0.9291, 0.9292, 0.9293, 0.9294, 0.9295, 0.9296, 0.9297, 0.9298, 0.9299, 0.9299999999999999, 0.9301, 0.9302, 0.9303, 0.9304, 0.9305, 0.9306, 0.9307, 0.9308, 0.9309000000000001, 0.931, 0.9311, 0.9312, 0.9313, 0.9314, 0.9315, 0.9316, 0.9317, 0.9318, 0.9319, 0.9319999999999999, 0.9321, 0.9322, 0.9323, 0.9324, 0.9325, 0.9326, 0.9327, 0.9328, 0.9329000000000001, 0.933, 0.9331, 0.9332, 0.9333, 0.9334, 0.9335, 0.9336, 0.9337, 0.9338, 0.9339, 0.9339999999999999, 0.9340999999999999, 0.9342, 0.9343, 0.9344, 0.9345, 0.9346, 0.9347, 0.9348, 0.9349000000000001, 0.935, 0.9351, 0.9352, 0.9353, 0.9354, 0.9355, 0.9356, 0.9357, 0.9358, 0.9359, 0.9359999999999999, 0.9360999999999999, 0.9362, 0.9363, 0.9364, 0.9365, 0.9366, 0.9367, 0.9368, 0.9369000000000001, 0.937, 0.9371, 0.9372, 0.9373, 0.9374, 0.9375, 0.9376, 0.9377, 0.9378, 0.9379, 0.938, 0.9380999999999999, 0.9382, 0.9383, 0.9384, 0.9385, 0.9386, 0.9387, 0.9388, 0.9389000000000001, 0.9390000000000001, 0.9391, 0.9392, 0.9393, 0.9394, 0.9395, 0.9396, 0.9397, 0.9398, 0.9399, 0.94, 0.9400999999999999, 0.9402, 0.9403, 0.9404, 0.9405, 0.9406, 0.9407, 0.9408, 0.9409, 0.9410000000000001, 0.9411, 0.9412, 0.9413, 0.9414, 0.9415, 0.9416, 0.9417, 0.9418, 0.9419, 0.942, 0.9420999999999999, 0.9422, 0.9423, 0.9424, 0.9425, 0.9426, 0.9427, 0.9428, 0.9429, 0.9430000000000001, 0.9431, 0.9432, 0.9433, 0.9434, 0.9435, 0.9436, 0.9437, 0.9438, 0.9439, 0.944, 0.9440999999999999, 0.9442, 0.9443, 0.9444, 0.9445, 0.9446, 0.9447, 0.9448, 0.9449, 0.9450000000000001, 0.9451, 0.9452, 0.9453, 0.9454, 0.9455, 0.9456, 0.9457, 0.9458, 0.9459, 0.946, 0.9460999999999999, 0.9462, 0.9463, 0.9464, 0.9465, 0.9466, 0.9467, 0.9468, 0.9469, 0.9470000000000001, 0.9471, 0.9472, 0.9473, 0.9474, 0.9475, 0.9476, 0.9477, 0.9478, 0.9479, 0.948, 0.9480999999999999, 0.9481999999999999, 0.9483, 0.9484, 0.9485, 0.9486, 0.9487, 0.9488, 0.9489, 0.9490000000000001, 0.9491, 0.9492, 0.9493, 0.9494, 0.9495, 0.9496, 0.9497, 0.9498, 0.9499, 0.95, 0.9501, 0.9501999999999999, 0.9502999999999999, 0.9504, 0.9505, 0.9506, 0.9507, 0.9508, 0.9509000000000001, 0.9510000000000001, 0.9511000000000001, 0.9512, 0.9513, 0.9514, 0.9515, 0.9516, 0.9517, 0.9518, 0.9519, 0.952, 0.9521, 0.9521999999999999, 0.9522999999999999, 0.9524, 0.9525, 0.9526, 0.9527, 0.9528, 0.9529000000000001, 0.9530000000000001, 0.9531000000000001, 0.9532, 0.9533, 0.9534, 0.9535, 0.9536, 0.9537, 0.9538, 0.9539, 0.954, 0.9541, 0.9541999999999999, 0.9542999999999999, 0.9544, 0.9545, 0.9546, 0.9547, 0.9548, 0.9549000000000001, 0.9550000000000001, 0.9551000000000001, 0.9552, 0.9553, 0.9554, 0.9555, 0.9556, 0.9557, 0.9558, 0.9559, 0.956, 0.9561, 0.9561999999999999, 0.9562999999999999, 0.9564, 0.9565, 0.9566, 0.9567, 0.9568, 0.9569000000000001, 0.9570000000000001, 0.9571000000000001, 0.9572, 0.9573, 0.9574, 0.9575, 0.9576, 0.9577, 0.9578, 0.9579, 0.958, 0.9581, 0.9581999999999999, 0.9582999999999999, 0.9584, 0.9585, 0.9586, 0.9587, 0.9588, 0.9589000000000001, 0.9590000000000001, 0.9591000000000001, 0.9592, 0.9593, 0.9594, 0.9595, 0.9596, 0.9597, 0.9598, 0.9599, 0.96, 0.9601, 0.9601999999999999, 0.9602999999999999, 0.9603999999999999, 0.9605, 0.9606, 0.9607, 0.9608, 0.9609, 0.9610000000000001, 0.9611000000000001, 0.9612, 0.9613, 0.9614, 0.9615, 0.9616, 0.9617, 0.9618, 0.9619, 0.962, 0.9621, 0.9621999999999999, 0.9622999999999999, 0.9623999999999999, 0.9625, 0.9626, 0.9627, 0.9628, 0.9629, 0.9630000000000001, 0.9631000000000001, 0.9632000000000001, 0.9633, 0.9634, 0.9635, 0.9636, 0.9637, 0.9638, 0.9639, 0.964, 0.9641, 0.9642, 0.9642999999999999, 0.9643999999999999, 0.9645, 0.9646, 0.9647, 0.9648, 0.9649, 0.9650000000000001, 0.9651000000000001, 0.9652000000000001, 0.9653, 0.9654, 0.9655, 0.9656, 0.9657, 0.9658, 0.9659, 0.966, 0.9661, 0.9662, 0.9662999999999999, 0.9663999999999999, 0.9665, 0.9666, 0.9667, 0.9668, 0.9669, 0.9670000000000001, 0.9671000000000001, 0.9672000000000001, 0.9673, 0.9674, 0.9675, 0.9676, 0.9677, 0.9678, 0.9679, 0.968, 0.9681, 0.9682, 0.9682999999999999, 0.9683999999999999, 0.9685, 0.9686, 0.9687, 0.9688, 0.9689, 0.9690000000000001, 0.9691000000000001, 0.9692000000000001, 0.9693, 0.9694, 0.9695, 0.9696, 0.9697, 0.9698, 0.9699, 0.97, 0.9701, 0.9702, 0.9702999999999999, 0.9703999999999999, 0.9705, 0.9706, 0.9707, 0.9708, 0.9709, 0.9710000000000001, 0.9711000000000001, 0.9712000000000001, 0.9713, 0.9714, 0.9715, 0.9716, 0.9717, 0.9718, 0.9719, 0.972, 0.9721, 0.9722, 0.9722999999999999, 0.9723999999999999, 0.9725, 0.9726, 0.9727, 0.9728, 0.9729, 0.9730000000000001, 0.9731000000000001, 0.9732000000000001, 0.9733, 0.9734, 0.9735, 0.9736, 0.9737, 0.9738, 0.9739, 0.974, 0.9741, 0.9742, 0.9742999999999999, 0.9743999999999999, 0.9745, 0.9746, 0.9747, 0.9748, 0.9749, 0.9750000000000001, 0.9751000000000001, 0.9752000000000001, 0.9753000000000001, 0.9754, 0.9755, 0.9756, 0.9757, 0.9758, 0.9759, 0.976, 0.9761, 0.9762, 0.9763, 0.9763999999999999, 0.9764999999999999, 0.9766, 0.9767, 0.9768, 0.9769, 0.977, 0.9771000000000001, 0.9772000000000001, 0.9773000000000001, 0.9774, 0.9775, 0.9776, 0.9777, 0.9778, 0.9779, 0.978, 0.9781, 0.9782, 0.9783, 0.9783999999999999, 0.9784999999999999, 0.9786, 0.9787, 0.9788, 0.9789, 0.979, 0.9791000000000001, 0.9792000000000001, 0.9793000000000001, 0.9794, 0.9795, 0.9796, 0.9797, 0.9798, 0.9799, 0.98, 0.9801, 0.9802, 0.9803, 0.9803999999999999, 0.9804999999999999, 0.9806, 0.9807, 0.9808, 0.9809, 0.981, 0.9811000000000001, 0.9812000000000001, 0.9813000000000001, 0.9814, 0.9815, 0.9816, 0.9817, 0.9818, 0.9819, 0.982, 0.9821, 0.9822, 0.9823, 0.9823999999999999, 0.9824999999999999, 0.9826, 0.9827, 0.9828, 0.9829, 0.983, 0.9831000000000001, 0.9832000000000001, 0.9833000000000001, 0.9834, 0.9835, 0.9836, 0.9837, 0.9838, 0.9839, 0.984, 0.9841, 0.9842, 0.9843, 0.9843999999999999, 0.9844999999999999, 0.9846, 0.9847, 0.9848, 0.9849, 0.985, 0.9851000000000001, 0.9852000000000001, 0.9853000000000001, 0.9854, 0.9855, 0.9856, 0.9857, 0.9858, 0.9859, 0.986, 0.9861, 0.9862, 0.9863, 0.9863999999999999, 0.9864999999999999, 0.9866, 0.9867, 0.9868, 0.9869, 0.987, 0.9871000000000001, 0.9872000000000001, 0.9873000000000001, 0.9874, 0.9875, 0.9876, 0.9877, 0.9878, 0.9879, 0.988, 0.9881, 0.9882, 0.9883, 0.9884, 0.9884999999999999, 0.9886, 0.9887, 0.9888, 0.9889, 0.989, 0.9891000000000001, 0.9892000000000001, 0.9893000000000001, 0.9894000000000001, 0.9895, 0.9896, 0.9897, 0.9898, 0.9899, 0.99, 0.9901, 0.9902, 0.9903, 0.9904, 0.9904999999999999, 0.9906, 0.9907, 0.9908, 0.9909, 0.991, 0.9911000000000001, 0.9912000000000001, 0.9913000000000001, 0.9914000000000001, 0.9915, 0.9916, 0.9917, 0.9918, 0.9919, 0.992, 0.9921, 0.9922, 0.9923, 0.9924, 0.9924999999999999, 0.9925999999999999, 0.9927, 0.9928, 0.9929, 0.993, 0.9931, 0.9932000000000001, 0.9933000000000001, 0.9934000000000001, 0.9935, 0.9936, 0.9937, 0.9938, 0.9939, 0.994, 0.9941, 0.9942, 0.9943, 0.9944, 0.9944999999999999, 0.9945999999999999, 0.9947, 0.9948, 0.9949, 0.995, 0.9951, 0.9952000000000001, 0.9953000000000001, 0.9954000000000001, 0.9955, 0.9956, 0.9957, 0.9958, 0.9959, 0.996, 0.9961, 0.9962, 0.9963, 0.9964, 0.9964999999999999, 0.9965999999999999, 0.9967, 0.9968, 0.9969, 0.997, 0.9971, 0.9972000000000001, 0.9973000000000001, 0.9974000000000001, 0.9975, 0.9976, 0.9977, 0.9978, 0.9979, 0.998, 0.9981, 0.9982, 0.9983, 0.9984, 0.9984999999999999, 0.9985999999999999, 0.9987, 0.9988, 0.9989, 0.999, 0.9991, 0.9992000000000001, 0.9993000000000001, 0.9994000000000001, 0.9995, 0.9996, 0.9997, 0.9998, 0.9999, 1.0};
    for(int k = 0; k < trials+1; k++){
      probcuts[k] = probcuts_temp[k];
    }
  }

  if(opt == 1){
    //Float_t probcuts_temp[trials+1] = {0.95, 0.95001, 0.95002, 0.95003, 0.95004, 0.95005, 0.95006, 0.95007, 0.95008, 0.95009, 0.9501, 0.95011, 0.95012, 0.95013, 0.95014, 0.95015, 0.95016, 0.95017, 0.95018, 0.95019, 0.9502, 0.95021, 0.95022, 0.95023, 0.95024, 0.95025, 0.95026, 0.95027, 0.95028, 0.95029, 0.95030, 0.95031, 0.95032, 0.95033, 0.95034, 0.95035, 0.95036, 0.95037, 0.95038, 0.95039, 0.9504, 0.95041, 0.95042, 0.95043, 0.95044, 0.95045, 0.95046, 0.95047, 0.95048, 0.95049, 0.95050, 0.95051, 0.95052, 0.95053, 0.95054, 0.95055, 0.95056, 0.95057, 0.95058, 0.95059, 0.9506, 0.95061, 0.95062, 0.95063, 0.95064, 0.95065, 0.95066, 0.95067, 0.95068, 0.95069, 0.9507, 0.95071, 0.95072, 0.95073, 0.95074, 0.95075, 0.95076, 0.95077, 0.95078, 0.95079, 0.9508, 0.95081, 0.95082, 0.95083, 0.95084, 0.95085, 0.95086, 0.95087, 0.95088, 0.95089, 0.9509, 0.95091, 0.95092, 0.95093, 0.95094, 0.95095, 0.95096, 0.95097, 0.95098, 0.95099, 0.951, 0.95101, 0.95102, 0.95103, 0.95104, 0.95105, 0.95106, 0.95107, 0.95108, 0.95109, 0.9511, 0.95111, 0.95112, 0.95113, 0.95114, 0.95115, 0.95116, 0.95117, 0.95118, 0.95119, 0.9512, 0.95121, 0.95122, 0.95123, 0.95124, 0.95125, 0.95126, 0.95127, 0.95128, 0.95129, 0.95130, 0.95131, 0.95132, 0.95133, 0.95134, 0.95135, 0.95136, 0.95137, 0.95138, 0.95139, 0.9514, 0.95141, 0.95142, 0.95143, 0.95144, 0.95145, 0.95146, 0.95147, 0.95148, 0.95149, 0.9515, 0.95151, 0.95152, 0.95153, 0.95154, 0.95155, 0.95156, 0.95157, 0.95158, 0.95159, 0.9516, 0.95161, 0.95162, 0.95163, 0.95164, 0.95165, 0.95166, 0.95167, 0.95168, 0.95169, 0.9517, 0.95171, 0.95172, 0.95173, 0.95174, 0.95175, 0.95176, 0.95177, 0.95178, 0.95179, 0.9518, 0.95181, 0.95182, 0.95183, 0.95184, 0.95185, 0.95186, 0.95187, 0.95188, 0.95189, 0.9519, 0.95191, 0.95192, 0.95193, 0.95194, 0.95195, 0.95196, 0.95197, 0.95198, 0.95199, 0.952, 0.95201, 0.95202, 0.95203, 0.95204, 0.95205, 0.95206, 0.95207, 0.95208, 0.95209, 0.9521, 0.95211, 0.95212, 0.95213, 0.95214, 0.95215, 0.95216, 0.95217, 0.95218, 0.95219, 0.9522, 0.95221, 0.95222, 0.95223, 0.95224, 0.95225, 0.95226, 0.95227, 0.95228, 0.95229, 0.95230, 0.95231, 0.95232, 0.95233, 0.95234, 0.95235, 0.95236, 0.95237, 0.95238, 0.95239, 0.9524, 0.95241, 0.95242, 0.95243, 0.95244, 0.95245, 0.95246, 0.95247, 0.95248, 0.95249, 0.95250, 0.95251, 0.95252, 0.95253, 0.95254, 0.95255, 0.95256, 0.95257, 0.95258, 0.95259, 0.9526, 0.95261, 0.95262, 0.95263, 0.95264, 0.95265, 0.95266, 0.95267, 0.95268, 0.95269, 0.9527, 0.95271, 0.95272, 0.95273, 0.95274, 0.95275, 0.95276, 0.95277, 0.95278, 0.95279, 0.9528, 0.95281, 0.95282, 0.95283, 0.95284, 0.95285, 0.95286, 0.95287, 0.95288, 0.95289, 0.9529, 0.95291, 0.95292, 0.95293, 0.95294, 0.95295, 0.95296, 0.95297, 0.95298, 0.95299, 0.953, 0.95301, 0.95302, 0.95303, 0.95304, 0.95305, 0.95306, 0.95307, 0.95308, 0.95309, 0.9531, 0.95311, 0.95312, 0.95313, 0.95314, 0.95315, 0.95316, 0.95317, 0.95318, 0.95319, 0.9532, 0.95321, 0.95322, 0.95323, 0.95324, 0.95325, 0.95326, 0.95327, 0.95328, 0.95329, 0.95330, 0.95331, 0.95332, 0.95333, 0.95334, 0.95335, 0.95336, 0.95337, 0.95338, 0.95339, 0.9534, 0.95341, 0.95342, 0.95343, 0.95344, 0.95345, 0.95346, 0.95347, 0.95348, 0.95349, 0.9535, 0.95351, 0.95352, 0.95353, 0.95354, 0.95355, 0.95356, 0.95357, 0.95358, 0.95359, 0.9536, 0.95361, 0.95362, 0.95363, 0.95364, 0.95365, 0.95366, 0.95367, 0.95368, 0.95369, 0.9537, 0.95371, 0.95372, 0.95373, 0.95374, 0.95375, 0.95376, 0.95377, 0.95378, 0.95379, 0.9538, 0.95381, 0.95382, 0.95383, 0.95384, 0.95385, 0.95386, 0.95387, 0.95388, 0.95389, 0.9539, 0.95391, 0.95392, 0.95393, 0.95394, 0.95395, 0.95396, 0.95397, 0.95398, 0.95399, 0.954, 0.95401, 0.95402, 0.95403, 0.95404, 0.95405, 0.95406, 0.95407, 0.95408, 0.95409, 0.9541, 0.95411, 0.95412, 0.95413, 0.95414, 0.95415, 0.95416, 0.95417, 0.95418, 0.95419, 0.9542, 0.95421, 0.95422, 0.95423, 0.95424, 0.95425, 0.95426, 0.95427, 0.95428, 0.95429, 0.95430, 0.95431, 0.95432, 0.95433, 0.95434, 0.95435, 0.95436, 0.95437, 0.95438, 0.95439, 0.9544, 0.95441, 0.95442, 0.95443, 0.95444, 0.95445, 0.95446, 0.95447, 0.95448, 0.95449, 0.95450, 0.95451, 0.95452, 0.95453, 0.95454, 0.95455, 0.95456, 0.95457, 0.95458, 0.95459, 0.9546, 0.95461, 0.95462, 0.95463, 0.95464, 0.95465, 0.95466, 0.95467, 0.95468, 0.95469, 0.9547, 0.95471, 0.95472, 0.95473, 0.95474, 0.95475, 0.95476, 0.95477, 0.95478, 0.95479, 0.9548, 0.95481, 0.95482, 0.95483, 0.95484, 0.95485, 0.95486, 0.95487, 0.95488, 0.95489, 0.9549, 0.95491, 0.95492, 0.95493, 0.95494, 0.95495, 0.95496, 0.95497, 0.95498, 0.95499, 0.95500, 0.95501, 0.95502, 0.95503, 0.95504, 0.95505, 0.95506, 0.95507, 0.95508, 0.95509, 0.9551, 0.95511, 0.95512, 0.95513, 0.95514, 0.95515, 0.95516, 0.95517, 0.95518, 0.95519, 0.9552, 0.95521, 0.95522, 0.95523, 0.95524, 0.95525, 0.95526, 0.95527, 0.95528, 0.95529, 0.95530, 0.95531, 0.95532, 0.95533, 0.95534, 0.95535, 0.95536, 0.95537, 0.95538, 0.95539, 0.9554, 0.95541, 0.95542, 0.95543, 0.95544, 0.95545, 0.95546, 0.95547, 0.95548, 0.95549, 0.9555, 0.95551, 0.95552, 0.95553, 0.95554, 0.95555, 0.95556, 0.95557, 0.95558, 0.95559, 0.9556, 0.95561, 0.95562, 0.95563, 0.95564, 0.95565, 0.95566, 0.95567, 0.95568, 0.95569, 0.9557, 0.95571, 0.95572, 0.95573, 0.95574, 0.95575, 0.95576, 0.95577, 0.95578, 0.95579, 0.9558, 0.95581, 0.95582, 0.95583, 0.95584, 0.95585, 0.95586, 0.95587, 0.95588, 0.95589, 0.95590, 0.95591, 0.95592, 0.95593, 0.95594, 0.95595, 0.95596, 0.95597, 0.95598, 0.95599, 0.956, 0.95601, 0.95602, 0.95603, 0.95604, 0.95605, 0.95606, 0.95607, 0.95608, 0.95609, 0.9561, 0.95611, 0.95612, 0.95613, 0.95614, 0.95615, 0.95616, 0.95617, 0.95618, 0.95619, 0.9562, 0.95621, 0.95622, 0.95623, 0.95624, 0.95625, 0.95626, 0.95627, 0.95628, 0.95629, 0.95630, 0.95631, 0.95632, 0.95633, 0.95634, 0.95635, 0.95636, 0.95637, 0.95638, 0.95639, 0.9564, 0.95641, 0.95642, 0.95643, 0.95644, 0.95645, 0.95646, 0.95647, 0.95648, 0.95649, 0.95650, 0.95651, 0.95652, 0.95653, 0.95654, 0.95655, 0.95656, 0.95657, 0.95658, 0.95659, 0.9566, 0.95661, 0.95662, 0.95663, 0.95664, 0.95665, 0.95666, 0.95667, 0.95668, 0.95669, 0.9567, 0.95671, 0.95672, 0.95673, 0.95674, 0.95675, 0.95676, 0.95677, 0.95678, 0.95679, 0.9568, 0.95681, 0.95682, 0.95683, 0.95684, 0.95685, 0.95686, 0.95687, 0.95688, 0.95689, 0.9569, 0.95691, 0.95692, 0.95693, 0.95694, 0.95695, 0.95696, 0.95697, 0.95698, 0.95699, 0.957, 0.95701, 0.95702, 0.95703, 0.95704, 0.95705, 0.95706, 0.95707, 0.95708, 0.95709, 0.9571, 0.95711, 0.95712, 0.95713, 0.95714, 0.95715, 0.95716, 0.95717, 0.95718, 0.95719, 0.9572, 0.95721, 0.95722, 0.95723, 0.95724, 0.95725, 0.95726, 0.95727, 0.95728, 0.95729, 0.95730, 0.95731, 0.95732, 0.95733, 0.95734, 0.95735, 0.95736, 0.95737, 0.95738, 0.95739, 0.9574, 0.95741, 0.95742, 0.95743, 0.95744, 0.95745, 0.95746, 0.95747, 0.95748, 0.95749, 0.9575, 0.95751, 0.95752, 0.95753, 0.95754, 0.95755, 0.95756, 0.95757, 0.95758, 0.95759, 0.9576, 0.95761, 0.95762, 0.95763, 0.95764, 0.95765, 0.95766, 0.95767, 0.95768, 0.95769, 0.9577, 0.95771, 0.95772, 0.95773, 0.95774, 0.95775, 0.95776, 0.95777, 0.95778, 0.95779, 0.9578, 0.95781, 0.95782, 0.95783, 0.95784, 0.95785, 0.95786, 0.95787, 0.95788, 0.95789, 0.9579, 0.95791, 0.95792, 0.95793, 0.95794, 0.95795, 0.95796, 0.95797, 0.95798, 0.95799, 0.958, 0.95801, 0.95802, 0.95803, 0.95804, 0.95805, 0.95806, 0.95807, 0.95808, 0.95809, 0.9581, 0.95811, 0.95812, 0.95813, 0.95814, 0.95815, 0.95816, 0.95817, 0.95818, 0.95819, 0.9582, 0.95821, 0.95822, 0.95823, 0.95824, 0.95825, 0.95826, 0.95827, 0.95828, 0.95829, 0.95830, 0.95831, 0.95832, 0.95833, 0.95834, 0.95835, 0.95836, 0.95837, 0.95838, 0.95839, 0.9584, 0.95841, 0.95842, 0.95843, 0.95844, 0.95845, 0.95846, 0.95847, 0.95848, 0.95849, 0.95850, 0.95851, 0.95852, 0.95853, 0.95854, 0.95855, 0.95856, 0.95857, 0.95858, 0.95859, 0.9586, 0.95861, 0.95862, 0.95863, 0.95864, 0.95865, 0.95866, 0.95867, 0.95868, 0.95869, 0.9587, 0.95871, 0.95872, 0.95873, 0.95874, 0.95875, 0.95876, 0.95877, 0.95878, 0.95879, 0.9588, 0.95881, 0.95882, 0.95883, 0.95884, 0.95885, 0.95886, 0.95887, 0.95888, 0.95889, 0.9589, 0.95891, 0.95892, 0.95893, 0.95894, 0.95895, 0.95896, 0.95897, 0.95898, 0.95899, 0.959, 0.95901, 0.95902, 0.95903, 0.95904, 0.95905, 0.95906, 0.95907, 0.95908, 0.95909, 0.9591, 0.95911, 0.95912, 0.95913, 0.95914, 0.95915, 0.95916, 0.95917, 0.95918, 0.95919, 0.9592, 0.95921, 0.95922, 0.95923, 0.95924, 0.95925, 0.95926, 0.95927, 0.95928, 0.95929, 0.95930, 0.95931, 0.95932, 0.95933, 0.95934, 0.95935, 0.95936, 0.95937, 0.95938, 0.95939, 0.9594, 0.95941, 0.95942, 0.95943, 0.95944, 0.95945, 0.95946, 0.95947, 0.95948, 0.95949, 0.9595, 0.95951, 0.95952, 0.95953, 0.95954, 0.95955, 0.95956, 0.95957, 0.95958, 0.95959, 0.9596, 0.95961, 0.95962, 0.95963, 0.95964, 0.95965, 0.95966, 0.95967, 0.95968, 0.95969, 0.9597, 0.95971, 0.95972, 0.95973, 0.95974, 0.95975, 0.95976, 0.95977, 0.95978, 0.95979, 0.9598, 0.95981, 0.95982, 0.95983, 0.95984, 0.95985, 0.95986, 0.95987, 0.95988, 0.95989, 0.9599, 0.95991, 0.95992, 0.95993, 0.95994, 0.95995, 0.95996, 0.95997, 0.95998, 0.95999, 0.96, 0.96001, 0.96002, 0.96003, 0.96004, 0.96005, 0.96006, 0.96007, 0.96008, 0.96009, 0.9601, 0.96011, 0.96012, 0.96013, 0.96014, 0.96015, 0.96016, 0.96017, 0.96018, 0.96019, 0.9602, 0.96021, 0.96022, 0.96023, 0.96024, 0.96025, 0.96026, 0.96027, 0.96028, 0.96029, 0.96030, 0.96031, 0.96032, 0.96033, 0.96034, 0.96035, 0.96036, 0.96037, 0.96038, 0.96039, 0.9604, 0.96041, 0.96042, 0.96043, 0.96044, 0.96045, 0.96046, 0.96047, 0.96048, 0.96049, 0.96050, 0.96051, 0.96052, 0.96053, 0.96054, 0.96055, 0.96056, 0.96057, 0.96058, 0.96059, 0.9606, 0.96061, 0.96062, 0.96063, 0.96064, 0.96065, 0.96066, 0.96067, 0.96068, 0.96069, 0.9607, 0.96071, 0.96072, 0.96073, 0.96074, 0.96075, 0.96076, 0.96077, 0.96078, 0.96079, 0.9608, 0.96081, 0.96082, 0.96083, 0.96084, 0.96085, 0.96086, 0.96087, 0.96088, 0.96089, 0.9609, 0.96091, 0.96092, 0.96093, 0.96094, 0.96095, 0.96096, 0.96097, 0.96098, 0.96099, 0.961, 0.96101, 0.96102, 0.96103, 0.96104, 0.96105, 0.96106, 0.96107, 0.96108, 0.96109, 0.9611, 0.96111, 0.96112, 0.96113, 0.96114, 0.96115, 0.96116, 0.96117, 0.96118, 0.96119, 0.9612, 0.96121, 0.96122, 0.96123, 0.96124, 0.96125, 0.96126, 0.96127, 0.96128, 0.96129, 0.96130, 0.96131, 0.96132, 0.96133, 0.96134, 0.96135, 0.96136, 0.96137, 0.96138, 0.96139, 0.9614, 0.96141, 0.96142, 0.96143, 0.96144, 0.96145, 0.96146, 0.96147, 0.96148, 0.96149, 0.9615, 0.96151, 0.96152, 0.96153, 0.96154, 0.96155, 0.96156, 0.96157, 0.96158, 0.96159, 0.9616, 0.96161, 0.96162, 0.96163, 0.96164, 0.96165, 0.96166, 0.96167, 0.96168, 0.96169, 0.9617, 0.96171, 0.96172, 0.96173, 0.96174, 0.96175, 0.96176, 0.96177, 0.96178, 0.96179, 0.9618, 0.96181, 0.96182, 0.96183, 0.96184, 0.96185, 0.96186, 0.96187, 0.96188, 0.96189, 0.9619, 0.96191, 0.96192, 0.96193, 0.96194, 0.96195, 0.96196, 0.96197, 0.96198, 0.96199, 0.962, 0.96201, 0.96202, 0.96203, 0.96204, 0.96205, 0.96206, 0.96207, 0.96208, 0.96209, 0.9621, 0.96211, 0.96212, 0.96213, 0.96214, 0.96215, 0.96216, 0.96217, 0.96218, 0.96219, 0.9622, 0.96221, 0.96222, 0.96223, 0.96224, 0.96225, 0.96226, 0.96227, 0.96228, 0.96229, 0.96230, 0.96231, 0.96232, 0.96233, 0.96234, 0.96235, 0.96236, 0.96237, 0.96238, 0.96239, 0.9624, 0.96241, 0.96242, 0.96243, 0.96244, 0.96245, 0.96246, 0.96247, 0.96248, 0.96249, 0.96250, 0.96251, 0.96252, 0.96253, 0.96254, 0.96255, 0.96256, 0.96257, 0.96258, 0.96259, 0.9626, 0.96261, 0.96262, 0.96263, 0.96264, 0.96265, 0.96266, 0.96267, 0.96268, 0.96269, 0.9627, 0.96271, 0.96272, 0.96273, 0.96274, 0.96275, 0.96276, 0.96277, 0.96278, 0.96279, 0.9628, 0.96281, 0.96282, 0.96283, 0.96284, 0.96285, 0.96286, 0.96287, 0.96288, 0.96289, 0.9629, 0.96291, 0.96292, 0.96293, 0.96294, 0.96295, 0.96296, 0.96297, 0.96298, 0.96299, 0.963, 0.96301, 0.96302, 0.96303, 0.96304, 0.96305, 0.96306, 0.96307, 0.96308, 0.96309, 0.9631, 0.96311, 0.96312, 0.96313, 0.96314, 0.96315, 0.96316, 0.96317, 0.96318, 0.96319, 0.9632, 0.96321, 0.96322, 0.96323, 0.96324, 0.96325, 0.96326, 0.96327, 0.96328, 0.96329, 0.96330, 0.96331, 0.96332, 0.96333, 0.96334, 0.96335, 0.96336, 0.96337, 0.96338, 0.96339, 0.9634, 0.96341, 0.96342, 0.96343, 0.96344, 0.96345, 0.96346, 0.96347, 0.96348, 0.96349, 0.9635, 0.96351, 0.96352, 0.96353, 0.96354, 0.96355, 0.96356, 0.96357, 0.96358, 0.96359, 0.9636, 0.96361, 0.96362, 0.96363, 0.96364, 0.96365, 0.96366, 0.96367, 0.96368, 0.96369, 0.9637, 0.96371, 0.96372, 0.96373, 0.96374, 0.96375, 0.96376, 0.96377, 0.96378, 0.96379, 0.9638, 0.96381, 0.96382, 0.96383, 0.96384, 0.96385, 0.96386, 0.96387, 0.96388, 0.96389, 0.9639, 0.96391, 0.96392, 0.96393, 0.96394, 0.96395, 0.96396, 0.96397, 0.96398, 0.96399, 0.964, 0.96401, 0.96402, 0.96403, 0.96404, 0.96405, 0.96406, 0.96407, 0.96408, 0.96409, 0.9641, 0.96411, 0.96412, 0.96413, 0.96414, 0.96415, 0.96416, 0.96417, 0.96418, 0.96419, 0.9642, 0.96421, 0.96422, 0.96423, 0.96424, 0.96425, 0.96426, 0.96427, 0.96428, 0.96429, 0.96430, 0.96431, 0.96432, 0.96433, 0.96434, 0.96435, 0.96436, 0.96437, 0.96438, 0.96439, 0.9644, 0.96441, 0.96442, 0.96443, 0.96444, 0.96445, 0.96446, 0.96447, 0.96448, 0.96449, 0.96450, 0.96451, 0.96452, 0.96453, 0.96454, 0.96455, 0.96456, 0.96457, 0.96458, 0.96459, 0.9646, 0.96461, 0.96462, 0.96463, 0.96464, 0.96465, 0.96466, 0.96467, 0.96468, 0.96469, 0.9647, 0.96471, 0.96472, 0.96473, 0.96474, 0.96475, 0.96476, 0.96477, 0.96478, 0.96479, 0.9648, 0.96481, 0.96482, 0.96483, 0.96484, 0.96485, 0.96486, 0.96487, 0.96488, 0.96489, 0.9649, 0.96491, 0.96492, 0.96493, 0.96494, 0.96495, 0.96496, 0.96497, 0.96498, 0.96499, 0.965, 0.96501, 0.96502, 0.96503, 0.96504, 0.96505, 0.96506, 0.96507, 0.96508, 0.96509, 0.9651, 0.96511, 0.96512, 0.96513, 0.96514, 0.96515, 0.96516, 0.96517, 0.96518, 0.96519, 0.9652, 0.96521, 0.96522, 0.96523, 0.96524, 0.96525, 0.96526, 0.96527, 0.96528, 0.96529, 0.96530, 0.96531, 0.96532, 0.96533, 0.96534, 0.96535, 0.96536, 0.96537, 0.96538, 0.96539, 0.9654, 0.96541, 0.96542, 0.96543, 0.96544, 0.96545, 0.96546, 0.96547, 0.96548, 0.96549, 0.9655, 0.96551, 0.96552, 0.96553, 0.96554, 0.96555, 0.96556, 0.96557, 0.96558, 0.96559, 0.9656, 0.96561, 0.96562, 0.96563, 0.96564, 0.96565, 0.96566, 0.96567, 0.96568, 0.96569, 0.9657, 0.96571, 0.96572, 0.96573, 0.96574, 0.96575, 0.96576, 0.96577, 0.96578, 0.96579, 0.9658, 0.96581, 0.96582, 0.96583, 0.96584, 0.96585, 0.96586, 0.96587, 0.96588, 0.96589, 0.9659, 0.96591, 0.96592, 0.96593, 0.96594, 0.96595, 0.96596, 0.96597, 0.96598, 0.96599, 0.966, 0.96601, 0.96602, 0.96603, 0.96604, 0.96605, 0.96606, 0.96607, 0.96608, 0.96609, 0.9661, 0.96611, 0.96612, 0.96613, 0.96614, 0.96615, 0.96616, 0.96617, 0.96618, 0.96619, 0.9662, 0.96621, 0.96622, 0.96623, 0.96624, 0.96625, 0.96626, 0.96627, 0.96628, 0.96629, 0.96630, 0.96631, 0.96632, 0.96633, 0.96634, 0.96635, 0.96636, 0.96637, 0.96638, 0.96639, 0.9664, 0.96641, 0.96642, 0.96643, 0.96644, 0.96645, 0.96646, 0.96647, 0.96648, 0.96649, 0.96650, 0.96651, 0.96652, 0.96653, 0.96654, 0.96655, 0.96656, 0.96657, 0.96658, 0.96659, 0.9666, 0.96661, 0.96662, 0.96663, 0.96664, 0.96665, 0.96666, 0.96667, 0.96668, 0.96669, 0.9667, 0.96671, 0.96672, 0.96673, 0.96674, 0.96675, 0.96676, 0.96677, 0.96678, 0.96679, 0.9668, 0.96681, 0.96682, 0.96683, 0.96684, 0.96685, 0.96686, 0.96687, 0.96688, 0.96689, 0.9669, 0.96691, 0.96692, 0.96693, 0.96694, 0.96695, 0.96696, 0.96697, 0.96698, 0.96699, 0.967, 0.96701, 0.96702, 0.96703, 0.96704, 0.96705, 0.96706, 0.96707, 0.96708, 0.96709, 0.9671, 0.96711, 0.96712, 0.96713, 0.96714, 0.96715, 0.96716, 0.96717, 0.96718, 0.96719, 0.9672, 0.96721, 0.96722, 0.96723, 0.96724, 0.96725, 0.96726, 0.96727, 0.96728, 0.96729, 0.96730, 0.96731, 0.96732, 0.96733, 0.96734, 0.96735, 0.96736, 0.96737, 0.96738, 0.96739, 0.9674, 0.96741, 0.96742, 0.96743, 0.96744, 0.96745, 0.96746, 0.96747, 0.96748, 0.96749, 0.9675, 0.96751, 0.96752, 0.96753, 0.96754, 0.96755, 0.96756, 0.96757, 0.96758, 0.96759, 0.9676, 0.96761, 0.96762, 0.96763, 0.96764, 0.96765, 0.96766, 0.96767, 0.96768, 0.96769, 0.9677, 0.96771, 0.96772, 0.96773, 0.96774, 0.96775, 0.96776, 0.96777, 0.96778, 0.96779, 0.9678, 0.96781, 0.96782, 0.96783, 0.96784, 0.96785, 0.96786, 0.96787, 0.96788, 0.96789, 0.9679, 0.96791, 0.96792, 0.96793, 0.96794, 0.96795, 0.96796, 0.96797, 0.96798, 0.96799, 0.968, 0.96801, 0.96802, 0.96803, 0.96804, 0.96805, 0.96806, 0.96807, 0.96808, 0.96809, 0.9681, 0.96811, 0.96812, 0.96813, 0.96814, 0.96815, 0.96816, 0.96817, 0.96818, 0.96819, 0.9682, 0.96821, 0.96822, 0.96823, 0.96824, 0.96825, 0.96826, 0.96827, 0.96828, 0.96829, 0.96830, 0.96831, 0.96832, 0.96833, 0.96834, 0.96835, 0.96836, 0.96837, 0.96838, 0.96839, 0.9684, 0.96841, 0.96842, 0.96843, 0.96844, 0.96845, 0.96846, 0.96847, 0.96848, 0.96849, 0.96850, 0.96851, 0.96852, 0.96853, 0.96854, 0.96855, 0.96856, 0.96857, 0.96858, 0.96859, 0.9686, 0.96861, 0.96862, 0.96863, 0.96864, 0.96865, 0.96866, 0.96867, 0.96868, 0.96869, 0.9687, 0.96871, 0.96872, 0.96873, 0.96874, 0.96875, 0.96876, 0.96877, 0.96878, 0.96879, 0.9688, 0.96881, 0.96882, 0.96883, 0.96884, 0.96885, 0.96886, 0.96887, 0.96888, 0.96889, 0.9689, 0.96891, 0.96892, 0.96893, 0.96894, 0.96895, 0.96896, 0.96897, 0.96898, 0.96899, 0.969, 0.96901, 0.96902, 0.96903, 0.96904, 0.96905, 0.96906, 0.96907, 0.96908, 0.96909, 0.9691, 0.96911, 0.96912, 0.96913, 0.96914, 0.96915, 0.96916, 0.96917, 0.96918, 0.96919, 0.9692, 0.96921, 0.96922, 0.96923, 0.96924, 0.96925, 0.96926, 0.96927, 0.96928, 0.96929, 0.96930, 0.96931, 0.96932, 0.96933, 0.96934, 0.96935, 0.96936, 0.96937, 0.96938, 0.96939, 0.9694, 0.96941, 0.96942, 0.96943, 0.96944, 0.96945, 0.96946, 0.96947, 0.96948, 0.96949, 0.9695, 0.96951, 0.96952, 0.96953, 0.96954, 0.96955, 0.96956, 0.96957, 0.96958, 0.96959, 0.9696, 0.96961, 0.96962, 0.96963, 0.96964, 0.96965, 0.96966, 0.96967, 0.96968, 0.96969, 0.9697, 0.96971, 0.96972, 0.96973, 0.96974, 0.96975, 0.96976, 0.96977, 0.96978, 0.96979, 0.9698, 0.96981, 0.96982, 0.96983, 0.96984, 0.96985, 0.96986, 0.96987, 0.96988, 0.96989, 0.9699, 0.96991, 0.96992, 0.96993, 0.96994, 0.96995, 0.96996, 0.96997, 0.96998, 0.96999, 0.97, 0.97001, 0.97002, 0.97003, 0.97004, 0.97005, 0.97006, 0.97007, 0.97008, 0.97009, 0.9701, 0.97011, 0.97012, 0.97013, 0.97014, 0.97015, 0.97016, 0.97017, 0.97018, 0.97019, 0.9702, 0.97021, 0.97022, 0.97023, 0.97024, 0.97025, 0.97026, 0.97027, 0.97028, 0.97029, 0.97030, 0.97031, 0.97032, 0.97033, 0.97034, 0.97035, 0.97036, 0.97037, 0.97038, 0.97039, 0.9704, 0.97041, 0.97042, 0.97043, 0.97044, 0.97045, 0.97046, 0.97047, 0.97048, 0.97049, 0.97050, 0.97051, 0.97052, 0.97053, 0.97054, 0.97055, 0.97056, 0.97057, 0.97058, 0.97059, 0.9706, 0.97061, 0.97062, 0.97063, 0.97064, 0.97065, 0.97066, 0.97067, 0.97068, 0.97069, 0.9707, 0.97071, 0.97072, 0.97073, 0.97074, 0.97075, 0.97076, 0.97077, 0.97078, 0.97079, 0.9708, 0.97081, 0.97082, 0.97083, 0.97084, 0.97085, 0.97086, 0.97087, 0.97088, 0.97089, 0.9709, 0.97091, 0.97092, 0.97093, 0.97094, 0.97095, 0.97096, 0.97097, 0.97098, 0.97099, 0.971, 0.97101, 0.97102, 0.97103, 0.97104, 0.97105, 0.97106, 0.97107, 0.97108, 0.97109, 0.9711, 0.97111, 0.97112, 0.97113, 0.97114, 0.97115, 0.97116, 0.97117, 0.97118, 0.97119, 0.9712, 0.97121, 0.97122, 0.97123, 0.97124, 0.97125, 0.97126, 0.97127, 0.97128, 0.97129, 0.97130, 0.97131, 0.97132, 0.97133, 0.97134, 0.97135, 0.97136, 0.97137, 0.97138, 0.97139, 0.9714, 0.97141, 0.97142, 0.97143, 0.97144, 0.97145, 0.97146, 0.97147, 0.97148, 0.97149, 0.9715, 0.97151, 0.97152, 0.97153, 0.97154, 0.97155, 0.97156, 0.97157, 0.97158, 0.97159, 0.9716, 0.97161, 0.97162, 0.97163, 0.97164, 0.97165, 0.97166, 0.97167, 0.97168, 0.97169, 0.9717, 0.97171, 0.97172, 0.97173, 0.97174, 0.97175, 0.97176, 0.97177, 0.97178, 0.97179, 0.9718, 0.97181, 0.97182, 0.97183, 0.97184, 0.97185, 0.97186, 0.97187, 0.97188, 0.97189, 0.9719, 0.97191, 0.97192, 0.97193, 0.97194, 0.97195, 0.97196, 0.97197, 0.97198, 0.97199, 0.972, 0.97201, 0.97202, 0.97203, 0.97204, 0.97205, 0.97206, 0.97207, 0.97208, 0.97209, 0.9721, 0.97211, 0.97212, 0.97213, 0.97214, 0.97215, 0.97216, 0.97217, 0.97218, 0.97219, 0.9722, 0.97221, 0.97222, 0.97223, 0.97224, 0.97225, 0.97226, 0.97227, 0.97228, 0.97229, 0.97230, 0.97231, 0.97232, 0.97233, 0.97234, 0.97235, 0.97236, 0.97237, 0.97238, 0.97239, 0.9724, 0.97241, 0.97242, 0.97243, 0.97244, 0.97245, 0.97246, 0.97247, 0.97248, 0.97249, 0.97250, 0.97251, 0.97252, 0.97253, 0.97254, 0.97255, 0.97256, 0.97257, 0.97258, 0.97259, 0.9726, 0.97261, 0.97262, 0.97263, 0.97264, 0.97265, 0.97266, 0.97267, 0.97268, 0.97269, 0.9727, 0.97271, 0.97272, 0.97273, 0.97274, 0.97275, 0.97276, 0.97277, 0.97278, 0.97279, 0.9728, 0.97281, 0.97282, 0.97283, 0.97284, 0.97285, 0.97286, 0.97287, 0.97288, 0.97289, 0.9729, 0.97291, 0.97292, 0.97293, 0.97294, 0.97295, 0.97296, 0.97297, 0.97298, 0.97299, 0.973, 0.97301, 0.97302, 0.97303, 0.97304, 0.97305, 0.97306, 0.97307, 0.97308, 0.97309, 0.9731, 0.97311, 0.97312, 0.97313, 0.97314, 0.97315, 0.97316, 0.97317, 0.97318, 0.97319, 0.9732, 0.97321, 0.97322, 0.97323, 0.97324, 0.97325, 0.97326, 0.97327, 0.97328, 0.97329, 0.97330, 0.97331, 0.97332, 0.97333, 0.97334, 0.97335, 0.97336, 0.97337, 0.97338, 0.97339, 0.9734, 0.97341, 0.97342, 0.97343, 0.97344, 0.97345, 0.97346, 0.97347, 0.97348, 0.97349, 0.9735, 0.97351, 0.97352, 0.97353, 0.97354, 0.97355, 0.97356, 0.97357, 0.97358, 0.97359, 0.9736, 0.97361, 0.97362, 0.97363, 0.97364, 0.97365, 0.97366, 0.97367, 0.97368, 0.97369, 0.9737, 0.97371, 0.97372, 0.97373, 0.97374, 0.97375, 0.97376, 0.97377, 0.97378, 0.97379, 0.9738, 0.97381, 0.97382, 0.97383, 0.97384, 0.97385, 0.97386, 0.97387, 0.97388, 0.97389, 0.9739, 0.97391, 0.97392, 0.97393, 0.97394, 0.97395, 0.97396, 0.97397, 0.97398, 0.97399, 0.974, 0.97401, 0.97402, 0.97403, 0.97404, 0.97405, 0.97406, 0.97407, 0.97408, 0.97409, 0.9741, 0.97411, 0.97412, 0.97413, 0.97414, 0.97415, 0.97416, 0.97417, 0.97418, 0.97419, 0.9742, 0.97421, 0.97422, 0.97423, 0.97424, 0.97425, 0.97426, 0.97427, 0.97428, 0.97429, 0.97430, 0.97431, 0.97432, 0.97433, 0.97434, 0.97435, 0.97436, 0.97437, 0.97438, 0.97439, 0.9744, 0.97441, 0.97442, 0.97443, 0.97444, 0.97445, 0.97446, 0.97447, 0.97448, 0.97449, 0.97450, 0.97451, 0.97452, 0.97453, 0.97454, 0.97455, 0.97456, 0.97457, 0.97458, 0.97459, 0.9746, 0.97461, 0.97462, 0.97463, 0.97464, 0.97465, 0.97466, 0.97467, 0.97468, 0.97469, 0.9747, 0.97471, 0.97472, 0.97473, 0.97474, 0.97475, 0.97476, 0.97477, 0.97478, 0.97479, 0.9748, 0.97481, 0.97482, 0.97483, 0.97484, 0.97485, 0.97486, 0.97487, 0.97488, 0.97489, 0.9749, 0.97491, 0.97492, 0.97493, 0.97494, 0.97495, 0.97496, 0.97497, 0.97498, 0.97499, 0.975, 0.97501, 0.97502, 0.97503, 0.97504, 0.97505, 0.97506, 0.97507, 0.97508, 0.97509, 0.9751, 0.97511, 0.97512, 0.97513, 0.97514, 0.97515, 0.97516, 0.97517, 0.97518, 0.97519, 0.9752, 0.97521, 0.97522, 0.97523, 0.97524, 0.97525, 0.97526, 0.97527, 0.97528, 0.97529, 0.9753, 0.97531, 0.97532, 0.97533, 0.97534, 0.97535, 0.97536, 0.97537, 0.97538, 0.97539, 0.9754, 0.97541, 0.97542, 0.97543, 0.97544, 0.97545, 0.97546, 0.97547, 0.97548, 0.97549, 0.9755, 0.97551, 0.97552, 0.97553, 0.97554, 0.97555, 0.97556, 0.97557, 0.97558, 0.97559, 0.9756, 0.97561, 0.97562, 0.97563, 0.97564, 0.97565, 0.97566, 0.97567, 0.97568, 0.97569, 0.9757, 0.97571, 0.97572, 0.97573, 0.97574, 0.97575, 0.97576, 0.97577, 0.97578, 0.97579, 0.9758, 0.97581, 0.97582, 0.97583, 0.97584, 0.97585, 0.97586, 0.97587, 0.97588, 0.97589, 0.9759, 0.97591, 0.97592, 0.97593, 0.97594, 0.97595, 0.97596, 0.97597, 0.97598, 0.97599, 0.976, 0.97601, 0.97602, 0.97603, 0.97604, 0.97605, 0.97606, 0.97607, 0.97608, 0.97609, 0.9761, 0.97611, 0.97612, 0.97613, 0.97614, 0.97615, 0.97616, 0.97617, 0.97618, 0.97619, 0.9762, 0.97621, 0.97622, 0.97623, 0.97624, 0.97625, 0.97626, 0.97627, 0.97628, 0.97629, 0.9763, 0.97631, 0.97632, 0.97633, 0.97634, 0.97635, 0.97636, 0.97637, 0.97638, 0.97639, 0.9764, 0.97641, 0.97642, 0.97643, 0.97644, 0.97645, 0.97646, 0.97647, 0.97648, 0.97649, 0.97650, 0.97651, 0.97652, 0.97653, 0.97654, 0.97655, 0.97656, 0.97657, 0.97658, 0.97659, 0.9766, 0.97661, 0.97662, 0.97663, 0.97664, 0.97665, 0.97666, 0.97667, 0.97668, 0.97669, 0.9767, 0.97671, 0.97672, 0.97673, 0.97674, 0.97675, 0.97676, 0.97677, 0.97678, 0.97679, 0.9768, 0.97681, 0.97682, 0.97683, 0.97684, 0.97685, 0.97686, 0.97687, 0.97688, 0.97689, 0.9769, 0.97691, 0.97692, 0.97693, 0.97694, 0.97695, 0.97696, 0.97697, 0.97698, 0.97699, 0.977, 0.97701, 0.97702, 0.97703, 0.97704, 0.97705, 0.97706, 0.97707, 0.97708, 0.97709, 0.9771, 0.97711, 0.97712, 0.97713, 0.97714, 0.97715, 0.97716, 0.97717, 0.97718, 0.97719, 0.9772, 0.97721, 0.97722, 0.97723, 0.97724, 0.97725, 0.97726, 0.97727, 0.97728, 0.97729, 0.9773, 0.97731, 0.97732, 0.97733, 0.97734, 0.97735, 0.97736, 0.97737, 0.97738, 0.97739, 0.9774, 0.97741, 0.97742, 0.97743, 0.97744, 0.97745, 0.97746, 0.97747, 0.97748, 0.97749, 0.9775, 0.97751, 0.97752, 0.97753, 0.97754, 0.97755, 0.97756, 0.97757, 0.97758, 0.97759, 0.9776, 0.97761, 0.97762, 0.97763, 0.97764, 0.97765, 0.97766, 0.97767, 0.97768, 0.97769, 0.9777, 0.97771, 0.97772, 0.97773, 0.97774, 0.97775, 0.97776, 0.97777, 0.97778, 0.97779, 0.9778, 0.97781, 0.97782, 0.97783, 0.97784, 0.97785, 0.97786, 0.97787, 0.97788, 0.97789, 0.9779, 0.97791, 0.97792, 0.97793, 0.97794, 0.97795, 0.97796, 0.97797, 0.97798, 0.97799, 0.978, 0.97801, 0.97802, 0.97803, 0.97804, 0.97805, 0.97806, 0.97807, 0.97808, 0.97809, 0.9781, 0.97811, 0.97812, 0.97813, 0.97814, 0.97815, 0.97816, 0.97817, 0.97818, 0.97819, 0.9782, 0.97821, 0.97822, 0.97823, 0.97824, 0.97825, 0.97826, 0.97827, 0.97828, 0.97829, 0.9783, 0.97831, 0.97832, 0.97833, 0.97834, 0.97835, 0.97836, 0.97837, 0.97838, 0.97839, 0.9784, 0.97841, 0.97842, 0.97843, 0.97844, 0.97845, 0.97846, 0.97847, 0.97848, 0.97849, 0.97850, 0.97851, 0.97852, 0.97853, 0.97854, 0.97855, 0.97856, 0.97857, 0.97858, 0.97859, 0.9786, 0.97861, 0.97862, 0.97863, 0.97864, 0.97865, 0.97866, 0.97867, 0.97868, 0.97869, 0.9787, 0.97871, 0.97872, 0.97873, 0.97874, 0.97875, 0.97876, 0.97877, 0.97878, 0.97879, 0.9788, 0.97881, 0.97882, 0.97883, 0.97884, 0.97885, 0.97886, 0.97887, 0.97888, 0.97889, 0.9789, 0.97891, 0.97892, 0.97893, 0.97894, 0.97895, 0.97896, 0.97897, 0.97898, 0.97899, 0.979, 0.97901, 0.97902, 0.97903, 0.97904, 0.97905, 0.97906, 0.97907, 0.97908, 0.97909, 0.9791, 0.97911, 0.97912, 0.97913, 0.97914, 0.97915, 0.97916, 0.97917, 0.97918, 0.97919, 0.9792, 0.97921, 0.97922, 0.97923, 0.97924, 0.97925, 0.97926, 0.97927, 0.97928, 0.97929, 0.9793, 0.97931, 0.97932, 0.97933, 0.97934, 0.97935, 0.97936, 0.97937, 0.97938, 0.97939, 0.9794, 0.97941, 0.97942, 0.97943, 0.97944, 0.97945, 0.97946, 0.97947, 0.97948, 0.97949, 0.9795, 0.97951, 0.97952, 0.97953, 0.97954, 0.97955, 0.97956, 0.97957, 0.97958, 0.97959, 0.9796, 0.97961, 0.97962, 0.97963, 0.97964, 0.97965, 0.97966, 0.97967, 0.97968, 0.97969, 0.9797, 0.97971, 0.97972, 0.97973, 0.97974, 0.97975, 0.97976, 0.97977, 0.97978, 0.97979, 0.9798, 0.97981, 0.97982, 0.97983, 0.97984, 0.97985, 0.97986, 0.97987, 0.97988, 0.97989, 0.9799, 0.97991, 0.97992, 0.97993, 0.97994, 0.97995, 0.97996, 0.97997, 0.97998, 0.97999, 0.98, 0.98001, 0.98002, 0.98003, 0.98004, 0.98005, 0.98006, 0.98007, 0.98008, 0.98009, 0.9801, 0.98011, 0.98012, 0.98013, 0.98014, 0.98015, 0.98016, 0.98017, 0.98018, 0.98019, 0.9802, 0.98021, 0.98022, 0.98023, 0.98024, 0.98025, 0.98026, 0.98027, 0.98028, 0.98029, 0.9803, 0.98031, 0.98032, 0.98033, 0.98034, 0.98035, 0.98036, 0.98037, 0.98038, 0.98039, 0.9804, 0.98041, 0.98042, 0.98043, 0.98044, 0.98045, 0.98046, 0.98047, 0.98048, 0.98049, 0.98050, 0.98051, 0.98052, 0.98053, 0.98054, 0.98055, 0.98056, 0.98057, 0.98058, 0.98059, 0.9806, 0.98061, 0.98062, 0.98063, 0.98064, 0.98065, 0.98066, 0.98067, 0.98068, 0.98069, 0.9807, 0.98071, 0.98072, 0.98073, 0.98074, 0.98075, 0.98076, 0.98077, 0.98078, 0.98079, 0.9808, 0.98081, 0.98082, 0.98083, 0.98084, 0.98085, 0.98086, 0.98087, 0.98088, 0.98089, 0.9809, 0.98091, 0.98092, 0.98093, 0.98094, 0.98095, 0.98096, 0.98097, 0.98098, 0.98099, 0.981, 0.98101, 0.98102, 0.98103, 0.98104, 0.98105, 0.98106, 0.98107, 0.98108, 0.98109, 0.9811, 0.98111, 0.98112, 0.98113, 0.98114, 0.98115, 0.98116, 0.98117, 0.98118, 0.98119, 0.9812, 0.98121, 0.98122, 0.98123, 0.98124, 0.98125, 0.98126, 0.98127, 0.98128, 0.98129, 0.9813, 0.98131, 0.98132, 0.98133, 0.98134, 0.98135, 0.98136, 0.98137, 0.98138, 0.98139, 0.9814, 0.98141, 0.98142, 0.98143, 0.98144, 0.98145, 0.98146, 0.98147, 0.98148, 0.98149, 0.9815, 0.98151, 0.98152, 0.98153, 0.98154, 0.98155, 0.98156, 0.98157, 0.98158, 0.98159, 0.9816, 0.98161, 0.98162, 0.98163, 0.98164, 0.98165, 0.98166, 0.98167, 0.98168, 0.98169, 0.9817, 0.98171, 0.98172, 0.98173, 0.98174, 0.98175, 0.98176, 0.98177, 0.98178, 0.98179, 0.9818, 0.98181, 0.98182, 0.98183, 0.98184, 0.98185, 0.98186, 0.98187, 0.98188, 0.98189, 0.9819, 0.98191, 0.98192, 0.98193, 0.98194, 0.98195, 0.98196, 0.98197, 0.98198, 0.98199, 0.982, 0.98201, 0.98202, 0.98203, 0.98204, 0.98205, 0.98206, 0.98207, 0.98208, 0.98209, 0.9821, 0.98211, 0.98212, 0.98213, 0.98214, 0.98215, 0.98216, 0.98217, 0.98218, 0.98219, 0.9822, 0.98221, 0.98222, 0.98223, 0.98224, 0.98225, 0.98226, 0.98227, 0.98228, 0.98229, 0.9823, 0.98231, 0.98232, 0.98233, 0.98234, 0.98235, 0.98236, 0.98237, 0.98238, 0.98239, 0.9824, 0.98241, 0.98242, 0.98243, 0.98244, 0.98245, 0.98246, 0.98247, 0.98248, 0.98249, 0.98250, 0.98251, 0.98252, 0.98253, 0.98254, 0.98255, 0.98256, 0.98257, 0.98258, 0.98259, 0.9826, 0.98261, 0.98262, 0.98263, 0.98264, 0.98265, 0.98266, 0.98267, 0.98268, 0.98269, 0.9827, 0.98271, 0.98272, 0.98273, 0.98274, 0.98275, 0.98276, 0.98277, 0.98278, 0.98279, 0.9828, 0.98281, 0.98282, 0.98283, 0.98284, 0.98285, 0.98286, 0.98287, 0.98288, 0.98289, 0.9829, 0.98291, 0.98292, 0.98293, 0.98294, 0.98295, 0.98296, 0.98297, 0.98298, 0.98299, 0.983, 0.98301, 0.98302, 0.98303, 0.98304, 0.98305, 0.98306, 0.98307, 0.98308, 0.98309, 0.9831, 0.98311, 0.98312, 0.98313, 0.98314, 0.98315, 0.98316, 0.98317, 0.98318, 0.98319, 0.9832, 0.98321, 0.98322, 0.98323, 0.98324, 0.98325, 0.98326, 0.98327, 0.98328, 0.98329, 0.9833, 0.98331, 0.98332, 0.98333, 0.98334, 0.98335, 0.98336, 0.98337, 0.98338, 0.98339, 0.9834, 0.98341, 0.98342, 0.98343, 0.98344, 0.98345, 0.98346, 0.98347, 0.98348, 0.98349, 0.9835, 0.98351, 0.98352, 0.98353, 0.98354, 0.98355, 0.98356, 0.98357, 0.98358, 0.98359, 0.9836, 0.98361, 0.98362, 0.98363, 0.98364, 0.98365, 0.98366, 0.98367, 0.98368, 0.98369, 0.9837, 0.98371, 0.98372, 0.98373, 0.98374, 0.98375, 0.98376, 0.98377, 0.98378, 0.98379, 0.9838, 0.98381, 0.98382, 0.98383, 0.98384, 0.98385, 0.98386, 0.98387, 0.98388, 0.98389, 0.9839, 0.98391, 0.98392, 0.98393, 0.98394, 0.98395, 0.98396, 0.98397, 0.98398, 0.98399, 0.984, 0.98401, 0.98402, 0.98403, 0.98404, 0.98405, 0.98406, 0.98407, 0.98408, 0.98409, 0.9841, 0.98411, 0.98412, 0.98413, 0.98414, 0.98415, 0.98416, 0.98417, 0.98418, 0.98419, 0.9842, 0.98421, 0.98422, 0.98423, 0.98424, 0.98425, 0.98426, 0.98427, 0.98428, 0.98429, 0.9843, 0.98431, 0.98432, 0.98433, 0.98434, 0.98435, 0.98436, 0.98437, 0.98438, 0.98439, 0.9844, 0.98441, 0.98442, 0.98443, 0.98444, 0.98445, 0.98446, 0.98447, 0.98448, 0.98449, 0.98450, 0.98451, 0.98452, 0.98453, 0.98454, 0.98455, 0.98456, 0.98457, 0.98458, 0.98459, 0.9846, 0.98461, 0.98462, 0.98463, 0.98464, 0.98465, 0.98466, 0.98467, 0.98468, 0.98469, 0.9847, 0.98471, 0.98472, 0.98473, 0.98474, 0.98475, 0.98476, 0.98477, 0.98478, 0.98479, 0.9848, 0.98481, 0.98482, 0.98483, 0.98484, 0.98485, 0.98486, 0.98487, 0.98488, 0.98489, 0.9849, 0.98491, 0.98492, 0.98493, 0.98494, 0.98495, 0.98496, 0.98497, 0.98498, 0.98499, 0.985, 0.98501, 0.98502, 0.98503, 0.98504, 0.98505, 0.98506, 0.98507, 0.98508, 0.98509, 0.9851, 0.98511, 0.98512, 0.98513, 0.98514, 0.98515, 0.98516, 0.98517, 0.98518, 0.98519, 0.9852, 0.98521, 0.98522, 0.98523, 0.98524, 0.98525, 0.98526, 0.98527, 0.98528, 0.98529, 0.9853, 0.98531, 0.98532, 0.98533, 0.98534, 0.98535, 0.98536, 0.98537, 0.98538, 0.98539, 0.9854, 0.98541, 0.98542, 0.98543, 0.98544, 0.98545, 0.98546, 0.98547, 0.98548, 0.98549, 0.9855, 0.98551, 0.98552, 0.98553, 0.98554, 0.98555, 0.98556, 0.98557, 0.98558, 0.98559, 0.9856, 0.98561, 0.98562, 0.98563, 0.98564, 0.98565, 0.98566, 0.98567, 0.98568, 0.98569, 0.9857, 0.98571, 0.98572, 0.98573, 0.98574, 0.98575, 0.98576, 0.98577, 0.98578, 0.98579, 0.9858, 0.98581, 0.98582, 0.98583, 0.98584, 0.98585, 0.98586, 0.98587, 0.98588, 0.98589, 0.9859, 0.98591, 0.98592, 0.98593, 0.98594, 0.98595, 0.98596, 0.98597, 0.98598, 0.98599, 0.986, 0.98601, 0.98602, 0.98603, 0.98604, 0.98605, 0.98606, 0.98607, 0.98608, 0.98609, 0.9861, 0.98611, 0.98612, 0.98613, 0.98614, 0.98615, 0.98616, 0.98617, 0.98618, 0.98619, 0.9862, 0.98621, 0.98622, 0.98623, 0.98624, 0.98625, 0.98626, 0.98627, 0.98628, 0.98629, 0.9863, 0.98631, 0.98632, 0.98633, 0.98634, 0.98635, 0.98636, 0.98637, 0.98638, 0.98639, 0.9864, 0.98641, 0.98642, 0.98643, 0.98644, 0.98645, 0.98646, 0.98647, 0.98648, 0.98649, 0.98650, 0.98651, 0.98652, 0.98653, 0.98654, 0.98655, 0.98656, 0.98657, 0.98658, 0.98659, 0.9866, 0.98661, 0.98662, 0.98663, 0.98664, 0.98665, 0.98666, 0.98667, 0.98668, 0.98669, 0.9867, 0.98671, 0.98672, 0.98673, 0.98674, 0.98675, 0.98676, 0.98677, 0.98678, 0.98679, 0.9868, 0.98681, 0.98682, 0.98683, 0.98684, 0.98685, 0.98686, 0.98687, 0.98688, 0.98689, 0.9869, 0.98691, 0.98692, 0.98693, 0.98694, 0.98695, 0.98696, 0.98697, 0.98698, 0.98699, 0.987, 0.98701, 0.98702, 0.98703, 0.98704, 0.98705, 0.98706, 0.98707, 0.98708, 0.98709, 0.9871, 0.98711, 0.98712, 0.98713, 0.98714, 0.98715, 0.98716, 0.98717, 0.98718, 0.98719, 0.9872, 0.98721, 0.98722, 0.98723, 0.98724, 0.98725, 0.98726, 0.98727, 0.98728, 0.98729, 0.9873, 0.98731, 0.98732, 0.98733, 0.98734, 0.98735, 0.98736, 0.98737, 0.98738, 0.98739, 0.9874, 0.98741, 0.98742, 0.98743, 0.98744, 0.98745, 0.98746, 0.98747, 0.98748, 0.98749, 0.9875, 0.98751, 0.98752, 0.98753, 0.98754, 0.98755, 0.98756, 0.98757, 0.98758, 0.98759, 0.9876, 0.98761, 0.98762, 0.98763, 0.98764, 0.98765, 0.98766, 0.98767, 0.98768, 0.98769, 0.9877, 0.98771, 0.98772, 0.98773, 0.98774, 0.98775, 0.98776, 0.98777, 0.98778, 0.98779, 0.9878, 0.98781, 0.98782, 0.98783, 0.98784, 0.98785, 0.98786, 0.98787, 0.98788, 0.98789, 0.9879, 0.98791, 0.98792, 0.98793, 0.98794, 0.98795, 0.98796, 0.98797, 0.98798, 0.98799, 0.988, 0.98801, 0.98802, 0.98803, 0.98804, 0.98805, 0.98806, 0.98807, 0.98808, 0.98809, 0.9881, 0.98811, 0.98812, 0.98813, 0.98814, 0.98815, 0.98816, 0.98817, 0.98818, 0.98819, 0.9882, 0.98821, 0.98822, 0.98823, 0.98824, 0.98825, 0.98826, 0.98827, 0.98828, 0.98829, 0.9883, 0.98831, 0.98832, 0.98833, 0.98834, 0.98835, 0.98836, 0.98837, 0.98838, 0.98839, 0.9884, 0.98841, 0.98842, 0.98843, 0.98844, 0.98845, 0.98846, 0.98847, 0.98848, 0.98849, 0.98850, 0.98851, 0.98852, 0.98853, 0.98854, 0.98855, 0.98856, 0.98857, 0.98858, 0.98859, 0.9886, 0.98861, 0.98862, 0.98863, 0.98864, 0.98865, 0.98866, 0.98867, 0.98868, 0.98869, 0.9887, 0.98871, 0.98872, 0.98873, 0.98874, 0.98875, 0.98876, 0.98877, 0.98878, 0.98879, 0.9888, 0.98881, 0.98882, 0.98883, 0.98884, 0.98885, 0.98886, 0.98887, 0.98888, 0.98889, 0.9889, 0.98891, 0.98892, 0.98893, 0.98894, 0.98895, 0.98896, 0.98897, 0.98898, 0.98899, 0.989, 0.98901, 0.98902, 0.98903, 0.98904, 0.98905, 0.98906, 0.98907, 0.98908, 0.98909, 0.9891, 0.98911, 0.98912, 0.98913, 0.98914, 0.98915, 0.98916, 0.98917, 0.98918, 0.98919, 0.9892, 0.98921, 0.98922, 0.98923, 0.98924, 0.98925, 0.98926, 0.98927, 0.98928, 0.98929, 0.9893, 0.98931, 0.98932, 0.98933, 0.98934, 0.98935, 0.98936, 0.98937, 0.98938, 0.98939, 0.9894, 0.98941, 0.98942, 0.98943, 0.98944, 0.98945, 0.98946, 0.98947, 0.98948, 0.98949, 0.9895, 0.98951, 0.98952, 0.98953, 0.98954, 0.98955, 0.98956, 0.98957, 0.98958, 0.98959, 0.9896, 0.98961, 0.98962, 0.98963, 0.98964, 0.98965, 0.98966, 0.98967, 0.98968, 0.98969, 0.9897, 0.98971, 0.98972, 0.98973, 0.98974, 0.98975, 0.98976, 0.98977, 0.98978, 0.98979, 0.9898, 0.98981, 0.98982, 0.98983, 0.98984, 0.98985, 0.98986, 0.98987, 0.98988, 0.98989, 0.9899, 0.98991, 0.98992, 0.98993, 0.98994, 0.98995, 0.98996, 0.98997, 0.98998, 0.98999, 0.99, 0.99001, 0.99002, 0.99003, 0.99004, 0.99005, 0.99006, 0.99007, 0.99008, 0.99009, 0.9901, 0.99011, 0.99012, 0.99013, 0.99014, 0.99015, 0.99016, 0.99017, 0.99018, 0.99019, 0.9902, 0.99021, 0.99022, 0.99023, 0.99024, 0.99025, 0.99026, 0.99027, 0.99028, 0.99029, 0.9903, 0.99031, 0.99032, 0.99033, 0.99034, 0.99035, 0.99036, 0.99037, 0.99038, 0.99039, 0.9904, 0.99041, 0.99042, 0.99043, 0.99044, 0.99045, 0.99046, 0.99047, 0.99048, 0.99049, 0.99050, 0.99051, 0.99052, 0.99053, 0.99054, 0.99055, 0.99056, 0.99057, 0.99058, 0.99059, 0.9906, 0.99061, 0.99062, 0.99063, 0.99064, 0.99065, 0.99066, 0.99067, 0.99068, 0.99069, 0.9907, 0.99071, 0.99072, 0.99073, 0.99074, 0.99075, 0.99076, 0.99077, 0.99078, 0.99079, 0.9908, 0.99081, 0.99082, 0.99083, 0.99084, 0.99085, 0.99086, 0.99087, 0.99088, 0.99089, 0.9909, 0.99091, 0.99092, 0.99093, 0.99094, 0.99095, 0.99096, 0.99097, 0.99098, 0.99099, 0.991, 0.99101, 0.99102, 0.99103, 0.99104, 0.99105, 0.99106, 0.99107, 0.99108, 0.99109, 0.9911, 0.99111, 0.99112, 0.99113, 0.99114, 0.99115, 0.99116, 0.99117, 0.99118, 0.99119, 0.9912, 0.99121, 0.99122, 0.99123, 0.99124, 0.99125, 0.99126, 0.99127, 0.99128, 0.99129, 0.9913, 0.99131, 0.99132, 0.99133, 0.99134, 0.99135, 0.99136, 0.99137, 0.99138, 0.99139, 0.9914, 0.99141, 0.99142, 0.99143, 0.99144, 0.99145, 0.99146, 0.99147, 0.99148, 0.99149, 0.9915, 0.99151, 0.99152, 0.99153, 0.99154, 0.99155, 0.99156, 0.99157, 0.99158, 0.99159, 0.9916, 0.99161, 0.99162, 0.99163, 0.99164, 0.99165, 0.99166, 0.99167, 0.99168, 0.99169, 0.9917, 0.99171, 0.99172, 0.99173, 0.99174, 0.99175, 0.99176, 0.99177, 0.99178, 0.99179, 0.9918, 0.99181, 0.99182, 0.99183, 0.99184, 0.99185, 0.99186, 0.99187, 0.99188, 0.99189, 0.9919, 0.99191, 0.99192, 0.99193, 0.99194, 0.99195, 0.99196, 0.99197, 0.99198, 0.99199, 0.992, 0.99201, 0.99202, 0.99203, 0.99204, 0.99205, 0.99206, 0.99207, 0.99208, 0.99209, 0.9921, 0.99211, 0.99212, 0.99213, 0.99214, 0.99215, 0.99216, 0.99217, 0.99218, 0.99219, 0.9922, 0.99221, 0.99222, 0.99223, 0.99224, 0.99225, 0.99226, 0.99227, 0.99228, 0.99229, 0.9923, 0.99231, 0.99232, 0.99233, 0.99234, 0.99235, 0.99236, 0.99237, 0.99238, 0.99239, 0.9924, 0.99241, 0.99242, 0.99243, 0.99244, 0.99245, 0.99246, 0.99247, 0.99248, 0.99249, 0.99250, 0.99251, 0.99252, 0.99253, 0.99254, 0.99255, 0.99256, 0.99257, 0.99258, 0.99259, 0.9926, 0.99261, 0.99262, 0.99263, 0.99264, 0.99265, 0.99266, 0.99267, 0.99268, 0.99269, 0.9927, 0.99271, 0.99272, 0.99273, 0.99274, 0.99275, 0.99276, 0.99277, 0.99278, 0.99279, 0.9928, 0.99281, 0.99282, 0.99283, 0.99284, 0.99285, 0.99286, 0.99287, 0.99288, 0.99289, 0.9929, 0.99291, 0.99292, 0.99293, 0.99294, 0.99295, 0.99296, 0.99297, 0.99298, 0.99299, 0.993, 0.99301, 0.99302, 0.99303, 0.99304, 0.99305, 0.99306, 0.99307, 0.99308, 0.99309, 0.9931, 0.99311, 0.99312, 0.99313, 0.99314, 0.99315, 0.99316, 0.99317, 0.99318, 0.99319, 0.9932, 0.99321, 0.99322, 0.99323, 0.99324, 0.99325, 0.99326, 0.99327, 0.99328, 0.99329, 0.9933, 0.99331, 0.99332, 0.99333, 0.99334, 0.99335, 0.99336, 0.99337, 0.99338, 0.99339, 0.9934, 0.99341, 0.99342, 0.99343, 0.99344, 0.99345, 0.99346, 0.99347, 0.99348, 0.99349, 0.9935, 0.99351, 0.99352, 0.99353, 0.99354, 0.99355, 0.99356, 0.99357, 0.99358, 0.99359, 0.9936, 0.99361, 0.99362, 0.99363, 0.99364, 0.99365, 0.99366, 0.99367, 0.99368, 0.99369, 0.9937, 0.99371, 0.99372, 0.99373, 0.99374, 0.99375, 0.99376, 0.99377, 0.99378, 0.99379, 0.9938, 0.99381, 0.99382, 0.99383, 0.99384, 0.99385, 0.99386, 0.99387, 0.99388, 0.99389, 0.9939, 0.99391, 0.99392, 0.99393, 0.99394, 0.99395, 0.99396, 0.99397, 0.99398, 0.99399, 0.994, 0.99401, 0.99402, 0.99403, 0.99404, 0.99405, 0.99406, 0.99407, 0.99408, 0.99409, 0.9941, 0.99411, 0.99412, 0.99413, 0.99414, 0.99415, 0.99416, 0.99417, 0.99418, 0.99419, 0.9942, 0.99421, 0.99422, 0.99423, 0.99424, 0.99425, 0.99426, 0.99427, 0.99428, 0.99429, 0.9943, 0.99431, 0.99432, 0.99433, 0.99434, 0.99435, 0.99436, 0.99437, 0.99438, 0.99439, 0.9944, 0.99441, 0.99442, 0.99443, 0.99444, 0.99445, 0.99446, 0.99447, 0.99448, 0.99449, 0.99450, 0.99451, 0.99452, 0.99453, 0.99454, 0.99455, 0.99456, 0.99457, 0.99458, 0.99459, 0.9946, 0.99461, 0.99462, 0.99463, 0.99464, 0.99465, 0.99466, 0.99467, 0.99468, 0.99469, 0.9947, 0.99471, 0.99472, 0.99473, 0.99474, 0.99475, 0.99476, 0.99477, 0.99478, 0.99479, 0.9948, 0.99481, 0.99482, 0.99483, 0.99484, 0.99485, 0.99486, 0.99487, 0.99488, 0.99489, 0.9949, 0.99491, 0.99492, 0.99493, 0.99494, 0.99495, 0.99496, 0.99497, 0.99498, 0.99499, 0.995, 0.99501, 0.99502, 0.99503, 0.99504, 0.99505, 0.99506, 0.99507, 0.99508, 0.99509, 0.9951, 0.99511, 0.99512, 0.99513, 0.99514, 0.99515, 0.99516, 0.99517, 0.99518, 0.99519, 0.9952, 0.99521, 0.99522, 0.99523, 0.99524, 0.99525, 0.99526, 0.99527, 0.99528, 0.99529, 0.9953, 0.99531, 0.99532, 0.99533, 0.99534, 0.99535, 0.99536, 0.99537, 0.99538, 0.99539, 0.9954, 0.99541, 0.99542, 0.99543, 0.99544, 0.99545, 0.99546, 0.99547, 0.99548, 0.99549, 0.9955, 0.99551, 0.99552, 0.99553, 0.99554, 0.99555, 0.99556, 0.99557, 0.99558, 0.99559, 0.9956, 0.99561, 0.99562, 0.99563, 0.99564, 0.99565, 0.99566, 0.99567, 0.99568, 0.99569, 0.9957, 0.99571, 0.99572, 0.99573, 0.99574, 0.99575, 0.99576, 0.99577, 0.99578, 0.99579, 0.9958, 0.99581, 0.99582, 0.99583, 0.99584, 0.99585, 0.99586, 0.99587, 0.99588, 0.99589, 0.9959, 0.99591, 0.99592, 0.99593, 0.99594, 0.99595, 0.99596, 0.99597, 0.99598, 0.99599, 0.996, 0.99601, 0.99602, 0.99603, 0.99604, 0.99605, 0.99606, 0.99607, 0.99608, 0.99609, 0.9961, 0.99611, 0.99612, 0.99613, 0.99614, 0.99615, 0.99616, 0.99617, 0.99618, 0.99619, 0.9962, 0.99621, 0.99622, 0.99623, 0.99624, 0.99625, 0.99626, 0.99627, 0.99628, 0.99629, 0.9963, 0.99631, 0.99632, 0.99633, 0.99634, 0.99635, 0.99636, 0.99637, 0.99638, 0.99639, 0.9964, 0.99641, 0.99642, 0.99643, 0.99644, 0.99645, 0.99646, 0.99647, 0.99648, 0.99649, 0.99650, 0.99651, 0.99652, 0.99653, 0.99654, 0.99655, 0.99656, 0.99657, 0.99658, 0.99659, 0.9966, 0.99661, 0.99662, 0.99663, 0.99664, 0.99665, 0.99666, 0.99667, 0.99668, 0.99669, 0.9967, 0.99671, 0.99672, 0.99673, 0.99674, 0.99675, 0.99676, 0.99677, 0.99678, 0.99679, 0.9968, 0.99681, 0.99682, 0.99683, 0.99684, 0.99685, 0.99686, 0.99687, 0.99688, 0.99689, 0.9969, 0.99691, 0.99692, 0.99693, 0.99694, 0.99695, 0.99696, 0.99697, 0.99698, 0.99699, 0.997, 0.99701, 0.99702, 0.99703, 0.99704, 0.99705, 0.99706, 0.99707, 0.99708, 0.99709, 0.9971, 0.99711, 0.99712, 0.99713, 0.99714, 0.99715, 0.99716, 0.99717, 0.99718, 0.99719, 0.9972, 0.99721, 0.99722, 0.99723, 0.99724, 0.99725, 0.99726, 0.99727, 0.99728, 0.99729, 0.9973, 0.99731, 0.99732, 0.99733, 0.99734, 0.99735, 0.99736, 0.99737, 0.99738, 0.99739, 0.9974, 0.99741, 0.99742, 0.99743, 0.99744, 0.99745, 0.99746, 0.99747, 0.99748, 0.99749, 0.9975, 0.99751, 0.99752, 0.99753, 0.99754, 0.99755, 0.99756, 0.99757, 0.99758, 0.99759, 0.9976, 0.99761, 0.99762, 0.99763, 0.99764, 0.99765, 0.99766, 0.99767, 0.99768, 0.99769, 0.9977, 0.99771, 0.99772, 0.99773, 0.99774, 0.99775, 0.99776, 0.99777, 0.99778, 0.99779, 0.9978, 0.99781, 0.99782, 0.99783, 0.99784, 0.99785, 0.99786, 0.99787, 0.99788, 0.99789, 0.9979, 0.99791, 0.99792, 0.99793, 0.99794, 0.99795, 0.99796, 0.99797, 0.99798, 0.99799, 0.998, 0.99801, 0.99802, 0.99803, 0.99804, 0.99805, 0.99806, 0.99807, 0.99808, 0.99809, 0.9981, 0.99811, 0.99812, 0.99813, 0.99814, 0.99815, 0.99816, 0.99817, 0.99818, 0.99819, 0.9982, 0.99821, 0.99822, 0.99823, 0.99824, 0.99825, 0.99826, 0.99827, 0.99828, 0.99829, 0.9983, 0.99831, 0.99832, 0.99833, 0.99834, 0.99835, 0.99836, 0.99837, 0.99838, 0.99839, 0.9984, 0.99841, 0.99842, 0.99843, 0.99844, 0.99845, 0.99846, 0.99847, 0.99848, 0.99849, 0.99850, 0.99851, 0.99852, 0.99853, 0.99854, 0.99855, 0.99856, 0.99857, 0.99858, 0.99859, 0.9986, 0.99861, 0.99862, 0.99863, 0.99864, 0.99865, 0.99866, 0.99867, 0.99868, 0.99869, 0.9987, 0.99871, 0.99872, 0.99873, 0.99874, 0.99875, 0.99876, 0.99877, 0.99878, 0.99879, 0.9988, 0.99881, 0.99882, 0.99883, 0.99884, 0.99885, 0.99886, 0.99887, 0.99888, 0.99889, 0.9989, 0.99891, 0.99892, 0.99893, 0.99894, 0.99895, 0.99896, 0.99897, 0.99898, 0.99899, 0.999, 0.99901, 0.99902, 0.99903, 0.99904, 0.99905, 0.99906, 0.99907, 0.99908, 0.99909, 0.9991, 0.99911, 0.99912, 0.99913, 0.99914, 0.99915, 0.99916, 0.99917, 0.99918, 0.99919, 0.9992, 0.99921, 0.99922, 0.99923, 0.99924, 0.99925, 0.99926, 0.99927, 0.99928, 0.99929, 0.9993, 0.99931, 0.99932, 0.99933, 0.99934, 0.99935, 0.99936, 0.99937, 0.99938, 0.99939, 0.9994, 0.99941, 0.99942, 0.99943, 0.99944, 0.99945, 0.99946, 0.99947, 0.99948, 0.99949, 0.9995, 0.99951, 0.99952, 0.99953, 0.99954, 0.99955, 0.99956, 0.99957, 0.99958, 0.99959, 0.9996, 0.99961, 0.99962, 0.99963, 0.99964, 0.99965, 0.99966, 0.99967, 0.99968, 0.99969, 0.9997, 0.99971, 0.99972, 0.99973, 0.99974, 0.99975, 0.99976, 0.99977, 0.99978, 0.99979, 0.9998, 0.99981, 0.99982, 0.99983, 0.99984, 0.99985, 0.99986, 0.99987, 0.99988, 0.99989, 0.9999, 0.99991, 0.99992, 0.99993, 0.99994, 0.99995, 0.99996, 0.99997, 0.99998, 0.99999, 1.0};
    Float_t probcuts_temp[trials+1] = {0.97, 0.97001, 0.97002, 0.97003, 0.97004, 0.97005, 0.97006, 0.97007, 0.97008, 0.97009, 0.9701, 0.97011, 0.97012, 0.97013, 0.97014, 0.97015, 0.97016, 0.97017, 0.97018, 0.97019, 0.9702, 0.97021, 0.97022, 0.97023, 0.97024, 0.97025, 0.97026, 0.97027, 0.97028, 0.97029, 0.97030, 0.97031, 0.97032, 0.97033, 0.97034, 0.97035, 0.97036, 0.97037, 0.97038, 0.97039, 0.9704, 0.97041, 0.97042, 0.97043, 0.97044, 0.97045, 0.97046, 0.97047, 0.97048, 0.97049, 0.97050, 0.97051, 0.97052, 0.97053, 0.97054, 0.97055, 0.97056, 0.97057, 0.97058, 0.97059, 0.9706, 0.97061, 0.97062, 0.97063, 0.97064, 0.97065, 0.97066, 0.97067, 0.97068, 0.97069, 0.9707, 0.97071, 0.97072, 0.97073, 0.97074, 0.97075, 0.97076, 0.97077, 0.97078, 0.97079, 0.9708, 0.97081, 0.97082, 0.97083, 0.97084, 0.97085, 0.97086, 0.97087, 0.97088, 0.97089, 0.9709, 0.97091, 0.97092, 0.97093, 0.97094, 0.97095, 0.97096, 0.97097, 0.97098, 0.97099, 0.971, 0.97101, 0.97102, 0.97103, 0.97104, 0.97105, 0.97106, 0.97107, 0.97108, 0.97109, 0.9711, 0.97111, 0.97112, 0.97113, 0.97114, 0.97115, 0.97116, 0.97117, 0.97118, 0.97119, 0.9712, 0.97121, 0.97122, 0.97123, 0.97124, 0.97125, 0.97126, 0.97127, 0.97128, 0.97129, 0.97130, 0.97131, 0.97132, 0.97133, 0.97134, 0.97135, 0.97136, 0.97137, 0.97138, 0.97139, 0.9714, 0.97141, 0.97142, 0.97143, 0.97144, 0.97145, 0.97146, 0.97147, 0.97148, 0.97149, 0.9715, 0.97151, 0.97152, 0.97153, 0.97154, 0.97155, 0.97156, 0.97157, 0.97158, 0.97159, 0.9716, 0.97161, 0.97162, 0.97163, 0.97164, 0.97165, 0.97166, 0.97167, 0.97168, 0.97169, 0.9717, 0.97171, 0.97172, 0.97173, 0.97174, 0.97175, 0.97176, 0.97177, 0.97178, 0.97179, 0.9718, 0.97181, 0.97182, 0.97183, 0.97184, 0.97185, 0.97186, 0.97187, 0.97188, 0.97189, 0.9719, 0.97191, 0.97192, 0.97193, 0.97194, 0.97195, 0.97196, 0.97197, 0.97198, 0.97199, 0.972, 0.97201, 0.97202, 0.97203, 0.97204, 0.97205, 0.97206, 0.97207, 0.97208, 0.97209, 0.9721, 0.97211, 0.97212, 0.97213, 0.97214, 0.97215, 0.97216, 0.97217, 0.97218, 0.97219, 0.9722, 0.97221, 0.97222, 0.97223, 0.97224, 0.97225, 0.97226, 0.97227, 0.97228, 0.97229, 0.97230, 0.97231, 0.97232, 0.97233, 0.97234, 0.97235, 0.97236, 0.97237, 0.97238, 0.97239, 0.9724, 0.97241, 0.97242, 0.97243, 0.97244, 0.97245, 0.97246, 0.97247, 0.97248, 0.97249, 0.97250, 0.97251, 0.97252, 0.97253, 0.97254, 0.97255, 0.97256, 0.97257, 0.97258, 0.97259, 0.9726, 0.97261, 0.97262, 0.97263, 0.97264, 0.97265, 0.97266, 0.97267, 0.97268, 0.97269, 0.9727, 0.97271, 0.97272, 0.97273, 0.97274, 0.97275, 0.97276, 0.97277, 0.97278, 0.97279, 0.9728, 0.97281, 0.97282, 0.97283, 0.97284, 0.97285, 0.97286, 0.97287, 0.97288, 0.97289, 0.9729, 0.97291, 0.97292, 0.97293, 0.97294, 0.97295, 0.97296, 0.97297, 0.97298, 0.97299, 0.973, 0.97301, 0.97302, 0.97303, 0.97304, 0.97305, 0.97306, 0.97307, 0.97308, 0.97309, 0.9731, 0.97311, 0.97312, 0.97313, 0.97314, 0.97315, 0.97316, 0.97317, 0.97318, 0.97319, 0.9732, 0.97321, 0.97322, 0.97323, 0.97324, 0.97325, 0.97326, 0.97327, 0.97328, 0.97329, 0.97330, 0.97331, 0.97332, 0.97333, 0.97334, 0.97335, 0.97336, 0.97337, 0.97338, 0.97339, 0.9734, 0.97341, 0.97342, 0.97343, 0.97344, 0.97345, 0.97346, 0.97347, 0.97348, 0.97349, 0.9735, 0.97351, 0.97352, 0.97353, 0.97354, 0.97355, 0.97356, 0.97357, 0.97358, 0.97359, 0.9736, 0.97361, 0.97362, 0.97363, 0.97364, 0.97365, 0.97366, 0.97367, 0.97368, 0.97369, 0.9737, 0.97371, 0.97372, 0.97373, 0.97374, 0.97375, 0.97376, 0.97377, 0.97378, 0.97379, 0.9738, 0.97381, 0.97382, 0.97383, 0.97384, 0.97385, 0.97386, 0.97387, 0.97388, 0.97389, 0.9739, 0.97391, 0.97392, 0.97393, 0.97394, 0.97395, 0.97396, 0.97397, 0.97398, 0.97399, 0.974, 0.97401, 0.97402, 0.97403, 0.97404, 0.97405, 0.97406, 0.97407, 0.97408, 0.97409, 0.9741, 0.97411, 0.97412, 0.97413, 0.97414, 0.97415, 0.97416, 0.97417, 0.97418, 0.97419, 0.9742, 0.97421, 0.97422, 0.97423, 0.97424, 0.97425, 0.97426, 0.97427, 0.97428, 0.97429, 0.97430, 0.97431, 0.97432, 0.97433, 0.97434, 0.97435, 0.97436, 0.97437, 0.97438, 0.97439, 0.9744, 0.97441, 0.97442, 0.97443, 0.97444, 0.97445, 0.97446, 0.97447, 0.97448, 0.97449, 0.97450, 0.97451, 0.97452, 0.97453, 0.97454, 0.97455, 0.97456, 0.97457, 0.97458, 0.97459, 0.9746, 0.97461, 0.97462, 0.97463, 0.97464, 0.97465, 0.97466, 0.97467, 0.97468, 0.97469, 0.9747, 0.97471, 0.97472, 0.97473, 0.97474, 0.97475, 0.97476, 0.97477, 0.97478, 0.97479, 0.9748, 0.97481, 0.97482, 0.97483, 0.97484, 0.97485, 0.97486, 0.97487, 0.97488, 0.97489, 0.9749, 0.97491, 0.97492, 0.97493, 0.97494, 0.97495, 0.97496, 0.97497, 0.97498, 0.97499, 0.975, 0.97501, 0.97502, 0.97503, 0.97504, 0.97505, 0.97506, 0.97507, 0.97508, 0.97509, 0.9751, 0.97511, 0.97512, 0.97513, 0.97514, 0.97515, 0.97516, 0.97517, 0.97518, 0.97519, 0.9752, 0.97521, 0.97522, 0.97523, 0.97524, 0.97525, 0.97526, 0.97527, 0.97528, 0.97529, 0.9753, 0.97531, 0.97532, 0.97533, 0.97534, 0.97535, 0.97536, 0.97537, 0.97538, 0.97539, 0.9754, 0.97541, 0.97542, 0.97543, 0.97544, 0.97545, 0.97546, 0.97547, 0.97548, 0.97549, 0.9755, 0.97551, 0.97552, 0.97553, 0.97554, 0.97555, 0.97556, 0.97557, 0.97558, 0.97559, 0.9756, 0.97561, 0.97562, 0.97563, 0.97564, 0.97565, 0.97566, 0.97567, 0.97568, 0.97569, 0.9757, 0.97571, 0.97572, 0.97573, 0.97574, 0.97575, 0.97576, 0.97577, 0.97578, 0.97579, 0.9758, 0.97581, 0.97582, 0.97583, 0.97584, 0.97585, 0.97586, 0.97587, 0.97588, 0.97589, 0.9759, 0.97591, 0.97592, 0.97593, 0.97594, 0.97595, 0.97596, 0.97597, 0.97598, 0.97599, 0.976, 0.97601, 0.97602, 0.97603, 0.97604, 0.97605, 0.97606, 0.97607, 0.97608, 0.97609, 0.9761, 0.97611, 0.97612, 0.97613, 0.97614, 0.97615, 0.97616, 0.97617, 0.97618, 0.97619, 0.9762, 0.97621, 0.97622, 0.97623, 0.97624, 0.97625, 0.97626, 0.97627, 0.97628, 0.97629, 0.9763, 0.97631, 0.97632, 0.97633, 0.97634, 0.97635, 0.97636, 0.97637, 0.97638, 0.97639, 0.9764, 0.97641, 0.97642, 0.97643, 0.97644, 0.97645, 0.97646, 0.97647, 0.97648, 0.97649, 0.97650, 0.97651, 0.97652, 0.97653, 0.97654, 0.97655, 0.97656, 0.97657, 0.97658, 0.97659, 0.9766, 0.97661, 0.97662, 0.97663, 0.97664, 0.97665, 0.97666, 0.97667, 0.97668, 0.97669, 0.9767, 0.97671, 0.97672, 0.97673, 0.97674, 0.97675, 0.97676, 0.97677, 0.97678, 0.97679, 0.9768, 0.97681, 0.97682, 0.97683, 0.97684, 0.97685, 0.97686, 0.97687, 0.97688, 0.97689, 0.9769, 0.97691, 0.97692, 0.97693, 0.97694, 0.97695, 0.97696, 0.97697, 0.97698, 0.97699, 0.977, 0.97701, 0.97702, 0.97703, 0.97704, 0.97705, 0.97706, 0.97707, 0.97708, 0.97709, 0.9771, 0.97711, 0.97712, 0.97713, 0.97714, 0.97715, 0.97716, 0.97717, 0.97718, 0.97719, 0.9772, 0.97721, 0.97722, 0.97723, 0.97724, 0.97725, 0.97726, 0.97727, 0.97728, 0.97729, 0.9773, 0.97731, 0.97732, 0.97733, 0.97734, 0.97735, 0.97736, 0.97737, 0.97738, 0.97739, 0.9774, 0.97741, 0.97742, 0.97743, 0.97744, 0.97745, 0.97746, 0.97747, 0.97748, 0.97749, 0.9775, 0.97751, 0.97752, 0.97753, 0.97754, 0.97755, 0.97756, 0.97757, 0.97758, 0.97759, 0.9776, 0.97761, 0.97762, 0.97763, 0.97764, 0.97765, 0.97766, 0.97767, 0.97768, 0.97769, 0.9777, 0.97771, 0.97772, 0.97773, 0.97774, 0.97775, 0.97776, 0.97777, 0.97778, 0.97779, 0.9778, 0.97781, 0.97782, 0.97783, 0.97784, 0.97785, 0.97786, 0.97787, 0.97788, 0.97789, 0.9779, 0.97791, 0.97792, 0.97793, 0.97794, 0.97795, 0.97796, 0.97797, 0.97798, 0.97799, 0.978, 0.97801, 0.97802, 0.97803, 0.97804, 0.97805, 0.97806, 0.97807, 0.97808, 0.97809, 0.9781, 0.97811, 0.97812, 0.97813, 0.97814, 0.97815, 0.97816, 0.97817, 0.97818, 0.97819, 0.9782, 0.97821, 0.97822, 0.97823, 0.97824, 0.97825, 0.97826, 0.97827, 0.97828, 0.97829, 0.9783, 0.97831, 0.97832, 0.97833, 0.97834, 0.97835, 0.97836, 0.97837, 0.97838, 0.97839, 0.9784, 0.97841, 0.97842, 0.97843, 0.97844, 0.97845, 0.97846, 0.97847, 0.97848, 0.97849, 0.97850, 0.97851, 0.97852, 0.97853, 0.97854, 0.97855, 0.97856, 0.97857, 0.97858, 0.97859, 0.9786, 0.97861, 0.97862, 0.97863, 0.97864, 0.97865, 0.97866, 0.97867, 0.97868, 0.97869, 0.9787, 0.97871, 0.97872, 0.97873, 0.97874, 0.97875, 0.97876, 0.97877, 0.97878, 0.97879, 0.9788, 0.97881, 0.97882, 0.97883, 0.97884, 0.97885, 0.97886, 0.97887, 0.97888, 0.97889, 0.9789, 0.97891, 0.97892, 0.97893, 0.97894, 0.97895, 0.97896, 0.97897, 0.97898, 0.97899, 0.979, 0.97901, 0.97902, 0.97903, 0.97904, 0.97905, 0.97906, 0.97907, 0.97908, 0.97909, 0.9791, 0.97911, 0.97912, 0.97913, 0.97914, 0.97915, 0.97916, 0.97917, 0.97918, 0.97919, 0.9792, 0.97921, 0.97922, 0.97923, 0.97924, 0.97925, 0.97926, 0.97927, 0.97928, 0.97929, 0.9793, 0.97931, 0.97932, 0.97933, 0.97934, 0.97935, 0.97936, 0.97937, 0.97938, 0.97939, 0.9794, 0.97941, 0.97942, 0.97943, 0.97944, 0.97945, 0.97946, 0.97947, 0.97948, 0.97949, 0.9795, 0.97951, 0.97952, 0.97953, 0.97954, 0.97955, 0.97956, 0.97957, 0.97958, 0.97959, 0.9796, 0.97961, 0.97962, 0.97963, 0.97964, 0.97965, 0.97966, 0.97967, 0.97968, 0.97969, 0.9797, 0.97971, 0.97972, 0.97973, 0.97974, 0.97975, 0.97976, 0.97977, 0.97978, 0.97979, 0.9798, 0.97981, 0.97982, 0.97983, 0.97984, 0.97985, 0.97986, 0.97987, 0.97988, 0.97989, 0.9799, 0.97991, 0.97992, 0.97993, 0.97994, 0.97995, 0.97996, 0.97997, 0.97998, 0.97999, 0.98, 0.98001, 0.98002, 0.98003, 0.98004, 0.98005, 0.98006, 0.98007, 0.98008, 0.98009, 0.9801, 0.98011, 0.98012, 0.98013, 0.98014, 0.98015, 0.98016, 0.98017, 0.98018, 0.98019, 0.9802, 0.98021, 0.98022, 0.98023, 0.98024, 0.98025, 0.98026, 0.98027, 0.98028, 0.98029, 0.9803, 0.98031, 0.98032, 0.98033, 0.98034, 0.98035, 0.98036, 0.98037, 0.98038, 0.98039, 0.9804, 0.98041, 0.98042, 0.98043, 0.98044, 0.98045, 0.98046, 0.98047, 0.98048, 0.98049, 0.98050, 0.98051, 0.98052, 0.98053, 0.98054, 0.98055, 0.98056, 0.98057, 0.98058, 0.98059, 0.9806, 0.98061, 0.98062, 0.98063, 0.98064, 0.98065, 0.98066, 0.98067, 0.98068, 0.98069, 0.9807, 0.98071, 0.98072, 0.98073, 0.98074, 0.98075, 0.98076, 0.98077, 0.98078, 0.98079, 0.9808, 0.98081, 0.98082, 0.98083, 0.98084, 0.98085, 0.98086, 0.98087, 0.98088, 0.98089, 0.9809, 0.98091, 0.98092, 0.98093, 0.98094, 0.98095, 0.98096, 0.98097, 0.98098, 0.98099, 0.981, 0.98101, 0.98102, 0.98103, 0.98104, 0.98105, 0.98106, 0.98107, 0.98108, 0.98109, 0.9811, 0.98111, 0.98112, 0.98113, 0.98114, 0.98115, 0.98116, 0.98117, 0.98118, 0.98119, 0.9812, 0.98121, 0.98122, 0.98123, 0.98124, 0.98125, 0.98126, 0.98127, 0.98128, 0.98129, 0.9813, 0.98131, 0.98132, 0.98133, 0.98134, 0.98135, 0.98136, 0.98137, 0.98138, 0.98139, 0.9814, 0.98141, 0.98142, 0.98143, 0.98144, 0.98145, 0.98146, 0.98147, 0.98148, 0.98149, 0.9815, 0.98151, 0.98152, 0.98153, 0.98154, 0.98155, 0.98156, 0.98157, 0.98158, 0.98159, 0.9816, 0.98161, 0.98162, 0.98163, 0.98164, 0.98165, 0.98166, 0.98167, 0.98168, 0.98169, 0.9817, 0.98171, 0.98172, 0.98173, 0.98174, 0.98175, 0.98176, 0.98177, 0.98178, 0.98179, 0.9818, 0.98181, 0.98182, 0.98183, 0.98184, 0.98185, 0.98186, 0.98187, 0.98188, 0.98189, 0.9819, 0.98191, 0.98192, 0.98193, 0.98194, 0.98195, 0.98196, 0.98197, 0.98198, 0.98199, 0.982, 0.98201, 0.98202, 0.98203, 0.98204, 0.98205, 0.98206, 0.98207, 0.98208, 0.98209, 0.9821, 0.98211, 0.98212, 0.98213, 0.98214, 0.98215, 0.98216, 0.98217, 0.98218, 0.98219, 0.9822, 0.98221, 0.98222, 0.98223, 0.98224, 0.98225, 0.98226, 0.98227, 0.98228, 0.98229, 0.9823, 0.98231, 0.98232, 0.98233, 0.98234, 0.98235, 0.98236, 0.98237, 0.98238, 0.98239, 0.9824, 0.98241, 0.98242, 0.98243, 0.98244, 0.98245, 0.98246, 0.98247, 0.98248, 0.98249, 0.98250, 0.98251, 0.98252, 0.98253, 0.98254, 0.98255, 0.98256, 0.98257, 0.98258, 0.98259, 0.9826, 0.98261, 0.98262, 0.98263, 0.98264, 0.98265, 0.98266, 0.98267, 0.98268, 0.98269, 0.9827, 0.98271, 0.98272, 0.98273, 0.98274, 0.98275, 0.98276, 0.98277, 0.98278, 0.98279, 0.9828, 0.98281, 0.98282, 0.98283, 0.98284, 0.98285, 0.98286, 0.98287, 0.98288, 0.98289, 0.9829, 0.98291, 0.98292, 0.98293, 0.98294, 0.98295, 0.98296, 0.98297, 0.98298, 0.98299, 0.983, 0.98301, 0.98302, 0.98303, 0.98304, 0.98305, 0.98306, 0.98307, 0.98308, 0.98309, 0.9831, 0.98311, 0.98312, 0.98313, 0.98314, 0.98315, 0.98316, 0.98317, 0.98318, 0.98319, 0.9832, 0.98321, 0.98322, 0.98323, 0.98324, 0.98325, 0.98326, 0.98327, 0.98328, 0.98329, 0.9833, 0.98331, 0.98332, 0.98333, 0.98334, 0.98335, 0.98336, 0.98337, 0.98338, 0.98339, 0.9834, 0.98341, 0.98342, 0.98343, 0.98344, 0.98345, 0.98346, 0.98347, 0.98348, 0.98349, 0.9835, 0.98351, 0.98352, 0.98353, 0.98354, 0.98355, 0.98356, 0.98357, 0.98358, 0.98359, 0.9836, 0.98361, 0.98362, 0.98363, 0.98364, 0.98365, 0.98366, 0.98367, 0.98368, 0.98369, 0.9837, 0.98371, 0.98372, 0.98373, 0.98374, 0.98375, 0.98376, 0.98377, 0.98378, 0.98379, 0.9838, 0.98381, 0.98382, 0.98383, 0.98384, 0.98385, 0.98386, 0.98387, 0.98388, 0.98389, 0.9839, 0.98391, 0.98392, 0.98393, 0.98394, 0.98395, 0.98396, 0.98397, 0.98398, 0.98399, 0.984, 0.98401, 0.98402, 0.98403, 0.98404, 0.98405, 0.98406, 0.98407, 0.98408, 0.98409, 0.9841, 0.98411, 0.98412, 0.98413, 0.98414, 0.98415, 0.98416, 0.98417, 0.98418, 0.98419, 0.9842, 0.98421, 0.98422, 0.98423, 0.98424, 0.98425, 0.98426, 0.98427, 0.98428, 0.98429, 0.9843, 0.98431, 0.98432, 0.98433, 0.98434, 0.98435, 0.98436, 0.98437, 0.98438, 0.98439, 0.9844, 0.98441, 0.98442, 0.98443, 0.98444, 0.98445, 0.98446, 0.98447, 0.98448, 0.98449, 0.98450, 0.98451, 0.98452, 0.98453, 0.98454, 0.98455, 0.98456, 0.98457, 0.98458, 0.98459, 0.9846, 0.98461, 0.98462, 0.98463, 0.98464, 0.98465, 0.98466, 0.98467, 0.98468, 0.98469, 0.9847, 0.98471, 0.98472, 0.98473, 0.98474, 0.98475, 0.98476, 0.98477, 0.98478, 0.98479, 0.9848, 0.98481, 0.98482, 0.98483, 0.98484, 0.98485, 0.98486, 0.98487, 0.98488, 0.98489, 0.9849, 0.98491, 0.98492, 0.98493, 0.98494, 0.98495, 0.98496, 0.98497, 0.98498, 0.98499, 0.985, 0.98501, 0.98502, 0.98503, 0.98504, 0.98505, 0.98506, 0.98507, 0.98508, 0.98509, 0.9851, 0.98511, 0.98512, 0.98513, 0.98514, 0.98515, 0.98516, 0.98517, 0.98518, 0.98519, 0.9852, 0.98521, 0.98522, 0.98523, 0.98524, 0.98525, 0.98526, 0.98527, 0.98528, 0.98529, 0.9853, 0.98531, 0.98532, 0.98533, 0.98534, 0.98535, 0.98536, 0.98537, 0.98538, 0.98539, 0.9854, 0.98541, 0.98542, 0.98543, 0.98544, 0.98545, 0.98546, 0.98547, 0.98548, 0.98549, 0.9855, 0.98551, 0.98552, 0.98553, 0.98554, 0.98555, 0.98556, 0.98557, 0.98558, 0.98559, 0.9856, 0.98561, 0.98562, 0.98563, 0.98564, 0.98565, 0.98566, 0.98567, 0.98568, 0.98569, 0.9857, 0.98571, 0.98572, 0.98573, 0.98574, 0.98575, 0.98576, 0.98577, 0.98578, 0.98579, 0.9858, 0.98581, 0.98582, 0.98583, 0.98584, 0.98585, 0.98586, 0.98587, 0.98588, 0.98589, 0.9859, 0.98591, 0.98592, 0.98593, 0.98594, 0.98595, 0.98596, 0.98597, 0.98598, 0.98599, 0.986, 0.98601, 0.98602, 0.98603, 0.98604, 0.98605, 0.98606, 0.98607, 0.98608, 0.98609, 0.9861, 0.98611, 0.98612, 0.98613, 0.98614, 0.98615, 0.98616, 0.98617, 0.98618, 0.98619, 0.9862, 0.98621, 0.98622, 0.98623, 0.98624, 0.98625, 0.98626, 0.98627, 0.98628, 0.98629, 0.9863, 0.98631, 0.98632, 0.98633, 0.98634, 0.98635, 0.98636, 0.98637, 0.98638, 0.98639, 0.9864, 0.98641, 0.98642, 0.98643, 0.98644, 0.98645, 0.98646, 0.98647, 0.98648, 0.98649, 0.98650, 0.98651, 0.98652, 0.98653, 0.98654, 0.98655, 0.98656, 0.98657, 0.98658, 0.98659, 0.9866, 0.98661, 0.98662, 0.98663, 0.98664, 0.98665, 0.98666, 0.98667, 0.98668, 0.98669, 0.9867, 0.98671, 0.98672, 0.98673, 0.98674, 0.98675, 0.98676, 0.98677, 0.98678, 0.98679, 0.9868, 0.98681, 0.98682, 0.98683, 0.98684, 0.98685, 0.98686, 0.98687, 0.98688, 0.98689, 0.9869, 0.98691, 0.98692, 0.98693, 0.98694, 0.98695, 0.98696, 0.98697, 0.98698, 0.98699, 0.987, 0.98701, 0.98702, 0.98703, 0.98704, 0.98705, 0.98706, 0.98707, 0.98708, 0.98709, 0.9871, 0.98711, 0.98712, 0.98713, 0.98714, 0.98715, 0.98716, 0.98717, 0.98718, 0.98719, 0.9872, 0.98721, 0.98722, 0.98723, 0.98724, 0.98725, 0.98726, 0.98727, 0.98728, 0.98729, 0.9873, 0.98731, 0.98732, 0.98733, 0.98734, 0.98735, 0.98736, 0.98737, 0.98738, 0.98739, 0.9874, 0.98741, 0.98742, 0.98743, 0.98744, 0.98745, 0.98746, 0.98747, 0.98748, 0.98749, 0.9875, 0.98751, 0.98752, 0.98753, 0.98754, 0.98755, 0.98756, 0.98757, 0.98758, 0.98759, 0.9876, 0.98761, 0.98762, 0.98763, 0.98764, 0.98765, 0.98766, 0.98767, 0.98768, 0.98769, 0.9877, 0.98771, 0.98772, 0.98773, 0.98774, 0.98775, 0.98776, 0.98777, 0.98778, 0.98779, 0.9878, 0.98781, 0.98782, 0.98783, 0.98784, 0.98785, 0.98786, 0.98787, 0.98788, 0.98789, 0.9879, 0.98791, 0.98792, 0.98793, 0.98794, 0.98795, 0.98796, 0.98797, 0.98798, 0.98799, 0.988, 0.98801, 0.98802, 0.98803, 0.98804, 0.98805, 0.98806, 0.98807, 0.98808, 0.98809, 0.9881, 0.98811, 0.98812, 0.98813, 0.98814, 0.98815, 0.98816, 0.98817, 0.98818, 0.98819, 0.9882, 0.98821, 0.98822, 0.98823, 0.98824, 0.98825, 0.98826, 0.98827, 0.98828, 0.98829, 0.9883, 0.98831, 0.98832, 0.98833, 0.98834, 0.98835, 0.98836, 0.98837, 0.98838, 0.98839, 0.9884, 0.98841, 0.98842, 0.98843, 0.98844, 0.98845, 0.98846, 0.98847, 0.98848, 0.98849, 0.98850, 0.98851, 0.98852, 0.98853, 0.98854, 0.98855, 0.98856, 0.98857, 0.98858, 0.98859, 0.9886, 0.98861, 0.98862, 0.98863, 0.98864, 0.98865, 0.98866, 0.98867, 0.98868, 0.98869, 0.9887, 0.98871, 0.98872, 0.98873, 0.98874, 0.98875, 0.98876, 0.98877, 0.98878, 0.98879, 0.9888, 0.98881, 0.98882, 0.98883, 0.98884, 0.98885, 0.98886, 0.98887, 0.98888, 0.98889, 0.9889, 0.98891, 0.98892, 0.98893, 0.98894, 0.98895, 0.98896, 0.98897, 0.98898, 0.98899, 0.989, 0.98901, 0.98902, 0.98903, 0.98904, 0.98905, 0.98906, 0.98907, 0.98908, 0.98909, 0.9891, 0.98911, 0.98912, 0.98913, 0.98914, 0.98915, 0.98916, 0.98917, 0.98918, 0.98919, 0.9892, 0.98921, 0.98922, 0.98923, 0.98924, 0.98925, 0.98926, 0.98927, 0.98928, 0.98929, 0.9893, 0.98931, 0.98932, 0.98933, 0.98934, 0.98935, 0.98936, 0.98937, 0.98938, 0.98939, 0.9894, 0.98941, 0.98942, 0.98943, 0.98944, 0.98945, 0.98946, 0.98947, 0.98948, 0.98949, 0.9895, 0.98951, 0.98952, 0.98953, 0.98954, 0.98955, 0.98956, 0.98957, 0.98958, 0.98959, 0.9896, 0.98961, 0.98962, 0.98963, 0.98964, 0.98965, 0.98966, 0.98967, 0.98968, 0.98969, 0.9897, 0.98971, 0.98972, 0.98973, 0.98974, 0.98975, 0.98976, 0.98977, 0.98978, 0.98979, 0.9898, 0.98981, 0.98982, 0.98983, 0.98984, 0.98985, 0.98986, 0.98987, 0.98988, 0.98989, 0.9899, 0.98991, 0.98992, 0.98993, 0.98994, 0.98995, 0.98996, 0.98997, 0.98998, 0.98999, 0.99, 0.99001, 0.99002, 0.99003, 0.99004, 0.99005, 0.99006, 0.99007, 0.99008, 0.99009, 0.9901, 0.99011, 0.99012, 0.99013, 0.99014, 0.99015, 0.99016, 0.99017, 0.99018, 0.99019, 0.9902, 0.99021, 0.99022, 0.99023, 0.99024, 0.99025, 0.99026, 0.99027, 0.99028, 0.99029, 0.9903, 0.99031, 0.99032, 0.99033, 0.99034, 0.99035, 0.99036, 0.99037, 0.99038, 0.99039, 0.9904, 0.99041, 0.99042, 0.99043, 0.99044, 0.99045, 0.99046, 0.99047, 0.99048, 0.99049, 0.99050, 0.99051, 0.99052, 0.99053, 0.99054, 0.99055, 0.99056, 0.99057, 0.99058, 0.99059, 0.9906, 0.99061, 0.99062, 0.99063, 0.99064, 0.99065, 0.99066, 0.99067, 0.99068, 0.99069, 0.9907, 0.99071, 0.99072, 0.99073, 0.99074, 0.99075, 0.99076, 0.99077, 0.99078, 0.99079, 0.9908, 0.99081, 0.99082, 0.99083, 0.99084, 0.99085, 0.99086, 0.99087, 0.99088, 0.99089, 0.9909, 0.99091, 0.99092, 0.99093, 0.99094, 0.99095, 0.99096, 0.99097, 0.99098, 0.99099, 0.991, 0.99101, 0.99102, 0.99103, 0.99104, 0.99105, 0.99106, 0.99107, 0.99108, 0.99109, 0.9911, 0.99111, 0.99112, 0.99113, 0.99114, 0.99115, 0.99116, 0.99117, 0.99118, 0.99119, 0.9912, 0.99121, 0.99122, 0.99123, 0.99124, 0.99125, 0.99126, 0.99127, 0.99128, 0.99129, 0.9913, 0.99131, 0.99132, 0.99133, 0.99134, 0.99135, 0.99136, 0.99137, 0.99138, 0.99139, 0.9914, 0.99141, 0.99142, 0.99143, 0.99144, 0.99145, 0.99146, 0.99147, 0.99148, 0.99149, 0.9915, 0.99151, 0.99152, 0.99153, 0.99154, 0.99155, 0.99156, 0.99157, 0.99158, 0.99159, 0.9916, 0.99161, 0.99162, 0.99163, 0.99164, 0.99165, 0.99166, 0.99167, 0.99168, 0.99169, 0.9917, 0.99171, 0.99172, 0.99173, 0.99174, 0.99175, 0.99176, 0.99177, 0.99178, 0.99179, 0.9918, 0.99181, 0.99182, 0.99183, 0.99184, 0.99185, 0.99186, 0.99187, 0.99188, 0.99189, 0.9919, 0.99191, 0.99192, 0.99193, 0.99194, 0.99195, 0.99196, 0.99197, 0.99198, 0.99199, 0.992, 0.99201, 0.99202, 0.99203, 0.99204, 0.99205, 0.99206, 0.99207, 0.99208, 0.99209, 0.9921, 0.99211, 0.99212, 0.99213, 0.99214, 0.99215, 0.99216, 0.99217, 0.99218, 0.99219, 0.9922, 0.99221, 0.99222, 0.99223, 0.99224, 0.99225, 0.99226, 0.99227, 0.99228, 0.99229, 0.9923, 0.99231, 0.99232, 0.99233, 0.99234, 0.99235, 0.99236, 0.99237, 0.99238, 0.99239, 0.9924, 0.99241, 0.99242, 0.99243, 0.99244, 0.99245, 0.99246, 0.99247, 0.99248, 0.99249, 0.99250, 0.99251, 0.99252, 0.99253, 0.99254, 0.99255, 0.99256, 0.99257, 0.99258, 0.99259, 0.9926, 0.99261, 0.99262, 0.99263, 0.99264, 0.99265, 0.99266, 0.99267, 0.99268, 0.99269, 0.9927, 0.99271, 0.99272, 0.99273, 0.99274, 0.99275, 0.99276, 0.99277, 0.99278, 0.99279, 0.9928, 0.99281, 0.99282, 0.99283, 0.99284, 0.99285, 0.99286, 0.99287, 0.99288, 0.99289, 0.9929, 0.99291, 0.99292, 0.99293, 0.99294, 0.99295, 0.99296, 0.99297, 0.99298, 0.99299, 0.993, 0.99301, 0.99302, 0.99303, 0.99304, 0.99305, 0.99306, 0.99307, 0.99308, 0.99309, 0.9931, 0.99311, 0.99312, 0.99313, 0.99314, 0.99315, 0.99316, 0.99317, 0.99318, 0.99319, 0.9932, 0.99321, 0.99322, 0.99323, 0.99324, 0.99325, 0.99326, 0.99327, 0.99328, 0.99329, 0.9933, 0.99331, 0.99332, 0.99333, 0.99334, 0.99335, 0.99336, 0.99337, 0.99338, 0.99339, 0.9934, 0.99341, 0.99342, 0.99343, 0.99344, 0.99345, 0.99346, 0.99347, 0.99348, 0.99349, 0.9935, 0.99351, 0.99352, 0.99353, 0.99354, 0.99355, 0.99356, 0.99357, 0.99358, 0.99359, 0.9936, 0.99361, 0.99362, 0.99363, 0.99364, 0.99365, 0.99366, 0.99367, 0.99368, 0.99369, 0.9937, 0.99371, 0.99372, 0.99373, 0.99374, 0.99375, 0.99376, 0.99377, 0.99378, 0.99379, 0.9938, 0.99381, 0.99382, 0.99383, 0.99384, 0.99385, 0.99386, 0.99387, 0.99388, 0.99389, 0.9939, 0.99391, 0.99392, 0.99393, 0.99394, 0.99395, 0.99396, 0.99397, 0.99398, 0.99399, 0.994, 0.99401, 0.99402, 0.99403, 0.99404, 0.99405, 0.99406, 0.99407, 0.99408, 0.99409, 0.9941, 0.99411, 0.99412, 0.99413, 0.99414, 0.99415, 0.99416, 0.99417, 0.99418, 0.99419, 0.9942, 0.99421, 0.99422, 0.99423, 0.99424, 0.99425, 0.99426, 0.99427, 0.99428, 0.99429, 0.9943, 0.99431, 0.99432, 0.99433, 0.99434, 0.99435, 0.99436, 0.99437, 0.99438, 0.99439, 0.9944, 0.99441, 0.99442, 0.99443, 0.99444, 0.99445, 0.99446, 0.99447, 0.99448, 0.99449, 0.99450, 0.99451, 0.99452, 0.99453, 0.99454, 0.99455, 0.99456, 0.99457, 0.99458, 0.99459, 0.9946, 0.99461, 0.99462, 0.99463, 0.99464, 0.99465, 0.99466, 0.99467, 0.99468, 0.99469, 0.9947, 0.99471, 0.99472, 0.99473, 0.99474, 0.99475, 0.99476, 0.99477, 0.99478, 0.99479, 0.9948, 0.99481, 0.99482, 0.99483, 0.99484, 0.99485, 0.99486, 0.99487, 0.99488, 0.99489, 0.9949, 0.99491, 0.99492, 0.99493, 0.99494, 0.99495, 0.99496, 0.99497, 0.99498, 0.99499, 0.995, 0.99501, 0.99502, 0.99503, 0.99504, 0.99505, 0.99506, 0.99507, 0.99508, 0.99509, 0.9951, 0.99511, 0.99512, 0.99513, 0.99514, 0.99515, 0.99516, 0.99517, 0.99518, 0.99519, 0.9952, 0.99521, 0.99522, 0.99523, 0.99524, 0.99525, 0.99526, 0.99527, 0.99528, 0.99529, 0.9953, 0.99531, 0.99532, 0.99533, 0.99534, 0.99535, 0.99536, 0.99537, 0.99538, 0.99539, 0.9954, 0.99541, 0.99542, 0.99543, 0.99544, 0.99545, 0.99546, 0.99547, 0.99548, 0.99549, 0.9955, 0.99551, 0.99552, 0.99553, 0.99554, 0.99555, 0.99556, 0.99557, 0.99558, 0.99559, 0.9956, 0.99561, 0.99562, 0.99563, 0.99564, 0.99565, 0.99566, 0.99567, 0.99568, 0.99569, 0.9957, 0.99571, 0.99572, 0.99573, 0.99574, 0.99575, 0.99576, 0.99577, 0.99578, 0.99579, 0.9958, 0.99581, 0.99582, 0.99583, 0.99584, 0.99585, 0.99586, 0.99587, 0.99588, 0.99589, 0.9959, 0.99591, 0.99592, 0.99593, 0.99594, 0.99595, 0.99596, 0.99597, 0.99598, 0.99599, 0.996, 0.99601, 0.99602, 0.99603, 0.99604, 0.99605, 0.99606, 0.99607, 0.99608, 0.99609, 0.9961, 0.99611, 0.99612, 0.99613, 0.99614, 0.99615, 0.99616, 0.99617, 0.99618, 0.99619, 0.9962, 0.99621, 0.99622, 0.99623, 0.99624, 0.99625, 0.99626, 0.99627, 0.99628, 0.99629, 0.9963, 0.99631, 0.99632, 0.99633, 0.99634, 0.99635, 0.99636, 0.99637, 0.99638, 0.99639, 0.9964, 0.99641, 0.99642, 0.99643, 0.99644, 0.99645, 0.99646, 0.99647, 0.99648, 0.99649, 0.99650, 0.99651, 0.99652, 0.99653, 0.99654, 0.99655, 0.99656, 0.99657, 0.99658, 0.99659, 0.9966, 0.99661, 0.99662, 0.99663, 0.99664, 0.99665, 0.99666, 0.99667, 0.99668, 0.99669, 0.9967, 0.99671, 0.99672, 0.99673, 0.99674, 0.99675, 0.99676, 0.99677, 0.99678, 0.99679, 0.9968, 0.99681, 0.99682, 0.99683, 0.99684, 0.99685, 0.99686, 0.99687, 0.99688, 0.99689, 0.9969, 0.99691, 0.99692, 0.99693, 0.99694, 0.99695, 0.99696, 0.99697, 0.99698, 0.99699, 0.997, 0.99701, 0.99702, 0.99703, 0.99704, 0.99705, 0.99706, 0.99707, 0.99708, 0.99709, 0.9971, 0.99711, 0.99712, 0.99713, 0.99714, 0.99715, 0.99716, 0.99717, 0.99718, 0.99719, 0.9972, 0.99721, 0.99722, 0.99723, 0.99724, 0.99725, 0.99726, 0.99727, 0.99728, 0.99729, 0.9973, 0.99731, 0.99732, 0.99733, 0.99734, 0.99735, 0.99736, 0.99737, 0.99738, 0.99739, 0.9974, 0.99741, 0.99742, 0.99743, 0.99744, 0.99745, 0.99746, 0.99747, 0.99748, 0.99749, 0.9975, 0.99751, 0.99752, 0.99753, 0.99754, 0.99755, 0.99756, 0.99757, 0.99758, 0.99759, 0.9976, 0.99761, 0.99762, 0.99763, 0.99764, 0.99765, 0.99766, 0.99767, 0.99768, 0.99769, 0.9977, 0.99771, 0.99772, 0.99773, 0.99774, 0.99775, 0.99776, 0.99777, 0.99778, 0.99779, 0.9978, 0.99781, 0.99782, 0.99783, 0.99784, 0.99785, 0.99786, 0.99787, 0.99788, 0.99789, 0.9979, 0.99791, 0.99792, 0.99793, 0.99794, 0.99795, 0.99796, 0.99797, 0.99798, 0.99799, 0.998, 0.99801, 0.99802, 0.99803, 0.99804, 0.99805, 0.99806, 0.99807, 0.99808, 0.99809, 0.9981, 0.99811, 0.99812, 0.99813, 0.99814, 0.99815, 0.99816, 0.99817, 0.99818, 0.99819, 0.9982, 0.99821, 0.99822, 0.99823, 0.99824, 0.99825, 0.99826, 0.99827, 0.99828, 0.99829, 0.9983, 0.99831, 0.99832, 0.99833, 0.99834, 0.99835, 0.99836, 0.99837, 0.99838, 0.99839, 0.9984, 0.99841, 0.99842, 0.99843, 0.99844, 0.99845, 0.99846, 0.99847, 0.99848, 0.99849, 0.99850, 0.99851, 0.99852, 0.99853, 0.99854, 0.99855, 0.99856, 0.99857, 0.99858, 0.99859, 0.9986, 0.99861, 0.99862, 0.99863, 0.99864, 0.99865, 0.99866, 0.99867, 0.99868, 0.99869, 0.9987, 0.99871, 0.99872, 0.99873, 0.99874, 0.99875, 0.99876, 0.99877, 0.99878, 0.99879, 0.9988, 0.99881, 0.99882, 0.99883, 0.99884, 0.99885, 0.99886, 0.99887, 0.99888, 0.99889, 0.9989, 0.99891, 0.99892, 0.99893, 0.99894, 0.99895, 0.99896, 0.99897, 0.99898, 0.99899, 0.999, 0.99901, 0.99902, 0.99903, 0.99904, 0.99905, 0.99906, 0.99907, 0.99908, 0.99909, 0.9991, 0.99911, 0.99912, 0.99913, 0.99914, 0.99915, 0.99916, 0.99917, 0.99918, 0.99919, 0.9992, 0.99921, 0.99922, 0.99923, 0.99924, 0.99925, 0.99926, 0.99927, 0.99928, 0.99929, 0.9993, 0.99931, 0.99932, 0.99933, 0.99934, 0.99935, 0.99936, 0.99937, 0.99938, 0.99939, 0.9994, 0.99941, 0.99942, 0.99943, 0.99944, 0.99945, 0.99946, 0.99947, 0.99948, 0.99949, 0.9995, 0.99951, 0.99952, 0.99953, 0.99954, 0.99955, 0.99956, 0.99957, 0.99958, 0.99959, 0.9996, 0.99961, 0.99962, 0.99963, 0.99964, 0.99965, 0.99966, 0.99967, 0.99968, 0.99969, 0.9997, 0.99971, 0.99972, 0.99973, 0.99974, 0.99975, 0.99976, 0.99977, 0.99978, 0.99979, 0.9998, 0.99981, 0.99982, 0.99983, 0.99984, 0.99985, 0.99986, 0.99987, 0.99988, 0.99989, 0.9999, 0.99991, 0.99992, 0.99993, 0.99994, 0.99995, 0.99996, 0.99997, 0.99998, 0.99999, 1.0};
    for(int k = 0; k < trials+1; k++){
      probcuts[k] = probcuts_temp[k];
    }
  }

  if(opt == 3){
    Float_t probcuts_temp[trials+1] = {0.5,  0.5001,  0.5002,  0.5003,  0.5004,  0.5005,  0.5006,  0.5007,  0.5008,  0.5009,  0.501,  0.5011,  0.5012,  0.5013,  0.5014,  0.5015,  0.5016,  0.5017,  0.5018,  0.5019,  0.502,  0.5021,  0.5022,  0.5023,  0.5024,  0.5025,  0.5026,  0.5027,  0.5028,  0.5029,  0.503,  0.5031,  0.5032,  0.5033,  0.5034,  0.5035,  0.5036,  0.5037,  0.5038,  0.5039,  0.504,  0.5041,  0.5042,  0.5043,  0.5044,  0.5045,  0.5046,  0.5047,  0.5048,  0.5049,  0.505,  0.5051,  0.5052,  0.5053,  0.5054,  0.5055,  0.5056,  0.5057,  0.5058,  0.5059,  0.506,  0.5061,  0.5062,  0.5063,  0.5064,  0.5065,  0.5066,  0.5067,  0.5068,  0.5069,  0.507,  0.5071,  0.5072,  0.5073,  0.5074,  0.5075,  0.5076,  0.5077,  0.5078,  0.5079,  0.508,  0.5081,  0.5082,  0.5083,  0.5084,  0.5085,  0.5086,  0.5087,  0.5088,  0.5089,  0.509,  0.5091,  0.5092,  0.5093,  0.5094,  0.5095,  0.5096,  0.5097,  0.5098,  0.5099,  0.51,  0.5101,  0.5102,  0.5103,  0.5104,  0.5105,  0.5106,  0.5107,  0.5108,  0.5109,  0.511,  0.5111,  0.5112,  0.5113,  0.5114,  0.5115,  0.5116,  0.5117,  0.5118,  0.5119,  0.512,  0.5121,  0.5122,  0.5123,  0.5124,  0.5125,  0.5126,  0.5127,  0.5128,  0.5129,  0.513,  0.5131,  0.5132,  0.5133,  0.5134,  0.5135,  0.5136,  0.5137,  0.5138,  0.5139,  0.514,  0.5141,  0.5142,  0.5143,  0.5144,  0.5145,  0.5146,  0.5147,  0.5148,  0.5149,  0.515,  0.5151,  0.5152,  0.5153,  0.5154,  0.5155,  0.5156,  0.5157,  0.5158,  0.5159,  0.516,  0.5161,  0.5162,  0.5163,  0.5164,  0.5165,  0.5166,  0.5167,  0.5168,  0.5169,  0.517,  0.5171,  0.5172,  0.5173,  0.5174,  0.5175,  0.5176,  0.5177,  0.5178,  0.5179,  0.518,  0.5181,  0.5182,  0.5183,  0.5184,  0.5185,  0.5186,  0.5187,  0.5188,  0.5189,  0.519,  0.5191,  0.5192,  0.5193,  0.5194,  0.5195,  0.5196,  0.5197,  0.5198,  0.5199,  0.52,  0.5201,  0.5202,  0.5203,  0.5204,  0.5205,  0.5206,  0.5207,  0.5208,  0.5209,  0.521,  0.5211,  0.5212,  0.5213,  0.5214,  0.5215,  0.5216,  0.5217,  0.5218,  0.5219,  0.522,  0.5221,  0.5222,  0.5223,  0.5224,  0.5225,  0.5226,  0.5227,  0.5228,  0.5229,  0.523,  0.5231,  0.5232,  0.5233,  0.5234,  0.5235,  0.5236,  0.5237,  0.5238,  0.5239,  0.524,  0.5241,  0.5242,  0.5243,  0.5244,  0.5245,  0.5246,  0.5247,  0.5248,  0.5249,  0.525,  0.5251,  0.5252,  0.5253,  0.5254,  0.5255,  0.5256,  0.5257,  0.5258,  0.5259,  0.526,  0.5261,  0.5262,  0.5263,  0.5264,  0.5265,  0.5266,  0.5267,  0.5268,  0.5269,  0.527,  0.5271,  0.5272,  0.5273,  0.5274,  0.5275,  0.5276,  0.5277,  0.5278,  0.5279,  0.528,  0.5281,  0.5282,  0.5283,  0.5284,  0.5285,  0.5286,  0.5287,  0.5288,  0.5289,  0.529,  0.5291,  0.5292,  0.5293,  0.5294,  0.5295,  0.5296,  0.5297,  0.5298,  0.5299,  0.53,  0.5301,  0.5302,  0.5303,  0.5304,  0.5305,  0.5306,  0.5307,  0.5308,  0.5309,  0.531,  0.5311,  0.5312,  0.5313,  0.5314,  0.5315,  0.5316,  0.5317,  0.5318,  0.5319,  0.532,  0.5321,  0.5322,  0.5323,  0.5324,  0.5325,  0.5326,  0.5327,  0.5328,  0.5329,  0.533,  0.5331,  0.5332,  0.5333,  0.5334,  0.5335,  0.5336,  0.5337,  0.5338,  0.5339,  0.534,  0.5341,  0.5342,  0.5343,  0.5344,  0.5345,  0.5346,  0.5347,  0.5348,  0.5349,  0.535,  0.5351,  0.5352,  0.5353,  0.5354,  0.5355,  0.5356,  0.5357,  0.5358,  0.5359,  0.536,  0.5361,  0.5362,  0.5363,  0.5364,  0.5365,  0.5366,  0.5367,  0.5368,  0.5369,  0.537,  0.5371,  0.5372,  0.5373,  0.5374,  0.5375,  0.5376,  0.5377,  0.5378,  0.5379,  0.538,  0.5381,  0.5382,  0.5383,  0.5384,  0.5385,  0.5386,  0.5387,  0.5388,  0.5389,  0.539,  0.5391,  0.5392,  0.5393,  0.5394,  0.5395,  0.5396,  0.5397,  0.5398,  0.5399,  0.54,  0.5401,  0.5402,  0.5403,  0.5404,  0.5405,  0.5406,  0.5407,  0.5408,  0.5409,  0.541,  0.5411,  0.5412,  0.5413,  0.5414,  0.5415,  0.5416,  0.5417,  0.5418,  0.5419,  0.542,  0.5421,  0.5422,  0.5423,  0.5424,  0.5425,  0.5426,  0.5427,  0.5428,  0.5429,  0.543,  0.5431,  0.5432,  0.5433,  0.5434,  0.5435,  0.5436,  0.5437,  0.5438,  0.5439,  0.544,  0.5441,  0.5442,  0.5443,  0.5444,  0.5445,  0.5446,  0.5447,  0.5448,  0.5449,  0.545,  0.5451,  0.5452,  0.5453,  0.5454,  0.5455,  0.5456,  0.5457,  0.5458,  0.5459,  0.546,  0.5461,  0.5462,  0.5463,  0.5464,  0.5465,  0.5466,  0.5467,  0.5468,  0.5469,  0.547,  0.5471,  0.5472,  0.5473,  0.5474,  0.5475,  0.5476,  0.5477,  0.5478,  0.5479,  0.548,  0.5481,  0.5482,  0.5483,  0.5484,  0.5485,  0.5486,  0.5487,  0.5488,  0.5489,  0.549,  0.5491,  0.5492,  0.5493,  0.5494,  0.5495,  0.5496,  0.5497,  0.5498,  0.5499,  0.55,  0.5501,  0.5502,  0.5503,  0.5504,  0.5505,  0.5506,  0.5507,  0.5508,  0.5509,  0.551,  0.5511,  0.5512,  0.5513,  0.5514,  0.5515,  0.5516,  0.5517,  0.5518,  0.5519,  0.552,  0.5521,  0.5522,  0.5523,  0.5524,  0.5525,  0.5526,  0.5527,  0.5528,  0.5529,  0.553,  0.5531,  0.5532,  0.5533,  0.5534,  0.5535,  0.5536,  0.5537,  0.5538,  0.5539,  0.554,  0.5541,  0.5542,  0.5543,  0.5544,  0.5545,  0.5546,  0.5547,  0.5548,  0.5549,  0.555,  0.5551,  0.5552,  0.5553,  0.5554,  0.5555,  0.5556,  0.5557,  0.5558,  0.5559,  0.556,  0.5561,  0.5562,  0.5563,  0.5564,  0.5565,  0.5566,  0.5567,  0.5568,  0.5569,  0.557,  0.5571,  0.5572,  0.5573,  0.5574,  0.5575,  0.5576,  0.5577,  0.5578,  0.5579,  0.558,  0.5581,  0.5582,  0.5583,  0.5584,  0.5585,  0.5586,  0.5587,  0.5588,  0.5589,  0.559,  0.5591,  0.5592,  0.5593,  0.5594,  0.5595,  0.5596,  0.5597,  0.5598,  0.5599,  0.56,  0.5601,  0.5602,  0.5603,  0.5604,  0.5605,  0.5606,  0.5607,  0.5608,  0.5609,  0.561,  0.5611,  0.5612,  0.5613,  0.5614,  0.5615,  0.5616,  0.5617,  0.5618,  0.5619,  0.562,  0.5621,  0.5622,  0.5623,  0.5624,  0.5625,  0.5626,  0.5627,  0.5628,  0.5629,  0.563,  0.5631,  0.5632,  0.5633,  0.5634,  0.5635,  0.5636,  0.5637,  0.5638,  0.5639,  0.564,  0.5641,  0.5642,  0.5643,  0.5644,  0.5645,  0.5646,  0.5647,  0.5648,  0.5649,  0.565,  0.5651,  0.5652,  0.5653,  0.5654,  0.5655,  0.5656,  0.5657,  0.5658,  0.5659,  0.566,  0.5661,  0.5662,  0.5663,  0.5664,  0.5665,  0.5666,  0.5667,  0.5668,  0.5669,  0.567,  0.5671,  0.5672,  0.5673,  0.5674,  0.5675,  0.5676,  0.5677,  0.5678,  0.5679,  0.568,  0.5681,  0.5682,  0.5683,  0.5684,  0.5685,  0.5686,  0.5687,  0.5688,  0.5689,  0.569,  0.5691,  0.5692,  0.5693,  0.5694,  0.5695,  0.5696,  0.5697,  0.5698,  0.5699,  0.57,  0.5701,  0.5702,  0.5703,  0.5704,  0.5705,  0.5706,  0.5707,  0.5708,  0.5709,  0.571,  0.5711,  0.5712,  0.5713,  0.5714,  0.5715,  0.5716,  0.5717,  0.5718,  0.5719,  0.572,  0.5721,  0.5722,  0.5723,  0.5724,  0.5725,  0.5726,  0.5727,  0.5728,  0.5729,  0.573,  0.5731,  0.5732,  0.5733,  0.5734,  0.5735,  0.5736,  0.5737,  0.5738,  0.5739,  0.574,  0.5741,  0.5742,  0.5743,  0.5744,  0.5745,  0.5746,  0.5747,  0.5748,  0.5749,  0.575,  0.5751,  0.5752,  0.5753,  0.5754,  0.5755,  0.5756,  0.5757,  0.5758,  0.5759,  0.576,  0.5761,  0.5762,  0.5763,  0.5764,  0.5765,  0.5766,  0.5767,  0.5768,  0.5769,  0.577,  0.5771,  0.5772,  0.5773,  0.5774,  0.5775,  0.5776,  0.5777,  0.5778,  0.5779,  0.578,  0.5781,  0.5782,  0.5783,  0.5784,  0.5785,  0.5786,  0.5787,  0.5788,  0.5789,  0.579,  0.5791,  0.5792,  0.5793,  0.5794,  0.5795,  0.5796,  0.5797,  0.5798,  0.5799,  0.58,  0.5801,  0.5802,  0.5803,  0.5804,  0.5805,  0.5806,  0.5807,  0.5808,  0.5809,  0.581,  0.5811,  0.5812,  0.5813,  0.5814,  0.5815,  0.5816,  0.5817,  0.5818,  0.5819,  0.582,  0.5821,  0.5822,  0.5823,  0.5824,  0.5825,  0.5826,  0.5827,  0.5828,  0.5829,  0.583,  0.5831,  0.5832,  0.5833,  0.5834,  0.5835,  0.5836,  0.5837,  0.5838,  0.5839,  0.584,  0.5841,  0.5842,  0.5843,  0.5844,  0.5845,  0.5846,  0.5847,  0.5848,  0.5849,  0.585,  0.5851,  0.5852,  0.5853,  0.5854,  0.5855,  0.5856,  0.5857,  0.5858,  0.5859,  0.586,  0.5861,  0.5862,  0.5863,  0.5864,  0.5865,  0.5866,  0.5867,  0.5868,  0.5869,  0.587,  0.5871,  0.5872,  0.5873,  0.5874,  0.5875,  0.5876,  0.5877,  0.5878,  0.5879,  0.588,  0.5881,  0.5882,  0.5883,  0.5884,  0.5885,  0.5886,  0.5887,  0.5888,  0.5889,  0.589,  0.5891,  0.5892,  0.5893,  0.5894,  0.5895,  0.5896,  0.5897,  0.5898,  0.5899,  0.59,  0.5901,  0.5902,  0.5903,  0.5904,  0.5905,  0.5906,  0.5907,  0.5908,  0.5909,  0.591,  0.5911,  0.5912,  0.5913,  0.5914,  0.5915,  0.5916,  0.5917,  0.5918,  0.5919,  0.592,  0.5921,  0.5922,  0.5923,  0.5924,  0.5925,  0.5926,  0.5927,  0.5928,  0.5929,  0.593,  0.5931,  0.5932,  0.5933,  0.5934,  0.5935,  0.5936,  0.5937,  0.5938,  0.5939,  0.594,  0.5941,  0.5942,  0.5943,  0.5944,  0.5945,  0.5946,  0.5947,  0.5948,  0.5949,  0.595,  0.5951,  0.5952,  0.5953,  0.5954,  0.5955,  0.5956,  0.5957,  0.5958,  0.5959,  0.596,  0.5961,  0.5962,  0.5963,  0.5964,  0.5965,  0.5966,  0.5967,  0.5968,  0.5969,  0.597,  0.5971,  0.5972,  0.5973,  0.5974,  0.5975,  0.5976,  0.5977,  0.5978,  0.5979,  0.598,  0.5981,  0.5982,  0.5983,  0.5984,  0.5985,  0.5986,  0.5987,  0.5988,  0.5989,  0.599,  0.5991,  0.5992,  0.5993,  0.5994,  0.5995,  0.5996,  0.5997,  0.5998,  0.5999,  0.6,  0.6001,  0.6002,  0.6003,  0.6004,  0.6005,  0.6006,  0.6007,  0.6008,  0.6009,  0.601,  0.6011,  0.6012,  0.6013,  0.6014,  0.6015,  0.6016,  0.6017,  0.6018,  0.6019,  0.602,  0.6021,  0.6022,  0.6023,  0.6024,  0.6025,  0.6026,  0.6027,  0.6028,  0.6029,  0.603,  0.6031,  0.6032,  0.6033,  0.6034,  0.6035,  0.6036,  0.6037,  0.6038,  0.6039,  0.604,  0.6041,  0.6042,  0.6043,  0.6044,  0.6045,  0.6046,  0.6047,  0.6048,  0.6049,  0.605,  0.6051,  0.6052,  0.6053,  0.6054,  0.6055,  0.6056,  0.6057,  0.6058,  0.6059,  0.606,  0.6061,  0.6062,  0.6063,  0.6064,  0.6065,  0.6066,  0.6067,  0.6068,  0.6069,  0.607,  0.6071,  0.6072,  0.6073,  0.6074,  0.6075,  0.6076,  0.6077,  0.6078,  0.6079,  0.608,  0.6081,  0.6082,  0.6083,  0.6084,  0.6085,  0.6086,  0.6087,  0.6088,  0.6089,  0.609,  0.6091,  0.6092,  0.6093,  0.6094,  0.6095,  0.6096,  0.6097,  0.6098,  0.6099,  0.61,  0.6101,  0.6102,  0.6103,  0.6104,  0.6105,  0.6106,  0.6107,  0.6108,  0.6109,  0.611,  0.6111,  0.6112,  0.6113,  0.6114,  0.6115,  0.6116,  0.6117,  0.6118,  0.6119,  0.612,  0.6121,  0.6122,  0.6123,  0.6124,  0.6125,  0.6126,  0.6127,  0.6128,  0.6129,  0.613,  0.6131,  0.6132,  0.6133,  0.6134,  0.6135,  0.6136,  0.6137,  0.6138,  0.6139,  0.614,  0.6141,  0.6142,  0.6143,  0.6144,  0.6145,  0.6146,  0.6147,  0.6148,  0.6149,  0.615,  0.6151,  0.6152,  0.6153,  0.6154,  0.6155,  0.6156,  0.6157,  0.6158,  0.6159,  0.616,  0.6161,  0.6162,  0.6163,  0.6164,  0.6165,  0.6166,  0.6167,  0.6168,  0.6169,  0.617,  0.6171,  0.6172,  0.6173,  0.6174,  0.6175,  0.6176,  0.6177,  0.6178,  0.6179,  0.618,  0.6181,  0.6182,  0.6183,  0.6184,  0.6185,  0.6186,  0.6187,  0.6188,  0.6189,  0.619,  0.6191,  0.6192,  0.6193,  0.6194,  0.6195,  0.6196,  0.6197,  0.6198,  0.6199,  0.62,  0.6201,  0.6202,  0.6203,  0.6204,  0.6205,  0.6206,  0.6207,  0.6208,  0.6209,  0.621,  0.6211,  0.6212,  0.6213,  0.6214,  0.6215,  0.6216,  0.6217,  0.6218,  0.6219,  0.622,  0.6221,  0.6222,  0.6223,  0.6224,  0.6225,  0.6226,  0.6227,  0.6228,  0.6229,  0.623,  0.6231,  0.6232,  0.6233,  0.6234,  0.6235,  0.6236,  0.6237,  0.6238,  0.6239,  0.624,  0.6241,  0.6242,  0.6243,  0.6244,  0.6245,  0.6246,  0.6247,  0.6248,  0.6249,  0.625,  0.6251,  0.6252,  0.6253,  0.6254,  0.6255,  0.6256,  0.6257,  0.6258,  0.6259,  0.626,  0.6261,  0.6262,  0.6263,  0.6264,  0.6265,  0.6266,  0.6267,  0.6268,  0.6269,  0.627,  0.6271,  0.6272,  0.6273,  0.6274,  0.6275,  0.6276,  0.6277,  0.6278,  0.6279,  0.628,  0.6281,  0.6282,  0.6283,  0.6284,  0.6285,  0.6286,  0.6287,  0.6288,  0.6289,  0.629,  0.6291,  0.6292,  0.6293,  0.6294,  0.6295,  0.6296,  0.6297,  0.6298,  0.6299,  0.63,  0.6301,  0.6302,  0.6303,  0.6304,  0.6305,  0.6306,  0.6307,  0.6308,  0.6309,  0.631,  0.6311,  0.6312,  0.6313,  0.6314,  0.6315,  0.6316,  0.6317,  0.6318,  0.6319,  0.632,  0.6321,  0.6322,  0.6323,  0.6324,  0.6325,  0.6326,  0.6327,  0.6328,  0.6329,  0.633,  0.6331,  0.6332,  0.6333,  0.6334,  0.6335,  0.6336,  0.6337,  0.6338,  0.6339,  0.634,  0.6341,  0.6342,  0.6343,  0.6344,  0.6345,  0.6346,  0.6347,  0.6348,  0.6349,  0.635,  0.6351,  0.6352,  0.6353,  0.6354,  0.6355,  0.6356,  0.6357,  0.6358,  0.6359,  0.636,  0.6361,  0.6362,  0.6363,  0.6364,  0.6365,  0.6366,  0.6367,  0.6368,  0.6369,  0.637,  0.6371,  0.6372,  0.6373,  0.6374,  0.6375,  0.6376,  0.6377,  0.6378,  0.6379,  0.638,  0.6381,  0.6382,  0.6383,  0.6384,  0.6385,  0.6386,  0.6387,  0.6388,  0.6389,  0.639,  0.6391,  0.6392,  0.6393,  0.6394,  0.6395,  0.6396,  0.6397,  0.6398,  0.6399,  0.64,  0.6401,  0.6402,  0.6403,  0.6404,  0.6405,  0.6406,  0.6407,  0.6408,  0.6409,  0.641,  0.6411,  0.6412,  0.6413,  0.6414,  0.6415,  0.6416,  0.6417,  0.6418,  0.6419,  0.642,  0.6421,  0.6422,  0.6423,  0.6424,  0.6425,  0.6426,  0.6427,  0.6428,  0.6429,  0.643,  0.6431,  0.6432,  0.6433,  0.6434,  0.6435,  0.6436,  0.6437,  0.6438,  0.6439,  0.644,  0.6441,  0.6442,  0.6443,  0.6444,  0.6445,  0.6446,  0.6447,  0.6448,  0.6449,  0.645,  0.6451,  0.6452,  0.6453,  0.6454,  0.6455,  0.6456,  0.6457,  0.6458,  0.6459,  0.646,  0.6461,  0.6462,  0.6463,  0.6464,  0.6465,  0.6466,  0.6467,  0.6468,  0.6469,  0.647,  0.6471,  0.6472,  0.6473,  0.6474,  0.6475,  0.6476,  0.6477,  0.6478,  0.6479,  0.648,  0.6481,  0.6482,  0.6483,  0.6484,  0.6485,  0.6486,  0.6487,  0.6488,  0.6489,  0.649,  0.6491,  0.6492,  0.6493,  0.6494,  0.6495,  0.6496,  0.6497,  0.6498,  0.6499,  0.65,  0.6501,  0.6502,  0.6503,  0.6504,  0.6505,  0.6506,  0.6507,  0.6508,  0.6509,  0.651,  0.6511,  0.6512,  0.6513,  0.6514,  0.6515,  0.6516,  0.6517,  0.6518,  0.6519,  0.652,  0.6521,  0.6522,  0.6523,  0.6524,  0.6525,  0.6526,  0.6527,  0.6528,  0.6529,  0.653,  0.6531,  0.6532,  0.6533,  0.6534,  0.6535,  0.6536,  0.6537,  0.6538,  0.6539,  0.654,  0.6541,  0.6542,  0.6543,  0.6544,  0.6545,  0.6546,  0.6547,  0.6548,  0.6549,  0.655,  0.6551,  0.6552,  0.6553,  0.6554,  0.6555,  0.6556,  0.6557,  0.6558,  0.6559,  0.656,  0.6561,  0.6562,  0.6563,  0.6564,  0.6565,  0.6566,  0.6567,  0.6568,  0.6569,  0.657,  0.6571,  0.6572,  0.6573,  0.6574,  0.6575,  0.6576,  0.6577,  0.6578,  0.6579,  0.658,  0.6581,  0.6582,  0.6583,  0.6584,  0.6585,  0.6586,  0.6587,  0.6588,  0.6589,  0.659,  0.6591,  0.6592,  0.6593,  0.6594,  0.6595,  0.6596,  0.6597,  0.6598,  0.6599,  0.66,  0.6601,  0.6602,  0.6603,  0.6604,  0.6605,  0.6606,  0.6607,  0.6608,  0.6609,  0.661,  0.6611,  0.6612,  0.6613,  0.6614,  0.6615,  0.6616,  0.6617,  0.6618,  0.6619,  0.662,  0.6621,  0.6622,  0.6623,  0.6624,  0.6625,  0.6626,  0.6627,  0.6628,  0.6629,  0.663,  0.6631,  0.6632,  0.6633,  0.6634,  0.6635,  0.6636,  0.6637,  0.6638,  0.6639,  0.664,  0.6641,  0.6642,  0.6643,  0.6644,  0.6645,  0.6646,  0.6647,  0.6648,  0.6649,  0.665,  0.6651,  0.6652,  0.6653,  0.6654,  0.6655,  0.6656,  0.6657,  0.6658,  0.6659,  0.666,  0.6661,  0.6662,  0.6663,  0.6664,  0.6665,  0.6666,  0.6667,  0.6668,  0.6669,  0.667,  0.6671,  0.6672,  0.6673,  0.6674,  0.6675,  0.6676,  0.6677,  0.6678,  0.6679,  0.668,  0.6681,  0.6682,  0.6683,  0.6684,  0.6685,  0.6686,  0.6687,  0.6688,  0.6689,  0.669,  0.6691,  0.6692,  0.6693,  0.6694,  0.6695,  0.6696,  0.6697,  0.6698,  0.6699,  0.67,  0.6701,  0.6702,  0.6703,  0.6704,  0.6705,  0.6706,  0.6707,  0.6708,  0.6709,  0.671,  0.6711,  0.6712,  0.6713,  0.6714,  0.6715,  0.6716,  0.6717,  0.6718,  0.6719,  0.672,  0.6721,  0.6722,  0.6723,  0.6724,  0.6725,  0.6726,  0.6727,  0.6728,  0.6729,  0.673,  0.6731,  0.6732,  0.6733,  0.6734,  0.6735,  0.6736,  0.6737,  0.6738,  0.6739,  0.674,  0.6741,  0.6742,  0.6743,  0.6744,  0.6745,  0.6746,  0.6747,  0.6748,  0.6749,  0.675,  0.6751,  0.6752,  0.6753,  0.6754,  0.6755,  0.6756,  0.6757,  0.6758,  0.6759,  0.676,  0.6761,  0.6762,  0.6763,  0.6764,  0.6765,  0.6766,  0.6767,  0.6768,  0.6769,  0.677,  0.6771,  0.6772,  0.6773,  0.6774,  0.6775,  0.6776,  0.6777,  0.6778,  0.6779,  0.678,  0.6781,  0.6782,  0.6783,  0.6784,  0.6785,  0.6786,  0.6787,  0.6788,  0.6789,  0.679,  0.6791,  0.6792,  0.6793,  0.6794,  0.6795,  0.6796,  0.6797,  0.6798,  0.6799,  0.68,  0.6801,  0.6802,  0.6803,  0.6804,  0.6805,  0.6806,  0.6807,  0.6808,  0.6809,  0.681,  0.6811,  0.6812,  0.6813,  0.6814,  0.6815,  0.6816,  0.6817,  0.6818,  0.6819,  0.682,  0.6821,  0.6822,  0.6823,  0.6824,  0.6825,  0.6826,  0.6827,  0.6828,  0.6829,  0.683,  0.6831,  0.6832,  0.6833,  0.6834,  0.6835,  0.6836,  0.6837,  0.6838,  0.6839,  0.684,  0.6841,  0.6842,  0.6843,  0.6844,  0.6845,  0.6846,  0.6847,  0.6848,  0.6849,  0.685,  0.6851,  0.6852,  0.6853,  0.6854,  0.6855,  0.6856,  0.6857,  0.6858,  0.6859,  0.686,  0.6861,  0.6862,  0.6863,  0.6864,  0.6865,  0.6866,  0.6867,  0.6868,  0.6869,  0.687,  0.6871,  0.6872,  0.6873,  0.6874,  0.6875,  0.6876,  0.6877,  0.6878,  0.6879,  0.688,  0.6881,  0.6882,  0.6883,  0.6884,  0.6885,  0.6886,  0.6887,  0.6888,  0.6889,  0.689,  0.6891,  0.6892,  0.6893,  0.6894,  0.6895,  0.6896,  0.6897,  0.6898,  0.6899,  0.69,  0.6901,  0.6902,  0.6903,  0.6904,  0.6905,  0.6906,  0.6907,  0.6908,  0.6909,  0.691,  0.6911,  0.6912,  0.6913,  0.6914,  0.6915,  0.6916,  0.6917,  0.6918,  0.6919,  0.692,  0.6921,  0.6922,  0.6923,  0.6924,  0.6925,  0.6926,  0.6927,  0.6928,  0.6929,  0.693,  0.6931,  0.6932,  0.6933,  0.6934,  0.6935,  0.6936,  0.6937,  0.6938,  0.6939,  0.694,  0.6941,  0.6942,  0.6943,  0.6944,  0.6945,  0.6946,  0.6947,  0.6948,  0.6949,  0.695,  0.6951,  0.6952,  0.6953,  0.6954,  0.6955,  0.6956,  0.6957,  0.6958,  0.6959,  0.696,  0.6961,  0.6962,  0.6963,  0.6964,  0.6965,  0.6966,  0.6967,  0.6968,  0.6969,  0.697,  0.6971,  0.6972,  0.6973,  0.6974,  0.6975,  0.6976,  0.6977,  0.6978,  0.6979,  0.698,  0.6981,  0.6982,  0.6983,  0.6984,  0.6985,  0.6986,  0.6987,  0.6988,  0.6989,  0.699,  0.6991,  0.6992,  0.6993,  0.6994,  0.6995,  0.6996,  0.6997,  0.6998,  0.6999,  0.7,  0.7001,  0.7002,  0.7003,  0.7004,  0.7005,  0.7006,  0.7007,  0.7008,  0.7009,  0.701,  0.7011,  0.7012,  0.7013,  0.7014,  0.7015,  0.7016,  0.7017,  0.7018,  0.7019,  0.702,  0.7021,  0.7022,  0.7023,  0.7024,  0.7025,  0.7026,  0.7027,  0.7028,  0.7029,  0.703,  0.7031,  0.7032,  0.7033,  0.7034,  0.7035,  0.7036,  0.7037,  0.7038,  0.7039,  0.704,  0.7041,  0.7042,  0.7043,  0.7044,  0.7045,  0.7046,  0.7047,  0.7048,  0.7049,  0.705,  0.7051,  0.7052,  0.7053,  0.7054,  0.7055,  0.7056,  0.7057,  0.7058,  0.7059,  0.706,  0.7061,  0.7062,  0.7063,  0.7064,  0.7065,  0.7066,  0.7067,  0.7068,  0.7069,  0.707,  0.7071,  0.7072,  0.7073,  0.7074,  0.7075,  0.7076,  0.7077,  0.7078,  0.7079,  0.708,  0.7081,  0.7082,  0.7083,  0.7084,  0.7085,  0.7086,  0.7087,  0.7088,  0.7089,  0.709,  0.7091,  0.7092,  0.7093,  0.7094,  0.7095,  0.7096,  0.7097,  0.7098,  0.7099,  0.71,  0.7101,  0.7102,  0.7103,  0.7104,  0.7105,  0.7106,  0.7107,  0.7108,  0.7109,  0.711,  0.7111,  0.7112,  0.7113,  0.7114,  0.7115,  0.7116,  0.7117,  0.7118,  0.7119,  0.712,  0.7121,  0.7122,  0.7123,  0.7124,  0.7125,  0.7126,  0.7127,  0.7128,  0.7129,  0.713,  0.7131,  0.7132,  0.7133,  0.7134,  0.7135,  0.7136,  0.7137,  0.7138,  0.7139,  0.714,  0.7141,  0.7142,  0.7143,  0.7144,  0.7145,  0.7146,  0.7147,  0.7148,  0.7149,  0.715,  0.7151,  0.7152,  0.7153,  0.7154,  0.7155,  0.7156,  0.7157,  0.7158,  0.7159,  0.716,  0.7161,  0.7162,  0.7163,  0.7164,  0.7165,  0.7166,  0.7167,  0.7168,  0.7169,  0.717,  0.7171,  0.7172,  0.7173,  0.7174,  0.7175,  0.7176,  0.7177,  0.7178,  0.7179,  0.718,  0.7181,  0.7182,  0.7183,  0.7184,  0.7185,  0.7186,  0.7187,  0.7188,  0.7189,  0.719,  0.7191,  0.7192,  0.7193,  0.7194,  0.7195,  0.7196,  0.7197,  0.7198,  0.7199,  0.72,  0.7201,  0.7202,  0.7203,  0.7204,  0.7205,  0.7206,  0.7207,  0.7208,  0.7209,  0.721,  0.7211,  0.7212,  0.7213,  0.7214,  0.7215,  0.7216,  0.7217,  0.7218,  0.7219,  0.722,  0.7221,  0.7222,  0.7223,  0.7224,  0.7225,  0.7226,  0.7227,  0.7228,  0.7229,  0.723,  0.7231,  0.7232,  0.7233,  0.7234,  0.7235,  0.7236,  0.7237,  0.7238,  0.7239,  0.724,  0.7241,  0.7242,  0.7243,  0.7244,  0.7245,  0.7246,  0.7247,  0.7248,  0.7249,  0.725,  0.7251,  0.7252,  0.7253,  0.7254,  0.7255,  0.7256,  0.7257,  0.7258,  0.7259,  0.726,  0.7261,  0.7262,  0.7263,  0.7264,  0.7265,  0.7266,  0.7267,  0.7268,  0.7269,  0.727,  0.7271,  0.7272,  0.7273,  0.7274,  0.7275,  0.7276,  0.7277,  0.7278,  0.7279,  0.728,  0.7281,  0.7282,  0.7283,  0.7284,  0.7285,  0.7286,  0.7287,  0.7288,  0.7289,  0.729,  0.7291,  0.7292,  0.7293,  0.7294,  0.7295,  0.7296,  0.7297,  0.7298,  0.7299,  0.73,  0.7301,  0.7302,  0.7303,  0.7304,  0.7305,  0.7306,  0.7307,  0.7308,  0.7309,  0.731,  0.7311,  0.7312,  0.7313,  0.7314,  0.7315,  0.7316,  0.7317,  0.7318,  0.7319,  0.732,  0.7321,  0.7322,  0.7323,  0.7324,  0.7325,  0.7326,  0.7327,  0.7328,  0.7329,  0.733,  0.7331,  0.7332,  0.7333,  0.7334,  0.7335,  0.7336,  0.7337,  0.7338,  0.7339,  0.734,  0.7341,  0.7342,  0.7343,  0.7344,  0.7345,  0.7346,  0.7347,  0.7348,  0.7349,  0.735,  0.7351,  0.7352,  0.7353,  0.7354,  0.7355,  0.7356,  0.7357,  0.7358,  0.7359,  0.736,  0.7361,  0.7362,  0.7363,  0.7364,  0.7365,  0.7366,  0.7367,  0.7368,  0.7369,  0.737,  0.7371,  0.7372,  0.7373,  0.7374,  0.7375,  0.7376,  0.7377,  0.7378,  0.7379,  0.738,  0.7381,  0.7382,  0.7383,  0.7384,  0.7385,  0.7386,  0.7387,  0.7388,  0.7389,  0.739,  0.7391,  0.7392,  0.7393,  0.7394,  0.7395,  0.7396,  0.7397,  0.7398,  0.7399,  0.74,  0.7401,  0.7402,  0.7403,  0.7404,  0.7405,  0.7406,  0.7407,  0.7408,  0.7409,  0.741,  0.7411,  0.7412,  0.7413,  0.7414,  0.7415,  0.7416,  0.7417,  0.7418,  0.7419,  0.742,  0.7421,  0.7422,  0.7423,  0.7424,  0.7425,  0.7426,  0.7427,  0.7428,  0.7429,  0.743,  0.7431,  0.7432,  0.7433,  0.7434,  0.7435,  0.7436,  0.7437,  0.7438,  0.7439,  0.744,  0.7441,  0.7442,  0.7443,  0.7444,  0.7445,  0.7446,  0.7447,  0.7448,  0.7449,  0.745,  0.7451,  0.7452,  0.7453,  0.7454,  0.7455,  0.7456,  0.7457,  0.7458,  0.7459,  0.746,  0.7461,  0.7462,  0.7463,  0.7464,  0.7465,  0.7466,  0.7467,  0.7468,  0.7469,  0.747,  0.7471,  0.7472,  0.7473,  0.7474,  0.7475,  0.7476,  0.7477,  0.7478,  0.7479,  0.748,  0.7481,  0.7482,  0.7483,  0.7484,  0.7485,  0.7486,  0.7487,  0.7488,  0.7489,  0.749,  0.7491,  0.7492,  0.7493,  0.7494,  0.7495,  0.7496,  0.7497,  0.7498,  0.7499,  0.75,  0.7501,  0.7502,  0.7503,  0.7504,  0.7505,  0.7506,  0.7507,  0.7508,  0.7509,  0.751,  0.7511,  0.7512,  0.7513,  0.7514,  0.7515,  0.7516,  0.7517,  0.7518,  0.7519,  0.752,  0.7521,  0.7522,  0.7523,  0.7524,  0.7525,  0.7526,  0.7527,  0.7528,  0.7529,  0.753,  0.7531,  0.7532,  0.7533,  0.7534,  0.7535,  0.7536,  0.7537,  0.7538,  0.7539,  0.754,  0.7541,  0.7542,  0.7543,  0.7544,  0.7545,  0.7546,  0.7547,  0.7548,  0.7549,  0.755,  0.7551,  0.7552,  0.7553,  0.7554,  0.7555,  0.7556,  0.7557,  0.7558,  0.7559,  0.756,  0.7561,  0.7562,  0.7563,  0.7564,  0.7565,  0.7566,  0.7567,  0.7568,  0.7569,  0.757,  0.7571,  0.7572,  0.7573,  0.7574,  0.7575,  0.7576,  0.7577,  0.7578,  0.7579,  0.758,  0.7581,  0.7582,  0.7583,  0.7584,  0.7585,  0.7586,  0.7587,  0.7588,  0.7589,  0.759,  0.7591,  0.7592,  0.7593,  0.7594,  0.7595,  0.7596,  0.7597,  0.7598,  0.7599,  0.76,  0.7601,  0.7602,  0.7603,  0.7604,  0.7605,  0.7606,  0.7607,  0.7608,  0.7609,  0.761,  0.7611,  0.7612,  0.7613,  0.7614,  0.7615,  0.7616,  0.7617,  0.7618,  0.7619,  0.762,  0.7621,  0.7622,  0.7623,  0.7624,  0.7625,  0.7626,  0.7627,  0.7628,  0.7629,  0.763,  0.7631,  0.7632,  0.7633,  0.7634,  0.7635,  0.7636,  0.7637,  0.7638,  0.7639,  0.764,  0.7641,  0.7642,  0.7643,  0.7644,  0.7645,  0.7646,  0.7647,  0.7648,  0.7649,  0.765,  0.7651,  0.7652,  0.7653,  0.7654,  0.7655,  0.7656,  0.7657,  0.7658,  0.7659,  0.766,  0.7661,  0.7662,  0.7663,  0.7664,  0.7665,  0.7666,  0.7667,  0.7668,  0.7669,  0.767,  0.7671,  0.7672,  0.7673,  0.7674,  0.7675,  0.7676,  0.7677,  0.7678,  0.7679,  0.768,  0.7681,  0.7682,  0.7683,  0.7684,  0.7685,  0.7686,  0.7687,  0.7688,  0.7689,  0.769,  0.7691,  0.7692,  0.7693,  0.7694,  0.7695,  0.7696,  0.7697,  0.7698,  0.7699,  0.77,  0.7701,  0.7702,  0.7703,  0.7704,  0.7705,  0.7706,  0.7707,  0.7708,  0.7709,  0.771,  0.7711,  0.7712,  0.7713,  0.7714,  0.7715,  0.7716,  0.7717,  0.7718,  0.7719,  0.772,  0.7721,  0.7722,  0.7723,  0.7724,  0.7725,  0.7726,  0.7727,  0.7728,  0.7729,  0.773,  0.7731,  0.7732,  0.7733,  0.7734,  0.7735,  0.7736,  0.7737,  0.7738,  0.7739,  0.774,  0.7741,  0.7742,  0.7743,  0.7744,  0.7745,  0.7746,  0.7747,  0.7748,  0.7749,  0.775,  0.7751,  0.7752,  0.7753,  0.7754,  0.7755,  0.7756,  0.7757,  0.7758,  0.7759,  0.776,  0.7761,  0.7762,  0.7763,  0.7764,  0.7765,  0.7766,  0.7767,  0.7768,  0.7769,  0.777,  0.7771,  0.7772,  0.7773,  0.7774,  0.7775,  0.7776,  0.7777,  0.7778,  0.7779,  0.778,  0.7781,  0.7782,  0.7783,  0.7784,  0.7785,  0.7786,  0.7787,  0.7788,  0.7789,  0.779,  0.7791,  0.7792,  0.7793,  0.7794,  0.7795,  0.7796,  0.7797,  0.7798,  0.7799,  0.78,  0.7801,  0.7802,  0.7803,  0.7804,  0.7805,  0.7806,  0.7807,  0.7808,  0.7809,  0.781,  0.7811,  0.7812,  0.7813,  0.7814,  0.7815,  0.7816,  0.7817,  0.7818,  0.7819,  0.782,  0.7821,  0.7822,  0.7823,  0.7824,  0.7825,  0.7826,  0.7827,  0.7828,  0.7829,  0.783,  0.7831,  0.7832,  0.7833,  0.7834,  0.7835,  0.7836,  0.7837,  0.7838,  0.7839,  0.784,  0.7841,  0.7842,  0.7843,  0.7844,  0.7845,  0.7846,  0.7847,  0.7848,  0.7849,  0.785,  0.7851,  0.7852,  0.7853,  0.7854,  0.7855,  0.7856,  0.7857,  0.7858,  0.7859,  0.786,  0.7861,  0.7862,  0.7863,  0.7864,  0.7865,  0.7866,  0.7867,  0.7868,  0.7869,  0.787,  0.7871,  0.7872,  0.7873,  0.7874,  0.7875,  0.7876,  0.7877,  0.7878,  0.7879,  0.788,  0.7881,  0.7882,  0.7883,  0.7884,  0.7885,  0.7886,  0.7887,  0.7888,  0.7889,  0.789,  0.7891,  0.7892,  0.7893,  0.7894,  0.7895,  0.7896,  0.7897,  0.7898,  0.7899,  0.79,  0.7901,  0.7902,  0.7903,  0.7904,  0.7905,  0.7906,  0.7907,  0.7908,  0.7909,  0.791,  0.7911,  0.7912,  0.7913,  0.7914,  0.7915,  0.7916,  0.7917,  0.7918,  0.7919,  0.792,  0.7921,  0.7922,  0.7923,  0.7924,  0.7925,  0.7926,  0.7927,  0.7928,  0.7929,  0.793,  0.7931,  0.7932,  0.7933,  0.7934,  0.7935,  0.7936,  0.7937,  0.7938,  0.7939,  0.794,  0.7941,  0.7942,  0.7943,  0.7944,  0.7945,  0.7946,  0.7947,  0.7948,  0.7949,  0.795,  0.7951,  0.7952,  0.7953,  0.7954,  0.7955,  0.7956,  0.7957,  0.7958,  0.7959,  0.796,  0.7961,  0.7962,  0.7963,  0.7964,  0.7965,  0.7966,  0.7967,  0.7968,  0.7969,  0.797,  0.7971,  0.7972,  0.7973,  0.7974,  0.7975,  0.7976,  0.7977,  0.7978,  0.7979,  0.798,  0.7981,  0.7982,  0.7983,  0.7984,  0.7985,  0.7986,  0.7987,  0.7988,  0.7989,  0.799,  0.7991,  0.7992,  0.7993,  0.7994,  0.7995,  0.7996,  0.7997,  0.7998,  0.7999,  0.8,  0.8001,  0.8002,  0.8003,  0.8004,  0.8005,  0.8006,  0.8007,  0.8008,  0.8009,  0.801,  0.8011,  0.8012,  0.8013,  0.8014,  0.8015,  0.8016,  0.8017,  0.8018,  0.8019,  0.802,  0.8021,  0.8022,  0.8023,  0.8024,  0.8025,  0.8026,  0.8027,  0.8028,  0.8029,  0.803,  0.8031,  0.8032,  0.8033,  0.8034,  0.8035,  0.8036,  0.8037,  0.8038,  0.8039,  0.804,  0.8041,  0.8042,  0.8043,  0.8044,  0.8045,  0.8046,  0.8047,  0.8048,  0.8049,  0.805,  0.8051,  0.8052,  0.8053,  0.8054,  0.8055,  0.8056,  0.8057,  0.8058,  0.8059,  0.806,  0.8061,  0.8062,  0.8063,  0.8064,  0.8065,  0.8066,  0.8067,  0.8068,  0.8069,  0.807,  0.8071,  0.8072,  0.8073,  0.8074,  0.8075,  0.8076,  0.8077,  0.8078,  0.8079,  0.808,  0.8081,  0.8082,  0.8083,  0.8084,  0.8085,  0.8086,  0.8087,  0.8088,  0.8089,  0.809,  0.8091,  0.8092,  0.8093,  0.8094,  0.8095,  0.8096,  0.8097,  0.8098,  0.8099,  0.81,  0.8101,  0.8102,  0.8103,  0.8104,  0.8105,  0.8106,  0.8107,  0.8108,  0.8109,  0.811,  0.8111,  0.8112,  0.8113,  0.8114,  0.8115,  0.8116,  0.8117,  0.8118,  0.8119,  0.812,  0.8121,  0.8122,  0.8123,  0.8124,  0.8125,  0.8126,  0.8127,  0.8128,  0.8129,  0.813,  0.8131,  0.8132,  0.8133,  0.8134,  0.8135,  0.8136,  0.8137,  0.8138,  0.8139,  0.814,  0.8141,  0.8142,  0.8143,  0.8144,  0.8145,  0.8146,  0.8147,  0.8148,  0.8149,  0.815,  0.8151,  0.8152,  0.8153,  0.8154,  0.8155,  0.8156,  0.8157,  0.8158,  0.8159,  0.816,  0.8161,  0.8162,  0.8163,  0.8164,  0.8165,  0.8166,  0.8167,  0.8168,  0.8169,  0.817,  0.8171,  0.8172,  0.8173,  0.8174,  0.8175,  0.8176,  0.8177,  0.8178,  0.8179,  0.818,  0.8181,  0.8182,  0.8183,  0.8184,  0.8185,  0.8186,  0.8187,  0.8188,  0.8189,  0.819,  0.8191,  0.8192,  0.8193,  0.8194,  0.8195,  0.8196,  0.8197,  0.8198,  0.8199,  0.82,  0.8201,  0.8202,  0.8203,  0.8204,  0.8205,  0.8206,  0.8207,  0.8208,  0.8209,  0.821,  0.8211,  0.8212,  0.8213,  0.8214,  0.8215,  0.8216,  0.8217,  0.8218,  0.8219,  0.822,  0.8221,  0.8222,  0.8223,  0.8224,  0.8225,  0.8226,  0.8227,  0.8228,  0.8229,  0.823,  0.8231,  0.8232,  0.8233,  0.8234,  0.8235,  0.8236,  0.8237,  0.8238,  0.8239,  0.824,  0.8241,  0.8242,  0.8243,  0.8244,  0.8245,  0.8246,  0.8247,  0.8248,  0.8249,  0.825,  0.8251,  0.8252,  0.8253,  0.8254,  0.8255,  0.8256,  0.8257,  0.8258,  0.8259,  0.826,  0.8261,  0.8262,  0.8263,  0.8264,  0.8265,  0.8266,  0.8267,  0.8268,  0.8269,  0.827,  0.8271,  0.8272,  0.8273,  0.8274,  0.8275,  0.8276,  0.8277,  0.8278,  0.8279,  0.828,  0.8281,  0.8282,  0.8283,  0.8284,  0.8285,  0.8286,  0.8287,  0.8288,  0.8289,  0.829,  0.8291,  0.8292,  0.8293,  0.8294,  0.8295,  0.8296,  0.8297,  0.8298,  0.8299,  0.83,  0.8301,  0.8302,  0.8303,  0.8304,  0.8305,  0.8306,  0.8307,  0.8308,  0.8309,  0.831,  0.8311,  0.8312,  0.8313,  0.8314,  0.8315,  0.8316,  0.8317,  0.8318,  0.8319,  0.832,  0.8321,  0.8322,  0.8323,  0.8324,  0.8325,  0.8326,  0.8327,  0.8328,  0.8329,  0.833,  0.8331,  0.8332,  0.8333,  0.8334,  0.8335,  0.8336,  0.8337,  0.8338,  0.8339,  0.834,  0.8341,  0.8342,  0.8343,  0.8344,  0.8345,  0.8346,  0.8347,  0.8348,  0.8349,  0.835,  0.8351,  0.8352,  0.8353,  0.8354,  0.8355,  0.8356,  0.8357,  0.8358,  0.8359,  0.836,  0.8361,  0.8362,  0.8363,  0.8364,  0.8365,  0.8366,  0.8367,  0.8368,  0.8369,  0.837,  0.8371,  0.8372,  0.8373,  0.8374,  0.8375,  0.8376,  0.8377,  0.8378,  0.8379,  0.838,  0.8381,  0.8382,  0.8383,  0.8384,  0.8385,  0.8386,  0.8387,  0.8388,  0.8389,  0.839,  0.8391,  0.8392,  0.8393,  0.8394,  0.8395,  0.8396,  0.8397,  0.8398,  0.8399,  0.84,  0.8401,  0.8402,  0.8403,  0.8404,  0.8405,  0.8406,  0.8407,  0.8408,  0.8409,  0.841,  0.8411,  0.8412,  0.8413,  0.8414,  0.8415,  0.8416,  0.8417,  0.8418,  0.8419,  0.842,  0.8421,  0.8422,  0.8423,  0.8424,  0.8425,  0.8426,  0.8427,  0.8428,  0.8429,  0.843,  0.8431,  0.8432,  0.8433,  0.8434,  0.8435,  0.8436,  0.8437,  0.8438,  0.8439,  0.844,  0.8441,  0.8442,  0.8443,  0.8444,  0.8445,  0.8446,  0.8447,  0.8448,  0.8449,  0.845,  0.8451,  0.8452,  0.8453,  0.8454,  0.8455,  0.8456,  0.8457,  0.8458,  0.8459,  0.846,  0.8461,  0.8462,  0.8463,  0.8464,  0.8465,  0.8466,  0.8467,  0.8468,  0.8469,  0.847,  0.8471,  0.8472,  0.8473,  0.8474,  0.8475,  0.8476,  0.8477,  0.8478,  0.8479,  0.848,  0.8481,  0.8482,  0.8483,  0.8484,  0.8485,  0.8486,  0.8487,  0.8488,  0.8489,  0.849,  0.8491,  0.8492,  0.8493,  0.8494,  0.8495,  0.8496,  0.8497,  0.8498,  0.8499,  0.85,  0.8501,  0.8502,  0.8503,  0.8504,  0.8505,  0.8506,  0.8507,  0.8508,  0.8509,  0.851,  0.8511,  0.8512,  0.8513,  0.8514,  0.8515,  0.8516,  0.8517,  0.8518,  0.8519,  0.852,  0.8521,  0.8522,  0.8523,  0.8524,  0.8525,  0.8526,  0.8527,  0.8528,  0.8529,  0.853,  0.8531,  0.8532,  0.8533,  0.8534,  0.8535,  0.8536,  0.8537,  0.8538,  0.8539,  0.854,  0.8541,  0.8542,  0.8543,  0.8544,  0.8545,  0.8546,  0.8547,  0.8548,  0.8549,  0.855,  0.8551,  0.8552,  0.8553,  0.8554,  0.8555,  0.8556,  0.8557,  0.8558,  0.8559,  0.856,  0.8561,  0.8562,  0.8563,  0.8564,  0.8565,  0.8566,  0.8567,  0.8568,  0.8569,  0.857,  0.8571,  0.8572,  0.8573,  0.8574,  0.8575,  0.8576,  0.8577,  0.8578,  0.8579,  0.858,  0.8581,  0.8582,  0.8583,  0.8584,  0.8585,  0.8586,  0.8587,  0.8588,  0.8589,  0.859,  0.8591,  0.8592,  0.8593,  0.8594,  0.8595,  0.8596,  0.8597,  0.8598,  0.8599,  0.86,  0.8601,  0.8602,  0.8603,  0.8604,  0.8605,  0.8606,  0.8607,  0.8608,  0.8609,  0.861,  0.8611,  0.8612,  0.8613,  0.8614,  0.8615,  0.8616,  0.8617,  0.8618,  0.8619,  0.862,  0.8621,  0.8622,  0.8623,  0.8624,  0.8625,  0.8626,  0.8627,  0.8628,  0.8629,  0.863,  0.8631,  0.8632,  0.8633,  0.8634,  0.8635,  0.8636,  0.8637,  0.8638,  0.8639,  0.864,  0.8641,  0.8642,  0.8643,  0.8644,  0.8645,  0.8646,  0.8647,  0.8648,  0.8649,  0.865,  0.8651,  0.8652,  0.8653,  0.8654,  0.8655,  0.8656,  0.8657,  0.8658,  0.8659,  0.866,  0.8661,  0.8662,  0.8663,  0.8664,  0.8665,  0.8666,  0.8667,  0.8668,  0.8669,  0.867,  0.8671,  0.8672,  0.8673,  0.8674,  0.8675,  0.8676,  0.8677,  0.8678,  0.8679,  0.868,  0.8681,  0.8682,  0.8683,  0.8684,  0.8685,  0.8686,  0.8687,  0.8688,  0.8689,  0.869,  0.8691,  0.8692,  0.8693,  0.8694,  0.8695,  0.8696,  0.8697,  0.8698,  0.8699,  0.87,  0.8701,  0.8702,  0.8703,  0.8704,  0.8705,  0.8706,  0.8707,  0.8708,  0.8709,  0.871,  0.8711,  0.8712,  0.8713,  0.8714,  0.8715,  0.8716,  0.8717,  0.8718,  0.8719,  0.872,  0.8721,  0.8722,  0.8723,  0.8724,  0.8725,  0.8726,  0.8727,  0.8728,  0.8729,  0.873,  0.8731,  0.8732,  0.8733,  0.8734,  0.8735,  0.8736,  0.8737,  0.8738,  0.8739,  0.874,  0.8741,  0.8742,  0.8743,  0.8744,  0.8745,  0.8746,  0.8747,  0.8748,  0.8749,  0.875,  0.8751,  0.8752,  0.8753,  0.8754,  0.8755,  0.8756,  0.8757,  0.8758,  0.8759,  0.876,  0.8761,  0.8762,  0.8763,  0.8764,  0.8765,  0.8766,  0.8767,  0.8768,  0.8769,  0.877,  0.8771,  0.8772,  0.8773,  0.8774,  0.8775,  0.8776,  0.8777,  0.8778,  0.8779,  0.878,  0.8781,  0.8782,  0.8783,  0.8784,  0.8785,  0.8786,  0.8787,  0.8788,  0.8789,  0.879,  0.8791,  0.8792,  0.8793,  0.8794,  0.8795,  0.8796,  0.8797,  0.8798,  0.8799,  0.88,  0.8801,  0.8802,  0.8803,  0.8804,  0.8805,  0.8806,  0.8807,  0.8808,  0.8809,  0.881,  0.8811,  0.8812,  0.8813,  0.8814,  0.8815,  0.8816,  0.8817,  0.8818,  0.8819,  0.882,  0.8821,  0.8822,  0.8823,  0.8824,  0.8825,  0.8826,  0.8827,  0.8828,  0.8829,  0.883,  0.8831,  0.8832,  0.8833,  0.8834,  0.8835,  0.8836,  0.8837,  0.8838,  0.8839,  0.884,  0.8841,  0.8842,  0.8843,  0.8844,  0.8845,  0.8846,  0.8847,  0.8848,  0.8849,  0.885,  0.8851,  0.8852,  0.8853,  0.8854,  0.8855,  0.8856,  0.8857,  0.8858,  0.8859,  0.886,  0.8861,  0.8862,  0.8863,  0.8864,  0.8865,  0.8866,  0.8867,  0.8868,  0.8869,  0.887,  0.8871,  0.8872,  0.8873,  0.8874,  0.8875,  0.8876,  0.8877,  0.8878,  0.8879,  0.888,  0.8881,  0.8882,  0.8883,  0.8884,  0.8885,  0.8886,  0.8887,  0.8888,  0.8889,  0.889,  0.8891,  0.8892,  0.8893,  0.8894,  0.8895,  0.8896,  0.8897,  0.8898,  0.8899,  0.89,  0.8901,  0.8902,  0.8903,  0.8904,  0.8905,  0.8906,  0.8907,  0.8908,  0.8909,  0.891,  0.8911,  0.8912,  0.8913,  0.8914,  0.8915,  0.8916,  0.8917,  0.8918,  0.8919,  0.892,  0.8921,  0.8922,  0.8923,  0.8924,  0.8925,  0.8926,  0.8927,  0.8928,  0.8929,  0.893,  0.8931,  0.8932,  0.8933,  0.8934,  0.8935,  0.8936,  0.8937,  0.8938,  0.8939,  0.894,  0.8941,  0.8942,  0.8943,  0.8944,  0.8945,  0.8946,  0.8947,  0.8948,  0.8949,  0.895,  0.8951,  0.8952,  0.8953,  0.8954,  0.8955,  0.8956,  0.8957,  0.8958,  0.8959,  0.896,  0.8961,  0.8962,  0.8963,  0.8964,  0.8965,  0.8966,  0.8967,  0.8968,  0.8969,  0.897,  0.8971,  0.8972,  0.8973,  0.8974,  0.8975,  0.8976,  0.8977,  0.8978,  0.8979,  0.898,  0.8981,  0.8982,  0.8983,  0.8984,  0.8985,  0.8986,  0.8987,  0.8988,  0.8989,  0.899,  0.8991,  0.8992,  0.8993,  0.8994,  0.8995,  0.8996,  0.8997,  0.8998,  0.8999,  0.9,  0.9001,  0.9002,  0.9003,  0.9004,  0.9005,  0.9006,  0.9007,  0.9008,  0.9009,  0.901,  0.9011,  0.9012,  0.9013,  0.9014,  0.9015,  0.9016,  0.9017,  0.9018,  0.9019,  0.902,  0.9021,  0.9022,  0.9023,  0.9024,  0.9025,  0.9026,  0.9027,  0.9028,  0.9029,  0.903,  0.9031,  0.9032,  0.9033,  0.9034,  0.9035,  0.9036,  0.9037,  0.9038,  0.9039,  0.904,  0.9041,  0.9042,  0.9043,  0.9044,  0.9045,  0.9046,  0.9047,  0.9048,  0.9049,  0.905,  0.9051,  0.9052,  0.9053,  0.9054,  0.9055,  0.9056,  0.9057,  0.9058,  0.9059,  0.906,  0.9061,  0.9062,  0.9063,  0.9064,  0.9065,  0.9066,  0.9067,  0.9068,  0.9069,  0.907,  0.9071,  0.9072,  0.9073,  0.9074,  0.9075,  0.9076,  0.9077,  0.9078,  0.9079,  0.908,  0.9081,  0.9082,  0.9083,  0.9084,  0.9085,  0.9086,  0.9087,  0.9088,  0.9089,  0.909,  0.9091,  0.9092,  0.9093,  0.9094,  0.9095,  0.9096,  0.9097,  0.9098,  0.9099,  0.91,  0.9101,  0.9102,  0.9103,  0.9104,  0.9105,  0.9106,  0.9107,  0.9108,  0.9109,  0.911,  0.9111,  0.9112,  0.9113,  0.9114,  0.9115,  0.9116,  0.9117,  0.9118,  0.9119,  0.912,  0.9121,  0.9122,  0.9123,  0.9124,  0.9125,  0.9126,  0.9127,  0.9128,  0.9129,  0.913,  0.9131,  0.9132,  0.9133,  0.9134,  0.9135,  0.9136,  0.9137,  0.9138,  0.9139,  0.914,  0.9141,  0.9142,  0.9143,  0.9144,  0.9145,  0.9146,  0.9147,  0.9148,  0.9149,  0.915,  0.9151,  0.9152,  0.9153,  0.9154,  0.9155,  0.9156,  0.9157,  0.9158,  0.9159,  0.916,  0.9161,  0.9162,  0.9163,  0.9164,  0.9165,  0.9166,  0.9167,  0.9168,  0.9169,  0.917,  0.9171,  0.9172,  0.9173,  0.9174,  0.9175,  0.9176,  0.9177,  0.9178,  0.9179,  0.918,  0.9181,  0.9182,  0.9183,  0.9184,  0.9185,  0.9186,  0.9187,  0.9188,  0.9189,  0.919,  0.9191,  0.9192,  0.9193,  0.9194,  0.9195,  0.9196,  0.9197,  0.9198,  0.9199,  0.92,  0.9201,  0.9202,  0.9203,  0.9204,  0.9205,  0.9206,  0.9207,  0.9208,  0.9209,  0.921,  0.9211,  0.9212,  0.9213,  0.9214,  0.9215,  0.9216,  0.9217,  0.9218,  0.9219,  0.922,  0.9221,  0.9222,  0.9223,  0.9224,  0.9225,  0.9226,  0.9227,  0.9228,  0.9229,  0.923,  0.9231,  0.9232,  0.9233,  0.9234,  0.9235,  0.9236,  0.9237,  0.9238,  0.9239,  0.924,  0.9241,  0.9242,  0.9243,  0.9244,  0.9245,  0.9246,  0.9247,  0.9248,  0.9249,  0.925,  0.9251,  0.9252,  0.9253,  0.9254,  0.9255,  0.9256,  0.9257,  0.9258,  0.9259,  0.926,  0.9261,  0.9262,  0.9263,  0.9264,  0.9265,  0.9266,  0.9267,  0.9268,  0.9269,  0.927,  0.9271,  0.9272,  0.9273,  0.9274,  0.9275,  0.9276,  0.9277,  0.9278,  0.9279,  0.928,  0.9281,  0.9282,  0.9283,  0.9284,  0.9285,  0.9286,  0.9287,  0.9288,  0.9289,  0.929,  0.9291,  0.9292,  0.9293,  0.9294,  0.9295,  0.9296,  0.9297,  0.9298,  0.9299,  0.93,  0.9301,  0.9302,  0.9303,  0.9304,  0.9305,  0.9306,  0.9307,  0.9308,  0.9309,  0.931,  0.9311,  0.9312,  0.9313,  0.9314,  0.9315,  0.9316,  0.9317,  0.9318,  0.9319,  0.932,  0.9321,  0.9322,  0.9323,  0.9324,  0.9325,  0.9326,  0.9327,  0.9328,  0.9329,  0.933,  0.9331,  0.9332,  0.9333,  0.9334,  0.9335,  0.9336,  0.9337,  0.9338,  0.9339,  0.934,  0.9341,  0.9342,  0.9343,  0.9344,  0.9345,  0.9346,  0.9347,  0.9348,  0.9349,  0.935,  0.9351,  0.9352,  0.9353,  0.9354,  0.9355,  0.9356,  0.9357,  0.9358,  0.9359,  0.936,  0.9361,  0.9362,  0.9363,  0.9364,  0.9365,  0.9366,  0.9367,  0.9368,  0.9369,  0.937,  0.9371,  0.9372,  0.9373,  0.9374,  0.9375,  0.9376,  0.9377,  0.9378,  0.9379,  0.938,  0.9381,  0.9382,  0.9383,  0.9384,  0.9385,  0.9386,  0.9387,  0.9388,  0.9389,  0.939,  0.9391,  0.9392,  0.9393,  0.9394,  0.9395,  0.9396,  0.9397,  0.9398,  0.9399,  0.94,  0.9401,  0.9402,  0.9403,  0.9404,  0.9405,  0.9406,  0.9407,  0.9408,  0.9409,  0.941,  0.9411,  0.9412,  0.9413,  0.9414,  0.9415,  0.9416,  0.9417,  0.9418,  0.9419,  0.942,  0.9421,  0.9422,  0.9423,  0.9424,  0.9425,  0.9426,  0.9427,  0.9428,  0.9429,  0.943,  0.9431,  0.9432,  0.9433,  0.9434,  0.9435,  0.9436,  0.9437,  0.9438,  0.9439,  0.944,  0.9441,  0.9442,  0.9443,  0.9444,  0.9445,  0.9446,  0.9447,  0.9448,  0.9449,  0.945,  0.9451,  0.9452,  0.9453,  0.9454,  0.9455,  0.9456,  0.9457,  0.9458,  0.9459,  0.946,  0.9461,  0.9462,  0.9463,  0.9464,  0.9465,  0.9466,  0.9467,  0.9468,  0.9469,  0.947,  0.9471,  0.9472,  0.9473,  0.9474,  0.9475,  0.9476,  0.9477,  0.9478,  0.9479,  0.948,  0.9481,  0.9482,  0.9483,  0.9484,  0.9485,  0.9486,  0.9487,  0.9488,  0.9489,  0.949,  0.9491,  0.9492,  0.9493,  0.9494,  0.9495,  0.9496,  0.9497,  0.9498,  0.9499,  0.95,  0.9501,  0.9502,  0.9503,  0.9504,  0.9505,  0.9506,  0.9507,  0.9508,  0.9509,  0.951,  0.9511,  0.9512,  0.9513,  0.9514,  0.9515,  0.9516,  0.9517,  0.9518,  0.9519,  0.952,  0.9521,  0.9522,  0.9523,  0.9524,  0.9525,  0.9526,  0.9527,  0.9528,  0.9529,  0.953,  0.9531,  0.9532,  0.9533,  0.9534,  0.9535,  0.9536,  0.9537,  0.9538,  0.9539,  0.954,  0.9541,  0.9542,  0.9543,  0.9544,  0.9545,  0.9546,  0.9547,  0.9548,  0.9549,  0.955,  0.9551,  0.9552,  0.9553,  0.9554,  0.9555,  0.9556,  0.9557,  0.9558,  0.9559,  0.956,  0.9561,  0.9562,  0.9563,  0.9564,  0.9565,  0.9566,  0.9567,  0.9568,  0.9569,  0.957,  0.9571,  0.9572,  0.9573,  0.9574,  0.9575,  0.9576,  0.9577,  0.9578,  0.9579,  0.958,  0.9581,  0.9582,  0.9583,  0.9584,  0.9585,  0.9586,  0.9587,  0.9588,  0.9589,  0.959,  0.9591,  0.9592,  0.9593,  0.9594,  0.9595,  0.9596,  0.9597,  0.9598,  0.9599,  0.96,  0.9601,  0.9602,  0.9603,  0.9604,  0.9605,  0.9606,  0.9607,  0.9608,  0.9609,  0.961,  0.9611,  0.9612,  0.9613,  0.9614,  0.9615,  0.9616,  0.9617,  0.9618,  0.9619,  0.962,  0.9621,  0.9622,  0.9623,  0.9624,  0.9625,  0.9626,  0.9627,  0.9628,  0.9629,  0.963,  0.9631,  0.9632,  0.9633,  0.9634,  0.9635,  0.9636,  0.9637,  0.9638,  0.9639,  0.964,  0.9641,  0.9642,  0.9643,  0.9644,  0.9645,  0.9646,  0.9647,  0.9648,  0.9649,  0.965,  0.9651,  0.9652,  0.9653,  0.9654,  0.9655,  0.9656,  0.9657,  0.9658,  0.9659,  0.966,  0.9661,  0.9662,  0.9663,  0.9664,  0.9665,  0.9666,  0.9667,  0.9668,  0.9669,  0.967,  0.9671,  0.9672,  0.9673,  0.9674,  0.9675,  0.9676,  0.9677,  0.9678,  0.9679,  0.968,  0.9681,  0.9682,  0.9683,  0.9684,  0.9685,  0.9686,  0.9687,  0.9688,  0.9689,  0.969,  0.9691,  0.9692,  0.9693,  0.9694,  0.9695,  0.9696,  0.9697,  0.9698,  0.9699,  0.97,  0.9701,  0.9702,  0.9703,  0.9704,  0.9705,  0.9706,  0.9707,  0.9708,  0.9709,  0.971,  0.9711,  0.9712,  0.9713,  0.9714,  0.9715,  0.9716,  0.9717,  0.9718,  0.9719,  0.972,  0.9721,  0.9722,  0.9723,  0.9724,  0.9725,  0.9726,  0.9727,  0.9728,  0.9729,  0.973,  0.9731,  0.9732,  0.9733,  0.9734,  0.9735,  0.9736,  0.9737,  0.9738,  0.9739,  0.974,  0.9741,  0.9742,  0.9743,  0.9744,  0.9745,  0.9746,  0.9747,  0.9748,  0.9749,  0.975,  0.9751,  0.9752,  0.9753,  0.9754,  0.9755,  0.9756,  0.9757,  0.9758,  0.9759,  0.976,  0.9761,  0.9762,  0.9763,  0.9764,  0.9765,  0.9766,  0.9767,  0.9768,  0.9769,  0.977,  0.9771,  0.9772,  0.9773,  0.9774,  0.9775,  0.9776,  0.9777,  0.9778,  0.9779,  0.978,  0.9781,  0.9782,  0.9783,  0.9784,  0.9785,  0.9786,  0.9787,  0.9788,  0.9789,  0.979,  0.9791,  0.9792,  0.9793,  0.9794,  0.9795,  0.9796,  0.9797,  0.9798,  0.9799,  0.98,  0.9801,  0.9802,  0.9803,  0.9804,  0.9805,  0.9806,  0.9807,  0.9808,  0.9809,  0.981,  0.9811,  0.9812,  0.9813,  0.9814,  0.9815,  0.9816,  0.9817,  0.9818,  0.9819,  0.982,  0.9821,  0.9822,  0.9823,  0.9824,  0.9825,  0.9826,  0.9827,  0.9828,  0.9829,  0.983,  0.9831,  0.9832,  0.9833,  0.9834,  0.9835,  0.9836,  0.9837,  0.9838,  0.9839,  0.984,  0.9841,  0.9842,  0.9843,  0.9844,  0.9845,  0.9846,  0.9847,  0.9848,  0.9849,  0.985,  0.9851,  0.9852,  0.9853,  0.9854,  0.9855,  0.9856,  0.9857,  0.9858,  0.9859,  0.986,  0.9861,  0.9862,  0.9863,  0.9864,  0.9865,  0.9866,  0.9867,  0.9868,  0.9869,  0.987,  0.9871,  0.9872,  0.9873,  0.9874,  0.9875,  0.9876,  0.9877,  0.9878,  0.9879,  0.988,  0.9881,  0.9882,  0.9883,  0.9884,  0.9885,  0.9886,  0.9887,  0.9888,  0.9889,  0.989,  0.9891,  0.9892,  0.9893,  0.9894,  0.9895,  0.9896,  0.9897,  0.9898,  0.9899,  0.99,  0.9901,  0.9902,  0.9903,  0.9904,  0.9905,  0.9906,  0.9907,  0.9908,  0.9909,  0.991,  0.9911,  0.9912,  0.9913,  0.9914,  0.9915,  0.9916,  0.9917,  0.9918,  0.9919,  0.992,  0.9921,  0.9922,  0.9923,  0.9924,  0.9925,  0.9926,  0.9927,  0.9928,  0.9929,  0.993,  0.9931,  0.9932,  0.9933,  0.9934,  0.9935,  0.9936,  0.9937,  0.9938,  0.9939,  0.994,  0.9941,  0.9942,  0.9943,  0.9944,  0.9945,  0.9946,  0.9947,  0.9948,  0.9949,  0.995,  0.9951,  0.9952,  0.9953,  0.9954,  0.9955,  0.9956,  0.9957,  0.9958,  0.9959,  0.996,  0.9961,  0.9962,  0.9963,  0.9964,  0.9965,  0.9966,  0.9967,  0.9968,  0.9969,  0.997,  0.9971,  0.9972,  0.9973,  0.9974,  0.9975,  0.9976,  0.9977,  0.9978,  0.9979,  0.998,  0.9981,  0.9982,  0.9983,  0.9984,  0.9985,  0.9986,  0.9987,  0.9988,  0.9989,  0.999,  0.9991,  0.9992,  0.9993,  0.9994,  0.9995,  0.9996,  0.9997,  0.9998,  0.9999, 1.0};
    for(int k = 0; k < trials+1; k++){
      probcuts[k] = probcuts_temp[k];
    }
  }

}



void calculate_Dsbackground_MLindep(){

  TH1F* htamu_pr = (TH1F*)ftamuDspr->Get("hTAMUcent");
  TH1F* htamu_fd = (TH1F*)ftamuDsfd->Get("hTAMUcent");

  TH1F* hfonll_pr = (TH1F*)ffonllDs->Get("hDsPhipitoKkpipred_central");
  TH1F* hfonll_fd = (TH1F*)ffonllDs->Get("hDsPhipitoKkpifromBpred_central_corr");
  TH1F* hfonllmin_pr = (TH1F*)ffonllDs->Get("hDsPhipitoKkpipred_min");
  TH1F* hfonllmin_fd = (TH1F*)ffonllDs->Get("hDsPhipitoKkpifromBpred_min_corr");
  TH1F* hfonllmax_pr = (TH1F*)ffonllDs->Get("hDsPhipitoKkpipred_max");
  TH1F* hfonllmax_fd = (TH1F*)ffonllDs->Get("hDsPhipitoKkpifromBpred_max_corr");

  TH1F* h_gen_pr = (TH1F*)feffDs->Get("h_gen_pr");
  TH1F* h_gen_fd = (TH1F*)feffDs->Get("h_gen_fd");
  TH1F* h_presel_pr = (TH1F*)feffDs->Get("h_presel_pr");
  TH1F* h_presel_fd = (TH1F*)feffDs->Get("h_presel_fd");

  TH1F* h_match_pr[nptbins];
  TH1F* h_match_fd[nptbins];
  for(int i = 0; i <  nptbins; i++){
    h_match_pr[i] = (TH1F*)feffDsmatchBs->Get(Form("h_match_pr%d%d",ptbins[i], ptbins[i+1]));
    h_match_fd[i] = (TH1F*)feffDsmatchBs->Get(Form("h_match_fd%d%d",ptbins[i], ptbins[i+1]));
  }

  TH1F* hscaling_pr = new TH1F("hscaling_pr", "Prompt: FONLL-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{pr. D_{s}} / Ev", hfonll_pr->GetNbinsX()-1, hfonll_pr->GetBinLowEdge(1), hfonll_pr->GetBinLowEdge(hfonll_pr->GetNbinsX()));
  TH1F* hscalingmin_pr = new TH1F("hscalingmin_pr", "Prompt: FONLL_{min}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{pr. D_{s}} / Ev", hfonll_pr->GetNbinsX()-1, hfonll_pr->GetBinLowEdge(1), hfonll_pr->GetBinLowEdge(hfonll_pr->GetNbinsX()));
  TH1F* hscalingmax_pr = new TH1F("hscalingmax_pr", "Prompt: FONLL_{max}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{pr. D_{s}} / Ev", hfonll_pr->GetNbinsX()-1, hfonll_pr->GetBinLowEdge(1), hfonll_pr->GetBinLowEdge(hfonll_pr->GetNbinsX()));

  TH1F* hscaling_fd = new TH1F("hscaling_fd", "Feed-down: FONLL-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{fd. D_{s}} / Ev", hfonll_fd->GetNbinsX()-1, hfonll_fd->GetBinLowEdge(1), hfonll_fd->GetBinLowEdge(hfonll_fd->GetNbinsX()));
  TH1F* hscalingmin_fd = new TH1F("hscalingmin_fd", "Feed-down: FONLL_{min}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{fd. D_{s}} / Ev", hfonll_fd->GetNbinsX()-1, hfonll_fd->GetBinLowEdge(1), hfonll_fd->GetBinLowEdge(hfonll_fd->GetNbinsX()));
  TH1F* hscalingmax_fd = new TH1F("hscalingmax_fd", "Feed-down: FONLL_{max}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{fd. D_{s}} / Ev", hfonll_fd->GetNbinsX()-1, hfonll_fd->GetBinLowEdge(1), hfonll_fd->GetBinLowEdge(hfonll_fd->GetNbinsX()));

  Double_t testNPDs12;

  for(Int_t i = 1; i <= hfonll_pr->GetNbinsX(); i++){
    Double_t pt = hfonll_pr->GetBinCenter(i);
    Double_t dpt = hfonll_pr->GetBinWidth(i);
    Double_t yfonll = hfonll_pr->GetBinContent(i);
    Double_t yfonllmin = hfonllmin_pr->GetBinContent(i);
    Double_t yfonllmax = hfonllmax_pr->GetBinContent(i);
    Double_t ytamu = htamu_pr->GetBinContent(htamu_pr->FindBin(pt));
    Double_t dy = 1;
    Double_t acceffDs;
    if(h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt)) > 0){
      acceffDs = h_presel_pr->GetBinContent(h_presel_pr->FindBin(pt)) / h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt));
    } else {
      acceffDs = 0;
    }

    Double_t signal = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonll * ytamu;
    Double_t signalmin = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonllmin * ytamu;
    Double_t signalmax = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonllmax * ytamu;
    hscaling_pr->SetBinContent(i, signal);
    hscalingmin_pr->SetBinContent(i, signalmin);
    hscalingmax_pr->SetBinContent(i, signalmax);

    if(pt >= 1 && pt < 2){
      //Test to non-prompt Ds ITS3 analysis, bin [1-2]
      //https://indico.cern.ch/event/893570/contributions/3769353/attachments/1999626/3337235/Ds_ITS3_estimates_2020_03_06.pdf#page=7
      //  Acc x Eff Bs TTree: 0.00565959
      //  Acc x Eff NP Ds: 0.002
      //    Exp Ev Fabrizio: "Indeed I always write 8, but I use 7.7"
      testNPDs12 += (0.002/0.00565959) * signal*7700000000.;
    }
  }

  for(Int_t i = 1; i <= hfonll_fd->GetNbinsX(); i++){
    Double_t pt = hfonll_fd->GetBinCenter(i);
    Double_t dpt = hfonll_fd->GetBinWidth(i);
    Double_t yfonll = hfonll_fd->GetBinContent(i);
    Double_t yfonllmin = hfonllmin_fd->GetBinContent(i);
    Double_t yfonllmax = hfonllmax_fd->GetBinContent(i);
    Double_t ytamu = htamu_fd->GetBinContent(htamu_fd->FindBin(pt));
    Double_t dy = 1;
    Double_t acceffDs;
    if(h_gen_fd->GetBinContent(h_gen_fd->FindBin(pt)) > 0){
      acceffDs = h_presel_fd->GetBinContent(h_presel_fd->FindBin(pt)) / h_gen_fd->GetBinContent(h_gen_fd->FindBin(pt));
    } else {
      acceffDs = 0;
    }

    Double_t signal = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonll * ytamu;
    Double_t signalmin = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonllmin * ytamu;
    Double_t signalmax = 2 * dpt * dy * BRDs * acceffDs * TAA2 * yfonllmax * ytamu;
    hscaling_fd->SetBinContent(i, signal);
    hscalingmin_fd->SetBinContent(i, signalmin);
    hscalingmax_fd->SetBinContent(i, signalmax);

    if(pt >= 1 && pt < 2){
      //Test to non-prompt Ds ITS3 analysis, bin [1-2]
      //https://indico.cern.ch/event/893570/contributions/3769353/attachments/1999626/3337235/Ds_ITS3_estimates_2020_03_06.pdf#page=7
      //  Acc x Eff Bs TTree: 0.00523262
      //  Acc x Eff NP Ds: 0.006
      //    Exp Ev Fabrizio: "Indeed I always write 8, but I use 7.7"
      testNPDs12 += (0.006/0.00523262) * signal*7700000000.;
    }
  }
  cout << "\nCrossCheck expected Ds signal for pT [1-2] with NP Ds analysis: " << endl;
  cout << "  NP: S = f(Sgnf, S/B) = ((70*70)*0.9 + 70*70)/0.9 = 10344.444" << endl;
  cout << "  Bs: S (scaled to match efficiency NP analysis) = " << testNPDs12 << endl;
  cout << "    Should be very close (03/05/20)!" << endl << endl;

  TH1F* h_expsig_pr = new TH1F("h_expsig_pr", ";#it{p}_{T}(B_{s}) (GeV/it{c});Expected pr. signal per event", nptbins, ptbinsfl);
  TH1F* h_expsigmin_pr = new TH1F("h_expsigmin_pr", ";#it{p}_{T}(B_{s}) (GeV/it{c});Min. Expected pr. signal per event", nptbins, ptbinsfl);
  TH1F* h_expsigmax_pr = new TH1F("h_expsigmax_pr", ";#it{p}_{T}(B_{s}) (GeV/it{c});Max. Expected pr. signal per event", nptbins, ptbinsfl);
  TH1F* h_expsig_fd = new TH1F("h_expsig_fd", ";#it{p}_{T}(B_{s}) (GeV/it{c});Expected fd. signal per event", nptbins, ptbinsfl);
  TH1F* h_expsigmin_fd = new TH1F("h_expsigmin_fd", ";#it{p}_{T}(B_{s}) (GeV/it{c});Min. Expected fd. signal per event", nptbins, ptbinsfl);
  TH1F* h_expsigmax_fd = new TH1F("h_expsigmax_fd", ";#it{p}_{T}(B_{s}) (GeV/it{c});Max. Expected fd. signal per event", nptbins, ptbinsfl);

  for(int i = 0; i < nptbins; i++){
    Double_t pt = 0.5 * (ptbins[i] + ptbins[i+1]);

    Double_t sigpr = 0;
    Double_t sigprmin = 0;
    Double_t sigprmax = 0;
    Double_t sigpr_nonmatch = 0;
    Double_t sigprmin_nonmatch = 0;
    Double_t sigprmax_nonmatch = 0;
    for(int j = 1; j <= h_gen_pr->GetNbinsX(); j++){
      sigpr += h_match_pr[i]->GetBinContent(j) * hscaling_pr->GetBinContent(j);
      sigprmin += h_match_pr[i]->GetBinContent(j) * hscalingmin_pr->GetBinContent(j);
      sigprmax += h_match_pr[i]->GetBinContent(j) * hscalingmax_pr->GetBinContent(j);

      if(hscaling_pr->GetBinCenter(j) >= ptbins[i] && hscaling_pr->GetBinCenter(j) < ptbins[i+1]){
        sigpr_nonmatch += hscaling_pr->GetBinContent(j);
        sigprmin_nonmatch += hscalingmin_pr->GetBinContent(j);
        sigprmax_nonmatch += hscalingmax_pr->GetBinContent(j);
      }
    }
    h_expsig_pr->SetBinContent(h_expsig_pr->FindBin(pt), sigpr);
    h_expsigmin_pr->SetBinContent(h_expsigmin_pr->FindBin(pt), sigprmin);
    h_expsigmax_pr->SetBinContent(h_expsigmax_pr->FindBin(pt), sigprmax);

    Double_t sigfd = 0;
    Double_t sigfdmin = 0;
    Double_t sigfdmax = 0;
    Double_t sigfd_nonmatch = 0;
    Double_t sigfdmin_nonmatch = 0;
    Double_t sigfdmax_nonmatch = 0;
    for(int j = 1; j <= h_gen_fd->GetNbinsX(); j++){
      sigfd += h_match_fd[i]->GetBinContent(j) * hscaling_fd->GetBinContent(j);
      sigfdmin += h_match_fd[i]->GetBinContent(j) * hscalingmin_fd->GetBinContent(j);
      sigfdmax += h_match_fd[i]->GetBinContent(j) * hscalingmax_fd->GetBinContent(j);

      if(hscaling_fd->GetBinCenter(j) >= ptbins[i] && hscaling_fd->GetBinCenter(j) < ptbins[i+1]){
        sigfd_nonmatch += hscaling_fd->GetBinContent(j);
        sigfdmin_nonmatch += hscalingmin_fd->GetBinContent(j);
        sigfdmax_nonmatch += hscalingmax_fd->GetBinContent(j);
      }
    }
    h_expsig_fd->SetBinContent(h_expsig_fd->FindBin(pt), sigfd);
    h_expsigmin_fd->SetBinContent(h_expsigmin_fd->FindBin(pt), sigfdmin);
    h_expsigmax_fd->SetBinContent(h_expsigmax_fd->FindBin(pt), sigfdmax);
  }

  fExpDsSignal->cd();
  h_expsig_pr->Write();
  h_expsigmin_pr->Write();
  h_expsigmax_pr->Write();
  h_expsig_fd->Write();
  h_expsigmin_fd->Write();
  h_expsigmax_fd->Write();

}
