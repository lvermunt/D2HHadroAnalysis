R__LOAD_LIBRARY(aliphysics/AliHFInvMassMultiTrialFit1_cxx.so)
R__LOAD_LIBRARY(aliphysics/AliHFInvMassFitter1_cxx.so)

enum{kParaCent, kParaLow, kParaUp, kFit, kFit2};
TString optname[5] = {"kParaCent", "kFit", "kParaLow", "kParaUp", "kFit2"};

const int nptbins = 6;
Int_t ptlimits[nptbins+1] = {0, 2, 4, 8, 12, 16, 24};

void MultiTrialFit(TH1D* hmass, Int_t ptbin, Int_t poisson, Int_t opt, Double_t signal, Double_t width, Double_t mean);
Double_t AliHFFitter(TH1F* hmass, Int_t bkgfunc, Double_t mass, Double_t sigma);
void SetStyleHisto(TH1F *h);

void systematic_rawyield(TString inputfile, Int_t ptbin = 1, Int_t opt = kParaCent){
  TFile* fin = new TFile(inputfile.Data());

  TF1* fbkg = (TF1*)fin->Get(Form("fbkgscaled_%d",ptbin));
  //TF1* fbkg2 = (TF1*)fin->Get(Form("fbkg2scaled_%d",ptbin));
  TF1* fcent = (TF1*)fin->Get(Form("fcentscaled_%d",ptbin));
  TF1* fup = (TF1*)fin->Get(Form("fupscaled_%d",ptbin));
  TF1* flow = (TF1*)fin->Get(Form("flowscaled_%d",ptbin));
  TF1* fdspr1 = (TF1*)fin->Get(Form("fdspr1scaled_%d",ptbin));
  TF1* fdsfdbzero1 = (TF1*)fin->Get(Form("fdsfdbzero1scaled_%d",ptbin));
  TF1* fdsfdbplus1 = (TF1*)fin->Get(Form("fdsfdbplus1scaled_%d",ptbin));
  TF1* fdsfdbs1 = (TF1*)fin->Get(Form("fdsfdbs1scaled_%d",ptbin));
  TF1* fdsfdlambdab1 = (TF1*)fin->Get(Form("fdsfdlambdab1scaled_%d",ptbin));
  TF1* fsig = (TF1*)fin->Get(Form("fsigfull_%d",ptbin));

  fbkg->SetName("fbkg");
  //fbkg2->SetName("fbkg2");
  fcent->SetName("fcent");
  fup->SetName("fup");
  flow->SetName("flow");
  fdspr1->SetName("fdspr1");
  fdsfdbzero1->SetName("fdsfdbzero1");
  fdsfdbplus1->SetName("fdsfdbplus1");
  fdsfdbs1->SetName("fdsfdbs1");
  fdsfdlambdab1->SetName("fdsfdlambdab1");
  fsig->SetName("fsig");

  TH1F* hmass = (TH1F*)fin->Get(Form("histo_invmass_%d",ptbin));

  TH1F* hmass2 = (TH1F*)hmass->Clone(Form("histo_invmass2_%d",ptbin));
  hmass2->Reset("ICEMS");
  hmass2->FillRandom("fcent", 10 * hmass->Integral(hmass->FindBin(fcent->GetXmin()), hmass->FindBin(fcent->GetXmax())));

  TF1* fnewexpo = new TF1("fnewexpo","expo",fcent->GetXmin(), fcent->GetXmax());
  hmass2->Fit("fnewexpo", "R,E,+");
  TF1* fnewpol2;
  if(ptbin < 3) fnewpol2 = new TF1("fnewpol2","pol2",fcent->GetXmin(), fcent->GetXmax());
  else fnewpol2 = new TF1("fnewpol2","pol1",fcent->GetXmin(), fcent->GetXmax());
  hmass2->Fit("fnewpol2", "R,E,+");
  hmass->Draw("same ep");

  Double_t signalexp = (1./ hmass->GetBinWidth(1)) * fsig->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1));
  Double_t meanexp = fsig->GetParameter(1);
  Double_t sigmaexp = fsig->GetParameter(2);

  Int_t ntrials = 25;
  TH1F* hRMSfit = new TH1F("hRMSfit",Form("%d < #it{p}_{T} < %d GeV/#it{c} (fit);trial;RMS (%%)",ptlimits[ptbin], ptlimits[ptbin+1]),3*ntrials,-0.5,(3*ntrials)-0.5);
  TH1F* hMeanfit = new TH1F("hMeanfit",Form("%d < #it{p}_{T} < %d GeV/#it{c} (fit);trial;Mean",ptlimits[ptbin], ptlimits[ptbin+1]),3*ntrials,-0.5,(3*ntrials)-0.5);
  TH1F* hMeanfit2 = new TH1F("hMeanfit2",Form("%d < #it{p}_{T} < %d GeV/#it{c} (fit);Mean MultiFitter;Entries",ptlimits[ptbin], ptlimits[ptbin+1]),(int)signalexp,0,2*signalexp);
  TH1F* hRMSbc = new TH1F("hRMSbc",Form("%d < #it{p}_{T} < %d GeV/#it{c} (BC 3#sigma);trial;RMS (%%)",ptlimits[ptbin], ptlimits[ptbin+1]),3*ntrials,-0.5,(3*ntrials)-0.5);
  SetStyleHisto(hRMSfit);
  SetStyleHisto(hRMSbc);

  for(opt = kParaCent; opt <= kParaUp; opt++){
    for(int j = 0; j < ntrials; j++){
      hmass->Reset("ICEMS");
      if(opt == kParaCent) hmass->FillRandom("fcent",(1./ hmass->GetBinWidth(1)) * fcent->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      //if(opt == kFit) hmass->FillRandom("fbkg",(1./ hmass->GetBinWidth(1)) * fbkg->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      //if(opt == kFit) hmass->FillRandom("fcent",(1./ hmass->GetBinWidth(1)) * fcent->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      //if(opt == kParaLow) hmass->FillRandom("flow",(1./ hmass->GetBinWidth(1)) * flow->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      //if(opt == kParaUp) hmass->FillRandom("fup",(1./ hmass->GetBinWidth(1)) * fup->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      if(opt == kParaLow) hmass->FillRandom("fnewexpo",(1./ hmass->GetBinWidth(1)) * fcent->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      if(opt == kParaUp) hmass->FillRandom("fnewpol2",(1./ hmass->GetBinWidth(1)) * fcent->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));

      hmass->FillRandom("fdspr1",(1./ hmass->GetBinWidth(1)) * fdspr1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      hmass->FillRandom("fdsfdbzero1",(1./ hmass->GetBinWidth(1)) * fdsfdbzero1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      hmass->FillRandom("fdsfdbplus1",(1./ hmass->GetBinWidth(1)) * fdsfdbplus1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      hmass->FillRandom("fdsfdbs1",(1./ hmass->GetBinWidth(1)) * fdsfdbs1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      hmass->FillRandom("fdsfdlambdab1",(1./ hmass->GetBinWidth(1)) * fdsfdlambdab1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
      hmass->FillRandom("fsig",(1./ hmass->GetBinWidth(1)) * fsig->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));

      MultiTrialFit((TH1D*)hmass, ptbin, j, opt, signalexp, sigmaexp, meanexp);

      TCanvas* ctest = new TCanvas();
      ctest->cd();
      TH1F* hmass_plot = (TH1F*)hmass->Clone(Form("hmass_plot_%d",j));
      hmass_plot->Draw("same ep");
    }
  }

  for(opt = kParaCent; opt <= kParaUp; opt++){
    for(int j = 0; j < ntrials; j++){
      TFile* fmulti = new TFile(Form("aliphysics/multitrialBs/Output_Multitrial_%d_%d_Poission%dOpt%s.root",ptlimits[ptbin],ptlimits[ptbin+1],j,optname[opt].Data()));
      TH1F* fHistoRawYieldDistAll = (TH1F*)fmulti->Get("hRawYieldDistAll");
      TH1F* fHistoRawYieldDistBinC1All = (TH1F*)fmulti->Get("hRawYieldDistBinC1All");

      if(fHistoRawYieldDistAll->GetMean() > 0){
        hRMSfit->SetBinContent(opt*ntrials + j+1, (fHistoRawYieldDistAll->GetRMS())/fHistoRawYieldDistAll->GetMean() * 100);
        hRMSfit->SetBinError(opt*ntrials + j+1, (fHistoRawYieldDistAll->GetRMSError())/fHistoRawYieldDistAll->GetMean() * 100);

        hMeanfit->SetBinContent(opt*ntrials + j+1, fHistoRawYieldDistAll->GetMean());
        hMeanfit->SetBinError(opt*ntrials + j+1, fHistoRawYieldDistAll->GetMeanError());
        hMeanfit2->Fill(fHistoRawYieldDistAll->GetMean());

        hRMSbc->SetBinContent(opt*ntrials + j+1, (fHistoRawYieldDistBinC1All->GetRMS())/fHistoRawYieldDistBinC1All->GetMean() * 100);
        hRMSbc->SetBinError(opt*ntrials + j+1, (fHistoRawYieldDistBinC1All->GetRMSError())/fHistoRawYieldDistBinC1All->GetMean() * 100);
      }
    }
  }
  hRMSfit->SetLineColor(kBlue+1);
  hMeanfit->SetLineColor(kRed+1);
  hMeanfit2->SetLineColor(kRed+1);
  hRMSbc->SetLineColor(kGreen+2);
  hRMSfit->SetMarkerColor(kBlue+1);
  hMeanfit->SetMarkerColor(kRed+1);
  hMeanfit2->SetMarkerColor(kRed+1);
  hRMSbc->SetMarkerColor(kGreen+2);
  hRMSfit->SetMarkerStyle(20);
  hMeanfit->SetMarkerStyle(20);
  hMeanfit2->SetMarkerStyle(20);
  hRMSbc->SetMarkerStyle(20);

  TLine* l1 = new TLine(ntrials-0.5,0,ntrials-0.5,20);
  TLine* l2 = new TLine(2*ntrials-0.5,0,2*ntrials-0.5,20);
  TLine* l3 = new TLine(3*ntrials-0.5,0,3*ntrials-0.5,20);
  l1->SetLineStyle(2); l2->SetLineStyle(2); l3->SetLineStyle(2);

  TLine* l1_2 = new TLine(ntrials-0.5,0,ntrials-0.5,2*signalexp);
  TLine* l2_2 = new TLine(2*ntrials-0.5,0,2*ntrials-0.5,2*signalexp);
  TLine* l3_2 = new TLine(3*ntrials-0.5,0,3*ntrials-0.5,2*signalexp);
  l1_2->SetLineStyle(2); l2_2->SetLineStyle(2); l3_2->SetLineStyle(2);

  TCanvas* csummary = new TCanvas("csummary","",800,400);
  csummary->Divide(2);
  csummary->cd(1);
  hRMSfit->Draw("ep");
  hRMSfit->GetYaxis()->SetRangeUser(0,20);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  l1->Draw();
  l2->Draw();
  //l3->Draw();
  csummary->cd(2);
  hMeanfit->Draw("ep");
  hMeanfit->GetYaxis()->SetRangeUser(0,2*signalexp);
  //hRMSbc->Draw("ep");
  //hRMSbc->GetYaxis()->SetRangeUser(0,20);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  l1_2->Draw();
  l2_2->Draw();
  //l3_2->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18); info1.SetTextColor(kRed+1);
  info1.DrawLatex(0.25, 0.8, Form("RMS = %.1f (%.1f%%)",hMeanfit2->GetRMS(), 100*hMeanfit2->GetRMS()/hMeanfit2->GetMean()));

  csummary->SaveAs(Form("rawyielsystematic_Bs_ITS2_v4_ptbin%d%d.eps",ptlimits[ptbin],ptlimits[ptbin+1]));
}

void MultiTrialFit(TH1D* hmass, Int_t ptbin, Int_t poisson, Int_t opt, Double_t signal, Double_t width, Double_t mean=5.366){

  TGaxis::SetMaxDigits(3);
  Int_t ptlimsLow = ptlimits[ptbin];
  Int_t ptlimsHigh = ptlimits[ptbin+1];

  Double_t signalRef = signal;
  Double_t sigmaMC = width;
  Double_t mass = mean;

  AliHFInvMassMultiTrialFit1* MultiTrialFit = new AliHFInvMassMultiTrialFit1();

  MultiTrialFit->AddInvMassFitSaveAsFormat("png");
  MultiTrialFit->SetDrawIndividualFits(kFALSE);

  Double_t minMassStep[4]={5.066,5.116,5.166,5.216};
  Double_t maxMassStep[4]={5.666,5.616,5.566,5.516};
  MultiTrialFit->ConfigureLowLimFitSteps(4,minMassStep);
  MultiTrialFit->ConfigureUpLimFitSteps(4,maxMassStep);

  Double_t nSigmasBC[1]={3.0};
  MultiTrialFit->ConfigurenSigmaBinCSteps(1,nSigmasBC);

  Int_t rebinStep[1] = {1};
  MultiTrialFit->ConfigureRebinSteps(1,rebinStep);
  MultiTrialFit->SetNumOfFirstBinSteps(rebinStep[0]);

  MultiTrialFit->SetUseDoubleGausSignal(kFALSE);
  MultiTrialFit->SetMass( mass );
  MultiTrialFit->SetSigmaGaussMC( sigmaMC );
  //Use Default LL fit
  //MultiTrialFit->SetUseChi2Fit();

  MultiTrialFit->SetUseExpoBackground(kTRUE);
  MultiTrialFit->SetUseLinBackground(kTRUE);//kFALSE);
  MultiTrialFit->SetUsePol2Background(kTRUE);
  MultiTrialFit->SetUsePol3Background(kFALSE);
  MultiTrialFit->SetUsePol4Background(kFALSE);
  MultiTrialFit->SetUsePol5Background(kFALSE);
  MultiTrialFit->SetUsePowerLawBackground(kFALSE);
  MultiTrialFit->SetUsePowerLawTimesExpoBackground(kFALSE);

  MultiTrialFit->SetSigmaMCVariation(0.15);
  MultiTrialFit->SetUseFixSigUpFreeMean(kFALSE);
  MultiTrialFit->SetUseFixSigDownFreeMean(kFALSE);
  MultiTrialFit->SetUseFreeS(kFALSE);//kTRUE);
  MultiTrialFit->SetUseFixSigFreeMean(kTRUE);
  MultiTrialFit->SetUseFixedMeanFreeS(kFALSE);//kTRUE);
  MultiTrialFit->SetUseFixSigFixMean(kTRUE);

  MultiTrialFit->SetpTMin(ptlimsLow);
  MultiTrialFit->SetpTMax(ptlimsHigh);
  MultiTrialFit->SetMaxChi2(2.);
  MultiTrialFit->SetSigmaRef( sigmaMC );
  MultiTrialFit->SetRawYieldRef( signalRef );

  TCanvas* cFit = new TCanvas();
  TPad* pad = (TPad*) cFit->cd();
  MultiTrialFit->DoMultiTrials(hmass, pad);

  TString outputname = Form("aliphysics/multitrialBs/Output_Multitrial_%d_%d_Poission%dOpt%s.root",ptlimsLow,ptlimsHigh,poisson,optname[opt].Data());
  MultiTrialFit->SaveToRoot(outputname.Data());

  TCanvas* c1 = new TCanvas(Form("c1_%d_%d_Poission%dOpt%s",ptlimsLow,ptlimsHigh,poisson,optname[opt].Data()),"",10,10,960,540);
  MultiTrialFit->DrawHistos(c1);

  c1->SaveAs(Form("aliphysics/multitrialBs/MultiFit_Bs_%d_%d_Poission%dOpt%s.eps",ptlimsLow,ptlimsHigh,poisson,optname[opt].Data()));
  TFile* foutput = new TFile(outputname.Data(),"UPDATE");
  foutput->cd();
  c1->Write();
  foutput->Close();
}

void systematic_mean_rawyield(TString inputfile, Int_t ptbin = 1, Int_t its = 3){
  TFile* fin = new TFile(inputfile.Data());

  TF1* fbkg = (TF1*)fin->Get(Form("fbkgscaled_%d",ptbin));
  //TF1* fbkg2 = (TF1*)fin->Get(Form("fbkg2scaled_%d",ptbin));
  TF1* fcent = (TF1*)fin->Get(Form("fcentscaled_%d",ptbin));
  TF1* fup = (TF1*)fin->Get(Form("fupscaled_%d",ptbin));
  TF1* flow = (TF1*)fin->Get(Form("flowscaled_%d",ptbin));
  TF1* fdspr1 = (TF1*)fin->Get(Form("fdspr1scaled_%d",ptbin));
  TF1* fdsfdbzero1 = (TF1*)fin->Get(Form("fdsfdbzero1scaled_%d",ptbin));
  TF1* fdsfdbplus1 = (TF1*)fin->Get(Form("fdsfdbplus1scaled_%d",ptbin));
  TF1* fdsfdbs1 = (TF1*)fin->Get(Form("fdsfdbs1scaled_%d",ptbin));
  TF1* fdsfdlambdab1 = (TF1*)fin->Get(Form("fdsfdlambdab1scaled_%d",ptbin));
  TF1* fsig = (TF1*)fin->Get(Form("fsigfull_%d",ptbin));

  fbkg->SetName("fbkg");
  //fbkg2->SetName("fbkg2");
  fcent->SetName("fcent");
  fup->SetName("fup");
  flow->SetName("flow");
  fdspr1->SetName("fdspr1");
  fdsfdbzero1->SetName("fdsfdbzero1");
  fdsfdbplus1->SetName("fdsfdbplus1");
  fdsfdbs1->SetName("fdsfdbs1");
  fdsfdlambdab1->SetName("fdsfdlambdab1");
  fsig->SetName("fsig");

  Bool_t is_expo = kTRUE;
  if(fcent->GetNpar() < 3) is_expo = kFALSE;

  TH1F* hmass = (TH1F*)fin->Get(Form("histo_invmass_%d",ptbin));

  Double_t signalexp = (1./ hmass->GetBinWidth(1)) * fsig->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1));
  Double_t meanexp = fsig->GetParameter(1);
  Double_t sigmaexp = fsig->GetParameter(2);

  Int_t ntrials = 1000;
  TH1F* hMeanfit0 = new TH1F("hMeanfit0",Form("%d < #it{p}_{T} < %d GeV/#it{c} ITS%d (expo fit);trial;Yield",ptlimits[ptbin], ptlimits[ptbin+1], its),ntrials,-0.5,ntrials-0.5);
  TH1F* hMeanfit2 = new TH1F("hMeanfit2",Form("%d < #it{p}_{T} < %d GeV/#it{c} ITS%d (pol fit);trial;Yield",ptlimits[ptbin], ptlimits[ptbin+1], its),ntrials,-0.5,ntrials-0.5);
  TH1F* hResidual = new TH1F("hResidual",Form("%d < #it{p}_{T} < %d GeV/#it{c} ITS%d (expo - pol);trial;Residual",ptlimits[ptbin], ptlimits[ptbin+1], its),ntrials,-0.5,ntrials-0.5);
  TH1F* hResidual2 = new TH1F("hResidual2",Form("%d < #it{p}_{T} < %d GeV/#it{c} ITS%d (expo - pol);Residual (%%);Entries",ptlimits[ptbin], ptlimits[ptbin+1], its),100, -40, 40);
  TH1F* hResidual3 = new TH1F("hResidual3",Form("%d < #it{p}_{T} < %d GeV/#it{c} ITS%d (expo - pol);Residual;Entries",ptlimits[ptbin], ptlimits[ptbin+1], its),2000, -1000, 1000);

  SetStyleHisto(hMeanfit0);
  SetStyleHisto(hMeanfit2);
  SetStyleHisto(hResidual);
  SetStyleHisto(hResidual2);

  for(int j = 0; j < ntrials; j++){
    hmass->Reset("ICEMS");
    hmass->FillRandom("fcent",(1./ hmass->GetBinWidth(1)) * fcent->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
    hmass->FillRandom("fdspr1",(1./ hmass->GetBinWidth(1)) * fdspr1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
    hmass->FillRandom("fdsfdbzero1",(1./ hmass->GetBinWidth(1)) * fdsfdbzero1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
    hmass->FillRandom("fdsfdbplus1",(1./ hmass->GetBinWidth(1)) * fdsfdbplus1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
    hmass->FillRandom("fdsfdbs1",(1./ hmass->GetBinWidth(1)) * fdsfdbs1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
    hmass->FillRandom("fdsfdlambdab1",(1./ hmass->GetBinWidth(1)) * fdsfdlambdab1->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));
    hmass->FillRandom("fsig",(1./ hmass->GetBinWidth(1)) * fsig->Integral(hmass->GetBinLowEdge(1), hmass->GetBinLowEdge(hmass->GetNbinsX()+1)));

    Double_t yield0 = AliHFFitter(hmass, 0, meanexp, sigmaexp);
    Double_t yield2 = AliHFFitter(hmass, 2, meanexp, sigmaexp);

    hMeanfit0->SetBinContent(j+1, yield0);
    hMeanfit2->SetBinContent(j+1, yield2);
    hMeanfit0->SetBinError(j+1, 0.0001);
    hMeanfit2->SetBinError(j+1, 0.0001);

    hResidual->SetBinContent(j+1, yield0 - yield2);
    hResidual->SetBinError(j+1, 0.0001);
    hResidual3->Fill(yield0 - yield2);

    Double_t yieldref = yield0;
    if(!is_expo) yieldref = yield2;
    hResidual2->Fill(100*(yield0 - yield2)/yieldref);
  }
  hMeanfit0->SetLineColor(kBlue+1);
  hMeanfit2->SetLineColor(kGreen+2);
  hResidual->SetLineColor(kRed+1);
  hResidual2->SetLineColor(kRed+1);
  hMeanfit0->SetMarkerColor(kBlue+1);
  hMeanfit2->SetMarkerColor(kGreen+2);
  hResidual->SetMarkerColor(kRed+1);
  hResidual2->SetMarkerColor(kRed+1);
  hMeanfit0->SetMarkerStyle(20);
  hMeanfit2->SetMarkerStyle(20);
  hResidual->SetMarkerStyle(20);
  hResidual2->SetMarkerStyle(20);

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18); info1.SetTextColor(kBlack);

  TCanvas* csummary = new TCanvas("csummary","",800,800);
  csummary->Divide(2,2);
  csummary->cd(1);
  hMeanfit0->Draw("ep");
  hMeanfit0->GetYaxis()->SetRangeUser(0,2*signalexp);
  hMeanfit0->GetYaxis()->SetMaxDigits(2);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  csummary->cd(2);
  hMeanfit2->Draw("ep");
  hMeanfit2->GetYaxis()->SetRangeUser(0,2*signalexp);
  hMeanfit2->GetYaxis()->SetMaxDigits(2);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  csummary->cd(3);
  hResidual->Draw("ep");
  hResidual->GetYaxis()->SetRangeUser(-0.5*signalexp,0.5*signalexp);
  hResidual->GetYaxis()->SetMaxDigits(2);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  info1.DrawLatex(0.14, 0.82, Form("Mean = %.1f", hResidual3->GetMean()));
  info1.DrawLatex(0.14, 0.76, Form("RMS = %.1f", hResidual3->GetRMS()));
  csummary->cd(4);
  hResidual2->Draw("hist");
  hResidual2->GetYaxis()->SetMaxDigits(2);
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  info1.DrawLatex(0.14, 0.82, Form("Mean = %.1f%%", hResidual2->GetMean()));
  info1.DrawLatex(0.14, 0.76, Form("RMS = %.1f%%", hResidual2->GetRMS()));


  csummary->SaveAs(Form("tempfig/mean_yieldsystematic_Bs_ITS%d_ptbin%d%d.eps",ptlimits[ptbin],ptlimits[ptbin+1],its));
}

Double_t AliHFFitter(TH1F* hmass, Int_t bkgfunc, Double_t mass, Double_t sigma){

  AliHFInvMassFitter1 *fitter = new AliHFInvMassFitter1(hmass, 5.066, 5.666, bkgfunc, 0);
  fitter->SetUseLikelihoodFit();
  fitter->SetInitialGaussianMean(mass);
  fitter->SetFixGaussianSigma(sigma);
  fitter->MassFitter(0);

  return fitter->GetRawYield();
}

void SetStyleHisto(TH1F *h){

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
  h->GetXaxis()->SetTitleOffset(0.88);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetNdivisions(505);

}
