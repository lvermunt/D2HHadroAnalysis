/**************************************************************************
 * Copyright(c) 2008-2019, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TMath.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TPaveText.h>
#include "AliHFInvMassFitter1.h"
#include "AliHFInvMassMultiTrialFit1.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFInvMassMultiTrialFit1);
/// \endcond


//_________________________________________________________________________
AliHFInvMassMultiTrialFit1::AliHFInvMassMultiTrialFit1() :
TNamed(),
fNumOfRebinSteps(4),
fRebinSteps(0x0),
fNumOfFirstBinSteps(1),
fNumOfLowLimFitSteps(6),
fLowLimFitSteps(0x0),
fNumOfUpLimFitSteps(6),
fUpLimFitSteps(0x0),
fNumOfnSigmaBinCSteps(11),
fnSigmaBinCSteps(0x0),
fnSigmaForBkgEval(3),
fSigmaGausMC(0.010),
fSigmaMCVariation(0.15),
fMassD(1.86484),
fSuffix(""),
fFitOption(0),
fUseExpoBkg(kTRUE),
fUseLinBkg(kTRUE),
fUsePol2Bkg(kTRUE),
fUsePol3Bkg(kTRUE),
fUsePol4Bkg(kTRUE),
fUsePol5Bkg(kFALSE),
fUsePowLawBkg(kFALSE),
fUsePowLawTimesExpoBkg(kFALSE),
fUseNoBackgroundOnlySignal(kFALSE),
fUseFixSigUpFreeMean(kTRUE),
fUseFixSigDownFreeMean(kTRUE),
fUseFreeS(kTRUE),
fUseFixedMeanFreeS(kTRUE),
fUseFixSigFreeMean(kTRUE),
fUseFixSigFixMean(kTRUE),
fUseSecondPeak(kFALSE),
fMassSecondPeak(1.86958),
fSigmaSecondPeak(0.01),
fFixMassSecondPeak(kFALSE),
fFixSigmaSecondPeak(kFALSE),
fRangeGaus2sigma(kFALSE),
fSaveBkgVal(kFALSE),
fDrawIndividualFits(kFALSE),
fHistoRawYieldDistAll(0x0),
fHistoRawYieldTrialAll(0x0),
fHistoSigmaTrialAll(0x0),
fHistoMeanTrialAll(0x0),
fHistoChi2TrialAll(0x0),
fHistoSignifTrialAll(0x0),
fHistoBkgTrialAll(0x0),
fHistoBkgInBinEdgesTrialAll(0x0),
fHistoRawYieldDistBinC0All(0x0),
fHistoRawYieldTrialBinC0All(0x0),
fHistoRawYieldDistBinC1All(0x0),
fHistoRawYieldTrialBinC1All(0x0),
fHistoRawYieldDistBinC0All_2(0x0),
fHistoRawYieldTrialBinC0All_2(0x0),
fHistoRawYieldDistBinC1All_2(0x0),
fHistoRawYieldTrialBinC1All_2(0x0),
fHistoRawYieldDist(0x0),
fHistoRawYieldTrial(0x0),
fHistoSigmaTrial(0x0),
fHistoMeanTrial(0x0),
fHistoChi2Trial(0x0),
fHistoSignifTrial(0x0),
fHistoBkgTrial(0x0),
fHistoBkgInBinEdgesTrial(0x0),
fHistoRawYieldDistBinC0(0x0),
fHistoRawYieldTrialBinC0(0x0),
fHistoRawYieldDistBinC1(0x0),
fHistoRawYieldTrialBinC1(0x0),
fHistoRawYieldDistBinC0_2(0x0),
fHistoRawYieldTrialBinC0_2(0x0),
fHistoRawYieldDistBinC1_2(0x0),
fHistoRawYieldTrialBinC1_2(0x0),
fhTemplRefl(0x0),
fhTemplSign(0x0),
fFixRefloS(1.),
fNtupleMultiTrials(0x0),
fMinYieldGlob(0),
fMaxYieldGlob(0),
fMinYieldGlobBC(0),
fMaxYieldGlobBC(0),
fMinYieldGlobBC_2(0),
fMaxYieldGlobBC_2(0),
fMassFitters()
{
  // constructor
  Int_t rebinStep[4]={3,4,5,6};
  Double_t minMassStep[6]={1.68,1.70,1.72,1.74,1.76,1.78};
  Double_t maxMassStep[6]={2.06,2.04,2.02,2.00,1.98,1.96};
  Double_t nSigmasBC[11]={2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0};
  ConfigureRebinSteps(4,rebinStep);
  ConfigureLowLimFitSteps(6,minMassStep);
  ConfigureUpLimFitSteps(6,maxMassStep);
  ConfigurenSigmaBinCSteps(11,nSigmasBC);
}

//________________________________________________________________________
AliHFInvMassMultiTrialFit1::~AliHFInvMassMultiTrialFit1(){
  // destructor
  delete [] fRebinSteps;
  delete [] fLowLimFitSteps;
  delete [] fUpLimFitSteps;
  if(fhTemplRefl) delete fhTemplRefl;
  if(fhTemplSign) delete fhTemplSign;
  for (auto fitter : fMassFitters) delete fitter;
}

//________________________________________________________________________
Bool_t AliHFInvMassMultiTrialFit1::CreateHistos(){
  // creates output histograms
  
  Int_t kNBkgFuncCasesTurnedOn = 0;
  for(Int_t typeb=0; typeb<kNBkgFuncCases; typeb++){
    if(typeb==kExpoBkg && !fUseExpoBkg) continue;
    if(typeb==kLinBkg && !fUseLinBkg) continue;
    if(typeb==kPol2Bkg && !fUsePol2Bkg) continue;
    if(typeb==kPol3Bkg && !fUsePol3Bkg) continue;
    if(typeb==kPol4Bkg && !fUsePol4Bkg) continue;
    if(typeb==kPol5Bkg && !fUsePol5Bkg) continue;
    if(typeb==kPowBkg && !fUsePowLawBkg) continue;
    if(typeb==kPowTimesExpoBkg && !fUsePowLawTimesExpoBkg) continue;
    if(typeb==kNoBkgOnlySig && !fUseNoBackgroundOnlySignal) continue;
    kNBkgFuncCasesTurnedOn++;
  }
  Int_t kNFitConfCasesTurnedOn = 0;
  for(Int_t igs=0; igs<kNFitConfCases; igs++){
    if (igs==kFixSigUpFreeMean && !fUseFixSigUpFreeMean) continue;
    if (igs==kFixSigDownFreeMean && !fUseFixSigDownFreeMean) continue;
    if (igs==kFreeSigFixMean  && !fUseFixedMeanFreeS) continue;
    if (igs==kFreeSigFreeMean  && !fUseFreeS) continue;
    if (igs==kFixSigFreeMean  && !fUseFixSigFreeMean) continue;
    if (igs==kFixSigFixMean   && !fUseFixSigFixMean) continue;
    cout << igs << "  (" << kNFitConfCases <<","<< kFixSigDownFreeMean <<","<< kFreeSigFixMean <<","<<kFreeSigFreeMean <<","<<kFixSigFreeMean <<","<< kFixSigFixMean<< endl;
    kNFitConfCasesTurnedOn++;
  }
  
  Int_t nCasesXaxisnonCst=kNBkgFuncCasesTurnedOn*kNFitConfCasesTurnedOn;
  const Int_t nCases=kNBkgFuncCases*kNFitConfCases;
  
  Double_t maxAllhistograms = 100000;
  Int_t nbinsAllhistograms = 100000;
  if(fUseNoBackgroundOnlySignal){ maxAllhistograms = 100000; nbinsAllhistograms = 50000; nCasesXaxisnonCst=1; }
  const Int_t nCasesXaxis = (const Int_t) nCasesXaxisnonCst;
  
  TString funcBkg[kNBkgFuncCases]={"PowLaw","PowLawExpo""Expo","Lin","Pol2","Pol3","Pol4","Pol5","noBkg"};
  TString gausSig[kNFitConfCases]={"FixedS","FixedSp20","FixedSm20","FreeS","FixedMeanFixedS","FixedMeanFreeS"};
  
  Int_t totTrials=fNumOfRebinSteps*fNumOfFirstBinSteps*fNumOfLowLimFitSteps*fNumOfUpLimFitSteps;
  fHistoRawYieldDistAll = new TH1F(Form("hRawYieldDistAll%s",fSuffix.Data()),"  ; Raw Yield",nbinsAllhistograms,0.,maxAllhistograms);
  fHistoRawYieldTrialAll = new TH1F(Form("hRawYieldTrialAll%s",fSuffix.Data())," ; Trial # ; Raw Yield",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  fHistoSigmaTrialAll = new TH1F(Form("hSigmaTrialAll%s",fSuffix.Data())," ; Trial # ; Sigma (GeV/c^{2})",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  fHistoMeanTrialAll = new TH1F(Form("hMeanTrialAll%s",fSuffix.Data())," ; Trial # ; Mean (GeV/c^{2})",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  fHistoChi2TrialAll = new TH1F(Form("hChi2TrialAll%s",fSuffix.Data()),"  ; Trial # ; #chi^{2}",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  fHistoSignifTrialAll = new TH1F(Form("hSignifTrialAll%s",fSuffix.Data()),"  ; Trial # ; Significance",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  if(fSaveBkgVal) {
    fHistoBkgTrialAll = new TH1F(Form("hBkgTrialAll%s",fSuffix.Data()),"  ; Background",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
    fHistoBkgInBinEdgesTrialAll = new TH1F(Form("hBkgInBinEdgesTrialAll%s",fSuffix.Data()),"  ; Background in bin edges",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  }
  
  
  fHistoRawYieldDistBinC0All = new TH1F(Form("hRawYieldDistBinC0All%s",fSuffix.Data()),"  ; Raw Yield (bin count)",nbinsAllhistograms,0.,maxAllhistograms);
  fHistoRawYieldTrialBinC0All = new TH1F(Form("hRawYieldTrialBinC0All%s",fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  fHistoRawYieldDistBinC1All = new TH1F(Form("hRawYieldDistBinC1All%s",fSuffix.Data()),"  ; Raw Yield (bin count)",nbinsAllhistograms,0.,maxAllhistograms);
  fHistoRawYieldTrialBinC1All = new TH1F(Form("hRawYieldTrialBinC1All%s",fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  
  fHistoRawYieldDistBinC0All_2 = new TH1F(Form("hRawYieldDistBinC0All%s_2",fSuffix.Data()),"  ; Raw Yield (bin count)",nbinsAllhistograms,0.,maxAllhistograms);
  fHistoRawYieldTrialBinC0All_2 = new TH1F(Form("hRawYieldTrialBinC0All%s_2",fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  fHistoRawYieldDistBinC1All_2 = new TH1F(Form("hRawYieldDistBinC1All%s_2",fSuffix.Data()),"  ; Raw Yield (bin count)",nbinsAllhistograms,0.,maxAllhistograms);
  fHistoRawYieldTrialBinC1All_2 = new TH1F(Form("hRawYieldTrialBinC1All%s_2",fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",nCasesXaxis*totTrials,-0.5,nCasesXaxis*totTrials-0.5);
  
  fHistoRawYieldDist = new TH1F*[nCases];
  fHistoRawYieldTrial = new TH1F*[nCases];
  fHistoSigmaTrial = new TH1F*[nCases];
  fHistoMeanTrial = new TH1F*[nCases];
  fHistoChi2Trial = new TH1F*[nCases];
  fHistoSignifTrial = new TH1F*[nCases];
  if(fSaveBkgVal) {
    fHistoBkgTrial = new TH1F*[nCases];
    fHistoBkgInBinEdgesTrial = new TH1F*[nCases];
  }
  /*
   fHistoRawYieldDistBinC0 = new TH1F*[nCases];
   fHistoRawYieldTrialBinC0 = new TH2F*[nCases];
   fHistoRawYieldDistBinC1 = new TH1F*[nCases];
   fHistoRawYieldTrialBinC1 = new TH2F*[nCases];
   
   fHistoRawYieldDistBinC0_2 = new TH1F*[nCases];
   fHistoRawYieldTrialBinC0_2 = new TH2F*[nCases];
   fHistoRawYieldDistBinC1_2 = new TH1F*[nCases];
   fHistoRawYieldTrialBinC1_2 = new TH2F*[nCases];
   
   for(Int_t ib=0; ib<kNBkgFuncCases; ib++){
   for(Int_t igs=0; igs<kNFitConfCases; igs++){
   Int_t theCase;
   if(ib == 6 && igs == 0) theCase = 1*2+0;
   else if(ib == 6 && igs == 3) theCase = 0*2+0;
   else if(ib == 7 && igs == 0) theCase = 1*2+1;
   else if(ib == 7 && igs == 3) theCase = 0*2+1;
   else{ theCase = -99; }
   //Int_t theCase=igs*kNBkgFuncCases+ib;
   fHistoRawYieldDist[theCase]=new TH1F(Form("hRawYieldDist%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Raw Yield",2500,0.,5000.);
   fHistoRawYieldDistBinC0[theCase]=new TH1F(Form("hRawYieldDistBinC0%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Raw Yield (bin count)",2500,0.,5000.);
   fHistoRawYieldDistBinC1[theCase]=new TH1F(Form("hRawYieldDistBinC1%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Raw Yield (bin count)",2500,0.,5000.);
   fHistoRawYieldDistBinC0_2[theCase]=new TH1F(Form("hRawYieldDistBinC0%s%s%s_2",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Raw Yield (bin count)",2500,0.,5000.);
   fHistoRawYieldDistBinC1_2[theCase]=new TH1F(Form("hRawYieldDistBinC1%s%s%s_2",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Raw Yield (bin count)",2500,0.,5000.);
   fHistoRawYieldTrial[theCase]=new TH1F(Form("hRawYieldTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Raw Yield",totTrials,-0.5,totTrials-0.5);
   fHistoRawYieldTrialBinC0[theCase]=new TH2F(Form("hRawYieldTrialBinC0%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);
   fHistoRawYieldTrialBinC1[theCase]=new TH2F(Form("hRawYieldTrialBinC1%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);
   fHistoRawYieldTrialBinC0_2[theCase]=new TH2F(Form("hRawYieldTrialBinC0%s%s%s_2",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);
   fHistoRawYieldTrialBinC1_2[theCase]=new TH2F(Form("hRawYieldTrialBinC1%s%s%s_2",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Range for count ; Raw Yield (bin count)",totTrials,-0.5,totTrials-0.5,fNumOfnSigmaBinCSteps,-0.5,fNumOfnSigmaBinCSteps-0.5);
   fHistoSigmaTrial[theCase]=new TH1F(Form("hSigmaTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Sigma (GeV/c^{2})",totTrials,-0.5,totTrials-0.5);
   fHistoMeanTrial[theCase]=new TH1F(Form("hMeanTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Mean (GeV/c^{2})",totTrials,-0.5,totTrials-0.5);
   fHistoChi2Trial[theCase]=new TH1F(Form("hChi2Trial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; #chi^{2}",totTrials,-0.5,totTrials-0.5);
   fHistoSignifTrial[theCase]=new TH1F(Form("hSignifTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data())," ; Trial # ; Significance",totTrials,-0.5,totTrials-0.5);
   if(fSaveBkgVal) {
   fHistoBkgTrial[theCase] = new TH1F(Form("hBkgTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Background",totTrials,-0.5,totTrials-0.5);
   fHistoBkgInBinEdgesTrial[theCase] = new TH1F(Form("hBkgInBinEdgesTrial%s%s%s",funcBkg[ib].Data(),gausSig[igs].Data(),fSuffix.Data()),"  ; Background in bin edges",totTrials,-0.5,totTrials-0.5);
   }
   
   fHistoChi2Trial[theCase]->SetMarkerStyle(7);
   fHistoSignifTrial[theCase]->SetMarkerStyle(7);
   fHistoSigmaTrial[theCase]->SetMarkerStyle(7);
   fHistoMeanTrial[theCase]->SetMarkerStyle(7);
   if(fSaveBkgVal) {
   fHistoBkgTrial[theCase]->SetMarkerStyle(7);
   fHistoBkgInBinEdgesTrial[theCase]->SetMarkerStyle(7);
   }
   
   }
   }
   */
  fHistoRawYieldTrialBinC1All->SetMarkerColor(kCyan);
  fHistoRawYieldTrialBinC1All_2->SetMarkerColor(kBlack);
  fHistoRawYieldTrialBinC1All->SetLineColor(kGreen+2);
  fHistoRawYieldTrialBinC1All_2->SetLineColor(kOrange+1);
  fHistoRawYieldDistBinC1All->SetMarkerColor(kCyan);
  fHistoRawYieldDistBinC1All_2->SetMarkerColor(kBlack);
  fHistoRawYieldDistBinC1All->SetLineColor(kGreen+2);
  fHistoRawYieldDistBinC1All_2->SetLineColor(kOrange+1);
  
  fHistoRawYieldTrialBinC1All_2->GetYaxis()->SetTitle("raw yield");
  fHistoRawYieldTrialBinC1All_2->GetXaxis()->SetTitle("Trial #");
  fHistoRawYieldTrialBinC1All_2->SetTitle(Form("%0.1f < #it{p}_{T} < %0.1f (GeV/#it{c})",ptmin,ptmax));
  fHistoRawYieldTrialBinC1All_2->SetMarkerStyle(20);
  fHistoRawYieldTrialBinC1All_2->SetLineWidth(1);
  fHistoRawYieldTrialBinC1All_2->SetMarkerSize(0.5);
  fHistoRawYieldTrialBinC1All_2->SetLineColor(kOrange+1);
  fHistoRawYieldTrialBinC1All_2->SetMarkerColor(kOrange+1);
  
  fHistoRawYieldTrialBinC1All->GetYaxis()->SetTitle("raw yield");
  fHistoRawYieldTrialBinC1All->GetXaxis()->SetTitle("Trial #");
  fHistoRawYieldTrialBinC1All->SetTitle(Form("%0.1f < #it{p}_{T} < %0.1f (GeV/#it{c})",ptmin,ptmax));
  fHistoRawYieldTrialBinC1All->SetMarkerStyle(20);
  fHistoRawYieldTrialBinC1All->SetLineWidth(1);
  fHistoRawYieldTrialBinC1All->SetMarkerSize(0.5);
  fHistoRawYieldTrialBinC1All->SetLineColor(kGreen+2);
  fHistoRawYieldTrialBinC1All->SetMarkerColor(kGreen+2);
  
  fHistoRawYieldTrialAll->GetYaxis()->SetTitle("raw yield");
  fHistoRawYieldTrialAll->GetXaxis()->SetTitle("Trial #");
  fHistoRawYieldTrialAll->SetTitle(Form("%0.1f < #it{p}_{T} < %0.1f (GeV/#it{c})",ptmin,ptmax));
  fHistoRawYieldTrialAll->SetMarkerStyle(20);
  fHistoRawYieldTrialAll->SetLineWidth(1);
  fHistoRawYieldTrialAll->SetMarkerSize(0.5);
  fHistoRawYieldTrialAll->SetLineColor(kBlue+1);
  fHistoRawYieldTrialAll->SetMarkerColor(kBlack);
  
  fHistoSigmaTrialAll->GetYaxis()->SetRangeUser(0.001,0.05);
  fHistoSigmaTrialAll->GetYaxis()->SetTitle("width (GeV/c^{2})");
  fHistoSigmaTrialAll->GetXaxis()->SetTitle("Trial #");
  fHistoSigmaTrialAll->SetTitle(Form("%0.1f < #it{p}_{T} < %0.1f (GeV/c)",ptmin,ptmax));
  fHistoSigmaTrialAll->SetMarkerStyle(20);
  fHistoSigmaTrialAll->SetLineWidth(1);
  fHistoSigmaTrialAll->SetMarkerSize(0.5);
  fHistoSigmaTrialAll->SetLineColor(kBlue+1);
  fHistoSigmaTrialAll->SetMarkerColor(kBlack);
  
  fHistoMeanTrialAll->GetYaxis()->SetRangeUser(5.35, 5.42);
  fHistoMeanTrialAll->GetYaxis()->SetTitle("mean (GeV/c^{2})");
  fHistoMeanTrialAll->GetXaxis()->SetTitle("Trial #");
  fHistoMeanTrialAll->SetTitle(Form("%0.1f < #it{p}_{T} < %0.1f (GeV/c)",ptmin,ptmax));
  fHistoMeanTrialAll->SetMarkerStyle(20);
  fHistoMeanTrialAll->SetLineWidth(1);
  fHistoMeanTrialAll->SetMarkerSize(0.5);
  fHistoMeanTrialAll->SetLineColor(kBlue+1);
  fHistoMeanTrialAll->SetMarkerColor(kBlack);
  
  fHistoChi2TrialAll->Draw("p");
  fHistoChi2TrialAll->GetYaxis()->SetRangeUser(0,maxchisquare);
  fHistoChi2TrialAll->GetYaxis()->SetTitle("#chi^{2}");
  fHistoChi2TrialAll->GetXaxis()->SetTitle("Trial #");
  fHistoChi2TrialAll->SetTitle(Form("%0.1f < #it{p}_{T} < %0.1f (GeV/#it{c})",ptmin,ptmax));
  fHistoChi2TrialAll->SetMarkerStyle(20);
  fHistoChi2TrialAll->SetLineWidth(1);
  fHistoChi2TrialAll->SetMarkerSize(0.5);
  fHistoChi2TrialAll->SetLineColor(kBlue+1);
  fHistoChi2TrialAll->SetMarkerColor(kBlack);
  
  fHistoRawYieldDistBinC1All->GetXaxis()->SetTitle("raw yield");
  fHistoRawYieldDistBinC1All->GetYaxis()->SetTitle("Entries");
  fHistoRawYieldDistBinC1All->SetTitle(Form("%0.1f < #it{p}_{T} < %0.1f (GeV/c)",ptmin,ptmax));
  fHistoRawYieldDistBinC1All->SetFillStyle(3004);
  fHistoRawYieldDistBinC1All->SetLineWidth(2);
  fHistoRawYieldDistBinC1All->SetLineColor(kGreen+2);
  fHistoRawYieldDistBinC1All->SetFillColor(kGreen+2);
  
  fHistoRawYieldDistBinC1All_2->GetXaxis()->SetTitle("raw yield");
  fHistoRawYieldDistBinC1All_2->GetYaxis()->SetTitle("Entries");
  fHistoRawYieldDistBinC1All_2->SetTitle(Form("%0.1f < #it{p}_{T} < %0.1f (GeV/c)",ptmin,ptmax));
  fHistoRawYieldDistBinC1All_2->SetFillStyle(3004);
  fHistoRawYieldDistBinC1All_2->SetLineWidth(2);
  fHistoRawYieldDistBinC1All_2->SetLineColor(kOrange+1);
  fHistoRawYieldDistBinC1All_2->SetFillColor(kOrange+1);
  
  fHistoRawYieldDistAll->GetXaxis()->SetTitle("raw yield");
  fHistoRawYieldDistAll->GetYaxis()->SetTitle("Entries");
  fHistoRawYieldDistAll->SetTitle(Form("%0.1f < #it{p}_{T} < %0.1f (GeV/c)",ptmin,ptmax));
  fHistoRawYieldDistAll->SetFillStyle(3004);
  fHistoRawYieldDistAll->SetLineWidth(2);
  fHistoRawYieldDistAll->SetLineColor(kBlue+1);
  fHistoRawYieldDistAll->SetFillColor(kBlue+1);
  
  fNtupleMultiTrials = new TNtuple(Form("ntuMultiTrial%s",fSuffix.Data()),Form("ntuMultiTrial%s",fSuffix.Data()),"rebin:firstb:minfit:maxfit:bkgfunc:confsig:confmean:chi2:signif:mean:emean:sigma:esigma:rawy:erawy",128000);
  fNtupleMultiTrials->SetDirectory(nullptr);
  return kTRUE;
  
}

//________________________________________________________________________
Bool_t AliHFInvMassMultiTrialFit1::DoMultiTrials(TH1D* hInvMassHisto, TPad* thePad){
  // perform the multiple fits
  
  Bool_t hOK=CreateHistos();
  if(!hOK) return kFALSE;
  
  Int_t itrial=0;
  Int_t types=0;
  if(fUseDoubleGausSignal) types = 1;
  Int_t itrialBC=0;
  Int_t totTrials=fNumOfRebinSteps*fNumOfFirstBinSteps*fNumOfLowLimFitSteps*fNumOfUpLimFitSteps;
  
  fMinYieldGlob=999999.;
  fMaxYieldGlob=0.;
  fMinYieldGlobBC=999999.;
  fMaxYieldGlobBC=0.;
  fMinYieldGlobBC_2=999999.;
  fMaxYieldGlobBC_2=0.;
  Float_t xnt[15];
  
  Int_t kNBkgFuncCasesTurnedOn = 0;
  for(Int_t typeb=0; typeb<kNBkgFuncCases; typeb++){
    if(typeb==kExpoBkg && !fUseExpoBkg) continue;
    if(typeb==kLinBkg && !fUseLinBkg) continue;
    if(typeb==kPol2Bkg && !fUsePol2Bkg) continue;
    if(typeb==kPol3Bkg && !fUsePol3Bkg) continue;
    if(typeb==kPol4Bkg && !fUsePol4Bkg) continue;
    if(typeb==kPol5Bkg && !fUsePol5Bkg) continue;
    if(typeb==kPowBkg && !fUsePowLawBkg) continue;
    if(typeb==kPowTimesExpoBkg && !fUsePowLawTimesExpoBkg) continue;
    if(typeb==kNoBkgOnlySig && !fUseNoBackgroundOnlySignal) continue;
    kNBkgFuncCasesTurnedOn++;
  }
  
  for(Int_t ir=0; ir<fNumOfRebinSteps; ir++){
    Int_t rebin=fRebinSteps[ir];
    for(Int_t iFirstBin=1; iFirstBin<=fNumOfFirstBinSteps; iFirstBin++) {
      TH1F* hRebinned=0x0;
      if(fNumOfFirstBinSteps==1) hRebinned=(TH1F*)AliVertexingHFUtils::RebinHisto(hInvMassHisto,rebin,-1);
      else hRebinned=(TH1F*)AliVertexingHFUtils::RebinHisto(hInvMassHisto,rebin,iFirstBin);
      for(Int_t iMinMass=0; iMinMass<fNumOfLowLimFitSteps; iMinMass++){
        Double_t minMassForFit=fLowLimFitSteps[iMinMass];
        Double_t hmin=TMath::Max(minMassForFit,hRebinned->GetBinLowEdge(2));
        for(Int_t iMaxMass=0; iMaxMass<fNumOfUpLimFitSteps; iMaxMass++){
          Double_t maxMassForFit=fUpLimFitSteps[iMaxMass];
          Double_t hmax=TMath::Min(maxMassForFit,hRebinned->GetBinLowEdge(hRebinned->GetNbinsX()));
          ++itrial;
          for(Int_t typeb=0; typeb<kNBkgFuncCases; typeb++){
            if(typeb==kExpoBkg && !fUseExpoBkg) continue;
            if(typeb==kLinBkg && !fUseLinBkg) continue;
            if(typeb==kPol2Bkg && !fUsePol2Bkg) continue;
            if(typeb==kPol3Bkg && !fUsePol3Bkg) continue;
            if(typeb==kPol4Bkg && !fUsePol4Bkg) continue;
            if(typeb==kPol5Bkg && !fUsePol5Bkg) continue;
            if(typeb==kPowBkg && !fUsePowLawBkg) continue;
            if(typeb==kPowTimesExpoBkg && !fUsePowLawTimesExpoBkg) continue;
            if(typeb==kNoBkgOnlySig && !fUseNoBackgroundOnlySignal) continue;
            for(Int_t igs=0; igs<kNFitConfCases; igs++){
              if (igs==kFixSigUpFreeMean && !fUseFixSigUpFreeMean) continue;
              if (igs==kFixSigDownFreeMean && !fUseFixSigDownFreeMean) continue;
              if (igs==kFreeSigFixMean  && !fUseFixedMeanFreeS) continue;
              if (igs==kFreeSigFreeMean  && !fUseFreeS) continue;
              if (igs==kFixSigFreeMean  && !fUseFixSigFreeMean) continue;
              if (igs==kFixSigFixMean   && !fUseFixSigFixMean) continue;
              
              Int_t theCase;
              if(kNBkgFuncCasesTurnedOn == 3){
                if(typeb == 2 && igs == 0) theCase = 1*3+0;
                else if(typeb == 2 && igs == 3) theCase = 2*3+0;
                else if(typeb == 2 && igs == 5) theCase = 3*3+0;
                else if(typeb == 2 && igs == 4) theCase = 0*3+0;
                else if(typeb == 0 && igs == 0) theCase = 1*3+1;
                else if(typeb == 0 && igs == 3) theCase = 2*3+1;
                else if(typeb == 0 && igs == 5) theCase = 3*3+1;
                else if(typeb == 0 && igs == 4) theCase = 0*3+1;
                else if(typeb == 1 && igs == 0) theCase = 1*3+2;
                else if(typeb == 1 && igs == 3) theCase = 2*3+2;
                else if(typeb == 1 && igs == 5) theCase = 3*3+2;
                else if(typeb == 1 && igs == 4) theCase = 0*3+2;
                else if(typeb == 6 && igs == 0) theCase = 1*3+2;
                else if(typeb == 6 && igs == 3) theCase = 2*3+2;
                else if(typeb == 6 && igs == 5) theCase = 3*3+2;
                else if(typeb == 6 && igs == 4) theCase = 0*3+2;
                else{ theCase = -99; cout << "Different: " << typeb << " "  << igs << endl; }
              } else {
                if(typeb == 0 && igs == 0) theCase = 1*2+0;
                else if(typeb == 0 && igs == 3) theCase = 2*2+0;
                else if(typeb == 0 && igs == 5) theCase = 3*2+0;
                else if(typeb == 0 && igs == 4) theCase = 0*2+0;
                else if(typeb == 2 && igs == 0) theCase = 1*2+1;
                else if(typeb == 2 && igs == 3) theCase = 2*2+1;
                else if(typeb == 2 && igs == 5) theCase = 3*2+1;
                else if(typeb == 2 && igs == 4) theCase = 0*2+1;
                else if(typeb == 1 && igs == 0) theCase = 1*2+1;
                else if(typeb == 1 && igs == 3) theCase = 2*2+1;
                else if(typeb == 1 && igs == 5) theCase = 3*2+1;
                else if(typeb == 1 && igs == 5) theCase = 0*2+1;
                else if(typeb == 6 && igs == 0) theCase = 1*2+1;
                else if(typeb == 6 && igs == 3) theCase = 2*2+1;
                else if(typeb == 6 && igs == 5) theCase = 3*2+1;
                else if(typeb == 6 && igs == 4) theCase = 0*2+1;
                else{ theCase = -99; cout << "Different: " << typeb << " "  << igs << endl; }
              }
            
              //Int_t theCase=igs*kNBkgFuncCases+typeb;
              Int_t globBin=itrial+theCase*totTrials;
              for(Int_t j=0; j<15; j++) xnt[j]=0.;
              
              Bool_t mustDeleteFitter = kTRUE;
              AliHFInvMassFitter1*  fitter=0x0;
              if(typeb==kExpoBkg){
                fitter=new AliHFInvMassFitter1(hRebinned, hmin, hmax, AliHFInvMassFitter1::kExpo, types);
              }else if(typeb==kLinBkg){
                fitter=new AliHFInvMassFitter1(hRebinned, hmin, hmax, AliHFInvMassFitter1::kLin, types);
              }else if(typeb==kPol2Bkg){
                fitter=new AliHFInvMassFitter1(hRebinned, hmin, hmax, AliHFInvMassFitter1::kPol2, types);
              }else if(typeb==kPowBkg){
                fitter=new AliHFInvMassFitter1(hRebinned, hmin, hmax, AliHFInvMassFitter1::kPow, types);
              }else if(typeb==kPowTimesExpoBkg){
                fitter=new AliHFInvMassFitter1(hRebinned, hmin, hmax, AliHFInvMassFitter1::kPowEx, types);
              }else if(typeb==kNoBkgOnlySig){
                fitter=new AliHFInvMassFitter1(hRebinned, hmin, hmax, AliHFInvMassFitter1::kNoBk, types);
                //if(fRangeGaus2sigma) fitter->SetRangeGaus2sigma(kTRUE);
              }else{
                fitter=new AliHFInvMassFitter1(hRebinned, hmin, hmax, 6, types);
                if(typeb==kPol3Bkg) fitter->SetPolDegreeForBackgroundFit(3);
                if(typeb==kPol4Bkg) fitter->SetPolDegreeForBackgroundFit(4);
                if(typeb==kPol5Bkg) fitter->SetPolDegreeForBackgroundFit(5);
              }
              // D0 Reflection
              if(fhTemplRefl && fhTemplSign){
                TH1F *hReflModif=(TH1F*)AliVertexingHFUtils::AdaptTemplateRangeAndBinning(fhTemplRefl,hRebinned,minMassForFit,maxMassForFit);
                TH1F *hSigModif=(TH1F*)AliVertexingHFUtils::AdaptTemplateRangeAndBinning(fhTemplSign,hRebinned,minMassForFit,maxMassForFit);
                TH1F* hrfl=fitter->SetTemplateReflections(hReflModif,"2gaus",minMassForFit,maxMassForFit);
                if(fFixRefloS>0){
                  Double_t fixSoverRefAt=fFixRefloS*(hReflModif->Integral(hReflModif->FindBin(minMassForFit*1.0001),hReflModif->FindBin(maxMassForFit*0.999))/hSigModif->Integral(hSigModif->FindBin(minMassForFit*1.0001),hSigModif->FindBin(maxMassForFit*0.999)));
                  fitter->SetFixReflOverS(fixSoverRefAt);
                }
                delete hReflModif;
                delete hSigModif;
              }
              if(fUseSecondPeak){
                fitter->IncludeSecondGausPeak(fMassSecondPeak, fFixMassSecondPeak, fSigmaSecondPeak, fFixSigmaSecondPeak);
              }
              if(fFitOption==1) fitter->SetUseChi2Fit();
              fitter->SetInitialGaussianMean(fMassD);
              fitter->SetInitialGaussianSigma(fSigmaGausMC);
              xnt[0]=rebin;
              xnt[1]=iFirstBin;
              xnt[2]=minMassForFit;
              xnt[3]=maxMassForFit;
              xnt[4]=typeb;
              xnt[6]=0;
              if(igs==kFixSigFreeMean){
                fitter->SetFixGaussianSigma(fSigmaGausMC);
                xnt[5]=1;
              }else if(igs==kFixSigUpFreeMean){
                fitter->SetFixGaussianSigma(fSigmaGausMC*(1.+fSigmaMCVariation));
                xnt[5]=2;
              }else if(igs==kFixSigDownFreeMean){
                fitter->SetFixGaussianSigma(fSigmaGausMC*(1.-fSigmaMCVariation));
                xnt[5]=3;
              }else if(igs==kFreeSigFreeMean){
                xnt[5]=0;
              }else if(igs==kFixSigFixMean){
                fitter->SetFixGaussianSigma(fSigmaGausMC);
                fitter->SetFixGaussianMean(fMassD);
                xnt[5]=1;
                xnt[6]=1;
              }else if(igs==kFreeSigFixMean){
                fitter->SetFixGaussianMean(fMassD);
                xnt[5]=0;
                xnt[6]=1;
              }
              Bool_t out=kFALSE;
              Double_t chisq=-1.;
              Double_t sigma=0.;
              Double_t esigma=0.;
              Double_t pos=.0;
              Double_t epos=.0;
              Double_t ry=.0;
              Double_t ery=.0;
              Double_t significance=0.;
              Double_t erSignif=0.;
              Double_t bkg=0.;
              Double_t erbkg=0.;
              Double_t bkgBEdge=0;
              Double_t erbkgBEdge=0;
              TF1* fB1=0x0;
              if(typeb<kNBkgFuncCases){
                printf("****** START FIT OF HISTO %s WITH REBIN %d FIRST BIN %d MASS RANGE %f-%f BACKGROUND FIT FUNCTION=%d CONFIG SIGMA/MEAN=%d GLOBIN=%d\n",hInvMassHisto->GetName(),rebin,iFirstBin,minMassForFit,maxMassForFit,typeb,igs,globBin);
                out=fitter->MassFitter(0);
                chisq=fitter->GetReducedChiSquare();
                fitter->Significance(fnSigmaForBkgEval,significance,erSignif);
                sigma=fitter->GetSigma();
                pos=fitter->GetMean();
                esigma=fitter->GetSigmaUncertainty();
                if(esigma<0.00001) esigma=0.000001;
                epos=fitter->GetMeanUncertainty();
                if(epos<0.00001) epos=0.000001;
                ry=fitter->GetRawYield();
                ery=fitter->GetRawYieldError();
                fB1=fitter->GetBackgroundFullRangeFunc();
                fitter->Background(fnSigmaForBkgEval,bkg,erbkg);
                Double_t minval = hInvMassHisto->GetXaxis()->GetBinLowEdge(hInvMassHisto->FindBin(pos-fnSigmaForBkgEval*sigma));
                Double_t maxval = hInvMassHisto->GetXaxis()->GetBinUpEdge(hInvMassHisto->FindBin(pos+fnSigmaForBkgEval*sigma));
                fitter->Background(minval,maxval,bkgBEdge,erbkgBEdge);
                if(out && fDrawIndividualFits && thePad){
                  thePad->Clear();
                  fitter->DrawHere(thePad, fnSigmaForBkgEval);
                  fMassFitters.push_back(fitter);
                  mustDeleteFitter = kFALSE;
                  //for (auto format : fInvMassFitSaveAsFormats) {
                  //thePad->SaveAs(Form("FitOutput_%s_Trial%d.%s",hInvMassHisto->GetName(),globBin, format.c_str()));
                  thePad->SaveAs(Form("FitOutput_%s_Trial%d.png",hInvMassHisto->GetName(),globBin));
                  //}
                }
              }
              // else{
              //   out=DoFitWithPol3Bkg(hRebinned,hmin,hmax,igs);
              //   if(out && thePad){
              // 	thePad->Clear();
              // 	hRebinned->Draw();
              // 	TF1* fSB=(TF1*)hRebinned->GetListOfFunctions()->FindObject("fSB");
              // 	fB1=new TF1("fB1","[0]+[1]*x+[2]*x*x+[3]*x*x*x",hmin,hmax);
              // 	for(Int_t j=0; j<4; j++) fB1->SetParameter(j,fSB->GetParameter(3+j));
              // 	fB1->SetLineColor(2);
              // 	fB1->Draw("same");
              // 	fSB->SetLineColor(4);
              // 	fSB->Draw("same");
              // 	thePad->Update();
              // 	chisq=fSB->GetChisquare()/fSB->GetNDF();;
              // 	sigma=fSB->GetParameter(2);
              // 	esigma=fSB->GetParError(2);
              // 	if(esigma<0.00001) esigma=0.0001;
              // 	pos=fSB->GetParameter(1);
              // 	epos=fSB->GetParError(1);
              // 	if(epos<0.00001) epos=0.0001;
              // 	ry=fSB->GetParameter(0)/hRebinned->GetBinWidth(1);
              // 	ery=fSB->GetParError(0)/hRebinned->GetBinWidth(1);
              //   }
              // }
              xnt[7]=chisq;
              if(out && chisq>0. && chisq < maxchisquare && sigma>0.7*fSigmaGausMC && sigma<1.3*fSigmaGausMC){
                xnt[8]=significance;
                xnt[9]=pos;
                xnt[10]=epos;
                xnt[11]=sigma;
                xnt[12]=esigma;
                xnt[13]=ry;
                xnt[14]=ery;
                fHistoRawYieldDistAll->Fill(ry);
                fHistoRawYieldTrialAll->SetBinContent(globBin,ry);
                fHistoRawYieldTrialAll->SetBinError(globBin,ery);
                fHistoSigmaTrialAll->SetBinContent(globBin,sigma);
                fHistoSigmaTrialAll->SetBinError(globBin,esigma);
                fHistoMeanTrialAll->SetBinContent(globBin,pos);
                fHistoMeanTrialAll->SetBinError(globBin,epos);
                fHistoChi2TrialAll->SetBinContent(globBin,chisq);
                fHistoChi2TrialAll->SetBinError(globBin,0.00001);
                fHistoSignifTrialAll->SetBinContent(globBin,significance);
                fHistoSignifTrialAll->SetBinError(globBin,erSignif);
                if(fSaveBkgVal) {
                  fHistoBkgTrialAll->SetBinContent(globBin,bkg);
                  fHistoBkgTrialAll->SetBinError(globBin,erbkg);
                  fHistoBkgInBinEdgesTrialAll->SetBinContent(globBin,bkgBEdge);
                  fHistoBkgInBinEdgesTrialAll->SetBinError(globBin,erbkgBEdge);
                }
                
                if(ry<fMinYieldGlob) fMinYieldGlob=ry;
                if(ry>fMaxYieldGlob) fMaxYieldGlob=ry;
                /*
                 fHistoRawYieldDist[theCase]->Fill(ry);
                 fHistoRawYieldTrial[theCase]->SetBinContent(itrial,ry);
                 fHistoRawYieldTrial[theCase]->SetBinError(itrial,ery);
                 fHistoSigmaTrial[theCase]->SetBinContent(itrial,sigma);
                 fHistoSigmaTrial[theCase]->SetBinError(itrial,esigma);
                 fHistoMeanTrial[theCase]->SetBinContent(itrial,pos);
                 fHistoMeanTrial[theCase]->SetBinError(itrial,epos);
                 fHistoChi2Trial[theCase]->SetBinContent(itrial,chisq);
                 fHistoChi2Trial[theCase]->SetBinError(itrial,0.00001);
                 fHistoSignifTrial[theCase]->SetBinContent(itrial,significance);
                 fHistoSignifTrial[theCase]->SetBinError(itrial,erSignif);
                 if(fSaveBkgVal) {
                 fHistoBkgTrial[theCase]->SetBinContent(itrial,bkg);
                 fHistoBkgTrial[theCase]->SetBinError(itrial,erbkg);
                 fHistoBkgInBinEdgesTrial[theCase]->SetBinContent(itrial,bkgBEdge);
                 fHistoBkgInBinEdgesTrial[theCase]->SetBinError(itrial,erbkgBEdge);
                 }
                 */
                TF1* fBkg=fitter->GetBackgroundRecalcFunc();
                for(Int_t iStepBC=0; iStepBC<fNumOfnSigmaBinCSteps; iStepBC++){
                  Double_t minMassBC=pos-fnSigmaBinCSteps[iStepBC]*sigma;
                  Double_t maxMassBC=pos+fnSigmaBinCSteps[iStepBC]*sigma;
                  if(minMassBC>minMassForFit &&
                     maxMassBC<maxMassForFit &&
                     minMassBC>(hRebinned->GetXaxis()->GetXmin()) &&
                     maxMassBC<(hRebinned->GetXaxis()->GetXmax())){
                    Double_t cnts0 = 0;
                    Double_t ecnts0 = 0;
                    Double_t cnts1 = 0;
                    Double_t ecnts1 = 0;
                    cnts0=fitter->GetRawYieldBinCounting(ecnts0,minMassBC,maxMassBC,0);
                    cnts1=fitter->GetRawYieldBinCounting(ecnts1,minMassBC,maxMassBC,1);
                    //          for(Int_t iMB=hRebinned->FindBin(minMassBC*1.001); iMB<=hRebinned->FindBin(maxMassBC*0.999); iMB++){
                    //              Double_t bkg=fBkg ? fBkg->Eval(hRebinned->GetBinCenter(iMB)) : 0;
                    //              cnts0 += (hRebinned->GetBinContent(iMB)-bkg); //Dirty fix, not correct
                    //              ecnts0 += (hRebinned->GetBinContent(iMB));  //Dirty fix, not correct
                    //              cnts1 += (hRebinned->GetBinContent(iMB)-bkg);
                    //              ecnts1 += (hRebinned->GetBinContent(iMB));
                    //         }
                    //          ecnts0=TMath::Sqrt(ecnts0);
                    //          ecnts1=TMath::Sqrt(ecnts1);
                    
                    ++itrialBC;
                    if(iStepBC == 0){
                      fHistoRawYieldDistBinC0All->Fill(cnts0);
                      fHistoRawYieldTrialBinC0All->SetBinContent(globBin,cnts0);
                      fHistoRawYieldTrialBinC0All->SetBinError(globBin,ecnts0);
                      //fHistoRawYieldTrialBinC0[theCase]->SetBinContent(itrial,iStepBC+1,cnts0);
                      //fHistoRawYieldTrialBinC0[theCase]->SetBinError(itrial,iStepBC+1,ecnts0);
                      //fHistoRawYieldDistBinC0[theCase]->Fill(cnts0);
                      if(cnts1<fMinYieldGlobBC) fMinYieldGlobBC=cnts1;
                      if(cnts1>fMaxYieldGlobBC) fMaxYieldGlobBC=cnts1;
                      fHistoRawYieldDistBinC1All->Fill(cnts1);
                      fHistoRawYieldTrialBinC1All->SetBinContent(globBin,cnts1);
                      fHistoRawYieldTrialBinC1All->SetBinError(globBin,ecnts1);
                      //fHistoRawYieldTrialBinC1[theCase]->SetBinContent(itrial,iStepBC+1,cnts1);
                      //fHistoRawYieldTrialBinC1[theCase]->SetBinError(itrial,iStepBC+1,ecnts1);
                      //fHistoRawYieldDistBinC1[theCase]->Fill(cnts1);
                    } else if(iStepBC == 1) {
                      fHistoRawYieldDistBinC0All_2->Fill(cnts0);
                      fHistoRawYieldTrialBinC0All_2->SetBinContent(globBin,cnts0);
                      fHistoRawYieldTrialBinC0All_2->SetBinError(globBin,ecnts0);
                      //fHistoRawYieldTrialBinC0_2[theCase]->SetBinContent(itrial,iStepBC+1,cnts0);
                      //fHistoRawYieldTrialBinC0_2[theCase]->SetBinError(itrial,iStepBC+1,ecnts0);
                      //fHistoRawYieldDistBinC0_2[theCase]->Fill(cnts0);
                      if(cnts1<fMinYieldGlobBC_2) fMinYieldGlobBC_2=cnts1;
                      if(cnts1>fMaxYieldGlobBC_2) fMaxYieldGlobBC_2=cnts1;
                      fHistoRawYieldDistBinC1All_2->Fill(cnts1);
                      fHistoRawYieldTrialBinC1All_2->SetBinContent(globBin,cnts1);
                      fHistoRawYieldTrialBinC1All_2->SetBinError(globBin,ecnts1);
                      //fHistoRawYieldTrialBinC1_2[theCase]->SetBinContent(itrial,iStepBC+1,cnts1);
                      //fHistoRawYieldTrialBinC1_2[theCase]->SetBinError(itrial,iStepBC+1,ecnts1);
                      //fHistoRawYieldDistBinC1_2[theCase]->Fill(cnts1);
                    } else { cout << "No more options build for bincounting" << endl; }
                  }
                }
              }
              if (mustDeleteFitter) delete fitter;
              fNtupleMultiTrials->Fill(xnt);
            }
          }
        }
      }
      delete hRebinned;
    }
  }
  return kTRUE;
}

//________________________________________________________________________
void AliHFInvMassMultiTrialFit1::SaveToRoot(TString fileName, TString option) const{
  // save histos in a root file for further analysis
  
  const Int_t nCases=kNBkgFuncCases*kNFitConfCases;
  //const Int_t nCases=2*2;
  TFile outHistos(fileName.Data(),option.Data());
  if (outHistos.IsZombie()) {
    Printf("Could not open file '%s'!", fileName.Data());
    return;
  }
  outHistos.cd();
  fHistoRawYieldTrialAll->Write();
  fHistoSigmaTrialAll->Write();
  fHistoMeanTrialAll->Write();
  fHistoChi2TrialAll->Write();
  fHistoSignifTrialAll->Write();
  if(fSaveBkgVal) {
    fHistoBkgTrialAll->Write();
    fHistoBkgInBinEdgesTrialAll->Write();
  }
  fHistoRawYieldDistAll->Write();
  fHistoRawYieldDistBinC0All->Write();
  fHistoRawYieldTrialBinC0All->Write();
  fHistoRawYieldDistBinC1All->Write();
  fHistoRawYieldTrialBinC1All->Write();
  if(fNumOfnSigmaBinCSteps>1){
    fHistoRawYieldDistBinC0All_2->Write();
    fHistoRawYieldTrialBinC0All_2->Write();
    fHistoRawYieldDistBinC1All_2->Write();
    fHistoRawYieldTrialBinC1All_2->Write();
  }
  /*
   for(Int_t ic=0; ic<nCases; ic++){
   fHistoRawYieldTrial[ic]->Write();
   fHistoSigmaTrial[ic]->Write();
   fHistoMeanTrial[ic]->Write();
   fHistoChi2Trial[ic]->Write();
   fHistoSignifTrial[ic]->Write();
   if(fSaveBkgVal) {
   fHistoBkgTrial[ic]->Write();
   fHistoBkgInBinEdgesTrial[ic]->Write();
   }
   fHistoRawYieldTrialBinC0[ic]->Write();
   fHistoRawYieldDistBinC0[ic]->Write();
   fHistoRawYieldTrialBinC1[ic]->Write();
   fHistoRawYieldDistBinC1[ic]->Write();
   if(fNumOfnSigmaBinCSteps>1){
   fHistoRawYieldTrialBinC0_2[ic]->Write();
   fHistoRawYieldDistBinC0_2[ic]->Write();
   fHistoRawYieldTrialBinC1_2[ic]->Write();
   fHistoRawYieldDistBinC1_2[ic]->Write();
   }
   }
   */
  
  fNtupleMultiTrials->SetDirectory(&outHistos);
  fNtupleMultiTrials->Write();
  outHistos.Close();
}

//________________________________________________________________________
void AliHFInvMassMultiTrialFit1::DrawHistos(TCanvas* cry) const{
  // draw histos
  
  Int_t kNBkgFuncCasesTurnedOn = 0;
  for(Int_t typeb=0; typeb<kNBkgFuncCases; typeb++){
    if(typeb==kExpoBkg && !fUseExpoBkg) continue;
    if(typeb==kLinBkg && !fUseLinBkg) continue;
    if(typeb==kPol2Bkg && !fUsePol2Bkg) continue;
    if(typeb==kPol3Bkg && !fUsePol3Bkg) continue;
    if(typeb==kPol4Bkg && !fUsePol4Bkg) continue;
    if(typeb==kPol5Bkg && !fUsePol5Bkg) continue;
    if(typeb==kPowBkg && !fUsePowLawBkg) continue;
    if(typeb==kPowTimesExpoBkg && !fUsePowLawTimesExpoBkg) continue;
    if(typeb==kNoBkgOnlySig && !fUseNoBackgroundOnlySignal) continue;
    kNBkgFuncCasesTurnedOn++;
  }
  Int_t kNFitConfCasesTurnedOn = 0;
  for(Int_t igs=0; igs<kNFitConfCases; igs++){
    if (igs==kFixSigUpFreeMean && !fUseFixSigUpFreeMean) continue;
    if (igs==kFixSigDownFreeMean && !fUseFixSigDownFreeMean) continue;
    if (igs==kFreeSigFixMean  && !fUseFixedMeanFreeS) continue;
    if (igs==kFreeSigFreeMean  && !fUseFreeS) continue;
    if (igs==kFixSigFreeMean  && !fUseFixSigFreeMean) continue;
    if (igs==kFixSigFixMean   && !fUseFixSigFixMean) continue;
    kNFitConfCasesTurnedOn++;
  }
  
  Int_t totTrials=fNumOfRebinSteps*fNumOfFirstBinSteps*fNumOfLowLimFitSteps*fNumOfUpLimFitSteps;
  gStyle->SetOptStat(0);
  TLine* lSigma = new TLine(0,SigmaRef,kNBkgFuncCasesTurnedOn*kNFitConfCasesTurnedOn*totTrials-0.5,SigmaRef);
  if(fUseNoBackgroundOnlySignal) lSigma = new TLine(0,SigmaRef,totTrials-0.5,SigmaRef);
  lSigma->SetLineWidth(2);
  lSigma->SetLineColor(kRed);
  
  TLine* lRawYield = new TLine(0,RawYieldRef,kNBkgFuncCasesTurnedOn*kNFitConfCasesTurnedOn*totTrials-0.5,RawYieldRef);
  if(fUseNoBackgroundOnlySignal) lRawYield = new TLine(0,RawYieldRef,totTrials-0.5,RawYieldRef);
  lRawYield->SetLineWidth(2);
  lRawYield->SetLineColor(kRed);
  
  Double_t fMassD = 2.286;

  TLine* lPDG = new TLine(0,fMassD,kNBkgFuncCasesTurnedOn*kNFitConfCasesTurnedOn*totTrials-0.5,fMassD);
  lPDG->SetLineWidth(2);
  lPDG->SetLineColor(kMagenta);
  
  cry->Divide(3,2);
  cry->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  fHistoSigmaTrialAll->Draw();
  lSigma->Draw();
  cry->cd(5);
  gPad->SetTickx();
  gPad->SetTicky();
  fHistoMeanTrialAll->Draw();
  lPDG->Draw();
  cry->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  cry->cd(3)->SetTopMargin(0.12);
  fHistoChi2TrialAll->Draw("p");
  if(fUseNoBackgroundOnlySignal) fHistoChi2TrialAll->SetMaximum(150);
  cry->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  
  Double_t rawmin = 1.;
  Double_t rawmax = -1.;
  rawmin = RawYieldRef*(1-0.6);
  rawmax = RawYieldRef*(1+0.6);
  if(fUseNoBackgroundOnlySignal){ rawmin = RawYieldRef*(1-0.15); rawmax = RawYieldRef*(1+0.15);}
  
  TLegend* leg = new TLegend(0.5,0.65,0.89,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.043);
  leg->AddEntry(lSigma,"Central value","l");
  leg->AddEntry(fHistoRawYieldTrialAll,"Fit method","lpe");
  if(fNumOfnSigmaBinCSteps>0) leg->AddEntry(fHistoRawYieldTrialBinC1All,Form("Bin counting (%0.1f#sigma)",fnSigmaBinCSteps[0]),"lpe");
  if(fNumOfnSigmaBinCSteps>1) leg->AddEntry(fHistoRawYieldTrialBinC1All_2,Form("Bin counting (%0.1f#sigma)",fnSigmaBinCSteps[1]),"lpe");
  
  cry->cd(1)->SetTopMargin(0.12);
  fHistoRawYieldTrialBinC1All->Draw();
  if(fNumOfnSigmaBinCSteps>1) fHistoRawYieldTrialBinC1All_2->Draw("same");
  fHistoRawYieldTrialAll->Draw("same");
  lRawYield->Draw();
  fHistoRawYieldTrialBinC1All->GetYaxis()->SetRangeUser(rawmin-0.25*rawmin,rawmax+0.35*rawmax);
  
  TLine** lBkg = new TLine*[(const int)(kNBkgFuncCasesTurnedOn*kNFitConfCasesTurnedOn)];
  for(int ii = 1; ii < kNBkgFuncCasesTurnedOn*kNFitConfCasesTurnedOn; ii++){
    lBkg[ii] = new TLine(ii*totTrials,rawmin-0.2*rawmin,ii*totTrials,rawmax+0.3*rawmax);
    lBkg[ii]->SetLineColor(kGray+2);
    lBkg[ii]->SetLineStyle(2);
    lBkg[ii]->Draw();
  }
  leg->Draw("same");
  
  cry->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  cry->cd(2)->SetTopMargin(0.12);;
  fHistoRawYieldDistBinC1All->GetXaxis()->SetRangeUser(rawmin, rawmax);
  fHistoRawYieldDistBinC1All->Draw();
  if(fNumOfnSigmaBinCSteps>1){
    fHistoRawYieldDistBinC1All_2->Draw();
    fHistoRawYieldDistBinC1All->Draw("same");
  }
  fHistoRawYieldDistAll->Draw("same");
  Int_t rebin = fHistoRawYieldDistAll->GetMaximumBin()/250 + 1;
  fHistoRawYieldDistAll->Rebin(rebin);
  fHistoRawYieldDistBinC1All->Rebin(rebin);
  fHistoRawYieldDistBinC1All_2->Rebin(rebin);
  fHistoRawYieldDistBinC1All->GetYaxis()->SetRangeUser(0.,fHistoRawYieldDistAll->GetMaximum()*1.5);
  fHistoRawYieldDistBinC1All->GetXaxis()->SetRangeUser(rawmin-0.25*rawmin,rawmax+0.35*rawmax);
  
  TLine* lRawRef = new TLine(RawYieldRef,0,RawYieldRef,fHistoRawYieldDistAll->GetMaximum()*1.5);
  lRawRef->SetLineColor(kRed);
  lRawRef->SetLineWidth(2);
  lRawRef->Draw("same");
  
  TPaveText* stats;
  TPaveText* statsbc;
  TPaveText* statsbc_2;
  TPaveText* statsdiffbc;
  
  stats=new TPaveText(0.15,0.65,0.44,0.85,"NDC");
  stats->SetTextSize(0.043);
  stats->SetFillColor(0);
  stats->SetFillStyle(0);
  stats->SetBorderSize(0);
  stats->SetTextFont(42);
  stats->SetTextColor(kBlue+1);
  stats->AddText(Form("mean = %0.1f",fHistoRawYieldDistAll->GetMean()));
//  stats->AddText(Form("RMS = %0.1f (%0.1f%%)",fHistoRawYieldDistAll->GetRMS(),(fHistoRawYieldDistAll->GetRMS())/RawYieldRef * 100));
  stats->AddText(Form("RMS = %0.1f (%0.1f%%)",fHistoRawYieldDistAll->GetRMS(),(fHistoRawYieldDistAll->GetRMS())/fHistoRawYieldDistAll->GetMean() * 100));
  stats->AddText(Form("#frac{max-min}{#sqrt{12}} = %.1f (%0.1f%%)",(fMaxYieldGlob-fMinYieldGlob)/TMath::Sqrt(12), ((fMaxYieldGlob-fMinYieldGlob)/TMath::Sqrt(12) )/RawYieldRef * 100 ));
  
  statsbc=new TPaveText(0.6,0.65-0.35*0,0.89,0.85-0.35*0,"NDC");
  statsbc->SetTextSize(0.043);
  statsbc->SetFillColor(0);
  statsbc->SetFillStyle(0);
  statsbc->SetBorderSize(0);
  statsbc->SetTextFont(42);
  statsbc->SetTextColor(kGreen+2);
  statsbc->AddText(Form("mean = %0.1f",fHistoRawYieldDistBinC1All->GetMean()));
//  statsbc->AddText(Form("RMS = %0.1f (%0.1f%%)",fHistoRawYieldDistBinC1All->GetRMS(),(fHistoRawYieldDistBinC1All->GetRMS())/RawYieldRef * 100));
  statsbc->AddText(Form("RMS = %0.1f (%0.1f%%)",fHistoRawYieldDistBinC1All->GetRMS(),(fHistoRawYieldDistBinC1All->GetRMS())/fHistoRawYieldDistBinC1All->GetMean() * 100));
  statsbc->AddText(Form("#frac{max-min}{#sqrt{12}} = %0.1f (%0.1f%%)",(fMaxYieldGlobBC-fMinYieldGlobBC)/TMath::Sqrt(12), ((fMaxYieldGlobBC-fMinYieldGlobBC)/TMath::Sqrt(12) )/RawYieldRef * 100 ));
  
  statsbc_2=new TPaveText(0.6,0.65-0.35*1,0.89,0.85-0.35*1,"NDC");
  statsbc_2->SetTextSize(0.043);
  statsbc_2->SetFillColor(0);
  statsbc_2->SetFillStyle(0);
  statsbc_2->SetBorderSize(0);
  statsbc_2->SetTextFont(42);
  statsbc_2->SetTextColor(kOrange+1);
  statsbc_2->AddText(Form("mean = %0.1f",fHistoRawYieldDistBinC1All_2->GetMean()));
  statsbc_2->AddText(Form("RMS = %0.1f (%0.1f%%)",fHistoRawYieldDistBinC1All_2->GetRMS(),(fHistoRawYieldDistBinC1All_2->GetRMS())/RawYieldRef * 100));
  statsbc_2->AddText(Form("#frac{max-min}{#sqrt{12}} = %0.1f (%0.1f%%)",(fMaxYieldGlobBC_2-fMinYieldGlobBC_2)/TMath::Sqrt(12), ((fMaxYieldGlobBC_2-fMinYieldGlobBC_2)/TMath::Sqrt(12) )/RawYieldRef * 100 ));
  
  statsdiffbc=new TPaveText(0.6,0.65-0.35*1,0.89,0.85-0.35*1,"NDC");
  statsdiffbc->SetTextSize(0.04);
  statsdiffbc->SetFillColor(0);
  statsdiffbc->SetFillStyle(0);
  statsdiffbc->SetBorderSize(0);
  statsdiffbc->SetTextFont(42);
  statsdiffbc->SetTextColor(kOrange+1);
  statsdiffbc->AddText(Form("Diff mean (BC/fit) = (%0.1f%%)",fHistoRawYieldDistBinC1All->GetMean()/fHistoRawYieldDistAll->GetMean() * 100. - 100));
  
  stats->Draw("same");
  if(fNumOfnSigmaBinCSteps>0) statsbc->Draw("same");
  if(fNumOfnSigmaBinCSteps>1) statsbc_2->Draw("same");
  //if(/*fUseNoBackgroundOnlySignal &&*/ fNumOfnSigmaBinCSteps<= 1) statsdiffbc->Draw("same");
  
}
//________________________________________________________________________
Bool_t AliHFInvMassMultiTrialFit1::DoFitWithPol3Bkg(TH1F* histoToFit, Double_t  hmin, Double_t  hmax,
                                                    Int_t iCase){
  //
  
  TH1F *hCutTmp=(TH1F*)histoToFit->Clone("hCutTmp");
  for(Int_t ib=1; ib<=hCutTmp->GetNbinsX(); ib++){
    Double_t xc=hCutTmp->GetBinCenter(ib);
    if(xc>(fMassD-5.*fSigmaGausMC) && xc<(fMassD+5.*fSigmaGausMC)){
      hCutTmp->SetBinContent(ib,0.);
      hCutTmp->SetBinError(ib,0.);
    }
  }
  
  hCutTmp->Fit("pol2","E0","",hmin,hmax);
  TF1* f2=(TF1*)hCutTmp->GetListOfFunctions()->FindObject("pol2");
  TF1* f3=new TF1("myPol3","pol3");
  for(Int_t i=0; i<3;i++) f3->SetParameter(i,f2->GetParameter(i));
  hCutTmp->Fit(f3,"E0","",hmin,hmax);
  Double_t quickCount=0.;
  for(Int_t ib=1; ib<=histoToFit->GetNbinsX(); ib++){
    Double_t xc=hCutTmp->GetBinCenter(ib);
    if(xc>(fMassD-3.*fSigmaGausMC) && xc<(fMassD+3.*fSigmaGausMC)){
      quickCount+=(histoToFit->GetBinContent(ib)-f3->Eval(xc));
    }
  }
  TF1* fSB=new TF1("fSB","[0]*1./(TMath::Sqrt(2.*TMath::Pi())*[2])*TMath::Exp(-(x-[1])*(x-[1])/(2.*[2]*[2]))+[3]+[4]*x+[5]*x*x+[6]*x*x*x");
  fSB->SetParameter(0,quickCount);
  fSB->SetParameter(1,fMassD);
  fSB->SetParameter(2,fSigmaGausMC);
  for(Int_t j=0; j<4; j++) fSB->SetParameter(j+3,f3->GetParameter(j));
  if(iCase==0) fSB->FixParameter(2,fSigmaGausMC);
  else if(iCase==1) fSB->FixParameter(2,fSigmaGausMC*(1.+fSigmaMCVariation));
  else if(iCase==2) fSB->FixParameter(2,fSigmaGausMC*(1.-fSigmaMCVariation));
  else if(iCase==4){
    fSB->FixParameter(1,fMassD);
    fSB->FixParameter(2,fSigmaGausMC);
  } else if(iCase==5){
    fSB->FixParameter(1,fMassD);
  }
  histoToFit->Fit(fSB,"ME0","",hmin,hmax);
  // quality cuts
  if(fSB->GetParError(0)<0.01*fSB->GetParameter(0)) return kFALSE;
  if(fSB->GetParError(0)>0.6*fSB->GetParameter(0)) return kFALSE;
  
  delete hCutTmp;
  return kTRUE;
}
//__________________________________________________________________________________
void AliHFInvMassMultiTrialFit1::SetTemplatesForReflections(const TH1F *hr, const TH1F *hs) {
  /// signal and reflection templates
  if(fhTemplSign) delete fhTemplSign;
  if(fhTemplRefl) delete fhTemplRefl;
  fhTemplRefl=(TH1F*)hr->Clone("hTemplRefl");
  fhTemplSign=(TH1F*)hs->Clone("hTemplSign");
  return;
}

