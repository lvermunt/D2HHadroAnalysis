const Int_t colors[] = {kAzure-2, kRed+1, kGreen+3, kMagenta-4, kGray+3};
const Int_t colorsFill[] = {kAzure-8, kRed-6, kGreen-5, kMagenta-3, kGray+1};
const Int_t colorsFillLight[] = {kAzure-9, kRed-9, kGreen-8, kMagenta-3, kGray+1};
Int_t fillstyle = 3145;

void SetStyleHisto(TH1D *h);
Int_t GetEmptyMarker(Int_t m);
void SetErrorXTGraph(TGraphAsymmErrors* gr, Double_t errx);
TGraph* extract_TAMU(TString filn);
void RescaleForFiducialAcceptance(TH1D* h);

const Int_t nptbins = 6;
Float_t ptbinsfl[nptbins+1] = {0, 2, 4, 8, 12, 16, 24};

const Int_t nptbins2 = 4;
Float_t ptbinsfl2[nptbins+1] = {0, 4, 8, 12, 18};

const Int_t nptbinsits2 = 2;
Float_t ptbinsflits2[nptbins+1] = {4, 8, 12};
Double_t ptbinsdbits2[nptbins+1] = {4, 8, 12};

const Int_t nptbinsits3 = 4;
Float_t ptbinsflits3[nptbins+1] = {2, 4, 8, 12, 16};
Double_t ptbinsdbits3[nptbins+2] = {2, 4, 8, 12, 16, 24};

Double_t yieldsyst2[nptbinsits3] = {0., 0.06, 0.06, 0.};
Double_t yieldsyst3[nptbinsits3] = {0.06, 0.06, 0.06, 0.06};
Double_t selsyst[nptbinsits3] = {0.08, 0.06, 0.04, 0.04};
Double_t pidsyst[nptbinsits3] = {0., 0., 0., 0.};
Double_t ptshapesyst[nptbinsits3] = {0.01, 0.04, 0.02, 0.01};
Double_t trackingsyst[nptbinsits3] = {0.082, 0.082, 0.082, 0.082};
Double_t brunc = 0.083;
Double_t taaunc = 0.019;

Double_t centTAMU[nptbinsits3] = {2.15415, 1.59321, 0.885321, 0.605254};

void drawLowpTBin(Int_t opt = 1){
  TFile* feff2 = new TFile("theory/pythia_Bs_ITS2.root");
  TFile* feff3 = new TFile("theory/pythia_Bs_ITS3.root");
if(opt == 1){
  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}",nptbins,ptbinsfl);
  SetStyleHisto(href);

  TH1D* hgen = (TH1D*)feff2->Get("h_gen_pr");
  TH1D* hreco2 = (TH1D*)feff2->Get("h_sel_pr");
  TH1D* hreco3 = (TH1D*)feff3->Get("h_sel_pr");

  SetStyleHisto(hgen);
  SetStyleHisto(hreco2);
  SetStyleHisto(hreco3);

  hgen->SetLineColor(colors[4]);
  hgen->SetMarkerColor(colors[4]);
  hgen->SetMarkerStyle(20);

  hreco2->SetLineColor(colors[0]);
  hreco2->SetMarkerColor(colors[0]);
  hreco2->SetMarkerStyle(20);

  hreco3->SetLineColor(colors[1]);
  hreco3->SetMarkerColor(colors[1]);
  hreco3->SetMarkerStyle(20);

  TCanvas* ceff = new TCanvas("ceff", "ceff", 450, 400);
  ceff->cd();
  href->Draw("");
  href->GetXaxis()->SetRangeUser(0,5.);
  href->GetYaxis()->SetRangeUser(1,100000.);
  hgen->Draw("same hist");
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  hreco2->Draw("same hist");
  hreco3->Draw("same hist");

  TLegend* leg = new TLegend(0.36, 0.24, 0.78, 0.42, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(24); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hgen, "Generated", "lm");
  leg->AddEntry(hreco2, "Reco. ITS2", "lm");
  leg->AddEntry(hreco3, "Reco. ITS3", "lm");
  leg->Draw();
} else {
  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});(Acc #times eff.) #times 2#it{y}_{fid}",nptbins,ptbinsfl);
  SetStyleHisto(href);

  TH1D* hgen = (TH1D*)feff2->Get("h_gen_pr");
  TH1D* hreco2 = (TH1D*)feff2->Get("h_sel_pr");
  TH1D* hreco3 = (TH1D*)feff3->Get("h_sel_pr");

  hreco2->Divide(hreco2,hgen,1,1,"B");
  hreco3->Divide(hreco3,hgen,1,1,"B");

  SetStyleHisto(hreco2);
  SetStyleHisto(hreco3);

  hreco2->SetLineColor(colors[0]);
  hreco2->SetMarkerColor(colors[0]);
  hreco2->SetMarkerStyle(20);

  hreco3->SetLineColor(colors[1]);
  hreco3->SetMarkerColor(colors[1]);
  hreco3->SetMarkerStyle(20);

  TCanvas* ceff = new TCanvas("ceff", "ceff", 450, 400);
  ceff->cd();
  href->Draw("");
  href->GetYaxis()->SetRangeUser(0.0001,1.);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  hreco2->Draw("same hist");
  hreco3->Draw("same hist");

  TLegend* leg = new TLegend(0.36, 0.24, 0.78, 0.36, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(24); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hreco2, "Reco. ITS2", "lm");
  leg->AddEntry(hreco3, "Reco. ITS3", "lm");
  leg->Draw();
}
}

void drawEfficiency(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_1505_v2/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_1505_v2/"){
  TString patheff2 = path2 + "signal_background_efficiency.root";
  TString patheff3 = path3 + "signal_background_efficiency.root";

  TFile* feff2 = new TFile(patheff2.Data());
  TFile* feff3 = new TFile(patheff3.Data());

  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});Acceptance #times Efficiency",nptbins2,ptbinsfl2);
  SetStyleHisto(href);
  href->GetYaxis()->SetTitleOffset(1.);

  TH1D* heff2 = (TH1D*)feff2->Get("hEffPythia");
  RescaleForFiducialAcceptance(heff2);
  heff2->SetBinContent(1,0.00001);
  heff2->SetBinError(1,0.);
  heff2->SetBinContent(6,0.00001);
  heff2->SetBinError(6,0.);
  TH1D* heff3 = (TH1D*)feff3->Get("hEffPythia");
  RescaleForFiducialAcceptance(heff3);
  heff3->SetBinContent(1,0.00001);
  heff3->SetBinError(1,0.);
  heff3->SetBinContent(6,0.00001);
  heff3->SetBinError(6,0.);

  heff2->SetTitle(";#it{p}_{T} (GeV/#it{c});Acceptance #times Efficiency");
  heff3->SetTitle(";#it{p}_{T} (GeV/#it{c});Acceptance #times Efficiency");

  heff2->SetLineColor(colors[0]);
  heff2->SetMarkerColor(colors[0]);
  heff2->SetMarkerStyle(20);

  heff3->SetLineColor(colors[1]);
  heff3->SetMarkerColor(colors[1]);
  heff3->SetMarkerStyle(21);

  SetStyleHisto(heff2);
  SetStyleHisto(heff3);

  TH1D* hempty2=(TH1D*)heff2->Clone("hempty2");
  hempty2->SetMarkerStyle(GetEmptyMarker(heff2->GetMarkerStyle()));
  hempty2->SetLineColor(colors[0]);
  hempty2->SetMarkerColor(1);
  TH1D* hempty3=(TH1D*)heff3->Clone("hempty3");
  hempty3->SetMarkerStyle(GetEmptyMarker(heff3->GetMarkerStyle()));
  hempty3->SetLineColor(colors[1]);
  hempty3->SetMarkerColor(1);

  TCanvas* ceff = new TCanvas("ceff", "ceff", 450, 400);
  ceff->cd();
  href->Draw("");
  href->GetYaxis()->SetRangeUser(0.005,1.);
  heff2->Draw("same ep");
  hempty2->Draw("same");
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  heff3->Draw("same ep");
  hempty3->Draw("same");

  TLegend* leg = new TLegend(0.65, 0.24, 0.85, 0.38, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(22); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(heff2, "ITS2", "pm");
  leg->AddEntry(heff3, "ITS3", "pm");
  leg->Draw();

  TLegend* leg2CE=(TLegend*)leg->Clone();
  TList* lCE=leg2CE->GetListOfPrimitives();
  for(Int_t j=0; j<lCE->GetEntries(); j++){
    TLegendEntry* e1=(TLegendEntry*)lCE->At(j);
    TString label=e1->GetLabel();
    if(label.Contains("ITS2")) e1->SetObject(hempty2);
    if(label.Contains("ITS3")) e1->SetObject(hempty3);
    e1->SetLabel(" ");
  }
  leg2CE->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(20);
  TLatex info2; info2.SetNDC(); info2.SetTextFont(43); info2.SetTextSize(17);
  info1.DrawLatex(0.14, 0.84, "ALICE Upgrade Projection");
  info2.DrawLatex(0.14, 0.78, "0^{ }#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info2.DrawLatex(0.14, 0.725, "#font[12]{L}_{int} = 10 nb^{-1}");
  info2.DrawLatex(0.71, 0.84, "B_{s}^{0} #rightarrow D_{s}^{#minus} #pi^{#plus}");
}

void drawBackground(Bool_t allbins = kFALSE, TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_1505_v2/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_1505_v2/"){
  TString pathbkg2 = path2 + "signal_background_efficiency.root";
  TString pathbkg3 = path3 + "signal_background_efficiency.root";

  TFile* fbkg2 = new TFile(pathbkg2.Data());
  TFile* fbkg3 = new TFile(pathbkg3.Data());

  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fbkg2->Get("grBkgFONLL");
  TGraphAsymmErrors* grBkgHIJING2 = (TGraphAsymmErrors*)fbkg2->Get("grBkgHIJING");

  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fbkg3->Get("grBkgFONLL");
  TGraphAsymmErrors* grBkgHIJING3 = (TGraphAsymmErrors*)fbkg3->Get("grBkgHIJING");

  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});#it{B} (3#sigma)",nptbins,ptbinsfl);
  SetStyleHisto(href);

  grBkgFONLL2->SetLineColor(colors[0]);
  grBkgFONLL2->SetMarkerColor(colors[0]);
  grBkgFONLL2->SetMarkerStyle(20);
  grBkgHIJING2->SetFillStyle(0); grBkgHIJING2->SetFillColor(kWhite); grBkgHIJING2->SetLineWidth(2); grBkgHIJING2->SetLineColor(colors[0]);
  //grBkgHIJING2->SetFillColor(colorsFill[0]); grBkgHIJING2->SetLineWidth(2); grBkgHIJING2->SetLineColor(colorsFill[0]); grBkgHIJING2->SetFillStyle(fillstyle);
  SetErrorXTGraph(grBkgHIJING2, 0.75);

  grBkgFONLL3->SetLineColor(colors[1]);
  grBkgFONLL3->SetMarkerColor(colors[1]);
  grBkgFONLL3->SetMarkerStyle(21);
  grBkgHIJING3->SetFillStyle(0); grBkgHIJING3->SetFillColor(kWhite); grBkgHIJING3->SetLineWidth(2); grBkgHIJING3->SetLineColor(colors[1]);
  //grBkgHIJING3->SetFillColor(colorsFill[1]); grBkgHIJING3->SetLineWidth(2); grBkgHIJING3->SetLineColor(colorsFill[1]); grBkgHIJING3->SetFillStyle(fillstyle);
  SetErrorXTGraph(grBkgHIJING3, 0.75);

  if(!allbins){
    grBkgFONLL2->RemovePoint(5);
    grBkgFONLL2->RemovePoint(4);
    grBkgFONLL2->RemovePoint(1);
    grBkgFONLL3->RemovePoint(5);
  }
  //grBkgHIJING2->RemovePoint(5);
  //grBkgHIJING2->RemovePoint(4);
  //grBkgHIJING2->RemovePoint(1);
  grBkgFONLL2->RemovePoint(0);
  grBkgHIJING2->RemovePoint(0);

  //grBkgHIJING3->RemovePoint(5);
  grBkgFONLL3->RemovePoint(0);
  grBkgHIJING3->RemovePoint(0);

  TH1D* hempty2=(TH1D*)grBkgFONLL2->Clone("hempty2");
  hempty2->SetMarkerStyle(GetEmptyMarker(grBkgFONLL2->GetMarkerStyle()));
  hempty2->SetLineColor(colors[0]);
  hempty2->SetMarkerColor(1);
  TH1D* hempty3=(TH1D*)grBkgFONLL3->Clone("hempty3");
  hempty3->SetMarkerStyle(GetEmptyMarker(grBkgFONLL3->GetMarkerStyle()));
  hempty3->SetLineColor(colors[1]);
  hempty3->SetMarkerColor(1);

  //grBkgHIJING2->SetPointEYhigh(2, 2*grBkgHIJING2->GetErrorYlow(1));
  //grBkgHIJING2->SetPointEYhigh(3, grBkgHIJING2->GetErrorYlow(2));

  //grBkgHIJING3->SetPointEYhigh(2, grBkgHIJING3->GetErrorYlow(1));
  //grBkgHIJING3->SetPointEYlow(3, 0.75*grBkgHIJING3->GetErrorYhigh(2));
  //grBkgHIJING3->SetPointEYhigh(3, 0.75*grBkgHIJING3->GetErrorYhigh(2));
  //grBkgHIJING3->SetPointEYhigh(4, grBkgHIJING3->GetErrorYlow(3));

  TCanvas* cbkg = new TCanvas("cbkg", "cbkg", 450, 400);
  cbkg->cd();
  href->Draw("");
  href->GetYaxis()->SetRangeUser(500,2000000);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  grBkgFONLL2->Draw("same ep");
  grBkgHIJING2->Draw("same2");
  hempty2->Draw("same p");
  grBkgFONLL3->Draw("same ep");
  grBkgHIJING3->Draw("same2");
  hempty3->Draw("same p");

  TLegend* leg = new TLegend(0.16, 0.24, 0.58, 0.36, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(24); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(grBkgFONLL2, "ITS2", "plm");
  leg->AddEntry(grBkgFONLL3, "ITS3", "plm");
  leg->Draw();

  TLegend* leg2CE=(TLegend*)leg->Clone();
  TList* lCE=leg2CE->GetListOfPrimitives();
  for(Int_t j=0; j<lCE->GetEntries(); j++){
    TLegendEntry* e1=(TLegendEntry*)lCE->At(j);
    TString label=e1->GetLabel();
    if(label.Contains("ITS2")) e1->SetObject(hempty2);
    if(label.Contains("ITS3")) e1->SetObject(hempty3);
    e1->SetLabel(" ");
  }
  leg2CE->Draw();

  TGraphAsymmErrors* grforleg = (TGraphAsymmErrors*)grBkgFONLL2->Clone("grforleg");
  grforleg->SetLineColor(colors[4]);
  TLegend* legsyst = new TLegend(0.16, 0.14, 0.58, 0.22, 0, "NDC");
  legsyst->SetTextFont(43); legsyst->SetTextSize(14); legsyst->SetFillColor(0); legsyst->SetLineColor(0);
  legsyst->AddEntry(grforleg,"Syst. from FONLL D-meson predictions","e");
  legsyst->AddEntry(grforleg,"Syst. from bkg. shape parametrisation","f");
  legsyst->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  info1.DrawLatex(0.485, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.485, 0.785, "Centrality 0^{ }#font[122]{-}10%");
  info1.DrawLatex(0.485, 0.73, "#font[12]{L}_{int} = 10 nb^{-1}");
}

void drawSignal(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_1505_v2/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_1505_v2/"){
  TString pathsig2 = path2 + "signal_background_efficiency.root";
  TString pathsig3 = path3 + "signal_background_efficiency.root";

  TFile* fsig2 = new TFile(pathsig2.Data());
  TFile* fsig3 = new TFile(pathsig3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsig2->Get("grSigPythia");
  TGraphAsymmErrors* grSigFONLL2 = (TGraphAsymmErrors*)fsig2->Get("grSigFONLL");

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsig3->Get("grSigPythia");
  TGraphAsymmErrors* grSigFONLL3 = (TGraphAsymmErrors*)fsig3->Get("grSigFONLL");

  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});#it{S} (3#sigma)",nptbins,ptbinsfl);
  SetStyleHisto(href);

  grSigPYTHIA2->SetLineColor(colors[0]);
  grSigPYTHIA2->SetMarkerColor(colors[0]);
  grSigPYTHIA2->SetMarkerStyle(20);
  grSigFONLL2->SetFillStyle(0); grSigFONLL2->SetFillColor(kWhite); grSigFONLL2->SetLineWidth(2); grSigFONLL2->SetLineColor(colors[0]);
  //grSigFONLL2->SetFillColor(colorsFill[0]); grSigFONLL2->SetLineWidth(2); grSigFONLL2->SetLineColor(colorsFill[0]); grSigFONLL2->SetFillStyle(fillstyle);
  SetErrorXTGraph(grSigFONLL2, 0.75);

  grSigPYTHIA3->SetLineColor(colors[1]);
  grSigPYTHIA3->SetMarkerColor(colors[1]);
  grSigPYTHIA3->SetMarkerStyle(21);
  grSigFONLL3->SetFillStyle(0); grSigFONLL3->SetFillColor(kWhite); grSigFONLL3->SetLineWidth(2); grSigFONLL3->SetLineColor(colors[1]);
  //grSigFONLL3->SetFillColor(colorsFill[1]); grSigFONLL3->SetLineWidth(2); grSigFONLL3->SetLineColor(colorsFill[1]); grSigFONLL3->SetFillStyle(fillstyle);
  SetErrorXTGraph(grSigFONLL3, 0.75);

  grSigPYTHIA2->RemovePoint(0);
  grSigFONLL2->RemovePoint(0);
  grSigPYTHIA3->RemovePoint(0);
  grSigFONLL3->RemovePoint(0);

  TH1D* hempty2=(TH1D*)grSigPYTHIA2->Clone("hempty2");
  hempty2->SetMarkerStyle(GetEmptyMarker(grSigPYTHIA2->GetMarkerStyle()));
  hempty2->SetLineColor(colors[0]);
  hempty2->SetMarkerColor(1);
  TH1D* hempty3=(TH1D*)grSigPYTHIA3->Clone("hempty3");
  hempty3->SetMarkerStyle(GetEmptyMarker(grSigPYTHIA3->GetMarkerStyle()));
  hempty3->SetLineColor(colors[1]);
  hempty3->SetMarkerColor(1);

  TCanvas* csig = new TCanvas("csig", "csig", 450, 400);
  csig->cd();
  href->Draw("");
  href->GetYaxis()->SetRangeUser(50,10000);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  grSigPYTHIA2->Draw("same ep");
  grSigFONLL2->Draw("same2");
  hempty2->Draw("same p");
  grSigPYTHIA3->Draw("same ep");
  grSigFONLL3->Draw("same2");
  hempty3->Draw("same p");

  TLegend* leg = new TLegend(0.16, 0.24, 0.58, 0.36, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(24); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(grSigPYTHIA2, "ITS2", "plm");
  leg->AddEntry(grSigPYTHIA3, "ITS3", "plm");
  leg->Draw();

  TLegend* leg2CE=(TLegend*)leg->Clone();
  TList* lCE=leg2CE->GetListOfPrimitives();
  for(Int_t j=0; j<lCE->GetEntries(); j++){
    TLegendEntry* e1=(TLegendEntry*)lCE->At(j);
    TString label=e1->GetLabel();
    if(label.Contains("ITS2")) e1->SetObject(hempty2);
    if(label.Contains("ITS3")) e1->SetObject(hempty3);
    e1->SetLabel(" ");
  }
  leg2CE->Draw();

  TGraphAsymmErrors* grforleg = (TGraphAsymmErrors*)grSigPYTHIA2->Clone("grforleg");
  grforleg->SetLineColor(colors[4]);
  TLegend* legsyst = new TLegend(0.16, 0.14, 0.58, 0.22, 0, "NDC");
  legsyst->SetTextFont(43); legsyst->SetTextSize(14); legsyst->SetFillColor(0); legsyst->SetLineColor(0);
  legsyst->AddEntry(grforleg,"Syst. from stat. unc. efficiency","e");
  legsyst->AddEntry(grforleg,"Syst. from FONLL B-meson predictions","f");
  legsyst->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  info1.DrawLatex(0.485, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.485, 0.785, "Centrality 0^{ }#font[122]{-}10%");
  info1.DrawLatex(0.485, 0.73, "#font[12]{L}_{int} = 10 nb^{-1}");
}

void drawSignalOverBackground(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_1505_v2/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_1505_v2/"){
  TString pathsoverb2 = path2 + "signal_background_efficiency.root";
  TString pathsoverb3 = path3 + "signal_background_efficiency.root";

  TFile* fsoverb2 = new TFile(pathsoverb2.Data());
  TFile* fsoverb3 = new TFile(pathsoverb3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsoverb2->Get("grSigPythia");
  TGraphAsymmErrors* grSigFONLL2 = (TGraphAsymmErrors*)fsoverb2->Get("grSigFONLLrel");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsoverb2->Get("grBkgFONLL");
  TGraphAsymmErrors* grBkgHIJING2 = (TGraphAsymmErrors*)fsoverb2->Get("grBkgHIJINGrel");
  TGraphAsymmErrors* grSoverBkgFONLL2 = new TGraphAsymmErrors(nptbinsits2);
  TGraphAsymmErrors* grSoverBkgMC2 = new TGraphAsymmErrors(nptbinsits2);

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsoverb3->Get("grSigPythia");
  TGraphAsymmErrors* grSigFONLL3 = (TGraphAsymmErrors*)fsoverb3->Get("grSigFONLLrel");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsoverb3->Get("grBkgFONLL");
  TGraphAsymmErrors* grBkgHIJING3 = (TGraphAsymmErrors*)fsoverb3->Get("grBkgHIJINGrel");
  TGraphAsymmErrors* grSoverBkgFONLL3 = new TGraphAsymmErrors(nptbinsits3);
  TGraphAsymmErrors* grSoverBkgMC3 = new TGraphAsymmErrors(nptbinsits3);

  grBkgFONLL2->RemovePoint(5);
  grSigPYTHIA2->RemovePoint(5);
  //grBkgFONLL2->RemovePoint(4);
  //grSigPYTHIA2->RemovePoint(4);
  //grBkgFONLL2->RemovePoint(1);
  //grSigPYTHIA2->RemovePoint(1);
  grBkgFONLL2->RemovePoint(0);
  grSigPYTHIA2->RemovePoint(0);

  grSigFONLL2->RemovePoint(5);
  grBkgHIJING2->RemovePoint(5);
  //grSigFONLL2->RemovePoint(4);
  //grBkgHIJING2->RemovePoint(4);
  //grSigFONLL2->RemovePoint(1);
  //grBkgHIJING2->RemovePoint(1);
  grSigFONLL2->RemovePoint(0);
  grBkgHIJING2->RemovePoint(0);

  grBkgFONLL3->RemovePoint(5);
  grSigPYTHIA3->RemovePoint(5);
  grBkgFONLL3->RemovePoint(0);
  grSigPYTHIA3->RemovePoint(0);

  grSigFONLL3->RemovePoint(5);
  grBkgHIJING3->RemovePoint(5);
  grSigFONLL3->RemovePoint(0);
  grBkgHIJING3->RemovePoint(0);

  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});#it{S} / #it{B}",nptbins2,ptbinsfl2);
  SetStyleHisto(href);
  href->GetYaxis()->SetTitleOffset(1.);

  TH1D* hsoverb2 = new TH1D("hsoverb2", ";#it{p}_{T} (GeV/#it{c});#it{S} / #it{B}",nptbinsits3,ptbinsflits3);
  TH1D* hsoverb3 = new TH1D("hsoverb3", ";#it{p}_{T} (GeV/#it{c});#it{S} / #it{B}",nptbinsits3,ptbinsflits3);

  Double_t soverb = 0;
  Double_t systupfonll = 0;
  Double_t systdownfonll = 0;
  Double_t systupmc = 0;
  Double_t systdownmc = 0;
  for(int j = 0; j < nptbinsits3; j++){
    soverb = grSigPYTHIA2->GetY()[j] / grBkgFONLL2->GetY()[j];
    hsoverb2->SetBinContent(j+1, soverb);
    hsoverb2->SetBinError(j+1, 0.00001);

    //if(j < 2){
      systupfonll = soverb * grSigFONLL2->GetErrorYhigh(j); //syst FONLL background is negligible wrt FONLL signal
      systdownfonll = soverb * grSigFONLL2->GetErrorYlow(j); //syst FONLL background is negligible wrt FONLL signal
      systupmc = soverb * grBkgHIJING2->GetErrorYhigh(j); //syst MC signal is negligible wrt MC background
      systdownmc = soverb * grBkgHIJING2->GetErrorYlow(j); //syst MC signal is negligible wrt MC background
      grSoverBkgFONLL2->SetPoint(j, hsoverb2->GetBinCenter(j+1), soverb);
      grSoverBkgFONLL2->SetPointError(j, 0.75, 0.75, systdownfonll, systupfonll);
      grSoverBkgMC2->SetPoint(j, hsoverb2->GetBinCenter(j+1), soverb);
      grSoverBkgMC2->SetPointError(j, 0.75, 0.75, systdownmc, systupmc);
     //}

    soverb = grSigPYTHIA3->GetY()[j] / grBkgFONLL3->GetY()[j];
    hsoverb3->SetBinContent(j+1, soverb);
    hsoverb3->SetBinError(j+1, 0.00001);

    systupfonll = soverb * grSigFONLL3->GetErrorYhigh(j); //syst FONLL background is negligible wrt FONLL signal
    systdownfonll = soverb * grSigFONLL3->GetErrorYlow(j); //syst FONLL background is negligible wrt FONLL signal
    systupmc = soverb * grBkgHIJING3->GetErrorYhigh(j); //syst MC signal is negligible wrt MC background
    systdownmc = soverb * grBkgHIJING3->GetErrorYlow(j); //syst MC signal is negligible wrt MC background
    grSoverBkgFONLL3->SetPoint(j, hsoverb3->GetBinCenter(j+1), soverb);
    grSoverBkgFONLL3->SetPointError(j, 0.75, 0.75, systdownfonll, systupfonll);
    grSoverBkgMC3->SetPoint(j, hsoverb3->GetBinCenter(j+1), soverb);
    grSoverBkgMC3->SetPointError(j, 0.75, 0.75, systdownmc, systupmc);
  }
  SetStyleHisto(hsoverb2);
  SetStyleHisto(hsoverb3);

  hsoverb2->SetLineColor(colors[0]);
  hsoverb2->SetMarkerColor(colors[0]);
  hsoverb2->SetMarkerStyle(20);

  hsoverb3->SetLineColor(colors[1]);
  hsoverb3->SetMarkerColor(colors[1]);
  hsoverb3->SetMarkerStyle(21);

  grSoverBkgFONLL2->SetLineColor(colors[0]);
  grSoverBkgFONLL2->SetMarkerColor(colors[0]);
  grSoverBkgFONLL2->SetMarkerStyle(20);
  grSoverBkgFONLL2->SetFillColor(colorsFill[0]); grSoverBkgFONLL2->SetLineWidth(2); grSoverBkgFONLL2->SetLineColor(colorsFill[0]); grSoverBkgFONLL2->SetFillStyle(fillstyle);
  SetErrorXTGraph(grSoverBkgFONLL2, 0.3);
  TGraphAsymmErrors* grSoverBkgFONLL2_empty =(TGraphAsymmErrors*)grSoverBkgFONLL2->Clone("grSoverBkgFONLL2_empty");
  grSoverBkgFONLL2_empty->SetFillStyle(0); grSoverBkgFONLL2_empty->SetLineColor(colors[0]); grSoverBkgFONLL2_empty->SetLineWidth(1);

  grSoverBkgMC2->SetLineColor(colors[0]);
  grSoverBkgMC2->SetMarkerColor(colors[0]);
  grSoverBkgMC2->SetMarkerStyle(20);
  grSoverBkgMC2->SetFillStyle(0); grSoverBkgMC2->SetFillColor(kWhite); grSoverBkgMC2->SetLineWidth(2); grSoverBkgMC2->SetLineColor(colors[0]);
  SetErrorXTGraph(grSoverBkgMC2, 0.5);

  grSoverBkgFONLL3->SetLineColor(colors[1]);
  grSoverBkgFONLL3->SetMarkerColor(colors[1]);
  grSoverBkgFONLL3->SetMarkerStyle(21);
  grSoverBkgFONLL3->SetFillColor(colorsFill[1]); grSoverBkgFONLL3->SetLineWidth(2); grSoverBkgFONLL3->SetLineColor(colorsFill[1]); grSoverBkgFONLL3->SetFillStyle(fillstyle);
  SetErrorXTGraph(grSoverBkgFONLL3, 0.3);
  TGraphAsymmErrors* grSoverBkgFONLL3_empty =(TGraphAsymmErrors*)grSoverBkgFONLL3->Clone("grSoverBkgFONLL3_empty");
  grSoverBkgFONLL3_empty->SetFillStyle(0); grSoverBkgFONLL3_empty->SetLineColor(colors[1]); grSoverBkgFONLL3_empty->SetLineWidth(1);

  grSoverBkgMC3->SetLineColor(colors[1]);
  grSoverBkgMC3->SetMarkerColor(colors[1]);
  grSoverBkgMC3->SetMarkerStyle(21);
  grSoverBkgMC3->SetFillStyle(0); grSoverBkgMC3->SetFillColor(kWhite); grSoverBkgMC3->SetLineWidth(2); grSoverBkgMC3->SetLineColor(colors[1]);
  SetErrorXTGraph(grSoverBkgMC3, 0.5);

  TH1D* hempty2=(TH1D*)hsoverb2->Clone("hempty2");
  TH1D* hempty2_v2=(TH1D*)hsoverb2->Clone("hempty2_v2");
  TH1D* hempty2_v3=(TH1D*)hsoverb2->Clone("hempty2_v3");
  hempty2->SetMarkerStyle(GetEmptyMarker(hsoverb2->GetMarkerStyle()));
  hempty2_v2->SetMarkerStyle(GetEmptyMarker(hsoverb2->GetMarkerStyle()));
  hempty2->SetLineColor(colors[0]);
  hempty2->SetMarkerColor(1);
  hempty2_v3->SetMarkerColor(colorsFill[0]);
  TH1D* hempty3=(TH1D*)hsoverb3->Clone("hempty3");
  hempty3->SetMarkerStyle(GetEmptyMarker(hsoverb3->GetMarkerStyle()));
  hempty3->SetLineColor(colors[1]);
  hempty3->SetMarkerColor(1);

  TGraphAsymmErrors* grSoverBkgFONLL2_v2 = (TGraphAsymmErrors*)grSoverBkgFONLL2->Clone("grSoverBkgFONLL2_v2");
  TGraphAsymmErrors* grSoverBkgFONLL2_v2_empty = (TGraphAsymmErrors*)grSoverBkgFONLL2_empty->Clone("grSoverBkgFONLL2_v2_empty");
  TGraphAsymmErrors* grSoverBkgMC2_v2 = (TGraphAsymmErrors*)grSoverBkgMC2->Clone("grSoverBkgMC2_v2");
  grSoverBkgFONLL2->RemovePoint(3);
  grSoverBkgFONLL2->RemovePoint(0);
  grSoverBkgFONLL2_empty->RemovePoint(3);
  grSoverBkgFONLL2_empty->RemovePoint(0);
  grSoverBkgMC2->RemovePoint(3);
  grSoverBkgMC2->RemovePoint(0);
  grSoverBkgFONLL2_v2->RemovePoint(2);
  grSoverBkgFONLL2_v2->RemovePoint(1);
  grSoverBkgFONLL2_v2_empty->RemovePoint(2);
  grSoverBkgFONLL2_v2_empty->RemovePoint(1);
  grSoverBkgMC2_v2->RemovePoint(2);
  grSoverBkgMC2_v2->RemovePoint(1);

  hempty2_v2->SetBinContent(2,-1);
  hempty2_v2->SetBinError(2,0.);
  hempty2_v2->SetBinContent(3,-1);
  hempty2_v2->SetBinError(3,0.);
  hempty2_v3->SetBinContent(2,-1);
  hempty2_v3->SetBinError(2,0.);
  hempty2_v3->SetBinContent(3,-1);
  hempty2_v3->SetBinError(3,0.);

  TCanvas* csoverb = new TCanvas("csoverb", "csoverb", 450, 400);
  csoverb->cd();
  href->Draw("");
  href->GetYaxis()->SetRangeUser(0.0005,0.195);
  grSoverBkgFONLL2->Draw("same2");
  grSoverBkgFONLL2_empty->Draw("same2");
  grSoverBkgFONLL2_v2->Draw("same2");
  grSoverBkgFONLL2_v2_empty->Draw("same2");
  hsoverb2->Draw("same ep");
  grSoverBkgMC2->Draw("same2");
  grSoverBkgMC2_v2->Draw("same2");
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  hempty2->Draw("same p");
  hempty2_v3->Draw("same p");
  hempty2_v2->Draw("same p");
  grSoverBkgFONLL3->Draw("same2");
  grSoverBkgFONLL3_empty->Draw("same2");
  hsoverb3->Draw("same ep");
  grSoverBkgMC3->Draw("same2");
  hempty3->Draw("same p");

  TLegend* leg = new TLegend(0.5647321, 0.25, 0.7633929, 0.39, 0, "NDC");
  //TLegend* leg = new TLegend(0.155, 0.52, 0.35, 0.66, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(22); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hsoverb2, "ITS2", "pm");
  leg->AddEntry(hsoverb3, "ITS3", "pm");
  leg->Draw();

  TLegend* leg2CE=(TLegend*)leg->Clone();
  TList* lCE=leg2CE->GetListOfPrimitives();
  for(Int_t j=0; j<lCE->GetEntries(); j++){
    TLegendEntry* e1=(TLegendEntry*)lCE->At(j);
    TString label=e1->GetLabel();
    if(label.Contains("ITS2")) e1->SetObject(hempty2);
    if(label.Contains("ITS3")) e1->SetObject(hempty3);
    e1->SetLabel(" ");
  }
  leg2CE->Draw();

  TGraphAsymmErrors* grSoverBkgFONLL_leg =(TGraphAsymmErrors*)grSoverBkgFONLL2->Clone("grSoverBkgFONLL_leg");
  grSoverBkgFONLL_leg->SetFillColor(colorsFill[4]); grSoverBkgFONLL_leg->SetLineColor(colorsFill[4]);
  TGraphAsymmErrors* grSoverBkgMC_leg =(TGraphAsymmErrors*)grSoverBkgMC2->Clone("grSoverBkgMC_leg");
  grSoverBkgMC_leg->SetLineColor(colors[4]);
  TLegend* legsyst = new TLegend(0.13, 0.125, 0.40, 0.205, 0, "NDC");
  legsyst->SetTextFont(43); legsyst->SetTextSize(14); legsyst->SetFillColor(0); legsyst->SetLineColor(0);
  legsyst->AddEntry(grSoverBkgFONLL_leg, "Unc. signal estimation","f");
  legsyst->AddEntry(grSoverBkgMC_leg, "Unc. background estimation","f");
  legsyst->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(20);
  TLatex info2; info2.SetNDC(); info2.SetTextFont(43); info2.SetTextSize(17);
  info1.DrawLatex(0.14, 0.84, "ALICE Upgrade Projection");
  info2.DrawLatex(0.14, 0.78, "0^{ }#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info2.DrawLatex(0.14, 0.725, "#font[12]{L}_{int} = 10 nb^{-1}");
  info2.DrawLatex(0.71, 0.84, "B_{s}^{0} #rightarrow D_{s}^{#minus} #pi^{#plus}");
}

void drawSignificance(Bool_t allbins = kFALSE, TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_1505_v2/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_1505_v2/"){
  TString pathsignf2 = path2 + "signal_background_efficiency.root";
  TString pathsignf3 = path3 + "signal_background_efficiency.root";

  TFile* fsignf2 = new TFile(pathsignf2.Data());
  TFile* fsignf3 = new TFile(pathsignf3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsignf2->Get("grSigPythia");
  TGraphAsymmErrors* grSigFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grSigFONLL");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgFONLL");
  TGraphAsymmErrors* grBkgHIJING2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgHIJING");
  TGraphAsymmErrors* grSgnfFONLL2 = new TGraphAsymmErrors(nptbinsits2);
  TGraphAsymmErrors* grSgnfMC2 = new TGraphAsymmErrors(nptbinsits2);

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsignf3->Get("grSigPythia");
  TGraphAsymmErrors* grSigFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grSigFONLL");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgFONLL");
  TGraphAsymmErrors* grBkgHIJING3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgHIJING");
  TGraphAsymmErrors* grSgnfFONLL3 = new TGraphAsymmErrors(nptbinsits3);
  TGraphAsymmErrors* grSgnfMC3 = new TGraphAsymmErrors(nptbinsits3);

  if(!allbins){
    grBkgFONLL2->RemovePoint(5);
    grSigPYTHIA2->RemovePoint(5);
    //grBkgFONLL2->RemovePoint(4);
    //grSigPYTHIA2->RemovePoint(4);
    //grBkgFONLL2->RemovePoint(1);
    //grSigPYTHIA2->RemovePoint(1);
    grBkgFONLL2->RemovePoint(0);
    grSigPYTHIA2->RemovePoint(0);

    grBkgHIJING2->RemovePoint(5);
    grSigFONLL2->RemovePoint(5);
    //grBkgHIJING2->RemovePoint(4);
    //grSigFONLL2->RemovePoint(4);
    //grBkgHIJING2->RemovePoint(1);
    //grSigFONLL2->RemovePoint(1);
    grBkgHIJING2->RemovePoint(0);
    grSigFONLL2->RemovePoint(0);

    grBkgFONLL3->RemovePoint(5);
    grSigPYTHIA3->RemovePoint(5);
    grBkgFONLL3->RemovePoint(0);
    grSigPYTHIA3->RemovePoint(0);

    grBkgHIJING3->RemovePoint(5);
    grSigFONLL3->RemovePoint(5);
    grBkgHIJING3->RemovePoint(0);
    grSigFONLL3->RemovePoint(0);
  }

  TH1D* href;
  if(!allbins) href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});Expected significance",nptbins2,ptbinsfl2);
  else  href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});Expected significance",nptbins,ptbinsfl);
  SetStyleHisto(href);

  TH1D* hsignf2;
  TH1D* hsignf3;
  if(allbins){
    hsignf2 = new TH1D("hsignf2", ";#it{p}_{T} (GeV/#it{c});Expected significance",nptbins,ptbinsfl);
    hsignf3 = new TH1D("hsignf3", ";#it{p}_{T} (GeV/#it{c});Expected significance",nptbins,ptbinsfl);
  } else {
    hsignf2 = new TH1D("hsignf2", ";#it{p}_{T} (GeV/#it{c});Expected significance",nptbinsits3,ptbinsflits3);
    hsignf3 = new TH1D("hsignf3", ";#it{p}_{T} (GeV/#it{c});Expected significance",nptbinsits3,ptbinsflits3);
  }

  Double_t significance = 0;
  Double_t errSigSq = 0;
  Double_t errBkgSq = 0;
  Double_t sigPlusBkg = 0;
  Double_t background = 0;
  Double_t signal = 0;

  Double_t systupfonll = 0;
  Double_t systdownfonll = 0;
  Double_t systupmc = 0;
  Double_t systdownmc = 0;

  Int_t nforloop = nptbinsits3;
  if(allbins) nforloop = nptbins;
  for(int j = 0; j < nforloop; j++){
    signal = grSigPYTHIA2->GetY()[j];
    background = grBkgFONLL2->GetY()[j];
    sigPlusBkg = grSigPYTHIA2->GetY()[j] + grBkgFONLL2->GetY()[j];
    significance = grSigPYTHIA2->GetY()[j] / TMath::Sqrt( grSigPYTHIA2->GetY()[j] + grBkgFONLL2->GetY()[j]);

    hsignf2->SetBinContent(j+1, significance);
    hsignf2->SetBinError(j+1, 0.00001);
    //if(j < 2){
      //syst FONLL background is negligible wrt FONLL signal
      errSigSq = grSigFONLL2->GetErrorYhigh(j) * grSigFONLL2->GetErrorYhigh(j); errBkgSq = 0 * 0;
      systupfonll = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
      errSigSq = grSigFONLL2->GetErrorYlow(j) * grSigFONLL2->GetErrorYlow(j); errBkgSq = 0 * 0;
      systdownfonll = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

      //syst MC signal is negligible wrt MC background
      errSigSq = 0 * 0; errBkgSq = grBkgHIJING2->GetErrorYhigh(j) * grBkgHIJING2->GetErrorYhigh(j);
      systupmc = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
      errSigSq = 0 * 0; errBkgSq = grBkgHIJING2->GetErrorYlow(j) * grBkgHIJING2->GetErrorYlow(j);
      systdownmc = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

      grSgnfFONLL2->SetPoint(j, hsignf2->GetBinCenter(j+1), significance);
      grSgnfFONLL2->SetPointError(j, 0.75, 0.75, systdownfonll, systupfonll);
      grSgnfMC2->SetPoint(j, hsignf2->GetBinCenter(j+1), significance);
      grSgnfMC2->SetPointError(j, 0.75, 0.75, systdownmc, systupmc);
    //}

    signal = grSigPYTHIA3->GetY()[j];
    background = grBkgFONLL3->GetY()[j];
    sigPlusBkg = grSigPYTHIA3->GetY()[j] + grBkgFONLL3->GetY()[j];
    significance = grSigPYTHIA3->GetY()[j] / TMath::Sqrt( grSigPYTHIA3->GetY()[j] + grBkgFONLL3->GetY()[j]);

    hsignf3->SetBinContent(j+1, significance);
    hsignf3->SetBinError(j+1, 0.00001);

    //syst FONLL background is negligible wrt FONLL signal
    errSigSq = grSigFONLL3->GetErrorYhigh(j) * grSigFONLL3->GetErrorYhigh(j); errBkgSq = 0 * 0;
    systupfonll = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
    errSigSq = grSigFONLL3->GetErrorYlow(j) * grSigFONLL3->GetErrorYlow(j); errBkgSq = 0 * 0;
    systdownfonll = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

    //syst MC signal is negligible wrt MC background
    errSigSq = 0 * 0; errBkgSq = grBkgHIJING3->GetErrorYhigh(j) * grBkgHIJING3->GetErrorYhigh(j);
    systupmc = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
    errSigSq = 0 * 0; errBkgSq = grBkgHIJING3->GetErrorYlow(j) * grBkgHIJING3->GetErrorYlow(j);
    systdownmc = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

    grSgnfFONLL3->SetPoint(j, hsignf3->GetBinCenter(j+1), significance);
    grSgnfFONLL3->SetPointError(j, 0.75, 0.75, systdownfonll, systupfonll);
    grSgnfMC3->SetPoint(j, hsignf3->GetBinCenter(j+1), significance);
    grSgnfMC3->SetPointError(j, 0.75, 0.75, systdownmc, systupmc);
  }
  SetStyleHisto(hsignf2);
  SetStyleHisto(hsignf3);

  hsignf2->SetLineColor(colors[0]);
  hsignf2->SetMarkerColor(colors[0]);
  hsignf2->SetMarkerStyle(20);

  hsignf3->SetLineColor(colors[1]);
  hsignf3->SetMarkerColor(colors[1]);
  hsignf3->SetMarkerStyle(21);

  grSgnfFONLL2->SetLineColor(colors[0]);
  grSgnfFONLL2->SetMarkerColor(colors[0]);
  grSgnfFONLL2->SetMarkerStyle(20);
  grSgnfFONLL2->SetFillColor(colorsFill[0]); grSgnfFONLL2->SetLineWidth(2); grSgnfFONLL2->SetLineColor(colorsFill[0]); grSgnfFONLL2->SetFillStyle(fillstyle);
  SetErrorXTGraph(grSgnfFONLL2, 0.3);
  TGraphAsymmErrors* grSgnfFONLL2_empty =(TGraphAsymmErrors*)grSgnfFONLL2->Clone("grSgnfFONLL2_empty");
  grSgnfFONLL2_empty->SetFillStyle(0); grSgnfFONLL2_empty->SetLineColor(colors[0]); grSgnfFONLL2_empty->SetLineWidth(1);

  grSgnfMC2->SetLineColor(colors[0]);
  grSgnfMC2->SetMarkerColor(colors[0]);
  grSgnfMC2->SetMarkerStyle(20);
  grSgnfMC2->SetFillStyle(0); grSgnfMC2->SetFillColor(kWhite); grSgnfMC2->SetLineWidth(2); grSgnfMC2->SetLineColor(colors[0]);
  SetErrorXTGraph(grSgnfMC2, 0.5);

  grSgnfFONLL3->SetLineColor(colors[1]);
  grSgnfFONLL3->SetMarkerColor(colors[1]);
  grSgnfFONLL3->SetMarkerStyle(21);
  grSgnfFONLL3->SetFillColor(colorsFill[1]); grSgnfFONLL3->SetLineWidth(2); grSgnfFONLL3->SetLineColor(colorsFill[1]); grSgnfFONLL3->SetFillStyle(fillstyle);
  SetErrorXTGraph(grSgnfFONLL3, 0.3);
  TGraphAsymmErrors* grSgnfFONLL3_empty =(TGraphAsymmErrors*)grSgnfFONLL3->Clone("grSgnfFONLL3_empty");
  grSgnfFONLL3_empty->SetFillStyle(0); grSgnfFONLL3_empty->SetLineColor(colors[1]); grSgnfFONLL3_empty->SetLineWidth(1);

  grSgnfMC3->SetLineColor(colors[1]);
  grSgnfMC3->SetMarkerColor(colors[1]);
  grSgnfMC3->SetMarkerStyle(21);
  grSgnfMC3->SetFillStyle(0); grSgnfMC3->SetFillColor(kWhite); grSgnfMC3->SetLineWidth(2); grSgnfMC3->SetLineColor(colors[1]);
  SetErrorXTGraph(grSgnfMC3, 0.5);

  TH1D* hempty2=(TH1D*)hsignf2->Clone("hempty2");
  TH1D* hempty2_v2=(TH1D*)hsignf2->Clone("hempty2_v2");
  TH1D* hempty2_v3=(TH1D*)hsignf2->Clone("hempty2_v3");
  hempty2->SetMarkerStyle(GetEmptyMarker(hsignf2->GetMarkerStyle()));
  hempty2_v2->SetMarkerStyle(GetEmptyMarker(hsignf2->GetMarkerStyle()));
  hempty2->SetLineColor(colors[0]);
  hempty2->SetMarkerColor(1);
  hempty2_v3->SetMarkerColor(colorsFill[0]);
  TH1D* hempty3=(TH1D*)hsignf3->Clone("hempty3");
  hempty3->SetMarkerStyle(GetEmptyMarker(hsignf3->GetMarkerStyle()));
  hempty3->SetLineColor(colors[1]);
  hempty3->SetMarkerColor(1);

  TGraphAsymmErrors* grSgnfFONLL2_v2 = (TGraphAsymmErrors*)grSgnfFONLL2->Clone("grSgnfFONLL2_v2");
  TGraphAsymmErrors* grSgnfFONLL2_v2_empty = (TGraphAsymmErrors*)grSgnfFONLL2_empty->Clone("grSgnfFONLL2_v2_empty");
  TGraphAsymmErrors* grSgnfMC2_v2 = (TGraphAsymmErrors*)grSgnfMC2->Clone("grSgnfMC2_v2");
  grSgnfFONLL2->RemovePoint(3);
  grSgnfFONLL2->RemovePoint(0);
  grSgnfFONLL2_empty->RemovePoint(3);
  grSgnfFONLL2_empty->RemovePoint(0);
  grSgnfMC2->RemovePoint(3);
  grSgnfMC2->RemovePoint(0);
  grSgnfFONLL2_v2->RemovePoint(2);
  grSgnfFONLL2_v2->RemovePoint(1);
  grSgnfFONLL2_v2_empty->RemovePoint(2);
  grSgnfFONLL2_v2_empty->RemovePoint(1);
  grSgnfMC2_v2->RemovePoint(2);
  grSgnfMC2_v2->RemovePoint(1);

  hempty2_v2->SetBinContent(2,-1);
  hempty2_v2->SetBinError(2,0.);
  hempty2_v2->SetBinContent(3,-1);
  hempty2_v2->SetBinError(3,0.);
  hempty2_v3->SetBinContent(2,-1);
  hempty2_v3->SetBinError(2,0.);
  hempty2_v3->SetBinContent(3,-1);
  hempty2_v3->SetBinError(3,0.);

  TLine* l;
  l = new TLine(0,3,18,3);
  l->SetLineStyle(9);

  TCanvas* csignf = new TCanvas("csignf", "csignf", 450, 400);
  csignf->cd();
  href->Draw("");
  href->GetYaxis()->SetRangeUser(0.,15);
  if(!allbins) href->GetYaxis()->SetRangeUser(0.,20);
  //if(!allbins) gPad->SetLogy();
  //if(!allbins) href->GetYaxis()->SetRangeUser(0.4,40.);
  gPad->SetTickx();
  gPad->SetTicky();
  if(!allbins) l->Draw();
  if(!allbins) grSgnfFONLL2->Draw("same2");
  if(!allbins) grSgnfFONLL2_empty->Draw("same2");
  if(!allbins) grSgnfFONLL2_v2->Draw("same2");
  if(!allbins) grSgnfFONLL2_v2_empty->Draw("same2");
  hsignf2->Draw("same ep");
  if(!allbins) grSgnfMC2->Draw("same2");
  if(!allbins) grSgnfMC2_v2->Draw("same2");
  hempty2->Draw("same p");
  if(!allbins) hempty2_v3->Draw("same2");
  if(!allbins) hempty2_v2->Draw("same2");
  if(!allbins) grSgnfFONLL3->Draw("same2");
  if(!allbins) grSgnfFONLL3_empty->Draw("same2");
  hsignf3->Draw("same ep");
  if(!allbins) grSgnfMC3->Draw("same2");
  hempty3->Draw("same p");

  TLegend* leg = new TLegend(0.5647321,0.608,0.7633929,0.7493333, 0, "NDC");
  //TLegend* leg = new TLegend(0.155, 0.52, 0.35, 0.66, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(22); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hsignf2, "ITS2", "pm");
  leg->AddEntry(hsignf3, "ITS3", "pm");
  leg->Draw();

  TLegend* leg2CE=(TLegend*)leg->Clone();
  TList* lCE=leg2CE->GetListOfPrimitives();
  for(Int_t j=0; j<lCE->GetEntries(); j++){
    TLegendEntry* e1=(TLegendEntry*)lCE->At(j);
    TString label=e1->GetLabel();
    if(label.Contains("ITS2")) e1->SetObject(hempty2);
    if(label.Contains("ITS3")) e1->SetObject(hempty3);
    e1->SetLabel(" ");
  }
  leg2CE->Draw();

  TGraphAsymmErrors* grSgnfFONLL_leg =(TGraphAsymmErrors*)grSgnfFONLL2->Clone("grSgnfFONLL_leg");
  grSgnfFONLL_leg->SetFillColor(colorsFill[4]); grSgnfFONLL_leg->SetLineColor(colorsFill[4]);
  TGraphAsymmErrors* grSgnfMC_leg =(TGraphAsymmErrors*)grSgnfMC2->Clone("grSgnfMC_leg");
  grSgnfMC_leg->SetLineColor(colors[4]);
  TLegend* legsyst = new TLegend(0.4642857,0.512,0.734375,0.592, 0, "NDC");
  legsyst->SetTextFont(43); legsyst->SetTextSize(13); legsyst->SetFillColor(0); legsyst->SetLineColor(0);
  //legsyst->AddEntry(grSgnfMC_leg,"Unc. on background estimation (MC)","f");
  //legsyst->AddEntry(grSgnfFONLL_leg,"Unc. on signal estimation (FONLL)","f");
  legsyst->AddEntry(grSgnfFONLL_leg,"Unc. signal estimation","f");
  legsyst->AddEntry(grSgnfMC_leg,"Unc. background estimation","f");

  if(!allbins) legsyst->Draw();

  if(!allbins){
    TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(20);
    TLatex info2; info2.SetNDC(); info2.SetTextFont(43); info2.SetTextSize(17);
    info1.DrawLatex(0.14, 0.84, "ALICE Upgrade Projection");
    info2.DrawLatex(0.14, 0.78, "0^{ }#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
    info2.DrawLatex(0.14, 0.725, "#font[12]{L}_{int} = 10 nb^{-1}");
    info2.DrawLatex(0.71, 0.84, "B_{s}^{0} #rightarrow D_{s}^{#minus} #pi^{#plus}");
  } else {
    TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
    info1.DrawLatex(0.155, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
    info1.DrawLatex(0.155, 0.785, "Centrality 0^{ }#font[122]{-}10%");
    info1.DrawLatex(0.155, 0.73, "#font[12]{L}_{int} = 10 nb^{-1}");
  }
}

void drawSignificanceRatio(Bool_t allbins = kFALSE, TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_1505_v2/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_1505_v2/"){
  TString pathsignf2 = path2 + "signal_background_efficiency.root";
  TString pathsignf3 = path3 + "signal_background_efficiency.root";

  TFile* fsignf2 = new TFile(pathsignf2.Data());
  TFile* fsignf3 = new TFile(pathsignf3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsignf2->Get("grSigPythia");
  TGraphAsymmErrors* grSigFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grSigFONLL");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgFONLL");
  TGraphAsymmErrors* grBkgHIJING2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgHIJING");
  TGraphAsymmErrors* grSgnfFONLL2 = new TGraphAsymmErrors(nptbinsits2);
  TGraphAsymmErrors* grSgnfMC2 = new TGraphAsymmErrors(nptbinsits2);

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsignf3->Get("grSigPythia");
  TGraphAsymmErrors* grSigFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grSigFONLL");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgFONLL");
  TGraphAsymmErrors* grBkgHIJING3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgHIJING");
  TGraphAsymmErrors* grSgnfFONLL3 = new TGraphAsymmErrors(nptbinsits3);
  TGraphAsymmErrors* grSgnfMC3 = new TGraphAsymmErrors(nptbinsits3);

  TH1D* href;
  if(allbins) href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});Expected significance ratio ITS3 / ITS2",nptbins,ptbinsfl);
  else href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});Expected significance ratio ITS3 / ITS2",nptbins2,ptbinsfl2);
  SetStyleHisto(href);

  TH1D* hsignf2;
  TH1D* hsignf3;
  hsignf2 = new TH1D("hsignf2", ";#it{p}_{T} (GeV/#it{c});Expected significance ratio ITS3 / ITS2",nptbins,ptbinsfl);
  hsignf3 = new TH1D("hsignf3", ";#it{p}_{T} (GeV/#it{c});Expected significance ratio ITS3 / ITS2",nptbins,ptbinsfl);

  Double_t significance = 0;
  Double_t errSigSq = 0;
  Double_t errBkgSq = 0;
  Double_t sigPlusBkg = 0;
  Double_t background = 0;
  Double_t signal = 0;

  Double_t systupfonll = 0;
  Double_t systdownfonll = 0;
  Double_t systupmc = 0;
  Double_t systdownmc = 0;

  Int_t nforloop = nptbinsits3;
  Int_t nforstart = 0;
  if(allbins){ nforloop = nptbins; nforstart = 1;}
  for(int j = 0; j < nptbins; j++){
    signal = grSigPYTHIA2->GetY()[j];
    background = grBkgFONLL2->GetY()[j];
    sigPlusBkg = grSigPYTHIA2->GetY()[j] + grBkgFONLL2->GetY()[j];
    significance = grSigPYTHIA2->GetY()[j] / TMath::Sqrt( grSigPYTHIA2->GetY()[j] + grBkgFONLL2->GetY()[j]);

    hsignf2->SetBinContent(j+1, significance);
    hsignf2->SetBinError(j+1, 0.00001);

    //syst FONLL background is negligible wrt FONLL signal
    errSigSq = grSigFONLL2->GetErrorYhigh(j) * grSigFONLL2->GetErrorYhigh(j); errBkgSq = 0 * 0;
    systupfonll = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
    errSigSq = grSigFONLL2->GetErrorYlow(j) * grSigFONLL2->GetErrorYlow(j); errBkgSq = 0 * 0;
    systdownfonll = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

    //syst MC signal is negligible wrt MC background
    errSigSq = 0 * 0; errBkgSq = grBkgHIJING2->GetErrorYhigh(j) * grBkgHIJING2->GetErrorYhigh(j);
    systupmc = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
    errSigSq = 0 * 0; errBkgSq = grBkgHIJING2->GetErrorYlow(j) * grBkgHIJING2->GetErrorYlow(j);
    systdownmc = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

    grSgnfFONLL2->SetPoint(j, hsignf2->GetBinCenter(j+1), significance);
    grSgnfFONLL2->SetPointError(j, 0.75, 0.75, systdownfonll, systupfonll);
    grSgnfMC2->SetPoint(j, hsignf2->GetBinCenter(j+1), significance);
    grSgnfMC2->SetPointError(j, 0.75, 0.75, systdownmc, systupmc);

    signal = grSigPYTHIA3->GetY()[j];
    background = grBkgFONLL3->GetY()[j];
    sigPlusBkg = grSigPYTHIA3->GetY()[j] + grBkgFONLL3->GetY()[j];
    significance = grSigPYTHIA3->GetY()[j] / TMath::Sqrt( grSigPYTHIA3->GetY()[j] + grBkgFONLL3->GetY()[j]);

    hsignf3->SetBinContent(j+1, significance);
    hsignf3->SetBinError(j+1, 0.00001);

    //syst FONLL background is negligible wrt FONLL signal
    errSigSq = grSigFONLL3->GetErrorYhigh(j) * grSigFONLL3->GetErrorYhigh(j); errBkgSq = 0 * 0;
    systupfonll = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
    errSigSq = grSigFONLL3->GetErrorYlow(j) * grSigFONLL3->GetErrorYlow(j); errBkgSq = 0 * 0;
    systdownfonll = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

    //syst MC signal is negligible wrt MC background
    errSigSq = 0 * 0; errBkgSq = grBkgHIJING3->GetErrorYhigh(j) * grBkgHIJING3->GetErrorYhigh(j);
    systupmc = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
    errSigSq = 0 * 0; errBkgSq = grBkgHIJING3->GetErrorYlow(j) * grBkgHIJING3->GetErrorYlow(j);
    systdownmc = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);

    grSgnfFONLL3->SetPoint(j, hsignf3->GetBinCenter(j+1), significance);
    grSgnfFONLL3->SetPointError(j, 0.75, 0.75, systdownfonll, systupfonll);
    grSgnfMC3->SetPoint(j, hsignf3->GetBinCenter(j+1), significance);
    grSgnfMC3->SetPointError(j, 0.75, 0.75, systdownmc, systupmc);
  }
  SetStyleHisto(hsignf2);
  SetStyleHisto(hsignf3);

  hsignf3->Divide(hsignf2);
  for(int j = 0; j < nptbins; j++) hsignf3->SetBinError(j+1, 0.00001);

  for(int j = 0; j < nptbins; j++){
    Double_t systfonllup = TMath::Abs(grSgnfFONLL2->GetErrorYhigh(j)/grSgnfFONLL2->GetY()[j] - grSgnfFONLL3->GetErrorYhigh(j)/grSgnfFONLL3->GetY()[j]);
    Double_t systfonlldown = TMath::Abs(grSgnfFONLL2->GetErrorYlow(j)/grSgnfFONLL2->GetY()[j] - grSgnfFONLL3->GetErrorYlow(j)/grSgnfFONLL3->GetY()[j]);
    cout << "FONLL: " << hsignf3->GetBinCenter(j+1) << " " << systfonlldown << " " << systfonllup << endl;

    Double_t systmcup = TMath::Abs(grSgnfMC2->GetErrorYhigh(j)/grSgnfMC2->GetY()[j] - grSgnfMC3->GetErrorYhigh(j)/grSgnfMC3->GetY()[j]);
    Double_t systmcdown = TMath::Abs(grSgnfMC2->GetErrorYlow(j)/grSgnfMC2->GetY()[j] - grSgnfMC3->GetErrorYlow(j)/grSgnfMC3->GetY()[j]);
    cout << "MC: " << hsignf3->GetBinCenter(j+1) << " " << systmcdown << " " << systmcup << endl;
    if(j == 0){ systmcup = 0.03; systmcdown = 0.03;}

    Double_t systup = TMath::Sqrt(systfonllup*systfonllup + systmcup*systmcup);
    Double_t systdown = TMath::Sqrt(systfonlldown*systfonlldown + systmcdown*systmcdown);

    grSgnfMC2->SetPoint(j, hsignf3->GetBinCenter(j+1), hsignf3->GetBinContent(j+1));
    grSgnfMC2->SetPointError(j, 0.5, 0.5, systdown*hsignf3->GetBinContent(j+1), systup*hsignf3->GetBinContent(j+1));
  }

  if(allbins){
    //hsignf3 = (TH1D*)hsignf3->Rebin(5,hsignf3->GetName(),ptbinsdbits3);
    //grSgnfMC2->RemovePoint(0);
  } else {
    hsignf3 = (TH1D*)hsignf3->Rebin(4,hsignf3->GetName(),ptbinsdbits3);
    grSgnfMC2->RemovePoint(5);
    //grSgnfMC2->RemovePoint(4);
    //grSgnfMC2->RemovePoint(1);
    grSgnfMC2->RemovePoint(0);
  }

  hsignf3->SetLineColor(colors[4]);
  hsignf3->SetMarkerColor(colors[4]);
  hsignf3->SetMarkerStyle(20);

  grSgnfMC2->SetLineColor(colors[4]);
  grSgnfMC2->SetMarkerColor(colors[4]);
  grSgnfMC2->SetMarkerStyle(20);
  grSgnfMC2->SetFillStyle(0); grSgnfMC2->SetFillColor(kWhite); grSgnfMC2->SetLineWidth(2); grSgnfMC2->SetLineColor(colors[4]);
  SetErrorXTGraph(grSgnfMC2, 0.5);

  TH1F* hopen = (TH1F*)hsignf3->Clone("hopen");
  hopen->SetMarkerStyle(GetEmptyMarker(hsignf3->GetMarkerStyle()));
  hopen->SetBinContent(2,-1);
  hopen->SetBinContent(3,-1);
  hopen->SetBinError(2,0.);
  hopen->SetBinError(3,0.);
  TH1F* hopen2 = (TH1F*)hopen->Clone("hopen2");
  hopen2->SetMarkerStyle(hsignf3->GetMarkerStyle());
  hopen2->SetMarkerColor(kWhite);

  TCanvas* csignfratio = new TCanvas("csignfratio", "csignfratio", 450, 400);
  csignfratio->cd();
  href->Draw("");
  gPad->SetTickx();
  gPad->SetTicky();
  href->GetYaxis()->SetRangeUser(0.,3.5);
  grSgnfMC2->Draw("same2");
  hsignf3->Draw("same ep");
  hopen2->Draw("same p");
  hopen->Draw("same p");

  if(!allbins){
    TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(20);
    TLatex info2; info2.SetNDC(); info2.SetTextFont(43); info2.SetTextSize(17);
    info1.DrawLatex(0.14, 0.84, "ALICE Upgrade Projection");
    info2.DrawLatex(0.14, 0.78, "0^{ }#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
    info2.DrawLatex(0.14, 0.725, "#font[12]{L}_{int} = 10 nb^{-1}");
    info2.DrawLatex(0.71, 0.84, "B_{s}^{0} #rightarrow D_{s}^{#minus} #pi^{#plus}");
  } else {
    TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
    info1.DrawLatex(0.285, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
    info1.DrawLatex(0.285, 0.785, "Centrality 0^{ }#font[122]{-}10%");
    info1.DrawLatex(0.285, 0.73, "#font[12]{L}_{int} = 10 nb^{-1}");
    info1.DrawLatex(0.71, 0.84, "B_{s}^{0} #rightarrow D_{s}^{#minus} #pi^{#plus}");
  }

  TLatex info2; info2.SetNDC(); info2.SetTextFont(43); info2.SetTextSize(12);
  info2.DrawLatex(0.14, 0.29, "Open markers: Only in reach for ITS3");

  TLine* l;
  if(allbins) l = new TLine(0,1,24,1);
  else l = new TLine(0,1,18,1);
  l->SetLineStyle(9);
  l->Draw();

  TLegend* legsyst = new TLegend(0.13, 0.13, 0.40, 0.17, 0, "NDC");
  legsyst->SetTextFont(43); legsyst->SetTextSize(14); legsyst->SetFillColor(0); legsyst->SetLineColor(0);
  legsyst->AddEntry(grSgnfMC2,"Unc. on signal and background estimation","f");
  legsyst->Draw();
}

void drawRelStatError(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_1505_v2/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_1505_v2/"){
  TString pathsignf2 = path2 + "signal_background_efficiency.root";
  TString pathsignf3 = path3 + "signal_background_efficiency.root";

  TFile* fsignf2 = new TFile(pathsignf2.Data());
  TFile* fsignf3 = new TFile(pathsignf3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsignf2->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgFONLL");

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsignf3->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgFONLL");

  grBkgFONLL2->RemovePoint(5);
  grSigPYTHIA2->RemovePoint(5);
  grBkgFONLL2->RemovePoint(4);
  grSigPYTHIA2->RemovePoint(4);
  grBkgFONLL2->RemovePoint(1);
  grSigPYTHIA2->RemovePoint(1);
  grBkgFONLL2->RemovePoint(0);
  grSigPYTHIA2->RemovePoint(0);

  grBkgFONLL3->RemovePoint(5);
  grSigPYTHIA3->RemovePoint(5);
  grBkgFONLL3->RemovePoint(0);
  grSigPYTHIA3->RemovePoint(0);

  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbins2,ptbinsfl2);
  SetStyleHisto(href);

  TH1D* hsignf2 = new TH1D("hsignf2", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbinsits2,ptbinsflits2);
  TH1D* hsyst2 = new TH1D("hsyst2", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbinsits2,ptbinsflits2);
  TH1D* hsignf3 = new TH1D("hsignf3", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbinsits3,ptbinsflits3);
  TH1D* hsyst = new TH1D("hsyst", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbinsits3,ptbinsflits3);

  for(int j = 0; j < nptbinsits3; j++){
    Double_t syst = TMath::Sqrt(yieldsyst3[j]*yieldsyst3[j] + selsyst[j]*selsyst[j] + pidsyst[j]*pidsyst[j] + ptshapesyst[j]*ptshapesyst[j] + trackingsyst[j]*trackingsyst[j] + taaunc*taaunc);
    Double_t syst2 = TMath::Sqrt(yieldsyst2[j]*yieldsyst2[j] + selsyst[j]*selsyst[j] + pidsyst[j]*pidsyst[j] + ptshapesyst[j]*ptshapesyst[j] + trackingsyst[j]*trackingsyst[j] + taaunc*taaunc);

    if(j == 0 || j == 3){
      //hsignf2->SetBinContent(j+1, 0);
      //hsignf2->SetBinError(j+1, 0.00001);
    } else {
      hsignf2->SetBinContent(j, 1. / (grSigPYTHIA2->GetY()[j-1] / TMath::Sqrt(grSigPYTHIA2->GetY()[j-1] + grBkgFONLL2->GetY()[j-1])));
      hsignf2->SetBinError(j, 0.00001);
      hsyst2->SetBinContent(j, syst2);
      hsyst2->SetBinError(j, 0.00001);
    }
    cout << "CHECK: " << syst << "      " << yieldsyst3[j] << " " << selsyst[j] << " " << pidsyst[j] << " " << ptshapesyst[j] << " " << trackingsyst[j] << " " << taaunc << endl;

    hsignf3->SetBinContent(j+1, 1. / (grSigPYTHIA3->GetY()[j] / TMath::Sqrt(grSigPYTHIA3->GetY()[j] + grBkgFONLL3->GetY()[j])));
    hsignf3->SetBinError(j+1, 0.00001);

    hsyst->SetBinContent(j+1, syst);
    hsyst->SetBinError(j+1, 0.00001);
  }
  SetStyleHisto(hsignf2);
  SetStyleHisto(hsignf3);
  SetStyleHisto(hsyst);
  SetStyleHisto(hsyst2);

  hsignf2->SetLineColor(colors[0]);
  hsignf2->SetMarkerColor(colors[0]);
  hsignf2->SetMarkerStyle(20);
  hsignf2->SetLineWidth(3);

  hsignf3->SetLineColor(colors[1]);
  hsignf3->SetMarkerColor(colors[1]);
  hsignf3->SetMarkerStyle(21);
  hsignf3->SetLineWidth(3);

  hsyst2->SetLineColor(kGray+3);
  hsyst2->SetMarkerColor(kGray+3);
  hsyst2->SetMarkerStyle(21);
  hsyst2->SetLineWidth(2);
  hsyst2->SetLineStyle(9);

  hsyst->SetLineColor(kGray+3);
  hsyst->SetMarkerColor(kGray+3);
  hsyst->SetMarkerStyle(21);
  hsyst->SetLineWidth(3);
  hsyst->SetLineStyle(9);

  TCanvas* csignf = new TCanvas("csignf", "csignf", 450, 400);
  csignf->cd();
  href->Draw("");
  href->GetYaxis()->SetRangeUser(0.,0.55);
  gPad->SetTickx();
  gPad->SetTicky();
  //hsyst2->Draw("same hist");
  hsyst->Draw("same hist");
  hsignf2->Draw("same hist");
  hsignf3->Draw("same hist");

  TLegend* leg = new TLegend(0.18, 0.52, 0.605, 0.67, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(15); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hsignf2, "stat. unc. ITS2", "lm");
  leg->AddEntry(hsignf3, "stat. unc. ITS3", "lm");
  leg->AddEntry(hsyst, "syst. unc.", "lm");
  leg->Draw();

  //TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(20);
  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  TLatex info2; info2.SetNDC(); info2.SetTextFont(43); info2.SetTextSize(17);
  //info1.DrawLatex(0.14, 0.84, "ALICE Upgrade Projection");
  //info2.DrawLatex(0.14, 0.78, "0^{ }#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  //info2.DrawLatex(0.14, 0.725, "#font[12]{L}_{int} = 10 nb^{-1}");
  info1.DrawLatex(0.14, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.14, 0.785, "Centrality 0^{ }#font[122]{-}10%");
  info1.DrawLatex(0.14, 0.73, "#font[12]{L}_{int} = 10 nb^{-1}");
  info2.DrawLatex(0.71, 0.84, "B_{s}^{0} #rightarrow D_{s}^{#minus} #pi^{#plus}");
}

void drawRAA(Bool_t wnonstrange = kTRUE, TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_1505_v2/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_1505_v2/"){
  TString pathsignf2 = path2 + "signal_background_efficiency.root";
  TString pathsignf3 = path3 + "signal_background_efficiency.root";

  TFile* fsignf2 = new TFile(pathsignf2.Data());
  TFile* fsignf3 = new TFile(pathsignf3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsignf2->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgFONLL");

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsignf3->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgFONLL");

  grBkgFONLL2->RemovePoint(5);
  grSigPYTHIA2->RemovePoint(5);
  grBkgFONLL2->RemovePoint(4);
  grSigPYTHIA2->RemovePoint(4);
  grBkgFONLL2->RemovePoint(1);
  grSigPYTHIA2->RemovePoint(1);
  grBkgFONLL2->RemovePoint(0);
  grSigPYTHIA2->RemovePoint(0);

  grBkgFONLL3->RemovePoint(5);
  grSigPYTHIA3->RemovePoint(5);
  grBkgFONLL3->RemovePoint(0);
  grSigPYTHIA3->RemovePoint(0);

  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",nptbins2,ptbinsfl2);
  TH1D* hraa2 = new TH1D("hraa2", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",nptbinsits2,ptbinsflits2);
  TH1D* hraa3 = new TH1D("hraa3", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbinsits3,ptbinsflits3);
  TGraphErrors* hsyst = new TGraphErrors(nptbinsits3);

  TGraph* grtamubs = (TGraph*)extract_TAMU("theory/input_RAA_TAMU_Bs.txt");
  TGraph* grtamub = (TGraph*)extract_TAMU("theory/input_RAA_TAMU_B.txt");
  grtamubs->SetLineColor(kRed+1);
  grtamubs->SetLineStyle(9);
  grtamubs->SetLineWidth(2);
  grtamub->SetLineColor(kAzure-7);
  grtamub->SetLineStyle(5);
  grtamub->SetLineWidth(2);

  for(int j = 0; j < nptbins2; j++){
    Double_t syst = TMath::Sqrt(yieldsyst3[j]*yieldsyst3[j] + selsyst[j]*selsyst[j] + pidsyst[j]*pidsyst[j] + ptshapesyst[j]*ptshapesyst[j] + trackingsyst[j]*trackingsyst[j] + taaunc*taaunc);
    Double_t pt = ptbinsflits3[j] + 0.5 * (ptbinsflits3[j+1] - ptbinsflits3[j]);

    if(j == 0 || j == 3){
      //hsignf2->SetBinContent(j+1, 0);
      //hsignf2->SetBinError(j+1, 0.00001);
    } else {
      hraa2->SetBinContent(j, grtamubs->Eval(pt));
      hraa2->SetBinError(j, grtamubs->Eval(pt) / (grSigPYTHIA2->GetY()[j-1] / TMath::Sqrt(grSigPYTHIA2->GetY()[j-1] + grBkgFONLL2->GetY()[j-1])));
    }
    hraa3->SetBinContent(j+1, grtamubs->Eval(pt));
    hraa3->SetBinError(j+1, grtamubs->Eval(pt) / (grSigPYTHIA3->GetY()[j] / TMath::Sqrt(grSigPYTHIA3->GetY()[j] + grBkgFONLL3->GetY()[j])));

    hsyst->SetPoint(j, pt, grtamubs->Eval(pt));
    hsyst->SetPointError(j, 0.5, grtamubs->Eval(pt) * syst);
  }
  SetStyleHisto(href);
  SetStyleHisto(hraa2);
  SetStyleHisto(hraa3);

  hraa2->SetLineColor(kGray+1);
  hraa2->SetMarkerColor(kGray+1);
  hraa2->SetMarkerStyle(25);

  hraa3->SetLineColor(kBlack);
  hraa3->SetMarkerColor(kBlack);
  hraa3->SetMarkerStyle(20);
  hraa3->SetMarkerSize(0.92);

  hsyst->SetLineColor(kBlack);
  hsyst->SetMarkerColor(kBlack);
  hsyst->SetFillStyle(0); hsyst->SetFillColor(kWhite); hsyst->SetLineWidth(2); hsyst->SetLineColor(kBlack);

  TLine* l = new TLine(0,1,18,1);
  l->SetLineStyle(9);

  TCanvas* csignf = new TCanvas("csignf", "csignf", 450, 400);
  csignf->cd();
  href->Draw();
  l->Draw();
  href->GetYaxis()->SetRangeUser(0.,3.8);
  grtamubs->Draw("same l");
  if(wnonstrange) grtamub->Draw("same l");
  hraa2->Draw("same ep");
  gPad->SetTickx();
  gPad->SetTicky();
  hraa3->Draw("same ep");
  hsyst->Draw("same2");

  TLegend* leg = new TLegend(0.4631696,0.6053333,0.8671875,0.7173333, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(18); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hraa2, "ITS2", "pm");
  leg->AddEntry(hraa3, "ITS3", "pm");
  leg->Draw();

  TLegend* leg2 = new TLegend(0.4631696,0.412,0.8671875,0.58, 0, "NDC");
  if(!wnonstrange) leg2 = new TLegend(0.4631696,0.468,0.8671875,0.58, 0, "NDC");
  leg2->SetTextFont(43); leg2->SetTextSize(18); leg2->SetFillColor(0); leg2->SetFillStyle(0); leg2->SetLineColor(0);
  leg2->SetHeader("TAMU");
  leg2->AddEntry(grtamubs, "B_{s}^{0}", "l");
  if(wnonstrange) leg2->AddEntry(grtamub, "B (non-strange)", "l");
  leg2->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(20);
  TLatex info2; info2.SetNDC(); info2.SetTextFont(43); info2.SetTextSize(17);
  info1.DrawLatex(0.14, 0.84, "ALICE Upgrade Projection");
  info2.DrawLatex(0.14, 0.78, "0^{ }#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info2.DrawLatex(0.14, 0.725, "#font[12]{L}_{int} = 10 nb^{-1}");
  info2.DrawLatex(0.71, 0.84, "B_{s}^{0} #rightarrow D_{s}^{#minus} #pi^{#plus}");

}

void drawRAA_asNPDs(Bool_t wnonstrange = kTRUE, TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_1505_v2/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_1505_v2/"){
  TString pathsignf2 = path2 + "signal_background_efficiency.root";
  TString pathsignf3 = path3 + "signal_background_efficiency.root";

  TFile* fsignf2 = new TFile(pathsignf2.Data());
  TFile* fsignf3 = new TFile(pathsignf3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsignf2->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgFONLL");

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsignf3->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgFONLL");

  grBkgFONLL2->RemovePoint(5);
  grSigPYTHIA2->RemovePoint(5);
  grBkgFONLL2->RemovePoint(4);
  grSigPYTHIA2->RemovePoint(4);
  grBkgFONLL2->RemovePoint(1);
  grSigPYTHIA2->RemovePoint(1);
  grBkgFONLL2->RemovePoint(0);
  grSigPYTHIA2->RemovePoint(0);

  grBkgFONLL3->RemovePoint(5);
  grSigPYTHIA3->RemovePoint(5);
  grBkgFONLL3->RemovePoint(0);
  grSigPYTHIA3->RemovePoint(0);

  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",nptbins,ptbinsfl);
  TH1D* hraa2 = new TH1D("hraa2", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",nptbinsits2,ptbinsflits2);
  TH1D* hraa3 = new TH1D("hraa3", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbinsits3,ptbinsflits3);
  TGraphErrors* hsyst = new TGraphErrors(nptbinsits3);

  TGraph* grtamubs = (TGraph*)extract_TAMU("theory/input_RAA_TAMU_Bs.txt");
  TGraph* grtamub = (TGraph*)extract_TAMU("theory/input_RAA_TAMU_B.txt");
  grtamubs->SetLineColor(kRed+1);
  grtamubs->SetLineStyle(9);
  grtamubs->SetLineWidth(3);
  grtamub->SetLineColor(kAzure+4);
  grtamub->SetLineStyle(5);
  grtamub->SetLineWidth(3);

  for(int j = 0; j < nptbins2; j++){
    Double_t syst = TMath::Sqrt(yieldsyst3[j]*yieldsyst3[j] + selsyst[j]*selsyst[j] + pidsyst[j]*pidsyst[j] + ptshapesyst[j]*ptshapesyst[j] + trackingsyst[j]*trackingsyst[j] + taaunc*taaunc);
    Double_t pt = ptbinsflits3[j] + 0.5 * (ptbinsflits3[j+1] - ptbinsflits3[j]);

    if(j == 0 || j == 3){
      //hsignf2->SetBinContent(j+1, 0);
      //hsignf2->SetBinError(j+1, 0.00001);
    } else {
      hraa2->SetBinContent(j, grtamubs->Eval(pt));
      hraa2->SetBinError(j, grtamubs->Eval(pt) / (grSigPYTHIA2->GetY()[j-1] / TMath::Sqrt(grSigPYTHIA2->GetY()[j-1] + grBkgFONLL2->GetY()[j-1])));
    }
    hraa3->SetBinContent(j+1, grtamubs->Eval(pt));
    hraa3->SetBinError(j+1, grtamubs->Eval(pt) / (grSigPYTHIA3->GetY()[j] / TMath::Sqrt(grSigPYTHIA3->GetY()[j] + grBkgFONLL3->GetY()[j])));

    hsyst->SetPoint(j, pt, grtamubs->Eval(pt));
    hsyst->SetPointError(j, 0.5, grtamubs->Eval(pt) * syst);
  }
  SetStyleHisto(href);
  SetStyleHisto(hraa2);
  SetStyleHisto(hraa3);

  hraa2->SetLineColor(kGray+1);
  hraa2->SetMarkerColor(kGray+1);
  hraa2->SetMarkerStyle(25);

  hraa3->SetLineColor(kBlack);
  hraa3->SetMarkerColor(kBlack);
  hraa3->SetMarkerStyle(20);
  //hraa3->SetMarkerSize(0.92);

  hsyst->SetLineColor(kBlack);
  hsyst->SetMarkerColor(kBlack);
  hsyst->SetFillStyle(0); hsyst->SetFillColor(kWhite); hsyst->SetLineWidth(2); hsyst->SetLineColor(kBlack);

  TLine* l = new TLine(0,1,20,1);
  l->SetLineStyle(9);
  l->SetLineWidth(1);
  l->SetLineColor(kGray+1);

  gStyle->SetPadRightMargin(0.035);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.035);
  gStyle->SetTitleSize(0.05, "xyz");
  gStyle->SetLabelSize(0.045, "xyz");
  gStyle->SetTitleOffset(1., "x");
  gStyle->SetTitleOffset(1.4, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  TCanvas* csignf = new TCanvas("csignf", "csignf", 500, 500);
  TH1F* hFrame = (TH1F*)csignf->DrawFrame(0., 0.01, 20., 3.75, ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}");
  hFrame->GetYaxis()->SetDecimals();
  l->Draw();
  grtamubs->Draw("same l");
  if(wnonstrange) grtamub->Draw("same l");
  hsyst->Draw("same2");
  hraa2->Draw("same ep");
  hraa3->Draw("same ep");

  TLegend* leg = new TLegend(0.55,0.62,0.9,0.8, 0, "NDC");
  leg->SetTextSize(0.045); leg->SetFillStyle(0); leg->SetBorderSize(0);
  leg->SetHeader("B_{s}^{0} #rightarrow D_{s}^{#minus} #pi^{#plus}");
  leg->AddEntry(hraa2, "ITS2", "p");
  leg->AddEntry(hraa3, "ITS3", "p");
  leg->Draw();

  TLegend* leg2 = new TLegend(0.55,0.50,0.9,0.62, 0, "NDC");
  if(wnonstrange) leg2 = new TLegend(0.55,0.44,0.9,0.62, 0, "NDC");
  leg2->SetTextSize(0.045); leg2->SetFillStyle(0); leg2->SetBorderSize(0);
  leg2->SetHeader("TAMU");
  leg2->AddEntry(grtamubs, "B_{s}^{0}", "l");
  if(wnonstrange) leg2->AddEntry(grtamub, "B (non-strange)", "l");
  leg2->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(42); info1.SetTextSize(0.045); info1.SetTextColor(kBlack);
  TLatex info2; info2.SetNDC(); info2.SetTextFont(43); info2.SetTextSize(17);
  info1.DrawLatex(0.18, 0.90, "ALICE Upgrade Projection");
  info1.DrawLatex(0.18, 0.84, "0#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.18, 0.78, "#font[12]{L}_{int} = 10 nb^{-1}");
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
  h->GetXaxis()->SetTitleOffset(0.88);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetNdivisions(505);

}

Int_t GetEmptyMarker(Int_t m){
  if(m==20) return 24;
  if(m==21) return 25;
  if(m==22) return 26;
  if(m==23) return 32;
  if(m==29) return 24;
  if(m==33) return 27;
  if(m==34) return 28;

  return -99;
}

void SetErrorXTGraph(TGraphAsymmErrors* gr, Double_t errx){

  for(int i = 0; i < gr->GetN(); i++){
    gr->SetPointEXhigh(i, errx);
    gr->SetPointEXlow(i, errx);
  }
}

TGraph* extract_TAMU(TString filn){

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
    g->SetPoint(iPt++,pt,raa);
  }

  g->SetLineColor(kOrange+1);
  g->SetLineWidth(3);
  g->SetLineStyle(1);

  fclose(f);

  return g;
}

void RescaleForFiducialAcceptance(TH1D* h){
  for (Int_t ibin = 1; ibin<=h->GetNbinsX(); ibin++){
    Double_t pt = h->GetBinCenter(ibin);
    Double_t yfid=0.8;
    if(pt<5){
      yfid=-0.2/15*pt*pt+1.9/15*pt+0.5;
    }
    printf("Fiducial acceptance: bin %d (pt = %f) yfid= %f\n",ibin,pt,yfid);
    Double_t aey=h->GetBinContent(ibin);
    Double_t erraey=h->GetBinError(ibin);
    h->SetBinContent(ibin,aey/(2*yfid));
    h->SetBinError(ibin,erraey/(2*yfid));
  }
}