const Int_t colors[] = {kAzure-2, kRed+1, kGreen+3, kMagenta-4, kBlack};
const Int_t colorsFill[] = {kAzure-8, kRed-6, kGreen-5, kMagenta-3, kGray+1};
const Int_t colorsFillLight[] = {kAzure-9, kRed-9, kGreen-8, kMagenta-3, kGray+1};
Int_t fillstyle = 3145;

void SetStyleHisto(TH1D *h);
Int_t GetEmptyMarker(Int_t m);
void SetErrorXTGraph(TGraphAsymmErrors* gr, Double_t errx);
TGraph* extract_TAMU(TString filn);

const Int_t nptbins = 5;
Float_t ptbinsfl[nptbins+1] = {0, 4, 8, 12, 16, 24};

const Int_t nptbins2 = 4;
Float_t ptbinsfl2[nptbins+1] = {0, 4, 8, 12, 16};

Double_t yieldsyst[nptbins2] = {0.1, 0.1, 0.1, 0.1};
Double_t selsyst[nptbins2] = {0.08, 0.06, 0.04, 0.04};
Double_t pidsyst[nptbins2] = {0., 0., 0., 0.};
Double_t ptshapesyst[nptbins2] = {0.05, 0.03, 0.02, 0.};
Double_t trackingsyst[nptbins2] = {0.082, 0.082, 0.082, 0.082};
Double_t brunc = 0.1;

Double_t centTAMU[nptbins2] = {2.15415, 1.59321, 0.885321, 0.605254};

void drawEfficiency(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_0505/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_0505/"){
  TString patheff2 = path2 + "signal_background_efficiency.root";
  TString patheff3 = path3 + "signal_background_efficiency.root";

  TFile* feff2 = new TFile(patheff2.Data());
  TFile* feff3 = new TFile(patheff3.Data());

  TH1D* heff2 = (TH1D*)feff2->Get("hEffPythia");
  TH1D* heff3 = (TH1D*)feff3->Get("hEffPythia");

  heff2->SetTitle(";#it{p}_{T} (GeV/#it{c});(Acc #times eff.) #times 2#it{y}_{fid}");
  heff3->SetTitle(";#it{p}_{T} (GeV/#it{c});(Acc #times eff.) #times 2#it{y}_{fid}");

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
  heff2->Draw("ep");
  heff2->GetYaxis()->SetRangeUser(0.001,1.);
  hempty2->Draw("same");
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  heff3->Draw("same ep");
  hempty3->Draw("same");

  TLegend* leg = new TLegend(0.16, 0.24, 0.58, 0.36, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(24); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(heff2, "ITS2", "lpm");
  leg->AddEntry(heff3, "ITS3", "lpm");
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

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  info1.DrawLatex(0.155, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.155, 0.785, "Centrality 0^{ }#font[122]{-}10%");
  info1.DrawLatex(0.155, 0.73, "#it{L}_{int} = 10 nb^{-1}");

}

void drawBackground(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_0505/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_0505/"){
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

  TH1D* hempty2=(TH1D*)grBkgFONLL2->Clone("hempty2");
  hempty2->SetMarkerStyle(GetEmptyMarker(grBkgFONLL2->GetMarkerStyle()));
  hempty2->SetLineColor(colors[0]);
  hempty2->SetMarkerColor(1);
  TH1D* hempty3=(TH1D*)grBkgFONLL3->Clone("hempty3");
  hempty3->SetMarkerStyle(GetEmptyMarker(grBkgFONLL3->GetMarkerStyle()));
  hempty3->SetLineColor(colors[1]);
  hempty3->SetMarkerColor(1);

  grBkgHIJING2->SetPointEYhigh(0, 3*grBkgHIJING2->GetErrorYlow(0));
  grBkgHIJING2->SetPointEYhigh(1, 2*grBkgHIJING2->GetErrorYlow(1));
  grBkgHIJING2->SetPointEYhigh(2, grBkgHIJING2->GetErrorYlow(2));

  grBkgHIJING3->SetPointEYhigh(0, 2*grBkgHIJING3->GetErrorYlow(0));
  grBkgHIJING3->SetPointEYhigh(1, grBkgHIJING3->GetErrorYlow(1));
  grBkgHIJING3->SetPointEYlow(2, 0.75*grBkgHIJING3->GetErrorYhigh(2));
  grBkgHIJING3->SetPointEYhigh(2, 0.75*grBkgHIJING3->GetErrorYhigh(2));
  grBkgHIJING3->SetPointEYhigh(3, grBkgHIJING3->GetErrorYlow(3));

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
  info1.DrawLatex(0.485, 0.73, "#it{L}_{int} = 10 nb^{-1}");
}

void drawSignal(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_0505/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_0505/"){
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
  legsyst->AddEntry(grforleg,"Syst. from FONLL D-meson predictions","f");
  legsyst->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  info1.DrawLatex(0.485, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.485, 0.785, "Centrality 0^{ }#font[122]{-}10%");
  info1.DrawLatex(0.485, 0.73, "#it{L}_{int} = 10 nb^{-1}");
}

void drawSignalOverBackground(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_0505/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_0505/"){
  TString pathsoverb2 = path2 + "signal_background_efficiency.root";
  TString pathsoverb3 = path3 + "signal_background_efficiency.root";

  TFile* fsoverb2 = new TFile(pathsoverb2.Data());
  TFile* fsoverb3 = new TFile(pathsoverb3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsoverb2->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsoverb2->Get("grBkgFONLL");

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsoverb3->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsoverb3->Get("grBkgFONLL");

  TH1D* hsoverb2 = new TH1D("hsoverb2", ";#it{p}_{T} (GeV/#it{c});#it{S} (3#sigma) / #it{B} (3#sigma)",nptbins,ptbinsfl);
  TH1D* hsoverb3 = new TH1D("hsoverb3", ";#it{p}_{T} (GeV/#it{c});#it{S} (3#sigma) / #it{B} (3#sigma)",nptbins,ptbinsfl);

  for(int j = 0; j < nptbins; j++){
    cout << grSigPYTHIA2->GetY()[j] / grBkgFONLL2->GetY()[j] << endl;
    hsoverb2->SetBinContent(j+1, grSigPYTHIA2->GetY()[j] / grBkgFONLL2->GetY()[j]);
    hsoverb2->SetBinError(j+1, 0.00001);
    cout << grSigPYTHIA3->GetY()[j] / grBkgFONLL3->GetY()[j] << endl;
    hsoverb3->SetBinContent(j+1, grSigPYTHIA3->GetY()[j] / grBkgFONLL3->GetY()[j]);
    hsoverb3->SetBinError(j+1, 0.00001);
  }
  SetStyleHisto(hsoverb2);
  SetStyleHisto(hsoverb3);

  hsoverb2->SetLineColor(colors[0]);
  hsoverb2->SetMarkerColor(colors[0]);
  hsoverb2->SetMarkerStyle(20);

  hsoverb3->SetLineColor(colors[1]);
  hsoverb3->SetMarkerColor(colors[1]);
  hsoverb3->SetMarkerStyle(21);

  TH1D* hempty2=(TH1D*)hsoverb2->Clone("hempty2");
  hempty2->SetMarkerStyle(GetEmptyMarker(hsoverb2->GetMarkerStyle()));
  hempty2->SetLineColor(colors[0]);
  hempty2->SetMarkerColor(1);
  TH1D* hempty3=(TH1D*)hsoverb3->Clone("hempty3");
  hempty3->SetMarkerStyle(GetEmptyMarker(hsoverb3->GetMarkerStyle()));
  hempty3->SetLineColor(colors[1]);
  hempty3->SetMarkerColor(1);

  TCanvas* csoverb = new TCanvas("csoverb", "csoverb", 450, 400);
  csoverb->cd();
  hsoverb2->Draw("ep");
  hsoverb2->GetYaxis()->SetRangeUser(0.0001,0.5);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLogy();
  hempty2->Draw("same p");
  hsoverb3->Draw("same ep");
  hempty3->Draw("same p");

  TLegend* leg = new TLegend(0.16, 0.18, 0.58, 0.30, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(24); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hsoverb2, "ITS2", "plm");
  leg->AddEntry(hsoverb3, "ITS3", "plm");
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

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  info1.DrawLatex(0.155, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.155, 0.785, "Centrality 0^{ }#font[122]{-}10%");
  info1.DrawLatex(0.155, 0.73, "#it{L}_{int} = 10 nb^{-1}");
}

void drawSignificance(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_0505/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_0505/"){
  TString pathsignf2 = path2 + "signal_background_efficiency.root";
  TString pathsignf3 = path3 + "signal_background_efficiency.root";

  TFile* fsignf2 = new TFile(pathsignf2.Data());
  TFile* fsignf3 = new TFile(pathsignf3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsignf2->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgFONLL");

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsignf3->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgFONLL");

  TH1D* hsignf2 = new TH1D("hsignf2", ";#it{p}_{T} (GeV/#it{c});Expected significance",nptbins,ptbinsfl);
  TH1D* hsignf3 = new TH1D("hsignf3", ";#it{p}_{T} (GeV/#it{c});Expected significance",nptbins,ptbinsfl);

  for(int j = 0; j < nptbins; j++){
    Double_t fs020 = 2 * 18.83 / 23.26;
    Double_t fb020 = 2 * 309.7 / 357.3;
    cout << fs020 << " " << fb020 << endl;

    cout << grSigPYTHIA2->GetY()[j] / TMath::Sqrt( grSigPYTHIA2->GetY()[j] + grBkgFONLL2->GetY()[j]) << "    ";
    cout << fs020 * grSigPYTHIA2->GetY()[j] / TMath::Sqrt( fs020 * grSigPYTHIA2->GetY()[j] + fb020 * grBkgFONLL2->GetY()[j]) << endl;
    hsignf2->SetBinContent(j+1, fs020 * grSigPYTHIA2->GetY()[j] / TMath::Sqrt( fs020 * grSigPYTHIA2->GetY()[j] + fb020 * grBkgFONLL2->GetY()[j]));
    hsignf2->SetBinError(j+1, 0.00001);
    cout << fs020 * grSigPYTHIA3->GetY()[j] / TMath::Sqrt( fs020 * grSigPYTHIA3->GetY()[j] + fb020 * grBkgFONLL3->GetY()[j]) << endl;
    hsignf3->SetBinContent(j+1, fs020 * grSigPYTHIA3->GetY()[j] / TMath::Sqrt( fs020 * grSigPYTHIA3->GetY()[j] + fb020 * grBkgFONLL3->GetY()[j]));
    hsignf3->SetBinError(j+1, 0.00001);
  }
  SetStyleHisto(hsignf2);
  SetStyleHisto(hsignf3);

  hsignf2->SetLineColor(colors[0]);
  hsignf2->SetMarkerColor(colors[0]);
  hsignf2->SetMarkerStyle(20);

  hsignf3->SetLineColor(colors[1]);
  hsignf3->SetMarkerColor(colors[1]);
  hsignf3->SetMarkerStyle(21);

  TH1D* hempty2=(TH1D*)hsignf2->Clone("hempty2");
  hempty2->SetMarkerStyle(GetEmptyMarker(hsignf2->GetMarkerStyle()));
  hempty2->SetLineColor(colors[0]);
  hempty2->SetMarkerColor(1);
  TH1D* hempty3=(TH1D*)hsignf3->Clone("hempty3");
  hempty3->SetMarkerStyle(GetEmptyMarker(hsignf3->GetMarkerStyle()));
  hempty3->SetLineColor(colors[1]);
  hempty3->SetMarkerColor(1);

  TCanvas* csignf = new TCanvas("csignf", "csignf", 450, 400);
  csignf->cd();
  hsignf2->Draw("ep");
  hsignf2->GetYaxis()->SetRangeUser(0.,13.);
  gPad->SetTickx();
  gPad->SetTicky();
  hempty2->Draw("same p");
  hsignf3->Draw("same ep");
  hempty3->Draw("same p");

  TLegend* leg = new TLegend(0.475, 0.6, 0.895, 0.72, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(24); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hsignf2, "ITS2", "plm");
  leg->AddEntry(hsignf3, "ITS3", "plm");
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

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  info1.DrawLatex(0.155, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.155, 0.785, "Centrality 0^{ }#font[122]{-}20%");
  info1.DrawLatex(0.155, 0.73, "#it{L}_{int} = 10 nb^{-1}");
}

void drawSignificanceRatio(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_0505/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_0505/"){
  TString pathsignf2 = path2 + "signal_background_efficiency.root";
  TString pathsignf3 = path3 + "signal_background_efficiency.root";

  TFile* fsignf2 = new TFile(pathsignf2.Data());
  TFile* fsignf3 = new TFile(pathsignf3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsignf2->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgFONLL");

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsignf3->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgFONLL");

  TH1D* hsignf2 = new TH1D("hsignf2", ";#it{p}_{T} (GeV/#it{c});Expected significance ratio ITS3 / ITS2",nptbins,ptbinsfl);
  TH1D* hsignf3 = new TH1D("hsignf3", ";#it{p}_{T} (GeV/#it{c});Expected significance ratio ITS3 / ITS2",nptbins,ptbinsfl);

  for(int j = 0; j < nptbins; j++){
    cout << grSigPYTHIA2->GetY()[j] / TMath::Sqrt( grSigPYTHIA2->GetY()[j] + grBkgFONLL2->GetY()[j]) << endl;
    hsignf2->SetBinContent(j+1, grSigPYTHIA2->GetY()[j] / TMath::Sqrt( grSigPYTHIA2->GetY()[j] + grBkgFONLL2->GetY()[j]));
    hsignf2->SetBinError(j+1, 0.00001);
    cout << grSigPYTHIA3->GetY()[j] / TMath::Sqrt( grSigPYTHIA3->GetY()[j] + grBkgFONLL3->GetY()[j]) << endl;
    hsignf3->SetBinContent(j+1, grSigPYTHIA3->GetY()[j] / TMath::Sqrt( grSigPYTHIA3->GetY()[j] + grBkgFONLL3->GetY()[j]));
    hsignf3->SetBinError(j+1, 0.00001);
  }
  SetStyleHisto(hsignf2);
  SetStyleHisto(hsignf3);

  hsignf3->Divide(hsignf2);
  for(int j = 0; j < nptbins; j++) hsignf3->SetBinError(j+1, 0.00001);

  hsignf3->SetLineColor(colors[4]);
  hsignf3->SetMarkerColor(colors[4]);
  hsignf3->SetMarkerStyle(20);

  TCanvas* csignfratio = new TCanvas("csignfratio", "csignfratio", 450, 400);
  csignfratio->cd();
  hsignf3->Draw("ep");
  hsignf3->GetYaxis()->SetRangeUser(0.,3.);
  gPad->SetTickx();
  gPad->SetTicky();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  info1.DrawLatex(0.485, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.485, 0.785, "Centrality 0^{ }#font[122]{-}10%");
  info1.DrawLatex(0.485, 0.73, "#it{L}_{int} = 10 nb^{-1}");

  TLine* l = new TLine(0,1,24,1);
  l->SetLineStyle(2);
  l->Draw();
}

void drawRelStatError(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_0505/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_0505/"){
  TString pathsignf2 = path2 + "signal_background_efficiency.root";
  TString pathsignf3 = path3 + "signal_background_efficiency.root";

  TFile* fsignf2 = new TFile(pathsignf2.Data());
  TFile* fsignf3 = new TFile(pathsignf3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsignf2->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgFONLL");

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsignf3->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgFONLL");

  TH1D* hsignf2 = new TH1D("hsignf2", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbins2,ptbinsfl2);
  TH1D* hsignf3 = new TH1D("hsignf3", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbins2,ptbinsfl2);
  TH1D* hsyst = new TH1D("hsyst", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbins2,ptbinsfl2);

  for(int j = 0; j < nptbins2; j++){
    Double_t fs020 = 1;//2 * 18.83 / 23.26;
    Double_t fb020 = 1;//2 * 309.7 / 357.3;

    Double_t syst = TMath::Sqrt(yieldsyst[j]*yieldsyst[j] + selsyst[j]*selsyst[j] + pidsyst[j]*pidsyst[j] + ptshapesyst[j]*ptshapesyst[j] + trackingsyst[j]*trackingsyst[j]);

    if(j == 0 || j == 3){
      //hsignf2->SetBinContent(j+1, 0);
      //hsignf2->SetBinError(j+1, 0.00001);
    } else {
      hsignf2->SetBinContent(j+1, 1. / (fs020 * grSigPYTHIA2->GetY()[j] / TMath::Sqrt( fs020 * grSigPYTHIA2->GetY()[j] + fb020 * grBkgFONLL2->GetY()[j])));
      hsignf2->SetBinError(j+1, 0.00001);
    }
    hsignf3->SetBinContent(j+1, 1. / (fs020 * grSigPYTHIA3->GetY()[j] / TMath::Sqrt( fs020 * grSigPYTHIA3->GetY()[j] + fb020 * grBkgFONLL3->GetY()[j])));
    hsignf3->SetBinError(j+1, 0.00001);

    hsyst->SetBinContent(j+1, syst);
    hsyst->SetBinError(j+1, 0.00001);
  }
  SetStyleHisto(hsignf2);
  SetStyleHisto(hsignf3);
  SetStyleHisto(hsyst);

  hsignf2->SetLineColor(colors[0]);
  hsignf2->SetMarkerColor(colors[0]);
  hsignf2->SetMarkerStyle(20);
  hsignf2->SetLineWidth(3);

  hsignf3->SetLineColor(colors[1]);
  hsignf3->SetMarkerColor(colors[1]);
  hsignf3->SetMarkerStyle(21);
  hsignf3->SetLineWidth(3);

  hsyst->SetLineColor(kGray+3);
  hsyst->SetMarkerColor(kGray+3);
  hsyst->SetMarkerStyle(21);
  hsyst->SetLineWidth(3);
  hsyst->SetLineStyle(9);

  TCanvas* csignf = new TCanvas("csignf", "csignf", 450, 400);
  csignf->cd();
  hsignf2->Draw("hist");
  hsignf2->GetYaxis()->SetRangeUser(0.,0.55);
  gPad->SetTickx();
  gPad->SetTicky();
  hsignf3->Draw("same hist");
  hsyst->Draw("same hist");

  TLegend* leg = new TLegend(0.475, 0.57, 0.895, 0.75, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(22); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hsignf2, "stat. unc. ITS2", "lm");
  leg->AddEntry(hsignf3, "stat. unc. ITS3", "lm");
  leg->AddEntry(hsyst, "syst. unc.", "lm");
  leg->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  info1.DrawLatex(0.155, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.155, 0.785, "Centrality 0^{ }#font[122]{-}10%");
  info1.DrawLatex(0.155, 0.73, "#it{L}_{int} = 10 nb^{-1}");
}

void drawRAA(TString path2 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS2_0505/", TString path3 = "~/cernbox/Analyses/ML/input/DerivedResults/finalinput_ITS3_0505/"){
  TString pathsignf2 = path2 + "signal_background_efficiency.root";
  TString pathsignf3 = path3 + "signal_background_efficiency.root";

  TFile* fsignf2 = new TFile(pathsignf2.Data());
  TFile* fsignf3 = new TFile(pathsignf3.Data());

  TGraphAsymmErrors* grSigPYTHIA2 = (TGraphAsymmErrors*)fsignf2->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL2 = (TGraphAsymmErrors*)fsignf2->Get("grBkgFONLL");

  TGraphAsymmErrors* grSigPYTHIA3 = (TGraphAsymmErrors*)fsignf3->Get("grSigPythia");
  TGraphAsymmErrors* grBkgFONLL3 = (TGraphAsymmErrors*)fsignf3->Get("grBkgFONLL");

  TH1D* href = new TH1D("href", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",nptbins,ptbinsfl);
  TH1D* hraa2 = new TH1D("hraa2", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",nptbins2,ptbinsfl2);
  TH1D* hraa3 = new TH1D("hraa3", ";#it{p}_{T} (GeV/#it{c});Relative uncertainties",nptbins2,ptbinsfl2);
  TGraphErrors* hsyst = new TGraphErrors(nptbins2);

  for(int j = 0; j < nptbins2; j++){
    Double_t fs020 = 1;//2 * 18.83 / 23.26;
    Double_t fb020 = 1;//2 * 309.7 / 357.3;

    Double_t syst = TMath::Sqrt(yieldsyst[j]*yieldsyst[j] + selsyst[j]*selsyst[j] + pidsyst[j]*pidsyst[j] + ptshapesyst[j]*ptshapesyst[j] + trackingsyst[j]*trackingsyst[j]);

    if(j == 0 || j == 3){
      //hsignf2->SetBinContent(j+1, 0);
      //hsignf2->SetBinError(j+1, 0.00001);
    } else {
      hraa2->SetBinContent(j+1, centTAMU[j]);
      hraa2->SetBinError(j+1, centTAMU[j] / (fs020 * grSigPYTHIA2->GetY()[j] / TMath::Sqrt( fs020 * grSigPYTHIA2->GetY()[j] + fb020 * grBkgFONLL2->GetY()[j])));
    }
    hraa3->SetBinContent(j+1, centTAMU[j]);
    hraa3->SetBinError(j+1, centTAMU[j] / (fs020 * grSigPYTHIA3->GetY()[j] / TMath::Sqrt( fs020 * grSigPYTHIA3->GetY()[j] + fb020 * grBkgFONLL3->GetY()[j])));

    hsyst->SetPoint(j, ptbinsfl[j] + 0.5 * (ptbinsfl[j+1] - ptbinsfl[j]), centTAMU[j]);
    hsyst->SetPointError(j, 0.5, centTAMU[j] * syst);
  }
  SetStyleHisto(href);
  SetStyleHisto(hraa2);
  SetStyleHisto(hraa3);

  hraa2->SetLineColor(colors[0]);
  hraa2->SetMarkerColor(colors[0]);
  hraa2->SetMarkerStyle(20);
  hraa2->SetLineWidth(3);

  hraa3->SetLineColor(colors[1]);
  hraa3->SetMarkerColor(colors[1]);
  hraa3->SetMarkerStyle(21);
  hraa3->SetLineWidth(3);

  hsyst->SetLineColor(kGray+2);
  hsyst->SetMarkerColor(kGray+2);
  hsyst->SetMarkerStyle(21);
  hsyst->SetFillStyle(0); hsyst->SetFillColor(kWhite); hsyst->SetLineWidth(1); hsyst->SetLineColor(kGray+2);

  TCanvas* csignf = new TCanvas("csignf", "csignf", 450, 400);
  csignf->cd();
  href->Draw();
  href->GetYaxis()->SetRangeUser(0.,3);
  TGraph* grtamu = (TGraph*)extract_TAMU("theory/input_RAA_TAMU_Bs.txt");
  grtamu->Draw("same l");
  hraa2->Draw("same ep");
  gPad->SetTickx();
  gPad->SetTicky();
  hraa3->Draw("same ep");
  hsyst->Draw("same2");

  TLegend* leg = new TLegend(0.475, 0.55, 0.895, 0.65, 0, "NDC");
  leg->SetTextFont(43); leg->SetTextSize(20); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(hraa2, "ITS2", "plm");
  leg->AddEntry(hraa3, "ITS3", "plm");
  leg->Draw();

  TLegend* leg2 = new TLegend(0.475, 0.43, 0.895, 0.48, 0, "NDC");
  leg2->SetTextFont(43); leg2->SetTextSize(20); leg2->SetFillColor(0); leg2->SetFillStyle(0); leg2->SetLineColor(0);
  leg2->AddEntry(grtamu, "TAMU", "l");
  leg2->Draw();

  TLatex info1; info1.SetNDC(); info1.SetTextFont(43); info1.SetTextSize(18);
  info1.DrawLatex(0.35, 0.84, "Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.5 TeV");
  info1.DrawLatex(0.35, 0.785, "Centrality 0^{ }#font[122]{-}10%");
  info1.DrawLatex(0.35, 0.73, "#it{L}_{int} = 10 nb^{-1}");

  TLine* l = new TLine(0,1,24,1);
  l->SetLineStyle(2);
  l->Draw();
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