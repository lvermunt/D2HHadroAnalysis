void SetStyleHisto(TH1D *h);

void drawInjectedBkg(TString filename, Int_t ptmin, Int_t ptmax, Int_t iscan){

  TGaxis::SetMaxDigits(2);

  TFile* f = new TFile(filename);

  TString suffix = Form("pt_cand%d_%d_%d", ptmin, ptmax, iscan);
  TF1* fdspr2 = (TF1*)f->Get(Form("fbkg_dspr%s", suffix.Data()));
  TF1* fdsfdbzero2 = (TF1*)f->Get(Form("fbkg_dsfdbzero%s", suffix.Data()));
  TF1* fdsfdbplus2 = (TF1*)f->Get(Form("fbkg_dsfdbplus%s", suffix.Data()));
  TF1* fdsfdlambdab2 = (TF1*)f->Get(Form("fbkg_dsfdlambdab%s", suffix.Data()));
  TF1* fdsfdbs2 = (TF1*)f->Get(Form("fbkg_dsfdbs%s", suffix.Data()));
  fdspr2->SetLineWidth(2);
  fdsfdbzero2->SetLineWidth(2);
  fdsfdbplus2->SetLineWidth(2);
  fdsfdlambdab2->SetLineWidth(2);
  fdsfdbs2->SetLineWidth(2);

  TF1* fdspr1 = (TF1*)f->Get(Form("fbkg_dspr1%s", suffix.Data()));
  TF1* fdsfdbzero1 = (TF1*)f->Get(Form("fbkg_dsfdbzero1%s", suffix.Data()));
  TF1* fdsfdbplus1 = (TF1*)f->Get(Form("fbkg_dsfdbplus1%s", suffix.Data()));
  TF1* fdsfdlambdab1 = (TF1*)f->Get(Form("fbkg_dsfdlambdab1%s", suffix.Data()));
  TF1* fdsfdbs1 = (TF1*)f->Get(Form("fbkg_dsfdbs1%s", suffix.Data()));
  fdspr1->SetLineWidth(2);
  fdsfdbzero1->SetLineWidth(2);
  fdsfdbplus1->SetLineWidth(2);
  fdsfdlambdab1->SetLineWidth(2);
  fdsfdbs1->SetLineWidth(2);

  TH1D *href = (TH1D*)f->Get(Form("h_norm_dspr%s", suffix.Data()));
  href->Reset("ICEMS");
  href->SetTitle(Form("%.1f < #it{p}_{T} < %.1f (GeV/#it{c})", (float)ptmin, (float)ptmax));
  href->SetMinimum(0.);
  href->SetMaximum(0.05);
  SetStyleHisto(href);

  TCanvas* cpol1 = new TCanvas(Form("cpol1_%d_%d_%d",ptmin,ptmax,iscan), "cpol1", 450, 400);
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
  leg->SetHeader("Injected D_{s}^{+} with HIJING #pi^{-}");
  leg->SetTextFont(43); leg->SetTextSize(16); leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineColor(0);
  leg->AddEntry(fdspr1, "(Prompt)   #scale[0.6]{ }D_{s} + #pi", "l");
  leg->AddEntry(fdsfdbzero1, "(B^{0} #rightarrow X +) D_{s} + #pi", "l");
  leg->AddEntry(fdsfdbplus1, "(B^{+} #rightarrow X +) D_{s} + #pi", "l");
  leg->AddEntry(fdsfdbs1, "(B_{s} #rightarrow X +) D_{s} + #pi", "l");
  leg->AddEntry(fdsfdlambdab1, "(#Lambda_{b} #rightarrow X +) D_{s} + #pi", "l");
  leg->Draw();
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
