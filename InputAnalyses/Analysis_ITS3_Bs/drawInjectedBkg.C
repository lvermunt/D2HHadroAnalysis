void SetStyleHisto(TH1D *h);

//https://indico.cern.ch/event/873138/contributions/3823930/attachments/2019548/3376384/D2H_FONLL_predictions_feeddown.pdf
Double_t ff_bplus = 0.09 * 0.412 / (0.117 * 0.412 + 0.09 * 0.412 + 0.93 * 0.088 + 0.011 * 0.089);
Double_t ff_bzero = 0.117 * 0.412 / (0.117 * 0.412 + 0.09 * 0.412 + 0.93 * 0.088 + 0.011 * 0.089);
Double_t ff_lambdab = 0.011 * 0.089 / (0.117 * 0.412 + 0.09 * 0.412 + 0.93 * 0.088 + 0.011 * 0.089);
Double_t ff_bs = 0.93 * 0.088 / (0.117 * 0.412 + 0.09 * 0.412 + 0.93 * 0.088 + 0.011 * 0.089);

void drawInjectedBkg(TString path, Int_t ptmin, Int_t ptmax, Int_t iscan){

  TGaxis::SetMaxDigits(2);

  TString filename = path + "/masshisto_bkgshape.root";
  TFile* f = new TFile(filename);

  TString filename2 = path + "/expected_Ds_signal_perevent.root";
  TFile* f2 = new TFile(filename2);
  TH1F* h_expsig_pr = (TH1F*)f2->Get("h_expsig_pr");
  TH1F* h_expsigmin_pr = (TH1F*)f2->Get("h_expsigmin_pr");
  TH1F* h_expsigmax_pr = (TH1F*)f2->Get("h_expsigmax_pr");
  TH1F* h_expsig_fd = (TH1F*)f2->Get("h_expsig_fd");
  TH1F* h_expsigmin_fd = (TH1F*)f2->Get("h_expsigmin_fd");
  TH1F* h_expsigmax_fd = (TH1F*)f2->Get("h_expsigmax_fd");
  Double_t exp_pr_cent_perev = h_expsig_pr->GetBinContent(h_expsig_pr->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_pr_min_perev = h_expsigmin_pr->GetBinContent(h_expsigmin_pr->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_pr_max_perev = h_expsigmax_pr->GetBinContent(h_expsigmax_pr->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_fd_cent_perev = h_expsig_fd->GetBinContent(h_expsig_fd->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_fd_min_perev = h_expsigmin_fd->GetBinContent(h_expsigmin_fd->FindBin(0.5 * (ptmin + ptmax)));
  Double_t exp_fd_max_perev = h_expsigmax_fd->GetBinContent(h_expsigmax_fd->FindBin(0.5 * (ptmin + ptmax)));

  TString suffixpt = Form("pt_cand%d_%d", ptmin, ptmax);
  TString suffix = Form("pt_cand%d_%d_%d", ptmin, ptmax, iscan);

  TH1F* hUniqueDs = (TH1F*)f->Get(Form("hNormUniqueDs%s",suffixpt.Data()));
  TH1F* hdspr = (TH1F*)f->Get(Form("hmass_DsPr%s", suffix.Data()));
  TH1F* hdsfdbplus = (TH1F*)f->Get(Form("hmass_DsFDBplus%s", suffix.Data()));
  TH1F* hdsfdbzero = (TH1F*)f->Get(Form("hmass_DsFDBzero%s", suffix.Data()));
  TH1F* hdsfdlambdab = (TH1F*)f->Get(Form("hmass_DsFDLambdab%s", suffix.Data()));
  TH1F* hdsfdbs = (TH1F*)f->Get(Form("hmass_DsFDBs%s", suffix.Data()));
  Double_t effBs_dspr = hdspr->GetEntries()/hUniqueDs->GetBinContent(1);
  Double_t effBs_dsfdbplus = hdsfdbplus->GetEntries()/hUniqueDs->GetBinContent(2) * ff_bplus;
  Double_t effBs_dsfdbzero = hdsfdbzero->GetEntries()/hUniqueDs->GetBinContent(3) * ff_bzero;
  Double_t effBs_dsfdlambdab = hdsfdlambdab->GetEntries()/hUniqueDs->GetBinContent(4) * ff_lambdab;
  Double_t effBs_dsfdbs = hdsfdbs->GetEntries()/hUniqueDs->GetBinContent(5) * ff_bs;
  cout << "Efficiency Ds -> Bs: " << effBs_dspr << " " << effBs_dsfdbplus << " " << effBs_dsfdbzero << " " << effBs_dsfdlambdab << " " << effBs_dsfdbs << endl;

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
  cout << "Expected (pr. Ds) + pi -> Bs cand: " << exp_Bspr_cent_perev << " + " << exp_Bspr_max_perev - exp_Bspr_cent_perev << " - " << exp_Bspr_cent_perev - exp_Bspr_min_perev << endl;
  cout << "Expected (B+ fd. Ds) + pi -> Bs cand: " << exp_Bsfdbplus_cent_perev << " + " << exp_Bsfdbplus_max_perev - exp_Bsfdbplus_cent_perev << " - " << exp_Bsfdbplus_cent_perev - exp_Bsfdbplus_min_perev << endl;
  cout << "Expected (B0 fd. Ds) + pi -> Bs cand: " << exp_Bsfdbzero_cent_perev << " + " << exp_Bsfdbzero_max_perev - exp_Bsfdbzero_cent_perev << " - " << exp_Bsfdbzero_cent_perev - exp_Bsfdbzero_min_perev << endl;
  cout << "Expected (Lb fd. Ds) + pi -> Bs cand: " << exp_Bsfdlambdab_cent_perev << " + " << exp_Bsfdlambdab_max_perev - exp_Bsfdlambdab_cent_perev << " - " << exp_Bsfdlambdab_cent_perev - exp_Bsfdlambdab_min_perev << endl;
  cout << "Expected (Bs fd. Ds) + pi -> Bs cand: " << exp_Bsfdbs_cent_perev << " + " << exp_Bsfdbs_max_perev - exp_Bsfdbs_cent_perev << " - " << exp_Bsfdbs_cent_perev - exp_Bsfdbs_min_perev << endl;

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
