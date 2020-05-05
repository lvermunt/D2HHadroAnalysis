const Int_t nptbins = 5;
Int_t ptlims[nptbins+1] = {0, 4, 8, 12, 16, 24};
Float_t ptlims_fl[nptbins+1] = {0, 4, 8, 12, 16, 24};

Double_t BRDs = 0.0227;
Double_t errBRDs = 0.0008;
Double_t TAA = 23.07 * 1e-9; //mb^-1 -> pb^-1 in which we put FONLL


void scaling_injDsHIJINGpion(TString path, TString outname){

  TFile* ftamu_pr = new TFile("theory/input_RAA_TAMU_prDs.root");
  TFile* ftamu_fd = new TFile("theory/input_RAA_TAMU_fdDs.root");
  TH1F* htamu_pr = (TH1F*)ftamu_pr->Get("hTAMUcent");
  TH1F* htamu_fd = (TH1F*)ftamu_fd->Get("hTAMUcent");

  TFile* ffonll = new TFile("fonll/DmesonLcPredictions_sqrt5500_y05_pythia8_FFee_BRPDG.root");
  TH1F* hfonll_pr = (TH1F*)ffonll->Get("hDsPhipitoKkpipred_central");
  TH1F* hfonll_fd = (TH1F*)ffonll->Get("hDsPhipitoKkpifromBpred_central_corr");
  TH1F* hfonllmin_pr = (TH1F*)ffonll->Get("hDsPhipitoKkpipred_min");
  TH1F* hfonllmin_fd = (TH1F*)ffonll->Get("hDsPhipitoKkpifromBpred_min_corr");
  TH1F* hfonllmax_pr = (TH1F*)ffonll->Get("hDsPhipitoKkpipred_max");
  TH1F* hfonllmax_fd = (TH1F*)ffonll->Get("hDsPhipitoKkpifromBpred_max_corr");

  TFile* feff = new TFile(Form("%s/efficiency_ITS3_Ds_050_50MeVbins.root",path.Data()));
  TH1F* h_gen_pr = (TH1F*)feff->Get("h_gen_pr");
  TH1F* h_gen_fd = (TH1F*)feff->Get("h_gen_fd");
  TH1F* h_presel_pr = (TH1F*)feff->Get("h_presel_pr");
  TH1F* h_presel_fd = (TH1F*)feff->Get("h_presel_fd");

  TH1F* h_match_pr[nptbins];
  TH1F* h_match_fd[nptbins];
  TFile* feff_match = new TFile(Form("%s/matchDspT_toBspt_ITS3corr_050_50MeVbins.root",path.Data()));
  for(int i = 0; i <  nptbins; i++){
    h_match_pr[i] = (TH1F*)feff_match->Get(Form("h_match_pr%d%d",ptlims[i], ptlims[i+1]));
    h_match_fd[i] = (TH1F*)feff_match->Get(Form("h_match_fd%d%d",ptlims[i], ptlims[i+1]));
  }
  Double_t test12 = 0;

  TH1F* hscaling_pr = new TH1F("hscaling_pr", "Prompt: FONLL-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{pr. D_{s}} / Ev", hfonll_pr->GetNbinsX()-1, hfonll_pr->GetBinLowEdge(1), hfonll_pr->GetBinLowEdge(hfonll_pr->GetNbinsX()));
  TH1F* hscalingmin_pr = new TH1F("hscalingmin_pr", "Prompt: FONLL_{min}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{pr. D_{s}} / Ev", hfonll_pr->GetNbinsX()-1, hfonll_pr->GetBinLowEdge(1), hfonll_pr->GetBinLowEdge(hfonll_pr->GetNbinsX()));
  TH1F* hscalingmax_pr = new TH1F("hscalingmax_pr", "Prompt: FONLL_{max}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{pr. D_{s}} / Ev", hfonll_pr->GetNbinsX()-1, hfonll_pr->GetBinLowEdge(1), hfonll_pr->GetBinLowEdge(hfonll_pr->GetNbinsX()));
  for(Int_t i = 1; i <= hfonll_pr->GetNbinsX(); i++){
    Double_t pt = hfonll_pr->GetBinCenter(i);
    Double_t dpt = hfonll_pr->GetBinWidth(i);
    Double_t yfonll = hfonll_pr->GetBinContent(i);
    Double_t yfonllmin = hfonllmin_pr->GetBinContent(i);
    Double_t yfonllmax = hfonllmax_pr->GetBinContent(i);
    Double_t ytamu = htamu_pr->GetBinContent(htamu_pr->FindBin(pt));
    Double_t dy = 1;
    Double_t acceffDs;
    if(h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt)) > 0) acceffDs = h_presel_pr->GetBinContent(h_presel_pr->FindBin(pt)) / h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt));
    else acceffDs = 0;

    Double_t signal = 2 * dpt * dy * BRDs * acceffDs * TAA * yfonll * ytamu;
    Double_t signalmin = 2 * dpt * dy * BRDs * acceffDs * TAA * yfonllmin * ytamu;
    Double_t signalmax = 2 * dpt * dy * BRDs * acceffDs * TAA * yfonllmax * ytamu;
    hscaling_pr->SetBinContent(i, signal);
    hscalingmin_pr->SetBinContent(i, signalmin);
    hscalingmax_pr->SetBinContent(i, signalmax);

    if(pt >= 1 && pt < 2){
      //test12 += signal*8000000000.;
      test12 += (0.002/0.00565959) * signal*8000000000.;
    }
  }

  TH1F* hscaling_fd = new TH1F("hscaling_fd", "Feed-down: FONLL-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{fd. D_{s}} / Ev", hfonll_fd->GetNbinsX()-1, hfonll_fd->GetBinLowEdge(1), hfonll_fd->GetBinLowEdge(hfonll_fd->GetNbinsX()));
  TH1F* hscalingmin_fd = new TH1F("hscalingmin_fd", "Feed-down: FONLL_{min}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{fd. D_{s}} / Ev", hfonll_fd->GetNbinsX()-1, hfonll_fd->GetBinLowEdge(1), hfonll_fd->GetBinLowEdge(hfonll_fd->GetNbinsX()));
  TH1F* hscalingmax_fd = new TH1F("hscalingmax_fd", "Feed-down: FONLL_{max}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{fd. D_{s}} / Ev", hfonll_fd->GetNbinsX()-1, hfonll_fd->GetBinLowEdge(1), hfonll_fd->GetBinLowEdge(hfonll_fd->GetNbinsX()));
  for(Int_t i = 1; i <= hfonll_fd->GetNbinsX(); i++){
    Double_t pt = hfonll_fd->GetBinCenter(i);
    Double_t dpt = hfonll_fd->GetBinWidth(i);
    Double_t yfonll = hfonll_fd->GetBinContent(i);
    Double_t yfonllmin = hfonllmin_fd->GetBinContent(i);
    Double_t yfonllmax = hfonllmax_fd->GetBinContent(i);
    Double_t ytamu = htamu_fd->GetBinContent(htamu_fd->FindBin(pt));
    Double_t dy = 1;
    Double_t acceffDs;
    if(h_gen_fd->GetBinContent(h_gen_fd->FindBin(pt)) > 0) acceffDs = h_presel_fd->GetBinContent(h_presel_fd->FindBin(pt)) / h_gen_fd->GetBinContent(h_gen_fd->FindBin(pt));
    else acceffDs = 0;

    Double_t signal = 2 * dpt * dy * BRDs * acceffDs * TAA * yfonll * ytamu;
    Double_t signalmin = 2 * dpt * dy * BRDs * acceffDs * TAA * yfonllmin * ytamu;
    Double_t signalmax = 2 * dpt * dy * BRDs * acceffDs * TAA * yfonllmax * ytamu;
    hscaling_fd->SetBinContent(i, signal);
    hscalingmin_fd->SetBinContent(i, signalmin);
    hscalingmax_fd->SetBinContent(i, signalmax);

    if(pt >= 1 && pt < 2){
      //test12 += signal*8000000000.;
      test12 += (0.006/0.00523262) * signal*8000000000.;
    }
  }
  //cout << test12  << endl;

  TH1F* h_expsig_pr = new TH1F("h_expsig_pr", ";#it{p}_{T}(B_{s}) (GeV/it{c});Expected pr. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsigmin_pr = new TH1F("h_expsigmin_pr", ";#it{p}_{T}(B_{s}) (GeV/it{c});Min. Expected pr. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsigmax_pr = new TH1F("h_expsigmax_pr", ";#it{p}_{T}(B_{s}) (GeV/it{c});Max. Expected pr. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsig_fd = new TH1F("h_expsig_fd", ";#it{p}_{T}(B_{s}) (GeV/it{c});Expected fd. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsigmin_fd = new TH1F("h_expsigmin_fd", ";#it{p}_{T}(B_{s}) (GeV/it{c});Min. Expected fd. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsigmax_fd = new TH1F("h_expsigmax_fd", ";#it{p}_{T}(B_{s}) (GeV/it{c});Max. Expected fd. signal per event", nptbins, ptlims_fl);
  for(int i = 0; i < nptbins; i++){
    Double_t pt = 0.5 * (ptlims[i] + ptlims[i+1]);

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

      if(hscaling_pr->GetBinCenter(j) >= ptlims[i] && hscaling_pr->GetBinCenter(j) < ptlims[i+1]){
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

      if(hscaling_fd->GetBinCenter(j) >= ptlims[i] && hscaling_fd->GetBinCenter(j) < ptlims[i+1]){
        sigfd_nonmatch += hscaling_fd->GetBinContent(j);
        sigfdmin_nonmatch += hscalingmin_fd->GetBinContent(j);
        sigfdmax_nonmatch += hscalingmax_fd->GetBinContent(j);
      }
    }
    h_expsig_fd->SetBinContent(h_expsig_fd->FindBin(pt), sigfd);
    h_expsigmin_fd->SetBinContent(h_expsigmin_fd->FindBin(pt), sigfdmin);
    h_expsigmax_fd->SetBinContent(h_expsigmax_fd->FindBin(pt), sigfdmax);

    //cout << "   " << 8000000000*sigpr << " " << 8000000000*sigfd << " " << 8000000000*sigpr + 8000000000*sigfd << endl;
    //cout << "      " << 8000000000*sigpr_nonmatch << " " << 8000000000*sigfd_nonmatch << " " << 8000000000*sigpr_nonmatch + 8000000000*sigfd_nonmatch << endl;
    //cout << "         " << (8000000000*sigpr + 8000000000*sigfd) / (8000000000*sigpr_nonmatch + 8000000000*sigfd_nonmatch) << endl;

  }

  TCanvas* c1 = new TCanvas();
  c1->cd();
  hscaling_pr->Draw();
  hscalingmin_pr->SetLineStyle(2);
  hscalingmin_pr->Draw("same");
  hscalingmax_pr->SetLineStyle(2);
  hscalingmax_pr->Draw("same");
  hscaling_fd->SetLineColor(kRed);
  hscalingmin_fd->SetLineColor(kRed);
  hscalingmax_fd->SetLineColor(kRed);
  hscalingmin_fd->SetLineStyle(2);
  hscalingmax_fd->SetLineStyle(2);
  hscaling_fd->Draw("same");
  hscalingmin_fd->Draw("same");
  hscalingmax_fd->Draw("same");

  TCanvas* c2 = new TCanvas();
  c2->cd();
  h_expsig_pr->Draw();
  h_expsigmin_pr->SetLineStyle(2);
  h_expsigmin_pr->Draw("same");
  h_expsigmax_pr->SetLineStyle(2);
  h_expsigmax_pr->Draw("same");
  h_expsig_fd->SetLineColor(kRed);
  h_expsigmin_fd->SetLineColor(kRed);
  h_expsigmax_fd->SetLineColor(kRed);
  h_expsigmin_fd->SetLineStyle(2);
  h_expsigmax_fd->SetLineStyle(2);
  h_expsig_fd->Draw("same");
  h_expsigmin_fd->Draw("same");
  h_expsigmax_fd->Draw("same");

  outname = path + outname;
  TFile* fout = new TFile(outname.Data(),"RECREATE");
  fout->cd();
  h_expsig_pr->Write();
  h_expsigmin_pr->Write();
  h_expsigmax_pr->Write();
  h_expsig_fd->Write();
  h_expsigmin_fd->Write();
  h_expsigmax_fd->Write();
  fout->Close();

  //cout << test << " " << testeffpr/(double)nbinssum << " " << testefffd/(double)nbinssum << endl;
  //13311 0.00565959 0.00523262 (for 1 < pt < 2) (scaling by efficiency difference with Fabrizio)
  //vs
  //((70*70)*0.9 + 70*70)/0.9  (from looking at plots https://indico.cern.ch/event/893570/contributions/3769353/attachments/1999626/3337235/Ds_ITS3_estimates_2020_03_06.pdf#page=7)
  //(double) 10344.444
  //So close, think it is ok.

  //Removed code above, see below snippets
  //  Double_t acceffDs;
  //  if(h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt)) > 0) acceffDs = h_presel_pr->GetBinContent(h_presel_pr->FindBin(pt)) / h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt));
  //  else acceffDs = 0;

  // Double_t signal = 2 * dpt * dy * BRDs * acceffDs * TAA * yfonll * ytamu;

  //  if(pt >= 1 && pt < 2){
  //    testeffpr += acceffDs;
  //    test += (0.002/0.00565959) * signal*8000000000.;
  //    nbinssum++;
  //  }
  //

}


/*


void scaling_injDsHIJINGpion(TString path){

  TFile* ftamu_pr = new TFile("theory/input_RAA_TAMU_prDs.root");
  TFile* ftamu_fd = new TFile("theory/input_RAA_TAMU_fdDs.root");
  TH1F* htamu_pr = (TH1F*)ftamu_pr->Get("hTAMUcent");
  TH1F* htamu_fd = (TH1F*)ftamu_fd->Get("hTAMUcent");

  TFile* ffonll = new TFile("fonll/DmesonLcPredictions_sqrt5500_y05_pythia8_FFee_BRPDG.root");
  TH1F* hfonll_pr = (TH1F*)ffonll->Get("hDsPhipitoKkpipred_central");
  TH1F* hfonll_fd = (TH1F*)ffonll->Get("hDsPhipitoKkpifromBpred_central_corr");
  TH1F* hfonllmin_pr = (TH1F*)ffonll->Get("hDsPhipitoKkpipred_min");
  TH1F* hfonllmin_fd = (TH1F*)ffonll->Get("hDsPhipitoKkpifromBpred_min_corr");
  TH1F* hfonllmax_pr = (TH1F*)ffonll->Get("hDsPhipitoKkpipred_max");
  TH1F* hfonllmax_fd = (TH1F*)ffonll->Get("hDsPhipitoKkpifromBpred_max_corr");

  TFile* feff = new TFile(Form("%s/efficiency_ITS2_Ds_050_50MeVbins.root",path.Data()));
  TH1F* h_gen_pr = (TH1F*)feff->Get("h_gen_pr");
  TH1F* h_gen_fd = (TH1F*)feff->Get("h_gen_fd");
  TH1F* h_presel_pr = (TH1F*)feff->Get("h_presel_pr");
  TH1F* h_presel_fd = (TH1F*)feff->Get("h_presel_fd");

  TH1F* h_match_pr[nptbins];
  TH1F* h_match_fd[nptbins];
  TFile* feff_match = new TFile(Form("%s/matchDspT_toBspt_ITS2_050_50MeVbins.root",path.Data()));
  for(int i = 0; i <  nptbins; i++){
    h_match_pr[i] = (TH1F*)feff_match->Get(Form("h_match_pr%d%d",ptlims[i], ptlims[i+1]));
    h_match_fd[i] = (TH1F*)feff_match->Get(Form("h_match_fd%d%d",ptlims[i], ptlims[i+1]));
  }

  TH1F* hscaling_pr = new TH1F("hscaling_pr", "Prompt: FONLL-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{pr. D_{s}} / Ev", hfonll_pr->GetNbinsX()-1, hfonll_pr->GetBinLowEdge(1), hfonll_pr->GetBinLowEdge(hfonll_pr->GetNbinsX()));
  TH1F* hscalingmin_pr = new TH1F("hscalingmin_pr", "Prompt: FONLL_{min}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{pr. D_{s}} / Ev", hfonll_pr->GetNbinsX()-1, hfonll_pr->GetBinLowEdge(1), hfonll_pr->GetBinLowEdge(hfonll_pr->GetNbinsX()));
  TH1F* hscalingmax_pr = new TH1F("hscalingmax_pr", "Prompt: FONLL_{max}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{pr. D_{s}} / Ev", hfonll_pr->GetNbinsX()-1, hfonll_pr->GetBinLowEdge(1), hfonll_pr->GetBinLowEdge(hfonll_pr->GetNbinsX()));
  for(Int_t i = 1; i <= hfonll_pr->GetNbinsX(); i++){
    Double_t pt = hfonll_pr->GetBinCenter(i);
    Double_t dpt = hfonll_pr->GetBinWidth(i);
    Double_t yfonll = hfonll_pr->GetBinContent(i);
    Double_t yfonllmin = hfonllmin_pr->GetBinContent(i);
    Double_t yfonllmax = hfonllmax_pr->GetBinContent(i);
    Double_t ytamu = htamu_pr->GetBinContent(htamu_pr->FindBin(pt));
    Double_t dy = 1;

    Double_t signal = 2 * dpt * dy * BRDs * TAA * yfonll * ytamu;
    Double_t signalmin = 2 * dpt * dy * BRDs * TAA * yfonllmin * ytamu;
    Double_t signalmax = 2 * dpt * dy * BRDs * TAA * yfonllmax * ytamu;
    hscaling_pr->SetBinContent(i, signal);
    hscalingmin_pr->SetBinContent(i, signalmin);
    hscalingmax_pr->SetBinContent(i, signalmax);
  }

  TH1F* hscaling_fd = new TH1F("hscaling_fd", "Feed-down: FONLL-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{fd. D_{s}} / Ev", hfonll_fd->GetNbinsX()-1, hfonll_fd->GetBinLowEdge(1), hfonll_fd->GetBinLowEdge(hfonll_fd->GetNbinsX()));
  TH1F* hscalingmin_fd = new TH1F("hscalingmin_fd", "Feed-down: FONLL_{min}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{fd. D_{s}} / Ev", hfonll_fd->GetNbinsX()-1, hfonll_fd->GetBinLowEdge(1), hfonll_fd->GetBinLowEdge(hfonll_fd->GetNbinsX()));
  TH1F* hscalingmax_fd = new TH1F("hscalingmax_fd", "Feed-down: FONLL_{max}-TAMU Ds scaling per event (without Bs eff);#it{p}_{T} (GeV/#it{c});Expected #it{N}_{fd. D_{s}} / Ev", hfonll_fd->GetNbinsX()-1, hfonll_fd->GetBinLowEdge(1), hfonll_fd->GetBinLowEdge(hfonll_fd->GetNbinsX()));
  for(Int_t i = 1; i <= hfonll_fd->GetNbinsX(); i++){
    Double_t pt = hfonll_fd->GetBinCenter(i);
    Double_t dpt = hfonll_fd->GetBinWidth(i);
    Double_t yfonll = hfonll_fd->GetBinContent(i);
    Double_t yfonllmin = hfonllmin_fd->GetBinContent(i);
    Double_t yfonllmax = hfonllmax_fd->GetBinContent(i);
    Double_t ytamu = htamu_fd->GetBinContent(htamu_fd->FindBin(pt));
    Double_t dy = 1;

    Double_t signal = 2 * dpt * dy * BRDs * TAA * yfonll * ytamu;
    Double_t signalmin = 2 * dpt * dy * BRDs * TAA * yfonllmin * ytamu;
    Double_t signalmax = 2 * dpt * dy * BRDs * TAA * yfonllmax * ytamu;
    hscaling_fd->SetBinContent(i, signal);
    hscalingmin_fd->SetBinContent(i, signalmin);
    hscalingmax_fd->SetBinContent(i, signalmax);
  }

  TH1F* h_expsig_pr = new TH1F("h_expsig_pr", ";#it{p}_{T}(B_{s}) (GeV/it{c});Expected pr. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsigmin_pr = new TH1F("h_expsigmin_pr", ";#it{p}_{T}(B_{s}) (GeV/it{c});Min. Expected pr. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsigmax_pr = new TH1F("h_expsigmax_pr", ";#it{p}_{T}(B_{s}) (GeV/it{c});Max. Expected pr. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsig_fd = new TH1F("h_expsig_fd", ";#it{p}_{T}(B_{s}) (GeV/it{c});Expected fd. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsigmin_fd = new TH1F("h_expsigmin_fd", ";#it{p}_{T}(B_{s}) (GeV/it{c});Min. Expected fd. signal per event", nptbins, ptlims_fl);
  TH1F* h_expsigmax_fd = new TH1F("h_expsigmax_fd", ";#it{p}_{T}(B_{s}) (GeV/it{c});Max. Expected fd. signal per event", nptbins, ptlims_fl);
  for(int i = 0; i < nptbins; i++){
    Double_t pt = 0.5 * (ptlims[i] + ptlims[i+1]);

    Double_t acceffDs_pr_match = 0;
    Double_t acceffDs_pr_nonmatch = 0;
    Int_t nbins_nonmatch = 0;
    for(int j = 1; j <= h_gen_pr->GetNbinsX(); j++){
      if(h_gen_pr->GetBinContent(j) > 0){
        acceffDs_pr_match += h_match_pr[i]->GetBinContent(j) * h_presel_pr->GetBinContent(j) / h_gen_pr->GetBinContent(j);
        if(hscaling_pr->GetBinCenter(j) >= ptlims[i] && hscaling_pr->GetBinCenter(j) < ptlims[i+1]){
          acceffDs_pr_nonmatch += h_presel_pr->GetBinContent(j) / h_gen_pr->GetBinContent(j);
          nbins_nonmatch++;
        }
      }
    }
    acceffDs_pr_nonmatch /= ((double)nbins_nonmatch);

    Double_t sigpr = 0;
    Double_t sigprmin = 0;
    Double_t sigprmax = 0;
    Double_t sigpr_nonmatch = 0;
    Double_t sigprmin_nonmatch = 0;
    Double_t sigprmax_nonmatch = 0;
    for(int j = 1; j <= h_gen_pr->GetNbinsX(); j++){
      if(hscaling_pr->GetBinCenter(j) >= ptlims[i] && hscaling_pr->GetBinCenter(j) < ptlims[i+1]){
        sigpr += acceffDs_pr_match * hscaling_pr->GetBinContent(j);
        sigprmin += acceffDs_pr_match * hscalingmin_pr->GetBinContent(j);
        sigprmax += acceffDs_pr_match * hscalingmax_pr->GetBinContent(j);

        sigpr_nonmatch += acceffDs_pr_nonmatch * hscaling_pr->GetBinContent(j);
        sigprmin_nonmatch += acceffDs_pr_nonmatch * hscalingmin_pr->GetBinContent(j);
        sigprmax_nonmatch += acceffDs_pr_nonmatch * hscalingmax_pr->GetBinContent(j);
      }
    }
    h_expsig_pr->SetBinContent(h_expsig_pr->FindBin(pt), 8000000000*sigpr);
    h_expsigmin_pr->SetBinContent(h_expsigmin_pr->FindBin(pt), 8000000000*sigprmin);
    h_expsigmax_pr->SetBinContent(h_expsigmax_pr->FindBin(pt), 8000000000*sigprmax);


    Double_t acceffDs_fd_match = 0;
    Double_t acceffDs_fd_nonmatch = 0;
    nbins_nonmatch = 0;
    for(int j = 1; j <= h_gen_fd->GetNbinsX(); j++){
      if(h_gen_fd->GetBinContent(j) > 0){
        acceffDs_fd_match += h_match_fd[i]->GetBinContent(j) * h_presel_fd->GetBinContent(j) / h_gen_fd->GetBinContent(j);
        if(hscaling_fd->GetBinCenter(j) >= ptlims[i] && hscaling_fd->GetBinCenter(j) < ptlims[i+1]){
          acceffDs_fd_nonmatch += h_presel_fd->GetBinContent(j) / h_gen_fd->GetBinContent(j);
          nbins_nonmatch++;
        }
      }
    }
    acceffDs_fd_nonmatch /= ((double)nbins_nonmatch);

    Double_t sigfd = 0;
    Double_t sigfdmin = 0;
    Double_t sigfdmax = 0;
    Double_t sigfd_nonmatch = 0;
    Double_t sigfdmin_nonmatch = 0;
    Double_t sigfdmax_nonmatch = 0;
    for(int j = 1; j <= h_gen_fd->GetNbinsX(); j++){
      if(hscaling_fd->GetBinCenter(j) >= ptlims[i] && hscaling_fd->GetBinCenter(j) < ptlims[i+1]){
        sigfd += acceffDs_fd_match * hscaling_fd->GetBinContent(j);
        sigfdmin += acceffDs_fd_match * hscalingmin_fd->GetBinContent(j);
        sigfdmax += acceffDs_fd_match * hscalingmax_fd->GetBinContent(j);

        sigfd_nonmatch += acceffDs_fd_nonmatch * hscaling_fd->GetBinContent(j);
        sigfdmin_nonmatch += acceffDs_fd_nonmatch * hscalingmin_fd->GetBinContent(j);
        sigfdmax_nonmatch += acceffDs_fd_nonmatch * hscalingmax_fd->GetBinContent(j);
      }
    }
    h_expsig_fd->SetBinContent(h_expsig_fd->FindBin(pt), 8000000000*sigfd);
    h_expsigmin_fd->SetBinContent(h_expsigmin_fd->FindBin(pt), 8000000000*sigfdmin);
    h_expsigmax_fd->SetBinContent(h_expsigmax_fd->FindBin(pt), 8000000000*sigfdmax);

    cout << "   " << 8000000000*sigpr << " " << 8000000000*sigfd << " " << 8000000000*sigpr + 8000000000*sigfd << endl;
    cout << "      " << 8000000000*sigpr_nonmatch << " " << 8000000000*sigfd_nonmatch << " " << 8000000000*sigpr_nonmatch + 8000000000*sigfd_nonmatch << endl;
    cout << "         " << (8000000000*sigpr + 8000000000*sigfd) / (8000000000*sigpr_nonmatch + 8000000000*sigfd_nonmatch) << endl;
  }

  TCanvas* c1 = new TCanvas();
  c1->cd();
  hscaling_pr->Draw();
  hscalingmin_pr->SetLineStyle(2);
  hscalingmin_pr->Draw("same");
  hscalingmax_pr->SetLineStyle(2);
  hscalingmax_pr->Draw("same");
  hscaling_fd->SetLineColor(kRed);
  hscalingmin_fd->SetLineColor(kRed);
  hscalingmax_fd->SetLineColor(kRed);
  hscalingmin_fd->SetLineStyle(2);
  hscalingmax_fd->SetLineStyle(2);
  hscaling_fd->Draw("same");
  hscalingmin_fd->Draw("same");
  hscalingmax_fd->Draw("same");

  TCanvas* c2 = new TCanvas();
  c2->cd();
  h_expsig_pr->Draw();
  h_expsigmin_pr->SetLineStyle(2);
  h_expsigmin_pr->Draw("same");
  h_expsigmax_pr->SetLineStyle(2);
  h_expsigmax_pr->Draw("same");
  h_expsig_fd->SetLineColor(kRed);
  h_expsigmin_fd->SetLineColor(kRed);
  h_expsigmax_fd->SetLineColor(kRed);
  h_expsigmin_fd->SetLineStyle(2);
  h_expsigmax_fd->SetLineStyle(2);
  h_expsig_fd->Draw("same");
  h_expsigmin_fd->Draw("same");
  h_expsigmax_fd->Draw("same");

  //cout << test << " " << testeffpr/(double)nbinssum << " " << testefffd/(double)nbinssum << endl;
  //13311 0.00565959 0.00523262 (for 1 < pt < 2) (scaling by efficiency difference with Fabrizio)
  //vs
  //((70*70)*0.9 + 70*70)/0.9  (from looking at plots https://indico.cern.ch/event/893570/contributions/3769353/attachments/1999626/3337235/Ds_ITS3_estimates_2020_03_06.pdf#page=7)
  //(double) 10344.444
  //So close, think it is ok.

  //Removed code above, see below snippets
  //  Double_t acceffDs;
  //  if(h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt)) > 0) acceffDs = h_presel_pr->GetBinContent(h_presel_pr->FindBin(pt)) / h_gen_pr->GetBinContent(h_gen_pr->FindBin(pt));
  //  else acceffDs = 0;

  // Double_t signal = 2 * dpt * dy * BRDs * acceffDs * TAA * yfonll * ytamu;

  //  if(pt >= 1 && pt < 2){
  //    testeffpr += acceffDs;
  //    test += (0.002/0.00565959) * signal*8000000000.;
  //    nbinssum++;
  //  }
  //

}
*/