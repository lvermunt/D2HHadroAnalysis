Double_t fbtoB = 0.407; //http://pdg.lbl.gov/2019/reviews/rpp2018-rev-b-meson-prod-decay.pdf, table 85.1
Double_t fbtoBUnc = 0.007;
Double_t fLHCbBBs = 2*0.122; //https://journals.aps.org/prd/pdf/10.1103/PhysRevD.100.031102, factor 2 because B0 + B+
Double_t fLHCbBBsUnc = 2*0.006;
Double_t relSystFF = fbtoBUnc/fbtoB;
Double_t relSystLHCb = fLHCbBBsUnc/fLHCbBBs;

const Int_t nptbins = 5;
Float_t ptbins[nptbins+1] = {0, 4, 8, 12, 16, 24};

void extract_fonll(TString filnam, TString foutput){

  if(filnam=="") return 0x0;
  FILE* infil=fopen(filnam.Data(),"r");
  Char_t line[200];

  for(Int_t il=0; il<18; il++){
    fgets(line,200,infil);
    if(strstr(line,"central")) break;
  }
  Float_t ptmin,ptmax,csc,csmin,csmax,dum;
  fscanf(infil,"%f",&ptmin);

  TH1F* hfonll = new TH1F("hFONLLcent","Central FONLL (multiplied by F(b->B) and f(B/Bs));#it{p}_{T} (GeV/#it{c});",nptbins,ptbins);
  TH1F* hfonllmin = new TH1F("hFONLLmin","Minimum FONLL (multiplied by F(b->B) and f(B/Bs));#it{p}_{T} (GeV/#it{c});",nptbins,ptbins);
  TH1F* hfonllmax = new TH1F("hFONLLmax","Maximum FONLL (multiplied by F(b->B) and f(B/Bs));#it{p}_{T} (GeV/#it{c});",nptbins,ptbins);

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
    
    hfonllmin->SetBinContent(hfonllmin->FindBin(ptmin+0.1), csc - errtotdw);
    hfonll->SetBinContent(hfonll->FindBin(ptmin+0.1), csc);
    hfonllmax->SetBinContent(hfonllmax->FindBin(ptmin+0.1), csc + errtotup);

    ptmin=ptmax;
  }

  fclose(infil);

  TFile* fout = new TFile(foutput.Data(),"RECREATE");
  fout->cd();
  hfonll->Write();
  hfonllmin->Write();
  hfonllmax->Write();
  fout->Close();
}

void extract_TAMU(TString filn, TString foutput){

  TH1F* htamu = new TH1F("hTAMUcent","Central TAMU;",nptbins,ptbins);

  for(int j = 0; j < nptbins; j++){
    FILE* f=fopen(filn.Data(),"r");
    Float_t pt, raa;
    Double_t meanRAA = 0.;
    Int_t ncount = 0;
    while(!feof(f)){
      fscanf(f,"%f %f\n",&pt,&raa);

      if(pt > ptbins[j] && pt <= ptbins[j+1]){
        meanRAA += raa;
        ncount++;
      }
    }
    Double_t tamucent = meanRAA / ((double)ncount);
    fclose(f);

    htamu->SetBinContent(htamu->FindBin(ptbins[j]+0.1), tamucent);
  }
  TFile* fout = new TFile(foutput.Data(),"RECREATE");
  fout->cd();
  htamu->Write();
  fout->Close();

}

