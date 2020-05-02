void extract_TAMU_fd(TString filn, TString foutput){

  TH1F* htamu = new TH1F("hTAMUcent","Central TAMU;",125,0,25);

  FILE* f=fopen(filn.Data(),"r");
  Char_t line[200];
  for(Int_t il=0; il<5; il++){
    fgets(line,200,f);
    if(strstr(line,"GeV")) break;
  }

  Float_t pt, dummy, raa;
  while(!feof(f)){
    fscanf(f,"%f %f %f\n",&pt,&dummy,&raa);

    htamu->SetBinContent(htamu->FindBin(pt), raa);
  }
  fclose(f);

  TFile* fout = new TFile(foutput.Data(),"RECREATE");
  fout->cd();
  htamu->Write();
  fout->Close();
}

void extract_TAMU_pr(TString filn, TString foutput){

  TH1F* htamu = new TH1F("hTAMUcent","Central TAMU;",50,0,50);
  TH1F* htamumin = new TH1F("hTAMUmin","Minimum TAMU;",50,0,50);
  TH1F* htamumax = new TH1F("hTAMUmax","Maximum TAMU;",50,0,50);

  FILE* f=fopen(filn.Data(),"r");
  Char_t line[200];
  for(Int_t il=0; il<5; il++){
    fgets(line,200,f);
    if(strstr(line,"GeV")) break;
  }

  Float_t pt, raamin, raamax;
  while(!feof(f)){
    fscanf(f,"%f %f %f\n",&pt,&raamin,&raamax);

    htamumin->SetBinContent(htamu->FindBin(pt), raamin);
    htamumax->SetBinContent(htamu->FindBin(pt), raamax);
    htamu->SetBinContent(htamu->FindBin(pt), 0.5*(raamax + raamin));
  }
  fclose(f);

  TFile* fout = new TFile(foutput.Data(),"RECREATE");
  fout->cd();
  htamu->Write();
  htamumin->Write();
  htamumax->Write();
  fout->Close();
}
