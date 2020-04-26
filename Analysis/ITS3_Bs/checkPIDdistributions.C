void checkPIDdistributions(TString filename, Bool_t its2 = kTRUE){
	
  TFile* f = new TFile(filename.Data());

  TDirectoryFile* dir;
  if(its2) dir = (TDirectoryFile*)f->Get("PWGHF_TreeCreatorITS2");
  else dir = (TDirectoryFile*)f->Get("PWGHF_TreeCreatorITS3");

  TTree* tr = (TTree*)dir->Get("tree_Bs");

  float nsigTPC_Pi_0, nsigTPC_K_0, nsigTOF_Pi_0, nsigTOF_K_0;
  float nsigTPC_Pi_1, nsigTPC_K_1, nsigTOF_Pi_1, nsigTOF_K_1;
  float nsigTPC_Pi_2, nsigTPC_K_2, nsigTOF_Pi_2, nsigTOF_K_2;
  float nsigTPC_Pi_3, nsigTPC_K_3, nsigTOF_Pi_3, nsigTOF_K_3;

  float pt_prong0, pt_prong1, pt_prong2, pt_prong3;

  tr->SetBranchAddress("nsigTPC_Pi_0",&nsigTPC_Pi_0);
  tr->SetBranchAddress("nsigTPC_K_0",&nsigTPC_K_0);
  tr->SetBranchAddress("nsigTOF_Pi_0",&nsigTOF_Pi_0);
  tr->SetBranchAddress("nsigTOF_K_0",&nsigTOF_K_0);

  tr->SetBranchAddress("nsigTPC_Pi_1",&nsigTPC_Pi_1);
  tr->SetBranchAddress("nsigTPC_K_1",&nsigTPC_K_1);
  tr->SetBranchAddress("nsigTOF_Pi_1",&nsigTOF_Pi_1);
  tr->SetBranchAddress("nsigTOF_K_1",&nsigTOF_K_1);

  tr->SetBranchAddress("nsigTPC_Pi_2",&nsigTPC_Pi_2);
  tr->SetBranchAddress("nsigTPC_K_2",&nsigTPC_K_2);
  tr->SetBranchAddress("nsigTOF_Pi_2",&nsigTOF_Pi_2);
  tr->SetBranchAddress("nsigTOF_K_2",&nsigTOF_K_2);

  tr->SetBranchAddress("nsigTPC_Pi_3",&nsigTPC_Pi_3);
  tr->SetBranchAddress("nsigTPC_K_3",&nsigTPC_K_3);
  tr->SetBranchAddress("nsigTOF_Pi_3",&nsigTOF_Pi_3);
  tr->SetBranchAddress("nsigTOF_K_3",&nsigTOF_K_3);

  tr->SetBranchAddress("pt_prong0",&pt_prong0);
  tr->SetBranchAddress("pt_prong1",&pt_prong1);
  tr->SetBranchAddress("pt_prong2",&pt_prong2);
  tr->SetBranchAddress("pt_prong3",&pt_prong3);

  TH1F* hnsigTPC_Pi_0 = new TH1F("hnsigTPC_Pi_0","hnsigTPC_Pi_0",500,-5,5);
  TH1F* hnsigTPC_K_0 = new TH1F("hnsigTPC_K_0","hnsigTPC_Pi_0",500,-5,5);
  TH1F* hnsigTOF_Pi_0 = new TH1F("hnsigTOF_Pi_0","hnsigTOF_Pi_0",500,-5,5);
  TH1F* hnsigTOF_K_0 = new TH1F("hnsigTOF_K_0","hnsigTOF_K_0",500,-5,5);
  TH1F* hnsigTPC_Pi_1 = new TH1F("hnsigTPC_Pi_1","hnsigTPC_Pi_1",500,-5,5);
  TH1F* hnsigTPC_K_1 = new TH1F("hnsigTPC_K_1","hnsigTPC_K_1",500,-5,5);
  TH1F* hnsigTOF_Pi_1 = new TH1F("hnsigTOF_Pi_1","hnsigTOF_Pi_1",500,-5,5);
  TH1F* hnsigTOF_K_1 = new TH1F("hnsigTOF_K_1","hnsigTOF_K_1",500,-5,5);
  TH1F* hnsigTPC_Pi_2 = new TH1F("hnsigTPC_Pi_2","hnsigTPC_Pi_2",500,-5,5);
  TH1F* hnsigTPC_K_2 = new TH1F("hnsigTPC_K_2","hnsigTPC_K_2",500,-5,5);
  TH1F* hnsigTOF_Pi_2 = new TH1F("hnsigTOF_Pi_2","hnsigTOF_Pi_2",500,-5,5);
  TH1F* hnsigTOF_K_2 = new TH1F("hnsigTOF_K_2","hnsigTOF_K_2",500,-5,5);
  TH1F* hnsigTPC_Pi_3 = new TH1F("hnsigTPC_Pi_3","hnsigTPC_Pi_3",500,-5,5);
  TH1F* hnsigTPC_K_3 = new TH1F("hnsigTPC_K_3","hnsigTPC_K_3",500,-5,5);
  TH1F* hnsigTOF_Pi_3 = new TH1F("hnsigTOF_Pi_3","hnsigTOF_Pi_3",500,-5,5);
  TH1F* hnsigTOF_K_3 = new TH1F("hnsigTOF_K_3","hnsigTOF_K_3",500,-5,5);

  TH2F* h2nsigTPC_Pi_0 = new TH2F("h2nsigTPC_Pi_0","hnsigTPC_Pi_0",50,0,5,100,-5,5);
  TH2F* h2nsigTPC_K_0 = new TH2F("h2nsigTPC_K_0","hnsigTPC_Pi_0",50,0,5,100,-5,5);
  TH2F* h2nsigTOF_Pi_0 = new TH2F("h2nsigTOF_Pi_0","hnsigTOF_Pi_0",50,0,5,100,-5,5);
  TH2F* h2nsigTOF_K_0 = new TH2F("h2nsigTOF_K_0","hnsigTOF_K_0",50,0,5,100,-5,5);
  TH2F* h2nsigTPC_Pi_1 = new TH2F("h2nsigTPC_Pi_1","hnsigTPC_Pi_1",50,0,5,100,-5,5);
  TH2F* h2nsigTPC_K_1 = new TH2F("h2nsigTPC_K_1","hnsigTPC_K_1",50,0,5,100,-5,5);
  TH2F* h2nsigTOF_Pi_1 = new TH2F("h2nsigTOF_Pi_1","hnsigTOF_Pi_1",50,0,5,100,-5,5);
  TH2F* h2nsigTOF_K_1 = new TH2F("h2nsigTOF_K_1","hnsigTOF_K_1",50,0,5,100,-5,5);
  TH2F* h2nsigTPC_Pi_2 = new TH2F("h2nsigTPC_Pi_2","hnsigTPC_Pi_2",50,0,5,100,-5,5);
  TH2F* h2nsigTPC_K_2 = new TH2F("h2nsigTPC_K_2","hnsigTPC_K_2",50,0,5,100,-5,5);
  TH2F* h2nsigTOF_Pi_2 = new TH2F("h2nsigTOF_Pi_2","hnsigTOF_Pi_2",50,0,5,100,-5,5);
  TH2F* h2nsigTOF_K_2 = new TH2F("h2nsigTOF_K_2","hnsigTOF_K_2",50,0,5,100,-5,5);
  TH2F* h2nsigTPC_Pi_3 = new TH2F("h2nsigTPC_Pi_3","hnsigTPC_Pi_3",50,0,5,100,-5,5);
  TH2F* h2nsigTPC_K_3 = new TH2F("h2nsigTPC_K_3","hnsigTPC_K_3",50,0,5,100,-5,5);
  TH2F* h2nsigTOF_Pi_3 = new TH2F("h2nsigTOF_Pi_3","hnsigTOF_Pi_3",50,0,5,100,-5,5);
  TH2F* h2nsigTOF_K_3 = new TH2F("h2nsigTOF_K_3","hnsigTOF_K_3",50,0,5,100,-5,5);

  Int_t entr = tr->GetEntriesFast();
  for(int i = 0; i < entr; i++){
    //if(i > 1000000) continue;
    if(i % 100000 == 0) cout << "Analysing entry: " << i << " of " << entr << endl;
    tr->GetEntry(i);
    hnsigTPC_Pi_0->Fill(nsigTPC_Pi_0);
    hnsigTPC_K_0->Fill(nsigTPC_K_0);
    hnsigTOF_Pi_0->Fill(nsigTOF_Pi_0);
    hnsigTOF_K_0->Fill(nsigTOF_K_0);

    hnsigTPC_Pi_1->Fill(nsigTPC_Pi_1);
    hnsigTPC_K_1->Fill(nsigTPC_K_1);
    hnsigTOF_Pi_1->Fill(nsigTOF_Pi_1);
    hnsigTOF_K_1->Fill(nsigTOF_K_1);

    hnsigTPC_Pi_2->Fill(nsigTPC_Pi_2);
    hnsigTPC_K_2->Fill(nsigTPC_K_2);
    hnsigTOF_Pi_2->Fill(nsigTOF_Pi_2);
    hnsigTOF_K_2->Fill(nsigTOF_K_2);

    hnsigTPC_Pi_3->Fill(nsigTPC_Pi_3);
    hnsigTPC_K_3->Fill(nsigTPC_K_3);
    hnsigTOF_Pi_3->Fill(nsigTOF_Pi_3);
    hnsigTOF_K_3->Fill(nsigTOF_K_3);

    h2nsigTPC_Pi_0->Fill(pt_prong0, nsigTPC_Pi_0);
    h2nsigTPC_K_0->Fill(pt_prong0, nsigTPC_K_0);
    h2nsigTOF_Pi_0->Fill(pt_prong0, nsigTOF_Pi_0);
    h2nsigTOF_K_0->Fill(pt_prong0, nsigTOF_K_0);

    h2nsigTPC_Pi_1->Fill(pt_prong1, nsigTPC_Pi_1);
    h2nsigTPC_K_1->Fill(pt_prong1, nsigTPC_K_1);
    h2nsigTOF_Pi_1->Fill(pt_prong1, nsigTOF_Pi_1);
    h2nsigTOF_K_1->Fill(pt_prong1, nsigTOF_K_1);

    h2nsigTPC_Pi_2->Fill(pt_prong2, nsigTPC_Pi_2);
    h2nsigTPC_K_2->Fill(pt_prong2, nsigTPC_K_2);
    h2nsigTOF_Pi_2->Fill(pt_prong2, nsigTOF_Pi_2);
    h2nsigTOF_K_2->Fill(pt_prong2, nsigTOF_K_2);

    h2nsigTPC_Pi_3->Fill(pt_prong3, nsigTPC_Pi_3);
    h2nsigTPC_K_3->Fill(pt_prong3, nsigTPC_K_3);
    h2nsigTOF_Pi_3->Fill(pt_prong3, nsigTOF_Pi_3);
    h2nsigTOF_K_3->Fill(pt_prong3, nsigTOF_K_3);
  }

  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);
  c1->cd(1); hnsigTPC_Pi_0->Draw("hist"); c1->cd(2); hnsigTPC_K_0->Draw("hist"); 
  c1->cd(3); hnsigTOF_Pi_0->Draw("hist"); c1->cd(4); hnsigTOF_K_0->Draw("hist"); 

  TCanvas* c2 = new TCanvas("c2","c2",800,800);
  c2->Divide(2,2);
  c2->cd(1); hnsigTPC_Pi_1->Draw("hist"); c2->cd(2); hnsigTPC_K_1->Draw("hist"); 
  c2->cd(3); hnsigTOF_Pi_1->Draw("hist"); c2->cd(4); hnsigTOF_K_1->Draw("hist"); 

  TCanvas* c3 = new TCanvas("c3","c3",800,800);
  c3->Divide(2,2);
  c3->cd(1); hnsigTPC_Pi_2->Draw("hist"); c3->cd(2); hnsigTPC_K_2->Draw("hist");
  c3->cd(3); hnsigTOF_Pi_2->Draw("hist"); c3->cd(4); hnsigTOF_K_2->Draw("hist");

  TCanvas* c4 = new TCanvas("c4","c4",800,800);
  c4->Divide(2,2);
  c4->cd(1); hnsigTPC_Pi_3->Draw("hist"); c4->cd(2); hnsigTPC_K_3->Draw("hist");
  c4->cd(3); hnsigTOF_Pi_3->Draw("hist"); c4->cd(4); hnsigTOF_K_3->Draw("hist");

  TCanvas* c21 = new TCanvas("c21","c21",800,800);
  c21->Divide(2,2);
  c21->cd(1); h2nsigTPC_Pi_0->Draw("COLZ"); c21->cd(2); h2nsigTPC_K_0->Draw("COLZ");
  c21->cd(3); h2nsigTOF_Pi_0->Draw("COLZ"); c21->cd(4); h2nsigTOF_K_0->Draw("COLZ");

  TCanvas* c22 = new TCanvas("c22","c22",800,800);
  c22->Divide(2,2);
  c22->cd(1); h2nsigTPC_Pi_1->Draw("COLZ"); c22->cd(2); h2nsigTPC_K_1->Draw("COLZ");
  c22->cd(3); h2nsigTOF_Pi_1->Draw("COLZ"); c22->cd(4); h2nsigTOF_K_1->Draw("COLZ");

  TCanvas* c23 = new TCanvas("c23","c23",800,800);
  c23->Divide(2,2);
  c23->cd(1); h2nsigTPC_Pi_2->Draw("COLZ"); c23->cd(2); h2nsigTPC_K_2->Draw("COLZ");
  c23->cd(3); h2nsigTOF_Pi_2->Draw("COLZ"); c23->cd(4); h2nsigTOF_K_2->Draw("COLZ");

  TCanvas* c24 = new TCanvas("c24","c24",800,800);
  c24->Divide(2,2);
  c24->cd(1); h2nsigTPC_Pi_3->Draw("COLZ"); c24->cd(2); h2nsigTPC_K_3->Draw("COLZ");
  c24->cd(3); h2nsigTOF_Pi_3->Draw("COLZ"); c24->cd(4); h2nsigTOF_K_3->Draw("COLZ");

}