#include <Riostream.h>
#include <TFile.h>
#include <AliRDHFCutsDplustoKpipi.h>
#include <TClonesArray.h>
#include <TParameter.h>

/*
 whichCuts=0, nameCuts="DplustoKpipiFilteringCuts"
 whichCuts=1, nameCuts="DplustoKpipiAnalysisCuts"
 */

AliRDHFCutsDplustoKpipi *makeInputCutsDplustoKpipi(Int_t whichCuts=0, TString nameCuts="DplustoKpipiFilteringCuts", Float_t minc=0.,Float_t maxc=100., Bool_t isMC=kFALSE, Int_t OptPreSelect = 1, Int_t TPCClsPID = 50, Bool_t PIDcorrection=kTRUE, Bool_t usePID=kTRUE)
{
  
  cout << "\n\033[1;31m--Warning (08/06/20)--\033[0m\n";
  cout << "  Don't blindly trust these cuts." << endl;
  cout << "  Relatively old and never tested." << endl;
  cout << "\033[1;31m----------------------\033[0m\n\n";

  AliRDHFCutsDplustoKpipi* cuts=new AliRDHFCutsDplustoKpipi();
  cuts->SetName(nameCuts.Data());
  cuts->SetTitle(nameCuts.Data());
  
  //UPDATE 21/06/19, use the same track quality cuts for filtering and analysis cuts
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //Should not use SetMinNClustersTPC anymore, not well described in MC
  //Two lines below replace this cut (for value 70)
  //  esdTrackCuts->SetMinNClustersTPC(80);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetMinNCrossedRowsTPC(70);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.6,1.e10);
  esdTrackCuts->SetMaxDCAToVertexXY(1.);
  esdTrackCuts->SetMaxDCAToVertexZ(1.);
  esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0060*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
  cuts->AddTrackCuts(esdTrackCuts);
  
  cuts->SetScaleNormDLxyBypOverPt(kFALSE);
  cuts->SetRemoveTrackletOutliers(kFALSE);
  //TODO: Should this be false or true for ITS Upgrade?
  cuts->SetUseTrackSelectionWithFilterBits(kFALSE);
  
  //UPDATE 08/06/20, Add cut on TPC clusters for PID (similar to geometrical cut)
  cuts->SetMinNumTPCClsForPID(TPCClsPID);

  if(whichCuts==0){
    const Int_t nptbins=2;
    Float_t* ptbins;
    ptbins=new Float_t[nptbins+1];
    
    ptbins[0]=2.;
    ptbins[1]=5.;
    ptbins[2]=1.e+9;
    
    const Int_t nvars=14;
    Float_t** anacutsval;
    anacutsval=new Float_t*[nvars];
    for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
    
    anacutsval[0][0]=0.2;
    anacutsval[1][0]=0.6;
    anacutsval[2][0]=0.6;
    anacutsval[3][0]=0.;
    anacutsval[4][0]=0.;
    anacutsval[5][0]=0.;
    anacutsval[6][0]=0.035;
    anacutsval[7][0]=0.08;
    anacutsval[8][0]=0.0;
    anacutsval[9][0]=0.98;
    anacutsval[10][0]=0.0;
    anacutsval[11][0]=10000000000.;
    anacutsval[12][0]=5.;
    anacutsval[13][0]=0.99;
    
    anacutsval[0][1]=0.25;
    anacutsval[1][1]=0.6;
    anacutsval[2][1]=0.6;
    anacutsval[3][1]=0.;
    anacutsval[4][1]=0.;
    anacutsval[5][1]=0.01;
    anacutsval[6][1]=0.05;
    anacutsval[7][1]=0.05;
    anacutsval[8][1]=0.0;
    anacutsval[9][1]=0.95;
    anacutsval[10][1]=0.0;
    anacutsval[11][1]=10000000000.;
    anacutsval[12][1]=3.;
    anacutsval[13][1]=0.;
    
    cuts->SetPtBins(nptbins+1,ptbins);
    cuts->SetCuts(nvars,nptbins,anacutsval);
    
    cuts->SetMinPtCandidate(2.);
    cuts->SetMaxPtCandidate(1000.);
    
    Bool_t pidflag=usePID;
    cuts->SetUsePID(pidflag);
    if(pidflag) cout<<"PID is used for filtering cuts"<<endl;
    else cout<<"PID is not used for filtering cuts"<<endl;
  }
  else if(whichCuts==1){
    
    const Int_t nptbins=15;
    Float_t* ptbins;
    ptbins=new Float_t[nptbins+1];
    
    ptbins[0]=2.;
    ptbins[1]=3.;
    ptbins[2]=4.;
    ptbins[3]=5.;
    ptbins[4]=6.;
    ptbins[5]=7.;
    ptbins[6]=8.;
    ptbins[7]=9.;
    ptbins[8]=10.;
    ptbins[9]=11.;
    ptbins[10]=12.;
    ptbins[11]=14.;
    ptbins[12]=16.;
    ptbins[13]=24.;
    ptbins[14]=36.;
    ptbins[15]=50.;
    
    const Int_t nvars=14;
    Float_t** anacutsval;
    anacutsval=new Float_t*[nvars];
    for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
    
    
    Int_t ic=0;//minv
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[ic][ipt]=0.2;
    }
    
    ic=1;//ptK
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[1][ipt]=0.6;
    }
    
    ic=2;//ptPi
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[ic][ipt]=0.6;
    }
    
    ic=3;//d0K
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[ic][ipt]=0.;
    }
    ic=4;//d0Pi
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[ic][ipt]=0.;
    }
    ic=5;//dist12
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[ic][ipt]=0.;
    }
    
    ic=6;//sigvert
    anacutsval[ic][0]=0.020;//2.0-3.0
    anacutsval[ic][1]=0.020;//3.0-4.0
    anacutsval[ic][2]=0.020;//4.0-5.0
    anacutsval[ic][3]=0.020;//5.0-6.0
    anacutsval[ic][4]=0.020;//6.0-7.0
    anacutsval[ic][5]=0.024;//7.0-8.0
    anacutsval[ic][6]=0.024;//8.0-9.0
    anacutsval[ic][7]=0.024;//9.0-10.0
    anacutsval[ic][8]=0.024;//10.0-11.0
    anacutsval[ic][9]=0.024;//11.0-12.0
    anacutsval[ic][10]=0.024;//12.0-14.0
    anacutsval[ic][11]=0.024;//14.0-16.0
    anacutsval[ic][12]=0.024;//16.0-24.0
    anacutsval[ic][13]=0.034;//24.0-36.0
    anacutsval[ic][14]=0.040;//36.0-50.0
    
    ic=7;//declen
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[ic][ipt]=0.14;
    }
    anacutsval[ic][0]=0.13;//2.0-3.0
    anacutsval[ic][1]=0.13;//3.0-4.0
    anacutsval[ic][2]=0.13;//4.0-5.0
    anacutsval[ic][3]=0.13;//5.0-6.0
    anacutsval[ic][10]=0.16;//12.0-14.0
    anacutsval[ic][11]=0.16;//14.0-16.0
    anacutsval[ic][12]=0.18;//16.0-24.0
    anacutsval[ic][13]=0.20;//24.0-36.0
    anacutsval[ic][14]=0.20;//36.0-50.0
    
    ic=8;//pM
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[ic][ipt]=0.0;
    }
    
    //cosp
    ic=9;
    anacutsval[ic][0]=0.998;//2.0-3.0
    anacutsval[ic][1]=0.998;//3.0-4.0
    anacutsval[ic][2]=0.998;//4.0-5.0
    anacutsval[ic][3]=0.997;//5.0-6.0
    anacutsval[ic][4]=0.997;//6.0-7.0
    anacutsval[ic][5]=0.995;//7.0-8.0
    anacutsval[ic][6]=0.992;//8.0-9.0
    anacutsval[ic][7]=0.992;//9.0-10.0
    anacutsval[ic][8]=0.992;//10.0-11.0
    anacutsval[ic][9]=0.992;//11.0-12.0
    anacutsval[ic][10]=0.990;//12.0-14.0
    anacutsval[ic][11]=0.990;//14.0-16.0
    anacutsval[ic][12]=0.990;//16.0-24.0
    anacutsval[ic][13]=0.990;//24.0-36.0
    anacutsval[ic][14]=0.990;//36.0-50.0
    
    ic=10;//sumd02
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[ic][ipt]=0.0;
    }
    
    ic=11;//dca
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      anacutsval[ic][ipt]=10000000000.;
    }
    
    ic=12;//ndlXY
    anacutsval[ic][0]=15.;//2.0-3.0
    anacutsval[ic][1]=14.;//3.0-4.0
    anacutsval[ic][2]=12.;//4.0-5.0
    anacutsval[ic][3]=12.;//5.0-6.0
    anacutsval[ic][4]=12.;//6.0-7.0
    anacutsval[ic][5]=12.;//7.0-8.0
    anacutsval[ic][6]=12.;//8.0-9.0
    anacutsval[ic][7]=12.;//9.0-10.0
    anacutsval[ic][8]=10.;//10.0-11.0
    anacutsval[ic][9]=10.;//11.0-12.0
    anacutsval[ic][10]=10.;//12.0-14.0
    anacutsval[ic][11]=10.;//14.0-16.0
    anacutsval[ic][12]=10.;//16.0-24.0
    anacutsval[ic][13]=10.;//24.0-36.0
    anacutsval[ic][14]=10.;//36.0-50.0
    
    ic=13;//cospXY
    anacutsval[ic][0]=0.999;//2.0-3.0
    anacutsval[ic][1]=0.999;//3.0-4.0
    anacutsval[ic][2]=0.999;//4.0-5.0
    anacutsval[ic][3]=0.998;//5.0-6.0
    anacutsval[ic][4]=0.998;//6.0-7.0
    anacutsval[ic][5]=0.996;//7.0-8.0
    anacutsval[ic][6]=0.994;//8.0-9.0
    anacutsval[ic][7]=0.994;//9.0-10.0
    anacutsval[ic][8]=0.994;//10.0-11.0
    anacutsval[ic][9]=0.994;//11.0-12.0
    anacutsval[ic][10]=0.990;//12.0-14.0
    anacutsval[ic][11]=0.990;//14.0-16.0
    anacutsval[ic][12]=0.990;//16.0-24.0
    anacutsval[ic][13]=0.990;//24.0-36.0
    anacutsval[ic][14]=0.990;//36.0-50.0
    
    Float_t *d0cutsval=new Float_t[nptbins];
    for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0
      d0cutsval[ipt]=60;
    }
    d0cutsval[0]=90;
    d0cutsval[12]=80;
    d0cutsval[13]=80;
    d0cutsval[14]=80;
    
    Float_t *d0d0expcutsval=new Float_t[nptbins];
    for(Int_t ipt=0;ipt<nptbins;ipt++){ //d0d0exp
      d0d0expcutsval[ipt]=2.0;
    }
    d0d0expcutsval[13]=2.5;
    d0d0expcutsval[14]=2.5;
    
    cuts->SetPtBins(nptbins+1,ptbins);
    cuts->SetCuts(nvars,nptbins,anacutsval);
    cuts->Setd0Cut(nptbins,d0cutsval);
    cuts->Setd0MeasMinusExpCut(nptbins,d0d0expcutsval);
    
    //        cuts->SetUseStrongPid(3);
    //        cuts->SetMaxPStrongPidK(1);
    //        cuts->SetMaxPStrongPidpi(1);
    cuts->SetUseImpParProdCorrCut(kFALSE);
    
    cuts->SetMinPtCandidate(2.);
    cuts->SetMaxPtCandidate(50.);
    
    Bool_t pidflag=usePID;
    cuts->SetUsePID(pidflag);
    if(pidflag) cout<<"PID is used for analysis cuts"<<endl;
    else cout<<"PID is not used for analysis cuts"<<endl;
  }
  
  //UPDATE 08/06/20: PreSelect, acting before FillRecoCasc.
  //NOTE: actual function not implemented for all HF hadrons yet (please check)
  cuts->SetUsePreSelect(OptPreSelect);

  //Do not recalculate the vertex
  cuts->SetRemoveDaughtersFromPrim(kFALSE); //activate for pp
  
  //Temporary PID fix for 2018 PbPb (only to be used on data)
  if(!isMC && PIDcorrection) cuts->EnableNsigmaDataDrivenCorrection(kTRUE, AliAODPidHF::kPbPb010);

  //event selection
  cuts->SetUsePhysicsSelection(kFALSE);
  cuts->SetTriggerClass("");
  if(!isMC) cuts->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kCentral);
  else      cuts->SetTriggerMask(AliVEvent::kAny);
  cuts->SetTriggerMask(AliVEvent::kAny);
  cuts->SetOptPileup(AliRDHFCuts::kNoPileupSelection);
  cuts->SetMaxVtxZ(10.);
  cuts->SetCutOnzVertexSPD(3);
  
  cout<<"This is the object I'm going to save:"<<endl;
  cuts->SetName(nameCuts.Data());
  cuts->SetTitle(nameCuts.Data());
  cuts->PrintAll();
  
  return cuts;
}
