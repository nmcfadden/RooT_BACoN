void readNtuple(){
 
  TString filename = "baconRunAna34.340000.root";//
  
  TFile *f = new TFile(filename);
  TTree *t1 = (TTree*)f->Get("ntupleRun");
  //Float_t irun,ievent,integral,integralOff,QMax,noiseQMax,ChargeMaxSig3,ChargeMaxOffSig3,ChargeMaxSig5,ChargeMaxOffSig5,baseline,Sdev;
  Float_t irun,ievent,pmtNum,sigma,charge,baseline,Sdev,startTime,peakWidth,qMax,nhits,qMaxTime;
 
  t1->SetBranchAddress("irun",&irun);
  t1->SetBranchAddress("ientry",&ievent);
  t1->SetBranchAddress("pmtNum",&pmtNum);
  t1->SetBranchAddress("sigma",&sigma);
  t1->SetBranchAddress("charge",&charge);
  t1->SetBranchAddress("startTime",&startTime);
  t1->SetBranchAddress("peakWidth",&peakWidth);
  t1->SetBranchAddress("qMaxTime",&qMaxTime);
  t1->SetBranchAddress("qMax",&qMax);
  t1->SetBranchAddress("nhits",&nhits);

  Double_t qStart = 9.1e-3;
  Double_t qWindow = 2.4e-3;


  TH1F * htotal = new TH1F("htotal","",100,0,6e-10);
  TH1F * h1 = new TH1F("h1","",100,0,6e-10);
  TH1F * h2 = new TH1F("h2","",100,0,6e-10);
  TH1F * h3 = new TH1F("h3","",100,0,6e-10);
  TH1F * h4 = new TH1F("h4","",100,0,6e-10);
  TH1F * h5 = new TH1F("h5","",100,0,6e-10);
  TH1F * h6 = new TH1F("h6","",100,0,6e-10);
   
  Int_t Nentries = (Int_t)t1->GetEntries();

  for(int i = 0; i < Nentries; i++){
    t1->GetEntry(i);
    if(pmtNum != 1) continue;
    if(sigma > 3.3 && sigma < 3.1) continue;
    //if(sigma > 3.1) continue;
    htotal->Fill(charge);
    
  }

  //TF1  *f1 = new TF1("myfit","[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]) ) + [3]*exp(-0.5*(x-[4])*(x-[4])/([5]*[5]) + [6]*exp(-0.5*(x-[7])*(x-[7])/([8]*[8])",0.7e-10,3.8e-10);
  TF1  *f1 = new TF1("myfit","gaus(0)+gaus(3)+gaus(6)+gaus(9)",0.5e-10,5.5e-10);
  f1->SetParameter(0,1.4e4);
  f1->SetParLimits(0,.5e4,5e4);

  f1->SetParameter(1,1.05e-10);
  f1->SetParLimits(1,.8e-10,1.2e-10);

  f1->SetParameter(2,6.9e-11);
  f1->SetParLimits(2,1e-11,2e-10);

  f1->SetParameter(3,1e4);
  f1->SetParLimits(3,1e3,2e4);

  f1->SetParameter(4,2.1e-10);
  f1->SetParLimits(4,1.8e-10,2.2e-10);

  f1->SetParameter(5,.68e-10);
  f1->SetParLimits(5,1e-11,2e-10);

  f1->SetParameter(6,6e3);
  f1->SetParLimits(6,1e3,1e4);

  f1->SetParameter(7,3.15e-10);
  f1->SetParLimits(7,2.8e-10,3.2e-10);

  f1->SetParameter(8,.9e-10);
  f1->SetParLimits(8,1e-11,2e-10);

  f1->SetParameter(9,1e3);
  f1->SetParLimits(9,.5e3,4e3);

  f1->SetParameter(10,4e-10);
  f1->SetParLimits(10,3.5e-10,4.2e-10);

  f1->SetParameter(11,.9e-10);
  f1->SetParLimits(11,1e-11,2e-10);
/*
  f1->SetParameter(12,1.9e4);
  f1->SetParLimits(12,1e4,4e4);

  f1->SetParameter(13,-1.7e-12);
  f1->SetParLimits(13,-2e-12,2e12);

  f1->SetParameter(14,6e-11);
  f1->SetParLimits(14,1e-12,2e-10);

*/
  htotal->Fit("myfit","","",7e-11,5.5e-10);
  htotal->Draw();
  TF1 *pe1 = new TF1("pe1","gaus(0)");
  pe1->SetParameter(0,f1->GetParameter(0));
  pe1->SetParameter(1,f1->GetParameter(1));
  pe1->SetParameter(2,f1->GetParameter(2));
  pe1->Draw("same");

  TF1 *pe2 = new TF1("pe2","gaus(0)");
  pe2->SetParameter(0,f1->GetParameter(3));
  pe2->SetParameter(1,f1->GetParameter(4));
  pe2->SetParameter(2,f1->GetParameter(5));
  pe2->Draw("same");

  TF1 *pe3 = new TF1("pe3","gaus(0)");
  pe3->SetParameter(0,f1->GetParameter(6));
  pe3->SetParameter(1,f1->GetParameter(7));
  pe3->SetParameter(2,f1->GetParameter(8));
  pe3->Draw("same");

  TF1 *pe4 = new TF1("pe4","gaus(0)");
  pe4->SetParameter(0,f1->GetParameter(9));
  pe4->SetParameter(1,f1->GetParameter(10));
  pe4->SetParameter(2,f1->GetParameter(11));
  pe4->Draw("same");
} 
