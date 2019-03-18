//#include "anaEvents.hh"

//anaEvents::ananEvents( Int_t irun){
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm>    // std::sort

void tripletFit(){
 
  Int_t irunStart = 34;

  TFile *outfile = new TFile("PulseAnalysis_Results.root","recreate");

    
    TString filename; filename.Form("PulseAnalysis_%i.root",irunStart);
    TFile *f = new TFile(filename);
    cout<<"opening file "<<filename<<endl;
    if(f->IsZombie()){
      cout<<"cannot open file"<<endl;
    }
    TTree *t1 = (TTree*)f->Get("ntuplePulse");

    Float_t irun,ientry,pmtNum,sigma,nhits,charge,baseline,Sdev,startTime,peakWidth,vMax,vMaxTime,T0,sigmaP,sigmaQ; 

    t1->SetBranchAddress("irun",&irun);
    t1->SetBranchAddress("ientry",&ientry);
    t1->SetBranchAddress("pmtNum",&pmtNum);
    t1->SetBranchAddress("nhits",&nhits);
    t1->SetBranchAddress("charge",&charge);
    t1->SetBranchAddress("startTime",&startTime);
    t1->SetBranchAddress("vMax",&vMax);
    t1->SetBranchAddress("vMaxTime",&vMaxTime);
    t1->SetBranchAddress("T0",&T0);
    t1->SetBranchAddress("Sdev",&Sdev);
    t1->SetBranchAddress("peakWidth",&peakWidth);
    t1->SetBranchAddress("baseline",&baseline);

    Int_t NEntries = (Int_t)t1->GetEntries();
 
    TCanvas * c = new TCanvas("c1","c1");
    TH1D * hHits = new TH1D("Pmt_hit","PMT_arrivial_time_weighted_by_peak_height",250,0e-9,4e-6);
    hHits->Sumw2();
    TH1F * hqMax = new TH1F("QMax","Max_Charge_per_event_" ,100,0,1e-8);
    TH1F * hPull = new TH1F("pull","",hHits->GetNbinsX(),-20,20);
    TH1F * hError = new TH1F("error","error",100,0,.1); 
    TH1F * hArrivalTime = new TH1F("ArrivalTime","ArrivalTime",200,0,4e-6);

    outfile->cd();
    
    Int_t currentID = 0,pastID = -1;
    Double_t peakMax = -9999, maxTime = -9999;
    std::vector<Double_t> deltaQ,Q;
    deltaQ.resize(hHits->GetNbinsX());Q.resize(hHits->GetNbinsX());
    for(int i = 0; i < NEntries; i++){
      t1->GetEntry(i);
      if(charge < 0) continue;
      //if(pmtNum != 1) continue;
      //if(sigma != 3) continue;
      //if(ientry < 295 && irun == 25) continue;
      //if(nhits < 3) continue;
      //if( Sdev*peakWidth/charge > 1 || Sdev*peakWidth/charge <= 0) continue;
      //if(Sdev/sigmaP > 1) continue;
      //if(sigmaQ/charge > 10) continue;
      //if(charge < 2e-11) continue;
      if(vMax < 1.e-3) continue;
            
      currentID = ientry;
      
      if(pastID != currentID){
        if(currentID != 0) hqMax->Fill(peakMax);
        peakMax = charge;
        pastID = currentID;
      }
      else{
        if(peakMax < vMax) peakMax = charge;
      }
      if(startTime-T0 > maxTime) maxTime = startTime-T0;
      Int_t hitBin = hHits->FindBin(startTime);
      Double_t err;
      err = sqrt(charge*1.43e-10 + std::pow(Sdev*peakWidth,2) );
      deltaQ[hitBin] += 1./(err*err);
      Q[hitBin]      += (charge/(err*err) );
      hArrivalTime->Fill(startTime-T0);
      //hHits->Fill(startTime);
      //hHits->SetBinContent(hitBin,hHits->GetBinContent(hitBin)+charge);//std::floor(qMax/Sdev) );
      //hHits->SetBinError(hitBin,hHits->GetBinError(hitBin)+ Sdev*peakWidth*Sdev*peakWidth);
      //hHits->SetBinContent(hitBin,hHits->GetBinContent(hitBin)+qMax);
      //hHits->SetBinError(hitBin,hHits->GetBinError(hitBin)+ Sdev*Sdev);
    }
    for( int i = 0; i < hHits->GetNbinsX(); i++){
      if(Q[i] == 0 ) continue;
      Q[i] = Q[i]/deltaQ[i];
      deltaQ[i] = 1./std::sqrt(deltaQ[i]);
      hHits->SetBinContent(i,Q[i]);
      hHits->SetBinError(i,deltaQ[i]);
      hError->Fill(deltaQ[i]/Q[i] );
      hArrivalTime->SetBinError(i,std::sqrt(hArrivalTime->GetBinContent(i)) );
    }
    /*  
       TF1  *f1 = new TF1("myfit","[0]*exp(-x*[1]) + [2]*exp(-x*[3])+ [4]",0.,maxTime);
       //f1->SetParameter(1,1./6e-9);
       f1->SetParLimits(1,1./40e-9,1./3e-9);
       f1->SetParLimits(0,0,10);
       f1->SetParLimits(3,1./1.e-6,1./2.25e-7);
       //f1->SetParameter(3,1./.1e-6);
       //f1->SetParameter(4,2.2e-3);
       hHits->Fit("myfit","","",0,maxTime);
       if(f1->GetParameter(1) != 0 && f1->GetParameter(3) != 0){
         Double_t t_shirt_error = f1->GetParError(1)/(f1->GetParameter(1)*f1->GetParameter(1) );
         Double_t t_long_error = f1->GetParError(3)/(f1->GetParameter(3)*f1->GetParameter(3) );
         std::cout<<"T_shirt "<<1/f1->GetParameter(1)<<" +/- "<< t_shirt_error<<", T_long "<<1/f1->GetParameter(3)<<" +/- "<<t_long_error<<std::endl;
       }
    */ 
    ///*
    TF1 *f1 = new TF1("myfit","[0]*exp(-x*[1])+[2]",1e-9,maxTime);
    f1->SetParameter(1,1./1e-6);
    hHits->Fit("myfit","","",60e-9,maxTime);
    Double_t t_long_error = f1->GetParError(1)/(f1->GetParameter(1)*f1->GetParameter(1) );
    std::cout<<"Run = "<<"  "<<1/f1->GetParameter(1)<<" +/- "<< t_long_error<<" Max Bin "<<hqMax->GetBinCenter(hqMax->GetMaximumBin() )<<" ChiSquare "<< 
      f1->GetChisquare()<<std::endl;
    //*/
    
    hHits->GetYaxis()->SetTitle("Arrival time weighted with Peak Height");
    //hHits->GetYaxis()->SetTitle("integrated charge");
    hHits->GetXaxis()->SetTitle("time (s)");
    gStyle->SetOptFit(1111);
    hHits->Draw();
    c->SetLogy();
    hHits->Write();
    hError->Write();
    //    hError->Draw();
    for(int i = 0 ; i < hHits->GetNbinsX();i++){
      if(deltaQ[i] == 0) continue;
      Double_t time = hHits->GetBinCenter(i);
      Double_t func = f1->Eval(time);
//      hPull->SetBinContent(i,(Q[i]-func)/deltaQ[i]);
        hPull->Fill((Q[i]-func)/deltaQ[i] );
    }
    hPull->Write();
    hArrivalTime->Write();
    //TCanvas * c2 = new TCanvas("c2","c2");
    //hqMax->Draw();
  //new TBrowser();
  
outfile->Write();
}
