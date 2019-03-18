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

void tripletFit400ns(){
 
  Int_t irunStart = 34;

  TFile *outfile = new TFile("PulseAnalysis_Results.root","recreate");
  Double_t time [] = {0,0.885,1.181,3.013,4.267,5.267,6.893/*,7.881*/,7.989,8.840,9.840,10.254,10.99,11.80,12.766,13.24,14.01,14.246,15.22,16.11,16.83,17.25,17.31,18.05,18.52,19.03,19.55,19.99,20.5,21.13,21.55,33.03,33.5};
  Int_t Ntimes = sizeof(time)/sizeof(Double_t) - 1;
  std::vector<Double_t> triplet,tripletErr;
    
    //TString filename; filename.Form("PulseAnalysis_%i.root",irunStart);
    //TString filename = "PulseAnalysisLAr_1_1.root";
    //TString filename = "Out_run_3.root";
    //TString filename = "Out_run_4.root";
    //TString filename = "Out_run_5.root";
    //TString filename = "Out_run_6.root";
  for(int z =3;z<=35;z++){
    if( z == 10) continue;
    if( z == 26) continue; 
    //TString filename = TString("dataSet1/Out_run_")+to_string(z)+TString(".root");
    //TString filename = TString("dataSet2/Out_run_")+to_string(z)+TString(".root");
    //TString filename = TString("dataSet3/Out_run_")+to_string(z)+TString(".root");
    //TString filename = "SimAnaResults.124103-20181206.root;
    //TString filename = TString("Out_run_")+to_string(z)+TString(".root");
    TString filename = TString("SimAnaResults.baconRun_10kBins_400ns_10mV_")+to_string(z)+TString(".root");
    TFile *f = new TFile(filename);
    cout<<"Run "<<z<<" opening file "<<filename<<endl;
    if(f->IsZombie()){
      cout<<"cannot open file"<<endl;
    }
    
    Float_t Sdev,baseline,integral,deltaT,irun;
    TTree *t2 = (TTree*)f->Get("ntupleEvent");
    
    t2->SetBranchAddress("irun",&irun);
    t2->SetBranchAddress("Sdev",&Sdev);
    t2->SetBranchAddress("baseline",&baseline);
    t2->SetBranchAddress("integral",&integral);
    t2->SetBranchAddress("deltaT",&deltaT);

    //intialize branch
    t2->GetEntry(0);
    
    TTree *t1 = (TTree*)f->Get("ntuplePulse");

    Float_t ientry,pmtNum,sigma,nhits,charge,startTime,peakWidth,vMax,vMaxTime,T0,sigmaP,sigmaQ; 

    t1->SetBranchAddress("ientry",&ientry);
    t1->SetBranchAddress("pmtNum",&pmtNum);
    t1->SetBranchAddress("nhits",&nhits);
    t1->SetBranchAddress("charge",&charge);
    t1->SetBranchAddress("startTime",&startTime);
    t1->SetBranchAddress("vMax",&vMax);
    t1->SetBranchAddress("vMaxTime",&vMaxTime);
    t1->SetBranchAddress("T0",&T0);
    //t1->SetBranchAddress("Sdev",&Sdev);
    t1->SetBranchAddress("peakWidth",&peakWidth);
    //t1->SetBranchAddress("baseline",&baseline);

    Int_t NEntries = (Int_t)t1->GetEntries();
    TString runStr = to_string(z);
    TCanvas * c = new TCanvas(TString("c1")+runStr,TString("c1") +runStr);
    TH1D * hHits = new TH1D("Pmt_hit"+runStr,"PMT_arrivial_time_weighted_by_peak_height"+runStr,200,0e-9,4e-6);
    hHits->Sumw2();
    TH1F * hqMax = new TH1F("QMax"+runStr,"Max_Charge_per_event_"+runStr ,100,0,1e-8);
    TH1F * hTotalCharge = new TH1F("TotalCharge"+runStr,"Total_Charge_per_event_"+runStr ,250,0,500);
    TH1F * hPull = new TH1F("pull"+runStr,""+runStr,hHits->GetNbinsX(),-20,20);
    TH1F * hError = new TH1F("error"+runStr,"error"+runStr,200,0,.1); 
    TH1F * hArrivalTime = new TH1F("ArrivalTime weighted by charge"+runStr,"ArrivalTime weighted by charge"+runStr,200,0,4e-6);
    TH1F * hSinglePhoton = new TH1F("SinglePhoton"+runStr,"SinglePhoton"+runStr,200,-5e-11,20e-10);

   
    TF1 *fSinglePhoton;
    if(filename.Contains("Out_run")){
    cout<<"Fitting Out Run"<<endl;  
    fSinglePhoton = new TF1("SinglePE","gaus(0) + gaus(3)",-5e-11,20e-10);
    Double_t p0_lo = 19, p0_up = 1e5;
    fSinglePhoton->SetParameter(0,(p0_lo+p0_up)/2);
    fSinglePhoton->SetParLimits(0,p0_lo,p0_up);
    Double_t p1_lo = 0, p1_up = 1e-10;
    fSinglePhoton->SetParameter(1,(p1_lo+p1_up)/2);
    fSinglePhoton->SetParLimits(1,p1_lo,p1_up);
    Double_t p2_lo = 1e-12, p2_up = 1e-11;
    fSinglePhoton->SetParameter(2,(p2_lo+p2_up)/2);
    fSinglePhoton->SetParLimits(2,p2_lo,p2_up);
    Double_t p3_lo = 10, p3_up = 1e5;
    fSinglePhoton->SetParameter(3,(p3_lo+p3_up)/2);
    fSinglePhoton->SetParLimits(3,p3_lo,p3_up);
    Double_t p4_lo = 2.5e-11, p4_up = 7e-11;
    fSinglePhoton->SetParameter(4,4.2e-11);
    fSinglePhoton->SetParLimits(4,p4_lo,p4_up);
    Double_t p5_lo = 1e-11, p5_up = 10e-11;
    fSinglePhoton->SetParameter(5,(p5_lo+p5_up)/2);
    fSinglePhoton->SetParLimits(5,p5_lo,p5_up);
    }
    else{
      fSinglePhoton = new TF1("SinglePE","gaus(0) ",2e-10,8e-10);
      Double_t p3_lo = 10, p3_up = 1e6;
      fSinglePhoton->SetParameter(0,(p3_lo+p3_up)/2);
      fSinglePhoton->SetParLimits(0,p3_lo,p3_up);
      Double_t p4_lo = 3.0e-10, p4_up = 8e-10;
      fSinglePhoton->SetParameter(1,4.2e-11);
      fSinglePhoton->SetParLimits(1,p4_lo,p4_up);
      Double_t p5_lo = 1e-11, p5_up = 10e-10;
      fSinglePhoton->SetParameter(2,(p5_lo+p5_up)/2);
      fSinglePhoton->SetParLimits(2,p5_lo,p5_up);
    }

    outfile->cd();

    for(int i = 0; i < NEntries; i++){
      t1->GetEntry(i);
      if(startTime > 3e-6)
        hSinglePhoton->Fill(charge);
    }
    if(filename.Contains("Out_run"))
      hSinglePhoton->Fit(fSinglePhoton,"Q","",-10e-12,75e-12);
    else
      hSinglePhoton->Fit(fSinglePhoton,"Q","",3e-10,8e-10);
    /*
    if(fSinglePhoton->GetParameter(4) == p4_lo ){
      fSinglePhoton->SetParLimits(4,p4_lo-p4_lo/2.,p4_up+p4_up/2.);
      hSinglePhoton->Fit(fSinglePhoton,"Q","",-10e-12,75e-12);
    }
    */
    Double_t singlePE = 0;
    Double_t singlePEerr = 0;
    if(fSinglePhoton->GetChisquare() && filename.Contains("Out_run")){
      singlePE = fSinglePhoton->GetParameter(4);
      singlePEerr = fSinglePhoton->GetParError(4);
    }
    else if(fSinglePhoton->GetChisquare()){
      singlePE = fSinglePhoton->GetParameter(1);
      singlePEerr = fSinglePhoton->GetParError(1);
    }
    cout<<"SinglePE = "<<singlePE<<"+/-"<<singlePEerr<<endl;
    Int_t currentID = 0,pastID = -1;
    Double_t peakMax = -9999, maxTime = -9999,totalCharge = 0;
    std::vector<Double_t> deltaQ,Q,dQNorm;
    deltaQ.resize(hHits->GetNbinsX());Q.resize(hHits->GetNbinsX());dQNorm.resize(hHits->GetNbinsX());
    for(int i = 0; i < NEntries; i++){
      t1->GetEntry(i);
      //if(charge < cut) continue;
      //if(charge < singlePE-singlePEerr) continue;
      //if(pmtNum != 1) continue;
      //if(sigma != 3) continue;
      //if(ientry < 295 && irun == 25) continue;
      //if(nhits < 10) continue;
      //if( Sdev*peakWidth/charge > 1 || Sdev*peakWidth/charge <= 0) continue;
      //if(Sdev/sigmaP > 1) continue;
      //if(sigmaQ/charge > 10) continue;
      //if(charge < 2e-11) continue;
      if(vMax < 0) continue;
      if(peakWidth < 0) {
        cout<<"Negative Peak Width"<<endl;
        continue;
      }
            
      currentID = ientry;
      
      if(pastID != currentID){
        if(currentID != 0) hqMax->Fill(peakMax);
        peakMax = charge;
        pastID = currentID;
        if(singlePE){
          hTotalCharge->Fill(deltaT*integral/singlePE);
        }
        else{
          hTotalCharge->Fill(deltaT*integral);
        
        }
        totalCharge = 0;
        t2->GetEntry(currentID);
      }
      else{
        if(peakMax < vMax) peakMax = charge;
      }
      if(startTime > maxTime) maxTime = startTime;
      Int_t hitBin = hHits->FindBin(startTime);
      Double_t err;
      if(Sdev == 0)Sdev = 1.e-3;
      //err like sqrt(N)^2 + uncertainty of pulse width ^2
      //err = sqrt(charge*4.93e-11 + std::pow(Sdev*peakWidth,2) );
      //err = sqrt(charge*singlePE + std::pow(Sdev*peakWidth,2) );
      err = sqrt(peakWidth/deltaT)*Sdev*deltaT;
      deltaQ[hitBin] +=err*err;
      //if(hitBin == deltaQ.size() - 1) cout<<deltaQ[hitBin]<<" "<<sqrt(peakWidth/deltaT)<<" "<<Sdev<<endl;
      dQNorm[hitBin] += 1./err;
      //deltaQ[hitBin] += charge*charge;//1./(err*err);
      //deltaQ[hitBin] +=err*err;
      Q[hitBin]      += (charge/err) ;
      hArrivalTime->Fill(startTime,charge);
      totalCharge+=charge;
      //hError->Fill(err/charge);
      //hHits->Fill(startTime);
      //hHits->SetBinContent(hitBin,hHits->GetBinContent(hitBin)+charge);//std::floor(qMax/Sdev) );
      //hHits->SetBinError(hitBin,hHits->GetBinError(hitBin)+ Sdev*peakWidth*Sdev*peakWidth);
      //hHits->SetBinContent(hitBin,hHits->GetBinContent(hitBin)+qMax);
      //hHits->SetBinError(hitBin,hHits->GetBinError(hitBin)+ Sdev*Sdev);
    }
    for( int i = 0; i < hHits->GetNbinsX(); i++){
      //if(Q[i] == 0 ) continue;
      Q[i] = Q[i]/dQNorm[i];
      if (dQNorm[i] != 0 && singlePE != 0){
        hHits->SetBinContent(i,Q[i]/singlePE);
        hHits->SetBinError(i,1./(dQNorm[i]) );
      }
      //cout<<"hArrivalTime "<<hArrivalTime->GetBinCenter(i+1)<<" "<<Q[i]<<"+/- "<<1./(dQNorm[i])<<endl;
      //hError->Fill(sqrt(deltaQ[i])/hArrivalTime->GetBinContent(i+1));
      hError->SetBinContent(i+1,sqrt(deltaQ[i])/hArrivalTime->GetBinContent(i+1));
      //hArrivalTime->SetBinError(i+1,deltaQ[i]);
      //hArrivalTime->SetBinError(i+1,std::sqrt(deltaQ[i]));
      //cout<<"hArrivalTime "<<hArrivalTime->GetBinCenter(i+1)<<" "<<hArrivalTime->GetBinContent(i+1)<<" +/- "<<sqrt(deltaQ[i])<<endl;
      //hArrivalTime->SetBinError(i,std::sqrt(hArrivalTime->GetBinContent(i)) );
    }

    ///*
    TH1D *hChi = new TH1D("Chi/NDF"+runStr,"Chi/NDF"+runStr,200,0,4e-6);
    TH1D *hTriplet = new TH1D("Triplet"+runStr,"Triplet"+runStr,200,0,4e-6);
    Double_t bestChi = 9999999,bestTrip = 0,bestErr,bestStart;
    for(int i = 0;i <40;i++){
      //TF1 *f1 = new TF1("myfit","[0]*exp(-x*[1])+[2]",1e-9,maxTime);
      TF1 *f1 = new TF1("myfit","[0]*exp(-x*[1])",hArrivalTime->GetMaximumBin(),maxTime);
      f1->SetParameter(1,1./1e-6);
      Double_t shiftTime = 100e-9; 
      Double_t startFitTime = hArrivalTime->GetBinCenter(hArrivalTime->GetMaximumBin())+shiftTime*(i+1);
      if(maxTime-startFitTime < 2e-6 ) break;
      //cout<<"start Time "<<startFitTime<<", maxTime "<<maxTime<<endl;
      hArrivalTime->Fit("myfit","Q","",startFitTime,maxTime);
      Double_t t_long_error = f1->GetParError(1)/(f1->GetParameter(1)*f1->GetParameter(1) );
      //cout<<f1->GetParameter(1)<<" "<<f1->GetParError(1)<<endl;
      //std::cout<<"Run = "<<"  "<<1/f1->GetParameter(1)<<" +/- "<< t_long_error<<" Max Bin "<<hqMax->GetBinCenter(hqMax->GetMaximumBin() )<<" ChiSquare "<< 
      // f1->GetChisquare()<<"/"<<f1->GetNDF()<<std::endl;
      hChi->Fill(startFitTime,f1->GetChisquare()/f1->GetNDF());
      hTriplet->Fill(startFitTime,1/f1->GetParameter(1));
      if(fabs(f1->GetChisquare()/f1->GetNDF() - 1) < fabs(bestChi-1) ){
        bestChi = f1->GetChisquare()/f1->GetNDF();
        bestTrip = 1/f1->GetParameter(1);
        bestErr = t_long_error;
        bestStart = startFitTime;
      }
      //*/
    }
    triplet.push_back(bestTrip);
    tripletErr.push_back(bestErr);
    cout<<bestTrip<<" +/- "<<bestErr<<" Chi/ndf = "<<bestChi<<" startTime = "<<bestStart<<endl;
    //hHits->GetYaxis()->SetTitle("Arrival time weighted with Peak Height");
    //hHits->GetYaxis()->SetTitle("integrated charge");
    hHits->GetXaxis()->SetTitle("time (s)");
    gStyle->SetOptFit(1111);
    hArrivalTime->Draw();
    c->SetLogy();
    hHits->Write();
    hError->Write();
    /*
    //    hError->Draw();
    for(int i = 0 ; i < hHits->GetNbinsX();i++){
      if(deltaQ[i] == 0) continue;
      Double_t time = hHits->GetBinCenter(i);
      Double_t func = f1->Eval(time);
//      hPull->SetBinContent(i,(Q[i]-func)/deltaQ[i]);
        hPull->Fill((Q[i]-func)/deltaQ[i] );
    }
    */

    hSinglePhoton->Write();
    //hChi->Write();
    //hTriplet->Write();
    hTotalCharge->Scale(1./hTotalCharge->GetEntries());
    hTotalCharge->Write();
    hPull->Write();
    hArrivalTime->Write();
    delete c;
  }
 
  ///*
  TH1D * hTripletvTime = new TH1D("Triplet_versus_Time","Triplet_versus_Time",Ntimes,time);
  for(int i = 0; i < Ntimes; i++){
    hTripletvTime->SetBinContent(i+1,triplet[i]);
    hTripletvTime->SetBinError(i+1,tripletErr[i]);
  }
  hTripletvTime->Draw();
  //hTripletvTime->Write();
  //*/
    //TCanvas * c2 = new TCanvas("c2","c2");
    //hqMax->Draw();
  //new TBrowser();
  
outfile->Write();
}
