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
#include <string>

Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}



TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");

   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }


   //his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
   his->Fit(FunName,"Q");   // fit within specified range, use ParLimits, plot, do not report
   
   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Seaches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);

      if (l < lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}

void tripletFit1us(){
 
  Int_t irunStart = 10000;
  Int_t irunStop  = 10784;
  bool simFlag = false;

  TString outFileName = TString("1usTripletFit");
  if( irunStart == 100 && irunStop == 272) outFileName = outFileName + TString("_100ppmN2.root");
  else if( irunStart == 3000 && irunStop < 4000) outFileName = outFileName + TString("_30ShortRuns.root");
  else if( irunStart == 4000 && irunStop < 5000) outFileName = outFileName + TString("_PulseTrigger.root");
  else if (irunStart == 10000 && irunStop < 11000) outFileName = outFileName + TString("_20ppmN2.root");
  else outFileName = outFileName + TString("Temp.root");
  TFile *outfile = new TFile(outFileName ,"recreate");
  //Double_t time [] = {0,0.885,1.181,3.013,4.267,5.267,6.893/*,7.881*/,7.989,8.840,9.840,10.254,10.99,11.80,12.766,13.24,14.01,14.246,15.22,16.11,16.83,17.25,17.31,18.05,18.52,19.03,19.55,19.99,20.5,21.13,21.55,33.03,33.5};
  std::vector<Double_t> monthVector = {0,31,59,90,120,151,181,212,243,273,304,334,365};
  std::vector<Double_t> time,totalCharge,totalChargeErr,promptTime,promptTimeError,intermediateTime,intermediateTimeError,longTime,longTimeError,singleCharge,singleChargeError,A_s,A_l;
  //time.push_back(0);
  std::vector<Double_t> triplet,tripletErr;
  std::vector<Double_t> runNumber;
  Double_t globalStartTime = -999;  
    //TString filename; filename.Form("PulseAnalysis_%i.root",irunStart);
    //TString filename = "PulseAnalysisLAr_1_1.root";
    //TString filename = "Out_run_3.root";
    //TString filename = "Out_run_4.root";
    //TString filename = "Out_run_5.root";
    //TString filename = "Out_run_6.root";
  bool cutDebug = false;
  for(int z =irunStart;z<=irunStop;z++){

    if( z == 10117) continue;
    if( z == 10274) continue;
    if( z == 10330) continue;
    if( z == 10353) continue;
    if( z == 10358) continue;
    if( z == 5002) continue;
    if( z == 2046) continue;
    if( z == 185) continue;
    if( z == 159) continue;
    //these runs are for cobalt
    if( z >= 262 && z <= 271) continue;
    if( z == 107) continue;
    if( z == 108) continue;
    if( z == 109) continue;
    /*
    if( z == 107) continue;
    if( z == 108) continue; 
    if( z == 109) continue; 
    */
    //TString filename = TString("dataSet1/Out_run_")+to_string(z)+TString(".root");
    //TString filename = TString("dataSet2/Out_run_")+to_string(z)+TString(".root");
    //TString filename = TString("dataSet3/Out_run_")+to_string(z)+TString(".root");
    //TString filename = "SimAnaResults.124103-20181206.root;
    //TString filename = TString("Out_run_")+to_string(z)+TString(".root");
    TString filename;
    if(z>= 100 && z < 2000)
      filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_10mV_")+to_string(z)+TString(".root");
    else if(z >= 2000 && z <4000)
      filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_20mV_div_-30mV_thresh_")+to_string(z)+TString(".root");
    else if(z >= 4000 && z <6000)
      filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_20mV_div_-6mV_thresh_")+to_string(z)+TString(".root");
    else if(z>=10000)
      filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_")+to_string(z)+TString(".root");
    else{ 
      filename = "processedData/SimAnaResults.simEvents_20190308_10000.root";
      simFlag = true;
      z = irunStop;
    }
    TFile *f = new TFile(filename);
    cout<<"Run "<<z<<" opening file "<<filename<<endl;
    if(f->IsZombie()){
      cout<<"cannot open file"<<endl;
    }
   
    TNamed * eventInfo;

    f->GetObject("eventInfo",eventInfo);
    if(eventInfo == NULL){
      cout<<"Error..Warning..Fatal.."<<endl;
      cout<<"Run "<<z<<" was not properly processed"<<endl;
      cout<<"Error..Warning..Fatal.."<<endl;
      if(!simFlag) continue;
      else eventInfo = new TNamed("title","2019-02-22 13:57 0.03999996185");

    }
    runNumber.push_back(z);
    string fullTitle = eventInfo->GetTitle();
    //cout<<"eventInfo is "<<eventInfo->GetTitle()<<endl;
    ///*
    double DateYear = stod(fullTitle.substr(0,4));
    double DateMonth = stod(fullTitle.substr(5,2));
    double DateDay = stod(fullTitle.substr(8,2));
    double EndTimeHour = stod(fullTitle.substr(11,2));
    double EndTimeMinute = stod(fullTitle.substr(14,2));
    double RunTime = stod(fullTitle.substr(25) );//,fullTitle.size() - 1);
    double globalTimeStop = DateYear*365.+ monthVector[DateMonth] + DateDay + EndTimeHour/24. +EndTimeMinute/(24.*60.) + RunTime/(3600.*24.);
    double globalTimeStart = DateYear*365.+ monthVector[DateMonth] + DateDay + EndTimeHour/24. +EndTimeMinute/(24.*60.);
    double GetterStartTime = 0;
    if(z< 1000){
      //019-1-11, endTime 11:27//run 126
      GetterStartTime = 2019.*365.+monthVector[1]+11.+11./24.+10/(24.*60);
    }
    else if(z >= 10000){
      GetterStartTime = 2019.*365.+monthVector[2]+15.+12./24.+41/(24.*60);
    }
    if(globalStartTime < 0) globalStartTime = globalTimeStart;
    time.push_back(globalTimeStart-GetterStartTime);
    time.push_back(globalTimeStop-GetterStartTime);
    cout<<"Time Since run zero "<<time[time.size() -2]<<", date "<<DateYear<<"-"<<DateMonth<<"-"<<DateDay<<", endTime "<<EndTimeHour<<":"<<EndTimeMinute<<", runTime "<<RunTime<<endl;

    Float_t Sdev,baseline,integral,deltaT,irun,tMax= 0,vMaxEvent;
    TTree *t2 = (TTree*)f->Get("ntupleEvent");
    t2->SetBranchAddress("irun",&irun);
    t2->SetBranchAddress("Sdev",&Sdev);
    t2->SetBranchAddress("baseline",&baseline);
    t2->SetBranchAddress("integral",&integral);
    t2->SetBranchAddress("deltaT",&deltaT);
    t2->SetBranchAddress("tMax",&tMax);
    t2->SetBranchAddress("vMax",&vMaxEvent);
    //intialize branch
    Int_t NEvents = (Int_t)t2->GetEntries();
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
    TH1D * hHits = new TH1D("Pmt_hit"+runStr,"PMT_arrivial_time_weighted_by_peak_height"+runStr,200,0e-9,10e-6);
    TH1D * hNhits = new TH1D("N_hits"+runStr,"N_hits",200,0,200);
    hHits->Sumw2();
    TH1F * hTotalCharge;
    TH1F * hPromptTime;
    TH1F * hIntermediateTime;
    TH1F * hLongTime; 
    if(z <10000){
      hTotalCharge = new TH1F("TotalCharge"+runStr,"Total_Charge_per_event_"+runStr ,150,100,700);
      hPromptTime = new TH1F("PromptTime"+runStr,"PromptTime"+runStr,100,20,200);
      hIntermediateTime = new TH1F("IntermediateTime"+runStr,"IntermediateTime"+runStr,150,0,500);
      hLongTime = new TH1F("LongTime"+runStr,"LongTime"+runStr,150,100,700);
    }
    else{
      hTotalCharge = new TH1F("TotalCharge"+runStr,"Total_Charge_per_event_"+runStr ,150,0,450);
      hPromptTime = new TH1F("PromptTime"+runStr,"PromptTime"+runStr,100,0,150);
      hIntermediateTime = new TH1F("IntermediateTime"+runStr,"IntermediateTime"+runStr,150,0,450);
      hLongTime = new TH1F("LongTime"+runStr,"LongTime"+runStr,150,0,300);
    }
    TH1F * hArrivalTime = new TH1F("ArrivalTime weighted by charge"+runStr,"ArrivalTime weighted by charge"+runStr,400,0,10e-6);
    hArrivalTime->Sumw2();
    TH1F * hArrivalTimeUnweighted = new TH1F("ArrivalTime Unweighted"+runStr,"ArrivalTime Unweighted"+runStr,200,0,10e-6);
    TH1F * hSinglePhoton = new TH1F("SinglePhoton"+runStr,"SinglePhoton"+runStr,200,0,2e-10);
   
    Double_t singlePEThreshold = 0.99, sumSPE = 0,startTimeSPE = 7e-6;
    /*
    //Calculate where the starting point is for SPE
    //integrate last 5% of points
    Double_t sumArrivalTime = 0;
    for(int i = 0; i < NEntries; i++){
      t1->GetEntry(i);
      if(vMax < 1e-3){
        if(cutDebug)cout<<"ientry "<<ientry<<", vMax < 0 cut, vMax = "<<vMax<<endl;
        continue;
      }
      hArrivalTime->Fill(startTime,charge);
      sumArrivalTime += charge;
    }
    //cout<<"Sum is "<<sumArrivalTime<<endl;
    for(int i = 0; i < hArrivalTime->GetNbinsX(); i++){
      sumSPE += hArrivalTime->GetBinContent(i);
      if(sumSPE/(sumArrivalTime) >singlePEThreshold){
        startTimeSPE = hArrivalTime->GetBinCenter(i);
        cout<<"Found start time for SPE at "<<startTimeSPE<<endl;
        break;
      }
    }
     */
    hArrivalTime->Reset();
    TF1 *fSinglePhoton = new TF1("SinglePE","landau(0)",0,2e-9);
    outfile->cd();
    Int_t currentID = 0,pastID = -1;
    std::vector<Double_t> eventSkip;eventSkip.resize(NEvents);
    std::vector<Double_t> maxPeakTime;maxPeakTime.resize(NEvents);
    Double_t peakMax = -9999, peakTime = -9999,maxTime = -9999,Integral = 0,promptCharge = 0,intermediateCharge = 0,longCharge = 0;
    Double_t IntegralRun = 0,promptChargeRun = 0,longChargeRun = 0;
    Double_t tMaxCut = 1.25e-6,tMinCut = 0e-6,vMaxEventCut = 6e-3,vMinCut = 3e-3,peakWidthCut = 0;//25e-9;
    if(z<4000 && z >=2000) tMaxCut = 1.2e-6;
    for(int i = 0; i < NEntries; i++){
      t1->GetEntry(i);
      t2->GetEntry(ientry);
      ///*
      //pulse cuts go here
      if(vMax < vMinCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", vMax < 0 cut, vMax = "<<vMax<<endl;
        continue;
      }
      //one sigma cut on charge
      if(peakWidth < peakWidthCut) {
        if(cutDebug)cout<<"ientry "<<ientry<<", Negative Peak Width"<<endl;
        continue;
      }
      //avoid events that start late
      if(tMax > tMaxCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", cut on tMax "<<tMax<<endl;
        continue;
      }
      if(vMaxEvent < vMaxEventCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", cut on vMaxEvent "<<vMaxEvent<<endl;
        continue;
      }
      //*/
      currentID = ientry;
      //end of the event fill 
      if(pastID != currentID && currentID != 0){
        pastID = currentID;
        if(currentID != 0){
          t2->GetEntry(currentID);
        }
      }
      hNhits->Fill(nhits);
      if(startTime > startTimeSPE) 
        hSinglePhoton->Fill(charge);
    }
    hSinglePhoton->Fit(fSinglePhoton,"Q","",0.02e-9,.1e-9);
    Double_t singlePE = 0;
    Double_t singlePEerr = 0;
    singlePE = fSinglePhoton->GetParameter(1);
    singleCharge.push_back(singlePE);
    //singlePE = 2.273e-11;
    singlePEerr = fSinglePhoton->GetParameter(2);
    cout<<"SinglePE = "<<singlePE<<"+/-"<<singlePEerr<<" chi/ndf "<<fSinglePhoton->GetChisquare()/fSinglePhoton->GetNDF()<<endl;
    singleChargeError.push_back(singlePEerr);
    //taken from low triplet runs at late times
    singlePE = 5.37e-11;//fSinglePhoton->GetParameter(1);
    singlePEerr = 1.67e-11;//fSinglePhoton->GetParameter(2);
    //TF1 *hGaus = new TF1("Gaus","gaus(0)",0,200);
    //hNhits->Fit(hGaus,"Q","",0,200);

    //Double_t meanNhits = hGaus->GetParameter(1);
    //Double_t sigmaNhits = hGaus->GetParameter(2);
    std::vector<Double_t> deltaQ,Q,dQNorm;
    deltaQ.resize(hHits->GetNbinsX());Q.resize(hHits->GetNbinsX());dQNorm.resize(hHits->GetNbinsX());
    //cout<<hGaus->GetChisquare()/hGaus->GetNDF()<<endl;
    for(int i = 0; i < NEntries; i++){
      t1->GetEntry(i);
      t2->GetEntry(ientry);
      ///*
      //pulse cuts go here
      if(vMax < vMinCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", vMax < 0 cut, vMax = "<<vMax<<endl;
        continue;
      }
      //one sigma cut on charge
      if(charge < singlePE-2*singlePEerr){
        if(cutDebug)cout<<"ientry "<<ientry<<", one sigma cut on charge... "<<charge<<", less than "<<singlePE-singlePEerr<<endl;
        continue;
      }
      /*
      if(nhits < meanNhits - sigmaNhits*2 && hGaus->GetChisquare()/hGaus->GetNDF() < 1000){
        if(cutDebug)cout<<"ientry "<<ientry<<", nhits cut... "<<nhits<<" less than "<< meanNhits - sigmaNhits*2<<endl;
        continue;
      }
      */
      if(peakWidth < peakWidthCut) {
        if(cutDebug)cout<<"ientry "<<ientry<<", Negative Peak Width"<<endl;
        continue;
      }
      //avoid events that start late
      if(tMax > tMaxCut && tMax < tMinCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", cut on tMax "<<tMax<<endl;
        continue;
      }
      if(vMaxEvent < vMaxEventCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", cut on vMaxEvent "<<vMaxEvent<<endl;
        continue;
      }
      //*/

      currentID = ientry;
      //End of event Fill,
      if(pastID != currentID && currentID != 0){
        pastID = currentID;
        if(cutDebug)cout<<"ientry "<<ientry<<", tMax "<<tMax<<", vMaxEvent "<<vMaxEvent<<endl;
        if(singlePE && promptCharge > 0 && longCharge > 0){
          Integral = promptCharge+intermediateCharge+longCharge;
          hTotalCharge->Fill(Integral/singlePE);
          hPromptTime->Fill(promptCharge/singlePE);
          hIntermediateTime->Fill(intermediateCharge/singlePE);
          hLongTime->Fill(longCharge/singlePE);
        }
        Integral = 0;
        promptCharge = 0;
        intermediateCharge = 0;
        longCharge = 0;
      }
      //*/
      if(startTime > maxTime) maxTime = startTime;
      Int_t hitBin = hHits->FindBin(startTime-tMax);
      Double_t err;
      if(Sdev == 0)Sdev = 1.e-3;
      //err like sqrt(N)^2 + uncertainty of pulse width ^2
      //err = sqrt(charge*singlePE + std::pow(Sdev*peakWidth,2) );
      err = (charge);
      deltaQ[hitBin] +=charge*charge;
      dQNorm[hitBin] += 1./err;
      Q[hitBin]      += (charge/err) ;
      hArrivalTimeUnweighted->Fill(startTime,vMax);
      //hArrivalTime->Fill(startTime-tMax,charge);
      hArrivalTime->Fill(startTime-tMax,charge);
      Integral+=charge;
      IntegralRun += charge;
      if(startTime > 900e-9 - 50e-9 && startTime < 900e-9 + 100e-9){
        promptCharge += charge;
        promptChargeRun += charge;
      }
      else if( startTime > 900e-9 +100e-9 && startTime < 900e-9 +200e-9){
        intermediateCharge += charge;
      }
      else if( startTime > 900e-9 + 100e-9){
        longCharge   += charge;
        longChargeRun += charge;
      }
      if(i == NEntries -1){
        if(singlePE){
          hTotalCharge->Fill(Integral/singlePE);
          hPromptTime->Fill(promptCharge/singlePE);
          hIntermediateTime->Fill(intermediateCharge/singlePE);
          hLongTime->Fill(longCharge/singlePE);
        }
        Integral = 0;
        promptCharge = 0;
        intermediateCharge = 0;
        longCharge = 0;
      }
      //*/
    }
    
    ///*
    //TH1D *hChi = new TH1D("Chi/NDF"+runStr,"Chi/NDF"+runStr,200,0,10e-6);
    //TH1D *hTriplet = new TH1D("Triplet"+runStr,"Triplet"+runStr,200,0,10e-6);
    Double_t bestChi = 9999999,bestTrip = 0,bestErrTrip = 0,bestStart,bestSinglet = 0,bestErrSinglet = 0;
    TF1 *f1;
    /*
    for(int i = 10;i >=0;i--){
      //TF1 *f1 = new TF1("myfit","[0]*exp(-x*[1])+[2]",1e-9,maxTime);
      Double_t shiftTime = 50e-9; 
      if(z>=2000)
        f1 = new TF1("myfit","([0]/[1])*exp(-x/[1]) + ([2]/[3])*exp(-x/[3])",hArrivalTime->GetBinCenter(hArrivalTime->GetMaximumBin())+shiftTime*(i),maxTime);
      else
        f1 = new TF1("myfit","[0]*exp(-x*[1]) + [2]*exp(-x*[3]) + [4]",hArrivalTime->GetBinCenter(hArrivalTime->GetMaximumBin())+shiftTime*(i),maxTime);
      f1->SetParLimits(0,0,1e-10);
      f1->SetParameter(1,1e-6);
      f1->SetParLimits(1,0.5e-6,2e-6);
      f1->SetParLimits(2,0,1e-10);
      f1->SetParameter(3,1e-7);
      f1->SetParLimits(3,5e-9,5e-7);
      if(z < 2000)
        f1->SetParLimits(4,0,1);
      
         Double_t shiftTime = 100e-9; 
         Double_t startFitTime = hArrivalTime->GetBinCenter(hArrivalTime->GetMaximumBin())+shiftTime*(i+1);
         
      Double_t startFitTime = hArrivalTime->GetBinCenter(hArrivalTime->GetMaximumBin())+shiftTime*(i);
      //if(maxTime-startFitTime < 8e-6 ) break;
      //cout<<"start Time "<<startFitTime<<", maxTime "<<maxTime<<endl;
      hArrivalTime->Fit("myfit","Q","",startFitTime,maxTime-100e-9);
      Double_t t_short_error = f1->GetParError(1)/(f1->GetParameter(1)*f1->GetParameter(1) );
      Double_t t_long_error = f1->GetParError(3)/(f1->GetParameter(3)*f1->GetParameter(3) );
      //cout<<f1->GetParameter(1)<<" "<<f1->GetParError(1)<<endl;
      //std::cout<<"\t Run = "<<"start Time "<<startFitTime<<", maxTime "<<maxTime<<"..."<<1/f1->GetParameter(3)<<" +/- "<< t_long_error<<", Chisquared/ndf "<< f1->GetChisquare()/f1->GetNDF()<<endl;
      // f1->GetChisquare()<<"/"<<f1->GetNDF()<<std::endl;
      //hChi->Fill(startFitTime,f1->GetChisquare()/f1->GetNDF());
      //hTriplet->Fill(startFitTime,1/f1->GetParameter(1));
      if(fabs(f1->GetChisquare()/f1->GetNDF() - 1) < fabs(bestChi-1) ){
        if(t_long_error/(1/f1->GetParameter(3)) > 1.) continue;
        bestChi = f1->GetChisquare()/f1->GetNDF();
        bestTrip = 1/f1->GetParameter(3);
        bestSingle = 1/f1->GetParameter(1);
        bestErrTrip = t_long_error;
        bestErrSinglet = t_short_error;
        bestStart = startFitTime;
      }
    }
    */

    Double_t startFitTime = hArrivalTime->GetBinCenter(hArrivalTime->GetMaximumBin());
    if(z>1000)
      f1 = new TF1("myfit","([0]/[1])*exp(-x/[1]) + ([2]/[3])*exp(-x/[3])",startFitTime,maxTime);
    else
      f1 = new TF1("myfit","([0]/[1])*exp(-x/[1]) + ([2]/[3])*exp(-x/[3])",startFitTime,maxTime);
    f1->SetParLimits(0,0,1e-10);
    f1->SetParameter(1,1e-6);
    f1->SetParLimits(1,4e-7,2e-6);
    f1->SetParLimits(2,0,1e-10);
    f1->SetParameter(3,1e-7);
    f1->SetParLimits(3,5e-9,5e-7);
    f1->SetParLimits(4,0,1);

    //hArrivalTime->Fit("myfit","Q","",0e-9,maxTime-1000e-9);
    cout<<"startFitTime "<<startFitTime<<", stopTime "<<maxTime-1000e-9<<endl;
    hArrivalTime->Fit("myfit","Q","",startFitTime,maxTime-1000e-9);
    int i =1;
    do{
      cout<<"Refitting...old Triplet is "<<f1->GetParameter(1);
      hArrivalTime->Fit("myfit","Q","",startFitTime + 10e-9*i,maxTime-1000e-9*(i+1));
      cout<<"...new Triplet is "<<f1->GetParameter(1)<<", Chi/ndf = "<<f1->GetChisquare()/f1->GetNDF()<<", err "<<f1->GetParError(1)/f1->GetParameter(1)<<endl;
      i++;
      if(i> 10){
        cout<<"warning...fit did not get below chiSquare/ndf = 2"<<endl;
        break;
      }
    }
    while(f1->GetChisquare()/f1->GetNDF() > 2 || f1->GetParError(1)/f1->GetParameter(1) > 0.1);
    if(z >= 110 && z <= 140) hArrivalTime->Fit("myfit","Q","",50e-9,2000e-9);
    bestTrip = f1->GetParameter(1);
    bestErrTrip = f1->GetParError(1);
    bestSinglet = f1->GetParameter(3);
    bestErrSinglet = f1->GetParError(3);
    gStyle->SetOptFit(1111); 
    triplet.push_back(bestTrip);
    //tripletErr.push_back(bestErrTrip+bestTrip*0.0272821);
    tripletErr.push_back(bestErrTrip);
    cout<<"Triplet = "<<bestTrip<<" +/- "<<bestErrTrip<<", Singlet = "<<bestSinglet<<"+/-"<<bestErrSinglet<<" Chi/ndf = "<<f1->GetChisquare()/f1->GetNDF()<<" ratio = "<<f1->GetParameter(0)/f1->GetParameter(2)<<endl;
    A_s.push_back(f1->GetParameter(0));
    A_l.push_back(f1->GetParameter(2));
    //hHits->GetYaxis()->SetTitle("Arrival time weighted with Peak Height");
    //hHits->GetYaxis()->SetTitle("integrated charge");
    //hArrivalTime->Draw();
    //c->SetLogy();
    hHits->Write();
    /*
    hPromptTime->Fit("landau","Q","");
    promptTime.push_back(hPromptTime->GetFunction("landau")->GetParameter(1));
    promptTimeError.push_back(hPromptTime->GetFunction("landau")->GetParameter(2));
    hLongTime->Fit("landau","Q","");
    longTime.push_back(hLongTime->GetFunction("landau")->GetParameter(1));
    longTimeError.push_back(hLongTime->GetFunction("landau")->GetParameter(2));
    hPromptTime->Write();
    hLongTime->Write();
    */
    /*
    //    hError->Draw();
    for(int i = 0 ; i < hHits->GetNbinsX();i++){
      if(deltaQ[i] == 0) continue;
      Double_t time = hHits->GetBinCenter(i);
      Double_t func = f1->Eval(time);
    }
    */

    hSinglePhoton->Write();
    hNhits->Write();
    //hChi->Write();
    //hTriplet->Write();
    //hTotalCharge->Scale(1./hTotalCharge->GetEntries());

    /*
    Double_t totalChargeMax = hTotalCharge->GetBinCenter(hTotalCharge->GetMaximumBin());
    Double_t totalChargeMaxErr = hTotalCharge->GetBinError(hTotalCharge->GetMaximumBin());
    TF1 *hGausCharge = new TF1("GausCharge","landau(0)",0,500);
    hGausCharge->SetParLimits(0,0,1000);
    hGausCharge->SetParLimits(1,1,500);
    hGausCharge->SetParLimits(3,0,1000);
    hTotalCharge->Fit(hGausCharge,"Q","",0,500);
    Double_t sigmaCharge = hGausCharge->GetParameter(1);
    Double_t sigmaErr = hGausCharge->GetParameter(2);
    //cout<<" SigmaCharge "<< sigmaCharge<<"+/-"<<sigmaErr<<", ChargeMax "<<totalChargeMax<<"+/-"<<totalChargeMaxErr<<endl;
    if(0 > sigmaCharge) {
      sigmaCharge = totalChargeMax;
      sigmaErr = totalChargeMaxErr;
    }
    */
    Double_t sigmaCharge = 0;
    Double_t sigmaErr = 0;
    Double_t promptErr = 0;
    Double_t intermediateTimeErr = 0;
    Double_t longChargeErr = 0;
    /*
    // Setting fit range and start values
    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    fr[0]=0.3*hTotalCharge->GetMean();
    fr[1]=3.0*hTotalCharge->GetMean();

    //pllo[0]=0.5; pllo[1]=5.0; pllo[2]=1.0; pllo[3]=0.4;
    //plhi[0]=5.0; plhi[1]=50.0; plhi[2]=1000000.0; plhi[3]=5.0;
    pllo[0]=5.; pllo[1]=5.0; pllo[2]=1.0; pllo[3]=0.4;
    plhi[0]=100.0; plhi[1]=500.0; plhi[2]=1000000.0; plhi[3]=100.0;
    //"Width","MP","Area","GSigma")
    sv[0]=10; sv[1]=150.0; sv[2]=10.0; sv[3]=3.0;
    //sv[0]=1.8; sv[1]=20.0; sv[2]=50000.0; sv[3]=3.0;


    Double_t chisqr;
    Int_t    ndf;
    //totalCharge fit
    TF1 *fitsnr = langaufit(hTotalCharge,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    
    Double_t SNRPeak, SNRFWHM;
    langaupro(fp,SNRPeak,SNRFWHM);

    TF1 *fLandau = new TF1("LandauCharge","landau(0)",fr[0],fr[1]);
    sigmaCharge = SNRPeak;//fitsnr->GetParameter(1);
    sigmaErr = fitsnr->GetParameter(0);
    if(chisqr/ndf > 5) {
      hTotalCharge->Fit(fLandau,"Q","",fr[0],fr[1]);
      sigmaCharge = fLandau->GetParameter(1);
      //sigmaErr = fLandau->GetParError(1);
      sigmaErr = fLandau->GetParameter(2);
    }
    totalCharge.push_back(sigmaCharge);
    totalChargeErr.push_back(sigmaErr);


    //promptTime fit
    fr[0]=0.3*hPromptTime->GetMean();
    fr[1]=3.0*hPromptTime->GetMean();
    //"Width","MP","Area","GSigma")
    pllo[0]=1.; pllo[1]=1.0; pllo[2]=1.0; pllo[3]=0.4;
    plhi[0]=50.0; plhi[1]=200.0; plhi[2]=1000000.0; plhi[3]=100.0;
    sv[0]=10; sv[1]=50.0; sv[2]=10.0; sv[3]=3.0;
    fitsnr = langaufit(hPromptTime,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    langaupro(fp,SNRPeak,SNRFWHM);
    promptCharge = fitsnr->GetParameter(1);
    promptErr = fitsnr->GetParError(1);
    promptErr = fitsnr->GetParameter(0);
    promptTime.push_back(SNRPeak);//promptCharge);
    promptTimeError.push_back(promptErr);


    //intermediate fit
    fitsnr = langaufit(hIntermediateTime,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    intermediateCharge = fitsnr->GetParameter(1);
    intermediateTimeErr = fitsnr->GetParError(1);
    intermediateTime.push_back(fitsnr->GetParameter(1));
    intermediateTimeError.push_back(fitsnr->GetParError(1));


    //longTime fit
    fr[0]=0.3*hLongTime->GetMean();
    fr[1]=3.0*hLongTime->GetMean();
    //"Width","MP","Area","GSigma")
    pllo[0]=1.; pllo[1]=10.0; pllo[2]=1.0; pllo[3]=0.4;
    plhi[0]=100.0; plhi[1]=700.0; plhi[2]=1000000.0; plhi[3]=200.0;
    sv[0]=10; sv[1]=100.0; sv[2]=10.0; sv[3]=3.0;
    fitsnr = langaufit(hLongTime,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    langaupro(fp,SNRPeak,SNRFWHM);

    if(chisqr/ndf > 5){
      if(z<4000){
        hLongTime->Fit("gaus","Q","");
        longTime.push_back(hLongTime->GetFunction("gaus")->GetParameter(1));
        //longTimeError.push_back(fitsnr->GetParError(1));
        longTimeError.push_back(hLongTime->GetFunction("gaus")->GetParameter(2));
      }
      else if( z >=4000 && z<5000){
        hLongTime->Fit(fLandau,"Q","",fr[0],fr[1]);
        longTime.push_back(fLandau->GetParameter(1));
        //sigmaErr = fLandau->GetParError(1);
        longTimeError.push_back(fLandau->GetParameter(2));
      }
    }
    else{
      longTime.push_back(SNRPeak);//fitsnr->GetParameter(1));
      //longTimeError.push_back(fitsnr->GetParError(1));
      longTimeError.push_back(fitsnr->GetParameter(0) );
    }
    */
    sigmaCharge  = hArrivalTime->IntegralAndError(1,hArrivalTime->GetNbinsX(),sigmaErr,"");
    sigmaCharge /= NEvents*singlePE;
    //totalCharge.push_back(sigmaCharge);
    sigmaErr/= NEvents*singlePE;
    totalChargeErr.push_back(sigmaErr);
    promptCharge = hArrivalTime->IntegralAndError(1,hArrivalTime->FindBin(50e-9),promptErr,"");
    promptCharge/= NEvents*singlePE;
    //promptTime.push_back(promptCharge);
    promptErr/= NEvents*singlePE;
    promptTimeError.push_back(promptErr);
    longCharge   = hArrivalTime->IntegralAndError(hArrivalTime->FindBin(50e-9),hArrivalTime->GetNbinsX(),longChargeErr,"");
    longCharge/= NEvents*singlePE;
    //longTime.push_back(longCharge);
    longChargeErr/= NEvents*singlePE;
    longTimeError.push_back(longChargeErr);
    
    totalCharge.push_back(IntegralRun/(NEvents*singlePE));
    promptTime.push_back(promptChargeRun/(NEvents*singlePE));
    longTime.push_back(longChargeRun/(NEvents*singlePE));


    //cout<<"Total Charge "<<sigmaCharge<<"+/-"<<sigmaErr<<", Singlet "<<promptCharge<<"+/-"<<promptErr<<", Intermediate "<<intermediateCharge<<"+/-"<<intermediateTimeErr<<", Triplet "<<longCharge<<"+/-"<<longChargeErr<<endl;

    cout<<"Total Charge "<<IntegralRun/(NEvents*singlePE)<<"+/-"<<sigmaErr<<", Singlet "<<promptChargeRun/(NEvents*singlePE)<<"+/-"<<promptErr<<", Triplet "<<longChargeRun/(NEvents*singlePE)<<"+/-"<<longChargeErr<<endl;
    
    hTotalCharge->Write();
    hPromptTime->Write();
    hIntermediateTime->Write();
    hLongTime->Write();
    hArrivalTimeUnweighted->Write();
    hArrivalTime->Write();
  }

  std::sort(time.begin(),time.end());
  TH1D * hTripletvTime = new TH1D("Triplet_versus_Time","Triplet_versus_Time",time.size()-1,&(time[0]));
  TH1D * hChargevTime = new TH1D("Charge_versus_Time","Total_Photo_Electron_versus_Time",time.size()-1,&(time[0]));
  TH1D * hPromtvTime = new TH1D("Prompt_versus_Time","Promt_versus_Time",time.size()-1,&(time[0]));
  //TH1D * hIntermediatevTime = new TH1D("Intermediate_versus_Time","Intermediate_versus_Time",time.size()-1,&(time[0]));
  TH1D * hLongvTime = new TH1D("Long_versus_Time","Long_versus_Time",time.size()-1,&(time[0]));
  TH1D * hSinglePEvTime = new TH1D("SinglePE_versus_Time","SinglePE_versus_Time",time.size()-1,&(time[0]));
  Double_t tripMean = 0, tripErr = 0;
  Double_t SPEMean = 0, SPEErr = 0;
  for(int i = 0; i < time.size(); i++){
    if((i-1)%2 == 0){
      tripMean += triplet[(i-1)/2.];
      SPEMean += singleCharge[(i-1)/2.];
      SPEErr +=  singleChargeError[(i-1)/2];
      hTripletvTime->SetBinContent(i,triplet[(i-1)/2.]);
      hTripletvTime->SetBinError(i,tripletErr[(i-1)/2.]);
      hChargevTime->SetBinContent(i,totalCharge[(i-1)/2]);
      hChargevTime->SetBinError(i,totalChargeErr[(i-1)/2]);
      hPromtvTime->SetBinContent(i,promptTime[(i-1)/2]);
      hPromtvTime->SetBinError(i,promptTimeError[(i-1)/2]);
      //hIntermediatevTime->SetBinContent(i,intermediateTime[(i-1)/2]);
      //hIntermediatevTime->SetBinError(i,intermediateTimeError[(i-1)/2]);
      hLongvTime->SetBinContent(i,longTime[(i-1)/2]);
      hLongvTime->SetBinError(i,longTimeError[(i-1)/2]);
      hSinglePEvTime->SetBinContent(i,singleCharge[(i-1)/2]);
      hSinglePEvTime->SetBinError(i,singleChargeError[(i-1)/2]);
    }
  }
  tripMean /= triplet.size();
  SPEMean /= singleCharge.size();
  SPEErr /= singleChargeError.size();
  std::vector<Double_t> timeMod, tripletMod, tripletModErr;
  Double_t startTime = 0,stopTime = 0;
  for(int i = 0; i < runNumber.size(); i++){
    if(runNumber[i] == 161) startTime = time[2*i];
    if(runNumber[i] == 166) stopTime = time[2*i];
    if(runNumber[i] < 161){ 
      timeMod.push_back(time[2*i]);
      timeMod.push_back(time[2*i + 1]);
      tripletMod.push_back(triplet[i]);
      tripletModErr.push_back(tripletErr[i]);
      //cout<<"Erasing "<<runNumber[i]<<" time "<<time[2*i]<<"..."<<time[2*i+1]<<" triplet "<<triplet[i]<<endl;
    }
    else if (runNumber[i] > 165){
      timeMod.push_back(time[2*i] - (stopTime-startTime));
      timeMod.push_back(time[2*i + 1] - (stopTime-startTime));
      //cout<<"time "<<time[2*i]<<", time + 1 "<<time[2*i+1]<<", delta time "<<(stopTime-startTime)<<endl;
      tripletMod.push_back(triplet[i]);
      tripletModErr.push_back(tripletErr[i]);
    }
  }
  TH1D * hTripletvTimeModified = new TH1D("Triplet_versus_TimeModified","Triplet_versus_TimeModified",timeMod.size()-1,&(timeMod[0]));
  for(int i = 0; i < timeMod.size(); i++){
    if((i-1)%2 == 0){
      //cout<<timeMod[i-1]<<"..."<<timeMod[i]<<" tripelt "<<tripletMod[(i-1)/2]<<endl;
      hTripletvTimeModified->SetBinContent(i,tripletMod[(i-1)/2.]);
      hTripletvTimeModified->SetBinError(i,tripletModErr[(i-1)/2.]);
    }
  }

  TF1 *fFilter = new TF1("Filter Rate","[0]/(1+[1]*[0]*100*exp(-(x)/[2]))",0,30);
  fFilter->SetParameter(0,1.6e-6);
  fFilter->SetParLimits(0,1.e-6,2e-6);
  fFilter->SetParameter(1,40);
  fFilter->SetParLimits(1,10,1e15);
  fFilter->SetParameter(2,3);
  fFilter->SetParLimits(2,0,100);

  hTripletvTime->SetMarkerStyle(4);
  hTripletvTime->SetMarkerColor(4);
  hTripletvTime->Draw();
  hTripletvTimeModified->Fit(fFilter,"Q","",0,17);
  hTripletvTimeModified->SetMarkerColor(2);
  hTripletvTimeModified->SetMarkerStyle(5);
  if(irunStart == 100)
    hTripletvTimeModified->Draw("same");
  TCanvas * c2 = new TCanvas("c2","c2");
  c2->cd();
  hChargevTime->Draw();
  hPromtvTime->Draw("same");
  hPromtvTime->SetLineColor(2);
  hLongvTime->Draw("same");
  hLongvTime->SetLineColor(3);
  //hIntermediatevTime->Draw("same");
  //hIntermediatevTime->SetLineColor(4);
  
  std::vector<Double_t> concentration;
  for(int i = 0; i < tripletMod.size();i++){
    Double_t val = (100e-6*exp(-(timeMod[2*i]-10.5)/fFilter->GetParameter(2)));
    if(timeMod[2*i]-10.5 > 0){
      concentration.push_back( val );    
      //cout<<"time "<<timeMod[2*i]<<", N2 "<<concentration[i]<<" triplet "<<tripletMod[i]<<endl;
    }
    else{ 
      tripletMod.erase(tripletMod.begin()+i);
    }
  }
  std::sort(concentration.begin(),concentration.end());
  TH1D * hTripletModvConcentration = new TH1D("Triplet_versus_Nitrogen_Concentration","Triplet_versus_Nitrogen_Concentration",concentration.size()-1,&(concentration[0]));
  for(int i = 0; i < tripletMod.size();i++){
    hTripletModvConcentration->SetBinContent(i,tripletMod[tripletMod.size() -i -1]/1.622e-6);
    hTripletModvConcentration->SetBinError(i,tripletModErr[tripletMod.size() -i -1]/1.622e-6);
  }
  
  TCanvas * cA = new TCanvas("cA","cA");
  cA->cd();
  if(irunStart == 100)
    hTripletModvConcentration->Draw(); 

  TH1D * hTripletvCharge = new TH1D("Tripet_versus_Photo_Electron_Yield","Triplet_v_P.E.",1000,400e-9,1600e-9);
  TH1D * hTripletvRunNumber = new TH1D("Tripet_versus_RunNumber","Triplet_v_RunNumber.",runNumber.size()-1,&(runNumber[0]));
  TH1D * hTotalChargevRunNumber = new TH1D("Total_Charge_versus_RunNumber","Total_Charge_versus_RunNumber",runNumber.size()-1,&(runNumber[0]));
  TH1D * hTripletChargevRunNumber = new TH1D("Triplet_Charge_versus_RunNumber","Triplet_Charge_versus_RunNumber",runNumber.size()-1,&(runNumber[0]));
  TH1D * hSingletChargevRunNumber = new TH1D("Singlet_Charge_versus_RunNumber","Singlet_Charge_versus_RunNumber",runNumber.size()-1,&(runNumber[0]));
  Double_t survial = 0;
  for(int i = 0; i < triplet.size(); i++){
    Int_t tripletBin = hTripletvCharge->FindBin(triplet[i]);
    hTripletvCharge->SetBinContent(tripletBin,totalCharge[i]);
    hTripletvCharge->SetBinError(tripletBin,totalChargeErr[i]);
    hTripletvRunNumber->SetBinContent(i,triplet[i]);
    hTripletvRunNumber->SetBinError(i,tripletErr[i]);
    hTotalChargevRunNumber->SetBinContent(i,totalCharge[i]);
    hTripletChargevRunNumber->SetBinContent(i,longTime[i]);
    hSingletChargevRunNumber->SetBinContent(i,promptTime[i]);
    tripErr += (tripMean-triplet[i])*(tripMean-triplet[i]);
    survial += A_s[i] +A_l[i];
  }
  tripErr = std::sqrt(tripErr/triplet.size());
  cout<<"mean triplet "<<tripMean<<", error "<<tripErr<<", fraction "<<tripErr/tripMean<<endl;
  cout<<"mean amplitude "<<survial/triplet.size()<<endl;
  cout<<"mean SPE "<<SPEMean<<", error "<<SPEErr<<endl;
  TCanvas * c3 = new TCanvas("c3","c3");
  c3->cd();
  //hTripletvCharge->Draw();
  hTotalChargevRunNumber->Draw();
  hTripletChargevRunNumber->Draw("same");
  hSingletChargevRunNumber->Draw("same");


  TCanvas * c4 = new TCanvas("c4","c4");
  c4->cd();
  hTripletvRunNumber->Draw();

  
outfile->Write();
}
