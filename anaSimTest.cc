#include "anaSim.hh"
//

anaSim::anaSim(Int_t z,Int_t nDer){
   //time 12:36:26 133626
  ///date 24/12/1997 19971224
  //TString outFileName = TString("SimAnaResults.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root");
  TString fileDir = "/home/nmcfadde/RooT/PMT/rootData/";
  //for(int z = 183; z <= 184; z++){
  //if(z == 10) continue;
  //if(z == 26) continue;
  //TString fileName = "simEvents_NEvents_1000_nPhotons_50.root";
  //TString fileName = TString("baconRun_10kBins_400ns_10mV_")+to_string(z)+TString(".root");
  //TString fileName = TString("baconRun_10kBins_1us_10mV_")+to_string(z)+TString(".root");
  TString fileName;
  if(z >= 1000 && z <2000)
    fileName = TString("baconRun_10kBins_1us_20mV_muon_")+to_string(z)+TString(".root");
  else if(z>=2000 && z < 4000)
    fileName = TString("baconRun_10kBins_1us_20mV_div_-30mV_thresh_") + to_string(z) +TString(".root");
  else if(z>=4000 && z < 6000)
    fileName = TString("baconRun_10kBins_1us_20mV_div_-6mV_thresh_") + to_string(z) +TString(".root");
  else if(z>=10000)
    fileName = TString("baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_")+ to_string(z) +TString(".root");
  else if (z >= 100 && z <=272)
    fileName = TString("baconRun_10kBins_1us_10mV_")+ to_string(z) +TString(".root");
  else
    fileName = "simEvents_20190418_10000.root";//"simEvents_20190308_10000.root";//simEvents_20190226_10000.root";//"simEvents_20190227_100.root";//"simEvents_20190226_10000.root";
  TString outFileName = TString("/home/nmcfadde/RooT/PMT/processedData/SimAnaResults.")+fileName;
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());

  TNtuple *ntuplePulse = new TNtuple("ntuplePulse","ntuplePulse","ientry:pmtNum:nhits:charge:startTime:peakWidth:T0:V0:vMax:vMaxTime:Sdev:baseline");
  //TNtuple *ntupleEvent = new TNtuple("ntupleEvent","ntupleEvent","irun:ientry:pmtNum:Sdev:baseline:integral:deltaT:deltaV:vMax:tMax");
  TNtuple *ntupleEvent = new TNtuple("ntupleEvent","ntupleEvent","irun:ientry:Sdev:baseline:integral:deltaT:tMax:vMax");
  TNtuple *ntupleSim = new TNtuple("ntupleSim","ntupleSim","irun:ientry:nMatch:deltaT:vMaxFound:truthVMax:charge");
  TNtuple *ntupleSimEvent = new TNtuple("ntupleSimEvent","ntupleSimEvent","irun:ientry:TotalCharge:nFound:nMiss:nNoise:truthN");  

//  TString fileName = "pmtEvents.100.100photons.root";
  
  TFile *fin = new TFile(fileDir+fileName,"readonly");
  
  if(fin->IsZombie()) {
    printf(" looking for file %s%s\n",fileDir.Data(),fileName.Data());
    return;
  }
  else 
    printf(" looking for file %s%s\n",fileDir.Data(),fileName.Data());

  // get pmtTree from file 
  TTree * pmtTree = new TTree();
  fin->GetObject("pmtTree",pmtTree);
  Long64_t nentries = pmtTree->GetEntries();
  cout << " number of entries is " << nentries << endl;

  // set up memory for reading
  simEvent = new TPmtSimulation();
  pmtEvent = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtSimulation", &simEvent);
  cout<<"Simulation branch set"<<endl;
  pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);
  cout<<"Event branch set"<<endl;

  TDatime time;

  TNamed * eventInfo;

  fin->GetObject("eventInfo",eventInfo);
  if(eventInfo != NULL)
    cout<<"eventInfo is "<<eventInfo->GetTitle()<<endl;
  //switch to output file
  outFile->cd();
  
  TH1D *hTruthPulsesMissed = new TH1D("MissedPulses","MissedPulses",10000,0,10e-6);
  TH1D *hNoisePulsesFound = new TH1D("NoisePulsesFound","NoisePulsesFound",10000,0,10e-6);
  
  for(Long64_t ientry = 0; ientry < nentries; ientry++){
    pmtTree->GetEntry(ientry);
    signal.resize(1);derivative.resize(signal.size());integral.resize(derivative.size());
   
    if(z == 107 || z == 108 || z == 109){
      for(int i = 0; i < pmtEvent->time.size();i++){
        pmtEvent->time[i] = 1e-9*i*10000./ pmtEvent->time.size();
      }
    }
    //simEvent is not intialized if it is NULL
    if(z == 0 && simEvent->Nphotons == 0) continue;
    std::vector<Double_t> simTimes = simEvent->startTime;
    std::vector<Double_t> simCharge = simEvent->q;
    std::sort(simTimes.begin(),simTimes.end());//,std::greater<Double_t>());
    //cout<<"starting event"<<endl;
    signal[0] = pmtEvent->volt1;
    NEvents = pmtEvent->volt1.size();
    deltaT = (pmtEvent->time[NEvents -1] - pmtEvent->time[0])/(NEvents-1);
    if(ientry == 0) cout<<"time bin size is "<<deltaT<<endl;  
    peakFindingDebug = false;
    ///*
    if( ientry == 5){
      if(eventInfo != NULL)
        eventInfo->Write();
      break;
    }
    peakFindingDebug = true;
    Int_t skipBin = 4;
    NHistograms = skipBin;
    //if(skipBin < ientry ) break;
    if( skipBin != ientry) continue;
    hSum[0] = new TH1D(TString("Sum_")+to_string(0),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
    hSumFresh[0] = new TH1D(TString("SumVirgin_")+to_string(0),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
    //*/


    TH1D* hPeakFinding[1];
    Int_t nInt = 1;
    Double_t sum = 0;
    for(int j = 0; j < signal.size();j++){
      hSignal[j] = new TH1D(TString("Signal_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      hIntegral[j] = new TH1D(TString("Integral_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      if(ientry == 0){
        hSum[j] = new TH1D(TString("SumModified_")+to_string(j),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
        hSumFresh[j] = new TH1D(TString("SumVirgin_")+to_string(j),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      }
      hPeakFinding[j] = new TH1D(TString("PeakFinding")+to_string(ientry),TString("PeakFinding")+to_string(ientry),NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);

      std::vector<Double_t> integral;
      for(int i = 0; i < signal[j].size(); i++){
        //signal[j][i] = signal[j][i] - fPoly[j]->Eval(pmtEvent->time[i]);
        //hSum[j]->SetBinContent(i+1,hSum[j]->GetBinContent(i+1)-signal[j][i]);
        sum += signal[j][i];
        integral.push_back(sum/nInt);
      }
      //signal[j] = TrapFilter(integral,40,5,j,ientry);
      /*
      if(signal[j].size() == 10000)
        signal[j] = Derivative(integral,nInt);
      else nInt = 1;
        signal[j] = Derivative(integral,nInt);
      */
      //signal[j] = DownSampler(signal[j],5);
      if(signal[j].size() != 10000) nDer = 2;
      derivative[j] = Derivative(signal[j],nDer);
      hDerivative[j] = new TH1D(TString("derivative_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
    }
    sum = 0;
    //cout<<"finished first loop"<<endl;
    for(int j = 0; j < derivative.size(); j++){
      for(int i = 0; i < derivative[j].size(); i++){
        //signal[j][i] /= (double) nInt;
        sum += signal[j][i];
        hSignal[j]->SetBinContent(i,signal[j][i]);
        hPeakFinding[j]->SetBinContent(i+1,0);
        hDerivative[j]->SetBinContent(i+1,derivative[j][i]);
        hIntegral[j]->SetBinContent(i+1,sum*deltaT);
      }
    }
    if(peakFindingDebug) hSignal[0]->Draw();
    //cout<<"setting derivative"<<endl;
    //Peak Finding
    std::vector< std::vector<Int_t> > peakTime;peakTime.resize(signal.size());
    std::vector< std::vector<Double_t> > peakTimeDoubleStart;peakTimeDoubleStart.resize(signal.size());
    std::vector< std::vector<Double_t> > peakVoltsStart;peakVoltsStart.resize(signal.size());
    std::vector< std::vector<Double_t> > peakTimeDoubleStop;peakTimeDoubleStop.resize(signal.size());
    std::vector< std::vector<Double_t> > peakVoltsStop;peakVoltsStop.resize(signal.size());
    Int_t maxPeakWidth = 0;
    for(int i = 0; i < signal.size(); i++){
      //cout<<"peakFinding"<<" derSize "<<derivative[i].size()<<" "<<ientry<<endl;
      std::vector<Int_t> pTime = PeakFinding(signal[i],derivative[i],i,ientry);
      if(pTime.size() == 0) continue;
      //cout<<"Culling the small and weak"<<endl;
      //Cull the small pulses
      for(int j = 0; j< pTime.size() -1 ;j +=2){
        Int_t startBin = pTime[j],stopBin = pTime[j+1];
        Double_t peakWidth = pmtEvent->time[stopBin] - pmtEvent->time[startBin] ;
        if( (stopBin-startBin) > maxPeakWidth) maxPeakWidth = (stopBin-startBin);
        if(signal[j].size() == 500) minPeakWidth = deltaT;//One bin 
        
        if(peakFindingDebug) cout<<"start "<<startBin<<" "<<pmtEvent->time[startBin]<<", stop "<<stopBin<<" "<< pmtEvent->time[stopBin]<<endl;
        //cout<<"\t peakWidth "<<peakWidth<<endl;
        if(peakWidth > minPeakWidth && peakWidth < maxPeakWidth) {
          peakTime[i].push_back(pTime[j]);
          peakVoltsStart[i].push_back(signal[i][pTime[j]]);
          peakTimeDoubleStart[i].push_back(pmtEvent->time[startBin]);
          //cout<<signal[i][pTime[j]]<<" "<< pTime[j]<<endl;
          peakTime[i].push_back(pTime[j+1]);
          peakVoltsStop[i].push_back(signal[i][pTime[j+1]]);
          peakTimeDoubleStop[i].push_back(pmtEvent->time[stopBin]);
          //cout<<"\t"<<signal[i][pTime[j+1]]<<" "<< pTime[j+1]<<endl;
        }
      }
    }
    TGraph *tFindingStart = new TGraph(peakTime[0].size(),&(peakTimeDoubleStart[0][0]),&(peakVoltsStart[0][0]));
    if(peakFindingDebug){
      tFindingStart->SetMarkerColor(3);
      tFindingStart->Draw("same*");
    }
    TGraph *tFindingStop = new TGraph(peakTime[0].size(),&(peakTimeDoubleStop[0][0]),&(peakVoltsStop[0][0]));
    if(peakFindingDebug){
      tFindingStop->SetMarkerStyle(2);
      tFindingStop->SetMarkerColor(2);
      tFindingStop->Draw("same*");
    }
    //cout<<"finished peakFinding"<<endl;
    //2N + 1 moving window for WMA
    Int_t movingWindow = 300;
    if(signal[0].size() == 500)
      movingWindow = 10;
    std::vector<Double_t> BaSeLiNe = BaselineWMA(signal[0],peakTime[0],ientry,movingWindow);
    /*
    for(int i = 0; i < signal[0].size();i++){
      if(peakTime[0].empty()) continue;
      hSum[0]->SetBinContent(i+1,hSum[0]->GetBinContent(i+1)-signal[0][i]+BaSeLiNe[i]);
      hSumFresh[0]->SetBinContent(i+1,hSumFresh[0]->GetBinContent(i+1)-signal[0][i]);
    }
    */
    //filling ntuple
    for(int j = 0; j < signal.size(); j++){
      Double_t vMaxEvent = 0,tMaxEvent = 0;
      TString sTitle = TString("");
      Double_t baseline = CalculateMean(signal[j],0,signal[j].size()*0.05);//500);
      Double_t Sdev =  CalculateSdev(signal[j],0,signal[j].size()*0.05,baseline);
      //Fill event information
      //ntupleEvent->Fill(z,ientry,j,Sdev,baseline,-sum,deltaT,peakTime.size()/2);
      
      if(peakTime[j].empty()) continue;
      Double_t T0=0,V0 = 0;
      Double_t nFound = 0,nMiss = 0,nNoise = 0,totalCharge = 0;
      //cout<<"entry "<<ientry<<endl;
      for(int i = 0; i < peakTime[j].size()-1; i += 2){
        Double_t vMax = 0,charge = 0,vMaxTime = 0,startTime,stopTime,peakWidth;
        Int_t startBin = peakTime[j][i],stopBin = peakTime[j][i+1];
        //Int_t startBin = hSignal[j]->FindBin(400e-9);
        //Int_t stopBin  = hSignal[j]->FindBin(500e-9);
        startTime = pmtEvent->time[startBin];
        stopTime = pmtEvent->time[stopBin];
        peakWidth = stopTime - startTime;
        for(int k = startBin; k < stopBin; k++){
          //if(signal[j][k]-BaSeLiNe[k] > 0) continue;
          charge += signal[j][k]-BaSeLiNe[k];
          hSum[0]->SetBinContent(k+1,hSum[0]->GetBinContent(k+1)-signal[j][k]+BaSeLiNe[i]);
          hSumFresh[0]->SetBinContent(k+1,hSumFresh[0]->GetBinContent(k+1)-signal[j][k]);
          if(std::fabs(signal[j][k]-BaSeLiNe[k]) > std::fabs(vMax) ){
            vMax = signal[j][k]-BaSeLiNe[k];
            vMaxTime = pmtEvent->time[k];
          }
          hPeakFinding[j]->SetBinContent(k,signal[j][k]-BaSeLiNe[k]);
        }
        if(i == 0) {
          T0 = startTime;
          V0 = vMax;
        }
        ntuplePulse->Fill(ientry,j,peakTime[j].size()/2,-charge*deltaT,startTime,peakWidth,T0,V0,-vMax,vMaxTime,0,0);
        //sTitle += TString("Charge_")+to_string(-1e9*deltaT*charge)+TString("_start/stop_")+to_string(1.e6*startTime)+TString("/")+to_string(1.e6*stopTime)+TString("_vMax_")+to_string(-vMax)+TString("_");
        if(std::fabs(vMax) > std::fabs(vMaxEvent) ){
          vMaxEvent = vMax;
          tMaxEvent = vMaxTime;
        }
        Double_t simTimeWindow = 10e-9,deltaSimTime = 9999,truthTime = 0,truthCharge = 0.05,truthVMax = 0;
        Double_t nMatch = 0;
        for(int i = 0 ; i < simTimes.size();i++){
            if(std::fabs(vMaxTime-simTimes[i]) < simTimeWindow){
              if(deltaSimTime > vMaxTime-simTimes[i]) deltaSimTime = vMaxTime-simTimes[i];
              //deltaSimTime = vMaxTime-simTimes[i];
              truthTime = simTimes[i];
              truthCharge = simCharge[i];
              Int_t vBin = hSignal[j]->FindBin(truthTime);
              nMatch++;
              simTimes.erase(simTimes.begin()+i);
              i--;
              if(peakFindingDebug){
                cout<<"charge "<<charge/-truthCharge<<", peakTime "<<vMaxTime<<", truth Time "<<truthTime<<", nMatch "<<nMatch<<endl;
              }
            }
        }
        totalCharge += charge;
        if(nMatch > 0) nFound+=nMatch;
        else{
          nNoise++;
          hNoisePulsesFound->Fill(vMaxTime);
        }
        ntupleSim->Fill(z,ientry,nMatch,deltaSimTime,vMax,truthVMax,charge/-truthCharge);//,(Double_t)nFound/simTimes.size());
      }
      if(nFound < (Double_t) simEvent->startTime.size()) nMiss = (Double_t) simEvent->startTime.size() - nFound;
      if(z == 0){
        ntupleSimEvent->Fill(z,ientry,-totalCharge/((Double_t)simEvent->q[0]* simEvent->startTime.size()),nFound/(Double_t)simEvent->startTime.size(),nMiss,nNoise,(Double_t) simEvent->startTime.size());
        for(int i = 0; i < simTimes.size(); i++) hTruthPulsesMissed->Fill(simTimes[i]);
      }
      //ntupleEvent->Fill(z,ientry,j,Sdev,baseline,-sum,deltaT,peakTime.size()/2,-vMaxEvent,tMaxEvent);
      sort(simEvent->startTime.begin(),simEvent->startTime.end());
      if(z == 0) sTitle += TString("M_") + to_string(nFound)+TString("/G")+to_string(simEvent->startTime.size())+TString("_N/") + to_string(nNoise);
      else sTitle += TString("tMax_")+to_string(tMaxEvent*1e6)+TString("_vMax_")+to_string(vMaxEvent)+TString("_totalCharge_")+to_string(totalCharge);
      ntupleEvent->Fill(z,ientry,Sdev,baseline,-sum,deltaT,tMaxEvent,-vMaxEvent,totalCharge);//-vMaxEvent,tMaxEvent);
      for(int i = 0; i < simTimes.size(); i++) sTitle +=TString("_") + to_string(1e6*simTimes[i]);
      hSignal[j]->SetTitle(sTitle);
    }
//    cout<<"finished filling ntuple"<<endl;
    hPeakFinding[0]->SetLineColor(2);
    if(peakFindingDebug) hPeakFinding[0]->Draw("same");
    //###################//
    //End of Event Clean-Up//
    //###################//
    /*
    TH1 *hm =0;
    TVirtualFFT::SetTransform(0);
    hm = hSignal[0]->FFT(hm, "MAG");
    hm->SetTitle("Magnitude of the 1st transform");
    hm->Draw();
    */
    if(ientry > NHistograms){
      for(int i = 0; i < signal.size(); i++){
        delete hSignal[i];
        delete hDerivative[i];
        delete hIntegral[i];
        delete hPeakFinding[i];
      }
    }
    signal.clear();
    derivative.clear();
    if(ientry == nentries -1 && eventInfo != NULL){
      eventInfo->Write();
      cout<<"event info written "<<endl;
    }
    if(ientry%100 == 0) cout<<"completed "<<ientry<<" events"<<endl;
  }
  //hSumFresh[0]->Fit("expo","","",2e-6,7e-6);
  //hSum[0]->Fit("expo","","",2e-6,7e-6);
  cout<<"root -l "<<outFileName<<endl;
  //ntuplePulse->Write();
  outFile->Write();
  //}
  //cout<<"root -l "<<outFileName<<endl;
  cout<<"you did it! your code ran with crashing"<<endl;
}
/*
std::vector<Double_t> anaSim::Derivative(std::vector<Double_t> sig,Int_t N){
  std::vector<Double_t> diff;
  if(N%2 != 0) N++;
  
  for(int i = 0; i < sig.size(); i++){
    Double_t front= 0, back = 0;
    Double_t nFront = 0, nBack = 0;

    for(int j = i;j<i+N;j++){
      if(j >= sig.size()) continue;
      front += sig[j];
      //cout<<front<<" "<<sig[j]<<" "<<j<<endl;
      nFront++;
    }
    for(int j = i;j>i-N;j--){
      if(j < 0) continue;
      back  += sig[j];
      nBack++;
    }
    if(nFront == 0 || nBack == 0){
      diff.push_back(0);
    }
    //Double subtraction leaving tiny remaining difference
    else if( log(std::fabs(front/nFront-back/nBack))  < -10)
      diff.push_back(0);
    else
      diff.push_back(front/nFront-back/nBack);
  }
  
  return diff;
}
*/

std::vector<Double_t> anaSim::DownSampler(std::vector<Double_t> sig,Int_t N){

  if(N%2 == 0) N+1;
  std::vector<Double_t> downSample;
  cout<<"starting down sampler "<<endl;
  for(int i = 0; i < sig.size(); ){
    Int_t arr [] = {N,i,(Int_t)sig.size()-1-i};
    Int_t window = *std::min_element(arr,arr+3);
    Double_t val = 0;
      for(int j = 1; j < window ; j++){
        val += (sig[i+j] + sig[i-j])/window;
      }
      downSample.push_back(window);
      if( window == 0) i++;
      else i += window;
  }
  cout<<downSample.size()<<endl;
  return downSample;
}

std::vector<Double_t> anaSim::Derivative(std::vector<Double_t> sig,Int_t N){
  std::vector<Double_t> diff;
  //diff.push_back(0);
  Double_t leftSum = 0, rightSum = 0;
  //if(N<2) N = 2;
  for(int i = 0; i < sig.size(); i++){
    Double_t der = 0;
    Int_t arr [] = {N,i,(Int_t)sig.size()-1-i};
    Int_t window = *std::min_element(arr,arr+3);
    if( i <= N && 2*i <= sig.size() - 1){
      leftSum  = leftSum + sig[i - 1];
      rightSum = rightSum - sig[i] + sig[2*i-1] +sig[2*i]; 
    }
    else if (i > N && i + N <= sig.size() - 1){
      leftSum = leftSum + sig[i-1] - sig[i-1-N];
      rightSum = rightSum - sig[i] + sig[i+N];
    }
    else if( i + N > sig.size() -1 && 2*i > sig.size() -1){
      leftSum = leftSum + sig[i-1] - sig[2*i-sig.size()-1]- sig[2*i-sig.size()];
      rightSum = rightSum -sig[i];
    }
    else{ cout<<"you done messed the if else statements"<<endl;}
    diff.push_back((rightSum-leftSum)/(1+window));
    /*
    for(int j = 1; j < window ; j++){
      der += (sig[i+j] - sig[i-j]);
    }
    diff.push_back(der);
    */
  }
  return diff;
}

Double_t anaSim::RMSCalculator(std::vector<Double_t> vec,Int_t ientry){
  Double_t rms = 0,sum = 0;
  TH1D * h1 = new TH1D(TString("rmsHistogram_")+to_string(ientry),TString("rmsHistogram_")+to_string(ientry),500,-BubbleSort(vec)[vec.size() -1],BubbleSort(vec)[vec.size() -1]);
  for(int i = 0; i < vec.size(); i++){
    h1->Fill(vec[i]);
  }
  Int_t centerBin = h1->GetNbinsX()/2;
  h1->SetBinContent(centerBin,sqrt(h1->GetBinContent(centerBin)*(h1->GetBinContent(centerBin+1)+h1->GetBinContent(centerBin-1))/2.));
  for(int i = 0; i < h1->GetNbinsX(); i++){
    h1->SetBinContent(i+1,TMath::Exp(h1->GetBinContent(i+1)/vec.size())-1);
  }
  TF1 * fGaus = new TF1("fit","gaus(0)");
  h1->Fit(fGaus,"Q");
  rms = fGaus->GetParameter(2);
  if(ientry > NHistograms) delete h1;
  return rms;
}
//Weighted Moving Average
std::vector<Double_t> anaSim::BaselineWMA(std::vector<Double_t> sig,std::vector<Int_t> peaks,Int_t ientry,Int_t N){
  std::vector<Double_t> Baseline;
  if(peaks.size() == 0){
    cout<<"BaselineWMA::Fatal, signal size is 0"<<endl;
    return Baseline;
  }
  Int_t ipeak = 1;
  Int_t startPeak = peaks[0];
  Int_t stopPeak = -1;//peaks[ipeak+1];
  std::vector<Double_t> weight;
  Int_t windowCounter = 0;
  //set peak weights
  for(int i = 0; i <sig.size();i++){
    if(i <= startPeak && i >= stopPeak){
      weight.push_back(startPeak - stopPeak -1);
      if(windowCounter > N) N = windowCounter;
      windowCounter = 0;
    }
    else {
      weight.push_back(1.e-6);
      windowCounter++;
     }
    if(i == peaks[ipeak]){
      ipeak +=2;
      if(ipeak -1 == peaks.size()){
        startPeak = sig.size();
      }
      else startPeak = peaks[ipeak-1];
      stopPeak = peaks[ipeak-2];
    }
    //cout<<"Bin "<<i<<", startPeak "<<startPeak<<", stopPeak "<<stopPeak<<", weight "<<weight[i]<<", window "<<window[i]<<", ipeak "<<ipeak<<endl;
  }
  if(N < 100)
    N*=1.5;
//  N += 10;
  //cout<<N<<endl;
  TH1D * f1 = new TH1D(TString("WMA")+to_string(ientry),TString("WMA")+to_string(ientry)+TString("_N_")+to_string(N),sig.size(),0,pmtEvent->time[NEvents-1]);
  /*
  //calculate WMA
  for(int i = 0; i < sig.size(); i++){
    Int_t hiWindow = std::min(i+N,(Int_t)sig.size() -1);
    Int_t loWindow = std::max(0,i-N);
    Double_t top = 0, bot = 0;
    f2->SetBinContent(i+1,weight[i]);
    //cout<<"loWindow "<<loWindow<<", hiWindow "<<hiWindow<<endl;
    for(int j = loWindow;j<hiWindow;j++){
      top +=sig[j]*weight[j]*(1+cos((j-i)*pi/N));
      bot +=weight[j]*(1+cos((j-i)*pi/N));
      //cout<<"\t top "<<top<<", bot "<<bot<<endl;
    }
    Baseline.push_back(top/bot);
    f1->SetBinContent(i+1,top/bot);
  }
  */
  ///*
  //obnoxious way of making things faster
  //i.e. run times does not depend on how big the window, N, is. 
  //It only depends on total number of points
  Double_t K_sw = 0;
  Double_t C_sw = 0;
  Double_t S_sw = 0;
  Double_t K_w = 0;
  Double_t C_w = 0;
  Double_t S_w = 0;
  Double_t kc = cos(pi/N);
  Double_t ks = sin(pi/N);
  Int_t resetBin = 100, lowBinMax = 15;
  Double_t lowBinWeight = N;
  if(sig.size() == 500) resetBin = 20;
  for(int i = 0;i< sig.size();i++){
    Int_t hiWindow = std::min(i+N,(Int_t)sig.size() -1);
    Int_t loWindow = std::max(0,i-N);
    if(weight[i] <= lowBinMax && weight[i] >= 1 && sig[i] > 0 ) weight[i] *= lowBinWeight;
    if(i%resetBin ==  0){
      K_sw = 0;C_sw = 0;S_sw = 0;K_w = 0;C_w = 0;S_w = 0;
      for(int j = loWindow;j<hiWindow;j++){
        if(weight[j] <= lowBinMax && weight[j] >= 1 && sig[j] > 0 ) weight[j] *= lowBinWeight;
        K_sw += sig[j]*weight[j];
        C_sw += sig[j]*weight[j]*cos((j-i)*pi/N);
        S_sw += sig[j]*weight[j]*sin((j-i)*pi/N);
        K_w += weight[j];
        C_w += weight[j]*cos((j-i)*pi/N);
        S_w += weight[j]*sin((j-i)*pi/N);
      }
    }
    else{
      if(i - N <= 0 && i + N <= sig.size() - 1){
        K_sw = K_sw + sig[i+N]*weight[i+N];
        Double_t tempC_sw = C_sw;
        C_sw = kc*C_sw + ks*S_sw - sig[i+N]*weight[i+N];
        S_sw = kc*S_sw - ks*tempC_sw ;
        
        K_w = K_w + weight[i+N];
        Double_t tempC_w = C_w;
        C_w = kc*C_w + ks*S_w - weight[i+N];
        S_w = kc*S_w - ks*tempC_w ;
      }
      else if(i - N > 0 && i + N <= sig.size() - 1){
        K_sw = K_sw - sig[i-1-N]*weight[i-1-N] + sig[i+N]*weight[i+N];
        Double_t tempC_sw = C_sw;
        C_sw = kc*C_sw + ks*S_sw + kc*sig[i-1-N]*weight[i-1-N] - sig[i+N]*weight[i+N];
        S_sw = kc*S_sw - ks*tempC_sw - ks*sig[i-1-N]*weight[i-1-N];
        
        K_w = K_w - weight[i-1-N] + weight[i+N];
        Double_t tempC_w = C_w;
        C_w = kc*C_w + ks*S_w + kc*weight[i-1-N] - weight[i+N];
        S_w = kc*S_w - ks*tempC_w - ks*weight[i-1-N];
      }
      else if(i - N > 0 && i + N > sig.size() - 1){
        K_sw = K_sw - sig[i-1-N]*weight[i-1-N];
        Double_t tempC_sw = C_sw;
        C_sw = kc*C_sw + ks*S_sw + kc*sig[i-1-N]*weight[i-1-N];
        S_sw = kc*S_sw - ks*tempC_sw - ks*sig[i-1-N]*weight[i-1-N];

        K_w = K_w - weight[i-1-N];
        Double_t tempC_w = C_w;
        C_w = kc*C_w + ks*S_w + kc*weight[i-1-N];
        S_w = kc*S_w - ks*tempC_w-ks*weight[i-1-N];
      }
      else if(i - N <=0 && i + N > sig.size() - 1){
        K_sw = K_sw;
        Double_t tempC_sw = C_sw;
        C_sw = kc*C_sw + ks*S_sw;
        S_sw = kc*S_sw - ks*tempC_sw;

        K_w = K_w - weight[i-1-N];
        Double_t tempC_w = C_w;
        C_w = kc*C_w + ks*S_w;
        S_w = kc*S_w - ks*tempC_w;
      }
      else{ cout<<"you messed up the logic on the moving baseline"<<endl;}
    }
    //cout<<"ientry"<<ientry<<" "<<i<<"...K_sw "<<K_sw<<", C_sw "<<C_sw<<", S_sw "<<S_sw<<", K_w "<<K_w<<", C_w "<<C_w<<"..."<< (K_sw+C_sw)/(K_w+C_w)<<endl;
    f1->SetBinContent(i+1,(K_sw+C_sw)/(K_w+C_w));
    Baseline.push_back((K_sw+C_sw)/(K_w+C_w));
  }
  //*/
  f1->SetLineColor(6);
  if(peakFindingDebug) f1->Draw("same");
  //f1->Write();
  if(ientry > NHistograms){
    delete f1;
  }

  return Baseline;
}


std::vector<Double_t> anaSim::Integral(std::vector<Double_t> sig,Int_t pmtNum, Int_t ientry){
  std::vector<Double_t> integ;
  Double_t sum = 0;
  for(int i = 0; i < sig.size(); i++){
    sum += sig[i];
    integ.push_back(sum);
  }
  return integ;
}
std::vector<Int_t> anaSim::PeakFinding(std::vector<Double_t> sig, std::vector<Double_t> diff,Int_t pmtNum,Int_t ientry){
  std::vector<Int_t> peakTime;
  std::vector<Int_t> peakTimeTrimmed;
  //-1 ll, 0, llhh, 1 hh
  std::vector<Int_t>  peakFlag;

  Int_t rmsBin = sig.size()*0.05;//500;
  //Double_t rmsFit = RMSCalculator(diff,ientry);
  Double_t rmsMean = CalculateMean(diff);
  Double_t rmsSdevlo = CalculateSdev(diff,0,rmsBin,rmsMean);
  Double_t rmsSdevhi = CalculateSdev(diff,diff.size()-rmsBin,diff.size()-1,rmsMean);
  Double_t rmsSdevAll = CalculateSdev(diff,0,diff.size()-1,rmsMean);
  Double_t arr [] = {rmsSdevlo,rmsSdevhi,rmsSdevAll};
  Double_t rmsSdev = *std::min_element(arr,arr+3);
  //Double_t rmsSdev = CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  Double_t rms =rmsSdev;// min(rmsFit,rmsSdev);
  Double_t loRMS = 3.5*rms;//CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  Double_t midRMS = 2.5*rms;
  Double_t hiRMS = 3.5*rms;//CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  //cout<<"RMS = "<<rms<<endl;
  hDerivative[pmtNum]->SetTitle(TString("loThreshold_")+to_string(loRMS)+TString("hiThreshold_")+to_string(hiRMS)+TString("_mean_")+to_string(CalculateMean(diff)) );
  hDerivative[pmtNum]->SetLineColor(8);
  Bool_t low1 = false,low2 = false, mid1 = false, mid2 = false,high1 = false, high2 = false, ll = false, hh = false;

  Int_t LowStartBin = -1,LowStopBin = -1,MidStartBin = -1, MidStopBin = -1, HighStartBin = -1,HighStopBin = -1;
  Int_t zeroCrossing = 0;
  for(int i = 1; i < diff.size(); i++){

    if(diff[i] <= -loRMS && diff[i-1] > -loRMS){
      LowStartBin = i;
      low1 =  true;
    }
    if(diff[i] >= -loRMS && diff[i-1] < -loRMS ){
      LowStopBin = i;
      low2 = true;
    }
    if(diff[i] >= midRMS && diff[i-1] < midRMS ){
      MidStartBin = i;
      mid1 = true;
    }
    if(diff[i] <= midRMS && diff[i-1] > midRMS ){
      MidStopBin = i;
      mid2 = true;
    }
    if(diff[i] >= hiRMS && diff[i-1] < hiRMS ){
      HighStartBin = i;
      high1 = true;
    }
    if(diff[i] <= hiRMS && diff[i-1] > hiRMS ){
      HighStopBin = i;
      high2 = true;
    }
    //sometimes noise will trigger a lolo then a long time later a hihi will 
    //be found resulting in a long noise peak that forms a lohi
    //if baseline crosses zero more than onces between lo hi, don't count it
    //resets if a hi && lo were found
    if( ((diff[i] < 0 && diff[i-1] > 0) || (diff[i] > 0 && diff[i-1] < 0))){
      zeroCrossing++;
    }
    if(low1 && low2) {
      if(peakFindingDebug){
        cout<<"low "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[LowStopBin]<<" "<<diff[i]<<" "<<-rms<<endl;
        cout<<"\thi "<<pmtEvent->time[HighStartBin]<<" "<<pmtEvent->time[HighStopBin]<<endl;
      }

      peakTime.push_back(LowStartBin);
      peakTime.push_back(LowStopBin);
      peakFlag.push_back(-1);
      low1 = false; low2 = false;
      if(hh)
        ll = false;
      else{
        ll = true;
        //reset zero crossing if a new lolo peak is found
        zeroCrossing = 0;
      }
    }
    
     
    if(high1 && high2){
      if(peakFindingDebug){
        cout<<"hi "<<pmtEvent->time[HighStartBin]<<" "<<pmtEvent->time[HighStopBin]<<" "<<diff[i]<<" "<<rms<<endl;
        cout<<"\tlow "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[LowStopBin]<<endl;
      }
      peakTime.push_back(HighStartBin);
      peakTime.push_back(HighStopBin);
      peakFlag.push_back(1);
      high1 = false;high2 = false;
      if(ll)
        hh = true;
      else 
        hh = false;
    }
    else if(mid1 && mid2 && ll){
      if(peakFindingDebug){
        cout<<"mid "<<pmtEvent->time[MidStartBin]<<" "<<pmtEvent->time[MidStopBin]<<" "<<diff[i]<<" "<<rms<<endl;
        cout<<"\tlow "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[LowStopBin]<<endl;
      }
      peakTime.push_back(MidStartBin);
      peakTime.push_back(MidStopBin);
      peakFlag.push_back(1);
      mid1 = false;mid2 = false;
      if(ll)
        hh = true;
      else 
        hh = false;
    }

    if( ll && hh) {
      if(zeroCrossing>1){
        zeroCrossing = 0;
        ll = false;
        hh = false;
        low1 = false; low2 = false;
        mid1 = false;mid2 = false;
        continue;
      }
      if(peakFindingDebug) cout<<"llhh "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[HighStopBin]<<" zeroCrossing "<<zeroCrossing<<endl;
      peakFlag.erase(peakFlag.end()-2,peakFlag.end());
      peakFlag.push_back(0);
      peakTime.erase(peakTime.end() - 4,peakTime.end() );
      peakTime.push_back(LowStartBin);
      peakTime.push_back(MidStopBin);
      low1 = false; low2 = false;
      mid1 = false;mid2 = false;
      ll = false;hh = false;
      //LowStartBin = -1;LowStopBin = -1;HighStartBin = -1;HighStopBin = -1;
    }
   
  }
  //find local minimum around found peaks
  //also merge peaks
  for(int i = 0; i < peakTime.size(); i +=2){
    Int_t lowBin = peakTime[i];
    Int_t hiBin = peakTime[i+1];
    Double_t startTime = pmtEvent->time[lowBin];
    Double_t stopTime = pmtEvent->time[hiBin];
    /*
    if( peakFlag[i/2] == 1 && (hiBin-lowBin)*deltaT < 1.5*minPeakWidth){
      if(peakFindingDebug) cout<<"Removing UP-TYPE small pulse start "<<pmtEvent->time[lowBin]<<", stop "<<pmtEvent->time[hiBin]<<endl;
      continue;
    }
    */
    if( (hiBin-lowBin)*deltaT < minPeakWidth){
      if(peakFindingDebug) cout<<"Removing small pulse start "<<pmtEvent->time[lowBin]<<", stop "<<pmtEvent->time[hiBin]<<endl;
      continue;
    }
    //*/
    Double_t peakThresh =1*rms;
    Int_t shiftBin = sig.size()*1e-3;
    if(shiftBin == 0) shiftBin = 1;
    shiftBin = 1;
    ///*
    for(int k = lowBin; k > 1; k--){
      if( ((diff[k] < 0 && diff[k-1] > 0) || (diff[k] > 0 && diff[k-1] < 0))){
        zeroCrossing++;
      }
      //#1 crossing
      if( (peakFlag[i/2] != 1) && (diff[k] < -peakThresh && diff[k-1] > -peakThresh)){
      //if((diff[k] < peakThresh && diff[k-1] > peakThresh) ){
      //if( (fabs(diff[k]) < peakThresh && fabs(diff[k-1]) > peakThresh)){
        if(peakFindingDebug)cout<<"extending low bin from "<<pmtEvent->time[lowBin]<<", to "<<pmtEvent->time[k]<<", peakFlag "<<peakFlag[i/2]<<endl;
        if(peakFindingDebug)cout<<"\tdiff[k] "<<diff[k]<<", diff[k-1] "<<diff[k-1]<<" peakThresh "<<peakThresh<<endl;
        if(sig[k] > 0)
          lowBin = k + shiftBin;//10;
        else
          lowBin = k;
        break;
      }
      //#1 crossing
      else if( (peakFlag[i/2] == 1) && ( (diff[k] < -peakThresh && diff[k-1] > -peakThresh)) ){
      //if( (diff[k] < peakThresh && diff[k-1] > peakThresh)){
      //if( (fabs(diff[k]) > peakThresh && fabs(diff[k-1]) < peakThresh)){
        if(peakFindingDebug)cout<<"extending low bin from "<<pmtEvent->time[lowBin]<<", to "<<pmtEvent->time[k]<<", peakFlag "<<peakFlag[i/2]<<endl;
        if(peakFindingDebug)cout<<"\tdiff[k] "<<diff[k]<<", diff[k-1] "<<diff[k-1]<<" peakThresh "<<peakThresh<<endl;
        if(sig[k] > 0)
          lowBin = k - shiftBin;//10;
        else
          lowBin = k;
        break;
      }
    }
    //keep maximum for baseline fitting
    for(int k = hiBin; k < diff.size()-1; k++){
      //#4 crossing
      if((diff[k+1] < 0*peakThresh && diff[k] > 0*peakThresh) ){
        if(peakFindingDebug)cout<<"extending hi bin from "<<pmtEvent->time[hiBin]<<", to "<<pmtEvent->time[k]<<", peakFlag "<<peakFlag[i/2]<<endl;
        if(peakFindingDebug)cout<<"\tdiff[k+1] "<<diff[k+1]<<", diff[k] "<<diff[k]<<" peakThresh "<<peakThresh<<endl;
        if(sig[k] > 0)
          hiBin = k - shiftBin;//10;
        else
          hiBin = k;
        break;
      }
    }
    //*/
    //cull identical pulses
    if(peakTimeTrimmed.size() > 0 && lowBin == peakTimeTrimmed[peakTimeTrimmed.size()-2] && hiBin == peakTimeTrimmed[peakTimeTrimmed.size()-1]){
      if(peakFindingDebug) cout<<"Culling identical pulses..("<<pmtEvent->time[lowBin]<<","<<pmtEvent->time[hiBin]<<")..and.."<<pmtEvent->time[peakTime[i-2]]<<","<<pmtEvent->time[peakTime[i-1]]<<endl;
      continue;
    }
    //cull pulses inside other pulses
    if(peakTimeTrimmed.size() > 0 && lowBin < peakTimeTrimmed[peakTimeTrimmed.size()-1]){
      if(peakFindingDebug)
        cout<<"Pulse ("<<pmtEvent->time[lowBin]<<","<<pmtEvent->time[hiBin]<<") starts inside this pulse ("<<pmtEvent->time[peakTimeTrimmed[peakTimeTrimmed.size()-2]]<<","<<pmtEvent->time[peakTimeTrimmed[peakTimeTrimmed.size()-1]]<<")"<<endl;
      if(sig[peakTimeTrimmed[peakTimeTrimmed.size()-1]] > 0)
        lowBin = peakTimeTrimmed[peakTimeTrimmed.size()-1]+1;
      else
        lowBin = peakTimeTrimmed[peakTimeTrimmed.size()-1];
      //continue;
    }

    if(peakFindingDebug) cout<<"#### Pushing back lowBinTime "<<pmtEvent->time[lowBin]<<", highBin "<<pmtEvent->time[hiBin]<<" ####"<<endl;
    peakTimeTrimmed.push_back(lowBin);
    peakTimeTrimmed.push_back(hiBin );
  }

  return peakTimeTrimmed;
  //return peakTime;
}

//calc baseline and exclude region where peaks are found
std::pair<Double_t,Double_t> anaSim::BaselineSubtraction(std::vector<Double_t> sig, std::vector<Int_t> weight, Int_t pmtNum, Int_t ientry){
  Int_t startBin = 0, stopBin = 0;
  
  for(int i = 0; i < weight.size() - 1; i += 2){
    startBin = weight[i];
    stopBin = weight[i+1];
    for(int k = 0; k< sig.size(); k++){
      if(k >= startBin && k <= stopBin)
        sig.erase(sig.begin()+k);
    } 
  }
  //std::sort(sig.begin(),sig.end());

  return std::make_pair(CalculateMean(sig),CalculateSdev(sig,0,sig.size(),CalculateMean(sig)) );
}
//calc moving bassline
//P. Funk, G. Funk ect...
std::vector<Double_t> anaSim::BaselineSubtraction(std::vector<Double_t> sig, Int_t weight, Int_t pmtNum, Int_t ientry){
  std::vector<Double_t> baseline;
  std::vector<Double_t> gaus;
  //weight should be even
  if(weight%2 == 1) weight++;

  Double_t sigma = weight;
  Double_t norm = 0;
  for(int i = -weight; i <= weight; i++){
    Double_t val = exp(-i*i/(2*sigma*sigma));
//    gaus.push_back((1/sqrt(sigma*sigma*3.1415*2))*exp(-i*i/(2*sigma*sigma)));
    gaus.push_back(val);
    norm += val;
  }
  std::vector<Double_t> median;
  //for(int i = sig.size() - 1; i >= 0; i--){
  for(int i = 0; i < sig.size(); i++){
    std::vector<Double_t> localMedian;
    for(int j = i-weight; j <= i+weight; j++){
      if(j < 0) continue;
      if(j >= sig.size() ) continue;
      localMedian.push_back(sig[j]);  
    }
    //sorts vector from lowest to highest
    std::sort(localMedian.begin(),localMedian.end());
    //middle value is the median of the window per Local Median method
    //tuned value to 0.8 to better match long/large peak baselines
    //N. McFadden 
    median.push_back(localMedian[(int)localMedian.size()*.5]);
    //median.push_back(localMedian[localMedian.size()*.2]);
  }
  for(int i = 0; i < median.size(); i++){
  //for(int i = sig.size() - 1; i >= 0; i--){
    Double_t val = 0;
      for(int j = i-weight; j <= i+weight; j++){
        if(j < 0) continue;
        if(j >= median.size() ) continue;
        val += median[j]*gaus[j+weight-i]/norm;
        //cout<<"\t"<<median[j]<<" "<<gaus[j+weight-i]/norm<<" "<<i<<" "<<j<<" "<<j+weight-i<<endl;;
      }
      baseline.push_back(val);
  }
  //hBaseline[pmtNum] = new TH1D(TString("Baseline_Channel_")+to_string(pmtNum)+TString("_event_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
  for(int i = 0; i < baseline.size();i++){
    //Int_t timeBin = hBaseline[pmtNum]->FindBin(pmtEvent->time[i]);
    //hBaseline[pmtNum]->SetBinContent(timeBin,-baseline[i]);
  }
  return baseline;
}
//calculate baseline of whole waveform
std::pair<Double_t,Double_t> anaSim::BaselineSubtraction(std::vector<Double_t> sig, Int_t pmtNum, Int_t ientry){
  std::sort(sig.begin(),sig.end());
  return std::make_pair(sig[sig.size()/2],sig[sig.size()*.68]);
}

// seems like it is some type of smoothing filter...
std::vector<Double_t> anaSim::TrigFilter(std::vector<Double_t> sig, Int_t N, Int_t pmtNum, Int_t ientry){
  
  if(N%2 != 0) N++;
  std::vector<Double_t> vec(sig.size(),0); 
  for(int i = 0; i < sig.size(); i++){
    Double_t topSum = 0,botSum = 0;
    Int_t window = min(i,N);
    for(int j = i-window/2; j <=i+window/2; j++){
      if(j < 0) continue;
      if(j >= sig.size() ) continue;
      topSum += sig[j]*(1+std::cos((j-i)*2*pi/N));
      botSum += (1+std::cos((j-i)*2*pi/N));
    }
    if(botSum == 0) continue;
    vec[i] = topSum/botSum;
  }

  return vec;
}

Double_t anaSim::CalculateMean(std::vector<Double_t> vec){
  Double_t mean = 0;
  for(int i = 0; i< vec.size(); i++) mean += vec[i];
  return mean/vec.size();
}

Double_t anaSim::CalculateMean(std::vector<Int_t> vec){
  Double_t mean = 0;
  for(int i = 0; i< vec.size(); i++) mean += vec[i];
  return mean/vec.size();
}
Double_t anaSim::CalculateMean(std::vector<Double_t> vec,Int_t startBin, Int_t stopBin){
  Double_t mean = 0;
  Double_t N = (stopBin-startBin);
  for(int i = startBin; i< stopBin; i++) mean += vec[i];
  return mean/N;
}

Double_t anaSim::CalculateSdev(std::vector<Double_t> vec,Int_t startBin, Int_t stopBin, Double_t mean){
  Double_t sdev = 0;
  for(int i = startBin; i < stopBin; i++){
    sdev += (mean-vec[i])*(mean-vec[i]);
  }
  return std::sqrt(sdev/(stopBin-startBin));
}

std::vector<Double_t> anaSim::BubbleSort(std::vector<Double_t> A){
  int i, j, N = A.size();

  for (i = 0; i < N; i++){
    for (j = i + 1; j < N; j++){
      if (A[j] < A[i]){
        Double_t buff;
        buff = A[i];
        A[i] = A[j];
        A[j] = buff;
      }
    }
  }
  return A;
}

std::vector<Double_t> anaSim::TrapFilter(std::vector<Double_t> sig, Int_t ramp, Int_t flat,Int_t pmtNum, Int_t ientry){
  std::vector<Double_t> filter;
  if(flat%2 != 0) flat++;
  for(int i = 0; i < sig.size(); i++){
    Double_t val1 = 0,val2 = 0;
    Int_t minArr [] = {ramp,i,(Int_t)sig.size() - 1-i};
    Int_t min = *std::min_element(minArr,minArr+3);
    if(i < ramp + flat || i > sig.size() - ramp - flat - 1){
      for(int j = 1; j < min ;j ++){
        val1 += sig[i-j];
        val2 += sig[i+j];
      }
    }
    
    else{
      for(int j = 0; j < min ; j++){
        val1 += sig[i-j-flat];
        val2 += sig[i+j+flat];
      }
    }
    filter.push_back((val2-val1)/(min+1));
  }
    /*
    for(int j = i; j < ramp+i; j++){
      if(j >= sig.size()) continue;
      val1 += sig[j];
    }
    val1 /= ramp;
    for(int j = i +ramp +flat; j < i+ 2*ramp + flat; j++){
      if(j >= sig.size()) continue;
      val2 += sig[j];
    }
    val2 /= ramp;
    filter.push_back(val2-val1); 
  }
  */
  /*
  hFilter[pmtNum] = new TH1D(TString("Trap_Filter_Channel_")+to_string(pmtNum)+TString("_event_")+to_string(ientry)+TString("_flat_")+to_string(flat)+TString("_ramp_")+to_string(ramp),"",filter.size(),0,filter.size());
  for(int i = 0; i < filter.size(); i++){
    hFilter[pmtNum]->SetBinContent(i,filter[i]);
  }
  */

  return filter;
}
