#include "anaSimBackup.hh"
//

anaSimBackup::anaSimBackup(Int_t z){
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
  if(z >= 1000 && z <=1002)
    fileName = TString("baconRun_10kBins_1us_20mV_muon_perpendicular_")+to_string(z)+TString(".root");
  else if(z >= 1003 && z <=1007)
    fileName = TString("baconRun_10kBins_1us_20mV_muon_parallel_")+to_string(z)+TString(".root");
  else if(z>=2000 && z < 4000)
    fileName = TString("baconRun_10kBins_1us_20mV_div_-30mV_thresh_") + to_string(z) +TString(".root");
  else if(z>=4000 && z < 6000)
    fileName = TString("baconRun_10kBins_1us_20mV_div_-6mV_thresh_") + to_string(z) +TString(".root");
  else if(z>=10000)
    fileName = TString("baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_")+ to_string(z) +TString(".root");
  else if (z >= 100 && z <=272)
    fileName = TString("baconRun_10kBins_1us_10mV_")+ to_string(z) +TString(".root");
  else
    fileName = "simEvents_20190308_10000.root";//simEvents_20190226_10000.root";//"simEvents_20190227_100.root";//"simEvents_20190226_10000.root";
  TString outFileName = TString("/home/nmcfadde/RooT/PMT/processedData/SimAnaResults.")+fileName;
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());
 
  TNtuple *ntuplePulse = new TNtuple("ntuplePulse","ntuplePulse","ientry:pmtNum:nhits:charge:startTime:peakWidth:T0:vMax:vMaxTime:Sdev:baseline");
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

    signal[0] = pmtEvent->volt1;
    NEvents = pmtEvent->volt1.size();
    deltaT = (pmtEvent->time[NEvents -1] - pmtEvent->time[0])/(NEvents-1);
    if(ientry == 0) cout<<"time bin size is "<<deltaT<<endl;
    peakFindingDebug = false;
    /*
    if( ientry == 1){
      if(eventInfo != NULL)
        eventInfo->Write();
      break;
    }
    peakFindingDebug = true;
    Int_t skipBin = 0;
    NHistograms = skipBin;
    //if(skipBin < ientry ) break;
    if( skipBin != ientry) continue;
    hSum[0] = new TH1D(TString("Sum_")+to_string(0),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
    hSumFresh[0] = new TH1D(TString("SumVirgin_")+to_string(0),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
    */
    //cout<<"starting event"<<endl;

    //TH1F * hBaseFit[1];
    TH1D* hPeakFinding[1];
    Int_t nDer = 15, nInt =10;
    Double_t sum = 0;
    //cout<<"Filling signal"<<endl;
    for(int j = 0; j < signal.size();j++){
      hRawSignal[j] = new TH1D(TString("RawSignal_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]); 
      hSignal[j] = new TH1D(TString("Signal_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      hIntegral[j] = new TH1D(TString("Integral_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      if(ientry == 0) hSum[j] = new TH1D(TString("Sum_")+to_string(j),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      //hBaseFit[j] = new TH1F(TString("baseFit_")+to_string(j)+TString("_")+to_string(ientry),TString("baseFit")+to_string(j),NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      hPeakFinding[j] = new TH1D(TString("PeakFinding")+to_string(ientry),TString("PeakFinding")+to_string(ientry),NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      /*
      for(int i = 0; i < signal[j].size(); i++){
        hBaseFit[j]->SetBinContent(i+1,signal[j][i]);
      }
      */
      //cout<<signal[j].size()<<endl;
      //signal[j] = TrigFilter(signal[j],100,j,ientry);
      //cout<<signal[j].size()<<endl;
      //hBaseFit[j]->Fit(fPoly[j],"Q");
      std::vector<Double_t> integral;
      for(int i = 0; i < signal[j].size(); i++){
        //signal[j][i] = signal[j][i] - fPoly[j]->Eval(pmtEvent->time[i]);
        hRawSignal[j]->SetBinContent(i+1,signal[j][i]);
        //hSum[j]->SetBinContent(i+1,hSum[j]->GetBinContent(i+1)-signal[j][i]);
        sum += signal[j][i];
        integral.push_back(sum/nInt);
      }
      //signal[j] = TrapFilter(integral,40,5,j,ientry);
      if(signal[j].size() == 10000)
        signal[j] = Derivative(integral,nInt);
      else nInt = 1;
        signal[j] = Derivative(integral,nInt);
      //signal[j] = DownSampler(signal[j],5);
      if(signal[j].size() != 10000) nDer = 2;
      derivative[j] = Derivative(signal[j],nDer);
      hDerivative[j] = new TH1D(TString("derivative_")+to_string(j)+TString("_")+to_string(ientry),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
    }
    sum = 0;
    //cout<<"Setting derivative"<<endl;
    for(int j = 0; j < derivative.size(); j++){
      for(int i = 0; i < derivative[j].size(); i++){
        Int_t timeBin = hDerivative[j]->FindBin(pmtEvent->time[i]);
        //signal[j][i] /= (double) nInt;
        sum += signal[j][i];
        hSignal[j]->SetBinContent(i,signal[j][i]);
        hPeakFinding[j]->SetBinContent(i+1,0);
        hDerivative[j]->SetBinContent(i+1,derivative[j][i]);
        hIntegral[j]->SetBinContent(i+1,sum*deltaT);
      }
    }
    //hSignal[0]->Draw();
    //cout<<"Peak Finding"<<endl;
    //Peak Finding
    std::vector< std::vector<Int_t> > peakTime;peakTime.resize(signal.size());
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
          peakTime[i].push_back(pTime[j+1]);
        }
      }
    }
    //cout<<"finished peakFinding"<<endl;
    //2N + 1 moving window for WMA
    Int_t movingWindow = 100;
    if(signal[0].size() == 500)
      movingWindow = 10;
    std::vector<Double_t> BaSeLiNe = BaselineWMA(signal[0],peakTime[0],ientry,movingWindow);
    for(int i = 0; i < signal[0].size();i++){
      if(peakTime[0].empty()) continue;
      hSum[0]->SetBinContent(i+1,hSum[0]->GetBinContent(i+1)-signal[0][i]+BaSeLiNe[i]);
    }
    //cout<<"Filling ntuple"<<endl;
    //filling ntuple
    for(int j = 0; j < signal.size(); j++){
      Double_t vMaxEvent = 0,tMaxEvent = 0;
      TString sTitle = TString("");
      Double_t baseline = CalculateMean(signal[j],0,signal[j].size()*0.05);//500);
      Double_t Sdev =  CalculateSdev(signal[j],0,signal[j].size()*0.05,baseline);
      //Fill event information
      //ntupleEvent->Fill(z,ientry,j,Sdev,baseline,-sum,deltaT,peakTime.size()/2);
      
      if(peakTime[j].empty()) continue;
      Double_t nFound = 0,nMiss = 0,nNoise = 0,totalCharge = 0;
      Double_t T0=0;
      //cout<<"Finding peaks charge and voltage "<<ientry<<endl;
      for(int i = 0; i < peakTime[j].size()-1; i += 2){
        Double_t vMax = 0,charge = 0,vMaxTime = 0,startTime,stopTime,peakWidth;
        Int_t startBin = peakTime[j][i],stopBin = peakTime[j][i+1];
        //Int_t startBin = hSignal[j]->FindBin(400e-9);
        //Int_t stopBin  = hSignal[j]->FindBin(500e-9);
        startTime = pmtEvent->time[startBin];
        if(i == 0) T0 = startTime;
        stopTime = pmtEvent->time[stopBin];
        peakWidth = stopTime - startTime;
        for(int k = startBin; k < stopBin; k++){
          if(signal[j][k] - BaSeLiNe[k] > 0) continue;
          charge += signal[j][k]-BaSeLiNe[k];
          if(std::fabs(signal[j][k]) > std::fabs(vMax) ){
            vMax = signal[j][k]-BaSeLiNe[k];
            vMaxTime = pmtEvent->time[k];
          }
          hPeakFinding[j]->SetBinContent(k,signal[j][k]-BaSeLiNe[k]);
        }
        if(peakWidth < 0){
          cout<<"Peak Width is negative, irun = "<<z<<", entry = "<<ientry<<", peakWidth = "<<peakWidth<<", start/stop "<<startTime<<"/"<<stopTime<<", startBin/stopBin "<<startBin<<"/"<<stopBin<<endl;
          cout<<"\tstart "<<pmtEvent->time[startBin]<<", stop "<< pmtEvent->time[stopBin]<<" vmax "<< -vMax<<" charge "<<charge*deltaT*-1e9<<endl;
        }
        if(peakFindingDebug) cout<<"start "<<pmtEvent->time[startBin]<<", stop "<<  pmtEvent->time[stopBin]<<" vmax "<< -vMax<<" charge "<<charge*deltaT*-1e9<<endl;
        ntuplePulse->Fill(ientry,j,peakTime[j].size()/2,-charge*deltaT,startTime,peakWidth,T0,-vMax,vMaxTime,0,0);
        sTitle += TString("Charge_")+to_string(-1e9*deltaT*charge)+TString("_start/stop_")+to_string(1.e6*startTime)+TString("/")+to_string(1.e6*stopTime)+TString("_vMax_")+to_string(-vMax)+TString("_");
        if(std::fabs(vMax) > std::fabs(vMaxEvent) ){
          vMaxEvent = vMax;
          tMaxEvent = vMaxTime;
        }
        //cout<<"Setting sim data"<<endl;
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
      sTitle += TString("M_") + to_string(nFound)+TString("/G")+to_string(simEvent->startTime.size());
      ntupleEvent->Fill(z,ientry,Sdev,baseline,-sum,deltaT,tMaxEvent,-vMaxEvent);//-vMaxEvent,tMaxEvent);
      for(int i = 0; i < simTimes.size(); i++) sTitle +=TString("_") + to_string(1e6*simTimes[i]);
      hSignal[j]->SetTitle(sTitle);
    }
//    cout<<"finished filling ntuple"<<endl;
    hPeakFinding[0]->SetLineColor(2);
    //hPeakFinding[0]->Draw("same");
    //###################//
    //End of Event Clean-Up//
    //###################//

    if(ientry > NHistograms){
      for(int i = 0; i < signal.size(); i++){
        delete hSignal[i];
        delete hDerivative[i];
        delete hRawSignal[i];
        //delete hBaseFit[i];
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
  cout<<"root -l "<<outFileName<<endl;
  //ntuplePulse->Write();
  outFile->Write();
  //}
  //cout<<"root -l "<<outFileName<<endl;
  cout<<"you did it! your code ran with crashing"<<endl;
}
/*
std::vector<Double_t> anaSimBackup::Derivative(std::vector<Double_t> sig,Int_t N){
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

std::vector<Double_t> anaSimBackup::DownSampler(std::vector<Double_t> sig,Int_t N){

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

std::vector<Double_t> anaSimBackup::Derivative(std::vector<Double_t> sig,Int_t N){
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

Double_t anaSimBackup::RMSCalculator(std::vector<Double_t> vec,Int_t ientry){
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
std::vector<Double_t> anaSimBackup::BaselineWMA(std::vector<Double_t> sig,std::vector<Int_t> peaks,Int_t ientry,Int_t N){
  std::vector<Double_t> Baseline;
  if(peaks.size() == 0) return Baseline;
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
  if(N > 100)
    N/=1.5;
//  N += 10;
  //cout<<N<<endl;
  TH1D * f1 = new TH1D(TString("WMA")+to_string(ientry),TString("WMA")+to_string(ientry)+TString("_N_")+to_string(N),sig.size(),0,pmtEvent->time[NEvents-1]);
  TH1D * f2 = new TH1D(TString("Weight")+to_string(ientry),TString("Weight")+to_string(ientry),sig.size(),0,pmtEvent->time[NEvents-1]);
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
  Int_t resetBin = 100;
  if(sig.size() == 500) resetBin = 20;
  for(int i = 0;i< sig.size();i++){
    Int_t hiWindow = std::min(i+N,(Int_t)sig.size() -1);
    Int_t loWindow = std::max(0,i-N);
    if(i%resetBin ==  0){
      K_sw = 0;C_sw = 0;S_sw = 0;K_w = 0;C_w = 0;S_w = 0;
      for(int j = loWindow;j<hiWindow;j++){
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
  //f1->Draw("same");
  //f1->Write();
  if(ientry > NHistograms){
    delete f2;
    delete f1;
  }

  return Baseline;
}


std::vector<Double_t> anaSimBackup::Integral(std::vector<Double_t> sig,Int_t pmtNum, Int_t ientry){
  std::vector<Double_t> integ;
  Double_t sum = 0;
  for(int i = 0; i < sig.size(); i++){
    sum += sig[i];
    integ.push_back(sum);
  }
  return integ;
}
std::vector<Int_t> anaSimBackup::PeakFinding(std::vector<Double_t> sig, std::vector<Double_t> diff,Int_t pmtNum,Int_t ientry){
  std::vector<Int_t> peakTime;
  std::vector<Int_t> peakTimeTrimmed;
  /*
  std::vector<Double_t> sort = diff;
  std::sort(sort.begin(),sort.end());
  Double_t RMS = 5*std::fabs(sort[sort.size()*.68]);
  hDerivative[pmtNum]->SetTitle(TString("Threshold_")+to_string(RMS));
  */
  //Int_t rmsBin = hDerivative[pmtNum]->FindBin(150e-9);
  Int_t rmsBin = sig.size()*0.05;//500;//sig.size();
  //Double_t rmsFit = RMSCalculator(diff,ientry);
  Double_t rmsSdev = CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  Double_t rms =rmsSdev;// min(rmsFit,rmsSdev);
  Double_t loRMS = 3.5*rms;//CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  Double_t hiRMS = 3.5*rms;//CalculateSdev(diff,0,rmsBin,CalculateMean(diff));
  //cout<<"RMS = "<<rms<<endl;
  hDerivative[pmtNum]->SetTitle(TString("loThreshold_")+to_string(loRMS)+TString("hiThreshold_")+to_string(hiRMS)+TString("_mean_")+to_string(CalculateMean(diff)) );
  Bool_t low1 = false,low2 = false,high1 = false, high2 = false, ll = false, hh = false;

  Int_t LowStartBin = -1,LowStopBin = -1,HighStartBin = -1,HighStopBin = -1;
  Int_t zeroCrossing = 0;
  for(int i = 1; i < diff.size(); i++){
    if(diff[i] <= -loRMS && diff[i-1] > -loRMS){
      LowStartBin = i;
      low1 =  true;
      //cout<<pmtEvent->time[LowStartBin]<<"  "<<diff[i]<<" "<<diff[i-1]<<endl;
    }
    if(diff[i] >= -loRMS && diff[i-1] < -loRMS ){
      LowStopBin = i;
      low2 = true;
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
      //cout<<"low "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[LowStopBin]<<" "<<diff[i]<<" "<<-RMS<<endl;
      //cout<<"\thi "<<pmtEvent->time[HighStartBin]<<" "<<pmtEvent->time[HighStopBin]<<endl;
      /*
      for(int k = LowStartBin; k > 1; k--){
        if(diff[k] < 0 && diff[k-1] > 0){
          //cout<<"LowStartBin "<<LowStartBin<<" k "<<k<<endl;
          LowStartBin = k;
          break;
        }
      }
      for(int k = LowStopBin; k < diff.size(); k++){
        if(diff[k] < 0 && diff[k-1] > 0){
          //cout<<"LowStopBin "<<LowStopBin<<" k "<<k<<endl;
          LowStopBin = k;
          break;
        }
      }
      */
      peakTime.push_back(LowStartBin);
      peakTime.push_back(LowStopBin);
      low1 = false; low2 = false;
      if(hh)
        ll = false;
      else{
        ll = true;
        //reset zero crossing if a new lolo peak is found
        zeroCrossing = 0;
      }
      //LowStartBin = -1;LowStopBin = -1;
    }
    
    if(high1 && high2){
      //cout<<"hi "<<pmtEvent->time[HighStartBin]<<" "<<pmtEvent->time[HighStopBin]<<" "<<diff[i]<<" "<<RMS<<endl;
      //cout<<"\tlow "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[LowStopBin]<<endl;
      /*
      for(int k = HighStartBin; k > 1; k--){
        if(diff[k] < 0 && diff[k-1] > 0){
          //cout<<"HighStartBin "<<HighStartBin<<" k "<<k<<endl;
          HighStartBin = k;
          break;
        }
      }
      for(int k = HighStopBin; k < diff.size(); k++){
        if(diff[k] < 0 && diff[k-1] > 0){
          //cout<<"HighStopBin "<<HighStopBin<<" k "<<k<<endl;
          HighStopBin = k;
          break;
        }
      }
      */
      peakTime.push_back(HighStartBin);
      peakTime.push_back(HighStopBin);
      high1 = false;high2 = false;
      if(ll)
        hh = true;
      else 
        hh = false;
      //HighStartBin = -1;HighStopBin = -1;
    }
    
    if( ll && hh) {
      //cout<<"llhh "<<pmtEvent->time[LowStartBin]<<" "<<pmtEvent->time[HighStopBin]<<endl;
      if(zeroCrossing>1) continue;
      zeroCrossing = 0;
      peakTime.erase(peakTime.end() - 4,peakTime.end() );
      peakTime.push_back(LowStartBin);
      peakTime.push_back(HighStopBin);
      low1 = false; low2 = false;
      high1 = false;high2 = false;
      ll = false;hh = false;
      //LowStartBin = -1;LowStopBin = -1;HighStartBin = -1;HighStopBin = -1;
    }
   
  }
  //find local minimum around found peaks
  //also merge peaks
  for(int i = 0; i < peakTime.size(); i +=2){
    Int_t lowBin = peakTime[i];
    Int_t hiBin = peakTime[i+1];
    Int_t shiftBin = sig.size()*1e-3;
    if(shiftBin == 0) shiftBin = 1;
    ///*
    for(int k = lowBin; k > 1; k--){
      if(diff[k] < 0 && diff[k-1] > 0){
        //cout<<"lowBin "<<lowBin<<" k "<<k<<endl;
        if(sig[k] > 0)
          lowBin = k + shiftBin;//10;
        else
          lowBin = k;
        break;
      }
    }
    //keep maximum for baseline fitting
    for(int k = hiBin; k < diff.size(); k++){
      //cout<<"hiBin "<<hiBin<<" k "<<k<<" diff size "<<diff.size()<<endl;
      if(diff[k] < 0 && diff[k-1] > 0){
        //cout<<"hiBin "<<hiBin<<" k "<<k<<endl;
        if(sig[k] > 0)
          hiBin = k - shiftBin;//10;
        else
          hiBin = k;
        break;
      }
    }
    //*/
    ///*
    //cull identical pulses
    if(i > 1 && lowBin == peakTime[i-2] && hiBin == peakTime[i-1]){
      continue;
    }
    //cull pulses inside other pulses
    if(lowBin < peakTime[i-1]){
      continue;
    }
    /*
    Double_t peakWidth = (hiBin-lowBin)*deltaT;
    if(peakWidth < minPeakWidth && peakWidth > maxPeakWidth)
      continue;
    */
   
    //*/
    //cout<<"peakTimeTrimmed.size() "<<peakTimeTrimmed.size()<<", i "<<i<<", lowBin "<<lowBin<<"/"<<peakTime[i]<<", hiBin "<<hiBin<<"/"<<peakTime[i+1]<<endl;
    peakTimeTrimmed.push_back(lowBin);
    peakTimeTrimmed.push_back(hiBin );
  }

  return peakTimeTrimmed;
}

//calc baseline and exclude region where peaks are found
std::pair<Double_t,Double_t> anaSimBackup::BaselineSubtraction(std::vector<Double_t> sig, std::vector<Int_t> weight, Int_t pmtNum, Int_t ientry){
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
std::vector<Double_t> anaSimBackup::BaselineSubtraction(std::vector<Double_t> sig, Int_t weight, Int_t pmtNum, Int_t ientry){
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
std::pair<Double_t,Double_t> anaSimBackup::BaselineSubtraction(std::vector<Double_t> sig, Int_t pmtNum, Int_t ientry){
  std::sort(sig.begin(),sig.end());
  return std::make_pair(sig[sig.size()/2],sig[sig.size()*.68]);
}

// seems like it is some type of smoothing filter...
std::vector<Double_t> anaSimBackup::TrigFilter(std::vector<Double_t> sig, Int_t N, Int_t pmtNum, Int_t ientry){
  
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

Double_t anaSimBackup::CalculateMean(std::vector<Double_t> vec){
  Double_t mean = 0;
  for(int i = 0; i< vec.size(); i++) mean += vec[i];
  return mean/vec.size();
}

Double_t anaSimBackup::CalculateMean(std::vector<Int_t> vec){
  Double_t mean = 0;
  for(int i = 0; i< vec.size(); i++) mean += vec[i];
  return mean/vec.size();
}
Double_t anaSimBackup::CalculateMean(std::vector<Double_t> vec,Int_t startBin, Int_t stopBin){
  Double_t mean = 0;
  Double_t N = (stopBin-startBin);
  for(int i = startBin; i< stopBin; i++) mean += vec[i];
  return mean/N;
}

Double_t anaSimBackup::CalculateSdev(std::vector<Double_t> vec,Int_t startBin, Int_t stopBin, Double_t mean){
  Double_t sdev = 0;
  for(int i = startBin; i < stopBin; i++){
    sdev += (mean-vec[i])*(mean-vec[i]);
  }
  return std::sqrt(sdev/(stopBin-startBin));
}

std::vector<Double_t> anaSimBackup::BubbleSort(std::vector<Double_t> A){
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

std::vector<Double_t> anaSimBackup::TrapFilter(std::vector<Double_t> sig, Int_t ramp, Int_t flat,Int_t pmtNum, Int_t ientry){
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
