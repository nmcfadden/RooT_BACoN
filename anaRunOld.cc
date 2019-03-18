////////////////////////////////////////////////////////
#include "anaRun.hh"


anaRun::anaRun(Int_t irun)
{

  Double_t irunStart = 31, irunStop = 31,irunName = irunStart + .01*irunStop;

  TString outFileName ; outFileName.Form("baconRunAna%f.root",irunName);
  TFile *outfile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());
  
  ntupleRun = new TNtuple("ntupleRun","ntupleRun","irun:ientry:pmtNum:sigma:nhits:charge:startTime:peakWidth:qMax:qMaxTime:T0");
  ntupleEvent = new TNtuple("ntupleEvent","ntupleEvent","irun:ientry:pmtNum:time:volt");
  ntupleCal = new TNtuple("ntupleCal","ntupleCal","irun:ientry:pmtNum:tStart:tStop:tWindow:charge:zeroNoise:baseline:Sdev");


  for(irun = irunStart; irun<irunStop +1; irun++){
    // open ouput file and make some histograms
    TString fileName; fileName.Form("rootData/baconRun_%i.root",irun);
    printf(" looking for file %s\n",fileName.Data());
    TFile *fin = new TFile(fileName,"readonly");
    if(fin->IsZombie()) {
      printf(" couldnt open file %s\n",fileName.Data());
      irun++;
      if(irun == irunStop)
        return;
    }
    else 
      printf("  found file %s \n",fileName.Data() );

    // get pmtTree from file 
    fin->GetObject("pmtTree",pmtTree);
    Long64_t nentries = pmtTree->GetEntries();
    cout << " number of entries is " << nentries << endl;

    // set up memory for reading
    pmtEvent = new TPmtEvent();
    pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);
    
    //switch to output file
    outfile->cd();

    //Event Vectors

    Int_t movingN1 = 0;
    Int_t movingN2 = 0;
    Double_t runningNoise1 = 0;
    Double_t runningNoise2 = 0;
    std::vector<Double_t> runningNoise(2);
    std::vector<Int_t> movingN(2);


    // loop over entries 
    for (UInt_t ientry=0; ientry<nentries; ientry++) {
      pmtTree->GetEntry(ientry);
      //if(ientry != 2) continue;
      if(pmtEvent->time.size() == 0) continue;
      NEvents = pmtEvent->time.size();


      Double_t nameDouble = double(ientry)+.01*double(irun);
      TString name; name.Form("%f", nameDouble );
      
      if( ientry == 0){ 
        hPMTSum1 = new TH1D("PMT_Sum_0" + name,"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
        hPMTSum2 = new TH1D("PMT_Sum_1" + name,"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      }

      
      //Be careful using Histograms, they get deleted after 1000 events!
      TH1D * hPMTSignal1 = new TH1D("PMTSignal1_" + name + TString("_baseline_subtracted"),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      TH1D * hPMTSignal2 = new TH1D("PMTSignal2_" + name + TString("_baseline_subtracted"),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
//      TH1D * hPMTSignalFiltered1 = new TH1D("PMTSignal1_" + name + TString("_Filtered"),"",NEvents-1,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
//      TH1D * hPMTSignalFiltered2 = new TH1D("PMTSignal2_" + name + TString("_Filtered"),"",NEvents-1,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      
      TH1D * hBaseline1 = new TH1D("Baseline1_"+ name,"",NEvents-1,-50e-3,50e-3);
      TH1D * hBaseline2 = new TH1D("Baseline2_"+ name,"",NEvents-1,-50e-3,50e-3);
      
      //TH1D * hNormline1 = new TH1D("Normline1_"+ name,"",NEvents-1,-1e-3,5e-3);
      //TH1D * hNormline2 = new TH1D("Normline2_"+ name,"",NEvents-1,-1e-3,5e-3);

      // initialize fft 
      fFFT = TVirtualFFT::FFT(1, &NEvents, "R2C M K");
      fInverseFFT = TVirtualFFT::FFT(1, &NEvents, "C2R M K");
      std::vector<Double_t>  sdigi1,ddigi1,fftdigi1;
      std::vector<Double_t>  sdigi2,ddigi2,fftdigi2;
      Double_t baseline1 = 0,sum1 = 0,T1 = 0,T2 = 0;
      Double_t baseline2 = 0,sum2 = 0;
      deltaT = 0;

      Int_t ZeroTime = hPMTSignal1->FindBin(hPMTSignal1->GetBinCenter(0));
      Int_t StopTime = hPMTSignal1->FindBin(-1.2e-6);//-200e-9);
      Int_t totTime = std::fabs(StopTime-ZeroTime);
        
      //For cal ntuple
      Int_t tStart = hPMTSignal1->FindBin(0e-9);
      Int_t tStop = hPMTSignal1->FindBin(500e-9);
      
      if(ientry == 0)  cout<<"ZeroTime = "<<hPMTSignal1->GetBinCenter(ZeroTime)<<", StopTime = "<<hPMTSignal1->GetBinCenter(StopTime)<<endl;

      /*
      for (UInt_t is=0; is<NEvents; ++is) {
        if(is > unsigned(ZeroTime) && is < unsigned(StopTime) ) {
          baseline1 += pmtEvent->volt1[is];
          baseline2 += pmtEvent->volt2[is];
        }
        if(is > 0)
          deltaT += std::fabs(pmtEvent->time[is] - pmtEvent->time[is - 1]);
      }*/

      
      std::vector<std::vector<Double_t> > ddigi;ddigi.resize(movingN.size() );
      
      for(Int_t ievent = 0; ievent < NEvents; ievent++){
        Double_t volt1 = pmtEvent->volt1[ievent];
        Double_t volt2 = pmtEvent->volt2[ievent];
        //double digi1 = -1.0*(double(volt1)-baseline1);
        //double digi2 = -1.0*(double(volt2)-baseline2);
        double digi1 = -1.0*(double(volt1));
        double digi2 = -1.0*(double(volt2));
        //if(digi1 < 0) digi1 = 0;
        //if(digi2 < 0) digi2 = 0;
        ddigi[0].push_back(digi1);
        ddigi[1].push_back(digi2);
        
        hBaseline1->Fill(-volt1);
        hBaseline2->Fill(-volt2);

        if(ievent > 0)
          deltaT += std::fabs(pmtEvent->time[ievent] - pmtEvent->time[ievent - 1]);
      }
      
      deltaT /= (NEvents - 1);
      minLength = int(backwardLengthTimeThreshold/deltaT);
      forwardHalfLength = int(forwardLengthTimeThreshold/deltaT);

      if(ientry == 0 ){ 
        cout<<"deltaT = "<<deltaT<<", Minimum Peak width = "
          <<deltaT*minLength<<", minHalfLength "<<minLength<<
          ", Max Peak width = "<<deltaT*forwardHalfLength<<", forwardHalfLength "<<forwardHalfLength<<endl;
      }

      std::vector<Double_t> baseline; baseline.resize(ddigi.size() );
      std::vector<Double_t> Sdev;Sdev.resize( ddigi.size() );
      
      //Q is quite mode! very important
      hBaseline1->Fit("gaus","Q");
      if(!hBaseline1->GetFunction("gaus") ) continue;
      baseline[0] = hBaseline1->GetFunction("gaus")->GetParameter(1);
      Sdev[0] = hBaseline1->GetFunction("gaus")->GetParameter(2);
      delete hBaseline1;

      hBaseline2->Fit("gaus","Q");
      if(!hBaseline2->GetFunction("gaus") ) continue;
      baseline[1] = hBaseline2->GetFunction("gaus")->GetParameter(1);
      Sdev[1] = hBaseline2->GetFunction("gaus")->GetParameter(2);
      delete hBaseline2;

      for(int j = 0; j < ddigi.size();j++){
        for(int i = 0; i < ddigi[j].size(); i++){
          ddigi[j][i] = ddigi[j][i]-baseline[j];// + Sdev[j];
          //ntupleEvent->Fill(irun,ientry,j,pmtEvent->time[i],ddigi[j][i]);

          if(j == 0) {
            hPMTSignal1->SetBinContent(i+1,ddigi[j][i]);
          }
          else{
            hPMTSignal2->SetBinContent(i+1,ddigi[j][i]);
          }

        }
      }
     
      Int_t nbins = NEvents;
      for(int j = 0; j < 2; j++){
        Double_t * source = new Double_t[ddigi[0].size()];
        Double_t * dest = new Double_t[ddigi[0].size()];
        Double_t * response = new Double_t[ddigi[0].size()];

        for(int i = 0; i < ddigi[0].size(); i++) source[i] = ddigi[j][i];
        Int_t nhits = spec->SearchHighRes(source, dest, nbins, 8, 2, kFALSE, 3, kFALSE, 3);
        Double_t * peaks = spec->GetPositionX();
        for(int i = 0; i < nbins; i++) ddigi[j][i] = dest[i];
      }
      /*
      for(int i = 0; i < ddigi.size(); i++){
        std::vector<std::complex<double> > fftPair = FFT(i,ientry,ddigi[i]);
        std::vector<Double_t> vecTemp = inverseFFT(i,ientry,fftPair,sum);
      }
      */    
     /* 
      for(int j = 0 ; j <ddigi.size(); j++){
        for(int i = 0; i < 4; i++){ 
          ddigi[j] = MovingAverageFilter(ddigi[j],5);
        }
      }
      for(int j = 0; j < ddigi.size(); j++){
        for(int i = 0; i < ddigi[j].size() ; i++){
          if(j == 0)
            hPMTSignalFiltered1->SetBinContent(i,ddigi[j][i]);
          else
            hPMTSignalFiltered2->SetBinContent(i,ddigi[j][i]);
        }
      }
      */
      /*
      if(ientry > 1000 ){
        delete hPMTSignal1;
        delete hPMTSignal2;
       // delete hPMTSignalFiltered1;
       // delete hPMTSignalFiltered2;
      }
      */
      std::vector<Double_t> tempPmtNum,tempSigma,tempNhits,tempQsum,tempStartTime,tempPeakWidth,tempqHit,tempPeakMaxTime;
      T0.clear();
      Double_t T0_temp = 1;
      for(int pmtNum = 0 ; pmtNum < ddigi.size();pmtNum++){
        for(int sigma = 2 ; sigma < 3; sigma++){
          Double_t minDev = 0*(sigma - 2)*Sdev[pmtNum], maxDev = 1*sigma*Sdev[pmtNum];
          std::vector<Int_t> peakTime = findPeaks(ddigi[pmtNum],maxDev,minDev,tStart,tStop );
          Int_t nhits = findHits( peakTime, ddigi[pmtNum]);
          

          if(ientry < 50 && pmtNum == 0 && 0 ){
            cout<<"\t\t"<<sigma<<"-SIGMA THRESHOLD == "<<maxDev<<endl;
            for(int v = 0; v < peakTime.size() ; v++){
              if(v > 0)
                if( (peakTime[v] - peakTime[v-1]) != 1) cout<<"\t new peak found"<<endl;
              cout<<ientry<<" = ientry, "<<sigma*Sdev[pmtNum]<<" = "<<sigma<<"-Sigma, "<<pmtEvent->time[peakTime[v]]<<" = time,"<<ddigi[pmtNum][peakTime[v] ]<<" = volts"<<endl;
            }
          }
          for(int i = 0; i < nhits; i++)
            ntupleRun->Fill(irun,ientry,pmtNum,sigma,nhits,qSum[i]*deltaT, startTime[i] ,peakWidth[i],qhitMax[i],peakMaxTime[i],T0_temp);
          if(nhits != 1 && pmtNum == 0) delete hPMTSignal1;
          if(nhits != 1 && pmtNum == 1) delete hPMTSignal2;
         /*
          for(int i = 0; i < nhits; i++){
            tempPmtNum.push_back(pmtNum);
            tempSigma.push_back(sigma);
            tempNhits.push_back(nhits);
            tempQsum.push_back(qSum[i]*deltaT);
            tempStartTime.push_back(startTime[i]);
            tempPeakWidth.push_back(peakWidth[i]);
            tempqHit.push_back(qhitMax[i]);
            tempPeakMaxTime.push_back(peakMaxTime[i]);
            T0.push_back(startTime[i]);
            //printf("ientry %i, sigma %e, qMax %e, qMaxTime %e \n",ientry,sigma*Sdev[pmtNum],qhitMax[i],peakMaxTime[i]);
          }
        }
        T0 = BubbleSort(T0);
        */
        //for(int i = 0; i < tempPmtNum.size(); i ++){
          //ntupleRun->Fill(irun,ientry,tempPmtNum[i],tempSigma[i],tempNhits[i],tempQsum[i], tempStartTime[i]-T0[0],tempPeakWidth[i],tempqHit[i],tempPeakMaxTime[i],T0[0]);
          //cout<<"Event "<<ientry<<", StartTime "<<tempStartTime[i]<< " T0 "<<T0[0]<<endl;
        }
        tempPmtNum.clear();tempSigma.clear();tempNhits.clear();tempQsum.clear();tempStartTime.clear();tempPeakWidth.clear();tempqHit.clear();tempPeakMaxTime.clear();T0.clear();
      }
      
      //cout<<"Start bin = "<<tStart<<", Start time = "<< pmtEvent->time[tStart]<<", Stop bin = "<<tStop<<", Stop time = "<< pmtEvent->time[tStop]<<endl;
      Double_t charge = 0,zeroNoise = 0;
      for(int j = 0; j < ddigi.size(); j++){
        for(int i = tStart ; i < tStop ; i++){
            charge += ddigi[j][i]*deltaT;
        }
        for(int i = ddigi[j].size() - std::fabs(tStart-tStop)-1 ; i < ddigi[j].size() ; i++){
          zeroNoise += ddigi[j][i]*deltaT;
        }
          ntupleCal->Fill(irun,ientry,j,tStart,tStop,std::fabs(tStart-tStop),charge,zeroNoise,baseline[j],Sdev[j]);
          charge = 0;
          zeroNoise = 0;
      }
      if(ientry%5000 == 0) cout<<ientry<<" Events Processed"<<endl;

      if(ientry == 1000) break;
      //if(ientry>25) break;
      //if(ientry > 50) break;
      //if(ientry>12364) break;//run_16 good data set
      if(ientry>40000) break;
    }
    //hPMTSum1->Scale(1./movingN[0]);
    //hPMTSum2->Scale(1./movingN[1]);
    cout<<movingN[0]<<" events processed out of "<<nentries<<endl;
    cout<<"\t WARNING SYSTEM IS IN CRITICAL CONDITION! \n\t NEEDS MORE CHEESE FOR MOUSE!"<<endl;
  }//irun loop

  outfile->Write();

return;
}

std::vector<std::complex<double> > anaRun::FFT(Int_t ipmt,Int_t ievent,std::vector<Double_t> signal)
{
  std::vector<std::complex<double> > VectorComplex;
  for(int is =0; is<NEvents; ++is) {
    if(ipmt==0) fFFT->SetPoint(is, signal[is]);
    else fFFT->SetPoint(is, signal[is]);
  }
  fFFT->Transform();
  if(ipmt == 0){
    hfft = new TH1D(Form("FFTChannelOne_%i",ievent),Form("FFT Channel One %i",ievent),NEvents/2,0,NEvents/2);
    //hfft = new TH1D(Form("FFTChannelOne_%i",ievent),Form("FFT Channel One %i",ievent),NEvents,0,NEvents);
  }
  else if(ipmt == 1){
    hfft = new TH1D(Form("FFTChannelTwo_%i",ievent),Form("FFT Channel Two %i",ievent),NEvents/2,0,NEvents/2);
    //hfft = new TH1D(Form("FFTChannelTwo_%i",ievent),Form("FFT Channel Two %i",ievent),NEvents,0,NEvents);
  }

  std::vector<Double_t> realVec,imVec;
  for (int i = 0; i<NEvents; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im); //.real or .imag accessors
    VectorComplex.push_back(c);
    // skip first bin which is pedestal
    //if( i> 0) hfft->SetBinContent(i+1,hfft->GetBinContent(i+1)+std::abs(c));
    hfft->SetBinContent(i,hfft->GetBinContent(i)+std::abs(c));
    realVec.push_back( VectorComplex[i].real());
    imVec.push_back(VectorComplex[i].imag() );
  }
  /*
  realVec = MovingAverageFilter(realVec,5);
  imVec   = MovingAverageFilter(imVec,5);
  for(int i = 0; i < NEvents; i++){
    std::complex<double> c(realVec[i],imVec[i] );
    VectorComplex[i] = c;
  }
  */

return VectorComplex;
}

std::vector< Double_t > anaRun::inverseFFT(Int_t ipmt,Int_t ievent,std::vector<std::complex<double> > VectorComplex,std::vector<Double_t> sum)
{
  std::vector<Double_t > Signal;
  for(int is =0; is<NEvents; ++is) {
    if(ipmt==0) fInverseFFT->SetPoint(is, VectorComplex[is].real(),VectorComplex[is].imag() );
    else fInverseFFT->SetPoint( is, VectorComplex[is].real(),VectorComplex[is].imag() );
  }
  fInverseFFT->Transform();
  if(ipmt == 0){
    hIfft = new TH1D(Form("inverseFFTChannelOne_%i",ievent),Form("InverseFFT_One_%i",ievent),NEvents,pmtEvent->time[0],pmtEvent->time[NEvents - 1]);
    //hfft = new TH1D(Form("FFTChannelOne_%i",ievent),Form("FFT Channel One %i",ievent),NEvents,0,NEvents);
  }
  else if(ipmt == 1){
    hIfft = new TH1D(Form("inverseFFTChannelTwo_%i",ievent),Form("InverseFFT_Two_%i",ievent),NEvents,pmtEvent->time[0],pmtEvent->time[NEvents - 1]);
    //hfft = new TH1D(Form("FFTChannelTwo_%i",ievent),Form("FFT Channel Two %i",ievent),NEvents,0,NEvents);
  }
  Double_t norm = 0;
  for (int i = 0; i<NEvents; ++i) {
    double rl = fInverseFFT->GetPointReal(i);
    norm += rl;
    Signal.push_back(rl);
  }
  for(int i = 0; i < NEvents; i++){
    Signal[i] = Signal[i]*(sum[ipmt]/norm);
    hIfft->SetBinContent(i,hIfft->GetBinContent(i) + Signal[i]);
  }

return Signal;
}

std::vector<Int_t> anaRun::findPeaks(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold,Int_t startBin,Int_t stopBin) 
{
  // Produces a list of peaks above the threshold
  std::vector<Int_t> peakTime;
  Int_t klow=0;
  Int_t khigh=0;
  Int_t kover=0;
  Int_t vsize = Int_t(v.size());
  Double_t fMin = 9999,bMin = 9999;
  
  //printf(" findPeaks \n");
  //for(Int_t  ibin=0; ibin<= vsize; ++ibin ) {
  for(Int_t  ibin=startBin; ibin< stopBin; ++ibin ) {
    if( v[ibin]>threshold ){//&& (v[ibin] > v[ibin - 1] ) && (v[ibin] < v[ibin + 1] ) ){// starting possible new hit and only find hits on a positive slope
      // consider this a "seed" and find full hit
      klow=ibin;
      for(Int_t k=ibin-1; k>=max(0,ibin-minLength); --k) {
        if(bMin > v[k]){
          bMin = v[k];
          klow = k;
        }
        if(v[k]<sthreshold){
          klow=k; 
          break;
        }
        if(k == max(0,ibin-minLength) ){
          klow=k;
        }
        //klow=k;
      }
      khigh=ibin;
      for(Int_t k=ibin+1; k<=min(ibin+forwardHalfLength,vsize); ++k) {
        if(fMin > v[k]){
          fMin = v[k];
          khigh = k;
        }
        if(v[k]<sthreshold){
          khigh = k; 
          break;
        }
        if(k == min(ibin+forwardHalfLength,vsize)){
          khigh=k;
        }
      }
      kover = khigh-klow+1;
      // found good pulse
      if(kover>minLength) { 
        for(Int_t k=klow ; k<= khigh; ++k){
          peakTime.push_back(k);
          fMin = 9999;
          bMin = 9999;
          //cout<<"hi time "<<pmtEvent->time[k]<<" volts "<<v[k]<<" "<<k<<", peakTime "<<peakTime[k-klow]<<endl;
        }
        //for(Int_t i = 0; i < peakTime.size(); i++) cout<<"\t"<<peakTime[i]<<endl;
      }
      // skip to end of sthreshold search 
      ibin=khigh + minLength;
    }
  }
   
  return peakTime;
}

Int_t anaRun::findHits( std::vector<Int_t> peakTime, std::vector<Double_t> ddigi) 
{

  if(peakTime.size()<1) return 0;
  std::vector<Int_t> hitTime;
  std::vector<std::vector<Int_t> > hitList;
  UInt_t nlast = peakTime.size()-1;
  for(UInt_t it=nlast; it>0; --it) {
    bool makeHit=false;
    //if at the end of the hit or end of the peakTime and hit is large enough
    if(peakTime[it]-peakTime[it-1]!=1 ||(it==1&&hitTime.size()>=minLength)) makeHit=true;
    hitTime.push_back(peakTime[it]);
//cout<<peakTime[it]<<endl;
    if(makeHit) {
      //cout<<"hit made "<<peakTime[it]<<" "<<peakTime[it-1]<<endl;
      //hit list has all time windows for each hit in an event
      hitList.push_back(hitTime);
      if(hitTime.size()<minLength) printf(" WARNING:: saving list %zu size %zu \n",hitList.size(),hitTime.size());
      //for(int r = 0; r < hitTime.size() ; r++) cout<<"\t"<<hitTime[r]<<endl;
      hitTime.clear();
      continue;
    }
    //hitTime.push_back(peakTime[it]);
  }

  Int_t nhits=0;
  //Double_t qhitMax = 0;
  qhitMax.clear();
  peakMaxTime.clear();
  peakBin.clear();
  qSum.clear();
  startTime.clear();
  peakWidth.clear();
  qhitMax.resize(hitList.size() );
  peakBin.resize(hitList.size() );
  peakMaxTime.resize(hitList.size() );
  qSum.resize(hitList.size() );
  startTime.resize(hitList.size() );
  peakWidth.resize(hitList.size() );

  for(UInt_t il=0; il<hitList.size(); ++il) {
    hitTime=hitList[il];
    if(hitTime.size()*deltaT > 2*forwardLengthTimeThreshold){
  //    cout<<"Found peak with width "<<hitTime.size()*deltaT<<", max peak size is "<<2*forwardLengthTimeThreshold<<endl;
      continue;
    }
   // printf(" %i hitTime.Size %i \n ",il,hitTime.size());
    Double_t qhit=0;
    UInt_t peakt=0;
    Double_t qpeak=0;
    //Double_t qUnhit=0;
    Double_t qUnpeak=0;
    for(UInt_t ih=0; ih<hitTime.size(); ++ih) {
//      printf(" \t ih = %i time  %i sample %f  \n ",ih,hitTime[ih],ddigi[hitTime[ih]]);
      if(ddigi[hitTime[ih]]>qpeak) {
        peakt=hitTime[ih];
        qpeak = ddigi[hitTime[ih]];
      }
      qSum[il]+=ddigi[hitTime[ih]];
    }
    qhitMax[il] = qpeak;
    peakMaxTime[il] = pmtEvent->time[peakt];
    peakBin[il] = peakt;

    startTime[il] = pmtEvent->time[hitTime[hitTime.size() - 1] ];
    peakWidth[il] = pmtEvent->time[hitTime[0] ] - startTime[il];
    ++nhits;
  }

  //qhitMax = BubbleSort(qhitMax);
  //peakMaxTime = BubbleSort(peakMaxTime);

  return  hitList.size();
}
std::vector<Double_t> anaRun::BubbleSort(std::vector<Double_t> A){
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

std::vector<Double_t> anaRun::SimpleLowPassFilter(std::vector<Double_t> signal, Double_t alpha){
  std::vector<Double_t> FilteredSignal;
  //Double_t alpha = 0.992105;//100 MHz cut off frequency
  for(int i = 0; i < signal.size(); i++){
    if(i == 0) FilteredSignal.push_back(alpha*signal[i]);
    else{
      FilteredSignal.push_back(alpha*signal[i] + (1.-alpha)*FilteredSignal[i-1] );
    }
  }
    
  //cout<<"Filter Size "<<FilteredSignal.size()<<", Ddigi size "<<signal.size()<<endl;
  return FilteredSignal;
}

std::vector<Double_t> anaRun::SimpleHighPassFilter(std::vector<Double_t> signal, Double_t alpha){
  std::vector<Double_t> FilteredSignal;
  //Double_t alpha = 0.992105;//100 MHz cut off frequency
  for(int i = 0; i < signal.size(); i++){
    if(i == 0) FilteredSignal.push_back(signal[i]);
    else{
      FilteredSignal.push_back( alpha*( signal[i]-signal[i-1] + FilteredSignal[i-1] )   );
    }
  }
    
  //cout<<"Filter Size "<<FilteredSignal.size()<<", Ddigi size "<<signal.size()<<endl;
  return FilteredSignal;

}

std::vector<Double_t> anaRun::MovingAverageFilter(std::vector<Double_t> signal,Int_t N)
{
  std::vector<Double_t> filter;
  Int_t N2 = std::floor(N/2);
  for(int i = N2; i < unsigned(signal.size())-N2; i++){
    Double_t sum = 0;
    for(int j = i-N2; j <= i+N2; j++){
      sum += signal[j];
    }
    sum /= N;
    filter.push_back(sum);
  }
  
  for(int i = 0; i < N2 ; i++){
    std::vector<Double_t>::iterator it;
    it = filter.begin();
    filter.insert(it,0.);
    filter.push_back(0);
  }
  return filter;
}



