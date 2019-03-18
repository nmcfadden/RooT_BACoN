////////////////////////////////////////////////////////
#include "anaRun.hh"


anaRun::anaRun(Int_t irun)
{

  Double_t irunStart = 38, irunStop = irunStart,irunName = irunStart + .01*irunStop;

  TString outFileName ; outFileName.Form("baconRunAna%f.root",irunName);
  TFile *outfile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());
  
  ntupleRun = new TNtuple("ntupleRun","ntupleRun","irun:ientry:pmtNum:sigma:nhits:charge:startTime:peakWidth:qMax:qMaxTime:nphotons:T0:Sdev:sigmaP:sigmaQ");
  ntupleEvent = new TNtuple("ntupleEvent","ntupleEvent","irun:ientry:pmtNum:time:volt");
  ntupleCal = new TNtuple("ntupleCal","ntupleCal","irun:ientry:pmtNum:tStart:tStop:peakMax:charge:zeroNoise:baseline:Sdev:Chi");


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
    for (ientry=0; ientry<nentries; ientry++) {
      pmtTree->GetEntry(ientry);
      if(irun == 32)   if(ientry < 593) continue;
//      if(ientry > 81) break;
  //    if(ientry != 81) continue;// || ientry != 45 || ientry != 69 || ientry != 81) continue;
      if(pmtEvent->time.size() == 0) continue;
      NEvents = pmtEvent->time.size();


      Double_t nameDouble = double(ientry)+.01*double(irun);
      TString name; name.Form("%f", nameDouble );
      
      //hPMTSum1 = new TH1D("PMT_Sum_0" + name,"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      //hPMTSum2 = new TH1D("PMT_Sum_1" + name,"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);

      
      //Be careful using Histograms, they get deleted after 1000 events!
      TH1D * hPMTSignal1 = new TH1D("PMTSignal1_" + name + TString("_baseline_subtracted"),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      TH1D * hPMTSignal2 = new TH1D("PMTSignal2_" + name + TString("_baseline_subtracted"),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      //TH1D * hPMTSignalFiltered1 = new TH1D("PMTSignal1_" + name + TString("_Filtered"),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);
      //TH1D * hPMTSignalFiltered2 = new TH1D("PMTSignal2_" + name + TString("_Filtered"),"",NEvents,pmtEvent->time[0],pmtEvent->time[NEvents-1]);

            
      TH1D * hBaseline1 = new TH1D("Baseline1_"+ name,"",NEvents,-50e-3,50e-3);
      TH1D * hBaseline2 = new TH1D("Baseline2_"+ name,"",NEvents,-50e-3,50e-3);
      
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

        
      //For cal ntuple
      //Int_t tStart = hPMTSignal1->FindBin(-50e-9);
      //Int_t tStop = hPMTSignal1->FindBin(50e-9);
      
      Int_t tStart = 0;//hPMTSignal1->FindBin(0e-9);
      Int_t tStop = NEvents;//hPMTSignal1->FindBin(500e-9);
      Int_t tStartZeroWidth = hPMTSignal1->FindBin(0e-9);
      Int_t tStopZeroWidth = hPMTSignal1->FindBin(50e-9);
      
      if(ientry == 0)  cout<<"StartTime = "<<hPMTSignal1->GetBinCenter(tStart)<<", StopTime = "<<hPMTSignal1->GetBinCenter(tStop)<<endl;

      
      std::vector<std::vector<Double_t> > ddigi;ddigi.resize(movingN.size() );
      Double_t vMin1 = 9999, vMin2 = 9999,vMax1 = -9999,vMax2 = -9999;

      for(Int_t ievent = 0; ievent < NEvents; ievent++){
        Double_t volt1 = pmtEvent->volt1[ievent];
        Double_t volt2 = pmtEvent->volt2[ievent];
        double digi1 = -1.0*(double(volt1));
        double digi2 = -1.0*(double(volt2));
        ddigi[0].push_back(digi1);
        ddigi[1].push_back(digi2);
       
        if(ievent > hPMTSignal1->FindBin(-500e-9) && ievent < hPMTSignal1->FindBin(-32e-9) ){
            hBaseline1->Fill(digi1);
            hBaseline2->Fill(digi2);
        }
        /*
        if(ievent > hPMTSignal1->FindBin(1e-6) && ievent < hPMTSignal1->FindBin(2.5e-6) ){
            hBaseline1->Fill(digi1);
            hBaseline2->Fill(digi2);
        } 
        */
        if(vMin1 > digi1) vMin1 = digi1;
        if(vMin2 > digi2) vMin2 = digi2;
        if(vMax1 < digi1) vMax1 = digi1;
        if(vMax2 < digi2) vMax2 = digi2;

        if(ievent > 0)
          deltaT += std::fabs(pmtEvent->time[ievent] - pmtEvent->time[ievent - 1]);
      }
      
      deltaT /= (NEvents - 1);
      //minLength = int(backwardLengthTimeThreshold/deltaT);
      backwardLengthTimeThreshold = backwardLengthTimeThresholdDouble/deltaT;
      minLength = int(minLengthDouble/deltaT);
      forwardLengthTimeThreshold = int(forwardLengthTimeThresholdDouble/deltaT);

      if(ientry == 0 ){ 
        cout<<"deltaT = "<<deltaT<<", Minimum Peak width = "
          <<minLength<<", backwardLengthTime "<<backwardLengthTimeThreshold<<", forwardLengthTimeThreshold "<<forwardLengthTimeThreshold<<endl;
      }

      std::vector<Double_t> baseline; baseline.resize(ddigi.size() );
      std::vector<Double_t> Chi;Chi.resize(ddigi.size() );
      std::vector<Double_t> Sdev;Sdev.resize( ddigi.size() );
     
      TF1 * myfit1 = new TF1("myfit1","gaus(0)",-50e-3,50e-3);
      //Q is quite mode! very important
      hBaseline1->Fit("myfit1","Q","",vMin1,vMax1);
      baseline[0] = myfit1->GetParameter(1);
      Chi[0] = myfit1->GetChisquare()/myfit1->GetNDF();
      Sdev[0] = myfit1->GetParameter(2);
      delete hBaseline1;
      if(Chi[0] > 15 || std::fabs(baseline[0]) > 3.e-3){
        delete hPMTSignal1;
        delete hPMTSignal2;
        delete hBaseline2;
        continue; 
      }
      TF1 * myfit2 = new TF1("myfit2","gaus(0)",-50e-3,50e-3);
      hBaseline2->Fit("myfit2","Q","",vMin2,vMax2);
      //if(!hBaseline2->GetFunction("gaus") ) continue;
      baseline[1] = myfit2->GetParameter(1);
      Chi[1] = myfit2->GetChisquare()/myfit2->GetNDF();
      Sdev[1] = myfit2->GetParameter(2);
      delete hBaseline2;
      if(Chi[1]>15 || std::fabs(baseline[1]) > 3.e-3 ){
        delete hPMTSignal1;
        delete hPMTSignal2;
        continue; 
      }
      for(int j = 0; j < ddigi.size();j++){
        for(int i = 0; i < ddigi[j].size(); i++){
          ddigi[j][i] = ddigi[j][i]-baseline[j];// + Sdev[j];
          //ntupleEvent->Fill(irun,ientry,j,pmtEvent->time[i],ddigi[j][i]);
          Double_t time = pmtEvent->time[i];
          Int_t timeBin = hPMTSignal1->FindBin(time);
          if(j == 0) {
            hPMTSignal1->SetBinContent(timeBin,ddigi[j][i]);
          }
          else{
            hPMTSignal2->SetBinContent(timeBin,ddigi[j][i]);
          }

        }
      }
      

      std::vector<Double_t> tempPmtNum,tempSigma,tempNhits,tempQsum,tempStartTime,tempPeakWidth,tempqHit,tempPeakMaxTime,tempNphotons;
      T0.clear();
      Double_t T0_temp = 1;
      std::vector<bool> nohitFlag;nohitFlag.resize(2);
      for(int pmtNum = 0 ; pmtNum < ddigi.size();pmtNum++){
        for(double sigma = 3 ; sigma < 4; sigma++){
          Double_t minDev = 1*Sdev[pmtNum], maxDev = sigma*Sdev[pmtNum];
          std::vector<Int_t> peakTime = findPeaks(ddigi[pmtNum],maxDev,minDev,tStart,tStop);
          Int_t nhits = findHits( peakTime, ddigi[pmtNum],maxDev);
          
          if(nhits == 0) nohitFlag[pmtNum] = true;
          else nohitFlag[pmtNum] = false;

          if(pmtNum == 1 && 0){
            cout<<"\t\t"<<sigma<<"-SIGMA THRESHOLD == "<<maxDev<<endl;
            for(int v = 0; v < peakTime.size() ; v++){
              if(v > 0)
                if( (peakTime[v] - peakTime[v-1]) != 1) cout<<"\t new peak found"<<endl;
              cout<<ientry<<" = ientry, "<<sigma*Sdev[pmtNum]<<" = "<<nhits<<"-nhits, "<<pmtEvent->time[peakTime[v]]<<" = time,"<<ddigi[pmtNum][peakTime[v] ]<<" = volts"<<endl;
            }
          }
          if(pmtNum == 1 && 0){
            for(int i = 0; i < nhits; i++)
              cout<<"ientry "<<ientry<<", qMax "<<qhitMax[i]<<", time "<<peakMaxTime[i]<<", start "<< startTime[i]<<", end "<<startTime[i]+peakWidth[i] <<endl;
          }
          Double_t nphotons = 0;bool writeFlag = false;
          for(int i = 0; i < nhits; i++) nphotons += qSum[i]*deltaT/1.02e-10;
          bool chargeFlag = false; int chargeBin = 0;
          for(int i = 0; i < nhits; i++){
            /*
            qSum[i] = 0;
            for(int j = peakBin[i] - 16e-9/deltaT; j <= peakBin[i] + 24e-9/deltaT; j++){
              if( j < 0 ) continue;
              if( j > (ddigi[pmtNum].size()-1) ) continue;
              qSum[i] += ddigi[pmtNum][j];
            }
            */
            //ntupleRun->Fill(irun,ientry,pmtNum,sigma,nhits,qSum[i]*deltaT, startTime[i] ,peakWidth[i],qhitMax[i],peakMaxTime[i],nphotons);
            tempPmtNum.push_back(pmtNum);tempSigma.push_back(sigma);tempQsum.push_back(qSum[i]);tempStartTime.push_back(startTime[i]);
            tempPeakWidth.push_back(peakWidth[i]);tempqHit.push_back(qhitMax[i]);tempPeakMaxTime.push_back(peakMaxTime[i]);
            if(qSum[i] < 0 ){
              chargeFlag = true;
              chargeBin = i;
            }
            if (1. < (tempPeakWidth[i]*Sdev[tempPmtNum[i]]) / sqrt(tempQsum[i]*deltaT/10.4e9) ) writeFlag = true;
          }
          if(pmtNum == 0){
            tempNhits.push_back(nhits);
            tempNphotons.push_back(nphotons);
          }
          else{
            tempNhits[0] += nhits;
            tempNphotons[0] += nphotons;
          }

          if(qSum.size() > 0 && sigma == 3){
          //if(chargeFlag){
            if(pmtNum == 0)
              hPMTSignal1->SetTitle(TString("PMT_")+to_string(pmtNum)+TString(".")+to_string(ientry)+TString("_charge_")+to_string(deltaT*qSum[chargeBin]*1e10)+TString("_peak_")+to_string(qhitMax[chargeBin])+
                TString("peakBin_")+to_string(peakBin[chargeBin])+TString("_start_")+to_string(startTime[chargeBin]*1e6)+TString("_stop_")+to_string((startTime[chargeBin]+peakWidth[chargeBin])*1e6 ) + TString("_nhits_")+to_string(nhits) );
            else
              hPMTSignal2->SetTitle(TString("PMT_")+to_string(pmtNum)+TString(".")+to_string(ientry)+TString("_charge_")+to_string(deltaT*qSum[chargeBin]*1e10)+TString("_peak_")+to_string(qhitMax[chargeBin]) +
                TString("peakBin_")+to_string(peakBin[chargeBin])+TString("_start_")+to_string(startTime[chargeBin]*1e6)+TString("_stop_")+to_string( (startTime[chargeBin]+peakWidth[chargeBin])*1e6 ) + TString("_nhits_")+to_string(nhits) );
            singlePhoton++;
          }
          else if(ientry < 1000 && sigma == 3){
            if(pmtNum == 0)
              delete hPMTSignal1;
            else
              delete hPMTSignal2;
          }
          qhitMax.clear();
          peakMaxTime.clear();
          peakBin.clear();
          qSum.clear();
          startTime.clear();
          peakWidth.clear();
        }
        for(int i = 0; i < tempPmtNum.size() ; i++){
          ntupleRun->Fill(irun,ientry,tempPmtNum[i],tempSigma[i],tempNhits[0],tempQsum[i]*deltaT, tempStartTime[i] ,tempPeakWidth[i],tempqHit[i],tempPeakMaxTime[i],tempNphotons[0],BubbleSort(tempStartTime)[0],tempPeakWidth[i]*Sdev[tempPmtNum[i]],sqrt(tempQsum[i]*deltaT/10.4e9),sqrt(tempQsum[i]*deltaT/10.4e9 + std::pow(Sdev[tempPmtNum[i]]*tempPeakWidth[i],2) ));
        }
        tempPmtNum.clear();tempSigma.clear();tempNhits.clear();tempQsum.clear();tempStartTime.clear();tempPeakWidth.clear();tempqHit.clear();tempPeakMaxTime.clear();tempNphotons.clear();T0.clear();
      }
     
      
      //cout<<"Start bin = "<<tStart<<", Start time = "<< pmtEvent->time[tStart]<<", Stop bin = "<<tStop<<", Stop time = "<< pmtEvent->time[tStop]<<endl;
      Double_t charge = 0,zeroNoise = 0,peakVal = -999;
        Int_t peakBin = 0; 
      for(int j = 0; j <2 ; j++){
        
        for(int i = tStartZeroWidth ; i < tStopZeroWidth ; i++){
          charge += ddigi[j][i]*deltaT;
        }

        for(int i = 0 ; i < 40e-9/deltaT ; i++){
          zeroNoise += ddigi[j][i]*deltaT;
        }
        
          ntupleCal->Fill(irun,ientry,j,tStart,tStop,peakVal,charge,zeroNoise,baseline[j],Sdev[j],Chi[j]);
          peakVal = -999;
          charge = 0;
          zeroNoise = 0;
      }
      

      ///*
      if(ientry >= 1000 ){
        delete hPMTSignal1;
        delete hPMTSignal2;
       // delete hPMTSignalFiltered1;
       // delete hPMTSignalFiltered2;
      }
      //*/
      //
      if(ientry%5000 == 0) cout<<ientry<<" Events Processed"<<endl;

     // if(ientry == 100) break;
      //if(ientry>25) break;
      //if(ientry > 50) break;
      //if(ientry>12364) break;//run_16 good data set
      //if(ientry>40000) break;
    }
    cout<<singlePhoton<<" events processed out of "<<nentries<<endl;
    cout<<"\t WARNING SYSTEM IS IN CRITICAL CONDITION! \n\t NEEDS MORE CHEESE FOR MOUSE!"<<endl;
  }//irun loop

  outfile->Write();

return;
}


std::vector<Int_t> anaRun::findPeaks(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold,Int_t startBin,Int_t stopBin) 
{
  // Produces a list of peaks above the threshold
  std::vector<Int_t> peakTime;
  Int_t klow=0;
  Int_t khigh=0;
  Int_t kover=0;
  Int_t vsize = Int_t(v.size())-1;
  Double_t fMin = 9999,bMin = 9999;
  Int_t klast = 0;
  
  //printf(" findPeaks \n");
  //for(Int_t  ibin=0; ibin<= vsize; ++ibin ) {
  for(Int_t  ibin=startBin; ibin< stopBin; ++ibin ) {
    if( v[ibin]>threshold ){//&& (v[ibin] > v[ibin - 1] ) && (v[ibin] < v[ibin + 1] ) ){// starting possible new hit and only find hits on a positive slope
      // consider this a "seed" and find full hit
      klow=ibin;
      for(Int_t k=ibin-1; k>=std::max(0,ibin-backwardLengthTimeThreshold); --k) {
        
        if(bMin > v[k]){
          //bMin = v[k];
          //klow = k;
        }
        
        if(v[k]<=sthreshold){
          klow=k; 
          break;
        }
        if(k == std::max(0,ibin-backwardLengthTimeThreshold) ){
        //  klow=k;
        }
        klow=k;
      }
      khigh=ibin;
      for(Int_t k=ibin+1; k<=min(ibin+forwardLengthTimeThreshold,vsize); ++k) {
        /*
        if(fMin > v[k]){
          fMin = v[k];
          khigh = k;
        }
        */
        if(v[k]<=sthreshold){
          khigh = k; 
          break;
        }
        if(k == min(ibin+forwardLengthTimeThreshold,vsize)){
          //khigh=k;
        }
        khigh=k;
      }
      kover = khigh-klow+1;
      if(klast < klow) klast = klow;
      // found good pulse
      if(kover>minLength) { 
        for(Int_t k=klast ; k<= khigh; ++k){
          peakTime.push_back(k);
          fMin = 9999;
          bMin = 9999;
          //cout<<"hi time "<<pmtEvent->time[k]<<" volts "<<v[k]<<" "<<k<<", peakTime "<<peakTime[k-klow]<<endl;
        }
        klast = khigh;
        //for(Int_t i = 0; i < peakTime.size(); i++) cout<<"\t"<<peakTime[i]<<endl;
      }
      // skip to end of sthreshold search 
      ibin=khigh;// + backwardLengthTimeThreshold;
    }
  }
   
  return peakTime;
}

Int_t anaRun::findHits( std::vector<Int_t> peakTime, std::vector<Double_t> ddigi,Double_t sigma) 
{

  if(peakTime.size()<1) return 0;
  std::vector<Int_t> hitTime;
  std::vector<std::vector<Int_t> > hitList;
  UInt_t nlast = peakTime.size()-1;
  for(Int_t it=nlast; it>0; --it) {
    if(ddigi[it] > .082) return 0;//for run 34 only
    bool makeHit=false;
    //if at the end of the hit or end of the peakTime and hit is large enough
    Double_t fDelta1 = 0,fDelta2 = 0,bDelta1 = 0,bDelta2 = 0;
    if( (it+1) < peakTime.size() && (peakTime[it+1]) < ddigi.size())
      fDelta1 = ddigi[peakTime[it+1]]-ddigi[peakTime[it]];

    if( ( (it+2) < peakTime.size()) &&(peakTime[it+2]) < ddigi.size())
      fDelta2 = ddigi[peakTime[it+2]]-ddigi[peakTime[it+1]];

    if( (it-1)>=0 && (peakTime[it-1]) >= 0 )
      bDelta1 = ddigi[peakTime[it]] -ddigi[peakTime[it-1]];

    if( (it-2)>=0 && (peakTime[it-2]) >= 0 )
      bDelta2 = ddigi[peakTime[it-1] ] - ddigi[peakTime[it-2] ];
    
    if(peakTime[it]-peakTime[it-1]!=1 ||(it==1&&hitTime.size()>=minLength)) makeHit=true;
    
    
    else if( (bDelta1 <= 0) && (bDelta2 < 0) && (fDelta1 > 0)  && (fDelta2 >= 0) &&hitTime.size()>=minLength ){
      makeHit=true;
      //cout<<pmtEvent->time[peakTime[it]]<<" bDelta2 "<<bDelta2<<" bDelta1 "<<bDelta1<<" fDelta1 "<<fDelta1<<endl;
    }
    
    else if( (bDelta1 < 0) && (bDelta2 <= 0) && (fDelta1 >= 0)  && (fDelta2 > 0) &&hitTime.size()>=minLength ){
      makeHit=true;
      //cout<<pmtEvent->time[peakTime[it]]<<" bDelta2 "<<bDelta2<<" bDelta1 "<<bDelta1<<" fDelta1 "<<fDelta1<<endl;
    }
    else if( (bDelta1 < 0) && (fDelta1 >= 0)  && (fDelta2 > 0) && -bDelta1 > sigma && hitTime.size()>=minLength ){
      makeHit=true;
      //cout<<pmtEvent->time[peakTime[it]]<<" bDelta2 "<<bDelta2<<" bDelta1 "<<bDelta1<<" fDelta1 "<<fDelta1<<endl;
    }
    else if( (bDelta1 <= 0) && (bDelta2 < 0)  && (fDelta1 > 0) && fDelta1 > sigma && hitTime.size()>=minLength ){
      makeHit=true;
      //cout<<pmtEvent->time[peakTime[it]]<<" bDelta2 "<<bDelta2<<" bDelta1 "<<bDelta1<<" fDelta1 "<<fDelta1<<endl;
    }
 
    
    hitTime.push_back(peakTime[it]);
//cout<<peakTime[it]<<endl;
    if(makeHit) {
      if(hitTime.size() <= minLength){
//        hitTime.pop_back(); 
        continue;
      }
      //hit list has all time windows for each hit in an event
      if(it == 1) hitTime.push_back(peakTime[it-1]);
      hitList.push_back(hitTime);
      //if(hitTime.size()<minLength) printf(" WARNING:: saving list %zu size %f \n",hitList.size(),deltaT*hitTime.size());
      //for(int r = 0; r < hitTime.size() ; r++) cout<<"\t"<<hitTime[r]<<endl;
      //cout<<"new hit"<<endl;
      hitTime.clear();
      //new hit starts a min
      if(peakTime[it] -peakTime[it-1] == 1) it++;
      continue;
    }
    //hitTime.push_back(peakTime[it]);
  }

  Int_t nhits=0;
 
  for(UInt_t il=0; il<hitList.size(); ++il) {
    hitTime=hitList[il];
    if(hitTime.size() > forwardLengthTimeThreshold){
      //cout<<"Found peak with width "<<hitTime.size()<<", max peak size is "<<forwardLengthTimeThreshold<<endl;
      continue;
    }
   // printf(" %i hitTime.Size %i \n ",il,hitTime.size());
    Double_t qhit=0;
    UInt_t peakt=0;
    Double_t qpeak=0,qsum = 0;
    Double_t qUnpeak=0;
    for(UInt_t ih=0; ih<hitTime.size(); ++ih) {
      //printf(" \t ih = %i time(us)  %f sample %f  \n ",hitTime[ih],1e6*pmtEvent->time[hitTime[ih] ],ddigi[hitTime[ih]]);
      if(ddigi[hitTime[ih]]>qpeak) {
        peakt=hitTime[ih];
        qpeak = ddigi[hitTime[ih]];
      }
      qsum += ddigi[hitTime[ih]];
    }
    if (qpeak < sigma){
      continue;
    }
    if(pmtEvent->time[hitTime[0] ] - pmtEvent->time[hitTime[hitTime.size() - 1] ] > forwardLengthTimeThreshold*deltaT) continue;
    qSum.push_back(qsum);
    qhitMax.push_back(qpeak);
    peakMaxTime.push_back(pmtEvent->time[peakt]);
    peakBin.push_back(peakt);
    startTime.push_back(pmtEvent->time[hitTime[hitTime.size() - 1] ] );
    peakWidth.push_back( pmtEvent->time[hitTime[0] ] - pmtEvent->time[hitTime[hitTime.size() - 1] ] );
    ++nhits;
  }


  return  nhits;
}

std::vector<std::complex<double> > anaRun::FFT(Int_t ipmt,Int_t ievent,std::vector<Double_t> signal)
{
  /*
     for(int i = 0; i < ddigi.size(); i++){
     std::vector<std::complex<double> > fftPair = FFT(i,ientry,ddigi[i]);
     std::vector<Double_t> vecTemp = inverseFFT(i,ientry,fftPair,sum);
     }
  */  
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



