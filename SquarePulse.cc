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

void SquarePulse(){
  TString outFileName = "SquarePulseOut.root";
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());

  TRandom2 rand;
  TDatime time;
  rand.SetSeed(time.GetTime());
  Double_t gausMean = 0,gausSigma = .1;
  Int_t Npulses = 1e5;
  Double_t deltaT = 4e-10;
  Int_t nBins =100,pulseWidth = 40, signalStart = 50, noiseStart = 10;
  

  TH1D * hSig = new TH1D("SignalArea","SignalArea",nBins,35,45);
  TH1D * hNoise = new TH1D("NoiseArea","NoiseArea",nBins,-20*gausSigma,20*gausSigma);

  for(int i =0; i < Npulses ;i++){
    TH1D * hPulse = new TH1D(TString("SquarePulse")+to_string(i),TString("SquarePulse")+to_string(i),nBins,0,nBins*deltaT);
    Double_t sumSig = 0, sumNoise = 0;
    for(int j = 0; j < nBins;j++){
      Double_t noise = gausMean+gausSigma*sqrt(-2.0*log(rand.Rndm()))*cos(2*TMath::Pi()*rand.Rndm()) ;
      if(j >=signalStart && j< signalStart+pulseWidth){
        hPulse->SetBinContent(j,1);
      }
      hPulse->SetBinContent(j,hPulse->GetBinContent(j)+noise);
      if(j >= noiseStart && j < noiseStart + pulseWidth) sumNoise += hPulse->GetBinContent(j);
      if(j >= signalStart && j < signalStart + pulseWidth) sumSig += hPulse->GetBinContent(j);
    }
    hSig->Fill(sumSig);
    hNoise->Fill(sumNoise);
    if(i > 100) delete hPulse;
    //hPulse->Draw();
  }
}
