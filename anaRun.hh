//////////////////////////////////////////////////////////
//  M.Gold August 2017
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>//includes std::pair, std::make_pair
#include <valarray>
//
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm>    // std::sort
#include "TSpectrum.h"
//
#include "TPmtEvent.hxx"

typedef std::complex<double> Complex;


class anaRun {
public :
  enum {NPMT=2};
  enum {MAXSAMPLES=10000};
  //pretty sure it is the min/max number of bins 
  //enum {minLength=100,maxHalfLength=400};
  Int_t minLength;
  UInt_t ientry;
  Double_t minLengthDouble = 24e-9;
  Double_t backwardLengthTimeThresholdDouble = 32e-9;
  Int_t    backwardLengthTimeThreshold;
  Int_t    forwardLengthTimeThreshold;
  Double_t forwardLengthTimeThresholdDouble = 100e-9; //halfMax time window (looks N seconds left and right)
  Double_t deltaT = 0,deltaV = 999,singlePhoton = 0;
  std::vector<Double_t> startTime;
  std::vector<Double_t> peakWidth;
  std::vector<Double_t> qhitMax;
  std::vector<Double_t> peakMaxTime;
  std::vector<Double_t> peakBin;
  std::vector<Double_t> qSum;
  std::vector<Double_t> T0;

 
  anaRun(Int_t irun=0);
  virtual ~anaRun(){;}

  TTree* pmtTree;
  TPmtEvent* pmtEvent;
  

  TSpectrum * spec = new TSpectrum(20000);
  //std::vector<Int_t> findMaxPeak(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); 
  //std::vector<Int_t> findPeaks(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); 
  //Int_t findHits(Int_t ipmt, Double_t sum,  std::vector<Int_t> peakTime, std::vector<Double_t> ddigi); 

  std::vector<Int_t> findPeaks(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold,Int_t startBin, Int_t stopBin);
  Int_t findHits( std::vector<Int_t> peakTime, std::vector<Double_t> ddigi,Double_t sigma);
  std::vector<Double_t> SimpleFilter(std::vector<Double_t> in);
  std::vector<Double_t> SimpleLowPassFilter(std::vector<Double_t> in, Double_t a);
  std::vector<Double_t> SimpleHighPassFilter(std::vector<Double_t> in, Double_t a);
  std::vector<Double_t> MovingAverageFilter(std::vector<Double_t> signal,Int_t N);
  std::vector<Double_t> BubbleSort(std::vector<Double_t> A);

  //Double_t convolve(Double_t *x, Double_t *par);
  Int_t NEvents;
  /// The fft class to take the fourier transform.
  TVirtualFFT *fFFT;
  /// The fft class to take the inverse fourier transform.
  TVirtualFFT *fInverseFFT;
  //std::vector<Double_t> FFTFilter(Int_t ipmt,Int_t ievent,Double_t norm);
  //TH1D* FFTFilter(Int_t ipmt,Int_t ievent);
  //std::pair<TH1D*,std::vector<Double_t> > FFTFilter(Int_t ipmt,Int_t ievent,Double_t norm,std::vector<Double_t> signal);
  std::vector<std::complex<double> > FFT(Int_t ipmt,Int_t ievent,std::vector<Double_t> signal);
  std::vector<Double_t > inverseFFT(Int_t ipmt,Int_t ievent, std::vector<std::complex<double> > VectorComplex,std::vector<Double_t> sum);
  // histogram pointers
  TNtuple *ntupleRun;
  TNtuple *ntupleEvent;
  TNtuple *ntupleCal;
  TH1D* hSamples[NPMT];  // include RF
  TH1D* hFFT[NPMT];
  TH1D* hHitQ[NPMT];
  TH1D* hNHits[NPMT];

  TH1D* hPMTSum1;
  TH1D* hPMTSum2;
  
  TH1D* hNoise;
  TH1D* hFFTNoise;
  TH1D* hBase;
  TH1D* hChargeWindow;
  TH1D* hChargePeak;
  TH1D* hPeakWidth;
  TH1D* hBaseline;
  TH1D* hBaselineFFT;
  TH1D* hfft;
  TH1D* hIfft;
  TH1D* hHits;
  TH1D* hDeltaV;
  TH1D* hStartTime;
  TH1D* hQMax;
  TH1D* hChargePeakOff;
};

