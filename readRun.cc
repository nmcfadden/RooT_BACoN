////////////////////////////////////////////////////////
#include<stdio.h>
#include <string.h>
#include <dirent.h> 
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
// local 
#include "TPmtEvent.hxx"
TTree *pmtTree;
TPmtEvent *pmtEvent;

Int_t eventCounter = 0;
TString date,localTime;
Double_t elapseTime;

void getEventsInDirectory(std::string directory, std::vector<Int_t> &eventNumbers)
{
  DIR *dir;
  struct dirent *dp;
  //class stat st;
    string tag("event");
    printf(" looking in directory %s \n",directory.c_str());
    dir = opendir(directory.c_str());
    while ((dp = readdir(dir)) != NULL) {
        string file_name = dp->d_name;
        if ( strstr(file_name.c_str(), "event" )==NULL ) continue;
        if ( strstr(file_name.c_str(), ".txt" )== NULL ) continue;
        TString fname(file_name.c_str());
        Int_t ievent = TString(fname( fname.First("_")+1, fname.First(".") - fname.First("_"))).Atoi();
        eventNumbers.push_back(ievent);
    }
    closedir(dir);
    eventNumbers.pop_back();
    printf(" found %i event files \n",Int_t(eventNumbers.size()));
} // GetFilesInDirectory

int readEvent(Int_t ievent, TString fileName, Int_t irun){
  //cout << " reading file " << fileName << endl;

  //pmtEvent->clear();
  //pmtEvent->event=ievent;
  //printf(" looking for file %s \n",fileName.Data());

  Int_t nlines = 0;
  Double_t time,volt1,volt2;
  std::vector<Double_t> timeVec,volt1Vec,volt2Vec;
  ifstream in;
  in.open(fileName);
  if(!in.is_open() ) {
    printf(" cannot open file %s \n",fileName.Data());
    return 0;
  }
  while (1) {
    if(nlines == 0){
      in >> date >> localTime >> elapseTime ;
    }
    else{
      if(irun < 5)
        in >> time >> volt1 >> volt2 ;
      else
        in >> time >> volt1 ;
      if(!in.good()) break;
      timeVec.push_back(time);
      volt1Vec.push_back(volt1);
      //volt2Vec.push_back(volt2);
      /*
      pmtEvent->time.push_back(time);
      pmtEvent->volt1.push_back(volt1);
      pmtEvent->volt2.push_back(volt2);
      */
    }
    nlines++;
  }
  //printf(" have read %i lines and %i entries \n",nlines,pmtEvent->time.size());
  in.close();
  
  //cout<<date<<" "<<localTime<<" "<<elapseTime<<endl;
  if(ievent != 0 && 
      pmtEvent->volt1[0] == volt1Vec[0]     && pmtEvent->volt1[10] == volt1Vec[10]     && 
      pmtEvent->volt1[100] == volt1Vec[100] && pmtEvent->volt1[110] == volt1Vec[110] &&
      pmtEvent->volt1[150] == volt1Vec[150] && pmtEvent->volt1[160] == volt1Vec[160] &&
      pmtEvent->volt1[200] == volt1Vec[200] && pmtEvent->volt1[210] == volt1Vec[210] 
    ){
    //if(pmtEvent->volt1[0] == volt1Vec[0]) cout<<"same volts1, "<<ievent<<endl;
    //if(pmtEvent->volt2[0] == volt2Vec[0]) cout<<"same volts2, "<<ievent<<endl;
    pmtEvent->volt1 = volt1Vec;
    //pmtEvent->volt2 = volt2Vec;
    pmtEvent->time = timeVec;
  }
  else{
    pmtEvent->volt1 = volt1Vec;
    //pmtEvent->volt2 = volt2Vec;
    pmtEvent->time = timeVec;
    eventCounter++;
    pmtTree->Fill();
  }
  pmtEvent->event=eventCounter;
  //pmtEvent->event=ievent;
  pmtEvent->clear();

  return nlines;
    
}

void readRun(Int_t irun)
{
  //for(irun  = 170; irun <=185;irun++){
    if(irun == 26) return;
  // open ouput file and make some histograms
  //TString outFileName; outFileName.Form("rootData/LAr_%i.root",irun);
  //TString outFileName; outFileName.Form("rootData/recirRoot_%i.root",irun);
  //TString outFileName = TString("cosmic_1_11.27.2018_10000_events.root");
  //TString outFileName = TString("rootData/run_6.root");
  TString outFileName; 
  if(irun >= 1000 && irun <= 1002)
    outFileName.Form("rootData/baconRun_10kBins_1us_20mV_muon_perpendicular_%i.root",irun);
  else if(irun < 2000 && irun >1002)
    outFileName.Form("rootData/baconRun_10kBins_1us_20mV_muon_parallel_%i.root",irun);
  else if(irun >= 4000 && irun <5000)
    outFileName.Form("rootData/baconRun_10kBins_1us_20mV_div_-6mV_thresh_%i.root",irun);
  else if(irun >=5000 && irun < 6000)
    outFileName.Form("rootData/baconRun_10kBins_1us_20mV_div_-6mV_thresh_%i.root",irun);
  else if(irun >= 10000 )
    outFileName.Form("rootData/baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_%i.root",irun);
  TFile *outfile = new TFile(outFileName,"recreate");
  outfile->cd();
  printf(" opening output file %s \n",outFileName.Data());

  // ttree
  pmtTree = new TTree("pmtTree","pmtTree");
  pmtEvent  = new TPmtEvent();
  pmtTree->Branch("pmtEvent",&pmtEvent);

  // get list of files
  TString dirName;
  dirName.Form("rawData/run_%i",irun);
  //dirName = TString("rawData/cosmic_1_11.27.2018_10000_events");
  //dirName = TString("rawData/run_6");
  //dirName = TString("rawData/background");
  std::string directory(dirName.Data());

  std::vector<Int_t> eventNumbers;
  getEventsInDirectory(directory,eventNumbers);
  Int_t nlines=0;
  for( unsigned ievent =0; ievent < eventNumbers.size() ; ++ievent ) {
    TString fname;
    fname.Form("event_%i.txt",ievent);
//    fname.Form("backup_%i.txt",ievent);
    //int ievent = TString(fname( fname.First("_")+1, fname.First(".") - fname.First("_"))).Atoi();
    TString fullFileName = dirName + TString("/")+fname;
    if(ievent%100==0) cout << ievent << " tree size  " << pmtTree->GetEntries()  << endl;
    nlines += readEvent(ievent,fullFileName,irun);
  }
  TString eventInfoString = date + TString("-") + localTime + TString("_runTime_") + to_string(elapseTime);
  TNamed eventInfo(TString("eventInfo"),eventInfoString);
  eventInfo.Write();
  printf(" total of lines %i total number of events is %i \n",nlines,int(pmtTree->GetEntries()));
  outfile->Write();
  //}
}
