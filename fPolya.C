// test polya distributions
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TKey.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"

using namespace TMath;
TString htag;
std::vector<TH1D*> hlist;
static double xmax=20;
static double xmin=0;
static double nsum;
// single polya
static double fpolya(double *x, double *par)
{
  double y = x[0]/par[1];
  double e = par[0];
  double norm = Gamma(e,xmax/par[1]*e)-Gamma(e,xmin/par[1]*e);   
  double f = par[2]*e*pow(e*y,e-1)*Exp(-e*y)/Gamma(e)/par[1]/norm;
  return f; 
}
// polya + exp 
static double fpolyaE(double *x, double *par)
{
  double y = x[0]/par[1];
  double e = par[0];
  double g = Gamma(par[0]);
  double A = 1. - xmin/par[2];
  double B = 1. - xmax/par[2];
  double norm = Gamma(e,xmax/par[1]*e)-Gamma(e,xmin/par[1]*e);
  // prevent numerical errors
  if(xmin/par[2]>1E-5) Exp(-xmin/par[2]);
  if(xmax/par[2]>1E-5) Exp(-xmax/par[2]);
  double norme = A-B;
  double C = 1.0-e*y;
  if(e*y>1E-5) C=Exp(-e*y);
  double sig  = par[3]*e*pow(e*y,e-1)*C/g/par[1]/norm;
  double bg = (1-par[3])*exp(-x[0]/par[2])/par[2]/norme;
  return par[4]*(sig+bg); 
}

/*****************************************************
   read all histograms 
*******************************************************/

void reading(TDirectory *fdir) 
{
  TObject *obj=NULL;
  TKey* key;
  TIter nextkey(fdir->GetListOfKeys());
  //htag = TString("QhitNoBeam");
  htag = TString("QhitOff");
  while ( (key = (TKey*)nextkey()) )  {
    fdir->cd();
    obj = key->ReadObj();
    if(obj->IsA()->InheritsFrom("TH1D")) { //case of TH1 or TProfile
      TH1D* h1 = dynamic_cast<TH1D*>(obj);
      if( TString(h1->GetName()).Contains(htag) ) hlist.push_back(h1);
    } 
  }
}
/*****8*********
  read ntuple
******8*********/
void readNtuple(Double_t tag){
  tag = tag +.01*tag; 
  TString filename; filename.Form("baconRunAna%f.root",tag);
  TFile *f = new TFile(filename);
  if(f->IsZombie()){
    return;
  }

  for(int i = 0; i < hlist.size();i++){
    hlist[i] = new TH1D(filename+TString("_pmt_")+to_string(i),filename+TString("_pmt_")+to_string(i),200,-1e-10,5e-10);
  }

  TTree *t1 = (TTree*)f->Get("ntupleRun");

  Float_t irun,ientry,pmtNum,sigma,nhits,charge,baseline,Sdev,startTime,qMax,qMaxTime,T0; 

  t1->SetBranchAddress("irun",&irun);
  t1->SetBranchAddress("ientry",&ientry);
  t1->SetBranchAddress("pmtNum",&pmtNum);
  t1->SetBranchAddress("sigma",&sigma);
  t1->SetBranchAddress("nhits",&nhits);
  t1->SetBranchAddress("charge",&charge);
  t1->SetBranchAddress("startTime",&startTime);
  t1->SetBranchAddress("qMax",&qMax);
  t1->SetBranchAddress("qMaxTime",&qMaxTime);
  t1->SetBranchAddress("T0",&T0);
    
  Int_t NEntries = (Int_t)t1->GetEntries();

  for(int i = 0; i < NEntries; i++){
    t1->GetEntry(i);
    if(sigma != 3) continue;
    hlist[pmtNum]->Fill(charge);
  }
}

void fPolya(Double_t tag=34)
{
  enum {NPMT=2};

  hlist.resize(NPMT);
  readNtuple(tag);
  for(UInt_t ih=0; ih<hlist.size(); ++ih) printf(" %i %s \n",ih,hlist[ih]->GetName());

  float en=1e-10;
  float qn=1e6;

  TF1 *fp[NPMT];
  for(int i=0; i<NPMT; ++i ) {
    fp[i] = new TF1(Form("poyaE-%i",i),fpolyaE,xmin,xmax,5);
    fp[i]->SetNpx(1000); // numb points for function
    fp[i]->SetParName(0,"primary charge en");
    fp[i]->SetParName(1,"gain Qn");
    fp[i]->SetParName(2,"exp const");
    fp[i]->SetParName(3,"sig frac");
    fp[i]->SetParName(4,"norm");
    fp[i]->SetParameter(0,en);
    fp[i]->SetParameter(1,qn);
    fp[i]->SetParameter(2,10);
    fp[i]->FixParameter(3,1.0);
    fp[i]->SetParameter(4,1.0);
    fp[i]->SetLineColor(kRed);
  }
  TCanvas *cpolyaE = new TCanvas("polyaE","polyaE");
  fp[0]->Print();
  TCanvas *cpolya = new TCanvas("polya","polya");
  //fp[0]->Draw();
  double epar=fp[0]->GetParameter(0); 
  double qpar=fp[0]->GetParameter(1); 
  double pnorm = Gamma(epar,xmax/qpar*epar)-Gamma(epar,xmin/qpar*epar); // from lower limit of xmin
  double econst=fp[0]->GetParameter(2); 
  double norme = Exp(-xmin/econst) - Exp(-xmax/econst);
  printf(" range (%f,%f) exp norm is %f and polya norm is %f \n",xmin,xmax,norme,pnorm );
  float spole=fp[0]->GetHistogram()->Integral("width");
  printf(" bin width %f pol integrates in the range (%f,%f) to  %e\n",fp[0]->GetHistogram()->GetBinWidth(0),xmin,xmax,spole);
 
  TString canTitle;
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  for(UInt_t ih=1; ih<hlist.size(); ++ih) { 
    double hintegral = hlist[ih]->Integral();
    fp[ih]->SetParameter(4,hintegral);
    hlist[ih]->SetNormFactor(hintegral);
    printf(" \n *********** fitting %s to pmt %i ************ \n",htag.Data(),ih);
    hlist[ih]->Fit(fp[ih],"R+","",4.4e-11,4.5e-10);
    canTitle.Form("%s-pmt%i",htag.Data(),ih);
    TCanvas *c0 = new TCanvas(canTitle,canTitle);
    gPad->SetLogy();
    hlist[ih]->SetLineColor(kBlack);
    hlist[ih]->Draw();
    c0->Update();
    c0->Print(".pdf");
  }

  cout << " fits " << htag << endl;
  for(int ih=0; ih<NPMT; ++ih ) {
    double width = hlist[ih]->GetBinWidth(0); // all same
    double qmax = hlist[ih]->GetBinLowEdge(hlist[ih]->GetMaximumBin()) + 0.5*width;
    double qfit=fp[ih]->GetParameter(1); 
    double qfite=fp[ih]->GetParError(1); 
    printf(" %i max bin %.3f qfit %.3f +/- %.3f \n",ih,qmax,qfit,qfite);
  }

  cout << " // fits " << htag << endl;
  if(htag.Contains("NoBeam")) {
      for(int ih=0; ih<NPMT; ++ih ) printf(" fitNoBeam[%i]=  %.3f ;\n",ih,fp[ih]->GetParameter(1));
  } else {
      for(int ih=0; ih<NPMT; ++ih ) printf(" fitOff[%i]=  %.3f ;\n",ih,fp[ih]->GetParameter(1));
  }
      
}

