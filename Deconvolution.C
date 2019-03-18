// Example to illustrate deconvolution function (class TSpectrum).
// To execute this example, do
// root > .x Deconvolution.C
#include "TSpectrum.h"
void Deconvolution() {
   Int_t i;
   Int_t nbins = 256;
   Double_t xmin  = 0;
   Double_t xmax  = nbins;
   Double_t * source = new Double_t[nbins];
   Double_t * response = new Double_t[nbins];
   TH1F *h = new TH1F("h","Deconvolution",nbins,xmin,xmax);
   TH1F *d = new TH1F("d","",nbins,xmin,xmax);
   TFile *f = new TFile("/usr/local/bin/root-6.10.06/tutorials/spectrum/TSpectrum.root");
   h=(TH1F*) f->Get("decon1;1");
   TFile *fr = new TFile("/usr/local/bin/root-6.10.06/tutorials/spectrum/TSpectrum.root");
   d=(TH1F*) fr->Get("decon_response;1");
   for (i = 0; i < nbins; i++) source[i]=h->GetBinContent(i + 1);
   for (i = 0; i < nbins; i++) response[i]=d->GetBinContent(i + 1);
   //TCanvas *Decon1 = gROOT->GetListOfCanvases()->FindObject("Decon1");
   //if (!Decon1) Decon1 = new TCanvas("Decon1","Decon1",10,10,1000,700);
   TCanvas *Decon1 =  new TCanvas("Decon1","Decon1",10,10,1000,700);
   h->Draw("L");
   TSpectrum *s = new TSpectrum();
   s->Deconvolution(source,response,256,1000,1,1);
   for (i = 0; i < nbins; i++) d->SetBinContent(i + 1,source[i]);
   d->SetLineColor(kRed);
   d->Draw("SAME L");
}
