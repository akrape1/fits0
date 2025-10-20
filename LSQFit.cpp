#include "TRandom2.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGClient.h"
#include "TStyle.h"


#include <iostream>
using namespace std;

using TMath::Log;

//parms
const double xmin=1;
const double xmax=20;
const int npoints=12;
const double sigma=0.2;

double f(double x){
  const double a=0.5;
  const double b=1.3;
  const double c=0.5;
  return a+b*Log(x)+c*Log(x)*Log(x);
}

void getX(double *x){
  double step=(xmax-xmin)/npoints;
  for (int i=0; i<npoints; i++){
    x[i]=xmin+i*step;
  }
}

void getY(const double *x, double *y, double *ey){
  static TRandom2 tr(0);
  for (int i=0; i<npoints; i++){
    y[i]=f(x[i])+tr.Gaus(0,sigma);
    ey[i]=sigma;
  }
}

/*
//this beast won't do anything when I run my script so get commented silly method
void leastsq(){
  double x[npoints];
  double y[npoints];
  double ey[npoints];
  getX(x);
  getY(x,y,ey);
  auto tg = new TGraphErrors(npoints,x,y,0,ey);
  tg->Draw("alp");
}
*/

//I run all of my stuff directly in ROOT so a lot of the skeleton script i threw out
void LSQFit() {
  gStyle->SetOptStat(0);
  const int nexperiments = 1000;

  double x[npoints], y[npoints], ey[npoints];
  getX(x);

  TH2F *hab = new TH2F("hab", "Parameter b vs a;a;b", 60, 0.3, 0.7, 60, 1.1, 1.5);
  TH2F *hac = new TH2F("hac", "Parameter c vs a;a;c", 60, 0.3, 0.7, 60, 0.3, 0.7);
  TH2F *hbc = new TH2F("hbc", "Parameter c vs b;b;c", 60, 1.1, 1.5, 60, 0.3, 0.7);
  TH1F *hchi_red = new TH1F("hchi", "Reduced #chi^{2};#chi^{2}_{red};Frequency", 80, 0, 3);

  TH1F *ha = new TH1F("ha", "Parameter a;a;Frequency", 80, 0.3, 0.7);
  TH1F *hb = new TH1F("hb", "Parameter b;b;Frequency", 80, 1.1, 1.5);
  TH1F *hc = new TH1F("hc", "Parameter c;c;Frequency", 80, 0.3, 0.7);
  TH1F *hchi = new TH1F("hchi", "#chi^{2};#chi^{2};Frequency", 80, 0, 30);

  for (int iexp = 0; iexp < nexperiments; iexp++) {
    getY(x, y, ey);

    TMatrixD A(npoints, 3);
    TVectorD Y(npoints);
    TMatrixD W(npoints, npoints);

    for (int i = 0; i < npoints; i++) {
      double lx = Log(x[i]);
      A(i, 0) = 1.0;
      A(i, 1) = lx;
      A(i, 2) = lx * lx;
      Y(i) = y[i];
      W(i, i) = 1.0 / (ey[i] * ey[i]);
    }

    TMatrixD AT(TMatrixD::kTransposed, A);
    TMatrixD ATW = AT * W;
    TMatrixD ATA = ATW * A;
    TMatrixD cov(ATA);
    cov.Invert();

    TVectorD rhs = ATW * Y;
    TVectorD best = cov * rhs;

    double a = best[0];
    double b = best[1];
    double c = best[2];

    double chi2 = 0.0;
    for (int i = 0; i < npoints; i++) {
      double yfit = a + b * Log(x[i]) + c * Log(x[i]) * Log(x[i]);
      chi2 += pow((y[i] - yfit) / ey[i], 2);
    }
    double chi2_reduced = chi2 / (npoints - 3);

    hab->Fill(a, b);
    hac->Fill(a, c);
    hbc->Fill(b, c);
    hchi_red->Fill(chi2_reduced);

    ha->Fill(a);
    hb->Fill(b);
    hc->Fill(c);
    hchi->Fill(chi2);
  }

  TCanvas *c = new TCanvas("c", "Chi2 Linear Fit Results", 1000, 800);
  c->Divide(2, 2);
  c->cd(1); hab->Draw("colz");
  c->cd(2); hac->Draw("colz");
  c->cd(3); hbc->Draw("colz");
  c->cd(4); hchi_red->Draw();
  c->SaveAs("redChi_root.png");

  TCanvas *c2 = new TCanvas("c2", "Parameter and Chi2 Distributions", 1000, 800);
  c2->Divide(2, 2);
  c2->cd(1); ha->Draw();
  c2->cd(2); hb->Draw();
  c2->cd(3); hc->Draw();
  c2->cd(4); hchi->Draw();
  c2->SaveAs("chi_root.png");
}
