//------------------------------------------
// Calculate Bayesian limit with uniform 
// prior for counting experiment with 
// a background and its uncertainty, as well
// as an uncertainty in the acceptance.
//
// All uncertinties are taken to be (truncated) gaussians.
// If the acceptance uncertainty is too large (>15-20%?), 
// the results should not be trusted.
//
//
// The function returns the 95% CL upper limit
// on the number of events (N_UL) so that N_UL
// can be directly plugged into the equation 
// for the upper limit on the cross-section:
//
// sigma_UL = N_Ul / (Acc*eff*lumi)
//
// It also plots the pdf and shows you the location
// of the UL.  You may need to change Nmin and Nmax 
// if the pdf is not "well contained" within the plot.
//
// Use from inside root
// root> .L ccBayes.C+
// root> ccBayes(....)
// 
//                               Claudio 16 July 2014
//-------------------------------------------------
#include <iostream>
#include "TCanvas.h"
#include "TLine.h"
#include "TMath.h"
#include "TH1D.h"
#include "TRandom3.h"

void     ccBayes(Int_t Nobs,       // number of observed events
                 Double_t Nbg,     // expected number of BG events
                 Double_t errbg,   // uncertainty on the above
                 Double_t erreff,  // fractional uncertainty on Acc*eff*lumi
                 Double_t Nmin=0., // Min value over which PDF is integrated
         Double_t Nmax=10.,// Max value over which PDF is integrated
                 Int_t Ntries=10000, // Number of MC smearings
                 Double_t step=0.01) // quantization of the pdf
{


  // Book the histogram that will hold the unnormalized PDF
  Int_t nbins = (Int_t) ((Nmax - Nmin)/step);
  TH1D* pdf = new TH1D("pdf","pdf",nbins,Nmin,Nmax);
  Double_t stepsize = (Nmax - Nmin)/nbins; 

  // Initialize the random number generator
  TRandom3 *random = new TRandom3();
  
  // Loop over MC smearings
  for (Int_t i=0; i<Ntries; i++) {

    // Get the BG, truncate gaussian at zero
    Double_t thisBG = random->Gaus(Nbg,errbg);
    if (thisBG < 0) continue;

    // Get the efficiency, truncate gaussian at zero
    Double_t thisEff = random->Gaus(1.,erreff);
    if (thisEff < 0) continue;

    // Fill the pdf
    Double_t N = Nmin - 0.5*stepsize;
    for (Int_t j=0; j<nbins; j++) {
      N = N + stepsize;              // center of j-th bin of the PDF
      Double_t thisN  = N/thisEff;   // efficiency correction
      Double_t mu = thisN + thisBG;
      Double_t pdfvalue = TMath::Exp(-mu) * TMath::Power(mu,(Double_t) Nobs);
      pdf->Fill(N,pdfvalue);
    }
  }

  // Now scan the pdf and find the 95% CL
  Double_t total = pdf->GetSum();
  Double_t integral=0.0;
  for (Int_t i=nbins; i>0; i--) {
    integral = integral + pdf->GetBinContent(i);
    if (integral > 0.05*total) {
      Double_t limit = pdf->GetBinLowEdge(i)+pdf->GetBinWidth(i);
      std::cout << "95% CL limit is " << limit << std::endl;
      TCanvas* c = new TCanvas();
      c->Draw();
      pdf->Draw();
      TLine* l = new TLine(limit,0.,limit,pdf->GetMaximum());
      l->Draw();
      break;
    }
  }
  
  delete pdf;

}
