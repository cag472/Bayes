//------------------------------------------
// Calculate Bayesian limit with uniform 
// prior for counting experiment with 
// a background and its uncertainty, as well
// as an uncertainty in the acceptance.
//
// All uncertinties are taken to be (truncated) gaussians
// or lognormals.
// If the acceptance uncertainty is too large (>15-20%?), 
// the results using gaussian should not be trusted.
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
//
// Added a switch (and warnng!) to move on to Gaussian
// statistics when appropriate   Claudio 6 August 2014
//-----------------------------------------------------
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
                 Double_t Nmax=10., // Max value over which PDF is integrated
         Bool_t gaussian=true,   // Gaussian nuisances 
                 Int_t Ntries=10000, // Number of MC smearings
                 Double_t step=0.01,  // quantization of the pdf
         Bool_t usePois=true) // are you sure you want Poisson?
{


  // Warnings
  if (gaussian && erreff>0.15) {
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Your efficiency uncertainty is kind of high" << std::endl;
    std::cout << "Maybe you should consider lognormal nuisances" << std::endl; 
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
 }
  if (usePois && (Nobs > 20 || Nbg > 20)) {
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Your counts/backgrounds are high enough that" << std::endl;
    std::cout << "you should consider switching to gaussian" << std::endl;
    std::cout << "which can be achieved by usePois=false." << std::endl;
    std::cout << "If the program crashes or returns nonsense or" << std::endl;
    std::cout << "returns nothing, that is probably the reason." << std::endl;
    std::cout << " " << std::endl;
    std::cout << " " << std::endl;
  }

  // Initialize the random number generator
  TRandom3 *random = new TRandom3();

  // Book the histogram that will hold the unnormalized PDF
  Int_t nbins = (Int_t) ((Nmax - Nmin)/step);
  TH1D* pdf = new TH1D("pdf","pdf",nbins,Nmin,Nmax);
  Double_t stepsize = (Nmax - Nmin)/nbins; 

  // Book and fill histograms for the lognormal nuisance PDFs 
  // (technically truncated at 6 sigmas)
  TH1D* bPrior;
  TH1D* ePrior;
  if (errbg > 0) {
    bPrior = new TH1D("bPrior","bPrior log-normal",3000,
              TMath::Max(0.,Nbg-6*errbg), Nbg+6*errbg);
    for(Int_t i=0; i<1000000; i++) {
      bPrior->Fill(Nbg*(TMath::Power(1+errbg/Nbg,random->Gaus(0.,1.))));
    }
  }
  if (erreff>0) {
    ePrior = new TH1D("ePrior","ePrior log-normal",3000,
              TMath::Max(0.,1.-6*erreff), 1.+6*erreff);
    for(Int_t i=0; i<1000000; i++) {
      ePrior->Fill(TMath::Power(1+erreff, random->Gaus(0.,1.)));
    }
  }
  
  // Loop over MC smearings
  for (Int_t i=0; i<Ntries; i++) {

    // Get the BG, truncate gaussian at zero
    Double_t thisBG = Nbg;
    if (gaussian) {
    thisBG = random->Gaus(Nbg,errbg);
    if (thisBG < 0) continue;
    } else if (errbg > 0) {
      thisBG = bPrior->GetRandom();
    }

    // Get the efficiency, truncate gaussian at zero
    Double_t thisEff = 1.0;
    if (gaussian) {
      thisEff = random->Gaus(1.,erreff);
      if (thisEff < 0) continue;
    } else if (erreff > 0) {
      thisEff = ePrior->GetRandom();
    }

    // Fill the pdf
    Double_t N = Nmin - 0.5*stepsize;
    for (Int_t j=0; j<nbins; j++) {
      N = N + stepsize;              // center of j-th bin of the PDF
      Double_t thisN  = N/thisEff;   // efficiency correction
      Double_t mu = thisN + thisBG;
      Double_t pdfvalue;
      if (usePois) {
    pdfvalue = TMath::Exp(-mu) * TMath::Power(mu,(Double_t) Nobs);
      } else {
    // Dispense with niceties of factors of pi..they dont matter
    pdfvalue = (1./mu) * TMath::Exp(-(mu-Nobs)*(mu-Nobs)/(2*mu));
      }
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
  
  // Comment out the deletion in case one wants to replot it or save it
  // delete pdf;

}

