#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1I.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TLegend.h>
#include "TLegendEntry.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include <TObjArray.h>
#include <TString.h>
#include <vector>
#include "TF2.h"
#include "TH2.h"
#include <fstream>
#include <TGraph.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

template <typename T>

  std::string NumberToString ( T Number )
  {
     std::ostringstream ss;
     ss << Number;
     return ss.str();
  }


//void create2DPulseShape(const char * hName, bool normalize, string  partition)
void post_Analyzer_V0_angular_correlation()
{

   int delta_phi_min = -4;
   int delta_phi_max = 4;
   int delta_phi_steps = 100;
   int delta_eta_min = -10;
   int delta_eta_max = 10;
   int delta_eta_steps = 100;


   //TFile* f1 = new TFile("../Results/PulseShape_Hole_06-06-18_10:16/PulseShape.root");
   TFile* f1 = new TFile("/user/jdeclerc/Analysis/SexaQuark/CMSSW_8_0_30/src/SexaQAnalysis/analysed_SingleMuon_LambdaKshort_correlation_allFilters.root");
   f1->cd("Analyzer_V0_angular_correlation");

   TH2F* h2_delta_phi_delta_eta = (TH2F*)gDirectory->Get("ZeroBias_check_angular_corr_h_L0_Ks_delta_phi_delta_eta");

 


   TH2F* h2_delta_phi_delta_eta_normalised_on_delta_eta_0 = new TH2F("h2_delta_phi_delta_eta_normalised_on_delta_eta_0","h2_delta_phi_delta_eta_normalised_on_delta_eta_0;delta phi; delta eta", 100, -4,4, 100, -10,10);
  // TH2F* h2_delta_phi_delta_eta_normalised_on_delta_eta_0 = new TH2F("h2_delta_phi_delta_eta_normalised_on_delta_eta_0","h2_delta_phi_delta_eta_normalised_on_delta_eta_0;delta phi; delta eta", 100, h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_phi_min), h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_phi_max), 100, h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_eta_min), h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_eta_max));

   for(int i_delta_phi = h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_phi_min); i_delta_phi <= h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_phi_max); i_delta_phi++ ){

	int mean_delta_eta = 1+(h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_eta_max) - h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_eta_min))/2;
	//cout << "getting first delta eta 0 point " << i_delta_phi << " " << h2_delta_phi_delta_eta->GetBinContent(i_delta_phi,mean_delta_eta) << endl;
	Double_t n_at_delta_eta_0 = h2_delta_phi_delta_eta->GetBinContent(i_delta_phi,mean_delta_eta);

	for(int i_delta_eta = h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_eta_min); i_delta_eta <= h2_delta_phi_delta_eta->GetXaxis()->FindBin(delta_eta_max); i_delta_eta++  ){

		Double_t n_at_delta_eta_delta_phi = h2_delta_phi_delta_eta->GetBinContent(i_delta_phi, i_delta_eta);
	//	cout << "n_at_delta_eta_0 "  << n_at_delta_eta_0 << endl;
		if(n_at_delta_eta_0 > 0){
			//cout << "getting other point: " << i_delta_phi << " " << i_delta_eta << " " << n_at_delta_eta_delta_phi/n_at_delta_eta_0 << endl;
		//	h2_delta_phi_delta_eta_normalised_on_delta_eta_0->SetBinContent(i_delta_phi, i_delta_eta, n_at_delta_eta_delta_phi/n_at_delta_eta_0);
		//	cout << (double)delta_phi_min + (double)i_delta_phi*((double)delta_phi_max-(double)delta_phi_min)/(double)delta_phi_steps << " " << (double)delta_eta_min + (double)i_delta_eta*((double)delta_eta_max-(double)delta_eta_min)/(double)delta_eta_steps << endl;

			int bin = h2_delta_phi_delta_eta_normalised_on_delta_eta_0->FindBin((double)delta_phi_min + (double)i_delta_phi*((double)delta_phi_max-(double)delta_phi_min)/(double)delta_phi_steps, (double)delta_eta_min + (double)i_delta_eta*((double)delta_eta_max-(double)delta_eta_min)/(double)delta_eta_steps);
 
			h2_delta_phi_delta_eta_normalised_on_delta_eta_0->SetBinContent(bin, n_at_delta_eta_delta_phi/n_at_delta_eta_0);
		}
		else {cout << "eta =  0 entries is 0 at " << i_delta_phi << endl; }
	}
   }

   TFile *OutputFile = new TFile("post_Analyzer_V0_angular_correlation_SingleMuon_allFilters.root","RECREATE");
   //TCanvas *c1 = new TCanvas("c1","c1",600,400);
   gStyle->SetOptStat(kFALSE);

   //h2_delta_phi_delta_eta_normalised_on_delta_eta_0->Draw();
   OutputFile->cd();   
   h2_delta_phi_delta_eta_normalised_on_delta_eta_0->Write();


}

