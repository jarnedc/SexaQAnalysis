//example of macro illustrating how to superimpose two histograms
//script to overlap multiple histograms from multiple partitions.
//to run the script do for example in an interactive root session: .x merge_over_partitions.cc++("h_ADC_offset", true, "barrel") if you want to have the h_ADC_offset histograms merged and normalised from all barrel parts of the detector, if you want only barrel type "barrel" for the partition parameter and if you want only disk type "disk". 
//where h_V125 is the histo you want to merge from all partitions from the ../Results/xxxx/fit_plots_xxxx.root files 
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
using namespace std;

void mergePlotsAnalyzer()
{
   //example of macro illustrating how to superimpose two histograms
   //with different scales in the "same" pad.
   // To see the output of this macro, click begin_html <a href="gif/twoscales.gif" >here</a> end_html
   //Author: Rene Brun

   TFile *OutputFile = new TFile("Combine_Analyzer_plots.root","RECREATE");
   TFile *hfile1 = new TFile("../MET.root");
   TFile *hfile2 = new TFile("../ZeroBias.root"); 
	
   vector<string> v_histo;
   v_histo.push_back("analyzer/WjetsMC_rCandMass_L0");
   v_histo.push_back("analyzer/WjetsMC_rCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors");
   v_histo.push_back("analyzer/WjetsMC_sCandMass_L0");
   v_histo.push_back("analyzer/WjetsMC_rCandMass_antiL0");
   v_histo.push_back("analyzer/WjetsMC_rCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors");
   v_histo.push_back("analyzer/WjetsMC_rCandMass_with_dxy_smaller_0p2_and_dr_daughters_smaller_0p8");
   v_histo.push_back("analyzer/WjetsMC_sCand_delta_phi");
   v_histo.push_back("analyzer/WjetsMC_sCand_delta_phi_dz(PCA_PV0)_below_1mm");
   v_histo.push_back("analyzer/WjetsMC_sCand_delta_eta_dz(PCA_PV0)_above_1mm");
   v_histo.push_back("analyzer/WjetsMC_sCand_delta_eta");
   v_histo.push_back("analyzer/WjetsMC_sCand_delta_R");
   v_histo.push_back("analyzer/WjetsMC_sCand_delta_R_dz(PCA_PV0)_below_1mm");
   v_histo.push_back("analyzer/WjetsMC_sCand_dxy(sCandVtx_beamspot)_signed_anti_L0");
   
   
   for(int i = 0; i < v_histo.size(); i++  ){
	   TCanvas *c1 = new TCanvas(v_histo[i].c_str(),v_histo[i].c_str(),600,400);
	   gStyle->SetOptStat(kFALSE);
	 //  auto legend = new TLegend(0.1,0.7,0.48,0.9);
	   
	   TH1F *h1 = hfile1->Get(v_histo[i].c_str()); 
	   TH1F *h2 = hfile2->Get(v_histo[i].c_str());

           


	   h1->Scale(1/h1->GetEntries());
	   h2->Scale(1/h2->GetEntries());

	   int binmax1 = h1->GetMaximumBin();
	   int binmax2 = h2->GetMaximumBin();

	   double binmax1content = h1-> GetBin(binmax1);
	   double binmax2content = h2-> GetBin(binmax2);

	   double binmaxcontent = binmax1content;
 	   if(binmax1content < binmax2content) binmaxcontent = binmax2content;
           
	   h1->GetYaxis()->SetRangeUser(0,binmaxcontent);
	   h2->GetYaxis()->SetRangeUser(0,binmaxcontent);

	   h1->Draw();
	   h1->SetLineColor(kRed);
	   h2->Draw("same");
	   h2->SetLineColor(kBlue);
	   h1->GetYaxis()->SetRangeUser(0,0.01);
	   h2->GetYaxis()->SetRangeUser(0,0.01);
	   c1->Update();
	   
   OutputFile->cd();
   c1->Write();

   }

   //OutputFile->cd();
   //c1->Write();

}

