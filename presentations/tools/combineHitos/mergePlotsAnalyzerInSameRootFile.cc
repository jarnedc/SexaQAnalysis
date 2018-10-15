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

void mergePlotsAnalyzerInSameRootFile()
{
   //to run: .x mergePlotsAnalyzer.cc++()

   TFile *OutputFile = new TFile("Combine_Analyzer_plots.root","RECREATE");
   
 
   TFile *hfile1 = new TFile("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/Combined_by_run151018/Combined_by_run/combined_2017.root");

   TH1F *h1 = (TH1F*)hfile1->Get("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/Angular_Correlation_Daughters_S/Angular_Correlation_Daughters_S_kinematics_S_in_correlation/h_s_candidates_lxy_signed_in_delta_phi_delta_eta_corr");
   TH1F *h2 = (TH1F*)hfile1->Get("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/Angular_Correlation_Daughters_S/Angular_Correlation_Daughters_S_kinematics_S_in_peak/h_s_candidates_lxy_signed_in_peak");
   TH1F *h3 = (TH1F*)hfile1->Get("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/Angular_Correlation_Daughters_S/Angular_Correlation_Daughters_S_kinematics_S_outside_correlation/h_s_candidates_lxy_signed_outside_delta_phi_delta_eta_corr");
   //TH1F *h4 = (TH1F*)hfile1->Get("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/PCA_R/h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_then_30_smaller_or_equal_to_40");
   //TH1F *h5 = (TH1F*)hfile1->Get("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/PCA_R/h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_then_40");
   //v_histos.push_back("_sCand_L0_delta_phi");
   //v_histos.push_back("_sCand_antiL0_delta_phi");
   
	   TCanvas *c1 = new TCanvas("c1","c1",600,400);
	   gStyle->SetOptStat(kFALSE);



	   int binmax1 = h1->GetMaximumBin();
	   int binmax2 = h2->GetMaximumBin();
	   int binmax3 = h3->GetMaximumBin();
	   //int binmax4 = h4->GetMaximumBin();
	   double binmax1content = h1-> GetBin(binmax1);
	   double binmax2content = h2-> GetBin(binmax2);
	   double binmax3content = h3-> GetBin(binmax3);
	   //double binmax4content = h4-> GetBin(binmax4);

	   double binmaxcontent = binmax1content;
 	   if(binmax1content < binmax2content) binmaxcontent = binmax2content; 
	   if(binmaxcontent < binmax3content) binmaxcontent = binmax3content;
	
           int maxNEntries = h1->GetEntries();
	   if(h1->GetEntries() > h2->GetEntries()) maxNEntries = h2->GetEntries();
	   if(h2->GetEntries() > h3->GetEntries()) maxNEntries = h3->GetEntries();
          
	   h1->SetLineColor(kRed);
	   h1->DrawNormalized("E1l");
	   h2->SetLineColor(kBlue);
	   h2->DrawNormalized("sameE1l");
	   h3->SetLineColor(kGreen);
	   h3->DrawNormalized("sameE1l");
	   //h4->SetLineColor(kBlack);
	   //h4->DrawNormalized("sameE1l");
	   //h5->SetLineColor(kOrange);
	   //h5->DrawNormalized("sameE1l");
	   
	   h1->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   h2->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   h3->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   //h4->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   //h5->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);

	   c1->Update();
	   TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
   	  /* legend->AddEntry(h1,data[0].c_str(),"l");
   	   legend->AddEntry(h2,data[1].c_str(),"l");
   	   legend->AddEntry(h3,data[2].c_str(),"l");
   	   legend->AddEntry(h4,data[3].c_str(),"l");
   	  */ 
	   legend->AddEntry(h1,"in correlation (fabs(delta_phi) < 0.1 && fabs(delta_eta) > 0.5 && fabs(delta_eta) < 3)","l");
   	   legend->AddEntry(h2,"in peak (fabs(delta_phi < 0.1 && fabs(delta_eta) < 0.1))","l");
   	   legend->AddEntry(h3,"others","l");
   	   //legend->AddEntry(h4,"30 < nPVs <= 40","l");
   	   //legend->AddEntry(h5,"nPVs > 40","l");
   	   legend->Draw();
   	   legend->Draw();
           OutputFile->cd();
           c1->Write();


   //OutputFile->cd();
   //c1->Write();

}

