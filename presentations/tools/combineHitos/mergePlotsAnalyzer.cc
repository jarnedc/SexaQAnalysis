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
   //to run: .x mergePlotsAnalyzer.cc++()

   TFile *OutputFile = new TFile("Combine_Analyzer_plots.root","RECREATE");
   
 
   TFile *hfile1 = new TFile("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/Combined_by_run/2016_runB.root");
   TFile *hfile2 = new TFile("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/Combined_by_run/2016_runG.root");

   TFile *hfile3 = new TFile("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/Combined_by_run/2017_runB.root");
   TFile *hfile4 = new TFile("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/Combined_by_run/2017_runC.root");
   TFile *hfile5 = new TFile("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/Combined_by_run/2017_runD.root");
   TFile *hfile6 = new TFile("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/Combined_by_run/2017_runE.root");

   TFile *hfile7 = new TFile("/user/jdeclerc/Analysis/SexaQuark/CMSSW_9_4_7/src/SexaQAnalysis/V0_angular_correlation_analysis/bashRunAnalyzer/Results/Combined_by_run/2017_runF.root");

   vector<string> v_histos;
   v_histos.push_back("h_r_candidates_lxy_signed_delta_phi_between_1_and_2p5_error_lxy_smaller_0p1_after_LambdaKshortVertexFilter");
   //v_histos.push_back("h_Lambda_lxy_error");
   //v_histos.push_back("h_Lambda_lxy_signed");
   //v_histos.push_back("_sCand_L0_delta_phi");
   //v_histos.push_back("_sCand_antiL0_delta_phi");
   
   for(int i = 0; i < (int)v_histos.size(); i++  ){
	   TCanvas *c1 = new TCanvas(v_histos[i].c_str(),v_histos[i].c_str(),600,400);
	   gStyle->SetOptStat(kFALSE);


	   TH1F *h1 = (TH1F*)hfile1->Get(("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/dist_beamspot_S/"+v_histos[i]).c_str()); 
	   TH1F *h2 = (TH1F*)hfile2->Get(("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/dist_beamspot_S/"+v_histos[i]).c_str()); 
	   TH1F *h3 = (TH1F*)hfile3->Get(("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/dist_beamspot_S/"+v_histos[i]).c_str()); 
	   TH1F *h4 = (TH1F*)hfile4->Get(("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/dist_beamspot_S/"+v_histos[i]).c_str()); 
	   TH1F *h5 = (TH1F*)hfile5->Get(("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/dist_beamspot_S/"+v_histos[i]).c_str()); 
	   TH1F *h6 = (TH1F*)hfile6->Get(("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/dist_beamspot_S/"+v_histos[i]).c_str()); 
	   TH1F *h7 = (TH1F*)hfile7->Get(("Analyzer_V0_angular_correlation/LambdaKshortVertexFilter/dist_beamspot_S/"+v_histos[i]).c_str()); 
	   //TH1F *h8 = (TH1F*)hfile8->Get(("Analyzer_V0_angular_correlation/OriginalV0s/displacement/Lambda/"+v_histos[i]).c_str()); 

	   int binmax1 = h1->GetMaximumBin();
	   int binmax2 = h2->GetMaximumBin();
	   int binmax3 = h3->GetMaximumBin();
	   int binmax4 = h4->GetMaximumBin();
	   double binmax1content = h1-> GetBin(binmax1);
	   double binmax2content = h2-> GetBin(binmax2);
	   double binmax3content = h3-> GetBin(binmax3);
	   double binmax4content = h4-> GetBin(binmax4);

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
	   h4->SetLineColor(kBlack);
	   h4->DrawNormalized("sameE1l");
	   h5->SetLineColor(kYellow);
	   h5->DrawNormalized("sameE1l");
	   h6->SetLineColor(kCyan);
	   h6->DrawNormalized("sameE1l");
	   h7->SetLineColor(kMagenta);
	   h7->DrawNormalized("sameE1l");
	  // h8->SetLineColor(kOrange);
	  // h8->DrawNormalized("sameE1l");
	   
	   h1->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   h2->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   h3->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   h4->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   h5->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   h6->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	   h7->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);
	  // h8->GetYaxis()->SetRangeUser(0,binmaxcontent/maxNEntries*20);

	   c1->Update();
	   TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
   	  /* legend->AddEntry(h1,data[0].c_str(),"l");
   	   legend->AddEntry(h2,data[1].c_str(),"l");
   	   legend->AddEntry(h3,data[2].c_str(),"l");
   	   legend->AddEntry(h4,data[3].c_str(),"l");
   	  */ 
	   legend->AddEntry(h1,"2016 RunB","l");
	   legend->AddEntry(h2,"2016 RunG","l");
	   legend->AddEntry(h3,"2017 RunB","l");
	   legend->AddEntry(h4,"2017 RunC","l");
	   legend->AddEntry(h5,"2017 RunD","l");
   	   legend->AddEntry(h6,"2017 RunE","l");
   	   legend->AddEntry(h7,"2017 RunF","l");
   	   //legend->AddEntry(h8,"2017 RunF","l");
   	   legend->Draw();
   	   legend->Draw();
           OutputFile->cd();
           c1->Write();

   }

   //OutputFile->cd();
   //c1->Write();

}

