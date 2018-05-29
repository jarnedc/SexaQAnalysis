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
   
   vector<string> data;
   data.push_back("MET");
   data.push_back("SingleMuon");
   
   TFile *hfile1 = new TFile(("../Data_analysis/Results/"+data[0]+".root").c_str());
   TFile *hfile2 = new TFile(("../Data_analysis/Results/"+data[1]+".root").c_str());
	
   vector<string> v_histos;
   v_histos.push_back("_sCand_delta_phi");
   v_histos.push_back("_sCand_L0_delta_phi");
   v_histos.push_back("_sCand_antiL0_delta_phi");
   
   for(int i = 0; i < (int)v_histos.size(); i++  ){
	   TCanvas *c1 = new TCanvas(v_histos[i].c_str(),v_histos[i].c_str(),600,400);
	   gStyle->SetOptStat(kFALSE);
	 //  auto legend = new TLegend(0.1,0.7,0.48,0.9);
	   cout << ("analyzer/"+data[0]+v_histos[i]).c_str() << endl;
	   cout << ("analyzer/"+data[1]+v_histos[i]).c_str() << endl; 
	   TH1F *h1 = (TH1F*)hfile1->Get(("analyzer/"+data[0]+v_histos[i]).c_str()); 
	   TH1F *h2 = (TH1F*)hfile2->Get(("analyzer/"+data[1]+v_histos[i]).c_str()); 

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

	   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   	   legend->AddEntry(h1,data[0].c_str(),"l");
   	   legend->AddEntry(h2,data[1].c_str(),"l");
   	   legend->Draw();
	   
           OutputFile->cd();
           c1->Write();

   }

   //OutputFile->cd();
   //c1->Write();

}

