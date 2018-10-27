//example of macro illustrating how to superimpose two histograms
//script to overlap multiple histograms from multiple partitions.
//to run the script do for example in an interactive root session: .x merge_over_partitions.cc++("h_ADC_offset", true, "barrel") if you want to have the h_ADC_offset histograms merged and normalised from all barrel parts of the detector, if you want only barrel type "barrel" for the partition parameter and if you want only disk type "disk". 
//where h_V125 is the histo you want to merge from all partitions from the ../Results/xxxx/fit_plots_xxxx.root files 
/*#include "TStyle.h"
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
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <TTreeReaderArray.h>
using namespace std;
*/


#include "TStyle.h"
#include "TGaxis.h"

#include <TTree.h>
#include "TFile.h"
//#include "TTreeReader.h"
//#include "TTreeReaderValue.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include "TLegendEntry.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include <TObjArray.h>
#include <TString.h>
#include <vector>
#include "TF2.h"
#include <fstream>
#include <TMultiGraph.h>
#include <iostream>
#include <string>
#include <TPaveLabel.h>

//#include "TBranch.h"
//#include "TTreeReaderArray.h"

using namespace std;
void mergeFlatTrees()
{
   //to run: .x mergePlotsAnalyzer.cc++()

  TCanvas *myCanvas = new TCanvas();

 
  ifstream runList("listFlatTrees.txt");
  TFile* output_f = TFile::Open("combined_FlatTrees.root","RECREATE");
  output_f->cd(); 

  ofstream nevents_file;
  nevents_file.open("nevents_file.txt");

  int totalNumberOfEvents = 0;
  int grand_totalNumberOfEvents = 0;


  TPad *pad_plot=new TPad("pad_plot","a transparent pad",0,0,.5,1);
  pad_plot->SetFillStyle(4000);
  pad_plot->Draw();


  TPad *pad_text=new TPad("pad_text","a transparent pad",0.5,0,1,1);
  pad_text->SetFillStyle(4000);
  pad_text->Draw();

  std::size_t pos = 0;       	  
  std::string for_legend;

  string line;
  int i_file = 0;
  while(std::getline(runList, line)){
 // for(unsigned int i_file = 0; i_file < v_files.size(); i_file++){
 	std::string sline = line.c_str();
	//cout << line.c_str()  << endl;
	//if(sline.find("NEW_RUN_AND_OR_YEAR") != string::npos){
	//if(strcmp(line,"NEW_RUN_AND_OR_YEAR") == 0){
	if(sline.length()  < 30){

  		output_f->cd(); 
		myCanvas->Write(for_legend.c_str());
		i_file = 0;
  		cout << "Total number of events, "  << totalNumberOfEvents << endl;
		nevents_file << "Total number of events, "  << totalNumberOfEvents << '\n';
		totalNumberOfEvents = 0;
		continue;
	}

        pos = line.find("combined") + 18;       	  
	for_legend = line.substr(pos);

	pad_plot->cd();	
 	TFile *f = new TFile(line.c_str());
 	f->cd("tree");
  	TTree *MyTree = (TTree*)f->Get("tree/SexaQAnalysis");
  	if(i_file  == 0)MyTree->Draw("nPV");
  	else MyTree->Draw("nPV", "", "same");
	cout << "NEvents " << for_legend.c_str() << ", " << MyTree->GetEntries() << endl;
	nevents_file << "NEvents " << for_legend.c_str() << ", " << MyTree->GetEntries() << '\n';
	((TH1F*)(gPad->GetListOfPrimitives()->At(i_file)))->SetLineColor(2+i_file);
	totalNumberOfEvents  = totalNumberOfEvents + MyTree->GetEntries();
	grand_totalNumberOfEvents = grand_totalNumberOfEvents + MyTree->GetEntries();	
	pad_text->cd();
	
 
	TPaveLabel *leg1 = new TPaveLabel(0,0.1+0.1*i_file,0.9,0.2+0.1*i_file,for_legend.c_str());
	leg1->SetFillColor(2+i_file);
	leg1->Draw(); 

        i_file++;

  } 

  cout << "Grand Total number of events, "  << grand_totalNumberOfEvents << endl;
  nevents_file << "Grand Total number of events, "  << grand_totalNumberOfEvents << '\n';

  nevents_file.close();
  
 
 
  output_f->Close();
   //OutputFile->cd();
   //c1->Write();

}

