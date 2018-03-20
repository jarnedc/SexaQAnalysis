#include "./headers/exportplots.h"

void exportplots_MCAnalysis_sames(){ 
	
	//string category = "c123"; //c1, c2, c3, c123, c1c2c3
	string inputfile_c1 = string("../rootfiles/HQ_plots_minbias_10M_events_skimmed_fixedGEANT_c1.root");
	string inputfile_c2 = string("../rootfiles/HQ_plots_minbias_10M_events_skimmed_fixedGEANT_c2.root");
	string inputfile_c3 = string("../rootfiles/HQ_plots_minbias_10M_events_skimmed_fixedGEANT_c3.root");
	string inputfile_c123 = string("../rootfiles/HQ_plots_minbias_10M_events_skimmed_fixedGEANT_c123.root");

	string outputfolder = string("../plots/plots_MCAnalysis/c1c2c3/"); //output

	TFile *infile_c1 = new TFile(inputfile_c1.c_str());
	TFile *infile_c2 = new TFile(inputfile_c2.c_str());
	TFile *infile_c3 = new TFile(inputfile_c3.c_str());
	TFile *infile_c123 = new TFile(inputfile_c123.c_str());


	TIter nextkey1(infile_c1->GetListOfKeys());
	TIter nextkey2(infile_c2->GetListOfKeys());
	TIter nextkey3(infile_c3->GetListOfKeys());
	TIter nextkey123(infile_c123->GetListOfKeys());
	int i=0;
	while(TKey *key1 = (TKey*)nextkey1()) {

		cout << "i: " << i << endl;
		++i;

		TKey *key2 = (TKey*)nextkey2();
		TKey *key3 = (TKey*)nextkey3();
		TKey *key123 = (TKey*)nextkey123();
		

		/*
		string name1 = key1->ReadObj()->GetName();
		string title1 = key1->ReadObj()->GetTitle();

		cout << title1 <<endl;

		string name2 = key2->ReadObj()->GetName();
		string title2 = key2->ReadObj()->GetTitle();

		cout << title2 <<endl;
		*/

		if (key1->ReadObj()->IsA()->InheritsFrom("TH1") && !(key1->ReadObj()->IsA()->InheritsFrom("TH2")) && !(key1->ReadObj()->IsA()->InheritsFrom("TEfficiency"))) {      
			
			gStyle->SetOptStat(111111);

			TH1* histo1;  
			histo1 = (TH1*)key1->ReadObj();  
			TH1* histo2;  
			histo2 = (TH1*)key2->ReadObj();
			TH1* histo3;  
			histo3 = (TH1*)key3->ReadObj();
			TH1* histo123;  
			histo123 = (TH1*)key123->ReadObj(); 
			
			string name = histo1->GetName();
			string filename = outputfolder + name + string(".pdf");
			
			TCanvas *c1 = new TCanvas("c1","c1");
			c1->cd();
			c1->SetLogy(1);
			histo1->SetMinimum(0.5);
			histo2->SetMinimum(0.5);
			histo3->SetMinimum(0.5);
			histo123->SetMinimum(0.5);

			Float_t w = 0.18; //width fraction
			Float_t h = 0.18; //height fraction
			Float_t x = 0.82;
			Float_t y = 0.75;

			histo123->SetLineColor(4);
			histo123->Draw();
			gPad->Update();
			TPaveStats *stats123 = (TPaveStats*)histo123->GetListOfFunctions()->FindObject("stats");
			stats123->SetName("c123");
			histo123->SetName("c123");
			stats123->SetY1NDC(y-h);
			stats123->SetY2NDC(y);
			stats123->SetX1NDC(x);
			stats123->SetX2NDC(x+w);
			stats123->SetTextColor(4);

			
			histo1->SetLineColor(1);
			histo1->Draw("sames");
			gPad->Update();
			TPaveStats *stats1 = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
			stats1->SetName("c1");
			histo1->SetName("c1");
			stats1->SetY1NDC(y);
			stats1->SetY2NDC(y+h);
			stats1->SetX1NDC(x-w);
			stats1->SetX2NDC(x);
			stats1->SetTextColor(1);

			
			histo2->SetLineColor(2);
			histo2->Draw("sames");
			gPad->Update();
			TPaveStats *stats2 = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
			stats2->SetName("c2");
			histo2->SetName("c2");
			stats2->SetY1NDC(y);
			stats2->SetY2NDC(y+h);
			stats2->SetX1NDC(x);
			stats2->SetX2NDC(x+w);
			stats2->SetTextColor(2);

			
			histo3->SetLineColor(3);
			histo3->Draw("sames");
			gPad->Update();
			TPaveStats *stats3 = (TPaveStats*)histo3->GetListOfFunctions()->FindObject("stats");
			stats3->SetName("c3");
			histo3->SetName("c3");
			stats3->SetY1NDC(y-h);
			stats3->SetY2NDC(y);
			stats3->SetX1NDC(x-w);
			stats3->SetX2NDC(x);
			stats3->SetTextColor(3);

			


			histo1->SetTitle(TString("c1c2c3_")+TString(histo1->GetTitle()));


			c1->SaveAs(filename.c_str(),"pdf"); //c_str() converts sting to const char*                                                                                                                                                                             
			cout << histo1->GetName() << endl;  
		}



	}
	
}
