#include "./headers/exportplots.h"

void exportplots_MCAnalysis(TString category){ //root -l -q 'exportplots_MCAnalysis.C+("c123")'
	
	//string category = "c123"; //c1, c2, c3, c123, c1c2c3

	TString inputfile = TString("../rootfiles/HQ_plots_minbias_10M_events_skimmed_fixedGEANT")+TString("_")+category+TString(".root");
	TString outputfolder = TString("../plots/plots_MCAnalysis/")+category+TString("/");

	TFile *infile = new TFile(inputfile);

	int i=0;
	
	TIter nextkey(gDirectory->GetListOfKeys());

	while(TKey *key = (TKey*)nextkey()) {
		cout << "i: " << i << endl;
		++i;
		
		
		TString name = key->ReadObj()->GetName();
		TString title = key->ReadObj()->GetTitle();

		//cout << "name: " << "plots/"+title+".pdf" << endl;																																									
                                                                                                                                                                      
		if (key->ReadObj()->IsA()->InheritsFrom("TH1") && !(key->ReadObj()->IsA()->InheritsFrom("TH2")) && !(key->ReadObj()->IsA()->InheritsFrom("TEfficiency"))) {      
			
			gStyle->SetOptStat(111111);
			TH1* histo;  
			histo = (TH1*)key->ReadObj();   
			
			name = histo->GetName();
			histo->SetTitle(category+TString("_")+TString(histo->GetTitle()));
			TString filename = outputfolder + name + TString(".pdf");
			
			TCanvas *c1 = new TCanvas("c1","c1");
			c1->cd();
			c1->SetLogy(1);
			histo->SetMinimum(0.5);


			histo->Draw();		
			c1->SaveAs(filename,"pdf"); //c_str() converts sting to const char*                                                                                                                                                                             
			cout << histo->GetName() << endl;                                                                                                                                           
			                                                                                                                                                                                                                                                                                                                                        																																														  
		} else if (key->ReadObj()->IsA()->InheritsFrom("TH2")) { 
			gStyle->SetStatX(0.9);                
			// Set x-position (fraction of pad size)
			//gStyle->SetStatW(0.4);                
			// Set width of stat-box (fraction of pad size)
			//gStyle->SetStatH(0.2);                
			// Set height of stat-box (fraction of pad size)
			gStyle->SetOptStat(111111);     
			TH2* histo;                                                                                                                                                  
			histo = (TH2*)key->ReadObj();
			
			name = histo->GetName();
			histo->SetTitle(category+TString("_")+TString(histo->GetTitle()));

			TString filename = outputfolder + name + TString(".pdf");
			
			TCanvas *c1 = new TCanvas("c1","c1");
			c1->cd();
			c1->SetLogy(0);
			c1->SetLogz(1);
			histo->Draw("colz");		
			c1->SaveAs(filename,"pdf"); //.c_str() converts sting to const char*                                                                                                                                                                             
			cout << histo->GetName() << endl;                                                                                                                                                                              		                                                                                                                                                                                                                                                                                                                                         
																																															  
		} else if (key->ReadObj()->IsA()->InheritsFrom("TEfficiency")) {      
			TEfficiency* histo;                                                                                                                                                  
			histo = (TEfficiency*)key->ReadObj(); 
			
			name = histo->GetName();
			histo->SetTitle(category+TString("_")+TString(histo->GetTitle()));
			TString filename = outputfolder + name + TString(".pdf");
			
			TCanvas *c1 = new TCanvas("c1","c1");
			c1->cd();
			histo->Draw();		
			c1->SaveAs(filename,"pdf"); //c_str() converts sting to const char*                                                                                                                                                                             
			cout << histo->GetName() << endl;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
																																															  
		}
				

		

	}//while
}



