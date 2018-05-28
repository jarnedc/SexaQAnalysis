#include "./headers/exportplots.h"

void exportplots_new_setup(string name){ //"SingleMuonData" of "WjetsMC" of "ZeroBiasData"
	
	
	string inputfile = string("../../rootfiles/analysis_")+name+string(".root");
	string outputfolder = string("../../plots/plots_new_setup/plots_")+name+string("/");

	TFile *infile = new TFile(inputfile.c_str());

	int i=0;
	
	TIter nextkey(gDirectory->GetListOfKeys());

	while(TKey *key2 = (TKey*)nextkey()) {



		if ( key2->ReadObj()->IsA()->InheritsFrom("TDirectory") ){
			cout <<"test"<<endl;

			std::string name = key2->ReadObj()->GetName();
			std::cout << "name: " << name << std::endl;
			TDirectory *dir = (TDirectory*)key2->ReadObj();
			dir->cd();
			TIter nextkey(gDirectory->GetListOfKeys());


			while(TKey *key = (TKey*)nextkey()) {
				cout << "i: " << i << endl;
				++i;
				
				
				string name = key->ReadObj()->GetName();
				string title = key->ReadObj()->GetTitle();

				//cout << "name: " << "plots/"+title+".pdf" << endl;																																									
		                                                                                                                                                                      
				if (key->ReadObj()->IsA()->InheritsFrom("TH1") && !(key->ReadObj()->IsA()->InheritsFrom("TH2")) && !(key->ReadObj()->IsA()->InheritsFrom("TEfficiency"))) {      
					gStyle->SetOptStat(111111);
					TH1* histo;  
					histo = (TH1*)key->ReadObj();   
					
					name = histo->GetName();
					string filename = outputfolder + name + string(".pdf");
					
					TCanvas *c1 = new TCanvas("c1","c1");
					c1->cd();
					c1->SetLogy(0);

					histo->Draw();		
					c1->SaveAs(filename.c_str(),"pdf"); //c_str() converts sting to const char*                                                                                                                                                                             
					cout << histo->GetName() << endl;                                                                                                                                           
					                                                                                                                                                                                                                                                                                                                                        																																														  
				} else if (key->ReadObj()->IsA()->InheritsFrom("TH2")) {
					gStyle->SetStatX(0.9);
					gStyle->SetOptStat(1111);      
					TH2* histo;                                                                                                                                                  
					histo = (TH2*)key->ReadObj();
					
					name = histo->GetName();
					string filename = outputfolder + name + string(".pdf");
					
					TCanvas *c1 = new TCanvas("c1","c1");
					c1->cd();
					c1->SetLogy(0);
					c1->SetLogz(1);
					histo->Draw("colz");		
					c1->SaveAs(filename.c_str(),"pdf"); //c_str() converts sting to const char*                                                                                                                                                                             
					cout << histo->GetName() << endl;                                                                                                                                                                              		                                                                                                                                                                                                                                                                                                                                         
																																																	  
				} else if (key->ReadObj()->IsA()->InheritsFrom("TEfficiency")) {      
					TEfficiency* histo;                                                                                                                                                  
					histo = (TEfficiency*)key->ReadObj(); 
					
					name = histo->GetName();
					string filename = outputfolder + name + string(".pdf");
					
					TCanvas *c1 = new TCanvas("c1","c1");
					c1->cd();
					histo->Draw();		
					c1->SaveAs(filename.c_str(),"pdf"); //c_str() converts sting to const char*                                                                                                                                                                             
					cout << histo->GetName() << endl;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
																																																	  
				}
						

				

			}//while
		}
	}
}









/* //OLD CODE

void exportplots(){
	
	//read TFile
	char inputfile[] = "HQ_plots_minbias_10M_events_skimmed_fixedGEANT.root";

	TFile *infile = new TFile(inputfile);

	//loop counter(s)
	int i=0;
	int j=0;
	
	TIter nextkey(gDirectory->GetListOfKeys());

	while(TKey *key = (TKey*)nextkey()) {
		std::cout << "i: " << i << std::endl;
		++i;
		//obj = key->ReadObj();
		cout<<typeid(key->ReadObj()).name()<<endl;
		
		std::string name = obj->GetName();
		std::cout << "name: " << name << std::endl;
		if ( obj->IsA()->InheritsFrom("TDirectory") ){
																																																	  TIter nextkey(gDirectory->GetListOfKeys());                                                                                                                                                     
					 while (TKey *key == (TKey*)nextkey()) {                                                                                                                                                                
					 obj = key->ReadObj();                                                                                                                                                                         
					 if (obj->IsA()->InheritsFrom("TH1")) {                                                                                                                                                        
					 h = (TH1*)obj;                                                                                                                                                                              
					 std::cout << h->GetName() << std::endl;                                                                                                                                                     
					 std::cout << "j: " << j << std::endl;                                                                                                                                                         
					 ++j;                                                                                                                                                                                          
				 }                                                                                                                                                                                               
			 }
		}

	}//while
}
*/
