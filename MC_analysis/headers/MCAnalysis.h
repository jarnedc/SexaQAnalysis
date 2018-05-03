#ifndef MCAnalysis_h
#define MCAnalysis_h

#include <map>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <iostream>
#include <sstream>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <math.h>
#include <string.h>
#include <boost/progress.hpp>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"


///GLOBAL VARIABLES

//particle pdgids
const Int_t Kshort_pdgid=310; //Kshort0
const Int_t Kshort_pdgid2=311; //Ks
const Int_t Kshort_pdgid3=130; //Klong0
const Int_t Lambda_pdgid=3122;
const Int_t Proton_pdgid=2212;
const Int_t Neutron_pdgid=2112;
const Int_t Piplus_pdgid=211;
const Int_t Xi1530_pdgid=3324; //Xi*0
const Int_t Xi1530_pdgid2=3314; //Xi*-

//Vector of pdgids for dividing in the three categories
// https://docs.google.com/spreadsheets/d/1dWtpiBBY6rUssBeWuRazb-ggjGvidOhGDALXyHHThLY/edit#gid=0

vector<Int_t> baryonMesonDecays_pdgids={3112, 3222, 3312, 3322, 421, 511, 521, 531, 5122, 5232 }; //Ks or L0 is product of mother being a meson/baryon and decaying
vector<Int_t> materialInteraction_pdgids={130, 211, 321, 2112, 2212, 3122, 310}; //Ks or L0 is product of mother parttaking in a material interaction
vector<Int_t> promptProduction_pdgids={311,333,10331,9010221,1,2,3,4,5,21,443,445,1103,2101,2103,2203,3101,3103,3114,3201,3203,3212,3214,3224,3303,4103,4122,4132,4232,20443,100441}; //Ks or L0 is product of mother being promptly produced during the collision

//c1 = prompt production
//c2 = material interaction
//c3 = decays

//constants
const Float_t pi=3.1415927;
const Float_t Mproton=0.938272; //GeV

//testvariables
Int_t testvar0=0;   Int_t testvar1=0;	Int_t testvar2=0;	Int_t testvar3=0;	Int_t testvar4=0;	Int_t testvar5=0;	
Int_t testvar6=0;	Int_t testvar7=0;	Int_t testvar8=0;	Int_t testvar9=0;	Int_t testvar10=0;


//other variables
Float_t reco_eta;
Float_t reco_rapidity;
Float_t reco_phi;
Float_t reco_pt;

//average amounts of particles
Float_t nr_gen_L0=0;
Float_t nr_gen_Ks=0;
Float_t nr_gen_L0_decay=0;
Float_t nr_gen_Ks_decay=0;
Float_t nr_gen_L0_decay_category=0;
Float_t nr_gen_Ks_decay_category=0;
Float_t nr_gen_L0_category=0;
Float_t nr_gen_Ks_category=0;
Float_t nr_reco_L0=0;
Float_t nr_reco_Ks=0;

//A bunch of Lorentz vectors

TLorentzVector lorentz;
TLorentzVector lorentz2;

TLorentzVector pion;

//deltaR and eta cut related variables
Float_t dR=0; //delta R
Float_t dRmax=0.02; //delta R cut value for calculating reco efficiencies
Float_t max_eta=2.5; //Cut on eta for the Lambda0s, Kshorts and their daughters, before their reco efficiencies are calculated

Float_t prompt_max_dr=10.0; //1 cm cut on prompt category

//Stuff related to vertex positions
TVector3 DV; //decay vertex position of some particle
TVector3 PV; //primary vertex position

TVector3 CV; //creation vertex position of some particle
TVector3 momCV; //creation vertex position of the mother of some particle
TVector3 momDV; //decay vertex position of the mother of some particle
TVector3 d1CV; //creation vertex position of the first daughter of some particle
TVector3 d2CV; //creation vertex position of the second daughter of some particle

TVector3 delta_DV_PV;
TVector3 delta_DV_momCV; 
TVector3 delta_PV_momCV;
TVector3 delta_CV_momCV;
TVector3 delta_DV_momDV; 
TVector3 delta_CV_PV;	

TH1F *temp; //temp TH1F histo for extracting 'passed' and 'total' histograms from TEfficiency histograms later on




///FUNCTION DEFINITIONS


void pdgid_amount_histo(TH1F *histogram, Int_t option){

	cout<<endl<<histogram->GetName()<<" content:"<<endl;
	float *array=histogram->fArray;
	Int_t Nelem=histogram->fN;
	cout <<"sizeofarray: "<<Nelem<<endl;

	if(option==0){
		cout<<"pdgid: "<<" amount: "<<endl;
		for(Int_t i=0; i < Nelem; i++){
			float value = array[i];
	    	if(value!=0) cout <<i-1<<" "<< value << endl;
		}
	} else if(option==1){
		cout<<"pdgid: "<<endl;
		for(Int_t i=0; i < Nelem; i++){
			float value = array[i];
	    	if(value!=0) cout <<i-1<< endl;
		}
	} else if(option==2){
		cout<<"amount: "<<endl;
		for(Int_t i=0; i < Nelem; i++){
			float value = array[i];
	    	if(value!=0) cout << value << endl;
		}
	}
	cout<<endl;

}


bool in_pdgid_list(Int_t pdgid, vector<Int_t> vec){ //absolute value is taken

	if (std::find(vec.begin(), vec.end(), abs(pdgid)) != vec.end()){ // found value in vec

		return true;

	} else return false;

}

Int_t get_category(Int_t mom_pdgid, Float_t dr){ //returns the category a given gen particle belongs to: 1, 2 or 3. 

	//think about grandmom pdgid potentially

	mom_pdgid = abs(mom_pdgid);

	if(in_pdgid_list(mom_pdgid, promptProduction_pdgids)){ 

		if(dr < prompt_max_dr){
			return 1;
		} else return 2;

	} else if(in_pdgid_list(mom_pdgid, materialInteraction_pdgids)){ 
		return 2;
	} else if(in_pdgid_list(mom_pdgid, baryonMesonDecays_pdgids)){
		return 3;
	} else {
		return 1; //put everything else in category 1 (promptly produced at collision)
	}
	return 0;
}

bool check_category(Int_t mom_pdgid, Float_t dr, Int_t Int_category){

	if(Int_category==123){
		return true;
	} else {

		if(get_category(mom_pdgid, dr)==Int_category) return true;
		return false;
	}

}


Float_t pythagoras2D(Float_t x, Float_t y){
	return sqrt(pow(x,2)+pow(y,2));	
}
Float_t pythagoras3D(Float_t x, Float_t y, Float_t z){
	return sqrt(pow(x,2)+pow(y,2)+pow(z,2));	
}

/*Float_t deltaR(Float_t e1, Float_t e2, Float_t p1, Float_t p2){ //e = eta; p = phi
	auto dp=std::abs(p1-p2); if (dp>3.1415927) dp-=2.0*3.1415927;  
    return std::sqrt((e1-e2)*(e1-e2) + dp*dp);
	
}
Float_t xy_norm(TVector3 *v){ //call as xy_norm(&vectorname)
	return sqrt(pow(v->X(),2)+pow(v->Y(),2));
}*/

Int_t cout_pdgid_name(Int_t pdgid){ //if pdgid is known: print out name of associated particle and return 1
									//if pdgid is not in switch choices: return 0
	
	switch (abs(pdgid)) {
        case Kshort_pdgid: cout << "Kshort 310"<<endl; return 1;
			break; 
		case Kshort_pdgid2: cout << "Kshort 311"<<endl; return 1;
			break; 
		case Lambda_pdgid: cout << "Lambda"<<endl; return 1;
			break; 
		case Proton_pdgid: cout << "Proton"<<endl; return 1;
			break; 
		case Neutron_pdgid: cout << "Neutron"<<endl; return 1;
			break; 
		case Piplus_pdgid: cout << "Pi+-"<<endl; return 1;
			break; 

		default: cout << pdgid <<endl; return 0;
    }
	
}

//saving all histos in the HQ_plots_... file and potentially saving them as pdf
void writeHistos(TFile*& outf, std::map<string, TH1F*> histos, std::map<string, TH2F*> histos2d, std::map<string, TEfficiency*> histoseff){
       
    for (std::pair<std::string, TH1F*> element : histos) {
      std::string word = element.first;
      TH1F* histo = element.second;
      outf->cd();
      //histo->SetName(histo->GetTitle());
      histo->Write();        
    }

    for (std::pair<std::string, TH2F*> element2 : histos2d) {
      std::string word = element2.first;
      TH2F* histo = element2.second;
      outf->cd();
      //histo->SetName(histo->GetTitle());
      histo->Write();     
    }
        
    for (std::pair<std::string, TEfficiency*> element3 : histoseff) {
      std::string word = element3.first;
      TEfficiency* histo = element3.second;
      outf->cd();
      //histo->SetName(histo->GetTitle());
      histo->Write();     
    }

}

#endif
