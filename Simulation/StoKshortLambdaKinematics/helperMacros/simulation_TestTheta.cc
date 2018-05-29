/// \author Jarne De Clercq
// simple simulation for the kinematics of an S annihilating on a neutron and decaying to a Ks and a Lambda

//notations: all variables are in the detector rest frame, the once that are in the S+n rest frame are noted with "star", energy unit = GeV
//to run this for example: root -l -> in interactive mode:

//.x simulation.cc++("Xi_TSalis_polarised", "Xi1820polarised.root", "Xi1820", 100000, "polarisedJ32M12" )
//.x simulation.cc++("Xi_TSalis_unpolarised", "Xi1820unpolarised.root", 100000, "unpolarised" )

//.x simulation.cc++("S_ptTSalis", "S.root", "S", 100000, "unpolarised" )

#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"
#include "TEfficiency.h"


int simulation_TestTheta( string rootFileName, int nIterations)
{
	
    TRandom *random = new TRandom3();

    TH1F *h_theta_1 =  new TH1F(("h_theta1"), ("h_theta1"), 200, -5, 5);
    TH1F *h_Cos_theta_1 =  new TH1F(("h_Cos_theta_1"), ("h_Cos_theta_1"), 200, -5, 5);
    TH1F *h_delta_theta_1_2 =  new TH1F(("h_delta_theta_1_2"), ("h_delta_theta_1_2"),200,-5,5);

    int i = 0; 
 
    while(i<nIterations){

            //******************settings*****************************************//

          Double_t u = random->Uniform(-1.,1.);
          Double_t theta_1 = TMath::ACos(u);   
          Double_t delta_theta_1_2 = TMath::Pi() - 2*theta_1; 
          
          h_theta_1->Fill(theta_1);
          h_Cos_theta_1->Fill(TMath::Cos(theta_1));
          h_delta_theta_1_2->Fill(delta_theta_1_2);


          i++;
    }

 TFile *top = new TFile(rootFileName.c_str(),"RECREATE");
 TDirectory *dir_settings = top->mkdir("settings");
 dir_settings->cd();    
 
 h_theta_1->Write();
 h_Cos_theta_1->Write();
 h_delta_theta_1_2->Write();

 return 0;
    
}

