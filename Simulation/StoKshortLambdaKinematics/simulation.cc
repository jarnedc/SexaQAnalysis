/// \author Jarne De Clercq
// simple simulation for the kinematics of an S annihilating on a neutron and decaying to a Ks and a Lambda

//notations: all variables are in the detector rest frame, the once that are in the S+n rest frame are noted with "star", energy unit = GeV

#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"


void simulation()
{
	
    TRandom *random = new TRandom3();


    //constants
    Double_t M_S = 1.8; //S particle mass
    Double_t m_Ks = 0.497611;//Kshort mass
    Double_t m_l = 1.115683;//Lambda mass
    Double_t m_n = 0.9395654133;//Lambda mass
    Double_t M_S_n = M_S + m_n;//Lambda mass

    //plots
    TH1F *h_delta_phi_Ks_l = new TH1F("h_delta_phi_Ks_l", "h_delta_phi_Ks_l", 100, 0, 7);
    TH1F *h_delta_theta_Ks_l = new TH1F("h_delta_theta_Ks_l", "h_delta_theta_Ks_l", 100, 0, 4);
    
    int i = 0; 
   
    while(i<1000){

	//******************calculate the S particle parameters******************************************//
	//given
	//p_S the size of the 3 momentum of the S particle, should become a distribution at some point, p_t goes exponential?
	Double_t p_S = 2.;
	//theta of the S, should also be a distribution at some point
	//from https://arxiv.org/pdf/1002.0621.pdf, page 7 the eta distribution is relatively flat
	Double_t eta = random->Uniform(-2,2);
	Double_t theta_S = 2*TMath::ATan(TMath::Exp(-eta));
	//Double_t theta_S = TMath::Pi()*1./3.; 
	//phi of the S, should be a uniform distribution in 0 to 2*pi
	//Double_t phi_S = TMath::Pi();
	Double_t phi_S = random->Uniform(2.*TMath::Pi());	
	//p3_S this is the 3 vector of the S particle can be calculated from p_t_S and theta_S
	TVector3 p3_S(0.,0.,0.); 
	p3_S.SetMagThetaPhi(p_S,theta_S,phi_S); 
	cout << "momenta S particle: " << p3_S.Px() << " " << p3_S.Py() << " " << p3_S.Pz() << endl;
	//******************calculate the S particle parameters******************************************//
			
	//the S hits on a neutron and transfers its momentum to the S+n system.
	
	//******************calculate the S+n particle parameters******************************************//
	//assume the S momentum is transfered to the S+n system, then the momenta are the same
	TVector3 p3_S_n =  p3_S;	
	Double_t E_S = pow(M_S_n*M_S_n+p_S*p_S,0.5);
  	cout << "Energy of the S and n system: " << E_S << endl;	
	cout << "-----------------------------"<< endl;
	//******************calculate the S+n particle parameters******************************************//
	
	//The S + n then decays to a Ks and Lambda. 
	//This is a 2 body decay (http://www.helsinki.fi/~www_sefo/phenomenology/Schlippe_relativistic_kinematics.pdf)
	
	//******************calculate the Ks and Lambda parameter****************************************//
	//E_ks_star and p_Ks_star energy and momentum of the Ks in the rest frame of the S + n
	Double_t E_ks_star = (M_S_n*M_S_n+m_Ks*m_Ks-m_l*m_l)/(2*M_S_n);
	Double_t p_Ks_star = pow(E_ks_star,2) - m_Ks*m_Ks;
	//E_l_star and p_l_star energy and momentum of the l in the rest frame of the S + n
	Double_t E_l_star = (M_S_n*M_S_n+m_l*m_l-m_Ks*m_Ks)/(2*M_S_n);
	Double_t p_l_star = pow(E_l_star,2) - m_l*m_l;
	cout << "momentum of Ks and Lambda in the rest frame of S+n, should be the same...: " << p_Ks_star << " " << p_l_star << endl;
	//theta of the Ks, should be uniform distribution from 0 to pi at some point
	//Double_t theta_Ks_star = TMath::Pi()*0.5;
	Double_t theta_Ks_star = random->Uniform(TMath::Pi());	
	//phi of the Ks, should be uniform distribution from 0 to 2*pi at some point
	//Double_t phi_Ks_star = TMath::Pi()*0.5;
	Double_t phi_Ks_star = random->Uniform(2.*TMath::Pi());	
	//theta of the Lambda, back to back with Ks
	Double_t theta_l_star = TMath::Pi()-theta_Ks_star;
	//phi of the Lambda, back to back with Ks
	Double_t phi_l_star = TMath::Pi()+phi_Ks_star;
	//p4_Ks_star 4 momentum vector of the Ks 
	TLorentzVector p4_Ks_star(0.,0.,0.,0.);
	cout << "theta_Ks_star " << theta_Ks_star << endl;
	p4_Ks_star.SetPtEtaPhiE(p_Ks_star*TMath::Sin(theta_Ks_star),-TMath::Log(TMath::Tan(theta_Ks_star/2.)),phi_Ks_star,E_ks_star);
	cout << "p4_Ks_star: " << p4_Ks_star.Pt() << " " << p4_Ks_star.Eta() << " " << p4_Ks_star.Phi() << " " << p4_Ks_star.E() << endl;
	//p4_Ks_star.SetRho(p_Ks_star); 
        //p4_Ks_star.SetTheta(theta_Ks_star);
        //p4_Ks_star.SetPhi(phi_Ks_star); 
        //p4_Ks_star.SetE(E_ks_star); 
	//p3_l_star 3 momentum vector of the l
	TLorentzVector p4_l_star(0.,0.,0.,0.);
	cout << "theta_l_star " << theta_l_star << endl;
	p4_l_star.SetPtEtaPhiE(p_l_star*TMath::Sin(theta_l_star),-TMath::Log(TMath::Tan(theta_l_star/2.)),phi_l_star,E_l_star);	
	cout << "p4_l_star: " << p4_l_star.Pt() << " " << p4_l_star.Eta() << " " << p4_l_star.Phi() << " " << p4_l_star.E() << endl;
	//p4_l_star.SetRho(p_l_star); 
        //p4_l_star.SetTheta(theta_l_star);
        //p4_l_star.SetPhi(phi_l_star);	
        //p4_l_star.SetE(E_l_star);	
	//******************calculate the Ks and Lambda parameter****************************************//

	//you have now the vector of the Ks and the Lambda in the rest frame of the S+n. Now will boost these vectors to the reference
	//frame where the S+n is actually moving with momentum p3_S, this is the detector reference frame
	
	
	//******************boost to the reference frame of the detector****************************************//
	//boosting the p4_Ks_star and p4_l_star allong the S+n system momentum	
	cout << "-----------------------------"<< endl;
	cout << "before boost: "<< endl;
	cout << "p4_Ks_star: " << p4_Ks_star.Px() << " " << p4_Ks_star.Py() << " " << p4_Ks_star.Pz() <<  " "<< p4_Ks_star.E() << endl;
	cout << "p4_l_star: " << p4_l_star.Px() << " " << p4_l_star.Py() << " " << p4_l_star.Pz() << " " << p4_Ks_star.E() << endl;
	Double_t delta_p = p4_Ks_star.Rho()-p4_l_star.Rho();
        Double_t delta_phi = p4_Ks_star.Phi()-p4_l_star.Phi();
        Double_t delta_theta = p4_Ks_star.Theta()-p4_l_star.Theta();
        cout << "delta_p " << delta_p << " delta_phi " << delta_phi << " delta_theta " << delta_theta << endl;	


	p4_Ks_star.Boost(p3_S_n.x()/E_S,p3_S_n.y()/E_S,p3_S_n.z()/E_S); //beta = p/E
	p4_l_star.Boost(p3_S_n.x()/E_S,p3_S_n.y()/E_S,p3_S_n.z()/E_S); //beta = p/E


	cout << "-----------------------------"<< endl;
	cout << "after boost: " << endl;
	cout << "p4_Ks_star: " << p4_Ks_star.Px() << " " << p4_Ks_star.Py() << " " << p4_Ks_star.Pz() <<  " "<< p4_Ks_star.E() << endl;
        cout << "p4_l_star: " << p4_l_star.Px() << " " << p4_l_star.Py() << " " << p4_l_star.Pz() << " " << p4_Ks_star.E() << endl;
	delta_p = p4_Ks_star.Rho()-p4_l_star.Rho();
	delta_phi = p4_Ks_star.Phi()-p4_l_star.Phi();
	delta_theta = p4_Ks_star.Theta()-p4_l_star.Theta();
	cout << "delta_p " << delta_p << " delta_phi " << delta_phi << " delta_theta " << delta_theta << endl;
	//******************boost to the reference frame of the detector****************************************//
	i++;    
	cout << "********************************************************************" << endl; 
   	h_delta_phi_Ks_l->Fill(delta_phi);
   	h_delta_theta_Ks_l->Fill(delta_theta);
   }//end for loop
 TFile *MyFile = new TFile("Simulation.root","NEW");
 h_delta_phi_Ks_l->Write();
 h_delta_theta_Ks_l->Write();
 MyFile->Write();
}
