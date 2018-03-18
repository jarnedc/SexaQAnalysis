/// \author Jarne De Clercq
// simple simulation for the kinematics of an S annihilating on a neutron and decaying to a Ks and a Lambda

//notations: all variables are in the detector rest frame, the once that are in the S+n rest frame are noted with "star", energy unit = GeV

#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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
    Double_t M_S_n = M_S + m_n;//neutron + S mass

    //TSalis distribution for the S momentum (https://arxiv.org/pdf/1401.4835.pdf)
    Double_t A = 1.;//factor in front of the TSalis function = g*V/(2*pi)^3, not important actually as it will not change the shape of the distribution, just normalisation
    Double_t q = 1.15; //q from the TSalis function 1.15 for CMS at 7TeV according to https://arxiv.org/pdf/1401.4835.pdf
    Double_t T = 0.080;//T in the TSalis function. 80MeV according to https://arxiv.org/pdf/1401.4835.pdf 
    TF1 *fTSalisSpt = new TF1("fTSalisSpt","[0]*x*pow(x*x+[1]*[1],0.5)*pow(1+([2]-1)*pow(x*x+[1]*[1],0.5)/[3],-[2]/([2]-1))",0,10);
    fTSalisSpt->SetParameter(0,A);
    fTSalisSpt->SetParameter(1,M_S);
    fTSalisSpt->SetParameter(2,q);
    fTSalisSpt->SetParameter(3,T);

    //plots
    //S partilce
    TH1F *h_S_pt = new TH1F("h_S_pt", "h_S_pt", 500, 0,10);
    TH1F *h_S_p = new TH1F("h_S_p", "h_S_p", 500, 0,10);
    TH1F *h_theta_S = new TH1F("h_theta_S", "h_theta_S", 100, -10,10);
    TH2F *h_S_p_xy = new TH2F("h_S_p_xy", "h_S_p_xy", 1000, -10, 10,1000, -10, 10);
    TH2F *h_S_p_xz = new TH2F("h_S_p_xz", "h_S_p_xz", 1000, -10, 10,1000, -10, 10);
    TH2F *h_S_p_yz = new TH2F("h_S_p_yz", "h_S_p_yz", 1000, -10, 10,1000, -10, 10);
    TH1F *h_E_S = new TH1F("h_E_S", "h_E_S", 200, 0, 20);
    //Ks in the star frame
    TH1F *h_E_Ks_star = new TH1F("h_E_Ks_star", "h_E_Ks_star", 200, 0, 20);
    TH1F *h_m_Ks_star = new TH1F("m_Ks_star","m_Ks_star", 200, 0,20);
    TH1F *h_E_l_star = new TH1F("h_E_l_star", "h_E_l_star", 200, 0, 20);
    TH1F *h_m_l_star = new TH1F("m_l_star","m_l_star", 200, 0,20);
    TH1F *h_p_Ks_star = new TH1F("h_p_Ks_star", "h_p_Ks_star", 200, 0, 20);
    TH1F *h_p_l_star = new TH1F("h_p_l_star", "h_p_l_star", 200, 0, 20);
    TH2F *h_Ks_star_p_xy = new TH2F("h_Ks_star_p_xy", "h_Ks_star_p_xy", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_star_p_xz = new TH2F("h_Ks_star_p_xz", "h_Ks_star_p_xz", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_star_p_yz = new TH2F("h_Ks_star_p_yz", "h_Ks_star_p_yz", 1000, -10, 10,1000, -10, 10);
    TH2F *h_l_star_p_xy = new TH2F("h_l_star_p_xy", "h_l_star_p_xy", 1000, -10, 10,1000, -10, 10);
    TH2F *h_l_star_p_xz = new TH2F("h_l_star_p_xz", "h_l_star_p_xz", 1000, -10, 10,1000, -10, 10);
    TH2F *h_l_star_p_yz = new TH2F("h_l_star_p_yz", "h_l_star_p_yz", 1000, -10, 10,1000, -10, 10);

    TH2F *h_Ks_l_star_p_x = new TH2F("h_Ks_l_star_p_x", "h_Ks_l_star_p_x", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_l_star_p_y = new TH2F("h_Ks_l_star_p_y", "h_Ks_l_star_p_y", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_l_star_p_z = new TH2F("h_Ks_l_star_p_z", "h_Ks_l_star_p_z", 1000, -10, 10,1000, -10, 10);

    TH2F *h_Ks_l_p_x = new TH2F("h_Ks_l_p_x", "h_Ks_l_p_x", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_l_p_y = new TH2F("h_Ks_l_p_y", "h_Ks_l_p_y", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_l_p_z = new TH2F("h_Ks_l_p_z", "h_Ks_l_p_z", 1000, -10, 10,1000, -10, 10);

    TH1F *h_delta_phi_Ks_l_star = new TH1F("h_delta_phi_Ks_l_star", "h_delta_phi_Ks_l_star", 200, -7, 7);
    TH1F *h_delta_theta_Ks_l_star = new TH1F("h_delta_theta_Ks_l_star", "h_delta_theta_Ks_l_star", 200, -4, 4);
 
    TH1F *h_delta_phi_Ks_l = new TH1F("h_delta_phi_Ks_l", "h_delta_phi_Ks_l", 100, -7, 7);
    TH1F *h_delta_theta_Ks_l = new TH1F("h_delta_theta_Ks_l", "h_delta_theta_Ks_l", 200, -4, 4);
    
    int i = 0; 
    bool verbose = false;  
 
    while(i<100000){

	//******************settings*****************************************//
    //pt distribution of the S particle:  https://arxiv.org/pdf/1401.4835.pdf
    //Double_t pt_S = fTSalisSpt->GetRandom(); 
    Double_t pt_S = 2; 
    //if(pt_S<2.)continue;
	//from https://arxiv.org/pdf/1002.0621.pdf, page 7 the eta distribution is relatively flat
    Double_t eta_S = random->Uniform(-2,2);
    //phi of the S, should be a uniform distribution in 0 to 2*pi
	Double_t phi_S = random->Uniform(2.*TMath::Pi());	

    //theta of the Ks, should be uniform distribution from 0 to pi at some point
	Double_t theta_Ks_star = random->Uniform(TMath::Pi());	
	//phi of the Ks, should be uniform distribution from 0 to 2*pi at some point
	Double_t phi_Ks_star = random->Uniform(2.*TMath::Pi());	

	//******************end settings*****************************************//

	if(verbose)cout << "---------------------------------------------------------------------------------" << endl;
	//******************calculate the S particle parameters******************************************//
	h_S_pt->Fill(pt_S);
	Double_t theta_S = 2*TMath::ATan(TMath::Exp(-eta_S));
	Double_t p_S = pt_S/TMath::Sin(theta_S);
    h_S_p->Fill(p_S);
    h_theta_S->Fill(theta_S);
	//p3_S this is the 3 vector of the S particle can be calculated from p_t_S and theta_S
	TVector3 p3_S(0.,0.,0.); 
	p3_S.SetMagThetaPhi(p_S,theta_S,phi_S); 
	h_S_p_xy->Fill(p3_S.Px(),p3_S.Py());
	h_S_p_xz->Fill(p3_S.Px(),p3_S.Pz());
	h_S_p_yz->Fill(p3_S.Py(),p3_S.Pz());
	//******************end calculate the S particle parameters******************************************//
		
	//the S hits on a neutron and transfers its momentum to the S+n system.
	
	//******************calculate the S+n particle parameters******************************************//
	//assume the S momentum is transfered to the S+n system, then the momenta are the same
	TVector3 p3_S_n =  p3_S;	
	Double_t E_S = pow(M_S_n*M_S_n+p_S*p_S,0.5);
	h_E_S->Fill(E_S);
	//******************end calculate the S+n particle parameters******************************************//
	
	//The S + n then decays to a Ks and Lambda. 
	//This is a 2 body decay (http://www.helsinki.fi/~www_sefo/phenomenology/Schlippe_relativistic_kinematics.pdf)
	
	//******************calculate the Ks and Lambda parameter****************************************//
	//E_Ks_star and p_Ks_star energy and momentum of the Ks in the rest frame of the S + n
	Double_t E_Ks_star = (M_S_n*M_S_n+m_Ks*m_Ks-m_l*m_l)/(2*M_S_n);
	Double_t p_Ks_star = pow(E_Ks_star*E_Ks_star - m_Ks*m_Ks,0.5);
	h_E_Ks_star->Fill(E_Ks_star);
	h_p_Ks_star->Fill(p_Ks_star);
	Double_t m_Ks_star = pow(E_Ks_star*E_Ks_star-p_Ks_star*p_Ks_star,0.5);
	h_m_Ks_star->Fill(m_Ks_star);
	//E_l_star and p_l_star energy and momentum of the l in the rest frame of the S + n
	Double_t E_l_star = (M_S_n*M_S_n+m_l*m_l-m_Ks*m_Ks)/(2*M_S_n);
	Double_t p_l_star = pow(E_l_star*E_l_star - m_l*m_l,0.5);
	h_E_l_star->Fill(E_l_star);
    h_p_l_star->Fill(p_l_star);
    Double_t m_l_star = pow(E_l_star*E_l_star-p_l_star*p_l_star,0.5);
	h_m_l_star->Fill(m_l_star);
	//p4_Ks_star 4 momentum vector of the Ks 
	TLorentzVector p4_Ks_star(0.,0.,0.,0.);
	p4_Ks_star.SetPtEtaPhiE(p_Ks_star*TMath::Sin(theta_Ks_star),-TMath::Log(TMath::Tan(theta_Ks_star/2.)),phi_Ks_star,E_Ks_star);	
	TLorentzVector p4_l_star(0.,0.,0.,0.);
	p4_l_star.SetPxPyPzE(-p4_Ks_star.Px(),-p4_Ks_star.Py(),-p4_Ks_star.Pz(),E_l_star);
	if(verbose)cout << "theta of the Ks: " << p4_Ks_star.Theta() << endl;
	if(verbose)cout << "theta of the l: " << p4_l_star.Theta() << endl;

	if(verbose)cout << "phi of the Ks: " << p4_Ks_star.Phi() << endl;
	if(verbose)cout << "phi of the l: " << p4_l_star.Phi() << endl;

 	if(verbose)cout << "delta phi: " << p4_Ks_star.Phi()-p4_l_star.Phi() << endl;
	if(verbose)cout << "delta theta: " << p4_Ks_star.Theta()-p4_l_star.Theta() << endl;

	h_Ks_star_p_xy->Fill(p4_Ks_star.Px(),p4_Ks_star.Py());
    h_Ks_star_p_xz->Fill(p4_Ks_star.Px(),p4_Ks_star.Pz());
    h_Ks_star_p_yz->Fill(p4_Ks_star.Py(),p4_Ks_star.Pz());
	
	h_l_star_p_xy->Fill(p4_l_star.Px(),p4_l_star.Py());
    h_l_star_p_xz->Fill(p4_l_star.Px(),p4_l_star.Pz());
    h_l_star_p_yz->Fill(p4_l_star.Py(),p4_l_star.Pz());

	h_Ks_l_star_p_x->Fill(p4_Ks_star.Px(),p4_l_star.Px());
    h_Ks_l_star_p_y->Fill(p4_Ks_star.Py(),p4_l_star.Py());
    h_Ks_l_star_p_z->Fill(p4_Ks_star.Pz(),p4_l_star.Pz());
	h_delta_phi_Ks_l_star->Fill(p4_Ks_star.Phi()-p4_l_star.Phi());
   	h_delta_theta_Ks_l_star->Fill(p4_Ks_star.Theta()-p4_l_star.Theta());

	//******************end calculate the Ks and Lambda parameter****************************************//

	//you have now the vector of the Ks and the Lambda in the rest frame of the S+n. Now will boost these vectors to the reference
	//frame where the S+n is actually moving with momentum p3_S, this is the detector reference frame
	
	
	//******************boost to the reference frame of the detector****************************************//
	//boosting the p4_Ks_star and p4_l_star allong the S+n system momentum	
	p4_Ks_star.Boost(p3_S_n.x()/E_S,p3_S_n.y()/E_S,p3_S_n.z()/E_S); //beta = p/E
	p4_l_star.Boost(p3_S_n.x()/E_S,p3_S_n.y()/E_S,p3_S_n.z()/E_S); //beta = p/E

	if(verbose)cout << "After boost" << endl;
	if(verbose)cout << "theta of the Ks: " << p4_Ks_star.Theta() << endl;
	if(verbose)cout << "theta of the l: " << p4_l_star.Theta() << endl;

	if(verbose)cout << "phi of the Ks: " << p4_Ks_star.Phi() << endl;
	if(verbose)cout << "phi of the l: " << p4_l_star.Phi() << endl;

 	if(verbose)cout << "delta phi: " << p4_Ks_star.Phi()-p4_l_star.Phi() << endl;
	if(verbose)cout << "delta theta: " << p4_Ks_star.Theta()-p4_l_star.Theta() << endl; 	

        h_Ks_l_p_x->Fill(p4_Ks_star.Px(),p4_l_star.Px());
        h_Ks_l_p_y->Fill(p4_Ks_star.Py(),p4_l_star.Py());
        h_Ks_l_p_z->Fill(p4_Ks_star.Pz(),p4_l_star.Pz());


	//******************end boost to the reference frame of the detector****************************************//
	i++;
    
    Double_t delta_phi_Ks_l = p4_Ks_star.Phi()-p4_l_star.Phi();
    if(delta_phi_Ks_l < -TMath::Pi()) delta_phi_Ks_l = delta_phi_Ks_l + TMath::Pi();
    else if(delta_phi_Ks_l > TMath::Pi()) delta_phi_Ks_l = delta_phi_Ks_l - TMath::Pi();

   	h_delta_phi_Ks_l->Fill(delta_phi_Ks_l);
   	h_delta_theta_Ks_l->Fill(p4_Ks_star.Theta()-p4_l_star.Theta());

    
	//******************check the above calculation by calculating the invariant mass of the S****************************************//
    //Double_t M_S_n_check = m_l_star*m_l_star + m_Ks_star*m_Ks_star + 2*(p4_Ks_star.E()*p4_l_star.E()-p4_Ks_star*p4_l_star);
    Double_t M_S_n_check = pow((p4_Ks_star+p4_l_star)*(p4_Ks_star+p4_l_star),0.5);
    if(fabs(M_S_n-M_S_n_check)>0.001) cout << "Something is wrong in your calculation. Invariant mass of the Ks and Lambda does not match m_S + m_neutron" << endl;
    //******************end check the above calculation by calculating the invariant mass of the S****************************************//
	
	i++;
    }//end for loop
 

 TFile *MyFile = new TFile("Simulation.root","RECREATE");
 h_S_pt->Write();
 h_S_p->Write();
 h_theta_S->Write();
 h_S_p_xy->Write();
 h_S_p_xz->Write();
 h_S_p_yz->Write();
 h_E_S->Write();
 
 h_E_Ks_star->Write();
 h_m_Ks_star->Write();
 h_E_l_star->Write();
 h_m_l_star->Write();
 h_p_Ks_star->Write();
 h_p_l_star->Write();
 h_Ks_star_p_xy->Write();
 h_Ks_star_p_xz->Write();
 h_Ks_star_p_yz->Write();
 h_l_star_p_xy->Write();
 h_l_star_p_xz->Write();
 h_l_star_p_yz->Write();

 h_Ks_l_star_p_x->Write();
 h_Ks_l_star_p_y->Write();
 h_Ks_l_star_p_z->Write();
 
 h_Ks_l_p_x->Write();
 h_Ks_l_p_y->Write();
 h_Ks_l_p_z->Write();
 
 h_delta_phi_Ks_l_star->Write();
 h_delta_theta_Ks_l_star->Write();

 h_delta_phi_Ks_l->Write();
 h_delta_theta_Ks_l->Write();
 MyFile->Write();
 MyFile->Write();
}
