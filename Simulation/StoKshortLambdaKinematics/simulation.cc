/// \author Jarne De Clercq
// simple simulation for the kinematics of an S annihilating on a neutron and decaying to a Ks and a Lambda

//notations: all variables are in the detector rest frame, the once that are in the S+n rest frame are noted with "star", energy unit = GeV
//to run this for example: root -l -> in interactive mode (and change the name of particle variable below):

//.x simulation.cc++("Xi_TSalis_polarised", "Xi1820polarised.root",  100000, "polarisedJ32M12" )
//.x simulation.cc++("Xi_TSalis_unpolarised", "Xi1820unpolarised.root", 100000, "unpolarised" )

//.x simulation.cc++("S_ptTSalis", "S.root", 100000, "unpolarised" )

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
    

    //settings

    string particle = "S";
    //constants
    Double_t m_Ks = 0.497611;//Kshort mass
    Double_t m_l = 1.115683;//Lambda mass

    //cuts
    Double_t max_eta = 2.5;
    Double_t min_pt_Ks = 0.9;
    Double_t min_pt_l = 1.5;

    //plots
    //settings
    TH1F *h_u = new TH1F("h_u","h_u",1000,0,5);
    TH1F *h_cos_theta_Ks = new TH1F("h_cos_theta_Ks","h_cos_theta_Ks",1000,-2,2);
    TH2F *h2_cos_theta_Ks = new TH2F("h2_cos_theta_Ks","h2_cos_theta_Ks",100,-5,5,1000,-2,2);
    TH2F *h2_theta_Ks = new TH2F("h2_theta_Ks","h2_theta_Ks",10000,-5,5,1000,-2.*TMath::Pi(),2.*TMath::Pi());
    TH1F *h_J32_to_M12_theta = new TH1F("h_J32_to_M12_theta","h_J32_to_M12_theta",100,-5,5);
    TH1F *h_sqrt = new TH1F("h_sqrt","h_sqrt",10000,-5,5);
    TH2F *h2_sqrt = new TH2F("h2_sqrt","h2_sqrt;u;sqrt((u-1)/3)",1000,0,5,1000,-5,5);
    TH1F *h_theta_Ks = new TH1F("h_theta_Ks","h_theta_Ks",1000,-2*TMath::Pi(),2*TMath::Pi());

    //start particle plots
    TH1F *h_S_pt = new TH1F("h_Start_pt", "h_Start_pt", 50, 0,10);
    TH1F *h_S_p = new TH1F(("h_"+particle+"_p").c_str(), ("h_"+particle+"_p").c_str(), 500, 0,10);
    TH1F *h_theta_S = new TH1F("h_theta_S", "h_theta_S", 100, -10,10);
    TH1F *h_eta_S = new TH1F("h_eta_S", "h_eta_S", 1200, -12,12);
    TH1F *h_S_p_x = new TH1F(("h_"+particle+"_p_x").c_str(), ("h_"+particle+"_p_x").c_str(), 1000, -10, 10);
    TH2F *h2_S_p_x = new TH2F(("h2_"+particle+"_p_x").c_str(), ("h2_"+particle+"_p_x").c_str(),1000,-2.*TMath::Pi(),2.*TMath::Pi(), 1000, -10, 10);
    TH1F *h_S_p_y = new TH1F(("h_"+particle+"_p_y").c_str(), ("h_"+particle+"_p_y").c_str(), 1000, -10, 10);
    TH2F *h2_S_p_y = new TH2F(("h2_"+particle+"_p_y").c_str(), ("h2_"+particle+"_p_y").c_str(), 1000, -2.*TMath::Pi(),2.*TMath::Pi(), 1000, -10, 10);
    TH1F *h_S_p_z = new TH1F(("h_"+particle+"_p_z").c_str(), ("h_"+particle+"_p_z").c_str(), 1000, -10, 10);
    TH2F *h2_S_p_z = new TH2F(("h2_"+particle+"_p_z").c_str(), ("h2_"+particle+"_p_z").c_str(), 1000, -2.*TMath::Pi(),2.*TMath::Pi(), 1000, -10, 10);
    TH2F *h2_eta_S_p_z = new TH2F(("h2_eta_"+particle+"_p_z").c_str(), ("h2_eta_"+particle+"_p_z").c_str(), 1000, -2.*TMath::Pi(),2.*TMath::Pi(), 1000, -10, 10);
    TH1F *h_S_P = new TH1F(("h_"+particle+"_P").c_str(), ("h_"+particle+"_P").c_str(), 1000, -10, 10);
    TH2F *h_S_p_xy = new TH2F(("h_"+particle+"_p_xy").c_str(), ("h_"+particle+"_p_xy;px;py").c_str(), 1000, -10, 10,1000, -10, 10);
    TH2F *h_S_p_xz = new TH2F(("h_"+particle+"_p_xz").c_str(), ("h_"+particle+"_p_xz;px;pz").c_str(), 1000, -10, 10,1000, -10, 10);
    TH2F *h_S_p_yz = new TH2F(("h_"+particle+"_p_yz").c_str(), ("h_"+particle+"_p_yz;py;pz").c_str(), 1000, -10, 10,1000, -10, 10);

    TH1F *h_E_S_n = new TH1F(("h_E_"+particle+"_n").c_str(), ("h_E_"+particle+"_n").c_str(), 2000, 0, 20);
    TH1F *h_M_S_n = new TH1F(("h_M_"+particle+"_n").c_str(), ("h_M_"+particle+"_n").c_str(), 200, 0, 20);


    //Ks and l in the star frame
    TH1F *h_E_Ks_star = new TH1F("h_E_Ks_star", "h_E_Ks_star", 200, 0, 20);
    TH1F *h_m_Ks_star = new TH1F("m_Ks_star","m_Ks_star", 200, 0,20);
    TH1F *h_E_l_star = new TH1F("h_E_l_star", "h_E_l_star", 200, 0, 20);
    TH1F *h_m_l_star = new TH1F("m_l_star","m_l_star", 200, 0,20);
    TH1F *h_p_Ks_star = new TH1F("h_p_Ks_star", "h_p_Ks_star", 200, 0, 20);
    TH1F *h_p_l_star = new TH1F("h_p_l_star", "h_p_l_star", 200, 0, 20);

    TH1F * h_Ks_star_p_x =  new TH1F("h_Ks_star_p_x", "h_Ks_star_p_x", 200, -3, 3);
    TH1F * h_Ks_star_p_y =  new TH1F("h_Ks_star_p_y", "h_Ks_star_p_y", 200, -3, 3);
    TH1F * h_Ks_star_p_z =  new TH1F("h_Ks_star_p_z", "h_Ks_star_p_z", 200, -3, 3);
    TH1F * h_Ks_star_p =  new TH1F("h_Ks_star_p", "h_Ks_star_p", 200, -3, 3);
    
    TH1F * h_l_star_p_x =  new TH1F("h_l_star_p_x", "h_l_star_p_x", 200, -3, 3);
    TH1F * h_l_star_p_y =  new TH1F("h_l_star_p_y", "h_l_star_p_y", 200, -3, 3);
    TH1F * h_l_star_p_z =  new TH1F("h_l_star_p_z", "h_l_star_p_z", 200, -3, 3);
    TH1F * h_l_star_p =  new TH1F("h_l_star_p", "h_l_star_p", 200, -3, 3);

   
    TH2F *h_Ks_star_p_xy_norm = new TH2F("h_Ks_star_p_xy_norm", "h_Ks_star_p_xy_norm", 1000, -1.1, 1.1,1000, -1.1, 1.1);
    TH2F *h_Ks_star_p_xz_norm = new TH2F("h_Ks_star_p_xz_norm", "h_Ks_star_p_xz_norm", 1000, -1.1, 1.1,1000, -1.1, 1.1);
    TH2F *h_Ks_star_p_yz_norm = new TH2F("h_Ks_star_p_yz_norm", "h_Ks_star_p_yz_norm", 1000, -1.1, 1.1,1000, -1.1, 1.1);
    TH2F *h_l_star_p_xy = new TH2F("h_l_star_p_xy", "h_l_star_p_xy", 1000, -10, 10,1000, -10, 10);
    TH2F *h_l_star_p_xz = new TH2F("h_l_star_p_xz", "h_l_star_p_xz", 1000, -10, 10,1000, -10, 10);
    TH2F *h_l_star_p_yz = new TH2F("h_l_star_p_yz", "h_l_star_p_yz", 1000, -10, 10,1000, -10, 10);

    TH2F *h_Ks_l_star_p_x = new TH2F("h_Ks_l_star_p_x", "h_Ks_l_star_p_x", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_l_star_p_y = new TH2F("h_Ks_l_star_p_y", "h_Ks_l_star_p_y", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_l_star_p_z = new TH2F("h_Ks_l_star_p_z", "h_Ks_l_star_p_z", 1000, -10, 10,1000, -10, 10);

    TH1F *h_delta_phi_Ks_l_star = new TH1F("h_delta_phi_Ks_l_star", "h_delta_phi_Ks_l_star", 100, -TMath::Pi(), TMath::Pi());
    TH1F *h_delta_theta_Ks_l_star = new TH1F("h_delta_theta_Ks_l_star", "h_delta_theta_Ks_l_star", 100, -TMath::Pi(), TMath::Pi());
    TH2F *h2_delta_theta_Ks_l_star = new TH2F("h2_delta_theta_Ks_l_star", "h2_delta_theta_Ks_l_star", 100, -2, 2, 200, -4, 4);
    TH2F *h2_delta_eta_Ks_l_star = new TH2F("h2_delta_eta_Ks_l_star", "h2_delta_eta_Ks_l_star", 100, -2, 2, 200, -4, 4);
    TH1F *h_sum_theta_Ks_l_star = new TH1F("h_sum_theta_Ks_l_star", "h_sum_theta_Ks_l_star", 200, -4, 4);

    Double_t p_xyz_range = 6.; //1.1
    TH2F *h_Ks_star_p_xy = new TH2F("h_Ks_star_p_xy", "h_Ks_star_p_xy", 1000, -p_xyz_range, p_xyz_range,1000, -p_xyz_range, p_xyz_range);
    TH1D *h_Ks_star_p_xy_proj_x =  new TH1D("h_Ks_star_p_xy_proj_x", "h_Ks_star_p_xy_proj_x", 200, -20, 20);
    TH2F *h_Ks_star_p_xz = new TH2F("h_Ks_star_p_xz", "h_Ks_star_p_xz", 1000, -p_xyz_range, p_xyz_range,1000, -p_xyz_range, p_xyz_range);
    TH2F *h_Ks_star_p_yz = new TH2F("h_Ks_star_p_yz", "h_Ks_star_p_yz", 1000, -p_xyz_range, p_xyz_range,1000, -p_xyz_range, p_xyz_range);


    //after boost
    TH1F *h_Ks_p_x = new TH1F("h_Ks_p_x", "h_Ks_p_x", 1000, -10, 10);
    TH1F *h_Ks_p_y = new TH1F("h_Ks_p_y", "h_Ks_p_y", 1000, -10, 10);
    TH1F *h_Ks_p_z = new TH1F("h_Ks_p_z", "h_Ks_p_z", 1000, -10, 10);
    TH1F *h_Ks_pt = new TH1F("h_Ks_pt", "h_Ks_pt", 1000, -10, 10);

    TH1F *h_l_p_x = new TH1F("h_l_p_x", "h_l_p_x", 1000, -10, 10);
    TH1F *h_l_p_y = new TH1F("h_l_p_y", "h_l_p_y", 1000, -10, 10);
    TH1F *h_l_p_z = new TH1F("h_l_p_z", "h_l_p_z", 1000, -10, 10);
    TH1F *h_l_pt = new TH1F("h_l_pt", "h_l_pt", 1000, -10, 10);

    TH2F *h_Ks_l_p_x = new TH2F("h_Ks_l_p_x", "h_Ks_l_p_x", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_l_p_y = new TH2F("h_Ks_l_p_y", "h_Ks_l_p_y", 1000, -10, 10,1000, -10, 10);
    TH2F *h_Ks_l_p_z = new TH2F("h_Ks_l_p_z", "h_Ks_l_p_z", 1000, -10, 10,1000, -10, 10);

    TH1F *h_p_l_trans = new TH1F("h_p_l_trans", "h_p_l_trans; Lambda transversal momentum along Start particle direction", 200, -5, 5);
    TH1F *h_p_Ks_trans = new TH1F("h_p_Ks_trans", "h_p_Ks_trans; Ks transversal momentum along Start particle direction", 200, -5, 5);
    TH1F *h_p_l_long = new TH1F("h_p_l_long", "h_p_l_long; Lambda longitudinal momentum along Start particle direction", 200, -5, 5);
    TH1F *h_p_Ks_long = new TH1F("h_p_Ks_long", "h_p_Ks_long; Ks longitudnal momentum along Start particle direction", 200, -5, 5);

    TH2F *h2_p_Ks_l_trans = new TH2F("h2_p_Ks_l_trans","h2_p_Ks_l_trans;transversal momentum Ks wrt Start particle;transversal momentum Lambda wrt Start particle",1000,-1,1,1000,-1,1);
    TH2F *h2_p_Ks_l_long = new TH2F("h2_p_Ks_l_long","h2_p_Ks_l_long;longitudnal momentum Ks wrt Start particle;longitudnal momentum Lambda wrt Start particle",100,-10,10,100,-10,10);

    TH1F *h_delta_phi_Ks_l = new TH1F("h_delta_phi_Ks_l", "h_delta_phi_Ks_l", 100, -4, 4);
    TH1F *h_delta_theta_Ks_l = new TH1F("h_delta_theta_Ks_l", "h_delta_theta_Ks_l", 100, -TMath::Pi(), TMath::Pi());
    TH1F *h_delta_eta_Ks_l = new TH1F("h_delta_eta_Ks_l", "h_delta_eta_Ks_l", 100, -TMath::Pi(), TMath::Pi());
    TH1F *h_delta_R_Ks_l = new TH1F("h_delta_R_Ks_l", "h_delta_R_Ks_l", 200, 0, 10);
   
   
    TH1F *h_S_pt_before_cuts = new TH1F("h_S_pt_before_cuts", "h_S_pt_before_cuts; pT(S) (GeV)", 200, 0, 10);
    TH1F *h_S_eta_before_cuts = new TH1F("h_S_eta_before_cuts", "h_S_eta_before_cuts; eta(S)", 200, -15, 15);
    TH1F *h_Ks_pt_before_cuts = new TH1F("h_Ks_pt_before_cuts", "h_Ks_pt_before_cuts; pT(Ks) (GeV)", 200, 0, 10);
    TH1F *h_Ks_eta_before_cuts = new TH1F("h_Ks_eta_before_cuts", "h_Ks_eta_before_cuts; eta(Ks)", 200, -15, 15);
    TH1F *h_l_pt_before_cuts = new TH1F("h_l_pt_before_cuts", "h_l_pt_before_cuts;pT(L) (GeV)", 200, 0, 10);
    TH1F *h_l_eta_before_cuts = new TH1F("h_l_eta_before_cuts", "h_l_eta_before_cuts; eta(L)", 200, -15, 15);

    TH1F *h_S_pt_after_cuts = new TH1F("h_S_pt_after_cuts", "h_S_pt_after_cuts; pT(S) (GeV)", 200, 0, 10);
    TH1F *h_S_p_after_cuts = new TH1F("h_S_p_after_cuts", "h_S_p_after_cuts; p(S) (GeV)", 200, 0, 10);
    TH1F *h_S_eta_after_cuts = new TH1F("h_S_eta_after_cuts", "h_S_eta_after_cuts; eta(S)", 200, -15, 15);
    TH1F *h_Ks_pt_after_cuts = new TH1F("h_Ks_pt_after_cuts", "h_Ks_pt_after_cuts; pT(Ks) (GeV)", 200, 0, 10);
    TH1F *h_Ks_eta_after_cuts = new TH1F("h_Ks_eta_after_cuts", "h_Ks_eta_after_cuts; eta(Ks)", 200, -15, 15);
    TH1F *h_l_pt_after_cuts = new TH1F("h_l_pt_after_cuts", "h_l_pt_after_cuts;pT(L) (GeV)", 200, 0, 10);
    TH1F *h_l_eta_after_cuts = new TH1F("h_l_eta_after_cuts", "h_l_eta_after_cuts; eta(L)", 200, -15, 15);

    Double_t delta_phi_detla_theta_range = 4.; //2
    TH2F *h_pt_S_delta_phi_Ks_l = new TH2F(("h_pt_"+particle+"_delta_phi_Ks_l").c_str(),("h_pt_"+particle+"_delta_phi_Ks_l;pt_"+particle+";delta_phi_Ks_l").c_str(),100,0,10,100,-delta_phi_detla_theta_range,delta_phi_detla_theta_range);
    TH2F *h_p_S_delta_phi_Ks_l = new TH2F(("h_p_"+particle+"_delta_phi_Ks_l").c_str(),("h_p_"+particle+"_delta_phi_Ks_l").c_str(),100,0,10,100,-7,7);
    TH2F *h_pt_S_delta_theta_Ks_l = new TH2F(("h_pt_"+particle+"_delta_theta_Ks_l").c_str(),("h_pt_"+particle+"_delta_theta_Ks_l;pt_"+particle+";delta_theta_Ks_l").c_str(),100,0,10,100,-delta_phi_detla_theta_range,delta_phi_detla_theta_range);
    TH2F *h_p_S_delta_theta_Ks_l = new TH2F(("h_p_"+particle+"_delta_theta_Ks_l").c_str(),("h_p_"+particle+"_delta_theta_Ks_l").c_str(),100,0,10,100,-7,7);
    TH2F *h_delta_phi_delta_theta_Ks_l = new TH2F("h_delta_phi_delta_theta_Ks_l","h_delta_phi_delta_theta_Ks_l;delta_phi_Ks_l;delta_theta_Ks_l",100,-delta_phi_detla_theta_range,delta_phi_detla_theta_range,100,-delta_phi_detla_theta_range,delta_phi_detla_theta_range);
    TH2F *h_delta_phi_delta_eta_Ks_l = new TH2F("h_delta_phi_delta_eta_Ks_l","h_delta_phi_delta_eta_Ks_l;delta_phi_Ks_l;delta_eta_Ks_l",100,-delta_phi_detla_theta_range,delta_phi_detla_theta_range,100,-delta_phi_detla_theta_range,delta_phi_detla_theta_range);
 
    Double_t ArmPod_alpha_range = 5.; //1
    TH2F *h2_ArmPod = new TH2F("h2_ArmPod","h2_ArmPod;alpha;pT(Ks,Lambda)",100,-ArmPod_alpha_range,ArmPod_alpha_range,100,0,5);
    TH2F *h2_ArmPod_Ks = new TH2F("h2_ArmPod_Ks","h2_ArmPod_Ks;alpha;pT(Ks)",100,-ArmPod_alpha_range,ArmPod_alpha_range,100,0,5);
    TH2F *h2_ArmPod_Lambda = new TH2F("h2_ArmPod_Lambda","h2_ArmPod_Lambda;alpha;pT(Lambda)",100,-ArmPod_alpha_range,ArmPod_alpha_range,100,0,5);

    TH1F *h_M_Start_n_check =  new TH1F(("h_M_"+particle+"_n_check").c_str(), ("h_M_"+particle+"_n_check").c_str(), 200, 0, 20);    
    TH1F *h_M_Start_n_Inv_Mass_calc =  new TH1F(("h_M_"+particle+"_n_Inv_Mass_calc").c_str(), ("h_M_"+particle+"_n_Inv_Mass_calc").c_str(),200,0,20);

    

    TH2F *h2_check_intEff_Ks = new TH2F("h2_check_intEff_Ks","h2_check_intEff_KS;pT;pTeff",1000,0,10,1000,0,1);
    TH2F *h2_check_intEff_L0 = new TH2F("h2_check_intEff_L0","h2_check_intEff_L0;pT;pTeff",1000,0,10,1000,0,1);
    
    TEfficiency* hEff_Ks = new TEfficiency("hEff_Ks","hEff_Ks;pT_S;Ks overall efficiency",50,0,10);
    TEfficiency* hEff_L0 = new TEfficiency("hEff_L0","hEff_L0;pT_S;L0 overall efficiency",50,0,10);
    TEfficiency* hEff_L0_Ks = new TEfficiency("hEff_L0_Ks","hEff_L0_Ks;pT_S;L0 & Ks overall efficiency",50,0,10);





void fillHistosSettings(Double_t u, Double_t theta_Ks_star){
    h_u->Fill(u);
    h_theta_Ks->Fill(theta_Ks_star);
    h2_theta_Ks->Fill(u,theta_Ks_star);
    h_cos_theta_Ks->Fill(TMath::Cos(theta_Ks_star));
    h2_cos_theta_Ks->Fill(u,TMath::Cos(theta_Ks_star));
    h_J32_to_M12_theta->Fill(1+3*TMath::Cos(theta_Ks_star)*TMath::Cos(theta_Ks_star));
}

TF1* TSalis(Double_t M_Start){
    //TSalis distribution for the S momentum (https://arxiv.org/pdf/1401.4835.pdf)
    Double_t A = 1.;//factor in front of the TSalis function = g*V/(2*pi)^3, not important actually as it will not change the shape of the distribution, just normalisation
    Double_t q = 1.15; //q from the TSalis function 1.15 for CMS at 7TeV according to https://arxiv.org/pdf/1401.4835.pdf
    Double_t T = 0.080;//T in the TSalis function. 80MeV according to https://arxiv.org/pdf/1401.4835.pdf 
    TF1 *fTSalisSpt = new TF1("fTSalisSpt","[0]*x*pow(x*x+[1]*[1],0.5)*pow(1+([2]-1)*pow(x*x+[1]*[1],0.5)/[3],-[2]/([2]-1))",0,10);
    fTSalisSpt->SetParameter(0,A);
    fTSalisSpt->SetParameter(1,M_Start);
    fTSalisSpt->SetParameter(2,q);
    fTSalisSpt->SetParameter(3,T);

    return fTSalisSpt;
}

TF1* EtaDistrHadr(){
    //eta distribution from: http://iopscience.iop.org/article/10.1209/0295-5075/95/61001/pdf --> distribution for 14TeV in figure 14 fitted with a 6th order polynomial (had to read the values from the plot, put them in excel and do the fit)

    TF1 *fEtaDistrHadr = new TF1("fEtaDistrHadr","-0.0009*pow(x,5)+0.0221*pow(x,4)-0.1614*pow(x,3)+0.2221*pow(x,2)+0.2747*x+6.6822",0.,10.);

    return fEtaDistrHadr;
}

Double_t return_M_Start(string particle){
        Double_t M_Start =-99999.;

        if(particle == "S"){
            cout << "Running for " << particle << endl;
            M_Start = 2.;
        }
        else if(particle == "Xi1820"){
            cout << "Running for " << particle << endl;
            M_Start = 1.82; //Xi particle mass
        }
        else{
            cout << "particle not known" << endl;
            return 1;
        }
        return M_Start;
}

Double_t return_M_n(string particle){
        Double_t m_n =-99999.;
        
        if(particle == "S"){
            m_n = 0.9395654133;
        }
        else if(particle == "Xi1820"){
            m_n = 0.;
        }
        else{
            cout << "particle not known" << endl;
            return 1;
        }
        return m_n;
}

TEfficiency* getRecoEfficiencyKs_L0(string plot){
//loading TEfficiencies
    TFile *efficiency_file = new TFile("./reco_efficiency/HQ_plots_minbias_10M_events_skimmed_fixedGEANT_c123.root");
    TEfficiency *TEff = (TEfficiency*)efficiency_file->Get(plot.c_str());
    return TEff;
}

Double_t return_argument_Acos_theta_Ks_angle(TRandom *random, string polarization){
    Double_t u = -9999;
    Double_t theta_Ks_star = -9999;   
    if(polarization == "unpolarised"){
       Double_t u = random->Uniform(-1.,1.);
       theta_Ks_star = TMath::ACos(u);  
       return u; 
    }

    //theta of the Ks, distributed non isotropically for a polarised decay
    else if(polarization == "polarisedJ32M12"){
        Double_t Damocles = random->Uniform(0.,1.);
        u = random->Uniform(1.,4.);
        Double_t sqrt = pow(((u-1.)/3.),0.5);
        if(Damocles > 0.5) sqrt = -sqrt;
        theta_Ks_star = TMath::ACos(sqrt);

        h_sqrt->Fill(sqrt);
        h2_sqrt->Fill(u,sqrt);
        return sqrt;
    }
    else {
        cout << "I do not know the polarisation";
        return -9999;
    }

}

void fillHistosSParticle (Double_t pt_S, Double_t p_S, Double_t theta_S, Double_t phi_S, TVector3 p3_S, Double_t eta_S){

    h_S_pt->Fill(pt_S);
    h_S_p->Fill(p_S);
    h_eta_S->Fill(eta_S);
    h_theta_S->Fill(theta_S);
    h_S_p_x->Fill(p3_S.Px());
    h2_S_p_x->Fill(phi_S,p3_S.Px());
    h_S_p_y->Fill(p3_S.Py());
    h2_S_p_y->Fill(phi_S,p3_S.Py());
    h_S_p_z->Fill(p3_S.Pz());
    h2_S_p_z->Fill(theta_S,p3_S.Pz());
    h2_eta_S_p_z->Fill(eta_S,p3_S.Pz());
    h_S_P->Fill(p3_S.Mag());

    h_S_p_xy->Fill(p3_S.Px(),p3_S.Py());
    h_S_p_xz->Fill(p3_S.Px(),p3_S.Pz());
    h_S_p_yz->Fill(p3_S.Py(),p3_S.Pz());

}

void fillHistosS_n(Double_t E_S_n, Double_t M_Start_n){

    h_E_S_n->Fill(E_S_n);
    h_M_S_n->Fill(M_Start_n);

}


 void fillHistosStar(Double_t E_Ks_star, Double_t p_Ks_star, Double_t m_Ks_star, Double_t E_l_star, Double_t p_l_star, Double_t m_l_star, TLorentzVector p4_Ks_star, TLorentzVector p4_l_star, Double_t argument_cos){
    h_E_Ks_star->Fill(E_Ks_star);
    h_p_Ks_star->Fill(p_Ks_star);

    h_m_Ks_star->Fill(m_Ks_star);   

    h_E_l_star->Fill(E_l_star);
    h_p_l_star->Fill(p_l_star);  

    h_m_l_star->Fill(m_l_star);


    h_Ks_star_p_x->Fill(p4_Ks_star.Px());
    h_Ks_star_p_y->Fill(p4_Ks_star.Py());
    h_Ks_star_p_z->Fill(p4_Ks_star.Pz());
    h_Ks_star_p->Fill(p4_Ks_star.P());

    h_l_star_p_x->Fill(p4_l_star.Px());
    h_l_star_p_y->Fill(p4_l_star.Py());
    h_l_star_p_z->Fill(p4_l_star.Pz());
    h_l_star_p->Fill(p4_l_star.P());

    h_Ks_star_p_xy->Fill(p4_Ks_star.Px(),p4_Ks_star.Py());
    h_Ks_star_p_xz->Fill(p4_Ks_star.Px(),p4_Ks_star.Pz());
    h_Ks_star_p_yz->Fill(p4_Ks_star.Py(),p4_Ks_star.Pz());
    

    h_Ks_star_p_xy_norm->Fill(p4_Ks_star.Px()/p4_Ks_star.P(),p4_Ks_star.Py()/p4_Ks_star.P());
    h_Ks_star_p_xz_norm->Fill(p4_Ks_star.Px()/p4_Ks_star.P(),p4_Ks_star.Pz()/p4_Ks_star.P());
    h_Ks_star_p_yz_norm->Fill(p4_Ks_star.Py()/p4_Ks_star.P(),p4_Ks_star.Pz()/p4_Ks_star.P());
    
    h_l_star_p_xy->Fill(p4_l_star.Px(),p4_l_star.Py());
    h_l_star_p_xz->Fill(p4_l_star.Px(),p4_l_star.Pz());
    h_l_star_p_yz->Fill(p4_l_star.Py(),p4_l_star.Pz());

    h_Ks_l_star_p_x->Fill(p4_Ks_star.Px(),p4_l_star.Px());
    h_Ks_l_star_p_y->Fill(p4_Ks_star.Py(),p4_l_star.Py());
    h_Ks_l_star_p_z->Fill(p4_Ks_star.Pz(),p4_l_star.Pz());
    h_delta_phi_Ks_l_star->Fill(p4_Ks_star.Phi()-p4_l_star.Phi());
    h_delta_theta_Ks_l_star->Fill(p4_Ks_star.Theta()-p4_l_star.Theta());
    h2_delta_theta_Ks_l_star->Fill(argument_cos,p4_Ks_star.Theta()-p4_l_star.Theta());
    h2_delta_eta_Ks_l_star->Fill(argument_cos,p4_Ks_star.Eta()-p4_l_star.Eta());
    h_sum_theta_Ks_l_star->Fill(p4_Ks_star.Theta()+p4_l_star.Theta());
 }

 void  fillHistosAfterBoost(Double_t pt_S, Double_t eta_S, TLorentzVector p4_Ks_star, TLorentzVector p4_l_star){

    h_S_pt_before_cuts->Fill(pt_S);
    h_S_eta_before_cuts->Fill(eta_S);
    h_Ks_pt_before_cuts->Fill(p4_Ks_star.Pt());
    h_Ks_eta_before_cuts->Fill(p4_Ks_star.Eta());
    h_l_pt_before_cuts->Fill(p4_l_star.Pt());
    h_l_eta_before_cuts->Fill(p4_l_star.Eta());     

}

void fillHistosAfterCuts(TLorentzVector p4_Ks_star, TLorentzVector p4_l_star, Double_t pt_S, TVector3 p3_S, Double_t eta_S, Double_t delta_phi_Ks_l, Double_t delta_theta_Ks_l, Double_t delta_eta_Ks_l, Double_t delta_R_Ks_l){
    h_Ks_p_x->Fill(p4_Ks_star.Px());
    h_Ks_p_y->Fill(p4_Ks_star.Py());
    h_Ks_p_z->Fill(p4_Ks_star.Pz());
    h_Ks_pt->Fill(p4_Ks_star.Pt());

    h_l_p_x->Fill(p4_l_star.Px());
    h_l_p_y->Fill(p4_l_star.Py());
    h_l_p_z->Fill(p4_l_star.Pz());
    h_l_pt->Fill(p4_l_star.Pt());

    h_Ks_l_p_x->Fill(p4_Ks_star.Px(),p4_l_star.Px());
    h_Ks_l_p_y->Fill(p4_Ks_star.Py(),p4_l_star.Py());
    h_Ks_l_p_z->Fill(p4_Ks_star.Pz(),p4_l_star.Pz());


    h_delta_phi_Ks_l->Fill(delta_phi_Ks_l);
    h_delta_theta_Ks_l->Fill(delta_theta_Ks_l);
    h_delta_eta_Ks_l->Fill(delta_eta_Ks_l);
    h_delta_R_Ks_l->Fill(delta_R_Ks_l);

    h_pt_S_delta_phi_Ks_l->Fill(pt_S,delta_phi_Ks_l);
    h_p_S_delta_phi_Ks_l->Fill(p3_S.Mag(),delta_phi_Ks_l);
    h_pt_S_delta_theta_Ks_l->Fill(pt_S,delta_theta_Ks_l);
    h_p_S_delta_theta_Ks_l->Fill(p3_S.Mag(),delta_theta_Ks_l);
    h_delta_phi_delta_theta_Ks_l->Fill(delta_phi_Ks_l,delta_theta_Ks_l);
    h_delta_phi_delta_eta_Ks_l->Fill(delta_phi_Ks_l,delta_eta_Ks_l);

    h_S_pt_after_cuts->Fill(pt_S);
    h_S_p_after_cuts->Fill(p3_S.Mag());
    h_S_eta_after_cuts->Fill(eta_S);
    h_Ks_pt_after_cuts->Fill(p4_Ks_star.Pt());
    h_Ks_eta_after_cuts->Fill(p4_Ks_star.Eta());
    h_l_pt_after_cuts->Fill(p4_l_star.Pt());
    h_l_eta_after_cuts->Fill(p4_l_star.Eta());
 }

 void fillArmenterosHist(Double_t p_l_trans, Double_t p_Ks_trans, Double_t p_l_long, Double_t p_Ks_long, Double_t alpha_ArmPod){


    h_p_l_trans -> Fill(p_l_trans);
    h_p_Ks_trans -> Fill(p_Ks_trans);
    h_p_l_long -> Fill(p_l_long);
    h_p_Ks_long -> Fill(p_Ks_long);
    h2_p_Ks_l_trans ->Fill(p_Ks_trans, p_l_trans);
    h2_p_Ks_l_long ->Fill(p_Ks_long, p_l_long);

    h2_ArmPod_Lambda->Fill(alpha_ArmPod,p_l_trans); //make Armenteros-Podolanski Plot for Lambda
    h2_ArmPod_Ks->Fill(alpha_ArmPod,p_Ks_trans); // make Armenteros-Podolanski Plot for Ks
    h2_ArmPod->Fill(alpha_ArmPod,p_l_trans); //combine Lambda and Kshort
    h2_ArmPod->Fill(alpha_ArmPod,p_Ks_trans); //combine Lambda and Kshort

}


int simulation(string outputDir, string rootFileName, int nIterations, string polarization)
{
    //some hitso's that need to know about "particle"


    TRandom *random = new TRandom3();
    
    TFile *top = new TFile(rootFileName.c_str(),"RECREATE");
    
    TDirectory *dir_settings = top->mkdir("settings");
    TDirectory *dir_StartParticle = top->mkdir("StartParticle");
    TDirectory *dir_StarFrame = top->mkdir("StarFrame");
    TDirectory *dir_DecayProductsLabBeforeCuts = top->mkdir("DecayProductsLabBeforeCuts");
    TDirectory *dir_DecayProductsLabAfterCuts = top->mkdir("DecayProductsLabAfterCuts");
    TDirectory *dir_Arementeros = top->mkdir("dir_Arementeros");

    //getting proper mass of start particle (Xi1820 or S mass) and the neutron mass (put to 0 when there is no collission)
    Double_t M_Start = return_M_Start(particle);
    Double_t m_n = return_M_n(particle);
    //loading TSalis fucntion
    TF1 *fTSalisSpt = TSalis(M_Start);
    //loading eta distributions (one for positive eta and one for negative eta)
    TF1 *fEtaDistrHadr = EtaDistrHadr();
    //getting reco efficiency of the Ks and Lambda in function of pt    
    TEfficiency *TEff_Ks_pt = getRecoEfficiencyKs_L0("Ks_recoeff_gen_Pt");    
    TEfficiency *TEff_L0_pt = getRecoEfficiencyKs_L0("L0_recoeff_gen_Pt");    

    int i = 0; 
    int numerator_eff = 0; 
    bool verbose = false;  
    
    //start the main loop where you generate particles
    while(i<nIterations){

	//******************start conditions*****************************************//

    //pt distribution of the S particle:  https://arxiv.org/pdf/1401.4835.pdf
    Double_t pt_S = fTSalisSpt->GetRandom(); 

	//from https://arxiv.org/pdf/1002.0621.pdf, page 14 the eta distribution is relatively flat
    Double_t eta_S = random->Uniform(-2.5,2.5);
    //Double_t eta_S = fEtaDistrHadr->GetRandom(0.,10.);
    //if(random->Uniform(0,1)<0.5) eta_S = -eta_S;
    

    //phi of the S, should be a uniform distribution in 0 to 2*pi
	Double_t phi_S = random->Uniform(2.*TMath::Pi());	

    //theta of the Ks, distributed according to a cosine for a isotropic distribution or differently for polarised particles
    Double_t argument_cos = return_argument_Acos_theta_Ks_angle(random, polarization);
    Double_t theta_Ks_star = TMath::ACos(argument_cos);
    
    //phi of the Ks, should be uniform distribution from 0 to 2*pi
	Double_t phi_Ks_star = random->Uniform(2.*TMath::Pi());	

    fillHistosSettings(argument_cos, theta_Ks_star);
   
	//******************end start conditions*****************************************//

	//******************calculate the S particle parameters******************************************//
	
    Double_t theta_S = 2*TMath::ATan(TMath::Exp(-eta_S));

	Double_t p_S = pt_S/TMath::Sin(theta_S);

	TVector3 p3_S(0.,0.,0.); 

	p3_S.SetMagThetaPhi(p_S,theta_S,phi_S); 

    
    fillHistosSParticle(pt_S, p_S, theta_S, phi_S, p3_S, eta_S);

	//******************end calculate the S particle parameters******************************************//
		
	//the S hits on a neutron and transfers its momentum to the S+n system.
	
	//******************calculate the S+n particle parameters******************************************//
	//assume the S momentum is transfered to the S+n system, then the momenta are the same
	TVector3 p3_S_n =  p3_S;	
  
    Double_t E_S = pow(p3_S*p3_S+M_Start*M_Start,0.5);
  
    Double_t E_S_n = E_S + m_n;
  
    Double_t M_Start_n = pow(E_S_n*E_S_n-p3_S_n*p3_S_n,0.5);

    fillHistosS_n(E_S_n, M_Start_n);
	//******************end calculate the S+n particle parameters******************************************//
	
	//The S + n then decays to a Ks and Lambda. 
	//This is a 2 body decay (http://www.helsinki.fi/~www_sefo/phenomenology/Schlippe_relativistic_kinematics.pdf)
	
	//******************calculate the Ks and Lambda parameter in the S+n rest frame****************************************//
	//E_Ks_star and p_Ks_star energy and momentum of the Ks in the rest frame of the S + n
	Double_t E_Ks_star = (M_Start_n*M_Start_n+m_Ks*m_Ks-m_l*m_l)/(2*M_Start_n);
	
    Double_t p_Ks_star = pow(E_Ks_star*E_Ks_star - m_Ks*m_Ks,0.5);
    
    if(verbose) cout << "2nd way of calculating p_Ks_star:  " << pow((M_Start_n*M_Start_n-(m_Ks-m_l)*(m_Ks-m_l))*(M_Start_n*M_Start_n-(m_Ks+m_l)*(m_Ks+m_l)),0.5)/(2.*M_Start_n) << endl;

	Double_t m_Ks_star = pow(E_Ks_star*E_Ks_star-p_Ks_star*p_Ks_star,0.5);

	//E_l_star and p_l_star energy and momentum of the l in the rest frame of the S + n
	Double_t E_l_star = (M_Start_n*M_Start_n+m_l*m_l-m_Ks*m_Ks)/(2*M_Start_n);

	Double_t p_l_star = pow(E_l_star*E_l_star - m_l*m_l,0.5);

    Double_t m_l_star = pow(E_l_star*E_l_star-p_l_star*p_l_star,0.5);
	//p4_Ks_star 4 momentum vector of the Ks 
	TLorentzVector p4_Ks_star(0.,0.,0.,0.);

	p4_Ks_star.SetPtEtaPhiE(p_Ks_star*TMath::Sin(theta_Ks_star),-TMath::Log(TMath::Tan(theta_Ks_star/2.)),phi_Ks_star,E_Ks_star);	

	TLorentzVector p4_l_star(0.,0.,0.,0.);

	p4_l_star.SetPxPyPzE(-p4_Ks_star.Px(),-p4_Ks_star.Py(),-p4_Ks_star.Pz(),E_l_star);

    fillHistosStar( E_Ks_star,  p_Ks_star,  m_Ks_star,  E_l_star,  p_l_star,  m_l_star,  p4_Ks_star,  p4_l_star,  argument_cos);


	//******************end calculate the Ks and Lambda parameter in the S+n rest frame****************************************//

	//you have now the vector of the Ks and the Lambda in the rest frame of the S+n. Now will boost these vectors to the reference
	//frame where the S+n is actually moving with momentum p3_S, this is the detector reference frame
	
	
	//******************boost to the reference frame of the detector****************************************//
	//boosting the p4_Ks_star and p4_l_star allong the S+n system momentum	
	p4_Ks_star.Boost(p3_S_n.x()/E_S_n,p3_S_n.y()/E_S_n,p3_S_n.z()/E_S_n); //beta = p/E
	p4_l_star.Boost(p3_S_n.x()/E_S_n,p3_S_n.y()/E_S_n,p3_S_n.z()/E_S_n); //beta = p/E

    fillHistosAfterBoost(pt_S, eta_S, p4_Ks_star, p4_l_star);

    //******************end boost to the reference frame of the detector****************************************//

    //******************calculate efficiencies (pt dependence)**************************************************//

    Double_t binSize = 0.2;
    Double_t binOffset = 0.1;

    int binKsEff = TEff_Ks_pt->FindFixBin(p4_Ks_star.Pt());
    Double_t intEffKs = TEff_Ks_pt->GetEfficiency(binKsEff);
    Double_t lower_error_intEffKs = TEff_Ks_pt->GetEfficiencyErrorLow(binKsEff);

    h2_check_intEff_Ks->Fill(p4_Ks_star.Pt(),intEffKs); 


    int binL0Eff = TEff_L0_pt->FindFixBin(p4_l_star.Pt());
    Double_t intEffL0 = TEff_L0_pt->GetEfficiency(binL0Eff);
    Double_t lower_error_intEffL0 = TEff_L0_pt->GetEfficiencyErrorLow(binL0Eff);
    
    h2_check_intEff_L0->Fill(p4_l_star.Pt(),intEffL0); 

    //*****************end calculate efficiencies**************************************************************//    

    //*****************doing cuts******************************************************************************//

    //efficiency of the reconstruction of a Ks and a Lambda
    Double_t L0_reco_eff = random->Uniform(0.,1.);
    Double_t Ks_reco_eff = random->Uniform(0.,1.);

    bool KsReconstructed = false;
    if(p4_Ks_star.Pt()< 3.){
        if(Ks_reco_eff <= intEffKs) KsReconstructed = true;
    }
    else{ 
        if(Ks_reco_eff <= 0.4) KsReconstructed = true;
    }
    bool Ks_survives_cuts =  (p4_Ks_star.Pt() > min_pt_Ks) && (abs(p4_Ks_star.Eta()) <= max_eta)  && KsReconstructed;


    hEff_Ks->Fill(Ks_survives_cuts, pt_S);
    
    bool L0Reconstructed = false;
    if(p4_l_star.Pt()< 4.){
        if(L0_reco_eff <= intEffL0) L0Reconstructed = true;
    }
    else{ 
        if(L0_reco_eff <= 0.55) L0Reconstructed = true;
    }
    bool l_survives_cuts = (p4_l_star.Pt() > min_pt_l) && (abs(p4_l_star.Eta()) <= max_eta) &&  L0Reconstructed;


    hEff_L0->Fill(l_survives_cuts, pt_S);

    hEff_L0_Ks->Fill(l_survives_cuts && Ks_survives_cuts, pt_S);


    if(!Ks_survives_cuts) continue;
    if(!l_survives_cuts) continue;

    Double_t delta_phi_Ks_l = p4_Ks_star.DeltaPhi(p4_l_star);
    Double_t delta_theta_Ks_l =p4_Ks_star.Theta()-p4_l_star.Theta();
    Double_t delta_eta_Ks_l =p4_Ks_star.Eta()-p4_l_star.Eta();
    Double_t delta_R_Ks_l = p4_Ks_star.DeltaR(p4_l_star);

    fillHistosAfterCuts( p4_Ks_star,  p4_l_star,  pt_S,  p3_S,  eta_S, delta_phi_Ks_l, delta_theta_Ks_l, delta_eta_Ks_l, delta_R_Ks_l);
  
    //*****************end doing cuts**************************************************************************//



    //*******************************making Armenteros-Podolanski Plot**************************************************************//
    //http://www.star.bnl.gov/~gorbunov/main/node48.html
   // TVector3 p3_S_norm(0.,0.,0.); 
   // p3_S_norm.SetXYZ(p3_S.Px()/p3_S.Mag(),p3_S.Py()/p3_S.Mag(),p3_S.Pz()/p3_S.Mag());
    TVector3 p3_l(0.,0.,0.);
    p3_l.SetXYZ(p4_l_star.Px(),p4_l_star.Py(),p4_l_star.Pz());

    TVector3 p3_Ks(0.,0.,0.);
    p3_Ks.SetXYZ(p4_Ks_star.Px(),p4_Ks_star.Py(),p4_Ks_star.Pz());
    
    Double_t angle_S_l  = p3_S.Angle(p3_l); //angle between S and Lambda in the lab frame
    Double_t angle_S_Ks  = p3_S.Angle(p3_Ks); //angle between S and Kshort in the lab frame
    Double_t p_l_long = p3_l.Mag()*TMath::Cos(angle_S_l); //longitudnal momentum of Lambda along S direction in the lab frame
    Double_t p_Ks_long = p3_Ks.Mag()*TMath::Cos(angle_S_Ks); //longitudnal momentum of Ks along S direction in the lab frame
    Double_t p_l_trans = p3_l.Mag()*TMath::Sin(angle_S_l); //transversal momentum of Lambda along S direction in the lab frame
    Double_t p_Ks_trans = p3_Ks.Mag()*TMath::Sin(angle_S_Ks); //transversal momentum of Ks along S direction in the lab frame

    Double_t alpha_ArmPod = (p_Ks_long-p_l_long)/(p_Ks_long+p_l_long); //alpha parameter in the Armenteros-Podolanski Plot take Ks as pos and Lambda as neg (convention)

    fillArmenterosHist(p_l_trans, p_Ks_trans, p_l_long, p_Ks_long, alpha_ArmPod);

	//******************check the above calculation by calculating the invariant mass of the S****************************************//
    TLorentzVector p4_sum_Ks_l = p4_Ks_star+p4_l_star;

    Double_t M_Start_n_check = p4_sum_Ks_l.M();

    Double_t approx_Inv_Mass =  M_Start + m_n + m_n*p3_S*p3_S/(4*M_Start*(M_Start+m_n)) ;
    //cout << "Inv mass Ks and Lambda 2: " << approx_Inv_Mass << endl;
    h_M_Start_n_check->Fill(M_Start_n_check);
    h_M_Start_n_Inv_Mass_calc->Fill(approx_Inv_Mass);
    h_M_Start_n_Inv_Mass_calc->SetFillColor(42);

    //******************end check the above calculation by calculating the invariant mass of the S****************************************//

	if(abs(delta_phi_Ks_l)>1 && abs(delta_phi_Ks_l)<2.5) numerator_eff++;

	i++;

    }//end for loop

 cout << "the efficiency on the signal as a result of the delta_phi cut on background: " << (double)numerator_eff/(double)i << endl; 

 //*****************************CALCULATE THE RECONSTRUCTION EFFICIENCY FOR ONE S*******************************************
 Double_t numerator = 0;
 Double_t denumenator = 0;
 for(int pt_bin = 0; pt_bin <= 50; pt_bin++){
    Double_t n_S_pt = h_S_pt->GetBinContent(pt_bin);
    Double_t effL0 = hEff_L0->GetEfficiency(pt_bin);
    Double_t effKs = hEff_Ks->GetEfficiency(pt_bin);

    numerator = numerator + n_S_pt*effL0*effKs;
    denumenator = denumenator + n_S_pt;
 }

 cout << "The efficiency we have all been waiting for: " << numerator/denumenator << endl;
 //*************************************************************************************************************************



 TCanvas * c1 = new TCanvas("c", "c");
 
 dir_settings->cd();

 TEff_Ks_pt->Write();

 h_cos_theta_Ks->Write();
 h2_cos_theta_Ks->Write();
 h2_theta_Ks->Write();
 h_J32_to_M12_theta->Write();
 h_sqrt->Write();
 h2_sqrt->Write();

 h_theta_Ks->Write();
 h_u->Write();

 



 dir_StartParticle->cd(); 

 h_S_pt->Write();
 h_S_pt->Draw();
 c1->SaveAs((outputDir+"/h_"+particle+"_pt.pdf").c_str(),"pdf");

 h_S_p->Write();
 h_theta_S->Write();
 h_eta_S->Write();
    
 h_S_p_x->Write();
 h_S_p_x->Draw();
 c1->SaveAs((outputDir+"/h_"+particle+"_p_x.pdf").c_str(),"pdf");
 h2_S_p_x->Write();
 h_S_p_y->Write();
 h_S_p_y->Draw();
 c1->SaveAs((outputDir+"/h_"+particle+"_p_y.pdf").c_str(),"pdf");
 h2_S_p_y->Write();
 h_S_p_z->Write();
 h_S_p_z->Draw();
 c1->SaveAs((outputDir+"/h_"+particle+"_p_z.pdf").c_str(),"pdf");
 h2_S_p_z->Write();
 h2_eta_S_p_z->Write();
 h_S_P->Write(); 

 h_S_p_xy->Write();
 h_S_p_xy->Draw("colz");
 c1->SaveAs((outputDir+"/h_"+particle+"_p_xy.png").c_str(),"png");
 
 h_S_p_xz->Write();
 h_S_p_xz->Draw("colz");
 c1->SaveAs((outputDir+"/h_"+particle+"_p_xz.png").c_str(),"png");
 
 h_S_p_yz->Write();
 h_S_p_yz->Draw("colz");
 c1->SaveAs((outputDir+"/h_"+particle+"_p_yz.png").c_str(),"png");

 h_E_S_n->Write();
 h_E_S_n->Draw();
 c1->SaveAs((outputDir+"/h_E_"+particle+"_n.pdf").c_str(),"pdf");
 h_M_S_n->Write();
 h_M_S_n->Draw();
 c1->SaveAs((outputDir+"/h_M_"+particle+"_n.pdf").c_str(),"pdf");
 





 dir_StarFrame->cd();    

 h_E_Ks_star->Write();
 h_E_Ks_star->Draw();
 c1->SaveAs((outputDir+"/h_E_Ks_star.pdf").c_str(),"pdf");
 
 h_m_Ks_star->Write();
 h_m_Ks_star->Draw();
 c1->SaveAs((outputDir+"/h_m_Ks_star.pdf").c_str(),"pdf");
 
 h_E_l_star->Write();
 h_E_l_star->Draw();
 c1->SaveAs((outputDir+"/h_E_l_star.pdf").c_str(),"pdf");
 
 h_m_l_star->Write();
 h_m_l_star->Draw();
 c1->SaveAs((outputDir+"/h_m_l_star.pdf").c_str(),"pdf");
 
 h_p_Ks_star->Write();
 h_p_Ks_star->Draw();
 c1->SaveAs((outputDir+"/h_p_Ks_star.pdf").c_str(),"pdf");
 
 h_p_l_star->Write();
 h_p_l_star->Draw();
 c1->SaveAs((outputDir+"/h_p_l_star.pdf").c_str(),"pdf");
 
 h_Ks_star_p_x->Write();
 h_Ks_star_p_x->Draw();
 c1->SaveAs((outputDir+"/h_Ks_star_p_x.pdf").c_str(),"pdf");
 
 h_Ks_star_p_y->Write();
 h_Ks_star_p_y->Draw();
 c1->SaveAs((outputDir+"/h_Ks_star_p_y.pdf").c_str(),"pdf");
 
 h_Ks_star_p_z->Write();
 h_Ks_star_p_z->Draw();
 c1->SaveAs((outputDir+"/h_Ks_star_p_z.pdf").c_str(),"pdf");
 
 h_Ks_star_p->Write();
 h_Ks_star_p->Draw();
 c1->SaveAs((outputDir+"/h_Ks_star_p.pdf").c_str(),"pdf");
 

 h_l_star_p_x->Write();
 h_l_star_p_x->Draw();
 c1->SaveAs((outputDir+"/h_l_star_p_x.pdf").c_str(),"pdf");
 
 h_l_star_p_y->Write();
 h_l_star_p_y->Draw();
 c1->SaveAs((outputDir+"/h_l_star_p_y.pdf").c_str(),"pdf");
 
 h_l_star_p_z->Write();
 h_l_star_p_z->Draw();
 c1->SaveAs((outputDir+"/h_l_star_p_z.pdf").c_str(),"pdf");
  
 h_l_star_p->Write();
 h_l_star_p->Draw();
 c1->SaveAs((outputDir+"/h_l_star_p.pdf").c_str(),"pdf");
 
 h_Ks_star_p_xy->Write();
 h_Ks_star_p_xy->Draw("colz");
 c1->SaveAs((outputDir+"/h_Ks_star_p_xy.png").c_str(),"png");
 
 h_Ks_star_p_xz->Write();
 h_Ks_star_p_xz->Draw("colz");
 c1->SaveAs((outputDir+"/h_Ks_star_p_xz.png").c_str(),"png");
 
 h_Ks_star_p_yz->Write();
 h_Ks_star_p_yz->Draw("colz");
 c1->SaveAs((outputDir+"/h_Ks_star_p_yz.png").c_str(),"png");
 cout << "writing projection in X " << endl;
 h_Ks_star_p_xy_proj_x = h_Ks_star_p_xy->ProjectionX();
 h_Ks_star_p_xy_proj_x->Write();

 h_Ks_star_p_xy_norm->Write();
 h_Ks_star_p_xy_norm->Draw("colz");
 c1->SaveAs((outputDir+"/h_Ks_star_p_xy_norm.png").c_str(),"png");
 
 h_Ks_star_p_xz_norm->Write();
 h_Ks_star_p_xz_norm->Draw("colz");
 c1->SaveAs((outputDir+"/h_Ks_star_p_xz_norm.png").c_str(),"png");
 
 h_Ks_star_p_yz_norm->Write();
 h_Ks_star_p_yz_norm->Draw("colz");
 c1->SaveAs((outputDir+"/h_Ks_star_p_yz_norm.png").c_str(),"png");

  h_l_star_p_xy->Write();
 h_l_star_p_xz->Write();
 h_l_star_p_yz->Write();

 h_Ks_l_star_p_x->Write();
 h_Ks_l_star_p_y->Write();
 h_Ks_l_star_p_z->Write();

h_delta_phi_Ks_l_star->Write();
 h_delta_phi_Ks_l_star->Draw();
 c1->SaveAs((outputDir+"/h_delta_phi_Ks_l_star.pdf").c_str(),"pdf");

 h_delta_theta_Ks_l_star->Write();
 h_delta_theta_Ks_l_star->Draw();
 c1->SaveAs((outputDir+"/h_delta_theta_Ks_l_star.pdf").c_str(),"pdf");
 
 h2_delta_theta_Ks_l_star->Write();
 h2_delta_theta_Ks_l_star->Draw();
 c1->SaveAs((outputDir+"/h2_delta_theta_Ks_l_star.pdf").c_str(),"pdf");
 
 h2_delta_eta_Ks_l_star->Write();
 h2_delta_eta_Ks_l_star->Draw();
 c1->SaveAs((outputDir+"/h2_delta_eta_Ks_l_star.pdf").c_str(),"pdf");
 
  h_sum_theta_Ks_l_star->Write();





 dir_DecayProductsLabBeforeCuts->cd();    

  h_S_pt_before_cuts->Write();
  h_S_eta_before_cuts->Write();
  h_Ks_pt_before_cuts->Write();
  h_Ks_eta_before_cuts->Write();
  h_l_pt_before_cuts->Write();
  h_l_eta_before_cuts->Write();





dir_DecayProductsLabAfterCuts->cd();

  h_S_pt_after_cuts->Write();
  h_S_p_after_cuts->Write();
  h_S_eta_after_cuts->Write();
  h_Ks_pt_after_cuts->Write();
  h_Ks_eta_after_cuts->Write();
  h_l_pt_after_cuts->Write();
  h_l_eta_after_cuts->Write();

 h_Ks_l_p_x->Write();
 h_Ks_l_p_y->Write();
 h_Ks_l_p_z->Write();

 h_Ks_p_x->Write();
 h_Ks_p_x->Draw();
 c1->SaveAs((outputDir+"/h_Ks_p_x.pdf").c_str(),"pdf");
 
 h_Ks_p_y->Write();
 h_Ks_p_y->Draw();
 c1->SaveAs((outputDir+"/h_Ks_p_y.pdf").c_str(),"pdf");
 
 h_Ks_p_z->Write();
 h_Ks_p_z->Draw();
 c1->SaveAs((outputDir+"/h_Ks_p_z.pdf").c_str(),"pdf");
 
 h_Ks_pt->Write();
 h_Ks_pt->Draw();
 c1->SaveAs((outputDir+"/h_Ks_pt.pdf").c_str(),"pdf");
 
 h_l_p_x->Write();
 h_l_p_x->Draw();
 c1->SaveAs((outputDir+"/h_l_p_x.pdf").c_str(),"pdf");
 
 h_l_p_y->Write();
 h_l_p_y->Draw();
 c1->SaveAs((outputDir+"/h_l_p_y.pdf").c_str(),"pdf");
 
 h_l_p_z->Write();
 h_l_p_z->Draw();
 c1->SaveAs((outputDir+"/h_l_p_z.pdf").c_str(),"pdf");
 
 h_l_pt->Write();
 h_l_pt->Draw();
 c1->SaveAs((outputDir+"/h_l_pt.pdf").c_str(),"pdf");

 Double_t norm1 = h_delta_phi_Ks_l->GetEntries();
 //if(norm1) h_delta_phi_Ks_l->Scale(1/norm1);
 h_delta_phi_Ks_l->Write();
 h_delta_phi_Ks_l->Draw();
 c1->SaveAs((outputDir+"/h_delta_phi_Ks_l.pdf").c_str(),"pdf");

 Double_t norm2 = h_delta_theta_Ks_l->GetEntries();
 //if(norm2) h_delta_theta_Ks_l->Scale(1/norm2); 
 h_delta_theta_Ks_l->Write();
 h_delta_theta_Ks_l->Draw();
 c1->SaveAs((outputDir+"/h_delta_theta_Ks_l.pdf").c_str(),"pdf");
 
 Double_t norm3 = h_delta_eta_Ks_l->GetEntries();
 //if(norm3) h_delta_eta_Ks_l->Scale(1/norm3); 
 h_delta_eta_Ks_l->Write();
 h_delta_eta_Ks_l->Draw();
 c1->SaveAs((outputDir+"/h_delta_eta_Ks_l.pdf").c_str(),"pdf");
 
 Double_t norm4 = h_delta_R_Ks_l->GetEntries();
 //if(norm4) h_delta_R_Ks_l->Scale(1/norm4); 
 h_delta_R_Ks_l->Write();
 h_delta_R_Ks_l->Draw();
 c1->SaveAs((outputDir+"/h_delta_R_Ks_l.pdf").c_str(),"pdf");
 
 h_pt_S_delta_phi_Ks_l->Write();
 h_pt_S_delta_phi_Ks_l->Draw("colz");
 c1->SaveAs((outputDir+"/h_pt_"+particle+"_delta_phi_Ks_l.png").c_str(),"png");
 
 h_p_S_delta_phi_Ks_l->Write();
 h_p_S_delta_phi_Ks_l->Draw("colz");
 c1->SaveAs((outputDir+"/h_p_"+particle+"_delta_phi_Ks_l.png").c_str(),"png");
 
 h_pt_S_delta_theta_Ks_l->Write();
 h_pt_S_delta_theta_Ks_l->Draw("colz");
 c1->SaveAs((outputDir+"/h_pt_"+particle+"_delta_theta_Ks_l.png").c_str(),"png");
 
 h_p_S_delta_theta_Ks_l->Write();
 h_p_S_delta_theta_Ks_l->Draw("colz");
 c1->SaveAs((outputDir+"/h_p_"+particle+"_delta_theta_Ks_l.png").c_str(),"png");
 
 h_delta_phi_delta_theta_Ks_l->Write();
 h_delta_phi_delta_theta_Ks_l->Draw("colz");
 c1->SaveAs((outputDir+"/h_delta_phi_delta_theta_Ks_l.png").c_str(),"png");

 h_delta_phi_delta_eta_Ks_l->Write();
 h_delta_phi_delta_eta_Ks_l->Draw("colz");
 c1->SaveAs((outputDir+"/h_delta_phi_delta_eta_Ks_l.png").c_str(),"png");

 h_M_Start_n_check->Write();
 h_M_Start_n_check->Draw();
 c1->SaveAs((outputDir+"/h_M_"+particle+"_n_check.pdf").c_str(),"pdf");
 
 h_M_Start_n_Inv_Mass_calc->Write();
 h_M_Start_n_Inv_Mass_calc->Draw("same");
 c1->SaveAs((outputDir+"/h_M_"+particle+"_n_Inv_Mass_calc.pdf").c_str(),"pdf");
 
 h2_check_intEff_Ks->Write(); 
 h2_check_intEff_L0->Write();  

 hEff_Ks->Write();
 hEff_L0->Write();
 hEff_L0_Ks->Write();





dir_Arementeros->cd();

 h2_ArmPod->Write();
 h2_ArmPod->Draw("colz");
 c1->SaveAs((outputDir+"/h2_ArmPod.png").c_str(),"png"); 

 h2_ArmPod_Ks->Write();
 h2_ArmPod_Ks->Draw("colz");
 c1->SaveAs((outputDir+"/h2_ArmPod_Ks.png").c_str(),"png"); 

 h2_ArmPod_Lambda->Write();
 h2_ArmPod_Lambda->Draw("colz");
 c1->SaveAs((outputDir+"/h2_ArmPod_Lambda.png").c_str(),"png"); 

 EtaDistrHadr()->Draw();
 
  h_p_l_trans->Write();
 h_p_Ks_trans->Write();
 h_p_l_long->Write();
 h_p_Ks_long->Write();

 h2_p_Ks_l_trans ->Write();
 h2_p_Ks_l_long ->Write(); 



 top->Write();
 return 0;
}

