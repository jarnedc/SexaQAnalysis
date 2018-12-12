#include "../interface/Analyzer_SIM_Sexaq.h"
#include <typeinfo>


Analyzer_SIM_Sexaq::Analyzer_SIM_Sexaq(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  m_genParticlesTag_GEN(pset.getParameter<edm::InputTag>("genCollection_GEN")),
  m_genParticlesTag_SIM_GEANT(pset.getParameter<edm::InputTag>("genCollection_SIM_GEANT")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),

  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN)),
  m_genParticlesToken_SIM_GEANT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_SIM_GEANT)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag))

{

}


void Analyzer_SIM_Sexaq::beginJob() {
  

    //for the GEN particles
    TFileDirectory dir_genParticles = m_fs->mkdir("genParticles");
    histos_th1f[b+"h_genParticles_pdgid"]= dir_genParticles.make<TH1F>(b+"h_genParticles_pdgid", b+"h_genParticles_pdgid; pdgId of Gen particles; pdgId ",6000,-3000,3000);


    TFileDirectory dir_genParticles_antiS = dir_genParticles.mkdir("genParticles_antiS");
    //for the GEN particles: only the S
    histos_th1f[b+"h_genParticles_antiS_pdgid"]= dir_genParticles_antiS.make<TH1F>(b+"h_genParticles_antiS_pdgid", b+"h_genParticles_antiS_pdgid; pdgId of anti S Gen particles; pdgId ",10,-1020000020-5,-1020000020+5);
    histos_th1f[b+"h_genParticles_antiS_pt"]= dir_genParticles_antiS.make<TH1F>(b+"h_genParticles_antiS_pt", b+"h_genParticles_antiS_pt; pt of anti S Gen particles; pt (GeV) ",200,0,20);
    histos_th1f[b+"h_genParticles_antiS_n_daughters"]= dir_genParticles_antiS.make<TH1F>(b+"h_genParticles_antiS_n_daughters", b+"h_genParticles_antiS_n_daughters; #daughters of anti S Gen particles ",20,0,20);
  
    //for the SIM-GEANT particles 
    //S
    TFileDirectory dir_simgeantParticles = m_fs->mkdir("simgeantParticles");
    TFileDirectory dir_simgeantParticles_antiS = dir_simgeantParticles.mkdir("antiS");
    histos_th1f[b+"h_simgeantParticles_antiS_pdgid"]= dir_simgeantParticles_antiS.make<TH1F>(b+"h_simgeantParticles_antiS_pdgid", b+"h_simgeantParticles_antiS_pdgid; pdgid anti-S(GeV); ",10,-1020000020-5,-1020000020+5);
    histos_th1f[b+"h_simgeantParticles_antiS_pt"]= dir_simgeantParticles_antiS.make<TH1F>(b+"h_simgeantParticles_antiS_pt", b+"h_simgeantParticles_antiS_pt; pt anti-S(GeV); ", 200, 0,20);
    histos_th1f[b+"h_simgeantParticles_antiS_p"]= dir_simgeantParticles_antiS.make<TH1F>(b+"h_simgeantParticles_antiS_p", b+"h_simgeantParticles_antiS_p; p anti-S(GeV); ", 2000, 0, 200);
    histos_th2f[b+"h2_simgeantParticles_antiS_p_eta"]= dir_simgeantParticles_antiS.make<TH2F>(b+"h2_simgeantParticles_antiS_p_eta", b+"h2_simgeantParticles_antiS_p_eta; p anti-S (GeV); eta anti-S ", 2000, 0, 200, 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_mass"]= dir_simgeantParticles_antiS.make<TH1F>(b+"h_simgeantParticles_antiS_mass", b+"h_simgeantParticles_antiS_mass; mass anti-S (GeV); ", 1000, 0,10);
    histos_th1f[b+"h_simgeantParticles_antiS_mass_from_Ep"]= dir_simgeantParticles_antiS.make<TH1F>(b+"h_simgeantParticles_antiS_mass_from_Ep", b+"h_simgeantParticles_antiS_mass_from_Ep; mass anti-S (GeV); ", 1000, 0,10);
    histos_th1f[b+"h_simgeantParticles_antiS_eta"]= dir_simgeantParticles_antiS.make<TH1F>(b+"h_simgeantParticles_antiS_eta", b+"h_simgeantParticles_antiS_eta; #eta(anti-S); ", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_phi"]= dir_simgeantParticles_antiS.make<TH1F>(b+"h_simgeantParticles_antiS_phi", b+"h_simgeantParticles_antiS_phi; #Phi(anti-S)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_vxy"]= dir_simgeantParticles_antiS.make<TH1F>(b+"h_simgeantParticles_antiS_vxy", b+"h_simgeantParticles_antiS_vxy; vxy (anti-S) (cm)", 1000, 0, 100);
    //for the S with no daughters
    TFileDirectory dir_simgeantParticles_antiS_no_daughters = dir_simgeantParticles_antiS.mkdir("no_daughters");
    histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_pt"]= dir_simgeantParticles_antiS_no_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_no_daughters_pt", b+"h_simgeantParticles_antiS_no_daughters_pt; pt anti-S(GeV); ", 200, 0,20);
    histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_p"]= dir_simgeantParticles_antiS_no_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_no_daughters_p", b+"h_simgeantParticles_antiS_no_daughters_p; p anti-S(GeV); ", 2000, 0, 200);
    histos_th2f[b+"h2_simgeantParticles_antiS_no_daughters_p_eta"]= dir_simgeantParticles_antiS_no_daughters.make<TH2F>(b+"h2_simgeantParticles_antiS_no_daughters_p_eta", b+"h2_simgeantParticles_antiS_no_daughters_p_eta; p anti-S (GeV); eta anti-S ", 2000, 0, 200, 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_eta"]= dir_simgeantParticles_antiS_no_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_no_daughters_eta", b+"h_simgeantParticles_antiS_no_daughters_eta; #eta(anti-S); ", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_phi"]= dir_simgeantParticles_antiS_no_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_no_daughters_phi", b+"h_simgeantParticles_antiS_no_daughters_phi; #Phi(anti-S)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_vxy"]= dir_simgeantParticles_antiS_no_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_no_daughters_vxy", b+"h_simgeantParticles_antiS_no_daughters_vxy; vxy (anti-S) (cm)", 1000, 0, 100);
    //for the S with 1 daughter
    TFileDirectory dir_simgeantParticles_antiS_1_daughter = dir_simgeantParticles_antiS.mkdir("1_daughter");
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_pt"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_pt", b+"h_simgeantParticles_antiS_1_daughter_pt; pt anti-S(GeV); ", 200, 0,20);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_p"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_p", b+"h_simgeantParticles_antiS_1_daughter_p; p anti-S(GeV); ", 2000, 0, 200);
    histos_th2f[b+"h2_simgeantParticles_antiS_1_daughter_p_eta"]= dir_simgeantParticles_antiS_1_daughter.make<TH2F>(b+"h2_simgeantParticles_antiS_1_daughter_p_eta", b+"h2_simgeantParticles_antiS_1_daughter_p_eta; p anti-S (GeV); eta anti-S ", 2000, 0, 200, 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_eta"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_eta", b+"h_simgeantParticles_antiS_1_daughter_eta; #eta(anti-S); ", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_phi"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_phi", b+"h_simgeantParticles_antiS_1_daughter_phi; #Phi(anti-S)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_vxy"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_vxy", b+"h_simgeantParticles_antiS_1_daughter_vxy; vxy (anti-S) (cm)", 1000, 0, 100);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_pdgid_daughter"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_pdgid_daughter", b+"h_simgeantParticles_antiS_1_daughter_pdgid_daughter; pdgid daughter anti-S", 8000, -4000, 4000);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_pt_daughter"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_pt_daughter", b+"h_simgeantParticles_antiS_1_daughter_pt_daughter; pt daughter (GeV)", 100, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_p_daughter"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_p_daughter", b+"h_simgeantParticles_antiS_1_daughter_p_daughter; p daughter (GeV)", 200, 0, 20);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_eta_daughter"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_eta_daughter", b+"h_simgeantParticles_antiS_1_daughter_eta_daughter; eta daughter ", 100, -4, 4);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_phi_daughter"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_phi_daughter", b+"h_simgeantParticles_antiS_1_daughter_phi_daughter; phi daughter (rad)", 100, -4, 4);
    histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_vxy_daughter"]= dir_simgeantParticles_antiS_1_daughter.make<TH1F>(b+"h_simgeantParticles_antiS_1_daughter_vxy_daughter", b+"h_simgeantParticles_antiS_1_daughter_vxy_daughter; vxy daughter (cm)", 1000, 0, 100);
    histos_th2f[b+"h_simgeantParticles_antiS_1_daughter_vzvx_daughter"]= dir_simgeantParticles_antiS_1_daughter.make<TH2F>(b+"h_simgeantParticles_antiS_1_daughter_vzvx_daughter", b+"h_simgeantParticles_antiS_1_daughter_vzvx_daughter; pdgid daughter anti-S", 2000, -500, 500, 2000, -500, 500);
    //for the S with 2 daughters
    TFileDirectory dir_simgeantParticles_antiS_2_daughters = dir_simgeantParticles_antiS.mkdir("2_daughters");
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pt"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_pt", b+"h_simgeantParticles_antiS_2_daughters_pt; pt anti-S(GeV); ", 200, 0,20);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_p"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_p", b+"h_simgeantParticles_antiS_2_daughters_p; p anti-S(GeV); ", 2000, 0, 200);
    histos_th2f[b+"h2_simgeantParticles_antiS_2_daughters_p_eta"]= dir_simgeantParticles_antiS_2_daughters.make<TH2F>(b+"h2_simgeantParticles_antiS_2_daughters_p_eta", b+"h2_simgeantParticles_antiS_2_daughters_p_eta; p anti-S (GeV); eta anti-S ", 2000, 0, 200, 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_eta"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_eta", b+"h_simgeantParticles_antiS_2_daughters_eta; #eta(anti-S); ", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_phi"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_phi", b+"h_simgeantParticles_antiS_2_daughters_phi; #Phi(anti-S)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_vxy"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_vxy", b+"h_simgeantParticles_antiS_2_daughters_vxy; vxy (anti-S) (cm)", 1000, 0, 100);
	
    //for the daughters of the S
    TFileDirectory dir_simgeantParticles_antiS_daughters = dir_simgeantParticles_antiS.mkdir("daughters");
    histos_th1f[b+"h_simgeantParticles_antiS_n_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_n_daughters", b+"h_simgeantParticles_antiS_n_daughters; #daughters of anti S Gen particles; ",20,0,20);
    histos_th1f[b+"h_simgeantParticles_antiS_pdgid_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_pdgid_daughters", b+"h_simgeantParticles_antiS_pdgid_daughters; pdgid of anti-S daughters", 8000, -4000, 4000);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_delta_phi_daughters", b+"h_simgeantParticles_antiS_delta_phi_daughters; #Delta Phi(Ks,L) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters_S_pt_larger_1GeV"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_delta_phi_daughters_S_pt_larger_1GeV", b+"h_simgeantParticles_antiS_delta_phi_daughters_S_pt_larger_1GeV; #Delta Phi(Ks,L) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_eta_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_delta_eta_daughters", b+"h_simgeantParticles_antiS_delta_eta_daughters; #Delta eta(Ks,L) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_R_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_delta_R_daughters", b+"h_simgeantParticles_antiS_delta_R_daughters; #Delta R(Ks,L)", 100, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_inv_mass_antiS_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_inv_mass_antiS_daughters", b+"h_simgeantParticles_antiS_inv_mass_antiS_daughters; inv mass anti-S (GeV)", 2000, 1.7, 1.9);
    histos_th1f[b+"h_simgeantParticles_antiS_openings_angle_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_openings_angle_daughters", b+"h_simgeantParticles_antiS_openings_angle_daughters; openings angle", 200, -10, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_pt_antiS_openings_angle_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH2F>(b+"h_simgeantParticles_antiS_pt_antiS_openings_angle_daughters", b+"h_simgeantParticles_antiS_pt_antiS_openings_angle_daughters;pt(S)(GeV); openings angle", 100, 0, 20, 200, -10, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_p_antiS_openings_angle_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH2F>(b+"h_simgeantParticles_antiS_p_antiS_openings_angle_daughters", b+"h_simgeantParticles_antiS_p_antiS_openings_angle_daughters; p(S)(GeV);openings angle", 100, 0, 100, 200, -10, 10);
    //with pt cuts on the Ks (pt > 0.9GeV) and the L(pt > 1.5GeV) 
   TFileDirectory dir_simgeantParticles_antiS_daughters_pt_cuts = dir_simgeantParticles_antiS.mkdir("pt_cuts");
   histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters_pt_cuts"]= dir_simgeantParticles_antiS_daughters_pt_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_delta_phi_daughters_pt_cuts", b+"h_simgeantParticles_antiS_delta_phi_daughters_pt_cuts; #Delta Phi(Ks,L) (rad)", 200, -10, 10); 
    histos_th1f[b+"h_simgeantParticles_antiS_delta_eta_daughters_pt_cuts"]= dir_simgeantParticles_antiS_daughters_pt_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_delta_eta_daughters_pt_cuts", b+"h_simgeantParticles_antiS_delta_eta_daughters_pt_cuts; #Delta eta(Ks,L) (rad)", 200, -10, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_delta_eta_delta_phi_daughters_pt_cuts"]= dir_simgeantParticles_antiS_daughters_pt_cuts.make<TH2F>(b+"h_simgeantParticles_antiS_delta_eta_delta_phi_daughters_pt_cuts", b+"h_simgeantParticles_antiS_delta_eta_delta_phi_daughters_pt_cuts; #Delta eta(Ks,L) (rad),  #Delta phi(Ks,L) (rad)", 200, -10, 10, 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_R_daughters_pt_cuts"]= dir_simgeantParticles_antiS_daughters_pt_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_delta_R_daughters_pt_cuts", b+"h_simgeantParticles_antiS_delta_R_daughters_pt_cuts; #Delta R(Ks,L)", 100, 0, 10);
    //with pt cuts on the Ks (pt > 0.9GeV) and the L(pt > 1.5GeV) and eta cuts on for the daughters of the Ks and Lambda to have abs(eta) < 2.5 
   TFileDirectory dir_simgeantParticles_antiS_daughters_pt_and_eta_and_delta_phi_cuts = dir_simgeantParticles_antiS_daughters_pt_cuts.mkdir("pt_and_eta_and_delta_phi_cuts");
   //antiS particle with  pt cuts on the Ks (pt > 0.9GeV) and the L(pt > 1.5GeV) and eta cuts on for the daughters of the Ks and Lambda to have abs(eta) < 2.5
   TFileDirectory dir_simgeantParticles_antiS_antiS_daughters_pt_and_eta_and_delta_phi_cuts = dir_simgeantParticles_antiS_daughters_pt_and_eta_and_delta_phi_cuts.mkdir("antiS");
   histos_th1f[b+"h_simgeantParticles_antiS_n_surv_background_cuts"]= dir_simgeantParticles_antiS_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_n_surv_background_cuts", b+"h_simgeantParticles_antiS_n_surv_background_cuts; # anti-S surv back cuts (1=survive)", 2, 0, 2); 
   histos_th1f[b+"h_simgeantParticles_antiS_p_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_antiS_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_p_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_p_pt_and_eta_and_delta_phi_cuts; p anti-S (GeV)", 200, 0, 10); 
   histos_th1f[b+"h_simgeantParticles_antiS_pt_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_antiS_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_pt_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_pt_pt_and_eta_and_delta_phi_cuts; pt anti-S (GeV)", 200, 0, 10); 
   histos_th1f[b+"h_simgeantParticles_antiS_phi_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_antiS_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_phi_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_phi_pt_and_eta_and_delta_phi_cuts; phi anti-S ", 200, -4, 4); 
   histos_th1f[b+"h_simgeantParticles_antiS_eta_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_antiS_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_eta_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_eta_pt_and_eta_and_delta_phi_cuts; eta anti-S ", 200, -4, 4); 
   //Ks particle with  pt cuts on the Ks (pt > 0.9GeV) and the L(pt > 1.5GeV) and eta cuts on for the daughters of the Ks and Lambda to have abs(eta) < 2.5
   TFileDirectory dir_simgeantParticles_Ks_antiS_daughters_pt_and_eta_and_delta_phi_cuts = dir_simgeantParticles_antiS_daughters_pt_and_eta_and_delta_phi_cuts.mkdir("Ks");
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_p_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_p_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_Ks_p_pt_and_eta_and_delta_phi_cuts; p Ks (GeV) ", 200, 0, 10); 
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_pt_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_pt_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_Ks_pt_pt_and_eta_and_delta_phi_cuts; pt Ks (GeV)", 200, 0, 10); 
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_phi_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_phi_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_Ks_phi_pt_and_eta_and_delta_phi_cuts; phi Ks ", 200, -4, 4); 
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_eta_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_eta_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_Ks_eta_pt_and_eta_and_delta_phi_cuts; eta Ks ", 200, -4, 4); 
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_vxy_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_vxy_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_Ks_vxy_pt_and_eta_and_delta_phi_cuts; vxy Ks (cm) ", 2000, 0, 200); 
   //L particle with  pt cuts on the Ks (pt > 0.9GeV) and the L(pt > 1.5GeV) and eta cuts on for the daughters of the Ks and Lambda to have abs(eta) < 2.5
   TFileDirectory dir_simgeantParticles_L_antiS_daughters_pt_and_eta_and_delta_phi_cuts = dir_simgeantParticles_antiS_daughters_pt_and_eta_and_delta_phi_cuts.mkdir("L");
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_p_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_L_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_p_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_L_p_pt_and_eta_and_delta_phi_cuts; p L (GeV) ", 200, 0, 10); 
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_pt_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_L_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_pt_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_L_pt_pt_and_eta_and_delta_phi_cuts; pt L (GeV)", 200, 0, 10); 
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_phi_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_L_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_phi_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_L_phi_pt_and_eta_and_delta_phi_cuts; phi L ", 200, -4, 4); 
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_eta_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_L_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_eta_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_L_eta_pt_and_eta_and_delta_phi_cuts; eta L ", 200, -4, 4); 
   histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_vxy_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_L_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_vxy_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_daughter_L_vxy_pt_and_eta_and_delta_phi_cuts; vxy L (cm) ", 2000, 0, 200); 
   //correlation Ks and L particle with  pt cuts on the Ks (pt > 0.9GeV) and the L(pt > 1.5GeV) and eta cuts on for the daughters of the Ks and Lambda to have abs(eta) < 2.5
   TFileDirectory dir_simgeantParticles_Ks_and_L_corr_antiS_daughters_pt_and_eta_and_delta_phi_cuts = dir_simgeantParticles_antiS_daughters_pt_and_eta_and_delta_phi_cuts.mkdir("Ks_and_L_corr");
   histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_delta_phi_daughters_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_delta_phi_daughters_pt_and_eta_and_delta_phi_cuts; #Delta Phi(Ks,L) (rad)", 200, -10, 10); 
    histos_th1f[b+"h_simgeantParticles_antiS_delta_eta_daughters_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_delta_eta_daughters_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_delta_eta_daughters_pt_and_eta_and_delta_phi_cuts; #Delta eta(Ks,L) (rad)", 200, -10, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_delta_eta_delta_phi_daughters_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH2F>(b+"h_simgeantParticles_antiS_delta_eta_delta_phi_daughters_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_delta_eta_delta_phi_daughters_pt_and_eta_and_delta_phi_cuts; #Delta eta(Ks,L) (rad),  #Delta phi(Ks,L) (rad)", 200, -10, 10, 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_R_daughters_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_delta_R_daughters_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_delta_R_daughters_pt_and_eta_and_delta_phi_cuts; #Delta R(Ks,L)", 100, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_openings_angle_daughters_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_openings_angle_daughters_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_openings_angle_daughters_pt_and_eta_and_delta_phi_cuts; openings angle(Ks,L)", 100, 0, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_pt_corr_daughters_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH2F>(b+"h_simgeantParticles_antiS_pt_corr_daughters_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_pt_corr_daughters_pt_and_eta_and_delta_phi_cuts; pt(Ks) (rad); pt(L) (rad)", 200, 0, 10, 200, 0, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_p_corr_daughters_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH2F>(b+"h_simgeantParticles_antiS_p_corr_daughters_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_p_corr_daughters_pt_and_eta_and_delta_phi_cuts; p(Ks) (rad); p(L) (rad)", 200, 0, 10, 200, 0, 10);

    //antiS particle with  pt cuts on the Ks (pt > 0.9GeV) and the L(pt > 1.5GeV) and eta cuts on for the daughters of the Ks and Lambda to have abs(eta) < 2.5
    TFileDirectory dir_simgeantParticles_Ks_and_L_corr_Ks_COM_antiS_daughters_pt_and_eta_and_delta_phi_cuts = dir_simgeantParticles_antiS_daughters_pt_and_eta_and_delta_phi_cuts.mkdir("Ks_and_L_corr_COM");
    histos_th2f[b+"h_simgeantParticles_antiS_px_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_Ks_COM_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH2F>(b+"h_simgeantParticles_antiS_px_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_px_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts; px Ks (GeV); px L (GeV)", 2000, -100, 100, 2000, -100, 100);
    histos_th2f[b+"h_simgeantParticles_antiS_py_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_Ks_COM_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH2F>(b+"h_simgeantParticles_antiS_py_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_py_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts; py Ks (GeV); py L (GeV)", 2000, -100, 100, 2000, -100, 100);
    histos_th2f[b+"h_simgeantParticles_antiS_pz_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_Ks_COM_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH2F>(b+"h_simgeantParticles_antiS_pz_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_pz_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts; pz Ks (GeV); pz L (GeV)", 2000, -100, 100, 2000, -100, 100);
    histos_th2f[b+"h_simgeantParticles_antiS_p_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_Ks_COM_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH2F>(b+"h_simgeantParticles_antiS_p_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_p_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts; p Ks (GeV); p L (GeV)", 2000, -100, 100, 2000, -100, 100);
    histos_th2f[b+"h_simgeantParticles_antiS_pt_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]= dir_simgeantParticles_Ks_and_L_corr_Ks_COM_antiS_daughters_pt_and_eta_and_delta_phi_cuts.make<TH2F>(b+"h_simgeantParticles_antiS_pt_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts", b+"h_simgeantParticles_antiS_pt_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts; pt Ks (GeV); pt L (GeV)", 2000, -100, 100, 2000, -100, 100);
   
    //Ks daughters of anti-S
    TFileDirectory dir_simgeantParticles_antiS_daughters_Ks = dir_simgeantParticles_antiS_daughters.mkdir("Ks");
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_pdgid"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_pdgid", b+"h_simgeantParticles_antiS_daughter_Ks_pdgid; Ks pdgid", 3, 309,311);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_pt"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_pt", b+"h_simgeantParticles_antiS_daughter_Ks_pt; pt(Ks) (GeV)", 200, 0, 20);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_p"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_p", b+"h_simgeantParticles_antiS_daughter_Ks_p; p(Ks) (GeV)", 200, 0, 20);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_mass"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_mass", b+"h_simgeantParticles_antiS_daughter_Ks_mass; m(Ks) (GeV)", 100, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_eta"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_eta", b+"h_simgeantParticles_antiS_daughter_Ks_eta; #eta(Ks)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_phi"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_phi", b+"h_simgeantParticles_antiS_daughter_Ks_phi; #Phi(Ks) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_vxy"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_vxy", b+"h_simgeantParticles_antiS_daughter_Ks_vxy; vxy(Ks) (cm)", 1000, 0, 100);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_Ks_vxvy"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_Ks_vxvy", b+"h_simgeantParticles_antiS_daughter_Ks_vxvy; vx(Ks)(cm); vy(Ks)(cm)", 2000, -500, 500, 2000, -500, 500);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_Ks_vzvx"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_Ks_vzvx", b+"h_simgeantParticles_antiS_daughter_Ks_vzvx; vz(Ks)(cm); vx(Ks)(cm)", 2000, -500, 500, 2000, -500, 500);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_Ks_vzvy"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_Ks_vzvy", b+"h_simgeantParticles_antiS_daughter_Ks_vzvy; vz(Ks)(cm); vy(Ks)(cm) ", 2000, -500, 500, 2000, -500, 500);

    //Ks daughters of anti-S
    TFileDirectory dir_simgeantParticles_antiS_daughters_L = dir_simgeantParticles_antiS_daughters.mkdir("L");
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_pdgid"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_pdgid", b+"h_simgeantParticles_antiS_daughter_L_pdgid; L pdgid", 3, -3123, -3121);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_pt"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_pt", b+"h_simgeantParticles_antiS_daughter_L_pt; pt(L) (GeV)", 200, 0, 20);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_p"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_p", b+"h_simgeantParticles_antiS_daughter_L_p; p(L) (GeV)", 200, 0, 20);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_mass"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_mass", b+"h_simgeantParticles_antiS_daughter_L_mass; m(L) (GeV)", 100, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_eta"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_eta", b+"h_simgeantParticles_antiS_daughter_L_eta; #eta(L)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_phi"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_phi", b+"h_simgeantParticles_antiS_daughter_L_phi; #Phi(L) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_vxy"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_vxy", b+"h_simgeantParticles_antiS_daughter_L_vxy; vxy(L) (cm)", 1000, 0, 100);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_L_vxvy"]= dir_simgeantParticles_antiS_daughters_L.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_L_vxvy", b+"h_simgeantParticles_antiS_daughter_L_vxvy; vx(L)(cm); vy(L)(cm)", 2000, -500, 500, 2000, -500, 500);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_L_vzvx"]= dir_simgeantParticles_antiS_daughters_L.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_L_vzvx", b+"h_simgeantParticles_antiS_daughter_L_vzvx; vz(L)(cm); vx(L)(cm)", 2000, -500, 500, 2000, -500, 500);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_L_vzvy"]= dir_simgeantParticles_antiS_daughters_L.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_L_vzvy", b+"h_simgeantParticles_antiS_daughter_L_vzvy; vz(L)(cm); vy(L)(cm)", 2000, -500, 500, 2000, -500, 500);






   //For the reconstructed S particles:
    TFileDirectory dir_LambdaKshortVertexFilter = m_fs->mkdir("LambdaKshortVertexFilter");
    TFileDirectory dir_LambdaKshortVertexFilter_S = dir_LambdaKshortVertexFilter.mkdir("S");
    histos_th1f[b+"h_LambdaKshortVertexFilter_S_mass"]= dir_LambdaKshortVertexFilter_S.make<TH1F>(b+"h_LambdaKshortVertexFilter_S_mass", b+"h_LambdaKshortVertexFilter_S_mass; S mass (GeV)", 2000, -100, 100);
    histos_th1f[b+"h_LambdaKshortVertexFilter_Sn_mass"]= dir_LambdaKshortVertexFilter_S.make<TH1F>(b+"h_LambdaKshortVertexFilter_Sn_mass", b+"h_LambdaKshortVertexFilter_Sn_mass; S+n mass (GeV)", 2000, -100, 100);

}


void Analyzer_SIM_Sexaq::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


 
 
  //gen particles 
  edm::Handle<vector<reco::GenParticle>> h_genParticles_GEN;
  iEvent.getByToken(m_genParticlesToken_GEN, h_genParticles_GEN);
 
  //SIM GEANT particles
  edm::Handle<vector<reco::GenParticle>> h_genParticles_SIM_GEANT;
  iEvent.getByToken(m_genParticlesToken_SIM_GEANT, h_genParticles_SIM_GEANT);

  //lambdaKshortVertexFilter sexaquark candidates
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);
 
  //look at the pure GEN particles
  if(h_genParticles_GEN.isValid() )
  {
     for(unsigned int i = 0; i < h_genParticles_GEN->size(); ++i){
       	 histos_th1f[b+"h_genParticles_pdgid"]->Fill(h_genParticles_GEN->at(i).pdgId());

	//only looking at the anti-S
	if(h_genParticles_GEN->at(i).pdgId() == -1020000020){
		
        	histos_th1f[b+"h_genParticles_antiS_pdgid"]->Fill(h_genParticles_GEN->at(i).pdgId());
        	histos_th1f[b+"h_genParticles_antiS_pt"]->Fill(h_genParticles_GEN->at(i).pt());
		
		//look at the daughters of the S particle
		//number of daughters
		histos_th1f[b+"h_genParticles_antiS_n_daughters"]->Fill(h_genParticles_GEN->at(i).numberOfDaughters());

	}
    }
  }

  //look ath the SIM_GEANT particles
  if(h_genParticles_SIM_GEANT.isValid())
  {

      for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i){
	if(h_genParticles_SIM_GEANT->at(i).pdgId() == -1020000020){
		const reco::GenParticle antiS =  h_genParticles_SIM_GEANT->at(i);
		histos_th1f[b+"h_simgeantParticles_antiS_pdgid"]->Fill(antiS.pdgId());
                histos_th1f[b+"h_simgeantParticles_antiS_pt"]->Fill(antiS.pt());
                histos_th1f[b+"h_simgeantParticles_antiS_p"]->Fill(antiS.p());
                histos_th2f[b+"h2_simgeantParticles_antiS_p_eta"]->Fill(antiS.p(), antiS.eta());
                histos_th1f[b+"h_simgeantParticles_antiS_mass"]->Fill(antiS.mass());
                histos_th1f[b+"h_simgeantParticles_antiS_mass_from_Ep"]->Fill(pow(antiS.energy()*antiS.energy()-antiS.p()*antiS.p(),0.5));
                histos_th1f[b+"h_simgeantParticles_antiS_eta"]->Fill(antiS.eta());
                histos_th1f[b+"h_simgeantParticles_antiS_phi"]->Fill(antiS.phi());
                Double_t vx_antiS = antiS.vx();
                Double_t vy_antiS = antiS.vy();
                histos_th1f[b+"h_simgeantParticles_antiS_vxy"]->Fill(pow(vx_antiS*vx_antiS+vy_antiS*vy_antiS,0.5));


		unsigned int n_daughters_antiS = antiS.numberOfDaughters();
                histos_th1f[b+"h_simgeantParticles_antiS_n_daughters"]->Fill(n_daughters_antiS);

		if(n_daughters_antiS == 0){
			
			histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_pt"]->Fill(antiS.pt());
			histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_p"]->Fill(antiS.p());
			histos_th2f[b+"h2_simgeantParticles_antiS_no_daughters_p_eta"]->Fill(antiS.p(), antiS.eta());
			histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_eta"]->Fill(antiS.eta());
			histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_phi"]->Fill(antiS.phi());
			histos_th1f[b+"h_simgeantParticles_antiS_no_daughters_vxy"]->Fill(pow(vx_antiS*vx_antiS+vy_antiS*vy_antiS,0.5));
	
		}
		
		if(n_daughters_antiS == 1){
		        cout << "FOUND AN ANTI-S ONLY DECAYING TO ONE PARTICLE...: "<< antiS.daughter(0)->pdgId() << endl;	

			TLorentzVector p4n(0,0,0,0);
			TLorentzVector p4antiS(0,0,0,0);
			TLorentzVector p4Daughter(0,0,0,0);

			p4n.SetPxPyPzE(0,0,0,0.939565);
			p4antiS.SetPxPyPzE(antiS.px(),antiS.py(),antiS.pz(),antiS.energy());

			TLorentzVector p4antiSn = p4antiS + p4n;

			p4Daughter.SetPxPyPzE(antiS.daughter(0)->px(),antiS.daughter(0)->py(),antiS.daughter(0)->pz(),antiS.daughter(0)->energy());
			cout << "S+n energy: "  << p4antiSn.Energy() << endl;
			cout << "daughter energy: " << p4Daughter.Energy() << endl;
			
			//check if there is anything special with the anti-S that only have 1 daughter 
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_pt"]->Fill(antiS.pt());
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_p"]->Fill(antiS.p());
			histos_th2f[b+"h2_simgeantParticles_antiS_1_daughter_p_eta"]->Fill(antiS.p(), antiS.eta());
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_eta"]->Fill(antiS.eta());
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_phi"]->Fill(antiS.phi());
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_vxy"]->Fill(pow(vx_antiS*vx_antiS+vy_antiS*vy_antiS,0.5));
			
			//check the parameters of the single daugther
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_pdgid_daughter"]->Fill(antiS.daughter(0)->pdgId());
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_pt_daughter"]->Fill(antiS.daughter(0)->pt());
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_p_daughter"]->Fill(antiS.daughter(0)->p());
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_eta_daughter"]->Fill(antiS.daughter(0)->eta());
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_phi_daughter"]->Fill(antiS.daughter(0)->phi());
			histos_th1f[b+"h_simgeantParticles_antiS_1_daughter_vxy_daughter"]->Fill(pow(antiS.daughter(0)->vx()*antiS.daughter(0)->vx()+antiS.daughter(0)->vy()*antiS.daughter(0)->vy(),0.5));
			histos_th2f[b+"h_simgeantParticles_antiS_1_daughter_vzvx_daughter"]->Fill(antiS.daughter(0)->vz(),antiS.daughter(0)->vx());
	
		}
		
		if(n_daughters_antiS == 2){
			
		        cout << "FOUND AN ANTI-S DECAYING TO TWO PARTICLES...: " << endl;	
		        TLorentzVector p4n(0,0,0,0);
                        TLorentzVector p4antiS(0,0,0,0);
                        TLorentzVector p4Daughter0(0,0,0,0);
                        TLorentzVector p4Daughter1(0,0,0,0);

                        p4n.SetPxPyPzE(0,0,0,0.939565);
                        p4antiS.SetPxPyPzE(antiS.px(),antiS.py(),antiS.pz(),antiS.energy());

                        TLorentzVector p4antiSn = p4antiS + p4n;

                        p4Daughter0.SetPxPyPzE(antiS.daughter(0)->px(),antiS.daughter(0)->py(),antiS.daughter(0)->pz(),antiS.daughter(0)->energy());
                        p4Daughter1.SetPxPyPzE(antiS.daughter(1)->px(),antiS.daughter(1)->py(),antiS.daughter(1)->pz(),antiS.daughter(1)->energy());
                        cout << "S+n energy: "  << p4antiSn.Energy() << endl;
                        cout << "daughter0 energy: " << p4Daughter0.Energy() << endl;
                        cout << "daughter1 energy: " << p4Daughter1.Energy() << endl;

			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pt"]->Fill(antiS.pt());
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_p"]->Fill(antiS.p());
			histos_th2f[b+"h2_simgeantParticles_antiS_2_daughters_p_eta"]->Fill(antiS.p(), antiS.eta());
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_eta"]->Fill(antiS.eta());
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_phi"]->Fill(antiS.phi());
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_vxy"]->Fill(pow(vx_antiS*vx_antiS+vy_antiS*vy_antiS,0.5));
	
		}

		if(n_daughters_antiS == 2){

			const reco::Candidate *daughter0 =  antiS.daughter(0);
			const reco::Candidate *daughter1 = antiS.daughter(1);
			//finding out which daughter is which
			const reco::Candidate *Ks = NULL;
			const reco::Candidate *L = NULL;

			if(daughter0->pdgId() == 310) Ks = daughter0;
			else if(daughter0->pdgId() == -3122) L = daughter0;
			else cout << "DAUGHTER0 FROM THE anti-S IS NOT KS OR LAMBDA" << endl;
			
			if(daughter1->pdgId() == 310) Ks = daughter1;
			else if(daughter1->pdgId() == -3122) L = daughter1;
			else cout << "DAUGHTER1 FROM THE anti-S IS NOT KS OR LAMBDA" << endl;
			 //calculating some basics about the daughters
			 Double_t delta_phi_Ks_L = reco::deltaPhi(Ks->phi(), L->phi());
			 Double_t delta_eta_Ks_L = Ks->eta() - L->eta();
			 Double_t delta_R_Ks_L = pow(delta_phi_Ks_L*delta_phi_Ks_L+delta_eta_Ks_L*delta_eta_Ks_L,0.5);

			 Double_t E_Ks = Ks->energy();
			 Double_t E_L = L->energy();
			 math::XYZVector p_Ks = Ks->momentum();
			 math::XYZVector p_L = L->momentum();
			 
			 math::XYZVector p_daughters = Ks->momentum() + L->momentum();
			 Double_t p_daughters_size = pow(p_daughters.X()*p_daughters.X()+p_daughters.Y()*p_daughters.Y()+p_daughters.Z()*p_daughters.Z(),0.5);
			 Double_t invMass_S = pow(pow(E_Ks+E_L-0.939565,2)-pow(p_daughters_size,2),0.5);
			 Double_t opening_angle = TMath::ACos((p_Ks.X()*p_L.X()+p_Ks.Y()*p_L.Y()+p_Ks.Z()*p_L.Z())/(Ks->p()*L->p()));
			 
			 histos_th1f[b+"h_simgeantParticles_antiS_inv_mass_antiS_daughters"]->Fill(invMass_S);
			 histos_th1f[b+"h_simgeantParticles_antiS_openings_angle_daughters"]->Fill(opening_angle); 
			 histos_th2f[b+"h_simgeantParticles_antiS_p_antiS_openings_angle_daughters"]->Fill(antiS.p(), opening_angle); 
			 histos_th2f[b+"h_simgeantParticles_antiS_pt_antiS_openings_angle_daughters"]->Fill(antiS.pt(), opening_angle); 
		 	 histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters"]->Fill(delta_phi_Ks_L);
			 histos_th1f[b+"h_simgeantParticles_antiS_delta_eta_daughters"]->Fill(delta_eta_Ks_L);
			 histos_th1f[b+"h_simgeantParticles_antiS_delta_R_daughters"]->Fill(delta_R_Ks_L);
			 //just check wether the delta phi is smaller when the S has higher momentum 
			 if(antiS.pt()>1)histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters_S_pt_larger_1GeV"]->Fill(delta_phi_Ks_L);
			
  			 //define the pt cuts noramlly done when running the skimming 
			 bool KsSurvPtCut = false;
			 bool LSurvPtCut = false;
			 if(Ks->pt() > 0.9) KsSurvPtCut = true;
			 if(L->pt() > 1.5) LSurvPtCut = true;
			 
			 if(KsSurvPtCut && LSurvPtCut){
				histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters_pt_cuts"]->Fill(delta_phi_Ks_L);
                         	histos_th1f[b+"h_simgeantParticles_antiS_delta_eta_daughters_pt_cuts"]->Fill(delta_eta_Ks_L);
                         	histos_th2f[b+"h_simgeantParticles_antiS_delta_eta_delta_phi_daughters_pt_cuts"]->Fill(delta_eta_Ks_L, delta_phi_Ks_L);
                         	histos_th1f[b+"h_simgeantParticles_antiS_delta_R_daughters_pt_cuts"]->Fill(delta_R_Ks_L);
			 }

			 //check wether the daughters of the Ks and the Lambda are within tracker accepatance. This makes them potentially reconstructible.
			 int n_daughters_Ks =  Ks->numberOfDaughters();
			 int n_daughters_L =  L->numberOfDaughters();
			 bool KsandLDaughtersSurvEtaCut = false;
			 if(n_daughters_Ks == 2 && n_daughters_L == 2){
				if(fabs(Ks->daughter(0)->eta()) < 2.5 && fabs(Ks->daughter(1)->eta()) && fabs(L->daughter(0)->eta()) && fabs(L->daughter(1)->eta())) KsandLDaughtersSurvEtaCut = true;
			 }

			 bool antiSSurvDeltaPhiCut = false;
			 if(fabs(delta_phi_Ks_L)>0.5) antiSSurvDeltaPhiCut = true;
			
			//boost back to the S+n reference frame. You have the S particle so still to add the S+n particle (S+n is the COM frame). When you do this on the RECO data the particle you reconstruct is actually the S+n, so there you do not need to add the n.  
			TLorentzVector p4n(0,0,0,0);
			TLorentzVector p4antiS(0,0,0,0);
			TLorentzVector p4Ks(0,0,0,0);
			TLorentzVector p4L(0,0,0,0);

			p4n.SetPxPyPzE(0,0,0,0.939565);
			p4antiS.SetPxPyPzE(antiS.px(),antiS.py(),antiS.pz(),antiS.energy());
			TLorentzVector p4antiSn = p4antiS + p4n;
			p4Ks.SetPxPyPzE(Ks->px(),Ks->py(),Ks->pz(),Ks->energy());
			p4L.SetPxPyPzE(L->px(),L->py(),L->pz(),L->energy());
			p4Ks.Boost(-p4antiSn.BoostVector());
			p4L.Boost(-p4antiSn.BoostVector());
			TLorentzVector p4_Ks_boosted_COM = p4Ks;
			TLorentzVector p4_L_boosted_COM = p4L;
			p4Ks.Boost(p4antiSn.BoostVector());
                        p4L.Boost(p4antiSn.BoostVector());

			 if(KsSurvPtCut && LSurvPtCut && KsandLDaughtersSurvEtaCut && antiSSurvDeltaPhiCut){
			
				if(fabs(delta_phi_Ks_L) > 0.5 && fabs(delta_eta_Ks_L) < 2 && openings_angle(Ks->momentum(), L->momentum()) > 0.5 && openings_angle(Ks->momentum(), L->momentum()) < 2 && delta_R_Ks_L < 3 && p4_L_boosted_COM.Pz() > -p4_Ks_boosted_COM.Pz() - 0.3 && p4_L_boosted_COM.Pz() < -p4_Ks_boosted_COM.Pz() + 0.3){
					 histos_th1f[b+"h_simgeantParticles_antiS_n_surv_background_cuts"]->Fill(1);
				}
				else{
					 histos_th1f[b+"h_simgeantParticles_antiS_n_surv_background_cuts"]->Fill(0);
				}
				cout << "fabs(delta_phi_Ks_L) " << fabs(delta_phi_Ks_L) << endl;
				cout << "fabs(delta_eta_Ks_L) " << fabs(delta_eta_Ks_L) << endl;
				cout << "openings_angle(Ks->momentum(), L->momentum()) " << openings_angle(Ks->momentum(), L->momentum()) << endl;
				cout << "delta_R_Ks_L " <<  delta_R_Ks_L << endl;
				//antiS particle				
				histos_th1f[b+"h_simgeantParticles_antiS_p_pt_and_eta_and_delta_phi_cuts"]->Fill(antiS.p());
				histos_th1f[b+"h_simgeantParticles_antiS_pt_pt_and_eta_and_delta_phi_cuts"]->Fill(antiS.pt());
				histos_th1f[b+"h_simgeantParticles_antiS_phi_pt_and_eta_and_delta_phi_cuts"]->Fill(antiS.phi());
				histos_th1f[b+"h_simgeantParticles_antiS_eta_pt_and_eta_and_delta_phi_cuts"]->Fill(antiS.eta());
				//Ks
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_p_pt_and_eta_and_delta_phi_cuts"]->Fill(Ks->p());
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_pt_pt_and_eta_and_delta_phi_cuts"]->Fill(Ks->pt());
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_phi_pt_and_eta_and_delta_phi_cuts"]->Fill(Ks->phi());
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_eta_pt_and_eta_and_delta_phi_cuts"]->Fill(Ks->eta());
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_vxy_pt_and_eta_and_delta_phi_cuts"]->Fill(pow(Ks->vx()*Ks->vx()+Ks->vy()*Ks->vy(),0.5));	
				//L
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_p_pt_and_eta_and_delta_phi_cuts"]->Fill(L->p());
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_pt_pt_and_eta_and_delta_phi_cuts"]->Fill(L->pt());
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_phi_pt_and_eta_and_delta_phi_cuts"]->Fill(L->phi());
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_eta_pt_and_eta_and_delta_phi_cuts"]->Fill(L->eta());
				histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_vxy_pt_and_eta_and_delta_phi_cuts"]->Fill(pow(L->vx()*L->vx()+L->vy()*L->vy(),0.5));
				//Ks and eta correlation
				histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters_pt_and_eta_and_delta_phi_cuts"]->Fill(delta_phi_Ks_L);
                                histos_th1f[b+"h_simgeantParticles_antiS_delta_eta_daughters_pt_and_eta_and_delta_phi_cuts"]->Fill(delta_eta_Ks_L);
                                histos_th2f[b+"h_simgeantParticles_antiS_delta_eta_delta_phi_daughters_pt_and_eta_and_delta_phi_cuts"]->Fill(delta_eta_Ks_L, delta_phi_Ks_L);
                                histos_th1f[b+"h_simgeantParticles_antiS_delta_R_daughters_pt_and_eta_and_delta_phi_cuts"]->Fill(delta_R_Ks_L);
                                histos_th1f[b+"h_simgeantParticles_antiS_openings_angle_daughters_pt_and_eta_and_delta_phi_cuts"]->Fill(openings_angle(Ks->momentum(), L->momentum()));
                                histos_th2f[b+"h_simgeantParticles_antiS_pt_corr_daughters_pt_and_eta_and_delta_phi_cuts"]->Fill(Ks->pt(),L->pt());
                                histos_th2f[b+"h_simgeantParticles_antiS_p_corr_daughters_pt_and_eta_and_delta_phi_cuts"]->Fill(Ks->p(),L->p());
				//Ks and L boosted back to S+n reference frame
                                histos_th2f[b+"h_simgeantParticles_antiS_px_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]->Fill(p4_Ks_boosted_COM.Px(),p4_L_boosted_COM.Px());
                                histos_th2f[b+"h_simgeantParticles_antiS_py_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]->Fill(p4_Ks_boosted_COM.Py(),p4_L_boosted_COM.Py());
                                histos_th2f[b+"h_simgeantParticles_antiS_pz_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]->Fill(p4_Ks_boosted_COM.Pz(),p4_L_boosted_COM.Pz());
                                histos_th2f[b+"h_simgeantParticles_antiS_p_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]->Fill(p4_Ks_boosted_COM.P(),p4_L_boosted_COM.P());
                                histos_th2f[b+"h_simgeantParticles_antiS_pt_corr_daughters_COM_pt_and_eta_and_delta_phi_cuts"]->Fill(p4_Ks_boosted_COM.Pt(),p4_L_boosted_COM.Pt());
					
			 }
		
		//reference kinematics for all the Ks and anti lambda	
		for(unsigned int j = 0; j < n_daughters_antiS; j++ ){
                    int pdgId_antiS_daughter = h_genParticles_SIM_GEANT->at(i).daughter(j)->pdgId();
                    histos_th1f[b+"h_simgeantParticles_antiS_pdgid_daughters"]->Fill(pdgId_antiS_daughter);
                    
                    //if daughter is Ks
                    if(pdgId_antiS_daughter == 310){
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_pdgid"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->pdgId());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->pt());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_p"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->p());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_mass"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->mass());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->eta());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->phi());
                        Double_t vx_Ks = h_genParticles_SIM_GEANT->at(i).daughter(j)->vx();
                        Double_t vy_Ks = h_genParticles_SIM_GEANT->at(i).daughter(j)->vy();
                        Double_t vz_Ks = h_genParticles_SIM_GEANT->at(i).daughter(j)->vz();
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_vxy"]->Fill(pow(vx_antiS*vx_antiS+vy_antiS*vy_antiS,0.5));
                        histos_th2f[b+"h_simgeantParticles_antiS_daughter_Ks_vxvy"]->Fill(vx_Ks,vy_Ks);
                        histos_th2f[b+"h_simgeantParticles_antiS_daughter_Ks_vzvx"]->Fill(vz_Ks,vx_Ks);
                        histos_th2f[b+"h_simgeantParticles_antiS_daughter_Ks_vzvy"]->Fill(vz_Ks,vy_Ks);
                    }
                    //if daugther is L
                    if(pdgId_antiS_daughter == -3122){
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_pdgid"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->pdgId());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->pt());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_p"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->p());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_mass"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->mass());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->eta());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->phi());
                        Double_t vx_L = h_genParticles_SIM_GEANT->at(i).daughter(j)->vx();
                        Double_t vy_L = h_genParticles_SIM_GEANT->at(i).daughter(j)->vy();
                        Double_t vz_L = h_genParticles_SIM_GEANT->at(i).daughter(j)->vz();
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_vxy"]->Fill(pow(vx_antiS*vx_antiS+vy_antiS*vy_antiS,0.5));
                        histos_th2f[b+"h_simgeantParticles_antiS_daughter_L_vxvy"]->Fill(vx_L,vy_L);
                        histos_th2f[b+"h_simgeantParticles_antiS_daughter_L_vzvx"]->Fill(vz_L,vx_L);
                        histos_th2f[b+"h_simgeantParticles_antiS_daughter_L_vzvy"]->Fill(vz_L,vy_L);

                    }
                }		
		
	}// if(n_daughters_antiS == 2)
      }//if(h_genParticles_SIM_GEANT->at(i).pdgId() == -1020000020)

    }//for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i)
  }//if(h_genParticles_SIM_GEANT.isValid())
 

  if(h_sCands.isValid()) {
        unsigned int n_sCands = h_sCands->size();
        for (unsigned int i = 0; i < n_sCands; ++i) { 
		//only look at the candidates with a charge == 1 (this is the charge of the proton), so you are looking at S candidates, not the anti-S candidates, so you are looking at background. 
		if(h_sCands->at(i).charge() == -1)continue;
			//the reconstructed particle is actually the S+n particle so have to subtract the neutron from the S candidate p4
			Double_t Sn_mass = h_sCands->at(i).mass();
			TLorentzVector p4n(0,0,0,0);
			p4n.SetPxPyPzE(0,0,0,0.939565);
			TLorentzVector p4Sn(0,0,0,0);
			p4Sn.SetPxPyPzE(h_sCands->at(i).px(),h_sCands->at(i).py(),h_sCands->at(i).pz(),h_sCands->at(i).energy());
			TLorentzVector p4S = p4Sn-p4n;
			Double_t S_mass = p4S.M();

			histos_th1f[b+"h_LambdaKshortVertexFilter_S_mass"]->Fill(S_mass);
			histos_th1f[b+"h_LambdaKshortVertexFilter_Sn_mass"]->Fill(Sn_mass);
			
		
		
			
	}//for (unsigned int i = 0; i < n_sCands; ++i)
  }//if(h_sCands.isValid())

} //end of analyzer

double Analyzer_SIM_Sexaq::openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2){
  double opening_angle = TMath::ACos((momentum1.Dot(momentum2))/(pow(momentum1.Mag2()*momentum2.Mag2(),0.5)));
  return opening_angle;
}

void Analyzer_SIM_Sexaq::endJob()
{
}

Analyzer_SIM_Sexaq::~Analyzer_SIM_Sexaq()
{
}


DEFINE_FWK_MODULE(Analyzer_SIM_Sexaq);
