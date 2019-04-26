#include "../interface/Analyzer_SIM_Sexaq.h"
#include <typeinfo>


Analyzer_SIM_Sexaq::Analyzer_SIM_Sexaq(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  m_bsTag(pset.getParameter<edm::InputTag>("beamspot")),
  m_offlinePVTag(pset.getParameter<edm::InputTag>("offlinePV")),
  m_genParticlesTag_GEN(pset.getParameter<edm::InputTag>("genCollection_GEN")),
  m_genParticlesTag_SIM_GEANT(pset.getParameter<edm::InputTag>("genCollection_SIM_GEANT")),
  m_generalTracksTag(pset.getParameter<edm::InputTag>("generalTracksCollection")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_V0KsTag(pset.getParameter<edm::InputTag>("V0KsCollection")),
  m_V0LTag(pset.getParameter<edm::InputTag>("V0LCollection")),

  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_offlinePVToken    (consumes<vector<reco::Vertex>>(m_offlinePVTag)),
  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN)),
  m_genParticlesToken_SIM_GEANT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_SIM_GEANT)),
  m_generalTracksToken(consumes<vector<reco::Track> >(m_generalTracksTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_V0KsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0KsTag)),
  m_V0LToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0LTag))


{

}


void Analyzer_SIM_Sexaq::beginJob() {
 
    //for the tracks
    TFileDirectory dir_generalTracks = m_fs->mkdir("generalTracks");
    histos_th1f[b+"h_tracks_deltaR"]= dir_generalTracks.make<TH1F>(b+"h_tracks_deltaR", b+"h_tracks_deltaR; deltaR(track, daughter of V0) ",1000,0,10);
    histos_th1f[b+"h_tracks_deltaR_grandDaughtersAntiS"]= dir_generalTracks.make<TH1F>(b+"h_tracks_deltaR_grandDaughtersAntiS", b+"h_tracks_deltaR_grandDaughtersAntiS; deltaR(track, grandDaughter of antiS) ",1000,0,10);
    histos_th1f[b+"h_tracks_deltaR_min"]= dir_generalTracks.make<TH1F>(b+"h_tracks_deltaR_min", b+"h_tracks_deltaR_min; deltaR min(track, daughter of V0) ",1000,0,10);
    histos_th1f[b+"h_tracks_deltaR_min_grandDaughtersAntiS"]= dir_generalTracks.make<TH1F>(b+"h_tracks_deltaR_min_grandDaughtersAntiS", b+"h_tracks_deltaR_min_grandDaughtersAntiS; deltaR min(track, daughter of V0) ",1000,0,10);
   histos_teff[b+"tracks_reco_eff_daughters_antiS"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_daughters_antiS",b+"tracks_reco_eff_daughters_antiS; 0 = tracks no daughters of antiS, 1 = tracks daughters of the antiS",2,0,2);
   histos_teff[b+"tracks_reco_eff_daughters_antiS_pt"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_daughters_antiS_pt",b+"tracks_reco_eff_daughters_antiS_pt; pt of V0 daughters (GeV)",100,0,10);
   histos_teff[b+"tracks_reco_eff_daughters_antiS_eta"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_daughters_antiS_eta",b+"tracks_reco_eff_daughters_antiS_eta; eta of V0 daughters",100,-4,4);
   histos_teff[b+"tracks_reco_eff_daughters_antiS_vxy"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_daughters_antiS_vxy",b+"tracks_reco_eff_daughters_antiS_vxy; vxy of V0 daughters (cm)",200,0,100);
   histos_teff[b+"tracks_reco_eff_daughters_antiS_surv_track_cuts_pt"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_daughters_antiS_surv_track_cuts_pt",b+"tracks_reco_eff_daughters_antiS_surv_track_cuts_pt; pt of V0 daughters (GeV)",100,0,10);
   histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_no_daughters_antiS_pt",b+"tracks_reco_eff_no_daughters_antiS_pt; pt of V0 daughters (GeV)",100,0,10);
   histos_teff[b+"tracks_reco_eff_no_daughters_antiS_eta"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_no_daughters_antiS_eta",b+"tracks_reco_eff_no_daughters_antiS_eta; eta of V0 daughters",100,-4,4);
   histos_teff[b+"tracks_reco_eff_no_daughters_antiS_vxy"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_no_daughters_antiS_vxy",b+"tracks_reco_eff_no_daughters_antiS_vxy; vxy of V0 daughters (cm)",200,0,100);
   histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5",b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5; pt of V0 daughters (GeV)",100,0,10);

    //for the GEN particles
    TFileDirectory dir_genParticles = m_fs->mkdir("genParticles");
    histos_th1f[b+"h_genParticles_pdgid"]= dir_genParticles.make<TH1F>(b+"h_genParticles_pdgid", b+"h_genParticles_pdgid; pdgId of Gen particles; pdgId ",6000,-3000,3000);
    histos_th1f[b+"h_genParticles_antiS_mass_check"]= dir_genParticles.make<TH1F>(b+"h_genParticles_antiS_mass_check", b+"h_genParticles_antiS_mass_check; mass S ",200,-10,10);
    histos_th1f[b+"h_genParticles_Xi_mass_check"]= dir_genParticles.make<TH1F>(b+"h_genParticles_Xi_mass_check", b+"h_genParticles_Xi_mass_check; mass Xi ",200,-10,10);
    histos_th1f[b+"h_genParticles_antiS_angular_matching_mass_check_GEN"]= dir_genParticles.make<TH1F>(b+"h_genParticles_antiS_angular_matching_mass_check_GEN", b+"h_genParticles_antiS_angular_matching_mass_check_GEN; mass S ",200,-10,10);
    histos_th1f[b+"h_genParticles_Xi_mass_angular_matching_check_GEN"]= dir_genParticles.make<TH1F>(b+"h_genParticles_Xi_mass_angular_matching_check_GEN", b+"h_genParticles_Xi_mass_angular_matching_check_GEN; mass Xi ",200,-10,10);
    histos_th1f[b+"h_genParticles_antiS_angular_matching_mass_check_RECO"]= dir_genParticles.make<TH1F>(b+"h_genParticles_antiS_angular_matching_mass_check_RECO", b+"h_genParticles_antiS_angular_matching_mass_check_RECO; mass S ",200,-10,10);
    histos_th1f[b+"h_genParticles_Xi_mass_angular_matching_check_RECO"]= dir_genParticles.make<TH1F>(b+"h_genParticles_Xi_mass_angular_matching_check_RECO", b+"h_genParticles_Xi_mass_angular_matching_check_RECO; mass Xi ",200,-10,10);




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
    histos_th1f[b+"h_simgeantParticles_antiS_status"]= dir_simgeantParticles_antiS.make<TH1F>(b+"h_simgeantParticles_antiS_status", b+"h_simgeantParticles_antiS_status; anti-S status; ",10, 0, 10);
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
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_lxy_V0s"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_lxy_V0s", b+"h_simgeantParticles_antiS_2_daughters_lxy_V0s; lxy (V0 vertex, beamspot) (cm)", 2000, 0, 200);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_vz_V0s"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_vz_V0s", b+"h_simgeantParticles_antiS_2_daughters_vz_V0s; vz (V0 vertex) (cm)", 3000, -300, 300);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_deltaPhi_V0s"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_deltaPhi_V0s", b+"h_simgeantParticles_antiS_2_daughters_deltaPhi_V0s; #Delta Phi (Ks, Lambda) (cm)", 100, -4, 4);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_deltaEta_V0s"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_deltaEta_V0s", b+"h_simgeantParticles_antiS_2_daughters_deltaEta_V0s; #Delta Eta (Ks, Lambda) (cm)", 100, -4, 4);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_deltaR_V0s"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_deltaR_V0s", b+"h_simgeantParticles_antiS_2_daughters_deltaR_V0s; #Delta R (Ks, Lambda) (cm)", 80, 0, 8);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pdgId_granddaughters"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_pdgId_granddaughters", b+"h_simgeantParticles_antiS_2_daughters_pdgId_granddaughters; pdgId antiS granddaughters (cm)", 20000, -10000, 10000);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_lxy_granddaughters"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_lxy_granddaughters", b+"h_simgeantParticles_antiS_2_daughters_lxy_granddaughters; lxy antiS granddaughters", 1000, 0, 100);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pt_Ksdaughters"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_pt_Ksdaughters", b+"h_simgeantParticles_antiS_2_daughters_pt_Ksdaughters; pt Ks granddaughters", 100, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pt_AntiLdaughters"]= dir_simgeantParticles_antiS_2_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_2_daughters_pt_AntiLdaughters", b+"h_simgeantParticles_antiS_2_daughters_pt_AntiLdaughters; pt antiL granddaughters", 100, 0, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_n_daughters_Ks_n_daughters_L"]= dir_simgeantParticles_antiS_2_daughters.make<TH2F>(b+"h_simgeantParticles_antiS_n_daughters_Ks_n_daughters_L", b+"h_simgeantParticles_antiS_n_daughters_Ks_n_daughters_L; number of daughters of the Ks; number of daughters of the L;", 20, 0, 20, 20, 0, 20);
    histos_th2f[b+"h_simgeantParticles_antiS_2_daughters_pdgId_Ks_daughters_absdiff_pdgId_L_daughters_absdiff_pdgId"]= dir_simgeantParticles_antiS_2_daughters.make<TH2F>(b+"h_simgeantParticles_antiS_2_daughters_pdgId_Ks_daughters_absdiff_pdgId_L_daughters_absdiff_pdgId", b+"h_simgeantParticles_antiS_2_daughters_pdgId_Ks_daughters_absdiff_pdgId_L_daughters_absdiff_pdgId; absdiff of the Ks daughters pdgID; absdiff of the L daughters pdgID;", 3000, 0, 3000, 3000, 0, 3000);
	
    //check on RECO level how many daughters of the antiS are reconstructed
    TFileDirectory dir_simgeantParticles_antiS_2_daughters_V0s = dir_simgeantParticles_antiS_2_daughters.mkdir("V0s");
    histos_th2f[b+"h_V0_nantiSLDaughtersReconstucted_nantiSKsDaughtersReconstucted"]= dir_simgeantParticles_antiS_2_daughters_V0s.make<TH2F>(b+"h_V0_nantiSLDaughtersReconstucted_nantiSKsDaughtersReconstucted", b+"h_V0_nantiSLDaughtersReconstucted_nantiSKsDaughtersReconstucted; matched V0 RECO Ks with GEN Ks antiS daughter;  matched V0 RECO L with GEN L antiS daughter ", 10, 0, 10, 10, 0, 10);
	
    //for the daughters of the S
    TFileDirectory dir_simgeantParticles_antiS_daughters = dir_simgeantParticles_antiS.mkdir("daughters");
    histos_th1f[b+"h_simgeantParticles_antiS_n_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_n_daughters", b+"h_simgeantParticles_antiS_n_daughters; #daughters of anti S Gen particles; ",20,0,20);
    histos_th1f[b+"h_simgeantParticles_antiS_pdgid_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_pdgid_daughters", b+"h_simgeantParticles_antiS_pdgid_daughters; pdgid of anti-S daughters", 8000, -4000, 4000);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_delta_phi_daughters", b+"h_simgeantParticles_antiS_delta_phi_daughters; #Delta Phi(Ks,L) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters_S_pt_larger_1GeV"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_delta_phi_daughters_S_pt_larger_1GeV", b+"h_simgeantParticles_antiS_delta_phi_daughters_S_pt_larger_1GeV; #Delta Phi(Ks,L) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_eta_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_delta_eta_daughters", b+"h_simgeantParticles_antiS_delta_eta_daughters; #Delta eta(Ks,L) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_delta_R_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_delta_R_daughters", b+"h_simgeantParticles_antiS_delta_R_daughters; #Delta R(Ks,L)", 100, 0, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_pt_corr_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH2F>(b+"h_simgeantParticles_antiS_pt_corr_daughters", b+"h_simgeantParticles_antiS_pt_corr_daughters; pt(Ks); pt(L)", 100, 0, 10, 100, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_inv_mass_antiS_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_inv_mass_antiS_daughters", b+"h_simgeantParticles_antiS_inv_mass_antiS_daughters; inv mass anti-S (GeV)", 2000, 1.7, 1.9);
    histos_th1f[b+"h_simgeantParticles_antiS_openings_angle_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH1F>(b+"h_simgeantParticles_antiS_openings_angle_daughters", b+"h_simgeantParticles_antiS_openings_angle_daughters; openings angle", 200, -10, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_pt_antiS_openings_angle_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH2F>(b+"h_simgeantParticles_antiS_pt_antiS_openings_angle_daughters", b+"h_simgeantParticles_antiS_pt_antiS_openings_angle_daughters;pt(S)(GeV); openings angle", 100, 0, 20, 200, -10, 10);
    histos_th2f[b+"h_simgeantParticles_antiS_p_antiS_openings_angle_daughters"]= dir_simgeantParticles_antiS_daughters.make<TH2F>(b+"h_simgeantParticles_antiS_p_antiS_openings_angle_daughters", b+"h_simgeantParticles_antiS_p_antiS_openings_angle_daughters; p(S)(GeV);openings angle", 100, 0, 100, 200, -10, 10);
    //with pt cuts on the Ks (pt > 0.9GeV) and the L(pt > 1.5GeV) 
   TFileDirectory dir_simgeantParticles_antiS_daughters_pt_cuts = dir_simgeantParticles_antiS.mkdir("pt_cuts");
   histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]= dir_simgeantParticles_antiS_daughters_pt_cuts.make<TH1F>(b+"h_simgeantParticles_antiS_cut_flow", b+"h_simgeantParticles_antiS_cut_flow; ", 15, 0, 15); 
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
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_status"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_status", b+"h_simgeantParticles_antiS_daughter_Ks_status; status", 10, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_pt"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_pt", b+"h_simgeantParticles_antiS_daughter_Ks_pt; pt(Ks) (GeV)", 200, 0, 20);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_p"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_p", b+"h_simgeantParticles_antiS_daughter_Ks_p; p(Ks) (GeV)", 200, 0, 20);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_mass"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_mass", b+"h_simgeantParticles_antiS_daughter_Ks_mass; m(Ks) (GeV)", 100, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_eta"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_eta", b+"h_simgeantParticles_antiS_daughter_Ks_eta; #eta(Ks)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_phi"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_phi", b+"h_simgeantParticles_antiS_daughter_Ks_phi; #Phi(Ks) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_vxy"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_vxy", b+"h_simgeantParticles_antiS_daughter_Ks_vxy; vxy(Ks) (cm)", 1000, 0, 100);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_Ks_vxvy"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_Ks_vxvy", b+"h_simgeantParticles_antiS_daughter_Ks_vxvy; vx(Ks)(cm); vy(Ks)(cm)", 2000, -500, 500, 2000, -500, 500);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_Ks_vzvx"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_Ks_vzvx", b+"h_simgeantParticles_antiS_daughter_Ks_vzvx; vz(Ks)(cm); vx(Ks)(cm)", 2000, -500, 500, 2000, -500, 500);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_Ks_vzvy"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_Ks_vzvy", b+"h_simgeantParticles_antiS_daughter_Ks_vzvy; vz(Ks)(cm); vy(Ks)(cm) ", 2000, -500, 500, 2000, -500, 500);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_n_daughters"]= dir_simgeantParticles_antiS_daughters_Ks.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_Ks_n_daughters", b+"h_simgeantParticles_antiS_daughter_Ks_n_daughters; number of daughters", 10, 0, 10);

    //Ks daughters of anti-S
    TFileDirectory dir_simgeantParticles_antiS_daughters_L = dir_simgeantParticles_antiS_daughters.mkdir("L");
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_pdgid"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_pdgid", b+"h_simgeantParticles_antiS_daughter_L_pdgid; L pdgid", 3, -3123, -3121);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_status"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_status", b+"h_simgeantParticles_antiS_daughter_L_status; status", 10, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_pt"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_pt", b+"h_simgeantParticles_antiS_daughter_L_pt; pt(L) (GeV)", 200, 0, 20);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_p"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_p", b+"h_simgeantParticles_antiS_daughter_L_p; p(L) (GeV)", 200, 0, 20);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_mass"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_mass", b+"h_simgeantParticles_antiS_daughter_L_mass; m(L) (GeV)", 100, 0, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_eta"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_eta", b+"h_simgeantParticles_antiS_daughter_L_eta; #eta(L)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_phi"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_phi", b+"h_simgeantParticles_antiS_daughter_L_phi; #Phi(L) (rad)", 200, -10, 10);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_vxy"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_vxy", b+"h_simgeantParticles_antiS_daughter_L_vxy; vxy(L) (cm)", 1000, 0, 100);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_L_vxvy"]= dir_simgeantParticles_antiS_daughters_L.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_L_vxvy", b+"h_simgeantParticles_antiS_daughter_L_vxvy; vx(L)(cm); vy(L)(cm)", 2000, -500, 500, 2000, -500, 500);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_L_vzvx"]= dir_simgeantParticles_antiS_daughters_L.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_L_vzvx", b+"h_simgeantParticles_antiS_daughter_L_vzvx; vz(L)(cm); vx(L)(cm)", 2000, -500, 500, 2000, -500, 500);
    histos_th2f[b+"h_simgeantParticles_antiS_daughter_L_vzvy"]= dir_simgeantParticles_antiS_daughters_L.make<TH2F>(b+"h_simgeantParticles_antiS_daughter_L_vzvy", b+"h_simgeantParticles_antiS_daughter_L_vzvy; vz(L)(cm); vy(L)(cm)", 2000, -500, 500, 2000, -500, 500);
    histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_n_daughters"]= dir_simgeantParticles_antiS_daughters_L.make<TH1F>(b+"h_simgeantParticles_antiS_daughter_L_n_daughters", b+"h_simgeantParticles_antiS_daughter_L_n_daughters; number of daughters", 10, 0, 10);


   //V0 particles
   TFileDirectory dir_V0s = m_fs->mkdir("V0s");
   //Ks 
   //kinematics of the GEN Ks
   TFileDirectory dir_V0s_GEN_Ks_kinematics = dir_V0s.mkdir("GEN_Ks_kinematics");
   histos_th1f[b+"GEN_Ks_pt"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_pt", b+"GEN_Ks_pt; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_Ks_pt_daughter_antiS"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_pt_daughter_antiS", b+"GEN_Ks_pt_daughter_antiS; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_Ks_daughter_antiS_pt_of_Ks_daughters"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_daughter_antiS_pt_of_Ks_daughters", b+"GEN_Ks_daughter_antiS_pt_of_Ks_daughters; pt(GeV)", 100, 0, 10);
   histos_th2f[b+"GEN_Ks_daughter_antiS_min_max_pt_of_Ks_daughters"]= dir_V0s_GEN_Ks_kinematics.make<TH2F>(b+"GEN_Ks_daughter_antiS_min_max_pt_of_Ks_daughters", b+"GEN_Ks_daughter_antiS_min_max_pt_of_Ks_daughters; daughter Ks min pt(GeV); daughter Ks max pt(GeV)", 100, 0, 10, 100, 0, 10);
   histos_th1f[b+"GEN_Ks_eta"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_eta", b+"GEN_Ks_eta; eta", 100, -4, 4);
   histos_th1f[b+"GEN_Ks_phi"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_phi", b+"GEN_Ks_phi; phi(rad)", 200, -10, 10);
   histos_th1f[b+"GEN_Ks_lxy"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_lxy", b+"GEN_Ks_lxy; lxy(cm)", 1000, 0, 100);
   histos_th2f[b+"GEN_Ks_lxy_pt"]= dir_V0s_GEN_Ks_kinematics.make<TH2F>(b+"GEN_Ks_lxy_pt", b+"GEN_Ks_lxy_pt; lxy(cm); pt(GeV)", 1000, 0, 100, 100, 0, 10);
   histos_th1f[b+"GEN_Ks_vz"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_vz", b+"GEN_Ks_vz; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"GEN_Ks_ndaughters"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_ndaughters", b+"GEN_Ks_ndaughters; #daughters", 10, 0, 10);
   histos_th1f[b+"GEN_Ks_status"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_status", b+"GEN_Ks_status; status", 10, 0, 10);
   histos_th2f[b+"GEN_Ks_status_n_daughters"]= dir_V0s_GEN_Ks_kinematics.make<TH2F>(b+"GEN_Ks_status_n_daughters", b+"GEN_Ks_status_n_daughters; status; n daughters", 10, 0, 10, 10, 0, 10);

   //kinematics of the GEN with status 2
   TFileDirectory dir_V0s_GEN_Ks_kinematics_status1 = dir_V0s.mkdir("GEN_Ks_kinematics_status1");
   histos_th1f[b+"GEN_Ks_pt_status1"]= dir_V0s_GEN_Ks_kinematics_status1.make<TH1F>(b+"GEN_Ks_pt_status1", b+"GEN_Ks_pt_status1; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_Ks_eta_status1"]= dir_V0s_GEN_Ks_kinematics_status1.make<TH1F>(b+"GEN_Ks_eta_status1", b+"GEN_Ks_eta_status1; eta", 100, -4, 4);
   histos_th1f[b+"GEN_Ks_phi_status1"]= dir_V0s_GEN_Ks_kinematics_status1.make<TH1F>(b+"GEN_Ks_phi_status1", b+"GEN_Ks_phi_status1; phi(rad)", 200, -10, 10);
   histos_th1f[b+"GEN_Ks_lxy_status1"]= dir_V0s_GEN_Ks_kinematics_status1.make<TH1F>(b+"GEN_Ks_lxy_status1", b+"GEN_Ks_lxy_status1; lxy(cm)", 1000, 0, 100);
   histos_th2f[b+"GEN_Ks_lxy_pt_status1"]= dir_V0s_GEN_Ks_kinematics_status1.make<TH2F>(b+"GEN_Ks_lxy_pt_status1", b+"GEN_Ks_lxy_pt_status1; lxy(cm); pt(GeV)", 1000, 0, 100, 100, 0, 10);
   histos_th1f[b+"GEN_Ks_vz_status1"]= dir_V0s_GEN_Ks_kinematics_status1.make<TH1F>(b+"GEN_Ks_vz_status1", b+"GEN_Ks_vz_status1; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"GEN_Ks_ndaughters_status1"]= dir_V0s_GEN_Ks_kinematics_status1.make<TH1F>(b+"GEN_Ks_ndaughters_status1", b+"GEN_Ks_ndaughters_status1; #daughters", 10, 0, 10);
   histos_th1f[b+"GEN_Ks_status_status1"]= dir_V0s_GEN_Ks_kinematics_status1.make<TH1F>(b+"GEN_Ks_status_status1", b+"GEN_Ks_status_status1; status", 10, 0, 10);
   histos_th2f[b+"GEN_Ks_status_n_daughters_status1"]= dir_V0s_GEN_Ks_kinematics_status1.make<TH2F>(b+"GEN_Ks_status_n_daughters_status1", b+"GEN_Ks_status_n_daughters_status1; status; n daughters", 10, 0, 10, 10, 0, 10);


   //kinematics of the GEN with status 8
   TFileDirectory dir_V0s_GEN_Ks_kinematics_status8 = dir_V0s.mkdir("GEN_Ks_kinematics_status8");
   histos_th1f[b+"GEN_Ks_pt_status8"]= dir_V0s_GEN_Ks_kinematics_status8.make<TH1F>(b+"GEN_Ks_pt_status8", b+"GEN_Ks_pt_status8; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_Ks_eta_status8"]= dir_V0s_GEN_Ks_kinematics_status8.make<TH1F>(b+"GEN_Ks_eta_status8", b+"GEN_Ks_eta_status8; eta", 100, -4, 4);
   histos_th1f[b+"GEN_Ks_phi_status8"]= dir_V0s_GEN_Ks_kinematics_status8.make<TH1F>(b+"GEN_Ks_phi_status8", b+"GEN_Ks_phi_status8; phi(rad)", 200, -10, 10);
   histos_th1f[b+"GEN_Ks_lxy_status8"]= dir_V0s_GEN_Ks_kinematics_status8.make<TH1F>(b+"GEN_Ks_lxy_status8", b+"GEN_Ks_lxy_status8; lxy(cm)", 1000, 0, 100);
   histos_th2f[b+"GEN_Ks_lxy_pt_status8"]= dir_V0s_GEN_Ks_kinematics_status8.make<TH2F>(b+"GEN_Ks_lxy_pt_status8", b+"GEN_Ks_lxy_pt_status8; lxy(cm); pt(GeV)", 1000, 0, 100, 100, 0, 10);
   histos_th1f[b+"GEN_Ks_vz_status8"]= dir_V0s_GEN_Ks_kinematics_status8.make<TH1F>(b+"GEN_Ks_vz_status8", b+"GEN_Ks_vz_status8; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"GEN_Ks_ndaughters_status8"]= dir_V0s_GEN_Ks_kinematics_status8.make<TH1F>(b+"GEN_Ks_ndaughters_status8", b+"GEN_Ks_ndaughters_status8; #daughters", 10, 0, 10);
   histos_th1f[b+"GEN_Ks_status_status8"]= dir_V0s_GEN_Ks_kinematics_status8.make<TH1F>(b+"GEN_Ks_status_status8", b+"GEN_Ks_status_status8; status", 10, 0, 10);
   histos_th2f[b+"GEN_Ks_status_n_daughters_status8"]= dir_V0s_GEN_Ks_kinematics_status8.make<TH2F>(b+"GEN_Ks_status_n_daughters_status8", b+"GEN_Ks_status_n_daughters_status8; status; n daughters", 10, 0, 10, 10, 0, 10);

   //kinematics of the RECO Ks
   TFileDirectory dir_V0s_RECO_Ks_kinematics = dir_V0s.mkdir("RECO_Ks_kinematics");
   histos_th1f[b+"V0s_Ks_pt"]= dir_V0s_RECO_Ks_kinematics.make<TH1F>(b+"V0s_Ks_pt", b+"V0s_Ks_pt; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_eta"]= dir_V0s_RECO_Ks_kinematics.make<TH1F>(b+"V0s_Ks_eta", b+"V0s_Ks_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_phi"]= dir_V0s_RECO_Ks_kinematics.make<TH1F>(b+"V0s_Ks_phi", b+"V0s_Ks_phi; phi(rad)", 200, -10, 10);
   histos_th1f[b+"V0s_Ks_lxy"]= dir_V0s_RECO_Ks_kinematics.make<TH1F>(b+"V0s_Ks_lxy", b+"V0s_Ks_lxy; lxy(cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_vz"]= dir_V0s_RECO_Ks_kinematics.make<TH1F>(b+"V0s_Ks_vz", b+"V0s_Ks_vz; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_Ks_ndaughters"]= dir_V0s_RECO_Ks_kinematics.make<TH1F>(b+"V0s_Ks_ndaughters", b+"V0s_Ks_ndaughters; #daughters", 10, 0, 10);

   //matched Ks to GEN, with Ks not a daughter of the S
   TFileDirectory dir_V0s_matched_Ks = dir_V0s.mkdir("matched_Ks_with_daughters_in_tracker");
   histos_th1f[b+"V0s_Ks_reconstructed"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed", b+"V0s_Ks_reconstructed; 0 = NO GEN  Ks (with daughters in tk accepatance) particle matching a RECO Ks V0", 2, 0, 2);
   histos_th1f[b+"V0s_Ks_deltaR_Ks_GEN_RECO"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_deltaR_Ks_GEN_RECO", b+"V0s_Ks_deltaR_Ks_GEN_RECO; deltaR(Ks RECO, KS GEN)", 1000, 0, 10);
   histos_th1f[b+"V0s_Ks_deltaR_Ks_GEN_RECO_status1"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_deltaR_Ks_GEN_RECO_status1", b+"V0s_Ks_deltaR_Ks_GEN_RECO_status1; deltaR(Ks RECO, KS GEN)", 1000, 0, 10);
   histos_th1f[b+"V0s_Ks_deltaR_Ks_GEN_RECO_status8"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_deltaR_Ks_GEN_RECO_status8", b+"V0s_Ks_deltaR_Ks_GEN_RECO_status8; deltaR(Ks RECO, KS GEN)", 1000, 0, 10);
   histos_th1f[b+"V0s_Ks_deltaR_Ks_GEN_RECO_daughter_antiS"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_deltaR_Ks_GEN_RECO_daughter_antiS", b+"V0s_Ks_deltaR_Ks_GEN_RECO_daughter_antiS; deltaR(Ks RECO, KS GEN daughter antiS)", 1000, 0, 10);

   histos_th1f[b+"V0s_Ks_reconstructed_GEN_pt"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_pt", b+"V0s_Ks_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_p"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_p", b+"V0s_Ks_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_eta"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_eta", b+"V0s_Ks_reconstructed_GEN_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_phi"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_phi", b+"V0s_Ks_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_vxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_vxy", b+"V0s_Ks_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_lxy", b+"V0s_Ks_reconstructed_GEN_lxy; lxy (beamspot, Ks creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_decay_lxy", b+"V0s_Ks_reconstructed_GEN_decay_lxy; lxy(beamspot, Ks decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_decay_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_decay_vz", b+"V0s_Ks_reconstructed_GEN_decay_vz; vz(beamspot, Ks decay vertex) (cm)", 1000, -400, 400);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_dxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_dxy", b+"V0s_Ks_reconstructed_GEN_dxy; dxy (beamspot, Ks) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_vz", b+"V0s_Ks_reconstructed_GEN_vz; vz(beamspot, Ks creation vertex) (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_status"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_status", b+"V0s_Ks_reconstructed_GEN_status; status", 10, 0, 10);

   //matched Ks to GEN, with Ks a daughter of the S
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_pt"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_pt", b+"V0s_Ks_daughterS_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_p"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_p", b+"V0s_Ks_daughterS_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_eta"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_eta", b+"V0s_Ks_daughterS_reconstructed_GEN_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_phi"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_phi", b+"V0s_Ks_daughterS_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_vxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_vxy", b+"V0s_Ks_daughterS_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_lxy", b+"V0s_Ks_daughterS_reconstructed_GEN_lxy; lxy (beamspot, Ks creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_decay_lxy", b+"V0s_Ks_daughterS_reconstructed_GEN_decay_lxy; lxy(beamspot, Ks decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_decay_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_decay_vz", b+"V0s_Ks_daughterS_reconstructed_GEN_decay_vz; vz(beamspot, Ks decay vertex) (cm)", 1000, -400, 400);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_dxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_dxy", b+"V0s_Ks_daughterS_reconstructed_GEN_dxy; dxy (beamspot, Ks) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_vz", b+"V0s_Ks_daughterS_reconstructed_GEN_vz; vz(beamspot, Ks creation vertex) (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_status"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_status", b+"V0s_Ks_daughterS_reconstructed_GEN_status; status", 10, 0, 10);

   //non matched Ks to GEN, with the Ks not a daughter of the S
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_pt"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_pt", b+"V0s_Ks_non_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_p"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_p", b+"V0s_Ks_non_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_eta"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_eta", b+"V0s_Ks_non_reconstructed_GEN_eta; eta ", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_phi"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_phi", b+"V0s_Ks_non_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_vxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_vxy", b+"V0s_Ks_non_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_lxy", b+"V0s_Ks_non_reconstructed_GEN_lxy; lxy(beamspot, Ks creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_decay_lxy", b+"V0s_Ks_non_reconstructed_GEN_decay_lxy; lxy (beamspot, Ks decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_decay_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_decay_vz", b+"V0s_Ks_non_reconstructed_GEN_decay_vz; vz (beamspot, Ks decay vertex) (cm)", 1000, -400, 400);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_dxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_dxy", b+"V0s_Ks_non_reconstructed_GEN_dxy; dxy(beamspot, Ks) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_vz", b+"V0s_Ks_non_reconstructed_GEN_vz; vz(beamspot, Ks creation vertex) (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_status"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_status", b+"V0s_Ks_non_reconstructed_GEN_status; status", 10, 0, 10);

   //non matched Ks to GEN, with the Ks a daughter of the S
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_pt"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_pt", b+"V0s_Ks_daughterS_non_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_p"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_p", b+"V0s_Ks_daughterS_non_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_eta"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_eta", b+"V0s_Ks_daughterS_non_reconstructed_GEN_eta; eta ", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_phi"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_phi", b+"V0s_Ks_daughterS_non_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_vxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_vxy", b+"V0s_Ks_daughterS_non_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_lxy", b+"V0s_Ks_daughterS_non_reconstructed_GEN_lxy; lxy(beamspot, Ks creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_lxy", b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_lxy; lxy (beamspot, Ks decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_vz", b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_vz; vz (beamspot, Ks decay vertex) (cm)", 1000, -400, 400);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_dxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_dxy", b+"V0s_Ks_daughterS_non_reconstructed_GEN_dxy; dxy(beamspot, Ks) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_vz", b+"V0s_Ks_daughterS_non_reconstructed_GEN_vz; vz(beamspot, Ks creation vertex) (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_status"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_status", b+"V0s_Ks_daughterS_non_reconstructed_GEN_status; status", 10, 0, 10);

   histos_teff[b+"V0_Ks_reconstructed_pt"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_pt",b+"V0_Ks_reconstructed_pt; Ks pt(GeV); Ks reconstruction eff",100,0,10);
   histos_teff[b+"V0_Ks_reconstructed_p"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_p",b+"V0_Ks_reconstructed_p; Ks p(GeV);Ks reconstruction eff",300,0,30);
   histos_teff[b+"V0_Ks_reconstructed_vxy"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_vxy",b+"V0_Ks_reconstructed_vxy; vxy(cm);Ks reconstruction eff",100,0,100);
   histos_teff[b+"V0_Ks_reconstructed_lxy"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_lxy",b+"V0_Ks_reconstructed_lxy; lxy (beamspot, Ks creation vertex)(cm);Ks reconstruction eff",100,0,100);
   histos_teff[b+"V0_Ks_reconstructed_decay_lxy"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_decay_lxy",b+"V0_Ks_reconstructed_decay_lxy; lxy (beamspot, Ks decay vertex)(cm);Ks reconstruction eff",200,0,100);
   histos_teff[b+"V0_Ks_reconstructed_decay_vz"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_decay_vz",b+"V0_Ks_reconstructed_decay_vz; vz (beamspot, Ks decay vertex)(cm);Ks reconstruction eff",100,-400,400);
   histos_teff[b+"V0_Ks_reconstructed_dxy"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_dxy",b+"V0_Ks_reconstructed_dxy; dxy (beamspot, Ks)(cm);Ks reconstruction eff",400,-100,100);
   histos_teff[b+"V0_Ks_reconstructed_vz"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_vz",b+"V0_Ks_reconstructed_vz; vz (beamspot, Ks creation)(cm);Ks reconstruction eff",400,-400,400);
   histos_teff[b+"V0_Ks_reconstructed_eta"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_eta",b+"V0_Ks_reconstructed_eta; eta;Ks reconstruction eff",100,-4,4);
   histos_teff[b+"V0_Ks_reconstructed_phi"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_phi",b+"V0_Ks_reconstructed_phi; phi(rad);Ks reconstruction eff",100,-4,4);
   histos_teff[b+"V0_Ks_reconstructed_status"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_status",b+"V0_Ks_reconstructed_status; status;Ks reconstruction eff",10,0,10);


				
   //Lambda
   //kinematics of the GEN L
   TFileDirectory dir_V0s_GEN_L_kinematics = dir_V0s.mkdir("GEN_L_kinematics");
   histos_th1f[b+"GEN_L_pt"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_pt", b+"GEN_L_pt; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_L_pt_daughter_antiS"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_pt_daughter_antiS", b+"GEN_L_pt_daughter_antiS; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_L_daughter_antiS_pt_of_L_daughters"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_daughter_antiS_pt_of_L_daughters", b+"GEN_L_daughter_antiS_pt_of_L_daughters; pt(GeV)", 100, 0, 10);
   histos_th2f[b+"GEN_L_daughter_antiS_min_max_pt_of_L_daughters"]= dir_V0s_GEN_L_kinematics.make<TH2F>(b+"GEN_L_daughter_antiS_min_max_pt_of_L_daughters", b+"GEN_L_daughter_antiS_min_max_pt_of_L_daughters; min pt daughter L(GeV); max pt daughter L(GeV)", 100, 0, 10,100,0,10);
   histos_th1f[b+"GEN_L_eta"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_eta", b+"GEN_L_eta; eta", 100, -4, 4);
   histos_th1f[b+"GEN_L_phi"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_phi", b+"GEN_L_phi; phi(rad)", 200, -10, 10);
   histos_th1f[b+"GEN_L_lxy"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_lxy", b+"GEN_L_lxy; lxy(cm)", 1000, 0, 100);
   histos_th2f[b+"GEN_L_lxy_pt"]= dir_V0s_GEN_L_kinematics.make<TH2F>(b+"GEN_L_lxy_pt", b+"GEN_L_lxy_pt; lxy(cm); pt(GeV)", 1000, 0, 100, 100, 0, 10);
   histos_th1f[b+"GEN_L_vz"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_vz", b+"GEN_L_vz; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"GEN_L_ndaughters"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_ndaughters", b+"GEN_L_ndaughters; #daughters", 10, 0, 10);
   histos_th1f[b+"GEN_L_status"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_status", b+"GEN_L_status; status", 10, 0, 10);
   histos_th2f[b+"GEN_L_status_n_daughters"]= dir_V0s_GEN_L_kinematics.make<TH2F>(b+"GEN_L_status_n_daughters", b+"GEN_L_status_n_daughters; status; n daughters", 10, 0, 10, 10, 0, 10);

   //kinematics of the GEN with status 2
   TFileDirectory dir_V0s_GEN_L_kinematics_status1 = dir_V0s.mkdir("GEN_L_kinematics_status1");
   histos_th1f[b+"GEN_L_pt_status1"]= dir_V0s_GEN_L_kinematics_status1.make<TH1F>(b+"GEN_L_pt_status1", b+"GEN_L_pt_status1; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_L_eta_status1"]= dir_V0s_GEN_L_kinematics_status1.make<TH1F>(b+"GEN_L_eta_status1", b+"GEN_L_eta_status1; eta", 100, -4, 4);
   histos_th1f[b+"GEN_L_phi_status1"]= dir_V0s_GEN_L_kinematics_status1.make<TH1F>(b+"GEN_L_phi_status1", b+"GEN_L_phi_status1; phi(rad)", 200, -10, 10);
   histos_th1f[b+"GEN_L_lxy_status1"]= dir_V0s_GEN_L_kinematics_status1.make<TH1F>(b+"GEN_L_lxy_status1", b+"GEN_L_lxy_status1; lxy(cm)", 1000, 0, 100);
   histos_th2f[b+"GEN_L_lxy_pt_status1"]= dir_V0s_GEN_L_kinematics_status1.make<TH2F>(b+"GEN_L_lxy_pt_status1", b+"GEN_L_lxy_pt_status1; lxy(cm); pt(GeV)", 1000, 0, 100, 100,0,10);
   histos_th1f[b+"GEN_L_vz_status1"]= dir_V0s_GEN_L_kinematics_status1.make<TH1F>(b+"GEN_L_vz_status1", b+"GEN_L_vz_status1; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"GEN_L_ndaughters_status1"]= dir_V0s_GEN_L_kinematics_status1.make<TH1F>(b+"GEN_L_ndaughters_status1", b+"GEN_L_ndaughters_status1; #daughters", 10, 0, 10);
   histos_th1f[b+"GEN_L_status_status1"]= dir_V0s_GEN_L_kinematics_status1.make<TH1F>(b+"GEN_L_status_status1", b+"GEN_L_status_status1; status", 10, 0, 10);
   histos_th2f[b+"GEN_L_status_n_daughters_status1"]= dir_V0s_GEN_L_kinematics_status1.make<TH2F>(b+"GEN_L_status_n_daughters_status1", b+"GEN_L_status_n_daughters_status1; status; n daughters", 10, 0, 10, 10, 0, 10);


   //kinematics of the GEN with status 8
   TFileDirectory dir_V0s_GEN_L_kinematics_status8 = dir_V0s.mkdir("GEN_L_kinematics_status8");
   histos_th1f[b+"GEN_L_pt_status8"]= dir_V0s_GEN_L_kinematics_status8.make<TH1F>(b+"GEN_L_pt_status8", b+"GEN_L_pt_status8; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_L_eta_status8"]= dir_V0s_GEN_L_kinematics_status8.make<TH1F>(b+"GEN_L_eta_status8", b+"GEN_L_eta_status8; eta", 100, -4, 4);
   histos_th1f[b+"GEN_L_phi_status8"]= dir_V0s_GEN_L_kinematics_status8.make<TH1F>(b+"GEN_L_phi_status8", b+"GEN_L_phi_status8; phi(rad)", 200, -10, 10);
   histos_th1f[b+"GEN_L_lxy_status8"]= dir_V0s_GEN_L_kinematics_status8.make<TH1F>(b+"GEN_L_lxy_status8", b+"GEN_L_lxy_status8; lxy(cm)", 1000, 0, 100);
   histos_th2f[b+"GEN_L_lxy_pt_status8"]= dir_V0s_GEN_L_kinematics_status8.make<TH2F>(b+"GEN_L_lxy_pt_status8", b+"GEN_L_lxy_pt_status8; lxy(cm); pt(GeV)", 1000, 0, 100, 100,0,10);
   histos_th1f[b+"GEN_L_vz_status8"]= dir_V0s_GEN_L_kinematics_status8.make<TH1F>(b+"GEN_L_vz_status8", b+"GEN_L_vz_status8; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"GEN_L_ndaughters_status8"]= dir_V0s_GEN_L_kinematics_status8.make<TH1F>(b+"GEN_L_ndaughters_status8", b+"GEN_L_ndaughters_status8; #daughters", 10, 0, 10);
   histos_th1f[b+"GEN_L_status_status8"]= dir_V0s_GEN_L_kinematics_status8.make<TH1F>(b+"GEN_L_status_status8", b+"GEN_L_status_status8; status", 10, 0, 10);
   histos_th2f[b+"GEN_L_status_n_daughters_status8"]= dir_V0s_GEN_L_kinematics_status8.make<TH2F>(b+"GEN_L_status_n_daughters_status8", b+"GEN_L_status_n_daughters_status8; status; n daughters", 10, 0, 10, 10, 0, 10);

   //kinematics of the RECO L
   TFileDirectory dir_V0s_RECO_L_kinematics = dir_V0s.mkdir("RECO_L_kinematics");
   histos_th1f[b+"V0s_L_pt"]= dir_V0s_RECO_L_kinematics.make<TH1F>(b+"V0s_L_pt", b+"V0s_L_pt; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_eta"]= dir_V0s_RECO_L_kinematics.make<TH1F>(b+"V0s_L_eta", b+"V0s_L_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_L_phi"]= dir_V0s_RECO_L_kinematics.make<TH1F>(b+"V0s_L_phi", b+"V0s_L_phi; phi(rad)", 200, -10, 10);
   histos_th1f[b+"V0s_L_lxy"]= dir_V0s_RECO_L_kinematics.make<TH1F>(b+"V0s_L_lxy", b+"V0s_L_lxy; lxy(cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_vz"]= dir_V0s_RECO_L_kinematics.make<TH1F>(b+"V0s_L_vz", b+"V0s_L_vz; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_L_ndaughters"]= dir_V0s_RECO_L_kinematics.make<TH1F>(b+"V0s_L_ndaughters", b+"V0s_L_ndaughters; #daughters", 10, 0, 10);

   //matched L to GEN
   TFileDirectory dir_V0s_matched_L = dir_V0s.mkdir("matched_L_with_daughters_in_tracker");
   histos_th1f[b+"V0s_L_reconstructed"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed", b+"V0s_L_reconstructed; 0 = NO GEN L (with daughters in tk accepatance)  particle matching a RECO L V0", 2, 0, 2);
   histos_th1f[b+"V0s_L_deltaR_L_GEN_RECO"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_deltaR_L_GEN_RECO", b+"V0s_L_deltaR_L_GEN_RECO; deltaR(L RECO, L GEN)", 1000, 0, 10);
   histos_th1f[b+"V0s_L_deltaR_L_GEN_RECO_status1"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_deltaR_L_GEN_RECO_status1", b+"V0s_L_deltaR_L_GEN_RECO_status1; deltaR(L RECO, L GEN)", 1000, 0, 10);
   histos_th1f[b+"V0s_L_deltaR_L_GEN_RECO_status8"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_deltaR_L_GEN_RECO_status8", b+"V0s_L_deltaR_L_GEN_RECO_status8; deltaR(L RECO, L GEN)", 1000, 0, 10);
   histos_th1f[b+"V0s_L_deltaR_L_GEN_RECO_daughter_antiS"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_deltaR_L_GEN_RECO_daughter_antiS", b+"V0s_L_deltaR_L_GEN_RECO_daughter_antiS; deltaR(L RECO, L GEN)", 1000, 0, 10);
   //matched L to GEN which are not daughters of the S
   histos_th1f[b+"V0s_L_reconstructed_GEN_pt"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_pt", b+"V0s_L_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_reconstructed_GEN_p"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_p", b+"V0s_L_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_reconstructed_GEN_eta"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_eta", b+"V0s_L_reconstructed_GEN_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_L_reconstructed_GEN_phi"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_phi", b+"V0s_L_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_L_reconstructed_GEN_vxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_vxy", b+"V0s_L_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_reconstructed_GEN_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_lxy", b+"V0s_L_reconstructed_GEN_lxy; lxy (beamspot, L creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_decay_lxy", b+"V0s_L_reconstructed_GEN_decay_lxy; lxy (beamspot, L decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_reconstructed_GEN_decay_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_decay_vz", b+"V0s_L_reconstructed_GEN_decay_vz; vz (beamspot, L decay vertex) (cm)", 1000, -400, 400);
   histos_th1f[b+"V0s_L_reconstructed_GEN_dxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_dxy", b+"V0s_L_reconstructed_GEN_dxy; dxy (beamspot, L) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_L_reconstructed_GEN_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_vz", b+"V0s_L_reconstructed_GEN_vz; vz(beamspot, #Lambda creation vertex) (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_L_reconstructed_GEN_status"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_status", b+"V0s_L_reconstructed_GEN_status; status", 10, 0, 10);

   //matched L to GEN which are daughters of the S
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_pt"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_pt", b+"V0s_L_daughterS_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_p"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_p", b+"V0s_L_daughterS_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_eta"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_eta", b+"V0s_L_daughterS_reconstructed_GEN_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_phi"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_phi", b+"V0s_L_daughterS_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_vxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_vxy", b+"V0s_L_daughterS_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_lxy", b+"V0s_L_daughterS_reconstructed_GEN_lxy; lxy (beamspot, L creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_decay_lxy", b+"V0s_L_daughterS_reconstructed_GEN_decay_lxy; lxy (beamspot, L decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_decay_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_decay_vz", b+"V0s_L_daughterS_reconstructed_GEN_decay_vz; vz (beamspot, L decay vertex) (cm)", 1000, -400, 400);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_dxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_dxy", b+"V0s_L_daughterS_reconstructed_GEN_dxy; dxy (beamspot, L) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_vz", b+"V0s_L_daughterS_reconstructed_GEN_vz; vz(beamspot, #Lambda creation vertex) (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_status"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_status", b+"V0s_L_daughterS_reconstructed_GEN_status; status", 10, 0, 10);

   //non-matched L to GEN which are not daughters of the S
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_pt"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_pt", b+"V0s_L_non_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_p"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_p", b+"V0s_L_non_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_eta"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_eta", b+"V0s_L_non_reconstructed_GEN_eta; eta ", 100, -4, 4);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_phi"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_phi", b+"V0s_L_non_reconstructed_GEN_phi; phi(rad) ", 100, -4, 4);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_vxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_vxy", b+"V0s_L_non_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_lxy", b+"V0s_L_non_reconstructed_GEN_lxy; lxy (beamspot, L creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_decay_lxy", b+"V0s_L_non_reconstructed_GEN_decay_lxy; lxy (beamspot, L decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_decay_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_decay_vz", b+"V0s_L_non_reconstructed_GEN_decay_vz; vz (beamspot, L decay vertex) (cm)", 1000, -400, 400);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_dxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_dxy", b+"V0s_L_non_reconstructed_GEN_dxy; dxy (beamspot, L) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_vz", b+"V0s_L_non_reconstructed_GEN_vz; vz(beamspot, #Lambda creation vertex) (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_status"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_status", b+"V0s_L_non_reconstructed_GEN_status; status", 10, 0, 10);

   //non-matched L to GEN which are daughters of the S
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_pt"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_pt", b+"V0s_L_daughterS_non_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_p"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_p", b+"V0s_L_daughterS_non_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_eta"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_eta", b+"V0s_L_daughterS_non_reconstructed_GEN_eta; eta ", 100, -4, 4);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_phi"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_phi", b+"V0s_L_daughterS_non_reconstructed_GEN_phi; phi(rad) ", 100, -4, 4);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_vxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_vxy", b+"V0s_L_daughterS_non_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_lxy", b+"V0s_L_daughterS_non_reconstructed_GEN_lxy; lxy (beamspot, L creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_decay_lxy", b+"V0s_L_daughterS_non_reconstructed_GEN_decay_lxy; lxy (beamspot, L decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_decay_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_decay_vz", b+"V0s_L_daughterS_non_reconstructed_GEN_decay_vz; vz (beamspot, L decay vertex) (cm)", 1000, -400, 400);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_dxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_dxy", b+"V0s_L_daughterS_non_reconstructed_GEN_dxy; dxy (beamspot, L) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_vz", b+"V0s_L_daughterS_non_reconstructed_GEN_vz; vz(beamspot, #Lambda creation vertex) (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_status"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_status", b+"V0s_L_daughterS_non_reconstructed_GEN_status; status", 10, 0, 10);


   histos_teff[b+"V0_L_reconstructed_pt"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_pt",b+"V0_L_reconstructed_pt; #Lambda pt(GeV);#Lambda reconstruction eff",100,0,10);
   histos_teff[b+"V0_L_reconstructed_p"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_p",b+"V0_L_reconstructed_p; #Lambda p(GeV);#Lambda reconstruction eff",300,0,30);
   histos_teff[b+"V0_L_reconstructed_vxy"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_vxy",b+"V0_L_reconstructed_vxy; vxy(cm);#Lambda reconstruction eff",100,0,100);
   histos_teff[b+"V0_L_reconstructed_lxy"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_lxy",b+"V0_L_reconstructed_lxy; lxy (beamspot, L creation vertex)(cm);#Lambda reconstruction eff",100,0,100);
   histos_teff[b+"V0_L_reconstructed_decay_lxy"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_decay_lxy",b+"V0_L_reconstructed_decay_lxy; lxy (beamspot, L decay vertex)(cm);#Lambda reconstruction eff",200,0,100);
   histos_teff[b+"V0_L_reconstructed_decay_vz"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_decay_vz",b+"V0_L_reconstructed_decay_vz; vz (beamspot, L decay vertex)(cm);#Lambda reconstruction eff",100,-400,400);
   histos_teff[b+"V0_L_reconstructed_dxy"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_dxy",b+"V0_L_reconstructed_dxy; dxy (beamspot, L)(cm);#Lambda reconstruction eff",400,-100,100);
   histos_teff[b+"V0_L_reconstructed_vz"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_vz",b+"V0_L_reconstructed_vz; vz (beamspot, #Lambda creation vertex)(cm);#Lambda reconstruction eff",400,-400,400);
   histos_teff[b+"V0_L_reconstructed_eta"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_eta",b+"V0_L_reconstructed_eta; eta;#Lambda reconstruction eff",100,-4,4);
   histos_teff[b+"V0_L_reconstructed_phi"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_phi",b+"V0_L_reconstructed_phi; phi(rad);#Lambda reconstruction eff",100,-4,4);
   histos_teff[b+"V0_L_reconstructed_status"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_status",b+"V0_L_reconstructed_status; status;#Lambda reconstruction eff",10,0,10);

   //For the reconstructed S particles:
    TFileDirectory dir_LambdaKshortVertexFilter = m_fs->mkdir("LambdaKshortVertexFilter");
    TFileDirectory dir_LambdaKshortVertexFilter_S = dir_LambdaKshortVertexFilter.mkdir("S");
    histos_th1f[b+"h_LambdaKshortVertexFilter_S_mass"]= dir_LambdaKshortVertexFilter_S.make<TH1F>(b+"h_LambdaKshortVertexFilter_S_mass", b+"h_LambdaKshortVertexFilter_S_mass; S mass (GeV)", 2000, -100, 100);
    histos_th1f[b+"h_LambdaKshortVertexFilter_Sn_mass"]= dir_LambdaKshortVertexFilter_S.make<TH1F>(b+"h_LambdaKshortVertexFilter_Sn_mass", b+"h_LambdaKshortVertexFilter_Sn_mass; S+n mass (GeV)", 2000, -100, 100);
    histos_th1f[b+"h_LambdaKshortVertexFilter_S_mass_with_displacement_larger_1p2"]= dir_LambdaKshortVertexFilter_S.make<TH1F>(b+"h_LambdaKshortVertexFilter_S_mass_with_displacement_larger_1p2", b+"h_LambdaKshortVertexFilter_S_mass_with_displacement_larger_1p2; S mass (GeV)", 2000, -100, 100);
    histos_th2f[b+"h_LambdaKshortVertexFilter_S_vx_vy"]= dir_LambdaKshortVertexFilter_S.make<TH2F>(b+"h_LambdaKshortVertexFilter_S_vx_vy", b+"h_LambdaKshortVertexFilter_S_vx_vy; S vx (cm); S vy (cm)", 2000, -100, 100, 2000, -100, 100);

}


void Analyzer_SIM_Sexaq::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


//  cout << "------------------------------------------------------------------------------------------------------------------------------------------------------------------"<< endl; 
 
  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  //iEvent.getByToken(m_bsToken, h_bs);

  //primary vertex
  edm::Handle<vector<reco::Vertex>> h_offlinePV;
  iEvent.getByToken(m_offlinePVToken, h_offlinePV);

  //gen particles 
  edm::Handle<vector<reco::GenParticle>> h_genParticles_GEN;
  iEvent.getByToken(m_genParticlesToken_GEN, h_genParticles_GEN);
 
  //SIM GEANT particles
  edm::Handle<vector<reco::GenParticle>> h_genParticles_SIM_GEANT;
  iEvent.getByToken(m_genParticlesToken_SIM_GEANT, h_genParticles_SIM_GEANT);

  //General tracks particles
  edm::Handle<vector<reco::Track>> h_generalTracks;
  iEvent.getByToken(m_generalTracksToken, h_generalTracks);

  //lambdaKshortVertexFilter sexaquark candidates
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);

  //V0 Kshorts
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks;
  iEvent.getByToken(m_V0KsToken, h_V0Ks);

  //V0 Lambdas
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L;
  iEvent.getByToken(m_V0LToken, h_V0L);


  //beamspot
  TVector3 beamspot(0,0,0);
  if(h_bs.isValid()){ 	

  	const reco::BeamSpot* theBeamSpot = h_bs.product();
  	math::XYZPoint referencePos(theBeamSpot->position());
	//beamspot
	double bx_x = h_bs->x0(); 
	double bx_y = h_bs->y0(); 
	double bx_z = h_bs->z0();
	beamspot.SetXYZ(bx_x,bx_y,bx_z);
	
  }
  TVector3 FirstOfflinePV(0.,0.,0.);
  if(h_offlinePV.isValid()){ FirstOfflinePV.SetX(h_offlinePV->at(0).x()); FirstOfflinePV.SetY(h_offlinePV->at(0).y()); FirstOfflinePV.SetZ(h_offlinePV->at(0).z());}

  //to store the GEN Ks and L particles that have a RECO counterpart
  vector<reco::GenParticle> v_gen_Ks;
  vector<reco::GenParticle> v_gen_L;
  //look at all the GEN Ks and Lambda, split up between Ks and Lambda with antiS and no antiS mother, and check the RECO eff of the pi+-, proton. Is there a difference between Ks and Lambda from and not from the antiS???
 if(h_genParticles_SIM_GEANT.isValid() && h_generalTracks.isValid()){
      for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i){//loop all gen particles
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() != 2) continue;
		if(fabs(h_genParticles_SIM_GEANT->at(i).pdgId())!=310&&fabs(h_genParticles_SIM_GEANT->at(i).pdgId())!=3122) continue;
		for(unsigned int d = 0; d < 2; ++d){//loop both daughters of the GEN V0 particles
			const reco::Candidate * V0_daug = h_genParticles_SIM_GEANT->at(i).daughter(d);	
			Double_t V0_daug_phi = V0_daug->phi();
			Double_t V0_daug_eta = V0_daug->eta();
			if(fabs(V0_daug_eta)>2.5)continue;
			if(V0_daug->charge()==0)continue;
			//if(V0_daug->pt()<0.5)continue;
			Double_t V0_daug_vxy = pow(V0_daug->vx()*V0_daug->vx()+V0_daug->vy()+V0_daug->vy(),0.5);
			Double_t Delta_R_V0_daug0_track_min = 999;
			//const reco::Track * matchingTrack = nullptr;
			for(unsigned int t = 0; t < h_generalTracks->size(); ++t){//loop over all the tracks and find the minimal deltaR
				Double_t track_phi = h_generalTracks->at(t).phi();
				Double_t track_eta = h_generalTracks->at(t).eta();
				Double_t Delta_R_V0_daug0_track = pow(pow(V0_daug_phi-track_phi,2)+pow(V0_daug_eta-track_eta,2),0.5);
				if(Delta_R_V0_daug0_track < Delta_R_V0_daug0_track_min) {
					Delta_R_V0_daug0_track_min = Delta_R_V0_daug0_track;
					//matchingTrack = &h_generalTracks->at(t);
				}
				histos_th1f[b+"h_tracks_deltaR"]->Fill(Delta_R_V0_daug0_track);
			 	if(h_genParticles_SIM_GEANT->at(i).numberOfMothers()>0)if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId()==-1020000020)histos_th1f[b+"h_tracks_deltaR_grandDaughtersAntiS"]->Fill(Delta_R_V0_daug0_track);	
			}
			histos_th1f[b+"h_tracks_deltaR_min"]->Fill(Delta_R_V0_daug0_track_min);
			if(h_genParticles_SIM_GEANT->at(i).numberOfMothers()>0)if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId()==-1020000020)histos_th1f[b+"h_tracks_deltaR_min_grandDaughtersAntiS"]->Fill(Delta_R_V0_daug0_track_min);
			if(Delta_R_V0_daug0_track_min<0.1){//matched a track with a daughter of a Ks or Lambda
				if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId()==-1020000020){//matched and the daughter is coming from an antiS
					histos_teff[b+"tracks_reco_eff_daughters_antiS"]->Fill(true,1);						
					histos_teff[b+"tracks_reco_eff_daughters_antiS_pt"]->Fill(true,V0_daug->pt());
					histos_teff[b+"tracks_reco_eff_daughters_antiS_eta"]->Fill(true,V0_daug->eta());
					histos_teff[b+"tracks_reco_eff_daughters_antiS_vxy"]->Fill(true,V0_daug_vxy);						
					cout << "found an antiS V0 daughter track" << endl;
				

				 	//now check if this specific track survives the track cuts in the V0Fitter:
					   // fill vectors of TransientTracks and TrackRefs after applying preselection cuts
					    /*  const reco::Track * tmpTrack = matchingTrack;
					      double ipsigXY = std::abs(tmpTrack->dxy(*theBeamSpot)/tmpTrack->dxyError());
					      double ipsigZ = std::abs(tmpTrack->dz(referencePos)/tmpTrack->dzError());
					      if (tmpTrack->normalizedChi2() < 10 && tmpTrack->numberOfValidHits() >= 7 &&
						  tmpTrack->pt() > 0.35 && ipsigXY > 2 && ipsigZ > -1) {
						 //reco::TrackRef tmpRef(theTrackHandle, std::distance(theTrackCollection->begin(), iTk));
						 //theTrackRefs.push_back(std::move(tmpRef));
						 //reco::TransientTrack tmpTransient(*tmpRef, theMagneticField);
						 //theTransTracks.push_back(std::move(tmpTransient));
						histos_teff[b+"tracks_reco_eff_daughters_antiS_surv_track_cuts_pt"]->Fill(true,V0_daug->pt());
						cout << "found an antiS V0 daughter track which survives the track cuts in the V0Fitter" << endl;
					      }
					      else{
						histos_teff[b+"tracks_reco_eff_daughters_antiS_surv_track_cuts_pt"]->Fill(false,V0_daug->pt());
					      }*/

				}
				else{//mathed, but no daughter of antiS
					histos_teff[b+"tracks_reco_eff_daughters_antiS"]->Fill(true,0);						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt"]->Fill(true,V0_daug->pt());						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_eta"]->Fill(true,V0_daug->eta());						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_vxy"]->Fill(true,V0_daug_vxy);						
					if(V0_daug_vxy<3.5) histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5"]->Fill(true,V0_daug->pt());						
				}

			}
			else{//no matched track found
				if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId()==-1020000020){
					histos_teff[b+"tracks_reco_eff_daughters_antiS"]->Fill(false,1);						
					histos_teff[b+"tracks_reco_eff_daughters_antiS_pt"]->Fill(false,V0_daug->pt());						
					histos_teff[b+"tracks_reco_eff_daughters_antiS_eta"]->Fill(false,V0_daug->eta());						
					histos_teff[b+"tracks_reco_eff_daughters_antiS_vxy"]->Fill(false,V0_daug_vxy);						
				}
				else{
					histos_teff[b+"tracks_reco_eff_daughters_antiS"]->Fill(false,0);						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt"]->Fill(false,V0_daug->pt());						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_eta"]->Fill(false,V0_daug->eta());						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_vxy"]->Fill(false,V0_daug_vxy);						
					if(V0_daug_vxy<3.5) histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5"]->Fill(false,V0_daug->pt());						
				}
			}

		}//loop over daughters of the Lambda or Ks
      }//for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i)
 }//if(h_genParticles_SIM_GEANT.isValid())
 
  //look at the pure GEN particles
 if(h_genParticles_GEN.isValid() )
 {
	 //    cout << "-------------------------------------THE GEN PARTICLES:--------------------------------------------" << endl;
	 for(unsigned int i = 0; i < h_genParticles_GEN->size(); ++i){

		 histos_th1f[b+"h_genParticles_pdgid"]->Fill(h_genParticles_GEN->at(i).pdgId());
		 //	cout << "h_genParticles_GEN->at(i).pdgId()" << h_genParticles_GEN->at(i).pdgId() << endl;
		 //only looking at the anti-S
		 //cout << h_genParticles_GEN->at(i).pdgId() << " " << h_genParticles_GEN->at(i).px() << " " << h_genParticles_GEN->at(i).py() << " " << h_genParticles_GEN->at(i).pz()  << " " << h_genParticles_GEN->at(i).vx()  << " " << h_genParticles_GEN->at(i).vy() << " " << h_genParticles_GEN->at(i).vz() << endl;


		 if(h_genParticles_GEN->at(i).pdgId() == -1020000020){ //Xi: 13324; antiS: -1020000020

			 const reco::GenParticle antiS_gen =  h_genParticles_GEN->at(i);

			 histos_th1f[b+"h_genParticles_antiS_pdgid"]->Fill(antiS_gen.pdgId());
			 histos_th1f[b+"h_genParticles_antiS_pt"]->Fill(antiS_gen.pt());

			 //look at the daughters of the S particle
			 //number of daughters
			 histos_th1f[b+"h_genParticles_antiS_n_daughters"]->Fill(antiS_gen.numberOfDaughters());

			 //try to match a reco particle with the generated s particle
			 Double_t antiS_gen_phi = antiS_gen.phi();
			 Double_t antiS_gen_eta = antiS_gen.eta();
			 if(h_sCands.isValid()) {
				 unsigned int n_sCands = h_sCands->size();
				 for (unsigned int j = 0; j < n_sCands; ++j) {

					 const reco::GenParticle antiSn_reco = h_sCands->at(j);
					 TLorentzVector p4_n(0,0,0,0.939565);
					 TLorentzVector p4_reco_antiSn(antiSn_reco.px(),antiSn_reco.py(),antiSn_reco.pz(),antiSn_reco.energy());
					 TLorentzVector p4_reco_antiS = p4_reco_antiSn-p4_n;



					 Double_t antiS_reco_phi = p4_reco_antiS.Phi();
					 Double_t antiS_reco_eta = p4_reco_antiS.Eta();
					 Double_t delta_phi_antiS_reco_gen = reco::deltaPhi(antiS_reco_phi,antiS_gen_phi); 
					 Double_t delta_eta_antiS_reco_gen = antiS_gen_eta - antiS_reco_eta;
					 Double_t delta_R_antiS_reco_gen = pow(delta_phi_antiS_reco_gen*delta_phi_antiS_reco_gen+delta_eta_antiS_reco_gen*delta_eta_antiS_reco_gen,0.5);

					 if(delta_R_antiS_reco_gen  < 100){
						 cout << "found a GEN and RECO antiS:" << endl;
						 cout << "--coordinates daughters--" << endl;
						 cout << "RECO anti S: daughter Ks vx, vy, vz: " << h_sCands->at(j).daughter(1)->vx() << ", " << h_sCands->at(j).daughter(1)->vy() << "," << h_sCands->at(j).daughter(1)->vz() << endl;
						 cout << "RECO anti S: daughter Ks px, py, pz: " << h_sCands->at(j).daughter(1)->px() << ", " << h_sCands->at(j).daughter(1)->py() << "," << h_sCands->at(j).daughter(1)->pz() << endl;
						 cout << "RECO anti S: daughter L vx, vy, vz: " << h_sCands->at(j).daughter(0)->vx() << ", " << h_sCands->at(j).daughter(0)->vy() << "," << h_sCands->at(j).daughter(0)->vz() << endl;
						 cout << "RECO anti S: daughter L px, py, pz: " << h_sCands->at(j).daughter(0)->px() << ", " << h_sCands->at(j).daughter(0)->py() << "," << h_sCands->at(j).daughter(0)->pz() << endl;
						 cout << "--directions antiS--" << endl;
						 cout << "GEN antiS phi, eta: " << antiS_gen_phi << ", " << antiS_gen_eta << endl;
						 cout << "RECO antiS phi, eta: " << antiS_reco_phi << ", " << antiS_reco_eta << endl;
						 cout << "RECO antiSn phi, eta: " << antiSn_reco.phi() << ", " << antiSn_reco.eta() << endl;
						 cout << "--vertices antiS--" << endl;
						 cout << "GEN antiS vx, vy, vz: " << antiS_gen.vx() << ", " << antiS_gen.vy() << ","  << antiS_gen.vz() << endl;
						 cout << "RECO antiSn vx, vy, vz: " << antiSn_reco.vx() << ", " << antiSn_reco.vy() << ", " << antiSn_reco.vz() << endl;
						 cout << "--momenta antiS--" << endl;
						 cout << "GEN antiS px, py, pz: " << antiS_gen.px() << ", " << antiS_gen.py() << ", " << antiS_gen.pz() << endl;
						 cout << "RECO antiS px, py, pz: " << p4_reco_antiS.Px() << ", " << p4_reco_antiS.Py() << ", " << p4_reco_antiS.Pz()  << endl;
						 cout << "RECO antiSn px, py, pz: " << antiSn_reco.px() << ", " << antiSn_reco.py() << ", "  << antiSn_reco.pz()  << endl;
						 cout << "--energy antiS--" << endl;
						 cout << "GEN antiS E: " << antiS_gen.energy()  << endl;
						 cout << "RECO antiS E: " << p4_reco_antiS.Energy() << endl;
						 cout << "RECO antiSn E: " << antiSn_reco.energy()  << endl;
						 cout << "--masses antiS--" << endl;
						 cout << "GEN antiS M (antiS_gen.mass()): " << antiS_gen.mass()  << endl;
						 cout << "RECO antiS M (p4_reco_antiS.M()): " << p4_reco_antiS.M() << endl;
						 Double_t energy_term = h_sCands->at(j).daughter(0)->energy()+h_sCands->at(j).daughter(1)->energy()-0.939565;
						 TVector3 momentum_term(h_sCands->at(j).daughter(0)->px()+h_sCands->at(j).daughter(1)->px(), h_sCands->at(j).daughter(0)->py()+h_sCands->at(j).daughter(1)->py(), h_sCands->at(j).daughter(0)->pz()+h_sCands->at(j).daughter(1)->pz());
						 cout << "RECO antiS inv M: " << pow(pow(energy_term,2)-momentum_term.Mag2(),0.5) << endl; 
						 cout << "RECO antiSn M (antiSn_reco.mass()): " << antiSn_reco.mass() << endl;
						 Double_t energy_term_no_neutron = h_sCands->at(j).daughter(0)->energy()+h_sCands->at(j).daughter(1)->energy();
						 cout << "RECO antiSn inv M: " << pow(pow(energy_term_no_neutron,2)-momentum_term.Mag2(),0.5) << endl; 
						 cout << "-----------------------------------------------------------------" << endl;
					 }
				 }
			 }

		 }

		 //for the generated KS
		 if(h_genParticles_GEN->at(i).pdgId() == 310){
			 const reco::GenParticle Ks_gen =  h_genParticles_GEN->at(i);
			 Double_t Ks_gen_phi = Ks_gen.phi();
			 Double_t Ks_gen_eta = Ks_gen.eta();
			 if(h_sCands.isValid()) {
				 unsigned int n_sCands = h_sCands->size();
				 for (unsigned int j = 0; j < n_sCands; ++j) {
					 Double_t Ks_reco_phi = h_sCands->at(j).daughter(1)->phi();
					 Double_t Ks_reco_eta = h_sCands->at(j).daughter(1)->eta();
					 Double_t delta_phi_Ks_reco_gen = reco::deltaPhi(Ks_reco_phi,Ks_gen_phi);
					 Double_t delta_eta_Ks_reco_gen = Ks_gen_eta - Ks_reco_eta;
					 Double_t delta_R_Ks_reco_gen = pow(delta_phi_Ks_reco_gen*delta_phi_Ks_reco_gen+delta_eta_Ks_reco_gen*delta_eta_Ks_reco_gen,0.5);
					 if(delta_R_Ks_reco_gen  < 0.1){
						 v_gen_Ks.push_back(Ks_gen);
						 cout << "found a GEN and RECO Ks overlapping" << endl;
						 cout << "GEN Ks vx, vy, vz: " << Ks_gen.vx() << ", " << Ks_gen.vy() << "," << Ks_gen.vz() << endl;
						 cout << "RECO Ks vx, vy, vz: " << h_sCands->at(j).daughter(1)->vx() << ", " << h_sCands->at(j).daughter(1)->vy() << "," << h_sCands->at(j).daughter(1)->vz() << endl;
						 cout << "--directions Ks--" << endl;
						 cout << "GEN Ks phi, eta: " << Ks_gen_phi << ", " << Ks_gen_eta << endl;
						 cout << "RECO Ks phi, eta: " << Ks_reco_phi << ", " << Ks_reco_eta << endl;
						 cout << "--momenta Ks--" << endl;
						 cout << "GEN Ks px, py, pz: " << Ks_gen.px() << ", " << Ks_gen.py() << "," << Ks_gen.pz() << endl;
						 cout << "RECO Ks px, py, pz: " << h_sCands->at(j).daughter(1)->px() << ", " << h_sCands->at(j).daughter(1)->py() << "," << h_sCands->at(j).daughter(1)->pz() << endl;
						 cout << "--Enery Ks--" << endl; 
						 cout << "GEN Ks E: " << Ks_gen.energy()  << endl;
						 cout << "RECO Ks E: " << h_sCands->at(j).daughter(1)->energy() << endl;
						 cout << "--Mass Ks--" << endl;
						 cout << "GEN Ks Mass: " << Ks_gen.mass() << endl;
						 cout << "RECO Ks Mass: " << h_sCands->at(j).daughter(1)->mass() << endl;
						 cout << "-----------------------------------------------------------------" << endl;
					 }
				 }
			 }
		 }

		 //for the generated anti Lambda
		 if(h_genParticles_GEN->at(i).pdgId() == -3122){
			 const reco::GenParticle L_gen =  h_genParticles_GEN->at(i);
			 Double_t L_gen_phi = L_gen.phi();
			 Double_t L_gen_eta = L_gen.eta();
			 if(h_sCands.isValid()) {
				 unsigned int n_sCands = h_sCands->size();
				 for (unsigned int j = 0; j < n_sCands; ++j) {
					 Double_t L_reco_phi = h_sCands->at(j).daughter(0)->phi();
					 Double_t L_reco_eta = h_sCands->at(j).daughter(0)->eta();
					 Double_t delta_phi_L_reco_gen = reco::deltaPhi(L_reco_phi,L_gen_phi);
					 Double_t delta_eta_L_reco_gen = L_gen_eta - L_reco_eta;
					 Double_t delta_R_L_reco_gen = pow(delta_phi_L_reco_gen*delta_phi_L_reco_gen+delta_eta_L_reco_gen*delta_eta_L_reco_gen,0.5);
					 if(delta_R_L_reco_gen  < 0.1){
						 v_gen_L.push_back(L_gen);
						 cout << "found a GEN and RECO L overlapping" << endl;
						 cout << "GEN L vx, vy, vz: " << L_gen.vx() << ", " << L_gen.vy() << "," << L_gen.vz() << endl;
						 cout << "RECO L vx, vy, vz: " << h_sCands->at(j).daughter(0)->vx() << ", " << h_sCands->at(j).daughter(0)->vy() << "," << h_sCands->at(j).daughter(0)->vz() << endl;
						 cout << "--directions L--" << endl;
						 cout << "GEN L phi, eta: " << L_gen_phi << ", " << L_gen_eta << endl;
						 cout << "RECO L phi, eta: " << L_reco_phi << ", " << L_reco_eta << endl;
						 cout << "--momenta L--" << endl;
						 cout << "GEN L px, py, pz: " << L_gen.px() << ", " << L_gen.py() << "," << L_gen.pz() << endl;
						 cout << "RECO L px, py, pz: " << h_sCands->at(j).daughter(0)->px() << ", " << h_sCands->at(j).daughter(0)->py() << "," << h_sCands->at(j).daughter(0)->pz() << endl;
						 cout << "--Energy L--" << endl;
						 cout << "GEN L E: " << L_gen.energy()  << endl;
						 cout << "RECO L E: " << h_sCands->at(j).daughter(0)->energy() << endl;
						 cout << "--Mass L--" << endl;
						 cout << "GEN L Mass: " << L_gen.mass() << endl;
						 cout << "RECO L Mass: " << h_sCands->at(j).daughter(0)->mass() << endl;
						 cout << "-----------------------------------------------------------------" << endl;
					 }
				 }
			 }
		 }



	 }
 }


 //  cout << "-----------------------------------------------------------------" << endl;
 //  cout << "number of properly (only when there is an S reconsturcted in this event and the angular separation of the S daughter Ks with a RECO Ks is small) reconstructed Ks: " << v_gen_Ks.size() << endl;
 //  cout << "number of properly (only when there is an S reconsturcted in this event and the angular separation of the S daughter Ks with a RECO L is small) reconstructed L: " << v_gen_L.size() << endl;
 for(unsigned int i = 0; i < v_gen_Ks.size(); i++){
	 for(unsigned int j = 0; j < v_gen_L.size(); j++){
		 Double_t energy_term = v_gen_Ks[i].energy()+v_gen_L[j].energy()-0.939565;
		 TVector3 momentum_term(v_gen_Ks[i].px()+v_gen_L[j].px(), v_gen_Ks[i].py()+v_gen_L[j].py(), v_gen_Ks[i].pz()+v_gen_L[j].pz());
		 cout << "GEN antiS inv M: " << pow(pow(energy_term,2)-momentum_term.Mag2(),0.5) << endl; 

		 //		 math::XYZVector p_daughters = v_gen_Ks[i].momentum() + v_gen_L[j].momentum();
		 ///		 Double_t p_daughters_size = pow(p_daughters.X()*p_daughters.X()+p_daughters.Y()*p_daughters.Y()+p_daughters.Z()*p_daughters.Z(),0.5);
		 //		 Double_t invMass_S = pow(pow(v_gen_Ks[i].energy()+v_gen_L[j].energy()-0.939565,2)-pow(p_daughters_size,2),0.5);
		 //		 cout << "GEN antiS invMass_S: " << invMass_S << endl;
	 }
 }
 // cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;












 //calculate the invariant mass of any Lambda Kshort combination:
 /*  if(h_genParticles_GEN.isValid() )
     {
     for(unsigned int i = 0; i < h_genParticles_GEN->size(); ++i){
     if(h_genParticles_GEN->at(i).pdgId() == -3122){
     for(unsigned int j = 0; j < h_genParticles_GEN->size(); ++j){
     if(h_genParticles_GEN->at(j).pdgId() == 310){
     math::XYZVector p_daughters = h_genParticles_GEN->at(i).momentum() + h_genParticles_GEN->at(j).momentum();
     Double_t p_daughters_size = pow(p_daughters.X()*p_daughters.X()+p_daughters.Y()*p_daughters.Y()+p_daughters.Z()*p_daughters.Z(),0.5);
     Double_t invMass_S = pow(pow(h_genParticles_GEN->at(i).energy()+h_genParticles_GEN->at(j).energy()-0.939565,2)-pow(p_daughters_size,2),0.5);
     cout << "GEN antiS invMass_S: " << invMass_S << endl;	
     histos_th1f[b+"h_genParticles_antiS_mass_check"]->Fill(invMass_S);
     Double_t invMass_Xi = pow(pow(h_genParticles_GEN->at(i).energy()+h_genParticles_GEN->at(j).energy(),2)-pow(p_daughters_size,2),0.5);
 //cout << "GEN Xi inv mass: " << invMass_Xi << endl;
 histos_th1f[b+"h_genParticles_Xi_mass_check"]->Fill(invMass_Xi);

 }
 }
 }
 }
 }
 */

 /*
    if(h_genParticles_GEN.isValid() && h_sCands.isValid())
    {

    for(unsigned int i = 0; i < h_genParticles_GEN->size(); ++i){
    if(h_genParticles_GEN->at(i).pdgId() == -3122){

    for(unsigned int j = 0; j < h_genParticles_GEN->size(); ++j){
    if(h_genParticles_GEN->at(j).pdgId() == 310){


    for(unsigned int k = 0; k < h_sCands->size(); ++k){

    const reco::GenParticle GEN_Lambda = h_genParticles_GEN->at(i);
    const reco::GenParticle GEN_Ks = h_genParticles_GEN->at(j);

    Double_t Delta_phi_Lambda = reco::deltaPhi(h_sCands->at(k).daughter(0)->phi(), GEN_Lambda.phi());
    Double_t Delta_eta_Lambda = h_sCands->at(k).daughter(0)->eta() -  GEN_Lambda.eta();
    Double_t Delta_R_Lambda = pow(Delta_phi_Lambda*Delta_phi_Lambda+Delta_eta_Lambda*Delta_eta_Lambda,0.5);

    Double_t Delta_phi_Ks = reco::deltaPhi(h_sCands->at(k).daughter(1)->phi(), GEN_Ks.phi());
    Double_t Delta_eta_Ks = h_sCands->at(k).daughter(1)->eta() -  GEN_Ks.eta();
    Double_t Delta_R_Ks = pow(Delta_phi_Ks*Delta_phi_Ks+Delta_eta_Ks*Delta_eta_Ks,0.5);


    if(Delta_R_Lambda > 0.1 || Delta_R_Ks > 0.1) continue;
    cout << "-------Angular separation daughters----------"<<endl;
    cout << "Delta_R_Lambda " << Delta_R_Lambda << endl;
    cout << "Delta_R_Ks " << Delta_R_Ks << endl; 

    cout << "-----matched Ks--------" << endl;
    cout << "GEN KS phi, eta " << GEN_Ks.phi() << " " << GEN_Ks.eta() << endl;	 
    cout << "RECO KS phi, eta " << h_sCands->at(k).daughter(1)->phi() << " " << h_sCands->at(k).daughter(1)->eta() << endl;

    cout << "-----matched Lambda--------" << endl;
    cout << "GEN Lambda phi, eta " << GEN_Lambda.phi() << " " << GEN_Lambda.eta() << endl;	
    cout << "RECO Lambda phi, eta " << h_sCands->at(k).daughter(0)->phi() << " " << h_sCands->at(k).daughter(0)->eta() << endl;

    const reco::Candidate * GEN_antiS = GEN_Ks.mother(0);
    cout << "-----matched S mother Ks--------" << endl;
    cout << "GEN antiS mother Ks phi, eta " << GEN_antiS->phi() << " " << GEN_antiS->eta() << endl;	
    cout << "RECO antiS mother Ks phi, eta " << h_sCands->at(k).phi() << " " << h_sCands->at(k).eta() << endl;

    const reco::Candidate * GEN_antiS2 = GEN_Lambda.mother(0);
    cout << "-----matched S2 mother L--------" << endl;
    cout << "GEN antiS mother L phi, eta " << GEN_antiS2->phi() << " " << GEN_antiS2->eta() << endl;	 
    cout << "RECO antiS mother L phi, eta " << h_sCands->at(k).phi() << " " << h_sCands->at(k).eta() << endl;

    cout << "-----Invariant masses GEN--------" << endl;
    math::XYZVector p_daughters_GEN = h_genParticles_GEN->at(i).momentum() + h_genParticles_GEN->at(j).momentum();
    Double_t p_daughters_size_GEN = pow(p_daughters_GEN.X()*p_daughters_GEN.X()+p_daughters_GEN.Y()*p_daughters_GEN.Y()+p_daughters_GEN.Z()*p_daughters_GEN.Z(),0.5);
    Double_t invMass_S_GEN = pow(pow(h_genParticles_GEN->at(i).energy()+h_genParticles_GEN->at(j).energy()-0.939565,2)-pow(p_daughters_size_GEN,2),0.5);
    cout << "GEN antiS invMass_S: " << invMass_S_GEN << endl;	
    histos_th1f[b+"h_genParticles_antiS_angular_matching_mass_check_GEN"]->Fill(invMass_S_GEN);


    Double_t invMass_Xi_GEN = pow(pow(h_genParticles_GEN->at(i).energy()+h_genParticles_GEN->at(j).energy(),2)-pow(p_daughters_size_GEN,2),0.5);
    cout << "GEN Xi inv mass: " << invMass_Xi_GEN << endl;
    histos_th1f[b+"h_genParticles_Xi_mass_angular_matching_check_GEN"]->Fill(invMass_Xi_GEN);

    cout << "-----Invariant masses RECO--------" << endl;
    math::XYZVector p_daughters_RECO = h_sCands->at(k).daughter(0)->momentum() + h_sCands->at(k).daughter(1)->momentum();
    Double_t p_daughters_size_RECO = pow(p_daughters_RECO.X()*p_daughters_RECO.X()+p_daughters_RECO.Y()*p_daughters_RECO.Y()+p_daughters_RECO.Z()*p_daughters_RECO.Z(),0.5);
    Double_t invMass_S_RECO = pow(pow(h_sCands->at(k).daughter(0)->energy()+h_sCands->at(k).daughter(1)->energy()-0.939565,2)-pow(p_daughters_size_RECO,2),0.5);
    cout << "RECO antiS invMass_S: " << invMass_S_RECO << endl;	
    histos_th1f[b+"h_genParticles_antiS_angular_matching_mass_check_RECO"]->Fill(invMass_S_RECO);

    Double_t invMass_Xi_RECO = pow(pow(h_sCands->at(k).daughter(0)->energy()+h_sCands->at(k).daughter(1)->energy(),2)-pow(p_daughters_size_RECO,2),0.5); 
    cout << "RECO Xi inv mass: " << invMass_Xi_RECO << endl;
    histos_th1f[b+"h_genParticles_Xi_mass_angular_matching_check_RECO"]->Fill(invMass_Xi_RECO);
    }
}					
}
}
}
}
*/

//print all the h_genParticles_SIM_GEANT particles
/*  if(h_genParticles_GEN.isValid() )
    {
    cout << "------------------------------------------------THE GEN PARTICLES:----------------------------------------" << endl;
     for(unsigned int i = 0; i < h_genParticles_GEN->size(); ++i){
	cout << h_genParticles_GEN->at(i).status() << " " << h_genParticles_GEN->at(i).pdgId() << " " << h_genParticles_GEN->at(i).px() << " " << h_genParticles_GEN->at(i).py() << " " << h_genParticles_GEN->at(i).pz()  << " " << h_genParticles_GEN->at(i).vx()  << " " << h_genParticles_GEN->at(i).vy() << " " << h_genParticles_GEN->at(i).vz() << endl;
     }
  }

  //print all the h_genParticles_SIM_GEANT particles
  if(h_genParticles_SIM_GEANT.isValid() )
  {
     cout << "------------------------------------------------THE GEN+GEANT PARTICLES:----------------------------------------" << endl;
     for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i){
	cout << h_genParticles_SIM_GEANT->at(i).status() << " " << h_genParticles_SIM_GEANT->at(i).pdgId() << " " << h_genParticles_SIM_GEANT->at(i).px() << " " << h_genParticles_SIM_GEANT->at(i).py() << " " << h_genParticles_SIM_GEANT->at(i).pz()  << " " << h_genParticles_SIM_GEANT->at(i).vx()  << " " << h_genParticles_SIM_GEANT->at(i).vy() << " " << h_genParticles_SIM_GEANT->at(i).vz() << endl;
     }
  }
*/
/*  if(h_genParticles_SIM_GEANT.isValid() )
  {
     cout << "------------------------------------------------THE GEN+GEANT PARTICLES:----------------------------------------" << endl;
     for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i){

	if(h_genParticles_SIM_GEANT->at(i).pdgId() == -1020000020 || h_genParticles_SIM_GEANT->at(i).pdgId() == 310  || h_genParticles_SIM_GEANT->at(i).pdgId() == -3122){

		cout << h_genParticles_SIM_GEANT->at(i).pdgId() << " " << h_genParticles_SIM_GEANT->at(i).px() << " " << h_genParticles_SIM_GEANT->at(i).py() << " " << h_genParticles_SIM_GEANT->at(i).pz()  << " " << h_genParticles_SIM_GEANT->at(i).vx()  << " " << h_genParticles_SIM_GEANT->at(i).vy() << " " << h_genParticles_SIM_GEANT->at(i).vz() << endl;

	}
	
	if(h_genParticles_SIM_GEANT->at(i).pdgId() == -3122){
		Double_t vx_L = h_genParticles_SIM_GEANT->at(i).vx();
		Double_t vy_L = h_genParticles_SIM_GEANT->at(i).vy();
		Double_t vz_L = h_genParticles_SIM_GEANT->at(i).vz();
		for(unsigned int j = 0; j < h_genParticles_SIM_GEANT->size(); ++j){
			if(h_genParticles_SIM_GEANT->at(j).pdgId() == 310){
				Double_t vx_Ks = h_genParticles_SIM_GEANT->at(j).vx();
				Double_t vy_Ks = h_genParticles_SIM_GEANT->at(j).vy();
				Double_t vz_Ks = h_genParticles_SIM_GEANT->at(j).vz();
				if(vx_L == vx_Ks && vy_L == vy_Ks && vz_L == vz_Ks){
					cout << "Found a Ks and Lambda with same point of origin" << endl;
					cout << h_genParticles_SIM_GEANT->at(i).pdgId() << " " << h_genParticles_SIM_GEANT->at(i).px() << " " << h_genParticles_SIM_GEANT->at(i).py() << " " << h_genParticles_SIM_GEANT->at(i).pz()  << " " << h_genParticles_SIM_GEANT->at(i).vx()  << " " << h_genParticles_SIM_GEANT->at(i).vy() << " " << h_genParticles_SIM_GEANT->at(i).vz() << endl;
					cout << h_genParticles_SIM_GEANT->at(j).pdgId() << " " << h_genParticles_SIM_GEANT->at(j).px() << " " << h_genParticles_SIM_GEANT->at(j).py() << " " << h_genParticles_SIM_GEANT->at(j).pz()  << " " << h_genParticles_SIM_GEANT->at(j).vx()  << " " << h_genParticles_SIM_GEANT->at(j).vy() << " " << h_genParticles_SIM_GEANT->at(j).vz() << endl;
					
					 Double_t energy_term = h_genParticles_SIM_GEANT->at(i).energy()+h_genParticles_SIM_GEANT->at(j).energy()-0.939565;
					 TVector3 momentum_term(h_genParticles_SIM_GEANT->at(i).px()+h_genParticles_SIM_GEANT->at(j).px(), h_genParticles_SIM_GEANT->at(i).py()+h_genParticles_SIM_GEANT->at(j).py(), h_genParticles_SIM_GEANT->at(i).pz()+h_genParticles_SIM_GEANT->at(j).pz());
					 Double_t invMassS = pow(pow(energy_term,2)-momentum_term.Mag2(),0.5);
					 cout << "inv mass of S calculated from the Ks and Lambda with same point of origin "<< invMassS << endl;
				}
			}
		}
	}

        if(h_genParticles_SIM_GEANT->at(i).pdgId() == -1020000020 && h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2){
		 cout << "antiS with 2 daughters found in the h_genParticles_SIM_GEANT" << endl;
		 Double_t energy_term = h_genParticles_SIM_GEANT->at(i).daughter(0)->energy()+h_genParticles_SIM_GEANT->at(i).daughter(1)->energy()-0.939565;
		 TVector3 momentum_term(h_genParticles_SIM_GEANT->at(i).daughter(0)->px()+h_genParticles_SIM_GEANT->at(i).daughter(1)->px(), h_genParticles_SIM_GEANT->at(i).daughter(0)->py()+h_genParticles_SIM_GEANT->at(i).daughter(1)->py(), h_genParticles_SIM_GEANT->at(i).daughter(0)->pz()+h_genParticles_SIM_GEANT->at(i).daughter(1)->pz());
		 Double_t invMassS = pow(pow(energy_term,2)-momentum_term.Mag2(),0.5);
//		 cout << "SIM antiS invMass from h_genParticles_SIM_GEANT: " << invMassS << endl;
		 if(invMassS > 1.81 || invMassS < 1.79) cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!The invMass calculated from the h_genParticles_SIM_GEANT particles is not within the expected range" << endl;
	}
     }

  }
*/
 

/*
  if(h_genParticles_GEN.isValid() )
  {
     cout << "------------------------------------------------THE GEN PARTICLES:----------------------------------------" << endl;
     for(unsigned int i = 0; i < h_genParticles_GEN->size(); ++i){

	if(h_genParticles_GEN->at(i).pdgId() == -1020000020 || h_genParticles_GEN->at(i).pdgId() == 310  || h_genParticles_GEN->at(i).pdgId() == -3122){

		cout << h_genParticles_GEN->at(i).pdgId() << " " << h_genParticles_GEN->at(i).px() << " " << h_genParticles_GEN->at(i).py() << " " << h_genParticles_GEN->at(i).pz()  << " " << h_genParticles_GEN->at(i).vx()  << " " << h_genParticles_GEN->at(i).vy() << " " << h_genParticles_GEN->at(i).vz() << endl;

	}
	
	if(h_genParticles_GEN->at(i).pdgId() == -3122){
		Double_t vx_L = h_genParticles_GEN->at(i).vx();
		Double_t vy_L = h_genParticles_GEN->at(i).vy();
		Double_t vz_L = h_genParticles_GEN->at(i).vz();
		for(unsigned int j = 0; j < h_genParticles_GEN->size(); ++j){
			if(h_genParticles_GEN->at(j).pdgId() == 310){
				Double_t vx_Ks = h_genParticles_GEN->at(j).vx();
				Double_t vy_Ks = h_genParticles_GEN->at(j).vy();
				Double_t vz_Ks = h_genParticles_GEN->at(j).vz();
				if(vx_L == vx_Ks && vy_L == vy_Ks && vz_L == vz_Ks){
					cout << "Found a Ks and Lambda with same point of origin" << endl;
					cout << h_genParticles_GEN->at(i).pdgId() << " " << h_genParticles_GEN->at(i).px() << " " << h_genParticles_GEN->at(i).py() << " " << h_genParticles_GEN->at(i).pz()  << " " << h_genParticles_GEN->at(i).vx()  << " " << h_genParticles_GEN->at(i).vy() << " " << h_genParticles_GEN->at(i).vz() << endl;
					cout << h_genParticles_GEN->at(j).pdgId() << " " << h_genParticles_GEN->at(j).px() << " " << h_genParticles_GEN->at(j).py() << " " << h_genParticles_GEN->at(j).pz()  << " " << h_genParticles_GEN->at(j).vx()  << " " << h_genParticles_GEN->at(j).vy() << " " << h_genParticles_GEN->at(j).vz() << endl;
					
					 Double_t energy_term = h_genParticles_GEN->at(i).energy()+h_genParticles_GEN->at(j).energy()-0.939565;
					 TVector3 momentum_term(h_genParticles_GEN->at(i).px()+h_genParticles_GEN->at(j).px(), h_genParticles_GEN->at(i).py()+h_genParticles_GEN->at(j).py(), h_genParticles_GEN->at(i).pz()+h_genParticles_GEN->at(j).pz());
					 Double_t invMassS = pow(pow(energy_term,2)-momentum_term.Mag2(),0.5);
					 cout << "inv mass of S calculated from the Ks and Lambda with same point of origin "<< invMassS << endl;
				}
			}
		}
	}

        if(h_genParticles_GEN->at(i).pdgId() == -1020000020 && h_genParticles_GEN->at(i).numberOfDaughters() == 2){
		 cout << "antiS with 2 daughters found in the h_genParticles_GEN" << endl;
		 Double_t energy_term = h_genParticles_GEN->at(i).daughter(0)->energy()+h_genParticles_GEN->at(i).daughter(1)->energy()-0.939565;
		 TVector3 momentum_term(h_genParticles_GEN->at(i).daughter(0)->px()+h_genParticles_GEN->at(i).daughter(1)->px(), h_genParticles_GEN->at(i).daughter(0)->py()+h_genParticles_GEN->at(i).daughter(1)->py(), h_genParticles_GEN->at(i).daughter(0)->pz()+h_genParticles_GEN->at(i).daughter(1)->pz());
		 Double_t invMassS = pow(pow(energy_term,2)-momentum_term.Mag2(),0.5);
//		 cout << "SIM antiS invMass from h_genParticles_GEN: " << invMassS << endl;
		 if(invMassS > 1.81 || invMassS < 1.79) cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!The invMass calculated from the h_genParticles_GEN particles is not within the expected range" << endl;
	}
     }

  }

*/



  bool genAntiSThisEventGoodDaughters = false;

  if(h_genParticles_SIM_GEANT.isValid()){
      for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i){
	if(h_genParticles_SIM_GEANT->at(i).pdgId()==-1020000020){
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() != 2) continue;
		Double_t vx_L = h_genParticles_SIM_GEANT->at(i).daughter(0)->vx();
		Double_t vy_L = h_genParticles_SIM_GEANT->at(i).daughter(0)->vy();
		Double_t vxy_L = pow(vx_L*vx_L+vy_L*vy_L,0.5);
		Double_t vz_L = h_genParticles_SIM_GEANT->at(i).daughter(0)->vz();
		Double_t vx_Ks = h_genParticles_SIM_GEANT->at(i).daughter(1)->vx();
		Double_t vy_Ks = h_genParticles_SIM_GEANT->at(i).daughter(1)->vy();
		Double_t vz_Ks = h_genParticles_SIM_GEANT->at(i).daughter(1)->vz();
		Double_t vxy_Ks = pow(vx_Ks*vx_Ks+vy_Ks*vy_Ks,0.5);

		if(vxy_L < 3 && vxy_Ks < 3 && fabs(vz_L) < 10 && fabs(vz_Ks) < 10) genAntiSThisEventGoodDaughters = true; 
		if(genAntiSThisEventGoodDaughters){
			
			cout<<  "--------------------------genParticlesPlusGEANT all the antiS and it's daughters--------------------------" << endl;
			cout << "genParticlesPlusGEANT antiS:" << endl;
			cout << "pdgId " << h_genParticles_SIM_GEANT->at(i).pdgId() << endl;
			cout << "status " << h_genParticles_SIM_GEANT->at(i).status() << endl;
			cout << "phi, eta " <<  h_genParticles_SIM_GEANT->at(i).phi() << " " <<  h_genParticles_SIM_GEANT->at(i).eta()   << endl;
			cout << "px, py, pz " <<  h_genParticles_SIM_GEANT->at(i).px() << " " <<  h_genParticles_SIM_GEANT->at(i).py() << " " <<  h_genParticles_SIM_GEANT->at(i).pz()   << endl;
			cout << "vx, vy, vz " << h_genParticles_SIM_GEANT->at(i).vx() << " " <<  h_genParticles_SIM_GEANT->at(i).vy() << " " <<  h_genParticles_SIM_GEANT->at(i).vz()   << endl;

			cout << "genParticlesPlusGEANT Daughters of the antiS: " << endl;
			cout << "pdgId " << h_genParticles_SIM_GEANT->at(i).daughter(0)->pdgId() << endl;
			cout << "status " << h_genParticles_SIM_GEANT->at(i).daughter(0)->status() << endl;
			cout << "phi, eta " <<  h_genParticles_SIM_GEANT->at(i).daughter(0)->phi() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(0)->eta()   << endl;
			cout << "px, py, pz " <<  h_genParticles_SIM_GEANT->at(i).daughter(0)->px() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(0)->py() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(0)->pz()   << endl;
			cout << "vx, vy, vz " << h_genParticles_SIM_GEANT->at(i).daughter(0)->vx() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(0)->vy() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(0)->vz()   << endl;
			cout << "#daughters of the daughter " << h_genParticles_SIM_GEANT->at(i).daughter(0)->numberOfDaughters() << endl;
			cout << "pdgId " << h_genParticles_SIM_GEANT->at(i).daughter(1)->pdgId() << endl;
			cout << "status " << h_genParticles_SIM_GEANT->at(i).daughter(1)->status() << endl;
			cout << "phi, eta " <<  h_genParticles_SIM_GEANT->at(i).daughter(1)->phi() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(1)->eta()   << endl;
			cout << "px, py, pz " <<  h_genParticles_SIM_GEANT->at(i).daughter(1)->px() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(1)->py() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(1)->pz()   << endl;
			cout << "vx, vy, vz " << h_genParticles_SIM_GEANT->at(i).daughter(1)->vx() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(1)->vy() << " " <<  h_genParticles_SIM_GEANT->at(i).daughter(1)->vz()   << endl;
			cout << "#daughters of the daughter " << h_genParticles_SIM_GEANT->at(i).daughter(1)->numberOfDaughters() << endl;


		}
	}
      }
  }

  //just print all the h_genParticles_SIM_GEANT particles
  if(h_genParticles_SIM_GEANT.isValid() && genAntiSThisEventGoodDaughters){
	cout<<  "--------------------------genParticlesPlusGEANT all  Ks and Lambdas--------------------------" << endl;
      for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i){
		//if(h_genParticles_SIM_GEANT->at(i).status()!=8)continue;
		if(fabs(h_genParticles_SIM_GEANT->at(i).pdgId())!=310 && fabs(h_genParticles_SIM_GEANT->at(i).pdgId())!=3122)continue;
		cout << "genParticlesPlusGEANT particle #" << i << endl;
		cout << "pdgId " << h_genParticles_SIM_GEANT->at(i).pdgId() << endl;
		cout << "status " << h_genParticles_SIM_GEANT->at(i).status() << endl;
		cout << "phi, eta " <<  h_genParticles_SIM_GEANT->at(i).phi() << " " <<  h_genParticles_SIM_GEANT->at(i).eta()   << endl;
		cout << "px, py, pz " <<  h_genParticles_SIM_GEANT->at(i).px() << " " <<  h_genParticles_SIM_GEANT->at(i).py() << " " <<  h_genParticles_SIM_GEANT->at(i).pz()   << endl;
		cout << "vx, vy, vz " <<  h_genParticles_SIM_GEANT->at(i).vx() << " " <<  h_genParticles_SIM_GEANT->at(i).vy() << " " <<  h_genParticles_SIM_GEANT->at(i).vz()   << endl;
	}
  }


  //just print all the GEN particles
  if(h_genParticles_GEN.isValid() && genAntiSThisEventGoodDaughters){
	cout<<  "--------------------------genParticles_GEN all  Ks and Lambdas--------------------------" << endl;
      for(unsigned int i = 0; i < h_genParticles_GEN->size(); ++i){
		if(fabs(h_genParticles_GEN->at(i).pdgId())!=310 && fabs(h_genParticles_GEN->at(i).pdgId())!=3122)continue;
		cout << "genParticles particle #" << i << endl;
		cout << "pdgId " << h_genParticles_GEN->at(i).pdgId() << endl;
		cout << "status " << h_genParticles_GEN->at(i).status() << endl;
		cout << "phi, eta " <<  h_genParticles_GEN->at(i).phi() << " " <<  h_genParticles_GEN->at(i).eta()   << endl;
		cout << "px, py, pz " <<  h_genParticles_GEN->at(i).px() << " " <<  h_genParticles_GEN->at(i).py() << " " <<  h_genParticles_GEN->at(i).pz()   << endl;
		cout << "vx, vy, vz " <<  h_genParticles_GEN->at(i).vx() << " " <<  h_genParticles_GEN->at(i).vy() << " " <<  h_genParticles_GEN->at(i).vz()   << endl;
	}
  }
	
  if(!h_genParticles_SIM_GEANT.isValid()) cout << "h_genParticles_SIM_GEANT not valid" << endl;

  //the kinematics of the GenParticlesPlusGeant Ks
  if(h_genParticles_SIM_GEANT.isValid())
  {
    for(unsigned int k = 0; k < h_genParticles_SIM_GEANT->size(); k++){
        if(fabs(h_genParticles_SIM_GEANT->at(k).pdgId())==310){
			histos_th1f[b+"GEN_Ks_pt"]->Fill(h_genParticles_SIM_GEANT->at(k).pt());				
			if(h_genParticles_SIM_GEANT->at(k).numberOfMothers()>0)if(h_genParticles_SIM_GEANT->at(k).mother()->pdgId()==-1020000020){
				histos_th1f[b+"GEN_Ks_pt_daughter_antiS"]->Fill(h_genParticles_SIM_GEANT->at(k).pt());				
				if(h_genParticles_SIM_GEANT->at(k).numberOfDaughters() == 2){
					Double_t pt_Ks_daughter1 = h_genParticles_SIM_GEANT->at(k).daughter(0)->pt();
					Double_t pt_Ks_daughter2 = h_genParticles_SIM_GEANT->at(k).daughter(1)->pt();
					Double_t min_pt_Ks_daughter = std::min(pt_Ks_daughter1,pt_Ks_daughter2);
					Double_t max_pt_Ks_daughter = std::max(pt_Ks_daughter1,pt_Ks_daughter2);
                                        histos_th1f[b+"GEN_Ks_daughter_antiS_pt_of_Ks_daughters"]->Fill(pt_Ks_daughter1);
                                        histos_th1f[b+"GEN_Ks_daughter_antiS_pt_of_Ks_daughters"]->Fill(pt_Ks_daughter2);
                                        histos_th2f[b+"GEN_Ks_daughter_antiS_min_max_pt_of_Ks_daughters"]->Fill(min_pt_Ks_daughter,max_pt_Ks_daughter);
                                }
			}
			histos_th1f[b+"GEN_Ks_eta"]->Fill(h_genParticles_SIM_GEANT->at(k).eta());				
			histos_th1f[b+"GEN_Ks_phi"]->Fill(h_genParticles_SIM_GEANT->at(k).phi());				
			TVector3 xyz_Ks(h_genParticles_SIM_GEANT->at(k).vx(),h_genParticles_SIM_GEANT->at(k).vy(),h_genParticles_SIM_GEANT->at(k).vz());
			Double_t lxy_Ks = lxy(beamspot,xyz_Ks);
			histos_th1f[b+"GEN_Ks_lxy"]->Fill(lxy_Ks);
			histos_th2f[b+"GEN_Ks_lxy_pt"]->Fill(lxy_Ks,h_genParticles_SIM_GEANT->at(k).pt());
			histos_th1f[b+"GEN_Ks_vz"]->Fill(h_genParticles_SIM_GEANT->at(k).vz());
			histos_th1f[b+"GEN_Ks_ndaughters"]->Fill(h_genParticles_SIM_GEANT->at(k).numberOfDaughters());
			histos_th1f[b+"GEN_Ks_status"]->Fill(h_genParticles_SIM_GEANT->at(k).status());
			histos_th2f[b+"GEN_Ks_status_n_daughters"]->Fill(h_genParticles_SIM_GEANT->at(k).status(),h_genParticles_SIM_GEANT->at(k).numberOfDaughters());
		
	}
     }
  }
  //the kinematics of the GenParticlesPlusGeant Ks with status 2
  if(h_genParticles_SIM_GEANT.isValid())
  {
    for(unsigned int k = 0; k < h_genParticles_SIM_GEANT->size(); k++){
        if(fabs(h_genParticles_SIM_GEANT->at(k).pdgId())==310 && h_genParticles_SIM_GEANT->at(k).status()==1){
			histos_th1f[b+"GEN_Ks_pt_status1"]->Fill(h_genParticles_SIM_GEANT->at(k).pt());				
			histos_th1f[b+"GEN_Ks_eta_status1"]->Fill(h_genParticles_SIM_GEANT->at(k).eta());				
			histos_th1f[b+"GEN_Ks_phi_status1"]->Fill(h_genParticles_SIM_GEANT->at(k).phi());				
			TVector3 xyz_Ks(h_genParticles_SIM_GEANT->at(k).vx(),h_genParticles_SIM_GEANT->at(k).vy(),h_genParticles_SIM_GEANT->at(k).vz());
			Double_t lxy_Ks = lxy(beamspot,xyz_Ks);
			histos_th1f[b+"GEN_Ks_lxy_status1"]->Fill(lxy_Ks);
			histos_th2f[b+"GEN_Ks_lxy_pt_status1"]->Fill(lxy_Ks, h_genParticles_SIM_GEANT->at(k).pt());
			histos_th1f[b+"GEN_Ks_vz_status1"]->Fill(h_genParticles_SIM_GEANT->at(k).vz());
			histos_th1f[b+"GEN_Ks_ndaughters_status1"]->Fill(h_genParticles_SIM_GEANT->at(k).numberOfDaughters());
			histos_th1f[b+"GEN_Ks_status_status1"]->Fill(h_genParticles_SIM_GEANT->at(k).status());
			histos_th2f[b+"GEN_Ks_status_n_daughters_status1"]->Fill(h_genParticles_SIM_GEANT->at(k).status(),h_genParticles_SIM_GEANT->at(k).numberOfDaughters());
		
	}
     }
  }

  //the kinematics of the GenParticlesPlusGeant Ks with status 8
  if(h_genParticles_SIM_GEANT.isValid())
  {
    for(unsigned int k = 0; k < h_genParticles_SIM_GEANT->size(); k++){
        if(fabs(h_genParticles_SIM_GEANT->at(k).pdgId())==310 && h_genParticles_SIM_GEANT->at(k).status()==8){
			histos_th1f[b+"GEN_Ks_pt_status8"]->Fill(h_genParticles_SIM_GEANT->at(k).pt());				
			histos_th1f[b+"GEN_Ks_eta_status8"]->Fill(h_genParticles_SIM_GEANT->at(k).eta());				
			histos_th1f[b+"GEN_Ks_phi_status8"]->Fill(h_genParticles_SIM_GEANT->at(k).phi());				
			TVector3 xyz_Ks(h_genParticles_SIM_GEANT->at(k).vx(),h_genParticles_SIM_GEANT->at(k).vy(),h_genParticles_SIM_GEANT->at(k).vz());
			Double_t lxy_Ks = lxy(beamspot,xyz_Ks);
			histos_th1f[b+"GEN_Ks_lxy_status8"]->Fill(lxy_Ks);
			histos_th2f[b+"GEN_Ks_lxy_pt_status8"]->Fill(lxy_Ks,h_genParticles_SIM_GEANT->at(k).pt());
			histos_th1f[b+"GEN_Ks_vz_status8"]->Fill(h_genParticles_SIM_GEANT->at(k).vz());
			histos_th1f[b+"GEN_Ks_ndaughters_status8"]->Fill(h_genParticles_SIM_GEANT->at(k).numberOfDaughters());
			histos_th1f[b+"GEN_Ks_status_status8"]->Fill(h_genParticles_SIM_GEANT->at(k).status());
			histos_th2f[b+"GEN_Ks_status_n_daughters_status8"]->Fill(h_genParticles_SIM_GEANT->at(k).status(),h_genParticles_SIM_GEANT->at(k).numberOfDaughters());
		
	}
     }
  }
  //the kinematics of the RECO Ks
  if(h_V0Ks.isValid()){
	for(unsigned int k = 0; k<h_V0Ks->size(); k++){
		histos_th1f[b+"V0s_Ks_pt"]->Fill(h_V0Ks->at(k).pt());				
		histos_th1f[b+"V0s_Ks_eta"]->Fill(h_V0Ks->at(k).eta());				
		histos_th1f[b+"V0s_Ks_phi"]->Fill(h_V0Ks->at(k).phi());				
		TVector3 xyz_Ks(h_V0Ks->at(k).vx(),h_V0Ks->at(k).vy(),h_V0Ks->at(k).vz());
		Double_t lxy_Ks = lxy(beamspot,xyz_Ks);
		histos_th1f[b+"V0s_Ks_lxy"]->Fill(lxy_Ks);
		histos_th1f[b+"V0s_Ks_vz"]->Fill(h_V0Ks->at(k).vz());
		histos_th1f[b+"V0s_Ks_ndaughters"]->Fill(h_V0Ks->at(k).numberOfDaughters());
		
	}
  }

  //the kinematics of the GenParticlesPlusGeant L
  if(h_genParticles_SIM_GEANT.isValid())
  {
    for(unsigned int l = 0; l < h_genParticles_SIM_GEANT->size(); l++){
        if(fabs(h_genParticles_SIM_GEANT->at(l).pdgId())==3122){
			histos_th1f[b+"GEN_L_pt"]->Fill(h_genParticles_SIM_GEANT->at(l).pt());				
			if(h_genParticles_SIM_GEANT->at(l).numberOfMothers()>0)if(h_genParticles_SIM_GEANT->at(l).mother()->pdgId()==-1020000020){
				histos_th1f[b+"GEN_L_pt_daughter_antiS"]->Fill(h_genParticles_SIM_GEANT->at(l).pt());				
				if(h_genParticles_SIM_GEANT->at(l).numberOfDaughters() == 2){
					Double_t pt_L_daughter1 = h_genParticles_SIM_GEANT->at(l).daughter(0)->pt();
					Double_t pt_L_daughter2 = h_genParticles_SIM_GEANT->at(l).daughter(1)->pt();
					Double_t min_pt_L_daughter = std::min(pt_L_daughter1,pt_L_daughter2);
					Double_t max_pt_L_daughter = std::max(pt_L_daughter1,pt_L_daughter2);
                                        histos_th1f[b+"GEN_L_daughter_antiS_pt_of_L_daughters"]->Fill(pt_L_daughter1);
                                        histos_th1f[b+"GEN_L_daughter_antiS_pt_of_L_daughters"]->Fill(pt_L_daughter2);
                                        histos_th2f[b+"GEN_L_daughter_antiS_min_max_pt_of_L_daughters"]->Fill(min_pt_L_daughter,max_pt_L_daughter);

				}
			}
			histos_th1f[b+"GEN_L_eta"]->Fill(h_genParticles_SIM_GEANT->at(l).eta());				
			histos_th1f[b+"GEN_L_phi"]->Fill(h_genParticles_SIM_GEANT->at(l).phi());				
			TVector3 xyz_L(h_genParticles_SIM_GEANT->at(l).vx(),h_genParticles_SIM_GEANT->at(l).vy(),h_genParticles_SIM_GEANT->at(l).vz());
			Double_t lxy_L = lxy(beamspot,xyz_L);
			histos_th1f[b+"GEN_L_lxy"]->Fill(lxy_L);
			histos_th2f[b+"GEN_L_lxy_pt"]->Fill(lxy_L,h_genParticles_SIM_GEANT->at(l).pt());
			histos_th1f[b+"GEN_L_vz"]->Fill(h_genParticles_SIM_GEANT->at(l).vz());
			histos_th1f[b+"GEN_L_ndaughters"]->Fill(h_genParticles_SIM_GEANT->at(l).numberOfDaughters());
			histos_th1f[b+"GEN_L_status"]->Fill(h_genParticles_SIM_GEANT->at(l).status());
			histos_th2f[b+"GEN_L_status_n_daughters"]->Fill(h_genParticles_SIM_GEANT->at(l).status(),h_genParticles_SIM_GEANT->at(l).numberOfDaughters());
		
	}
     }
  }
  //the kinematics of the GenParticlesPlusGeant L of status 2
  if(h_genParticles_SIM_GEANT.isValid())
  {
    for(unsigned int l = 0; l < h_genParticles_SIM_GEANT->size(); l++){
        if(fabs(h_genParticles_SIM_GEANT->at(l).pdgId())==3122&& h_genParticles_SIM_GEANT->at(l).status()==1){
			histos_th1f[b+"GEN_L_pt_status1"]->Fill(h_genParticles_SIM_GEANT->at(l).pt());				
			histos_th1f[b+"GEN_L_eta_status1"]->Fill(h_genParticles_SIM_GEANT->at(l).eta());				
			histos_th1f[b+"GEN_L_phi_status1"]->Fill(h_genParticles_SIM_GEANT->at(l).phi());				
			TVector3 xyz_L(h_genParticles_SIM_GEANT->at(l).vx(),h_genParticles_SIM_GEANT->at(l).vy(),h_genParticles_SIM_GEANT->at(l).vz());
			Double_t lxy_L = lxy(beamspot,xyz_L);
			histos_th1f[b+"GEN_L_lxy_status1"]->Fill(lxy_L);
			histos_th2f[b+"GEN_L_lxy_pt_status1"]->Fill(lxy_L,h_genParticles_SIM_GEANT->at(l).pt());
			histos_th1f[b+"GEN_L_vz_status1"]->Fill(h_genParticles_SIM_GEANT->at(l).vz());
			histos_th1f[b+"GEN_L_ndaughters_status1"]->Fill(h_genParticles_SIM_GEANT->at(l).numberOfDaughters());
			histos_th1f[b+"GEN_L_status_status1"]->Fill(h_genParticles_SIM_GEANT->at(l).status());
			histos_th2f[b+"GEN_L_status_n_daughters_status1"]->Fill(h_genParticles_SIM_GEANT->at(l).status(),h_genParticles_SIM_GEANT->at(l).numberOfDaughters());
		
	}
     }
  }

  //the kinematics of the GenParticlesPlusGeant L with status 8
  if(h_genParticles_SIM_GEANT.isValid())
  {
    for(unsigned int l = 0; l < h_genParticles_SIM_GEANT->size(); l++){
        if(fabs(h_genParticles_SIM_GEANT->at(l).pdgId())==3122&& h_genParticles_SIM_GEANT->at(l).status()==8){
			histos_th1f[b+"GEN_L_pt_status8"]->Fill(h_genParticles_SIM_GEANT->at(l).pt());				
			histos_th1f[b+"GEN_L_eta_status8"]->Fill(h_genParticles_SIM_GEANT->at(l).eta());				
			histos_th1f[b+"GEN_L_phi_status8"]->Fill(h_genParticles_SIM_GEANT->at(l).phi());				
			TVector3 xyz_L(h_genParticles_SIM_GEANT->at(l).vx(),h_genParticles_SIM_GEANT->at(l).vy(),h_genParticles_SIM_GEANT->at(l).vz());
			Double_t lxy_L = lxy(beamspot,xyz_L);
			histos_th1f[b+"GEN_L_lxy_status8"]->Fill(lxy_L);
			histos_th2f[b+"GEN_L_lxy_pt_status8"]->Fill(lxy_L,h_genParticles_SIM_GEANT->at(l).pt());
			histos_th1f[b+"GEN_L_vz_status8"]->Fill(h_genParticles_SIM_GEANT->at(l).vz());
			histos_th1f[b+"GEN_L_ndaughters_status8"]->Fill(h_genParticles_SIM_GEANT->at(l).numberOfDaughters());
			histos_th1f[b+"GEN_L_status_status8"]->Fill(h_genParticles_SIM_GEANT->at(l).status());
			histos_th2f[b+"GEN_L_status_n_daughters_status8"]->Fill(h_genParticles_SIM_GEANT->at(l).status(),h_genParticles_SIM_GEANT->at(l).numberOfDaughters());
		
	}
     }
  }
  //the kinematics of the RECO L
  if(h_V0L.isValid()){
	for(unsigned int l = 0; l<h_V0L->size(); l++){
		histos_th1f[b+"V0s_L_pt"]->Fill(h_V0L->at(l).pt());				
		histos_th1f[b+"V0s_L_eta"]->Fill(h_V0L->at(l).eta());				
		histos_th1f[b+"V0s_L_phi"]->Fill(h_V0L->at(l).phi());				
		TVector3 xyz_L(h_V0L->at(l).vx(),h_V0L->at(l).vy(),h_V0L->at(l).vz());
		Double_t lxy_L = lxy(beamspot,xyz_L);
		histos_th1f[b+"V0s_L_lxy"]->Fill(lxy_L);
		histos_th1f[b+"V0s_L_vz"]->Fill(h_V0L->at(l).vz());
		histos_th1f[b+"V0s_L_ndaughters"]->Fill(h_V0L->at(l).numberOfDaughters());
		
	}
  }


  //look ath the SIM_GEANT, find the Ks and the Lambdas and check wether they got reconstructed
  if(h_genParticles_SIM_GEANT.isValid())
  {

    for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); i++){

        //if the gen particle is a Ks
        if(fabs(h_genParticles_SIM_GEANT->at(i).pdgId())==310){
		
                Double_t GEN_Ks_phi = h_genParticles_SIM_GEANT->at(i).phi();
                Double_t GEN_Ks_eta = h_genParticles_SIM_GEANT->at(i).eta();
		Double_t GEN_Ks_vx = h_genParticles_SIM_GEANT->at(i).vx();
		Double_t GEN_Ks_vy = h_genParticles_SIM_GEANT->at(i).vy();
		Double_t GEN_Ks_vz = h_genParticles_SIM_GEANT->at(i).vz();
		Double_t GEN_Ks_vxy = pow(GEN_Ks_vx*GEN_Ks_vx+GEN_Ks_vy*GEN_Ks_vy,0.5);
		TVector3 GEN_Ks_vxyz(GEN_Ks_vx,GEN_Ks_vy,GEN_Ks_vz);
		Double_t GEN_Ks_lxy = lxy(FirstOfflinePV,GEN_Ks_vxyz);
		//Double_t GEN_Ks_lxy = lxy(beamspot,GEN_Ks_vxyz);

		Double_t GEN_Ks_decay_lxy = 999;
		Double_t GEN_Ks_decay_vz = 999;
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2){
			Double_t GEN_Ks_decay_vx = h_genParticles_SIM_GEANT->at(i).daughter(0)->vx();
			Double_t GEN_Ks_decay_vy = h_genParticles_SIM_GEANT->at(i).daughter(0)->vy();
			GEN_Ks_decay_vz = h_genParticles_SIM_GEANT->at(i).daughter(0)->vz();

			TVector3 GEN_Ks_decay_vxyz(GEN_Ks_decay_vx,GEN_Ks_decay_vy,GEN_Ks_decay_vz);
			GEN_Ks_decay_lxy = lxy(beamspot,GEN_Ks_decay_vxyz);
		}

		TVector3 GEN_Ks_p(h_genParticles_SIM_GEANT->at(i).px(),h_genParticles_SIM_GEANT->at(i).py(),h_genParticles_SIM_GEANT->at(i).pz());
		Double_t GEN_Ks_dxy = dxy_signed_line_point(GEN_Ks_vxyz, GEN_Ks_p, beamspot);

		bool KsDaughtersInTrackerAcceptance = false;
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2){
			if(fabs(h_genParticles_SIM_GEANT->at(i).daughter(0)->eta()) < 2.5 && fabs(h_genParticles_SIM_GEANT->at(i).daughter(1)->eta()) < 2.5) KsDaughtersInTrackerAcceptance = true;
		}

		bool KsDaughtersAreChargedPions = false;
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2){
			if(  (h_genParticles_SIM_GEANT->at(i).daughter(0)->pdgId() == 211 && h_genParticles_SIM_GEANT->at(i).daughter(1)->pdgId() == -211 ) || (h_genParticles_SIM_GEANT->at(i).daughter(0)->pdgId() == -211 && h_genParticles_SIM_GEANT->at(i).daughter(1)->pdgId() == 211 )     ) KsDaughtersAreChargedPions = true; 
		}

                if(h_V0Ks.isValid() && KsDaughtersInTrackerAcceptance && KsDaughtersAreChargedPions){
                        for (unsigned int k = 0; k < h_V0Ks->size(); ++k) {
                                Double_t RECO_Ks_phi = h_V0Ks->at(k).phi();
                                Double_t RECO_Ks_eta = h_V0Ks->at(k).eta();
                                Double_t deltaR_Ks_GEN_RECO = deltaR(GEN_Ks_phi,GEN_Ks_eta,RECO_Ks_phi,RECO_Ks_eta);
                                histos_th1f[b+"V0s_Ks_deltaR_Ks_GEN_RECO"]->Fill(deltaR_Ks_GEN_RECO);
				if(h_genParticles_SIM_GEANT->at(i).status()==1)histos_th1f[b+"V0s_Ks_deltaR_Ks_GEN_RECO_status1"]->Fill(deltaR_Ks_GEN_RECO);
				if(h_genParticles_SIM_GEANT->at(i).status()==8)histos_th1f[b+"V0s_Ks_deltaR_Ks_GEN_RECO_status8"]->Fill(deltaR_Ks_GEN_RECO);
				if(h_genParticles_SIM_GEANT->at(i).numberOfMothers()>0)if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId()==-1020000020)histos_th1f[b+"V0s_Ks_deltaR_Ks_GEN_RECO_daughter_antiS"]->Fill(deltaR_Ks_GEN_RECO);
				bool matched = false;
                                if(deltaR_Ks_GEN_RECO<0.02){
					//plot the properties of the Ks which are not daughters of the antiS and which get reconstructed
					if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId() != -1020000020){
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_vxy"]->Fill(GEN_Ks_vxy);
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_lxy"]->Fill(GEN_Ks_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_Ks_reconstructed_GEN_decay_lxy"]->Fill(GEN_Ks_decay_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_Ks_reconstructed_GEN_decay_vz"]->Fill(GEN_Ks_decay_vz);
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_dxy"]->Fill(GEN_Ks_dxy);
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_status"]->Fill(h_genParticles_SIM_GEANT->at(i).status());
					}
					else{//plot the properties of the Ks which are daughters of the antiS and which get reconstructed
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_vxy"]->Fill(GEN_Ks_vxy);
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_lxy"]->Fill(GEN_Ks_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_decay_lxy"]->Fill(GEN_Ks_decay_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_decay_vz"]->Fill(GEN_Ks_decay_vz);
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_dxy"]->Fill(GEN_Ks_dxy);
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_status"]->Fill(h_genParticles_SIM_GEANT->at(i).status());
					}
					histos_th1f[b+"V0s_Ks_reconstructed"]->Fill(1);
					matched = true;				
                                }
                                else{

					//plot the properties of the Ks which are not daughters of the antiS and which do not get reconstructed
					if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId() != -1020000020){
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_vxy"]->Fill(GEN_Ks_vxy);
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_lxy"]->Fill(GEN_Ks_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_decay_lxy"]->Fill(GEN_Ks_decay_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_decay_vz"]->Fill(GEN_Ks_decay_vz);
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_dxy"]->Fill(GEN_Ks_dxy);
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_status"]->Fill(h_genParticles_SIM_GEANT->at(i).status());
					}
					else{//plot the properties of the Ks which are daughters of the antiS and which do not get reconstructed
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_vxy"]->Fill(GEN_Ks_vxy);
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_lxy"]->Fill(GEN_Ks_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_lxy"]->Fill(GEN_Ks_decay_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_vz"]->Fill(GEN_Ks_decay_vz);
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_dxy"]->Fill(GEN_Ks_dxy);
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_status"]->Fill(h_genParticles_SIM_GEANT->at(i).status());
					}
					histos_th1f[b+"V0s_Ks_reconstructed"]->Fill(0);

					
                                }
				histos_teff[b+"V0_Ks_reconstructed_pt"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).pt());
				histos_teff[b+"V0_Ks_reconstructed_p"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).p());
				histos_teff[b+"V0_Ks_reconstructed_vxy"]->Fill(matched,GEN_Ks_vxy);
				histos_teff[b+"V0_Ks_reconstructed_lxy"]->Fill(matched,GEN_Ks_lxy);
				if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_teff[b+"V0_Ks_reconstructed_decay_lxy"]->Fill(matched,GEN_Ks_decay_lxy);
				if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_teff[b+"V0_Ks_reconstructed_decay_vz"]->Fill(matched,GEN_Ks_decay_vz);
				histos_teff[b+"V0_Ks_reconstructed_dxy"]->Fill(matched,GEN_Ks_dxy);
				histos_teff[b+"V0_Ks_reconstructed_vz"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).vz());
				histos_teff[b+"V0_Ks_reconstructed_eta"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).eta());
				histos_teff[b+"V0_Ks_reconstructed_phi"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).phi());
				histos_teff[b+"V0_Ks_reconstructed_status"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).status());
				
                        }
                }
        }


        //if the gen particle is a L
        if(fabs(h_genParticles_SIM_GEANT->at(i).pdgId())==3122 ){
                Double_t GEN_L_phi = h_genParticles_SIM_GEANT->at(i).phi();
                Double_t GEN_L_eta = h_genParticles_SIM_GEANT->at(i).eta();
		Double_t GEN_L_vx = h_genParticles_SIM_GEANT->at(i).vx();
		Double_t GEN_L_vy = h_genParticles_SIM_GEANT->at(i).vy();
		Double_t GEN_L_vz = h_genParticles_SIM_GEANT->at(i).vz();
		Double_t GEN_L_vxy = pow(GEN_L_vx*GEN_L_vx+GEN_L_vy*GEN_L_vy,0.5);
		TVector3 GEN_L_vxyz(GEN_L_vx,GEN_L_vy,GEN_L_vz);
		//Double_t GEN_L_lxy = lxy(beamspot,GEN_L_vxyz);
		Double_t GEN_L_lxy = lxy(FirstOfflinePV,GEN_L_vxyz);

		Double_t GEN_L_decay_lxy = 999;
		Double_t GEN_L_decay_vz = 999;
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2){
			Double_t GEN_L_decay_vx = h_genParticles_SIM_GEANT->at(i).daughter(0)->vx();
			Double_t GEN_L_decay_vy = h_genParticles_SIM_GEANT->at(i).daughter(0)->vy();
			GEN_L_decay_vz = h_genParticles_SIM_GEANT->at(i).daughter(0)->vz();

			TVector3 GEN_L_decay_vxyz(GEN_L_decay_vx,GEN_L_decay_vy,GEN_L_decay_vz);
			GEN_L_decay_lxy = lxy(beamspot,GEN_L_decay_vxyz);
		}	

		TVector3 GEN_L_p(h_genParticles_SIM_GEANT->at(i).px(),h_genParticles_SIM_GEANT->at(i).py(),h_genParticles_SIM_GEANT->at(i).pz());
		Double_t GEN_L_dxy = dxy_signed_line_point(GEN_L_vxyz, GEN_L_p, beamspot);

	
		bool LDaughtersInTrackerAcceptance = false;
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2){
			if(fabs(h_genParticles_SIM_GEANT->at(i).daughter(0)->eta()) < 2.5 && fabs(h_genParticles_SIM_GEANT->at(i).daughter(1)->eta()) < 2.5) LDaughtersInTrackerAcceptance = true;
		}
	
		bool LDaughtersAreChargedPions = false;
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2){
			if(  (h_genParticles_SIM_GEANT->at(i).daughter(0)->pdgId() == -2212 && h_genParticles_SIM_GEANT->at(i).daughter(1)->pdgId() == 211 ) || (h_genParticles_SIM_GEANT->at(i).daughter(0)->pdgId() == 211 && h_genParticles_SIM_GEANT->at(i).daughter(1)->pdgId() == -2212 )     ) LDaughtersAreChargedPions = true; 
		}

                if(h_V0L.isValid() && LDaughtersInTrackerAcceptance && LDaughtersAreChargedPions){
                        for (unsigned int l = 0; l < h_V0L->size(); ++l) {
                                Double_t RECO_L_phi = h_V0L->at(l).phi();
                                Double_t RECO_L_eta = h_V0L->at(l).eta();
                                Double_t deltaR_L_GEN_RECO = deltaR(GEN_L_phi,GEN_L_eta,RECO_L_phi,RECO_L_eta);
                                histos_th1f[b+"V0s_L_deltaR_L_GEN_RECO"]->Fill(deltaR_L_GEN_RECO);
				if(h_genParticles_SIM_GEANT->at(i).status()==1)histos_th1f[b+"V0s_L_deltaR_L_GEN_RECO_status1"]->Fill(deltaR_L_GEN_RECO);
				if(h_genParticles_SIM_GEANT->at(i).status()==8)histos_th1f[b+"V0s_L_deltaR_L_GEN_RECO_status8"]->Fill(deltaR_L_GEN_RECO);
				if(h_genParticles_SIM_GEANT->at(i).numberOfMothers()>0)if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId()==-1020000020)histos_th1f[b+"V0s_L_deltaR_L_GEN_RECO_daughter_antiS"]->Fill(deltaR_L_GEN_RECO);
				bool matched = false;
                                if(deltaR_L_GEN_RECO<0.02){
					//plot the properties of the L which are not daughters of the antiS and which get reconstructed
					if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId() != -1020000020){
						histos_th1f[b+"V0s_L_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_L_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_L_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_L_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_L_reconstructed_GEN_vxy"]->Fill(GEN_L_vxy);
						histos_th1f[b+"V0s_L_reconstructed_GEN_lxy"]->Fill(GEN_L_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_L_reconstructed_GEN_decay_lxy"]->Fill(GEN_L_decay_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_L_reconstructed_GEN_decay_vz"]->Fill(GEN_L_decay_vz);
						histos_th1f[b+"V0s_L_reconstructed_GEN_dxy"]->Fill(GEN_L_dxy);
						histos_th1f[b+"V0s_L_reconstructed_GEN_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
						histos_th1f[b+"V0s_L_reconstructed_GEN_status"]->Fill(h_genParticles_SIM_GEANT->at(i).status());
					}
					//plot the properties of the L which are daughters of the antiS and which get reconstructed
					else{
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_vxy"]->Fill(GEN_L_vxy);
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_lxy"]->Fill(GEN_L_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_decay_lxy"]->Fill(GEN_L_decay_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_decay_vz"]->Fill(GEN_L_decay_vz);
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_dxy"]->Fill(GEN_L_dxy);
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_status"]->Fill(h_genParticles_SIM_GEANT->at(i).status());
					}
					histos_th1f[b+"V0s_L_reconstructed"]->Fill(1);
					matched = true;	
                                }
                                else{
					//plot the properties of the L which are not daughters of the antiS and which do not get reconstructed
					if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId() != -1020000020){
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_vxy"]->Fill(GEN_L_vxy);
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_lxy"]->Fill(GEN_L_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_L_non_reconstructed_GEN_decay_lxy"]->Fill(GEN_L_decay_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_L_non_reconstructed_GEN_decay_vz"]->Fill(GEN_L_decay_vz);
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_dxy"]->Fill(GEN_L_dxy);
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_status"]->Fill(h_genParticles_SIM_GEANT->at(i).status());
					}
					else{//plot the properties of the L which are daughters of the antiS and which do not get reconstructed
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_vxy"]->Fill(GEN_L_vxy);
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_lxy"]->Fill(GEN_L_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_decay_lxy"]->Fill(GEN_L_decay_lxy);
						if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_decay_vz"]->Fill(GEN_L_decay_vz);
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_dxy"]->Fill(GEN_L_dxy);
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_status"]->Fill(h_genParticles_SIM_GEANT->at(i).status());
					}
					histos_th1f[b+"V0s_L_reconstructed"]->Fill(0);
					
                                }

				histos_teff[b+"V0_L_reconstructed_pt"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).pt());
				histos_teff[b+"V0_L_reconstructed_p"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).p());
				histos_teff[b+"V0_L_reconstructed_vxy"]->Fill(matched,GEN_L_vxy);
				histos_teff[b+"V0_L_reconstructed_lxy"]->Fill(matched,GEN_L_lxy);
				if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_teff[b+"V0_L_reconstructed_decay_lxy"]->Fill(matched,GEN_L_decay_lxy);
				if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2)histos_teff[b+"V0_L_reconstructed_decay_vz"]->Fill(matched,GEN_L_decay_vz);
				histos_teff[b+"V0_L_reconstructed_dxy"]->Fill(matched,GEN_L_dxy);
				histos_teff[b+"V0_L_reconstructed_vz"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).vz());
				histos_teff[b+"V0_L_reconstructed_eta"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).eta());
				histos_teff[b+"V0_L_reconstructed_phi"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).phi());
				histos_teff[b+"V0_L_reconstructed_status"]->Fill(matched,h_genParticles_SIM_GEANT->at(i).status());
				
                        }
                }
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
		histos_th1f[b+"h_simgeantParticles_antiS_status"]->Fill(antiS.status());
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
		       // cout << "FOUND AN ANTI-S ONLY DECAYING TO ONE PARTICLE...: "<< antiS.daughter(0)->pdgId() << endl;	

			TLorentzVector p4n(0,0,0,0);
			TLorentzVector p4antiS(0,0,0,0);
			TLorentzVector p4Daughter(0,0,0,0);

			p4n.SetPxPyPzE(0,0,0,0.939565);
			p4antiS.SetPxPyPzE(antiS.px(),antiS.py(),antiS.pz(),antiS.energy());

			TLorentzVector p4antiSn = p4antiS + p4n;

			p4Daughter.SetPxPyPzE(antiS.daughter(0)->px(),antiS.daughter(0)->py(),antiS.daughter(0)->pz(),antiS.daughter(0)->energy());
			//cout << "S+n energy: "  << p4antiSn.Energy() << endl;
			//cout << "daughter energy: " << p4Daughter.Energy() << endl;
			
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
			
		        //cout << "FOUND AN ANTI-S DECAYING TO TWO PARTICLES...: " << endl;	
		        TLorentzVector p4n(0,0,0,0);
                        TLorentzVector p4antiS(0,0,0,0);
                        TLorentzVector p4Daughter0(0,0,0,0);
                        TLorentzVector p4Daughter1(0,0,0,0);

                        p4n.SetPxPyPzE(0,0,0,0.939565);
                        p4antiS.SetPxPyPzE(antiS.px(),antiS.py(),antiS.pz(),antiS.energy());

                        TLorentzVector p4antiSn = p4antiS + p4n;

                        p4Daughter0.SetPxPyPzE(antiS.daughter(0)->px(),antiS.daughter(0)->py(),antiS.daughter(0)->pz(),antiS.daughter(0)->energy());
                        p4Daughter1.SetPxPyPzE(antiS.daughter(1)->px(),antiS.daughter(1)->py(),antiS.daughter(1)->pz(),antiS.daughter(1)->energy());
                        //cout << "S+n energy: "  << p4antiSn.Energy() << endl;
                        //cout << "daughter0 energy: " << p4Daughter0.Energy() << endl;
                        //cout << "daughter1 energy: " << p4Daughter1.Energy() << endl;

			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pt"]->Fill(antiS.pt());
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_p"]->Fill(antiS.p());
			histos_th2f[b+"h2_simgeantParticles_antiS_2_daughters_p_eta"]->Fill(antiS.p(), antiS.eta());
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_eta"]->Fill(antiS.eta());
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_phi"]->Fill(antiS.phi());
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_vxy"]->Fill(pow(vx_antiS*vx_antiS+vy_antiS*vy_antiS,0.5));
			TVector3 vertex_V0(antiS.daughter(0)->vx(),antiS.daughter(0)->vy(),antiS.daughter(0)->vz());
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_lxy_V0s"]->Fill(lxy(vertex_V0,beamspot));
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_vz_V0s"]->Fill(antiS.daughter(0)->vz());

			double deltaPhiV0s = reco::deltaPhi(antiS.daughter(0)->phi(),antiS.daughter(1)->phi());
			double deltaEtaV0s = antiS.daughter(0)->eta() - antiS.daughter(1)->eta();
			double deltaRV0s = pow(deltaPhiV0s*deltaPhiV0s+deltaEtaV0s*deltaEtaV0s,0.5);
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_deltaPhi_V0s"]->Fill(deltaPhiV0s);
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_deltaEta_V0s"]->Fill(deltaEtaV0s);
			histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_deltaR_V0s"]->Fill(deltaRV0s);

			if(antiS.daughter(0)->numberOfDaughters() ==2 && antiS.daughter(1)->numberOfDaughters() == 2){
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pdgId_granddaughters"]->Fill(antiS.daughter(0)->daughter(0)->pdgId());	
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pdgId_granddaughters"]->Fill(antiS.daughter(0)->daughter(1)->pdgId());	
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pdgId_granddaughters"]->Fill(antiS.daughter(1)->daughter(0)->pdgId());	
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pdgId_granddaughters"]->Fill(antiS.daughter(1)->daughter(1)->pdgId());	

				TVector3 vertex_V0_Ks_daug0(antiS.daughter(0)->daughter(0)->vx(),antiS.daughter(0)->daughter(0)->vy(),antiS.daughter(0)->daughter(0)->vz());
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_lxy_granddaughters"]->Fill(lxy(vertex_V0_Ks_daug0,beamspot));
				TVector3 vertex_V0_Ks_daug1(antiS.daughter(0)->daughter(1)->vx(),antiS.daughter(0)->daughter(1)->vy(),antiS.daughter(0)->daughter(1)->vz());
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_lxy_granddaughters"]->Fill(lxy(vertex_V0_Ks_daug1,beamspot));
				TVector3 vertex_V0_antiL_daug0(antiS.daughter(1)->daughter(0)->vx(),antiS.daughter(1)->daughter(0)->vy(),antiS.daughter(1)->daughter(0)->vz());
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_lxy_granddaughters"]->Fill(lxy(vertex_V0_antiL_daug0,beamspot));
				TVector3 vertex_V0_antiL_daug1(antiS.daughter(1)->daughter(1)->vx(),antiS.daughter(1)->daughter(1)->vy(),antiS.daughter(1)->daughter(1)->vz());
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_lxy_granddaughters"]->Fill(lxy(vertex_V0_antiL_daug1,beamspot));

				
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pt_Ksdaughters"]->Fill(antiS.daughter(0)->daughter(0)->pt());
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pt_Ksdaughters"]->Fill(antiS.daughter(0)->daughter(1)->pt());
				
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pt_AntiLdaughters"]->Fill(antiS.daughter(1)->daughter(0)->pt());
				histos_th1f[b+"h_simgeantParticles_antiS_2_daughters_pt_AntiLdaughters"]->Fill(antiS.daughter(1)->daughter(1)->pt());
				
			}
	
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
			
			//plotting the pdgids of the V0 daughters	
			histos_th2f[b+"h_simgeantParticles_antiS_n_daughters_Ks_n_daughters_L"]->Fill(Ks->numberOfDaughters(), L->numberOfDaughters());
			if(Ks->numberOfDaughters()==2 && L->numberOfDaughters()==2){
				histos_th2f[b+"h_simgeantParticles_antiS_2_daughters_pdgId_Ks_daughters_absdiff_pdgId_L_daughters_absdiff_pdgId"]->Fill(fabs(Ks->daughter(0)->pdgId()-Ks->daughter(1)->pdgId()), fabs(L->daughter(0)->pdgId()-L->daughter(1)->pdgId()));
			}
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
                         histos_th2f[b+"h_simgeantParticles_antiS_pt_corr_daughters"]->Fill(Ks->pt(),L->pt());
			 //just check wether the delta phi is smaller when the S has higher momentum 
			 if(antiS.pt()>1)histos_th1f[b+"h_simgeantParticles_antiS_delta_phi_daughters_S_pt_larger_1GeV"]->Fill(delta_phi_Ks_L);
			





			//now check for the daughters (Ks and L) of the GEN antiS if you can find a matching V0 candidate

			int nantiSKsDaughtersReconstucted = 0;
			if(h_V0Ks.isValid() ){
				cout << "-----RECO V0 Ks-----" << endl;
				for (unsigned int k = 0; k < h_V0Ks->size(); ++k) {
					cout << "RECO VO Ks #" << k << endl;
					cout << "phi, eta " << h_V0Ks->at(k).phi() << " " << h_V0Ks->at(k).eta()   << endl;
					cout << "px, py, pz " << h_V0Ks->at(k).px() << " " << h_V0Ks->at(k).py() << " " << h_V0Ks->at(k).pz()   << endl;
					cout << "vx, vy, vz " << h_V0Ks->at(k).vx() << " " << h_V0Ks->at(k).vy() << " " << h_V0Ks->at(k).vz()   << endl;
					Double_t RECO_V0_Ks_phi = h_V0Ks->at(k).phi();
					Double_t RECO_V0_Ks_eta = h_V0Ks->at(k).eta();
					Double_t deltaR_GEN_RECO_Ks = deltaR(RECO_V0_Ks_phi, RECO_V0_Ks_eta, Ks->phi(), Ks->eta());
					if(deltaR_GEN_RECO_Ks<0.1){
						cout << "found a RECO V0 Ks overlapping with a GEN Ks dauhgter of the antiS" << endl;
						nantiSKsDaughtersReconstucted++;
					}
				}//(unsigned int i = 0; i < h_V0Ks->size(); ++i)
			}//(h_V0Ks.isValid())

			int nantiSLDaughtersReconstucted = 0;
			if(h_V0L.isValid() ){
				cout << "-----RECO V0 L-----" << endl;
				for (unsigned int l = 0; l < h_V0L->size(); ++l) {
					cout << "RECO VO L #" << l << endl;
					cout << "phi, eta " << h_V0L->at(l).phi() << " " << h_V0L->at(l).eta()   << endl;
					cout << "px, py, pz " << h_V0L->at(l).px() << " " << h_V0L->at(l).py() << " " << h_V0L->at(l).pz()   << endl;
					cout << "vx, vy, vz " << h_V0L->at(l).vx() << " " << h_V0L->at(l).vy() << " " << h_V0L->at(l).vz()   << endl;
					Double_t RECO_V0_L_phi = h_V0L->at(l).phi();
					Double_t RECO_V0_L_eta = h_V0L->at(l).eta();
					Double_t deltaR_GEN_RECO_L = deltaR(RECO_V0_L_phi, RECO_V0_L_eta, L->phi(), L->eta());
					if(deltaR_GEN_RECO_L<0.1){
						cout << "found a RECO V0 L overlapping with a GEN L dauhgter of the antiS" << endl;
						nantiSLDaughtersReconstucted++;		
					}
				}//(unsigned int i = 0; i < h_V0L->size(); ++i)
			}
				
			histos_th2f[b+"h_V0_nantiSLDaughtersReconstucted_nantiSKsDaughtersReconstucted"]->Fill(nantiSKsDaughtersReconstucted, nantiSLDaughtersReconstucted);








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
			 bool KsDaughtersSurvEtaCut = false;
			 bool LDaughtersSurvEtaCut = false;
			 if(n_daughters_Ks == 2){
				if(fabs(Ks->daughter(0)->eta()) < 2.5 && fabs(Ks->daughter(1)->eta()) < 2.5 ) KsDaughtersSurvEtaCut = true;
			 }
			 if( n_daughters_L == 2){
				if(fabs(L->daughter(0)->eta()) < 2.5 && fabs(L->daughter(1)->eta()) < 2.5) LDaughtersSurvEtaCut = true;
			 }


			 bool antiSSurvDeltaPhiCut = false;
			 if(fabs(delta_phi_Ks_L)>0.5) antiSSurvDeltaPhiCut = true;

			//make some plots to See which cuts are the difficult ones to survive
			histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(0);
			if(KsSurvPtCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(1);
			if(LSurvPtCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(2);
			if(KsSurvPtCut && LSurvPtCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(3);
			if(KsDaughtersSurvEtaCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(4);
			if(LDaughtersSurvEtaCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(5);
			if(KsDaughtersSurvEtaCut && LDaughtersSurvEtaCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(6);
			if(KsSurvPtCut && KsDaughtersSurvEtaCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(7);
			if(LSurvPtCut && LDaughtersSurvEtaCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(8);
			if(antiSSurvDeltaPhiCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(9);
			if(KsSurvPtCut && LSurvPtCut && KsDaughtersSurvEtaCut && LDaughtersSurvEtaCut && antiSSurvDeltaPhiCut) histos_th1f[b+"h_simgeantParticles_antiS_cut_flow"]->Fill(10);
			
			
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

			 if(KsSurvPtCut && LSurvPtCut && KsDaughtersSurvEtaCut && LDaughtersSurvEtaCut && antiSSurvDeltaPhiCut){
			
				if(fabs(delta_phi_Ks_L) > 0.5 && fabs(delta_eta_Ks_L) < 2 && openings_angle(Ks->momentum(), L->momentum()) > 0.5 && openings_angle(Ks->momentum(), L->momentum()) < 2 && delta_R_Ks_L < 3 && p4_L_boosted_COM.Pz() > -p4_Ks_boosted_COM.Pz() - 0.3 && p4_L_boosted_COM.Pz() < -p4_Ks_boosted_COM.Pz() + 0.3){
					 histos_th1f[b+"h_simgeantParticles_antiS_n_surv_background_cuts"]->Fill(1);
				}
				else{
					 histos_th1f[b+"h_simgeantParticles_antiS_n_surv_background_cuts"]->Fill(0);
				}
				//cout << "fabs(delta_phi_Ks_L) " << fabs(delta_phi_Ks_L) << endl;
				//cout << "fabs(delta_eta_Ks_L) " << fabs(delta_eta_Ks_L) << endl;
				//cout << "openings_angle(Ks->momentum(), L->momentum()) " << openings_angle(Ks->momentum(), L->momentum()) << endl;
				//cout << "delta_R_Ks_L " <<  delta_R_Ks_L << endl;
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
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_status"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->status());
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
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_Ks_n_daughters"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->numberOfDaughters());
                    }
                    //if daugther is L
                    if(pdgId_antiS_daughter == -3122){
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_pdgid"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->pdgId());
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_status"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->status());
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
                        histos_th1f[b+"h_simgeantParticles_antiS_daughter_L_n_daughters"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(j)->numberOfDaughters());

                    }
                }		

		
		//now will start matching these PLUSGEANT particles to RECO particles
		if(h_sCands.isValid()) {
			unsigned int n_sCands = h_sCands->size();
        		for (unsigned int i = 0; i < n_sCands; ++i) {
				Double_t RECO_antiS_phi = h_sCands->at(i).phi();
				Double_t RECO_antiS_eta = h_sCands->at(i).eta();
				Double_t RECO_L_phi = h_sCands->at(i).daughter(0)->phi();
				Double_t RECO_L_eta = h_sCands->at(i).daughter(0)->eta();
				Double_t RECO_Ks_phi = h_sCands->at(i).daughter(1)->phi();
				Double_t RECO_Ks_eta = h_sCands->at(i).daughter(1)->eta();

				bool antiS_isReconstructed = false;
				bool L_isReconstructed = false;
				bool Ks_isReconstructed = false;
				
				Double_t deltaR_GEN_RECO_antiS = deltaR(RECO_antiS_phi, RECO_antiS_eta, antiS.phi(), antiS.eta());
				Double_t deltaR_GEN_RECO_L = deltaR(RECO_L_phi, RECO_L_eta, L->phi(), L->eta());
				Double_t deltaR_GEN_RECO_Ks = deltaR(RECO_Ks_phi, RECO_Ks_eta, Ks->phi(), Ks->eta());

				if(deltaR_GEN_RECO_antiS<0.1)antiS_isReconstructed=true;
				if(deltaR_GEN_RECO_L<0.1)L_isReconstructed=true;
				if(deltaR_GEN_RECO_Ks<0.1)Ks_isReconstructed=true;

				if(antiS_isReconstructed){
					cout << "an antiS got reconstructed" << endl;
					cout << "GEN antiS phi, eta: " << antiS.phi() << ", " << antiS.eta() << endl;
					cout << "RECO antiS phi, eta: "<< RECO_antiS_phi << ", " << RECO_antiS_eta  <<  endl;
				}
				if(L_isReconstructed){
					cout << "a Lambda as daughter of the antiS got reconstructed" << endl;
					cout << "GEN L phi, eta: " << L->phi() << ", " << L->eta() << endl;
					cout << "RECO L phi, eta: " << RECO_L_phi << ", " << RECO_L_eta << endl;
				}
				if(Ks_isReconstructed){
					cout << "a Ks as daughter of the antiS got reconstructed" << endl;
					cout << "GEN Ks phi, eta: " << Ks->phi() << ", " << Ks->eta() << endl;
					cout << "RECO Ks phi, eta: " << RECO_Ks_phi << ", " << RECO_Ks_eta << endl;
				}

				//Ks and L are available as the GEN Ks and Lambda		
				//antiS is available as the GEN antiS
						
			}//end for (unsigned int i = 0; i < n_sCands; ++i) 
		}//end h_sCands.isValid()

	
	}// if(n_daughters_antiS == 2)
      }//if(h_genParticles_SIM_GEANT->at(i).pdgId() == -1020000020)

    }//for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i)
  }//if(h_genParticles_SIM_GEANT.isValid())
 

if(h_sCands.isValid()) {
        unsigned int n_sCands = h_sCands->size();
        for (unsigned int i = 0; i < n_sCands; ++i) { 
			//the reconstructed particle is actually the S+n particle so have to subtract the neutron from the S candidate p4
			Double_t Sn_mass = h_sCands->at(i).mass();
			TLorentzVector p4n(0,0,0,0);
			p4n.SetPxPyPzE(0,0,0,0.939565);
			TLorentzVector p4Sn(0,0,0,0);
			p4Sn.SetPxPyPzE(h_sCands->at(i).px(),h_sCands->at(i).py(),h_sCands->at(i).pz(),h_sCands->at(i).energy());
			TLorentzVector p4S = p4Sn-p4n;
			Double_t S_mass = p4S.M();
			cout << "+++++++++++++++++++++++++++++++++++" << endl;
			cout << "using 4 momenta Sn_mass: " << Sn_mass << endl;
			cout << "using 4 momenta S_mass: " << S_mass << endl;
			histos_th1f[b+"h_LambdaKshortVertexFilter_S_mass"]->Fill(S_mass);
			histos_th1f[b+"h_LambdaKshortVertexFilter_Sn_mass"]->Fill(Sn_mass);
			histos_th2f[b+"h_LambdaKshortVertexFilter_S_vx_vy"]->Fill(h_sCands->at(i).vx(), h_sCands->at(i).vy());

			TVector3 S_vertex(h_sCands->at(i).vx(),h_sCands->at(i).vy(),h_sCands->at(i).vz());			
			Double_t lxy_S_b = lxy(S_vertex,beamspot);
			
			if(lxy_S_b > 1.2){		
				histos_th1f[b+"h_LambdaKshortVertexFilter_S_mass_with_displacement_larger_1p2"]->Fill(S_mass);
			}
			
	}//for (unsigned int i = 0; i < n_sCands; ++i)
  }//if(h_sCands.isValid())

} //end of analyzer

double Analyzer_SIM_Sexaq::openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2){
  double opening_angle = TMath::ACos((momentum1.Dot(momentum2))/(pow(momentum1.Mag2()*momentum2.Mag2(),0.5)));
  return opening_angle;
}

double Analyzer_SIM_Sexaq::deltaR(double phi1, double eta1, double phi2, double eta2){
	double deltaPhi = reco::deltaPhi(phi1,phi2);
	double deltaEta = eta1-eta2;
	return pow(deltaPhi*deltaPhi+deltaEta*deltaEta,0.5);
}


double Analyzer_SIM_Sexaq::lxy(TVector3 v1, TVector3 v2){
	double x1 = v1.X();
	double x2 = v2.X();
	double y1 = v1.Y();
	double y2 = v2.Y();
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}


TVector3 Analyzer_SIM_Sexaq::PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point){
   //first move the vector along the line to the starting point of Point_line
   double normalise = sqrt(Vector_along_line.X()*Vector_along_line.X()+Vector_along_line.Y()*Vector_along_line.Y()+Vector_along_line.Z()*Vector_along_line.Z());
   TVector3 n(Vector_along_line.X()/normalise,Vector_along_line.Y()/normalise,Vector_along_line.Z()/normalise);
   TVector3 a = Point_line;
   TVector3 p = Point;

   //see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line (Vector formulation)
   TVector3 vector_PCA = (a-p)-((a-p)*n)*n;
   return vector_PCA ;
}

double Analyzer_SIM_Sexaq::dxy_signed_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in){

  //looking at XY, so put the Z component to 0 first
  TVector3 Point_line(Point_line_in.X(),Point_line_in.Y(),0.);
  TVector3 Vector_along_line(Vector_along_line_in.X(), Vector_along_line_in.Y(),0.);
  TVector3 Point(Point_in.X(), Point_in.Y(), 0.);

  TVector3 shortest_distance = PCA_line_point(Point_line,  Vector_along_line, Point);
  double dxy_signed_line_point = sqrt(shortest_distance.X()*shortest_distance.X()+shortest_distance.Y()*shortest_distance.Y());

  TVector3 displacement = Point_line - Point; 
  if(displacement*Vector_along_line<0)dxy_signed_line_point = -dxy_signed_line_point;

  return dxy_signed_line_point;
}

void Analyzer_SIM_Sexaq::endJob()
{
}

Analyzer_SIM_Sexaq::~Analyzer_SIM_Sexaq()
{
}


DEFINE_FWK_MODULE(Analyzer_SIM_Sexaq);
