#include "../interface/Analyzer_Steven.h"
#include <typeinfo>


Analyzer_Steven::Analyzer_Steven(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  m_bsTag(pset.getParameter<edm::InputTag>("beamspot")),
  m_genParticlesTag_GEN(pset.getParameter<edm::InputTag>("genCollection_GEN")),
  m_genParticlesTag_SIM_GEANT(pset.getParameter<edm::InputTag>("genCollection_SIM_GEANT")),
  m_generalTracksTag(pset.getParameter<edm::InputTag>("generalTracksCollection")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_V0KsTag(pset.getParameter<edm::InputTag>("V0KsCollection")),
  m_V0LTag(pset.getParameter<edm::InputTag>("V0LCollection")),

  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN)),
  m_genParticlesToken_SIM_GEANT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_SIM_GEANT)),
  m_generalTracksToken(consumes<vector<reco::Track> >(m_generalTracksTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_V0KsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0KsTag)),
  m_V0LToken(consumes<vector<reco::VertexCompositeCandidate> >(m_V0LTag))


{

}


void Analyzer_Steven::beginJob() {
 
    //for the tracks
    TFileDirectory dir_generalTracks = m_fs->mkdir("generalTracks");
    histos_th1f[b+"h_tracks_deltaR"]= dir_generalTracks.make<TH1F>(b+"h_tracks_deltaR", b+"h_tracks_deltaR; deltaR(track, daughter of V0) ",1000,0,10);
    histos_th1f[b+"h_tracks_deltaR_min"]= dir_generalTracks.make<TH1F>(b+"h_tracks_deltaR_min", b+"h_tracks_deltaR_min; deltaR(track, daughter of V0) ",1000,0,10);
   histos_teff[b+"tracks_reco_eff_daughters_antiS"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_daughters_antiS",b+"tracks_reco_eff_daughters_antiS; 0 = tracks no daughters of antiS, 1 = tracks daughters of the antiS",2,0,2);
   histos_teff[b+"tracks_reco_eff_daughters_antiS_pt"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_daughters_antiS_pt",b+"tracks_reco_eff_daughters_antiS_pt; pt of V0 daughters (GeV)",100,0,10);
   histos_teff[b+"tracks_reco_eff_daughters_antiS_vxy"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_daughters_antiS_vxy",b+"tracks_reco_eff_daughters_antiS_vxy; vxy of V0 daughters (cm)",200,0,100);
   histos_teff[b+"tracks_reco_eff_daughters_antiS_surv_track_cuts_pt"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_daughters_antiS_surv_track_cuts_pt",b+"tracks_reco_eff_daughters_antiS_surv_track_cuts_pt; pt of V0 daughters (GeV)",100,0,10);
   histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_no_daughters_antiS_pt",b+"tracks_reco_eff_no_daughters_antiS_pt; pt of V0 daughters (GeV)",100,0,10);
   histos_teff[b+"tracks_reco_eff_no_daughters_antiS_vxy"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_no_daughters_antiS_vxy",b+"tracks_reco_eff_no_daughters_antiS_vxy; vxy of V0 daughters (cm)",200,0,100);
   histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5"] = dir_generalTracks.make<TEfficiency>(b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5",b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5; pt of V0 daughters (GeV)",100,0,10);



  


   //V0 particles
   TFileDirectory dir_V0s = m_fs->mkdir("V0s");
   //Ks 
   //kinematics of the GEN Ks
   TFileDirectory dir_V0s_GEN_Ks_kinematics = dir_V0s.mkdir("GEN_Ks_kinematics");
   histos_th1f[b+"GEN_Ks_pt"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_pt", b+"GEN_Ks_pt; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_Ks_eta"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_eta", b+"GEN_Ks_eta; eta", 100, -4, 4);
   histos_th1f[b+"GEN_Ks_phi"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_phi", b+"GEN_Ks_phi; phi(rad)", 200, -10, 10);
   histos_th1f[b+"GEN_Ks_lxy"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_lxy", b+"GEN_Ks_lxy; lxy(cm)", 1000, 0, 100);
   histos_th2f[b+"GEN_Ks_lxy_pt"]= dir_V0s_GEN_Ks_kinematics.make<TH2F>(b+"GEN_Ks_lxy_pt", b+"GEN_Ks_lxy_pt; lxy(cm); pt(GeV)", 1000, 0, 100, 100, 0, 10);
   histos_th1f[b+"GEN_Ks_vz"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_vz", b+"GEN_Ks_vz; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"GEN_Ks_ndaughters"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_ndaughters", b+"GEN_Ks_ndaughters; #daughters", 10, 0, 10);
   histos_th1f[b+"GEN_Ks_status"]= dir_V0s_GEN_Ks_kinematics.make<TH1F>(b+"GEN_Ks_status", b+"GEN_Ks_status; status", 10, 0, 10);
   histos_th2f[b+"GEN_Ks_status_n_daughters"]= dir_V0s_GEN_Ks_kinematics.make<TH2F>(b+"GEN_Ks_status_n_daughters", b+"GEN_Ks_status_n_daughters; status; n daughters", 10, 0, 10, 10, 0, 10);


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

   histos_th1f[b+"V0s_Ks_reconstructed_GEN_pt"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_pt", b+"V0s_Ks_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_p"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_p", b+"V0s_Ks_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_eta"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_eta", b+"V0s_Ks_reconstructed_GEN_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_phi"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_phi", b+"V0s_Ks_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_vxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_vxy", b+"V0s_Ks_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_lxy", b+"V0s_Ks_reconstructed_GEN_lxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_decay_lxy", b+"V0s_Ks_reconstructed_GEN_decay_lxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_dxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_dxy", b+"V0s_Ks_reconstructed_GEN_dxy; dxy (beamspot, Ks) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_vz", b+"V0s_Ks_reconstructed_GEN_vz; vz (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_Ks_reconstructed_GEN_status"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_reconstructed_GEN_status", b+"V0s_Ks_reconstructed_GEN_status; status", 10, 0, 10);

   //matched Ks to GEN, with Ks a daughter of the S
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_pt"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_pt", b+"V0s_Ks_daughterS_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_p"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_p", b+"V0s_Ks_daughterS_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_eta"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_eta", b+"V0s_Ks_daughterS_reconstructed_GEN_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_phi"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_phi", b+"V0s_Ks_daughterS_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_vxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_vxy", b+"V0s_Ks_daughterS_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_lxy", b+"V0s_Ks_daughterS_reconstructed_GEN_lxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_decay_lxy", b+"V0s_Ks_daughterS_reconstructed_GEN_decay_lxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_dxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_dxy", b+"V0s_Ks_daughterS_reconstructed_GEN_dxy; dxy (beamspot, Ks) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_vz", b+"V0s_Ks_daughterS_reconstructed_GEN_vz; vz (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_status"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_reconstructed_GEN_status", b+"V0s_Ks_daughterS_reconstructed_GEN_status; status", 10, 0, 10);

   //non matched Ks to GEN, with the Ks not a daughter of the S
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_pt"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_pt", b+"V0s_Ks_non_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_p"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_p", b+"V0s_Ks_non_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_eta"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_eta", b+"V0s_Ks_non_reconstructed_GEN_eta; eta ", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_phi"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_phi", b+"V0s_Ks_non_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_vxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_vxy", b+"V0s_Ks_non_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_lxy", b+"V0s_Ks_non_reconstructed_GEN_lxy; lxy(beamspot, Ks creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_decay_lxy", b+"V0s_Ks_non_reconstructed_GEN_decay_lxy; lxy (beamspot, Ks decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_dxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_dxy", b+"V0s_Ks_non_reconstructed_GEN_dxy; dxy(beamspot, Ks) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_vz", b+"V0s_Ks_non_reconstructed_GEN_vz; vz (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_status"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_non_reconstructed_GEN_status", b+"V0s_Ks_non_reconstructed_GEN_status; status", 10, 0, 10);

   //non matched Ks to GEN, with the Ks a daughter of the S
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_pt"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_pt", b+"V0s_Ks_daughterS_non_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_p"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_p", b+"V0s_Ks_daughterS_non_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_eta"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_eta", b+"V0s_Ks_daughterS_non_reconstructed_GEN_eta; eta ", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_phi"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_phi", b+"V0s_Ks_daughterS_non_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_vxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_vxy", b+"V0s_Ks_daughterS_non_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_lxy", b+"V0s_Ks_daughterS_non_reconstructed_GEN_lxy; lxy(beamspot, Ks creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_lxy", b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_lxy; lxy (beamspot, Ks decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_dxy"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_dxy", b+"V0s_Ks_daughterS_non_reconstructed_GEN_dxy; dxy(beamspot, Ks) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_vz"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_vz", b+"V0s_Ks_daughterS_non_reconstructed_GEN_vz; vz (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_status"]= dir_V0s_matched_Ks.make<TH1F>(b+"V0s_Ks_daughterS_non_reconstructed_GEN_status", b+"V0s_Ks_daughterS_non_reconstructed_GEN_status; status", 10, 0, 10);

   histos_teff[b+"V0_Ks_reconstructed_pt"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_pt",b+"V0_Ks_reconstructed_pt; pt(GeV)",100,0,10);
   histos_teff[b+"V0_Ks_reconstructed_p"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_p",b+"V0_Ks_reconstructed_p; p(GeV)",300,0,30);
   histos_teff[b+"V0_Ks_reconstructed_vxy"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_vxy",b+"V0_Ks_reconstructed_vxy; vxy(cm)",100,0,100);
   histos_teff[b+"V0_Ks_reconstructed_lxy"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_lxy",b+"V0_Ks_reconstructed_lxy; lxy (beamspot, Ks creation vertex)(cm)",1000,0,100);
   histos_teff[b+"V0_Ks_reconstructed_decay_lxy"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_decay_lxy",b+"V0_Ks_reconstructed_decay_lxy; lxy (beamspot, Ks decay vertex)(cm)",200,0,100);
   histos_teff[b+"V0_Ks_reconstructed_dxy"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_dxy",b+"V0_Ks_reconstructed_dxy; dxy (beamspot, Ks)(cm)",400,-100,100);
   histos_teff[b+"V0_Ks_reconstructed_vz"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_vz",b+"V0_Ks_reconstructed_vz; vz(cm)",200,-100,100);
   histos_teff[b+"V0_Ks_reconstructed_eta"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_eta",b+"V0_Ks_reconstructed_eta; eta",100,-4,4);
   histos_teff[b+"V0_Ks_reconstructed_phi"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_phi",b+"V0_Ks_reconstructed_phi; phi(rad)",100,-4,4);
   histos_teff[b+"V0_Ks_reconstructed_status"] = dir_V0s_matched_Ks.make<TEfficiency>(b+"V0_Ks_reconstructed_status",b+"V0_Ks_reconstructed_status; status",10,0,10);


				
   //Lambda
   //kinematics of the GEN L
   TFileDirectory dir_V0s_GEN_L_kinematics = dir_V0s.mkdir("GEN_L_kinematics");
   histos_th1f[b+"GEN_L_pt"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_pt", b+"GEN_L_pt; pt(GeV)", 100, 0, 10);
   histos_th1f[b+"GEN_L_eta"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_eta", b+"GEN_L_eta; eta", 100, -4, 4);
   histos_th1f[b+"GEN_L_phi"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_phi", b+"GEN_L_phi; phi(rad)", 200, -10, 10);
   histos_th1f[b+"GEN_L_lxy"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_lxy", b+"GEN_L_lxy; lxy(cm)", 1000, 0, 100);
   histos_th2f[b+"GEN_L_lxy_pt"]= dir_V0s_GEN_L_kinematics.make<TH2F>(b+"GEN_L_lxy_pt", b+"GEN_L_lxy_pt; lxy(cm); pt(GeV)", 1000, 0, 100, 100, 0, 10);
   histos_th1f[b+"GEN_L_vz"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_vz", b+"GEN_L_vz; vz(cm)", 2000, -100, 100);
   histos_th1f[b+"GEN_L_ndaughters"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_ndaughters", b+"GEN_L_ndaughters; #daughters", 10, 0, 10);
   histos_th1f[b+"GEN_L_status"]= dir_V0s_GEN_L_kinematics.make<TH1F>(b+"GEN_L_status", b+"GEN_L_status; status", 10, 0, 10);
   histos_th2f[b+"GEN_L_status_n_daughters"]= dir_V0s_GEN_L_kinematics.make<TH2F>(b+"GEN_L_status_n_daughters", b+"GEN_L_status_n_daughters; status; n daughters", 10, 0, 10, 10, 0, 10);


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
   //matched L to GEN which are not daughters of the S
   histos_th1f[b+"V0s_L_reconstructed_GEN_pt"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_pt", b+"V0s_L_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_reconstructed_GEN_p"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_p", b+"V0s_L_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_reconstructed_GEN_eta"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_eta", b+"V0s_L_reconstructed_GEN_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_L_reconstructed_GEN_phi"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_phi", b+"V0s_L_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_L_reconstructed_GEN_vxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_vxy", b+"V0s_L_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_reconstructed_GEN_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_lxy", b+"V0s_L_reconstructed_GEN_lxy; lxy (beamspot, L creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_decay_lxy", b+"V0s_L_reconstructed_GEN_decay_lxy; lxy (beamspot, L decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_reconstructed_GEN_dxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_dxy", b+"V0s_L_reconstructed_GEN_dxy; dxy (beamspot, L) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_L_reconstructed_GEN_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_vz", b+"V0s_L_reconstructed_GEN_vz; vz (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_L_reconstructed_GEN_status"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_reconstructed_GEN_status", b+"V0s_L_reconstructed_GEN_status; status", 10, 0, 10);

   //matched L to GEN which are daughters of the S
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_pt"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_pt", b+"V0s_L_daughterS_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_p"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_p", b+"V0s_L_daughterS_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_eta"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_eta", b+"V0s_L_daughterS_reconstructed_GEN_eta; eta", 100, -4, 4);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_phi"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_phi", b+"V0s_L_daughterS_reconstructed_GEN_phi; phi (rad)", 100, -4, 4);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_vxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_vxy", b+"V0s_L_daughterS_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_lxy", b+"V0s_L_daughterS_reconstructed_GEN_lxy; lxy (beamspot, L creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_decay_lxy", b+"V0s_L_daughterS_reconstructed_GEN_decay_lxy; lxy (beamspot, L decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_dxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_dxy", b+"V0s_L_daughterS_reconstructed_GEN_dxy; dxy (beamspot, L) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_vz", b+"V0s_L_daughterS_reconstructed_GEN_vz; vz (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_status"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_reconstructed_GEN_status", b+"V0s_L_daughterS_reconstructed_GEN_status; status", 10, 0, 10);

   //non-matched L to GEN which are not daughters of the S
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_pt"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_pt", b+"V0s_L_non_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_p"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_p", b+"V0s_L_non_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_eta"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_eta", b+"V0s_L_non_reconstructed_GEN_eta; eta ", 100, -4, 4);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_phi"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_phi", b+"V0s_L_non_reconstructed_GEN_phi; phi(rad) ", 100, -4, 4);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_vxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_vxy", b+"V0s_L_non_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_lxy", b+"V0s_L_non_reconstructed_GEN_lxy; lxy (beamspot, L creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_decay_lxy", b+"V0s_L_non_reconstructed_GEN_decay_lxy; lxy (beamspot, L decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_dxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_dxy", b+"V0s_L_non_reconstructed_GEN_dxy; dxy (beamspot, L) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_vz", b+"V0s_L_non_reconstructed_GEN_vz; vz (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_L_non_reconstructed_GEN_status"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_non_reconstructed_GEN_status", b+"V0s_L_non_reconstructed_GEN_status; status", 10, 0, 10);

   //non-matched L to GEN which are daughters of the S
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_pt"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_pt", b+"V0s_L_daughterS_non_reconstructed_GEN_pt; pt (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_p"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_p", b+"V0s_L_daughterS_non_reconstructed_GEN_p; p (GeV)", 100, 0, 10);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_eta"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_eta", b+"V0s_L_daughterS_non_reconstructed_GEN_eta; eta ", 100, -4, 4);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_phi"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_phi", b+"V0s_L_daughterS_non_reconstructed_GEN_phi; phi(rad) ", 100, -4, 4);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_vxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_vxy", b+"V0s_L_daughterS_non_reconstructed_GEN_vxy; vxy (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_lxy", b+"V0s_L_daughterS_non_reconstructed_GEN_lxy; lxy (beamspot, L creation vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_decay_lxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_decay_lxy", b+"V0s_L_daughterS_non_reconstructed_GEN_decay_lxy; lxy (beamspot, L decay vertex) (cm)", 1000, 0, 100);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_dxy"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_dxy", b+"V0s_L_daughterS_non_reconstructed_GEN_dxy; dxy (beamspot, L) (cm)", 2000, -100, 100);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_vz"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_vz", b+"V0s_L_daughterS_non_reconstructed_GEN_vz; vz (cm)", 20000, -1000, 1000);
   histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_status"]= dir_V0s_matched_L.make<TH1F>(b+"V0s_L_daughterS_non_reconstructed_GEN_status", b+"V0s_L_daughterS_non_reconstructed_GEN_status; status", 10, 0, 10);


   histos_teff[b+"V0_L_reconstructed_pt"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_pt",b+"V0_L_reconstructed_pt; pt(GeV)",100,0,10);
   histos_teff[b+"V0_L_reconstructed_p"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_p",b+"V0_L_reconstructed_p; p(GeV)",300,0,30);
   histos_teff[b+"V0_L_reconstructed_vxy"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_vxy",b+"V0_L_reconstructed_vxy; vxy(cm)",100,0,100);
   histos_teff[b+"V0_L_reconstructed_lxy"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_lxy",b+"V0_L_reconstructed_lxy; lxy (beamspot, L creation vertex)(cm)",1000,0,100);
   histos_teff[b+"V0_L_reconstructed_decay_lxy"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_decay_lxy",b+"V0_L_reconstructed_decay_lxy; lxy (beamspot, L decay vertex)(cm)",200,0,100);
   histos_teff[b+"V0_L_reconstructed_dxy"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_dxy",b+"V0_L_reconstructed_dxy; dxy (beamspot, L)(cm)",400,-100,100);
   histos_teff[b+"V0_L_reconstructed_vz"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_vz",b+"V0_L_reconstructed_vz; vz(cm)",200,-100,100);
   histos_teff[b+"V0_L_reconstructed_eta"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_eta",b+"V0_L_reconstructed_eta; eta",100,-4,4);
   histos_teff[b+"V0_L_reconstructed_phi"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_phi",b+"V0_L_reconstructed_phi; phi(rad)",100,-4,4);
   histos_teff[b+"V0_L_reconstructed_status"] = dir_V0s_matched_L.make<TEfficiency>(b+"V0_L_reconstructed_status",b+"V0_L_reconstructed_status; status",10,0,10);

   //For the reconstructed S particles:
    TFileDirectory dir_LambdaKshortVertexFilter = m_fs->mkdir("LambdaKshortVertexFilter");
    TFileDirectory dir_LambdaKshortVertexFilter_S = dir_LambdaKshortVertexFilter.mkdir("S");
    histos_th1f[b+"h_LambdaKshortVertexFilter_S_mass"]= dir_LambdaKshortVertexFilter_S.make<TH1F>(b+"h_LambdaKshortVertexFilter_S_mass", b+"h_LambdaKshortVertexFilter_S_mass; S mass (GeV)", 2000, -100, 100);
    histos_th1f[b+"h_LambdaKshortVertexFilter_Sn_mass"]= dir_LambdaKshortVertexFilter_S.make<TH1F>(b+"h_LambdaKshortVertexFilter_Sn_mass", b+"h_LambdaKshortVertexFilter_Sn_mass; S+n mass (GeV)", 2000, -100, 100);
    histos_th1f[b+"h_LambdaKshortVertexFilter_S_mass_with_displacement_larger_1p5"]= dir_LambdaKshortVertexFilter_S.make<TH1F>(b+"h_LambdaKshortVertexFilter_S_mass_with_displacement_larger_1p5", b+"h_LambdaKshortVertexFilter_S_mass_with_displacement_larger_1p5; S mass (GeV)", 2000, -100, 100);
    histos_th2f[b+"h_LambdaKshortVertexFilter_S_vx_vy"]= dir_LambdaKshortVertexFilter_S.make<TH2F>(b+"h_LambdaKshortVertexFilter_S_vx_vy", b+"h_LambdaKshortVertexFilter_S_vx_vy; S vx (cm); S vy (cm)", 2000, -100, 100, 2000, -100, 100);

}


void Analyzer_Steven::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


//  cout << "------------------------------------------------------------------------------------------------------------------------------------------------------------------"<< endl; 
 
  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

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
  const reco::BeamSpot* theBeamSpot = h_bs.product();
  math::XYZPoint referencePos(theBeamSpot->position());
  if(h_bs.isValid()){ 	

	//beamspot
	double bx_x = h_bs->x0(); 
	double bx_y = h_bs->y0(); 
	double bx_z = h_bs->z0();
	beamspot.SetXYZ(bx_x,bx_y,bx_z);
	
  }

 
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
			const reco::Track * matchingTrack = nullptr;
			for(unsigned int t = 0; t < h_generalTracks->size(); ++t){//loop over all the tracks and find the minimal deltaR
				Double_t track_phi = h_generalTracks->at(t).phi();
				Double_t track_eta = h_generalTracks->at(t).eta();
				Double_t Delta_R_V0_daug0_track = pow(pow(V0_daug_phi-track_phi,2)+pow(V0_daug_eta-track_eta,2),0.5);
				if(Delta_R_V0_daug0_track < Delta_R_V0_daug0_track_min) {
					Delta_R_V0_daug0_track_min = Delta_R_V0_daug0_track;
					matchingTrack = &h_generalTracks->at(t);
				}
				histos_th1f[b+"h_tracks_deltaR"]->Fill(Delta_R_V0_daug0_track);
				
			}
			histos_th1f[b+"h_tracks_deltaR_min"]->Fill(Delta_R_V0_daug0_track_min);
			if(Delta_R_V0_daug0_track_min<0.1){//matched a track with a daughter of a Ks or Lambda
				if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId()==-1020000020){//and the daughter is coming from an antiS
					histos_teff[b+"tracks_reco_eff_daughters_antiS"]->Fill(true,1);						
					histos_teff[b+"tracks_reco_eff_daughters_antiS_pt"]->Fill(true,V0_daug->pt());
					histos_teff[b+"tracks_reco_eff_daughters_antiS_vxy"]->Fill(true,V0_daug_vxy);						
					cout << "found an antiS V0 daughter track" << endl;
				

				 	//now check if this specific track survives the track cuts in the V0Fitter:
					   // fill vectors of TransientTracks and TrackRefs after applying preselection cuts
					      const reco::Track * tmpTrack = matchingTrack;
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
					      }

				}
				else{
					histos_teff[b+"tracks_reco_eff_daughters_antiS"]->Fill(true,0);						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt"]->Fill(true,V0_daug->pt());						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_vxy"]->Fill(true,V0_daug_vxy);						
					if(V0_daug_vxy<3.5) histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5"]->Fill(true,V0_daug->pt());						
				}

			}
			else{//no matched track found
				if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId()==-1020000020){
					histos_teff[b+"tracks_reco_eff_daughters_antiS"]->Fill(false,1);						
					histos_teff[b+"tracks_reco_eff_daughters_antiS_pt"]->Fill(false,V0_daug->pt());						
					histos_teff[b+"tracks_reco_eff_daughters_antiS_vxy"]->Fill(false,V0_daug_vxy);						
				}
				else{
					histos_teff[b+"tracks_reco_eff_daughters_antiS"]->Fill(false,0);						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt"]->Fill(false,V0_daug->pt());						
					histos_teff[b+"tracks_reco_eff_no_daughters_antiS_vxy"]->Fill(false,V0_daug_vxy);						
					if(V0_daug_vxy<3.5) histos_teff[b+"tracks_reco_eff_no_daughters_antiS_pt_vxy_smaller_3p5"]->Fill(false,V0_daug->pt());						
				}
			}

		}//loop over daughters of the Lambda or Ks
      }//for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i)
 }//if(h_genParticles_SIM_GEANT.isValid())
 

  //the kinematics of the GenParticlesPlusGeant Ks
  if(h_genParticles_SIM_GEANT.isValid())
  {
    for(unsigned int k = 0; k < h_genParticles_SIM_GEANT->size(); k++){
        if(fabs(h_genParticles_SIM_GEANT->at(k).pdgId())==310){
			histos_th1f[b+"GEN_Ks_pt"]->Fill(h_genParticles_SIM_GEANT->at(k).pt());				
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
        if(fabs(h_genParticles_SIM_GEANT->at(l).pdgId())==3112){
			histos_th1f[b+"GEN_L_pt"]->Fill(h_genParticles_SIM_GEANT->at(l).pt());				
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
		Double_t GEN_Ks_lxy = lxy(beamspot,GEN_Ks_vxyz);

		Double_t GEN_Ks_decay_lxy = 0;
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2){
			Double_t GEN_Ks_decay_vx = h_genParticles_SIM_GEANT->at(i).daughter(0)->vx();
			Double_t GEN_Ks_decay_vy = h_genParticles_SIM_GEANT->at(i).daughter(0)->vy();
			Double_t GEN_Ks_decay_vz = h_genParticles_SIM_GEANT->at(i).daughter(0)->vz();

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
				bool matched = false;
                                if(deltaR_Ks_GEN_RECO<0.1){
					//plot the properties of the Ks which are not daughters of the antiS and which get reconstructed
					if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId() != -1020000020){
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_vxy"]->Fill(GEN_Ks_vxy);
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_lxy"]->Fill(GEN_Ks_lxy);
						histos_th1f[b+"V0s_Ks_reconstructed_GEN_decay_lxy"]->Fill(GEN_Ks_decay_lxy);
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
						histos_th1f[b+"V0s_Ks_daughterS_reconstructed_GEN_decay_lxy"]->Fill(GEN_Ks_decay_lxy);
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
						histos_th1f[b+"V0s_Ks_non_reconstructed_GEN_decay_lxy"]->Fill(GEN_Ks_decay_lxy);
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
						histos_th1f[b+"V0s_Ks_daughterS_non_reconstructed_GEN_decay_lxy"]->Fill(GEN_Ks_decay_lxy);
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
				histos_teff[b+"V0_Ks_reconstructed_decay_lxy"]->Fill(matched,GEN_Ks_decay_lxy);
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
		Double_t GEN_L_lxy = lxy(beamspot,GEN_L_vxyz);

		Double_t GEN_L_decay_lxy = 0;
		if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2){
			Double_t GEN_L_decay_vx = h_genParticles_SIM_GEANT->at(i).daughter(0)->vx();
			Double_t GEN_L_decay_vy = h_genParticles_SIM_GEANT->at(i).daughter(0)->vy();
			Double_t GEN_L_decay_vz = h_genParticles_SIM_GEANT->at(i).daughter(0)->vz();

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
				bool matched = false;
                                if(deltaR_L_GEN_RECO<0.1){
					//plot the properties of the L which are not daughters of the antiS and which get reconstructed
					if(h_genParticles_SIM_GEANT->at(i).mother()->pdgId() != -1020000020){
						histos_th1f[b+"V0s_L_reconstructed_GEN_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
						histos_th1f[b+"V0s_L_reconstructed_GEN_p"]->Fill(h_genParticles_SIM_GEANT->at(i).p());
						histos_th1f[b+"V0s_L_reconstructed_GEN_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
						histos_th1f[b+"V0s_L_reconstructed_GEN_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
						histos_th1f[b+"V0s_L_reconstructed_GEN_vxy"]->Fill(GEN_L_vxy);
						histos_th1f[b+"V0s_L_reconstructed_GEN_lxy"]->Fill(GEN_L_lxy);
						histos_th1f[b+"V0s_L_reconstructed_GEN_decay_lxy"]->Fill(GEN_L_decay_lxy);
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
						histos_th1f[b+"V0s_L_daughterS_reconstructed_GEN_decay_lxy"]->Fill(GEN_L_decay_lxy);
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
						histos_th1f[b+"V0s_L_non_reconstructed_GEN_decay_lxy"]->Fill(GEN_L_decay_lxy);
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
						histos_th1f[b+"V0s_L_daughterS_non_reconstructed_GEN_decay_lxy"]->Fill(GEN_L_decay_lxy);
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
				histos_teff[b+"V0_L_reconstructed_decay_lxy"]->Fill(matched,GEN_L_decay_lxy);
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
			
			if(lxy_S_b > 1.5){		
				histos_th1f[b+"h_LambdaKshortVertexFilter_S_mass_with_displacement_larger_1p5"]->Fill(S_mass);
			}
			
	}//for (unsigned int i = 0; i < n_sCands; ++i)
  }//if(h_sCands.isValid())

} //end of analyzer

double Analyzer_Steven::openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2){
  double opening_angle = TMath::ACos((momentum1.Dot(momentum2))/(pow(momentum1.Mag2()*momentum2.Mag2(),0.5)));
  return opening_angle;
}

double Analyzer_Steven::deltaR(double phi1, double eta1, double phi2, double eta2){
	double deltaPhi = reco::deltaPhi(phi1,phi2);
	double deltaEta = eta1-eta2;
	return pow(deltaPhi*deltaPhi+deltaEta*deltaEta,0.5);
}


double Analyzer_Steven::lxy(TVector3 v1, TVector3 v2){
	double x1 = v1.X();
	double x2 = v2.X();
	double y1 = v1.Y();
	double y2 = v2.Y();
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}


TVector3 Analyzer_Steven::PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point){
   //first move the vector along the line to the starting point of Point_line
   double normalise = sqrt(Vector_along_line.X()*Vector_along_line.X()+Vector_along_line.Y()*Vector_along_line.Y()+Vector_along_line.Z()*Vector_along_line.Z());
   TVector3 n(Vector_along_line.X()/normalise,Vector_along_line.Y()/normalise,Vector_along_line.Z()/normalise);
   TVector3 a = Point_line;
   TVector3 p = Point;

   //see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line (Vector formulation)
   TVector3 vector_PCA = (a-p)-((a-p)*n)*n;
   return vector_PCA ;
}

double Analyzer_Steven::dxy_signed_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in){

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

void Analyzer_Steven::endJob()
{
}

Analyzer_Steven::~Analyzer_Steven()
{
}


DEFINE_FWK_MODULE(Analyzer_Steven);
