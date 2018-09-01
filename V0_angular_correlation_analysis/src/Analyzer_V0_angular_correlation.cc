#include "SexaQAnalysis/V0_angular_correlation_analysis/interface/Analyzer_V0_angular_correlation.h"
#include <typeinfo>


Analyzer_V0_angular_correlation::Analyzer_V0_angular_correlation(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  //m_rCandsTag(pset.getParameter<edm::InputTag>("resonCandidates")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_KshortsTag(pset.getParameter<edm::InputTag>("kshortCollection")),
  m_LambdasTag(pset.getParameter<edm::InputTag>("lambdaCollection")),
  m_LambdasLambdaKshortFilterTag(pset.getParameter<edm::InputTag>("lambdaCollectionlambdaKshortFilter")),
  m_KshortsLambdaKshortFilterTag(pset.getParameter<edm::InputTag>("kshortCollectionlambdaKshortFilter")),
  m_sCollectionMassFilterTag(pset.getParameter<edm::InputTag>("sCollectionMassFilter")),
  m_rCollectionMassFilterTag(pset.getParameter<edm::InputTag>("rCollectionMassFilter")),
  m_nPVsTag(pset.getParameter<edm::InputTag>("nPVsCollection")),
  m_nelectronsTag(pset.getParameter<edm::InputTag>("nelectronsCollection")),
  m_njetsTag(pset.getParameter<edm::InputTag>("njetsCollection")),
  m_nkshortsTag(pset.getParameter<edm::InputTag>("nkshortsCollection")),
  m_nlambdasTag(pset.getParameter<edm::InputTag>("nlambdasCollection")),
  m_nmuonsTag(pset.getParameter<edm::InputTag>("nmuonsCollection")),
  m_ntracksTag(pset.getParameter<edm::InputTag>("ntracksCollection")),
  m_HTTag(pset.getParameter<edm::InputTag>("HTCollection")),
  m_TKHTTag(pset.getParameter<edm::InputTag>("TKHTCollection")),
  m_TwoTopJetsTag(pset.getParameter<edm::InputTag>("TwoTopJetsCollection")),
  m_METTag(pset.getParameter<edm::InputTag>("METCollection")),
  m_TKMETTag(pset.getParameter<edm::InputTag>("TKMETCollection")),
  //m_rCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_rCandsTag)),
  //m_sCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_sCandsTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_KshortsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_KshortsTag)),
  m_LambdasToken(consumes<vector<reco::VertexCompositeCandidate> >(m_LambdasTag)),
  m_LambdasLambdaKshortFilterToken(consumes<edm::PtrVector<reco::Candidate >>(m_LambdasLambdaKshortFilterTag)),
  m_KshortsLambdaKshortFilterToken(consumes<edm::PtrVector<reco::Candidate >>(m_KshortsLambdaKshortFilterTag)),
  m_sCollectionMassFilterToken(consumes<vector<reco::VertexCompositePtrCandidate>>(m_sCollectionMassFilterTag)),
  m_rCollectionMassFilterToken(consumes<vector<reco::VertexCompositePtrCandidate>>(m_rCollectionMassFilterTag)),
  m_nPVsToken(consumes<vector<int>>(m_nPVsTag)),
  m_nelectronsToken(consumes<vector<int>>(m_nelectronsTag)),
  m_njetsToken(consumes<vector<int>>(m_njetsTag)),
  m_nkshortsToken(consumes<vector<int>>(m_nkshortsTag)),
  m_nlambdasToken(consumes<vector<int>>(m_nlambdasTag)),
  m_nmuonsToken(consumes<vector<int>>(m_nmuonsTag)),
  m_ntracksToken(consumes<vector<int>>(m_ntracksTag)),
  m_HTToken(consumes<vector<double>>(m_HTTag)),
  m_TKHTToken(consumes<vector<double>>(m_TKHTTag)),
  m_TwoTopJetsToken(consumes<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>>(m_TwoTopJetsTag)),
  m_METToken(consumes<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>>(m_METTag)),
  m_TKMETToken(consumes<vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag>>>(m_TKMETTag))


{

}


void Analyzer_V0_angular_correlation::beginJob() {
  
   //for scalar
    histos_th1f[b+"h_h_nPVs"]= m_fs->make<TH1F>(b+"h_h_nPVs",b+"h_h_nPVs; nPVs",1000,0,1000);
    histos_th1f[b+"h_h_nelectrons"]= m_fs->make<TH1F>(b+"h_h_nelectrons",b+"h_h_nelectrons; nelectrons",1000,0,1000);
    histos_th1f[b+"h_h_njets"]= m_fs->make<TH1F>(b+"h_h_njets",b+"h_h_njets; njets",1000,0,1000);
    histos_th1f[b+"h_h_nkshorts"]= m_fs->make<TH1F>(b+"h_h_nkshorts",b+"h_h_nkshorts; nkshorts",1000,0,1000);
    histos_th1f[b+"h_h_nlambdas"]= m_fs->make<TH1F>(b+"h_h_nlambdas",b+"h_h_nlambdas; nlambdas",1000,0,1000);
    histos_th1f[b+"h_h_nmuons"]= m_fs->make<TH1F>(b+"h_h_nmuons",b+"h_h_nmuons; nmuons",1000,0,1000);
    histos_th1f[b+"h_h_ntracks"]= m_fs->make<TH1F>(b+"h_h_ntracks",b+"h_h_ntracks; ntracks",1000,0,1000);
  
    histos_th1f[b+"h_h_TKMET_pT"]= m_fs->make<TH1F>(b+"h_h_TKMET_pT",b+"h_h_TKMET_pT; ",2000,-1000,1000);
    histos_th1f[b+"h_h_MET_pT"]= m_fs->make<TH1F>(b+"h_h_MET_pT",b+"h_h_MET_pT; ",2000,-1000,1000);
    histos_th1f[b+"h_h_TwoTopJets_pT_0"]= m_fs->make<TH1F>(b+"h_h_TwoTopJets_pT_0",b+"h_h_TwoTopJets_pT_0; ",2000,-1000,1000);
    histos_th1f[b+"h_h_TwoTopJets_pT_1"]= m_fs->make<TH1F>(b+"h_h_TwoTopJets_pT_1",b+"h_h_TwoTopJets_pT_1; ",2000,-1000,1000);
    histos_th1f[b+"h_h_HT"]= m_fs->make<TH1F>(b+"h_h_HT",b+"h_h_HT; ",5000,0,5000);
    histos_th1f[b+"h_h_TKHT"]= m_fs->make<TH1F>(b+"h_h_TKHT",b+"h_h_TKHT; ",5000,0,5000);

    
  
//  histos_th1f[b+"nPV"]      = m_fs->make<TH1F>(b+"nPV",     a+" Number of PV; #PVs",60,0.,60);
   //masses from the original kshort and lambda collection
   histos_th1f[b+"h_Lambda_mass"] = m_fs->make<TH1F>(b+"h_Lambda_mass",b+"h_Lambda_mass; Lambda mass from original V0 collection",4000,0,20);   
   histos_th1f[b+"h_Kshort_mass"] = m_fs->make<TH1F>(b+"h_Kshort_mass",b+"h_Kshort_mass; Kshort mass from original V0 collection",4000,0,20);   
   
   //masses of Ks and Lambda after the LambdaKshortVertexFilter
   histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass_all"] = m_fs->make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass_all", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass_all; Lambda mass after LambdaKshortVertexFilter",4000, 0,20);
   histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass_all"] = m_fs->make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass_all", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass_all; Kshort mass after LambdaKshortVertexFilter",4000, 0,20);
   histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass"] = m_fs->make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass; Lambda mass after LambdaKshortVertexFilter for S candidate",4000, 0,20);
   histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass"] = m_fs->make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass; Kshort mass after LambdaKshortVertexFilter for S candidate",4000, 0,20);
   histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass; Lambda mass after LambdaKshortVertexFilter for anti S candidate",4000, 0,20);
   histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass; Kshort mass after LambdaKshortVertexFilter for anti S candidate",4000, 0,20);
  
  //other kinematic quantities after the LambdaKshortVertexFilter 
	//Lambda
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_pt_all"]= m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_pt_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_pt_all; Lambda pt after the LambdaKshortVertexFilter; pT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_p_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_p_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_p_all; Lambda p after the LambdaKshortVertexFilter; p (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_energy_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_energy_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_energy_all; Lambda energy after the LambdaKshortVertexFilter; E (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_et_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_et_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_et_all; Lambda transversal E after the LambdaKshortVertexFilter; eT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_mt_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_mt_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_mt_all; Lambda transversal mass after the LambdaKshortVertexFilter; mT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_phi_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_phi_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_phi_all; Lambda phi after the LambdaKshortVertexFilter; phi (rad) ",8000,-4,4);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_theta_all"]= m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_theta_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_theta_all; Lambda theta after the LambdaKshortVertexFilter; theta (rad) ",800,-4,4);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_eta_all"]= m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_eta_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_eta_all; Lambda eta after the LambdaKshortVertexFilter; eta ",800,-10,10);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_all; Lambda vx after the LambdaKshortVertexFilter; vx (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_all"]= m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_all; Lambda vy after the LambdaKshortVertexFilter; vy (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vz_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vz_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vz_all; Lambda vz after the LambdaKshortVertexFilter; vz (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vertexNormalizedChi2_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vertexNormalizedChi2_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vertexNormalizedChi2_all; Lambda vertex normalised Chi2 after the LambdaKshortVertexFilter ",1000, 0,100);
	//Kshort
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_pt_all"]= m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_pt_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_pt_all; Kshort pt after the KshortKshortVertexFilter; pT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_p_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_p_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_p_all; Kshort p after the LambdaKshortVertexFilter; p (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_energy_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_energy_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_energy_all; Kshort energy after the LambdaKshortVertexFilter; E (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_et_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_et_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_et_all; Kshort transversal E after the LambdaKshortVertexFilter; eT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_mt_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_mt_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_mt_all; Kshort transversal mass after the LambdaKshortVertexFilter; mT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_phi_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_phi_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_phi_all; Kshort phi after the LambdaKshortVertexFilter; phi (rad) ",8000,-4,4);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_theta_all"]= m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_theta_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_theta_all; Kshort theta after the LambdaKshortVertexFilter; theta (rad) ",800,-4,4);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_eta_all"]= m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_eta_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_eta_all; Kshort eta after the LambdaKshortVertexFilter; eta ",800,-10,10);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_all; Kshort vx after the LambdaKshortVertexFilter; vx (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_all"]= m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_all; Kshort vy after the LambdaKshortVertexFilter; vy (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vz_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vz_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vz_all; Kshort vz after the LambdaKshortVertexFilter; vz (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vertexNormalizedChi2_all"] = m_fs->make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vertexNormalizedChi2_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vertexNormalizedChi2_all; Kshort vertex normalised Chi2 after the LambdaKshortVertexFilter ",1000, 0,100);
   
   //for angular correlation between Ks and L0 produced in the filters
   histos_th1f[b+"h_L0_Ks_delta_phi_Filters_Collection"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_phi_Filters_Collection",b+"h_L0_Ks_delta_phi_Filters_Collection; delta phi",1000,-4,4); 
   histos_th1f[b+"h_L0_Ks_delta_eta_Filters_Collection"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_eta_Filters_Collection",b+"h_L0_Ks_delta_eta_Filters_Collection; delta eta",4000,-10,10); 
   histos_th1f[b+"h_L0_Ks_delta_R_Filters_Collection"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_R_Filters_Collection",b+"h_L0_Ks_delta_R_Filters_Collection; delta R",3000,0,30); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection"]= m_fs->make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection",b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection; delta phi; delta_eta",100,-4,4, 100, -10, 10); 

   //for angular correlation between Ks and L0 
   histos_th1f[b+"h_L0_Ks_delta_phi"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_phi",b+"h_L0_Ks_delta_phi; delta phi",1000,-4,4); 
   histos_th1f[b+"h_L0_Ks_delta_eta"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_eta",b+"h_L0_Ks_delta_eta; delta eta",4000,-10,10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta"]= m_fs->make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta",b+"h_L0_Ks_delta_phi_delta_eta; delta phi; delta_eta",100,-4,4, 100, -10, 10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_no_cent"]= m_fs->make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta_no_cent",b+"h_L0_Ks_delta_phi_delta_eta_no_cent; delta phi; delta_eta",1000,-4,4, 4000, -10, 10); 
  
   //for angular correlation between Ks and L0 from different events
   histos_th1f[b+"h_L0_Ks_delta_phi_prev_and_current"]= m_fs->make<TH1F>(b+"h_L0_delta_phi_prev_and_current",b+"h_L0_delta_phi_prev_and_current; delta phi prev event and current event;",1000,-4,4); 
   histos_th1f[b+"h_L0_Ks_delta_eta_prev_and_current"]= m_fs->make<TH1F>(b+"h_L0_delta_eta_prev_and_current",b+"h_L0_delta_eta_prev_and_current; delta eta prev event and current event;",4000,-10,10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]= m_fs->make<TH2F>(b+"h_L0_delta_phi_delta_eta_prev_and_current",b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current; delta phi prev event and current event; delta_eta previous event and current event",1000,-4,4, 4000, -10, 10); 
   
   //for angular coorelation between the daughters of the S candidates
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_phi"]= m_fs->make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_phi",b+"h_Sdaughters_L0_Ks_delta_phi; delta phi;",1000,-4,4); 
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_eta"]= m_fs->make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_eta",b+"h_Sdaughters_L0_Ks_delta_eta; delta eta;",4000,-10,10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta"]= m_fs->make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta; delta phi; delta eta",1000,-4,4, 4000, -10, 10); 
  
   //for the s candidates after the massfilter
    
   histos_th1f[b+"h_s_candidates_mass_after_massFilter"]= m_fs->make<TH1F>(b+"h_s_candidates_mass_after_massFilter",b+"h_s_candidates_mass_after_massFilter; S candidate mass (GeV);",40000,-200,200); 
   histos_th1f[b+"h_s_candidates_pt_after_massFilter"] =  m_fs->make<TH1F>(b+"h_s_candidates_pt_after_massFilter",b+"h_s_candidates_pt_after_massFilter; S candidate pt (GeV);",40000,-200,200);
   histos_th1f[b+"h_s_candidates_p_after_massFilter"] =  m_fs->make<TH1F>(b+"h_s_candidates_p_after_massFilter",b+"h_s_candidates_p_after_massFilter; S candidate p (GeV);",40000,-200,200);
   histos_th1f[b+"h_s_candidates_energy_after_massFilter"]  =  m_fs->make<TH1F>(b+"h_s_candidates_energy_after_massFilter",b+"h_s_candidates_energy_after_massFilter; S candidate E (GeV);",20000,0,200);
   histos_th1f[b+"h_s_candidates_et_after_massFilter"] =  m_fs->make<TH1F>(b+"h_s_candidates_et_after_massFilter",b+"h_s_candidates_et_after_massFilter; S candidate Et (GeV);",20000,0,200);
   histos_th1f[b+"h_s_candidates_mt_after_massFilter"] = m_fs->make<TH1F>(b+"h_s_candidates_mt_after_massFilter",b+"h_s_candidates_mt_after_massFilter; S candidate mt (GeV);",40000,-200,200);
   histos_th1f[b+"h_s_candidates_phi_after_massFilter"] = m_fs->make<TH1F>(b+"h_s_candidates_phi_after_massFilter",b+"h_s_candidates_phi_after_massFilter; S candidate phi (rad);",8000,-4,4);
   histos_th1f[b+"h_s_candidates_theta_after_massFilter"] = m_fs->make<TH1F>(b+"h_s_candidates_theta_after_massFilter",b+"h_s_candidates_theta_after_massFilter; S candidate theta (rad);",8000,-4,4);
   histos_th1f[b+"h_s_candidates_eta_after_massFilter"] = m_fs->make<TH1F>(b+"h_s_candidates_eta_after_massFilter",b+"h_s_candidates_eta_after_massFilter; S candidate eta;",800,-10,10);
   histos_th1f[b+"h_s_candidates_vx_after_massFilter"] = m_fs->make<TH1F>(b+"h_s_candidates_vx_after_massFilter",b+"h_s_candidates_vx_after_massFilter; S candidate vx (cm);",2000,-100,100);
   histos_th1f[b+"h_s_candidates_vy_after_massFilter"] = m_fs->make<TH1F>(b+"h_s_candidates_vy_after_massFilter",b+"h_s_candidates_vy_after_massFilter; S candidate vy (cm);",2000,-100,100);
   histos_th1f[b+"h_s_candidates_vz_after_massFilter"] = m_fs->make<TH1F>(b+"h_s_candidates_vz_after_massFilter",b+"h_s_candidates_vz_after_massFilter; S candidate vz (cm);",2000,-100,100);
   histos_th1f[b+"h_s_candidates_vertexNormalizedChi2_after_massFilter"] = m_fs->make<TH1F>(b+"h_s_candidates_vertexNormalizedChi2_after_massFilter",b+"h_s_candidates_vertexNormalizedChi2_after_massFilter; S candidate vertex normalised Chi2;",1000,0,100);
		
   histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter"]= m_fs->make<TH1F>(b+"h_anti_s_candidates_mass_after_massFilter",b+"h_anti_s_candidates_mass_after_massFilter; anti S candidate mass (GeV);",40000,-200,200); 
   histos_th1f[b+"h_anti_s_candidates_pt_after_massFilter"] =  m_fs->make<TH1F>(b+"h_anti_s_candidates_pt_after_massFilter",b+"h_anti_s_candidates_pt_after_massFilter; anti S candidate pt (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_p_after_massFilter"] =  m_fs->make<TH1F>(b+"h_anti_s_candidates_p_after_massFilter",b+"h_anti_s_candidates_p_after_massFilter; anti S candidate p (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_energy_after_massFilter"]  =  m_fs->make<TH1F>(b+"h_anti_s_candidates_energy_after_massFilter",b+"h_anti_s_candidates_energy_after_massFilter; anti S candidate E (GeV);",20000,0,200);
   histos_th1f[b+"h_anti_s_candidates_et_after_massFilter"] =  m_fs->make<TH1F>(b+"h_anti_s_candidates_et_after_massFilter",b+"h_anti_s_candidates_et_after_massFilter; anti S candidate Et (GeV);",20000,0,200);
   histos_th1f[b+"h_anti_s_candidates_mt_after_massFilter"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_mt_after_massFilter",b+"h_anti_s_candidates_mt_after_massFilter; anti S candidate mt (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_phi_after_massFilter"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_phi_after_massFilter",b+"h_anti_s_candidates_phi_after_massFilter; anti S candidate phi (rad);",8000,-4,4);
   histos_th1f[b+"h_anti_s_candidates_theta_after_massFilter"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_theta_after_massFilter",b+"h_anti_s_candidates_theta_after_massFilter; anti S candidate theta (rad);",8000,-4,4);
   histos_th1f[b+"h_anti_s_candidates_eta_after_massFilter"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_eta_after_massFilter",b+"h_anti_s_candidates_eta_after_massFilter; anti S candidate eta;",800,-10,10);
   histos_th1f[b+"h_anti_s_candidates_vx_after_massFilter"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_vx_after_massFilter",b+"h_anti_s_candidates_vx_after_massFilter; anti S candidate vx (cm);",2000,-100,100);
   histos_th1f[b+"h_anti_s_candidates_vy_after_massFilter"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_vy_after_massFilter",b+"h_anti_s_candidates_vy_after_massFilter; anti S candidate vy (cm);",2000,-100,100);
   histos_th1f[b+"h_anti_s_candidates_vz_after_massFilter"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_vz_after_massFilter",b+"h_anti_s_candidates_vz_after_massFilter; anti S candidate vz (cm);",2000,-100,100);
   histos_th1f[b+"h_anti_s_candidates_vertexNormalizedChi2_after_massFilter"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_vertexNormalizedChi2_after_massFilter",b+"h_anti_s_candidates_vertexNormalizedChi2_after_massFilter; anti S candidate vertex normalised Chi2;",1000,0,100);
		

   //for the daughters of the s candidates after the massfilter
/*
   histos_th1f[b+"h_s_candidates_mass_after_massFilter_Lambda_mass"] = m_fs->make<TH1F>(b+"h_s_candidates_mass_after_massFilter_Lambda_mass",b+"h_s_candidates_mass_after_massFilter_Lambda_mass; Lambda mass in the S decay chain (GeV);",40000,-200,200);
   histos_th1f[b+"h_s_candidates_mass_after_massFilter_Kshort_mass"] = m_fs->make<TH1F>(b+"h_s_candidates_mass_after_massFilter_Kshort_mass",b+"h_s_candidates_mass_after_massFilter_Kshort_mass; Kshort mass in the S decay chain (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter_Lambda_mass"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_mass_after_massFilter_Lambda_mass",b+"h_anti_s_candidates_mass_after_massFilter_Lambda_mass; Lambda mass in the anti S decay chain (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter_Kshort_mass"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_mass_after_massFilter_Kshort_mass",b+"h_anti_s_candidates_mass_after_massFilter_Kshort_mass; Kshort mass in the anti S decay chain (GeV);",40000,-200,200);
*/
   //for the r candidates after the massfilter
    
   histos_th1f[b+"h_r_candidates_mass_after_massFilter"]= m_fs->make<TH1F>(b+"h_r_candidates_mass_after_massFilter",b+"h_r_candidates_mass_after_massFilter; r candidate mass (GeV);",40000,-200,200); 
   histos_th1f[b+"h_anti_r_candidates_mass_after_massFilter"]= m_fs->make<TH1F>(b+"h_anti_r_candidates_mass_after_massFilter",b+"h_anti_r_candidates_mass_after_massFilter; anti r candidate mass (GeV);",40000,-200,200); 
    
}


void Analyzer_V0_angular_correlation::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  //std::cout << "START LOOPING OVER THIS EVENT" << std::endl;


  // resonance candidates
  
  //edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_rCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
  //iEvent.getByToken(m_rCandsToken, h_rCands);

  // sexaquark candidates
  
  //edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_sCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
  //iEvent.getByToken(m_sCandsToken, h_sCands);
  
  //scalars
  edm::Handle <vector<int> > h_nPVs;
  iEvent.getByToken(m_nPVsToken, h_nPVs);

  edm::Handle <vector<int> > h_nelectrons;
  iEvent.getByToken(m_nelectronsToken, h_nelectrons);

  edm::Handle <vector<int> > h_njets;
  iEvent.getByToken(m_njetsToken, h_njets);
  
  edm::Handle <vector<int> > h_nkshorts;
  iEvent.getByToken(m_nkshortsToken, h_nkshorts);

  edm::Handle <vector<int> > h_nlambdas;
  iEvent.getByToken(m_nlambdasToken, h_nlambdas);

  edm::Handle <vector<int> > h_nmuons;
  iEvent.getByToken(m_nmuonsToken, h_nmuons);

  edm::Handle <vector<int> > h_ntracks;
  iEvent.getByToken(m_ntracksToken, h_ntracks);


  edm::Handle <vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > > h_TKMET;
  iEvent.getByToken(m_TKMETToken, h_TKMET);

  edm::Handle <vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >h_MET;
  iEvent.getByToken(m_METToken, h_MET);

  edm::Handle <vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > h_TwoTopJets;
  iEvent.getByToken(m_TwoTopJetsToken, h_TwoTopJets);

  edm::Handle <vector<double> > h_HT;
  iEvent.getByToken(m_HTToken, h_HT);

  edm::Handle <vector<double> > h_TKHT;
  iEvent.getByToken(m_TKHTToken, h_TKHT);

  //lambdaKshortVertexFilter sexaquark candidates
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);
  
  //reco Kshorts V0
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_Kshorts;
  iEvent.getByToken(m_KshortsToken, h_Kshorts);
  
  //reco Lambdas V0
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_Lambdas;
  iEvent.getByToken(m_LambdasToken, h_Lambdas);

  //filtered Lambdas
  edm::Handle <edm::PtrVector<reco::Candidate > > h_Lambdas_LambdaKshortFilter;
  iEvent.getByToken(m_LambdasLambdaKshortFilterToken, h_Lambdas_LambdaKshortFilter);

  //filtered Kshorts
  edm::Handle <edm::PtrVector<reco::Candidate > > h_Kshorts_LambdaKshortFilter;
  iEvent.getByToken(m_KshortsLambdaKshortFilterToken, h_Kshorts_LambdaKshortFilter);

  //s candidates stored in the sMassFilter
  edm::Handle <vector<reco::VertexCompositePtrCandidate> > h_s_MassFilter;
  iEvent.getByToken(m_sCollectionMassFilterToken, h_s_MassFilter);

  //r candidates stored in the sMassFilter
  edm::Handle <vector<reco::VertexCompositePtrCandidate> > h_r_MassFilter;
  iEvent.getByToken(m_rCollectionMassFilterToken, h_r_MassFilter);


  //do stuff with the scalars produced in the filters
  std::cout << "the scalars" << std::endl;

 if(h_nPVs.isValid() && h_nPVs->size() > 0)  		 	histos_th1f[b+"h_h_nPVs"]->Fill(h_nPVs->at(0));
 if(h_nelectrons.isValid() && h_nelectrons->size() > 0 )	histos_th1f[b+"h_h_nelectrons"]->Fill(h_nelectrons->at(0));
 if(h_njets.isValid() && h_njets->size() > 0)			histos_th1f[b+"h_h_njets"]->Fill(h_njets->at(0));
 if(h_nkshorts.isValid() && h_nkshorts->size() > 0)		histos_th1f[b+"h_h_nkshorts"]->Fill(h_nkshorts->at(0));
 if(h_nlambdas.isValid() && h_nlambdas->size() > 0)		histos_th1f[b+"h_h_nlambdas"]->Fill(h_nlambdas->at(0));
 if(h_nmuons.isValid() && h_nmuons->size() > 0)			histos_th1f[b+"h_h_nmuons"]->Fill(h_nmuons->at(0));
// if(h_ntracks.isValid() && h_ntracks->size() > 0)		histos_th1f[b+"h_h_tracks"]->Fill(h_ntracks->at(0));

 if(h_TKMET.isValid() && h_TKMET->size() > 0){			histos_th1f[b+"h_h_TKMET_pT"]->Fill(pow(h_TKMET->at(0).X()*h_TKMET->at(0).X()+h_TKMET->at(0).Y()*h_TKMET->at(0).Y(),0.5)); cout << "TKMET x and y component: " << h_TKMET->at(0).X() << " " << h_TKMET->at(0).Y() << endl; }
 if(h_MET.isValid() && h_MET->size() > 0)			histos_th1f[b+"h_h_MET_pT"]->Fill(h_MET->at(0).pt());
 if(h_TwoTopJets.isValid() && h_TwoTopJets->size() > 0)		histos_th1f[b+"h_h_TwoTopJets_pT_0"]->Fill(h_TwoTopJets->at(0).pt()); 
 if(h_TwoTopJets.isValid() && h_TwoTopJets->size() > 0)		histos_th1f[b+"h_h_TwoTopJets_pT_1"]->Fill(h_TwoTopJets->at(1).pt()); 
 if(h_HT.isValid() && h_HT->size() > 0)				histos_th1f[b+"h_h_HT"]->Fill(h_HT->at(0));
 if(h_TKHT.isValid() && h_TKHT->size() > 0)			histos_th1f[b+"h_h_TKHT"]->Fill(h_TKHT->at(0));

 
  std::cout << "LambdaKshortfilter" << std::endl;
  //do stuff with the collections produced in the filters
  if(h_Lambdas_LambdaKshortFilter.isValid() && h_Kshorts_LambdaKshortFilter.isValid()){
  std::cout << "------------------------------------------" << std::endl;
  for(unsigned int l = 0; l < h_Lambdas_LambdaKshortFilter->size(); ++l){
	for(unsigned int k = 0; k < h_Kshorts_LambdaKshortFilter->size(); ++k){
		double phi1 = (*h_Kshorts_LambdaKshortFilter)[k]->phi();
                double phi2 = (*h_Lambdas_LambdaKshortFilter)[l]->phi();
                double delta_phi = reco::deltaPhi(phi1, phi2);

                double eta1 = (*h_Kshorts_LambdaKshortFilter)[k]->eta();
                double eta2 = (*h_Lambdas_LambdaKshortFilter)[l]->eta();
                double delta_eta = eta1-eta2;

		double delta_R = sqrt(delta_phi*delta_phi+delta_eta*delta_eta);

                histos_th1f[b+"h_L0_Ks_delta_phi_Filters_Collection"]->Fill(delta_phi);
                histos_th1f[b+"h_L0_Ks_delta_eta_Filters_Collection"]->Fill(delta_eta);
                histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection"]->Fill(delta_phi,delta_eta);	
                histos_th1f[b+"h_L0_Ks_delta_R_Filters_Collection"]->Fill(delta_R);	
		
   	}
  }
  }
/*
  for(unsigned int l = 0; l < h_Lambdas_LambdaKshortFilter->size(); ++l){
        for(unsigned int k = 0; k < h_Kshorts_LambdaKshortFilter->size(); ++k){
	
     }
   }
*/

   //std::cout << "looping over the lambdas and kshorts together" << std::endl;
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STANDARD PLOT: ONLY THE CURRENT EVENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //h_lambdas and h_kshorts are the original lambda and kshort collections
  std::cout << "running on h_Lambdas" << std::endl;
  if(h_Lambdas.isValid()){
	for (unsigned int l = 0; l < h_Lambdas->size(); ++l) {
		histos_th1f[b+"h_Lambda_mass"]->Fill(h_Lambdas->at(l).mass());
	}
  }


  std::cout << "running on h_Kshorts" << std::endl;
  if(h_Kshorts.isValid()){
	for (unsigned int k = 0; k < h_Kshorts->size(); ++k) {
		histos_th1f[b+"h_Kshort_mass"]->Fill(h_Kshorts->at(k).mass());
	}
  }

  std::cout << "the lambdas and kshorts" << std::endl;
 if(h_Lambdas.isValid() && h_Kshorts.isValid()){
 for (unsigned int l = 0; l < h_Lambdas->size(); ++l) { //loop over reco V0 Lambdas
	for (unsigned int k = 0; k < h_Kshorts->size(); ++k) { //loop over reco Kshorts	
		//to check the sharp peak at 0 check if the kshort and lambdas are exactly the same sometimes
		//std::cout << "---------------------" << std::endl;
		//std::cout << setprecision(15) << "kshort momenta: " << h_Kshorts->at(k).px() << " " << h_Kshorts->at(k).py() << " " << h_Kshorts->at(k).pz() << std::endl;
		//std::cout << setprecision(15) << "lambda momenta: " << h_Lambdas->at(l).px() << " " << h_Lambdas->at(l).py() << " " << h_Lambdas->at(l).pz() << std::endl;
		//for debugging: check that the lambda and kshort which you save are really different
/*		 if(h_Lambdas->at(l).px() == h_Kshorts->at(k).px()){
			std::cout << "Lambda number " << l << " momenta: "<<  h_Lambdas->at(l).px()  << ", " << h_Lambdas->at(l).py() << ", " <<  h_Lambdas->at(l).pz() << std::endl;
			std::cout << "Kshort number " << k << " momenta: "<<  h_Kshorts->at(k).px()  << ", " << h_Kshorts->at(k).py() << ", " <<  h_Kshorts->at(k).pz() << std::endl;
		 }
*/
  		double phi1 = h_Kshorts->at(k).phi();
		double phi2 = h_Lambdas->at(l).phi();
		double delta_phi = reco::deltaPhi(phi1, phi2);
				
		double eta1 = h_Kshorts->at(k).eta();
		double eta2 = h_Lambdas->at(l).eta();
		double delta_eta = eta1-eta2;

		histos_th1f[b+"h_L0_Ks_delta_phi"]->Fill(delta_phi);
		histos_th1f[b+"h_L0_Ks_delta_eta"]->Fill(delta_eta);
		histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta"]->Fill(delta_phi,delta_eta);
		if( abs(delta_phi) > 0.1 || abs(delta_eta) > 0.1 )histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_no_cent"]->Fill(delta_phi,delta_eta);
	}
 }
 }


  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COMPARE THE PREVIOUS EVENT WITH THE CURRENT EVENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  bool emptyVectors = false;
  //std::cout << "---------------------" << std::endl;
  //std::cout << "going to check if the vectors are empty" << std::endl;
  if(v_L0_phi.empty() && v_L0_eta.empty() && v_Ks_phi.empty() && v_Ks_eta.empty()) emptyVectors = true;//this means they have been cleared in the previous event or it is the first event

 //if(emptyVectors) std::cout << "vectors empty" << std::endl;  
 //else std::cout << "vectors not empty" << std::endl;  
  
   //std::cout << "getting the size of the lambda and kshort collections" << std::endl;
  if(h_Lambdas.isValid()){
  for (unsigned int i = 0; i < h_Lambdas->size(); ++i) { //loop over reco V0 Lambdas
  	if(emptyVectors){//just save the data if there is no data in these vectrors, make the histogram in the next event
		v_L0_phi.push_back(h_Lambdas->at(i).phi());	
		v_L0_eta.push_back(h_Lambdas->at(i).eta());	
	} 
	else{//now make the histo, also using the data from the previous event
		for(unsigned int i_prev = 0; i_prev < v_Ks_phi.size(); i_prev++){//loop over all the data from the previous event
			double phi1 = v_Ks_phi[i_prev];
			double phi2 = h_Lambdas->at(i).phi();
			double delta_phi = reco::deltaPhi(phi1, phi2);
			
			double eta1 = v_Ks_eta[i_prev];
			double eta2 = h_Lambdas->at(i).eta();
			double delta_eta = eta1-eta2;
			histos_th1f[b+"h_L0_Ks_delta_phi_prev_and_current"]->Fill(delta_phi);
			histos_th1f[b+"h_L0_Ks_delta_eta_prev_and_current"]->Fill(delta_eta);
			histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]->Fill(delta_phi,delta_eta);
		}
	} 
  }//end loop over lambda
  }
  //std::cout << "loop over Lambdas done" << std::endl; 
  
  if(h_Kshorts.isValid()){
  for (unsigned int i = 0; i < h_Kshorts->size(); ++i) { //loop over reco Kshorts
 	 if(emptyVectors){
                v_Ks_phi.push_back(h_Kshorts->at(i).phi());
                v_Ks_eta.push_back(h_Kshorts->at(i).eta());
         } 
 
	 else{//now make the histo, also using the data from the previous event
		for(unsigned int i_prev = 0; i_prev < v_L0_phi.size(); i_prev++){//loop over all the data from the previous event
			double phi1 = v_L0_phi[i_prev];
			double phi2 = h_Kshorts->at(i).phi();
			double delta_phi = reco::deltaPhi(phi1, phi2);
			
			double eta1 = v_L0_eta[i_prev];
			double eta2 = h_Kshorts->at(i).eta();
			double delta_eta = eta1-eta2;

			histos_th1f[b+"h_L0_Ks_delta_phi_prev_and_current"]->Fill(delta_phi);
			histos_th1f[b+"h_L0_Ks_delta_eta_prev_and_current"]->Fill(delta_eta);
			histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]->Fill(delta_phi,delta_eta);
		}
	}
  }// end loop over Kshort
 }

 // std::cout << "loop over Kshorts done" << std::endl; 

  if(!emptyVectors){
	v_L0_phi.clear();
	v_L0_eta.clear();
	v_Ks_phi.clear();
	v_Ks_eta.clear();
  } 
 // std::cout << "vectors cleared" << std::endl; 


 //!!!!!!!!!!!!!!!!!!!!!!!!!!!make the same distribution but for the S candidates (if they are there, they are only there if they ran through the 2nd filter)!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  std::cout << "the scands" << std::endl;
  if(h_sCands.isValid()) {
	unsigned int n_sCands = h_sCands->size();
	for (unsigned int i = 0; i < n_sCands; ++i) { //loop over reco Kshorts	
			double phi1 = h_sCands->at(i).daughter(0)->phi();
			double phi2 = h_sCands->at(i).daughter(1)->phi();
			double delta_phi = reco::deltaPhi(phi1, phi2);
			
			double eta1 = h_sCands->at(i).daughter(0)->eta();
			double eta2 = h_sCands->at(i).daughter(1)->eta();
			double delta_eta = eta1-eta2;
   			
			histos_th1f[b+"h_Sdaughters_L0_Ks_delta_phi"]->Fill(delta_phi);
                        histos_th1f[b+"h_Sdaughters_L0_Ks_delta_eta"]->Fill(delta_eta);
                        histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta"]->Fill(delta_phi,delta_eta);
			
			//masses of the daughters
			histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass_all"]->Fill(h_sCands->at(i).daughter(0)->mass());
                        histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass_all"]->Fill(h_sCands->at(i).daughter(1)->mass());	
			if(h_sCands->at(i).charge()==1){
			              histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass"]->Fill(h_sCands->at(i).daughter(0)->mass());
			              histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass"]->Fill(h_sCands->at(i).daughter(1)->mass());
			}
			if(h_sCands->at(i).charge()==-1){
			              histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass"]->Fill(h_sCands->at(i).daughter(0)->mass());
			              histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass"]->Fill(h_sCands->at(i).daughter(1)->mass());
			}
			//other kinematic variables of the daughters
			//Lambda
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_pt_all"]->Fill(h_sCands->at(i).daughter(0)->pt());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_p_all"]->Fill(h_sCands->at(i).daughter(0)->p());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_energy_all"]->Fill(h_sCands->at(i).daughter(0)->energy());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_et_all"]->Fill(h_sCands->at(i).daughter(0)->et());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_mt_all"]->Fill(h_sCands->at(i).daughter(0)->mt());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_phi_all"]->Fill(h_sCands->at(i).daughter(0)->phi());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_theta_all"]->Fill(h_sCands->at(i).daughter(0)->theta());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_eta_all"]->Fill(h_sCands->at(i).daughter(0)->eta());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_all"]->Fill(h_sCands->at(i).daughter(0)->vx());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_all"]->Fill(h_sCands->at(i).daughter(0)->vy());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vz_all"]->Fill(h_sCands->at(i).daughter(0)->vz());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vertexNormalizedChi2_all"]->Fill(h_sCands->at(i).daughter(0)->vertexNormalizedChi2());
			//Kshort
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_pt_all"]->Fill(h_sCands->at(i).daughter(1)->pt());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_p_all"]->Fill(h_sCands->at(i).daughter(1)->p());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_energy_all"]->Fill(h_sCands->at(i).daughter(1)->energy());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_et_all"]->Fill(h_sCands->at(i).daughter(1)->et());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_mt_all"]->Fill(h_sCands->at(i).daughter(1)->mt());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_phi_all"]->Fill(h_sCands->at(i).daughter(1)->phi());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_theta_all"]->Fill(h_sCands->at(i).daughter(1)->theta());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_eta_all"]->Fill(h_sCands->at(i).daughter(1)->eta());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_all"]->Fill(h_sCands->at(i).daughter(1)->vx());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_all"]->Fill(h_sCands->at(i).daughter(1)->vy());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vz_all"]->Fill(h_sCands->at(i).daughter(1)->vz());
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vertexNormalizedChi2_all"]->Fill(h_sCands->at(i).daughter(1)->vertexNormalizedChi2());
		}

  }//end of sCands present

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!make some distributions for the S and r candidates
  std::cout << "the rMassFilter" << std::endl;
  if((h_r_MassFilter).isValid()){
  	for (unsigned int r = 0; r < (*h_r_MassFilter).size(); ++r) {
		if((*h_r_MassFilter)[r].charge() == 1) histos_th1f[b+"h_r_candidates_mass_after_massFilter"]->Fill((*h_r_MassFilter)[r].mass());
		if((*h_r_MassFilter)[r].charge() == -1) histos_th1f[b+"h_anti_r_candidates_mass_after_massFilter"]->Fill((*h_r_MassFilter)[r].mass());
	}
  }
 
  std::cout << "the sMassFilter" << std::endl;
  if((h_s_MassFilter).isValid()){
  	for (unsigned int s = 0; s < (*h_s_MassFilter).size(); ++s) {
		
		//mass plots of the S candidates. The charge is the charge of the proton in the decay chain, so a antiLambda would give an antiproton
		std::cout << "getting the charge and the mass" << std::endl;
		if((*h_s_MassFilter)[s].charge()==1){
			histos_th1f[b+"h_s_candidates_mass_after_massFilter"]->Fill((*h_s_MassFilter)[s].mass());
			histos_th1f[b+"h_s_candidates_pt_after_massFilter"]->Fill((*h_s_MassFilter)[s].pt());
			histos_th1f[b+"h_s_candidates_p_after_massFilter"]->Fill((*h_s_MassFilter)[s].p());
			histos_th1f[b+"h_s_candidates_energy_after_massFilter"]->Fill((*h_s_MassFilter)[s].energy());
			histos_th1f[b+"h_s_candidates_et_after_massFilter"]->Fill((*h_s_MassFilter)[s].et());
			histos_th1f[b+"h_s_candidates_mt_after_massFilter"]->Fill((*h_s_MassFilter)[s].mt());
			histos_th1f[b+"h_s_candidates_phi_after_massFilter"]->Fill((*h_s_MassFilter)[s].phi());
			histos_th1f[b+"h_s_candidates_theta_after_massFilter"]->Fill((*h_s_MassFilter)[s].theta());
			histos_th1f[b+"h_s_candidates_eta_after_massFilter"]->Fill((*h_s_MassFilter)[s].eta());
			histos_th1f[b+"h_s_candidates_vx_after_massFilter"]->Fill((*h_s_MassFilter)[s].vx());
			histos_th1f[b+"h_s_candidates_vy_after_massFilter"]->Fill((*h_s_MassFilter)[s].vy());
			histos_th1f[b+"h_s_candidates_vz_after_massFilter"]->Fill((*h_s_MassFilter)[s].vz());
			histos_th1f[b+"h_s_candidates_vertexNormalizedChi2_after_massFilter"]->Fill((*h_s_MassFilter)[s].vertexNormalizedChi2());
			
			
		}
		if((*h_s_MassFilter)[s].charge()==-1){
			histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter"]->Fill((*h_s_MassFilter)[s].mass());
			histos_th1f[b+"h_anti_s_candidates_pt_after_massFilter"]->Fill((*h_s_MassFilter)[s].pt());
			histos_th1f[b+"h_anti_s_candidates_p_after_massFilter"]->Fill((*h_s_MassFilter)[s].p());
			histos_th1f[b+"h_anti_s_candidates_energy_after_massFilter"]->Fill((*h_s_MassFilter)[s].energy());
			histos_th1f[b+"h_anti_s_candidates_et_after_massFilter"]->Fill((*h_s_MassFilter)[s].et());
			histos_th1f[b+"h_anti_s_candidates_mt_after_massFilter"]->Fill((*h_s_MassFilter)[s].mt());
			histos_th1f[b+"h_anti_s_candidates_phi_after_massFilter"]->Fill((*h_s_MassFilter)[s].phi());
			histos_th1f[b+"h_anti_s_candidates_theta_after_massFilter"]->Fill((*h_s_MassFilter)[s].theta());
			histos_th1f[b+"h_anti_s_candidates_eta_after_massFilter"]->Fill((*h_s_MassFilter)[s].eta());
			histos_th1f[b+"h_anti_s_candidates_vx_after_massFilter"]->Fill((*h_s_MassFilter)[s].vx());
			histos_th1f[b+"h_anti_s_candidates_vy_after_massFilter"]->Fill((*h_s_MassFilter)[s].vy());
			histos_th1f[b+"h_anti_s_candidates_vz_after_massFilter"]->Fill((*h_s_MassFilter)[s].vz());
			histos_th1f[b+"h_anti_s_candidates_vertexNormalizedChi2_after_massFilter"]->Fill((*h_s_MassFilter)[s].vertexNormalizedChi2());
		
		}
		
		//THE BELOW DOES NOT WORK: I CANNOT ACCESS THE DAUGHTERS OF THESE S CANDIDATES, I GET A product not found exception. ALTHOUGH THE DAUGHTERS ARE PUT INSIDE IN THE massfilter...
		//some plots for the S particle daughters: the first added daughter is the Lamda, the second is the Kshort (look at LambdaKshortVertexFilter.cc for this)
/*		if((*h_s_MassFilter)[s].charge()==1){
			const reco::Candidate *daughter0 = (*h_s_MassFilter)[s].daughter(0); 
			histos_th1f[b+"h_s_candidates_mass_after_massFilter_Lambda_mass"]->Fill((*h_s_MassFilter)[s].daughter(0)->mass());
			histos_th1f[b+"h_s_candidates_mass_after_massFilter_Kshort_mass"]->Fill((*h_s_MassFilter)[s].daughter(1)->mass());
		}
		if((*h_s_MassFilter)[s].charge()==-1){
			histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter_Lambda_mass"]->Fill((*h_s_MassFilter)[s].daughter(0)->mass());
			histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter_Kshort_mass"]->Fill((*h_s_MassFilter)[s].daughter(1)->mass());
		}
*/
	}
  }


} //end of analyzer


void Analyzer_V0_angular_correlation::endJob()
{
}

Analyzer_V0_angular_correlation::~Analyzer_V0_angular_correlation()
{
}


DEFINE_FWK_MODULE(Analyzer_V0_angular_correlation);
