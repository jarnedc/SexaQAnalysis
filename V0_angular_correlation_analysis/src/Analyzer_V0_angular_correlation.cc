#include "SexaQAnalysis/V0_angular_correlation_analysis/interface/Analyzer_V0_angular_correlation.h"
#include <typeinfo>


Analyzer_V0_angular_correlation::Analyzer_V0_angular_correlation(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  //m_rCandsTag(pset.getParameter<edm::InputTag>("resonCandidates")),
  m_bsTag    (pset.getParameter<edm::InputTag>("beamspot")),
  m_genParticlesTag(pset.getParameter<edm::InputTag>("genCollection")),
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
  //m_TKMETTag(pset.getParameter<edm::InputTag>("TKMETCollection")),
  //m_rCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_rCandsTag)),
  //m_sCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_sCandsTag)),
  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_genParticlesToken(consumes<vector<reco::GenParticle> >(m_genParticlesTag)),
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
    TFileDirectory dir_InitialProducer = m_fs->mkdir("InitialProducer");
    histos_th1f[b+"h_h_nPVs"]= dir_InitialProducer.make<TH1F>(b+"h_h_nPVs",b+"h_h_nPVs; nPVs",1000,0,1000);
    //histos_th1f[b+"h_h_nPVs"]= dir_InitialProducer.make<TH1F>(b+"h_h_nPVs",b+"h_h_nPVs; nPVs",1000,0,1000);
    histos_th1f[b+"h_h_nelectrons"]= dir_InitialProducer.make<TH1F>(b+"h_h_nelectrons",b+"h_h_nelectrons; nelectrons",1000,0,1000);
    histos_th1f[b+"h_h_njets"]= dir_InitialProducer.make<TH1F>(b+"h_h_njets",b+"h_h_njets; njets",1000,0,1000);
    histos_th1f[b+"h_h_nkshorts"]= dir_InitialProducer.make<TH1F>(b+"h_h_nkshorts",b+"h_h_nkshorts; nkshorts",1000,0,1000);
    histos_th1f[b+"h_h_nlambdas"]= dir_InitialProducer.make<TH1F>(b+"h_h_nlambdas",b+"h_h_nlambdas; nlambdas",1000,0,1000);
    histos_th1f[b+"h_h_nmuons"]= dir_InitialProducer.make<TH1F>(b+"h_h_nmuons",b+"h_h_nmuons; nmuons",1000,0,1000);
    histos_th1f[b+"h_h_ntracks"]= dir_InitialProducer.make<TH1F>(b+"h_h_ntracks",b+"h_h_ntracks; ntracks",1000,0,1000);
  
    histos_th1f[b+"h_h_TKMET_pT"]= dir_InitialProducer.make<TH1F>(b+"h_h_TKMET_pT",b+"h_h_TKMET_pT; ",2000,-1000,1000);
    histos_th1f[b+"h_h_MET_pT"]= dir_InitialProducer.make<TH1F>(b+"h_h_MET_pT",b+"h_h_MET_pT; ",2000,-1000,1000);
    histos_th1f[b+"h_h_TwoTopJets_pT_0"]= dir_InitialProducer.make<TH1F>(b+"h_h_TwoTopJets_pT_0",b+"h_h_TwoTopJets_pT_0; ",2000,-1000,1000);
    histos_th1f[b+"h_h_TwoTopJets_pT_1"]= dir_InitialProducer.make<TH1F>(b+"h_h_TwoTopJets_pT_1",b+"h_h_TwoTopJets_pT_1; ",2000,-1000,1000);
    histos_th1f[b+"h_h_HT"]= dir_InitialProducer.make<TH1F>(b+"h_h_HT",b+"h_h_HT; ",5000,0,5000);
    histos_th1f[b+"h_h_TKHT"]= dir_InitialProducer.make<TH1F>(b+"h_h_TKHT",b+"h_h_TKHT; ",5000,0,5000);

    //for beamspot
    TFileDirectory dir_Beamspot = m_fs->mkdir("Beamspot");
    histos_th1f[b+"h_beamspot_vx"]= dir_Beamspot.make<TH1F>(b+"h_beamspot_vx",b+"h_beamspot_vx; vx beamspot (cm)",1000,-10,10);
    histos_th1f[b+"h_beamspot_vy"]= dir_Beamspot.make<TH1F>(b+"h_beamspot_vy",b+"h_beamspot_vy; vy beamspot (cm)",1000,-10,10);
    histos_th1f[b+"h_beamspot_vx_std_dev"]= dir_Beamspot.make<TH1F>(b+"h_beamspot_vx_std_dev",b+"h_beamspot_vx_std_dev; vx beamspot std dev (cm)",1000,0,10);
    histos_th1f[b+"h_beamspot_vy_std_dev"]= dir_Beamspot.make<TH1F>(b+"h_beamspot_vy_std_dev",b+"h_beamspot_vy_std_dev; vy beamspot (cm)",1000,0,10);
    histos_th1f[b+"h_beamspot_lxy"]= dir_Beamspot.make<TH1F>(b+"h_beamspot_lxy",b+"h_beamspot_lxy; lxy beamspot ((0,0,0) to beamspot) (cm)",1000,0,10);
    histos_th1f[b+"h_beamspot_lxy_std_dev"]= dir_Beamspot.make<TH1F>(b+"h_beamspot_lxy_std_dev",b+"h_beamspot_lxy_std_dev; lxy beamspot ((0,0,0) to beamspot) std dev (cm)",1000,0,10);

    //for the gen particles
    TFileDirectory dir_genParticles = m_fs->mkdir("genParticles");
    histos_th1f[b+"h_genParticles_pdgid"]= dir_genParticles.make<TH1F>(b+"h_genParticles_pdgid", b+"h_genParticles_pdgid; pdgId of Gen particles; pdgId ",3000,-15000,15000);
    //kinematics of the generated Xi
    TFileDirectory dir_genXi1820 = dir_genParticles.mkdir("Xi1820");
	histos_th1f[b+"h_genParticles_Xi1820_pt_all"]= dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_pt_all", b+"h_genParticles_Xi1820_pt_all; #Xi(1820) pt ; pT (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Xi1820_p_all"] = dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_p_all", b+"h_genParticles_Xi1820_p_all; #Xi(1820) p ; p (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Xi1820_energy_all"] = dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_energy_all", b+"h_genParticles_Xi1820_energy_all; #Xi(1820) energy ; E (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Xi1820_et_all"] = dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_et_all", b+"h_genParticles_Xi1820_et_all; #Xi(1820) transversal E ; eT (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Xi1820_mt_all"] = dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_mt_all", b+"h_genParticles_Xi1820_mt_all; #Xi(1820) transversal mass ; mT (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Xi1820_phi_all"] = dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_phi_all", b+"h_genParticles_Xi1820_phi_all; #Xi(1820) phi ; phi (rad) ",8000,-4,4);
	histos_th1f[b+"h_genParticles_Xi1820_theta_all"]= dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_theta_all", b+"h_genParticles_Xi1820_theta_all; #Xi(1820) theta; theta (rad) ",800,-4,4);
	histos_th1f[b+"h_genParticles_Xi1820_eta_all"]= dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_eta_all", b+"h_genParticles_Xi1820_eta_all; #Xi(1820) eta; eta ",800,-10,10);
	histos_th1f[b+"h_genParticles_Xi1820_vx_all"] = dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_vx_all", b+"h_genParticles_Xi1820_vx_all; #Xi(1820) vx; vx (cm) ",2000,-100,100);
	histos_th1f[b+"h_genParticles_Xi1820_vy_all"]= dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_vy_all", b+"h_genParticles_Xi1820_vy_all; #Xi(1820) vy; vy (cm) ",2000,-100,100);
	histos_th1f[b+"h_genParticles_Xi1820_vz_all"] = dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_vz_all", b+"h_genParticles_Xi1820_vz_all; #Xi(1820) vz; vz (cm) ",2000,-100,100);
	histos_th1f[b+"h_genParticles_Xi1820_vertexNormalizedChi2_all"] = dir_genXi1820.make<TH1F>(b+"h_genParticles_Xi1820_vertexNormalizedChi2_all", b+"h_genParticles_Xi1820_vertexNormalizedChi2_all; #Xi(1820) vertex normalised Chi2  ",1000, 0,100);
    //kinematics of the generated Lambda
    TFileDirectory dir_genLambda = dir_genParticles.mkdir("Lambda");
	histos_th1f[b+"h_genParticles_Lambda_pt_all"]= dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_pt_all", b+"h_genParticles_Lambda_pt_all; Lambda pt ; pT (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Lambda_p_all"] = dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_p_all", b+"h_genParticles_Lambda_p_all; Lambda p ; p (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Lambda_energy_all"] = dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_energy_all", b+"h_genParticles_Lambda_energy_all; Lambda energy ; E (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Lambda_et_all"] = dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_et_all", b+"h_genParticles_Lambda_et_all; Lambda transversal E ; eT (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Lambda_mt_all"] = dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_mt_all", b+"h_genParticles_Lambda_mt_all; Lambda transversal mass ; mT (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Lambda_phi_all"] = dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_phi_all", b+"h_genParticles_Lambda_phi_all; Lambda phi ; phi (rad) ",8000,-4,4);
	histos_th1f[b+"h_genParticles_Lambda_theta_all"]= dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_theta_all", b+"h_genParticles_Lambda_theta_all; Lambda theta ; theta (rad) ",800,-4,4);
	histos_th1f[b+"h_genParticles_Lambda_eta_all"]= dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_eta_all", b+"h_genParticles_Lambda_eta_all; Lambda eta ; eta ",800,-10,10);
	histos_th1f[b+"h_genParticles_Lambda_vx_all"] = dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_vx_all", b+"h_genParticles_Lambda_vx_all; Lambda vx ; vx (cm) ",2000,-100,100);
	histos_th1f[b+"h_genParticles_Lambda_vy_all"]= dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_vy_all", b+"h_genParticles_Lambda_vy_all; Lambda vy ; vy (cm) ",2000,-100,100);
	histos_th1f[b+"h_genParticles_Lambda_vz_all"] = dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_vz_all", b+"h_genParticles_Lambda_vz_all; Lambda vz ; vz (cm) ",2000,-100,100);
	histos_th1f[b+"h_genParticles_Lambda_vertexNormalizedChi2_all"] = dir_genLambda.make<TH1F>(b+"h_genParticles_Lambda_vertexNormalizedChi2_all", b+"h_genParticles_Lambda_vertexNormalizedChi2_all; Lambda vertex normalised Chi2  ",1000, 0,100);

    //kinematics of the generated Kshort
    TFileDirectory dir_genKshort = dir_genParticles.mkdir("Kshort");
	histos_th1f[b+"h_genParticles_Kshort_pt_all"]= dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_pt_all", b+"h_genParticles_Kshort_pt_all; Kshort pt ; pT (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Kshort_p_all"] = dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_p_all", b+"h_genParticles_Kshort_p_all; Kshort p ; p (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Kshort_energy_all"] = dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_energy_all", b+"h_genParticles_Kshort_energy_all; Kshort energy ; E (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Kshort_et_all"] = dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_et_all", b+"h_genParticles_Kshort_et_all; Kshort transversal E ; eT (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Kshort_mt_all"] = dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_mt_all", b+"h_genParticles_Kshort_mt_all; Kshort transversal mass ; mT (GeV) ",4000,0,20);
	histos_th1f[b+"h_genParticles_Kshort_phi_all"] = dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_phi_all", b+"h_genParticles_Kshort_phi_all; Kshort phi ; phi (rad) ",8000,-4,4);
	histos_th1f[b+"h_genParticles_Kshort_theta_all"]= dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_theta_all", b+"h_genParticles_Kshort_theta_all; Kshort theta ; theta (rad) ",800,-4,4);
	histos_th1f[b+"h_genParticles_Kshort_eta_all"]= dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_eta_all", b+"h_genParticles_Kshort_eta_all; Kshort eta ; eta ",800,-10,10);
	histos_th1f[b+"h_genParticles_Kshort_vx_all"] = dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_vx_all", b+"h_genParticles_Kshort_vx_all; Kshort vx ; vx (cm) ",2000,-100,100);
	histos_th1f[b+"h_genParticles_Kshort_vy_all"]= dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_vy_all", b+"h_genParticles_Kshort_vy_all; Kshort vy ; vy (cm) ",2000,-100,100);
	histos_th1f[b+"h_genParticles_Kshort_vz_all"] = dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_vz_all", b+"h_genParticles_Kshort_vz_all; Kshort vz ; vz (cm) ",2000,-100,100);
	histos_th1f[b+"h_genParticles_Kshort_vertexNormalizedChi2_all"] = dir_genKshort.make<TH1F>(b+"h_genParticles_Kshort_vertexNormalizedChi2_all", b+"h_genParticles_Kshort_vertexNormalizedChi2_all; Kshort vertex normalised Chi2  ",1000, 0,100);




   //masses from the original kshort and lambda collection in the original V0 collections
    TFileDirectory dir_OriginalV0s = m_fs->mkdir("OriginalV0s");
    TFileDirectory dir_OriginalV0s_masses = dir_OriginalV0s.mkdir("masses");
   histos_th1f[b+"h_Lambda_mass"] = dir_OriginalV0s_masses.make<TH1F>(b+"h_Lambda_mass",b+"h_Lambda_mass; Lambda mass from original V0 collection",3000,0,3);   
   histos_th1f[b+"h_Kshort_mass"] = dir_OriginalV0s_masses.make<TH1F>(b+"h_Kshort_mass",b+"h_Kshort_mass; Kshort mass from original V0 collection",3000,0,3);   
   //comparison between the gen particles and the original V0s
    TFileDirectory dir_comparison_gen_OriginalV0s = dir_OriginalV0s.mkdir("Comparison with genererated");

   //Lambda
    TFileDirectory dir_comparison_gen_OriginalV0s_Lambda = dir_comparison_gen_OriginalV0s.mkdir("Lambda");
   histos_th1f[b+"h_genPartices_Lambda_origV0_pdgId"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_pdgId",b+"h_genPartices_Lambda_origV0_pdgId; pdgId;",30000,-15000,15000);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_pt"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_pt",b+"h_genPartices_Lambda_origV0_diff_pt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_mt"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_mt",b+"h_genPartices_Lambda_origV0_diff_mt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_phi"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_phi",b+"h_genPartices_Lambda_origV0_diff_phi; diff(phi) (gen-reco)[rad];",800,-4,4);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_eta"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_eta",b+"h_genPartices_Lambda_origV0_diff_eta; diff(eta) (gen-reco);",200,-10,10);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vx"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_vx",b+"h_genPartices_Lambda_origV0_diff_vx; diff(vx) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vx_error_vx_smaller_0.01"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_vx_error_vx_smaller_0.01",b+"h_genPartices_Lambda_origV0_diff_vx_error_vx_smaller_0.01; diff(vx) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vy"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_vy",b+"h_genPartices_Lambda_origV0_diff_vy; diff(vy) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vy_error_vy_smaller_0.01"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_vy_error_vy_smaller_0.01",b+"h_genPartices_Lambda_origV0_diff_vy_error_vy_smaller_0.01; diff(vy) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vz"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_vz",b+"h_genPartices_Lambda_origV0_diff_vz; diff(vz) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vxy"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_vxy",b+"h_genPartices_Lambda_origV0_diff_vxy; diff(vy) |gen-reco|[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vxy_error_vx_and_vy_smaller_0.01"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH1F>(b+"h_genPartices_Lambda_origV0_diff_vxy_error_vx_and_vy_smaller_0.01",b+"h_genPartices_Lambda_origV0_diff_vxy_error_vx_and_vy_smaller_0.01; diff(vy) |gen-reco|[cm];",10000,-50,50);
   histos_th2f[b+"h_genPartices_Lambda_origV0_pt_diff_vx"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH2F>(b+"h_genPartices_Lambda_origV0_pt_diff_vx",b+"h_genPartices_Lambda_origV0_pt_diff_vx; gen particle pt [GeV]; abs(diff(vx)) |gen-reco|[cm];",5000, 0,50, 5000, 0, 50);
   histos_th2f[b+"h_genPartices_Lambda_origV0_diff_vx_error_vx"] = dir_comparison_gen_OriginalV0s_Lambda.make<TH2F>(b+"h_genPartices_Lambda_origV0_diff_vx_error_vx",b+"h_genPartices_Lambda_origV0_diff_vx_error_vx; abs(diff(vx)) |gen-reco|[cm]; error vx reconstructed [cm];",5000, 0,50, 5000, 0, 50);

   //Kshort
    TFileDirectory dir_comparison_gen_OriginalV0s_Kshort = dir_comparison_gen_OriginalV0s.mkdir("Kshort");
   histos_th1f[b+"h_genPartices_Kshort_origV0_pdgId"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_pdgId",b+"h_genPartices_Kshort_origV0_pdgId; pdgId;",30000,-15000,15000);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_pt"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_pt",b+"h_genPartices_Kshort_origV0_diff_pt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_mt"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_mt",b+"h_genPartices_Kshort_origV0_diff_mt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_phi"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_phi",b+"h_genPartices_Kshort_origV0_diff_phi; diff(phi) (gen-reco)[rad];",800,-4,4);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_eta"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_eta",b+"h_genPartices_Kshort_origV0_diff_eta; diff(eta) (gen-reco);",200,-10,10);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vx"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_vx",b+"h_genPartices_Kshort_origV0_diff_vx; diff(vx) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vx_error_vx_smaller_0.01"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_vx_error_vx_smaller_0.01","h_genPartices_Kshort_origV0_diff_vx_error_vx_smaller_0.01; diff(vx) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vy"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_vy",b+"h_genPartices_Kshort_origV0_diff_vy; diff(vy) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vy_error_vy_smaller_0.01"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_vy_error_vy_smaller_0.01","h_genPartices_Kshort_origV0_diff_vy_error_vy_smaller_0.01; diff(vy) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vz"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_vz",b+"h_genPartices_Kshort_origV0_diff_vz; diff(vz) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vxy"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_vxy",b+"h_genPartices_Kshort_origV0_diff_vxy; diff(vxy) |gen-reco|[cm];",10000,-50,50);

   histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vxy_error_vx_and_vy_smaller_0.01"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH1F>(b+"h_genPartices_Kshort_origV0_diff_vxy_error_vx_and_vy_smaller_0.01",b+"h_genPartices_Kshort_origV0_diff_vxy_error_vx_and_vy_smaller_0.01; diff(vxy) |gen-reco|[cm];",10000,-50,50);
   histos_th2f[b+"h_genPartices_Kshort_origV0_pt_diff_vx"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH2F>(b+"h_genPartices_Kshort_origV0_pt_diff_vx",b+"h_genPartices_Kshort_origV0_pt_diff_vx; gen particle pt [GeV]; abs(diff(vx)) |gen-reco|[cm];",5000, 0,50, 5000, 0, 50);
   histos_th2f[b+"h_genPartices_Kshort_origV0_diff_vx_error_vx"] = dir_comparison_gen_OriginalV0s_Kshort.make<TH2F>(b+"h_genPartices_Kshort_origV0_diff_vx_error_vx",b+"h_genPartices_Kshort_origV0_diff_vx_error_vx; abs(diff(vx)) |gen-reco|[cm]; error vx reconstructed [cm];",5000, 0,50, 5000, 0, 50);

			
   //for angular correlation between Ks and L0 in the original V0 collections 
    TFileDirectory dir_Angular_correlation = dir_OriginalV0s.mkdir("Angular_correlation");
   histos_th1f[b+"h_L0_Ks_delta_phi"]= dir_Angular_correlation.make<TH1F>(b+"h_L0_Ks_delta_phi",b+"h_L0_Ks_delta_phi; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]",1000,-4,4); 
   histos_th1f[b+"h_L0_Ks_delta_eta"]= dir_Angular_correlation.make<TH1F>(b+"h_L0_Ks_delta_eta",b+"h_L0_Ks_delta_eta; delta eta",4000,-10,10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta"]= dir_Angular_correlation.make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta",b+"h_L0_Ks_delta_phi_delta_eta; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta_eta",100,-4,4, 100, -10, 10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_no_cent"]= dir_Angular_correlation.make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta_no_cent",b+"h_L0_Ks_delta_phi_delta_eta_no_cent; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta_eta",1000,-4,4, 4000, -10, 10); 

    //vertex and error on vertex for the Ks and the Lambda
    TFileDirectory dir_displacement = dir_OriginalV0s.mkdir("displacement");
    TFileDirectory dir_displacement_Lambda = dir_displacement.mkdir("Lambda");    
    histos_th1f[b+"h_Lambda_lxy"] = dir_displacement_Lambda.make<TH1F>(b+"h_Lambda_lxy",b+"h_Lambda_lxy; lxy Lambda from original V0 collection (cm)",4000,-20,20);   
    histos_th1f[b+"h_Lambda_lxy_error"] = dir_displacement_Lambda.make<TH1F>(b+"h_Lambda_lxy_error",b+"h_Lambda_lxy_error; error lxy Lambda from original V0 collection (cm)",4000,-20,20);   
    histos_th1f[b+"h_Lambda_lxy_signed"] = dir_displacement_Lambda.make<TH1F>(b+"h_Lambda_lxy_signed",b+"h_Lambda_lxy_signed; signed lxy Lambda from original V0 collection (cm)",4000,-20,20);   
    

    TFileDirectory dir_displacement_Kshort = dir_displacement.mkdir("Kshort");
    histos_th1f[b+"h_Kshort_lxy"] = dir_displacement_Kshort.make<TH1F>(b+"h_Kshort_lxy",b+"h_Kshort_lxy; lxy Kshort from original V0 collection (cm)",4000,-20,20);   
    histos_th1f[b+"h_Kshort_lxy_error"] = dir_displacement_Kshort.make<TH1F>(b+"h_Kshort_lxy_error",b+"h_Kshort_lxy_error; error lxy Kshort from original V0 collection (cm)",4000,-20,20);   
    histos_th1f[b+"h_Kshort_lxy_signed"] = dir_displacement_Kshort.make<TH1F>(b+"h_Kshort_lxy_signed",b+"h_Kshort_lxy_signed; signed lxy Kshort from original V0 collection (cm)",4000,-20,20);   
    

   //for angular correlation between Ks and L0 from different events, from the original V0 collections
    TFileDirectory dir_Angular_correlation_cross_events = dir_OriginalV0s.mkdir("Angular_correlation_cross_events");
   histos_th1f[b+"h_L0_Ks_delta_phi_prev_and_current"]= dir_Angular_correlation_cross_events.make<TH1F>(b+"h_L0_Ks_delta_phi_prev_and_current",b+"h_L0_Ks_delta_phi_prev_and_current; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad] prev event and current event;",1000,-4,4); 
   histos_th1f[b+"h_L0_Ks_delta_eta_prev_and_current"]= dir_Angular_correlation_cross_events.make<TH1F>(b+"h_L0_Ks_delta_eta_prev_and_current",b+"h_L0_Ks_delta_eta_prev_and_current; delta eta prev event and current event;",4000,-10,10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]= dir_Angular_correlation_cross_events.make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current",b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad] prev event and current event; delta_eta previous event and current event",1000,-4,4, 4000, -10, 10); 
   
    //masses using the LambdaKshortVertexFilter results 
    TFileDirectory dir_LambdaKshortVertexFilter = m_fs->mkdir("LambdaKshortVertexFilter");
    TFileDirectory dir_masses_S = dir_LambdaKshortVertexFilter.mkdir("masses S");
    
    histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter; S mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm; S mass (GeV)",20000, 0,20);
    histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_error_smaller_0.01cm"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_error_smaller_0.01cm", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_error_smaller_0.01cm; S mass (GeV)",20000, 0,20);
    histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; S mass (GeV)",20000, 0,20);
    histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5; S mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; S mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_s_candidates_vertex_beamspot_dist_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_vertex_beamspot_dist_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm", b+"h_s_candidates_vertex_beamspot_dist_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; lxy (cm)",20000, 0,20);
    histos_th1f[b+"h_s_candidates_error_vertex_beamspot_dist_error_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_error_vertex_beamspot_dist_error_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm", b+"h_s_candidates_error_vertex_beamspot_dist_error_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; error lxy (cm)",20000, 0,20);
    histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8; S mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"] = dir_masses_S.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; S mass [GeV]",20000, 0,20);

    TFileDirectory dir_masses_anti_S = dir_LambdaKshortVertexFilter.mkdir("masses anti-S (BLINDED!!!)");
    
    histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter"] = dir_masses_anti_S.make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter; S mass after LambdaKshortVertexFilter",20000, 0,20);
    histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"] = dir_masses_anti_S.make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; anti S mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5"] = dir_masses_anti_S.make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5; anti S mass [GeV]",20000, 0,20);

    histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"] = dir_masses_anti_S.make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; anti S mass [GeV]",20000, 0,20);

    histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8"] = dir_masses_anti_S.make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8; anti S mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"] = dir_masses_anti_S.make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; anti S mass [GeV]",20000, 0,20);

    //masses of all the reconances
    TFileDirectory dir_masses_all_R = dir_LambdaKshortVertexFilter.mkdir("masses_all_r");
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8; R mass [GeV]",20000, 0,20);
    //for the r candidate, you made some cuts. Now investigate these cuts:    
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_for_delta_R_slammer_0.8"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_for_delta_R_slammer_0.8", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_for_delta_R_slammer_0.8; lxy(S interaction vertex, beamspot)", 5000, 0, 50);
    histos_th1f[b+"h_all_r_candidates_delta_R"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_delta_R", b+"h_all_r_candidates_delta_R; #Delta R(K_{S}^{0},#Lambda^{0})", 5000, 0, 50);
    histos_th1f[b+"h_all_r_candidates_r_vertex_beamspot_distance"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_r_vertex_beamspot_distance", b+"h_all_r_candidates_r_vertex_beamspot_distance; lxy(S interaction vertex, beamspot)", 5000, 0, 50);
    histos_th1f[b+"h_all_r_candidates_r_vertex_beamspot_distance_error"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_r_vertex_beamspot_distance_error", b+"h_all_r_candidates_r_vertex_beamspot_distance_error; error lxy(S interaction vertex, beamspot)", 5000, 0, 50);
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_error_r_vertex_beamspot_distance_for_delta_R_slammer_0.8"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_error_r_vertex_beamspot_distance_for_delta_R_slammer_0.8", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_error_r_vertex_beamspot_distance_for_delta_R_slammer_0.8;error lxy(S interaction vertex, beamspot)", 5000, 0, 50);
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_error_r_vertex_beamspot_distance_for_r_vertex_beamspot_distance_smaller_0.1"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_error_r_vertex_beamspot_distance_for_r_vertex_beamspot_distance_smaller_0.1", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_error_r_vertex_beamspot_distance_for_r_vertex_beamspot_distance_smaller_0.1;error lxy(S interaction vertex, beamspot)", 5000, 0, 50);
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_for_r_vertex_beamspot_distance_smaller_0.1"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_for_r_vertex_beamspot_distance_smaller_0.1", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_for_r_vertex_beamspot_distance_smaller_0.1;#Delta R(K_{S}^{0},#Lambda^{0})", 2000, 0, 20);
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_for_error_r_vertex_beamspot_distance_smaller_0.01"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_for_error_r_vertex_beamspot_distance_smaller_0.01", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_for_error_r_vertex_beamspot_distance_smaller_0.01;lxy(S interaction vertex, beamspot)", 5000, 0, 50);
    histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_for_error_r_vertex_beamspot_distance_smaller_0.01"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_for_error_r_vertex_beamspot_distance_smaller_0.01", b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_for_error_r_vertex_beamspot_distance_smaller_0.01;#Delta R(K_{S}^{0},#Lambda^{0})", 2000, 0, 20);
    
   histos_th1f[b+"h_all_r_candidates_mass_PCA_xy_smaller_0p1"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_PCA_xy_smaller_0p1", b+"h_all_r_candidates_mass_PCA_xy_smaller_0p1; R Mass",20000, -20,20);
   histos_th1f[b+"h_all_r_candidates_mass_PCA_z_smaller_0p05"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_PCA_z_smaller_0p05", b+"h_all_r_candidates_mass_PCA_z_smaller_0p05; R Mass",20000, -20,20);
   histos_th1f[b+"h_all_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8", b+"h_all_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8; R Mass",20000, -20,20);
   histos_th1f[b+"h_all_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8", b+"h_all_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8; R Mass",20000, -20,20);
   histos_th1f[b+"h_all_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8_PCA_xy_smaller_0p5"] = dir_masses_all_R.make<TH1F>(b+"h_all_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8_PCA_xy_smaller_0p5", b+"h_all_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8_PCA_xy_smaller_0p5; R Mass",20000, -20,20);
    //masses of the reconances decaying to  protons
    TFileDirectory dir_masses_R = dir_LambdaKshortVertexFilter.mkdir("masses_r");
    histos_th1f[b+"h_r_candidates_mass_after_LambdaKshortVertexFilter"] = dir_masses_R.make<TH1F>(b+"h_r_candidates_mass_after_LambdaKshortVertexFilter", b+"h_r_candidates_mass_after_LambdaKshortVertexFilter; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm"] = dir_masses_R.make<TH1F>(b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm", b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8"] = dir_masses_R.make<TH1F>(b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8", b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8"] = dir_masses_R.make<TH1F>(b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8", b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8"] = dir_masses_R.make<TH1F>(b+"h_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8", b+"h_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8"] = dir_masses_R.make<TH1F>(b+"h_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8", b+"h_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8; R mass [GeV]",20000, 0,20);
    
    //masses of the reconances decaying to anti-proton 
    TFileDirectory dir_masses_anti_R = dir_LambdaKshortVertexFilter.mkdir("masses_anti_r");
    histos_th1f[b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter"] = dir_masses_anti_R.make<TH1F>(b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter", b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm"] = dir_masses_anti_R.make<TH1F>(b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm", b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8"] = dir_masses_anti_R.make<TH1F>(b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8", b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8"] = dir_masses_anti_R.make<TH1F>(b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8", b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_anti_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8"] = dir_masses_R.make<TH1F>(b+"h_anti_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8", b+"h_anti_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8; R mass [GeV]",20000, 0,20);
    histos_th1f[b+"h_anti_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8"] = dir_masses_R.make<TH1F>(b+"h_anti_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8", b+"h_anti_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8; R mass [GeV]",20000, 0,20);
    

   //PCA 
    TFileDirectory dir_PCA_R = dir_LambdaKshortVertexFilter.mkdir("PCA_R");
   histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter", b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter; dxy [cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dxy_signed_after_LambdaKshortVertexFilter"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dxy_signed_after_LambdaKshortVertexFilter", b+"h_r_candidates_dxy_signed_after_LambdaKshortVertexFilter; dxy signed[cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dz_after_LambdaKshortVertexFilter", b+"h_r_candidates_dz_after_LambdaKshortVertexFilter; dz [cm]",20000, -20,20);
   histos_th2f[b+"h_r_candidates_dxy_dz_after_LambdaKshortVertexFilter"] = dir_PCA_R.make<TH2F>(b+"h_r_candidates_dxy_dz_after_LambdaKshortVertexFilter", b+"h_r_candidates_dxy_dz_after_LambdaKshortVertexFilter;dxy [cm]; dz [cm]",4000, -20,20, 4000, -20,20);
   histos_th2f[b+"h_r_candidates_signed_dxy_dz_after_LambdaKshortVertexFilter"] = dir_PCA_R.make<TH2F>(b+"h_r_candidates_signed_dxy_dz_after_LambdaKshortVertexFilter", b+"h_r_candidates_signed_dxy_dz_after_LambdaKshortVertexFilter; signed dxy [cm]; dz [cm]",4000, -20,20, 4000, -20,20);
   
   histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_smaller_11"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_smaller_11", b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_smaller_11; dxy [cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_smaller_11"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_smaller_11", b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_smaller_11; dz [cm]",20000, -20,20);
   
   //pile-up dependence of the dxy and dz 
   histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_10_smaller_or_equal_to_20"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_10_smaller_or_equal_to_20", b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_10_smaller_or_equal_to_20; dxy [cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_10_smaller_or_equal_to_20"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_10_smaller_or_equal_to_20", b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_10_smaller_or_equal_to_20; dz [cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_20_smaller_or_equal_to_30"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_20_smaller_or_equal_to_30", b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_20_smaller_or_equal_to_30; dxy [cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_20_smaller_or_equal_to_30"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_20_smaller_or_equal_to_30", b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_20_smaller_or_equal_to_30; dz [cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_30_smaller_or_equal_to_40"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_30_smaller_or_equal_to_40", b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_30_smaller_or_equal_to_40; dxy [cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_30_smaller_or_equal_to_40"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_30_smaller_or_equal_to_40", b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_30_smaller_or_equal_to_40; dz [cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_40"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_40", b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_40; dxy [cm]",20000, -20,20);
   histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_40"] = dir_PCA_R.make<TH1F>(b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_40", b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_40; dz [cm]",20000, -20,20);

   //masses of Ks and Lambda after the LambdaKshortVertexFilter
    TFileDirectory dir_masses_daugthers = dir_LambdaKshortVertexFilter.mkdir("masses daughters");
   histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass_all"] = dir_masses_daugthers.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass_all", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass_all; Lambda mass after LambdaKshortVertexFilter",4000, 0,20);
   histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass_all"] = dir_masses_daugthers.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass_all", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass_all; Kshort mass after LambdaKshortVertexFilter",4000, 0,20);
   histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass"] = dir_masses_daugthers.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass; Lambda mass after LambdaKshortVertexFilter for S candidate",4000, 0,20);
   histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass"] = dir_masses_daugthers.make<TH1F>(b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass", b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass; Kshort mass after LambdaKshortVertexFilter for S candidate",4000, 0,20);
   histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass"] = dir_masses_daugthers.make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Lambda_mass; Lambda mass after LambdaKshortVertexFilter for anti S candidate",4000, 0,20);
   histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass"] = dir_masses_daugthers.make<TH1F>(b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass", b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_Kshort_mass; Kshort mass after LambdaKshortVertexFilter for anti S candidate",4000, 0,20);
  
  //other kinematic quantities after the LambdaKshortVertexFilter 
	//Lambda
    	TFileDirectory dir_Lambda_Kinematics = dir_LambdaKshortVertexFilter.mkdir("Lambda_Kinematics");
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_pt_all"]= dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_pt_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_pt_all; Lambda pt after the LambdaKshortVertexFilter; pT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_p_all"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_p_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_p_all; Lambda p after the LambdaKshortVertexFilter; p (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_energy_all"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_energy_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_energy_all; Lambda energy after the LambdaKshortVertexFilter; E (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_et_all"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_et_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_et_all; Lambda transversal E after the LambdaKshortVertexFilter; eT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_mt_all"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_mt_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_mt_all; Lambda transversal mass after the LambdaKshortVertexFilter; mT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_phi_all"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_phi_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_phi_all; Lambda phi after the LambdaKshortVertexFilter; phi (rad) ",8000,-4,4);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_theta_all"]= dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_theta_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_theta_all; Lambda theta after the LambdaKshortVertexFilter; theta (rad) ",800,-4,4);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_eta_all"]= dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_eta_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_eta_all; Lambda eta after the LambdaKshortVertexFilter; eta ",800,-10,10);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_all"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_all; Lambda vx after the LambdaKshortVertexFilter; vx (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_all"]= dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_all; Lambda vy after the LambdaKshortVertexFilter; vy (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vz_all"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vz_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vz_all; Lambda vz after the LambdaKshortVertexFilter; vz (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vertexNormalizedChi2_all"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vertexNormalizedChi2_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vertexNormalizedChi2_all; Lambda vertex normalised Chi2 after the LambdaKshortVertexFilter ",1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_std_dev"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_std_dev", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_std_dev;Lambda x vertex std dev (cm)",1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_std_dev"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_std_dev", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_std_dev;Lambda y vertex std dev (cm)",1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_lxy"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_lxy", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_lxy;Lambda lxy(beamspot and Lambda vertex) (cm)",1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_lxy_std_dev"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_lxy_std_dev", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_lxy_std_dev;Lambda lxy(beamspot and Lambda vertex) std_dev (cm)",1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_dxy"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_dxy", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_dxy; Lambda dxy (cm)",10000, -100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_dz"] = dir_Lambda_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_dz", b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_dz; Lambda dz (cm)",10000, -100,100);


	//Kshort
    	TFileDirectory dir_Kshort_Kinematics = dir_LambdaKshortVertexFilter.mkdir("Kshort_Kinematics");
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_pt_all"]= dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_pt_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_pt_all; Kshort pt after the LambdaKshortVertexFilter; pT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_p_all"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_p_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_p_all; Kshort p after the LambdaKshortVertexFilter; p (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_energy_all"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_energy_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_energy_all; Kshort energy after the LambdaKshortVertexFilter; E (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_et_all"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_et_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_et_all; Kshort transversal E after the LambdaKshortVertexFilter; eT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_mt_all"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_mt_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_mt_all; Kshort transversal mass after the LambdaKshortVertexFilter; mT (GeV) ",4000,0,20);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_phi_all"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_phi_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_phi_all; Kshort phi after the LambdaKshortVertexFilter; phi (rad) ",8000,-4,4);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_theta_all"]= dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_theta_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_theta_all; Kshort theta after the LambdaKshortVertexFilter; theta (rad) ",800,-4,4);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_eta_all"]= dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_eta_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_eta_all; Kshort eta after the LambdaKshortVertexFilter; eta ",800,-10,10);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_all"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_all; Kshort vx after the LambdaKshortVertexFilter; vx (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_all"]= dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_all; Kshort vy after the LambdaKshortVertexFilter; vy (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vz_all"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vz_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vz_all; Kshort vz after the LambdaKshortVertexFilter; vz (cm) ",2000,-100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vertexNormalizedChi2_all"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vertexNormalizedChi2_all", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vertexNormalizedChi2_all; Kshort vertex normalised Chi2 after the LambdaKshortVertexFilter ",1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_std_dev"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_std_dev", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_std_dev;Kshort x vertex std dev(cm)" ,1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_std_dev"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_std_dev", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_std_dev;Kshort y vertex std dev (cm)",1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_lxy"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_lxy", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_lxy;Kshort lxy(beamspot and Kshort vertex) (cm)",1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_lxy_std_dev"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_lxy_std_dev", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_lxy_std_dev;Kshort lxy(beamspot and Kshort vertex) std_dev (cm)",1000, 0,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_dxy"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_dxy", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_dxy; Kshort dxy",10000, -100,100);
	histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_dz"] = dir_Kshort_Kinematics.make<TH1F>(b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_dz", b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_dz; Kshort dz",10000, -100,100);

  
   //comparison between the reconstructed lambda and kshort and the gen gen particles
 
   TFileDirectory dir_comparison_generated_LambdaKshortVertexFilter = dir_LambdaKshortVertexFilter.mkdir("Comparison with generator particles");
   TFileDirectory dir_comparison_generated_LambdaKshortVertexFilter_Lambda = dir_comparison_generated_LambdaKshortVertexFilter.mkdir("Lambda");
      
   histos_th1f[b+"h_genPartices_Lambda_diff_pdgId"] = dir_comparison_generated_LambdaKshortVertexFilter_Lambda.make<TH1F>(b+"h_genPartices_Lambda_diff_pdgId",b+"h_genPartices_Lambda_diff_pdgId; pdgid as a check;",15000,-30000,30000);
   histos_th1f[b+"h_genPartices_Lambda_diff_pt"] = dir_comparison_generated_LambdaKshortVertexFilter_Lambda.make<TH1F>(b+"h_genPartices_Lambda_diff_pt",b+"h_genPartices_Lambda_diff_pt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Lambda_diff_mt"] = dir_comparison_generated_LambdaKshortVertexFilter_Lambda.make<TH1F>(b+"h_genPartices_Lambda_diff_mt",b+"h_genPartices_Lambda_diff_mt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Lambda_diff_phi"] = dir_comparison_generated_LambdaKshortVertexFilter_Lambda.make<TH1F>(b+"h_genPartices_Lambda_diff_phi",b+"h_genPartices_Lambda_diff_phi; diff(phi) (gen-reco)[rad];",800,-4,4);
   histos_th1f[b+"h_genPartices_Lambda_diff_eta"] = dir_comparison_generated_LambdaKshortVertexFilter_Lambda.make<TH1F>(b+"h_genPartices_Lambda_diff_eta",b+"h_genPartices_Lambda_diff_eta; diff(eta) (gen-reco);",200,-10,10);
   histos_th1f[b+"h_genPartices_Lambda_diff_vx"] = dir_comparison_generated_LambdaKshortVertexFilter_Lambda.make<TH1F>(b+"h_genPartices_Lambda_diff_vx",b+"h_genPartices_Lambda_diff_vx; diff(vx) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Lambda_diff_vy"] = dir_comparison_generated_LambdaKshortVertexFilter_Lambda.make<TH1F>(b+"h_genPartices_Lambda_diff_vy",b+"h_genPartices_Lambda_diff_vy; diff(vy) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Lambda_diff_vz"] = dir_comparison_generated_LambdaKshortVertexFilter_Lambda.make<TH1F>(b+"h_genPartices_Lambda_diff_vz",b+"h_genPartices_Lambda_diff_vz; diff(vz) (gen-reco)[cm];",10000,-50,50);

   TFileDirectory dir_comparison_generated_LambdaKshortVertexFilter_Kshort = dir_comparison_generated_LambdaKshortVertexFilter.mkdir("Kshort");
      
   histos_th1f[b+"h_genPartices_Kshort_diff_pdgId"] = dir_comparison_generated_LambdaKshortVertexFilter_Kshort.make<TH1F>(b+"h_genPartices_Kshort_diff_pdgId",b+"h_genPartices_Kshort_diff_pdgId; pdgid as a check;",15000,-30000,30000);
   histos_th1f[b+"h_genPartices_Kshort_diff_pt"] = dir_comparison_generated_LambdaKshortVertexFilter_Kshort.make<TH1F>(b+"h_genPartices_Kshort_diff_pt",b+"h_genPartices_Kshort_diff_pt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Kshort_diff_mt"] = dir_comparison_generated_LambdaKshortVertexFilter_Kshort.make<TH1F>(b+"h_genPartices_Kshort_diff_mt",b+"h_genPartices_Kshort_diff_mt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Kshort_diff_phi"] = dir_comparison_generated_LambdaKshortVertexFilter_Kshort.make<TH1F>(b+"h_genPartices_Kshort_diff_phi",b+"h_genPartices_Kshort_diff_phi; diff(phi) (gen-reco)[rad];",800,-4,4);
   histos_th1f[b+"h_genPartices_Kshort_diff_eta"] = dir_comparison_generated_LambdaKshortVertexFilter_Kshort.make<TH1F>(b+"h_genPartices_Kshort_diff_eta",b+"h_genPartices_Kshort_diff_eta; diff(eta) (gen-reco);",200,-10,10);
   histos_th1f[b+"h_genPartices_Kshort_diff_vx"] = dir_comparison_generated_LambdaKshortVertexFilter_Kshort.make<TH1F>(b+"h_genPartices_Kshort_diff_vx",b+"h_genPartices_Kshort_diff_vx; diff(vx) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Kshort_diff_vy"] = dir_comparison_generated_LambdaKshortVertexFilter_Kshort.make<TH1F>(b+"h_genPartices_Kshort_diff_vy",b+"h_genPartices_Kshort_diff_vy; diff(vy) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Kshort_diff_vz"] = dir_comparison_generated_LambdaKshortVertexFilter_Kshort.make<TH1F>(b+"h_genPartices_Kshort_diff_vz",b+"h_genPartices_Kshort_diff_vz; diff(vz) (gen-reco)[cm];",10000,-50,50);

   //for distance between beamspot and S candidate vertex 
   TFileDirectory dir_dist_beamspot_S = dir_LambdaKshortVertexFilter.mkdir("dist_beamspot_S");
   histos_th1f[b+"h_S_vtx_distance_to_beamspot"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot",b+"h_S_vtx_distance_to_beamspot; lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_S_vtx_distance_to_beamspot_error"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_error",b+"h_S_vtx_distance_to_beamspot_error; error lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   //try to see if the error on the vertex is dependent on the #PVs

   histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_smaller_11"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_for_nPVs_smaller_11",b+"h_S_vtx_distance_to_beamspot_for_nPVs_smaller_11; lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_10_smaller_or_equal_to_20"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_10_smaller_or_equal_to_20",b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_10_smaller_or_equal_to_20; lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_20_smaller_or_equal_to_30"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_20_smaller_or_equal_to_30",b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_20_smaller_or_equal_to_30; lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_30_smaller_or_equal_to_40"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_30_smaller_or_equal_to_40",b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_30_smaller_or_equal_to_40; lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_40"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_40",b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_40; lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);

   histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_smaller_11"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_smaller_11",b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_smaller_11; error lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_10_smaller_or_equal_to_20"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_10_smaller_or_equal_to_20",b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_10_smaller_or_equal_to_20; error lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_20_smaller_or_equal_to_30"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_20_smaller_or_equal_to_30",b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_20_smaller_or_equal_to_30; error lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_30_smaller_or_equal_to_40"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_30_smaller_or_equal_to_40",b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_30_smaller_or_equal_to_40; error lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_40"] = dir_dist_beamspot_S.make<TH1F>(b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_40",b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_40; error lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);

   histos_th2f[b+"h_delta_phi_S_vtx_distance_to_beamspot"]= dir_dist_beamspot_S.make<TH2F>(b+"h_delta_phi_S_vtx_distance_to_beamspot",b+"h_delta_phi_S_vtx_distance_to_beamspot; #Delta #Phi [rad]; lxy(S interaction vertex,  beamspot) [cm];",4000, -4, 4,3000,0,30); 
 //  histos_th3f[b+"h_delta_phi_delta_eta_S_vtx_distance_to_beamspot"]= dir_dist_beamspot_S.make<TH3F>(b+"h_delta_phi_delta_eta_S_vtx_distance_to_beamspot",b+"h_delta_phi_delta_eta_S_vtx_distance_to_beamspot; #Delta #Phi [rad]; #Delta #eta; lxy(S interaction vertex,  beamspot) [cm];",4000, -4, 4, 4000, -4, 4,3000,0,30); 
   histos_th2f[b+"h_delta_phi_S_vtx_distance_to_beamspot_vertex_beamspot_dist_error_smaller_0.01cm"]= dir_dist_beamspot_S.make<TH2F>(b+"h_delta_phi_S_vtx_distance_to_beamspot_vertex_beamspot_dist_error_smaller_0.01cm",b+"h_delta_phi_S_vtx_distance_to_beamspot_vertex_beamspot_dist_error_smaller_0.01cm; #Delta #Phi [rad]; lxy(S interaction vertex,  beamspot) [cm];",4000, -4, 4,3000,0,30); 
 //  histos_th3f[b+"h_delta_phi_delta_eta_S_vtx_distance_to_beamspot_vertex_beamspot_dist_error_smaller_0.01cm"]= dir_dist_beamspot_S.make<TH3F>(b+"h_delta_phi_delta_eta_S_vtx_distance_to_beamspot_vertex_beamspot_dist_error_smaller_0.01cm",b+"h_delta_phi_delta_eta_S_vtx_distance_to_beamspot_vertex_beamspot_dist_error_smaller_0.01cm; #Delta #Phi [rad]; #Delta #eta; lxy(S interaction vertex,  beamspot) [cm];",4000, -4, 4, 4000, -4, 4,3000,0,30); 
   histos_th2f[b+"h_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error"]= dir_dist_beamspot_S.make<TH2F>(b+"h_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error",b+"h_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error; lxy(S interaction vertex,  beamspot) [cm]; errror lxy(S interaction vertex,  beamspot) [cm];",3000, 0, 30,3000,0,30); 
   histos_th2f[b+"h_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error_for_delta_phi_between_1_and_2.5"]= dir_dist_beamspot_S.make<TH2F>(b+"h_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error_for_delta_phi_between_1_and_2.5",b+"h_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error_for_delta_phi_between_1_and_2.5; lxy(S interaction vertex,  beamspot) [cm]; errror lxy(S interaction vertex,  beamspot) [cm];",3000, 0, 30,3000,0,30); 
   histos_th2f[b+"h_delta_phi_S_vtx_distance_to_beamspot_significance"]= dir_dist_beamspot_S.make<TH2F>(b+"h_delta_phi_S_vtx_distance_to_beamspot_significance",b+"h_delta_phi_S_vtx_distance_to_beamspot_significance; #Delta #Phi [rad]; lxy(S interaction vertex,  beamspot) significance;",4000, -4, 4,3000,0,30); 
   histos_th2f[b+"h_S_vtx_distance_to_beamspot_vx_vy"]= dir_dist_beamspot_S.make<TH2F>(b+"h_S_vtx_distance_to_beamspot_vx_vy",b+"h_S_vtx_distance_to_beamspot_vx_vy; vx (distance beamspot to S annihilation vertex); vy (distance beamspot to S annihilation vertex)",1000, -50, 50,1000, -50, 50); 

   TFileDirectory dir_dist_beamspot_S_cut_on_decay_vertex = dir_dist_beamspot_S.mkdir("dist_beamspot_S_cut_on_decay_vertex");
   histos_th2f[b+"h_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_dist_beamspot_S_cut_on_decay_vertex.make<TH2F>(b+"h_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; vx (distance beamspot to S annihilation vertex); vy (distance beamspot to S annihilation vertex)",1000, -50, 50,1000, -50, 50); 
   histos_th2f[b+"h_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm_delta_phi_between_1_and_2.5"]= dir_dist_beamspot_S_cut_on_decay_vertex.make<TH2F>(b+"h_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm_delta_phi_between_1_and_2.5",b+"h_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm_delta_phi_between_1_and_2.5; vx (distance beamspot to S annihilation vertex); vy (distance beamspot to S annihilation vertex)",1000, -50, 50,1000, -50, 50); 

   //for distance between beamspot and anti-S candidate vertex 
   TFileDirectory dir_dist_beamspot_anti_S = dir_LambdaKshortVertexFilter.mkdir("dist_beamspot_anti_S_BLINDED");
   histos_th1f[b+"h_anti_S_vtx_distance_to_beamspot"] = dir_dist_beamspot_anti_S.make<TH1F>(b+"h_anti_S_vtx_distance_to_beamspot",b+"h_anti_S_vtx_distance_to_beamspot; lxy(anti S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th1f[b+"h_anti_S_vtx_distance_to_beamspot_error"] = dir_dist_beamspot_S.make<TH1F>(b+"h_anti_S_vtx_distance_to_beamspot_error",b+"h_anti_S_vtx_distance_to_beamspot_error; error lxy(S interaction vertex,  beamspot) [cm];",3000,0,30);
   histos_th2f[b+"h_delta_phi_anti_S_vtx_distance_to_beamspot"]= dir_dist_beamspot_anti_S.make<TH2F>(b+"h_delta_phi_anti_S_vtx_distance_to_beamspot",b+"h_delta_phi_anti_S_vtx_distance_to_beamspot; #Delta #Phi [rad]; lxy(anti S interaction vertex,  beamspot) [cm];",4000, -4, 4,3000,0,30); 
   histos_th2f[b+"h_anti_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error"]= dir_dist_beamspot_anti_S.make<TH2F>(b+"h_anti_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error",b+"h_anti_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error; lxy(anti S interaction vertex,  beamspot) [cm]; errror lxy(anti S interaction vertex,  beamspot) [cm];",3000, 0, 30,3000,0,30); 
   histos_th2f[b+"h_delta_phi_anti_S_vtx_distance_to_beamspot_significance"]= dir_dist_beamspot_anti_S.make<TH2F>(b+"h_delta_phi_anti_S_vtx_distance_to_beamspot_significance",b+"h_delta_phi_anti_S_vtx_distance_to_beamspot_significance; #Delta #Phi [rad]; lxy(anti S interaction vertex,  beamspot) significance;",4000, -4, 4,3000,0,30); 
   histos_th2f[b+"h_anti_S_vtx_distance_to_beamspot_vx_vy"]= dir_dist_beamspot_anti_S.make<TH2F>(b+"h_anti_S_vtx_distance_to_beamspot_vx_vy",b+"h_anti_S_vtx_distance_to_beamspot_vx_vy; vx (distance beamspot to anti S annihilation vertex); vy (distance beamspot to anti S annihilation vertex)",1000, -50, 50,1000, -50, 50); 

   TFileDirectory dir_dist_beamspot_anti_S_cut_on_decay_vertex = dir_dist_beamspot_anti_S.mkdir("dist_beamspot_anti_S_cut_on_decay_vertex");
   histos_th2f[b+"h_anti_S_vtx_distance_to_beamspot_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_dist_beamspot_anti_S_cut_on_decay_vertex.make<TH2F>(b+"h_anti_S_vtx_distance_to_beamspot_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_anti_S_vtx_distance_to_beamspot_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; vx (distance beamspot to anti S annihilation vertex); vy (distance beamspot to anti S annihilation vertex)",1000, -50, 50,1000, -50, 50); 
   histos_th2f[b+"h_anti_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm_delta_phi_between_1_and_2.5"]= dir_dist_beamspot_anti_S_cut_on_decay_vertex.make<TH2F>(b+"h_anti_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm_delta_phi_between_1_and_2.5",b+"h_anti_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm_delta_phi_between_1_and_2.5; vx (distance beamspot to anti S annihilation vertex); vy (distance beamspot to anti S annihilation vertex)",1000, -50, 50,1000, -50, 50); 


   //for angular coorelation between the daughters of the S candidates
   TFileDirectory dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S = dir_LambdaKshortVertexFilter.mkdir("Angular_Correlation_Daughters_S");
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_phi"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S.make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_phi",b+"h_Sdaughters_L0_Ks_delta_phi; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad];",1000,-4,4); 
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_eta"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S.make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_eta",b+"h_Sdaughters_L0_Ks_delta_eta; delta eta;",4000,-10,10); 
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_R"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S.make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_R",b+"h_Sdaughters_L0_Ks_delta_R; delta R;",1000,0,10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 

   TFileDirectory dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex = dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S.mkdir("Angular_Correlation_Daughters_S_cut_on_decay_vertices");
   			
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_phi_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_phi_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_Sdaughters_L0_Ks_delta_phi_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad];",1000,-4,4); 
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; h_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm;",4000,-10,10); 
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_R_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_R_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_Sdaughters_L0_Ks_delta_R_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; delta R;",1000,0,10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_0.5cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_0.5cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_0.5cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.0cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.0cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.0cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.5cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.5cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.5cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.01cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.01cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.02cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.02cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.02cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.03cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.03cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.03cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.04cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.04cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.04cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.06cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.06cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.06cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.08cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.08cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.08cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.1cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.1cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.1cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.04cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.04cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.04cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_0.5cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_0.5cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_0.5cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.0cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.0cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.0cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.5cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.5cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.5cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.9cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.9cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.9cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.01cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.01cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.02cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.02cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.02cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.03cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.03cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.03cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.04cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.04cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.04cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.06cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.06cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.06cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.08cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.08cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.08cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.1cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.1cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.1cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.9cm_and_error_larger_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S_cut_decay_vertex.make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.9cm_and_error_larger_0.01cm",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.9cm_and_error_larger_0.01cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 



   TFileDirectory dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation = dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S.mkdir("Angular_Correlation_Daughters_S_kinematics_S_in_correlation");
 
   histos_th1f[b+"h_s_candidates_mass_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_mass_in_delta_phi_delta_eta_corr",b+"h_s_candidates_mass_in_delta_phi_delta_eta_corr; mass (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_pt_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_pt_in_delta_phi_delta_eta_corr",b+"h_s_candidates_pt_in_delta_phi_delta_eta_corr; pt (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_p_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_p_in_delta_phi_delta_eta_corr",b+"h_s_candidates_p_in_delta_phi_delta_eta_corr; p (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_energy_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_energy_in_delta_phi_delta_eta_corr",b+"h_s_candidates_energy_in_delta_phi_delta_eta_corr; energy (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_et_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_et_in_delta_phi_delta_eta_corr",b+"h_s_candidates_et_in_delta_phi_delta_eta_corr; et (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_mt_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_mt_in_delta_phi_delta_eta_corr",b+"h_s_candidates_mt_in_delta_phi_delta_eta_corr; mt (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_phi_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_phi_in_delta_phi_delta_eta_corr",b+"h_s_candidates_phi_in_delta_phi_delta_eta_corr; #Phi (rad);",4000,-4,4); 
   histos_th1f[b+"h_s_candidates_theta_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_theta_in_delta_phi_delta_eta_corr",b+"h_s_candidates_theta_in_delta_phi_delta_eta_corr; #theta (rad);",4000,-4,4); 
   histos_th1f[b+"h_s_candidates_eta_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_eta_in_delta_phi_delta_eta_corr",b+"h_s_candidates_eta_in_delta_phi_delta_eta_corr; #eta (rad);",4000,-4,50); 
   histos_th1f[b+"h_s_candidates_vx_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_vx_in_delta_phi_delta_eta_corr",b+"h_s_candidates_vx_in_delta_phi_delta_eta_corr; vx (cm);",1000,-50,50); 
   histos_th1f[b+"h_s_candidates_vy_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_vy_in_delta_phi_delta_eta_corr",b+"h_s_candidates_vy_in_delta_phi_delta_eta_corr; vy (GeV);",1000,-50,50); 
   histos_th1f[b+"h_s_candidates_vz_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_vz_in_delta_phi_delta_eta_corr",b+"h_s_candidates_vz_in_delta_phi_delta_eta_corr; vz (GeV);",1000,-50,50); 
   histos_th1f[b+"h_s_candidates_vertexNormalizedChi2_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_vertexNormalizedChi2_in_delta_phi_delta_eta_corr",b+"h_s_candidates_vertexNormalizedChi2_in_delta_phi_delta_eta_corr; Chi2/ndof;",1000,0,50); 
   histos_th1f[b+"h_s_candidates_lxy_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_lxy_in_delta_phi_delta_eta_corr",b+"h_s_candidates_lxy_in_delta_phi_delta_eta_corr; lxy (cm);",1000,0,20); 
   histos_th1f[b+"h_s_candidates_lxy_signed_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_lxy_signed_in_delta_phi_delta_eta_corr",b+"h_s_candidates_lxy_signed_in_delta_phi_delta_eta_corr; signed lxy (cm);",1000,-20,20); 
   histos_th1f[b+"h_s_candidates_lxy_signed_in_delta_phi_delta_eta_corr_error_lxy_smaller_0p1"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_lxy_signed_in_delta_phi_delta_eta_corr_error_lxy_smaller_0p1",b+"h_s_candidates_lxy_signed_in_delta_phi_delta_eta_corr_error_lxy_smaller_0p1; signed lxy (cm);",1000,-20,20); 
   histos_th1f[b+"h_s_candidates_error_lxy_in_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_correlation.make<TH1F>(b+"h_s_candidates_error_lxy_in_delta_phi_delta_eta_corr",b+"h_s_candidates_error_lxy_in_delta_phi_delta_eta_corr; error lxy (cm);",1000,0,20); 
   TFileDirectory dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak = dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S.mkdir("Angular_Correlation_Daughters_S_kinematics_S_in_peak");
   histos_th1f[b+"h_s_candidates_mass_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_mass_in_peak",b+"h_s_candidates_mass_in_peak; Mass (GeV);",1000,0,20); 
   histos_th1f[b+"h_s_candidates_pt_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_pt_in_peak",b+"h_s_candidates_pt_in_peak; pT (GeV);",1000,0,20); 
   histos_th1f[b+"h_s_candidates_phi_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_phi_in_peak",b+"h_s_candidates_phi_in_peak; #Phi(Rad);",1000,-4,4); 
   histos_th1f[b+"h_s_candidates_eta_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_eta_in_peak",b+"h_s_candidates_eta_in_peak;#eta ;",1000,-4,4); 
   histos_th1f[b+"h_s_candidates_vx_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_vx_in_peak",b+";h_s_candidates_vx_in_peak ;vx(cm)",1000,-20,20); 
   histos_th1f[b+"h_s_candidates_vy_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_vy_in_peak",b+";h_s_candidates_vy_in_peak ;vy(cm)",1000,-20,20); 
   histos_th1f[b+"h_s_candidates_vz_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_vz_in_peak",b+";h_s_candidates_vz_in_peak ;vz(cm)",1000,-20,20); 
   histos_th1f[b+"h_s_candidates_lxy_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_lxy_in_peak",b+"h_s_candidates_lxy_in_peak; lxy (cm);",1000,0,20); 
   histos_th1f[b+"h_s_candidates_lxy_signed_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_lxy_signed_in_peak",b+"h_s_candidates_lxy_signed_in_peak;lxy signed (cm) ;",1000,-20,20); 
   histos_th1f[b+"h_s_candidates_lxy_signed_in_peak_error_lxy_smaller_0p1"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_lxy_signed_in_peak_error_lxy_smaller_0p1",b+"h_s_candidates_lxy_signed_in_peak_error_lxy_smaller_0p1;lxy signed (cm) ;",1000,-20,20); 
   histos_th1f[b+"h_s_candidates_error_lxy_in_peak"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_in_peak.make<TH1F>(b+"h_s_candidates_error_lxy_in_peak",b+"h_s_candidates_error_lxy_in_peak; error lxy (cm);",1000,-20,20); 

   TFileDirectory dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation = dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_S.mkdir("Angular_Correlation_Daughters_S_kinematics_S_outside_correlation");
 
   histos_th1f[b+"h_s_candidates_mass_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_mass_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_mass_outside_delta_phi_delta_eta_corr; mass (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_pt_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_pt_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_pt_outside_delta_phi_delta_eta_corr; pt (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_p_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_p_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_p_outside_delta_phi_delta_eta_corr; p (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_energy_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_energy_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_energy_outside_delta_phi_delta_eta_corr; energy (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_et_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_et_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_et_outside_delta_phi_delta_eta_corr; et (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_mt_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_mt_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_mt_outside_delta_phi_delta_eta_corr; mt (GeV);",1000,0,50); 
   histos_th1f[b+"h_s_candidates_phi_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_phi_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_phi_outside_delta_phi_delta_eta_corr; #Phi (rad);",4000,-4,4); 
   histos_th1f[b+"h_s_candidates_theta_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_theta_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_theta_outside_delta_phi_delta_eta_corr; #theta (rad);",4000,-4,4); 
   histos_th1f[b+"h_s_candidates_eta_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_eta_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_eta_outside_delta_phi_delta_eta_corr; #eta (rad);",4000,-10,10); 
   histos_th1f[b+"h_s_candidates_vx_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_vx_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_vx_outside_delta_phi_delta_eta_corr; vx (cm);",1000,-50,50); 
   histos_th1f[b+"h_s_candidates_vy_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_vy_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_vy_outside_delta_phi_delta_eta_corr; vy (GeV);",1000,-50,50); 
   histos_th1f[b+"h_s_candidates_vz_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_vz_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_vz_outside_delta_phi_delta_eta_corr; vz (GeV);",1000,-50,50); 
   histos_th1f[b+"h_s_candidates_vertexNormalizedChi2_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_vertexNormalizedChi2_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_vertexNormalizedChi2_outside_delta_phi_delta_eta_corr; Chi2/ndof;",1000,0,50); 
   histos_th1f[b+"h_s_candidates_lxy_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_lxy_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_lxy_outside_delta_phi_delta_eta_corr; lxy (cm);",1000,0,20); 
   histos_th1f[b+"h_s_candidates_lxy_signed_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_lxy_signed_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_lxy_signed_outside_delta_phi_delta_eta_corr; lxy signed (cm);",1000,-20,20); 
   histos_th1f[b+"h_s_candidates_lxy_signed_outside_delta_phi_delta_eta_corr_error_lxy_smaller_0p1"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_lxy_signed_outside_delta_phi_delta_eta_corr_error_lxy_smaller_0p1",b+"h_s_candidates_lxy_signed_outside_delta_phi_delta_eta_corr_error_lxy_smaller_0p1; lxy signed (cm);",1000,-20,20); 
   histos_th1f[b+"h_s_candidates_error_lxy_outside_delta_phi_delta_eta_corr"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_kinematics_S_outside_correlation.make<TH1F>(b+"h_s_candidates_error_lxy_outside_delta_phi_delta_eta_corr",b+"h_s_candidates_error_lxy_outside_delta_phi_delta_eta_corr; error lxy (cm);",1000,0,20); 


   //for angular coorelation between the daughters of the anti-S candidates
   TFileDirectory dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S = dir_LambdaKshortVertexFilter.mkdir("Angular_Correlation_Daughters_anti_S_BLINDED!!!!");
   histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_phi"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S.make<TH1F>(b+"h_anti_Sdaughters_L0_Ks_delta_phi",b+"h_anti_Sdaughters_L0_Ks_delta_phi; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad];",1000,-4,4); 
   histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_eta"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S.make<TH1F>(b+"h_anti_Sdaughters_L0_Ks_delta_eta",b+"h_anti_Sdaughters_L0_Ks_delta_eta; delta eta;",4000,-10,10); 
   histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_R"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S.make<TH1F>(b+"h_anti_Sdaughters_L0_Ks_delta_R",b+"h_anti_Sdaughters_L0_Ks_delta_R; delta R;",1000,0,10); 
   histos_th2f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S.make<TH2F>(b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta",b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S.make<TH2F>(b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent",b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 

   TFileDirectory dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S_cut_decay_vertex = dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S.mkdir("Angular_Correlation_Daughters_anti_S_cut_on_decay_vertices_BLINDED!!!!");
   			
   histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S_cut_decay_vertex.make<TH1F>(b+"h_anti_Sdaughters_L0_Ks_delta_phi_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_anti_Sdaughters_L0_Ks_delta_phi_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad];",1000,-4,4); 
   histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S_cut_decay_vertex.make<TH1F>(b+"h_anti_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_anti_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; h_anti_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm;",4000,-10,10); 
   histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_R_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S_cut_decay_vertex.make<TH1F>(b+"h_anti_Sdaughters_L0_Ks_delta_R_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_anti_Sdaughters_L0_Ks_delta_R_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; delta R;",1000,0,10); 
   histos_th2f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S_cut_decay_vertex.make<TH2F>(b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 
   histos_th2f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]= dir_LambdaKshortVertexFilter_Angular_Correlation_Daughters_anti_S_cut_decay_vertex.make<TH2F>(b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm",b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta eta",1000,-4,4, 4000, -10, 10); 

			

   //for angular correlation between Ks and L0 produced in the LambdaKshortFilter
   TFileDirectory dir_LambdaKshortFilter = m_fs->mkdir("LambdaKshortFilter");
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_dz_Filters_Collection"]= dir_LambdaKshortFilter.make<TH2F>(b+"h_L0_Ks_delta_phi_delta_dz_Filters_Collection",b+"h_L0_Ks_delta_phi_delta_dz_Filters_Collection; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; #Delta(dz(K_{S}^{0}),dz(#Lambda^{0}))",1000,-4,4, 1000, -50,50); 
   histos_th1f[b+"h_L0_Ks_delta_phi_Filters_Collection"]= dir_LambdaKshortFilter.make<TH1F>(b+"h_L0_Ks_delta_phi_Filters_Collection",b+"h_L0_Ks_delta_phi_Filters_Collection; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]",1000,-4,4); 
   histos_th1f[b+"h_L0_Ks_delta_eta_Filters_Collection"]= dir_LambdaKshortFilter.make<TH1F>(b+"h_L0_Ks_delta_eta_Filters_Collection",b+"h_L0_Ks_delta_eta_Filters_Collection; delta eta",4000,-10,10); 
   histos_th1f[b+"h_L0_Ks_delta_R_Filters_Collection"]= dir_LambdaKshortFilter.make<TH1F>(b+"h_L0_Ks_delta_R_Filters_Collection",b+"h_L0_Ks_delta_R_Filters_Collection; delta R",3000,0,30); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection"]= dir_LambdaKshortFilter.make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection",b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta_eta",100,-4,4, 100, -10, 10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection_no_cent"]= dir_LambdaKshortFilter.make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection_no_cent",b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection_no_cent; #Delta #Phi (K_{S}^{0},#Lambda^{0}) [rad]; delta_eta",100,-4,4, 100, -10, 10); 

   
   //for the s candidates after the massfilter
   TFileDirectory dir_sMassFilter = m_fs->mkdir("sMassFilter"); 
   TFileDirectory dir_s_candidates = dir_sMassFilter.mkdir("S_candidates");
   histos_th1f[b+"h_s_candidates_mass_after_massFilter"]= dir_s_candidates.make<TH1F>(b+"h_s_candidates_mass_after_massFilter",b+"h_s_candidates_mass_after_massFilter; S candidate mass (GeV);",40000,-200,200); 
   histos_th1f[b+"h_s_candidates_pt_after_massFilter"] =  dir_s_candidates.make<TH1F>(b+"h_s_candidates_pt_after_massFilter",b+"h_s_candidates_pt_after_massFilter; S candidate pt (GeV);",40000,-200,200);
   histos_th1f[b+"h_s_candidates_p_after_massFilter"] =  dir_s_candidates.make<TH1F>(b+"h_s_candidates_p_after_massFilter",b+"h_s_candidates_p_after_massFilter; S candidate p (GeV);",40000,-200,200);
   histos_th1f[b+"h_s_candidates_energy_after_massFilter"]  =  dir_s_candidates.make<TH1F>(b+"h_s_candidates_energy_after_massFilter",b+"h_s_candidates_energy_after_massFilter; S candidate E (GeV);",20000,0,200);
   histos_th1f[b+"h_s_candidates_et_after_massFilter"] =  dir_s_candidates.make<TH1F>(b+"h_s_candidates_et_after_massFilter",b+"h_s_candidates_et_after_massFilter; S candidate Et (GeV);",20000,0,200);
   histos_th1f[b+"h_s_candidates_mt_after_massFilter"] = dir_s_candidates.make<TH1F>(b+"h_s_candidates_mt_after_massFilter",b+"h_s_candidates_mt_after_massFilter; S candidate mt (GeV);",40000,-200,200);
   histos_th1f[b+"h_s_candidates_phi_after_massFilter"] = dir_s_candidates.make<TH1F>(b+"h_s_candidates_phi_after_massFilter",b+"h_s_candidates_phi_after_massFilter; S candidate phi (rad);",8000,-4,4);
   histos_th1f[b+"h_s_candidates_theta_after_massFilter"] = dir_s_candidates.make<TH1F>(b+"h_s_candidates_theta_after_massFilter",b+"h_s_candidates_theta_after_massFilter; S candidate theta (rad);",8000,-4,4);
   histos_th1f[b+"h_s_candidates_eta_after_massFilter"] = dir_s_candidates.make<TH1F>(b+"h_s_candidates_eta_after_massFilter",b+"h_s_candidates_eta_after_massFilter; S candidate eta;",800,-10,10);
   histos_th1f[b+"h_s_candidates_vx_after_massFilter"] = dir_s_candidates.make<TH1F>(b+"h_s_candidates_vx_after_massFilter",b+"h_s_candidates_vx_after_massFilter; S candidate vx (cm);",2000,-100,100);
   histos_th1f[b+"h_s_candidates_vy_after_massFilter"] = dir_s_candidates.make<TH1F>(b+"h_s_candidates_vy_after_massFilter",b+"h_s_candidates_vy_after_massFilter; S candidate vy (cm);",2000,-100,100);
   histos_th1f[b+"h_s_candidates_vz_after_massFilter"] = dir_s_candidates.make<TH1F>(b+"h_s_candidates_vz_after_massFilter",b+"h_s_candidates_vz_after_massFilter; S candidate vz (cm);",2000,-100,100);
   histos_th1f[b+"h_s_candidates_vertexNormalizedChi2_after_massFilter"] = dir_s_candidates.make<TH1F>(b+"h_s_candidates_vertexNormalizedChi2_after_massFilter",b+"h_s_candidates_vertexNormalizedChi2_after_massFilter; S candidate vertex normalised Chi2;",1000,0,100);
   histos_th1f[b+"h_s_candidates_vertexNdof_after_massFilter"] = dir_s_candidates.make<TH1F>(b+"h_s_candidates_vertexNdof_after_massFilter",b+"h_s_candidates_vertexNdof_after_massFilter; S candidate vertex Ndof;",100,0,100);

   TFileDirectory dir_anti_s_candidates = dir_sMassFilter.mkdir("Anti-S_candidates_(BLINDED!!!!!)");
   histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter"]= dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_mass_after_massFilter",b+"h_anti_s_candidates_mass_after_massFilter; anti S candidate mass (GeV);",40000,-200,200); 
   histos_th1f[b+"h_anti_s_candidates_pt_after_massFilter"] =  dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_pt_after_massFilter",b+"h_anti_s_candidates_pt_after_massFilter; anti S candidate pt (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_p_after_massFilter"] =  dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_p_after_massFilter",b+"h_anti_s_candidates_p_after_massFilter; anti S candidate p (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_energy_after_massFilter"]  =  dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_energy_after_massFilter",b+"h_anti_s_candidates_energy_after_massFilter; anti S candidate E (GeV);",20000,0,200);
   histos_th1f[b+"h_anti_s_candidates_et_after_massFilter"] =  dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_et_after_massFilter",b+"h_anti_s_candidates_et_after_massFilter; anti S candidate Et (GeV);",20000,0,200);
   histos_th1f[b+"h_anti_s_candidates_mt_after_massFilter"] = dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_mt_after_massFilter",b+"h_anti_s_candidates_mt_after_massFilter; anti S candidate mt (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_phi_after_massFilter"] = dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_phi_after_massFilter",b+"h_anti_s_candidates_phi_after_massFilter; anti S candidate phi (rad);",8000,-4,4);
   histos_th1f[b+"h_anti_s_candidates_theta_after_massFilter"] = dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_theta_after_massFilter",b+"h_anti_s_candidates_theta_after_massFilter; anti S candidate theta (rad);",8000,-4,4);
   histos_th1f[b+"h_anti_s_candidates_eta_after_massFilter"] = dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_eta_after_massFilter",b+"h_anti_s_candidates_eta_after_massFilter; anti S candidate eta;",800,-10,10);
   histos_th1f[b+"h_anti_s_candidates_vx_after_massFilter"] = dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_vx_after_massFilter",b+"h_anti_s_candidates_vx_after_massFilter; anti S candidate vx (cm);",2000,-100,100);
   histos_th1f[b+"h_anti_s_candidates_vy_after_massFilter"] = dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_vy_after_massFilter",b+"h_anti_s_candidates_vy_after_massFilter; anti S candidate vy (cm);",2000,-100,100);
   histos_th1f[b+"h_anti_s_candidates_vz_after_massFilter"] = dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_vz_after_massFilter",b+"h_anti_s_candidates_vz_after_massFilter; anti S candidate vz (cm);",2000,-100,100);
   histos_th1f[b+"h_anti_s_candidates_vertexNormalizedChi2_after_massFilter"] = dir_anti_s_candidates.make<TH1F>(b+"h_anti_s_candidates_vertexNormalizedChi2_after_massFilter",b+"h_anti_s_candidates_vertexNormalizedChi2_after_massFilter; anti S candidate vertex normalised Chi2;",1000,0,100);
   //for the comparison between the reconstructed S particle and the gen particle

   TFileDirectory dir_comparison_generated = dir_sMassFilter.mkdir("Comparison_with_generator_particles");
   
   histos_th1f[b+"h_genPartices_Scands_pdgId"] = dir_comparison_generated.make<TH1F>(b+"h_genPartices_Scands_pdgId",b+"h_genPartices_Scands_pdgId; pdgId;",30000,-15000,15000);
   histos_th1f[b+"h_genPartices_Scands_mass"] = dir_comparison_generated.make<TH1F>(b+"h_genPartices_Scands_mass",b+"h_genPartices_Scands_mass; mass [GeV];",2000,0,20);
   histos_th1f[b+"h_genPartices_Scands_diff_pt"] = dir_comparison_generated.make<TH1F>(b+"h_genPartices_Scands_diff_pt",b+"h_genPartices_Scands_diff_pt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Scands_diff_mt"] = dir_comparison_generated.make<TH1F>(b+"h_genPartices_Scands_diff_mt",b+"h_genPartices_Scands_diff_mt; diff(pt) (gen-reco)[GeV];",2000,-10,10);
   histos_th1f[b+"h_genPartices_Scands_diff_phi"] = dir_comparison_generated.make<TH1F>(b+"h_genPartices_Scands_diff_phi",b+"h_genPartices_Scands_diff_phi; diff(phi) (gen-reco)[rad];",8000,-4,4);
   histos_th1f[b+"h_genPartices_Scands_diff_eta"] = dir_comparison_generated.make<TH1F>(b+"h_genPartices_Scands_diff_eta",b+"h_genPartices_Scands_diff_eta; diff(eta) (gen-reco);",20000,-10,10);
   histos_th1f[b+"h_genPartices_Scands_diff_vx"] = dir_comparison_generated.make<TH1F>(b+"h_genPartices_Scands_diff_vx",b+"h_genPartices_Scands_diff_vx; diff(vx) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Scands_diff_vy"] = dir_comparison_generated.make<TH1F>(b+"h_genPartices_Scands_diff_vy",b+"h_genPartices_Scands_diff_vy; diff(vy) (gen-reco)[cm];",10000,-50,50);
   histos_th1f[b+"h_genPartices_Scands_diff_vz"] = dir_comparison_generated.make<TH1F>(b+"h_genPartices_Scands_diff_vz",b+"h_genPartices_Scands_diff_vz; diff(vz) (gen-reco)[cm];",10000,-50,50);
   histos_th2f[b+"h_genPartices_Scands_pt_diff_vx"] = dir_comparison_generated.make<TH2F>(b+"h_genPartices_Scands_pt_diff_vx",b+"h_genPartices_Scands_pt_diff_vx;gen particle pt; diff(vx) (gen-reco)[cm];",1000,0,10,10000,-50,50);
   histos_th2f[b+"h_genPartices_Scands_pt_diff_vy"] = dir_comparison_generated.make<TH2F>(b+"h_genPartices_Scands_pt_diff_vy",b+"h_genPartices_Scands_pt_diff_vy;gen particle pt; diff(vy) (gen-reco)[cm];",1000,0,10,10000,-50,50);

   //for the daughters of the s candidates after the massfilter
/*
   histos_th1f[b+"h_s_candidates_mass_after_massFilter_Lambda_mass"] = m_fs->make<TH1F>(b+"h_s_candidates_mass_after_massFilter_Lambda_mass",b+"h_s_candidates_mass_after_massFilter_Lambda_mass; Lambda mass in the S decay chain (GeV);",40000,-200,200);
   histos_th1f[b+"h_s_candidates_mass_after_massFilter_Kshort_mass"] = m_fs->make<TH1F>(b+"h_s_candidates_mass_after_massFilter_Kshort_mass",b+"h_s_candidates_mass_after_massFilter_Kshort_mass; Kshort mass in the S decay chain (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter_Lambda_mass"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_mass_after_massFilter_Lambda_mass",b+"h_anti_s_candidates_mass_after_massFilter_Lambda_mass; Lambda mass in the anti S decay chain (GeV);",40000,-200,200);
   histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter_Kshort_mass"] = m_fs->make<TH1F>(b+"h_anti_s_candidates_mass_after_massFilter_Kshort_mass",b+"h_anti_s_candidates_mass_after_massFilter_Kshort_mass; Kshort mass in the anti S decay chain (GeV);",40000,-200,200);
*/
   //for the r candidates after the massfilter
    
   TFileDirectory dir_rMassFilter = m_fs->mkdir("rMassFilter"); 
   histos_th1f[b+"h_all_r_candidates_mass_after_massFilter"]= dir_rMassFilter.make<TH1F>(b+"h_all_r_candidates_mass_after_massFilter",b+"h_all_r_candidates_mass_after_massFilter; all r candidate mass (GeV);",40000,-200,200); 
   histos_th1f[b+"h_r_candidates_mass_after_massFilter"]= dir_rMassFilter.make<TH1F>(b+"h_r_candidates_mass_after_massFilter",b+"h_r_candidates_mass_after_massFilter; r candidate mass (GeV);",40000,-200,200); 
   histos_th1f[b+"h_anti_r_candidates_mass_after_massFilter"]= dir_rMassFilter.make<TH1F>(b+"h_anti_r_candidates_mass_after_massFilter",b+"h_anti_r_candidates_mass_after_massFilter; anti r candidate mass (GeV);",40000,-200,200); 
    
}


void Analyzer_V0_angular_correlation::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


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


 // edm::Handle <vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > > h_TKMET;
 // iEvent.getByToken(m_TKMETToken, h_TKMET);

  edm::Handle <vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >h_MET;
  iEvent.getByToken(m_METToken, h_MET);

  edm::Handle <vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > h_TwoTopJets;
  iEvent.getByToken(m_TwoTopJetsToken, h_TwoTopJets);

  edm::Handle <vector<double> > h_HT;
  iEvent.getByToken(m_HTToken, h_HT);

  edm::Handle <vector<double> > h_TKHT;
  iEvent.getByToken(m_TKHTToken, h_TKHT);

  //gen particles 
  edm::Handle<vector<reco::GenParticle>> h_genParticles;
  iEvent.getByToken(m_genParticlesToken, h_genParticles);
  
  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

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

 if(h_nPVs.isValid() && h_nPVs->size() > 0)  		 	histos_th1f[b+"h_h_nPVs"]->Fill(h_nPVs->at(0));
 if(h_nelectrons.isValid() && h_nelectrons->size() > 0 )	histos_th1f[b+"h_h_nelectrons"]->Fill(h_nelectrons->at(0));
 if(h_njets.isValid() && h_njets->size() > 0)			histos_th1f[b+"h_h_njets"]->Fill(h_njets->at(0));
 if(h_nkshorts.isValid() && h_nkshorts->size() > 0)		histos_th1f[b+"h_h_nkshorts"]->Fill(h_nkshorts->at(0));
 if(h_nlambdas.isValid() && h_nlambdas->size() > 0)		histos_th1f[b+"h_h_nlambdas"]->Fill(h_nlambdas->at(0));
 if(h_nmuons.isValid() && h_nmuons->size() > 0)			histos_th1f[b+"h_h_nmuons"]->Fill(h_nmuons->at(0));
// if(h_ntracks.isValid() && h_ntracks->size() > 0)		histos_th1f[b+"h_h_tracks"]->Fill(h_ntracks->at(0));

 if(h_MET.isValid() && h_MET->size() > 0)			histos_th1f[b+"h_h_MET_pT"]->Fill(h_MET->at(0).pt());
 if(h_TwoTopJets.isValid() && h_TwoTopJets->size() > 0)		histos_th1f[b+"h_h_TwoTopJets_pT_0"]->Fill(h_TwoTopJets->at(0).pt()); 
 if(h_TwoTopJets.isValid() && h_TwoTopJets->size() > 0)		histos_th1f[b+"h_h_TwoTopJets_pT_1"]->Fill(h_TwoTopJets->at(1).pt()); 
 if(h_HT.isValid() && h_HT->size() > 0)				histos_th1f[b+"h_h_HT"]->Fill(h_HT->at(0));
 if(h_TKHT.isValid() && h_TKHT->size() > 0)			histos_th1f[b+"h_h_TKHT"]->Fill(h_TKHT->at(0));


  //get info on the beamspot:

  //beamspot
  double bx_x = h_bs->x0(); 
  double bx_y = h_bs->y0(); 
  double bx_z = h_bs->z0(); 
  TVector3 beamspot(bx_x,bx_y,bx_x);
  
  double bx_x_std_dev = h_bs->x0Error(); 
  double bx_y_std_dev = h_bs->y0Error(); 
  double bx_z_std_dev = h_bs->z0Error(); 


  if(h_bs.isValid()){ 	

	double bx_x = h_bs->x0(); 
	double bx_y = h_bs->y0(); 
	double bx_z = h_bs->z0(); 

	double bx_x_std_dev = h_bs->x0Error(); 
	double bx_y_std_dev = h_bs->y0Error(); 
	double bx_z_std_dev = h_bs->z0Error(); 
	
        double lxy_beamspot = sqrt(bx_x*bx_x+bx_y*bx_y);
	double lxy_beamspot_std_dev = std_dev_lxy(bx_x, bx_y, bx_x_std_dev*bx_x_std_dev, bx_y_std_dev*bx_y_std_dev, 0, 0, 0, 0);
	histos_th1f[b+"h_beamspot_vx"]->Fill(bx_x);
	histos_th1f[b+"h_beamspot_vy"]->Fill(bx_y);
	histos_th1f[b+"h_beamspot_vx_std_dev"]->Fill(bx_x_std_dev);
	histos_th1f[b+"h_beamspot_vy_std_dev"]->Fill(bx_y_std_dev);
	histos_th1f[b+"h_beamspot_lxy"]->Fill(lxy_beamspot);
	histos_th1f[b+"h_beamspot_lxy_std_dev"]->Fill(lxy_beamspot_std_dev);

  }
 
  //do stuff with the collections produced in the LambdaKshort filters
  if(h_Lambdas_LambdaKshortFilter.isValid() && h_Kshorts_LambdaKshortFilter.isValid()){
  for(unsigned int l = 0; l < h_Lambdas_LambdaKshortFilter->size(); ++l){

	double phi2 = (*h_Lambdas_LambdaKshortFilter)[l]->phi();
	double eta2 = (*h_Lambdas_LambdaKshortFilter)[l]->eta();


	//vertex	
	double vx_Lambda = h_Lambdas->at(l).vx(); 
	double vy_Lambda = h_Lambdas->at(l).vy(); 
	double vz_Lambda = h_Lambdas->at(l).vz();
	TVector3 Lambda_v(vx_Lambda,vy_Lambda,vz_Lambda);

	//momenta
	double px_Lambda = h_Lambdas->at(l).px(); 
	double py_Lambda = h_Lambdas->at(l).py(); 
	double pz_Lambda = h_Lambdas->at(l).pz();
	TVector3 Lambda_p(px_Lambda,py_Lambda,pz_Lambda);

	
	double z_PCA_Lambda = z_PCA_line_point(Lambda_v, Lambda_p, beamspot);

	for(unsigned int k = 0; k < h_Kshorts_LambdaKshortFilter->size(); ++k){

		//vertex	
		double vx_Kshort = h_Kshorts->at(k).vx(); 
		double vy_Kshort = h_Kshorts->at(k).vy(); 
		double vz_Kshort = h_Kshorts->at(k).vz();
		TVector3 Kshort_v(vx_Kshort,vy_Kshort,vz_Kshort);

		//momenta
		double px_Kshort = h_Kshorts->at(k).px(); 
		double py_Kshort = h_Kshorts->at(k).py(); 
		double pz_Kshort = h_Kshorts->at(k).pz();
		TVector3 Kshort_p(px_Kshort,py_Kshort,pz_Kshort);

		double z_PCA_Kshort = z_PCA_line_point(Kshort_v, Kshort_p, beamspot);

		double delta_z_PCA = fabs(z_PCA_Lambda-z_PCA_Kshort);
	
		double phi1 = (*h_Kshorts_LambdaKshortFilter)[k]->phi();
                double delta_phi = reco::deltaPhi(phi1, phi2);

                double eta1 = (*h_Kshorts_LambdaKshortFilter)[k]->eta();
                double delta_eta = eta1-eta2;

		double delta_R = sqrt(delta_phi*delta_phi+delta_eta*delta_eta);

                histos_th2f[b+"h_L0_Ks_delta_phi_delta_dz_Filters_Collection"]->Fill(delta_phi,delta_z_PCA);

                histos_th1f[b+"h_L0_Ks_delta_phi_Filters_Collection"]->Fill(delta_phi);
                histos_th1f[b+"h_L0_Ks_delta_eta_Filters_Collection"]->Fill(delta_eta);
                histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection"]->Fill(delta_phi,delta_eta);	
                if( abs(delta_phi) > 0.1 || abs(delta_eta) > 0.1 ) histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection_no_cent"]->Fill(delta_phi,delta_eta);	
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
  cout << "original Lambda V0" << endl;
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STANDARD PLOT: ONLY THE CURRENT EVENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //h_lambdas and h_kshorts are the original lambda and kshort collections
  if(h_Lambdas.isValid()){
	for (unsigned int l = 0; l < h_Lambdas->size(); ++l) {
		histos_th1f[b+"h_Lambda_mass"]->Fill(h_Lambdas->at(l).mass());
	
		//vertex	
		double vx_Lambda = h_Lambdas->at(l).vx(); 
		double vy_Lambda = h_Lambdas->at(l).vy(); 
		double vz_Lambda = h_Lambdas->at(l).vz();
		TVector3 Lambda_v(vx_Lambda,vy_Lambda,vz_Lambda);

		//momenta
		double px_Lambda = h_Lambdas->at(l).px(); 
		double py_Lambda = h_Lambdas->at(l).py(); 
		double pz_Lambda = h_Lambdas->at(l).pz();
		TVector3 Lambda_p(px_Lambda,py_Lambda,pz_Lambda);

		//vertex error
		double vx_Lambda_var = h_Lambdas->at(l).vertexCovariance(0,0);
		double vy_Lambda_var = h_Lambdas->at(l).vertexCovariance(1,1);
		double vz_Lambda_var = h_Lambdas->at(l).vertexCovariance(2,2);
	
		//displacement
		double lxy_Lambda_b = lxy(Lambda_v, beamspot);
                double lxy_Lambda_b_std_dev = std_dev_lxy(vx_Lambda, vy_Lambda, vx_Lambda_var, vy_Lambda_var, bx_x, bx_y, pow(bx_x_std_dev,2) , pow(bx_y_std_dev,2));
                double lxy_Lambda_b_signed =  lxy_signed(Lambda_v, beamspot, Lambda_p);
		
	
		histos_th1f[b+"h_Lambda_lxy"]->Fill(lxy_Lambda_b);	
		histos_th1f[b+"h_Lambda_lxy_error"]->Fill(lxy_Lambda_b_std_dev);	
		histos_th1f[b+"h_Lambda_lxy_signed"]->Fill(lxy_Lambda_b_signed);	

		if(h_genParticles.isValid() )
		{
			//the 1th gen particle is always the Labda particle 
			histos_th1f[b+"h_genPartices_Lambda_origV0_pdgId"]->Fill(h_genParticles->at(1).pdgId());
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_pt"]->Fill(h_genParticles->at(1).pt()-h_Lambdas->at(l).pt());
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_mt"]->Fill(h_genParticles->at(1).mt()-h_Lambdas->at(l).mt());
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_phi"]->Fill(h_genParticles->at(1).phi()-h_Lambdas->at(l).phi());
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_eta"]->Fill(h_genParticles->at(1).eta()-h_Lambdas->at(l).eta());


			double diff_vx = h_genParticles->at(1).daughter(0)->vx()-h_Lambdas->at(l).vx();
			double diff_vy = h_genParticles->at(1).daughter(0)->vy()-h_Lambdas->at(l).vy();
			double diff_vz = h_genParticles->at(1).daughter(0)->vz()-h_Lambdas->at(l).vz();

			double error_vx = h_Lambdas->at(l).vertexCovariance(0,0);
			double error_vy = h_Lambdas->at(l).vertexCovariance(1,1);
			
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vx"]->Fill(diff_vx);
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vy"]->Fill(diff_vy);
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vz"]->Fill(diff_vz);

			histos_th2f[b+"h_genPartices_Lambda_origV0_pt_diff_vx"]->Fill(h_genParticles->at(1).pt() , diff_vx);
			histos_th2f[b+"h_genPartices_Lambda_origV0_diff_vx_error_vx"]->Fill(abs(diff_vx), h_Lambdas->at(l).vertexCovariance(0,0));
			
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vx"]->Fill(diff_vx);
			if(error_vx<0.01) histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vx_error_vx_smaller_0.01"]->Fill(diff_vx);
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vy"]->Fill(diff_vy);
			if(error_vy<0.01) histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vy_error_vy_smaller_0.01"]->Fill(diff_vy);
			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vz"]->Fill(diff_vz);

			histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vxy"]->Fill(sqrt(diff_vx*diff_vx+diff_vy*diff_vy));
			if(error_vx<0.01 && error_vy<0.01)histos_th1f[b+"h_genPartices_Lambda_origV0_diff_vxy_error_vx_and_vy_smaller_0.01"]->Fill(sqrt(diff_vx*diff_vx+diff_vy*diff_vy));
			
		}	

	}
  }
  cout << "original Kshort V0" << endl;

  if(h_Kshorts.isValid()){
	for (unsigned int k = 0; k < h_Kshorts->size(); ++k) {
		histos_th1f[b+"h_Kshort_mass"]->Fill(h_Kshorts->at(k).mass());
		
	
		//vertex	
		double vx_Kshort = h_Kshorts->at(k).vx(); 
		double vy_Kshort = h_Kshorts->at(k).vy(); 
		double vz_Kshort = h_Kshorts->at(k).vz();
		TVector3 Kshort_v(vx_Kshort,vy_Kshort,vz_Kshort);

		//momenta
		double px_Kshort = h_Kshorts->at(k).px(); 
		double py_Kshort = h_Kshorts->at(k).py(); 
		double pz_Kshort = h_Kshorts->at(k).pz();
		TVector3 Kshort_p(px_Kshort,py_Kshort,pz_Kshort);

		//vertex error
		double vx_Kshort_var = h_Kshorts->at(k).vertexCovariance(0,0);
		double vy_Kshort_var = h_Kshorts->at(k).vertexCovariance(1,1);
		double vz_Kshort_var = h_Kshorts->at(k).vertexCovariance(2,2);
	
		//displacement
		double lxy_Kshort_b = lxy(Kshort_v, beamspot);
                double lxy_Kshort_b_std_dev = std_dev_lxy(vx_Kshort, vy_Kshort, vx_Kshort_var, vy_Kshort_var, bx_x, bx_y, pow(bx_x_std_dev,2) , pow(bx_y_std_dev,2));
                double lxy_Kshort_b_signed =  lxy_signed(Kshort_v, beamspot, Kshort_p);
		
	
		histos_th1f[b+"h_Kshort_lxy"]->Fill(lxy_Kshort_b);	
		histos_th1f[b+"h_Kshort_lxy_error"]->Fill(lxy_Kshort_b_std_dev);	
		histos_th1f[b+"h_Kshort_lxy_signed"]->Fill(lxy_Kshort_b_signed);	

		if(h_genParticles.isValid() )
		{
			//the 2th gen particle is always the Kshort particle which corresponds to the S here
			histos_th1f[b+"h_genPartices_Kshort_origV0_pdgId"]->Fill(h_genParticles->at(2).pdgId());
			histos_th1f[b+"h_genPartices_Kshort_origV0_diff_pt"]->Fill(h_genParticles->at(2).pt()-h_Kshorts->at(k).pt());
			histos_th1f[b+"h_genPartices_Kshort_origV0_diff_mt"]->Fill(h_genParticles->at(2).mt()-h_Kshorts->at(k).mt());
			histos_th1f[b+"h_genPartices_Kshort_origV0_diff_phi"]->Fill(h_genParticles->at(2).phi()-h_Kshorts->at(k).phi());
			histos_th1f[b+"h_genPartices_Kshort_origV0_diff_eta"]->Fill(h_genParticles->at(2).eta()-h_Kshorts->at(k).eta());

			double diff_vx = h_genParticles->at(2).daughter(0)->vx()-h_Kshorts->at(k).vx();
			double diff_vy = h_genParticles->at(2).daughter(0)->vy()-h_Kshorts->at(k).vy();
			double diff_vz = h_genParticles->at(2).daughter(0)->vz()-h_Kshorts->at(k).vz();

			double error_vx = h_Kshorts->at(k).vertexCovariance(0,0);
			double error_vy = h_Kshorts->at(k).vertexCovariance(1,1);
			
			histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vx"]->Fill(diff_vx);
			if(error_vx<0.01) histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vx_error_vx_smaller_0.01"]->Fill(diff_vx);
			histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vy"]->Fill(diff_vy);
			if(error_vy<0.01) histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vy_error_vy_smaller_0.01"]->Fill(diff_vy);
			histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vz"]->Fill(diff_vz);

			histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vxy"]->Fill(sqrt(diff_vx*diff_vx+diff_vy*diff_vy));
			if(error_vx<0.01 && error_vy<0.01)histos_th1f[b+"h_genPartices_Kshort_origV0_diff_vxy_error_vx_and_vy_smaller_0.01"]->Fill(sqrt(diff_vx*diff_vx+diff_vy*diff_vy));
			
			histos_th2f[b+"h_genPartices_Kshort_origV0_pt_diff_vx"]->Fill(h_genParticles->at(2).pt() , diff_vx);
			histos_th2f[b+"h_genPartices_Kshort_origV0_diff_vx_error_vx"]->Fill(abs(diff_vx), h_Kshorts->at(k).vertexCovariance(0,0));
			
		}	

	}
  }

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
  if(v_L0_phi.empty() && v_L0_eta.empty() && v_Ks_phi.empty() && v_Ks_eta.empty()) emptyVectors = true;//this means they have been cleared in the previous event or it is the first event

  
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


  if(!emptyVectors){
	v_L0_phi.clear();
	v_L0_eta.clear();
	v_Ks_phi.clear();
	v_Ks_eta.clear();
  } 


 //make some distributions for the gen particles
 //
 if(h_genParticles.isValid() ){
	for(unsigned int i = 0; i < h_genParticles->size(); ++i){
			int pdgId = h_genParticles->at(i).pdgId();
			histos_th1f[b+"h_genParticles_pdgid"]->Fill(pdgId) ;			
			//bool validXiDecayParticle = false; 
			//if(pdgID == 13324 || pdgID == 3122 ||pdgID == -311 ||pdgID == -211 ||pdgID == 2212 ||pdgID == 211 || pdgID == -211 )validXiDecayParticle = true;
			string particle_name ="";
			if(pdgId == 13324) particle_name = "Xi1820";
			else if(pdgId == 3122) particle_name = "Lambda";
			else if(pdgId == -310) particle_name = "Kshort";
			if(pdgId == 13324 || pdgId == 3122 ||pdgId == -310){
				histos_th1f[b+"h_genParticles_"+particle_name+"_pt_all"]->Fill(h_genParticles->at(i).pt());
				histos_th1f[b+"h_genParticles_"+particle_name+"_p_all"]->Fill(h_genParticles->at(i).p());
				histos_th1f[b+"h_genParticles_"+particle_name+"_energy_all"]->Fill(h_genParticles->at(i).energy());
				histos_th1f[b+"h_genParticles_"+particle_name+"_et_all"]->Fill(h_genParticles->at(i).et());
				histos_th1f[b+"h_genParticles_"+particle_name+"_mt_all"]->Fill(h_genParticles->at(i).mt());
				histos_th1f[b+"h_genParticles_"+particle_name+"_phi_all"]->Fill(h_genParticles->at(i).phi());
				histos_th1f[b+"h_genParticles_"+particle_name+"_theta_all"]->Fill(h_genParticles->at(i).theta());
				histos_th1f[b+"h_genParticles_"+particle_name+"_eta_all"]->Fill(h_genParticles->at(i).eta());
				histos_th1f[b+"h_genParticles_"+particle_name+"_vx_all"]->Fill(h_genParticles->at(i).daughter(0)->vx());
				histos_th1f[b+"h_genParticles_"+particle_name+"_vy_all"]->Fill(h_genParticles->at(i).daughter(0)->vy());
				histos_th1f[b+"h_genParticles_"+particle_name+"_vz_all"]->Fill(h_genParticles->at(i).daughter(0)->vz());
				histos_th1f[b+"h_genParticles_"+particle_name+"_vertexNormalizedChi2_all"]->Fill(h_genParticles->at(i).vertexNormalizedChi2());
			}
	}
 }
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!make the same distribution but for the S candidates (if they are there, they are only there if they ran through the 2nd filter)!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(h_sCands.isValid()) {
	unsigned int n_sCands = h_sCands->size();
	for (unsigned int i = 0; i < n_sCands; ++i) { //loop over reco Kshorts	
			//calculate the inv mass of the S: still have to subtract the neutron mass
			double px = h_sCands->at(i).px(); 
			double py = h_sCands->at(i).py();
			double pz = h_sCands->at(i).pz();
			TVector3 S_p(px,py,pz);
			double E = h_sCands->at(i).energy();
			    
			double neutronMass = 0.9395654133; //GeV
			TLorentzVector rCand_4momentum(px, py, pz, E);
			TLorentzVector neutron_4momentum(0,0,0,neutronMass);
			TLorentzVector sCand_4momentum = rCand_4momentum - neutron_4momentum; //h_sCand actually holds the resonance candidate --> need to substract neutron 4-momentum
			double M_Xi = h_sCands->at(i).mass();
			double M_S = sCand_4momentum.M();

			//angles
			double phi1 = h_sCands->at(i).daughter(0)->phi();
			double phi2 = h_sCands->at(i).daughter(1)->phi();
			double delta_phi = reco::deltaPhi(phi1, phi2);
			
			double eta1 = h_sCands->at(i).daughter(0)->eta();
			double eta2 = h_sCands->at(i).daughter(1)->eta();
			double delta_eta = eta1-eta2;
			
			double delta_R = sqrt(delta_phi*delta_phi+delta_eta*delta_eta);

			//vertices
			double vx_S_candidate = h_sCands->at(i).vx(); 
			double vy_S_candidate = h_sCands->at(i).vy(); 
			double vz_S_candidate = h_sCands->at(i).vz();
			TVector3 S_v(vx_S_candidate,vy_S_candidate,vz_S_candidate);

			double vx_S_candidate_var = h_sCands->at(i).vertexCovariance(0,0);
			double vy_S_candidate_var = h_sCands->at(i).vertexCovariance(1,1);
			double vz_S_candidate_var = h_sCands->at(i).vertexCovariance(2,2);

			//Lambda
			double vx_Lambda = h_sCands->at(i).daughter(0)->vx(); 
			double vy_Lambda = h_sCands->at(i).daughter(0)->vy(); 
			double vz_Lambda = h_sCands->at(i).daughter(0)->vz();
			TVector3 Lambda_v(vx_Lambda,vy_Lambda,vz_Lambda);

			double px_Lambda = h_sCands->at(i).daughter(0)->px(); 
			double py_Lambda = h_sCands->at(i).daughter(0)->py(); 
			double pz_Lambda = h_sCands->at(i).daughter(0)->pz();
			TVector3 Lambda_p(px_Lambda,py_Lambda,pz_Lambda);

			/*double vx_Lambda_var = h_sCands->at(i).daughter(0)->vertexCovariance(0,0);
			double vy_Lambda_var = h_sCands->at(i).daughter(0)->vertexCovariance(1,1);
			double vz_Lambda_var = h_sCands->at(i).daughter(0)->vertexCovariance(2,2);
			*/
			//Kshort
			double vx_Kshort = h_sCands->at(i).daughter(1)->vx(); 
			double vy_Kshort = h_sCands->at(i).daughter(1)->vy(); 
			double vz_Kshort = h_sCands->at(i).daughter(1)->vz();
			TVector3 Kshort_v(vx_Kshort,vy_Kshort,vz_Kshort);

			double px_Kshort = h_sCands->at(i).daughter(1)->px(); 
			double py_Kshort = h_sCands->at(i).daughter(1)->py(); 
			double pz_Kshort = h_sCands->at(i).daughter(1)->pz();
			TVector3 Kshort_p(px_Kshort,py_Kshort,pz_Kshort);

			/*double vx_Kshort_var = h_sCands->at(i).daughter(1)->vertexCovariance(0,0);
			double vy_Kshort_var = h_sCands->at(i).daughter(1)->vertexCovariance(1,1);
			double vz_Kshort_var = h_sCands->at(i).daughter(1)->vertexCovariance(2,2);
			*/
			//xy distance between beamspot and S vertex   		
			double lxy_S_b = lxy(S_v, beamspot);
			double lxy_S_b_std_dev = std_dev_lxy(vx_S_candidate, vy_S_candidate, vx_S_candidate_var, vy_S_candidate_var, bx_x, bx_y, pow(bx_x_std_dev,2) , pow(bx_y_std_dev,2));
			double lxy_S_b_signed =  lxy_signed(S_v, beamspot, S_p);
	
			//xy distance between beamspot and Ks decay
			double lxy_Kshort_b = lxy(Kshort_v, beamspot);
			//double lxy_Kshort_b_std_dev = std_dev_lxy(vx_Kshort, vy_Kshort, vx_Kshort_var, vy_Kshort_var, bx_x, bx_y, pow(bx_x_std_dev,2) , pow(bx_y_std_dev,2));

			//xy distance between beamspot and Lambda decay
			double lxy_Lambda_b = lxy(Lambda_v, beamspot);
			//double lxy_Lambda_b_std_dev = std_dev_lxy(vx_Lambda, vy_Lambda, vx_Lambda_var, vy_Lambda_var, bx_x, bx_y, pow(bx_x_std_dev,2) , pow(bx_y_std_dev,2));

			double significance_lxy_S_b = lxy_S_b/lxy_S_b_std_dev;


			double xy_PCA_S_b = xy_PCA_line_point(S_v, S_p, beamspot);
			double xy_signed_PCA_S_b = xy_signed_PCA_line_point(S_v, S_p, beamspot);
			double xy_PCA_Kshort_b = xy_PCA_line_point(Kshort_v, Kshort_p, beamspot);
			double xy_PCA_Lambda_b = xy_PCA_line_point(Lambda_v, Lambda_p, beamspot);
			
			double z_PCA_S_b = z_PCA_line_point(S_v, S_p, beamspot);
			double z_PCA_Kshort_b = z_PCA_line_point(Kshort_v, Kshort_p, beamspot);
			double z_PCA_Lambda_b = z_PCA_line_point(Lambda_v, Lambda_p, beamspot);
			

			//impact parameters
			//std::cout << "starting to get the dxy" << std::endl;
			//const reco::Candidate* S_candidate_daughter = h_sCands->at(i).daughter(0);
			//const reco::Track* S_daughter_track =S_candidate_daughter->bestTrack();	
			//const reco::Candidate* S_candidate = h_sCands->at(i);
			//std::cout << "got the mother" << std::endl;
			//const reco::Track* S_candidate_track =S_candidate->bestTrack();	
			//std::cout << "got the track" << std::endl;
			//std::cout << "S candidate transv impact parameter (dxy): " << S_candidate_track->dxy() << std::endl;
			//if(S_daughter_track)std::cout << "S candidate daughter pt : " << S_daughter_track->dxy() << std::endl;

			//angular correlations of the daughter:
			//for the S candidates
		        if(h_sCands->at(i).charge() == 1){	
				histos_th1f[b+"h_Sdaughters_L0_Ks_delta_phi"]->Fill(delta_phi);
				histos_th1f[b+"h_Sdaughters_L0_Ks_delta_eta"]->Fill(delta_eta);
				histos_th1f[b+"h_Sdaughters_L0_Ks_delta_R"]->Fill(delta_R);
				histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta"]->Fill(delta_phi,delta_eta);
				if( abs(delta_phi) > 0.1 || abs(delta_eta) > 0.1 )histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent"]->Fill(delta_phi,delta_eta);
				//make angular correlation plots for the particles which survive the constraint on the decay vertex
				if(lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01){
					
					histos_th1f[b+"h_Sdaughters_L0_Ks_delta_phi_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(delta_phi);
					histos_th1f[b+"h_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(delta_eta);
					histos_th1f[b+"h_Sdaughters_L0_Ks_delta_R_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(delta_R);
					histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(delta_phi,delta_eta);
				
				}
/*
				if(lxy_S_b > 0.5) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_0.5cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b > 1.0) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.0cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b > 1.5) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.5cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b > 1.9) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b_std_dev < 0.01) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.01cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b_std_dev < 0.02) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.02cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b_std_dev < 0.03) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.03cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b_std_dev < 0.04) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.04cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b_std_dev < 0.06) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.06cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b_std_dev < 0.08) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.08cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b_std_dev < 0.1) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_smaller_0.1cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b < 0.5) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_0.5cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b < 1.0) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.0cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b < 1.5) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.5cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b < 1.9) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.9cm"]->Fill(delta_phi,delta_eta);
                                if(lxy_S_b_std_dev > 0.01) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.01cm"]->Fill(delta_phi,delta_eta);
                                if(lxy_S_b_std_dev > 0.02) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.02cm"]->Fill(delta_phi,delta_eta);
                                if(lxy_S_b_std_dev > 0.03) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.03cm"]->Fill(delta_phi,delta_eta);
                                if(lxy_S_b_std_dev > 0.04) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.04cm"]->Fill(delta_phi,delta_eta);
                                if(lxy_S_b_std_dev > 0.06) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.06cm"]->Fill(delta_phi,delta_eta);
                                if(lxy_S_b_std_dev > 0.08) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.08cm"]->Fill(delta_phi,delta_eta);
                                if(lxy_S_b_std_dev > 0.1) histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_error_larger_0.1cm"]->Fill(delta_phi,delta_eta);
				if(lxy_S_b < 1.9 && lxy_S_b_std_dev > 0.01)histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_smaller_1.9cm_and_error_larger_0.01cm"]->Fill(delta_phi,delta_eta);
*/
				cout << "1" << endl;
				//look at the kinematics of the S candidates in the reconance
				if(fabs(delta_phi) < 0.1 && fabs(delta_eta) > 0.5 && fabs(delta_eta) < 3){
       		                        histos_th1f[b+"h_s_candidates_mass_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).mass());
					histos_th1f[b+"h_s_candidates_pt_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).pt());
					histos_th1f[b+"h_s_candidates_p_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).p());
					histos_th1f[b+"h_s_candidates_energy_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).energy());
					histos_th1f[b+"h_s_candidates_et_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).et());
					histos_th1f[b+"h_s_candidates_mt_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).mt());
					histos_th1f[b+"h_s_candidates_phi_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).phi());
					histos_th1f[b+"h_s_candidates_theta_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).theta());
					histos_th1f[b+"h_s_candidates_eta_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).eta());
					histos_th1f[b+"h_s_candidates_vx_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).vx());
					histos_th1f[b+"h_s_candidates_vy_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).vy());
					histos_th1f[b+"h_s_candidates_vz_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).vz());
					histos_th1f[b+"h_s_candidates_vertexNormalizedChi2_in_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).vertexNormalizedChi2());
					histos_th1f[b+"h_s_candidates_lxy_in_delta_phi_delta_eta_corr"]->Fill(lxy_S_b);
					histos_th1f[b+"h_s_candidates_lxy_signed_in_delta_phi_delta_eta_corr"]->Fill(lxy_S_b_signed);
					if(lxy_S_b_std_dev<0.1)histos_th1f[b+"h_s_candidates_lxy_signed_in_delta_phi_delta_eta_corr_error_lxy_smaller_0p1"]->Fill(lxy_S_b_signed);
					histos_th1f[b+"h_s_candidates_error_lxy_in_delta_phi_delta_eta_corr"]->Fill(lxy_S_b_std_dev);
				}
				//select the ones in the 0,0 peak
				else if(fabs(delta_phi) < 0.1 && fabs(delta_eta) < 0.1){
				cout << "2" << endl;
					histos_th1f[b+"h_s_candidates_mass_in_peak"]->Fill(h_sCands->at(i).mass());
					histos_th1f[b+"h_s_candidates_pt_in_peak"]->Fill(h_sCands->at(i).pt());
					histos_th1f[b+"h_s_candidates_phi_in_peak"]->Fill(h_sCands->at(i).phi());
					histos_th1f[b+"h_s_candidates_eta_in_peak"]->Fill(h_sCands->at(i).eta());
					histos_th1f[b+"h_s_candidates_vx_in_peak"]->Fill(h_sCands->at(i).vx());
					histos_th1f[b+"h_s_candidates_vy_in_peak"]->Fill(h_sCands->at(i).vy());
					histos_th1f[b+"h_s_candidates_vz_in_peak"]->Fill(h_sCands->at(i).vz());
					histos_th1f[b+"h_s_candidates_lxy_in_peak"]->Fill(lxy_S_b);
                                        histos_th1f[b+"h_s_candidates_lxy_signed_in_peak"]->Fill(lxy_S_b_signed);
                                        if(lxy_S_b_std_dev<0.1)histos_th1f[b+"h_s_candidates_lxy_signed_in_peak_error_lxy_smaller_0p1"]->Fill(lxy_S_b_signed);
                                        histos_th1f[b+"h_s_candidates_error_lxy_in_peak"]->Fill(lxy_S_b_std_dev);
				}
				else{
				cout << "3" << endl;
       		                        histos_th1f[b+"h_s_candidates_mass_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).mass());
					histos_th1f[b+"h_s_candidates_pt_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).pt());
					histos_th1f[b+"h_s_candidates_p_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).p());
					histos_th1f[b+"h_s_candidates_energy_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).energy());
					histos_th1f[b+"h_s_candidates_et_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).et());
					histos_th1f[b+"h_s_candidates_mt_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).mt());
					histos_th1f[b+"h_s_candidates_phi_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).phi());
					histos_th1f[b+"h_s_candidates_theta_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).theta());
					histos_th1f[b+"h_s_candidates_eta_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).eta());
					histos_th1f[b+"h_s_candidates_vx_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).vx());
					histos_th1f[b+"h_s_candidates_vy_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).vy());
					histos_th1f[b+"h_s_candidates_vz_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).vz());
					histos_th1f[b+"h_s_candidates_vertexNormalizedChi2_outside_delta_phi_delta_eta_corr"]->Fill(h_sCands->at(i).vertexNormalizedChi2());
					histos_th1f[b+"h_s_candidates_lxy_outside_delta_phi_delta_eta_corr"]->Fill(lxy_S_b);
					histos_th1f[b+"h_s_candidates_lxy_signed_outside_delta_phi_delta_eta_corr"]->Fill(lxy_S_b_signed);
					if(lxy_S_b_std_dev<0.1)histos_th1f[b+"h_s_candidates_lxy_signed_outside_delta_phi_delta_eta_corr_error_lxy_smaller_0p1"]->Fill(lxy_S_b_signed);
					histos_th1f[b+"h_s_candidates_error_lxy_outside_delta_phi_delta_eta_corr"]->Fill(lxy_S_b_std_dev);

				}

			}
			cout << "3.1" << endl;
			//for the anti-S candidates
			if(h_sCands->at(i).charge() == -1){

				histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_phi"]->Fill(delta_phi);
				histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_eta"]->Fill(delta_eta);
				histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_R"]->Fill(delta_R);
				histos_th2f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta"]->Fill(delta_phi,delta_eta);
				if( abs(delta_phi) > 0.1 || abs(delta_eta) > 0.1 )histos_th2f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent"]->Fill(delta_phi,delta_eta);
				//make angular correlation plots for the particles which survive the constraint on the decay vertex
				if(lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01){
					
					histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(delta_phi);
					histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(delta_eta);
					histos_th1f[b+"h_anti_Sdaughters_L0_Ks_delta_R_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(delta_R);
					histos_th2f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(delta_phi,delta_eta);
					if( abs(delta_phi) > 0.1 || abs(delta_eta) > 0.1 )histos_th2f[b+"h_anti_Sdaughters_L0_Ks_delta_phi_delta_eta_no_cent_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(delta_phi,delta_eta);
				
				}

			}
			cout << "3.2" << endl;
		        //distances to the beamspot
		        //for the S candidates
		        if(h_sCands->at(i).charge() == 1){	
				histos_th1f[b+"h_S_vtx_distance_to_beamspot"]->Fill(lxy_S_b);
				histos_th1f[b+"h_S_vtx_distance_to_beamspot_error"]->Fill(lxy_S_b_std_dev);
				//look for #PV dependence on the error of the reconstructed vertex
				 cout << "3.2.1" << endl;
				double n_PVs = h_nPVs->at(0);
				if(n_PVs <= 10) {
					histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_smaller_11"]->Fill(lxy_S_b_std_dev);
					histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_smaller_11"]->Fill(lxy_S_b_std_dev);
				}
				if(n_PVs > 10 && n_PVs <= 20){
				 	histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_10_smaller_or_equal_to_20"]->Fill(lxy_S_b_std_dev);
				 	histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_10_smaller_or_equal_to_20"]->Fill(lxy_S_b_std_dev);
				}
				if(n_PVs > 20 && n_PVs <= 30){ 
					histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_20_smaller_or_equal_to_30"]->Fill(lxy_S_b_std_dev);
					histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_20_smaller_or_equal_to_30"]->Fill(lxy_S_b_std_dev);
				}
				 cout << "3.2.2" << endl;
				if(n_PVs > 30 && n_PVs <= 40){
					histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_30_smaller_or_equal_to_40"]->Fill(lxy_S_b_std_dev);
					histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_30_smaller_or_equal_to_40"]->Fill(lxy_S_b_std_dev);
				}
				if(n_PVs > 40){
					 histos_th1f[b+"h_S_vtx_distance_to_beamspot_for_nPVs_larger_than_40"]->Fill(lxy_S_b_std_dev);
					 histos_th1f[b+"h_S_vtx_distance_to_beamspot_error_for_nPVs_larger_than_40"]->Fill(lxy_S_b_std_dev);
				}
				
				histos_th2f[b+"h_delta_phi_S_vtx_distance_to_beamspot"]->Fill(delta_phi, lxy_S_b);
				 cout << "3.2.2.1" << endl;
				//histos_th3f[b+"h_delta_phi_delta_eta_S_vtx_distance_to_beamspot"]->Fill(delta_phi, delta_eta, lxy_S_b);
				 cout << "3.2.3" << endl;
				
				if(lxy_S_b_std_dev < 0.01){
					histos_th2f[b+"h_delta_phi_S_vtx_distance_to_beamspot_vertex_beamspot_dist_error_smaller_0.01cm"]->Fill(delta_phi, lxy_S_b);
				//	histos_th3f[b+"h_delta_phi_delta_eta_S_vtx_distance_to_beamspot_vertex_beamspot_dist_error_smaller_0.01cm"]->Fill(delta_phi, delta_eta, lxy_S_b);

				}
				 cout << "3.2.4" << endl;
				histos_th2f[b+"h_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error"]->Fill(lxy_S_b, lxy_S_b_std_dev);
				if(1 < abs(delta_phi) && abs(delta_phi) < 2.5)histos_th2f[b+"h_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error_for_delta_phi_between_1_and_2.5"]->Fill(lxy_S_b, lxy_S_b_std_dev);
				histos_th2f[b+"h_delta_phi_S_vtx_distance_to_beamspot_significance"]->Fill(delta_phi, significance_lxy_S_b);
				histos_th2f[b+"h_S_vtx_distance_to_beamspot_vx_vy"]->Fill(vx_S_candidate,vy_S_candidate);
				 cout << "3.2.5" << endl;
				if(lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01){
					histos_th2f[b+"h_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(vx_S_candidate,vy_S_candidate);
				}
				if(lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01 && 1 < abs(delta_phi) && abs(delta_phi) < 2.5) histos_th2f[b+"h_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm_delta_phi_between_1_and_2.5"]->Fill(vx_S_candidate,vy_S_candidate);
				 cout << "3.2.6" << endl;
			}
			cout << "3.3" << endl;
		        //for the anti-S candidates
		        if(h_sCands->at(i).charge() == -1){	
				histos_th1f[b+"h_anti_S_vtx_distance_to_beamspot"]->Fill(lxy_S_b);
				histos_th1f[b+"h_anti_S_vtx_distance_to_beamspot_error"]->Fill(lxy_S_b_std_dev);
				histos_th2f[b+"h_delta_phi_anti_S_vtx_distance_to_beamspot"]->Fill(delta_phi, lxy_S_b);
				histos_th2f[b+"h_anti_S_vtx_distance_to_beamspot_S_vtx_distance_to_beamspot_error"]->Fill(lxy_S_b, lxy_S_b_std_dev);
				histos_th2f[b+"h_anti_S_vtx_distance_to_beamspot_vx_vy"]->Fill(vx_S_candidate,vy_S_candidate);
				if(lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01){
					histos_th2f[b+"h_anti_S_vtx_distance_to_beamspot_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(vx_S_candidate,vy_S_candidate);
				}
				if(lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01 && 1 < abs(delta_phi) && abs(delta_phi) < 2.5) histos_th2f[b+"h_anti_S_vx_vy_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm_delta_phi_between_1_and_2.5"]->Fill(vx_S_candidate,vy_S_candidate);
			}

			cout << "3.4" << endl;
			//masses of the S candidates
			if(h_sCands->at(i).charge() == 1){
				histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter"]->Fill(M_S);
				if(lxy_S_b > 1.9) histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm"]->Fill(M_S);
				if(lxy_S_b_std_dev < 0.01) histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_error_smaller_0.01cm"]->Fill(M_S);
				if(lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01) histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(M_S);
				if(abs(delta_phi)<2.5 && abs(delta_phi)>1) histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5"]->Fill(M_S);
				if(abs(delta_phi)<2.5 && abs(delta_phi)>1 && lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01){
					 //std::cout << "EventID which survive the abs(delta_phi)<2.5 && abs(delta_phi)>1 && lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01 cuts: " << iEvent.GetEventID() << std::endl;
					 histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(M_S);
					 histos_th1f[b+"h_s_candidates_vertex_beamspot_dist_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(lxy_S_b);
					 histos_th1f[b+"h_s_candidates_error_vertex_beamspot_dist_error_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(lxy_S_b_std_dev);
					
				}
				if(delta_R > 0.8) histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8"]->Fill(M_S);
				if(delta_R > 0.8 && lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01) histos_th1f[b+"h_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(M_S);
			}
			cout << "3.5" << endl;
		
			//masses of the anti S candidates
			if(h_sCands->at(i).charge() == -1){
				histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter"]->Fill(M_S);
				if(lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01) histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(M_S);
				if(abs(delta_phi)<2.5 && abs(delta_phi)>1) histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5"]->Fill(M_S);
				if(abs(delta_phi)<2.5 && abs(delta_phi)>1 && lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01) histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_phi_daughters_between_1_and_2.5_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(M_S);
				if(delta_R > 0.8) histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8"]->Fill(M_S);
				if(delta_R > 0.8 && lxy_S_b > 1.9 && lxy_S_b_std_dev < 0.01) histos_th1f[b+"h_anti_s_candidates_mass_after_LambdaKshortVertexFilter_delta_R_larger_0.8_and_S_vertex_beamspot_dist_larger_1.9cm_and_error_smaller_0.01cm"]->Fill(M_S);
			}

				cout << "5" << endl;
			//masses of the reconance candidate (Xi)
			histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter"]->Fill(M_Xi);
			if(lxy_S_b < 0.1 && lxy_S_b_std_dev < 0.01)histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm"]->Fill(M_Xi);
			if(delta_R < 0.8)histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8"]->Fill(M_Xi);
			if(lxy_S_b < 0.1 && lxy_S_b_std_dev < 0.01 && delta_R < 0.8)histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8"]->Fill(M_Xi);
		
			//look at the cut by cut for the Xi candidate and look at the same time at the distribution of the other variables:
			histos_th1f[b+"h_all_r_candidates_delta_R"]->Fill(delta_R);
			histos_th1f[b+"h_all_r_candidates_r_vertex_beamspot_distance"]->Fill(lxy_S_b);
			histos_th1f[b+"h_all_r_candidates_r_vertex_beamspot_distance_error"]->Fill(lxy_S_b_std_dev);
			
			if(delta_R < 0.8){
				histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_for_delta_R_slammer_0.8"]->Fill(lxy_S_b);
				histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_error_r_vertex_beamspot_distance_for_delta_R_slammer_0.8"]->Fill(lxy_S_b_std_dev);
			}
			if(lxy_S_b < 0.1){
				histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_error_r_vertex_beamspot_distance_for_r_vertex_beamspot_distance_smaller_0.1"]->Fill(lxy_S_b_std_dev);
				histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_for_r_vertex_beamspot_distance_smaller_0.1"]->Fill(delta_R);
			}
			if(lxy_S_b_std_dev < 0.01){
				histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_for_error_r_vertex_beamspot_distance_smaller_0.01"]->Fill(lxy_S_b);
				histos_th1f[b+"h_all_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_for_error_r_vertex_beamspot_distance_smaller_0.01"]->Fill(delta_R);
			}
			
			if(h_sCands->at(i).charge()==1){

				histos_th1f[b+"h_r_candidates_mass_after_LambdaKshortVertexFilter"]->Fill(M_Xi);
				if(lxy_S_b < 0.1 && lxy_S_b_std_dev < 0.01)histos_th1f[b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm"]->Fill(M_Xi);
				if(delta_R < 0.8)histos_th1f[b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8"]->Fill(M_Xi);
				if(lxy_S_b < 0.1 && lxy_S_b_std_dev < 0.01 && delta_R < 0.8)histos_th1f[b+"h_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8"]->Fill(M_Xi);
				if(z_PCA_S_b < 0.05 && delta_R < 0.8){
                                	histos_th1f[b+"h_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8"]->Fill(M_Xi);
                        	}
                        	if(z_PCA_S_b < 0.01 && delta_R < 0.8){
                                	histos_th1f[b+"h_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8"]->Fill(M_Xi);
                        	}	
			}	

			if(h_sCands->at(i).charge()==-1){

				histos_th1f[b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter"]->Fill(M_Xi);
				if(lxy_S_b < 0.1 && lxy_S_b_std_dev < 0.01)histos_th1f[b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm"]->Fill(M_Xi);
				if(delta_R < 0.8)histos_th1f[b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_delta_R_smaller_0.8"]->Fill(M_Xi);
				if(lxy_S_b < 0.1 && lxy_S_b_std_dev < 0.01 && delta_R < 0.8)histos_th1f[b+"h_anti_r_candidates_mass_after_LambdaKshortVertexFilter_r_vertex_beamspot_distance_smaller_0.1_and_error_smaller_0.01cm_delta_R_smaller_0.8"]->Fill(M_Xi);
				if(z_PCA_S_b < 0.05 && delta_R < 0.8){
                                	histos_th1f[b+"h_anti_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8"]->Fill(M_Xi);
                        	}
                        	if(z_PCA_S_b < 0.01 && delta_R < 0.8){
                                	histos_th1f[b+"h_anti_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8"]->Fill(M_Xi);
                        	}	
			
			}

			//PCA for the reconance candidate:
			histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter"]->Fill(xy_PCA_S_b);
			histos_th1f[b+"h_r_candidates_dxy_signed_after_LambdaKshortVertexFilter"]->Fill(xy_signed_PCA_S_b);
			histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter"]->Fill(z_PCA_S_b);
			histos_th2f[b+"h_r_candidates_dxy_dz_after_LambdaKshortVertexFilter"]->Fill(xy_PCA_S_b ,z_PCA_S_b);
			histos_th2f[b+"h_r_candidates_signed_dxy_dz_after_LambdaKshortVertexFilter"]->Fill(xy_signed_PCA_S_b ,z_PCA_S_b);
			


			double n_PVs = h_nPVs->at(0);
			if(n_PVs <= 10) {
				histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_smaller_11"]->Fill(xy_PCA_S_b);
				histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_smaller_11"]->Fill(z_PCA_S_b);
			}
			if(n_PVs > 10 && n_PVs <= 20){
				histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_10_smaller_or_equal_to_20"]->Fill(xy_PCA_S_b);
				histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_10_smaller_or_equal_to_20"]->Fill(z_PCA_S_b);
			}
			if(n_PVs > 20 && n_PVs <= 30){ 
				histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_20_smaller_or_equal_to_30"]->Fill(xy_PCA_S_b);
				histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_20_smaller_or_equal_to_30"]->Fill(z_PCA_S_b);
			}
			if(n_PVs > 30 && n_PVs <= 40){
				histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_30_smaller_or_equal_to_40"]->Fill(xy_PCA_S_b);
				histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_30_smaller_or_equal_to_40"]->Fill(z_PCA_S_b);
			}
			if(n_PVs > 40){
				 histos_th1f[b+"h_r_candidates_dxy_after_LambdaKshortVertexFilter_for_nPVs_larger_than_40"]->Fill(xy_PCA_S_b);
				 histos_th1f[b+"h_r_candidates_dz_after_LambdaKshortVertexFilter_for_nPVs_larger_than_40"]->Fill(z_PCA_S_b);
			}
			


			if(xy_PCA_S_b < 0.1){
				histos_th1f[b+"h_all_r_candidates_mass_PCA_xy_smaller_0p1"]->Fill(M_Xi);
			}
			if(z_PCA_S_b < 0.05){
				histos_th1f[b+"h_all_r_candidates_mass_PCA_z_smaller_0p05"]->Fill(M_Xi);
			}
			if(z_PCA_S_b < 0.05 && delta_R < 0.8){
				histos_th1f[b+"h_all_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8"]->Fill(M_Xi);
			}	
			if(z_PCA_S_b < 0.01 && delta_R < 0.8){
				histos_th1f[b+"h_all_r_candidates_mass_PCA_z_smaller_0p01_and_delta_R_smaller_0p8"]->Fill(M_Xi);
			}	
			if(fabs(z_PCA_S_b) < 0.05 && delta_R < 0.8 && fabs(xy_signed_PCA_S_b) < 0.5){
				histos_th1f[b+"h_all_r_candidates_mass_PCA_z_smaller_0p05_and_delta_R_smaller_0p8_PCA_xy_smaller_0p5"]->Fill(M_Xi);
			}	
			//masses of the daughters
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
			//histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vx_std_dev"]->Fill(sqrt(vx_Lambda_var));
			//histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_vy_std_dev"]->Fill(sqrt(vy_Lambda_var));
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_lxy"]->Fill(lxy_Lambda_b);
			//histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_lxy_std_dev"]->Fill(lxy_Lambda_b_std_dev);
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_dxy"]->Fill(xy_PCA_Lambda_b);
                        histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Lambda_dz"]->Fill(z_PCA_Lambda_b);
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
			//histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vx_std_dev"]->Fill(sqrt(vx_Kshort_var));
			//histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_vy_std_dev"]->Fill(sqrt(vy_Kshort_var));
			histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_lxy"]->Fill(lxy_Kshort_b);
			//histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_lxy_std_dev"]->Fill(lxy_Kshort_b_std_dev);
                        histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_dxy"]->Fill(xy_PCA_Kshort_b);
                        histos_th1f[b+"h_s_candidates_after_LambdaKshortVertexFilter_Kshort_dz"]->Fill(z_PCA_Kshort_b);

			if(h_genParticles.isValid()  )
			{
				histos_th1f[b+"h_genPartices_Lambda_diff_pdgId"]->Fill(h_genParticles->at(1).pdgId());
				histos_th1f[b+"h_genPartices_Lambda_diff_pt"]->Fill(h_genParticles->at(1).pt()-h_sCands->at(i).daughter(0)->pt());
				histos_th1f[b+"h_genPartices_Lambda_diff_mt"]->Fill(h_genParticles->at(1).mt()-h_sCands->at(i).daughter(0)->mt());
				histos_th1f[b+"h_genPartices_Lambda_diff_phi"]->Fill(h_genParticles->at(1).phi()-h_sCands->at(i).daughter(0)->phi());
				histos_th1f[b+"h_genPartices_Lambda_diff_eta"]->Fill(h_genParticles->at(1).eta()-h_sCands->at(i).daughter(0)->eta());
				histos_th1f[b+"h_genPartices_Lambda_diff_vx"]->Fill(h_genParticles->at(1).daughter(0)->vx()-h_sCands->at(i).daughter(0)->vx());
				histos_th1f[b+"h_genPartices_Lambda_diff_vy"]->Fill(h_genParticles->at(1).daughter(0)->vy()-h_sCands->at(i).daughter(0)->vy());
				histos_th1f[b+"h_genPartices_Lambda_diff_vz"]->Fill(h_genParticles->at(1).daughter(0)->vz()-h_sCands->at(i).daughter(0)->vz());
				
				histos_th1f[b+"h_genPartices_Kshort_diff_pdgId"]->Fill(h_genParticles->at(2).pdgId());
				histos_th1f[b+"h_genPartices_Kshort_diff_pt"]->Fill(h_genParticles->at(2).pt()-h_sCands->at(i).daughter(1)->pt());
				histos_th1f[b+"h_genPartices_Kshort_diff_mt"]->Fill(h_genParticles->at(2).mt()-h_sCands->at(i).daughter(1)->mt());
				histos_th1f[b+"h_genPartices_Kshort_diff_phi"]->Fill(h_genParticles->at(2).phi()-h_sCands->at(i).daughter(1)->phi());
				histos_th1f[b+"h_genPartices_Kshort_diff_eta"]->Fill(h_genParticles->at(2).eta()-h_sCands->at(i).daughter(1)->eta());
				histos_th1f[b+"h_genPartices_Kshort_diff_vx"]->Fill(h_genParticles->at(2).daughter(0)->vx()-h_sCands->at(i).daughter(1)->vx());
				histos_th1f[b+"h_genPartices_Kshort_diff_vy"]->Fill(h_genParticles->at(2).daughter(0)->vy()-h_sCands->at(i).daughter(1)->vy());
				histos_th1f[b+"h_genPartices_Kshort_diff_vz"]->Fill(h_genParticles->at(2).daughter(0)->vz()-h_sCands->at(i).daughter(1)->vz());
			}	
	}
  }//end of sCands present

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!make some distributions for the S and r candidates
  if((h_r_MassFilter).isValid()){
  	for (unsigned int r = 0; r < (*h_r_MassFilter).size(); ++r) {
		histos_th1f[b+"h_all_r_candidates_mass_after_massFilter"]->Fill((*h_r_MassFilter)[r].mass());
		if((*h_r_MassFilter)[r].charge() == 1) histos_th1f[b+"h_r_candidates_mass_after_massFilter"]->Fill((*h_r_MassFilter)[r].mass());
		if((*h_r_MassFilter)[r].charge() == -1) histos_th1f[b+"h_anti_r_candidates_mass_after_massFilter"]->Fill((*h_r_MassFilter)[r].mass());
	}
  }
 
  if((h_s_MassFilter).isValid()){
  	for (unsigned int s = 0; s < (*h_s_MassFilter).size(); ++s) {
		
		//mass plots of the S candidates. The charge is the charge of the proton in the decay chain, so a antiLambda would give an antiproton
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
			histos_th1f[b+"h_s_candidates_vertexNdof_after_massFilter"]->Fill((*h_s_MassFilter)[s].vertexNdof());
			
			
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

		if(h_genParticles.isValid() )
		{
			//the 0th gen particle is always the Xi particle which corresponds to the S here
			histos_th1f[b+"h_genPartices_Scands_pdgId"]->Fill(h_genParticles->at(0).pdgId());
			histos_th1f[b+"h_genPartices_Scands_mass"]->Fill(h_genParticles->at(0).mass());
			histos_th1f[b+"h_genPartices_Scands_diff_pt"]->Fill(h_genParticles->at(0).pt()-(*h_s_MassFilter)[s].pt());
			histos_th1f[b+"h_genPartices_Scands_diff_mt"]->Fill(h_genParticles->at(0).mt()-(*h_s_MassFilter)[s].mt());
			histos_th1f[b+"h_genPartices_Scands_diff_phi"]->Fill(h_genParticles->at(0).phi()-(*h_s_MassFilter)[s].phi());
			histos_th1f[b+"h_genPartices_Scands_diff_eta"]->Fill(h_genParticles->at(0).eta()-(*h_s_MassFilter)[s].eta());
			histos_th1f[b+"h_genPartices_Scands_diff_vx"]->Fill(h_genParticles->at(0).daughter(0)->vx()-(*h_s_MassFilter)[s].vx());
			histos_th1f[b+"h_genPartices_Scands_diff_vy"]->Fill(h_genParticles->at(0).daughter(0)->vy()-(*h_s_MassFilter)[s].vy());
			histos_th1f[b+"h_genPartices_Scands_diff_vz"]->Fill(h_genParticles->at(0).daughter(0)->vz()-(*h_s_MassFilter)[s].vz());
			histos_th2f[b+"h_genPartices_Scands_pt_diff_vx"]->Fill(h_genParticles->at(0).pt(),abs(h_genParticles->at(0).daughter(0)->vx()-(*h_s_MassFilter)[s].vx()));
			histos_th2f[b+"h_genPartices_Scands_pt_diff_vy"]->Fill(h_genParticles->at(0).pt(),abs(h_genParticles->at(0).daughter(0)->vy()-(*h_s_MassFilter)[s].vy()));
		}	

		
		//THE BELOW DOES NOT WORK: I CANNOT ACCESS THE DAUGHTERS OF THESE S CANDIDATES, I GET A product not found exception. ALTHOUGH THE DAUGHTERS ARE PUT INSIDE IN THE massfilter...
		//some plots for the S particle daughters: the first added daughter is the Lambda, the second is the Kshort (look at LambdaKshortVertexFilter.cc for this)
/*		if((*h_s_MassFilter)[s].charge()==1){
			const reco::Candidate *daughter0 = (*h_s_MassFilter)[s].daughter(0); 
//			histos_th1f[b+"h_s_candidates_mass_after_massFilter_Lambda_mass"]->Fill((*h_s_MassFilter)[s].daughter(0)->mass());
//			histos_th1f[b+"h_s_candidates_mass_after_massFilter_Kshort_mass"]->Fill((*h_s_MassFilter)[s].daughter(1)->mass());
		}
		if((*h_s_MassFilter)[s].charge()==-1){
			histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter_Lambda_mass"]->Fill((*h_s_MassFilter)[s].daughter(0)->mass());
			histos_th1f[b+"h_anti_s_candidates_mass_after_massFilter_Kshort_mass"]->Fill((*h_s_MassFilter)[s].daughter(1)->mass());
		}
*/
	}
  }


} //end of analyzer

double Analyzer_V0_angular_correlation::lxy(TVector3 v1, TVector3 v2){
	double x1 = v1.X();
	double x2 = v2.X();
	double y1 = v1.Y();
	double y2 = v2.Y();
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}


double Analyzer_V0_angular_correlation::lxy_signed(TVector3 particle_vertex, TVector3 beamspot, TVector3 particle_direction){
	double x1 = particle_vertex.X();
	double x2 = beamspot.X();
	double y1 = particle_vertex.Y();
	double y2 = beamspot.Y();
	double lxy_signed = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
	TVector3 displacement = particle_vertex - beamspot;
	
	if(displacement*particle_direction < 0) lxy_signed = -lxy_signed;
	return lxy_signed;
}

double Analyzer_V0_angular_correlation::std_dev_lxy(double vx, double vy, double vx_var, double vy_var, double bx_x, double bx_y, double bx_x_var, double bx_y_var){
	
	double lxy_std_dev_nominator = pow(vx-bx_x,2)*(vx_var+bx_x_var) + pow(vy-bx_y,2)*(vy_var+bx_y_var);
	double lxy_std_dev_denominator = pow(vx-bx_x,2) + pow(vy-bx_y,2);
	double lxy_b_std_dev = sqrt(lxy_std_dev_nominator/lxy_std_dev_denominator);
	return lxy_b_std_dev;

}

TVector3 Analyzer_V0_angular_correlation::PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point){
   double normalise = sqrt(Vector_along_line.X()*Vector_along_line.X()+Vector_along_line.Y()*Vector_along_line.Y()+Vector_along_line.Z()*Vector_along_line.Z());
   TVector3 n(Vector_along_line.X()/normalise,Vector_along_line.Y()/normalise,Vector_along_line.Z()/normalise);
   TVector3 a = Point_line;
   TVector3 p = Point;

   //see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line (Vector formulation)
   return (a-p)-((a-p)*n)*n;
}

double Analyzer_V0_angular_correlation::xy_PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point){
  TVector3 shortest_distance = PCA_line_point(Point_line,  Vector_along_line, Point);
  return sqrt(shortest_distance.X()*shortest_distance.X()+shortest_distance.Y()*shortest_distance.Y());
}


double Analyzer_V0_angular_correlation::xy_signed_PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point){
  TVector3 shortest_distance = PCA_line_point(Point_line,  Vector_along_line, Point);
  double xy_signed_PCA_line_point = sqrt(shortest_distance.X()*shortest_distance.X()+shortest_distance.Y()*shortest_distance.Y());
  
  TVector3 displacement = Point_line - Point; 
  if(displacement*Vector_along_line<0)xy_signed_PCA_line_point = -xy_signed_PCA_line_point;
  return xy_signed_PCA_line_point;
}

double Analyzer_V0_angular_correlation::z_PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point){
  TVector3 shortest_distance = PCA_line_point(Point_line,  Vector_along_line, Point);
  return shortest_distance.Z();
}

void Analyzer_V0_angular_correlation::endJob()
{
}

Analyzer_V0_angular_correlation::~Analyzer_V0_angular_correlation()
{
}


DEFINE_FWK_MODULE(Analyzer_V0_angular_correlation);
