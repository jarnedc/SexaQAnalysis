#include "SexaQAnalysis/Data_analysis/interface/Analyzer.h"
#include <typeinfo>


Analyzer::Analyzer(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  m_bsTag    (pset.getParameter<edm::InputTag>("beamspot")),
  m_vertexTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  //m_rCandsTag(pset.getParameter<edm::InputTag>("resonCandidates")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_KshortsTag(pset.getParameter<edm::InputTag>("lambdaCollection")),
  m_LambdasTag(pset.getParameter<edm::InputTag>("kshortCollection")),
  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_vertexToken(consumes<vector<reco::Vertex> >(m_vertexTag)),
  //m_rCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_rCandsTag)),
  //m_sCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_sCandsTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_KshortsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_KshortsTag)),
  m_LambdasToken(consumes<vector<reco::VertexCompositeCandidate> >(m_LambdasTag))


{

}


void Analyzer::beginJob() {
  
 
  
  histos_th1f[b+"nPV"]      = m_fs->make<TH1F>(b+"nPV",     a+" Number of PV; #PVs",60,0.,60);
  histos_th1f[b+"n_sCand"]  = m_fs->make<TH1F>(b+"n_sCand",     a+" Number of Kshort-Lambda vertex fits; #S candidates",8,0.,8);
  
  
  ///MASSES
  
  histos_th1f[b+"rCandMass_L0"] = m_fs->make<TH1F>(b+"rCandMass_L0",b+"rCandMass_L0; mass (GeV)",1000,0,10.);
  histos_th1f[b+"rCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors"]=m_fs->make<TH1F>(b+"rCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors",b+"rCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_L0"] = m_fs->make<TH1F>(b+"sCandMass_L0",b+"sCandMass_L0; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_L0_absDeltaPhi_between_1_and_2.5_cut"] = m_fs->make<TH1F>(b+"sCandMass_L0_absDeltaPhi_between_1_and_2.5_cut",b+"sCandMass_L0_absDeltaPhi_between_1_and_2.5_cut; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors"] = m_fs->make<TH1F>(b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors",b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut"] = m_fs->make<TH1F>(b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut",b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut; mass (GeV)",1000,0,10.);
    
  histos_th1f[b+"rCandMass_antiL0"] = m_fs->make<TH1F>(b+"rCandMass_antiL0",b+"rCandMass_antiL0; mass (GeV)",1000,0,10.);
  histos_th1f[b+"rCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors"]=m_fs->make<TH1F>(b+"rCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors",b+"rCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_antiL0"] = m_fs->make<TH1F>(b+"sCandMass_antiL0",b+"sCandMass_antiL0; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_antiL0_absDeltaPhi_between_1_and_2.5_cut"] = m_fs->make<TH1F>(b+"sCandMass_antiL0_absDeltaPhi_between_1_and_2.5_cut",b+"sCandMass_antiL0_absDeltaPhi_between_1_and_2.5_cut; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors"] = m_fs->make<TH1F>(b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors",b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut"] = m_fs->make<TH1F>(b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut",b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut; mass (GeV)",1000,0,10.);
  
  //with fermi motion of neutron
  histos_th1f[b+"rCandMass_L0_Fermi"] = m_fs->make<TH1F>(b+"rCandMass_L0_Fermi",b+"rCandMass_L0_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"rCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi"]=m_fs->make<TH1F>(b+"rCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi",b+"rCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_L0_Fermi"] = m_fs->make<TH1F>(b+"sCandMass_L0_Fermi",b+"sCandMass_L0_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_L0_absDeltaPhi_between_1_and_2.5_cut_Fermi"] = m_fs->make<TH1F>(b+"sCandMass_L0_absDeltaPhi_between_1_and_2.5_cut_Fermi",b+"sCandMass_L0_absDeltaPhi_between_1_and_2.5_cut_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi"] = m_fs->make<TH1F>(b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi",b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut_Fermi"] = m_fs->make<TH1F>(b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut_Fermi",b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut_Fermi; mass (GeV)",1000,0,10.);
    
  histos_th1f[b+"rCandMass_antiL0_Fermi"] = m_fs->make<TH1F>(b+"rCandMass_antiL0_Fermi",b+"rCandMass_antiL0_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"rCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi"]=m_fs->make<TH1F>(b+"rCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi",b+"rCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_antiL0_Fermi"] = m_fs->make<TH1F>(b+"sCandMass_antiL0_Fermi",b+"sCandMass_antiL0_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_antiL0_absDeltaPhi_between_1_and_2.5_cut_Fermi"] = m_fs->make<TH1F>(b+"sCandMass_antiL0_absDeltaPhi_between_1_and_2.5_cut_Fermi",b+"sCandMass_antiL0_absDeltaPhi_between_1_and_2.5_cut_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi"] = m_fs->make<TH1F>(b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi",b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi; mass (GeV)",1000,0,10.);
  histos_th1f[b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut_Fermi"] = m_fs->make<TH1F>(b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut_Fermi",b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut_Fermi; mass (GeV)",1000,0,10.);
  
 
  ///rCAND DISTIBUTIONS
  //rCAND mass distributions with certain cuts
  histos_th1f[b+"rCandMass_L0_with_dxy_smaller_0.2_and_dR_daughters_smaller_0.8"] = m_fs->make<TH1F>(b+"rCandMass_L0_with_dxy_smaller_0.2_and_dR_daughters_smaller_0.8",b+"rCandMass_L0_with_dxy_smaller_0.2_and_dR_daughters_smaller_0.8; mass (GeV)",1000,0,7);
  histos_th1f[b+"rCandMass_antiL0_with_dxy_smaller_0.2_and_dR_daughters_smaller_0.8"] = m_fs->make<TH1F>(b+"rCandMass_antiL0_with_dxy_smaller_0.2_and_dR_daughters_smaller_0.8",b+"rCandMass_antiL0_with_dxy_smaller_0.2_and_dR_daughters_smaller_0.8; mass (GeV)",1000,0,7);

  //rCAND mass distributions with certain cuts
  histos_th1f[b+"rCandMass_with_dxy_smaller_0p2_and_dr_daughters_smaller_0p8"] = m_fs->make<TH1F>(b+"rCandMass_with_dxy_smaller_0p2_and_dr_daughters_smaller_0p8",b+"rCandMass_with_dxy_smaller_0p2_and_dr_daughters_smaller_0p8; mass (GeV)",1000,-10.,10.);
  histos_th1f[b+"rCandMass_with_dxy_smaller_0p1_and_dr_daughters_smaller_0p8"] = m_fs->make<TH1F>(b+"rCandMass_with_dxy_smaller_0p1_and_dr_daughters_smaller_0p8",b+"rCandMass_with_dxy_smaller_0p1_and_dr_daughters_smaller_0p8; mass (GeV)",1000,-10.,10.);
  histos_th1f[b+"rCandMass_with_dxy_smaller_0p05_and_dr_daughters_smaller_0p8"] = m_fs->make<TH1F>(b+"rCandMass_with_dxy_smaller_0p05_and_dr_daughters_smaller_0p8",b+"rCandMass_with_dxy_smaller_0p05_and_dr_daughters_smaller_0p8; mass (GeV)",1000,-10.,10.);
  histos_th1f[b+"rCandMass_with_dxy_smaller_0p01_and_dr_daughters_smaller_0p8"] = m_fs->make<TH1F>(b+"rCandMass_with_dxy_smaller_0p01_and_dr_daughters_smaller_0p8",b+"rCandMass_with_dxy_smaller_0p01_and_dr_daughters_smaller_0p8; mass (GeV)",1000,-10.,10.);
  //histos_th1f[b+"rCand_mass_below_2GeV_Eta"] = m_fs->make<TH1F>(b+"rCand_mass_below_2GeV_Eta",b+"rCand_mass_below_2GeV_Eta; pseudorapidity",1000,-5.,5.);
  
  ///sCAND DISTIBUTIONS
  
  histos_th1f[b+"sCand_delta_phi"] = m_fs->make<TH1F>(b+"sCand_delta_phi",b+"sCand_delta_phi; delta Phi",1000,-4.,4.);
  histos_th1f[b+"sCand_delta_phi_dz(PCA_PV0)_below_1mm"] = m_fs->make<TH1F>(b+"sCand_delta_phi_dz(PCA_PV0)_below_1mm",b+"sCand_delta_phi_dz(PCA_PV0)_below_1mm; delta Phi",1000,-4.,4.);
  histos_th1f[b+"sCand_delta_phi_dz(PCA_PV0)_above_1mm"] = m_fs->make<TH1F>(b+"sCand_delta_phi_dz(PCA_PV0)_above_1mm",b+"sCand_delta_phi_dz(PCA_PV0)_above_1mm; delta Phi",1000,-4.,4.);
    
  histos_th1f[b+"sCand_delta_eta_dz(PCA_PV0)_below_1mm"] = m_fs->make<TH1F>(b+"sCand_delta_eta_dz(PCA_PV0)_below_1mm",b+"sCand_delta_eta_dz(PCA_PV0)_below_1mm; delta eta",1000,-5.,5.);
  histos_th1f[b+"sCand_delta_eta_dz(PCA_PV0)_above_1mm"] = m_fs->make<TH1F>(b+"sCand_delta_eta_dz(PCA_PV0)_above_1mm",b+"sCand_delta_eta_dz(PCA_PV0)_above_1mm; delta eta",1000,-5.,5.);
  histos_th1f[b+"sCand_delta_eta"] = m_fs->make<TH1F>(b+"sCand_delta_eta",b+"sCand_delta_eta; delta eta",1000,-5.,5.);
  
  histos_th1f[b+"sCand_delta_R"] = m_fs->make<TH1F>(b+"sCand_delta_R",b+"sCand_delta_R; delta R",1000,0,6);
  
  histos_th1f[b+"sCand_delta_R_dz(PCA_PV0)_below_1mm"] = m_fs->make<TH1F>(b+"sCand_delta_R_dz(PCA_PV0)_below_1mm",b+"sCand_delta_R_dz(PCA_PV0)_below_1mm; delta R",1000,0,6);
  histos_th1f[b+"sCand_delta_R_dz(PCA_PV0)_above_1mm"] = m_fs->make<TH1F>(b+"sCand_delta_R_dz(PCA_PV0)_above_1mm",b+"sCand_delta_R_dz(PCA_PV0)_above_1mm; delta R",1000,0,6);

  
  histos_th1f[b+"sCand_dxy(sCandVtx_bs)"]= m_fs->make<TH1F>(b+"sCand_dxy(sCandVtx_bs)",b+"sCand_dxy(sCandVtx_bs); distance (cm)",1000,0,2.5);
  
  histos_th1f[b+"sCand_edxy(sCandVtx_bs)"]= m_fs->make<TH1F>(b+"sCand_edxy(sCandVtx_bs)",b+"sCand_edxy(sCandVtx_bs); standard deviation on distance (cm)",1000,0,1);
  histos_th2f[b+"sCand_dxy(sCandVtx_bs)_edxy(sCandVtx_bs)"]= m_fs->make<TH2F>(b+"sCand_dxy(sCandVtx_bs)_edxy(sCandVtx_bs)",b+"sCand_dxy(sCandVtx_bs)_edxy(sCandVtx_bs); distance (cm); standard deviation on distance (cm)",1000,0,4, 1000, 0, 2);
  
  histos_th1f[b+"sCand_edxy(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut"]= m_fs->make<TH1F>(b+"sCand_edxy(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut",b+"sCand_edxy(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut; standard deviation on distance (cm)",1000,0,1);
  histos_th2f[b+"sCand_dxy(sCandVtx_bs)_edxy(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut"]= m_fs->make<TH2F>(b+"sCand_dxy(sCandVtx_bs)_edxy(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut",b+"sCand_dxy(sCandVtx_bs)_edxy(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut; distance (cm); standard deviation on distance (cm)",1000,0,4, 1000, 0, 2);

  histos_th1f[b+"sCand_edz(sCandVtx_bs)"]= m_fs->make<TH1F>(b+"sCand_edz(sCandVtx_bs)",b+"sCand_edz(sCandVtx_bs); standard deviation on distance (cm)",1000,0,1);
  histos_th1f[b+"sCand_edz(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut"]= m_fs->make<TH1F>(b+"sCand_edz(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut",b+"sCand_edz(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut; standard deviation on distance (cm)",1000,0,1);


  histos_th1f[b+"sCand_dxyz(sCandVtx_bs)"]= m_fs->make<TH1F>(b+"sCand_dxyz(sCandVtx_bs)",b+"sCand_dxyz(sCandVtx_bs); distance (cm)",1000,0,30);
  histos_th1f[b+"sCand_dxy_signif_(sCandVtx_bs)"]= m_fs->make<TH1F>(b+"sCand_dxy_signif_(sCandVtx_bs)",b+"sCand_dxy_signif_(sCandVtx_bs); significance",1000,0,5);
  histos_th1f[b+"sCand_dxyz_signif_(sCandVtx_bs)"]= m_fs->make<TH1F>(b+"sCand_dxyz_signif_(sCandVtx_bs)",b+"sCand_dxyz_signif_(sCandVtx_bs); significance",1000,0,30);
  
  
  ///sCAND DAUGHTER DISTRIBUTIONS
  
  histos_th1f[b+"sCand_daughter_L0_pt"] = m_fs->make<TH1F>(b+"sCand_daughter_L0_pt",b+"sCand_daughter_L0_pt; pt (GeV)",1000,0,10.);
  histos_th1f[b+"sCand_daughter_Ks_pt"] = m_fs->make<TH1F>(b+"sCand_daughter_Ks_pt",b+"sCand_daughter_Ks_pt; pt (GeV)",1000,0,10.);
  
  histos_th1f[b+"sCand_daughter_L0_eta"] = m_fs->make<TH1F>(b+"sCand_daughter_L0_eta",b+"sCand_daughter_L0_eta; pseudorapidity",1000,-5,5.);
  histos_th1f[b+"sCand_daughter_Ks_eta"] = m_fs->make<TH1F>(b+"sCand_daughter_Ks_eta",b+"sCand_daughter_Ks_eta; pseudorapidity",1000,-5,5.);

    
  histos_th1f[b+"sCand_daughters_m"] = m_fs->make<TH1F>(b+"sCand_daughters_mass",b+"sCand_daughters_mass; mass (GeV)",1000,0,10.);
  
  //SCATTERPLOTS
  
  histos_th2f[b+"scatterplot_sCand_xy_pos_wrt_origin"]    = m_fs->make<TH2F>(b+"scatterplot_sCand_xy_pos_wrt_origin", b+"scatterplot_sCand_xy_pos_wrt_origin; x position (cm); y position (cm)",1000,-100.,100.,1000,-100.,100.);
  histos_th2f[b+"scatterplot_sCand_rz_pos_wrt_origin"]    = m_fs->make<TH2F>(b+"scatterplot_sCand_rz_pos_wrt_origin", b+"scatterplot_sCand_rz_pos_wrt_origin;z position (cm); r=dxy position (cm) ",1000,-100,100.,1000,0,100.);

  histos_th2f[b+"scatterplot_sCand_xy_pos_wrt_bs"]    = m_fs->make<TH2F>(b+"scatterplot_sCand_xy_pos_wrt_bs", b+"scatterplot_sCand_xy_pos_wrt_bs; x position (cm); y position (cm)",1000,-100.,100.,1000,-100.,100.);
  histos_th2f[b+"scatterplot_sCand_rz_pos_wrt_bs"]    = m_fs->make<TH2F>(b+"scatterplot_sCand_rz_pos_wrt_bs", b+"scatterplot_sCand_rz_pos_wrt_bs;z position (cm); r=dxy position (cm)",1000,-100,100.,1000,0,100.);

  histos_th2f[b+"scatterplot_PV0_xy_pos"]    = m_fs->make<TH2F>(b+"scatterplot_PV0_xy_pos", b+"scatterplot_PV0_xy_pos; x position (cm); y position (cm)",1000,-3.,3.,1000,-3.,3.);

  histos_th2f[b+"scatterplot_bs_xy_pos"]    = m_fs->make<TH2F>(b+"scatterplot_bs_xy_pos", b+"scatterplot_bs_xy_pos; x position (cm); y position (cm)",1000,-0.5,0.5,1000,-0.5,0.5);
  histos_th2f[b+"scatterplot_sCand_PCA_to_bs_xy_pos_wrt_bs"]    = m_fs->make<TH2F>(b+"scatterplot_sCand_PCA_to_bs_xy_pos_wrt_bs", b+"scatterplot_sCand_PCA_to_bs_xy_pos_wrt_bs; x position (cm); y position (cm)",1000,-10.,10.,1000,-10.,10.);
  
  ///DELTA PHI INVESTIGATION
  
  histos_th2f[b+"sCand_delta_phi_dxy(PCA_beamspot)_signed"]= m_fs->make<TH2F>(b+"sCand_delta_phi_dxy(PCA_beamspot)_signed",b+"sCand_delta_phi_dxy(PCA_beamspot)_signed; delta Phi; distance (cm)",1000,-4.,4., 1000, 0, 20);
  histos_th2f[b+"sCand_L0_delta_phi_dxy(sCandVtx_beamspot)_signed"]= m_fs->make<TH2F>(b+"sCand_L0_delta_phi_dxy(sCandVtx_beamspot)_signed",b+"sCand_L0_delta_phi_dxy(sCandVtx_beamspot)_signed; delta Phi; distance (cm)",1000,-4.,4., 1000, 0, 30);
  histos_th2f[b+"sCand_antiL0_delta_phi_dxy(sCandVtx_beamspot)_signed"]= m_fs->make<TH2F>(b+"sCand_antiL0_delta_phi_dxy(sCandVtx_beamspot)_signed",b+"sCand_antiL0_delta_phi_dxy(sCandVtx_beamspot)_signed; delta Phi; distance (cm)",1000,-4.,4., 1000, 0, 30);
  histos_th2f[b+"sCand_L0_delta_phi_dxy(sCandVtx_beamspot)_signed_zoom"]= m_fs->make<TH2F>(b+"sCand_L0_delta_phi_dxy(sCandVtx_beamspot)_signed_zoom",b+"sCand_L0_delta_phi_dxy(sCandVtx_beamspot)_signed_zoom; delta Phi; distance (cm)",1000,-4.,4., 1000, 0, 5);
  histos_th2f[b+"sCand_antiL0_delta_phi_dxy(sCandVtx_beamspot)_signed_zoom"]= m_fs->make<TH2F>(b+"sCand_antiL0_delta_phi_dxy(sCandVtx_beamspot)_signed_zoom",b+"sCand_antiL0_delta_phi_dxy(sCandVtx_beamspot)_signed_zoom; delta Phi; distance (cm)",1000,-4.,4., 1000, 0, 5);

  ///Vxy PLOTS
  
  histos_th1f[b+"sCand_dxy(sCandVtx_beamspot)_signed_L0"] = m_fs->make<TH1F>(b+"sCand_dxy(sCandVtx_beamspot)_signed_L0",b+"sCand_dxy(sCandVtx_beamspot)_signed_L0; distance (cm)",1000,-2.5,2.5);
  histos_th1f[b+"sCand_dxy(sCandVtx_beamspot)_signed_anti_L0"] = m_fs->make<TH1F>(b+"sCand_dxy(sCandVtx_beamspot)_signed_anti_L0",b+"sCand_dxy(sCandVtx_beamspot)_signed_anti_L0; distance (cm)",1000,-2.5,2.5);

  histos_th1f[b+"sCand_dxy(sCandVtx_beamspot)_signed_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy(sCandVtx_beamspot)_signed_L0_deltaPhiCut",b+"sCand_dxy(sCandVtx_beamspot)_signed_L0_deltaPhiCut; distance (cm)",1000,-2.5,2.5);
  histos_th1f[b+"sCand_dxy(sCandVtx_beamspot)_signed_anti_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy(sCandVtx_beamspot)_signed_anti_L0_deltaPhiCut",b+"sCand_dxy(sCandVtx_beamspot)_signed_anti_L0_deltaPhiCut; distance (cm)",1000,-2.5,2.5);

  ///dxy PCA to beamspot 
  
  histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0"] = m_fs->make<TH1F>(b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0",b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0; distance (cm)",1000,-0.5,0.5);
  histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0"] = m_fs->make<TH1F>(b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0",b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0; distance (cm)",1000,-0.5,0.5);
  histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0"] = m_fs->make<TH1F>(b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0",b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0; distance (cm)",1000,-0.5,0.5);
  histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0"] = m_fs->make<TH1F>(b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0",b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0; distance (cm)",1000,-0.5,0.5);

  histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut",b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut; distance (cm)",1000,-0.5,0.5);
  histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut",b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut; distance (cm)",1000,-0.5,0.5);
  histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut",b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut; distance (cm)",1000,-0.5,0.5);
  histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut",b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut; distance (cm)",1000,-0.5,0.5);
  
  ///dz PCA to PV(0)
  
  histos_th1f[b+"sCand_dz(PCA_PV0)"] = m_fs->make<TH1F>(b+"sCand_dz(PCA_PV0)",b+"sCand_dz(PCA_PV0); PCA to PV(0) z-distance (cm)",1000,-0.5,0.5);

  
  ///dxy PCA significance PLOTS 
  histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0"] = m_fs->make<TH1F>(b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0",b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0; significance (dxy / sigma_dxy)", 1000, -20, 20);
  histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0"] = m_fs->make<TH1F>(b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0",b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0; significance (dxy / sigma_dxy)", 1000, -20, 20);
  histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0"] = m_fs->make<TH1F>(b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0",b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0; significance (dxy / sigma_dxy)",1000, -20, 20);
  histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0"] = m_fs->make<TH1F>(b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0",b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0; significance (dxy / sigma_dxy)", 1000, -20, 20);

  histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut",b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut; significance (dxy / sigma_dxy)", 1000, -20, 20);
  histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut",b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut; significance (dxy / sigma_dxy)", 1000, -20, 20);
  histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut",b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut; significance (dxy / sigma_dxy)", 1000, -20, 20);
  histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH1F>(b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut",b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut; significance (dxy / sigma_dxy)", 1000, -20, 20);

  histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_L0"] = m_fs->make<TH2F>(b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_L0",b+"sCand_dxy(PCA_beamspot)_2D_Vxy_between_5mm_18mm_L0; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 20);
  histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_anti_L0"] = m_fs->make<TH2F>(b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_anti_L0",b+"sCand_dxy(PCA_beamspot)_2D_Vxy_between_5mm_18mm_anti_L0; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 20);
  histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_L0"] = m_fs->make<TH2F>(b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_L0",b+"sCand_dxy(PCA_beamspot)_2D_Vxy_over_18mm_L0; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 20);
  histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_anti_L0"] = m_fs->make<TH2F>(b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_anti_L0",b+"sCand_dxy(PCA_beamspot)_2D_Vxy_over_18mm_anti_L0; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 20);

  histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_L0_deltaPhiCut"] = m_fs->make<TH2F>(b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_L0_deltaPhiCut",b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_L0_deltaPhiCut; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 20);
  histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH2F>(b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut",b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 20);
  histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_L0_deltaPhiCut"] = m_fs->make<TH2F>(b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_L0_deltaPhiCut",b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_L0_deltaPhiCut; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 20);
  histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH2F>(b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_anti_L0_deltaPhiCut",b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_anti_L0_deltaPhiCut; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 20);



  ///RECO LAMBDA PLOTS
  
  histos_th1f[b+"L0_eta"] = m_fs->make<TH1F>(b+"L0_eta",b+"L0_eta; pseudorapidity",1000,-4,4.);
  histos_th1f[b+"L0_phi"] = m_fs->make<TH1F>(b+"L0_phi",b+"L0_phi; phi",1000,-5,5.);
  histos_th1f[b+"L0_pt"] = m_fs->make<TH1F>(b+"L0_pt",b+"L0_pt; pt (GeV)",1000,0,10.);
    
  histos_th1f[b+"L0_dxy(PCA_beamspot)"] = m_fs->make<TH1F>(b+"L0_dxy(PCA_beamspot)",b+"L0_dxy(PCA_beamspot); distance (cm)", 1000, 0, 0.5);

  histos_th1f[b+"L0_dxy_signif_(PCA_beamspot)"]= m_fs->make<TH1F>(b+"L0_dxy_signif_(PCA_beamspot)",b+"L0_dxy_signif_(PCA_beamspot); significance (dxy / sigma_dxy)", 1000, 0, 20);
    
  ///RECO KSHORT PLOTS
  
  histos_th1f[b+"Ks_eta"] = m_fs->make<TH1F>(b+"Ks_eta",b+"Ks_eta; pseudorapidity",1000,-4,4.);
  histos_th1f[b+"Ks_phi"] = m_fs->make<TH1F>(b+"Ks_phi",b+"Ks_phi; phi",1000,-5,5.);
  histos_th1f[b+"Ks_pt"] = m_fs->make<TH1F>(b+"Ks_pt",b+"Ks_pt; pt (GeV)",1000,0,10.);  
    
  histos_th1f[b+"Ks_dxy(PCA_beamspot)"] = m_fs->make<TH1F>(b+"Ks_dxy(PCA_beamspot)",b+"Ks_dxy(PCA_beamspot); distance (cm)", 1000, 0, 0.5);

  histos_th1f[b+"Ks_dxy_signif_(PCA_beamspot)"]= m_fs->make<TH1F>(b+"Ks_dxy_signif_(PCA_beamspot)",b+"Ks_dxy_signif_(PCA_beamspot); significance (dxy / sigma_dxy)", 1000, 0, 20);
    

}


void Analyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

  edm::Handle<vector<reco::Vertex> > h_vert; //Primary vertices
  iEvent.getByToken(m_vertexToken, h_vert);


  // resonance candidates
  
  //edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_rCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
  //iEvent.getByToken(m_rCandsToken, h_rCands);

  // sexaquark candidates
  
  //edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_sCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
  //iEvent.getByToken(m_sCandsToken, h_sCands);
  
  //lambdaKshortVertexFilter sexaquark candidates
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);
  
  //reco Kshorts V0
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_Kshorts;
  iEvent.getByToken(m_KshortsToken, h_Kshorts);
  
  //reco Lambdas V0
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_Lambdas;
  iEvent.getByToken(m_LambdasToken, h_Lambdas);

  // Check validity
  if(!h_bs.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_bsTag << " ... skip entry !" << endl;
    return;
  }

  if(!h_vert.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_vertexTag << " ... skip entry !" << endl;
    return;
  }
/*
  if(!h_rCands.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_rCandsTag << " ... skip entry !" << endl;
    return;
  }
*/

  if(!h_sCands.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_sCandsTag << " ... skip entry !" << endl;
    return;
  }
  
  if(!h_Lambdas.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_LambdasTag << " ... skip entry !" << endl;
    return;
  }
  
  if(!h_Kshorts.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_KshortsTag << " ... skip entry !" << endl;
    return;
  }

  // GLOBAL EVENT INFORMATIONS //
  m_nRun   = iEvent.id().run();
  m_nLumi  = iEvent.luminosityBlock();
  m_nEvent = iEvent.id().event();
  Int_t nPV = h_vert->size();
  Int_t n_sCands = h_sCands->size();
  
  unsigned int n_Lambdas = h_Lambdas->size(); //Jarne: this is the original V0 collection
  unsigned int n_Kshorts = h_Kshorts->size(); //Jarne: this is the original V0 collection
  //Int_t n_rCands = h_rCands->size();
 

  //std::cout<<m_nRun<<"\t"<<m_nLumi<<"\t"<<m_nEvent<<std::endl;

  histos_th1f[b+"nPV"]->Fill(nPV);
  histos_th1f[b+"n_sCand"]->Fill(h_sCands->size());    
  //histos_th1f["n_rCand"]->Fill(h_rCands->size());
  if (h_vert->size() == 0) return;  //only look at events with reconstructed L0-Ks vertices //Jarne: why exactly do you do this?
  
  
  //beamspot position
  float bs_x = h_bs->x0(); 
  float bs_y = h_bs->y0();
  float bs_z = h_bs->z0();
  
  
  TVector3 bs_pos(bs_x,bs_y,bs_z); 
  TVector2 bs_pos2D(bs_x,bs_y);    
  histos_th2f[b+"scatterplot_bs_xy_pos"]->Fill(bs_x,bs_y);
  
  //PV(0) position

  TVector3 PV0_pos(h_vert->at(0).x(),h_vert->at(0).y(),h_vert->at(0).z()); 
  TVector2 PV0_pos2D(h_vert->at(0).x(),h_vert->at(0).y()); 
  
  histos_th2f[b+"scatterplot_PV0_xy_pos"]->Fill(h_vert->at(0).x(), h_vert->at(0).y());

//@@@@@ rCAND = sCAND with rCAND 4 momentum = sum of sCAND daughter's 4 momenta !!!!! 


  for (unsigned int i = 0; i < n_Lambdas; ++i) { //loop over reco V0 Lambdas
    
  histos_th1f[b+"L0_eta"]->Fill(h_Lambdas->at(i).eta());
  histos_th1f[b+"L0_phi"]->Fill(h_Lambdas->at(i).phi());
  histos_th1f[b+"L0_pt"]->Fill(h_Lambdas->at(i).pt());

  
  
    
  TVector2 momentum(h_Lambdas->at(i).px(), h_Lambdas->at(i).py());
  TVector2 position(h_Lambdas->at(i).vx(), h_Lambdas->at(i).vy());
    
  //Extrapolate vertex momentum in direction of the unit momentum vector to the Point of Closest Approach with the beamspot in 2D (x,y)

  float alpha = ((momentum * bs_pos2D) - (momentum * position)) / (momentum * momentum);

  TVector2 PCA2D = position + (alpha * momentum);

  TVector2 PCA_bs2D = PCA2D - bs_pos2D;

  float dxy_PCA_bs = sqrt(pow(PCA_bs2D.X(),2)+pow(PCA_bs2D.Y(),2));

  TVector2 ptcl_bs2D = position - bs_pos2D; //vector between particle (Kshort or Lambda) and beamspot
    
  //float dxy_ptcl_bs = sqrt(pow(ptcl_bs2D.X(),2)+pow(ptcl_bs2D.Y(),2));
    
  float dxy_PCA_bs_signed = dxy_PCA_bs * signum(ptcl_bs2D * momentum); //2D dot product
  
  histos_th1f[b+"L0_dxy(PCA_beamspot)"] ->Fill(dxy_PCA_bs_signed);
  
  TMatrixD CovMx2D(2,2); //2D xy covariance matrix
        CovMx2D(0,0)=h_Lambdas->at(i).vertexCovariance(0,0); //elements are Sigma(i,j)²
        CovMx2D(0,1)=h_Lambdas->at(i).vertexCovariance(0,1);
        CovMx2D(1,0)=h_Lambdas->at(i).vertexCovariance(1,0);
        CovMx2D(1,1)=h_Lambdas->at(i).vertexCovariance(1,1);
    
        TVectorD EValues;
        TMatrixD EVectors = CovMx2D.EigenVectors(EValues);
    
  for (Int_t i = 0; i < EValues.GetNrows(); ++i) { 
    TVectorD EVector(TMatrixTColumn_const<double>(EVectors, i));  
  }

  float MinEValue = ( EValues(0)<EValues(1) ) ? EValues(0) : EValues(1);
  //float dxy_signif_PCA_bs= abs(dxy_PCA_bs_signed)/sqrt(abs(MinEValue)); //significance = dxy / sigma_dxy ≃ dxy / sqrt(min(eigenvalues(particle vtx 2D Cov matrix)))
    float dxy_signif_PCA_bs_signed = dxy_PCA_bs_signed/sqrt(abs(MinEValue)); //significance = dxy / sigma_dxy ≃ dxy / sqrt(min(eigenvalues(particle vtx 2D Cov matrix)))
    
    histos_th1f[b+"L0_dxy_signif_(PCA_beamspot)"] ->Fill(dxy_signif_PCA_bs_signed);
   
    
    
  }
  
  for (unsigned int i = 0; i < n_Kshorts; ++i) { //loop over reco Kshorts
    
  histos_th1f[b+"Ks_eta"]->Fill(h_Kshorts->at(i).eta());
  histos_th1f[b+"Ks_phi"]->Fill(h_Kshorts->at(i).phi());
  histos_th1f[b+"Ks_pt"]->Fill(h_Kshorts->at(i).pt());

  
  
    
  TVector2 momentum(h_Kshorts->at(i).px(), h_Kshorts->at(i).py());
  TVector2 position(h_Kshorts->at(i).vx(), h_Kshorts->at(i).vy());
    
  //Extrapolate vertex momentum in direction of the unit momentum vector to the Point of Closest Approach with the beamspot in 2D (x,y)

  float alpha = ((momentum * bs_pos2D) - (momentum * position)) / (momentum * momentum);

  TVector2 PCA2D = position + (alpha * momentum);

  TVector2 PCA_bs2D = PCA2D - bs_pos2D;

  float dxy_PCA_bs = sqrt(pow(PCA_bs2D.X(),2)+pow(PCA_bs2D.Y(),2));

  TVector2 ptcl_bs2D = position - bs_pos2D; //vector between particle (Kshort or Lambda) and beamspot
    
  //float dxy_ptcl_bs = sqrt(pow(ptcl_bs2D.X(),2)+pow(ptcl_bs2D.Y(),2));
    
  float dxy_PCA_bs_signed = dxy_PCA_bs * signum(ptcl_bs2D * momentum); //2D dot product
  
  histos_th1f[b+"Ks_dxy(PCA_beamspot)"] ->Fill(dxy_PCA_bs_signed);
  
  TMatrixD CovMx2D(2,2); //2D xy covariance matrix
        CovMx2D(0,0)=h_Kshorts->at(i).vertexCovariance(0,0); //elements are Sigma(i,j)²
        CovMx2D(0,1)=h_Kshorts->at(i).vertexCovariance(0,1);
        CovMx2D(1,0)=h_Kshorts->at(i).vertexCovariance(1,0);
        CovMx2D(1,1)=h_Kshorts->at(i).vertexCovariance(1,1);
    
        TVectorD EValues;
        TMatrixD EVectors = CovMx2D.EigenVectors(EValues);
    
  for (Int_t i = 0; i < EValues.GetNrows(); ++i) { 
    TVectorD EVector(TMatrixTColumn_const<double>(EVectors, i));  
  }

  float MinEValue = ( EValues(0)<EValues(1) ) ? EValues(0) : EValues(1);
  //float dxy_signif_PCA_bs= abs(dxy_PCA_bs_signed)/sqrt(abs(MinEValue)); //significance = dxy / sigma_dxy ≃ dxy / sqrt(min(eigenvalues(particle vtx 2D Cov matrix)))
    float dxy_signif_PCA_bs_signed = dxy_PCA_bs_signed/sqrt(abs(MinEValue)); //significance = dxy / sigma_dxy ≃ dxy / sqrt(min(eigenvalues(particle vtx 2D Cov matrix)))
    
    histos_th1f[b+"Ks_dxy_signif_(PCA_beamspot)"] ->Fill(dxy_signif_PCA_bs_signed);
    
    
  }
  
  

  for (unsigned int i = 0; i < h_sCands->size(); ++i) { //loop over S candidates
    
    int proton_charge = h_sCands->at(i).charge();   //in the case the sCand charge (which is the proton charge) is (-1)+1, the sCand daughter includes an (anti)Lambda

    
    //sCand daughter distributions
   
    histos_th1f[b+"sCand_daughter_L0_pt"]->Fill(h_sCands->at(i).daughter(0)->pt()); 
    histos_th1f[b+"sCand_daughter_Ks_pt"]->Fill(h_sCands->at(i).daughter(1)->pt()); 
  
    histos_th1f[b+"sCand_daughter_L0_eta"]->Fill(h_sCands->at(i).daughter(0)->eta());
    histos_th1f[b+"sCand_daughter_Ks_eta"]->Fill(h_sCands->at(i).daughter(1)->eta());

    histos_th1f[b+"sCand_daughters_m"] ->Fill(h_sCands->at(i).daughter(0)->mass());   
    histos_th1f[b+"sCand_daughters_m"] ->Fill(h_sCands->at(i).daughter(1)->mass());   
    
    
    
    //cout <<h_sCands->at(i).vertexChi2() <<endl;

    //sCand vertex position and errors
    float x  = h_sCands->at(i).vx();
    float y  = h_sCands->at(i).vy();
    float z  = h_sCands->at(i).vz();
    float xe = h_sCands->at(i).vertexCovariance(0,0); //Watch out: these are variances, not standard deviations
    float ye = h_sCands->at(i).vertexCovariance(1,1);
    float ze = h_sCands->at(i).vertexCovariance(2,2);
    float edxy = sqrt(xe+ye); //Jarne: should this not be: sqrt(xe*xe+ye*ye); Florian: no because xe and ye are variances, not standard deviations

    
    TVector2 cand_pos2D(x,y);
    TVector3 cand_pos(x,y,z);
    
    float delta_phi = reco::deltaPhi(h_sCands->at(i).daughter(0)->phi(), h_sCands->at(i).daughter(1)->phi());

    float delta_eta = h_sCands->at(i).daughter(0)->eta() - h_sCands->at(i).daughter(1)->eta();

    float delta_R = deltaR(h_sCands->at(i).daughter(0)->eta(), h_sCands->at(i).daughter(0)->phi(), h_sCands->at(i).daughter(1)->eta(), h_sCands->at(i).daughter(1)->phi());
  
    
    histos_th1f[b+"sCand_dxy(sCandVtx_bs)"]->Fill(sqrt(pow(x - bs_x, 2) + pow(y - bs_y, 2)));
    histos_th1f[b+"sCand_dxyz(sCandVtx_bs)"]->Fill(sqrt(pow(x - bs_x, 2) + pow(y - bs_y, 2) + pow(z - bs_z, 2)) ); 
    histos_th1f[b+"sCand_dxy_signif_(sCandVtx_bs)"]->Fill(sqrt(pow(x - bs_x, 2) + pow(y - bs_y, 2))/sqrt(xe+ye)); 
    histos_th1f[b+"sCand_dxyz_signif_(sCandVtx_bs)"]->Fill(sqrt(pow(x - bs_x, 2) + pow(y - bs_y, 2) + pow(z - bs_z, 2))/sqrt(xe+ye+ze) );
  
  
	
    histos_th1f[b+"sCand_edxy(sCandVtx_bs)"]->Fill(edxy); //Here I assume the error on the beamspot position to be negligible compared to the error on the sCandVtx position
    histos_th1f[b+"sCand_edz(sCandVtx_bs)"]->Fill(sqrt(ze)); //ze = variance --> sqrt(ze) = std
    histos_th2f[b+"sCand_dxy(sCandVtx_bs)_edxy(sCandVtx_bs)"]->Fill(sqrt(pow(x - bs_x, 2) + pow(y - bs_y, 2)), edxy);
    
    if(1<abs(delta_phi) && abs(delta_phi)<2.5){
		histos_th1f[b+"sCand_edxy(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut"]->Fill(edxy); //Here I assume the error on the beamspot position to be negligible compared to the error on the sCandVtx position
		histos_th2f[b+"sCand_dxy(sCandVtx_bs)_edxy(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut"]->Fill(sqrt(pow(x - bs_x, 2) + pow(y - bs_y, 2)), edxy);
	    histos_th1f[b+"sCand_edz(sCandVtx_bs)_absDeltaPhi_between_1_and_2.5_cut"]->Fill(sqrt(ze)); //ze = variance --> sqrt(ze) = std
	}
    
    //momentum
    float px = h_sCands->at(i).px(); 
    float py = h_sCands->at(i).py();
    float pz = h_sCands->at(i).pz();
    TVector2 cand_mom2D(px,py);
    TVector3 cand_mom(px,py,pz);
    
    float E = h_sCands->at(i).energy();
    float neutronMass = 0.9395654133; //GeV
    TLorentzVector rCand_4momentum(px, py, pz, E);
    TLorentzVector neutron_4momentum(0,0,0,neutronMass);
    TLorentzVector sCand_4momentum = rCand_4momentum - neutron_4momentum; //h_sCand actually holds the resonance candidate --> need to substract neutron 4-momentum
    float rCand_mass = h_sCands->at(i).mass();
    //cout <<"mass ok: "<< (h_sCands->at(i).mass() - rCand_4momentum.M() )<<endl;
    float sCand_mass = sCand_4momentum.M();

     
	//Fermi motion: https://fias.uni-frankfurt.de/~svogel/lecture_ws_2011_12/slides_bratkovskaya_2.pdf
    //fermi motion momentum of nucleon approximately 250MeV -->Generate random momentum of neutron between 0 and 250MeV

    float FermiPx = RandomFloat(-0.250/sqrt(3), 0.250/sqrt(3)); //Fermi momentum uniformly distributed between 0 and 250/sqrt(3) MeV --> total Fermi momentum maximally 250MeV
    float FermiPy = RandomFloat(-0.250/sqrt(3), 0.250/sqrt(3));
    float FermiPz = RandomFloat(-0.250/sqrt(3), 0.250/sqrt(3));
    
    TLorentzVector neutron_4momentum_Fermi(FermiPx,FermiPy,FermiPz,neutronMass);
    TLorentzVector sCand_4momentum_Fermi = rCand_4momentum - neutron_4momentum_Fermi;
    float sCand_mass_Fermi = sCand_4momentum_Fermi.M();
    
      
    
    TMatrixD CovMx2D(2,2); //2D xy covariance matrix
    CovMx2D(0,0)=h_sCands->at(i).vertexCovariance(0,0); //elements are Sigma(i,j)²
    CovMx2D(0,1)=h_sCands->at(i).vertexCovariance(0,1);
    CovMx2D(1,0)=h_sCands->at(i).vertexCovariance(1,0);
    CovMx2D(1,1)=h_sCands->at(i).vertexCovariance(1,1);
    
    //////
    //CovMx2D.Print();
    
    
     
    //calculate eigenvalues and eigenvectors
    //used for significance in PCA plots
    //http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
    TVectorD EValues;
    TMatrixD EVectors = CovMx2D.EigenVectors(EValues);
    
    
  for (Int_t i = 0; i < EValues.GetNrows(); ++i) { 
    
    TVectorD EVector(TMatrixTColumn_const<double>(EVectors, i)); 
    //cout << "eigen-value " << i << " is " << EValues(i);
    //cout << " with eigen-vector ";  EVector.Print(); 
    
  }
  //cout<<endl;
  //}
  
  
  float MinEValue = ( EValues(0)<EValues(1) ) ? EValues(0) : EValues(1);
    
 
	
  if(proton_charge==1) histos_th1f[b+"rCandMass_L0"]->Fill(rCand_mass);
  if(proton_charge==-1) histos_th1f[b+"rCandMass_antiL0"]->Fill(rCand_mass);

  float dxy = sqrt(pow(x,2)+pow(y,2));
  float dxy_bs = sqrt(pow(x-bs_x,2)+pow(y-bs_y,2)); //xy distance to the beamspot

  if (dxy_bs > 2-3*edxy && edxy < .1 && sqrt(ze) < .1){
    if(proton_charge==1) histos_th1f[b+"rCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors"]->Fill(rCand_mass); 
    if(proton_charge==-1) histos_th1f[b+"rCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors"]->Fill(rCand_mass); 

  }
  

  //Extrapolate vertex momentum in direction of the unit momentum vector to the Point of Closest Approach with the beamspot in 2D (x,y)

  float alpha = ((cand_mom2D * bs_pos2D) - (cand_mom2D * cand_pos2D)) / (cand_mom2D * cand_mom2D);

  TVector2 PCA2D = cand_pos2D + (alpha * cand_mom2D);

  TVector2 PCA_bs2D = PCA2D - bs_pos2D;

  float dxy_PCA_bs = sqrt(pow(PCA_bs2D.X(),2)+pow(PCA_bs2D.Y(),2));

  TVector2 cand_bs2D = cand_pos2D - bs_pos2D; //vector between sCand and beamspot
    
  float dxy_cand_bs = sqrt(pow(cand_bs2D.X(),2)+pow(cand_bs2D.Y(),2));
    

  float dxy_PCA_bs_signed = dxy_PCA_bs * signum(cand_bs2D * cand_mom2D); //2D dot product

    //cout << signum(cand_bs * cand_mom) << endl;
    


    //Extrapolate vertex momentum in direction of the unit momentum vector to the Point of Closest Approach with the beamspot in 3D
    //Here the PCA is first calculated in 3D and then projected onto the 2D (x,y) plane, resulting in a 2D xy-distance
    //Not used in thesis
    
    /*
    TVector3 target = bs_pos; //Calculate PCA wrt to this position
    float alpha = ((cand_mom * target) - (cand_mom * cand_pos)) / (cand_mom * cand_mom); //This is calculated in 3D
    TVector3 PCA = cand_pos + (alpha * cand_mom);
    TVector3 PCA_target = PCA - target;
    float dxy_PCA_bs = sqrt(pow(PCA_target.X(),2)+pow(PCA_target.Y(),2));
    TVector3 cand_bs = cand_pos - target; //vector between sCand and beamspot
    
    float dxy_cand_bs_3D = sqrt(pow(cand_bs.X(),2)+pow(cand_bs.Y(),2));
    
    float dxy_PCA_bs_3D_signed = dxy_PCA_bs_3D * signum(cand_bs.XYvector() * cand_mom.XYvector()); //2D dot product
    */


  //Extrapolate vertex momentum in direction of the unit momentum vector to the Point of Closest Approach with the PV(0) in z direction
  //PCA first calculated in 3D and then projected to 1D z-distance
    
  TVector3 target = PV0_pos; //Calculate PCA wrt to this position

  alpha = ((cand_mom * target) - (cand_mom * cand_pos)) / (cand_mom * cand_mom); //This is calculated in 3D

  TVector3 PCA = cand_pos + (alpha * cand_mom);

  TVector3 PCA_target = PCA - target;

  float dz_PCA_PV0 = PCA_target.Z();
    
  histos_th1f[b+"sCand_dz(PCA_PV0)"]->Fill(dz_PCA_PV0);
    
  histos_th2f[b+"scatterplot_sCand_xy_pos_wrt_bs"] ->Fill(x-bs_x,y-bs_y);
  histos_th2f[b+"scatterplot_sCand_rz_pos_wrt_bs"] ->Fill(z-bs_z,sqrt(pow(x-bs_x,2)+pow(y-bs_y,2)) );
    
  histos_th2f[b+"scatterplot_sCand_xy_pos_wrt_origin"] ->Fill(x,y);
  histos_th2f[b+"scatterplot_sCand_rz_pos_wrt_origin"] ->Fill(z, sqrt(pow(x,2)+pow(y,2)) );

  histos_th2f[b+"scatterplot_sCand_PCA_to_bs_xy_pos_wrt_bs"]   ->Fill(PCA2D.X()-bs_x, PCA2D.Y()-bs_y);
  
 
  histos_th1f[b+"sCand_delta_R"]->Fill(delta_R);
  if(abs(dz_PCA_PV0)>0.1) histos_th1f[b+"sCand_delta_R_dz(PCA_PV0)_above_1mm"]->Fill(delta_R);
  if(abs(dz_PCA_PV0)<0.1) histos_th1f[b+"sCand_delta_R_dz(PCA_PV0)_below_1mm"]->Fill(delta_R);
  histos_th1f[b+"sCand_delta_phi"]->Fill(delta_phi);
  if(abs(dz_PCA_PV0)>0.1) histos_th1f[b+"sCand_delta_phi_dz(PCA_PV0)_above_1mm"]->Fill(delta_phi);
  if(abs(dz_PCA_PV0)<0.1) histos_th1f[b+"sCand_delta_phi_dz(PCA_PV0)_below_1mm"]->Fill(delta_phi);
  histos_th1f[b+"sCand_delta_eta"]->Fill(delta_eta);
  if(abs(dz_PCA_PV0)>0.1) histos_th1f[b+"sCand_delta_eta_dz(PCA_PV0)_above_1mm"]->Fill(delta_eta);
  if(abs(dz_PCA_PV0)<0.1) histos_th1f[b+"sCand_delta_eta_dz(PCA_PV0)_below_1mm"]->Fill(delta_eta);
  
  
  
  histos_th2f[b+"sCand_delta_phi_dxy(PCA_beamspot)_signed"]->Fill(delta_phi, dxy_PCA_bs_signed);
  if(proton_charge == 1) histos_th2f[b+"sCand_L0_delta_phi_dxy(sCandVtx_beamspot)_signed"]->Fill(delta_phi, dxy_cand_bs * signum(cand_bs2D * cand_mom2D));
  if(proton_charge == -1)   histos_th2f[b+"sCand_antiL0_delta_phi_dxy(sCandVtx_beamspot)_signed"]->Fill(delta_phi, dxy_cand_bs * signum(cand_bs2D * cand_mom2D));
  if(proton_charge == 1) histos_th2f[b+"sCand_L0_delta_phi_dxy(sCandVtx_beamspot)_signed_zoom"]->Fill(delta_phi, dxy_cand_bs * signum(cand_bs2D * cand_mom2D));
  if(proton_charge == -1)   histos_th2f[b+"sCand_antiL0_delta_phi_dxy(sCandVtx_beamspot)_signed_zoom"]->Fill(delta_phi, dxy_cand_bs * signum(cand_bs2D * cand_mom2D));

  ///sCAND MASS PLOTS WITH AND WITHOUT CUTS
  //Fermi motion: https://fias.uni-frankfurt.de/~svogel/lecture_ws_2011_12/slides_bratkovskaya_2.pdf
  //fermi motion momentum of nucleon approximately 250MeV -->Generate random momentum of neutron between 0 and 250MeV
  
  if(proton_charge==1) histos_th1f[b+"sCandMass_L0"]->Fill(sCand_mass);
  if(proton_charge==-1) histos_th1f[b+"sCandMass_antiL0"]->Fill(sCand_mass);
  if(proton_charge==1) histos_th1f[b+"sCandMass_L0_Fermi"]->Fill(sCand_mass_Fermi);
  if(proton_charge==-1) histos_th1f[b+"sCandMass_antiL0_Fermi"]->Fill(sCand_mass_Fermi);
  
  if(1<abs(delta_phi) && abs(delta_phi)<2.5){
    if(proton_charge==1) histos_th1f[b+"sCandMass_L0_absDeltaPhi_between_1_and_2.5_cut"]->Fill(sCand_mass);
    if(proton_charge==-1) histos_th1f[b+"sCandMass_antiL0_absDeltaPhi_between_1_and_2.5_cut"]->Fill(sCand_mass);
    if(proton_charge==1) histos_th1f[b+"sCandMass_L0_absDeltaPhi_between_1_and_2.5_cut_Fermi"]->Fill(sCand_mass_Fermi);
    if(proton_charge==-1) histos_th1f[b+"sCandMass_antiL0_absDeltaPhi_between_1_and_2.5_cut_Fermi"]->Fill(sCand_mass_Fermi);
  }
	

	
	
	
	
     
    //make rCandidate mass plot with cuts on the PCA (should be coming from the beam spot) and dR < a given value, from the simulation we saw that for the rCandidate, this could be the Xi1820 for example the has dR ~< 0.8
    if(fabs(dxy_cand_bs) < 0.2 && delta_R < 0.8){
	histos_th1f[b+"rCandMass_with_dxy_smaller_0p2_and_dr_daughters_smaller_0p8"]->Fill(h_sCands->at(i).mass());
    }     
    
     if(fabs(dxy_cand_bs) < 0.1 && delta_R < 0.8){
	histos_th1f[b+"rCandMass_with_dxy_smaller_0p1_and_dr_daughters_smaller_0p8"]->Fill(h_sCands->at(i).mass());
    }      
     if(fabs(dxy_cand_bs) < 0.05 && delta_R < 0.8){
	histos_th1f[b+"rCandMass_with_dxy_smaller_0p05_and_dr_daughters_smaller_0p8"]->Fill(h_sCands->at(i).mass());
    }  
     if(fabs(dxy_cand_bs) < 0.01 && delta_R < 0.8){
	histos_th1f[b+"rCandMass_with_dxy_smaller_0p01_and_dr_daughters_smaller_0p8"]->Fill(h_sCands->at(i).mass());
    } 
     //PCA AND SIGNIFICANCE PLOTS
    
  
  
  if (dxy_bs > 2-3*edxy && edxy < .1 && sqrt(ze) < .1){
    if(proton_charge==1) histos_th1f[b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors"]->Fill(sCand_mass); 
    if(proton_charge==-1) histos_th1f[b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors"]->Fill(sCand_mass); 
    if(proton_charge==1) histos_th1f[b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi"]->Fill(sCand_mass_Fermi); 
    if(proton_charge==-1) histos_th1f[b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_Fermi"]->Fill(sCand_mass_Fermi); 

  }
  
  if (dxy_bs > 2-3*edxy && edxy < .1 && sqrt(ze) < .1 && 1<abs(delta_phi) && abs(delta_phi)<2.5){
    if(proton_charge==1) histos_th1f[b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut"]->Fill(sCand_mass); 
    if(proton_charge==-1) histos_th1f[b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut"]->Fill(sCand_mass); 
	if(proton_charge==1) histos_th1f[b+"sCandMass_L0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut_Fermi"]->Fill(sCand_mass_Fermi); 
    if(proton_charge==-1) histos_th1f[b+"sCandMass_antiL0_with_dxy(sCandVtx_bs)_over_2cm_cut_and_conditions_on_errors_and_absDeltaPhi_between_1_and_2.5_cut_Fermi"]->Fill(sCand_mass_Fermi); 
  
  }
     
  //make rCandidate mass plot with cuts on the PCA (should be coming from the beam spot) and dR < a given value, from the simulation we saw that for the rCandidate, this could be the Xi1820 for example the has dR ~< 0.8
  if(fabs(dxy_cand_bs) < 0.2 && delta_R < 0.8){
    if(proton_charge==1) histos_th1f[b+"rCandMass_L0_with_dxy_smaller_0.2_and_dR_daughters_smaller_0.8"]->Fill(h_sCands->at(i).mass());
    if(proton_charge==-1) histos_th1f[b+"rCandMass_antiL0_with_dxy_smaller_0.2_and_dR_daughters_smaller_0.8"]->Fill(h_sCands->at(i).mass());

  }     
    
         
    
    ///PCA AND SIGNIFICANCE PLOTS
    
    float dxy_signif_PCA_bs= abs(dxy_PCA_bs_signed)/sqrt(abs(MinEValue)); //significance = dxy / sigma_dxy ≃ dxy / sqrt(min(eigenvalues(sCand vtx 2D Cov matrix)))
    float dxy_signif_PCA_bs_signed = dxy_PCA_bs_signed/sqrt(abs(MinEValue)); //significance = dxy / sigma_dxy ≃ dxy / sqrt(min(eigenvalues(sCand vtx 2D Cov matrix)))
   
    
   
    
    
    
    
    //in the case the sCand charge (which is the proton charge) is (-1)+1, the sCand daughter includes an (anti)Lambda
    
    
    
    if(0.5<dxy_cand_bs && dxy_cand_bs<1.8) {
    
      if(proton_charge == 1) histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0"] ->Fill(dxy_PCA_bs_signed);
      
     //In these 2D plots the distance is signed but the significance not (because it would waste space since only the 1st and 3 quadrants would be filled since both distance and significance would have the same sign)
      if(proton_charge == 1) histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_L0"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs); 
      if(proton_charge == -1) histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_anti_L0"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
      
      //Here the significance is signed because it is useful
      if(proton_charge == 1) histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0"] ->Fill(dxy_signif_PCA_bs_signed); 
      if(proton_charge == -1) histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0"] ->Fill(dxy_signif_PCA_bs_signed);
    
    } else if(1.8 < dxy_cand_bs) {
    
      if(proton_charge == 1) histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0"] ->Fill(dxy_PCA_bs_signed);
      
      
      if(proton_charge == 1) histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_L0"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
      if(proton_charge == -1) histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_anti_L0"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
   
      if(proton_charge == 1) histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0"] ->Fill(dxy_signif_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0"] ->Fill(dxy_signif_PCA_bs_signed);

    }
    
    if(proton_charge == 1) histos_th1f[b+"sCand_dxy(sCandVtx_beamspot)_signed_L0"]->Fill(dxy_cand_bs * signum(cand_bs2D * cand_mom2D));
    if(proton_charge == -1) histos_th1f[b+"sCand_dxy(sCandVtx_beamspot)_signed_anti_L0"]->Fill(dxy_cand_bs * signum(cand_bs2D * cand_mom2D));
    
    
    //Same plots but with a delta phi cut
    if(abs(delta_phi)<0.5 || abs(delta_phi) > (TMath::Pi() - 0.5)) continue;
    
  if(0.5<dxy_cand_bs && dxy_cand_bs<1.8) {
    
      if(proton_charge == 1) histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed);
      
     
      if(proton_charge == 1) histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
      if(proton_charge == -1) histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
      
      if(proton_charge == 1) histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut"] ->Fill(dxy_signif_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_signif_PCA_bs_signed);
    
    } else if(1.8 < dxy_cand_bs) {
    
      if(proton_charge == 1) histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f[b+"sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed);
      
      
      if(proton_charge == 1) histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
      if(proton_charge == -1) histos_th2f[b+"sCand_dxy_signif_(PCA_beamspot)_2D_Vxy_over_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
      
      if(proton_charge == 1) histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut"] ->Fill(dxy_signif_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f[b+"sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_signif_PCA_bs_signed);
    
    
    }
    
    if(proton_charge == 1) histos_th1f[b+"sCand_dxy(sCandVtx_beamspot)_signed_L0_deltaPhiCut"]->Fill(dxy_cand_bs * signum(cand_bs2D * cand_mom2D));
    if(proton_charge == -1) histos_th1f[b+"sCand_dxy(sCandVtx_beamspot)_signed_anti_L0_deltaPhiCut"]->Fill(dxy_cand_bs * signum(cand_bs2D * cand_mom2D));

  }
 



}


void Analyzer::endJob()
{
}

Analyzer::~Analyzer()
{
}


DEFINE_FWK_MODULE(Analyzer);
