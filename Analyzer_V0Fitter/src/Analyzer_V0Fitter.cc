#include "../interface/Analyzer_V0Fitter.h"
#include <typeinfo>

//#include "V0Fitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <typeinfo>
#include <memory>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

// pdg mass constants
namespace {
   const double piMass = 0.13957018;
   const double piMassSquared = piMass*piMass;
   const double protonMass = 0.938272046;
   const double protonMassSquared = protonMass*protonMass;
   const double kShortMass = 0.497614;
   const double lambdaMass = 1.115683;
}

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;
//edm::InputTag m_genParticlesTag_SIM_GEANT;
//edm::EDGetTokenT<std::vector<reco::GenParticle>> m_genParticlesToken_SIM_GEANT;


bool studyOnlySGrandDaughterTracks = false;
//Analyzer_V0Fitter::Analyzer_V0Fitter(const edm::ParameterSet& theParameters, edm::ConsumesCollector && iC){
Analyzer_V0Fitter::Analyzer_V0Fitter(edm::ParameterSet const& theParameters):

   vertexFitter_(theParameters.getParameter<bool>("vertexFitter")),
   useRefTracks_(theParameters.getParameter<bool>("useRefTracks")),
   
   // whether to reconstruct KShorts
   doKShorts_(theParameters.getParameter<bool>("doKShorts")),
   // whether to reconstruct Lambdas
   doLambdas_(theParameters.getParameter<bool>("doLambdas")),

   // cuts on initial track selection
   tkChi2Cut_(theParameters.getParameter<double>("tkChi2Cut")),
   tkNHitsCut_(theParameters.getParameter<int>("tkNHitsCut")),
   tkPtCut_(theParameters.getParameter<double>("tkPtCut")),
   tkIPSigXYCut_(theParameters.getParameter<double>("tkIPSigXYCut")),
   tkIPSigZCut_(theParameters.getParameter<double>("tkIPSigZCut")),
   
   // cuts on vertex
   vtxChi2Cut_(theParameters.getParameter<double>("vtxChi2Cut")),
   vtxDecaySigXYCut_(theParameters.getParameter<double>("vtxDecaySigXYCut")),
   vtxDecaySigXYZCut_(theParameters.getParameter<double>("vtxDecaySigXYZCut")),
   // miscellaneous cuts
   tkDCACut_(theParameters.getParameter<double>("tkDCACut")),
   mPiPiCut_(theParameters.getParameter<double>("mPiPiCut")),
   innerHitPosCut_(theParameters.getParameter<double>("innerHitPosCut")),
   cosThetaXYCut_(theParameters.getParameter<double>("cosThetaXYCut")),
   cosThetaXYZCut_(theParameters.getParameter<double>("cosThetaXYZCut")),
   // cuts on the V0 candidate mass
   kShortMassCut_(theParameters.getParameter<double>("kShortMassCut")),
   lambdaMassCut_(theParameters.getParameter<double>("lambdaMassCut")),
   useVertex_(theParameters.getParameter<bool>("useVertex")),


   m_tag_tracks(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm")),
   m_tag_beamspot(theParameters.getParameter<edm::InputTag>("beamSpot")),
   m_tag_vertices(theParameters.getParameter<edm::InputTag>("vertices")),
   m_tag_genParticlesTag_SIM_GEANT(theParameters.getParameter<edm::InputTag>("genCollection_SIM_GEANT")),
   m_tag_TrackingParticlesTag(theParameters.getParameter<edm::InputTag>("TrackingParticleCollection")),
   
   token_tracks(consumes<reco::TrackCollection>(m_tag_tracks)),
   token_beamspot(consumes<reco::BeamSpot>(m_tag_beamspot)),
   token_vertices(consumes<std::vector<reco::Vertex>>(m_tag_vertices)),
   token_genParticlesTag_SIM_GEANT(consumes<std::vector<reco::GenParticle>>(m_tag_genParticlesTag_SIM_GEANT)),
   token_TrackingParticles(consumes<std::vector<TrackingParticle>>(m_tag_TrackingParticlesTag))

{

   token_associatorByHits= consumes<reco::TrackToTrackingParticleAssociator>(edm::InputTag("associators"));

   token_recoTracks      = consumes<edm::View<reco::Track> >(theParameters.getParameter<edm::InputTag>( "trackRecoAlgorithm" ));

}


void Analyzer_V0Fitter::beginJob() {
 
    //histos_th2f[b+"h_LambdaKshortVertexFilter_S_vx_vy"]= dir_LambdaKshortVertexFilter_S.make<TH2F>(b+"h_LambdaKshortVertexFilter_S_vx_vy", b+"h_LambdaKshortVertexFilter_S_vx_vy; S vx (cm); S vy (cm)", 2000, -100, 100, 2000, -100, 100);
      TFileDirectory dir_general = m_fs->mkdir("general");
      histos_th1f[b+"h_chargeParticles"] = dir_general.make<TH1F>(b+"h_chargeParticles",b+"h_chargeParticles; charge of particles which I try to match to track",100,-50,50);
      histos_th1f[b+"h_delta_eta_pt_GEN_RECO_gen_particle"] = dir_general.make<TH1F>(b+"h_delta_eta_pt_GEN_RECO_gen_particle",b+"h_delta_eta_pt_GEN_RECO_gen_particle; min quadratic delta(eta,pt) between  GEN particle and  track",1000,0,10);
      histos_th1f[b+"h_delta_eta_GEN_RECO_gen_particle"] = dir_general.make<TH1F>(b+"h_delta_eta_GEN_RECO_gen_particle",b+"h_delta_eta_GEN_RECO_gen_particle; delta(eta) between  GEN particle and  track",1000,0,10);
      histos_th1f[b+"h_delta_pt_GEN_RECO_gen_particle"] = dir_general.make<TH1F>(b+"h_delta_pt_GEN_RECO_gen_particle",b+"h_delta_pt_GEN_RECO_gen_particle; delta(pt) between  GEN particle and  track",1000,0,10);
      histos_th1f[b+"h_delta_pt_GEN_RECO_gen_particle_small_delta_eta"] = dir_general.make<TH1F>(b+"h_delta_pt_GEN_RECO_gen_particle_small_delta_eta",b+"h_delta_pt_GEN_RECO_gen_particle_small_delta_eta; delta(pt) between  GEN particle and track for delta_eta < 0.01",1000,0,10);
      histos_th1f[b+"h_delta_pt_over_pt_GEN_RECO_gen_particle"] = dir_general.make<TH1F>(b+"h_delta_pt_over_pt_GEN_RECO_gen_particle",b+"h_delta_pt_over_pt_GEN_RECO_gen_particle; delta(pt) between  GEN particle and  track/pt GEN",1000,0,10);
      histos_th1f[b+"h_delta_pt_over_pt_GEN_RECO_gen_particle_small_delta_eta"] = dir_general.make<TH1F>(b+"h_delta_pt_over_pt_GEN_RECO_gen_particle_small_delta_eta",b+"h_delta_pt_over_pt_GEN_RECO_gen_particle_small_delta_eta; delta(pt) between  GEN particle and  track/pt GEN for delta_eta < 0.01",1000,0,10);
      histos_th1f[b+"h_deltaR_GEN_RECO_gen_particle"] = dir_general.make<TH1F>(b+"h_deltaR_GEN_RECO_gen_particle",b+"h_deltaR_GEN_RECO_gen_particle; deltaR between  GEN particle and  track",1000,0,10);
      histos_th1f[b+"h_delta_eta_pt_min_GEN_RECO_gen_particle"] = dir_general.make<TH1F>(b+"h_delta_eta_pt_min_GEN_RECO_gen_particle",b+"h_delta_eta_pt_min_GEN_RECO_gen_particle; min quadratic delta(eta,pt) between  GEN particle and  track",1000,0,10);
      histos_th1f[b+"h_deltaR_min_GEN_RECO_gen_particle"] = dir_general.make<TH1F>(b+"h_deltaR_min_GEN_RECO_gen_particle",b+"h_deltaR_min_GEN_RECO_gen_particle; min deltaR between  GEN particle and  track",1000,0,10);
      histos_th1f[b+"h_deltaR_min_GEN_RECO_gen_particle_deltaPhi"] = dir_general.make<TH1F>(b+"h_deltaR_min_GEN_RECO_gen_particle_deltaPhi",b+"h_deltaR_min_GEN_RECO_gen_particle_deltaPhi; deltaPhi for min deltaR between  GEN particle and  track",1000,-10,10);
      histos_th1f[b+"h_deltaR_min_GEN_RECO_gen_particle_deltaEta"] = dir_general.make<TH1F>(b+"h_deltaR_min_GEN_RECO_gen_particle_deltaEta",b+"h_deltaR_min_GEN_RECO_gen_particle_deltaEta; deltaEta for min deltaR between  GEN particle and  track",1000,-10,10);
      histos_th1f[b+"h_deltaEta_min_GEN_RECO_gen_particle"] = dir_general.make<TH1F>(b+"h_deltaEta_min_GEN_RECO_gen_particle",b+"h_deltaEta_min_GEN_RECO_gen_particle; min deltaEta between  GEN particle and  track",1000,-10,10);
      histos_th1f[b+"h_deltaPhi_min_GEN_RECO_gen_particle"] = dir_general.make<TH1F>(b+"h_deltaPhi_min_GEN_RECO_gen_particle",b+"h_deltaPhi_min_GEN_RECO_gen_particle; min deltaPhi between  GEN particle and  track",1000,-10,10);
      histos_teff[b+"heff_gen_particle_pt"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_pt",b+"heff_gen_particle_pt; GEN particle pt (GeV) ; Track reconstruction efficiency",100,0,10);
      histos_teff[b+"heff_gen_particle_pt_no_daughter_antiS"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_pt_no_daughter_antiS",b+"heff_gen_particle_pt_no_daughter_antiS; GEN particle no daugther of antiS pt (GeV) ; Track reconstruction efficiency",100,0,10);
      histos_teff[b+"heff_gen_particle_pt_REF_no_daughter_antiS"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_pt_REF_no_daughter_antiS",b+"heff_gen_particle_pt_REF_no_daughter_antiS; GEN particle no daugther of antiS pt (GeV) ; Track reconstruction efficiency",100,0,10);
      histos_teff[b+"heff_gen_particle_pt_REF"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_pt_REF",b+"heff_gen_particle_pt_REF; GEN particle pt (GeV); ; Track reconstruction efficiency",100,0,10);
      histos_teff[b+"heff_gen_particle_pt_REF_ultraTight_lxy_and_vz_cut"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_pt_REF_ultraTight_lxy_and_vz_cut",b+"heff_gen_particle_pt_REF_ultraTight_lxy_and_vz_cut; GEN particle pt (GeV); ; Track reconstruction efficiency",100,0,10);
      histos_teff[b+"heff_gen_particle_pt_REF_ultraTight_lxy_and_vz_cut_and_daughter_displaced"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_pt_REF_ultraTight_lxy_and_vz_cut_and_daughter_displaced",b+"heff_gen_particle_pt_REF_ultraTight_lxy_and_vz_cut_and_daughter_displaced; GEN particle pt (GeV); ; Track reconstruction efficiency",100,0,10);
      histos_teff[b+"heff_gen_particle_pt_REF_status1"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_pt_REF_status1",b+"heff_gen_particle_pt_REF_status1; GEN particle pt (GeV); ; Track reconstruction efficiency",100,0,10);
      histos_teff[b+"heff_gen_particle_vz_deltaR0p01"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_deltaR0p01",b+"heff_gen_particle_vz_deltaR0p01; RECO eff of  particle vs vz (cm)",800,-400,400);
      histos_teff[b+"heff_gen_particle_vz_deltaR0p02"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_deltaR0p02",b+"heff_gen_particle_vz_deltaR0p02; RECO eff of  particle vs vz (cm)",800,-400,400);
     histos_teff[b+"heff_gen_particle_vz_deltaR0p04"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_deltaR0p04",b+"heff_gen_particle_vz_deltaR0p04; RECO eff of  particle vs vz (cm)",800,-400,400);
      histos_teff[b+"heff_gen_particle_vz_deltaR0p05"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_deltaR0p05",b+"heff_gen_particle_vz_deltaR0p05; RECO eff of  particle vs vz (cm)",800,-400,400);
      histos_teff[b+"heff_gen_particle_vz"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz",b+"heff_gen_particle_vz; GEN particle creation vertex vz (cm); Track reconstruction efficiency",800,-400,400);
      histos_teff[b+"heff_gen_particle_vz_lxy_smaller_0p1"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_lxy_smaller_0p1",b+"heff_gen_particle_vz_lxy_smaller_0p1; RECO eff of  particle vs vz with lxy < 0.1 cm (cm)",800,-400,400);
      histos_teff[b+"heff_gen_particle_vz_lxy_0p1_to_0p2"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_lxy_0p1_to_0p2",b+"heff_gen_particle_vz_lxy_0p1_to_0p2; RECO eff of  particle vs vz with 0.1 =< lxy < 0.2 cm (cm)",800,-400,400);
      histos_teff[b+"heff_gen_particle_vz_lxy_0p2_to_0p5"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_lxy_0p2_to_0p5",b+"heff_gen_particle_vz_lxy_0p2_to_0p5; RECO eff of  particle vs vz with 0.2 =< lxy < 0.5 cm (cm)",800,-400,400);
      histos_teff[b+"heff_gen_particle_vz_lxy_0p5_to_1"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_lxy_0p5_to_1",b+"heff_gen_particle_vz_lxy_0p5_to_1; RECO eff of  particle vs vz with 0.5 =< lxy < 1 cm (cm)",800,-400,400);
      histos_teff[b+"heff_gen_particle_vz_lxy_1_to_5"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_lxy_1_to_5",b+"heff_gen_particle_vz_lxy_1_to_5; RECO eff of  particle vs vz with 1 =< lxy < 5 cm (cm)",800,-400,400);
      histos_teff[b+"heff_gen_particle_vz_lxy_larger_5"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_vz_lxy_larger_5",b+"heff_gen_particle_vz_lxy_larger_5; RECO eff of  particle vs vz with lxy > 5  cm (cm)",800,-400,400);

      histos_th1f[b+"h_trackMatched_large_vz_dxy"] = dir_general.make<TH1F>(b+"h_trackMatched_large_vz_dxy",b+"h_trackMatched_large_vz_dxy; dxy(cm)",2000,-10,10);
      histos_th1f[b+"h_trackMatched_large_vz_dz"] = dir_general.make<TH1F>(b+"h_trackMatched_large_vz_dz",b+"h_trackMatched_large_vz_dz; dz(cm)",1000,0,10);
      histos_th1f[b+"h_trackMatched_large_vz_lxy"] = dir_general.make<TH1F>(b+"h_trackMatched_large_vz_lxy",b+"h_trackMatched_large_vz_lxy; lxy(cm)",1000,0,10);
      histos_th1f[b+"h_trackNotMatched_large_vz_dxy"] = dir_general.make<TH1F>(b+"h_trackNotMatched_large_vz_dxy",b+"h_trackNotMatched_large_vz_dxy; dxy(cm)",2000,-10,10);
      histos_th1f[b+"h_trackNotMatched_large_vz_dz"] = dir_general.make<TH1F>(b+"h_trackNotMatched_large_vz_dz",b+"h_trackNotMatched_large_vz_dz; dz(cm)",1000,0,10);
      histos_th1f[b+"h_trackNotMatched_large_vz_lxy"] = dir_general.make<TH1F>(b+"h_trackNotMatched_large_vz_lxy",b+"h_trackNotMatched_large_vz_lxy; lxy(cm)",1000,0,10);
      //efficiency plots for REF plots
      histos_teff[b+"heff_REF_tracks_pt"] = dir_general.make<TEfficiency>(b+"heff_REF_tracks_pt",b+"heff_REF_tracks_pt; GEN particle pt(GeV);  Track reconstruction efficiency",100,0,10);
      histos_teff[b+"heff_REF_tracks_pt_only_lxy_cut"] = dir_general.make<TEfficiency>(b+"heff_REF_tracks_pt_only_lxy_cut",b+"heff_REF_tracks_pt_only_lxy_cut; pt(GeV)",100,0,10);
      histos_teff[b+"heff_REF_tracks_pt_only_eta_cut"] = dir_general.make<TEfficiency>(b+"heff_REF_tracks_pt_only_eta_cut",b+"heff_REF_tracks_pt_only_eta_cut; pt(GeV)",100,0,10);
      histos_teff[b+"heff_REF_tracks_pt_loose_lxy"] = dir_general.make<TEfficiency>(b+"heff_REF_tracks_pt_loose_lxy",b+"heff_REF_tracks_pt_loose_lxy; pt(GeV)",100,0,10);
      histos_teff[b+"heff_REF_tracks_pt_medium_lxy"] = dir_general.make<TEfficiency>(b+"heff_REF_tracks_pt_medium_lxy",b+"heff_REF_tracks_pt_medium_lxy; pt(GeV)",100,0,10);
      histos_teff[b+"heff_REF_tracks_pt_tight_lxy"] = dir_general.make<TEfficiency>(b+"heff_REF_tracks_pt_tight_lxy",b+"heff_REF_tracks_pt_tight_lxy; pt(GeV)",100,0,10);
      histos_teff[b+"heff_REF_tracks_pt_tighter_lxy"] = dir_general.make<TEfficiency>(b+"heff_REF_tracks_pt_tighter_lxy",b+"heff_REF_tracks_pt_tighter_lxy; pt(GeV)",100,0,10);
      histos_teff[b+"heff_REF_tracks_eta"] = dir_general.make<TEfficiency>(b+"heff_REF_tracks_eta",b+"heff_REF_tracks_eta; GEN particle #eta; Track reconstruction efficiency",100,-4,4);
      histos_teff[b+"heff_REF_tracks_lxy"] = dir_general.make<TEfficiency>(b+"heff_REF_tracks_lxy",b+"heff_REF_tracks_lxy; GEN particle lxy(beamspot, GEN creation vertex)(cm); Track reconstruction efficiency",400,0,400);
 
      histos_teff[b+"heff_gen_particle_lxy"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_lxy",b+"heff_gen_particle_lxy;  GEN particle lxy(beamspot, GEN creation vertex)(cm); Track reconstruction efficiency",800,0,400);
     
      histos_th1f[b+"h_trackMatched_gen_pt"] = dir_general.make<TH1F>(b+"h_trackMatched_gen_pt",b+"h_trackMatched_gen_pt; pt(GeV)",1000,0,10);
      histos_th2f[b+"h_gen_pt_Delta_R_track_min"] = dir_general.make<TH2F>(b+"h_gen_pt_Delta_R_track_min",b+"h_gen_pt_Delta_R_track_min; pt(GeV); deltaR(RECO,GEN)",100,0,10,1000,0,10);
      histos_th2f[b+"h_gen_lxy_Delta_R_track_min"] = dir_general.make<TH2F>(b+"h_gen_lxy_Delta_R_track_min",b+"h_gen_lxy_Delta_R_track_min; lxy(cm); deltaR(RECO,GEN)",1000,0,100,1000,0,10);
      histos_th2f[b+"h_gen_lxy_Delta_R_track_min_deltaPhi"] = dir_general.make<TH2F>(b+"h_gen_lxy_Delta_R_track_min_deltaPhi",b+"h_gen_lxy_Delta_R_track_min_deltaPhi; lxy(cm); deltaPhi for min deltaR(RECO,GEN)",1000,0,100,1000,-10,10);
      histos_th2f[b+"h_gen_lxy_Delta_R_track_min_deltaEta"] = dir_general.make<TH2F>(b+"h_gen_lxy_Delta_R_track_min_deltaEta",b+"h_gen_lxy_Delta_R_track_min_deltaEta; lxy(cm); deltaEta for min deltaR(RECO,GEN)",1000,0,100,1000,-10,10);
      histos_th2f[b+"h_gen_lxy_track_min_deltaEta"] = dir_general.make<TH2F>(b+"h_gen_lxy_track_min_deltaEta",b+"h_gen_lxy_track_min_deltaEta; lxy(cm);min  deltaEta (RECO,GEN)",1000,0,100,1000,-10,10);
      histos_th2f[b+"h_gen_lxy_track_min_deltaPhi"] = dir_general.make<TH2F>(b+"h_gen_lxy_track_min_deltaPhi",b+"h_gen_lxy_track_min_deltaPhi; lxy(cm);min  deltaPhi (RECO,GEN)",1000,0,100,1000,-10,10);
      histos_th1f[b+"h_gen_lxy_status_1"] = dir_general.make<TH1F>(b+"h_gen_lxy_status_1",b+"h_gen_lxy_status_1; status 1 particles lxy(cm)", 400000,0,400);
      histos_th1f[b+"h_gen_vz_status_1"] = dir_general.make<TH1F>(b+"h_gen_vz_status_1",b+"h_gen_vz_status_1; status 1 particles vz(cm)", 8000,-400,400);
      histos_th1f[b+"h_gen_lxy_status_1_ZeroZeroZero"] = dir_general.make<TH1F>(b+"h_gen_lxy_status_1_ZeroZeroZero",b+"h_gen_lxy_status_1_ZeroZeroZero; status 1 particles lxy(cm) to ZeroZeroZero", 4000,0,400);
      histos_th1f[b+"h_gen_lxy_status_8"] = dir_general.make<TH1F>(b+"h_gen_lxy_status_8",b+"h_gen_lxy_status_8; status 8 particles lxy(cm)", 4000,0,400);
      histos_th1f[b+"h_gen_vz_status_8"] = dir_general.make<TH1F>(b+"h_gen_vz_status_8",b+"h_gen_vz_status_8; status 8 particles vz(cm)", 8000,-400,400);
      histos_th1f[b+"h_gen_lxy_status_other"] = dir_general.make<TH1F>(b+"h_gen_lxy_status_other",b+"h_gen_lxy_status_other; status != 1 and !=8  particles lxy(cm)", 4000,0,400);
      histos_th2f[b+"h_gen_lxyz_Delta_R_track_min"] = dir_general.make<TH2F>(b+"h_gen_lxyz_Delta_R_track_min",b+"h_gen_lxyz_Delta_R_track_min; lxyz(cm); deltaR(RECO,GEN)",5000,0,500,1000,0,10);
      histos_th1f[b+"h_trackMatched_gen_eta"] = dir_general.make<TH1F>(b+"h_trackMatched_gen_eta",b+"h_trackMatched_gen_eta; eta",100,-4,4);
      histos_th1f[b+"h_trackMatched_gen_phi"] = dir_general.make<TH1F>(b+"h_trackMatched_gen_phi",b+"h_trackMatched_gen_phi; phi(rad)",100,-4,4);
      histos_th1f[b+"h_trackMatched_gen_vxy"] = dir_general.make<TH1F>(b+"h_trackMatched_gen_vxy",b+"h_trackMatched_gen_vxy; vxy(cm)",1000,0,100);
      histos_th1f[b+"h_trackMatched_gen_lxy"] = dir_general.make<TH1F>(b+"h_trackMatched_gen_lxy",b+"h_trackMatched_gen_lxy; lxy(beamspot, creation vertex)(cm)",1000,0,100);
      histos_th1f[b+"h_trackMatched_gen_dxy"] = dir_general.make<TH1F>(b+"h_trackMatched_gen_dxy",b+"h_trackMatched_gen_dxy; dxy(cm)",1000,0,100);
      histos_th1f[b+"h_trackMatched_gen_vz"] = dir_general.make<TH1F>(b+"h_trackMatched_gen_vz",b+"h_trackMatched_gen_vz; vz(cm)",8000,-400,400);
      histos_th2f[b+"h_trackMatched_gen_vz_lxy"] = dir_general.make<TH2F>(b+"h_trackMatched_gen_vz_lxy",b+"h_trackMatched_gen_vz_lxy; vz(cm);lxy(cm)",800,-400,400,400,0,400);

     histos_teff[b+"heff_gen_particle_with_mother_2_daughters_and_mother_V0_pt"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_with_mother_2_daughters_and_mother_V0_pt",b+"heff_gen_particle_with_mother_2_daughters_and_mother_V0_pt; RECO eff of a particle with a mother that has 2 daughters and the mother is a V0 vs pt",1000,0,10);
      histos_teff[b+"heff_gen_particle_with_grandMother_antiS_pt"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_with_grandMother_antiS_pt",b+"heff_gen_particle_with_grandMother_antiS_pt; RECO eff of the GEN with grandmother antiS wrt pt (GeV)",1000,0,10);
      histos_teff[b+"heff_gen_particle_with_grandMother_antiS_vz"] = dir_general.make<TEfficiency>(b+"heff_gen_particle_with_grandMother_antiS_vz",b+"heff_gen_particle_with_grandMother_antiS_vz; RECO eff of the GEN with grandmother antiS wrt vz (cm)",80,-400,400);


     histos_th1f[b+"h_grandDaughterOfAntiS_gen_pt"] = dir_general.make<TH1F>(b+"h_grandDaughterOfAntiS_gen_pt",b+"h_grandDaughterOfAntiS_gen_pt; grandDaughter of antiS pt",1000,0,10);
     histos_th1f[b+"h_grandDaughterOfAntiS_gen_eta"] = dir_general.make<TH1F>(b+"h_grandDaughterOfAntiS_gen_eta",b+"h_grandDaughterOfAntiS_gen_eta; grandDaughter of antiS eta",100,-4,4);
     histos_th1f[b+"h_grandDaughterOfAntiS_gen_phi"] = dir_general.make<TH1F>(b+"h_grandDaughterOfAntiS_gen_phi",b+"h_grandDaughterOfAntiS_gen_phi; grandDaughter of antiS phi",100,-4,4);
     histos_th1f[b+"h_grandDaughterOfAntiS_gen_vxy"] = dir_general.make<TH1F>(b+"h_grandDaughterOfAntiS_gen_vxy",b+"h_grandDaughterOfAntiS_gen_vxy; grandDaughter of antiS vxy (cm)",1000,0,100);
     histos_th1f[b+"h_grandDaughterOfAntiS_gen_lxy"] = dir_general.make<TH1F>(b+"h_grandDaughterOfAntiS_gen_lxy",b+"h_grandDaughterOfAntiS_gen_lxy; grandDaughter of antiS lxy (cm)",1000,0,100);
     histos_th1f[b+"h_grandDaughterOfAntiS_gen_dxy"] = dir_general.make<TH1F>(b+"h_grandDaughterOfAntiS_gen_dxy",b+"h_grandDaughterOfAntiS_gen_dxy; grandDaughter of antiS dxy (cm)",2000,-100,100);
     histos_th1f[b+"h_grandDaughterOfAntiS_gen_dz"] = dir_general.make<TH1F>(b+"h_grandDaughterOfAntiS_gen_dz",b+"h_grandDaughterOfAntiS_gen_dz; grandDaughter of antiS dz (cm)",2000,-100,100);
     histos_th1f[b+"h_grandDaughterOfAntiS_gen_vz"] = dir_general.make<TH1F>(b+"h_grandDaughterOfAntiS_gen_vz",b+"h_grandDaughterOfAntiS_gen_vz; grandDaughter of antiS creation vertex vz (cm)",8000,-400,400);
     histos_th2f[b+"h_grandDaughterOfAntiS_gen_vz_lxy"] = dir_general.make<TH2F>(b+"h_grandDaughterOfAntiS_gen_vz_lxy",b+"h_grandDaughterOfAntiS_gen_vz_lxy;grandDaughter of antiS vz (cm); grandDaughter of antiS lxy(cm)",800,-400,400,400, 0,400);

     histos_th1f[b+"h_noGrandDaughterOfAntiS_gen_vz"] = dir_general.make<TH1F>(b+"h_noGrandDaughterOfAntiS_gen_vz",b+"h_noGrandDaughterOfAntiS_gen_vz; no grandDaughter of antiS vz (cm)",8000,-400,400);


      histos_th1f[b+"h_deltaR_GEN_RECO_gen_particle_with_grand_mother_antiS"] = dir_general.make<TH1F>(b+"h_deltaR_GEN_RECO_gen_particle_with_grand_mother_antiS",b+"h_deltaR_GEN_RECO_gen_particle_with_grand_mother_antiS; deltaR between GEN particles with as grandmother the antiS",1000,0,10);
      histos_th1f[b+"h_deltaR_V0_daughter_GEN_RECO_grandDaughterS"] = dir_general.make<TH1F>(b+"h_deltaR_V0_daughter_GEN_RECO_grandDaughterS",b+"h_deltaR_V0_daughter_GEN_RECO_grandDaughterS; deltaR between  antiSGrandaughter GEN particle and  track",1000,0,10);

      histos_th1f[b+"h_nTracks"] = dir_general.make<TH1F>(b+"h_Tracks",b+"h_nTracks; number of antiS in events",1000,0,1000);
      histos_th1f[b+"h_nAntiSEvent"] = dir_general.make<TH1F>(b+"h_nAntiSEvent",b+"h_nAntiSEvent; number of antiS in events",5,0,5);
      histos_th1f[b+"h_nDaughtersAntiS"] = dir_general.make<TH1F>(b+"h_nDaughtersAntiS",b+"h_nDaughtersAntiS; number of antiS Daughters in events",10,0,10);
      histos_th1f[b+"h_nGrandDaughtersAntiS"] = dir_general.make<TH1F>(b+"h_nGrandDaughtersAntiS",b+"h_nGrandDaughtersAntiS; number of antiS grandDaughters in events",20,0,20);
      histos_th1f[b+"h_nAntiSWithReconstructableGrandDaughters"] = dir_general.make<TH1F>(b+"h_nAntiSWithReconstructableGrandDaughters",b+"h_nAntiSWithReconstructableGrandDaughters; number of antiS with reconstructable granddaughters (this is 4 tracks with all the expected charge and all in tracker acceptance)",2,0,2);


      TFileDirectory dir_cuts_single_tracks = m_fs->mkdir("cuts_single_track");
      histos_th1f[b+"h_deltaR_V0_daughter_GEN_RECO"] = dir_cuts_single_tracks.make<TH1F>(b+"h_deltaR_V0_daughter_GEN_RECO",b+"h_deltaR_V0_daughter_GEN_RECO; h_deltaR_V0_daughter_GEN_RECO",1000,0,10);
      histos_th1f[b+"h_ReconstructableGrandDaughtersMatchedToTrack"] = dir_cuts_single_tracks.make<TH1F>(b+"h_ReconstructableGrandDaughtersMatchedToTrack",b+"h_ReconstructableGrandDaughtersMatchedToTrack; #reconstructable granddaughters that were effectively reconstructed to a track (only the ones with 4 could possibly lead to the RECO of an antiS)",10,0,10);

      histos_th1f[b+"h_track_normalizedChi2"] = dir_cuts_single_tracks.make<TH1F>(b+"h_track_normalizedChi2",b+"h_track_normalizedChi2; h_track_normalizedChi2",100,0,100);
      histos_th1f[b+"h_track_numberOfValidHits"] = dir_cuts_single_tracks.make<TH1F>(b+"h_track_numberOfValidHits",b+"h_track_numberOfValidHits; h_track_numberOfValidHits",50,0,50);
      histos_th1f[b+"h_track_pt"] = dir_cuts_single_tracks.make<TH1F>(b+"h_track_pt",b+"h_track_pt; h_track_pt",1000,0,100);
      histos_th1f[b+"h_ipsigXY"] = dir_cuts_single_tracks.make<TH1F>(b+"h_ipsigXY",b+"h_ipsigXY; h_ipsigXY",200,-100,100);
      histos_th1f[b+"h_ipsigZ"] = dir_cuts_single_tracks.make<TH1F>(b+"h_ipsigZ",b+"h_ipsigZ;h_ipsigZ",200,-100,100);
      histos_th1f[b+"h_goodTrackSelection"] = dir_cuts_single_tracks.make<TH1F>(b+"h_goodTrackSelection",b+"h_goodTrackSelection; good track = 1",2,0,2);

      TFileDirectory dir_cuts_two_tracks = m_fs->mkdir("cuts_two_tracks");
      histos_th1f[b+"h_oppositeCharge"] = dir_cuts_two_tracks.make<TH1F>(b+"h_oppositeCharge",b+"h_oppositeCharge;1 = tracks opposite charge",2,0,2);
      histos_th1f[b+"h_impactPointTSCP_valid"] = dir_cuts_two_tracks.make<TH1F>(b+"h_impactPointTSCP_valid",b+"h_impactPointTSCP_valid;1 = valid",2,0,2);
      histos_th1f[b+"h_cApp_status"] = dir_cuts_two_tracks.make<TH1F>(b+"h_cApp_status",b+"h_cApp_status;1 = status ok",2,0,2);
      histos_th1f[b+"h_dca"] = dir_cuts_two_tracks.make<TH1F>(b+"h_dca",b+"h_dca; dca",200,-100,100);
      histos_th1f[b+"h_POCA_sensVolume_xy"] = dir_cuts_two_tracks.make<TH1F>(b+"h_POCA_sensVolume_xy",b+"h_POCA_sensVolume_xy; POCA xy",100,0,100);
      histos_th1f[b+"h_POCA_sensVolume_z"] = dir_cuts_two_tracks.make<TH1F>(b+"h_POCA_sensVolume_z",b+"h_POCA_sensVolume_z; POCA z",100,0,100);
      histos_th1f[b+"h_TSCP_valid"] = dir_cuts_two_tracks.make<TH1F>(b+"h_TSCP_valid",b+"h_TSCP_valid;TSCP valid",2,0,2);
      histos_th1f[b+"h_TSCP_directions"] = dir_cuts_two_tracks.make<TH1F>(b+"h_TSCP_directions",b+"h_TSCP_directions;1 = opposite direction",2,0,2);
      histos_th1f[b+"h_massPiPi"] = dir_cuts_two_tracks.make<TH1F>(b+"h_massPiPi",b+"h_massPiPi; PiPiMass",100,0,100);
      histos_th1f[b+"h_theRecoVertexValid"] = dir_cuts_two_tracks.make<TH1F>(b+"h_theRecoVertexValid",b+"h_theRecoVertexValid; theRecoVertexValid = 1",2,0,2);
      histos_th1f[b+"h_vtxNormalizedChi2"] = dir_cuts_two_tracks.make<TH1F>(b+"h_vtxNormalizedChi2",b+"h_vtxNormalizedChi2; h_vtxNormalizedChi2",100,0,100);
      histos_th1f[b+"h_distMagXY_over_sigmaDistMagXY"] = dir_cuts_two_tracks.make<TH1F>(b+"h_distMagXY_over_sigmaDistMagXY",b+"h_distMagXY_over_sigmaDistMagXY; h_distMagXY_over_sigmaDistMagXY",100,0,100);
      histos_th1f[b+"h_distMagXYZ_over_sigmaDistMagXYZ"] = dir_cuts_two_tracks.make<TH1F>(b+"h_distMagXYZ_over_sigmaDistMagXYZ",b+"h_distMagXYZ_over_sigmaDistMagXYZ; h_distMagXYZ_over_sigmaDistMagXYZ",100,0,100);
      histos_th1f[b+"h_vertexRadius_within_inner_track1"] = dir_cuts_two_tracks.make<TH1F>(b+"h_vertexRadius_within_inner_track1",b+"h_vertexRadius_within_inner_track1; Valid = 1",2,0,2);
      histos_th1f[b+"h_vertexRadius_within_inner_track2"] = dir_cuts_two_tracks.make<TH1F>(b+"h_vertexRadius_within_inner_track2",b+"h_vertexRadius_within_inner_track2; Valid = 1",2,0,2);
      histos_th1f[b+"h_RefTracksValid"] = dir_cuts_two_tracks.make<TH1F>(b+"h_RefTracksValid",b+"h_RefTracksValid; Valid = 1",2,0,2);
      histos_th1f[b+"h_trajPlus_trajMins_valid"] = dir_cuts_two_tracks.make<TH1F>(b+"h_trajPlus_trajMins_valid",b+"h_trajPlus_trajMins_valid; Valid = 1",2,0,2);
      histos_th1f[b+"h_angleXY"] = dir_cuts_two_tracks.make<TH1F>(b+"h_angleXY",b+"h_angleXY; h_angleXY",200,-100,100);
      histos_th1f[b+"h_angleXYZ"] = dir_cuts_two_tracks.make<TH1F>(b+"h_angleXYZ",b+"h_angleXYZ; h_angleXYZ",200,-100,100);
      histos_th1f[b+"h_V0_formed"] = dir_cuts_two_tracks.make<TH1F>(b+"h_V0_formed",b+"h_V0_formed; 1=Ks formed, 2=Lambda formed, 3=Lambdabar formed",4,0,4);
      
      
}


void Analyzer_V0Fitter::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

/*
  //beamspot
  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

  //gen particles 
  edm::Handle<vector<reco::GenParticle>> h_genParticles_GEN;
  iEvent.getByToken(m_genParticlesToken_GEN, h_genParticles_GEN);
 
  //SIM GEANT particles
*/
  edm::Handle<vector<reco::GenParticle>> h_genParticles_SIM_GEANT;
  iEvent.getByToken(token_genParticlesTag_SIM_GEANT, h_genParticles_SIM_GEANT);

  edm::Handle<reco::TrackCollection> theTrackHandle;
  iEvent.getByToken(token_tracks, theTrackHandle);
  //if (!theTrackHandle->size()) return;
  const reco::TrackCollection* theTrackCollection = theTrackHandle.product();  

  edm::Handle<edm::View<reco::Track> > trackCollectionH;
  iEvent.getByToken(token_recoTracks, trackCollectionH);
  const edm::View<reco::Track> & recoTrackCollection = *(trackCollectionH.product());

  edm::Handle<vector<TrackingParticle>> h_trackingParticles;
  iEvent.getByToken(token_TrackingParticles, h_trackingParticles);


  edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
  iEvent.getByToken(token_associatorByHits, theAssociator);
  
  if(!trackCollectionH.isValid()) cout << "trackCollectionH is not valid" << endl;
  if(!h_trackingParticles.isValid()) cout << "h_trackingParticles is not valid" << endl;
  if(!theAssociator.isValid()) cout << "theAssociator collection is not valid" << endl;
  
  reco::RecoToSimCollection  recSimColl;

  if(trackCollectionH.isValid() && h_trackingParticles.isValid()){
	cout << "trackCollectionH.isValid() && h_trackingParticles.isValid()" << endl;
	cout << "trackCollectionH size: " << trackCollectionH->size() << endl;
	cout << "h_trackingParticles size: " << h_trackingParticles->size() << endl;
	if(theAssociator.isValid()){
		cout << "theAssociator is valid" << endl;
  		reco::RecoToSimCollection  recSimColl = theAssociator->associateRecoToSim (trackCollectionH, h_trackingParticles);
	}
  }
  cout << "recSimColl size: " << recSimColl.size() << endl;

  for(unsigned int i = 0; i < trackCollectionH->size(); i++){
	  edm::RefToBase<reco::Track> track_test(trackCollectionH, i);
	  if(recSimColl.find(track_test) != recSimColl.end()) {
		  std::vector<std::pair<TrackingParticleRef, double> > tp =  recSimColl[track_test];
		  if(tp.size()!=0) cout << "association quality: " << tp.begin()->second << endl;
	  }
	  else{
		cout<< "reco::Track #" << int(i) << " with pt = " << track_test->pt()
					     << " NOT associated to any TrackingParticle" << "\n";
	  }
  }
  
  

/*  for (unsigned int i = 0; i < h_trackingParticles->size(); ++i){
	for (unsigned int j = 0; j < h_trackingParticles->size(); ++j) {
		TrackingParticleRef tp (trackingParticleHandle, j);
		std::cout << Associations<edm::RefToBase<reco::Track> >("Track",  tp, recSimColl);
	}
  }
*/
/*
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
*/
   std::cout << "------------------------------------------------------------------NEW EVENT----------------------------------------------------------------------------------------"<< std::endl;
   using std::vector;

//SIM GEANT particles
//  edm::Handle<vector<reco::GenParticle>> h_genParticles_SIM_GEANT;
//  iEvent.getByToken(m_genParticlesToken_SIM_GEANT, h_genParticles_SIM_GEANT);

   //edm::Handle<vector<reco::GenParticle> > h_genParticles_SIM_GEANT;
   //iEvent.getByToken(token_genParticlesTag_SIM_GEANT, h_genParticles_SIM_GEANT);
  // if (!theGenPartilclesPlusGeantHandle->size()) return;
  // const std::vector<reco::GenParticle>* theGenPartilclesPlusGeantCollection = theGenPartilclesPlusGeantHandle.product();   




   //matching of GEN to reco tracks based on hit info
  /* for(View<reco::Track>::size_type i_track=0; i_track<trackCollectionH->size(); ++i_track){
   	RefToBase<reco::Track> track(trackCollectionH, i_track);
   	std::vector<std::pair<TrackingParticleRef, double> > tp;
   	if(recSimColl.find(track) != recSimColl.end()){
     		tp = recSimColl[track];
     		if (tp.size()!=0) {
              		TrackingParticleRef tpr = tp.begin()->first;
              		double associationQuality = tp.begin()->second;
              		cout << "associationQuality " << associationQuality << endl; 
          	}
        }
   }
*/


   edm::Handle<reco::BeamSpot> theBeamSpotHandle;
   iEvent.getByToken(token_beamspot, theBeamSpotHandle);
   const reco::BeamSpot* theBeamSpot = theBeamSpotHandle.product();
   math::XYZPoint referencePos(theBeamSpot->position());
   std::cout << "--V0Fitter.cc--: location of the referencePos (i.e. the beamspot): " << referencePos.x() << ", " << referencePos.y() << ", " << referencePos.z() << std::endl; 
   reco::Vertex referenceVtx;
   if (useVertex_) {
      edm::Handle<std::vector<reco::Vertex>> vertices;
      iEvent.getByToken(token_vertices, vertices);
      referenceVtx = vertices->at(0);
      referencePos = referenceVtx.position();
   }

   edm::ESHandle<MagneticField> theMagneticFieldHandle;
   iSetup.get<IdealMagneticFieldRecord>().get(theMagneticFieldHandle);
   const MagneticField* theMagneticField = theMagneticFieldHandle.product();




   //loop over the theGenPartilclesPlusGeantCollection and check for an antiS, check for it to have 2x2 daughters in tk accepatance and with the correct pdgid, if found, print the parameters of these 4 particles
   bool antiSWithReconstructableKsGrandDaughters = false;
   bool antiSWithReconstructableLGrandDaughters = false;
   bool antiSWithReconstructableGrandDaughters = false;
   vector<vector<double>> V0_daughters_phi_eta;
   int nAntiSInEvent = 0;
   for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i){
	if(h_genParticles_SIM_GEANT->at(i).pdgId() == -1020000020){
		std::cout << "FOUND AN ANTIS IN THE GENPARTICLESPLUSGEANT" << std::endl;
		nAntiSInEvent++;
	}
	else continue;
 	cout << "h_genParticles_SIM_GEANT->at(i).numberOfDaughters()" << endl;	
 	cout << h_genParticles_SIM_GEANT->at(i).numberOfDaughters() << endl;	
	//check how m of the 

	histos_th1f[b+"h_nDaughtersAntiS"]->Fill(h_genParticles_SIM_GEANT->at(i).numberOfDaughters());	    
	if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() == 2 ){
	    histos_th1f[b+"h_nGrandDaughtersAntiS"]->Fill(h_genParticles_SIM_GEANT->at(i).daughter(0)->numberOfDaughters() + h_genParticles_SIM_GEANT->at(i).daughter(1)->numberOfDaughters());	    
	    if(h_genParticles_SIM_GEANT->at(i).daughter(0)->numberOfDaughters() == 2 && h_genParticles_SIM_GEANT->at(i).daughter(1)->numberOfDaughters() == 2){
	 		const reco::Candidate *daughter0 =  h_genParticles_SIM_GEANT->at(i).daughter(0);
                        const reco::Candidate *daughter1 = h_genParticles_SIM_GEANT->at(i).daughter(1);
                        //finding out which daughter is which
                        const reco::Candidate *Ks = NULL;
                        const reco::Candidate *L = NULL;
                        if(daughter0->pdgId() == 310) Ks = daughter0;
                        else if(daughter0->pdgId() == -3122) L = daughter0;
                        else std::cout << "DAUGHTER0 FROM THE anti-S IS NOT KS OR LAMBDA" << std::endl;
                        
                        if(daughter1->pdgId() == 310) Ks = daughter1;
                        else if(daughter1->pdgId() == -3122) L = daughter1;
                        else std::cout << "DAUGHTER1 FROM THE anti-S IS NOT KS OR LAMBDA" << std::endl;	

			bool KsDaughtersInTrackerAcceptance = false;
			if(fabs(Ks->daughter(0)->eta()) < 2.5 && fabs(Ks->daughter(1)->eta()) < 2.5) KsDaughtersInTrackerAcceptance = true;

			bool KsDaughtersAreChargedPions = false;
			if(  (Ks->daughter(0)->pdgId() == 211 && Ks->daughter(1)->pdgId() == -211 ) || (Ks->daughter(0)->pdgId() == -211 && Ks->daughter(1)->pdgId() == 211 )     ) KsDaughtersAreChargedPions = true;
		        bool LDaughtersInTrackerAcceptance = false;
                        if(fabs(L->daughter(0)->eta()) < 2.5 && fabs(L->daughter(1)->eta()) < 2.5) LDaughtersInTrackerAcceptance = true;

			bool LDaughtersAreChargedPions = false;
			if(  (L->daughter(0)->pdgId() == -2212 && L->daughter(1)->pdgId() == 211 ) || (L->daughter(0)->pdgId() == 211 && L->daughter(1)->pdgId() == -2212 )     ) LDaughtersAreChargedPions = true;
			if(KsDaughtersInTrackerAcceptance && KsDaughtersAreChargedPions) antiSWithReconstructableKsGrandDaughters = true;
			if(LDaughtersInTrackerAcceptance && LDaughtersAreChargedPions) antiSWithReconstructableLGrandDaughters = true;
			if(antiSWithReconstructableKsGrandDaughters && antiSWithReconstructableLGrandDaughters) antiSWithReconstructableGrandDaughters = true;

			if(antiSWithReconstructableGrandDaughters){
				//std::cout << "FOUND AN ANTIS IN THE GENPARTICLESPLUSGEANT WITH DAUGHTERS OF V0S IN TRACKER ACCEPTANCE AND DECAYING TO CHARGED PRODUCTS" << std::endl;
				histos_th1f[b+"h_nAntiSWithReconstructableGrandDaughters"]->Fill(1);

				std::cout << "the daughters of the V0s: " << std::endl;
				std::cout << "Ks daug0 phi, eta: " << Ks->daughter(0)->phi() << " " << Ks->daughter(0)->eta()  << std::endl;
				std::cout << "Ks daug1 phi, eta: " << Ks->daughter(1)->phi() << " " << Ks->daughter(1)->eta()  << std::endl;
				std::cout << "L daug0 phi, eta: " << L->daughter(0)->phi() << " " << L->daughter(0)->eta()  << std::endl;
				std::cout << "L daug1 phi, eta: " << L->daughter(1)->phi() << " " << L->daughter(1)->eta()  << std::endl;
				   
				vector<double> Ks_daug0_phi_eta;
				vector<double> Ks_daug1_phi_eta;
				vector<double> L_daug0_phi_eta;
				vector<double> L_daug1_phi_eta;
				
				Ks_daug0_phi_eta.push_back(Ks->daughter(0)->phi());Ks_daug0_phi_eta.push_back(Ks->daughter(0)->eta());
				Ks_daug1_phi_eta.push_back(Ks->daughter(1)->phi());Ks_daug1_phi_eta.push_back(Ks->daughter(1)->phi());
				L_daug0_phi_eta.push_back(L->daughter(0)->phi());L_daug0_phi_eta.push_back(L->daughter(0)->eta());
				L_daug1_phi_eta.push_back(L->daughter(1)->phi());L_daug1_phi_eta.push_back(L->daughter(1)->eta());
				V0_daughters_phi_eta.push_back(Ks_daug0_phi_eta);
				V0_daughters_phi_eta.push_back(Ks_daug1_phi_eta);
				V0_daughters_phi_eta.push_back(L_daug0_phi_eta);
				V0_daughters_phi_eta.push_back(L_daug1_phi_eta);
			}
			else  histos_th1f[b+"h_nAntiSWithReconstructableGrandDaughters"]->Fill(0);
			
			//additional cuts on the daughters to find an nice event to look at in fireworks: 4 granddaughters also pt > 1 GeV and 4 granddaughters abs(vz) < 10cm
			bool granddaughtersPtHigh = false;
			bool granddaughtersLPtHigh = false;
			bool granddaughtersKsPtHigh = false;
			bool granddaughtersVzSmall = false;
			bool granddaughtersLVzSmall = false;
			bool granddaughtersKsVzSmall = false;
			if(L->daughter(0)->pt() > 1 && L->daughter(1)->pt() > 1 ) granddaughtersLPtHigh = true;
			if(Ks->daughter(0)->pt() > 1 && Ks->daughter(1)->pt() > 1 ) granddaughtersKsPtHigh = true;
			if(granddaughtersLPtHigh && granddaughtersKsPtHigh ) granddaughtersPtHigh = true;
			if(fabs(L->daughter(0)->vz()) < 10 && fabs(L->daughter(1)->vz()) < 10 ) granddaughtersLVzSmall = true;
			if(fabs(Ks->daughter(0)->vz()) < 10 && fabs(Ks->daughter(1)->vz()) < 10  ) granddaughtersKsVzSmall = true;
			if(granddaughtersLVzSmall && granddaughtersKsVzSmall ) granddaughtersVzSmall = true;

			if(antiSWithReconstructableGrandDaughters) cout << "found an antiS with magic1 granddaughters (=4x charged, 4x|eta|<2.5)" << std::endl;
			if(antiSWithReconstructableGrandDaughters && granddaughtersPtHigh) cout << "found an antiS with magic2 granddaughters (=4x charged, 4x|eta|<2.5, 4x pt > 1GeV)" << std::endl;
			if(antiSWithReconstructableGrandDaughters && granddaughtersVzSmall) cout << "found an antiS with magic3 granddaughters (=4x charged, 4x|eta|<2.5, 4x |vz|<10cm)" << std::endl;
			if(antiSWithReconstructableGrandDaughters && granddaughtersPtHigh && granddaughtersVzSmall) cout << "found an antiS with magic4 granddaughters (=4x charged, 4x|eta|<2.5, 4x pt > 1GeV, 4x|vz|<10cm)" << std::endl;

			if(antiSWithReconstructableKsGrandDaughters) cout << "found a Ks with magic11 daughters (=2x charged, 2x|eta|<2.5)" << std::endl;
			if(antiSWithReconstructableKsGrandDaughters && granddaughtersKsPtHigh) cout << "found a Ks with magic12 daughters (=2x charged, 2x|eta|<2.5, 2x pt > 1GeV)" << std::endl;
			if(antiSWithReconstructableKsGrandDaughters && granddaughtersKsVzSmall) cout << "found a Ks with magic13 daughters (=2x charged, 2x|eta|<2.5, 2x |vz|<10cm)" << std::endl;
			if(antiSWithReconstructableKsGrandDaughters && granddaughtersKsPtHigh && granddaughtersKsVzSmall) cout << "found a Ks with magic14 daughters (=2x charged, 2x|eta|<2.5, 2x pt > 1GeV, 2x|vz|<10cm)" << std::endl;

			if(antiSWithReconstructableLGrandDaughters) cout << "found a L with magic21 daughters (=2x charged, 2x|eta|<2.5)" << std::endl;
			if(antiSWithReconstructableLGrandDaughters && granddaughtersLPtHigh) cout << "found a L with magic22 daughters (=2x charged, 2x|eta|<2.5, 2x pt > 1GeV)" << std::endl;
			if(antiSWithReconstructableLGrandDaughters && granddaughtersLVzSmall) cout << "found a L with magic23 daughters (=2x charged, 2x|eta|<2.5, 2x |vz|<10cm)" << std::endl;
			if(antiSWithReconstructableLGrandDaughters && granddaughtersLPtHigh && granddaughtersLVzSmall) cout << "found a L with magic24 daughters (=2x charged, 2x|eta|<2.5, 2x pt > 1GeV, 2x|vz|<10cm)" << std::endl;
			
			
			


	    }//Ks and Lambda enough daughters
	}//antiS enough daughters
	
   }//end of loop over h_genParticles_SIM_GEANT 

   histos_th1f[b+"h_nTracks"]->Fill(theTrackCollection->size());
   histos_th1f[b+"h_nAntiSEvent"]->Fill(nAntiSInEvent);

   std::vector<reco::TrackRef> theTrackRefs;
   std::vector<reco::TransientTrack> theTransTracks;


 
   //check the general track mathing
   if(h_genParticles_SIM_GEANT.isValid()){
   for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i){

	Double_t delta_eta_pt_min = 999;
	Double_t Delta_R_track_min = 999;
	Double_t Delta_R_track_min_phi = 999;
	Double_t Delta_R_track_min_eta = 999;
	Double_t Delta_eta_track_min = 999;
	Double_t Delta_phi_track_min = 999;
	const reco::Track* bestTrack = NULL;

//	if(h_genParticles_SIM_GEANT->at(i).charge() == 0 || fabs(h_genParticles_SIM_GEANT->at(i).eta()) > 2.5)continue;
	if(h_genParticles_SIM_GEANT->at(i).charge() == 0)continue;

	Double_t pt_genParticle = h_genParticles_SIM_GEANT->at(i).pt();
	Double_t phi_genParticle = h_genParticles_SIM_GEANT->at(i).phi();
	cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	cout << "REF postion X,Y (xREF,yREF): "  << referencePos.X() << " , " << referencePos.Y() << endl;
	cout << "GEN postion X,Y (x1,y1): "  << h_genParticles_SIM_GEANT->at(i).vx() << " , " << h_genParticles_SIM_GEANT->at(i).vy() << endl;
	cout << "GEN postion pX,pY (p1x,p1y): "  << h_genParticles_SIM_GEANT->at(i).px() << " , " << h_genParticles_SIM_GEANT->at(i).py() << endl;
	cout << "GEN particle pt: "  <<  h_genParticles_SIM_GEANT->at(i).pt() << endl;
	cout << "normally you use the following phi for the GEN particle (phi1): " << phi_genParticle << endl;
 
	Double_t phi_genParticleCorrected = Analyzer_V0Fitter::phiInTrackRefPosition(h_genParticles_SIM_GEANT->at(i).charge(), h_genParticles_SIM_GEANT->at(i).px(), h_genParticles_SIM_GEANT->at(i).py(), h_genParticles_SIM_GEANT->at(i).vx(), h_genParticles_SIM_GEANT->at(i).vy(), referencePos.X(), referencePos.Y());
	cout << "but you should use this corrected one (phi2), which should match the phi for tracks: " <<  phi_genParticleCorrected << endl;
	//phi_genParticle = phi_genParticleCorrected;

	Double_t eta_genParticle = h_genParticles_SIM_GEANT->at(i).eta();
	for (reco::TrackCollection::const_iterator iTk = theTrackCollection->begin(); iTk != theTrackCollection->end(); ++iTk) {
		const reco::Track* tmpTrack = &(*iTk);
		Double_t ptTrack = tmpTrack->pt();
		Double_t phiTrack = tmpTrack->phi();
		Double_t etaTrack = tmpTrack->eta();

		Double_t deltaPt = ptTrack - pt_genParticle; 
		Double_t deltaEta = etaTrack - eta_genParticle;
		Double_t deltaPhi = reco::deltaPhi(phiTrack,phi_genParticle);

		Double_t delta_eta_pt = pow(pow(deltaEta,2)+pow(deltaPt,2),0.5);
		Double_t deltaR_track = pow(pow(deltaPhi,2)+pow(deltaEta,2),0.5);

		histos_th1f[b+"h_delta_eta_pt_GEN_RECO_gen_particle"]->Fill(delta_eta_pt);
		histos_th1f[b+"h_delta_eta_GEN_RECO_gen_particle"]->Fill(fabs(deltaEta));
		histos_th1f[b+"h_delta_pt_GEN_RECO_gen_particle"]->Fill(fabs(deltaPt));
		histos_th1f[b+"h_delta_pt_over_pt_GEN_RECO_gen_particle"]->Fill(fabs(deltaPt)/pt_genParticle);
		if(deltaEta<0.01) histos_th1f[b+"h_delta_pt_GEN_RECO_gen_particle_small_delta_eta"]->Fill(fabs(deltaPt));
		if(deltaEta<0.01) histos_th1f[b+"h_delta_pt_over_pt_GEN_RECO_gen_particle_small_delta_eta"]->Fill(fabs(deltaPt)/pt_genParticle);
		histos_th1f[b+"h_deltaR_GEN_RECO_gen_particle"]->Fill(deltaR_track);


		if(deltaR_track<Delta_R_track_min){
			Delta_R_track_min=deltaR_track;
			Delta_R_track_min_phi = deltaPhi;
			Delta_R_track_min_eta = deltaEta;
			bestTrack = tmpTrack;
		}
		if(delta_eta_pt < delta_eta_pt_min){
			delta_eta_pt_min = delta_eta_pt;
		}

		if(fabs(deltaEta) < Delta_eta_track_min){
			Delta_eta_track_min = fabs(deltaEta);
		}
		
		if(fabs(deltaPhi) < Delta_phi_track_min){
			Delta_phi_track_min = fabs(deltaPhi);
		}
		

	}//end loop over track collection
	if(bestTrack)cout << "best track reference position: " << bestTrack->vx() << ", " << bestTrack->vy() << "," << bestTrack->vz() << endl;

	bool genParticleGrandDaughterOfAntiS = false;
	if(h_genParticles_SIM_GEANT->at(i).numberOfMothers()>0){
		if(h_genParticles_SIM_GEANT->at(i).mother(0)->numberOfMothers()>0){
			if(h_genParticles_SIM_GEANT->at(i).mother(0)->mother(0)->pdgId()==-1020000020){
				genParticleGrandDaughterOfAntiS = true;
			}
		}//check if mother has mother
	}//check of mother


	TVector3 xyz(h_genParticles_SIM_GEANT->at(i).vx(),h_genParticles_SIM_GEANT->at(i).vy(),h_genParticles_SIM_GEANT->at(i).vz());		
	TVector3 p(h_genParticles_SIM_GEANT->at(i).px(),h_genParticles_SIM_GEANT->at(i).py(),h_genParticles_SIM_GEANT->at(i).pz());
	TVector3 refPos(referencePos.X(),referencePos.Y(),referencePos.Z());
	TVector3 ZeroZeroZero(0,0,0);

	histos_th1f[b+"h_chargeParticles"]->Fill(h_genParticles_SIM_GEANT->at(i).charge());
	histos_th1f[b+"h_delta_eta_pt_min_GEN_RECO_gen_particle"]->Fill(delta_eta_pt_min);
	histos_th1f[b+"h_deltaR_min_GEN_RECO_gen_particle"]->Fill(Delta_R_track_min);
	histos_th1f[b+"h_deltaR_min_GEN_RECO_gen_particle_deltaPhi"]->Fill(Delta_R_track_min_phi);
	histos_th1f[b+"h_deltaR_min_GEN_RECO_gen_particle_deltaEta"]->Fill(Delta_R_track_min_eta);
	histos_th1f[b+"h_deltaEta_min_GEN_RECO_gen_particle"]->Fill(Delta_eta_track_min);
	histos_th1f[b+"h_deltaPhi_min_GEN_RECO_gen_particle"]->Fill(Delta_phi_track_min);
	//below plots to check if DeltaR is by construction larger for displaced particles
	histos_th2f[b+"h_gen_pt_Delta_R_track_min"]->Fill(h_genParticles_SIM_GEANT->at(i).pt(), Delta_R_track_min);
	histos_th2f[b+"h_gen_lxy_Delta_R_track_min"]->Fill(lxy(refPos, xyz), Delta_R_track_min);
	histos_th2f[b+"h_gen_lxy_Delta_R_track_min_deltaPhi"]->Fill(lxy(refPos, xyz), Delta_R_track_min_phi);
	histos_th2f[b+"h_gen_lxy_Delta_R_track_min_deltaEta"]->Fill(lxy(refPos, xyz), Delta_R_track_min_eta);
	histos_th2f[b+"h_gen_lxy_track_min_deltaEta"]->Fill(lxy(refPos, xyz), Delta_eta_track_min);
	histos_th2f[b+"h_gen_lxy_track_min_deltaPhi"]->Fill(lxy(refPos, xyz), Delta_phi_track_min);
	if(h_genParticles_SIM_GEANT->at(i).status() == 1){
		histos_th1f[b+"h_gen_lxy_status_1"]->Fill(lxy(refPos, xyz));
		histos_th1f[b+"h_gen_lxy_status_1_ZeroZeroZero"]->Fill(lxy(ZeroZeroZero, xyz));
		histos_th1f[b+"h_gen_vz_status_1"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
	
	}
	else if(h_genParticles_SIM_GEANT->at(i).status() == 8){
		histos_th1f[b+"h_gen_lxy_status_8"]->Fill(lxy(refPos, xyz));
		histos_th1f[b+"h_gen_vz_status_8"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
	}
	else histos_th1f[b+"h_gen_lxy_status_other"]->Fill(lxy(refPos, xyz));
	histos_th2f[b+"h_gen_lxyz_Delta_R_track_min"]->Fill(pow(lxy(refPos, xyz)*lxy(refPos, xyz)+h_genParticles_SIM_GEANT->at(i).vz()*h_genParticles_SIM_GEANT->at(i).vz(),0.5), Delta_R_track_min);
	bool trackMatched1 = false;
	bool trackMatched2 = false;
	bool trackMatched = false;
	bool trackMatched4 = false;
	bool trackMatched5 = false;
	if(Delta_R_track_min<0.01) trackMatched1 = true;
	if(Delta_R_track_min<0.02) trackMatched2 = true;
	if(Delta_R_track_min<0.03) trackMatched = true;
	if(Delta_R_track_min<0.04) trackMatched4 = true;
	if(Delta_R_track_min<0.05) trackMatched5 = true;

	histos_teff[b+"heff_gen_particle_pt"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).pt());
	if(lxy(refPos,xyz) < 3.5 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5) histos_teff[b+"heff_gen_particle_pt_REF"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).pt());
	if(lxy(refPos,xyz) < 0.005 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5 && fabs(h_genParticles_SIM_GEANT->at(i).vz())<1) histos_teff[b+"heff_gen_particle_pt_REF_ultraTight_lxy_and_vz_cut"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).pt());

	//if particle has a daughter check if the daughter originates far enough: this makes the mother particle going outward of you put a max cut on the lxy(refPos,xyz)
	//by default set decayvertex very large, this is basically the case for particles with no daughter
	TVector3 decayVertexGenParticle(9999,9999,9999);
	if(h_genParticles_SIM_GEANT->at(i).numberOfDaughters() != 0){ 
		decayVertexGenParticle.SetX(h_genParticles_SIM_GEANT->at(i).daughter(0)->vx());
		decayVertexGenParticle.SetY(h_genParticles_SIM_GEANT->at(i).daughter(0)->vy());
		decayVertexGenParticle.SetZ(h_genParticles_SIM_GEANT->at(i).daughter(0)->vz());
	}
 
	if(lxy(refPos,xyz) < 0.005 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5 && fabs(h_genParticles_SIM_GEANT->at(i).vz())<1 && lxy(refPos,decayVertexGenParticle) > 80) histos_teff[b+"heff_gen_particle_pt_REF_ultraTight_lxy_and_vz_cut_and_daughter_displaced"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).pt());


	if(lxy(refPos,xyz) < 3.5 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5 && h_genParticles_SIM_GEANT->at(i).status() == 1) histos_teff[b+"heff_gen_particle_pt_REF_status1"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).pt());

	if(!genParticleGrandDaughterOfAntiS)histos_teff[b+"heff_gen_particle_pt_no_daughter_antiS"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).pt());
	if(!genParticleGrandDaughterOfAntiS && lxy(refPos,xyz) < 3.5 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5)histos_teff[b+"heff_gen_particle_pt_REF_no_daughter_antiS"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).pt());

	histos_teff[b+"heff_gen_particle_vz_deltaR0p01"]->Fill(trackMatched1,h_genParticles_SIM_GEANT->at(i).vz());
	histos_teff[b+"heff_gen_particle_vz_deltaR0p02"]->Fill(trackMatched2,h_genParticles_SIM_GEANT->at(i).vz());
	histos_teff[b+"heff_gen_particle_vz"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).vz());
	histos_teff[b+"heff_gen_particle_vz_deltaR0p04"]->Fill(trackMatched4,h_genParticles_SIM_GEANT->at(i).vz());
	histos_teff[b+"heff_gen_particle_vz_deltaR0p05"]->Fill(trackMatched5,h_genParticles_SIM_GEANT->at(i).vz());
	histos_teff[b+"heff_gen_particle_lxy"]->Fill(trackMatched,lxy(refPos, xyz));

	//check if the efficiency over vz is different for particles with displacement
	if(lxy(refPos, xyz)<0.1)histos_teff[b+"heff_gen_particle_vz_lxy_smaller_0p1"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).vz());
	if(lxy(refPos, xyz)>=0.1 && lxy(refPos, xyz)<0.2)histos_teff[b+"heff_gen_particle_vz_lxy_0p1_to_0p2"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).vz());
	if(lxy(refPos, xyz)>=0.2 && lxy(refPos, xyz)<0.5)histos_teff[b+"heff_gen_particle_vz_lxy_0p2_to_0p5"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).vz());
	if(lxy(refPos, xyz)>=0.5 && lxy(refPos, xyz)<1)histos_teff[b+"heff_gen_particle_vz_lxy_0p5_to_1"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).vz());
	if(lxy(refPos, xyz)>=1 && lxy(refPos, xyz)<5)histos_teff[b+"heff_gen_particle_vz_lxy_1_to_5"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).vz());
	if(lxy(refPos, xyz)>=5)histos_teff[b+"heff_gen_particle_vz_lxy_larger_5"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).vz());

	//check for the particles which do get reconstructed at high vz what is special about them, compared to particles at large vz which do not get reconstructed
	if(trackMatched && fabs(h_genParticles_SIM_GEANT->at(i).vz())>20){
		histos_th1f[b+"h_trackMatched_large_vz_dxy"]->Fill(dxy_signed_line_point(xyz, p, refPos));	
		histos_th1f[b+"h_trackMatched_large_vz_dz"]->Fill(dz_line_point(xyz, p, refPos));	
		histos_th1f[b+"h_trackMatched_large_vz_lxy"]->Fill(lxy(refPos, xyz));	
	}
	if(!trackMatched && fabs(h_genParticles_SIM_GEANT->at(i).vz())>20){
		histos_th1f[b+"h_trackNotMatched_large_vz_dxy"]->Fill(dxy_signed_line_point(xyz, p, refPos));
		histos_th1f[b+"h_trackNotMatched_large_vz_dz"]->Fill(dz_line_point(xyz, p, refPos));	
		histos_th1f[b+"h_trackNotMatched_large_vz_lxy"]->Fill(lxy(refPos, xyz));	
	}

	//for a subset of "reference" (low lxy, eta in tracker acceptance, high enough pt,...) tracks make the efficiency plots
	if(lxy(refPos, xyz)<3.5 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5)histos_teff[b+"heff_REF_tracks_pt"]->Fill(trackMatched, h_genParticles_SIM_GEANT->at(i).pt());
	if(lxy(refPos, xyz)<3.5)histos_teff[b+"heff_REF_tracks_pt_only_lxy_cut"]->Fill(trackMatched, h_genParticles_SIM_GEANT->at(i).pt());
	if(fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5)histos_teff[b+"heff_REF_tracks_pt_only_eta_cut"]->Fill(trackMatched, h_genParticles_SIM_GEANT->at(i).pt());

	if(lxy(refPos, xyz)<0.5 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5)histos_teff[b+"heff_REF_tracks_pt_loose_lxy"]->Fill(trackMatched, h_genParticles_SIM_GEANT->at(i).pt());
	if(lxy(refPos, xyz)<0.1 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5)histos_teff[b+"heff_REF_tracks_pt_medium_lxy"]->Fill(trackMatched, h_genParticles_SIM_GEANT->at(i).pt());
	if(lxy(refPos, xyz)<0.01 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5)histos_teff[b+"heff_REF_tracks_pt_tight_lxy"]->Fill(trackMatched, h_genParticles_SIM_GEANT->at(i).pt());
	if(lxy(refPos, xyz)<0.005 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5)histos_teff[b+"heff_REF_tracks_pt_tighter_lxy"]->Fill(trackMatched, h_genParticles_SIM_GEANT->at(i).pt());
	if(lxy(refPos, xyz)<3.5 && h_genParticles_SIM_GEANT->at(i).pt() > 0.9)histos_teff[b+"heff_REF_tracks_eta"]->Fill(trackMatched, h_genParticles_SIM_GEANT->at(i).eta());
	if(h_genParticles_SIM_GEANT->at(i).pt() > 0.9 && fabs(h_genParticles_SIM_GEANT->at(i).eta())<2.5)histos_teff[b+"heff_REF_tracks_lxy"]->Fill(trackMatched, lxy(refPos, xyz));
	//now for all the tracks you were able to match to a GEN particle with charge plot the properties of the GEN particle linked to this track. 
	if(trackMatched){

		/*
		cout << "found a matched track: " << endl;
		cout << "gen particle lxyz: " << pow(lxy(refPos, xyz)*lxy(refPos, xyz)+h_genParticles_SIM_GEANT->at(i).vz()*h_genParticles_SIM_GEANT->at(i).vz(),0.5) << endl;
		cout << "gen particle lxy: " << lxy(refPos, xyz) << endl;
		cout << "gen particl vz: " << h_genParticles_SIM_GEANT->at(i).vz() << endl;
		cout << "Delta_R_track_min: " << Delta_R_track_min << endl;
		cout << "gen particle vx, vy, vz: " << h_genParticles_SIM_GEANT->at(i).vx() << ", " << h_genParticles_SIM_GEANT->at(i).vy() << ", " << h_genParticles_SIM_GEANT->at(i).vz() << endl;
		cout << "reco track   vx, vy, vz: " << bestTrack->vx() << ", " << bestTrack->vy() << ", " << bestTrack->vz() << endl;
		cout << "gen particle eta, phi " << h_genParticles_SIM_GEANT->at(i).phi() << ", " << h_genParticles_SIM_GEANT->at(i).eta() << endl; 
		cout << "reco track   eta, phi " << bestTrack->phi() << ", " << bestTrack->eta() << endl; 
		cout << "----------------------------------------------------------" << endl;
		*/
		histos_th1f[b+"h_trackMatched_gen_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
		histos_th1f[b+"h_trackMatched_gen_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
		histos_th1f[b+"h_trackMatched_gen_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
		Double_t vx = h_genParticles_SIM_GEANT->at(i).vx();
		Double_t vy = h_genParticles_SIM_GEANT->at(i).vy();
		histos_th1f[b+"h_trackMatched_gen_vxy"]->Fill(pow(vx*vx+vy*vy,0.5));
		histos_th1f[b+"h_trackMatched_gen_lxy"]->Fill(lxy(refPos, xyz));
		histos_th1f[b+"h_trackMatched_gen_dxy"]->Fill(dxy_signed_line_point(xyz, p, refPos)); 	
		histos_th1f[b+"h_trackMatched_gen_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz()); 	
		histos_th2f[b+"h_trackMatched_gen_vz_lxy"]->Fill(lxy(refPos, xyz),h_genParticles_SIM_GEANT->at(i).vz()); 	
	}

	//check the extras you were asking for in the previous script. You were also asking for the particle to be a daughter of a V0 (fabs(310) or fabs(3122)) and this V0 to have 2 daughters
	if(h_genParticles_SIM_GEANT->at(i).numberOfMothers()>0){
		if((fabs(h_genParticles_SIM_GEANT->at(i).mother(0)->pdgId())==3122 || fabs(h_genParticles_SIM_GEANT->at(i).mother(0)->pdgId()==310)) && h_genParticles_SIM_GEANT->at(i).mother(0)->numberOfDaughters()==2 ){
			histos_teff[b+"heff_gen_particle_with_mother_2_daughters_and_mother_V0_pt"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).pt());
		}
	}




	//now for the gen particles which have an AntiS as grandmother check if they are antiProtons or charged pions 
	if(genParticleGrandDaughterOfAntiS && (h_genParticles_SIM_GEANT->at(i).pdgId() == 211 || h_genParticles_SIM_GEANT->at(i).pdgId() == -211 || h_genParticles_SIM_GEANT->at(i).pdgId() == -2212)){

		 Double_t vx = h_genParticles_SIM_GEANT->at(i).vx();
                 Double_t vy = h_genParticles_SIM_GEANT->at(i).vy();
		 Double_t dz = dz_line_point(xyz, p, refPos);
       		 histos_th1f[b+"h_deltaR_GEN_RECO_gen_particle_with_grand_mother_antiS"]->Fill(Delta_R_track_min);
		 histos_th1f[b+"h_grandDaughterOfAntiS_gen_pt"]->Fill(h_genParticles_SIM_GEANT->at(i).pt());
		 histos_th1f[b+"h_grandDaughterOfAntiS_gen_eta"]->Fill(h_genParticles_SIM_GEANT->at(i).eta());
		 histos_th1f[b+"h_grandDaughterOfAntiS_gen_phi"]->Fill(h_genParticles_SIM_GEANT->at(i).phi());
		 histos_th1f[b+"h_grandDaughterOfAntiS_gen_vxy"]->Fill(pow(vx*vx+vy*vy,0.5));
		 histos_th1f[b+"h_grandDaughterOfAntiS_gen_lxy"]->Fill(lxy(refPos, xyz));
		 histos_th1f[b+"h_grandDaughterOfAntiS_gen_dxy"]->Fill(dxy_signed_line_point(xyz, p, refPos)); 
		 histos_th1f[b+"h_grandDaughterOfAntiS_gen_dz"]->Fill(dz); 
		 histos_th1f[b+"h_grandDaughterOfAntiS_gen_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz()); 
		 histos_th2f[b+"h_grandDaughterOfAntiS_gen_vz_lxy"]->Fill(lxy(refPos, xyz), h_genParticles_SIM_GEANT->at(i).vz()); 

		 histos_teff[b+"heff_gen_particle_with_grandMother_antiS_pt"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).pt());
		 histos_teff[b+"heff_gen_particle_with_grandMother_antiS_vz"]->Fill(trackMatched,h_genParticles_SIM_GEANT->at(i).vz());
		 
		 	
	}
	else if(!genParticleGrandDaughterOfAntiS ){
		histos_th1f[b+"h_noGrandDaughterOfAntiS_gen_vz"]->Fill(h_genParticles_SIM_GEANT->at(i).vz());
	}
   }//end loop over h_genParticles_SIM_GEANT
   }//end if(h_genParticles_SIM_GEANT.isValid())











   //check if there is any antiS granddaughter which can be linked to a track
   for(unsigned int i = 0; i < h_genParticles_SIM_GEANT->size(); ++i){
	if(h_genParticles_SIM_GEANT->at(i).pdgId() == -1020000020){
		for(unsigned int j = 0; j < h_genParticles_SIM_GEANT->at(i).numberOfDaughters(); j++){
			for(unsigned int k = 0; k < h_genParticles_SIM_GEANT->at(i).daughter(j)->numberOfDaughters(); k++){
				const reco::Candidate * antiSGrandDaughter = h_genParticles_SIM_GEANT->at(i).daughter(j)->daughter(k);
				if(antiSGrandDaughter->charge() != 0 && fabs(antiSGrandDaughter->eta()) < 2.5){
					Double_t phi_antiSGrandDaughter = antiSGrandDaughter->phi();
					Double_t eta_antiSGrandDaughter = antiSGrandDaughter->eta();
					for (reco::TrackCollection::const_iterator iTk = theTrackCollection->begin(); iTk != theTrackCollection->end(); ++iTk) {
						const reco::Track* tmpTrack = &(*iTk);
						Double_t phiTrack = tmpTrack->phi();
						Double_t etaTrack = tmpTrack->eta();
						//TODO
						Double_t deltaR_track_grandDaughter = pow(pow( reco::deltaPhi(phi_antiSGrandDaughter,phiTrack),2)+pow(eta_antiSGrandDaughter-etaTrack,2),0.5);
						histos_th1f[b+"h_deltaR_V0_daughter_GEN_RECO_grandDaughterS"]->Fill(deltaR_track_grandDaughter);					
					}//end loop over track collection
				}//end if granddaughter is charged and in tracker accepatance
			}//end loop over antiS granddaughters		
		}//end loop over antiS daughters
	}//end if antiS	
   }//end loop over h_genParticles_SIM_GEANT


   // fill vectors of TransientTracks and TrackRefs after applying preselection cuts
   for (reco::TrackCollection::const_iterator iTk = theTrackCollection->begin(); iTk != theTrackCollection->end(); ++iTk) {
      const reco::Track* tmpTrack = &(*iTk);
      double ipsigXY = std::abs(tmpTrack->dxy(*theBeamSpot)/tmpTrack->dxyError());
      if (useVertex_) ipsigXY = std::abs(tmpTrack->dxy(referencePos)/tmpTrack->dxyError());
      double ipsigZ = std::abs(tmpTrack->dz(referencePos)/tmpTrack->dzError());

      //check if the track here matches a track which is a granddaughter of the antiS
      bool trackIsGrandDaughterOfAntiS = false;
      int nReconstructableGrandDaughtersMatchedToTrack = 0;
      if(antiSWithReconstructableGrandDaughters){
	      for(unsigned int i_V0_daughter = 0; i_V0_daughter < V0_daughters_phi_eta.size(); i_V0_daughter++){
		Double_t phi_V0_daughter = V0_daughters_phi_eta[i_V0_daughter][0];
		Double_t eta_V0_daughter = V0_daughters_phi_eta[i_V0_daughter][1];
		//TODO
		Double_t deltaR_V0_daughter_GEN_RECO = pow(pow(reco::deltaPhi(tmpTrack->phi(),phi_V0_daughter),2)+pow(tmpTrack->eta()-eta_V0_daughter,2),0.5);
		histos_th1f[b+"h_deltaR_V0_daughter_GEN_RECO"]->Fill(deltaR_V0_daughter_GEN_RECO);
		if(deltaR_V0_daughter_GEN_RECO < 0.1){
			nReconstructableGrandDaughtersMatchedToTrack++;
			trackIsGrandDaughterOfAntiS = true;
		}
		
	      }
      } 
      histos_th1f[b+"h_ReconstructableGrandDaughtersMatchedToTrack"]->Fill(nReconstructableGrandDaughtersMatchedToTrack);
      //use the switch to only study the grandaughter of antiS particles and check if the track could be matched to an antiS granddaighter
      if(specialGate(studyOnlySGrandDaughterTracks,trackIsGrandDaughterOfAntiS)){      
	      histos_th1f[b+"h_track_normalizedChi2"]->Fill(tmpTrack->normalizedChi2());
	      histos_th1f[b+"h_track_numberOfValidHits"]->Fill(tmpTrack->numberOfValidHits());
	      histos_th1f[b+"h_track_pt"]->Fill(tmpTrack->pt());
	      histos_th1f[b+"h_ipsigXY"]->Fill(ipsigXY);
	      histos_th1f[b+"h_ipsigZ"]->Fill(ipsigZ);
      }
      
      if (tmpTrack->normalizedChi2() < tkChi2Cut_ && tmpTrack->numberOfValidHits() >= tkNHitsCut_ &&
          tmpTrack->pt() > tkPtCut_ && ipsigXY > tkIPSigXYCut_ && ipsigZ > tkIPSigZCut_) {
         reco::TrackRef tmpRef(theTrackHandle, std::distance(theTrackCollection->begin(), iTk));
         theTrackRefs.push_back(std::move(tmpRef));
         reco::TransientTrack tmpTransient(*tmpRef, theMagneticField);
         theTransTracks.push_back(std::move(tmpTransient));
	 if(specialGate(studyOnlySGrandDaughterTracks,trackIsGrandDaughterOfAntiS)) histos_th1f[b+"h_goodTrackSelection"]->Fill(1);	 
      }
      else{
      	if(specialGate(studyOnlySGrandDaughterTracks,trackIsGrandDaughterOfAntiS)) histos_th1f[b+"h_goodTrackSelection"]->Fill(0);
      }
   }
   // good tracks have now been selected for vertexing

   // loop over tracks and vertex good charged track pairs
   for (unsigned int trdx1 = 0; trdx1 < theTrackRefs.size(); ++trdx1) {
      Double_t phi_track1 = theTrackRefs[trdx1]->phi();
      Double_t eta_track1 = theTrackRefs[trdx1]->eta();
      //std::cout << "track1 phi, eta: " << phi_track1 << " " << eta_track1 << std::endl;
      int foundKsDaugMatchingStep1 = 0;
      int foundLDaugMatchingStep1 = 0;
      //now check if this track matches a daughter of a V0
      if(antiSWithReconstructableGrandDaughters){
	      for(unsigned int i_V0_daughter = 0; i_V0_daughter < V0_daughters_phi_eta.size(); i_V0_daughter++){
		Double_t phi_V0_daughter = V0_daughters_phi_eta[i_V0_daughter][0];
		Double_t eta_V0_daughter = V0_daughters_phi_eta[i_V0_daughter][1];
		//TODO
		Double_t deltaR_V0_daughter_GEN_RECO = pow(pow( reco::deltaPhi(phi_track1,phi_V0_daughter),2)+pow(eta_track1-eta_V0_daughter,2),0.5);
		if(deltaR_V0_daughter_GEN_RECO < 0.1){
			if(i_V0_daughter == 0 || i_V0_daughter == 1)  foundKsDaugMatchingStep1++;
			if(i_V0_daughter == 2 || i_V0_daughter == 3)  foundLDaugMatchingStep1++;
		}
		
	      }
      } 

   for (unsigned int trdx2 = trdx1 + 1; trdx2 < theTrackRefs.size(); ++trdx2) {

      reco::TrackRef positiveTrackRef;
      reco::TrackRef negativeTrackRef;
      reco::TransientTrack* posTransTkPtr = nullptr;
      reco::TransientTrack* negTransTkPtr = nullptr;

      Double_t phi_track2 = theTrackRefs[trdx2]->phi();
      Double_t eta_track2 = theTrackRefs[trdx2]->eta();
      //std::cout << "track2 phi, eta: " << phi_track2 << " " << eta_track2 << std::endl;
      //now check if this track matches a daughter of a V0
      int foundKsDaugMatchingStep2 = foundKsDaugMatchingStep1;
      int foundLDaugMatchingStep2 = foundLDaugMatchingStep1;
      if(antiSWithReconstructableGrandDaughters){
	      for(unsigned int i_V0_daughter = 0; i_V0_daughter < V0_daughters_phi_eta.size(); i_V0_daughter++){
		Double_t phi_V0_daughter = V0_daughters_phi_eta[i_V0_daughter][0];
		Double_t eta_V0_daughter = V0_daughters_phi_eta[i_V0_daughter][1];
		//TODO
		Double_t deltaR_V0_daughter_GEN_RECO = pow(pow( reco::deltaPhi(phi_track2,phi_V0_daughter),2)+pow(eta_track2-eta_V0_daughter,2),0.5);
		if(deltaR_V0_daughter_GEN_RECO < 0.1){
			if(i_V0_daughter == 0 || i_V0_daughter == 1)  foundKsDaugMatchingStep2++;
			if(i_V0_daughter == 2 || i_V0_daughter == 3)  foundLDaugMatchingStep2++;
		}
		
	      }
      }
      if(foundKsDaugMatchingStep2==2)std::cout << "found 2 daughters of the GEN Ks from the antiS matching with a track " << std::endl;
      if(foundLDaugMatchingStep2==2)std::cout << "found 2 daughters of the GEN L from the antiS matching with a track " << std::endl;

      bool found2DaugMatching = false;
      if(foundKsDaugMatchingStep2==2 || foundLDaugMatchingStep2 == 2)found2DaugMatching = true;

      if (theTrackRefs[trdx1]->charge() < 0. && theTrackRefs[trdx2]->charge() > 0.) {
         negativeTrackRef = theTrackRefs[trdx1];
         positiveTrackRef = theTrackRefs[trdx2];
         negTransTkPtr = &theTransTracks[trdx1];
         posTransTkPtr = &theTransTracks[trdx2];
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_oppositeCharge"]->Fill(1);	 
      } else if (theTrackRefs[trdx1]->charge() > 0. && theTrackRefs[trdx2]->charge() < 0.) {
         negativeTrackRef = theTrackRefs[trdx2];
         positiveTrackRef = theTrackRefs[trdx1];
         negTransTkPtr = &theTransTracks[trdx2];
         posTransTkPtr = &theTransTracks[trdx1];
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_oppositeCharge"]->Fill(1);	 
      } else {
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_oppositeCharge"]->Fill(0);	 
         continue;
      }
      if(found2DaugMatching) std::cout << "tracks have opposite charge" << std::endl;
      // measure distance between tracks at their closest approach
      if (!posTransTkPtr->impactPointTSCP().isValid() || !negTransTkPtr->impactPointTSCP().isValid()){
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_impactPointTSCP_valid"]->Fill(0);	 
	 continue;
      }
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_impactPointTSCP_valid"]->Fill(1);	 

      FreeTrajectoryState const & posState = posTransTkPtr->impactPointTSCP().theState();
      FreeTrajectoryState const & negState = negTransTkPtr->impactPointTSCP().theState();
      ClosestApproachInRPhi cApp;
      cApp.calculate(posState, negState);
      if (!cApp.status()){
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_cApp_status"]->Fill(0);
	 continue;
      }
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_cApp_status"]->Fill(1);
      float dca = std::abs(cApp.distance());
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_dca"]->Fill(dca);
      if (dca > tkDCACut_) continue;

      if(found2DaugMatching) std::cout << "tracks PCA cut" << std::endl;
      // the POCA should at least be in the sensitive volume
      GlobalPoint cxPt = cApp.crossingPoint();
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_POCA_sensVolume_xy"]->Fill(sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()));
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_POCA_sensVolume_z"]->Fill(std::abs(cxPt.z()));
      if (sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.) continue;
      if(found2DaugMatching) std::cout << "tracks PCA sensitive volume cut" << std::endl;

      // the tracks should at least point in the same quadrant
      TrajectoryStateClosestToPoint const & posTSCP = posTransTkPtr->trajectoryStateClosestToPoint(cxPt);
      TrajectoryStateClosestToPoint const & negTSCP = negTransTkPtr->trajectoryStateClosestToPoint(cxPt);
      if (!posTSCP.isValid() || !negTSCP.isValid()){ 
	if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_TSCP_valid"]->Fill(0);
	continue;
      }
      histos_th1f[b+"h_TSCP_valid"]->Fill(1);	
      if (posTSCP.momentum().dot(negTSCP.momentum())  < 0){
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_TSCP_directions"]->Fill(0);	
	 continue;
      }
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_TSCP_directions"]->Fill(1);	
     

      if(found2DaugMatching) std::cout << "quadrant cut" << std::endl;
     
      // calculate mPiPi
      double totalE = sqrt(posTSCP.momentum().mag2() + piMassSquared) + sqrt(negTSCP.momentum().mag2() + piMassSquared);
      double totalESq = totalE*totalE;
      double totalPSq = (posTSCP.momentum() + negTSCP.momentum()).mag2();
      double mass = sqrt(totalESq - totalPSq);
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_massPiPi"]->Fill(mass);
      if (mass > mPiPiCut_) continue;
      if(found2DaugMatching) std::cout << "quadrant cut" << std::endl;

      // Fill the vector of TransientTracks to send to KVF
      std::vector<reco::TransientTrack> transTracks;
      transTracks.reserve(2);
      transTracks.push_back(*posTransTkPtr);
      transTracks.push_back(*negTransTkPtr);

      // create the vertex fitter object and vertex the tracks
      TransientVertex theRecoVertex;
      if (vertexFitter_) {
         KalmanVertexFitter theKalmanFitter(useRefTracks_ == 0 ? false : true);
         theRecoVertex = theKalmanFitter.vertex(transTracks);
      } else if (!vertexFitter_) {
         useRefTracks_ = false;
         AdaptiveVertexFitter theAdaptiveFitter;
         theRecoVertex = theAdaptiveFitter.vertex(transTracks);
      }
      if (!theRecoVertex.isValid()){
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_theRecoVertexValid"]->Fill(0);
	 continue;
      }	
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_theRecoVertexValid"]->Fill(1);
      if(found2DaugMatching) std::cout << "vertex valid cut" << std::endl;
     
      reco::Vertex theVtx = theRecoVertex;
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_vtxNormalizedChi2"]->Fill(theVtx.normalizedChi2()); 
      if (theVtx.normalizedChi2() > vtxChi2Cut_) continue;
      if(found2DaugMatching) std::cout << "chi2 cut" << std::endl;
      GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());

      // 2D decay significance
      SMatrixSym3D totalCov = theBeamSpot->rotatedCovariance3D() + theVtx.covariance();
      if (useVertex_) totalCov = referenceVtx.covariance() + theVtx.covariance();
      SVector3 distVecXY(vtxPos.x()-referencePos.x(), vtxPos.y()-referencePos.y(), 0.);
      double distMagXY = ROOT::Math::Mag(distVecXY);
      double sigmaDistMagXY = sqrt(ROOT::Math::Similarity(totalCov, distVecXY)) / distMagXY;
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_distMagXY_over_sigmaDistMagXY"]->Fill(distMagXY/sigmaDistMagXY);
      if (distMagXY/sigmaDistMagXY < vtxDecaySigXYCut_) continue;
      //if(found2DaugMatching) std::cout << "2D significance cut" << std::endl;
      std::cout << "2D significance cut" << std::endl;

      // 3D decay significance
      SVector3 distVecXYZ(vtxPos.x()-referencePos.x(), vtxPos.y()-referencePos.y(), vtxPos.z()-referencePos.z());
      double distMagXYZ = ROOT::Math::Mag(distVecXYZ);
      double sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_distMagXYZ_over_sigmaDistMagXYZ"]->Fill(distMagXY/sigmaDistMagXY);
      if (distMagXYZ/sigmaDistMagXYZ < vtxDecaySigXYZCut_) continue;
//      if(found2DaugMatching) std::cout << "3D significance cut" << std::endl;
      std::cout << "3D significance cut" << std::endl;
     cout << "0" << endl;

      // make sure the vertex radius is within the inner track hit radius
      if (innerHitPosCut_ > 0. && positiveTrackRef->innerOk()) {
	 cout << "1" << endl;
         reco::Vertex::Point posTkHitPos = positiveTrackRef->innerPosition();
	  cout << "2" << endl;
         double posTkHitPosD2 =  (posTkHitPos.x()-referencePos.x())*(posTkHitPos.x()-referencePos.x()) +
            (posTkHitPos.y()-referencePos.y())*(posTkHitPos.y()-referencePos.y());
	 cout << "3" << endl;
         if (sqrt(posTkHitPosD2) < (distMagXY - sigmaDistMagXY*innerHitPosCut_)){
	 	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_vertexRadius_within_inner_track1"]->Fill(0);
		 continue;
	 }
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_vertexRadius_within_inner_track1"]->Fill(1);
         cout << "4" << endl;
      }
      //if(found2DaugMatching) std::cout << "inner track radius cut1" << std::endl;
      std::cout << "inner track radius cut1" << std::endl;
      if (innerHitPosCut_ > 0. && negativeTrackRef->innerOk()) {
         reco::Vertex::Point negTkHitPos = negativeTrackRef->innerPosition();
         double negTkHitPosD2 = (negTkHitPos.x()-referencePos.x())*(negTkHitPos.x()-referencePos.x()) +
            (negTkHitPos.y()-referencePos.y())*(negTkHitPos.y()-referencePos.y());
         if (sqrt(negTkHitPosD2) < (distMagXY - sigmaDistMagXY*innerHitPosCut_)){
	 	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_vertexRadius_within_inner_track2"]->Fill(0);
		 continue;
	 }
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_vertexRadius_within_inner_track2"]->Fill(1);
      }
      //if(found2DaugMatching) std::cout << "inner track radius cut2" << std::endl;
      std::cout << "inner track radius cut2" << std::endl;
      
      std::auto_ptr<TrajectoryStateClosestToPoint> trajPlus;
      std::auto_ptr<TrajectoryStateClosestToPoint> trajMins;
      std::vector<reco::TransientTrack> theRefTracks;
      if (theRecoVertex.hasRefittedTracks()) {
         theRefTracks = theRecoVertex.refittedTracks();
      }

      if (useRefTracks_ && theRefTracks.size() > 1) {
         reco::TransientTrack* thePositiveRefTrack = 0;
         reco::TransientTrack* theNegativeRefTrack = 0;
         for (std::vector<reco::TransientTrack>::iterator iTrack = theRefTracks.begin(); iTrack != theRefTracks.end(); ++iTrack) {
            if (iTrack->track().charge() > 0.) {
               thePositiveRefTrack = &*iTrack;
            } else if (iTrack->track().charge() < 0.) {
               theNegativeRefTrack = &*iTrack;
            }
         }
         if (thePositiveRefTrack == 0 || theNegativeRefTrack == 0){
	 	if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_RefTracksValid"]->Fill(0);
		 continue;
	 }
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_RefTracksValid"]->Fill(1);
         trajPlus.reset(new TrajectoryStateClosestToPoint(thePositiveRefTrack->trajectoryStateClosestToPoint(vtxPos)));
         trajMins.reset(new TrajectoryStateClosestToPoint(theNegativeRefTrack->trajectoryStateClosestToPoint(vtxPos)));
      } else {
         trajPlus.reset(new TrajectoryStateClosestToPoint(posTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
         trajMins.reset(new TrajectoryStateClosestToPoint(negTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
      }

      if (trajPlus.get() == 0 || trajMins.get() == 0 || !trajPlus->isValid() || !trajMins->isValid()){
	 if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_trajPlus_trajMins_valid"]->Fill(0);
	 continue;
      }
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_trajPlus_trajMins_valid"]->Fill(1);
 

      GlobalVector positiveP(trajPlus->momentum());
      GlobalVector negativeP(trajMins->momentum());
      GlobalVector totalP(positiveP + negativeP);

      // 2D pointing angle
      double dx = theVtx.x()-referencePos.x();
      double dy = theVtx.y()-referencePos.y();
      double px = totalP.x();
      double py = totalP.y();
      double angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_angleXY"]->Fill(angleXY);
      if (angleXY < cosThetaXYCut_) continue;
      //if(found2DaugMatching) std::cout << "2D pointing angle cut" << std::endl;
      std::cout << "2D pointing angle cut" << std::endl;

      // 3D pointing angle
      double dz = theVtx.z()-referencePos.z();
      double pz = totalP.z();
      double angleXYZ = (dx*px+dy*py+dz*pz)/(sqrt(dx*dx+dy*dy+dz*dz)*sqrt(px*px+py*py+pz*pz));
      if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_angleXYZ"]->Fill(angleXYZ);
      if (angleXYZ < cosThetaXYZCut_) continue;
      //if(found2DaugMatching) std::cout << "3D pointing angle cut" << std::endl;
      std::cout << "3D pointing angle cut" << std::endl;

      // calculate total energy of V0 3 ways: assume it's a kShort, a Lambda, or a LambdaBar.
      double piPlusE = sqrt(positiveP.mag2() + piMassSquared);
      double piMinusE = sqrt(negativeP.mag2() + piMassSquared);
      double protonE = sqrt(positiveP.mag2() + protonMassSquared);
      double antiProtonE = sqrt(negativeP.mag2() + protonMassSquared);
      double kShortETot = piPlusE + piMinusE;
      double lambdaEtot = protonE + piMinusE;
      double lambdaBarEtot = antiProtonE + piPlusE;

      // Create momentum 4-vectors for the 3 candidate types
      const reco::Particle::LorentzVector kShortP4(totalP.x(), totalP.y(), totalP.z(), kShortETot);
      const reco::Particle::LorentzVector lambdaP4(totalP.x(), totalP.y(), totalP.z(), lambdaEtot);
      const reco::Particle::LorentzVector lambdaBarP4(totalP.x(), totalP.y(), totalP.z(), lambdaBarEtot);

      reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
      const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
      double vtxChi2(theVtx.chi2());
      double vtxNdof(theVtx.ndof());

      // Create the VertexCompositeCandidate object that will be stored in the Event
      reco::VertexCompositeCandidate* theKshort = nullptr;
      reco::VertexCompositeCandidate* theLambda = nullptr;
      reco::VertexCompositeCandidate* theLambdaBar = nullptr;

      if (doKShorts_) {
         theKshort = new reco::VertexCompositeCandidate(0, kShortP4, vtx, vtxCov, vtxChi2, vtxNdof);
      }
      if (doLambdas_) {
         if (positiveP.mag2() > negativeP.mag2()) {
            theLambda = new reco::VertexCompositeCandidate(0, lambdaP4, vtx, vtxCov, vtxChi2, vtxNdof);
         } else {
            theLambdaBar = new reco::VertexCompositeCandidate(0, lambdaBarP4, vtx, vtxCov, vtxChi2, vtxNdof);
         }
      }

      // Create daughter candidates for the VertexCompositeCandidates
      reco::RecoChargedCandidate thePiPlusCand(
         1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), piPlusE), vtx);
      thePiPlusCand.setTrack(positiveTrackRef);
      
      reco::RecoChargedCandidate thePiMinusCand(
         -1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), piMinusE), vtx);
      thePiMinusCand.setTrack(negativeTrackRef);
      
      reco::RecoChargedCandidate theProtonCand(
         1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), protonE), vtx);
      theProtonCand.setTrack(positiveTrackRef);

      reco::RecoChargedCandidate theAntiProtonCand(
         -1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), antiProtonE), vtx);
      theAntiProtonCand.setTrack(negativeTrackRef);

      AddFourMomenta addp4;
      // Store the daughter Candidates in the VertexCompositeCandidates if they pass mass cuts
      if (doKShorts_) {
         theKshort->addDaughter(thePiPlusCand);
         theKshort->addDaughter(thePiMinusCand);
         theKshort->setPdgId(310);
         addp4.set(*theKshort);
         if (theKshort->mass() < kShortMass + kShortMassCut_ && theKshort->mass() > kShortMass - kShortMassCut_) {
            //theKshorts.push_back(std::move(*theKshort));
      	    if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_V0_formed"]->Fill(1);
	    if(foundKsDaugMatchingStep2==2){
		std::cout << "put a Ks in the event which is linked to an antiS" << std::endl;
		std::cout << "daug0 of the Lambda phi, eta: "<< V0_daughters_phi_eta[0][0] << " " << V0_daughters_phi_eta[0][1] << std::endl;
		std::cout << "daug1 of the Lambda phi, eta: "<< V0_daughters_phi_eta[1][0] << " " << V0_daughters_phi_eta[1][1] << std::endl;
		std::cout << "daug0 of the Lambda phi, eta: "<< V0_daughters_phi_eta[2][0] << " " << V0_daughters_phi_eta[2][1] << std::endl;
		std::cout << "daug1 of the Lambda phi, eta: "<< V0_daughters_phi_eta[3][0] << " " << V0_daughters_phi_eta[3][1] << std::endl;
      		std::cout << "track1 phi, eta: " << phi_track1 << " " << eta_track1 << std::endl;
      		std::cout << "track2 phi, eta: " << phi_track2 << " " << eta_track2 << std::endl;
	    }
         }
      }
      if (doLambdas_ && theLambda) {
         theLambda->addDaughter(theProtonCand);
         theLambda->addDaughter(thePiMinusCand);
         theLambda->setPdgId(3122);
         addp4.set( *theLambda );
         if (theLambda->mass() < lambdaMass + lambdaMassCut_ && theLambda->mass() > lambdaMass - lambdaMassCut_) {
            //theLambdas.push_back(std::move(*theLambda));
      	    if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_V0_formed"]->Fill(2);
	    if(foundLDaugMatchingStep2==2){
		std::cout << "put a Lambda in the event which is linked to an antiS" << std::endl;
		std::cout << "daug0 of the Lambda phi, eta: "<< V0_daughters_phi_eta[0][0] << " " << V0_daughters_phi_eta[0][1] << std::endl;
		std::cout << "daug1 of the Lambda phi, eta: "<< V0_daughters_phi_eta[1][0] << " " << V0_daughters_phi_eta[1][1] << std::endl;
		std::cout << "daug0 of the Lambda phi, eta: "<< V0_daughters_phi_eta[2][0] << " " << V0_daughters_phi_eta[2][1] << std::endl;
		std::cout << "daug1 of the Lambda phi, eta: "<< V0_daughters_phi_eta[3][0] << " " << V0_daughters_phi_eta[3][1] << std::endl;
      		std::cout << "track1 phi, eta: " << phi_track1 << " " << eta_track1 << std::endl;
      		std::cout << "track2 phi, eta: " << phi_track2 << " " << eta_track2 << std::endl;
	    }
         }
      } else if (doLambdas_ && theLambdaBar) {
         theLambdaBar->addDaughter(theAntiProtonCand);
         theLambdaBar->addDaughter(thePiPlusCand);
         theLambdaBar->setPdgId(-3122);
         addp4.set(*theLambdaBar);
         if (theLambdaBar->mass() < lambdaMass + lambdaMassCut_ && theLambdaBar->mass() > lambdaMass - lambdaMassCut_) {
            //theLambdas.push_back(std::move(*theLambdaBar));
    	    if(specialGate(studyOnlySGrandDaughterTracks,found2DaugMatching))histos_th1f[b+"h_V0_formed"]->Fill(3);
	    if(foundLDaugMatchingStep2==2){
		std::cout << "put a Lambdabar in the event which is linked to an antiS" << std::endl;
		std::cout << "daug0 of the Lambda phi, eta: "<< V0_daughters_phi_eta[0][0] << " " << V0_daughters_phi_eta[0][1] << std::endl;
		std::cout << "daug1 of the Lambda phi, eta: "<< V0_daughters_phi_eta[1][0] << " " << V0_daughters_phi_eta[1][1] << std::endl;
		std::cout << "daug0 of the Lambda phi, eta: "<< V0_daughters_phi_eta[2][0] << " " << V0_daughters_phi_eta[2][1] << std::endl;
		std::cout << "daug1 of the Lambda phi, eta: "<< V0_daughters_phi_eta[3][0] << " " << V0_daughters_phi_eta[3][1] << std::endl;
      		std::cout << "track1 phi, eta: " << phi_track1 << " " << eta_track1 << std::endl;
      		std::cout << "track2 phi, eta: " << phi_track2 << " " << eta_track2 << std::endl;
	    }
         }
      }

      delete theKshort;
      delete theLambda;
      delete theLambdaBar;
      theKshort = theLambda = theLambdaBar = nullptr;

    }
  } 

} //end of analyzer

double Analyzer_V0Fitter::openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2){
  double opening_angle = TMath::ACos((momentum1.Dot(momentum2))/(pow(momentum1.Mag2()*momentum2.Mag2(),0.5)));
  return opening_angle;
}

double Analyzer_V0Fitter::deltaR(double phi1, double eta1, double phi2, double eta2){
	double deltaPhi = reco::deltaPhi(phi1,phi2);
	double deltaEta = eta1-eta2;
	return pow(deltaPhi*deltaPhi+deltaEta*deltaEta,0.5);
}


double Analyzer_V0Fitter::lxy(TVector3 v1, TVector3 v2){
	double x1 = v1.X();
	double x2 = v2.X();
	double y1 = v1.Y();
	double y2 = v2.Y();
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
}


TVector3 Analyzer_V0Fitter::PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point){
   //first move the vector along the line to the starting point of Point_line
   double normalise = sqrt(Vector_along_line.X()*Vector_along_line.X()+Vector_along_line.Y()*Vector_along_line.Y()+Vector_along_line.Z()*Vector_along_line.Z());
   TVector3 n(Vector_along_line.X()/normalise,Vector_along_line.Y()/normalise,Vector_along_line.Z()/normalise);
   TVector3 a = Point_line;
   TVector3 p = Point;

   //see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line (Vector formulation)
   TVector3 vector_PCA = (a-p)-((a-p)*n)*n;
   return vector_PCA ;
}

double Analyzer_V0Fitter::dxy_signed_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in){

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

double Analyzer_V0Fitter::phiInTrackRefPosition(int charge, double p1x, double p1y, double x1, double y1, double xREF, double yREF){

    double pt = pow(p1x*p1x+p1y*p1y,0.5);
    cout << "pt "<<  pt << endl;
    cout << "charge "<<  charge << endl;
    double r = pt*10/(3.8);
    cout << "r " << r << endl;
    


    double phi1 =  TMath::ACos(p1x/sqrt(p1x*p1x+p1y*p1y));
    if(p1y < 0){
        phi1 =  - phi1;
    }
    cout << "phi1 "<< phi1 << endl;
    cout << "x1 "<< x1 << endl;
    cout << "y1 "<< y1 << endl;

    //check if particle points towards zero zero or away
    double inproductCoordMom = (x1-xREF)*p1x+(y1-yREF)*p1y; 
    bool partilcePointsToZeroZero = false;
    if(inproductCoordMom < 0){
        partilcePointsToZeroZero = true;
    }	

    double x0 = 0;
    double y0 = 0;

    double x0_small = x1 - sqrt( pow(TMath::Tan(phi1),2) *r*r /  (1+(pow(TMath::Tan(phi1),2))) );
    double x0_large = x1 + sqrt( pow(TMath::Tan(phi1),2) *r*r /  (1+(pow(TMath::Tan(phi1),2))) );


    //with the minus signs
    double y0_small = y1 - r*pow(1-pow(TMath::Tan(phi1),2)/(1+pow(TMath::Tan(phi1),2)),0.5);
    //with the pos signs
    double y0_large = y1 + r*pow(1-pow(TMath::Tan(phi1),2)/(1+pow(TMath::Tan(phi1),2)),0.5);

    cout << "partilcePointsToZeroZero " << partilcePointsToZeroZero << endl;

    //pos charged particles:
    if(charge == 1){
        //1st quarter
        if(phi1 >= 0 && phi1 < TMath::Pi()/2){
            x0 = x0_large;
            y0 = y0_small;
	}
        //2nd quarter
        else if(phi1 >=  TMath::Pi()/2 && phi1 <  TMath::Pi()){
            x0 = x0_large;
            y0 = y0_large;
	}
        //3rd quarter
        else if(phi1 > - TMath::Pi() && phi1 <= - TMath::Pi()/2){
            x0 = x0_small;
            y0 = y0_large;
	}
        //4th quarter
        else if(phi1 < 0 && phi1 >= - TMath::Pi()/2){
            x0 = x0_small;
            y0 = y0_small;
	}
    }
    //pos charged particles:
    if(charge == -1){
        //1st quarter
        if(phi1 >= 0 && phi1 <  TMath::Pi()/2){
            x0 = x0_small;
            y0 = y0_large;
	}
        //2nd quarter
        else if(phi1 >=  TMath::Pi()/2 && phi1 <  TMath::Pi()){
            x0 = x0_small;
            y0 = y0_small;
	}
        //3rd quarter
        else if(phi1 > - TMath::Pi() && phi1 <= - TMath::Pi()/2){
            x0 = x0_large;
            y0 = y0_small;
	}
        //4th quarter
        else if(phi1 < 0 && phi1 >= - TMath::Pi()/2){
            x0 = x0_large;
            y0 = y0_large;
	}
    }	
/*
    //look ath the 1st quarter
    if(x1>=0 && y1>=0){
        if((!partilcePointsToZeroZero && charge < 0) || (partilcePointsToZeroZero && charge > 0)){
            x0 = x0_small;
            y0 = y0_large;
	}
        else{
            x0 = x0_large;
            y0 = y0_small;
	}
    }		
    //look ath the 2nd quarter
    else if(x1<=0 && y1>=0){
        if((!partilcePointsToZeroZero && charge < 0) || (partilcePointsToZeroZero && charge > 0)){
            x0 = x0_small;
            y0 = y0_small;
	}
        else{
            x0 = x0_large;
            y0 = y0_large;
	}
    }
    //look ath the 3rd quarter
    else if(x1<=0 && y1<=0){
        if((!partilcePointsToZeroZero && charge < 0) || (partilcePointsToZeroZero && charge > 0)){
            x0 = x0_large;
            y0 = y0_small;
	}
        else{
            x0 = x0_small;
            y0 = y0_large;
	}
    }
    //look ath the 4th quarter
    else{
        if((!partilcePointsToZeroZero && charge < 0) || (partilcePointsToZeroZero && charge > 0)){
            x0 = x0_large;
            y0 = y0_large;
	}
        else{
            x0 = x0_small;
            y0 = y0_small;
	}
    }
*/
    cout << "x0_small "<< x0_small << endl;
    cout << "x0_large "<< x0_large << endl;
    cout << "y0_small "<< y0_small << endl;
    cout << "y0_large "<< y0_large << endl;

    cout << "x0 " << x0 << endl;
    cout << "y0 " << y0 << endl;


     //calculate where the circle (center x0,y0 and radius r) and the line connecting the REF point and x0,y0 cross
    double x2_1 = r*pow(1+pow(yREF-y1,2)/pow(xREF-x0,2),-0.5)+x0;
    double x2_2 = -r*pow(1+pow(yREF-y1,2)/pow(xREF-x0,2),-0.5)+x0;
    
    double y2_1 = y0 + (yREF-y0)*(x2_1-x0)/(xREF-x0);
    double y2_2 = y0 + (yREF-y0)*(x2_2-x0)/(xREF-x0);

    cout << "x2_1, y2_1 " << x2_1 << ", " << y2_1 << endl;
    cout << "x2_2, y2_2 " << x2_2 << ", " << y2_2 << endl;
    cout << "as a check the difference between this should be 2*r: " << sqrt(pow(x2_1-x2_2,2)+pow(y2_1-y2_2,2)) << endl;

    //now calculate which one of the point2 is closest to the REF point
    double distREFx2_1_y2_1 = sqrt(pow(xREF-x2_1,2)+pow(yREF-y2_1,2));
    double distREFx2_2_y2_2 = sqrt(pow(xREF-x2_2,2)+pow(yREF-y2_2,2));

    double distREFX2Y2 = distREFx2_1_y2_1;
    if(distREFx2_2_y2_2<distREFx2_1_y2_1){
	distREFX2Y2 = distREFx2_2_y2_2;
    }

    double distREFX1Y1 = sqrt(pow(x1-xREF,2)+pow(y1-yREF,2));

    //if the dist between REF and point1 (the GEN ref point) is smaller than the distance between the REF and the point2 (the track ref point calculated) then you can just use the phi1 for the phi2
    double phi2 = phi1;	
    cout << "distREFX1Y1 " << distREFX1Y1 << endl;
    cout << "distREFX2Y2 " << distREFX2Y2 << endl;
    //if however the dist between REF and point1 is bigger than the distance between REF and point2, then you need to use the one calculated below
    if(distREFX1Y1 > distREFX2Y2){

	    double phi0 = TMath::ACos((x0-xREF)/sqrt((x0-xREF)*(x0-xREF)+(y0-yREF)*(y0-yREF)));

	    if(y0 < 0){
		phi0 = -phi0;
	    }
	    cout << "phi0" <<  phi0 << endl;

	    if(charge == 1){
		phi2 = phi0 + TMath::Pi()/2;
	    }
	    else if(charge == -1){
		phi2 = phi0 - TMath::Pi()/2;
	    }

	    if(phi2 > TMath::Pi()){
		phi2 = phi2 - 2*TMath::Pi();
	    }
	    if(phi2 < TMath::Pi()){
		phi2 = phi2 + 2*TMath::Pi();
	    }
    }

    return phi2;
}


bool Analyzer_V0Fitter::specialGate(bool flag1, bool flag2){
   bool returnBool = false;
   if(!flag1 && !flag2) returnBool =  true;
   if(!flag1 && flag2) returnBool = true;
   if(flag1 && !flag2) returnBool = false;
   if(flag1 && flag2) returnBool = true;
   return returnBool;
}

double Analyzer_V0Fitter::dz_line_point(TVector3 Point_line_in, TVector3 Vector_along_line_in, TVector3 Point_in){
  //looking at Z, so put the XY component to 0 first
//  TVector3 Point_line(0.,0., Point_line_in.Z());
//  TVector3 Vector_along_line(0.,0., Vector_along_line_in.Z());
//  TVector3 Point( 0., 0., Point_in.Z());

//  TVector3 shortest_distance = PCA_line_point(Point_line,  Vector_along_line, Point);
//  return shortest_distance.Z();
  Double_t vz = Point_line_in.Z();
  Double_t vx = Point_line_in.X();
  Double_t vy = Point_line_in.Y();
  Double_t px = Vector_along_line_in.X();
  Double_t py = Vector_along_line_in.Y();
  Double_t pz = Vector_along_line_in.Z();
  Double_t pt = sqrt(px*px+py*py);
  return (vz - Point_in.Z()) - ((vx - Point_in.X()) * px + (vy - Point_in.Y()) * py) / pt * pz / pt;

}

void Analyzer_V0Fitter::endJob()
{
}

Analyzer_V0Fitter::~Analyzer_V0Fitter()
{
}


DEFINE_FWK_MODULE(Analyzer_V0Fitter);
