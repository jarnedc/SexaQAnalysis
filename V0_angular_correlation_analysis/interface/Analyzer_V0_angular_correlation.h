#ifndef My_azimuthal_MC_h
#define My_azimuthal_MC_h
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TProfile.h"
#include <TMath.h>
#include <math.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TString.h>

#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
using namespace edm;
using namespace std;

#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//vectors to save data accross events
    static vector<double> v_L0_phi;
    static vector<double> v_L0_eta;
    static vector<double> v_Ks_phi;
    static vector<double> v_Ks_eta;


int signum(float x){
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

float RandomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

  
class Analyzer_V0_angular_correlation : public edm::EDAnalyzer
 {
  public:
    explicit Analyzer_V0_angular_correlation(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~Analyzer_V0_angular_correlation();
    TVector3 PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double dxy_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double dxy_signed_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double dxy_signed_line_point2(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double dz_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double lxy(TVector3 v1, TVector3 v2);
    double std_dev_lxy(double vx, double vy, double vx_var, double vy_var, double bx_x, double bx_y, double bx_x_var, double bx_y_var);
    double lxy_signed(TVector3 particle_vertex, TVector3 beamspot, TVector3 particle_direction);
  private:
    //---- configurable parameters --------
    bool m_isData;

    int m_nEvent, m_nRun, m_nLumi;
    
    edm::Service<TFileService> m_fs;
 
    edm::InputTag m_bsTag;
    edm::InputTag m_vertexTag;
    edm::InputTag m_genParticlesTag;
    //edm::InputTag m_rCandsTag;
    edm::InputTag m_sCandsTag;
    
    edm::InputTag m_KshortsTag;
    edm::InputTag m_LambdasTag;
   
    edm::InputTag m_LambdasLambdaKshortFilterTag;
    edm::InputTag m_KshortsLambdaKshortFilterTag;
    
    edm::InputTag m_sCollectionMassFilterTag;
    edm::InputTag m_rCollectionMassFilterTag;
 

    edm::InputTag m_nPVsTag;
    edm::InputTag m_nelectronsTag;
    edm::InputTag m_njetsTag;
    edm::InputTag m_nkshortsTag;
    edm::InputTag m_nlambdasTag;
    edm::InputTag m_nmuonsTag;
    edm::InputTag m_ntracksTag;
    edm::InputTag m_HTTag;
    edm::InputTag m_TKHTTag;
    edm::InputTag m_TwoTopJetsTag;
    edm::InputTag m_METTag;
    edm::InputTag m_TKMETTag;
 
    edm::EDGetTokenT<reco::BeamSpot> m_bsToken;
    edm::EDGetTokenT<vector<reco::Vertex> > m_vertexToken;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken;
    //edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate> > m_rCandsToken;
    //edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate> > m_sCandsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_sCandsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_KshortsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_LambdasToken;
    edm::EDGetTokenT<edm::PtrVector<reco::Candidate > > m_LambdasLambdaKshortFilterToken;
    edm::EDGetTokenT<edm::PtrVector<reco::Candidate > > m_KshortsLambdaKshortFilterToken;
    edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate> > m_sCollectionMassFilterToken;
    edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate> > m_rCollectionMassFilterToken;
    
    //from the initialproducer 
    edm::EDGetTokenT<vector<int> > m_nPVsToken;
    edm::EDGetTokenT<vector<int> > m_nelectronsToken;
    edm::EDGetTokenT<vector<int> > m_njetsToken;
    edm::EDGetTokenT<vector<int> > m_nkshortsToken;
    edm::EDGetTokenT<vector<int> > m_nlambdasToken;
    edm::EDGetTokenT<vector<int> > m_nmuonsToken;
    edm::EDGetTokenT<vector<int> > m_ntracksToken;
    edm::EDGetTokenT<vector<double> > m_HTToken;
    edm::EDGetTokenT<vector<double> > m_TKHTToken;
    edm::EDGetTokenT<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > m_TwoTopJetsToken;
    edm::EDGetTokenT<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > m_METToken;
    edm::EDGetTokenT<vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > > m_TKMETToken;

    
    


    int verbose=1;
    
    TString a = ""; //"WjetsMC" "ZeroBias" "MET" "MinBiasMC" "SingleMuon"
    TString b = a + "";
    
    //--------- Histogram Declaration --------------------//
    std::map<TString, TH1F *> histos_th1f;
    std::map<TString, TH2F *> histos_th2f;
    //std::map<TString, TH3F *> histos_th3f;
    std::map<TString, TProfile *> histos_TProfile;

     };

#endif

