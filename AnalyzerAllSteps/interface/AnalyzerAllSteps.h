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
#include "TEfficiency.h"
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


  
class AnalyzerAllSteps : public edm::EDAnalyzer
 {
  public:
    explicit AnalyzerAllSteps(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~AnalyzerAllSteps();
    double openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2);
    double deltaR(double phi1, double eta1, double phi2, double eta2);
    double lxy(TVector3 v1, TVector3 v2);
    TVector3 PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double dxy_signed_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double std_dev_lxy(double vx, double vy, double vx_var, double vy_var, double bx_x, double bx_y, double bx_x_var, double bx_y_var);

    void FillHistosGENAntiS(const reco::Candidate * , TVector3 beamspot);
    void FillHistosGENKsNonAntiS(const reco::Candidate *, TVector3 beamspot );
    void FillHistosGENAntiLambdaNonAntiS(const reco::Candidate *,TVector3 beamspot );
    void FillHistosGENKsAntiS(const reco::Candidate *, TVector3 beamspot );
    void FillHistosGENAntiLambdaAntiS(const reco::Candidate *,TVector3 beamspot );
    void FillHistoAllGranddaughtersCorrect(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L, int eventId);
    void FillHistosRECOKs(const reco::VertexCompositeCandidate * RECOKs, TVector3 beamspot);
    void FillHistosRECOAntiLambda(const reco::VertexCompositeCandidate * RECOAntiLambda, TVector3 beamspot);
    void FillHistosRECOAntiS(const reco::VertexCompositeCandidate * RECOAntiS, TVector3 beamspot,TVector3 beamspotVariance, edm::Handle<vector<reco::GenParticle>>  h_genParticles, int eventId);
    void RecoEvaluationKsNonAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks,  TVector3 beamspot);
    void RecoEvaluationKsAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0Ks,  TVector3 beamspot);
    void RecoEvaluationAntiLambdaNonAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L,  TVector3 beamspot);
    void RecoEvaluationAntiLambdaAntiS(const reco::Candidate  * genParticle, edm::Handle<vector<reco::VertexCompositeCandidate> > h_V0L,  TVector3 beamspot);

    const int pdgIdAntiS = -1020000020;
    const int pdgIdKs = 310;
    const int pdgIdAntiLambda = -3122;
    const int pdgIdAntiProton = -2212;
    const int pdgIdPosPion = 211;
    const int pdgIdNegPion = -211;


  private:
    //---- configurable parameters --------
    bool m_isData;
    int m_nEvent, m_nRun, m_nLumi;
    
    edm::Service<TFileService> m_fs;
 
    edm::InputTag m_bsTag;
    edm::InputTag m_offlinePVTag;
    edm::InputTag m_genParticlesTag_GEN;
    edm::InputTag m_genParticlesTag_SIM_GEANT;
    edm::InputTag m_generalTracksTag;
    edm::InputTag m_sCandsTag;
    edm::InputTag m_V0KsTag;
    edm::InputTag m_V0LTag;

    edm::EDGetTokenT<reco::BeamSpot> m_bsToken;
    edm::EDGetTokenT<vector<reco::Vertex>> m_offlinePVToken;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_GEN;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_SIM_GEANT;
    edm::EDGetTokenT<vector<reco::Track>> m_generalTracksToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_sCandsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0KsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0LToken;
    
    int verbose=1;
    
    TString a = ""; //"WjetsMC" "ZeroBias" "MET" "MinBiasMC" "SingleMuon"
    TString b = a + "";
    
    //--------- Histogram Declaration --------------------//
    std::map<TString, TH1F *> histos_th1f;
    std::map<TString, TH2F *> histos_th2f;
    std::map<TString, TEfficiency *> histos_teff;
    std::map<TString, TProfile *> histos_TProfile;

     };

#endif
