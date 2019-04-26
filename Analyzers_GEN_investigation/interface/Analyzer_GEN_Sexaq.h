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


  
class Analyzer_GEN_Sexaq : public edm::EDAnalyzer
 {
  public:
    explicit Analyzer_GEN_Sexaq(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~Analyzer_GEN_Sexaq();
    double openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2);
  private:
    //---- configurable parameters --------
    bool m_isData;

    int m_nEvent, m_nRun, m_nLumi;
    
    edm::Service<TFileService> m_fs;
 
    edm::InputTag m_genParticlesTag_GEN;
    edm::InputTag m_genParticlesTag_SIM;
    edm::InputTag m_genParticlesTag_PLUSGEANT;
    edm::InputTag m_genParticlesTag_HLT;
    edm::InputTag m_sCandsTag;

    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_GEN;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_SIM;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_PLUSGEANT;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_HLT;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_sCandsToken; 


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

