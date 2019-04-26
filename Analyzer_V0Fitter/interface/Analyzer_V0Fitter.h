#ifndef My_azimuthal_MC_h
#define My_azimuthal_MC_h
 
//#ifndef RECOVERTEX__V0_FITTER_H
//#define RECOVERTEX__V0_FITTER_H

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
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

//#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"
//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include <SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h>
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
 
class Analyzer_V0Fitter : public edm::EDAnalyzer
 {
  public:
    explicit Analyzer_V0Fitter(edm::ParameterSet const& theParams);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~Analyzer_V0Fitter();

    double openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2);
    double deltaR(double phi1, double eta1, double phi2, double eta2);
    double dz_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double lxy(TVector3 v1, TVector3 v2);
    TVector3 PCA_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double dxy_signed_line_point(TVector3 Point_line, TVector3 Vector_along_line, TVector3 Point);
    double phiInTrackRefPosition(int charge, double p1x, double p1y, double x1, double y1, double xREF, double yREF);
    bool specialGate(bool flag1, bool flag2);
  private:
    //---- configurable parameters --------
    //bool m_isData;

    //int m_nEvent, m_nRun, m_nLumi;
    
    bool vertexFitter_;
      bool useRefTracks_;
      bool doKShorts_;
      bool doLambdas_;

      // cuts on initial track selection
      double tkChi2Cut_;
      int tkNHitsCut_;
      double tkPtCut_;
      double tkIPSigXYCut_;
      double tkIPSigZCut_;
      // cuts on the vertex
      double vtxChi2Cut_;
      double vtxDecaySigXYCut_;
      double vtxDecaySigXYZCut_;
      // miscellaneous cuts
      double tkDCACut_;
      double mPiPiCut_;
      double innerHitPosCut_;
      double cosThetaXYCut_;
      double cosThetaXYZCut_;
      // cuts on the V0 candidate mass
      double kShortMassCut_;
      double lambdaMassCut_;
      bool useVertex_;

    edm::Service<TFileService> m_fs;
 
/*    edm::InputTag m_bsTag;
    edm::InputTag m_genParticlesTag_GEN;
    edm::InputTag m_genParticlesTag_SIM_GEANT;
    edm::InputTag m_generalTracksTag;
    edm::InputTag m_sCandsTag;
    edm::InputTag m_V0KsTag;
    edm::InputTag m_V0LTag;

    edm::EDGetTokenT<reco::BeamSpot> m_bsToken;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_GEN;
    edm::EDGetTokenT<vector<reco::GenParticle>> m_genParticlesToken_SIM_GEANT;
    edm::EDGetTokenT<vector<reco::Track>> m_generalTracksToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_sCandsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0KsToken;
    edm::EDGetTokenT<vector<reco::VertexCompositeCandidate> > m_V0LToken;
*/    
    int verbose=1;
    
    TString a = ""; //"WjetsMC" "ZeroBias" "MET" "MinBiasMC" "SingleMuon"
    TString b = a + "";
    
    //--------- Histogram Declaration --------------------//
    std::map<TString, TH1F *> histos_th1f;
    std::map<TString, TH2F *> histos_th2f;
    std::map<TString, TEfficiency *> histos_teff;
    //std::map<TString, TH3F *> histos_th3f;
    std::map<TString, TProfile *> histos_TProfile;

      edm::InputTag m_tag_tracks;
      edm::InputTag m_tag_beamspot;
      edm::InputTag m_tag_vertices;
      edm::InputTag m_tag_genParticlesTag_SIM_GEANT;
      edm::InputTag m_tag_TrackingParticlesTag;

      edm::EDGetTokenT<reco::TrackCollection> token_tracks;
      edm::EDGetTokenT<reco::BeamSpot> token_beamspot;
      edm::EDGetTokenT<std::vector<reco::Vertex>> token_vertices;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> token_genParticlesTag_SIM_GEANT;
      edm::EDGetTokenT<std::vector<TrackingParticle>> token_TrackingParticles;
      edm::EDGetTokenT<edm::View<reco::Track> >     token_recoTracks;
      edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> token_associatorByHits;
     };

#endif

