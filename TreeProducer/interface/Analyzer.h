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
#include "TFile.h"
#include "TProfile.h"
#include <TMath.h>

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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


  
class Analyzer : public edm::EDAnalyzer
 {
  public:
    explicit Analyzer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~Analyzer();

  private:
    //---- configurable parameters --------
    bool m_isData;

    int m_nEvent, m_nRun, m_nLumi;
    
    edm::Service<TFileService> m_fs ;
    
    edm::InputTag m_vertexCollectionTag;
    edm::InputTag m_trackCollectionTag;
    edm::InputTag m_partonsTag;
  
    edm::EDGetTokenT<vector<reco::Vertex> > m_vertexCollectionToken;
    edm::EDGetTokenT<vector<reco::Track> > m_trackCollectionToken;
    edm::EDGetTokenT<vector<reco::GenParticle> >  m_partons;
    int verbose=1;
    
    //--------- Histogram Declaration --------------------//
    // Vertices
    std::map<TString, TH1F*> histos_th1f;
    //TH1F *num_of_VtxGood;



 };


 namespace reco {
  template<typename T>
  class RecoPtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }
  };
}

#endif

