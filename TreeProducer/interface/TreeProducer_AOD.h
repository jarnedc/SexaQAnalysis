// -*- C++ -*-
// Package:    TreeProducer_AOD
// Class:      TreeProducer_AOD
///**\class TreeProducer_AOD TreeProducer_AOD.cc HexaAnalysis/TreeProducer/src/TreeProducer_AOD.cc
// Description: EDAnalyzer produce flat trees from AOD for HexaAnalysis

// C++ lib
#include <vector>

// ROOT
#include "TTree.h"
#include "TMatrixD.h"
//#include "TLorentzVector.h"
#include "TPRegexp.h"

// CMSSW standard lib
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// CMSSW specific lib
//#include "DataFormats/Common/interface/ValueMap.h"
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

// others
using namespace std;
int verbose=1;

//
// class declaration
//

class TreeProducer_AOD : public edm::EDAnalyzer {
 public:
  explicit TreeProducer_AOD(const edm::ParameterSet&);
  ~TreeProducer_AOD();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static bool ptSorter(const reco::Track & i, const reco::Track & j);

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void Init();

  // ----------member data ---------------------------
  edm::InputTag _trigResultsTag;
  edm::InputTag _vertexCollectionTag;
  edm::InputTag _trackCollectionTag;

  edm::EDGetTokenT<edm::TriggerResults> _trigResultsToken;
  edm::EDGetTokenT<vector<reco::Vertex> > _vertexCollectionToken;
  edm::EDGetTokenT<vector<reco::Track> > _trackCollectionToken;
  bool _isData;
  edm::EDGetToken  m_partons;

  HLTConfigProvider hltConfig_;
  HLTPrescaleProvider hltPrescale_;

  std::vector<std::string> triggerNames_;
  std::vector<unsigned int> triggerIndex_;

//   GlobalPoint vertexPosition;

  // Tree and its branches
  TTree* _tree;

  // Global quantities
  int _nEvent, _nRun, _nLumi, _nTrack, _nTrack_stored;

  // Vertices
  int _vtx_N, _vtx_N_stored;
  std::vector<int> _vtx_ndof, _vtx_nTracks;
  std::vector<int> _vtx_tracksSize;
  std::vector<bool> _vtx_isValid, _vtx_isFake;
  std::vector<double> _vtx_x, _vtx_y, _vtx_z;
  std::vector<double> _vtx_normalizedChi2, _vtx_d0;
  std::vector< std::vector<double> > _vtx_covariance;

	//Tracks
	std::vector<int> _track_fromPV, _track_Nhits, _track_NpixHits, _track_purity, _track_ndof;
	std::vector<double> _track_eta, _track_pt, _track_px, _track_py, _track_pz, _track_phi, _track_ptError, _track_dxy, _track_d0, _track_dzError, _track_dz, _track_normalizedChi2;
  std::vector< std::vector<double> > _track_covariance;
  std::vector<int> _track_charge;

  //GenParticles
  std::vector<double> _genp_px;
  std::vector<double> _genp_py;
  std::vector<double> _genp_pz;
  std::vector<double> _genp_pt;
  std::vector<double> _genp_p;
  std::vector<double> _genp_eta;
  std::vector<double> _genp_phi;
  std::vector<double> _genp_mass;
  std::vector<double> _genp_energy;
  std::vector<int> _genp_charge;
  std::vector<int> _genp_pdgid;
  std::vector<int> _genp_status;
  
  //Trigger
  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
  int _singlejet_450;
  //prescales
  double  _pswgt_singlejet_450;
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

//
// constants, enums and typedefs
//
namespace reco {
  typedef std::vector<Track> TrackCollection;
  typedef edm::Ref<TrackCollection> TrackRef;
}
//
// static data member definitions
//
