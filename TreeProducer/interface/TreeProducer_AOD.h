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
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// CMSSW specific lib
//#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

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
  edm::InputTag _lambdaKshortCollectionTag;

  edm::EDGetTokenT<vector<reco::Vertex> > _vertexCollectionToken;
  edm::EDGetTokenT<vector<reco::Track> > _trackCollectionToken;
  edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate> > _lambdaKshortCollectionToken;
  bool _isData;
  edm::EDGetToken  m_partons;

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
  std::vector<int> _track_Nhits, _track_NpixHits, _track_purity, _track_ndof;
  std::vector<double> _track_eta, _track_pt, _track_px, _track_py, _track_pz, _track_x, _track_y, _track_z, _track_phi, _track_ptError, _track_dxy, _track_d0, _track_dzError, _track_dz, _track_normalizedChi2;
  std::vector< std::vector<double> > _track_covariance;
  std::vector<int> _track_charge;

  // Lambdas
  int _lambda_N;
  std::vector<int> _lambda_ndof, _lambda_d1ch, _lambda_d2ch;
  std::vector<double> _lambda_normalizedChi2, _lambda_m;
  std::vector<double> _lambda_x, _lambda_y, _lambda_z;
  std::vector<double> _lambda_px, _lambda_py, _lambda_pz;
  std::vector<double> _lambda_d1px, _lambda_d1py, _lambda_d1pz;
  std::vector<double> _lambda_d2px, _lambda_d2py, _lambda_d2pz;

  // Kshorts
  int _kshort_N;
  std::vector<int> _kshort_ndof;
  std::vector<double> _kshort_normalizedChi2, _kshort_m;
  std::vector<double> _kshort_x, _kshort_y, _kshort_z;
  std::vector<double> _kshort_px, _kshort_py, _kshort_pz;
  std::vector<double> _kshort_d1px, _kshort_d1py, _kshort_d1pz;
  std::vector<double> _kshort_d2px, _kshort_d2py, _kshort_d2pz;

  //GenParticles
  std::vector<double> _genp_x;
  std::vector<double> _genp_y;
  std::vector<double> _genp_z;
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
  std::vector<int> _genp_mom;
  std::vector<int> _genp_m2;
  std::vector<int> _genp_d1;
  std::vector<int> _genp_d2;
  
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
