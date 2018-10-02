#include "SexaQAnalysis/TreeProducer/interface/TreeProducer_AOD.h"


//
// constructors and destructor
//

TreeProducer_AOD::TreeProducer_AOD(edm::ParameterSet const& pset):
trackCollectionTag_(pset.getParameter<edm::InputTag>("TrackCollection")),
lambdaCollectionTag_(pset.getParameter<edm::InputTag>("lambdaCollection")),
kshortCollectionTag_(pset.getParameter<edm::InputTag>("kshortCollection")),
offlinePrimaryVerticesCollectionTag_(pset.getParameter<edm::InputTag>("offlinePrimaryVerticesCollection")),
ak4PFJetsCollectionTag_(pset.getParameter<edm::InputTag>("ak4PFJetsCollection")),
muonsCollectionTag_(pset.getParameter<edm::InputTag>("muonsCollection")),
electronsCollectionTag_(pset.getParameter<edm::InputTag>("electronsCollection")),
METCollectionTag_(pset.getParameter<edm::InputTag>("METCollection")),
tracksCollectionToken_(consumes<std::vector<reco::Track> >(trackCollectionTag_)),
lambdaCollectionToken_(consumes<std::vector<reco::VertexCompositeCandidate> >(lambdaCollectionTag_)),
kshortCollectionToken_(consumes<std::vector<reco::VertexCompositeCandidate> >(kshortCollectionTag_)),
offlinePrimaryVerticesCollectionToken_(consumes<std::vector<reco::Vertex> >(offlinePrimaryVerticesCollectionTag_)),
ak4PFJetsCollectionToken_(consumes<std::vector<reco::PFJet> >(ak4PFJetsCollectionTag_)),
muonsCollectionToken_(consumes<std::vector<edm::FwdPtr<reco::PFCandidate> > >(muonsCollectionTag_)),
electronsCollectionToken_(consumes<std::vector<edm::FwdPtr<reco::PFCandidate> > >(electronsCollectionTag_)),
METCollectionToken_(consumes<std::vector<reco::PFMET> >(METCollectionTag_))
{
}

TreeProducer_AOD::~TreeProducer_AOD()
{
}

//
// member functions
//

// ------------ method called for each event  ------------
void
TreeProducer_AOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Initialize branches
  Init();

  // HANDLES //
  // Get collections
/*  edm::Handle<vector<reco::GenParticle> > H_partons;
  if(!_isData){
    iEvent.getByToken(m_partons, H_partons);
    if(!H_partons.isValid()) {
      if(verbose>0) cout << "Missing collection during TreeProducer_AOD: genParticles ... skip entry !" << endl;
      return;
    }
  }
*/
  //ntracks part
  edm::Handle<std::vector<reco::Track >> h_tracks;
  iEvent.getByToken(tracksCollectionToken_, h_tracks);
  if(!h_tracks.isValid()) {
      std::cout << "Missing collection during TreeProducer_AOD : " << trackCollectionTag_ << " ... skip entry !" << std::endl;
  }
  _nTrack = h_tracks->size();
/*  edm::Handle<vector<reco::VertexCompositePtrCandidate> > H_lk;
  if(_isData){
   iEvent.getByToken(_lambdaKshortCollectionToken , H_lk);
   if(!H_lk.isValid()) {
     if(verbose>0) cout << "Missing collection during TreeProducer_AOD: " << _lambdaKshortCollectionTag << " ... skip entry !" << endl;
     return;
   }
  }
*/

/*  edm::Handle<vector<reco::Vertex> > H_vert;
  iEvent.getByToken(_vertexCollectionToken, H_vert);
  if(!H_vert.isValid()) {
    if(verbose>0) cout << "Missing collection during TreeProducer_AOD: " << _vertexCollectionTag << " ... skip entry !" << endl;
    return;
  }

  edm::Handle<vector<reco::Track> > H_track;
  iEvent.getByToken(_trackCollectionToken , H_track);
  if(!H_track.isValid()) {
    if(verbose>0) cout << "Missing collection during TreeProducer_AOD: " << _trackCollectionTag << " ... skip entry !" << endl;
    return;
  }*/
/*
  edm::Handle<vector<reco::VertexCompositeCandidate> > H_S;
  if(_isData){
   iEvent.getByToken(_sCollectionToken , H_S);
   if(!H_S.isValid()) {
    if(verbose>0) cout << "Missing collection during TreeProducer_AOD : " << _sCollectionTag << " ... skip entry !" << endl;
    return;
   }
  }
  // throw away events on data withouts - for MC we check gen
  if (_isData && H_S->size() == 0) return; // only use events with at least one s
*/

/*  edm::Handle<vector<reco::Track> > H_S_tracks;
  if(_isData){
   iEvent.getByToken(_sTracksCollectionToken , H_S_tracks);
   if(!H_S_tracks.isValid()) {
    if(verbose>0) cout << "Missing collection during TreeProducer_AOD : " << _sTracksCollectionTag << " ... skip entry !" << endl;
    return;
   }
  }
  // throw away events on data withouts - for MC we check gen
  if (_isData && H_S_tracks->size() == 0) return; // only use events with at least one s
*/
  // GLOBAL EVENT INFORMATIONS //

  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();
  _nEvent = iEvent.id().event();




// counting of lambdas and kaons is now a little trickier, since we need to
// check overlaps to avoid doublecounting (otherwise #lambdas = #kaons by
// construction)
//  _lambda_N = H_lambda->size();
//  _kshort_N = H_kshort->size();

/*
  for(std::vector<reco::VertexCompositeCandidate>::const_iterator lambda = H_lambda->begin(); lambda != H_lambda->end(); ++lambda){
    if (lambda->numberOfDaughters() != 2) continue;
    for(std::vector<reco::VertexCompositeCandidate>::const_iterator kshort = H_kshort->begin(); kshort != H_kshort->end(); ++kshort){
      if (kshort->numberOfDaughters() != 2) continue;
      std::vector<reco::Track> v0Tracks;
std::cout << "Lambda: " << lambda->bestTrack()->charge() << " " << lambda->bestTrack()->pt()  << std::endl;
std::cout << "Kshort: " << kshort->bestTrack()->charge() << " " << kshort->bestTrack()->pt() << std::endl;
      v0Tracks.push_back(*lambda->bestTrack());
      v0Tracks.push_back(*kshort->bestTrack());
      KalmanVertexFitter fitter;
      TransientVertex vtx = fitter.vertex(v0Tracks);
std::cout << "Vertex: " << vtx.position().x() << " " << vtx.position().y() << " " << vtx.position().z() << std::endl;
    }
  }
*/

  _tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
TreeProducer_AOD::beginJob()
{
        // Initialize when class is created
        edm::Service<TFileService> fs ;
        _tree = fs->make <TTree>("SexaQAnalysis","tree");

        // Declare tree's branches
        // Event
        _tree->Branch("nEvent",&_nEvent,"nEvent/I");
        _tree->Branch("nRun",&_nRun,"nRun/I");
        _tree->Branch("nLumi",&_nLumi,"nLumi/I");
        //
        // Vertices
        _tree->Branch("nTrack",&_nTrack,"nTrack/I");
       
}

// ------------ method called once each job just after ending the event loop  ------------
void
TreeProducer_AOD::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
TreeProducer_AOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

// ------------ method called when ending the processing of a run  ------------
void
TreeProducer_AOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
TreeProducer_AOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
TreeProducer_AOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeProducer_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
TreeProducer_AOD::Init()
{
  _nEvent = _nRun = _nLumi = 0;

  //Tracks
  _nTrack = 0;

}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer_AOD);
