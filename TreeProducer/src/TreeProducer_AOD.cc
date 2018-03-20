#include "SexaQAnalysis/TreeProducer/interface/TreeProducer_AOD.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"


//
// constructors and destructor
//
TreeProducer_AOD::TreeProducer_AOD(const edm::ParameterSet& pset):
  _vertexCollectionTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  _trackCollectionTag(pset.getParameter<edm::InputTag>("trackCollection")),
  _lambdaKshortCollectionTag(pset.getParameter<edm::InputTag>("lambdaKshortCollection")),
  _sCollectionTag(pset.getParameter<edm::InputTag>("sCollection")),
  _sTracksCollectionTag(pset.getParameter<edm::InputTag>("sTrackCollection")),
  _vertexCollectionToken(consumes<vector<reco::Vertex> >(_vertexCollectionTag)),
  _trackCollectionToken(consumes<vector<reco::Track> >(_trackCollectionTag)),
  _lambdaKshortCollectionToken(consumes<vector<reco::VertexCompositePtrCandidate> >(_lambdaKshortCollectionTag)),
  _sCollectionToken(consumes<vector<reco::VertexCompositeCandidate> >(_sCollectionTag)),
  _sTracksCollectionToken(consumes<vector<reco::Track> >(_sTracksCollectionTag)),
  _isData(pset.getUntrackedParameter<bool>("isData")),
  m_partons(consumes<vector<reco::GenParticle> >(pset.getParameter<edm::InputTag>("genCollection")))
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
  edm::Handle<vector<reco::GenParticle> > H_partons;
  if(!_isData){
    iEvent.getByToken(m_partons, H_partons);
    if(!H_partons.isValid()) {
      if(verbose>0) cout << "Missing collection during TreeProducer_AOD: genParticles ... skip entry !" << endl;
      return;
    }
  }


  edm::Handle<vector<reco::VertexCompositePtrCandidate> > H_lk;
  if(_isData){
   iEvent.getByToken(_lambdaKshortCollectionToken , H_lk);
   if(!H_lk.isValid()) {
     if(verbose>0) cout << "Missing collection during TreeProducer_AOD: " << _lambdaKshortCollectionTag << " ... skip entry !" << endl;
     return;
   }
  }


  edm::Handle<vector<reco::Vertex> > H_vert;
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
  }

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


  edm::Handle<vector<reco::Track> > H_S_tracks;
  if(_isData){
   iEvent.getByToken(_sTracksCollectionToken , H_S_tracks);
   if(!H_S_tracks.isValid()) {
    if(verbose>0) cout << "Missing collection during TreeProducer_AOD : " << _sTracksCollectionTag << " ... skip entry !" << endl;
    return;
   }
  }
  // throw away events on data withouts - for MC we check gen
  if (_isData && H_S_tracks->size() == 0) return; // only use events with at least one s

  // GLOBAL EVENT INFORMATIONS //

  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();
  _nEvent = iEvent.id().event();


  // GENPARTICLES //

  if(!_isData){
    bool hasLambda0 = false;
    for (vector<reco::GenParticle>::const_iterator theparton = H_partons->begin(); theparton != H_partons->end(); ++theparton){
      // ask for a lambda
      //  -> within tracker acceptance (otherwise at least one decay product
      //     not in tracker) ; we could tighten this, asking p an pi in tracker
      //  -> and with some minimum pt (m(Lambda)-m(p)-m(pi) is tiny.
      if (abs(theparton->pdgId())==3122 &&
          abs(theparton->eta())<2.5 && theparton->pt()>0.5) {
        hasLambda0=true;
        break;
      }
    }
    if (!hasLambda0) return;
    for (vector<reco::GenParticle>::const_iterator theparton = H_partons->begin(); theparton != H_partons->end(); ++theparton){
      _genp_x.push_back(theparton->vx());
      _genp_y.push_back(theparton->vy());
      _genp_z.push_back(theparton->vz());
      _genp_px.push_back(theparton->px());
      _genp_py.push_back(theparton->py());
      _genp_pz.push_back(theparton->pz());
      _genp_p.push_back(theparton->p());
      _genp_pt.push_back(theparton->pt());
      _genp_eta.push_back(theparton->eta());
      _genp_phi.push_back(theparton->phi());
      _genp_mass.push_back(theparton->mass());
      _genp_energy.push_back(theparton->energy());
      _genp_charge.push_back(theparton->charge());
      _genp_pdgid.push_back(theparton->pdgId());
      _genp_status.push_back(theparton->status());
      _genp_mom.push_back(theparton->numberOfMothers()>0 ? theparton->motherRef(0).key() : 0);
      _genp_m2.push_back(theparton->numberOfMothers()>1 ? theparton->motherRef(1).key() : 0);
      _genp_d1.push_back(theparton->numberOfDaughters()>0 ? theparton->daughterRef(0).key() : 0);
      _genp_d2.push_back(theparton->numberOfDaughters()>1 ? theparton->daughterRef(1).key() : 0);
    }///for partons
  }



  // VERTICES //

  for( std::vector<reco::Vertex>::const_iterator PV = H_vert->begin(); PV != H_vert->end(); ++PV){
    _vtx_normalizedChi2.push_back(PV->normalizedChi2());
    _vtx_ndof.push_back(PV->ndof());
    _vtx_nTracks.push_back(PV->nTracks());
    _vtx_tracksSize.push_back(PV->tracksSize());
    _vtx_d0.push_back(PV->position().Rho());
    _vtx_x.push_back(PV->x());
    _vtx_y.push_back(PV->y());
    _vtx_z.push_back(PV->z());
    _vtx_isValid.push_back(PV->isValid());
    _vtx_isFake.push_back(PV->isFake());
    std::vector<double> cov;
    for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
        cov.push_back(PV->covariance(i,j));
      }
    }
    _vtx_covariance.push_back(cov);
  } // for loop on primary vertices
  _vtx_N = H_vert->size();


  // TRACKS //
  for (vector<reco::Track>::const_iterator theTrack = H_track->begin(); theTrack != H_track->end(); ++theTrack){
    _track_purity.push_back(theTrack->highPurity);
    _track_Nhits.push_back(theTrack->numberOfValidHits());
    _track_NpixHits.push_back(theTrack->hitPattern().numberOfValidPixelHits());
    _track_pt.push_back(theTrack->pt());
    _track_px.push_back(theTrack->px());
    _track_py.push_back(theTrack->py());
    _track_pz.push_back(theTrack->pz());
    _track_x.push_back(theTrack->vx());
    _track_y.push_back(theTrack->vy());
    _track_z.push_back(theTrack->vz());
    _track_eta.push_back(theTrack->eta());
    _track_phi.push_back(theTrack->phi());
    _track_normalizedChi2.push_back(theTrack->normalizedChi2());
    _track_ndof.push_back(theTrack->ndof());
    _track_ptError.push_back(theTrack->ptError());
    _track_dzError.push_back(theTrack->dzError());
    _track_dz.push_back(theTrack->dz());
    _track_dxy.push_back(theTrack->dxy());
    _track_d0.push_back(theTrack->d0());
    _track_charge.push_back(theTrack->charge());
    std::vector<double> cov;
    for(int k=0; k<4; k++){
      for(int j=0; j<4; j++){
        cov.push_back(theTrack->covariance(k,j));
      }
    }
    _track_covariance.push_back(cov);
  }
  _nTrack = H_track->size();

// S //
  if(_isData){
  for(std::vector<reco::VertexCompositeCandidate>::const_iterator s = H_S->begin(); s != H_S->end(); ++s){
    if (s->numberOfDaughters() == 2) {
      _S_normalizedChi2.push_back(s->vertexNormalizedChi2());
      _S_ndof.push_back(s->vertexNdof());
      _S_x.push_back(s->vx());
      _S_y.push_back(s->vy());
      _S_z.push_back(s->vz());
      _S_px.push_back(s->px());
      _S_py.push_back(s->py());
      _S_pz.push_back(s->pz());
      _S_m.push_back(s->mass());
      _S_d1px.push_back(s->daughter(0)->px());
      _S_d1py.push_back(s->daughter(0)->py());
      _S_d1pz.push_back(s->daughter(0)->pz());
      _S_d2px.push_back(s->daughter(1)->px());
      _S_d2py.push_back(s->daughter(1)->py());
      _S_d2pz.push_back(s->daughter(1)->pz());
    } else {
      cout << "%%%s vertex with " << s->numberOfDaughters() << " daughters" << std::endl;
      continue;
    }
  } // for loop on s
  _S_N = H_S->size();
  }


  // S TRACKS //
  if(_isData){
  for (vector<reco::Track>::const_iterator theTrack = H_S_tracks->begin(); theTrack != H_S_tracks->end(); ++theTrack){
    _Strack_purity.push_back(theTrack->highPurity);
    _Strack_Nhits.push_back(theTrack->numberOfValidHits());
    _Strack_NpixHits.push_back(theTrack->hitPattern().numberOfValidPixelHits());
    _Strack_pt.push_back(theTrack->pt());
    _Strack_px.push_back(theTrack->px());
    _Strack_py.push_back(theTrack->py());
    _Strack_pz.push_back(theTrack->pz());
    _Strack_x.push_back(theTrack->vx());
    _Strack_y.push_back(theTrack->vy());
    _Strack_z.push_back(theTrack->vz());
    _Strack_eta.push_back(theTrack->eta());
    _Strack_phi.push_back(theTrack->phi());
    _Strack_normalizedChi2.push_back(theTrack->normalizedChi2());
    _Strack_ndof.push_back(theTrack->ndof());
    _Strack_ptError.push_back(theTrack->ptError());
    _Strack_dzError.push_back(theTrack->dzError());
    _Strack_dz.push_back(theTrack->dz());
    _Strack_dxy.push_back(theTrack->dxy());
    _Strack_d0.push_back(theTrack->d0());
    _Strack_charge.push_back(theTrack->charge());
    std::vector<double> cov;
    for(int k=0; k<4; k++){
      for(int j=0; j<4; j++){
        cov.push_back(theTrack->covariance(k,j));
      }
    }
    _Strack_covariance.push_back(cov);
  }
  _SnTrack = H_S_tracks->size();
  }
/*
  // LAMBDAS //

  for(vector<reco::VertexCompositePtrCandidate>::const_iterator lk = H_lk->begin(); lk != H_lk->end(); ++lk){
    const reco::VertexCompositeCandidate * lambda = dynamic_cast<const reco::VertexCompositeCandidate *>(lk->daughter(0));
    if (lambda->numberOfDaughters() == 2) {
      _lambda_normalizedChi2.push_back(lambda->vertexNormalizedChi2());
      _lambda_ndof.push_back(lambda->vertexNdof());
      _lambda_x.push_back(lambda->vx());
      _lambda_y.push_back(lambda->vy());
      _lambda_z.push_back(lambda->vz());
      _lambda_px.push_back(lambda->px());
      _lambda_py.push_back(lambda->py());
      _lambda_pz.push_back(lambda->pz());
      _lambda_m.push_back(lambda->mass());
      _lambda_d1px.push_back(lambda->daughter(0)->px());
      _lambda_d1py.push_back(lambda->daughter(0)->py());
      _lambda_d1pz.push_back(lambda->daughter(0)->pz());
      _lambda_d1ch.push_back(lambda->daughter(0)->charge());
      _lambda_d2px.push_back(lambda->daughter(1)->px());
      _lambda_d2py.push_back(lambda->daughter(1)->py());
      _lambda_d2pz.push_back(lambda->daughter(1)->pz());
      _lambda_d2ch.push_back(lambda->daughter(1)->charge());
    } else {
      cout << "%%% Lambda vertex with " << lambda->numberOfDaughters() << " daughters" << std::endl;
      continue;
    }

    const reco::VertexCompositeCandidate * kshort = dynamic_cast<const reco::VertexCompositeCandidate *>(lk->daughter(1));
    if (kshort->numberOfDaughters() == 2) {
      _kshort_normalizedChi2.push_back(kshort->vertexNormalizedChi2());
      _kshort_ndof.push_back(kshort->vertexNdof());
      _kshort_x.push_back(kshort->vx());
      _kshort_y.push_back(kshort->vy());
      _kshort_z.push_back(kshort->vz());
      _kshort_px.push_back(kshort->px());
      _kshort_py.push_back(kshort->py());
      _kshort_pz.push_back(kshort->pz());
      _kshort_m.push_back(kshort->mass());
      _kshort_d1px.push_back(kshort->daughter(0)->px());
      _kshort_d1py.push_back(kshort->daughter(0)->py());
      _kshort_d1pz.push_back(kshort->daughter(0)->pz());
      _kshort_d2px.push_back(kshort->daughter(1)->px());
      _kshort_d2py.push_back(kshort->daughter(1)->py());
      _kshort_d2pz.push_back(kshort->daughter(1)->pz());
    } else {
      cout << "%%% kshort vertex with " << kshort->numberOfDaughters() << " daughters" << std::endl;
      continue;
    }
  } // for loop on lambda-kshort pairs
*/

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
        _tree = fs->make <TTree>("HexaQAnalysis","tree");

        // Declare tree's branches
        // Event
        _tree->Branch("nEvent",&_nEvent,"nEvent/I");
        _tree->Branch("nRun",&_nRun,"nRun/I");
        _tree->Branch("nLumi",&_nLumi,"nLumi/I");
        //
        // Vertices
        _tree->Branch("vtx_N",&_vtx_N,"vtx_N/I");
        _tree->Branch("vtx_N_stored",&_vtx_N_stored,"vtx_N_stored/I");
        _tree->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2);
        _tree->Branch("vtx_ndof",&_vtx_ndof);
        _tree->Branch("vtx_nTracks",&_vtx_nTracks);
        _tree->Branch("vtx_d0",&_vtx_d0);
        _tree->Branch("vtx_x",&_vtx_x);
        _tree->Branch("vtx_y",&_vtx_y);
        _tree->Branch("vtx_z",&_vtx_z);
        _tree->Branch("vtx_covariance",&_vtx_covariance);
        _tree->Branch("vtx_isFake",&_vtx_isFake);
        _tree->Branch("vtx_isValid",&_vtx_isValid);
        //
	// Tracks
        _tree->Branch("nTrack_stored",&_nTrack_stored,"nTrack_stored/I");
        _tree->Branch("nTrack",&_nTrack,"nTrack/I");
        _tree->Branch("track_pt",&_track_pt);
        _tree->Branch("track_px",&_track_px);
        _tree->Branch("track_py",&_track_py);
        _tree->Branch("track_pz",&_track_pz);
        _tree->Branch("track_x",&_track_x);
        _tree->Branch("track_y",&_track_y);
        _tree->Branch("track_z",&_track_z);
        _tree->Branch("track_eta",&_track_eta);
	_tree->Branch("track_phi",&_track_phi);
	_tree->Branch("track_normalizedChi2",&_track_normalizedChi2);
	_tree->Branch("track_ndof",&_track_ndof);
	_tree->Branch("track_ptError",&_track_ptError);
	_tree->Branch("track_dzError",&_track_dzError);
	_tree->Branch("track_dz",&_track_dz);
	_tree->Branch("track_purity",&_track_purity);
	_tree->Branch("track_nhits",&_track_Nhits);
	_tree->Branch("track_nPixHits",&_track_NpixHits);
        _tree->Branch("track_d0",&_track_d0);
        _tree->Branch("track_charge",&_track_charge);
        _tree->Branch("track_dxy",&_track_dxy);
	_tree->Branch("track_covariance",&_track_covariance);
	//
	// Lambdas
	_tree->Branch("nLambda",&_lambda_N,"nLambda/I");
	_tree->Branch("lambda_normalizedChi2",&_lambda_normalizedChi2);
	_tree->Branch("lambda_ndof",&_lambda_ndof);
	_tree->Branch("lambda_x",&_lambda_x);
	_tree->Branch("lambda_y",&_lambda_y);
	_tree->Branch("lambda_z",&_lambda_z);
	_tree->Branch("lambda_px",&_lambda_px);
	_tree->Branch("lambda_py",&_lambda_py);
	_tree->Branch("lambda_pz",&_lambda_pz);
	_tree->Branch("lambda_m",&_lambda_m);
	_tree->Branch("lambda_d1px",&_lambda_d1px);
	_tree->Branch("lambda_d1py",&_lambda_d1py);
	_tree->Branch("lambda_d1pz",&_lambda_d1pz);
	_tree->Branch("lambda_d1ch",&_lambda_d1ch);
	_tree->Branch("lambda_d2px",&_lambda_d2px);
	_tree->Branch("lambda_d2py",&_lambda_d2py);
	_tree->Branch("lambda_d2pz",&_lambda_d2pz);
	_tree->Branch("lambda_d2ch",&_lambda_d2ch);
	//
	// Kshorts
	_tree->Branch("nKshort",&_kshort_N,"nKshort/I");
	_tree->Branch("kshort_normalizedChi2",&_kshort_normalizedChi2);
	_tree->Branch("kshort_ndof",&_kshort_ndof);
	_tree->Branch("kshort_x",&_kshort_x);
	_tree->Branch("kshort_y",&_kshort_y);
	_tree->Branch("kshort_z",&_kshort_z);
	_tree->Branch("kshort_px",&_kshort_px);
	_tree->Branch("kshort_py",&_kshort_py);
	_tree->Branch("kshort_pz",&_kshort_pz);
	_tree->Branch("kshort_m",&_kshort_m);
	_tree->Branch("kshort_d1px",&_kshort_d1px);
	_tree->Branch("kshort_d1py",&_kshort_d1py);
	_tree->Branch("kshort_d1pz",&_kshort_d1pz);
	_tree->Branch("kshort_d2px",&_kshort_d2px);
	_tree->Branch("kshort_d2py",&_kshort_d2py);
	_tree->Branch("kshort_d2pz",&_kshort_d2pz);

	// s
	_tree->Branch("nS",&_S_N,"nKshort/I");
	_tree->Branch("S_normalizedChi2",&_S_normalizedChi2);
	_tree->Branch("S_ndof",&_S_ndof);
	_tree->Branch("S_x",&_S_x);
	_tree->Branch("S_y",&_S_y);
	_tree->Branch("S_z",&_S_z);
	_tree->Branch("S_px",&_S_px);
	_tree->Branch("S_py",&_S_py);
	_tree->Branch("S_pz",&_S_pz);
	_tree->Branch("S_m",&_S_m);
	_tree->Branch("S_d1px",&_S_d1px);
	_tree->Branch("S_d1py",&_S_d1py);
	_tree->Branch("S_d1pz",&_S_d1pz);
	_tree->Branch("S_d2px",&_S_d2px);
	_tree->Branch("S_d2py",&_S_d2py);
	_tree->Branch("S_d2pz",&_S_d2pz);

    ///GenParticles
	_tree->Branch("gen_x",&_genp_x);
	_tree->Branch("gen_y",&_genp_y);
	_tree->Branch("gen_z",&_genp_z);
	_tree->Branch("gen_px",&_genp_px);
	_tree->Branch("gen_py",&_genp_py);
	_tree->Branch("gen_pz",&_genp_pz);
        _tree->Branch("gen_p",&_genp_p);
        _tree->Branch("gen_pt",&_genp_pt);
	_tree->Branch("gen_eta",&_genp_eta);
	_tree->Branch("gen_phi",&_genp_phi);
	_tree->Branch("gen_mass",&_genp_mass);
        _tree->Branch("gen_energy",&_genp_energy);
        _tree->Branch("gen_charge",&_genp_charge);
        _tree->Branch("gen_pdgid",&_genp_pdgid);
        _tree->Branch("gen_status",&_genp_status);
        _tree->Branch("gen_mom",&_genp_mom);
        _tree->Branch("gen_m2",&_genp_m2);
        _tree->Branch("gen_d1",&_genp_d1);
        _tree->Branch("gen_d2",&_genp_d2);
 
        ///S tracks
       	// Tracks
        _tree->Branch("SnTrack_stored",&_SnTrack_stored,"SnTrack_stored/I");
        _tree->Branch("SnTrack",&_SnTrack,"SnTrack/I");
        _tree->Branch("Strack_pt",&_Strack_pt);
        _tree->Branch("Strack_px",&_Strack_px);
        _tree->Branch("Strack_py",&_Strack_py);
        _tree->Branch("Strack_pz",&_Strack_pz);
        _tree->Branch("Strack_x",&_Strack_x);
        _tree->Branch("Strack_y",&_Strack_y);
        _tree->Branch("Strack_z",&_Strack_z);
        _tree->Branch("Strack_eta",&_Strack_eta);
	_tree->Branch("Strack_phi",&_Strack_phi);
	_tree->Branch("Strack_normalizedChi2",&_Strack_normalizedChi2);
	_tree->Branch("Strack_ndof",&_Strack_ndof);
	_tree->Branch("Strack_ptError",&_Strack_ptError);
	_tree->Branch("Strack_dzError",&_Strack_dzError);
	_tree->Branch("Strack_dz",&_Strack_dz);
	_tree->Branch("Strack_purity",&_Strack_purity);
	_tree->Branch("Strack_nhits",&_Strack_Nhits);
	_tree->Branch("Strack_nPixHits",&_Strack_NpixHits);
        _tree->Branch("Strack_d0",&_Strack_d0);
        _tree->Branch("Strack_charge",&_Strack_charge);
        _tree->Branch("Strack_dxy",&_Strack_dxy);
	_tree->Branch("Strack_covariance",&_Strack_covariance);
       
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

  // Vertices
  _vtx_N = 0;
  _vtx_N_stored = 0;
  _vtx_normalizedChi2.clear();
  _vtx_ndof.clear();
  _vtx_nTracks.clear();
  _vtx_d0.clear();
  _vtx_x.clear();
  _vtx_y.clear();
  _vtx_z.clear();
  _vtx_covariance.clear();
  _vtx_isFake.clear();
  _vtx_isValid.clear();


  //Tracks
  _nTrack = 0;
  _nTrack_stored = 0;
  _track_eta.clear();
  _track_ndof.clear();
  _track_Nhits.clear();
  _track_normalizedChi2.clear();
  _track_NpixHits.clear();
  _track_phi.clear();
  _track_pt.clear();
  _track_px.clear();
  _track_py.clear();
  _track_pz.clear();  
  _track_x.clear();
  _track_y.clear();
  _track_z.clear();
  _track_ptError.clear();
  _track_dzError.clear();
  _track_dz.clear();
  _track_d0.clear();
  _track_charge.clear();
  _track_covariance.clear();
  _track_dxy.clear();
  _track_purity.clear();
  
  //Lambdas
  _lambda_N = 0;
  _lambda_normalizedChi2.clear();
  _lambda_ndof.clear();
  _lambda_x.clear();
  _lambda_y.clear();
  _lambda_z.clear();
  _lambda_px.clear();
  _lambda_py.clear();
  _lambda_pz.clear();
  _lambda_m.clear();
  _lambda_d1px.clear();
  _lambda_d1py.clear();
  _lambda_d1pz.clear();
  _lambda_d1ch.clear();
  _lambda_d2px.clear();
  _lambda_d2py.clear();
  _lambda_d2pz.clear();
  _lambda_d2ch.clear();
  
  //Kshorts
  _kshort_N = 0;
  _kshort_normalizedChi2.clear();
  _kshort_ndof.clear();
  _kshort_x.clear();
  _kshort_y.clear();
  _kshort_z.clear();
  _kshort_px.clear();
  _kshort_py.clear();
  _kshort_pz.clear();
  _kshort_m.clear();
  _kshort_d1px.clear();
  _kshort_d1py.clear();
  _kshort_d1pz.clear();
  _kshort_d2px.clear();
  _kshort_d2py.clear();
  _kshort_d2pz.clear();
  
  //S
  _S_N = 0;
  _S_normalizedChi2.clear();
  _S_ndof.clear();
  _S_x.clear();
  _S_y.clear();
  _S_z.clear();
  _S_px.clear();
  _S_py.clear();
  _S_pz.clear();
  _S_m.clear();
  _S_d1px.clear();
  _S_d1py.clear();
  _S_d1pz.clear();
  _S_d2px.clear();
  _S_d2py.clear();
  _S_d2pz.clear();

  //GenParticles
  _genp_x.clear();
  _genp_y.clear();
  _genp_z.clear();
  _genp_px.clear();
  _genp_py.clear();
  _genp_pz.clear();
  _genp_p.clear();
  _genp_pt.clear();  
  _genp_eta.clear();
  _genp_phi.clear();
  _genp_mass.clear();
  _genp_energy.clear();
  _genp_charge.clear();
  _genp_pdgid.clear();
  _genp_status.clear();
  _genp_mom.clear();
  _genp_m2.clear();
  _genp_d1.clear();
  _genp_d2.clear();

    //Tracks
  _SnTrack = 0;
  _SnTrack_stored = 0;
  _Strack_eta.clear();
  _Strack_ndof.clear();
  _Strack_Nhits.clear();
  _Strack_normalizedChi2.clear();
  _Strack_NpixHits.clear();
  _Strack_phi.clear();
  _Strack_pt.clear();
  _Strack_px.clear();
  _Strack_py.clear();
  _Strack_pz.clear();  
  _Strack_x.clear();
  _Strack_y.clear();
  _Strack_z.clear();
  _Strack_ptError.clear();
  _Strack_dzError.clear();
  _Strack_dz.clear();
  _Strack_d0.clear();
  _Strack_charge.clear();
  _Strack_covariance.clear();
  _Strack_dxy.clear();
  _Strack_purity.clear();

}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer_AOD);
