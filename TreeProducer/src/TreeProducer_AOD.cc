#include "HexaAnalysis/TreeProducer/interface/TreeProducer_AOD.h"

//
// constructors and destructor
//
TreeProducer_AOD::TreeProducer_AOD(const edm::ParameterSet& pset):
  _trigResultsTag(pset.getParameter<edm::InputTag>("triggerResults")),
  _genjetCollectionTag(pset.getParameter<edm::InputTag>("genjetCollection")),
  _vertexCollectionTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  _trackCollectionTag(pset.getParameter<edm::InputTag>("trackCollection")),
  _trigResultsToken(consumes<edm::TriggerResults>(_trigResultsTag)),
  _genjetCollectionToken(consumes<vector<reco::GenJet> >(_genjetCollectionTag)),
  _vertexCollectionToken(consumes<vector<reco::Vertex> >(_vertexCollectionTag)),
  _trackCollectionToken(consumes<vector<reco::Track> >(_trackCollectionTag)),
  _isData(pset.getUntrackedParameter<bool>("isData")),
  hltPrescale_(pset, consumesCollector(), *this)
{
 triggerNames_ = pset.getParameter<std::vector<std::string> > ("triggerName");


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
  edm::Handle<edm::TriggerResults> H_trig;//, H_trig1, H_trig2;
  iEvent.getByToken(_trigResultsToken, H_trig);

  edm::Handle<vector<reco::Vertex> > H_vert;
  iEvent.getByToken(_vertexCollectionToken, H_vert);

  edm::Handle<vector<reco::GenJet> > H_genjets;
  iEvent.getByToken(_genjetCollectionToken , H_genjets);

  edm::Handle<vector<reco::Track> > H_track;
  iEvent.getByToken(_trackCollectionToken , H_track);

  // Check validity
  if(!H_trig.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _trigResultsTag << " ... skip entry !" << endl;
    return;
  }

  if(!H_vert.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _vertexCollectionTag << " ... skip entry !" << endl;
    return;
  }

  if(!H_track.isValid()) {
    if(verbose>0) cout << "Missing collection : " << _trackCollectionTag << " ... skip entry !" << endl;
    return;
  }

  // GLOBAL EVENT INFORMATIONS //
  _nRun   = iEvent.id().run();
  _nLumi  = iEvent.luminosityBlock();
  _nEvent = iEvent.id().event();

  // VERTICES //
  UInt_t vtx_counter=0;
  _vtx_N = H_vert->size();
  _vtx_N_stored = nV;

  // select the primary vertex as the one with higest sum of (pt)^2 of tracks
  PrimaryVertexSorter PVSorter;
  std::vector<reco::Vertex> sortedVertices = PVSorter.sortedList( *(H_vert.product()) );

  for( std::vector<reco::Vertex>::const_iterator PV = sortedVertices.begin(); PV != sortedVertices.end(); ++PV){
    _vtx_normalizedChi2.push_back(PV->normalizedChi2());
    _vtx_ndof.push_back(PV->ndof());
    _vtx_nTracks.push_back(PV->tracksSize());
//     _vtx_nTracks.push_back(PV->nTracks());
    _vtx_d0.push_back(PV->position().Rho());
    _vtx_x.push_back(PV->x());
    _vtx_y.push_back(PV->y());
    _vtx_z.push_back(PV->z());

    vtx_counter++;
    if(vtx_counter >= nV) break;
  } // for loop on primary vertices

  // TRACKS //
	vector<reco::TrackRef> trackRef;
	for (vector<reco::Track>::const_iterator theTrack = H_track->begin(); theTrack != H_track->end(); ++theTrack){
		reco::TrackRef ref(H_track, theTrack - H_track->begin());
			trackRef.push_back(ref);
	}

	reco::RecoPtSorter<reco::TrackRef> trackSorter;
	std::sort( trackRef.begin(), trackRef.end(), trackSorter);

  UInt_t iT=0;
	for (size_t i = 0; i < trackRef.size(); i++) {
		_track_purity.push_back(trackRef[i]->highPurity);
		_track_Nhits.push_back(trackRef[i]->numberOfValidHits());
		_track_NpixHits.push_back(trackRef[i]->hitPattern().numberOfValidPixelHits());
		if (trackRef[i]->vz() == 0) _track_fromPV.push_back(1);
		else _track_fromPV.push_back(0);
		_track_pt.push_back(trackRef[i]->pt());
		_track_eta.push_back(trackRef[i]->eta());
		_track_phi.push_back(trackRef[i]->phi());
		_track_normalizedChi2.push_back(trackRef[i]->normalizedChi2());
		_track_ndof.push_back(trackRef[i]->ndof());
		_track_ptError.push_back(trackRef[i]->ptError());
		_track_dzError.push_back(trackRef[i]->dzError());
		_track_dz.push_back(trackRef[i]->dz());
		_track_dxy.push_back(trackRef[i]->dxy());
		_track_d0.push_back(trackRef[i]->d0());
		iT++ ;
    if(iT>=nT) break;
	}

  _nTrack = iT;
  _nTrack_stored = nT;

  UInt_t iGJ=0;
  //
	if(!_isData){
    for (vector<reco::GenJet>::const_iterator theGenJet = H_genjets->begin(); theGenJet != H_genjets->end(); ++theGenJet){
      if (theGenJet->pt() > 20){
        //Area
        _genjet_area.push_back(theGenJet->jetArea());

        // Vertex
        _genjet_vx.push_back(theGenJet->vx());
        _genjet_vy.push_back(theGenJet->vy());
        _genjet_vz.push_back(theGenJet->vz());

        // Kinematics
        _genjet_pt.push_back(theGenJet->pt());
        _genjet_eta.push_back(theGenJet->eta());
        _genjet_phi.push_back(theGenJet->phi());
        _genjet_e.push_back(theGenJet->energy());
        _genjet_m.push_back(theGenJet->mass());

        // Energy fractions
        double efrac = 0;
        for (size_t i = 0; i < theGenJet->numberOfDaughters(); i++){
          const reco::Candidate * constituent = theGenJet->daughter(i);
          if (constituent->charge() != 0) efrac += constituent->energy()/theGenJet->energy();
        }
        _genjet_efrac_ch.push_back(efrac);
        iGJ++ ;
      }
      if(iGJ>=nGJ) break;
    }
    _nGenJet = iGJ;
    _nGenJet_stored = nGJ;
  }

  //TRIGGER//
  if (triggerPathsMap[triggerPathsVector[0]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[0]])) _dijet_170_0p1 = 1;
  if (triggerPathsMap[triggerPathsVector[1]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[1]])) _dijet_220_0p3 = 1;
  if (triggerPathsMap[triggerPathsVector[2]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[2]])) _dijet_330_0p5 = 1;
  if (triggerPathsMap[triggerPathsVector[3]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[3]])) _dijet_430 = 1;
  if (triggerPathsMap[triggerPathsVector[4]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[4]])) _dijet_170 = 1;
  if (triggerPathsMap[triggerPathsVector[5]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[5]])) _singlejet_170_0p1 = 1;
  if (triggerPathsMap[triggerPathsVector[6]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[6]])) _singlejet_450 = 1;
  if (triggerPathsMap[triggerPathsVector[7]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[7]])) _singlejet_500 = 1;
  if (triggerPathsMap[triggerPathsVector[8]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[8]])) _isomu24 = 1;
  if (triggerPathsMap[triggerPathsVector[9]] != -1 && H_trig->accept(triggerPathsMap[triggerPathsVector[9]])) _isomu27 = 1;

  //prescales
//   const edm::TriggerNames &trignames = iEvent.triggerNames(*H_trig);
// 	for (size_t i = 0; i < H_trig->size(); i++) {
// 			if (trignames.triggerName(i).find("HLT_DiCentralPFJet170_v") != string::npos) _pswgt_dijet_170 = H_prescale->getPrescaleForIndex(i);
// 			if (trignames.triggerName(i).find("HLT_SingleCentralPFJet170_CFMax0p1_v") != string::npos) _pswgt_singlejet_170_0p1 = H_prescale->getPrescaleForIndex(i);
// 	}


// sanity check
assert(H_trig->size() == hltConfig_.size());
//------ loop over all trigger names ---------
for(unsigned itrig=0;itrig<triggerNames_.size() ;itrig++) {
  int preL1(-1);
  int preHLT(-1);

  if (triggerIndex_[itrig] < hltConfig_.size()) {
    ///In detail get prescale info from hltConfig_
    std::pair<std::vector<std::pair<std::string,int> >,int> detailedPrescaleInfo = hltPrescale_.prescaleValuesInDetail(iEvent, iSetup, triggerNames_[itrig]);
    preHLT = detailedPrescaleInfo.second ;

    // save l1 prescale values in standalone vector
    std::vector <int> l1prescalevals;
    for( size_t varind = 0; varind < detailedPrescaleInfo.first.size(); varind++ ){
      l1prescalevals.push_back(detailedPrescaleInfo.first.at(varind).second);
    }

    //find and save minimum l1 prescale of any ORed L1 that seeds the HLT
    std::vector<int>::iterator result = std::min_element(std::begin(l1prescalevals), std::end(l1prescalevals));
    size_t minind = std::distance(std::begin(l1prescalevals), result);
    // sometimes there are no L1s associated with a HLT. In that case, this branch stores -1 for the l1prescale
    preL1 = minind < l1prescalevals.size() ? l1prescalevals.at(minind) : -1 ;//commented for 76X

    //prescales
    if(triggerNames_[itrig].find("HLT_DiCentralPFJet170_v") != std::string::npos)   _pswgt_dijet_170 = preL1 * preHLT;
    if(triggerNames_[itrig].find("HLT_SingleCentralPFJet170_CFMax0p1_v") != std::string::npos)   _pswgt_singlejet_170_0p1 = preL1 * preHLT;

  }
  }

  _tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
TreeProducer_AOD::beginJob()
{
  // Initialize when class is created
  edm::Service<TFileService> fs ;
  _tree = fs->make <TTree>("SimpAnalysis","tree");

  // Declare tree's branches
  //_tree->Branch("",&,"");
  //
  // Event
  _tree->Branch("nEvent",&_nEvent,"nEvent/I");
  _tree->Branch("nRun",&_nRun,"nRun/I");
  _tree->Branch("nLumi",&_nLumi,"nLumi/I");
  //
  // Vertices
  _tree->Branch("vtx_N",&_vtx_N,"vtx_N/I");
  _tree->Branch("vtx_N_stored",&_vtx_N_stored,"vtx_N_stored/I");
  _tree->Branch("vtx_normalizedChi2",&_vtx_normalizedChi2,"vtx_normalizedChi2[vtx_N_stored]/D");
  _tree->Branch("vtx_ndof",&_vtx_ndof,"vtx_ndof[vtx_N_stored]/I");
  _tree->Branch("vtx_nTracks",&_vtx_nTracks,"vtx_nTracks[vtx_N_stored]/I");
  _tree->Branch("vtx_d0",&_vtx_d0,"vtx_d0[vtx_N_stored]/D");
  _tree->Branch("vtx_x",&_vtx_x,"vtx_x[vtx_N_stored]/D");
  _tree->Branch("vtx_y",&_vtx_y,"vtx_y[vtx_N_stored]/D");
  _tree->Branch("vtx_z",&_vtx_z,"vtx_z[vtx_N_stored]/D");
  //
	// Tracks
	_tree->Branch("nTrack_stored",&_nTrack_stored,"nTrack_stored/I");
	_tree->Branch("nTrack",&_nTrack,"nTrack/I");

	_tree->Branch("track_pt",&_track_pt,"track_pt[nTrack_stored]/D");
	_tree->Branch("track_eta",&_track_eta,"track_eta[nTrack_stored]/D");
	_tree->Branch("track_phi",&_track_phi,"track_phi[nTrack_stored]/D");

	_tree->Branch("track_normalizedChi2",&_track_normalizedChi2,"track_normalizedChi2[nTrack_stored]/D");
	_tree->Branch("track_ndof",&_track_ndof,"track_ndof[nTrack_stored]/I");
	_tree->Branch("track_ptError",&_track_ptError,"track_ptError[nTrack_stored]/D");
	_tree->Branch("track_dzError",&_track_dzError,"track_dzError[nTrack_stored]/D");
	_tree->Branch("track_dz",&_track_dz,"track_dz[nTrack_stored]/D");
	_tree->Branch("track_fromPV",&_track_fromPV,"track_fromPV[nTrack_stored]/I");
	_tree->Branch("track_purity",&_track_purity,"track_purity[nTrack_stored]/I");
	_tree->Branch("track_nhits",&_track_Nhits,"track_nhits[nTrack_stored]/I");
	_tree->Branch("track_nPixHits",&_track_NpixHits,"track_nPixHits[nTrack_stored]/I");
	_tree->Branch("track_d0",&_track_d0,"track_d0[nTrack_stored]/D");
	_tree->Branch("track_dxy",&_track_dxy,"track_dxy[nTrack_stored]/D");

  // GenJets
  _tree->Branch("nGenJet_stored",&_nGenJet_stored,"nGenJet_stored/I");
  _tree->Branch("nGenJet",&_nGenJet,"nGenJet/I");
  //
  _tree->Branch("genjetArea",&_genjet_area,"genjetArea[nGenJet_stored]/D");
  //
  _tree->Branch("genjet_vx",&_genjet_vx,"genjet_vx[nGenJet_stored]/D");
  _tree->Branch("genjet_vy",&_genjet_vy,"genjet_vy[nGenJet_stored]/D");
  _tree->Branch("genjet_vz",&_genjet_vz,"genjet_vz[nGenJet_stored]/D");
  //
  _tree->Branch("genjet_eta",&_genjet_eta,"genjet_eta[nGenJet_stored]/D");
  _tree->Branch("genjet_phi",&_genjet_phi,"genjet_phi[nGenJet_stored]/D");
  _tree->Branch("genjet_pt",&_genjet_pt,"genjet_pt[nGenJet_stored]/D");
  _tree->Branch("genjet_e",&_genjet_e,"genjet_e[nGenJet_stored]/D");
  _tree->Branch("genjet_m",&_genjet_m,"genjet_m[nGenJet_stored]/D");
  //
  _tree->Branch("genjet_efrac_ch", &_genjet_efrac_ch, "genjet_efrac_ch[nGenJet_stored]/D");

  //Trigger
  _tree->Branch("HLT_DiCentralPFJet170_CFMax0p1", &_dijet_170_0p1);
  _tree->Branch("HLT_DiCentralPFJet220_CFMax0p3", &_dijet_220_0p3);
  _tree->Branch("HLT_DiCentralPFJet330_CFMax0p5", &_dijet_330_0p5);
  _tree->Branch("HLT_DiCentralPFJet430", &_dijet_430);
  _tree->Branch("HLT_DiCentralPFJet170", &_dijet_170);
  _tree->Branch("HLT_SingleCentralPFJet170_CFMax0p1", &_singlejet_170_0p1);
  _tree->Branch("HLT_PFJet500", &_singlejet_500);
  _tree->Branch("HLT_PFJet450", &_singlejet_450);
  _tree->Branch("HLT_IsoMu24", &_isomu24);
  _tree->Branch("HLT_IsoMu27", &_isomu27);

  //prescales
  _tree->Branch("pswgt_dijet_170", &_pswgt_dijet_170, "pswgt_dijet_170/D");
  _tree->Branch("pswgt_singlejet_170_0p1", &_pswgt_singlejet_170_0p1, "pswgt_singlejet_170_0p1/D");

  //MET filters
  _tree->Branch("Flag_HBHENoiseFilter", &_HBHENoiseFlag);
  _tree->Branch("Flag_HBHENoiseIsoFilter", &_HBHENoiseIsoFlag);
  _tree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &_ECALFlag);
  _tree->Branch("Flag_goodVertices", &_vertexFlag);
  _tree->Branch("Flag_eeBadScFilter", &_eeFlag);
  _tree->Branch("Flag_globalTightHalo2016Filter", &_beamhaloFlag);
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


  triggerPathsVector.push_back("HLT_DiCentralPFJet170_CFMax0p1_v");
  triggerPathsVector.push_back("HLT_DiCentralPFJet220_CFMax0p3_v");
  triggerPathsVector.push_back("HLT_DiCentralPFJet330_CFMax0p5_v");
  triggerPathsVector.push_back("HLT_DiCentralPFJet430_v");
  triggerPathsVector.push_back("HLT_DiCentralPFJet170_v");
  triggerPathsVector.push_back("HLT_SingleCentralPFJet170_CFMax0p1_v");
  triggerPathsVector.push_back("HLT_PFJet450_v");
  triggerPathsVector.push_back("HLT_PFJet500_v");
  triggerPathsVector.push_back("HLT_IsoMu24_v");
  triggerPathsVector.push_back("HLT_IsoMu27_v");

  bool changedConfig = false;
  hltConfig_.init(iRun, iSetup, _trigResultsTag.process(), changedConfig);

  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    triggerPathsMap[triggerPathsVector[i]] = -1;
  }

  for(size_t i = 0; i < triggerPathsVector.size(); i++){
    TPRegexp pattern(triggerPathsVector[i]);
    for(size_t j = 0; j < hltConfig_.triggerNames().size(); j++){
      std::string pathName = hltConfig_.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
        triggerPathsMap[triggerPathsVector[i]] = j;
      }
    }
  }


  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,"HLT",changed) && hltPrescale_.init(iRun, iSetup, "HLT", changed) ) {
    if (changed) {
      // check if trigger names in (new) config
      cout<<"New trigger menu found !!!"<<endl;
      triggerIndex_.clear();
      const unsigned int n(hltConfig_.size());
      for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
        triggerIndex_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
        cout<<triggerNames_[itrig]<<" "<<triggerIndex_[itrig]<<" ";
        if (triggerIndex_[itrig] >= n)
          cout<<"does not exist in the current menu"<<endl;
        else
          cout<<"exists"<<endl;
      }
      cout << "Available TriggerNames are: " << endl;
      if (true)
        hltConfig_.dump("Triggers");
    }
  }
  else {
    cout << "ProcessedTreeProducerBTag::analyze:"
         << " config extraction failure with process name "
         << "HLT" << endl;
}

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

  //Trigger
  _dijet_170_0p1 = 0;
  _dijet_220_0p3 = 0;
  _dijet_330_0p5 = 0;
  _dijet_430 = 0;
  _dijet_170 = 0;
  _singlejet_170_0p1 = 0;
  _singlejet_450 = 0;
  _singlejet_500 = 0;
  _isomu24 = 0;
  _isomu27 = 0;
  //prescales
  _pswgt_dijet_170 = 1;
  _pswgt_singlejet_170_0p1 = 1;

  //MET filters
  _HBHENoiseFlag = 0;
  _HBHENoiseIsoFlag = 0;
  _ECALFlag = 0;
  _vertexFlag = 0;
  _eeFlag = 0;
  _beamhaloFlag = 0;


  // Vertices
  _vtx_N = 0;
  _vtx_N_stored = 0;
  //for(UInt_t iv=0;iv<nV;iv++) {
    _vtx_normalizedChi2.clear();
    _vtx_ndof.clear();
    _vtx_nTracks.clear();
    _vtx_d0.clear();
    _vtx_x.clear();
    _vtx_y.clear();
    _vtx_z.clear();
      //}

  //Tracks
  _nTrack = 0;
	_nTrack_stored = 0;
  //for(UInt_t it=0;it<nT;it++) {
		_track_eta.clear();
		_track_fromPV.clear();
		_track_ndof.clear();
		_track_Nhits.clear();
		_track_normalizedChi2.clear();
		_track_NpixHits.clear();
		_track_phi.clear();
		_track_pt.clear();
		_track_ptError.clear();
		_track_dzError.clear();
		_track_dz.clear();
		_track_d0.clear();
		_track_dxy.clear();
		_track_purity.clear();
	//}

  //GenJets
  _nGenJet = 0;
	_nGenJet_stored = 0;
  //for(UInt_t i=0 ; i<nGJ ; i++) {
    _genjet_vx.clear();
    _genjet_vy.clear();
    _genjet_vz.clear();
    _genjet_area.clear();
    _genjet_eta.clear();
    _genjet_phi.clear();
    _genjet_pt.clear();
    _genjet_e.clear();
    _genjet_m.clear();
    _genjet_efrac_ch.clear();
	//}

}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeProducer_AOD);
