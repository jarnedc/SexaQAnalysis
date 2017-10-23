#include "HexaAnalysis/TreeProducer/interface/Analyzer.h"


Analyzer::Analyzer(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  m_vertexCollectionTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  m_trackCollectionTag(pset.getParameter<edm::InputTag>("trackCollection")),
  m_partonsTag(pset.getParameter<edm::InputTag>("genCollection")),
  m_vertexCollectionToken(consumes<vector<reco::Vertex> >(m_vertexCollectionTag)),
  m_trackCollectionToken(consumes<vector<reco::Track> >(m_trackCollectionTag)),
  m_partons(consumes<vector<reco::GenParticle> >(m_partonsTag))
{

  m_isData       = pset.getUntrackedParameter<bool>             ("isData");

}


void Analyzer::beginJob()
{
  // Initialize when class is created
  //_tree = fs->make <TTree>("HexaQAnalysis","tree");
  histos_th1f["num_of_Vtx"] = m_fs->make<TH1F>("num_of_Vtx","num_of_Vtx",100,0.,100.);
  histos_th1f["track_pt"] = m_fs->make<TH1F>("track_pt","track_pt",5000,0.,5000.);
  histos_th1f["gen_pt"] = m_fs->make<TH1F>("gen_pt","gen_pt",5000,0.,5000.);

}

void Analyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
{
  // HANDLES //
  // Get collections
  //edm::Handle<edm::TriggerResults> H_trig;//, H_trig1, H_trig2;
  //iEvent.getByToken(_trigResultsToken, H_trig);

  edm::Handle<vector<reco::Vertex> > H_vert;
  iEvent.getByToken(m_vertexCollectionToken, H_vert);

  edm::Handle<vector<reco::Track> > H_track;
  iEvent.getByToken(m_trackCollectionToken , H_track);

  edm::Handle<vector<reco::GenParticle> > H_partons;
  iEvent.getByToken(m_partons, H_partons);

  // Check validity
  if(!H_vert.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_vertexCollectionTag << " ... skip entry !" << endl;
    return;
  }

  if(!H_track.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_trackCollectionTag << " ... skip entry !" << endl;
    return;
  }

  if(!H_partons.isValid()) {
    if(verbose>0) cout << "Missing collection : "<< m_partonsTag <<"... skip entry !" << endl;
    return;
  }

  // GLOBAL EVENT INFORMATIONS //
  m_nRun   = iEvent.id().run();
  m_nLumi  = iEvent.luminosityBlock();
  m_nEvent = iEvent.id().event();

  //std::cout<<m_nRun<<"\t"<<m_nLumi<<"\t"<<m_nEvent<<std::endl;


  histos_th1f["num_of_Vtx"]->Fill(H_vert->size());


  // TRACKS //
  vector<reco::TrackRef> trackRef;
  for (vector<reco::Track>::const_iterator theTrack = H_track->begin(); theTrack != H_track->end(); ++theTrack){
          reco::TrackRef ref(H_track, theTrack - H_track->begin());
                  trackRef.push_back(ref);
  }

  reco::RecoPtSorter<reco::TrackRef> trackSorter;
  std::sort( trackRef.begin(), trackRef.end(), trackSorter);

  for (size_t i = 0; i < trackRef.size(); i++) {
          if(trackRef[i]->pt()<0.5) continue;
          histos_th1f["track_pt"]->Fill(trackRef[i]->pt());
  }

  if(!m_isData){
            for (vector<reco::GenParticle>::const_iterator thepartons = H_partons->begin();
            thepartons != H_partons->end(); ++thepartons){

                if(thepartons->status()!=2) continue;
                histos_th1f["gen_pt"]->Fill(thepartons->pt());

            }///for partons
  }


}


void Analyzer::endJob()
{

}

Analyzer::~Analyzer()
{
}


DEFINE_FWK_MODULE(Analyzer);
