#include "SexaQAnalysis/TreeProducer/interface/Analyzer.h"


Analyzer::Analyzer(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  m_bsTag    (pset.getParameter<edm::InputTag>("beamspot")),
  m_vertexTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  m_rCandsTag(pset.getParameter<edm::InputTag>("resonCandidates")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_vertexToken(consumes<vector<reco::Vertex> >(m_vertexTag)),
  m_rCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_rCandsTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_sCandsTag))
{
}


void Analyzer::beginJob() {
  histos_th1f["nVtx"]      = m_fs->make<TH1F>("nVtx",     "",100,0.,100.);
  histos_th2f["vtx_rz"]    = m_fs->make<TH2F>("vtx_rz",   "",800,-200.,200.,200,0.,100.);
  histos_th2f["vtx_xy"]    = m_fs->make<TH2F>("vtx_xy",   "",400,-100.,100.,400,-100.,100.);
  histos_th1f["vtx_chi2"]  = m_fs->make<TH1F>("vtx_chi2", "",100,0.,10.);
  histos_th1f["dxy"]       = m_fs->make<TH1F>("dxy",      "",200,0.,50.);
  histos_th1f["dr"]        = m_fs->make<TH1F>("dr",       "",200,0.,50.);
  histos_th1f["sdxy"]      = m_fs->make<TH1F>("sdxy",     "",100,0.,100.);
  histos_th1f["sdr"]       = m_fs->make<TH1F>("sdr",      "",100,0.,100.);
  histos_th1f["rCandMass"] = m_fs->make<TH1F>("rCandMass","",1000,0.,5.);
  histos_th1f["sCandMass"] = m_fs->make<TH1F>("sCandMass","",1000,0.,5.);
}


void Analyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

  edm::Handle<vector<reco::Vertex> > h_vert;
  iEvent.getByToken(m_vertexToken, h_vert);

  // resonance candidates
  edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_rCands;
  iEvent.getByToken(m_rCandsToken, h_rCands);

  // sexaquark candidates
  edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);

  // Check validity
  if(!h_bs.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_bsTag << " ... skip entry !" << endl;
    return;
  }

  if(!h_vert.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_vertexTag << " ... skip entry !" << endl;
    return;
  }

  if(!h_rCands.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_rCandsTag << " ... skip entry !" << endl;
    return;
  }

  if(!h_sCands.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_sCandsTag << " ... skip entry !" << endl;
    return;
  }

  // GLOBAL EVENT INFORMATIONS //
  m_nRun   = iEvent.id().run();
  m_nLumi  = iEvent.luminosityBlock();
  m_nEvent = iEvent.id().event();

  //std::cout<<m_nRun<<"\t"<<m_nLumi<<"\t"<<m_nEvent<<std::endl;

  histos_th1f["nVtx"]->Fill(h_vert->size());
  if (h_vert->size() == 0) return;
  float vx = h_bs->x0();
  float vy = h_bs->y0();
  float vz = h_bs->z0();

  for (unsigned int i = 0; i < h_rCands->size(); ++i) {
    float x  = h_rCands->at(i).vx();
    float y  = h_rCands->at(i).vy();
    float z  = h_rCands->at(i).vz();
    float xe = h_rCands->at(i).vertexCovariance(0,0);
    float ye = h_rCands->at(i).vertexCovariance(1,1);
    float ze = h_rCands->at(i).vertexCovariance(2,2);
    histos_th2f["vtx_xy"]  ->Fill(x,y);
    histos_th2f["vtx_rz"]  ->Fill(z,sqrt(pow(x,2)+pow(y,2)));
//    histos_th1f["vtx_chi2"]->Fill(10*sqrt(xe+ye));
    histos_th1f["dxy"] ->Fill(sqrt(pow(x-vx,2)+pow(y-vy,2)));
    histos_th1f["dr"]  ->Fill(sqrt(pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2)));
    histos_th1f["sdxy"]->Fill(sqrt((pow(x-vx,2)+pow(y-vy,2))/(xe+ye)));
    histos_th1f["sdr"] ->Fill(sqrt((pow(x-vx,2)+pow(y-vy,2)+pow(z-vz,2))/(xe+ye+ze)));

    float dxy = sqrt(pow(x-vx,2)+pow(y-vy,2));
    float edxy = sqrt(xe+ye);
    if (dxy < 0.1 && 
        dxy < 3*edxy &&
        edxy < .1 &&
        sqrt(ze) < .1)
      histos_th1f["rCandMass"]->Fill(h_rCands->at(i).mass());
  }

  for (unsigned int i = 0; i < h_sCands->size(); ++i) {
    float x  = h_sCands->at(i).vx();
    float y  = h_sCands->at(i).vy();
    float xe = h_sCands->at(i).vertexCovariance(0,0);
    float ye = h_sCands->at(i).vertexCovariance(1,1);
    float ze = h_sCands->at(i).vertexCovariance(2,2);

/* to distinguish lambda from anti-lambda, one can use this info:
    std::cout << h_sCands->at(i).daughterPtr(0)->mass() << " "
              << h_sCands->at(i).daughterPtr(0)->daughter(0)->mass() << " "
              << h_sCands->at(i).daughterPtr(0)->daughter(0)->charge()
              << std::endl;
*/

    float dxy = sqrt(pow(x,2)+pow(y,2));
    float edxy = sqrt(xe+ye);
    if (dxy > 2-3*edxy &&
        edxy < .1 &&
        sqrt(ze) < .1)
      histos_th1f["sCandMass"]->Fill(h_sCands->at(i).mass());
  }


}


void Analyzer::endJob()
{

}

Analyzer::~Analyzer()
{
}


DEFINE_FWK_MODULE(Analyzer);
