#include "SexaQAnalysis/Data_analysis/interface/Analyzer.h"
#include <typeinfo>


Analyzer::Analyzer(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  m_bsTag    (pset.getParameter<edm::InputTag>("beamspot")),
  m_vertexTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  m_rCandsTag(pset.getParameter<edm::InputTag>("resonCandidates")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_vertexToken(consumes<vector<reco::Vertex> >(m_vertexTag)),
  //m_rCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_rCandsTag)),
  //m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag))
  m_rCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_rCandsTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_sCandsTag))
{
}


void Analyzer::beginJob() {
  histos_th1f["nPV"]      = m_fs->make<TH1F>("nPV",     "Number of PV; #PVs",100,0.,100.);
  histos_th1f["n_sCand"]  = m_fs->make<TH1F>("n_sCand",     "Number of sCands; #sCands",100,0.,100.);
  histos_th1f["n_rCand"]  = m_fs->make<TH1F>("n_rCand",     "Number of rCands; #rCands",100,0.,100.);
  histos_th2f["vtx_rz"]   = m_fs->make<TH2F>("vtx_rz",   "vtx_rz",800,-200.,200.,200,0.,100.);
  histos_th2f["vtx_xy"]   = m_fs->make<TH2F>("vtx_xy",   "vtx_xy",400,-100.,100.,400,-100.,100.);
  histos_th1f["vtx_chi2"] = m_fs->make<TH1F>("vtx_chi2", "vtx_chi2",100,0.,10.);
  histos_th1f["dxy"]      = m_fs->make<TH1F>("dxy",      "dxy",200,0.,50.);
  histos_th1f["dr"]       = m_fs->make<TH1F>("dr",       "dr",200,0.,50.);
  histos_th1f["sdxy"]     = m_fs->make<TH1F>("sdxy",     "sdxy",100,0.,100.);
  histos_th1f["sdr"]      = m_fs->make<TH1F>("sdr",      "sdr",100,0.,100.);
  
  //nr L0 and Ks per PV TProfiles
  
  //histos_TProfile["nL0_per_PV_vs_nPV"] = m_fs->make<TProfile>("nL0_per_PV_vs_nPV","nL0_per_PV_vs_nPV; nPV; nL0 per PV",30, 0.,30., 1000, 0, 1);
  //histos_TProfile["nKs_per_PV_vs_nPV"] = m_fs->make<TProfile>("nKs_per_PV_vs_nPV","nKs_per_PV_vs_nPV; nPV; nKs per PV",30, 0.,30., 1000, 0, 1);
 
  histos_TProfile["nL0Ks_per_PV_vs_nPV"] = m_fs->make<TProfile>("nL0Ks_per_PV_vs_nPV","nL0Ks_per_PV_vs_nPV; nPV; numberof L0 Ks fits per PV",45, 0.,45., 0., 0.25);
  histos_TProfile["nL0Ks_vs_nPV"] = m_fs->make<TProfile>("nL0Ks_vs_nPV","nL0Ks_vs_nPV; nPV; number of L0 Ks fits",45, 0.,45., 0., 10.);


  
  //MASSES
  
  histos_th1f["rCandMass"] = m_fs->make<TH1F>("rCandMass","rCandMass; mass (GeV)",1000,-10.,10.);
  histos_th1f["sCandMass"] = m_fs->make<TH1F>("sCandMass","sCandMass; mass (GeV)",1000,-10.,10.);
  
  //rCAND DISTIBUTIONS
  
  histos_th1f["rCand_mass_below_2GeV_Eta"] = m_fs->make<TH1F>("rCand_mass_below_2GeV_Eta","rCand_mass_below_2GeV_Eta",1000,-5.,5.);
  
  //sCAND DISTIBUTIONS
  
  histos_th1f["sCand_delta_phi"] = m_fs->make<TH1F>("sCand_delta_phi","sCand_delta_phi; delta Phi",1000,-5.,5.);
  histos_th1f["sCand_delta_R"] = m_fs->make<TH1F>("sCand_delta_R","sCand_delta_R; delta R",1000,0,10.);
    histos_th1f["sCand_delta_eta"] = m_fs->make<TH1F>("sCand_delta_eta","sCand_delta_eta; delta eta",1000,-5.,5.);



  
  //SCATTERPLOTS
  
  histos_th2f["scatterplot_sCand_pos"]    = m_fs->make<TH2F>("scatterplot_sCand_pos", "scatterplot_sCand_pos",1000,-100.,100.,1000,-100.,100.);
  histos_th2f["scatterplot_bs_pos"]    = m_fs->make<TH2F>("scatterplot_bs_pos", "scatterplot_bs_pos",1000,-100.,100.,1000,-100.,100.);
  histos_th2f["scatterplot_PCA_pos"]    = m_fs->make<TH2F>("scatterplot_PCA_pos", "scatterplot_PCA_pos",1000,-100.,100.,1000,-100.,100.);
	
  //DELTA PHI INVESTIGATION
  
  histos_th2f["sCand_delta_phi_dxy(PCA_beamspot)_signed"]= m_fs->make<TH2F>("sCand_delta_phi_dxy(PCA_beamspot)_signed","sCand_delta_phi_dxy(PCA_beamspot)_signed; delta Phi; distance (cm)",1000,-5.,5., 1000, 0, 20);
  histos_th2f["sCand_delta_phi_dxy(sCandVtx_beamspot)_signed"]= m_fs->make<TH2F>("sCand_delta_phi_dxy(sCandVtx_beamspot)_signed","sCand_delta_phi_dxy(sCandVtx_beamspot)_signed; delta Phi; distance (cm)",1000,-5.,5., 1000, 0, 20);
  
  //Vxy PLOTS
  
  histos_th1f["data_sCand_dxy(sCandVtx_beamspot)_signed_L0"] = m_fs->make<TH1F>("data_sCand_dxy(sCandVtx_beamspot)_signed_L0","data_sCand_dxy(sCandVtx_beamspot)_signed_L0; distance (cm)",1000,-2.5,2.5);
  histos_th1f["data_sCand_dxy(sCandVtx_beamspot)_signed_anti_L0"] = m_fs->make<TH1F>("data_sCand_dxy(sCandVtx_beamspot)_signed_anti_L0","data_sCand_dxy(sCandVtx_beamspot)_signed_anti_L0; distance (cm)",1000,-2.5,2.5);

  histos_th1f["data_sCand_dxy(sCandVtx_beamspot)_signed_L0_deltaPhiCut"] = m_fs->make<TH1F>("data_sCand_dxy(sCandVtx_beamspot)_signed_L0_deltaPhiCut","data_sCand_dxy(sCandVtx_beamspot)_signed_L0_deltaPhiCut; distance (cm)",1000,-2.5,2.5);
  histos_th1f["data_sCand_dxy(sCandVtx_beamspot)_signed_anti_L0_deltaPhiCut"] = m_fs->make<TH1F>("data_sCand_dxy(sCandVtx_beamspot)_signed_anti_L0_deltaPhiCut","data_sCand_dxy(sCandVtx_beamspot)_signed_anti_L0_deltaPhiCut; distance (cm)",1000,-2.5,2.5);

  //dxy PCA PLOTS
  
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0; distance (cm)",1000,-0.5,0.5);

  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut; distance (cm)",1000,-0.5,0.5);
  
  //dxy PCA PLOTS with negative variance cut
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_PosVar"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_PosVar","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_PosVar; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_PosVar"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_PosVar","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_PosVar; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_PosVar"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_PosVar","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_PosVar; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_PosVar"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_PosVar","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_PosVar; distance (cm)",1000,-0.5,0.5);

  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut_PosVar"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut_PosVar","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut_PosVar; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut_PosVar"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut_PosVar","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut_PosVar; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut_PosVar"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut_PosVar","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut_PosVar; distance (cm)",1000,-0.5,0.5);
  histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut_PosVar"] = m_fs->make<TH1F>("data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut_PosVar","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut_PosVar; distance (cm)",1000,-0.5,0.5);

  //dxy PCA significance PLOTS (also with negative variance cut)
  histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0"] = m_fs->make<TH2F>("data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 10^7);
  histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0"] = m_fs->make<TH2F>("data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0","data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 10^7);
  histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0"] = m_fs->make<TH2F>("data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 10^7);
  histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0"] = m_fs->make<TH2F>("data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0","data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 10^7);

  histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut"] = m_fs->make<TH2F>("data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut","data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 10^7);
  histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH2F>("data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut","data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 10^7);
  histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut"] = m_fs->make<TH2F>("data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut","data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 10^7);
  histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut"] = m_fs->make<TH2F>("data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut","data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut; distance (cm);significance (dxy / sigma_dxy)",1000,-0.5,0.5, 1000, 0, 10^7);

}


void Analyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {

  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

  edm::Handle<vector<reco::Vertex> > h_vert; //Primary vertices
  iEvent.getByToken(m_vertexToken, h_vert);

  // resonance candidates
  //edm::Handle<vector<reco::VertexCompositeCandidate> > h_rCands;
  edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_rCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
  iEvent.getByToken(m_rCandsToken, h_rCands);

  // sexaquark candidates
  //edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_sCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
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
  Int_t nPV = h_vert->size();
  Int_t n_sCands = h_sCands->size();
  Int_t n_rCands = h_rCands->size();
 

  //std::cout<<m_nRun<<"\t"<<m_nLumi<<"\t"<<m_nEvent<<std::endl;

  histos_th1f["nPV"]->Fill(nPV);
  histos_th1f["n_sCand"]->Fill(h_sCands->size());    
  histos_th1f["n_rCand"]->Fill(h_rCands->size());
  if (h_vert->size() == 0) return;  //only look at events with reconstructed L0-Ks vertices
  
  histos_TProfile["nL0Ks_per_PV_vs_nPV"] ->Fill(nPV, Float_t(h_rCands->size()) / Float_t(nPV));
  histos_TProfile["nL0Ks_vs_nPV"]  ->Fill(nPV, Float_t(h_rCands->size()));
  
  //cout<<Float_t(h_rCands->size()) / Float_t(nPV)<<endl;
  //cout<<h_rCands->size()<<endl;


  //beamspot position
  float bs_x = h_bs->x0(); 
  float bs_y = h_bs->y0();
  float bs_z = h_bs->z0();
  TVector3 bs_pos(bs_x,bs_y,bs_z); 
  TVector2 bs_pos2D(bs_x,bs_y);    
  histos_th2f["scatterplot_bs_pos"]->Fill(bs_x,bs_y);


  for (unsigned int i = 0; i < h_rCands->size(); ++i) { //loop over resonance candidates
    float x  = h_rCands->at(i).vx();
    float y  = h_rCands->at(i).vy();
    float z  = h_rCands->at(i).vz();
    float xe = h_rCands->at(i).vertexCovariance(0,0);
    float ye = h_rCands->at(i).vertexCovariance(1,1);
    float ze = h_rCands->at(i).vertexCovariance(2,2);
    
    //cout<<"Cov rCand: "<<h_rCands->at(i).vertexCovariance(0,0)<<" "<<h_rCands->at(i).vertexCovariance(0,1)<<" "<<h_rCands->at(i).vertexCovariance(1,0)<<" "<<h_rCands->at(i).vertexCovariance(1,1)<<endl;
    histos_th2f["vtx_xy"]  ->Fill(x,y);
    histos_th2f["vtx_rz"]  ->Fill(z,sqrt(pow(x,2)+pow(y,2)));
//    histos_th1f["vtx_chi2"]->Fill(10*sqrt(xe+ye));
    histos_th1f["dxy"] ->Fill(sqrt(pow(x-bs_x,2)+pow(y-bs_y,2)));
    histos_th1f["dr"]  ->Fill(sqrt(pow(x-bs_x,2)+pow(y-bs_y,2)+pow(z-bs_z,2)));
    histos_th1f["sdxy"]->Fill(sqrt((pow(x-bs_x,2)+pow(y-bs_y,2))/(xe+ye)));
    histos_th1f["sdr"] ->Fill(sqrt((pow(x-bs_x,2)+pow(y-bs_y,2)+pow(z-bs_z,2))/(xe+ye+ze)));
    
    if(h_rCands->at(i).mass() < 2.) histos_th1f["rCand_mass_below_2GeV_Eta"]->Fill(h_rCands->at(i).eta()); //to be compared to gen Xi(1530) eta distri

    float dxy = sqrt(pow(x-bs_x,2)+pow(y-bs_y,2));
    float edxy = sqrt(xe+ye);
    if (dxy < 0.1 && 
        dxy < 3*edxy &&
        edxy < .1 &&
        sqrt(ze) < .1)
      histos_th1f["rCandMass"]->Fill(h_rCands->at(i).mass());
  }
  

  for (unsigned int i = 0; i < h_sCands->size(); ++i) { //loop over S candidates
	  
     //cout <<h_sCands->at(i).vertexChi2() <<endl;

	//sCand vertex position and errors
    float x  = h_sCands->at(i).vx();
    float y  = h_sCands->at(i).vy();
    //float z  = h_sCands->at(i).vz();
    TVector2 cand_pos(x,y);

    float xe = h_sCands->at(i).vertexCovariance(0,0);
    float ye = h_sCands->at(i).vertexCovariance(1,1);
    float ze = h_sCands->at(i).vertexCovariance(2,2);
    
    TMatrixD CovMx2D(2,2); //2D xy covariance matrix
    CovMx2D(0,0)=h_sCands->at(i).vertexCovariance(0,0); //elements are Sigma(i,j)²
    CovMx2D(0,1)=h_sCands->at(i).vertexCovariance(0,1);
    CovMx2D(1,0)=h_sCands->at(i).vertexCovariance(1,0);
    CovMx2D(1,1)=h_sCands->at(i).vertexCovariance(1,1);
    
    
    
    //CovMx2D.Print();
    
    
     
    //calculate eigenvalues and eigenvectors
    //used for significance in PCA plots
    //http://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/
    TVectorD EValues;
    TMatrixD EVectors = CovMx2D.EigenVectors(EValues);
    
    //if(!(CovMx2D(1,1)<0 || CovMx2D(0,0)<0)){ //cut negative sigma_y^2 variances
		
	for (Int_t i = 0; i < EValues.GetNrows(); ++i) { 
		
		TVectorD EVector(TMatrixTColumn_const<double>(EVectors, i)); 
		//cout << "eigen-value " << i << " is " << EValues(i);
		//cout << " with eigen-vector "; 	EVector.Print(); 
		
	}
	//cout<<endl;
	//}
	
	
	float MinEValue = ( EValues(0)<EValues(1) ) ? EValues(0) : EValues(1);
	//cout <<"MinEValue: "<<MinEValue<<endl;
	//if(!(CovMx2D(1,1)<0 || CovMx2D(0,0)<0)) cout <<"EValues: "<<EValues(0)<<" "<<EValues(1)<<" min: "<<MinEValue<<endl;

    
    

    //momentum
    float px = h_sCands->at(i).px(); 
    float py = h_sCands->at(i).py();
    //float pz = h_sCands->at(i).pz();
    TVector2 cand_mom(px,py);

    
    //to distinguish lambda from anti-lambda, one can use this info:
    //std::cout << "Lambda = 1115; "<<h_sCands->at(i).daughterPtr(0)->mass() << " "
              //<< "Proton = 938; "<<h_sCands->at(i).daughterPtr(0)->daughterPtr(0)->mass() << " "
              //<< h_sCands->at(i).daughterPtr(0)->daughterPtr(0)->charge()
             // << std::endl;
 

    float dxy = sqrt(pow(x,2)+pow(y,2));
    float edxy = sqrt(xe+ye);
    if (dxy > 2-3*edxy && edxy < .1 && sqrt(ze) < .1){
      histos_th1f["sCandMass"]->Fill(h_sCands->at(i).mass());
	}
	

    //Extrapolate vertex momentum in direction of the unit momentum vector to the Point of Closest Approach with the beamspot

    float alpha = ((cand_mom * bs_pos2D) - (cand_mom * cand_pos)) / (cand_mom * cand_mom);

    TVector2 PCA = cand_pos + (alpha * cand_mom);

    TVector2 PCA_bs = PCA - bs_pos2D;

    float dxy_PCA_bs = sqrt(pow(PCA_bs.X(),2)+pow(PCA_bs.Y(),2));

    TVector2 cand_bs = cand_pos - bs_pos2D; //vector between sCand and beamspot
    
    float dxy_cand_bs = sqrt(pow(cand_bs.X(),2)+pow(cand_bs.Y(),2));
    

    float dxy_PCA_bs_signed = dxy_PCA_bs * signum(cand_bs * cand_mom); //2D dot product

    //cout << signum(cand_bs * cand_mom) << endl;
    
    
    
    cout << h_sCands->at(i).daughterPtr(0)->Eta() << endl;
    //int proton_charge = h_sCands->at(i).daughterPtr(0)->daughter(0)->charge(); //This granddaughterPtr of the sCand is the proton.
                                        //in the case the proton charge is (-1)+1, the sCand daughterPtr includes an (anti)Lambda
    
    /*
    histos_th2f["scatterplot_sCand_pos"] ->Fill(x,y);

	histos_th2f["scatterplot_PCA_pos"]   ->Fill(PCA.X(), PCA.Y());
	
	float delta_phi = reco::deltaPhi(h_sCands->at(i).daughterPtr(0)->phi(), h_sCands->at(i).daughterPtr(1)->phi());

	float delta_eta = h_sCands->at(i).daughterPtr(0)->eta() - h_sCands->at(i).daughterPtr(1)->eta();
	
	histos_th1f["sCand_delta_phi"]->Fill(delta_phi);
	histos_th1f["sCand_delta_R"]->Fill(deltaR(h_sCands->at(i).daughterPtr(0)->eta(), h_sCands->at(i).daughterPtr(1)->eta(), h_sCands->at(i).daughterPtr(0)->phi(), h_sCands->at(i).daughterPtr(1)->phi()));
	histos_th1f["sCand_delta_eta"]->Fill(delta_eta);
	
	
	
	histos_th2f["sCand_delta_phi_dxy(PCA_beamspot)_signed"]->Fill(delta_phi, dxy_PCA_bs_signed);
	histos_th2f["sCand_delta_phi_dxy(sCandVtx_beamspot)_signed"]->Fill(delta_phi, dxy_cand_bs * signum(cand_bs * cand_mom));

    
    
    ///PCA AND SIGNIFICANCE PLOTS
    
    float dxy_signif_PCA_bs= abs(dxy_PCA_bs_signed)/sqrt(abs(MinEValue)); //significance = dxy / sigma_dxy ≃ dxy / sqrt(min(eigenvalues(sCand vtx 2D Cov matrix)))
    float dxy_signif_PCA_bs_signed = dxy_PCA_bs_signed/sqrt(abs(MinEValue)); //significance = dxy / sigma_dxy ≃ dxy / sqrt(min(eigenvalues(sCand vtx 2D Cov matrix)))
    //cut away instances when the x or y variance of the sCand position is negative
    bool posVar = true;
    if(CovMx2D(1,1)<0 || CovMx2D(0,0)<0){cout<<"neg variance"<<endl; posVar = false; }//negative variance is bad
    
    
    
    
    if(0.5<dxy_cand_bs && dxy_cand_bs<1.8) {
		
      if(proton_charge == 1) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0"] ->Fill(dxy_PCA_bs_signed);
      
      if(proton_charge == 1 && posVar) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_PosVar"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1 && posVar) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_PosVar"] ->Fill(dxy_PCA_bs_signed);
      
	  if(proton_charge == 1 && posVar) histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs); 
      if(proton_charge == -1 && posVar) histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
    
    } else if(1.8 < dxy_cand_bs) {
		
      if(proton_charge == 1) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0"] ->Fill(dxy_PCA_bs_signed);
      
      if(proton_charge == 1 && posVar) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_PosVar"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1 && posVar) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_PosVar"] ->Fill(dxy_PCA_bs_signed);
      
      if(proton_charge == 1 && posVar) histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
      if(proton_charge == -1 && posVar) histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
    }
    
    if(proton_charge == 1) histos_th1f["data_sCand_dxy(sCandVtx_beamspot)_signed_L0"]->Fill(dxy_cand_bs * signum(cand_bs * cand_mom));
    if(proton_charge == -1) histos_th1f["data_sCand_dxy(sCandVtx_beamspot)_signed_anti_L0"]->Fill(dxy_cand_bs * signum(cand_bs * cand_mom));
    
    
    //Same plots but with a delta phi cut
    if(abs(delta_phi)<0.5 || abs(delta_phi) > (TMath::Pi() - 0.5)) continue;
    
	if(0.5<dxy_cand_bs && dxy_cand_bs<1.8) {
		
      if(proton_charge == 1) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed);
      
      if(proton_charge == 1 && posVar) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut_PosVar"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1 && posVar) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut_PosVar"] ->Fill(dxy_PCA_bs_signed);
      
      if(proton_charge == 1 && posVar) histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
      if(proton_charge == -1 && posVar) histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_between_5mm_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
    
    } else if(1.8 < dxy_cand_bs) {
		
      if(proton_charge == 1) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed);
      
      if(proton_charge == 1 && posVar) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut_PosVar"] ->Fill(dxy_PCA_bs_signed);
      if(proton_charge == -1 && posVar) histos_th1f["data_sCand_dxy(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut_PosVar"] ->Fill(dxy_PCA_bs_signed);
      
      if(proton_charge == 1 && posVar) histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
      if(proton_charge == -1 && posVar) histos_th2f["data_sCand_dxy_signif_(PCA_beamspot)_Vxy_over_18mm_anti_L0_deltaPhiCut"] ->Fill(dxy_PCA_bs_signed, dxy_signif_PCA_bs);
    }
    
    if(proton_charge == 1) histos_th1f["data_sCand_dxy(sCandVtx_beamspot)_signed_L0_deltaPhiCut"]->Fill(dxy_cand_bs * signum(cand_bs * cand_mom));
    if(proton_charge == -1) histos_th1f["data_sCand_dxy(sCandVtx_beamspot)_signed_anti_L0_deltaPhiCut"]->Fill(dxy_cand_bs * signum(cand_bs * cand_mom));
*/
  }
  



}


void Analyzer::endJob()
{

}

Analyzer::~Analyzer()
{
}


DEFINE_FWK_MODULE(Analyzer);
