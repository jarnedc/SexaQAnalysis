#include "SexaQAnalysis/V0_angular_correlation_analysis/interface/Analyzer_V0_angular_correlation.h"
#include <typeinfo>


Analyzer_V0_angular_correlation::Analyzer_V0_angular_correlation(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  m_bsTag    (pset.getParameter<edm::InputTag>("beamspot")),
  m_vertexTag(pset.getParameter<edm::InputTag>("vertexCollection")),
  //m_rCandsTag(pset.getParameter<edm::InputTag>("resonCandidates")),
  //m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_KshortsTag(pset.getParameter<edm::InputTag>("lambdaCollection")),
  m_LambdasTag(pset.getParameter<edm::InputTag>("kshortCollection")),
  m_bsToken    (consumes<reco::BeamSpot>(m_bsTag)),
  m_vertexToken(consumes<vector<reco::Vertex> >(m_vertexTag)),
  //m_rCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_rCandsTag)),
  //m_sCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_sCandsTag)),
  //m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_KshortsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_KshortsTag)),
  m_LambdasToken(consumes<vector<reco::VertexCompositeCandidate> >(m_LambdasTag))


{

}


void Analyzer_V0_angular_correlation::beginJob() {
  
 
  
//  histos_th1f[b+"nPV"]      = m_fs->make<TH1F>(b+"nPV",     a+" Number of PV; #PVs",60,0.,60);
  
  
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta"]= m_fs->make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta",b+"h_L0_Ks_delta_phi_delta_eta; delta phi; delta_eta",100,-4,4, 100, -4, 4); 
  
  
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]= m_fs->make<TH2F>(b+"h_L0_delta_phi_delta_eta_prev_and_current",b+"h_L0_delta_phi_delta_eta; delta phi prev event and current event; delta_eta previous event and current event",100,-4,4, 100, -4, 4); 
   
    

}


void Analyzer_V0_angular_correlation::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  edm::Handle<reco::BeamSpot> h_bs;
  iEvent.getByToken(m_bsToken, h_bs);

  edm::Handle<vector<reco::Vertex> > h_vert; //Primary vertices
  iEvent.getByToken(m_vertexToken, h_vert);


  // resonance candidates
  
  //edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_rCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
  //iEvent.getByToken(m_rCandsToken, h_rCands);

  // sexaquark candidates
  
  //edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_sCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
  //iEvent.getByToken(m_sCandsToken, h_sCands);
  
  //lambdaKshortVertexFilter sexaquark candidates
 // edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
 // iEvent.getByToken(m_sCandsToken, h_sCands);
  
  //reco Kshorts V0
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_Kshorts;
  iEvent.getByToken(m_KshortsToken, h_Kshorts);
  
  //reco Lambdas V0
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_Lambdas;
  iEvent.getByToken(m_LambdasToken, h_Lambdas);

  // Check validity
  if(!h_bs.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_bsTag << " ... skip entry !" << endl;
    return;
  }

  if(!h_vert.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_vertexTag << " ... skip entry !" << endl;
    return;
  }
/*
  if(!h_rCands.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_rCandsTag << " ... skip entry !" << endl;
    return;
  }
*/

 // if(!h_sCands.isValid()) {
 //   if(verbose>0) cout << "Missing collection : " << m_sCandsTag << " ... skip entry !" << endl;
 //   return;
 // }
  
  if(!h_Lambdas.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_LambdasTag << " ... skip entry !" << endl;
    return;
  }
  
  if(!h_Kshorts.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_KshortsTag << " ... skip entry !" << endl;
    return;
  }

   unsigned int n_Lambdas = h_Lambdas->size(); //Jarne: this is the original V0 collection
   unsigned int n_Kshorts = h_Kshorts->size();


 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STANDARD PLOT: ONLY THE CURRENT EVENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 for (unsigned int i = 0; i < n_Lambdas; ++i) { //loop over reco V0 Lambdas
	for (unsigned int i = 0; i < n_Kshorts; ++i) { //loop over reco Kshorts	
		double phi1 = h_Kshorts->at(i).phi();
		double phi2 = h_Lambdas->at(i).phi();
		double delta_phi = reco::deltaPhi(phi1, phi2);
		
		double eta1 = h_Kshorts->at(i).eta();
		double eta2 = h_Lambdas->at(i).eta();
		double delta_eta = eta1-eta2;

		histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta"]->Fill(delta_phi,delta_eta);
	}
 }


  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COMPARE THE PREVIOUS EVENT WITH THE CURRENT EVENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  bool emptyVectors = false;
  if(v_L0_phi.empty() && v_L0_eta.empty() && v_Ks_phi.empty() && v_Ks_eta.empty()) emptyVectors = true;//this means they have been cleared in the previous event or it is the first event

  for (unsigned int i = 0; i < n_Lambdas; ++i) { //loop over reco V0 Lambdas
  	if(emptyVectors){//just save the data if there is no data in these vectrors, make the histogram in the next event
		v_L0_phi.push_back(h_Lambdas->at(i).phi());	
		v_L0_eta.push_back(h_Lambdas->at(i).eta());	
	} 
	else{//now make the histo, also using the data from the previous event
		for(unsigned int i_prev = 0; i_prev < v_Ks_phi.size(); i_prev++){//loop over all the data from the previous event
			double phi1 = v_Ks_phi[i_prev];
			double phi2 = h_Lambdas->at(i).phi();
			double delta_phi = reco::deltaPhi(phi1, phi2);
			
			double eta1 = v_Ks_eta[i_prev];
			double eta2 = h_Lambdas->at(i).eta();
			double delta_eta = eta1-eta2;
			histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]->Fill(delta_phi,delta_eta);
		}
	} 
  }//end loop over lambda

 if(emptyVectors) std::cout << "vectors empty" << std::endl;  
 else std::cout << "vectors not empty" << std::endl;  
  
  for (unsigned int i = 0; i < n_Kshorts; ++i) { //loop over reco Kshorts
 	 if(emptyVectors){
                v_Ks_phi.push_back(h_Kshorts->at(i).phi());
                v_Ks_eta.push_back(h_Kshorts->at(i).eta());
         } 
 
	 else{//now make the histo, also using the data from the previous event
		for(unsigned int i_prev = 0; i_prev < v_L0_phi.size(); i_prev++){//loop over all the data from the previous event
			double phi1 = v_L0_phi[i_prev];
			double phi2 = h_Kshorts->at(i).phi();
			double delta_phi = reco::deltaPhi(phi1, phi2);
			
			double eta1 = v_L0_eta[i_prev];
			double eta2 = h_Kshorts->at(i).eta();
			double delta_eta = eta1-eta2;

			histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]->Fill(delta_phi,delta_eta);
		}
	}
  }// end loop over Kshort
  
  if(!emptyVectors){
	v_L0_phi.clear();
	v_L0_eta.clear();
	v_Ks_phi.clear();
	v_Ks_eta.clear();
  }  


} //end of analyzer


void Analyzer_V0_angular_correlation::endJob()
{
}

Analyzer_V0_angular_correlation::~Analyzer_V0_angular_correlation()
{
}


DEFINE_FWK_MODULE(Analyzer_V0_angular_correlation);
