#include "SexaQAnalysis/V0_angular_correlation_analysis/interface/Analyzer_V0_angular_correlation.h"
#include <typeinfo>


Analyzer_V0_angular_correlation::Analyzer_V0_angular_correlation(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  //m_rCandsTag(pset.getParameter<edm::InputTag>("resonCandidates")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),
  m_KshortsTag(pset.getParameter<edm::InputTag>("kshortCollection")),
  m_LambdasTag(pset.getParameter<edm::InputTag>("lambdaCollection")),
  m_LambdasLambdaKshortFilterTag(pset.getParameter<edm::InputTag>("lambdaCollectionlambdaKshortFilter")),
  m_KshortsLambdaKshortFilterTag(pset.getParameter<edm::InputTag>("kshortCollectionlambdaKshortFilter")),
  m_sCollectionMassFilterTag(pset.getParameter<edm::InputTag>("sCollectionMassFilter")),
  m_rCollectionMassFilterTag(pset.getParameter<edm::InputTag>("rCollectionMassFilter")),
  m_nPVsTag(pset.getParameter<edm::InputTag>("nPVsCollection")),
  m_nelectronsTag(pset.getParameter<edm::InputTag>("nelectronsCollection")),
  m_njetsTag(pset.getParameter<edm::InputTag>("njetsCollection")),
  m_nkshortsTag(pset.getParameter<edm::InputTag>("nkshortsCollection")),
  m_nlambdasTag(pset.getParameter<edm::InputTag>("nlambdasCollection")),
  m_nmuonsTag(pset.getParameter<edm::InputTag>("nmuonsCollection")),
  m_ntracksTag(pset.getParameter<edm::InputTag>("ntracksCollection")),
  //m_rCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_rCandsTag)),
  //m_sCandsToken(consumes<vector<reco::VertexCompositePtrCandidate> >(m_sCandsTag)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag)),
  m_KshortsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_KshortsTag)),
  m_LambdasToken(consumes<vector<reco::VertexCompositeCandidate> >(m_LambdasTag)),
  m_LambdasLambdaKshortFilterToken(consumes<edm::PtrVector<reco::Candidate >>(m_LambdasLambdaKshortFilterTag)),
  m_KshortsLambdaKshortFilterToken(consumes<edm::PtrVector<reco::Candidate >>(m_KshortsLambdaKshortFilterTag)),
  m_sCollectionMassFilterToken(consumes<vector<reco::VertexCompositePtrCandidate>>(m_sCollectionMassFilterTag)),
  m_rCollectionMassFilterToken(consumes<vector<reco::VertexCompositePtrCandidate>>(m_rCollectionMassFilterTag)),
  m_nPVsToken(consumes<vector<int>>(m_nPVsTag)),
  m_nelectronsToken(consumes<vector<int>>(m_nelectronsTag)),
  m_njetsToken(consumes<vector<int>>(m_njetsTag)),
  m_nkshortsToken(consumes<vector<int>>(m_nkshortsTag)),
  m_nlambdasToken(consumes<vector<int>>(m_nlambdasTag)),
  m_nmuonsToken(consumes<vector<int>>(m_nmuonsTag)),
  m_ntracksToken(consumes<vector<int>>(m_ntracksTag))


{

}


void Analyzer_V0_angular_correlation::beginJob() {
  
   //for scalar
    histos_th1f[b+"h_h_nPVs"]= m_fs->make<TH1F>(b+"h_h_nPVs",b+"h_h_nPVs; nPVs",1000,0,1000);
    histos_th1f[b+"h_h_nelectrons"]= m_fs->make<TH1F>(b+"h_h_nelectrons",b+"h_h_nelectrons; nelectrons",1000,0,1000);
    histos_th1f[b+"h_h_njets"]= m_fs->make<TH1F>(b+"h_h_njets",b+"h_h_njets; njets",1000,0,1000);
    histos_th1f[b+"h_h_nkshorts"]= m_fs->make<TH1F>(b+"h_h_nkshorts",b+"h_h_nkshorts; nkshorts",1000,0,1000);
    histos_th1f[b+"h_h_nlambdas"]= m_fs->make<TH1F>(b+"h_h_nlambdas",b+"h_h_nlambdas; nlambdas",1000,0,1000);
    histos_th1f[b+"h_h_nmuons"]= m_fs->make<TH1F>(b+"h_h_nmuons",b+"h_h_nmuons; nmuons",1000,0,1000);
    histos_th1f[b+"h_h_ntracks"]= m_fs->make<TH1F>(b+"h_h_ntracks",b+"h_h_ntracks; ntracks",1000,0,1000);

  
//  histos_th1f[b+"nPV"]      = m_fs->make<TH1F>(b+"nPV",     a+" Number of PV; #PVs",60,0.,60);
  
   //for angular correlation between Ks and L0 produced in the filters
   histos_th1f[b+"h_L0_Ks_delta_phi_Filters_Collection"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_phi_Filters_Collection",b+"h_L0_Ks_delta_phi_Filters_Collection; delta phi",1000,-4,4); 
   histos_th1f[b+"h_L0_Ks_delta_eta_Filters_Collection"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_eta_Filters_Collection",b+"h_L0_Ks_delta_eta_Filters_Collection; delta eta",4000,-10,10); 
   histos_th1f[b+"h_L0_Ks_delta_R_Filters_Collection"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_R_Filters_Collection",b+"h_L0_Ks_delta_R_Filters_Collection; delta R",3000,0,30); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection"]= m_fs->make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection",b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection; delta phi; delta_eta",100,-4,4, 100, -10, 10); 

   //for angular correlation between Ks and L0 
   histos_th1f[b+"h_L0_Ks_delta_phi"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_phi",b+"h_L0_Ks_delta_phi; delta phi",1000,-4,4); 
   histos_th1f[b+"h_L0_Ks_delta_eta"]= m_fs->make<TH1F>(b+"h_L0_Ks_delta_eta",b+"h_L0_Ks_delta_eta; delta eta",4000,-10,10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta"]= m_fs->make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta",b+"h_L0_Ks_delta_phi_delta_eta; delta phi; delta_eta",100,-4,4, 100, -10, 10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_no_cent"]= m_fs->make<TH2F>(b+"h_L0_Ks_delta_phi_delta_eta_no_cent",b+"h_L0_Ks_delta_phi_delta_eta_no_cent; delta phi; delta_eta",1000,-4,4, 4000, -10, 10); 
  
   //for angular correlation between Ks and L0 from different events
   histos_th1f[b+"h_L0_Ks_delta_phi_prev_and_current"]= m_fs->make<TH1F>(b+"h_L0_delta_phi_prev_and_current",b+"h_L0_delta_phi_prev_and_current; delta phi prev event and current event;",1000,-4,4); 
   histos_th1f[b+"h_L0_Ks_delta_eta_prev_and_current"]= m_fs->make<TH1F>(b+"h_L0_delta_eta_prev_and_current",b+"h_L0_delta_eta_prev_and_current; delta eta prev event and current event;",4000,-10,10); 
   histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]= m_fs->make<TH2F>(b+"h_L0_delta_phi_delta_eta_prev_and_current",b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current; delta phi prev event and current event; delta_eta previous event and current event",1000,-4,4, 4000, -10, 10); 
   
   //for angular coorelation between the daughters of the S candidates
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_phi"]= m_fs->make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_phi",b+"h_Sdaughters_L0_Ks_delta_phi; delta phi;",1000,-4,4); 
   histos_th1f[b+"h_Sdaughters_L0_Ks_delta_eta"]= m_fs->make<TH1F>(b+"h_Sdaughters_L0_Ks_delta_eta",b+"h_Sdaughters_L0_Ks_delta_eta; delta eta;",4000,-10,10); 
   histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta"]= m_fs->make<TH2F>(b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta",b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta; delta phi; delta eta",1000,-4,4, 4000, -10, 10); 
  
   //for the s candidates after the massfilter
    
   histos_th1f[b+"h_s_candidates_mass_after_massFilter"]= m_fs->make<TH1F>(b+"h_s_candidates_mass_after_massFilter",b+"h_s_candidates_mass_after_massFilter; r candidate mass (GeV);",40000,-200,200); 
    

   //for the r candidates after the massfilter
    
   histos_th1f[b+"h_r_candidates_mass_after_massFilter"]= m_fs->make<TH1F>(b+"h_r_candidates_mass_after_massFilter",b+"h_r_candidates_mass_after_massFilter; r candidate mass (GeV);",40000,-200,200); 
    
}


void Analyzer_V0_angular_correlation::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  //std::cout << "START LOOPING OVER THIS EVENT" << std::endl;


  // resonance candidates
  
  //edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_rCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
  //iEvent.getByToken(m_rCandsToken, h_rCands);

  // sexaquark candidates
  
  //edm::Handle<vector<reco::VertexCompositePtrCandidate> > h_sCands; //https://github.com/cms-sw/cmssw/blob/master/DataFormats/Candidate/interface/VertexCompositePtrCandidate.h
  //iEvent.getByToken(m_sCandsToken, h_sCands);
  
  //scalars
  edm::Handle <vector<int> > h_nPVs;
  iEvent.getByToken(m_nPVsToken, h_nPVs);

  edm::Handle <vector<int> > h_nelectrons;
  iEvent.getByToken(m_nelectronsToken, h_nelectrons);

  edm::Handle <vector<int> > h_njets;
  iEvent.getByToken(m_njetsToken, h_njets);
  
  edm::Handle <vector<int> > h_nkshorts;
  iEvent.getByToken(m_nkshortsToken, h_nkshorts);

  edm::Handle <vector<int> > h_nlambdas;
  iEvent.getByToken(m_nlambdasToken, h_nlambdas);

  edm::Handle <vector<int> > h_nmuons;
  iEvent.getByToken(m_nmuonsToken, h_nmuons);

  edm::Handle <vector<int> > h_ntracks;
  iEvent.getByToken(m_ntracksToken, h_ntracks);


  //lambdaKshortVertexFilter sexaquark candidates
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);
  
  //reco Kshorts V0
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_Kshorts;
  iEvent.getByToken(m_KshortsToken, h_Kshorts);
  
  //reco Lambdas V0
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_Lambdas;
  iEvent.getByToken(m_LambdasToken, h_Lambdas);

  //filtered Lambdas
  edm::Handle <edm::PtrVector<reco::Candidate > > h_Lambdas_LambdaKshortFilter;
  iEvent.getByToken(m_LambdasLambdaKshortFilterToken, h_Lambdas_LambdaKshortFilter);

  //filtered Kshorts
  edm::Handle <edm::PtrVector<reco::Candidate > > h_Kshorts_LambdaKshortFilter;
  iEvent.getByToken(m_KshortsLambdaKshortFilterToken, h_Kshorts_LambdaKshortFilter);

  //s candidates stored in the sMassFilter
  edm::Handle <vector<reco::VertexCompositePtrCandidate> > h_s_MassFilter;
  iEvent.getByToken(m_sCollectionMassFilterToken, h_s_MassFilter);

  //r candidates stored in the sMassFilter
  edm::Handle <vector<reco::VertexCompositePtrCandidate> > h_r_MassFilter;
  iEvent.getByToken(m_rCollectionMassFilterToken, h_r_MassFilter);



/*
  if(!h_rCands.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_rCandsTag << " ... skip entry !" << endl;
    return;
  }
*/
/*
  if(!h_Lambdas.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_LambdasTag << " ... skip entry !" << endl;
    return;
  }
  
  if(!h_Kshorts.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_KshortsTag << " ... skip entry !" << endl;
    return;
  }
  if(!h_Lambdas_LambdaKshortFilter.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_LambdasLambdaKshortFilterTag << " ... skip entry !" << endl;
    return;
  }
  if(!h_Kshorts_LambdaKshortFilter.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_KshortsLambdaKshortFilterTag << " ... skip entry !" << endl;
    return;
  }

  if(!h_s_MassFilter.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_sCollectionMassFilterTag << " ... skip entry !" << endl;
    return;
  }

  if(!h_r_MassFilter.isValid()) {
    if(verbose>0) cout << "Missing collection : " << m_rCollectionMassFilterTag << " ... skip entry !" << endl;
    return;
  }
 */ 
  
  //do stuff with the scalars produced in the filters
  std::cout << "the scalars" << std::endl;

 if(h_nPVs.isValid() && h_nPVs->size() > 0)  		 	histos_th1f[b+"h_h_nPVs"]->Fill(h_nPVs->at(0));
 if(h_nelectrons.isValid() && h_nelectrons->size() > 0 )	histos_th1f[b+"h_h_nelectrons"]->Fill(h_nelectrons->at(0));
 if(h_njets.isValid() && h_njets->size() > 0)			histos_th1f[b+"h_h_njets"]->Fill(h_njets->at(0));
 if(h_nkshorts.isValid() && h_nkshorts->size() > 0)		histos_th1f[b+"h_h_nkshorts"]->Fill(h_nkshorts->at(0));
 if(h_nlambdas.isValid() && h_nlambdas->size() > 0)		histos_th1f[b+"h_h_nlambdas"]->Fill(h_nlambdas->at(0));
 if(h_nmuons.isValid() && h_nmuons->size() > 0)			histos_th1f[b+"h_h_nmuons"]->Fill(h_nmuons->at(0));
// if(h_ntracks.isValid() && h_ntracks->size() > 0)		histos_th1f[b+"h_h_tracks"]->Fill(h_ntracks->at(0));


 
  std::cout << "LambdaKshortfilter" << std::endl;
  //do stuff with the collections produced in the filters
  if(h_Lambdas_LambdaKshortFilter.isValid() && h_Kshorts_LambdaKshortFilter.isValid()){
  std::cout << "------------------------------------------" << std::endl;
  for(unsigned int l = 0; l < h_Lambdas_LambdaKshortFilter->size(); ++l){
	for(unsigned int k = 0; k < h_Kshorts_LambdaKshortFilter->size(); ++k){
		double phi1 = (*h_Kshorts_LambdaKshortFilter)[k]->phi();
                double phi2 = (*h_Lambdas_LambdaKshortFilter)[l]->phi();
                double delta_phi = reco::deltaPhi(phi1, phi2);

                double eta1 = (*h_Kshorts_LambdaKshortFilter)[k]->eta();
                double eta2 = (*h_Lambdas_LambdaKshortFilter)[l]->eta();
                double delta_eta = eta1-eta2;

		double delta_R = sqrt(delta_phi*delta_phi+delta_eta*delta_eta);

                histos_th1f[b+"h_L0_Ks_delta_phi_Filters_Collection"]->Fill(delta_phi);
                histos_th1f[b+"h_L0_Ks_delta_eta_Filters_Collection"]->Fill(delta_eta);
                histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_Filters_Collection"]->Fill(delta_phi,delta_eta);	
                histos_th1f[b+"h_L0_Ks_delta_R_Filters_Collection"]->Fill(delta_R);	
		
   	}
  }
  }
/*
  for(unsigned int l = 0; l < h_Lambdas_LambdaKshortFilter->size(); ++l){
        for(unsigned int k = 0; k < h_Kshorts_LambdaKshortFilter->size(); ++k){
	
     }
   }
*/

   //std::cout << "looping over the lambdas and kshorts together" << std::endl;
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STANDARD PLOT: ONLY THE CURRENT EVENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  std::cout << "the lambdas and kshorts" << std::endl;
 if(h_Lambdas.isValid() && h_Kshorts.isValid()){
 for (unsigned int l = 0; l < h_Lambdas->size(); ++l) { //loop over reco V0 Lambdas
	for (unsigned int k = 0; k < h_Kshorts->size(); ++k) { //loop over reco Kshorts	
		//to check the sharp peak at 0 check if the kshort and lambdas are exactly the same sometimes
		//std::cout << "---------------------" << std::endl;
		//std::cout << setprecision(15) << "kshort momenta: " << h_Kshorts->at(k).px() << " " << h_Kshorts->at(k).py() << " " << h_Kshorts->at(k).pz() << std::endl;
		//std::cout << setprecision(15) << "lambda momenta: " << h_Lambdas->at(l).px() << " " << h_Lambdas->at(l).py() << " " << h_Lambdas->at(l).pz() << std::endl;
		//for debugging: check that the lambda and kshort which you save are really different
/*		 if(h_Lambdas->at(l).px() == h_Kshorts->at(k).px()){
			std::cout << "Lambda number " << l << " momenta: "<<  h_Lambdas->at(l).px()  << ", " << h_Lambdas->at(l).py() << ", " <<  h_Lambdas->at(l).pz() << std::endl;
			std::cout << "Kshort number " << k << " momenta: "<<  h_Kshorts->at(k).px()  << ", " << h_Kshorts->at(k).py() << ", " <<  h_Kshorts->at(k).pz() << std::endl;
		 }
*/
  		double phi1 = h_Kshorts->at(k).phi();
		double phi2 = h_Lambdas->at(l).phi();
		double delta_phi = reco::deltaPhi(phi1, phi2);
				
		double eta1 = h_Kshorts->at(k).eta();
		double eta2 = h_Lambdas->at(l).eta();
		double delta_eta = eta1-eta2;

		histos_th1f[b+"h_L0_Ks_delta_phi"]->Fill(delta_phi);
		histos_th1f[b+"h_L0_Ks_delta_eta"]->Fill(delta_eta);
		histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta"]->Fill(delta_phi,delta_eta);
		if( abs(delta_phi) > 0.1 || abs(delta_eta) > 0.1 )histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_no_cent"]->Fill(delta_phi,delta_eta);
	}
 }
 }


  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COMPARE THE PREVIOUS EVENT WITH THE CURRENT EVENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  bool emptyVectors = false;
  //std::cout << "---------------------" << std::endl;
  //std::cout << "going to check if the vectors are empty" << std::endl;
  if(v_L0_phi.empty() && v_L0_eta.empty() && v_Ks_phi.empty() && v_Ks_eta.empty()) emptyVectors = true;//this means they have been cleared in the previous event or it is the first event

 //if(emptyVectors) std::cout << "vectors empty" << std::endl;  
 //else std::cout << "vectors not empty" << std::endl;  
  
   //std::cout << "getting the size of the lambda and kshort collections" << std::endl;
  if(h_Lambdas.isValid()){
  for (unsigned int i = 0; i < h_Lambdas->size(); ++i) { //loop over reco V0 Lambdas
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
			histos_th1f[b+"h_L0_Ks_delta_phi_prev_and_current"]->Fill(delta_phi);
			histos_th1f[b+"h_L0_Ks_delta_eta_prev_and_current"]->Fill(delta_eta);
			histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]->Fill(delta_phi,delta_eta);
		}
	} 
  }//end loop over lambda
  }
  //std::cout << "loop over Lambdas done" << std::endl; 
  
  if(h_Kshorts.isValid()){
  for (unsigned int i = 0; i < h_Kshorts->size(); ++i) { //loop over reco Kshorts
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

			histos_th1f[b+"h_L0_Ks_delta_phi_prev_and_current"]->Fill(delta_phi);
			histos_th1f[b+"h_L0_Ks_delta_eta_prev_and_current"]->Fill(delta_eta);
			histos_th2f[b+"h_L0_Ks_delta_phi_delta_eta_prev_and_current"]->Fill(delta_phi,delta_eta);
		}
	}
  }// end loop over Kshort
 }

 // std::cout << "loop over Kshorts done" << std::endl; 

  if(!emptyVectors){
	v_L0_phi.clear();
	v_L0_eta.clear();
	v_Ks_phi.clear();
	v_Ks_eta.clear();
  } 
 // std::cout << "vectors cleared" << std::endl; 


 //!!!!!!!!!!!!!!!!!!!!!!!!!!!make the same distribution but for the S candidates (if they are there, they are only there if they ran through the 2nd filter)!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  std::cout << "the scands" << std::endl;
  if(h_sCands.isValid()) {
	unsigned int n_sCands = h_sCands->size();
	for (unsigned int i = 0; i < n_sCands; ++i) { //loop over reco Kshorts	
			double phi1 = h_sCands->at(i).daughter(0)->phi();
			double phi2 = h_sCands->at(i).daughter(1)->phi();
			double delta_phi = reco::deltaPhi(phi1, phi2);
			
			double eta1 = h_sCands->at(i).daughter(0)->eta();
			double eta2 = h_sCands->at(i).daughter(1)->eta();
			double delta_eta = eta1-eta2;
   			
			histos_th1f[b+"h_Sdaughters_L0_Ks_delta_phi"]->Fill(delta_phi);
                        histos_th1f[b+"h_Sdaughters_L0_Ks_delta_eta"]->Fill(delta_eta);
                        histos_th2f[b+"h_Sdaughters_L0_Ks_delta_phi_delta_eta"]->Fill(delta_phi,delta_eta);
		}

  }//end of sCands present

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!make some distributions for the S and r candidates
  std::cout << "the rMassFilter" << std::endl;
  if((h_r_MassFilter).isValid()){
  for (unsigned int r = 0; r < (*h_r_MassFilter).size(); ++r) {
	
	histos_th1f[b+"h_r_candidates_mass_after_massFilter"]->Fill((*h_r_MassFilter)[r].mass());
  }
  }
 
  std::cout << "the sMassFilter" << std::endl;
  if((h_s_MassFilter).isValid()){
  for (unsigned int s = 0; s < (*h_s_MassFilter).size(); ++s) {
	
	histos_th1f[b+"h_s_candidates_mass_after_massFilter"]->Fill((*h_s_MassFilter)[s].mass());
  }
  }


} //end of analyzer


void Analyzer_V0_angular_correlation::endJob()
{
}

Analyzer_V0_angular_correlation::~Analyzer_V0_angular_correlation()
{
}


DEFINE_FWK_MODULE(Analyzer_V0_angular_correlation);
