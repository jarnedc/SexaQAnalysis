
#include "SexaQAnalysis/Skimming/plugins/OverlapFilter.h"


OverlapFilter::OverlapFilter(edm::ParameterSet const& pset):
  lambdaCollectionTag_(pset.getParameter<edm::InputTag>("lambdaCollection")),
  kshortCollectionTag_(pset.getParameter<edm::InputTag>("kshortCollection")),
  genCollectionTag_   (pset.getParameter<edm::InputTag>("genCollection")),
  isData_       (pset.getParameter<bool>  ("isData")),
  minNrLambda_  (pset.getParameter<unsigned int>("minNrLambda")),
  minNrKshort_  (pset.getParameter<unsigned int>("minNrKshort")),
  prescaleFalse_(pset.getParameter<unsigned int>("prescaleFalse"))
{
  lambdaCollectionToken_ = consumes<std::vector<reco::VertexCompositeCandidate> >(lambdaCollectionTag_);
  kshortCollectionToken_ = consumes<std::vector<reco::VertexCompositeCandidate> >(kshortCollectionTag_);
  genCollectionToken_    = consumes<std::vector<reco::GenParticle> >             (genCollectionTag_);
  nreject_ = 0;
  produces<reco::CandidatePtrVector>("kshort");
  produces<reco::CandidatePtrVector>("lambda");
}


bool OverlapFilter::filter(edm::Event & iEvent, edm::EventSetup const & iSetup)
{
  auto kshorts = std::make_unique<reco::CandidatePtrVector>();
  auto lambdas = std::make_unique<reco::CandidatePtrVector>();

  std::vector<int> goodLambdas;
  std::vector<int> goodKshorts;
  // select on reco lambdas and kaons
  edm::Handle<std::vector<reco::VertexCompositeCandidate> > h_lambda;
  iEvent.getByToken(lambdaCollectionToken_, h_lambda);
  edm::Handle<std::vector<reco::VertexCompositeCandidate> > h_kshort;
  iEvent.getByToken(kshortCollectionToken_ , h_kshort);

  if(!h_lambda.isValid()) {
    std::cout << "Missing collection during OverlapFilter : " << lambdaCollectionTag_ << " ... skip entry !" << std::endl;
    return false;
  }
  if(!h_kshort.isValid()) {
    std::cout << "Missing collection during OverlapFilter : " << kshortCollectionTag_ << " ... skip entry !" << std::endl;
    return false;
  }
  //std::cout << "OverlapFilter: Starting a new event" << std::endl;
  // select on reco lambdas and kaons
  if (isData_) {
    //std::cout << "OverlapFilter: Starting a new event" << std::endl;
    for (unsigned int l = 0; l < h_lambda->size(); ++l) {
        edm::Ptr<reco::VertexCompositeCandidate> lptr(h_lambda,l);
        lambdas->push_back(std::move(lptr));
	goodLambdas.push_back(l);
    }

  // select the kshorts passing kinematic cuts and non-overlapping with lambdas
    for (unsigned int k = 0; k < h_kshort->size(); ++k) {
        edm::Ptr<reco::VertexCompositeCandidate> kptr(h_kshort,k);
	// check for overlaps with the lambdas, and keep the lambda in case
        bool overlap = false;
        for (auto lptr : *lambdas) {
          for (unsigned int li = 0; li < lptr->numberOfDaughters() && !overlap; ++li) {
            for (unsigned int ki = 0; ki < kptr->numberOfDaughters() && !overlap; ++ki) {
	      if (lptr->daughter(li)->px() == kptr->daughter(ki)->px() &&
	          lptr->daughter(li)->py() == kptr->daughter(ki)->py() &&
	          lptr->daughter(li)->pz() == kptr->daughter(ki)->pz()) {
                overlap = true;
	      }
	    }
          }
          if (overlap){
		 //std::cout << "LambdaKshortFilter: OVERLAP FOUND" << std::endl;
		 break;
	  }	
        }
        if (!overlap){ 
		kshorts->push_back(std::move(kptr));
		goodKshorts.push_back(k);
	}
      //  kshorts->push_back(std::move(kptr));
    }

  // if not data, then select on gen particles
  } else {

    // read out genparticles
    edm::Handle<std::vector<reco::GenParticle> > h_genparts;
    iEvent.getByToken(genCollectionToken_, h_genparts);
    if(!h_genparts.isValid()) {
      std::cout << "Missing collection : " << genCollectionTag_ << " ... skip entry !" << std::endl;
      return false;
    }



  }

  // get the vector sizes before they disappear when putting in the event
  unsigned int nl = lambdas->size(), nk = kshorts->size();

  iEvent.put(std::move(lambdas),"lambda");
  iEvent.put(std::move(kshorts),"kshort");
/*  if(kshorts != NULL){
  for(int k = 0; k<(int)kshorts->size(); k++){
        std::cout << "OverlapFilter: kshort momenta " << k << "px,py,pz: " << (*kshorts)[k]->px() << "," << (*kshorts)[k]->py() << "," << (*kshorts)[k]->pz() << std::endl; 
  }
  }

  if(lambdas != NULL){
  for(int l = 0; l<(int)lambdas->size(); l++){
        std::cout << "OverlapFilter: lambda momenta " << l << "px,py,pz: " << (*lambdas)[l]->px() << "," << (*lambdas)[l]->py() << "," << (*lambdas)[l]->pz() << std::endl; 
  }
  }
*/
  //std::cout << "OverlapFilter: n lambdas put in events " << nl << std::endl;
  //std::cout << "OverlapFilter: n kshorts put in events " << nk << std::endl;

  // throw away events on data without sufficient lambdas or kshorts
  if (nl < minNrLambda_ || nk < minNrKshort_) {
    ++nreject_;
    return (prescaleFalse_ ? !(nreject_ % prescaleFalse_) : false);
  }
  std::cout << "OverlapFilter: n lambdas put in events " << nl << std::endl;
  std::cout << "OverlapFilter: n kshorts put in events " << nk << std::endl;
  for(unsigned int l = 0; l < goodLambdas.size(); l++){
        std::cout << "LambdaKshortFilter: momenta of the surviving lambda: " << h_lambda->at(goodLambdas[l]).px() << "," << h_lambda->at(goodLambdas[l]).py() << "," << h_lambda->at(goodLambdas[l]).pz() << std::endl;
  }
  for(unsigned int k = 0; k < goodKshorts.size(); k++){
        std::cout << "LambdaKshortFilter: momenta of the surviving kshort: " << h_kshort->at(goodKshorts[k]).px() << "," << h_kshort->at(goodKshorts[k]).py() << "," << h_kshort->at(goodKshorts[k]).pz() << std::endl;
  }

  // if we reach here there's a sufficient number of good lambdas and kshorts
 
  //for debugging: check that the lambda and kshort which you save are really different

/*  for(unsigned int l = 0; l < lambdas->size(); ++l){
 	std::cout << "put a lambda" << std::endl;
  }
  for(unsigned int k = 0; k < kshorts->size(); ++k){
	std::cout << "put a kshort" << std::endl;
  }

  for(unsigned int l = 0; l < lambdas->size(); ++l){
  	for(unsigned int k = 0; k < kshorts->size(); ++k){
	 if((*lambdas)[l]->px() == (*kshorts)[k]->px()){
  		std::cout << "------------------------------------------" << std::endl;
      		std::cout << "Lambda number " << l << " momenta: "<< (*lambdas)[l]->px()  << ", " << (*lambdas)[l]->py() << ", " <<  (*lambdas)[l]->pz() << std::endl;
      		std::cout << "Kshort number " << k << " momenta: "<<  (*kshorts)[k]->px() << ", " << (*kshorts)[k]->py() << ", " <<  (*kshorts)[k]->pz() << std::endl;
         }
  	}
  }
*/


   return true;

}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(OverlapFilter);
