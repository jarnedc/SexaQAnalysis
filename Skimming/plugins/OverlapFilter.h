#ifndef OverlapFilter_h
#define OverlapFilter_h
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <vector>

  
class OverlapFilter : public edm::EDFilter {

  public:

    explicit OverlapFilter(edm::ParameterSet const& cfg);
    virtual ~OverlapFilter() {}
    virtual bool filter(edm::Event & iEvent, edm::EventSetup const & iSetup);

  private:
  
    edm::InputTag lambdaCollectionTag_;
    edm::InputTag kshortCollectionTag_;
    edm::InputTag genCollectionTag_;
    edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > lambdaCollectionToken_;
    edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > kshortCollectionToken_;
    edm::EDGetTokenT<std::vector<reco::GenParticle> >              genCollectionToken_;
    bool isData_;
    unsigned int minNrLambda_,   minNrKshort_;
    unsigned int prescaleFalse_, nreject_;

};


#endif
