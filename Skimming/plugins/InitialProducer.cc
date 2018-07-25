// -*- C++ -*-
//
// Package:    SexaQAnalysis/Skimming
// Class:      InitialProducer
// 
/**\class InitialProducer InitialProducer.cc SexaQAnalysis/Skimming/plugins/InitialProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  localusers user
//         Created:  Tue, 24 Jul 2018 10:26:06 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>

//
// class declaration
//

class InitialProducer : public edm::stream::EDProducer<> {
   public:
      explicit InitialProducer(const edm::ParameterSet&);
      ~InitialProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      //************remove 2 lines below
      edm::InputTag trackCollectionTag_;
      edm::EDGetTokenT<std::vector<reco::Track> > tracksCollectionToken_;
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::InputTag src_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
//InitialProducer::InitialProducer(const edm::ParameterSet& iConfig)
InitialProducer::InitialProducer(edm::ParameterSet const& pset)
//***************delte the line below
//:trackCollectionTag_(iConfig.getParameter<edm::InputTag>("TrackCollection"))
:trackCollectionTag_(pset.getParameter<edm::InputTag>("TrackCollection"))
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  //***********uncomment 2 lines below
//   src_  = iConfig.getParameter<edm::InputTag>( "src" );
//   produces<int>( "ntracks" ).setBranchAlias( "ntracks"); 
  //***********delete the 2 lines below
   tracksCollectionToken_ = consumes<std::vector<reco::Track> >(trackCollectionTag_);
   produces<std::vector<int>>("ntracks");
}


InitialProducer::~InitialProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
InitialProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco; 
   using namespace std;
   
   //get information on the tracks
 //*********uncomment the block below
 /*  edm::Handle<reco::TrackCollection> Tracks;
   iEvent.getByLabel( src_, Tracks );
   if(!Tracks.isValid()) {
      std::cout << "Missing collection during InitialProducer : ... skip entry !" << std::endl;
   }
   auto ntracks = std::make_unique<std::vector<int>>();
   ntracks->push_back((int)Tracks->size());
   
   iEvent.put(std::move(ntracks), "ntracks");
*/

  //********delete the block below
  edm::Handle<std::vector<reco::Track >> h_tracks;
  iEvent.getByToken(tracksCollectionToken_, h_tracks);
  if(!h_tracks.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << trackCollectionTag_ << " ... skip entry !" << std::endl;
  }

   auto ntracks = std::make_unique<std::vector<int>>();
   ntracks->push_back((int)h_tracks->size());
   
   iEvent.put(std::move(ntracks), "ntracks");

   

/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::unique_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(std::move(pOut));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
InitialProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
InitialProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
InitialProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
InitialProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
InitialProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
InitialProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
InitialProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(InitialProducer);
