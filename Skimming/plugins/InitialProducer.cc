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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
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
      edm::InputTag lambdaCollectionTag_;
      edm::InputTag kshortCollectionTag_;
      edm::InputTag offlinePrimaryVerticesCollectionTag_;
      edm::InputTag ak4PFJetsCollectionTag_;
      edm::InputTag muonsCollectionTag_;
      edm::InputTag electronsCollectionTag_;
      edm::InputTag MHTCollectionTag_;
      edm::InputTag METCollectionTag_;
      edm::EDGetTokenT<std::vector<reco::Track> > tracksCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > lambdaCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > kshortCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > offlinePrimaryVerticesCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::TrackJet> > ak4PFJetsCollectionToken_;
      edm::EDGetTokenT<std::vector<reco::Muon> > muonsCollectionToken_;
      edm::EDGetTokenT<edm::ValueMap<edm::Ptr<reco::PFCandidate> >> electronsCollectionToken_;
      edm::EDGetTokenT<std::vector<l1extra::L1EtMissParticle> > MHTCollectionToken_;
      edm::EDGetTokenT<std::vector<l1extra::L1EtMissParticle>  > METCollectionToken_;
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
InitialProducer::InitialProducer(edm::ParameterSet const& pset):
trackCollectionTag_(pset.getParameter<edm::InputTag>("TrackCollection")),
lambdaCollectionTag_(pset.getParameter<edm::InputTag>("lambdaCollection")),
kshortCollectionTag_(pset.getParameter<edm::InputTag>("kshortCollection")),
offlinePrimaryVerticesCollectionTag_(pset.getParameter<edm::InputTag>("offlinePrimaryVerticesCollection")),
ak4PFJetsCollectionTag_(pset.getParameter<edm::InputTag>("ak4PFJetsCollection")),
muonsCollectionTag_(pset.getParameter<edm::InputTag>("muonsCollection")),
electronsCollectionTag_(pset.getParameter<edm::InputTag>("electronsCollection")),
MHTCollectionTag_(pset.getParameter<edm::InputTag>("MHTCollection")),
METCollectionTag_(pset.getParameter<edm::InputTag>("METCollection"))
{
   //register your products
   tracksCollectionToken_ = consumes<std::vector<reco::Track> >(trackCollectionTag_);
   lambdaCollectionToken_ = consumes<std::vector<reco::VertexCompositeCandidate> >(lambdaCollectionTag_);
   kshortCollectionToken_ = consumes<std::vector<reco::VertexCompositeCandidate> >(kshortCollectionTag_);
   offlinePrimaryVerticesCollectionToken_ = consumes<std::vector<reco::Vertex> >(offlinePrimaryVerticesCollectionTag_);
   ak4PFJetsCollectionToken_ = consumes<std::vector<reco::TrackJet> >(ak4PFJetsCollectionTag_);
   muonsCollectionToken_ = consumes<std::vector<reco::Muon> >(muonsCollectionTag_);
   electronsCollectionToken_ = consumes<edm::ValueMap<edm::Ptr<reco::PFCandidate> > >(electronsCollectionTag_);
   MHTCollectionToken_ = consumes<std::vector<l1extra::L1EtMissParticle> >(METCollectionTag_);
   METCollectionToken_ = consumes<std::vector<l1extra::L1EtMissParticle> >(MHTCollectionTag_);
   produces<std::vector<int>>("ntracks");
   produces<std::vector<int>>("nlambdas");
   produces<std::vector<int>>("nkshorts");
   produces<std::vector<int>>("nPVs");
   produces<std::vector<int>>("njets");
   produces<std::vector<reco::Particle::LorentzVector>>("TwoTopJets");
   produces<std::vector<int>>("nmuons");
   produces<std::vector<int>>("nelectrons");
   produces<std::vector<int>>("MHT");
   produces<std::vector<int>>("MET");
  
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
   
  //ntracks part
  edm::Handle<std::vector<reco::Track >> h_tracks;
  iEvent.getByToken(tracksCollectionToken_, h_tracks);
  if(!h_tracks.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << trackCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto ntracks = std::make_unique<std::vector<int>>();
  ntracks->push_back((int)h_tracks->size());
  iEvent.put(std::move(ntracks), "ntracks");

  //nlambdas part
  edm::Handle<std::vector<reco::VertexCompositeCandidate> > h_lambdas;
  iEvent.getByToken(lambdaCollectionToken_, h_lambdas);
  if(!h_lambdas.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << lambdaCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nlambdas = std::make_unique<std::vector<int>>();
  nlambdas->push_back((int)h_lambdas->size());
  iEvent.put(std::move(nlambdas), "nlambdas");

  //nkshorts part
  edm::Handle<std::vector<reco::VertexCompositeCandidate> > h_kshorts;
  iEvent.getByToken(kshortCollectionToken_, h_kshorts);
  if(!h_kshorts.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << kshortCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nkshorts = std::make_unique<std::vector<int>>();
  nkshorts->push_back((int)h_kshorts->size());
  iEvent.put(std::move(nkshorts), "nkshorts");

  //nPVs part
  edm::Handle<std::vector<reco::Vertex> > h_PVs;
  iEvent.getByToken(offlinePrimaryVerticesCollectionToken_, h_PVs);
  if(!h_PVs.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << offlinePrimaryVerticesCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nPVs = std::make_unique<std::vector<int>>();
  nPVs->push_back((int)h_PVs->size());
  iEvent.put(std::move(nPVs), "nPVs");

  //njets part
  edm::Handle<std::vector<reco::TrackJet> > h_jets;
  iEvent.getByToken(ak4PFJetsCollectionToken_, h_jets);
  if(!h_jets.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << ak4PFJetsCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto njets = std::make_unique<std::vector<int>>();
  njets->push_back((int)h_jets->size());
  iEvent.put(std::move(njets), "njets");
  //select the 2 jets with highest momentum
  auto TwoTopJets = std::make_unique<std::vector<reco::TrackJet>>();
  reco::TrackJet dummyTrackJet; 
  TwoTopJets->push_back(dummyTrackJet); 
  TwoTopJets->push_back(dummyTrackJet);
  for (unsigned int j = 0; j < h_jets->size(); ++j) {

	double thisJetMomentum = h_jets->at(j).p();
	if(thisJetMomentum > TwoTopJets->at(0).p()) TwoTopJets->at(0)  = h_jets->at(j);
	else if(thisJetMomentum > TwoTopJets->at(1).p()) TwoTopJets->at(1)  = h_jets->at(j);
 
  }
  auto highestMomentJetsLorentzVectors = std::make_unique<std::vector<reco::Particle::LorentzVector>>();
  highestMomentJetsLorentzVectors->push_back(TwoTopJets->at(0).p4());
  highestMomentJetsLorentzVectors->push_back(TwoTopJets->at(1).p4());
  iEvent.put(std::move(highestMomentJetsLorentzVectors), "TwoTopJets");

  //nmuons part
  edm::Handle<std::vector<reco::Muon> > h_muons;
  iEvent.getByToken(muonsCollectionToken_, h_muons);
  if(!h_muons.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << muonsCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nmuons = std::make_unique<std::vector<int>>();
  nmuons->push_back((int)h_muons->size());
  iEvent.put(std::move(nmuons), "nmuons");

  //nelectrons part
  edm::Handle<edm::ValueMap<edm::Ptr<reco::PFCandidate> > > h_electrons;
  iEvent.getByToken(electronsCollectionToken_, h_electrons);
  if(!h_electrons.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << electronsCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto nelectrons = std::make_unique<std::vector<int>>();
  nelectrons->push_back((int)h_electrons->size());
  iEvent.put(std::move(nelectrons), "nelectrons");

  //MHT part
  edm::Handle<vector<l1extra::L1EtMissParticle> > h_MHT;
  iEvent.getByToken(MHTCollectionToken_, h_MHT);
  if(!h_MHT.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << MHTCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto MHT = std::make_unique<std::vector<int>>();
  MHT->push_back((int)h_MHT->at(0).etMiss());
  iEvent.put(std::move(MHT), "MHT");

  //MET part
  edm::Handle<vector<l1extra::L1EtMissParticle> > h_MET;
  iEvent.getByToken(METCollectionToken_, h_MET);
  if(!h_MET.isValid()) {
      std::cout << "Missing collection during InitialProducer : " << METCollectionTag_ << " ... skip entry !" << std::endl;
  }
  auto MET = std::make_unique<std::vector<int>>();
  MET->push_back((int)h_MET->at(0).etMiss()); 
  iEvent.put(std::move(MET), "MET");

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
