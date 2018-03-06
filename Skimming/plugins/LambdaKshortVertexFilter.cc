#include "LambdaKshortVertexFilter.h"
using namespace reco;
using namespace edm;
using namespace std;

LambdaKshortVertexFilter::LambdaKshortVertexFilter(edm::ParameterSet const& pset):
  //collections
  lambdaCollectionTag_		(pset.getParameter<edm::InputTag>("lambdaCollection")),
  kshortCollectionTag_		(pset.getParameter<edm::InputTag>("kshortCollection")),
  //parameters
  maxchi2ndofVertexFit_  	(pset.getParameter<double>("maxchi2ndofVertexFit"))
{
  //collections
  lambdaCollectionToken_ = consumes<reco::CandidatePtrVector>(lambdaCollectionTag_);
  kshortCollectionToken_ = consumes<reco::CandidatePtrVector>(kshortCollectionTag_);
  //producer
  produces<std::vector<reco::VertexCompositePtrCandidate> >("");
}


//the real filter
bool LambdaKshortVertexFilter::filter(edm::Event & iEvent, edm::EventSetup const & iSetup)
{ 

  // initialize the transient track builder
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  auto sParticles = std::make_unique<std::vector<reco::VertexCompositePtrCandidate> >();

  // collections
  edm::Handle<reco::CandidatePtrVector> h_lambda;
  iEvent.getByToken(lambdaCollectionToken_, h_lambda);
  edm::Handle<reco::CandidatePtrVector> h_kshort;
  iEvent.getByToken(kshortCollectionToken_ , h_kshort);
  //check all the above collections and return false if any of them is invalid
  if (!allCollectionValid(h_lambda, h_kshort)) return false;

  std::vector<RefCountedKinematicParticle> lambdaKinFitted, kshortKinFitted;
  std::vector<unsigned int> lambdaIdx, kshortIdx;

  // loop over all the lambdas in an event
  for (unsigned int l = 0; l < h_lambda->size(); ++l) {
    //get the daughters from the Lambdas
    const Candidate * V0LambdasDaughter1 = (*h_lambda)[l]->daughter(0);
    const Candidate * V0LambdasDaughter2 = (*h_lambda)[l]->daughter(1);
    //get the tracks corresponding to these daughters
    const Track * TrackV0LambdasDaughterProton = 0;
    const Track * TrackV0LambdasDaughterPion = 0;
    //look at the mass assigned to the Candidate to check which track is the proton track and which one is the pion track
    if(fabs(V0LambdasDaughter1->mass() - proton_mass) < fabs(V0LambdasDaughter2->mass() - proton_mass)) {
      TrackV0LambdasDaughterProton = V0LambdasDaughter1->bestTrack();
      TrackV0LambdasDaughterPion = V0LambdasDaughter2->bestTrack();
    }
    else {
      TrackV0LambdasDaughterProton = V0LambdasDaughter2->bestTrack();
      TrackV0LambdasDaughterPion = V0LambdasDaughter1->bestTrack();
    }
    //get the ttracks
    TransientTrack TTrackV0LambdasDaughterProton = (*theB).build(TrackV0LambdasDaughterProton);
    TransientTrack TTrackV0LambdasDaughterPion = (*theB).build(TrackV0LambdasDaughterPion);
    //now do a kinfit on the two transient tracks from the lambda
    RefCountedKinematicTree LambdaTree = KinfitTwoTTracks(TTrackV0LambdasDaughterPion, TTrackV0LambdasDaughterProton, charged_pi_mass, charged_pi_mass_sigma, proton_mass, proton_mass_sigma, LambdaMass, LambdaMassSigma);
    //check if the LambdaTree is not a null pointer, isValid and is not empty. i.e.: the fit succeeded
    if(!checkRefCountedKinematicTree(LambdaTree)){cout << "Lambda tree not succesfully build" << endl; return false;}
    //get the Lambda particle from the tree
    lambdaKinFitted.push_back(getTopParticleFromTree(LambdaTree));
    lambdaIdx.push_back(l);
  }


  // loop over all kaons in the event
  for(unsigned int k = 0; k < h_kshort->size(); ++k){
    //get the daughters from the Kshorts
    const Candidate * V0KaonsDaughter1 = (*h_kshort)[k]->daughter(0);
    const Candidate * V0KaonsDaughter2 = (*h_kshort)[k]->daughter(1);
    //get the tracks corresponding to these daughters
    const Track * TrackV0KaonsDaughter1 = V0KaonsDaughter1->bestTrack();
    const Track * TrackV0KaonsDaughter2 = V0KaonsDaughter2->bestTrack();
    //get the ttracks
    TransientTrack TTrackV0KaonsDaughter1 = (*theB).build(TrackV0KaonsDaughter1);
    TransientTrack TTrackV0KaonsDaughter2 = (*theB).build(TrackV0KaonsDaughter2);
    //now do a kinfit to the two transient tracks
    RefCountedKinematicTree KshortTree = KinfitTwoTTracks(TTrackV0KaonsDaughter1, TTrackV0KaonsDaughter2, charged_pi_mass, charged_pi_mass_sigma, charged_pi_mass, charged_pi_mass_sigma, KshortMass, KshortMassSigma);
    if(!checkRefCountedKinematicTree(KshortTree)){cout << "Kshort tree not succesfully build" << endl; return false;}
    //get the Kshort particle from the tree
    kshortKinFitted.push_back(getTopParticleFromTree(KshortTree));
    kshortIdx.push_back(k);
  }

  for (unsigned int l = 0; l < lambdaKinFitted.size(); ++l) {
    for(unsigned int k = 0; k < kshortKinFitted.size(); ++k){
      vector<RefCountedKinematicParticle> daughtersS; 
      daughtersS.push_back(lambdaKinFitted.at(l));
      daughtersS.push_back(kshortKinFitted.at(k));
      //fit the S daughters to a common vertex
      KinematicParticleVertexFitter kpvFitter;
      RefCountedKinematicTree STree = kpvFitter.fit(daughtersS);
      if(!checkRefCountedKinematicTree(STree)){cout << "S tree not succesfully build" << endl; return false;};
      STree->movePointerToTheTop();
      RefCountedKinematicParticle Sparticle = STree->currentParticle();
      const reco::TransientTrack SparticleTrasientTrack = Sparticle->refittedTransientTrack();
      const reco::Track SparticleTrack = SparticleTrasientTrack.track();
      RefCountedKinematicVertex STreeVertex = STree->currentDecayVertex();

      //cut on some things now related to the Sparticle: the chi2 of the vertex fit, the the vertex location compared to the closest primary vertex, the fact if the reconstructed S momentum points to the primary vertex. 	
      if(STreeVertex->chiSquared()/STreeVertex->degreesOfFreedom() > maxchi2ndofVertexFit_)continue;

      //now that you passed all these cuts can start making the Sparticle candidate as a VertexCompositeCandidate
      const reco::Particle::LorentzVector SparticleP(Sparticle->currentState().globalMomentum().x(), Sparticle->currentState().globalMomentum().y(), Sparticle->currentState().globalMomentum().z(),sqrt(pow(Sparticle->currentState().globalMomentum().x()+Sparticle->currentState().globalMomentum().y()+Sparticle->currentState().globalMomentum().z(),2) + pow(Sparticle->currentState().kinematicParameters().mass(),2)));
      //reco::VertexCompositeCandidate* theSparticle = new reco::VertexCompositeCandidate(0,SparticleP, STreeVertex, STreeVertex.covariance(), STreeVertex->chiSquared(), STreeVertex->degreesOfFreedom());
      //Charge SCharge= 0;
      Point STreeVertexPoint(STreeVertex->position().x(),STreeVertex->position().y(),STreeVertex->position().z());
      //reco::CompositeCandidate* theSparticleCompositeCandidate = new reco::CompositeCandidate(0.,SparticleP, STreeVertexPoint, 0, 0, true, "Sparticle");
      //reco::VertexCompositeCandidate* theSparticleVertexCompositeCandidate = new reco::VertexCompositeCandidate(*theSparticleCompositeCandidate);
      reco::VertexCompositePtrCandidate theSparticleVertexCompositePtrCandidate(0, SparticleP, STreeVertexPoint);
      //adding daughters to the Sparticle
      theSparticleVertexCompositePtrCandidate.addDaughter((*h_lambda)[l]);
      theSparticleVertexCompositePtrCandidate.addDaughter((*h_kshort)[k]);
      // set the vertex info
      theSparticleVertexCompositePtrCandidate.setChi2AndNdof(STreeVertex->chiSquared(),STreeVertex->degreesOfFreedom());
      theSparticleVertexCompositePtrCandidate.setCovariance(STreeVertex->error().matrix());
      // adding Sparticles to the event
      sParticles->push_back(std::move(theSparticleVertexCompositePtrCandidate));
    }//end loop over kshort
  }//end loop over lambda
  
  int ns = sParticles->size();
  iEvent.put(std::move(sParticles)); 
  return (ns > 0);
}


//check the validness of all collections needed in the filter
bool LambdaKshortVertexFilter::allCollectionValid(edm::Handle<reco::CandidatePtrVector> h_lambda,edm::Handle<reco::CandidatePtrVector> h_kshort){
  if(!h_lambda.isValid()) {
      std::cout << "Missing collection : " << lambdaCollectionTag_   << " ... skip entry !" << std::endl;
      return false;
  }
  else if(!h_kshort.isValid()) {
      std::cout << "Missing collection : " << kshortCollectionTag_   << " ... skip entry !" << std::endl;
      return false;
  }
  else {
      return true;
  }
} 


//check if a RefCountedKinematicTree is not a nullpointer isValis and !isEmpty
bool LambdaKshortVertexFilter::checkRefCountedKinematicTree(RefCountedKinematicTree Tree){
  if(Tree){
      if(Tree->isValid() && !Tree->isEmpty()) { return true;}
      else{return false;}
  }
  else{return false;} 
}


//kinematic fit of two ttracks with the mass of each tracks, it's sigma and the mass of the parent
RefCountedKinematicTree LambdaKshortVertexFilter::KinfitTwoTTracks(TransientTrack ttrack1, TransientTrack ttrack2, ParticleMass trackMass1, float trackMassSigma1, ParticleMass trackMass2, float trackMassSigma2, ParticleMass combinedMass, float combinedMassSigma){ 
  //making particles
  vector<RefCountedKinematicParticle> daughters;
  KinematicParticleFactoryFromTransientTrack pFactory;
  //assumption 1: only needed when trackMass1 == trackMass2, this is the case for example for the Kshort where both daughters are charged pions and have the same mass
  daughters.push_back(pFactory.particle (ttrack1,trackMass1,chi,ndf,trackMassSigma1));
  daughters.push_back(pFactory.particle (ttrack2,trackMass2,chi,ndf,trackMassSigma2));
  //creating the vertex fitter
  KinematicParticleVertexFitter kpvFitter;
  //creating the particle fitter
  KinematicParticleFitter csFitter;
  // creating the mass constraint
  KinematicConstraint * ParentMassConstr = new MassKinematicConstraint(combinedMass,combinedMassSigma);
  //reconstructing a Kshort decay tree
  RefCountedKinematicTree ParentTreeSequential = kpvFitter.fit(daughters);
  //update the tree with a constrained fit:
  if(!ParentTreeSequential->isEmpty()){
    ParentTreeSequential = csFitter.fit(ParentMassConstr,ParentTreeSequential);
    if(!ParentTreeSequential->isEmpty()){
    }
    else{
      cout << "second fit in KinfitTwoTTracks failed" << endl;
      return NULL; 
    }
  }
  else{
   cout << "first fit in KinfitTwoTTracks failed" << endl;
   return NULL; //return a null pointer 
  }

  return ParentTreeSequential;
}


RefCountedKinematicParticle LambdaKshortVertexFilter::getTopParticleFromTree(RefCountedKinematicTree Tree){
  Tree->movePointerToTheTop();
  return  Tree->currentParticle();
}


RefCountedKinematicVertex LambdaKshortVertexFilter::returnVertexFromTree(const RefCountedKinematicTree& myTree) const
{
  if (!myTree->isValid()) {
    cout <<"Tree is invalid. Fit failed.\n";
    return 0;
  }
  //accessing the tree components, move pointer to top
  myTree->movePointerToTheTop();
  RefCountedKinematicVertex dec_vertex = myTree->currentDecayVertex();
  return dec_vertex;
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(LambdaKshortVertexFilter);

