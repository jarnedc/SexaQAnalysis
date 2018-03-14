#include "LambdaKshortVertexFilter.h"
using namespace reco;
using namespace edm;
using namespace std;

LambdaKshortVertexFilter::LambdaKshortVertexFilter(edm::ParameterSet const& pset):
  //collections
  lambdaCollectionTag_		(pset.getParameter<edm::InputTag>("lambdaCollection")),
  kshortCollectionTag_		(pset.getParameter<edm::InputTag>("kshortCollection")),
  //genCollectionTag_   (pset.getParameter<edm::InputTag>("genCollection")),
  //parameters
  maxchi2ndofVertexFit_  	(pset.getParameter<double>("maxchi2ndofVertexFit"))
{
  //collections
  lambdaCollectionToken_ = consumes<reco::CandidatePtrVector>(lambdaCollectionTag_);
  kshortCollectionToken_ = consumes<reco::CandidatePtrVector>(kshortCollectionTag_);
  //genCollectionToken_    = consumes<std::vector<reco::GenParticle> > (genCollectionTag_);
  //producer
  produces<std::vector<reco::Track> >("sParticlesTracks");
  produces<std::vector<reco::VertexCompositeCandidate> >("sParticles");
}


//the real filter
bool LambdaKshortVertexFilter::filter(edm::Event & iEvent, edm::EventSetup const & iSetup)
{ 
  // initialize the transient track builder
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  //these are for the producer
  auto sParticlesTracks = std::make_unique<std::vector<reco::Track>>();
  auto sParticles = std::make_unique<std::vector<reco::VertexCompositeCandidate> >();


  // collections
  edm::Handle<reco::CandidatePtrVector> h_lambda;
  iEvent.getByToken(lambdaCollectionToken_, h_lambda);
  edm::Handle<reco::CandidatePtrVector> h_kshort;
  iEvent.getByToken(kshortCollectionToken_ , h_kshort);
 
  //check all the above collections and return false if any of them is invalid
  if (!allCollectionValid(h_lambda, h_kshort)) return false;

  //to save the results from the kinfit
  std::vector<int> chargeProton;
  std::vector<RefCountedKinematicParticle> lambdaKinFitted, kshortKinFitted;
  std::vector<RefCountedKinematicVertex> lambdaKinFittedVertex, kshortKinFittedVertex;

  
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
    //save the charge of the proton: if pos then we know the lamda was a particle if neg we know it was an antilambda
    chargeProton.push_back(TrackV0LambdasDaughterProton->charge());
    //get the ttracks
    TransientTrack TTrackV0LambdasDaughterProton = (*theB).build(TrackV0LambdasDaughterProton);
    TransientTrack TTrackV0LambdasDaughterPion = (*theB).build(TrackV0LambdasDaughterPion);
    //now do a kinfit on the two transient tracks from the lambda
    RefCountedKinematicTree LambdaTree = KinfitTwoTTracks(TTrackV0LambdasDaughterPion, TTrackV0LambdasDaughterProton, charged_pi_mass, charged_pi_mass_sigma, proton_mass, proton_mass_sigma, LambdaMass, LambdaMassSigma);
    //check if the LambdaTree is not a null pointer, isValid and is not empty. i.e.: the fit succeeded
    if(!checkRefCountedKinematicTree(LambdaTree)){cout << "Lambda tree not succesfully build" << endl; return false;}
    //get the Lambda particle from the tree and save it to a vector
    LambdaTree->movePointerToTheTop();
    lambdaKinFitted.push_back(LambdaTree->currentParticle());
    lambdaKinFittedVertex.push_back(LambdaTree->currentDecayVertex());
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
    RefCountedKinematicTree  KshortTree = KinfitTwoTTracks(TTrackV0KaonsDaughter1, TTrackV0KaonsDaughter2, charged_pi_mass, charged_pi_mass_sigma, charged_pi_mass, charged_pi_mass_sigma, KshortMass, KshortMassSigma);
    if(!checkRefCountedKinematicTree(KshortTree)){cout << "Kshort tree not succesfully build" << endl; return false;}
    //get the Kshort particle from the tree and put in a vector
   KshortTree->movePointerToTheTop(); 
   kshortKinFitted.push_back(KshortTree->currentParticle());
   kshortKinFittedVertex.push_back(KshortTree->currentDecayVertex());

  }

  for (unsigned int l = 0; l < lambdaKinFitted.size(); ++l) {
    for(unsigned int k = 0; k < kshortKinFitted.size(); ++k){
      //put the supposed daughters of the S togethers
      vector<RefCountedKinematicParticle> daughtersS; 
      daughtersS.push_back(lambdaKinFitted.at(l));
      daughtersS.push_back(kshortKinFitted.at(k));
      //fit the S daughters to a common vertex
      KinematicParticleVertexFitter kpvFitter;
      RefCountedKinematicTree STree = kpvFitter.fit(daughtersS);
      //check succesfulness of the kpvFitter
      if(!checkRefCountedKinematicTree(STree)){cout << "S tree not succesfully build" << endl; return false;};
      //get the S particle
      STree->movePointerToTheTop();
      RefCountedKinematicParticle Sparticle = STree->currentParticle();
      //get the decay vertex of the S
      RefCountedKinematicVertex STreeVertex = STree->currentDecayVertex();

      //cut on some things now related to the Sparticle: the chi2 of the vertex fit, the vertex location compared to the closest primary vertex, the fact if the reconstructed S momentum points to the primary vertex. 	
      if(STreeVertex->chiSquared()/STreeVertex->degreesOfFreedom() > maxchi2ndofVertexFit_)continue;
      //get the track associated with this particle
      const reco::TransientTrack SparticleTrasientTrack = Sparticle->refittedTransientTrack();
      const reco::Track SparticleTrack = SparticleTrasientTrack.track();
      //get parameters from this track to put them back in later with a higher momentum...(no better way to do this?)
      const Vector & StrackArtificialMomentum =  SparticleTrack.momentum()*1000.;
      const reco::Track SparticleTrackArtificial = Track(SparticleTrack.chi2(),SparticleTrack.ndof(),SparticleTrack.referencePoint(),StrackArtificialMomentum, SparticleTrack.charge(), SparticleTrack.covariance(), SparticleTrack.algo(), reco::TrackBase::undefQuality); 
      //now that you passed cuts you can start making the Sparticle candidate as a VertexCompositeCandidate
      //momentum
      const reco::Particle::LorentzVector SparticleP(Sparticle->currentState().globalMomentum().x(), Sparticle->currentState().globalMomentum().y(), Sparticle->currentState().globalMomentum().z(),sqrt(pow(Sparticle->currentState().globalMomentum().x()+Sparticle->currentState().globalMomentum().y()+Sparticle->currentState().globalMomentum().z(),2) + pow(Sparticle->currentState().kinematicParameters().mass(),2)));
      //decay vertex
      Point STreeVertexPoint(STreeVertex->position().x(),STreeVertex->position().y(),STreeVertex->position().z()); 

      //will use the charge in the VertexCompositeCandidate constructor to indicate if in the decay there is an antiproton present if the antiproton
      //create the S as VertexCompositeCandidate
      reco::VertexCompositeCandidate theSparticleVertexCompositeCandidate(chargeProton[l], SparticleP, STreeVertexPoint);
     theSparticleVertexCompositeCandidate.setCovariance(STreeVertex->error().matrix());
     theSparticleVertexCompositeCandidate.setChi2AndNdof(STreeVertex->chiSquared(),STreeVertex->degreesOfFreedom());
     

      //making daughters to the Sparticle
      //Lambda      
      //momentum
      const reco::Particle::LorentzVector LambdaDaughterP(lambdaKinFitted.at(l)->currentState().globalMomentum().x(), lambdaKinFitted.at(l)->currentState().globalMomentum().y(), lambdaKinFitted.at(l)->currentState().globalMomentum().z(),sqrt(pow(lambdaKinFitted.at(l)->currentState().globalMomentum().x()+lambdaKinFitted.at(l)->currentState().globalMomentum().y()+lambdaKinFitted.at(l)->currentState().globalMomentum().z(),2) + pow(lambdaKinFitted.at(l)->currentState().kinematicParameters().mass(),2)));
      //vertex
      const Point LambdaDaughterPoint(lambdaKinFittedVertex.at(l)->position().x(),lambdaKinFittedVertex.at(l)->position().y(),lambdaKinFittedVertex.at(l)->position().z());
     //daughter
     LeafCandidate LambdaDaughter(0,  LambdaDaughterP, LambdaDaughterPoint);
      //Kshort
      //momentum
      const reco::Particle::LorentzVector KshortDaughterP(kshortKinFitted.at(k)->currentState().globalMomentum().x(), kshortKinFitted.at(k)->currentState().globalMomentum().y(), kshortKinFitted.at(k)->currentState().globalMomentum().z(),sqrt(pow(kshortKinFitted.at(k)->currentState().globalMomentum().x()+kshortKinFitted.at(k)->currentState().globalMomentum().y()+kshortKinFitted.at(k)->currentState().globalMomentum().z(),2) + pow(kshortKinFitted.at(k)->currentState().kinematicParameters().mass(),2)));
     //vertex
     const Point KshortDaughterPoint(kshortKinFittedVertex.at(k)->position().x(),kshortKinFittedVertex.at(k)->position().y(),kshortKinFittedVertex.at(k)->position().z()); 
     //daughter
     LeafCandidate KshortDaughter(0,  KshortDaughterP, KshortDaughterPoint);
     //add daughters to the S
     theSparticleVertexCompositeCandidate.addDaughter(LambdaDaughter);
     theSparticleVertexCompositeCandidate.addDaughter(KshortDaughter);
     

       //adding SparticlesTracks to the event
      sParticlesTracks->push_back(std::move(SparticleTrackArtificial));
       //adding Sparticles to the event
      sParticles->push_back(std::move(theSparticleVertexCompositeCandidate));
    }//end loop over kshort
  }//end loop over lambda




  
  int ns = sParticles->size();
  iEvent.put(std::move(sParticlesTracks),"sParticlesTracks");
  iEvent.put(std::move(sParticles),"sParticles"); 
  return (ns > 0);

}//end filter


//check the validness of all collections needed in the filter
bool LambdaKshortVertexFilter::allCollectionValid(edm::Handle<reco::CandidatePtrVector> h_lambda,edm::Handle<reco::CandidatePtrVector> h_kshort){
  if(!h_lambda.isValid()) {
      std::cout << "Missing collection during LambdaKshortVertexFilter: " << lambdaCollectionTag_   << " ... skip entry !" << std::endl;
      return false;
  }
  else if(!h_kshort.isValid()) {
      std::cout << "Missing collection during LambdaKshortVertexFilter: " << kshortCollectionTag_   << " ... skip entry !" << std::endl;
      return false;
  }
/*  else if(!h_beamspot.isValid()) {
      std::cout << "Missing collection : " << beamspotCollectionTag_   << " ... skip entry !" << std::endl;
      return false;
  }*/
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

