#include "../interface/Analyzer_GEN_Sexaq.h"
#include <typeinfo>


Analyzer_GEN_Sexaq::Analyzer_GEN_Sexaq(edm::ParameterSet const& pset):
  m_isData(pset.getUntrackedParameter<bool>("isData")),
  m_genParticlesTag_GEN(pset.getParameter<edm::InputTag>("genCollection_GEN")),
  m_genParticlesTag_SIM(pset.getParameter<edm::InputTag>("genCollection_SIM")),
  m_genParticlesTag_PLUSGEANT(pset.getParameter<edm::InputTag>("genCollection_PLUSGEANT")),
  m_genParticlesTag_HLT(pset.getParameter<edm::InputTag>("genCollection_HLT")),
  m_sCandsTag(pset.getParameter<edm::InputTag>("sexaqCandidates")),

  m_genParticlesToken_GEN(consumes<vector<reco::GenParticle> >(m_genParticlesTag_GEN)),
  m_genParticlesToken_SIM(consumes<vector<reco::GenParticle> >(m_genParticlesTag_SIM)),
  m_genParticlesToken_PLUSGEANT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_PLUSGEANT)),
  m_genParticlesToken_HLT(consumes<vector<reco::GenParticle> >(m_genParticlesTag_HLT)),
  m_sCandsToken(consumes<vector<reco::VertexCompositeCandidate> >(m_sCandsTag))

{

}


void Analyzer_GEN_Sexaq::beginJob() {
}


void Analyzer_GEN_Sexaq::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {


  cout << "------------------------------------------------------------------------------------------------------------------------------------------------------------------"<< endl; 
 
  //gen particles 
  edm::Handle<vector<reco::GenParticle>> h_genParticles_GEN;
  iEvent.getByToken(m_genParticlesToken_GEN, h_genParticles_GEN);
 
  //gen particles 
  edm::Handle<vector<reco::GenParticle>> h_genParticles_SIM;
  iEvent.getByToken(m_genParticlesToken_SIM, h_genParticles_SIM);
 
  //gen particles 
  edm::Handle<vector<reco::GenParticle>> h_genParticles_PLUSGEANT;
  iEvent.getByToken(m_genParticlesToken_PLUSGEANT, h_genParticles_PLUSGEANT);
 
  //gen particles 
  edm::Handle<vector<reco::GenParticle>> h_genParticles_HLT;
  iEvent.getByToken(m_genParticlesToken_HLT, h_genParticles_HLT);
 
  //lambdaKshortVertexFilter sexaquark candidates
  edm::Handle<vector<reco::VertexCompositeCandidate> > h_sCands;
  iEvent.getByToken(m_sCandsToken, h_sCands);

 
  //look at the GEN particles
  if(h_genParticles_GEN.isValid() )
  {
     cout << "-------------------------------------THE GEN PARTICLES:--------------------------------------------" << endl;
     for(unsigned int i = 0; i < h_genParticles_GEN->size(); ++i){
	if(h_genParticles_GEN->at(i).pdgId() == -1020000020 || h_genParticles_GEN->at(i).pdgId() == -3122 || h_genParticles_GEN->at(i).pdgId() == 310){ //Xi: 13324; antiS: -1020000020
		cout << h_genParticles_GEN->at(i).pdgId()  << " " << h_genParticles_GEN->at(i).px() << " " << h_genParticles_GEN->at(i).py()  << " " << h_genParticles_GEN->at(i).pz() << " " << h_genParticles_GEN->at(i).vx()	 << " " << h_genParticles_GEN->at(i).vy() << " " << h_genParticles_GEN->at(i).vz() << endl;
		//cout << "ndaughters: " << h_genParticles_GEN->at(i).numberOfDaughters() << endl;
        }
     }
  }

  //look at the SIM particles
  if(h_genParticles_SIM.isValid() )
  {
     cout << "-------------------------------------THE SIM PARTICLES:--------------------------------------------" << endl;
     for(unsigned int i = 0; i < h_genParticles_SIM->size(); ++i){
	if(h_genParticles_SIM->at(i).pdgId() == -1020000020 || h_genParticles_SIM->at(i).pdgId() == -3122 || h_genParticles_SIM->at(i).pdgId() == 310){ //Xi: 13324; antiS: -1020000020
		cout << h_genParticles_SIM->at(i).pdgId()  << " " << h_genParticles_SIM->at(i).px() << " " << h_genParticles_SIM->at(i).py()  << " " << h_genParticles_SIM->at(i).pz() << " " << h_genParticles_SIM->at(i).vx()	 << " " << h_genParticles_SIM->at(i).vy() << " " << h_genParticles_SIM->at(i).vz() << endl;
   		//cout << "ndaughters: " << h_genParticles_SIM->at(i).numberOfDaughters() << endl;
        }
     }
     for(unsigned int i = 0; i < h_genParticles_SIM->size(); ++i){
	if(h_genParticles_SIM->at(i).pdgId() == -3122){
		for(unsigned int j = 0; j < h_genParticles_SIM->size(); ++j){
			if(h_genParticles_SIM->at(j).pdgId() == 310){
				 Double_t energy_term = h_genParticles_SIM->at(i).energy() + h_genParticles_SIM->at(j).energy()-0.939565;
				 TVector3 momentum_term(h_genParticles_SIM->at(i).px()+h_genParticles_SIM->at(j).px(), h_genParticles_SIM->at(i).py()+h_genParticles_SIM->at(j).py(), h_genParticles_SIM->at(i).pz()+h_genParticles_SIM->at(j).pz());
				 //cout << "SIM antiS inv M: " << pow(pow(energy_term,2)-momentum_term.Mag2(),0.5) << endl;
			}
		}
	}
    }
  }

  //look at the PLUSGEANT particles
  if(h_genParticles_PLUSGEANT.isValid() )
  {
     cout << "-------------------------------------THE PLUSGEANT PARTICLES:--------------------------------------------" << endl;
     for(unsigned int i = 0; i < h_genParticles_PLUSGEANT->size(); ++i){
	if(h_genParticles_PLUSGEANT->at(i).pdgId() == -1020000020 || h_genParticles_PLUSGEANT->at(i).pdgId() == -3122 || h_genParticles_PLUSGEANT->at(i).pdgId() == 310){ //Xi: 13324; antiS: -1020000020
		cout << h_genParticles_PLUSGEANT->at(i).pdgId()  << " " << h_genParticles_PLUSGEANT->at(i).px() << " " << h_genParticles_PLUSGEANT->at(i).py()  << " " << h_genParticles_PLUSGEANT->at(i).pz() << " " << h_genParticles_PLUSGEANT->at(i).vx()	 << " " << h_genParticles_PLUSGEANT->at(i).vy() << " " << h_genParticles_PLUSGEANT->at(i).vz() << endl;
		//cout << "ndaughters: " << h_genParticles_PLUSGEANT->at(i).numberOfDaughters() << endl;
        }
	if(h_genParticles_PLUSGEANT->at(i).pdgId() == -1020000020 && h_genParticles_PLUSGEANT->at(i).numberOfDaughters() == 2){
		cout << "daughters of the PLUSGEANT antiS: " << endl;
		cout << h_genParticles_PLUSGEANT->at(i).daughter(0)->pdgId() << " " << h_genParticles_PLUSGEANT->at(i).daughter(0)->px() << " "  << h_genParticles_PLUSGEANT->at(i).daughter(0)->py() << " " << h_genParticles_PLUSGEANT->at(i).daughter(0)->pz() << " " << h_genParticles_PLUSGEANT->at(i).daughter(0)->vx()<< " " << h_genParticles_PLUSGEANT->at(i).daughter(0)->vy() << " " << h_genParticles_PLUSGEANT->at(i).daughter(0)->vz() << endl;
		cout << h_genParticles_PLUSGEANT->at(i).daughter(1)->pdgId() << " " << h_genParticles_PLUSGEANT->at(i).daughter(1)->px() << " "  << h_genParticles_PLUSGEANT->at(i).daughter(1)->py() << " " << h_genParticles_PLUSGEANT->at(i).daughter(1)->pz() << " " << h_genParticles_PLUSGEANT->at(i).daughter(1)->vx()<< " " << h_genParticles_PLUSGEANT->at(i).daughter(1)->vy() << " " << h_genParticles_PLUSGEANT->at(i).daughter(1)->vz() << endl;
		cout << "end daughters of the PLUSGEANT antiS: " << endl;
	}
	if(h_genParticles_PLUSGEANT->at(i).pdgId() == -1020000020 && h_genParticles_PLUSGEANT->at(i).numberOfDaughters() == 2){
		 Double_t energy_term = h_genParticles_PLUSGEANT->at(i).daughter(0)->energy() + h_genParticles_PLUSGEANT->at(i).daughter(1)->energy()-0.939565;
                 TVector3 momentum_term(h_genParticles_PLUSGEANT->at(i).daughter(0)->px()+h_genParticles_PLUSGEANT->at(i).daughter(1)->px(), h_genParticles_PLUSGEANT->at(i).daughter(0)->py()+h_genParticles_PLUSGEANT->at(i).daughter(1)->py(), h_genParticles_PLUSGEANT->at(i).daughter(0)->pz()+h_genParticles_PLUSGEANT->at(i).daughter(1)->pz());
                 //cout << "PLUSGEANT antiS inv M: " << pow(pow(energy_term,2)-momentum_term.Mag2(),0.5) << endl;
	}
     }
  }

  //look at the HLT particles
  if(h_genParticles_HLT.isValid() )
  {
     cout << "-------------------------------------THE HLT PARTICLES:--------------------------------------------" << endl;
     for(unsigned int i = 0; i < h_genParticles_HLT->size(); ++i){
	if(h_genParticles_HLT->at(i).pdgId() == -1020000020 || h_genParticles_HLT->at(i).pdgId() == -3122 || h_genParticles_HLT->at(i).pdgId() == 310){ //Xi: 13324; antiS: -1020000020
		cout << h_genParticles_HLT->at(i).pdgId()  << " " << h_genParticles_HLT->at(i).px() << " " << h_genParticles_HLT->at(i).py()  << " " << h_genParticles_HLT->at(i).pz() << " " << h_genParticles_HLT->at(i).vx()	 << " " << h_genParticles_HLT->at(i).vy() << " " << h_genParticles_HLT->at(i).vz() << endl;
		//cout << "ndaughters: " << h_genParticles_HLT->at(i).numberOfDaughters() << endl 
        }
     }
     for(unsigned int i = 0; i < h_genParticles_HLT->size(); ++i){
	if(h_genParticles_HLT->at(i).pdgId() == -3122){
		for(unsigned int j = 0; j < h_genParticles_HLT->size(); ++j){
			if(h_genParticles_HLT->at(j).pdgId() == 310){
				 Double_t energy_term = h_genParticles_HLT->at(i).energy() + h_genParticles_HLT->at(j).energy()-0.939565;
				 TVector3 momentum_term(h_genParticles_HLT->at(i).px()+h_genParticles_HLT->at(j).px(), h_genParticles_HLT->at(i).py()+h_genParticles_HLT->at(j).py(), h_genParticles_HLT->at(i).pz()+h_genParticles_HLT->at(j).pz());
			//	 cout << "HLT antiS inv M: " << pow(pow(energy_term,2)-momentum_term.Mag2(),0.5) << endl;
			}
		}
	}
     }
   
  }

  //look at the HLT particles
  if(h_sCands.isValid() )
  {
     cout << "-------------------------------------THE RECO PARTICLES:--------------------------------------------" << endl;
     for(unsigned int i = 0; i < h_sCands->size(); ++i){
		cout << h_sCands->at(i).pdgId()  << " " << h_sCands->at(i).px() << " " << h_sCands->at(i).py()  << " " << h_sCands->at(i).pz() << " " << h_sCands->at(i).vx()	 << " " << h_sCands->at(i).vy() << " " << h_sCands->at(i).vz() << endl;
     }
  }

} //end of analyzer

double Analyzer_GEN_Sexaq::openings_angle(reco::Candidate::Vector momentum1, reco::Candidate::Vector momentum2){
  double opening_angle = TMath::ACos((momentum1.Dot(momentum2))/(pow(momentum1.Mag2()*momentum2.Mag2(),0.5)));
  return opening_angle;
}

void Analyzer_GEN_Sexaq::endJob()
{
}

Analyzer_GEN_Sexaq::~Analyzer_GEN_Sexaq()
{
}


DEFINE_FWK_MODULE(Analyzer_GEN_Sexaq);
