// -*- C++ -*-
//
// Package:    ttHGenAnalyzer
// Class:      ttHGenAnalyzer
//
// Original Author:  Zhengcheng Tao
//         Created:  Mon, 26 Jun 2017 18:20:01 GMT
//
//

// system include files
#include <memory>
#include <iostream>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "RecoTauTag/TauTagTools/interface/GeneratorTau.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace std;

class ttHGenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ttHGenAnalyzer(const edm::ParameterSet&);
      ~ttHGenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	// member function
	
	void printDaughters(const reco::Candidate& p, unsigned int& index,
						int mother_index, bool kinematics, bool recursive)
	{
		
		int ndaughters = p.numberOfDaughters();
		
		for (int j=0; j < ndaughters; ++j) {
			const reco::Candidate *daug = p.daughter(j);
			index++;
			cout << index << "\t" << daug->pdgId() << "\t" << daug->status()
				 << "\t" << mother_index << "\t" << daug->numberOfDaughters();
			if (kinematics) {
				cout << "\t" << daug->pt() << "\t" << daug->eta() << "\t"
					 << daug->phi() << "\t" << daug->energy()<< "\t"
					 << daug->mass() << endl;
			} else
				cout << endl;
			
			if (daug->status() != 1 and recursive)
				printDaughters(*daug, index, index, kinematics, true);
		}
	}

	void printDecayChain(const reco::GenParticle& p, bool kinematics=false)
	{

		int ndaughters = p.numberOfDaughters();
		unsigned int index = 0;
		
		cout << "index" << "\t" << "pdgId" << "\t" << "status" << "\t"
			 << "mother" << "\t" << "nDaughters" << "\t" << "isHardProcess"
			 << "\t" << "isLastCopy";
		if (kinematics) {
			cout << "\t" << "pt" << "\t" << "eta" << "\t" << "phi" << "\t"
				 << "energy" << "\t" << "mass" <<endl;
		} else
			cout << endl;

		cout << index << "\t" << p.pdgId() << "\t" << p.status() << "\t"
			 << "\t" << ndaughters << "\t" << p.isHardProcess() << "\t"
			 << p.isLastCopy();
		if (kinematics) {
			cout << "\t" << p.pt() << "\t" << p.eta() << "\t" << p.phi()
				 << "\t" << p.energy() << "\t" << p.mass() << endl;
		}
		else
			cout << endl;

		if (p.numberOfDaughters()<1) {
			cout << "The particle is stable." << endl;
			return;
		}
		else {
			const reco::Candidate * cand = &p;
 			printDaughters(*cand, index, 0, kinematics, true);
		}
	}
	
	const reco::GenParticle& getLastCopy(const reco::GenParticle& particle,
										 const vector<reco::GenParticle>& collection)
	{	
		if (particle.isLastCopy()) {
			return particle;
		}
		else {
			assert(particle.numberOfDaughters()>0);

			const reco::GenParticle* copy = 0;
			for (size_t i=0; i<particle.numberOfDaughters(); ++i) {
				const reco::GenParticle& daug =
					collection[(particle.daughterRef(i)).key()];
				
				if (daug.pdgId()==particle.pdgId()) {
					copy = &daug;
					break;
				}
			}

			assert(copy->pdgId()==particle.pdgId());
			return getLastCopy(*copy, collection);
		}
	}
	
	float Upsilon(const TLorentzVector& p4_leadtrack,
				  const TLorentzVector& p4_visible) 
	{
		// (2*E_ldgtrk - E_vis) / E_vis
		return (2*p4_leadtrack.Energy()/p4_visible.Energy() - 1);
	}

	float TauDecayEnergyAsymmetry(const TLorentzVector& p4_charged,
								  const TLorentzVector& p4_neutral)
	{
		// (E_charged - E_neutral) / (E_charged + E_neutral) of tau decay products
		return (p4_charged.Energy() - p4_neutral.Energy()) / (p4_charged.Energy() + p4_neutral.Energy());
	}

	float lambda(float a, float b, float c)
	{
		return a*a + b*b + c*c - 2*a*b - 2*b*c - 2*a*c;
	}
	
	float CosPsi(const TLorentzVector& p1, const TLorentzVector& p2,
				 const TLorentzVector& p3, float mass = 0.139)
	// lab frame, 3 prongs only
	{
		// angle between the direction of the normal to the plane defined by the three pions and the direction of flight of a1
		float s = (p1+p2+p3).M2();
		float s12 = (p1+p2).M2();
		float s23 = (p2+p3).M2();
		float s13 = (p1+p3).M2();

		TVector3 v1 = p1.Vect();
		TVector3 v2 = p2.Vect();
		TVector3 v3 = p3.Vect();

		float m2 = mass*mass;
		
		float numerator = v1.Dot(v2.Cross(v3)) / (v1+v2+v3).Mag();
		float denominator = TMath::Sqrt(-lambda(lambda(s,s12,m2),lambda(s,s13,m2),
												lambda(s,s23,m2)));

		return 8*s*numerator/denominator;
	}
	
private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;
	
	// ----------member data ---------------------------
	edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;
	//edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
	//edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;

	int boson_pdgid_;
	
	// Output file
	edm::Service<TFileService> fs_;

	// event variables to store
	int decayMode_plus_;
	int decayMode_minus_;
	// lab frame
	float evis_plus_lab_;  // visble energy / (boson mass / 2)
	float evis_minus_lab_;
	float upsilon_plus_lab_;  // (E_charged - E_neutral) / (E_charged + E_neutral)
	float upsilon_minus_lab_;
	float cosTauTheta_plus_lab_;  // angle between tau and visible decay products
	float cosTauTheta_minus_lab_;
	float x1ThreeProngs_plus_lab_; // E1/(E1+E2+E3) of 3 pions
	float x1ThreeProngs_minus_lab_;
	float x2ThreeProngs_plus_lab_;
	float x2ThreeProngs_minus_lab_;
	float x3ThreeProngs_plus_lab_;
	float x3ThreeProngs_minus_lab_;
	// Boson rest frame
	float evis_plus_bRF_;
	float evis_minus_bRF_;
	float upsilon_plus_bRF_;
	float upsilon_minus_bRF_;
	float cosTauTheta_plus_bRF_;
	float cosTauTheta_minus_bRF_;
	float x1ThreeProngs_plus_bRF_; // E1/(E1+E2+E3) of 3 pions
	float x1ThreeProngs_minus_bRF_;
	float x2ThreeProngs_plus_bRF_;
	float x2ThreeProngs_minus_bRF_;
	float x3ThreeProngs_plus_bRF_;
	float x3ThreeProngs_minus_bRF_;
	// visible tau pair rest frame
	float evis_plus_visRF_;
	float evis_minus_visRF_;
	float upsilon_plus_visRF_;
	float upsilon_minus_visRF_;
	//float cosTauTheta_plus_visRF_;
	//float cosTauTheta_minus_visRF_;
	//float x1ThreeProngs_plus_visRF_; // E1/(E1+E2+E3) of 3 pions
	//float x1ThreeProngs_minus_visRF_;
	//float x2ThreeProngs_plus_visRF_;
	//float x2ThreeProngs_minus_visRF_;
	//float x3ThreeProngs_plus_visRF_;
	//float x3ThreeProngs_minus_visRF_;

	// angle between visible daughter and tau direction of flight in tau plus frame
	float cos_plus_; 
	// angle between visible daughter and tau direction of flight in tau minus frame
	float cos_minus_;

	// 3-prong decay mode
	float cosPsi_plus_;
	float cosPsi_minus_;
	
	// Tree
	TTree* eventTree_;
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
ttHGenAnalyzer::ttHGenAnalyzer(const edm::ParameterSet& iConfig):
	boson_pdgid_ (iConfig.getParameter<int>("boson"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

   prunedGenToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("pruned"));
   //prunedGenToken_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"));
   //packedGenToken_ = consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"));

   // Setup tree
   eventTree_ = fs_->make<TTree>("eventTree", "Event tree");
   eventTree_->Branch("decayMode_plus", &decayMode_plus_);
   eventTree_->Branch("decayMode_minus", &decayMode_minus_);
   eventTree_->Branch("evis_plus_lab", &evis_plus_lab_);
   eventTree_->Branch("evis_minus_lab", &evis_minus_lab_);
   eventTree_->Branch("upsilon_plus_lab", &upsilon_plus_lab_);
   eventTree_->Branch("upsilon_minus_lab", &upsilon_minus_lab_);
   eventTree_->Branch("cosTauTheta_plus_lab", &cosTauTheta_plus_lab_);
   eventTree_->Branch("cosTauTheta_minus_lab", &cosTauTheta_minus_lab_);
   eventTree_->Branch("x1ThreeProngs_plus_lab", &x1ThreeProngs_plus_lab_);
   eventTree_->Branch("x1ThreeProngs_minus_lab", &x1ThreeProngs_minus_lab_);
   eventTree_->Branch("x2ThreeProngs_plus_lab", &x2ThreeProngs_plus_lab_);
   eventTree_->Branch("x2ThreeProngs_minus_lab", &x2ThreeProngs_minus_lab_);
   eventTree_->Branch("x3ThreeProngs_plus_lab", &x3ThreeProngs_plus_lab_);
   eventTree_->Branch("x3ThreeProngs_minus_lab", &x3ThreeProngs_minus_lab_);
   eventTree_->Branch("evis_plus_bRF", &evis_plus_bRF_);
   eventTree_->Branch("evis_minus_bRF", &evis_minus_bRF_);
   eventTree_->Branch("upsilon_plus_bRF", &upsilon_plus_bRF_);
   eventTree_->Branch("upsilon_minus_bRF", &upsilon_minus_bRF_);
   eventTree_->Branch("cosTauTheta_plus_bRF", &cosTauTheta_plus_bRF_);
   eventTree_->Branch("cosTauTheta_minus_bRF", &cosTauTheta_minus_bRF_);
   eventTree_->Branch("x1ThreeProngs_plus_bRF", &x1ThreeProngs_plus_bRF_);
   eventTree_->Branch("x1ThreeProngs_minus_bRF", &x1ThreeProngs_minus_bRF_);
   eventTree_->Branch("x2ThreeProngs_plus_bRF", &x2ThreeProngs_plus_bRF_);
   eventTree_->Branch("x2ThreeProngs_minus_bRF", &x2ThreeProngs_minus_bRF_);
   eventTree_->Branch("x3ThreeProngs_plus_bRF", &x3ThreeProngs_plus_bRF_);
   eventTree_->Branch("x3ThreeProngs_minus_bRF", &x3ThreeProngs_minus_bRF_);
   eventTree_->Branch("evis_plus_visRF", &evis_plus_visRF_);
   eventTree_->Branch("evis_minus_visRF", &evis_minus_visRF_);
   eventTree_->Branch("upsilon_plus_visRF", &upsilon_plus_visRF_);
   eventTree_->Branch("upsilon_minus_visRF", &upsilon_minus_visRF_);
   //eventTree_->Branch("cosTauTheta_plus_visRF", &cosTauTheta_plus_visRF_);
   //eventTree_->Branch("cosTauTheta_minus_visRF", &cosTauTheta_minus_visRF_);
   //eventTree_->Branch("x1ThreeProngs_plus_visRF", &x1ThreeProngs_plus_visRF_);
   //eventTree_->Branch("x1ThreeProngs_minus_visRF", &x1ThreeProngs_minus_visRF_);
   //eventTree_->Branch("x2ThreeProngs_plus_visRF", &x2ThreeProngs_plus_visRF_);
   //eventTree_->Branch("x2ThreeProngs_minus_visRF", &x2ThreeProngs_minus_visRF_);
   //eventTree_->Branch("x3ThreeProngs_plus_visRF", &x3ThreeProngs_plus_visRF_);
   //eventTree_->Branch("x3ThreeProngs_minus_visRF", &x3ThreeProngs_minus_visRF_);
   eventTree_->Branch("cos_plus", &cos_plus_);
   eventTree_->Branch("cos_minus", &cos_minus_);
   eventTree_->Branch("cosPsi_plus", &cosPsi_plus_);
   eventTree_->Branch("cosPsi_minus", &cosPsi_minus_);
}


ttHGenAnalyzer::~ttHGenAnalyzer(){}


//
// member functions
//

// ------------ method called for each event  ------------
void
ttHGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Handle<edm::View<reco::GenParticle> > pruned;
   Handle<reco::GenParticleCollection> pruned;
   iEvent.getByToken(prunedGenToken_,pruned);
   
   //Handle<edm::View<pat::PackedGenParticle> > packed;
   //iEvent.getByToken(packedGenToken_,packed);

   const reco::GenParticle* top = 0;
   const reco::GenParticle* antitop = 0;
   const reco::GenParticle* boson = 0;
   vector<GeneratorTau> taus_i;  // taus from higgs decay before FSR
   vector<GeneratorTau> taus;  // last copy of taus from higgs decay

   for (const auto & p : (*pruned)) {
	     
	   if (p.pdgId()==6) {  // top
		   if (not p.isLastCopy()) continue;  // p.isHardProcess() for first copy
		   top = &p;
	   }
	   
	   if (p.pdgId()==-6) { // antitop
		   if (not p.isLastCopy()) continue;  // p.isHardProcess() for first copy
		   antitop = &p;
	   }
	   
	   //if (p.pdgId()==23 or p.pdgId()==25) {  // higgs or Z
	   if (p.pdgId()==boson_pdgid_) {
		   
		   if (not p.isLastCopy()) continue;  // p.isHardProcess() for first copy
		   boson = &p;
		   
		   //printDecayChain(*boson);
		   
		   auto daughterRefs = p.daughterRefVector();	   
		   assert(daughterRefs.size()==2);
		   
		   // loop over daughters
		   for (auto & dauRef : daughterRefs) {
			   const reco::GenParticle& daug = (*pruned)[dauRef.key()];

			   // only looks for Higgs decaying to taus
			   if ( abs(daug.pdgId()) != 15 ) continue;

			   GeneratorTau tau_i(daug);
			   tau_i.init();
			   taus_i.push_back(tau_i);
			   
			   const reco::GenParticle& lasttau = getLastCopy(daug, *pruned);
			   GeneratorTau tau(lasttau);
			   tau.init();
			   taus.push_back(tau);
		   }
		   
	   }
	   	   
	   // print decay chain
	   //cout << "------------------------------------------------------" << endl;
	   //printDecayChain(p, false);
	   
   } // end of genParticles loop
   
   //assert(boson);
   if (!boson) return;
   
   if (taus.size()==0) return;
   assert(taus.size()==2);

   if (taus_i.size()==0) return;
   assert(taus_i.size()==2);
   
   // sort by pt
   //sort(taus.begin(), taus.end(), [] (GeneratorTau& t1, GeneratorTau& t2) {return t1.pt() > t2.pt();});

   TLorentzVector tauplusp4, visplusp4, chargedplusp4, neutralplusp4;
   TLorentzVector tauminusp4, visminusp4, chargedminusp4, neutralminusp4;
   TLorentzVector pion1tauplusp4, pion2tauplusp4, pion3tauplusp4;
   TLorentzVector pion1tauminusp4, pion2tauminusp4, pion3tauminusp4;
   GeneratorTau::tauDecayModeEnum decaymode_plus =
	   GeneratorTau::tauDecayModeEnum::kUndefined;
   GeneratorTau::tauDecayModeEnum decaymode_minus =
	   GeneratorTau::tauDecayModeEnum::kUndefined;
   
   for (const auto & tau : taus) {
	   assert(tau.isLastCopy());
	   //printDecayChain(tau);
	   
	   LorentzVector vis = tau.getVisibleFourVector();
	   auto chargedpions = tau.getGenChargedPions();
	   LorentzVector charged;
	   for (const auto & pion : chargedpions)
		   charged += pion->p4();
	   LorentzVector neutral = vis - charged;

	   vector<LorentzVector> pions_sorted;
	   if (tau.getDecayType()==GeneratorTau::tauDecayModeEnum::kThreeProng0pi0) {
		   assert(chargedpions.size()==3);
		   // assign pion1 and pion2 to the two indistinguishable pions
		   int i_distinguishable = -1;
		   for (unsigned int i=0; i<3; ++i) {
			   if (chargedpions[i]->charge() * tau.charge() < 0)
				   i_distinguishable = i;
			   else
				   pions_sorted.push_back(chargedpions[i]->p4());
			   
		   }
		   assert(i_distinguishable>=0 and i_distinguishable<3);
		   pions_sorted.push_back(chargedpions[i_distinguishable]->p4());
		   assert(pions_sorted.size()==3);
		   /*
		   if (neutral.energy()!=0) {
			   auto npis = tau.getGenNeutralPions();
			   cout << "number of neutral pions : " << npis.size() << endl;
			   auto stabledaugs = tau.getStableDecayProducts();
			   for (const auto& d: stabledaugs)
				   cout << d->pdgId() << endl;
		   }
		   */	   
	   }
	   
	   if (tau.charge()>0) { // tau+
		   decaymode_plus = tau.getDecayType();
		   tauplusp4.SetPtEtaPhiM(tau.pt(),tau.eta(),tau.phi(),tau.mass());
		   visplusp4.SetPtEtaPhiM(vis.pt(),vis.eta(),vis.phi(),vis.mass()); 
		   chargedplusp4.SetPtEtaPhiM(charged.pt(),charged.eta(),charged.phi(),charged.mass());
		   neutralplusp4.SetPtEtaPhiM(neutral.pt(),neutral.eta(),neutral.phi(),neutral.mass());
		   
		   if (pions_sorted.size() > 0) { // 3 prongs
			   pion1tauplusp4.SetPtEtaPhiM(pions_sorted[0].pt(),pions_sorted[0].eta(),pions_sorted[0].phi(),pions_sorted[0].mass());
			   pion2tauplusp4.SetPtEtaPhiM(pions_sorted[1].pt(),pions_sorted[1].eta(),pions_sorted[1].phi(),pions_sorted[1].mass());
			   pion3tauplusp4.SetPtEtaPhiM(pions_sorted[2].pt(),pions_sorted[2].eta(),pions_sorted[2].phi(),pions_sorted[2].mass());
		   }		   
	   }
	   else {  // tau-
		   decaymode_minus = tau.getDecayType();
		   tauminusp4.SetPtEtaPhiM(tau.pt(),tau.eta(),tau.phi(),tau.mass());
		   visminusp4.SetPtEtaPhiM(vis.pt(),vis.eta(),vis.phi(),vis.mass()); 
		   chargedminusp4.SetPtEtaPhiM(charged.pt(),charged.eta(),charged.phi(),charged.mass());
		   neutralminusp4.SetPtEtaPhiM(neutral.pt(),neutral.eta(),neutral.phi(),neutral.mass());
		   
		   if (pions_sorted.size() > 0) { // 3 prongs
			   pion1tauminusp4.SetPtEtaPhiM(pions_sorted[0].pt(),pions_sorted[0].eta(),pions_sorted[0].phi(),pions_sorted[0].mass());
			   pion2tauminusp4.SetPtEtaPhiM(pions_sorted[1].pt(),pions_sorted[1].eta(),pions_sorted[1].phi(),pions_sorted[1].mass());
			   pion3tauminusp4.SetPtEtaPhiM(pions_sorted[2].pt(),pions_sorted[2].eta(),pions_sorted[2].phi(),pions_sorted[2].mass());
		   }
	   }
   }

   assert(decaymode_plus!=GeneratorTau::tauDecayModeEnum::kUndefined and
		  decaymode_minus!=GeneratorTau::tauDecayModeEnum::kUndefined);

   decayMode_plus_ = decaymode_plus;
   decayMode_minus_ = decaymode_minus;

   // Lab frame
   TLorentzVector xp4;
   xp4.SetPtEtaPhiM(boson->pt(),boson->eta(),boson->phi(),boson->mass());
   
   evis_plus_lab_ = 2*visplusp4.Energy()/xp4.M();
   evis_minus_lab_ = 2*visminusp4.Energy()/xp4.M();
   upsilon_plus_lab_ = TauDecayEnergyAsymmetry(chargedplusp4, neutralplusp4);
   upsilon_minus_lab_ = TauDecayEnergyAsymmetry(chargedminusp4, neutralminusp4);
   cosTauTheta_plus_lab_ = TMath::Cos(visplusp4.Angle(tauplusp4.Vect()));
   cosTauTheta_minus_lab_ = TMath::Cos(visminusp4.Angle(tauminusp4.Vect()));
   if (decayMode_plus_ == GeneratorTau::tauDecayModeEnum::kThreeProng0pi0) {
	   float E1 = pion1tauplusp4.Energy();  float E2 = pion2tauplusp4.Energy();
	   float E3 = pion3tauplusp4.Energy();
	   float Ea = (pion1tauplusp4+pion2tauplusp4+pion3tauplusp4).Energy();
	   x1ThreeProngs_plus_lab_ = E1/Ea;
	   x2ThreeProngs_plus_lab_ = E2/Ea;
	   x3ThreeProngs_plus_lab_ = E3/Ea;
	   cosPsi_plus_ = CosPsi(pion1tauplusp4, pion2tauplusp4, pion3tauplusp4);
   }
   else {
	   x1ThreeProngs_plus_lab_ = -9;
	   x2ThreeProngs_plus_lab_ = -9;
	   x3ThreeProngs_plus_lab_ = -9;
	   cosPsi_plus_ = -9;
   }
   if (decayMode_minus_ == GeneratorTau::tauDecayModeEnum::kThreeProng0pi0) {
	   float E1 = pion1tauminusp4.Energy();  float E2 = pion2tauminusp4.Energy();
	   float E3 = pion3tauminusp4.Energy();
	   float Ea = (pion1tauminusp4+pion2tauminusp4+pion3tauminusp4).Energy();
	   x1ThreeProngs_minus_lab_ = E1/Ea;
	   x2ThreeProngs_minus_lab_ = E2/Ea;
	   x3ThreeProngs_minus_lab_ = E3/Ea;
	   cosPsi_minus_ = CosPsi(pion1tauminusp4, pion2tauminusp4, pion3tauminusp4);
   }
   else {
	   x1ThreeProngs_minus_lab_ = -9;
	   x2ThreeProngs_minus_lab_ = -9;
	   x3ThreeProngs_minus_lab_ = -9;
	   cosPsi_minus_ = -9;
   }
   
   // Boost to boson rest frame
   TLorentzVector taupx(tauplusp4), vispx(visplusp4);
   TLorentzVector chargedpx(chargedplusp4), neutralpx(neutralplusp4);
   TLorentzVector taumx(tauminusp4), vismx(visminusp4);
   TLorentzVector chargedmx(chargedminusp4), neutralmx(neutralminusp4);
   
   TVector3 xboost = xp4.BoostVector();
   xp4.Boost(-xboost);

   taupx.Boost(-xboost);
   vispx.Boost(-xboost);
   chargedpx.Boost(-xboost);
   neutralpx.Boost(-xboost);
   taumx.Boost(-xboost);
   vismx.Boost(-xboost);
   chargedmx.Boost(-xboost);
   neutralmx.Boost(-xboost);

   pion1tauplusp4.Boost(-xboost);
   pion2tauplusp4.Boost(-xboost);
   pion3tauplusp4.Boost(-xboost);
   pion1tauminusp4.Boost(-xboost);
   pion2tauminusp4.Boost(-xboost);
   pion3tauminusp4.Boost(-xboost);

   evis_plus_bRF_ = 2*vispx.Energy()/xp4.M();
   evis_minus_bRF_ = 2*vismx.Energy()/xp4.M();
   upsilon_plus_bRF_ = TauDecayEnergyAsymmetry(chargedpx, neutralpx);
   upsilon_minus_bRF_ = TauDecayEnergyAsymmetry(chargedmx, neutralmx);
   cosTauTheta_plus_bRF_ = TMath::Cos(vispx.Angle(taupx.Vect()));
   cosTauTheta_minus_bRF_ = TMath::Cos(vismx.Angle(taumx.Vect()));
   if (decayMode_plus_ == GeneratorTau::tauDecayModeEnum::kThreeProng0pi0) {
	   float E1 = pion1tauplusp4.Energy();  float E2 = pion2tauplusp4.Energy();
	   float E3 = pion3tauplusp4.Energy();
	   float Ea = (pion1tauplusp4+pion2tauplusp4+pion3tauplusp4).Energy();
	   x1ThreeProngs_plus_bRF_ = E1/Ea;
	   x2ThreeProngs_plus_bRF_ = E2/Ea;
	   x3ThreeProngs_plus_bRF_ = E3/Ea;
   }
   else {
	   x1ThreeProngs_plus_bRF_ = -9;
	   x2ThreeProngs_plus_bRF_ = -9;
	   x3ThreeProngs_plus_bRF_ = -9;
   }
   if (decayMode_minus_ == GeneratorTau::tauDecayModeEnum::kThreeProng0pi0) {
	   float E1 = pion1tauminusp4.Energy();  float E2 = pion2tauminusp4.Energy();
	   float E3 = pion3tauminusp4.Energy();
	   float Ea = (pion1tauminusp4+pion2tauminusp4+pion3tauminusp4).Energy();
	   x1ThreeProngs_minus_bRF_ = E1/Ea;
	   x2ThreeProngs_minus_bRF_ = E2/Ea;
	   x3ThreeProngs_minus_bRF_ = E3/Ea;
   }
   else {
	   x1ThreeProngs_minus_bRF_ = -9;
	   x2ThreeProngs_minus_bRF_ = -9;
	   x3ThreeProngs_minus_bRF_ = -9;
   }

   // Boost to tau+ rest frame
   TVector3 tpboost = taupx.BoostVector();
   vispx.Boost(-tpboost);
   cos_plus_ = TMath::Cos(vispx.Angle(tpboost));
      
   // Boost to tau- rest frame
   TVector3 tmboost = taumx.BoostVector();
   vismx.Boost(-tmboost);
   cos_minus_ = TMath::Cos(vismx.Angle(tmboost));

   // Boost to visible tau pair rest frame
   TVector3 visboost = (visplusp4+visminusp4).BoostVector();
   visplusp4.Boost(-visboost);
   visminusp4.Boost(-visboost);
   chargedplusp4.Boost(-visboost);
   neutralplusp4.Boost(-visboost);
   chargedminusp4.Boost(-visboost);
   neutralminusp4.Boost(-visboost);

   evis_plus_visRF_ = 2*visplusp4.Energy()/xp4.M();
   evis_minus_visRF_ = 2*visminusp4.Energy()/xp4.M();
   upsilon_plus_visRF_ = (chargedplusp4.Energy()-neutralplusp4.Energy())/(chargedplusp4.Energy()+neutralplusp4.Energy());
   upsilon_minus_visRF_ = (chargedminusp4.Energy()-neutralminusp4.Energy())/(chargedminusp4.Energy()+neutralminusp4.Energy());

   //h_taudecay_Evis_correlation_visRF[mode]->Fill(zmv, zpv);

   eventTree_->Fill();
   
   if (false) {
	   printDecayChain(*top);
	   printDecayChain(*antitop);
	   printDecayChain(*boson);
	   
	   for (const auto & tau_i : taus_i) 
		   printDecayChain(tau_i);
   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void 
ttHGenAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ttHGenAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ttHGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ttHGenAnalyzer);
