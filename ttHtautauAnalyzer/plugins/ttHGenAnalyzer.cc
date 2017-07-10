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
	float evis_plus_lab_;
	float evis_minus_lab_;
	float upsilon_plus_lab_;
	float upsilon_minus_lab_;
	// Boson rest frame
	float evis_plus_bRF_;
	float evis_minus_bRF_;
	float upsilon_plus_bRF_;
	float upsilon_minus_bRF_;
	// visible tau pair rest frame
	float evis_plus_visRF_;
	float evis_minus_visRF_;
	float upsilon_plus_visRF_;
	float upsilon_minus_visRF_;

	float cos_plus_;
	float cos_minus_;
	
	// Tree
	TTree* eventTree_;

	/*
	// histograms
	// decay mode: [0]:pi+/pi-, [1]:pi+/rho-, [2]:pi+/l-, [3]:rho+/rho-, [4]:rho+/l-
	// not include a1 yet
	TH2F* h_taudecay_costheta_correlation[5];
	TH2F* h_taudecay_Evis_correlation[5];
	TH2F* h_taudecay_Evis_correlation_visRF[5];
	*/
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
   eventTree_->Branch("evis_plus_bRF", &evis_plus_bRF_);
   eventTree_->Branch("evis_minus_bRF", &evis_minus_bRF_);
   eventTree_->Branch("upsilon_plus_bRF", &upsilon_plus_bRF_);
   eventTree_->Branch("upsilon_minus_bRF", &upsilon_minus_bRF_);
   eventTree_->Branch("evis_plus_visRF", &evis_plus_visRF_);
   eventTree_->Branch("evis_minus_visRF", &evis_minus_visRF_);
   eventTree_->Branch("upsilon_plus_visRF", &upsilon_plus_visRF_);
   eventTree_->Branch("upsilon_minus_visRF", &upsilon_minus_visRF_);
   eventTree_->Branch("cos_plus", &cos_plus_);
   eventTree_->Branch("cos_minus", &cos_minus_);
   
   /*
   // histograms
   h_taudecay_costheta_correlation[0] = fs_->make<TH2F>("h_taudecay_costheta_pipi","1prong0pi0 + 1prong0pi0",20, -1, 1, 20, -1, 1);
   h_taudecay_costheta_correlation[1] = fs_->make<TH2F>("h_taudecay_costheta_pirho","1prong0pi0 + 1prong1pi0",20, -1, 1, 20, -1, 1);
   h_taudecay_costheta_correlation[2] = fs_->make<TH2F>("h_taudecay_costheta_pil","1prong0pi0 + e/mu",20, -1, 1, 20, -1, 1);
   h_taudecay_costheta_correlation[3] = fs_->make<TH2F>("h_taudecay_costheta_rhorho","1prong1pi0 + 1prong1pi0",20, -1, 1, 20, -1, 1);
   h_taudecay_costheta_correlation[4] = fs_->make<TH2F>("h_taudecay_costheta_rhol","1prong1pi0 + e/mu",20, -1, 1, 20, -1, 1);
   // set axis labels
   for (int i=0; i<5; ++i) {
	   h_taudecay_costheta_correlation[i]->GetXaxis()->SetTitle("cos#theta^{-}");
	   h_taudecay_costheta_correlation[i]->GetYaxis()->SetTitle("cos#theta^{+}");
   }
   
   h_taudecay_Evis_correlation[0] = fs_->make<TH2F>("h_taudecay_Evis_pipi","1prong0pi0 + 1prong0pi0 (Boson rest frame)",20, 0, 1, 20, 0, 1);
   h_taudecay_Evis_correlation[1] = fs_->make<TH2F>("h_taudecay_Evis_pirho","1prong0pi0 + 1prong1pi0 (Boson rest frame)",20, 0, 1, 20, 0, 1);
   h_taudecay_Evis_correlation[2] = fs_->make<TH2F>("h_taudecay_Evis_pil","1prong0pi0 + e/mu (Boson rest frame)",20, 0, 1, 20, 0, 1);
   h_taudecay_Evis_correlation[3] = fs_->make<TH2F>("h_taudecay_Evis_rhorho","1prong1pi0 + 1prong1pi0 (Boson rest frame)",20, 0, 1, 20, 0, 1);
   h_taudecay_Evis_correlation[4] = fs_->make<TH2F>("h_taudecay_Evis_rhol","1prong1pi0 + e/mu (Boson rest frame)",20, 0, 1, 20, 0, 1);
   // set axis labels
   for (int i=0; i<5; ++i) {
	   h_taudecay_Evis_correlation[i]->GetXaxis()->SetTitle("2E_{vis}^{-}/M_{X}");
	   h_taudecay_Evis_correlation[i]->GetYaxis()->SetTitle("2E_{vis}^{+}/M_{X}");
   }

   h_taudecay_Evis_correlation_visRF[0] = fs_->make<TH2F>("h_taudecay_Evis_pipi_visRF","1prong0pi0 + 1prong0pi0 (visible tau pair rest frame)",20, 0, 1, 20, 0, 1);
   h_taudecay_Evis_correlation_visRF[1] = fs_->make<TH2F>("h_taudecay_Evis_pirho_visRF","1prong0pi0 + 1prong1pi0 (visible tau pair rest frame)",20, 0, 1, 20, 0, 1);
   h_taudecay_Evis_correlation_visRF[2] = fs_->make<TH2F>("h_taudecay_Evis_pil_visRF","1prong0pi0 + e/mu (visible tau pair rest frame)",20, 0, 1, 20, 0, 1);
   h_taudecay_Evis_correlation_visRF[3] = fs_->make<TH2F>("h_taudecay_Evis_rhorho_visRF","1prong1pi0 + 1prong1pi0 (visible tau pair rest frame)",20, 0, 1, 20, 0, 1);
   h_taudecay_Evis_correlation_visRF[4] = fs_->make<TH2F>("h_taudecay_Evis_rhol_visRF","1prong1pi0 + e/mu (visible tau pair rest frame)",20, 0, 1, 20, 0, 1);
   // set axis labels
   for (int i=0; i<5; ++i) {
	   h_taudecay_Evis_correlation_visRF[i]->GetXaxis()->SetTitle("2E_{vis}^{-}/M_{X}");
	   h_taudecay_Evis_correlation_visRF[i]->GetYaxis()->SetTitle("2E_{vis}^{+}/M_{X}");
   } 
   */
   
}


ttHGenAnalyzer::~ttHGenAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


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
   GeneratorTau::tauDecayModeEnum decaymode_plus =
	   GeneratorTau::tauDecayModeEnum::kUndefined;
   GeneratorTau::tauDecayModeEnum decaymode_minus =
	   GeneratorTau::tauDecayModeEnum::kUndefined;
   
   for (const auto & tau : taus) {
	   assert(tau.isLastCopy());
	   //printDecayChain(tau);
	   
	   LorentzVector vis = tau.getVisibleFourVector();
	   vector<LorentzVector> chargedpions = tau.getChargedPions();
	   LorentzVector charged;
	   for (auto & pion : chargedpions)
		   charged += pion;
	   LorentzVector neutral = vis - charged;
	   
	   if (tau.charge()>0) { // tau+
		   tauplusp4.SetPtEtaPhiM(tau.pt(), tau.eta(), tau.phi(), tau.mass());
		   visplusp4.SetPtEtaPhiM(vis.pt(), vis.eta(), vis.phi(), vis.mass());
		   chargedplusp4.SetPtEtaPhiM(charged.pt(), charged.eta(), charged.phi(), charged.mass());
		   neutralplusp4.SetPtEtaPhiM(neutral.pt(), neutral.eta(), neutral.phi(), neutral.mass());
		   decaymode_plus = tau.getDecayType();
	   }
	   else {  // tau-
		   tauminusp4.SetPtEtaPhiM(tau.pt(), tau.eta(), tau.phi(), tau.mass());
		   visminusp4.SetPtEtaPhiM(vis.pt(), vis.eta(), vis.phi(), vis.mass());
		   chargedminusp4.SetPtEtaPhiM(charged.pt(), charged.eta(), charged.phi(), charged.mass());
		   neutralminusp4.SetPtEtaPhiM(neutral.pt(), neutral.eta(), neutral.phi(), neutral.mass());
		   decaymode_minus = tau.getDecayType();
	   }
   }

   assert(decaymode_plus!=GeneratorTau::tauDecayModeEnum::kUndefined and
		  decaymode_minus!=GeneratorTau::tauDecayModeEnum::kUndefined);
   /*
   int mode = -1;
   if (decaymode_plus==GeneratorTau::tauDecayModeEnum::kOneProng0pi0 and
	   decaymode_minus==GeneratorTau::tauDecayModeEnum::kOneProng0pi0) {
	   mode = 0;
   }
   else if (decaymode_plus==GeneratorTau::tauDecayModeEnum::kOneProng0pi0 and
			decaymode_minus==GeneratorTau::tauDecayModeEnum::kOneProng1pi0) {
	   mode = 1;
   }
   else if (decaymode_plus==GeneratorTau::tauDecayModeEnum::kOneProng0pi0 and
			(decaymode_minus==GeneratorTau::tauDecayModeEnum::kElectron or
			 decaymode_minus==GeneratorTau::tauDecayModeEnum::kMuon)) {
	   mode = 2;
   }
   else if (decaymode_plus==GeneratorTau::tauDecayModeEnum::kOneProng1pi0 and
			decaymode_minus==GeneratorTau::tauDecayModeEnum::kOneProng1pi0) {
	   mode = 3;
   }
   else if (decaymode_plus==GeneratorTau::tauDecayModeEnum::kOneProng1pi0 and
			(decaymode_minus==GeneratorTau::tauDecayModeEnum::kElectron or
			 decaymode_minus==GeneratorTau::tauDecayModeEnum::kMuon)) {
	   mode = 4;
   }
   else
   	   return;

   assert(mode > -1);
   */

   decayMode_plus_ = decaymode_plus;
   decayMode_minus_ = decaymode_minus;

   // Lab frame
   TLorentzVector xp4;
   xp4.SetPtEtaPhiM(boson->pt(), boson->eta(), boson->phi(), boson->mass());
   
   evis_plus_lab_ = 2*visplusp4.Energy()/xp4.M();
   evis_minus_lab_ = 2*visminusp4.Energy()/xp4.M();
   upsilon_plus_lab_ = (chargedplusp4.Energy()-neutralplusp4.Energy())/(chargedplusp4.Energy()+neutralplusp4.Energy());
   upsilon_minus_lab_ = (chargedminusp4.Energy()-neutralminusp4.Energy())/(chargedminusp4.Energy()+neutralminusp4.Energy());
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

   evis_plus_bRF_ = 2*vispx.Energy()/xp4.M();
   evis_minus_bRF_ = 2*vismx.Energy()/xp4.M();
   upsilon_plus_bRF_ = (chargedpx.Energy()-neutralpx.Energy())/(chargedpx.Energy()+neutralpx.Energy());
   upsilon_minus_bRF_ = (chargedmx.Energy()-neutralmx.Energy())/(chargedmx.Energy()+neutralmx.Energy());
   
   //h_taudecay_Evis_correlation[mode]->Fill(zm, zp);

   // Boost to tau+ rest frame
   TVector3 tpboost = taupx.BoostVector();
   vispx.Boost(-tpboost);
   cos_plus_ = TMath::Cos(vispx.Angle(tpboost));
      
   // Boost to tau- rest frame
   TVector3 tmboost = taumx.BoostVector();
   vismx.Boost(-tmboost);
   cos_minus_ = TMath::Cos(vismx.Angle(tmboost));

   //h_taudecay_costheta_correlation[mode]->Fill(cosm, cosp);

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
