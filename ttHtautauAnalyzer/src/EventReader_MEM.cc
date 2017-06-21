
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/EventReader_MEM.h"

using namespace std;

// constructor
EventReader_MEM::EventReader_MEM(const TString& filename,const int64_t startEvent,
								 int64_t maxNbrEventsToRead, bool loose,
								 const TString& treename)
{	
	infile_ = new TFile(filename, "read");
	intree_ = (TTree*)infile_->Get(treename);

	evNtuple_.set_branch_address(intree_);
	
	startevent_ = startEvent;
	nevents_ = maxNbrEventsToRead;

	assert(startevent_>=0 and nevents_>=0);
	
	if (startevent_+nevents_ >= intree_->GetEntries()) {
		startevent_ = intree_->GetEntries()-1;
		nevents_ = 0;
	}	
	else if (startevent_+nevents_ > intree_->GetEntries()) {
		nevents_ = intree_->GetEntries() - startevent_;
	}

	looseSel_ = loose;
}
// destructor
EventReader_MEM::~EventReader_MEM()
{
	infile_->Close();
}


void EventReader_MEM::fillEventList(IntegrationMsg_t* integration, int64_t nEvents)
{
	int64_t k = 0; // Number of event stored
	
	for (int64_t i=startevent_; i < startevent_+nevents_; ++i) {
		intree_->GetEntry(i);

		// apply additional selection here if needed
		
		// build four momentum
		auto leptons = evNtuple_.buildLeptons(looseSel_);
		auto taus = evNtuple_.buildFourVectorTaus(looseSel_);
		auto jets = evNtuple_.buildFourVectorJets(looseSel_);
		auto jets_btag = evNtuple_.buildFourVectorBtagJets(looseSel_);
		auto met = evNtuple_.buildFourVectorMET();
		
		std::vector<TLorentzVector> jets_untag;
		set_difference(jets.begin(),jets.end(),jets_btag.begin(),jets_btag.end(),
					   inserter(jets_untag, jets_untag.end()),
					   [] (TLorentzVector l1, TLorentzVector l2)
					   {return l1.Pt() > l2.Pt();});

		assert(leptons.size()>1);
		assert(taus.size()>0);
		assert(jets_btag.size()>1);

		// Fill event list
		// leptons
		integration[k].evLep1_4P_[0] = leptons[0].p4().Px();
		integration[k].evLep1_4P_[1] = leptons[0].p4().Py(); 
		integration[k].evLep1_4P_[2] = leptons[0].p4().Pz();
		integration[k].evLep1_4P_[3] = leptons[0].p4().E();
		integration[k].lepton1_Type_ = leptons[0].pdgId();

		integration[k].evLep2_4P_[0] = leptons[1].p4().Px();
		integration[k].evLep2_4P_[1] = leptons[1].p4().Py(); 
		integration[k].evLep2_4P_[2] = leptons[1].p4().Pz();
		integration[k].evLep2_4P_[3] = leptons[1].p4().E();
		integration[k].lepton2_Type_ = leptons[1].pdgId();
		// tau
		integration[k].evHadSys_Tau_4P_[0] = taus[0].Px();
		integration[k].evHadSys_Tau_4P_[1] = taus[0].Py();
		integration[k].evHadSys_Tau_4P_[2] = taus[0].Pz();
		integration[k].evHadSys_Tau_4P_[3] = taus[0].E();
		integration[k].HadtauDecayMode_ =  evNtuple_.tau_decayMode->at(0);
		// btagged jets
		integration[k].evBJet1_4P_[0] = jets_btag[0].Px();
		integration[k].evBJet1_4P_[1] = jets_btag[0].Py();
		integration[k].evBJet1_4P_[2] = jets_btag[0].Pz();
		integration[k].evBJet1_4P_[3] = jets_btag[0].E();
		integration[k].evBJet2_4P_[0] = jets_btag[1].Px();
		integration[k].evBJet2_4P_[1] = jets_btag[1].Py();
		integration[k].evBJet2_4P_[2] = jets_btag[1].Pz();
		integration[k].evBJet2_4P_[3] = jets_btag[1].E();
		// untagged jets
		integration[k].n_lightJets_ = min(10,int(jets_untag.size()) );
		for( int j=0; j<integration[k].n_lightJets_; j++){
			integration[k].evJets_4P_[j][0] = jets_untag[j].Px();
			integration[k].evJets_4P_[j][1] = jets_untag[j].Py();
			integration[k].evJets_4P_[j][2] = jets_untag[j].Pz();
			integration[k].evJets_4P_[j][3] = jets_untag[j].E();
		}
		if (jets_untag.size()<10) {
			for (unsigned int j=jets_untag.size(); j<10; j++) {
				integration[k].evJets_4P_[j][0] = 0;
				integration[k].evJets_4P_[j][1] = 0;
				integration[k].evJets_4P_[j][2] = 0;
				integration[k].evJets_4P_[j][3] = 0;
			}
		}

		if (jets_untag.size()<=1) {
			integration[k].integration_type_ = 1;
			integration[k].evJet1_4P_[0] = 0;
			integration[k].evJet1_4P_[1] = 0;
			integration[k].evJet1_4P_[2] = 0;
			integration[k].evJet1_4P_[3] = 0;
			integration[k].evJet2_4P_[0] = 0;
			integration[k].evJet2_4P_[1] = 0;
			integration[k].evJet2_4P_[2] = 0;
			integration[k].evJet2_4P_[3] = 0;
		}
		else {
			TLorentzVector jetw1 = jets_untag[0];
			TLorentzVector jetw2 = jets_untag[1];
			
			double mW = 80.4;
			double delta_Mjj = fabs( (jetw1+jetw2).M() - mW );
			
			for (unsigned int j1=0; j1 < jets_untag.size()-1; j1++) {
				for (unsigned int j2=j1+1; j2 < jets_untag.size(); j2++) {
					if (fabs((jets_untag[j1]+jets_untag[j2]).M()-mW)<delta_Mjj) {
						jetw1 = jets_untag[j1];
						jetw2 = jets_untag[j2];
						delta_Mjj = fabs((jetw1+jetw2).M() - mW);
					}
				}
			}

			if ((jetw1+jetw2).M()>60 and (jetw1+jetw2).M()<100) {
				integration[k].integration_type_ = 0;
				integration[k].evJet1_4P_[0] = jetw1.Px();
				integration[k].evJet1_4P_[1] = jetw1.Py();
				integration[k].evJet1_4P_[2] = jetw1.Pz();
				integration[k].evJet1_4P_[3] = jetw1.E();
				integration[k].evJet2_4P_[0] = jetw2.Px();
				integration[k].evJet2_4P_[1] = jetw2.Py();
				integration[k].evJet2_4P_[2] = jetw2.Pz();
				integration[k].evJet2_4P_[3] = jetw2.E();
			}
			else {
				integration[k].integration_type_ = 1;
				integration[k].evJet1_4P_[0] = 0;
				integration[k].evJet1_4P_[1] = 0;
				integration[k].evJet1_4P_[2] = 0;
				integration[k].evJet1_4P_[3] = 0;
				integration[k].evJet2_4P_[0] = 0;
				integration[k].evJet2_4P_[1] = 0;
				integration[k].evJet2_4P_[2] = 0;
				integration[k].evJet2_4P_[3] = 0;
			}
		}

		integration[k].evRecoMET4P_[0] = met.Px();
		integration[k].evRecoMET4P_[1] = met.Py();
		integration[k].evRecoMET4P_[2] = met.Pz();
		integration[k].evRecoMET4P_[3] = met.E();

		// MET covariance Matrix component (index order 11, 01, 10, 00 )
		double det = evNtuple_.METCov00 * evNtuple_.METCov11 - evNtuple_.METCov01 * evNtuple_.METCov10;
		integration[k].evV_[0] =   evNtuple_.METCov11 / det;
		integration[k].evV_[1] = - evNtuple_.METCov01 / det;
		integration[k].evV_[2] = - evNtuple_.METCov10 / det;
		integration[k].evV_[3] =   evNtuple_.METCov00 / det;
		
		k++;
		
	} // end of event loop
}
