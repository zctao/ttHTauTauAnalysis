#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/EventSelector.h"

// member functions
bool EventSelector::pass_lepton_number(const std::vector<miniLepton>& lep_fakeable,
									   const std::vector<miniLepton>& lep_tight)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	// at least 2 fakeable leptons
	// no more than 2 tight leptons
	
	if (lep_fakeable.size() >= 2 and lep_tight.size() <= 2)
		return true;
	
	if (debug_) {
		std::cout << "FAIL lepton number requirement" << std::endl;
		std::cout << "nLepFakeable : " << lep_fakeable.size() << std::endl;
		std::cout << "nLepTight : " << lep_tight.size() << std::endl;
	}
	return false;
	
}

bool EventSelector::pass_tau_number(int ntaus)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);
	
	// at least 1 selected tau	
	if (ntaus >= 1)
		return true;
	
	if (debug_) {
		std::cout << "FAIL tau number requirement" << std::endl;
		std::cout << "nTau : " << ntaus << std::endl;
	}
	return false;
}

bool EventSelector::pass_lepton_pt(const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);
	
	assert(leps.size() >= 2);
	// lepton pt
	float minpt_ldg = 25.;
	float minpt_subldg = abs(leps[1].pdgId())==11 ? 15. : 10.;

	if (leps[0].pt() > minpt_ldg and leps[1].pt() > minpt_subldg)
		return true;

	if (debug_) {
		std::cout << "FAIL lepton pT cut" << std::endl;
		std::cout << "lep0 pt id : " << leps[0].pt() << " " << leps[0].pdgId() << std::endl;
		std::cout << "lep1 pt id : " << leps[1].pt() << " " << leps[1].pdgId() << std::endl;
	}
	return false;
}

bool EventSelector::pass_pairMass_veto(const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);
	
	// veto two leptons with invariant mass  < 12 GeV
	for (auto it = leps.begin(); it != leps.end()-1; ++it) {
		for (auto it2 = it+1; it2 != leps.end(); ++it2) {
	if ( (it->p4() + it2->p4()).M() < 12. ) {
		if (debug_) {
			std::cout << "FAIL any pair of loose leptons has invariant mass >= 12 GeV" << std::endl;
			std::cout << "pt eta phi mass : " << it->pt() << " " << it->eta()
					  << " " << it->phi() << " " << it->mass() << std::endl;
			std::cout << "pt eta phi mass : " << it2->pt() << " " << it2->eta()
					  << " " << it2->phi() << " " << it2->mass() << std::endl;
		}
		return false;
	}
		} // end of inner loop
	} // end of outerloop

	return true;
}

bool EventSelector::pass_tight_charge(const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	// tight charge
	if (leps[0].passTightCharge() and leps[1].passTightCharge())
		return true;

	if (debug_) {
		std::cout << "FAIL tight charge" << std::endl;
		std::cout << "lep0 " << leps[0].pdgId() << leps[0].passTightCharge()
				  << std::endl;
		std::cout << "lep1 " << leps[1].pdgId() << leps[1].passTightCharge()
				  << std::endl;
	}
	return false;
}

bool EventSelector::pass_Zmass_veto(const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	// Zmass veto: 91.2 +/- 10
	// ee only
	assert(leps.size() >= 2);
	if (not (abs(leps[0].pdgId())==11 and abs(leps[1].pdgId())==11) ) 
		return true;

	double eeInvMass = (leps[0].p4()+leps[1].p4()).M();

	if (eeInvMass < (91.2 - 10.0) or eeInvMass > (91.2 + 10.0))
		return true;

	if (debug_) {
		std::cout << "FAIL Z mass veto" << std::endl;
		std::cout << "ee pair invariant mass : " << eeInvMass << std::endl;
	}
	return false;
}

bool EventSelector::pass_metLD(float metLD, const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	// metLD cut
	// ee only
	assert(leps.size() >= 2);
	if (not (abs(leps[0].pdgId())==11 and abs(leps[1].pdgId())==11) ) 
		return true;

	if (metLD > 0.2)
		return true;

	if (debug_) {
		std::cout << "FAIL metLD cut" << std::endl;
		std::cout << "metLD : " << metLD << std::endl;
	}
	return false;
}

bool EventSelector::pass_jet_number(int njets)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	if (njets >= 3)
		return true;

	if (debug_) {
		std::cout << "FAIL number of jets requirement" << std::endl;
		std::cout << "njets : " << njets << std::endl;
	}
	return false;
}

bool EventSelector::pass_btag_number(int nbtags_loose, int nbtags_medium)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	if (nbtags_loose >= 2 or nbtags_medium >= 1)
		return true;

	if (debug_) {
		std::cout << "FAIL number of btags requirement" << std::endl;
		std::cout << "loose btags : " << nbtags_loose << std::endl;
		std::cout << "medium btags : " << nbtags_medium << std::endl;
	}
	return false;
}

bool EventSelector::pass_lepton_charge(int lep0Charge, int lep1Charge)
{
	// 2lss1tau category only
	assert(anaType_==Analyze_2lss1tau);

	bool samesign = lep0Charge * lep1Charge > 0;
	bool pass = (selType_==Control_2los1tau) ? (!samesign) : samesign;	
	if (pass)
		return true;

	if (debug_) {
		std::cout << "FAIL lepton charge requirement" << std::endl;
		std::cout << "lep charges : " << lep0Charge << " " << lep1Charge
				  << std::endl;
	}
	return false;

}

bool EventSelector::pass_tau_charge(int tauCharge,
									const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only
	assert(anaType_==Analyze_2lss1tau);

	assert(leps.size() >= 2);

	if (selType_==Control_2los1tau) {
		assert(leps[0].charge()*leps[1].charge() < 0);
		if (leps[0].charge()==tauCharge and abs(leps[0].pdgId())==11)
			return true;
		if (leps[1].charge()==tauCharge and abs(leps[1].pdgId())==11)
			return true;
	}
	else {
		assert(leps[0].charge()*leps[1].charge() > 0);
		if (tauCharge * leps[0].charge() < 0)
			return true;
	}

	if (debug_) {
		std::cout << "FAIL tau charge requirement" << std::endl;
		std::cout << "tau charge : " << tauCharge << std::endl;
		std::cout << "lep charges : " << leps[0].charge() << " "
				  << leps[1].charge() << std::endl;
	}
	return false;
}

bool EventSelector::pass_lepton_ID(bool lep0IsTight, bool lep1IsTight)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);
	
	// signal region: two leading leptons are tight
	bool passLepWP = lep0IsTight and lep1IsTight;

	if (selType_==Control_1lfakeable)
		passLepWP = not passLepWP;

	if (passLepWP)
		return true;

	if (debug_) {
		std::cout << "FAIL lepton WP requirement" << std::endl;
		std::cout << "passtight?  " << lep0IsTight << " " << lep1IsTight
				  << std::endl;
	}
	return false;
}

bool EventSelector::pass_lep_mc_match(const std::vector<miniLepton>& leps)
{
	// only for MC samples
	
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	assert(leps.size() >= 2);

	bool matchGenLeps = true;
	int lep_cnt = 0;
	
	for (const auto & l : leps) {
		// for electron
		int isPrompt = abs(l.pdgId())==11 ? 1 : 2;
		int isPromptTauDecay = abs(l.pdgId())==11 ? 3 : 4;
		
		if (not (l.MCMatchType()==isPrompt or l.MCMatchType()==isPromptTauDecay))
			matchGenLeps = false;

		// match only the two leading leptons
		if (++lep_cnt > 1) break;
	}

	if (matchGenLeps)
		return true;

	if (debug_) {
		std::cout << "FAIL lepton MC match" << std::endl;
		std::cout << "MC match type : " << leps[0].MCMatchType() << " "
				  << leps[1].MCMatchType() << std::endl;
	}
	return false;
}

bool EventSelector::pass_tau_mc_match(const pat::Tau& tau)
{
	// only for MC samples

	int mtype = -1;
	if (tau.hasUserInt("MCMatchType"))
		mtype = tau.userInt("MCMatchType");
	else {
		std::cerr << "ERROR: MC matching has not been done yet!" << std::endl;
		return false;
	}

	if (mtype==1 or mtype==2 or mtype==3 or mtype==4 or mtype==5)
		return true;

	if (debug_) {
		std::cout << "FAIL tau MC match" << std::endl;
		std::cout << "MC match type : " << mtype << std::endl;
	}
	return false;
	
}
