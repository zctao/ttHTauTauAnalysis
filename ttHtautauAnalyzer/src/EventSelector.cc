#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/EventSelector.h"

// member functions
bool EventSelector::pass_lepton_number(const std::vector<miniLepton>& lep_fakeable,
									   const std::vector<miniLepton>& lep_tight)
{
	if (debug_) {
		std::cout << "nLepFakeable : " << lep_fakeable.size() << std::endl;
		std::cout << "nLepTight : " << lep_tight.size() << std::endl;
	}

	if (anaType_==Analyze_2lss1tau) {
		// at least 2 fakeable leptons and no more than 2 tight leptons
		if (lep_fakeable.size() >= 2 and lep_tight.size() <= 2)
			return true;
	}
	else if (anaType_==Analyze_1l2tau) {
		// at least 1 fakeable lepton and no more than 1 tight lepton
		if (lep_fakeable.size() >= 1 and lep_tight.size() <= 1)
			return true;
	}
	//else if (anaType_==Analyze_3l1tau) {
	//
	//}
	else
		std::cout << "Analysis type not available!!" << std::endl;
	
	if (debug_) {
		std::cout << "FAIL lepton number requirement" << std::endl;
	}
	return false;
	
}

bool EventSelector::pass_lepton_pt(const std::vector<miniLepton>& leps)
{
	if (debug_) {
		for (const auto& lep : leps)
			std::cout << "lep pt id" << lep.pt()<<" "<< lep.pdgId() << std::endl;
	}
	
	if (anaType_==Analyze_2lss1tau) {
		assert(leps.size() >= 2);
		// lepton pt
		float minpt_ldg = 25.;
		float minpt_subldg = abs(leps[1].pdgId())==11 ? 15. : 10.;

		if (leps[0].pt() > minpt_ldg and leps[1].pt() > minpt_subldg)
			return true;
	}
	else if (anaType_==Analyze_1l2tau) {
		assert(leps.size() >= 1);
		// lepton pt
		float minpt = abs(leps[0].pdgId())==11 ? 30. : 25.;

		if (leps[0].pt() > minpt)
			return true;
	}
	//else if (anaType_==Analyze_3l1tau) {
	//
	//}
	else
		std::cout << "Analysis type not available!!" << std::endl;

	if (debug_) {
		std::cout << "FAIL lepton pT cut" << std::endl;
	}
	return false;
}

bool EventSelector::pass_lepton_charge(int lep0Charge, int lep1Charge)
{
	// 2lss1tau category only
	assert(anaType_==Analyze_2lss1tau);

	if (debug_) {
		std::cout << "lep charges : " << lep0Charge << " " << lep1Charge
				  << std::endl;
	}
	
	bool samesign = lep0Charge * lep1Charge > 0;
	bool pass = (selType_==Control_2los1tau) ? (!samesign) : samesign;	
	if (pass)
		return true;

	if (debug_) {
		std::cout << "FAIL lepton charge requirement" << std::endl;
	}
	return false;

}

bool EventSelector::pass_tight_charge(const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	// tight charge
	if (debug_) {
		std::cout << "lep0 id tightCharge : " << leps[0].pdgId() << " "
				  << leps[0].passTightCharge()
				  << std::endl;
		std::cout << "lep1 id tightCharge : " << leps[1].pdgId() << " "
				  << leps[1].passTightCharge()
				  << std::endl;
	}
	
	if (leps[0].passTightCharge() and leps[1].passTightCharge())
		return true;

	if (debug_) {
		std::cout << "FAIL tight charge" << std::endl;
	}
	return false;
}

bool EventSelector::pass_pairMass_veto(const std::vector<miniLepton>& leps)
{
	assert(leps.size()>0);
	
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

	if (debug_) {
		std::cout << "ee pair invariant mass : " << eeInvMass << std::endl;
	}
	
	if (eeInvMass < (91.2 - 10.0) or eeInvMass > (91.2 + 10.0))
		return true;

	if (debug_) {
		std::cout << "FAIL Z mass veto" << std::endl;
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

	if (debug_)
		std::cout << "metLD : " << metLD << std::endl;
	
	if (metLD > 0.2)
		return true;

	if (debug_) {
		std::cout << "FAIL metLD cut" << std::endl;
	}
	return false;
}

bool EventSelector::pass_lepton_ID(bool lep0IsTight, bool lep1IsTight,
								   bool lep2IsTight)
{
	if (debug_) {
		std::cout << "lepID passtight?  " << lep0IsTight;
		if (anaType_!=Analyze_1l2tau)
			std::cout << " " << lep1IsTight;
		if (anaType_==Analyze_3l1tau)
			std::cout << " " << lep2IsTight;
		std::cout << std::endl;
	}

	if (anaType_==Analyze_2lss1tau) {
		// signal region: two leading leptons are tight
		bool passLepWP = lep0IsTight and lep1IsTight;

		if (selType_==Control_fake_2lss1tau) passLepWP = not passLepWP;
		if (selType_==Loose_2lss1tau) passLepWP = true;

		if (passLepWP) return true;
	}
	else if (anaType_==Analyze_1l2tau) {
		// signal region: leading lepton is tight
		bool passLepWP = lep0IsTight;

		//if (selType_==Control_fake_1l2tau) passLepWP = not lep0IsTight;

		if (passLepWP) return true;
	}
	//else if (anaType_==Analyze_3l1tau) {
	//
	//}
	else
		std::cout << "Analysis type not available!!" << std::endl;	

	if (debug_) {
		std::cout << "FAIL lepton WP requirement" << std::endl;
	}
	return false;
}

bool EventSelector::pass_tau_number(int ntaus)
{
	if (debug_) {
		std::cout << "nTau : " << ntaus << std::endl;
	}

	if (anaType_==Analyze_2lss1tau) {
		// at least 1 selected tau
		if (ntaus >= 1) return true;
	}
	else if (anaType_==Analyze_1l2tau) {
		if (ntaus >= 2) return true;
	}
	//else if (anaType_==Analyze_3l1tau) {
	//
	//}
	else
		std::cout << "Analysis type not available!!" << std::endl;	
	
	if (debug_) {
		std::cout << "FAIL tau number requirement" << std::endl;
	}
	return false;
}

bool EventSelector::pass_tau_charge(int tauCharge,
									const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only
	assert(anaType_==Analyze_2lss1tau);

	assert(leps.size() >= 2);

	if (debug_) {
		std::cout << "tau charge : " << tauCharge << std::endl;
		std::cout << "lep charges : " << leps[0].charge() << " "
				  << leps[1].charge() << std::endl;
	}
	
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
	}
	return false;
}

bool EventSelector::pass_taupair_charge(int tau0charge, int tau1charge)
{
	assert(anaType_==Analyze_1l2tau);
	assert(selType_==Signal_1l2tau or selType_==Control_fake_1l2tau);
	
	if (debug_) {
		std::cout << "tau charge : " << tau0charge << " " << tau1charge
				  << std::endl;
	}
		
	bool opposite = tau0charge * tau1charge < 0;

	if (opposite) return true;
	
	if (debug_) {
		std::cout << "FAIL tau pair charge requirement" << std::endl;
	}
	return false;
	
}

bool EventSelector::pass_tau_ID(int ntau_tight)
{
	// for 1l2tau only
	assert(anaType_==Analyze_1l2tau);

	if (debug_) {
		std::cout << "number of tight taus : " << ntau_tight << std::endl;
	}

	bool passTauWP = ntau_tight >= 2;

	if (passTauWP) return true;

	if (debug_) {
		std::cout << "FAIL tau WP requirement" << std::endl;
	}
	return false;
}

bool EventSelector::pass_lep_tau_ID(bool lepIsTight, int ntau_tight)
{
	// for 1l2tau only
	assert(anaType_==Analyze_1l2tau);

	// signal region: the leading lepton is tight
	// and at least 2 taus are tight
	bool passWP = pass_lepton_ID(lepIsTight) and pass_tau_ID(ntau_tight);

	if (selType_==Control_fake_1l2tau)
		passWP = not passWP;

	return passWP;
}

bool EventSelector::pass_jet_number(int njets)
{
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
	if (debug_) {
		std::cout << "loose btags : " << nbtags_loose << std::endl;
		std::cout << "medium btags : " << nbtags_medium << std::endl;
	}
	
	if (nbtags_loose >= 2 or nbtags_medium >= 1)
		return true;

	if (debug_) {
		std::cout << "FAIL number of btags requirement" << std::endl;
	}
	return false;
}

bool EventSelector::pass_lep_mc_match(const miniLepton& lep)
{
	// for MC only

	if (debug_)
		std::cout << "MC match type : " << lep.MCMatchType() << std::endl;

	int isPrompt = abs(lep.pdgId())==11 ? 1 : 2;
	int isPromptTauDecay = abs(lep.pdgId())==11 ? 3 : 4;

	if (lep.MCMatchType()==isPrompt or lep.MCMatchType()==isPromptTauDecay)
		return true;

	if (debug_) {
		std::cout << "FAIL lepton MC match" << std::endl;
	}
	return false;
}

bool EventSelector::pass_lep_mc_match(const std::vector<miniLepton>& leps,
									  int nleps)
{	
	bool matchGenLeps = true;
	int lep_cnt = 0;
	
	for (const auto & l : leps) {
		
		if (not pass_lep_mc_match(l))
			matchGenLeps = false;

		// match only the leading leptons
		if (++lep_cnt >= nleps) break;
	}

	if (matchGenLeps)
		return true;

	if (debug_) {
		std::cout << "FAIL lepton MC match" << std::endl;
	}
	return false;
}

bool EventSelector::pass_tau_mc_match(const pat::Tau& tau)
{
	// for MC only

	int mtype = -1;
	if (tau.hasUserInt("MCMatchType"))
		mtype = tau.userInt("MCMatchType");
	else {
		std::cerr << "ERROR: MC matching has not been done yet!" << std::endl;
		return false;
	}

	if (debug_) {
		std::cout << "MC match type : " << mtype << std::endl;
	}
	
	if (mtype==1 or mtype==2 or mtype==3 or mtype==4 or mtype==5)
		return true;
	
	if (debug_) {
		std::cout << "FAIL tau MC match" << std::endl;
	}
	return false;
}

bool EventSelector::pass_tau_mc_match(const std::vector<pat::Tau>& taus,
									  int ntaus)
{
	bool matchGenTaus = true;
	int tau_cnt = 0;

	for (const auto & tau : taus) {
		
		if (not pass_tau_mc_match(tau))
			matchGenTaus = false;
		
		// match only the leading taus
		if (++tau_cnt >= ntaus) break;
	}

	if (matchGenTaus)
		return true;

	if (debug_) {
		std::cout << "FAIL tau MC match" << std::endl;
	}
	return false;
}
