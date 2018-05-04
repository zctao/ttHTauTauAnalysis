#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/EventSelector.h"

// member functions

void EventSelector::fill_cutflow(TH1* h, int ibin, const char* name)
{
	assert(h);
	assert(ibin>0);

	if (ibin > h->GetNbinsX()) {
		std::cout << "WARNING : ibin " << ibin
				  << " exceeds number of bins of histogram" << std::endl;
		return;
	}
	
	if (h->GetBinContent(ibin)==0) { // first time filling this bin
		// set bin label first
		h->GetXaxis()->SetBinLabel(ibin, name);
	}
	
	//assert(h->GetXaxis()->GetBinLabel(ibin)).c_str()==name.c_str());
	h->AddBinContent(ibin);

	return;
}

bool EventSelector::pass_hlt_paths(Analysis_types anatype,
								   TriggerHelper* const trighelper,
								   unsigned int triggerBits)
{
	bool passhlt = false;
	if (anatype == Analyze_1l2tau) {
		passhlt =
			trighelper->pass_single_lep_triggers(triggerBits) or
			trighelper->pass_leptau_cross_triggers(triggerBits);
	}
	else if (anatype == Analyze_2lss1tau) {
		passhlt =
			trighelper->pass_single_lep_triggers(triggerBits) or
			trighelper->pass_dilep_triggers(triggerBits);
	}
	else if (anatype == Analyze_3l1tau) {
		passhlt =
			trighelper->pass_single_lep_triggers(triggerBits) or
			trighelper->pass_dilep_triggers(triggerBits) or
			trighelper->pass_trilep_triggers(triggerBits);
	}
	else if (anatype == Analyze_2l2tau) {
		passhlt =
			trighelper->pass_single_lep_triggers(triggerBits) or
			trighelper->pass_dilep_triggers(triggerBits) or
			trighelper->pass_leptau_cross_triggers(triggerBits);
	}

	if (not passhlt and verbose_) {
		std::cout << "FAIL HLT" << std::endl;
		std::cout << "trigger bits : " << triggerBits << std::endl;
	}
	
	return passhlt;
}

bool EventSelector::pass_hlt_match(Analysis_types anatype,
								   TriggerHelper* const trighelper,
								   unsigned int triggerBits,
								   int nElectron, int nMuon)
{
	bool pass = false;

	bool pass_e = trighelper->pass_single_e_triggers(triggerBits);
	bool pass_m = trighelper->pass_single_m_triggers(triggerBits);
	bool pass_2e = trighelper->pass_2e_triggers(triggerBits);
	bool pass_2m = trighelper->pass_2m_triggers(triggerBits);
	bool pass_em = trighelper->pass_em_triggers(triggerBits);
	bool pass_etau = trighelper->pass_etau_triggers(triggerBits);
	bool pass_mtau = trighelper->pass_mtau_triggers(triggerBits);
	
	if (anatype == Analyze_1l2tau) {		
		pass = ((pass_e or pass_etau) and nElectron>0) or
			((pass_m or pass_mtau) and nMuon>0);
	}
	else if (anatype == Analyze_2lss1tau) {
		pass = (pass_e and nElectron>0) or (pass_2e and nElectron>1) or
			(pass_m and nMuon>0) or (pass_2m and nMuon>1) or
			(pass_em and nElectron>0 and nMuon>0);
	}
	else if (anatype == Analyze_3l1tau) {
		bool pass_3e = trighelper->pass_3e_triggers(triggerBits);
		bool pass_m2e = trighelper->pass_m_2e_triggers(triggerBits);
		bool pass_2me = trighelper->pass_2m_e_triggers(triggerBits);
		bool pass_3m = trighelper->pass_3m_triggers(triggerBits);
		pass = (pass_e and nElectron>0) or (pass_2e and nElectron>1) or
			(pass_m and nMuon>0) or (pass_2m and nMuon>1) or
			(pass_em and nElectron>0 and nMuon>0) or
			(pass_3e and nElectron>2) or (pass_3m and nMuon>2) or
			(pass_m2e and nElectron>1 and nMuon>0) or
			(pass_2me and nElectron>0 and nMuon>1);
	}
	else if (anatype == Analyze_2l2tau) {
		pass = ((pass_e or pass_etau) and nElectron>0) or
			((pass_m or pass_mtau) and nMuon>0) or
			(pass_2e and nElectron>1) or (pass_2m and nMuon>1) or
			(pass_em and nElectron>0 and nMuon>0);
	}

	if (not pass and verbose_) {
		std::cout << "FAIL to match the number of offline objects to HLT paths" << std::endl;
		std::cout << "trigger bits : " << triggerBits << std::endl;
		std::cout << "nElectron = " << nElectron << " ";
		std::cout << "nMuon = " << nMuon << std::endl;
	}
	
	return pass;
}

/////////////////////////////////
// ttH l+tau inclusive
/////////////////////////////////
bool EventSelector::pass_ttH_ltau_inclusive_selection(
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& fakeableTaus,
	int njets, TH1* h_cutflow)
{
	if (verbose_) std::cout << "start inclusive event selection" << std::endl;
	
	int ibin = 1;
	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");

	/////////////////////////////////
	// lepton and tau number
	if (verbose_) {
		std::cout << "nFakeableLeps = " << fakeableLeps.size() << std::endl;
		std::cout << "nFakeableTaus = " << fakeableTaus.size() << std::endl;
	}
	// nlep + ntau >= 3
	if (fakeableLeps.size()+fakeableTaus.size() >= 3) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "nlep+ntau>2");
	}
	else {
		if (verbose_) std::cout << "FAIL nlepton+ntau >= 3" << std::endl;
		return false;
	}

	// nlep >=1
	if (fakeableLeps.size() > 0) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "nlep>0");
	}
	else {
		if (verbose_) std::cout << "FAIL nlep >= 1" << std::endl;
		return false;
	}

	/*
	// leading lep pt > 20
	assert(fakeableLeps.size() > 0);
	if (fakeableLeps[0].conept() >= 20.) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep0_pt>20");
	}
	else {
		if (verbose_) std::cout << "FAIL leading lepton pt cut" << std::endl;
		return false;
	}
	*/

	/////////////////////////////////
	// At least 2 jets
	if (verbose_) std::cout << "nJets = " << njets << std::endl;
	if (njets >= 2) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "njet>1");
	}
	else {
		if (verbose_) std::cout << "FAIL jet number requirement" << std::endl;
		return false;
	}

	/*
	/////////////////////////////////
	// Dilepton mass of any loose lepton pair > 12 GeV
	bool passMll = pass_pairMass_veto(looseLeps);
	if (passMll) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Mll>12GeV");
	}
	else
		return false;
	*/

	
	/////////////////////////////////
	if (verbose_) std::cout << "PASSED ttH inclusive event selection!" << std::endl;

	return true;
}

/////////////////////////////////
bool EventSelector::pass_extra_event_selection(
    Analysis_types anatype, Selection_types seltype,
    std::vector<miniLepton> const * const fakeableLeps,
	std::vector<miniTau> const * const selectedTaus,
	std::vector<miniTau> const * const fakeableTaus)
{
	// assume already passed inclusive selection

	bool pass=false;

	if (verbose_) std::cout << std::endl;
	
	if (anatype == Analyze_1l2tau) {		
		if (seltype == Signal_1l2tau)
			pass = pass_1l2tau_SR_selection(*fakeableLeps, *selectedTaus);
		else if (seltype == Control_fake_1l2tau)
			pass = pass_1l2tau_FakeAR_selection(*fakeableLeps, *selectedTaus,
												*fakeableTaus);
	}
	else if (anatype == Analyze_2lss1tau) {
		if (seltype == Signal_2lss1tau)
			pass = pass_2lss1tau_SR_selection(*fakeableLeps, *selectedTaus);
		else if (seltype == Control_fake_2lss1tau)
			pass = pass_2lss1tau_FakeAR_selection(*fakeableLeps, *selectedTaus);
		else if (seltype == Control_2los1tau)
			pass = pass_2lss1tau_FlipAR_selection(*fakeableLeps, *selectedTaus);
    }
	else if (anatype == Analyze_3l1tau) {
		if (seltype == Signal_3l1tau)
			pass = pass_3l1tau_SR_selection(*fakeableLeps, *selectedTaus);
		else if (seltype == Control_fake_3l1tau)
			pass = pass_3l1tau_FakeAR_selection(*fakeableLeps, *selectedTaus);
	}
	else if (anatype == Analyze_2l2tau) {
		if (seltype == Signal_2l2tau)
			pass = pass_2l2tau_SR_selection(*fakeableLeps, *selectedTaus);
		else if (seltype == Control_fake_2l2tau)
			pass = pass_2l2tau_FakeAR_selection(*fakeableLeps, *selectedTaus,
												*fakeableTaus);
	}
	
	return pass;
}

/////////////////////////////////
bool EventSelector::pass_full_event_selection(
    Analysis_types anatype, Selection_types seltype,
	const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniLepton>& tightLeps,
	const std::vector<miniTau>& fakeableTaus,
	const std::vector<miniTau>& selectedTaus,
	int njets, int nbtags_loose, int nbtags_medium, float metLD,
	TH1* h_cutflow)
{
	bool pass=false;
	if (verbose_) std::cout << std::endl;

	if (anatype == Analyze_1l2tau) {
		if (not pass_1l2tau_inclusive_selection(looseLeps, fakeableLeps, tightLeps,
												fakeableTaus, njets, nbtags_loose,
												nbtags_medium, h_cutflow) )
			return false;

		if (seltype == Signal_1l2tau)
			pass = pass_1l2tau_SR_selection(fakeableLeps, selectedTaus);
		else if (seltype == Control_fake_1l2tau)
			pass = pass_1l2tau_FakeAR_selection(fakeableLeps, selectedTaus,
												fakeableTaus);
	}
	else if (anatype == Analyze_2lss1tau) {
		if (seltype == Control_ttW) {
			pass = pass_ttW_CR_selection(looseLeps, fakeableLeps, tightLeps,
										 njets, nbtags_loose, nbtags_medium,
										 metLD, h_cutflow);
		}
		else {
			if (not pass_2l1tau_inclusive_selection(looseLeps, fakeableLeps,
													tightLeps, fakeableTaus, njets,
													nbtags_loose, nbtags_medium,
													metLD, h_cutflow) )
				return false;
			
			if (seltype == Signal_2lss1tau)
				pass = pass_2lss1tau_SR_selection(fakeableLeps, selectedTaus);
			else if (seltype == Control_fake_2lss1tau)
				pass = pass_2lss1tau_FakeAR_selection(fakeableLeps, selectedTaus);
			else if (seltype == Control_2los1tau)
				pass = pass_2lss1tau_FlipAR_selection(fakeableLeps, selectedTaus);
		}
	}
	else if (anatype == Analyze_3l1tau) {
		if (seltype == Control_ttZ) {
			pass = pass_ttZ_CR_selection(looseLeps, fakeableLeps, tightLeps,
										 njets, nbtags_medium, metLD, h_cutflow);
		}
		else {	
			if (not pass_3l1tau_inclusive_selection(looseLeps, fakeableLeps,
													fakeableTaus, njets,
													nbtags_loose, nbtags_medium,
													metLD, h_cutflow) )
				return false;

			if (seltype == Signal_3l1tau)
				pass = pass_3l1tau_SR_selection(fakeableLeps, selectedTaus);
			else if (seltype == Control_fake_3l1tau)
				pass = pass_3l1tau_FakeAR_selection(fakeableLeps, selectedTaus);
		}
	}
	else if (anatype == Analyze_2l2tau) {
		if (not pass_2l2tau_inclusive_selection(looseLeps, fakeableLeps,
												fakeableTaus, njets, nbtags_loose,
												nbtags_medium, metLD, h_cutflow) )
			return false;

		if (seltype == Signal_2l2tau)
			pass = pass_2l2tau_SR_selection(fakeableLeps, selectedTaus);
		else if (seltype == Control_fake_2l2tau)
			pass = pass_2l2tau_FakeAR_selection(fakeableLeps, selectedTaus,
												fakeableTaus);
	}
	else {
		std::cout << "Unsupported analysis type. Return false." << std::endl;
		return false;
	}

	return pass;
}

/////////////////////////////////
// 1l2tau
/////////////////////////////////
bool EventSelector::pass_1l2tau_inclusive_selection(
	const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniLepton>& tightLeps,
	const std::vector<miniTau>& fakeableTaus,
	int njets, int nbtags_loose, int nbtags_medium,
	TH1* h_cutflow)
{	
	if (verbose_) std::cout << "start event selection: 1l2tau" << std::endl;

	int ibin = 1;
	
	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");
	
	/////////////////////////////////
	// At least 1 fakeable lepton
	// and no more than 1 tight lepton
	if (verbose_) {
		std::cout << "nFakeableLeptons = "<< fakeableLeps.size()
				  << "  nTightLeptons = " << tightLeps.size() << std::endl;
	}
	//
	bool passLepNumber = fakeableLeps.size() > 0 and tightLeps.size() < 2;
	if (passLepNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep num");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// lepton pt
	if (verbose_) {
		std::cout << "pdgID conept pt : " << fakeableLeps[0].pdgId() << " "
				  << fakeableLeps[0].conept() << " " << fakeableLeps[0].pt()
				  << std::endl;
	}
	
	float minpt = abs(fakeableLeps[0].pdgId())==11 ? 30. : 25.;
	bool passLepPt = fakeableLeps[0].conept() >= minpt;
	if (passLepPt) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep pt");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton pT requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// lepton eta (to match trigger)
	if (verbose_) {
		std::cout << "pdgID eta : " << fakeableLeps[0].pdgId() << " "
				  << fakeableLeps[0].eta() << std::endl;
	}

	bool passLepEta = fabs(fakeableLeps[0].eta()) < 2.1;
	if (passLepEta) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep eta");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton eta requirement" << std::endl;
		return false;
	}
	
	/////////////////////////////////
	// At least 2 fakeable taus
	if (verbose_) std::cout << "nFakeableTaus = " << fakeableTaus.size() << std::endl;
	
	bool passTauNumber = fakeableTaus.size() >= 2;
	if (passTauNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tau num");
	}
	else {
		if (verbose_) std::cout << "FAIL tau number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// tau pt
	if (verbose_) {
		std::cout << "Two leading tau pT : " << fakeableTaus[0].pt() << " "
				  << fakeableTaus[1].pt() << std::endl;
	}
	
	bool passTauPt = fakeableTaus[0].pt() >= 30. and fakeableTaus[1].pt() >= 20.;  // FIXME
	if (passTauPt) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tau pt");
	}
	else {
		if (verbose_) std::cout << "FAIL tau pT requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// At least 3 selected jets
	if (verbose_) std::cout << "nJets = " << njets << std::endl;
	
	bool passJetNumber = njets >= 3;
	if (passJetNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "njets>=3");
	}
	else {
		if (verbose_) std::cout << "FAIL jet number requirement" << std::endl;
		return false;
	}
	
	/////////////////////////////////
	// At least 2 loose btag or 1 medium btag
	if (verbose_) {
		std::cout << "nbtags loose : " << nbtags_loose << std::endl;
		std::cout << "nbtags medium : " << nbtags_medium << std::endl;
	}
	
	bool passBTagNumber = (nbtags_loose >= 2) or (nbtags_medium >= 1);
	if (passBTagNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "btag num");
	}
	else {
		if (verbose_) std::cout << "FAIL btag number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// Dilepton mass of any loose lepton pair > 12 GeV
	bool passMll = pass_pairMass_veto(looseLeps);
	if (passMll) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Mll>12GeV");
	}
	else
		return false;
	
	/////////////////////////////////
	if (verbose_)
		std::cout << "PASSED inclusive 1l2tau event selection!" << std::endl;

	return true;
}

bool EventSelector::pass_1l2tau_tightID(const std::vector<miniLepton>& fakeableLeps,
										const std::vector<miniTau>& tightTaus)
{	
	// SR: tight lepton; both taus pass tight selection (VTight MVA)
	if (verbose_) {
		std::cout << "number of tight taus : " << tightTaus.size() << std::endl;
		std::cout << "lep istight : " << fakeableLeps[0].passTightSel() << std::endl;
	}
	
	if (tightTaus.size()<2) {
		return false;
	}

	assert(fakeableLeps.size()>0);
	if (not fakeableLeps[0].passTightSel())
		return false;

	return true;
	
	
	//if (verbose_) {
	//	std::cout << "IsTight lep tau1 tau2 : " << fakeableLeps[0].passTightSel()
	//			  << " " << selectedTaus[0].passTightSel()
	//			  << " " << selectedTaus[1].passTightSel() << std::endl;
	//	std::cout << "tauID WP tau1 tau2 : " << selectedTaus[0].tauIDMVAWPindex()
	//			  << " " << selectedTaus[1].tauIDMVAWPindex() << std::endl;
	//}
}

bool EventSelector::pass_1l2tau_charge(const std::vector<miniTau>& taus)
{
	assert(taus.size()>1);

	if (verbose_)
		std::cout << "tau charges : " << taus[0].charge() << " "
				  << taus[1].charge() << std::endl;

	// SR: the two taus have opposite charge
	bool passCharge = (taus[0].charge() * taus[1].charge() < 0);
	return passCharge;
}

bool EventSelector::pass_1l2tau_SR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& tightTaus)
{	
	// event is assumed to already pass 1l2tau inclusive selection
	if (verbose_) std::cout << "start 1l2tau signal region selection" << std::endl;
	
	///////////////////////////////
	// lepton and tau ID
	if (not pass_1l2tau_tightID(fakeableLeps, tightTaus) ) {
		if (verbose_) std::cout << "FAIL lepton and tau ID requirement" << std::endl;
		return false;
	}
		
	///////////////////////////////
	// tau charges
	assert(tightTaus.size()>1);
	if ( not pass_1l2tau_charge(tightTaus) ) {
		if (verbose_) std::cout << "FAIL tau charge requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// MC Matching
	if (isMC_) {
		assert(fakeableLeps.size()>0 and tightTaus.size()>1);
		bool passMCMatch = fakeableLeps[0].isGenMatched() and
			tightTaus[0].isGenMatched() and tightTaus[1].isGenMatched();
		if (not passMCMatch) {
			if (verbose_) {
				std::cout << "isGenMatched lep tau1 tau2 : "
						  << fakeableLeps[0].isGenMatched() << " "
						  << tightTaus[0].isGenMatched() << " "
						  << tightTaus[1].isGenMatched() << std::endl;
				std::cout << "mcMatchType lep tau1 tau2 : "
						  << fakeableLeps[0].MCMatchType() << " "
						  << tightTaus[0].MCMatchType() << " "
						  << tightTaus[1].MCMatchType() << std::endl;
				std::cout << "FAIL MC Matching" << std::endl;
			}
			return false;
		}
	}
	
	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 1l2tau signal region selection!" << std::endl;
	
	return true;
}

bool EventSelector::pass_1l2tau_FakeAR_selection(
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& tightTaus,
	const std::vector<miniTau>& fakeableTaus)
{
	assert(not looseselection_);
	
	if ( pass_1l2tau_tightID(fakeableLeps, tightTaus) ) {
		if (verbose_) std::cout << "FAIL lepton and tau ID WP" << std::endl;
		return false;
	}
	if ( not pass_1l2tau_charge(fakeableTaus) ) {
		if (verbose_) std::cout << "FAIL charge requirement" << std::endl;
		return false;
	}
	return true;
}

bool EventSelector::pass_1l2tau_CR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& tightTaus)
{
	assert(not looseselection_);
	
	if ( not pass_1l2tau_tightID(fakeableLeps, tightTaus) ) return false;
	if ( pass_1l2tau_charge(tightTaus) ) return false;
	return true;

	// MC match?
}

bool EventSelector::pass_1l2tau_FakeARCR_selection(
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& tightTaus,
	const std::vector<miniTau>& fakeableTaus)
{
	assert(not looseselection_);
	
	if ( pass_1l2tau_tightID(fakeableLeps, tightTaus) ) return false;
	if ( pass_1l2tau_charge(fakeableTaus) ) return false;
    return true;
}

/////////////////////////////////
// 2l1tau
/////////////////////////////////

bool EventSelector::pass_2l_generic_selection(
    const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniLepton>& tightLeps,
	int njets, int nbtags_loose, int nbtags_medium, float metLD,
    int& ibin, TH1* h_cutflow)
{
	if (verbose_) std::cout << "start event selection: generic 2l" << std::endl;

	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");

	/////////////////////////////////
	// At least 2 fakeable leptons and no more than 2 tight leptons
	if (verbose_) {
		std::cout << "nFakeableLeptons = "<< fakeableLeps.size()
				  << "  nTightLeptons = " << tightLeps.size() << std::endl;
	}
	bool passLepNumber = fakeableLeps.size() >=2  and tightLeps.size() <= 2;
	if (passLepNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep num");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// lepton pt
	if (verbose_) {
		std::cout << "lep1 pdgid conept pt : " << fakeableLeps[0].pdgId() << " "
				  << fakeableLeps[0].conept() << " " << fakeableLeps[0].pt()
				  << std::endl;
		std::cout << "lep2 pdgid conept pt : " << fakeableLeps[1].pdgId() << " "
				  << fakeableLeps[1].conept() << " " << fakeableLeps[1].pt()
				  << std::endl;
	}
	float minpt = 25.;
	float minpt2 = abs(fakeableLeps[1].pdgId())==11 ? 15. : 10.;
	bool passLepPt =
		fakeableLeps[0].conept() >= minpt and fakeableLeps[1].conept() >= minpt2;
	if (passLepPt) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep pt");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton pT requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// muon, electron tight charge
	if (verbose_) {
		std::cout << "lep1 pdgid tigtcharge : " << fakeableLeps[0].pdgId() << " "
				  << fakeableLeps[0].passTightCharge() << std::endl;
		std::cout << "lep2 pdgid tigtcharge : " << fakeableLeps[1].pdgId() << " "
				  << fakeableLeps[1].passTightCharge() << std::endl;
	}
	bool passTightCharge =
		fakeableLeps[0].passTightCharge() and fakeableLeps[1].passTightCharge();
	if (passTightCharge) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tight charge");
	}
	else {
		if (verbose_) std::cout << "FAIL tight charge requirement" << std::endl;
		return false;
	}
	
	/////////////////////////////////
	// metLD > 0.2 (Pair of electrons only)
	if (abs(fakeableLeps[0].pdgId())==11 and abs(fakeableLeps[1].pdgId())==11) {
		if (verbose_) std::cout << "metLD : " << metLD << std::endl; 
		if (metLD > 0.2) {
			if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "metLD>0.2");
		}
		else {
			if (verbose_) std::cout << "FAIL metLD > 0.2" << std::endl;
			return false;
		}
	}

	/////////////////////////////////
	// Z mass veto: 91.2 +/- 10 (ee only)
	if ( pass_Zmass_veto(fakeableLeps, false, true) ) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Zmass veto");
	}
	else
		return false;

	/////////////////////////////////
	// At least 2 selected jets
	if (verbose_) std::cout << "nJets = " << njets << std::endl;
		
	bool passJetNumber = njets >= 2;
	if (passJetNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "jet num");
	}
	else {
		if (verbose_) std::cout << "FAIL jet number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// At least 2 loose btag or 1 medium btag
	if (verbose_) {
		std::cout << "nbtags loose : " << nbtags_loose << std::endl;
		std::cout << "nbtags medium : " << nbtags_medium << std::endl;
	}
	
	bool passBTagNumber = (nbtags_loose >= 2) or (nbtags_medium >= 1);
	if (passBTagNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "btag num");
	}
	else {
		if (verbose_) std::cout << "FAIL btag number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// Dilepton mass of any loose lepton pair > 12 GeV
	bool passMll = pass_pairMass_veto(looseLeps);
	if (passMll) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Mll>12GeV");
	}
	else
		return false;
	
	/////////////////////////////////
	if (verbose_) std::cout << "PASSED generic 2l event selection!" << std::endl;

	return true;
}

bool EventSelector::pass_2ltight_ss_selection(
    const std::vector<miniLepton>& tightLeps, int njets)
{
	int dummy = 0;
	return pass_2ltight_ss_selection(tightLeps, njets, dummy);
}

bool EventSelector::pass_2ltight_ss_selection(
    const std::vector<miniLepton>& tightLeps,
	int njets, int& ibin, TH1* h_cutflow)
{
	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");

	/////////////////////////////////
	// At least 2 tight leptons
	if (tightLeps.size() >= 2) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep num");
	}
	else
		return false;

	/////////////////////////////////
	// same sign
	if (tightLeps[0].charge() * tightLeps[1].charge() > 0) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "ss");
	}
	else
		return false;

	/////////////////////////////////
	// At least 2 selected jets
	bool passJetNumber = njets >= 2;
	if (passJetNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "jet num");
	}
	else
		return false;

	return true;
}

bool EventSelector::pass_2l1tau_inclusive_selection(
    const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniLepton>& tightLeps,
	const std::vector<miniTau>& fakeableTaus,
	int njets, int nbtags_loose, int nbtags_medium, float metLD,
    TH1* h_cutflow)
{
	if (verbose_) std::cout << "start event selection: 2l1tau" << std::endl;

	int ibin = 1;
	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");

	bool passes2lGenericSel =
		pass_2l_generic_selection(looseLeps, fakeableLeps, tightLeps,
								  njets, nbtags_loose, nbtags_medium, metLD,
								  ibin, h_cutflow);

	if (not passes2lGenericSel) {
		if (verbose_) std::cout << "FAIL generic 2l selection" << std::endl;
		return false;
	}

	/////////////////////////////////
	// At least 1 fakeable tau
	if (verbose_) std::cout << "nTaus = " << fakeableTaus.size() << std::endl;
	bool passTauNumber = fakeableTaus.size() > 0;
	if (passTauNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tau num");
	}
	else {
		if (verbose_)
			std::cout << "FAIL fakeable tau number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// At least 3 selected jets
	if (verbose_) std::cout << "nJets = " << njets << std::endl;
	
	bool passJetNumber = njets >= 3;
	if (passJetNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "njet>=3");
	}
	else {
		if (verbose_) std::cout << "FAIL jet number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 2l1tau inclusive selection!" << std::endl;

	return true;
}

bool EventSelector::pass_2lss1tau_tightLepID(const std::vector<miniLepton>& leptons)
{	
	assert(leptons.size()>1);

	if (looseselection_) return true;
	
	if (verbose_) {
		std::cout << "isTight lep1 lep2 : " << leptons[0].passTightSel()
				  << " " << leptons[1].passTightSel() << std::endl;
	}
	// SR: both leptons are tight
	bool passID = leptons[0].passTightSel() and leptons[1].passTightSel();
	return passID;
}

bool EventSelector::pass_2lss1tau_2lss(const std::vector<miniLepton>& leptons)
{
	assert(leptons.size()>1);
	if (verbose_) {
		std::cout << "charge lep1 lep2 : " << leptons[0].charge() << " "
				  << leptons[1].charge() << std::endl;
	}
	// SR: two leptons same sign
	bool pass2lss = leptons[0].charge() * leptons[1].charge() > 0;
	return pass2lss;
}

bool EventSelector::pass_2lss1tau_taucharge(const miniTau& tau, const miniLepton& lep)
{
	if (verbose_) {
		std::cout << "tau charge : " << tau.charge() << std::endl;
		std::cout << "lep charge : " << lep.charge() << std::endl;
	}
	// SR: opposite sign between tau and lepton
	bool passTauCharge = tau.charge() * lep.charge() < 0;
	return passTauCharge;
}

bool EventSelector::pass_2lss1tau_tauNumber(const std::vector<miniTau>& selectedTaus)
{
	// SR: at least one tau pass Loose MVA ID
	if (verbose_)
		std::cout << "number of taus : " << selectedTaus.size() << std::endl;
	
	if (selectedTaus.size() <= 0) {
		if (verbose_) std::cout << "FAIL tau number requirement" << std::endl;
		return false;
	}
	
	// At most one tau passing Medium MVA (avoid overlap with 2l2tau)
	int nvtighttau = 0;
	for (const auto & tau : selectedTaus) {
		if (not looseselection_) assert(tau.passMVAID("L"));
		if (tau.passMVAID("M")) nvtighttau++;
	}
	if ( nvtighttau > 1 ) {
		if (verbose_) std::cout << "FAIL: more than one Medium taus" << std::endl;
		return false;
	}

	return true;
}

bool EventSelector::pass_2lss1tau_SR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& selectedTaus)
{
	// event is assumed to already pass 2l1tau inclusive selection
	if (verbose_) std::cout << "start 2lss1tau signal region selection" << std::endl;

	///////////////////////////////
	// At least one tau pass medium MVA ID
	// and at most one tau passing VTight MVA (avoid overlap with 2l2tau)
	if ( not pass_2lss1tau_tauNumber(selectedTaus) )
		return false;
	
	///////////////////////////////
	// lepton ID
	if ( not pass_2lss1tau_tightLepID(fakeableLeps) ) {
		if (verbose_) std::cout << "FAIL tight lepton requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// same sign leptons
	if ( not pass_2lss1tau_2lss(fakeableLeps) ) {
		if (verbose_) std::cout << "FAIL lepton same sign requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// tau and leptons opposite sign
	assert(selectedTaus.size()>0);
	if ( not pass_2lss1tau_taucharge(selectedTaus[0], fakeableLeps[0]) ) {
		if (verbose_) std::cout << "FAIL tau charge requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// MC Matching
	if (isMC_) {
		assert(fakeableLeps.size()>1);
		bool passMCMatch = fakeableLeps[0].isGenMatched() and
			fakeableLeps[1].isGenMatched();
		if (not passMCMatch) {
			if (verbose_) {
				std::cout << "mcMatchType lep1 lep2 : "
						  << fakeableLeps[0].MCMatchType() << " "
						  << fakeableLeps[1].MCMatchType() << std::endl;
				std::cout << "FAIL MC Matching" << std::endl;
			}
			return false;
		}
	}
	
	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 2lss1tau signal region selection!" << std::endl;
	
	return true;
}

bool EventSelector::pass_2lss1tau_FakeAR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& selectedTaus)
{
	assert(not looseselection_);
	
	if (not pass_2lss1tau_tauNumber(selectedTaus))
		return false;
	
	assert(fakeableLeps.size()>1 and selectedTaus.size()>0);

	if (not pass_2lss1tau_2lss(fakeableLeps)) {
		if (verbose_) std::cout << "FAIL lepton same sign requirement" << std::endl;
		return false;
	}

	if (not pass_2lss1tau_taucharge(selectedTaus[0], fakeableLeps[0])) {
		if (verbose_) std::cout << "FAIL tau charge requirement" << std::endl;
		return false;
	}

	if (pass_2lss1tau_tightLepID(fakeableLeps)) {
		if (verbose_) std::cout << "FAIL lepton ID requirement" << std::endl;
		return false;
	}

	if (verbose_) std::cout << "PASSED 2lss1tau Fake AR selection!" << std::endl;
	return true;
}

bool EventSelector::pass_2lss1tau_FlipAR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& selectedTaus)
{
	assert(not looseselection_);
	
	if (not pass_2lss1tau_tauNumber(selectedTaus))
		return false;
	
	assert(fakeableLeps.size()>1 and selectedTaus.size()>0);

	bool lep1IsElectron = abs(fakeableLeps[0].pdgId())==11;
	bool lep2IsElectron = abs(fakeableLeps[1].pdgId())==11;

	if (verbose_ and not (lep1IsElectron or lep2IsElectron)) {
		std::cout << "FAIL: neither leptons are electrons" << std::endl;
	}
	
	// The lepton that has the same sign as tau has to be a electron
	// And the charge flip rate is only applied to this electron
	bool passTauCharge = false;
	if (lep1IsElectron and fakeableLeps[0].charge()== selectedTaus[0].charge())
		passTauCharge = true;
	if (lep2IsElectron and fakeableLeps[1].charge()== selectedTaus[0].charge())
		passTauCharge = true;

	if (verbose_ and not passTauCharge) {
		std::cout << "lep1 isElectron charge : " << lep1IsElectron << " "
				  << fakeableLeps[0].charge() << std::endl;
		std::cout << "lep2 isElectron charge : " << lep2IsElectron << " "
				  << fakeableLeps[1].charge() << std::endl;
		std::cout << "tau charge : " << selectedTaus[0].charge() << std::endl;
		std::cout << "FAIL tau charge requirement" << std::endl;
	}
	
	return ( (lep1IsElectron or lep2IsElectron) and
			 pass_2lss1tau_tightLepID(fakeableLeps) and
			 not pass_2lss1tau_2lss(fakeableLeps) and
			 passTauCharge );
}

// ttW control regions
bool EventSelector::pass_ttW_CR_selection(
    const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniLepton>& tightLeps,
	int njets, int nbtags_loose, int nbtags_medium, float metLD,
    TH1* h_cutflow)
{
	if (verbose_)
		std::cout << "start event selection: ttW control region" << std::endl;

	int ibin = 1;
	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");

	bool pass2lGenericSel =
		pass_2l_generic_selection(looseLeps, fakeableLeps, tightLeps,
								  njets, nbtags_loose, nbtags_medium, metLD,
								  ibin, h_cutflow);
	// need to loose tight lep number cut
	
	if (not pass2lGenericSel) {
		if (verbose_) std::cout << "FAIL generic 2l selection" << std::endl;
		return false;
	}

	// Extend Z mass veto to SFOS lepton pairs
	if (pass_Zmass_veto(fakeableLeps, true, false)) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Zmass veto (sfos)");
	}
	
	bool pass2ltightssSel =
		pass_2ltight_ss_selection(tightLeps, njets, ibin, h_cutflow);
	if (not pass2ltightssSel) {
		if (verbose_)
			std::cout << "FAIL 2l tight and same sign selection" << std::endl;
		return false;
	}

	// At least 2 jets passing medium WP b-tag
	if (nbtags_medium >= 2) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "btag num");
	}
	else {
		if (verbose_)
			std::cout << "FAIL to have at least 2 medium btags" << std::endl;
		return false;
	}

	if (verbose_) std::cout << "PASSED ttW control region selection!" << std::endl;
	return true;
}

/////////////////////////////////
// 3l1tau
/////////////////////////////////
bool EventSelector::pass_3l_generic_selection(
    const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	int njets, float metLD, int& ibin, TH1* h_cutflow)
{
	if (verbose_) std::cout << "start event selection: generic 3l" << std::endl;

	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");
	
	/////////////////////////////////
	// At least 3 fakeable leptons
	if (verbose_)
		std::cout << "nFakeableLeptons = " << fakeableLeps.size() << std::endl;
	bool passLepNumber = fakeableLeps.size() >= 3;
	if (passLepNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep num");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// lepton pt
	if (verbose_) {
		std::cout << "lep1 pdgid conept pt : " << fakeableLeps[0].pdgId() << " "
				  << fakeableLeps[0].conept() << " " << fakeableLeps[0].pt()
				  << std::endl;
		std::cout << "lep2 pdgid conept pt : " << fakeableLeps[1].pdgId() << " "
				  << fakeableLeps[1].conept() << " " << fakeableLeps[1].pt()
				  << std::endl;
		std::cout << "lep3 pdgid conept pt : " << fakeableLeps[2].pdgId() << " "
				  << fakeableLeps[2].conept() << " " << fakeableLeps[2].pt()
				  << std::endl;
	}
	bool passLepPt = fakeableLeps[0].conept() >= 20. and
		fakeableLeps[1].conept() >= 10. and fakeableLeps[2].conept() >= 10.;
	if (passLepPt) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep pt");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton pT requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// metLD 
	if ( pass_metLD_3l(metLD, fakeableLeps, njets) ) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "metLD cut");
	}
	else
		return false;

	/////////////////////////////////
	// At least 2 selected jets
	if (verbose_) std::cout << "nJets = " << njets << std::endl;
		
	bool passJetNumber = njets >= 2;
	if (passJetNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "jet num");
	}
	else {
		if (verbose_) std::cout << "FAIL jet number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// Dilepton mass of any loose lepton pair > 12 GeV
	bool passMll = pass_pairMass_veto(looseLeps);
	if (passMll) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Mll>12GeV");
	}
	else
		return false;

	/////////////////////////////////
	if (verbose_) std::cout << "PASSED generic 3l event selection!" << std::endl;

	return true;
}

bool EventSelector::pass_3l1tau_inclusive_selection(
    const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& fakeableTaus,
	int njets, int nbtags_loose, int nbtags_medium, float metLD,
    TH1* h_cutflow)
{
	if (verbose_) std::cout << "start event selcetion: 3l1tau" << std::endl;

	int ibin = 1;
	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");

	bool passes3lGenericSel =
		pass_3l_generic_selection(looseLeps, fakeableLeps, njets, metLD,
								  ibin, h_cutflow);

	if (not passes3lGenericSel) {
		if (verbose_) std::cout << "FAIL generic 3l selection" << std::endl;
		return false;
	}

	/////////////////////////////////
	// At least 1 fakeable tau
	if (verbose_) std::cout << "nTaus = " << fakeableTaus.size() << std::endl;
	bool passTauNumber = fakeableTaus.size() > 0;
	if (passTauNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tau num");
	}
	else {
		if (verbose_) std::cout << "FAIL tau number requirement" << std::endl;
		return false;
	}
	
	/////////////////////////////////
	// Z mass veto: 91.2 +/- 10 GeV (SFOS)
	if ( pass_Zmass_veto(fakeableLeps, true, false) ) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Zmass veto");
	}
	else
		return false;

	/////////////////////////////////
	// At least 2 loose btag or 1 medium btag
	if (verbose_) {
		std::cout << "nbtags loose : " << nbtags_loose << std::endl;
		std::cout << "nbtags medium : " << nbtags_medium << std::endl;
	}
	
	bool passBTagNumber = (nbtags_loose >= 2) or (nbtags_medium >= 1);
	if (passBTagNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "btag num");
	}
	else {
		if (verbose_) std::cout << "FAIL btag number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 3l1tau inclusive selection!" << std::endl;

	return true;
}

bool EventSelector::pass_3l1tau_tightID(const std::vector<miniLepton>& leptons)
{	
	assert(leptons.size()>2);

	if (looseselection_) return true;
	
	if (verbose_) {
		std::cout << "isTight lep1 lep2 lep3 : " << leptons[0].passTightSel()
				  << " " << leptons[1].passTightSel() << " "
				  << leptons[2].passTightSel() << std::endl;
	}
	// SR: all three leptons are tight
	return (leptons[0].passTightSel() and leptons[1].passTightSel() and
			leptons[2].passTightSel());
}

bool EventSelector::pass_3l1tau_charge(
    const std::vector<miniLepton>& leptons,
	const std::vector<miniTau>& taus)
{
	assert(leptons.size()>2 and taus.size()>0);

	if (verbose_) {
		std::cout << "charge lep1 lep2 lep3 : " << leptons[0].charge() << " "
				  << leptons[1].charge() << " " <<  leptons[2].charge() << std::endl;
		std::cout << "charge tau : " << taus[0].charge() << std::endl;
	}
	// SR: charge sum of the three leptons and tau is zero
	int chargesum = leptons[0].charge() + leptons[1].charge() + leptons[2].charge()
		+ taus[0].charge();
	return chargesum==0;
}

bool EventSelector::pass_3l1tau_tauNumber(const std::vector<miniTau>& selectedTaus)
{
	assert(not looseselection_);
	
	if (verbose_)
		std::cout << "number of taus : " << selectedTaus.size() << std::endl;
	
	for (const auto & tau : selectedTaus)
		if (not looseselection_) assert(tau.passMVAID("L"));

	return selectedTaus.size() > 0;
}

bool EventSelector::pass_3l1tau_SR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& selectedTaus)
{	
	// event is assumed to already pass 3l1tau inclusive selection
	if (verbose_) std::cout << "start 3l1tau signal region selection" << std::endl;

	///////////////////////////////
	// tau number
    if ( not pass_3l1tau_tauNumber(selectedTaus) ) {
		if (verbose_) std::cout << "FAIL tau number requirement" << std::endl;
		return false;
	}
	
	///////////////////////////////
	// lepton ID
	if ( not pass_3l1tau_tightID(fakeableLeps) ) {
		if (verbose_) std::cout << "FAIL tight lepton requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// charge sum
	if ( not pass_3l1tau_charge(fakeableLeps, selectedTaus) ) {
		if (verbose_) std::cout << "FAIL charge sum" << std::endl;
		return false;
	}

	///////////////////////////////
	// MC Matching
	if (isMC_) {
		assert(fakeableLeps.size()>2);
		bool passMCMatch = fakeableLeps[0].isGenMatched() and
			fakeableLeps[1].isGenMatched() and fakeableLeps[2].isGenMatched();
		if (not passMCMatch) {
			if (verbose_) {
				std::cout << "mcMatchType lep1 lep2 lep3 : "
						  << fakeableLeps[0].MCMatchType() << " "
						  << fakeableLeps[1].MCMatchType() << " "
						  << fakeableLeps[2].MCMatchType() << std::endl;
				std::cout << "FAIL MC Matching" << std::endl;
			}
			return false;
		}
	}
	
	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 3l1tau signal region selection!" << std::endl;
	
	return true;	
}

bool EventSelector::pass_3l1tau_FakeAR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& selectedTaus)
{
	assert(not looseselection_);
	
	if (not pass_3l1tau_tauNumber(selectedTaus)) {
		if (verbose_) std::cout << "FAIL tau number requirement" << std::endl;
		return false;
	}

	if (pass_3l1tau_tightID(fakeableLeps)) {
		if (verbose_) std::cout << "FAIL lepton ID requirement" << std::endl;
		return false;
	}

	if (not pass_3l1tau_charge(fakeableLeps, selectedTaus)) {
		if (verbose_) std::cout << "FAIL charge sum" << std::endl;
		return false;
	}

	if (verbose_) std::cout << "PASSED 2l2tau Fake AR selection!" << std::endl;
	return true;
}

bool EventSelector::pass_3l1tau_CR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& selectedTaus)
{
	assert(not looseselection_);
	
	if (not pass_3l1tau_tauNumber(selectedTaus))
		return false;
	
	return ( pass_3l1tau_tightID(fakeableLeps) and 
			 not pass_3l1tau_charge(fakeableLeps, selectedTaus));
}

bool EventSelector::pass_3l1tau_FakeARCR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& selectedTaus)
{
	assert(not looseselection_);
	
	if (not pass_3l1tau_tauNumber(selectedTaus))
		return false;
	
	return ( not pass_3l1tau_tightID(fakeableLeps) and 
			 not pass_3l1tau_charge(fakeableLeps, selectedTaus));
}

// ttZ control region
bool EventSelector::pass_ttZ_CR_selection(
    const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniLepton>& tightLeps,
    int njets, int nbtags_medium, float metLD,
	TH1* h_cutflow)
{
	if (verbose_)
		std::cout << "start event selection: ttZ control region" << std::endl;

	int ibin = 1;
	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");

	bool passes3lGenericSel =
		pass_3l_generic_selection(looseLeps, fakeableLeps, njets, metLD, ibin,
								  h_cutflow);
	if (not passes3lGenericSel) {
		if (verbose_) std::cout << "FAIL generic 3l selection" << std::endl;
		return false;
	}

	// At least 3 tight leptons
	if (verbose_) std::cout << "nTightLeptons = " << tightLeps.size() << std::endl;
	if (tightLeps.size() >= 3) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, ">=3 tight lep");
	}
	else {
		if (verbose_)
			std::cout << "FAIL to have at least 3 tight leptons" << std::endl;
		return false;
	}

	// tighter lepton pT cut
	assert(tightLeps.size()>2);
	bool passpt = tightLeps[0].conept()>25. and tightLeps[2].conept()>10. and
		( (tightLeps[1].conept()>15. and abs(tightLeps[1].pdgId())==11) or
		  (tightLeps[1].conept()>10. and abs(tightLeps[1].pdgId())==13) );
	if (passpt) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep pt");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton pt cuts" << std::endl;
		return false;
	}
	
	// At least 2 jets passing medium WP b-tag
	if (verbose_) std::cout << "n medium btags : " << nbtags_medium << std::endl;
	if (nbtags_medium >= 2) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "btag num");
	}
	else {
		if (verbose_)
			std::cout << "FAIL to have at least 2 medium btags" << std::endl;
		return false;
	}

	// At least 1 same flavour, opposite charge lepton pair passes the Z mass window cut: |mll - mZ| < 10 GeV 
	if ( not pass_Zmass_veto(tightLeps, true, false) ) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Zmass");
	}
	else
		return false;

	if (verbose_) std::cout << "PASSED ttZ control region selection!" << std::endl;

	return true;
}

/*
TODO
bool EventSelector::pass_3l_inclusive_CR_selection()
{

}

bool EventSelector::pass_3l_WZ_CR_selection()
{

}
*/

/////////////////////////////////
// 2l2tau
/////////////////////////////////
bool EventSelector::pass_2l2tau_inclusive_selection(
	const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& fakeableTaus,
	int njets, int nbtags_loose, int nbtags_medium, float metLD,
	TH1* h_cutflow)
{
	if (verbose_) std::cout << "start event selection: 2l2tau" << std::endl;

	int ibin = 1;

	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");

	/////////////////////////////////
	// At least 2 fakeable lepton
	if (verbose_) {
		std::cout << "nFakeableLeptons = " << fakeableLeps.size() << std::endl;
	}
	//
	bool passLepNumber = fakeableLeps.size() > 1;
	if (passLepNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep num");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// lepton pt
	if (verbose_) {
		std::cout << "lep1 pdgid conept pt : " << fakeableLeps[0].pdgId() << " "
				  << fakeableLeps[0].conept() << " " << fakeableLeps[0].pt()
				  << std::endl;
		std::cout << "lep2 pdgid conept pt : " << fakeableLeps[1].pdgId() << " "
				  << fakeableLeps[1].conept() << " " << fakeableLeps[1].pt()
				  << std::endl;
	}
    float minpt = 25.;
	float minpt2 = abs(fakeableLeps[1].pdgId())==11 ? 15. : 10.;
	bool passLepPt =
		fakeableLeps[0].conept() >= minpt and fakeableLeps[1].conept() >= minpt2;
	if (passLepPt) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep pt");
	}
	else {
		if (verbose_) std::cout << "FAIL lepton pT requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// At least 2 fakeable taus
	if (verbose_) std::cout << "nFakeableTaus = " << fakeableTaus.size() << std::endl;
	if (fakeableTaus.size()>1) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tau num");
	}
	else {
		if (verbose_) std::cout << "FAIL tau number requirement" << std::endl;
		return false;
	}
	
	/////////////////////////////////
	// Dilepton mass of any loose lepton pair > 12 GeV
	bool passMll = pass_pairMass_veto(looseLeps);
	if (passMll) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Mll>12GeV");
	}
	else
		return false;

	/////////////////////////////////
	// Z mass veto: 91.2 +/- 10 GeV (SFOS)
	if ( pass_Zmass_veto(fakeableLeps, true, false) ) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "Zmass veto");
	}
	else
		return false;

	/////////////////////////////////
	// metLD
	if ( pass_metLD_3l(metLD, fakeableLeps, njets) ) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "metLD cut");
	}
	else
		return false;

	/////////////////////////////////
	// At least 2 jets
	if (verbose_) std::cout << "nJets = " << njets << std::endl;
	bool passJetNumber = njets >= 2;
	if (passJetNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "jet num");
	}
	else {
		if (verbose_) std::cout << "FAIL jet number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// At least 2 loose btag or 1 medium btag
	if (verbose_) {
		std::cout << "nbtags loose : " << nbtags_loose << std::endl;
		std::cout << "nbtags medium : " << nbtags_medium << std::endl;
	}
	
	bool passBTagNumber = (nbtags_loose >= 2) or (nbtags_medium >= 1);
	if (passBTagNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "btag num");
	}
	else {
		if (verbose_) std::cout << "FAIL btag number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 2l2tau inclusive selection!" << std::endl;

	return true;
}

bool EventSelector::pass_2l2tau_tightID(const std::vector<miniLepton>& fakeableLeps,
										const std::vector<miniTau>& tightTaus)
{
	// SR: the 2 leading leptons are tight and >=2 tau pass tight selection
	assert(fakeableLeps.size()>1);
	if (verbose_) {
		std::cout << "lep1 istight: " << fakeableLeps[0].passTightSel() << std::endl;
		std::cout << "lep2 istight: " << fakeableLeps[1].passTightSel() << std::endl;
		std::cout << "number of tight taus: " << tightTaus.size() << std::endl;
	}

	if (tightTaus.size()<2) return false;

	if (not (fakeableLeps[0].passTightSel() and fakeableLeps[1].passTightSel()))
		return false;

	return true;
}

bool EventSelector::pass_2l2tau_charge(const std::vector<miniLepton>& leptons,
									   const std::vector<miniTau>& taus)
{
	assert(leptons.size()>1 and taus.size()>1);

	if (verbose_) {
		std::cout << "lep charges: " << leptons[0].charge() << " "
				  << leptons[1].charge() << std::endl;
		std::cout << "tau charges: " << taus[0].charge() << " " << taus[1].charge()
				  << std::endl;
	}
	
	// SR: sum of the charge == 0
	bool passCharge = (leptons[0].charge()+leptons[1].charge()+taus[0].charge()+
					   taus[1].charge())==0;
	return passCharge;
}

bool EventSelector::pass_2l2tau_SR_selection(
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& tightTaus)
{
	// event is assumed to already pass 2l2tau inclusive selection
	if (verbose_) std::cout << "start 2l2tau signal region selection" << std::endl;

	///////////////////////////////
	// lepton and tau ID
	if (not pass_2l2tau_tightID(fakeableLeps, tightTaus)) {
		if (verbose_) std::cout << "FAIL lepton and tau ID requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// charges
	if (not pass_2l2tau_charge(fakeableLeps, tightTaus)) {
		if (verbose_) std::cout << "FAIL charge sum requirement" << std::endl;
		return false;
	}

	// MC Matching
	if (isMC_) {
		assert(fakeableLeps.size()>1 and tightTaus.size()>1);
		bool passMCMatch =
			fakeableLeps[0].isGenMatched() and fakeableLeps[1].isGenMatched() and
			tightTaus[0].isGenMatched() and tightTaus[1].isGenMatched();
		if (not passMCMatch) {
			if (verbose_) {
				std::cout << "isGenMatched lep1 lep2 tau1 tau2: "
						  << fakeableLeps[0].isGenMatched() << " "
						  << fakeableLeps[1].isGenMatched() << " "
						  << tightTaus[0].isGenMatched() << " "
						  << tightTaus[1].isGenMatched() << std::endl;
				std::cout << "mcMatchType lep1 lep2 tau1 tau2: "
						  << fakeableLeps[0].MCMatchType() << " "
						  << fakeableLeps[1].MCMatchType() << " "
						  << tightTaus[0].MCMatchType() << " "
						  << tightTaus[1].MCMatchType() << std::endl;
				std::cout << "FAIL MC Matching" << std::endl;
			}
			return false;
		}	
	}

	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 2l2tau signal region selection!" << std::endl;

	return true;
}

bool EventSelector::pass_2l2tau_FakeAR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& tightTaus,
	const std::vector<miniTau>& fakeableTaus)
{
	assert(not looseselection_);
	if (pass_2l2tau_tightID(fakeableLeps, tightTaus)) {
		if (verbose_) std::cout << "FAIL lepton and tau ID requirement" << std::endl;
		return false;
	}
	
	if (not pass_2l2tau_charge(fakeableLeps, fakeableTaus)) {
		if (verbose_) std::cout << "FAIL charge sum requirement" << std::endl;
		return false;
	}

	if (verbose_) std::cout << "PASSED 2l2tau Fake AR selection!" << std::endl;
	return true;
}

bool EventSelector::pass_2l2tau_CR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& tightTaus)
{
	assert(not looseselection_);

	if (not pass_2l2tau_tightID(fakeableLeps, tightTaus)) return false;
	if (pass_2l2tau_charge(fakeableLeps, tightTaus)) return false;
	return true;

	// MC match?
}

bool EventSelector::pass_2l2tau_FakeARCR_selection(
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& fakeableTaus,
	const std::vector<miniTau>& tightTaus)
{
	assert(not looseselection_);

	if (pass_2l2tau_tightID(fakeableLeps, tightTaus)) return false;
	if (pass_2l2tau_charge(fakeableLeps, fakeableTaus)) return false;
	return true;
}

bool EventSelector::pass_pairMass_veto(const std::vector<miniLepton>& leps)
{
	if (leps.size()<2) return true;
	
	// veto two leptons with invariant mass  < 12 GeV
	for (auto it = leps.begin(); it != leps.end()-1; ++it) {
		for (auto it2 = it+1; it2 != leps.end(); ++it2) {
	if ( (it->p4() + it2->p4()).M() < 12. ) {
		if (verbose_) {
			std::cout << "FAIL Mll >= 12 GeV" << std::endl;
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

bool EventSelector::pass_Zmass_veto(const std::vector<miniLepton>& leps,
									bool opposite_sign, bool electron_only)
{
	// Zmass veto: 91.2 +/- 10	
	if (leps.size() < 2) return true;

	for (auto it = leps.begin(); it != leps.end()-1; ++it) {
		for (auto it2 = it+1; it2 != leps.end(); ++it2) {
			bool skip = false;
			
			if (electron_only) {
				if (not (abs(it->pdgId())==11 and abs(it2->pdgId())==11)) skip = true;
			}
			
			if (opposite_sign) {  // SFOS
				if (it->pdgId() + it2->pdgId() != 0) skip = true;
			}

			if (skip) continue;
			
			double invMass = (it->p4() + it2->p4()).M();
			if (invMass > (91.2 - 10.0) and invMass < (91.2 + 10.0)) {
				if (verbose_) {
					std::cout << "FAIL Z mass veto" << std::endl;
					std::cout << "lepton pair invariant mass : " << invMass
							  << std::endl;
				}
				return false;
			}
		}
	}

	return true;
}

bool EventSelector::pass_metLD_3l(float metLD, const std::vector<miniLepton>& leps,
								  int njets)
{
	// metLD cut
	assert(leps.size() > 1);

	if (verbose_)
		std::cout << "metLD : " << metLD << std::endl;
	
	if (njets >= 4)
		return true;
	else { // njets < 4
		bool sfos = false;
		for (auto it = leps.begin(); it != leps.end()-1; ++it) {
			for (auto it2 = it+1; it2 != leps.end(); ++it2) {
				if (it->pdgId() + it2->pdgId() == 0) sfos = true;
			}
		}
			
		double cut = sfos ? 0.3 : 0.2;
		if (metLD > cut) return true;
		
		if (verbose_)
			std::cout << "njets: " << njets << " sfos: " << sfos << std::endl;
	}

	if (verbose_) {
		std::cout << "FAIL metLD cut" << std::endl;
	}
	return false;
}


//////////////////////////
// Deprecated methods
//////////////////////////
bool EventSelector::pass_lepton_number(const std::vector<miniLepton>& lep_fakeable,
									   const std::vector<miniLepton>& lep_tight)
{
	if (verbose_) {
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
	else if (anaType_==Analyze_3l1tau) {
		// at least 3 fakeable leptons
		if (lep_fakeable.size() >= 3)
			return true;
	}
	else
		std::cout << "Analysis type not available!!" << std::endl;
	
	if (verbose_) {
		std::cout << "FAIL lepton number requirement" << std::endl;
	}
	return false;
	
}

bool EventSelector::pass_lepton_pt(const std::vector<miniLepton>& leps)
{
	if (verbose_) {
		for (const auto& lep : leps)
			std::cout << "lep pt id " << lep.pt()<<" "<< lep.pdgId() << std::endl;
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
		float minpt = abs(leps[0].pdgId())==11 ? 25. : 20.;

		if (leps[0].pt() > minpt)
			return true;
	}
	else if (anaType_==Analyze_3l1tau) {
		assert(leps.size() >= 3);
		// lepton pt
		if (leps[0].pt() > 20. and leps[1].pt() > 10. and leps[2].pt() > 10.)
			return true;
	}
	else
		std::cout << "Analysis type not available!!" << std::endl;

	if (verbose_) {
		std::cout << "FAIL lepton pT cut" << std::endl;
	}
	return false;
}

bool EventSelector::pass_lepton_charge(int lep0Charge, int lep1Charge)
{
	// 2lss1tau category only
	assert(anaType_==Analyze_2lss1tau);

	if (verbose_) {
		std::cout << "lep charges : " << lep0Charge << " " << lep1Charge
				  << std::endl;
	}
	
	bool samesign = lep0Charge * lep1Charge > 0;
	bool pass = (selType_==Control_2los1tau) ? (!samesign) : samesign;	
	if (pass)
		return true;

	if (verbose_) {
		std::cout << "FAIL lepton charge requirement" << std::endl;
	}
	return false;

}

bool EventSelector::pass_tight_charge(const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	// tight charge
	if (verbose_) {
		std::cout << "lep0 id tightCharge : " << leps[0].pdgId() << " "
				  << leps[0].passTightCharge()
				  << std::endl;
		std::cout << "lep1 id tightCharge : " << leps[1].pdgId() << " "
				  << leps[1].passTightCharge()
				  << std::endl;
	}
	
	if (leps[0].passTightCharge() and leps[1].passTightCharge())
		return true;

	if (verbose_) {
		std::cout << "FAIL tight charge" << std::endl;
	}
	return false;
}

bool EventSelector::pass_Zmass_veto(const std::vector<miniLepton>& leps)
{
	// Zmass veto: 91.2 +/- 10	
	if (leps.size() < 2) return true;

	for (auto it = leps.begin(); it != leps.end()-1; ++it) {
		for (auto it2 = it+1; it2 != leps.end(); ++it2) {
			bool skip = false;
			if (anaType_==Analyze_2lss1tau) { // ee only
				if (not (abs(it->pdgId())==11 and abs(it2->pdgId())==11) )
					skip = true;
			}
			else if (anaType_==Analyze_3l1tau) {  // SFOS
				if (it->pdgId() + it2->pdgId() != 0)
					skip = true;
			}

			if (skip) continue;
			
			double invMass = (it->p4() + it2->p4()).M();
			if (invMass > (91.2 - 10.0) and invMass < (91.2 + 10.0)) {
				if (verbose_) {
					std::cout << "FAIL Z mass veto" << std::endl;
					std::cout << "lepton pair invariant mass : " << invMass
							  << std::endl;
				}
				return false;
			}
		}
	}

	return true;
}

bool EventSelector::pass_metLD(float metLD, const std::vector<miniLepton>& leps,
							   int njets)
{
	// metLD cut
	assert(leps.size() >= 2);

	if (verbose_)
		std::cout << "metLD : " << metLD << std::endl;
	
	if (anaType_==Analyze_2lss1tau) {  //ee only
		if (not (abs(leps[0].pdgId())==11 and abs(leps[1].pdgId())==11) )
			return true;

		if (metLD > 0.2) return true;
	}
	else if (anaType_==Analyze_3l1tau) {
		if (njets >= 4)
			return true;
		else { // njets < 4
			bool sfos = false;
			for (auto it = leps.begin(); it != leps.end()-1; ++it) {
				for (auto it2 = it+1; it2 != leps.end(); ++it2) {
					if (it->pdgId() + it2->pdgId() == 0) sfos = true;
				}
			}
			
			double cut = sfos ? 0.3 : 0.2;
			if (metLD > cut) return true;

			if (verbose_)
				std::cout << "njets: " << njets << " sfos: " << sfos << std::endl;
		}
	}

	if (verbose_) {
		std::cout << "FAIL metLD cut" << std::endl;
	}
	return false;
}

/*
{
	// 2lss1tau category only for now
	assert(anaType_==Analyze_2lss1tau);

	// metLD cut
	// ee only
	assert(leps.size() >= 2);
	if (not (abs(leps[0].pdgId())==11 and abs(leps[1].pdgId())==11) ) 
		return true;

	if (verbose_)
		std::cout << "metLD : " << metLD << std::endl;
	
	if (metLD > 0.2)
		return true;

	if (verbose_) {
		std::cout << "FAIL metLD cut" << std::endl;
	}
	return false;
}*/

bool EventSelector::pass_lepton_ID(bool lep0IsTight, bool lep1IsTight,
								   bool lep2IsTight)
{
	if (verbose_) {
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

		if (passLepWP) return true;
	}
	else if (anaType_==Analyze_1l2tau) {
		// signal region: leading lepton is tight
		bool passLepWP = lep0IsTight;

		//if (selType_==Control_fake_1l2tau) passLepWP = not lep0IsTight;

		if (passLepWP) return true;
	}
	else if (anaType_==Analyze_3l1tau) {
		// signal region: the three leading leptons are tight
		bool passLepWP = lep0IsTight and lep1IsTight and lep2IsTight;

		if (selType_==Control_fake_3l1tau) passLepWP = not passLepWP;

		if (passLepWP) return true;
	}
	else
		std::cout << "Analysis type not available!!" << std::endl;	

	if (verbose_) {
		std::cout << "FAIL lepton WP requirement" << std::endl;
	}
	return false;
}

bool EventSelector::pass_tau_number(int ntaus)
{
	if (verbose_) {
		std::cout << "nTau : " << ntaus << std::endl;
	}

	if (anaType_==Analyze_2lss1tau or anaType_==Analyze_3l1tau) {
		// at least 1 selected tau
		if (ntaus >= 1) return true;
	}
	else if (anaType_==Analyze_1l2tau) {
		if (ntaus >= 2) return true;
	}
	else
		std::cout << "Analysis type not available!!" << std::endl;	
	
	if (verbose_) {
		std::cout << "FAIL tau number requirement" << std::endl;
	}
	return false;
}

bool EventSelector::pass_charge_sum(int tauCharge,
									const std::vector<miniLepton>& leps)
{
	// 3l1tau category only
	assert(anaType_==Analyze_3l1tau);
	assert(leps.size() >= 3);

	if (verbose_) {
		std::cout << "tau charge : " << tauCharge << std::endl;
		std::cout << "lepton charges : " << leps[0].charge() << " "
				  << leps[1].charge() << " " << leps[2].charge() << std::endl;
	}

	int chargesum = tauCharge + leps[0].charge() + leps[1].charge()
		+leps[2].charge();

	if (chargesum == 0) return true;

	if (verbose_) {
		std::cout << "FAIL charge sum requirement" << std::endl;
	}
	return false;
}

bool EventSelector::pass_tau_charge(int tauCharge,
									const std::vector<miniLepton>& leps)
{
	// 2lss1tau category only
	assert(anaType_==Analyze_2lss1tau);

	assert(leps.size() >= 2);

	if (verbose_) {
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

	if (verbose_) {
		std::cout << "FAIL tau charge requirement" << std::endl;
	}
	return false;
}

bool EventSelector::pass_taupair_charge(int tau0charge, int tau1charge)
{
	// 1l2tau category only
	assert(anaType_==Analyze_1l2tau);
	assert(selType_==Signal_1l2tau or selType_==Control_fake_1l2tau or
		   selType_==Loose_1l2tau);
	
	if (verbose_) {
		std::cout << "tau charge : " << tau0charge << " " << tau1charge
				  << std::endl;
	}
		
	bool opposite = tau0charge * tau1charge < 0;

	if (opposite) return true;
	
	if (verbose_) {
		std::cout << "FAIL tau pair charge requirement" << std::endl;
	}
	return false;
	
}

bool EventSelector::pass_tau_ID(int ntau_tight)
{
	// for 1l2tau only
	assert(anaType_==Analyze_1l2tau);

	if (verbose_) {
		std::cout << "number of tight taus : " << ntau_tight << std::endl;
	}

	bool passTauWP = ntau_tight >= 2;

	if (passTauWP) return true;

	if (verbose_) {
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
	if (anaType_==Analyze_3l1tau) {
		if (njets >= 2) return true;
	}
	else {
		if (njets >= 3) return true;
	}

	if (verbose_) {
		std::cout << "FAIL number of jets requirement" << std::endl;
		std::cout << "njets : " << njets << std::endl;
	}
	return false;
}

bool EventSelector::pass_btag_number(int nbtags_loose, int nbtags_medium)
{
	if (verbose_) {
		std::cout << "loose btags : " << nbtags_loose << std::endl;
		std::cout << "medium btags : " << nbtags_medium << std::endl;
	}
	
	if (nbtags_loose >= 2 or nbtags_medium >= 1)
		return true;

	if (verbose_) {
		std::cout << "FAIL number of btags requirement" << std::endl;
	}
	return false;
}

bool EventSelector::pass_lep_mc_match(const miniLepton& lep)
{
	// for MC only

	if (verbose_)
		std::cout << "MC match type : " << lep.MCMatchType() << std::endl;

	int isPrompt = abs(lep.pdgId())==11 ? 1 : 2;
	int isPromptTauDecay = abs(lep.pdgId())==11 ? 3 : 4;

	if (lep.MCMatchType()==isPrompt or lep.MCMatchType()==isPromptTauDecay)
		return true;

	if (verbose_) {
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

	if (verbose_) {
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

	if (verbose_) {
		std::cout << "MC match type : " << mtype << std::endl;
	}
	
	if (mtype==1 or mtype==2 or mtype==3 or mtype==4 or mtype==5)
		return true;
	
	if (verbose_) {
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

	if (verbose_) {
		std::cout << "FAIL tau MC match" << std::endl;
	}
	return false;
}
