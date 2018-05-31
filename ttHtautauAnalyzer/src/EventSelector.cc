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
	else if (anatype == Analyze_2lss1tau or anatype == Analyze_2lss) {
		passhlt =
			trighelper->pass_single_lep_triggers(triggerBits) or
			trighelper->pass_dilep_triggers(triggerBits);
	}
	else if (anatype == Analyze_3l1tau or anatype == Analyze_3l) {
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
	else if (anatype == Analyze_2lss1tau or anatype == Analyze_2lss) {
		pass = (pass_e and nElectron>0) or (pass_2e and nElectron>1) or
			(pass_m and nMuon>0) or (pass_2m and nMuon>1) or
			(pass_em and nElectron>0 and nMuon>0);
	}
	else if (anatype == Analyze_3l1tau or anatype == Analyze_3l) {
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

bool EventSelector::pass_hlt_and_filters(Analysis_types anatype,
										 TriggerHelper* const trighelper,
										 unsigned int triggerBits,
										 int nElectron, int nMuon,
										 unsigned int filterBits, bool isdata)
{
	if (not trighelper->pass_filters(filterBits, isdata)) {
		if (verbose_) {
			std::cout << "FAIL MET filters" << std::endl;
			std::cout << "filter bits : " << filterBits << std::endl;
		}
		return false;
	}
	
	if (not pass_hlt_paths(anatype, trighelper, triggerBits))
		return false;

	if (not pass_hlt_match(anatype, trighelper, triggerBits, nElectron, nMuon))
		return false;

	if (verbose_) std::cout << "PASSED HLT and MET filters!" << std::endl;
		
	return true;
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

	// nlep >=2 or (nlep==1 and ntau>=2)
	bool passLepTauNum = (fakeableLeps.size() == 1 and fakeableTaus.size() >= 2)
		or fakeableLeps.size() >= 2;

	if (passLepTauNum) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "nlep>1||(nlep==1&&ntau>1)");
	}
	else {
		if (verbose_)
			std::cout << "FAIL lepton and tau number requirement" << std::endl;
		return false;
	}

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
	bool pass = false;
	if (verbose_) std::cout << std::endl;

	if (anatype == Analyze_1l2tau) {
		if (not pass_1l2tau_inclusive_selection(looseLeps, fakeableLeps, tightLeps,
												fakeableTaus, njets, nbtags_loose,
												nbtags_medium, h_cutflow) )
			return false;

		pass = pass_1l2tau_selection(seltype, fakeableLeps, fakeableTaus);
	}
	else if (anatype == Analyze_2lss1tau) {
		if (not pass_2l1tau_inclusive_selection(looseLeps, fakeableLeps,
												tightLeps, fakeableTaus, njets,
												nbtags_loose, nbtags_medium, metLD,
												h_cutflow) )
			return false;

		pass = pass_2lss1tau_selection(seltype, fakeableLeps, selectedTaus);
	}
	else if (anatype == Analyze_3l1tau) {
		if (not pass_3l1tau_inclusive_selection(looseLeps, fakeableLeps, fakeableTaus,
												njets, nbtags_loose, nbtags_medium,
												metLD, h_cutflow) )
			return false;

		pass = pass_3l1tau_selection(seltype, fakeableLeps, selectedTaus);
	}
	else if (anatype == Analyze_2l2tau) {
		if (not pass_2l2tau_inclusive_selection(looseLeps, fakeableLeps, fakeableTaus,
												njets, nbtags_loose, nbtags_medium,
												metLD, h_cutflow) )
			return false;

		pass = pass_2l2tau_selection(seltype, fakeableLeps, fakeableTaus);
	}
	else if (anatype == Analyze_2lss) {
		if (seltype==Control_ttW or seltype==Control_FakeAR_ttW or
			seltype==Control_FlipAR_ttW) {
			pass = pass_ttW_CR_selection(seltype, looseLeps, fakeableLeps, tightLeps,
										 selectedTaus, njets, nbtags_loose,
										 nbtags_medium, metLD, h_cutflow);
		}
		else {  // not yet
			assert(0);
		}
	}
	else if (anatype == Analyze_3l) {
		if (seltype==Control_ttZ or seltype==Control_FakeAR_ttZ) {
			pass = pass_ttZ_CR_selection(seltype, looseLeps, fakeableLeps, tightLeps,
										 selectedTaus, njets, nbtags_loose,
										 nbtags_medium, metLD, h_cutflow);
		}
		else if (seltype==Control_WZ or seltype==Control_FakeAR_WZ) { // coming soon?
			assert(0);
		}
		else {
			assert(0);
		}
	}
	else {
		std::cerr << "Analysis type is not supported. Return false." << std::endl;
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
	
	bool passTauPt = fakeableTaus[0].pt() >= 30. and fakeableTaus[1].pt() >= 20.; 
	if (passTauPt) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tau pt");
	}
	else {
		if (verbose_) std::cout << "FAIL tau pT requirement" << std::endl;
		return false;
	}

	/*
	/////////////////////////////////
	// tau eta (to match trigger)
	if (verbose_) {
		std::cout << "Two leading tau eta : " << fakeableTaus[0].eta() << " "
				  << fakeableTaus[1].eta() << std::endl;
	}

	bool passTauEta =
		fabs(fakeableTaus[0].eta())<2.1 and fabs(fakeableTaus[1].eta())<2.1;
	if (passTauEta) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tau eta");
	}
	else {
		if (verbose_) std::cout << "FAIL tau eta requirement" << std::endl;
		return false;
	}
	*/
	
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
/*
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
*/
// check leading two fakeable tau ID instead of counting number of tight taus
bool EventSelector::pass_1l2tau_tightID(const std::vector<miniLepton>& fakeableLeps,
										const std::vector<miniTau>& fakeableTaus)
{
	// SR: tight lepton; both taus pass tight selection (VTight MVA)
	if (verbose_) {
		std::cout << "Istight lep : " << fakeableLeps[0].passTightSel() << std::endl;
		std::cout << "Istight tau1 tau2 : " << fakeableTaus[0].passTightSel() << " "
				  << fakeableTaus[1].passTightSel() << std::endl;
	}

	assert(fakeableLeps.size()>0 and fakeableTaus.size()>1);

	bool passTightID = fakeableLeps[0].passTightSel() and
		fakeableTaus[0].passTightSel() and fakeableTaus[1].passTightSel();

	return passTightID;
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

bool EventSelector::pass_1l2tau_selection(Selection_types seltype,
										  const std::vector<miniLepton>& fakeableLeps,
										  const std::vector<miniTau>& fakeableTaus)
{
	// event is assumed to already pass 1l2tau inclusive selection
	if (verbose_) std::cout << "start 1l2tau selection" << std::endl;

	///////////////////////////////
	// lepton and tau ID
	bool pass_ID = pass_1l2tau_tightID(fakeableLeps, fakeableTaus);
	// invert the cut for fake application region
	if (seltype==Application_Fake_1l2tau or seltype==Control_FakeAR_1l2tau)
		pass_ID = not pass_ID;

	if (not pass_ID) {
		if (verbose_) std::cout << "FAIL lepton and tau ID requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// tau charges
	bool pass_charge = pass_1l2tau_charge(fakeableTaus);
	// invert the cut for control region
	if (seltype==Control_1l2tau or seltype==Control_FakeAR_1l2tau)
		pass_charge = not pass_charge;

	if (not pass_charge) {
		if (verbose_) std::cout << "FAIL tau charge requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// MC Matching
	if (isMC_ and (seltype==Signal_1l2tau or seltype==Control_1l2tau)) {// CR too?
		bool passMCMatch = pass_MCMatch(fakeableLeps, 1, fakeableTaus, 2);
		if (not passMCMatch) {
			if (verbose_) std::cout << "FAIL MC Matching" << std::endl;
			return false;
		}
	}

	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 1l2tau selection!" << std::endl;

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
	if (pass_lepton_pT_2l(fakeableLeps)) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep pt");
	}
	else
		return false;

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
	if ( pass_Zmass_veto(looseLeps, false, true) ) {
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

bool EventSelector::pass_2l_tightLepID(const std::vector<miniLepton>& leptons)
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

bool EventSelector::pass_2l_2lSS(const std::vector<miniLepton>& leptons)
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

bool EventSelector::pass_2l1tau_leptauOS(const miniTau& tau, const miniLepton& lep)
{
	if (verbose_) {
		std::cout << "tau charge : " << tau.charge() << std::endl;
		std::cout << "lep charge : " << lep.charge() << std::endl;
	}
	// SR: opposite sign between tau and lepton
	bool passTauCharge = tau.charge() * lep.charge() < 0;
	return passTauCharge;
}

bool EventSelector::pass_2l1tau_tauNumber(const std::vector<miniTau>& selectedTaus)
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

bool EventSelector::pass_2lss1tau_selection(Selection_types seltype,
    const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniTau>& selectedTaus)
{
	// event is assumed to already pass 2l1tau inclusive selection
	if (verbose_) std::cout << "start 2lss1tau signal region selection" << std::endl;

	///////////////////////////////
	// At least one tau pass Loose MVA ID
	// and at most one tau passing Medium MVA (avoid overlap with 2l2tau)
	if ( not pass_2l1tau_tauNumber(selectedTaus) ) return false;

	///////////////////////////////
	// lepton ID
	bool pass_lepID = pass_2l_tightLepID(fakeableLeps);
	// invert the cut for fake application region
	if (seltype==Application_Fake_2lss1tau or seltype==Control_FakeAR_2lss1tau)
		pass_lepID = not pass_lepID;

	if (not pass_lepID) {
		if (verbose_) std::cout << "FAIL lepton ID requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// lepton charge
	bool pass_lepCharge = pass_2l_2lSS(fakeableLeps);
	// invert the cut for charge flip application region
	if (seltype==Application_Flip_2lss1tau or seltype==Control_FlipAR_2lss1tau) {
		pass_lepCharge = not pass_lepCharge;
	}
	
	if (not pass_lepCharge) {
		if (verbose_) std::cout << "FAIL lepton charge requirement" << std::endl;
		return false;
	}

	assert(fakeableLeps.size()>1);
	bool lep1IsElectron = abs(fakeableLeps[0].pdgId())==11;
	bool lep2IsElectron = abs(fakeableLeps[1].pdgId())==11;
	if (seltype==Application_Flip_2lss1tau or seltype==Control_FlipAR_2lss1tau) {
		// at least one electron
		if ( not (lep1IsElectron or lep2IsElectron) ) {
			if (verbose_)
				std::cout << "FAIL: neither leptons are electrons" << std::endl;
			return false;
		}
	}
	
	///////////////////////////////
	// tau charge
	bool pass_tauCharge = false;
	assert(selectedTaus.size()>0);

	if (verbose_) {
		std::cout << "lep1 isElectron charge : " << lep1IsElectron << " "
				  << fakeableLeps[0].charge() << std::endl;
		std::cout << "lep2 isElectron charge : " << lep2IsElectron << " "
				  << fakeableLeps[1].charge() << std::endl;
		std::cout << "tau charge : " << selectedTaus[0].charge() << std::endl;
	}
	
	if (seltype==Signal_2lss1tau or seltype==Application_Fake_2lss1tau) {
		pass_tauCharge = pass_2l1tau_leptauOS(selectedTaus[0], fakeableLeps[0]);
	}
	else if (seltype==Control_2lss1tau or seltype==Control_FakeAR_2lss1tau) {
		pass_tauCharge = not pass_2l1tau_leptauOS(selectedTaus[0], fakeableLeps[0]);
	}
	else if (seltype==Application_Flip_2lss1tau) {
		// For charge flip region, the lepton that has the same sign as tau has to be
		// a electron, and the charge flip rate is only applied to this electron
		if (lep1IsElectron and fakeableLeps[0].charge()==selectedTaus[0].charge())
			pass_tauCharge = true;
		if (lep2IsElectron and fakeableLeps[1].charge()==selectedTaus[0].charge())
			pass_tauCharge = true;
	}
	else if (seltype==Control_FlipAR_2lss1tau) {
		if (lep1IsElectron and fakeableLeps[0].charge()!=selectedTaus[0].charge())
			pass_tauCharge = true;
		if (lep2IsElectron and fakeableLeps[1].charge()!=selectedTaus[0].charge())
			pass_tauCharge = true;
		// Need to be careful to pick the right electron to apply charge flip rate
	}

	if (not pass_tauCharge) {
		if (verbose_) std::cout << "FAIL tau charge requirement" << std::endl;
		return false;
	}
	
	///////////////////////////////
	// MC Matching
	if (isMC_ and (seltype==Signal_2lss1tau or seltype==Control_2lss1tau)) {
		bool passMCMatch = pass_MCMatch(fakeableLeps, 2, selectedTaus, 0);
		if (not passMCMatch) {
			if (verbose_) std::cout << "FAIL MC Matching" << std::endl;
			return false;
		}
	}

	///////////////////////////////
	if (verbose_) std::cout << "PASSED 2lss1tau selection!" << std::endl;

	return true;
}

// 2lss signal region
//bool EventSelector::pass_2lss_selection() {}

// ttW control region
bool EventSelector::pass_ttW_CR_selection(Selection_types seltype,
    const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniLepton>& tightLeps,
	const std::vector<miniTau>& tightTaus,
	int njets, int nbtags_loose, int nbtags_medium, float metLD,
	TH1* h_cutflow)
{
	if (verbose_) std::cout << "start ttW control region selection" << std::endl;

	int ibin = 1;
	if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "total");

	bool passes2lGenericSel =
		pass_2l_generic_selection(looseLeps, fakeableLeps, tightLeps,
								  njets, nbtags_loose, nbtags_medium, metLD,
								  ibin, h_cutflow);

	if (not passes2lGenericSel) {
		if (verbose_) std::cout << "FAIL generic 2l selection" << std::endl;
		return false;
	}

	/////////////////////////////////
	// raise pT cut to 15 GeV for sub-leading lepton
	// for all 2lss categories
	assert(fakeableLeps.size()>1);
	if (verbose_) {
		std::cout << "lep2 pdgid conept pt : " << fakeableLeps[1].pdgId() << " "
				  << fakeableLeps[1].conept() << " " << fakeableLeps[1].pt()
				  << std::endl;;
	}
	bool passTightLepPt = fakeableLeps[1].conept() >= 15.;
	if (not passTightLepPt) {
		if (verbose_) std::cout << "FAIL lepton pT cut" << std::endl;
		return false;
	}
	
	/////////////////////////////////
	// No tau passing "byLooseIsolationMVArun2v1DBdR03oldDMwLT"
	if (verbose_) std::cout << "ntightTaus = " << tightTaus.size() << std::endl;
	bool passTauNumber = (tightTaus.size()==0);
	if (passTauNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tau num");
	}
	else {
		if (verbose_) std::cout << "FAIL tight tau number requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// jet number
	if (verbose_) std::cout << "nJets = " << njets << std::endl;
	
	// Exactly 3 selected jets
	bool passJetNumber = (njets==3);
	if (passJetNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "jet num");
	}
	else {
		if (verbose_) std::cout << "FAIL jet number requirement" << std::endl;
		return false;
	}

	// At least 4 selected jets for 2lss selection
	
	/////////////////////////////////
	// lepton ID
	bool pass_lepID = pass_2l_tightLepID(fakeableLeps);
	// invert the cut for fake application region
	if (seltype==Control_FakeAR_ttW)   // Application_Fake_2lss
		pass_lepID = not pass_lepID;
	
	if (not pass_lepID) {
		if (verbose_) std::cout << "FAIL lepton ID requirement" << std::endl;
		return false;
	}

	/////////////////////////////////
	// lepton charge
	bool pass_lepCharge = pass_2l_2lSS(fakeableLeps);
	// invert the cut for charge flip region
	if (seltype==Control_FlipAR_ttW)   // Application_Flip
		pass_lepCharge = not pass_lepCharge;

	if (not pass_lepCharge) {
		if (verbose_) std::cout << "FAIL lepton charge requirement" << std::endl;
		return false;
	}

	assert(fakeableLeps.size()>1);
    bool lep1IsElectron = abs(fakeableLeps[0].pdgId())==11;
	bool lep2IsElectron = abs(fakeableLeps[1].pdgId())==11;
	if (seltype==Control_FlipAR_ttW) {
		// at least one electron
		if ( not (lep1IsElectron or lep2IsElectron) ) {
			if (verbose_)
				std::cout << "FAIL: neither leptons are electrons" << std::endl;
			return false;
		}
	}

	///////////////////////////////
	// MC Matching
	if (isMC_ and (seltype==Control_ttW) ) {  // Signal_2lss
		bool passMCMatch = pass_MCMatch(fakeableLeps, 2, tightTaus, 0);
		if (not passMCMatch) {
			if (verbose_) std::cout << "FAIL MC Matching" << std::endl;
			return false;
		}
	}

	/*
	////////////////
	
	// At least 2 jets passing medium WP b-tag
	if (nbtags_medium >= 2) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "btag num");
	}
	else {
		if (verbose_)
			std::cout << "FAIL to have at least 2 medium btags" << std::endl;
		return false;
	}
	*/
	
	/////////////////////////////////
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
	if ( pass_Zmass_veto(looseLeps, true, false) ) {
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

bool EventSelector::pass_3l_tightLepID(const std::vector<miniLepton>& leptons)
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

bool EventSelector::pass_3l_charge(const std::vector<miniLepton>& leptons)
{
	assert(leptons.size()>2);

	if (verbose_) {
		std::cout << "charge lep1 lep2 lep3 : " << leptons[0].charge() << " "
				  << leptons[1].charge() << " " <<  leptons[2].charge() << std::endl;
	}
	// SR: charge sum of the three leptons == +/- 1
	int chargesum = leptons[0].charge() + leptons[1].charge() + leptons[2].charge();
	return abs(chargesum)==1;
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

bool EventSelector::pass_3l1tau_selection(Selection_types seltype,
										  const std::vector<miniLepton>& fakeableLeps,
										  const std::vector<miniTau>& selectedTaus)
{
	// event is assumed to already pass 3l1tau inclusive selection
	if (verbose_) std::cout << "start 3l1tau selection" << std::endl;

	///////////////////////////////
	// tau number
	if ( not pass_3l1tau_tauNumber(selectedTaus) ) {
		if (verbose_) std::cout << "FAIL tau number requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// lepton ID
	bool pass_lepID = pass_3l_tightLepID(fakeableLeps);
	// invert the cut for fake application region
	if (seltype==Application_Fake_3l1tau or seltype==Control_FakeAR_3l1tau)
		pass_lepID = not pass_lepID;

	if (not pass_lepID) {
		if (verbose_) std::cout << "FAIL lepton ID requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// charge sum
	bool pass_charge = pass_3l1tau_charge(fakeableLeps, selectedTaus);
	// invert the cut for control region
	if (seltype==Control_3l1tau or seltype==Control_FakeAR_3l1tau)
		pass_charge = not pass_charge;

	if (not pass_charge) {
		if (verbose_) std::cout << "FAIL charge sum requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// MC Matching
	if (isMC_ and (seltype==Signal_3l1tau or seltype==Control_3l1tau)) { // CR too?
		bool passMCMatch = pass_MCMatch(fakeableLeps, 3, selectedTaus, 0);
		if (not passMCMatch) {
			if (verbose_) std::cout << "FAIL MC Matching" << std::endl;
			return false;
		}
	}

	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 3l1tau selection!" << std::endl;
	
	return true;
}

// ttZ control region
bool EventSelector::pass_ttZ_CR_selection(Selection_types seltype,
    const std::vector<miniLepton>& looseLeps,
	const std::vector<miniLepton>& fakeableLeps,
	const std::vector<miniLepton>& tightLeps,
	const std::vector<miniTau>& tightTaus,
    int njets, int nbtags_loose, int nbtags_medium, float metLD,
	TH1* h_cutflow)
{
	if (verbose_) std::cout << "start ttZ control region selection" << std::endl;

	int ibin = 1;
	if (h_cutflow and ibin==1) fill_cutflow(h_cutflow, ibin++, "total");

	bool passes3lGenericSel =
		pass_3l_generic_selection(looseLeps, fakeableLeps, njets, metLD, ibin,
								  h_cutflow);
	if (not passes3lGenericSel) {
		if (verbose_) std::cout << "FAIL generic 3l selection" << std::endl;
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
	
	/* only in 3l
	///////////////////////////////
	// FIXME: Necessary?
	// No tau passing "byLooseIsolationMVArun2v1DBdR03oldDMwLT"
	if (verbose_) std::cout << "ntightTaus = " << tightTaus.size() << std::endl;
	bool passTauNumber = (tightTaus.size()==0);
	if (passTauNumber) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "tau num");
	}
	else {
		if (verbose_) std::cout << "FAIL tight tau number requirement" << std::endl;
		return false;
	}
	///////////////////////////////
	*/
	
	///////////////////////////////
	// lepton ID
	bool pass_lepID = pass_3l_tightLepID(fakeableLeps);
	// invert the cut for fake application region
	if (seltype==Control_FakeAR_ttZ)
		pass_lepID = not pass_lepID;

	if (not pass_lepID) {
		if (verbose_) std::cout << "FAIL lepton ID requirement" << std::endl;
		return false;
	}
	
	///////////////////////////////
	// lepton charge
	bool pass_charge = pass_3l_charge(fakeableLeps);
	if (not pass_charge) {
		if (verbose_) std::cout << "FAIL charge sum" << std::endl;
		return false;
	}

	// tighter lepton pT cut
	bool passLepPt_tighter = fakeableLeps[1].conept() >= 15.;
	if (not passLepPt_tighter) {
		if (verbose_) std::cout << "FAIL lepton pT cut" << std::endl;
		return false;
	}

	/* only in 3l
	///////////////////////////////
	// No more than 3 tight leptons
	if (tightLeps.size() > 3) {
		if (verbose_)
			std::cout << "FAIL due to more than 3 tight leptons" << std::endl;
		return false;
	}
	*/
	
	///////////////////////////////
	// MC Matching
	if (isMC_ and (seltype==Control_ttZ)) {
		bool passMCMatch = pass_MCMatch(fakeableLeps, 3, tightTaus, 0);
		if (not passMCMatch) {
			if (verbose_) std::cout << "FAIL MC Matching" << std::endl;
			return false;
		}
	}

	/* only in 3l
	///////////////////////////////
	// veto H->ZZ*->4l
	if (not pass_HZZ4l_veto(looseLeps)) {
		if (verbose_) std::cout << "FAIL H->ZZ*->4l veto" << std::endl;
		return false;
	}
	*/

	///////////////////////////////
	// reverted Z mass veto: 91.2 +/- 10 GeV (SFOS)
	if (pass_Zmass_veto(looseLeps, true, false)) {
		if (verbose_) std::cout << "FAIL reverted Z mass veto" << std::endl;
		return false;
	}

	/////////////////////////////////
	if (verbose_) std::cout << "PASSED ttZ control region selection!" << std::endl;

	return true;
}

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
	if (pass_lepton_pT_2l(fakeableLeps)) {
		if (h_cutflow) fill_cutflow(h_cutflow, ibin++, "lep pt");
	}
	else 
		return false;

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
	if ( pass_Zmass_veto(looseLeps, true, false) ) {
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
/*
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
*/
bool EventSelector::pass_2l2tau_tightID(const std::vector<miniLepton>& fakeableLeps,
										const std::vector<miniTau>& fakeableTaus)
{
	// SR: 2 leading leptons are tight and 2 leading taus are tight
	assert(fakeableLeps.size()>1 and fakeableTaus.size()>1);
	if (verbose_) {
		std::cout << "Istight lep1 lep2 : " << fakeableLeps[0].passTightSel() << " "
				  << fakeableLeps[1].passTightSel() << std::endl;
		std::cout << "Istight tau1 tau2 : " << fakeableTaus[0].passTightSel() << " "
				  << fakeableTaus[1].passTightSel() << std::endl;
	}

	bool passTightID =
		fakeableLeps[0].passTightSel() and fakeableLeps[1].passTightSel() and
		fakeableTaus[0].passTightSel() and fakeableTaus[1].passTightSel();

	return passTightID;
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

bool EventSelector::pass_2l2tau_selection(Selection_types seltype,
										  const std::vector<miniLepton>& fakeableLeps,
										  const std::vector<miniTau>& fakeableTaus)
{
	// event is assumed to already pass 2l2tau inclusive selection
	if (verbose_) std::cout << "start 2l2tau selection" << std::endl;

	///////////////////////////////
	// lepton and tau ID
	bool pass_ID = pass_2l2tau_tightID(fakeableLeps, fakeableTaus);
	// invert the cut for fake application region
	if (seltype==Application_Fake_2l2tau or seltype==Control_FakeAR_2l2tau)
		pass_ID = not pass_ID;

	if (not pass_ID) {
		if (verbose_) std::cout << "FAIL lepton and tau ID requirement" << std::endl;
		return false;
	}

	///////////////////////////////
	// charges
	bool pass_charge = pass_2l2tau_charge(fakeableLeps, fakeableTaus);
	// invert the cut for control region
	if (seltype==Control_2l2tau or seltype==Control_FakeAR_2l2tau)
		pass_charge = not pass_charge;
	
	if (not pass_charge) {
		if (verbose_) std::cout << "FAIL charge sum requirement" << std::endl;
		return false;
	}

	// MC Matching
	if (isMC_ and (seltype==Signal_2l2tau or seltype==Control_2l2tau)) {
		bool passMCMatch = pass_MCMatch(fakeableLeps, 2, fakeableTaus, 2);
		if (not passMCMatch) {
			if (verbose_) std::cout << "FAIL MC Matching" << std::endl;
			return false;
		}
	}

	/////////////////////////////////
	if (verbose_) std::cout << "PASSED 2l2tau selection!" << std::endl;

	return true;
}


bool EventSelector::pass_MCMatch(const std::vector<miniLepton>& leptons, int nLep,
								 const std::vector<miniTau>& taus, int nTau)
{
	assert((int)leptons.size() >= nLep);
	assert((int)taus.size() >= nTau);

	int ilep = 0, itau = 0;
	bool passMCMatch = true;
	
	for (const auto & lep : leptons) {
		if (ilep >= nLep) break;
		passMCMatch = passMCMatch and lep.isGenMatched();
		ilep++;
		
		if (verbose_) {
			std::cout << "lep" << ilep << " isGenMatched mcMatchType : "
					  << lep.isGenMatched() << " " << lep.MCMatchType() << std::endl;
		}
	}
	
	for (const auto & tau : taus) {
		if (itau >= nTau) break;
		passMCMatch = passMCMatch and tau.isGenMatched();
		itau++;

		if (verbose_) {
			std::cout << "tau" << itau << " isGenMatched mcMatchType : "
					  << tau.isGenMatched() << " " << tau.MCMatchType() << std::endl;
		}
	}

	return passMCMatch;
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
				    std::cout << "pt eta phi mass : " << it->pt() << " " << it->eta()
							  << " " << it->phi() << " " << it->mass() << std::endl;
					std::cout << "pt eta phi mass : " << it2->pt() << " " << it2->eta()
							  << " " << it2->phi() << " " << it2->mass() << std::endl;
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

bool EventSelector::pass_lepton_pT_2l(const std::vector<miniLepton>& leptons)
{
	if (leptons.size()<2) {
		if (verbose_) std::cout << "Not enough leptons!!" << std::endl;
		return false;
	}

	if (verbose_) {
		std::cout << "lep1 pdgid conept pt : " << leptons[0].pdgId() << " "
				  << leptons[0].conept() << " " << leptons[0].pt() << std::endl;
		std::cout << "lep2 pdgid conept pt : " << leptons[1].pdgId() << " "
				  << leptons[1].conept() << " " << leptons[1].pt() << std::endl;
	}
    float minpt = 25.;
	float minpt2 = abs(leptons[1].pdgId())==11 ? 15. : 10.;
	bool passLepPt =
		leptons[0].conept() >= minpt and leptons[1].conept() >= minpt2;
	if (passLepPt)
		return true;
	else {
		if (verbose_) std::cout << "FAIL lepton pT requirement" << std::endl;
		return false;
	}
}

bool EventSelector::pass_HZZ4l_veto(const std::vector<miniLepton>& leptons)
{
	// Veto events containing two SFOS pairs of preselected leptons, with mass m_4l < 140 GeV
	if (leptons.size()<4) return true;

	for (auto l1 = leptons.begin(); l1 != leptons.end()-3; ++l1) {
		for (auto l2 = l1+1; l2 != leptons.end()-2; ++l2) {
			for (auto l3=l2+1; l3 != leptons.end()-1; ++l3) {
				for (auto l4=l3+1; l4 != leptons.end(); ++l4) {
					// check if two SFOS pairs
					bool IsTwoSFOS = (l1->pdgId()+l2->pdgId()+l3->pdgId()+l4->pdgId() ==0);
					float m4l = (l1->p4()+l2->p4()+l3->p4()+l4->p4()).M();
					if (IsTwoSFOS and m4l<140.) return false;
				}
			}
		}
	}

	return true;
}
