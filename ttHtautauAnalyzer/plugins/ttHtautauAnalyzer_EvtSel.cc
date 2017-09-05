#ifndef ttHtautauAnalyzer_EvtSel_cc
#define ttHtautauAnalyzer_EvtSel_cc

#include "ttHTauTauAnalysis/ttHtautauAnalyzer/plugins/ttHtautauAnalyzer.h"

bool ttHtautauAnalyzer::pass_event_sel_2lss1tau (
    const std::vector<miniLepton>& lep_loose,
    const std::vector<miniLepton>& lep_fakeable,
	const std::vector<miniLepton>& lep_tight,
	const std::vector<pat::Tau>& taus,
    const int njets, const int nbtags_loose, const int nbtags_medium,
    const float metLD)
{
	if (debug_) std::cout << "start event selection: 2lss1tau" << std::endl;

	int ibin = 1;
	if (doCutflow_) fill_CutFlow(ibin++,"total");
	
	//////////////////////////
	// at least 2 fakeable leptons and no more than 2 tight leptons
	if (evt_selector_->pass_lepton_number(lep_fakeable, lep_tight)) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep num");
	}
	else
		return false;

	//////////////////////////
	// at least 1 selected tau
	// "byMediumIsolationMVArun2v1DBdR03oldDMwLT"
	if (evt_selector_->pass_tau_number(taus.size())) {
		if (doCutflow_) fill_CutFlow(ibin++,"tau num");
	}
	else
		return false;

	//////////////////////////
	// lepton pt
	if (evt_selector_->pass_lepton_pt(lep_fakeable)) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep pt");
	}
	else
		return false;

	//////////////////////////
	// tight charge
	if (evt_selector_->pass_tight_charge(lep_fakeable)) {
		if (doCutflow_) fill_CutFlow(ibin++,"tight charge");
	}
	else
		return false;

	//////////////////////////
	// veto two loose leptons with invariant mass < 12 GeV
	if (evt_selector_->pass_pairMass_veto(lep_loose)) {
		if (doCutflow_) fill_CutFlow(ibin++,"Mll<12GeV");
	}
	else
		return false;
	
	//////////////////////////
	// Z mass veto: 91.2 +/- 10 (ee only)
	if (evt_selector_->pass_Zmass_veto(lep_fakeable)) {
		if (doCutflow_) fill_CutFlow(ibin++,"Zmass veto");
	}
	else
		return false;

	//////////////////////////
	// metLD cut (ee only)
	if (evt_selector_->pass_metLD(metLD, lep_fakeable)) {
		if (doCutflow_) fill_CutFlow(ibin++,"metLD");
	}
	else
		return false;

	//////////////////////////
	// number of jets
	if (evt_selector_->pass_jet_number(njets)) {
		if (doCutflow_) fill_CutFlow(ibin++,"jet num");
	}
	else
		return false;
	
	//////////////////////////
	// number of btags
	if (evt_selector_->pass_btag_number(nbtags_loose, nbtags_medium)) {
		fill_CutFlow(ibin++,"btag num");
	}
	else
		return false;
	
	//////////////////////////
	// lepton charge
	assert(lep_fakeable.size() >= 2);
	if (evt_selector_->pass_lepton_charge(lep_fakeable[0].charge(),
										  lep_fakeable[1].charge())
		) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep charge");
	}
	else
		return false;

	//////////////////////////
	// tau charge
	/*
	assert(taus.size() >= 1);
	if (evt_selector_->pass_tau_charge(taus[0].charge(), lep_fakeable)) {
		if (doCutflow_) fill_CutFlow(ibin++,"tau charge");
	}
	else
		return false;
	*/
	// for signal region, opposite sign between tau and either lepton
	// for control region, the lepton that has the same sign as tau has to be a electron, and the charge flip rate is only applied to this electron
	// To save computing time, additional requirement on tau charge applied after ntuple production for signal or control region
	// save the selection flag in the ntuple

	//////////////////////////
	// lepton WP
	if (evt_selector_->pass_lepton_ID(lep_fakeable[0].passTightSel(),
									  lep_fakeable[1].passTightSel())
		) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep WP");
	}
	else
		return false;

	//////////////////////////
	// MC truth matching
	if (not isdata_) {
		/*
		if (evt_selector_->pass_lep_mc_match(lep_fakeable)) {
			if (doCutflow_) fill_CutFlow(ibin++,"lep MC");
		}
		else
			return false;
		*/
		/*
		assert(taus.size() >= 1);
		if (evt_selector_->pass_tau_mc_match(taus[0])) {
			if (doCutflow_) fill_CutFlow(ibin++,"tau MC");
		}
		else
			return false;
		*/
	}

	//////////////////////////
	if (debug_) std::cout << "PASSED event selection!" << std::endl;

	return true;
}

bool ttHtautauAnalyzer::pass_event_sel_1l2tau (
	const std::vector<miniLepton>& lep_loose,
    const std::vector<miniLepton>& lep_fakeable,
	const std::vector<miniLepton>& lep_tight,
	const std::vector<pat::Tau>& taus_fakeable,
	const std::vector<pat::Tau>& taus,
    const int njets, const int nbtags_loose, const int nbtags_medium)
{
	if (debug_) std::cout << "start event selection: 1l2tau" << std::endl;

	int ibin = 1;
	if (doCutflow_) fill_CutFlow(ibin++,"total");

	//////////////////////////
	// at least 1 fakeable lepton and no more than 1 tight lepton
	if (evt_selector_->pass_lepton_number(lep_fakeable, lep_tight)) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep num");
	}
	else
		return false;

	//////////////////////////
	// lepton pt
	if (evt_selector_->pass_lepton_pt(lep_fakeable)) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep pt");
	}
	else
		return false;

	//////////////////////////
	// lepton eta cut (to match trigger)
	assert(lep_fakeable.size()>0);
	if (abs(lep_fakeable[0].eta())<2.1) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep eta");
	}
	else
		return false;

	//////////////////////////
	// veto two loose leptons with invariant mass < 12 GeV
	if (evt_selector_->pass_pairMass_veto(lep_loose)) {
		if (doCutflow_) fill_CutFlow(ibin++,"Mll<12GeV");
	}
	else
		return false;
	
	//////////////////////////
	// at least 2 fakeable taus
	if (evt_selector_->pass_tau_number(taus_fakeable.size())) {
		if (doCutflow_) fill_CutFlow(ibin++,"tau num");
	}
	else
		return false;

	//////////////////////////
	// tau pt
	assert(taus_fakeable.size()>0);
	if (taus_fakeable[0].pt()>30) {
		if (doCutflow_) fill_CutFlow(ibin++,"tau pt");
	}
	else
		return false;

	//////////////////////////
	// lepton and taus WP
	assert(lep_fakeable.size()>0);
	if (evt_selector_->pass_lep_tau_ID(lep_fakeable[0].passTightSel(),
									   taus.size())) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep&tau WP");
	}
	else
		return false;

	//////////////////////////
	// tau charge
	/*
	bool pass_taupair_charge = false;
	if (selType_==Signal_1l2tau) {
		assert(taus.size()>1);
		pass_taupair_charge =
			evt_selector_->pass_taupair_charge(taus[0].charge(),taus[1].charge());
	}
	else if (selType_==Control_fake_1l2tau) {
		assert(taus_fakeable.size()>1);
		pass_taupair_charge =
			evt_selector_->pass_taupair_charge(taus_fakeable[0].charge(),
											   taus_fakeable[1].charge());
	}
	else {
		std::cout << "Selection type not available!!" << std::endl;
		return false;
	}

	if (pass_taupair_charge) {
		if (doCutflow_) fill_CutFlow(ibin++,"tau charge");
	}
	else
		return false;
	*/
	// Requirement on tau charge is applied after ntuple production for signal or control region
	
	//////////////////////////
	// number of jets
	if (evt_selector_->pass_jet_number(njets)) {
		if (doCutflow_) fill_CutFlow(ibin++, "jet num");
	}
	else
		return false;

	//////////////////////////
	// number of btags
	if (evt_selector_->pass_btag_number(nbtags_loose, nbtags_medium)) {
		fill_CutFlow(ibin++, "btag num");
	}
	else
		return false;

	//////////////////////////
	// MC truth matching
	if (not isdata_ and selType_==Signal_1l2tau) {
		if (evt_selector_->pass_lep_mc_match(lep_fakeable[0])) {
			fill_CutFlow(ibin++, "lep MC");
		}
		else
			return false;

		assert(taus.size()>=2);
		if (evt_selector_->pass_tau_mc_match(taus)) {
			fill_CutFlow(ibin++, "tau MC");
		}
		else
			return false;
	}

	//////////////////////////
	if (debug_) std::cout << "PASSED event selection!" << std::endl;

	return true;
}											  

bool ttHtautauAnalyzer::pass_event_sel_3l1tau (
    const std::vector<miniLepton>& lep_loose,
    const std::vector<miniLepton>& lep_fakeable,
	const std::vector<miniLepton>& lep_tight,
	const std::vector<pat::Tau>& taus,
    const int njets, const int nbtags_loose, const int nbtags_medium,
    const float metLD)
{
	if (debug_) std::cout << "start event selection: 3l1tau" << std::endl;

	int ibin = 1;
	if (doCutflow_) fill_CutFlow(ibin++,"total");

	//////////////////////////
	// at least 3 fakeable leptons
	if (evt_selector_->pass_lepton_number(lep_fakeable, lep_tight)) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep num");
	}
	else
		return false;

	//////////////////////////
	// lepton pt
	if (evt_selector_->pass_lepton_pt(lep_fakeable)) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep pt");
	}
	else
		return false;

	//////////////////////////
	// veto two loose leptons with invariant mass < 12 GeV
	if (evt_selector_->pass_pairMass_veto(lep_loose)) {
		if (doCutflow_) fill_CutFlow(ibin++,"Mll<12GeV");
	}
	else
		return false;

	//////////////////////////
	// Z mass veto: 91.2 +/- 10 (SFOS)
	if (evt_selector_->pass_Zmass_veto(lep_fakeable)) {
		if (doCutflow_) fill_CutFlow(ibin++,"Zmass veto");
	}
	else
		return false;

	//////////////////////////
	// lepton WP
	if (evt_selector_->pass_lepton_ID(lep_fakeable[0].passTightSel(),
									  lep_fakeable[1].passTightSel(),
									  lep_fakeable[2].passTightSel())) {
		if (doCutflow_) fill_CutFlow(ibin++,"lep WP");
	}
	else
		return false;
	
	//////////////////////////
	// at least 1 selected tau
	// "byMediumIsolationMVArun2v1DBdR03oldDMwLT"
	if (evt_selector_->pass_tau_number(taus.size())) {
		if (doCutflow_) fill_CutFlow(ibin++,"tau num");
	}
	else
		return false;

	//////////////////////////
	// number of jets
	if (evt_selector_->pass_jet_number(njets)) {
		if (doCutflow_) fill_CutFlow(ibin++,"jet num");
	}
	else
		return false;

	//////////////////////////
	// number of btags
	if (evt_selector_->pass_btag_number(nbtags_loose, nbtags_medium)) {
		fill_CutFlow(ibin++,"btag num");
	}
	else
		return false;

	//////////////////////////
	// metLD cut
	if (evt_selector_->pass_metLD(metLD, lep_fakeable, njets)) {
		if (doCutflow_) fill_CutFlow(ibin++,"metLD");
	}
	else
		return false;
	
	//////////////////////////
	// Charge sum
	// do not apply cut here
	// save the charge sum in ntuple
	/*
	if (evt_selector_->pass_charge_sum(taus[0].charge(), lep_fakeable)) {
		if (doCutflow_) fill_CutFlow(ibin++,"charge sum");
	}
	else
		return false;
	*/
	
	//////////////////////////
	// MC truth matching
	if (not isdata_) {
		// do not apply cut here
		// save the MC truth matching info in ntuple
	}

	//////////////////////////
	if (debug_) std::cout << "PASSED event selection!" << std::endl;

	return true;
}

void ttHtautauAnalyzer::fill_CutFlow(int ibin, const char* name)
{
	assert(h_CutFlow_);

	if (not firstpass_)
		h_CutFlow_->GetXaxis()->SetBinLabel(ibin, name);

	if (ibin > h_CutFlow_->GetNbinsX()) {
		std::cout << "WARNING : ibin " << ibin
				  << " exceeds number of bins of histogram" << std::endl;
	}
	
	h_CutFlow_->Fill(ibin);

	return;
}

#endif
