#ifndef eventNtuple_cc
#define eventNtuple_cc

#include "../interface/eventNtuple.h"

std::vector<miniLepton> eventNtuple::buildLeptons(bool loose)
{
	std::vector<miniLepton> leptons;

	for (unsigned int l = 0; l<mu_pt->size(); ++l) {
		// if not loose selection require to at least pass fakeable id
		if (!loose and !(mu_isfakeablesel->at(l))) continue;
		TLorentzVector mu;
		mu.SetPtEtaPhiE(mu_pt->at(l),mu_eta->at(l),mu_phi->at(l),mu_E->at(l));
		miniLepton lep(mu, mu_conept->at(l), -13*mu_charge->at(l),
					   mu_charge->at(l), true, mu_isfakeablesel->at(l),
					   mu_ismvasel->at(l));
		leptons.push_back(lep);
	}

	for (unsigned int l = 0; l<ele_pt->size(); ++l) {
		// if not loose selection require to at least pass fakeable id
		if (!loose and !(ele_isfakeablesel->at(l))) continue;
		TLorentzVector ele;
		ele.SetPtEtaPhiE(ele_pt->at(l),ele_eta->at(l),ele_phi->at(l),ele_E->at(l));
		miniLepton lep(ele, ele_conept->at(l), -11*ele_charge->at(l),
					   ele_charge->at(l), true, ele_isfakeablesel->at(l),
					   ele_ismvasel->at(l));
		leptons.push_back(lep);
	}
	// sort by conept
	sort(leptons.begin(), leptons.end(), [] (miniLepton l1, miniLepton l2)
		 {return l1.conept()>l2.conept();} );

	return leptons;
}

std::vector<TLorentzVector> eventNtuple::buildFourVectorLeps(bool loose)
{
	auto leptons = buildLeptons(loose);
	
	std::vector<TLorentzVector> lepsP4;
	for (const auto & lep : leptons)
		lepsP4.push_back(lep.p4());

	return lepsP4;
}

std::vector<miniTau> eventNtuple::buildTaus(bool loose)
{
	std::vector<miniTau> taus;
		
	for (unsigned int t = 0; t < tau_pt->size(); ++t) {
		// require tight tau if not loose selection
		if (!loose and !(tau_idSelection->at(t))) continue;
		TLorentzVector tauP4;
		tauP4.SetPtEtaPhiE(tau_pt->at(t),tau_eta->at(t),tau_phi->at(t),tau_E->at(t));
		miniTau tau(tauP4,tau_charge->at(t),tau_decayMode->at(t),
					tau_idPreselection->at(t), tau_idSelection->at(t));
		if (tau_mcMatchType->size()>0)
			tau.set_MCMatchType(tau_mcMatchType->at(t));

		// Tau decay products
		// charged hadron
		std::vector<TLorentzVector> chP4;		
		for (unsigned int ch = 0; ch<(tau_signalChargedHadrCands_pt->at(t)).size();++ch) {
			TLorentzVector chargedhadr;
			chargedhadr.SetPtEtaPhiE((tau_signalChargedHadrCands_pt->at(t))[ch],
									 (tau_signalChargedHadrCands_eta->at(t))[ch],
									 (tau_signalChargedHadrCands_phi->at(t))[ch],
									 (tau_signalChargedHadrCands_E->at(t))[ch]);
			chP4.push_back(chargedhadr);
		}			
		tau.set_signalChargedHadrCands(chP4);

		// gamma (pi0)
		std::vector<TLorentzVector> gmP4;
		for (unsigned int g=0; g<(tau_signalGammaCands_pt->at(t)).size(); ++g) {
			TLorentzVector gamma;
			gamma.SetPtEtaPhiE((tau_signalGammaCands_pt->at(t))[g],
							   (tau_signalGammaCands_eta->at(t))[g],
							   (tau_signalGammaCands_phi->at(t))[g],
							   (tau_signalGammaCands_E->at(t))[g]);
			gmP4.push_back(gamma);
		}
		tau.set_signalGammaCands(gmP4);

		// neutral hadron
		std::vector<TLorentzVector> nhP4;
		for (unsigned int nh=0; nh<(tau_signalNeutrHadrCands_pt->at(t)).size(); ++nh) {
			TLorentzVector neutralhadr;
			neutralhadr.SetPtEtaPhiE((tau_signalNeutrHadrCands_pt->at(t))[nh],
									 (tau_signalNeutrHadrCands_eta->at(t))[nh],
									 (tau_signalNeutrHadrCands_phi->at(t))[nh],
									 (tau_signalNeutrHadrCands_E->at(t))[nh]);
		    nhP4.push_back(neutralhadr);
		}
		tau.set_signalNeutrHadrCands(nhP4);

		taus.push_back(tau);
	}
	// should be already sorted by pt
	return taus;
}

std::vector<TLorentzVector> eventNtuple::buildFourVectorTaus(std::vector<int>& decaymode, bool loose)
{
	
	std::vector<TLorentzVector> tausP4;

	decaymode.clear();
	
	for (unsigned int t = 0; t < tau_pt->size(); ++t) {
		// require tight tau if not loose selection
		if (!loose and !(tau_idSelection->at(t))) continue;

		TLorentzVector tau;
		tau.SetPtEtaPhiE(tau_pt->at(t),tau_eta->at(t),tau_phi->at(t),tau_E->at(t));
		tausP4.push_back(tau);
		decaymode.push_back(tau_decayMode->at(t));
	}
	// should be already sorted by pt
	return tausP4;
}

std::vector<TLorentzVector> eventNtuple::buildFourVectorTaus(bool loose)
{
	std::vector<int> dummy;
	return buildFourVectorTaus(dummy, loose);
}

std::vector<TLorentzVector> eventNtuple::buildFourVectorTauDaugsCharged(bool loose)
{
	std::vector<TLorentzVector> tauDaugsChargedP4;

	for (unsigned int t = 0; t < tau_pt->size(); ++t) {
		// require tight tau if not loose selection
		if (!loose and !(tau_idSelection->at(t))) continue;

		TLorentzVector chargedDaugsP4;
		
		// loop over tau charged daughters and add up four vectors
		for (unsigned int ch=0; ch<(tau_signalChargedHadrCands_pt->at(t)).size(); ++ch) {
			TLorentzVector chargedhadr;
			chargedhadr.SetPtEtaPhiE((tau_signalChargedHadrCands_pt->at(t))[ch],
									 (tau_signalChargedHadrCands_eta->at(t))[ch],
									 (tau_signalChargedHadrCands_phi->at(t))[ch],
									 (tau_signalChargedHadrCands_E->at(t))[ch]);
			chargedDaugsP4 += chargedhadr;
		}

		tauDaugsChargedP4.push_back(chargedDaugsP4);
	}

	return tauDaugsChargedP4;
}

std::vector<TLorentzVector> eventNtuple::buildFourVectorTauDaugsNeutral(bool loose)
{
	std::vector<TLorentzVector> tauDaugsNeutralP4;

	for (unsigned int t = 0; t < tau_pt->size(); ++t) {
		// require tight tau if not loose selection
		if (!loose and !(tau_idSelection->at(t))) continue;

		TLorentzVector neutralDaugsP4;

		// loop over tau neutral daughters and add up four vectors
		for (unsigned int g=0; g<(tau_signalGammaCands_pt->at(t)).size(); ++g) {
			TLorentzVector gamma;
			gamma.SetPtEtaPhiE((tau_signalGammaCands_pt->at(t))[g],
							   (tau_signalGammaCands_eta->at(t))[g],
							   (tau_signalGammaCands_phi->at(t))[g],
							   (tau_signalGammaCands_E->at(t))[g]);
			neutralDaugsP4 += gamma;
		}

		for (unsigned int nh=0; nh<(tau_signalNeutrHadrCands_pt->at(t)).size(); ++nh) {
			TLorentzVector neutralhadr;
			neutralhadr.SetPtEtaPhiE((tau_signalNeutrHadrCands_pt->at(t))[nh],
									 (tau_signalNeutrHadrCands_eta->at(t))[nh],
									 (tau_signalNeutrHadrCands_phi->at(t))[nh],
									 (tau_signalNeutrHadrCands_E->at(t))[nh]);
			neutralDaugsP4 += neutralhadr;
		}
		

		tauDaugsNeutralP4.push_back(neutralDaugsP4);
	}

	return tauDaugsNeutralP4;
}

std::vector<TLorentzVector> eventNtuple::buildFourVectorJets(bool loose)
{
	std::vector<TLorentzVector> jetsP4;
	
	for (unsigned int j = 0; j < jet_pt->size(); ++j) {
		TLorentzVector jet;
		jet.SetPtEtaPhiE(jet_pt->at(j),jet_eta->at(j),jet_phi->at(j),jet_E->at(j));
		jetsP4.push_back(jet);
	}
	// should be already sorted by pt
	return jetsP4;
}

std::vector<TLorentzVector> eventNtuple::buildFourVectorBtagJets(bool loose)
{
	std::vector<TLorentzVector> btagsP4;
	btagsP4.reserve(2);

	int icsv0 = -1; int icsv1 = -1;
	float csv0 = -233.; float csv1 = -233;
	
	for (unsigned int j = 0; j < jet_pt->size(); ++j) {
		if (jet_csv->at(j) > csv0) {
			csv1 = csv0;
			icsv1 = icsv0;
			csv0 = jet_csv->at(j);
			icsv0 = j;
		}
		else if (jet_csv->at(j) > csv1) {
			csv1 = jet_csv->at(j);
			icsv1 = j;
		}
	}
	
	TLorentzVector bj0, bj1;
	bj0.SetPtEtaPhiE(jet_pt->at(icsv0),jet_eta->at(icsv0),jet_phi->at(icsv0),
					 jet_E->at(icsv0));
	assert(jet_csv->at(icsv0)==csv0);
	bj1.SetPtEtaPhiE(jet_pt->at(icsv1),jet_eta->at(icsv1),jet_phi->at(icsv1),
					 jet_E->at(icsv1));
	assert(jet_csv->at(icsv1)==csv1);
	btagsP4.push_back(bj0);
	btagsP4.push_back(bj1);

	return btagsP4;
}

TLorentzVector eventNtuple::buildFourVectorMET()
{
	TLorentzVector met;
	met.SetPtEtaPhiM(PFMET, 0, PFMETphi, 0);
	return met;
}

void eventNtuple::initialize()
{
	// event variables
	//std::cout << "eventNtuple::initialize(): event ids" << std::endl;
	run = 0;
	ls = 0;
    nEvent = 0;

	//std::cout << "eventNtuple::initialize(): weights" << std::endl;
	event_weight = -9999.;
	PU_weight = -9999.;
	MC_weight = -9999.;
	bTagSF_weight = -9999.;
	leptonSF_weight = -9999.;
	tauSF_weight = -9999.;
	triggerSF_weight = -9999.;
	FR_weight = -9999.;

	/////////////////////////////
	// systematics
	MC_weight_scale_muF0p5 = -9999.;
	MC_weight_scale_muF2 = -9999.;
	MC_weight_scale_muR0p5 = -9999.;
	MC_weight_scale_muR2 = -9999.;
	btagSF_weight_LFUp = -9999.;
	btagSF_weight_LFDown = -9999.;
	btagSF_weight_HFUp = -9999.;
	btagSF_weight_HFDown = -9999.;
	btagSF_weight_HFStats1Up = -9999.;
	btagSF_weight_HFStats1Down = -9999.;
	btagSF_weight_HFStats2Up = -9999.;
	btagSF_weight_HFStats2Down = -9999.;
	btagSF_weight_LFStats1Up = -9999.;
	btagSF_weight_LFStats1Down = -9999.;
	btagSF_weight_LFStats2Up = -9999.;
	btagSF_weight_LFStats2Down = -9999.;
	btagSF_weight_cErr1Up = -9999.;
	btagSF_weight_cErr1Down = -9999.;
	btagSF_weight_cErr2Up = -9999.;
	btagSF_weight_cErr2Down = -9999.;
	/////////////////////////////

	HiggsDecayType = -9999;  // Higgs decay product pdgId

    lepCategory = -9999;    // 0: mumu; 1: ee; 2:emu
    btagCategory = -9999;   // 0: loose; 1: medium (>=2 medium btags)

	// pileup
	npuTrue = -9999.;
	npuInTime = -9999.;

	pvx = -9999.;
	pvy = -9999.;
	pvz = -9999.;
	
	// triggers
	// triggerBits_single_e = 0;
	// triggerBits_single_mu = 0;
	// triggerBits_double_e = 0;
	// triggerBits_double_mu = 0;
	// triggerBits_elemu = 0;
	triggerBits = 0;
	
	filterBits = 0;

	nBadMuons = -9999;

	// event selection flag
	passTauCharge = -9999;
	isGenMatchedLep = -9999;
	isGenMatchedTau  = -9999;
	
	// event level MVA variables
	//MVA_2lss_ttV;
	//MVA_2lss_ttbar;
	//MT_met_lep0;
	//mindr_lep0_jet;
	//mindr_lep1_jet;
	//lep0_conept;
	//lep1_conept;
	//avg_dr_jet;

	//ibin;  // bin index in 1D BDT shape template

    n_presel_mu = -9999;
    n_fakeable_mu = -9999;
    n_mvasel_mu = -9999;
    n_presel_ele = -9999;
	n_fakeable_ele = -9999;
	n_mvasel_ele = -9999;
	n_presel_tau = -9999;
	n_tau = -9999;
	n_presel_jet = -9999;
	n_jet = -9999;
	n_btag_medium = -9999;
	n_btag_loose = -9999;

	// muons
	//std::cout << "eventNtuple::initialize(): muons" << std::endl;
	mu_pt->clear();
	mu_conept->clear();
	mu_eta->clear();
	mu_phi->clear();
	mu_E->clear();
	mu_charge->clear();
	mu_dxy->clear();
	mu_dz->clear();
	mu_isfakeablesel->clear();
	mu_ismvasel->clear();
	mu_jetNDauChargedMVASel->clear();
	mu_miniRelIso->clear();
	mu_miniIsoCharged->clear();
	mu_miniIsoNeutral->clear();
	mu_jetPtRel->clear();
	mu_jetPtRatio->clear();
	mu_jetCSV->clear();
	mu_sip3D->clear();
	mu_segmentCompatibility->clear();
	mu_leptonMVA->clear();
	mu_mediumID->clear();
	mu_dpt_div_pt->clear();
	mu_mcMatchType->clear();
	mu_isPFMuon->clear();

	// electrons
	//std::cout << "eventNtuple::initialize(): electrons" << std::endl;
	ele_pt->clear();
	ele_conept->clear();
	ele_eta->clear();
	ele_phi->clear();
	ele_E->clear();
	ele_charge->clear();
	ele_dxy->clear();
	ele_dz->clear();
	ele_isfakeablesel->clear();
	ele_ismvasel->clear();
	ele_jetNDauChargedMVASel->clear();
	ele_miniRelIso->clear();
	ele_miniIsoCharged->clear();
	ele_miniIsoNeutral->clear();
	ele_jetPtRel->clear();
	ele_jetPtRatio->clear();
	ele_jetCSV->clear();
	ele_sip3D->clear();
	ele_ntMVAeleID->clear();
	ele_leptonMVA->clear();
	ele_isChargeConsistent->clear();
	ele_passesConversionVeto->clear();
	ele_nMissingHits->clear();
	ele_mcMatchType->clear();
	
	// taus
	//std::cout << "eventNtuple::initialize(): taus" << std::endl;
	tau_pt->clear();
	tau_eta->clear();
	tau_phi->clear();
	tau_E->clear();
	tau_charge->clear();
	tau_dxy->clear();
	tau_dz->clear();
	tau_decayMode->clear();
	tau_decayModeFindingOldDMs->clear();
	tau_decayModeFindingNewDMs->clear();
	tau_byCombinedIsolationDeltaBetaCorr3Hits->clear();
	tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->clear();
	tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->clear();
	tau_byTightCombinedIsolationDeltaBetaCorr3Hits->clear();
	//tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03->clear();
	//tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03->clear();
	//tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03->clear();
	tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->clear();
	tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->clear();
	tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->clear();
	tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->clear();
	tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->clear();
	tau_againstMuonLoose3->clear();
	tau_againstMuonTight3->clear();
	tau_againstElectronVLooseMVA6->clear();
	tau_againstElectronLooseMVA6->clear();
	tau_againstElectronMediumMVA6->clear();
	tau_againstElectronTightMVA6->clear();
	tau_idPreselection->clear();
	tau_idSelection->clear();
	tau_mcMatchType->clear();
	tau_ecalEnergy->clear();
	tau_hcalEnergy->clear();
	//tau_isPFTau->clear();
	//tau_isCaloTau->clear();
	///
	tau_signalChargedHadrCands_pt->clear();
	tau_signalChargedHadrCands_eta->clear();
	tau_signalChargedHadrCands_phi->clear();
	tau_signalChargedHadrCands_E->clear();
	tau_signalNeutrHadrCands_pt->clear();
	tau_signalNeutrHadrCands_eta->clear();
	tau_signalNeutrHadrCands_phi->clear();
	tau_signalNeutrHadrCands_E->clear();
	tau_signalGammaCands_pt->clear();
	tau_signalGammaCands_eta->clear();
	tau_signalGammaCands_phi->clear();
	tau_signalGammaCands_E->clear();
	
	// jets
	//std::cout << "eventNtuple::initialize(): jets" << std::endl;
	jet_pt->clear();
	jet_eta->clear();
	jet_phi->clear();
	jet_E->clear();
	jet_csv->clear();
	jet_flavor->clear();

	// met
	//std::cout << "eventNtuple::initialize(): met" << std::endl;
	PFMET = -9999.;
	PFMETphi = -9999.;
	MHT = -9999.;
	metLD = -9999.;
	METSignificance = -9999.;
	METCov00 = -9999.;
	METCov10 = -9999.;
	METCov01 = -9999.;
	METCov11 = -9999.;
}

void eventNtuple::setup_branches(TTree* tree)
{
	tree->Branch("run", &run);
	tree->Branch("ls", &ls);
	tree->Branch("nEvent", &nEvent);
	tree->Branch("event_weight", &event_weight);
	tree->Branch("PU_weight", &PU_weight);
	tree->Branch("MC_weight", &MC_weight);
	tree->Branch("bTagSF_weight", &bTagSF_weight);
	tree->Branch("leptonSF_weight", &leptonSF_weight);
	tree->Branch("tauSF_weight", &tauSF_weight);
	tree->Branch("triggerSF_weight", &triggerSF_weight);
	tree->Branch("FR_weight", &FR_weight);
	tree->Branch("MC_weight_scale_muF0p5", &MC_weight_scale_muF0p5);
	tree->Branch("MC_weight_scale_muF2", &MC_weight_scale_muF2);
	tree->Branch("MC_weight_scale_muR0p5", &MC_weight_scale_muR0p5);
	tree->Branch("MC_weight_scale_muR2", &MC_weight_scale_muR2);
	tree->Branch("btagSF_weight_LFUp", &btagSF_weight_LFUp);
	tree->Branch("btagSF_weight_LFDown", &btagSF_weight_LFDown);
	tree->Branch("btagSF_weight_HFUp", &btagSF_weight_HFUp);
	tree->Branch("btagSF_weight_HFDown", &btagSF_weight_HFDown);
	tree->Branch("btagSF_weight_HFStats1Up", &btagSF_weight_HFStats1Up);
	tree->Branch("btagSF_weight_HFStats1Down", &btagSF_weight_HFStats1Down);
	tree->Branch("btagSF_weight_HFStats2Up", &btagSF_weight_HFStats2Up);
	tree->Branch("btagSF_weight_HFStats2Down", &btagSF_weight_HFStats2Down);
	tree->Branch("btagSF_weight_LFStats1Up", &btagSF_weight_LFStats1Up);
	tree->Branch("btagSF_weight_LFStats1Down", &btagSF_weight_LFStats1Down);
	tree->Branch("btagSF_weight_LFStats2Up", &btagSF_weight_LFStats2Up);
	tree->Branch("btagSF_weight_LFStats2Down", &btagSF_weight_LFStats2Down);
	tree->Branch("btagSF_weight_cErr1Up", &btagSF_weight_cErr1Up);
	tree->Branch("btagSF_weight_cErr1Down", &btagSF_weight_cErr1Down);
	tree->Branch("btagSF_weight_cErr2Up", &btagSF_weight_cErr2Up);
	tree->Branch("btagSF_weight_cErr2Down", &btagSF_weight_cErr2Down);
	tree->Branch("HiggsDecayType", &HiggsDecayType);
	tree->Branch("lepCategory", &lepCategory);
	tree->Branch("btagCategory", &btagCategory);
	tree->Branch("npuTrue", &npuTrue);
	tree->Branch("npuInTime", &npuInTime);
	tree->Branch("pvx", &pvx);
	tree->Branch("pvy", &pvy);
	tree->Branch("pvz", &pvz);
	//tree->Branch("triggerBits_single_e", &triggerBits_single_e);
	//tree->Branch("triggerBits_single_mu", &triggerBits_single_mu);
	//tree->Branch("triggerBits_double_e", &triggerBits_double_e);
	//tree->Branch("triggerBits_double_mu", &triggerBits_double_mu);
	//tree->Branch("triggerBits_elemu", &triggerBits_elemu);
	tree->Branch("triggerBits", &triggerBits);
	tree->Branch("filterBits", &filterBits);
	tree->Branch("nBadMuons", &nBadMuons);
	
	tree->Branch("passTauCharge", &passTauCharge);
	tree->Branch("isGenMatchedLep", &isGenMatchedLep);
	tree->Branch("isGenMatchedTau", &isGenMatchedTau);

	tree->Branch("n_presel_mu", &n_presel_mu);
	tree->Branch("n_fakeable_mu", &n_fakeable_mu);
	tree->Branch("n_mvasel_mu", &n_mvasel_mu);
	tree->Branch("n_presel_ele", &n_presel_ele);
	tree->Branch("n_fakeable_ele", &n_fakeable_ele);
	tree->Branch("n_mvasel_ele", &n_mvasel_ele);
	tree->Branch("n_presel_tau", &n_presel_tau);
	tree->Branch("n_tau", &n_tau);
	tree->Branch("n_presel_jet", &n_presel_jet);
	tree->Branch("n_jet", &n_jet);
	tree->Branch("n_btag_loose", &n_btag_loose);
	tree->Branch("n_btag_medium", &n_btag_medium);

	tree->Branch("mu_pt",                   &mu_pt);
	tree->Branch("mu_conept",               &mu_conept);
	tree->Branch("mu_eta",                  &mu_eta);
	tree->Branch("mu_phi",                  &mu_phi);
	tree->Branch("mu_E",                    &mu_E);
	tree->Branch("mu_charge",               &mu_charge);
	tree->Branch("mu_jetNDauChargedMVASel", &mu_jetNDauChargedMVASel);
	tree->Branch("mu_miniRelIso",           &mu_miniRelIso);
	tree->Branch("mu_miniIsoCharged",       &mu_miniIsoCharged);
	tree->Branch("mu_miniIsoNeutral",       &mu_miniIsoNeutral);
	tree->Branch("mu_jetPtRel",             &mu_jetPtRel);
	tree->Branch("mu_jetPtRatio",           &mu_jetPtRatio);
	tree->Branch("mu_jetCSV",               &mu_jetCSV);
	tree->Branch("mu_sip3D",                &mu_sip3D);
	tree->Branch("mu_dxy",                  &mu_dxy);
	tree->Branch("mu_dz",                   &mu_dz);
	tree->Branch("mu_segmentCompatibility", &mu_segmentCompatibility);
	tree->Branch("mu_leptonMVA",            &mu_leptonMVA);
	tree->Branch("mu_mediumID",             &mu_mediumID);
	tree->Branch("mu_dpt_div_pt",           &mu_dpt_div_pt);
	tree->Branch("mu_ismvasel",             &mu_ismvasel);
	tree->Branch("mu_isfakeablesel",        &mu_isfakeablesel);
	tree->Branch("mu_mcMatchType",          &mu_mcMatchType);
	tree->Branch("mu_isPFMuon",             &mu_isPFMuon);

	tree->Branch("ele_pt",                   &ele_pt);
	tree->Branch("ele_conept",               &ele_conept);
	tree->Branch("ele_eta",                  &ele_eta);
	tree->Branch("ele_phi",                  &ele_phi);
	tree->Branch("ele_E",                    &ele_E);
	tree->Branch("ele_charge",               &ele_charge);
	tree->Branch("ele_jetNDauChargedMVASel", &ele_jetNDauChargedMVASel);
	tree->Branch("ele_miniRelIso",           &ele_miniRelIso);
	tree->Branch("ele_miniIsoCharged",       &ele_miniIsoCharged);
	tree->Branch("ele_miniIsoNeutral",       &ele_miniIsoNeutral);
	tree->Branch("ele_jetPtRel",             &ele_jetPtRel);
	tree->Branch("ele_jetPtRatio",           &ele_jetPtRatio);
	tree->Branch("ele_jetCSV",               &ele_jetCSV);
	tree->Branch("ele_sip3D",                &ele_sip3D);
	tree->Branch("ele_dxy",                  &ele_dxy);
	tree->Branch("ele_dz",                   &ele_dz);
	tree->Branch("ele_ntMVAeleID",           &ele_ntMVAeleID);
	tree->Branch("ele_leptonMVA",            &ele_leptonMVA);
	tree->Branch("ele_isChargeConsistent",   &ele_isChargeConsistent);
	tree->Branch("ele_passesConversionVeto", &ele_passesConversionVeto);
	tree->Branch("ele_nMissingHits",         &ele_nMissingHits);
	tree->Branch("ele_ismvasel",             &ele_ismvasel);
	tree->Branch("ele_isfakeablesel",        &ele_isfakeablesel);
	tree->Branch("ele_mcMatchType",          &ele_mcMatchType);

	tree->Branch("tau_pt", &tau_pt);
	tree->Branch("tau_eta", &tau_eta);
	tree->Branch("tau_phi", &tau_phi);
	tree->Branch("tau_E", &tau_E);
	tree->Branch("tau_charge", &tau_charge);
	tree->Branch("tau_dxy", &tau_dxy);
	tree->Branch("tau_dz", &tau_dz);
	tree->Branch("tau_decayMode", &tau_decayMode);
	tree->Branch("tau_decayModeFindingOldDMs", &tau_decayModeFindingOldDMs);
	tree->Branch("tau_decayModeFindingNewDMs", &tau_decayModeFindingNewDMs);
	tree->Branch("tau_byCombinedIsolationDeltaBetaCorr3Hits", &tau_byCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
	tree->Branch("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
	//tree->Branch("tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03", &tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03);
	//tree->Branch("tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03", &tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03);
	//tree->Branch("tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03", &tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03);
	tree->Branch("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT", &tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byTightIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
	tree->Branch("tau_againstMuonLoose3", &tau_againstMuonLoose3);
	tree->Branch("tau_againstMuonTight3", &tau_againstMuonTight3);
	tree->Branch("tau_againstElectronVLooseMVA6", &tau_againstElectronVLooseMVA6);
	tree->Branch("tau_againstElectronLooseMVA6", &tau_againstElectronLooseMVA6);
	tree->Branch("tau_againstElectronMediumMVA6", &tau_againstElectronMediumMVA6);
	tree->Branch("tau_againstElectronTightMVA6", &tau_againstElectronTightMVA6);
	tree->Branch("tau_idPreselection", &tau_idPreselection);
	tree->Branch("tau_idSelection", &tau_idSelection);
	tree->Branch("tau_mcMatchType", &tau_mcMatchType);
	tree->Branch("tau_ecalEnergy", &tau_ecalEnergy);
	tree->Branch("tau_hcalEnergy", &tau_hcalEnergy);
	//tree->Branch("tau_isPFTau", &tau_isPFTau);
	//tree->Branch("tau_isCaloTau", &tau_isCaloTau);
	///
	tree->Branch("tau_signalChargedHadrCands_pt", &tau_signalChargedHadrCands_pt);
	tree->Branch("tau_signalChargedHadrCands_eta", &tau_signalChargedHadrCands_eta);
	tree->Branch("tau_signalChargedHadrCands_phi", &tau_signalChargedHadrCands_phi);
	tree->Branch("tau_signalChargedHadrCands_E", &tau_signalChargedHadrCands_E);
	tree->Branch("tau_signalNeutrHadrCands_pt", &tau_signalNeutrHadrCands_pt);
	tree->Branch("tau_signalNeutrHadrCands_eta", &tau_signalNeutrHadrCands_eta);
	tree->Branch("tau_signalNeutrHadrCands_phi", &tau_signalNeutrHadrCands_phi);
	tree->Branch("tau_signalNeutrHadrCands_E", &tau_signalNeutrHadrCands_E);
	tree->Branch("tau_signalGammaCands_pt", &tau_signalGammaCands_pt);
	tree->Branch("tau_signalGammaCands_eta", &tau_signalGammaCands_eta);
	tree->Branch("tau_signalGammaCands_phi", &tau_signalGammaCands_phi);
	tree->Branch("tau_signalGammaCands_E", &tau_signalGammaCands_E);
	///
	
	tree->Branch("jet_pt",     &jet_pt);
	tree->Branch("jet_eta",    &jet_eta);
	tree->Branch("jet_phi",    &jet_phi);
	tree->Branch("jet_E",      &jet_E);
	tree->Branch("jet_csv",    &jet_csv);
	tree->Branch("jet_flavor", &jet_flavor);

	tree->Branch("PFMET", &PFMET);
	tree->Branch("PFMETphi", &PFMETphi);
	tree->Branch("METSignificance", &METSignificance);
	tree->Branch("METCov00", &METCov00);
	tree->Branch("METCov01", &METCov01);
	tree->Branch("METCov10", &METCov10);
	tree->Branch("METCov11", &METCov11);
	tree->Branch("MHT", &MHT);
	tree->Branch("metLD", &metLD);
}

void eventNtuple::set_branch_address(TTree* tree)
{
	tree->SetBranchAddress("run", &run);
	tree->SetBranchAddress("ls", &ls);
	tree->SetBranchAddress("nEvent", &nEvent);
	tree->SetBranchAddress("event_weight", &event_weight);
	tree->SetBranchAddress("PU_weight", &PU_weight);
	tree->SetBranchAddress("MC_weight", &MC_weight);
	tree->SetBranchAddress("bTagSF_weight", &bTagSF_weight);
	tree->SetBranchAddress("leptonSF_weight", &leptonSF_weight);
	tree->SetBranchAddress("tauSF_weight", &tauSF_weight);
	tree->SetBranchAddress("triggerSF_weight", &triggerSF_weight);
	tree->SetBranchAddress("FR_weight", &FR_weight);
	tree->SetBranchAddress("MC_weight_scale_muF0p5", &MC_weight_scale_muF0p5);
	tree->SetBranchAddress("MC_weight_scale_muF2", &MC_weight_scale_muF2);
	tree->SetBranchAddress("MC_weight_scale_muR0p5", &MC_weight_scale_muR0p5);
	tree->SetBranchAddress("MC_weight_scale_muR2", &MC_weight_scale_muR2);
	tree->SetBranchAddress("btagSF_weight_LFUp", &btagSF_weight_LFUp);
	tree->SetBranchAddress("btagSF_weight_LFDown", &btagSF_weight_LFDown);
	tree->SetBranchAddress("btagSF_weight_HFUp", &btagSF_weight_HFUp);
	tree->SetBranchAddress("btagSF_weight_HFDown", &btagSF_weight_HFDown);
	tree->SetBranchAddress("btagSF_weight_HFStats1Up", &btagSF_weight_HFStats1Up);
	tree->SetBranchAddress("btagSF_weight_HFStats1Down", &btagSF_weight_HFStats1Down);
	tree->SetBranchAddress("btagSF_weight_HFStats2Up", &btagSF_weight_HFStats2Up);
	tree->SetBranchAddress("btagSF_weight_HFStats2Down", &btagSF_weight_HFStats2Down);
	tree->SetBranchAddress("btagSF_weight_LFStats1Up", &btagSF_weight_LFStats1Up);
	tree->SetBranchAddress("btagSF_weight_LFStats1Down", &btagSF_weight_LFStats1Down);
	tree->SetBranchAddress("btagSF_weight_LFStats2Up", &btagSF_weight_LFStats2Up);
	tree->SetBranchAddress("btagSF_weight_LFStats2Down", &btagSF_weight_LFStats2Down);
	tree->SetBranchAddress("btagSF_weight_cErr1Up", &btagSF_weight_cErr1Up);
	tree->SetBranchAddress("btagSF_weight_cErr1Down", &btagSF_weight_cErr1Down);
	tree->SetBranchAddress("btagSF_weight_cErr2Up", &btagSF_weight_cErr2Up);
	tree->SetBranchAddress("btagSF_weight_cErr2Down", &btagSF_weight_cErr2Down);
	tree->SetBranchAddress("HiggsDecayType", &HiggsDecayType);
	tree->SetBranchAddress("lepCategory", &lepCategory);
	tree->SetBranchAddress("btagCategory", &btagCategory);
	tree->SetBranchAddress("npuTrue", &npuTrue);
	tree->SetBranchAddress("npuInTime", &npuInTime);
	tree->SetBranchAddress("pvx", &pvx);
	tree->SetBranchAddress("pvy", &pvy);
	tree->SetBranchAddress("pvz", &pvz);
	//tree->SetBranchAddress("triggerBits_single_e", &triggerBits_single_e);
	//tree->SetBranchAddress("triggerBits_single_mu", &triggerBits_single_mu);
	//tree->SetBranchAddress("triggerBits_double_e", &triggerBits_double_e);
	//tree->SetBranchAddress("triggerBits_double_mu", &triggerBits_double_mu);
	//tree->SetBranchAddress("triggerBits_elemu", &triggerBits_elemu);
	tree->SetBranchAddress("triggerBits", &triggerBits);
	tree->SetBranchAddress("filterBits", &filterBits);
	tree->SetBranchAddress("nBadMuons", &nBadMuons);

	tree->SetBranchAddress("passTauCharge", &passTauCharge);
	tree->SetBranchAddress("isGenMatchedLep", &isGenMatchedLep);
	tree->SetBranchAddress("isGenMatchedTau", &isGenMatchedTau);
	
	tree->SetBranchAddress("n_presel_mu", &n_presel_mu);
	tree->SetBranchAddress("n_fakeable_mu", &n_fakeable_mu);
	tree->SetBranchAddress("n_mvasel_mu", &n_mvasel_mu);
	tree->SetBranchAddress("n_presel_ele", &n_presel_ele);
	tree->SetBranchAddress("n_fakeable_ele", &n_fakeable_ele);
	tree->SetBranchAddress("n_mvasel_ele", &n_mvasel_ele);
	tree->SetBranchAddress("n_presel_tau", &n_presel_tau);
	tree->SetBranchAddress("n_tau", &n_tau);
	tree->SetBranchAddress("n_presel_jet", &n_presel_jet);
	tree->SetBranchAddress("n_jet", &n_jet);
	tree->SetBranchAddress("n_btag_loose", &n_btag_loose);
	tree->SetBranchAddress("n_btag_medium", &n_btag_medium);

	tree->SetBranchAddress("mu_pt",                   &mu_pt);
	tree->SetBranchAddress("mu_conept",               &mu_conept);
	tree->SetBranchAddress("mu_eta",                  &mu_eta);
	tree->SetBranchAddress("mu_phi",                  &mu_phi);
	tree->SetBranchAddress("mu_E",                    &mu_E);
	tree->SetBranchAddress("mu_charge",               &mu_charge);
	tree->SetBranchAddress("mu_jetNDauChargedMVASel", &mu_jetNDauChargedMVASel);
	tree->SetBranchAddress("mu_miniRelIso",           &mu_miniRelIso);
	tree->SetBranchAddress("mu_miniIsoCharged",       &mu_miniIsoCharged);
	tree->SetBranchAddress("mu_miniIsoNeutral",       &mu_miniIsoNeutral);
	tree->SetBranchAddress("mu_jetPtRel",             &mu_jetPtRel);
	tree->SetBranchAddress("mu_jetPtRatio",           &mu_jetPtRatio);
	tree->SetBranchAddress("mu_jetCSV",               &mu_jetCSV);
	tree->SetBranchAddress("mu_sip3D",                &mu_sip3D);
	tree->SetBranchAddress("mu_dxy",                  &mu_dxy);
	tree->SetBranchAddress("mu_dz",                   &mu_dz);
	tree->SetBranchAddress("mu_segmentCompatibility", &mu_segmentCompatibility);
	tree->SetBranchAddress("mu_leptonMVA",            &mu_leptonMVA);
	tree->SetBranchAddress("mu_mediumID",             &mu_mediumID);
	tree->SetBranchAddress("mu_dpt_div_pt",           &mu_dpt_div_pt);
	tree->SetBranchAddress("mu_ismvasel",             &mu_ismvasel);
	tree->SetBranchAddress("mu_isfakeablesel",        &mu_isfakeablesel);
	tree->SetBranchAddress("mu_mcMatchType",          &mu_mcMatchType);
	tree->SetBranchAddress("mu_isPFMuon",             &mu_isPFMuon);

	tree->SetBranchAddress("ele_pt",                   &ele_pt);
	tree->SetBranchAddress("ele_conept",               &ele_conept);
	tree->SetBranchAddress("ele_eta",                  &ele_eta);
	tree->SetBranchAddress("ele_phi",                  &ele_phi);
	tree->SetBranchAddress("ele_E",                    &ele_E);
	tree->SetBranchAddress("ele_charge",               &ele_charge);
	tree->SetBranchAddress("ele_jetNDauChargedMVASel", &ele_jetNDauChargedMVASel);
	tree->SetBranchAddress("ele_miniRelIso",           &ele_miniRelIso);
	tree->SetBranchAddress("ele_miniIsoCharged",       &ele_miniIsoCharged);
	tree->SetBranchAddress("ele_miniIsoNeutral",       &ele_miniIsoNeutral);
	tree->SetBranchAddress("ele_jetPtRel",             &ele_jetPtRel);
	tree->SetBranchAddress("ele_jetPtRatio",           &ele_jetPtRatio);
	tree->SetBranchAddress("ele_jetCSV",               &ele_jetCSV);
	tree->SetBranchAddress("ele_sip3D",                &ele_sip3D);
	tree->SetBranchAddress("ele_dxy",                  &ele_dxy);
	tree->SetBranchAddress("ele_dz",                   &ele_dz);
	tree->SetBranchAddress("ele_ntMVAeleID",           &ele_ntMVAeleID);
	tree->SetBranchAddress("ele_leptonMVA",            &ele_leptonMVA);
	tree->SetBranchAddress("ele_isChargeConsistent",   &ele_isChargeConsistent);
	tree->SetBranchAddress("ele_passesConversionVeto", &ele_passesConversionVeto);
	tree->SetBranchAddress("ele_nMissingHits",         &ele_nMissingHits);
	tree->SetBranchAddress("ele_ismvasel",             &ele_ismvasel);
	tree->SetBranchAddress("ele_isfakeablesel",        &ele_isfakeablesel);
	tree->SetBranchAddress("ele_mcMatchType",          &ele_mcMatchType);

	tree->SetBranchAddress("tau_pt", &tau_pt);
	tree->SetBranchAddress("tau_eta", &tau_eta);
	tree->SetBranchAddress("tau_phi", &tau_phi);
	tree->SetBranchAddress("tau_E", &tau_E);
	tree->SetBranchAddress("tau_charge", &tau_charge);
	tree->SetBranchAddress("tau_dxy", &tau_dxy);
	tree->SetBranchAddress("tau_dz", &tau_dz);
	tree->SetBranchAddress("tau_decayMode", &tau_decayMode);
	tree->SetBranchAddress("tau_decayModeFindingOldDMs", &tau_decayModeFindingOldDMs);
	tree->SetBranchAddress("tau_decayModeFindingNewDMs", &tau_decayModeFindingNewDMs);
	tree->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorr3Hits", &tau_byCombinedIsolationDeltaBetaCorr3Hits);
	tree->SetBranchAddress("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
	tree->SetBranchAddress("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
	tree->SetBranchAddress("tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
	//tree->SetBranchAddress("tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03", &tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03);
	//tree->SetBranchAddress("tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03", &tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03);
	//tree->SetBranchAddress("tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03", &tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03);
	tree->SetBranchAddress("tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT);
	tree->SetBranchAddress("tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT", &tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
	tree->SetBranchAddress("tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT", &tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
	tree->SetBranchAddress("tau_byTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byTightIsolationMVArun2v1DBdR03oldDMwLT);
	tree->SetBranchAddress("tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT", &tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
	tree->SetBranchAddress("tau_againstMuonLoose3", &tau_againstMuonLoose3);
	tree->SetBranchAddress("tau_againstMuonTight3", &tau_againstMuonTight3);
	tree->SetBranchAddress("tau_againstElectronVLooseMVA6", &tau_againstElectronVLooseMVA6);
	tree->SetBranchAddress("tau_againstElectronLooseMVA6", &tau_againstElectronLooseMVA6);
	tree->SetBranchAddress("tau_againstElectronMediumMVA6", &tau_againstElectronMediumMVA6);
	tree->SetBranchAddress("tau_againstElectronTightMVA6", &tau_againstElectronTightMVA6);
	tree->SetBranchAddress("tau_idPreselection", &tau_idPreselection);
	tree->SetBranchAddress("tau_idSelection", &tau_idSelection);
	tree->SetBranchAddress("tau_mcMatchType", &tau_mcMatchType);
	tree->SetBranchAddress("tau_ecalEnergy", &tau_ecalEnergy);
	tree->SetBranchAddress("tau_hcalEnergy", &tau_hcalEnergy);
	//tree->SetBranchAddress("tau_isPFTau", &tau_isPFTau);
	//tree->SetBranchAddress("tau_isCaloTau", &tau_isCaloTau);
	
	tree->SetBranchAddress("tau_signalChargedHadrCands_pt", &tau_signalChargedHadrCands_pt);
	tree->SetBranchAddress("tau_signalChargedHadrCands_eta", &tau_signalChargedHadrCands_eta);
	tree->SetBranchAddress("tau_signalChargedHadrCands_phi", &tau_signalChargedHadrCands_phi);
	tree->SetBranchAddress("tau_signalChargedHadrCands_E", &tau_signalChargedHadrCands_E);
	tree->SetBranchAddress("tau_signalNeutrHadrCands_pt", &tau_signalNeutrHadrCands_pt);
	tree->SetBranchAddress("tau_signalNeutrHadrCands_eta", &tau_signalNeutrHadrCands_eta);
	tree->SetBranchAddress("tau_signalNeutrHadrCands_phi", &tau_signalNeutrHadrCands_phi);
	tree->SetBranchAddress("tau_signalNeutrHadrCands_E", &tau_signalNeutrHadrCands_E);
	tree->SetBranchAddress("tau_signalGammaCands_pt", &tau_signalGammaCands_pt);
	tree->SetBranchAddress("tau_signalGammaCands_eta", &tau_signalGammaCands_eta);
	tree->SetBranchAddress("tau_signalGammaCands_phi", &tau_signalGammaCands_phi);
	tree->SetBranchAddress("tau_signalGammaCands_E", &tau_signalGammaCands_E);
	
	tree->SetBranchAddress("jet_pt",     &jet_pt);
	tree->SetBranchAddress("jet_eta",    &jet_eta);
	tree->SetBranchAddress("jet_phi",    &jet_phi);
	tree->SetBranchAddress("jet_E",      &jet_E);
	tree->SetBranchAddress("jet_csv",    &jet_csv);
	tree->SetBranchAddress("jet_flavor", &jet_flavor);

	tree->SetBranchAddress("PFMET", &PFMET);
	tree->SetBranchAddress("PFMETphi", &PFMETphi);
	tree->SetBranchAddress("METSignificance", &METSignificance);
	tree->SetBranchAddress("METCov00", &METCov00);
	tree->SetBranchAddress("METCov01", &METCov01);
	tree->SetBranchAddress("METCov10", &METCov10);
	tree->SetBranchAddress("METCov11", &METCov11);
	tree->SetBranchAddress("MHT", &MHT);
	tree->SetBranchAddress("metLD", &metLD);
}

#endif
