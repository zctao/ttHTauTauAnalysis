#ifndef ttHtautauAnalyzer_Ntuple_cc
#define ttHtautauAnalyzer_Ntuple_cc

#include "ttHTauTauAnalysis/ttHtautauAnalyzer/plugins/ttHtautauAnalyzer.h"

void ttHtautauAnalyzer::write_ntuple_muons(const std::vector<pat::Muon>& muons)
{
	for (auto & mu : muons) {
		evNtuple_.mu_pt->push_back(mu.pt());
		evNtuple_.mu_conept->push_back(mu.userFloat("conePt"));
		evNtuple_.mu_eta->push_back(mu.eta());
		evNtuple_.mu_phi->push_back(mu.phi());
		evNtuple_.mu_E->push_back(mu.energy());
		evNtuple_.mu_charge->push_back(mu.charge());
		evNtuple_.mu_dxy->push_back(mu.userFloat("dxy"));
		evNtuple_.mu_dz->push_back(mu.userFloat("dz"));
		evNtuple_.mu_isfakeablesel->push_back(mu.userInt("isFakeable"));
		evNtuple_.mu_ismvasel->push_back(mu.userInt("isTight"));
		evNtuple_.mu_jetNDauChargedMVASel->push_back(mu.userFloat("nearestJetNDauCharged"));
		evNtuple_.mu_miniRelIso->push_back(mu.userFloat("miniRelIso"));
		evNtuple_.mu_miniIsoCharged->push_back(mu.userFloat("miniAbsIsoCharged"));
		evNtuple_.mu_miniIsoNeutral->push_back(mu.userFloat("miniAbsIsoNeutralcorr"));
		evNtuple_.mu_jetPtRel->push_back(mu.userFloat("nearestJetPtRel"));
		evNtuple_.mu_jetPtRatio->push_back(mu.userFloat("nearestJetPtRatio"));
		evNtuple_.mu_jetCSV->push_back(std::max(mu.userFloat("nearestJetCsv"), 0.f));
		evNtuple_.mu_sip3D->push_back(mu.userFloat("sip3D"));
		evNtuple_.mu_segmentCompatibility->push_back(mu.segmentCompatibility());
		evNtuple_.mu_leptonMVA->push_back(mu.userFloat("leptonMVA"));
		evNtuple_.mu_pfRelIso04->push_back(mu.userFloat("relIso04"));
		evNtuple_.mu_istightcharge->push_back(mu.userInt("isTightCharge"));
		//
		evNtuple_.mu_mediumID->push_back(mu.isMediumMuon());
		//
		evNtuple_.mu_dpt_div_pt->push_back(mu.muonBestTrack()->ptError()/
								 mu.muonBestTrack()->pt());
		if (!isdata_)
			evNtuple_.mu_mcMatchType->push_back(mu.userInt("MCMatchType"));
		evNtuple_.mu_isPFMuon->push_back(mu.isPFMuon());
	}
}

void ttHtautauAnalyzer::write_ntuple_electrons(const std::vector<pat::Electron>& electrons)
{
	for (auto & ele : electrons) {
		evNtuple_.ele_pt->push_back(ele.pt());
		evNtuple_.ele_conept->push_back(ele.userFloat("conePt"));
		evNtuple_.ele_eta->push_back(ele.eta());
		evNtuple_.ele_phi->push_back(ele.phi());
		evNtuple_.ele_E->push_back(ele.energy());
		evNtuple_.ele_charge->push_back(ele.charge());
		evNtuple_.ele_dxy->push_back(ele.userFloat("dxy"));
		evNtuple_.ele_dz->push_back(ele.userFloat("dz"));
		evNtuple_.ele_isfakeablesel->push_back(ele.userInt("isFakeable"));
		evNtuple_.ele_ismvasel->push_back(ele.userInt("isTight"));
		evNtuple_.ele_jetNDauChargedMVASel->
			push_back(ele.userFloat("nearestJetNDauCharged"));
		evNtuple_.ele_miniRelIso->push_back(ele.userFloat("miniRelIso"));
		evNtuple_.ele_miniIsoCharged->push_back(ele.userFloat("miniAbsIsoCharged"));
		evNtuple_.ele_miniIsoNeutral->push_back(ele.userFloat("miniAbsIsoNeutralcorr"));
		evNtuple_.ele_jetPtRel->push_back(ele.userFloat("nearestJetPtRel"));
		evNtuple_.ele_jetPtRatio->push_back(ele.userFloat("nearestJetPtRatio"));
		evNtuple_.ele_jetCSV->push_back(std::max(ele.userFloat("nearestJetCsv"), 0.f));
		evNtuple_.ele_sip3D->push_back(ele.userFloat("sip3D"));
		evNtuple_.ele_ntMVAeleID->push_back(ele.userFloat("eleMvaId"));
		evNtuple_.ele_leptonMVA->push_back(ele.userFloat("leptonMVA"));
		evNtuple_.ele_pfRelIso04->push_back(ele.userFloat("relIso04"));
		evNtuple_.ele_istightcharge->push_back(ele.userInt("isTightCharge"));
		evNtuple_.ele_isChargeConsistent->push_back(ele.isGsfCtfScPixChargeConsistent() + ele.isGsfScPixChargeConsistent() > 1);
		evNtuple_.ele_passesConversionVeto->push_back(ele.passConversionVeto());
		evNtuple_.ele_nMissingHits->push_back(ele.userFloat("numMissingHits"));
		if (!isdata_)
			evNtuple_.ele_mcMatchType->push_back(ele.userInt("MCMatchType"));
	}
}

void ttHtautauAnalyzer::write_ntuple_taus(const std::vector<pat::Tau>& taus)
{
	for (auto & tau : taus) {
		evNtuple_.tau_pt->push_back(tau.pt());
		evNtuple_.tau_eta->push_back(tau.eta());
		evNtuple_.tau_phi->push_back(tau.phi());
		evNtuple_.tau_E->push_back(tau.energy());
		evNtuple_.tau_charge->push_back(tau.charge());
		evNtuple_.tau_dxy->push_back(tau.userFloat("dxy"));
		evNtuple_.tau_dz->push_back(tau.userFloat("dz"));
		evNtuple_.tau_decayMode->push_back(tau.decayMode());
		evNtuple_.tau_decayModeFindingOldDMs->push_back(tau.tauID("decayModeFinding"));
		evNtuple_.tau_decayModeFindingNewDMs->push_back(tau.tauID("decayModeFindingNewDMs"));
		evNtuple_.tau_byCombinedIsolationDeltaBetaCorr3Hits->push_back(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
		evNtuple_.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
		evNtuple_.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
		evNtuple_.tau_byTightCombinedIsolationDeltaBetaCorr3Hits->push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
		//evNtuple_.tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03->push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03"));
		//evNtuple_.tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03->push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03"));
		//evNtuple_.tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03->push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3HitsdR03"));
		evNtuple_.tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->push_back(tau.tauID("byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		evNtuple_.tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->push_back(tau.tauID("byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		evNtuple_.tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->push_back(tau.tauID("byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		evNtuple_.tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->push_back(tau.tauID("byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		evNtuple_.tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->push_back(tau.tauID("byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017"));
		evNtuple_.tau_byIsolationMVArun2v1DBdR03oldDMwLTraw->push_back(tau.tauID("byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017"));
		evNtuple_.tau_againstMuonLoose3->push_back(tau.tauID("againstMuonLoose3"));
		evNtuple_.tau_againstMuonTight3->push_back(tau.tauID("againstMuonTight3"));
		evNtuple_.tau_againstElectronVLooseMVA6->push_back(tau.tauID("againstElectronVLooseMVA6"));
		evNtuple_.tau_againstElectronLooseMVA6->push_back(tau.tauID("againstElectronLooseMVA6"));
		evNtuple_.tau_againstElectronMediumMVA6->push_back(tau.tauID("againstElectronMediumMVA6"));
		evNtuple_.tau_againstElectronTightMVA6->push_back(tau.tauID("againstElectronTightMVA6"));
		//evNtuple_.tau_idPreselection->push_back(tau.userFloat("idPreselection"));
		//evNtuple_.tau_idSelection->push_back(tau.userFloat("idSelection"));
		evNtuple_.tau_idPreselection->push_back(tau.userInt("isLoose"));
		evNtuple_.tau_idSelection->push_back(tau.userInt("isTight"));
		if (!isdata_)
			evNtuple_.tau_mcMatchType->push_back(tau.userInt("MCMatchType"));
		//evNtuple_.tau_isPFTau->push_back(tau.isPFTau());
		//evNtuple_.tau_isCaloTau->push_back(tau.isCaloTau());
		
		evNtuple_.tau_ecalEnergy->push_back(tau.ecalEnergy());
		evNtuple_.tau_hcalEnergy->push_back(tau.hcalEnergy());

		// tau decay sub-structure
		std::vector<float> hpts;
		std::vector<float> hetas;
		std::vector<float> hphis;
		std::vector<float> hEs;
		std::vector<int> hcharges;
		for (const auto & h : tau.signalChargedHadrCands()) {
			hpts.push_back(h->pt());
			hetas.push_back(h->eta());
			hphis.push_back(h->phi());
			hEs.push_back(h->energy());
			hcharges.push_back(h->charge());
		}
		evNtuple_.tau_signalChargedHadrCands_pt->push_back(hpts);
		evNtuple_.tau_signalChargedHadrCands_eta->push_back(hetas);
		evNtuple_.tau_signalChargedHadrCands_phi->push_back(hphis);
		evNtuple_.tau_signalChargedHadrCands_E->push_back(hEs);
		evNtuple_.tau_signalChargedHadrCands_charge->push_back(hcharges);

		std::vector<float> npts;
		std::vector<float> netas;
		std::vector<float> nphis;
		std::vector<float> nEs;
		for (const auto & n : tau.signalNeutrHadrCands()) {
			npts.push_back(n->pt());
			netas.push_back(n->eta());
			nphis.push_back(n->phi());
			nEs.push_back(n->energy());
		}
		evNtuple_.tau_signalNeutrHadrCands_pt->push_back(npts);
		evNtuple_.tau_signalNeutrHadrCands_eta->push_back(netas);
		evNtuple_.tau_signalNeutrHadrCands_phi->push_back(nphis);
		evNtuple_.tau_signalNeutrHadrCands_E->push_back(nEs);

		std::vector<float> gpts;
		std::vector<float> getas;
		std::vector<float> gphis;
		std::vector<float> gEs;
		for (const auto & g : tau.signalGammaCands()) {
			gpts.push_back(g->pt());
			getas.push_back(g->eta());
			gphis.push_back(g->phi());
			gEs.push_back(g->energy());
		}
		evNtuple_.tau_signalGammaCands_pt->push_back(gpts);
		evNtuple_.tau_signalGammaCands_eta->push_back(getas);
		evNtuple_.tau_signalGammaCands_phi->push_back(gphis);
		evNtuple_.tau_signalGammaCands_E->push_back(gEs);
 	}
	
}

void ttHtautauAnalyzer::write_ntuple_jets(const std::vector<pat::Jet>& jets)
{
	for (auto & jet : jets) {
		evNtuple_.jet_pt->push_back(jet.pt());
		evNtuple_.jet_eta->push_back(jet.eta());
		evNtuple_.jet_phi->push_back(jet.phi());
		evNtuple_.jet_E->push_back(jet.energy());
		evNtuple_.jet_csv->push_back(getJetCSV(jet));
		evNtuple_.jet_deepcsv->push_back(getJetCSV(jet));
		evNtuple_.jet_deepCvsB->push_back(getJetDeepCvsB(jet));
		evNtuple_.jet_deepCvsL->push_back(getJetDeepCvsL(jet));
		evNtuple_.jet_csvv2->push_back(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
		evNtuple_.jet_ctagCvsB->push_back(jet.bDiscriminator("pfCombinedCvsBJetTags"));
		evNtuple_.jet_ctagCvsL->push_back(jet.bDiscriminator("pfCombinedCvsLJetTags"));
		evNtuple_.jet_flavor->push_back(jet.hadronFlavour());
		evNtuple_.jet_qgLikelihood->push_back(jet.userFloat("qgLikelihood"));
		evNtuple_.jet_axis2->push_back(jet.userFloat("axis2"));
		evNtuple_.jet_ptD->push_back(jet.userFloat("ptD"));
		evNtuple_.jet_mult->push_back(jet.userInt("mult"));
		evNtuple_.jet_jesUnc->push_back(jet.userFloat("jesUnc"));
	}
}

/*
void ttHtautauAnalyzer::write_ntuple_bTagSF(const std::vector<pat::Jet>& jets)
{
	evNtuple_.bTagSF_weight =
		sf_helper_->Get_EvtCSVWeight(jets,"NA");
	evNtuple_.btagSF_weight_LFUp =
		sf_helper_->Get_EvtCSVWeight(jets,"LFUp");
	evNtuple_.btagSF_weight_LFDown =
		sf_helper_->Get_EvtCSVWeight(jets,"LFDown");
	evNtuple_.btagSF_weight_HFUp =
		sf_helper_->Get_EvtCSVWeight(jets,"HFUp");
	evNtuple_.btagSF_weight_HFDown =
		sf_helper_->Get_EvtCSVWeight(jets,"HFDown");
	evNtuple_.btagSF_weight_HFStats1Up =
		sf_helper_->Get_EvtCSVWeight(jets,"HFStats1Up");
	evNtuple_.btagSF_weight_HFStats1Down =
		sf_helper_->Get_EvtCSVWeight(jets,"HFStats1Down");
	evNtuple_.btagSF_weight_HFStats2Up =
		sf_helper_->Get_EvtCSVWeight(jets,"HFStats2Up");
	evNtuple_.btagSF_weight_HFStats2Down =
		sf_helper_->Get_EvtCSVWeight(jets,"HFStats2Down");
	evNtuple_.btagSF_weight_LFStats1Up =
		sf_helper_->Get_EvtCSVWeight(jets,"LFStats1Up");
	evNtuple_.btagSF_weight_LFStats1Down =
		sf_helper_->Get_EvtCSVWeight(jets,"LFStats1Down");
	evNtuple_.btagSF_weight_LFStats2Up =
		sf_helper_->Get_EvtCSVWeight(jets,"LFStats2Up");
	evNtuple_.btagSF_weight_LFStats2Down =
		sf_helper_->Get_EvtCSVWeight(jets,"LFStats2Down");
	evNtuple_.btagSF_weight_cErr1Up =
		sf_helper_->Get_EvtCSVWeight(jets,"cErr1Up");
	evNtuple_.btagSF_weight_cErr1Down =
		sf_helper_->Get_EvtCSVWeight(jets,"cErr1Down");
	evNtuple_.btagSF_weight_cErr2Up =
		sf_helper_->Get_EvtCSVWeight(jets,"cErr2Up");
	evNtuple_.btagSF_weight_cErr2Down =
		sf_helper_->Get_EvtCSVWeight(jets,"cErr2Down");
}

void ttHtautauAnalyzer::write_ntuple_leptonSF(const std::vector<miniLepton>& leps)
{	
	size_t nleps = 0;
	if (anaType_==Analyze_1l2tau)
		nleps = 1;
	else if (anaType_==Analyze_2lss1tau)
		nleps = 2;
	else if (anaType_==Analyze_3l1tau)
		nleps = 3;
	
	assert(leps.size() >= nleps);

	float lepSF = 1.;
	for (size_t ilep = 0; ilep < nleps; ilep++) {
		lepSF *= sf_helper_->Get_LeptonIDSF(leps[ilep]);
	}

	evNtuple_.leptonSF_weight = lepSF;
}

void ttHtautauAnalyzer::write_ntuple_tauSF(const std::vector<pat::Tau>& taus)
{
	size_t ntaus = 0;
	if (anaType_==Analyze_1l2tau)
		ntaus = 2;
	else if (anaType_==Analyze_2lss1tau)
		ntaus = 1;
	else if (anaType_==Analyze_3l1tau)
		ntaus = 1;

	assert(taus.size() >= ntaus);

	float tauSF = 1.;

	for (size_t itau = 0; itau < ntaus; itau++) {
		bool isGenMatched = evt_selector_->pass_tau_mc_match(taus[itau]);
		tauSF *= sf_helper_->Get_TauIDSF(taus[itau], isGenMatched);
	}

	evNtuple_.tauSF_weight = tauSF;
}

void ttHtautauAnalyzer::write_ntuple_tauSF(const std::vector<miniTau>& taus)
{
	size_t ntaus = 0;
	if (anaType_==Analyze_1l2tau)
		ntaus = 2;
	else if (anaType_==Analyze_2lss1tau)
		ntaus = 1;
	else if (anaType_==Analyze_3l1tau)
		ntaus = 1;

	//assert(taus.size() >= ntaus);

	float tauSF = 1.;
	for (size_t itau = 0; itau < ntaus; itau++) {
		tauSF *= sf_helper_->Get_TauIDSF(taus[itau], taus[itau].isGenMatched());
	}

	evNtuple_.tauSF_weight = tauSF;
}

void ttHtautauAnalyzer::write_ntuple_triggerSF(int lepCategory)
{
	if (anaType_==Analyze_2lss1tau)
		evNtuple_.triggerSF_weight = sf_helper_->Get_HLTSF_2l1tau(lepCategory);
	else if (anaType_==Analyze_3l1tau)
		evNtuple_.triggerSF_weight = sf_helper_->Get_HLTSF_3l1tau();
}

void ttHtautauAnalyzer::write_ntuple_triggerSF(const miniLepton& lep,
											   const std::vector<pat::Tau>& taus,
											   bool LTriggered, bool XTriggered)
{
	assert(anaType_==Analyze_1l2tau);
	//if (selType_==Signal_1l2tau)
	evNtuple_.triggerSF_weight =
		sf_helper_->Get_HLTSF_1l2tau(lep, taus, LTriggered, XTriggered);
}

void ttHtautauAnalyzer::write_ntuple_triggerSF(const miniLepton& lep,
											   const std::vector<miniTau>& taus,
											   bool LTriggered, bool XTriggered)
{
	assert(anaType_==Analyze_1l2tau);
	evNtuple_.triggerSF_weight =
		sf_helper_->Get_HLTSF_1l2tau(lep, taus, LTriggered, XTriggered);
}

void ttHtautauAnalyzer::write_ntuple_frweight(const std::vector<miniLepton>& leps,
											  const std::vector<pat::Tau>& taus)
{
	if (selType_==Application_Flip_2lss1tau)
		evNtuple_.FR_weight = sf_helper_->Get_ChargeFlipWeight(leps, taus);
	else
		evNtuple_.FR_weight = sf_helper_->Get_FR_weight(leps, taus);
}

void ttHtautauAnalyzer::write_ntuple_frweight(const std::vector<miniLepton>& leps,
											  const std::vector<miniTau>& taus)
{
	evNtuple_.FR_weight = sf_helper_->Get_FR_weight(leps, taus);
}
*/

#endif
