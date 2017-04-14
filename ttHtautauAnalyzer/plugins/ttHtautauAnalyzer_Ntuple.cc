#ifndef ttHtautauAnalyzer_Ntuple_cc
#define ttHtautauAnalyzer_Ntuple_cc

#include "ttHTauTauAnalysis/ttHtautauAnalyzer/plugins/ttHtautauAnalyzer.h"


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
	assert(anaType_==Analyze_2lss1tau); // only 2lss1tau for now

	assert(leps.size() >= 2);

	float lepSF = 1.;
	lepSF *= sf_helper_->Get_LeptonIDSF(leps[0]);
	lepSF *= sf_helper_->Get_LeptonIDSF(leps[1]);

	evNtuple_.leptonSF_weight = lepSF;
}

void ttHtautauAnalyzer::write_ntuple_tauSF(const pat::Tau& tau, bool isGenMatched)
{
	evNtuple_.tauSF_weight = sf_helper_->Get_TauIDSF(tau, isGenMatched);
}

void ttHtautauAnalyzer::write_ntuple_triggerSF(int lepCategory)
{
	evNtuple_.triggerSF_weight = sf_helper_->Get_HLTSF(lepCategory);
}

void ttHtautauAnalyzer::write_ntuple_frweight(const std::vector<miniLepton>& leps,
											  const std::vector<pat::Tau>& taus)
{
	assert(leps.size() >= 2);
	assert(leps[0].passFakeableSel() and leps[1].passFakeableSel());
	assert(taus.size() >= 1);
	
	if (selType_==Control_1lfakeable) {
		float f1 = sf_helper_->Get_FakeRate(leps[0]);
		float f2 = sf_helper_->Get_FakeRate(leps[1]);

		float F1 = leps[0].passTightSel() ? -1. : f1/(1-f1);
		float F2 = leps[1].passTightSel() ? -1. : f2/(1-f2);

		if (debug_) {
			std::cout << "lep0 pdgid passTight? : " << leps[0].pdgId() << " "
					  << leps[0].passTightSel() << std::endl;
			std::cout << "f1 F1 : " << f1 << " " << F1 << std::endl;
			std::cout << "lep1 pdgid passTight? : " << leps[1].pdgId() << " "
					  << leps[1].passTightSel() << std::endl;
			std::cout << "f2 F2 : " << f2 << " " << F2 << std::endl;
		}

		evNtuple_.FR_weight = -1. * F1 * F2;
		if (debug_) std::cout << "FR_weight : " << -1. * F1 * F2 << std::endl;
		
	}
	else if (selType_==Control_2los1tau) {
		assert(anaType_==Analyze_2lss1tau);
		float P1_misCharge =
			sf_helper_->Get_EleChargeMisIDProb(leps[0], taus[0].charge());
		float P2_misCharge =
			sf_helper_->Get_EleChargeMisIDProb(leps[1], taus[0].charge());

		// only one of the above two can be non-zero
		//assert(P1_misCharge*P2_misCharge==0.);

		if (debug_) {
			std::cout << "tau charge : " << taus[0].charge() << std::endl;
			std::cout << "lep0 pt conept eta : " << leps[0].pt() << " "
					  << leps[0].conept() << " " << leps[0].eta() << std::endl;
			std::cout << "p1_mischarge : " << P1_misCharge << std::endl;
			std::cout << "lep1 pt conept eta : " << leps[1].pt() << " "
					  << leps[1].conept() << " " << leps[1].eta() << std::endl;
			std::cout << "p2_mischarge : " << P2_misCharge << std::endl;
		}

		evNtuple_.FR_weight = P1_misCharge + P2_misCharge;
	}
}

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
		evNtuple_.mu_miniRelIso->push_back(mu.userFloat("miniIso"));
		evNtuple_.mu_miniIsoCharged->push_back(mu.userFloat("miniAbsIsoCharged"));
		evNtuple_.mu_miniIsoNeutral->push_back(mu.userFloat("miniAbsIsoNeutralcorr"));
		evNtuple_.mu_jetPtRel->push_back(mu.userFloat("nearestJetPtRel"));
		evNtuple_.mu_jetPtRatio->push_back(mu.userFloat("nearestJetPtRatio"));
		evNtuple_.mu_jetCSV->push_back(mu.userFloat("nearestJetCsv"));
		evNtuple_.mu_sip3D->push_back(mu.userFloat("sip3D"));
		evNtuple_.mu_segmentCompatibility->push_back(mu.segmentCompatibility());
		evNtuple_.mu_leptonMVA->push_back(mu.userFloat("leptonMVA"));
		//
		evNtuple_.mu_mediumID->push_back(mu.isMediumMuon());
		//
		evNtuple_.mu_dpt_div_pt->push_back(mu.muonBestTrack()->ptError()/
								 mu.muonBestTrack()->pt());
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
		evNtuple_.ele_miniRelIso->push_back(ele.userFloat("miniIso"));
		evNtuple_.ele_miniIsoCharged->push_back(ele.userFloat("miniAbsIsoCharged"));
		evNtuple_.ele_miniIsoNeutral->push_back(ele.userFloat("miniAbsIsoNeutralcorr"));
		evNtuple_.ele_jetPtRel->push_back(ele.userFloat("nearestJetPtRel"));
		evNtuple_.ele_jetPtRatio->push_back(ele.userFloat("nearestJetPtRatio"));
		evNtuple_.ele_jetCSV->push_back(ele.userFloat("nearestJetCsv"));
		evNtuple_.ele_sip3D->push_back(ele.userFloat("sip3D"));
		evNtuple_.ele_ntMVAeleID->push_back(ele.userFloat("eleMvaId"));
		evNtuple_.ele_leptonMVA->push_back(ele.userFloat("leptonMVA"));
		evNtuple_.ele_isChargeConsistent->push_back(ele.isGsfCtfScPixChargeConsistent() + ele.isGsfScPixChargeConsistent() > 1);
		evNtuple_.ele_passesConversionVeto->push_back(ele.passConversionVeto());
		evNtuple_.ele_nMissingHits->push_back(ele.userFloat("numMissingHits"));
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
		evNtuple_.tau_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03->push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3HitsdR03"));
		evNtuple_.tau_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03->push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3HitsdR03"));
		evNtuple_.tau_byTightCombinedIsolationDeltaBetaCorr3HitsdR03->push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3HitsdR03"));
		evNtuple_.tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->push_back(tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"));
		evNtuple_.tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->push_back(tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT"));
		evNtuple_.tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->push_back(tau.tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT"));
		evNtuple_.tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->push_back(tau.tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT"));
		evNtuple_.tau_againstMuonLoose3->push_back(tau.tauID("againstMuonLoose3"));
		evNtuple_.tau_againstMuonTight3->push_back(tau.tauID("againstMuonTight3"));
		evNtuple_.tau_againstElectronVLooseMVA6->push_back(tau.tauID("againstElectronVLooseMVA6"));
		evNtuple_.tau_againstElectronLooseMVA6->push_back(tau.tauID("againstElectronLooseMVA6"));
		evNtuple_.tau_againstElectronMediumMVA6->push_back(tau.tauID("againstElectronMediumMVA6"));
		evNtuple_.tau_againstElectronTightMVA6->push_back(tau.tauID("againstElectronTightMVA6"));
		evNtuple_.tau_mcMatchType->push_back(tau.userInt("MCMatchType"));
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
		evNtuple_.jet_flavor->push_back(jet.hadronFlavour());
	}
}

#endif
