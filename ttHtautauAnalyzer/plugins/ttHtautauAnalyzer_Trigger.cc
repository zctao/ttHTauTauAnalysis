#ifndef ttHtautauAnalyzer_Trigger_cc
#define ttHtautauAnalyzer_Trigger_cc

#include "ttHTauTauAnalysis/ttHtautauAnalyzer/plugins/ttHtautauAnalyzer.h"

void ttHtautauAnalyzer::Set_up_triggerNames()
{
	trigger_names_no_ver_ = {
		// single electron triggers
		"HLT_Ele27_WPTight_Gsf_v",
		"HLT_Ele27_eta2p1_WPLoose_Gsf",
		"HLT_Ele25_eta2p1_WPTight_Gsf",
		// single muon triggers
		"HLT_IsoMu24_v",
		"HLT_IsoTkMu24_v",
		"HLT_IsoMu22_eta2p1_v",
		"HLT_IsoTkMu22_eta2p1_v",
		"HLT_IsoMu22_v",
		"HLT_IsoTkMu22_v",
		// di-lepton triggers
		"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
		//
		"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",
		"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
		//
		"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
	};
}

void ttHtautauAnalyzer::Set_up_filterNames()
{
	filter_names_ = {
		"Flag_HBHENoiseFilter",
		"Flag_HBHENoiseIsoFilter",
		"Flag_EcalDeadCellTriggerPrimitiveFilter",
		"Flag_goodVertices",
		"Flag_eeBadScFilter",
		"Flag_globalTightHalo2016Filter",
		"Flag_muonBadTrackFilter",
		"Flag_chargedHadronTrackResolutionFilter"
	};
}

void ttHtautauAnalyzer::Add_trigger_versionNumber()
{
	for (const auto & name : trigger_names_no_ver_) {
		
		bool foundPath = false;
		
		for (const auto & triggername : hlt_config_.triggerNames()) {
			if (triggername.substr(0,name.size()) == name) {
				trigger_names_.push_back(triggername);
				foundPath = true;
				break;
			}
		}

		if (not foundPath)
			std::cerr << "WARNING!! Cannot find HLT path "<< name << std::endl;
	}
}

unsigned int ttHtautauAnalyzer::getTriggerBits(edm::Handle<edm::TriggerResults> triggerResults, std::vector<std::string>& names, HLTConfigProvider& config)
{
	unsigned int bits = 0;
	unsigned int iname = 0;

	for (const auto & n : names) {
		unsigned int index = config.triggerIndex(n);

		if (index >= triggerResults->size()) {
			std::cerr << "Failed to find " << n << std::endl;
			continue;
		}

		if (triggerResults->accept(index))
			bits |= 1 << iname;

		++iname;
	}

	return bits;
}

#endif
