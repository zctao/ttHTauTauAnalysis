#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/TriggerHelper.h"

const std::vector<std::string> TriggerHelper::hlt_paths_2l1tau_ = {
	// single electron triggers
	"HLT_Ele27_WPTight_Gsf_v",
	"HLT_Ele27_eta2p1_WPLoose_Gsf_v",
	"HLT_Ele25_eta2p1_WPTight_Gsf_v",
	// single muon triggers
	"HLT_IsoMu24_v",
	"HLT_IsoTkMu24_v",
	"HLT_IsoMu22_eta2p1_v",
	"HLT_IsoTkMu22_eta2p1_v",
	"HLT_IsoMu22_v",
	"HLT_IsoTkMu22_v",
	// di-lepton triggers
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
	//
	"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",
	"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
	//
	"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_1l2tau_ = {
	// single lepton triggers
	"HLT_Ele25_eta2p1_WPTight_Gsf_v",
	"HLT_IsoMu22_eta2p1_v",
	"HLT_IsoTkMu22_eta2p1_v",
	"HLT_IsoMu22_v",
	"HLT_IsoTkMu22_v",
	// lepton+tau cross triggers
	"HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v",
	"HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v",
	"HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v",
	"HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v"
};

const std::vector<std::string> TriggerHelper::filter_paths_ = {
	"Flag_HBHENoiseFilter",
	"Flag_HBHENoiseIsoFilter",
	"Flag_EcalDeadCellTriggerPrimitiveFilter",
	"Flag_goodVertices",
	"Flag_eeBadScFilter",
	"Flag_globalTightHalo2016Filter",
	"Flag_muonBadTrackFilter",
	"Flag_chargedHadronTrackResolutionFilter"
};

TriggerHelper::TriggerHelper(Analysis_types analysis, bool verbose)
{
	anaType_ = analysis;
	verbose_ = verbose;
}

void TriggerHelper::add_trigger_version_number(HLTConfigProvider& hlt_config)
{
	if (anaType_ == Analyze_2lss1tau)
		add_trigger_version_number(hlt_config, hlt_paths_2l1tau_);
	else if (anaType_ == Analyze_1l2tau)
		add_trigger_version_number(hlt_config, hlt_paths_1l2tau_);
	//else if (anaType_ == Analyze_3l1tau)
	//	add_trigger_version_number();
}											   
// overload
void TriggerHelper::add_trigger_version_number(
	     HLTConfigProvider& hlt_config, const std::vector<std::string>& hlt_paths)
{
	hlt_paths_version_.clear();
	hlt_paths_version_.reserve(hlt_paths.size());
	
	for (const auto & name : hlt_paths) {
		
		bool foundPath = false;

		for (const auto & triggername : hlt_config.triggerNames()) {
			if (triggername.substr(0,name.size()) == name) {
				hlt_paths_version_.push_back(triggername);
				foundPath = true;
				break;
			}
		} // hlt_config:triggerNames()

		if ((not foundPath) and verbose_)
			std::cerr << "WARNING!!<< Cannot find HLT path " << name << std::endl;
		
	} // hlt_paths
}

unsigned int TriggerHelper::get_trigger_bits(edm::Handle<edm::TriggerResults> triggerResults, HLTConfigProvider& config)
{
	//assert(hlt_paths_version_.size()==hlt_paths_.size());
	return encode_bits(triggerResults, config, hlt_paths_version_);
}

unsigned int TriggerHelper::get_filter_bits(edm::Handle<edm::TriggerResults> filterResults, HLTConfigProvider& config)
{
	return encode_bits(filterResults, config, filter_paths_);
}

unsigned int TriggerHelper::encode_bits(edm::Handle<edm::TriggerResults> results, HLTConfigProvider& config, const std::vector<std::string>& vnames)
{
	unsigned int bits = 0;
	unsigned int iname = 0;
	
	for (const auto & name : vnames) {
		unsigned int index = config.triggerIndex(name);

		if ((index >= results->size()) and verbose_) {
			std::cerr << "Failed to find " << name << std::endl;
			++iname;
			continue;
		}
		if (results->accept(index)) bits |= 1 << iname;

		++iname;
	}

	return bits;
}
