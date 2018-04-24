#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/TriggerHelper.h"

const std::vector<std::string> TriggerHelper::hlt_paths_e_ = {
	"HLT_Ele32_WPTight_Gsf_v",
	"HLT_Ele35_WPTight_Gsf_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_m_ = {
	"HLT_IsoMu24_v",
	"HLT_IsoMu27_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_2e_ = {
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_2m_ = {
	//"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_em_ = {
	"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
	"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
	"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_mtau_ = {
	"HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_etau_ = {
	"HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_3e_ = {
	"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_m2e_ = {
	"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_2me_ = {
	"HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v"
};

const std::vector<std::string> TriggerHelper::hlt_paths_3m_ = {
	"HLT_TripleMu_12_10_5_v"
};

// Filters
const std::vector<std::string> TriggerHelper::filter_paths_ = {
	"Flag_goodVertices",
	"Flag_globalTightHalo2016Filter",
	"Flag_HBHENoiseFilter",
	"Flag_HBHENoiseIsoFilter",
	"Flag_EcalDeadCellTriggerPrimitiveFilter",
	"Flag_BadPFMuonFilter",
	"Flag_BadChargedCandidateFilter",
	"Flag_eeBadScFilter",
	"Flag_ecalBadCalibFilter"
};

TriggerHelper::TriggerHelper(Analysis_types analysis, bool verbose)
{
	anaType_ = analysis;
	verbose_ = verbose;
	
	hlt_paths_.clear();

	if (anaType_ == Analyze_1l2tau)
		set_up_paths_1l2tau();
	else if (anaType_ == Analyze_2lss1tau)
		set_up_paths_2lss1tau();
	else if (anaType_ == Analyze_3l1tau)
		set_up_paths_3l1tau();
	else if (anaType_ == Analyze_2l2tau)
		set_up_paths_2l2tau();
}

void TriggerHelper::set_up_paths_1l2tau()
{
	assert(anaType_ == Analyze_1l2tau);
	unsigned int totalsize = hlt_paths_e_.size() + hlt_paths_m_.size() +
		hlt_paths_mtau_.size() + hlt_paths_etau_.size();
	hlt_paths_.reserve(totalsize);

	bitmask_e_ = add_paths(hlt_paths_e_, totalsize);
	bitmask_m_ = add_paths(hlt_paths_m_, totalsize);
	bitmask_mtau_ = add_paths(hlt_paths_mtau_, totalsize);
	bitmask_etau_ = add_paths(hlt_paths_etau_, totalsize);
	
	bitmask_2e_ = 0;
	bitmask_2m_ = 0;
	bitmask_em_ = 0;
	bitmask_3e_ = 0;
	bitmask_m2e_ = 0;
	bitmask_2me_ = 0;
	bitmask_3m_ = 0;

	assert(hlt_paths_.size()==totalsize);
}

void TriggerHelper::set_up_paths_2lss1tau()
{
	assert(anaType_ == Analyze_2lss1tau);
	unsigned int totalsize = hlt_paths_e_.size() + hlt_paths_m_.size() +
		hlt_paths_2e_.size() + hlt_paths_2m_.size() + hlt_paths_em_.size();
	hlt_paths_.reserve(totalsize);

	bitmask_e_ = add_paths(hlt_paths_e_, totalsize);
	bitmask_m_ = add_paths(hlt_paths_m_, totalsize);
	bitmask_2e_ = add_paths(hlt_paths_2e_, totalsize);
	bitmask_2m_ = add_paths(hlt_paths_2m_, totalsize);
	bitmask_em_ = add_paths(hlt_paths_em_, totalsize);

	bitmask_mtau_ = 0;
	bitmask_etau_ = 0;
	bitmask_3e_ = 0;
	bitmask_m2e_ = 0;
	bitmask_2me_ = 0;
	bitmask_3m_ = 0;

	assert(hlt_paths_.size()==totalsize);
}

void TriggerHelper::set_up_paths_3l1tau()
{
	assert(anaType_ == Analyze_3l1tau);
	unsigned int totalsize = hlt_paths_e_.size() + hlt_paths_m_.size() +
		hlt_paths_2e_.size() + hlt_paths_2m_.size() + hlt_paths_em_.size() +
		hlt_paths_3e_.size() + hlt_paths_m2e_.size() + hlt_paths_2me_.size() +
		hlt_paths_3m_.size();
	hlt_paths_.reserve(totalsize);

	bitmask_e_ = add_paths(hlt_paths_e_, totalsize);
	bitmask_m_ = add_paths(hlt_paths_m_, totalsize);
	bitmask_2e_ = add_paths(hlt_paths_2e_, totalsize);
	bitmask_2m_ = add_paths(hlt_paths_2m_, totalsize);
	bitmask_em_ = add_paths(hlt_paths_em_, totalsize);
    bitmask_3e_ = add_paths(hlt_paths_3e_, totalsize);
	bitmask_m2e_ = add_paths(hlt_paths_m2e_, totalsize);
	bitmask_2me_ = add_paths(hlt_paths_2me_, totalsize);
	bitmask_3m_ = add_paths(hlt_paths_3m_, totalsize);

	bitmask_mtau_ = 0;
	bitmask_etau_ = 0;

	assert(hlt_paths_.size()==totalsize);
}

void TriggerHelper::set_up_paths_2l2tau()
{
	assert(anaType_ == Analyze_2l2tau);
	unsigned int totalsize = hlt_paths_e_.size() + hlt_paths_m_.size() +
		hlt_paths_2e_.size() + hlt_paths_2m_.size() + hlt_paths_em_.size() +
		hlt_paths_mtau_.size() + hlt_paths_etau_.size();
	hlt_paths_.reserve(totalsize);

	bitmask_e_ = add_paths(hlt_paths_e_, totalsize);
	bitmask_m_ = add_paths(hlt_paths_m_, totalsize);
	bitmask_2e_ = add_paths(hlt_paths_2e_, totalsize);
	bitmask_2m_ = add_paths(hlt_paths_2m_, totalsize);
	bitmask_em_ = add_paths(hlt_paths_em_, totalsize);
	bitmask_mtau_ = add_paths(hlt_paths_mtau_, totalsize);
	bitmask_etau_ = add_paths(hlt_paths_etau_, totalsize);
	
	bitmask_3e_ = 0;
	bitmask_m2e_ = 0;
	bitmask_2me_ = 0;
	bitmask_3m_ = 0;

	assert(hlt_paths_.size()==totalsize);
}

unsigned int TriggerHelper::add_paths(const std::vector<std::string>& pathsToAdd,
									  unsigned int totalHLTPathSize)
{
	if (hlt_paths_.size() >= totalHLTPathSize) {
		std::cout << "hlt_paths_ size = " << hlt_paths_.size();
		std::cout << " max size = " << totalHLTPathSize << std::endl;
		std::cout << "Cannot add paths : ";
		for (const auto & path : pathsToAdd)
			std::cout << path << " ";
		std::cout << std::endl;
		return 0;
	}

	hlt_paths_.insert(hlt_paths_.end(), pathsToAdd.begin(), pathsToAdd.end());
	unsigned int bitmask =
		((1<<pathsToAdd.size())-1) << (totalHLTPathSize - hlt_paths_.size());

	return bitmask;
}

void TriggerHelper::add_trigger_version_number(HLTConfigProvider& hlt_config)
{
	hlt_paths_version_.clear();
	hlt_paths_version_.reserve(hlt_paths_.size());
	
	for (const auto & name : hlt_paths_) {
		
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
	// make sure add_trigger_version_number() is called before this
	if (hlt_paths_version_.size()!=hlt_paths_.size()) {
		std::cout << "hlt_paths_ : " << std::endl;
		for (auto p : hlt_paths_) {
			std::cout << p << std::endl;
		}
		std::cout << "hlt_paths_version_ : " << std::endl;
		for (auto p : hlt_paths_version_) {
			std::cout << p << std::endl;
		}
	}
	
	assert(hlt_paths_version_.size()==hlt_paths_.size());
	
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

		if (index >= results->size()) {
			if (verbose_) std::cerr << "Failed to find " << name << std::endl;
		}
		else {
			if (results->accept(index)) bits |= 1 << iname;
		}
		
		++iname;
	}

	return bits;
}

bool TriggerHelper::pass_leptau_cross_triggers(unsigned int triggerBits)
{
	return pass_etau_triggers(triggerBits) or pass_mtau_triggers(triggerBits);
}

bool TriggerHelper::pass_single_lep_triggers(unsigned int triggerBits)
{
	return pass_single_e_triggers(triggerBits) or pass_single_m_triggers(triggerBits);
}

bool TriggerHelper::pass_dilep_triggers(unsigned int triggerBits)
{
	return pass_2e_triggers(triggerBits) or pass_2m_triggers(triggerBits) or
		pass_em_triggers(triggerBits);
}

bool TriggerHelper::pass_trilep_triggers(unsigned int triggerBits)
{
	return pass_3e_triggers(triggerBits) or pass_m_2e_triggers(triggerBits) or
		pass_2m_e_triggers(triggerBits) or pass_3m_triggers(triggerBits);
}
