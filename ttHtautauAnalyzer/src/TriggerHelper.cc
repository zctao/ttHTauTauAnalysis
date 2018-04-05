#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/TriggerHelper.h"

const std::vector<std::string> TriggerHelper::hlt_paths_e_ = {
    "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele35_WPTight_Gsf_v"};

const std::vector<std::string> TriggerHelper::hlt_paths_m_ = {"HLT_IsoMu24_v",
                                                              "HLT_IsoMu27_v"};

const std::vector<std::string> TriggerHelper::hlt_paths_2l_ = {
    // di e
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    // di mu
    //"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v4",
    // e mu
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"};

const std::vector<std::string> TriggerHelper::hlt_paths_ltau_ = {
    "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v",
    "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"};

const std::vector<std::string> TriggerHelper::hlt_paths_3l_ = {
    "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v",
    "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v",
    "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v", "HLT_TripleMu_12_10_5_v"};

// Filters
const std::vector<std::string> TriggerHelper::filter_paths_ = {
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_goodVertices",
    "Flag_eeBadScFilter",
    "Flag_globalTightHalo2016Filter",
    "Flag_muonBadTrackFilter",
    "Flag_chargedHadronTrackResolutionFilter"};

TriggerHelper::TriggerHelper(Analysis_types analysis, bool verbose)
{
    anaType_ = analysis;
    verbose_ = verbose;

    hlt_paths_.clear();

    if (anaType_ == Analyze_2lss1tau) {
        hlt_paths_.reserve(hlt_paths_e_.size() + hlt_paths_m_.size() +
                           hlt_paths_2l_.size());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_e_.begin(),
                          hlt_paths_e_.end());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_m_.begin(),
                          hlt_paths_m_.end());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_2l_.begin(),
                          hlt_paths_2l_.end());

        bitmask_e_ = ((1 << hlt_paths_e_.size()) - 1)
                     << (hlt_paths_m_.size() + hlt_paths_2l_.size());
        bitmask_m_ = ((1 << hlt_paths_m_.size()) - 1) << hlt_paths_2l_.size();
        bitmask_2l_ = (1 << hlt_paths_2l_.size()) - 1;
        bitmask_ltau_ = 0;
        bitmask_3l_ = 0;
    }

    if (anaType_ == Analyze_1l2tau) {
        hlt_paths_.reserve(hlt_paths_e_.size() + hlt_paths_m_.size() +
                           hlt_paths_ltau_.size());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_e_.begin(),
                          hlt_paths_e_.end());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_m_.begin(),
                          hlt_paths_m_.end());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_ltau_.begin(),
                          hlt_paths_ltau_.end());

        bitmask_e_ = ((1 << hlt_paths_e_.size()) - 1)
                     << (hlt_paths_m_.size() + hlt_paths_ltau_.size());
        bitmask_m_ = ((1 << hlt_paths_m_.size()) - 1)
                     << (hlt_paths_ltau_.size());
        bitmask_ltau_ = (1 << hlt_paths_ltau_.size()) - 1;
        bitmask_2l_ = 0;
        bitmask_3l_ = 0;
    }

    if (anaType_ == Analyze_3l1tau) {
        hlt_paths_.reserve(hlt_paths_e_.size() + hlt_paths_m_.size() +
                           hlt_paths_2l_.size() + hlt_paths_3l_.size());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_e_.begin(),
                          hlt_paths_e_.end());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_m_.begin(),
                          hlt_paths_m_.end());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_2l_.begin(),
                          hlt_paths_2l_.end());
        hlt_paths_.insert(hlt_paths_.end(), hlt_paths_3l_.begin(),
                          hlt_paths_3l_.end());

        bitmask_e_ = ((1 << hlt_paths_e_.size()) - 1)
                     << (hlt_paths_m_.size() + hlt_paths_2l_.size() +
                         hlt_paths_3l_.size());
        bitmask_m_ = ((1 << hlt_paths_m_.size()) - 1)
                     << (hlt_paths_2l_.size() + hlt_paths_3l_.size());
        bitmask_2l_ = ((1 << hlt_paths_2l_.size()) - 1)
                      << (hlt_paths_3l_.size());
        bitmask_3l_ = (1 << hlt_paths_3l_.size()) - 1;
        bitmask_ltau_ = 0;
    }
}

void TriggerHelper::add_trigger_version_number(HLTConfigProvider &hlt_config)
{
    hlt_paths_version_.clear();
    hlt_paths_version_.reserve(hlt_paths_.size());

    for (const auto &name : hlt_paths_) {

        bool foundPath = false;

        for (const auto &triggername : hlt_config.triggerNames()) {
            if (triggername.substr(0, name.size()) == name) {
                hlt_paths_version_.push_back(triggername);
                foundPath = true;
                break;
            }
        } // hlt_config:triggerNames()

        if ((not foundPath) and verbose_)
            std::cerr << "WARNING!!<< Cannot find HLT path " << name
                      << std::endl;

    } // hlt_paths
}

unsigned int
TriggerHelper::get_trigger_bits(edm::Handle<edm::TriggerResults> triggerResults,
                                HLTConfigProvider &config)
{
    // assert(hlt_paths_version_.size()==hlt_paths_.size());
    return encode_bits(triggerResults, config, hlt_paths_version_);
}

unsigned int
TriggerHelper::get_filter_bits(edm::Handle<edm::TriggerResults> filterResults,
                               HLTConfigProvider &config)
{
    return encode_bits(filterResults, config, filter_paths_);
}

unsigned int
TriggerHelper::encode_bits(edm::Handle<edm::TriggerResults> results,
                           HLTConfigProvider &config,
                           const std::vector<std::string> &vnames)
{
    unsigned int bits = 0;
    unsigned int iname = 0;

    for (const auto &name : vnames) {
        unsigned int index = config.triggerIndex(name);

        if (index >= results->size()) {
            if (verbose_)
                std::cerr << "Failed to find " << name << std::endl;
        } else {
            if (results->accept(index))
                bits |= 1 << iname;
        }

        ++iname;
    }

    return bits;
}

bool TriggerHelper::pass_leptau_cross_triggers(unsigned int triggerBits)
{
    if (anaType_ == Analyze_1l2tau)
        return triggerBits & bitmask_ltau_;
    else
        return false;
}

bool TriggerHelper::pass_single_e_triggers(unsigned int triggerBits)
{
    // assert(anaType_==Analyze_2lss1tau);
    return triggerBits & bitmask_e_;
}

bool TriggerHelper::pass_single_m_triggers(unsigned int triggerBits)
{
    // assert(anaType_==Analyze_2lss1tau);
    return triggerBits & bitmask_m_;
}

bool TriggerHelper::pass_single_lep_triggers(unsigned int triggerBits)
{
    // assert(anaType_==Analyze_1l2tau);
    return pass_single_e_triggers(triggerBits) or
           pass_single_m_triggers(triggerBits);
}

bool TriggerHelper::pass_dilep_triggers(unsigned int triggerBits)
{
    assert(anaType_ == Analyze_2lss1tau or anaType_ == Analyze_3l1tau);
    return triggerBits & bitmask_2l_;
}

bool TriggerHelper::pass_trilep_triggers(unsigned int triggerBits)
{
    assert(anaType_ == Analyze_3l1tau);
    return triggerBits & bitmask_3l_;
}
