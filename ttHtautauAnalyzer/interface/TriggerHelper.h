#ifndef TriggerHelper_h
#define TriggerHelper_h

#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/Types_enum.h"

#include <vector>
#include <string>
#include <iostream>

class TriggerHelper
{	
 public:
	
	// constructor and destructor
	TriggerHelper(Analysis_types, bool verbose = true);
	~TriggerHelper(){};

	// member functions
	void add_trigger_version_number(HLTConfigProvider&);
	unsigned int get_trigger_bits(edm::Handle<edm::TriggerResults>,
								  HLTConfigProvider&);
	unsigned int get_filter_bits(edm::Handle<edm::TriggerResults>,
								  HLTConfigProvider&);

	bool pass_leptau_cross_triggers(unsigned int);
	bool pass_single_lep_triggers(unsigned int);
	bool pass_single_e_triggers(unsigned int);
	bool pass_single_m_triggers(unsigned int);
	bool pass_dilep_triggers(unsigned int);
	
	//bool pass_trigger(std::string&);  // Todo
	//bool pass_filter(std::string&);   // Todo
	//void dump_hlt_paths(HLTConfigProvider&);  // Todo
	
 protected:

	unsigned int encode_bits(edm::Handle<edm::TriggerResults>,
							 HLTConfigProvider&, const std::vector<std::string>&);

	static const std::vector<std::string> filter_paths_;
	
	static const std::vector<std::string> hlt_paths_e_2l1tau_;
	static const std::vector<std::string> hlt_paths_m_2l1tau_;
	static const std::vector<std::string> hlt_paths_2l_2l1tau_;
	static const std::vector<std::string> hlt_paths_l_1l2tau_;
	static const std::vector<std::string> hlt_paths_x_1l2tau_;
	std::vector<std::string> hlt_paths_;
	std::vector<std::string> hlt_paths_version_;

	Analysis_types anaType_;
	bool verbose_;

	unsigned int bitmask_e_2l1tau_;
	unsigned int bitmask_m_2l1tau_;
	unsigned int bitmask_2l_2l1tau_;
	unsigned int bitmask_l_1l2tau_;
	unsigned int bitmask_x_1l2tau_;
};

#endif
