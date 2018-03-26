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
	bool pass_trilep_triggers(unsigned int);
	
	//bool pass_trigger(std::string&);  // Todo
	//bool pass_filter(std::string&);   // Todo
	//void dump_hlt_paths(HLTConfigProvider&);  // Todo
	
 protected:

	unsigned int encode_bits(edm::Handle<edm::TriggerResults>,
							 HLTConfigProvider&, const std::vector<std::string>&);

	static const std::vector<std::string> filter_paths_;

	static const std::vector<std::string> hlt_paths_e_;
	static const std::vector<std::string> hlt_paths_m_;
	static const std::vector<std::string> hlt_paths_2l_;
	static const std::vector<std::string> hlt_paths_ltau_;
	static const std::vector<std::string> hlt_paths_3l_;
	std::vector<std::string> hlt_paths_;
	std::vector<std::string> hlt_paths_version_;

	Analysis_types anaType_;
	bool verbose_;

	unsigned int bitmask_e_;
	unsigned int bitmask_m_;
	unsigned int bitmask_2l_;
	unsigned int bitmask_ltau_;
	unsigned int bitmask_3l_;
};

#endif
