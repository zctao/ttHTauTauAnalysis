#ifndef TriggerHelper_h
#define TriggerHelper_h

#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <vector>
#include <string>
#include <iostream>

class TriggerHelper
{	
 public:
	
	// constructor and destructor
	TriggerHelper(bool verbose = true);
	~TriggerHelper(){};

	// member functions
	void add_trigger_version_number(HLTConfigProvider&);
	unsigned int get_trigger_bits(edm::Handle<edm::TriggerResults>,
								  HLTConfigProvider&);
	unsigned int get_filter_bits(edm::Handle<edm::TriggerResults>,
								  HLTConfigProvider&);
	bool pass_trigger(std::string&);  // Todo
	bool pass_filter(std::string&);   // Todo
	void dump_hlt_paths(HLTConfigProvider&);  // Todo
	
 protected:

	unsigned int encode_bits(edm::Handle<edm::TriggerResults>,
							 HLTConfigProvider&, const std::vector<std::string>&);
	
	static const std::vector<std::string> hlt_paths_;
	static const std::vector<std::string> filter_paths_;

	std::vector<std::string> hlt_paths_version_;

	bool verbose_;
};

#endif
