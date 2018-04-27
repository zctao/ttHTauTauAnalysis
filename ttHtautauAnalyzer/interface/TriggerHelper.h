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
	void set_up_paths_1l2tau();
	void set_up_paths_2lss1tau();
	void set_up_paths_3l1tau();
	void set_up_paths_2l2tau();
	void set_up_paths_all();
	void add_trigger_version_number(HLTConfigProvider&);
	unsigned int get_trigger_bits(edm::Handle<edm::TriggerResults>,
								  HLTConfigProvider&);
	unsigned int get_filter_bits(edm::Handle<edm::TriggerResults>,
								  HLTConfigProvider&);
	unsigned int add_paths(const std::vector<std::string>&, unsigned int);
	
	bool pass_etau_triggers(unsigned int tbits) {return tbits & bitmask_etau_;}
	bool pass_mtau_triggers(unsigned int tbits) {return tbits & bitmask_mtau_;}
	bool pass_leptau_cross_triggers(unsigned int);
	bool pass_single_e_triggers(unsigned int tbits) {return tbits & bitmask_e_;}
	bool pass_single_m_triggers(unsigned int tbits) {return tbits & bitmask_m_;}
	bool pass_single_lep_triggers(unsigned int);
	bool pass_2e_triggers(unsigned int tbits) {return tbits & bitmask_2e_;}
	bool pass_2m_triggers(unsigned int tbits) {return tbits & bitmask_2m_;}
	bool pass_em_triggers(unsigned int tbits) {return tbits & bitmask_em_;}
	bool pass_dilep_triggers(unsigned int);
	bool pass_3e_triggers(unsigned int tbits) {return tbits & bitmask_3e_;}
	bool pass_m_2e_triggers(unsigned int tbits) {return tbits & bitmask_m2e_;}
	bool pass_2m_e_triggers(unsigned int tbits) {return tbits & bitmask_2me_;}
	bool pass_3m_triggers(unsigned int tbits) {return tbits & bitmask_3m_;}
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
	static const std::vector<std::string> hlt_paths_2e_;
	static const std::vector<std::string> hlt_paths_2m_;
	static const std::vector<std::string> hlt_paths_em_;
	static const std::vector<std::string> hlt_paths_mtau_;
	static const std::vector<std::string> hlt_paths_etau_;
	static const std::vector<std::string> hlt_paths_3e_;
	static const std::vector<std::string> hlt_paths_m2e_;
	static const std::vector<std::string> hlt_paths_2me_;
	static const std::vector<std::string> hlt_paths_3m_;
	std::vector<std::string> hlt_paths_;
	std::vector<std::string> hlt_paths_version_;

	Analysis_types anaType_;
	bool verbose_;

	unsigned int bitmask_e_;
	unsigned int bitmask_m_;
	unsigned int bitmask_2e_;
	unsigned int bitmask_2m_;
	unsigned int bitmask_em_;
	unsigned int bitmask_mtau_;
	unsigned int bitmask_etau_;
	unsigned int bitmask_3e_;
	unsigned int bitmask_m2e_;
	unsigned int bitmask_2me_;
	unsigned int bitmask_3m_;
};

#endif
