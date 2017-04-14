#ifndef EventSelector_h
#define EventSelector_h

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
#endif

#include "Types_enum.h"
#include "miniLepton.h"

class EventSelector
{
 public:
	// constructor and destructor
	EventSelector(Analysis_types anatype, Selection_types seltype, bool debug) {
		anaType_ = anatype;
		selType_ = seltype;
		debug_ = debug;
	}
	   
	
	~EventSelector(){};

	// member functions
	bool pass_lepton_number(const std::vector<miniLepton>&,
							const std::vector<miniLepton>&);
	bool pass_tau_number(int);
	bool pass_lepton_pt(const std::vector<miniLepton>&);
	bool pass_pairMass_veto(const std::vector<miniLepton>&);
	bool pass_tight_charge(const std::vector<miniLepton>&);
	bool pass_Zmass_veto(const std::vector<miniLepton>&);
	bool pass_metLD(float, const std::vector<miniLepton>&);
	bool pass_jet_number(int);
	bool pass_btag_number(int, int);
	bool pass_lepton_charge(int, int);
	bool pass_tau_charge(int, const std::vector<miniLepton>&);
	bool pass_lepton_ID(bool, bool);
	bool pass_lep_mc_match(const std::vector<miniLepton>&);
	bool pass_tau_mc_match(const pat::Tau&);
	
 protected:
	Analysis_types  anaType_;
	Selection_types selType_;
	bool debug_;
};

#endif
