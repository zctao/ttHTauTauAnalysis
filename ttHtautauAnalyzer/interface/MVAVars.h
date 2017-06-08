#ifndef MVAVars_h
#define MVAVars_h

//#include "DataFormats/Math/interface/deltaR.h"
#include "miniLepton.h"
#include "TLorentzVector.h"
#include "TMVA/Reader.h"

#include <vector>
#include <cmath>
#include <algorithm>

class MVAVars
{
 public:

	MVAVars(const std::vector<miniLepton>&, const std::vector<TLorentzVector>&,
			const std::vector<TLorentzVector>&, float, float, float);
	MVAVars(){};
	~MVAVars(){};

	void set_up_tmva_reader();
	
	//int nJet() const {return nJet_;}
	float nJet() const {return nJet_;}
	float mindr_lep0_jet() const {return mindr_lep0_jet_;}
	float mindr_lep1_jet() const {return mindr_lep1_jet_;}
	float avg_dr_jet() const {return avg_dr_jet_;}
	float max_lep_eta() const {return max_lep_eta_;}
	float met() const {return met_;}
	float mht() const {return mht_;}
	float mT_met_lep0() const {return mT_met_lep0_;}
	float lep0_conept() const {return lep0_conept_;}
	float lep1_conept() const {return lep1_conept_;}
	//bool isGenMatched() const {return isGenMatched_;}
	float dr_leps() const {return dr_leps_;}
	float tau_pt() const {return tau_pt_;}
	float dr_lep0_tau() const {return dr_lep0_tau_;}
	float dr_lep1_tau() const {return dr_lep1_tau_;}
	float mvis_lep0_tau() const {return mvis_lep0_tau_;}
	float mvis_lep1_tau() const {return mvis_lep1_tau_;}
	void compute_all_variables(const std::vector<miniLepton>&,
							   const std::vector<TLorentzVector>&,
							   const std::vector<TLorentzVector>&,
							   float, float, float);
	float BDT_2lss1tau_ttV();
	float BDT_2lss1tau_ttbar();
	
 private:

	//int nJet_;
	float nJet_;
	float mindr_lep0_jet_;
	float mindr_lep1_jet_;
	float avg_dr_jet_;
	float max_lep_eta_;
	float met_;
	float mht_;
	float mT_met_lep0_;
	float lep0_conept_;
	float lep1_conept_;
	//bool isGenMatched_;
	float dr_leps_;
	float tau_pt_;
	float dr_lep0_tau_;
	float dr_lep1_tau_;
	float mvis_lep0_tau_;
	float mvis_lep1_tau_;

	float compute_average_dr(const std::vector<TLorentzVector>&);
	float compute_mindr(const TLorentzVector&, const std::vector<TLorentzVector>&);
	
	TMVA::Reader *reader_2lss1tau_ttV;
	TMVA::Reader *reader_2lss1tau_ttbar;
};

#endif
