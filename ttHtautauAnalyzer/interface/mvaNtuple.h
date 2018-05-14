#ifndef mvaNtuple_h
#define mvaNtuple_h

#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"

#include "miniLepton.h"
#include "miniTau.h"
#include "miniJet.h"
#include "Types_enum.h"
#include "MVAEvaluator.h"

#include "HadTop/HTT_kinfit/interface/HadTopKinFit.h"

class mvaNtuple
{
 public:

    mvaNtuple(Analysis_types anaType, bool doSystematics,
			  const std::string& version="2017", bool doHTT=false, bool eval=false,
			  bool control_region=false) {
		anatype_ = anaType;
		dosystematics_ = doSystematics;
		version_ = version;
		doHTT_ = doHTT;
		evalMVA_ = eval;
		control_ = control_region;
		
		if (doHTT) {
			const std::string tf_fileName = std::string(getenv("CMSSW_BASE")) +
				"/src/HadTop/HTT_kinfit/data/TF_jets.root";
			kinFit_ = new HadTopKinFit(1, tf_fileName);
		}

		if (doHTT_ or evalMVA_) {
			mva_eval_ = new MVAEvaluator(false);
		}
	}
	
	~mvaNtuple(){};

	//void set_branch_address(TTree*);
	void setup_branches(TTree*);
	void compute_mva_variables(const std::vector<miniLepton>&,
							   const std::vector<miniTau>&,
							   const std::vector<TLorentzVector>&,
							   float, float, float, int, int);
	void compute_mva_variables(const std::vector<miniLepton>&,
							   const std::vector<miniTau>&,
							   const std::vector<TLorentzVector>&,
							   float, float, float, int, int,
							   const std::vector<TLorentzVector>&);
	void compute_tauDecay_variables(const std::vector<miniTau>&, bool test=false);
	void compute_HTT_input_variables(const miniJet&, const miniJet&, const miniJet&);
	void compute_HTT(const std::vector<miniJet>&);
	
	void assign_four_momentum(const std::vector<miniLepton>&,
							  const std::vector<miniTau>&);
	int count_electrons(const std::vector<miniLepton>&);
	int count_muons(const std::vector<miniLepton>&);
	
	float compute_average_dr(const std::vector<TLorentzVector>&);
	float compute_average_dr(const std::vector<TLorentzVector>&,
							 const std::vector<TLorentzVector>&);
	float compute_max_dr(const std::vector<TLorentzVector>&);
	float compute_max_dr(const std::vector<TLorentzVector>&,
						 const std::vector<TLorentzVector>&);
	float compute_min_dr(const TLorentzVector&, const std::vector<TLorentzVector>&);
	float compute_min_dr(const std::vector<TLorentzVector>&,
						 const std::vector<TLorentzVector>&);
	float compute_cosThetaS(const TLorentzVector&);
	float compute_mT_lep(const miniLepton&, float, float);
	float compute_max_lep_eta(const std::vector<miniLepton>&);
	float compute_upsilon(const miniTau&);
	float compute_upsilon_pt(const miniTau&);
	float compute_cosPsi(const miniTau&, float mass=0.139);
	float compute_cosPsi(const TLorentzVector&, const TLorentzVector&,
						 const TLorentzVector&, float mass=0.139);
	float compute_mll(const std::vector<miniLepton>&);
	void evaluate_BDTs();

	//////////////////////////////
	//// variables
	//////////////////////////////

	// event ID
	unsigned int run;
	unsigned int lumi;
	unsigned long long nEvent;

	// generic variables
	int nEle;
	int nMu;
	int nTau;
	float nJet;
	float nbtags_medium;
	float nbtags_loose;
	float met;
	float metLD;

	// kinematic variables
	float mll;  // dilepton mass
	float lep0_conept;
	float lep1_conept;
	float lep2_conept;
	float lep0_eta;
	float lep1_eta;
	float lep2_eta;
	float lep0_phi;
	float lep1_phi;
	float lep2_phi;
	float lep3_phi;
	float lep0_E;
	float lep1_E;
	float lep2_E;
	float tau0_pt;
	float tau1_pt;
	float tau0_eta;
	float tau1_eta;
	float tau0_phi;
	float tau1_phi;
	float tau0_E;
	float tau1_E;
	
	// Other MVA input variables
	float mindr_lep0_jet;
	float mindr_lep1_jet;
	float mindr_lep2_jet;
	float mindr_tau0_jet;
	float mindr_tau1_jet;
	float avg_dr_jet;
	float max_lep_eta;
	float mT_met_lep0;
	float mT_met_lep1;
	float costS_tau;
	float dr_leps;
	float dr_lep0_tau;
	float dr_lep1_tau;
	float dr_lep_tau_ss;
	float dr_lep_tau_lead;
	float dr_lep_tau_sublead;
	float mvis_lep0_tau;
	float mvis_lep1_tau;
	float dr_taus;
	float mTauTauVis;
	float tt_pt;
	float max_dr_jet;
	float HT;
	float mht;
	float mbb;
	float is_OS;
	float min_dr_lep_jet;
	float mindr_tau_jet;
	float max_dr_lep_tau;
	float min_dr_lep_tau;
	float avg_dr_lep_tau;

	int tau0_decaymode;
	int tau1_decaymode;
	float tau0_easym;
	float tau1_easym;

	int tau0_tightWP;
	int tau1_tightWP;
	float tau0_ldgtrkpt;
	float tau1_ldgtrkpt;
	float tau0_ldgtrkE;
	float tau1_ldgtrkE;

	int taup_decaymode;
	int taum_decaymode;
	float taup_E;
	float taum_E;
	float evisTaus_diff;
	float evisTaus_sum;
	float evisTausAsym;
	float taup_easym;
	float taum_easym;
	float taup_cosPsi;
	float taum_cosPsi;

	float taup_pt;
	float taum_pt;
	float taup_ldgtrkpt;
	float taum_ldgtrkpt;
	float taup_ldgtrkE;
	float taum_ldgtrkE;
	int taup_tightWP;
	int taum_tightWP;
	float taup_upsilon;  // using pt
	float taum_upsilon;  // using pt
	
	// event weights
	float event_weight;
	float event_weight_thu_shape_x1Up;
	float event_weight_thu_shape_x1Down;
	float event_weight_thu_shape_y1Up;
	float event_weight_thu_shape_y1Down;
	float event_weight_btag_LFUp;
	float event_weight_btag_LFDown;
	float event_weight_btag_HFUp;
	float event_weight_btag_HFDown;
	float event_weight_btag_HFStats1Up;
	float event_weight_btag_HFStats1Down;
	float event_weight_btag_HFStats2Up;
	float event_weight_btag_HFStats2Down;
	float event_weight_btag_LFStats1Up;
	float event_weight_btag_LFStats1Down;
	float event_weight_btag_LFStats2Up;
	float event_weight_btag_LFStats2Down;
	float event_weight_btag_cErr1Up;
	float event_weight_btag_cErr1Down;
	float event_weight_btag_cErr2Up;
	float event_weight_btag_cErr2Down;
	float event_weight_FRjt_normUp;
	float event_weight_FRjt_normDown;
	float event_weight_FRjt_shapeUp;
	float event_weight_FRjt_shapeDown;
	float event_weight_FRe_normUp;
	float event_weight_FRe_normDown;
	float event_weight_FRe_ptUp;
	float event_weight_FRe_ptDown;
	float event_weight_FRe_bUp;
	float event_weight_FRe_bDown;
	float event_weight_FRe_ecUp;
	float event_weight_FRe_ecDown;
	float event_weight_FRm_normUp;
	float event_weight_FRm_normDown;
	float event_weight_FRm_ptUp;
	float event_weight_FRm_ptDown;
	float event_weight_FRm_bUp;
	float event_weight_FRm_bDown;
	float event_weight_FRm_ecUp;
	float event_weight_FRm_ecDown;

	double xsection_weight;
	double xsection_weight_gen;
	
	// selection flags
	int isGenMatchedTau;
	int HiggsDecayType;

	// TEST
	// Four vector product of all combinations
	double pp1_pp1;
	double pp1_pp2;
	double pp1_pp3;
	double pp1_pm1;
	double pp1_pm2;
	double pp1_pm3;
	double pp2_pp2;
	double pp2_pp3;
	double pp2_pm1;
	double pp2_pm2;
	double pp2_pm3;
	double pp3_pp3;
	double pp3_pm1;
	double pp3_pm2;
	double pp3_pm3;
	double pm1_pm1;
	double pm1_pm2;
	double pm1_pm3;
	double pm2_pm2;
	double pm2_pm3;
	double pm3_pm3;

	// mem likelihood
	float mem_LR;
	
	// Hj_tagger
	float Hj_tagger;

	// HTT and input variables
	float HTT;
	float HadTop_pt;
	
	float CSV_b;
	float qg_Wj2;
	float pT_bWj1Wj2;
	float pT_Wj2;
	float m_Wj1Wj2;
	float nllKinFit;
	float pT_b_o_kinFit_pT_b;

	// BDT output
	//float mva_ttbar;
	//float mva_ttV;

	float mva_1l2tau_BDT1;
	float mva_1l2tau_BDT2;
	float mva_2lss1tau_BDT1;
	float mva_2lss1tau_BDT2;
	float mva_2lss1tau_BDT3;
	float mva_2lss1tau_BDT4;
	float mva_2lss1tau_BDT5;
	float mva_2lss1tau_BDT6;
	float mva_3l1tau_BDT1;
	float mva_3l1tau_BDT2;
	float mva_3l1tau_BDT3;
	float mva_3l1tau_BDT4;
	float mva_2l2tau_BDT1;
	float mva_2l2tau_BDT2;
	float mva_2l2tau_BDT3;
	float mva_2l2tau_BDT4;
	
 protected:

	Analysis_types anatype_;
	bool dosystematics_;
	std::string version_;

	bool doHTT_;
	bool evalMVA_;
	HadTopKinFit *kinFit_;
	MVAEvaluator *mva_eval_;

	bool control_;

	float lam(float, float, float);
	int count_leptons(const std::vector<miniLepton>&, int);
};

#endif
