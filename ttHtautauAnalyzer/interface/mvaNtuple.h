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
	void compute_all_variables(const std::vector<miniLepton>&,
							   const std::vector<miniTau>&,
							   const std::vector<miniJet>&,
							   float, float, float, int, int);
	void compute_mva_variables(const std::vector<miniLepton>&,
							   const std::vector<miniTau>&,
							   const std::vector<TLorentzVector>&,
							   float, float, float, int, int);
	void compute_mva_variables(const std::vector<miniLepton>&,
							   const std::vector<miniTau>&,
							   const std::vector<TLorentzVector>&,
							   float, float, float, int, int,
							   const std::vector<TLorentzVector>&);
	void compute_mva_variables(const std::vector<miniLepton>&,
							   const std::vector<miniTau>&,
							   const std::vector<miniJet>&,
							   float, float, float, int, int);
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
	int nEle = -9999;
	int nMu = -9999;
	int nTau = -9999;
	float nJet = -9999.;
	float nbtags_medium = -9999.;
	float nbtags_loose = -9999.;
	float met = -9999.;
	float metLD = -9999.;

	// kinematic variables
	float mll = -9999.;  // dilepton mass
	float lep0_conept = -9999.;
	float lep1_conept = -9999.;
	float lep2_conept = -9999.;
	float lep0_eta = -9999.;
	float lep1_eta = -9999.;
	float lep2_eta = -9999.;
	float lep0_phi = -9999.;
	float lep1_phi = -9999.;
	float lep2_phi = -9999.;
	float lep3_phi = -9999.;
	float lep0_E = -9999.;
	float lep1_E = -9999.;
	float lep2_E = -9999.;
	float tau0_pt = -9999.;
	float tau1_pt = -9999.;
	float tau0_eta = -9999.;
	float tau1_eta = -9999.;
	float tau0_phi = -9999.;
	float tau1_phi = -9999.;
	float tau0_E = -9999.;
	float tau1_E = -9999.;
	
	// Other MVA input variables
	float mindr_lep0_jet = -9999.;
	float mindr_lep1_jet = -9999.;
	float mindr_lep2_jet = -9999.;
	float mindr_tau0_jet = -9999.;
	float mindr_tau1_jet = -9999.;
	float avg_dr_jet = -9999.;
	float max_lep_eta = -9999.;
	float mT_met_lep0 = -9999.;
	float mT_met_lep1 = -9999.;
	float costS_tau = -9999.;
	float dr_leps = -9999.;
	float dr_lep0_tau = -9999.;
	float dr_lep1_tau = -9999.;
	float dr_lep_tau_ss = -9999.;
	float dr_lep_tau_lead = -9999.;
	float dr_lep_tau_sublead = -9999.;
	float mvis_lep0_tau = -9999.;
	float mvis_lep1_tau = -9999.;
	float dr_taus = -9999.;
	float mTauTauVis = -9999.;
	float tt_pt = -9999.;
	float max_dr_jet = -9999.;
	float HT = -9999.;
	float mht = -9999.;
	float mbb = -9999.;
	float is_OS = -9999.;
	float min_dr_lep_jet = -9999.;
	float mindr_tau_jet = -9999.;
	float max_dr_lep_tau = -9999.;
	float min_dr_lep_tau = -9999.;
	float avg_dr_lep_tau = -9999.;

	int tau0_decaymode = -9999;
	int tau1_decaymode = -9999;
	float tau0_easym = -9999.;
	float tau1_easym = -9999.;

	int tau0_tightWP = -9999;
	int tau1_tightWP = -9999;
	float tau0_ldgtrkpt = -9999.;
	float tau1_ldgtrkpt = -9999.;
	float tau0_ldgtrkE = -9999.;
	float tau1_ldgtrkE = -9999.;

	int taup_decaymode = -9999;
	int taum_decaymode = -9999;
	float taup_E = -9999.;
	float taum_E = -9999.;
	float evisTaus_diff = -9999.;
	float evisTaus_sum = -9999.;
	float evisTausAsym = -9999.;
	float taup_easym = -9999.;
	float taum_easym = -9999.;
	float taup_cosPsi = -9999.;
	float taum_cosPsi = -9999.;

	float taup_pt = -9999.;
	float taum_pt = -9999.;
	float taup_ldgtrkpt = -9999.;
	float taum_ldgtrkpt = -9999.;
	float taup_ldgtrkE = -9999.;
	float taum_ldgtrkE = -9999.;
	int taup_tightWP = -9999;
	int taum_tightWP = -9999;
	float taup_upsilon = -9999.;  // using pt
	float taum_upsilon = -9999.;  // using pt
	
	// event weights
	float event_weight = -9999.;
	float event_weight_thu_shape_x1Up = -9999.;
	float event_weight_thu_shape_x1Down = -9999.;
	float event_weight_thu_shape_y1Up = -9999.;
	float event_weight_thu_shape_y1Down = -9999.;
	float event_weight_btag_LFUp = -9999.;
	float event_weight_btag_LFDown = -9999.;
	float event_weight_btag_HFUp = -9999.;
	float event_weight_btag_HFDown = -9999.;
	float event_weight_btag_HFStats1Up = -9999.;
	float event_weight_btag_HFStats1Down = -9999.;
	float event_weight_btag_HFStats2Up = -9999.;
	float event_weight_btag_HFStats2Down = -9999.;
	float event_weight_btag_LFStats1Up = -9999.;
	float event_weight_btag_LFStats1Down = -9999.;
	float event_weight_btag_LFStats2Up = -9999.;
	float event_weight_btag_LFStats2Down = -9999.;
	float event_weight_btag_cErr1Up = -9999.;
	float event_weight_btag_cErr1Down = -9999.;
	float event_weight_btag_cErr2Up = -9999.;
	float event_weight_btag_cErr2Down = -9999.;
	float event_weight_FRjt_normUp = -9999.;
	float event_weight_FRjt_normDown = -9999.;
	float event_weight_FRjt_shapeUp = -9999.;
	float event_weight_FRjt_shapeDown = -9999.;
	float event_weight_FRe_normUp = -9999.;
	float event_weight_FRe_normDown = -9999.;
	float event_weight_FRe_ptUp = -9999.;
	float event_weight_FRe_ptDown = -9999.;
	float event_weight_FRe_beUp = -9999.;
	float event_weight_FRe_beDown = -9999.;
	float event_weight_FRe_bUp = -9999.;
	float event_weight_FRe_bDown = -9999.;
	float event_weight_FRe_ecUp = -9999.;
	float event_weight_FRe_ecDown = -9999.;
	float event_weight_FRm_normUp = -9999.;
	float event_weight_FRm_normDown = -9999.;
	float event_weight_FRm_ptUp = -9999.;
	float event_weight_FRm_ptDown = -9999.;
	float event_weight_FRm_beUp = -9999.;
	float event_weight_FRm_beDown = -9999.;
	float event_weight_FRm_bUp = -9999.;
	float event_weight_FRm_bDown = -9999.;
	float event_weight_FRm_ecUp = -9999.;
	float event_weight_FRm_ecDown = -9999.;

	// scale factors
	float pu_weight = -9999.;
	float mc_weight = -9999.;
	float btag_sf = -9999.;
	float lepid_sf = -9999.;
	float tauid_sf = -9999.;
	float hlt_sf = -9999.;

	double xsection_weight = -9999.;
	double xsection_weight_gen = -9999.;
	
	// selection flags
	int isGenMatchedLep = -9999;
	int isGenMatchedTau = -9999;
	int HiggsDecayType = -9999;

	// TEST
	// Four vector product of all combinations
	double pp1_pp1 = -9999.;
	double pp1_pp2 = -9999.;
	double pp1_pp3 = -9999.;
	double pp1_pm1 = -9999.;
	double pp1_pm2 = -9999.;
	double pp1_pm3 = -9999.;
	double pp2_pp2 = -9999.;
	double pp2_pp3 = -9999.;
	double pp2_pm1 = -9999.;
	double pp2_pm2 = -9999.;
	double pp2_pm3 = -9999.;
	double pp3_pp3 = -9999.;
	double pp3_pm1 = -9999.;
	double pp3_pm2 = -9999.;
	double pp3_pm3 = -9999.;
	double pm1_pm1 = -9999.;
	double pm1_pm2 = -9999.;
	double pm1_pm3 = -9999.;
	double pm2_pm2 = -9999.;
	double pm2_pm3 = -9999.;
	double pm3_pm3 = -9999.;

	// mem likelihood
	float mem_LR = -9999.;
	
	// Hj_tagger
	float Hj_tagger = -9999.;

	// HTT and input variables
	float HTT = -9999.;
	float HadTop_pt = -9999.;
	
	float CSV_b = -9999.;
	float qg_Wj2 = -9999.;
	float pT_bWj1Wj2 = -9999.;
	float pT_Wj2 = -9999.;
	float m_Wj1Wj2 = -9999.;
	float nllKinFit = -9999.;
	float pT_b_o_kinFit_pT_b = -9999.;

	// BDT output
	//float mva_ttbar;
	//float mva_ttV;

	float mva_output = -9999.;
	
	// out-dated
	float mva_1l2tau_BDT1 = -9999.;
	float mva_1l2tau_BDT2 = -9999.;
	float mva_2lss1tau_BDT1 = -9999.;
	float mva_2lss1tau_BDT2 = -9999.;
	float mva_2lss1tau_BDT3 = -9999.;
	float mva_2lss1tau_BDT4 = -9999.;
	float mva_2lss1tau_BDT5 = -9999.;
	float mva_2lss1tau_BDT6 = -9999.;
	float mva_3l1tau_BDT1 = -9999.;
	float mva_3l1tau_BDT2 = -9999.;
	float mva_3l1tau_BDT3 = -9999.;
	float mva_3l1tau_BDT4 = -9999.;
	float mva_2l2tau_BDT1 = -9999.;
	float mva_2l2tau_BDT2 = -9999.;
	float mva_2l2tau_BDT3 = -9999.;
	float mva_2l2tau_BDT4 = -9999.;
	
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
