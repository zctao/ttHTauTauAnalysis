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
#include "Types_enum.h"

class mvaNtuple
{
 public:

    mvaNtuple(Analysis_types anaType, bool doSystematics,
			  const std::string& version="2017") :
	anatype_(anaType),dosystematics_(doSystematics),
		version_(version){};
	
	~mvaNtuple(){};

	//void set_branch_address(TTree*);
	void setup_branches(TTree*);
	void compute_variables(const std::vector<miniLepton>&,
						   const std::vector<miniTau>&,
						   const std::vector<TLorentzVector>&,
						   float, float, float, int, int);

	void compute_tauDecay_variables(const std::vector<miniTau>&, bool test=false);
	float compute_average_dr(const std::vector<TLorentzVector>&);
	float compute_max_dr(const std::vector<TLorentzVector>&);
	float compute_min_dr(const TLorentzVector&, const std::vector<TLorentzVector>&);
	float compute_cosThetaS(const TLorentzVector&);
	float compute_mT_lep(const miniLepton&, float, float);
	float compute_max_lep_eta(const std::vector<miniLepton>&);
	float compute_upsilon(const miniTau&);
	float compute_upsilon_pt(const miniTau&);
	float compute_cosPsi(const miniTau&, float mass=0.139);
	float compute_cosPsi(const TLorentzVector&, const TLorentzVector&,
						 const TLorentzVector&, float mass=0.139);

	//////////////////////////////
	//// variables
	//////////////////////////////

	// event ID
	unsigned int run;
	unsigned int lumi;
	unsigned long long nEvent;
	
	// MVA variables
	int nJet;
	float mindr_lep0_jet;
	float mindr_lep1_jet;
	float mindr_lep2_jet;
	float mindr_tau0_jet;
	float mindr_tau1_jet;
	float avg_dr_jet;
	float max_lep_eta;
	float met;
	float mht;
	float mT_met_lep0;
	float lep0_conept;
	float lep1_conept;
	float lep2_conept;
	float costS_tau;
	float dr_leps;
	float tau0_pt;
	float tau1_pt;
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
	int nbtags_medium;
	int nbtags_loose;

	int tau0_decaymode;
	int tau1_decaymode;
	float tau0_E;
	float tau1_E;
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
	
 protected:

	Analysis_types anatype_;
	bool dosystematics_;
	std::string version_;

	float lam(float, float, float);
};

#endif
