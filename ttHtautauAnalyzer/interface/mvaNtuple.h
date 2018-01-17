#ifndef mvaNtuple_h
#define mvaNtuple_h

#include "TTree.h"

#include "Types_enum.h"

class mvaNtuple
{
 public:

	mvaNtuple(Analysis_types anaType, bool evaluate, bool doSystematics):
	anatype_(anaType),evaluate_(evaluate), dosystematics_(doSystematics){};
	
	~mvaNtuple(){};

	//void set_branch_address(TTree*);
	void setup_branches(TTree*);

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
	float avg_dr_jet;
	float max_lep_eta;
	float met;
	float mht;
	float mT_met_lep0;
	float lep0_conept;
	float lep1_conept;
	float lep2_conept;
	float dr_leps;
	float tau0_pt;
	float tau1_pt;
	float dr_lep0_tau;
	float dr_lep1_tau;
	float mvis_lep0_tau;
	float mvis_lep1_tau;
	float tt_deltaR;
	float tt_mvis;
	float tt_pt;
	float max_dr_jet;
	float HT;
	int ntags;
	int ntags_loose;

	int tau0_decaymode;
	int tau1_decaymode;
	float tau0_E;
	float tau1_E;
	float tau0_upsilon;
	float tau1_upsilon;

	int taup_decaymode;
	int taum_decaymode;
	float taup_E;
	float taum_E;
	float taup_upsilon;
	float taum_upsilon;

	float mva_ttV;
	float mva_ttbar;

	// Todo
	//float memLR_ttV;
	//float memLR_ttbar;
	//tree_in->AddFriend();
	
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

	// selection flags
	int isGenMatchedTau;
	int HiggsDecayType;

 protected:

	Analysis_types anatype_;
	bool evaluate_;
	bool dosystematics_;
	
};

#endif
