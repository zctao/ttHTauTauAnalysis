#ifndef mvaNtuple_cc
#define mvaNtuple_cc

#include "../interface/mvaNtuple.h"

void mvaNtuple::setup_branches(TTree* tree)
{
	tree->Branch("event_weight", &event_weight);
	if (dosystematics_) {
		tree->Branch("event_weight_thu_shape_x1Up",&event_weight_thu_shape_x1Up);
		tree->Branch("event_weight_thu_shape_x1Down",&event_weight_thu_shape_x1Down);
		tree->Branch("event_weight_thu_shape_y1Up",&event_weight_thu_shape_y1Up);
		tree->Branch("event_weight_thu_shape_y1Down",&event_weight_thu_shape_y1Down);
		tree->Branch("event_weight_btag_LFUp",&event_weight_btag_LFUp);
		tree->Branch("event_weight_btag_LFDown",&event_weight_btag_LFDown);
		tree->Branch("event_weight_btag_HFUp",&event_weight_btag_HFUp);
		tree->Branch("event_weight_btag_HFDown",&event_weight_btag_HFDown);
		tree->Branch("event_weight_btag_HFStats1Up",&event_weight_btag_HFStats1Up);
		tree->Branch("event_weight_btag_HFStats1Down",&event_weight_btag_HFStats1Down);
		tree->Branch("event_weight_btag_HFStats2Up",&event_weight_btag_HFStats2Up);
		tree->Branch("event_weight_btag_HFStats2Down",&event_weight_btag_HFStats2Down);
		tree->Branch("event_weight_btag_LFStats1Up",&event_weight_btag_LFStats1Up);
		tree->Branch("event_weight_btag_LFStats1Down",&event_weight_btag_LFStats1Down);
		tree->Branch("event_weight_btag_LFStats2Up",&event_weight_btag_LFStats2Up);
		tree->Branch("event_weight_btag_LFStats2Down",&event_weight_btag_LFStats2Down);
		tree->Branch("event_weight_btag_cErr1Up",&event_weight_btag_cErr1Up);
		tree->Branch("event_weight_btag_cErr1Down",&event_weight_btag_cErr1Down);
		tree->Branch("event_weight_btag_cErr2Up",&event_weight_btag_cErr2Up);
		tree->Branch("event_weight_btag_cErr2Down",&event_weight_btag_cErr2Down);
		tree->Branch("event_weight_FRjt_normUp",&event_weight_FRjt_normUp);
		tree->Branch("event_weight_FRjt_normDown",&event_weight_FRjt_normDown);
		tree->Branch("event_weight_FRjt_shapeUp",&event_weight_FRjt_shapeUp);
		tree->Branch("event_weight_FRjt_shapeDown",&event_weight_FRjt_shapeDown);
		tree->Branch("event_weight_FRe_normUp",&event_weight_FRe_normUp);
		tree->Branch("event_weight_FRe_normDown",&event_weight_FRe_normDown);
		tree->Branch("event_weight_FRe_ptUp",&event_weight_FRe_ptUp);
		tree->Branch("event_weight_FRe_ptDown",&event_weight_FRe_ptDown);
		tree->Branch("event_weight_FRe_bUp",&event_weight_FRe_bUp);
		tree->Branch("event_weight_FRe_bDown",&event_weight_FRe_bDown);
		tree->Branch("event_weight_FRe_ecUp",&event_weight_FRe_ecUp);
		tree->Branch("event_weight_FRe_ecDown",&event_weight_FRe_ecDown);
		tree->Branch("event_weight_FRm_normUp",&event_weight_FRm_normUp);
		tree->Branch("event_weight_FRm_normDown",&event_weight_FRm_normDown);
		tree->Branch("event_weight_FRm_ptUp",&event_weight_FRm_ptUp);
		tree->Branch("event_weight_FRm_ptDown",&event_weight_FRm_ptDown);
		tree->Branch("event_weight_FRm_bUp",&event_weight_FRm_bUp);
		tree->Branch("event_weight_FRm_bDown",&event_weight_FRm_bDown);
		tree->Branch("event_weight_FRm_ecUp",&event_weight_FRm_ecUp);
		tree->Branch("event_weight_FRm_ecDown",&event_weight_FRm_ecDown);
	}
	
	tree->Branch("nJet", &nJet);
	tree->Branch("avg_dr_jet", &avg_dr_jet);
	
	if (anatype_ == Analyze_2lss1tau) {		
		tree->Branch("mindr_lep0_jet", &mindr_lep0_jet);
		tree->Branch("mindr_lep1_jet", &mindr_lep1_jet);		
		tree->Branch("max_lep_eta", &max_lep_eta);
		tree->Branch("met", &met);
		tree->Branch("mht", &mht);
		tree->Branch("mT_met_lep0", &mT_met_lep0);
		tree->Branch("lep0_conept", &lep0_conept);
		tree->Branch("lep1_conept", &lep1_conept);
		tree->Branch("dr_leps", &dr_leps);
		tree->Branch("tau_pt", &tau0_pt);
		tree->Branch("dr_lep0_tau", &dr_lep0_tau);
		tree->Branch("dr_lep1_tau", &dr_lep1_tau);
		tree->Branch("mvis_lep0_tau", &mvis_lep0_tau);
		tree->Branch("mvis_lep1_tau", &mvis_lep1_tau);
		tree->Branch("tau0_decaymode", &tau0_decaymode);
		tree->Branch("tau0_E", &tau0_E);
		tree->Branch("tau0_upsilon", &tau0_upsilon);
	}
	else if (anatype_ == Analyze_1l2tau) {
		tree->Branch("ht", &HT);
		tree->Branch("tt_deltaR", &tt_deltaR);
		tree->Branch("tt_mvis", &tt_mvis);
		tree->Branch("tt_sumpt", &tt_pt);
		tree->Branch("max_dr_jet", &max_dr_jet);
		tree->Branch("tau0_pt", &tau0_pt);
		tree->Branch("tau1_pt", &tau1_pt);
		tree->Branch("ntags", &ntags);
		tree->Branch("ntags_loose", &ntags_loose);
		
		tree->Branch("taup_decaymode", &taup_decaymode);
		tree->Branch("taum_decaymode", &taum_decaymode);
		tree->Branch("taup_E", &taup_E);
		tree->Branch("taum_E", &taum_E);
		tree->Branch("taup_upsilon", &taup_upsilon);
		tree->Branch("taum_upsilon", &taum_upsilon);
	}
	else if (anatype_ == Analyze_3l1tau) {
		tree->Branch("max_lep_eta", &max_lep_eta);
		tree->Branch("mindr_lep0_jet", &mindr_lep0_jet);
		tree->Branch("mindr_lep1_jet", &mindr_lep1_jet);
		tree->Branch("mindr_lep2_jet", &mindr_lep2_jet);
		tree->Branch("mT_met_lep0", &mT_met_lep0);
		tree->Branch("lep0_conept", &lep0_conept);
		tree->Branch("lep1_conept", &lep1_conept);
		tree->Branch("lep2_conept", &lep2_conept);
		tree->Branch("tau0_decaymode", &tau0_decaymode);
		tree->Branch("tau0_E", &tau0_E);
		tree->Branch("tau0_upsilon", &tau0_upsilon);
	}

	tree->Branch("isGenMatchedTau", &isGenMatchedTau);
	tree->Branch("HiggsDecayType", &HiggsDecayType);

	if (evaluate_) {
		tree->Branch("mva_ttbar", &mva_ttbar);
		tree->Branch("mva_ttV", &mva_ttV);
	}
}

#endif
