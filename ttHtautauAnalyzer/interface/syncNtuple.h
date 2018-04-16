#ifndef syncNtuple_h
#define syncNtuple_h

#include "TTree.h"

class syncNtuple
{
 public:

	syncNtuple(){};
	~syncNtuple(){};

	void initialize();
	void set_up_branches(TTree*);

		/// variables
	// event variables
	unsigned long long nEvent;
	int ls;   // luminosity section number
	int run;  // run number

	int n_presel_mu;
	int n_fakeablesel_mu;
	int n_mvasel_mu;
	int n_presel_ele;
	int n_fakeablesel_ele;
	int n_mvasel_ele;
	int n_presel_tau;
	//int n_tau;
	int n_presel_jet;

	// muons
	float mu1_pt;
	float mu1_conept;
	float mu1_eta;
	float mu1_phi;
	float mu1_E;
	int   mu1_charge;
	int   mu1_jetNDauChargedMVASel;
	float mu1_miniRelIso;
	float mu1_miniIsoCharged;
	float mu1_miniIsoNeutral;
	float mu1_jetPtRel;
	float mu1_jetPtRatio;
	float mu1_jetCSV;
	float mu1_sip3D;
	float mu1_dxy;
	float mu1_dxyAbs;
	float mu1_dz;
	float mu1_segmentCompatibility;
	float mu1_leptonMVA;
	float mu1_mediumID;
	float mu1_dpt_div_pt;
	int   mu1_ismvasel;
	int   mu1_isfakeablesel;
	float mu1_PFRelIso04;
	//int   mu1_mcMatchType;
	//int   mu1_isPFMuon;
	float mu2_pt;
	float mu2_conept;
	float mu2_eta;
	float mu2_phi;
	float mu2_E;
	int   mu2_charge;
	int   mu2_jetNDauChargedMVASel;
	float mu2_miniRelIso;
	float mu2_miniIsoCharged;
	float mu2_miniIsoNeutral;
	float mu2_jetPtRel;
	float mu2_jetPtRatio;
	float mu2_jetCSV;
	float mu2_sip3D;
	float mu2_dxy;
	float mu2_dxyAbs;
	float mu2_dz;
	float mu2_segmentCompatibility;
	float mu2_leptonMVA;
	float mu2_mediumID;
	float mu2_dpt_div_pt;
	int   mu2_ismvasel;
	int   mu2_isfakeablesel;
	float mu2_PFRelIso04;
	//int   mu2_mcMatchType;
	//int   mu2_isPFMuon;
	
	// electrons
	float ele1_pt;
	float ele1_conept;
	float ele1_eta;
	float ele1_phi;
	float ele1_E;
	int   ele1_charge;
	int   ele1_jetNDauChargedMVASel;
	float ele1_miniRelIso;
	float ele1_miniIsoCharged;
	float ele1_miniIsoNeutral;
	float ele1_jetPtRel;
	float ele1_jetPtRatio;
	float ele1_jetCSV;
	float ele1_sip3D;
	float ele1_dxy;
	float ele1_dxyAbs;
	float ele1_dz;
	float ele1_ntMVAeleID;
	float ele1_leptonMVA;
	int   ele1_isChargeConsistent;
	int   ele1_passesConversionVeto;
	int   ele1_nMissingHits;
	int   ele1_ismvasel;
	int   ele1_isfakeablesel;
	float ele1_PFRelIso04;
	//int   ele1_mcMatchType;
	float ele2_pt;
	float ele2_conept;
	float ele2_eta;
	float ele2_phi;
	float ele2_E;
	int   ele2_charge;
	int   ele2_jetNDauChargedMVASel;
	float ele2_miniRelIso;
	float ele2_miniIsoCharged;
	float ele2_miniIsoNeutral;
	float ele2_jetPtRel;
	float ele2_jetPtRatio;
	float ele2_jetCSV;
	float ele2_sip3D;
	float ele2_dxy;
	float ele2_dxyAbs;
	float ele2_dz;
	float ele2_ntMVAeleID;
	float ele2_leptonMVA;
	int   ele2_isChargeConsistent;
	int   ele2_passesConversionVeto;
	int   ele2_nMissingHits;
	int   ele2_ismvasel;
	int   ele2_isfakeablesel;
	float ele2_PFRelIso04;
	//int   ele2_mcMatchType;
	
	// taus
	float tau1_pt;
	float tau1_eta;
	float tau1_phi;
	float tau1_E;
	int   tau1_charge;
	float tau1_dxy;
	float tau1_dz;
	//int   tau1_decayMode;
	int   tau1_decayModeFindingOldDMs;
	int   tau1_decayModeFindingNewDMs;
	float tau1_byCombinedIsolationDeltaBetaCorr3Hits;
	int   tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits;
	int   tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits;
	int   tau1_byTightCombinedIsolationDeltaBetaCorr3Hits;
	int   tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03;
	int   tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03;
	int   tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03;
	int   tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
	int   tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
	int   tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT;
	int   tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
	float tau1_rawMVArun2v1DBdR03oldDMwLT;
	int   tau1_againstMuonLoose3;
	int   tau1_againstMuonTight3;
	int   tau1_againstElectronVLooseMVA6;
	int   tau1_againstElectronLooseMVA6;
	int   tau1_againstElectronMediumMVA6;
	int   tau1_againstElectronTightMVA6;
	//int   tau1_mcMatchType;
	float tau2_pt;
	float tau2_eta;
	float tau2_phi;
	float tau2_E;
	int   tau2_charge;
	float tau2_dxy;
	float tau2_dz;
	//int   tau2_decayMode;
	int   tau2_decayModeFindingOldDMs;
	int   tau2_decayModeFindingNewDMs;
	float tau2_byCombinedIsolationDeltaBetaCorr3Hits;
	int   tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits;
	int   tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits;
	int   tau2_byTightCombinedIsolationDeltaBetaCorr3Hits;
	int   tau2_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03;
	int   tau2_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03;
	int   tau2_byTightCombinedIsolationDeltaBetaCorr3HitsdR03;
	int   tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
	int   tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
	int   tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT;
	int   tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT;
	float tau2_rawMVArun2v1DBdR03oldDMwLT;
	int   tau2_againstMuonLoose3;
	int   tau2_againstMuonTight3;
	int   tau2_againstElectronVLooseMVA6;
	int   tau2_againstElectronLooseMVA6;
	int   tau2_againstElectronMediumMVA6;
	int   tau2_againstElectronTightMVA6;
	//int   tau2_mcMatchType;
	
	// jets
	float jet1_pt;
	float jet1_eta;
	float jet1_phi;
	float jet1_E;
	float jet1_CSV;
	float jet2_pt;
	float jet2_eta;
	float jet2_phi;
	float jet2_E;
	float jet2_CSV;
	float jet3_pt;
	float jet3_eta;
	float jet3_phi;
	float jet3_E;
	float jet3_CSV;
	float jet4_pt;
	float jet4_eta;
	float jet4_phi;
	float jet4_E;
	float jet4_CSV;
	
	// MET
	float PFMET;
	float PFMETphi;
	float MHT;
	float metLD;
	//float METSignificance;
	//float METCov00;
	//float METCov10;
	//float METCov01;
	//float METCov11;

	// event-level MVA input variables
	int isGenMatched;	
	float lep1_conept;
	float lep2_conept;
	float lep3_conept;
	float mindr_lep1_jet;
	float mindr_lep2_jet;
	float mindr_lep3_jet;
	float mindr_tau1_jet;
	float mindr_tau2_jet;
	float avg_dr_jet;
	float avg_dr_lep_tau;
	float max_dr_lep_tau;
	float mindr_tau_jet;
	float min_dr_lep_tau;
	float min_dr_lep_jet;
	float dr_leps;
	float dr_taus;
	float dR_lep_tau_ss;
	float dr_lep1_tau;
	float dr_lep2_tau;
	float max_lep_eta;
	float mT_lep1;
	float mT_lep2;
	float mTauTauVis;
	float mTauTauVis1;
	float mTauTauVis2;
	float mbb;
	float mbb_loose;
	float cosThetaS_hadTau;
	float HTT;
	float HadTop_pt;
	float Hj_tagger;
	int nBJetLoose;

	// event-level MVA output variables
	float mvaOutput_plainKin_ttV;
	float mvaOutput_plainKin_ttbar;
	float mvaOutput_1l_2tau_HTT_SUM_VT;
	float mvaOutput_2l_2tau_plainKin_1B_VT;
	float mvaOutput_2l_2tau_plainKin_SUM_VT;
	float mvaOutput_2lss_ttV;
	float mvaOutput_2lss_ttbar;
	float mvaOutput_2lss_1tau_plainKin_ttbar;
	float mvaOutput_2lss_1tau_plainKin_ttV;
	float mvaOutput_2lss_1tau_plainKin_1B_M;
	float mvaOutput_2lss_1tau_plainKin_SUM_M;
	float mvaOutput_2lss_1tau_HTT_SUM_M;
	float mvaOutput_2lss_1tau_HTTMEM_SUM_M;
	float mvaOutput_3l_ttV;
	float mvaOutput_3l_ttbar;
	float mvaOutput_3l_1tau_plainKin_SUM_M;
	float mvaOutput_3l_1tau_plainKin_1B_M;

	// event weights
	//float event_weight;
	float FR_weight;
	float triggerSF_weight;//hltSF;
	float leptonSF_weight;
	float tauSF_weight;
	float bTagSF_weight; //csv_weight;
	float PU_weight;
	float MC_weight;
	
	// MEM variables
	float Integral_ttH;
	float Integral_ttZ;
	float Integral_ttZ_Zll;
	float Integral_ttbar;
	int   Integration_type;
	float memOutput_LR;
   	
	//ClassDef(CU_ttH_EDA_Ntuple,1);
};

#endif
