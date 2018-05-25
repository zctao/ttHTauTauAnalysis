#ifndef SFHelper_h
#define SFHelper_h

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#endif
#include "miniLepton.h"
#include "miniTau.h"
#include "miniJet.h"

#include "Types_enum.h"
// Root
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

#include <iostream>
#include <algorithm>
#include <assert.h>

class SFHelper
{
 public:
	// constructor and destructor
	SFHelper(Analysis_types, Selection_types, bool, bool debug=false);
	~SFHelper();

	// member functions
	float Get_HLTSF(const std::vector<miniLepton>&, const std::vector<miniTau>&,
					bool, bool);
	float Get_HLTSF_2l1tau(int);
	float Get_HLTSF_3l1tau();
	float Get_HLTSF_1l2tau(const miniLepton&, const std::vector<miniTau>&, bool,
						   bool);
	float Get_HLTSF_1l2tau(float, float, int, float, float, int, float, float,
						   int, bool, bool);
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	float Get_HLTSF_1l2tau(const miniLepton&, const std::vector<pat::Tau>&, bool,
						   bool);
#endif
	float Compute_trig_eff_OR_1l2tau(float, float, float, float, bool, bool);
	float Get_trig_eff_lepLeg_crossTrigger(float, float, int, bool);
	float Get_trig_eff_tauLeg_crossTrigger(float, float, int, int, bool);
	float Get_trig_eff_singleLep(float, float, int, bool);
	
	float Get_PUWeight(int);
	
	//float Get_MCWeight();
	
	float Get_FakeRate_lep(float,float,bool,bool,TString syst="NA");
	float Get_FakeRate_lep(const miniLepton&, TString syst="NA");
	float Get_FakeRate_tau(float,float,TString);
	//float Get_FR_weight(float,float,bool,bool,bool,float,float,bool,bool,bool,
	//					TString syst="NA");
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	float Get_FakeRate_tau(const pat::Tau&, TString);
	float Get_FR_weight(const std::vector<miniLepton>&,
						const std::vector<pat::Tau>&, TString syst="NA");
#endif
	float Get_FR_weight(const std::vector<miniLepton>&,
						const std::vector<miniTau>&, TString syst="NA");

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	float Get_EvtCSVWeight(const std::vector<pat::Jet> &, const std::string &);
	float Get_ChargeFlipWeight(const std::vector<miniLepton>&,
							   const std::vector<pat::Tau>&);
#endif
	float Get_EvtCSVWeight(const std::vector<miniJet>&, const std::string&);
	float Get_JetCSVWeight(const miniJet&, std::string);
	
	float Get_ChargeFlipWeight(const std::vector<miniLepton>&, int);
	float Get_EleChargeMisIDProb(float,float,int,int);
	float Get_EleChargeMisIDProb(const miniLepton&, int);

	float Get_LeptonIDSF(float,float,bool,bool,bool);
	float Get_LeptonIDSF(const miniLepton&);
	float Get_LeptonSF_loose(float,float,bool,bool);
	float Get_LeptonSF_tight_vs_loose(float,float,bool,bool);
	float Get_LeptonSF_loose(const miniLepton&);
	float Get_LeptonSF_tight_vs_loose(const miniLepton&);
	float Get_LeptonIDSF_weight(const std::vector<miniLepton>&);

	float Get_TauIDSF(float,float,bool,TString syst="NA");
	float Get_TauIDSF(const miniTau&, TString syst="NA");
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	float Get_TauIDSF(const pat::Tau&, bool, TString syst="NA");
#endif
	float Get_TauIDSF_weight(const std::vector<miniTau>&, TString syst="NA");
		
	// utilities
	float read2DHist(TH2*, float, float);
	float evalTGraph(TGraphAsymmErrors*, float);
	float readTGraph(TGraphAsymmErrors*, float);
	float readTF(TF1*, float);

	inline Analysis_types getAnaType() const {return _analysis;}
	inline Selection_types getSelType() const {return _selection;}
	
 private:

	bool _isdata;
	bool _debug;
	Analysis_types _analysis;
	Selection_types _selection;

	//float _hlt_sf;      // trigger scale factor
	//float _lepIDEff_sf; // lepton ID efficiency scale factor
	//float _csv_weight;  // btag csv re-weighting
	//float _pu_weight;   // pileup weight
	//float _tauID_sf;    // tau ID scale factor
	//float _mc_weight;   // gen weight
	//float _FR_weight;   // fake rate weight;
	//float _QF_weight;   // charge flip weight;

	// Fake lepton rate 
	TFile *file_fr_lep;
	TH2F *h_fakerate_el;
	TH2F *h_fakerate_mu;
	// systematics
	TH2F *h_fakerate_el_normUp;
	TH2F *h_fakerate_el_normDown;
	TH2F *h_fakerate_el_ptUp;
	TH2F *h_fakerate_el_ptDown;
	TH2F *h_fakerate_el_beUp;
	TH2F *h_fakerate_el_beDown;
	//TH2F *h_fakerate_el_bUp;
	//TH2F *h_fakerate_el_bDown;
	//TH2F *h_fakerate_el_ecUp;
	//TH2F *h_fakerate_el_ecDown;
	TH2F *h_fakerate_mu_normUp;
	TH2F *h_fakerate_mu_normDown;
	TH2F *h_fakerate_mu_ptUp;
	TH2F *h_fakerate_mu_ptDown;
	TH2F *h_fakerate_mu_beUp;
	TH2F *h_fakerate_mu_beDown;
	//TH2F *h_fakerate_mu_bUp;
	//TH2F *h_fakerate_mu_bDown;
	//TH2F *h_fakerate_mu_ecUp;
	//TH2F *h_fakerate_mu_ecDown;

	// Jet to tau fake rate
	TFile* file_fr_tau;
	TGraphAsymmErrors *g_fakerate_tau_mvaM_etaL_mc;
	TGraphAsymmErrors *g_fakerate_tau_mvaM_etaH_mc;
	TF1 *f_fakerate_tau_mvaM_etaL_ratio;
	TF1 *f_fakerate_tau_mvaM_etaH_ratio;
	// systematics
	TF1 *f_fakerate_tau_mvaM_etaL_ratio_normUp;
	TF1 *f_fakerate_tau_mvaM_etaL_ratio_normDown;
	TF1 *f_fakerate_tau_mvaM_etaL_ratio_shapeUp;
	TF1 *f_fakerate_tau_mvaM_etaL_ratio_shapeDown;
	TF1 *f_fakerate_tau_mvaM_etaH_ratio_normUp;
	TF1 *f_fakerate_tau_mvaM_etaH_ratio_normDown;
	TF1 *f_fakerate_tau_mvaM_etaH_ratio_shapeUp;
	TF1 *f_fakerate_tau_mvaM_etaH_ratio_shapeDown;

	// PU weight
	TFile* file_puweight;
	TH1F *h_puweight;
	
	// Lepton ID scale factor lookup tables
	TFile* file_recoToLoose_leptonSF_mu1;
	TFile* file_recoToLoose_leptonSF_mu2;
	TFile* file_recoToLoose_leptonSF_mu3;
	TFile* file_recoToLoose_leptonSF_mu4;
	TH2F *h_recoToLoose_leptonSF_mu1;
	TH2F *h_recoToLoose_leptonSF_mu2;
	TH2F *h_recoToLoose_leptonSF_mu3;
	TGraphAsymmErrors *h_recoToLoose_leptonSF_mu4;
	
	TFile* file_recoToLoose_leptonSF_el;
	TFile* file_recoToLoose_leptonSF_gsf;
	TH2F *h_recoToLoose_leptonSF_el1;
	TH2F *h_recoToLoose_leptonSF_el2;
	TH2F *h_recoToLoose_leptonSF_el3;
	TH2F *h_recoToLoose_leptonSF_gsf;
	
	TFile* file_looseToTight_leptonSF_mu_2lss;
	TFile* file_looseToTight_leptonSF_el_2lss;
	TFile* file_looseToTight_leptonSF_mu_3l;
	TFile* file_looseToTight_leptonSF_el_3l;
	TH2F *h_looseToTight_leptonSF_mu_2lss;
	TH2F *h_looseToTight_leptonSF_el_2lss;
	TH2F *h_looseToTight_leptonSF_mu_3l;
	TH2F *h_looseToTight_leptonSF_el_3l;

	// Electron Charge MisID
	TFile *file_eleMisCharge;
	TH2F *h_chargeMisId;

	// trigger scale factors
	// 1l2tau
	TFile *file_Mu22OR_eff;
	TGraphAsymmErrors* g_Mu_ZMassEtaLt0p9_MC;
	TGraphAsymmErrors* g_Mu_ZMassEta0p9to1p2_MC;
	TGraphAsymmErrors* g_Mu_ZMassEta1p2to2p1_MC;
	TGraphAsymmErrors* g_Mu_ZMassEtaLt0p9_Data;
	TGraphAsymmErrors* g_Mu_ZMassEta0p9to1p2_Data;
	TGraphAsymmErrors* g_Mu_ZMassEta1p2to2p1_Data;	

	TFile *file_Ele25WPTight_eff;
	TGraphAsymmErrors* g_Ele_ZMassEtaLt1p48_MC;
	TGraphAsymmErrors* g_Ele_ZMassEta1p48to2p1_MC;
	TGraphAsymmErrors* g_Ele_ZMassEtaGt2p1_MC;
	TGraphAsymmErrors* g_Ele_ZMassEtaLt1p48_Data;
	TGraphAsymmErrors* g_Ele_ZMassEta1p48to2p1_Data;
	TGraphAsymmErrors* g_Ele_ZMassEtaGt2p1_Data;

	TFile *file_Mu19leg_eff;
	TGraphAsymmErrors* g_Muleg_ZMassEtaLt0p9_MC;
	TGraphAsymmErrors* g_Muleg_ZMassEta0p9to1p2_MC;
	TGraphAsymmErrors* g_Muleg_ZMassEta1p2to2p1_MC;
	TGraphAsymmErrors* g_Muleg_ZMassEtaLt0p9_Data;
	TGraphAsymmErrors* g_Muleg_ZMassEta0p9to1p2_Data;
	TGraphAsymmErrors* g_Muleg_ZMassEta1p2to2p1_Data;

	TFile *file_Ele24leg_eff;
	TGraphAsymmErrors* g_Eleleg_ZMassEtaLt1p48_MC;
	TGraphAsymmErrors* g_Eleleg_ZMassEta1p48to2p1_MC;
	TGraphAsymmErrors* g_Eleleg_ZMassEtaGt2p1_MC;
	TGraphAsymmErrors* g_Eleleg_ZMassEtaLt1p48_Data;
	TGraphAsymmErrors* g_Eleleg_ZMassEta1p48to2p1_Data;
	TGraphAsymmErrors* g_Eleleg_ZMassEtaGt2p1_Data;

	TFile *file_mt_Tauleg_eff;
	TGraphAsymmErrors* g_mt_mc_genuine_barrel_TightIso;
	TGraphAsymmErrors* g_mt_mc_genuine_endcap_TightIso;
	TGraphAsymmErrors* g_mt_data_genuine_barrel_TightIso;
	TGraphAsymmErrors* g_mt_data_genuine_endcap_TightIso;
	
	TFile *file_et_Tauleg_eff;
	TGraphAsymmErrors* g_et_mc_genuine_barrel_TightIso;
	TGraphAsymmErrors* g_et_mc_genuine_endcap_TightIso;
	TGraphAsymmErrors* g_et_data_genuine_barrel_TightIso_dm0;
	TGraphAsymmErrors* g_et_data_genuine_barrel_TightIso_dm1;
	TGraphAsymmErrors* g_et_data_genuine_barrel_TightIso_dm10;
	TGraphAsymmErrors* g_et_data_genuine_endcap_TightIso_dm0;
	TGraphAsymmErrors* g_et_data_genuine_endcap_TightIso_dm1;
	TGraphAsymmErrors* g_et_data_genuine_endcap_TightIso_dm10;
	
	// CSV
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	BTagCalibrationReader* BTagCaliReader;
	//std::string Btag_sysList[16] =
	//	{"LFUp","LFDown","HFUp","HFDown",
	//	 "HFStats1Up","HFStats1Down","HFStats2Up","HFStats2Down",
	//	 "LFStats1Up","LFStats1Down","LFStats2Up","LFStats2Down",
	//	 "cErr1Up","cErr1Down","cErr2Up","cErr2Down"};
	void Set_up_BTagCalibration_Readers();
	void Delete_BTagCalibration_Readers();
#endif
	
	void Set_up_FakeRate_Lut(TString tauIDWP="dR03mvaTight");
	void Set_up_TauSF_Lut(TString tauIDWP="dR03mvaMedium");
	void Set_up_ChargeMisID_Lut();
	void Set_up_PUWeight_hist();
	void Set_up_LeptonSF_Lut();
	void Set_up_triggerSF_Lut();

	void Delete_FakeRate_Lut();
	void Delete_TauSF_Lut();
	void Delete_ChargeMisID_Lut();
	void Delete_PUWeight_hist();
	void Delete_LeptonSF_Lut();
	void Delete_triggerSF_Lut();


	//float Get_LeptonSF_loose(float,float,bool,bool);
	//float Get_LeptonSF_tight_vs_loose(float,float,bool,bool);
	//float Get_LeptonSF_loose(const miniLepton&);
	//float Get_LeptonSF_tight_vs_loose(const miniLepton&);
};

#endif
