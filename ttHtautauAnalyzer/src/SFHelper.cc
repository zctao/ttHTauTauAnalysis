#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/SFHelper.h"

// constructor
SFHelper::SFHelper(Analysis_types analysis, Selection_types selection,
				   std::string samplename, bool isdata, bool debug)
{
	_analysis = analysis;
	_selection = selection;
	_isdata = isdata;
	_debug = debug;

	TString tauIDWP = (analysis==Analyze_1l2tau or analysis==Analyze_2l2tau) ?
	    "dR03mvaMedium" : "dR03mvaLoose";

	Set_up_TauFR_Lut(tauIDWP);
	
	if (not _isdata) {
		//Set_up_PUWeight_hist();
		Set_up_PUWeights(samplename);
		Set_up_LeptonSF_Lut();
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
		Set_up_BTagCalibration_Readers();
#endif
		if (_analysis == Analyze_1l2tau)
			Set_up_triggerSF_Lut();
	}

	if (_selection == Application_Fake_2lss1tau or
		_selection == Application_Fake_1l2tau or
		_selection == Application_Fake_3l1tau or
		_selection == Application_Fake_2l2tau or
		_selection == Control_FakeAR_1l2tau or
		_selection == Control_FakeAR_2lss1tau or
		_selection == Control_FakeAR_3l1tau or
		_selection == Control_FakeAR_2l2tau or
		_selection == Control_FakeAR_ttW or
		_selection == Control_FakeAR_ttZ) {
		Set_up_FakeRate_Lut(tauIDWP);
	}

	if (_selection == Application_Flip_2lss1tau or
		_selection == Control_FlipAR_2lss1tau or
		_selection == Control_FlipAR_ttW) {
		Set_up_ChargeMisID_Lut();
	}
}

SFHelper::~SFHelper()
{
	Delete_TauFR_Lut();
	
	if (not _isdata) {	
		//Delete_PUWeight_hist();
		Delete_LeptonSF_Lut();
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
		Delete_BTagCalibration_Readers();
#endif
		if (_analysis == Analyze_1l2tau)
			Delete_triggerSF_Lut();
	}

	if (_selection == Application_Fake_2lss1tau or
		_selection == Application_Fake_1l2tau or
		_selection == Application_Fake_3l1tau or
		_selection == Application_Fake_2l2tau or
		_selection == Control_FakeAR_1l2tau or
		_selection == Control_FakeAR_2lss1tau or
		_selection == Control_FakeAR_3l1tau or
		_selection == Control_FakeAR_2l2tau or
		_selection == Control_FakeAR_ttW or
		_selection == Control_FakeAR_ttZ) {
		Delete_FakeRate_Lut();
	}

	if (_selection == Application_Flip_2lss1tau or
		_selection == Control_FlipAR_2lss1tau or
		_selection == Control_FlipAR_ttW) {
		Delete_ChargeMisID_Lut();
	}
}

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
void SFHelper::Set_up_BTagCalibration_Readers()
{
	const std::string base =
		std::string(getenv("CMSSW_BASE")) +  "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/";
	
	//BTagCalibration calib_csvv2("csvv2", base + "CSVv2_Moriond17_B_H.csv");
	//BTagCalibration calib_csvv2("csvv2", base + "CSVv2_94XSF_V1_B_F.csv");
	BTagCalibration calib_deepcsv("deepcsv", base + "DeepCSV_94XSF_V1_B_F.csv");
	
	BTagCaliReader = new BTagCalibrationReader(
	     BTagEntry::OP_RESHAPING, // operating point
		 "central",
		 {"up_jes", "down_jes", "up_lf", "down_lf", "up_hf", "down_hf",
				 "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2",
				 "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2",
				 "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"}
											   );
	
	BTagCaliReader->load(calib_deepcsv,BTagEntry::FLAV_B,"iterativefit");
	BTagCaliReader->load(calib_deepcsv,BTagEntry::FLAV_C,"iterativefit");
	BTagCaliReader->load(calib_deepcsv,BTagEntry::FLAV_UDSG,"iterativefit");
}
#endif

void SFHelper::Set_up_TauFR_Lut(TString tauIDWP/*"dR03mvaMedium"*/)
{
	//taus fake rate
	file_fr_tau = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/FR_tau_2017_v1.root").c_str(), "read");
	
	g_fakerate_tau_mvaM_etaL_mc = (TGraphAsymmErrors*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEtaLt1_5/jetToTauFakeRate_mc_hadTaus_pt");
	g_fakerate_tau_mvaM_etaH_mc = (TGraphAsymmErrors*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEta1_5to9_9/jetToTauFakeRate_mc_hadTaus_pt");
	f_fakerate_tau_mvaM_etaL_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEtaLt1_5/fitFunction_data_div_mc_hadTaus_pt");
	f_fakerate_tau_mvaM_etaH_ratio = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEta1_5to9_9/fitFunction_data_div_mc_hadTaus_pt");

	// systematics
	f_fakerate_tau_mvaM_etaL_ratio_normUp = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEtaLt1_5/fitFunction_data_div_mc_hadTaus_pt_par1Up");
	f_fakerate_tau_mvaM_etaL_ratio_normDown = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEtaLt1_5/fitFunction_data_div_mc_hadTaus_pt_par1Down");
	f_fakerate_tau_mvaM_etaL_ratio_shapeUp = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEtaLt1_5/fitFunction_data_div_mc_hadTaus_pt_par2Up");
	f_fakerate_tau_mvaM_etaL_ratio_shapeDown = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEtaLt1_5/fitFunction_data_div_mc_hadTaus_pt_par2Down");
	f_fakerate_tau_mvaM_etaH_ratio_normUp = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEta1_5to9_9/fitFunction_data_div_mc_hadTaus_pt_par1Up");
	f_fakerate_tau_mvaM_etaH_ratio_normDown = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEta1_5to9_9/fitFunction_data_div_mc_hadTaus_pt_par1Down");
	f_fakerate_tau_mvaM_etaH_ratio_shapeUp = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEta1_5to9_9/fitFunction_data_div_mc_hadTaus_pt_par2Up");
	f_fakerate_tau_mvaM_etaH_ratio_shapeDown = (TF1*) file_fr_tau->Get("jetToTauFakeRate/"+tauIDWP+"/absEta1_5to9_9/fitFunction_data_div_mc_hadTaus_pt_par2Down");
}

/*
void SFHelper::Set_up_PUWeight_hist()
{
	//file_puweight = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/data/PU_weights/PU_weights_2016_271036_284044.root").c_str(), "read");
	file_puweight = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/pu_weight/PU_weights_MCSummer2016_271036_284044.root").c_str(), "read");
	
	h_puweight = (TH1F*) file_puweight->Get("h_ratio_data_MC");
}
*/

const std::map<std::string, std::string> SFHelper::PUSampleNameMap_ = {
	{"DYJets_M10to50","DYJetsToLL_M-10to50"},
	{"DYJets_M50","DYJetsToLL_M-50"},
	{"DYJets_M50_ext","DYJetsToLL_M-50_ext1"},
	{"ggHZZ4l","GluGluHToZZTo4L"},
	{"ST_sLep","ST_s-channel_4f_leptonDecays"},
	{"ST_sLep_psw","ST_s-channel_4f_leptonDecays"}, // FIXME
	{"ST_tTbar","ST_t-channel_antitop_4f_inclusiveDecays"},
	{"ST_tT","ST_t-channel_top_4f_inclusiveDecays"},
	{"ST_tWTbar","ST_tW_antitop_5f_inclusiveDecays"},
	{"ST_tWT","ST_tW_top_5f_inclusiveDecays"},
	{"ST_tWll","ST_tWll"},
	{"tHq","THQ"},
	{"tHW","THW"},
	{"TTGJets","TTGJets"},
	{"ttHJetTobb","ttHJetTobb_M125_amcatnlo"},
	{"ttHJetToNonbb","ttHJetToNonbb_M125_amcatnlo"},
	{"ttHToNonbb","ttHToNonbb_M125_powheg"},
	{"TTToDiLep","TTTo2L2Nu"},
	{"TTToDiLep_psw","TTTo2L2Nu_PSweights"},
	{"TTToHad","TTToHadronic"},
	{"TTToHad_psw","TTToHadronic_PSweights"},
	{"TTToSemiLep","TTToSemiLeptonic"},
	{"TTToSemiLep_psw","TTToSemiLeptonic_PSweights"},
	{"TTTT","TTTT"},
	{"TTW","TTWJetsToLNu"},
	{"TTW_psw","TTWJetsToLNu_PSweights"},
	{"TTWW","TTWW"},
	{"TTZ_M1to10","TTZToLL_M-1to10"},
	{"TTZ","TTZToLL_M10"},
	{"tZq","tZq_ll_4f"},
	{"VHToNonbb","VHToNonbb_M125"},
	{"WJets","WJetsToLNu"},
	{"WpWpJJ","WpWpJJ_EWK_QCD"},
	{"WW","WWTo2L2Nu"},
	//{"WWds","WWTo2L2Nu_DoubleScattering"}, // FIXME
	{"WWTo2L2Nuds","WWTo2L2Nu_DoubleScattering"},
	{"WWW","WWW_4F"},
	{"WWZ","WWZ_4F"},
	{"WZ","WZTo3LNu"}, // FIXME
	{"WZTo3LNu","WZTo3LNu"},
	{"WZZ","WZZ"},
	{"ZZ","ZZTo4L"},
	{"ZZ_ext","ZZTo4L_ext1"},
	{"ZZZ","ZZZ"}
};		

void SFHelper::Set_up_PUWeights(const std::string& samplename)
{
	std::string fname_pu_data = std::string(getenv("CMSSW_BASE"))+"/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/pu_weight/PileupData_ReRecoJSON_Full2017.root";
    std::string fname_pu_mc = std::string(getenv("CMSSW_BASE"))+"/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/pu_weight/pileup_2017.root";

	if (not PUSampleNameMap_.count(samplename)) {
		std::cerr << "sample name : " << samplename << " not exsit in the PU name map"
				  << std::endl;
		assert(0);
	}
	
	LumiWeights_ = edm::LumiReWeighting(fname_pu_mc, fname_pu_data,
										PUSampleNameMap_.at(samplename), "pileup");
}

void SFHelper::Set_up_FakeRate_Lut(TString tauIDWP/*"dR03mvaTight"*/)
{
	// electrons and muons
	//file_fr_lep = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/FR_lep_ttH_mva090_2017_CERN_2018May29.root").c_str(),"read");
	file_fr_lep = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/FR_data_ttH_mva.root").c_str(),"read");
	assert(file_fr_lep->IsOpen());
	
	//h_fakerate_el = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb");
	h_fakerate_el = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC");
	h_fakerate_mu = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb");

	// systematics
	h_fakerate_el_normUp = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_up");
	h_fakerate_el_normDown = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_down");
	h_fakerate_el_ptUp = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_pt1");
	h_fakerate_el_ptDown = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_pt2");
	h_fakerate_el_beUp = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_be1");
	h_fakerate_el_beDown = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_be2");
	//h_fakerate_el_bUp = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_b1");
	//h_fakerate_el_bDown = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_b2");
	//h_fakerate_el_ecUp = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_ec1");
	//h_fakerate_el_ecDown = (TH2F*) file_fr_lep->Get("FR_mva090_el_data_comb_NC_ec2");
	h_fakerate_mu_normUp = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_up");
	h_fakerate_mu_normDown = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_down");
	h_fakerate_mu_ptUp = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_pt1");
	h_fakerate_mu_ptDown = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_pt2");
	h_fakerate_mu_beUp = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_be1");
	h_fakerate_mu_beDown = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_be2");
	//h_fakerate_mu_bUp = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_b1");
	//h_fakerate_mu_bDown = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_b2");
	//h_fakerate_mu_ecUp = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_ec1");
	//h_fakerate_mu_ecDown = (TH2F*) file_fr_lep->Get("FR_mva090_mu_data_comb_ec2");

	h_fakerate_el_TT = (TH2F*) file_fr_lep->Get("FR_mva090_el_TT");
	h_fakerate_el_QCD = (TH2F*) file_fr_lep->Get("FR_mva090_el_QCD_NC");
	h_fakerate_mu_TT = (TH2F*) file_fr_lep->Get("FR_mva090_mu_TT");
	h_fakerate_mu_QCD = (TH2F*) file_fr_lep->Get("FR_mva090_mu_QCD");
}

void SFHelper::Set_up_ChargeMisID_Lut()
{
	file_eleMisCharge = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/ElectronChargeMisIdRates_2017.root").c_str(), "read");
	h_chargeMisId = (TH2D*) file_eleMisCharge->Get("eChargeMisIdRates");
}

void SFHelper::Set_up_LeptonSF_Lut()
{
	//// loose vs reco
	/// muons
	// tracking efficiency SF (pT>10)
	file_recoToLoose_leptonSF_mu1 = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/trk_Eff.root").c_str(),"read");
	g_recoToLoose_leptonSF_mu1 = (TGraphAsymmErrors*)(file_recoToLoose_leptonSF_mu1->Get("ratio_eff_eta3_dr030e030_corr"));

	// loose ID
	// pt>30
	file_recoToLoose_leptonSF_mu2 = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/muon_RunBCDEF_SF_ID.root").c_str(),"read");
	h_recoToLoose_leptonSF_mu2 = (TH2F*)(file_recoToLoose_leptonSF_mu2->Get("NUM_LooseID_DEN_genTracks_pt_abseta"));
	// pt< 30
	file_recoToLoose_leptonSF_mu3 = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/muon_RunBCDEF_SF_ID_JPsi.root").c_str(),"read");
	h_recoToLoose_leptonSF_mu3 = (TH2F*)(file_recoToLoose_leptonSF_mu3->Get("NUM_LooseID_DEN_genTracks_pt_abseta"));
	
	// muon vertex and isolation cut
	file_recoToLoose_leptonSF_mu4 = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/scaleFactors.root").c_str(),"read");
	h_recoToLoose_leptonSF_mu4 = (TH2D*)(file_recoToLoose_leptonSF_mu4->Get("NUM_ttHLoo_DEN_LooseID"));
	
	/// for electrons
	// reco pt > 20
	file_recoToLoose_leptonSF_gsf = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root").c_str(),"read");
	h_recoToLoose_leptonSF_gsf = (TH2F*)(file_recoToLoose_leptonSF_gsf->Get("EGamma_SF2D"));
	
	// reco pt < 20
	file_recoToLoose_leptonSF_gsf_low = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root").c_str(),"read");
	h_recoToLoose_leptonSF_gsf_low = (TH2F*)(file_recoToLoose_leptonSF_gsf_low->Get("EGamma_SF2D"));

	/*
	// out of date
	// loose vs reco
	file_recoToLoose_leptonSF_el = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/el_scaleFactors_Moriond17.root").c_str(),"read");
	h_recoToLoose_leptonSF_el1 = (TH2F*)(file_recoToLoose_leptonSF_el->Get("GsfElectronToMVAVLooseFOIDEmuTightIP2D"));
	h_recoToLoose_leptonSF_el2 = (TH2F*)(file_recoToLoose_leptonSF_el->Get("MVAVLooseElectronToMini4"));
	h_recoToLoose_leptonSF_el3 = (TH2F*)(file_recoToLoose_leptonSF_el->Get("MVAVLooseElectronToConvVetoIHit1"));
	*/
	///////////////////////////////////////////////////
	
	//// tight vs loose
	if (_analysis == Analyze_2lss1tau or _analysis == Analyze_2lss) {
		// tight charge
		/// for muon
		file_looseToTight_leptonSF_mu_2lss = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/lepMVAEffSF_m_2lss.root").c_str(),"read");
		h_looseToTight_leptonSF_mu_2lss = (TH2F*)(file_looseToTight_leptonSF_mu_2lss->Get("sf"));
		
		/// for electron
		file_looseToTight_leptonSF_el_2lss = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/lepMVAEffSF_e_2lss.root").c_str(),"read");
		h_looseToTight_leptonSF_el_2lss = (TH2F*)(file_looseToTight_leptonSF_el_2lss->Get("sf"));
	}
	
	if (_analysis == Analyze_3l1tau or _analysis == Analyze_3l or
		_analysis == Analyze_2l2tau or _analysis == Analyze_1l2tau) {
		// no tight charge
		/// for muon
		file_looseToTight_leptonSF_mu_3l = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/lepMVAEffSF_m_3l.root").c_str(),"read");
		h_looseToTight_leptonSF_mu_3l = (TH2F*)(file_looseToTight_leptonSF_mu_3l->Get("sf"));
		
		/// for electron
		file_looseToTight_leptonSF_el_3l = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/lepton_sf/lepMVAEffSF_e_3l.root").c_str(),"read");
		h_looseToTight_leptonSF_el_3l = (TH2F*)(file_looseToTight_leptonSF_el_3l->Get("sf"));
	}
}

void SFHelper::Set_up_triggerSF_Lut()
{
	assert(_analysis == Analyze_1l2tau);

	// Single Muon triggers
	file_Mu_SingleMu_hlt_eff = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/trigger_sf/Muon_IsoMu24orIsoMu27_eff.root").c_str(),"read");
	assert(file_Mu_SingleMu_hlt_eff->IsOpen());
	g_Mu_ZMassEtaLt0p9_MC = (TGraphAsymmErrors*)file_Mu_SingleMu_hlt_eff->Get("ZMassEtaLt0p9_MC");
	g_Mu_ZMassEta0p9to1p2_MC = (TGraphAsymmErrors*)file_Mu_SingleMu_hlt_eff->Get("ZMassEta0p9to1p2_MC");
	g_Mu_ZMassEta1p2to2p1_MC = (TGraphAsymmErrors*)file_Mu_SingleMu_hlt_eff->Get("ZMassEta1p2to2p1_MC");
	g_Mu_ZMassEtaGt2p1_MC = (TGraphAsymmErrors*)file_Mu_SingleMu_hlt_eff->Get("ZMassEtaGt2p1_MC");
	g_Mu_ZMassEtaLt0p9_Data = (TGraphAsymmErrors*)file_Mu_SingleMu_hlt_eff->Get("ZMassEtaLt0p9_Data");
	g_Mu_ZMassEta0p9to1p2_Data = (TGraphAsymmErrors*)file_Mu_SingleMu_hlt_eff->Get("ZMassEta0p9to1p2_Data");
	g_Mu_ZMassEta1p2to2p1_Data = (TGraphAsymmErrors*)file_Mu_SingleMu_hlt_eff->Get("ZMassEta1p2to2p1_Data");
	g_Mu_ZMassEtaGt2p1_Data = (TGraphAsymmErrors*)file_Mu_SingleMu_hlt_eff->Get("ZMassEtaGt2p1_Data");

	// Single Eletron triggers
	file_Ele_SingleElE_hlt_eff = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/trigger_sf/Electron_Ele32orEle35_eff.root").c_str(),"read");
	assert(file_Ele_SingleElE_hlt_eff->IsOpen());
	g_Ele_ZMassEtaLt1p48_MC = (TGraphAsymmErrors*)file_Ele_SingleElE_hlt_eff->Get("ZMassEtaLt1p48_MC");
	g_Ele_ZMassEta1p48to2p1_MC = (TGraphAsymmErrors*)file_Ele_SingleElE_hlt_eff->Get("ZMassEta1p48to2p1_MC");
	g_Ele_ZMassEtaLt1p48_Data = (TGraphAsymmErrors*)file_Ele_SingleElE_hlt_eff->Get("ZMassEtaLt1p48_Data");
	g_Ele_ZMassEta1p48to2p1_Data = (TGraphAsymmErrors*)file_Ele_SingleElE_hlt_eff->Get("ZMassEta1p48to2p1_Data");

	// Muon leg of MuTau cross triggers
	file_Mu_MuTau_hlt_eff = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/trigger_sf/Muon_MuTau_IsoMu20_eff.root").c_str(),"read");
	assert(file_Mu_MuTau_hlt_eff->IsOpen());
	g_Muleg_ZMassEtaLt0p9_MC = (TGraphAsymmErrors*)file_Mu_MuTau_hlt_eff->Get("ZMassEtaLt0p9_MC");
	g_Muleg_ZMassEta0p9to1p2_MC = (TGraphAsymmErrors*)file_Mu_MuTau_hlt_eff->Get("ZMassEta0p9to1p2_MC");
	g_Muleg_ZMassEta1p2to2p1_MC = (TGraphAsymmErrors*)file_Mu_MuTau_hlt_eff->Get("ZMassEta1p2to2p1_MC");
	g_Muleg_ZMassEtaGt2p1_MC = (TGraphAsymmErrors*)file_Mu_MuTau_hlt_eff->Get("ZMassEtaGt2p1_MC");
	g_Muleg_ZMassEtaLt0p9_Data = (TGraphAsymmErrors*)file_Mu_MuTau_hlt_eff->Get("ZMassEtaLt0p9_Data");
	g_Muleg_ZMassEta0p9to1p2_Data = (TGraphAsymmErrors*)file_Mu_MuTau_hlt_eff->Get("ZMassEta0p9to1p2_Data");
	g_Muleg_ZMassEta1p2to2p1_Data = (TGraphAsymmErrors*)file_Mu_MuTau_hlt_eff->Get("ZMassEta1p2to2p1_Data");
	g_Muleg_ZMassEtaGt2p1_Data = (TGraphAsymmErrors*)file_Mu_MuTau_hlt_eff->Get("ZMassEtaGt2p1_Data");

	// Electron leg of EleTau cross triggers
	file_Ele_EleTau_hlt_eff = new TFile((std::string(getenv("CMSSW_BASE")) + "/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/trigger_sf/Electron_EleTau_Ele24_eff.root").c_str(),"read");
	assert(file_Ele_EleTau_hlt_eff->IsOpen());
	g_Eleleg_ZMassEtaLt1p48_MC = (TGraphAsymmErrors*)file_Ele_EleTau_hlt_eff->Get("ZMassEtaLt1p48_MC");
	g_Eleleg_ZMassEta1p48to2p1_MC = (TGraphAsymmErrors*)file_Ele_EleTau_hlt_eff->Get("ZMassEta1p48to2p1_MC");
	g_Eleleg_ZMassEtaLt1p48_Data = (TGraphAsymmErrors*)file_Ele_EleTau_hlt_eff->Get("ZMassEtaLt1p48_Data");
	g_Eleleg_ZMassEta1p48to2p1_Data = (TGraphAsymmErrors*)file_Ele_EleTau_hlt_eff->Get("ZMassEta1p48to2p1_Data");

	// tau leg of lepton-tau cross triggers
	tauTrigSFhelper = new TauTriggerSFs2017(std::string(getenv("CMSSW_BASE")) + "/src/TriggerSF/TauTriggerSFs2017/data/tauTriggerEfficiencies2017.root","medium");
	
}

void SFHelper::Delete_FakeRate_Lut()
{
	file_fr_lep->Close();	
	delete file_fr_lep;
}

void SFHelper::Delete_TauFR_Lut()
{
	file_fr_tau->Close();		
	delete file_fr_tau;
}

void SFHelper::Delete_ChargeMisID_Lut()
{
	file_eleMisCharge->Close();
	delete file_eleMisCharge;
}

void SFHelper::Delete_PUWeight_hist()
{
	//file_puweight->Close();
	//delete file_puweight;
	
	//file_pu_data->Close();
	//file_pu_mc->Close();
	//delete file_pu_data;
	//delete file_pu_mc;
}

void SFHelper::Delete_LeptonSF_Lut()
{
	file_recoToLoose_leptonSF_mu1->Close();
	file_recoToLoose_leptonSF_mu2->Close();
	file_recoToLoose_leptonSF_mu3->Close();
	file_recoToLoose_leptonSF_mu4->Close();
	
	//file_recoToLoose_leptonSF_el->Close();
	file_recoToLoose_leptonSF_gsf->Close();
	file_recoToLoose_leptonSF_gsf_low->Close();
	
	delete file_recoToLoose_leptonSF_mu1;
	delete file_recoToLoose_leptonSF_mu2;
	delete file_recoToLoose_leptonSF_mu3;
	delete file_recoToLoose_leptonSF_mu4;
	
	//delete file_recoToLoose_leptonSF_el;
	delete file_recoToLoose_leptonSF_gsf;
	delete file_recoToLoose_leptonSF_gsf_low;
	
	if (_analysis == Analyze_2lss1tau or _analysis == Analyze_2lss) {
		file_looseToTight_leptonSF_mu_2lss->Close();
		file_looseToTight_leptonSF_el_2lss->Close();
		
		delete file_looseToTight_leptonSF_mu_2lss;
		delete file_looseToTight_leptonSF_el_2lss;
	}
	if (_analysis == Analyze_3l1tau or _analysis == Analyze_3l or
		_analysis == Analyze_2l2tau or _analysis == Analyze_1l2tau) {
		file_looseToTight_leptonSF_mu_3l->Close();
		file_looseToTight_leptonSF_el_3l->Close();
		
		delete file_looseToTight_leptonSF_mu_3l;
		delete file_looseToTight_leptonSF_el_3l;
	}
}

void SFHelper::Delete_triggerSF_Lut()
{
	file_Mu_SingleMu_hlt_eff->Close();
	file_Ele_SingleElE_hlt_eff->Close();
	file_Mu_MuTau_hlt_eff->Close();
	file_Ele_EleTau_hlt_eff->Close();
	
	delete file_Mu_SingleMu_hlt_eff;
	delete file_Ele_SingleElE_hlt_eff;
	delete file_Mu_MuTau_hlt_eff;
	delete file_Ele_EleTau_hlt_eff;

	delete tauTrigSFhelper;
}

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
void SFHelper::Delete_BTagCalibration_Readers()
{
	delete BTagCaliReader;
}
#endif

float SFHelper::Get_HLTSF(const std::vector<miniLepton>& leptons,
						 const std::vector<miniTau>& taus,
						 bool hlt1LTriggered, bool hltXTriggered)
{
	if (_analysis==Analyze_2lss1tau or _analysis==Analyze_2lss or
		_analysis==Analyze_2l2tau) {
		assert(leptons.size()>1);
		return Get_HLTSF_2l(leptons);
	}
	else if (_analysis==Analyze_1l2tau) {
		assert(leptons.size()>0);
		return Get_HLTSF_1l2tau(leptons[0], taus, hlt1LTriggered, hltXTriggered);
	}
	else if (_analysis==Analyze_3l1tau or _analysis==Analyze_3l) {
		return Get_HLTSF_3l();
	}
	else {
		std::cout << "SFHelper::Get_HLTSF() : WARNING! analysis type not supported."
				  << std::endl;
		return 1.;
	}
}

float SFHelper::Get_HLTSF_2l(const std::vector<miniLepton>& leptons, TString syst)
{
	assert(_analysis==Analyze_2lss1tau or _analysis==Analyze_2lss or
		   _analysis==Analyze_2l2tau);
	assert(leptons.size()>1);
	// assume leptons are already sorted in descending order by pT

	float hltsf = 1.;
	
	// function of leading lepton pt only for now
	if (abs(leptons[0].pdgId())==13 and abs(leptons[1].pdgId())==13) { // 2mu
		if (leptons[0].conept() < 35.) { // pT < 35
			hltsf = 0.972;  // +/- 0.006
			if (syst=="triggerUp") hltsf += 0.006;
			if (syst=="triggerDown") hltsf -= 0.006;
		}
		else { // pT >= 35
			hltsf = 0.994;  // +/- 0.001
			if (syst=="triggerUp") hltsf += 0.001;
			if (syst=="triggerDown") hltsf -= 0.001;
		}
	}
	else if (abs(leptons[0].pdgId())==11 and abs(leptons[1].pdgId())==11) { // 2ele
		if (leptons[0].conept() < 30.) { // pT < 30
			hltsf = 0.937;  // +/- 0.027
			if (syst=="triggerUp") hltsf += 0.027;
			if (syst=="triggerDown") hltsf -= 0.027;
		}
		else { // pT >= 30
			hltsf = 0.991;  // +/- 0.002
			if (syst=="triggerUp") hltsf += 0.002;
			if (syst=="triggerDown") hltsf -= 0.002;
		}
	}
	else { // ele+mu
		if (leptons[0].conept() < 35.) {
			hltsf = 0.952;  // +/- 0.008
			if (syst=="triggerUp") hltsf += 0.008;
			if (syst=="triggerDown") hltsf -= 0.008;
		}
		else if (leptons[0].conept() < 50.) {
			hltsf = 0.983;  // +/- 0.003
			if (syst=="triggerUp") hltsf += 0.003;
			if (syst=="triggerDown") hltsf -= 0.003;
		}
		else {
			hltsf = 1.0;  // +/- 0.001
			if (syst=="triggerUp") hltsf += 0.001;
			if (syst=="triggerDown") hltsf -= 0.001;
		}
	}

	return hltsf;
}

float SFHelper::Get_HLTSF_3l()
{
	assert(_analysis==Analyze_3l1tau or _analysis==Analyze_3l);

	return 1.;  // +/- 0.05
}

float SFHelper::Get_HLTSF_1l2tau(const miniLepton& lepton,
								 const std::vector<miniTau>& taus,
								 bool passHLT1l, bool passHLT1l1tau)
{
	assert(_analysis==Analyze_1l2tau);
	assert(taus.size()>1);

	// trigger efficiency for single lepton trigger
	float eff_L_data = Get_trig_eff_singleLep(lepton.pt(), lepton.eta(),
											  lepton.pdgId(), true);
	float eff_L_mc = Get_trig_eff_singleLep(lepton.pt(), lepton.eta(),
											lepton.pdgId(), false);
	
 	// trigger efficiency for lepton leg of lep tau cross trigger
	float eff_l_data = Get_trig_eff_lepLeg_crossTrigger(lepton.pt(), lepton.eta(),
														lepton.pdgId(), true);
	float eff_l_mc = Get_trig_eff_lepLeg_crossTrigger(lepton.pt(), lepton.eta(),
													  lepton.pdgId(), false);

	// trigger efficiency for tau leg of lep tau cross trigger
	float eff_t0_data, eff_t0_mc, eff_t1_data, eff_t1_mc;
	if (abs(lepton.pdgId())==11) { // ele+tau
		eff_t0_data = tauTrigSFhelper->getETauEfficiencyData(taus[0].pt(),taus[0].eta(),taus[0].phi());
		eff_t0_mc = tauTrigSFhelper->getETauEfficiencyMC(taus[0].pt(),taus[0].eta(),taus[0].phi());
		eff_t1_data = tauTrigSFhelper->getETauEfficiencyData(taus[1].pt(),taus[1].eta(),taus[1].phi());
		eff_t1_mc = tauTrigSFhelper->getETauEfficiencyMC(taus[1].pt(),taus[1].eta(),taus[1].phi());
	}
	else { // mu+tau
		eff_t0_data = tauTrigSFhelper->getMuTauEfficiencyData(taus[0].pt(),taus[0].eta(),taus[0].phi());
		eff_t0_mc = tauTrigSFhelper->getMuTauEfficiencyMC(taus[0].pt(),taus[0].eta(),taus[0].phi());
		eff_t1_data = tauTrigSFhelper->getMuTauEfficiencyData(taus[1].pt(),taus[1].eta(),taus[1].phi());
		eff_t1_mc = tauTrigSFhelper->getMuTauEfficiencyMC(taus[1].pt(),taus[1].eta(),taus[1].phi());
	}

	// compute trigger efficiency based on triggered HLT paths
	float eff_data = Compute_trig_eff_OR_1l2tau(eff_L_data,eff_l_data,eff_t0_data,
												eff_t1_data, passHLT1l,
												passHLT1l1tau);
	float eff_mc = Compute_trig_eff_OR_1l2tau(eff_L_mc, eff_l_mc, eff_t0_mc,
											  eff_t1_mc, passHLT1l,passHLT1l1tau);
	
	// scale factor
	return std::min( eff_data/std::max(1.e-6f, eff_mc), 1.e+1f);

}

/*
float SFHelper::Get_HLTSF_1l2tau(float lepPt, float lepEta, int lepPdgId,
								 float tau0Pt, float tau0Eta, int decaymode0,
								 float tau1Pt, float tau1Eta, int decaymode1,
								 bool passHLT1l, bool passHLT1l1tau)
{
	
	assert(abs(lepPdgId)==11 or abs(lepPdgId)==13);
	
	// trigger efficiency for tau leg of lep tau cross trigger
	float eff_t0_data =
		Get_trig_eff_tauLeg_crossTrigger(tau0Pt,tau0Eta, decaymode0, lepPdgId, 1);
	float eff_t0_mc =
		Get_trig_eff_tauLeg_crossTrigger(tau0Pt,tau0Eta, decaymode0, lepPdgId, 0);
	float eff_t1_data =
		Get_trig_eff_tauLeg_crossTrigger(tau1Pt,tau1Eta, decaymode1, lepPdgId, 1);
	float eff_t1_mc =
		Get_trig_eff_tauLeg_crossTrigger(tau1Pt,tau1Eta, decaymode1, lepPdgId, 0);	
}
*/

float SFHelper::Compute_trig_eff_OR_1l2tau(float eff_L, float eff_l,
										   float eff_tau0, float eff_tau1,
										   bool passHLT1l, bool passHLT1l1tau)
{
	if (passHLT1l and !passHLT1l1tau) { // only the single lepton trigger fires
		return std::max(1.e-2f,
				   eff_L-std::min(eff_L,eff_l)*(1-(1-eff_tau0)*(1-eff_tau1)));
	}
	else if (!passHLT1l and passHLT1l1tau) { // only cross trigger fires
		return std::max(1.e-2f, (eff_l-eff_L)*(1-(1-eff_tau0)*(1-eff_tau1)));
	}
	else if (passHLT1l and passHLT1l1tau) { // both triggers fire
		return std::min(eff_L, eff_l)*(1-(1-eff_tau0)*(1-eff_tau1));
	}

	return 0.;
}

float SFHelper::Get_trig_eff_lepLeg_crossTrigger(float pt, float eta, int pdgid,
												 bool isdata)
{
	assert(_analysis==Analyze_1l2tau);

	float eff = 0.;
	
	if (abs(pdgid)==11) {  // Electron
		if (fabs(eta)<1.48) {
			if (isdata)
				eff = readTGraph(g_Eleleg_ZMassEtaLt1p48_Data, pt);
			else
				eff = readTGraph(g_Eleleg_ZMassEtaLt1p48_MC, pt);
		}
		else if (fabs(eta)<2.1) {
			if (isdata)
				eff = readTGraph(g_Eleleg_ZMassEta1p48to2p1_Data, pt);
			else
				eff = readTGraph(g_Eleleg_ZMassEta1p48to2p1_MC, pt);
		}
		else {
			// shouldn't be here
			assert(0);
			//if (isdata)
			//	eff = readTGraph(g_Eleleg_ZMassEtaGt2p1_Data, pt);
			//else
			//	eff = readTGraph(g_Eleleg_ZMassEtaGt2p1_MC, pt);
		}
	}
	else if (abs(pdgid)==13) {  // Muon
		if (fabs(eta)<0.9) {
			if (isdata)
				eff = readTGraph(g_Muleg_ZMassEtaLt0p9_Data, pt);
			else
				eff = readTGraph(g_Muleg_ZMassEtaLt0p9_MC, pt);
		}
		else if (fabs(eta)<1.2) {
			if (isdata)
				eff = readTGraph(g_Muleg_ZMassEta0p9to1p2_Data, pt);
			else
				eff = readTGraph(g_Muleg_ZMassEta0p9to1p2_MC, pt);
		}
		else if (fabs(eta)<2.1) {
			if (isdata)
				eff = readTGraph(g_Muleg_ZMassEta1p2to2p1_Data, pt);
			else
				eff = readTGraph(g_Muleg_ZMassEta1p2to2p1_MC, pt);
		}
		else // shouldn't be here
			assert(0);
	}
	else {
		std::cout << "Get_trig_eff_lepLeg_crossTrigger: Oops..." << std::endl;
		assert(0);
	}

	return eff;
}
/*
float SFHelper::Get_trig_eff_tauLeg_crossTrigger(float pt, float eta,
												 int decaymode, int lepPdgid,
												 bool isdata)
{
	assert(_analysis==Analyze_1l2tau);

	float eff = 0.;

	if (abs(lepPdgid)==11) {  // e tau cross trigger
		if (fabs(eta)<1.479) {  // barrel
			if (isdata) {
				if (decaymode==0)
					eff = evalTGraph(g_et_data_genuine_barrel_TightIso_dm0, pt);
				else if (decaymode==1)
					eff = evalTGraph(g_et_data_genuine_barrel_TightIso_dm1, pt);
				else if (decaymode==10)
					eff = evalTGraph(g_et_data_genuine_barrel_TightIso_dm10, pt);
			}
			else
				eff = evalTGraph(g_et_mc_genuine_barrel_TightIso, pt);
		}
		else {  // endcap
			if (isdata) {
				if (decaymode==0)
					eff = evalTGraph(g_et_data_genuine_endcap_TightIso_dm0, pt);
				else if (decaymode==1)
					eff = evalTGraph(g_et_data_genuine_endcap_TightIso_dm1, pt);
				else if (decaymode==10)
					eff = evalTGraph(g_et_data_genuine_endcap_TightIso_dm10, pt);
			}
			else
				eff = evalTGraph(g_et_mc_genuine_endcap_TightIso, pt);

		}
	}
	else if (abs(lepPdgid)==13) {  // mu tau cross trigger
		if (fabs(eta)<1.479) {  // barrel
			if (isdata)
				eff = evalTGraph(g_mt_data_genuine_barrel_TightIso, pt);
			else
				eff = evalTGraph(g_mt_mc_genuine_barrel_TightIso, pt);
		}
		else {  // endcap
			if (isdata)
				eff = evalTGraph(g_mt_data_genuine_endcap_TightIso, pt);
			else
				eff = evalTGraph(g_mt_mc_genuine_endcap_TightIso, pt);
		}
	}
	else
		assert(0);
	
	return eff;
}
*/
float SFHelper::Get_trig_eff_singleLep(float pt, float eta, int pdgid, bool isdata)
{
	assert(_analysis==Analyze_1l2tau);

	float eff = 0.;

	if (abs(pdgid)==11) {  // Electron
		if (fabs(eta)<1.48) {
			if (isdata)
				eff = readTGraph(g_Ele_ZMassEtaLt1p48_Data, pt);
			else
				eff = readTGraph(g_Ele_ZMassEtaLt1p48_MC, pt);
		}
		else if (fabs(eta)<2.1) {
			if (isdata)
				eff = readTGraph(g_Ele_ZMassEta1p48to2p1_Data, pt);
			else
				eff = readTGraph(g_Ele_ZMassEta1p48to2p1_MC, pt);
		}
		else {
			// shouldn't be here
			assert(0);
			//if (isdata)
			//	eff = readTGraph(g_Ele_ZMassEtaGt2p1_Data, pt);
			//else
			//	eff = readTGraph(g_Ele_ZMassEtaGt2p1_MC, pt);
		}
	}
	else if (abs(pdgid)==13) {  // Muon
		if (fabs(eta)<0.9) {
			if (isdata)
				eff = readTGraph(g_Mu_ZMassEtaLt0p9_Data, pt);
			else
				eff = readTGraph(g_Mu_ZMassEtaLt0p9_MC, pt);
		}
		else if (fabs(eta)<1.2) {
			if (isdata)
				eff = readTGraph(g_Mu_ZMassEta0p9to1p2_Data, pt);
			else
				eff = readTGraph(g_Mu_ZMassEta0p9to1p2_MC, pt);
		}
		else if (fabs(eta)<2.1) {
			if (isdata)
				eff = readTGraph(g_Mu_ZMassEta1p2to2p1_Data, pt);
			else
				eff = readTGraph(g_Mu_ZMassEta1p2to2p1_MC, pt);
		}
		else // shouldn't be here
			assert(0);
	}
	else {
		std::cout << "Get_trig_eff_singleLep: Oops..." << std::endl;
		assert(0);
	}

	return eff;
}

//#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
float SFHelper::Get_LeptonIDSF_weight(const std::vector<miniLepton>& leptons)
{
	size_t nleps = 0;
	if (_analysis==Analyze_1l2tau)
		nleps = 1;
	else if (_analysis==Analyze_2lss1tau or _analysis==Analyze_2lss or
			 _analysis==Analyze_2l2tau)
		nleps = 2;
	else if (_analysis==Analyze_3l1tau or _analysis==Analyze_3l)
		nleps = 3;

	assert(leptons.size() >= nleps);

	float lepSF = 1.;
	for (size_t ilep = 0; ilep < nleps; ilep++) {
		lepSF *= Get_LeptonIDSF(leptons[ilep]);
	}

	return lepSF;
}

float SFHelper::Get_LeptonIDSF(const miniLepton& lepton)
{
	assert(not _isdata);
	
	assert(lepton.passLooseSel());

	float sf = 1.;
	sf *= Get_LeptonSF_loose(lepton);
	// SF fakeable to loose is currently assumed to be 1.
	if (lepton.passTightSel())
		sf *= Get_LeptonSF_tight_vs_loose(lepton);

	return sf;
}
//#endif
float SFHelper::Get_LeptonIDSF(float lepPt, float lepEta, bool isEle, bool isMu,
							   bool isTight)
{
	assert(not _isdata);
	
	float sf = 1.;
	sf *= Get_LeptonSF_loose(lepPt, lepEta, isEle, isMu);
	// SF fakeable to loose is currently assumed to be 1.
	if (isTight)
		sf *= Get_LeptonSF_tight_vs_loose(lepPt, lepEta, isEle, isMu);

	return sf;
}
//#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
float SFHelper::Get_LeptonSF_loose(const miniLepton& lepton)
{
	return Get_LeptonSF_loose(lepton.pt(), lepton.eta(),
							  abs(lepton.pdgId())==11,
							  abs(lepton.pdgId())==13);
}
//#endif
float SFHelper::Get_LeptonSF_loose(float lepPt, float lepEta,
								   bool isEle, bool isMu)
{
	float sf =1.;
	if (isMu) {
		sf *= evalTGraph(g_recoToLoose_leptonSF_mu1, lepEta);
		
		if (lepPt > 30)
			sf *= read2DHist(h_recoToLoose_leptonSF_mu2, lepPt, fabs(lepEta));
		else
			sf *= read2DHist(h_recoToLoose_leptonSF_mu3, lepPt, fabs(lepEta));
		
		sf *= read2DHist(h_recoToLoose_leptonSF_mu4, lepPt, fabs(lepEta));
	}
	else if (isEle) {
		if (lepPt > 20) {
			sf *= read2DHist(h_recoToLoose_leptonSF_gsf, lepEta, lepPt);
		}
		else {
			sf *= read2DHist(h_recoToLoose_leptonSF_gsf_low, lepEta, lepPt);
		}
		////////////////////////
		// UPDATE NEEDED for electron
		////////////////////////
	}
	
	return sf;
	
}

//#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
float SFHelper::Get_LeptonSF_tight_vs_loose(const miniLepton& lepton)
{
	assert(lepton.passTightSel());

	return Get_LeptonSF_tight_vs_loose(lepton.pt(), lepton.eta(),
									   abs(lepton.pdgId())==11,
									   abs(lepton.pdgId())==13);
}
//#endif
float SFHelper::Get_LeptonSF_tight_vs_loose(float lepPt, float lepEta,
											bool isEle, bool isMu)
{
	float sf = 1.;
	
	if (_analysis == Analyze_2lss1tau or _analysis == Analyze_2lss) {
		if (isMu) {
			sf = read2DHist(h_looseToTight_leptonSF_mu_2lss,
							lepPt, std::abs(lepEta));
		}
		else if (isEle) {
			sf = read2DHist(h_looseToTight_leptonSF_el_2lss,
							lepPt, std::abs(lepEta));
		}
	}
	else if (_analysis == Analyze_3l1tau or _analysis == Analyze_3l or
			 _analysis == Analyze_2l2tau or _analysis == Analyze_1l2tau) {
		if (isMu) {
			sf = read2DHist(h_looseToTight_leptonSF_mu_3l,
							lepPt, std::abs(lepEta));
		}
		else if (isEle) {
			sf = read2DHist(h_looseToTight_leptonSF_el_3l,
							lepPt, std::abs(lepEta));
		}
	}

	return sf;
}

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
float SFHelper::Get_EvtCSVWeight(const std::vector<pat::Jet>& jets,
								 const std::string& sys)
{
	std::vector<miniJet> minijets;
	
	for (const auto & j : jets) {
		TLorentzVector jp4;
		jp4.SetPtEtaPhiE(j.pt(), j.eta(), j.phi(), j.energy());
		//double csv = j.bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll");
		double csv = j.bDiscriminator("pfDeepCSVJetTags:probb")
			+ j.bDiscriminator("pfDeepCSVJetTags:probb");
		
		miniJet mj(jp4, csv, j.hadronFlavour());
		minijets.push_back(mj);
	}

	return Get_EvtCSVWeight(minijets, sys);
}
#endif

float SFHelper::Get_EvtCSVWeight(const std::vector<miniJet>& jets,
								 const std::string& sys)
{
	assert(not _isdata);

	double weight_evt = 1.;

	for (const auto & j : jets) {
		double w = Get_JetCSVWeight(j, sys);
		if (w!=0) weight_evt *= w;
	}

	return weight_evt;
}

float SFHelper::Get_JetCSVWeight(const miniJet& jet, std::string sys/*pass by copy*/)
{
	double pt = jet.pt();
	double eta = jet.eta();
	//double csv = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	//double csv = jet.bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll");	
	//int flavor = jet.hadronFlavour();

	double csv = jet.csv();
	int flavor = jet.flavor();
	
	BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;	
	if ( std::abs(flavor) == 5 )
		jf = BTagEntry::FLAV_B;
	else if ( std::abs(flavor) == 4 )
		jf = BTagEntry::FLAV_C;

	double weight_jet = 1.;

	if (sys == "JESUp" and jf != BTagEntry::FLAV_C)
		sys = "up_jes";
	else if (sys == "JESDown" and jf != BTagEntry::FLAV_C)
		sys = "down_jes";
	else if (sys == "LFUp" and jf == BTagEntry::FLAV_B)
		sys = "up_lf";
	else if (sys == "LFDown" and jf == BTagEntry::FLAV_B)
		sys = "down_lf";
	else if (sys == "HFStats1Up" and jf == BTagEntry::FLAV_B)
		sys = "up_hfstats1";
	else if (sys == "HFStats1Down" and jf == BTagEntry::FLAV_B)
		sys = "down_hfstats1";
	else if (sys == "HFStats2Up" and jf == BTagEntry::FLAV_B)
		sys = "up_hfstats2";
	else if (sys == "HFStats2Down" and jf == BTagEntry::FLAV_B)
		sys = "down_hfstats2";
	else if (sys == "HFUp" and jf == BTagEntry::FLAV_UDSG)
		sys = "up_hf";
	else if (sys == "HFDown" and jf == BTagEntry::FLAV_UDSG)
		sys = "down_hf";
	else if (sys == "LFStats1Up" and jf == BTagEntry::FLAV_UDSG)
		sys = "up_lfstats1";
	else if (sys == "LFStats1Down" and jf == BTagEntry::FLAV_UDSG)
		sys = "down_lfstats1";
	else if (sys == "LFStats2Up" and jf == BTagEntry::FLAV_UDSG)
		sys = "up_lfstats2";
	else if (sys == "LFStats2Down" and jf == BTagEntry::FLAV_UDSG)
		sys = "down_lfstats2";
	else if (sys == "cErr1Up" and jf == BTagEntry::FLAV_C)
		sys = "up_cferr1";
	else if (sys == "cErr1Down" and jf == BTagEntry::FLAV_C)
		sys = "down_cferr1";
	else if (sys == "cErr2Up" and jf == BTagEntry::FLAV_C)
		sys = "up_cferr2";
	else if (sys == "cErr2Down" and jf == BTagEntry::FLAV_C)
		sys = "down_cferr2";
	else
		sys = "central";

	weight_jet = BTagCaliReader->eval_auto_bounds(sys, jf, std::abs(eta), pt, csv);
	
	if ((sys == "central" or sys=="up_jes" or sys=="down_jes") and weight_jet==0.) {
		std::cerr << "SFHelper::Get_JetCSVWeight" << std::endl;
		std::cerr << "sys jetflavor |eta| pt csv : " << sys << " " << jf << " "
				  << std::abs(eta) << " " << pt << " " << csv << std::endl;
	}
	//assert(weight_jet > 0.);

	return weight_jet;
}

float SFHelper::Get_PUWeight(int nPU)
{
	/*
	assert(not _isdata);
	
	int xbin = h_puweight->FindBin(nPU);
	return h_puweight->GetBinContent(xbin);
	*/
	return LumiWeights_.weight(nPU);
}

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
float SFHelper::Get_TauIDSF(const pat::Tau& tau, bool isGenMatched, TString syst)
{
	return Get_TauIDSF(tau.pt(), tau.eta(), isGenMatched, syst);
}
#endif

float SFHelper::Get_TauIDSF(const miniTau& tau, TString syst)
{
	return Get_TauIDSF(tau.pt(), tau.eta(), tau.isGenMatched(), syst);
}

float SFHelper::Get_TauIDSF(float tauPt, float tauEta, bool isGenMatched, TString syst)
{
	assert(not _isdata);

	// tau ID efficiency data/MC scale factor
	// https://indico.cern.ch/event/719250/contributions/2971854/attachments/1635435/2609013/tauid_recommendations2017.pdf
	float tauEff_sf = 0.89; // dR03mvaLoose or dR03mvaMedium

	// tau ID fake rate data/MC scale factor
	float tauFR_sf = fabs(tauEta) < 1.497 ?
	    readTF(f_fakerate_tau_mvaM_etaL_ratio, tauPt) :
		readTF(f_fakerate_tau_mvaM_etaH_ratio, tauPt);

	if (syst=="FRjt_normUp") {
		tauFR_sf = fabs(tauEta) < 1.497 ?
			readTF(f_fakerate_tau_mvaM_etaL_ratio_normUp, tauPt) :
			readTF(f_fakerate_tau_mvaM_etaH_ratio_normUp, tauPt);
	}
	else if (syst=="FRjt_normDown") {
		tauFR_sf = fabs(tauEta) < 1.497 ?
			readTF(f_fakerate_tau_mvaM_etaL_ratio_normDown, tauPt) :
			readTF(f_fakerate_tau_mvaM_etaH_ratio_normDown, tauPt);
	}
	else if (syst=="FRjt_shapeUp") {
		tauFR_sf = fabs(tauEta) < 1.497 ?
			readTF(f_fakerate_tau_mvaM_etaL_ratio_shapeUp, tauPt) :
			readTF(f_fakerate_tau_mvaM_etaH_ratio_shapeUp, tauPt);
	}
	else if (syst=="FRjt_shapeDown") {
		tauFR_sf = fabs(tauEta) < 1.497 ?
			readTF(f_fakerate_tau_mvaM_etaL_ratio_shapeDown, tauPt) :
			readTF(f_fakerate_tau_mvaM_etaH_ratio_shapeDown, tauPt);
	}

	return (isGenMatched ? tauEff_sf : tauFR_sf);
}

float SFHelper::Get_TauIDSF_weight(const std::vector<miniTau>& taus, TString syst)
{
	size_t ntaus = 1;
	if (_analysis==Analyze_1l2tau or _analysis==Analyze_2l2tau)
		ntaus = 2;
	else if (_analysis==Analyze_2lss or _analysis==Analyze_3l)
		ntaus = 0;
	// TODO: ttW,ttZ control regions ?

	assert(taus.size()>=ntaus);

	float tauSF = 1.;

	for (size_t itau = 0; itau < ntaus; itau++) {
		tauSF *= Get_TauIDSF(taus[itau],syst);
	}

	return tauSF;
}

/*float SFHelper::Get_MCWeight()
{
	return 0.;
	}*/

float SFHelper::Get_FakeRate_lep(const miniLepton& lepton, TString syst)
{
	bool isEle = abs(lepton.pdgId())==11;
	bool isMu = abs(lepton.pdgId())==13;
	
	return Get_FakeRate_lep(lepton.conept(), lepton.eta(), isEle, isMu, syst);
}

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
float SFHelper::Get_FakeRate_tau(const pat::Tau& tau, TString syst)
{
	return Get_FakeRate_tau(tau.pt(),tau.eta(), syst);
}

float SFHelper::Get_FR_weight(const std::vector<miniLepton>& leps,
							  const std::vector<pat::Tau>& taus,
							  TString syst)
{
	std::vector<miniTau> mtaus;
	bool addDaughters = false;
	for (const auto & tau : taus) {
		miniTau mt(tau, addDaughters);
		mtaus.push_back(mt);
	}

	return Get_FR_weight(leps, mtaus, syst);
}

float SFHelper::Get_ChargeFlipWeight(const std::vector<miniLepton>& leps,
									 const std::vector<pat::Tau>& taus)
{
	assert(leps.size() >= 2);
	assert(leps[0].passFakeableSel() and leps[1].passFakeableSel());
	assert(taus.size() >= 1);

	return Get_ChargeFlipWeight(leps, taus[0].charge());
}
#endif

/*
float SFHelper::Get_FR_weight(const miniLepton& lep1, const miniLepton& lep2,
							  TString syst)
{
	bool lep1IsEle = abs(lep1.pdgId()) == 11;
	bool lep1IsMu = abs(lep1.pdgId()) == 13;
	bool lep2IsEle = abs(lep2.pdgId()) == 11;
	bool lep2IsMu = abs(lep2.pdgId()) == 13;

	return Get_FR_weight(
						 lep1.conept(),lep1.eta(),lep1IsEle,lep1IsMu,
						 lep1.passTightSel(),
						 lep2.conept(),lep2.eta(),lep2IsEle,lep2IsMu,
						 lep2.passTightSel(),
						 syst
						 );
}
*/

float SFHelper::Get_FakeRate_lep(float lepConePt, float lepEta,
							 bool isEle, bool isMuon, TString syst)
{
	float fakerate = 0;

	if (syst=="FRe_normUp" and isEle) {
		fakerate = read2DHist(h_fakerate_el_normUp, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRe_normDown" and isEle) {
		fakerate = read2DHist(h_fakerate_el_normDown, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRe_ptUp" and isEle) {
		fakerate = read2DHist(h_fakerate_el_ptUp, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRe_ptDown" and isEle) {
		fakerate = read2DHist(h_fakerate_el_ptDown, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRe_beUp" and isEle) {
		fakerate = read2DHist(h_fakerate_el_beUp, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRe_beDown" and isEle) {
		fakerate = read2DHist(h_fakerate_el_beDown, lepConePt, std::abs(lepEta));
	}
	//else if (syst=="FRe_bUp" and isEle) {
	//	fakerate = read2DHist(h_fakerate_el_bUp, lepConePt, std::abs(lepEta));
	//}
	//else if (syst=="FRe_bDown" and isEle) {
	//	fakerate = read2DHist(h_fakerate_el_bDown, lepConePt, std::abs(lepEta));
	//}
	//else if (syst=="FRe_ecUp" and isEle) {
	//	fakerate = read2DHist(h_fakerate_el_ecUp, lepConePt, std::abs(lepEta));
	//}
	//else if (syst=="FRe_ecDown" and isEle) {
	//	fakerate = read2DHist(h_fakerate_el_ecDown, lepConePt, std::abs(lepEta));
	//}
	else if (syst=="FRm_normUp" and isMuon) {
		fakerate = read2DHist(h_fakerate_mu_normUp, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRm_normDown" and isMuon) {
		fakerate = read2DHist(h_fakerate_mu_normDown, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRm_ptUp" and isMuon) {
		fakerate = read2DHist(h_fakerate_mu_ptUp, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRm_ptDown" and isMuon) {
		fakerate = read2DHist(h_fakerate_mu_ptDown, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRm_beUp" and isMuon) {
		fakerate = read2DHist(h_fakerate_mu_beUp, lepConePt, std::abs(lepEta));
	}
	else if (syst=="FRm_beDown" and isMuon) {
		fakerate = read2DHist(h_fakerate_mu_beDown, lepConePt, std::abs(lepEta));
	}
	//else if (syst=="FRm_bUp" and isMuon) {
	//	fakerate = read2DHist(h_fakerate_mu_bUp, lepConePt, std::abs(lepEta));
	//}
	//else if (syst=="FRm_bDown" and isMuon) {
	//	fakerate = read2DHist(h_fakerate_mu_bDown, lepConePt, std::abs(lepEta));
	//}
	//else if (syst=="FRm_ecUp" and isMuon) {
	//	fakerate = read2DHist(h_fakerate_mu_ecUp, lepConePt, std::abs(lepEta));
	//}
	//else if (syst=="FRm_ecDown" and isMuon) {
	//	fakerate = read2DHist(h_fakerate_mu_ecDown, lepConePt, std::abs(lepEta));
	//}
	else if ((syst=="FR_TT" or syst=="FR_el_TT_mu_QCD") and isEle) {
		fakerate = read2DHist(h_fakerate_el_TT, lepConePt, std::abs(lepEta));
	}
	else if ((syst=="FR_QCD" or syst=="FR_el_QCD_mu_TT") and isEle) {
		fakerate = read2DHist(h_fakerate_el_QCD, lepConePt, std::abs(lepEta));
	}
	else if ((syst=="FR_TT" or syst=="FR_el_QCD_mu_TT") and isMuon) {
		fakerate = read2DHist(h_fakerate_mu_TT, lepConePt, std::abs(lepEta));
	}
	else if ((syst=="FR_QCD" or syst=="FR_el_TT_mu_QCD") and isMuon) {
		fakerate = read2DHist(h_fakerate_mu_QCD, lepConePt, std::abs(lepEta));
	}
	//
	else if (isEle)
		fakerate = read2DHist(h_fakerate_el, lepConePt, std::abs(lepEta));
	else if (isMuon)
		fakerate = read2DHist(h_fakerate_mu, lepConePt, std::abs(lepEta));
	
	if (lepConePt < 10.) return 0.;
	
	return fakerate;
}

float SFHelper::Get_FakeRate_tau(float tauPt, float tauEta, TString syst)
{
	float fr_mc = 0;
	float ratio = 0;

	if (std::abs(tauEta) < 1.479) {
		fr_mc = readTGraph(g_fakerate_tau_mvaM_etaL_mc, tauPt);
	}
	else {
		fr_mc = readTGraph(g_fakerate_tau_mvaM_etaH_mc, tauPt);
	}
	
	if (syst=="FRjt_normUp") {
		if (std::abs(tauEta) < 1.479) {
			ratio = readTF(f_fakerate_tau_mvaM_etaL_ratio_normUp, tauPt);
		}
		else {
			ratio = readTF(f_fakerate_tau_mvaM_etaH_ratio_normUp, tauPt);
		}
	}
	else if (syst=="FRjt_normDown") {
		if (std::abs(tauEta) < 1.479) {
			ratio = readTF(f_fakerate_tau_mvaM_etaL_ratio_normDown, tauPt);
		}
		else {
			ratio = readTF(f_fakerate_tau_mvaM_etaH_ratio_normDown, tauPt);
		}
	}
	else if (syst=="FRjt_shapeUp") {
		if (std::abs(tauEta) < 1.479) {
			ratio = readTF(f_fakerate_tau_mvaM_etaL_ratio_shapeUp, tauPt);
		}
		else {
			ratio = readTF(f_fakerate_tau_mvaM_etaH_ratio_shapeUp, tauPt);
		}
	}
	else if (syst=="FRjt_shapeDown") {
		if (std::abs(tauEta) < 1.479) {
			ratio = readTF(f_fakerate_tau_mvaM_etaL_ratio_shapeDown, tauPt);
		}
		else {
			ratio = readTF(f_fakerate_tau_mvaM_etaH_ratio_shapeDown, tauPt);
		}
	}
	else {		
		if (std::abs(tauEta) < 1.479) {
			ratio = readTF(f_fakerate_tau_mvaM_etaL_ratio, tauPt);
		}
		else {
			ratio = readTF(f_fakerate_tau_mvaM_etaH_ratio, tauPt);
		}
	}
	
	return fr_mc * ratio;
}

float SFHelper::Get_FR_weight(const std::vector<miniLepton>& leps,
							  const std::vector<miniTau>& taus,
							  TString syst)
{
	float F1=-1., F2=-1., F3=-1.,F4=-1.;
	float f1=0., f2=0., f3=0., f4=0.;
	float FR_weight = 0.;

	if (_analysis==Analyze_1l2tau) {
		assert(leps.size() >= 1);
		assert(taus.size() >= 2);

		f1 = Get_FakeRate_lep(leps[0], syst);
		f2 = Get_FakeRate_tau(taus[0].pt(), taus[0].eta(), syst);
		f3 = Get_FakeRate_tau(taus[1].pt(), taus[1].eta(), syst);

		F1 = leps[0].passTightSel() ? -1. : f1/(1.-f1);
		F2 = taus[0].passTightSel() ? -1. : f2/(1.-f2);
		F3 = taus[1].passTightSel() ? -1. : f3/(1.-f3);

		if (_debug) {
			std::cout << "lep pdgid passTight? : " << leps[0].pdgId() << " "
					  << leps[0].passTightSel() << std::endl;
			std::cout << "f1 F1 : " << f1 << " " << F1 << std::endl;
			std::cout << "tau0 passTight? : " << taus[0].passTightSel() << std::endl;
			std::cout << "f2 F2 : " << f2 << " " << F2 << std::endl;
			std::cout << "tau1 passTight? : " << taus[1].passTightSel() << std::endl;
			std::cout << "f3 F3 : " << f3 << " " << F3 << std::endl;
		}

		FR_weight = F1 * F2 * F3;
		if (_debug) std::cout << "FR_weight : " << F1 * F2 * F3 << std::endl;
	}
	else if (_analysis==Analyze_2lss1tau or _analysis==Analyze_2lss) {
		assert(leps.size() >= 2);
		assert(leps[0].passFakeableSel() and leps[1].passFakeableSel());

		f1 = Get_FakeRate_lep(leps[0], syst);
		f2 = Get_FakeRate_lep(leps[1], syst);

		F1 = leps[0].passTightSel() ? -1. : f1/(1.-f1);
		F2 = leps[1].passTightSel() ? -1. : f2/(1.-f2);

		if (_debug) {
			std::cout << "lep0 pdgid passTight? : " << leps[0].pdgId() << " "
					  << leps[0].passTightSel() << std::endl;
			std::cout << "f1 F1 : " << f1 << " " << F1 << std::endl;
			std::cout << "lep1 pdgid passTight? : " << leps[1].pdgId() << " "
					  << leps[1].passTightSel() << std::endl;
			std::cout << "f2 F2 : " << f2 << " " << F2 << std::endl;
		}

		FR_weight = -1. * F1 * F2;
		if (_debug) std::cout << "FR_weight : " << -1 * F1 * F2 << std::endl;
	}
	else if (_analysis==Analyze_3l1tau or _analysis==Analyze_3l) {
		assert(leps.size() >= 3);
		assert(leps[0].passFakeableSel() and leps[1].passFakeableSel() and
			   leps[2].passFakeableSel());

		f1 = Get_FakeRate_lep(leps[0], syst);
		f2 = Get_FakeRate_lep(leps[1], syst);
		f3 = Get_FakeRate_lep(leps[2], syst);

		F1 = leps[0].passTightSel() ? -1. : f1/(1.-f1);
		F2 = leps[1].passTightSel() ? -1. : f2/(1.-f2);
		F3 = leps[2].passTightSel() ? -1. : f2/(1.-f3);

		if (_debug) {
			std::cout << "lep0 pdgid passTight? : " << leps[0].pdgId() << " "
					  << leps[0].passTightSel() << std::endl;
			std::cout << "f1 F1 : " << f1 << " " << F1 << std::endl;
			std::cout << "lep1 pdgid passTight? : " << leps[1].pdgId() << " "
					  << leps[1].passTightSel() << std::endl;
			std::cout << "f2 F2 : " << f2 << " " << F2 << std::endl;
			std::cout << "lep2 pdgid passTight? : " << leps[2].pdgId() << " "
					  << leps[2].passTightSel() << std::endl;
			std::cout << "f3 F3 : " << f3 << " " << F3 << std::endl;
		}

		FR_weight = F1 * F2 * F3;
		if (_debug) std::cout << "FR_weight : " << F1 * F2 * F3 << std::endl;
	}
	else if (_analysis==Analyze_2l2tau) {
		assert(leps.size() >= 2);
		assert(taus.size() >= 2);
		
		f1 = Get_FakeRate_lep(leps[0], syst);
		f2 = Get_FakeRate_lep(leps[1], syst);
		f3 = Get_FakeRate_tau(taus[0].pt(), taus[0].eta(), syst);
		f4 = Get_FakeRate_tau(taus[1].pt(), taus[1].eta(), syst);
		
		F1 = leps[0].passTightSel() ? -1. : f1/(1.-f1);
		F2 = leps[1].passTightSel() ? -1. : f2/(1.-f2);
		F3 = taus[0].passTightSel() ? -1. : f3/(1.-f3);
		F4 = taus[1].passTightSel() ? -1. : f4/(1.-f4);
		
		if (_debug) {
			std::cout << "lep0 pdgid passTight? : " << leps[0].pdgId() << " "
					  << leps[0].passTightSel() << std::endl;
			std::cout << "f1 F1 : " << f1 << " " << F1 << std::endl;
			std::cout << "lep1 pdgid passTight? : " << leps[1].pdgId() << " "
					  << leps[1].passTightSel() << std::endl;
			std::cout << "f2 F2 : " << f2 << " " << F2 << std::endl;
			std::cout << "tau0 passTight? : " << taus[0].passTightSel()
					  << std::endl;
			std::cout << "f3 F3 : " << f3 << " " << F3 << std::endl;
			std::cout << "tau1 passTight? : " << taus[1].passTightSel()
					  << std::endl;
			std::cout << "f4 F4 : " << f4 << " " << F4 << std::endl;
		}

		FR_weight = -1 * F1 * F2 * F3 * F4;
		if (_debug) std::cout << "FR_weight : " << (-1*F1*F2*F3*F4) << std::endl;
	}	

	return FR_weight;
}

/*
float SFHelper::Get_FR_weight(float lep1ConePt, float lep1Eta, bool lep1IsEle,
							  bool lep1IsMu, bool lep1IsTight,
							  float lep2ConePt, float lep2Eta, bool lep2IsEle,
							  bool lep2IsMu, bool lep2IsTight,
							  TString syst)
{
	float f1 = Get_FakeRate(lep1ConePt, lep1Eta, lep1IsEle, lep1IsMu, syst);
	float f2 = Get_FakeRate(lep2ConePt, lep2Eta, lep2IsEle, lep2IsMu, syst);
	float F1 = lep1IsTight ? -1. : f1/(1.-f1);
	float F2 = lep2IsTight ? -1. : f2/(1.-f2);

	float weight = -1. * F1 * F2;

	if (lep1IsTight and lep2IsTight) {
		assert((lep1IsEle and lep2IsEle) or (lep1IsMu and lep2IsMu));
		weight = 0;
	}
	
	return weight;
}*/

float SFHelper::Get_ChargeFlipWeight(const std::vector<miniLepton>& leps,
									 int tauCharge)
{
	assert(leps.size() >= 2);
	assert(leps[0].passFakeableSel() and leps[1].passFakeableSel());

	float P1_misCharge =
		Get_EleChargeMisIDProb(leps[0], tauCharge);
	float P2_misCharge =
	    Get_EleChargeMisIDProb(leps[1], tauCharge);

	if (_debug) {
		std::cout << "tau charge : " << tauCharge << std::endl;
		std::cout << "lep0 pt conept eta charge : " << leps[0].pt() << " "
				  << leps[0].conept() << " " << leps[0].eta() << " "
				  << leps[0].charge() << std::endl;
		std::cout << "p1_mischarge : " << P1_misCharge << std::endl;
		std::cout << "lep1 pt conept eta charge : " << leps[1].pt() << " "
				  << leps[1].conept() << " " << leps[1].eta() << " "
				  << leps[1].charge() << std::endl;
		std::cout << "p2_mischarge : " << P2_misCharge << std::endl;
	}

	return P1_misCharge + P2_misCharge;
}

float SFHelper::Get_EleChargeMisIDProb(const miniLepton& lepton, int tauCharge)
{
	// muon
	if (abs(lepton.pdgId())==13) return 0.;
	
	// electron
	assert(abs(lepton.pdgId())==11);

	if (_selection==Application_Flip_2lss1tau or _selection==Control_FlipAR_ttW) {
		// only apply the charge flip rate to the electron that is same sign as tau
		// due to the tau charge requirement in signal region
		if (lepton.charge() * tauCharge < 0) return 0;
	}
	else {
		assert(_selection==Control_FlipAR_2lss1tau);
		// only apply to the electron that is opposite sign as tau
		if (lepton.charge() * tauCharge > 0) return 0;
	}
	
	return read2DHist(h_chargeMisId, lepton.pt(), std::abs(lepton.eta()));
}

float SFHelper::read2DHist(TH2* h2d, float x, float y)
{
	TAxis* xaxis = h2d->GetXaxis();
	int nbinx = xaxis->GetNbins();
	int xbin = xaxis->FindBin(x);
    if (xbin < 1) xbin = 1;
	if (xbin > nbinx) xbin = nbinx;
	
	TAxis* yaxis = h2d->GetYaxis();
	int nbiny = yaxis->GetNbins();
	int ybin = yaxis->FindBin(y);
    if (ybin < 1) ybin = 1;
	if (ybin > nbiny) ybin = nbiny;
	
	float result = h2d->GetBinContent(xbin, ybin);

	//result += errVar * h2d->GetBinError(xbin, ybin);

	return result;
}

float SFHelper::evalTGraph(TGraphAsymmErrors* graph, float x)
{
	float x1 = std::max(float(graph->GetXaxis()->GetXmin()+1e-5),
						std::min(float(graph->GetXaxis()->GetXmax()-1e-5), x)
						);
	return graph->Eval(x1);
}

float SFHelper::readTGraph(TGraphAsymmErrors* graph, float x)
{
	int nPoints = graph -> GetN();
	double xp, yp = -1.;
	for (int i=0; i<nPoints; ++i) {
		int ip = graph -> GetPoint(i,xp,yp);
		float xerr_h = graph->GetErrorXhigh(i);
		float xerr_l = graph->GetErrorXlow(i);

		assert(xerr_h-xerr_l>=0);
		
		if (ip != -1) {
			if (x>xp-xerr_l and x<xp+xerr_h) return yp;
		}
	}
	return yp;
}

float SFHelper::readTF(TF1* f, float x)
{
	//float x1 = std::max(float(f->GetXaxis()->GetXmin()),
	//					std::min(float(f->GetXaxis()->GetXmax()), x)
	//					);
	return f->Eval(x);
}
