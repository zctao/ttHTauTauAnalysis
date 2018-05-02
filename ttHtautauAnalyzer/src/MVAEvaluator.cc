#include "../interface/MVAEvaluator.h"

using namespace std;

const TString MVAEvaluator::data_directory_ = (string(getenv("CMSSW_BASE"))+"/src/ttHTauTauAnalysis/ttHtautauAnalyzer/dataFiles/mva/").c_str();

MVAEvaluator::MVAEvaluator()
{
	setup_tmva_reader_HTT();
	setup_tmva_reader_1l2tau_BDT1();
	setup_tmva_reader_1l2tau_BDT2();
	setup_tmva_reader_2lss1tau_BDT1();
	setup_tmva_reader_2lss1tau_BDT2();
	setup_tmva_reader_2lss1tau_BDT3();
	setup_tmva_reader_2lss1tau_BDT4();
	setup_tmva_reader_2lss1tau_BDT5();
	setup_tmva_reader_3l1tau_BDT1();
	setup_tmva_reader_3l1tau_BDT2();
	setup_tmva_reader_3l1tau_BDT3();
	setup_tmva_reader_2l2tau_BDT1();
	setup_tmva_reader_2l2tau_BDT2();
	setup_tmva_reader_2l2tau_BDT3();
}

void MVAEvaluator::setup_tmva_reader_HTT()
{
	const TString weights =
		data_directory_ + "HadTopTagger_resolved_XGB_CSV_sort_withKinFit.xml";
	reader_HTT_ = new TMVA::Reader("!Color:!Silent");

	reader_HTT_ -> AddVariable("CSV_b", &(inputVars_HTT_[0]));
	reader_HTT_ -> AddVariable("qg_Wj2", &(inputVars_HTT_[1]));
	reader_HTT_ -> AddVariable("pT_bWj1Wj2", &(inputVars_HTT_[2]));
	reader_HTT_ -> AddVariable("m_Wj1Wj2", &(inputVars_HTT_[3]));
	reader_HTT_ -> AddVariable("nllKinFit", &(inputVars_HTT_[4]));
	reader_HTT_ -> AddVariable("pT_b_o_kinFit_pT_b", &(inputVars_HTT_[5]));
	reader_HTT_ -> AddVariable("pT_Wj2", &(inputVars_HTT_[6]));

	reader_HTT_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_HTT(float Vars[7])
{
	for (size_t i=0; i<7; ++i)
		inputVars_HTT_[i] = Vars[i];

	return reader_HTT_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_1l2tau_BDT1()
{
	const TString weights =
		data_directory_ + "1l2tau/1l_2tau_XGB_plainKin_evtLevelTT_TTH_13Var.xml";
	reader_1l2tau_BDT1_ = new TMVA::Reader("!Color:!Silent");

	reader_1l2tau_BDT1_ ->AddVariable("avg_dr_jet", &(inputVars_1l2tau_BDT1_[0]));
	reader_1l2tau_BDT1_ ->AddVariable("dr_taus", &(inputVars_1l2tau_BDT1_[1]));
	reader_1l2tau_BDT1_ ->AddVariable("ptmiss", &(inputVars_1l2tau_BDT1_[2]));
	reader_1l2tau_BDT1_ ->AddVariable("mT_lep", &(inputVars_1l2tau_BDT1_[3]));
	reader_1l2tau_BDT1_ ->AddVariable("nJet", &(inputVars_1l2tau_BDT1_[4]));
	reader_1l2tau_BDT1_ ->AddVariable("mTauTauVis", &(inputVars_1l2tau_BDT1_[5]));
	reader_1l2tau_BDT1_ ->AddVariable("mindr_lep_jet", &(inputVars_1l2tau_BDT1_[6]));
	reader_1l2tau_BDT1_ ->AddVariable("mindr_tau1_jet", &(inputVars_1l2tau_BDT1_[7]));
	reader_1l2tau_BDT1_ ->AddVariable("mindr_tau2_jet", &(inputVars_1l2tau_BDT1_[8]));
	reader_1l2tau_BDT1_ ->AddVariable("dr_lep_tau_lead", &(inputVars_1l2tau_BDT1_[9]));
	reader_1l2tau_BDT1_ ->AddVariable("nBJetLoose", &(inputVars_1l2tau_BDT1_[10]));
	reader_1l2tau_BDT1_ ->AddVariable("tau1_pt", &(inputVars_1l2tau_BDT1_[11]));
	reader_1l2tau_BDT1_ ->AddVariable("tau2_pt", &(inputVars_1l2tau_BDT1_[12]));

	reader_1l2tau_BDT1_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_1l2tau_BDT1(float Vars[13])
{
	for (size_t i=0; i<13; ++i)
		inputVars_1l2tau_BDT1_[i] = Vars[i];

	return reader_1l2tau_BDT1_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_1l2tau_BDT2()
{
	const TString weights =
		data_directory_ + "1l2tau/1l_2tau_XGB_HTT_evtLevelSUM_TTH_VT_17Var.xml";
	reader_1l2tau_BDT2_ = new TMVA::Reader("!Color:!Silent");

	reader_1l2tau_BDT2_ -> AddVariable("avg_dr_jet", &(inputVars_1l2tau_BDT2_[0]));
	reader_1l2tau_BDT2_ -> AddVariable("dr_taus", &(inputVars_1l2tau_BDT2_[1]));
	reader_1l2tau_BDT2_ -> AddVariable("ptmiss", &(inputVars_1l2tau_BDT2_[2]));
	reader_1l2tau_BDT2_ -> AddVariable("lep_conePt", &(inputVars_1l2tau_BDT2_[3]));
	reader_1l2tau_BDT2_ -> AddVariable("mT_lep", &(inputVars_1l2tau_BDT2_[4]));
	reader_1l2tau_BDT2_ -> AddVariable("mTauTauVis", &(inputVars_1l2tau_BDT2_[5]));
	reader_1l2tau_BDT2_ -> AddVariable("mindr_lep_jet", &(inputVars_1l2tau_BDT2_[6]));
	reader_1l2tau_BDT2_ -> AddVariable("mindr_tau1_jet", &(inputVars_1l2tau_BDT2_[7]));
	reader_1l2tau_BDT2_ -> AddVariable("mindr_tau2_jet", &(inputVars_1l2tau_BDT2_[8]));
	reader_1l2tau_BDT2_ -> AddVariable("dr_lep_tau_ss", &(inputVars_1l2tau_BDT2_[9]));
	reader_1l2tau_BDT2_ -> AddVariable("dr_lep_tau_lead", &(inputVars_1l2tau_BDT2_[10]));
	reader_1l2tau_BDT2_ -> AddVariable("costS_tau", &(inputVars_1l2tau_BDT2_[11]));
	reader_1l2tau_BDT2_ -> AddVariable("nBJetLoose", &(inputVars_1l2tau_BDT2_[12]));
	reader_1l2tau_BDT2_ -> AddVariable("tau1_pt", &(inputVars_1l2tau_BDT2_[13]));
	reader_1l2tau_BDT2_ -> AddVariable("tau2_pt", &(inputVars_1l2tau_BDT2_[14]));
	reader_1l2tau_BDT2_ -> AddVariable("HTT", &(inputVars_1l2tau_BDT2_[15]));
	reader_1l2tau_BDT2_ -> AddVariable("HadTop_pt", &(inputVars_1l2tau_BDT2_[16]));
	
	reader_1l2tau_BDT2_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_1l2tau_BDT2(float Vars[17])
{
	for (size_t i=0; i<17; ++i)
		inputVars_1l2tau_BDT2_[i] = Vars[i];

	return reader_1l2tau_BDT2_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_2lss1tau_BDT1()
{
	const TString weights =
		data_directory_ + "2lss1tau/2lss_1tau_XGB_plainKin_evtLevelTTV_TTH_15Var.xml";
	reader_2lss1tau_BDT1_ = new TMVA::Reader("!Color:!Silent");

	reader_2lss1tau_BDT1_ -> AddVariable("avg_dr_jet", &(inputVars_2lss1tau_BDT1_[0]));
	reader_2lss1tau_BDT1_ -> AddVariable("dr_lep1_tau", &(inputVars_2lss1tau_BDT1_[1]));
	reader_2lss1tau_BDT1_ -> AddVariable("dr_leps", &(inputVars_2lss1tau_BDT1_[2]));
	reader_2lss1tau_BDT1_ -> AddVariable("lep1_conePt", &(inputVars_2lss1tau_BDT1_[3]));
	reader_2lss1tau_BDT1_ -> AddVariable("lep2_conePt", &(inputVars_2lss1tau_BDT1_[4]));
	reader_2lss1tau_BDT1_ -> AddVariable("mT_lep1", &(inputVars_2lss1tau_BDT1_[5]));
	reader_2lss1tau_BDT1_ -> AddVariable("mT_lep2", &(inputVars_2lss1tau_BDT1_[6]));
	reader_2lss1tau_BDT1_ -> AddVariable("mTauTauVis1", &(inputVars_2lss1tau_BDT1_[7]));
	reader_2lss1tau_BDT1_ -> AddVariable("mTauTauVis2", &(inputVars_2lss1tau_BDT1_[8]));
	reader_2lss1tau_BDT1_ -> AddVariable("mindr_lep1_jet", &(inputVars_2lss1tau_BDT1_[9]));
	reader_2lss1tau_BDT1_ -> AddVariable("mindr_lep2_jet", &(inputVars_2lss1tau_BDT1_[10]));
	reader_2lss1tau_BDT1_ -> AddVariable("mindr_tau_jet", &(inputVars_2lss1tau_BDT1_[11]));
	reader_2lss1tau_BDT1_ -> AddVariable("ptmiss", &(inputVars_2lss1tau_BDT1_[12]));
	reader_2lss1tau_BDT1_ -> AddVariable("max_lep_eta", &(inputVars_2lss1tau_BDT1_[13]));
	reader_2lss1tau_BDT1_ -> AddVariable("tau_pt", &(inputVars_2lss1tau_BDT1_[14]));
	
	reader_2lss1tau_BDT1_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_2lss1tau_BDT1(float Vars[15])
{
	for (size_t i=0; i<15; ++i)
		inputVars_2lss1tau_BDT1_[i] = Vars[i];

	return reader_2lss1tau_BDT1_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_2lss1tau_BDT2()
{
	const TString weights =
		data_directory_ + "2lss1tau/2lss_1tau_XGB_plainKin_evtLevelTT_TTH_16Var.xml";
	reader_2lss1tau_BDT2_ = new TMVA::Reader("!Color:!Silent");

	reader_2lss1tau_BDT2_ -> AddVariable("avg_dr_jet", &(inputVars_2lss1tau_BDT2_[0]));
	reader_2lss1tau_BDT2_ -> AddVariable("dr_lep1_tau", &(inputVars_2lss1tau_BDT2_[1]));
	reader_2lss1tau_BDT2_ -> AddVariable("dr_lep2_tau", &(inputVars_2lss1tau_BDT2_[2]));
	reader_2lss1tau_BDT2_ -> AddVariable("dr_leps", &(inputVars_2lss1tau_BDT2_[3]));
	reader_2lss1tau_BDT2_ -> AddVariable("lep1_conePt", &(inputVars_2lss1tau_BDT2_[4]));
	reader_2lss1tau_BDT2_ -> AddVariable("lep2_conePt", &(inputVars_2lss1tau_BDT2_[5]));
	reader_2lss1tau_BDT2_ -> AddVariable("mT_lep2", &(inputVars_2lss1tau_BDT2_[6]));
	reader_2lss1tau_BDT2_ -> AddVariable("mTauTauVis1", &(inputVars_2lss1tau_BDT2_[7]));
	reader_2lss1tau_BDT2_ -> AddVariable("mTauTauVis2", &(inputVars_2lss1tau_BDT2_[8]));
	reader_2lss1tau_BDT2_ -> AddVariable("mbb", &(inputVars_2lss1tau_BDT2_[9]));
	reader_2lss1tau_BDT2_ -> AddVariable("mindr_lep1_jet", &(inputVars_2lss1tau_BDT2_[10]));
	reader_2lss1tau_BDT2_ -> AddVariable("mindr_lep2_jet", &(inputVars_2lss1tau_BDT2_[11]));
	reader_2lss1tau_BDT2_ -> AddVariable("mindr_tau_jet", &(inputVars_2lss1tau_BDT2_[12]));
	reader_2lss1tau_BDT2_ -> AddVariable("nJet", &(inputVars_2lss1tau_BDT2_[13]));
	reader_2lss1tau_BDT2_ -> AddVariable("ptmiss", &(inputVars_2lss1tau_BDT2_[14]));
	reader_2lss1tau_BDT2_ -> AddVariable("tau_pt", &(inputVars_2lss1tau_BDT2_[15]));

	reader_2lss1tau_BDT2_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_2lss1tau_BDT2(float Vars[16])
{
	for (size_t i=0; i<16; ++i)
		inputVars_2lss1tau_BDT2_[i] = Vars[i];

	return reader_2lss1tau_BDT2_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_2lss1tau_BDT3()
{
	const TString weights =
		data_directory_ + "2lss1tau/2lss_1tau_XGB_plainKin_evtLevelSUM_TTH_M_18Var.xml";
	reader_2lss1tau_BDT3_ = new TMVA::Reader("!Color:!Silent");

	reader_2lss1tau_BDT3_ -> AddVariable("avg_dr_jet", &(inputVars_2lss1tau_BDT3_[0]));
	reader_2lss1tau_BDT3_ -> AddVariable("dr_lep1_tau", &(inputVars_2lss1tau_BDT3_[1]));
	reader_2lss1tau_BDT3_ -> AddVariable("dr_lep2_tau", &(inputVars_2lss1tau_BDT3_[2]));
	reader_2lss1tau_BDT3_ -> AddVariable("dr_leps", &(inputVars_2lss1tau_BDT3_[3]));
	reader_2lss1tau_BDT3_ -> AddVariable("lep1_conePt", &(inputVars_2lss1tau_BDT3_[4]));
	reader_2lss1tau_BDT3_ -> AddVariable("lep2_conePt", &(inputVars_2lss1tau_BDT3_[5]));
	reader_2lss1tau_BDT3_ -> AddVariable("mT_lep1", &(inputVars_2lss1tau_BDT3_[6]));
	reader_2lss1tau_BDT3_ -> AddVariable("mT_lep2", &(inputVars_2lss1tau_BDT3_[7]));
	reader_2lss1tau_BDT3_ -> AddVariable("mTauTauVis1", &(inputVars_2lss1tau_BDT3_[8]));
	reader_2lss1tau_BDT3_ -> AddVariable("mTauTauVis2", &(inputVars_2lss1tau_BDT3_[9]));
	reader_2lss1tau_BDT3_ -> AddVariable("max_lep_eta", &(inputVars_2lss1tau_BDT3_[10]));
	reader_2lss1tau_BDT3_ -> AddVariable("mbb", &(inputVars_2lss1tau_BDT3_[11]));
	reader_2lss1tau_BDT3_ -> AddVariable("mindr_lep1_jet", &(inputVars_2lss1tau_BDT3_[12]));
	reader_2lss1tau_BDT3_ -> AddVariable("mindr_lep2_jet", &(inputVars_2lss1tau_BDT3_[13]));
	reader_2lss1tau_BDT3_ -> AddVariable("mindr_tau_jet", &(inputVars_2lss1tau_BDT3_[14]));
	reader_2lss1tau_BDT3_ -> AddVariable("nJet", &(inputVars_2lss1tau_BDT3_[15]));
	reader_2lss1tau_BDT3_ -> AddVariable("ptmiss", &(inputVars_2lss1tau_BDT3_[16]));
	reader_2lss1tau_BDT3_ -> AddVariable("tau_pt", &(inputVars_2lss1tau_BDT3_[17]));

	reader_2lss1tau_BDT3_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_2lss1tau_BDT3(float Vars[18])
{
	for (size_t i=0; i<18; ++i)
		inputVars_2lss1tau_BDT3_[i] = Vars[i];

	return reader_2lss1tau_BDT3_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_2lss1tau_BDT4()
{
	const TString weights =
		data_directory_ + "2lss1tau/2lss_1tau_XGB_HTT_evtLevelSUM_TTH_M_19Var.xml";
	reader_2lss1tau_BDT4_ = new TMVA::Reader("!Color:!Silent");

	reader_2lss1tau_BDT4_ -> AddVariable("avg_dr_jet", &(inputVars_2lss1tau_BDT4_[0]));
	reader_2lss1tau_BDT4_ -> AddVariable("dr_lep1_tau", &(inputVars_2lss1tau_BDT4_[1]));
	reader_2lss1tau_BDT4_ -> AddVariable("dr_lep2_tau", &(inputVars_2lss1tau_BDT4_[2]));
	reader_2lss1tau_BDT4_ -> AddVariable("dr_leps", &(inputVars_2lss1tau_BDT4_[3]));
	reader_2lss1tau_BDT4_ -> AddVariable("lep2_conePt", &(inputVars_2lss1tau_BDT4_[4]));
	reader_2lss1tau_BDT4_ -> AddVariable("mT_lep1", &(inputVars_2lss1tau_BDT4_[5]));
	reader_2lss1tau_BDT4_ -> AddVariable("mT_lep2", &(inputVars_2lss1tau_BDT4_[6]));
	reader_2lss1tau_BDT4_ -> AddVariable("mTauTauVis2", &(inputVars_2lss1tau_BDT4_[7]));
	reader_2lss1tau_BDT4_ -> AddVariable("max_lep_eta", &(inputVars_2lss1tau_BDT4_[8]));
	reader_2lss1tau_BDT4_ -> AddVariable("mbb", &(inputVars_2lss1tau_BDT4_[9]));
	reader_2lss1tau_BDT4_ -> AddVariable("mindr_lep1_jet", &(inputVars_2lss1tau_BDT4_[10]));
	reader_2lss1tau_BDT4_ -> AddVariable("mindr_lep2_jet", &(inputVars_2lss1tau_BDT4_[11]));
	reader_2lss1tau_BDT4_ -> AddVariable("mindr_tau_jet", &(inputVars_2lss1tau_BDT4_[12]));
	reader_2lss1tau_BDT4_ -> AddVariable("nJet", &(inputVars_2lss1tau_BDT4_[13]));
	reader_2lss1tau_BDT4_ -> AddVariable("ptmiss", &(inputVars_2lss1tau_BDT4_[14]));
	reader_2lss1tau_BDT4_ -> AddVariable("tau_pt", &(inputVars_2lss1tau_BDT4_[15]));
	reader_2lss1tau_BDT4_ -> AddVariable("HTT", &(inputVars_2lss1tau_BDT4_[16]));
	reader_2lss1tau_BDT4_ -> AddVariable("Hj_tagger", &(inputVars_2lss1tau_BDT4_[17]));
	reader_2lss1tau_BDT4_ -> AddVariable("HadTop_pt", &(inputVars_2lss1tau_BDT4_[18]));

	reader_2lss1tau_BDT4_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_2lss1tau_BDT4(float Vars[19])
{
	for (size_t i=0; i<19; ++i)
		inputVars_2lss1tau_BDT4_[i] = Vars[i];

	return reader_2lss1tau_BDT4_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_2lss1tau_BDT5()
{
	const TString weights =
		data_directory_ + "2lss1tau/2lss_1tau_XGB_HTTMEM_evtLevelSUM_TTH_M_20Var.xml";
	reader_2lss1tau_BDT5_ = new TMVA::Reader("!Color:!Silent");

	reader_2lss1tau_BDT5_ -> AddVariable("avg_dr_jet", &(inputVars_2lss1tau_BDT5_[0]));
	reader_2lss1tau_BDT5_ -> AddVariable("dr_lep1_tau", &(inputVars_2lss1tau_BDT5_[1]));
	reader_2lss1tau_BDT5_ -> AddVariable("dr_lep2_tau", &(inputVars_2lss1tau_BDT5_[2]));
	reader_2lss1tau_BDT5_ -> AddVariable("dr_leps", &(inputVars_2lss1tau_BDT5_[3]));
	reader_2lss1tau_BDT5_ -> AddVariable("lep2_conePt", &(inputVars_2lss1tau_BDT5_[4]));
	reader_2lss1tau_BDT5_ -> AddVariable("mT_lep1", &(inputVars_2lss1tau_BDT5_[5]));
	reader_2lss1tau_BDT5_ -> AddVariable("mT_lep2", &(inputVars_2lss1tau_BDT5_[6]));
	reader_2lss1tau_BDT5_ -> AddVariable("mTauTauVis2", &(inputVars_2lss1tau_BDT5_[7]));
	reader_2lss1tau_BDT5_ -> AddVariable("max_lep_eta", &(inputVars_2lss1tau_BDT5_[8]));
	reader_2lss1tau_BDT5_ -> AddVariable("mbb", &(inputVars_2lss1tau_BDT5_[9]));
	reader_2lss1tau_BDT5_ -> AddVariable("mindr_lep1_jet", &(inputVars_2lss1tau_BDT5_[10]));
	reader_2lss1tau_BDT5_ -> AddVariable("mindr_lep2_jet", &(inputVars_2lss1tau_BDT5_[11]));
	reader_2lss1tau_BDT5_ -> AddVariable("mindr_tau_jet", &(inputVars_2lss1tau_BDT5_[12]));
	reader_2lss1tau_BDT5_ -> AddVariable("nJet", &(inputVars_2lss1tau_BDT5_[13]));
	reader_2lss1tau_BDT5_ -> AddVariable("ptmiss", &(inputVars_2lss1tau_BDT5_[14]));
	reader_2lss1tau_BDT5_ -> AddVariable("tau_pt", &(inputVars_2lss1tau_BDT5_[15]));
	reader_2lss1tau_BDT5_ -> AddVariable("memOutput_LR", &(inputVars_2lss1tau_BDT5_[16]));
	reader_2lss1tau_BDT5_ -> AddVariable("HTT", &(inputVars_2lss1tau_BDT5_[17]));
	reader_2lss1tau_BDT5_ -> AddVariable("Hj_tagger", &(inputVars_2lss1tau_BDT5_[18]));
	reader_2lss1tau_BDT5_ -> AddVariable("HadTop_pt", &(inputVars_2lss1tau_BDT5_[19]));
	
	reader_2lss1tau_BDT5_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_2lss1tau_BDT5(float Vars[20])
{
	for (size_t i=0; i<20; ++i)
		inputVars_2lss1tau_BDT5_[i] = Vars[i];

	return reader_2lss1tau_BDT5_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_3l1tau_BDT1()
{
	const TString weights =
		data_directory_ + "3l1tau/3l_1tau_XGB_plainKin_evtLevelTTV_TTH_13Var.xml";
	reader_3l1tau_BDT1_ = new TMVA::Reader("!Color:!Silent");

	reader_3l1tau_BDT1_ -> AddVariable("lep1_conePt", &(inputVars_3l1tau_BDT1_[0]));
	reader_3l1tau_BDT1_ -> AddVariable("lep2_conePt", &(inputVars_3l1tau_BDT1_[1]));
	reader_3l1tau_BDT1_ -> AddVariable("mindr_lep1_jet", &(inputVars_3l1tau_BDT1_[2]));
	reader_3l1tau_BDT1_ -> AddVariable("mindr_lep2_jet", &(inputVars_3l1tau_BDT1_[3]));
	reader_3l1tau_BDT1_ -> AddVariable("mT_lep2", &(inputVars_3l1tau_BDT1_[4]));
	reader_3l1tau_BDT1_ -> AddVariable("mT_lep1", &(inputVars_3l1tau_BDT1_[5]));
	reader_3l1tau_BDT1_ -> AddVariable("max_lep_eta", &(inputVars_3l1tau_BDT1_[6]));
	reader_3l1tau_BDT1_ -> AddVariable("avg_dr_jet", &(inputVars_3l1tau_BDT1_[7]));
	reader_3l1tau_BDT1_ -> AddVariable("ptmiss", &(inputVars_3l1tau_BDT1_[8]));
	reader_3l1tau_BDT1_ -> AddVariable("tau_pt", &(inputVars_3l1tau_BDT1_[9]));
	reader_3l1tau_BDT1_ -> AddVariable("dr_leps", &(inputVars_3l1tau_BDT1_[10]));
	reader_3l1tau_BDT1_ -> AddVariable("mTauTauVis1", &(inputVars_3l1tau_BDT1_[11]));
	reader_3l1tau_BDT1_ -> AddVariable("mTauTauVis2", &(inputVars_3l1tau_BDT1_[12]));
	
	reader_3l1tau_BDT1_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_3l1tau_BDT1(float Vars[13])
{
	for (size_t i=0; i<13; ++i)
		inputVars_3l1tau_BDT1_[i] = Vars[i];

	return reader_3l1tau_BDT1_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_3l1tau_BDT2()
{
	const TString weights =
		data_directory_ + "3l1tau/3l_1tau_XGB_plainKin_evtLevelTT_TTH_15Var.xml";
	reader_3l1tau_BDT2_ = new TMVA::Reader("!Color:!Silent");

	reader_3l1tau_BDT2_ -> AddVariable("mindr_lep1_jet", &(inputVars_3l1tau_BDT2_[0]));
	reader_3l1tau_BDT2_ -> AddVariable("mindr_lep2_jet", &(inputVars_3l1tau_BDT2_[1]));
	reader_3l1tau_BDT2_ -> AddVariable("mT_lep2", &(inputVars_3l1tau_BDT2_[2]));
	reader_3l1tau_BDT2_ -> AddVariable("mT_lep1", &(inputVars_3l1tau_BDT2_[3]));
	reader_3l1tau_BDT2_ -> AddVariable("max_lep_eta", &(inputVars_3l1tau_BDT2_[4]));
	reader_3l1tau_BDT2_ -> AddVariable("lep3_conePt", &(inputVars_3l1tau_BDT2_[5]));
	reader_3l1tau_BDT2_ -> AddVariable("mindr_lep3_jet", &(inputVars_3l1tau_BDT2_[6]));
	reader_3l1tau_BDT2_ -> AddVariable("mindr_tau_jet", &(inputVars_3l1tau_BDT2_[7]));
	reader_3l1tau_BDT2_ -> AddVariable("avg_dr_jet", &(inputVars_3l1tau_BDT2_[8]));
	reader_3l1tau_BDT2_ -> AddVariable("ptmiss", &(inputVars_3l1tau_BDT2_[9]));
	reader_3l1tau_BDT2_ -> AddVariable("tau_pt", &(inputVars_3l1tau_BDT2_[10]));
	reader_3l1tau_BDT2_ -> AddVariable("dr_leps", &(inputVars_3l1tau_BDT2_[11]));
	reader_3l1tau_BDT2_ -> AddVariable("mTauTauVis1", &(inputVars_3l1tau_BDT2_[12]));
	reader_3l1tau_BDT2_ -> AddVariable("mTauTauVis2", &(inputVars_3l1tau_BDT2_[13]));
	reader_3l1tau_BDT2_ -> AddVariable("mbb_loose", &(inputVars_3l1tau_BDT2_[14]));
	
	reader_3l1tau_BDT2_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_3l1tau_BDT2(float Vars[15])
{
	for (size_t i=0; i<15; ++i)
		inputVars_3l1tau_BDT2_[i] = Vars[i];

	return reader_3l1tau_BDT2_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_3l1tau_BDT3()
{
	const TString weights =
		data_directory_ + "3l1tau/3l_1tau_XGB_plainKin_evtLevelSUM_TTH_M_12Var.xml";
	reader_3l1tau_BDT3_ = new TMVA::Reader("!Color:!Silent");

	reader_3l1tau_BDT3_ -> AddVariable("lep1_conePt", &(inputVars_3l1tau_BDT3_[0]));
	reader_3l1tau_BDT3_ -> AddVariable("lep2_conePt", &(inputVars_3l1tau_BDT3_[1]));
	reader_3l1tau_BDT3_ -> AddVariable("mindr_lep1_jet", &(inputVars_3l1tau_BDT3_[2]));
	reader_3l1tau_BDT3_ -> AddVariable("max_lep_eta", &(inputVars_3l1tau_BDT3_[3]));
	reader_3l1tau_BDT3_ -> AddVariable("mindr_tau_jet", &(inputVars_3l1tau_BDT3_[4]));
	reader_3l1tau_BDT3_ -> AddVariable("ptmiss", &(inputVars_3l1tau_BDT3_[5]));
	reader_3l1tau_BDT3_ -> AddVariable("tau_pt", &(inputVars_3l1tau_BDT3_[6]));
	reader_3l1tau_BDT3_ -> AddVariable("dr_leps", &(inputVars_3l1tau_BDT3_[7]));
	reader_3l1tau_BDT3_ -> AddVariable("mTauTauVis1", &(inputVars_3l1tau_BDT3_[8]));
	reader_3l1tau_BDT3_ -> AddVariable("mTauTauVis2", &(inputVars_3l1tau_BDT3_[9]));
	reader_3l1tau_BDT3_ -> AddVariable("mbb_loose", &(inputVars_3l1tau_BDT3_[10]));
	reader_3l1tau_BDT3_ -> AddVariable("nJet", &(inputVars_3l1tau_BDT3_[11]));

	reader_3l1tau_BDT3_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_3l1tau_BDT3(float Vars[12])
{
	for (size_t i=0; i<12; ++i)
		inputVars_3l1tau_BDT3_[i] = Vars[i];

	return reader_3l1tau_BDT3_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_2l2tau_BDT1()
{
	const TString weights =
		data_directory_ + "2l2tau/2l_2tau_XGB_plainKin_evtLevelTTV_TTH_14Var.xml";
	reader_2l2tau_BDT1_ = new TMVA::Reader("!Color:!Silent");

	reader_2l2tau_BDT1_ -> AddVariable("mTauTauVis", &(inputVars_2l2tau_BDT1_[0]));
	reader_2l2tau_BDT1_ -> AddVariable("cosThetaS_hadTau", &(inputVars_2l2tau_BDT1_[1]));
	reader_2l2tau_BDT1_ -> AddVariable("lep1_conePt", &(inputVars_2l2tau_BDT1_[2]));
	reader_2l2tau_BDT1_ -> AddVariable("lep2_conePt", &(inputVars_2l2tau_BDT1_[3]));
	reader_2l2tau_BDT1_ -> AddVariable("mT_lep1", &(inputVars_2l2tau_BDT1_[4]));
	reader_2l2tau_BDT1_ -> AddVariable("mT_lep2", &(inputVars_2l2tau_BDT1_[5]));
	reader_2l2tau_BDT1_ -> AddVariable("dr_taus", &(inputVars_2l2tau_BDT1_[6]));
	reader_2l2tau_BDT1_ -> AddVariable("min_dr_lep_jet", &(inputVars_2l2tau_BDT1_[7]));
	reader_2l2tau_BDT1_ -> AddVariable("mindr_tau1_jet", &(inputVars_2l2tau_BDT1_[8]));
	reader_2l2tau_BDT1_ -> AddVariable("avg_dr_jet", &(inputVars_2l2tau_BDT1_[9]));
	reader_2l2tau_BDT1_ -> AddVariable("min_dr_lep_tau", &(inputVars_2l2tau_BDT1_[10]));
	reader_2l2tau_BDT1_ -> AddVariable("max_dr_lep_tau", &(inputVars_2l2tau_BDT1_[11]));
	reader_2l2tau_BDT1_ -> AddVariable("is_OS", &(inputVars_2l2tau_BDT1_[12]));
	reader_2l2tau_BDT1_ -> AddVariable("nJet", &(inputVars_2l2tau_BDT1_[13]));

	reader_2l2tau_BDT1_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_2l2tau_BDT1(float Vars[14])
{
	for (size_t i=0; i<14; ++i)
		inputVars_2l2tau_BDT1_[i] = Vars[i];

	return reader_2l2tau_BDT1_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_2l2tau_BDT2()
{
	const TString weights =
		data_directory_ + "2l2tau/2l_2tau_XGB_plainKin_evtLevelTT_TTH_11Var.xml";
	reader_2l2tau_BDT2_ = new TMVA::Reader("!Color:!Silent");

	reader_2l2tau_BDT2_ -> AddVariable("mTauTauVis", &(inputVars_2l2tau_BDT2_[0]));
	reader_2l2tau_BDT2_ -> AddVariable("cosThetaS_hadTau", &(inputVars_2l2tau_BDT2_[1]));
	reader_2l2tau_BDT2_ -> AddVariable("tau1_pt", &(inputVars_2l2tau_BDT2_[2]));
	reader_2l2tau_BDT2_ -> AddVariable("tau2_pt", &(inputVars_2l2tau_BDT2_[3]));
	reader_2l2tau_BDT2_ -> AddVariable("tau2_eta", &(inputVars_2l2tau_BDT2_[4]));
	reader_2l2tau_BDT2_ -> AddVariable("mindr_lep1_jet", &(inputVars_2l2tau_BDT2_[5]));
	reader_2l2tau_BDT2_ -> AddVariable("mT_lep1", &(inputVars_2l2tau_BDT2_[6]));
	reader_2l2tau_BDT2_ -> AddVariable("mindr_tau_jet", &(inputVars_2l2tau_BDT2_[7]));
	reader_2l2tau_BDT2_ -> AddVariable("max_dr_lep_tau", &(inputVars_2l2tau_BDT2_[8]));
	reader_2l2tau_BDT2_ -> AddVariable("is_OS", &(inputVars_2l2tau_BDT2_[9]));
	reader_2l2tau_BDT2_ -> AddVariable("nBJetLoose", &(inputVars_2l2tau_BDT2_[10]));

	reader_2l2tau_BDT2_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_2l2tau_BDT2(float Vars[11])
{
	for (size_t i=0; i<11; ++i)
		inputVars_2l2tau_BDT2_[i] = Vars[i];

	return reader_2l2tau_BDT2_ -> EvaluateMVA("BDT");
}

void MVAEvaluator::setup_tmva_reader_2l2tau_BDT3()
{
	const TString weights =
		data_directory_ + "2l2tau/2l_2tau_XGB_plainKin_evtLevelSUM_TTH_VT_13Var.xml";
	reader_2l2tau_BDT3_ = new TMVA::Reader("!Color:!Silent");

	reader_2l2tau_BDT3_ -> AddVariable("mTauTauVis", &(inputVars_2l2tau_BDT3_[0]));
	reader_2l2tau_BDT3_ -> AddVariable("cosThetaS_hadTau", &(inputVars_2l2tau_BDT3_[1]));
	reader_2l2tau_BDT3_ -> AddVariable("tau1_pt", &(inputVars_2l2tau_BDT3_[2]));
	reader_2l2tau_BDT3_ -> AddVariable("tau2_pt", &(inputVars_2l2tau_BDT3_[3]));
	reader_2l2tau_BDT3_ -> AddVariable("lep2_conePt", &(inputVars_2l2tau_BDT3_[4]));
	reader_2l2tau_BDT3_ -> AddVariable("mindr_lep1_jet", &(inputVars_2l2tau_BDT3_[5]));
	reader_2l2tau_BDT3_ -> AddVariable("mT_lep1", &(inputVars_2l2tau_BDT3_[6]));
	reader_2l2tau_BDT3_ -> AddVariable("mindr_tau_jet", &(inputVars_2l2tau_BDT3_[7]));
	reader_2l2tau_BDT3_ -> AddVariable("avg_dr_jet", &(inputVars_2l2tau_BDT3_[8]));
	reader_2l2tau_BDT3_ -> AddVariable("avr_dr_lep_tau", &(inputVars_2l2tau_BDT3_[9]));
	reader_2l2tau_BDT3_ -> AddVariable("dr_taus", &(inputVars_2l2tau_BDT3_[10]));
	reader_2l2tau_BDT3_ -> AddVariable("is_OS", &(inputVars_2l2tau_BDT3_[11]));
	reader_2l2tau_BDT3_ -> AddVariable("nBJetLoose", &(inputVars_2l2tau_BDT3_[12]));

	reader_2l2tau_BDT3_ -> BookMVA("BDT", weights);
}

float MVAEvaluator::evaluate_bdt_2l2tau_BDT3(float Vars[13])
{
	for (size_t i=0; i<13; ++i)
		inputVars_2l2tau_BDT3_[i] = Vars[i];

	return reader_2l2tau_BDT3_ -> EvaluateMVA("BDT");
}
