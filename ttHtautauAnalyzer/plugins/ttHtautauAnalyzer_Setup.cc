#ifndef ttHtautauAnalyzer_Setup_cc
#define ttHtautauAnalyzer_Setup_cc

#include "ttHTauTauAnalysis/ttHtautauAnalyzer/plugins/ttHtautauAnalyzer.h"

void ttHtautauAnalyzer::Set_up_AnaType(const std::string& anatype)
{
	if (anatype == "1l2tau") {
		anaType_ =  Analysis_types::Analyze_1l2tau;
	}
	else if (anatype == "2lss1tau") {
		anaType_ = Analysis_types::Analyze_2lss1tau;
	}
	else if (anatype == "3l1tau") {
		anaType_ = Analysis_types::Analyze_3l;
	}

	return;
}

void ttHtautauAnalyzer::Set_up_SelType(const std::string& seltype)
{
	if (seltype == "signal_2lss1tau") {
		selType_ = Selection_types::Signal_2lss1tau;
	}
	else if (seltype == "signal_1l2tau") {
		selType_ = Selection_types::Signal_1l2tau;
	}
	else if (seltype == "signal_3l") {
		selType_ = Selection_types::Signal_3l;
	}
	else if (seltype == "control_2los1tau") {
		selType_ = Selection_types::Control_2los1tau;
	}
	else if (seltype == "control_1lfakeable") {
		selType_ = Selection_types::Control_1lfakeable;
	}
	else {
		std::cerr << "Not valid selection region!" << std::endl;
		assert(0);
	}

	return;
}

void ttHtautauAnalyzer::Set_up_histograms()
{
	h_nProcessed_ = fs_->make<TH1I>("h_nProcessed","",1,0.5,1.5);
	h_SumGenWeight_ = fs_->make<TH1D>("h_SumGenWeight","",1,0.5,1.5);
	h_SumGenWeightxPU_ = fs_->make<TH1D>("h_SumGenWeightxPU","",1,0.5,1.5);

	// cut flow
	h_CutFlow_ = fs_->make<TH1D>("h_CutFlow","",15,0,15);
}

void ttHtautauAnalyzer::Set_up_tree()
{
	eventTree_ = fs_->make<TTree>("eventTree", "Event tree");

	evNtuple_.setup_branches(eventTree_);
}

#endif
