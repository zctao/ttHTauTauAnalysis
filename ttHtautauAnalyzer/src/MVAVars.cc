#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/MVAVars.h"

// constructor
MVAVars::MVAVars(const std::vector<miniLepton>& leptons, const std::vector<TLorentzVector>& taus, const std::vector<TLorentzVector>& jets, float met, float metphi, float mht)
{
	compute_all_variables(leptons, taus, jets, met, metphi, mht);
}

void MVAVars::compute_all_variables(const std::vector<miniLepton>& leptons, const std::vector<TLorentzVector>& taus, const std::vector<TLorentzVector>& jets, float met, float metphi, float mht)
{
	assert(leptons.size()>0 and taus.size()>0);
	
	nJet_ = jets.size();

	mindr_lep0_jet_ =compute_mindr(leptons[0].p4(),jets);
	mindr_lep1_jet_ = (leptons.size()>1) ? compute_mindr(leptons[1].p4(),jets) : -9999;

	avg_dr_jet_ = compute_average_dr(jets);

	max_lep_eta_ = (leptons.size()>1) ? std::max(fabs(leptons[0].eta()),fabs(leptons[1].eta())) : -9999.;

	met_ = met;
	mht_ = mht;

	mT_met_lep0_ =
		sqrt( 2*leptons[0].conept()*met*(1-cos(leptons[0].phi()-metphi)) );

	lep0_conept_ = leptons[0].conept();
	lep1_conept_ = (leptons.size()>1) ? leptons[1].conept() : -9999.;

	dr_leps_ = (leptons.size()>1) ? leptons[0].p4().DeltaR(leptons[1].p4()) : -9999.;

	tau_pt_ = taus[0].Pt();

	dr_lep0_tau_ = leptons[0].p4().DeltaR(taus[0]);
	dr_lep1_tau_ = (leptons.size()>1) ? leptons[1].p4().DeltaR(taus[0]) : -9999.;

	mvis_lep0_tau_ = (leptons[0].p4()+taus[0]).M();
	mvis_lep1_tau_ = (leptons.size()>1) ? (leptons[1].p4()+taus[0]).M() : -9999.;
}

float MVAVars::compute_mindr(const TLorentzVector& l, const std::vector<TLorentzVector>& vjs)
{
	if (vjs.size() < 1) return -9999.;

	float mindr = 666.;

	for (auto & j : vjs) {
		float dr = l.DeltaR(j);
		if (dr < mindr) mindr = dr;
	}

	return mindr;	
}

float MVAVars::compute_average_dr(const std::vector<TLorentzVector>& lvs)
{
	if (lvs.size()<2) return -9999.;

	float sum_dr = 0.;
	int ncomb = 0;

	for (auto i = lvs.begin(); i != lvs.end()-1; ++i) {
		for (auto j = i+1; j != lvs.end(); ++j) {
			++ncomb;
			sum_dr += i->DeltaR(*j);
		}
	}

	return sum_dr / ncomb;
}

void MVAVars::set_up_tmva_reader()
{
	/////////
	// ttV
	const TString weights_ttV = (std::string(getenv("CMSSW_BASE"))+"/src/ttHTauTauAnalysis/ttHtautauAnalyzer/data/tthMVA/2lss1tau_v2/mvaTTHvsTTV2lss1tau_BDTG.weights.xml").c_str();
	reader_2lss1tau_ttV = new TMVA::Reader("!Color:!Silent");

	reader_2lss1tau_ttV->AddVariable("mindr_lep1_jet", &mindr_lep0_jet_);
	reader_2lss1tau_ttV->AddVariable("mindr_lep2_jet", &mindr_lep1_jet_);
	reader_2lss1tau_ttV->AddVariable("avg_dr_jet", &avg_dr_jet_);
	reader_2lss1tau_ttV->AddVariable("TMath::Max(TMath::Abs(lep1_eta), TMath::Abs(lep2_eta))", &max_lep_eta_);
	//reader_2lss1tau_ttV->AddVariable("max_lep_eta", &max_lep_eta_);
	reader_2lss1tau_ttV->AddVariable("lep1_conePt", &lep0_conept_);
	reader_2lss1tau_ttV->AddVariable("lep2_conePt", &lep1_conept_);
	reader_2lss1tau_ttV->AddVariable("mT_lep1", &mT_met_lep0_);
	reader_2lss1tau_ttV->AddVariable("dr_leps", &dr_leps_);
	reader_2lss1tau_ttV->AddVariable("mTauTauVis1", &mvis_lep0_tau_);
	reader_2lss1tau_ttV->AddVariable("mTauTauVis2", &mvis_lep1_tau_);
	
	reader_2lss1tau_ttV->BookMVA("BDT", weights_ttV);

	/////////
	// ttbar
	const TString weights_ttbar = (std::string(getenv("CMSSW_BASE"))+"/src/ttHTauTauAnalysis/ttHtautauAnalyzer/data/tthMVA/2lss1tau_v2/mvaTTHvsTTbar2lss1tau_BDTG.weights.xml").c_str();
	reader_2lss1tau_ttbar = new TMVA::Reader("!Color:!Silent");

	reader_2lss1tau_ttbar->AddVariable("nJet", &nJet_);
	reader_2lss1tau_ttbar->AddVariable("mindr_lep1_jet", &mindr_lep0_jet_);
	//reader_2lss1tau_ttbar->AddVariable("mindr_lep2_jet", &mindr_lep1_jet_);
	reader_2lss1tau_ttbar->AddVariable("avg_dr_jet", &avg_dr_jet_);
	reader_2lss1tau_ttbar->AddVariable("TMath::Max(TMath::Abs(lep1_eta), TMath::Abs(lep2_eta))", &max_lep_eta_);
	//reader_2lss1tau_ttbar->AddVariable("max_lep_eta", &max_lep_eta_);
	//reader_2lss1tau_ttbar->AddVariable("lep1_conePt", &lep0_conept_);
	reader_2lss1tau_ttbar->AddVariable("lep2_conePt", &lep1_conept_);
	//reader_2lss1tau_ttbar->AddVariable("ptmiss", &met_);
	//reader_2lss1tau_ttbar->AddVariable("htmiss", &mht_);
	reader_2lss1tau_ttbar->AddVariable("dr_leps", &dr_leps_);
	reader_2lss1tau_ttbar->AddVariable("tau_pt", &tau_pt_);
	reader_2lss1tau_ttbar->AddVariable("dr_lep1_tau", &dr_lep0_tau_);
	
	reader_2lss1tau_ttbar->BookMVA("BDT", weights_ttbar);
}

float MVAVars::BDT_2lss1tau_ttV()
{
	return reader_2lss1tau_ttV->EvaluateMVA("BDT");
}

float MVAVars::BDT_2lss1tau_ttbar()
{
	return reader_2lss1tau_ttbar->EvaluateMVA("BDT");
}
