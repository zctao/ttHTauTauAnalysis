#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/MVAVars.h"

// constructor
//MVAVars::MVAVars(const std::vector<miniLepton>& leptons, const std::vector<TLorentzVector>& taus, const std::vector<TLorentzVector>& jets, float met, float metphi, float mht)
//{
//	compute_all_variables(leptons, taus, jets, met, metphi, mht);
//}

void MVAVars::compute_all_variables(const std::vector<miniLepton>& leptons, const std::vector<TLorentzVector>& taus, const std::vector<TLorentzVector>& jets, float met, float metphi, float mht, int ntags_loose)
{
	assert(leptons.size()>0 and taus.size()>0);

	nJet_ = jets.size();
	avg_dr_jet_ = compute_average_dr(jets);
	
	if (anaType_==Analyze_2lss1tau) {
	
		mindr_lep0_jet_ =compute_mindr(leptons[0].p4(),jets);
		mindr_lep1_jet_ = (leptons.size()>1) ? compute_mindr(leptons[1].p4(),jets) : -9999;
		
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
	else if (anaType_==Analyze_1l2tau) {
		assert(leptons.size()>0 and taus.size()>1);
		
		ht_ = mht; // FIXME

		dr_taus_ = taus[0].DeltaR(taus[1]);

		mvis_taus_ = (taus[0]+taus[1]).M();
		pt_taus_ = (taus[0]+taus[1]).Pt();

		max_dr_jet_ = compute_max_dr(jets);;

		tau0_pt_ = taus[0].Pt();
		tau1_pt_ = taus[1].Pt();

		ntags_loose_ = ntags_loose;
	}
	else if (anaType_==Analyze_3l1tau) {
		assert(leptons.size()>2 and taus.size()>0);
		
		max_lep_eta_ = std::max(std::max(fabs(leptons[0].eta()),fabs(leptons[1].eta())), fabs(leptons[2].eta()));
		
		mindr_lep0_jet_ = compute_mindr(leptons[0].p4(),jets);	
		mindr_lep1_jet_ = compute_mindr(leptons[1].p4(),jets);
		mindr_lep2_jet_ = compute_mindr(leptons[2].p4(),jets);
		
		mT_met_lep0_ =
			sqrt( 2*leptons[0].conept()*met*(1-cos(leptons[0].phi()-metphi)) );;
		
		lep0_conept_ = leptons[0].conept();
		lep1_conept_ = leptons[1].conept();
		lep2_conept_ = leptons[2].conept();

		mht_ = mht;
	}
}

void MVAVars::compute_taudecay_variables(const std::vector<TLorentzVector>& taus,
										 const std::vector<TLorentzVector>& chargedDaug,
										 const std::vector<TLorentzVector>& neutralDaug,
										 const std::vector<int> decaymode)
{
	assert(taus.size()>0);
	tau0_decaymode_ = decaymode[0];
	tau0_energy_ = taus[0].E();
	tau0_upsilon_ = (chargedDaug[0].E()-neutralDaug[0].E()) / (chargedDaug[0].E()+neutralDaug[0].E());
	
	if (anaType_==Analyze_1l2tau) {
		assert(taus.size()>1);
		tau1_decaymode_ = decaymode[1];
		tau1_energy_ = taus[1].E();
		tau1_upsilon_ = (chargedDaug[1].E()-neutralDaug[1].E()) / (chargedDaug[1].E()+neutralDaug[1].E());
	}
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

float MVAVars::compute_max_dr(const std::vector<TLorentzVector>& lvs)
{
	float max_dr = 0.;

	for (auto i = lvs.begin(); i != lvs.end()-1; ++i) {
		for (auto j = i+1; j != lvs.end(); ++j) {
		    float dr = i->DeltaR(*j);
			if (dr > max_dr) max_dr = dr;
		}
	}

	return max_dr;
}

void MVAVars::set_up_tmva_reader()
{
	const TString datadir = (std::string(getenv("CMSSW_BASE"))+"/src/ttHTauTauAnalysis/ttHtautauAnalyzer/data/").c_str();
	
	/////////////
	// 2lss1tau
	if (anaType_==Analyze_2lss1tau) {
		/// ttV
		const TString weights_ttV = datadir + "tthMVA/2lss1tau_v2/mvaTTHvsTTV2lss1tau_BDTG.weights.xml";
		reader_ttV = new TMVA::Reader("!Color:!Silent");

		reader_ttV->AddVariable("mindr_lep1_jet", &mindr_lep0_jet_);
		reader_ttV->AddVariable("mindr_lep2_jet", &mindr_lep1_jet_);
		reader_ttV->AddVariable("avg_dr_jet", &avg_dr_jet_);
		reader_ttV->AddVariable("TMath::Max(TMath::Abs(lep1_eta), TMath::Abs(lep2_eta))", &max_lep_eta_);
		//reader_ttV->AddVariable("max_lep_eta", &max_lep_eta_);
		reader_ttV->AddVariable("lep1_conePt", &lep0_conept_);
		reader_ttV->AddVariable("lep2_conePt", &lep1_conept_);
		reader_ttV->AddVariable("mT_lep1", &mT_met_lep0_);
		reader_ttV->AddVariable("dr_leps", &dr_leps_);
		reader_ttV->AddVariable("mTauTauVis1", &mvis_lep0_tau_);
		reader_ttV->AddVariable("mTauTauVis2", &mvis_lep1_tau_);
	
		reader_ttV->BookMVA("BDT", weights_ttV);
		
		//// ttbar
		const TString weights_ttbar = datadir + "tthMVA/2lss1tau_v2/mvaTTHvsTTbar2lss1tau_BDTG.weights.xml";
		reader_ttbar = new TMVA::Reader("!Color:!Silent");
		
		reader_ttbar->AddVariable("nJet", &nJet_);
		reader_ttbar->AddVariable("mindr_lep1_jet", &mindr_lep0_jet_);
		//reader_ttbar->AddVariable("mindr_lep2_jet", &mindr_lep1_jet_);
		reader_ttbar->AddVariable("avg_dr_jet", &avg_dr_jet_);
		reader_ttbar->AddVariable("TMath::Max(TMath::Abs(lep1_eta), TMath::Abs(lep2_eta))", &max_lep_eta_);
		//reader_ttbar->AddVariable("max_lep_eta", &max_lep_eta_);
		//reader_ttbar->AddVariable("lep1_conePt", &lep0_conept_);
		reader_ttbar->AddVariable("lep2_conePt", &lep1_conept_);
		//reader_ttbar->AddVariable("ptmiss", &met_);
		//reader_ttbar->AddVariable("htmiss", &mht_);
		reader_ttbar->AddVariable("dr_leps", &dr_leps_);
		reader_ttbar->AddVariable("tau_pt", &tau_pt_);
		reader_ttbar->AddVariable("dr_lep1_tau", &dr_lep0_tau_);
	
		reader_ttbar->BookMVA("BDT", weights_ttbar);
	}
	else if (anaType_==Analyze_1l2tau) {
		/// ttV
		const TString weights_ttV = datadir + "evtMVAWeights/1l2tau/sklearn_ttZ.xml";
		reader_ttV = new TMVA::Reader("!Color:!Silent");

		reader_ttV->AddVariable("ht", &ht_);
		reader_ttV->AddVariable("tt_deltaR", &dr_taus_);
		reader_ttV->AddVariable("tt_visiblemass", &mvis_taus_);
		reader_ttV->AddVariable("tt_sumpt", &pt_taus_);
		reader_ttV->AddVariable("jet_deltaRmax", &max_dr_jet_);
		reader_ttV->AddVariable("jet_deltaRavg", &avg_dr_jet_);
		reader_ttV->AddVariable("njets_inclusive", &nJet_);

		reader_ttV->BookMVA("BDT", weights_ttV);
		
		/// ttbar
		const TString weights_ttbar = datadir + "evtMVAWeights/1l2tau/sklearn_tt.xml";
		reader_ttbar = new TMVA::Reader("!Color:!Silent");

		reader_ttbar->AddVariable("ht", &ht_);
		reader_ttbar->AddVariable("tt_deltaR", &dr_taus_);
		reader_ttbar->AddVariable("tt_visiblemass", &mvis_taus_);
		reader_ttbar->AddVariable("tau1_pt", &tau0_pt_);
		reader_ttbar->AddVariable("tau2_pt", &tau1_pt_);
		reader_ttbar->AddVariable("jet_deltaRavg", &avg_dr_jet_);
		reader_ttbar->AddVariable("njets_inclusive", &nJet_);
		reader_ttbar->AddVariable("ntags_loose", &ntags_loose_);

		reader_ttbar->BookMVA("BDT", weights_ttbar);
	}
	else if (anaType_==Analyze_3l1tau) {
		/// ttV
		const TString weights_ttV = datadir + "evtMVAWeights/3l1tau/3l_ttV_BDTG.weights.xml";
		reader_ttV = new TMVA::Reader("!Color:!Silent");

		reader_ttV->AddVariable("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))", &max_lep_eta_);
		reader_ttV->AddVariable("nJet25_Recl", &nJet_);
		reader_ttV->AddVariable("mindr_lep1_jet", &mindr_lep0_jet_);
		reader_ttV->AddVariable("mindr_lep2_jet", &mindr_lep1_jet_);
		reader_ttV->AddVariable("MT_met_lep1", &mT_met_lep0_);
		reader_ttV->AddVariable("LepGood_conePt[iF_Recl[2]]", &lep2_conept_);
		reader_ttV->AddVariable("LepGood_conePt[iF_Recl[0]]", &lep0_conept_);
		// Spectators
		float iF0 = 0;
		float iF1 = 1;
		float iF2 = 2;
		reader_ttV->AddSpectator("iF_Recl[0]", &iF0);
		reader_ttV->AddSpectator("iF_Recl[1]", &iF1);
		reader_ttV->AddSpectator("iF_Recl[2]", &iF2);	

		reader_ttV->BookMVA("BDT", weights_ttV);
		
		/// ttbar
		const TString weights_ttbar = datadir + "evtMVAWeights/3l1tau/3l_ttbar_BDTG.weights.xml";
		reader_ttbar = new TMVA::Reader("!Color:!Silent");
		
		reader_ttbar->AddVariable("max(abs(LepGood_eta[iF_Recl[0]]),abs(LepGood_eta[iF_Recl[1]]))", &max_lep_eta_);
		reader_ttbar->AddVariable("nJet25_Recl", &nJet_);
		reader_ttbar->AddVariable("mindr_lep1_jet", &mindr_lep0_jet_);
		reader_ttbar->AddVariable("mindr_lep2_jet", &mindr_lep1_jet_);
		reader_ttbar->AddVariable("MT_met_lep1", &mT_met_lep0_);
		reader_ttbar->AddVariable("mhtJet25_Recl", &mht_);
		reader_ttbar->AddVariable("avg_dr_jet", &avg_dr_jet_);
		// Spectators
		float iF3 = 0;
		float iF4 = 1;
		float iF5 = 2;
		reader_ttbar->AddSpectator("iF_Recl[0]", &iF3);
		reader_ttbar->AddSpectator("iF_Recl[1]", &iF4);
		reader_ttbar->AddSpectator("iF_Recl[2]", &iF5);
		
		reader_ttbar->BookMVA("BDT", weights_ttbar);
	}
	else {
		std::cerr << "Analysis type is not supported!" << std::endl;
		assert(0);
	}
}

float MVAVars::BDT_ttV()
{
	return reader_ttV->EvaluateMVA("BDT");
}

float MVAVars::BDT_ttbar()
{
	return reader_ttbar->EvaluateMVA("BDT");
}
