#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/MVAVars.h"

// constructor
MVAVars::MVAVars(const std::vector<miniLepton>& leptons, const std::vector<TLorentzVector>& taus, const std::vector<TLorentzVector>& jets, float met, float metphi, float mht)
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
