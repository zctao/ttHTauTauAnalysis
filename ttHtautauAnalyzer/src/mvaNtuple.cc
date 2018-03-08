#ifndef mvaNtuple_cc
#define mvaNtuple_cc

#include "../interface/mvaNtuple.h"

void mvaNtuple::setup_branches(TTree* tree)
{
	tree->Branch("run", &run);
	tree->Branch("lumi", &lumi);
	tree->Branch("nEvent", &nEvent);
	
	tree->Branch("event_weight", &event_weight);
	tree->Branch("xsection_weight", &xsection_weight);
	tree->Branch("xsection_weight_gen", &xsection_weight);
		
	tree->Branch("isGenMatchedTau", &isGenMatchedTau);
	tree->Branch("HiggsDecayType", &HiggsDecayType);
	
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
		if (version_=="2016") { // 2016 analysis variables
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
			tree->Branch("tau0_easym", &tau0_easym);
		}
		else {
			std::cout << "WARNING: Version " << version_
					  << " for 2lss1tau category is unknown. Ntuple will be empty."
					  << std::endl;
		}
	}
	else if (anatype_ == Analyze_1l2tau) {
		if (version_=="2016") { // 2016 analysis variables
			tree->Branch("ht", &HT);
			tree->Branch("tt_deltaR", &dr_taus);
			tree->Branch("tt_mvis", &mTauTauVis);
			tree->Branch("tt_sumpt", &tt_pt);
			tree->Branch("max_dr_jet", &max_dr_jet);
			tree->Branch("tau0_pt", &tau0_pt);
			tree->Branch("tau1_pt", &tau1_pt);
			tree->Branch("nbtags", &nbtags_medium);
			tree->Branch("nbtags_loose", &nbtags_loose);
		}
		else if (version_=="2017") {
			//tree->Branch("avg_dr_jet", &avg_dr_jet);
			tree->Branch("costS_tau", &costS_tau); 
			tree->Branch("dr_taus", &dr_taus);
			tree->Branch("dr_lep_tau_lead", &dr_lep_tau_lead); 
			tree->Branch("dr_lep_tau_sublead", &dr_lep_tau_sublead);
			tree->Branch("dr_lep_tau_ss", &dr_lep_tau_ss);
			tree->Branch("lep_conePt", &lep0_conept);
			tree->Branch("mindr_lep_jet", &mindr_lep0_jet);
			tree->Branch("mindr_tau1_jet", &mindr_tau0_jet);
			tree->Branch("mindr_tau2_jet", &mindr_tau1_jet);
			tree->Branch("mTauTauVis", &mTauTauVis);
			tree->Branch("mT_lep", &mT_met_lep0);
			tree->Branch("ptmiss", &met);
			tree->Branch("ht", &HT);
			tree->Branch("tau1_pt", &tau0_pt);
			tree->Branch("tau2_pt", &tau1_pt);
			tree->Branch("nBjetLoose", &nbtags_loose);
			tree->Branch("nBjetMedium", &nbtags_medium);
			
			tree->Branch("taup_decaymode", &taup_decaymode);
			tree->Branch("taum_decaymode", &taum_decaymode);
			tree->Branch("taup_easym", &taup_easym);
			tree->Branch("taum_easym", &taum_easym);
			tree->Branch("taup_cosPsi", &taup_cosPsi);
			tree->Branch("taum_cosPsi", &taum_cosPsi);
			tree->Branch("evisTaus_diff", &evisTaus_diff);
			tree->Branch("evisTaus_sum", &evisTaus_sum);
			tree->Branch("evisTausAsym", &evisTausAsym);
		}
		else if (version_=="test") {
			tree->Branch("taup_decaymode", &taup_decaymode);
			tree->Branch("taum_decaymode", &taum_decaymode);
			tree->Branch("taup_easym", &taup_easym);
			tree->Branch("taum_easym", &taum_easym);
			tree->Branch("taup_cosPsi", &taup_cosPsi);
			tree->Branch("taum_cosPsi", &taum_cosPsi);
			tree->Branch("evisTaus_diff", &evisTaus_diff);
			tree->Branch("evisTaus_sum", &evisTaus_sum);
			tree->Branch("evisTausAsym", &evisTausAsym);
			
			tree->Branch("pp1_pp1", &pp1_pp1);
			tree->Branch("pp1_pp2", &pp1_pp2);
			tree->Branch("pp1_pp3", &pp1_pp3);
			tree->Branch("pp1_pm1", &pp1_pm1);
			tree->Branch("pp1_pm2", &pp1_pm2);
			tree->Branch("pp1_pm3", &pp1_pm3);
			tree->Branch("pp2_pp2", &pp2_pp2);
			tree->Branch("pp2_pp3", &pp2_pp3);
			tree->Branch("pp2_pm1", &pp2_pm1);
			tree->Branch("pp2_pm2", &pp2_pm2);
			tree->Branch("pp2_pm3", &pp2_pm3);
			tree->Branch("pp3_pp3", &pp3_pp3);
			tree->Branch("pp3_pm1", &pp3_pm1);
			tree->Branch("pp3_pm2", &pp3_pm2);
			tree->Branch("pp3_pm3", &pp3_pm3);
			tree->Branch("pm1_pm1", &pm1_pm1);
			tree->Branch("pm1_pm2", &pm1_pm2);
			tree->Branch("pm1_pm3", &pm1_pm3);
			tree->Branch("pm2_pm2", &pm2_pm2);
			tree->Branch("pm2_pm3", &pm2_pm3);
			tree->Branch("pm3_pm3", &pm3_pm3);
		}
		else {
			std::cout << "WARNING: Version " << version_
					  << " for 1l2tau category is unknown. Ntuple will be empty."
					  << std::endl;
		}
	}
	else if (anatype_ == Analyze_3l1tau) {
		if (version_=="2016") {
			tree->Branch("max_lep_eta", &max_lep_eta);
			tree->Branch("mindr_lep0_jet", &mindr_lep0_jet);
			tree->Branch("mindr_lep1_jet", &mindr_lep1_jet);
			tree->Branch("mindr_lep2_jet", &mindr_lep2_jet);
			tree->Branch("mT_met_lep0", &mT_met_lep0);
			tree->Branch("lep0_conept", &lep0_conept);
			tree->Branch("lep1_conept", &lep1_conept);
			tree->Branch("lep2_conept", &lep2_conept);
			//tree->Branch("tau0_decaymode", &tau0_decaymode);
			//tree->Branch("tau0_E", &tau0_E);
			//tree->Branch("tau0_easym", &tau0_easym);
		}
		else {
			std::cout << "WARNING: Version " << version_
					  << " for 3l1tau category is unknown. Ntuple will be empty."
					  << std::endl;
		}
	}
}

void mvaNtuple::compute_variables(const std::vector<miniLepton>& leptons,
								  const std::vector<miniTau>& taus,
								  const std::vector<TLorentzVector>& jets,
								  float MET, float METphi, float MHT,
								  int nbtagsloose, int nbtagsmedium)
{
	// common variables in all categories
	nJet = jets.size();
	avg_dr_jet = compute_average_dr(jets);
	
	switch(anatype_){
	case Analyze_1l2tau:
		assert(leptons.size()>0);
		assert(taus.size()>1);
		
		if (version_=="2016") {
			HT = mht;  // FIXME
			dr_taus = taus[0].p4().DeltaR(taus[1].p4());;
			mTauTauVis = (taus[0].p4()+taus[1].p4()).M();
			tt_pt = (taus[0].p4()+taus[1].p4()).Pt();
			max_dr_jet = compute_max_dr(jets);
			tau0_pt = taus[0].pt();
			tau1_pt = taus[1].pt();
			nbtags_medium = nbtagsloose;
			nbtags_loose = nbtagsmedium;
		}
		else if (version_=="2017") {
			costS_tau = compute_cosThetaS(taus[0].p4());
			dr_taus = taus[0].p4().DeltaR(taus[1].p4());
			dr_lep_tau_lead = leptons[0].p4().DeltaR(taus[0].p4());
			dr_lep_tau_sublead = leptons[0].p4().DeltaR(taus[1].p4());
			//assert(taus[0].charge()!=taus[1].charge());
			dr_lep_tau_ss = (leptons[0].charge()==taus[0].charge()) ? dr_lep_tau_lead : dr_lep_tau_sublead;
			lep0_conept = leptons[0].conept();
			mindr_lep0_jet = compute_min_dr(leptons[0].p4(),jets);
			mindr_tau0_jet = compute_min_dr(taus[0].p4(),jets);
			mindr_tau1_jet = compute_min_dr(taus[1].p4(),jets);
			mTauTauVis = (taus[0].p4()+taus[1].p4()).M();
			mT_met_lep0 = compute_mT_lep(leptons[0], MET, METphi);
			met = MET;
			HT = MHT; // FIXME
			tau0_pt = taus[0].pt();
			tau1_pt = taus[1].pt();
			nbtags_loose = nbtagsloose;
			nbtags_medium = nbtagsmedium;

			compute_tauDecay_variables(taus);

		}
		else if (version_=="test") {
			compute_tauDecay_variables(taus, true);
		}
		break;
	case Analyze_2lss1tau:
		assert(leptons.size()>1);
		assert(taus.size()>0);
		
		if (version_=="2016") {
			mindr_lep0_jet = compute_min_dr(leptons[0].p4(),jets);
			mindr_lep1_jet = compute_min_dr(leptons[1].p4(),jets);
			max_lep_eta = compute_max_lep_eta(leptons);
			met = MET;
			mht = MHT;
			mT_met_lep0 = compute_mT_lep(leptons[0], MET, METphi);
			lep0_conept = leptons[0].conept();;
			lep1_conept = leptons[1].conept();;
			dr_leps = leptons[0].p4().DeltaR(leptons[1].p4());
			tau0_pt = taus[0].pt();
			dr_lep0_tau = leptons[0].p4().DeltaR(taus[0].p4());
			dr_lep1_tau = leptons[1].p4().DeltaR(taus[0].p4());
			mvis_lep0_tau = (leptons[0].p4()+taus[0].p4()).M();
			mvis_lep1_tau = (leptons[1].p4()+taus[0].p4()).M();
			tau0_decaymode = taus[0].decaymode();
			tau0_E = taus[0].p4().E();
			tau0_easym = compute_upsilon(taus[0]);
		}
		else if (version_=="2017") {
			std::cout << "TO BE ADDED!!!!!!" << std::endl;
		}
		break;
	case Analyze_3l1tau:
		assert(leptons.size()>2);
		assert(taus.size()>0);
		
		if (version_=="2016") {
			max_lep_eta = compute_max_lep_eta(leptons);
			mindr_lep0_jet = compute_min_dr(leptons[0].p4(),jets);
			mindr_lep1_jet = compute_min_dr(leptons[1].p4(),jets);
			mindr_lep2_jet = compute_min_dr(leptons[2].p4(),jets);
			mT_met_lep0 = compute_mT_lep(leptons[0], MET, METphi);
			lep0_conept = leptons[0].conept();
			lep1_conept = leptons[1].conept();
			lep2_conept = leptons[2].conept();
		}
		else if (version_=="2017") {
			std::cout << "TO BE ADDED!!!!!!" << std::endl;
		}
		break;
	default:
		std::cout << "WARNING analysis type is not supported. Nothing is computed."
				  << std::endl;
		break;
	}
}

void mvaNtuple::compute_tauDecay_variables(const std::vector<miniTau>& taus,bool test)
{
	assert(taus.size()>1);
	//assert(taus[0].charge()!=taus[1].charge());
	int ip, im;
	if (taus[0].charge()>0) {ip=0; im=1;}
	else {ip=1; im=0;}
	
	taup_decaymode = taus[ip].decaymode();
	taum_decaymode = taus[im].decaymode();
	taup_easym = compute_upsilon(taus[ip]);
	taum_easym = compute_upsilon(taus[im]);
	taup_cosPsi = compute_cosPsi(taus[ip]);
	taum_cosPsi= compute_cosPsi(taus[im]);
	evisTaus_diff = (taus[im].p4().E()-taus[ip].p4().E());
	evisTaus_sum = (taus[im].p4()+taus[ip].p4()).E();
	evisTausAsym = evisTaus_diff / evisTaus_sum;

	if (test) {
		// LorentzVectors of tau decay products
		TLorentzVector pp1, pp2, pp3;
		TLorentzVector pm1, pm2, pm3;

		std::vector<TLorentzVector> taup_chargedHadrons =
			taus[ip].get_signalChargedHadrCands();
		std::vector<TLorentzVector> taum_chargedHadrons =
			taus[im].get_signalChargedHadrCands();

		std::cout << taup_chargedHadrons[0].Pt() << " "
				  << taup_chargedHadrons[0].Eta() << " "
				  << taup_chargedHadrons[0].Phi() << " "
				  << taup_chargedHadrons[0].M() << " "
				  << std::endl;
		
		assert(taup_chargedHadrons.size()>0);
		assert(taum_chargedHadrons.size()>0);
		/*
		if (taup_decaymode == 0) { // 1prong
			pp1 = taup_chargedHadrons[0];
		}
		else if (taup_decaymode == 1) { // 1prong1pi0
			pp1 = taup_chargedHadrons[0];
			pp2 = taus[ip].p4()-taup_chargedHadrons[0];
		}
		else if (taup_decaymode == 10) { // 3prongs
			assert(taup_chargedHadrons.size()>2);
			pp1 = taup_chargedHadrons[0];
			pp2 = taup_chargedHadrons[1];
			pp3 = taup_chargedHadrons[2];
		}

		if (taum_decaymode == 0) { // 1prong
			pm1 = taum_chargedHadrons[0];
		}
		else if (taum_decaymode == 1) { // 1prong1pi0
			pm1 = taum_chargedHadrons[0];
			pm2 = taus[im].p4()-taum_chargedHadrons[0];
		}
		else if (taum_decaymode == 10) { // 3prongs
			assert(taum_chargedHadrons.size()>2);
			pm1 = taum_chargedHadrons[0];
			pm2 = taum_chargedHadrons[1];
			pm3 = taum_chargedHadrons[2];
		}
		*/

		pp1_pp1 = pp1*pp1; //std::cout << pp1_pp1 << std::endl;
		pp1_pp2 = pp1*pp2; //std::cout << pp1_pp2 << std::endl;
		pp1_pp3 = pp1*pp3; //std::cout << pp1_pp3 << std::endl;
		pp1_pm1 = pp1*pm1; //std::cout << pp1_pm1 << std::endl;
		pp1_pm2 = pp1*pm2; //std::cout << pp1_pm2 << std::endl;
		pp1_pm3 = pp1*pm3; //std::cout << pp1_pm3 << std::endl;
		pp2_pp2 = pp2*pp2; //std::cout << pp2_pp2 << std::endl;
		pp2_pp3 = pp2*pp3; //std::cout << pp2_pp3 << std::endl;
		pp2_pm1 = pp2*pm1; //std::cout << pp2_pm1 << std::endl;
		pp2_pm2 = pp2*pm2; //std::cout << pp2_pm2 << std::endl;
		pp2_pm3 = pp2*pm3; //std::cout << pp2_pm3 << std::endl;
		pp3_pp3 = pp3*pp3; //std::cout << pp3_pp3 << std::endl;
		pp3_pm1 = pp3*pm1; //std::cout << pp3_pm1 << std::endl;
		pp3_pm2 = pp3*pm2; //std::cout << pp3_pm2 << std::endl;
		pp3_pm3 = pp3*pm3; //std::cout << pp3_pm3 << std::endl;
		pm1_pm1 = pm1*pm1; //std::cout << pm1_pm1 << std::endl;
		pm1_pm2 = pm1*pm2; //std::cout << pm1_pm2 << std::endl;
		pm1_pm3 = pm1*pm3; //std::cout << pm1_pm3 << std::endl;
		pm2_pm2 = pm2*pm2; //std::cout << pm2_pm2 << std::endl;
		pm2_pm3 = pm2*pm3; //std::cout << pm2_pm3 << std::endl;
		pm3_pm3 = pm3*pm3; //std::cout << pm3_pm3 << std::endl;
	}
}

float mvaNtuple::compute_average_dr(const std::vector<TLorentzVector>& lvs)
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

float mvaNtuple::compute_max_dr(const std::vector<TLorentzVector>& lvs)
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

float mvaNtuple::compute_min_dr(const TLorentzVector& l,
								const std::vector<TLorentzVector>& vjs)
{
	if (vjs.size() < 1) return -9999.;

	float mindr = 666.;

	for (auto & j : vjs) {
		float dr = l.DeltaR(j);
		if (dr < mindr) mindr = dr;
	}

	return mindr;
}

float mvaNtuple::compute_cosThetaS(const TLorentzVector& taup4_lead
								   /*,const TLorentzVector& taup4_sublead*/)
{
	// ????
	// https://github.com/HEP-KBFI/tth-htt/blob/master/bin/analyze_2l_2tau.cc#L104-L115
	
	//TLorentzVector ditaup4 = taup4_lead + taup4_sublead;
	//TLorentzVector tauboost = taup4_lead;

	return std::abs(taup4_lead.CosTheta());	
}

float mvaNtuple::compute_mT_lep(const miniLepton& lepton, float MET, float METphi)
{
	return sqrt(2*lepton.conept()*MET*(1-cos(lepton.phi()-METphi)));
}

float mvaNtuple::compute_max_lep_eta(const std::vector<miniLepton>& leptons)
{
	float max_eta = 0.;

	for (const auto & lep : leptons) {
		if (fabs(lep.eta())>max_eta) max_eta = fabs(lep.eta());
	}

	return max_eta;
}

float mvaNtuple::compute_upsilon(const miniTau& tau)
{
	// (2*E_ldgtrk - E_vis) / E_vis
	auto p4_ldgtrk = tau.leadtrackP4();
	auto p4_visible = tau.p4();
	return (2.*p4_ldgtrk.Energy()/p4_visible.Energy() - 1);
}

float mvaNtuple::compute_cosPsi(const miniTau& tau, float mass)
{
	std::vector<TLorentzVector> chargedhadrons = tau.get_signalChargedHadrCands();
	if (chargedhadrons.size()<3)
		return 1.;
	else
		// FIXME: directly compute from TLorentzVector?
		return compute_cosPsi(chargedhadrons[0], chargedhadrons[1],
							  chargedhadrons[2], mass);	
}

float mvaNtuple::compute_cosPsi(const TLorentzVector& p1, const TLorentzVector& p2,
							  const TLorentzVector& p3, float mass)
// lab frame, 3 prongs mode only
// angle between the direction of the normal to the plane defined by the three pions and the direction of flight of a1
{
	float s = (p1+p2+p3).M2();
	float s12 = (p1+p2).M2();
	float s23 = (p2+p3).M2();
	float s13 = (p1+p3).M2();

	TVector3 v1 = p1.Vect();
	TVector3 v2 = p2.Vect();
	TVector3 v3 = p3.Vect();

	float m2 = mass*mass;

	float numerator = v1.Dot(v2.Cross(v3)) / (v1+v2+v3).Mag();
	float denominator = TMath::Sqrt(-lam(lam(s,s12,m2),lam(s,s13,m2),
										  lam(s,s23,m2)));

	return 8*s*numerator/denominator;
}

float mvaNtuple::lam(float a, float b, float c)
{
	return a*a + b*b + c*c - 2*a*b - 2*b*c - 2*a*c;
}

#endif
