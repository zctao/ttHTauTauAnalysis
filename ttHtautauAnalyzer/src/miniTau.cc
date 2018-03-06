#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/miniTau.h"

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
miniTau::miniTau(const pat::Tau& tau, bool addDaughters)
{
	pt_ = tau.pt();
	eta_ = tau.eta();
	phi_ = tau.phi();
	mass_ = tau.mass();
	charge_ = tau.charge();
	//pdgid_ = tau.pdgId();

	if (tau.hasUserInt("isLoose"))
		isloose_ = tau.hasUserInt("isLoose");

	if (tau.hasUserInt("isTight"))
		istight_ = tau.userInt("isTight");

	if (tau.hasUserInt("MCMatchType"))
		mcmatchtype_ = tau.userInt("MCMatchType");

	if (addDaughters) {
		// TODO
	}
}
#endif

miniTau::miniTau(const TLorentzVector& t, int charge, int decaymode,
				 bool isloose, bool istight)//, int mcmatchtype)
{
	pt_ = t.Pt();
	eta_ = t.Eta();
	phi_ = t.Phi();
	mass_ = t.M();
	charge_ = charge;
	decaymode_ = decaymode;
	isloose_ = isloose;
	istight_ = istight;
	//mcmatchtype_ = mcmatchtype;
}

miniTau::miniTau(const TLorentzVector& t, int charge, int decaymode,
				 bool isloose, bool istight, //int mcmatchtype,
				 const std::vector<TLorentzVector>& signalChargedHadrCands,
				 const std::vector<TLorentzVector>& signalGammaCands,
				 const std::vector<TLorentzVector>& signalNeutrHadrCands)
{
	pt_ = t.Pt();
	eta_ = t.Eta();
	phi_ = t.Phi();
	mass_ = t.M();
	charge_ = charge;
	decaymode_ = decaymode;
	isloose_ = isloose;
	istight_ = istight;
	//mcmatchtype_ = mcmatchtype;
	signalChargedHadrCands_ = signalChargedHadrCands;
	signalGammaCands_ = signalGammaCands;
	signalNeutrHadrCands_ = signalNeutrHadrCands;
}

TLorentzVector miniTau::p4() const
{
	TLorentzVector t;
	t.SetPtEtaPhiM(pt_,eta_,phi_,mass_);
	return t;
}

TLorentzVector miniTau::chargedDaughtersP4() const
{
	TLorentzVector ChargedP4;

	for (const TLorentzVector & ch : signalChargedHadrCands_)
		ChargedP4 += ch;

	return ChargedP4;
}

TLorentzVector miniTau::neutralDaughtersP4() const
{
	TLorentzVector NeutralP4;

	for (const TLorentzVector & g : signalGammaCands_)
		NeutralP4 += g;
	for (const TLorentzVector & nh : signalNeutrHadrCands_)
		NeutralP4 += nh;

	return NeutralP4;
}

TLorentzVector miniTau::leadtrackP4() const
{
	TLorentzVector ldgtrkP4(0.,0.,0.,0.);
	for (const TLorentzVector& ch : signalChargedHadrCands_) {
		if (ch.Pt()>ldgtrkP4.Pt())
			ldgtrkP4 = ch;
	}

	return ldgtrkP4;
}

bool miniTau::isGenMatched() const
{
	return (mcmatchtype_==1 or mcmatchtype_==2 or mcmatchtype_==3 or
			mcmatchtype_==4 or mcmatchtype_==5);
}
