#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/miniLepton.h"

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
miniLepton::miniLepton(const pat::Electron& ele)
{
	pt_ = ele.pt();
	eta_ = ele.eta();
	phi_ = ele.phi();
	mass_ = ele.mass();
	charge_ = ele.charge();
	pdgid_ = ele.pdgId();

	if (ele.hasUserFloat("conePt"))
		conept_ = ele.userFloat("conePt");
	
	if (ele.hasUserInt("isTightCharge"))
		tightcharge_ = ele.userInt("isTightCharge");

	if (ele.hasUserInt("isLoose"))
		isloose_ = ele.hasUserInt("isLoose");

	if (ele.hasUserInt("isFakeable"))
		isfakeable_ = ele.userInt("isFakeable");

	if (ele.hasUserInt("isTight"))
		istight_ = ele.userInt("isTight");

	if (ele.hasUserInt("MCMatchType"))
		mcmatchtype_ = ele.userInt("MCMatchType");
}

miniLepton::miniLepton(const pat::Muon& mu)
{
	pt_ = mu.pt();
	eta_ = mu.eta();
	phi_ = mu.phi();
	mass_ = mu.mass();
	charge_ = mu.charge();
	pdgid_ = mu.pdgId();

	if (mu.hasUserFloat("conePt"))
		conept_ = mu.userFloat("conePt");
	
	if (mu.hasUserInt("isTightCharge"))
		tightcharge_ = mu.userInt("isTightCharge");

	if (mu.hasUserInt("isLoose"))
		isloose_ = mu.hasUserInt("isLoose");

	if (mu.hasUserInt("isFakeable"))
		isfakeable_ = mu.userInt("isFakeable");

	if (mu.hasUserInt("isTight"))
		istight_ = mu.userInt("isTight");

	if (mu.hasUserInt("MCMatchType"))
		mcmatchtype_ = mu.userInt("MCMatchType");
}
	
#endif

TLorentzVector miniLepton::p4() const
{
	TLorentzVector l;
	l.SetPtEtaPhiM(pt_, eta_, phi_, mass_);
	return l;
}
