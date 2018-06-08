#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/miniJet.h"

miniJet::miniJet(const TLorentzVector& j, float csv, float flavor, float qgL)
{
	pt_ = j.Pt();
	eta_ = j.Eta();
	phi_ = j.Phi();
	energy_ = j.Energy();
	csv_ = csv;
	flavor_ = flavor;
	qgLikelihood_ = qgL;
}

TLorentzVector miniJet::p4() const
{
	TLorentzVector j;
	j.SetPtEtaPhiE(pt_, eta_, phi_, energy_);
	return j;
}

float miniJet::mass() const
{
	return p4().M();
}

LorentzVector miniJet::lorentzvector() const
{
	LorentzVector l;
	l.SetPt(pt_);
	l.SetEta(eta_);
	l.SetPhi(phi_);
	l.SetM(mass());
	return l;
}

void miniJet::dump() const
{
	std::cout << "pt: " << pt() << " eta: " << eta() << " phi: " << phi()
			  << " energy: " << energy() << " csv: " << csv()
			  << " flavor: " << flavor() << " QGLikelihood: " << qgLikelihood()
			  << std::endl;;
}
