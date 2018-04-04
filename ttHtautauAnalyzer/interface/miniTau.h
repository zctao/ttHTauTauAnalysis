#ifndef miniTau_h
#define miniTau_h

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
#include "DataFormats/PatCandidates/interface/Tau.h"
#endif

#include "TLorentzVector.h"

class miniTau
{
 public:
	// constructor and destructor
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	miniTau(const pat::Tau&, bool addDaughters = false);
#endif

	miniTau(const TLorentzVector&, int, int, bool, bool, int mcmatchtype=-9999);
	miniTau(const TLorentzVector&, int, int, bool, bool, 
			const std::vector<TLorentzVector>&,const std::vector<TLorentzVector>&,
			const std::vector<TLorentzVector>&, int mcmatchtype=-9999);

	~miniTau(){};

	// member functions

	void set_pt(float ipt) {pt_ = ipt;}
	void set_eta(float ieta) {eta_ = ieta;}
	void set_phi(float iphi) {phi_ = iphi;}
	void set_mass(float imass) {mass_ = imass;}
	void set_charge(int icharge) {charge_ = icharge;}
	void set_decaymode(int idecaymode) {decaymode_ = idecaymode;}
	void set_MCMatchType(int imctype) {mcmatchtype_ = imctype;}
	//void set_pdgId(int id) {pdgid_ = id;}
	void set_signalChargedHadrCands(std::vector<TLorentzVector>& ch) {signalChargedHadrCands_ = ch;}
	void set_signalGammaCands(std::vector<TLorentzVector>& ga) {signalGammaCands_ = ga;}
	void set_signalNeutrHadrCands(std::vector<TLorentzVector>& nh) {signalNeutrHadrCands_ = nh;}
	void set_tauIDWPindex(int iWP) {tauIDMVAWP_ = iWP;}
	void set_tauIDWPindex(bool, bool, bool, bool);

	// TODO: check value was set before returning
	float pt() const {return pt_;}
	float eta() const {return eta_;}
	float phi() const {return phi_;}
	float mass() const {return mass_;}
	int charge() const {return charge_;}
	int decaymode() const {return decaymode_;}
	//int pdgId() const {return pdgid_;}
	bool passLooseSel() const {return isloose_;}
	bool passTightSel() const {return istight_;} // depend on analysis type
	int MCMatchType() const {return mcmatchtype_;}
	int tauIDMVAWPindex() const {assert(mvawp_set_); return tauIDMVAWP_;}
	bool passMVAID(char) const;
	bool isGenMatched() const;
	TLorentzVector p4() const;
	TLorentzVector chargedDaughtersP4() const;
	TLorentzVector neutralDaughtersP4() const;
	TLorentzVector leadtrackP4() const;

	void dump() const;

	std::vector<TLorentzVector> get_signalChargedHadrCands() const {return signalChargedHadrCands_;}
	std::vector<TLorentzVector> get_signalGammaCands() const {return signalGammaCands_;}
	std::vector<TLorentzVector> get_signalNeutrHadrCands() const {return signalNeutrHadrCands_;}

 protected:

	float pt_;
	float eta_;
	float phi_;
	float mass_;
	int charge_;
	int decaymode_;
	//float pdgid_;
	bool isloose_;
	bool istight_;
	int mcmatchtype_;
	// tau ID MVA work point: 0='Loose', 1='Medium', 2='Tight', 3='VTight'
	int tauIDMVAWP_;
	bool mvawp_set_;

	std::vector<TLorentzVector> signalChargedHadrCands_;
	std::vector<TLorentzVector> signalGammaCands_;
	std::vector<TLorentzVector> signalNeutrHadrCands_;
};

#endif
