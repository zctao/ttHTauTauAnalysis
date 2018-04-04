#ifndef miniLepton_h
#define miniLepton_h

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#endif

#include "TLorentzVector.h"

class miniLepton
{
 public:
	// constructor and destructor
#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
	miniLepton(const pat::Electron&);
	miniLepton(const pat::Muon&);
#endif
	miniLepton(const TLorentzVector&, float, int pdgid=-9999, int charge=-9999, 
			   bool isloose=false, bool isfakeable=false, bool istight=false,
			   bool tightcharge=false, int mcmatchtype=-9999);
	//miniLepton(float, float);
	
	~miniLepton(){};

	// member functions

	void set_pt(float ipt) {pt_ = ipt;}
	void set_conept(float iconept) {conept_ = iconept;}
	void set_eta(float ieta) {eta_ = ieta;}
	void set_phi(float iphi) {phi_ = iphi;}
	void set_mass(float imass) {mass_ = imass;}
	void set_charge(int icharge) {charge_ = icharge;}
	void set_pdgId(int id) {pdgid_ = id;}
	void set_MCMatchType(int imctype) {mcmatchtype_ = imctype;}

	// TODO: check value was set before returning
	float pt() const {return pt_;}
	float conept() const {return conept_;}
	float eta() const {return eta_;}
	float phi() const {return phi_;}
	float mass() const {return mass_;}
	int charge() const {return charge_;}
	int pdgId() const {return pdgid_;}
	bool passTightCharge() const {return tightcharge_;}
	bool passLooseSel() const {return isloose_;}
	bool passFakeableSel() const {return isfakeable_;}
	bool passTightSel() const {return istight_;}
	int MCMatchType() const {return mcmatchtype_;}
	bool isGenMatched() const;
	TLorentzVector p4() const;

	void dump() const;
		
 private:

	float pt_;
	float conept_;
	float eta_;
	float phi_;
	float mass_;
	int   charge_;
	int   pdgid_;
	bool  tightcharge_;
	bool  isloose_;
	bool  isfakeable_;
	bool  istight_;
	int   mcmatchtype_;
};

#endif
