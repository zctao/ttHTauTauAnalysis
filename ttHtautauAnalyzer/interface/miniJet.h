#ifndef miniJet_h
#define miniJet_h

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
#include "DataFormats/PatCandidates/interface/Jet.h"
#endif

#include "TLorentzVector.h"
#include <string>
#include <iostream>

class miniJet {
 public:
	// constructors and destructors	
	miniJet(const TLorentzVector&, float, float, float);
	
	~miniJet(){};
	
	// member functions
	float pt() const {return pt_;}
	float eta() const {return eta_;}
	float phi() const {return phi_;}
	float energy() const {return energy_;}
	float csv() const {return csv_;}
	float flavor() const {return flavor_;}
	float qgLikelihood() const {return qgLikelihood_;}

	TLorentzVector p4() const;
	void dump() const;
	
 protected:
	float pt_;
	float eta_;
	float phi_;
	float energy_;
	float csv_;
	float flavor_;
	float qgLikelihood_;
};

#endif
