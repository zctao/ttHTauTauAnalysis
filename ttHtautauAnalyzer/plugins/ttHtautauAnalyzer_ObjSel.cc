#ifndef ttHtautauAnalyzer_ObjSel_cc
#define ttHtautauAnalyzer_ObjSel_cc

#include "ttHTauTauAnalysis/ttHtautauAnalyzer/plugins/ttHtautauAnalyzer.h"

///////////////////
// leptons
template <typename T> float ttHtautauAnalyzer::ConePt(const T& lep)
{
	//return (lep.userFloat("leptonMVA") > 0.75) ? lep.pt() :
	//	(0.85 * lep.pt() / lep.userFloat("nearestJetPtRatio"));

	return lep.userFloat("correctedPt");
}
template float ttHtautauAnalyzer::ConePt<pat::Electron>(const pat::Electron&);
template float ttHtautauAnalyzer::ConePt<pat::Muon>(const pat::Muon&);

void ttHtautauAnalyzer::SortByConept(std::vector<pat::Electron>& electrons)
{
	std::sort(electrons.begin(), electrons.end(), [] (pat::Electron e1, pat::Electron e2) {return e1.userFloat("conePt") > e2.userFloat("conePt");});
}

void ttHtautauAnalyzer::SortByConept(std::vector<pat::Muon>& muons)
{
	std::sort(muons.begin(), muons.end(), [] (pat::Muon m1, pat::Muon m2) {return m1.userFloat("conePt") > m2.userFloat("conePt");});
}

void ttHtautauAnalyzer::SortByConept(std::vector<miniLepton>& leptons)
{
	std::sort(leptons.begin(), leptons.end(), [] (miniLepton l1, miniLepton l2) {return l1.conept() > l2.conept();});
}

// get selected collection
template <typename T> std::vector<T> ttHtautauAnalyzer::getSelectedLeptons(const std::vector<T>& leptons, const std::string& selection)
{
	std::vector<T> selected;
	for (auto & lep : leptons) {
	    if (lep.hasUserInt(selection)) {
			if (not lep.userInt(selection)) continue;
		}
		else {
			if (selection == "isLoose") {
				if (not isLooseID(lep)) continue;
			}
			else if (selection == "isFakeable") {
				if (not isFakeableID(lep)) continue;
			}
			else if (selection == "isTight") {
				if (not isTightID(lep)) continue;
			}
			else {
				std::cerr << "ERROR: not valid lepton selection. " << std::endl;
				std::cerr << "Avaiable selections are: "
						  << "isLoose, isFakebale and isTight" << std::endl;
				assert(0);
			}
		}
		selected.push_back(lep);
	}

	return selected;
}
template std::vector<pat::Muon> ttHtautauAnalyzer::getSelectedLeptons(const std::vector<pat::Muon>&, const std::string&);
template std::vector<pat::Electron> ttHtautauAnalyzer::getSelectedLeptons(const std::vector<pat::Electron>&, const std::string&);

// add ID related flags/variables
template <typename T> void ttHtautauAnalyzer::addIDFlags(std::vector<T>& leptons)
{
	for (auto & lep : leptons) {
		lep.addUserInt("isLoose", isLooseID(lep));
		lep.addUserInt("isFakeable", isFakeableID(lep));
		lep.addUserInt("isTight", isTightID(lep));
		lep.addUserInt("isTightCharge", isTightCharge(lep));
		lep.addUserFloat("conePt", ConePt(lep));
	}
}
template void ttHtautauAnalyzer::addIDFlags(std::vector<pat::Electron>&);
template void ttHtautauAnalyzer::addIDFlags(std::vector<pat::Muon>&);

// get matched gen particle and add mc match information
template <typename T>
std::vector<const reco::GenParticle*> ttHtautauAnalyzer::getGenLepMatchInfo
(std::vector<T>& leptons, edm::Handle<reco::GenParticleCollection> MC_particles)
{
  std::vector<const reco::GenParticle*> matched_collection;

  if (isdata_)
    return matched_collection;

  for (auto & lep : leptons) {
    // get matched gen particle
    auto matchedLep = getMatchedGenParticle(lep, *MC_particles);
    lep.addUserInt("MCMatchType", getMCMatchType(lep, matchedLep));
    lep.addUserInt("isGenPhotonMatched",
                   isGenPhotonMatched(lep, *MC_particles, false));
    lep.addUserInt("isPromptGenPhotonMatched",
                   isGenPhotonMatched(lep, *MC_particles, true));

    matched_collection.push_back(matchedLep);
  }

  return matched_collection;
}

template std::vector<const reco::GenParticle*> ttHtautauAnalyzer::getGenLepMatchInfo(std::vector<pat::Electron>&, edm::Handle<reco::GenParticleCollection>);
template std::vector<const reco::GenParticle*> ttHtautauAnalyzer::getGenLepMatchInfo(std::vector<pat::Muon>&, edm::Handle<reco::GenParticleCollection>);

///////////////////
// muons
bool ttHtautauAnalyzer::isLooseID(const pat::Muon& mu) const
{
	bool passKinematic = mu.pt() > 5. and std::abs(mu.eta()) < 2.4;
	
	return passKinematic and mu.userFloat("idPreselection") > 0.5;
}

bool ttHtautauAnalyzer::isFakeableID(const pat::Muon& mu) const
{
	return mu.userFloat("idFakeable") > 0.5;
}

bool ttHtautauAnalyzer::isTightID(const pat::Muon& mu) const 
{
	return mu.userFloat("idMVABased") > 0.5;
}

bool ttHtautauAnalyzer::isTightCharge(const pat::Muon& mu) const 
{
	bool tightcharge = false;

	if (mu.innerTrack().isAvailable())
		tightcharge = mu.innerTrack()->ptError()/mu.innerTrack()->pt() < 0.2;

	return tightcharge;
}

///////////////////
// electrons
bool ttHtautauAnalyzer::isLooseID(const pat::Electron& ele) const 
{
	bool passKinematic = ele.pt() > 7. and fabs(ele.eta()) < 2.5;
	//bool passPhotonVeto = (ele.passConversionVeto() and
	//					   ele.userFloat("numMissingHits") <= 1);

	return (passKinematic and ele.userFloat("idPreselection") > 0.5);
}

bool ttHtautauAnalyzer::isFakeableID(const pat::Electron& ele) const
{
	return (ele.userFloat("idFakeable") > 0.5 //and
			//ele.userFloat("numMissingHits") == 0 and
			//ele.passConversionVeto()
			);
}

bool ttHtautauAnalyzer::isTightID(const pat::Electron& ele) const
{
	return (ele.userFloat("idMVABased") > 0.5// and
			//ele.userFloat("numMissingHits") == 0 and
			//ele.passConversionVeto()
			);
}

bool ttHtautauAnalyzer::isTightCharge(const pat::Electron& ele) const
{
	bool tightcharge =
		ele.isGsfCtfScPixChargeConsistent()+ele.isGsfScPixChargeConsistent() > 1;
	
	//std::cout << "elel tight charge : " << ele.isGsfCtfScPixChargeConsistent()
	//		  << " " << ele.isGsfScPixChargeConsistent() << std::endl;
	
	return tightcharge;
}

///////////////////
// taus

std::vector<pat::Tau> ttHtautauAnalyzer::getPreselTaus(const std::vector<pat::Tau>& taus)
{
	std::vector<pat::Tau> preseltaus;
	for (auto & tau : taus) {

		if (tau.hasUserInt("isLoose")) {
			if (not tau.userInt("isLoose")) continue;
		}
		else {
			if (not isLooseID(tau)) continue;
		}
	
		preseltaus.push_back(tau);	
	}
	
	return preseltaus;
}

std::vector<pat::Tau> ttHtautauAnalyzer::getSelectedTaus(const std::vector<pat::Tau>& taus)
{
	std::vector<pat::Tau> seltaus;

	for (auto & tau : taus) {

		if (tau.hasUserInt("isTight")) {
			if (not tau.userInt("isTight")) continue;
		}
		else {
			if (not isTightID(tau)) continue;
		}

		seltaus.push_back(tau);
	}
	
	return seltaus;
}

// better to move this to leptonID package
bool ttHtautauAnalyzer::isLooseID(const pat::Tau& tau) const
{
    // 18 GeV pT cut is somewhat arbitrary. (Should be 20 GeV)
    // It accounts for the upward tau energy scale variation (should be <10%)
    // for systematic study
	bool passKinematic = tau.pt() > 18. and std::abs(tau.eta()) < 2.3;

	// tauID("byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017")
	bool passID = tau.userFloat("idPreselection") > 0.5;

	return (passKinematic and passID);
}

bool ttHtautauAnalyzer::isTightID(const pat::Tau& tau) const
{

	// loose ID is the prerequisite for tight ID
	if (tau.hasUserInt("isLoose")) {
		if (not tau.userInt("isLoose")) return false;
	}
	else {
		if (not isLooseID(tau)) return false;
	}

    // tauID("byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017")
    return tau.userFloat("idSelection") > 0.5;
}

void ttHtautauAnalyzer::addIDFlagsTau(std::vector<pat::Tau>& taus)
{
	for (auto & tau : taus) {
		tau.addUserInt("isLoose", isLooseID(tau));
		tau.addUserInt("isTight", isTightID(tau));
	}
}

std::vector<const reco::GenParticle*> ttHtautauAnalyzer::getGenTauMatchInfo
(std::vector<pat::Tau>& taus, edm::Handle<reco::GenParticleCollection> MC_particles)
{
  std::vector<const reco::GenParticle*> matched_collection;

  if (isdata_)
    return matched_collection;

  for (auto & tau : taus) {
    // get matched gen particle
    auto matchedTau = getMatchedGenParticle(tau, *MC_particles);
    tau.addUserInt("MCMatchType", getMCMatchType(tau, matchedTau));

    matched_collection.push_back(matchedTau);
  }

  return matched_collection;
}

///////////////////
// jets
bool ttHtautauAnalyzer::isTightJet(const pat::Jet& jet) const
{
	// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017#Preliminary_Recommendations_for
	
	bool istight = false;
	
	if (fabs(jet.eta()) <= 2.7) {
		istight =
			jet.neutralHadronEnergyFraction() < 0.90 and
			jet.neutralEmEnergyFraction() < 0.90 and
		    jet.numberOfDaughters() > 1;
		if (fabs(jet.eta()) <= 2.4) {
			istight = istight and
				jet.chargedHadronEnergyFraction() > 0 and
				jet.chargedMultiplicity() > 0;
		}																		  
	}
	else if (fabs(jet.eta()) <= 3.0) {
		istight =
			jet.neutralEmEnergyFraction() < 0.99 and
		    jet.neutralEmEnergyFraction() > 0.02 and
			jet.neutralMultiplicity() > 2;
	}
	else {
		istight =
			jet.neutralEmEnergyFraction() < 0.90 and
			jet.neutralHadronEnergyFraction() > 0.02 and
			jet.neutralMultiplicity() > 10;
	}

	return istight;
}

std::vector<pat::Jet> ttHtautauAnalyzer::getSelectedJets(const std::vector<pat::Jet>& jets)
{
	std::vector<pat::Jet> seljets;
	for (const auto & jet : jets) {
		float jetPt = jet.pt();
		// To account for upwards JEC variation for systematic study
		if (jet.hasUserFloat("jesUnc")) jetPt *= 1.+jet.userFloat("jesUnc");
		
		bool passKinematic = (jetPt > 25. and fabs(jet.eta()) < 2.4);

		bool passID = isTightJet(jet);

		if (passKinematic and passID)
			seljets.push_back(jet);
	}
	return seljets;
}

void ttHtautauAnalyzer::addJetCorrections(
    std::vector<pat::Jet>& input_jets, JetCorrectionUncertainty* jecUnc,
    JME::JetResolution* jetRes, JME::JetResolutionScaleFactor* jetRes_sf,
    double rho, const std::vector<reco::GenJet>& genJets,
    std::mt19937& random_generator)
{
  for (auto & jet : input_jets) {
    // JES
    if (jecUnc) {
      jecUnc->setJetEta(jet.eta());
      jecUnc->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
      float unc = jecUnc->getUncertainty(true);
      jet.addUserFloat("jesUnc",unc);
    }
    
    // JER
    if (jetRes and jetRes_sf) {
      double jet_res = jetRes->getResolution({{JME::Binning::JetPt, jet.pt()},
            {JME::Binning::JetEta, jet.eta()}, {JME::Binning::Rho, rho}});
      double jer_sf = jetRes_sf->getScaleFactor({{JME::Binning::JetEta, jet.eta()}});
      double jer_sf_up = jetRes_sf->getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::UP);
      double jer_sf_down = jetRes_sf->getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::DOWN);
      
      // Matched gen jet
      const reco::GenJet* genjet = getMatchedGenJet(jet, genJets, jet_res);
      
      // Smearing factors
      float smearFactor = getJERSmearFactor(jer_sf, jet_res, &jet, genjet,
                                            random_generator);
      float smearFactor_up = getJERSmearFactor(jer_sf_up, jet_res, &jet, genjet,
                                               random_generator);
      float smearFactor_down = getJERSmearFactor(jer_sf_down, jet_res, &jet,
                                               genjet, random_generator);
      jet.addUserFloat("jerSmearFactor", smearFactor);
      jet.addUserFloat("jerSmearFactor_up", smearFactor_up);
      jet.addUserFloat("jerSmearFactor_down", smearFactor_down);
    }
  }
}

float ttHtautauAnalyzer::getJERSmearFactor(double jerSF, double jetRes,
                                           const pat::Jet* jet,
                                           const reco::GenJet* genjet,
                                           std::mt19937& random_generator,
                                           bool debug)
// https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L247
{
  float factor = 1.;
  
  if (genjet) {
    //std::cout << "Found gen jet! jerSF = " << jerSF << std::endl;
    // Scaling method if there is a good gen jet matched to the reco jet
    float dPtOverPt = std::abs(jet->pt() - genjet->pt()) / jet->pt(); 
    factor = 1. + (jerSF - 1.) * dPtOverPt;
  }
  else {
    // Stochastic smearing if no matched gen jet
    if (jerSF >= 1.) {
      //std::cout << "Stochastic smearing. jerSF = " << jerSF << std::endl;
      double sigma = jetRes * std::sqrt(jerSF * jerSF - 1.);
      std::normal_distribution<> d(0, sigma);
      factor = 1. + d(random_generator);
    }
    else if (debug) {
      std::cout << "jerSF < 1 and no matched gen jet found. ";
      std::cout << "Cannot smear the jet. ";
      std::cout << "jet eta = " << jet->eta() << " jerSF = " << jerSF << std::endl;
    }
  }

  return max(factor, 0.f);
}

float ttHtautauAnalyzer::getJetCSV(const pat::Jet& jet)
{
	float csv = jet.bDiscriminator("pfDeepCSVJetTags:probb")
		+ jet.bDiscriminator("pfDeepCSVJetTags:probbb");

	if (isnan(csv))
		csv = -2.;
	
	return csv;

	/*
	float defaultvalue = -0.1;
	//float csv = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	//float csv = jet.bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll");
	if (std::isnan(csv)) return defaultvalue;
	if (csv > 1.) return 1.;
	if (csv < 0.) return defaultvalue;

	return csv;
	*/
}

float ttHtautauAnalyzer::getJetDeepCvsB(const pat::Jet& jet)
{
	float probC = jet.bDiscriminator("pfDeepCSVJetTags:probc");
	float probB = jet.bDiscriminator("pfDeepCSVJetTags:probb") +
		jet.bDiscriminator("pfDeepCSVJetTags:probbb");
	
	return probC/(probC+probB);
}

float ttHtautauAnalyzer::getJetDeepCvsL(const pat::Jet& jet)
{
	float probC = jet.bDiscriminator("pfDeepCSVJetTags:probc");
	float probL = jet.bDiscriminator("pfDeepCSVJetTags:probudsg");

	return probC/(probC+probL);
}

void ttHtautauAnalyzer::addJetQGLikelihood(std::vector<pat::Jet>& jets,
										   const edm::ValueMap<float>& QGL_valuemap)
{
	for (size_t i=0; i<jets.size(); ++i) {
		float qgLikelihood = QGL_valuemap.get(i);
		jets[i].addUserFloat("qgLikelihood", qgLikelihood);
	}
}

void ttHtautauAnalyzer::addJetQGTaggerInfo(std::vector<pat::Jet>& jets,
										   const edm::ValueMap<float>& QGL_valuemap,
										   const edm::ValueMap<float>& axis2_valuemap,
										   const edm::ValueMap<float>& ptD_valuemap,
										   const edm::ValueMap<int>& mult_valuemap)
{
	for (size_t i=0; i<jets.size(); ++i) {
		jets[i].addUserFloat("qgLikelihood", QGL_valuemap.get(i));
		jets[i].addUserFloat("axis2", axis2_valuemap.get(i));
		jets[i].addUserFloat("ptD", ptD_valuemap.get(i));
		jets[i].addUserInt("mult", mult_valuemap.get(i));
	}
}

/////////////////////
// MHT
float ttHtautauAnalyzer::getMHT(std::vector<pat::Muon>& muons,
								std::vector<pat::Electron>& electrons,
								std::vector<pat::Tau>& taus,
								std::vector<pat::Jet>& jets)
{
	float MHT_x = 0;
	float MHT_y = 0;

	for (const auto & mu : muons) {
		MHT_x -= mu.px();
		MHT_y -= mu.py();
	}

	for (const auto & ele : electrons) {
		MHT_x -= ele.px();
		MHT_y -= ele.py();
	}

	for (const auto & tau : taus) {
		MHT_x -= tau.px();
		MHT_y -= tau.py();
	}

	for (const auto & jet : jets) {
		MHT_x -= jet.px();
		MHT_y -= jet.py();
	}

	return sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
}

float ttHtautauAnalyzer::getMHT(std::vector<miniLepton>& leptons,
							    std::vector<pat::Tau>& taus,
								std::vector<pat::Jet>& jets)
{
	float MHT_x = 0;
	float MHT_y = 0;

	for (const auto & lep : leptons) {
		MHT_x -= lep.p4().Px();
		MHT_y -= lep.p4().Py();
	}

	for (const auto & tau : taus) {
		MHT_x -= tau.px();
		MHT_y -= tau.py();
	}

	for (const auto & jet : jets) {
		MHT_x -= jet.px();
		MHT_y -= jet.py();
	}

	return sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
}

/////////////////////
// MC truth matching
bool ttHtautauAnalyzer::isGenPhotonMatched(const pat::Electron& patEle, const std::vector<reco::GenParticle>& gen_particles, bool isPromptPhoton)
{
	for (const auto& gen : gen_particles) {
		if (gen.pdgId() != 22) continue;
		if (gen.status() != 1) continue;

		float dR = reco::deltaR(gen.eta(), gen.phi(), patEle.eta(), patEle.phi());
		if (dR > 0.3) continue;
		if (gen.pt() < 0.5*patEle.pt()) continue;

		if (isPromptPhoton) {
			if (not gen.isPromptFinalState()) continue;
		}

		// found match
		return true;
	}

	return false;
}
// dummy
bool ttHtautauAnalyzer::isGenPhotonMatched(const pat::Muon& patMuon, const std::vector<reco::GenParticle>& gen_particles, bool isPromptPhoton)
{
	return false;
}

const reco::GenParticle* ttHtautauAnalyzer::getMatchedGenParticle(const pat::Electron& patEle, const std::vector<reco::GenParticle>& gen_particles)
{
	const reco::GenParticle* out = NULL;
	float dRmin = 666.;
	
	// loop over genParticle collections to find match
	for (auto& gen : gen_particles) {
		if (abs(gen.pdgId()) == 11) {
			auto genStatus = gen.statusFlags();

			if (not (genStatus.isPrompt() or
					 genStatus.isDirectPromptTauDecayProduct()))
				continue;

			float dR = reco::deltaR(gen.eta(),gen.phi(),patEle.eta(),patEle.phi());
			if (dR > 0.2) continue;
			if (gen.pt() < 8.) continue;

			if (dR > dRmin) continue;  // find the closest in dR
			
			dRmin = dR;
			out = &gen;
		}
	}
	
	return out;
}

const reco::GenParticle* ttHtautauAnalyzer::getMatchedGenParticle(const pat::Muon& patMu, const std::vector<reco::GenParticle>& gen_particles)
{
	const reco::GenParticle* out = NULL;
	float dRmin = 666.;
	
	// loop over genParticle collections to find match
	for (auto& gen : gen_particles) {
		if (abs(gen.pdgId()) == 13) {
			
			auto genStatus = gen.statusFlags();

			if (not (genStatus.isPrompt() or
					 genStatus.isDirectPromptTauDecayProduct()))
				continue;

			float dR = reco::deltaR(gen.eta(),gen.phi(),patMu.eta(),patMu.phi());
			
			if (dR > 0.2) continue;
			if (gen.pt() < 8.) continue;

			if (dR > dRmin) continue;  // find the closest in dR

			dRmin = dR;
			out = &gen;
		}
	}

	return out;
}

const reco::GenParticle* ttHtautauAnalyzer::getMatchedGenParticle(const pat::Tau& patTau, const std::vector<reco::GenParticle>& gen_particles)
{
	const reco::GenParticle* out = NULL;
	float dRmin = 666.;
	
	// loop over genParticle collections to find match
	for (auto& gen : gen_particles) {
		if ( abs(gen.pdgId()) == 11 or abs(gen.pdgId()) == 13 ) {
			auto genStatus = gen.statusFlags();

			if (not (genStatus.isPrompt() or
					 genStatus.isDirectPromptTauDecayProduct()))
				continue;

			float dR = reco::deltaR(gen.eta(),gen.phi(),patTau.eta(),patTau.phi());

			if (dR > 0.2) continue;
			if (gen.pt() < 8.) continue;

		    if (dR > dRmin) continue;

			dRmin = dR;
			out = &gen;
		}

		if ( abs(gen.pdgId()) == 15 ) {
			auto genStatus = gen.statusFlags();

			if (not genStatus.isPrompt()) continue;

			reco::Candidate::LorentzVector visP4;
			for (unsigned int i = 0; i < gen.numberOfDaughters(); ++i) {
				auto daug = gen.daughter(i);
				int id = abs(daug->pdgId());
				if (id == 11 or id == 12 or id == 13 or id == 14 or id == 16)
					continue;
				visP4 += daug->p4();
			}

			float dR = reco::deltaR(visP4.eta(),visP4.phi(),patTau.eta(),patTau.phi());

			if (dR > 0.2) continue;
			// if (visP4.pt() < 15.) continue;
			if (gen.pt() < 15.) continue;

			if (dR > dRmin) continue;
			
			dRmin = dR;
			out = &gen;
		}
	}

	return out;
}

// Gen Jet matcihing for JER purpose
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
const reco::GenJet* ttHtautauAnalyzer::getMatchedGenJet
(const pat::Jet& jet, const std::vector<reco::GenJet>& genJets, double jetRes,
 double dR_max, double dPtMaxFactor)
{
  // dR_max = 0.4/2 by default for AK4 jets
  
  const reco::GenJet* matchedGenJet = nullptr;

  double dR_min = std::numeric_limits<double>::infinity();
  
  for (const auto & genjet : genJets) {
    double dR = reco::deltaR(jet, genjet);

    // pick the cloest genJet by dR
    if (dR > dR_min) continue;

    if (dR < dR_max) { // matched by dR
      // check if matched in Pt
      double dPt = std::abs(jet.pt() - genjet.pt());
      if (dPt > dPtMaxFactor * jet.pt() * jetRes) continue;

      // update using current best match
      dR_min = dR;
      matchedGenJet = &genjet;
    }
  }

  return matchedGenJet;
}

template <typename T>
int ttHtautauAnalyzer::getMCMatchType(const T& reco_particle, const reco::GenParticle* matchedGen)
{
	if (not matchedGen) return -1;

	auto genStatus = matchedGen->statusFlags();

	int mtype = -1;

	if (abs(matchedGen->pdgId()) == 11) {
		if (genStatus.isPrompt()) mtype = 1;
		if (genStatus.isDirectPromptTauDecayProduct()) mtype = 3;
	}
	else if (abs(matchedGen->pdgId()) == 13) {
		if (genStatus.isPrompt()) mtype = 2;
		if (genStatus.isDirectPromptTauDecayProduct()) mtype = 4;
	}
	else if (abs(matchedGen->pdgId()) == 15 and genStatus.isPrompt()) {
		mtype = 5;
	}
    
	return mtype;
}

template int ttHtautauAnalyzer::getMCMatchType<pat::Electron>(const pat::Electron&, const reco::GenParticle*);
template int ttHtautauAnalyzer::getMCMatchType<pat::Muon>(const pat::Muon&, const reco::GenParticle*);
template int ttHtautauAnalyzer::getMCMatchType<pat::Tau>(const pat::Tau&, const reco::GenParticle*);

int ttHtautauAnalyzer::HiggsDaughterPdgId(const std::vector<reco::GenParticle>& genParticles)
{
	for (auto & p : genParticles) {
		if (p.pdgId() != 25) continue;

		int ndaugs = p.numberOfDaughters();
		if (ndaugs != 2) continue;

		const reco::Candidate *d1 = p.daughter(0);
		const reco::Candidate *d2 = p.daughter(1);

		if ( abs(d1->pdgId()) != abs(d2->pdgId()) ) continue;

		//assert(p.statusFlags().isLastCopy());
		
		return d1->pdgId();
	}

	return -9999;
}

/////////////////////
// Vertex
bool ttHtautauAnalyzer::isGoodPV(const reco::Vertex& vertex)
{
    return (not vertex.isFake()) and (vertex.ndof() >= 4) and
		(fabs(vertex.z()) <= 24.) and fabs(vertex.position().Rho() > 2.);
}

reco::Vertex ttHtautauAnalyzer::getPrimaryVertex(edm::Handle<reco::VertexCollection> vertices)
{
	assert(vertices.isValid());

	reco::Vertex pv;

	for (reco::VertexCollection::const_iterator vtx = vertices->begin();
		 vtx != vertices->end(); ++vtx) {

		if (not isGoodPV(*vtx)) continue;

		pv = *vtx;
		break;
	}

	return pv;
}

/////////////////////
// Printout
void ttHtautauAnalyzer::dumpLeptons(const std::vector<miniLepton>& leptons)
{
	for (const auto & lep : leptons)
		lep.dump();
}

void ttHtautauAnalyzer::dumpTaus(const std::vector<miniTau>& taus)
{
	for (const auto & tau : taus)
		tau.dump();
}

void ttHtautauAnalyzer::dumpJets(const std::vector<pat::Jet>& jets)
{
	for (const auto & j : jets) {
		std::cout << " pt: " << j.pt() << " eta: " << j.eta() << " phi: " << j.phi()
				  << " energy: " << j.energy() << " csv: " << getJetCSV(j)
				  << " flavor: " << j.hadronFlavour() << std::endl;
	}
}

#endif
