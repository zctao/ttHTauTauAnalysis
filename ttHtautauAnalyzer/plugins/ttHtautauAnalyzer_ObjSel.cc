#ifndef ttHtautauAnalyzer_ObjSel_cc
#define ttHtautauAnalyzer_ObjSel_cc

#include "ttHTauTauAnalysis/ttHtautauAnalyzer/plugins/ttHtautauAnalyzer.h"

///////////////////
// leptons
template <typename T> float ttHtautauAnalyzer::ConePt(const T &lep)
{
    // return (lep.userFloat("leptonMVA") > 0.75) ? lep.pt() :
    //	(0.85 * lep.pt() / lep.userFloat("nearestJetPtRatio"));

    return lep.userFloat("correctedPt");
}
template float ttHtautauAnalyzer::ConePt<pat::Electron>(const pat::Electron &);
template float ttHtautauAnalyzer::ConePt<pat::Muon>(const pat::Muon &);

void ttHtautauAnalyzer::SortByConept(std::vector<pat::Electron> &electrons)
{
    std::sort(electrons.begin(), electrons.end(),
              [](pat::Electron e1, pat::Electron e2) {
                  return e1.userFloat("conePt") > e2.userFloat("conePt");
              });
}

void ttHtautauAnalyzer::SortByConept(std::vector<pat::Muon> &muons)
{
    std::sort(muons.begin(), muons.end(), [](pat::Muon m1, pat::Muon m2) {
        return m1.userFloat("conePt") > m2.userFloat("conePt");
    });
}

void ttHtautauAnalyzer::SortByConept(std::vector<miniLepton> &leptons)
{
    std::sort(leptons.begin(), leptons.end(), [](miniLepton l1, miniLepton l2) {
        return l1.conept() > l2.conept();
    });
}

// get selected collection
template <typename T>
std::vector<T>
ttHtautauAnalyzer::getSelectedLeptons(const std::vector<T> &leptons,
                                      const std::string &selection)
{
    std::vector<T> selected;
    for (auto &lep : leptons) {
        if (lep.hasUserInt(selection)) {
            if (not lep.userInt(selection))
                continue;
        } else {
            if (selection == "isLoose") {
                if (not isLooseID(lep))
                    continue;
            } else if (selection == "isFakeable") {
                if (not isFakeableID(lep))
                    continue;
            } else if (selection == "isTight") {
                if (not isTightID(lep))
                    continue;
            } else {
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
template std::vector<pat::Muon>
ttHtautauAnalyzer::getSelectedLeptons(const std::vector<pat::Muon> &,
                                      const std::string &);
template std::vector<pat::Electron>
ttHtautauAnalyzer::getSelectedLeptons(const std::vector<pat::Electron> &,
                                      const std::string &);

// add ID related flags/variables
template <typename T>
void ttHtautauAnalyzer::addIDFlags(
    std::vector<T> &leptons,
    edm::Handle<reco::GenParticleCollection> MC_particles)
{
    for (auto &lep : leptons) {
        lep.addUserInt("isLoose", isLooseID(lep));
        lep.addUserInt("isFakeable", isFakeableID(lep));
        lep.addUserInt("isTight", isTightID(lep));
        lep.addUserInt("isTightCharge", isTightCharge(lep));
        lep.addUserFloat("conePt", ConePt(lep));
        // MC matching
        if (not isdata_)
            lep.addUserInt("MCMatchType", getMCMatchType(lep, *MC_particles));
    }
}
template void
ttHtautauAnalyzer::addIDFlags(std::vector<pat::Electron> &,
                              edm::Handle<reco::GenParticleCollection>);
template void
ttHtautauAnalyzer::addIDFlags(std::vector<pat::Muon> &,
                              edm::Handle<reco::GenParticleCollection>);

///////////////////
// muons
bool ttHtautauAnalyzer::isLooseID(const pat::Muon &mu) const
{
    bool passKinematic = mu.pt() > 5. and std::abs(mu.eta()) < 2.4;

    return passKinematic and mu.userFloat("idPreselection") > 0.5;
}

bool ttHtautauAnalyzer::isFakeableID(const pat::Muon &mu) const
{
    return mu.userFloat("idFakeable") > 0.5;
}

bool ttHtautauAnalyzer::isTightID(const pat::Muon &mu) const
{
    return mu.userFloat("idMVABased") > 0.5;
}

bool ttHtautauAnalyzer::isTightCharge(const pat::Muon &mu) const
{
    bool tightcharge = false;

    if (mu.innerTrack().isAvailable())
        tightcharge = mu.innerTrack()->ptError() / mu.innerTrack()->pt() < 0.2;

    return tightcharge;
}

///////////////////
// electrons
bool ttHtautauAnalyzer::isLooseID(const pat::Electron &ele) const
{
    bool passKinematic = ele.pt() > 7. and fabs(ele.eta()) < 2.5;
    // bool passPhotonVeto = (ele.passConversionVeto() and
    //					   ele.userFloat("numMissingHits") <= 1);

    return (passKinematic and ele.userFloat("idPreselection") > 0.5);
}

bool ttHtautauAnalyzer::isFakeableID(const pat::Electron &ele) const
{
    return (ele.userFloat("idFakeable") > 0.5 // and
            // ele.userFloat("numMissingHits") == 0 and
            // ele.passConversionVeto()
            );
}

bool ttHtautauAnalyzer::isTightID(const pat::Electron &ele) const
{
    return (ele.userFloat("idMVABased") > 0.5 // and
            // ele.userFloat("numMissingHits") == 0 and
            // ele.passConversionVeto()
            );
}

bool ttHtautauAnalyzer::isTightCharge(const pat::Electron &ele) const
{
    bool tightcharge =
        ele.isGsfCtfScPixChargeConsistent() + ele.isGsfScPixChargeConsistent() >
        1;

    // std::cout << "elel tight charge : " <<
    // ele.isGsfCtfScPixChargeConsistent()
    //		  << " " << ele.isGsfScPixChargeConsistent() << std::endl;

    return tightcharge;
}

///////////////////
// taus

std::vector<pat::Tau>
ttHtautauAnalyzer::getPreselTaus(const std::vector<pat::Tau> &taus)
{
    std::vector<pat::Tau> preseltaus;
    for (auto &tau : taus) {

        if (tau.hasUserInt("isLoose")) {
            if (not tau.userInt("isLoose"))
                continue;
        } else {
            if (not isLooseID(tau))
                continue;
        }

        preseltaus.push_back(tau);
    }

    return preseltaus;
}

std::vector<pat::Tau>
ttHtautauAnalyzer::getSelectedTaus(const std::vector<pat::Tau> &taus)
{
    std::vector<pat::Tau> seltaus;

    for (auto &tau : taus) {

        if (tau.hasUserInt("isTight")) {
            if (not tau.userInt("isTight"))
                continue;
        } else {
            if (not isTightID(tau))
                continue;
        }

        seltaus.push_back(tau);
    }

    return seltaus;
}

std::vector<pat::Tau>
ttHtautauAnalyzer::getCorrectedTaus(const std::vector<pat::Tau> &input_taus,
                                    double tauES_Unc,
                                    const std::string &TESType)
{
    float shift = 0.;

    if (TESType == "tauESUp")
        shift = 1.;
    else if (TESType == "tauESDown")
        shift = -1.;
    else
        return input_taus;

    std::vector<pat::Tau> output_taus = input_taus;

    for (auto &tau : output_taus) {
        auto corrP4 = tau.p4();
        corrP4 *= 1 + shift * tauES_Unc;
        tau.setP4(corrP4);
    }

    return output_taus;
}

// better to move this to leptonID package
bool ttHtautauAnalyzer::isLooseID(const pat::Tau &tau) const
{
    bool passKinematic = tau.pt() > 20. and std::abs(tau.eta()) < 2.3;

    // tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT")
    bool passID = tau.userFloat("idPreselection") > 0.5;

    return (passKinematic and passID);
}

bool ttHtautauAnalyzer::isTightID(const pat::Tau &tau) const
{

    // loose ID is the prerequisite for tight ID
    if (tau.hasUserInt("isLoose")) {
        if (not tau.userInt("isLoose"))
            return false;
    } else {
        if (not isLooseID(tau))
            return false;
    }

    if (anaType_ == Analyze_1l2tau) {
        // tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT")
        return tau.tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT") > 0.5;
    } else {
        // tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT")
        return tau.userFloat("idSelection") > 0.5;
    }
}

void ttHtautauAnalyzer::addIDFlagsTau(std::vector<pat::Tau> &taus)
{
    for (auto &tau : taus) {
        tau.addUserInt("isLoose", isLooseID(tau));
        tau.addUserInt("isTight", isTightID(tau));
    }
}

///////////////////
// jets
bool ttHtautauAnalyzer::isTightJet(const pat::Jet &jet) const
{
    // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017#Preliminary_Recommendations_for

    bool istight = false;

    if (fabs(jet.eta()) <= 2.7) {
        istight = jet.neutralHadronEnergyFraction() < 0.90 and
                  jet.neutralEmEnergyFraction() < 0.90 and
                  jet.numberOfDaughters() > 1;
        if (fabs(jet.eta()) <= 2.4) {
            istight = istight and jet.chargedHadronEnergyFraction() > 0 and
                      jet.chargedMultiplicity() > 0;
        }
    } else if (fabs(jet.eta()) <= 3.0) {
        istight = jet.neutralEmEnergyFraction() < 0.99 and
                  jet.neutralEmEnergyFraction() > 0.02 and
                  jet.neutralMultiplicity() > 2;
    } else {
        istight = jet.neutralEmEnergyFraction() < 0.90 and
                  jet.neutralHadronEnergyFraction() > 0.02 and
                  jet.neutralMultiplicity() > 10;
    }

    return istight;
}

std::vector<pat::Jet>
ttHtautauAnalyzer::getSelectedJets(const std::vector<pat::Jet> &jets)
{
    std::vector<pat::Jet> seljets;
    for (const auto &jet : jets) {
        bool passKinematic = (jet.pt() > 25. and fabs(jet.eta()) < 2.4);

        bool passID = isTightJet(jet);

        /*
        bool looseID =
            jet.neutralHadronEnergyFraction() < 0.99 and
            jet.chargedEmEnergyFraction() < 0.99 and
            jet.neutralEmEnergyFraction() < 0.99 and
            jet.numberOfDaughters() > 1 and
            jet.chargedHadronEnergyFraction() > 0.0 and
            jet.chargedMultiplicity() > 0;
        */
        /*
        bool looseID =
            (jet.neutralHadronEnergyFraction()<0.99 and
             jet.neutralEmEnergyFraction()<0.99 and
             jet.chargedMultiplicity()+jet.neutralMultiplicity()>1
             )
            and
            ((std::abs(jet.eta())<=2.4 and
              jet.chargedHadronEnergyFraction()>0 and
              jet.chargedMultiplicity()>0 && jet.chargedEmEnergyFraction()<0.99)
             or std::abs(jet.eta())>2.4
             );
        */
        if (passKinematic and passID)
            seljets.push_back(jet);
    }
    return seljets;
}

std::vector<pat::Jet>
ttHtautauAnalyzer::getCorrectedJets(const std::vector<pat::Jet> &input_jets,
                                    JetCorrectionUncertainty *jecUnc,
                                    const std::string &JECType)
{
    float shift = 0.;

    if (JECType == "JESUp")
        shift = 1.;
    else if (JECType == "JESDown")
        shift = -1.;
    else
        return input_jets;

    std::vector<pat::Jet> output_jets = input_jets;

    for (auto &jet : output_jets) {
        jecUnc->setJetEta(jet.eta());
        jecUnc->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
        float unc = jecUnc->getUncertainty(true);
        float jes = 1 + shift * unc;
        jet.scaleEnergy(jes);
    }

    return output_jets;
}

float ttHtautauAnalyzer::getJetCSV(const pat::Jet &jet)
{
    float defaultvalue = -0.1;
    // float csv =
    // jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    float csv = jet.bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll");
    if (std::isnan(csv))
        return defaultvalue;
    if (csv > 1.)
        return 1.;
    if (csv < 0.)
        return defaultvalue;

    return csv;
}

/////////////////////
// MHT
float ttHtautauAnalyzer::getMHT(std::vector<pat::Muon> &muons,
                                std::vector<pat::Electron> &electrons,
                                std::vector<pat::Tau> &taus,
                                std::vector<pat::Jet> &jets)
{
    float MHT_x = 0;
    float MHT_y = 0;

    for (const auto &mu : muons) {
        MHT_x -= mu.px();
        MHT_y -= mu.py();
    }

    for (const auto &ele : electrons) {
        MHT_x -= ele.px();
        MHT_y -= ele.py();
    }

    for (const auto &tau : taus) {
        MHT_x -= tau.px();
        MHT_y -= tau.py();
    }

    for (const auto &jet : jets) {
        MHT_x -= jet.px();
        MHT_y -= jet.py();
    }

    return sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
}

float ttHtautauAnalyzer::getMHT(std::vector<miniLepton> &leptons,
                                std::vector<pat::Tau> &taus,
                                std::vector<pat::Jet> &jets)
{
    float MHT_x = 0;
    float MHT_y = 0;

    for (const auto &lep : leptons) {
        MHT_x -= lep.p4().Px();
        MHT_y -= lep.p4().Py();
    }

    for (const auto &tau : taus) {
        MHT_x -= tau.px();
        MHT_y -= tau.py();
    }

    for (const auto &jet : jets) {
        MHT_x -= jet.px();
        MHT_y -= jet.py();
    }

    return sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
}

/////////////////////
// MC truth matching

const reco::GenParticle *ttHtautauAnalyzer::getMatchedGenParticle(
    const pat::Muon &patMu, const std::vector<reco::GenParticle> &gen_particles)
{
    const reco::GenParticle *out = NULL;
    float dRmin = 666.;

    // loop over genParticle collections to find match
    for (auto &gen : gen_particles) {
        if (abs(gen.pdgId()) == 13) {

            auto genStatus = gen.statusFlags();

            if (not(genStatus.isPrompt() or
                    genStatus.isDirectPromptTauDecayProduct()))
                continue;

            float dR =
                reco::deltaR(gen.eta(), gen.phi(), patMu.eta(), patMu.phi());

            if (dR > 0.2)
                continue;
            if (gen.pt() < 8.)
                continue;

            if (dR > dRmin)
                continue; // find the closest in dR

            dRmin = dR;
            out = &gen;
        }
    }

    return out;
}

const reco::GenParticle *ttHtautauAnalyzer::getMatchedGenParticle(
    const pat::Electron &patEle,
    const std::vector<reco::GenParticle> &gen_particles)
{
    const reco::GenParticle *out = NULL;
    float dRmin = 666.;

    // loop over genParticle collections to find match
    for (auto &gen : gen_particles) {
        if (abs(gen.pdgId()) == 11) {
            auto genStatus = gen.statusFlags();

            if (not(genStatus.isPrompt() or
                    genStatus.isDirectPromptTauDecayProduct()))
                continue;

            float dR =
                reco::deltaR(gen.eta(), gen.phi(), patEle.eta(), patEle.phi());
            if (dR > 0.2)
                continue;
            if (gen.pt() < 8.)
                continue;

            if (dR > dRmin)
                continue; // find the closest in dR

            dRmin = dR;
            out = &gen;
        }
    }

    return out;
}

const reco::GenParticle *ttHtautauAnalyzer::getMatchedGenParticle(
    const pat::Tau &patTau, const std::vector<reco::GenParticle> &gen_particles)
{
    const reco::GenParticle *out = NULL;
    float dRmin = 666.;

    // loop over genParticle collections to find match
    for (auto &gen : gen_particles) {
        if (abs(gen.pdgId()) == 11 or abs(gen.pdgId()) == 13) {
            auto genStatus = gen.statusFlags();

            if (not(genStatus.isPrompt() or
                    genStatus.isDirectPromptTauDecayProduct()))
                continue;

            float dR =
                reco::deltaR(gen.eta(), gen.phi(), patTau.eta(), patTau.phi());

            if (dR > 0.2)
                continue;
            if (gen.pt() < 8.)
                continue;

            if (dR > dRmin)
                continue;

            dRmin = dR;
            out = &gen;
        }

        if (abs(gen.pdgId()) == 15) {
            auto genStatus = gen.statusFlags();

            if (not genStatus.isPrompt())
                continue;

            reco::Candidate::LorentzVector visP4;
            for (unsigned int i = 0; i < gen.numberOfDaughters(); ++i) {
                auto daug = gen.daughter(i);
                int id = abs(daug->pdgId());
                if (id == 11 or id == 12 or id == 13 or id == 14 or id == 16)
                    continue;
                visP4 += daug->p4();
            }

            float dR = reco::deltaR(visP4.eta(), visP4.phi(), patTau.eta(),
                                    patTau.phi());

            if (dR > 0.2)
                continue;
            // if (visP4.pt() < 15.) continue;
            if (gen.pt() < 15.)
                continue;

            if (dR > dRmin)
                continue;

            dRmin = dR;
            out = &gen;
        }
    }

    return out;
}

template <typename T>
int ttHtautauAnalyzer::getMCMatchType(
    const T &reco_particle, const std::vector<reco::GenParticle> &gen_particles)
{
    auto matchedGen = getMatchedGenParticle(reco_particle, gen_particles);

    if (not matchedGen)
        return -1;

    auto genStatus = matchedGen->statusFlags();

    int mtype = -1;

    if (abs(matchedGen->pdgId()) == 11) {
        if (genStatus.isPrompt())
            mtype = 1;
        if (genStatus.isDirectPromptTauDecayProduct())
            mtype = 3;
    } else if (abs(matchedGen->pdgId()) == 13) {
        if (genStatus.isPrompt())
            mtype = 2;
        if (genStatus.isDirectPromptTauDecayProduct())
            mtype = 4;
    } else if (abs(matchedGen->pdgId()) == 15 and genStatus.isPrompt()) {
        mtype = 5;
    }

    return mtype;
}

template int ttHtautauAnalyzer::getMCMatchType<pat::Electron>(
    const pat::Electron &, const std::vector<reco::GenParticle> &);
template int ttHtautauAnalyzer::getMCMatchType<pat::Muon>(
    const pat::Muon &, const std::vector<reco::GenParticle> &);
template int ttHtautauAnalyzer::getMCMatchType<pat::Tau>(
    const pat::Tau &, const std::vector<reco::GenParticle> &);

int ttHtautauAnalyzer::HiggsDaughterPdgId(
    const std::vector<reco::GenParticle> &genParticles)
{
    for (auto &p : genParticles) {
        if (p.pdgId() != 25)
            continue;

        int ndaugs = p.numberOfDaughters();
        if (ndaugs != 2)
            continue;

        const reco::Candidate *d1 = p.daughter(0);
        const reco::Candidate *d2 = p.daughter(1);

        if (abs(d1->pdgId()) != abs(d2->pdgId()))
            continue;

        // assert(p.statusFlags().isLastCopy());

        return d1->pdgId();
    }

    return -9999;
}

/////////////////////
// Vertex
bool ttHtautauAnalyzer::isGoodPV(const reco::Vertex &vertex)
{
    return (not vertex.isFake()) and (vertex.ndof() >= 4) and
           (fabs(vertex.z()) <= 24.) and fabs(vertex.position().Rho() > 2.);
}

reco::Vertex ttHtautauAnalyzer::getPrimaryVertex(
    edm::Handle<reco::VertexCollection> vertices)
{
    assert(vertices.isValid());

    reco::Vertex pv;

    for (reco::VertexCollection::const_iterator vtx = vertices->begin();
         vtx != vertices->end(); ++vtx) {

        if (not isGoodPV(*vtx))
            continue;

        pv = *vtx;
        break;
    }

    return pv;
}

/////////////////////
// Printout
void ttHtautauAnalyzer::dumpLeptons(const std::vector<miniLepton> &leptons)
{
    for (const auto &lep : leptons)
        lep.dump();
}

void ttHtautauAnalyzer::dumpTaus(const std::vector<miniTau> &taus)
{
    for (const auto &tau : taus)
        tau.dump();
}

void ttHtautauAnalyzer::dumpJets(const std::vector<pat::Jet> &jets)
{
    for (const auto &j : jets) {
        std::cout << " pt: " << j.pt() << " eta: " << j.eta()
                  << " phi: " << j.phi() << " energy: " << j.energy()
                  << " csv: " << getJetCSV(j)
                  << " flavor: " << j.hadronFlavour() << std::endl;
    }
}

#endif
