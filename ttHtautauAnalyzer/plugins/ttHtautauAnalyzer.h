#ifndef ttHtautauAnalyzer_h
#define ttHtautauAnalyzer_h

/// system include files
#include <memory>
#include <stdexcept>
#include <map>
#include <iostream>

/// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

/// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TAxis.h"
#include "TFile.h"
#include "TString.h"

/// event ntuple
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/eventNtuple.h"
/// scale helper
//#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/SFHelper.h"
/// trigger helper
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/TriggerHelper.h"
/// event selector
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/EventSelector.h"

//
// class declaration
//
class ttHtautauAnalyzer : public edm::EDAnalyzer
{
 public:
	
	explicit ttHtautauAnalyzer(const edm::ParameterSet&);
	~ttHtautauAnalyzer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
	
	void beginJob() override;
	void endJob() override;
	void beginRun(const edm::Run&, const edm::EventSetup&) override;
	void endRun(const edm::Run&, const edm::EventSetup&) override;
	void analyze(const edm::Event&, const edm::EventSetup&) override;

	/// member methods
	// Setup
	void Set_up_AnaType(const std::string&);
	void Set_up_SelType(const std::string&);
	void Set_up_histograms();
	void Set_up_tree();
	
	// Object ID
	template <typename T> float ConePt(const T&);
	template <typename T> std::vector<T> getSelectedLeptons(const std::vector<T>&,
															const std::string&);
	template <typename T> void addIDFlags(std::vector<T>&,
										  edm::Handle<reco::GenParticleCollection>);
	void addIDFlagsTau(std::vector<pat::Tau>&);
	bool isLooseID(const pat::Muon&) const;
	bool isLooseID(const pat::Electron&) const;
	bool isLooseID(const pat::Tau&) const;
	bool isFakeableID(const pat::Muon&) const;
	bool isFakeableID(const pat::Electron&) const;
	bool isTightID(const pat::Muon&) const;
	bool isTightID(const pat::Tau&) const;
	bool isTightID(const pat::Electron&) const;
	bool isTightCharge(const pat::Muon&) const;
	bool isTightCharge(const pat::Electron&) const;
	bool isTightJet(const pat::Jet&) const;
	std::vector<pat::Tau> getPreselTaus(const std::vector<pat::Tau>&);
	std::vector<pat::Tau> getSelectedTaus(const std::vector<pat::Tau>&);
	std::vector<pat::Tau> getCorrectedTaus(const std::vector<pat::Tau>&, double, const std::string&);
	std::vector<pat::Jet> getSelectedJets(const std::vector<pat::Jet>&);
	std::vector<pat::Jet> getCorrectedJets(const std::vector<pat::Jet>&, JetCorrectionUncertainty*, const std::string&);
	float getJetCSV(const pat::Jet&);
	void addJetQGLikelihood(std::vector<pat::Jet>&,
							const edm::ValueMap<float>&);
	float getMHT(std::vector<pat::Muon>&, std::vector<pat::Electron>&,
				 std::vector<pat::Tau>&, std::vector<pat::Jet>&);
	float getMHT(std::vector<miniLepton>&, std::vector<pat::Tau>&,
				 std::vector<pat::Jet>&);
	
	const reco::GenParticle* getMatchedGenParticle(const pat::Muon&, const std::vector<reco::GenParticle>&);
	const reco::GenParticle* getMatchedGenParticle(const pat::Electron&, const std::vector<reco::GenParticle>&);
	const reco::GenParticle* getMatchedGenParticle(const pat::Tau&, const std::vector<reco::GenParticle>&);
	template <typename T> int getMCMatchType(const T&, const std::vector<reco::GenParticle>&);

	int HiggsDaughterPdgId(const std::vector<reco::GenParticle>&);

	// vertices
	bool isGoodPV(const reco::Vertex&);
	reco::Vertex getPrimaryVertex(edm::Handle<reco::VertexCollection>);

	void dumpLeptons(const std::vector<miniLepton>&);
	void dumpTaus(const std::vector<miniTau>&);
	void dumpJets(const std::vector<pat::Jet>&);
	
	// Event selection
	EventSelector* evt_selector_;
	bool pass_event_sel_2lss1tau(const std::vector<miniLepton>&,
								 const std::vector<miniLepton>&,
								 const std::vector<miniLepton>&,
								 const std::vector<pat::Tau>&,
								 const int, const int, const int, const float);
	bool pass_event_sel_1l2tau(const std::vector<miniLepton>&,
							   const std::vector<miniLepton>&,
							   const std::vector<miniLepton>&,
							   const std::vector<pat::Tau>&,
							   const std::vector<pat::Tau>&,
							   const int, const int, const int);
	bool pass_event_sel_3l1tau(const std::vector<miniLepton>&,
							   const std::vector<miniLepton>&,
							   const std::vector<miniLepton>&,
							   const std::vector<pat::Tau>&,
							   const int, const int, const int, const float);
	void fill_CutFlow(int, const char*);

	// Write ntuple
	/*
	void write_ntuple_bTagSF(const std::vector<pat::Jet>&);
	void write_ntuple_leptonSF(const std::vector<miniLepton>&);
	void write_ntuple_tauSF(const std::vector<pat::Tau>&);
	void write_ntuple_tauSF(const std::vector<miniTau>&);
	void write_ntuple_triggerSF(int category = -1);
	void write_ntuple_triggerSF(const miniLepton&, const std::vector<pat::Tau>&,
								bool, bool);
	void write_ntuple_triggerSF(const miniLepton&, const std::vector<miniTau>&,
								bool, bool);
	void write_ntuple_frweight(const std::vector<miniLepton>&,
							   const std::vector<pat::Tau>&);
	void write_ntuple_frweight(const std::vector<miniLepton>&,
							   const std::vector<miniTau>&);
	*/
	void write_ntuple_muons(const std::vector<pat::Muon>&);
	void write_ntuple_electrons(const std::vector<pat::Electron>&);
	void write_ntuple_taus(const std::vector<pat::Tau>&);
	void write_ntuple_jets(const std::vector<pat::Jet>&);
	
	//Utilities
	template <typename T1, typename T2>
		std::vector<T1>
		removeOverlapdR(const std::vector<T1>& v1, const std::vector<T2>& v2, double dR = 0.02);
	template <typename T> void SortByPt(T&);
    void SortByConept(std::vector<pat::Electron>&);
	void SortByConept(std::vector<pat::Muon>&);
	void SortByConept(std::vector<miniLepton>&);
	
	/// member variables
	// Analysis type
	std::string config_analysis_type_;
	Analysis_types anaType_;	
	
	// event selection region
	std::string selection_region_;
    Selection_types selType_;
	bool looseSelection_;
	
	// sample
	//std::string sample_name_;
	TString sample_name_;
	
	// generic parameters
	int verbosity_;
	bool isdata_;
	bool debug_;
	//bool doSystematics_;
	bool doSync_;
	bool event_selection_off_;
	bool doCutflow_;
	bool firstpass_;
	unsigned long long event_count_;
	
	// triggers
	TriggerHelper* trig_helper_;
	bool dumpHLT_;
	std::string hltTag_;
	std::string filterTag_;
	//bool trigger_stats_;
	//bool hltcut_off_;
	HLTConfigProvider hlt_config_;
	HLTConfigProvider filter_config_;
	
	// tauES
	std::string TESType_;  // "tauESUp" or "tauESDown" or  "NA"
	double tauES_unc_;
	
	// JEC
	std::string JECType_;
	bool doJERsmear_;
	
	double genWeightSum_;
	double genWeightxPUSum_;
	
	TTree* eventTree_;
	eventNtuple evNtuple_;
	
	// helper class
	//SFHelper* sf_helper_;
	
	// tokens
	edm::EDGetTokenT<GenEventInfoProduct> geninfo_token_;
	edm::EDGetTokenT<LHEEventProduct> lheinfo_token_;
	edm::EDGetTokenT<edm::TriggerResults> trigger_token_;
	edm::EDGetTokenT<edm::TriggerResults> filter_token_;

	edm::EDGetTokenT<reco::VertexCollection> vtx_token_;
	edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pu_token_;
	edm::EDGetTokenT<reco::BeamSpot> bs_token_;
	edm::EDGetTokenT<double> rho_token_;

	edm::EDGetTokenT<pat::ElectronCollection> electrons_token_;
	edm::EDGetTokenT<pat::MuonCollection> muons_token_;
	edm::EDGetTokenT<pat::TauCollection> taus_token_;
	edm::EDGetTokenT<pat::JetCollection> jets_token_;
	edm::EDGetTokenT<pat::METCollection> mets_token_;
	edm::EDGetTokenT<reco::GenParticleCollection> genparticle_token_;
	edm::EDGetTokenT<reco::GenJetCollection> genjets_token_;
	edm::EDGetTokenT<edm::ValueMap<float>> qg_token_;

	//edm::EDGetTokenT<edm::View<reco::Muon>> badMuons_token_;
	//edm::EDGetTokenT<edm::View<reco::Muon>> clonedMuons_token_;

	/// Output file is opened/closed through CMS py config
	edm::Service<TFileService> fs_;
	
	// constants
	std::string sysList_[16] =
		{"LFUp","LFDown","HFUp","HFDown",
		 "HFStats1Up","HFStats1Down","HFStats2Up","HFStats2Down",
		 "LFStats1Up","LFStats1Down","LFStats2Up","LFStats2Down",
		 "cErr1Up","cErr1Down","cErr2Up","cErr2Down"};
	
	// CSV WP
	const double csv_loose_wp_;
	const double csv_medium_wp_;
	const double csv_tight_wp_;

	// histograms
	// event count
	TH1I* h_nProcessed_;
	TH1D* h_SumGenWeight_;
	TH1D* h_SumGenWeightxPU_;
	// cut flow
	TH1D* h_CutFlow_;
	// pileup
	TH1D* h_npuTrue_;
	TH1D* h_npuInTime_;
	// HLT and filters
	TH1I* h_hlt_;
	TH1I* h_flt_;
};

template <typename T1, typename T2>
	std::vector<T1>
	ttHtautauAnalyzer::removeOverlapdR(const std::vector<T1>& v1, const std::vector<T2>& v2, double dR)
{
	std::vector<T1> res;
	for (const auto& o1 : v1) {
		bool keep = true;
		for (const auto& o2 : v2) {
			if (reco::deltaR(o1.eta(),o1.phi(),o2.eta(),o2.phi()) < dR)
				keep = false;
		}
		if (keep) res.push_back(o1);
	}

	return res;
}

template <typename T> void ttHtautauAnalyzer::SortByPt( T& vec)
{
	std::sort(vec.begin(), vec.end(), [] (typename T::value_type a, typename T::value_type b) {return a.pt() > b.pt();});
}

#endif
