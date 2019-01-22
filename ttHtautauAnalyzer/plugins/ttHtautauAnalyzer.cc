// -*- C++ -*-
//
// Package:    ttHTauTauAnalysis/ttHtautauAnalyzer
// Class:      ttHtautauAnalyzer
// 
/**\class ttHtautauAnalyzer ttHtautauAnalyzer.cc ttHTauTauAnalysis/ttHtautauAnalyzer/plugins/ttHtautauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhengcheng Tao
//         Created:  Fri, 10 Mar 2017 17:15:13 GMT
//
//


// user include files
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/plugins/ttHtautauAnalyzer.h"

//
// constants, enums and typedefs
//
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/Types_enum.h"

//
// static data member definitions
//

//
// constructors and destructor
//
ttHtautauAnalyzer::ttHtautauAnalyzer(const edm::ParameterSet& iConfig):
	// Loose selection for more statistics
	looseSelection_ (iConfig.getParameter<bool>("looseSelection")),
	// Sample
	sample_name_ (iConfig.getParameter<string>("sample_name")),
	// Generic parameters
	verbosity_ (iConfig.getParameter<int>("verbosity")),
	isdata_ (iConfig.getParameter<bool>("using_collision_data")),
	debug_ (iConfig.getParameter<bool>("debug_mode")),
	//doSystematics_ (iConfig.getParameter<bool>("do_systematics")),
	doSync_ (iConfig.getParameter<bool>("do_sync")),
	event_selection_off_ (iConfig.getParameter<bool>("turn_off_event_sel")),
	doCutflow_ (iConfig.getParameter<bool>("doCutFlow")),
	// triggers
	dumpHLT_ (iConfig.getParameter<bool>("print_HLT_event_path")),
	hltTag_ (iConfig.getParameter<string>("HLT_config_tag")),
	filterTag_ (iConfig.getParameter<string>("filter_config_tag")),
	// tau energy scale uncerainty
	tauES_unc_ (iConfig.getParameter<double>("TauESUnc")),
	// CSV WP
	csv_loose_wp_ (iConfig.getParameter<double>("csv_loose_wp")),
	csv_medium_wp_ (iConfig.getParameter<double>("csv_medium_wp")),
	csv_tight_wp_ (iConfig.getParameter<double>("csv_tight_wp"))
{
	
	// register data access
	geninfo_token_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
	lheinfo_token_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
	trigger_token_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", iConfig.getParameter<std::string>("HLT_config_tag")));
	filter_token_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", iConfig.getParameter<std::string>("filter_config_tag")));
	vtx_token_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pv"));
	pu_token_ = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileup"));
	bs_token_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"));
	rho_token_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
	electrons_token_ = consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"));
	muons_token_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
	taus_token_ = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"));
	jets_token_ = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));
	mets_token_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"));
	genparticle_token_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedgen"));
	genjets_token_ = consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJets"));
	// quark-gluon likelihood
	qg_token_ = consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"));
	axis2_token_ = consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"));
	ptD_token_ = consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"));
	mult_token_ = consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult"));
	
	//badMuons_token_ = consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("badmu"));
	//clonedMuons_token_ = consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("clonemu"));

	// Trigger and filter
	// always include all provided trigger paths regardless of analysis type
	trig_helper_ = new TriggerHelper(Analysis_types::Analyze_inclusive, verbosity_);
	
	// histograms
	Set_up_histograms();
	// trees
	Set_up_tree();
	// event selection
	evt_selector_ = new EventSelector(debug_, not isdata_);
}


ttHtautauAnalyzer::~ttHtautauAnalyzer()
{
	
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

	delete trig_helper_;
	delete evt_selector_;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ttHtautauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if (debug_) {
		std::cout << "****************************" << std::endl;
		std::cout << "start of analyze() " << std::endl;
		std::cout << iEvent.id().run() << ":" << iEvent.id().luminosityBlock()
				  << ":" << iEvent.id().event() << std::endl;
	}
	
	using namespace edm;

	h_nProcessed_->Fill(1);
	++event_count_;
	
	//// set up handles
	Handle<edm::TriggerResults> triggerResults;
	iEvent.getByToken(trigger_token_,triggerResults);
	Handle<edm::TriggerResults> filterResults;
	iEvent.getByToken(filter_token_,filterResults);
	Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vtx_token_,vertices);
	Handle<std::vector<PileupSummaryInfo>> PU_info;
	iEvent.getByToken(pu_token_,PU_info);
	Handle<reco::BeamSpot> BS;
	iEvent.getByToken(bs_token_,BS);
	Handle<double> Rho;
	iEvent.getByToken(rho_token_,Rho);
	Handle<pat::ElectronCollection> electrons;
	iEvent.getByToken(electrons_token_,electrons);
	Handle<pat::MuonCollection> muons;
	iEvent.getByToken(muons_token_,muons);
	Handle<pat::TauCollection> taus;
	iEvent.getByToken(taus_token_,taus);
	Handle<pat::JetCollection> jets;
	iEvent.getByToken(jets_token_,jets);
	Handle<pat::METCollection> mets;
	iEvent.getByToken(mets_token_,mets);
	Handle<edm::ValueMap<float>> qgHandle;
	iEvent.getByToken(qg_token_, qgHandle);
	Handle<edm::ValueMap<float>> axis2Handle;
	iEvent.getByToken(axis2_token_, axis2Handle);
	Handle<edm::ValueMap<float>> ptDHandle;
	iEvent.getByToken(ptD_token_, ptDHandle);
	Handle<edm::ValueMap<int>> multHandle;
	iEvent.getByToken(mult_token_, multHandle);

	Handle<GenEventInfoProduct> event_gen_info;
	Handle<LHEEventProduct> event_lhe_info;
	Handle<reco::GenParticleCollection> MC_particles;
	Handle<reco::GenJetCollection> genJets;
	//Handle<View<reco::Muon>> badMuons;
	//Handle<View<reco::Muon>> clonedMuons;
	if (not isdata_) {
		iEvent.getByToken(geninfo_token_,event_gen_info);
		iEvent.getByToken(lheinfo_token_,event_lhe_info);
		iEvent.getByToken(genparticle_token_,MC_particles);		
		iEvent.getByToken(genjets_token_,genJets);
		//iEvent.getByToken(badMuons_token_,badMuons);
		//iEvent.getByToken(clonedMuons_token_,clonedMuons);
	}

	/////////////////////////////////////////
	// gen weights
	float mc_weight = -9999.;
	float mc_weight_scale_muF0p5 = -9999.;
	float mc_weight_scale_muF2 = -9999.;
	float mc_weight_scale_muR0p5 = -9999.;
	float mc_weight_scale_muR2 = -9999.;
	if (not isdata_) {
		float genWeight = event_gen_info.product()->weight();
		mc_weight = genWeight;// / std::abs(genWeight);

		if (event_lhe_info.isValid()) {
			// tHq, tHW samples: nominal weight is inverted
			// retrieve SM weight:
			if (sample_name_.EqualTo("thq",TString::ECaseCompare::kIgnoreCase)){
				mc_weight *= (event_lhe_info->weights()[893].wgt)/(event_lhe_info->originalXWGTUP());
			}
			else if (sample_name_.EqualTo("thw",TString::ECaseCompare::kIgnoreCase)){
				mc_weight *= (event_lhe_info->weights()[1091].wgt)/(event_lhe_info->originalXWGTUP());
			}

			if (event_lhe_info->weights().size() > 6) {
				mc_weight_scale_muF0p5 = // muF = 0.5 | muR = 1
					mc_weight * (event_lhe_info->weights()[2].wgt)/
					(event_lhe_info->originalXWGTUP());
			    mc_weight_scale_muF2 = // muF = 2 | muR = 1
				    mc_weight * (event_lhe_info->weights()[1].wgt)/
					(event_lhe_info->originalXWGTUP());
			    mc_weight_scale_muR0p5 = // muF = 1 | muR = 0.5
				    mc_weight * (event_lhe_info->weights()[6].wgt)/
					(event_lhe_info->originalXWGTUP());
			    mc_weight_scale_muR2 = // muF = 1 | muR = 2
				    mc_weight * (event_lhe_info->weights()[3].wgt)/
					(event_lhe_info->originalXWGTUP());
			}
		}
	}
	
	/////////////////////////////////////////
	// event counts before cuts
	if (not isdata_) {
		h_SumGenWeight_->Fill(1, mc_weight);
		//h_SumGenWeightxPU_->Fill(1, mc_weight * pu_weight);
	}

	/////////////////////////////////////////
	// Pile up
	float npuTrue = -9999.;
	float npuInTime = -9999.;
	//float pu_weight = -9999.;
	if (not isdata_) {
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		for (PVI = PU_info->begin(); PVI != PU_info->end(); ++PVI){
			int BX = PVI->getBunchCrossing();
			if (BX == 0) { // "0" is the in-time crossing, negative values are the early crossings, positive are late
				npuTrue = PVI -> getTrueNumInteractions();
				npuInTime = PVI -> getPU_NumInteractions();
				break;
			}
		}

		h_npuTrue_->Fill(npuTrue);
		h_npuTrue_genweight_->Fill(npuTrue, mc_weight);
		h_npuInTime_->Fill(npuInTime);
		h_npuInTime_genweight_->Fill(npuInTime, mc_weight);
		
		// pileup weight
	    //pu_weight = sf_helper_->Get_PUWeight(npuTrue);
	}
	
	/////////////////////////////////////////	
	// leptons
	std::vector<miniLepton> lep_loose;
	std::vector<miniLepton> lep_fakeable;
	std::vector<miniLepton> lep_tight;
	
	/////////////////////////////////////////
	// Muons
	if (debug_) std::cout << "muons" << std::endl;
	
	// loose
	std::vector<pat::Muon> mu_preselected =
		getSelectedLeptons(*muons, std::string("isLoose"));
	//SortByPt(mu_preselected);  // sort by pt
	// add ID flags to preselected muons
	addIDFlags(mu_preselected, MC_particles);
	SortByConept(mu_preselected);  // sort by conept

	unsigned int n_muon_fakeable = 0;
	unsigned int n_muon_tight = 0;

	for (const auto & mu : mu_preselected) {
		miniLepton l(mu);
		lep_loose.push_back(l);

		if (l.passFakeableSel()) {
			lep_fakeable.push_back(l);
			++n_muon_fakeable;
		}

		if (l.passTightSel()) {
			lep_tight.push_back(l);
			++n_muon_tight;
		}
	}
	assert(n_muon_tight<=n_muon_fakeable);
	
	if (debug_) {
		std::cout << "n_muon_loose : " << mu_preselected.size() << std::endl;
		std::cout << "n_muon_fakeable : " << n_muon_fakeable << std::endl;
		std::cout << "n_muon_tight : " << n_muon_tight << std::endl;
	}
	
	/////////////////////////////////////////
	// Electrons
	if (debug_) std::cout << "electrons" << std::endl;
	
	// loose
	std::vector<pat::Electron> ele_preselected =
		getSelectedLeptons(*electrons, std::string("isLoose"));
	// remove overlap with muons
	ele_preselected = removeOverlapdR(ele_preselected, mu_preselected, 0.05);
	//SortByPt(ele_preselected);   // sort by pt
	// add ID flags to preselected electrons
	addIDFlags(ele_preselected, MC_particles);
	SortByConept(ele_preselected);  // sort by conept

	unsigned int n_ele_fakeable = 0;
	unsigned int n_ele_tight = 0;
	
	for (const auto & ele : ele_preselected) {
		miniLepton l(ele);
		lep_loose.push_back(l);

		if (l.passFakeableSel()) {
			lep_fakeable.push_back(l);
			++n_ele_fakeable;
		}

		if (l.passTightSel()) {
			lep_tight.push_back(l);
			++n_ele_tight;
		}
	}	
	assert(n_ele_tight<=n_ele_fakeable);
	
	if (debug_) {
		std::cout << "n_electrons_loose : " << ele_preselected.size() << std::endl;
		std::cout << "n_electrons_fakeable : " << n_ele_fakeable << std::endl;
		std::cout << "n_electrons_tight : " << n_ele_tight << std::endl;
	}
	
	//sort leptons by conept
	SortByConept(lep_loose);	
	SortByConept(lep_fakeable);
	SortByConept(lep_tight);

	if (debug_) {
		std::cout << "n_lep_loose : " << lep_loose.size() << std::endl;
		std::cout << "n_lep_fakeable : " << lep_fakeable.size() << std::endl;
		std::cout << "n_lep_tight : " << lep_tight.size() << std::endl;
		std::cout << "loose leptons : " << std::endl;
		dumpLeptons(lep_loose);
	}
	
	/////////////////////////////////////////
	// Taus
	if (debug_) std::cout << "taus" << std::endl;

	// loose
	std::vector<pat::Tau> tau_loose = getPreselTaus(*taus);
	// remove overlap with muons and electrons
	tau_loose = removeOverlapdR(tau_loose, mu_preselected, 0.3);
	tau_loose = removeOverlapdR(tau_loose, ele_preselected, 0.3);
	addIDFlagsTau(tau_loose);
	// sort by pt
	SortByPt(tau_loose);
	// add mc match type if not data
	if (not isdata_) {
		for (auto & tau : tau_loose)
			tau.addUserInt("MCMatchType", getMCMatchType(tau, *MC_particles));
	}

	std::vector<miniTau> minitau_loose;
	std::vector<miniTau> minitau_tight;
	for (const auto & tau : tau_loose) {
		miniTau mt(tau);
		minitau_loose.push_back(mt);

		if (mt.passTightSel())
			minitau_tight.push_back(mt);
	}

	if (debug_) {
		std::cout << "n_tau_loose : " << minitau_loose.size() << std::endl;
		std::cout << "n_tau : " << minitau_tight.size() << std::endl;
		std::cout << "loose taus : " << std::endl;
		dumpTaus(minitau_loose);
	}
	
	/////////////////////////////////////////
	// Jets
	if (debug_) std::cout << "jets" << std::endl;

    std::vector<pat::Jet> jet_raw = *jets;
    
	if (not isdata_) { // correction factors for simulated jet energy

        // JES
        edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
		iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
		const JetCorrectorParameters & JetCorPar = (*JetCorParColl)["Uncertainty"];
		JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

        // JER
        JME::JetResolution resolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
        JME::JetResolutionScaleFactor res_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");

        // random number generator for smearing
        std::mt19937 random_generator;
        unsigned int runNum_uint = static_cast<unsigned int>(iEvent.id().run());
        unsigned int lumiNum_uint = static_cast<unsigned int>(iEvent.id().luminosityBlock());
        unsigned int evNum_uint = static_cast<unsigned int>(iEvent.id().event());
        unsigned int jet0eta = uint32_t(jet_raw.empty() ? 0 : jet_raw[0].eta()/0.01);
        std::uint32_t seed = jet0eta + 1 + (lumiNum_uint<<10) + (runNum_uint<<20) + evNum_uint;
        random_generator.seed(seed);
        
	    addJetCorrections(jet_raw, jecUnc, &resolution, &res_sf, *Rho, *genJets,
                          random_generator);

		delete jecUnc;
	}

	// QGLikelihood
	//addJetQGLikelihood(jet_raw, *qgHandle.product());
	addJetQGTaggerInfo(jet_raw, *qgHandle.product(), *axis2Handle.product(),
					   *ptDHandle.product(), *multHandle.product());
	
	// remove overlap
	//std::vector<pat::Jet> jet_no_lep = removeOverlapdR(jet_raw, lep_fakeable, 0.4);
	//std::vector<pat::Jet> jet_cleaned = removeOverlapdR(jet_no_lep, tau_loose, 0.4);

	// selected jets
	//std::vector<pat::Jet> jet_selected = getSelectedJets(jet_cleaned);
	std::vector<pat::Jet> jet_selected = getSelectedJets(jet_raw);

	if (debug_) {
		std::cout << "n_jets : " << jet_selected.size() << std::endl;
		dumpJets(jet_selected);
	}

	//////////////////
	// btags
	//std::vector<pat::Jet> jets_btag_loose;
	//std::vector<pat::Jet> jets_btag_medium;
	int n_btags_loose = 0;
	int n_btags_medium = 0;
	for (const auto& jet : jet_selected) {
		float csv = getJetCSV(jet);
		if (csv > csv_loose_wp_)
			n_btags_loose++;
		if (csv > csv_medium_wp_)
			n_btags_medium++;
	}

	if (debug_) {
		std::cout << "n_btags_loose : " << n_btags_loose << std::endl;
		std::cout << "n_btags_medium : " << n_btags_medium << std::endl;
	}

	/////////////////////////////////////////
	// MET
	if (debug_) std::cout << "mets" << std::endl;
	
	pat::MET pfMET = mets->front();
	float MET = pfMET.pt();
	// MHT
	float MHT = getMHT(lep_fakeable, tau_loose, jet_selected);
	float metLD = 0.00397 * MET + 0.00265 * MHT;

	/////////////////////////////////////////
	// Event selection
	/////////////////////////////////////////
    
	std::vector<miniLepton> *lep_selected = looseSelection_ ? &lep_loose : &lep_fakeable;

    // An inclusive selection that should cover all the region of interests
    // for event categories with hadronic taus
    bool pass_event_sel = evt_selector_ -> pass_ttH_ltau_inclusive_selection(
        *lep_selected, minitau_loose, jet_selected.size(), h_CutFlow_);

    if (not (pass_event_sel or event_selection_off_)) return;
    
	/////////////////////////////////////////
	// Write ntuple
	/////////////////////////////////////////

	//std::cout << "start writing ntuples" << std::endl;
	/// initialize ntuple
	evNtuple_.initialize();
	//std::cout << "done initialization" << std::endl;
	
	/// event id
	evNtuple_.run = iEvent.id().run();
	evNtuple_.ls = iEvent.id().luminosityBlock();
	evNtuple_.nEvent = iEvent.id().event();

	/// event variables
	if (not isdata_) {
		// pile up
		evNtuple_.npuTrue = npuTrue;
		evNtuple_.npuInTime = npuInTime;
		
		//
		if (sample_name_.Contains("ttH") or doSync_)
			evNtuple_.HiggsDecayType = HiggsDaughterPdgId(*MC_particles);

		//evNtuple_.nBadMuons = badMuons->size() + clonedMuons->size();
	}

	/// primary vertex
	reco::Vertex pv = getPrimaryVertex(vertices);
	evNtuple_.pvx = pv.x();
	evNtuple_.pvy = pv.y();
	evNtuple_.pvz = pv.z();
	
	evNtuple_.n_presel_mu = mu_preselected.size();
	evNtuple_.n_fakeable_mu = n_muon_fakeable;
	evNtuple_.n_mvasel_mu = n_muon_tight;
	evNtuple_.n_presel_ele = ele_preselected.size();
	evNtuple_.n_fakeable_ele = n_ele_fakeable;
	evNtuple_.n_mvasel_ele = n_ele_tight;
	evNtuple_.n_presel_tau = minitau_loose.size();
	evNtuple_.n_tau = minitau_tight.size();
	evNtuple_.n_jet = jet_selected.size();
	evNtuple_.n_btag_loose = n_btags_loose;
	evNtuple_.n_btag_medium = n_btags_medium;
	
	/// triggers and filters
	evNtuple_.triggerBits =
		trig_helper_->get_trigger_bits(triggerResults,hlt_config_);
	evNtuple_.filterBits =
		trig_helper_->get_filter_bits(filterResults,filter_config_);
	
	/// event weights
	if (not isdata_) {
		//evNtuple_.PU_weight = pu_weight;
		evNtuple_.MC_weight = mc_weight;
		evNtuple_.MC_weight_scale_muF0p5 = mc_weight_scale_muF0p5;
		evNtuple_.MC_weight_scale_muF2 = mc_weight_scale_muF2;
		evNtuple_.MC_weight_scale_muR0p5 = mc_weight_scale_muR0p5;
		evNtuple_.MC_weight_scale_muR2 = mc_weight_scale_muR2;
	}
	
	/// muons
	write_ntuple_muons(mu_preselected);
	
	/// electrons
	write_ntuple_electrons(ele_preselected);

	/// taus
	write_ntuple_taus(tau_loose);

	/// jets
	write_ntuple_jets(jet_selected);
	
	/// met
	write_ntuple_met(pfMET);
	evNtuple_.MHT = MHT;
	evNtuple_.metLD = metLD;

	////
	eventTree_ -> Fill();
	////
	
	firstpass_ = true;

	return;
}


// ------------ method called once each job just before starting event loop  ------------
void 
ttHtautauAnalyzer::beginJob()
{
	TH1::SetDefaultSumw2(true);
	TH2::SetDefaultSumw2(true);
	
	firstpass_ = false;

	event_count_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ttHtautauAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void ttHtautauAnalyzer::beginRun(const edm::Run &iRun, const edm::EventSetup &iSetup)
{
	// Update HLTConfigProvider for the new run
	bool hlt_config_changed = true; // init() updates this one
	
	if (!hlt_config_.init(iRun, iSetup, hltTag_, hlt_config_changed)) {
		std::cerr << "Warning, didn't find trigger process HLT,\t" << hltTag_
				  << std::endl;
		return;
	}

	if (hlt_config_changed) {
		std::cout << "New " << hltTag_ << " config has been loaded.\n";
		trig_helper_->add_trigger_version_number(hlt_config_);
	}

	//trig_helper_->dump_all_hlt_paths(hlt_config_);
	
	bool filter_config_changed = true; // init() updates this one
	
	if (!filter_config_.init(iRun, iSetup, filterTag_, filter_config_changed)) {
		std::cerr << "Warning, didn't find filter process HLT,\t" << filterTag_
				  << std::endl;
		return;
	}

	if (filter_config_changed) {
		std::cout << "New " << filterTag_ << " config has been loaded.\n";
	}

}

// ------------ method called when ending the processing of a run  ------------
void ttHtautauAnalyzer::endRun(const edm::Run &, const edm::EventSetup &)
{
	if (0/*dumpHLT_*/) {
		std::cout
			<< "***************************************************************"
			<< std::endl;
		std::cout << "  Summary for HLT: Total number of events = "
				  << event_count_ << std::endl;
		// FIXME
		std::cout
			<< "***************************************************************"
			<< std::endl;

		std::cout << "  Summary for Filters: Total number of events = "
				  << event_count_ << std::endl;

		std::cout
			<< "***************************************************************"
			<< std::endl;
	}

	std::cout
		<< "***************************************************************"
		<< std::endl;
	std::cout << "  Total number of events = " << event_count_ << std::endl;
	std::cout
		<< "***************************************************************"
		<< std::endl;
}

// ------------ method called when starting to processes a luminosity block
// ------------
/*
void ttHtautauAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock &,
const edm::EventSetup &)
{
}
*/

// ------------ method called when ending the processing of a luminosity block
// ------------
/*
void ttHtautauAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & lumi, const edm::EventSetup & iConfig)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ttHtautauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ttHtautauAnalyzer);
