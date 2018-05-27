// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>

#include "boost/program_options.hpp"
#include "../interface/Types_enum.h"
#include "../interface/eventNtuple.h"
#include "../interface/syncNtuple.h"
#include "../interface/miniLepton.h"
#include "../interface/miniTau.h"
#include "../interface/SFHelper.h"
#include "../interface/mvaNtuple.h"
#include "../interface/EventSelector.h"
#include "../interface/TriggerHelper.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <tuple>

TTree* makeSyncTree(const TString, const TString, const TString,
					Analysis_types anatype=Analyze_NA,
					Selection_types seltype=Selection_NA, bool evaluateMVA=false,
					bool doHTT=false, bool debug=false);


int main(int argc, char** argv)
{
	using namespace std;
	namespace po = boost::program_options;

	string dir, outname, infile, treename;
	bool makeObjectNtuple, make1l2tau, make2lss1tau, make3l1tau, make2l2tau;
	bool makeControl;
	bool evaluateMVA, doHTT;
	bool debug;
	
	po::options_description desc("Options");
	desc.add_options()
		("help,h", "produce help message")
		("directory,d", po::value<string>(&dir), "event ntuple directory")
		("output,o", po::value<string>(&outname), "output name")
		("input,i", po::value<string>(&infile)->default_value("output_sync_event_incl.root"))
		("treename,t", po::value<string>(&treename)->default_value("ttHtaus/eventTree"))
		("makeObjectNtuple", po::value<bool>(&makeObjectNtuple)->default_value(false))
		("make1l2tau", po::value<bool>(&make1l2tau)->default_value(false))
		("make2lss1tau", po::value<bool>(&make2lss1tau)->default_value(false))
		("make3l1tau", po::value<bool>(&make3l1tau)->default_value(false))
		("make2l2tau", po::value<bool>(&make2l2tau)->default_value(false))
		("makeControl", po::value<bool>(&makeControl)->default_value(false))
		("evaluateMVA", po::value<bool>(&evaluateMVA)->default_value(true))
		("doHTT", po::value<bool>(&doHTT)->default_value(true))
		("debug", po::value<bool>(&debug)->default_value(false));

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc << endl;
		return 1;
	}

	TString cdir = dir.c_str();

	TTree *synctree_obj = 0;
	TTree *synctree_1l2tau_sr = 0;
	TTree *synctree_1l2tau_fake = 0;
	TTree *synctree_2lss1tau_sr = 0;
	TTree *synctree_2lss1tau_fake = 0;
	TTree *synctree_2lss1tau_flip = 0;
	TTree *synctree_3l1tau_sr = 0;
	TTree *synctree_3l1tau_fake = 0;
	TTree *synctree_2l2tau_sr = 0;
	TTree *synctree_2l2tau_fake = 0;
	TTree *synctree_ttWctrl_sr = 0;
	TTree *synctree_ttWctrl_fake = 0;
	TTree *synctree_ttWctrl_flip = 0;
	TTree *synctree_ttZctrl_sr = 0;
	TTree *synctree_ttZctrl_fake = 0;

	if (makeObjectNtuple) {
		cout << "Object ntuple ... " << endl;
		synctree_obj = makeSyncTree(cdir+"output_sync.root","syncTree", treename);
		
		// event count
		cout << "number of events with at least 1 preselected muons : " << "\t"
			 << synctree_obj->GetEntries("n_presel_mu>0") << endl;
		cout << "number of events with at least 1 preselected electrons : " << "\t"
			 << synctree_obj->GetEntries("n_presel_ele>0") << endl;
		cout << "number of events with at least 1 preselected taus : " << "\t"
			 << synctree_obj->GetEntries("n_presel_tau>0") << endl;
		cout << "number of events with at least 1 preselected jets : " << "\t"
			 << synctree_obj->GetEntries("n_presel_jet>0") << endl;
	}

	if (make1l2tau) {
		cout << "1l2tau signal region ... " << endl;
		synctree_1l2tau_sr = makeSyncTree(cdir+infile.c_str(), "syncTree_1l2tau_SR",
										  treename, Analyze_1l2tau, Signal_1l2tau,
										  evaluateMVA, doHTT, debug);
		
		cout << "1l2tau fake extrapolation region ... " << endl;
		synctree_1l2tau_fake = makeSyncTree(cdir+infile.c_str(),
											"syncTree_1l2tau_Fake", treename, 
											Analyze_1l2tau, Application_Fake_1l2tau,
											evaluateMVA, doHTT, debug);
		
		// event count
		cout << "1l2tau : " << endl;
		cout << "SR : " << synctree_1l2tau_sr->GetEntries() << endl;
		cout << "Fake : " << synctree_1l2tau_fake->GetEntries() << endl;
	}

	if (make2lss1tau) {
		cout << "2lss1tau signal region ... " << endl;
		synctree_2lss1tau_sr =
			makeSyncTree(cdir+infile.c_str(), "syncTree_2lSS1tau_SR", treename, 
						 Analyze_2lss1tau, Signal_2lss1tau,
						 evaluateMVA, doHTT, debug);
		
		cout << "2lss1tau fake extrapolation region ... " << endl;
		synctree_2lss1tau_fake =
			makeSyncTree(cdir+infile.c_str(), "syncTree_2lSS1tau_Fake", treename, 
						 Analyze_2lss1tau, Application_Fake_2lss1tau,
						 evaluateMVA, doHTT, debug);

		cout << "2lss1tau charge flip region ... " << endl;
		synctree_2lss1tau_flip =
			makeSyncTree(cdir+infile.c_str(), "syncTree_2lSS1tau_Flip", treename, 
						 Analyze_2lss1tau, Application_Flip_2lss1tau,
						 evaluateMVA, doHTT, debug);
		
		// event count
		cout << "2lSS1tau : " << endl;
		cout << "SR : " << synctree_2lss1tau_sr->GetEntries() << endl;
		cout << "Fake : " << synctree_2lss1tau_fake->GetEntries() << endl;
		cout << "Flip : " << synctree_2lss1tau_flip->GetEntries() << endl;
	}

	if (make2l2tau) {
		cout << "2l2tau signal region ... " << endl;
		synctree_2l2tau_sr = makeSyncTree(cdir+infile.c_str(), "syncTree_2l2tau_SR",
										  treename, Analyze_2l2tau, Signal_2l2tau,
										  evaluateMVA, doHTT, debug);

		cout << "2l2tau fake extrapolation region ... " << endl;
		synctree_2l2tau_fake =
			makeSyncTree(cdir+infile.c_str(), "syncTree_2l2tau_Fake",
						 treename, Analyze_2l2tau, Application_Fake_2l2tau,
						 evaluateMVA, doHTT, debug);
		// event count
		cout << "2l2tau : " << endl;
		cout << "SR : " << synctree_2l2tau_sr->GetEntries() << endl;
		cout << "Fake : " << synctree_2l2tau_fake->GetEntries() << endl;
	}
	
	if (make3l1tau) {
		cout << "3l1tau signal region ... " << endl;
		synctree_3l1tau_sr = makeSyncTree(cdir+infile.c_str(), "syncTree_3l1tau_SR",
										  treename, Analyze_3l1tau, Signal_3l1tau,
										  evaluateMVA, doHTT, debug);

		cout << "3l1tau fake extrapolation region ... " << endl;
		synctree_3l1tau_fake =
			makeSyncTree(cdir+infile.c_str(), "syncTree_3l1tau_Fake",
						 treename, Analyze_3l1tau, Application_Fake_3l1tau,
						 evaluateMVA, doHTT, debug);
		// event count
		cout << "3l1tau : " << endl;
		cout << "SR : " << synctree_3l1tau_sr->GetEntries() << endl;
		cout << "Fake : " << synctree_3l1tau_fake->GetEntries() << endl;
	}

	if (makeControl) {
		cout << "ttW control region ... " << endl;
		synctree_ttWctrl_sr =
			makeSyncTree(cdir+infile.c_str(), "syncTree_ttWctrl_SR",
						 treename, Analyze_2lss, Control_ttW,
						 evaluateMVA, doHTT, debug);
		synctree_ttWctrl_fake =
			makeSyncTree(cdir+infile.c_str(), "syncTree_ttWctrl_Fake",
						 treename, Analyze_2lss, Control_FakeAR_ttW,
						 evaluateMVA, doHTT, debug);
		synctree_ttWctrl_flip =
			makeSyncTree(cdir+infile.c_str(), "syncTree_ttWctrl_Flip",
						 treename, Analyze_2lss, Control_FlipAR_ttW,
						 evaluateMVA, doHTT, debug);	
		cout << "ttZ control region ..." << endl;
		synctree_ttZctrl_sr =
			makeSyncTree(cdir+infile.c_str(), "syncTree_ttZctrl_SR",
						 treename, Analyze_3l, Control_ttZ,
						 evaluateMVA, doHTT, debug);
		synctree_ttZctrl_fake =
			makeSyncTree(cdir+infile.c_str(), "syncTree_ttZctrl_Fake",
						 treename, Analyze_3l, Control_FakeAR_ttZ,
						 evaluateMVA, doHTT, debug);
			
		// event count
		cout << "Control region : " << endl;
		cout << "ttW : " << synctree_ttWctrl_sr->GetEntries() << endl;
		cout << "ttW fakeAR : " << synctree_ttWctrl_fake->GetEntries() << endl;
		cout << "ttW flipAR : " << synctree_ttWctrl_flip->GetEntries() << endl;
		cout << "ttZ : " << synctree_ttZctrl_sr->GetEntries() << endl;
		cout << "ttZ fakeAR : " << synctree_ttZctrl_fake->GetEntries() << endl;
	}
	
	// create output file
	TFile* fileout = new TFile(outname.c_str(), "recreate");

	if (makeObjectNtuple) {
		synctree_obj->Write();
	}
	if (make1l2tau) {
		synctree_1l2tau_sr->Write();
		synctree_1l2tau_fake->Write();
	}
	if (make2lss1tau) {
		synctree_2lss1tau_sr->Write();
		synctree_2lss1tau_fake->Write();
		synctree_2lss1tau_flip->Write();
	}
	if (make3l1tau) {
		synctree_3l1tau_sr->Write();
		synctree_3l1tau_fake->Write();
	}
	if (make2l2tau) {
		synctree_2l2tau_sr->Write();
		synctree_2l2tau_fake->Write();
	}

	if (makeControl) {
		synctree_ttWctrl_sr->Write();
		synctree_ttWctrl_fake->Write();
		synctree_ttWctrl_flip->Write();
		synctree_ttZctrl_sr->Write();
		synctree_ttZctrl_fake->Write();
	}
	
	fileout->Close();

	cout << "output : " << outname.c_str() << endl;
	
	return 0;
}


TTree* makeSyncTree(const TString input_file, const TString treename,
					const TString intreename,
					Analysis_types anatype, Selection_types seltype,
					bool evaluateMVA, bool doHTT, bool debug)
{
	using namespace std;

	// open file and read tree
	cout << "Opening input file : " << input_file << endl;
	TFile* f_in = TFile::Open(input_file);
	TTree* tree_in = (TTree*)f_in->Get(intreename);

	eventNtuple evNtuple;
	evNtuple.set_branch_address(tree_in);

	// output tree
	TTree* tree_out = new TTree(treename, treename);
	tree_out -> SetDirectory(0);
	syncNtuple syncntuple;
	syncntuple.set_up_branches(tree_out);

	//////////////////////////////////////////////
	// Set up SFHelper
	SFHelper sf_helper(anatype, seltype, "ttHJetToNonbb_M125_amcatnlo", false, debug);

	//////////////////////////////////////////////
	// trigger Helper
	TriggerHelper trig_helper(Analysis_types::Analyze_inclusive, false);
	//TriggerHelper trig_helper(anatype, false);

	//////////////////////////////////////////////
	// mva ntuple
	mvaNtuple mvantuple(anatype, false, "2017", doHTT, evaluateMVA);

	//////////////////////////////////////////////
	// event selector
	EventSelector evt_selector(debug, true);
	
	//////////////////////////////////////////////
	
	// event loop
	int nEntries = tree_in->GetEntries();
	for (int i = 0; i < nEntries; ++i) {
		tree_in -> GetEntry(i);

		// reconstruct objects
		auto leptons_loose = evNtuple.buildLeptons('L');  // loose
	    auto leptons = evNtuple.buildLeptons('F');  // fakeable
		auto leptons_tight = evNtuple.buildLeptons('T'); // tight

		auto taus_fakeable = evNtuple.buildTaus(true, anatype); // fakeable
		auto taus_tight = evNtuple.buildTaus(false, anatype);  // tight

		//auto jets = evNtuple.buildJets();
		// Jet cleaning based on analysis type
		auto jets = evNtuple.buildCleanedJets(0.4, anatype, seltype, &leptons,
											  &taus_fakeable);

		int nbtags_loose, nbtags_medium;
		std::tie(nbtags_loose, nbtags_medium) = evNtuple.count_btags(jets);
		assert(nbtags_loose >= nbtags_medium);

		auto metp4 = evNtuple.buildFourVectorMET();
		float met = metp4.Pt();
		float mht = evNtuple.computeMHT(leptons, taus_fakeable, jets);
		float metld = 0.00397 * met + 0.00265 * mht;
		
		if (debug) {
			std::cout << std::endl;
			std::cout << "Event: " << evNtuple.run <<":"<< evNtuple.ls << ":"
					  << evNtuple.nEvent << std::endl;
		}
		
		if (anatype != Analyze_NA) {

			////////////////////////////////////////
			bool passEvtSel = evt_selector.pass_full_event_selection(
			    anatype, seltype, leptons_loose, leptons, leptons_tight,
				taus_fakeable, taus_tight, jets.size(), nbtags_loose,
				nbtags_medium, metld);
			
			if (not passEvtSel) continue;

			////////////////////////////////////////
			// HLT and MET filters
			int nfakeableEle = 0, nfkaeableMu = 0;
			for (const auto & lep : leptons) {
				if (abs(lep.pdgId())==11) nfakeableEle++;
				if (abs(lep.pdgId())==13) nfkaeableMu++;
			}

			bool isdata = false;
			if (not evt_selector.pass_hlt_and_filters(anatype, &trig_helper,
													  evNtuple.triggerBits,
													  nfakeableEle, nfkaeableMu,
													  evNtuple.filterBits, isdata))
				continue;
		}
		
		////////////////////////////////////////
		std::vector<miniTau> * taus = &taus_tight;
		//if (seltype==Application_Fake_1l2tau or seltype==Application_Fake_2l2tau)
		if (anatype==Analyze_1l2tau or anatype==Analyze_2l2tau)
			taus = &taus_fakeable;
		
		syncntuple.initialize();

		syncntuple.run = evNtuple.run;
		syncntuple.ls = evNtuple.ls;
		syncntuple.nEvent = evNtuple.nEvent;

		syncntuple.n_presel_mu = evNtuple.n_presel_mu;
		syncntuple.n_mvasel_mu = evNtuple.n_mvasel_mu;
		syncntuple.n_fakeablesel_mu = evNtuple.n_fakeable_mu;
		syncntuple.n_presel_ele = evNtuple.n_presel_ele;
		syncntuple.n_mvasel_ele = evNtuple.n_mvasel_ele;
		syncntuple.n_fakeablesel_ele = evNtuple.n_fakeable_ele;
		syncntuple.n_presel_tau = evNtuple.n_presel_tau;
		syncntuple.n_presel_jet = jets.size();

		// muons
		if (evNtuple.mu_pt->size()>0) {
			syncntuple.mu1_pt = evNtuple.mu_pt->at(0);
			syncntuple.mu1_conept = evNtuple.mu_conept->at(0);
			syncntuple.mu1_eta = evNtuple.mu_eta->at(0);
			syncntuple.mu1_phi = evNtuple.mu_phi->at(0);
			syncntuple.mu1_E = evNtuple.mu_E->at(0);
			syncntuple.mu1_charge = evNtuple.mu_charge->at(0);
			syncntuple.mu1_dxy = evNtuple.mu_dxy->at(0);
			syncntuple.mu1_dxyAbs = fabs(evNtuple.mu_dxy->at(0));
			syncntuple.mu1_dz = evNtuple.mu_dz->at(0);
			syncntuple.mu1_isfakeablesel = evNtuple.mu_isfakeablesel->at(0);
			syncntuple.mu1_ismvasel = evNtuple.mu_ismvasel->at(0);
			syncntuple.mu1_jetNDauChargedMVASel = evNtuple.mu_jetNDauChargedMVASel->at(0);
			syncntuple.mu1_miniRelIso = evNtuple.mu_miniRelIso->at(0);
			syncntuple.mu1_miniIsoCharged = evNtuple.mu_miniIsoCharged->at(0);
			syncntuple.mu1_miniIsoNeutral = evNtuple.mu_miniIsoNeutral->at(0);
			syncntuple.mu1_jetPtRel = evNtuple.mu_jetPtRel->at(0);
			syncntuple.mu1_jetPtRatio = evNtuple.mu_jetPtRatio->at(0);
			syncntuple.mu1_jetCSV = evNtuple.mu_jetCSV->at(0);
			syncntuple.mu1_sip3D = evNtuple.mu_sip3D->at(0);
			syncntuple.mu1_segmentCompatibility = evNtuple.mu_segmentCompatibility->at(0);
			syncntuple.mu1_leptonMVA = evNtuple.mu_leptonMVA->at(0);
			syncntuple.mu1_PFRelIso04 = evNtuple.mu_pfRelIso04->at(0);
			syncntuple.mu1_mediumID = evNtuple.mu_mediumID->at(0);
			syncntuple.mu1_dpt_div_pt = evNtuple.mu_dpt_div_pt->at(0);
		}
		if (evNtuple.mu_pt->size()>1) {
			syncntuple.mu2_pt = evNtuple.mu_pt->at(1);
			syncntuple.mu2_conept = evNtuple.mu_conept->at(1);
			syncntuple.mu2_eta = evNtuple.mu_eta->at(1);
			syncntuple.mu2_phi = evNtuple.mu_phi->at(1);
			syncntuple.mu2_E = evNtuple.mu_E->at(1);
			syncntuple.mu2_charge = evNtuple.mu_charge->at(1);
			syncntuple.mu2_dxy = evNtuple.mu_dxy->at(1);
			syncntuple.mu2_dxyAbs = fabs(evNtuple.mu_dxy->at(1));
			syncntuple.mu2_dz = evNtuple.mu_dz->at(1);
			syncntuple.mu2_isfakeablesel = evNtuple.mu_isfakeablesel->at(1);
			syncntuple.mu2_ismvasel = evNtuple.mu_ismvasel->at(1);
			syncntuple.mu2_jetNDauChargedMVASel = evNtuple.mu_jetNDauChargedMVASel->at(1);
			syncntuple.mu2_miniRelIso = evNtuple.mu_miniRelIso->at(1);
			syncntuple.mu2_miniIsoCharged = evNtuple.mu_miniIsoCharged->at(1);
			syncntuple.mu2_miniIsoNeutral = evNtuple.mu_miniIsoNeutral->at(1);
			syncntuple.mu2_jetPtRel = evNtuple.mu_jetPtRel->at(1);
			syncntuple.mu2_jetPtRatio = evNtuple.mu_jetPtRatio->at(1);
			syncntuple.mu2_jetCSV = evNtuple.mu_jetCSV->at(1);
			syncntuple.mu2_sip3D = evNtuple.mu_sip3D->at(1);
			syncntuple.mu2_segmentCompatibility = evNtuple.mu_segmentCompatibility->at(1);
			syncntuple.mu2_leptonMVA = evNtuple.mu_leptonMVA->at(1);
			syncntuple.mu2_PFRelIso04 = evNtuple.mu_pfRelIso04->at(1);
			syncntuple.mu2_mediumID = evNtuple.mu_mediumID->at(1);
			syncntuple.mu2_dpt_div_pt = evNtuple.mu_dpt_div_pt->at(1);
		}
		
		// electrons
		if (evNtuple.ele_pt->size()>0) {
			syncntuple.ele1_pt = evNtuple.ele_pt->at(0);
			syncntuple.ele1_conept = evNtuple.ele_conept->at(0);
			syncntuple.ele1_eta = evNtuple.ele_eta->at(0);
			syncntuple.ele1_phi = evNtuple.ele_phi->at(0);
			syncntuple.ele1_E = evNtuple.ele_E->at(0);
			syncntuple.ele1_charge = evNtuple.ele_charge->at(0);
			syncntuple.ele1_dxy = evNtuple.ele_dxy->at(0);
			syncntuple.ele1_dxyAbs = fabs(evNtuple.ele_dxy->at(0));
			syncntuple.ele1_dz = evNtuple.ele_dz->at(0);
			syncntuple.ele1_isfakeablesel = evNtuple.ele_isfakeablesel->at(0);
			syncntuple.ele1_ismvasel = evNtuple.ele_ismvasel->at(0);
			syncntuple.ele1_jetNDauChargedMVASel = evNtuple.ele_jetNDauChargedMVASel->at(0);
			syncntuple.ele1_miniRelIso = evNtuple.ele_miniRelIso->at(0);
			syncntuple.ele1_miniIsoCharged = evNtuple.ele_miniIsoCharged->at(0);
			syncntuple.ele1_miniIsoNeutral = evNtuple.ele_miniIsoNeutral->at(0);
			syncntuple.ele1_jetPtRel = evNtuple.ele_jetPtRel->at(0);
			syncntuple.ele1_jetPtRatio = evNtuple.ele_jetPtRatio->at(0);
			syncntuple.ele1_jetCSV = evNtuple.ele_jetCSV->at(0);
			syncntuple.ele1_sip3D = evNtuple.ele_sip3D->at(0);
			syncntuple.ele1_ntMVAeleID = evNtuple.ele_ntMVAeleID->at(0);
			syncntuple.ele1_leptonMVA = evNtuple.ele_leptonMVA->at(0);
			syncntuple.ele1_PFRelIso04 = evNtuple.ele_pfRelIso04->at(0);
			syncntuple.ele1_isChargeConsistent = evNtuple.ele_isChargeConsistent->at(0);
			syncntuple.ele1_passesConversionVeto = evNtuple.ele_passesConversionVeto->at(0);
			syncntuple.ele1_nMissingHits = evNtuple.ele_nMissingHits->at(0);
		}
		if (evNtuple.ele_pt->size()>1) {
			syncntuple.ele2_pt = evNtuple.ele_pt->at(1);
			syncntuple.ele2_conept = evNtuple.ele_conept->at(1);
			syncntuple.ele2_eta = evNtuple.ele_eta->at(1);
			syncntuple.ele2_phi = evNtuple.ele_phi->at(1);
			syncntuple.ele2_E = evNtuple.ele_E->at(1);
			syncntuple.ele2_charge = evNtuple.ele_charge->at(1);
			syncntuple.ele2_dxy = evNtuple.ele_dxy->at(1);
			syncntuple.ele2_dxyAbs = fabs(evNtuple.ele_dxy->at(1));
			syncntuple.ele2_dz = evNtuple.ele_dz->at(1);
			syncntuple.ele2_isfakeablesel = evNtuple.ele_isfakeablesel->at(1);
			syncntuple.ele2_ismvasel = evNtuple.ele_ismvasel->at(1);
			syncntuple.ele2_jetNDauChargedMVASel = evNtuple.ele_jetNDauChargedMVASel->at(1);
			syncntuple.ele2_miniRelIso = evNtuple.ele_miniRelIso->at(1);
			syncntuple.ele2_miniIsoCharged = evNtuple.ele_miniIsoCharged->at(1);
			syncntuple.ele2_miniIsoNeutral = evNtuple.ele_miniIsoNeutral->at(1);
			syncntuple.ele2_jetPtRel = evNtuple.ele_jetPtRel->at(1);
			syncntuple.ele2_jetPtRatio = evNtuple.ele_jetPtRatio->at(1);
			syncntuple.ele2_jetCSV = evNtuple.ele_jetCSV->at(1);
			syncntuple.ele2_sip3D = evNtuple.ele_sip3D->at(1);
			syncntuple.ele2_ntMVAeleID = evNtuple.ele_ntMVAeleID->at(1);
			syncntuple.ele2_leptonMVA = evNtuple.ele_leptonMVA->at(1);
			syncntuple.ele2_PFRelIso04 = evNtuple.ele_pfRelIso04->at(1);
			syncntuple.ele2_isChargeConsistent = evNtuple.ele_isChargeConsistent->at(1);
			syncntuple.ele2_passesConversionVeto = evNtuple.ele_passesConversionVeto->at(1);
			syncntuple.ele2_nMissingHits = evNtuple.ele_nMissingHits->at(1);
		}

		// taus
		if (evNtuple.tau_pt->size()>0) {
			syncntuple.tau1_pt = evNtuple.tau_pt->at(0);
			syncntuple.tau1_eta = evNtuple.tau_eta->at(0);
			syncntuple.tau1_phi = evNtuple.tau_phi->at(0);
			syncntuple.tau1_E = evNtuple.tau_E->at(0);
			syncntuple.tau1_charge = evNtuple.tau_charge->at(0);
			syncntuple.tau1_dxy = evNtuple.tau_dxy->at(0);
			syncntuple.tau1_dz = evNtuple.tau_dz->at(0);
			syncntuple.tau1_decayModeFindingOldDMs = evNtuple.tau_decayModeFindingOldDMs->at(0);
			syncntuple.tau1_decayModeFindingNewDMs = evNtuple.tau_decayModeFindingNewDMs->at(0);
			syncntuple.tau1_byCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byCombinedIsolationDeltaBetaCorr3Hits->at(0);
			syncntuple.tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(0);
			syncntuple.tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(0);
			syncntuple.tau1_byTightCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(0);
			syncntuple.tau1_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			syncntuple.tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			syncntuple.tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			syncntuple.tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			syncntuple.tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			syncntuple.tau1_rawMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byIsolationMVArun2v1DBdR03oldDMwLTraw->at(0);
			syncntuple.tau1_againstMuonLoose3 = evNtuple.tau_againstMuonLoose3->at(0);
			syncntuple.tau1_againstMuonTight3 = evNtuple.tau_againstMuonTight3->at(0);
			syncntuple.tau1_againstElectronVLooseMVA6 = evNtuple.tau_againstElectronVLooseMVA6->at(0);
			syncntuple.tau1_againstElectronLooseMVA6 = evNtuple.tau_againstElectronLooseMVA6->at(0);
			syncntuple.tau1_againstElectronMediumMVA6 = evNtuple.tau_againstElectronMediumMVA6->at(0);
			syncntuple.tau1_againstElectronTightMVA6 = evNtuple.tau_againstElectronTightMVA6->at(0);
		}
		if (evNtuple.tau_pt->size()>1) {
			syncntuple.tau2_pt = evNtuple.tau_pt->at(1);
			syncntuple.tau2_eta = evNtuple.tau_eta->at(1);
			syncntuple.tau2_phi = evNtuple.tau_phi->at(1);
			syncntuple.tau2_E = evNtuple.tau_E->at(1);
			syncntuple.tau2_charge = evNtuple.tau_charge->at(1);
			syncntuple.tau2_dxy = evNtuple.tau_dxy->at(1);
			syncntuple.tau2_dz = evNtuple.tau_dz->at(1);
			syncntuple.tau2_decayModeFindingOldDMs = evNtuple.tau_decayModeFindingOldDMs->at(1);
			syncntuple.tau2_decayModeFindingNewDMs = evNtuple.tau_decayModeFindingNewDMs->at(1);
			syncntuple.tau2_byCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byCombinedIsolationDeltaBetaCorr3Hits->at(1);
			syncntuple.tau2_byLooseCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(1);
			syncntuple.tau2_byMediumCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(1);
			syncntuple.tau2_byTightCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(1);
			syncntuple.tau2_byVLooseIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			syncntuple.tau2_byLooseIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			syncntuple.tau2_byMediumIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			syncntuple.tau2_byTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			syncntuple.tau2_byVTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			syncntuple.tau2_rawMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byIsolationMVArun2v1DBdR03oldDMwLTraw->at(1);
			syncntuple.tau2_againstMuonLoose3 = evNtuple.tau_againstMuonLoose3->at(1);
			syncntuple.tau2_againstMuonTight3 = evNtuple.tau_againstMuonTight3->at(1);
			syncntuple.tau2_againstElectronVLooseMVA6 = evNtuple.tau_againstElectronVLooseMVA6->at(1);
			syncntuple.tau2_againstElectronLooseMVA6 = evNtuple.tau_againstElectronLooseMVA6->at(1);
			syncntuple.tau2_againstElectronMediumMVA6 = evNtuple.tau_againstElectronMediumMVA6->at(1);
			syncntuple.tau2_againstElectronTightMVA6 = evNtuple.tau_againstElectronTightMVA6->at(1);
		}

		// jets
		/*
		if (evNtuple.jet_pt->size()>0) {
			syncntuple.jet1_pt = evNtuple.jet_pt->at(0);
			syncntuple.jet1_eta = evNtuple.jet_eta->at(0);
			syncntuple.jet1_phi = evNtuple.jet_phi->at(0);
			syncntuple.jet1_E = evNtuple.jet_E->at(0);
			syncntuple.jet1_CSV = evNtuple.jet_csv->at(0);
		}
		if (evNtuple.jet_pt->size()>1) {
			syncntuple.jet2_pt = evNtuple.jet_pt->at(1);
			syncntuple.jet2_eta = evNtuple.jet_eta->at(1);
			syncntuple.jet2_phi = evNtuple.jet_phi->at(1);
			syncntuple.jet2_E = evNtuple.jet_E->at(1);
			syncntuple.jet2_CSV = evNtuple.jet_csv->at(1);
		}
		if (evNtuple.jet_pt->size()>2) {
			syncntuple.jet3_pt = evNtuple.jet_pt->at(2);
			syncntuple.jet3_eta = evNtuple.jet_eta->at(2);
			syncntuple.jet3_phi = evNtuple.jet_phi->at(2);
			syncntuple.jet3_E = evNtuple.jet_E->at(2);
			syncntuple.jet3_CSV = evNtuple.jet_csv->at(2);
		}
		if (evNtuple.jet_pt->size()>3) {
			syncntuple.jet4_pt = evNtuple.jet_pt->at(3);
			syncntuple.jet4_eta = evNtuple.jet_eta->at(3);
			syncntuple.jet4_phi = evNtuple.jet_phi->at(3);
			syncntuple.jet4_E = evNtuple.jet_E->at(3);
			syncntuple.jet4_CSV = evNtuple.jet_csv->at(3);
		}
		*/
		if (jets.size()>0) {
			syncntuple.jet1_pt = jets[0].pt();
			syncntuple.jet1_eta = jets[0].eta();
			syncntuple.jet1_phi = jets[0].phi();
			syncntuple.jet1_E = jets[0].energy();
			syncntuple.jet1_CSV = jets[0].csv();
		}
		if (jets.size()>1) {
			syncntuple.jet2_pt = jets[1].pt();
			syncntuple.jet2_eta = jets[1].eta();
			syncntuple.jet2_phi = jets[1].phi();
			syncntuple.jet2_E = jets[1].energy();
			syncntuple.jet2_CSV = jets[1].csv();
		}
		if (jets.size()>2) {
			syncntuple.jet3_pt = jets[2].pt();
			syncntuple.jet3_eta = jets[2].eta();
			syncntuple.jet3_phi = jets[2].phi();
			syncntuple.jet3_E = jets[2].energy();
			syncntuple.jet3_CSV = jets[2].csv();
		}
		if (jets.size()>3) {
			syncntuple.jet4_pt = jets[3].pt();
			syncntuple.jet4_eta = jets[3].eta();
			syncntuple.jet4_phi = jets[3].phi();
			syncntuple.jet4_E = jets[3].energy();
			syncntuple.jet4_CSV = jets[3].csv();
		}

		
		// met
		syncntuple.PFMET = metp4.Pt();
		syncntuple.PFMETphi = metp4.Phi();
		syncntuple.MHT = mht;
		syncntuple.metLD = metld;

		if (anatype==Analyze_NA or seltype==Selection_NA or
			seltype==Control_ttW or seltype==Control_ttZ) {
			tree_out->Fill();
			continue;
		}

		// mva variables

		//assert(leptons.size()>0 and taus->size()>0);

		mvantuple.compute_all_variables(leptons, *taus, jets, syncntuple.PFMET,
										syncntuple.PFMETphi, syncntuple.MHT,
										nbtags_loose, nbtags_medium);
		
		syncntuple.nBJetLoose = nbtags_loose;
		
		if (anatype==Analyze_1l2tau) {
			assert(taus->size()>1);
			syncntuple.isGenMatched = leptons[0].isGenMatched() and
				taus->at(0).isGenMatched() and taus->at(1).isGenMatched();
			syncntuple.avg_dr_jet = mvantuple.avg_dr_jet;
			syncntuple.dr_taus = mvantuple.dr_taus;
			//mvantuple.met;
			syncntuple.lep1_conept = mvantuple.lep0_conept;
			syncntuple.mT_lep1 = mvantuple.mT_met_lep0;
			syncntuple.mTauTauVis = mvantuple.mTauTauVis;
			syncntuple.mindr_lep1_jet = mvantuple.mindr_lep0_jet;
			syncntuple.mindr_tau1_jet = mvantuple.mindr_tau0_jet;
			syncntuple.mindr_tau2_jet = mvantuple.mindr_tau1_jet;
			syncntuple.dr_lep1_tau = mvantuple.dr_lep_tau_lead;
			//mvantuple.nbtags_loose;
			//mvantuple.tau0_pt;
			//mvantuple.tau1_pt;
			syncntuple.dR_lep_tau_ss = mvantuple.dr_lep_tau_ss;
			syncntuple.cosThetaS_hadTau = mvantuple.costS_tau;
			if (doHTT) {
				syncntuple.HTT = mvantuple.HTT;
				syncntuple.HadTop_pt = mvantuple.HadTop_pt;
			}
			if (evaluateMVA) {
				syncntuple.mvaOutput_plainKin_ttbar =  mvantuple.mva_1l2tau_BDT1;
				syncntuple.mvaOutput_1l_2tau_HTT_SUM_VT = mvantuple.mva_1l2tau_BDT2;
			}
		}

		if (anatype==Analyze_2lss1tau) {
			assert(leptons.size()>1);
			syncntuple.isGenMatched = leptons[0].isGenMatched() and
				leptons[1].isGenMatched() and taus->at(0).isGenMatched();
			syncntuple.avg_dr_jet = mvantuple.avg_dr_jet;
			syncntuple.dr_lep1_tau = mvantuple.dr_lep0_tau;
			syncntuple.dr_lep2_tau = mvantuple.dr_lep1_tau;
			syncntuple.dr_leps = mvantuple.dr_leps;
			syncntuple.lep1_conept = mvantuple.lep0_conept;
			syncntuple.lep2_conept = mvantuple.lep1_conept;
			syncntuple.mT_lep1 = mvantuple.mT_met_lep0;
			syncntuple.mT_lep2 = mvantuple.mT_met_lep1;
			syncntuple.mTauTauVis1 = mvantuple.mvis_lep0_tau;
			syncntuple.mTauTauVis2 = mvantuple.mvis_lep1_tau;
			syncntuple.max_lep_eta = mvantuple.max_lep_eta;
			syncntuple.mbb = mvantuple.mbb;
			syncntuple.mindr_lep1_jet = mvantuple.mindr_lep0_jet;
			syncntuple.mindr_lep2_jet = mvantuple.mindr_lep1_jet;
			syncntuple.mindr_tau1_jet = mvantuple.mindr_tau0_jet;
			//mvantuple.nJet;
			//mvantuple.met;
			//mvantuple.tau0_pt;

			if (doHTT) {
				syncntuple.HTT = mvantuple.HTT;
				syncntuple.HadTop_pt = mvantuple.HadTop_pt;
			}
			// Hj_tagger
			// memOutput_LR
			if (evaluateMVA) {
				syncntuple.mvaOutput_2lss_1tau_plainKin_ttbar =
					mvantuple.mva_2lss1tau_BDT2;
				syncntuple.mvaOutput_2lss_1tau_plainKin_ttV =
					mvantuple.mva_2lss1tau_BDT1;
				syncntuple.mvaOutput_2lss_1tau_plainKin_1B_M =
					mvantuple.mva_2lss1tau_BDT6;
				syncntuple.mvaOutput_2lss_1tau_plainKin_SUM_M =
					mvantuple.mva_2lss1tau_BDT3;
				syncntuple.mvaOutput_2lss_1tau_HTT_SUM_M =
					mvantuple.mva_2lss1tau_BDT4;
				syncntuple.mvaOutput_2lss_1tau_HTTMEM_SUM_M =
					mvantuple.mva_2lss1tau_BDT5;
			}
		}

		if (anatype==Analyze_3l1tau) {
			assert(leptons.size()>2);
			syncntuple.isGenMatched = leptons[0].isGenMatched() and
				leptons[1].isGenMatched() and leptons[2].isGenMatched() and
				taus->at(0).isGenMatched();
			syncntuple.mTauTauVis1 = mvantuple.mvis_lep0_tau;
			syncntuple.mTauTauVis2 = mvantuple.mvis_lep1_tau;
			syncntuple.max_lep_eta = mvantuple.max_lep_eta;
			//mvantuple.met;
			//mvantuple.tau0_pt;
			syncntuple.dr_leps = mvantuple.dr_leps;
			syncntuple.mindr_lep1_jet = mvantuple.mindr_lep0_jet;
			syncntuple.mindr_lep2_jet = mvantuple.mindr_lep1_jet;
			syncntuple.mindr_lep3_jet = mvantuple.mindr_lep2_jet;
			syncntuple.mT_lep1 = mvantuple.mT_met_lep0;
			syncntuple.mT_lep2 = mvantuple.mT_met_lep1;
			syncntuple.avg_dr_jet = mvantuple.avg_dr_jet;
			syncntuple.mindr_tau1_jet = mvantuple.mindr_tau0_jet;
			syncntuple.mbb_loose = mvantuple.mbb;
			syncntuple.lep1_conept = mvantuple.lep0_conept;
			syncntuple.lep2_conept = mvantuple.lep1_conept;
			syncntuple.lep3_conept = mvantuple.lep2_conept;
			//mvantuple.nJet;
			if (evaluateMVA) {
				syncntuple.mvaOutput_plainKin_ttV = mvantuple.mva_3l1tau_BDT1;
				syncntuple.mvaOutput_plainKin_ttbar = mvantuple.mva_3l1tau_BDT2;
				syncntuple.mvaOutput_3l_1tau_plainKin_SUM_M =
					mvantuple.mva_3l1tau_BDT3;
				syncntuple.mvaOutput_3l_1tau_plainKin_1B_M =
					mvantuple.mva_3l1tau_BDT4;
			}
		}

		if (anatype==Analyze_2l2tau) {
			assert(taus->size()>1 and leptons.size()>1);
			syncntuple.mTauTauVis = mvantuple.mTauTauVis;
			syncntuple.cosThetaS_hadTau = mvantuple.costS_tau;
			syncntuple.min_dr_lep_jet = mvantuple.min_dr_lep_jet;
			syncntuple.mindr_tau_jet = mvantuple.mindr_tau_jet;
			syncntuple.mT_lep1 = mvantuple.mT_met_lep0;
			syncntuple.mT_lep2 = mvantuple.mT_met_lep1;	
			syncntuple.mindr_lep1_jet = mvantuple.mindr_lep0_jet;
			syncntuple.max_dr_lep_tau = mvantuple.max_dr_lep_tau;
			syncntuple.nBJetLoose = mvantuple.nbtags_loose;
			syncntuple.dr_taus = mvantuple.dr_taus;
			syncntuple.avg_dr_jet = mvantuple.avg_dr_jet;
			syncntuple.lep1_conept = mvantuple.lep0_conept;
			syncntuple.lep2_conept = mvantuple.lep1_conept;
			syncntuple.min_dr_lep_tau = mvantuple.min_dr_lep_tau;
			syncntuple.mindr_tau1_jet = mvantuple.mindr_tau0_jet;
			syncntuple.avg_dr_lep_tau = mvantuple.avg_dr_lep_tau;
			if (evaluateMVA) {
				syncntuple.mvaOutput_plainKin_ttV = mvantuple.mva_2l2tau_BDT1;
				syncntuple.mvaOutput_plainKin_ttbar = mvantuple.mva_2l2tau_BDT2;
				syncntuple.mvaOutput_2l_2tau_plainKin_1B_VT =
					mvantuple.mva_2l2tau_BDT4;
				syncntuple.mvaOutput_2l_2tau_plainKin_SUM_VT =
					mvantuple.mva_2l2tau_BDT3;
			}
		}

		// weights
		syncntuple.MC_weight = evNtuple.MC_weight;

		syncntuple.PU_weight = sf_helper.Get_PUWeight(evNtuple.npuTrue);

		syncntuple.bTagSF_weight = sf_helper.Get_EvtCSVWeight(jets,"NA");

		// FR_weight
		if (seltype==Application_Fake_1l2tau or seltype==Application_Fake_2lss1tau or
			seltype==Application_Fake_3l1tau or seltype==Application_Fake_2l2tau or
			seltype==Control_FakeAR_ttW or seltype==Control_FakeAR_ttZ)
			syncntuple.FR_weight = sf_helper.Get_FR_weight(leptons,*taus);
		else if (seltype==Application_Flip_2lss1tau) {
			syncntuple.FR_weight =
				sf_helper.Get_ChargeFlipWeight(leptons, taus->at(0).charge());
		}
		else if (seltype==Control_FlipAR_ttW) {
			syncntuple.FR_weight = sf_helper.Get_ChargeFlipWeight(leptons);
		}
		
		// leptonSF_weight
		// UPDATE NEEDED: loose vs reco
		syncntuple.leptonSF_weight = sf_helper.Get_LeptonIDSF_weight(leptons);
		
		// tauSF_weight
		syncntuple.tauSF_weight = sf_helper.Get_TauIDSF_weight(*taus);
			
		// triggerSF
		bool hlt1LTriggered =
				trig_helper.pass_single_lep_triggers(evNtuple.triggerBits);
		bool hltXTriggered =
			trig_helper.pass_leptau_cross_triggers(evNtuple.triggerBits);
		syncntuple.triggerSF_weight =
			sf_helper.Get_HLTSF(leptons, *taus, hlt1LTriggered, hltXTriggered);
		
		tree_out->Fill();
		
	} // end of event loop

	f_in->Close();

	return tree_out;
}
