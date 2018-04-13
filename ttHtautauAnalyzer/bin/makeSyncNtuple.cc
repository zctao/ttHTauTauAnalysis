// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

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

TTree* makeSyncTree(const TString, const TString, Analysis_types anatype=Analyze_NA,
					Selection_types seltype=Selection_NA, bool debug=false,
					const TString intreename="ttHtaus/eventTree");


int main(int argc, char** argv)
{
	using namespace std;
	namespace po = boost::program_options;

	string dir, outname;
	bool makeObjectNtuple, make1l2tau, make2lss1tau, make3l1tau;
	bool debug;
	
	po::options_description desc("Options");
	desc.add_options()
		("help,h", "produce help message")
		("directory,d", po::value<string>(&dir), "event ntuple directory")
		("output,o", po::value<string>(&outname), "output name")
		("makeObjectNtuple", po::value<bool>(&makeObjectNtuple)->default_value(false))
		("make1l2tau", po::value<bool>(&make1l2tau)->default_value(false))
		("make2lss1tau", po::value<bool>(&make2lss1tau)->default_value(false))
		("make3l1tau", po::value<bool>(&make3l1tau)->default_value(false))
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

	if (makeObjectNtuple) {
		cout << "Object ntuple ... " << endl;
		synctree_obj = makeSyncTree(cdir+"output_sync.root","syncTree");
		
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
		synctree_1l2tau_sr = makeSyncTree(cdir+"output_sync_event_1l2tau_incl.root",
										  "syncTree_1l2tau_SR",
										  Analyze_1l2tau, Signal_1l2tau, debug);
		
		cout << "1l2tau fake extrapolation region ... " << endl;
		synctree_1l2tau_fake = makeSyncTree(cdir+"output_sync_event_1l2tau_incl.root",
											"syncTree_1l2tau_Fake",
											Analyze_1l2tau, Control_fake_1l2tau, debug);
		
		// event count
		cout << "1l2tau : " << endl;
		cout << "SR : " << synctree_1l2tau_sr->GetEntries() << endl;
		cout << "Fake : " << synctree_1l2tau_fake->GetEntries() << endl;
	}

	if (make2lss1tau) {
		cout << "2lss1tau signal region ... " << endl;
		synctree_2lss1tau_sr =
			makeSyncTree(cdir+"output_sync_event_2lss1tau_incl.root",
						 "syncTree_2lSS1tau_SR",
						 Analyze_2lss1tau, Signal_2lss1tau, debug);
		
		cout << "2lss1tau fake extrapolation region ... " << endl;
		synctree_2lss1tau_fake =
			makeSyncTree(cdir+"output_sync_event_2lss1tau_incl.root",
						 "syncTree_2lSS1tau_Fake",
						 Analyze_2lss1tau, Control_fake_2lss1tau, debug);

		cout << "2lss1tau charge flip region ... " << endl;
		synctree_2lss1tau_flip =
			makeSyncTree(cdir+"output_sync_event_2lss1tau_incl.root",
						 "syncTree_2lSS1tau_Flip",
						 Analyze_2lss1tau, Control_2los1tau, debug);
		
		// event count
		cout << "2lSS1tau : " << endl;
		cout << "SR : " << synctree_2lss1tau_sr->GetEntries() << endl;
		cout << "Fake : " << synctree_2lss1tau_fake->GetEntries() << endl;
		cout << "Flip : " << synctree_2lss1tau_flip->GetEntries() << endl;
	}

	if (make3l1tau) {
		cout << "3l1tau signal region ... " << endl;
		synctree_3l1tau_sr = makeSyncTree(cdir+"output_sync_event_3l1tau_incl.root",
										  "syncTree_3l1tau_SR",
										  Analyze_3l1tau, Signal_3l1tau, debug);

		cout << "3l1tau fake extrapolation region ... " << endl;
		synctree_3l1tau_fake = makeSyncTree(cdir+"output_sync_event_3l1tau_incl.root",
											"syncTree_3l1tau_Fake",
											Analyze_3l1tau, Control_fake_3l1tau, debug);

		// event count
		cout << "3l1tau : " << endl;
		cout << "SR : " << synctree_3l1tau_sr->GetEntries() << endl;
		cout << "Fake : " << synctree_3l1tau_fake->GetEntries() << endl;
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
	
	fileout->Close();

	cout << "output : " << outname.c_str() << endl;
	
	return 0;
}


TTree* makeSyncTree(const TString input_file, const TString treename,
					Analysis_types anatype, Selection_types seltype, bool debug,
					const TString intreename)
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
	SFHelper sf_helper(anatype, seltype, false, debug);

	//////////////////////////////////////////////
	// trigger Helper
	TriggerHelper trig_helper(anatype, false);

	//////////////////////////////////////////////
	// mva ntuple
	mvaNtuple mvantuple(anatype, false, "2016");

	//////////////////////////////////////////////
	// event selector
	EventSelector evt_selector(debug, true);
	
	//////////////////////////////////////////////
	
	// event loop
	int nEntries = tree_in->GetEntries();
	for (int i = 0; i < nEntries; ++i) {
		tree_in -> GetEntry(i);

		// reconstruct objects
	    auto leptons = evNtuple.buildLeptons();  // fakeable
		auto taus_tight = evNtuple.buildTaus(); // tight taus
		// fakeable taus for 1l2tau fake AR
	    auto taus_fakeable = evNtuple.buildTaus(true);

		if (debug) {
			std::cout << std::endl;
			std::cout << "Event: " << evNtuple.run <<":"<< evNtuple.ls << ":"
					  << evNtuple.nEvent << std::endl;
		}

		if (anatype != Analyze_NA) {
		
			////////////////////////////////////////
			if (not evt_selector.pass_extra_event_selection(anatype, seltype,
															&leptons, &taus_tight,
															&taus_fakeable))
			continue;

			////////////////////////////////////////
			// HLT path
			if (not evt_selector.pass_hlt_paths(anatype, &trig_helper,
												evNtuple.triggerBits))
				continue;
		}
		
		////////////////////////////////////////
		std::vector<miniTau> * taus = &taus_tight;
		if (seltype==Control_fake_1l2tau)
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
		syncntuple.n_presel_jet = evNtuple.n_jet;

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

		syncntuple.ntags = evNtuple.n_btag_medium;
		syncntuple.ntags_loose = evNtuple.n_btag_loose;
		
		// met
		syncntuple.PFMET = evNtuple.PFMET;
		syncntuple.PFMETphi = evNtuple.PFMETphi;
		//
		//syncntuple.MHT = evNtuple.MHT;
		//syncntuple.metLD = evNtuple.metLD;
		syncntuple.MHT = evNtuple.computeMHT();
		syncntuple.metLD = 0.00397 * syncntuple.PFMET + 0.00265 * syncntuple.MHT;

		if (anatype==Analyze_NA or seltype==Selection_NA) {
			tree_out->Fill();
			continue;
		}

		// mva variables
		//auto leptons = evNtuple.buildLeptons();
		//auto taus = evNtuple.buildTaus(seltype==Control_fake_1l2tau);
		auto jets = evNtuple.buildFourVectorJets();

		assert(leptons.size()>0);
		assert(leptons.size()>0);
		assert(jets.size()>0);
		
		syncntuple.lep1_conept = leptons[0].conept();
		syncntuple.mindr_lep1_jet = mvantuple.compute_min_dr(leptons[0].p4(),jets);
		syncntuple.mindr_tau_jet = mvantuple.compute_min_dr(taus->at(0).p4(),jets);
		syncntuple.MT_met_lep1 = mvantuple.compute_mT_lep(leptons[0], syncntuple.PFMET, syncntuple.PFMETphi);
		syncntuple.avg_dr_jet = mvantuple.compute_average_dr(jets);
		syncntuple.max_dr_jet = mvantuple.compute_max_dr(jets);
		syncntuple.HT = syncntuple.MHT;
		syncntuple.dR_l0tau = leptons[0].p4().DeltaR(taus->at(0).p4());
		syncntuple.mvis_l0tau = (leptons[0].p4()+taus->at(0).p4()).M();
		
		if (leptons.size()>1) {
			syncntuple.lep2_conept = leptons[1].conept();
			syncntuple.mindr_lep2_jet = mvantuple.compute_min_dr(leptons[1].p4(),jets);
			syncntuple.mvis_l1tau = (leptons[1].p4()+taus->at(0).p4()).M();
			syncntuple.dR_l1tau = leptons[1].p4().DeltaR(taus->at(0).p4());
			syncntuple.dR_leps = leptons[0].p4().DeltaR(leptons[1].p4());
		}

		if (leptons.size()>2) {
			syncntuple.mindr_lep3_jet = mvantuple.compute_min_dr(leptons[2].p4(),jets);
			syncntuple.dR_l2tau = leptons[2].p4().DeltaR(taus->at(0).p4());
			syncntuple.MT_met_lep3 = mvantuple.compute_mT_lep(leptons[2], syncntuple.PFMET, syncntuple.PFMETphi);;
		}

		if (taus->size()>1) {
			syncntuple.tt_deltaR = taus->at(0).p4().DeltaR(taus->at(1).p4());
			syncntuple.tt_mvis = (taus->at(0).p4()+taus->at(1).p4()).M();
			syncntuple.tt_pt = (taus->at(0).p4()+taus->at(1).p4()).Pt();
		}		

		// weights
		syncntuple.PU_weight = evNtuple.PU_weight;
		syncntuple.MC_weight = evNtuple.MC_weight;
		syncntuple.bTagSF_weight = evNtuple.bTagSF_weight;

		// FR_weight
		if (seltype==Control_fake_1l2tau or seltype==Control_fake_2lss1tau or
			seltype==Control_fake_3l1tau)
			syncntuple.FR_weight = sf_helper.Get_FR_weight(leptons,*taus);
		else if (seltype==Control_2los1tau)
			syncntuple.FR_weight =
				sf_helper.Get_ChargeFlipWeight(leptons, taus->at(0).charge());
		
		// leptonSF_weight
		syncntuple.leptonSF_weight = sf_helper.Get_LeptonIDSF_weight(leptons);
		
		// tauSF_weight
		syncntuple.tauSF_weight = sf_helper.Get_TauIDSF_weight(*taus);
		
		// triggerSF_weight
		bool hlt1LTriggered =
				trig_helper.pass_single_lep_triggers(evNtuple.triggerBits);
		bool hltXTriggered =
			trig_helper.pass_leptau_cross_triggers(evNtuple.triggerBits);
		syncntuple.triggerSF_weight =
			sf_helper.Get_HLTSF(leptons, *taus, hlt1LTriggered, hltXTriggered);

		//syncntuple.isGenMatched =
			//evNtuple.isGenMatchedLep * evNtuple.isGenMatchedTau;
		
		tree_out->Fill();
		
	} // end of event loop

	f_in->Close();

	return tree_out;
}
