#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "../interface/Types_enum.h"
#include "../interface/eventNtuple.h"
#include "../interface/syncNtuple.h"
#include "../interface/MVAVars.h"
#include "../interface/miniLepton.h"
#include "../interface/SFHelper.h"

#include <iostream>
#include <vector>
#include <algorithm>

TTree* makeSyncTree(const TString, const TString, const bool, const bool,
					Analysis_types, Selection_types);

void makeSyncNtuples(const TString dir="~/nobackup/ttHTT_syncNtuple/80X/Test/",
					 const TString output_file="syncNtuple_event.root")
{
	using namespace std;

	gROOT->ProcessLine(".L ../src/miniLepton.cc+");
	gROOT->ProcessLine(".L ../src/eventNtuple.cc+");
	gROOT->ProcessLine(".L ../src/syncNtuple.cc+");
	gROOT->ProcessLine(".L ../src/MVAVars.cc+");
	gROOT->ProcessLine(".L ../src/SFHelper.cc+");
	
	cout << "Opening root files from directory " << dir << endl;

	cout << "1l2tau signal region ... " << endl;
	TTree* synctree_1l2tau_sr =
		makeSyncTree(dir+"output_sync_event_1l2tau_sr.root",
					 "syncTree_1l2tau_SR", true, true,
					 Analysis_types::Analyze_1l2tau,
					 Selection_types::Signal_1l2tau);

	cout << "1l2tau fake extrapolation region ... " << endl;
	TTree* synctree_1l2tau_fake =
		makeSyncTree(dir+"output_sync_event_1l2tau_fake.root",
					 "syncTree_1l2tau_Fake", false, true,
					 Analysis_types::Analyze_1l2tau,
					 Selection_types::Control_fake_1l2tau);

	cout << "2lss1tau signal region ... " << endl;
	TTree* synctree_2lss1tau_sr =
		makeSyncTree(dir+"output_sync_event_2lss1tau_sr.root",
					 "syncTree_2lSS1tau_SR", true, true,
					 Analysis_types::Analyze_2lss1tau,
					 Selection_types::Signal_2lss1tau);

	cout << "2lss1tau fake extrapolation region ... " << endl;
	TTree* synctree_2lss1tau_fake =
		makeSyncTree(dir+"output_sync_event_2lss1tau_fake.root",
					 "syncTree_2lSS1tau_Fake", false, true,
					 Analysis_types::Analyze_2lss1tau,
					 Selection_types::Control_fake_2lss1tau);

	cout << "2lss1tau charge flip region ... " << endl;
	TTree* synctree_2lss1tau_flip =
		makeSyncTree(dir+"output_sync_event_2lss1tau_flip.root",
					 "syncTree_2lSS1tau_Flip", false, true,
					 Analysis_types::Analyze_2lss1tau,
					 Selection_types::Control_2los1tau);
	
	cout << "3l1tau signal region ... " << endl;
	TTree* synctree_3l1tau_sr =
		makeSyncTree(dir+"output_sync_event_3l1tau_sr.root",
					 "syncTree_3l1tau_SR", true, true,
					 Analysis_types::Analyze_3l1tau,
					 Selection_types::Signal_3l1tau);

	cout << "3l1tau fake extrapolation region ... " << endl;
	TTree* synctree_3l1tau_fake =
		makeSyncTree(dir+"output_sync_event_3l1tau_fake.root",
					 "syncTree_3l1tau_Fake", false, true,
					 Analysis_types::Analyze_3l1tau,
					 Selection_types::Control_fake_3l1tau);
	
	// event count
	cout << "1l2tau : " << endl;
	cout << "SR : " << synctree_1l2tau_sr->GetEntries() << endl;
	cout << "Fake : " << synctree_1l2tau_fake->GetEntries() << endl;
	cout << "2lSS1tau : " << endl;
	cout << "SR : " << synctree_2lss1tau_sr->GetEntries() << endl;
	cout << "Fake : " << synctree_2lss1tau_fake->GetEntries() << endl;
	cout << "Flip : " << synctree_2lss1tau_flip->GetEntries() << endl;
	cout << "3l1tau : " << endl;
	cout << "SR : " << synctree_3l1tau_sr->GetEntries() << endl;
	cout << "Fake : " << synctree_3l1tau_fake->GetEntries() << endl;
	
	// create output file
	TFile* fileout = new TFile(output_file, "recreate");
	
	synctree_2lss1tau_sr->Write();
	synctree_2lss1tau_fake->Write();
	synctree_2lss1tau_flip->Write();
	synctree_1l2tau_sr->Write();
	synctree_1l2tau_fake->Write();
	synctree_3l1tau_sr->Write();
	synctree_3l1tau_fake->Write();
	
	fileout -> Close();
}

TTree* makeSyncTree(const TString input_file, const TString treename,
					const bool isSignalRegion, const bool evtSel,
					Analysis_types anaType, Selection_types selType)
{
	using namespace std;

	// open file and read tree
	cout << "Opening input file : " << input_file << endl;
	TFile* f_in = new TFile(input_file,"read");
	TTree* tree_in = (TTree*)f_in->Get("ttHtaus/eventTree");

	eventNtuple evNtuple;
	evNtuple.set_branch_address(tree_in);

	// output tree
	TTree* tree_out = new TTree(treename, treename);
	tree_out -> SetDirectory(0);
	syncNtuple syncntuple;
	syncntuple.set_up_branches(tree_out);

	//////////////////////////////////////////////
	// Set up SFHelper
	SFHelper sf_helper(anaType, selType, false);
	
	//////////////////////////////////////////////
	// Set up TMVA Reader
	MVAVars mvaVars(anaType);
	mvaVars.set_up_tmva_reader();
	
	//////////////////////////////////////////////
	
	// event loop
	int nEntries = tree_in->GetEntries();
	for (int i = 0; i < nEntries; ++i) {
		tree_in->GetEntry(i);

		// apply additional event selection
		if (evtSel) {
			if (anaType==Analyze_2lss1tau) {
				// MC match
				//if (isSignalRegion and (not evNtuple.isGenMatchedLep)) continue;
				// tau charge
				if (not evNtuple.passTauCharge) continue;
			}
			else if (anaType==Analyze_1l2tau) {
				// MC match
				if (isSignalRegion and
					(not (evNtuple.isGenMatchedLep and evNtuple.isGenMatchedTau))
					) continue;
				// tau pair charge
				if (not evNtuple.passTauCharge) continue;
			}
			else if (anaType==Analyze_3l1tau) {
				// MC match
				//if (isSignalRegion and (not evNtuple.isGenMatchedLep)) continue;
				// charge sum
				if (not evNtuple.passTauCharge) continue;
			}
			else
				cerr << "Analysis type not supported" << endl;
		}

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
			syncntuple.mu0_pt = evNtuple.mu_pt->at(0);
			syncntuple.mu0_conept = evNtuple.mu_conept->at(0);
			syncntuple.mu0_eta = evNtuple.mu_eta->at(0);
			syncntuple.mu0_phi = evNtuple.mu_phi->at(0);
			syncntuple.mu0_E = evNtuple.mu_E->at(0);
			syncntuple.mu0_charge = evNtuple.mu_charge->at(0);
			syncntuple.mu0_dxy = evNtuple.mu_dxy->at(0);
			syncntuple.mu0_dz = evNtuple.mu_dz->at(0);
			syncntuple.mu0_isfakeablesel = evNtuple.mu_isfakeablesel->at(0);
			syncntuple.mu0_ismvasel = evNtuple.mu_ismvasel->at(0);
			syncntuple.mu0_jetNDauChargedMVASel = evNtuple.mu_jetNDauChargedMVASel->at(0);
			syncntuple.mu0_miniRelIso = evNtuple.mu_miniRelIso->at(0);
			syncntuple.mu0_miniIsoCharged = evNtuple.mu_miniIsoCharged->at(0);
			syncntuple.mu0_miniIsoNeutral = evNtuple.mu_miniIsoNeutral->at(0);
			syncntuple.mu0_jetPtRel = evNtuple.mu_jetPtRel->at(0);
			syncntuple.mu0_jetPtRatio = evNtuple.mu_jetPtRatio->at(0);
			syncntuple.mu0_jetCSV = evNtuple.mu_jetCSV->at(0);
			syncntuple.mu0_sip3D = evNtuple.mu_sip3D->at(0);
			syncntuple.mu0_segmentCompatibility = evNtuple.mu_segmentCompatibility->at(0);
			syncntuple.mu0_leptonMVA = evNtuple.mu_leptonMVA->at(0);
			syncntuple.mu0_mediumID = evNtuple.mu_mediumID->at(0);
			syncntuple.mu0_dpt_div_pt = evNtuple.mu_dpt_div_pt->at(0);
		}
		if (evNtuple.mu_pt->size()>1) {
			syncntuple.mu1_pt = evNtuple.mu_pt->at(1);
			syncntuple.mu1_conept = evNtuple.mu_conept->at(1);
			syncntuple.mu1_eta = evNtuple.mu_eta->at(1);
			syncntuple.mu1_phi = evNtuple.mu_phi->at(1);
			syncntuple.mu1_E = evNtuple.mu_E->at(1);
			syncntuple.mu1_charge = evNtuple.mu_charge->at(1);
			syncntuple.mu1_dxy = evNtuple.mu_dxy->at(1);
			syncntuple.mu1_dz = evNtuple.mu_dz->at(1);
			syncntuple.mu1_isfakeablesel = evNtuple.mu_isfakeablesel->at(1);
			syncntuple.mu1_ismvasel = evNtuple.mu_ismvasel->at(1);
			syncntuple.mu1_jetNDauChargedMVASel = evNtuple.mu_jetNDauChargedMVASel->at(1);
			syncntuple.mu1_miniRelIso = evNtuple.mu_miniRelIso->at(1);
			syncntuple.mu1_miniIsoCharged = evNtuple.mu_miniIsoCharged->at(1);
			syncntuple.mu1_miniIsoNeutral = evNtuple.mu_miniIsoNeutral->at(1);
			syncntuple.mu1_jetPtRel = evNtuple.mu_jetPtRel->at(1);
			syncntuple.mu1_jetPtRatio = evNtuple.mu_jetPtRatio->at(1);
			syncntuple.mu1_jetCSV = evNtuple.mu_jetCSV->at(1);
			syncntuple.mu1_sip3D = evNtuple.mu_sip3D->at(1);
			syncntuple.mu1_segmentCompatibility = evNtuple.mu_segmentCompatibility->at(1);
			syncntuple.mu1_leptonMVA = evNtuple.mu_leptonMVA->at(1);
			syncntuple.mu1_mediumID = evNtuple.mu_mediumID->at(1);
			syncntuple.mu1_dpt_div_pt = evNtuple.mu_dpt_div_pt->at(1);
		}

		// electrons
		if (evNtuple.ele_pt->size()>0) {
			syncntuple.ele0_pt = evNtuple.ele_pt->at(0);
			syncntuple.ele0_conept = evNtuple.ele_conept->at(0);
			syncntuple.ele0_eta = evNtuple.ele_eta->at(0);
			syncntuple.ele0_phi = evNtuple.ele_phi->at(0);
			syncntuple.ele0_E = evNtuple.ele_E->at(0);
			syncntuple.ele0_charge = evNtuple.ele_charge->at(0);
			syncntuple.ele0_dxy = evNtuple.ele_dxy->at(0);
			syncntuple.ele0_dz = evNtuple.ele_dz->at(0);
			syncntuple.ele0_isfakeablesel = evNtuple.ele_isfakeablesel->at(0);
			syncntuple.ele0_ismvasel = evNtuple.ele_ismvasel->at(0);
			syncntuple.ele0_jetNDauChargedMVASel = evNtuple.ele_jetNDauChargedMVASel->at(0);
			syncntuple.ele0_miniRelIso = evNtuple.ele_miniRelIso->at(0);
			syncntuple.ele0_miniIsoCharged = evNtuple.ele_miniIsoCharged->at(0);
			syncntuple.ele0_miniIsoNeutral = evNtuple.ele_miniIsoNeutral->at(0);
			syncntuple.ele0_jetPtRel = evNtuple.ele_jetPtRel->at(0);
			syncntuple.ele0_jetPtRatio = evNtuple.ele_jetPtRatio->at(0);
			syncntuple.ele0_jetCSV = evNtuple.ele_jetCSV->at(0);
			syncntuple.ele0_sip3D = evNtuple.ele_sip3D->at(0);
			syncntuple.ele0_ntMVAeleID = evNtuple.ele_ntMVAeleID->at(0);
			syncntuple.ele0_leptonMVA = evNtuple.ele_leptonMVA->at(0);
			syncntuple.ele0_isChargeConsistent = evNtuple.ele_isChargeConsistent->at(0);
			syncntuple.ele0_passesConversionVeto = evNtuple.ele_passesConversionVeto->at(0);
			syncntuple.ele0_nMissingHits = evNtuple.ele_nMissingHits->at(0);
		}
		if (evNtuple.ele_pt->size()>1) {
			syncntuple.ele1_pt = evNtuple.ele_pt->at(1);
			syncntuple.ele1_conept = evNtuple.ele_conept->at(1);
			syncntuple.ele1_eta = evNtuple.ele_eta->at(1);
			syncntuple.ele1_phi = evNtuple.ele_phi->at(1);
			syncntuple.ele1_E = evNtuple.ele_E->at(1);
			syncntuple.ele1_charge = evNtuple.ele_charge->at(1);
			syncntuple.ele1_dxy = evNtuple.ele_dxy->at(1);
			syncntuple.ele1_dz = evNtuple.ele_dz->at(1);
			syncntuple.ele1_isfakeablesel = evNtuple.ele_isfakeablesel->at(1);
			syncntuple.ele1_ismvasel = evNtuple.ele_ismvasel->at(1);
			syncntuple.ele1_jetNDauChargedMVASel = evNtuple.ele_jetNDauChargedMVASel->at(1);
			syncntuple.ele1_miniRelIso = evNtuple.ele_miniRelIso->at(1);
			syncntuple.ele1_miniIsoCharged = evNtuple.ele_miniIsoCharged->at(1);
			syncntuple.ele1_miniIsoNeutral = evNtuple.ele_miniIsoNeutral->at(1);
			syncntuple.ele1_jetPtRel = evNtuple.ele_jetPtRel->at(1);
			syncntuple.ele1_jetPtRatio = evNtuple.ele_jetPtRatio->at(1);
			syncntuple.ele1_jetCSV = evNtuple.ele_jetCSV->at(1);
			syncntuple.ele1_sip3D = evNtuple.ele_sip3D->at(1);
			syncntuple.ele1_ntMVAeleID = evNtuple.ele_ntMVAeleID->at(1);
			syncntuple.ele1_leptonMVA = evNtuple.ele_leptonMVA->at(1);
			syncntuple.ele1_isChargeConsistent = evNtuple.ele_isChargeConsistent->at(1);
			syncntuple.ele1_passesConversionVeto = evNtuple.ele_passesConversionVeto->at(1);
			syncntuple.ele1_nMissingHits = evNtuple.ele_nMissingHits->at(1);
		}

		// taus
		if (evNtuple.tau_pt->size()>0) {
			syncntuple.tau0_pt = evNtuple.tau_pt->at(0);
			syncntuple.tau0_eta = evNtuple.tau_eta->at(0);
			syncntuple.tau0_phi = evNtuple.tau_phi->at(0);
			syncntuple.tau0_E = evNtuple.tau_E->at(0);
			syncntuple.tau0_charge = evNtuple.tau_charge->at(0);
			syncntuple.tau0_dxy = evNtuple.tau_dxy->at(0);
			syncntuple.tau0_dz = evNtuple.tau_dz->at(0);
			syncntuple.tau0_decayModeFindingOldDMs = evNtuple.tau_decayModeFindingOldDMs->at(0);
			syncntuple.tau0_decayModeFindingNewDMs = evNtuple.tau_decayModeFindingNewDMs->at(0);
			syncntuple.tau0_byCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byCombinedIsolationDeltaBetaCorr3Hits->at(0);
			syncntuple.tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(0);
			syncntuple.tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(0);
			syncntuple.tau0_byTightCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(0);
			syncntuple.tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			syncntuple.tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			syncntuple.tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			syncntuple.tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			syncntuple.tau0_againstMuonLoose3 = evNtuple.tau_againstMuonLoose3->at(0);
			syncntuple.tau0_againstMuonTight3 = evNtuple.tau_againstMuonTight3->at(0);
			syncntuple.tau0_againstElectronVLooseMVA6 = evNtuple.tau_againstElectronVLooseMVA6->at(0);
			syncntuple.tau0_againstElectronLooseMVA6 = evNtuple.tau_againstElectronLooseMVA6->at(0);
			syncntuple.tau0_againstElectronMediumMVA6 = evNtuple.tau_againstElectronMediumMVA6->at(0);
			syncntuple.tau0_againstElectronTightMVA6 = evNtuple.tau_againstElectronTightMVA6->at(0);
		}
		if (evNtuple.tau_pt->size()>1) {
			syncntuple.tau1_pt = evNtuple.tau_pt->at(1);
			syncntuple.tau1_eta = evNtuple.tau_eta->at(1);
			syncntuple.tau1_phi = evNtuple.tau_phi->at(1);
			syncntuple.tau1_E = evNtuple.tau_E->at(1);
			syncntuple.tau1_charge = evNtuple.tau_charge->at(1);
			syncntuple.tau1_dxy = evNtuple.tau_dxy->at(1);
			syncntuple.tau1_dz = evNtuple.tau_dz->at(1);
			syncntuple.tau1_decayModeFindingOldDMs = evNtuple.tau_decayModeFindingOldDMs->at(1);
			syncntuple.tau1_decayModeFindingNewDMs = evNtuple.tau_decayModeFindingNewDMs->at(1);
			syncntuple.tau1_byCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byCombinedIsolationDeltaBetaCorr3Hits->at(1);
			syncntuple.tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(1);
			syncntuple.tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(1);
			syncntuple.tau1_byTightCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(1);
			syncntuple.tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			syncntuple.tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			syncntuple.tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			syncntuple.tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			syncntuple.tau1_againstMuonLoose3 = evNtuple.tau_againstMuonLoose3->at(1);
			syncntuple.tau1_againstMuonTight3 = evNtuple.tau_againstMuonTight3->at(1);
			syncntuple.tau1_againstElectronVLooseMVA6 = evNtuple.tau_againstElectronVLooseMVA6->at(1);
			syncntuple.tau1_againstElectronLooseMVA6 = evNtuple.tau_againstElectronLooseMVA6->at(1);
			syncntuple.tau1_againstElectronMediumMVA6 = evNtuple.tau_againstElectronMediumMVA6->at(1);
			syncntuple.tau1_againstElectronTightMVA6 = evNtuple.tau_againstElectronTightMVA6->at(1);
		}

		// jets
		if (evNtuple.jet_pt->size()>0) {
			syncntuple.jet0_pt = evNtuple.jet_pt->at(0);
			syncntuple.jet0_eta = evNtuple.jet_eta->at(0);
			syncntuple.jet0_phi = evNtuple.jet_phi->at(0);
			syncntuple.jet0_E = evNtuple.jet_E->at(0);
			syncntuple.jet0_CSV = evNtuple.jet_csv->at(0);
		}
		if (evNtuple.jet_pt->size()>1) {
			syncntuple.jet1_pt = evNtuple.jet_pt->at(1);
			syncntuple.jet1_eta = evNtuple.jet_eta->at(1);
			syncntuple.jet1_phi = evNtuple.jet_phi->at(1);
			syncntuple.jet1_E = evNtuple.jet_E->at(1);
			syncntuple.jet1_CSV = evNtuple.jet_csv->at(1);
		}
		if (evNtuple.jet_pt->size()>2) {
			syncntuple.jet2_pt = evNtuple.jet_pt->at(2);
			syncntuple.jet2_eta = evNtuple.jet_eta->at(2);
			syncntuple.jet2_phi = evNtuple.jet_phi->at(2);
			syncntuple.jet2_E = evNtuple.jet_E->at(2);
			syncntuple.jet2_CSV = evNtuple.jet_csv->at(2);
		}
		if (evNtuple.jet_pt->size()>3) {
			syncntuple.jet3_pt = evNtuple.jet_pt->at(3);
			syncntuple.jet3_eta = evNtuple.jet_eta->at(3);
			syncntuple.jet3_phi = evNtuple.jet_phi->at(3);
			syncntuple.jet3_E = evNtuple.jet_E->at(3);
			syncntuple.jet3_CSV = evNtuple.jet_csv->at(3);
		}

		syncntuple.ntags = evNtuple.n_btag_medium;
		syncntuple.ntags_loose = evNtuple.n_btag_loose;
		
		// met
		syncntuple.PFMET = evNtuple.PFMET;
		syncntuple.PFMETphi = evNtuple.PFMETphi;
		//
		//syncntuple.MHT = evNtuple.MHT;
		//syncntuple.metLD = evNtuple.metLD;
		auto looseleptons = evNtuple.buildLeptons(true);
	    auto loosetaus = evNtuple.buildFourVectorTaus(true);
	    auto ljets = evNtuple.buildFourVectorJets(true);

		TLorentzVector mht;
		for (auto l : looseleptons)
			mht -= l.p4();
		for (auto t : loosetaus)
			mht -= t;
		for (auto j : ljets)
			mht -= j;
		syncntuple.MHT = mht.Pt();
		syncntuple.metLD = 0.00397 * syncntuple.PFMET + 0.00265 * syncntuple.MHT;

		// event-level MVA variables
		syncntuple.isGenMatched =
			evNtuple.isGenMatchedLep * evNtuple.isGenMatchedTau;

		if (not evtSel) {
			tree_out->Fill();
			continue;
		}
		
		// build object four momentum
		// loose muons and electrons, preselected taus, selected jets are stored
		// in the ntuple
		// by default, select leptons at least passing fakeable selection
		// and tight tau (fakeable tau for control_fake_1l2tau)
	    auto leptons = evNtuple.buildLeptons();
	    auto taus = evNtuple.buildFourVectorTaus(selType==Control_fake_1l2tau);
	    auto jets = evNtuple.buildFourVectorJets();
		
	    mvaVars.compute_all_variables(leptons, taus, jets, evNtuple.PFMET,
									  evNtuple.PFMETphi, evNtuple.MHT,
									  evNtuple.n_btag_loose);

		if (anaType==Analyze_2lss1tau) {
			syncntuple.lep0_conept = mvaVars.lep0_conept();
			syncntuple.lep1_conept = mvaVars.lep1_conept();
			syncntuple.mindr_lep0_jet = mvaVars.mindr_lep0_jet();
			syncntuple.mindr_lep1_jet = mvaVars.mindr_lep1_jet();
			syncntuple.MT_met_lep0 = mvaVars.mT_met_lep0();
			syncntuple.avg_dr_jet = mvaVars.avg_dr_jet();
			syncntuple.dr_leps = mvaVars.dr_leps();
			syncntuple.mvis_lep0_tau = mvaVars.mvis_lep0_tau();
			syncntuple.mvis_lep1_tau = mvaVars.mvis_lep1_tau();
			syncntuple.max_lep_eta = mvaVars.max_lep_eta();
			syncntuple.dr_lep0_tau = mvaVars.dr_lep0_tau();
			
			syncntuple.MVA_2lss_ttV = mvaVars.BDT_ttV();
			syncntuple.MVA_2lss_ttbar = mvaVars.BDT_ttbar();
			syncntuple.MVA_2lSS1tau_noMEM_ttV = mvaVars.BDT_ttV();
			syncntuple.MVA_2lSS1tau_noMEM_ttbar = mvaVars.BDT_ttbar();
		}
		else if (anaType==Analyze_1l2tau) {
			
			syncntuple.tt_deltaR = mvaVars.dr_taus();
			syncntuple.tt_mvis = mvaVars.mvis_taus();
			syncntuple.tt_pt = mvaVars.pt_taus();
			syncntuple.max_dr_jet = mvaVars.max_dr_jet();
			syncntuple.HT = mvaVars.ht();
			syncntuple.MVA_1l2tau_ttbar_v2 = mvaVars.BDT_ttbar();
			syncntuple.MVA_1l2tau_ttZ_v2 = mvaVars.BDT_ttV();
		}
		else if (anaType==Analyze_3l1tau) {
			syncntuple.lep0_conept = mvaVars.lep0_conept();
			syncntuple.lep1_conept = mvaVars.lep1_conept();
			syncntuple.mindr_lep0_jet = mvaVars.mindr_lep0_jet();
			syncntuple.mindr_lep1_jet = mvaVars.mindr_lep1_jet();
			syncntuple.max_lep_eta = mvaVars.max_lep_eta();
			syncntuple.MT_met_lep0 = mvaVars.mT_met_lep0();
			//syncntuple.mvis_l1tau;
			//syncntuple.dR_l0tau;
			//syncntuple.dR_l1tau;
			//syncntuple.dR_l2tau;
			//syncntuple.MT_met_lep2;				
			syncntuple.MVA_3l1tau_ttbar = mvaVars.BDT_ttbar();
			syncntuple.MVA_3l1tau_ttV = mvaVars.BDT_ttV();
		}
		
		// weights
		syncntuple.PU_weight = evNtuple.PU_weight;
		syncntuple.MC_weight = evNtuple.MC_weight;
		syncntuple.bTagSF_weight = evNtuple.bTagSF_weight;
		syncntuple.FR_weight = evNtuple.FR_weight;
		// update weights
		//syncntuple.leptonSF_weight = evNtuple.leptonSF_weight;
		syncntuple.leptonSF_weight =
			sf_helper.Get_LeptonIDSF(leptons[0]) *
			sf_helper.Get_LeptonIDSF(leptons[1]);
		/*
		cout << "event# " << evNtuple.nEvent << endl;
		cout << "leptonSF_weight : " << syncntuple.leptonSF_weight << endl;
		cout << "lep0 (loose_to_reco tight_to_loose) : ";
		cout << sf_helper.Get_LeptonIDSF(leptons[0]) << " (";
		cout << sf_helper.Get_LeptonSF_loose(leptons[0]) << "  ";
		cout << sf_helper.Get_LeptonSF_tight_vs_loose(leptons[0]) <<")" << endl;
		cout << "lep1 (loose_to_reco tight_to_loose) : ";
		cout << sf_helper.Get_LeptonIDSF(leptons[1]) << " (";
		cout << sf_helper.Get_LeptonSF_loose(leptons[1]) << "  ";
		cout << sf_helper.Get_LeptonSF_tight_vs_loose(leptons[1]) <<")" << endl;
		*/
		//syncntuple.tauSF_weight = evNtuple.tauSF_weight;
		syncntuple.tauSF_weight =
			sf_helper.Get_TauIDSF(taus[0].Pt(), taus[0].Eta(),
								  evNtuple.isGenMatchedTau);
		syncntuple.triggerSF_weight = evNtuple.triggerSF_weight;
		//syncntuple.triggerSF_weight	= sf_helper.Get_HLTSF(evNtuple.lepCategory);
		
		tree_out->Fill();
		
	} // end of tree loop

	f_in->Close();
	
	return tree_out;
	
}
