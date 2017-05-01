#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "../interface/eventNtuple.h"
#include "../interface/syncNtuple.h"

#include <iostream>
#include <vector>

void makeFlatNtuple(const TString input_file="../test/output_.root", const TString output_file="test.root", const bool evtSel = true)
{
	using namespace std;

	gROOT->ProcessLine(".L ../src/eventNtuple.cc+");
	gROOT->ProcessLine(".L ../src/syncNtuple.cc+");
	
	// open file and read tree
	cout << "Opening input file : " << input_file << endl;
	TFile* f_in = new TFile(input_file,"read");
	TTree* tree_in = (TTree*)f_in->Get("ttHtaus/eventTree");

	eventNtuple evNtuple;
	evNtuple.set_branch_address(tree_in);

	// output file
	cout << "Output file created: " << output_file << endl;
	TFile* f_out = new TFile(output_file, "recreate");
	// output tree
	TTree* tree_flat = new TTree("eventTree", "Sync Tree");	
	syncNtuple flatNtuple;
	flatNtuple.set_up_branches(tree_flat);

	//
	int nEntries = tree_in->GetEntries();
	for (int i = 0; i < nEntries; ++i) {
		tree_in->GetEntry(i);

		// apply additional event selection
		if (evtSel) {
			if (not evNtuple.isGenMatchedLep) continue;
			if (not evNtuple.passTauCharge) continue;
		}
		
		flatNtuple.initialize();

		flatNtuple.run = evNtuple.run;
		flatNtuple.ls = evNtuple.ls;
		flatNtuple.nEvent = evNtuple.nEvent;

		flatNtuple.event_weight = evNtuple.event_weight;
		flatNtuple.PU_weight = evNtuple.PU_weight;
		flatNtuple.MC_weight = evNtuple.MC_weight;
		flatNtuple.bTagSF_weight = evNtuple.bTagSF_weight;
		flatNtuple.leptonSF_weight = evNtuple.leptonSF_weight;
		flatNtuple.tauSF_weight = evNtuple.tauSF_weight;
		flatNtuple.triggerSF_weight = evNtuple.triggerSF_weight;
		flatNtuple.FR_weight = evNtuple.FR_weight;

		/*
		flatNtuple.MC_weight_scale_muF0p5 = evNtuple.MC_weight_scale_muF0p5;
		flatNtuple.MC_weight_scale_muF2 = evNtuple.MC_weight_scale_muF2;
		flatNtuple.MC_weight_scale_muR0p5 = evNtuple.MC_weight_scale_muR0p5;
		flatNtuple.MC_weight_scale_muR2 = evNtuple.MC_weight_scale_muR2;
		flatNtuple.btagSF_weight_LFUp = evNtuple.btagSF_weight_LFUp;
		flatNtuple.btagSF_weight_LFDown = evNtuple.btagSF_weight_LFDown;
		flatNtuple.btagSF_weight_HFUp = evNtuple.btagSF_weight_HFUp;
		flatNtuple.btagSF_weight_HFDown = evNtuple.btagSF_weight_HFDown;
		flatNtuple.btagSF_weight_HFStats1Up = evNtuple.btagSF_weight_HFStats1Up;
		flatNtuple.btagSF_weight_HFStats1Down = evNtuple.btagSF_weight_HFStats1Down;
		flatNtuple.btagSF_weight_HFStats2Up = evNtuple.btagSF_weight_HFStats2Up;
		flatNtuple.btagSF_weight_HFStats2Down = evNtuple.btagSF_weight_HFStats2Down;
		flatNtuple.btagSF_weight_LFStats1Up = evNtuple.btagSF_weight_LFStats1Up;
		flatNtuple.btagSF_weight_LFStats1Down = evNtuple.btagSF_weight_LFStats1Down;
		flatNtuple.btagSF_weight_LFStats2Up = evNtuple.btagSF_weight_LFStats2Up;
		flatNtuple.btagSF_weight_LFStats2Down = evNtuple.btagSF_weight_LFStats2Down;
		flatNtuple.btagSF_weight_cErr1Up = evNtuple.btagSF_weight_cErr1Up;
		flatNtuple.btagSF_weight_cErr1Down = evNtuple.btagSF_weight_cErr1Down;
		flatNtuple.btagSF_weight_cErr2Up = evNtuple.btagSF_weight_cErr2Up;
		flatNtuple.btagSF_weight_cErr2Down = evNtuple.btagSF_weight_cErr2Down;
		*/
		
		flatNtuple.n_presel_mu = evNtuple.n_presel_mu;
		flatNtuple.n_mvasel_mu = evNtuple.n_mvasel_mu;
		flatNtuple.n_fakeablesel_mu = evNtuple.n_fakeable_mu;
		flatNtuple.n_presel_ele = evNtuple.n_presel_ele;
		flatNtuple.n_mvasel_ele = evNtuple.n_mvasel_ele;
		flatNtuple.n_fakeablesel_ele = evNtuple.n_fakeable_ele;
		flatNtuple.n_presel_tau = evNtuple.n_presel_tau;
		flatNtuple.n_presel_jet = evNtuple.n_jet;
		
		// muons
		if (evNtuple.mu_pt->size()>0) {
			flatNtuple.mu0_pt = evNtuple.mu_pt->at(0);
			flatNtuple.mu0_conept = evNtuple.mu_conept->at(0);
			flatNtuple.mu0_eta = evNtuple.mu_eta->at(0);
			flatNtuple.mu0_phi = evNtuple.mu_phi->at(0);
			flatNtuple.mu0_E = evNtuple.mu_E->at(0);
			flatNtuple.mu0_charge = evNtuple.mu_charge->at(0);
			flatNtuple.mu0_dxy = evNtuple.mu_dxy->at(0);
			flatNtuple.mu0_dz = evNtuple.mu_dz->at(0);
			flatNtuple.mu0_isfakeablesel = evNtuple.mu_isfakeablesel->at(0);
			flatNtuple.mu0_ismvasel = evNtuple.mu_ismvasel->at(0);
			flatNtuple.mu0_jetNDauChargedMVASel = evNtuple.mu_jetNDauChargedMVASel->at(0);
			flatNtuple.mu0_miniRelIso = evNtuple.mu_miniRelIso->at(0);
			flatNtuple.mu0_miniIsoCharged = evNtuple.mu_miniIsoCharged->at(0);
			flatNtuple.mu0_miniIsoNeutral = evNtuple.mu_miniIsoNeutral->at(0);
			flatNtuple.mu0_jetPtRel = evNtuple.mu_jetPtRel->at(0);
			flatNtuple.mu0_jetPtRatio = evNtuple.mu_jetPtRatio->at(0);
			flatNtuple.mu0_jetCSV = evNtuple.mu_jetCSV->at(0);
			flatNtuple.mu0_sip3D = evNtuple.mu_sip3D->at(0);
			flatNtuple.mu0_segmentCompatibility = evNtuple.mu_segmentCompatibility->at(0);
			flatNtuple.mu0_leptonMVA = evNtuple.mu_leptonMVA->at(0);
			flatNtuple.mu0_mediumID = evNtuple.mu_mediumID->at(0);
		}
		if (evNtuple.mu_pt->size()>1) {
			flatNtuple.mu1_pt = evNtuple.mu_pt->at(1);
			flatNtuple.mu1_conept = evNtuple.mu_conept->at(1);
			flatNtuple.mu1_eta = evNtuple.mu_eta->at(1);
			flatNtuple.mu1_phi = evNtuple.mu_phi->at(1);
			flatNtuple.mu1_E = evNtuple.mu_E->at(1);
			flatNtuple.mu1_charge = evNtuple.mu_charge->at(1);
			flatNtuple.mu1_dxy = evNtuple.mu_dxy->at(1);
			flatNtuple.mu1_dz = evNtuple.mu_dz->at(1);
			flatNtuple.mu1_isfakeablesel = evNtuple.mu_isfakeablesel->at(1);
			flatNtuple.mu1_ismvasel = evNtuple.mu_ismvasel->at(1);
			flatNtuple.mu1_jetNDauChargedMVASel = evNtuple.mu_jetNDauChargedMVASel->at(1);
			flatNtuple.mu1_miniRelIso = evNtuple.mu_miniRelIso->at(1);
			flatNtuple.mu1_miniIsoCharged = evNtuple.mu_miniIsoCharged->at(1);
			flatNtuple.mu1_miniIsoNeutral = evNtuple.mu_miniIsoNeutral->at(1);
			flatNtuple.mu1_jetPtRel = evNtuple.mu_jetPtRel->at(1);
			flatNtuple.mu1_jetPtRatio = evNtuple.mu_jetPtRatio->at(1);
			flatNtuple.mu1_jetCSV = evNtuple.mu_jetCSV->at(1);
			flatNtuple.mu1_sip3D = evNtuple.mu_sip3D->at(1);
			flatNtuple.mu1_segmentCompatibility = evNtuple.mu_segmentCompatibility->at(1);
			flatNtuple.mu1_leptonMVA = evNtuple.mu_leptonMVA->at(1);
			flatNtuple.mu1_mediumID = evNtuple.mu_mediumID->at(1);
		}

		// electrons
		if (evNtuple.ele_pt->size()>0) {
			flatNtuple.ele0_pt = evNtuple.ele_pt->at(0);
			flatNtuple.ele0_conept = evNtuple.ele_conept->at(0);
			flatNtuple.ele0_eta = evNtuple.ele_eta->at(0);
			flatNtuple.ele0_phi = evNtuple.ele_phi->at(0);
			flatNtuple.ele0_E = evNtuple.ele_E->at(0);
			flatNtuple.ele0_charge = evNtuple.ele_charge->at(0);
			flatNtuple.ele0_dxy = evNtuple.ele_dxy->at(0);
			flatNtuple.ele0_dz = evNtuple.ele_dz->at(0);
			flatNtuple.ele0_isfakeablesel = evNtuple.ele_isfakeablesel->at(0);
			flatNtuple.ele0_ismvasel = evNtuple.ele_ismvasel->at(0);
			flatNtuple.ele0_jetNDauChargedMVASel = evNtuple.ele_jetNDauChargedMVASel->at(0);
			flatNtuple.ele0_miniRelIso = evNtuple.ele_miniRelIso->at(0);
			flatNtuple.ele0_miniIsoCharged = evNtuple.ele_miniIsoCharged->at(0);
			flatNtuple.ele0_miniIsoNeutral = evNtuple.ele_miniIsoNeutral->at(0);
			flatNtuple.ele0_jetPtRel = evNtuple.ele_jetPtRel->at(0);
			flatNtuple.ele0_jetPtRatio = evNtuple.ele_jetPtRatio->at(0);
			flatNtuple.ele0_jetCSV = evNtuple.ele_jetCSV->at(0);
			flatNtuple.ele0_sip3D = evNtuple.ele_sip3D->at(0);
			flatNtuple.ele0_ntMVAeleID = evNtuple.ele_ntMVAeleID->at(0);
			flatNtuple.ele0_leptonMVA = evNtuple.ele_leptonMVA->at(0);
			flatNtuple.ele0_isChargeConsistent = evNtuple.ele_isChargeConsistent->at(0);
			flatNtuple.ele0_passesConversionVeto = evNtuple.ele_passesConversionVeto->at(0);
			flatNtuple.ele0_nMissingHits = evNtuple.ele_nMissingHits->at(0);
		}
		if (evNtuple.ele_pt->size()>1) {
			flatNtuple.ele1_pt = evNtuple.ele_pt->at(1);
			flatNtuple.ele1_conept = evNtuple.ele_conept->at(1);
			flatNtuple.ele1_eta = evNtuple.ele_eta->at(1);
			flatNtuple.ele1_phi = evNtuple.ele_phi->at(1);
			flatNtuple.ele1_E = evNtuple.ele_E->at(1);
			flatNtuple.ele1_charge = evNtuple.ele_charge->at(1);
			flatNtuple.ele1_dxy = evNtuple.ele_dxy->at(1);
			flatNtuple.ele1_dz = evNtuple.ele_dz->at(1);
			flatNtuple.ele1_isfakeablesel = evNtuple.ele_isfakeablesel->at(1);
			flatNtuple.ele1_ismvasel = evNtuple.ele_ismvasel->at(1);
			flatNtuple.ele1_jetNDauChargedMVASel = evNtuple.ele_jetNDauChargedMVASel->at(1);
			flatNtuple.ele1_miniRelIso = evNtuple.ele_miniRelIso->at(1);
			flatNtuple.ele1_miniIsoCharged = evNtuple.ele_miniIsoCharged->at(1);
			flatNtuple.ele1_miniIsoNeutral = evNtuple.ele_miniIsoNeutral->at(1);
			flatNtuple.ele1_jetPtRel = evNtuple.ele_jetPtRel->at(1);
			flatNtuple.ele1_jetPtRatio = evNtuple.ele_jetPtRatio->at(1);
			flatNtuple.ele1_jetCSV = evNtuple.ele_jetCSV->at(1);
			flatNtuple.ele1_sip3D = evNtuple.ele_sip3D->at(1);
			flatNtuple.ele1_ntMVAeleID = evNtuple.ele_ntMVAeleID->at(1);
			flatNtuple.ele1_leptonMVA = evNtuple.ele_leptonMVA->at(1);
			flatNtuple.ele1_isChargeConsistent = evNtuple.ele_isChargeConsistent->at(1);
			flatNtuple.ele1_passesConversionVeto = evNtuple.ele_passesConversionVeto->at(1);
			flatNtuple.ele1_nMissingHits = evNtuple.ele_nMissingHits->at(1);
		}

		// taus
		if (evNtuple.tau_pt->size()>0) {
			flatNtuple.tau0_pt = evNtuple.tau_pt->at(0);
			flatNtuple.tau0_eta = evNtuple.tau_eta->at(0);
			flatNtuple.tau0_phi = evNtuple.tau_phi->at(0);
			flatNtuple.tau0_E = evNtuple.tau_E->at(0);
			flatNtuple.tau0_charge = evNtuple.tau_charge->at(0);
			flatNtuple.tau0_dxy = evNtuple.tau_dxy->at(0);
			flatNtuple.tau0_dz = evNtuple.tau_dz->at(0);
			flatNtuple.tau0_decayModeFindingOldDMs = evNtuple.tau_decayModeFindingOldDMs->at(0);
			flatNtuple.tau0_decayModeFindingNewDMs = evNtuple.tau_decayModeFindingNewDMs->at(0);
			flatNtuple.tau0_byCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byCombinedIsolationDeltaBetaCorr3Hits->at(0);
			flatNtuple.tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(0);
			flatNtuple.tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(0);
			flatNtuple.tau0_byTightCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(0);
			flatNtuple.tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			flatNtuple.tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			flatNtuple.tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			flatNtuple.tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(0);
			flatNtuple.tau0_againstMuonLoose3 = evNtuple.tau_againstMuonLoose3->at(0);
			flatNtuple.tau0_againstMuonTight3 = evNtuple.tau_againstMuonTight3->at(0);
			flatNtuple.tau0_againstElectronVLooseMVA6 = evNtuple.tau_againstElectronVLooseMVA6->at(0);
			flatNtuple.tau0_againstElectronLooseMVA6 = evNtuple.tau_againstElectronLooseMVA6->at(0);
			flatNtuple.tau0_againstElectronMediumMVA6 = evNtuple.tau_againstElectronMediumMVA6->at(0);
			flatNtuple.tau0_againstElectronTightMVA6 = evNtuple.tau_againstElectronTightMVA6->at(0);
		}
		if (evNtuple.tau_pt->size()>1) {
			flatNtuple.tau1_pt = evNtuple.tau_pt->at(1);
			flatNtuple.tau1_eta = evNtuple.tau_eta->at(1);
			flatNtuple.tau1_phi = evNtuple.tau_phi->at(1);
			flatNtuple.tau1_E = evNtuple.tau_E->at(1);
			flatNtuple.tau1_charge = evNtuple.tau_charge->at(1);
			flatNtuple.tau1_dxy = evNtuple.tau_dxy->at(1);
			flatNtuple.tau1_dz = evNtuple.tau_dz->at(1);
			flatNtuple.tau1_decayModeFindingOldDMs = evNtuple.tau_decayModeFindingOldDMs->at(1);
			flatNtuple.tau1_decayModeFindingNewDMs = evNtuple.tau_decayModeFindingNewDMs->at(1);
			flatNtuple.tau1_byCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byCombinedIsolationDeltaBetaCorr3Hits->at(1);
			flatNtuple.tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byLooseCombinedIsolationDeltaBetaCorr3Hits->at(1);
			flatNtuple.tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits->at(1);
			flatNtuple.tau1_byTightCombinedIsolationDeltaBetaCorr3Hits = evNtuple.tau_byTightCombinedIsolationDeltaBetaCorr3Hits->at(1);
			flatNtuple.tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			flatNtuple.tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			flatNtuple.tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byTightIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			flatNtuple.tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT = evNtuple.tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT->at(1);
			flatNtuple.tau1_againstMuonLoose3 = evNtuple.tau_againstMuonLoose3->at(1);
			flatNtuple.tau1_againstMuonTight3 = evNtuple.tau_againstMuonTight3->at(1);
			flatNtuple.tau1_againstElectronVLooseMVA6 = evNtuple.tau_againstElectronVLooseMVA6->at(1);
			flatNtuple.tau1_againstElectronLooseMVA6 = evNtuple.tau_againstElectronLooseMVA6->at(1);
			flatNtuple.tau1_againstElectronMediumMVA6 = evNtuple.tau_againstElectronMediumMVA6->at(1);
			flatNtuple.tau1_againstElectronTightMVA6 = evNtuple.tau_againstElectronTightMVA6->at(1);
		}

		// jets
		if (evNtuple.jet_pt->size()>0) {
			flatNtuple.jet0_pt = evNtuple.jet_pt->at(0);
			flatNtuple.jet0_eta = evNtuple.jet_eta->at(0);
			flatNtuple.jet0_phi = evNtuple.jet_phi->at(0);
			flatNtuple.jet0_E = evNtuple.jet_E->at(0);
			flatNtuple.jet0_CSV = evNtuple.jet_csv->at(0);
		}
		if (evNtuple.jet_pt->size()>1) {
			flatNtuple.jet1_pt = evNtuple.jet_pt->at(1);
			flatNtuple.jet1_eta = evNtuple.jet_eta->at(1);
			flatNtuple.jet1_phi = evNtuple.jet_phi->at(1);
			flatNtuple.jet1_E = evNtuple.jet_E->at(1);
			flatNtuple.jet1_CSV = evNtuple.jet_csv->at(1);
		}
		if (evNtuple.jet_pt->size()>2) {
			flatNtuple.jet2_pt = evNtuple.jet_pt->at(2);
			flatNtuple.jet2_eta = evNtuple.jet_eta->at(2);
			flatNtuple.jet2_phi = evNtuple.jet_phi->at(2);
			flatNtuple.jet2_E = evNtuple.jet_E->at(2);
			flatNtuple.jet2_CSV = evNtuple.jet_csv->at(2);
		}
		if (evNtuple.jet_pt->size()>3) {
			flatNtuple.jet3_pt = evNtuple.jet_pt->at(3);
			flatNtuple.jet3_eta = evNtuple.jet_eta->at(3);
			flatNtuple.jet3_phi = evNtuple.jet_phi->at(3);
			flatNtuple.jet3_E = evNtuple.jet_E->at(3);
			flatNtuple.jet3_CSV = evNtuple.jet_csv->at(3);
		}

		// met
		flatNtuple.PFMET = evNtuple.PFMET;
		flatNtuple.PFMETphi = evNtuple.PFMETphi;
		flatNtuple.MHT = evNtuple.MHT;
		flatNtuple.metLD = evNtuple.metLD;

		
		tree_flat->Fill();
	}

	delete tree_in;
	
	// write to output file
	f_out->Write();

	// event count
	cout << "number of events with at least 1 preselected muons : " << "\t"
		 << tree_flat->GetEntries("n_presel_mu>0") << endl;
	cout << "number of events with at least 1 preselected electrons : " << "\t"
		 << tree_flat->GetEntries("n_presel_ele>0") << endl;
	cout << "number of events with at least 1 preselected taus : " << "\t"
		 << tree_flat->GetEntries("n_presel_tau>0") << endl;
	cout << "number of events with at least 1 preselected jets : " << "\t"
		 << tree_flat->GetEntries("n_presel_jet>0") << endl;
		
	// close files
	f_out->Close();
	f_in->Close();
}
