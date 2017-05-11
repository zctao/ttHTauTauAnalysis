#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "../interface/eventNtuple.h"
#include "../interface/MVAVars.h"
#include "../interface/miniLepton.h"

#include <iostream>
#include <vector>
#include <algorithm>

void makeMVANtuple(const TString input_file, const TString output_file,
				   bool looseSelection)
{
	using namespace std;

	gROOT->ProcessLine(".L ../src/eventNtuple.cc+");
	gROOT->ProcessLine(".L ../src/miniLepton.cc+");
	gROOT->ProcessLine(".L ../src/MVAVars.cc+");

	// open file and read tree
	cout << "Opening root file : " << input_file << endl;
	TFile* f_in = new TFile(input_file, "read");
	TTree* tree_in = (TTree*)f_in->Get("ttHtaus/eventTree");

	eventNtuple evNtuple;
	evNtuple.set_branch_address(tree_in);

	// output
	cout << "Output file created : " << output_file << endl;
	TFile* f_out = new TFile(output_file, "recreate");

	TTree* tree_mva = new TTree("mva", "mva");

	float event_weight;
	int nJet;
	float mindr_lep0_jet;
	float mindr_lep1_jet;
	float avg_dr_jet;
	float max_lep_eta;
	float met;
	float mht;
	float mT_met_lep0;
	float lep0_conept;
	float lep1_conept;
	float dr_leps;
	float tau_pt;
	float dr_lep0_tau;
	float dr_lep1_tau;
	
	// set up branches
	tree_mva->Branch("event_weight", &event_weight);
	tree_mva->Branch("nJet", &nJet);
	tree_mva->Branch("mindr_lep0_jet", &mindr_lep0_jet);
	tree_mva->Branch("mindr_lep1_jet", &mindr_lep1_jet);
	tree_mva->Branch("avg_dr_jet", &avg_dr_jet);
	tree_mva->Branch("max_lep_eta", &max_lep_eta);
	tree_mva->Branch("met", &met);
	tree_mva->Branch("mht", &mht);
	tree_mva->Branch("mT_met_lep0", &mT_met_lep0);
	tree_mva->Branch("lep0_conept", &lep0_conept);
	tree_mva->Branch("lep1_conept", &lep1_conept);
	tree_mva->Branch("dr_leps", &dr_leps);
	tree_mva->Branch("tau_pt", &tau_pt);
	tree_mva->Branch("dr_lep0_tau", &dr_lep0_tau);
	tree_mva->Branch("dr_lep1_tau", &dr_lep1_tau);

	// loop over trees
	int nEntries = tree_in->GetEntries();
	for (int i = 0; i < nEntries; ++i) {
		tree_in->GetEntry(i);

		// apply tighter/additional selection if needed
		//TODO
		
		// build object four momentum
		// loose muons and electrons, preselected taus, selected jets are stored
		// in the ntuple

		/////////////////
		// leptons
		std::vector<miniLepton> leptons;
		
		for (unsigned int l = 0; l<evNtuple.mu_pt->size(); ++l) {
			// requrie at least fakeable id if not for training
			if (!looseSelection and !(evNtuple.mu_isfakeablesel->at(l))) continue;
			TLorentzVector mu;
			mu.SetPtEtaPhiE(evNtuple.mu_pt->at(l), evNtuple.mu_eta->at(l),
							evNtuple.mu_phi->at(l), evNtuple.mu_E->at(l));
			miniLepton lep(mu, evNtuple.mu_conept->at(l));
			leptons.push_back(lep);
		}

		for (unsigned int l = 0; l<evNtuple.ele_pt->size(); ++l) {
			// requrie at least fakeable id if not for training
			if (!looseSelection and !(evNtuple.ele_isfakeablesel->at(l))) continue;
			TLorentzVector ele;
			ele.SetPtEtaPhiE(evNtuple.ele_pt->at(l), evNtuple.ele_eta->at(l),
							 evNtuple.ele_phi->at(l), evNtuple.ele_E->at(l));
			miniLepton lep(ele, evNtuple.ele_conept->at(l));
			leptons.push_back(lep);
		}

		assert(leptons.size()>1);
		// sort by conept
		sort(leptons.begin(), leptons.end(), [] (miniLepton l1, miniLepton l2)
			 {return l1.conept()>l2.conept();} );

		/////////////////
		// taus
		std::vector<TLorentzVector> taus;
		
		for (unsigned int t = 0; t<evNtuple.tau_pt->size(); ++t) {
			// require tight tau if not for training
			if (!looseSelection and !(evNtuple.tau_idSelection->at(t))) continue;

			TLorentzVector tau;
			tau.SetPtEtaPhiE(evNtuple.tau_pt->at(t), evNtuple.tau_eta->at(t),
							 evNtuple.tau_phi->at(t), evNtuple.tau_E->at(t));
			taus.push_back(tau);
		}
		assert(taus.size()>0);
		
		/////////////////
		// jets
		std::vector<TLorentzVector> jets;

		for (unsigned int j = 0; j<evNtuple.jet_pt->size(); ++j) {
			TLorentzVector jet;
			jet.SetPtEtaPhiE(evNtuple.jet_pt->at(j), evNtuple.jet_eta->at(j),
							 evNtuple.jet_phi->at(j), evNtuple.jet_E->at(j));
			jets.push_back(jet);
		}

		MVAVars mvaVars(leptons, taus, jets, evNtuple.PFMET, evNtuple.PFMETphi,
						evNtuple.MHT);

		nJet = mvaVars.nJet();
		mindr_lep0_jet = mvaVars.mindr_lep0_jet();
		mindr_lep1_jet = mvaVars.mindr_lep1_jet();
		avg_dr_jet = mvaVars.avg_dr_jet();
		max_lep_eta = mvaVars.max_lep_eta();
		met = mvaVars.met();
		mht = mvaVars.mht();
		mT_met_lep0 = mvaVars.mT_met_lep0();
		lep0_conept = mvaVars.lep0_conept();
		lep1_conept = mvaVars.lep1_conept();
		dr_leps = mvaVars.dr_leps();
		tau_pt = mvaVars.tau_pt();
		dr_lep0_tau = mvaVars.dr_lep0_tau();
		dr_lep1_tau = mvaVars.dr_lep1_tau();

		event_weight = evNtuple.event_weight;

		tree_mva->Fill();
		
	} // end of tree loop

	delete tree_in;

	f_out->Write(); 
	
	// close files
	f_out->Close();
	f_in->Close();
}
