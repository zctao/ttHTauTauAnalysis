#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "../interface/eventNtuple.h"
#include "../interface/MVAVars.h"
#include "../interface/miniLepton.h"
#include "../interface/Types_enum.h"

#include <iostream>
#include <vector>
#include <algorithm>

#ifdef __ROOTCLING__
#pragma link C++ class std::vector < std::vector<int> >+;
#endif

void makeMVANtuple(const TString infile, const TString sample,
				   const TString analysis_type, bool looseSelection,
				   bool requireMCMatching, const TString outdir = "./",
				   const TString intree = "ttHtaus/eventTree")
{
	using namespace std;
	
	gROOT->ProcessLine(".L ../src/miniLepton.cc+");
	gROOT->ProcessLine(".L ../src/eventNtuple.cc+");
	gROOT->ProcessLine(".L ../src/MVAVars.cc+");

	// Set up analysis type
	Analysis_types anaType;
	if (analysis_type == "1l2tau")
		anaType = Analyze_1l2tau;
	else if (analysis_type == "2lss1tau")
		anaType = Analyze_2lss1tau;
	else if (analysis_type == "3l1tau")
		anaType = Analyze_3l1tau;
	else {
		cerr << "Analysis type not supported! Available types are: 1l2tau, 2lss1tau, 3l1tau" << endl;
		return;
	}	

	// Set up TMVA variables
	MVAVars mvaVars(anaType);
	
	// open file and read tree
	cout << "Opening root file : " << infile << endl;
	TFile* f_in = TFile::Open(infile);
	TTree* tree_in = (TTree*)f_in->Get(intree);

	eventNtuple evNtuple;
	evNtuple.set_branch_address(tree_in);

	// output
	const TString output_file = outdir+"mvaVars_"+sample+"_"+analysis_type+".root";
	cout << "Output file created : " << output_file << endl;
	TFile* f_out = new TFile(output_file, "recreate");

	TTree* tree_mva = new TTree("mva", "mva");

	// MVA variables
	float event_weight;
	int nJet;
	float mindr_lep0_jet;
	float mindr_lep1_jet;
	float mindr_lep2_jet;
	float avg_dr_jet;
	float max_lep_eta;
	float met;
	float mht;
	float mT_met_lep0;
	float lep0_conept;
	float lep1_conept;
	float lep2_conept;
	float dr_leps;
	float tau0_pt;
	float tau1_pt;
	float dr_lep0_tau;
	float dr_lep1_tau;
	float mvis_lep0_tau;
	float mvis_lep1_tau;
	float tt_deltaR;
	float tt_mvis;
	float tt_pt;
	float max_dr_jet;
	float HT;
	int ntags;
	int ntags_loose;

	int tau0_decaymode;
	int tau1_decaymode;
	float tau0_E;
	float tau1_E;
	float tau0_upsilon;
	float tau1_upsilon;
	
	// set up branches
	tree_mva->Branch("event_weight", &event_weight);
	tree_mva->Branch("nJet", &nJet);
	tree_mva->Branch("avg_dr_jet", &avg_dr_jet);

	tree_mva->Branch("tau0_decaymode", &tau0_decaymode);
	tree_mva->Branch("tau0_E", &tau0_E);
	tree_mva->Branch("tau0_upsilon", &tau0_upsilon);
	
	if (anaType == Analyze_2lss1tau) {		
		tree_mva->Branch("mindr_lep0_jet", &mindr_lep0_jet);
		tree_mva->Branch("mindr_lep1_jet", &mindr_lep1_jet);		
		tree_mva->Branch("max_lep_eta", &max_lep_eta);
		tree_mva->Branch("met", &met);
		tree_mva->Branch("mht", &mht);
		tree_mva->Branch("mT_met_lep0", &mT_met_lep0);
		tree_mva->Branch("lep0_conept", &lep0_conept);
		tree_mva->Branch("lep1_conept", &lep1_conept);
		tree_mva->Branch("dr_leps", &dr_leps);
		tree_mva->Branch("tau_pt", &tau0_pt);
		tree_mva->Branch("dr_lep0_tau", &dr_lep0_tau);
		tree_mva->Branch("dr_lep1_tau", &dr_lep1_tau);
		tree_mva->Branch("mvis_lep0_tau", &mvis_lep0_tau);
		tree_mva->Branch("mvis_lep1_tau", &mvis_lep1_tau);
	}
	else if (anaType == Analyze_1l2tau) {
		tree_mva->Branch("ht", &HT);
		tree_mva->Branch("tt_deltaR", &tt_deltaR);
		tree_mva->Branch("tt_mvis", &tt_mvis);
		tree_mva->Branch("tt_sumpt", &tt_pt);
		tree_mva->Branch("max_dr_jet", &max_dr_jet);
		tree_mva->Branch("tau0_pt", &tau0_pt);
		tree_mva->Branch("tau1_pt", &tau1_pt);
		tree_mva->Branch("ntags", &ntags);
		tree_mva->Branch("ntags_loose", &ntags_loose);
		
		tree_mva->Branch("tau1_decaymode", &tau1_decaymode);
		tree_mva->Branch("tau1_E", &tau1_E);
		tree_mva->Branch("tau1_upsilon", &tau1_upsilon);
	}
	else if (anaType == Analyze_3l1tau) {
		tree_mva->Branch("max_lep_eta", &max_lep_eta);
		tree_mva->Branch("mindr_lep0_jet", &mindr_lep0_jet);
		tree_mva->Branch("mindr_lep1_jet", &mindr_lep1_jet);
		tree_mva->Branch("mindr_lep2_jet", &mindr_lep2_jet);
		tree_mva->Branch("mT_met_lep0", &mT_met_lep0);
		tree_mva->Branch("lep0_conept", &lep0_conept);
		tree_mva->Branch("lep1_conept", &lep1_conept);
		tree_mva->Branch("lep2_conept", &lep2_conept);
	}
	
	// loop over trees
	int nEntries = tree_in->GetEntries();
	for (int i = 0; i < nEntries; ++i) {
		tree_in->GetEntry(i);

		// apply tighter/additional selection if needed
		if (anaType == Analyze_2lss1tau) {
			// tau charge
			if (not evNtuple.passTauCharge) continue;
		}
		else if (anaType == Analyze_1l2tau) {
			// mc match
			if (requireMCMatching and
				not(evNtuple.isGenMatchedLep and evNtuple.isGenMatchedTau))
				continue;
			// tau pair charge
			if (not evNtuple.passTauCharge) continue;
		}
		else if (anaType == Analyze_3l1tau) {
			if (not evNtuple.passTauCharge) continue;
		}

		//if (not looseSelection) {
			// TODO
			// Full selection with tight objects
		//}
		
		// build object four momentum
		// loose muons and electrons, preselected taus, selected jets are stored
		// in the ntuple

		/////////////////
		// leptons
		std::vector<miniLepton> leptons = evNtuple.buildLeptons(looseSelection);

		/////////////////
		// taus
		std::vector<int> taus_decaymode;
		std::vector<TLorentzVector> taus
			= evNtuple.buildFourVectorTaus(taus_decaymode, looseSelection);
		
		std::vector<TLorentzVector> taudaugs_charged
			= evNtuple.buildFourVectorTauDaugsCharged(looseSelection);
		std::vector<TLorentzVector> taudaugs_neutral
			= evNtuple.buildFourVectorTauDaugsNeutral(looseSelection);
		
		/////////////////
		// jets
		std::vector<TLorentzVector> jets
			= evNtuple.buildFourVectorJets(looseSelection);

		mvaVars.compute_all_variables(leptons, taus, jets, evNtuple.PFMET,
									  evNtuple.PFMETphi, evNtuple.MHT,
									  evNtuple.n_btag_loose);
		mvaVars.compute_taudecay_variables(taus,taudaugs_charged,taudaugs_neutral,
										   taus_decaymode);

		nJet = mvaVars.nJet();
		avg_dr_jet = mvaVars.avg_dr_jet();

		tau0_decaymode = mvaVars.tau0_decaymode();
		tau0_E = mvaVars.tau0_energy();
		tau0_upsilon = mvaVars.tau0_upsilon();
		
		if (anaType == Analyze_2lss1tau) {	
			
			mindr_lep0_jet = mvaVars.mindr_lep0_jet();
			mindr_lep1_jet = mvaVars.mindr_lep1_jet();	
			max_lep_eta = mvaVars.max_lep_eta();
			met = mvaVars.met();
			mht = mvaVars.mht();
			mT_met_lep0 = mvaVars.mT_met_lep0();
			lep0_conept = mvaVars.lep0_conept();
			lep1_conept = mvaVars.lep1_conept();
			dr_leps = mvaVars.dr_leps();
			tau0_pt = mvaVars.tau_pt();
			dr_lep0_tau = mvaVars.dr_lep0_tau();
			dr_lep1_tau = mvaVars.dr_lep1_tau();
			mvis_lep0_tau = mvaVars.mvis_lep0_tau();
			mvis_lep1_tau = mvaVars.mvis_lep1_tau();
		}
		else if (anaType == Analyze_1l2tau) {
			HT = mvaVars.ht();
			tt_deltaR = mvaVars.dr_taus();
			tt_mvis = mvaVars.mvis_taus();
			tt_pt = mvaVars.pt_taus();
			max_dr_jet = mvaVars.max_dr_jet();
			tau0_pt = mvaVars.tau0_pt();
			tau1_pt = mvaVars.tau1_pt();
			ntags = evNtuple.n_btag_medium;
			ntags_loose = evNtuple.n_btag_loose;

			tau1_decaymode = mvaVars.tau1_decaymode();
			tau1_E = mvaVars.tau1_energy();
			tau1_upsilon = mvaVars.tau1_upsilon();
		}
		else if (anaType == Analyze_3l1tau) {

		}

		event_weight = evNtuple.event_weight;

		tree_mva->Fill();
		
	} // end of tree loop

	delete tree_in;

	f_out->Write(); 
	
	// close files
	f_out->Close();
	f_in->Close();
}
