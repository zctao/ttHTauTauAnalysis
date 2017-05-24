#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "../interface/eventNtuple.h"
#include "../interface/miniLepton.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

void dumpEvents(const TString infile, const TString outfile = "EventDump.py",
				bool looseSelection = false, int maxEvents = -1)
{
	using namespace std;

	gROOT->ProcessLine(".L ../src/eventNtuple.cc+");
	gROOT->ProcessLine(".L ../src/miniLepton.cc+");

	// open file and read tree
	cout << "Opening root file : " << infile << endl;
	TFile* f_in = new TFile(infile, "read");
	TTree* tree_in = (TTree*)f_in->Get("ttHtaus/eventTree");

	eventNtuple evNtuple;
	evNtuple.set_branch_address(tree_in);

	// output
	ofstream eventDump;
    eventDump.open(outfile);
	cout << "Output file created : " << outfile << endl;

	eventDump << "\"\"\"Event list\n";
	eventDump << "Per event :\n";
	eventDump << " - leptons\n";
	eventDump << " - tau\n";
	eventDump << " - MET\n";
	eventDump << " - MET Cov. matrix\n";
	eventDump << " - 2 b-jets\n";
	eventDump << " - Untagged jets\n";
	eventDump << "\"\"\"\n\n";
	eventDump << "Events = [\n";

	// loop over events in tree
	int nEntries = tree_in->GetEntries();

	for (int i = 0; i < nEntries; ++i) {

		if (maxEvents>=0 and i>=maxEvents) break;

		tree_in->GetEntry(i);
		
		// apply tighter selection here if needed

		// build four momentums
		//auto leptons = evNtuple.buildFourVectorLeps(looseSelection);
		auto leptons = evNtuple.buildLeptons(looseSelection);
		auto taus = evNtuple.buildFourVectorTaus(looseSelection);
		auto jets = evNtuple.buildFourVectorJets(looseSelection);
		auto jets_btag = evNtuple.buildFourVectorBtagJets(looseSelection);
		std::vector<TLorentzVector> jets_untag;
		set_difference(jets.begin(),jets.end(),jets_btag.begin(),jets_btag.end(),
					   inserter(jets_untag, jets_untag.end()),
					   [] (TLorentzVector l1, TLorentzVector l2)
					   {return l1.Pt() > l2.Pt();});
		assert(leptons.size()>1);
		assert(taus.size()>0);
		assert(jets_btag.size()>1);
		
		eventDump << "{\n";
		// leptons
		eventDump << "\'_recolep_sel_px\': ["
				  << leptons[0].p4().Px()<< "," << leptons[1].p4().Px() << "],\n";
		eventDump << "\'_recolep_sel_py\': ["
				  << leptons[0].p4().Py()<< "," << leptons[1].p4().Py() << "],\n";
		eventDump << "\'_recolep_sel_pz\': ["
				  << leptons[0].p4().Pz()<< "," << leptons[1].p4().Pz() << "],\n";
		eventDump << "\'_recolep_sel_e\': ["
				  << leptons[0].p4().E() << "," << leptons[1].p4().E() << "],\n";
		eventDump << "\'_recolep_sel_pdg\': ["
				  << leptons[0].pdgId() << "," << leptons[1].pdgId() << "],\n";

		// tau
		eventDump << "\'_recotauh_sel_px\': [" << taus[0].Px() << "],\n";
		eventDump << "\'_recotauh_sel_py\': [" << taus[0].Py() << "],\n";
		eventDump << "\'_recotauh_sel_pz\': [" << taus[0].Pz() << "],\n";
		eventDump << "\'_recotauh_sel_e\': [" << taus[0].E() << "],\n";
		eventDump << "\'_recotauh_sel_decayMode\': ["
				  << evNtuple.tau_decayMode->at(0) << "],\n";

		// MET
		eventDump << "\'_PFMET\': " << evNtuple.PFMET
				  << ", \'_PFMET_phi\': " << evNtuple.PFMETphi << ",\n";
		eventDump << "\'_PFMETCov00\': " << evNtuple.METCov00 << ", "
				  << "\'_PFMETCov11\': " << evNtuple.METCov11 << ",\n";
		eventDump << "\'_PFMETCov01\': " << evNtuple.METCov01 << ", "
				  << "\'_PFMETCov10\': " << evNtuple.METCov10 << ",\n";

		// btagged jets
		eventDump << "\'_recoPFJet_btag_px\': [" << jets_btag[0].Px() << ","
				  << jets_btag[1].Px() << "],\n";
		eventDump << "\'_recoPFJet_btag_py\': [" << jets_btag[0].Py() << ","
				  << jets_btag[1].Py() << "],\n";
		eventDump << "\'_recoPFJet_btag_pz\': [" << jets_btag[0].Pz() << ","
				  << jets_btag[1].Pz() << "],\n";
		eventDump << "\'_recoPFJet_btag_e\': [" << jets_btag[0].E() << ","
				  << jets_btag[1].E() << "],\n";

		// untagged jets
		eventDump << "\'_n_recoPFJet_untag\': " << jets_untag.size() << ",\n";
		eventDump << "\'_recoPFJet_untag_px\': [";
		size_t ipx = 0;
		for (const auto & jet : jets_untag) {
			eventDump << jet.Px();
			++ipx;
			if (ipx < jets_untag.size()) eventDump << ",";
		}
		eventDump << "],\n";
		eventDump << "\'_recoPFJet_untag_py\': [";
		size_t ipy = 0;
		for (const auto & jet : jets_untag) {
			eventDump << jet.Py();
			++ipy;
			if (ipy < jets_untag.size()) eventDump << ",";
		}
		eventDump << "],\n";
		eventDump << "\'_recoPFJet_untag_pz\': [";
		size_t ipz = 0;
		for (const auto & jet : jets_untag) {
			eventDump << jet.Pz();
			++ipz;
			if (ipz < jets_untag.size()) eventDump << ",";
		}
		eventDump << "],\n";
		eventDump << "\'_recoPFJet_untag_e\': [";
		size_t ie = 0;
		for (const auto & jet : jets_untag) {
			eventDump << jet.E();
			++ie;
			if (ie < jets_untag.size()) eventDump << ",";
		}
		eventDump << "]\n";
		
		eventDump << "},\n";
		
	} // end of tree loop

	eventDump << "]";
	eventDump.close();		
}
