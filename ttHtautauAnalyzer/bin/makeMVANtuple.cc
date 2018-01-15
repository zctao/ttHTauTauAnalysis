// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>

#include "boost/program_options.hpp"

#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/miniLepton.h"
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/miniTau.h"
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/eventNtuple.h"
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/MVAVars.h"
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/Types_enum.h"
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/SFHelper.h"
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/TriggerHelper.h"
#include "ttHTauTauAnalysis/ttHtautauAnalyzer/interface/mvaNtuple.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>

//#pragma link C++ class std::vector < std::vector<int> >+;

Analysis_types getAnaType(std::string analysis_type)
{
	Analysis_types anaType = Analyze_NA;
		
	if (analysis_type=="1l2tau")
		anaType = Analyze_1l2tau;
	else if (analysis_type=="2lss1tau")
		anaType = Analyze_2lss1tau;
	else if (analysis_type=="3l1tau")
		anaType = Analyze_3l1tau;
	else {
		std::cerr << "Analysis type not supported! Available types are: 1l2tau, 2lss1tau, 3l1tau" << std::endl;
	}

	return anaType;
}

Selection_types getSelType(std::string selection_type)
{
	Selection_types selType = Selection_NA;
	
	if (selection_type=="signal_2lss1tau")
		selType = Signal_2lss1tau;
	else if (selection_type=="control_fake_2lss1tau")
		selType = Control_fake_2lss1tau;
	else if (selection_type=="control_2los1tau")
		selType = Control_2los1tau;
	else if (selection_type=="signal_1l2tau")
		selType = Signal_1l2tau;
	else if (selection_type=="control_fake_1l2tau")
		selType = Control_fake_1l2tau;
	else if (selection_type=="signal_3l1tau")
		selType = Signal_3l1tau;
	else if (selection_type=="control_fake_3l1tau")
		selType = Control_fake_3l1tau;
	else if (selection_type=="loose_1l2tau")
		selType = Loose_1l2tau;
	else if (selection_type=="loose_2lss1tau")
		selType = Loose_2lss1tau;
	else
		std::cerr << "Not valid selection type!" << std::endl;

	return selType;
}

int main(int argc, char** argv)
{
	using namespace std;
	namespace po = boost::program_options;

	string infile, outname, sample, analysis_type, selection_type, intree;
	bool requireMCMatching, evaluate, systematics, isdata;
	bool updateSF, looseSelection, setTreeWeight;

	po::options_description desc("Options");
	desc.add_options()
		("help,h", "produce help message")
		("infile,i", po::value<string>(&infile), "input file")
		//("outdir,o", po::value<string>(&outdir)->default_value("./"), "output directory")
		("outname,o", po::value<string>(&outname)->default_value("mvaNtuple.root"), "output file name")	
		//("sample,s", po::value<string>(&sample), "sample name")
		("anatype", po::value<string>(&analysis_type), "analysis type")
		("seltype", po::value<string>(&selection_type), "selection type")
		("treename,tn", po::value<string>(&intree)->default_value("ttHtaus/eventTree"))
		("isdata,d", po::value<bool>(&isdata)->default_value(false))
		("mc_matching,m", po::value<bool>(&requireMCMatching)->default_value(true))
		("evaluate,e", po::value<bool>(&evaluate)->default_value(false))
		("systematics,s", po::value<bool>(&systematics)->default_value(true))
		("update_sf,u", po::value<bool>(&updateSF)->default_value(false))
		("loose_selection,l", po::value<bool>(&looseSelection)->default_value(false))
		("tree_weight,w", po::value<bool>(&setTreeWeight)->default_value(false));
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc << endl;
		return 1;
	}

	//TString sample_ts = sample.c_str();
	//bool isdata = sample_ts.Contains("data");
	
	auto anaType = getAnaType(analysis_type);
	auto selType = getSelType(selection_type);

	// MVA variables
	MVAVars mvaVars(anaType);
	if (evaluate)
		mvaVars.set_up_tmva_reader();
	
	// Scale Factor Helper
	SFHelper sfhelper(anaType, selType, isdata);

	// Trigger Helper
	TriggerHelper trighelper(anaType, false);
	
	// open file and read tree
	cout << "Open root file : " << infile << endl;
	TFile* f_in = TFile::Open(infile.c_str());
	TTree* tree_in = (TTree*)f_in->Get(intree.c_str());

	eventNtuple evNtuple;
	evNtuple.set_branch_address(tree_in);

	// create output file
	//const string output_file = outdir+"mvaVars_"+sample+"_"+analysis_type+".root";
	const string output_file = outname;
	cout << "Output file created : " << output_file << endl;
	TFile* f_out = new TFile(output_file.c_str(), "recreate");
	TTree* tree_mva = new TTree("mva","mva");

	// mva ntuple
	mvaNtuple mvantuple(anaType, evaluate, systematics);
	mvantuple.setup_branches(tree_mva);
	
	// loop over events
	int nEntries = tree_in->GetEntries();
	for (int i = 0; i < nEntries; ++i) {
		tree_in->GetEntry(i);

		//////////////////////////////////////
		/// Additional selections
		//////////////////////////////////////

		// triggers/filters
		// also apply tighter/additional selection if needed
		if (anaType == Analyze_2lss1tau) {
			// tau charge
			if (not evNtuple.passTauCharge) continue;
		}
		else if (anaType == Analyze_1l2tau) {
			// mc match
			if (requireMCMatching and !isdata and
				not(evNtuple.isGenMatchedLep and evNtuple.isGenMatchedTau))
				continue;
			// tau pair charge
			if (not evNtuple.passTauCharge) continue;
		}
		else if (anaType == Analyze_3l1tau) {
			if (not evNtuple.passTauCharge) continue;
		}

		//////////////////////////////////////
		/// reconstruct object four momentums
		//////////////////////////////////////

		// loose muons and electrons, preselected taus, selected jets are stored
		// in the input ntuple

		/////////////////
		// leptons
		vector<miniLepton> leptons = evNtuple.buildLeptons(looseSelection);

		/////////////////
		// taus
		vector<miniTau> taus =
			evNtuple.buildTaus(looseSelection or selType==Control_fake_1l2tau);

		/////////////////
		// jets
		vector<TLorentzVector> jets
			= evNtuple.buildFourVectorJets(looseSelection);

		//////////////////////////////////////
		/// MVA variables
		//////////////////////////////////////

		mvaVars.compute_all_variables(leptons, taus, jets, evNtuple.PFMET,
									  evNtuple.PFMETphi, evNtuple.MHT,
									  evNtuple.n_btag_loose);

		mvantuple.nJet = mvaVars.nJet();
		mvantuple.avg_dr_jet = mvaVars.avg_dr_jet();

		if (anaType == Analyze_2lss1tau) {	
			
			mvantuple.mindr_lep0_jet = mvaVars.mindr_lep0_jet();
			mvantuple.mindr_lep1_jet = mvaVars.mindr_lep1_jet();	
			mvantuple.max_lep_eta = mvaVars.max_lep_eta();
			mvantuple.met = mvaVars.met();
			mvantuple.mht = mvaVars.mht();
			mvantuple.mT_met_lep0 = mvaVars.mT_met_lep0();
			mvantuple.lep0_conept = mvaVars.lep0_conept();
			mvantuple.lep1_conept = mvaVars.lep1_conept();
			mvantuple.dr_leps = mvaVars.dr_leps();
			mvantuple.tau0_pt = mvaVars.tau_pt();
			mvantuple.dr_lep0_tau = mvaVars.dr_lep0_tau();
			mvantuple.dr_lep1_tau = mvaVars.dr_lep1_tau();
			mvantuple.mvis_lep0_tau = mvaVars.mvis_lep0_tau();
			mvantuple.mvis_lep1_tau = mvaVars.mvis_lep1_tau();

			mvantuple.tau0_decaymode = mvaVars.tau0_decaymode();
			mvantuple.tau0_E = mvaVars.tau0_energy();
			mvantuple.tau0_upsilon = mvaVars.tau0_upsilon();
		}
		else if (anaType == Analyze_1l2tau) {
			mvantuple.HT = mvaVars.ht();
			mvantuple.tt_deltaR = mvaVars.dr_taus();
			mvantuple.tt_mvis = mvaVars.mvis_taus();
			mvantuple.tt_pt = mvaVars.pt_taus();
			mvantuple.max_dr_jet = mvaVars.max_dr_jet();
			mvantuple.tau0_pt = mvaVars.tau0_pt();
			mvantuple.tau1_pt = mvaVars.tau1_pt();
			mvantuple.ntags = evNtuple.n_btag_medium;
			mvantuple.ntags_loose = evNtuple.n_btag_loose;

			mvantuple.taup_decaymode = mvaVars.taup_decaymode();
			mvantuple.taum_decaymode = mvaVars.taum_decaymode();
			mvantuple.taup_E = mvaVars.taup_energy();
			mvantuple.taum_E = mvaVars.taum_energy();
			mvantuple.taup_upsilon = mvaVars.taup_upsilon();
			mvantuple.taum_upsilon = mvaVars.taum_upsilon();
		}
		else if (anaType == Analyze_3l1tau) {
			mvantuple.max_lep_eta = mvaVars.max_lep_eta();
			mvantuple.mindr_lep0_jet = mvaVars.mindr_lep0_jet();
			mvantuple.mindr_lep1_jet = mvaVars.mindr_lep1_jet();
			mvantuple.mindr_lep2_jet = mvaVars.mindr_lep2_jet();
			mvantuple.mT_met_lep0 = mvaVars.mT_met_lep0();
			mvantuple.lep0_conept = mvaVars.lep0_conept();
			mvantuple.lep1_conept = mvaVars.lep1_conept();
			mvantuple.lep2_conept = mvaVars.lep2_conept();

			mvantuple.tau0_decaymode = mvaVars.tau0_decaymode();
			mvantuple.tau0_E = mvaVars.tau0_energy();
			mvantuple.tau0_upsilon = mvaVars.tau0_upsilon();
		}

		if (evaluate) {
			mvantuple.mva_ttV = mvaVars.BDT_ttV();
			mvantuple.mva_ttbar = mvaVars.BDT_ttbar();
		}
		
		if (not isdata) { 
			mvantuple.isGenMatchedTau = evNtuple.isGenMatchedTau;
			mvantuple.HiggsDecayType = evNtuple.HiggsDecayType;
		}
		
		//////////////////////////////////////
		/// Event weights
		//////////////////////////////////////

		mvantuple.event_weight = evNtuple.event_weight;

		if (updateSF) {
			if (selType==Control_fake_1l2tau or selType==Control_fake_2lss1tau or
				selType==Control_fake_3l1tau) {
				mvantuple.event_weight =
					sfhelper.Get_FR_weight(leptons,taus);

				if (systematics) {
					if (selType==Control_fake_1l2tau) {
						mvantuple.event_weight_FRjt_normUp = sfhelper.Get_FR_weight(leptons,taus,"FRjt_normUp");
						mvantuple.event_weight_FRjt_normDown = sfhelper.Get_FR_weight(leptons,taus,"FRjt_normDown");
						mvantuple.event_weight_FRjt_shapeUp = sfhelper.Get_FR_weight(leptons,taus,"FRjt_shapeUp");
						mvantuple.event_weight_FRjt_shapeDown = sfhelper.Get_FR_weight(leptons,taus,"FRjt_shapeDown");
						
						// electron/muon fake rate systematics?
					}
					else { // fake 2lss1tau or 3l1tau
						mvantuple.event_weight_FRe_normUp =
							sfhelper.Get_FR_weight(leptons,taus,"FRe_normUp");
						mvantuple.event_weight_FRe_normDown =
							sfhelper.Get_FR_weight(leptons,taus,"FRe_normDown");
						mvantuple.event_weight_FRe_ptUp =
							sfhelper.Get_FR_weight(leptons,taus,"FRe_ptUp");
						mvantuple.event_weight_FRe_ptDown =
							sfhelper.Get_FR_weight(leptons,taus,"FRe_ptDown");
						mvantuple.event_weight_FRe_bUp =
							sfhelper.Get_FR_weight(leptons,taus,"FRe_bUp");
						mvantuple.event_weight_FRe_bDown =
							sfhelper.Get_FR_weight(leptons,taus,"FRe_bDown");
						mvantuple.event_weight_FRe_ecUp =
							sfhelper.Get_FR_weight(leptons,taus,"FRe_ecUp");
						mvantuple.event_weight_FRe_ecDown =
							sfhelper.Get_FR_weight(leptons,taus,"FRe_ecDown");
						mvantuple.event_weight_FRm_normUp =
							sfhelper.Get_FR_weight(leptons,taus,"FRm_normUp");
						mvantuple.event_weight_FRm_normDown =
							sfhelper.Get_FR_weight(leptons,taus,"FRm_normDown");
						mvantuple.event_weight_FRm_ptUp =
							sfhelper.Get_FR_weight(leptons,taus,"FRm_ptUp");
						mvantuple.event_weight_FRm_ptDown =
							sfhelper.Get_FR_weight(leptons,taus,"FRm_ptDown");
						mvantuple.event_weight_FRm_bUp =
							sfhelper.Get_FR_weight(leptons,taus,"FRm_bUp");
						mvantuple.event_weight_FRm_bDown =
							sfhelper.Get_FR_weight(leptons,taus,"FRm_bDown");
						mvantuple.event_weight_FRm_ecUp =
							sfhelper.Get_FR_weight(leptons,taus,"FRm_ecUp");
						mvantuple.event_weight_FRm_ecDown =
							sfhelper.Get_FR_weight(leptons,taus,"FRm_ecDown");
					}
				}
			}
			else if (selType==Control_2los1tau) {
				assert(taus.size()>0);
				mvantuple.event_weight =
					sfhelper.Get_ChargeFlipWeight(leptons, taus[0].charge());
			}
			else {
				float pu_weight = sfhelper.Get_PUWeight(evNtuple.npuTrue);
				float lepid_sf = sfhelper.Get_LeptonIDSF_weight(leptons);
				float tauid_sf = sfhelper.Get_TauIDSF_weight(taus);

				float hlt_sf = 1.;
				if (anaType==Analyze_1l2tau) {
					bool hlt1LTriggered =
						trighelper.pass_single_lep_triggers(evNtuple.triggerBits);
					bool hltXTriggered =
						trighelper.pass_leptau_cross_triggers(evNtuple.triggerBits);
					hlt_sf =
						sfhelper.Get_HLTSF_1l2tau(leptons[0],taus,hlt1LTriggered,
												  hltXTriggered);
				}
				else if (anaType==Analyze_2lss1tau) {
					hlt_sf =
						sfhelper.Get_HLTSF_2l1tau(evNtuple.lepCategory);
				}
				else if (anaType==Analyze_3l1tau) {
					hlt_sf = sfhelper.Get_HLTSF_3l1tau();
				}
				
				// not recomputable with current ntuple
				float mc_weight = evNtuple.MC_weight;
				float btag_sf = evNtuple.bTagSF_weight;

				mvantuple.event_weight =
					pu_weight * mc_weight * btag_sf * lepid_sf * tauid_sf * hlt_sf;
				if (systematics) {
					mvantuple.event_weight_thu_shape_x1Up = mvantuple.event_weight / mc_weight * evNtuple.MC_weight_scale_muF2;
					mvantuple.event_weight_thu_shape_x1Down = mvantuple.event_weight / mc_weight * evNtuple.MC_weight_scale_muF0p5;
					mvantuple.event_weight_thu_shape_y1Up = mvantuple.event_weight / mc_weight * evNtuple.MC_weight_scale_muR2;
					mvantuple.event_weight_thu_shape_y1Down = mvantuple.event_weight / mc_weight * evNtuple.MC_weight_scale_muR0p5;

					mvantuple.event_weight_btag_LFUp = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_LFUp;
					mvantuple.event_weight_btag_LFDown = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_LFDown;
					mvantuple.event_weight_btag_HFUp = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_HFUp;
					mvantuple.event_weight_btag_HFDown = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_HFDown;
					mvantuple.event_weight_btag_HFStats1Up = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_HFStats1Up;
					mvantuple.event_weight_btag_HFStats1Down = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_HFStats1Down;
					mvantuple.event_weight_btag_HFStats2Up = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_HFStats2Up;
					mvantuple.event_weight_btag_HFStats2Down = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_HFStats2Down;
					mvantuple.event_weight_btag_LFStats1Up = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_LFStats1Up;
					mvantuple.event_weight_btag_LFStats1Down = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_LFStats1Down;
					mvantuple.event_weight_btag_LFStats2Up = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_LFStats2Up;
					mvantuple.event_weight_btag_LFStats2Down = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_LFStats2Down;
					mvantuple.event_weight_btag_cErr1Up = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_cErr1Up;
					mvantuple.event_weight_btag_cErr1Down = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_cErr1Down;
					mvantuple.event_weight_btag_cErr2Up = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_cErr2Up;
					mvantuple.event_weight_btag_cErr2Down = mvantuple.event_weight / btag_sf * evNtuple.btagSF_weight_cErr2Down;

					if (anaType==Analyze_2lss1tau or anaType==Analyze_3l1tau) {
						mvantuple.event_weight_FRjt_normUp = mvantuple.event_weight / tauid_sf * sfhelper.Get_TauIDSF_weight(taus, "FRjt_normUp");
						mvantuple.event_weight_FRjt_normDown = mvantuple.event_weight / tauid_sf * sfhelper.Get_TauIDSF_weight(taus, "FRjt_normDown");
						mvantuple.event_weight_FRjt_shapeUp = mvantuple.event_weight / tauid_sf * sfhelper.Get_TauIDSF_weight(taus, "FRjt_shapeUp");
						mvantuple.event_weight_FRjt_shapeDown = mvantuple.event_weight / tauid_sf * sfhelper.Get_TauIDSF_weight(taus, "FRjt_shapeDown");
					}
					
				} // if (systematics)
				
			}
		} // if (updateSF)
		
		tree_mva->Fill();
		
	} // end of event loop

	if (setTreeWeight) {
		// Get SumGenWeights from histogram in the input file
		// Set output tree weight as 1./SumGenWeight
		auto h_SumGenWeight = (TH1D*)f_in->Get("ttHtaus/h_SumGenWeight");
		//auto h_SumGenWeight = (TH1D*)f_in->Get("ttHtaus/h_SumGenWeightxPU");
		double SumGenWeight = h_SumGenWeight->GetBinContent(1);
		tree_mva->SetWeight(1./SumGenWeight);
	}

	delete tree_in;

	f_out->Write();

	// close files
	f_out->Close();
	f_in->Close();
	
	return 0;
}
