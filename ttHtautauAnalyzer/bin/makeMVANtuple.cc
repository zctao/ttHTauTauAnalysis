// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TParameter.h>

#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"

#include "../interface/miniLepton.h"
#include "../interface/miniTau.h"
#include "../interface/miniJet.h"
#include "../interface/eventNtuple.h"
#include "../interface/Types_enum.h"
#include "../interface/SFHelper.h"
#include "../interface/TriggerHelper.h"
#include "../interface/EventSelector.h"
#include "../interface/mvaNtuple.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <string>
#include <tuple>

//#pragma link C++ class std::vector < std::vector<int> >+;

void Set_FR_weights(Selection_types, SFHelper&, eventNtuple&, mvaNtuple&,
					const std::vector<miniLepton>&, const std::vector<miniTau>&,
					bool);
void Set_SR_weights(Analysis_types,SFHelper&, eventNtuple&, mvaNtuple&,
					TriggerHelper&,
					const std::vector<miniLepton>&, const std::vector<miniTau>&,
					const std::vector<miniJet>&, bool, std::string);

int main(int argc, char** argv)
{
	using namespace std;
	namespace po = boost::program_options;

	string infile, outname, sample, analysis_type, selection_type, intree, version;
	string enCorr;
	string mem_filename, mem_treename;
	float xsection;
	bool systematics, isdata, evaluateMVA, genMatching;
	bool updateSF, setTreeWeight, addMEM;
	bool requireTrigger, looseSelection;
	//string tauWP;

	po::options_description desc("Options");
	desc.add_options()
		("help,h", "produce help message")
		("infile,i", po::value<string>(&infile), "input file")
		//("outdir,o", po::value<string>(&outdir)->default_value("./"), "output directory")
		("outname,o", po::value<string>(&outname)->default_value("mvaNtuple.root"), "output file name")	
		("sample,s", po::value<string>(&sample), "sample name")
		("anatype", po::value<string>(&analysis_type), "analysis type")
		("seltype", po::value<string>(&selection_type), "selection type")
		("version,v", po::value<string>(&version)->default_value("2017"), "analysis version")
		("correction,c", po::value<string>(&enCorr)->default_value("NA"), "Energy correction for tau and jets: NA, JESUp, JESDown, JERUp, JERDown, TESUp, TESDown")
		("xsection,x", po::value<float>(&xsection)->default_value(-99.),"cross section of the sample")
		("treename", po::value<string>(&intree)->default_value("ttHtaus/eventTree"))
		("isdata,d", po::value<bool>(&isdata)->default_value(false))
		("trigger,t", po::value<bool>(&requireTrigger)->default_value(true))
		("add_mem", po::value<bool>(&addMEM)->default_value(false))
		("mem_filename", po::value<string>(&mem_filename)->default_value("mem_output.root"))
		("mem_treename", po::value<string>(&mem_treename)->default_value("mem"))
		("gen_matching,m", po::value<bool>(&genMatching)->default_value(false))
		("evaluateMVA,e", po::value<bool>(&evaluateMVA)->default_value(true))
		("systematics", po::value<bool>(&systematics)->default_value(true))
		("update_sf,u", po::value<bool>(&updateSF)->default_value(true))
		("loose_selection,l", po::value<bool>(&looseSelection)->default_value(false))
		//("tauWP", po::value<string>(&tauWP)->default_value("-"))
		("tree_weight,w", po::value<bool>(&setTreeWeight)->default_value(false));
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		cout << desc << endl;
		return 1;
	}

	Analysis_types anaType = Types_enum::getAnaType(analysis_type);
	Selection_types selType = Types_enum::getSelType(selection_type);

	bool CR = boost::algorithm::contains(selection_type,"_ttW") or
		boost::algorithm::contains(selection_type,"_ttZ") or
		boost::algorithm::contains(selection_type,"_WZ");
	
	assert(enCorr=="NA" or enCorr=="JESUp" or enCorr=="JESDown" or
           enCorr=="JERUp" or enCorr=="JERDown" or
           enCorr=="TESUp" or enCorr=="TESDown");
	
	//////////////////////////////////////////////
	// open file and read tree
	cout << "Open root file : " << infile << endl;
	TFile* f_in = TFile::Open(infile.c_str());
	TTree* tree_in = (TTree*)f_in->Get(intree.c_str());

	eventNtuple evNtuple;
	evNtuple.set_branch_address(tree_in);

	//////////////////////////////////////////////
	// create output file
	const string output_file = outname;
	cout << "Output file created : " << output_file << endl;
	TFile* f_out = new TFile(output_file.c_str(), "recreate");
	TTree* tree_mva = new TTree("mva","mva");

	// mva ntuple
	bool doHTT = evaluateMVA;
	mvaNtuple mvantuple(anaType, systematics, version, doHTT, evaluateMVA, CR);
	mvantuple.setup_branches(tree_mva);

	//////////////////////////////////////////////
	// Trigger Helper
	TriggerHelper trighelper(Analysis_types::Analyze_inclusive, false);

	//////////////////////////////////////////////
	// Scale Factor Helper
	SFHelper sfhelper(anaType, selType, sample, isdata);

	//////////////////////////////////////////////
	// event selector
	EventSelector evt_selector(false, !isdata, genMatching, looseSelection);

	//////////////////////////////////////////////
	// total number of events processed
	auto h_nProcessed = (TH1D*)f_in->Get("ttHtaus/h_nProcessed");
	double nProcessed = h_nProcessed->GetBinContent(1);
	// sum of genWeights for all events processed
	auto h_SumGenWeight = (TH1D*)f_in->Get("ttHtaus/h_SumGenWeight");
	//auto h_SumGenWeight = (TH1D*)f_in->Get("ttHtaus/h_SumGenWeightxPU");
	double SumGenWeight = h_SumGenWeight->GetBinContent(1);
	if (h_SumGenWeight->GetEntries()==0) // empty h_SumGenWeight in data ntuples
		SumGenWeight = 1.;

	TParameter<double> xspar("xsection", xsection);
	TParameter<double> sumgenwpar("SumGenWeight", SumGenWeight);
	
	tree_mva->GetUserInfo()->Add(&xspar);
	tree_mva->GetUserInfo()->Add(&sumgenwpar);
	
	// loop over events
	int nEntries = tree_in->GetEntries();
	//std::cout << "nEntries : " << nEntries << std::endl;

	for (int i = 0; i < nEntries; ++i) {
		tree_in->GetEntry(i);

		// reject one LS: https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3067.html
		if (evNtuple.run==305366 and evNtuple.ls==935) continue;
		
		//////////////////////////////////////
		/// Reconstruct objects for event selection
		//////////////////////////////////////
		auto leptons_loose = evNtuple.buildLeptons('L');  // loose
		auto leptons_fakeable = evNtuple.buildLeptons('F');  // fakeable
		auto leptons_tight = evNtuple.buildLeptons('T');  // tight

		auto taus_fakeable = evNtuple.buildTaus(true, anaType, enCorr.c_str(), isdata);
		auto taus_tight = evNtuple.buildTaus(false, anaType, enCorr.c_str(), isdata);

		//auto jets = evNtuple.buildJets();
		// Jet cleaning based on analysis type
		auto jets = evNtuple.buildCleanedJets(0.4, anaType, selType,
											  &leptons_fakeable, &taus_fakeable,
											  enCorr.c_str(), true);

		int nbtags_loose, nbtags_medium;
		std::tie(nbtags_loose, nbtags_medium) = evNtuple.count_btags(jets);
		assert(nbtags_loose >= nbtags_medium);

        // FIXME: propagate jet or tau energy correction to MET
		auto metp4 = evNtuple.buildFourVectorMET();
		float met = metp4.Pt();
		float mht = evNtuple.computeMHT(leptons_fakeable, taus_fakeable, jets);
		float metld = 0.00397 * met + 0.00265 * mht;
		
		//////////////////////////////////////
		/// Event selections
		//////////////////////////////////////
		bool passEvtSel =  evt_selector.pass_full_event_selection(
			anaType, selType, leptons_loose, leptons_fakeable, leptons_tight,
			taus_fakeable, taus_tight, jets.size(), nbtags_loose,
			nbtags_medium, metld);

		if (not passEvtSel) continue;

		// trigger paths
		int nfakeableEle = mvantuple.count_electrons(leptons_fakeable);
		int nfakeableMu = mvantuple.count_muons(leptons_fakeable);
		//bool passHLT =
		//	evt_selector.pass_hlt_paths(anaType,&trighelper,evNtuple.triggerBits) and
		//	evt_selector.pass_hlt_match(anaType,&trighelper,evNtuple.triggerBits,
		//								nfakeableEle, nfakeableMu);
		//if (not passHLT and requireTrigger) continue;

		bool passHLTandFilters = evt_selector.pass_hlt_and_filters(
			anaType, &trighelper, evNtuple.triggerBits, nfakeableEle, nfakeableMu,
			evNtuple.filterBits, isdata);
		if (not passHLTandFilters and requireTrigger) continue;
		
		//////////////////////////////////////
		/// MVA variables
		//////////////////////////////////////

		const std::vector<miniLepton> *leptons = &leptons_fakeable;
		// TODO: if loose selection

		const std::vector<miniTau> *taus = &taus_tight;
		if (anaType==Analyze_1l2tau or anaType==Analyze_2l2tau)
            // since taus are included in the fake background estimation
			taus = &taus_fakeable;
		
		mvantuple.compute_all_variables(*leptons, *taus, jets,
										metp4.Pt(), metp4.Phi(), mht,
										nbtags_loose, nbtags_medium);
		
		//////////////////////////////////////
		/// Event ID
		//////////////////////////////////////
		mvantuple.run = evNtuple.run;
		mvantuple.lumi = evNtuple.ls;
		mvantuple.nEvent = evNtuple.nEvent;

		if (not isdata) {
			mvantuple.HiggsDecayType = evNtuple.HiggsDecayType;
			mvantuple.xsection_weight = xsection/nProcessed; // deprecated
			mvantuple.xsection_weight_gen = xsection/SumGenWeight; // deprecated
			mvantuple.isGenMatchedLep =
				evt_selector.pass_MCMatch_Leps(anaType, *leptons);
			mvantuple.isGenMatchedTau =
				evt_selector.pass_MCMatch_Taus(anaType, *taus);
			mvantuple.isGenPhotonMatched =
				evt_selector.is_MCMatch_Photon(anaType, *leptons);
		}
		
		//////////////////////////////////////
		/// Event weights and scale factors
		//////////////////////////////////////

		if (not updateSF) {
			mvantuple.event_weight = evNtuple.event_weight;
			tree_mva->Fill();
			continue;
		}

		if (isdata) {
			mvantuple.event_weight = 1;
		}
		else {
			// bad naming: should be Set_MC_weights
			Set_SR_weights(anaType, sfhelper, evNtuple, mvantuple, trighelper,
						   *leptons, *taus, jets, systematics, enCorr);
		}

		// if not signal region selection, overwrite event_weight with FF method
		if (selType==Application_Fake_1l2tau or selType==Application_Fake_2lss1tau or
			selType==Application_Fake_3l1tau or selType==Application_Fake_2l2tau or
			selType==Control_FakeAR_1l2tau or selType==Control_FakeAR_2lss1tau or
			selType==Control_FakeAR_3l1tau or selType==Control_FakeAR_2l2tau or
			selType==Control_FakeAR_ttW or selType==Control_FakeAR_ttZ or
			selType==Control_FakeAR_WZ) {
			// Fake application region weights
			Set_FR_weights(selType, sfhelper, evNtuple, mvantuple, *leptons, *taus,
						   systematics);
		}
		else if (selType==Application_Flip_2lss1tau or
				 selType==Control_FlipAR_2lss1tau) {
			// Charge flip weights
			mvantuple.event_weight =
				sfhelper.Get_ChargeFlipWeight(*leptons, taus->at(0).charge());
		}
		else if (selType==Control_FlipAR_ttW) {
			mvantuple.event_weight = sfhelper.Get_ChargeFlipWeight(*leptons);
		}
		
		tree_mva->Fill();
		
	} // end of event loop
	
	if (setTreeWeight and !isdata) {
		// Set output tree weight as 1./SumGenWeight
		tree_mva->SetWeight(1./SumGenWeight);
		//tree_mva->SetWeight(1./nProcessed);
	}

	delete tree_in;

	// Add mem tree as friend
	if (addMEM) {
		tree_mva->AddFriend(mem_treename.c_str(), mem_filename.c_str());
	}

	f_out->Write();

	// close files
	f_out->Close();
	f_in->Close();
	
	return 0;
}

void Set_FR_weights(Selection_types seltype, SFHelper& sfhelper,
					eventNtuple& evntuple, mvaNtuple& mvantuple,
					const std::vector<miniLepton>& leptons,
					const std::vector<miniTau>& taus, bool syst)
{
	assert(sfhelper.getSelType()==seltype);
	
	mvantuple.event_weight = sfhelper.Get_FR_weight(leptons, taus);

	if (syst) {
		mvantuple.event_weight_FRe_normUp =
			sfhelper.Get_FR_weight(leptons,taus,"FRe_normUp");
		mvantuple.event_weight_FRe_normDown =
			sfhelper.Get_FR_weight(leptons,taus,"FRe_normDown");
		mvantuple.event_weight_FRe_ptUp =
			sfhelper.Get_FR_weight(leptons,taus,"FRe_ptUp");
		mvantuple.event_weight_FRe_ptDown =
			sfhelper.Get_FR_weight(leptons,taus,"FRe_ptDown");
		mvantuple.event_weight_FRe_beUp =
			sfhelper.Get_FR_weight(leptons,taus,"FRe_beUp");
		mvantuple.event_weight_FRe_beDown =
			sfhelper.Get_FR_weight(leptons,taus,"FRe_beDown");
		//mvantuple.event_weight_FRe_bUp =
		//	sfhelper.Get_FR_weight(leptons,taus,"FRe_bUp");
		//mvantuple.event_weight_FRe_bDown =
		//	sfhelper.Get_FR_weight(leptons,taus,"FRe_bDown");
		//mvantuple.event_weight_FRe_ecUp =
		//	sfhelper.Get_FR_weight(leptons,taus,"FRe_ecUp");
		//mvantuple.event_weight_FRe_ecDown =
		//	sfhelper.Get_FR_weight(leptons,taus,"FRe_ecDown");
		mvantuple.event_weight_FRm_normUp =
			sfhelper.Get_FR_weight(leptons,taus,"FRm_normUp");
		mvantuple.event_weight_FRm_normDown =
			sfhelper.Get_FR_weight(leptons,taus,"FRm_normDown");
		mvantuple.event_weight_FRm_ptUp =
			sfhelper.Get_FR_weight(leptons,taus,"FRm_ptUp");
		mvantuple.event_weight_FRm_ptDown =
			sfhelper.Get_FR_weight(leptons,taus,"FRm_ptDown");
		mvantuple.event_weight_FRm_beUp =
			sfhelper.Get_FR_weight(leptons,taus,"FRm_beUp");
		mvantuple.event_weight_FRm_beDown =
			sfhelper.Get_FR_weight(leptons,taus,"FRm_beDown");
		//mvantuple.event_weight_FRm_bUp =
		//	sfhelper.Get_FR_weight(leptons,taus,"FRm_bUp");
		//mvantuple.event_weight_FRm_bDown =
		//	sfhelper.Get_FR_weight(leptons,taus,"FRm_bDown");
		//mvantuple.event_weight_FRm_ecUp =
		//	sfhelper.Get_FR_weight(leptons,taus,"FRm_ecUp");
		//mvantuple.event_weight_FRm_ecDown =
		//	sfhelper.Get_FR_weight(leptons,taus,"FRm_ecDown");

		mvantuple.event_weight_FR_TT =
			sfhelper.Get_FR_weight(leptons,taus,"FR_TT");
		mvantuple.event_weight_FR_QCD =
			sfhelper.Get_FR_weight(leptons,taus,"FR_QCD");
		mvantuple.event_weight_FR_el_QCD_mu_TT =
			sfhelper.Get_FR_weight(leptons,taus,"FR_el_QCD_mu_TT");
		mvantuple.event_weight_FR_el_TT_mu_QCD =
			sfhelper.Get_FR_weight(leptons,taus,"FR_el_TT_mu_QCD");
		
		if (seltype==Application_Fake_1l2tau or seltype==Application_Fake_2l2tau) {
			mvantuple.event_weight_FRjt_normUp =
				sfhelper.Get_FR_weight(leptons,taus,"FRjt_normUp");
			mvantuple.event_weight_FRjt_normDown =
				sfhelper.Get_FR_weight(leptons,taus,"FRjt_normDown");
			mvantuple.event_weight_FRjt_shapeUp =
				sfhelper.Get_FR_weight(leptons,taus,"FRjt_shapeUp");
			mvantuple.event_weight_FRjt_shapeDown =
				sfhelper.Get_FR_weight(leptons,taus,"FRjt_shapeDown");
		}
	}  // syst
}

void Set_SR_weights(Analysis_types anatype, SFHelper& sfhelper,
					eventNtuple& evntuple, mvaNtuple& mvantuple,
					TriggerHelper& trighelper,
					const std::vector<miniLepton>& leptons,
					const std::vector<miniTau>& taus,
					const std::vector<miniJet>& jets, bool syst,
					std::string JES)
{
	assert(sfhelper.getAnaType()==anatype);
	
	// signal region and MC samples

	// Gen weight
	float mc_weight = evntuple.MC_weight;
	
	// Pileup
	float pu_weight = sfhelper.Get_PUWeight(evntuple.npuTrue);
	
	// HLT scale factor
	bool hlt1LTriggered =
		trighelper.pass_single_lep_triggers(evntuple.triggerBits);
	bool hltXTriggered =
		trighelper.pass_leptau_cross_triggers(evntuple.triggerBits);
	float hlt_sf = sfhelper.Get_HLTSF(leptons, taus, hlt1LTriggered, hltXTriggered);

	// lepton ID scale factor
	// NEED UPDATE: electron loose vs reco
	float lepid_sf = sfhelper.Get_LeptonIDSF_weight(leptons);

	// tau ID scale factor
	float tauid_sf = sfhelper.Get_TauIDSF_weight(taus);
	
	// btag scale factor
	if (JES != "JESUp" and JES != "JESDown") JES = "NA";
	float btag_sf = sfhelper.Get_EvtCSVWeight(jets,JES);

	////////////////////
	mvantuple.event_weight =
		pu_weight * mc_weight * btag_sf * lepid_sf * tauid_sf * hlt_sf;

	mvantuple.pu_weight = pu_weight;
	mvantuple.mc_weight = mc_weight;
	mvantuple.btag_sf = btag_sf;
	mvantuple.lepid_sf = lepid_sf;
	mvantuple.tauid_sf = tauid_sf;
	mvantuple.hlt_sf = hlt_sf;
	
	////////////////////
	if (syst) { // systematic uncertainties
		// theoretical uncertainties
		mvantuple.event_weight_thu_shape_x1Up =
			mvantuple.event_weight / mc_weight * evntuple.MC_weight_scale_muF2;
		mvantuple.event_weight_thu_shape_x1Down =
			mvantuple.event_weight / mc_weight * evntuple.MC_weight_scale_muF0p5;
		mvantuple.event_weight_thu_shape_y1Up =
			mvantuple.event_weight / mc_weight * evntuple.MC_weight_scale_muR2;
		mvantuple.event_weight_thu_shape_y1Down =
			mvantuple.event_weight / mc_weight * evntuple.MC_weight_scale_muR0p5;

		// trigger
		if (anatype==Analyze_2lss1tau or anatype==Analyze_2lss or
			anatype==Analyze_2l2tau) {
			mvantuple.event_weight_triggerUp = mvantuple.event_weight / hlt_sf * (sfhelper.Get_HLTSF_2l(leptons, "triggerUp"));
			mvantuple.event_weight_triggerDown = mvantuple.event_weight / hlt_sf * (sfhelper.Get_HLTSF_2l(leptons, "triggerDown"));
		}
		
		// b-tagging
		mvantuple.event_weight_btag_LFUp = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "LFUp"));
		mvantuple.event_weight_btag_LFDown = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "LFDown"));
		mvantuple.event_weight_btag_HFUp = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "HFUp"));
		mvantuple.event_weight_btag_HFDown = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "HFDown"));
		mvantuple.event_weight_btag_HFStats1Up = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "HFStats1Up"));
		mvantuple.event_weight_btag_HFStats1Down = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "HFStats1Down"));
		mvantuple.event_weight_btag_HFStats2Up = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "HFStats2Up"));
		mvantuple.event_weight_btag_HFStats2Down = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "HFStats2Down"));
		mvantuple.event_weight_btag_LFStats1Up = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "LFStats1Up"));
		mvantuple.event_weight_btag_LFStats1Down = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "LFStats1Down"));
		mvantuple.event_weight_btag_LFStats2Up = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "LFStats2Up"));
		mvantuple.event_weight_btag_LFStats2Down = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "LFStats2Down"));
		mvantuple.event_weight_btag_cErr1Up = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "cErr1Up"));
		mvantuple.event_weight_btag_cErr1Down = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "cErr1Down"));
		mvantuple.event_weight_btag_cErr2Up = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "cErr2Up"));
		mvantuple.event_weight_btag_cErr2Down = mvantuple.event_weight / btag_sf * (sfhelper.Get_EvtCSVWeight(jets, "cErr2Down"));

		// tau ID
		if (anatype==Analyze_2lss1tau or anatype==Analyze_3l1tau) {
			mvantuple.event_weight_FRjt_normUp = mvantuple.event_weight / tauid_sf * sfhelper.Get_TauIDSF_weight(taus, "FRjt_normUp");
			mvantuple.event_weight_FRjt_normDown = mvantuple.event_weight / tauid_sf * sfhelper.Get_TauIDSF_weight(taus, "FRjt_normDown");
			mvantuple.event_weight_FRjt_shapeUp = mvantuple.event_weight / tauid_sf * sfhelper.Get_TauIDSF_weight(taus, "FRjt_shapeUp");
			mvantuple.event_weight_FRjt_shapeDown = mvantuple.event_weight / tauid_sf * sfhelper.Get_TauIDSF_weight(taus, "FRjt_shapeDown");
		}
		
		// lepton ID
		mvantuple.event_weight_lepEff_mulooseUp = mvantuple.event_weight / lepid_sf * sfhelper.Get_LeptonIDSF_weight(leptons, "lepEff_mulooseUp");
		mvantuple.event_weight_lepEff_mulooseDown = mvantuple.event_weight / lepid_sf * sfhelper.Get_LeptonIDSF_weight(leptons, "lepEff_mulooseDown");
		mvantuple.event_weight_lepEff_ellooseUp = mvantuple.event_weight / lepid_sf * sfhelper.Get_LeptonIDSF_weight(leptons, "lepEff_ellooseUp");
		mvantuple.event_weight_lepEff_ellooseDown = mvantuple.event_weight / lepid_sf * sfhelper.Get_LeptonIDSF_weight(leptons, "lepEff_ellooseDown");
		mvantuple.event_weight_lepEff_mutightUp = mvantuple.event_weight / lepid_sf * sfhelper.Get_LeptonIDSF_weight(leptons, "lepEff_mutightUp");
		mvantuple.event_weight_lepEff_mutightDown = mvantuple.event_weight / lepid_sf * sfhelper.Get_LeptonIDSF_weight(leptons, "lepEff_mutightDown");
		mvantuple.event_weight_lepEff_eltightUp = mvantuple.event_weight / lepid_sf * sfhelper.Get_LeptonIDSF_weight(leptons, "lepEff_eltightUp");
		mvantuple.event_weight_lepEff_eltightDown = mvantuple.event_weight / lepid_sf * sfhelper.Get_LeptonIDSF_weight(leptons, "lepEff_eltightDown");
	}
}
