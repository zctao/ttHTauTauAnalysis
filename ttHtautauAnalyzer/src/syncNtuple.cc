#include "../interface/syncNtuple.h"

void syncNtuple::initialize()
{
    // event variables
    nEvent = -9999;
    ls = -9999;
    run = -9999;

    n_presel_mu = -9999;
    n_mvasel_mu = -9999;
    n_fakeablesel_mu = -9999;
    n_presel_ele = -9999;
    n_mvasel_ele = -9999;
    n_fakeablesel_ele = -9999;
    n_presel_tau = -9999;
    // n_tau = -9999;
    n_presel_jet = -9999;

    // event_weight = -9999.;
    PU_weight = -9999.;
    MC_weight = -9999.;
    bTagSF_weight = -9999.;
    leptonSF_weight = -9999.;
    tauSF_weight = -9999.;
    triggerSF_weight = -9999.;
    FR_weight = -9999.;

    // muons
    mu0_pt = -9999.;
    mu0_conept = -9999.;
    mu0_eta = -9999.;
    mu0_phi = -9999.;
    mu0_E = -9999.;
    mu0_charge = -9999;
    mu0_jetNDauChargedMVASel = -9999;
    mu0_miniRelIso = -9999.;
    mu0_miniIsoCharged = -9999.;
    mu0_miniIsoNeutral = -9999.;
    mu0_jetPtRel = -9999.;
    mu0_jetPtRatio = -9999.;
    mu0_jetCSV = -9999.;
    mu0_sip3D = -9999.;
    mu0_dxy = -9999.;
    mu0_dz = -9999.;
    mu0_segmentCompatibility = -9999.;
    mu0_leptonMVA = -9999.;
    mu0_mediumID = -9999.;
    mu0_dpt_div_pt = -9999.;
    mu0_ismvasel = 0;      //-9999;
    mu0_isfakeablesel = 0; //-9999;
    // mu0_mcMatchType = -9999;
    // mu0_isPFMuon = -9999;
    mu1_pt = -9999.;
    mu1_conept = -9999.;
    mu1_eta = -9999.;
    mu1_phi = -9999.;
    mu1_E = -9999.;
    mu1_charge = -9999;
    mu1_jetNDauChargedMVASel = -9999;
    mu1_miniRelIso = -9999.;
    mu1_miniIsoCharged = -9999.;
    mu1_miniIsoNeutral = -9999.;
    mu1_jetPtRel = -9999.;
    mu1_jetPtRatio = -9999.;
    mu1_jetCSV = -9999.;
    mu1_sip3D = -9999.;
    mu1_dxy = -9999.;
    mu1_dz = -9999.;
    mu1_segmentCompatibility = -9999.;
    mu1_leptonMVA = -9999.;
    mu1_mediumID = -9999.;
    mu1_dpt_div_pt = -9999.;
    mu1_ismvasel = 0;      //-9999;
    mu1_isfakeablesel = 0; //-9999;
    // mu1_mcMatchType = -9999;
    // mu1_isPFMuon = -9999;

    // electrons
    ele0_pt = -9999.;
    ele0_conept = -9999.;
    ele0_eta = -9999.;
    ele0_phi = -9999.;
    ele0_E = -9999.;
    ele0_charge = -9999;
    ele0_jetNDauChargedMVASel = -9999;
    ele0_miniRelIso = -9999.;
    ele0_miniIsoCharged = -9999.;
    ele0_miniIsoNeutral = -9999.;
    ele0_jetPtRel = -9999.;
    ele0_jetPtRatio = -9999.;
    ele0_jetCSV = -9999.;
    ele0_sip3D = -9999.;
    ele0_dxy = -9999.;
    ele0_dz = -9999.;
    ele0_ntMVAeleID = -9999.;
    ele0_leptonMVA = -9999.;
    ele0_isChargeConsistent = 0;   //-9999;
    ele0_passesConversionVeto = 0; //-9999;
    ele0_nMissingHits = -9999;
    ele0_ismvasel = 0;      //-9999;
    ele0_isfakeablesel = 0; //-9999;
    // ele0_mcMatchType = -9999;
    ele1_pt = -9999.;
    ele1_conept = -9999.;
    ele1_eta = -9999.;
    ele1_phi = -9999.;
    ele1_E = -9999.;
    ele1_charge = -9999;
    ele1_jetNDauChargedMVASel = -9999;
    ele1_miniRelIso = -9999.;
    ele1_miniIsoCharged = -9999.;
    ele1_miniIsoNeutral = -9999.;
    ele1_jetPtRel = -9999.;
    ele1_jetPtRatio = -9999.;
    ele1_jetCSV = -9999.;
    ele1_sip3D = -9999.;
    ele1_dxy = -9999.;
    ele1_dz = -9999.;
    ele1_ntMVAeleID = -9999.;
    ele1_leptonMVA = -9999.;
    ele1_isChargeConsistent = 0;   //-9999;
    ele1_passesConversionVeto = 0; //-9999;
    ele1_nMissingHits = -9999;
    ele1_ismvasel = 0;      //-9999;
    ele1_isfakeablesel = 0; //-9999;
    // ele1_mcMatchType = -9999;

    // taus
    tau0_pt = -9999.;
    tau0_eta = -9999.;
    tau0_phi = -9999.;
    tau0_E = -9999.;
    tau0_charge = -9999;
    tau0_dxy = -9999.;
    tau0_dz = -9999.;
    // tau0_decayMode = -9999;
    tau0_decayModeFindingOldDMs = -9999;
    tau0_decayModeFindingNewDMs = -9999;
    tau0_byCombinedIsolationDeltaBetaCorr3Hits = -9999.;
    tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits = -9999;
    tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits = -9999;
    tau0_byTightCombinedIsolationDeltaBetaCorr3Hits = -9999;
    tau0_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
    tau0_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
    tau0_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
    tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT = -9999;
    tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT = -9999;
    tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
    tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
    tau0_againstMuonLoose3 = -9999;
    tau0_againstMuonTight3 = -9999;
    tau0_againstElectronVLooseMVA6 = -9999;
    tau0_againstElectronLooseMVA6 = -9999;
    tau0_againstElectronMediumMVA6 = -9999;
    tau0_againstElectronTightMVA6 = -9999;
    // tau0_mcMatchType = -9999;
    tau1_pt = -9999.;
    tau1_eta = -9999.;
    tau1_phi = -9999.;
    tau1_E = -9999.;
    tau1_charge = -9999;
    tau1_dxy = -9999.;
    tau1_dz = -9999.;
    // tau1_decayMode = -9999;
    tau1_decayModeFindingOldDMs = -9999;
    tau1_decayModeFindingNewDMs = -9999;
    tau1_byCombinedIsolationDeltaBetaCorr3Hits = -9999.;
    tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits = -9999;
    tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits = -9999;
    tau1_byTightCombinedIsolationDeltaBetaCorr3Hits = -9999;
    tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
    tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
    tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03 = -9999;
    tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT = -9999;
    tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT = -9999;
    tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
    tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT = -9999;
    tau1_againstMuonLoose3 = -9999;
    tau1_againstMuonTight3 = -9999;
    tau1_againstElectronVLooseMVA6 = -9999;
    tau1_againstElectronLooseMVA6 = -9999;
    tau1_againstElectronMediumMVA6 = -9999;
    tau1_againstElectronTightMVA6 = -9999;
    // tau1_mcMatchType = -9999;

    // jets
    jet0_pt = -9999.;
    jet0_eta = -9999.;
    jet0_phi = -9999.;
    jet0_E = -9999.;
    jet0_CSV = -9999.;
    jet1_pt = -9999.;
    jet1_eta = -9999.;
    jet1_phi = -9999.;
    jet1_E = -9999.;
    jet1_CSV = -9999.;
    jet2_pt = -9999.;
    jet2_eta = -9999.;
    jet2_phi = -9999.;
    jet2_E = -9999.;
    jet2_CSV = -9999.;
    jet3_pt = -9999.;
    jet3_eta = -9999.;
    jet3_phi = -9999.;
    jet3_E = -9999.;
    jet3_CSV = -9999.;

    // MET
    PFMET = -9999.;
    PFMETphi = -9999.;
    MHT = -9999.;
    metLD = -9999.;
    // METCov00 = -9999.;
    // METCov01 = -9999.;
    // METCov10 = -9999.;
    // METCov11 = -9999.;

    // event-level MVA variables
    isGenMatched = -9999;
    lep0_conept = -9999.;
    lep1_conept = -9999.;
    mindr_lep0_jet = -9999.;
    mindr_lep1_jet = -9999.;
    mindr_lep2_jet = -9999.;
    mindr_tau_jet = -9999.;
    MT_met_lep0 = -9999.;
    MT_met_lep2 = -9999.;
    avg_dr_jet = -9999.;

    dr_leps = -9999.;
    mvis_lep0_tau = -9999.;
    mvis_lep1_tau = -9999.;
    max_lep_eta = -9999.;
    dr_lep0_tau = -9999.;

    MVA_2lss_ttV = -9999.;
    MVA_2lss_ttbar = -9999.;
    tt_deltaR = -9999.;
    ntags = -9999;
    ntags_loose = -9999;
    tt_mvis = -9999.;
    tt_pt = -9999.;
    max_dr_jet = -9999.;
    HT = -9999.;
    MVA_1l2tau_ttbar = -9999.;
    MVA_1l2tau_ttbar_v2 = -9999.;
    MVA_1l2tau_ttZ_v2 = -9999.;
    MVA_1l2tau_2Dbin_v2 = -9999;
    mvis_l1tau = -9999.;
    dR_l0tau = -9999.;
    dR_l1tau = -9999.;
    dR_l2tau = -9999.;
    mT_lep2 = -9999.;
    MVA_3l1tau_ttbar = -9999.;
    MVA_3l1tau_ttV = -9999.;
    MVA_3l1tau_2Dbin = -9999;

    // MEM variables
    Integral_ttH = -9999.;
    Integral_ttZ = -9999.;
    Integral_ttZ_Zll = -9999.;
    Integral_ttbar = -9999.;
    Integration_type = -9999;
    MEM_LR = -9999.;
    dR_leps = -9999.;
    mvis_l0tau = -9999.;
    mT_lep0 = -9999.;
    MVA_2lSS1tau_noMEM_ttbar = -9999.;
    MVA_2lSS1tau_noMEM_ttV = -9999.;
    MVA_2lSS1tau_noMEM_2Dbin = -9999;
    MVA_2lSS1tau_MEM_ttbar = -9999.;
    MVA_2lSS1tau_MEM_ttV = -9999.;
    MVA_2lSS1tau_MEM_2Dbin = -9999;

    /*
    MC_weight_scale_muF0p5 = -9999.;
    MC_weight_scale_muF2 = -9999.;
    MC_weight_scale_muR0p5 = -9999.;
    MC_weight_scale_muR2 = -9999.;
    btagSF_weight_LFUp = -9999.;
    btagSF_weight_LFDown = -9999.;
    btagSF_weight_HFUp = -9999.;
    btagSF_weight_HFDown = -9999.;
    btagSF_weight_HFStats1Up = -9999.;
    btagSF_weight_HFStats1Down = -9999.;
    btagSF_weight_HFStats2Up = -9999.;
    btagSF_weight_HFStats2Down = -9999.;
    btagSF_weight_LFStats1Up = -9999.;
    btagSF_weight_LFStats1Down = -9999.;
    btagSF_weight_LFStats2Up = -9999.;
    btagSF_weight_LFStats2Down = -9999.;
    btagSF_weight_cErr1Up = -9999.;
    btagSF_weight_cErr1Down = -9999.;
    btagSF_weight_cErr2Up = -9999.;
    btagSF_weight_cErr2Down = -9999.;

    HiggsDecayType = -9999;
    lepCategory = -9999;
    btagCategory = -9999;
    npuTrue = -9999.;
    npuInTime = -9999.;

    pass_single_mu = -9999;
    pass_single_e = -9999;
    pass_double_mu = -9999;
    pass_double_e = -9999;
    pass_elemu = -9999;
    matchHLTPath = -9999;
    triggerBits = 0;
    filterBits = 0;

    nBadMuons = -9999;

    ibin = -9999;
    lepXtauCharge = -9999;

    n_jet25_recl = -9999;

    max_lep_eta = -9999.;
    */
}

void syncNtuple::set_up_branches(TTree *tree)
{
    // initialize ntuple
    this->initialize();

    // set up tree branches
    tree->Branch("nEvent", &nEvent);
    tree->Branch("ls", &ls);
    tree->Branch("run", &run);
    tree->Branch("n_presel_mu", &n_presel_mu);
    tree->Branch("n_mvasel_mu", &n_mvasel_mu);
    tree->Branch("n_fakeablesel_mu", &n_fakeablesel_mu);
    tree->Branch("n_presel_ele", &n_presel_ele);
    tree->Branch("n_mvasel_ele", &n_mvasel_ele);
    tree->Branch("n_fakeablesel_ele", &n_fakeablesel_ele);
    tree->Branch("n_presel_tau", &n_presel_tau);
    // tree->Branch("n_tau", &n_tau);
    tree->Branch("n_presel_jet", &n_presel_jet);
    // muons
    tree->Branch("mu0_pt", &mu0_pt);
    tree->Branch("mu0_conept", &mu0_conept);
    tree->Branch("mu0_eta", &mu0_eta);
    tree->Branch("mu0_phi", &mu0_phi);
    tree->Branch("mu0_E", &mu0_E);
    tree->Branch("mu0_charge", &mu0_charge);
    tree->Branch("mu0_jetNDauChargedMVASel", &mu0_jetNDauChargedMVASel);
    tree->Branch("mu0_miniRelIso", &mu0_miniRelIso);
    tree->Branch("mu0_miniIsoCharged", &mu0_miniIsoCharged);
    tree->Branch("mu0_miniIsoNeutral", &mu0_miniIsoNeutral);
    tree->Branch("mu0_jetPtRel", &mu0_jetPtRel);
    tree->Branch("mu0_jetPtRatio", &mu0_jetPtRatio);
    tree->Branch("mu0_jetCSV", &mu0_jetCSV);
    tree->Branch("mu0_sip3D", &mu0_sip3D);
    tree->Branch("mu0_dxy", &mu0_dxy);
    tree->Branch("mu0_dz", &mu0_dz);
    tree->Branch("mu0_segmentCompatibility", &mu0_segmentCompatibility);
    tree->Branch("mu0_leptonMVA", &mu0_leptonMVA);
    tree->Branch("mu0_mediumID", &mu0_mediumID);
    tree->Branch("mu0_dpt_div_pt", &mu0_dpt_div_pt);
    tree->Branch("mu0_ismvasel", &mu0_ismvasel);
    tree->Branch("mu0_isfakeablesel", &mu0_isfakeablesel);
    // tree->Branch("mu0_mcMatchType",          &mu0_mcMatchType);
    // tree->Branch("mu0_isPFMuon",             &mu0_isPFMuon);
    tree->Branch("mu1_pt", &mu1_pt);
    tree->Branch("mu1_conept", &mu1_conept);
    tree->Branch("mu1_eta", &mu1_eta);
    tree->Branch("mu1_phi", &mu1_phi);
    tree->Branch("mu1_E", &mu1_E);
    tree->Branch("mu1_charge", &mu1_charge);
    tree->Branch("mu1_jetNDauChargedMVASel", &mu1_jetNDauChargedMVASel);
    tree->Branch("mu1_miniRelIso", &mu1_miniRelIso);
    tree->Branch("mu1_miniIsoCharged", &mu1_miniIsoCharged);
    tree->Branch("mu1_miniIsoNeutral", &mu1_miniIsoNeutral);
    tree->Branch("mu1_jetPtRel", &mu1_jetPtRel);
    tree->Branch("mu1_jetPtRatio", &mu1_jetPtRatio);
    tree->Branch("mu1_jetCSV", &mu1_jetCSV);
    tree->Branch("mu1_sip3D", &mu1_sip3D);
    tree->Branch("mu1_dxy", &mu1_dxy);
    tree->Branch("mu1_dz", &mu1_dz);
    tree->Branch("mu1_segmentCompatibility", &mu1_segmentCompatibility);
    tree->Branch("mu1_leptonMVA", &mu1_leptonMVA);
    tree->Branch("mu1_mediumID", &mu1_mediumID);
    tree->Branch("mu1_dpt_div_pt", &mu1_dpt_div_pt);
    tree->Branch("mu1_ismvasel", &mu1_ismvasel);
    tree->Branch("mu1_isfakeablesel", &mu1_isfakeablesel);
    // tree->Branch("mu1_mcMatchType",          &mu1_mcMatchType);
    // tree->Branch("mu1_isPFMuon",             &mu1_isPFMuon);
    // electrons
    tree->Branch("ele0_pt", &ele0_pt);
    tree->Branch("ele0_conept", &ele0_conept);
    tree->Branch("ele0_eta", &ele0_eta);
    tree->Branch("ele0_phi", &ele0_phi);
    tree->Branch("ele0_E", &ele0_E);
    tree->Branch("ele0_charge", &ele0_charge);
    tree->Branch("ele0_jetNDauChargedMVASel", &ele0_jetNDauChargedMVASel);
    tree->Branch("ele0_miniRelIso", &ele0_miniRelIso);
    tree->Branch("ele0_miniIsoCharged", &ele0_miniIsoCharged);
    tree->Branch("ele0_miniIsoNeutral", &ele0_miniIsoNeutral);
    tree->Branch("ele0_jetPtRel", &ele0_jetPtRel);
    tree->Branch("ele0_jetPtRatio", &ele0_jetPtRatio);
    tree->Branch("ele0_jetCSV", &ele0_jetCSV);
    tree->Branch("ele0_sip3D", &ele0_sip3D);
    tree->Branch("ele0_dxy", &ele0_dxy);
    tree->Branch("ele0_dz", &ele0_dz);
    tree->Branch("ele0_ntMVAeleID", &ele0_ntMVAeleID);
    tree->Branch("ele0_leptonMVA", &ele0_leptonMVA);
    tree->Branch("ele0_isChargeConsistent", &ele0_isChargeConsistent);
    tree->Branch("ele0_passesConversionVeto", &ele0_passesConversionVeto);
    tree->Branch("ele0_nMissingHits", &ele0_nMissingHits);
    tree->Branch("ele0_ismvasel", &ele0_ismvasel);
    tree->Branch("ele0_isfakeablesel", &ele0_isfakeablesel);
    // tree->Branch("ele0_mcMatchType", &ele0_mcMatchType);
    tree->Branch("ele1_pt", &ele1_pt);
    tree->Branch("ele1_conept", &ele1_conept);
    tree->Branch("ele1_eta", &ele1_eta);
    tree->Branch("ele1_phi", &ele1_phi);
    tree->Branch("ele1_E", &ele1_E);
    tree->Branch("ele1_charge", &ele1_charge);
    tree->Branch("ele1_jetNDauChargedMVASel", &ele1_jetNDauChargedMVASel);
    tree->Branch("ele1_miniRelIso", &ele1_miniRelIso);
    tree->Branch("ele1_miniIsoCharged", &ele1_miniIsoCharged);
    tree->Branch("ele1_miniIsoNeutral", &ele1_miniIsoNeutral);
    tree->Branch("ele1_jetPtRel", &ele1_jetPtRel);
    tree->Branch("ele1_jetPtRatio", &ele1_jetPtRatio);
    tree->Branch("ele1_jetCSV", &ele1_jetCSV);
    tree->Branch("ele1_sip3D", &ele1_sip3D);
    tree->Branch("ele1_dxy", &ele1_dxy);
    tree->Branch("ele1_dz", &ele1_dz);
    tree->Branch("ele1_ntMVAeleID", &ele1_ntMVAeleID);
    tree->Branch("ele1_leptonMVA", &ele1_leptonMVA);
    tree->Branch("ele1_isChargeConsistent", &ele1_isChargeConsistent);
    tree->Branch("ele1_passesConversionVeto", &ele1_passesConversionVeto);
    tree->Branch("ele1_nMissingHits", &ele1_nMissingHits);
    tree->Branch("ele1_ismvasel", &ele1_ismvasel);
    tree->Branch("ele1_isfakeablesel", &ele1_isfakeablesel);
    // tree->Branch("ele1_mcMatchType", &ele1_mcMatchType);
    // taus
    tree->Branch("tau0_pt", &tau0_pt);
    tree->Branch("tau0_eta", &tau0_eta);
    tree->Branch("tau0_phi", &tau0_phi);
    tree->Branch("tau0_E", &tau0_E);
    tree->Branch("tau0_charge", &tau0_charge);
    tree->Branch("tau0_dxy", &tau0_dxy);
    tree->Branch("tau0_dz", &tau0_dz);
    // tree->Branch("tau0_decayMode", &tau0_decayMode);
    tree->Branch("tau0_decayModeFindingOldDMs", &tau0_decayModeFindingOldDMs);
    tree->Branch("tau0_decayModeFindingNewDMs", &tau0_decayModeFindingNewDMs);
    tree->Branch("tau0_byCombinedIsolationDeltaBetaCorr3Hits",
                 &tau0_byCombinedIsolationDeltaBetaCorr3Hits);
    tree->Branch("tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits",
                 &tau0_byLooseCombinedIsolationDeltaBetaCorr3Hits);
    tree->Branch("tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits",
                 &tau0_byMediumCombinedIsolationDeltaBetaCorr3Hits);
    tree->Branch("tau0_byTightCombinedIsolationDeltaBetaCorr3Hits",
                 &tau0_byTightCombinedIsolationDeltaBetaCorr3Hits);
    tree->Branch("tau0_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",
                 &tau0_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03);
    tree->Branch("tau0_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",
                 &tau0_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03);
    tree->Branch("tau0_byTightCombinedIsolationDeltaBetaCorr3HitsdR03",
                 &tau0_byTightCombinedIsolationDeltaBetaCorr3HitsdR03);
    tree->Branch("tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT",
                 &tau0_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
    tree->Branch("tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT",
                 &tau0_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
    tree->Branch("tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT",
                 &tau0_byTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree->Branch("tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT",
                 &tau0_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree->Branch("tau0_againstMuonLoose3", &tau0_againstMuonLoose3);
    tree->Branch("tau0_againstMuonTight3", &tau0_againstMuonTight3);
    tree->Branch("tau0_againstElectronVLooseMVA6",
                 &tau0_againstElectronVLooseMVA6);
    tree->Branch("tau0_againstElectronLooseMVA6",
                 &tau0_againstElectronLooseMVA6);
    tree->Branch("tau0_againstElectronMediumMVA6",
                 &tau0_againstElectronMediumMVA6);
    tree->Branch("tau0_againstElectronTightMVA6",
                 &tau0_againstElectronTightMVA6);
    // tree->Branch("tau0_mcMatchType", &tau0_mcMatchType);
    tree->Branch("tau1_pt", &tau1_pt);
    tree->Branch("tau1_eta", &tau1_eta);
    tree->Branch("tau1_phi", &tau1_phi);
    tree->Branch("tau1_E", &tau1_E);
    tree->Branch("tau1_charge", &tau1_charge);
    tree->Branch("tau1_dxy", &tau1_dxy);
    tree->Branch("tau1_dz", &tau1_dz);
    // tree->Branch("tau1_decayMode", &tau1_decayMode);
    tree->Branch("tau1_decayModeFindingOldDMs", &tau1_decayModeFindingOldDMs);
    tree->Branch("tau1_decayModeFindingNewDMs", &tau1_decayModeFindingNewDMs);
    tree->Branch("tau1_byCombinedIsolationDeltaBetaCorr3Hits",
                 &tau1_byCombinedIsolationDeltaBetaCorr3Hits);
    tree->Branch("tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits",
                 &tau1_byLooseCombinedIsolationDeltaBetaCorr3Hits);
    tree->Branch("tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits",
                 &tau1_byMediumCombinedIsolationDeltaBetaCorr3Hits);
    tree->Branch("tau1_byTightCombinedIsolationDeltaBetaCorr3Hits",
                 &tau1_byTightCombinedIsolationDeltaBetaCorr3Hits);
    tree->Branch("tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",
                 &tau1_byLooseCombinedIsolationDeltaBetaCorr3HitsdR03);
    tree->Branch("tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",
                 &tau1_byMediumCombinedIsolationDeltaBetaCorr3HitsdR03);
    tree->Branch("tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03",
                 &tau1_byTightCombinedIsolationDeltaBetaCorr3HitsdR03);
    tree->Branch("tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT",
                 &tau1_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
    tree->Branch("tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT",
                 &tau1_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
    tree->Branch("tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT",
                 &tau1_byTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree->Branch("tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT",
                 &tau1_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree->Branch("tau1_againstMuonLoose3", &tau1_againstMuonLoose3);
    tree->Branch("tau1_againstMuonTight3", &tau1_againstMuonTight3);
    tree->Branch("tau1_againstElectronVLooseMVA6",
                 &tau1_againstElectronVLooseMVA6);
    tree->Branch("tau1_againstElectronLooseMVA6",
                 &tau1_againstElectronLooseMVA6);
    tree->Branch("tau1_againstElectronMediumMVA6",
                 &tau1_againstElectronMediumMVA6);
    tree->Branch("tau1_againstElectronTightMVA6",
                 &tau1_againstElectronTightMVA6);
    // tree->Branch("tau1_mcMatchType", &tau1_mcMatchType);
    // jets
    tree->Branch("jet0_pt", &jet0_pt);
    tree->Branch("jet0_eta", &jet0_eta);
    tree->Branch("jet0_phi", &jet0_phi);
    tree->Branch("jet0_E", &jet0_E);
    tree->Branch("jet0_CSV", &jet0_CSV);
    tree->Branch("jet1_pt", &jet1_pt);
    tree->Branch("jet1_eta", &jet1_eta);
    tree->Branch("jet1_phi", &jet1_phi);
    tree->Branch("jet1_E", &jet1_E);
    tree->Branch("jet1_CSV", &jet1_CSV);
    tree->Branch("jet2_pt", &jet2_pt);
    tree->Branch("jet2_eta", &jet2_eta);
    tree->Branch("jet2_phi", &jet2_phi);
    tree->Branch("jet2_E", &jet2_E);
    tree->Branch("jet2_CSV", &jet2_CSV);
    tree->Branch("jet3_pt", &jet3_pt);
    tree->Branch("jet3_eta", &jet3_eta);
    tree->Branch("jet3_phi", &jet3_phi);
    tree->Branch("jet3_E", &jet3_E);
    tree->Branch("jet3_CSV", &jet3_CSV);
    // MET
    tree->Branch("PFMET", &PFMET);
    tree->Branch("PFMETphi", &PFMETphi);
    tree->Branch("MHT", &MHT);
    tree->Branch("metLD", &metLD);
    // tree->Branch("METSignificance", &METSignificance);
    // tree->Branch("METCov00", &METCov00);
    // tree->Branch("METCov01", &METCov01);
    // tree->Branch("METCov10", &METCov10);
    // tree->Branch("METCov11", &METCov11);
    // event weights
    // tree->Branch("event_weight", &event_weight);
    tree->Branch("PU_weight", &PU_weight);
    tree->Branch("MC_weight", &MC_weight);
    tree->Branch("bTagSF_weight", &bTagSF_weight);
    tree->Branch("leptonSF_weight", &leptonSF_weight);
    tree->Branch("tauSF_weight", &tauSF_weight);
    tree->Branch("triggerSF_weight", &triggerSF_weight);
    tree->Branch("FR_weight", &FR_weight);
    // additional event-level MVA variables
    tree->Branch("isGenMatched", &isGenMatched);
    tree->Branch("lep0_conept", &lep0_conept);
    tree->Branch("lep1_conept", &lep1_conept);
    tree->Branch("mindr_lep0_jet", &mindr_lep0_jet);
    tree->Branch("mindr_lep1_jet", &mindr_lep1_jet);
    tree->Branch("mindr_lep2_jet", &mindr_lep2_jet);
    tree->Branch("mindr_tau_jet", &mindr_tau_jet);
    tree->Branch("MT_met_lep0", &MT_met_lep0);
    tree->Branch("MT_met_lep2", &MT_met_lep2);
    tree->Branch("avg_dr_jet", &avg_dr_jet);

    tree->Branch("dr_leps", &dr_leps);
    tree->Branch("mvis_lep0_tau", &mvis_lep0_tau);
    tree->Branch("mvis_lep1_tau", &mvis_lep1_tau);
    tree->Branch("max_lep_eta", &max_lep_eta);
    tree->Branch("dr_lep0_tau", &dr_lep0_tau);

    tree->Branch("MVA_2lss_ttV", &MVA_2lss_ttV);
    tree->Branch("MVA_2lss_ttbar", &MVA_2lss_ttbar);
    tree->Branch("tt_deltaR", &tt_deltaR);
    tree->Branch("ntags", &ntags);
    tree->Branch("ntags_loose", &ntags_loose);
    tree->Branch("tt_mvis", &tt_mvis);
    tree->Branch("tt_pt", &tt_pt);
    tree->Branch("max_dr_jet", &max_dr_jet);
    tree->Branch("HT", &HT);
    tree->Branch("MVA_1l2tau_ttbar", &MVA_1l2tau_ttbar);
    tree->Branch("MVA_1l2tau_ttbar_v2", &MVA_1l2tau_ttbar_v2);
    tree->Branch("MVA_1l2tau_ttZ_v2", &MVA_1l2tau_ttZ_v2);
    tree->Branch("MVA_1l2tau_2Dbin_v2", &MVA_1l2tau_2Dbin_v2);
    tree->Branch("mvis_l1tau", &mvis_l1tau);
    tree->Branch("dR_l0tau", &dR_l0tau);
    tree->Branch("dR_l1tau", &dR_l1tau);
    tree->Branch("dR_l2tau", &dR_l2tau);
    // tree->Branch("mT_lep2", &mT_lep2);
    tree->Branch("MVA_3l1tau_ttbar", &MVA_3l1tau_ttbar);
    tree->Branch("MVA_3l1tau_ttV", &MVA_3l1tau_ttV);
    tree->Branch("MVA_3l1tau_2Dbin", &MVA_3l1tau_2Dbin);
    tree->Branch("Integral_ttH", &Integral_ttH);
    tree->Branch("Integral_ttZ", &Integral_ttZ);
    tree->Branch("Integral_ttZ_Zll", &Integral_ttZ_Zll);
    tree->Branch("Integral_ttbar", &Integral_ttbar);
    tree->Branch("Integration_type", &Integration_type);
    tree->Branch("MEM_LR", &MEM_LR);
    tree->Branch("dR_leps", &dR_leps);
    tree->Branch("mvis_l0tau", &mvis_l0tau);
    // tree->Branch("mT_lep0", &mT_lep0);
    tree->Branch("MVA_2lSS1tau_noMEM_ttbar", &MVA_2lSS1tau_noMEM_ttbar);
    tree->Branch("MVA_2lSS1tau_noMEM_ttV", &MVA_2lSS1tau_noMEM_ttV);
    tree->Branch("MVA_2lSS1tau_noMEM_2Dbin", &MVA_2lSS1tau_noMEM_2Dbin);
    tree->Branch("MVA_2lSS1tau_MEM_ttbar", &MVA_2lSS1tau_MEM_ttbar);
    tree->Branch("MVA_2lSS1tau_MEM_ttV", &MVA_2lSS1tau_MEM_ttV);
    tree->Branch("MVA_2lSS1tau_MEM_2Dbin", &MVA_2lSS1tau_MEM_2Dbin);

    /*
    tree->Branch("MC_weight_scale_muF0p5", &MC_weight_scale_muF0p5);
    tree->Branch("MC_weight_scale_muF2", &MC_weight_scale_muF2);
    tree->Branch("MC_weight_scale_muR0p5", &MC_weight_scale_muR0p5);
    tree->Branch("MC_weight_scale_muR2", &MC_weight_scale_muR2);
    tree->Branch("btagSF_weight_LFUp", &btagSF_weight_LFUp);
    tree->Branch("btagSF_weight_LFDown", &btagSF_weight_LFDown);
    tree->Branch("btagSF_weight_HFUp", &btagSF_weight_HFUp);
    tree->Branch("btagSF_weight_HFDown", &btagSF_weight_HFDown);
    tree->Branch("btagSF_weight_HFStats1Up", &btagSF_weight_HFStats1Up);
    tree->Branch("btagSF_weight_HFStats1Down", &btagSF_weight_HFStats1Down);
    tree->Branch("btagSF_weight_HFStats2Up", &btagSF_weight_HFStats2Up);
    tree->Branch("btagSF_weight_HFStats2Down", &btagSF_weight_HFStats2Down);
    tree->Branch("btagSF_weight_LFStats1Up", &btagSF_weight_LFStats1Up);
    tree->Branch("btagSF_weight_LFStats1Down", &btagSF_weight_LFStats1Down);
    tree->Branch("btagSF_weight_LFStats2Up", &btagSF_weight_LFStats2Up);
    tree->Branch("btagSF_weight_LFStats2Down", &btagSF_weight_LFStats2Down);
    tree->Branch("btagSF_weight_cErr1Up", &btagSF_weight_cErr1Up);
    tree->Branch("btagSF_weight_cErr1Down", &btagSF_weight_cErr1Down);
    tree->Branch("btagSF_weight_cErr2Up", &btagSF_weight_cErr2Up);
    tree->Branch("btagSF_weight_cErr2Down", &btagSF_weight_cErr2Down);
    tree->Branch("HiggsDecayType", &HiggsDecayType);
    tree->Branch("lepCategory", &lepCategory);
    tree->Branch("btagCategory", &btagCategory);
    tree->Branch("npuTrue", &npuTrue);
    tree->Branch("npuInTime", &npuInTime);
    tree->Branch("pass_single_e", &pass_single_e);
    tree->Branch("pass_single_mu", &pass_single_mu);
    tree->Branch("pass_double_e", &pass_double_e);
    tree->Branch("pass_double_mu", &pass_double_mu);
    tree->Branch("pass_elemu", &pass_elemu);
    tree->Branch("matchHLTPath", &matchHLTPath);
    tree->Branch("triggerBits", &triggerBits);
    tree->Branch("nBadMuons", &nBadMuons);
    tree->Branch("filterBits", &filterBits);
    tree->Branch("ibin", &ibin);
    tree->Branch("lepXtauCharge", &lepXtauCharge);
    */
}
