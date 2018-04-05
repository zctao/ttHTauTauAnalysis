#ifndef EventSelector_h
#define EventSelector_h

#if !defined(__ACLIC__) && !defined(__ROOTCLING__)
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
#endif

#include "TriggerHelper.h"
#include "Types_enum.h"
#include "miniLepton.h"
#include "miniTau.h"

class EventSelector
{
  public:
    // constructor and destructor
    EventSelector(Analysis_types anatype, Selection_types seltype, bool verbose,
                  bool isMC = true)
    {
        anaType_ = anatype;
        selType_ = seltype;
        verbose_ = verbose;
        isMC_ = isMC;
    }

    EventSelector(bool verbose, bool isMC, bool looseSelection = false)
    {
        verbose_ = verbose;
        isMC_ = isMC;
        looseselection_ = looseSelection;
        anaType_ = Analyze_NA;
        selType_ = Selection_NA;
    }

    ~EventSelector(){};

    // member functions
    void fill_cutflow(TH1 *, int ibin, const char *);
    bool pass_hlt_paths(Analysis_types, TriggerHelper *const, unsigned int);
    // bool pass_extra_event_selection(Analysis_types, Selection_types,
    //								const std::vector<miniLepton>&,
    //								const std::vector<miniTau>&);
    bool pass_extra_event_selection(
        Analysis_types, Selection_types, std::vector<miniLepton> const *const,
        std::vector<miniTau> const *const,
        std::vector<miniTau> const *const fakeabletau = 0);

    bool pass_1l2tau_inclusive_selection(const std::vector<miniLepton> &,
                                         const std::vector<miniLepton> &,
                                         const std::vector<miniLepton> &,
                                         const std::vector<miniTau> &, int, int,
                                         int, TH1 *h_cutflow = 0);
    bool pass_1l2tau_SR_selection(const std::vector<miniLepton> &,
                                  const std::vector<miniTau> &);
    bool pass_1l2tau_FakeAR_selection(const std::vector<miniLepton> &,
                                      const std::vector<miniTau> &,
                                      const std::vector<miniTau> &);
    bool pass_1l2tau_CR_selection(const std::vector<miniLepton> &,
                                  const std::vector<miniTau> &);
    bool pass_1l2tau_FakeARCR_selection(const std::vector<miniLepton> &,
                                        const std::vector<miniTau> &,
                                        const std::vector<miniTau> &);
    bool pass_1l2tau_tightID(const std::vector<miniLepton> &,
                             const std::vector<miniTau> &);
    bool pass_1l2tau_charge(const std::vector<miniTau> &);

    bool pass_2l_generic_selection(const std::vector<miniLepton> &,
                                   const std::vector<miniLepton> &,
                                   const std::vector<miniLepton> &, int, int,
                                   int, float, int &, TH1 *h_cutflow = 0);
    bool pass_2l1tau_inclusive_selection(const std::vector<miniLepton> &,
                                         const std::vector<miniLepton> &,
                                         const std::vector<miniLepton> &,
                                         const std::vector<miniTau> &, int, int,
                                         int, float, TH1 *h_cutflow = 0);
    bool pass_2lss1tau_SR_selection(const std::vector<miniLepton> &,
                                    const std::vector<miniTau> &);
    bool pass_2lss1tau_FakeAR_selection(const std::vector<miniLepton> &,
                                        const std::vector<miniTau> &);
    bool pass_2lss1tau_FlipAR_selection(const std::vector<miniLepton> &,
                                        const std::vector<miniTau> &);

    bool pass_2lss1tau_tightLepID(const std::vector<miniLepton> &);
    bool pass_2lss1tau_2lss(const std::vector<miniLepton> &);
    bool pass_2lss1tau_taucharge(const miniTau &, const miniLepton &);
    bool pass_2lss1tau_tauNumber(const std::vector<miniTau> &);
    // TODO
    bool pass_2lss1tau_FakeCR_selection();
    bool pass_2lss1tau_FakeARCR_selection();
    bool pass_2lss1tau_FlipCR_selection();
    bool pass_2lss1tau_FlipARCR_selection();
    //
    bool pass_2lss_inclusive_CR_selection();
    bool pass_2lss_ttW_CR_selection();
    bool pass_2lss_NJet_CR_selection();

    bool pass_3l_generic_selection(const std::vector<miniLepton> &,
                                   const std::vector<miniLepton> &, int, float,
                                   int &, TH1 *h_cutflow = 0);
    bool pass_3l1tau_inclusive_selection(const std::vector<miniLepton> &,
                                         const std::vector<miniLepton> &,
                                         const std::vector<miniTau> &, int, int,
                                         int, float, TH1 *h_cutflow = 0);
    bool pass_3l1tau_SR_selection(const std::vector<miniLepton> &,
                                  const std::vector<miniTau> &);
    bool pass_3l1tau_FakeAR_selection(const std::vector<miniLepton> &,
                                      const std::vector<miniTau> &);
    bool pass_3l1tau_CR_selection(const std::vector<miniLepton> &,
                                  const std::vector<miniTau> &);
    bool pass_3l1tau_FakeARCR_selection(const std::vector<miniLepton> &,
                                        const std::vector<miniTau> &);
    bool pass_3l1tau_tightID(const std::vector<miniLepton> &);
    bool pass_3l1tau_charge(const std::vector<miniLepton> &,
                            const std::vector<miniTau> &);
    bool pass_3l1tau_tauNumber(const std::vector<miniTau> &);

    // TODO
    bool pass_3l_inclusive_CR_selection();
    bool pass_3l_ttZ_CR_selection();
    bool pass_3l_WZ_CR_selection();

    bool pass_pairMass_veto(const std::vector<miniLepton> &);
    bool pass_Zmass_veto(const std::vector<miniLepton> &, bool, bool);
    bool pass_metLD_3l(float, const std::vector<miniLepton> &, int);

    //////////////////////////
    // Deprecated methods
    //////////////////////////
    bool pass_lepton_number(const std::vector<miniLepton> &,
                            const std::vector<miniLepton> &);
    bool pass_lepton_pt(const std::vector<miniLepton> &);
    bool pass_lepton_ID(bool, bool lep1IsTight = false,
                        bool lep2IsTight = false);
    bool pass_lepton_charge(int, int);

    bool pass_tight_charge(const std::vector<miniLepton> &);
    bool pass_Zmass_veto(const std::vector<miniLepton> &);
    bool pass_metLD(float, const std::vector<miniLepton> &, int njets = 0);
    bool pass_tau_number(int);
    bool pass_charge_sum(int, const std::vector<miniLepton> &);
    bool pass_tau_charge(int, const std::vector<miniLepton> &);
    bool pass_taupair_charge(int, int);
    bool pass_tau_ID(int);
    bool pass_lep_tau_ID(bool, int);
    bool pass_jet_number(int);
    bool pass_btag_number(int, int);
    bool pass_lep_mc_match(const miniLepton &);
    bool pass_lep_mc_match(const std::vector<miniLepton> &, int nleps = 2);
    bool pass_tau_mc_match(const pat::Tau &);
    bool pass_tau_mc_match(const std::vector<pat::Tau> &, int ntaus = 2);

  protected:
    Analysis_types anaType_;
    Selection_types selType_;
    bool verbose_;
    bool isMC_;
    bool looseselection_;
};

#endif
