#ifndef MVAVars_h
#define MVAVars_h

//#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TMVA/Reader.h"
#include "TMath.h"
#include "TVector3.h"
#include "Types_enum.h"
#include "miniLepton.h"
#include "miniTau.h"

#include <algorithm>
#include <cmath>
#include <vector>

class MVAVars
{
  public:
    // MVAVars(const std::vector<miniLepton>&, const
    // std::vector<TLorentzVector>&,
    //		const std::vector<TLorentzVector>&, float, float, float);
    MVAVars(Analysis_types type)
    {
        std::cout << "This class is deprecated. All functionalities were moved "
                     "to mvaNtuple class."
                  << std::endl;
        anaType_ = type;
    }
    ~MVAVars(){};

    void set_up_tmva_reader();

    // void compute_all_variables(const std::vector<miniLepton>&,
    ///						   const std::vector<TLorentzVector>&,
    //						   const std::vector<TLorentzVector>&,
    //						   float, float, float, int);
    void compute_all_variables(const std::vector<miniLepton> &,
                               const std::vector<miniTau> &,
                               const std::vector<TLorentzVector> &, float,
                               float, float, int);
    void compute_taudecay_variables(const std::vector<miniTau> &);
    void compute_taudecay_variables(const std::vector<TLorentzVector> &,
                                    const std::vector<TLorentzVector> &,
                                    const std::vector<TLorentzVector> &,
                                    const std::vector<int> decaymode);

    // int nJet() const {return nJet_;}
    float nJet() const { return nJet_; }
    float mindr_lep0_jet() const { return mindr_lep0_jet_; }
    float mindr_lep1_jet() const { return mindr_lep1_jet_; }
    float mindr_lep2_jet() const { return mindr_lep2_jet_; }
    float avg_dr_jet() const { return avg_dr_jet_; }
    float max_lep_eta() const { return max_lep_eta_; }
    float met() const { return met_; }
    float mht() const { return mht_; }
    float mT_met_lep0() const { return mT_met_lep0_; }
    float lep0_conept() const { return lep0_conept_; }
    float lep1_conept() const { return lep1_conept_; }
    float lep2_conept() const { return lep2_conept_; }
    // bool isGenMatched() const {return isGenMatched_;}
    float dr_leps() const { return dr_leps_; }
    float tau_pt() const { return tau_pt_; }
    float dr_lep0_tau() const { return dr_lep0_tau_; }
    float dr_lep1_tau() const { return dr_lep1_tau_; }
    float mvis_lep0_tau() const { return mvis_lep0_tau_; }
    float mvis_lep1_tau() const { return mvis_lep1_tau_; }
    float ht() const { return ht_; }
    float dr_taus() const { return dr_taus_; }
    float mvis_taus() const { return mvis_taus_; }
    float pt_taus() const { return pt_taus_; }
    float max_dr_jet() const { return max_dr_jet_; }
    float tau0_pt() const { return tau0_pt_; }
    float tau1_pt() const { return tau1_pt_; }

    float costS_tau() const { return costS_tau_; }
    float dr_lep_tau0() const { return dr_lep_tau0_; }
    float dr_lep_tau1() const { return dr_lep_tau1_; }
    float dr_lep_tau_ss() const { return dr_lep_tau_ss_; }
    float mindr_tau0_jet() const { return mindr_tau0_jet_; }
    float mindr_tau1_jet() const { return mindr_tau1_jet_; }
    float taup_cosPsi() const { return taup_cosPsi_; }
    float taum_cosPsi() const { return taum_cosPsi_; }
    float evistaus_diff() const { return taum_energy_ - taup_energy_; }
    float evistaus_sum() const { return taum_energy_ + taup_energy_; }

    int tau0_decaymode() const { return tau0_decaymode_; }
    int tau1_decaymode() const { return tau1_decaymode_; }
    float tau0_energy() const { return tau0_energy_; }
    float tau1_energy() const { return tau1_energy_; }
    float tau0_upsilon() const { return tau0_upsilon_; }
    float tau1_upsilon() const { return tau1_upsilon_; }
    int taup_decaymode() const { return taup_decaymode_; }
    int taum_decaymode() const { return taum_decaymode_; }
    float taup_energy() const { return taup_energy_; }
    float taum_energy() const { return taum_energy_; }
    float taup_upsilon() const { return taup_upsilon_; }
    float taum_upsilon() const { return taum_upsilon_; }

    float BDT_ttV();
    float BDT_ttbar();

  private:
    // int nJet_;
    float nJet_;
    float mindr_lep0_jet_;
    float mindr_lep1_jet_;
    float mindr_lep2_jet_;
    float avg_dr_jet_;
    float max_lep_eta_;
    float met_;
    float mht_;
    float mT_met_lep0_;
    float lep0_conept_;
    float lep1_conept_;
    float lep2_conept_;
    // bool isGenMatched_;
    float dr_leps_;
    float tau_pt_;
    float dr_lep0_tau_;
    float dr_lep1_tau_;
    float mvis_lep0_tau_;
    float mvis_lep1_tau_;

    float ht_;
    float dr_taus_;
    float mvis_taus_;
    float pt_taus_;
    float max_dr_jet_;
    float tau0_pt_;
    float tau1_pt_;
    float ntags_loose_;
    // float ntags_medium_;

    float costS_tau_;
    float dr_lep_tau0_;
    float dr_lep_tau1_;
    float dr_lep_tau_ss_;
    float mindr_tau0_jet_;
    float mindr_tau1_jet_;

    float taup_cosPsi_;
    float taum_cosPsi_;

    int tau0_decaymode_;
    int tau1_decaymode_;
    float tau0_energy_;
    float tau1_energy_;
    float tau0_upsilon_;
    float tau1_upsilon_;
    int taup_decaymode_;
    int taum_decaymode_;
    float taup_energy_;
    float taum_energy_;
    float taup_upsilon_;
    float taum_upsilon_;

    float compute_average_dr(const std::vector<TLorentzVector> &);
    float compute_max_dr(const std::vector<TLorentzVector> &);
    float compute_mindr(const TLorentzVector &,
                        const std::vector<TLorentzVector> &);
    float compute_cosThetaS(const TLorentzVector &);
    float compute_upsilon(const TLorentzVector &, const TLorentzVector &);
    // float compute_chargedEnergyAsym();
    float compute_cosPsi(const miniTau &, float mass = 0.139);
    float compute_cosPsi(const TLorentzVector &, const TLorentzVector &,
                         const TLorentzVector &, float mass = 0.139);
    float lam(float, float, float);

    Analysis_types anaType_;

    TMVA::Reader *reader_ttV;
    TMVA::Reader *reader_ttbar;
};

#endif
