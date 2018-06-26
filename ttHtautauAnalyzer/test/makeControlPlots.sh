#!/bin/bash

# example
plotDatacards.py "1l_2#tau CR" ~/nobackup/control_histograms/mTauTauVis_1l2tau_ctrl_41p53invfb.root mTauTauVis -c ttH TTW TTWW TTZ EWK Rares conversions fakes_data data_obs -v -o ~/public_html/Datacards/CR_plots/mTauTauVis_1l2tau_fakeCR.pdf -x "m_{#tau#tau}^{vis} [GeV]"

plotDatacards.py "1l_2#tau CR" ~/nobackup/control_histograms/nJet_1l2tau_ctrl_41p53invfb.root nJet -c ttH TTW TTWW TTZ EWK Rares conversions fakes_data data_obs -v -o ~/public_html/Datacards/CR_plots/nJet_1l2tau_fakeCR.pdf

plotDatacards.py "1l_2#tau CR" ~/nobackup/control_histograms/mvaOutput_1l2tau_ctrl_41p53invfb_binned.root x -c ttH TTW TTWW TTZ EWK Rares conversions fakes_data data_obs -v -o ~/public_html/Datacards/CR_plots/bdt_1l2tau_fakeCR.pdf 
