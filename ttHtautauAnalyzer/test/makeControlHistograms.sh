#!/bin/bash

variable=${1}
xmin=${2}
xmax=${3}
nbins=${4:-100}
mvaNtupleVersion=${5:-jun2018v3}

makeHistograms.py control_1l2tau $variable $xmin $xmax -b $nbins --channels ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs --luminosity 41.5 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --version $mvaNtupleVersion -o ${variable}_1l2tau_ctrl_41p53invfb.root

if [ "$variable" = "mvaOutput" ]; then
	binDatacards.py 1l2tau ${variable}_1l2tau_ctrl_41p53invfb.root -o ${variable}_1l2tau_ctrl_41p53invfb_binned.root -p -c ttH TTW TTWW TTZ EWK Rares conversions fakes_data data_obs --nosystematics
fi

#makeHistograms.py control_2lss1tau $variable $xmin $xmax --channels ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data flips_data data_obs --luminosity 41.5 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --version $mvaNtupleVersion -o ${variable}_2lss1tau_ctrl_41p53invfb.root

#makeHistograms.py control_3l1tau $variable $xmin $xmax --channels ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs --luminosity 41.5 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --version $mvaNtupleVersion -o ${variable}_3l1tau_ctrl_41p53invfb.root

#makeHistograms.py control_2l2tau $variable $xmin $xmax --channels ttH TTW TTWW TTZ EWK Rares tH VH conversions ggH fakes_data data_obs --luminosity 41.5 -d ../dataFiles/DatasetList_2017reMiniAODv2.csv -vvv --version $mvaNtupleVersion -o ${variable}_2l2tau_ctrl_41p53invfb.root 
